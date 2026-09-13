use cosmolkit_core::{
    __migration_conjugation::{ConjugationError, assign_conjugation, atom_has_conjugated_bond},
    AromaticityError, ValenceAssignment, ValenceError,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, Element};

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn bond(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
}

fn valence(explicit_valence: &[i32], implicit_hydrogens: &[i32]) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: explicit_valence.to_vec(),
        implicit_hydrogens: implicit_hydrogens.to_vec(),
    }
}

fn flags(topology: &TopologyBlock) -> Vec<bool> {
    topology.bonds.iter().map(Bond::is_conjugated).collect()
}

fn assert_only_flags_changed(
    source: &TopologyBlock,
    actual: &TopologyBlock,
    expected_flags: &[bool],
) {
    let mut expected = source.clone();
    for (bond, &flag) in expected.bonds.iter_mut().zip(expected_flags) {
        bond.set_conjugated(flag);
    }
    assert_eq!(actual, &expected);
    actual.validate().unwrap();
}

fn three_atom_motif(
    terminal: AtomSpec,
    terminal_explicit: i32,
    terminal_implicit: i32,
) -> (TopologyBlock, ValenceAssignment) {
    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            terminal,
        ],
        vec![bond(0, 1, BondOrder::Double), bond(1, 2, BondOrder::Single)],
    );
    (
        graph,
        valence(&[2, 3, terminal_explicit], &[2, 1, terminal_implicit]),
    )
}

#[test]
fn reset_is_exact_for_empty_isolated_aromatic_and_stale_disconnected_rows() {
    let empty = TopologyBlock::default();
    assert_eq!(
        assign_conjugation(&empty, &valence(&[], &[])).unwrap(),
        empty
    );

    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::N)
                .with_prop("row", "isolated")
                .unwrap(),
        ],
        vec![
            bond(0, 1, BondOrder::Aromatic)
                .with_aromatic(true)
                .with_conjugated(false),
            bond(2, 3, BondOrder::Single).with_conjugated(true),
        ],
    );
    let snapshot = graph.clone();
    let result = assign_conjugation(&graph, &valence(&[1, 1, 1, 1, 0], &[0; 5])).unwrap();
    assert_only_flags_changed(&snapshot, &result, &[true, false]);
    assert_eq!(graph, snapshot);
}

#[test]
fn issue_539_allyl_and_terminal_carbon_oxygen_nitrogen_rows_match_source() {
    let rows = [
        (AtomSpec::new(Element::C), 1, 3, vec![false, false], "C=C-C"),
        (AtomSpec::new(Element::O), 1, 1, vec![true, true], "C=C-O"),
        (AtomSpec::new(Element::N), 1, 2, vec![true, true], "C=C-N"),
        (
            AtomSpec::new(Element::N)
                .with_formal_charge(1)
                .with_explicit_hydrogens(3)
                .with_no_implicit(true),
            4,
            0,
            vec![false, false],
            "C=C-[NH3+]",
        ),
    ];
    for (terminal, explicit, implicit, expected, label) in rows {
        let (graph, assignment) = three_atom_motif(terminal, explicit, implicit);
        let snapshot = graph.clone();
        let result = assign_conjugation(&graph, &assignment).unwrap();
        assert_eq!(flags(&result), expected, "{label}");
        assert_only_flags_changed(&snapshot, &result, &expected);
        assert_eq!(graph, snapshot, "{label}");
    }

    let allyl = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C)
                .with_formal_charge(1)
                .with_explicit_hydrogens(2)
                .with_no_implicit(true),
        ],
        vec![bond(0, 1, BondOrder::Double), bond(1, 2, BondOrder::Single)],
    );
    let result = assign_conjugation(&allyl, &valence(&[2, 3, 3], &[2, 1, 0])).unwrap();
    assert_eq!(flags(&result), [true, true]);
}

#[test]
fn cyclopropenyl_and_tropylium_post_aromaticity_rows_remain_fully_conjugated() {
    for size in [3usize, 7] {
        let atoms = (0..size)
            .map(|index| {
                let spec = AtomSpec::new(Element::C).with_aromatic(true);
                if index + 1 == size {
                    spec.with_formal_charge(1)
                } else {
                    spec
                }
            })
            .collect();
        let bonds = (0..size)
            .map(|index| {
                bond(index, (index + 1) % size, BondOrder::Aromatic)
                    .with_aromatic(true)
                    .with_conjugated(false)
            })
            .collect();
        let graph = topology(atoms, bonds);
        let assignment = valence(&vec![3; size], &vec![0; size]);
        let result = assign_conjugation(&graph, &assignment).unwrap();
        assert_eq!(flags(&result), vec![true; size], "ring size {size}");
        assert_only_flags_changed(&graph, &result, &vec![true; size]);
    }
}

#[test]
fn methyl_and_fluoro_aromatic_attachment_rows_check_every_bond() {
    let mut methyl_atoms = vec![AtomSpec::new(Element::C)];
    methyl_atoms.extend((0..6).map(|_| AtomSpec::new(Element::C).with_aromatic(true)));
    let mut methyl_bonds = vec![bond(0, 1, BondOrder::Single).with_conjugated(true)];
    methyl_bonds.extend((0..6).map(|index| {
        bond(1 + index, 1 + (index + 1) % 6, BondOrder::Aromatic)
            .with_aromatic(true)
            .with_conjugated(false)
    }));
    let methyl = topology(methyl_atoms, methyl_bonds);
    let methyl_result = assign_conjugation(
        &methyl,
        &valence(&[1, 4, 3, 3, 3, 3, 3], &[3, 0, 1, 1, 1, 1, 1]),
    )
    .unwrap();
    assert_eq!(
        flags(&methyl_result),
        [false, true, true, true, true, true, true]
    );

    let fluoro = topology(
        vec![
            AtomSpec::new(Element::F),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_explicit_hydrogens(1),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::O),
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_explicit_hydrogens(1),
            AtomSpec::new(Element::C).with_aromatic(true),
            AtomSpec::new(Element::O),
        ],
        vec![
            bond(0, 1, BondOrder::Single).with_conjugated(true),
            bond(1, 2, BondOrder::Aromatic).with_aromatic(true),
            bond(2, 3, BondOrder::Aromatic).with_aromatic(true),
            bond(3, 4, BondOrder::Aromatic).with_aromatic(true),
            bond(4, 6, BondOrder::Aromatic).with_aromatic(true),
            bond(6, 7, BondOrder::Aromatic).with_aromatic(true),
            bond(7, 1, BondOrder::Aromatic).with_aromatic(true),
            bond(4, 5, BondOrder::Double),
            bond(7, 8, BondOrder::Double),
        ],
    );
    let fluoro_result =
        assign_conjugation(&fluoro, &valence(&[1, 4, 3, 3, 4, 2, 3, 4, 2], &[0; 9])).unwrap();
    assert_eq!(
        flags(&fluoro_result),
        [false, true, true, true, true, true, true, true, true]
    );
}

#[test]
fn issue_211_and_charged_hypervalence_candidate_gates_are_independent() {
    let (phosphorus, phosphorus_valence) = three_atom_motif(AtomSpec::new(Element::P), 3, 0);
    assert_eq!(
        flags(&assign_conjugation(&phosphorus, &phosphorus_valence).unwrap()),
        [false, false]
    );

    let (arsenic, arsenic_valence) = three_atom_motif(AtomSpec::new(Element::AS), 3, 0);
    assert_eq!(
        flags(&assign_conjugation(&arsenic, &arsenic_valence).unwrap()),
        [false, false]
    );

    let (silicon, silicon_valence) = three_atom_motif(AtomSpec::new(Element::SI), 1, 1);
    assert_eq!(
        flags(&assign_conjugation(&silicon, &silicon_valence).unwrap()),
        [true, true]
    );

    let (sulfur_degree_one, sulfur_degree_one_valence) =
        three_atom_motif(AtomSpec::new(Element::S), 1, 0);
    assert_eq!(
        flags(&assign_conjugation(&sulfur_degree_one, &sulfur_degree_one_valence).unwrap()),
        [true, true]
    );
    let (sulfur_degree_two, sulfur_degree_two_valence) =
        three_atom_motif(AtomSpec::new(Element::S).with_explicit_hydrogens(1), 2, 0);
    assert_eq!(
        flags(&assign_conjugation(&sulfur_degree_two, &sulfur_degree_two_valence).unwrap()),
        [false, false]
    );

    let (selenium_degree_one, selenium_degree_one_valence) =
        three_atom_motif(AtomSpec::new(Element::SE), 1, 0);
    assert_eq!(
        flags(&assign_conjugation(&selenium_degree_one, &selenium_degree_one_valence).unwrap()),
        [true, true]
    );

    let (neutral_hypervalent, neutral_assignment) =
        three_atom_motif(AtomSpec::new(Element::C), 5, 0);
    assert_eq!(
        flags(&assign_conjugation(&neutral_hypervalent, &neutral_assignment).unwrap()),
        [false, false]
    );
    let (charged_hypervalent, charged_assignment) =
        three_atom_motif(AtomSpec::new(Element::C).with_formal_charge(1), 5, 0);
    assert_eq!(
        flags(&assign_conjugation(&charged_hypervalent, &charged_assignment).unwrap()),
        [true, true]
    );
}

#[test]
fn substitution_bounds_and_multiple_bond_threshold_are_exact() {
    let ethene = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::C)],
        vec![bond(0, 1, BondOrder::Double)],
    );
    assert_eq!(
        flags(&assign_conjugation(&ethene, &valence(&[2, 2], &[2, 2])).unwrap()),
        [false]
    );

    for center_implicit in [0, 1] {
        let (graph, mut assignment) = three_atom_motif(AtomSpec::new(Element::O), 1, 1);
        assignment.implicit_hydrogens[1] = center_implicit;
        assert_eq!(
            flags(&assign_conjugation(&graph, &assignment).unwrap()),
            [true, true],
            "center substitutions {}",
            2 + center_implicit
        );
    }

    let (four_substitutions, mut assignment) = three_atom_motif(AtomSpec::new(Element::O), 1, 1);
    assignment.implicit_hydrogens[1] = 2;
    assert_eq!(
        flags(&assign_conjugation(&four_substitutions, &assignment).unwrap()),
        [false, false]
    );

    let single_single = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::O),
        ],
        vec![bond(0, 1, BondOrder::Single), bond(1, 2, BondOrder::Single)],
    );
    assert_eq!(
        flags(&assign_conjugation(&single_single, &valence(&[1, 2, 1], &[3, 2, 1])).unwrap()),
        [false, false]
    );
}

#[test]
fn concrete_bond_order_contributions_are_not_generalized() {
    for order in [
        BondOrder::OneAndHalf,
        BondOrder::Aromatic,
        BondOrder::Double,
        BondOrder::Triple,
    ] {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
            ],
            vec![
                bond(0, 1, order).with_aromatic(order == BondOrder::Aromatic),
                bond(1, 2, BondOrder::Single),
            ],
        );
        let result = assign_conjugation(&graph, &valence(&[2, 3, 1], &[1, 0, 1])).unwrap();
        assert_eq!(flags(&result), [true, true], "order {order:?}");
    }

    for order in [
        BondOrder::Single,
        BondOrder::Zero,
        BondOrder::Dative,
        BondOrder::DativeOne,
    ] {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
            ],
            vec![bond(0, 1, order), bond(1, 2, BondOrder::Single)],
        );
        let result = assign_conjugation(&graph, &valence(&[1, 2, 1], &[1, 0, 1])).unwrap();
        assert_eq!(flags(&result), [false, false], "order {order:?}");
    }

    for order in [BondOrder::DativeLeft, BondOrder::DativeRight] {
        let graph = topology(
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::O),
            ],
            vec![bond(0, 1, order), bond(1, 2, BondOrder::Single)],
        );
        assert!(matches!(
            assign_conjugation(&graph, &valence(&[1, 2, 1], &[1, 0, 1])),
            Err(ConjugationError::Aromaticity(AromaticityError::Valence(
                ValenceError::BadBondType {
                    bond: Some(bond_id),
                    order: actual_order,
                }
            ))) if bond_id == BondId::new(0) && actual_order == order
        ));
    }
}

#[test]
fn radical_and_charge_electron_gates_preserve_source_distinctions() {
    let (neutral_carbon_radical, assignment) =
        three_atom_motif(AtomSpec::new(Element::C).with_radical_electrons(1), 1, 1);
    assert_eq!(
        flags(&assign_conjugation(&neutral_carbon_radical, &assignment).unwrap()),
        [true, true]
    );

    let (charged_carbon_radical, assignment) = three_atom_motif(
        AtomSpec::new(Element::C)
            .with_radical_electrons(1)
            .with_formal_charge(1),
        1,
        1,
    );
    assert_eq!(
        flags(&assign_conjugation(&charged_carbon_radical, &assignment).unwrap()),
        [true, true]
    );

    let (nitrogen_radical, assignment) =
        three_atom_motif(AtomSpec::new(Element::N).with_radical_electrons(2), 1, 1);
    assert_eq!(
        flags(&assign_conjugation(&nitrogen_radical, &assignment).unwrap()),
        [true, true]
    );

    let (nitrogen_three_radicals, assignment) =
        three_atom_motif(AtomSpec::new(Element::N).with_radical_electrons(3), 1, 1);
    assert_eq!(
        flags(&assign_conjugation(&nitrogen_three_radicals, &assignment).unwrap()),
        [false, false]
    );
}

#[test]
fn invalid_inputs_report_exact_fields_without_mutating_sources() {
    let graph = topology(
        vec![AtomSpec::new(Element::C), AtomSpec::new(Element::O)],
        vec![bond(0, 1, BondOrder::Double).with_conjugated(true)],
    );
    let snapshot = graph.clone();

    assert!(matches!(
        assign_conjugation(&graph, &valence(&[2], &[0, 0])),
        Err(ConjugationError::ValenceAssignmentLength {
            field: "explicit_valence",
            actual: 1,
            expected: 2,
        })
    ));
    assert!(matches!(
        assign_conjugation(&graph, &valence(&[2, 2], &[0])),
        Err(ConjugationError::ValenceAssignmentLength {
            field: "implicit_hydrogens",
            actual: 1,
            expected: 2,
        })
    ));
    assert!(matches!(
        assign_conjugation(&graph, &valence(&[-1, 2], &[0, 0])),
        Err(ConjugationError::InvalidValenceRow {
            atom,
            field: "explicit_valence",
            value: -1,
        }) if atom == AtomId::new(0)
    ));
    assert!(matches!(
        assign_conjugation(&graph, &valence(&[2, 2], &[0, -1])),
        Err(ConjugationError::InvalidValenceRow {
            atom,
            field: "implicit_hydrogens",
            value: -1,
        }) if atom == AtomId::new(1)
    ));
    assert!(matches!(
        assign_conjugation(&graph, &valence(&[i32::MAX, 2], &[1, 0])),
        Err(ConjugationError::IntegerOverflow {
            field: "total valence"
        })
    ));
    assert!(matches!(
        atom_has_conjugated_bond(&graph, AtomId::new(2)),
        Err(ConjugationError::AtomOutOfRange {
            atom,
            atom_count: 2,
        }) if atom == AtomId::new(2)
    ));
    assert_eq!(graph, snapshot);

    let mut malformed = graph.clone();
    malformed.atoms[0] = Atom::from_spec(AtomId::new(9), AtomSpec::new(Element::C));
    assert!(matches!(
        assign_conjugation(&malformed, &valence(&[2, 2], &[0, 0])),
        Err(ConjugationError::InvalidTopology(_))
    ));
}

#[test]
fn repeated_assignment_predicate_and_complete_row_identity_are_deterministic() {
    let atom_zero = AtomSpec::new(Element::C)
        .with_prop("atom-key", "atom-value")
        .unwrap();
    let graph = topology(
        vec![
            atom_zero,
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::O),
        ],
        vec![
            bond(0, 1, BondOrder::Double)
                .with_prop("bond-key", "bond-value")
                .unwrap(),
            bond(1, 2, BondOrder::Single),
        ],
    );
    let assignment = valence(&[2, 3, 1], &[2, 1, 1]);
    let first = assign_conjugation(&graph, &assignment).unwrap();
    let second = assign_conjugation(&graph, &assignment).unwrap();
    let third = assign_conjugation(&first, &assignment).unwrap();
    assert_eq!(first, second);
    assert_eq!(first, third);
    assert_only_flags_changed(&graph, &first, &[true, true]);
    assert!(atom_has_conjugated_bond(&first, AtomId::new(0)).unwrap());
    assert!(atom_has_conjugated_bond(&first, AtomId::new(1)).unwrap());
    assert!(atom_has_conjugated_bond(&first, AtomId::new(2)).unwrap());

    let isolated = topology(vec![AtomSpec::new(Element::C)], vec![]);
    assert!(!atom_has_conjugated_bond(&isolated, AtomId::new(0)).unwrap());
}
