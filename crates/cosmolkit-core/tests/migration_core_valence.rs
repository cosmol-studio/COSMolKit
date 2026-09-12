use cosmolkit_core::{
    ValenceError, ValenceModel, ValenceParams, ValencePhase, assign_valence, bond_type_as_double,
    bond_valence_contrib, calculate_explicit_valence_from_parts, explicit_valence_for_atom,
    has_valence_violation, implicit_valence_for_atom,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn atom_spec(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).expect("modeled element"))
}

fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
    let atoms = atom_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect::<Vec<_>>();
    let bonds = bond_specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
        .collect::<Vec<_>>();
    TopologyBlock {
        adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
        atoms,
        bonds,
        ..TopologyBlock::default()
    }
}

fn star(center: AtomSpec, neighbor_count: usize, order: BondOrder) -> TopologyBlock {
    let mut atoms = vec![center];
    atoms.extend((0..neighbor_count).map(|_| atom_spec(9).with_no_implicit(true)));
    let bonds = (1..=neighbor_count)
        .map(|neighbor| BondSpec::new(AtomId::new(0), AtomId::new(neighbor), order))
        .collect();
    topology(atoms, bonds)
}

#[test]
fn bond_order_matrix_matches_rdkit_numeric_and_rejection_branches() {
    for (order, expected) in [
        (BondOrder::Unspecified, 0.0),
        (BondOrder::Ionic, 0.0),
        (BondOrder::Zero, 0.0),
        (BondOrder::Single, 1.0),
        (BondOrder::Double, 2.0),
        (BondOrder::Triple, 3.0),
        (BondOrder::Quadruple, 4.0),
        (BondOrder::Quintuple, 5.0),
        (BondOrder::Hextuple, 6.0),
        (BondOrder::OneAndHalf, 1.5),
        (BondOrder::TwoAndHalf, 2.5),
        (BondOrder::ThreeAndHalf, 3.5),
        (BondOrder::FourAndHalf, 4.5),
        (BondOrder::FiveAndHalf, 5.5),
        (BondOrder::Aromatic, 1.5),
        (BondOrder::Dative, 1.0),
        (BondOrder::DativeOne, 1.0),
        (BondOrder::Hydrogen, 0.0),
    ] {
        assert_eq!(bond_type_as_double(order).unwrap(), expected, "{order:?}");
    }
    for order in [
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
        BondOrder::ThreeCenter,
        BondOrder::Other,
    ] {
        assert_eq!(
            bond_type_as_double(order),
            Err(ValenceError::BadBondType { bond: None, order })
        );
    }
}

#[test]
fn bond_contribution_is_endpoint_and_dative_direction_sensitive() {
    let dative = Bond::from_spec(
        BondId::new(7),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::DativeOne),
    );
    assert_eq!(bond_valence_contrib(&dative, AtomId::new(0)).unwrap(), 0.0);
    assert_eq!(bond_valence_contrib(&dative, AtomId::new(1)).unwrap(), 1.0);
    assert_eq!(bond_valence_contrib(&dative, AtomId::new(2)).unwrap(), 0.0);

    let rejected = Bond::from_spec(
        BondId::new(9),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::ThreeCenter),
    );
    assert_eq!(
        bond_valence_contrib(&rejected, AtomId::new(1)),
        Err(ValenceError::BadBondType {
            bond: Some(BondId::new(9)),
            order: BondOrder::ThreeCenter,
        })
    );
}

#[test]
fn canonical_defaults_assign_in_stable_row_order_without_mutating_input() {
    assert_eq!(
        ValenceParams::default(),
        ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: true,
        }
    );
    let input = topology(
        vec![atom_spec(6), atom_spec(6), atom_spec(8)],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        ],
    );
    let before = input.clone();
    let assignment = assign_valence(&input, &ValenceParams::default()).unwrap();
    assert_eq!(assignment.explicit_valence, vec![1, 2, 1]);
    assert_eq!(assignment.implicit_hydrogens, vec![3, 2, 1]);
    assert_eq!(assignment.explicit_valence.len(), input.atoms.len());
    assert_eq!(assignment.implicit_hydrogens.len(), input.atoms.len());
    assert_eq!(input, before);
}

#[test]
fn canonical_atom_access_rejects_range_before_any_lookup() {
    let input = topology(vec![atom_spec(6)], vec![]);
    assert_eq!(
        explicit_valence_for_atom(&input, AtomId::new(3), true),
        Err(ValenceError::AtomOutOfRange {
            atom: AtomId::new(3),
            atom_count: 1,
        })
    );
    assert_eq!(
        has_valence_violation(&input, AtomId::new(4)),
        Err(ValenceError::AtomOutOfRange {
            atom: AtomId::new(4),
            atom_count: 1,
        })
    );
}

#[test]
fn canonical_assignment_preserves_structured_topology_error() {
    let bad_atom = Atom::from_spec(AtomId::new(1), atom_spec(6));
    let input = TopologyBlock {
        atoms: vec![bad_atom],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };
    assert_eq!(
        assign_valence(&input, &ValenceParams::default()),
        Err(ValenceError::InvalidTopology {
            source: TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            },
        })
    );
}

#[test]
fn raw_parts_helpers_reject_dangling_adjacency_bond_and_endpoint() {
    let atoms = vec![
        Atom::from_spec(AtomId::new(0), atom_spec(6)),
        Atom::from_spec(AtomId::new(1), atom_spec(6)),
        Atom::from_spec(AtomId::new(2), atom_spec(6)),
    ];
    let source_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    let adjacency = AdjacencyList::from_topology(3, std::slice::from_ref(&source_bond));
    assert_eq!(
        calculate_explicit_valence_from_parts(&atoms, &[], &adjacency, AtomId::new(0), true, false),
        Err(ValenceError::AdjacencyBondOutOfRange {
            atom: AtomId::new(0),
            bond: BondId::new(0),
            bond_count: 0,
        })
    );

    let different_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
    );
    assert_eq!(
        calculate_explicit_valence_from_parts(
            &atoms,
            &[different_bond],
            &adjacency,
            AtomId::new(0),
            true,
            false,
        ),
        Err(ValenceError::AdjacencyEndpointMismatch {
            atom: AtomId::new(0),
            neighbor_atom: 1,
            bond: BondId::new(0),
            begin: AtomId::new(0),
            end: AtomId::new(2),
        })
    );
}

#[test]
fn raw_parts_helpers_reject_neighbor_rows_outside_the_atom_slice() {
    let atoms = vec![Atom::from_spec(AtomId::new(0), atom_spec(6))];
    let bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
    );
    let adjacency = AdjacencyList::from_topology(2, std::slice::from_ref(&bond));
    assert_eq!(
        calculate_explicit_valence_from_parts(
            &atoms,
            &[bond],
            &adjacency,
            AtomId::new(0),
            true,
            false,
        ),
        Err(ValenceError::AdjacencyAtomOutOfRange {
            atom: AtomId::new(0),
            neighbor_atom: 1,
            atom_count: 1,
        })
    );
}

#[test]
fn implicit_special_cases_cover_no_implicit_dummy_hydrogen_and_radicals() {
    let no_implicit = topology(
        vec![
            atom_spec(6)
                .with_no_implicit(true)
                .with_explicit_hydrogens(2),
        ],
        vec![],
    );
    assert_eq!(
        implicit_valence_for_atom(&no_implicit, AtomId::new(0), None, true).unwrap(),
        0
    );

    let dummy = topology(vec![atom_spec(0)], vec![]);
    assert_eq!(
        implicit_valence_for_atom(&dummy, AtomId::new(0), None, true).unwrap(),
        0
    );

    for (charge, expected) in [(1, 0), (-1, 0), (0, 1)] {
        let hydrogen = topology(vec![atom_spec(1).with_formal_charge(charge)], vec![]);
        assert_eq!(
            implicit_valence_for_atom(&hydrogen, AtomId::new(0), Some(0), true).unwrap(),
            expected
        );
    }

    let radical = topology(vec![atom_spec(6).with_radical_electrons(2)], vec![]);
    assert_eq!(
        implicit_valence_for_atom(&radical, AtomId::new(0), Some(0), true).unwrap(),
        2
    );
}

#[test]
fn implicit_parameter_rejects_the_internal_negative_cache_sentinel() {
    let input = topology(vec![atom_spec(6)], vec![]);
    assert_eq!(
        implicit_valence_for_atom(&input, AtomId::new(0), Some(-1), true),
        Err(ValenceError::InvalidExplicitValenceInput {
            atom: AtomId::new(0),
            value: -1,
        })
    );
    assert_eq!(
        implicit_valence_for_atom(&input, AtomId::new(0), None, true).unwrap(),
        4
    );
}

#[test]
fn multivalence_exact_next_higher_unlimited_and_non_strict_exhaustion_are_distinct() {
    let sulfur_exact = star(atom_spec(16), 2, BondOrder::Single);
    assert_eq!(
        implicit_valence_for_atom(&sulfur_exact, AtomId::new(0), None, true).unwrap(),
        0
    );
    let sulfur_next = star(atom_spec(16), 3, BondOrder::Single);
    assert_eq!(
        implicit_valence_for_atom(&sulfur_next, AtomId::new(0), None, true).unwrap(),
        1
    );

    let transition_metal = topology(vec![atom_spec(26)], vec![]);
    assert_eq!(
        implicit_valence_for_atom(&transition_metal, AtomId::new(0), None, true).unwrap(),
        0
    );

    let overbonded = star(atom_spec(6), 5, BondOrder::Single);
    let error = explicit_valence_for_atom(&overbonded, AtomId::new(0), true).unwrap_err();
    assert!(matches!(
        error,
        ValenceError::InvalidValence {
            atom,
            atomic_number: 6,
            formal_charge: 0,
            phase: ValencePhase::Explicit,
            calculated: Some(5),
            reason: "greater than permitted",
            ..
        } if atom == AtomId::new(0)
    ));
    assert_eq!(
        assign_valence(
            &overbonded,
            &ValenceParams {
                model: ValenceModel::RdkitLike,
                strict: false,
            },
        )
        .unwrap()
        .implicit_hydrogens[0],
        0
    );
}

#[test]
fn aromatic_half_order_rounding_and_satisfied_branches_are_stable() {
    let mut atoms = (0..6)
        .map(|_| atom_spec(6).with_aromatic(true))
        .collect::<Vec<_>>();
    atoms[0] = atoms[0].clone().with_no_implicit(true);
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        .into_iter()
        .map(|(begin, end)| {
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true)
        })
        .collect();
    let benzene = topology(atoms, bonds);
    let assignment = assign_valence(&benzene, &ValenceParams::default()).unwrap();
    assert_eq!(assignment.explicit_valence, vec![3; 6]);
    assert_eq!(assignment.implicit_hydrogens, vec![0, 1, 1, 1, 1, 1]);
}

#[test]
fn hypervalent_anions_and_two_coordinate_hydride_remain_accepted() {
    for atomic_number in [15, 16, 33, 34] {
        let input = star(
            atom_spec(atomic_number).with_formal_charge(-1),
            5,
            BondOrder::Single,
        );
        assert!(!has_valence_violation(&input, AtomId::new(0)).unwrap());
    }

    let hydride = topology(
        vec![
            atom_spec(1).with_formal_charge(-1).with_no_implicit(true),
            atom_spec(1).with_no_implicit(true),
            atom_spec(1).with_no_implicit(true),
        ],
        vec![
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
        ],
    );
    assert_eq!(
        explicit_valence_for_atom(&hydride, AtomId::new(0), true).unwrap(),
        2
    );
    assert!(!has_valence_violation(&hydride, AtomId::new(0)).unwrap());
}

#[test]
fn violation_classification_covers_dummy_hydrogen_charge_and_period_change() {
    let dummy = topology(vec![atom_spec(0).with_formal_charge(100)], vec![]);
    assert!(!has_valence_violation(&dummy, AtomId::new(0)).unwrap());

    for charge in [-2, 2] {
        let hydrogen = topology(vec![atom_spec(1).with_formal_charge(charge)], vec![]);
        assert!(has_valence_violation(&hydrogen, AtomId::new(0)).unwrap());
    }

    let fluorine_crossing_period = topology(vec![atom_spec(9).with_formal_charge(-2)], vec![]);
    assert!(has_valence_violation(&fluorine_crossing_period, AtomId::new(0)).unwrap());
    let excessive_positive = topology(vec![atom_spec(6).with_formal_charge(7)], vec![]);
    assert!(has_valence_violation(&excessive_positive, AtomId::new(0)).unwrap());
}
