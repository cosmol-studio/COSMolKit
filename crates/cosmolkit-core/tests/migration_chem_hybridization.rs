use cosmolkit_core::{
    __migration_hybridization::{HybridizationError, assign_hybridization},
    ValenceAssignment,
};
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
use cosmolkit_types::{BondOrder, ChiralTag, Element, Hybridization};

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

fn valence(explicit: &[i32], implicit: &[i32]) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: explicit.to_vec(),
        implicit_hydrogens: implicit.to_vec(),
    }
}

fn isolated(spec: AtomSpec, explicit_valence: i32, implicit: i32) -> Hybridization {
    let graph = topology(vec![spec], vec![]);
    assign_hybridization(&graph, &valence(&[explicit_valence], &[implicit]))
        .unwrap()
        .values[0]
}

fn tagged_carbon(tag: ChiralTag, hydrogens: u8, explicit_valence: i32) -> Hybridization {
    isolated(
        AtomSpec::new(Element::C)
            .with_chiral_tag(tag)
            .with_explicit_hydrogens(hydrogens)
            .with_no_implicit(true),
        explicit_valence,
        0,
    )
}

fn zero_bond_star(tag: ChiralTag, degree: usize) -> Hybridization {
    let mut atoms = vec![
        AtomSpec::new(Element::C)
            .with_chiral_tag(tag)
            .with_no_implicit(true),
    ];
    atoms.extend((0..degree).map(|_| AtomSpec::new(Element::H).with_no_implicit(true)));
    let bonds = (0..degree)
        .map(|index| bond(0, index + 1, BondOrder::Zero))
        .collect();
    let graph = topology(atoms, bonds);
    let mut explicit = vec![0; degree + 1];
    let implicit = vec![0; degree + 1];
    explicit[0] = 0;
    assign_hybridization(&graph, &valence(&explicit, &implicit))
        .unwrap()
        .values[0]
}

#[test]
fn empty_dummy_and_all_chiral_tags_replace_stale_values_in_stable_order() {
    assert!(
        assign_hybridization(&TopologyBlock::default(), &valence(&[], &[]))
            .unwrap()
            .values
            .is_empty()
    );

    let tags = [
        (ChiralTag::Unspecified, Hybridization::Sp3),
        (ChiralTag::TetrahedralCw, Hybridization::Sp3),
        (ChiralTag::TetrahedralCcw, Hybridization::Sp3),
        (ChiralTag::Other, Hybridization::Sp3),
        (ChiralTag::Tetrahedral, Hybridization::Sp3),
        (ChiralTag::Allene, Hybridization::Sp3),
        (ChiralTag::SquarePlanar, Hybridization::Sp2d),
        (ChiralTag::TrigonalBipyramidal, Hybridization::Sp3d),
        (ChiralTag::Octahedral, Hybridization::Sp3d2),
    ];
    let mut specs = vec![
        AtomSpec::new(Element::DUMMY)
            .with_chiral_tag(ChiralTag::Octahedral)
            .with_hybridization(Hybridization::Other),
    ];
    specs.extend(tags.iter().map(|(tag, _)| {
        AtomSpec::new(Element::C)
            .with_chiral_tag(*tag)
            .with_explicit_hydrogens(4)
            .with_no_implicit(true)
            .with_hybridization(Hybridization::Other)
    }));
    let graph = topology(specs, vec![]);
    let snapshot = graph.clone();
    let assignment = valence(&[0, 4, 4, 4, 4, 4, 4, 4, 4, 4], &[0; 10]);
    let first = assign_hybridization(&graph, &assignment).unwrap();
    let second = assign_hybridization(&graph, &assignment).unwrap();
    let expected: Vec<_> = std::iter::once(Hybridization::Unspecified)
        .chain(tags.into_iter().map(|(_, expected)| expected))
        .collect();
    assert_eq!(first.values, expected);
    assert_eq!(first, second);
    assert_eq!(graph, snapshot);
}

#[test]
fn coordination_tag_degree_bounds_use_exact_fast_paths_and_fallbacks() {
    assert_eq!(
        tagged_carbon(ChiralTag::SquarePlanar, 2, 2),
        Hybridization::Sp2d
    );
    assert_eq!(
        tagged_carbon(ChiralTag::SquarePlanar, 4, 4),
        Hybridization::Sp2d
    );
    assert_eq!(
        tagged_carbon(ChiralTag::SquarePlanar, 1, 1),
        Hybridization::Sp
    );
    assert_eq!(
        zero_bond_star(ChiralTag::SquarePlanar, 5),
        Hybridization::Sp
    );

    assert_eq!(
        tagged_carbon(ChiralTag::TrigonalBipyramidal, 2, 2),
        Hybridization::Sp3d
    );
    assert_eq!(
        tagged_carbon(ChiralTag::TrigonalBipyramidal, 5, 5),
        Hybridization::Sp3d
    );
    assert_eq!(
        tagged_carbon(ChiralTag::TrigonalBipyramidal, 1, 1),
        Hybridization::Sp
    );
    assert_eq!(
        zero_bond_star(ChiralTag::TrigonalBipyramidal, 6),
        Hybridization::Sp
    );

    assert_eq!(
        tagged_carbon(ChiralTag::Octahedral, 2, 2),
        Hybridization::Sp3d2
    );
    assert_eq!(
        tagged_carbon(ChiralTag::Octahedral, 6, 6),
        Hybridization::Sp3d2
    );
    assert_eq!(
        tagged_carbon(ChiralTag::Octahedral, 1, 1),
        Hybridization::Sp
    );
    assert_eq!(zero_bond_star(ChiralTag::Octahedral, 7), Hybridization::Sp);

    for tag in [
        ChiralTag::Tetrahedral,
        ChiralTag::TetrahedralCw,
        ChiralTag::TetrahedralCcw,
    ] {
        assert_eq!(tagged_carbon(tag, 4, 4), Hybridization::Sp3);
        assert_eq!(tagged_carbon(tag, 3, 3), Hybridization::Sp2);
    }
}

#[test]
fn source_basic_carbon_nitrogen_oxygen_and_charged_rows_are_sp3() {
    let rows = [
        (AtomSpec::new(Element::C), 4, 4, 0),
        (AtomSpec::new(Element::N), 3, 3, 0),
        (AtomSpec::new(Element::O), 2, 2, 0),
        (AtomSpec::new(Element::C).with_formal_charge(-2), 2, 2, 0),
        (
            AtomSpec::new(Element::C)
                .with_formal_charge(-1)
                .with_explicit_hydrogens(1),
            3,
            2,
            0,
        ),
    ];
    for (center, center_valence, center_bonds, center_implicit) in rows {
        let mut atoms = vec![center];
        atoms.extend((0..center_bonds).map(|_| AtomSpec::new(Element::H).with_no_implicit(true)));
        let bonds = (0..center_bonds)
            .map(|index| bond(0, index + 1, BondOrder::Single))
            .collect();
        let graph = topology(atoms, bonds);
        let mut explicit = vec![1; center_bonds + 1];
        let mut implicit = vec![0; center_bonds + 1];
        explicit[0] = center_valence;
        implicit[0] = center_implicit;
        assert_eq!(
            assign_hybridization(&graph, &valence(&explicit, &implicit))
                .unwrap()
                .values[0],
            Hybridization::Sp3
        );
    }
}

#[test]
fn issue_276_count_four_uses_conjugation_only_at_degree_three_or_less() {
    let oxygen = AtomSpec::new(Element::O)
        .with_explicit_hydrogens(1)
        .with_no_implicit(true);
    for (conjugated, expected) in [(false, Hybridization::Sp3), (true, Hybridization::Sp2)] {
        let graph = topology(
            vec![AtomSpec::new(Element::C), oxygen.clone()],
            vec![bond(0, 1, BondOrder::Single).with_conjugated(conjugated)],
        );
        assert_eq!(
            assign_hybridization(&graph, &valence(&[1, 2], &[3, 0]))
                .unwrap()
                .values[1],
            expected
        );
    }

    let mut atoms = vec![AtomSpec::new(Element::C).with_no_implicit(true)];
    atoms.extend((0..4).map(|_| AtomSpec::new(Element::H).with_no_implicit(true)));
    let mut bonds: Vec<_> = (0..4)
        .map(|index| bond(0, index + 1, BondOrder::Single))
        .collect();
    bonds[0] = bonds[0].clone().with_conjugated(true);
    let graph = topology(atoms, bonds);
    assert_eq!(
        assign_hybridization(&graph, &valence(&[4, 1, 1, 1, 1], &[0; 5]))
            .unwrap()
            .values[0],
        Hybridization::Sp3
    );
}

#[test]
fn heavy_elements_use_degree_only_and_cover_every_orbital_mapping() {
    let expected = [
        Hybridization::S,
        Hybridization::S,
        Hybridization::Sp,
        Hybridization::Sp2,
        Hybridization::Sp3,
        Hybridization::Sp3d,
        Hybridization::Sp3d2,
        Hybridization::Unspecified,
    ];
    for (degree, expected) in expected.into_iter().enumerate() {
        assert_eq!(
            isolated(
                AtomSpec::new(Element::AC)
                    .with_explicit_hydrogens(degree as u8)
                    .with_no_implicit(true),
                degree as i32,
                0,
            ),
            expected
        );
    }
    assert_eq!(
        isolated(
            AtomSpec::new(Element::OG)
                .with_explicit_hydrogens(2)
                .with_no_implicit(true),
            0,
            0,
        ),
        Hybridization::Sp
    );
    assert_eq!(
        isolated(
            AtomSpec::new(Element::RA)
                .with_explicit_hydrogens(2)
                .with_no_implicit(true),
            0,
            0,
        ),
        Hybridization::Sp2
    );
}

#[test]
fn zero_and_all_dative_variants_adjust_only_the_source_or_donor_endpoint() {
    for order in [
        BondOrder::Dative,
        BondOrder::DativeOne,
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
    ] {
        let donor = topology(
            vec![
                AtomSpec::new(Element::H)
                    .with_explicit_hydrogens(1)
                    .with_no_implicit(true),
                AtomSpec::new(Element::H).with_no_implicit(true),
            ],
            vec![bond(0, 1, order)],
        );
        assert_eq!(
            assign_hybridization(&donor, &valence(&[0, 0], &[0, 0]))
                .unwrap()
                .values[0],
            Hybridization::S
        );

        let acceptor = topology(
            vec![
                AtomSpec::new(Element::H)
                    .with_explicit_hydrogens(1)
                    .with_no_implicit(true),
                AtomSpec::new(Element::H).with_no_implicit(true),
            ],
            vec![bond(1, 0, order)],
        );
        assert_eq!(
            assign_hybridization(&acceptor, &valence(&[0, 0], &[0, 0]))
                .unwrap()
                .values[0],
            Hybridization::Sp
        );
    }

    for (order, expected) in [
        (BondOrder::Zero, Hybridization::S),
        (BondOrder::Single, Hybridization::Sp),
    ] {
        let graph = topology(
            vec![
                AtomSpec::new(Element::H)
                    .with_explicit_hydrogens(1)
                    .with_no_implicit(true),
                AtomSpec::new(Element::H).with_no_implicit(true),
            ],
            vec![bond(0, 1, order)],
        );
        assert_eq!(
            assign_hybridization(&graph, &valence(&[0, 0], &[0, 0]))
                .unwrap()
                .values[0],
            expected
        );
    }
}

#[test]
fn explicit_and_implicit_hydrogens_charge_radicals_and_signed_division_are_exact() {
    assert_eq!(
        isolated(
            AtomSpec::new(Element::C)
                .with_explicit_hydrogens(4)
                .with_no_implicit(true),
            4,
            0,
        ),
        Hybridization::Sp3
    );
    assert_eq!(
        isolated(AtomSpec::new(Element::C), 0, 4),
        Hybridization::Sp3
    );

    for (radicals, expected) in [
        (0, Hybridization::Sp),
        (1, Hybridization::Sp),
        (2, Hybridization::Sp2),
        (3, Hybridization::Sp2),
    ] {
        assert_eq!(
            isolated(
                AtomSpec::new(Element::C)
                    .with_no_implicit(true)
                    .with_radical_electrons(radicals),
                0,
                0,
            ),
            expected
        );
    }

    assert_eq!(
        isolated(
            AtomSpec::new(Element::C)
                .with_formal_charge(1)
                .with_no_implicit(true),
            4,
            0,
        ),
        Hybridization::S
    );
    assert_eq!(
        isolated(
            AtomSpec::new(Element::N)
                .with_explicit_hydrogens(3)
                .with_no_implicit(true),
            3,
            0,
        ),
        Hybridization::Sp3
    );
}

#[test]
fn invalid_inputs_return_exact_structured_errors_without_mutation() {
    let graph = topology(
        vec![
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Other),
            AtomSpec::new(Element::O),
        ],
        vec![bond(0, 1, BondOrder::Single)],
    );
    let snapshot = graph.clone();
    assert!(matches!(
        assign_hybridization(&graph, &valence(&[1], &[3, 1])),
        Err(HybridizationError::ValenceAssignmentLength {
            field: "explicit_valence",
            actual: 1,
            expected: 2,
        })
    ));
    assert!(matches!(
        assign_hybridization(&graph, &valence(&[1, 1], &[3])),
        Err(HybridizationError::ValenceAssignmentLength {
            field: "implicit_hydrogens",
            actual: 1,
            expected: 2,
        })
    ));
    assert!(matches!(
        assign_hybridization(&graph, &valence(&[-1, 1], &[0, 0])),
        Err(HybridizationError::InvalidValenceRow {
            atom,
            field: "explicit_valence",
            value: -1,
        }) if atom == AtomId::new(0)
    ));
    assert!(matches!(
        assign_hybridization(&graph, &valence(&[1, 1], &[0, -1])),
        Err(HybridizationError::InvalidValenceRow {
            atom,
            field: "implicit_hydrogens",
            value: -1,
        }) if atom == AtomId::new(1)
    ));
    assert_eq!(graph, snapshot);

    let overflow_graph = topology(vec![AtomSpec::new(Element::C)], vec![]);
    assert!(matches!(
        assign_hybridization(&overflow_graph, &valence(&[i32::MAX], &[1])),
        Err(HybridizationError::IntegerOverflow {
            atom,
            field: "total valence",
        }) if atom == AtomId::new(0)
    ));

    let mut malformed = graph.clone();
    malformed.atoms[0] = Atom::from_spec(AtomId::new(9), AtomSpec::new(Element::C));
    assert!(matches!(
        assign_hybridization(&malformed, &valence(&[1, 1], &[3, 1])),
        Err(HybridizationError::InvalidTopology(_))
    ));
}
