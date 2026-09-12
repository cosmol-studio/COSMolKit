use cosmolkit_core::{
    CipRankError, ValenceAssignment, assign_atom_cip_ranks, refine_atom_cip_ranks_from_invariants,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, ChiralTag, Element};

fn atom(id: usize, element: Element) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(element))
}

fn atom_from_spec(id: usize, spec: AtomSpec) -> Atom {
    Atom::from_spec(AtomId::new(id), spec)
}

fn bond(id: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn valence(explicit: Vec<i32>, implicit: Vec<i32>) -> ValenceAssignment {
    ValenceAssignment {
        explicit_valence: explicit,
        implicit_hydrogens: implicit,
    }
}

fn zero_valence(atom_count: usize) -> ValenceAssignment {
    valence(vec![0; atom_count], vec![0; atom_count])
}

fn carbon_topology(atom_count: usize, bonds: Vec<Bond>) -> TopologyBlock {
    topology(
        (0..atom_count)
            .map(|index| atom(index, Element::C))
            .collect(),
        bonds,
    )
}

#[test]
fn empty_singleton_disconnected_symmetric_and_distinguished_rows_are_exact() {
    assert_eq!(
        assign_atom_cip_ranks(&TopologyBlock::default(), &zero_valence(0)),
        Ok(vec![])
    );

    let singleton = topology(vec![atom(0, Element::C)], vec![]);
    assert_eq!(
        assign_atom_cip_ranks(&singleton, &zero_valence(1)),
        Ok(vec![0])
    );

    let disconnected = carbon_topology(3, vec![]);
    assert_eq!(
        assign_atom_cip_ranks(&disconnected, &zero_valence(3)),
        Ok(vec![0, 0, 0])
    );

    let symmetric = carbon_topology(
        4,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(&symmetric, &valence(vec![1, 2, 2, 1], vec![3, 2, 2, 3])),
        Ok(vec![0, 1, 1, 0])
    );

    let distinguished = topology(
        vec![
            atom(0, Element::C),
            atom(1, Element::N),
            atom(2, Element::O),
        ],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&distinguished, &zero_valence(3)),
        Ok(vec![0, 1, 2])
    );
}

#[test]
fn rdkit_atom_ranking_source_cases_have_complete_rank_rows() {
    // RDKit molopstest.cpp: "Testing Atom Ranking", FC(Cl)(Br)C.
    let first = topology(
        vec![
            atom(0, Element::F),
            atom(1, Element::C),
            atom(2, Element::CL),
            atom(3, Element::BR),
            atom(4, Element::C),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 1, 3, BondOrder::Single),
            bond(3, 1, 4, BondOrder::Single),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(&first, &valence(vec![1, 4, 1, 1, 1], vec![0, 0, 0, 0, 3])),
        Ok(vec![2, 1, 3, 4, 0])
    );

    // RDKit molopstest.cpp: "Testing Atom Ranking", FC(Cl)(Br)C(F)(F)F.
    let second = topology(
        vec![
            atom(0, Element::F),
            atom(1, Element::C),
            atom(2, Element::CL),
            atom(3, Element::BR),
            atom(4, Element::C),
            atom(5, Element::F),
            atom(6, Element::F),
            atom(7, Element::F),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 1, 3, BondOrder::Single),
            bond(3, 1, 4, BondOrder::Single),
            bond(4, 4, 5, BondOrder::Single),
            bond(5, 4, 6, BondOrder::Single),
            bond(6, 4, 7, BondOrder::Single),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(&second, &valence(vec![1, 4, 1, 1, 4, 1, 1, 1], vec![0; 8])),
        Ok(vec![3, 1, 4, 5, 0, 2, 2, 2])
    );
}

#[test]
fn rdkit_issue_188_cases_preserve_complete_order_across_neighbor_shells() {
    // OC[C@H](C=C)C. The chiral tag is deliberately present but is not an
    // initial legacy CIP invariant; explicit and implicit H still refine ties.
    let first = topology(
        vec![
            atom(0, Element::O),
            atom(1, Element::C),
            atom_from_spec(
                2,
                AtomSpec::new(Element::C)
                    .with_chiral_tag(ChiralTag::TetrahedralCcw)
                    .with_explicit_hydrogens(1),
            ),
            atom(3, Element::C),
            atom(4, Element::C),
            atom(5, Element::C),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 2, 3, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Double),
            bond(4, 2, 5, BondOrder::Single),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(
            &first,
            &valence(vec![1, 2, 4, 3, 2, 1], vec![1, 2, 0, 1, 2, 3])
        ),
        Ok(vec![5, 4, 3, 2, 1, 0])
    );

    // CC(=N\N)/C=N/N. Bond directions are parser/stereo state and do not
    // participate in this detached primitive, so only topology is reproduced.
    let second = topology(
        vec![
            atom(0, Element::C),
            atom(1, Element::C),
            atom(2, Element::N),
            atom(3, Element::N),
            atom(4, Element::C),
            atom(5, Element::N),
            atom(6, Element::N),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Double),
            bond(2, 2, 3, BondOrder::Single),
            bond(3, 1, 4, BondOrder::Single),
            bond(4, 4, 5, BondOrder::Double),
            bond(5, 5, 6, BondOrder::Single),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(
            &second,
            &valence(vec![1, 4, 3, 1, 3, 3, 1], vec![3, 0, 0, 2, 1, 0, 2])
        ),
        Ok(vec![0, 2, 6, 4, 1, 5, 3])
    );
}

#[test]
fn provided_invariant_seeding_refines_ties_and_validates_row_values() {
    let chain = carbon_topology(
        3,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
        ],
    );
    let assignment = valence(vec![1, 2, 1], vec![3, 2, 3]);
    assert_eq!(
        refine_atom_cip_ranks_from_invariants(&chain, &assignment, &[9, 9, 9]),
        Ok(vec![0, 1, 0])
    );
    assert_eq!(
        refine_atom_cip_ranks_from_invariants(&chain, &assignment, &[7, 9, 7]),
        Ok(vec![0, 1, 0])
    );
    assert_eq!(
        refine_atom_cip_ranks_from_invariants(&chain, &assignment, &[1, 2]),
        Err(CipRankError::InvariantCount {
            actual: 2,
            atom_count: 3,
        })
    );
    for (atom_index, value) in [(0, i64::from(i32::MAX) + 1), (1, i64::from(i32::MIN) - 1)] {
        let mut invariants = vec![0; 3];
        invariants[atom_index] = value;
        assert_eq!(
            refine_atom_cip_ranks_from_invariants(&chain, &assignment, &invariants),
            Err(CipRankError::InvariantOutOfRange {
                atom: AtomId::new(atom_index),
                value,
            })
        );
    }
}

#[test]
fn isotope_and_atom_map_bit_packing_matches_source_wrap_rules() {
    let isotopes = topology(
        vec![
            atom_from_spec(0, AtomSpec::new(Element::C).with_isotope(11)),
            atom(1, Element::C),
            atom_from_spec(2, AtomSpec::new(Element::C).with_isotope(12)),
            atom_from_spec(3, AtomSpec::new(Element::C).with_isotope(1036)),
            atom_from_spec(4, AtomSpec::new(Element::C).with_isotope(13)),
        ],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&isotopes, &zero_valence(5)),
        Ok(vec![0, 1, 2, 2, 3])
    );

    let maps = topology(
        vec![
            atom(0, Element::C),
            atom_from_spec(1, AtomSpec::new(Element::C).with_atom_map(1023)),
            atom_from_spec(2, AtomSpec::new(Element::C).with_atom_map(0)),
            atom_from_spec(3, AtomSpec::new(Element::C).with_atom_map(5)),
            atom_from_spec(4, AtomSpec::new(Element::C).with_atom_map(i32::MAX as u32)),
        ],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&maps, &zero_valence(5)),
        Ok(vec![0, 0, 1, 2, 0])
    );

    let out_of_range = topology(
        vec![atom_from_spec(
            0,
            AtomSpec::new(Element::C).with_atom_map(i32::MAX as u32 + 1),
        )],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&out_of_range, &zero_valence(1)),
        Err(CipRankError::AtomMapOutOfRange {
            atom: AtomId::new(0),
            map_number: i32::MAX as u32 + 1,
        })
    );
}

#[test]
fn charge_aromaticity_and_stereo_are_excluded_but_hydrogens_refine_ties() {
    let excluded = topology(
        vec![
            atom_from_spec(0, AtomSpec::new(Element::C).with_formal_charge(-1)),
            atom_from_spec(1, AtomSpec::new(Element::C).with_aromatic(true)),
            atom_from_spec(
                2,
                AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
            ),
        ],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&excluded, &zero_valence(3)),
        Ok(vec![0, 0, 0])
    );

    let implicit = carbon_topology(2, vec![]);
    assert_eq!(
        assign_atom_cip_ranks(&implicit, &valence(vec![0, 0], vec![0, 1])),
        Ok(vec![0, 1])
    );
    let explicit = topology(
        vec![
            atom(0, Element::C),
            atom_from_spec(1, AtomSpec::new(Element::C).with_explicit_hydrogens(1)),
        ],
        vec![],
    );
    assert_eq!(
        assign_atom_cip_ranks(&explicit, &zero_valence(2)),
        Ok(vec![0, 1])
    );
}

#[test]
fn every_supported_bond_order_reaches_the_exact_multiplicity_table() {
    for order in [
        BondOrder::Unspecified,
        BondOrder::Ionic,
        BondOrder::Zero,
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Triple,
        BondOrder::Quadruple,
        BondOrder::Quintuple,
        BondOrder::Hextuple,
        BondOrder::OneAndHalf,
        BondOrder::TwoAndHalf,
        BondOrder::ThreeAndHalf,
        BondOrder::FourAndHalf,
        BondOrder::FiveAndHalf,
        BondOrder::Aromatic,
        BondOrder::DativeOne,
        BondOrder::Dative,
        BondOrder::Hydrogen,
    ] {
        let pair = carbon_topology(2, vec![bond(0, 0, 1, order)]);
        assert_eq!(
            assign_atom_cip_ranks(&pair, &zero_valence(2)),
            Ok(vec![0, 0]),
            "supported order {order:?}"
        );
    }
}

#[test]
fn every_unsupported_bond_order_reports_exact_bond_and_order() {
    for order in [
        BondOrder::DativeLeft,
        BondOrder::DativeRight,
        BondOrder::ThreeCenter,
        BondOrder::Other,
    ] {
        let pair = carbon_topology(2, vec![bond(0, 0, 1, order)]);
        assert_eq!(
            assign_atom_cip_ranks(&pair, &zero_valence(2)),
            Err(CipRankError::UnsupportedBondOrder {
                bond: BondId::new(0),
                order,
            })
        );
    }
}

fn phosphorus_case(degree: usize, reverse_double: bool) -> TopologyBlock {
    assert!((2..=4).contains(&degree));
    let mut atoms = vec![
        atom(0, Element::C),
        atom(1, Element::P),
        atom(2, Element::C),
    ];
    for index in 3..=degree {
        atoms.push(atom(
            index,
            if index == 3 { Element::F } else { Element::CL },
        ));
    }
    let mut bonds = vec![
        if reverse_double {
            bond(0, 1, 0, BondOrder::Double)
        } else {
            bond(0, 0, 1, BondOrder::Double)
        },
        bond(1, 2, 1, BondOrder::Single),
    ];
    for index in 3..=degree {
        bonds.push(bond(index - 1, 1, index, BondOrder::Single));
    }
    topology(atoms, bonds)
}

#[test]
fn phosphorus_degree_two_three_four_boundary_and_orientation_are_exact() {
    for degree in 2..=4 {
        let forward = phosphorus_case(degree, false);
        let reversed = phosphorus_case(degree, true);
        let forward_ranks = assign_atom_cip_ranks(&forward, &zero_valence(degree + 1)).unwrap();
        let reversed_ranks = assign_atom_cip_ranks(&reversed, &zero_valence(degree + 1)).unwrap();
        assert_eq!(forward_ranks, reversed_ranks, "degree {degree}");
        if degree == 2 {
            assert!(forward_ranks[0] > forward_ranks[2]);
        } else {
            assert!(forward_ranks[0] < forward_ranks[2]);
        }
    }
}

#[test]
fn rdkit_issue_3009911_aromatic_multiplicity_has_complete_rank_row() {
    // F[C@](O)(c1ccccc1)C(=C)CO, directly constructed from the pinned source
    // case. Aromatic flags are retained on both atom and bond rows.
    let aromatic_carbon = |id| atom_from_spec(id, AtomSpec::new(Element::C).with_aromatic(true));
    let molecule = topology(
        vec![
            atom(0, Element::F),
            atom_from_spec(
                1,
                AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
            ),
            atom(2, Element::O),
            aromatic_carbon(3),
            aromatic_carbon(4),
            aromatic_carbon(5),
            aromatic_carbon(6),
            aromatic_carbon(7),
            aromatic_carbon(8),
            atom(9, Element::C),
            atom(10, Element::C),
            atom(11, Element::C),
            atom(12, Element::O),
        ],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 1, 3, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Aromatic),
            bond(4, 4, 5, BondOrder::Aromatic),
            bond(5, 5, 6, BondOrder::Aromatic),
            bond(6, 6, 7, BondOrder::Aromatic),
            bond(7, 7, 8, BondOrder::Aromatic),
            bond(8, 1, 9, BondOrder::Single),
            bond(9, 9, 10, BondOrder::Double),
            bond(10, 9, 11, BondOrder::Single),
            bond(11, 11, 12, BondOrder::Single),
            bond(12, 8, 3, BondOrder::Aromatic),
        ],
    );
    assert_eq!(
        assign_atom_cip_ranks(
            &molecule,
            &valence(
                vec![1, 4, 1, 3, 3, 3, 3, 3, 3, 4, 2, 2, 1],
                vec![0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 2, 2, 1],
            )
        ),
        Ok(vec![10, 7, 9, 4, 3, 2, 1, 2, 3, 5, 0, 6, 8])
    );
}

#[test]
fn topology_and_valence_errors_are_rejected_before_ranking() {
    let valid = carbon_topology(2, vec![bond(0, 0, 1, BondOrder::Single)]);
    assert_eq!(
        assign_atom_cip_ranks(&valid, &valence(vec![0], vec![0, 0])),
        Err(CipRankError::ValenceRowCount {
            field: "explicit_valence",
            actual: 1,
            atom_count: 2,
        })
    );
    assert_eq!(
        assign_atom_cip_ranks(&valid, &valence(vec![0, 0], vec![0])),
        Err(CipRankError::ValenceRowCount {
            field: "implicit_hydrogens",
            actual: 1,
            atom_count: 2,
        })
    );
    assert_eq!(
        assign_atom_cip_ranks(&valid, &valence(vec![0, 0], vec![0, -2])),
        Err(CipRankError::NegativeImplicitHydrogen {
            atom: AtomId::new(1),
            value: -2,
        })
    );

    let mut invalid = valid.clone();
    invalid.atoms[0] = invalid.atoms[0].clone().with_id(AtomId::new(1));
    assert_eq!(
        assign_atom_cip_ranks(&invalid, &zero_valence(2)),
        Err(CipRankError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );
}

#[test]
fn sixteen_neighbors_fail_with_the_source_stride_boundary_fields() {
    let mut bonds = Vec::new();
    for leaf in 1..=16 {
        bonds.push(bond(leaf - 1, 0, leaf, BondOrder::Single));
    }
    let crowded = carbon_topology(17, bonds);
    assert_eq!(
        assign_atom_cip_ranks(&crowded, &zero_valence(17)),
        Err(CipRankError::TooManyNeighbors {
            atom: AtomId::new(0),
            degree: 16,
            maximum_supported: 15,
        })
    );
}

#[test]
fn repeated_calls_are_stable_and_leave_detached_inputs_unchanged() {
    let molecule = carbon_topology(
        5,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 1, 2, BondOrder::Double),
            bond(2, 2, 3, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Single),
        ],
    );
    let assignment = valence(vec![1, 3, 3, 2, 1], vec![3, 1, 1, 2, 3]);
    let topology_before = molecule.clone();
    let valence_before = assignment.clone();
    let first = assign_atom_cip_ranks(&molecule, &assignment).unwrap();
    let second = assign_atom_cip_ranks(&molecule, &assignment).unwrap();
    assert_eq!(first, second);
    assert_eq!(molecule, topology_before);
    assert_eq!(assignment, valence_before);
}
