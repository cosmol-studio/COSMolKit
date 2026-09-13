use cosmolkit_core::{
    __migration_property_cache::{PropertyCacheError, PropertyCacheParams, assign_property_cache},
    ValenceError, ValencePhase,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind,
    SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondOrder, Element};

fn atom(atomic_number: u8) -> AtomSpec {
    AtomSpec::new(Element::from_atomic_number(atomic_number).expect("modeled test element"))
}

fn bond(begin: usize, end: usize, order: BondOrder) -> BondSpec {
    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
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
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn star(center: AtomSpec, neighbor_count: usize) -> TopologyBlock {
    let mut atoms = vec![center];
    atoms.extend((0..neighbor_count).map(|_| atom(9).with_no_implicit(true)));
    let bonds = (1..=neighbor_count)
        .map(|neighbor| bond(0, neighbor, BondOrder::Single))
        .collect();
    topology(atoms, bonds)
}

#[test]
fn defaults_empty_topology_and_result_shape_are_exact() {
    assert_eq!(
        PropertyCacheParams::default(),
        PropertyCacheParams { strict: true }
    );

    let assignment =
        assign_property_cache(&TopologyBlock::default(), &PropertyCacheParams::default()).unwrap();
    assert!(assignment.valence().explicit_valence.is_empty());
    assert!(assignment.valence().implicit_hydrogens.is_empty());
    assert!(assignment.into_valence().explicit_valence.is_empty());
}

#[test]
fn stored_atom_order_is_stable_deterministic_and_not_chemically_sorted() {
    let input = topology(vec![atom(9), atom(6), atom(8), atom(7)], vec![]);
    let first = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();
    let second = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();

    assert_eq!(first, second);
    assert_eq!(first.valence.explicit_valence, vec![0, 0, 0, 0]);
    assert_eq!(first.valence.implicit_hydrogens, vec![1, 4, 2, 3]);
    assert_eq!(first.valence.explicit_valence.len(), input.atoms.len());
    assert_eq!(first.valence.implicit_hydrogens.len(), input.atoms.len());
}

#[test]
fn calculation_preserves_atoms_bonds_adjacency_stereo_and_sgroups() {
    let mut input = topology(
        vec![atom(6), atom(6), atom(8)],
        vec![bond(0, 1, BondOrder::Single), bond(1, 2, BondOrder::Single)],
    );
    input.atoms[0].set_prop("atom-note", "kept").unwrap();
    input.bonds[1].set_prop("bond-note", "kept").unwrap();
    input.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![AtomId::new(0), AtomId::new(1)])
            .with_bonds(vec![BondId::new(0)]),
    ];
    input.stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(1)],
            vec![BondId::new(1)],
        )
        .with_id(23),
    ];
    input.validate().unwrap();
    let snapshot = input.clone();

    let assignment = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();
    assert_eq!(assignment.valence.explicit_valence, vec![1, 2, 1]);
    assert_eq!(assignment.valence.implicit_hydrogens, vec![3, 2, 1]);
    assert_eq!(input, snapshot);
}

#[test]
fn dummy_explicit_h_no_implicit_radical_and_charge_rows_are_forwarded() {
    let input = topology(
        vec![
            atom(0),
            atom(6).with_explicit_hydrogens(2).with_no_implicit(true),
            atom(6).with_radical_electrons(2),
            atom(6).with_formal_charge(1),
        ],
        vec![],
    );
    let assignment = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();

    assert_eq!(assignment.valence.explicit_valence, vec![0, 2, 0, 0]);
    assert_eq!(assignment.valence.implicit_hydrogens, vec![0, 0, 2, 3]);
}

#[test]
fn aromatic_half_order_rows_follow_source_rounding() {
    let atoms = (0..6)
        .map(|_| atom(6).with_aromatic(true))
        .collect::<Vec<_>>();
    let bonds = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
        .into_iter()
        .map(|(begin, end)| bond(begin, end, BondOrder::Aromatic).with_aromatic(true))
        .collect();
    let input = topology(atoms, bonds);
    let assignment = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();

    assert_eq!(assignment.valence.explicit_valence, vec![3; 6]);
    assert_eq!(assignment.valence.implicit_hydrogens, vec![1; 6]);
}

#[test]
fn integral_fractional_zero_and_hydrogen_bond_orders_are_forwarded() {
    for (order, expected_explicit) in [
        (BondOrder::Unspecified, 0),
        (BondOrder::Ionic, 0),
        (BondOrder::Zero, 0),
        (BondOrder::Hydrogen, 0),
        (BondOrder::Single, 1),
        (BondOrder::OneAndHalf, 2),
        (BondOrder::Double, 2),
        (BondOrder::TwoAndHalf, 3),
        (BondOrder::Triple, 3),
        (BondOrder::Quadruple, 4),
        (BondOrder::Quintuple, 5),
        (BondOrder::Hextuple, 6),
    ] {
        let input = topology(vec![atom(6), atom(6)], vec![bond(0, 1, order)]);
        let assignment =
            assign_property_cache(&input, &PropertyCacheParams { strict: false }).unwrap();
        assert_eq!(
            assignment.valence.explicit_valence,
            vec![expected_explicit; 2],
            "{order:?}"
        );
    }
}

#[test]
fn dative_direction_and_source_no_op_bond_stage_have_no_bond_payload() {
    for order in [BondOrder::Dative, BondOrder::DativeOne] {
        let input = topology(vec![atom(7), atom(8)], vec![bond(0, 1, order)]);
        let assignment = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();
        assert_eq!(assignment.valence.explicit_valence, vec![0, 1]);
        assert_eq!(assignment.valence.implicit_hydrogens, vec![3, 1]);
    }
}

#[test]
fn rejected_directional_dative_orders_preserve_bond_id_and_order() {
    for order in [BondOrder::DativeLeft, BondOrder::DativeRight] {
        let input = topology(vec![atom(7), atom(8)], vec![bond(0, 1, order)]);
        assert_eq!(
            assign_property_cache(&input, &PropertyCacheParams::default()),
            Err(PropertyCacheError::Valence(ValenceError::BadBondType {
                bond: Some(BondId::new(0)),
                order,
            }))
        );
    }
}

#[test]
fn hypervalent_anions_keep_their_source_defined_assignments() {
    for atomic_number in [15, 16, 33, 34] {
        let input = star(atom(atomic_number).with_formal_charge(-1), 5);
        let assignment = assign_property_cache(&input, &PropertyCacheParams::default()).unwrap();
        assert_eq!(assignment.valence.explicit_valence, vec![5, 1, 1, 1, 1, 1]);
        assert_eq!(
            assignment.valence.implicit_hydrogens,
            vec![
                if matches!(atomic_number, 15 | 33) {
                    1
                } else {
                    0
                },
                0,
                0,
                0,
                0,
                0
            ]
        );
    }
}

#[test]
fn strict_error_and_non_strict_result_preserve_exact_source_fields() {
    let input = star(atom(6), 5);
    assert!(matches!(
        assign_property_cache(&input, &PropertyCacheParams::default()),
        Err(PropertyCacheError::Valence(ValenceError::InvalidValence {
            atom,
            atomic_number: 6,
            formal_charge: 0,
            phase: ValencePhase::Explicit,
            calculated: Some(5),
            reason: "greater than permitted",
            ref message,
        })) if atom == AtomId::new(0)
            && message == "Explicit valence for atom # 0 C, 5, is greater than permitted"
    ));

    let assignment = assign_property_cache(&input, &PropertyCacheParams { strict: false }).unwrap();
    assert_eq!(assignment.valence.explicit_valence, vec![5, 1, 1, 1, 1, 1]);
    assert_eq!(
        assignment.valence.implicit_hydrogens,
        vec![0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn malformed_topology_error_stays_typed_and_nested() {
    let input = TopologyBlock {
        atoms: vec![Atom::from_spec(AtomId::new(1), atom(6))],
        adjacency: AdjacencyList::from_topology(1, &[]),
        ..TopologyBlock::default()
    };

    assert_eq!(
        assign_property_cache(&input, &PropertyCacheParams::default()),
        Err(PropertyCacheError::Valence(ValenceError::InvalidTopology {
            source: TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            },
        }))
    );
}

#[test]
fn later_atom_failure_returns_no_prefix_and_does_not_mutate_input() {
    let input = topology(vec![atom(8), atom(1).with_formal_charge(2)], vec![]);
    let snapshot = input.clone();
    let result = assign_property_cache(&input, &PropertyCacheParams::default());

    assert!(matches!(
        result,
        Err(PropertyCacheError::Valence(ValenceError::InvalidValence {
            atom,
            atomic_number: 1,
            formal_charge: 2,
            phase: ValencePhase::Implicit,
            calculated: Some(0),
            reason: "unreasonable hydrogen formal charge",
            ref message,
        })) if atom == AtomId::new(1)
            && message == "Unreasonable formal charge on atom # 1."
    ));
    assert_eq!(input, snapshot);
}
