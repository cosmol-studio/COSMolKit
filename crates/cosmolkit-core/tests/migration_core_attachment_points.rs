#![cfg(feature = "op-contracts-strict")]

use cosmolkit_core::{
    __migration_hydrogens::{attachment_query_rows, place_terminal_attachment_coordinates},
    expand_attachment_points, rdkit_rb0,
};
use cosmolkit_model::{
    Atom, AtomId, AtomQueryPredicate, AtomSpec, Bond, BondId, BondOrder, BondQueryPredicate,
    BondSpec, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension, Element,
    Hybridization, QueryAtom, QueryNode, QueryStateError, QueryStateRef, TopologyBlock,
    TopologyMapping, remap_query_rows, remap_query_rows_with_appended,
};

#[test]
fn attachment_query_transport_retains_explicit_parent_and_new_null_query() {
    let old = attachment_topology(0, Hybridization::Unspecified);
    let new = attachment_topology(1, Hybridization::Unspecified);
    let old_atom = QueryAtom::from_parts(
        old.atoms[0].clone(),
        QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Any)),
        ]),
    );
    let state =
        QueryStateRef::try_for_topology(std::slice::from_ref(&old_atom), &[], &old).unwrap();
    let mapping = TopologyMapping::with_appended(1, 0, 1, 1);
    let (new_atom, new_bond) = attachment_query_rows(&new.atoms[1], &new.bonds[0], true);
    assert!(!new_atom.predicate_is_carrier_derived());
    assert_eq!(
        new_atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
    assert!(new_bond.predicate_is_carrier_derived());
    assert_eq!(
        new_bond.predicate(),
        &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
    );
    let (atoms, bonds) =
        remap_query_rows_with_appended(state, &new, &mapping, &[new_atom], &[new_bond]).unwrap();
    assert_eq!(atoms[0].predicate(), old_atom.predicate());
    assert!(!atoms[0].predicate_is_carrier_derived());
    assert_eq!(atoms[0].atom(), &new.atoms[0]);
    assert_eq!(atoms[1].atom(), &new.atoms[1]);
    assert_eq!(bonds[0].bond(), &new.bonds[0]);
    assert!(
        QueryStateRef::try_for_topology(&atoms, &bonds, &new)
            .unwrap()
            .atom_has_query(AtomId::new(1))
    );
    assert_eq!(old.atoms.len(), 1);
    assert_eq!(old_atom.atom(), &old.atoms[0]);
}

#[test]
fn attachment_query_transport_ordinary_rows_and_explicit_append_validation() {
    let old = attachment_topology(0, Hybridization::Unspecified);
    let mut new = attachment_topology(1, Hybridization::Unspecified);
    new.atoms[1].set_prop("_fromAttchpt", "2").unwrap();
    let old_atom = QueryAtom::from_carrier_parts(
        old.atoms[0].clone(),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
    );
    let state =
        QueryStateRef::try_for_topology(std::slice::from_ref(&old_atom), &[], &old).unwrap();
    let mapping = TopologyMapping::with_appended(1, 0, 1, 1);
    assert!(matches!(
        remap_query_rows(state, &new, &mapping),
        Err(QueryStateError::AppendedRow {
            entity: "atom",
            position: 1
        })
    ));
    let (new_atom, new_bond) = attachment_query_rows(&new.atoms[1], &new.bonds[0], false);
    assert!(new_atom.predicate_is_carrier_derived());
    assert_eq!(
        new_atom.predicate(),
        &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0))
    );
    let (atoms, bonds) = remap_query_rows_with_appended(
        state,
        &new,
        &mapping,
        &[new_atom.clone()],
        &[new_bond.clone()],
    )
    .unwrap();
    assert!(atoms.iter().all(QueryAtom::predicate_is_carrier_derived));
    assert_eq!(atoms[1].atom().prop("_fromAttchpt"), Some("2"));
    assert!(bonds[0].predicate_is_carrier_derived());
    assert!(matches!(
        remap_query_rows_with_appended(
            state,
            &new,
            &mapping,
            &[new_atom.clone(), new_atom],
            &[new_bond.clone()]
        ),
        Err(QueryStateError::AppendedRowCount {
            entity: "atom",
            actual: 2,
            expected: 1
        })
    ));
    let bad_bond = Bond::from_spec(
        BondId::new(0),
        BondSpec::new(AtomId::new(1), AtomId::new(0), BondOrder::Single),
    );
    let (_, bad_query_bond) = attachment_query_rows(&new.atoms[1], &bad_bond, false);
    assert!(matches!(
        remap_query_rows_with_appended(state, &new, &mapping, &atoms[1..], &[bad_query_bond]),
        Err(QueryStateError::BondEndpoints { position: 0, .. })
    ));
    assert_eq!(old.atoms[0].prop("_fromAttchpt"), None);
}

fn attachment_source(value: Option<i32>) -> TopologyBlock {
    let mut topology = attachment_topology(0, Hybridization::Unspecified);
    if let Some(value) = value {
        topology.atoms[0]
            .set_prop("molAttachPoint", value.to_string())
            .unwrap();
    }
    topology
}

#[test]
fn attachment_expansion_values_options_and_query_origins() {
    for value in [None, Some(0), Some(1), Some(2), Some(-1), Some(3)] {
        for add_as_queries in [false, true] {
            for add_coords in [false, true] {
                let source = attachment_source(value);
                let source_snapshot = source.clone();
                let original = CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(7, vec![[1.0, 2.0]])],
                    conformers_3d: vec![
                        Conformer3D::new(10, vec![[1.0, 2.0, -0.0]], false),
                        Conformer3D::new(11, vec![[2.0, 3.0, 0.5]], true),
                    ],
                    source_coordinate_dim: Some(CoordinateDimension::ThreeD),
                };
                let original_snapshot = original.clone();
                let result = expand_attachment_points(
                    source,
                    original,
                    None,
                    &[value],
                    add_as_queries,
                    add_coords,
                )
                .unwrap();
                let labels: &[&str] = match value {
                    Some(1) => &["1"],
                    Some(2) => &["2"],
                    Some(-1) => &["1", "2"],
                    _ => &[],
                };
                assert_eq!(result.topology.atoms.len(), 1 + labels.len());
                assert_eq!(result.topology.bonds.len(), labels.len());
                assert_eq!(result.appended_valence.len(), labels.len());
                assert_eq!(result.query_rows, None);
                assert_eq!(result.mapping.atoms().new_to_old()[0], Some(AtomId::new(0)));
                assert_eq!(result.mapping.atoms().new_to_old().len(), 1 + labels.len());
                assert_eq!(
                    result.warnings.len(),
                    usize::from(matches!(value, Some(0 | 3)))
                );
                assert_eq!(
                    result.topology.atoms[0].prop("molAttachPoint"),
                    if labels.is_empty() {
                        source_snapshot.atoms[0].prop("molAttachPoint")
                    } else {
                        None
                    }
                );
                for (offset, &label) in labels.iter().enumerate() {
                    let atom = &result.topology.atoms[offset + 1];
                    assert_eq!(atom.prop("_fromAttchpt"), Some(label));
                    assert_eq!(atom.atomic_number(), 0);
                    assert_eq!(result.topology.bonds[offset].begin(), AtomId::new(0));
                    assert_eq!(result.topology.bonds[offset].end(), AtomId::new(offset + 1));
                    assert_eq!(
                        result.appended_query_atoms[offset].predicate_is_carrier_derived(),
                        !add_as_queries
                    );
                    assert!(result.appended_query_bonds[offset].predicate_is_carrier_derived());
                    assert_eq!(result.appended_valence[offset].0, AtomId::new(offset + 1));
                    assert_eq!(result.appended_valence[offset].1, 1);
                }
                result.topology.validate().unwrap();
                result
                    .coordinates
                    .validate_for_atom_count(1 + labels.len())
                    .unwrap();
                assert_eq!(
                    result.coordinates.source_coordinate_dim,
                    original_snapshot.source_coordinate_dim
                );
                if !add_coords {
                    for conformer in &result.coordinates.conformers_3d {
                        assert!(
                            conformer.coordinates()[1..]
                                .iter()
                                .all(|row| *row == [0.0; 3])
                        );
                    }
                }
                assert_eq!(source_snapshot.atoms.len(), 1);
                assert_eq!(original_snapshot.conformers_3d[0].coordinates().len(), 1);
            }
        }
    }
}

#[test]
fn attachment_expansion_ck_coord_001_and_sequential_second_row() {
    let coordinates = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(10, vec![[1.0, 2.0, -0.0]], false),
            Conformer3D::new(11, vec![[2.0, 3.0, 0.5]], true),
        ],
        ..CoordinateBlock::default()
    };
    let result = expand_attachment_points(
        attachment_source(Some(-1)),
        coordinates.clone(),
        None,
        &[Some(-1)],
        true,
        true,
    )
    .unwrap();
    assert_eq!(
        result.coordinates.conformers_3d[0].coordinates()[0],
        [1.0, 2.0, -0.0]
    );
    assert_eq!(
        result.coordinates.conformers_3d[0].coordinates()[1],
        [2.0, 2.0, 0.0]
    );
    assert_eq!(
        result.coordinates.conformers_3d[0].coordinates()[2],
        [0.0, 2.0, -0.0]
    );
    let distance = rdkit_rb0(1) + rdkit_rb0(6);
    assert_eq!(
        result.coordinates.conformers_3d[1].coordinates()[1],
        [2.0, 3.0, 0.5 + distance]
    );
    let single = expand_attachment_points(
        attachment_source(Some(1)),
        coordinates,
        None,
        &[Some(1)],
        true,
        true,
    )
    .unwrap();
    assert_eq!(
        single.coordinates.conformers_3d[0].coordinates()[1],
        result.coordinates.conformers_3d[0].coordinates()[1]
    );
    assert_eq!(
        single.coordinates.conformers_3d[1].coordinates()[1],
        result.coordinates.conformers_3d[1].coordinates()[1]
    );
}

#[test]
fn attachment_expansion_preserves_explicit_query_and_rejects_invalid_inputs_atomically() {
    let source = attachment_source(Some(1));
    let old_query = QueryAtom::from_parts(
        source.atoms[0].clone(),
        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
    );
    let old_rows = [old_query.clone()];
    let state = QueryStateRef::try_for_topology(&old_rows, &[], &source).unwrap();
    let result = expand_attachment_points(
        source.clone(),
        CoordinateBlock::default(),
        Some(state),
        &[Some(1)],
        true,
        true,
    )
    .unwrap();
    let (atoms, bonds) = result.query_rows.unwrap();
    assert_eq!(atoms[0].predicate(), old_query.predicate());
    assert!(!atoms[0].predicate_is_carrier_derived());
    assert_eq!(atoms[0].atom().prop("molAttachPoint"), None);
    assert_eq!(
        atoms[1].predicate(),
        &QueryNode::predicate(AtomQueryPredicate::Any)
    );
    assert!(!atoms[1].predicate_is_carrier_derived());
    assert!(bonds[0].predicate_is_carrier_derived());
    assert_eq!(source.atoms[0].prop("molAttachPoint"), Some("1"));
    assert!(
        expand_attachment_points(
            source.clone(),
            CoordinateBlock::default(),
            None,
            &[],
            true,
            true,
        )
        .is_err()
    );
    let invalid_coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(1, vec![])],
        ..CoordinateBlock::default()
    };
    assert!(
        expand_attachment_points(
            source.clone(),
            invalid_coordinates,
            None,
            &[Some(1)],
            true,
            true,
        )
        .is_err()
    );
    assert_eq!(source.atoms[0].prop("molAttachPoint"), Some("1"));
}

fn attachment_topology(dummy_count: usize, hybridization: Hybridization) -> TopologyBlock {
    let mut atoms = vec![Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_hybridization(hybridization),
    )];
    let mut bonds = Vec::new();
    for index in 0..dummy_count {
        atoms.push(Atom::from_spec(
            AtomId::new(index + 1),
            AtomSpec::new(Element::DUMMY),
        ));
        bonds.push(Bond::from_spec(
            BondId::new(index),
            BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single),
        ));
    }
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

fn star_attachment_topology(neighbor_count: usize, hybridization: Hybridization) -> TopologyBlock {
    let mut atoms = vec![Atom::from_spec(
        AtomId::new(0),
        AtomSpec::new(Element::C).with_hybridization(hybridization),
    )];
    let mut bonds = Vec::new();
    for index in 0..neighbor_count {
        atoms.push(Atom::from_spec(
            AtomId::new(index + 1),
            AtomSpec::new(Element::C),
        ));
        bonds.push(Bond::from_spec(
            BondId::new(index),
            BondSpec::new(AtomId::new(0), AtomId::new(index + 1), BondOrder::Single),
        ));
    }
    let dummy = neighbor_count + 1;
    atoms.push(Atom::from_spec(
        AtomId::new(dummy),
        AtomSpec::new(Element::DUMMY),
    ));
    bonds.push(Bond::from_spec(
        BondId::new(neighbor_count),
        BondSpec::new(AtomId::new(0), AtomId::new(dummy), BondOrder::Single),
    ));
    TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
}

#[test]
fn attachment_coordinates_ck_coord_001_mixed_flags_preserve_xyz_and_order() {
    // RDKit 2026.03.1 AddHs.cpp::setTerminalAtomCoords case 1 uses a shared
    // dirVect and yields order-dependent rows. CK-COORD-001 deliberately
    // isolates each conformer's direction; raw RDKit rows are preserved in
    // IO-mol_post.md and are not treated as CK expectations here.
    let topology = attachment_topology(1, Hybridization::Unspecified);
    let original = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(3, vec![[8.0, 9.0], [0.0, 0.0]])],
        conformers_3d: vec![
            Conformer3D::new(10, vec![[1.0, 2.0, -0.0], [0.0; 3]], false),
            Conformer3D::new(11, vec![[2.0, 2.0, 0.0], [0.0; 3]], true),
            Conformer3D::new(12, vec![[4.0, 5.0, 7.0], [0.0; 3]], false),
        ],
        source_coordinate_dim: Some(CoordinateDimension::ThreeD),
    };
    let snapshot = original.clone();
    let result = place_terminal_attachment_coordinates(
        &topology,
        original,
        AtomId::new(1),
        AtomId::new(0),
        BondId::new(0),
    )
    .unwrap();
    let length = rdkit_rb0(1) + rdkit_rb0(6);
    assert_eq!(
        result.conformers_2d[0].coordinates(),
        &[[8.0, 9.0], [9.0, 9.0]]
    );
    assert_eq!(
        result.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, -0.0], [2.0, 2.0, 0.0]]
    );
    assert_eq!(
        result.conformers_3d[1].coordinates(),
        &[[2.0, 2.0, 0.0], [2.0, 2.0, length]]
    );
    assert_eq!(
        result.conformers_3d[2].coordinates(),
        &[[4.0, 5.0, 7.0], [5.0, 5.0, 7.0]]
    );
    assert_eq!(result.source_coordinate_dim, snapshot.source_coordinate_dim);
    assert_eq!(
        result
            .conformers_3d
            .iter()
            .map(|c| c.id())
            .collect::<Vec<_>>(),
        vec![10, 11, 12]
    );
    assert_eq!(snapshot.conformers_3d[2].coordinates()[0][2], 7.0);
    assert_ne!(result.conformers_3d[1].coordinates()[1], [3.1, 2.0, 1.1]);
    assert_ne!(result.conformers_3d[0].coordinates()[1], [2.0, 2.0, 1.0]);

    let reversed = CoordinateBlock {
        conformers_3d: snapshot.conformers_3d.iter().cloned().rev().collect(),
        ..snapshot.clone()
    };
    let reversed = place_terminal_attachment_coordinates(
        &topology,
        reversed,
        AtomId::new(1),
        AtomId::new(0),
        BondId::new(0),
    )
    .unwrap();
    for conformer in &result.conformers_3d {
        let reversed_row = reversed
            .conformers_3d
            .iter()
            .find(|row| row.id() == conformer.id())
            .unwrap();
        assert_eq!(conformer.coordinates(), reversed_row.coordinates());
        let isolated = place_terminal_attachment_coordinates(
            &topology,
            CoordinateBlock {
                conformers_3d: vec![
                    snapshot
                        .conformers_3d
                        .iter()
                        .find(|row| row.id() == conformer.id())
                        .unwrap()
                        .clone(),
                ],
                ..Default::default()
            },
            AtomId::new(1),
            AtomId::new(0),
            BondId::new(0),
        )
        .unwrap();
        assert_eq!(
            conformer.coordinates(),
            isolated.conformers_3d[0].coordinates()
        );
    }
}

#[test]
fn attachment_coordinates_degree_three_four_and_degenerate_source_branches() {
    // AddHs.cpp cases 3 and 4: the degree-three bisector uses the two
    // normalized away vectors; planar degree four chooses the outer-angle
    // bisector, while planar 3D uses the first nonzero cross product.
    let degree_three = star_attachment_topology(2, Hybridization::Sp2);
    let coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            1,
            vec![[0.0, 0.0], [-1.0, 0.0], [0.0, -1.0], [0.0; 2]],
        )],
        conformers_3d: vec![Conformer3D::new(
            2,
            vec![[0.0; 3], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0; 3]],
            true,
        )],
        ..Default::default()
    };
    let placed = place_terminal_attachment_coordinates(
        &degree_three,
        coordinates,
        AtomId::new(3),
        AtomId::new(0),
        BondId::new(2),
    )
    .unwrap();
    let diagonal = std::f64::consts::FRAC_1_SQRT_2;
    let length = rdkit_rb0(1) + rdkit_rb0(6);
    let planar = placed.conformers_2d[0].coordinates()[3];
    assert!((planar[0] - diagonal).abs() < 1e-12);
    assert!((planar[1] - diagonal).abs() < 1e-12);
    let spatial = placed.conformers_3d[0].coordinates()[3];
    assert!((spatial[0] - diagonal * length).abs() < 1e-12);
    assert!((spatial[1] - diagonal * length).abs() < 1e-12);
    assert_eq!(spatial[2], 0.0);

    let cancel = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            3,
            vec![[0.0; 2], [-1.0, 0.0], [1.0, 0.0], [0.0; 2]],
        )],
        ..Default::default()
    };
    let cancelled = place_terminal_attachment_coordinates(
        &degree_three,
        cancel,
        AtomId::new(3),
        AtomId::new(0),
        BondId::new(2),
    )
    .unwrap();
    assert_eq!(cancelled.conformers_2d[0].coordinates()[3], [0.0, 0.0]);

    let degree_four = star_attachment_topology(3, Hybridization::Sp3);
    let four = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            4,
            vec![[0.0; 2], [-1.0, 0.0], [0.0, -1.0], [1.0, 0.0], [0.0; 2]],
        )],
        conformers_3d: vec![Conformer3D::new(
            5,
            vec![
                [0.0; 3],
                [-1.0, 0.0, 0.0],
                [0.0, -1.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0; 3],
            ],
            true,
        )],
        ..Default::default()
    };
    let placed_four = place_terminal_attachment_coordinates(
        &degree_four,
        four,
        AtomId::new(4),
        AtomId::new(0),
        BondId::new(3),
    )
    .unwrap();
    assert_eq!(placed_four.conformers_2d[0].coordinates()[4], [0.0, 1.0]);
    assert_eq!(
        placed_four.conformers_3d[0].coordinates()[4],
        [0.0, 0.0, length]
    );
}

#[test]
fn attachment_coordinates_sequential_second_dummy_uses_first_as_neighbor() {
    // MolOps.cpp::expandAttachmentPoints appends -1's label 1 before label 2.
    // The second setTerminalAtomCoords sees parent degree two, so its default
    // hybridization branch points away from the first dummy.
    let topology = attachment_topology(2, Hybridization::Unspecified);
    let input = CoordinateBlock {
        conformers_3d: vec![
            Conformer3D::new(10, vec![[1.0, 2.0, -0.0], [0.0; 3], [0.0; 3]], false),
            Conformer3D::new(11, vec![[1.0, 2.0, -0.0], [0.0; 3], [0.0; 3]], true),
        ],
        ..Default::default()
    };
    let once = place_terminal_attachment_coordinates(
        &topology,
        input,
        AtomId::new(1),
        AtomId::new(0),
        BondId::new(0),
    )
    .unwrap();
    let twice = place_terminal_attachment_coordinates(
        &topology,
        once,
        AtomId::new(2),
        AtomId::new(0),
        BondId::new(1),
    )
    .unwrap();
    let length = rdkit_rb0(1) + rdkit_rb0(6);
    assert_eq!(
        twice.conformers_3d[0].coordinates(),
        &[[1.0, 2.0, -0.0], [2.0, 2.0, 0.0], [0.0, 2.0, 0.0]]
    );
    assert_eq!(
        twice.conformers_3d[1].coordinates(),
        &[[1.0, 2.0, -0.0], [1.0, 2.0, length], [1.0, 2.0, -length]]
    );
}

#[test]
fn attachment_coordinates_no_conformer_and_invalid_rows_fail_atomically() {
    let topology = attachment_topology(1, Hybridization::Unspecified);
    assert_eq!(
        place_terminal_attachment_coordinates(
            &topology,
            CoordinateBlock::default(),
            AtomId::new(1),
            AtomId::new(0),
            BondId::new(0),
        )
        .unwrap(),
        CoordinateBlock::default(),
    );
    let short = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(7, vec![[4.0, 5.0]])],
        ..Default::default()
    };
    let snapshot = short.clone();
    assert!(
        place_terminal_attachment_coordinates(
            &topology,
            short.clone(),
            AtomId::new(1),
            AtomId::new(0),
            BondId::new(0),
        )
        .is_err()
    );
    assert_eq!(short, snapshot);
    assert!(
        place_terminal_attachment_coordinates(
            &topology,
            CoordinateBlock::default(),
            AtomId::new(0),
            AtomId::new(0),
            BondId::new(0),
        )
        .is_err()
    );
}
