use std::cell::Cell;

use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomMapping, AtomSpec, Bond, BondId, BondMapping, BondOrder,
    BondSpec, BondStereo, Conformer2D, CoordinateBlock, Element, MappingValidationError,
    MoleculeProperties, SGroupAttachPoint, SGroupCState, SdfPropertyList, SdfPropertyListTarget,
    StereoGroup, StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
    TopologyBlock, TopologyMapping,
};

fn atom_id(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond_id(index: usize) -> BondId {
    BondId::new(index)
}

fn mapping(
    atom_old_to_new: Vec<Option<usize>>,
    atom_new_to_old: Vec<Option<usize>>,
    bond_old_to_new: Vec<Option<usize>>,
    bond_new_to_old: Vec<Option<usize>>,
) -> TopologyMapping {
    TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new
                .into_iter()
                .map(|id| id.map(atom_id))
                .collect(),
            new_to_old: atom_new_to_old
                .into_iter()
                .map(|id| id.map(atom_id))
                .collect(),
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new
                .into_iter()
                .map(|id| id.map(bond_id))
                .collect(),
            new_to_old: bond_new_to_old
                .into_iter()
                .map(|id| id.map(bond_id))
                .collect(),
        },
    }
}

fn assert_valid(mapping: &TopologyMapping, counts: (usize, usize, usize, usize)) {
    assert_eq!(
        mapping.validate_for_counts(counts.0, counts.1, counts.2, counts.3),
        Ok(())
    );
}

#[test]
fn identity_and_append_constructors_cover_independent_zero_and_nonzero_counts() {
    for (atoms, bonds) in [(0, 0), (3, 0), (0, 2), (3, 2)] {
        let identity = TopologyMapping::identity(atoms, bonds);
        assert_eq!(identity.atoms().old_to_new().len(), atoms);
        assert_eq!(identity.atoms().new_to_old().len(), atoms);
        assert_eq!(identity.bonds().old_to_new().len(), bonds);
        assert_eq!(identity.bonds().new_to_old().len(), bonds);
        assert_eq!(
            identity.atoms().old_to_new(),
            (0..atoms).map(|row| Some(atom_id(row))).collect::<Vec<_>>()
        );
        assert_eq!(
            identity.bonds().new_to_old(),
            (0..bonds).map(|row| Some(bond_id(row))).collect::<Vec<_>>()
        );
        assert_eq!(
            identity.retained_atom_indices(),
            (0..atoms).collect::<Vec<_>>()
        );
        assert_valid(&identity, (atoms, atoms, bonds, bonds));
    }

    for (old_atoms, old_bonds, added_atoms, added_bonds) in
        [(0, 0, 0, 0), (2, 0, 3, 0), (0, 2, 0, 3), (2, 1, 3, 2)]
    {
        let appended =
            TopologyMapping::with_appended(old_atoms, old_bonds, added_atoms, added_bonds);
        assert_eq!(
            appended.atoms().new_to_old(),
            (0..old_atoms)
                .map(|row| Some(atom_id(row)))
                .chain((0..added_atoms).map(|_| None))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            appended.bonds().new_to_old(),
            (0..old_bonds)
                .map(|row| Some(bond_id(row)))
                .chain((0..added_bonds).map(|_| None))
                .collect::<Vec<_>>()
        );
        assert_eq!(
            appended.retained_atom_indices(),
            (0..old_atoms).collect::<Vec<_>>()
        );
        assert_valid(
            &appended,
            (
                old_atoms,
                old_atoms + added_atoms,
                old_bonds,
                old_bonds + added_bonds,
            ),
        );
    }
}

#[test]
fn deletion_reorder_and_mixed_none_rows_validate_and_preserve_new_row_order() {
    let deleted = mapping(
        vec![Some(0), None, Some(1), None],
        vec![Some(0), Some(2)],
        vec![None, Some(0), None],
        vec![Some(1)],
    );
    assert_valid(&deleted, (4, 2, 3, 1));
    assert_eq!(deleted.retained_atom_indices(), vec![0, 2]);

    let reordered = mapping(
        vec![Some(1), Some(2), Some(0)],
        vec![Some(2), Some(0), Some(1)],
        vec![Some(1), Some(0)],
        vec![Some(1), Some(0)],
    );
    assert_valid(&reordered, (3, 3, 2, 2));
    assert_eq!(reordered.retained_atom_indices(), vec![2, 0, 1]);

    let mixed = mapping(
        vec![Some(1), None, Some(0)],
        vec![Some(2), Some(0), None],
        vec![None, Some(0)],
        vec![Some(1), None],
    );
    assert_valid(&mixed, (3, 3, 2, 2));
    assert_eq!(mixed.retained_atom_indices(), vec![2, 0]);
}

#[test]
fn every_atom_and_bond_direction_reports_exact_length_error_fields() {
    let cases = [
        (
            mapping(vec![], vec![Some(0)], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::Length {
                entity: "atom",
                direction: "old-to-new",
                actual: 0,
                expected: 1,
            },
        ),
        (
            mapping(vec![Some(0)], vec![], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::Length {
                entity: "atom",
                direction: "new-to-old",
                actual: 0,
                expected: 1,
            },
        ),
        (
            mapping(vec![], vec![], vec![], vec![Some(0)]),
            (0, 0, 1, 1),
            MappingValidationError::Length {
                entity: "bond",
                direction: "old-to-new",
                actual: 0,
                expected: 1,
            },
        ),
        (
            mapping(vec![], vec![], vec![Some(0)], vec![]),
            (0, 0, 1, 1),
            MappingValidationError::Length {
                entity: "bond",
                direction: "new-to-old",
                actual: 0,
                expected: 1,
            },
        ),
    ];

    for (candidate, counts, expected) in cases {
        assert_eq!(
            candidate.validate_for_counts(counts.0, counts.1, counts.2, counts.3),
            Err(expected)
        );
    }
}

#[test]
fn every_atom_and_bond_direction_reports_exact_range_error_fields() {
    let cases = [
        (
            mapping(vec![Some(2)], vec![None], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 2,
                target_count: 1,
            },
        ),
        (
            mapping(vec![None], vec![Some(2)], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::OutOfRange {
                entity: "atom",
                direction: "new-to-old",
                row: 0,
                mapped: 2,
                target_count: 1,
            },
        ),
        (
            mapping(vec![], vec![], vec![Some(3)], vec![None]),
            (0, 0, 1, 1),
            MappingValidationError::OutOfRange {
                entity: "bond",
                direction: "old-to-new",
                row: 0,
                mapped: 3,
                target_count: 1,
            },
        ),
        (
            mapping(vec![], vec![], vec![None], vec![Some(3)]),
            (0, 0, 1, 1),
            MappingValidationError::OutOfRange {
                entity: "bond",
                direction: "new-to-old",
                row: 0,
                mapped: 3,
                target_count: 1,
            },
        ),
    ];

    for (candidate, counts, expected) in cases {
        assert_eq!(
            candidate.validate_for_counts(counts.0, counts.1, counts.2, counts.3),
            Err(expected)
        );
    }
}

#[test]
fn every_atom_and_bond_direction_reports_exact_inverse_error_fields() {
    let cases = [
        (
            mapping(vec![Some(0)], vec![None], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "old-to-new",
                row: 0,
                mapped: 0,
            },
        ),
        (
            mapping(vec![None], vec![Some(0)], vec![], vec![]),
            (1, 1, 0, 0),
            MappingValidationError::InverseMismatch {
                entity: "atom",
                direction: "new-to-old",
                row: 0,
                mapped: 0,
            },
        ),
        (
            mapping(vec![], vec![], vec![Some(0)], vec![None]),
            (0, 0, 1, 1),
            MappingValidationError::InverseMismatch {
                entity: "bond",
                direction: "old-to-new",
                row: 0,
                mapped: 0,
            },
        ),
        (
            mapping(vec![], vec![], vec![None], vec![Some(0)]),
            (0, 0, 1, 1),
            MappingValidationError::InverseMismatch {
                entity: "bond",
                direction: "new-to-old",
                row: 0,
                mapped: 0,
            },
        ),
    ];

    for (candidate, counts, expected) in cases {
        assert_eq!(
            candidate.validate_for_counts(counts.0, counts.1, counts.2, counts.3),
            Err(expected)
        );
    }
}

#[test]
fn duplicate_targets_missing_reverse_links_and_error_precedence_are_deterministic() {
    let duplicate = mapping(vec![Some(0), Some(0)], vec![Some(0), None], vec![], vec![]);
    assert_eq!(
        duplicate.validate_for_counts(2, 2, 0, 0),
        Err(MappingValidationError::InverseMismatch {
            entity: "atom",
            direction: "old-to-new",
            row: 1,
            mapped: 0,
        })
    );

    let atom_before_bond = mapping(vec![], vec![Some(0)], vec![], vec![Some(9)]);
    assert_eq!(
        atom_before_bond.validate_for_counts(1, 1, 1, 1),
        Err(MappingValidationError::Length {
            entity: "atom",
            direction: "old-to-new",
            actual: 0,
            expected: 1,
        })
    );

    let direction_before_row = mapping(vec![Some(0), Some(3)], vec![None, Some(4)], vec![], vec![]);
    assert_eq!(
        direction_before_row.validate_for_counts(2, 2, 0, 0),
        Err(MappingValidationError::InverseMismatch {
            entity: "atom",
            direction: "old-to-new",
            row: 0,
            mapped: 0,
        })
    );
}

fn guarded_project(
    mapping: &TopologyMapping,
    counts: (usize, usize, usize, usize),
    coordinates: &mut CoordinateBlock,
    properties: &mut MoleculeProperties,
    projection_ran: &Cell<bool>,
) -> Result<(), MappingValidationError> {
    mapping.validate_for_counts(counts.0, counts.1, counts.2, counts.3)?;
    projection_ran.set(true);
    coordinates.remap_topology(&mapping.retained_atom_indices());
    properties.remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
    Ok(())
}

#[test]
fn validation_guard_rejects_invalid_some_before_dangerous_projection() {
    let invalid = mapping(
        vec![Some(0), None, None],
        vec![Some(0), Some(99)],
        vec![Some(0)],
        vec![Some(0)],
    );
    let original_coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            0,
            vec![[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]],
        )],
        ..Default::default()
    };
    let original_properties = MoleculeProperties::default()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atoms",
            vec![Some("a0".into()), None, Some("a2".into())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bonds",
            vec![Some("b0".into())],
        ));
    let mut coordinates = original_coordinates.clone();
    let mut properties = original_properties.clone();
    let projection_ran = Cell::new(false);

    assert_eq!(
        guarded_project(
            &invalid,
            (3, 2, 1, 1),
            &mut coordinates,
            &mut properties,
            &projection_ran,
        ),
        Err(MappingValidationError::OutOfRange {
            entity: "atom",
            direction: "new-to-old",
            row: 1,
            mapped: 99,
            target_count: 3,
        })
    );
    assert!(!projection_ran.get());
    assert_eq!(coordinates, original_coordinates);
    assert_eq!(properties, original_properties);
}

#[test]
fn valid_reorder_and_deletion_project_coordinates_and_properties_without_conflating_none() {
    let valid = mapping(
        vec![Some(1), None, Some(0)],
        vec![Some(2), Some(0)],
        vec![None, Some(0)],
        vec![Some(1)],
    );
    let mut coordinates = CoordinateBlock {
        conformers_2d: vec![Conformer2D::new(
            4,
            vec![[10.0, 11.0], [20.0, 21.0], [30.0, 31.0]],
        )],
        ..Default::default()
    };
    let mut properties = MoleculeProperties::default()
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Atom,
            "atoms",
            vec![Some("a0".into()), None, Some("a2".into())],
        ))
        .with_sdf_property_list(SdfPropertyList::new(
            SdfPropertyListTarget::Bond,
            "bonds",
            vec![Some("b0".into()), Some("b1".into())],
        ));
    let projection_ran = Cell::new(false);

    assert_eq!(
        guarded_project(
            &valid,
            (3, 2, 2, 1),
            &mut coordinates,
            &mut properties,
            &projection_ran,
        ),
        Ok(())
    );
    assert!(projection_ran.get());
    assert_eq!(
        coordinates.conformers_2d[0].coordinates(),
        &[[30.0, 31.0], [10.0, 11.0]]
    );
    let lists = properties.sdf_property_lists();
    assert_eq!(lists[0].values(), &[Some("a2".into()), Some("a0".into())]);
    assert_eq!(lists[1].values(), &[Some("b1".into())]);

    let append = TopologyMapping::with_appended(1, 0, 1, 0);
    assert_valid(&append, (1, 2, 0, 0));
    let mut appended_properties = MoleculeProperties::default().with_sdf_property_list(
        SdfPropertyList::new(SdfPropertyListTarget::Atom, "atoms", vec![None]),
    );
    appended_properties.remap_topology(append.atoms().new_to_old(), &[]);
    assert_eq!(
        appended_properties.sdf_property_lists()[0].values(),
        &[None, None]
    );
}

fn atom(index: usize) -> Atom {
    Atom::from_spec(atom_id(index), AtomSpec::new(Element::C))
}

#[test]
fn removal_mapping_covers_bonds_sgroups_stereo_and_adjacency_references() {
    let stereo_bond = Bond::from_spec(
        bond_id(0),
        BondSpec::new(atom_id(2), atom_id(3), BondOrder::Double)
            .with_stereo_atoms(atom_id(1), atom_id(4))
            .with_stereo(BondStereo::Cis),
    );
    let removed_bond = Bond::from_spec(
        bond_id(1),
        BondSpec::new(atom_id(0), atom_id(1), BondOrder::Single),
    );
    let bonds = vec![stereo_bond, removed_bond];
    let topology = TopologyBlock {
        atoms: (0..5).map(atom).collect(),
        adjacency: AdjacencyList::from_topology(5, &bonds),
        bonds,
        substance_groups: vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![atom_id(1), atom_id(4)])
                .with_bonds(vec![bond_id(0)])
                .with_parent_atoms(vec![atom_id(2)])
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: atom_id(3),
                    leaving_atom: Some(atom_id(4)),
                    label: Some("AP".into()),
                    order: Some(1),
                }])
                .with_cstates(vec![SGroupCState {
                    bond: bond_id(0),
                    vector: [1.0, 2.0],
                }]),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_atoms(vec![atom_id(2)])
                .with_parent(SubstanceGroupId::new(0)),
        ],
        stereo_groups: vec![StereoGroup::new(
            StereoGroupKind::Or,
            vec![atom_id(1), atom_id(4)],
            vec![bond_id(0)],
        )],
    };
    assert_eq!(topology.validate(), Ok(()));

    let mut edit = topology
        .begin_batch_edit()
        .expect("valid detached topology");
    edit.remove_atom(atom_id(0)).expect("valid atom removal");
    let (topology, topology_mapping) = edit.finish().expect("valid detached edit");
    assert_valid(&topology_mapping, (5, 4, 2, 1));
    assert_eq!(
        topology_mapping.atoms().old_to_new(),
        &[
            None,
            Some(atom_id(0)),
            Some(atom_id(1)),
            Some(atom_id(2)),
            Some(atom_id(3)),
        ]
    );
    assert_eq!(
        topology_mapping.atoms().new_to_old(),
        &[
            Some(atom_id(1)),
            Some(atom_id(2)),
            Some(atom_id(3)),
            Some(atom_id(4)),
        ]
    );
    assert_eq!(
        topology_mapping.bonds().old_to_new(),
        &[Some(bond_id(0)), None]
    );
    assert_eq!(topology_mapping.bonds().new_to_old(), &[Some(bond_id(0))]);
    assert_eq!(topology.validate(), Ok(()));
    assert_eq!(
        topology.bonds[0].stereo_atoms(),
        Some([atom_id(0), atom_id(3)])
    );
    assert_eq!(
        topology.substance_groups[0].atoms(),
        &[atom_id(0), atom_id(3)]
    );
    assert_eq!(
        topology.substance_groups[1].parent(),
        Some(SubstanceGroupId::new(0))
    );
    assert_eq!(topology.stereo_groups[0].atoms(), &[atom_id(0), atom_id(3)]);
    assert_eq!(
        topology.adjacency,
        AdjacencyList::from_topology(4, &topology.bonds)
    );
}
