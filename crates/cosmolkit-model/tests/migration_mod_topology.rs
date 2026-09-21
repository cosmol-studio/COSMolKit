use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, BondStereo, Element,
    NeighborRef, SGroupAttachPoint, SGroupCState, StereoGroup, StereoGroupKind, SubstanceGroup,
    SubstanceGroupId, SubstanceGroupKind, TopologyBlock, TopologyEditError,
    TopologyValidationError,
};

fn atom_id(index: usize) -> AtomId {
    AtomId::new(index)
}

fn bond_id(index: usize) -> BondId {
    BondId::new(index)
}

fn atom(index: usize) -> Atom {
    Atom::from_spec(atom_id(index), AtomSpec::new(Element::C))
}

fn bond(index: usize, begin: usize, end: usize) -> Bond {
    Bond::from_spec(
        bond_id(index),
        BondSpec::new(atom_id(begin), atom_id(end), BondOrder::Single),
    )
}

fn topology(atom_count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
    TopologyBlock::try_from_parts(
        (0..atom_count).map(atom).collect(),
        edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| bond(index, begin, end))
            .collect(),
        Vec::new(),
        Vec::new(),
    )
    .expect("test topology is valid")
}

fn neighbors(topology: &TopologyBlock, atom: usize) -> Vec<(usize, usize)> {
    topology
        .adjacency
        .neighbors_of(atom)
        .iter()
        .map(|neighbor| (neighbor.atom_index, neighbor.bond.index()))
        .collect()
}

#[test]
fn construction_preserves_empty_isolated_chain_branch_and_ring_order() {
    let empty = topology(0, &[]);
    assert!(empty.atoms.is_empty());
    assert!(empty.bonds.is_empty());
    assert!(empty.adjacency.neighbors_of(0).is_empty());

    let isolated = topology(3, &[]);
    assert_eq!(
        (0..3).map(|i| neighbors(&isolated, i)).collect::<Vec<_>>(),
        vec![vec![], vec![], vec![]]
    );

    let chain = topology(4, &[(0, 1), (1, 2), (2, 3)]);
    assert_eq!(
        (0..4).map(|i| neighbors(&chain, i)).collect::<Vec<_>>(),
        vec![
            vec![(1, 0)],
            vec![(0, 0), (2, 1)],
            vec![(1, 1), (3, 2)],
            vec![(2, 2)]
        ]
    );
    assert_eq!(
        chain
            .bonds
            .iter()
            .map(|b| (b.begin(), b.end()))
            .collect::<Vec<_>>(),
        vec![
            (atom_id(0), atom_id(1)),
            (atom_id(1), atom_id(2)),
            (atom_id(2), atom_id(3))
        ]
    );

    let branch = topology(4, &[(1, 3), (1, 0), (1, 2)]);
    assert_eq!(
        (0..4).map(|i| neighbors(&branch, i)).collect::<Vec<_>>(),
        vec![
            vec![(1, 1)],
            vec![(3, 0), (0, 1), (2, 2)],
            vec![(1, 2)],
            vec![(1, 0)]
        ]
    );

    let ring = topology(4, &[(2, 3), (0, 1), (3, 0), (1, 2)]);
    assert_eq!(
        (0..4).map(|i| neighbors(&ring, i)).collect::<Vec<_>>(),
        vec![
            vec![(1, 1), (3, 2)],
            vec![(0, 1), (2, 3)],
            vec![(3, 0), (1, 3)],
            vec![(2, 0), (0, 2)]
        ]
    );
    assert_eq!(
        ring.bonds.iter().map(|b| b.id()).collect::<Vec<_>>(),
        (0..4).map(bond_id).collect::<Vec<_>>()
    );
}

#[test]
fn validation_reports_row_endpoint_loop_and_stereo_errors_exactly() {
    assert_eq!(
        TopologyBlock::try_from_parts(vec![atom(1)], vec![], vec![], vec![]),
        Err(TopologyValidationError::AtomIdMismatch {
            position: 0,
            id: atom_id(1)
        })
    );

    let wrong_bond_id = bond(1, 0, 1);
    let block = TopologyBlock {
        atoms: vec![atom(0), atom(1)],
        adjacency: AdjacencyList::from_topology(2, std::slice::from_ref(&wrong_bond_id)),
        bonds: vec![wrong_bond_id],
        substance_groups: vec![],
        stereo_groups: vec![],
    };
    assert_eq!(
        block.validate(),
        Err(TopologyValidationError::BondIdMismatch {
            position: 0,
            id: bond_id(1)
        })
    );

    for (endpoint, bad_bond) in [("begin", bond(0, 2, 0)), ("end", bond(0, 0, 2))] {
        assert_eq!(
            TopologyBlock::try_from_parts(vec![atom(0)], vec![bad_bond], vec![], vec![]),
            Err(TopologyValidationError::BondEndpointOutOfRange {
                bond: bond_id(0),
                endpoint,
                atom: atom_id(2),
                atom_count: 1
            })
        );
    }

    assert_eq!(
        TopologyBlock::try_from_parts(vec![atom(0)], vec![bond(0, 0, 0)], vec![], vec![]),
        Err(TopologyValidationError::SelfLoopBond {
            bond: bond_id(0),
            atom: atom_id(0)
        })
    );

    let stereo_out_of_range = Bond::from_spec(
        bond_id(0),
        BondSpec::new(atom_id(0), atom_id(1), BondOrder::Double)
            .with_stereo(BondStereo::Cis)
            .with_stereo_atoms(atom_id(2), atom_id(0)),
    );
    assert_eq!(
        TopologyBlock::try_from_parts(
            vec![atom(0), atom(1)],
            vec![stereo_out_of_range],
            vec![],
            vec![]
        ),
        Err(TopologyValidationError::StereoAtomOutOfRange {
            bond: bond_id(0),
            begin: atom_id(2),
            end: atom_id(0),
            atom_count: 2
        })
    );

    let stereo_missing = Bond::from_spec(
        bond_id(0),
        BondSpec::new(atom_id(0), atom_id(1), BondOrder::Double).with_stereo(BondStereo::Trans),
    );
    assert_eq!(
        TopologyBlock::try_from_parts(vec![atom(0), atom(1)], vec![stereo_missing], vec![], vec![]),
        Err(TopologyValidationError::StereoAtomsRequired {
            bond: bond_id(0),
            stereo: BondStereo::Trans
        })
    );
}

#[test]
fn validation_reports_every_substance_group_reference_family() {
    let atoms = vec![atom(0), atom(1)];
    let bonds = vec![bond(0, 0, 1)];
    let invalid_atom = atom_id(2);
    let invalid_bond = bond_id(1);
    let group = |group: SubstanceGroup| {
        TopologyBlock::try_from_parts(atoms.clone(), bonds.clone(), vec![group], vec![])
    };

    assert_eq!(
        group(SubstanceGroup::new(
            SubstanceGroupId::new(1),
            SubstanceGroupKind::Data
        )),
        Err(TopologyValidationError::SubstanceGroupIdMismatch {
            position: 0,
            id: SubstanceGroupId::new(1)
        })
    );
    for invalid in
        [
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![invalid_atom]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_parent_atoms(vec![invalid_atom]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: invalid_atom,
                    leaving_atom: None,
                    label: None,
                    order: None,
                }]),
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: atom_id(0),
                    leaving_atom: Some(invalid_atom),
                    label: None,
                    order: None,
                }]),
        ]
    {
        assert_eq!(
            group(invalid),
            Err(TopologyValidationError::SubstanceGroupAtomOutOfRange {
                sgroup: SubstanceGroupId::new(0),
                atom: invalid_atom,
                atom_count: 2
            })
        );
    }
    for invalid in [
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_bonds(vec![invalid_bond]),
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data).with_cstates(vec![
            SGroupCState {
                bond: invalid_bond,
                vector: [1.0, 2.0, 3.0],
            },
        ]),
    ] {
        assert_eq!(
            group(invalid),
            Err(TopologyValidationError::SubstanceGroupBondOutOfRange {
                sgroup: SubstanceGroupId::new(0),
                bond: invalid_bond,
                bond_count: 1
            })
        );
    }
    assert_eq!(
        group(
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_parent(SubstanceGroupId::new(1))
        ),
        Err(TopologyValidationError::SubstanceGroupParentOutOfRange {
            sgroup: SubstanceGroupId::new(0),
            parent: SubstanceGroupId::new(1)
        })
    );
}

#[test]
fn validation_reports_stereo_group_duplicate_and_stale_adjacency_errors() {
    assert_eq!(
        TopologyBlock::try_from_parts(
            vec![atom(0), atom(1)],
            vec![bond(0, 0, 1)],
            vec![],
            vec![StereoGroup::new(
                StereoGroupKind::Or,
                vec![atom_id(2)],
                vec![]
            )],
        ),
        Err(TopologyValidationError::StereoGroupAtomOutOfRange {
            atom: atom_id(2),
            atom_count: 2
        })
    );
    assert_eq!(
        TopologyBlock::try_from_parts(
            vec![atom(0), atom(1)],
            vec![bond(0, 0, 1)],
            vec![],
            vec![StereoGroup::new(
                StereoGroupKind::Or,
                vec![],
                vec![bond_id(1)]
            )],
        ),
        Err(TopologyValidationError::StereoGroupBondOutOfRange {
            bond: bond_id(1),
            bond_count: 1
        })
    );
    assert_eq!(
        TopologyBlock::try_from_parts(
            vec![atom(0), atom(1), atom(2)],
            vec![bond(0, 0, 1), bond(0, 1, 2)],
            vec![],
            vec![]
        ),
        Err(TopologyValidationError::AdjacencyMismatch)
    );
    assert_eq!(
        TopologyBlock::try_from_parts(
            vec![atom(0), atom(1)],
            vec![bond(0, 0, 1), bond(1, 1, 0)],
            vec![],
            vec![]
        ),
        Err(TopologyValidationError::AdjacencyMismatch)
    );
    let stale = TopologyBlock {
        atoms: vec![atom(0), atom(1)],
        bonds: vec![bond(0, 0, 1)],
        adjacency: AdjacencyList::default(),
        substance_groups: vec![],
        stereo_groups: vec![],
    };
    assert_eq!(
        stale.validate(),
        Err(TopologyValidationError::AdjacencyMismatch)
    );
}

#[test]
fn batch_begin_ranges_duplicate_removal_abort_and_drop_preserve_source() {
    let invalid = TopologyBlock {
        atoms: vec![atom(1)],
        ..TopologyBlock::default()
    };
    assert_eq!(
        invalid.begin_batch_edit().unwrap_err(),
        TopologyEditError::InvalidSource(TopologyValidationError::AtomIdMismatch {
            position: 0,
            id: atom_id(1)
        })
    );

    let source = topology(3, &[(0, 1), (1, 2)]);
    let snapshot = source.clone();
    let mut edit = source.begin_batch_edit().unwrap();
    assert_eq!(
        edit.remove_atom(atom_id(3)),
        Err(TopologyEditError::AtomOutOfRange {
            atom: atom_id(3),
            atom_count: 3
        })
    );
    assert_eq!(
        edit.remove_bond(bond_id(2)),
        Err(TopologyEditError::BondOutOfRange {
            bond: bond_id(2),
            bond_count: 2
        })
    );
    edit.remove_atom(atom_id(1)).unwrap();
    edit.remove_atom(atom_id(1)).unwrap();
    edit.abort();
    assert_eq!(source, snapshot);

    let mut dropped = source.begin_batch_edit().unwrap();
    dropped.remove_bond(bond_id(0)).unwrap();
    drop(dropped);
    assert_eq!(source, snapshot);
}

#[test]
fn no_op_finish_is_identity_and_add_bond_rejects_invalid_rows() {
    let source = topology(3, &[(0, 1)]);
    let (result, mapping) = source.begin_batch_edit().unwrap().finish().unwrap();
    assert_eq!(result, source);
    assert_eq!(
        mapping.atoms().old_to_new(),
        &[Some(atom_id(0)), Some(atom_id(1)), Some(atom_id(2))]
    );
    assert_eq!(mapping.bonds().old_to_new(), &[Some(bond_id(0))]);
    mapping.validate_for_counts(3, 3, 1, 1).unwrap();

    let mut edit = source.begin_batch_edit().unwrap();
    assert_eq!(
        edit.add_bond(BondSpec::new(atom_id(0), atom_id(3), BondOrder::Single)),
        Err(TopologyEditError::AtomOutOfRange {
            atom: atom_id(3),
            atom_count: 3
        })
    );
    assert_eq!(
        edit.add_bond(BondSpec::new(atom_id(2), atom_id(2), BondOrder::Single)),
        Err(TopologyEditError::InvalidResult(
            TopologyValidationError::SelfLoopBond {
                bond: bond_id(1),
                atom: atom_id(2)
            }
        ))
    );
    assert_eq!(
        edit.add_bond(BondSpec::new(atom_id(1), atom_id(0), BondOrder::Single)),
        Err(TopologyEditError::DuplicateBond {
            begin: atom_id(1),
            end: atom_id(0)
        })
    );
}

#[test]
fn batch_atom_bond_mixed_multiple_and_all_removals_have_exact_mappings() {
    let source = topology(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);

    let mut atom_edit = source.begin_batch_edit().unwrap();
    atom_edit.remove_atom(atom_id(1)).unwrap();
    let (atom_result, atom_mapping) = atom_edit.finish().unwrap();
    assert_eq!(
        atom_result
            .bonds
            .iter()
            .map(|b| (b.begin(), b.end()))
            .collect::<Vec<_>>(),
        vec![(atom_id(1), atom_id(2)), (atom_id(2), atom_id(3))]
    );
    assert_eq!(
        atom_mapping.atoms().old_to_new(),
        &[
            Some(atom_id(0)),
            None,
            Some(atom_id(1)),
            Some(atom_id(2)),
            Some(atom_id(3))
        ]
    );
    assert_eq!(
        atom_mapping.bonds().old_to_new(),
        &[None, None, Some(bond_id(0)), Some(bond_id(1))]
    );
    atom_mapping.validate_for_counts(5, 4, 4, 2).unwrap();

    let mut bond_edit = source.begin_batch_edit().unwrap();
    bond_edit.remove_bond(bond_id(2)).unwrap();
    let (bond_result, bond_mapping) = bond_edit.finish().unwrap();
    assert_eq!(bond_result.bonds.len(), 3);
    assert_eq!(
        bond_mapping.bonds().old_to_new(),
        &[Some(bond_id(0)), Some(bond_id(1)), None, Some(bond_id(2))]
    );
    bond_mapping.validate_for_counts(5, 5, 4, 3).unwrap();

    let mut mixed = source.begin_batch_edit().unwrap();
    mixed.remove_atom(atom_id(3)).unwrap();
    mixed.remove_atom(atom_id(1)).unwrap();
    mixed.remove_bond(bond_id(0)).unwrap();
    let (mixed_result, mixed_mapping) = mixed.finish().unwrap();
    assert_eq!(mixed_result.atoms.len(), 3);
    assert!(mixed_result.bonds.is_empty());
    assert_eq!(
        mixed_mapping.atoms().new_to_old(),
        &[Some(atom_id(0)), Some(atom_id(2)), Some(atom_id(4))]
    );
    mixed_mapping.validate_for_counts(5, 3, 4, 0).unwrap();

    let mut all = source.begin_batch_edit().unwrap();
    for atom in (0..5).rev() {
        all.remove_atom(atom_id(atom)).unwrap();
    }
    let (empty, all_mapping) = all.finish().unwrap();
    assert!(empty.atoms.is_empty() && empty.bonds.is_empty());
    assert_eq!(
        all_mapping.atoms().old_to_new(),
        &[None, None, None, None, None]
    );
    assert_eq!(all_mapping.bonds().old_to_new(), &[None, None, None, None]);
    all_mapping.validate_for_counts(5, 0, 4, 0).unwrap();
}

#[test]
fn additions_during_batch_and_removal_of_appended_rows_are_mapped_explicitly() {
    let source = topology(2, &[(0, 1)]);
    let mut edit = source.begin_batch_edit().unwrap();
    let atom2 = edit.add_atom(AtomSpec::new(Element::N));
    let atom3 = edit.add_atom(AtomSpec::new(Element::O));
    let bond1 = edit
        .add_bond(BondSpec::new(atom_id(1), atom2, BondOrder::Single))
        .unwrap();
    let bond2 = edit
        .add_bond(BondSpec::new(atom2, atom3, BondOrder::Double))
        .unwrap();
    edit.remove_atom(atom3).unwrap();
    edit.remove_bond(bond2).unwrap();
    let (result, mapping) = edit.finish().unwrap();
    assert_eq!(result.atoms.len(), 3);
    assert_eq!(result.atoms[2].element(), Element::N);
    assert_eq!(result.bonds.len(), 2);
    assert_eq!(result.bonds[1].id(), bond1);
    assert_eq!(
        mapping.atoms().old_to_new(),
        &[Some(atom_id(0)), Some(atom_id(1))]
    );
    assert_eq!(
        mapping.atoms().new_to_old(),
        &[Some(atom_id(0)), Some(atom_id(1)), None]
    );
    assert_eq!(mapping.bonds().old_to_new(), &[Some(bond_id(0))]);
    assert_eq!(mapping.bonds().new_to_old(), &[Some(bond_id(0)), None]);
    mapping.validate_for_counts(2, 3, 1, 2).unwrap();
}

#[test]
fn invalid_stereo_reference_is_rejected_before_finish_and_does_not_change_source() {
    let source = topology(2, &[]);
    let snapshot = source.clone();
    let mut edit = source.begin_batch_edit().unwrap();
    assert_eq!(
        edit.add_bond(
            BondSpec::new(atom_id(0), atom_id(1), BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(atom_id(0), atom_id(9)),
        ),
        Err(TopologyEditError::InvalidResult(
            TopologyValidationError::StereoAtomOutOfRange {
                bond: bond_id(0),
                begin: atom_id(0),
                end: atom_id(9),
                atom_count: 2,
            }
        ))
    );
    let (result, mapping) = edit.finish().unwrap();
    assert_eq!(result, source);
    mapping.validate_for_counts(2, 2, 0, 0).unwrap();
    assert_eq!(source, snapshot);
}

#[test]
fn batch_remaps_retained_sgroups_and_drops_direct_and_transitive_dependents() {
    let mut source = topology(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    source.substance_groups = vec![
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![atom_id(2), atom_id(4)])
            .with_bonds(vec![bond_id(2)])
            .with_parent_atoms(vec![atom_id(3)])
            .with_attach_points(vec![SGroupAttachPoint {
                atom: atom_id(4),
                leaving_atom: Some(atom_id(3)),
                label: Some("AP".into()),
                order: Some(2),
            }])
            .with_cstates(vec![SGroupCState {
                bond: bond_id(2),
                vector: [3.0, 4.0, 5.0],
            }]),
        SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Superatom)
            .with_atoms(vec![atom_id(1)]),
        SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Monomer)
            .with_atoms(vec![atom_id(4)])
            .with_parent(SubstanceGroupId::new(1)),
    ];
    source.validate().unwrap();

    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(atom_id(1)).unwrap();
    let (result, _) = edit.finish().unwrap();
    assert_eq!(result.substance_groups.len(), 1);
    let retained = &result.substance_groups[0];
    assert_eq!(retained.id(), SubstanceGroupId::new(0));
    assert_eq!(retained.atoms(), &[atom_id(1), atom_id(3)]);
    assert_eq!(retained.bonds(), &[bond_id(0)]);
    assert_eq!(retained.parent_atoms(), &[atom_id(2)]);
    assert_eq!(retained.attach_points()[0].atom, atom_id(3));
    assert_eq!(retained.attach_points()[0].leaving_atom, Some(atom_id(2)));
    assert_eq!(retained.cstates()[0].bond, bond_id(0));
    result.validate().unwrap();
}

#[test]
fn typed_sgroup_bond_references_validate_and_preserve_order_through_compaction() {
    let mut source = topology(6, &[(0, 1), (1, 2), (2, 3), (3, 4)]);
    source.substance_groups = vec![
        SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_atoms(vec![atom_id(3), atom_id(4)])
        .with_bonds(vec![bond_id(2), bond_id(3)])
        .with_head_crossing_bonds(vec![bond_id(3), bond_id(2), bond_id(3)])
        .with_crossing_bond_correspondence(vec![bond_id(2), bond_id(3), bond_id(2)])
        .with_cstates(vec![SGroupCState {
            bond: bond_id(3),
            vector: [1.0, 2.0, 3.0],
        }]),
    ];
    source
        .validate()
        .expect("all typed references are in range");

    let mut compact = source.begin_batch_edit().expect("valid source");
    compact.remove_atom(atom_id(0)).expect("valid removal");
    let (compacted, mapping) = compact.finish().expect("surviving references remap");
    assert_eq!(
        mapping.bonds().old_to_new(),
        &[None, Some(bond_id(0)), Some(bond_id(1)), Some(bond_id(2))]
    );
    let group = &compacted.substance_groups[0];
    assert_eq!(group.bonds(), &[bond_id(1), bond_id(2)]);
    assert_eq!(
        group.head_crossing_bonds(),
        &[bond_id(2), bond_id(1), bond_id(2)]
    );
    assert_eq!(
        group.crossing_bond_correspondence(),
        &[bond_id(1), bond_id(2), bond_id(1)]
    );
    assert_eq!(group.cstates()[0].bond, bond_id(2));
    compacted.validate().expect("remapped topology is valid");

    let isolated_atom_order = [
        atom_id(5),
        atom_id(0),
        atom_id(1),
        atom_id(2),
        atom_id(3),
        atom_id(4),
    ];
    let (reordered, _) = source
        .reordered_atoms(&isolated_atom_order)
        .expect("atom-only reorder preserves the bond table");
    assert_eq!(
        reordered.substance_groups[0].head_crossing_bonds(),
        &[bond_id(3), bond_id(2), bond_id(3)]
    );
    assert_eq!(
        reordered.substance_groups[0].crossing_bond_correspondence(),
        &[bond_id(2), bond_id(3), bond_id(2)]
    );
}

#[test]
fn typed_sgroup_bond_reference_removal_drops_the_complete_group() {
    for group in [
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![atom_id(0)])
            .with_head_crossing_bonds(vec![bond_id(1)]),
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_atoms(vec![atom_id(0)])
            .with_crossing_bond_correspondence(vec![bond_id(1)]),
    ] {
        let mut source = topology(4, &[(0, 1), (1, 2), (2, 3)]);
        source.substance_groups = vec![group];
        source.validate().expect("reference begins valid");
        let mut edit = source.begin_batch_edit().expect("valid source");
        edit.remove_bond(bond_id(1)).expect("valid bond removal");
        let (result, _) = edit.finish().expect("group removal is atomic");
        assert!(result.substance_groups.is_empty());
        result.validate().expect("no stale reference survives");
    }
}

#[test]
fn typed_sgroup_bond_references_reject_each_out_of_range_category() {
    let invalid = bond_id(1);
    for group in [
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_head_crossing_bonds(vec![invalid]),
        SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_crossing_bond_correspondence(vec![invalid]),
    ] {
        let result = TopologyBlock::try_from_parts(
            vec![atom(0), atom(1)],
            vec![bond(0, 0, 1)],
            vec![group],
            vec![],
        );
        assert_eq!(
            result,
            Err(TopologyValidationError::SubstanceGroupBondOutOfRange {
                sgroup: SubstanceGroupId::new(0),
                bond: invalid,
                bond_count: 1,
            })
        );
    }
}

#[test]
fn enhanced_stereo_removal_retains_partial_members_in_order_and_drops_empty_group() {
    let mut source = topology(4, &[(0, 1), (1, 2), (2, 3)]);
    source.stereo_groups = vec![
        StereoGroup::new(
            StereoGroupKind::Or,
            vec![atom_id(0), atom_id(2), atom_id(3)],
            vec![bond_id(0), bond_id(2)],
        )
        .with_id(7),
        StereoGroup::new(StereoGroupKind::And, vec![atom_id(1)], vec![bond_id(0)]).with_id(8),
    ];
    source.validate().unwrap();
    let mut edit = source.begin_batch_edit().unwrap();
    edit.remove_atom(atom_id(1)).unwrap();
    let (result, _) = edit.finish().unwrap();
    assert_eq!(result.stereo_groups.len(), 1);
    assert_eq!(result.stereo_groups[0].id(), Some(7));
    assert_eq!(
        result.stereo_groups[0].atoms(),
        &[atom_id(0), atom_id(1), atom_id(2)]
    );
    assert_eq!(result.stereo_groups[0].bonds(), &[bond_id(0)]);
}

fn stereo_source() -> TopologyBlock {
    let atoms = (0..5).map(atom).collect();
    let bonds = vec![
        bond(0, 0, 1),
        Bond::from_spec(
            bond_id(1),
            BondSpec::new(atom_id(1), atom_id(2), BondOrder::Double)
                .with_stereo(BondStereo::Cis)
                .with_stereo_atoms(atom_id(0), atom_id(3)),
        ),
        bond(2, 2, 3),
    ];
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

#[test]
fn double_bond_stereo_is_remapped_when_retained_and_cleared_when_carrier_is_lost() {
    let source = stereo_source();
    let mut retain = source.begin_batch_edit().unwrap();
    retain.remove_atom(atom_id(4)).unwrap();
    let (retained, _) = retain.finish().unwrap();
    assert_eq!(retained.bonds[1].stereo(), BondStereo::Cis);
    assert_eq!(
        retained.bonds[1].stereo_atoms(),
        Some([atom_id(0), atom_id(3)])
    );

    let mut lose_atom = source.begin_batch_edit().unwrap();
    lose_atom.remove_atom(atom_id(0)).unwrap();
    let (atom_result, _) = lose_atom.finish().unwrap();
    assert_eq!(atom_result.bonds[0].stereo(), BondStereo::None);
    assert_eq!(atom_result.bonds[0].stereo_atoms(), None);

    let mut lose_bond = source.begin_batch_edit().unwrap();
    lose_bond.remove_bond(bond_id(2)).unwrap();
    let (bond_result, _) = lose_bond.finish().unwrap();
    assert_eq!(bond_result.bonds[1].stereo(), BondStereo::None);
    assert_eq!(bond_result.bonds[1].stereo_atoms(), None);
}

#[test]
fn reorder_remaps_all_references_preserves_bond_order_and_reports_every_permutation_error() {
    let atoms = vec![
        Atom::from_spec(atom_id(0), AtomSpec::new(Element::C)),
        Atom::from_spec(atom_id(1), AtomSpec::new(Element::N)),
        Atom::from_spec(atom_id(2), AtomSpec::new(Element::O)),
        Atom::from_spec(atom_id(3), AtomSpec::new(Element::S)),
    ];
    let bonds = vec![
        bond(0, 0, 1),
        Bond::from_spec(
            bond_id(1),
            BondSpec::new(atom_id(1), atom_id(2), BondOrder::Double)
                .with_stereo(BondStereo::Trans)
                .with_stereo_atoms(atom_id(0), atom_id(3)),
        ),
        bond(2, 2, 3),
    ];
    let sgroup = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
        .with_atoms(vec![atom_id(3), atom_id(0)])
        .with_bonds(vec![bond_id(2), bond_id(0)])
        .with_parent_atoms(vec![atom_id(1)])
        .with_attach_points(vec![SGroupAttachPoint {
            atom: atom_id(2),
            leaving_atom: Some(atom_id(0)),
            label: None,
            order: None,
        }])
        .with_cstates(vec![SGroupCState {
            bond: bond_id(1),
            vector: [1.0, 0.0, -0.0],
        }]);
    let source = TopologyBlock::try_from_parts(
        atoms,
        bonds,
        vec![sgroup],
        vec![StereoGroup::new(
            StereoGroupKind::Absolute,
            vec![atom_id(3), atom_id(1)],
            vec![bond_id(2), bond_id(0)],
        )],
    )
    .unwrap();
    let snapshot = source.clone();

    let (reordered, mapping) = source
        .reordered_atoms(&[atom_id(2), atom_id(0), atom_id(3), atom_id(1)])
        .unwrap();
    assert_eq!(
        reordered
            .atoms
            .iter()
            .map(Atom::element)
            .collect::<Vec<_>>(),
        vec![Element::O, Element::C, Element::S, Element::N]
    );
    assert_eq!(
        mapping.atoms().old_to_new(),
        &[
            Some(atom_id(1)),
            Some(atom_id(3)),
            Some(atom_id(0)),
            Some(atom_id(2))
        ]
    );
    assert_eq!(
        mapping.bonds().old_to_new(),
        &[Some(bond_id(0)), Some(bond_id(1)), Some(bond_id(2))]
    );
    mapping.validate_for_counts(4, 4, 3, 3).unwrap();
    assert_eq!(
        reordered
            .bonds
            .iter()
            .map(|b| (b.begin(), b.end()))
            .collect::<Vec<_>>(),
        vec![
            (atom_id(1), atom_id(3)),
            (atom_id(3), atom_id(0)),
            (atom_id(0), atom_id(2))
        ]
    );
    assert_eq!(
        reordered.bonds[1].stereo_atoms(),
        Some([atom_id(1), atom_id(2)])
    );
    let group = &reordered.substance_groups[0];
    assert_eq!(group.atoms(), &[atom_id(2), atom_id(1)]);
    assert_eq!(group.bonds(), &[bond_id(2), bond_id(0)]);
    assert_eq!(group.parent_atoms(), &[atom_id(3)]);
    assert_eq!(group.attach_points()[0].atom, atom_id(0));
    assert_eq!(group.attach_points()[0].leaving_atom, Some(atom_id(1)));
    assert_eq!(group.cstates()[0].bond, bond_id(1));
    assert_eq!(
        reordered.stereo_groups[0].atoms(),
        &[atom_id(2), atom_id(3)]
    );
    assert_eq!(
        reordered.stereo_groups[0].bonds(),
        &[bond_id(2), bond_id(0)]
    );
    assert_eq!(
        (0..4).map(|i| neighbors(&reordered, i)).collect::<Vec<_>>(),
        vec![
            vec![(3, 1), (2, 2)],
            vec![(3, 0)],
            vec![(0, 2)],
            vec![(1, 0), (0, 1)]
        ]
    );
    assert_eq!(source, snapshot);

    let (identity, identity_mapping) = source
        .reordered_atoms(&[atom_id(0), atom_id(1), atom_id(2), atom_id(3)])
        .unwrap();
    assert_eq!(identity, source);
    identity_mapping.validate_for_counts(4, 4, 3, 3).unwrap();
    assert_eq!(
        source.reordered_atoms(&[atom_id(0), atom_id(1)]),
        Err(TopologyEditError::PermutationLength {
            actual: 2,
            expected: 4
        })
    );
    assert_eq!(
        source.reordered_atoms(&[atom_id(0), atom_id(1), atom_id(4), atom_id(2)]),
        Err(TopologyEditError::PermutationAtomOutOfRange {
            position: 2,
            atom: atom_id(4),
            atom_count: 4
        })
    );
    assert_eq!(
        source.reordered_atoms(&[atom_id(0), atom_id(1), atom_id(1), atom_id(2)]),
        Err(TopologyEditError::PermutationDuplicateAtom {
            position: 2,
            atom: atom_id(1)
        })
    );
}

#[test]
fn no_bond_and_empty_batch_regressions_remain_valid() {
    let isolated = topology(3, &[]);
    let mut edit = isolated.begin_batch_edit().unwrap();
    edit.remove_atom(atom_id(1)).unwrap();
    let (result, mapping) = edit.finish().unwrap();
    assert_eq!(result.atoms.len(), 2);
    assert!(result.bonds.is_empty());
    mapping.validate_for_counts(3, 2, 0, 0).unwrap();

    let empty = TopologyBlock::default();
    let (result, mapping) = empty.begin_batch_edit().unwrap().finish().unwrap();
    assert_eq!(result, empty);
    mapping.validate_for_counts(0, 0, 0, 0).unwrap();
    assert_eq!(result.adjacency.neighbors_of(0), &[] as &[NeighborRef]);
}
