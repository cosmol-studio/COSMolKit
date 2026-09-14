use std::collections::BTreeSet;

use cosmolkit_core::{
    AtropisomerAssignment, AtropisomerBondUpdate, AtropisomerConformer, AtropisomerError,
    AtropisomerRejectionKind, RingFindType, RingInfo, RingSearchParams, atropisomer_carriers,
    cleanup_atropisomer_stereo_groups, detect_atropisomer_chirality,
    does_topology_have_atropisomers, find_sssr, stereo_group_atom_ids,
    wedge_bonds_from_atropisomers,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, StereoGroup,
    StereoGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Element, Hybridization};

fn atom(id: usize, hybridization: Hybridization) -> Atom {
    Atom::from_spec(
        AtomId::new(id),
        AtomSpec::new(Element::C).with_hybridization(hybridization),
    )
}

fn bond(
    id: usize,
    begin: usize,
    end: usize,
    order: BondOrder,
    direction: BondDirection,
    stereo: BondStereo,
) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
            .with_direction(direction)
            .with_stereo(stereo),
    )
}

fn topology_with_groups(
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    groups: Vec<StereoGroup>,
) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], groups).unwrap()
}

fn axial(left: BondDirection, right: BondDirection, axial_stereo: BondStereo) -> TopologyBlock {
    topology_with_groups(
        vec![
            atom(0, Hybridization::Sp3),
            atom(1, Hybridization::Sp2),
            atom(2, Hybridization::Sp2),
            atom(3, Hybridization::Sp3),
        ],
        vec![
            bond(0, 1, 0, BondOrder::Single, left, BondStereo::None),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                axial_stereo,
            ),
            bond(2, 2, 3, BondOrder::Single, right, BondStereo::None),
        ],
        vec![],
    )
}

fn sssr(topology: &TopologyBlock) -> RingInfo {
    find_sssr(topology, &RingSearchParams::default()).unwrap()
}

#[test]
fn carrier_query_rejects_zero_and_preserves_single_carriers_in_endpoint_order() {
    let no_carriers = topology_with_groups(
        vec![atom(0, Hybridization::Sp2), atom(1, Hybridization::Sp2)],
        vec![bond(
            0,
            0,
            1,
            BondOrder::Single,
            BondDirection::None,
            BondStereo::AtropCw,
        )],
        vec![],
    );
    assert_eq!(
        atropisomer_carriers(&no_carriers, BondId::new(0)).unwrap(),
        None
    );

    let one_each = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    let ends = atropisomer_carriers(&one_each, BondId::new(1))
        .unwrap()
        .unwrap();
    assert_eq!(ends[0].focus(), AtomId::new(1));
    assert_eq!(ends[0].carrier_bonds(), &[BondId::new(0)]);
    assert_eq!(ends[1].focus(), AtomId::new(2));
    assert_eq!(ends[1].carrier_bonds(), &[BondId::new(2)]);
}

#[test]
fn carrier_query_sorts_exactly_two_but_preserves_many_in_adjacency_order() {
    let topology = topology_with_groups(
        (0..8)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                7,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                3,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                2,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::AtropCw,
            ),
            bond(
                3,
                1,
                5,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                4,
                2,
                6,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                5,
                2,
                4,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
        ],
        vec![],
    );

    let ends = atropisomer_carriers(&topology, BondId::new(2))
        .unwrap()
        .unwrap();
    assert_eq!(ends[0].focus(), AtomId::new(1));
    assert_eq!(
        ends[0].carrier_bonds(),
        &[BondId::new(0), BondId::new(1), BondId::new(3)]
    );
    assert_eq!(ends[1].focus(), AtomId::new(2));
    assert_eq!(ends[1].carrier_bonds(), &[BondId::new(5), BondId::new(4)]);

    assert!(
        detect_atropisomer_chirality(&topology, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );
}

#[test]
fn carrier_query_reports_invalid_axial_bond_without_exposing_mutable_state() {
    let topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    let error = atropisomer_carriers(&topology, BondId::new(9)).unwrap_err();
    assert_eq!(
        error,
        AtropisomerError::AxialBondOutOfRange {
            bond: BondId::new(9),
            bond_count: 3,
        }
    );
}

#[test]
fn empty_predicate_and_detection_are_total_and_deterministic() {
    let empty = TopologyBlock::default();
    assert!(!does_topology_have_atropisomers(&empty));
    let first = detect_atropisomer_chirality(&empty, None).unwrap();
    let second = detect_atropisomer_chirality(&empty, None).unwrap();
    assert_eq!(first, second);
    assert!(first.bond_updates.is_empty());
    assert!(first.diagnostics.is_empty());

    assert!(does_topology_have_atropisomers(&axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    )));
}

#[test]
fn no_conformer_wedge_hash_pairs_assign_bond_cw_and_ccw() {
    let ccw = detect_atropisomer_chirality(
        &axial(
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondStereo::None,
        ),
        None,
    )
    .unwrap();
    assert_eq!(
        ccw.bond_updates,
        vec![AtropisomerBondUpdate {
            bond: BondId::new(1),
            stereo: BondStereo::AtropCcw,
        }]
    );

    let cw = detect_atropisomer_chirality(
        &axial(
            BondDirection::BeginDash,
            BondDirection::BeginWedge,
            BondStereo::None,
        ),
        None,
    )
    .unwrap();
    assert_eq!(cw.bond_updates[0].stereo, BondStereo::AtropCw);
}

#[test]
fn candidate_gates_reject_any_non_single_degree_and_non_sp2_controls() {
    let mut any = axial(
        BondDirection::BeginWedge,
        BondDirection::BeginDash,
        BondStereo::Any,
    );
    assert!(
        detect_atropisomer_chirality(&any, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );

    any.bonds[1].set_order(BondOrder::Double);
    assert!(
        detect_atropisomer_chirality(&any, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );

    let mut sulfinamide_like = axial(
        BondDirection::BeginWedge,
        BondDirection::BeginDash,
        BondStereo::None,
    );
    sulfinamide_like.atoms[1].set_hybridization(Hybridization::Sp3);
    assert!(
        detect_atropisomer_chirality(&sulfinamide_like, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );

    let too_high_degree = topology_with_groups(
        (0..6)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                0,
                BondOrder::Single,
                BondDirection::BeginWedge,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                2,
                2,
                3,
                BondOrder::Single,
                BondDirection::BeginDash,
                BondStereo::None,
            ),
            bond(
                3,
                1,
                4,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                4,
                1,
                5,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
        ],
        vec![],
    );
    assert!(
        detect_atropisomer_chirality(&too_high_degree, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );
}

#[test]
fn missing_unknown_and_inconsistent_carriers_do_not_fabricate_stereo() {
    let missing = topology_with_groups(
        vec![atom(0, Hybridization::Sp2), atom(1, Hybridization::Sp2)],
        vec![bond(
            0,
            0,
            1,
            BondOrder::Single,
            BondDirection::BeginWedge,
            BondStereo::None,
        )],
        vec![],
    );
    assert!(
        detect_atropisomer_chirality(&missing, None)
            .unwrap()
            .bond_updates
            .is_empty()
    );

    let unknown = topology_with_groups(
        (0..5)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                0,
                BondOrder::Single,
                BondDirection::BeginWedge,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                2,
                2,
                3,
                BondOrder::Single,
                BondDirection::BeginDash,
                BondStereo::None,
            ),
            bond(
                3,
                1,
                4,
                BondOrder::Single,
                BondDirection::Unknown,
                BondStereo::None,
            ),
        ],
        vec![],
    );
    let result = detect_atropisomer_chirality(&unknown, None).unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::UnknownCarrierDirection
    );

    let same_end = topology_with_groups(
        (0..5)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                0,
                BondOrder::Single,
                BondDirection::BeginWedge,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                2,
                2,
                3,
                BondOrder::Single,
                BondDirection::BeginDash,
                BondStereo::None,
            ),
            bond(
                3,
                1,
                4,
                BondOrder::Single,
                BondDirection::BeginWedge,
                BondStereo::None,
            ),
        ],
        vec![],
    );
    let result = detect_atropisomer_chirality(&same_end, None).unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::InconsistentDirections
    );
}

#[test]
fn two_dimensional_geometry_assigns_and_rejects_zero_axis() {
    let topology = axial(
        BondDirection::BeginWedge,
        BondDirection::None,
        BondStereo::None,
    );
    let conformer = Conformer2D::new(7, vec![[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let result =
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(&conformer)))
            .unwrap();
    assert_eq!(result.bond_updates[0].stereo, BondStereo::AtropCcw);

    let zero = Conformer2D::new(8, vec![[0.0, 1.0], [0.0, 0.0], [0.0, 0.0], [1.0, 1.0]]);
    let result =
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(&zero))).unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::ZeroLengthAxis
    );
}

#[test]
fn three_dimensional_geometry_covers_frame_branches_and_coplanarity() {
    let topology = axial(
        BondDirection::BeginWedge,
        BondDirection::None,
        BondStereo::None,
    );
    let along_x = Conformer3D::new(
        0,
        vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
        ],
        true,
    );
    let result =
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(&along_x)))
            .unwrap();
    assert_eq!(result.bond_updates[0].stereo, BondStereo::AtropCw);

    let along_z = Conformer3D::new(
        1,
        vec![
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        true,
    );
    assert_eq!(
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(&along_z)),)
            .unwrap()
            .bond_updates
            .len(),
        1
    );

    let coplanar = Conformer3D::new(
        2,
        vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
        ],
        true,
    );
    let result =
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(&coplanar)))
            .unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::CoplanarCarriers
    );
}

#[test]
fn two_carrier_projection_rejects_same_side_and_uses_collinear_fallback() {
    let topology = topology_with_groups(
        (0..6)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                0,
                BondOrder::Single,
                BondDirection::BeginWedge,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                2,
                2,
                3,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                3,
                1,
                4,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                4,
                2,
                5,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
        ],
        vec![],
    );
    let same_side = Conformer2D::new(
        0,
        vec![
            [0.0, 1.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 2.0],
            [1.0, -1.0],
        ],
    );
    let result =
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(&same_side)))
            .unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::SameSideCarriers
    );

    let fallback = Conformer3D::new(
        1,
        vec![
            [2.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ],
        true,
    );
    assert_eq!(
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(&fallback)),)
            .unwrap()
            .bond_updates
            .len(),
        1
    );
}

#[test]
fn cleanup_replaces_group_atoms_with_unique_atrop_bonds_and_preserves_other_groups() {
    let atrop_group = StereoGroup::new(
        StereoGroupKind::Or,
        vec![AtomId::new(1), AtomId::new(2)],
        vec![],
    )
    .with_id(17);
    let unaffected =
        StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(0)], vec![]).with_id(19);
    let mut topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    topology.stereo_groups = vec![atrop_group, unaffected.clone()];
    topology.validate().unwrap();

    let cleaned =
        cleanup_atropisomer_stereo_groups(&topology, &AtropisomerAssignment::default()).unwrap();
    assert_eq!(cleaned.groups[0].id(), Some(17));
    assert!(cleaned.groups[0].atoms().is_empty());
    assert_eq!(cleaned.groups[0].bonds(), &[BondId::new(1)]);
    assert_eq!(cleaned.groups[1], unaffected);
}

#[test]
fn stereo_group_atom_expansion_uses_existing_and_generated_wedges_once() {
    let topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    let group = StereoGroup::new(StereoGroupKind::Absolute, vec![], vec![BondId::new(1)]);
    let wedges =
        wedge_bonds_from_atropisomers(&topology, &sssr(&topology), None, &BTreeSet::new()).unwrap();
    let ids = stereo_group_atom_ids(&topology, &group, &wedges).unwrap();
    assert_eq!(ids.len(), 1);
    assert!(ids[0] == AtomId::new(1) || ids[0] == AtomId::new(2));
}

#[test]
fn wedge_generation_covers_no_conformer_2d_and_3d_direction_rules() {
    let topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    let rings = sssr(&topology);
    let no_conf = wedge_bonds_from_atropisomers(&topology, &rings, None, &BTreeSet::new()).unwrap();
    assert_eq!(no_conf.bond_updates.len(), 1);
    assert_eq!(no_conf.bond_updates[0].direction, BondDirection::BeginWedge);
    assert_eq!(no_conf.bond_updates[0].begin, AtomId::new(2));

    let two_d = Conformer2D::new(0, vec![[0.0, 1.0], [0.0, 0.0], [1.0, 0.0], [1.0, 1.0]]);
    let result = wedge_bonds_from_atropisomers(
        &topology,
        &rings,
        Some(AtropisomerConformer::TwoD(&two_d)),
        &BTreeSet::new(),
    )
    .unwrap();
    assert_eq!(result.bond_updates.len(), 1);

    let three_d = Conformer3D::new(
        1,
        vec![
            [0.0, 0.0, 1.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.0, -1.0],
        ],
        true,
    );
    let result = wedge_bonds_from_atropisomers(
        &topology,
        &rings,
        Some(AtropisomerConformer::ThreeD(&three_d)),
        &BTreeSet::new(),
    )
    .unwrap();
    assert_eq!(result.bond_updates.len(), 1);
    let update = result.bond_updates[0];
    let z_delta = three_d.coordinates()[update.end.index()][2]
        - three_d.coordinates()[update.begin.index()][2];
    assert_eq!(
        update.direction,
        if z_delta > 0.000_000_1 {
            BondDirection::BeginWedge
        } else {
            BondDirection::BeginDash
        }
    );
}

#[test]
fn two_dimensional_wedge_selection_prefers_ring_carriers_over_nonring() {
    let topology = topology_with_groups(
        (0..5)
            .map(|id| {
                atom(
                    id,
                    if id == 1 || id == 2 {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    },
                )
            })
            .collect(),
        vec![
            bond(
                0,
                1,
                0,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                1,
                1,
                2,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::AtropCcw,
            ),
            bond(
                2,
                2,
                3,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                3,
                0,
                4,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
            bond(
                4,
                4,
                1,
                BondOrder::Single,
                BondDirection::None,
                BondStereo::None,
            ),
        ],
        vec![],
    );
    let rings = sssr(&topology);
    let coordinates = Conformer2D::new(
        0,
        vec![
            [0.0, 1.0],
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, -1.0],
            [-1.0, -0.5],
        ],
    );
    let result = wedge_bonds_from_atropisomers(
        &topology,
        &rings,
        Some(AtropisomerConformer::TwoD(&coordinates)),
        &BTreeSet::new(),
    )
    .unwrap();
    assert_eq!(result.bond_updates.len(), 1);
    assert!(matches!(result.bond_updates[0].bond.index(), 0 | 4));
}

#[test]
fn occupied_conflicting_and_unknown_carriers_fail_closed() {
    let topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCcw,
    );
    let rings = sssr(&topology);
    let occupied = BTreeSet::from([BondId::new(0), BondId::new(2)]);
    let result = wedge_bonds_from_atropisomers(&topology, &rings, None, &occupied).unwrap();
    assert!(result.bond_updates.is_empty());
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::NoUsableWedgeBond
    );

    let conflict = axial(
        BondDirection::EndUpRight,
        BondDirection::None,
        BondStereo::AtropCcw,
    );
    let result =
        wedge_bonds_from_atropisomers(&conflict, &sssr(&conflict), None, &BTreeSet::new()).unwrap();
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::DirectionConflict
    );

    let unknown = axial(
        BondDirection::Unknown,
        BondDirection::None,
        BondStereo::AtropCcw,
    );
    let result =
        wedge_bonds_from_atropisomers(&unknown, &sssr(&unknown), None, &BTreeSet::new()).unwrap();
    assert_eq!(
        result.diagnostics[0].kind,
        AtropisomerRejectionKind::UnknownCarrierDirection
    );
}

#[test]
fn structured_validation_rejects_coordinates_rings_assignments_and_group_ids() {
    let topology = axial(
        BondDirection::None,
        BondDirection::None,
        BondStereo::AtropCw,
    );
    let short = Conformer2D::new(4, vec![[0.0, 0.0]]);
    assert!(matches!(
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::TwoD(&short)),),
        Err(AtropisomerError::InvalidCoordinates { .. })
    ));
    let false_3d = Conformer3D::new(5, vec![[0.0; 3]; 4], false);
    assert_eq!(
        detect_atropisomer_chirality(&topology, Some(AtropisomerConformer::ThreeD(&false_3d)),),
        Err(AtropisomerError::ConformerNotThreeDimensional { conformer: 5 })
    );
    assert_eq!(
        wedge_bonds_from_atropisomers(
            &topology,
            &RingInfo::new(RingFindType::Fast, 4, 3),
            None,
            &BTreeSet::new(),
        ),
        Err(AtropisomerError::RingInfoNotSssr)
    );
    assert_eq!(
        cleanup_atropisomer_stereo_groups(
            &topology,
            &AtropisomerAssignment {
                bond_updates: vec![AtropisomerBondUpdate {
                    bond: BondId::new(99),
                    stereo: BondStereo::AtropCw,
                }],
                diagnostics: vec![],
            },
        ),
        Err(AtropisomerError::AssignmentBondOutOfRange {
            bond: BondId::new(99),
            bond_count: 3,
        })
    );
    let bad_group = StereoGroup::new(StereoGroupKind::Absolute, vec![], vec![BondId::new(99)]);
    assert_eq!(
        stereo_group_atom_ids(&topology, &bad_group, &Default::default(),),
        Err(AtropisomerError::StereoGroupBondOutOfRange {
            bond: BondId::new(99),
            bond_count: 3,
        })
    );
}
