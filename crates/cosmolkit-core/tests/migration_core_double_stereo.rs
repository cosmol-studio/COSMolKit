use cosmolkit_core::{
    DoubleBondControl, DoubleBondStereoDescriptor, DoubleBondStereoError,
    DoubleBondStereoSpecified, RingFindType, RingInfo, assign_directional_double_bond_stereo,
    assign_double_bond_stereo_from_directions, clear_bond_directions, clear_single_bond_directions,
    double_bond_stereo_info, fast_find_rings, find_double_bond_stereo_atoms,
    has_stereo_bond_direction, is_double_bond_stereo_candidate, neighboring_directed_bond,
    opposite_stereo_bond_direction, set_double_bond_neighbor_directions,
    should_detect_double_bond_stereo, translate_ez_to_cis_trans, with_double_bond_stereo_reference,
};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer3D, CoordinateValidationError,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, Element};

fn atom(id: usize) -> Atom {
    Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C))
}

fn unknown_atom(id: usize) -> Atom {
    Atom::from_spec(
        AtomId::new(id),
        AtomSpec::new(Element::C).with_unknown_stereo(true),
    )
}

fn bond(id: usize, begin: usize, end: usize, order: BondOrder) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order),
    )
}

fn directed_bond(
    id: usize,
    begin: usize,
    end: usize,
    order: BondOrder,
    direction: BondDirection,
) -> Bond {
    Bond::from_spec(
        BondId::new(id),
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_direction(direction),
    )
}

fn stereo_bond(
    id: usize,
    begin: usize,
    end: usize,
    stereo: BondStereo,
    references: Option<[usize; 2]>,
) -> Bond {
    let mut spec =
        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Double).with_stereo(stereo);
    if let Some([begin_reference, end_reference]) = references {
        spec = spec.with_stereo_atoms(AtomId::new(begin_reference), AtomId::new(end_reference));
    }
    Bond::from_spec(BondId::new(id), spec)
}

fn topology_with_atoms(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

fn topology(atom_count: usize, bonds: Vec<Bond>) -> TopologyBlock {
    topology_with_atoms((0..atom_count).map(atom).collect(), bonds)
}

fn empty_rings(topology: &TopologyBlock) -> RingInfo {
    RingInfo::new(
        RingFindType::Fast,
        topology.atoms.len(),
        topology.bonds.len(),
    )
}

fn four_atom_chain(left: BondDirection, right: BondDirection, stereo: BondStereo) -> TopologyBlock {
    topology(
        4,
        vec![
            directed_bond(0, 0, 1, BondOrder::Single, left),
            stereo_bond(
                1,
                1,
                2,
                stereo,
                matches!(
                    stereo,
                    BondStereo::Cis | BondStereo::Trans | BondStereo::E | BondStereo::Z
                )
                .then_some([0, 3]),
            ),
            directed_bond(2, 2, 3, BondOrder::Single, right),
        ],
    )
}

fn ring(size: usize) -> TopologyBlock {
    let mut bonds = Vec::with_capacity(size);
    for index in 0..size {
        bonds.push(bond(
            index,
            index,
            (index + 1) % size,
            if index == 0 {
                BondOrder::Double
            } else {
                BondOrder::Single
            },
        ));
    }
    topology(size, bonds)
}

#[test]
fn direction_predicate_and_opposite_are_exact_and_structured() {
    for direction in [BondDirection::EndDownRight, BondDirection::EndUpRight] {
        assert!(has_stereo_bond_direction(direction));
    }
    for direction in [
        BondDirection::None,
        BondDirection::BeginWedge,
        BondDirection::BeginDash,
        BondDirection::EitherDouble,
        BondDirection::Unknown,
    ] {
        assert!(!has_stereo_bond_direction(direction));
    }
    assert_eq!(
        opposite_stereo_bond_direction(BondDirection::EndDownRight),
        Ok(BondDirection::EndUpRight)
    );
    assert_eq!(
        opposite_stereo_bond_direction(BondDirection::EndUpRight),
        Ok(BondDirection::EndDownRight)
    );
    assert_eq!(
        opposite_stereo_bond_direction(BondDirection::Unknown),
        Err(DoubleBondStereoError::InvalidDirection {
            direction: BondDirection::Unknown,
        })
    );
}

#[test]
fn ring_threshold_and_candidate_gates_match_source() {
    let acyclic = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::None);
    let acyclic_rings = empty_rings(&acyclic);
    assert_eq!(
        should_detect_double_bond_stereo(&acyclic, &acyclic_rings, BondId::new(1)),
        Ok(true)
    );
    assert_eq!(
        is_double_bond_stereo_candidate(&acyclic, &acyclic_rings, BondId::new(1)),
        Ok(true)
    );

    for (size, expected) in [(7, false), (8, true)] {
        let topology = ring(size);
        let rings = fast_find_rings(&topology).unwrap();
        assert_eq!(
            should_detect_double_bond_stereo(&topology, &rings, BondId::new(0)),
            Ok(expected)
        );
    }

    let any = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::Any);
    assert_eq!(
        is_double_bond_stereo_candidate(&any, &empty_rings(&any), BondId::new(1)),
        Ok(false)
    );
    let mut crossed = acyclic.clone();
    crossed.bonds[1].set_direction(BondDirection::EitherDouble);
    assert_eq!(
        is_double_bond_stereo_candidate(&crossed, &empty_rings(&crossed), BondId::new(1)),
        Ok(false)
    );
}

#[test]
fn ring_rank_and_conformer_rows_fail_before_assignment() {
    let source = four_atom_chain(
        BondDirection::EndDownRight,
        BondDirection::EndUpRight,
        BondStereo::None,
    );
    assert_eq!(
        assign_directional_double_bond_stereo(
            source.clone(),
            &[0, 1, 2, 3],
            &RingInfo::new(RingFindType::Fast, 3, 3),
        ),
        Err(DoubleBondStereoError::RingRowCount {
            dimension: "atom",
            actual: 3,
            expected: 4,
        })
    );
    assert_eq!(
        assign_directional_double_bond_stereo(source.clone(), &[0, 1], &empty_rings(&source),),
        Err(DoubleBondStereoError::RankCount {
            actual: 2,
            atom_count: 4,
        })
    );
    let bad_rows = Conformer3D::new(9, vec![[0.0; 3]; 3], true);
    assert_eq!(
        set_double_bond_neighbor_directions(source.clone(), &empty_rings(&source), Some(&bad_rows),),
        Err(DoubleBondStereoError::InvalidConformer(
            CoordinateValidationError::RowCount {
                dimension: "3D",
                conformer: 9,
                rows: 3,
                atom_count: 4,
            }
        ))
    );
    let nonfinite = Conformer3D::new(
        7,
        vec![[0.0; 3], [0.0; 3], [f64::NAN, 0.0, 0.0], [0.0; 3]],
        true,
    );
    assert!(matches!(
        set_double_bond_neighbor_directions(
            source.clone(),
            &empty_rings(&source),
            Some(&nonfinite),
        ),
        Err(DoubleBondStereoError::InvalidConformer(
            CoordinateValidationError::NonFiniteCoordinate {
                conformer: 7,
                atom: 2,
                axis: "x",
                ..
            }
        ))
    ));
    assert_eq!(source.bonds[1].stereo(), BondStereo::None);
}

#[test]
fn topology_and_lookup_errors_retain_exact_ids() {
    let topology = topology(2, vec![bond(0, 0, 1, BondOrder::Single)]);
    assert_eq!(
        neighboring_directed_bond(&topology, AtomId::new(2)),
        Err(DoubleBondStereoError::AtomOutOfRange {
            atom: AtomId::new(2),
            atom_count: 2,
        })
    );
    assert_eq!(
        should_detect_double_bond_stereo(&topology, &empty_rings(&topology), BondId::new(2)),
        Err(DoubleBondStereoError::BondOutOfRange {
            bond: BondId::new(2),
            bond_count: 1,
        })
    );
    assert_eq!(
        find_double_bond_stereo_atoms(&topology, BondId::new(0), &[0, 1]),
        Err(DoubleBondStereoError::NotDoubleBond {
            bond: BondId::new(0),
            order: BondOrder::Single,
        })
    );
    let mut invalid = topology.clone();
    invalid.atoms[0] = invalid.atoms[0].clone().with_id(AtomId::new(1));
    assert_eq!(
        neighboring_directed_bond(&invalid, AtomId::new(0)),
        Err(DoubleBondStereoError::InvalidTopology(
            TopologyValidationError::AtomIdMismatch {
                position: 0,
                id: AtomId::new(1),
            }
        ))
    );
}

#[test]
fn neighboring_directed_bond_uses_first_adjacency_entry_and_skips_double_bonds() {
    let topology = topology(
        4,
        vec![
            directed_bond(0, 0, 1, BondOrder::Double, BondDirection::EndDownRight),
            directed_bond(1, 0, 2, BondOrder::Single, BondDirection::EndUpRight),
            directed_bond(2, 0, 3, BondOrder::Single, BondDirection::EndDownRight),
        ],
    );
    assert_eq!(
        neighboring_directed_bond(&topology, AtomId::new(0)),
        Ok(Some(BondId::new(1)))
    );
}

#[test]
fn ez_translation_and_reference_discovery_cover_existing_unique_tie_and_errors() {
    assert_eq!(translate_ez_to_cis_trans(BondStereo::E), BondStereo::Trans);
    assert_eq!(translate_ez_to_cis_trans(BondStereo::Z), BondStereo::Cis);
    assert_eq!(translate_ez_to_cis_trans(BondStereo::Any), BondStereo::Any);

    let existing = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::E);
    assert_eq!(
        find_double_bond_stereo_atoms(&existing, BondId::new(1), &[0, 1, 2, 3]),
        Ok(Some([AtomId::new(0), AtomId::new(3)]))
    );

    let no_refs = topology(
        6,
        vec![
            bond(0, 2, 3, BondOrder::Double),
            bond(1, 2, 0, BondOrder::Single),
            bond(2, 2, 1, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Single),
            bond(4, 3, 5, BondOrder::Single),
        ],
    );
    let mut no_refs_e = no_refs.clone();
    no_refs_e.bonds[0].set_stereo(BondStereo::E).unwrap();
    assert_eq!(
        find_double_bond_stereo_atoms(&no_refs_e, BondId::new(0), &[1, 9, 0, 0, 8, 2]),
        Ok(Some([AtomId::new(1), AtomId::new(4)]))
    );
    assert_eq!(
        find_double_bond_stereo_atoms(&no_refs_e, BondId::new(0), &[9, 9, 0, 0, 8, 2]),
        Ok(None)
    );

    let none = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::None);
    assert_eq!(
        find_double_bond_stereo_atoms(&none, BondId::new(1), &[0, 1, 2, 3]),
        Err(DoubleBondStereoError::UndefinedStereo {
            bond: BondId::new(1),
            stereo: BondStereo::None,
        })
    );
}

#[test]
fn chemical_reference_validation_rejects_opposite_and_non_neighbors() {
    let opposite = topology(
        4,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            stereo_bond(1, 1, 2, BondStereo::E, Some([2, 3])),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        find_double_bond_stereo_atoms(&opposite, BondId::new(1), &[0, 1, 2, 3]),
        Err(DoubleBondStereoError::StereoReferenceIsOppositeEndpoint {
            bond: BondId::new(1),
            endpoint: "begin",
            reference: AtomId::new(2),
        })
    );
    let non_neighbor = topology(
        5,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            stereo_bond(1, 1, 2, BondStereo::E, Some([4, 3])),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    assert_eq!(
        find_double_bond_stereo_atoms(&non_neighbor, BondId::new(1), &[0, 1, 2, 3, 4]),
        Err(DoubleBondStereoError::StereoReferenceNotNeighbor {
            bond: BondId::new(1),
            endpoint: "begin",
            reference: AtomId::new(4),
        })
    );
}

#[test]
fn stereo_info_pads_implicit_slots_and_flips_for_second_control() {
    let simple = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::Z);
    let info = double_bond_stereo_info(&simple, BondId::new(1)).unwrap();
    assert_eq!(
        info.controlling_atoms,
        [
            DoubleBondControl::Atom(AtomId::new(0)),
            DoubleBondControl::Implicit,
            DoubleBondControl::Atom(AtomId::new(3)),
            DoubleBondControl::Implicit,
        ]
    );
    assert_eq!(info.specified, DoubleBondStereoSpecified::Specified);
    assert_eq!(info.descriptor, Some(DoubleBondStereoDescriptor::Cis));

    let flipped = topology(
        5,
        vec![
            stereo_bond(0, 2, 3, BondStereo::Cis, Some([1, 4])),
            bond(1, 2, 0, BondOrder::Single),
            bond(2, 2, 1, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Single),
        ],
    );
    let info = double_bond_stereo_info(&flipped, BondId::new(0)).unwrap();
    assert_eq!(info.descriptor, Some(DoubleBondStereoDescriptor::Trans));
}

#[test]
fn stereo_info_unknown_precedence_and_degree_errors_are_exact() {
    let unknown = topology_with_atoms(
        vec![atom(0), unknown_atom(1), atom(2), atom(3)],
        vec![
            bond(0, 0, 1, BondOrder::Single),
            stereo_bond(1, 1, 2, BondStereo::E, Some([0, 3])),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    let info = double_bond_stereo_info(&unknown, BondId::new(1)).unwrap();
    assert_eq!(info.specified, DoubleBondStereoSpecified::Unknown);
    assert_eq!(info.descriptor, None);

    let crowded = topology(
        6,
        vec![
            stereo_bond(0, 0, 1, BondStereo::None, None),
            bond(1, 0, 2, BondOrder::Single),
            bond(2, 0, 3, BondOrder::Single),
            bond(3, 0, 4, BondOrder::Single),
            bond(4, 1, 5, BondOrder::Single),
        ],
    );
    assert_eq!(
        double_bond_stereo_info(&crowded, BondId::new(0)),
        Err(DoubleBondStereoError::InvalidEndpointDegree {
            bond: BondId::new(0),
            endpoint: "begin",
            degree: 4,
        })
    );
}

#[test]
fn reference_selection_distinguishes_normal_cx_and_stored_endpoint_order() {
    let base = topology(
        6,
        vec![
            bond(0, 2, 3, BondOrder::Double),
            bond(1, 2, 0, BondOrder::Single),
            bond(2, 2, 1, BondOrder::Single),
            bond(3, 3, 4, BondOrder::Single),
            bond(4, 3, 5, BondOrder::Single),
        ],
    );
    let normal =
        with_double_bond_stereo_reference(base.clone(), BondId::new(0), BondStereo::Any, false)
            .unwrap();
    assert_eq!(
        normal.bonds[0].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(5)])
    );
    let cx =
        with_double_bond_stereo_reference(base, BondId::new(0), BondStereo::Any, true).unwrap();
    assert_eq!(
        cx.bonds[0].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(4)])
    );

    let reversed = topology(
        4,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 2, 1, BondOrder::Double),
            bond(2, 2, 3, BondOrder::Single),
        ],
    );
    let reversed =
        with_double_bond_stereo_reference(reversed, BondId::new(1), BondStereo::Any, false)
            .unwrap();
    assert_eq!(
        reversed.bonds[1].stereo_atoms(),
        Some([AtomId::new(3), AtomId::new(0)])
    );
}

#[test]
fn directional_perception_handles_reversal_ez_and_input_immutability() {
    let source = four_atom_chain(
        BondDirection::EndDownRight,
        BondDirection::EndDownRight,
        BondStereo::None,
    );
    let snapshot = source.clone();
    let result =
        assign_directional_double_bond_stereo(source, &[0, 1, 2, 3], &empty_rings(&snapshot))
            .unwrap();
    assert!(result.assigned_any);
    assert!(!result.has_unassigned);
    assert_eq!(result.topology.bonds[1].stereo(), BondStereo::E);
    assert_eq!(
        result.topology.bonds[1].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(3)])
    );
    assert_eq!(snapshot.bonds[1].stereo(), BondStereo::None);

    let z_source = four_atom_chain(
        BondDirection::EndDownRight,
        BondDirection::EndUpRight,
        BondStereo::None,
    );
    let z = assign_directional_double_bond_stereo(
        z_source.clone(),
        &[0, 1, 2, 3],
        &empty_rings(&z_source),
    )
    .unwrap();
    assert_eq!(z.topology.bonds[1].stereo(), BondStereo::Z);
}

#[test]
fn directional_perception_reports_missing_tied_unknown_and_conflicting_neighbors() {
    let missing = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::None);
    let missing_result = assign_directional_double_bond_stereo(
        missing.clone(),
        &[0, 1, 2, 3],
        &empty_rings(&missing),
    )
    .unwrap();
    assert!(missing_result.has_unassigned);
    assert!(!missing_result.assigned_any);

    let tied = topology(
        5,
        vec![
            directed_bond(0, 0, 2, BondOrder::Single, BondDirection::EndUpRight),
            bond(1, 1, 2, BondOrder::Single),
            bond(2, 2, 3, BondOrder::Double),
            directed_bond(3, 3, 4, BondOrder::Single, BondDirection::EndDownRight),
        ],
    );
    let tied_result =
        assign_directional_double_bond_stereo(tied.clone(), &[8, 8, 0, 0, 9], &empty_rings(&tied))
            .unwrap();
    assert!(tied_result.has_unassigned);

    let unknown = topology_with_atoms(
        vec![atom(0), unknown_atom(1), atom(2), atom(3)],
        vec![
            directed_bond(0, 0, 1, BondOrder::Single, BondDirection::EndDownRight),
            bond(1, 1, 2, BondOrder::Double),
            directed_bond(2, 2, 3, BondOrder::Single, BondDirection::EndUpRight),
        ],
    );
    let unknown_result = assign_directional_double_bond_stereo(
        unknown.clone(),
        &[0, 1, 2, 3],
        &empty_rings(&unknown),
    )
    .unwrap();
    assert_eq!(unknown_result.topology.bonds[1].stereo(), BondStereo::Any);

    let conflict = topology(
        6,
        vec![
            directed_bond(0, 0, 2, BondOrder::Single, BondDirection::EndUpRight),
            directed_bond(1, 1, 2, BondOrder::Single, BondDirection::EndUpRight),
            bond(2, 2, 3, BondOrder::Double),
            directed_bond(3, 3, 4, BondOrder::Single, BondDirection::EndDownRight),
            bond(4, 3, 5, BondOrder::Single),
        ],
    );
    let conflict_result = assign_directional_double_bond_stereo(
        conflict.clone(),
        &[1, 2, 0, 0, 3, 4],
        &empty_rings(&conflict),
    )
    .unwrap();
    assert!(conflict_result.assigned_any);
    assert_eq!(conflict_result.topology.bonds[2].stereo(), BondStereo::None);
    assert_eq!(
        conflict_result.topology.bonds[0].direction(),
        BondDirection::None
    );
    assert_eq!(
        conflict_result.topology.bonds[1].direction(),
        BondDirection::None
    );
}

#[test]
fn direction_to_absolute_stereo_handles_both_direction_relations() {
    let cis = assign_double_bond_stereo_from_directions(four_atom_chain(
        BondDirection::EndDownRight,
        BondDirection::EndUpRight,
        BondStereo::None,
    ))
    .unwrap();
    assert_eq!(cis.bonds[1].stereo(), BondStereo::Cis);
    assert_eq!(
        cis.bonds[1].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(3)])
    );

    let trans = assign_double_bond_stereo_from_directions(four_atom_chain(
        BondDirection::EndDownRight,
        BondDirection::EndDownRight,
        BondStereo::None,
    ))
    .unwrap();
    assert_eq!(trans.bonds[1].stereo(), BondStereo::Trans);
}

#[test]
fn cleanup_modes_preserve_unknown_metadata_and_requested_slashes() {
    let mut topology = topology(
        4,
        vec![
            directed_bond(0, 0, 1, BondOrder::Single, BondDirection::Unknown),
            directed_bond(1, 1, 2, BondOrder::Single, BondDirection::EndDownRight),
            directed_bond(2, 2, 3, BondOrder::Double, BondDirection::EitherDouble),
        ],
    );
    let single = clear_single_bond_directions(topology.clone(), true).unwrap();
    assert_eq!(single.bonds[0].direction(), BondDirection::None);
    assert_eq!(single.bonds[0].prop("_UnknownStereo"), Some("1"));
    assert_eq!(single.bonds[1].direction(), BondDirection::EndDownRight);
    assert_eq!(single.bonds[2].direction(), BondDirection::EitherDouble);

    let all_preserve_slashes = clear_bond_directions(topology.clone(), true).unwrap();
    assert_eq!(
        all_preserve_slashes.bonds[1].direction(),
        BondDirection::EndDownRight
    );
    assert_eq!(
        all_preserve_slashes.bonds[2].direction(),
        BondDirection::None
    );
    assert_eq!(
        all_preserve_slashes.bonds[2].prop("_UnknownStereo"),
        Some("1")
    );

    topology.bonds[1].set_direction(BondDirection::BeginWedge);
    let all = clear_bond_directions(topology, false).unwrap();
    assert!(
        all.bonds
            .iter()
            .all(|bond| bond.direction() == BondDirection::None)
    );
}

#[test]
fn neighbor_direction_generation_uses_existing_stereo_and_reference_flips() {
    let source = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::E);
    let result =
        set_double_bond_neighbor_directions(source.clone(), &empty_rings(&source), None).unwrap();
    assert!(has_stereo_bond_direction(result.bonds[0].direction()));
    assert!(has_stereo_bond_direction(result.bonds[2].direction()));
    assert_eq!(result.bonds[0].direction(), result.bonds[2].direction());
    let reconstructed = assign_double_bond_stereo_from_directions(result).unwrap();
    assert_eq!(reconstructed.bonds[1].stereo(), BondStereo::Trans);

    let reversed_refs = topology(
        6,
        vec![
            bond(0, 0, 1, BondOrder::Single),
            bond(1, 4, 1, BondOrder::Single),
            stereo_bond(2, 1, 2, BondStereo::E, Some([4, 5])),
            bond(3, 2, 3, BondOrder::Single),
            bond(4, 2, 5, BondOrder::Single),
        ],
    );
    let result = set_double_bond_neighbor_directions(
        reversed_refs.clone(),
        &empty_rings(&reversed_refs),
        None,
    )
    .unwrap();
    assert!(
        result
            .bonds
            .iter()
            .any(|bond| has_stereo_bond_direction(bond.direction()))
    );
}

#[test]
fn conformer_geometry_marks_linear_unknown_and_assigns_nonlinear_directions() {
    let source = four_atom_chain(BondDirection::None, BondDirection::None, BondStereo::None);
    let linear = Conformer3D::new(
        1,
        vec![
            [-1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ],
        true,
    );
    let linear_result =
        set_double_bond_neighbor_directions(source.clone(), &empty_rings(&source), Some(&linear))
            .unwrap();
    assert_eq!(linear_result.bonds[1].stereo(), BondStereo::Any);
    assert_eq!(
        linear_result.bonds[1].stereo_atoms(),
        Some([AtomId::new(0), AtomId::new(3)])
    );

    let bent = Conformer3D::new(
        2,
        vec![
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 1.0],
        ],
        true,
    );
    let bent_result =
        set_double_bond_neighbor_directions(source.clone(), &empty_rings(&source), Some(&bent))
            .unwrap();
    assert!(has_stereo_bond_direction(bent_result.bonds[0].direction()));
    assert!(has_stereo_bond_direction(bent_result.bonds[2].direction()));
}
