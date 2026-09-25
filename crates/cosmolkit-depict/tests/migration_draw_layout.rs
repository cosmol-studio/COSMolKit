// Private owner coverage without exposing unfinished layout helpers publicly.
#[path = "../src/embedded_frag.rs"]
mod embedded_frag;
#[path = "../src/geometry.rs"]
mod geometry;
#[path = "../src/nontetrahedral.rs"]
mod nontetrahedral;
#[path = "../src/templates.rs"]
mod templates;

use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicUsize, Ordering};

use cosmolkit_core::{RingFindType, RingInfo, symmetrize_sssr_with_options_from_parts};
use cosmolkit_depict::{Compute2DCoordinatesParams, DepictError, compute_2d_coordinates};
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, BondStereo, ChiralTag, Element,
    TopologyBlock,
};

#[test]
fn entry_output_defaults_cover_empty_single_chain_and_disconnected_rows() {
    let empty = topology(vec![], &[]);
    assert!(
        compute_2d_coordinates(&empty, &Default::default())
            .unwrap()
            .coordinates()
            .is_empty()
    );

    let single = topology(vec![AtomSpec::new(Element::C)], &[]);
    assert_eq!(
        compute_2d_coordinates(&single, &Default::default())
            .unwrap()
            .coordinates(),
        [[0.0, 0.0]]
    );

    let chain = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let source = chain.clone();
    let chain_xy = compute_2d_coordinates(&chain, &Default::default()).unwrap();
    assert_eq!(chain, source, "detached depiction must not mutate topology");
    assert_eq!(chain_xy.id(), 0);
    assert_eq!(
        chain_xy.coordinates(),
        [
            [0.0, 0.0],
            [1.299038105676658, 0.7499999999999998],
            [2.598076211353316, -6.661338147750939e-16],
        ]
    );

    let disconnected = topology(vec![AtomSpec::new(Element::C); 2], &[]);
    assert_eq!(
        compute_2d_coordinates(&disconnected, &Default::default())
            .unwrap()
            .coordinates(),
        [[0.0, 0.0], [1.0, 0.0]]
    );
}

#[test]
fn entry_output_coordinate_constraints_cover_empty_single_many_and_invalid() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let empty = compute_2d_coordinates(&graph, &Default::default()).unwrap();
    let explicit_empty = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            coordinate_map: BTreeMap::new(),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(explicit_empty, empty);

    let singleton = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            coordinate_map: BTreeMap::from([(1, [7.0, -3.0])]),
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(singleton.coordinates()[1], [7.0, -3.0]);

    let fixed = BTreeMap::from([(0, [-0.0, 2.0]), (1, [1.5, 2.0])]);
    let multiple = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            coordinate_map: fixed.clone(),
            ..Default::default()
        },
    )
    .unwrap();
    for (atom, expected) in fixed {
        assert_eq!(
            multiple.coordinates()[atom][0].to_bits(),
            expected[0].to_bits()
        );
        assert_eq!(
            multiple.coordinates()[atom][1].to_bits(),
            expected[1].to_bits()
        );
    }

    let error = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            coordinate_map: BTreeMap::from([(3, [0.0, 0.0])]),
            ..Default::default()
        },
    )
    .unwrap_err();
    assert!(matches!(
        error,
        DepictError::Fragment(
            cosmolkit_depict::Coordinate2DLayoutError::AtomIndexOutOfRange {
                atom: 3,
                atom_count: 3
            }
        )
    ));
}

#[test]
fn detached_typed_error_chain_preserves_topology_source() {
    let mut invalid = TopologyBlock::default();
    invalid
        .atoms
        .push(Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)));

    let error = compute_2d_coordinates(&invalid, &Default::default()).unwrap_err();
    let source = std::error::Error::source(&error)
        .expect("depiction topology errors must retain the typed topology source");
    assert!(matches!(
        source.downcast_ref::<cosmolkit_model::TopologyValidationError>(),
        Some(cosmolkit_model::TopologyValidationError::AtomIdMismatch {
            position: 0,
            id
        }) if *id == AtomId::new(1)
    ));
}

#[test]
fn detached_typed_error_chain_preserves_layout_fields() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let error = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            coordinate_map: BTreeMap::from([(3, [0.0, 0.0])]),
            ..Default::default()
        },
    )
    .unwrap_err();
    let source = std::error::Error::source(&error)
        .expect("layout errors must retain their public typed projection");
    assert!(matches!(
        source.downcast_ref::<cosmolkit_depict::Coordinate2DLayoutError>(),
        Some(
            cosmolkit_depict::Coordinate2DLayoutError::AtomIndexOutOfRange {
                atom: 3,
                atom_count: 3
            }
        )
    ));
}

#[test]
fn detached_typed_error_chain_preserves_sampling_safety_boundary() {
    // This is the deterministic CK boundary for an upstream unwritten
    // sampling-cost slot, not an RDKit-parity result or a replacement value.
    let disconnected = topology(vec![AtomSpec::new(Element::C); 2], &[]);
    let error = compute_2d_coordinates(
        &disconnected,
        &Compute2DCoordinatesParams {
            flips_per_sample: 1,
            samples: 1,
            sample_seed: 7,
            ..Default::default()
        },
    )
    .unwrap_err();
    let source = std::error::Error::source(&error)
        .expect("sampling boundary must be distinguishable as a typed CK error");
    assert!(matches!(
        source.downcast_ref::<cosmolkit_depict::Coordinate2DLayoutError>(),
        Some(
            cosmolkit_depict::Coordinate2DLayoutError::UndefinedSamplingDistance {
                first: 1,
                second: 0
            }
        )
    ));
}

#[test]
fn entry_output_options_cover_orientation_route_sampling_templates_and_clear_carrier() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (1, 2), (2, 3), (3, 4)],
    );
    let defaults = compute_2d_coordinates(&graph, &Default::default()).unwrap();
    let canonical = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            canonical_orientation: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_ne!(canonical.coordinates(), defaults.coordinates());
    for bond in &graph.bonds {
        let [left, right] = [bond.begin().index(), bond.end().index()];
        let distance = |rows: &[[f64; 2]]| {
            (rows[left][0] - rows[right][0]).hypot(rows[left][1] - rows[right][1])
        };
        assert!(
            (distance(canonical.coordinates()) - distance(defaults.coordinates())).abs() < 1e-12
        );
    }

    let routed = compute_2d_coordinates(
        &graph,
        &Compute2DCoordinatesParams {
            clear_existing_2d: false,
            force_rdkit: true,
            use_ring_templates: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        routed, defaults,
        "detached clear carrier and fixed RDKit route do not change XY"
    );

    let sampled_params = Compute2DCoordinatesParams {
        flips_per_sample: 1,
        samples: 2,
        sample_seed: 7,
        permute_degree_four: true,
        ..Default::default()
    };
    let sampled_a = compute_2d_coordinates(&graph, &sampled_params).unwrap();
    let sampled_b = compute_2d_coordinates(&graph, &sampled_params).unwrap();
    assert_eq!(
        sampled_a, sampled_b,
        "positive seed resets the shared source RNG"
    );
}

#[test]
fn entry_output_ring_and_stereo_inputs_complete_atom_ordered_finite_rows() {
    let ring = topology(
        vec![AtomSpec::new(Element::C); 6],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
    );
    let ring_xy = compute_2d_coordinates(
        &ring,
        &Compute2DCoordinatesParams {
            use_ring_templates: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(ring_xy.coordinates().len(), 6);
    assert!(
        ring_xy
            .coordinates()
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );

    let stereo = cis_trans_graph(
        4,
        &[
            (0, 1, BondOrder::Single, BondStereo::None, None),
            (1, 2, BondOrder::Double, BondStereo::E, Some((0, 3))),
            (2, 3, BondOrder::Single, BondStereo::None, None),
        ],
    );
    let stereo_xy = compute_2d_coordinates(&stereo, &Default::default()).unwrap();
    assert_eq!(stereo_xy.coordinates().len(), 4);
    assert!(
        stereo_xy
            .coordinates()
            .iter()
            .flatten()
            .all(|value| value.is_finite())
    );
    let left = stereo_xy.coordinates()[1];
    let right = stereo_xy.coordinates()[2];
    assert!(((left[0] - right[0]).hypot(left[1] - right[1]) - geometry::BOND_LEN).abs() < 1e-12);
}

#[test]
fn sampling_boundary_public_partial_fragments_do_not_gain_a_fallback_contract() {
    // Fixed RDKit 2026.03.1 happens to return coordinates for C.C and reaches
    // map::at for CC.CC, but both first read unwritten sampling-cost storage.
    // Those observations are retained in DRAW-layout.md, not promoted to a
    // portable expected value. CK must stop at its deterministic safety
    // boundary instead of zero-filling, filtering or skipping the sampler.
    for disconnected in [
        topology(vec![AtomSpec::new(Element::C); 2], &[]),
        topology(vec![AtomSpec::new(Element::C); 4], &[(0, 1), (2, 3)]),
    ] {
        let error = compute_2d_coordinates(
            &disconnected,
            &Compute2DCoordinatesParams {
                flips_per_sample: 1,
                samples: 1,
                sample_seed: 7,
                ..Default::default()
            },
        )
        .unwrap_err();
        assert!(matches!(
            error,
            DepictError::Fragment(
                cosmolkit_depict::Coordinate2DLayoutError::UndefinedSamplingDistance { .. }
            )
        ));
    }
}

fn collision_fragment<'a>(
    graph: &'a TopologyBlock,
    rings: &'a RingInfo,
    points: &[(usize, [f64; 2])],
) -> embedded_frag::EmbeddedFrag<'a> {
    let map = points.iter().copied().collect::<BTreeMap<_, _>>();
    embedded_frag::EmbeddedFrag::from_coord_map(graph, rings, &map).unwrap()
}

#[test]
fn random_engine_minstd_seed_and_zero_seed_continuation() {
    // RDGeneral/utils.cpp starts at 42 and reseeds only for a positive seed.
    let mut rng = embedded_frag::DepictRandom::with_seed(42);
    assert_eq!(rng.next_raw(), 2_027_382);
    rng.reset_if_positive(0);
    assert_eq!(rng.next_raw(), 1_226_992_407);
    rng.reset_if_positive(-1);
    assert_eq!(rng.next_raw(), 551_494_037);
    rng.reset_if_positive(7);
    assert_eq!(rng.next_raw(), 337_897);
    let mut modulo_zero = embedded_frag::DepictRandom::with_seed(2_147_483_647);
    assert_eq!(modulo_zero.next_raw(), 48_271);
}

#[test]
fn random_engine_boost_integer_bucket_and_rejection_order() {
    // Boost 1.85 uniform_int distribution discards a bucket result above the
    // requested range instead of introducing modulo bias.
    let mut rng = embedded_frag::DepictRandom::with_seed(7);
    assert_eq!(rng.next_index(5), 0);
    assert_eq!(rng.next_index(5), 2);
    assert_eq!(rng.next_index(5), 1);
    let mut rejection = embedded_frag::DepictRandom::with_seed(247_665_088);
    assert_eq!(rejection.next_index(5), 4);
    assert_eq!(rejection.next_raw(), 1_964_877_853);
    let mut wide = embedded_frag::DepictRandom::with_seed(7);
    assert_eq!(wide.next_index(i32::MAX as usize), 449_829_613);
    // The first two draws form an out-of-range concatenated result and are
    // rejected; draws three and four form the accepted value.
    assert_eq!(wide.next_raw(), 1_665_781_405);
}

#[test]
fn random_engine_boost_real_draw_shares_integer_stream() {
    let mut rng = embedded_frag::DepictRandom::with_seed(7);
    assert_eq!(
        rng.next_real().to_bits(),
        0.000_157_345_086_482_674_9_f64.to_bits()
    );
    assert_eq!(rng.next_index(5), 2);
    assert_eq!(
        rng.next_real().to_bits(),
        0.209_468_236_853_804_67_f64.to_bits()
    );
}

fn complete_sampling_fragment<'a>(
    graph: &'a TopologyBlock,
    rings: &'a RingInfo,
    points: &[(usize, [f64; 2])],
) -> embedded_frag::EmbeddedFrag<'a> {
    let mut fragment = collision_fragment(graph, rings, points);
    for atom in fragment.atoms.values_mut() {
        atom.fixed = false;
    }
    fragment
}

fn sampling_locations(fragment: &embedded_frag::EmbeddedFrag<'_>) -> Vec<(usize, [f64; 2])> {
    fragment
        .atoms
        .iter()
        .map(|(&atom_id, atom)| (atom_id, atom.loc))
        .collect()
}

#[test]
fn sampling_boundary_defined_zero_or_one_positive_count_performs_no_draw_or_move() {
    // EmbeddedFrag.cpp uses nPerSample=min(nt,nBondsPerSample), and its outer
    // loop is bounded by nSamples. Either zero therefore performs no draw.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let points = [
        (0, [-1.0, 1.0]),
        (1, [0.0, 0.0]),
        (2, [1.0, 0.0]),
        (3, [2.0, 1.0]),
    ];
    for (bonds_per_sample, samples) in [(0, 0), (1, 0), (0, 2)] {
        let mut fragment = complete_sampling_fragment(&graph, &rings, &points);
        let before = fragment.atoms.clone();
        fragment
            .random_sample_flips_and_permutations(
                bonds_per_sample,
                samples,
                20_000,
                None,
                0.0,
                false,
            )
            .unwrap();
        assert_eq!(fragment.atoms, before);
    }
}

#[test]
fn sampling_boundary_defined_connected_seed_uses_whole_bond_table_order_and_best_cost() {
    // With Boost seed 20000 and three candidates, the first uniform draw is
    // index 1. Source bond-table order therefore selects the central bond.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut fragment = complete_sampling_fragment(
        &graph,
        &rings,
        &[
            (0, [-1.0, 1.0]),
            (1, [0.0, 0.0]),
            (2, [1.0, 0.0]),
            (3, [2.0, 1.0]),
        ],
    );
    fragment
        .random_sample_flips_and_permutations(1, 1, 20_000, None, 0.0, false)
        .unwrap();
    assert_eq!(
        sampling_locations(&fragment),
        [
            (0, [-1.0, 1.0]),
            (1, [0.0, 0.0]),
            (2, [1.0, 0.0]),
            (3, [1.999_999_999_999_999_8, -1.0]),
        ]
    );
    assert!(fragment.atoms[&0].ccw);
    assert!(fragment.atoms[&1].ccw);
    assert!(!fragment.atoms[&2].ccw);
    assert!(!fragment.atoms[&3].ccw);
}

#[test]
fn random_sampling_degree_four_seed_uses_source_pair_and_coordinate_only_restore() {
    // With Boost seed 40000, uniform_int(0,4) selects the degree-four row and
    // the following real draw is <0.5, choosing source pair zero. The cardinal
    // geometry has equal cost after swapping those arms, so source restores
    // only coordinates while the two EmbeddedAtom ccw flags remain reflected.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (0, 2), (0, 3), (0, 4)],
    );
    let rings = source_rings(&graph);
    let points = [
        (0, [0.0, 0.0]),
        (1, [1.0, 0.0]),
        (2, [0.0, 1.0]),
        (3, [-1.0, 0.0]),
        (4, [0.0, -1.0]),
    ];
    let mut fragment = complete_sampling_fragment(&graph, &rings, &points);
    fragment
        .random_sample_flips_and_permutations(1, 1, 40_000, None, 0.0, true)
        .unwrap();
    assert_eq!(sampling_locations(&fragment), points);
    assert!(fragment.atoms[&0].ccw);
    assert!(!fragment.atoms[&1].ccw);
    assert!(!fragment.atoms[&2].ccw);
    assert!(fragment.atoms[&3].ccw);
    assert!(fragment.atoms[&4].ccw);
}

#[test]
fn sampling_boundary_ck_safety_error_reports_first_unwritten_pair_without_mutation() {
    // The pinned C++ allocation leaves this pair uninitialized. Rust must not
    // read it, assign a replacement number, or proceed to candidate sampling.
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = source_rings(&graph);
    let mut fragment =
        complete_sampling_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (1, [1.5, 0.0])]);
    let before = fragment.atoms.clone();
    assert_eq!(
        fragment.random_sample_flips_and_permutations(1, 1, 7, None, 0.0, false),
        Err(embedded_frag::FragmentError::UndefinedSamplingDistance {
            first: 2,
            second: 0,
        })
    );
    assert_eq!(fragment.atoms, before);
}

#[test]
fn sampling_boundary_source_at_lookup_error_remains_distinct_from_ck_safety_error() {
    // Whole-molecule candidate construction includes bond 1, but its atom 2
    // row is absent from this fragment. flipAboutBond fails at d_eatoms.at,
    // rather than filtering the candidate or fabricating a coordinate.
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = source_rings(&graph);
    let mut fragment =
        complete_sampling_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (1, [1.5, 0.0])]);
    let before = fragment.atoms.clone();
    assert_eq!(
        fragment.flip_about_bond(1, true),
        Err(embedded_frag::FragmentError::AtomNotEmbedded { atom: 2 })
    );
    assert_eq!(fragment.atoms, before);
}

#[test]
fn random_partial_flip_fixed_scan_source_default_inserts_recursive_rows() {
    // In the fixed-row branch the pinned source uses d_eatoms[endAtomId], so
    // the missing recursive-side atom is default inserted before side choice.
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = source_rings(&graph);
    let mut fragment = collision_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (1, [1.5, 0.0])]);
    fragment.atoms.get_mut(&1).unwrap().fixed = false;
    fragment.flip_about_bond(0, true).unwrap();
    assert_eq!(
        fragment.atoms.keys().copied().collect::<Vec<_>>(),
        vec![0, 1, 2]
    );
    let inserted = &fragment.atoms[&2];
    assert_eq!(inserted.aid, 0);
    assert_eq!(inserted.loc, [0.0, 0.0]);
    assert_eq!(inserted.normal, [0.0, 0.0]);
    assert_eq!(inserted.angle, -1.0);
    assert_eq!(inserted.nbr1, None);
    assert_eq!(inserted.nbr2, None);
    assert_eq!(inserted.cis_trans_nbr, None);
    assert!(inserted.ccw);
    assert_eq!(inserted.rot_dir, 0);
    assert!(inserted.neighs.is_empty());
    assert_eq!(inserted.density, -1.0);
    assert!(!inserted.fixed);
}

#[test]
fn random_partial_permutation_source_default_inserts_deeper_side_rows() {
    // permuteBonds uses operator[] in both reflection loops. All four direct
    // neighbor rows exist, but recursive atom 5 is absent and is inserted.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 6],
        &[(0, 1), (0, 2), (0, 3), (0, 4), (1, 5)],
    );
    let rings = source_rings(&graph);
    let mut fragment = complete_sampling_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.0, 1.0]),
            (3, [-1.0, 0.0]),
            (4, [0.0, -1.0]),
        ],
    );
    fragment.permute_bonds(0, 1, 2).unwrap();
    let inserted = &fragment.atoms[&5];
    assert_eq!(inserted.aid, 0);
    assert_eq!(inserted.loc, [8.326_672_684_688_674e-17, 0.0]);
    assert_eq!(inserted.normal, [0.0, 0.0]);
    assert!(!inserted.ccw);
    assert!(!inserted.fixed);
}

#[test]
fn orientation_shift_absent_and_empty_coordinate_maps_canonicalize_exact_xy() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let rings = source_rings(&graph);
    let input = geometry::PointMap::from([(0, [2.0, 2.0]), (1, [4.0, 4.0])]);
    // Pinned code normalizes [4,4] component-wise and then sums two
    // transformed products; independently evaluating sqrt(2) rounds one bit
    // differently from that source operation order.
    let expected = 1.414_213_562_373_095_f64;

    for coordinate_map_size in [None, Some(0)] {
        let mut fragments =
            vec![embedded_frag::EmbeddedFrag::from_coord_map(&graph, &rings, &input).unwrap()];
        embedded_frag::orient_and_shift_fragments(&mut fragments, true, coordinate_map_size);
        assert_eq!(fragments[0].atoms[&0].loc, [-expected, 0.0]);
        assert_eq!(fragments[0].atoms[&1].loc, [expected, 0.0]);
    }
}

#[test]
fn orientation_shift_nonempty_coordinate_map_suppresses_canonical_rotation() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let rings = source_rings(&graph);
    let input = geometry::PointMap::from([(0, [2.0, 2.0]), (1, [4.0, 4.0])]);
    let mut fragments =
        vec![embedded_frag::EmbeddedFrag::from_coord_map(&graph, &rings, &input).unwrap()];

    embedded_frag::orient_and_shift_fragments(&mut fragments, true, Some(2));

    assert_eq!(fragments[0].atoms[&0].loc, [2.0, 2.0]);
    assert_eq!(fragments[0].atoms[&1].loc, [4.0, 4.0]);
}

#[test]
fn orientation_shift_disconnected_fragments_follow_source_axis_and_spacing() {
    let graph = topology(vec![AtomSpec::new(Element::C); 6], &[]);
    let rings = source_rings(&graph);
    let mut fragments = vec![
        embedded_frag::EmbeddedFrag::from_coord_map(
            &graph,
            &rings,
            &geometry::PointMap::from([(0, [0.0, -1.0]), (1, [0.0, 1.0])]),
        )
        .unwrap(),
        embedded_frag::EmbeddedFrag::from_coord_map(
            &graph,
            &rings,
            &geometry::PointMap::from([(2, [-2.0, 0.0]), (3, [2.0, 0.0])]),
        )
        .unwrap(),
        embedded_frag::EmbeddedFrag::from_coord_map(
            &graph,
            &rings,
            &geometry::PointMap::from([(4, [0.0, -0.5]), (5, [0.0, 0.5])]),
        )
        .unwrap(),
    ];

    embedded_frag::orient_and_shift_fragments(&mut fragments, false, None);

    assert_eq!(fragments[0].atoms[&0].loc, [0.0, -1.0]);
    assert_eq!(fragments[0].atoms[&1].loc, [0.0, 1.0]);
    assert_eq!(fragments[1].atoms[&2].loc, [1.0, 0.0]);
    assert_eq!(fragments[1].atoms[&3].loc, [5.0, 0.0]);
    assert_eq!(fragments[2].atoms[&4].loc, [0.0, 2.0]);
    assert_eq!(fragments[2].atoms[&5].loc, [0.0, 3.0]);
}

#[test]
fn constraint_semantics_zero_one_many_seed_and_index_validation() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = source_rings(&graph);
    assert!(
        embedded_frag::seed_coordinate_constraints(&graph, &rings, None)
            .unwrap()
            .is_none()
    );
    assert!(
        embedded_frag::seed_coordinate_constraints(
            &graph,
            &rings,
            Some(&geometry::PointMap::new()),
        )
        .unwrap()
        .is_none()
    );
    let singleton = geometry::PointMap::from([(1, [8.0, 9.0])]);
    assert!(
        embedded_frag::seed_coordinate_constraints(&graph, &rings, Some(&singleton))
            .unwrap()
            .is_none()
    );

    let many = geometry::PointMap::from([(2, [-0.0, 4.0]), (0, [2.0, 3.0])]);
    let seeded = embedded_frag::seed_coordinate_constraints(&graph, &rings, Some(&many))
        .unwrap()
        .expect("multi-row coordinate map seeds the source fragment");
    assert_eq!(seeded.atoms.keys().copied().collect::<Vec<_>>(), vec![0, 2]);
    assert_eq!(seeded.atoms[&0].loc, [2.0, 3.0]);
    assert_eq!(seeded.atoms[&2].loc[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(seeded.atoms[&2].loc[1], 4.0);
    assert!(seeded.atoms.values().all(|atom| atom.fixed));

    let invalid = geometry::PointMap::from([(3, [0.0, 0.0])]);
    assert_eq!(
        embedded_frag::seed_coordinate_constraints(&graph, &rings, Some(&invalid)).unwrap_err(),
        embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 3,
            atom_count: 3,
        }
    );
}

#[test]
fn constraint_semantics_singleton_translates_all_rows_and_repeated_key_uses_last_value() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[]);
    let rings = source_rings(&graph);
    let points = geometry::PointMap::from([(0, [1.0, 2.0]), (1, [3.0, 4.0]), (2, [-1.0, -2.0])]);
    let mut fragments =
        vec![embedded_frag::EmbeddedFrag::from_coord_map(&graph, &rings, &points).unwrap()];
    let mut singleton = geometry::PointMap::new();
    singleton.insert(1, [90.0, 90.0]);
    singleton.insert(1, [8.0, 9.0]);
    embedded_frag::translate_single_coordinate_constraint(&graph, &mut fragments, Some(&singleton))
        .unwrap();
    assert_eq!(fragments[0].atoms[&0].loc, [6.0, 7.0]);
    assert_eq!(fragments[0].atoms[&1].loc, [8.0, 9.0]);
    assert_eq!(fragments[0].atoms[&2].loc, [4.0, 3.0]);

    let after_singleton = fragments[0].atoms.clone();
    embedded_frag::translate_single_coordinate_constraint(&graph, &mut fragments, None).unwrap();
    embedded_frag::translate_single_coordinate_constraint(
        &graph,
        &mut fragments,
        Some(&geometry::PointMap::new()),
    )
    .unwrap();
    let multiple = geometry::PointMap::from([(0, [0.0, 0.0]), (2, [2.0, 0.0])]);
    embedded_frag::translate_single_coordinate_constraint(&graph, &mut fragments, Some(&multiple))
        .unwrap();
    assert_eq!(fragments[0].atoms, after_singleton);

    let invalid = geometry::PointMap::from([(4, [0.0, 0.0])]);
    assert_eq!(
        embedded_frag::translate_single_coordinate_constraint(
            &graph,
            &mut fragments,
            Some(&invalid),
        ),
        Err(embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 4,
            atom_count: 3,
        })
    );
    assert_eq!(fragments[0].atoms, after_singleton);
}

#[test]
fn constraint_semantics_ring_template_threshold_and_fixed_seed_coordinates() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)],
    );
    let rings = source_rings(&graph);
    let rows = "[*]1-[*]2-[*]-1-[*]-2 |(10,0,;11,0,;11,1,;10,1,)|\n";
    let mut templates = template_match_registry(rows);

    let one = geometry::PointMap::from([(0, [99.0, 99.0])]);
    assert!(
        embedded_frag::seed_coordinate_constraints(&graph, &rings, Some(&one))
            .unwrap()
            .is_none()
    );
    let allowed =
        embedded_frag::embed_fused_systems(&graph, &rings, Some(&one), true, &mut templates)
            .unwrap();
    assert!(allowed[0].atoms.values().all(|atom| atom.fixed));
    assert_eq!(allowed[0].atoms[&0].loc, [10.0, 0.0]);

    let two = geometry::PointMap::from([(0, [-0.0, 5.0]), (1, [2.0, 5.0])]);
    let mut seeded = embedded_frag::seed_coordinate_constraints(&graph, &rings, Some(&two))
        .unwrap()
        .expect("two constraints seed before ring-system merge");
    let mut disabled =
        embedded_frag::embed_fused_systems(&graph, &rings, Some(&two), true, &mut templates)
            .unwrap();
    assert!(disabled[0].atoms.values().all(|atom| !atom.fixed));
    seeded.merge_frags_with_common(&mut disabled).unwrap();
    assert!(disabled.is_empty());
    assert_eq!(seeded.atoms[&0].loc[0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!(seeded.atoms[&0].loc[1], 5.0);
    assert_eq!(seeded.atoms[&1].loc, [2.0, 5.0]);
    assert!(seeded.atoms[&0].fixed);
    assert!(seeded.atoms[&1].fixed);
}

#[test]
fn collision_discovery_atom_threshold_hetero_factor_and_density() {
    // EmbeddedFrag.cpp::findCollisions: strict 0.70² threshold, independent
    // 1.3 hetero factors, and 1/d² density before factor scaling.
    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::N),
            AtomSpec::new(Element::O),
        ],
        &[],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [0.71, 0.0]),
            (2, [0.0, 0.79]),
            (3, [3.0, 3.0]),
        ],
    );
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, false), [(2, 0)]);
    let density_sum = frag.atoms.values().map(|row| row.density).sum::<f64>();
    assert_eq!(frag.total_density().to_bits(), density_sum.to_bits());
    assert!(frag.atoms[&0].density > 1.0 / (0.71 * 0.71));
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[]);
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (1, [0.0, 0.0])]);
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, false), [(1, 0)]);
    assert_eq!(frag.total_density(), 2000.0);
}

#[test]
fn collision_discovery_bond_crossing_selects_source_first_distance_tie() {
    // EmbeddedFrag.cpp::_findClosestPair compares (beg1,beg2),
    // (beg1,end2), (end1,beg2), (end1,end2) with strict less-than ties.
    let graph = topology(vec![AtomSpec::new(Element::C); 4], &[(0, 1), (2, 3)]);
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [-1.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.0, -1.0]),
            (3, [0.0, 1.0]),
        ],
    );
    let distance = frag.collision_distance_matrix().unwrap();
    assert!(frag.find_collisions(&distance, false).is_empty());
    assert_eq!(frag.find_collisions(&distance, true), [(0, 2)]);
}

#[test]
fn collision_discovery_rotatable_path_uses_interior_edges_only() {
    // DepictUtils.cpp::getRotatableBonds excludes endpoint edges, not
    // interior double bonds merely because of their bond order.
    let graph = cis_trans_graph(
        5,
        &[
            (0, 1, BondOrder::Single, BondStereo::None, None),
            (1, 2, BondOrder::Double, BondStereo::None, None),
            (2, 3, BondOrder::Single, BondStereo::None, None),
            (3, 4, BondOrder::Single, BondStereo::None, None),
        ],
    );
    let rings = source_rings(&graph);
    let frag = collision_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (4, [1.0, 0.0])]);
    assert_eq!(frag.rotatable_bonds_on_shortest_path(0, 4).unwrap(), [1, 2]);
    assert!(
        frag.rotatable_bonds_on_shortest_path(0, 2)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn collision_bond_flip_reflects_source_end_side_and_reverses_on_repeat() {
    // EmbeddedFrag.cpp::flipAboutBond chooses the smaller side, with an equal
    // split staying on the selected end. Reflection is an involution.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [-1.0, 1.0]),
            (1, [0.0, 0.0]),
            (2, [1.0, 0.0]),
            (3, [2.0, 1.0]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    frag.flip_about_bond(1, true).unwrap();
    assert_eq!(frag.atoms[&0].loc, [-1.0, 1.0]);
    assert!((frag.atoms[&3].loc[0] - 2.0).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[1] + 1.0).abs() < 1e-12);
    frag.flip_about_bond(1, true).unwrap();
    assert!((frag.atoms[&3].loc[0] - 2.0).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[1] - 1.0).abs() < 1e-12);
    frag.flip_about_bond(1, false).unwrap();
    assert!((frag.atoms[&0].loc[0] + 1.0).abs() < 1e-12);
    assert!((frag.atoms[&0].loc[1] + 1.0).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[0] - 2.0).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[1] - 1.0).abs() < 1e-12);
}

#[test]
fn collision_bond_flip_preserves_source_fixed_side_guard_and_ring_precondition() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [-1.0, 1.0]),
            (1, [0.0, 0.0]),
            (2, [1.0, 0.0]),
            (3, [2.0, 1.0]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    frag.atoms.get_mut(&3).unwrap().fixed = true;
    let before = frag
        .atoms
        .iter()
        .map(|(&id, atom)| (id, atom.loc))
        .collect::<Vec<_>>();
    frag.flip_about_bond(1, true).unwrap();
    assert_eq!(
        frag.atoms
            .iter()
            .map(|(&id, atom)| (id, atom.loc))
            .collect::<Vec<_>>(),
        before
    );
    let ring = topology(
        vec![AtomSpec::new(Element::C); 3],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let rings = source_rings(&ring);
    let mut frag = collision_fragment(
        &ring,
        &rings,
        &[(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])],
    );
    assert!(matches!(
        frag.flip_about_bond(0, true),
        Err(embedded_frag::FragmentError::CollisionBondInvalid { bond: 0 })
    ));
}

#[test]
fn collision_bond_flip_bounded_unresolved_disconnected_pair_keeps_coordinates() {
    // With no shortest path the source loop repeats at most 15 times but
    // cannot select a rotatable bond; it does not fabricate a connection.
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[]);
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (1, [0.1, 0.0])]);
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    frag.remove_collisions_bond_flip().unwrap();
    assert_eq!(frag.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(frag.atoms[&1].loc, [0.1, 0.0]);
}

#[test]
fn collision_bond_flip_repairs_first_collision_with_source_path_order() {
    // EmbeddedFrag.cpp::removeCollisionsBondFlip starts with the first pair
    // (larger embedded ID first), and accepts a strict collision-count drop.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (1, 2), (2, 3), (3, 4)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [2.0, 0.0]),
            (3, [2.0, 1.0]),
            (4, [0.1, 0.1]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, true), [(4, 0)]);
    frag.remove_collisions_bond_flip().unwrap();
    assert!(frag.find_collisions(&distance, true).is_empty());
    assert!((frag.atoms[&4].loc[0] - 3.9).abs() < 1e-12);
    assert!((frag.atoms[&4].loc[1] - 0.1).abs() < 1e-12);
}

#[test]
fn collision_bond_flip_path_excludes_source_stereo_and_ring_bonds() {
    let graph = cis_trans_graph(
        6,
        &[
            (0, 1, BondOrder::Single, BondStereo::None, None),
            (1, 2, BondOrder::Single, BondStereo::None, None),
            (2, 3, BondOrder::Double, BondStereo::E, Some((1, 4))),
            (3, 4, BondOrder::Single, BondStereo::None, None),
            (4, 5, BondOrder::Single, BondStereo::None, None),
        ],
    );
    let rings = source_rings(&graph);
    let frag = collision_fragment(&graph, &rings, &[(0, [0.0, 0.0]), (5, [1.0, 0.0])]);
    assert_eq!(frag.rotatable_bonds_on_shortest_path(0, 5).unwrap(), [1, 3]);

    let ring = topology(
        vec![AtomSpec::new(Element::C); 6],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 1), (4, 5)],
    );
    let rings = source_rings(&ring);
    let frag = collision_fragment(&ring, &rings, &[(0, [0.0, 0.0]), (5, [1.0, 0.0])]);
    assert!(
        frag.rotatable_bonds_on_shortest_path(0, 5)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn collision_bond_flip_near_tie_obeys_strict_computed_density() {
    // Geometric symmetry does not imply floating-point equality: reflectPoint
    // uses two transforms, so this near-tie can have a strictly lower computed
    // density. The source accepts that strict improvement, not an epsilon tie.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [2.0, 0.0]),
            (3, [0.1, 0.1]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, true), [(3, 0)]);
    let before_density = frag.total_density();
    frag.remove_collisions_bond_flip().unwrap();
    assert!(frag.total_density() < before_density);
    assert_eq!(frag.find_collisions(&distance, true), [(3, 0)]);
}

#[test]
fn collision_open_angle_two_degree_one_rows_rotate_around_own_neighbors() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 1.0]),
            (1, [0.0, 0.0]),
            (2, [1.0, 0.0]),
            (3, [1.0, 1.0]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = frag.collision_distance_matrix().unwrap();
    frag.open_angles(&distance, 0, 3).unwrap();
    let angle = 0.1222_f64;
    assert!((frag.atoms[&0].loc[0] + angle.sin()).abs() < 1e-12);
    assert!((frag.atoms[&0].loc[1] - angle.cos()).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[0] - (1.0 + angle.sin())).abs() < 1e-12);
    assert!((frag.atoms[&3].loc[1] - angle.cos()).abs() < 1e-12);
    assert_eq!(frag.atoms[&1].loc, [0.0, 0.0]);
    assert_eq!(frag.atoms[&2].loc, [1.0, 0.0]);
}

#[test]
fn collision_open_angle_one_terminal_row_uses_source_case_two_and_three() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (1, 2), (2, 3), (2, 4)],
    );
    let rings = source_rings(&graph);
    let points = [
        (0, [0.0, -1.0]),
        (1, [0.0, 0.0]),
        (2, [1.0, 0.0]),
        (3, [2.0, 1.0]),
        (4, [2.0, -1.0]),
    ];
    let mut case_two = collision_fragment(&graph, &rings, &points);
    for atom in case_two.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = case_two.collision_distance_matrix().unwrap();
    case_two.open_angles(&distance, 0, 2).unwrap();
    let expected = geometry::Transform2D::around([0.0, 0.0], 0.2444).transform_point([0.0, -1.0]);
    assert!((case_two.atoms[&0].loc[0] - expected[0]).abs() < 1e-12);
    assert!((case_two.atoms[&0].loc[1] - expected[1]).abs() < 1e-12);
    assert_eq!(case_two.atoms[&2].loc, [1.0, 0.0]);

    let mut case_three = collision_fragment(&graph, &rings, &points);
    for atom in case_three.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = case_three.collision_distance_matrix().unwrap();
    case_three.open_angles(&distance, 2, 0).unwrap();
    let expected = geometry::Transform2D::around([0.0, 0.0], -0.2444).transform_point([0.0, -1.0]);
    assert!((case_three.atoms[&0].loc[0] - expected[0]).abs() < 1e-12);
    assert!((case_three.atoms[&0].loc[1] - expected[1]).abs() < 1e-12);
    assert_eq!(case_three.atoms[&2].loc, [1.0, 0.0]);
}

#[test]
fn collision_open_angle_fixed_and_nonterminal_guards_preserve_rows() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let points = [
        (0, [0.0, 1.0]),
        (1, [0.0, 0.0]),
        (2, [1.0, 0.0]),
        (3, [1.0, 1.0]),
    ];
    let mut frag = collision_fragment(&graph, &rings, &points);
    let distance = frag.collision_distance_matrix().unwrap();
    let before = frag
        .atoms
        .iter()
        .map(|(&id, atom)| (id, atom.loc))
        .collect::<Vec<_>>();
    frag.open_angles(&distance, 0, 3).unwrap();
    frag.open_angles(&distance, 1, 2).unwrap();
    assert_eq!(
        frag.atoms
            .iter()
            .map(|(&id, atom)| (id, atom.loc))
            .collect::<Vec<_>>(),
        before
    );
}

#[test]
fn collision_open_angle_processes_initial_pair_snapshot_in_source_order() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let points = [
        (0, [0.0, 0.0]),
        (1, [1.0, 0.0]),
        (2, [0.1, 0.1]),
        (3, [0.2, 0.2]),
    ];
    let mut actual = collision_fragment(&graph, &rings, &points);
    let mut expected = collision_fragment(&graph, &rings, &points);
    for atom in actual.atoms.values_mut() {
        atom.fixed = false;
    }
    for atom in expected.atoms.values_mut() {
        atom.fixed = false;
    }
    let distance = expected.collision_distance_matrix().unwrap();
    let original_pairs = expected.find_collisions(&distance, false);
    assert_eq!(original_pairs, [(2, 0), (3, 0), (3, 2)]);
    for (first, second) in original_pairs {
        expected.open_angles(&distance, first, second).unwrap();
    }
    actual.remove_collisions_open_angles().unwrap();
    for (&id, row) in &expected.atoms {
        assert_eq!(actual.atoms[&id].loc, row.loc);
    }
}

#[test]
fn collision_shortening_contracts_only_unfixed_terminal_to_strict_threshold() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [2.0, 0.0]),
            (3, [0.1, 0.0]),
        ],
    );
    frag.atoms.get_mut(&0).unwrap().fixed = false;
    frag.remove_collisions_shorten_bonds().unwrap();
    assert!((frag.atoms[&0].loc[0] - 0.19).abs() < 1e-12);
    assert_eq!(frag.atoms[&3].loc, [0.1, 0.0]);
    assert_eq!(frag.atoms[&1].loc, [1.0, 0.0]);
}

#[test]
fn collision_shortening_skips_fixed_pair_and_disconnected_path() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[]);
    let rings = source_rings(&graph);
    let points = [(0, [0.0, 0.0]), (1, [0.1, 0.0])];
    let mut fixed = collision_fragment(&graph, &rings, &points);
    fixed.remove_collisions_shorten_bonds().unwrap();
    assert_eq!(fixed.atoms[&0].loc, points[0].1);
    assert_eq!(fixed.atoms[&1].loc, points[1].1);
    let mut disconnected = collision_fragment(&graph, &rings, &points);
    for atom in disconnected.atoms.values_mut() {
        atom.fixed = false;
    }
    disconnected.remove_collisions_shorten_bonds().unwrap();
    assert_eq!(disconnected.atoms[&0].loc, points[0].1);
    assert_eq!(disconnected.atoms[&1].loc, points[1].1);
}

#[test]
fn collision_shortening_ring_path_moves_inward_and_keeps_fixed_row() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3), (3, 0)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.1, 0.1]),
            (3, [0.0, 1.0]),
        ],
    );
    for atom in frag.atoms.values_mut() {
        atom.fixed = false;
    }
    frag.atoms.get_mut(&1).unwrap().fixed = true;
    let before = frag
        .atoms
        .iter()
        .map(|(&id, row)| (id, row.loc))
        .collect::<Vec<_>>();
    frag.remove_collisions_shorten_bonds().unwrap();
    assert_eq!(frag.atoms[&1].loc, before[1].1);
    assert!(
        frag.atoms
            .iter()
            .any(|(&id, row)| id != 1 && row.loc != before[id].1)
    );
}

#[test]
fn collision_repair_runs_flip_then_open_angle_then_shortening_in_source_order() {
    // RDDepictor.cpp::compute2DCoords runs the flip pass before the two
    // residual-collision passes. The first pass resolves this chain overlap.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (1, 2), (2, 3), (3, 4)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [2.0, 0.0]),
            (3, [2.0, 1.0]),
            (4, [0.1, 0.1]),
        ],
    );
    for row in frag.atoms.values_mut() {
        row.fixed = false;
    }
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, true), [(4, 0)]);
    frag.remove_collisions_bond_flip().unwrap();
    assert!(frag.find_collisions(&distance, true).is_empty());
    let after_flip = frag
        .atoms
        .iter()
        .map(|(&id, row)| (id, row.loc))
        .collect::<Vec<_>>();
    frag.remove_collisions_open_angles().unwrap();
    frag.remove_collisions_shorten_bonds().unwrap();
    assert_eq!(
        frag.atoms
            .iter()
            .map(|(&id, row)| (id, row.loc))
            .collect::<Vec<_>>(),
        after_flip
    );
}

#[test]
fn collision_repair_leaves_fixed_ring_overlap_unresolved_and_bounded() {
    // A complete embedded ring is a valid fragment. Fixed rows cannot be
    // moved by the residual passes, and ring edges cannot be flipped.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 3), (3, 0)],
    );
    let rings = source_rings(&graph);
    let mut frag = collision_fragment(
        &graph,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.1, 0.1]),
            (3, [0.0, 1.0]),
        ],
    );
    let distance = frag.collision_distance_matrix().unwrap();
    assert_eq!(frag.find_collisions(&distance, true), [(2, 0)]);
    frag.remove_collisions_bond_flip().unwrap();
    frag.remove_collisions_open_angles().unwrap();
    frag.remove_collisions_shorten_bonds().unwrap();
    assert_eq!(frag.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(frag.atoms[&2].loc, [0.1, 0.1]);
    assert_eq!(frag.find_collisions(&distance, true), [(2, 0)]);
}

#[test]
fn collision_repair_shortest_path_tie_cannot_flip_ring_or_stereo_bonds() {
    // The two equal routes through this cycle are both ring edges, so the
    // source shortest-path choice cannot make either interior edge flippable.
    let ring = topology(
        vec![AtomSpec::new(Element::C); 6],
        &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
    );
    let rings = source_rings(&ring);
    let mut frag = collision_fragment(
        &ring,
        &rings,
        &[
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [2.0, 0.0]),
            (3, [2.0, 1.0]),
            (4, [1.0, 1.0]),
            (5, [0.1, 0.1]),
        ],
    );
    for row in frag.atoms.values_mut() {
        row.fixed = false;
    }
    assert!(
        frag.rotatable_bonds_on_shortest_path(0, 3)
            .unwrap()
            .is_empty()
    );
    let before = frag
        .atoms
        .iter()
        .map(|(&id, row)| (id, row.loc))
        .collect::<Vec<_>>();
    frag.remove_collisions_bond_flip().unwrap();
    assert_eq!(
        frag.atoms
            .iter()
            .map(|(&id, row)| (id, row.loc))
            .collect::<Vec<_>>(),
        before
    );

    let stereo = cis_trans_graph(
        5,
        &[
            (0, 1, BondOrder::Single, BondStereo::None, None),
            (1, 2, BondOrder::Double, BondStereo::E, Some((0, 3))),
            (2, 3, BondOrder::Single, BondStereo::None, None),
            (3, 4, BondOrder::Single, BondStereo::None, None),
        ],
    );
    let rings = source_rings(&stereo);
    let frag = collision_fragment(&stereo, &rings, &[(0, [0.0, 0.0]), (4, [0.1, 0.1])]);
    assert_eq!(frag.rotatable_bonds_on_shortest_path(0, 4).unwrap(), [2]);
}

fn cis_trans_graph(
    count: usize,
    edges: &[(usize, usize, BondOrder, BondStereo, Option<(usize, usize)>)],
) -> TopologyBlock {
    let atoms = (0..count)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect();
    let bonds = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end, order, stereo, refs))| {
            let mut spec =
                BondSpec::new(AtomId::new(begin), AtomId::new(end), order).with_stereo(stereo);
            if let Some((left, right)) = refs {
                spec = spec.with_stereo_atoms(AtomId::new(left), AtomId::new(right));
            }
            Bond::from_spec(BondId::new(index), spec)
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
}

#[test]
fn cis_trans_e_z_cis_trans_seed_normals_and_substituent_ids() {
    // EmbeddedFrag.cpp:185-234 seeds only double-bond endpoints; Z/Cis puts
    // the second normal below the bond, E/Trans above it.
    for (stereo, end_normal, end_ccw) in [
        (BondStereo::Z, [0.0, -1.0], true),
        (BondStereo::Cis, [0.0, -1.0], true),
        (BondStereo::E, [0.0, 1.0], false),
        (BondStereo::Trans, [0.0, 1.0], false),
    ] {
        let graph = cis_trans_graph(
            4,
            &[
                (0, 1, BondOrder::Single, BondStereo::None, None),
                (1, 2, BondOrder::Double, stereo, Some((0, 3))),
                (2, 3, BondOrder::Single, BondStereo::None, None),
            ],
        );
        let rings = source_rings(&graph);
        let fragments = embedded_frag::embed_cis_trans_systems(&graph, &rings).unwrap();
        assert_eq!(fragments.len(), 1);
        let fragment = &fragments[0];
        assert_eq!(fragment.atoms.keys().copied().collect::<Vec<_>>(), [1, 2]);
        assert_eq!(fragment.atoms[&1].loc, [0.0, 0.0]);
        assert_eq!(fragment.atoms[&2].loc, [geometry::BOND_LEN, 0.0]);
        assert_eq!(fragment.atoms[&1].normal, [0.0, -1.0]);
        assert_eq!(fragment.atoms[&2].normal, end_normal);
        assert!(!fragment.atoms[&1].ccw);
        assert_eq!(fragment.atoms[&2].ccw, end_ccw);
        assert_eq!(fragment.atoms[&1].cis_trans_nbr, Some(0));
        assert_eq!(fragment.atoms[&2].cis_trans_nbr, Some(3));
        assert_eq!(fragment.atoms[&1].neighs, [0]);
        assert_eq!(fragment.atoms[&2].neighs, [3]);
        assert_eq!(fragment.attachment_points.len(), 2);
    }
}

#[test]
fn cis_trans_conjugated_seed_order_and_no_common_merge() {
    // RDDepictor.cpp::embedCisTransSystems scans bond-table order. The two
    // endpoint seeds are disjoint; expandEfrag joins them via mergeNoCommon
    // when it reaches the intervening single bond.
    let graph = cis_trans_graph(
        6,
        &[
            (0, 1, BondOrder::Single, BondStereo::None, None),
            (1, 2, BondOrder::Double, BondStereo::E, Some((0, 3))),
            (2, 3, BondOrder::Single, BondStereo::None, None),
            (3, 4, BondOrder::Double, BondStereo::Z, Some((2, 5))),
            (4, 5, BondOrder::Single, BondStereo::None, None),
        ],
    );
    let rings = source_rings(&graph);
    let mut fragments = embedded_frag::embed_cis_trans_systems(&graph, &rings).unwrap();
    assert_eq!(fragments.len(), 2);
    assert_eq!(
        fragments[0].atoms.keys().copied().collect::<Vec<_>>(),
        [1, 2]
    );
    assert_eq!(
        fragments[1].atoms.keys().copied().collect::<Vec<_>>(),
        [3, 4]
    );
    let mut first = fragments.remove(0);
    first
        .expand_fragment(&mut vec![0, 5], &mut fragments)
        .unwrap();
    assert!(fragments.is_empty());
    assert!(first.atoms.contains_key(&3));
    assert!(first.atoms.contains_key(&4));
    assert_eq!(first.atoms[&3].cis_trans_nbr, Some(2));
    assert_eq!(first.atoms[&4].cis_trans_nbr, Some(5));
}

#[test]
fn cis_trans_ring_bond_and_unspecified_double_bond_are_not_seeded() {
    // RDDepictor.cpp:298-318 excludes ring double bonds and stereo <= ANY.
    let ring = cis_trans_graph(
        3,
        &[
            (0, 1, BondOrder::Double, BondStereo::E, Some((2, 2))),
            (1, 2, BondOrder::Single, BondStereo::None, None),
            (2, 0, BondOrder::Single, BondStereo::None, None),
        ],
    );
    assert!(
        embedded_frag::embed_cis_trans_systems(&ring, &source_rings(&ring))
            .unwrap()
            .is_empty()
    );
    let plain = cis_trans_graph(2, &[(0, 1, BondOrder::Double, BondStereo::Any, None)]);
    assert!(
        embedded_frag::embed_cis_trans_systems(&plain, &source_rings(&plain))
            .unwrap()
            .is_empty()
    );
}

fn nontetra_graph(tag: ChiralTag, permutation: Option<u32>, ligands: &[Element]) -> TopologyBlock {
    let centre = match tag {
        ChiralTag::SquarePlanar => Element::PT,
        ChiralTag::TrigonalBipyramidal => Element::AS,
        ChiralTag::Octahedral => Element::CO,
        _ => Element::C,
    };
    let mut first = AtomSpec::new(centre).with_chiral_tag(tag);
    if let Some(value) = permutation {
        first = first.with_chiral_permutation(value);
    }
    let specs = std::iter::once(first)
        .chain(ligands.iter().copied().map(AtomSpec::new))
        .collect();
    let edges: Vec<_> = (1..=ligands.len()).map(|ligand| (0, ligand)).collect();
    topology(specs, &edges)
}

fn nontetra_points(graph: &TopologyBlock) -> geometry::PointMap {
    let rings = source_rings(graph);
    let fragments = nontetrahedral::embed_nontetrahedral_stereo(graph, &rings).unwrap();
    assert_eq!(fragments.len(), 1);
    assert!(fragments[0].atoms.values().all(|atom| atom.fixed));
    fragments[0]
        .atoms
        .iter()
        .map(|(&index, atom)| (index, atom.loc))
        .collect()
}

#[test]
fn nontetra_stereo_square_planar_sp1_source_quadrants_and_across() {
    // RDKit 2026.03.1: [Pt@SP1](F)(Cl)(Br)I, Compute2DCoords(canonOrient=False).
    let graph = nontetra_graph(
        ChiralTag::SquarePlanar,
        Some(1),
        &[Element::F, Element::CL, Element::BR, Element::I],
    );
    let d = 0.707107 * geometry::BOND_LEN;
    assert_eq!(
        nontetra_points(&graph),
        geometry::PointMap::from([
            (0, [0.0, 0.0]),
            (1, [d, d]),
            (2, [d, -d]),
            (3, [-d, -d]),
            (4, [-d, d]),
        ])
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_across_ligand(&graph, 0, 1),
        Some(3)
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_ideal_angle(&graph, 0, 1, 3),
        180.0
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_ideal_angle(&graph, 0, 1, 2),
        90.0
    );
}

#[test]
fn nontetra_stereo_trigonal_bipyramidal_tb1_axial_and_ranked_equatorial() {
    // RDKit 2026.03.1: [As@TB1](F)(Cl)(Br)(I)S; bond-order axial F/S.
    let graph = nontetra_graph(
        ChiralTag::TrigonalBipyramidal,
        Some(1),
        &[Element::F, Element::CL, Element::BR, Element::I, Element::S],
    );
    let e = 0.866025 * geometry::BOND_LEN;
    assert_eq!(
        nontetra_points(&graph),
        geometry::PointMap::from([
            (0, [0.0, 0.0]),
            (1, [0.0, 1.5]),
            (2, [-e, 0.75]),
            (3, [-e, -0.75]),
            (4, [1.5, 0.0]),
            (5, [0.0, -1.5]),
        ])
    );
    assert_eq!(
        cosmolkit_core::trigonal_bipyramidal_axial_ligand(&graph, 0, 1),
        Some(1)
    );
    assert_eq!(
        cosmolkit_core::trigonal_bipyramidal_axial_ligand(&graph, 0, -1),
        Some(5)
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_ideal_angle(&graph, 0, 1, 5),
        180.0
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_ideal_angle(&graph, 0, 1, 2),
        90.0
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_ideal_angle(&graph, 0, 2, 3),
        120.0
    );
}

#[test]
fn nontetra_stereo_octahedral_oh1_axial_across_and_source_ligand_order() {
    // RDKit 2026.03.1: [Co@OH1](F)(Cl)(Br)(I)(S)N; final depiction may rotate.
    let graph = nontetra_graph(
        ChiralTag::Octahedral,
        Some(1),
        &[
            Element::F,
            Element::CL,
            Element::BR,
            Element::I,
            Element::S,
            Element::N,
        ],
    );
    let e = 0.866025 * geometry::BOND_LEN;
    assert_eq!(
        nontetra_points(&graph),
        geometry::PointMap::from([
            (0, [0.0, 0.0]),
            (1, [0.0, -1.5]),
            (2, [e, -0.75]),
            (3, [-e, -0.75]),
            (4, [-e, 0.75]),
            (5, [e, 0.75]),
            (6, [0.0, 1.5]),
        ])
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_across_ligand(&graph, 0, 1),
        Some(6)
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_across_ligand(&graph, 0, 2),
        Some(4)
    );
}

#[test]
fn nontetra_stereo_missing_permutation_keeps_source_partial_state() {
    let graph = nontetra_graph(
        ChiralTag::SquarePlanar,
        None,
        &[Element::F, Element::CL, Element::BR, Element::I],
    );
    assert_eq!(
        cosmolkit_core::non_tetrahedral_across_ligand(&graph, 0, 1),
        None
    );
    let d = 0.707107 * geometry::BOND_LEN;
    let points = nontetra_points(&graph);
    assert_eq!(points[&1], [d, d]);
    assert_eq!(points[&2], [d, -d]);
    assert_eq!(points[&3], [-d, d]);
    assert_eq!(points[&4], [-d, d]);
    let untagged = topology(vec![AtomSpec::new(Element::C)], &[]);
    assert!(
        nontetrahedral::embed_nontetrahedral_stereo(&untagged, &source_rings(&untagged))
            .unwrap()
            .is_empty()
    );
}

fn topology(specs: Vec<AtomSpec>, edges: &[(usize, usize)]) -> TopologyBlock {
    let atoms = specs
        .into_iter()
        .enumerate()
        .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
        .collect();
    let bonds = edges
        .iter()
        .copied()
        .enumerate()
        .map(|(index, (begin, end))| {
            Bond::from_spec(
                BondId::new(index),
                BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
            )
        })
        .collect();
    TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).expect("valid test topology")
}

fn template_match_registry(rows: &str) -> templates::CoordinateTemplates {
    static NEXT: AtomicUsize = AtomicUsize::new(0);
    let path = std::env::temp_dir().join(format!(
        "cosmolkit-template-match-{}-{}.cxsmarts",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::write(&path, rows).expect("write owned template fixture");
    let mut registry = templates::CoordinateTemplates::default();
    let result = registry.set_ring_system_templates(&path);
    std::fs::remove_file(&path).expect("remove owned template fixture");
    result.expect("source-valid template fixture");
    registry
}

fn template_match_fragment<'a>(
    graph: &'a TopologyBlock,
    rings: &'a RingInfo,
) -> embedded_frag::EmbeddedFrag<'a> {
    embedded_frag::EmbeddedFrag::from_coord_map(graph, rings, &geometry::PointMap::new())
        .expect("empty source fragment")
}

fn source_rings(graph: &TopologyBlock) -> RingInfo {
    symmetrize_sssr_with_options_from_parts(
        graph.atoms.len(),
        &graph.bonds,
        &graph.adjacency,
        false,
        false,
    )
    .expect("source symmetrized rings")
}

#[test]
fn ring_embed_atom_overlap_groups_spiro_but_not_disconnected_rings() {
    // RingUtils::makeRingNeighborMap intersects atom rings, not bond rings.
    assert_eq!(
        embedded_frag::ring_systems(&[vec![0, 1, 2], vec![2, 3, 4], vec![5, 6, 7], vec![4, 8, 9],]),
        vec![vec![0, 1, 3], vec![2]]
    );
}

#[test]
fn ring_embed_single_ring_uses_source_polygon_and_cyclic_neighbors() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 3],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let rings = source_rings(&graph);
    let mut templates = template_match_registry("[*]1-[*]-[*]-1 |(9,9,;10,9,;9,10,)|\n");
    let fragments =
        embedded_frag::embed_fused_systems(&graph, &rings, None, true, &mut templates).unwrap();
    assert_eq!(fragments.len(), 1);
    assert_eq!(fragments[0].atoms.len(), 3);
    let ring: Vec<_> = rings.atom_rings()[0]
        .iter()
        .map(|atom| atom.index())
        .collect();
    let expected = geometry::embed_ring(&ring);
    for (index, &atom) in ring.iter().enumerate() {
        let state = &fragments[0].atoms[&atom];
        assert_eq!(state.loc, expected[&atom]);
        assert_eq!(
            state.nbr1,
            Some(ring[(index + ring.len() - 1) % ring.len()])
        );
        assert_eq!(state.nbr2, Some(ring[(index + 1) % ring.len()]));
        assert!(
            !state.fixed,
            "source ignores templates for rings of size <= 8"
        );
    }
}

#[test]
fn ring_embed_fused_spiro_and_bridged_complete_source_ring_systems() {
    let cases: [(&str, usize, &[(usize, usize)]); 3] = [
        ("fused", 4, &[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)]),
        (
            "spiro",
            5,
            &[(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)],
        ),
        (
            "bridged",
            7,
            &[
                (0, 1),
                (1, 2),
                (2, 3),
                (0, 4),
                (4, 5),
                (5, 3),
                (0, 6),
                (6, 3),
            ],
        ),
    ];
    for (name, count, edges) in cases {
        let graph = topology(vec![AtomSpec::new(Element::C); count], edges);
        let rings = source_rings(&graph);
        let actual_rings: Vec<Vec<usize>> = rings
            .atom_rings()
            .iter()
            .map(|ring| ring.iter().map(|atom| atom.index()).collect())
            .collect();
        let expected_rings: &[&[usize]] = match name {
            "fused" => &[&[0, 1, 2], &[3, 1, 2]],
            "spiro" => &[&[1, 0, 2], &[3, 0, 4]],
            "bridged" => &[&[1, 0, 6, 3, 2], &[4, 0, 6, 3, 5]],
            _ => unreachable!(),
        };
        assert_eq!(actual_rings, expected_rings, "{name} ring traversal");
        let mut templates = templates::CoordinateTemplates::default();
        let fragments =
            embedded_frag::embed_fused_systems(&graph, &rings, None, false, &mut templates)
                .unwrap_or_else(|error| panic!("{name}: {error:?}"));
        assert_eq!(fragments.len(), 1, "{name}");
        assert_eq!(fragments[0].atoms.len(), count, "{name}");
        assert!(
            fragments[0]
                .atoms
                .values()
                .all(|atom| atom.loc.iter().all(|coordinate| coordinate.is_finite())),
            "{name}"
        );
        let reference: &[[f64; 2]] = match name {
            "fused" => &[
                [0.8660254, 0.0],
                [-0.4330127, 0.75],
                [-0.4330127, -0.75],
                [-1.73205081, 0.0],
            ],
            "spiro" => &[
                [-0.4330127, 0.75],
                [0.8660254, 0.0],
                [-0.4330127, -0.75],
                [-1.73205081, 1.5],
                [-0.4330127, 2.25],
            ],
            "bridged" => &[
                [0.39429833, 1.21352549],
                [1.27597621, 0.0],
                [0.39429833, -1.21352549],
                [-1.03228644, -0.75],
                [-1.03228644, 1.67705098],
                [-1.91396432, 0.46352549],
                [-1.03228644, 0.75],
            ],
            _ => unreachable!(),
        };
        // Final RDKit translation/orientation may differ from this private
        // fragment, but source-preserving rigid transforms retain distances.
        for left in 0..count {
            for right in left + 1..count {
                let actual_left = fragments[0].atoms[&left].loc;
                let actual_right = fragments[0].atoms[&right].loc;
                let actual =
                    (actual_left[0] - actual_right[0]).hypot(actual_left[1] - actual_right[1]);
                let expected = (reference[left][0] - reference[right][0])
                    .hypot(reference[left][1] - reference[right][1]);
                assert!(
                    (actual - expected).abs() < 1e-6,
                    "{name} pair ({left},{right}) actual={actual} expected={expected} coordinates={:?}",
                    fragments[0].atoms
                );
            }
        }
    }
}

#[test]
fn ring_embed_first_matching_template_places_exact_xy_unless_two_atoms_constrained() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)],
    );
    let rings = source_rings(&graph);
    assert_eq!(rings.num_rings(), 2);
    let rows = "[*]1-[*]2-[*]-1-[*]-2 |(10,0,;11,0,;11,1,;10,1,)|\n\
[*]1-[*]2-[*]-1-[*]-2 |(20,0,;21,0,;21,1,;20,1,)|\n";
    let mut templates = template_match_registry(rows);
    let fragments =
        embedded_frag::embed_fused_systems(&graph, &rings, None, true, &mut templates).unwrap();
    assert_eq!(fragments.len(), 1);
    assert!(fragments[0].atoms.values().all(|atom| atom.fixed));
    assert_eq!(fragments[0].atoms[&0].loc, [10.0, 0.0]);
    assert_eq!(fragments[0].atoms[&1].loc, [11.0, 0.0]);
    assert_eq!(fragments[0].atoms[&2].loc, [11.0, 1.0]);
    assert_eq!(fragments[0].atoms[&3].loc, [10.0, 1.0]);

    let one = geometry::PointMap::from([(0, [99.0, 99.0])]);
    let allowed =
        embedded_frag::embed_fused_systems(&graph, &rings, Some(&one), true, &mut templates)
            .unwrap();
    assert!(allowed[0].atoms.values().all(|atom| atom.fixed));
    let two = geometry::PointMap::from([(0, [99.0, 99.0]), (1, [98.0, 99.0])]);
    let disabled =
        embedded_frag::embed_fused_systems(&graph, &rings, Some(&two), true, &mut templates)
            .unwrap();
    assert!(disabled[0].atoms.values().all(|atom| !atom.fixed));
    assert_ne!(disabled[0].atoms[&0].loc, [10.0, 0.0]);
}

#[test]
fn template_match_prefilters_size_bonds_rings_and_degree_before_query() {
    // EmbeddedFrag.cpp::matchToTemplate checks size bucket, ring-only bond
    // count, ring count and degree histogram in that order.
    let triangle = topology(
        vec![AtomSpec::new(Element::C); 3],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 3, 3);
    let mut fragment = template_match_fragment(&triangle, &rings);
    let mut square_templates =
        template_match_registry("[*]1-[*]-[*]-[*]-1 |(0,0,;1,0,;1,1,;0,1,)|\n");
    assert!(
        !fragment
            .match_to_template(&[0, 1, 2], 1, &mut square_templates)
            .unwrap()
    );
    assert!(fragment.atoms.is_empty());

    let mut triangle_templates = template_match_registry("[*]1-[*]-[*]-1 |(0,0,;1,0,;0,1,)|\n");
    assert!(
        !fragment
            .match_to_template(&[0, 1], 1, &mut triangle_templates)
            .unwrap()
    );
    let path = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let path_rings = RingInfo::new(RingFindType::Fast, 3, 2);
    let mut path_fragment = template_match_fragment(&path, &path_rings);
    assert!(
        !path_fragment
            .match_to_template(&[0, 1, 2], 1, &mut triangle_templates)
            .unwrap()
    );
    assert!(path_fragment.atoms.is_empty());
    assert!(
        !fragment
            .match_to_template(&[0, 1, 2], 2, &mut triangle_templates)
            .unwrap()
    );
    assert!(fragment.atoms.is_empty());

    let pendant = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 0), (0, 3)],
    );
    let pendant_rings = RingInfo::new(RingFindType::Fast, 4, 4);
    let mut pendant_fragment = template_match_fragment(&pendant, &pendant_rings);
    assert!(
        !pendant_fragment
            .match_to_template(&[0, 1, 2, 3], 1, &mut square_templates)
            .unwrap()
    );
    assert!(pendant_fragment.atoms.is_empty());
}

#[test]
fn template_match_sentinel_excludes_external_ring_and_maps_exact_xy() {
    // Source sets non-ring-system atomic numbers to 200; [!#200] therefore
    // cannot map the first disconnected triangle when the second is selected.
    let graph = topology(
        vec![AtomSpec::new(Element::C); 6],
        &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 6, 6);
    let mut fragment = template_match_fragment(&graph, &rings);
    let mut registry = template_match_registry("[!#200]1-[!#200]-[!#200]-1 |(0,0,;1,0,;0,1,)|\n");
    assert!(
        fragment
            .match_to_template(&[3, 4, 5], 1, &mut registry)
            .unwrap()
    );
    assert_eq!(
        fragment.atoms.keys().copied().collect::<Vec<_>>(),
        vec![3, 4, 5]
    );
    assert!(fragment.atoms.values().all(|atom| atom.fixed));
    assert_eq!(fragment.atoms[&3].loc, [0.0, 0.0]);
    assert_eq!(fragment.atoms[&4].loc, [1.0, 0.0]);
    assert_eq!(fragment.atoms[&5].loc, [0.0, 1.0]);
    assert!(fragment.attachment_points.is_empty());
}

#[test]
fn template_match_preserves_source_order_and_query_bond_predicates() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 3],
        &[(0, 1), (1, 2), (2, 0)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 3, 3);
    let mut fragment = template_match_fragment(&graph, &rings);
    let mut registry = template_match_registry(concat!(
        "[*]1=[*]-[*]-1 |(7,7,;8,7,;7,8,)|\n",
        "[!#200]1-[!#200]-[!#200]-1 |(0,0,;1,0,;0,1,)|\n",
        "[*]1-[*]-[*]-1 |(9,9,;10,9,;9,10,)|\n",
    ));
    assert!(
        fragment
            .match_to_template(&[0, 1, 2], 1, &mut registry)
            .unwrap()
    );
    assert_eq!(fragment.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(fragment.atoms[&1].loc, [1.0, 0.0]);
    assert_eq!(fragment.atoms[&2].loc, [0.0, 1.0]);
    let mut rejected = template_match_fragment(&graph, &rings);
    let mut only_double = template_match_registry("[*]1=[*]-[*]-1 |(7,7,;8,7,;7,8,)|\n");
    assert!(
        !rejected
            .match_to_template(&[0, 1, 2], 1, &mut only_double)
            .unwrap()
    );
    assert!(rejected.atoms.is_empty());
}

#[test]
fn template_match_rejects_first_stereo_geometry_then_accepts_next() {
    // EmbeddedFrag.cpp::checkStereoChemistry rejects the first cis-shaped
    // template for a source E bond, then matchToTemplate continues in row order.
    let atoms = (0..4)
        .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
        .collect();
    let edges = [(0, 1), (1, 2), (2, 3), (3, 0)];
    let bonds = edges
        .into_iter()
        .enumerate()
        .map(|(index, (begin, end))| {
            let mut spec = BondSpec::new(
                AtomId::new(begin),
                AtomId::new(end),
                if index == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                },
            );
            if index == 0 {
                spec = spec
                    .with_stereo(cosmolkit_model::BondStereo::E)
                    .with_stereo_atoms(AtomId::new(3), AtomId::new(2));
            }
            Bond::from_spec(BondId::new(index), spec)
        })
        .collect();
    let graph = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
    let rings = RingInfo::new(RingFindType::Fast, 4, 4);
    let mut fragment = template_match_fragment(&graph, &rings);
    let mut registry = template_match_registry(concat!(
        "[!#200]1=[!#200]-[!#200]-[!#200]-1 |(0,0,;1,0,;1,1,;0,1,)|\n",
        "[!#200]1=[!#200]-[!#200]-[!#200]-1 |(0,0,;1,0,;1,-1,;0,1,)|\n",
    ));
    assert!(
        fragment
            .match_to_template(&[0, 1, 2, 3], 1, &mut registry)
            .unwrap()
    );
    assert_eq!(fragment.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(fragment.atoms[&1].loc, [1.0, 0.0]);
    assert_eq!(fragment.atoms[&2].loc, [1.0, -1.0]);
    assert_eq!(fragment.atoms[&3].loc, [0.0, 1.0]);
    assert!(fragment.atoms.values().all(|atom| atom.fixed));
}

#[test]
fn template_match_refreshes_neighbor_and_attachment_state() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (2, 0), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 4, 4);
    let mut fragment = template_match_fragment(&graph, &rings);
    let mut registry = template_match_registry("[!#200]1-[!#200]-[!#200]-1 |(0,0,;1,0,;0,1,)|\n");
    assert!(
        fragment
            .match_to_template(&[0, 1, 2], 1, &mut registry)
            .unwrap()
    );
    assert_eq!(fragment.attachment_points, vec![0]);
    assert_eq!(fragment.atoms[&0].neighs, vec![3]);
    assert!(!fragment.atoms.contains_key(&3));
    assert!(fragment.atoms.values().all(|atom| atom.fixed));
}

#[test]
fn geometry_ranks_ring_and_affine_point_sources() {
    let ring = geometry::embed_ring(&[0, 1, 2, 3, 4, 5]);
    assert_eq!(ring.len(), 6);
    assert!((ring[&0][0] - geometry::BOND_LEN).abs() < 1e-12);
    assert!(ring[&0][1].abs() < 1e-12);
    let a = ring[&0];
    let b = ring[&1];
    assert!(((a[0] - b[0]).hypot(a[1] - b[1]) - geometry::BOND_LEN).abs() < 1e-12);

    let transform =
        geometry::Transform2D::from_point_pairs([0.0, 0.0], [1.0, 0.0], [2.0, 3.0], [2.0, 4.0]);
    let transformed = transform.transform_point([2.0, 4.0]);
    assert!((transformed[0] - 1.0).abs() < 1e-12);
    assert!(transformed[1].abs() < 1e-12);
    let mut points = geometry::PointMap::from([(0, [2.0, 3.0]), (1, [2.0, 4.0])]);
    geometry::transform_points(&mut points, transform);
    assert!(points[&0][0].abs() < 1e-12 && points[&0][1].abs() < 1e-12);
    assert!((points[&1][0] - 1.0).abs() < 1e-12 && points[&1][1].abs() < 1e-12);
    let degenerate =
        geometry::Transform2D::from_point_pairs([1.0, 1.0], [1.0, 1.0], [3.0, 4.0], [5.0, 6.0]);
    assert_eq!(degenerate.transform_point([7.0, 8.0]), [7.0, 8.0]);
}

#[test]
fn geometry_ranks_bisect_reflect_and_discarded_map_result() {
    assert_eq!(
        geometry::compute_bisect_point([0.0, 0.0], 0.0, [2.0, 0.0], [0.0, 2.0]),
        [1.0, 1.0]
    );
    assert_eq!(
        geometry::compute_bisect_point([0.0, 0.0], PI, [2.0, 0.0], [0.0, 2.0]),
        [1.0, 1.0]
    );
    assert_eq!(
        geometry::compute_bisect_point([0.0, 0.0], PI + 0.01, [2.0, 0.0], [0.0, 2.0]),
        [-1.0, -1.0]
    );
    let reflected = geometry::reflect_point([1.0, 2.0], [0.0, 0.0], [2.0, 0.0]);
    assert!((reflected[0] - 1.0).abs() < 1e-12);
    assert!((reflected[1] + 2.0).abs() < 1e-12);
    let original = geometry::PointMap::from([(0, [1.0, 2.0])]);
    let mut points = original.clone();
    geometry::reflect_points(&mut points, [0.0, 0.0], [2.0, 0.0]);
    // Pinned DepictUtils.cpp discards reflectPoint's return value here.
    assert_eq!(points, original);
}

#[test]
fn geometry_ranks_source_rank_property_precedence_and_hydrogen_order() {
    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::N)
                .with_prop("_CIPRank", "0")
                .unwrap(),
            AtomSpec::new(Element::O)
                .with_prop("_ChiralAtomRank", "1")
                .unwrap(),
            AtomSpec::new(Element::H),
        ],
        &[],
    );
    assert_eq!(geometry::atom_depict_rank(&graph, 3).unwrap(), 100_000);
    assert_eq!(
        geometry::rank_atoms_by_rank(&graph, &[0, 1, 2, 3], true).unwrap(),
        vec![1, 0, 2, 3]
    );
    assert_eq!(
        geometry::rank_atoms_by_rank(&graph, &[0, 1, 2, 3], false).unwrap(),
        vec![3, 2, 0, 1]
    );
    assert_eq!(
        geometry::rank_atoms_by_rank(&graph, &[4], true),
        Err(geometry::GeometryError::AtomIndexOutOfRange {
            atom: 4,
            atom_count: 4,
        })
    );
    let tied = topology(
        vec![
            AtomSpec::new(Element::C)
                .with_prop("_CIPRank", "1")
                .unwrap(),
            AtomSpec::new(Element::N)
                .with_prop("_CIPRank", "1")
                .unwrap(),
        ],
        &[],
    );
    assert_eq!(
        geometry::rank_atoms_by_rank(&tied, &[1, 0], true).unwrap(),
        vec![0, 1]
    );
}

#[test]
fn geometry_ranks_degree_four_reference_rotation_matches_source() {
    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::N),
            AtomSpec::new(Element::O),
            AtomSpec::new(Element::F),
        ],
        &[(0, 1), (0, 2), (0, 3), (0, 4)],
    );
    assert_eq!(
        geometry::set_neighbor_order(&graph, 0, &[1, 2, 3]).unwrap(),
        vec![1, 3, 2]
    );
    assert_eq!(
        geometry::set_neighbor_order(&graph, 0, &[1, 2, 3, 4]).unwrap(),
        vec![1, 3, 2, 4]
    );
    assert_eq!(
        geometry::set_neighbor_order(&graph, 0, &[1, 2]),
        Err(geometry::GeometryError::NotEnoughNeighbors { atom: 0, count: 3 })
    );
}

#[test]
fn fragment_seed_neighbors_singleton_and_empty_map_match_source() {
    let graph = topology(vec![AtomSpec::new(Element::C)], &[]);
    let rings = RingInfo::new(RingFindType::Fast, 1, 0);
    let seeded = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    assert_eq!(seeded.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(seeded.atoms[&0].normal, [1.0, 0.0]);
    assert_eq!(seeded.atoms[&0].angle, -1.0);
    assert!(!seeded.atoms[&0].fixed);
    assert!(seeded.attachment_points.is_empty());
    assert_eq!(seeded.find_neighbor(0).unwrap(), None);
    assert!(!seeded.done);
    let empty =
        embedded_frag::EmbeddedFrag::from_coord_map(&graph, &rings, &geometry::PointMap::new())
            .unwrap();
    assert!(empty.atoms.is_empty());
    assert!(empty.attachment_points.is_empty());
    assert!(!empty.done);
    assert_eq!(
        embedded_frag::EmbeddedFrag::from_single(1, &graph, &rings).unwrap_err(),
        embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 1,
            atom_count: 1
        }
    );
    assert_eq!(
        empty.find_neighbor(1),
        Err(embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 1,
            atom_count: 1
        })
    );
}

#[test]
fn fragment_seed_neighbors_rank_unembedded_rows_and_keep_map_points_fixed() {
    let graph = topology(
        vec![
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::H),
            AtomSpec::new(Element::O),
            AtomSpec::new(Element::N)
                .with_prop("_CIPRank", "0")
                .unwrap(),
        ],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let seeded = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    assert_eq!(seeded.atoms[&0].neighs, vec![3, 2, 1]);
    assert_eq!(seeded.attachment_points, vec![0]);

    let mapped = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [2.0, 3.0]), (2, [3.0, 3.0])]),
    )
    .unwrap();
    assert_eq!(mapped.atoms[&0].loc, [2.0, 3.0]);
    assert_eq!(mapped.atoms[&2].loc, [3.0, 3.0]);
    assert!(mapped.atoms.values().all(|atom| atom.fixed));
    assert_eq!(mapped.atoms[&0].neighs, vec![3, 1]);
    assert_eq!(mapped.atoms[&0].nbr1, Some(2));
    assert_eq!(mapped.atoms[&0].normal, [-0.0, 1.0]);
    assert_eq!(mapped.find_neighbor(1).unwrap(), Some(0));
    assert_eq!(mapped.attachment_points, vec![0]);
}

#[test]
fn fragment_seed_neighbors_two_and_three_done_neighbors_follow_source_angle_order() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 5],
        &[(0, 1), (0, 2), (0, 3), (0, 4)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let two = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    assert_eq!(two.atoms[&0].nbr1, Some(1));
    assert_eq!(two.atoms[&0].nbr2, Some(2));
    assert!((two.atoms[&0].angle - PI / 2.0).abs() < 1e-12);

    let three = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.0, 1.0]),
            (3, [-1.0, 0.0]),
        ]),
    )
    .unwrap();
    // The pinned `nbi3++` loop includes self-pairs, and stable angle sorting
    // chooses that zero-angle pair before the winning 180-degree pair.
    assert_eq!(three.atoms[&0].nbr1, Some(1));
    assert_eq!(three.atoms[&0].nbr2, Some(1));
    assert!((three.atoms[&0].angle - PI).abs() < 1e-12);
    assert_eq!(three.atoms[&0].rot_dir, -1);
    assert_eq!(three.atoms[&0].neighs, vec![4]);
}

#[test]
fn fragment_seed_neighbors_reject_invalid_map_and_coincident_embedded_points() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (0, 2)]);
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    assert_eq!(
        embedded_frag::EmbeddedFrag::from_coord_map(
            &graph,
            &rings,
            &geometry::PointMap::from([(3, [0.0, 0.0])]),
        )
        .unwrap_err(),
        embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 3,
            atom_count: 3
        }
    );
    assert_eq!(
        embedded_frag::EmbeddedFrag::from_coord_map(
            &graph,
            &rings,
            &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [0.0, 0.0])]),
        )
        .unwrap_err(),
        embedded_frag::FragmentError::CoincidentPoints
    );
}

#[test]
fn fragment_nonring_add_first_chain_and_branch_follow_source_angle_state() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (1, 2), (1, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    fragment.add_non_ring_atom(1, 0).unwrap();
    let first = fragment.atoms[&1].loc;
    assert!((first[0] - 0.75_f64.sqrt() * 1.5).abs() < 1e-12);
    assert!((first[1] - 0.75).abs() < 1e-12);
    assert_eq!(fragment.atoms[&0].nbr1, Some(1));
    assert!(fragment.atoms[&0].neighs.is_empty());
    fragment.add_non_ring_atom(2, 1).unwrap();
    let second = fragment.atoms[&2].loc;
    assert!((second[0] - 1.5 * 3.0_f64.sqrt()).abs() < 1e-12);
    assert!(second[1].abs() < 1e-12);
    assert_eq!(fragment.atoms[&1].nbr2, Some(2));
    assert!((fragment.atoms[&1].angle - 2.0 * PI / 3.0).abs() < 1e-12);
    fragment.add_non_ring_atom(3, 1).unwrap();
    let third = fragment.atoms[&3].loc;
    assert!((third[0] - first[0]).abs() < 1e-12);
    assert!((third[1] - 2.25).abs() < 1e-12);
    assert_eq!(fragment.atoms[&1].nbr2, Some(3));
    assert!(fragment.atoms[&1].neighs.is_empty());
}

#[test]
fn fragment_nonring_add_positive_angle_uses_source_rotation_and_density_tie() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut fragment = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    assert_eq!(fragment.atoms[&0].angle, PI / 2.0);
    assert_eq!(fragment.find_num_neigh([0.0, 0.0], 0.5), 1);
    assert_eq!(fragment.find_num_neigh([0.5, 0.0], 0.5), 0);
    fragment.add_non_ring_atom(3, 0).unwrap();
    let point = fragment.atoms[&3].loc;
    assert!((point[0] + 2.0_f64.sqrt() / 2.0).abs() < 1e-12);
    assert!((point[1] + 2.0_f64.sqrt() / 2.0).abs() < 1e-12);
    assert!((fragment.atoms[&0].angle - 5.0 * PI / 4.0).abs() < 1e-12);
    assert_eq!(fragment.atoms[&0].nbr2, Some(3));
    assert!(fragment.atoms[&0].neighs.is_empty());
}

#[test]
fn fragment_nonring_add_cis_trans_neighbor_reverses_only_nonmatching_direction() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (0, 2)]);
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut matching = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    matching.atoms.get_mut(&0).unwrap().cis_trans_nbr = Some(1);
    matching.add_non_ring_atom(1, 0).unwrap();
    let mut opposite = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    opposite.atoms.get_mut(&0).unwrap().cis_trans_nbr = Some(2);
    opposite.add_non_ring_atom(1, 0).unwrap();
    let a = matching.atoms[&1].loc;
    let b = opposite.atoms[&1].loc;
    assert!((a[0] + b[0]).abs() < 1e-12);
    assert!((a[1] - b[1]).abs() < 1e-12);
    assert_ne!(matching.atoms[&1].ccw, opposite.atoms[&1].ccw);
}

#[test]
fn fragment_nonring_add_rejects_preconditions_without_mutation() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let initial = fragment.atoms.clone();
    assert_eq!(
        fragment.add_non_ring_atom(2, 0),
        Err(embedded_frag::FragmentError::AtomIndexOutOfRange {
            atom: 2,
            atom_count: 2
        })
    );
    assert_eq!(
        fragment.add_non_ring_atom(1, 1),
        Err(embedded_frag::FragmentError::AtomNotEmbedded { atom: 1 })
    );
    assert_eq!(
        fragment.add_non_ring_atom(0, 0),
        Err(embedded_frag::FragmentError::AtomAlreadyEmbedded { atom: 0 })
    );
    assert_eq!(fragment.atoms, initial);
}

#[test]
fn fragment_merge_zero_common_bridge_preserves_both_fragment_rows() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut left = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut right = embedded_frag::EmbeddedFrag::from_single(1, &graph, &rings).unwrap();
    assert!(left.find_common_atoms(&right).is_empty());
    left.merge_no_common(&mut right, 0, 1).unwrap();
    assert_eq!(left.atoms.len(), 2);
    assert_eq!(left.find_common_atoms(&right), vec![0, 1]);
    let a = left.atoms[&0].loc;
    let b = left.atoms[&1].loc;
    assert!(((a[0] - b[0]).hypot(a[1] - b[1]) - geometry::BOND_LEN).abs() < 1e-12);
    assert_eq!(left.atoms[&1].aid, 0); // pinned operator= omits aid
}

#[test]
fn fragment_merge_one_common_extends_source_common_order_and_copies_state() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut left = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    let mut right = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (3, [-1.0, 0.0])]),
    )
    .unwrap();
    let mut common = left.find_common_atoms(&right);
    assert_eq!(common, vec![0]);
    left.merge_with_common(&mut right, &mut common).unwrap();
    assert_eq!(common, vec![0, 1]);
    assert_eq!(left.atoms.len(), 4);
    assert_eq!(left.atoms[&0].nbr1, right.atoms[&0].nbr1);
    assert_eq!(left.atoms[&0].nbr2, right.atoms[&0].nbr2);
    assert_eq!(left.atoms[&3].aid, 0);
    assert!(left.atoms[&3].fixed);
}

#[test]
fn fragment_merge_two_common_reflects_overcrowded_side_and_tracks_attachment() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut left = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    let mut right = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (3, [0.0, 1.0])]),
    )
    .unwrap();
    let mut common = left.find_common_atoms(&right);
    assert_eq!(common, vec![0, 1]);
    left.merge_with_common(&mut right, &mut common).unwrap();
    assert_eq!(left.atoms.len(), 4);
    assert!((left.atoms[&3].loc[0]).abs() < 1e-12);
    assert!((left.atoms[&3].loc[1] + 1.0).abs() < 1e-12);
    assert_ne!(left.atoms[&3].ccw, true);
    assert_eq!(left.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(left.atoms[&1].loc, [1.0, 0.0]);
}

#[test]
fn fragment_merge_three_common_uses_third_point_reflection() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (1, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut left = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    let mut right = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([
            (0, [0.0, 0.0]),
            (1, [1.0, 0.0]),
            (2, [0.0, -1.0]),
            (3, [1.0, -1.0]),
        ]),
    )
    .unwrap();
    let mut common = left.find_common_atoms(&right);
    assert_eq!(common, vec![0, 1, 2]);
    left.merge_with_common(&mut right, &mut common).unwrap();
    assert_eq!(left.atoms.len(), 4);
    assert!((left.atoms[&3].loc[0] - 1.0).abs() < 1e-12);
    assert!((left.atoms[&3].loc[1] - 1.0).abs() < 1e-12);
    assert_eq!(left.atoms[&3].aid, 0);
}

#[test]
fn fragment_merge_rejects_source_preconditions_without_mutating_inputs() {
    let graph = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let unrelated = topology(vec![AtomSpec::new(Element::C); 2], &[(0, 1)]);
    let rings = RingInfo::new(RingFindType::Fast, graph.atoms.len(), graph.bonds.len());
    let mut left = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut right = embedded_frag::EmbeddedFrag::from_single(1, &graph, &rings).unwrap();
    let before = left.atoms.clone();
    assert_eq!(
        left.merge_with_common(&mut right, &mut vec![]),
        Err(embedded_frag::FragmentError::NoCommonAtoms)
    );
    assert_eq!(left.atoms, before);
    let mut other_graph = embedded_frag::EmbeddedFrag::from_single(1, &unrelated, &rings).unwrap();
    assert_eq!(
        left.merge_no_common(&mut other_graph, 0, 1),
        Err(embedded_frag::FragmentError::MismatchedTopology)
    );
    assert_eq!(left.atoms, before);
}

#[test]
fn fragment_seed_expand_isolated_seed_has_no_attachment_work() {
    let graph = topology(vec![AtomSpec::new(Element::C)], &[]);
    let rings = RingInfo::new(RingFindType::Fast, 1, 0);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut nonring = Vec::new();
    let mut other = Vec::new();
    fragment.expand_fragment(&mut nonring, &mut other).unwrap();
    assert_eq!(fragment.atoms[&0].loc, [0.0, 0.0]);
    assert_eq!(fragment.atoms.len(), 1);
    assert!(fragment.attachment_points.is_empty());
}

#[test]
fn fragment_seed_expand_chain_consumes_nonring_rows_in_source_order() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = RingInfo::new(RingFindType::Fast, 3, 2);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut nonring = vec![2, 1];
    fragment
        .expand_fragment(&mut nonring, &mut Vec::new())
        .unwrap();
    assert!(nonring.is_empty());
    assert!(fragment.attachment_points.is_empty());
    assert_eq!(fragment.atoms.len(), 3);
    assert!((fragment.atoms[&1].loc[0] - 1.5 * (3.0_f64).sqrt() / 2.0).abs() < 1e-12);
    assert!((fragment.atoms[&1].loc[1] - 0.75).abs() < 1e-12);
    assert!((fragment.atoms[&2].loc[0] - 1.5 * (3.0_f64).sqrt()).abs() < 1e-12);
    assert!(fragment.atoms[&2].loc[1].abs() < 1e-12);
}

#[test]
fn fragment_seed_expand_branch_drains_attachment_and_preserves_bond_lengths() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 4, 3);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut nonring = vec![3, 1, 2];
    fragment
        .expand_fragment(&mut nonring, &mut Vec::new())
        .unwrap();
    assert!(nonring.is_empty());
    assert!(fragment.attachment_points.is_empty());
    assert_eq!(fragment.atoms.len(), 4);
    for aid in 1..4 {
        let point = fragment.atoms[&aid].loc;
        assert!((point[0].hypot(point[1]) - geometry::BOND_LEN).abs() < 1e-12);
    }
}

#[test]
fn fragment_seed_expand_bridge_searches_remaining_fragments() {
    let graph = topology(vec![AtomSpec::new(Element::C); 3], &[(0, 1), (1, 2)]);
    let rings = RingInfo::new(RingFindType::Fast, 3, 2);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut fragments = vec![embedded_frag::EmbeddedFrag::from_single(1, &graph, &rings).unwrap()];
    let mut nonring = vec![2];
    fragment
        .expand_fragment(&mut nonring, &mut fragments)
        .unwrap();
    assert!(fragments.is_empty());
    assert!(nonring.is_empty());
    assert_eq!(fragment.atoms.len(), 3);
    assert!(fragment.attachment_points.is_empty());
}

#[test]
fn fragment_seed_expand_shared_fragment_merges_before_queue_walk() {
    let graph = topology(
        vec![AtomSpec::new(Element::C); 4],
        &[(0, 1), (0, 2), (0, 3)],
    );
    let rings = RingInfo::new(RingFindType::Fast, 4, 3);
    let mut fragment = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 1.0])]),
    )
    .unwrap();
    let other = embedded_frag::EmbeddedFrag::from_coord_map(
        &graph,
        &rings,
        &geometry::PointMap::from([(0, [0.0, 0.0]), (3, [-1.0, 0.0])]),
    )
    .unwrap();
    let mut fragments = vec![other];
    fragment
        .expand_fragment(&mut Vec::new(), &mut fragments)
        .unwrap();
    assert!(fragments.is_empty());
    assert_eq!(fragment.atoms.len(), 4);
    assert!(fragment.attachment_points.is_empty());
    assert_eq!(fragment.atoms[&3].aid, 0); // pinned EmbeddedAtom assignment omits aid
}

#[test]
fn fragment_seed_expand_disconnected_fragment_is_left_for_component_handoff() {
    let graph = topology(vec![AtomSpec::new(Element::C); 4], &[(0, 1), (2, 3)]);
    let rings = RingInfo::new(RingFindType::Fast, 4, 2);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    let mut nonring = vec![3, 1];
    let mut fragments = vec![embedded_frag::EmbeddedFrag::from_single(2, &graph, &rings).unwrap()];
    fragment
        .expand_fragment(&mut nonring, &mut fragments)
        .unwrap();
    assert_eq!(fragment.atoms.len(), 2);
    assert_eq!(nonring, vec![3]);
    assert_eq!(fragments.len(), 1);
    assert_eq!(fragments[0].attachment_points, vec![2]);
    assert!(fragment.attachment_points.is_empty());
}

#[test]
fn fragment_seed_expand_reports_empty_attachment_invariant() {
    let graph = topology(vec![AtomSpec::new(Element::C)], &[]);
    let rings = RingInfo::new(RingFindType::Fast, 1, 0);
    let mut fragment = embedded_frag::EmbeddedFrag::from_single(0, &graph, &rings).unwrap();
    fragment.attachment_points.push(0);
    assert_eq!(
        fragment.expand_fragment(&mut Vec::new(), &mut Vec::new()),
        Err(embedded_frag::FragmentError::EmptyAttachment { atom: 0 })
    );
}
