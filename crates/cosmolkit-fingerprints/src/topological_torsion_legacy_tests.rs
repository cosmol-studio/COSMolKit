#![allow(non_snake_case)]
use super::*;
use crate::test_support::TestMolecule;
fn getTopologicalTorsionFingerprint(
    molecule: &TestMolecule,
    target: u32,
    from: Option<&[u32]>,
    ignore: Option<&[u32]>,
    inv: Option<&[u32]>,
    chiral: bool,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    legacy_topological_torsion_sparse_count(
        &molecule.input(),
        &LegacyTopologicalTorsionParams {
            torsion_atom_count: target,
            include_chirality: chiral,
            ..Default::default()
        },
        &TopologicalTorsionCall {
            from_atoms: from,
            ignore_atoms: ignore,
            custom_atom_invariants: inv,
            ..Default::default()
        },
    )
}
fn getHashedTopologicalTorsionFingerprint(
    molecule: &TestMolecule,
    size: u32,
    target: u32,
    from: Option<&[u32]>,
    ignore: Option<&[u32]>,
    inv: Option<&[u32]>,
    chiral: bool,
) -> Result<SparseCountFingerprint, TopologicalTorsionError> {
    legacy_topological_torsion_count(
        &molecule.input(),
        &LegacyTopologicalTorsionParams {
            fp_size: size,
            torsion_atom_count: target,
            include_chirality: chiral,
            ..Default::default()
        },
        &TopologicalTorsionCall {
            from_atoms: from,
            ignore_atoms: ignore,
            custom_atom_invariants: inv,
            ..Default::default()
        },
    )
}
fn getHashedTopologicalTorsionFingerprintAsBitVect(
    molecule: &TestMolecule,
    size: u32,
    target: u32,
    from: Option<&[u32]>,
    ignore: Option<&[u32]>,
    inv: Option<&[u32]>,
    entry: u32,
    chiral: bool,
) -> Result<Fingerprint, TopologicalTorsionError> {
    legacy_topological_torsion_bits(
        &molecule.input(),
        &LegacyTopologicalTorsionParams {
            fp_size: size,
            torsion_atom_count: target,
            include_chirality: chiral,
            bits_per_entry: entry,
        },
        &TopologicalTorsionCall {
            from_atoms: from,
            ignore_atoms: ignore,
            custom_atom_invariants: inv,
            ..Default::default()
        },
    )
}
fn total(fingerprint: &SparseCountFingerprint) -> i32 {
    fingerprint.nonzero_elements().values().sum()
}

#[test]
fn legacy_unfolded_ids_counts_and_compatibility_size_match_pinned_rdkit() {
    // Pinned RDKit test1.cpp::testGitHubIssue25 fixes both unfolded ids.
    let molecule = TestMolecule::from_smiles("CCCCO").expect("molecule");
    let legacy = getTopologicalTorsionFingerprint(&molecule, 4, None, None, None, false)
        .expect("legacy unfolded");

    assert_eq!(legacy.length(), (1_u64 << 36) - 1);
    assert_eq!(total(&legacy), 2);
    assert_eq!(
        legacy
            .nonzero_elements()
            .iter()
            .map(|(&id, &count)| (id, count))
            .collect::<Vec<_>>(),
        vec![(4_437_590_048, 1), (12_893_306_913, 1)]
    );

    let modern = topological_torsion_sparse_count(
        &molecule.input(),
        &TopologicalTorsionParams::default(),
        &TopologicalTorsionCall::default(),
        None,
    )
    .unwrap();
    assert_eq!(legacy.nonzero_elements(), modern.nonzero_elements());
    assert_eq!(legacy.length() + 1, modern.length());
}

#[test]
fn legacy_hashed_counts_match_exact_pinned_ids_and_modern_count_core() {
    // Pinned RDKit test1.cpp::testGitHubIssue25 fixes the 1000-bit ids.
    let molecule = TestMolecule::from_smiles("CCCCO").expect("molecule");
    let legacy =
        getHashedTopologicalTorsionFingerprint(&molecule, 1000, 4, None, None, None, false)
            .expect("legacy hashed");
    assert_eq!(legacy.length(), 1000);
    assert_eq!(total(&legacy), 2);
    assert_eq!(
        legacy
            .nonzero_elements()
            .iter()
            .map(|(&id, &count)| (id, count))
            .collect::<Vec<_>>(),
        vec![(24, 1), (288, 1)]
    );

    let modern = count_helper(
        &molecule.input(),
        &TopologicalTorsionParams {
            fp_size: 1000,
            ..Default::default()
        },
        &TopologicalTorsionParams::default().common().unwrap(),
        &TopologicalTorsionCall::default(),
        1000,
        None,
        TorsionCodeMode::Modern,
    )
    .unwrap();
    assert_eq!(legacy, modern);
}

#[test]
fn legacy_bit_vector_uses_exact_four_bit_and_non_four_bit_thresholds() {
    let molecule = TestMolecule::from_smiles("CCCCCCCCCCCC").expect("long chain");
    let invariants = vec![7; molecule.num_atoms()];
    let counts = getHashedTopologicalTorsionFingerprint(
        &molecule,
        16,
        4,
        None,
        None,
        Some(&invariants),
        false,
    )
    .expect("block counts");
    assert_eq!(counts.nonzero_elements().len(), 1);
    let (&block, &count) = counts.nonzero_elements().first_key_value().unwrap();
    assert_eq!(count, 9);

    let four = getHashedTopologicalTorsionFingerprintAsBitVect(
        &molecule,
        64,
        4,
        None,
        None,
        Some(&invariants),
        4,
        false,
    )
    .expect("four-bit thresholds");
    let four_base = u32::try_from(block * 4).unwrap();
    assert_eq!(
        four.on_bits(),
        vec![four_base, four_base + 1, four_base + 2, four_base + 3]
    );

    let six = getHashedTopologicalTorsionFingerprintAsBitVect(
        &molecule,
        96,
        4,
        None,
        None,
        Some(&invariants),
        6,
        false,
    )
    .expect("non-four-bit thresholds");
    let six_base = u32::try_from(block * 6).unwrap();
    assert_eq!(six.on_bits(), (six_base..six_base + 6).collect::<Vec<_>>());
}

#[test]
fn legacy_bit_vector_floors_non_divisible_block_sizes_and_keeps_tail_clear() {
    let molecule = TestMolecule::from_smiles("CCCCCCCCCCCC").expect("long chain");
    let invariants = vec![7; molecule.num_atoms()];
    let bits = getHashedTopologicalTorsionFingerprintAsBitVect(
        &molecule,
        67,
        4,
        None,
        None,
        Some(&invariants),
        4,
        false,
    )
    .expect("non-divisible size");

    assert_eq!(bits.n_bits(), 67);
    assert_eq!(bits.on_bits().len(), 4);
    assert!(bits.on_bits().iter().all(|&bit| bit < 64));
}

#[test]
fn legacy_selections_custom_invariants_and_chirality_use_the_shared_core() {
    let chain = TestMolecule::from_smiles("CCCCC").expect("chain");
    let invariants = [10, 20, 30, 40, 50];
    let full = getHashedTopologicalTorsionFingerprint(
        &chain,
        2048,
        4,
        None,
        None,
        Some(&invariants),
        false,
    )
    .expect("full");
    let rooted = getHashedTopologicalTorsionFingerprint(
        &chain,
        2048,
        4,
        Some(&[0]),
        None,
        Some(&invariants),
        false,
    )
    .expect("rooted");
    let ignored = getHashedTopologicalTorsionFingerprint(
        &chain,
        2048,
        4,
        None,
        Some(&[2]),
        Some(&invariants),
        false,
    )
    .expect("ignored");
    assert_eq!(total(&full), 2);
    assert_eq!(total(&rooted), 1);
    assert!(ignored.nonzero_elements().is_empty());

    let clockwise = TestMolecule::from_smiles("CC[C@H](F)Cl").expect("clockwise");
    let anticlockwise = TestMolecule::from_smiles("CC[C@@H](F)Cl").expect("anticlockwise");
    let clockwise_achiral =
        getHashedTopologicalTorsionFingerprint(&clockwise, 2048, 4, None, None, None, false)
            .unwrap();
    let anticlockwise_achiral =
        getHashedTopologicalTorsionFingerprint(&anticlockwise, 2048, 4, None, None, None, false)
            .unwrap();
    assert_eq!(clockwise_achiral, anticlockwise_achiral);

    let clockwise_chiral =
        getHashedTopologicalTorsionFingerprint(&clockwise, 2048, 4, None, None, None, true)
            .unwrap();
    let anticlockwise_chiral =
        getHashedTopologicalTorsionFingerprint(&anticlockwise, 2048, 4, None, None, None, true)
            .unwrap();
    assert_ne!(clockwise_chiral, anticlockwise_chiral);
}

fn raw_sanitized(smiles: &str) -> TestMolecule {
    let record =
        cosmolkit_smiles::parse_smiles(smiles, &cosmolkit_smiles::SmilesParseParams::default())
            .unwrap();
    let assignment = cosmolkit_core::sanitize_topology(
        &record.topology,
        &cosmolkit_core::SanitizeParams::default(),
    )
    .unwrap();
    TestMolecule {
        topology: assignment.topology,
        properties: record.properties,
        coordinates: record.coordinates,
        valence: assignment.final_valence.unwrap(),
        rings: assignment.final_rings.unwrap(),
    }
}
#[test]
fn legacy_conditional_stereo_raw_matches_independent_pinned_native_three_forms() {
    {
        let raw = raw_sanitized("F[C@H](Cl)CC");
        let input = raw.input();
        let original_topology = raw.topology.clone();
        let original_properties = raw.properties.clone();
        assert!(raw.properties.prop("_StereochemDone").is_none());
        assert!(
            raw.topology
                .atoms
                .iter()
                .all(|a| a.prop("_CIPCode").is_none())
        );
        let params = LegacyTopologicalTorsionParams {
            include_chirality: true,
            ..Default::default()
        };
        let call = TopologicalTorsionCall::default();
        let sparse = legacy_topological_torsion_sparse_count(&input, &params, &call).unwrap();
        let count = legacy_topological_torsion_count(&input, &params, &call).unwrap();
        let bits = legacy_topological_torsion_bits(&input, &params, &call).unwrap();
        assert_eq!(
            sparse
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(1101797589024, 1), (2201309216800, 1)]
        );
        assert_eq!(
            count
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(1076, 1), (1204, 1)]
        );
        assert_eq!(bits.on_bits(), vec![208, 720]);
        assert_eq!(raw.topology, original_topology);
        assert_eq!(raw.properties, original_properties);
    }
    {
        let raw = raw_sanitized("F[C@@H](Cl)CC");
        let input = raw.input();
        let original_topology = raw.topology.clone();
        let original_properties = raw.properties.clone();
        assert!(raw.properties.prop("_StereochemDone").is_none());
        assert!(
            raw.topology
                .atoms
                .iter()
                .all(|a| a.prop("_CIPCode").is_none())
        );
        let params = LegacyTopologicalTorsionParams {
            include_chirality: true,
            ..Default::default()
        };
        let call = TopologicalTorsionCall::default();
        let sparse = legacy_topological_torsion_sparse_count(&input, &params, &call).unwrap();
        let count = legacy_topological_torsion_count(&input, &params, &call).unwrap();
        let bits = legacy_topological_torsion_bits(&input, &params, &call).unwrap();
        assert_eq!(
            sparse
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(1103945072672, 1), (2203456700448, 1)]
        );
        assert_eq!(
            count
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(1269, 1), (1397, 1)]
        );
        assert_eq!(bits.on_bits(), vec![980, 1492]);
        assert_eq!(raw.topology, original_topology);
        assert_eq!(raw.properties, original_properties);
    }
    {
        let raw = raw_sanitized("N[C@H](C)CC");
        let input = raw.input();
        let original_topology = raw.topology.clone();
        let original_properties = raw.properties.clone();
        assert!(raw.properties.prop("_StereochemDone").is_none());
        assert!(
            raw.topology
                .atoms
                .iter()
                .all(|a| a.prop("_CIPCode").is_none())
        );
        let params = LegacyTopologicalTorsionParams {
            include_chirality: true,
            ..Default::default()
        };
        let call = TopologicalTorsionCall::default();
        let sparse = legacy_topological_torsion_sparse_count(&input, &params, &call).unwrap();
        let count = legacy_topological_torsion_count(&input, &params, &call).unwrap();
        let bits = legacy_topological_torsion_bits(&input, &params, &call).unwrap();
        assert_eq!(
            sparse
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(277163868192, 1), (552041775136, 1)]
        );
        assert_eq!(
            count
                .nonzero_elements()
                .iter()
                .map(|(&k, &v)| (k, v))
                .collect::<Vec<_>>(),
            vec![(1140, 1), (1940, 1)]
        );
        assert_eq!(bits.on_bits(), vec![464, 1616]);
        assert_eq!(raw.topology, original_topology);
        assert_eq!(raw.properties, original_properties);
    }
}
#[test]
fn legacy_stereo_done_presence_custom_invariants_and_modern_original_provider_control() {
    let mut raw = raw_sanitized("F[C@H](Cl)CC");
    let custom = [17, 18, 19, 20, 21];
    for done in [None, Some("0"), Some("1")] {
        if let Some(v) = done {
            raw.properties.set_prop("_StereochemDone", v).unwrap();
        }
        let input = raw.input();
        let prepared = crate::prepared::prepare_morgan_environment(
            input.topology,
            input.properties,
            input.valence,
            input.rings,
            true,
        )
        .unwrap();
        let prepared_input = AtomPairPreparedInput {
            topology: prepared.topology(),
            properties: prepared.properties(),
            ..input
        };
        let params = LegacyTopologicalTorsionParams {
            include_chirality: true,
            ..Default::default()
        };
        for inv in [None, Some(custom.as_slice())] {
            let call = TopologicalTorsionCall {
                custom_atom_invariants: inv,
                ..Default::default()
            };
            assert_eq!(
                legacy_topological_torsion_sparse_count(&input, &params, &call).unwrap(),
                legacy_topological_torsion_sparse_count(&prepared_input, &params, &call).unwrap()
            );
            assert_eq!(
                legacy_topological_torsion_count(&input, &params, &call).unwrap(),
                legacy_topological_torsion_count(&prepared_input, &params, &call).unwrap()
            );
            assert_eq!(
                legacy_topological_torsion_bits(&input, &params, &call).unwrap(),
                legacy_topological_torsion_bits(&prepared_input, &params, &call).unwrap()
            );
        }
        let modern_params = TopologicalTorsionParams {
            include_chirality: true,
            ..Default::default()
        };
        let call = TopologicalTorsionCall::default();
        let modern = topological_torsion_sparse_count(&input, &modern_params, &call, None).unwrap();
        // Modern source computes invariants from the original molecule even
        // when environments are conditionally prepared. Keep this control.
        assert_eq!(
            modern
                .nonzero_elements()
                .keys()
                .copied()
                .collect::<Vec<_>>(),
            vec![1099650105376, 2199161733152]
        );
        if done.is_some() {
            assert_eq!(
                legacy_topological_torsion_sparse_count(&input, &params, &call)
                    .unwrap()
                    .nonzero_elements(),
                modern.nonzero_elements()
            );
        }
    }
}
