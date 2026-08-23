use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_core::{
    AdditionalOutput, AtomPairFingerprintGenerator, AtomPairFingerprintParams,
    FingerprintFuncArguments, Molecule,
};

fn map(entries: &[(u64, i32)]) -> BTreeMap<u64, i32> {
    entries.iter().copied().collect()
}

fn bits(entries: &[u64]) -> BTreeSet<u64> {
    entries.iter().copied().collect()
}

fn all_atom_pair_output() -> AdditionalOutput {
    let mut output = AdditionalOutput::new();
    output.allocate_atom_counts();
    output.allocate_atom_to_bits();
    output.allocate_bit_info_map();
    output.allocate_atoms_per_bit();
    output
}

#[test]
fn default_generator_matches_rdkit_for_all_four_result_forms() {
    // Oracle: pinned RDKit 2026.03.1, GetAtomPairGenerator() on CCCO.
    let molecule = Molecule::from_smiles("CCCO").unwrap();
    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams::default())
        .expect("default AtomPair generator");

    let sparse_count = generator
        .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(sparse_count.size(), 1 << 23);
    assert_eq!(
        sparse_count.nonzero_elements(),
        &map(&[
            (558_113, 1),
            (558_114, 1),
            (558_145, 1),
            (1_590_307, 1),
            (1_590_337, 1),
            (1_590_338, 1),
        ])
    );

    let count = generator
        .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(count.size(), 2048);
    assert_eq!(
        count.nonzero_elements(),
        &map(&[(1310, 1), (1358, 2), (1375, 1), (1423, 1), (1692, 1)])
    );

    let sparse = generator
        .sparse_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(sparse.size(), 1 << 23);
    assert_eq!(
        sparse.on_bits(),
        &bits(&[
            7_918_712, 7_918_972, 7_920_240, 8_066_360, 8_066_361, 8_066_620
        ])
    );

    let explicit = generator
        .fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(explicit.n_bits(), 2048);
    assert_eq!(explicit.on_bits(), &[624, 1144, 1336, 1337, 1404, 1596]);
}

#[test]
fn custom_count_bounds_and_fold_collisions_match_rdkit_with_provenance() {
    // Oracle: pinned RDKit 2026.03.1 with countBounds=[1,3,5], fpSize=16.
    let molecule = Molecule::from_smiles("CCCCCC").unwrap();
    let params = AtomPairFingerprintParams {
        n_bits: 16,
        count_bounds: vec![1, 3, 5],
        ..Default::default()
    };
    let generator = AtomPairFingerprintGenerator::new(&params).unwrap();

    let sparse_count = generator
        .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        sparse_count.nonzero_elements(),
        &map(&[
            (541_733, 1),
            (558_113, 2),
            (558_114, 2),
            (558_115, 2),
            (558_116, 2),
            (558_145, 3),
            (558_146, 2),
            (558_147, 1),
        ])
    );

    let count = generator
        .count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        count.nonzero_elements(),
        &map(&[(11, 1), (12, 3), (13, 2), (14, 5), (15, 4)])
    );

    let sparse = generator
        .sparse_bit_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        sparse.on_bits(),
        &bits(&[
            1_747_161, 1_747_746, 1_747_941, 1_748_892, 1_749_087, 1_858_293, 1_858_482, 1_858_483,
            1_859_628,
        ])
    );

    let mut arguments = FingerprintFuncArguments {
        additional_output: Some(all_atom_pair_output()),
        ..Default::default()
    };
    let explicit = generator.fingerprint(&molecule, &mut arguments).unwrap();
    assert_eq!(explicit.on_bits(), &[0, 1, 2, 6, 7, 9, 10, 11]);

    let output = arguments.additional_output.unwrap();
    assert_eq!(output.atom_counts, Some(vec![5, 5, 5, 5, 5, 5]));
    assert_eq!(
        output.atom_to_bits,
        Some(vec![
            vec![0, 1, 2, 9, 10, 11],
            vec![0, 1, 2, 6, 7, 9, 10, 11],
            vec![0, 1, 2, 6, 7, 9, 10, 11],
            vec![0, 1, 2, 6, 7, 9, 10, 11],
            vec![0, 1, 2, 6, 7, 9, 10, 11],
            vec![0, 1, 2, 9, 10, 11],
        ])
    );
    let expected_pairs = BTreeMap::from([
        (
            0,
            vec![(0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (2, 5), (3, 4)],
        ),
        (
            1,
            vec![(0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (2, 5), (3, 4)],
        ),
        (
            2,
            vec![(0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (2, 5), (3, 4)],
        ),
        (6, vec![(1, 3), (1, 4), (2, 4)]),
        (7, vec![(1, 3), (1, 4), (2, 4)]),
        (9, vec![(0, 1), (0, 2), (0, 5), (3, 5), (4, 5)]),
        (10, vec![(0, 1), (0, 2), (0, 5), (3, 5), (4, 5)]),
        (11, vec![(0, 1), (0, 2), (0, 5), (3, 5), (4, 5)]),
    ]);
    assert_eq!(output.bit_info_map, Some(expected_pairs.clone()));
    assert_eq!(
        output.atoms_per_bit,
        Some(
            expected_pairs
                .into_iter()
                .map(|(bit, pairs)| {
                    (
                        bit,
                        pairs
                            .into_iter()
                            .map(|(first, second)| vec![first as usize, second as usize])
                            .collect(),
                    )
                })
                .collect()
        )
    );
}

#[test]
fn chirality_custom_invariants_and_atom_filters_match_rdkit() {
    // Oracle: pinned RDKit 2026.03.1 on C[C@H](O)F.
    let molecule = Molecule::from_smiles("C[C@H](O)F").unwrap();
    let chiral = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 64,
        use_chirality: true,
        count_simulation: false,
        ..Default::default()
    })
    .unwrap();
    let sparse = chiral
        .sparse_count_fingerprint(&molecule, &mut FingerprintFuncArguments::default())
        .unwrap();
    assert_eq!(
        sparse.nonzero_elements(),
        &map(&[
            (6_358_050, 1),
            (8_455_202, 1),
            (8_457_250, 1),
            (35_849_249, 1),
            (35_851_297, 1),
            (35_852_321, 1),
        ])
    );

    let restricted = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 64,
        min_distance: 1,
        max_distance: 2,
        use_chirality: true,
        count_simulation: false,
        ..Default::default()
    })
    .unwrap();
    let mut arguments = FingerprintFuncArguments {
        from_atoms: Some(vec![1]),
        ignore_atoms: Some(vec![3]),
        custom_atom_invariants: Some(vec![11, 22, 33, 44]),
        additional_output: Some(all_atom_pair_output()),
        ..Default::default()
    };
    let explicit = restricted.fingerprint(&molecule, &mut arguments).unwrap();
    assert_eq!(explicit.on_bits(), &[0, 37]);
    let output = arguments.additional_output.unwrap();
    assert_eq!(output.atom_counts, Some(vec![1, 2, 1, 0]));
    assert_eq!(
        output.atom_to_bits,
        Some(vec![vec![37], vec![37, 0], vec![0], vec![]])
    );
    assert_eq!(
        output.bit_info_map,
        Some(BTreeMap::from([(0, vec![(1, 2)]), (37, vec![(0, 1)])]))
    );
    assert_eq!(
        output.atoms_per_bit,
        Some(BTreeMap::from([
            (0, vec![vec![1, 2]]),
            (37, vec![vec![0, 1]]),
        ]))
    );
}

#[test]
fn three_dimensional_distance_selection_matches_rdkit() {
    // Oracle: pinned RDKit 2026.03.1 with use2D=false, min=2, max=3.
    let molecule = Molecule::from_mol_block(
        "\n     RDKit          3D\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    3.6000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0\n  2  3  1  0\nM  END\n",
    )
    .unwrap();
    let generator = AtomPairFingerprintGenerator::new(&AtomPairFingerprintParams {
        n_bits: 32,
        min_distance: 2,
        max_distance: 3,
        use_2d: false,
        count_simulation: false,
        conformer_id: 0,
        ..Default::default()
    })
    .unwrap();
    let arguments = || FingerprintFuncArguments {
        conf_id: 0,
        ..Default::default()
    };

    let sparse_count = generator
        .sparse_count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(
        sparse_count.nonzero_elements(),
        &map(&[(541_731, 1), (558_114, 1)])
    );
    let count = generator
        .count_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(count.nonzero_elements(), &map(&[(28, 1), (30, 1)]));
    let sparse = generator
        .sparse_bit_fingerprint(&molecule, &mut arguments())
        .unwrap();
    assert_eq!(sparse.on_bits(), &bits(&[6_173_982, 6_174_428]));
    let explicit = generator.fingerprint(&molecule, &mut arguments()).unwrap();
    assert_eq!(explicit.on_bits(), &[28, 30]);
}
