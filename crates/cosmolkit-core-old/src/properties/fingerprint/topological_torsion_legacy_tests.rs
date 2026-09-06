#![allow(deprecated)]

use super::{
    FingerprintError, FingerprintFuncArguments, TopologicalTorsionArguments,
    getHashedTopologicalTorsionFingerprint, getHashedTopologicalTorsionFingerprintAsBitVect,
    getTopologicalTorsionFingerprint, getTopologicalTorsionGenerator,
    getTopologicalTorsionGeneratorWithParams,
};
use crate::Molecule;

fn total(fingerprint: &super::SparseCountFingerprint) -> i32 {
    fingerprint.nonzero_elements().values().sum()
}

#[test]
fn legacy_unfolded_ids_counts_and_compatibility_size_match_pinned_rdkit() {
    // Pinned RDKit test1.cpp::testGitHubIssue25 fixes both unfolded ids.
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let legacy = getTopologicalTorsionFingerprint(&molecule, 4, None, None, None, false)
        .expect("legacy unfolded");

    assert_eq!(legacy.size(), (1_u64 << 36) - 1);
    assert_eq!(total(&legacy), 2);
    assert_eq!(
        legacy
            .nonzero_elements()
            .iter()
            .map(|(&id, &count)| (id, count))
            .collect::<Vec<_>>(),
        vec![(4_437_590_048, 1), (12_893_306_913, 1)]
    );

    let modern =
        getTopologicalTorsionGenerator(&TopologicalTorsionArguments::default(), None, true)
            .expect("modern generator")
            .getSparseCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
            .expect("modern unfolded");
    assert_eq!(legacy.nonzero_elements(), modern.nonzero_elements());
    assert_eq!(legacy.size() + 1, modern.size());
}

#[test]
fn legacy_hashed_counts_match_exact_pinned_ids_and_modern_count_core() {
    // Pinned RDKit test1.cpp::testGitHubIssue25 fixes the 1000-bit ids.
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    let legacy =
        getHashedTopologicalTorsionFingerprint(&molecule, 1000, 4, None, None, None, false)
            .expect("legacy hashed");
    assert_eq!(legacy.size(), 1000);
    assert_eq!(total(&legacy), 2);
    assert_eq!(
        legacy
            .nonzero_elements()
            .iter()
            .map(|(&id, &count)| (id, count))
            .collect::<Vec<_>>(),
        vec![(24, 1), (288, 1)]
    );

    let modern = getTopologicalTorsionGeneratorWithParams(
        false,
        4,
        None,
        true,
        1000,
        vec![1, 2, 4, 8],
        true,
    )
    .expect("modern generator")
    .getCountFingerprint(&molecule, &mut FingerprintFuncArguments::default())
    .expect("modern hashed");
    assert_eq!(legacy, modern);
}

#[test]
fn legacy_bit_vector_uses_exact_four_bit_and_non_four_bit_thresholds() {
    let molecule = Molecule::from_smiles("CCCCCCCCCCCC").expect("long chain");
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
    let four_base = usize::try_from(block * 4).unwrap();
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
    let six_base = usize::try_from(block * 6).unwrap();
    assert_eq!(six.on_bits(), (six_base..six_base + 6).collect::<Vec<_>>());
}

#[test]
fn legacy_bit_vector_floors_non_divisible_block_sizes_and_keeps_tail_clear() {
    let molecule = Molecule::from_smiles("CCCCCCCCCCCC").expect("long chain");
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
    let chain = Molecule::from_smiles("CCCCC").expect("chain");
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

    let clockwise = Molecule::from_smiles("CC[C@H](F)Cl").expect("clockwise");
    let anticlockwise = Molecule::from_smiles("CC[C@@H](F)Cl").expect("anticlockwise");
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

#[test]
fn legacy_invalid_parameters_return_structured_errors() {
    let molecule = Molecule::from_smiles("CCCCO").expect("molecule");
    assert!(matches!(
        getHashedTopologicalTorsionFingerprint(&molecule, 0, 4, None, None, None, false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        getHashedTopologicalTorsionFingerprintAsBitVect(
            &molecule, 64, 4, None, None, None, 0, false,
        ),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        getHashedTopologicalTorsionFingerprintAsBitVect(
            &molecule, 3, 4, None, None, None, 4, false,
        ),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        getTopologicalTorsionFingerprint(&molecule, 4, None, None, Some(&[1, 2]), false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        getHashedTopologicalTorsionFingerprint(&molecule, 2048, 8, None, None, None, false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}
