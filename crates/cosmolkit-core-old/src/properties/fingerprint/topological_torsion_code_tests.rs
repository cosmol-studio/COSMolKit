use super::{FingerprintError, get_topological_torsion_code, get_topological_torsion_hash};

#[test]
fn forward_and_reverse_paths_have_identical_codes_and_hashes() {
    let forward = [1, 2, 3, 4];
    let reverse = [4, 3, 2, 1];

    assert_eq!(
        get_topological_torsion_code(&forward, false).unwrap(),
        get_topological_torsion_code(&reverse, false).unwrap()
    );
    assert_eq!(
        get_topological_torsion_hash(&forward).unwrap(),
        get_topological_torsion_hash(&reverse).unwrap()
    );
    assert_eq!(get_topological_torsion_hash(&forward).unwrap(), 0x481a_d382);
}

#[test]
fn palindrome_keeps_forward_packing_order() {
    let palindrome = [7, 3, 3, 7];
    let expected = 7u64 | (3u64 << 9) | (3u64 << 18) | (7u64 << 27);

    assert_eq!(
        get_topological_torsion_code(&palindrome, false).unwrap(),
        expected
    );
}

#[test]
fn first_difference_after_equal_endpoints_selects_canonical_direction() {
    let path = [4, 9, 1, 4];
    let reversed = [4, 1, 9, 4];
    let expected = 4u64 | (1u64 << 9) | (9u64 << 18) | (4u64 << 27);

    assert_eq!(
        get_topological_torsion_code(&path, false).unwrap(),
        expected
    );
    assert_eq!(
        get_topological_torsion_code(&path, false).unwrap(),
        get_topological_torsion_code(&reversed, false).unwrap()
    );
}

#[test]
fn chirality_selects_eleven_bit_instead_of_nine_bit_fields() {
    let path = [1, 2, 3, 4];
    let expected_without_chirality = 1u64 | (2u64 << 9) | (3u64 << 18) | (4u64 << 27);
    let expected_with_chirality = 1u64 | (2u64 << 11) | (3u64 << 22) | (4u64 << 33);

    assert_eq!(
        get_topological_torsion_code(&path, false).unwrap(),
        expected_without_chirality
    );
    assert_eq!(
        get_topological_torsion_code(&path, true).unwrap(),
        expected_with_chirality
    );
}

#[test]
fn zero_path_codes_preserve_zero_packing_and_source_hashing() {
    let path = [0, 0, 0, 0];

    assert_eq!(get_topological_torsion_code(&path, false).unwrap(), 0);
    assert_eq!(get_topological_torsion_code(&path, true).unwrap(), 0);
    assert_eq!(get_topological_torsion_hash(&path).unwrap(), 0x484d_221a);
}

#[test]
fn maximum_defined_shift_widths_are_accepted() {
    let nine_bit_path = [1; 8];
    let eleven_bit_path = [1; 6];
    let expected_nine_bit = (0..8).fold(0u64, |packed, index| packed | (1u64 << (9 * index)));
    let expected_eleven_bit = (0..6).fold(0u64, |packed, index| packed | (1u64 << (11 * index)));

    assert_eq!(
        get_topological_torsion_code(&nine_bit_path, false).unwrap(),
        expected_nine_bit
    );
    assert_eq!(
        get_topological_torsion_code(&eleven_bit_path, true).unwrap(),
        expected_eleven_bit
    );
}

#[test]
fn source_undefined_empty_and_excessive_shifts_return_structured_errors() {
    assert!(matches!(
        get_topological_torsion_code(&[], false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        get_topological_torsion_hash(&[]),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        get_topological_torsion_code(&[1; 9], false),
        Err(FingerprintError::InvalidArguments { .. })
    ));
    assert!(matches!(
        get_topological_torsion_code(&[1; 7], true),
        Err(FingerprintError::InvalidArguments { .. })
    ));
}

#[test]
fn boost_hash_combine_overflow_matches_the_pinned_u32_result() {
    assert_eq!(
        get_topological_torsion_hash(&[u32::MAX; 4]).unwrap(),
        0x4841_0b19
    );
}
