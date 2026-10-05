//! Tests for the complete FP02 supplied-code numeric helpers.
//! Literal expectations derive from pinned RDKit integer algebra.
//! Empty/undefined-shift/nonterminating-path errors follow ROOT-approved COSMolKit safety policy.
use cosmolkit_fingerprints::{
    FingerprintError, atom_pair_code, topological_torsion_code, topological_torsion_hash,
};
fn pair_group(group: &str) {
    let mut observations = 0usize;
    for &(g, a, z, distance, chir, expected) in PAIR {
        if g == group {
            observations += 1;
            assert_eq!(
                atom_pair_code(a, z, distance, chir),
                Ok(expected),
                "{g}/{a}/{z}/{distance}/{chir}"
            );
        }
    }
    assert!(
        observations > 0,
        "frozen case group {group} must execute observations"
    );
}
fn code_group(group: &str) {
    let mut observations = 0usize;
    for &(g, codes, chir, expected) in CODE {
        if g == group {
            observations += 1;
            let before = codes.to_vec();
            assert_eq!(
                topological_torsion_code(codes, chir),
                Ok(expected),
                "{g}/{codes:?}/{chir}"
            );
            assert_eq!(codes, before.as_slice());
        }
    }
    assert!(
        observations > 0,
        "frozen case group {group} must execute observations"
    );
}
fn hash_group(group: &str) {
    let mut observations = 0usize;
    for &(g, codes, expected) in HASH {
        if g == group {
            observations += 1;
            let before = codes.to_vec();
            assert_eq!(
                topological_torsion_hash(codes),
                Ok(expected),
                "{g}/{codes:?}"
            );
            assert_eq!(codes, before.as_slice());
        }
    }
    assert!(
        observations > 0,
        "frozen case group {group} must execute observations"
    );
}

#[test]
fn packed_pair_all_source_distance_values_both_booleans() {
    pair_group("distance_0_through_30");
}

#[test]
fn packed_pair_minmax_precedes_truncation_and_is_unsigned() {
    pair_group("minmax_unsigned_before_truncate");
}

#[test]
fn packed_pair_equal_codes_and_full_width() {
    pair_group("equal_codes_full_width");
}

#[test]
fn packed_pair_all32_code_bits_without_mask() {
    pair_group("all32_bit_positions");
}

#[test]
fn packed_pair_chirality_changes_only_second_shift() {
    pair_group("chirality_and_highbit_overlap");
}

#[test]
fn packed_pair_original_upstream_scalar_ccccc_codes() {
    for &(a, z, d, e) in &[
        (33, 34, 1, 558113u32),
        (34, 33, 1, 558113),
        (33, 34, 2, 558114),
        (34, 33, 2, 558114),
    ] {
        assert_eq!(atom_pair_code(a, z, d, false), Ok(e));
    }
}

#[test]
fn packed_pair_source_precondition_strict_31() {
    for chir in [false, true] {
        for d in [31, 32, 0x80000000, u32::MAX] {
            assert_eq!(
                atom_pair_code(u32::MAX, 0, d, chir),
                Err(FingerprintError::PreconditionViolation {
                    what: "dist too long"
                })
            );
        }
    }
}

#[test]
fn packed_pair_distance_errors_do_not_mutate_input_tuple() {
    let input = (u32::MAX, 0u32, 31u32, true);
    let before = input;
    assert_eq!(
        atom_pair_code(input.0, input.1, input.2, input.3),
        Err(FingerprintError::PreconditionViolation {
            what: "dist too long"
        })
    );
    assert_eq!(input, before);
}

#[test]
fn packed_torsion_every_defined_length_both_booleans() {
    code_group("all_defined_lengths");
}

#[test]
fn packed_torsion_paired_end_outer_and_inner_decisions() {
    code_group("paired_end_order_branches");
}

#[test]
fn packed_torsion_equal_palindromes_keep_source_order() {
    code_group("palindrome_equal_pairs");
}

#[test]
fn packed_torsion_one_element_retains_all_u32_bits() {
    code_group("one_element_full_u32");
}

#[test]
fn packed_torsion_OR_overlap_is_not_masked_or_added() {
    code_group("OR_overlap_and_full_width");
}

#[test]
fn packed_torsion_last_defined_shifts63_and55() {
    code_group("last_defined_shift");
}

#[test]
fn packed_torsion_original_upstream_butane_literal() {
    assert_eq!(
        topological_torsion_code(&[32, 32, 32, 32], false),
        Ok(4303372320u64)
    );
}

#[test]
fn packed_torsion_empty_is_source_undefined_safe_proposal() {
    for chir in [false, true] {
        assert_eq!(
            topological_torsion_code(&[], chir),
            Err(FingerprintError::InvalidArguments {
                reason: "empty topological torsion code path is undefined in pinned source"
            })
        );
    }
}

#[test]
fn packed_torsion_zero_operands_do_not_define_out_of_width_shifts() {
    for (length, chir) in [(9, false), (10, false), (7, true), (8, true)] {
        for value in [0u32, u32::MAX] {
            let codes = vec![value; length];
            let before = codes.clone();
            assert_eq!(
                topological_torsion_code(&codes, chir),
                Err(FingerprintError::InvalidArguments {
                    reason: "topological torsion code executes a shift outside source uint64 width"
                })
            );
            assert_eq!(codes, before);
        }
    }
}

#[test]
fn packed_torsion_borrowed_codes_preserved_on_success_and_errors() {
    for codes in [
        vec![],
        vec![1],
        vec![1, 9, 2, 1],
        vec![0; 7],
        vec![u32::MAX; 9],
    ] {
        let before = codes.clone();
        for chir in [false, true] {
            let _ = topological_torsion_code(&codes, chir);
            assert_eq!(codes, before);
        }
    }
}

#[test]
fn packed_hash_one_element_wrap_and_unsigned_identity() {
    hash_group("one_element_u32_carry");
}

#[test]
fn packed_hash_reversal_and_paired_inner_order() {
    hash_group("paired_end_order_branches");
}

#[test]
fn packed_hash_palindrome_and_full_order_not_sort() {
    hash_group("palindrome_equal_pairs");
}

#[test]
fn packed_hash_ordinary_lengths_exceed_packing_limits() {
    hash_group("all_ordinary_lengths_no_packing_cap");
}

#[test]
fn packed_hash_all32_bits_and_combination_wrap() {
    hash_group("all32_bit_positions");
    hash_group("full_width_overflow");
}

#[test]
fn packed_hash_empty_is_source_undefined_safe_proposal() {
    assert_eq!(
        topological_torsion_hash(&[]),
        Err(FingerprintError::InvalidArguments {
            reason: "empty topological torsion hash path is undefined in pinned source"
        })
    );
}

#[test]
fn packed_hash_borrowed_codes_preserved() {
    for &(group, codes, expected) in HASH {
        let before = codes.to_vec();
        assert_eq!(topological_torsion_hash(codes), Ok(expected), "{group}");
        assert_eq!(topological_torsion_hash(codes), Ok(expected));
        assert_eq!(codes, before.as_slice());
    }
}

#[test]
fn packed_codes_exact_domain_signatures_existing_error() {
    let _: fn(u32, u32, u32, bool) -> Result<u32, FingerprintError> = atom_pair_code;
    let _: fn(&[u32], bool) -> Result<u64, FingerprintError> = topological_torsion_code;
    let _: fn(&[u32]) -> Result<u32, FingerprintError> = topological_torsion_hash;
}

const PAIR: &[(&str, u32, u32, u32, bool, u32)] = &[
    ("distance_0_through_30", 0u32, 0u32, 0, false, 0u32),
    ("distance_0_through_30", 1u32, 2u32, 0, false, 32800u32),
    ("distance_0_through_30", 33u32, 34u32, 0, false, 558112u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        0,
        false,
        4294950912u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 1, false, 1u32),
    ("distance_0_through_30", 1u32, 2u32, 1, false, 32801u32),
    ("distance_0_through_30", 33u32, 34u32, 1, false, 558113u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        1,
        false,
        4294950913u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 2, false, 2u32),
    ("distance_0_through_30", 1u32, 2u32, 2, false, 32802u32),
    ("distance_0_through_30", 33u32, 34u32, 2, false, 558114u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        2,
        false,
        4294950914u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 3, false, 3u32),
    ("distance_0_through_30", 1u32, 2u32, 3, false, 32803u32),
    ("distance_0_through_30", 33u32, 34u32, 3, false, 558115u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        3,
        false,
        4294950915u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 4, false, 4u32),
    ("distance_0_through_30", 1u32, 2u32, 4, false, 32804u32),
    ("distance_0_through_30", 33u32, 34u32, 4, false, 558116u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        4,
        false,
        4294950916u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 5, false, 5u32),
    ("distance_0_through_30", 1u32, 2u32, 5, false, 32805u32),
    ("distance_0_through_30", 33u32, 34u32, 5, false, 558117u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        5,
        false,
        4294950917u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 6, false, 6u32),
    ("distance_0_through_30", 1u32, 2u32, 6, false, 32806u32),
    ("distance_0_through_30", 33u32, 34u32, 6, false, 558118u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        6,
        false,
        4294950918u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 7, false, 7u32),
    ("distance_0_through_30", 1u32, 2u32, 7, false, 32807u32),
    ("distance_0_through_30", 33u32, 34u32, 7, false, 558119u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        7,
        false,
        4294950919u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 8, false, 8u32),
    ("distance_0_through_30", 1u32, 2u32, 8, false, 32808u32),
    ("distance_0_through_30", 33u32, 34u32, 8, false, 558120u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        8,
        false,
        4294950920u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 9, false, 9u32),
    ("distance_0_through_30", 1u32, 2u32, 9, false, 32809u32),
    ("distance_0_through_30", 33u32, 34u32, 9, false, 558121u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        9,
        false,
        4294950921u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 10, false, 10u32),
    ("distance_0_through_30", 1u32, 2u32, 10, false, 32810u32),
    ("distance_0_through_30", 33u32, 34u32, 10, false, 558122u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        10,
        false,
        4294950922u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 11, false, 11u32),
    ("distance_0_through_30", 1u32, 2u32, 11, false, 32811u32),
    ("distance_0_through_30", 33u32, 34u32, 11, false, 558123u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        11,
        false,
        4294950923u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 12, false, 12u32),
    ("distance_0_through_30", 1u32, 2u32, 12, false, 32812u32),
    ("distance_0_through_30", 33u32, 34u32, 12, false, 558124u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        12,
        false,
        4294950924u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 13, false, 13u32),
    ("distance_0_through_30", 1u32, 2u32, 13, false, 32813u32),
    ("distance_0_through_30", 33u32, 34u32, 13, false, 558125u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        13,
        false,
        4294950925u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 14, false, 14u32),
    ("distance_0_through_30", 1u32, 2u32, 14, false, 32814u32),
    ("distance_0_through_30", 33u32, 34u32, 14, false, 558126u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        14,
        false,
        4294950926u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 15, false, 15u32),
    ("distance_0_through_30", 1u32, 2u32, 15, false, 32815u32),
    ("distance_0_through_30", 33u32, 34u32, 15, false, 558127u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        15,
        false,
        4294950927u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 16, false, 16u32),
    ("distance_0_through_30", 1u32, 2u32, 16, false, 32816u32),
    ("distance_0_through_30", 33u32, 34u32, 16, false, 558128u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        16,
        false,
        4294950928u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 17, false, 17u32),
    ("distance_0_through_30", 1u32, 2u32, 17, false, 32817u32),
    ("distance_0_through_30", 33u32, 34u32, 17, false, 558129u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        17,
        false,
        4294950929u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 18, false, 18u32),
    ("distance_0_through_30", 1u32, 2u32, 18, false, 32818u32),
    ("distance_0_through_30", 33u32, 34u32, 18, false, 558130u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        18,
        false,
        4294950930u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 19, false, 19u32),
    ("distance_0_through_30", 1u32, 2u32, 19, false, 32819u32),
    ("distance_0_through_30", 33u32, 34u32, 19, false, 558131u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        19,
        false,
        4294950931u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 20, false, 20u32),
    ("distance_0_through_30", 1u32, 2u32, 20, false, 32820u32),
    ("distance_0_through_30", 33u32, 34u32, 20, false, 558132u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        20,
        false,
        4294950932u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 21, false, 21u32),
    ("distance_0_through_30", 1u32, 2u32, 21, false, 32821u32),
    ("distance_0_through_30", 33u32, 34u32, 21, false, 558133u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        21,
        false,
        4294950933u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 22, false, 22u32),
    ("distance_0_through_30", 1u32, 2u32, 22, false, 32822u32),
    ("distance_0_through_30", 33u32, 34u32, 22, false, 558134u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        22,
        false,
        4294950934u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 23, false, 23u32),
    ("distance_0_through_30", 1u32, 2u32, 23, false, 32823u32),
    ("distance_0_through_30", 33u32, 34u32, 23, false, 558135u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        23,
        false,
        4294950935u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 24, false, 24u32),
    ("distance_0_through_30", 1u32, 2u32, 24, false, 32824u32),
    ("distance_0_through_30", 33u32, 34u32, 24, false, 558136u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        24,
        false,
        4294950936u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 25, false, 25u32),
    ("distance_0_through_30", 1u32, 2u32, 25, false, 32825u32),
    ("distance_0_through_30", 33u32, 34u32, 25, false, 558137u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        25,
        false,
        4294950937u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 26, false, 26u32),
    ("distance_0_through_30", 1u32, 2u32, 26, false, 32826u32),
    ("distance_0_through_30", 33u32, 34u32, 26, false, 558138u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        26,
        false,
        4294950938u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 27, false, 27u32),
    ("distance_0_through_30", 1u32, 2u32, 27, false, 32827u32),
    ("distance_0_through_30", 33u32, 34u32, 27, false, 558139u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        27,
        false,
        4294950939u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 28, false, 28u32),
    ("distance_0_through_30", 1u32, 2u32, 28, false, 32828u32),
    ("distance_0_through_30", 33u32, 34u32, 28, false, 558140u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        28,
        false,
        4294950940u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 29, false, 29u32),
    ("distance_0_through_30", 1u32, 2u32, 29, false, 32829u32),
    ("distance_0_through_30", 33u32, 34u32, 29, false, 558141u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        29,
        false,
        4294950941u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 30, false, 30u32),
    ("distance_0_through_30", 1u32, 2u32, 30, false, 32830u32),
    ("distance_0_through_30", 33u32, 34u32, 30, false, 558142u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        30,
        false,
        4294950942u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        1u32,
        2u32,
        7,
        false,
        32807u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2u32,
        1u32,
        7,
        false,
        32807u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        33u32,
        34u32,
        7,
        false,
        558119u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        34u32,
        33u32,
        7,
        false,
        558119u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        4294967295u32,
        1u32,
        7,
        false,
        4294950951u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        1u32,
        4294967295u32,
        7,
        false,
        4294950951u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2147483648u32,
        2147483647u32,
        7,
        false,
        4294967271u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2147483647u32,
        2147483648u32,
        7,
        false,
        4294967271u32,
    ),
    ("equal_codes_full_width", 0u32, 0u32, 30, false, 30u32),
    ("equal_codes_full_width", 1u32, 1u32, 30, false, 16446u32),
    (
        "equal_codes_full_width",
        511u32,
        511u32,
        30,
        false,
        8388606u32,
    ),
    (
        "equal_codes_full_width",
        512u32,
        512u32,
        30,
        false,
        8405022u32,
    ),
    (
        "equal_codes_full_width",
        1023u32,
        1023u32,
        30,
        false,
        16777214u32,
    ),
    (
        "equal_codes_full_width",
        1024u32,
        1024u32,
        30,
        false,
        16810014u32,
    ),
    (
        "equal_codes_full_width",
        2047u32,
        2047u32,
        30,
        false,
        33554430u32,
    ),
    (
        "equal_codes_full_width",
        2048u32,
        2048u32,
        30,
        false,
        33619998u32,
    ),
    (
        "equal_codes_full_width",
        4294967295u32,
        4294967295u32,
        30,
        false,
        4294967294u32,
    ),
    ("all32_bit_positions", 0u32, 1u32, 3, false, 16387u32),
    ("all32_bit_positions", 1u32, 0u32, 3, false, 16387u32),
    ("all32_bit_positions", 1u32, 1u32, 3, false, 16419u32),
    ("all32_bit_positions", 0u32, 2u32, 3, false, 32771u32),
    ("all32_bit_positions", 2u32, 0u32, 3, false, 32771u32),
    ("all32_bit_positions", 2u32, 2u32, 3, false, 32835u32),
    ("all32_bit_positions", 0u32, 4u32, 3, false, 65539u32),
    ("all32_bit_positions", 4u32, 0u32, 3, false, 65539u32),
    ("all32_bit_positions", 4u32, 4u32, 3, false, 65667u32),
    ("all32_bit_positions", 0u32, 8u32, 3, false, 131075u32),
    ("all32_bit_positions", 8u32, 0u32, 3, false, 131075u32),
    ("all32_bit_positions", 8u32, 8u32, 3, false, 131331u32),
    ("all32_bit_positions", 0u32, 16u32, 3, false, 262147u32),
    ("all32_bit_positions", 16u32, 0u32, 3, false, 262147u32),
    ("all32_bit_positions", 16u32, 16u32, 3, false, 262659u32),
    ("all32_bit_positions", 0u32, 32u32, 3, false, 524291u32),
    ("all32_bit_positions", 32u32, 0u32, 3, false, 524291u32),
    ("all32_bit_positions", 32u32, 32u32, 3, false, 525315u32),
    ("all32_bit_positions", 0u32, 64u32, 3, false, 1048579u32),
    ("all32_bit_positions", 64u32, 0u32, 3, false, 1048579u32),
    ("all32_bit_positions", 64u32, 64u32, 3, false, 1050627u32),
    ("all32_bit_positions", 0u32, 128u32, 3, false, 2097155u32),
    ("all32_bit_positions", 128u32, 0u32, 3, false, 2097155u32),
    ("all32_bit_positions", 128u32, 128u32, 3, false, 2101251u32),
    ("all32_bit_positions", 0u32, 256u32, 3, false, 4194307u32),
    ("all32_bit_positions", 256u32, 0u32, 3, false, 4194307u32),
    ("all32_bit_positions", 256u32, 256u32, 3, false, 4202499u32),
    ("all32_bit_positions", 0u32, 512u32, 3, false, 8388611u32),
    ("all32_bit_positions", 512u32, 0u32, 3, false, 8388611u32),
    ("all32_bit_positions", 512u32, 512u32, 3, false, 8404995u32),
    ("all32_bit_positions", 0u32, 1024u32, 3, false, 16777219u32),
    ("all32_bit_positions", 1024u32, 0u32, 3, false, 16777219u32),
    (
        "all32_bit_positions",
        1024u32,
        1024u32,
        3,
        false,
        16809987u32,
    ),
    ("all32_bit_positions", 0u32, 2048u32, 3, false, 33554435u32),
    ("all32_bit_positions", 2048u32, 0u32, 3, false, 33554435u32),
    (
        "all32_bit_positions",
        2048u32,
        2048u32,
        3,
        false,
        33619971u32,
    ),
    ("all32_bit_positions", 0u32, 4096u32, 3, false, 67108867u32),
    ("all32_bit_positions", 4096u32, 0u32, 3, false, 67108867u32),
    (
        "all32_bit_positions",
        4096u32,
        4096u32,
        3,
        false,
        67239939u32,
    ),
    ("all32_bit_positions", 0u32, 8192u32, 3, false, 134217731u32),
    ("all32_bit_positions", 8192u32, 0u32, 3, false, 134217731u32),
    (
        "all32_bit_positions",
        8192u32,
        8192u32,
        3,
        false,
        134479875u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        16384u32,
        3,
        false,
        268435459u32,
    ),
    (
        "all32_bit_positions",
        16384u32,
        0u32,
        3,
        false,
        268435459u32,
    ),
    (
        "all32_bit_positions",
        16384u32,
        16384u32,
        3,
        false,
        268959747u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        32768u32,
        3,
        false,
        536870915u32,
    ),
    (
        "all32_bit_positions",
        32768u32,
        0u32,
        3,
        false,
        536870915u32,
    ),
    (
        "all32_bit_positions",
        32768u32,
        32768u32,
        3,
        false,
        537919491u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        65536u32,
        3,
        false,
        1073741827u32,
    ),
    (
        "all32_bit_positions",
        65536u32,
        0u32,
        3,
        false,
        1073741827u32,
    ),
    (
        "all32_bit_positions",
        65536u32,
        65536u32,
        3,
        false,
        1075838979u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        131072u32,
        3,
        false,
        2147483651u32,
    ),
    (
        "all32_bit_positions",
        131072u32,
        0u32,
        3,
        false,
        2147483651u32,
    ),
    (
        "all32_bit_positions",
        131072u32,
        131072u32,
        3,
        false,
        2151677955u32,
    ),
    ("all32_bit_positions", 0u32, 262144u32, 3, false, 3u32),
    ("all32_bit_positions", 262144u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        262144u32,
        262144u32,
        3,
        false,
        8388611u32,
    ),
    ("all32_bit_positions", 0u32, 524288u32, 3, false, 3u32),
    ("all32_bit_positions", 524288u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        524288u32,
        524288u32,
        3,
        false,
        16777219u32,
    ),
    ("all32_bit_positions", 0u32, 1048576u32, 3, false, 3u32),
    ("all32_bit_positions", 1048576u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        1048576u32,
        1048576u32,
        3,
        false,
        33554435u32,
    ),
    ("all32_bit_positions", 0u32, 2097152u32, 3, false, 3u32),
    ("all32_bit_positions", 2097152u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        2097152u32,
        2097152u32,
        3,
        false,
        67108867u32,
    ),
    ("all32_bit_positions", 0u32, 4194304u32, 3, false, 3u32),
    ("all32_bit_positions", 4194304u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        4194304u32,
        4194304u32,
        3,
        false,
        134217731u32,
    ),
    ("all32_bit_positions", 0u32, 8388608u32, 3, false, 3u32),
    ("all32_bit_positions", 8388608u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        8388608u32,
        8388608u32,
        3,
        false,
        268435459u32,
    ),
    ("all32_bit_positions", 0u32, 16777216u32, 3, false, 3u32),
    ("all32_bit_positions", 16777216u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        16777216u32,
        16777216u32,
        3,
        false,
        536870915u32,
    ),
    ("all32_bit_positions", 0u32, 33554432u32, 3, false, 3u32),
    ("all32_bit_positions", 33554432u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        33554432u32,
        33554432u32,
        3,
        false,
        1073741827u32,
    ),
    ("all32_bit_positions", 0u32, 67108864u32, 3, false, 3u32),
    ("all32_bit_positions", 67108864u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        67108864u32,
        67108864u32,
        3,
        false,
        2147483651u32,
    ),
    ("all32_bit_positions", 0u32, 134217728u32, 3, false, 3u32),
    ("all32_bit_positions", 134217728u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        134217728u32,
        134217728u32,
        3,
        false,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 268435456u32, 3, false, 3u32),
    ("all32_bit_positions", 268435456u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        268435456u32,
        268435456u32,
        3,
        false,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 536870912u32, 3, false, 3u32),
    ("all32_bit_positions", 536870912u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        536870912u32,
        536870912u32,
        3,
        false,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 1073741824u32, 3, false, 3u32),
    ("all32_bit_positions", 1073741824u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        1073741824u32,
        1073741824u32,
        3,
        false,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 2147483648u32, 3, false, 3u32),
    ("all32_bit_positions", 2147483648u32, 0u32, 3, false, 3u32),
    (
        "all32_bit_positions",
        2147483648u32,
        2147483648u32,
        3,
        false,
        3u32,
    ),
    (
        "chirality_and_highbit_overlap",
        511u32,
        512u32,
        29,
        false,
        8404989u32,
    ),
    (
        "chirality_and_highbit_overlap",
        512u32,
        1024u32,
        29,
        false,
        16793629u32,
    ),
    (
        "chirality_and_highbit_overlap",
        1023u32,
        1024u32,
        29,
        false,
        16809981u32,
    ),
    (
        "chirality_and_highbit_overlap",
        2047u32,
        2048u32,
        29,
        false,
        33619965u32,
    ),
    (
        "chirality_and_highbit_overlap",
        65535u32,
        65536u32,
        29,
        false,
        1075838973u32,
    ),
    (
        "chirality_and_highbit_overlap",
        131071u32,
        131072u32,
        29,
        false,
        2151677949u32,
    ),
    (
        "chirality_and_highbit_overlap",
        2147483647u32,
        4294967295u32,
        29,
        false,
        4294967293u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 0, true, 0u32),
    ("distance_0_through_30", 1u32, 2u32, 0, true, 131104u32),
    ("distance_0_through_30", 33u32, 34u32, 0, true, 2229280u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        0,
        true,
        4294901760u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 1, true, 1u32),
    ("distance_0_through_30", 1u32, 2u32, 1, true, 131105u32),
    ("distance_0_through_30", 33u32, 34u32, 1, true, 2229281u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        1,
        true,
        4294901761u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 2, true, 2u32),
    ("distance_0_through_30", 1u32, 2u32, 2, true, 131106u32),
    ("distance_0_through_30", 33u32, 34u32, 2, true, 2229282u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        2,
        true,
        4294901762u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 3, true, 3u32),
    ("distance_0_through_30", 1u32, 2u32, 3, true, 131107u32),
    ("distance_0_through_30", 33u32, 34u32, 3, true, 2229283u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        3,
        true,
        4294901763u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 4, true, 4u32),
    ("distance_0_through_30", 1u32, 2u32, 4, true, 131108u32),
    ("distance_0_through_30", 33u32, 34u32, 4, true, 2229284u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        4,
        true,
        4294901764u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 5, true, 5u32),
    ("distance_0_through_30", 1u32, 2u32, 5, true, 131109u32),
    ("distance_0_through_30", 33u32, 34u32, 5, true, 2229285u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        5,
        true,
        4294901765u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 6, true, 6u32),
    ("distance_0_through_30", 1u32, 2u32, 6, true, 131110u32),
    ("distance_0_through_30", 33u32, 34u32, 6, true, 2229286u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        6,
        true,
        4294901766u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 7, true, 7u32),
    ("distance_0_through_30", 1u32, 2u32, 7, true, 131111u32),
    ("distance_0_through_30", 33u32, 34u32, 7, true, 2229287u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        7,
        true,
        4294901767u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 8, true, 8u32),
    ("distance_0_through_30", 1u32, 2u32, 8, true, 131112u32),
    ("distance_0_through_30", 33u32, 34u32, 8, true, 2229288u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        8,
        true,
        4294901768u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 9, true, 9u32),
    ("distance_0_through_30", 1u32, 2u32, 9, true, 131113u32),
    ("distance_0_through_30", 33u32, 34u32, 9, true, 2229289u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        9,
        true,
        4294901769u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 10, true, 10u32),
    ("distance_0_through_30", 1u32, 2u32, 10, true, 131114u32),
    ("distance_0_through_30", 33u32, 34u32, 10, true, 2229290u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        10,
        true,
        4294901770u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 11, true, 11u32),
    ("distance_0_through_30", 1u32, 2u32, 11, true, 131115u32),
    ("distance_0_through_30", 33u32, 34u32, 11, true, 2229291u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        11,
        true,
        4294901771u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 12, true, 12u32),
    ("distance_0_through_30", 1u32, 2u32, 12, true, 131116u32),
    ("distance_0_through_30", 33u32, 34u32, 12, true, 2229292u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        12,
        true,
        4294901772u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 13, true, 13u32),
    ("distance_0_through_30", 1u32, 2u32, 13, true, 131117u32),
    ("distance_0_through_30", 33u32, 34u32, 13, true, 2229293u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        13,
        true,
        4294901773u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 14, true, 14u32),
    ("distance_0_through_30", 1u32, 2u32, 14, true, 131118u32),
    ("distance_0_through_30", 33u32, 34u32, 14, true, 2229294u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        14,
        true,
        4294901774u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 15, true, 15u32),
    ("distance_0_through_30", 1u32, 2u32, 15, true, 131119u32),
    ("distance_0_through_30", 33u32, 34u32, 15, true, 2229295u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        15,
        true,
        4294901775u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 16, true, 16u32),
    ("distance_0_through_30", 1u32, 2u32, 16, true, 131120u32),
    ("distance_0_through_30", 33u32, 34u32, 16, true, 2229296u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        16,
        true,
        4294901776u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 17, true, 17u32),
    ("distance_0_through_30", 1u32, 2u32, 17, true, 131121u32),
    ("distance_0_through_30", 33u32, 34u32, 17, true, 2229297u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        17,
        true,
        4294901777u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 18, true, 18u32),
    ("distance_0_through_30", 1u32, 2u32, 18, true, 131122u32),
    ("distance_0_through_30", 33u32, 34u32, 18, true, 2229298u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        18,
        true,
        4294901778u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 19, true, 19u32),
    ("distance_0_through_30", 1u32, 2u32, 19, true, 131123u32),
    ("distance_0_through_30", 33u32, 34u32, 19, true, 2229299u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        19,
        true,
        4294901779u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 20, true, 20u32),
    ("distance_0_through_30", 1u32, 2u32, 20, true, 131124u32),
    ("distance_0_through_30", 33u32, 34u32, 20, true, 2229300u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        20,
        true,
        4294901780u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 21, true, 21u32),
    ("distance_0_through_30", 1u32, 2u32, 21, true, 131125u32),
    ("distance_0_through_30", 33u32, 34u32, 21, true, 2229301u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        21,
        true,
        4294901781u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 22, true, 22u32),
    ("distance_0_through_30", 1u32, 2u32, 22, true, 131126u32),
    ("distance_0_through_30", 33u32, 34u32, 22, true, 2229302u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        22,
        true,
        4294901782u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 23, true, 23u32),
    ("distance_0_through_30", 1u32, 2u32, 23, true, 131127u32),
    ("distance_0_through_30", 33u32, 34u32, 23, true, 2229303u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        23,
        true,
        4294901783u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 24, true, 24u32),
    ("distance_0_through_30", 1u32, 2u32, 24, true, 131128u32),
    ("distance_0_through_30", 33u32, 34u32, 24, true, 2229304u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        24,
        true,
        4294901784u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 25, true, 25u32),
    ("distance_0_through_30", 1u32, 2u32, 25, true, 131129u32),
    ("distance_0_through_30", 33u32, 34u32, 25, true, 2229305u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        25,
        true,
        4294901785u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 26, true, 26u32),
    ("distance_0_through_30", 1u32, 2u32, 26, true, 131130u32),
    ("distance_0_through_30", 33u32, 34u32, 26, true, 2229306u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        26,
        true,
        4294901786u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 27, true, 27u32),
    ("distance_0_through_30", 1u32, 2u32, 27, true, 131131u32),
    ("distance_0_through_30", 33u32, 34u32, 27, true, 2229307u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        27,
        true,
        4294901787u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 28, true, 28u32),
    ("distance_0_through_30", 1u32, 2u32, 28, true, 131132u32),
    ("distance_0_through_30", 33u32, 34u32, 28, true, 2229308u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        28,
        true,
        4294901788u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 29, true, 29u32),
    ("distance_0_through_30", 1u32, 2u32, 29, true, 131133u32),
    ("distance_0_through_30", 33u32, 34u32, 29, true, 2229309u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        29,
        true,
        4294901789u32,
    ),
    ("distance_0_through_30", 0u32, 0u32, 30, true, 30u32),
    ("distance_0_through_30", 1u32, 2u32, 30, true, 131134u32),
    ("distance_0_through_30", 33u32, 34u32, 30, true, 2229310u32),
    (
        "distance_0_through_30",
        0u32,
        4294967295u32,
        30,
        true,
        4294901790u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        1u32,
        2u32,
        7,
        true,
        131111u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2u32,
        1u32,
        7,
        true,
        131111u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        33u32,
        34u32,
        7,
        true,
        2229287u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        34u32,
        33u32,
        7,
        true,
        2229287u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        4294967295u32,
        1u32,
        7,
        true,
        4294901799u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        1u32,
        4294967295u32,
        7,
        true,
        4294901799u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2147483648u32,
        2147483647u32,
        7,
        true,
        4294967271u32,
    ),
    (
        "minmax_unsigned_before_truncate",
        2147483647u32,
        2147483648u32,
        7,
        true,
        4294967271u32,
    ),
    ("equal_codes_full_width", 0u32, 0u32, 30, true, 30u32),
    ("equal_codes_full_width", 1u32, 1u32, 30, true, 65598u32),
    (
        "equal_codes_full_width",
        511u32,
        511u32,
        30,
        true,
        33505278u32,
    ),
    (
        "equal_codes_full_width",
        512u32,
        512u32,
        30,
        true,
        33570846u32,
    ),
    (
        "equal_codes_full_width",
        1023u32,
        1023u32,
        30,
        true,
        67076094u32,
    ),
    (
        "equal_codes_full_width",
        1024u32,
        1024u32,
        30,
        true,
        67141662u32,
    ),
    (
        "equal_codes_full_width",
        2047u32,
        2047u32,
        30,
        true,
        134217726u32,
    ),
    (
        "equal_codes_full_width",
        2048u32,
        2048u32,
        30,
        true,
        134283294u32,
    ),
    (
        "equal_codes_full_width",
        4294967295u32,
        4294967295u32,
        30,
        true,
        4294967294u32,
    ),
    ("all32_bit_positions", 0u32, 1u32, 3, true, 65539u32),
    ("all32_bit_positions", 1u32, 0u32, 3, true, 65539u32),
    ("all32_bit_positions", 1u32, 1u32, 3, true, 65571u32),
    ("all32_bit_positions", 0u32, 2u32, 3, true, 131075u32),
    ("all32_bit_positions", 2u32, 0u32, 3, true, 131075u32),
    ("all32_bit_positions", 2u32, 2u32, 3, true, 131139u32),
    ("all32_bit_positions", 0u32, 4u32, 3, true, 262147u32),
    ("all32_bit_positions", 4u32, 0u32, 3, true, 262147u32),
    ("all32_bit_positions", 4u32, 4u32, 3, true, 262275u32),
    ("all32_bit_positions", 0u32, 8u32, 3, true, 524291u32),
    ("all32_bit_positions", 8u32, 0u32, 3, true, 524291u32),
    ("all32_bit_positions", 8u32, 8u32, 3, true, 524547u32),
    ("all32_bit_positions", 0u32, 16u32, 3, true, 1048579u32),
    ("all32_bit_positions", 16u32, 0u32, 3, true, 1048579u32),
    ("all32_bit_positions", 16u32, 16u32, 3, true, 1049091u32),
    ("all32_bit_positions", 0u32, 32u32, 3, true, 2097155u32),
    ("all32_bit_positions", 32u32, 0u32, 3, true, 2097155u32),
    ("all32_bit_positions", 32u32, 32u32, 3, true, 2098179u32),
    ("all32_bit_positions", 0u32, 64u32, 3, true, 4194307u32),
    ("all32_bit_positions", 64u32, 0u32, 3, true, 4194307u32),
    ("all32_bit_positions", 64u32, 64u32, 3, true, 4196355u32),
    ("all32_bit_positions", 0u32, 128u32, 3, true, 8388611u32),
    ("all32_bit_positions", 128u32, 0u32, 3, true, 8388611u32),
    ("all32_bit_positions", 128u32, 128u32, 3, true, 8392707u32),
    ("all32_bit_positions", 0u32, 256u32, 3, true, 16777219u32),
    ("all32_bit_positions", 256u32, 0u32, 3, true, 16777219u32),
    ("all32_bit_positions", 256u32, 256u32, 3, true, 16785411u32),
    ("all32_bit_positions", 0u32, 512u32, 3, true, 33554435u32),
    ("all32_bit_positions", 512u32, 0u32, 3, true, 33554435u32),
    ("all32_bit_positions", 512u32, 512u32, 3, true, 33570819u32),
    ("all32_bit_positions", 0u32, 1024u32, 3, true, 67108867u32),
    ("all32_bit_positions", 1024u32, 0u32, 3, true, 67108867u32),
    (
        "all32_bit_positions",
        1024u32,
        1024u32,
        3,
        true,
        67141635u32,
    ),
    ("all32_bit_positions", 0u32, 2048u32, 3, true, 134217731u32),
    ("all32_bit_positions", 2048u32, 0u32, 3, true, 134217731u32),
    (
        "all32_bit_positions",
        2048u32,
        2048u32,
        3,
        true,
        134283267u32,
    ),
    ("all32_bit_positions", 0u32, 4096u32, 3, true, 268435459u32),
    ("all32_bit_positions", 4096u32, 0u32, 3, true, 268435459u32),
    (
        "all32_bit_positions",
        4096u32,
        4096u32,
        3,
        true,
        268566531u32,
    ),
    ("all32_bit_positions", 0u32, 8192u32, 3, true, 536870915u32),
    ("all32_bit_positions", 8192u32, 0u32, 3, true, 536870915u32),
    (
        "all32_bit_positions",
        8192u32,
        8192u32,
        3,
        true,
        537133059u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        16384u32,
        3,
        true,
        1073741827u32,
    ),
    (
        "all32_bit_positions",
        16384u32,
        0u32,
        3,
        true,
        1073741827u32,
    ),
    (
        "all32_bit_positions",
        16384u32,
        16384u32,
        3,
        true,
        1074266115u32,
    ),
    (
        "all32_bit_positions",
        0u32,
        32768u32,
        3,
        true,
        2147483651u32,
    ),
    (
        "all32_bit_positions",
        32768u32,
        0u32,
        3,
        true,
        2147483651u32,
    ),
    (
        "all32_bit_positions",
        32768u32,
        32768u32,
        3,
        true,
        2148532227u32,
    ),
    ("all32_bit_positions", 0u32, 65536u32, 3, true, 3u32),
    ("all32_bit_positions", 65536u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        65536u32,
        65536u32,
        3,
        true,
        2097155u32,
    ),
    ("all32_bit_positions", 0u32, 131072u32, 3, true, 3u32),
    ("all32_bit_positions", 131072u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        131072u32,
        131072u32,
        3,
        true,
        4194307u32,
    ),
    ("all32_bit_positions", 0u32, 262144u32, 3, true, 3u32),
    ("all32_bit_positions", 262144u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        262144u32,
        262144u32,
        3,
        true,
        8388611u32,
    ),
    ("all32_bit_positions", 0u32, 524288u32, 3, true, 3u32),
    ("all32_bit_positions", 524288u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        524288u32,
        524288u32,
        3,
        true,
        16777219u32,
    ),
    ("all32_bit_positions", 0u32, 1048576u32, 3, true, 3u32),
    ("all32_bit_positions", 1048576u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        1048576u32,
        1048576u32,
        3,
        true,
        33554435u32,
    ),
    ("all32_bit_positions", 0u32, 2097152u32, 3, true, 3u32),
    ("all32_bit_positions", 2097152u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        2097152u32,
        2097152u32,
        3,
        true,
        67108867u32,
    ),
    ("all32_bit_positions", 0u32, 4194304u32, 3, true, 3u32),
    ("all32_bit_positions", 4194304u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        4194304u32,
        4194304u32,
        3,
        true,
        134217731u32,
    ),
    ("all32_bit_positions", 0u32, 8388608u32, 3, true, 3u32),
    ("all32_bit_positions", 8388608u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        8388608u32,
        8388608u32,
        3,
        true,
        268435459u32,
    ),
    ("all32_bit_positions", 0u32, 16777216u32, 3, true, 3u32),
    ("all32_bit_positions", 16777216u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        16777216u32,
        16777216u32,
        3,
        true,
        536870915u32,
    ),
    ("all32_bit_positions", 0u32, 33554432u32, 3, true, 3u32),
    ("all32_bit_positions", 33554432u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        33554432u32,
        33554432u32,
        3,
        true,
        1073741827u32,
    ),
    ("all32_bit_positions", 0u32, 67108864u32, 3, true, 3u32),
    ("all32_bit_positions", 67108864u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        67108864u32,
        67108864u32,
        3,
        true,
        2147483651u32,
    ),
    ("all32_bit_positions", 0u32, 134217728u32, 3, true, 3u32),
    ("all32_bit_positions", 134217728u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        134217728u32,
        134217728u32,
        3,
        true,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 268435456u32, 3, true, 3u32),
    ("all32_bit_positions", 268435456u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        268435456u32,
        268435456u32,
        3,
        true,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 536870912u32, 3, true, 3u32),
    ("all32_bit_positions", 536870912u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        536870912u32,
        536870912u32,
        3,
        true,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 1073741824u32, 3, true, 3u32),
    ("all32_bit_positions", 1073741824u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        1073741824u32,
        1073741824u32,
        3,
        true,
        3u32,
    ),
    ("all32_bit_positions", 0u32, 2147483648u32, 3, true, 3u32),
    ("all32_bit_positions", 2147483648u32, 0u32, 3, true, 3u32),
    (
        "all32_bit_positions",
        2147483648u32,
        2147483648u32,
        3,
        true,
        3u32,
    ),
    (
        "chirality_and_highbit_overlap",
        511u32,
        512u32,
        29,
        true,
        33570813u32,
    ),
    (
        "chirality_and_highbit_overlap",
        512u32,
        1024u32,
        29,
        true,
        67125277u32,
    ),
    (
        "chirality_and_highbit_overlap",
        1023u32,
        1024u32,
        29,
        true,
        67141629u32,
    ),
    (
        "chirality_and_highbit_overlap",
        2047u32,
        2048u32,
        29,
        true,
        134283261u32,
    ),
    (
        "chirality_and_highbit_overlap",
        65535u32,
        65536u32,
        29,
        true,
        2097149u32,
    ),
    (
        "chirality_and_highbit_overlap",
        131071u32,
        131072u32,
        29,
        true,
        4194301u32,
    ),
    (
        "chirality_and_highbit_overlap",
        2147483647u32,
        4294967295u32,
        29,
        true,
        4294967293u32,
    ),
];
const CODE: &[(&str, &[u32], bool, u64)] = &[
    ("all_defined_lengths", &[1], false, 1u64),
    ("all_defined_lengths", &[0], false, 0u64),
    ("all_defined_lengths", &[1, 2], false, 1025u64),
    ("all_defined_lengths", &[0, 0], false, 0u64),
    ("all_defined_lengths", &[1, 2, 3], false, 787457u64),
    ("all_defined_lengths", &[0, 0, 0], false, 0u64),
    ("all_defined_lengths", &[1, 2, 3, 4], false, 537658369u64),
    ("all_defined_lengths", &[0, 0, 0, 0], false, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5],
        false,
        344135042049u64,
    ),
    ("all_defined_lengths", &[0, 0, 0, 0, 0], false, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5, 6],
        false,
        211450367575041u64,
    ),
    ("all_defined_lengths", &[0, 0, 0, 0, 0, 0], false, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5, 6, 7],
        false,
        126312239933948929u64,
    ),
    ("all_defined_lengths", &[0, 0, 0, 0, 0, 0, 0], false, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5, 6, 7, 8],
        false,
        126312239933948929u64,
    ),
    (
        "all_defined_lengths",
        &[0, 0, 0, 0, 0, 0, 0, 0],
        false,
        0u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 3, 4],
        false,
        537658369u64,
    ),
    (
        "paired_end_order_branches",
        &[4, 3, 2, 1],
        false,
        537658369u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 9, 2, 1],
        false,
        136578049u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 9, 1],
        false,
        136578049u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 7, 8, 2, 1],
        false,
        35322886620161u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 8, 7, 2, 1],
        false,
        35322886620161u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 3, 2, 4],
        false,
        537396737u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 3, 4],
        false,
        537658369u64,
    ),
    ("palindrome_equal_pairs", &[5, 3, 5], false, 1312261u64),
    ("palindrome_equal_pairs", &[3, 9, 9, 3], false, 405017091u64),
    ("palindrome_equal_pairs", &[0, 0, 0, 0], false, 0u64),
    (
        "palindrome_equal_pairs",
        &[4294967295, 1, 1, 4294967295],
        false,
        576460752303423487u64,
    ),
    ("palindrome_equal_pairs", &[7, 7], false, 3591u64),
    ("one_element_full_u32", &[0], false, 0u64),
    ("one_element_full_u32", &[1], false, 1u64),
    ("one_element_full_u32", &[511], false, 511u64),
    ("one_element_full_u32", &[512], false, 512u64),
    ("one_element_full_u32", &[2047], false, 2047u64),
    ("one_element_full_u32", &[2048], false, 2048u64),
    ("one_element_full_u32", &[2147483648], false, 2147483648u64),
    ("one_element_full_u32", &[4294967295], false, 4294967295u64),
    ("OR_overlap_and_full_width", &[512, 513], false, 262656u64),
    (
        "OR_overlap_and_full_width",
        &[2048, 2049],
        false,
        1051136u64,
    ),
    (
        "OR_overlap_and_full_width",
        &[4294967294, 4294967295],
        false,
        2199023255550u64,
    ),
    (
        "OR_overlap_and_full_width",
        &[0, 4294967295, 1],
        false,
        2199023255040u64,
    ),
    (
        "OR_overlap_and_full_width",
        &[1, 4294967295, 2, 3],
        false,
        2199023255041u64,
    ),
    ("last_defined_shift", &[0, 0, 0, 0, 0, 0, 0, 0], false, 0u64),
    (
        "last_defined_shift",
        &[0, 0, 0, 0, 0, 0, 0, 1],
        false,
        9223372036854775808u64,
    ),
    ("last_defined_shift", &[0, 0, 0, 0, 0, 0, 0, 2], false, 0u64),
    (
        "last_defined_shift",
        &[0, 0, 0, 0, 0, 0, 0, 4294967295],
        false,
        9223372036854775808u64,
    ),
    ("all_defined_lengths", &[1], true, 1u64),
    ("all_defined_lengths", &[0], true, 0u64),
    ("all_defined_lengths", &[1, 2], true, 4097u64),
    ("all_defined_lengths", &[0, 0], true, 0u64),
    ("all_defined_lengths", &[1, 2, 3], true, 12587009u64),
    ("all_defined_lengths", &[0, 0, 0], true, 0u64),
    ("all_defined_lengths", &[1, 2, 3, 4], true, 34372325377u64),
    ("all_defined_lengths", &[0, 0, 0, 0], true, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5],
        true,
        87995302547457u64,
    ),
    ("all_defined_lengths", &[0, 0, 0, 0, 0], true, 0u64),
    (
        "all_defined_lengths",
        &[1, 2, 3, 4, 5, 6],
        true,
        216260777416331265u64,
    ),
    ("all_defined_lengths", &[0, 0, 0, 0, 0, 0], true, 0u64),
    (
        "paired_end_order_branches",
        &[1, 2, 3, 4],
        true,
        34372325377u64,
    ),
    (
        "paired_end_order_branches",
        &[4, 3, 2, 1],
        true,
        34372325377u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 9, 2, 1],
        true,
        8627687425u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 9, 1],
        true,
        8627687425u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 7, 8, 2, 1],
        true,
        36064050139893761u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 8, 7, 2, 1],
        true,
        36064050139893761u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 3, 2, 4],
        true,
        34368133121u64,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 3, 4],
        true,
        34372325377u64,
    ),
    ("palindrome_equal_pairs", &[5, 3, 5], true, 20977669u64),
    (
        "palindrome_equal_pairs",
        &[3, 9, 9, 3],
        true,
        25807570947u64,
    ),
    ("palindrome_equal_pairs", &[0, 0, 0, 0], true, 0u64),
    (
        "palindrome_equal_pairs",
        &[4294967295, 1, 1, 4294967295],
        true,
        18446744069414584319u64,
    ),
    ("palindrome_equal_pairs", &[7, 7], true, 14343u64),
    ("one_element_full_u32", &[0], true, 0u64),
    ("one_element_full_u32", &[1], true, 1u64),
    ("one_element_full_u32", &[511], true, 511u64),
    ("one_element_full_u32", &[512], true, 512u64),
    ("one_element_full_u32", &[2047], true, 2047u64),
    ("one_element_full_u32", &[2048], true, 2048u64),
    ("one_element_full_u32", &[2147483648], true, 2147483648u64),
    ("one_element_full_u32", &[4294967295], true, 4294967295u64),
    ("OR_overlap_and_full_width", &[512, 513], true, 1051136u64),
    ("OR_overlap_and_full_width", &[2048, 2049], true, 4196352u64),
    (
        "OR_overlap_and_full_width",
        &[4294967294, 4294967295],
        true,
        8796093022206u64,
    ),
    (
        "OR_overlap_and_full_width",
        &[0, 4294967295, 1],
        true,
        8796093020160u64,
    ),
    (
        "OR_overlap_and_full_width",
        &[1, 4294967295, 2, 3],
        true,
        8796093020161u64,
    ),
    ("last_defined_shift", &[0, 0, 0, 0, 0, 0], true, 0u64),
    (
        "last_defined_shift",
        &[0, 0, 0, 0, 0, 1],
        true,
        36028797018963968u64,
    ),
    (
        "last_defined_shift",
        &[0, 0, 0, 0, 0, 2],
        true,
        72057594037927936u64,
    ),
    (
        "last_defined_shift",
        &[0, 0, 0, 0, 0, 4294967295],
        true,
        18410715276690587648u64,
    ),
];
const HASH: &[(&str, &[u32], u32)] = &[
    ("one_element_u32_carry", &[0], 2654435769u32),
    ("one_element_u32_carry", &[1], 2654435770u32),
    ("one_element_u32_carry", &[2], 2654435771u32),
    ("one_element_u32_carry", &[2147483647], 506952120u32),
    ("one_element_u32_carry", &[2147483648], 506952121u32),
    ("one_element_u32_carry", &[4294967295], 2654435768u32),
    ("paired_end_order_branches", &[1, 2], 3449077523u32),
    ("paired_end_order_branches", &[2, 1], 3449077523u32),
    ("paired_end_order_branches", &[1, 2, 3, 4], 1209717634u32),
    ("paired_end_order_branches", &[4, 3, 2, 1], 1209717634u32),
    ("paired_end_order_branches", &[1, 9, 2, 1], 1209717274u32),
    ("paired_end_order_branches", &[1, 2, 9, 1], 1209717274u32),
    (
        "paired_end_order_branches",
        &[1, 2, 7, 8, 2, 1],
        1762696026u32,
    ),
    (
        "paired_end_order_branches",
        &[1, 2, 8, 7, 2, 1],
        1762696026u32,
    ),
    ("paired_end_order_branches", &[1, 3, 2, 4], 1209702519u32),
    ("paired_end_order_branches", &[1, 2, 3, 4], 1209717634u32),
    ("palindrome_equal_pairs", &[5, 3, 5], 4216885398u32),
    ("palindrome_equal_pairs", &[3, 9, 9, 3], 1209184878u32),
    ("palindrome_equal_pairs", &[7, 7], 3449074160u32),
    (
        "palindrome_equal_pairs",
        &[4294967295, 1, 1, 4294967295],
        1212300213u32,
    ),
    ("palindrome_equal_pairs", &[0, 0, 0, 0], 1213014554u32),
    ("palindrome_equal_pairs", &[1, 0, 1], 4216901596u32),
    ("all_ordinary_lengths_no_packing_cap", &[0], 2654435769u32),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[4294967295],
        2654435768u32,
    ),
    ("all_ordinary_lengths_no_packing_cap", &[0], 2654435769u32),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1],
        3449077713u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[4294967295, 4294967295],
        3449077662u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0],
        3449077726u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2],
        4216857150u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[4294967295, 4294967295, 4294967295],
        4216860289u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0],
        4216856302u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3],
        1213084661u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[4294967295, 4294967295, 4294967295, 4294967295],
        1212222233u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0],
        1213014554u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4],
        2341875727u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[4294967295, 4294967295, 4294967295, 4294967295, 4294967295],
        2295040423u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0],
        2346609317u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5],
        759118222u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
        ],
        2072464454u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0],
        857200391u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5, 6],
        3563754028u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
        ],
        3848999439u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0, 0],
        1139074621u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5, 6, 7],
        966540135u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295,
        ],
        3611153396u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0, 0, 0],
        3951994037u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5, 6, 7, 8],
        707877181u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295,
        ],
        1950758977u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0, 0, 0, 0],
        1464530067u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
        522929772u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295,
        ],
        3096518729u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        3515723534u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
        1627643959u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295,
        ],
        1593840943u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        3565488559u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30,
        ],
        3049133390u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295,
        ],
        3618750256u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0,
        ],
        4051500288u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23,
            24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
            46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
        ],
        1601204595u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295, 4294967295,
            4294967295,
        ],
        3596870173u32,
    ),
    (
        "all_ordinary_lengths_no_packing_cap",
        &[
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0,
        ],
        4238431781u32,
    ),
    ("all32_bit_positions", &[0, 1, 4294967295], 4216857149u32),
    ("all32_bit_positions", &[0, 2, 4294967295], 4216857212u32),
    ("all32_bit_positions", &[0, 4, 4294967295], 4216857086u32),
    ("all32_bit_positions", &[0, 8, 4294967295], 4216856827u32),
    ("all32_bit_positions", &[0, 16, 4294967295], 4216857317u32),
    ("all32_bit_positions", &[0, 32, 4294967295], 4216899129u32),
    ("all32_bit_positions", &[0, 64, 4294967295], 4216901601u32),
    ("all32_bit_positions", &[0, 128, 4294967295], 4216864337u32),
    ("all32_bit_positions", &[0, 256, 4294967295], 4216905265u32),
    ("all32_bit_positions", &[0, 512, 4294967295], 4217220209u32),
    ("all32_bit_positions", &[0, 1024, 4294967295], 4217187825u32),
    ("all32_bit_positions", &[0, 2048, 4294967295], 4217252081u32),
    ("all32_bit_positions", &[0, 4096, 4294967295], 4217123569u32),
    ("all32_bit_positions", &[0, 8192, 4294967295], 4216325873u32),
    (
        "all32_bit_positions",
        &[0, 16384, 4294967295],
        4226756337u32,
    ),
    (
        "all32_bit_positions",
        &[0, 32768, 4294967295],
        4227825393u32,
    ),
    (
        "all32_bit_positions",
        &[0, 65536, 4294967295],
        4213169905u32,
    ),
    (
        "all32_bit_positions",
        &[0, 131072, 4294967295],
        4225670897u32,
    ),
    (
        "all32_bit_positions",
        &[0, 262144, 4294967295],
        4166721265u32,
    ),
    (
        "all32_bit_positions",
        &[0, 524288, 4294967295],
        4115799793u32,
    ),
    (
        "all32_bit_positions",
        &[0, 1048576, 4294967295],
        4283178737u32,
    ),
    (
        "all32_bit_positions",
        &[0, 2097152, 4294967295],
        2198871793u32,
    ),
    (
        "all32_bit_positions",
        &[0, 4194304, 4294967295],
        2332565233u32,
    ),
    (
        "all32_bit_positions",
        &[0, 8388608, 4294967295],
        2683838193u32,
    ),
    (
        "all32_bit_positions",
        &[0, 16777216, 4294967295],
        3164086001u32,
    ),
    (
        "all32_bit_positions",
        &[0, 33554432, 4294967295],
        2111315697u32,
    ),
    (
        "all32_bit_positions",
        &[0, 67108864, 4294967295],
        4233633521u32,
    ),
    (
        "all32_bit_positions",
        &[0, 134217728, 4294967295],
        4049084145u32,
    ),
    (
        "all32_bit_positions",
        &[0, 268435456, 4294967295],
        3210223345u32,
    ),
    (
        "all32_bit_positions",
        &[0, 536870912, 4294967295],
        3545767665u32,
    ),
    (
        "all32_bit_positions",
        &[0, 1073741824, 4294967295],
        190324465u32,
    ),
    (
        "all32_bit_positions",
        &[0, 2147483648, 4294967295],
        1532501745u32,
    ),
    (
        "full_width_overflow",
        &[4294967295, 2147483648, 1, 2147483647],
        1044527990u32,
    ),
    (
        "full_width_overflow",
        &[2147483647, 1, 2147483648, 4294967295],
        1044527990u32,
    ),
    ("full_width_overflow", &[0, 4294967295, 0], 4216856239u32),
    (
        "full_width_overflow",
        &[4294967295, 4294967294],
        3449076818u32,
    ),
    (
        "full_width_overflow",
        &[4294967294, 4294967295],
        3449076818u32,
    ),
];
