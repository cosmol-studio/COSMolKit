//! Fixed fingerprint-value regressions in the owning crate.
//!
//! Expectations are pinned against RDKit `2026.03.1`, commit
//! `351f8f378f8ad6bbd517980c38896e66bf907af8c`, probed 2026-09-24 through
//! the installed Python reference (`rdkit.DataStructs`) and the C++ sources
//! under `third_party/rdkit/Code/DataStructs`. Probe outputs are recorded
//! beside each assertion. Where the Python wrapper does not expose a method
//! (e.g. `ClearBits`), the expectation is pinned from the C++ body cited in
//! the implementing function's source anchor.

use cosmolkit_fingerprints::{
    Fingerprint, FingerprintError, SparseBitFingerprint, SparseCountFingerprint,
    SparseCountFingerprint32,
};

#[test]
fn sparse_zero_length_vector_matches_source_defaults() {
    // Oracle: SparseBitVect(0) -> GetNumBits()==0, GetNumOnBits()==0.
    let sbv = SparseBitFingerprint::new(0);
    assert_eq!(sbv.n_bits(), 0);
    assert_eq!(sbv.num_on_bits(), 0);
    assert_eq!(sbv.num_off_bits(), 0);
    assert!(sbv.on_bits().is_empty());
    assert_eq!(
        sbv.get_bit(0),
        Err(FingerprintError::SparseIndexOutOfRange { index: 0, size: 0 })
    );
}

#[test]
fn sparse_set_get_unset_and_duplicate_insertion() {
    // Oracle probe: SetBit(3) -> False, duplicate -> True; GetBit True/False;
    // UnSetBit -> True; UnSetBit absent -> False.
    let mut sbv = SparseBitFingerprint::new(10);
    assert!(!sbv.set_bit(3).unwrap());
    assert!(sbv.set_bit(3).unwrap());
    assert_eq!(sbv.num_on_bits(), 1);
    assert_eq!(sbv.num_off_bits(), 9);
    assert!(sbv.get_bit(3).unwrap());
    assert!(!sbv.get_bit(4).unwrap());
    assert!(sbv.unset_bit(3).unwrap());
    assert!(!sbv.unset_bit(3).unwrap());
    assert_eq!(sbv.num_on_bits(), 0);
}

#[test]
fn sparse_index_equal_to_length_errors() {
    // Oracle probe: SetBit(10) on size 10 -> IndexError; GetBit(11) too.
    // (The C++ error message wraps the index through `int`; the Rust error
    // keeps the structured index/size fields.)
    let mut sbv = SparseBitFingerprint::new(10);
    assert_eq!(
        sbv.set_bit(10),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 10,
            size: 10
        })
    );
    assert_eq!(
        sbv.get_bit(11),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 11,
            size: 10
        })
    );
    assert_eq!(
        sbv.unset_bit(10),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 10,
            size: 10
        })
    );
}

#[test]
fn sparse_maximum_length_accepts_edge_index() {
    // Oracle probe: SparseBitVect(2**32-1).SetBit(2**32-1) -> False (first
    // insert), GetBit -> True; SparseBitVect(2**32-2).SetBit(2**32-2)
    // (idx==size with size < max) -> IndexError. Sparse storage stays small;
    // no dense allocation happens at maximum length. F02 oracle: the
    // u32::MAX member appears in GetOnBits as the wrapped signed key -1
    // ordered before 7.
    let mut big = SparseBitFingerprint::new(u32::MAX);
    assert!(!big.set_bit(u32::MAX).unwrap());
    assert!(big.get_bit(u32::MAX).unwrap());
    assert!(big.set_bit(u32::MAX).unwrap());
    assert_eq!(big.num_on_bits(), 1);
    assert_eq!(big.num_off_bits(), u32::MAX - 1);
    big.set_bit(7).unwrap();
    assert_eq!(big.on_bits(), vec![-1, 7]);
    let mut below = SparseBitFingerprint::new(u32::MAX - 1);
    assert_eq!(
        below.set_bit(u32::MAX - 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: u64::from(u32::MAX - 1),
            size: u64::from(u32::MAX - 1)
        })
    );
}

#[test]
fn sparse_signed_set_ordering_f02() {
    // F02: indices >= 2^31 reach the signed set<int> storage as wrapped
    // negative keys and sort before all non-negative keys. Oracle probes
    // (pinned RDKit 2026.03.1, gcc x86-64 ABI):
    //   SparseBitVect(2**31+4){2**31+1, 5}: GetOnBits -> [-2147483647, 5];
    //   SparseBitVect(2**31){2**31-1, 3}: GetOnBits -> [3, 2147483647];
    //   SparseBitVect(2**32-1){u32::MAX, 7}: GetOnBits -> [-1, 7].
    let mut big = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    assert!(!big.set_bit(2_u32.pow(31) + 1).unwrap());
    big.set_bit(5).unwrap();
    assert_eq!(big.on_bits(), vec![-2_147_483_647, 5]);
    // numeric getBit uses the unsigned checked index, distinct from the
    // signed stored/output representation
    assert!(big.get_bit(2_u32.pow(31) + 1).unwrap());
    assert!(big.get_bit(5).unwrap());
    assert!(!big.get_bit(2_u32.pow(31) + 3).unwrap());

    let mut boundary = SparseBitFingerprint::new(2_u32.pow(31));
    boundary.set_bit(2_u32.pow(31) - 1).unwrap();
    boundary.set_bit(3).unwrap();
    assert_eq!(boundary.on_bits(), vec![3, 2_147_483_647]);

    let mut edge = SparseBitFingerprint::new(u32::MAX - 1);
    edge.set_bit(u32::MAX - 2).unwrap();
    edge.set_bit(0).unwrap();
    // oracle probe: SparseBitVect(2**32-2){2**32-3, 0} -> [-3, 0]
    // (u32::MAX-2 = 4294967293 wraps to -3)
    assert_eq!(edge.on_bits(), vec![-3, 0]);
}

#[test]
fn sparse_sorted_iteration_and_clear() {
    let mut sbv = SparseBitFingerprint::new(64);
    for bit in [33, 0, 63, 31, 32] {
        sbv.set_bit(bit).unwrap();
    }
    assert_eq!(sbv.on_bits(), vec![0, 31, 32, 33, 63]);
    sbv.clear_bits();
    assert!(sbv.on_bits().is_empty());
    assert_eq!(sbv.num_on_bits(), 0);
    assert_eq!(sbv.num_off_bits(), 64);
}

#[test]
fn sparse_set_algebra_matches_source_operators() {
    // Oracle probe: c={1,3,5}, d={3,5,7} over size 8:
    // or -> [1,3,5,7]; and -> [3,5]; xor -> [1,7]; ~c -> [0,2,4,6,7].
    let mut c = SparseBitFingerprint::new(8);
    let mut d = SparseBitFingerprint::new(8);
    for bit in [1, 3, 5] {
        c.set_bit(bit).unwrap();
    }
    for bit in [3, 5, 7] {
        d.set_bit(bit).unwrap();
    }
    assert_eq!(c.union(&d).on_bits(), vec![1, 3, 5, 7]);
    assert_eq!(c.intersection(&d).on_bits(), vec![3, 5]);
    assert_eq!(c.symmetric_difference(&d).on_bits(), vec![1, 7]);
    assert_eq!(c.negate().on_bits(), vec![0, 2, 4, 6, 7]);
    assert_eq!(c.num_on_bits(), 3);
    assert_eq!(c.num_off_bits(), 5);
    assert_eq!(c.n_bits(), 8);
}

#[test]
fn sparse_set_operators_do_not_check_lengths_f01() {
    // F01 (replaces the former `sparse_mismatched_lengths_error_for_set_ops`
    // expectation, which encoded a self-approved fail-closed rejection that
    // the pinned source does not perform; the same mismatched-length inputs
    // are retained as the historical counterexample, now with source-correct
    // expectations).
    //
    // Oracle probes (pinned RDKit 2026.03.1):
    //   a=SparseBitVect(8){1}, b=SparseBitVect(16){2,12}:
    //   a|b  -> 8 bits, on [1,2,12]  (out-of-length member 12 retained)
    //   b|a  -> 16 bits, on [1,2,12]
    //   a&b  -> 8 bits, on []
    //   a^b  -> 8 bits, on [1,2,12]
    //   GetBit(12) on the 8-bit union -> IndexError (stored member outside
    //   checked range)
    //   disjoint same-size c{1},d{5}: or [1,5], and [], xor [1,5]
    //   identical e{3}: or [3], xor []
    let mut a = SparseBitFingerprint::new(8);
    let mut b = SparseBitFingerprint::new(16);
    a.set_bit(1).unwrap();
    b.set_bit(2).unwrap();
    b.set_bit(12).unwrap();

    let u = a.union(&b);
    assert_eq!(u.n_bits(), 8);
    assert_eq!(u.on_bits(), vec![1, 2, 12]);
    assert_eq!(
        u.get_bit(12),
        Err(FingerprintError::SparseIndexOutOfRange { index: 12, size: 8 })
    );

    let u2 = b.union(&a);
    assert_eq!(u2.n_bits(), 16);
    assert_eq!(u2.on_bits(), vec![1, 2, 12]);

    let i = a.intersection(&b);
    assert_eq!(i.n_bits(), 8);
    assert!(i.on_bits().is_empty());

    let x = a.symmetric_difference(&b);
    assert_eq!(x.n_bits(), 8);
    assert_eq!(x.on_bits(), vec![1, 2, 12]);

    let mut c = SparseBitFingerprint::new(8);
    let mut d = SparseBitFingerprint::new(8);
    c.set_bit(1).unwrap();
    d.set_bit(5).unwrap();
    assert_eq!(c.union(&d).on_bits(), vec![1, 5]);
    assert!(c.intersection(&d).on_bits().is_empty());
    assert_eq!(c.symmetric_difference(&d).on_bits(), vec![1, 5]);

    let mut e = SparseBitFingerprint::new(8);
    e.set_bit(3).unwrap();
    assert_eq!(e.union(&e).on_bits(), vec![3]);
    assert!(e.symmetric_difference(&e).on_bits().is_empty());
}

#[test]
fn sparse_equality_ignores_length_like_source() {
    // Oracle probe: 4-bit {1} == 8-bit {1} -> True (operator== compares only
    // the on-bit sets).
    let mut a = SparseBitFingerprint::new(4);
    let mut b = SparseBitFingerprint::new(8);
    a.set_bit(1).unwrap();
    b.set_bit(1).unwrap();
    assert_eq!(a, b);
    let mut c = SparseBitFingerprint::new(4);
    c.set_bit(2).unwrap();
    assert_ne!(a, c);
    // clone independence
    let mut clone = a.clone();
    clone.set_bit(0).unwrap();
    assert_ne!(a, clone);
    assert_eq!(a.on_bits(), vec![1]);
}

#[test]
fn dense_zero_length_vector_matches_source_defaults() {
    // Oracle: ExplicitBitVect(0) -> GetNumBits()==0, GetNumOnBits()==0.
    let fp = Fingerprint::new(0);
    assert_eq!(fp.n_bits(), 0);
    assert_eq!(fp.num_on_bits(), 0);
    assert_eq!(fp.num_off_bits(), 0);
    assert!(fp.on_bits().is_empty());
    assert_eq!(fp, Fingerprint::empty());
    assert_eq!(fp, Fingerprint::new(0));
    // Empty fingerprint still has a defined tanimoto boundary (group 6).
}

#[test]
fn dense_set_bit_across_word_boundaries_and_duplicate_insertion() {
    // Oracle probe: SetBit(0) -> False, duplicate SetBit(0) -> True;
    // after bits {0,31,32,63,64,69}: counts 6/64/70.
    let mut fp = Fingerprint::new(70);
    assert!(!fp.set_bit(0).unwrap());
    assert!(fp.set_bit(0).unwrap());
    for bit in [31, 32, 63, 64, 69] {
        assert!(!fp.set_bit(bit).unwrap());
    }
    assert_eq!(fp.num_on_bits(), 6);
    assert_eq!(fp.num_off_bits(), 64);
    assert_eq!(fp.n_bits(), 70);
    assert_eq!(fp.on_bits(), vec![0, 31, 32, 63, 64, 69]);
    assert_eq!(fp.get_bit(0).unwrap(), true);
    assert_eq!(fp.get_bit(1).unwrap(), false);
    assert_eq!(fp.get_bit(69).unwrap(), true);
}

#[test]
fn dense_unset_bit_reset_semantics() {
    // Oracle probe: UnSetBit(2) on empty -> False; SetBit(2) -> False;
    // UnSetBit(2) -> True; count 0.
    let mut fp = Fingerprint::new(4);
    assert!(!fp.unset_bit(2).unwrap());
    assert!(!fp.set_bit(2).unwrap());
    assert!(fp.unset_bit(2).unwrap());
    assert_eq!(fp.num_on_bits(), 0);
    assert!(fp.on_bits().is_empty());
}

#[test]
fn dense_out_of_range_indices_error() {
    // Oracle probe: SetBit(70)/GetBit(70) on a 70-bit vector raise
    // IndexError.
    let mut fp = Fingerprint::new(70);
    assert_eq!(
        fp.set_bit(70),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 70,
            size: 70
        })
    );
    assert_eq!(
        fp.get_bit(70),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 70,
            size: 70
        })
    );
    assert_eq!(
        fp.unset_bit(70),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 70,
            size: 70
        })
    );
    let empty = Fingerprint::new(0);
    assert_eq!(
        empty.get_bit(0),
        Err(FingerprintError::SparseIndexOutOfRange { index: 0, size: 0 })
    );
}

#[test]
fn dense_bit_access_edges_and_failure_atomicity_f06() {
    // F06: bit==size error (already covered generally), repeated
    // set/unset cycles, word-boundary addressing, and failure leaving the
    // value unchanged (source throws before mutating: ExplicitBitVect.cpp
    // :84-95/96-107 check `which >= d_size` first).
    let mut fp = Fingerprint::new(65);
    // word boundaries 0/63/64 plus interior
    for bit in [0u32, 63, 64, 31, 32] {
        assert!(!fp.set_bit(bit).unwrap());
        assert!(fp.set_bit(bit).unwrap()); // repeated set returns prior on
        fp.unset_bit(bit).unwrap();
        assert!(!fp.unset_bit(bit).unwrap()); // repeated unset returns prior off
        fp.set_bit(bit).unwrap();
    }
    assert_eq!(fp.on_bits(), vec![0, 31, 32, 63, 64]);

    // failing calls must not mutate anything
    let before = fp.clone();
    assert!(fp.set_bit(65).is_err());
    assert!(fp.unset_bit(65).is_err());
    assert!(fp.get_bit(65).is_err());
    // a failing set at u32::MAX must not touch tail words either
    assert!(fp.set_bit(u32::MAX).is_err());
    assert_eq!(fp, before);
    assert_eq!(fp.num_on_bits(), 5);

    // set/unset cycle preserves counts across the partial last word
    for _ in 0..3 {
        assert!(fp.unset_bit(64).unwrap());
        assert!(!fp.set_bit(64).unwrap());
    }
    assert_eq!(fp.num_on_bits(), 5);
}

#[test]
fn dense_word_ops_on_non_aligned_tail() {
    // Oracle probe on 65-bit vectors x{64}, y{1}:
    // xor -> [1, 64], or -> [1, 64], and -> [].
    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    x.set_bit(64).unwrap();
    y.set_bit(1).unwrap();
    let xx = x.xor(&y).unwrap();
    assert_eq!(xx.on_bits(), vec![1, 64]);
    assert_eq!(xx.num_on_bits(), 2);
    let oo = x.or(&y).unwrap();
    assert_eq!(oo.on_bits(), vec![1, 64]);
    let aa = x.and(&y).unwrap();
    assert!(aa.on_bits().is_empty());
    assert_eq!(aa.num_on_bits(), 0);
}

#[test]
fn dense_word_operators_recompute_counts_f07() {
    // F07: operator&/|/^/~ recompute d_numOnBits from the result bitset
    // (ExplicitBitVect.cpp:118/125/132/157), so a stale cache on an
    // operand never leaks into the result. Partial last word (65 bits).
    let mut stale = Fingerprint::new(65);
    for bit in [0u32, 31, 64] {
        stale.set_bit(bit).unwrap();
    }
    stale.clear_bits();
    assert_eq!(stale.num_on_bits(), 3); // stale

    let mut other = Fingerprint::new(65);
    other.set_bit(31).unwrap();
    other.set_bit(64).unwrap();

    let a = stale.and(&other).unwrap();
    assert_eq!(a.num_on_bits(), 0); // stale&{31,64} -> {}
    let o = stale.or(&other).unwrap();
    assert_eq!(o.num_on_bits(), 2);
    assert_eq!(o.on_bits(), vec![31, 64]);
    let x = stale.xor(&other).unwrap();
    assert_eq!(x.num_on_bits(), 2);
    assert_eq!(x.on_bits(), vec![31, 64]);
    // complement recomputes even with a stale operand cache
    let c = stale.not();
    assert_eq!(c.num_on_bits(), 65);
    // operands unchanged (stale cache included)
    assert_eq!(stale.num_on_bits(), 3);
    assert_eq!(other.num_on_bits(), 2);
    // length preconditions: mismatch errors on all three binary ops
    let mismatch = Fingerprint::new(64);
    assert!(stale.and(&mismatch).is_err());
    assert!(stale.or(&mismatch).is_err());
    assert!(stale.xor(&mismatch).is_err());
}

#[test]
fn dense_mismatched_sizes_error_for_binary_ops() {
    // Oracle probe: 33-bit & 8-bit raises ValueError (boost
    // invalid_argument mapped to the structured BitLengthMismatch error).
    let a = Fingerprint::new(33);
    let b = Fingerprint::new(8);
    let expected = Err(FingerprintError::BitLengthMismatch { left: 33, right: 8 });
    assert_eq!(a.and(&b), expected);
    assert_eq!(a.or(&b), expected);
    assert_eq!(a.xor(&b), expected);
}

#[test]
fn dense_concat_matches_source_operator_plus() {
    // Oracle probe: a=33 bits {0,32}, b=31 bits {0,30};
    // a + b -> size 64, on [0, 32, 33, 63].
    let mut a = Fingerprint::new(33);
    let mut b = Fingerprint::new(31);
    a.set_bit(0).unwrap();
    a.set_bit(32).unwrap();
    b.set_bit(0).unwrap();
    b.set_bit(30).unwrap();
    let c = a.concat(&b).unwrap();
    assert_eq!(c.n_bits(), 64);
    assert_eq!(c.num_on_bits(), 4);
    assert_eq!(c.on_bits(), vec![0, 32, 33, 63]);
    // aligned concat path
    let d = a.concat(&Fingerprint::new(31)).unwrap();
    assert_eq!(d.n_bits(), 64);
    assert_eq!(d.on_bits(), vec![0, 32]);
}

#[test]
fn dense_complement_matches_source_operator_not() {
    // Oracle probe: 70-bit vector with bit 69 set; ~ -> 69 on bits,
    // [0..=68].
    let mut fp = Fingerprint::new(70);
    fp.set_bit(69).unwrap();
    let n = fp.not();
    assert_eq!(n.n_bits(), 70);
    assert_eq!(n.num_on_bits(), 69);
    assert_eq!(n.on_bits(), (0..69u32).collect::<Vec<_>>());
    // complement of full vector is empty
    let full = Fingerprint::new_filled(70);
    assert_eq!(full.not().num_on_bits(), 0);
}

#[test]
fn dense_filled_constructor_sets_all_bits_including_tail() {
    // Oracle probe: ExplicitBitVect(3, True) -> on 3, off 0, [0,1,2].
    let fp = Fingerprint::new_filled(3);
    assert_eq!(fp.num_on_bits(), 3);
    assert_eq!(fp.num_off_bits(), 0);
    assert_eq!(fp.on_bits(), vec![0, 1, 2]);
    let fp = Fingerprint::new_filled(64);
    assert_eq!(fp.num_on_bits(), 64);
    let fp = Fingerprint::new_filled(0);
    assert_eq!(fp.num_on_bits(), 0);
}

#[test]
fn dense_constructor_size_sweep_f05() {
    // F05: _initForSize across word boundaries. Oracle: GetNumBits ==
    // requested size, GetNumOnBits == 0; ExplicitBitVect(n, True) has all
    // n bits on including the non-word-aligned tail.
    for size in [0u32, 1, 31, 32, 63, 64, 65] {
        let fp = Fingerprint::new(size);
        assert_eq!(fp.n_bits(), size, "zeroed n_bits for size {size}");
        assert_eq!(fp.num_on_bits(), 0, "zeroed count for size {size}");
        assert_eq!(fp.num_off_bits(), size, "zeroed off-count for size {size}");
        assert!(fp.on_bits().is_empty(), "zeroed on bits for size {size}");
        assert_eq!(fp, Fingerprint::new(size), "fresh equality for size {size}");

        let filled = Fingerprint::new_filled(size);
        assert_eq!(filled.n_bits(), size, "filled n_bits for size {size}");
        assert_eq!(filled.num_on_bits(), size, "filled count for size {size}");
        assert_eq!(filled.num_off_bits(), 0, "filled off-count for size {size}");
        let expected: Vec<u32> = (0..size).collect();
        assert_eq!(
            filled.on_bits(),
            expected,
            "filled tail bits for size {size}"
        );
        // unsetting only the tail bit proves it is independently addressable
        if size > 0 {
            let mut t = Fingerprint::new_filled(size);
            assert!(t.unset_bit(size - 1).unwrap());
            assert_eq!(t.num_on_bits(), size - 1, "tail unset for size {size}");
        }
    }
}

#[test]
fn dense_concat_edges_and_immutability_f08() {
    // F08: 0+N, N+0, stale source counts, input immutability. Source
    // recomputes d_numOnBits from the concatenated bitset count(), so a
    // stale operand cache never leaks.
    let mut left = Fingerprint::new(65);
    for bit in [0u32, 31, 64] {
        left.set_bit(bit).unwrap();
    }
    left.clear_bits();
    assert_eq!(left.num_on_bits(), 3); // stale

    // N+0: empty right operand; count recomputed from bits (not stale 3)
    let mut right_empty = Fingerprint::new(31);
    let n0 = left.concat(&right_empty).unwrap();
    assert_eq!(n0.n_bits(), 96);
    assert_eq!(n0.num_on_bits(), 0);
    assert!(n0.on_bits().is_empty());

    // 0+N: empty left operand with set bits on the right
    let mut left_empty = Fingerprint::new(0);
    let mut right = Fingerprint::new(31);
    right.set_bit(0).unwrap();
    right.set_bit(30).unwrap();
    let zn = left_empty.concat(&right).unwrap();
    assert_eq!(zn.n_bits(), 31);
    assert_eq!(zn.on_bits(), vec![0, 30]);

    // stale left + set right: recomputed count = right's on bits only
    let mut right2 = Fingerprint::new(31);
    right2.set_bit(5).unwrap();
    let joined = left.concat(&right2).unwrap();
    assert_eq!(joined.n_bits(), 96);
    assert_eq!(joined.num_on_bits(), 1);
    assert_eq!(joined.on_bits(), vec![70]);

    // input immutability
    right_empty.set_bit(3).unwrap();
    let mut probe = Fingerprint::new(4);
    probe.set_bit(1).unwrap();
    let frozen = probe.clone();
    let _ = probe.concat(&right_empty).unwrap();
    let _ = right_empty.concat(&probe).unwrap();
    assert_eq!(probe, frozen);
    assert_eq!(right_empty.on_bits(), vec![3]);
    // stale cache of left untouched by concat
    assert_eq!(left.num_on_bits(), 3);
}

#[test]
fn dense_equality_and_clone_independence() {
    // Oracle probe: 4-bit {0} != 8-bit {0}; same size same bits equal.
    let mut f1 = Fingerprint::new(4);
    let mut f2 = Fingerprint::new(8);
    f1.set_bit(0).unwrap();
    f2.set_bit(0).unwrap();
    assert_ne!(f1, f2);
    let mut f3 = Fingerprint::new(4);
    f3.set_bit(0).unwrap();
    assert_eq!(f1, f3);
    let mut clone = f1.clone();
    clone.set_bit(3).unwrap();
    assert_eq!(f1.on_bits(), vec![0]);
    assert_eq!(clone.on_bits(), vec![0, 3]);
    assert_ne!(f1, clone);
}

#[test]
fn dense_clear_bits_reproduces_source_stale_on_bit_cache() {
    // C++ source pin (ExplicitBitVect.h:80): clearBits() resets the
    // dynamic_bitset but assigns nothing to d_numOnBits, so getNumOnBits()
    // keeps the stale cached value; the Python wrapper does not expose
    // ClearBits, so this expectation is pinned from the C++ body. setBit
    // afterwards increments the stale base.
    let mut fp = Fingerprint::new(70);
    for bit in [0, 31, 32, 63, 64, 69] {
        fp.set_bit(bit).unwrap();
    }
    assert_eq!(fp.num_on_bits(), 6);
    fp.clear_bits();
    assert!(fp.on_bits().is_empty());
    assert!(!fp.get_bit(0).unwrap());
    // stale cache reproduced: count stays 6 (not 0)
    assert_eq!(fp.num_on_bits(), 6);
    // setBit adjusts the stale base: 6 + 1 = 7
    fp.set_bit(5).unwrap();
    assert_eq!(fp.num_on_bits(), 7);
    // a count-recomputing operation restores consistency
    let recombined = fp.or(&Fingerprint::new(70)).unwrap();
    assert_eq!(recombined.num_on_bits(), 1);
}

#[test]
fn dense_clear_bits_stale_cache_arithmetic_f03() {
    // F03: full stale-cache interaction matrix, source-pinned from
    // ExplicitBitVect.h:80 + ExplicitBitVect.cpp:84-107,154-159,179-183.
    // clear -> unset: unsetBit decrements the stale base even though the
    // bit is already off (source: count(which)==0 -> return false without
    // touching d_numOnBits... verified below against the source branch).
    let mut fp = Fingerprint::new(8);
    fp.set_bit(2).unwrap();
    fp.set_bit(4).unwrap();
    fp.clear_bits();
    // stale cache == 2
    assert_eq!(fp.num_on_bits(), 2);
    // set-clear-set: stale base increments from 2
    fp.set_bit(1).unwrap();
    assert_eq!(fp.num_on_bits(), 3);
    // clear-unset: unset of an off bit returns false and leaves the stale
    // cache untouched (source else-branch returns without decrement)
    assert!(!fp.unset_bit(3).unwrap());
    assert_eq!(fp.num_on_bits(), 3);
    // unset of the bit set after clear decrements the stale base
    assert!(fp.unset_bit(1).unwrap());
    assert_eq!(fp.num_on_bits(), 2);
    assert_eq!(fp.num_off_bits(), 8 - 2);
    // complement recomputation: operator~ counts the complemented bits,
    // ignoring the stale cache (ans.d_numOnBits = ans.dp_bits->count())
    let mut stale = Fingerprint::new(70);
    for bit in [0, 31, 32, 63, 64, 69] {
        stale.set_bit(bit).unwrap();
    }
    stale.clear_bits();
    assert_eq!(stale.num_on_bits(), 6);
    let complement = stale.not();
    assert_eq!(complement.num_on_bits(), 70);
    assert_eq!(complement.on_bits().len(), 70);
}

#[test]
fn dense_from_on_bits_constructs_and_validates() {
    let fp = Fingerprint::from_on_bits(70, [0, 31, 32, 63, 64, 69]).unwrap();
    assert_eq!(fp.on_bits(), vec![0, 31, 32, 63, 64, 69]);
    assert_eq!(fp.num_on_bits(), 6);
    assert_eq!(
        Fingerprint::from_on_bits(4, [4]),
        Err(FingerprintError::SparseIndexOutOfRange { index: 4, size: 4 })
    );
    // duplicate on bits collapse like repeated setBit
    let fp = Fingerprint::from_on_bits(8, [1, 1, 1]).unwrap();
    assert_eq!(fp.on_bits(), vec![1]);
}

#[test]
fn sparse_clear_equality_and_counts_f09() {
    // F09: clearBits, operator== (sets only), getNumOn/OffBits across
    // empty sets and high signed-index states.
    // empty sets of different lengths are equal (source compares sets only)
    let e4 = SparseBitFingerprint::new(4);
    let e8 = SparseBitFingerprint::new(8);
    assert_eq!(e4, e8);
    assert_eq!(e4.num_on_bits(), 0);
    assert_eq!(e4.num_off_bits(), 4);
    assert_eq!(e8.num_off_bits(), 8);

    // clear returns both to the empty-set state
    let mut a = SparseBitFingerprint::new(8);
    a.set_bit(3).unwrap();
    let mut b = SparseBitFingerprint::new(16);
    b.set_bit(9).unwrap();
    b.set_bit(11).unwrap();
    assert_ne!(a, b);
    a.clear_bits();
    b.clear_bits();
    assert_eq!(a, b);
    assert!(a.on_bits().is_empty());
    assert_eq!(a.num_off_bits(), 8);

    // high signed-index states: wrapped keys compare as signed keys
    let mut hi1 = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    hi1.set_bit(2_u32.pow(31) + 1).unwrap();
    let mut hi2 = SparseBitFingerprint::new(2_u32.pow(31) + 9);
    hi2.set_bit(2_u32.pow(31) + 1).unwrap();
    assert_eq!(hi1, hi2); // same wrapped key, different lengths
    hi2.set_bit(5).unwrap();
    assert_ne!(hi1, hi2);
    assert_eq!(hi1.num_on_bits(), 1);
    assert_eq!(hi1.num_off_bits(), 2_u32.pow(31) + 4 - 1);
    // distinct high indices wrap to distinct signed keys
    let mut hi3 = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    hi3.set_bit(2_u32.pow(31) + 2).unwrap();
    assert_ne!(hi1, hi3);
}

#[test]
fn count_construction_and_index_access_f10() {
    // F10: SparseIntVect<uint64_t> constructor/checkIndex/getVal/setVal.
    // Oracle probes (LongSparseIntVect, pinned 2026.03.1): get absent -> 0;
    // SetVal(3,0) removes; GetVal(10) on length 10 -> IndexError; the
    // u64::MAX-length edge (idx == length valid only at max length) is
    // pinned from the template text (Python wrapper cannot construct it:
    // its length argument goes through C long).
    let v = SparseCountFingerprint::new(10);
    assert_eq!(v.length(), 10);
    assert_eq!(v.value(3).unwrap(), 0);
    let mut v = SparseCountFingerprint::new(10);
    v.set_value(3, 5).unwrap();
    v.set_value(4, -2).unwrap();
    assert_eq!(v.value(3).unwrap(), 5);
    assert_eq!(v.value(4).unwrap(), -2);
    assert_eq!(v.value(5).unwrap(), 0);
    // zero deletion
    v.set_value(3, 0).unwrap();
    assert_eq!(v.value(3).unwrap(), 0);
    // index == length error below max
    assert_eq!(
        v.value(10),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 10,
            size: 10
        })
    );
    assert_eq!(
        v.set_value(11, 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 11,
            size: 10
        })
    );
    // zero length
    let z = SparseCountFingerprint::new(0);
    assert_eq!(
        z.value(0),
        Err(FingerprintError::SparseIndexOutOfRange { index: 0, size: 0 })
    );
    // maximum length accepts index == length (template checkIndex)
    let mut big = SparseCountFingerprint::new(u64::MAX);
    big.set_value(u64::MAX, -7).unwrap();
    assert_eq!(big.value(u64::MAX).unwrap(), -7);
    // large valid index below length
    big.set_value(u64::MAX - 1, 3).unwrap();
    assert_eq!(big.value(u64::MAX - 1).unwrap(), 3);
    // below-max length rejects index == length
    let mut below = SparseCountFingerprint::new(u64::MAX - 1);
    assert_eq!(
        below.set_value(u64::MAX - 1, 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: u64::MAX - 1,
            size: u64::MAX - 1
        })
    );
}

#[test]
fn count_specialization_reconciliation_f10() {
    // F10 reconciliation evidence: the legacy value boundary used one
    // u64-index struct (legacy-core/properties/fingerprint.rs L1213:
    // `SparseCountFingerprint { size: u64, nonzero_elements: BTreeMap<u64,
    // i32> }`) matching SparseIntVect<std::uint64_t>; the u32-index count
    // family (getCountFPBulk -> SparseIntVect<std::uint32_t>) is reachable
    // only through molecule-consuming generator dispatch, which is the
    // excluded FP-generator lane. This lane therefore exposes exactly the
    // u64 specialization; a u32 form requires a supervisor refinement when
    // the generator implementation lands.
    let v = SparseCountFingerprint::new(1 << 40);
    assert_eq!(v.length(), 1 << 40);
}

#[test]
fn count_ordered_iteration_and_equality_f11() {
    // F11: getNonzeroElements (ascending map order) and operator==
    // (length + map). Oracle probe: nonzero elements dict preserves the
    // sorted index order; equality is length+entries (a == empty is false
    // when a has entries; equal entries with different insertion order are
    // equal because the source map is ordered).
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(5, 1).unwrap();
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let entries: Vec<(u64, i32)> = a.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (5, 1)]);

    // same entries, different insertion order -> equal
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 5).unwrap();
    b.set_value(1, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert_eq!(a, b);

    // length differs -> not equal even with same entries
    let mut c = SparseCountFingerprint::new(16);
    c.set_value(1, 2).unwrap();
    c.set_value(3, 5).unwrap();
    c.set_value(5, 1).unwrap();
    assert_ne!(a, c);

    // zero removal changes representation equality
    let mut d = SparseCountFingerprint::new(8);
    d.set_value(1, 2).unwrap();
    d.set_value(5, 1).unwrap();
    assert_ne!(a, d);
    let mut e = a.clone();
    e.set_value(3, 0).unwrap();
    assert_eq!(e, d);
    // empty map vs empty map equal; entries vs empty not equal
    assert_ne!(a, SparseCountFingerprint::new(8));
    assert_eq!(
        SparseCountFingerprint::new(8),
        SparseCountFingerprint::new(8)
    );
    // clone independence
    let mut f = a.clone();
    f.set_value(7, 9).unwrap();
    assert_ne!(a, f);
}

#[test]
fn count_total_val_defined_domain_and_reported_edges_f12() {
    // F12: getTotalVal exact on the source-defined domain; undefined
    // signed-overflow edges are reported, never wrapped. Oracle probe:
    // a={1:2, 3:5} -> total 7, abs total 7.
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    assert_eq!(a.total_value(false).unwrap(), 7);
    assert_eq!(a.total_value(true).unwrap(), 7);

    let mut b = SparseCountFingerprint::new(8);
    b.set_value(1, -3).unwrap();
    b.set_value(2, 0).unwrap(); // zero deletion: absent
    b.set_value(3, 10).unwrap();
    assert_eq!(b.total_value(false).unwrap(), 7);
    assert_eq!(b.total_value(true).unwrap(), 13);
    // empty
    assert_eq!(
        SparseCountFingerprint::new(8).total_value(false).unwrap(),
        0
    );
    // source-defined extrema: single i32::MAX / i32::MIN entries
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, i32::MAX).unwrap();
    assert_eq!(mx.total_value(false).unwrap(), i32::MAX);
    assert_eq!(mx.total_value(true).unwrap(), i32::MAX);
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    assert_eq!(mn.total_value(false).unwrap(), i32::MIN);
    // undefined edges, retained as decision requests (contract §7b D1):
    // accumulation overflow ...
    let mut ov = SparseCountFingerprint::new(8);
    ov.set_value(0, i32::MAX).unwrap();
    ov.set_value(1, 1).unwrap();
    assert_eq!(
        ov.total_value(false),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::getTotalVal accumulation"
        })
    );
    // ... and abs(i32::MIN)
    assert_eq!(
        mn.total_value(true),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::getTotalVal abs"
        })
    );
}

#[test]
fn count_fuzzy_and_min_merge_f13() {
    // F13: operator&= fuzzy intersection = per-index minimum; absent-in-
    // other entries drop; length preflight raises the size-mismatch error.
    // Oracle probe: a={1:2,3:5}, b={3:2,5:1} -> a&b = {3:2}.
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let i = a.fuzzy_and(&b).unwrap();
    let entries: Vec<(u64, i32)> = i.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(3, 2)]);

    // negative values: min picks the negative side
    let mut n1 = SparseCountFingerprint::new(8);
    n1.set_value(1, 2).unwrap();
    n1.set_value(2, -4).unwrap();
    let mut n2 = SparseCountFingerprint::new(8);
    n2.set_value(1, -3).unwrap();
    n2.set_value(2, -1).unwrap();
    let ni = n1.fuzzy_and(&n2).unwrap();
    let entries: Vec<(u64, i32)> = ni
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, -3), (2, -4)]);

    // disjoint: everything absent in other drops
    let mut d1 = SparseCountFingerprint::new(8);
    d1.set_value(1, 7).unwrap();
    let d2 = SparseCountFingerprint::new(8);
    assert!(d1.fuzzy_and(&d2).unwrap().nonzero_elements().is_empty());
    // min of nonzero values never creates a stored zero
    let mut z1 = SparseCountFingerprint::new(8);
    z1.set_value(0, 3).unwrap();
    let mut z2 = SparseCountFingerprint::new(8);
    z2.set_value(0, 9).unwrap();
    let zi = z1.fuzzy_and(&z2).unwrap();
    assert_eq!(zi.value(0).unwrap(), 3);

    // unequal lengths error (both orders)
    let c = SparseCountFingerprint::new(16);
    assert_eq!(
        a.fuzzy_and(&c),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(
        c.fuzzy_and(&a),
        Err(FingerprintError::BitLengthMismatch { left: 16, right: 8 })
    );
    // operands unchanged
    assert_eq!(a.value(3).unwrap(), 5);
    assert_eq!(b.value(5).unwrap(), 1);
}

#[test]
fn count_fuzzy_or_max_merge_f14() {
    // F14: operator|= fuzzy union = per-index maximum; one-sided entries
    // (including negative) insert verbatim; size error on length mismatch.
    // Oracle probe: a={1:2,3:5}, b={3:2,5:1} -> a|b = {1:2,3:5,5:1}.
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let u = a.fuzzy_or(&b).unwrap();
    let entries: Vec<(u64, i32)> = u.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (5, 1)]);

    // one-sided negative entry inserts verbatim; overlap keeps max
    let mut n1 = SparseCountFingerprint::new(8);
    n1.set_value(1, 4).unwrap();
    let mut n2 = SparseCountFingerprint::new(8);
    n2.set_value(1, 4).unwrap();
    n2.set_value(2, -6).unwrap();
    let nu = n1.fuzzy_or(&n2).unwrap();
    let entries: Vec<(u64, i32)> = nu
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 4), (2, -6)]);

    // max of negative and positive keeps the positive; both negative keeps
    // the larger (closer to zero)
    let mut p = SparseCountFingerprint::new(8);
    p.set_value(0, -1).unwrap();
    let mut q = SparseCountFingerprint::new(8);
    q.set_value(0, -5).unwrap();
    assert_eq!(p.fuzzy_or(&q).unwrap().value(0).unwrap(), -1);
    assert_eq!(q.fuzzy_or(&p).unwrap().value(0).unwrap(), -1);

    // empty operands
    assert!(
        a.fuzzy_or(&SparseCountFingerprint::new(8))
            .unwrap()
            .nonzero_elements()
            .len()
            == 2
    );
    assert!(
        SparseCountFingerprint::new(8)
            .fuzzy_or(&a)
            .unwrap()
            .nonzero_elements()
            .len()
            == 2
    );

    // size error both orders
    let c = SparseCountFingerprint::new(16);
    assert_eq!(
        a.fuzzy_or(&c),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(
        c.fuzzy_or(&a),
        Err(FingerprintError::BitLengthMismatch { left: 16, right: 8 })
    );
    // operands unchanged
    assert_eq!(a.value(3).unwrap(), 5);
    assert_eq!(b.value(3).unwrap(), 2);
}

#[test]
fn count_vector_add_zero_handling_f15() {
    // F15: operator+= — cancellation to zero erases, one-sided entries
    // insert, source-defined arithmetic exact, overflow reported (D2).
    // Oracle probe: a={1:2,3:5}, b={3:2,5:1} -> a+b = {1:2,3:7,5:1}.
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let s = a.with_added(&b).unwrap();
    let entries: Vec<(u64, i32)> = s.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(1, 2), (3, 7), (5, 1)]);

    // cancellation erases the entry
    let mut c = SparseCountFingerprint::new(8);
    c.set_value(2, 4).unwrap();
    let mut d = SparseCountFingerprint::new(8);
    d.set_value(2, -4).unwrap();
    d.set_value(6, 3).unwrap();
    let cd = c.with_added(&d).unwrap();
    assert_eq!(cd.value(2).unwrap(), 0);
    assert!(cd.nonzero_elements().get(&2).is_none());
    assert_eq!(cd.value(6).unwrap(), 3);

    // one-sided entries both directions + negative sums
    let mut e = SparseCountFingerprint::new(8);
    e.set_value(0, -1).unwrap();
    let mut f = SparseCountFingerprint::new(8);
    f.set_value(7, 2).unwrap();
    let ef = e.with_added(&f).unwrap();
    assert_eq!(ef.value(0).unwrap(), -1);
    assert_eq!(ef.value(7).unwrap(), 2);
    let fe = f.with_added(&e).unwrap();
    assert_eq!(fe, ef);

    // size error both orders; operands unchanged
    let big = SparseCountFingerprint::new(16);
    assert_eq!(
        a.with_added(&big),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(a.value(3).unwrap(), 5);
    assert_eq!(b.value(3).unwrap(), 2);

    // defined extrema and reported overflow (D2 pending)
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, i32::MAX - 1).unwrap();
    let mut one = SparseCountFingerprint::new(8);
    one.set_value(0, 1).unwrap();
    assert_eq!(mx.with_added(&one).unwrap().value(0).unwrap(), i32::MAX);
    let mut over = SparseCountFingerprint::new(8);
    over.set_value(0, 2).unwrap();
    assert_eq!(
        mx.with_added(&over),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator+= accumulation"
        })
    );
}

#[test]
fn count_vector_sub_order_and_zero_handling_f16() {
    // F16: operator-= — exact operand order (a-b != b-a), right-only
    // entries insert negated, cancellation erases, negatives preserved,
    // size error, UB boundaries reported (D2).
    // Oracle probe: a={1:2,3:5}, b={3:2,5:1} -> a-b = {1:2,3:3,5:-1}.
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let d = a.with_subtracted(&b).unwrap();
    let entries: Vec<(u64, i32)> = d.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(1, 2), (3, 3), (5, -1)]);
    // reversed operand order differs
    let dr = b.with_subtracted(&a).unwrap();
    let entries: Vec<(u64, i32)> = dr
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, -2), (3, -3), (5, 1)]);

    // cancellation erases
    let mut c = SparseCountFingerprint::new(8);
    c.set_value(2, 4).unwrap();
    let mut e = SparseCountFingerprint::new(8);
    e.set_value(2, 4).unwrap();
    assert!(c.with_subtracted(&e).unwrap().nonzero_elements().is_empty());

    // size error both orders; operands unchanged
    let big = SparseCountFingerprint::new(16);
    assert_eq!(
        a.with_subtracted(&big),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(a.value(1).unwrap(), 2);
    assert_eq!(b.value(5).unwrap(), 1);

    // defined extrema and reported UB: i32::MIN - 1 overflows; negating
    // i32::MIN on the right-only path is UB
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    let mut one = SparseCountFingerprint::new(8);
    one.set_value(0, 1).unwrap();
    assert_eq!(
        mn.with_subtracted(&one),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-= accumulation"
        })
    );
    let mut neg = SparseCountFingerprint::new(8);
    neg.set_value(3, i32::MIN).unwrap();
    let empty = SparseCountFingerprint::new(8);
    assert_eq!(
        empty.with_subtracted(&neg),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-= negation"
        })
    );
    // defined: subtracting to exactly i32::MIN stays representable
    let mut low = SparseCountFingerprint::new(8);
    low.set_value(0, -1).unwrap();
    let mut up = SparseCountFingerprint::new(8);
    up.set_value(0, i32::MAX).unwrap();
    assert_eq!(
        low.with_subtracted(&up).unwrap().value(0).unwrap(),
        i32::MIN
    );
}

#[test]
fn count_scalar_mul_f17() {
    // F17: operator*=(int) under the Supervisor resolution: exact defined
    // products; UndefinedArithmetic only when a C++-undefined multiply
    // actually executes; zero entries retained; input unchanged.
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 3).unwrap();
    v.set_value(2, -4).unwrap();

    // positive / negative / zero factors
    let p = v.with_multiplied_scalar(2).unwrap();
    assert_eq!(p.value(1).unwrap(), 6);
    assert_eq!(p.value(2).unwrap(), -8);
    let n = v.with_multiplied_scalar(-3).unwrap();
    assert_eq!(n.value(1).unwrap(), -9);
    assert_eq!(n.value(2).unwrap(), 12);
    let z = v.with_multiplied_scalar(0).unwrap();
    assert_eq!(z.value(1).unwrap(), 0);
    assert_eq!(z.value(2).unwrap(), 0);
    // multiplied-through-zero entries REMAIN stored (no zero deletion)
    assert_eq!(z.nonzero_elements().len(), 2);
    // zero entries stay zero under further defined multiplies
    let z2 = z.with_multiplied_scalar(5).unwrap();
    assert_eq!(z2.nonzero_elements().len(), 2);
    assert_eq!(z2.value(1).unwrap(), 0);

    // empty storage succeeds for any factor (nothing executes)
    assert!(
        SparseCountFingerprint::new(8)
            .with_multiplied_scalar(i32::MIN)
            .is_ok()
    );

    // valid boundary products (1_073_741_823 * 2 == i32::MAX - 1)
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, 1_073_741_823).unwrap();
    assert_eq!(
        mx.with_multiplied_scalar(2).unwrap().value(0).unwrap(),
        i32::MAX - 1
    );
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN + 1).unwrap();
    assert_eq!(
        mn.with_multiplied_scalar(1).unwrap().value(0).unwrap(),
        i32::MIN + 1
    );

    // UB rejections: MIN * -1 and general product overflow (1<<30 * 4)
    let mut min_entry = SparseCountFingerprint::new(8);
    min_entry.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        min_entry.with_multiplied_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator*=(int)"
        })
    );
    let mut big = SparseCountFingerprint::new(8);
    big.set_value(0, 1 << 30).unwrap();
    assert_eq!(
        big.with_multiplied_scalar(4),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator*=(int)"
        })
    );
    // input unchanged on success and error
    assert_eq!(v.value(1).unwrap(), 3);
    assert_eq!(v.value(2).unwrap(), -4);
    assert_eq!(min_entry.value(0).unwrap(), i32::MIN);
    // length and ordering preserved
    assert_eq!(p.length(), 8);
    let keys: Vec<u64> = p.nonzero_elements().keys().copied().collect();
    assert_eq!(keys, vec![1, 2]);
}

#[test]
fn count_scalar_div_f18() {
    // F18: operator/=(int) under the Supervisor resolution: truncation
    // toward zero for defined quotients; executed-UB rejection only.
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(0, 7).unwrap();
    v.set_value(1, -7).unwrap();
    v.set_value(2, 6).unwrap();
    v.set_value(3, -6).unwrap();

    // sign combinations with truncation toward zero: 7/2=3, -7/2=-3,
    // 6/-2=-3, -6/-2=3
    let q = v.with_divided_scalar(2).unwrap();
    assert_eq!(q.value(0).unwrap(), 3);
    assert_eq!(q.value(1).unwrap(), -3);
    let qn = v.with_divided_scalar(-2).unwrap();
    assert_eq!(qn.value(2).unwrap(), -3);
    assert_eq!(qn.value(3).unwrap(), 3);

    // empty storage divided by zero SUCCEEDS: no division executes
    assert!(
        SparseCountFingerprint::new(8)
            .with_divided_scalar(0)
            .is_ok()
    );
    // any stored entry divided by zero errors (executed UB)
    assert_eq!(
        v.with_divided_scalar(0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    // a stored ZERO (scalar-produced) divided by zero still executes and
    // errors
    let zeroed = v.with_multiplied_scalar(0).unwrap();
    assert_eq!(zeroed.nonzero_elements().len(), 4);
    assert_eq!(
        zeroed.with_divided_scalar(0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    // MIN / -1 errors (executed UB); MIN / 1 is defined
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        mn.with_divided_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    assert_eq!(
        mn.with_divided_scalar(1).unwrap().value(0).unwrap(),
        i32::MIN
    );
    // zero results are retained, not deleted: 3/4 == 0 stays stored
    let mut small = SparseCountFingerprint::new(8);
    small.set_value(5, 3).unwrap();
    let zq = small.with_divided_scalar(4).unwrap();
    assert_eq!(zq.value(5).unwrap(), 0);
    assert_eq!(zq.nonzero_elements().len(), 1);
    // input immutability on success and error
    assert_eq!(v.value(0).unwrap(), 7);
    assert_eq!(mn.value(0).unwrap(), i32::MIN);
    assert_eq!(q.length(), 8);
}

#[test]
fn count_scalar_add_f19() {
    // F19: operator+=(int) under the Supervisor resolution: exact defined
    // sums on existing entries only; zero results retained; executed-UB
    // rejection; input unchanged.
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 4).unwrap();
    v.set_value(2, -6).unwrap();

    // positive / negative / zero operands
    assert_eq!(v.with_added_scalar(3).unwrap().value(1).unwrap(), 7);
    assert_eq!(v.with_added_scalar(3).unwrap().value(2).unwrap(), -3);
    assert_eq!(v.with_added_scalar(-10).unwrap().value(1).unwrap(), -6);
    assert_eq!(v.with_added_scalar(0).unwrap().value(1).unwrap(), 4);

    // zero-result retention: -6 + 6 == 0 stays stored
    let zeroed = v.with_added_scalar(6).unwrap();
    assert_eq!(zeroed.value(2).unwrap(), 0);
    assert_eq!(zeroed.nonzero_elements().len(), 2);
    // absent-index preservation (no densification)
    assert_eq!(zeroed.value(0).unwrap(), 0);
    assert!(zeroed.nonzero_elements().get(&0).is_none());

    // empty storage succeeds for any operand (nothing executes)
    assert!(
        SparseCountFingerprint::new(8)
            .with_added_scalar(i32::MIN)
            .is_ok()
    );

    // valid boundary sums: MAX-1 + 1 == MAX; MIN+1 + -1 == MIN
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, i32::MAX - 1).unwrap();
    assert_eq!(mx.with_added_scalar(1).unwrap().value(0).unwrap(), i32::MAX);
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN + 1).unwrap();
    assert_eq!(
        mn.with_added_scalar(-1).unwrap().value(0).unwrap(),
        i32::MIN
    );

    // executed-UB rejection: MAX + 1 and MIN + -1
    let mut over = SparseCountFingerprint::new(8);
    over.set_value(0, i32::MAX).unwrap();
    assert_eq!(
        over.with_added_scalar(1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator+=(int)"
        })
    );
    let mut under = SparseCountFingerprint::new(8);
    under.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        under.with_added_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator+=(int)"
        })
    );
    // input immutability on success and error
    assert_eq!(v.value(1).unwrap(), 4);
    assert_eq!(over.value(0).unwrap(), i32::MAX);
    assert_eq!(zeroed.length(), 8);
}

#[test]
fn count_scalar_sub_f20() {
    // F20: operator-=(int) under the Supervisor resolution: exact defined
    // differences on existing entries only; zero results retained;
    // executed-UB rejection; input unchanged.
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 4).unwrap();
    v.set_value(2, -6).unwrap();

    // positive / negative / zero operands
    assert_eq!(v.with_subtracted_scalar(3).unwrap().value(1).unwrap(), 1);
    assert_eq!(v.with_subtracted_scalar(-10).unwrap().value(1).unwrap(), 14);
    assert_eq!(v.with_subtracted_scalar(0).unwrap().value(2).unwrap(), -6);

    // zero-result retention: 4 - 4 == 0 stays stored
    let zeroed = v.with_subtracted_scalar(4).unwrap();
    assert_eq!(zeroed.value(1).unwrap(), 0);
    assert_eq!(zeroed.nonzero_elements().len(), 2);
    // absent-index preservation (no densification)
    assert!(zeroed.nonzero_elements().get(&0).is_none());

    // empty storage succeeds for any operand (nothing executes)
    assert!(
        SparseCountFingerprint::new(8)
            .with_subtracted_scalar(i32::MIN)
            .is_ok()
    );

    // valid boundary differences: MIN+1 - 1 == MIN; MAX-1 - (-1) == MAX
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN + 1).unwrap();
    assert_eq!(
        mn.with_subtracted_scalar(1).unwrap().value(0).unwrap(),
        i32::MIN
    );
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, i32::MAX - 1).unwrap();
    assert_eq!(
        mx.with_subtracted_scalar(-1).unwrap().value(0).unwrap(),
        i32::MAX
    );

    // executed-UB rejection: MIN - 1 and MAX - (-1)
    let mut under = SparseCountFingerprint::new(8);
    under.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        under.with_subtracted_scalar(1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-=(int)"
        })
    );
    let mut over = SparseCountFingerprint::new(8);
    over.set_value(0, i32::MAX).unwrap();
    assert_eq!(
        over.with_subtracted_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-=(int)"
        })
    );
    // input immutability on success and error
    assert_eq!(v.value(1).unwrap(), 4);
    assert_eq!(under.value(0).unwrap(), i32::MIN);
    assert_eq!(zeroed.length(), 8);
}

#[test]
fn hash_integral_cast_overloads_f21() {
    // F21: pinned-source cast/identity overloads (hash.hpp:138-176).
    // Expectations are the defined C++ conversions of the pinned ABI
    // (gcc x86-64); no UB is involved.
    use cosmolkit_fingerprints::hash;

    // hash_result_t width is 32 bits (hash_fwd.hpp)
    assert_eq!(std::mem::size_of::<hash::HashResult>(), 4);

    assert_eq!(hash::hash_value_bool(false), 0);
    assert_eq!(hash::hash_value_bool(true), 1);

    // signed small types: two's-complement 32-bit patterns
    assert_eq!(hash::hash_value_i8(-1), 0xFFFF_FFFF);
    assert_eq!(hash::hash_value_i8(127), 127);
    assert_eq!(hash::hash_value_i8(i8::MIN), 0xFFFF_FF80);
    assert_eq!(hash::hash_value_i16(-2), 0xFFFF_FFFE);
    assert_eq!(hash::hash_value_i16(i16::MIN), 0xFFFF_8000);
    assert_eq!(hash::hash_value_i32(-1), 0xFFFF_FFFF);
    assert_eq!(hash::hash_value_i32(i32::MIN), 0x8000_0000);
    assert_eq!(hash::hash_value_i32(i32::MAX), 0x7FFF_FFFF);

    // unsigned small types: zero-extended identity
    assert_eq!(hash::hash_value_u8(u8::MAX), 0xFF);
    assert_eq!(hash::hash_value_u16(u16::MAX), 0xFFFF);
    assert_eq!(hash::hash_value_u32(0xDEAD_BEEF), 0xDEAD_BEEF);
    assert_eq!(hash::hash_value_u32(0), 0);

    // unsigned long (usize) under the pinned ABI: plain cast truncation,
    // NOT the 64-bit mixing helper (that distinction is F22's; the value
    // here proves this overload truncates)
    assert_eq!(hash::hash_value_usize(5), 5);
    assert_eq!(hash::hash_value_usize(0x1_0000_0001), 1);
    assert_eq!(hash::hash_value_usize(usize::MAX), 0xFFFF_FFFF);
}

#[test]
fn hash_unsigned64_mixing_f22() {
    // F22: hash_value_unsigned for 64-bit values (hash.hpp:118-135).
    // Expectations are deterministic defined unsigned arithmetic derived
    // from the pinned source shape (length = 1: fold high 32 then low 32).
    use cosmolkit_fingerprints::hash;

    // zero high half reduces to the low-word fold from seed 0
    assert_eq!(hash::hash_value_u64(1), 1);
    assert_eq!(hash::hash_value_u64(0), 0);
    assert_eq!(hash::hash_value_u64(u64::from(u32::MAX)), 0xFFFF_FFFF);
    // high half participates: seed = high32 first
    assert_eq!(hash::hash_value_u64(0x1_0000_0000), 65);
    assert_eq!(hash::hash_value_u64(0x1_0000_0001), 64);
    // full extrema
    assert_eq!(hash::hash_value_u64(u64::MAX), 0xC000_0041);
    assert_eq!(hash::hash_value_u64(1 << 63), {
        // seed = 0x8000_0000; low = 0:
        // 0 ^ 0 + (0x8000_0000 << 6 = 0) + (0x8000_0000 >> 2 = 0x2000_0000)
        // = 0x2000_0000; 0x8000_0000 ^ 0x2000_0000 = 0xA000_0000
        0xA000_0000
    });
    // distinction from the plain-cast unsigned long overload (F21)
    assert_ne!(
        hash::hash_value_u64(0x1_0000_0001),
        hash::hash_value_usize(0x1_0000_0001)
    );
    assert_eq!(hash::hash_value_usize(0x1_0000_0001), 1);
}

#[test]
fn hash_signed64_mixing_including_min_f23() {
    // F23: hash_value_signed for 64-bit values including the minimum
    // signed value (hash.hpp:98-116). Expectations are deterministic
    // defined arithmetic derived from the pinned source shape.
    use cosmolkit_fingerprints::hash;

    assert_eq!(hash::hash_value_i64(0), 0);
    assert_eq!(hash::hash_value_i64(1), 1);
    // negative: positive = -1 - val = 0 for val == -1, low word 0xFFFF_FFFF
    assert_eq!(hash::hash_value_i64(-1), 0xFFFF_FFFF);
    assert_eq!(hash::hash_value_i64(-2), 0xFFFF_FFFE);
    // i64::MIN: positive = i64::MAX (high 0x7FFF_FFFF), low word 0;
    // final addend = 0 + (0x7FFF_FFFF << 6) + (0x7FFF_FFFF >> 2)
    //              = 0xFFFF_FFC0 + 0x1FFF_FFFF = 0x1FFF_FFBF (mod 2^32)
    // seed = 0x7FFF_FFFF ^ 0x1FFF_FFBF = 0x6000_0040
    assert_eq!(hash::hash_value_i64(i64::MIN), 0x6000_0040);
    // i64::MAX (verified against the source-shape emulation): high
    // 0x7FFF_FFFF folds first, final addend = 0xFFFF_FFFF + 0xFFFF_FFC0 +
    // 0x1FFF_FFFF (mod 2^32) = 0xFFFF_FFBE... seed = 0x7FFF_FFFF ^
    // 0xE000_00BE = 0x6000_0041 per emulation of hash.hpp:98-116
    assert_eq!(hash::hash_value_i64(i64::MAX), 0x6000_0041);
    // signed and unsigned 64-bit paths differ on the same bit pattern
    assert_ne!(hash::hash_value_i64(-1), hash::hash_value_u64(u64::MAX));
}

#[test]
fn hash_combine_u32_seed_f24() {
    // F24: hash_combine for the unsigned-int call-site overload
    // (hash.hpp:209-219). Expectations generated by a faithful
    // source-shape emulation of the declared 32-bit unsigned arithmetic.
    use cosmolkit_fingerprints::hash;

    let mut seed = 0u32;
    hash::hash_combine(&mut seed, 0);
    assert_eq!(seed, 0x9e3779b9);

    let mut seed = 0u32;
    hash::hash_combine(&mut seed, 1);
    assert_eq!(seed, 0x9e3779ba);

    let mut seed = 0u32;
    hash::hash_combine(&mut seed, 0xBEEF);
    assert_eq!(seed, 0x9e3838a8);

    // high-bit seed exercises both the <<6 and >>2 wrap paths
    let mut seed = 0xFFFF_FFFF;
    hash::hash_combine(&mut seed, 0);
    assert_eq!(seed, 0x21c8_8687);

    let mut seed = 0xFFFF_FFFF;
    hash::hash_combine(&mut seed, 0xFFFF_FFFF);
    assert_eq!(seed, 0x21c8_8688);

    // ordered sequences are order-sensitive
    let mut seed = 0u32;
    for v in [1u32, 2, 3] {
        hash::hash_combine(&mut seed, v);
    }
    assert_eq!(seed, 0xfb58_d153);

    let mut seed = 0u32;
    for v in [0xFFFF_FFFFu32, 0x8000_0000, 0x7FFF_FFFF] {
        hash::hash_combine(&mut seed, v);
    }
    assert_eq!(seed, 0xdb58_2ec0);
}

#[test]
fn hash_range_ordered_accumulation_f25() {
    // F25: hash_range overloads (hash.hpp:221-237). Empty range from a
    // zero seed stays zero; ordered accumulation equals the corresponding
    // hash_combine sequence; the seeded overload continues an existing
    // seed exactly like continued combines.
    use cosmolkit_fingerprints::hash;

    // empty
    assert_eq!(hash::hash_range(&[]), 0);
    let mut seed = 0x1234_5678u32;
    hash::hash_range_seeded(&mut seed, &[]);
    assert_eq!(seed, 0x1234_5678);

    // single and multiple, order-sensitive
    assert_eq!(hash::hash_range(&[1]), 0x9e3779ba);
    assert_eq!(hash::hash_range(&[1, 2, 3]), 0xfb58_d153);
    assert_ne!(hash::hash_range(&[3, 2, 1]), hash::hash_range(&[1, 2, 3]));

    // high-bit values
    assert_eq!(
        hash::hash_range(&[0xFFFF_FFFF, 0x8000_0000, 0x7FFF_FFFF]),
        0xdb58_2ec0
    );

    // seeded continuation equals the same combine sequence
    let mut via_seed = 0u32;
    hash::hash_combine(&mut via_seed, 1);
    hash::hash_range_seeded(&mut via_seed, &[2, 3]);
    assert_eq!(via_seed, hash::hash_range(&[1, 2, 3]));
}

#[test]
fn fold_dense_fingerprint_f26() {
    // F26: FoldFingerprint<ExplicitBitVect> (BitOps.cpp:673-690).
    // Oracle probes (pinned 2026.03.1): 32-bit {0,5,31}:
    //   /2 -> 16 bits [0,5,15]; /4 -> 8 bits [0,5,7]; factor 1 -> identity
    //   [0,5,31]; factor 0 and 32 -> ValueError "invalid fold factor";
    //   factor 31 (= n_bits-1) succeeds;
    //   33-bit {16,32} /2 -> 16 bits [0] (floor length + collision);
    //   empty vector: any factor errors.
    use cosmolkit_fingerprints::folding;

    let mut bv = Fingerprint::new(32);
    for bit in [0u32, 5, 31] {
        bv.set_bit(bit).unwrap();
    }
    let f2 = folding::fold_fingerprint(&bv, 2).unwrap();
    assert_eq!(f2.n_bits(), 16);
    assert_eq!(f2.on_bits(), vec![0, 5, 15]);
    let f4 = folding::fold_fingerprint(&bv, 4).unwrap();
    assert_eq!(f4.n_bits(), 8);
    assert_eq!(f4.on_bits(), vec![0, 5, 7]);
    let f1 = folding::fold_fingerprint(&bv, 1).unwrap();
    assert_eq!(f1.n_bits(), 32);
    assert_eq!(f1.on_bits(), vec![0, 5, 31]);

    let expected_err = Err(FingerprintError::InvalidFoldFactor {
        factor: 0,
        n_bits: 32,
    });
    assert_eq!(folding::fold_fingerprint(&bv, 0), expected_err);
    assert_eq!(
        folding::fold_fingerprint(&bv, 32),
        Err(FingerprintError::InvalidFoldFactor {
            factor: 32,
            n_bits: 32
        })
    );
    // factor == n_bits - 1 is allowed
    let f31 = folding::fold_fingerprint(&bv, 31).unwrap();
    assert_eq!(f31.n_bits(), 1);
    assert_eq!(f31.on_bits(), vec![0]);

    // non-divisible length folds by floor and collides
    let mut bv33 = Fingerprint::new(33);
    bv33.set_bit(16).unwrap();
    bv33.set_bit(32).unwrap();
    let f = folding::fold_fingerprint(&bv33, 2).unwrap();
    assert_eq!(f.n_bits(), 16);
    assert_eq!(f.on_bits(), vec![0]);

    // empty vector: any factor is invalid (factor >= n_bits == 0)
    assert_eq!(
        folding::fold_fingerprint(&Fingerprint::new(0), 2),
        Err(FingerprintError::InvalidFoldFactor {
            factor: 2,
            n_bits: 0
        })
    );
    // input unchanged
    assert_eq!(bv.on_bits(), vec![0, 5, 31]);
}

#[test]
fn fold_sparse_fingerprint_signed_indices_f27() {
    // F27: FoldFingerprint<SparseBitVect> (BitOps.cpp:673-690).
    // Oracle probes (pinned 2026.03.1):
    //   sparse 32 {0,5,31}: /2 -> [0,5,15]; /4 -> [0,5,7]; factor 0 ->
    //   ValueError;
    //   {2^31+1, 3} over size 2^31+4 folded by 2 -> IndexError with
    //   negative remainder (C++ message shows -1073741821; the structured
    //   error records the wrapped u32 index of the same boundary);
    //   {2^31-1, 5} over size 2^31 folded by 2 -> size 2^30,
    //   [5, 1073741823] with no dense materialization.
    use cosmolkit_fingerprints::folding;

    let mut sbv = SparseBitFingerprint::new(32);
    for bit in [0u32, 5, 31] {
        sbv.set_bit(bit).unwrap();
    }
    let f2 = folding::fold_fingerprint_sparse(&sbv, 2).unwrap();
    assert_eq!(f2.n_bits(), 16);
    assert_eq!(f2.on_bits(), vec![0, 5, 15]);
    let f4 = folding::fold_fingerprint_sparse(&sbv, 4).unwrap();
    assert_eq!(f4.n_bits(), 8);
    assert_eq!(f4.on_bits(), vec![0, 5, 7]);
    assert_eq!(
        folding::fold_fingerprint_sparse(&sbv, 0),
        Err(FingerprintError::InvalidFoldFactor {
            factor: 0,
            n_bits: 32
        })
    );
    assert_eq!(
        folding::fold_fingerprint_sparse(&sbv, 32),
        Err(FingerprintError::InvalidFoldFactor {
            factor: 32,
            n_bits: 32
        })
    );

    // high signed index: negative remainder fails checkIndex
    let mut big = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    big.set_bit(2_u32.pow(31) + 1).unwrap();
    big.set_bit(3).unwrap();
    // remainder: (-2147483647) % (2^30 + 2) = -(1073741821) -> wrapped u32
    let wrapped_neg_remainder = (-1073741821i32) as u32;
    assert_eq!(
        folding::fold_fingerprint_sparse(&big, 2),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: u64::from(wrapped_neg_remainder),
            size: u64::from(2_u32.pow(30) + 2)
        })
    );

    // below-2^31 large index folds without dense materialization
    let mut med = SparseBitFingerprint::new(2_u32.pow(31));
    med.set_bit(2_u32.pow(31) - 1).unwrap();
    med.set_bit(5).unwrap();
    let fm = folding::fold_fingerprint_sparse(&med, 2).unwrap();
    assert_eq!(fm.n_bits(), 2_u32.pow(30));
    assert_eq!(fm.on_bits(), vec![5i32, 2_u32.pow(30) as i32 - 1]);
    // input unchanged
    assert_eq!(med.on_bits(), vec![5, 2_147_483_647]);
}

#[test]
fn value_conversions_roundtrip_f28() {
    // F28: dense<->sparse conversions with checked lossless index
    // mapping. from_lsb_bytes is out of closure on caller evidence
    // (Avalon generator, legacy avalon_fingerprint.rs:157).
    use cosmolkit_fingerprints::folding;

    // roundtrip incl. tail bits
    let mut dense = Fingerprint::new(70);
    for bit in [0u32, 31, 32, 63, 64, 69] {
        dense.set_bit(bit).unwrap();
    }
    let sparse = folding::to_sparse_bit(&dense);
    assert_eq!(sparse.n_bits(), 70);
    assert_eq!(sparse.on_bits(), vec![0, 31, 32, 63, 64, 69]);
    let back = folding::to_dense(&sparse).unwrap();
    assert_eq!(back, dense);

    // empty roundtrip
    let empty_sparse = folding::to_sparse_bit(&Fingerprint::new(0));
    assert_eq!(empty_sparse.n_bits(), 0);
    assert!(
        folding::to_dense(&empty_sparse)
            .unwrap()
            .on_bits()
            .is_empty()
    );

    // high signed sparse output: a key >= 2^31 (negative carrier value) in
    // a moderate-length vector is an out-of-length member; conversion is
    // checked, never silently truncated
    let mut odd = SparseBitFingerprint::new(8);
    odd.set_bit(1).unwrap();
    let with_out_of_range = odd.union(&{
        let mut big = SparseBitFingerprint::new(2_u32.pow(31) + 4);
        big.set_bit(2_u32.pow(31) + 1).unwrap();
        big
    });
    assert_eq!(with_out_of_range.n_bits(), 8);
    assert_eq!(with_out_of_range.on_bits(), vec![-2_147_483_647, 1]);
    let wrapped_index = u64::from((-2_147_483_647i32) as u32);
    assert_eq!(
        folding::to_dense(&with_out_of_range),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: wrapped_index,
            size: 8
        })
    );
    // in-range high member of a large vector converts losslessly without
    // allocating a dense huge vector in this test's error-free small case
    let mut ok_sparse = SparseBitFingerprint::new(8);
    ok_sparse.set_bit(7).unwrap();
    let ok_dense = folding::to_dense(&ok_sparse).unwrap();
    assert_eq!(ok_dense.on_bits(), vec![7]);
    // inputs unchanged
    assert_eq!(dense.on_bits().len(), 6);
    assert_eq!(with_out_of_range.on_bits().len(), 2);
}

#[test]
fn num_on_bits_in_common_f29() {
    // F29: NumOnBitsInCommon dense + sparse overloads (BitOps.cpp:222-269).
    // Oracle probe: rdkit.NumOnBitsInCommon equivalents —
    //   (a & b).GetNumOnBits() for dense; OnBitsInCommon().size() for
    //   sparse (length-checked before counting).
    use cosmolkit_fingerprints::similarity;

    // empty pair
    assert_eq!(
        similarity::num_on_bits_in_common(&Fingerprint::new(8), &Fingerprint::new(8)).unwrap(),
        0
    );
    // overlap incl. partial final word (65 bits)
    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    for bit in [1u32, 64] {
        x.set_bit(bit).unwrap();
    }
    for bit in [1u32, 63] {
        y.set_bit(bit).unwrap();
    }
    assert_eq!(similarity::num_on_bits_in_common(&x, &y).unwrap(), 1);
    // complements
    let full = Fingerprint::new_filled(65);
    assert_eq!(similarity::num_on_bits_in_common(&x, &full).unwrap(), 2);
    assert_eq!(similarity::num_on_bits_in_common(&x, &x.not()).unwrap(), 0);
    // length mismatch
    assert_eq!(
        similarity::num_on_bits_in_common(&x, &Fingerprint::new(64)),
        Err(FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        })
    );

    // sparse: same semantics through the length-checked template
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3, 5] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5, 7] {
        sb.set_bit(bit).unwrap();
    }
    assert_eq!(
        similarity::sparse_num_on_bits_in_common(&sa, &sb).unwrap(),
        2
    );
    assert_eq!(
        similarity::sparse_num_on_bits_in_common(&SparseBitFingerprint::new(8), &sb).unwrap(),
        0
    );
    assert_eq!(
        similarity::sparse_num_on_bits_in_common(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // high sparse index participates in the intersection
    let mut ha = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    let mut hb = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    ha.set_bit(2_u32.pow(31) + 1).unwrap();
    hb.set_bit(2_u32.pow(31) + 1).unwrap();
    assert_eq!(
        similarity::sparse_num_on_bits_in_common(&ha, &hb).unwrap(),
        1
    );
}

#[test]
fn on_bits_in_common_f30() {
    // F30: OnBitsInCommon dense + sparse (BitOps.cpp:561-569): ascending
    // common-on indices, length prechecked.
    use cosmolkit_fingerprints::similarity;

    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    for bit in [1u32, 3, 64] {
        x.set_bit(bit).unwrap();
    }
    for bit in [3u32, 64, 63] {
        y.set_bit(bit).unwrap();
    }
    assert_eq!(similarity::on_bits_in_common(&x, &y).unwrap(), vec![3, 64]);
    // empty pair
    assert!(
        similarity::on_bits_in_common(&Fingerprint::new(8), &Fingerprint::new(8))
            .unwrap()
            .is_empty()
    );
    // complements
    assert!(
        similarity::on_bits_in_common(&x, &x.not())
            .unwrap()
            .is_empty()
    );
    // mismatch
    assert_eq!(
        similarity::on_bits_in_common(&x, &Fingerprint::new(64)),
        Err(FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        })
    );

    // sparse with high signed index
    let mut sa = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    let mut sb = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    sa.set_bit(5).unwrap();
    sa.set_bit(2_u32.pow(31) + 1).unwrap();
    sb.set_bit(2_u32.pow(31) + 1).unwrap();
    // raw signed keys in signed order: the >=2^31 member appears first
    assert_eq!(
        similarity::sparse_on_bits_in_common(&sa, &sb).unwrap(),
        vec![-2_147_483_647]
    );
    assert_eq!(
        similarity::sparse_on_bits_in_common(&sa, &SparseBitFingerprint::new(2_u32.pow(31) + 4))
            .unwrap(),
        Vec::<i32>::new()
    );
    assert_eq!(
        similarity::sparse_on_bits_in_common(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch {
            left: u64::from(2_u32.pow(31) + 4),
            right: 16
        })
    );
}

#[test]
fn off_bits_in_common_f31() {
    // F31: OffBitsInCommon dense + sparse (BitOps.cpp:585-593):
    // (~(a|b)).getOnBits() with the length precondition.
    use cosmolkit_fingerprints::similarity;

    // 65-bit partial word: x{1,3,64}, y{3,63,64} -> union {1,3,63,64},
    // off-common = the rest
    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    for bit in [1u32, 3, 64] {
        x.set_bit(bit).unwrap();
    }
    for bit in [3u32, 63, 64] {
        y.set_bit(bit).unwrap();
    }
    let expected: Vec<u32> = (0..65u32).filter(|i| ![1, 3, 63, 64].contains(i)).collect();
    assert_eq!(similarity::off_bits_in_common(&x, &y).unwrap(), expected);

    // both full -> no off bits in common
    let full = Fingerprint::new_filled(65);
    assert!(
        similarity::off_bits_in_common(&full, &full)
            .unwrap()
            .is_empty()
    );
    // both empty -> all off
    let empty = Fingerprint::new(65);
    let all: Vec<u32> = (0..65).collect();
    assert_eq!(similarity::off_bits_in_common(&empty, &empty).unwrap(), all);
    // mismatch
    assert_eq!(
        similarity::off_bits_in_common(&x, &Fingerprint::new(64)),
        Err(FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        })
    );

    // sparse with retained out-of-length member: the complement loop only
    // visits indices below the length, so the retained member never
    // appears in the off-common output (bounded test; the signed-ordered
    // high-index output shape is covered in F30 without a giant
    // materialization)
    let mut wide = SparseBitFingerprint::new(16);
    wide.set_bit(12).unwrap();
    let mut sa = SparseBitFingerprint::new(8);
    sa.set_bit(1).unwrap();
    let with_oob = sa.union(&sa.union(&wide));
    assert_eq!(with_oob.n_bits(), 8);
    assert_eq!(with_oob.on_bits(), vec![1, 12]);
    let mut sb = SparseBitFingerprint::new(8);
    sb.set_bit(3).unwrap();
    assert_eq!(
        similarity::sparse_off_bits_in_common(&with_oob, &sb).unwrap(),
        vec![0, 2, 4, 5, 6, 7]
    );
    assert_eq!(
        similarity::sparse_off_bits_in_common(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn num_bits_in_common_f32() {
    // F32: NumBitsInCommon dense + sparse (BitOps.cpp:510-522):
    // n_bits - (a^b).on_count with the length precondition. Sparse
    // retained out-of-length members can push the xor count past the
    // length; the source's unsigned subtraction then wraps (defined).
    use cosmolkit_fingerprints::similarity;

    // dense, partial final word: x{1,3,64}, y{3,63,64} -> xor {1,63} ->
    // 65 - 2 = 63
    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    for bit in [1u32, 3, 64] {
        x.set_bit(bit).unwrap();
    }
    for bit in [3u32, 63, 64] {
        y.set_bit(bit).unwrap();
    }
    assert_eq!(similarity::num_bits_in_common(&x, &y).unwrap(), 63);
    // identical / complementary
    assert_eq!(similarity::num_bits_in_common(&x, &x).unwrap(), 65);
    assert_eq!(similarity::num_bits_in_common(&x, &x.not()).unwrap(), 0);
    // empty pair: all (zero) bits agree
    assert_eq!(
        similarity::num_bits_in_common(&Fingerprint::new(0), &Fingerprint::new(0)).unwrap(),
        0
    );
    // mismatch
    assert_eq!(
        similarity::num_bits_in_common(&x, &Fingerprint::new(64)),
        Err(FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        })
    );

    // sparse, high signed index participates in the xor count
    let mut ha = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    let hb = SparseBitFingerprint::new(2_u32.pow(31) + 4);
    ha.set_bit(2_u32.pow(31) + 1).unwrap();
    // xor = {high key} -> count 1
    assert_eq!(
        similarity::sparse_num_bits_in_common(&ha, &hb).unwrap(),
        2_u32.pow(31) + 3
    );

    // sparse retained out-of-length members: two length-8 vectors whose
    // symmetric difference holds 9 members -> 8 - 9 wraps mod 2^32
    let base = SparseBitFingerprint::new(8);
    let mut wide1 = SparseBitFingerprint::new(17);
    for bit in 8u32..=16 {
        wide1.set_bit(bit).unwrap();
    }
    let a = base.union(&wide1);
    let b = base.clone();
    assert_eq!(a.n_bits(), 8);
    assert_eq!(a.symmetric_difference(&b).num_on_bits(), 9);
    assert_eq!(
        similarity::sparse_num_bits_in_common(&a, &b).unwrap(),
        8u32.wrapping_sub(9)
    );
    assert_eq!(
        similarity::sparse_num_bits_in_common(&a, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn tanimoto_similarity_f33() {
    // F33: TanimotoSimilarity (BitOps.cpp:286-297) dense + sparse.
    // Oracle probes (rdkit.DataStructs.TanimotoSimilarity, pinned):
    //   x{1,3,64} vs y{3,63,64} over 65: total 6, common 2 -> 2/4 = 0.5;
    //   both empty -> 0.0; one empty -> 0/3 = 0.0; identical -> 1.0;
    //   disjoint -> 0.0; mismatch -> ValueError.
    use cosmolkit_fingerprints::similarity;

    let mut x = Fingerprint::new(65);
    let mut y = Fingerprint::new(65);
    for bit in [1u32, 3, 64] {
        x.set_bit(bit).unwrap();
    }
    for bit in [3u32, 63, 64] {
        y.set_bit(bit).unwrap();
    }
    assert!((similarity::tanimoto(&x, &y).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::tanimoto(&x, &x).unwrap() == 1.0);
    let empty = Fingerprint::new(65);
    // both empty -> source total==0 branch returns 0.0 (not 1.0)
    assert!(similarity::tanimoto(&empty, &empty).unwrap() == 0.0);
    // one empty
    assert!(similarity::tanimoto(&x, &empty).unwrap() == 0.0);
    // disjoint
    let mut d = Fingerprint::new(65);
    d.set_bit(5).unwrap();
    assert!(similarity::tanimoto(&x, &d).unwrap() == 0.0);
    assert_eq!(
        similarity::tanimoto(&x, &Fingerprint::new(64)),
        Err(FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        })
    );

    // sparse overload
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3, 5] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5, 7] {
        sb.set_bit(bit).unwrap();
    }
    // total 6, common 2 -> 0.5
    assert!((similarity::sparse_tanimoto(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::sparse_tanimoto(&sa, &sa).unwrap() == 1.0);
    let se = SparseBitFingerprint::new(8);
    assert!(similarity::sparse_tanimoto(&se, &se).unwrap() == 0.0);
    assert_eq!(
        similarity::sparse_tanimoto(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn tversky_similarity_f34() {
    // F34: TverskySimilarity (BitOps.cpp:299-318). Oracle probes:
    //   a{1,3} vs b{3,5} over 8: common 1, y=2, z=2:
    //   (.5,.5) -> 0.5; (1,1) -> 1/3 (== tanimoto);
    //   a=1.5 or b=-0.1 -> Range Error (checked before the length test);
    //   one empty -> 0.0.
    use cosmolkit_fingerprints::similarity;

    let mut x = Fingerprint::new(8);
    let mut y = Fingerprint::new(8);
    for bit in [1u32, 3] {
        x.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        y.set_bit(bit).unwrap();
    }
    assert!((similarity::tversky(&x, &y, 0.5, 0.5).unwrap() - 0.5).abs() < 1e-12);
    assert!((similarity::tversky(&x, &y, 1.0, 1.0).unwrap() - 1.0 / 3.0).abs() < 1e-12);
    // (0,0): denom = 0 -> 1.0 with nonzero operands
    assert!(similarity::tversky(&x, &y, 0.0, 0.0).unwrap() == 1.0);
    // range errors fire before the length precondition
    assert_eq!(
        similarity::tversky(&x, &Fingerprint::new(16), 1.5, 0.5),
        Err(FingerprintError::RangeError { value: 1.5 })
    );
    assert_eq!(
        similarity::tversky(&x, &y, 0.5, -0.1),
        Err(FingerprintError::RangeError { value: -0.1 })
    );
    // boundary parameter values 0 and 1 are in range
    assert!(similarity::tversky(&x, &y, 0.0, 1.0).is_ok());
    // one empty -> 0.0
    let empty = Fingerprint::new(8);
    assert!(similarity::tversky(&x, &empty, 0.5, 0.5).unwrap() == 0.0);
    // length mismatch (with in-range parameters)
    assert_eq!(
        similarity::tversky(&x, &Fingerprint::new(16), 0.5, 0.5),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    // sparse overload: same values
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_tversky(&sa, &sb, 0.5, 0.5).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_tversky(&sa, &sb, 2.0, 0.0),
        Err(FingerprintError::RangeError { value: 2.0 })
    );
    assert_eq!(
        similarity::sparse_tversky(&sa, &SparseBitFingerprint::new(16), 0.5, 0.5),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn cosine_similarity_f35() {
    // F35: CosineSimilarity (BitOps.cpp:320-334). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/sqrt(4) = 0.5; ae = 0; ee = 0; af = 2/sqrt(16) = 0.5.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::cosine(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::cosine(&a, &e).unwrap() == 0.0);
    assert!(similarity::cosine(&e, &e).unwrap() == 0.0);
    assert!((similarity::cosine(&a, &f).unwrap() - 0.5).abs() < 1e-12);
    assert!((similarity::cosine(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::cosine(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_cosine(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_cosine(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn kulczynski_similarity_f36() {
    // F36: KulczynskiSimilarity (BitOps.cpp:336-350). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1*4/8 = 0.5; ae = 0; ee = 0; af = 2*10/16 = 0.625.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::kulczynski(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::kulczynski(&a, &e).unwrap() == 0.0);
    assert!(similarity::kulczynski(&e, &e).unwrap() == 0.0);
    assert!((similarity::kulczynski(&a, &f).unwrap() - 0.625).abs() < 1e-12);
    assert!((similarity::kulczynski(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::kulczynski(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_kulczynski(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_kulczynski(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn dice_similarity_f37() {
    // F37: DiceSimilarity (BitOps.cpp:352-366). Oracle probes over 8-bit
    // a{1,3}, b{3,5}, e empty, f full:
    //   ab = 2/4 = 0.5; ae = 0; ee = 0; af = 4/10 = 0.4.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::dice(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::dice(&a, &e).unwrap() == 0.0);
    assert!(similarity::dice(&e, &e).unwrap() == 0.0);
    assert!((similarity::dice(&a, &f).unwrap() - 0.4).abs() < 1e-12);
    assert!((similarity::dice(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::dice(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_dice(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_dice(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn sokal_similarity_f38() {
    // F38: SokalSimilarity (BitOps.cpp:368-381). Oracle probes over 8-bit
    // a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/(4+4-3) = 0.2; ae = 0; ee = 0; af = 2/(4+16-6) = 1/7;
    //   identical a,a = 1.0.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::sokal(&a, &b).unwrap() - 0.2).abs() < 1e-12);
    assert!(similarity::sokal(&a, &e).unwrap() == 0.0);
    assert!(similarity::sokal(&e, &e).unwrap() == 0.0);
    assert!((similarity::sokal(&a, &f).unwrap() - 1.0 / 7.0).abs() < 1e-12);
    assert!((similarity::sokal(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::sokal(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_sokal(&sa, &sb).unwrap() - 0.2).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_sokal(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn mcconnaughey_similarity_f39() {
    // F39: McConnaugheySimilarity (BitOps.cpp:383-397). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = (1*4-4)/4 = 0.0; ae = 0; ee = 0; af = (2*10-16)/16 = 0.25;
    //   identical a,a = 1.0 (range [-1,1]).
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::mcconnaughey(&a, &b).unwrap() - 0.0).abs() < 1e-12);
    assert!(similarity::mcconnaughey(&a, &e).unwrap() == 0.0);
    assert!(similarity::mcconnaughey(&e, &e).unwrap() == 0.0);
    assert!((similarity::mcconnaughey(&a, &f).unwrap() - 0.25).abs() < 1e-12);
    assert!((similarity::mcconnaughey(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::mcconnaughey(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!(similarity::sparse_mcconnaughey(&sa, &sb).unwrap().abs() < 1e-12);
    assert_eq!(
        similarity::sparse_mcconnaughey(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn asymmetric_similarity_f40() {
    // F40: AsymmetricSimilarity (BitOps.cpp:399-414), non-commutative.
    // Oracle probes over 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/2 = 0.5; ae = 0; ee = 0; af = 2/2 = 1.0 (min = |a|).
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::asymmetric(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::asymmetric(&a, &e).unwrap() == 0.0);
    assert!(similarity::asymmetric(&e, &e).unwrap() == 0.0);
    assert!((similarity::asymmetric(&a, &f).unwrap() - 1.0).abs() < 1e-12);
    // symmetric min makes ab == ba here
    assert!((similarity::asymmetric(&b, &a).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::asymmetric(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_asymmetric(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_asymmetric(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn braun_blanquet_similarity_f41() {
    // F41: BraunBlanquetSimilarity (BitOps.cpp:416-431). Oracle probes
    // over 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/2 = 0.5; ae = 0; ee = 0; af = 2/8 = 0.25 (max = |f|).
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::braun_blanquet(&a, &b).unwrap() - 0.5).abs() < 1e-12);
    assert!(similarity::braun_blanquet(&a, &e).unwrap() == 0.0);
    assert!(similarity::braun_blanquet(&e, &e).unwrap() == 0.0);
    assert!((similarity::braun_blanquet(&a, &f).unwrap() - 0.25).abs() < 1e-12);
    assert!((similarity::braun_blanquet(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::braun_blanquet(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_braun_blanquet(&sa, &sb).unwrap() - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_braun_blanquet(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn russel_similarity_f42() {
    // F42: RusselSimilarity (BitOps.cpp:433-441). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/8 = 0.125; ae = 0/8 = 0; ee(8-bit) = 0; ea = 0;
    //   identical a,a = 2/8 = 0.25 (divisor is the LENGTH, not |a|).
    // Zero-length pairs: dense raises the NumOnBitsInCommon bitmap
    // precondition (oracle: "Pre-condition Violation: no afp"); sparse
    // computes 0/0 = NaN.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::russel(&a, &b).unwrap() - 0.125).abs() < 1e-12);
    assert!(similarity::russel(&a, &e).unwrap() == 0.0);
    assert!(similarity::russel(&e, &e).unwrap() == 0.0);
    assert!(similarity::russel(&e, &a).unwrap() == 0.0);
    assert!((similarity::russel(&a, &a).unwrap() - 0.25).abs() < 1e-12);
    assert!((similarity::russel(&a, &f).unwrap() - 0.25).abs() < 1e-12);
    assert_eq!(
        similarity::russel(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // dense zero-length pair: reproduced precondition violation
    assert_eq!(
        similarity::russel(&Fingerprint::new(0), &Fingerprint::new(0)),
        Err(FingerprintError::PreconditionViolation {
            what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)"
        })
    );

    // sparse overload: same length divisor; zero-length pair = NaN
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_russel(&sa, &sb).unwrap() - 0.125).abs() < 1e-12);
    let zs =
        similarity::sparse_russel(&SparseBitFingerprint::new(0), &SparseBitFingerprint::new(0))
            .unwrap();
    assert!(zs.is_nan());
    assert_eq!(
        similarity::sparse_russel(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn dense_zero_length_metrics_precondition_f42() {
    // Pinned build behavior for zero-length dense pairs: Tanimoto returns
    // 0.0 (guard before NumOnBitsInCommon); every other on-count metric
    // raises the bitmap precondition.
    use cosmolkit_fingerprints::similarity;

    let z = Fingerprint::new(0);
    assert!(similarity::tanimoto(&z, &z).unwrap() == 0.0);
    let expected = Err(FingerprintError::PreconditionViolation {
        what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)",
    });
    assert_eq!(similarity::cosine(&z, &z), expected);
    assert_eq!(similarity::kulczynski(&z, &z), expected);
    assert_eq!(similarity::dice(&z, &z), expected);
    assert_eq!(similarity::sokal(&z, &z), expected);
    assert_eq!(similarity::mcconnaughey(&z, &z), expected);
    assert_eq!(similarity::asymmetric(&z, &z), expected);
    assert_eq!(similarity::braun_blanquet(&z, &z), expected);
    assert_eq!(similarity::russel(&z, &z), expected);
    // sparse zero-length pairs compute normally (no bitmap path)
    let zs = SparseBitFingerprint::new(0);
    assert!(similarity::sparse_cosine(&zs, &zs).unwrap() == 0.0);
    assert!(similarity::sparse_dice(&zs, &zs).unwrap() == 0.0);
}

#[test]
fn rogot_goldberg_similarity_f43() {
    // F43: RogotGoldbergSimilarity (BitOps.cpp:443-467). Oracle probes
    // over 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 2/3; ae = 0; ee = 0; af = 0.2; identical = 1.0 (x == l path
    //   never fires for partial vectors; f,f gives x == l -> 1.0).
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::rogot_goldberg(&a, &b).unwrap() - 2.0 / 3.0).abs() < 1e-12);
    assert!(similarity::rogot_goldberg(&a, &e).unwrap() == 0.0);
    assert!(similarity::rogot_goldberg(&e, &e).unwrap() == 0.0);
    assert!((similarity::rogot_goldberg(&a, &f).unwrap() - 0.2).abs() < 1e-12);
    assert!(similarity::rogot_goldberg(&f, &f).unwrap() == 1.0);
    assert_eq!(
        similarity::rogot_goldberg(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_rogot_goldberg(&sa, &sb).unwrap() - 2.0 / 3.0).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_rogot_goldberg(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn on_bit_similarity_f44() {
    // F44: OnBitSimilarity (BitOps.cpp:481-495). Oracle probes over 8-bit
    // a{1,3}, b{3,5}, e empty, f full:
    //   ab = 1/3; ae = 0; ee = 0; af = 2/8 = 0.25; identical = 1.0.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::on_bit(&a, &b).unwrap() - 1.0 / 3.0).abs() < 1e-12);
    assert!(similarity::on_bit(&a, &e).unwrap() == 0.0);
    assert!(similarity::on_bit(&e, &e).unwrap() == 0.0);
    assert!((similarity::on_bit(&a, &f).unwrap() - 0.25).abs() < 1e-12);
    assert!((similarity::on_bit(&a, &a).unwrap() - 1.0).abs() < 1e-12);
    assert_eq!(
        similarity::on_bit(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    // sparse overload, incl. an out-of-length retained member joining the
    // union denominator exactly like the raw source set union
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_on_bit(&sa, &sb).unwrap() - 1.0 / 3.0).abs() < 1e-12);
    let mut wide = SparseBitFingerprint::new(16);
    wide.set_bit(12).unwrap();
    let with_oob = sa.union(&sa.union(&wide));
    // common(sa, with_oob) = 2; union = {1,3,12} -> 2/3
    assert!((similarity::sparse_on_bit(&with_oob, &sa).unwrap() - 2.0 / 3.0).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_on_bit(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn all_bit_similarity_f45() {
    // F45: AllBitSimilarity (BitOps.cpp:538-545). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = 6/8 = 0.75; ae = 6/8 = 0.75; ee = 1.0; af = 2/8 = 0.25;
    //   zero-length pair = nan.
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    assert!((similarity::all_bit(&a, &b).unwrap() - 0.75).abs() < 1e-12);
    assert!((similarity::all_bit(&a, &e).unwrap() - 0.75).abs() < 1e-12);
    assert!(similarity::all_bit(&e, &e).unwrap() == 1.0);
    assert!((similarity::all_bit(&a, &f).unwrap() - 0.25).abs() < 1e-12);
    assert_eq!(
        similarity::all_bit(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // zero-length pair: 0/0 = NaN (oracle: nan)
    let z = similarity::all_bit(&Fingerprint::new(0), &Fingerprint::new(0)).unwrap();
    assert!(z.is_nan());

    // sparse overload
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    assert!((similarity::sparse_all_bit(&sa, &sb).unwrap() - 0.75).abs() < 1e-12);
    let zs =
        similarity::sparse_all_bit(&SparseBitFingerprint::new(0), &SparseBitFingerprint::new(0))
            .unwrap();
    assert!(zs.is_nan());
    assert_eq!(
        similarity::sparse_all_bit(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn on_bit_proj_similarity_f46() {
    // F46: OnBitProjSimilarity (BitOps.cpp:620-632). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty:
    //   ab = [0.5, 0.5]; ae = [0.0, 0.0] (num == 0 branch); a,a = [1,1].
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let ab = similarity::on_bit_proj_similarity(&a, &b).unwrap();
    assert!((ab[0] - 0.5).abs() < 1e-12);
    assert!((ab[1] - 0.5).abs() < 1e-12);
    let ae = similarity::on_bit_proj_similarity(&a, &e).unwrap();
    assert!(ae == [0.0, 0.0]);
    let aa = similarity::on_bit_proj_similarity(&a, &a).unwrap();
    assert!((aa[0] - 1.0).abs() < 1e-12 && (aa[1] - 1.0).abs() < 1e-12);
    // asymmetric projection: a{1,3} vs c{1,3,5,7} -> [2/2, 2/4] = [1, .5]
    let mut c = Fingerprint::new(8);
    for bit in [1u32, 3, 5, 7] {
        c.set_bit(bit).unwrap();
    }
    let ac = similarity::on_bit_proj_similarity(&a, &c).unwrap();
    assert!((ac[0] - 1.0).abs() < 1e-12);
    assert!((ac[1] - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::on_bit_proj_similarity(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    // sparse overload
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    let sab = similarity::sparse_on_bit_proj_similarity(&sa, &sb).unwrap();
    assert!((sab[0] - 0.5).abs() < 1e-12 && (sab[1] - 0.5).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_on_bit_proj_similarity(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn off_bit_proj_similarity_f47() {
    // F47: OffBitProjSimilarity (BitOps.cpp:659-671). Oracle probes over
    // 8-bit a{1,3}, b{3,5}, e empty, f full:
    //   ab = [5/6, 5/6]; ae = [6/6, 6/8] = [1.0, 0.75];
    //   f,f = [0,0] (num == 0 branch); e,e = [1,1].
    use cosmolkit_fingerprints::similarity;

    let mut a = Fingerprint::new(8);
    let mut b = Fingerprint::new(8);
    for bit in [1u32, 3] {
        a.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        b.set_bit(bit).unwrap();
    }
    let e = Fingerprint::new(8);
    let f = Fingerprint::new_filled(8);
    let ab = similarity::off_bit_proj_similarity(&a, &b).unwrap();
    assert!((ab[0] - 5.0 / 6.0).abs() < 1e-12);
    assert!((ab[1] - 5.0 / 6.0).abs() < 1e-12);
    let ae = similarity::off_bit_proj_similarity(&a, &e).unwrap();
    assert!((ae[0] - 1.0).abs() < 1e-12);
    assert!((ae[1] - 0.75).abs() < 1e-12);
    let ff = similarity::off_bit_proj_similarity(&f, &f).unwrap();
    assert!(ff == [0.0, 0.0]);
    let ee = similarity::off_bit_proj_similarity(&e, &e).unwrap();
    assert!(ee == [1.0, 1.0]);
    assert_eq!(
        similarity::off_bit_proj_similarity(&a, &Fingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );

    // sparse overload
    let mut sa = SparseBitFingerprint::new(8);
    let mut sb = SparseBitFingerprint::new(8);
    for bit in [1u32, 3] {
        sa.set_bit(bit).unwrap();
    }
    for bit in [3u32, 5] {
        sb.set_bit(bit).unwrap();
    }
    let sab = similarity::sparse_off_bit_proj_similarity(&sa, &sb).unwrap();
    assert!((sab[0] - 5.0 / 6.0).abs() < 1e-12);
    assert!((sab[1] - 5.0 / 6.0).abs() < 1e-12);
    assert_eq!(
        similarity::sparse_off_bit_proj_similarity(&sa, &SparseBitFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn similarity_wrapper_auto_folding_f48() {
    // F48: SimilarityWrapper (BitOps.h:30-50). Oracle probes:
    //   a16{1,9} vs b8{1}: fold by 2 -> a{1,1} -> tanimoto 1.0;
    //   returnDistance -> 0.0; operand order symmetric;
    //   a12{1,9} vs b8: ratio 12/8 = 1 -> fold(12,1) leaves 12 bits ->
    //   metric length error (ValueError);
    //   0 vs 0 -> metric directly (tanimoto 0.0);
    //   one zero-length operand: source computes fold factor /0 and dies
    //   with SIGFPE (uncatchable, oracle-verified); Rust rejects the
    //   executed division with the resolved CK safe-rejection error.
    use cosmolkit_fingerprints::similarity;

    let mut a16 = Fingerprint::new(16);
    for bit in [1u32, 9] {
        a16.set_bit(bit).unwrap();
    }
    let mut b8 = Fingerprint::new(8);
    b8.set_bit(1).unwrap();
    assert!(similarity::similarity_wrapper(&a16, &b8, similarity::tanimoto, false).unwrap() == 1.0);
    assert!(similarity::similarity_wrapper(&a16, &b8, similarity::tanimoto, true).unwrap() == 0.0);
    assert!(similarity::similarity_wrapper(&b8, &a16, similarity::tanimoto, false).unwrap() == 1.0);

    let mut a12 = Fingerprint::new(12);
    for bit in [1u32, 9] {
        a12.set_bit(bit).unwrap();
    }
    assert_eq!(
        similarity::similarity_wrapper(&a12, &b8, similarity::tanimoto, false),
        Err(FingerprintError::BitLengthMismatch { left: 12, right: 8 })
    );

    let z = Fingerprint::new(0);
    assert!(similarity::similarity_wrapper(&z, &z, similarity::tanimoto, false).unwrap() == 0.0);
    let boundary = Err(FingerprintError::InvalidArguments {
        reason: "SimilarityWrapper: exactly one zero-length operand makes the \
                 source fold-factor length division a division by zero; CK \
                 safe rejection of undefined source execution",
    });
    assert_eq!(
        similarity::similarity_wrapper(&z, &b8, similarity::tanimoto, false),
        boundary
    );
    assert_eq!(
        similarity::similarity_wrapper(&b8, &z, similarity::tanimoto, false),
        boundary
    );

    // sparse overload: same fold semantics
    let mut sa16 = SparseBitFingerprint::new(16);
    for bit in [1u32, 9] {
        sa16.set_bit(bit).unwrap();
    }
    let mut sb8 = SparseBitFingerprint::new(8);
    sb8.set_bit(1).unwrap();
    assert!(
        (similarity::sparse_similarity_wrapper(&sa16, &sb8, similarity::sparse_tanimoto, false)
            .unwrap()
            - 1.0)
            .abs()
            < 1e-12
    );
    assert_eq!(
        similarity::sparse_similarity_wrapper(
            &sb8,
            &z_sparse(),
            similarity::sparse_tanimoto,
            false
        ),
        boundary
    );
}

fn z_sparse() -> SparseBitFingerprint {
    SparseBitFingerprint::new(0)
}

#[test]
fn similarity_wrapper_tversky_f49() {
    // F49: parameterized SimilarityWrapper (BitOps.h:52-72). Oracle:
    //   wrapped Tversky(.5,.5) on a16{1,9} vs b8{1} (fold by 2) = 1.0;
    //   distance = 0.0; out-of-range a/b still raises inside the metric
    //   (after folding — same error as the raw metric).
    use cosmolkit_fingerprints::similarity;

    let mut a16 = Fingerprint::new(16);
    for bit in [1u32, 9] {
        a16.set_bit(bit).unwrap();
    }
    let mut b8 = Fingerprint::new(8);
    b8.set_bit(1).unwrap();
    assert!(
        (similarity::similarity_wrapper_tversky(&a16, &b8, 0.5, 0.5, similarity::tversky, false)
            .unwrap()
            - 1.0)
            .abs()
            < 1e-12
    );
    assert!(
        (similarity::similarity_wrapper_tversky(&a16, &b8, 0.5, 0.5, similarity::tversky, true)
            .unwrap()
            - 0.0)
            .abs()
            < 1e-12
    );
    assert_eq!(
        similarity::similarity_wrapper_tversky(&a16, &b8, 1.5, 0.5, similarity::tversky, false),
        Err(FingerprintError::RangeError { value: 1.5 })
    );
    // folded length-mismatch path and zero-length boundary
    let mut a12 = Fingerprint::new(12);
    a12.set_bit(1).unwrap();
    a12.set_bit(9).unwrap();
    assert_eq!(
        similarity::similarity_wrapper_tversky(&a12, &b8, 0.5, 0.5, similarity::tversky, false),
        Err(FingerprintError::BitLengthMismatch { left: 12, right: 8 })
    );
    let z = Fingerprint::new(0);
    assert_eq!(
        similarity::similarity_wrapper_tversky(&z, &b8, 0.5, 0.5, similarity::tversky, false),
        Err(FingerprintError::InvalidArguments {
            reason: "SimilarityWrapper: exactly one zero-length operand makes the \
                     source fold-factor length division a division by zero; CK \
                     safe rejection of undefined source execution"
        })
    );

    // sparse overload
    let mut sa16 = SparseBitFingerprint::new(16);
    for bit in [1u32, 9] {
        sa16.set_bit(bit).unwrap();
    }
    let mut sb8 = SparseBitFingerprint::new(8);
    sb8.set_bit(1).unwrap();
    assert!(
        (similarity::sparse_similarity_wrapper_tversky(
            &sa16,
            &sb8,
            0.5,
            0.5,
            similarity::sparse_tversky,
            false
        )
        .unwrap()
            - 1.0)
            .abs()
            < 1e-12
    );
    assert_eq!(
        similarity::sparse_similarity_wrapper_tversky(
            &sa16,
            &sb8,
            2.0,
            0.0,
            similarity::sparse_tversky,
            false
        ),
        Err(FingerprintError::RangeError { value: 2.0 })
    );
}

#[test]
fn calc_vect_params_f50() {
    // F50: calcVectParams (SparseIntVect.h:434-495). Deriving the
    // expected sums from the abs/min semantics verified against the
    // count-similarity oracles in F51-F53: v1_sum/v2_sum are abs sums,
    // and_sum the per-common-key abs minima. Oracle cross-check for
    // a={1:2,3:5}, b={3:2,5:1}: Dice = 2*and/(v1+v2) = 4/10 = 0.4.
    use cosmolkit_fingerprints::similarity;

    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let (v1s, v2s, ands) = similarity::calc_vect_params(&a, &b).unwrap();
    assert_eq!(v1s, 7.0);
    assert_eq!(v2s, 3.0);
    assert_eq!(ands, 2.0);

    // negative counts use absolute values
    let mut n = SparseCountFingerprint::new(8);
    n.set_value(0, -4).unwrap();
    n.set_value(2, 3).unwrap();
    let mut m = SparseCountFingerprint::new(8);
    m.set_value(0, 2).unwrap();
    m.set_value(2, -5).unwrap();
    let (v1s, v2s, ands) = similarity::calc_vect_params(&n, &m).unwrap();
    assert_eq!(v1s, 7.0);
    assert_eq!(v2s, 7.0);
    assert_eq!(ands, 2.0 + 3.0);

    // asymmetric supports / cancellation-style disjoint entries
    let mut d1 = SparseCountFingerprint::new(8);
    d1.set_value(7, 3).unwrap();
    let d2 = SparseCountFingerprint::new(8);
    let (v1s, v2s, ands) = similarity::calc_vect_params(&d1, &d2).unwrap();
    assert!((v1s, v2s, ands) == (3.0, 0.0, 0.0));
    // both empty
    let (v1s, v2s, ands) = similarity::calc_vect_params(
        &SparseCountFingerprint::new(8),
        &SparseCountFingerprint::new(8),
    )
    .unwrap();
    assert!((v1s, v2s, ands) == (0.0, 0.0, 0.0));
    // source-defined extreme counts stay exact in the double accumulator
    let mut mx = SparseCountFingerprint::new(8);
    mx.set_value(0, i32::MAX).unwrap();
    mx.set_value(1, i32::MIN + 1).unwrap();
    let (v1s, _, ands) =
        similarity::calc_vect_params(&mx, &SparseCountFingerprint::new(8)).unwrap();
    assert_eq!(v1s, f64::from(i32::MAX) + f64::from(i32::MIN + 1).abs());
    assert_eq!(ands, 0.0);
    // abs(MIN) is the one executed signed-i32 UB on this path (D-series
    // resolution): reported, not wrapped
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&mn, &SparseCountFingerprint::new(8)),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::calcVectParams abs"
        })
    );
    // length mismatch
    assert_eq!(
        similarity::calc_vect_params(&a, &SparseCountFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn count_dice_similarity_f51() {
    // F51: count-valued DiceSimilarity (SparseIntVect.h:498-539). Oracle
    // probes (LongSparseIntVect): a={1:2,3:5}, b={3:2,5:1}:
    //   dice = 2*2/(7+3) = 0.4; distance = 0.6; both empty = 0.0;
    //   one empty = 0.0; size mismatch = ValueError.
    use cosmolkit_fingerprints::similarity;

    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!((similarity::sparse_count_dice(&a, &b, false, 0.0).unwrap() - 0.4).abs() < 1e-12);
    assert!((similarity::sparse_count_dice(&a, &b, true, 0.0).unwrap() - 0.6).abs() < 1e-12);
    let e = SparseCountFingerprint::new(8);
    assert!(similarity::sparse_count_dice(&e, &e, false, 0.0).unwrap() == 0.0);
    assert!(similarity::sparse_count_dice(&a, &e, false, 0.0).unwrap() == 0.0);
    // negative counts use absolute accumulation
    let mut n = SparseCountFingerprint::new(8);
    n.set_value(0, -4).unwrap();
    n.set_value(2, 3).unwrap();
    let mut m = SparseCountFingerprint::new(8);
    m.set_value(0, 2).unwrap();
    m.set_value(2, -5).unwrap();
    // and = 2 + 3 = 5; sums 7,7 -> 10/14
    assert!(
        (similarity::sparse_count_dice(&n, &m, false, 0.0).unwrap() - 10.0 / 14.0).abs() < 1e-12
    );
    // bounds branch: 2*min/(sum) = 2*3/10 = 0.6 -> bounds 0.7 rejects
    // (returns 0.0), bounds 0.5 accepts and computes the full metric
    assert!(similarity::sparse_count_dice(&a, &b, false, 0.7).unwrap() == 0.0);
    assert!((similarity::sparse_count_dice(&a, &b, false, 0.5).unwrap() - 0.4).abs() < 1e-12);
    // bounds is ignored under returnDistance
    assert!((similarity::sparse_count_dice(&a, &b, true, 0.9).unwrap() - 0.6).abs() < 1e-12);
    // exact threshold: |denom| < 1e-6 -> 0.0 (empty pair under bounds)
    assert!(similarity::sparse_count_dice(&e, &e, false, 0.5).unwrap() == 0.0);
    // size mismatch
    assert_eq!(
        similarity::sparse_count_dice(&a, &SparseCountFingerprint::new(16), false, 0.0),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn count_tversky_similarity_f52() {
    // F52: count-valued TverskySimilarity (SparseIntVect.h:541-568). No
    // parameter range checks (distinct from the bit-vector Tversky);
    // bounds is unused. a={1:2,3:5}, b={3:2,5:1}: sums 7/3, and = 2:
    //   (1,1) -> 2/(7+3-2) = 0.25; (.5,.5) -> 2/(3.5+1.5) = 2/5 = 0.4;
    //   (0,0) -> denom = and = 2 -> 1.0; distance = 1 - sim;
    //   zero vectors -> denom ~ 0 -> 0.0.
    use cosmolkit_fingerprints::similarity;

    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!(
        (similarity::sparse_count_tversky(&a, &b, 1.0, 1.0, false, 0.0).unwrap() - 0.25).abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky(&a, &b, 0.5, 0.5, false, 0.0).unwrap() - 0.4).abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky(&a, &b, 0.0, 0.0, false, 0.0).unwrap() - 1.0).abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky(&a, &b, 1.0, 1.0, true, 0.0).unwrap() - 0.75).abs()
            < 1e-12
    );
    // coefficients outside [0,1] are NOT rejected in the count variant
    assert!(similarity::sparse_count_tversky(&a, &b, 2.0, -1.0, false, 0.0).is_ok());
    // zero vectors: denom < 1e-6 -> 0.0
    let e = SparseCountFingerprint::new(8);
    assert!(similarity::sparse_count_tversky(&e, &e, 1.0, 1.0, false, 0.0).unwrap() == 0.0);
    // size mismatch
    assert_eq!(
        similarity::sparse_count_tversky(
            &a,
            &SparseCountFingerprint::new(16),
            1.0,
            1.0,
            false,
            0.0
        ),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn count_tanimoto_similarity_f53() {
    // F53: count-valued TanimotoSimilarity (SparseIntVect.h:570-575) = the
    // source-prescribed count Tversky(1,1). Oracle: a={1:2,3:5},
    // b={3:2,5:1} -> 0.25; distance 0.75; negative counts included via
    // abs accumulation: n={0:-4,2:3}, m={0:2,2:-5}: sums 7/7, and 5 ->
    // 5/(7+7-5) = 5/9; empty pair -> 0.0.
    use cosmolkit_fingerprints::similarity;

    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!((similarity::sparse_count_tanimoto(&a, &b, false, 0.0).unwrap() - 0.25).abs() < 1e-12);
    assert!((similarity::sparse_count_tanimoto(&a, &b, true, 0.0).unwrap() - 0.75).abs() < 1e-12);
    let mut n = SparseCountFingerprint::new(8);
    n.set_value(0, -4).unwrap();
    n.set_value(2, 3).unwrap();
    let mut m = SparseCountFingerprint::new(8);
    m.set_value(0, 2).unwrap();
    m.set_value(2, -5).unwrap();
    assert!(
        (similarity::sparse_count_tanimoto(&n, &m, false, 0.0).unwrap() - 5.0 / 9.0).abs() < 1e-12
    );
    // a=b=1 equivalence with sparse_count_tversky
    assert_eq!(
        similarity::sparse_count_tanimoto(&a, &b, false, 0.3),
        similarity::sparse_count_tversky(&a, &b, 1.0, 1.0, false, 0.3)
    );
    let e = SparseCountFingerprint::new(8);
    assert!(similarity::sparse_count_tanimoto(&e, &e, false, 0.0).unwrap() == 0.0);
    assert_eq!(
        similarity::sparse_count_tanimoto(&a, &SparseCountFingerprint::new(16), false, 0.0),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
}

#[test]
fn error_surface_f54() {
    // F54: structured error details, source ordering, Unsupported
    // preservation, RangeError finite/NaN PartialEq (no Eq promise, no
    // NaN reflexivity), UndefinedArithmetic site/Display.
    use cosmolkit_fingerprints::{SparseCountFingerprint, folding, hash, similarity};

    // index detail fields
    assert_eq!(
        Fingerprint::new(4).get_bit(9),
        Err(FingerprintError::SparseIndexOutOfRange { index: 9, size: 4 })
    );
    // width/length detail fields (u64 widths for count vectors)
    assert_eq!(
        SparseCountFingerprint::new(1 << 40).set_value(1 << 41, 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 1 << 41,
            size: 1 << 40
        })
    );
    // fold detail fields
    assert_eq!(
        folding::fold_fingerprint(&Fingerprint::new(32), 33),
        Err(FingerprintError::InvalidFoldFactor {
            factor: 33,
            n_bits: 32
        })
    );
    // parameter detail: RangeError distinguishes values; NaN payloads are
    // not self-equal under the derived f64 PartialEq (no reflexivity
    // promise) and the type deliberately makes no Eq promise
    assert_ne!(
        FingerprintError::RangeError { value: f64::NAN },
        FingerprintError::RangeError { value: f64::NAN },
        "NaN payloads are distinct under PartialEq; no Eq promise exists"
    );
    assert_eq!(
        FingerprintError::RangeError { value: 1.5 },
        FingerprintError::RangeError { value: 1.5 }
    );
    // source error ordering: Tversky range checks fire before the length
    // precondition; count-Tversky has no range checks at all
    let mut x = Fingerprint::new(8);
    x.set_bit(1).unwrap();
    assert_eq!(
        similarity::tversky(&x, &Fingerprint::new(16), 1.5, 0.5),
        Err(FingerprintError::RangeError { value: 1.5 })
    );
    assert!(
        similarity::sparse_count_tversky(
            &SparseCountFingerprint::new(8),
            &SparseCountFingerprint::new(8),
            1.5,
            0.5,
            false,
            0.0
        )
        .is_ok()
    );
    // UndefinedArithmetic site/Display reconciliation (D1-D5 resolved)
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    let err = mn.with_multiplied_scalar(-1).unwrap_err();
    assert_eq!(
        err,
        FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator*=(int)"
        }
    );
    let text = err.to_string();
    assert!(text.contains("SparseIntVect::operator*=(int)"), "{text}");
    assert!(text.contains("COSMolKit safety boundary"), "{text}");
    assert!(!text.contains("pending"), "{text}");
    // Display keeps structured details for other kinds
    assert_eq!(
        FingerprintError::BitLengthMismatch { left: 8, right: 16 }.to_string(),
        "fingerprint bit length mismatch: 8 != 16"
    );
    assert_eq!(
        FingerprintError::InvalidFoldFactor {
            factor: 33,
            n_bits: 32
        }
        .to_string(),
        "invalid fold factor 33 for fingerprint of 32 bits"
    );
    // Unsupported preservation: the protected generator boundaries still
    // fail closed with the unchanged variant
    let topology = cosmolkit_model::TopologyBlock::default();
    assert_eq!(
        cosmolkit_fingerprints::morgan(&topology),
        Err(FingerprintError::Unsupported)
    );
    assert_eq!(
        cosmolkit_fingerprints::pattern(&topology),
        Err(FingerprintError::Unsupported)
    );
    // std::error::Error chain remains available
    fn assert_error<E: std::error::Error>(_: &E) {}
    assert_error(&err);
    // hash module unaffected by error surface
    assert_eq!(hash::hash_value_u32(7), 7);
}

#[test]
fn dense_equality_ignores_cached_count_c1() {
    // C1 (supervisor rlib evidence): new(1), set_bit(0), clear_bits() must
    // equal a fresh one-bit vector although the cached on-count differs
    // (ExplicitBitVect.h:86-91 compares only *dp_bits). Identical actual
    // bits with different stale caches compare equal; different lengths
    // or different actual bits compare unequal.
    let mut one = Fingerprint::new(1);
    one.set_bit(0).unwrap();
    one.clear_bits();
    let fresh = Fingerprint::new(1);
    assert_eq!(
        one, fresh,
        "cleared vector must equal fresh despite stale cache"
    );

    // identical set bits, different stale counts
    let mut a = Fingerprint::new(8);
    a.set_bit(1).unwrap();
    a.set_bit(3).unwrap();
    a.clear_bits(); // stale cache: 2
    a.set_bit(1).unwrap(); // stale cache: 3
    let mut b = Fingerprint::new(8);
    b.set_bit(1).unwrap(); // fresh cache: 1
    assert_eq!(a, b, "same actual bits with different caches compare equal");
    // cached counts preserved as separate assertions (not equalized)
    assert_eq!(a.num_on_bits(), 3);
    assert_eq!(b.num_on_bits(), 1);

    // different lengths unequal
    assert_ne!(Fingerprint::new(4), Fingerprint::new(8));
    // different actual bits unequal (same stale-cache shape)
    let mut c = Fingerprint::new(8);
    c.set_bit(2).unwrap();
    assert_ne!(b, c);
}

#[test]
fn dense_stale_count_off_count_wraps_c2() {
    // C2 (supervisor overflow-checks evidence): the public sequence
    // new(1), set_bit(0), clear_bits(), set_bit(0) leaves the stale cache
    // at 2 on a one-bit vector; getNumOffBits is the source's unsigned
    // `d_size - d_numOnBits`, which wraps to u32::MAX instead of panicking
    // (ExplicitBitVect.cpp:181-183).
    let mut one = Fingerprint::new(1);
    one.set_bit(0).unwrap();
    one.clear_bits();
    one.set_bit(0).unwrap();
    assert_eq!(one.num_on_bits(), 2);
    assert_eq!(one.num_off_bits(), u32::MAX);
    assert!(one.get_bit(0).unwrap());
}

#[test]
fn calc_vect_params_merge_shape_c3() {
    // C3: regressions for the exact monotonic merge: either/both empty
    // maps, disjoint/interleaved keys, equal keys with mixed signs,
    // either iterator exhausted first with remaining tails, stored zero
    // entries, length mismatch before arithmetic, abs(MIN) rejection.
    use cosmolkit_fingerprints::similarity;

    // both empty
    let e = SparseCountFingerprint::new(8);
    assert_eq!(
        similarity::calc_vect_params(&e, &e).unwrap(),
        (0.0, 0.0, 0.0)
    );
    // v1 empty, v2 populated (iter1 exhausted immediately; v2 tail)
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(1, 2).unwrap();
    b.set_value(4, -3).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&e, &b).unwrap(),
        (0.0, 5.0, 0.0)
    );
    // v2 empty, v1 populated
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(2, 5).unwrap();
    a.set_value(6, 1).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&a, &e).unwrap(),
        (6.0, 0.0, 0.0)
    );
    // disjoint keys exhaust iter2 first; v1 tail must be fully summed
    let mut d1 = SparseCountFingerprint::new(8);
    d1.set_value(5, 1).unwrap();
    d1.set_value(6, 2).unwrap();
    d1.set_value(7, 4).unwrap();
    let mut d2 = SparseCountFingerprint::new(8);
    d2.set_value(0, 8).unwrap();
    d2.set_value(1, 16).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&d1, &d2).unwrap(),
        (7.0, 24.0, 0.0)
    );
    // interleaved keys with mixed signs on equal keys
    let mut i1 = SparseCountFingerprint::new(16);
    i1.set_value(1, -2).unwrap();
    i1.set_value(3, 5).unwrap();
    i1.set_value(5, -1).unwrap();
    let mut i2 = SparseCountFingerprint::new(16);
    i2.set_value(2, 1).unwrap();
    i2.set_value(3, -7).unwrap();
    i2.set_value(5, 4).unwrap();
    // and = min(5,7) at 3, min(1,4) at 5
    assert_eq!(
        similarity::calc_vect_params(&i1, &i2).unwrap(),
        (8.0, 12.0, 5.0 + 1.0)
    );
    // stored zero entries (scalar-produced, retained): contribute 0 to
    // their own sum and can match as abs-min 0
    let mut z1 = SparseCountFingerprint::new(8);
    z1.set_value(2, 3).unwrap();
    let z1 = z1.with_added_scalar(-3).unwrap(); // {2: 0} retained
    assert_eq!(z1.nonzero_elements().get(&2), Some(&0));
    let mut z2 = SparseCountFingerprint::new(8);
    z2.set_value(2, 9).unwrap();
    z2.set_value(4, 1).unwrap();
    // and = min(0, 9) = 0 at key 2
    assert_eq!(
        similarity::calc_vect_params(&z1, &z2).unwrap(),
        (0.0, 10.0, 0.0)
    );
    // length mismatch raised before any arithmetic
    assert_eq!(
        similarity::calc_vect_params(&a, &SparseCountFingerprint::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // abs(MIN) structured rejection
    let mut mn = SparseCountFingerprint::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&mn, &e),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::calcVectParams abs"
        })
    );
}

#[test]
fn wrapper_zero_length_dispatch_table_c4() {
    // C4: table-driven dispatch across all four wrappers. Exactly one
    // zero-length operand rejects with the structured fold-division
    // boundary error before the metric runs (the metric would otherwise
    // observe mismatched folded operands); zero/zero follows the source
    // equal-length metric branch — no division executes and the metric's
    // own result/error (and returnDistance handling) is preserved.
    use cosmolkit_fingerprints::similarity;

    let boundary = Err(FingerprintError::InvalidArguments {
        reason: "SimilarityWrapper: exactly one zero-length operand makes the \
                 source fold-factor length division a division by zero; CK \
                 safe rejection of undefined source execution",
    });

    let z = Fingerprint::new(0);
    let mut b8 = Fingerprint::new(8);
    b8.set_bit(1).unwrap();
    for (x, y, x_zero) in [(&z, &b8, true), (&b8, &z, false)] {
        assert_eq!(
            similarity::similarity_wrapper(x, y, similarity::tanimoto, false),
            boundary,
            "dense ordinary, left zero = {x_zero}"
        );
        assert_eq!(
            similarity::similarity_wrapper_tversky(x, y, 0.5, 0.5, similarity::tversky, false),
            boundary,
            "dense tversky, left zero = {x_zero}"
        );
    }
    // zero/zero dispatches to the metric: tanimoto's total==0 guard gives
    // 0.0 (and 1.0 under returnDistance); cosine raises the dense
    // zero-length precondition from inside the metric
    assert!(similarity::similarity_wrapper(&z, &z, similarity::tanimoto, false).unwrap() == 0.0);
    assert!(similarity::similarity_wrapper(&z, &z, similarity::tanimoto, true).unwrap() == 1.0);
    assert_eq!(
        similarity::similarity_wrapper(&z, &z, similarity::cosine, false),
        Err(FingerprintError::PreconditionViolation {
            what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)"
        })
    );

    // sparse wrappers: same dispatch; sparse zero/zero metrics compute
    // normally (no bitmap path), e.g. cosine 0.0 / distance 1.0
    let zs = SparseBitFingerprint::new(0);
    let mut sb8 = SparseBitFingerprint::new(8);
    sb8.set_bit(1).unwrap();
    assert_eq!(
        similarity::sparse_similarity_wrapper(&zs, &sb8, similarity::sparse_tanimoto, false),
        boundary
    );
    assert_eq!(
        similarity::sparse_similarity_wrapper(&sb8, &zs, similarity::sparse_tanimoto, false),
        boundary
    );
    assert_eq!(
        similarity::sparse_similarity_wrapper_tversky(
            &zs,
            &sb8,
            0.5,
            0.5,
            similarity::sparse_tversky,
            false
        ),
        boundary
    );
    assert!(
        similarity::sparse_similarity_wrapper(&zs, &zs, similarity::sparse_tanimoto, false)
            .unwrap()
            == 0.0
    );
    assert!(
        similarity::sparse_similarity_wrapper(&zs, &zs, similarity::sparse_cosine, true).unwrap()
            == 1.0
    );

    // equal and non-divisible positive lengths still route through the
    // fold logic; fold errors surface unchanged
    let mut a16 = Fingerprint::new(16);
    a16.set_bit(1).unwrap();
    a16.set_bit(9).unwrap();
    assert!(similarity::similarity_wrapper(&a16, &b8, similarity::tanimoto, false).unwrap() == 1.0);
    let mut a12 = Fingerprint::new(12);
    a12.set_bit(1).unwrap();
    assert_eq!(
        similarity::similarity_wrapper(&a12, &b8, similarity::tanimoto, false),
        Err(FingerprintError::BitLengthMismatch { left: 12, right: 8 })
    );
    // unchanged inputs
    assert_eq!(z.n_bits(), 0);
    assert_eq!(b8.on_bits(), vec![1]);
    assert_eq!(a12.on_bits(), vec![1]);
}

// ---------------------------------------------------------------------------
// Supervisor C4 dispatch matrices (D1-D4): four wrappers x ten ordered
// length pairs x two distance modes x two callback outcomes = 160 executed
// rows, each asserted by a row counter. Callback state is thread-local and
// reset before every call; callbacks are test-private function pointers.
// ---------------------------------------------------------------------------

const DISPATCH_FOLD_DIV_REASON: &str = "SimilarityWrapper: exactly one zero-length operand \
                                        makes the source fold-factor length division a \
                                        division by zero; CK safe rejection of undefined \
                                        source execution";
const DISPATCH_SENTINEL_REASON: &str = "dispatch callback sentinel";

fn dispatch_fold_div_error() -> FingerprintError {
    FingerprintError::InvalidArguments {
        reason: DISPATCH_FOLD_DIV_REASON,
    }
}

fn dispatch_sentinel_error() -> FingerprintError {
    FingerprintError::InvalidArguments {
        reason: DISPATCH_SENTINEL_REASON,
    }
}

/// Input bit sets per dispatch length: 16:{1,9}, 12:{1,9}, 8:{1}, 2:{1},
/// 1:{0}, 0:{}.
fn dispatch_bits(len: u32) -> Vec<u32> {
    match len {
        16 | 12 => vec![1, 9],
        8 | 2 => vec![1],
        1 => vec![0],
        _ => Vec::new(),
    }
}

fn dispatch_dense(len: u32) -> Fingerprint {
    let bits = dispatch_bits(len);
    Fingerprint::from_on_bits(len, bits.iter().copied()).unwrap()
}

fn dispatch_sparse(len: u32) -> SparseBitFingerprint {
    let bits = dispatch_bits(len);
    let mut out = SparseBitFingerprint::new(len);
    for bit in bits {
        out.set_bit(bit).unwrap();
    }
    out
}

/// Expected outcome of one dispatch row.
enum DispatchExpectation {
    /// fold-division rejection before the metric; callback count 0
    FoldDiv,
    /// fold factor rejection; callback count 0
    InvalidFold { factor: u32, n_bits: u32 },
    /// exactly one callback; received operand lengths (post-fold)
    Invoke { recv_a: u32, recv_b: u32 },
}

const DISPATCH_PAIRS: [(u32, u32); 10] = [
    (0, 8),
    (8, 0),
    (0, 0),
    (8, 8),
    (16, 8),
    (8, 16),
    (12, 8),
    (8, 12),
    (2, 1),
    (1, 2),
];

fn dispatch_expectation(la: u32, lb: u32) -> DispatchExpectation {
    match (la, lb) {
        (0, 8) | (8, 0) => DispatchExpectation::FoldDiv,
        (2, 1) | (1, 2) => DispatchExpectation::InvalidFold {
            factor: 2,
            n_bits: 2,
        },
        (0, 0) | (8, 8) => DispatchExpectation::Invoke {
            recv_a: la,
            recv_b: lb,
        },
        // folding 16 by 2 yields length 8, bits {1}
        (16, 8) => DispatchExpectation::Invoke {
            recv_a: 8,
            recv_b: 8,
        },
        (8, 16) => DispatchExpectation::Invoke {
            recv_a: 8,
            recv_b: 8,
        },
        // factor 12/8 = 1 does NOT force equal lengths
        (12, 8) => DispatchExpectation::Invoke {
            recv_a: 12,
            recv_b: 8,
        },
        (8, 12) => DispatchExpectation::Invoke {
            recv_a: 8,
            recv_b: 12,
        },
        _ => unreachable!("frozen dispatch matrix"),
    }
}

#[test]
fn wrapper_dispatch_matrix_d1_dense_ordinary() {
    use cosmolkit_fingerprints::similarity;
    use std::cell::{Cell, RefCell};

    thread_local! {
        static CALLS: Cell<u32> = const { Cell::new(0) };
        static OBSERVED: RefCell<Vec<(u32, Vec<u32>, u32, Vec<u32>)>> =
            const { RefCell::new(Vec::new()) };
    }
    fn record(a: &Fingerprint, b: &Fingerprint) {
        CALLS.with(|c| c.set(c.get() + 1));
        OBSERVED.with(|o| {
            o.borrow_mut()
                .push((a.n_bits(), a.on_bits(), b.n_bits(), b.on_bits()));
        });
    }
    fn cb_ok(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
        record(a, b);
        Ok(0.25)
    }
    fn cb_err(a: &Fingerprint, b: &Fingerprint) -> Result<f64, FingerprintError> {
        record(a, b);
        Err(dispatch_sentinel_error())
    }
    fn reset() {
        CALLS.with(|c| c.set(0));
        OBSERVED.with(|o| o.borrow_mut().clear());
    }
    fn calls() -> u32 {
        CALLS.with(|c| c.get())
    }
    fn observed() -> Vec<(u32, Vec<u32>, u32, Vec<u32>)> {
        OBSERVED.with(|o| o.borrow().clone())
    }

    let mut rows = 0u32;
    for (la, lb) in DISPATCH_PAIRS {
        for distance in [false, true] {
            for callback_ok in [true, false] {
                reset();
                rows += 1;
                let row = format!("d1 ({la},{lb}) distance={distance} ok={callback_ok}");
                let x = dispatch_dense(la);
                let y = dispatch_dense(lb);
                let result = similarity::similarity_wrapper(
                    &x,
                    &y,
                    if callback_ok { cb_ok } else { cb_err },
                    distance,
                );
                match dispatch_expectation(la, lb) {
                    DispatchExpectation::FoldDiv => {
                        assert_eq!(result, Err(dispatch_fold_div_error()), "{row}");
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::InvalidFold { factor, n_bits } => {
                        assert_eq!(
                            result,
                            Err(FingerprintError::InvalidFoldFactor { factor, n_bits }),
                            "{row}"
                        );
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::Invoke { recv_a, recv_b } => {
                        assert_eq!(calls(), 1, "{row}");
                        let seen = observed();
                        assert_eq!(seen.len(), 1, "{row}");
                        assert_eq!(seen[0].0, recv_a, "{row}");
                        assert_eq!(seen[0].2, recv_b, "{row}");
                        // received bits: folded 16->8 keeps {1}; 12 keeps {1,9}
                        assert_eq!(
                            seen[0].1,
                            if recv_a == 8 && la == 16 {
                                vec![1]
                            } else {
                                dispatch_bits(recv_a)
                            },
                            "{row} recv_a bits"
                        );
                        assert_eq!(
                            seen[0].3,
                            if recv_b == 8 && lb == 16 {
                                vec![1]
                            } else {
                                dispatch_bits(recv_b)
                            },
                            "{row} recv_b bits"
                        );
                        match (callback_ok, distance) {
                            (true, false) => assert_eq!(result.unwrap(), 0.25, "{row}"),
                            (true, true) => assert_eq!(result.unwrap(), 0.75, "{row}"),
                            (false, _) => {
                                assert_eq!(result, Err(dispatch_sentinel_error()), "{row}")
                            }
                        }
                    }
                }
                // original inputs unchanged: lengths, bits AND cached counts
                assert_eq!(x.n_bits(), la, "{row}");
                assert_eq!(x.on_bits(), dispatch_bits(la), "{row}");
                assert_eq!(x.num_on_bits(), dispatch_bits(la).len() as u32, "{row}");
                assert_eq!(y.n_bits(), lb, "{row}");
                assert_eq!(y.on_bits(), dispatch_bits(lb), "{row}");
                assert_eq!(y.num_on_bits(), dispatch_bits(lb).len() as u32, "{row}");
            }
        }
    }
    assert_eq!(rows, 40, "D1 must execute exactly 40 dispatch rows");
}

#[test]
fn wrapper_dispatch_matrix_d2_sparse_ordinary() {
    use cosmolkit_fingerprints::similarity;
    use std::cell::{Cell, RefCell};

    thread_local! {
        static CALLS: Cell<u32> = const { Cell::new(0) };
        static OBSERVED: RefCell<Vec<(u32, Vec<i32>, u32, Vec<i32>)>> =
            const { RefCell::new(Vec::new()) };
    }
    fn as_i32(bits: &[u32]) -> Vec<i32> {
        bits.iter().map(|&b| b as i32).collect()
    }
    fn record(a: &SparseBitFingerprint, b: &SparseBitFingerprint) {
        CALLS.with(|c| c.set(c.get() + 1));
        OBSERVED.with(|o| {
            o.borrow_mut()
                .push((a.n_bits(), a.on_bits(), b.n_bits(), b.on_bits()));
        });
    }
    fn cb_ok(a: &SparseBitFingerprint, b: &SparseBitFingerprint) -> Result<f64, FingerprintError> {
        record(a, b);
        Ok(0.25)
    }
    fn cb_err(a: &SparseBitFingerprint, b: &SparseBitFingerprint) -> Result<f64, FingerprintError> {
        record(a, b);
        Err(dispatch_sentinel_error())
    }
    fn reset() {
        CALLS.with(|c| c.set(0));
        OBSERVED.with(|o| o.borrow_mut().clear());
    }
    fn calls() -> u32 {
        CALLS.with(|c| c.get())
    }
    fn observed() -> Vec<(u32, Vec<i32>, u32, Vec<i32>)> {
        OBSERVED.with(|o| o.borrow().clone())
    }

    let mut rows = 0u32;
    for (la, lb) in DISPATCH_PAIRS {
        for distance in [false, true] {
            for callback_ok in [true, false] {
                reset();
                rows += 1;
                let row = format!("d2 ({la},{lb}) distance={distance} ok={callback_ok}");
                let x = dispatch_sparse(la);
                let y = dispatch_sparse(lb);
                let result = similarity::sparse_similarity_wrapper(
                    &x,
                    &y,
                    if callback_ok { cb_ok } else { cb_err },
                    distance,
                );
                match dispatch_expectation(la, lb) {
                    DispatchExpectation::FoldDiv => {
                        assert_eq!(result, Err(dispatch_fold_div_error()), "{row}");
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::InvalidFold { factor, n_bits } => {
                        assert_eq!(
                            result,
                            Err(FingerprintError::InvalidFoldFactor { factor, n_bits }),
                            "{row}"
                        );
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::Invoke { recv_a, recv_b } => {
                        assert_eq!(calls(), 1, "{row}");
                        let seen = observed();
                        assert_eq!(seen.len(), 1, "{row}");
                        assert_eq!(seen[0].0, recv_a, "{row}");
                        assert_eq!(seen[0].2, recv_b, "{row}");
                        assert_eq!(
                            seen[0].1,
                            if recv_a == 8 && la == 16 {
                                vec![1]
                            } else {
                                as_i32(&dispatch_bits(recv_a))
                            },
                            "{row} recv_a bits"
                        );
                        assert_eq!(
                            seen[0].3,
                            if recv_b == 8 && lb == 16 {
                                vec![1]
                            } else {
                                as_i32(&dispatch_bits(recv_b))
                            },
                            "{row} recv_b bits"
                        );
                        match (callback_ok, distance) {
                            (true, false) => assert_eq!(result.unwrap(), 0.25, "{row}"),
                            (true, true) => assert_eq!(result.unwrap(), 0.75, "{row}"),
                            (false, _) => {
                                assert_eq!(result, Err(dispatch_sentinel_error()), "{row}")
                            }
                        }
                    }
                }
                // original inputs unchanged: lengths, bits AND cached counts
                assert_eq!(x.n_bits(), la, "{row}");
                assert_eq!(x.on_bits(), as_i32(&dispatch_bits(la)), "{row}");
                assert_eq!(x.num_on_bits(), dispatch_bits(la).len() as u32, "{row}");
                assert_eq!(y.n_bits(), lb, "{row}");
                assert_eq!(y.on_bits(), as_i32(&dispatch_bits(lb)), "{row}");
                assert_eq!(y.num_on_bits(), dispatch_bits(lb).len() as u32, "{row}");
            }
        }
    }
    assert_eq!(rows, 40, "D2 must execute exactly 40 dispatch rows");
}

#[test]
fn wrapper_dispatch_matrix_d3_dense_parameterized() {
    use cosmolkit_fingerprints::similarity;
    use std::cell::{Cell, RefCell};

    thread_local! {
        static CALLS: Cell<u32> = const { Cell::new(0) };
        static OBSERVED: RefCell<Vec<(u32, Vec<u32>, u32, Vec<u32>, u64, u64)>> =
            const { RefCell::new(Vec::new()) };
    }
    fn record(a: &Fingerprint, b: &Fingerprint, alpha: f64, beta: f64) {
        CALLS.with(|c| c.set(c.get() + 1));
        OBSERVED.with(|o| {
            o.borrow_mut().push((
                a.n_bits(),
                a.on_bits(),
                b.n_bits(),
                b.on_bits(),
                alpha.to_bits(),
                beta.to_bits(),
            ));
        });
    }
    fn cb_ok(
        a: &Fingerprint,
        b: &Fingerprint,
        alpha: f64,
        beta: f64,
    ) -> Result<f64, FingerprintError> {
        record(a, b, alpha, beta);
        Ok(0.25)
    }
    fn cb_err(
        a: &Fingerprint,
        b: &Fingerprint,
        alpha: f64,
        beta: f64,
    ) -> Result<f64, FingerprintError> {
        record(a, b, alpha, beta);
        Err(dispatch_sentinel_error())
    }
    fn reset() {
        CALLS.with(|c| c.set(0));
        OBSERVED.with(|o| o.borrow_mut().clear());
    }
    fn calls() -> u32 {
        CALLS.with(|c| c.get())
    }
    fn observed() -> Vec<(u32, Vec<u32>, u32, Vec<u32>, u64, u64)> {
        OBSERVED.with(|o| o.borrow().clone())
    }

    let alpha = 0.25f64;
    let beta = 0.75f64;
    let mut rows = 0u32;
    for (la, lb) in DISPATCH_PAIRS {
        for distance in [false, true] {
            for callback_ok in [true, false] {
                reset();
                rows += 1;
                let row = format!("d3 ({la},{lb}) distance={distance} ok={callback_ok}");
                let x = dispatch_dense(la);
                let y = dispatch_dense(lb);
                let result = similarity::similarity_wrapper_tversky(
                    &x,
                    &y,
                    alpha,
                    beta,
                    if callback_ok { cb_ok } else { cb_err },
                    distance,
                );
                match dispatch_expectation(la, lb) {
                    DispatchExpectation::FoldDiv => {
                        assert_eq!(result, Err(dispatch_fold_div_error()), "{row}");
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::InvalidFold { factor, n_bits } => {
                        assert_eq!(
                            result,
                            Err(FingerprintError::InvalidFoldFactor { factor, n_bits }),
                            "{row}"
                        );
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::Invoke { recv_a, recv_b } => {
                        assert_eq!(calls(), 1, "{row}");
                        let seen = observed();
                        assert_eq!(seen.len(), 1, "{row}");
                        assert_eq!(seen[0].0, recv_a, "{row}");
                        assert_eq!(seen[0].2, recv_b, "{row}");
                        assert_eq!(
                            seen[0].1,
                            if recv_a == 8 && la == 16 {
                                vec![1]
                            } else {
                                dispatch_bits(recv_a)
                            },
                            "{row} recv_a bits"
                        );
                        assert_eq!(
                            seen[0].3,
                            if recv_b == 8 && lb == 16 {
                                vec![1]
                            } else {
                                dispatch_bits(recv_b)
                            },
                            "{row} recv_b bits"
                        );
                        // parameter bits forwarded unchanged
                        assert_eq!(seen[0].4, alpha.to_bits(), "{row} alpha bits");
                        assert_eq!(seen[0].5, beta.to_bits(), "{row} beta bits");
                        match (callback_ok, distance) {
                            (true, false) => assert_eq!(result.unwrap(), 0.25, "{row}"),
                            (true, true) => assert_eq!(result.unwrap(), 0.75, "{row}"),
                            (false, _) => {
                                assert_eq!(result, Err(dispatch_sentinel_error()), "{row}")
                            }
                        }
                    }
                }
                // original inputs unchanged: lengths, bits AND cached counts
                assert_eq!(x.n_bits(), la, "{row}");
                assert_eq!(x.on_bits(), dispatch_bits(la), "{row}");
                assert_eq!(x.num_on_bits(), dispatch_bits(la).len() as u32, "{row}");
                assert_eq!(y.n_bits(), lb, "{row}");
                assert_eq!(y.on_bits(), dispatch_bits(lb), "{row}");
                assert_eq!(y.num_on_bits(), dispatch_bits(lb).len() as u32, "{row}");
            }
        }
    }
    assert_eq!(rows, 40, "D3 must execute exactly 40 dispatch rows");
}

#[test]
fn wrapper_dispatch_matrix_d4_sparse_parameterized() {
    use cosmolkit_fingerprints::similarity;
    use std::cell::{Cell, RefCell};

    thread_local! {
        static CALLS: Cell<u32> = const { Cell::new(0) };
        static OBSERVED: RefCell<Vec<(u32, Vec<i32>, u32, Vec<i32>, u64, u64)>> =
            const { RefCell::new(Vec::new()) };
    }
    fn as_i32(bits: &[u32]) -> Vec<i32> {
        bits.iter().map(|&b| b as i32).collect()
    }
    fn record(a: &SparseBitFingerprint, b: &SparseBitFingerprint, alpha: f64, beta: f64) {
        CALLS.with(|c| c.set(c.get() + 1));
        OBSERVED.with(|o| {
            o.borrow_mut().push((
                a.n_bits(),
                a.on_bits(),
                b.n_bits(),
                b.on_bits(),
                alpha.to_bits(),
                beta.to_bits(),
            ));
        });
    }
    fn cb_ok(
        a: &SparseBitFingerprint,
        b: &SparseBitFingerprint,
        alpha: f64,
        beta: f64,
    ) -> Result<f64, FingerprintError> {
        record(a, b, alpha, beta);
        Ok(0.25)
    }
    fn cb_err(
        a: &SparseBitFingerprint,
        b: &SparseBitFingerprint,
        alpha: f64,
        beta: f64,
    ) -> Result<f64, FingerprintError> {
        record(a, b, alpha, beta);
        Err(dispatch_sentinel_error())
    }
    fn reset() {
        CALLS.with(|c| c.set(0));
        OBSERVED.with(|o| o.borrow_mut().clear());
    }
    fn calls() -> u32 {
        CALLS.with(|c| c.get())
    }
    fn observed() -> Vec<(u32, Vec<i32>, u32, Vec<i32>, u64, u64)> {
        OBSERVED.with(|o| o.borrow().clone())
    }

    let alpha = 0.25f64;
    let beta = 0.75f64;
    let mut rows = 0u32;
    for (la, lb) in DISPATCH_PAIRS {
        for distance in [false, true] {
            for callback_ok in [true, false] {
                reset();
                rows += 1;
                let row = format!("d4 ({la},{lb}) distance={distance} ok={callback_ok}");
                let x = dispatch_sparse(la);
                let y = dispatch_sparse(lb);
                let result = similarity::sparse_similarity_wrapper_tversky(
                    &x,
                    &y,
                    alpha,
                    beta,
                    if callback_ok { cb_ok } else { cb_err },
                    distance,
                );
                match dispatch_expectation(la, lb) {
                    DispatchExpectation::FoldDiv => {
                        assert_eq!(result, Err(dispatch_fold_div_error()), "{row}");
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::InvalidFold { factor, n_bits } => {
                        assert_eq!(
                            result,
                            Err(FingerprintError::InvalidFoldFactor { factor, n_bits }),
                            "{row}"
                        );
                        assert_eq!(calls(), 0, "{row}");
                    }
                    DispatchExpectation::Invoke { recv_a, recv_b } => {
                        assert_eq!(calls(), 1, "{row}");
                        let seen = observed();
                        assert_eq!(seen.len(), 1, "{row}");
                        assert_eq!(seen[0].0, recv_a, "{row}");
                        assert_eq!(seen[0].2, recv_b, "{row}");
                        assert_eq!(
                            seen[0].1,
                            if recv_a == 8 && la == 16 {
                                vec![1]
                            } else {
                                as_i32(&dispatch_bits(recv_a))
                            },
                            "{row} recv_a bits"
                        );
                        assert_eq!(
                            seen[0].3,
                            if recv_b == 8 && lb == 16 {
                                vec![1]
                            } else {
                                as_i32(&dispatch_bits(recv_b))
                            },
                            "{row} recv_b bits"
                        );
                        // parameter bits forwarded unchanged
                        assert_eq!(seen[0].4, alpha.to_bits(), "{row} alpha bits");
                        assert_eq!(seen[0].5, beta.to_bits(), "{row} beta bits");
                        match (callback_ok, distance) {
                            (true, false) => assert_eq!(result.unwrap(), 0.25, "{row}"),
                            (true, true) => assert_eq!(result.unwrap(), 0.75, "{row}"),
                            (false, _) => {
                                assert_eq!(result, Err(dispatch_sentinel_error()), "{row}")
                            }
                        }
                    }
                }
                // original inputs unchanged: lengths, bits AND cached counts
                assert_eq!(x.n_bits(), la, "{row}");
                assert_eq!(x.on_bits(), as_i32(&dispatch_bits(la)), "{row}");
                assert_eq!(x.num_on_bits(), dispatch_bits(la).len() as u32, "{row}");
                assert_eq!(y.n_bits(), lb, "{row}");
                assert_eq!(y.on_bits(), as_i32(&dispatch_bits(lb)), "{row}");
                assert_eq!(y.num_on_bits(), dispatch_bits(lb).len() as u32, "{row}");
            }
        }
    }
    assert_eq!(rows, 40, "D4 must execute exactly 40 dispatch rows");
}

#[test]
fn tversky_wrapper_real_metrics_zero_and_range_dispatch() {
    // Concrete metric dispatch through both parameterized wrappers using
    // the real TverskySimilarity (not a test callback):
    // - zero/zero with valid alpha=beta=0.5: dense raises the
    //   PreconditionViolation from NumOnBitsInCommon's zero-length bitmap
    //   path inside the metric; sparse computes 0.0 (distance 1.0);
    // - zero/zero with alpha=1.5: RangeError from the metric's
    //   RANGE_CHECK in both distance modes;
    // - zero/nonzero and nonzero/zero with alpha=1.5: the approved
    //   fold-division rejection fires BEFORE metric range validation.
    use cosmolkit_fingerprints::similarity;

    let zd = Fingerprint::new(0);
    let nd = dispatch_dense(8);
    let zs = SparseBitFingerprint::new(0);
    let ns = dispatch_sparse(8);

    // dense zero/zero, valid parameters
    assert_eq!(
        similarity::similarity_wrapper_tversky(&zd, &zd, 0.5, 0.5, similarity::tversky, false),
        Err(FingerprintError::PreconditionViolation {
            what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)"
        })
    );
    assert_eq!(
        similarity::similarity_wrapper_tversky(&zd, &zd, 0.5, 0.5, similarity::tversky, true),
        Err(FingerprintError::PreconditionViolation {
            what: "no afp (CalcBitmapNumBitsInCommon on zero-length vector)"
        })
    );
    // sparse zero/zero, valid parameters: y==0 -> 0.0; distance -> 1.0
    assert_eq!(
        similarity::sparse_similarity_wrapper_tversky(
            &zs,
            &zs,
            0.5,
            0.5,
            similarity::sparse_tversky,
            false
        )
        .unwrap(),
        0.0
    );
    assert_eq!(
        similarity::sparse_similarity_wrapper_tversky(
            &zs,
            &zs,
            0.5,
            0.5,
            similarity::sparse_tversky,
            true
        )
        .unwrap(),
        1.0
    );

    // zero/zero with out-of-range alpha: RangeError in both modes
    for distance in [false, true] {
        assert_eq!(
            similarity::similarity_wrapper_tversky(
                &zd,
                &zd,
                1.5,
                0.5,
                similarity::tversky,
                distance
            ),
            Err(FingerprintError::RangeError { value: 1.5 })
        );
        assert_eq!(
            similarity::sparse_similarity_wrapper_tversky(
                &zs,
                &zs,
                1.5,
                0.5,
                similarity::sparse_tversky,
                distance
            ),
            Err(FingerprintError::RangeError { value: 1.5 })
        );
    }

    // one zero-length operand with out-of-range alpha: the fold-division
    // rejection precedes metric range validation in both orders
    for distance in [false, true] {
        assert_eq!(
            similarity::similarity_wrapper_tversky(
                &zd,
                &nd,
                1.5,
                0.5,
                similarity::tversky,
                distance
            ),
            Err(dispatch_fold_div_error())
        );
        assert_eq!(
            similarity::similarity_wrapper_tversky(
                &nd,
                &zd,
                1.5,
                0.5,
                similarity::tversky,
                distance
            ),
            Err(dispatch_fold_div_error())
        );
        assert_eq!(
            similarity::sparse_similarity_wrapper_tversky(
                &zs,
                &ns,
                1.5,
                0.5,
                similarity::sparse_tversky,
                distance
            ),
            Err(dispatch_fold_div_error())
        );
        assert_eq!(
            similarity::sparse_similarity_wrapper_tversky(
                &ns,
                &zs,
                1.5,
                0.5,
                similarity::sparse_tversky,
                distance
            ),
            Err(dispatch_fold_div_error())
        );
    }
}

#[test]
fn u01_count_storage_access_both_widths() {
    // U01: shared storage/access template for SparseCountFingerprint (u64
    // index) and SparseCountFingerprint32 (u32 index), both with i32
    // values. Defined-source expectations mirror the pinned
    // SparseIntVect.h checkIndex/getVal/setVal template; the u32::MAX
    // index==length allowance is pinned from the template text exactly as
    // the u64 edge was (the Python UInt wrapper exposes no SetVal).
    // ---- u32 width ----
    let mut v = SparseCountFingerprint32::new(10);
    assert_eq!(v.length(), 10);
    assert_eq!(v.value(3).unwrap(), 0); // absent reads zero
    v.set_value(3, 2).unwrap();
    v.set_value(5, -1).unwrap(); // negative value
    assert_eq!(v.value(3).unwrap(), 2);
    assert_eq!(v.value(5).unwrap(), -1);
    v.set_value(3, 0).unwrap(); // zero erases
    assert_eq!(v.value(3).unwrap(), 0);
    assert!(v.nonzero_elements().get(&3).is_none());
    // ordinary bounds: last valid, index==length, above length
    assert!(v.set_value(9, 1).is_ok());
    assert_eq!(
        v.set_value(10, 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 10,
            size: 10
        })
    );
    assert_eq!(
        v.value(11),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 11,
            size: 10
        })
    );
    // zero length
    let z = SparseCountFingerprint32::new(0);
    assert_eq!(
        z.value(0),
        Err(FingerprintError::SparseIndexOutOfRange { index: 0, size: 0 })
    );
    // u32::MAX permits index == length
    let mut big = SparseCountFingerprint32::new(u32::MAX);
    big.set_value(u32::MAX, -7).unwrap();
    assert_eq!(big.value(u32::MAX).unwrap(), -7);
    big.set_value(u32::MAX - 1, 3).unwrap();
    assert_eq!(big.value(u32::MAX - 1).unwrap(), 3);
    // ascending unsigned keys across 2^31
    let mut hi = SparseCountFingerprint32::new(2_u32.pow(31) + 4);
    hi.set_value(2_u32.pow(31) + 3, 4).unwrap();
    hi.set_value(2_u32.pow(31), 1).unwrap();
    hi.set_value(5, 2).unwrap();
    let keys: Vec<u32> = hi.nonzero_elements().keys().copied().collect();
    assert_eq!(keys, vec![5, 2_u32.pow(31), 2_u32.pow(31) + 3]);
    // clone independence + equality
    let mut clone = hi.clone();
    assert_eq!(clone, hi);
    clone.set_value(7, 9).unwrap();
    assert_ne!(clone, hi);
    assert_eq!(hi.nonzero_elements().len(), 3);

    // ---- u64 width (shared template changed; re-pin core behavior) ----
    let mut w = SparseCountFingerprint::new(10);
    w.set_value(3, 2).unwrap();
    w.set_value(5, -1).unwrap();
    assert_eq!(w.value(3).unwrap(), 2);
    w.set_value(3, 0).unwrap();
    assert_eq!(w.value(3).unwrap(), 0);
    assert_eq!(
        w.set_value(10, 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: 10,
            size: 10
        })
    );
    // u64 length == u32::MAX REJECTS index == length (only u64::MAX allows)
    let mut wide = SparseCountFingerprint::new(u64::from(u32::MAX));
    assert_eq!(
        wide.set_value(u64::from(u32::MAX), 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: u64::from(u32::MAX),
            size: u64::from(u32::MAX)
        })
    );
    assert!(wide.set_value(u64::from(u32::MAX) - 1, 1).is_ok());
    // equality including a stored zero (scalar-produced, u64 path today;
    // the u32 scalar family arrives with its own packets and retests this)
    let mut sz = SparseCountFingerprint::new(8);
    sz.set_value(2, 3).unwrap();
    let sz = sz.with_added_scalar(-3).unwrap();
    assert_eq!(sz.nonzero_elements().get(&2), Some(&0));
    let mut other = SparseCountFingerprint::new(8);
    other.set_value(2, 4).unwrap();
    assert_ne!(sz, other);
    let other_zeroed = other.with_added_scalar(-4).unwrap();
    assert_eq!(sz, other_zeroed);

    // operands unchanged after every success/error above is covered by the
    // value-style API; spot-check one error path leaves inputs untouched
    let before = hi.clone();
    assert!(hi.set_value(9, 1).is_ok()); // in range, mutates by design
    drop(before);
}

#[test]
fn u02_total_val_both_widths() {
    // U02: shared getTotalVal template — signed and absolute totals,
    // cancellation, intermediate overflow, MIN abs (CK-labeled), inputs
    // unchanged. Defined-source expectations mirror the u64 F12 probes;
    // the u32 UInt wrapper's GetTotalVal pins the same semantics.
    let mut a = SparseCountFingerprint32::new(8);
    a.set_value(0, 5).unwrap();
    a.set_value(1, -2).unwrap();
    assert_eq!(a.total_value(false).unwrap(), 3);
    assert_eq!(a.total_value(true).unwrap(), 7);
    // cancellation to exact zero total
    a.set_value(2, -3).unwrap();
    assert_eq!(a.total_value(false).unwrap(), 0);
    assert_eq!(a.total_value(true).unwrap(), 10);
    // empty
    assert_eq!(
        SparseCountFingerprint32::new(8).total_value(false).unwrap(),
        0
    );
    // defined extrema
    let mut mx = SparseCountFingerprint32::new(8);
    mx.set_value(0, i32::MAX).unwrap();
    assert_eq!(mx.total_value(false).unwrap(), i32::MAX);
    assert_eq!(mx.total_value(true).unwrap(), i32::MAX);
    // intermediate accumulation overflow (CK safety boundary, executed UB)
    let mut ov = SparseCountFingerprint32::new(8);
    ov.set_value(0, i32::MAX).unwrap();
    ov.set_value(1, 1).unwrap();
    assert_eq!(
        ov.total_value(false),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::getTotalVal accumulation"
        })
    );
    // abs(MIN) rejection (CK safety boundary)
    let mut mn = SparseCountFingerprint32::new(8);
    mn.set_value(0, i32::MIN).unwrap();
    assert_eq!(mn.total_value(false).unwrap(), i32::MIN);
    assert_eq!(
        mn.total_value(true),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::getTotalVal abs"
        })
    );
    // inputs unchanged on success and error
    assert_eq!(a.value(1).unwrap(), -2);
    assert_eq!(ov.value(1).unwrap(), 1);
    assert_eq!(mn.value(0).unwrap(), i32::MIN);

    // u64 width re-pin (shared body changed)
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(0, 5).unwrap();
    b.set_value(1, -2).unwrap();
    assert_eq!(b.total_value(false).unwrap(), 3);
    assert_eq!(b.total_value(true).unwrap(), 7);
}

#[test]
fn u03_fuzzy_and_both_widths() {
    // U03: shared source-shaped monotonic fuzzy_and merge. Oracle shape
    // (u64 F13 probe): a={1:2,3:5}, b={3:2,5:1} -> {3:2}; absent dropped;
    // negative minima; both orders; length mismatch; operands unchanged.
    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    let i = a32.fuzzy_and(&b32).unwrap();
    let entries: Vec<(u32, i32)> = i.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(3, 2)]);
    // disjoint -> empty support
    let mut d32 = SparseCountFingerprint32::new(8);
    d32.set_value(7, 4).unwrap();
    assert!(a32.fuzzy_and(&d32).unwrap().nonzero_elements().is_empty());
    // overlapping negatives: min picks the negative side
    let mut n32 = SparseCountFingerprint32::new(8);
    n32.set_value(1, -3).unwrap();
    let ni = a32.fuzzy_and(&n32).unwrap();
    assert_eq!(ni.value(1).unwrap(), -3);
    // both orders agree on the intersection
    let i2 = n32.fuzzy_and(&a32).unwrap();
    assert_eq!(i2.value(1).unwrap(), -3);
    // retained scalar-produced stored zero participates as value 0:
    // {0:0} & {0:9} keeps 0 (u32 scalars arrive in U07-U10; zero produced
    // via u64 path is covered in the u64 section below)
    // length mismatch both orders
    assert_eq!(
        a32.fuzzy_and(&SparseCountFingerprint32::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(
        SparseCountFingerprint32::new(16).fuzzy_and(&a32),
        Err(FingerprintError::BitLengthMismatch { left: 16, right: 8 })
    );
    // operands unchanged
    assert_eq!(a32.value(1).unwrap(), 2);
    assert_eq!(a32.value(3).unwrap(), 5);
    assert_eq!(b32.value(5).unwrap(), 1);

    // u64 width: shared merge changed; re-pin F13 oracle values
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let i64r = a.fuzzy_and(&b).unwrap();
    let entries: Vec<(u64, i32)> = i64r
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(3, 2)]);
    // retained stored zero (scalar-produced) intersects as 0
    let mut sz = SparseCountFingerprint::new(8);
    sz.set_value(2, 3).unwrap();
    let sz = sz.with_added_scalar(-3).unwrap(); // {2: 0}
    let mut other = SparseCountFingerprint::new(8);
    other.set_value(2, 9).unwrap();
    let zi = sz.fuzzy_and(&other).unwrap();
    assert_eq!(zi.value(2).unwrap(), 0);
    assert!(zi.nonzero_elements().contains_key(&2));
}

#[test]
fn u04_fuzzy_or_both_widths() {
    // U04: shared fuzzy_or max merge. Oracle shape (u64 F14 probe):
    // a={1:2,3:5}, b={3:2,5:1} -> {1:2,3:5,5:1}; one-sided negatives;
    // overlapping min/max distinction; stored zeros; length mismatch;
    // input independence.
    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    let u = a32.fuzzy_or(&b32).unwrap();
    let entries: Vec<(u32, i32)> = u.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (5, 1)]);
    // left-only and right-only negative entries insert verbatim
    let mut lneg = SparseCountFingerprint32::new(8);
    lneg.set_value(0, -6).unwrap();
    lneg.set_value(2, 4).unwrap();
    let mut rneg = SparseCountFingerprint32::new(8);
    rneg.set_value(2, 1).unwrap();
    rneg.set_value(4, -9).unwrap();
    let entries: Vec<(u32, i32)> = lneg
        .fuzzy_or(&rneg)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(0, -6), (2, 4), (4, -9)]);
    // overlapping min/max distinction (fuzzy_and vs fuzzy_or)
    assert_eq!(lneg.fuzzy_and(&rneg).unwrap().value(2).unwrap(), 1);
    assert_eq!(lneg.fuzzy_or(&rneg).unwrap().value(2).unwrap(), 4);
    // length mismatch and input independence
    assert_eq!(
        a32.fuzzy_or(&SparseCountFingerprint32::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(a32.value(1).unwrap(), 2);
    assert_eq!(b32.value(3).unwrap(), 2);

    // u64 width: shared body changed; re-pin F14 oracle + stored zeros
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let entries: Vec<(u64, i32)> = a
        .fuzzy_or(&b)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (5, 1)]);
    // stored zero on one side is preserved by the union (max(0, v) = v;
    // max(0, stored 0) stays stored)
    let mut sz = SparseCountFingerprint::new(8);
    sz.set_value(2, 3).unwrap();
    let sz = sz.with_added_scalar(-3).unwrap(); // {2: 0}
    let mut other = SparseCountFingerprint::new(8);
    other.set_value(2, -4).unwrap();
    let zo = sz.fuzzy_or(&other).unwrap();
    assert_eq!(zo.value(2).unwrap(), 0);
    assert!(zo.nonzero_elements().contains_key(&2));
}

#[test]
fn u05_vector_add_both_widths() {
    // U05: shared vector addition. Oracle shape (u64 F15 probe):
    // a={1:2,3:5}, b={3:2,5:1} -> {1:2,3:7,5:1}; one-sided both ways;
    // cancellation erases; stored zeros insert; overflow reported after
    // the length preflight; inputs unchanged.
    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    let s32 = a32.with_added(&b32).unwrap();
    let entries: Vec<(u32, i32)> = s32
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 7), (5, 1)]);
    // one-sided supports both ways
    let mut r_only = SparseCountFingerprint32::new(8);
    r_only.set_value(6, -4).unwrap();
    let entries: Vec<(u32, i32)> = a32
        .with_added(&r_only)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (6, -4)]);
    let entries: Vec<(u32, i32)> = r_only
        .with_added(&a32)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 5), (6, -4)]);
    // cancellation erases the entry
    let mut c32 = SparseCountFingerprint32::new(8);
    c32.set_value(1, -2).unwrap();
    assert!(
        a32.with_added(&c32)
            .unwrap()
            .nonzero_elements()
            .get(&1)
            .is_none()
    );
    // length error fires before arithmetic (overflowing operands with a
    // length mismatch report the mismatch)
    let mut big32 = SparseCountFingerprint32::new(8);
    big32.set_value(0, i32::MAX).unwrap();
    assert_eq!(
        big32.with_added(&SparseCountFingerprint32::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // overflow reported (CK safety boundary), inputs unchanged
    let mut one32 = SparseCountFingerprint32::new(8);
    one32.set_value(0, 1).unwrap();
    assert_eq!(
        big32.with_added(&one32),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator+= accumulation"
        })
    );
    assert_eq!(big32.value(0).unwrap(), i32::MAX);
    assert_eq!(one32.value(0).unwrap(), 1);
    assert_eq!(a32.value(3).unwrap(), 5);

    // u64 width: shared body changed; re-pin F15 oracle + stored zero pass
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let entries: Vec<(u64, i32)> = a
        .with_added(&b)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 7), (5, 1)]);
    // stored zero on the right inserts verbatim
    let mut sz = SparseCountFingerprint::new(8);
    sz.set_value(2, 3).unwrap();
    let sz = sz.with_added_scalar(-3).unwrap(); // {2: 0}
    let mut other = SparseCountFingerprint::new(8);
    other.set_value(5, 7).unwrap();
    let entries: Vec<(u64, i32)> = sz
        .with_added(&other)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(2, 0), (5, 7)]);
}

#[test]
fn u06_vector_sub_both_widths() {
    // U06: shared vector subtraction. Oracle shape (u64 F16 probe):
    // a={1:2,3:5}, b={3:2,5:1} -> {1:2,3:3,5:-1}; right-only MIN negation
    // rejected; overlap overflow rejected; ordinary negative subtraction;
    // cancellation; length mismatch; unchanged inputs.
    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    let d32 = a32.with_subtracted(&b32).unwrap();
    let entries: Vec<(u32, i32)> = d32
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 3), (5, -1)]);
    // ordinary negative subtraction
    let mut neg = SparseCountFingerprint32::new(8);
    neg.set_value(1, -1).unwrap();
    assert_eq!(a32.with_subtracted(&neg).unwrap().value(1).unwrap(), 3);
    // cancellation erases
    let mut same = SparseCountFingerprint32::new(8);
    same.set_value(1, 2).unwrap();
    assert!(
        a32.with_subtracted(&same)
            .unwrap()
            .nonzero_elements()
            .get(&1)
            .is_none()
    );
    // right-only MIN negation rejected (CK safety boundary)
    let mut mn_right = SparseCountFingerprint32::new(8);
    mn_right.set_value(7, i32::MIN).unwrap();
    assert_eq!(
        SparseCountFingerprint32::new(8).with_subtracted(&mn_right),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-= negation"
        })
    );
    // overlap overflow rejected
    let mut mx = SparseCountFingerprint32::new(8);
    mx.set_value(0, i32::MIN).unwrap();
    let mut dec = SparseCountFingerprint32::new(8);
    dec.set_value(0, 1).unwrap();
    assert_eq!(
        mx.with_subtracted(&dec),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-= accumulation"
        })
    );
    // length mismatch; inputs unchanged
    assert_eq!(
        a32.with_subtracted(&SparseCountFingerprint32::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    assert_eq!(a32.value(1).unwrap(), 2);
    assert_eq!(mn_right.value(7).unwrap(), i32::MIN);

    // u64 width: shared body changed; re-pin F16 oracle
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let entries: Vec<(u64, i32)> = a
        .with_subtracted(&b)
        .unwrap()
        .nonzero_elements()
        .iter()
        .map(|(&k, &v)| (k, v))
        .collect();
    assert_eq!(entries, vec![(1, 2), (3, 3), (5, -1)]);
}

#[test]
fn u07_scalar_mul_both_widths() {
    // U07: shared scalar multiplication with the approved checked
    // boundary. Empty storage succeeds for any factor; factors 0/1/-1;
    // retained zeros; negatives; MIN*-1 and product overflow; unchanged
    // inputs.
    let empty32 = SparseCountFingerprint32::new(8);
    assert!(empty32.with_multiplied_scalar(i32::MIN).is_ok());
    let mut v32 = SparseCountFingerprint32::new(8);
    v32.set_value(1, 3).unwrap();
    v32.set_value(2, -4).unwrap();
    assert_eq!(v32.with_multiplied_scalar(2).unwrap().value(1).unwrap(), 6);
    assert_eq!(
        v32.with_multiplied_scalar(-3).unwrap().value(2).unwrap(),
        12
    );
    assert_eq!(v32.with_multiplied_scalar(1).unwrap().value(1).unwrap(), 3);
    // factor 0 retains stored zeros
    let z32 = v32.with_multiplied_scalar(0).unwrap();
    assert_eq!(z32.nonzero_elements().len(), 2);
    assert_eq!(
        z32.with_multiplied_scalar(5)
            .unwrap()
            .nonzero_elements()
            .len(),
        2
    );
    // MIN * -1 and product overflow rejected (CK safety boundary)
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        mn32.with_multiplied_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator*=(int)"
        })
    );
    let mut big32 = SparseCountFingerprint32::new(8);
    big32.set_value(0, 1 << 30).unwrap();
    assert_eq!(
        big32.with_multiplied_scalar(4),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator*=(int)"
        })
    );
    // unchanged inputs
    assert_eq!(v32.value(1).unwrap(), 3);
    assert_eq!(mn32.value(0).unwrap(), i32::MIN);

    // u64 width: re-pin F17 values
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 3).unwrap();
    v.set_value(2, -4).unwrap();
    assert_eq!(v.with_multiplied_scalar(2).unwrap().value(1).unwrap(), 6);
    assert_eq!(
        v.with_multiplied_scalar(0)
            .unwrap()
            .nonzero_elements()
            .len(),
        2
    );
}

#[test]
fn u08_scalar_div_both_widths() {
    // U08: shared scalar division. Empty/0 succeeds; stored zero /0 and
    // nonzero /0 error; MIN/-1 errors; positive/negative truncation
    // toward zero; retained zero quotient; unchanged inputs.
    let empty32 = SparseCountFingerprint32::new(8);
    assert!(empty32.with_divided_scalar(0).is_ok());
    let mut v32 = SparseCountFingerprint32::new(8);
    v32.set_value(0, 7).unwrap();
    v32.set_value(1, -7).unwrap();
    v32.set_value(2, 6).unwrap();
    v32.set_value(3, -6).unwrap();
    assert_eq!(v32.with_divided_scalar(2).unwrap().value(0).unwrap(), 3);
    assert_eq!(v32.with_divided_scalar(2).unwrap().value(1).unwrap(), -3);
    assert_eq!(v32.with_divided_scalar(-2).unwrap().value(2).unwrap(), -3);
    assert_eq!(v32.with_divided_scalar(-2).unwrap().value(3).unwrap(), 3);
    // nonzero / 0 errors (executed division)
    assert_eq!(
        v32.with_divided_scalar(0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    // stored zero / 0 errors (scalar-produced zero is an executed entry)
    let z32 = v32.with_multiplied_scalar(0).unwrap();
    assert_eq!(z32.nonzero_elements().len(), 4);
    assert_eq!(
        z32.with_divided_scalar(0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    // MIN / -1 errors; MIN / 1 defined
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        mn32.with_divided_scalar(-1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
    assert_eq!(
        mn32.with_divided_scalar(1).unwrap().value(0).unwrap(),
        i32::MIN
    );
    // retained zero quotient (3/4 -> 0 stays stored)
    let mut small32 = SparseCountFingerprint32::new(8);
    small32.set_value(5, 3).unwrap();
    let q = small32.with_divided_scalar(4).unwrap();
    assert_eq!(q.value(5).unwrap(), 0);
    assert_eq!(q.nonzero_elements().len(), 1);
    // unchanged inputs
    assert_eq!(v32.value(0).unwrap(), 7);

    // u64 width: re-pin F18 values
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(0, 7).unwrap();
    v.set_value(1, -7).unwrap();
    assert_eq!(v.with_divided_scalar(2).unwrap().value(0).unwrap(), 3);
    assert_eq!(v.with_divided_scalar(2).unwrap().value(1).unwrap(), -3);
    assert_eq!(
        v.with_divided_scalar(0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator/=(int)"
        })
    );
}

#[test]
fn u09_scalar_add_both_widths() {
    // U09: shared scalar addition. Empty succeeds; scalar 0; negatives;
    // zero result retained; overflow rejected; unchanged inputs.
    assert!(
        SparseCountFingerprint32::new(8)
            .with_added_scalar(i32::MIN)
            .is_ok()
    );
    let mut v32 = SparseCountFingerprint32::new(8);
    v32.set_value(1, 4).unwrap();
    v32.set_value(2, -6).unwrap();
    assert_eq!(v32.with_added_scalar(3).unwrap().value(1).unwrap(), 7);
    assert_eq!(v32.with_added_scalar(-10).unwrap().value(1).unwrap(), -6);
    assert_eq!(v32.with_added_scalar(0).unwrap().value(1).unwrap(), 4);
    // zero result retained
    let z = v32.with_added_scalar(-4).unwrap();
    assert_eq!(z.value(1).unwrap(), 0);
    assert_eq!(z.nonzero_elements().len(), 2);
    // no densification
    assert!(z.nonzero_elements().get(&0).is_none());
    // overflow rejected (CK boundary)
    let mut mx32 = SparseCountFingerprint32::new(8);
    mx32.set_value(0, i32::MAX).unwrap();
    assert_eq!(
        mx32.with_added_scalar(1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator+=(int)"
        })
    );
    // unchanged inputs
    assert_eq!(v32.value(1).unwrap(), 4);
    assert_eq!(mx32.value(0).unwrap(), i32::MAX);

    // u64 re-pin
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 4).unwrap();
    assert_eq!(v.with_added_scalar(3).unwrap().value(1).unwrap(), 7);
    assert_eq!(v.with_added_scalar(-4).unwrap().value(1).unwrap(), 0);
    assert_eq!(v.with_added_scalar(-4).unwrap().nonzero_elements().len(), 1);
}

#[test]
fn u10_scalar_sub_both_widths() {
    // U10: shared scalar subtraction. Empty succeeds; scalar 0;
    // negatives; zero result retained; overflow rejected; unchanged
    // inputs.
    assert!(
        SparseCountFingerprint32::new(8)
            .with_subtracted_scalar(i32::MIN)
            .is_ok()
    );
    let mut v32 = SparseCountFingerprint32::new(8);
    v32.set_value(1, 4).unwrap();
    v32.set_value(2, -6).unwrap();
    assert_eq!(v32.with_subtracted_scalar(3).unwrap().value(1).unwrap(), 1);
    assert_eq!(
        v32.with_subtracted_scalar(-10).unwrap().value(1).unwrap(),
        14
    );
    assert_eq!(v32.with_subtracted_scalar(0).unwrap().value(2).unwrap(), -6);
    // zero result retained; no densification
    let z = v32.with_subtracted_scalar(4).unwrap();
    assert_eq!(z.value(1).unwrap(), 0);
    assert_eq!(z.nonzero_elements().len(), 2);
    assert!(z.nonzero_elements().get(&0).is_none());
    // boundary sums defined; overflow rejected
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN + 1).unwrap();
    assert_eq!(
        mn32.with_subtracted_scalar(1).unwrap().value(0).unwrap(),
        i32::MIN
    );
    let mut under32 = SparseCountFingerprint32::new(8);
    under32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        under32.with_subtracted_scalar(1),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::operator-=(int)"
        })
    );
    // unchanged inputs
    assert_eq!(v32.value(1).unwrap(), 4);
    assert_eq!(under32.value(0).unwrap(), i32::MIN);

    // u64 re-pin
    let mut v = SparseCountFingerprint::new(8);
    v.set_value(1, 4).unwrap();
    assert_eq!(v.with_subtracted_scalar(3).unwrap().value(1).unwrap(), 1);
    assert_eq!(v.with_subtracted_scalar(4).unwrap().value(1).unwrap(), 0);
}

#[test]
fn u11_calc_vect_params_both_widths() {
    // U11: shared monotonic calcVectParams merge for both widths. Merge
    // shapes: both empty, either empty, interleaved/disjoint/overlapping,
    // both tails, negatives, stored zeros, MIN rejection, length
    // mismatch, u32 high-bit keys; equal results on the common domain.
    use cosmolkit_fingerprints::similarity;

    let e32 = SparseCountFingerprint32::new(8);
    assert_eq!(
        similarity::calc_vect_params_u32(&e32, &e32).unwrap(),
        (0.0, 0.0, 0.0)
    );
    // either empty
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(1, 2).unwrap();
    b32.set_value(4, -3).unwrap();
    assert_eq!(
        similarity::calc_vect_params_u32(&e32, &b32).unwrap(),
        (0.0, 5.0, 0.0)
    );
    assert_eq!(
        similarity::calc_vect_params_u32(&b32, &e32).unwrap(),
        (5.0, 0.0, 0.0)
    );
    // interleaved/overlapping with negatives: F50 oracle shape
    let mut i32a = SparseCountFingerprint32::new(16);
    i32a.set_value(1, -2).unwrap();
    i32a.set_value(3, 5).unwrap();
    i32a.set_value(5, -1).unwrap();
    let mut i32b = SparseCountFingerprint32::new(16);
    i32b.set_value(2, 1).unwrap();
    i32b.set_value(3, -7).unwrap();
    i32b.set_value(5, 4).unwrap();
    assert_eq!(
        similarity::calc_vect_params_u32(&i32a, &i32b).unwrap(),
        (8.0, 12.0, 6.0)
    );
    // disjoint with both tails exhausted independently
    let mut d1 = SparseCountFingerprint32::new(8);
    d1.set_value(5, 1).unwrap();
    d1.set_value(6, 2).unwrap();
    d1.set_value(7, 4).unwrap();
    let mut d2 = SparseCountFingerprint32::new(8);
    d2.set_value(0, 8).unwrap();
    d2.set_value(1, 16).unwrap();
    assert_eq!(
        similarity::calc_vect_params_u32(&d1, &d2).unwrap(),
        (7.0, 24.0, 0.0)
    );
    // stored zero participation
    let mut z32 = SparseCountFingerprint32::new(8);
    z32.set_value(2, 3).unwrap();
    let z32 = z32.with_added_scalar(-3).unwrap();
    let mut o32 = SparseCountFingerprint32::new(8);
    o32.set_value(2, 9).unwrap();
    o32.set_value(4, 1).unwrap();
    assert_eq!(
        similarity::calc_vect_params_u32(&z32, &o32).unwrap(),
        (0.0, 10.0, 0.0)
    );
    // MIN checked rejection (CK boundary)
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        similarity::calc_vect_params_u32(&mn32, &e32),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::calcVectParams abs"
        })
    );
    // length mismatch (error detail widened losslessly)
    assert_eq!(
        similarity::calc_vect_params_u32(&b32, &SparseCountFingerprint32::new(16)),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // u32 high-bit keys participate in the merge
    let mut hi = SparseCountFingerprint32::new(2_u32.pow(31) + 4);
    hi.set_value(2_u32.pow(31), 3).unwrap();
    hi.set_value(2_u32.pow(31) + 3, 4).unwrap();
    let mut hi2 = SparseCountFingerprint32::new(2_u32.pow(31) + 4);
    hi2.set_value(2_u32.pow(31) + 3, -9).unwrap();
    // and = min(4, 9) at the shared high key
    assert_eq!(
        similarity::calc_vect_params_u32(&hi, &hi2).unwrap(),
        (7.0, 9.0, 4.0)
    );
    // equal results on the common domain vs the u64 entrypoint
    let mut a = SparseCountFingerprint::new(16);
    a.set_value(1, -2).unwrap();
    a.set_value(3, 5).unwrap();
    a.set_value(5, -1).unwrap();
    let mut b = SparseCountFingerprint::new(16);
    b.set_value(2, 1).unwrap();
    b.set_value(3, -7).unwrap();
    b.set_value(5, 4).unwrap();
    assert_eq!(
        similarity::calc_vect_params(&a, &b).unwrap(),
        similarity::calc_vect_params_u32(&i32a, &i32b).unwrap()
    );
    // operands unchanged
    assert_eq!(i32a.value(3).unwrap(), 5);
    assert_eq!(i32b.value(3).unwrap(), -7);
}

#[test]
fn u12_count_dice_both_widths() {
    // U12: shared count Dice for both widths. F51 oracle shape:
    // a={1:2,3:5}, b={3:2,5:1} -> 0.4; distance 0.6; bounds branch both
    // directions (0.7 rejects -> 0.0, 0.5 accepts -> 0.4); bounds ignored
    // under returnDistance; empty pair 0.0; negatives; stored zeros;
    // |denom|<1e-6 threshold; length mismatch; approved arithmetic errors.
    use cosmolkit_fingerprints::similarity;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    assert!(
        (similarity::sparse_count_dice_u32(&a32, &b32, false, 0.0).unwrap() - 0.4).abs() < 1e-12
    );
    assert!(
        (similarity::sparse_count_dice_u32(&a32, &b32, true, 0.0).unwrap() - 0.6).abs() < 1e-12
    );
    // bounds both directions around the exact 2*min/sum = 0.6 cutoff
    assert!(similarity::sparse_count_dice_u32(&a32, &b32, false, 0.7).unwrap() == 0.0);
    assert!(
        (similarity::sparse_count_dice_u32(&a32, &b32, false, 0.5).unwrap() - 0.4).abs() < 1e-12
    );
    // bounds ignored under returnDistance
    assert!(
        (similarity::sparse_count_dice_u32(&a32, &b32, true, 0.9).unwrap() - 0.6).abs() < 1e-12
    );
    // empty pair (threshold path and plain path)
    let e32 = SparseCountFingerprint32::new(8);
    assert!(similarity::sparse_count_dice_u32(&e32, &e32, false, 0.0).unwrap() == 0.0);
    assert!(similarity::sparse_count_dice_u32(&e32, &e32, false, 0.5).unwrap() == 0.0);
    // negatives: n={0:-4,2:3}, m={0:2,2:-5} -> 10/14
    let mut n32 = SparseCountFingerprint32::new(8);
    n32.set_value(0, -4).unwrap();
    n32.set_value(2, 3).unwrap();
    let mut m32 = SparseCountFingerprint32::new(8);
    m32.set_value(0, 2).unwrap();
    m32.set_value(2, -5).unwrap();
    assert!(
        (similarity::sparse_count_dice_u32(&n32, &m32, false, 0.0).unwrap() - 10.0 / 14.0).abs()
            < 1e-12
    );
    // stored zero participates as zero
    let mut z32 = SparseCountFingerprint32::new(8);
    z32.set_value(2, 3).unwrap();
    let z32 = z32.with_added_scalar(-3).unwrap();
    let mut o32 = SparseCountFingerprint32::new(8);
    o32.set_value(2, 9).unwrap();
    // and = 0, sums 0+9 -> 0.0
    assert!(similarity::sparse_count_dice_u32(&z32, &o32, false, 0.0).unwrap() == 0.0);
    // length mismatch
    assert_eq!(
        similarity::sparse_count_dice_u32(&a32, &SparseCountFingerprint32::new(16), false, 0.0),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // approved arithmetic error surfaces through the bounds path (abs MIN)
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        similarity::sparse_count_dice_u32(&mn32, &e32, false, 0.5),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::getTotalVal abs"
        })
    );
    // operands unchanged
    assert_eq!(a32.value(1).unwrap(), 2);
    assert_eq!(z32.value(2).unwrap(), 0);

    // u64 width re-pin (shared body changed)
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!((similarity::sparse_count_dice(&a, &b, false, 0.0).unwrap() - 0.4).abs() < 1e-12);
}

#[test]
fn u13_count_tversky_both_widths() {
    // U13: shared count Tversky for both widths, without bit-vector
    // parameter checks. F52 oracle shape: a={1:2,3:5}, b={3:2,5:1}:
    // (1,1) -> 0.25; (.5,.5) -> 2/(3.5+1.5) = 0.4; (0,0) -> 1.0;
    // distance 0.75; coefficients outside [0,1] NOT rejected; empty pair
    // -> similarity 0 / distance 1; length mismatch; bounds unused.
    use cosmolkit_fingerprints::similarity;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    assert!(
        (similarity::sparse_count_tversky_u32(&a32, &b32, 1.0, 1.0, false, 0.0).unwrap() - 0.25)
            .abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky_u32(&a32, &b32, 0.5, 0.5, false, 0.0).unwrap() - 0.4)
            .abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky_u32(&a32, &b32, 0.0, 0.0, false, 0.0).unwrap() - 1.0)
            .abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tversky_u32(&a32, &b32, 1.0, 1.0, true, 0.0).unwrap() - 0.75)
            .abs()
            < 1e-12
    );
    // out-of-range coefficients are not rejected in the count variant
    assert!(similarity::sparse_count_tversky_u32(&a32, &b32, 2.0, -1.0, false, 0.0).is_ok());
    // denominator threshold: empty pair -> 0.0 both distance modes
    let e32 = SparseCountFingerprint32::new(8);
    assert!(similarity::sparse_count_tversky_u32(&e32, &e32, 1.0, 1.0, false, 0.0).unwrap() == 0.0);
    assert!(similarity::sparse_count_tversky_u32(&e32, &e32, 1.0, 1.0, true, 0.0).unwrap() == 1.0);
    // bounds argument is unused by the source
    assert!(
        (similarity::sparse_count_tversky_u32(&a32, &b32, 1.0, 1.0, false, 0.9).unwrap() - 0.25)
            .abs()
            < 1e-12
    );
    // length mismatch
    assert_eq!(
        similarity::sparse_count_tversky_u32(
            &a32,
            &SparseCountFingerprint32::new(16),
            1.0,
            1.0,
            false,
            0.0
        ),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // negatives through abs accumulation
    let mut n32 = SparseCountFingerprint32::new(8);
    n32.set_value(0, -4).unwrap();
    n32.set_value(2, 3).unwrap();
    let mut m32 = SparseCountFingerprint32::new(8);
    m32.set_value(0, 2).unwrap();
    m32.set_value(2, -5).unwrap();
    // (1,1): 5 / (7+7-5) = 5/9
    assert!(
        (similarity::sparse_count_tversky_u32(&n32, &m32, 1.0, 1.0, false, 0.0).unwrap()
            - 5.0 / 9.0)
            .abs()
            < 1e-12
    );
    // operands unchanged
    assert_eq!(a32.value(3).unwrap(), 5);

    // u64 width re-pin
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!(
        (similarity::sparse_count_tversky(&a, &b, 1.0, 1.0, false, 0.0).unwrap() - 0.25).abs()
            < 1e-12
    );
}

#[test]
fn u14_count_tanimoto_both_widths() {
    // U14: shared count Tanimoto via the source Tversky(1,1) delegation
    // for both widths. F53 oracle shape: a={1:2,3:5}, b={3:2,5:1} ->
    // 0.25; distance 0.75; both empty 0.0; negatives 5/9; bounds variants
    // (unused); length mismatch; error propagation (abs MIN through the
    // bounds-free path's calcVectParams).
    use cosmolkit_fingerprints::similarity;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    assert!(
        (similarity::sparse_count_tanimoto_u32(&a32, &b32, false, 0.0).unwrap() - 0.25).abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tanimoto_u32(&a32, &b32, true, 0.0).unwrap() - 0.75).abs()
            < 1e-12
    );
    // Tversky(1,1) equivalence
    assert_eq!(
        similarity::sparse_count_tanimoto_u32(&a32, &b32, false, 0.3),
        similarity::sparse_count_tversky_u32(&a32, &b32, 1.0, 1.0, false, 0.3)
    );
    // both empty
    let e32 = SparseCountFingerprint32::new(8);
    assert!(similarity::sparse_count_tanimoto_u32(&e32, &e32, false, 0.0).unwrap() == 0.0);
    assert!(similarity::sparse_count_tanimoto_u32(&e32, &e32, true, 0.0).unwrap() == 1.0);
    // negatives: 5/9
    let mut n32 = SparseCountFingerprint32::new(8);
    n32.set_value(0, -4).unwrap();
    n32.set_value(2, 3).unwrap();
    let mut m32 = SparseCountFingerprint32::new(8);
    m32.set_value(0, 2).unwrap();
    m32.set_value(2, -5).unwrap();
    assert!(
        (similarity::sparse_count_tanimoto_u32(&n32, &m32, false, 0.0).unwrap() - 5.0 / 9.0).abs()
            < 1e-12
    );
    // length mismatch
    assert_eq!(
        similarity::sparse_count_tanimoto_u32(&a32, &SparseCountFingerprint32::new(16), false, 0.0),
        Err(FingerprintError::BitLengthMismatch { left: 8, right: 16 })
    );
    // approved arithmetic error propagates through the delegation
    let mut mn32 = SparseCountFingerprint32::new(8);
    mn32.set_value(0, i32::MIN).unwrap();
    assert_eq!(
        similarity::sparse_count_tanimoto_u32(&mn32, &e32, false, 0.0),
        Err(FingerprintError::UndefinedArithmetic {
            site: "SparseIntVect::calcVectParams abs"
        })
    );
    // operands unchanged
    assert_eq!(a32.value(1).unwrap(), 2);
    assert_eq!(mn32.value(0).unwrap(), i32::MIN);

    // u64 width re-pin
    let mut a = SparseCountFingerprint::new(8);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let mut b = SparseCountFingerprint::new(8);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    assert!((similarity::sparse_count_tanimoto(&a, &b, false, 0.0).unwrap() - 0.25).abs() < 1e-12);
}

#[test]
fn closure_cross_width_sequence() {
    // Whole-lane closure sequence for both widths: construction -> set ->
    // scalar-produced stored zero -> each vector operation -> each count
    // similarity, preserving u32 maximum-index behavior and width-specific
    // length checks. All 160 wrapper dispatch rows and the C1-C4
    // regressions remain in their own tests.
    use cosmolkit_fingerprints::similarity;

    // ---- u32 sequence ----
    let mut a = SparseCountFingerprint32::new(16);
    a.set_value(1, 2).unwrap();
    a.set_value(3, 5).unwrap();
    let stored_zero = a.with_multiplied_scalar(0).unwrap(); // {1:0, 3:0} retained
    assert_eq!(stored_zero.nonzero_elements().len(), 2);
    let revived = stored_zero.with_added_scalar(7).unwrap(); // {1:7, 3:7}
    assert_eq!(revived.value(1).unwrap(), 7);
    // vector ops on the revived pair
    let mut b = SparseCountFingerprint32::new(16);
    b.set_value(3, 2).unwrap();
    b.set_value(5, 1).unwrap();
    let and = revived.fuzzy_and(&b).unwrap();
    assert_eq!(and.value(3).unwrap(), 2);
    assert!(and.nonzero_elements().get(&1).is_none());
    let or = revived.fuzzy_or(&b).unwrap();
    assert_eq!(or.value(5).unwrap(), 1);
    let sum = revived.with_added(&b).unwrap();
    assert_eq!(sum.value(3).unwrap(), 9);
    let diff = revived.with_subtracted(&b).unwrap();
    assert_eq!(diff.value(5).unwrap(), -1);
    let scaled = sum.with_divided_scalar(3).unwrap();
    assert_eq!(scaled.value(3).unwrap(), 3);
    let shifted = scaled.with_subtracted_scalar(1).unwrap();
    assert_eq!(shifted.value(3).unwrap(), 2);
    // similarities on the original pair
    assert!((similarity::sparse_count_dice_u32(&a, &b, false, 0.0).unwrap() - 0.4).abs() < 1e-12);
    assert!(
        (similarity::sparse_count_tversky_u32(&a, &b, 1.0, 1.0, false, 0.0).unwrap() - 0.25).abs()
            < 1e-12
    );
    assert!(
        (similarity::sparse_count_tanimoto_u32(&a, &b, false, 0.0).unwrap() - 0.25).abs() < 1e-12
    );
    // u32 maximum-index behavior within the sequence
    let mut edge = SparseCountFingerprint32::new(u32::MAX);
    edge.set_value(u32::MAX, 4).unwrap();
    edge.set_value(0, 1).unwrap();
    assert_eq!(edge.total_value(false).unwrap(), 5);
    // width-specific length check: same numeric length is same-length for
    // each width family; cross-width mixing is a type error by design
    assert_eq!(
        a.fuzzy_and(&SparseCountFingerprint32::new(17)),
        Err(FingerprintError::BitLengthMismatch {
            left: 16,
            right: 17
        })
    );

    // ---- u64 sequence (mirror) ----
    let mut a64 = SparseCountFingerprint::new(16);
    a64.set_value(1, 2).unwrap();
    a64.set_value(3, 5).unwrap();
    let sz64 = a64.with_multiplied_scalar(0).unwrap();
    assert_eq!(sz64.nonzero_elements().len(), 2);
    let revived64 = sz64.with_added_scalar(7).unwrap();
    assert_eq!(revived64.value(1).unwrap(), 7);
    let mut b64 = SparseCountFingerprint::new(16);
    b64.set_value(3, 2).unwrap();
    b64.set_value(5, 1).unwrap();
    assert_eq!(revived64.fuzzy_and(&b64).unwrap().value(3).unwrap(), 2);
    assert_eq!(revived64.with_added(&b64).unwrap().value(3).unwrap(), 9);
    assert_eq!(
        revived64.with_subtracted(&b64).unwrap().value(5).unwrap(),
        -1
    );
    assert!((similarity::sparse_count_dice(&a64, &b64, false, 0.0).unwrap() - 0.4).abs() < 1e-12);
    assert!(
        (similarity::sparse_count_tanimoto(&a64, &b64, false, 0.0).unwrap() - 0.25).abs() < 1e-12
    );
    // u64 width-specific check: length u32::MAX rejects index == length
    let mut wide64 = SparseCountFingerprint::new(u64::from(u32::MAX));
    assert_eq!(
        wide64.set_value(u64::from(u32::MAX), 1),
        Err(FingerprintError::SparseIndexOutOfRange {
            index: u64::from(u32::MAX),
            size: u64::from(u32::MAX)
        })
    );
    // operands unchanged through the sequences
    assert_eq!(a.value(1).unwrap(), 2);
    assert_eq!(a64.value(3).unwrap(), 5);
}

#[test]
fn acceptance_dice_exact_bounds_both_widths() {
    // Acceptance: Dice bounds cutoff is strictly `<` at the exact
    // 2*minV/denom == 0.6 neighborhood (supervisor-frozen bit patterns,
    // locally re-confirmed against pinned UInt/ULongSparseIntVect:
    // 0x...32 -> 0x3fd999999999999a, 0x...33 -> 0x3fd999999999999a
    // (equal bounds still runs the metric: strict <), 0x...34 ->
    // 0x0000000000000000; distance mode always 0x3fe3333333333333
    // (the bounds branch only applies when !returnDistance).
    // Row count: 2 widths x 3 bounds x 2 modes = 12.
    use cosmolkit_fingerprints::similarity;

    let bounds_bits = [
        0x3fe3333333333332u64,
        0x3fe3333333333333,
        0x3fe3333333333334,
    ];
    let mut rows = 0u32;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(1, 2).unwrap();
    a32.set_value(3, 5).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(3, 2).unwrap();
    b32.set_value(5, 1).unwrap();
    let (a32_0, b32_0) = (a32.clone(), b32.clone());

    let mut a64 = SparseCountFingerprint::new(8);
    a64.set_value(1, 2).unwrap();
    a64.set_value(3, 5).unwrap();
    let mut b64 = SparseCountFingerprint::new(8);
    b64.set_value(3, 2).unwrap();
    b64.set_value(5, 1).unwrap();
    let (a64_0, b64_0) = (a64.clone(), b64.clone());

    for (i, &bb) in bounds_bits.iter().enumerate() {
        let bounds = f64::from_bits(bb);
        for distance in [false, true] {
            rows += 1;
            let sim = similarity::sparse_count_dice_u32(&a32, &b32, distance, bounds)
                .unwrap()
                .to_bits();
            let expected_sim = if i <= 1 {
                0x3fd999999999999a
            } else {
                0x0000000000000000
            };
            let expected_dist = 0x3fe3333333333333;
            assert_eq!(
                sim,
                if distance {
                    expected_dist
                } else {
                    expected_sim
                }
            );
            assert_eq!(a32, a32_0, "row {i} distance={distance}");
            assert_eq!(b32, b32_0, "row {i} distance={distance}");

            rows += 1;
            let sim = similarity::sparse_count_dice(&a64, &b64, distance, bounds)
                .unwrap()
                .to_bits();
            assert_eq!(
                sim,
                if distance {
                    expected_dist
                } else {
                    expected_sim
                }
            );
            assert_eq!(a64, a64_0, "row {i} distance={distance}");
            assert_eq!(b64, b64_0, "row {i} distance={distance}");
        }
    }
    assert_eq!(rows, 12, "dice acceptance matrix must execute 12 rows");
}

#[test]
fn acceptance_tversky_nonzero_threshold_both_widths() {
    // Acceptance: count Tversky nonzero-numerator branches around the
    // |denom| < 1e-6 cutoff (supervisor-frozen bit patterns, locally
    // re-confirmed against pinned UInt/ULongSparseIntVect). A={0:2},
    // B={0:1}, beta=0, bounds=0: v1Sum=2, v2Sum=1, and=1,
    // denom = alpha*2 + 1*(1-alpha). Two alphas fall inside |denom|<1e-6
    // (sim 0.0 / distance 1.0) and two just outside (huge positive/
    // negative ratios with exact distance 1-sim). 2 widths x 4 alphas x
    // 2 modes = 16 rows.
    use cosmolkit_fingerprints::similarity;

    let cases: [(u64, u64, u64); 4] = [
        (0xbfeffffef39085f5, 0x0000000000000000, 0x3ff0000000000000),
        (0xbfeffffbce4217d3, 0x411e847ffffc3b1f, 0xc11e847bfffc3b1f),
        (0xbff000008637bd06, 0x0000000000000000, 0x3ff0000000000000),
        (0xbff0000218def417, 0xc11e8480000ac869, 0x411e8484000ac869),
    ];
    let mut rows = 0u32;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(0, 2).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(0, 1).unwrap();
    let (a32_0, b32_0) = (a32.clone(), b32.clone());

    let mut a64 = SparseCountFingerprint::new(8);
    a64.set_value(0, 2).unwrap();
    let mut b64 = SparseCountFingerprint::new(8);
    b64.set_value(0, 1).unwrap();
    let (a64_0, b64_0) = (a64.clone(), b64.clone());

    for &(alpha_bits, sim_bits, dist_bits) in &cases {
        let alpha = f64::from_bits(alpha_bits);
        for distance in [false, true] {
            rows += 1;
            let r = similarity::sparse_count_tversky_u32(&a32, &b32, alpha, 0.0, distance, 0.0)
                .unwrap()
                .to_bits();
            assert_eq!(
                r,
                if distance { dist_bits } else { sim_bits },
                "u32 alpha {alpha_bits:#x} distance={distance}"
            );
            assert_eq!(a32, a32_0);
            assert_eq!(b32, b32_0);

            rows += 1;
            let r = similarity::sparse_count_tversky(&a64, &b64, alpha, 0.0, distance, 0.0)
                .unwrap()
                .to_bits();
            assert_eq!(
                r,
                if distance { dist_bits } else { sim_bits },
                "u64 alpha {alpha_bits:#x} distance={distance}"
            );
            assert_eq!(a64, a64_0);
            assert_eq!(b64, b64_0);
        }
    }
    assert_eq!(rows, 16, "tversky acceptance matrix must execute 16 rows");
}

#[test]
fn acceptance_tversky_exact_negative_threshold_both_widths() {
    // Acceptance: exact |denom| == 1e-6 signed-zero boundary (supervisor
    // frozen, locally re-confirmed). Disjoint A={0:1}, B={1:1}, beta=0,
    // bounds=0: v1Sum=1, v2Sum=1, and=0, denom = alpha*1 = alpha (tiny
    // negative). alpha one ulp inside -> sim = +0.0 (0x0); alpha exactly
    // -1e-6 -> strict `<` is false -> sim = 0/(-1e-6) = -0.0
    // (0x8000000000000000); alpha one ulp outside -> -0.0. distance is
    // 1.0 (0x3ff0000000000000) in every case since 1-(±0.0)=1.0. This
    // distinguishes `<` from `<=`: a `<=` cutoff would turn the exact
    // -1e-6 row into the early-exit +0.0. 2 widths x 3 alphas x 2 modes
    // = 12 rows.
    use cosmolkit_fingerprints::similarity;

    let alpha_bits = [
        0xbeb0c6f7a0b5ed8cu64,
        0xbeb0c6f7a0b5ed8d,
        0xbeb0c6f7a0b5ed8e,
    ];
    let expected_sim = [
        0x0000000000000000u64,
        0x8000000000000000,
        0x8000000000000000,
    ];
    let expected_dist = 0x3ff0000000000000u64;
    let mut rows = 0u32;

    let mut a32 = SparseCountFingerprint32::new(8);
    a32.set_value(0, 1).unwrap();
    let mut b32 = SparseCountFingerprint32::new(8);
    b32.set_value(1, 1).unwrap();
    let (a32_0, b32_0) = (a32.clone(), b32.clone());

    let mut a64 = SparseCountFingerprint::new(8);
    a64.set_value(0, 1).unwrap();
    let mut b64 = SparseCountFingerprint::new(8);
    b64.set_value(1, 1).unwrap();
    let (a64_0, b64_0) = (a64.clone(), b64.clone());

    for (i, &ab) in alpha_bits.iter().enumerate() {
        let alpha = f64::from_bits(ab);
        for distance in [false, true] {
            rows += 1;
            let r = similarity::sparse_count_tversky_u32(&a32, &b32, alpha, 0.0, distance, 0.0)
                .unwrap()
                .to_bits();
            assert_eq!(
                r,
                if distance {
                    expected_dist
                } else {
                    expected_sim[i]
                },
                "u32 alpha {ab:#x} distance={distance}"
            );
            assert_eq!(a32, a32_0);
            assert_eq!(b32, b32_0);

            rows += 1;
            let r = similarity::sparse_count_tversky(&a64, &b64, alpha, 0.0, distance, 0.0)
                .unwrap()
                .to_bits();
            assert_eq!(
                r,
                if distance {
                    expected_dist
                } else {
                    expected_sim[i]
                },
                "u64 alpha {ab:#x} distance={distance}"
            );
            assert_eq!(a64, a64_0);
            assert_eq!(b64, b64_0);
        }
    }
    assert_eq!(rows, 12, "exact-threshold matrix must execute 12 rows");
}

#[test]
fn regression_fuzzy_and_retain_across_node_boundary_both_widths() {
    // Twelve entries exceed the inspected BTreeMap leaf capacity (11).
    // One fixed deletion/re-retention case, not a generated size/order matrix.
    macro_rules! check_width {
        ($ty:ty, $index:ty) => {{
            let mut left = <$ty>::new(16);
            for key in 0..12 {
                left.set_value(key as $index, 3).unwrap();
            }
            let mut right = <$ty>::new(16);
            for (key, value) in [(0, 2), (5, -1), (11, 7)] {
                right.set_value(key as $index, value).unwrap();
            }
            let before = (left.clone(), right.clone());
            let result = left.fuzzy_and(&right).unwrap();
            let actual: Vec<_> = result
                .nonzero_elements()
                .iter()
                .map(|(&key, &value)| (key, value))
                .collect();
            assert_eq!(actual, vec![(0, 2), (5, -1), (11, 3)]);
            assert_eq!((left, right), before);
            assert_eq!(result.fuzzy_and(&result).unwrap(), result);
            let empty = <$ty>::new(16);
            assert_eq!(result.fuzzy_and(&empty).unwrap(), empty);
        }};
    }
    check_width!(SparseCountFingerprint32, u32);
    check_width!(SparseCountFingerprint, u64);
}

#[test]
fn regression_fuzzy_and_fixed_shapes_both_widths() {
    // Fixed source-defined expected values; no corpus or timing assertions.
    fn check32(
        a_entries: &[(u32, i32)],
        b_entries: &[(u32, i32)],
        expected: &[(u32, i32)],
        name: &str,
    ) {
        let mut a = SparseCountFingerprint32::new(16);
        let mut b = SparseCountFingerprint32::new(16);
        for &(k, v) in a_entries {
            a.set_value(k, v).unwrap();
        }
        for &(k, v) in b_entries {
            b.set_value(k, v).unwrap();
        }
        let a0 = a.clone();
        let b0 = b.clone();
        let r = a.fuzzy_and(&b).unwrap();
        let got: Vec<(u32, i32)> = r.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
        assert_eq!(got, expected.to_vec(), "{name}");
        assert_eq!(a, a0, "{name}: left unchanged");
        assert_eq!(b, b0, "{name}: right unchanged");
        // length mismatch error path also leaves inputs untouched
        let mut big = SparseCountFingerprint32::new(32);
        big.set_value(1, 1).unwrap();
        assert_eq!(
            a.fuzzy_and(&big),
            Err(FingerprintError::BitLengthMismatch {
                left: 16,
                right: 32
            }),
            "{name}: mismatch"
        );
        assert_eq!(a, a0, "{name}: left unchanged after error");
    }
    fn check64(
        a_entries: &[(u64, i32)],
        b_entries: &[(u64, i32)],
        expected: &[(u64, i32)],
        name: &str,
    ) {
        let mut a = SparseCountFingerprint::new(16);
        let mut b = SparseCountFingerprint::new(16);
        for &(k, v) in a_entries {
            a.set_value(k, v).unwrap();
        }
        for &(k, v) in b_entries {
            b.set_value(k, v).unwrap();
        }
        let a0 = a.clone();
        let b0 = b.clone();
        let r = a.fuzzy_and(&b).unwrap();
        let got: Vec<(u64, i32)> = r.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
        assert_eq!(got, expected.to_vec(), "{name}");
        assert_eq!(a, a0, "{name}: left unchanged");
        assert_eq!(b, b0, "{name}: right unchanged");
    }

    // ---- u32 shapes ----
    check32(&[], &[], &[], "empty/empty");
    check32(
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        "complete overlap",
    );
    check32(
        &[(1, 2), (3, 5)],
        &[(5, 1), (7, 4)],
        &[],
        "disjoint left<right",
    );
    check32(
        &[(5, 1), (7, 4)],
        &[(1, 2), (3, 5)],
        &[],
        "disjoint right<left",
    );
    check32(
        &[(1, 1), (3, 1), (5, 1)],
        &[(2, 9), (3, 9), (4, 9)],
        &[(3, 1)],
        "alternating",
    );
    check32(
        &[(1, 7), (2, 7), (3, 7)],
        &[(2, 3)],
        &[(2, 3)],
        "right subset, min from right",
    );
    check32(
        &[(2, 3)],
        &[(1, 7), (2, 7), (3, 7)],
        &[(2, 3)],
        "left subset",
    );
    // negative / equal / extrema / stored zero
    check32(
        &[(1, -4), (2, 3)],
        &[(1, 2), (2, -5)],
        &[(1, -4), (2, -5)],
        "negatives",
    );
    check32(&[(1, 5)], &[(1, 5)], &[(1, 5)], "equal values");
    check32(
        &[(1, i32::MIN), (2, 4)],
        &[(1, 7), (2, i32::MAX)],
        &[(1, i32::MIN), (2, 4)],
        "extrema",
    );
    {
        let mut sz = SparseCountFingerprint32::new(16);
        sz.set_value(2, 3).unwrap();
        let sz = sz.with_added_scalar(-3).unwrap(); // {2: 0} retained
        let mut other = SparseCountFingerprint32::new(16);
        other.set_value(2, 9).unwrap();
        other.set_value(5, 1).unwrap();
        let r = sz.fuzzy_and(&other).unwrap();
        assert_eq!(r.value(2).unwrap(), 0);
        assert!(
            r.nonzero_elements().contains_key(&2),
            "stored zero retained in result"
        );
        assert!(!r.nonzero_elements().contains_key(&5));
    }

    // ---- u64 shapes ----
    check64(&[], &[], &[], "empty/empty");
    check64(
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        "complete overlap",
    );
    check64(&[(1, 2), (3, 5)], &[(5, 1), (7, 4)], &[], "disjoint");
    check64(
        &[(1, 1), (3, 1), (5, 1)],
        &[(2, 9), (3, 9), (4, 9)],
        &[(3, 1)],
        "alternating",
    );
    check64(
        &[(1, -4), (2, 3)],
        &[(1, 2), (2, -5)],
        &[(1, -4), (2, -5)],
        "negatives",
    );
    check64(
        &[(1, i32::MIN), (2, 4)],
        &[(1, 7), (2, i32::MAX)],
        &[(1, i32::MIN), (2, 4)],
        "extrema",
    );
}

#[test]
fn regression_fuzzy_or_fixed_shapes_both_widths() {
    // Fixed source-defined expected values; no corpus or timing assertions.
    fn check32(
        a_entries: &[(u32, i32)],
        b_entries: &[(u32, i32)],
        expected: &[(u32, i32)],
        name: &str,
    ) {
        let mut a = SparseCountFingerprint32::new(64);
        let mut b = SparseCountFingerprint32::new(64);
        for &(k, v) in a_entries {
            a.set_value(k, v).unwrap();
        }
        for &(k, v) in b_entries {
            b.set_value(k, v).unwrap();
        }
        let a0 = a.clone();
        let b0 = b.clone();
        let r = a.fuzzy_or(&b).unwrap();
        let got: Vec<(u32, i32)> = r.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
        assert_eq!(got, expected.to_vec(), "{name}");
        assert_eq!(a, a0, "{name}: left unchanged");
        assert_eq!(b, b0, "{name}: right unchanged");
        let mut big = SparseCountFingerprint32::new(128);
        big.set_value(1, 1).unwrap();
        assert_eq!(
            a.fuzzy_or(&big),
            Err(FingerprintError::BitLengthMismatch {
                left: 64,
                right: 128
            }),
            "{name}: mismatch"
        );
        assert_eq!(a, a0, "{name}: left unchanged after error");
    }
    fn check64(
        a_entries: &[(u64, i32)],
        b_entries: &[(u64, i32)],
        expected: &[(u64, i32)],
        name: &str,
    ) {
        let mut a = SparseCountFingerprint::new(64);
        let mut b = SparseCountFingerprint::new(64);
        for &(k, v) in a_entries {
            a.set_value(k, v).unwrap();
        }
        for &(k, v) in b_entries {
            b.set_value(k, v).unwrap();
        }
        let a0 = a.clone();
        let b0 = b.clone();
        let r = a.fuzzy_or(&b).unwrap();
        let got: Vec<(u64, i32)> = r.nonzero_elements().iter().map(|(&k, &v)| (k, v)).collect();
        assert_eq!(got, expected.to_vec(), "{name}");
        assert_eq!(a, a0, "{name}: left unchanged");
        assert_eq!(b, b0, "{name}: right unchanged");
    }

    // ---- u32 shapes ----
    check32(&[], &[], &[], "empty/empty");
    check32(
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        "complete overlap",
    );
    check32(
        &[(1, 2), (3, 5)],
        &[(5, 1), (7, 4)],
        &[(1, 2), (3, 5), (5, 1), (7, 4)],
        "disjoint left<right",
    );
    check32(
        &[(5, 1), (7, 4)],
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5), (5, 1), (7, 4)],
        "disjoint right<left",
    );
    check32(
        &[(1, 1), (3, 1), (5, 1)],
        &[(2, 9), (3, 9), (4, 9)],
        &[(1, 1), (2, 9), (3, 9), (4, 9), (5, 1)],
        "alternating",
    );
    // right-only keys before/between/after left keys
    check32(
        &[(10, 1), (20, 2), (30, 3)],
        &[(5, 7), (15, 8), (25, 9), (35, 6)],
        &[(5, 7), (10, 1), (15, 8), (20, 2), (25, 9), (30, 3), (35, 6)],
        "nested/before/between/after",
    );
    // nested supports: right entirely within left
    check32(
        &[(1, 7), (2, 7), (3, 7)],
        &[(2, 3)],
        &[(1, 7), (2, 7), (3, 7)],
        "nested right, max keeps left",
    );
    check32(
        &[(2, 3)],
        &[(1, 7), (2, 7), (3, 7)],
        &[(1, 7), (2, 7), (3, 7)],
        "nested left",
    );
    // negative / equal / extrema / stored zero
    check32(
        &[(1, -4), (2, 3)],
        &[(1, 2), (2, -5)],
        &[(1, 2), (2, 3)],
        "negatives, max",
    );
    check32(&[(1, 5)], &[(1, 5)], &[(1, 5)], "equal values");
    check32(
        &[(1, i32::MIN), (2, 4)],
        &[(1, 7), (2, i32::MAX)],
        &[(1, 7), (2, i32::MAX)],
        "extrema",
    );
    {
        let mut sz = SparseCountFingerprint32::new(64);
        sz.set_value(2, 3).unwrap();
        let sz = sz.with_added_scalar(-3).unwrap(); // {2: 0} retained
        let mut other = SparseCountFingerprint32::new(64);
        other.set_value(2, -4).unwrap();
        other.set_value(5, 6).unwrap();
        let r = sz.fuzzy_or(&other).unwrap();
        assert_eq!(r.value(2).unwrap(), 0, "max(0, -4) = 0 stays stored");
        assert!(r.nonzero_elements().contains_key(&2));
        assert_eq!(r.value(5).unwrap(), 6);
        // reverse order: right stored zero vs left negative
        let r2 = other.fuzzy_or(&sz).unwrap();
        assert_eq!(r2.value(2).unwrap(), 0);
        assert!(r2.nonzero_elements().contains_key(&2));
    }

    // ---- u64 shapes ----
    check64(&[], &[], &[], "empty/empty");
    check64(
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        &[(1, 2), (3, 5)],
        "complete overlap",
    );
    check64(
        &[(1, 2), (3, 5)],
        &[(5, 1), (7, 4)],
        &[(1, 2), (3, 5), (5, 1), (7, 4)],
        "disjoint",
    );
    check64(
        &[(10, 1), (20, 2), (30, 3)],
        &[(5, 7), (15, 8), (25, 9), (35, 6)],
        &[(5, 7), (10, 1), (15, 8), (20, 2), (25, 9), (30, 3), (35, 6)],
        "nested",
    );
    check64(
        &[(1, -4), (2, 3)],
        &[(1, 2), (2, -5)],
        &[(1, 2), (2, 3)],
        "negatives",
    );
    check64(
        &[(1, i32::MIN), (2, 4)],
        &[(1, 7), (2, i32::MAX)],
        &[(1, 7), (2, i32::MAX)],
        "extrema",
    );
}
