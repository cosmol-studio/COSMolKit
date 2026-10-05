//! Fixed explicit fingerprint construction and similarity value semantics.
#![cfg(feature = "cap-fingerprints")]
use cosmolkit::{BINDING_CONTRACT, Fingerprint, FingerprintError};
#[test]
fn constructor_sorts_deduplicates_and_retains_boundary_bits() {
    let source = vec![64, 0, 63, 64, 1];
    let value = Fingerprint::from_on_bits(65, source.clone()).unwrap();
    assert_eq!(value.n_bits(), 65);
    assert_eq!(value.on_bits(), [0, 1, 63, 64]);
    assert_eq!(source, [64, 0, 63, 64, 1]);
    let mut detached = value.on_bits();
    detached.clear();
    assert_eq!(value.on_bits(), [0, 1, 63, 64]);
}
#[test]
fn empty_width_and_first_invalid_index_validate_before_allocation() {
    let empty = Fingerprint::from_on_bits(0, []).unwrap();
    assert_eq!(empty.n_bits(), 0);
    assert!(empty.on_bits().is_empty());
    for (width, bits, index) in [
        (0, vec![0], 0),
        (7, vec![7, 8], 7),
        (u32::MAX, vec![u32::MAX], u32::MAX),
    ] {
        assert_eq!(
            Fingerprint::from_on_bits(width, bits).unwrap_err(),
            FingerprintError::SparseIndexOutOfRange {
                index: u64::from(index),
                size: u64::from(width)
            }
        );
    }
}
#[test]
fn exact_tanimoto_empty_intersection_and_equal_width_precondition() {
    let a = Fingerprint::from_on_bits(65, [0, 1, 64]).unwrap();
    let b = Fingerprint::from_on_bits(65, [1, 63, 64]).unwrap();
    assert_eq!(a.tanimoto(&b).unwrap().to_bits(), 0.5f64.to_bits());
    assert_eq!(b.tanimoto(&a).unwrap(), 0.5);
    assert_eq!(a.tanimoto(&a).unwrap(), 1.0);
    let zero = Fingerprint::from_on_bits(65, []).unwrap();
    assert_eq!(zero.tanimoto(&zero).unwrap(), 0.0);
    assert_eq!(a.tanimoto(&zero).unwrap(), 0.0);
    let mismatched = Fingerprint::from_on_bits(64, []).unwrap();
    assert_eq!(
        zero.tanimoto(&mismatched).unwrap_err(),
        FingerprintError::BitLengthMismatch {
            left: 65,
            right: 64
        }
    );
}
#[test]
fn real_constructor_and_similarity_contracts_are_declared() {
    for id in ["Fingerprint.from_on_bits", "Fingerprint.tanimoto"] {
        let row = BINDING_CONTRACT
            .iter()
            .find(|row| row.semantic_id == id)
            .unwrap();
        assert!(row.callable.is_some());
        assert_eq!(row.feature, "cap-fingerprints");
    }
}
