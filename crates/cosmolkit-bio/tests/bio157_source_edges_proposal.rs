//! Source boundary proposals; p1 review and ROOT acceptance remain pending.
use cosmolkit_bio::{
    BioChainId, BioResidueRow, BioRowSpan, BioSiftsUnpResidue, BioStructureError, BioTransform,
    EntityKind, PdbSeqId, ResidueInfoKind, ResidueName, ResidueSourceIds, auth_seq_id_to_label,
    label_seq_id_to_auth,
};

fn row(auth: i32, insertion: Option<u8>, label: Option<i32>) -> BioResidueRow {
    BioResidueRow::new(
        BioChainId::new(0),
        BioRowSpan::new(0, 0).unwrap(),
        ResidueName::from_ascii(b"ALA").unwrap(),
        ResidueInfoKind::Aa,
        EntityKind::Polymer,
        None,
        None,
        ResidueSourceIds::new(
            Some(PdbSeqId::new(auth, insertion)),
            label,
            None,
            None,
            None,
        )
        .unwrap(),
        BioSiftsUnpResidue::default(),
    )
}

#[test]
fn transform_nan_matrix_and_vector_follow_different_source_comparisons() {
    // Pinned Gemmi math.hpp: Mat33::approx uses fabs(delta)>epsilon;
    // Vec3::approx uses fabs(delta)<=epsilon. NaN makes both predicates false.
    let identity = BioTransform::identity();
    let mut matrix = *identity.matrix();
    matrix[0][0] = f64::NAN;
    assert!(BioTransform::new(matrix, [0.; 3]).approx(&identity, 0.));
    assert!(!BioTransform::new(*identity.matrix(), [f64::NAN, 0., 0.]).approx(&identity, 0.));
    assert!(!identity.approx(&identity, f64::NAN));
}

#[test]
fn transform_source_tolerance_is_componentwise_and_inclusive() {
    let identity = BioTransform::identity();
    // Each component equals epsilon: Vec3 length exceeds epsilon, but source
    // compares the three absolute components separately.
    let other = BioTransform::new(*identity.matrix(), [0.5; 3]);
    assert!(other.approx(&identity, 0.5));
    assert!(!other.approx(&identity, 0.499));
    assert!(!identity.approx(&identity, -1.));
}

#[test]
fn sequence_exact_match_preserves_insertion_and_case_insensitive_auth_identity() {
    // Pinned Gemmi model.hpp span exact-match paths and seqid.hpp SeqId equality.
    let a = row(10, Some(b'A'), Some(1));
    let b = row(10, Some(b'B'), Some(2));
    let c = row(20, None, Some(8));
    let span = [&a, &b, &c];
    assert_eq!(
        label_seq_id_to_auth(&span, Some(1)).unwrap(),
        Some(PdbSeqId::new(10, Some(b'A')))
    );
    assert_eq!(
        auth_seq_id_to_label(&span, Some(PdbSeqId::new(10, Some(b'a')))).unwrap(),
        Some(1)
    );
    assert_eq!(
        auth_seq_id_to_label(&span, Some(PdbSeqId::new(10, Some(b'B')))).unwrap(),
        Some(2)
    );
}

#[test]
fn sequence_interpolation_uses_later_tied_endpoint_and_clears_insertion() {
    // Source label lower_bound tie uses the later endpoint (strict <).
    let a = row(10, Some(b'A'), Some(2));
    let b = row(20, Some(b'B'), Some(8));
    let span = [&a, &b];
    assert_eq!(
        label_seq_id_to_auth(&span, Some(5)).unwrap(),
        Some(PdbSeqId::new(17, None))
    );
    assert_eq!(
        label_seq_id_to_auth(&span, Some(0)).unwrap(),
        Some(PdbSeqId::new(8, None))
    );
    assert_eq!(
        label_seq_id_to_auth(&span, Some(10)).unwrap(),
        Some(PdbSeqId::new(22, None))
    );
    assert_eq!(
        auth_seq_id_to_label(&span, Some(PdbSeqId::new(15, None))).unwrap(),
        Some(3)
    );
}

#[test]
fn sequence_missing_optional_number_propagates_and_empty_span_is_typed() {
    let a = row(10, None, Some(2));
    let span = [&a];
    assert_eq!(label_seq_id_to_auth(&span, None).unwrap(), None);
    assert_eq!(auth_seq_id_to_label(&span, None).unwrap(), None);
    assert!(matches!(
        label_seq_id_to_auth(&[], Some(2)),
        Err(BioStructureError::EmptyResidueSpan {
            operation: "label_seq_id_to_auth"
        })
    ));
    assert!(matches!(
        auth_seq_id_to_label(&[], None),
        Err(BioStructureError::EmptyResidueSpan {
            operation: "auth_seq_id_to_label"
        })
    ));
}
