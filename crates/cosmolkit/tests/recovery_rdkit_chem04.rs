use cosmolkit::{CipLabelerError, OperationError};
use std::error::Error;
#[test]
fn new_empty_pairlist_leaf_survives_existing_registered_error_envelope() {
    let e = OperationError::CipLabeler(CipLabelerError::EmptyPairListReference);
    assert_eq!(
        e.source().unwrap().downcast_ref::<CipLabelerError>(),
        Some(&CipLabelerError::EmptyPairListReference)
    );
    assert!(
        e.to_string()
            .contains("Cannot get a reference from an empty PairList")
    );
    let e = OperationError::AtomCode(cosmolkit_fingerprints::AtomCodeError::Cip(
        CipLabelerError::EmptyPairListReference,
    ));
    assert_eq!(
        e.source()
            .unwrap()
            .source()
            .unwrap()
            .downcast_ref::<CipLabelerError>(),
        Some(&CipLabelerError::EmptyPairListReference)
    );
}
#[test]
fn new_empty_pairlist_leaf_survives_existing_fingerprint_and_batch_envelopes() {
    use cosmolkit_batch::BatchRecordError;
    use cosmolkit_fingerprints::{AtomCodeError, AtomPairError, MorganError};
    let errors = [
        BatchRecordError::with_source(
            7,
            "atom_pair_fingerprint",
            cosmolkit::AtomPairReadError::Generator(AtomPairError::AtomCode(AtomCodeError::Cip(
                CipLabelerError::EmptyPairListReference,
            ))),
        ),
        BatchRecordError::with_source(
            7,
            "atom_pair_fingerprint",
            MorganError::AtomPair(Box::new(AtomPairError::AtomCode(AtomCodeError::Cip(
                CipLabelerError::EmptyPairListReference,
            )))),
        ),
    ];
    for e in errors {
        let mut cause: &(dyn Error + 'static) = &e;
        let mut found = false;
        let mut depth = 0;
        loop {
            if cause.downcast_ref::<CipLabelerError>()
                == Some(&CipLabelerError::EmptyPairListReference)
            {
                found = true;
            }
            depth += 1;
            match cause.source() {
                Some(next) => cause = next,
                None => break,
            }
        }
        assert!(found);
        assert_eq!(depth, 5);
        assert_eq!(e.index, 7);
        assert!(
            e.to_string()
                .contains("Cannot get a reference from an empty PairList")
        );
    }
}
