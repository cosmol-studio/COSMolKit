//! Delivery proposals against d892ec3 properties/batch.rs indexed scheduling.
use cosmolkit_batch::{BatchErrorMode, BatchRecordError, run_indexed, validate_record_errors};
use std::sync::atomic::{AtomicUsize, Ordering};

#[test]
fn ordered_completion_preserves_all_errors_and_configured_pool_width() {
    for jobs in [None, Some(1), Some(2), Some(4)] {
        let calls = AtomicUsize::new(0);
        let rows = run_indexed(127, jobs, Some(false), "proposal", |index| {
            calls.fetch_add(1, Ordering::SeqCst);
            assert_eq!(rayon::current_num_threads(), jobs.unwrap_or(1));
            if index % 7 == 0 {
                Err(BatchRecordError::with_source(
                    index,
                    "proposal.operation",
                    std::io::Error::new(std::io::ErrorKind::InvalidData, "original failure"),
                ))
            } else {
                Ok(index * 3)
            }
        })
        .unwrap();
        assert_eq!(calls.load(Ordering::SeqCst), 127);
        assert_eq!(rows.len(), 127);
        let mut errors = Vec::new();
        for (index, row) in rows.into_iter().enumerate() {
            if index % 7 == 0 {
                let error = row.unwrap_err();
                assert_eq!(error.index, index);
                assert!(std::error::Error::source(&error).is_some());
                errors.push(error);
            } else {
                assert_eq!(row.unwrap(), index * 3);
            }
        }
        validate_record_errors(errors.clone(), BatchErrorMode::KeepErrors).unwrap();
        let strict = validate_record_errors(errors, BatchErrorMode::Strict).unwrap_err();
        assert_eq!(strict.errors, 19);
        assert_eq!(
            strict
                .record_errors
                .iter()
                .map(|e| e.index)
                .collect::<Vec<_>>(),
            (0..127).step_by(7).collect::<Vec<_>>()
        );
    }
}
#[test]
fn invalid_configuration_fails_before_any_work_and_empty_work_still_validates() {
    let calls = AtomicUsize::new(0);
    assert!(
        run_indexed(4, Some(0), None, "proposal", |_| {
            calls.fetch_add(1, Ordering::SeqCst);
        })
        .is_err()
    );
    assert_eq!(calls.load(Ordering::SeqCst), 0);
    assert!(run_indexed(0, Some(0), None, "proposal", |_| ()).is_err());
    assert!(
        run_indexed(0, Some(1), Some(true), "proposal", |_| ())
            .unwrap()
            .is_empty()
    );
}
