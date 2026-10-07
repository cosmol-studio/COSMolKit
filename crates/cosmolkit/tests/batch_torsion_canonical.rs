#![cfg(all(
    feature = "cap-batch",
    feature = "cap-fingerprints",
    feature = "cap-smiles"
))]
use cosmolkit::*;
use std::error::Error as _;
use std::sync::{
    Arc,
    atomic::{AtomicUsize, Ordering},
};

fn original_batch() -> MoleculeBatch {
    MoleculeBatch::from_records(
        vec![
            BatchRecord::Molecule(Molecule::from_smiles("CCCCO").unwrap()),
            BatchRecord::Error(BatchError::new(
                1,
                "original-input",
                "retained invalid input",
            )),
            BatchRecord::Molecule(Molecule::from_smiles("CCCCC").unwrap()),
        ],
        BatchErrorMode::KeepErrors,
    )
    .unwrap()
    .with_parallel_jobs(Some(2))
    .unwrap()
    .with_progress_bar(Some(false))
}

#[test]
fn torsion_defaults_overrides_order_invalid_rows_and_progress_preserve_inputs() {
    let batch = original_batch();
    let before = batch.to_list();
    let options = TopologicalTorsionFingerprintParams::default();
    let expected = batch
        .records()
        .iter()
        .map(|row| match row {
            BatchRecord::Molecule(m) => Some(
                m.topological_torsion_fingerprint_with_params(&options, None)
                    .unwrap(),
            ),
            BatchRecord::Error(_) => None,
        })
        .collect::<Vec<_>>();
    assert_eq!(
        batch.fingerprint_topological_torsion_list().unwrap(),
        expected
    );
    for jobs in [1, 2] {
        let ticks = Arc::new(AtomicUsize::new(0));
        let worker = Arc::clone(&ticks);
        let params = BatchQueryParams {
            n_jobs: Some(jobs),
            progress_bar: Some(false),
            progress_callback: Some(Arc::new(move || {
                worker.fetch_add(1, Ordering::SeqCst);
            })),
        };
        assert_eq!(
            batch
                .fingerprint_topological_torsion_list_with_params(&options, &params)
                .unwrap(),
            expected
        );
        assert_eq!(ticks.load(Ordering::SeqCst), 3);
        assert_eq!(batch.to_list(), before);
        assert_eq!(batch.parallel_jobs(), Some(2));
        assert_eq!(batch.progress_bar(), Some(false));
    }
    assert_eq!(batch.errors()[0].index, 1);
    assert_eq!(batch.errors()[0].message, "retained invalid input");
}

#[test]
fn torsion_collects_every_typed_failure_in_original_order_without_skipping_work() {
    let batch = original_batch();
    let before = batch.to_list();
    let mut options = TopologicalTorsionFingerprintParams::default();
    options.from_atoms = Some(vec![5]);
    let ticks = Arc::new(AtomicUsize::new(0));
    let worker = Arc::clone(&ticks);
    let params = BatchQueryParams {
        n_jobs: Some(2),
        progress_bar: Some(false),
        progress_callback: Some(Arc::new(move || {
            worker.fetch_add(1, Ordering::SeqCst);
        })),
    };
    let error = batch
        .fingerprint_topological_torsion_list_with_params(&options, &params)
        .unwrap_err();
    assert_eq!(error.errors, 2);
    assert_eq!(ticks.load(Ordering::SeqCst), 3);
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        [0, 2]
    );
    for record in &error.record_errors {
        assert_eq!(record.operation, "batch.topological_torsion_fingerprint");
        let original = record
            .source()
            .unwrap()
            .downcast_ref::<TopologicalTorsionReadError>()
            .unwrap();
        assert!(matches!(
            original,
            TopologicalTorsionReadError::Generator(_)
        ));
        let reason = original
            .source()
            .and_then(|generator| generator.source())
            .and_then(|source| source.downcast_ref::<FingerprintError>())
            .expect("the public error chain retains the fingerprint selection reason");
        assert!(matches!(
            reason,
            FingerprintError::InvalidArguments {
                reason: "atom selection contains an atom index outside the molecule"
            }
        ));
    }
    assert_eq!(batch.to_list(), before);
}

#[test]
fn zero_threads_reject_before_work_and_empty_batches_keep_source_shape() {
    let batch = original_batch();
    let ticks = Arc::new(AtomicUsize::new(0));
    let worker = Arc::clone(&ticks);
    let params = BatchQueryParams {
        n_jobs: Some(0),
        progress_bar: Some(false),
        progress_callback: Some(Arc::new(move || {
            worker.fetch_add(1, Ordering::SeqCst);
        })),
    };
    let error = batch
        .fingerprint_topological_torsion_list_with_params(&Default::default(), &params)
        .unwrap_err();
    assert_eq!(error.record_errors[0].operation, "n_jobs");
    assert_eq!(ticks.load(Ordering::SeqCst), 0);
    let empty = MoleculeBatch::from_records(vec![], BatchErrorMode::Strict).unwrap();
    assert!(
        empty
            .fingerprint_topological_torsion_list()
            .unwrap()
            .is_empty()
    );
}

#[test]
fn oversized_fixed_bit_path_retains_source_empty_bits_and_original_invalid_row() {
    // The fixed-bit source route hashes environments into the configured size;
    // it does not call getResultSize, which belongs to sparse output selection.
    let batch = original_batch();
    let before = batch.to_list();
    let mut options = TopologicalTorsionFingerprintParams::default();
    options.generator.torsion_atom_count = u32::MAX;
    let ticks = Arc::new(AtomicUsize::new(0));
    let worker = Arc::clone(&ticks);
    let params = BatchQueryParams {
        n_jobs: Some(2),
        progress_bar: Some(false),
        progress_callback: Some(Arc::new(move || {
            worker.fetch_add(1, Ordering::SeqCst);
        })),
    };
    let values = batch
        .fingerprint_topological_torsion_list_with_params(&options, &params)
        .unwrap();
    let empty = Fingerprint::from_on_bits(2048, vec![]).unwrap();
    assert_eq!(values, [Some(empty.clone()), None, Some(empty)]);
    assert_eq!(ticks.load(Ordering::SeqCst), 3);
    assert_eq!(batch.to_list(), before);
}
