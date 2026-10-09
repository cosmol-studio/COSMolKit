#![cfg(all(feature = "cap-io", feature = "cap-batch"))]
// Author proposal: fixed p4 V8 supplier and batch source, not an accepted condition.
use cosmolkit::{
    BatchErrorMode, BatchRecord, MoleculeBatch, SdfDataset, SdfReadParams, SdfReader,
    SdfRecordStream,
};
use std::io::Cursor;
use std::sync::atomic::{AtomicU64, Ordering};

fn record(index: usize) -> String {
    format!(
        "row{index}\n  COSMolKit         2D\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n{:>10.4}{:>10.4}{:>10.4} C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n$$$$\n",
        index as f64, 0.0, 0.0
    )
}
fn input() -> String {
    [
        record(0),
        "bad1\n$$$$\n".into(),
        record(2),
        "bad3\n$$$$\n".into(),
        record(4),
    ]
    .concat()
}
fn params() -> SdfReadParams {
    SdfReadParams {
        sanitize: false,
        remove_hs: false,
        ..SdfReadParams::default()
    }
}
struct Fixture(std::path::PathBuf);
impl Fixture {
    fn new(text: &str) -> Self {
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let path = std::env::temp_dir().join(format!(
            "cosmolkit-p5-sdf4-{}-{}.sdf",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&path)
            .unwrap();
        std::io::Write::write_all(&mut file, text.as_bytes()).unwrap();
        Self(path)
    }
    fn path(&self) -> &str {
        self.0.to_str().unwrap()
    }
}
impl Drop for Fixture {
    fn drop(&mut self) {
        std::fs::remove_file(&self.0).unwrap();
    }
}
fn positions(batch: &MoleculeBatch) -> Vec<Option<f64>> {
    batch
        .records()
        .iter()
        .map(|record| match record {
            BatchRecord::Molecule(molecule) => Some(molecule.coordinates_2d().unwrap()[0][0]),
            BatchRecord::Error(_) => None,
        })
        .collect()
}

#[test]
fn dataset_keep_preserves_chunks_original_error_indices_jobs_and_duplicate_selection() {
    let file = Fixture::new(&input());
    let dataset = SdfDataset::open_with_params(file.path(), &params()).unwrap();
    assert_eq!(dataset.len(), 5);
    let mut batches = dataset
        .batches(3, BatchErrorMode::KeepErrors, Some(2), false)
        .unwrap();
    let first = batches.next_batch().unwrap().unwrap();
    assert_eq!(positions(&first), vec![Some(0.0), None, Some(2.0)]);
    assert_eq!(first.valid_mask(), vec![true, false, true]);
    assert_eq!(
        first.errors().iter().map(|e| e.index).collect::<Vec<_>>(),
        vec![1]
    );
    assert_eq!(first.errors()[0].operation, "read_sdf");
    assert_eq!(first.parallel_jobs(), Some(2));
    let second = batches.next_batch().unwrap().unwrap();
    assert_eq!(positions(&second), vec![None, Some(4.0)]);
    assert_eq!(
        second.errors().iter().map(|e| e.index).collect::<Vec<_>>(),
        vec![3]
    );
    assert!(batches.next_batch().unwrap().is_none());
    assert!(batches.next_batch().unwrap().is_none());
    let selected =
        MoleculeBatch::from_dataset_indices(&dataset, &[4, 0, 4], BatchErrorMode::Strict).unwrap();
    assert_eq!(positions(&selected), vec![Some(4.0), Some(0.0), Some(4.0)]);
}

#[test]
fn dataset_strict_consumes_whole_failed_chunk_and_continues_at_next_record() {
    let file = Fixture::new(&input());
    let dataset = SdfDataset::open_with_params(file.path(), &params()).unwrap();
    let mut batches = dataset
        .batches(4, BatchErrorMode::Strict, None, true)
        .unwrap();
    let error = batches.next_batch().unwrap_err();
    assert_eq!(error.errors, 2);
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1, 3]
    );
    assert!(
        error
            .record_errors
            .iter()
            .all(|e| e.operation == "read_sdf")
    );
    assert_eq!(
        positions(&batches.next_batch().unwrap().unwrap()),
        vec![Some(4.0)]
    );
    assert!(batches.next_batch().unwrap().is_none());
}

#[test]
fn forward_stream_strict_recovers_and_file_reader_reopens_with_same_keep_configuration() {
    // MolFromMolDataStream reads name/info/comments/counts before recovery.
    // Retain the original short headers: their $$$$ is consumed as info and
    // the next good record's prefix as comments/counts. Recovery then consumes
    // that good record through its delimiter, exactly as the pinned supplier.
    let mut stream = SdfRecordStream::with_params(Cursor::new(input().into_bytes()), params())
        .batches(4, BatchErrorMode::Strict, Some(3))
        .unwrap();
    let error = stream.next_batch().unwrap_err();
    assert_eq!(error.errors, 2);
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1, 2]
    );
    assert!(
        error
            .record_errors
            .iter()
            .all(|e| e.operation == "read_sdf" && std::error::Error::source(e).is_some())
    );
    assert!(stream.next_batch().unwrap().is_none());
    let original = Fixture::new(&input());
    let reader = SdfReader::open_with_params(original.path(), &params()).unwrap();
    for _ in 0..2 {
        let mut batches = reader
            .batches(2, BatchErrorMode::KeepErrors, Some(2))
            .unwrap();
        let first = batches.next_batch().unwrap().unwrap();
        assert_eq!(positions(&first), vec![Some(0.0), None]);
        assert_eq!(first.errors()[0].index, 1);
        assert_eq!(first.parallel_jobs(), Some(2));
        let last = batches.next_batch().unwrap().unwrap();
        assert_eq!(positions(&last), vec![None]);
        assert_eq!(last.errors()[0].index, 2);
        assert_eq!(last.parallel_jobs(), Some(2));
        assert!(batches.next_batch().unwrap().is_none());
    }
    // Full malformed-count headers fail before the unread delimiter, so all
    // original forward recovery, chunk, job and reopen assertions apply.
    let framed = [
        record(0),
        "bad1\n  COSMolKit\n\ninvalid\n$$$$\n".into(),
        record(2),
        "bad3\n  COSMolKit\n\ninvalid\n$$$$\n".into(),
        record(4),
    ]
    .concat();
    let mut stream = SdfRecordStream::with_params(Cursor::new(framed.as_bytes()), params())
        .batches(4, BatchErrorMode::Strict, Some(3))
        .unwrap();
    let error = stream.next_batch().unwrap_err();
    assert_eq!(
        error
            .record_errors
            .iter()
            .map(|e| e.index)
            .collect::<Vec<_>>(),
        vec![1, 3]
    );
    let last = stream.next_batch().unwrap().unwrap();
    assert_eq!(positions(&last), vec![Some(4.0)]);
    assert_eq!(last.parallel_jobs(), Some(3));
    assert!(stream.next_batch().unwrap().is_none());
    let file = Fixture::new(&framed);
    let reader = SdfReader::open_with_params(file.path(), &params()).unwrap();
    for _ in 0..2 {
        let mut batches = reader
            .batches(2, BatchErrorMode::KeepErrors, Some(2))
            .unwrap();
        assert_eq!(
            positions(&batches.next_batch().unwrap().unwrap()),
            vec![Some(0.0), None]
        );
        assert_eq!(batches.next_batch().unwrap().unwrap().errors()[0].index, 3);
        assert_eq!(
            positions(&batches.next_batch().unwrap().unwrap()),
            vec![Some(4.0)]
        );
        assert!(batches.next_batch().unwrap().is_none());
    }
}

#[test]
fn size_precedes_jobs_and_reader_io_failure_retains_typed_cause() {
    let file = Fixture::new(&record(0));
    let dataset = SdfDataset::open_with_params(file.path(), &params()).unwrap();
    let error = dataset
        .batches(0, BatchErrorMode::Strict, Some(0), false)
        .err()
        .unwrap();
    assert_eq!(error.record_errors[0].operation, "size");
    let error = dataset
        .batches(1, BatchErrorMode::Strict, Some(0), false)
        .err()
        .unwrap();
    assert_eq!(error.record_errors[0].operation, "n_jobs");
    let reader =
        SdfReader::open_with_params(&(file.path().to_owned() + ".missing"), &params()).unwrap();
    assert_eq!(
        reader
            .batches(0, BatchErrorMode::Strict, Some(0))
            .err()
            .unwrap()
            .record_errors[0]
            .operation,
        "size"
    );
    let error = reader
        .batches(1, BatchErrorMode::Strict, None)
        .err()
        .unwrap();
    assert_eq!(error.record_errors[0].operation, "SdfReader.open");
    assert!(std::error::Error::source(&error.record_errors[0]).is_some());
}
