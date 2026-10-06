//! Source-backed SDF batch scheduling over explicit detached values only.
use crate::BatchRecord;
use cosmolkit_io::{
    MolBlockRecord, MolPostParams, SdfDataReadParams, SdfGraphDataset, SdfGraphReader,
    SdfGraphRecord,
};
use rayon::prelude::*;
use std::io::BufRead;
use std::sync::Arc;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BatchErrorMode {
    #[default]
    Strict,
    KeepErrors,
}

#[derive(Debug, Clone)]
pub struct BatchRecordError {
    pub index: usize,
    pub operation: &'static str,
    pub message: String,
    source: Option<Arc<dyn std::error::Error + Send + Sync>>,
}
impl BatchRecordError {
    pub fn new(index: usize, operation: &'static str, message: impl Into<String>) -> Self {
        Self {
            index,
            operation,
            message: message.into(),
            source: None,
        }
    }
    pub fn with_source<E: std::error::Error + Send + Sync + 'static>(
        index: usize,
        operation: &'static str,
        source: E,
    ) -> Self {
        Self {
            index,
            operation,
            message: source.to_string(),
            source: Some(Arc::new(source)),
        }
    }
}
impl std::fmt::Display for BatchRecordError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} at {}: {}", self.operation, self.index, self.message)
    }
}
impl std::error::Error for BatchRecordError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.source
            .as_deref()
            .map(|e| e as &(dyn std::error::Error + 'static))
    }
}
#[derive(Debug, Clone, thiserror::Error)]
#[error("batch validation failed with {errors} errors{message_suffix}")]
pub struct BatchValidationError {
    pub errors: usize,
    pub record_errors: Vec<BatchRecordError>,
    message_suffix: String,
}
impl BatchValidationError {
    pub fn from_record_errors(record_errors: Vec<BatchRecordError>) -> Self {
        let message_suffix = record_errors
            .first()
            .map(|e| format!("; first error at {}: {}", e.operation, e.message))
            .unwrap_or_default();
        Self {
            errors: record_errors.len(),
            record_errors,
            message_suffix,
        }
    }
    pub fn parameter(name: &'static str, message: &'static str) -> Self {
        Self::from_record_errors(vec![BatchRecordError::new(0, name, message)])
    }
}
/// The common source error mode is evaluated only after the entire selected
/// chunk is consumed, retaining every failure in its original input order.
pub fn validate_record_errors(
    errors: Vec<BatchRecordError>,
    mode: BatchErrorMode,
) -> Result<(), BatchValidationError> {
    if mode == BatchErrorMode::Strict && !errors.is_empty() {
        Err(BatchValidationError::from_record_errors(errors))
    } else {
        Ok(())
    }
}
#[derive(Debug, Clone, Copy, Default)]
pub struct BatchReadParams {
    pub data: SdfDataReadParams,
    pub post: MolPostParams,
}
pub struct BatchProgressBar {
    inner: indicatif::ProgressBar,
}
impl BatchProgressBar {
    pub fn new(total: usize, message: &'static str) -> Self {
        let inner = indicatif::ProgressBar::with_draw_target(
            Some(total as u64),
            indicatif::ProgressDrawTarget::stderr_with_hz(20),
        );
        let style = indicatif::ProgressStyle::with_template(
            "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
        )
        .unwrap_or_else(|_| indicatif::ProgressStyle::default_bar());
        inner.set_style(style);
        inner.set_message(message);
        Self { inner }
    }
    pub fn inc(&self, n: u64) {
        self.inner.inc(n);
    }
    pub fn finish(&self) {
        self.inner.finish();
    }
}
pub fn finalize_sdf_record(
    record: SdfGraphRecord,
    post: MolPostParams,
    index: usize,
    operation: &'static str,
) -> Result<BatchRecord, BatchRecordError> {
    let mut record = record
        .finish_mol_post(post)
        .map_err(|e| BatchRecordError::with_source(index, operation, e))?;
    let post_state = record.take_post_state();
    match record.mol_block {
        MolBlockRecord::Concrete {
            topology,
            coordinates,
            properties,
        } => Ok(BatchRecord {
            topology,
            coordinates,
            properties,
            post_state,
        }),
        MolBlockRecord::Query(_) => Err(BatchRecordError::with_source(
            index,
            operation,
            cosmolkit_io::SdfReadError::QueryRecord,
        )),
    }
}
pub(super) use crate::scheduler::parallel;

pub fn read_sdf_text(
    text: &str,
    params: BatchReadParams,
    n_jobs: Option<usize>,
    progress_bar: bool,
) -> Result<Vec<Result<BatchRecord, BatchRecordError>>, BatchValidationError> {
    // Pinned legacy properties/batch.rs::read_sdf_records_from_str_with_params_and_options:
    // split all record strings, indexed par_iter+enumerate+collect, then validate.
    // Keep the same Rayon indexed ordering and one-record detached finalization;
    // no live Molecule, runtime cache or commit authority enters this owner.
    let records = cosmolkit_io::split_sdf_record_strings(text);
    let progress =
        progress_bar.then(|| BatchProgressBar::new(records.len(), "Reading SDF records"));
    let result = parallel(n_jobs, || {
        records
            .par_iter()
            .enumerate()
            .map(|(index, text)| {
                let operation = "batch.read_sdf_records_from_str";
                let result =
                    cosmolkit_io::read_sdf_graph_record_detached_with_params(text, params.data)
                        .map_err(|e| BatchRecordError::with_source(index, operation, e))
                        .and_then(|record| {
                            finalize_sdf_record(record, params.post, index, operation)
                        });
                if let Some(progress) = &progress {
                    progress.inc(1);
                }
                result
            })
            .collect()
    });
    if let Some(progress) = progress {
        progress.finish();
    }
    result
}
pub fn read_sdf_dataset(
    dataset: &SdfGraphDataset,
    params: BatchReadParams,
    n_jobs: Option<usize>,
    progress_bar: bool,
) -> Result<Vec<Result<BatchRecord, BatchRecordError>>, BatchValidationError> {
    let progress =
        progress_bar.then(|| BatchProgressBar::new(dataset.len(), "Reading SDF dataset"));
    let result = parallel(n_jobs, || {
        (0..dataset.len())
            .into_par_iter()
            .map(|index| {
                let operation = "batch.read_sdf_dataset";
                let result = dataset
                    .record_with_params(index, params.data)
                    .map_err(|e| BatchRecordError::with_source(index, operation, e))
                    .and_then(|record| finalize_sdf_record(record, params.post, index, operation));
                if let Some(progress) = &progress {
                    progress.inc(1);
                }
                result
            })
            .collect()
    });
    if let Some(progress) = progress {
        progress.finish();
    }
    result
}
pub fn read_sdf_reader<R: BufRead>(
    reader: R,
    params: BatchReadParams,
) -> Vec<Result<BatchRecord, BatchRecordError>> {
    // Source forward path is sequential even when the caller supplies n_jobs.
    let mut reader = SdfGraphReader::with_params(reader, params.data);
    let mut records = Vec::new();
    loop {
        let index = reader.records_consumed();
        match reader.next_record() {
            Ok(Some(record)) => {
                records.push(finalize_sdf_record(record, params.post, index, "read_sdf"))
            }
            Ok(None) => break,
            Err(error) => {
                records.push(Err(BatchRecordError::with_source(index, "read_sdf", error)))
            }
        }
    }
    records
}
