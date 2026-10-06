//! Checked live-value projection of the detached batch owner.
use crate::{BatchError, BatchErrorMode, BatchRecord, BatchValidationError, MoleculeBatch};
use crate::{Molecule, SdfDataset, SdfReadParams, SdfRecord, SdfRecordStream};
use std::io::{BufRead, BufReader};

fn params(params: &SdfReadParams) -> cosmolkit_batch::BatchReadParams {
    cosmolkit_batch::BatchReadParams {
        data: crate::sdf_supplier::data_params(params),
        post: cosmolkit_io::MolPostParams {
            sanitize: params.sanitize,
            remove_hs: params.remove_hydrogens,
            expand_attachment_points: params.expand_attachment_points,
        },
    }
}
impl MoleculeBatch {
    fn from_detached(
        records: Vec<Result<cosmolkit_batch::BatchRecord, BatchError>>,
        mode: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let records = records
            .into_iter()
            .enumerate()
            .map(|(index, record)| match record {
                Ok(record) => match Molecule::from_parsed_parts_with_derived_state(
                    record.topology,
                    record.coordinates,
                    record.properties,
                    record.post_state.valence,
                    record.post_state.rings,
                ) {
                    Ok(molecule) => BatchRecord::Molecule(molecule),
                    Err(e) => BatchRecord::Error(BatchError::with_source(index, "read_sdf", e)),
                },
                Err(e) => BatchRecord::Error(e),
            })
            .collect();
        Self::from_records(records, mode)
    }
    pub fn from_sdf_records(text: &str) -> Result<Self, BatchValidationError> {
        Self::from_sdf_records_with_params(
            text,
            &SdfReadParams::default(),
            BatchErrorMode::Strict,
            None,
        )
    }
    pub fn from_sdf_records_with_params(
        text: &str,
        read: &SdfReadParams,
        mode: BatchErrorMode,
        n_jobs: Option<usize>,
    ) -> Result<Self, BatchValidationError> {
        Self::from_detached(
            cosmolkit_batch::read_sdf_text(text, params(read), n_jobs, false)?,
            mode,
        )
    }
    pub fn read_sdf(path: &str) -> Result<Self, BatchValidationError> {
        Self::read_sdf_with_params(
            path,
            &SdfReadParams::default(),
            BatchErrorMode::Strict,
            None,
            false,
        )
    }
    pub fn read_sdf_with_params(
        path: &str,
        read: &SdfReadParams,
        mode: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> Result<Self, BatchValidationError> {
        if n_jobs == Some(0) {
            return Err(BatchValidationError::parameter(
                "n_jobs",
                "n_jobs must be >= 1",
            ));
        }
        let records = if progress_bar {
            let dataset = SdfDataset::open_with_params(path, read).map_err(|e| {
                BatchValidationError::from_record_errors(vec![BatchError::with_source(
                    0, "read_sdf", e,
                )])
            })?;
            cosmolkit_batch::read_sdf_dataset(
                dataset.detached_dataset(),
                params(read),
                n_jobs,
                true,
            )?
        } else {
            let file = crate::molecular_io::open_file(path).map_err(|e| {
                BatchValidationError::from_record_errors(vec![BatchError::with_source(
                    0, "read_sdf", e,
                )])
            })?;
            cosmolkit_batch::read_sdf_reader(BufReader::new(file), params(read))
        };
        Self::from_detached(records, mode)
    }
    pub fn from_dataset_indices(
        dataset: &SdfDataset,
        indices: &[usize],
        mode: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let records = indices
            .iter()
            .map(|&index| batch_record(dataset.record(index), index))
            .collect();
        Self::from_records(records, mode)
    }
}
fn batch_record(record: Result<SdfRecord, crate::SdfError>, index: usize) -> BatchRecord {
    match record.and_then(|record| record.molecule().cloned()) {
        Ok(molecule) => BatchRecord::Molecule(molecule),
        Err(error) => BatchRecord::Error(BatchError::with_source(index, "read_sdf", error)),
    }
}
pub struct SdfBatchIterator {
    dataset: SdfDataset,
    position: usize,
    size: usize,
    mode: BatchErrorMode,
    n_jobs: Option<usize>,
    progress: Option<cosmolkit_batch::BatchProgressBar>,
}
impl SdfDataset {
    pub fn batches(
        &self,
        size: usize,
        mode: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> Result<SdfBatchIterator, BatchValidationError> {
        if size == 0 {
            return Err(BatchValidationError::parameter("size", "size must be >= 1"));
        }
        if n_jobs == Some(0) {
            return Err(BatchValidationError::parameter(
                "n_jobs",
                "n_jobs must be >= 1",
            ));
        }
        Ok(SdfBatchIterator {
            dataset: self.clone(),
            position: 0,
            size,
            mode,
            n_jobs,
            progress: progress_bar
                .then(|| cosmolkit_batch::BatchProgressBar::new(self.len(), "read_sdf_batches")),
        })
    }
}
impl Iterator for SdfBatchIterator {
    type Item = Result<MoleculeBatch, BatchValidationError>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.position >= self.dataset.len() {
            if let Some(progress) = self.progress.take() {
                progress.finish();
            }
            return None;
        }
        let start = self.position;
        let end = self.dataset.len().min(start.saturating_add(self.size));
        self.position = end;
        let records = (start..end)
            .map(|index| {
                let record = batch_record(self.dataset.record(index), index);
                if let Some(progress) = &self.progress {
                    progress.inc(1);
                }
                record
            })
            .collect();
        Some(
            MoleculeBatch::from_records(records, self.mode)
                .and_then(|b| b.with_parallel_jobs(self.n_jobs)),
        )
    }
}
impl std::iter::FusedIterator for SdfBatchIterator {}
pub struct SdfReaderBatchIterator<R> {
    reader: SdfRecordStream<R>,
    size: usize,
    mode: BatchErrorMode,
    n_jobs: Option<usize>,
}
impl<R: BufRead> SdfRecordStream<R> {
    pub fn batches(
        self,
        size: usize,
        mode: BatchErrorMode,
        n_jobs: Option<usize>,
    ) -> Result<SdfReaderBatchIterator<R>, BatchValidationError> {
        if size == 0 {
            return Err(BatchValidationError::parameter("size", "size must be >= 1"));
        }
        if n_jobs == Some(0) {
            return Err(BatchValidationError::parameter(
                "n_jobs",
                "n_jobs must be >= 1",
            ));
        }
        Ok(SdfReaderBatchIterator {
            reader: self,
            size,
            mode,
            n_jobs,
        })
    }
}
impl<R: BufRead> Iterator for SdfReaderBatchIterator<R> {
    type Item = Result<MoleculeBatch, BatchValidationError>;
    fn next(&mut self) -> Option<Self::Item> {
        let mut records = Vec::with_capacity(self.size);
        for _ in 0..self.size {
            let index = self.reader.records_consumed();
            match self.reader.next_record() {
                Ok(Some(record)) => records.push(batch_record(Ok(record), index)),
                Ok(None) => break,
                Err(e) => records.push(batch_record(Err(e), index)),
            }
        }
        if records.is_empty() {
            None
        } else {
            Some(
                MoleculeBatch::from_records(records, self.mode)
                    .and_then(|b| b.with_parallel_jobs(self.n_jobs)),
            )
        }
    }
}

impl SdfBatchIterator {
    pub fn next_batch(&mut self) -> Result<Option<MoleculeBatch>, BatchValidationError> {
        self.next().transpose()
    }
}
impl<R: BufRead> SdfReaderBatchIterator<R> {
    pub fn next_batch(&mut self) -> Result<Option<MoleculeBatch>, BatchValidationError> {
        self.next().transpose()
    }
}

pub use cosmolkit_batch::{BatchExportParams, BatchExportReport};
impl MoleculeBatch {
    fn export_records(&self) -> Vec<cosmolkit_batch::SdfExportRecord<'_>> {
        self.records
            .iter()
            .map(|r| match r {
                BatchRecord::Molecule(m) => {
                    cosmolkit_batch::SdfExportRecord::Molecule(m.mol_write_input())
                }
                BatchRecord::Error(e) => cosmolkit_batch::SdfExportRecord::Error(e),
            })
            .collect()
    }
    pub fn to_sdf(&self, path: &str) -> Result<BatchExportReport, BatchValidationError> {
        self.to_sdf_with_params(path, &BatchExportParams::default(), None)
    }
    pub fn to_sdf_with_params(
        &self,
        path: &str,
        params: &BatchExportParams,
        report_path: Option<&str>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let path = crate::molecular_io::expanded_path(path).map_err(|e| {
            BatchValidationError::from_record_errors(vec![BatchError::with_source(
                0,
                "batch.write_sdf",
                e,
            )])
        })?;
        let params = BatchExportParams {
            n_jobs: params.n_jobs.or(self.n_jobs),
            progress_bar: params.progress_bar.or(self.progress_bar),
            ..*params
        };
        let report = cosmolkit_batch::export_sdf(&self.export_records(), &path, &params)?;
        Self::write_report(&report, report_path)?;
        Ok(report)
    }
    pub fn to_sdf_files(&self, directory: &str) -> Result<BatchExportReport, BatchValidationError> {
        self.to_sdf_files_with_params(directory, &BatchExportParams::default(), None, None)
    }
    pub fn to_sdf_files_with_params(
        &self,
        directory: &str,
        params: &BatchExportParams,
        filenames: Option<&[Option<String>]>,
        report_path: Option<&str>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let directory = crate::molecular_io::expanded_path(directory).map_err(|e| {
            BatchValidationError::from_record_errors(vec![BatchError::with_source(
                0,
                "batch.write_sdf_files",
                e,
            )])
        })?;
        let params = BatchExportParams {
            n_jobs: params.n_jobs.or(self.n_jobs),
            progress_bar: params.progress_bar.or(self.progress_bar),
            ..*params
        };
        let report = cosmolkit_batch::export_sdf_files(
            &self.export_records(),
            &directory,
            &params,
            filenames,
        )?;
        Self::write_report(&report, report_path)?;
        Ok(report)
    }
    fn write_report(
        report: &BatchExportReport,
        path: Option<&str>,
    ) -> Result<(), BatchValidationError> {
        if let Some(path) = path {
            let path = crate::molecular_io::expanded_path(path).map_err(|e| {
                BatchValidationError::from_record_errors(vec![BatchError::with_source(
                    0,
                    "write error report",
                    e,
                )])
            })?;
            cosmolkit_batch::write_export_report(&path, report)?;
        }
        Ok(())
    }
}
