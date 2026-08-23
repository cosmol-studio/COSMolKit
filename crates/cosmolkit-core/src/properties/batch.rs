#[cfg(test)]
use crate::draw::PreparedDrawMolecule;
use crate::io::molblock::{self, SdfFormat};
use crate::io::sdf::{
    SdfCoordinateMode, SdfDataset, SdfReadParams, SdfReader, read_sdf_from_str_with_params,
};
use crate::{Molecule, SmilesWriteParams};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::{BufRead, Write};
use std::path::{Component, Path, PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchErrorMode {
    Strict,
    KeepErrors,
}

impl BatchErrorMode {
    const fn raise_on_errors(self) -> bool {
        matches!(self, Self::Strict)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchRecordError {
    pub index: usize,
    pub operation: &'static str,
    pub message: String,
}

impl BatchRecordError {
    #[must_use]
    pub fn new(index: usize, operation: &'static str, message: impl Into<String>) -> Self {
        Self {
            index,
            operation,
            message: message.into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BatchRecord {
    Molecule(Molecule),
    Error(BatchRecordError),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MoleculeBatch {
    records: Vec<BatchRecord>,
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BatchExportReport {
    pub written: usize,
    pub skipped: usize,
    pub failed: usize,
    pub errors: Vec<BatchRecordError>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("batch validation failed with {errors} errors{message_suffix}")]
pub struct BatchValidationError {
    pub errors: usize,
    pub reason: Option<crate::UnsupportedFeatureError>,
    pub record_errors: Vec<BatchRecordError>,
    #[doc(hidden)]
    pub message_suffix: String,
}

impl BatchValidationError {
    pub fn from_record_errors(record_errors: Vec<BatchRecordError>) -> Self {
        let message_suffix = record_errors
            .first()
            .map(|error| format!("; first error at {}: {}", error.operation, error.message))
            .unwrap_or_default();
        Self {
            errors: record_errors.len(),
            reason: None,
            record_errors,
            message_suffix,
        }
    }

    pub fn unsupported(feature: &'static str, reason: &'static str) -> Self {
        Self {
            errors: 1,
            reason: Some(crate::UnsupportedFeatureError { feature, reason }),
            record_errors: Vec::new(),
            message_suffix: format!("; {feature}: {reason}"),
        }
    }

    pub fn simple(errors: usize, operation: &'static str, message: impl Into<String>) -> Self {
        let message = message.into();
        let record_errors = if errors == 0 {
            Vec::new()
        } else {
            vec![BatchRecordError::new(0, operation, message.clone())]
        };
        Self {
            errors,
            reason: None,
            record_errors,
            message_suffix: format!("; first error at {operation}: {message}"),
        }
    }
}

fn record_errors_from_records(records: &[BatchRecord]) -> Vec<BatchRecordError> {
    records
        .iter()
        .filter_map(|record| match record {
            BatchRecord::Molecule(_) => None,
            BatchRecord::Error(error) => Some(error.clone()),
        })
        .collect()
}

pub type BatchProgress<'a> = Option<&'a (dyn Fn() + Sync)>;

pub struct BatchProgressBar {
    inner: ProgressBar,
}

impl BatchProgressBar {
    #[must_use]
    pub fn new(total: usize, message: impl Into<String>) -> Self {
        let progress_bar = ProgressBar::with_draw_target(
            Some(total as u64),
            ProgressDrawTarget::stderr_with_hz(20),
        );
        let style = ProgressStyle::with_template(
            "{spinner:.green} {msg} [{elapsed_precise}] [{wide_bar:.cyan/blue}] {pos}/{len}",
        )
        .unwrap_or_else(|_| ProgressStyle::default_bar());
        progress_bar.set_style(style);
        progress_bar.set_message(message.into());
        Self {
            inner: progress_bar,
        }
    }

    #[must_use]
    pub fn callback(&self) -> Box<dyn Fn() + Sync + '_> {
        Box::new(|| self.inner.inc(1))
    }

    pub fn inc(&self, delta: u64) {
        self.inner.inc(delta);
    }

    pub fn finish(&self) {
        self.inner.finish();
    }
}

pub fn batch_progress_bar(total: usize, message: impl Into<String>) -> BatchProgressBar {
    BatchProgressBar::new(total, message)
}

fn tick_progress(progress: BatchProgress<'_>) {
    if let Some(progress) = progress {
        progress();
    }
}

fn with_progress_bar_for_option<R>(
    progress_bar: Option<bool>,
    total: usize,
    message: &'static str,
    f: impl FnOnce(BatchProgress<'_>) -> R,
) -> R {
    if progress_bar.unwrap_or(false) {
        let progress_bar = batch_progress_bar(total, message);
        let callback = progress_bar.callback();
        let result = f(Some(&*callback));
        progress_bar.finish();
        result
    } else {
        f(None)
    }
}

fn run_with_parallel_jobs_option<R: Send>(
    n_jobs: Option<usize>,
    f: impl FnOnce() -> R + Send,
) -> R {
    match n_jobs.map(|value| value.max(1)) {
        Some(n_jobs) => rayon::ThreadPoolBuilder::new()
            .num_threads(n_jobs)
            .build()
            .expect("batch rayon thread pool build must succeed")
            .install(f),
        None => f(),
    }
}

fn atom_pair_generator_for_batch(
    params: &crate::AtomPairFingerprintParams,
    operation: &'static str,
) -> Result<crate::AtomPairFingerprintGenerator, BatchValidationError> {
    crate::AtomPairFingerprintGenerator::new(params).map_err(|error| {
        BatchValidationError::simple(1, operation, format!("invalid configuration: {error}"))
    })
}

impl MoleculeBatch {
    #[must_use]
    pub fn new(records: Vec<BatchRecord>) -> Self {
        Self {
            records,
            n_jobs: None,
            progress_bar: None,
        }
    }

    #[must_use]
    pub fn with_parallel_jobs(mut self, n_jobs: Option<usize>) -> Self {
        self.n_jobs = n_jobs.map(|value| value.max(1));
        self
    }

    #[must_use]
    pub fn parallel_jobs(&self) -> Option<usize> {
        self.n_jobs
    }

    #[must_use]
    pub fn with_progress_bar(mut self, progress_bar: Option<bool>) -> Self {
        self.progress_bar = progress_bar;
        self
    }

    #[must_use]
    pub fn progress_bar(&self) -> Option<bool> {
        self.progress_bar
    }

    fn effective_progress_bar(&self, progress_bar: Option<bool>) -> bool {
        progress_bar.or(self.progress_bar).unwrap_or(false)
    }

    fn effective_parallel_jobs(&self, n_jobs: Option<usize>) -> Option<usize> {
        n_jobs.or(self.n_jobs).map(|value| value.max(1))
    }

    fn with_progress_bar_for<R>(
        &self,
        progress_bar: Option<bool>,
        total: usize,
        message: &'static str,
        f: impl FnOnce(BatchProgress<'_>) -> R,
    ) -> R {
        if self.effective_progress_bar(progress_bar) {
            let progress_bar = batch_progress_bar(total, message);
            let callback = progress_bar.callback();
            let result = f(Some(&*callback));
            progress_bar.finish();
            result
        } else {
            f(None)
        }
    }

    fn run_with_parallel_jobs<R: Send>(
        &self,
        n_jobs: Option<usize>,
        f: impl FnOnce() -> R + Send,
    ) -> R {
        match self.effective_parallel_jobs(n_jobs) {
            Some(n_jobs) => rayon::ThreadPoolBuilder::new()
                .num_threads(n_jobs)
                .build()
                .expect("batch rayon thread pool build must succeed")
                .install(f),
            None => f(),
        }
    }

    #[must_use]
    pub fn get(&self, index: usize) -> Option<&BatchRecord> {
        self.records.get(index)
    }

    pub fn iter(&self) -> impl Iterator<Item = &BatchRecord> {
        self.records.iter()
    }

    #[must_use]
    pub fn error_summary(&self) -> String {
        let errors = self.errors();
        if errors.is_empty() {
            return String::new();
        }
        let total = self.records.len();
        let valid = total - errors.len();
        let mut lines = Vec::with_capacity(errors.len() + 1);
        lines.push(format!(
            "Batch error summary: {} errors out of {} records ({} valid)",
            errors.len(),
            total,
            valid
        ));
        for error in &errors {
            lines.push(format!(
                "  [{}] {}: {}",
                error.index, error.operation, error.message
            ));
        }
        lines.join("\n")
    }

    #[must_use]
    pub fn error_report(&self) -> Vec<(usize, &'static str, String)> {
        self.errors()
            .into_iter()
            .map(|error| (error.index, error.operation, error.message))
            .collect()
    }

    #[must_use]
    pub fn from_smiles_list(smiles: &[String]) -> Self {
        Self {
            records: smiles
                .iter()
                .enumerate()
                .map(|(index, smiles)| match Molecule::from_smiles(smiles) {
                    Ok(molecule) => BatchRecord::Molecule(molecule),
                    Err(error) => BatchRecord::Error(BatchRecordError::new(
                        index,
                        "batch.from_smiles_list",
                        error.to_string(),
                    )),
                })
                .collect(),
            n_jobs: None,
            progress_bar: None,
        }
    }

    pub fn from_smiles_list_with_sanitize(
        smiles: &[String],
        sanitize: bool,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        Self::from_smiles_list_with_sanitize_and_options(smiles, sanitize, errors, None, None)
    }

    pub fn from_smiles_list_with_sanitize_and_options(
        smiles: &[String],
        sanitize: bool,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        with_progress_bar_for_option(progress_bar, smiles.len(), "Parsing SMILES", |progress| {
            let records = run_with_parallel_jobs_option(n_jobs, || {
                smiles
                    .par_iter()
                    .enumerate()
                    .map(|(index, smiles)| {
                        let out = match Molecule::from_smiles_with_sanitize(smiles, sanitize) {
                            Ok(molecule) => BatchRecord::Molecule(molecule),
                            Err(error) => BatchRecord::Error(BatchRecordError::new(
                                index,
                                "batch.from_smiles_list",
                                error.to_string(),
                            )),
                        };
                        tick_progress(progress);
                        out
                    })
                    .collect()
            });
            Self::from_records_with_mode(records, errors)
        })
    }

    pub fn read_sdf_records_from_str(
        sdf_text: &str,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_records_from_str_with_options(sdf_text, coordinate_mode, errors, None, None)
    }

    pub fn read_sdf_records_from_str_with_options(
        sdf_text: &str,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_records_from_str_with_params_and_options(
            sdf_text,
            SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
            errors,
            n_jobs,
            progress_bar,
        )
    }

    pub fn read_sdf_records_from_str_with_params_and_options(
        sdf_text: &str,
        params: SdfReadParams,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        let records = split_sdf_record_strings(sdf_text);
        with_progress_bar_for_option(
            progress_bar,
            records.len(),
            "Reading SDF records",
            |progress| {
                let batch_records = run_with_parallel_jobs_option(n_jobs, || {
                    records
                        .par_iter()
                        .enumerate()
                        .map(|(index, sdf)| {
                            let out = match read_sdf_from_str_with_params(sdf, params) {
                                Ok(record) => BatchRecord::Molecule(record.molecule),
                                Err(error) => BatchRecord::Error(BatchRecordError::new(
                                    index,
                                    "batch.read_sdf_records_from_str",
                                    error.to_string(),
                                )),
                            };
                            tick_progress(progress);
                            out
                        })
                        .collect()
                });
                Self::from_records_with_mode(batch_records, errors)
            },
        )
    }

    pub fn read_sdf_records_from_reader<R: BufRead>(
        reader: R,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_records_from_reader_with_options(reader, coordinate_mode, errors, None, None)
    }

    pub fn read_sdf_records_from_reader_with_options<R: BufRead>(
        reader: R,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_records_from_reader_with_params_and_options(
            reader,
            SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
            errors,
            n_jobs,
            progress_bar,
        )
    }

    pub fn read_sdf_records_from_reader_with_params_and_options<R: BufRead>(
        reader: R,
        params: SdfReadParams,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        let _ = n_jobs;
        with_progress_bar_for_option(progress_bar, 0, "Reading SDF records", |progress| {
            Self::read_sdf_records_from_reader_with_params_and_progress(
                reader, params, errors, progress,
            )
        })
    }

    pub fn read_sdf_records_from_reader_with_progress<R: BufRead>(
        reader: R,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_records_from_reader_with_params_and_progress(
            reader,
            SdfReadParams {
                coordinate_mode,
                ..Default::default()
            },
            errors,
            progress,
        )
    }

    pub fn read_sdf_records_from_reader_with_params_and_progress<R: BufRead>(
        reader: R,
        params: SdfReadParams,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        let mut reader = SdfReader::with_params(reader, params);
        let mut records = Vec::new();
        let mut index = 0usize;
        loop {
            match reader.next_record() {
                Ok(Some(record)) => {
                    records.push(BatchRecord::Molecule(record.molecule));
                    index += 1;
                    tick_progress(progress);
                }
                Ok(None) if reader.is_end() => break,
                Ok(None) => {
                    records.push(BatchRecord::Error(BatchRecordError::new(
                        index,
                        "read_sdf",
                        "SDF record returned no molecule before end of stream",
                    )));
                    index += 1;
                    tick_progress(progress);
                }
                Err(error) => {
                    records.push(BatchRecord::Error(BatchRecordError::new(
                        index,
                        "read_sdf",
                        error.to_string(),
                    )));
                    index += 1;
                    tick_progress(progress);
                }
            }
        }
        Self::from_records_with_mode(records, errors)
    }

    pub fn read_sdf_dataset_with_progress(
        dataset: &SdfDataset,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_dataset_with_params_and_progress(dataset, dataset.params(), errors, progress)
    }

    pub fn read_sdf_dataset_with_params_and_progress(
        dataset: &SdfDataset,
        params: SdfReadParams,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        Self::read_sdf_dataset_with_params_and_runtime(dataset, params, errors, progress, None)
    }

    pub fn read_sdf_dataset_with_params_and_options(
        dataset: &SdfDataset,
        params: SdfReadParams,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        with_progress_bar_for_option(
            progress_bar,
            dataset.len(),
            "Reading SDF dataset",
            |progress| {
                Self::read_sdf_dataset_with_params_and_runtime(
                    dataset, params, errors, progress, n_jobs,
                )
            },
        )
    }

    fn read_sdf_dataset_with_params_and_runtime(
        dataset: &SdfDataset,
        params: SdfReadParams,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Self, BatchValidationError> {
        let records = run_with_parallel_jobs_option(n_jobs, || {
            (0..dataset.len())
                .into_par_iter()
                .map(|index| {
                    let out = match dataset.record_with_params(index, params) {
                        Ok(record) => BatchRecord::Molecule(record.molecule),
                        Err(error) => BatchRecord::Error(BatchRecordError::new(
                            index,
                            "batch.read_sdf_dataset",
                            error.to_string(),
                        )),
                    };
                    tick_progress(progress);
                    out
                })
                .collect()
        });
        Self::from_records_with_mode(records, errors)
    }

    pub fn from_records_with_mode(
        records: Vec<BatchRecord>,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        if errors.raise_on_errors() {
            let record_errors = record_errors_from_records(&records);
            if !record_errors.is_empty() {
                return Err(BatchValidationError::from_record_errors(record_errors));
            }
        }
        Ok(Self {
            records,
            n_jobs: None,
            progress_bar: None,
        })
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.records.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    #[must_use]
    pub fn valid_mask(&self) -> Vec<bool> {
        self.records
            .iter()
            .map(|record| matches!(record, BatchRecord::Molecule(_)))
            .collect()
    }

    #[must_use]
    pub fn errors(&self) -> Vec<BatchRecordError> {
        self.records
            .iter()
            .filter_map(|record| match record {
                BatchRecord::Molecule(_) => None,
                BatchRecord::Error(error) => Some(error.clone()),
            })
            .collect()
    }

    #[must_use]
    pub fn valid_count(&self) -> usize {
        self.records
            .iter()
            .filter(|record| matches!(record, BatchRecord::Molecule(_)))
            .count()
    }

    #[must_use]
    pub fn invalid_count(&self) -> usize {
        self.len() - self.valid_count()
    }

    #[must_use]
    pub fn filter_valid(&self) -> Self {
        Self {
            records: self
                .records
                .iter()
                .filter(|record| matches!(record, BatchRecord::Molecule(_)))
                .cloned()
                .collect(),
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
        }
    }

    pub fn sanitize(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.sanitize_with_options(errors, None, None)
    }

    pub fn with_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.with_hydrogens_with_options(errors, None, None)
    }

    pub fn without_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.without_hydrogens_with_options(errors, None, None)
    }

    pub fn with_kekulized_bonds(
        &self,
        clear_aromatic_flags: bool,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.with_kekulized_bonds_with_options(clear_aromatic_flags, errors, None, None)
    }

    pub fn with_2d_coordinates(
        &self,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.with_2d_coordinates_with_options(errors, None, None)
    }

    pub fn with_2d_coordinates_with_params(
        &self,
        params: crate::With2DCoordinatesParams,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.with_2d_coordinates_with_params_and_options(params, errors, None, None)
    }

    pub fn sanitize_with_options(
        &self,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Sanitizing molecules",
            |progress| {
                self.transform_with_options(
                    "batch.sanitize",
                    errors,
                    Molecule::sanitize,
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn with_hydrogens_with_options(
        &self,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Adding hydrogens",
            |progress| {
                self.transform_with_options(
                    "batch.with_hydrogens",
                    errors,
                    Molecule::with_hydrogens,
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn without_hydrogens_with_options(
        &self,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Removing hydrogens",
            |progress| {
                self.transform_with_options(
                    "batch.without_hydrogens",
                    errors,
                    Molecule::without_hydrogens,
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn with_kekulized_bonds_with_options(
        &self,
        clear_aromatic_flags: bool,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Kekulizing molecules",
            |progress| {
                self.transform_with_options(
                    "batch.with_kekulized_bonds",
                    errors,
                    move |molecule| molecule.with_kekulized_bonds(clear_aromatic_flags),
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn with_2d_coordinates_with_options(
        &self,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing 2D coordinates",
            |progress| {
                self.transform_with_options(
                    "batch.with_2d_coordinates",
                    errors,
                    Molecule::with_2d_coordinates,
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn with_2d_coordinates_with_params_and_options(
        &self,
        params: crate::With2DCoordinatesParams,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Self, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing 2D coordinates",
            |progress| {
                self.transform_with_options(
                    "batch.with_2d_coordinates",
                    errors,
                    move |molecule| molecule.with_2d_coordinates_with_params(params),
                    progress,
                    n_jobs,
                )
            },
        )
    }

    pub fn with_hydrogens_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_with_options(
            "batch.with_hydrogens",
            errors,
            Molecule::with_hydrogens,
            progress,
            None,
        )
    }

    pub fn without_hydrogens_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_with_options(
            "batch.without_hydrogens",
            errors,
            Molecule::without_hydrogens,
            progress,
            None,
        )
    }

    pub fn sanitize_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_with_options("batch.sanitize", errors, Molecule::sanitize, progress, None)
    }

    pub fn with_kekulized_bonds_with_progress(
        &self,
        clear_aromatic_flags: bool,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_with_options(
            "batch.with_kekulized_bonds",
            errors,
            move |molecule| molecule.with_kekulized_bonds(clear_aromatic_flags),
            progress,
            None,
        )
    }

    pub fn with_2d_coordinates_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_with_options(
            "batch.with_2d_coordinates",
            errors,
            Molecule::with_2d_coordinates,
            progress,
            None,
        )
    }

    fn transform_with_options<F>(
        &self,
        operation: &'static str,
        errors: BatchErrorMode,
        transform: F,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Self, BatchValidationError>
    where
        F: Fn(&Molecule) -> Result<Molecule, crate::OperationError> + Sync + Send,
    {
        let records: Vec<BatchRecord> = self.run_with_parallel_jobs(n_jobs, || {
            self.records
                .par_iter()
                .enumerate()
                .map(|(index, record)| {
                    let out = match record {
                        BatchRecord::Molecule(molecule) => match transform(molecule) {
                            Ok(molecule) => BatchRecord::Molecule(molecule),
                            Err(error) => BatchRecord::Error(BatchRecordError::new(
                                index,
                                operation,
                                error.to_string(),
                            )),
                        },
                        BatchRecord::Error(error) => BatchRecord::Error(error.clone()),
                    };
                    tick_progress(progress);
                    out
                })
                .collect()
        });
        if errors.raise_on_errors() {
            let record_errors = record_errors_from_records(&records);
            if !record_errors.is_empty() {
                return Err(BatchValidationError::from_record_errors(record_errors));
            }
        }
        Ok(Self {
            records,
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
        })
    }

    pub fn to_smiles_list(
        &self,
        errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        self.to_smiles_list_with_options(errors, None, None)
    }

    pub fn to_smiles_list_with_params(
        &self,
        params: &SmilesWriteParams,
        errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        self.to_smiles_list_with_params_and_options(params, errors, None, None)
    }

    pub fn to_smiles_list_with_options(
        &self,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<String>, BatchValidationError> {
        let params = SmilesWriteParams::default();
        self.to_smiles_list_with_params_and_options(&params, errors, n_jobs, progress_bar)
    }

    pub fn to_smiles_list_with_params_and_options(
        &self,
        params: &SmilesWriteParams,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<String>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Writing SMILES",
            |progress| {
                let outcomes: Vec<Result<String, BatchRecordError>> =
                    self.run_with_parallel_jobs(n_jobs, || {
                        self.records
                            .par_iter()
                            .enumerate()
                            .map(|(index, record)| {
                                let out = match record {
                                    BatchRecord::Molecule(molecule) => {
                                        molecule.to_smiles_with_params(params).map_err(|error| {
                                            BatchRecordError::new(
                                                index,
                                                "to_smiles",
                                                error.to_string(),
                                            )
                                        })
                                    }
                                    BatchRecord::Error(error) => Err(error.clone()),
                                };
                                tick_progress(progress);
                                out
                            })
                            .collect()
                    });
                let mut results = Vec::with_capacity(outcomes.len());
                let mut record_errors = Vec::new();
                for outcome in outcomes {
                    match outcome {
                        Ok(smiles) => results.push(smiles),
                        Err(error) => {
                            record_errors.push(error);
                            results.push("?".to_string());
                        }
                    }
                }
                if errors.raise_on_errors() && !record_errors.is_empty() {
                    return Err(BatchValidationError::from_record_errors(record_errors));
                }
                Ok(results)
            },
        )
    }

    pub fn to_smiles_list_with_params_and_progress(
        &self,
        params: &SmilesWriteParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.to_smiles_optional_list_with_params_and_runtime(params, progress, None)
    }

    pub fn to_smiles_optional_list_with_params_and_options(
        &self,
        params: &SmilesWriteParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Writing SMILES",
            |progress| {
                self.to_smiles_optional_list_with_params_and_runtime(params, progress, n_jobs)
            },
        )
    }

    fn to_smiles_optional_list_with_params_and_runtime(
        &self,
        params: &SmilesWriteParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "to_smiles",
            |molecule| {
                molecule
                    .to_smiles_with_params(params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn dg_bounds_matrix_list_with_progress(
        &self,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.dg_bounds_matrix_list_with_runtime(progress, None)
    }

    pub fn dg_bounds_matrix_list_with_options(
        &self,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing distance-geometry bounds",
            |progress| self.dg_bounds_matrix_list_with_runtime(progress, n_jobs),
        )
    }

    fn dg_bounds_matrix_list_with_runtime(
        &self,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.dg_bounds_matrix",
            |molecule| {
                molecule
                    .dg_bounds_matrix()
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn atom_pair_sparse_count_fingerprint_list_with_progress(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.atom_pair_sparse_count_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn atom_pair_sparse_count_fingerprint_list_with_options(
        &self,
        params: &crate::AtomPairFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing sparse-count AtomPair fingerprints",
            |progress| {
                self.atom_pair_sparse_count_fingerprint_list_with_runtime(params, progress, n_jobs)
            },
        )
    }

    fn atom_pair_sparse_count_fingerprint_list_with_runtime(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        let generator =
            atom_pair_generator_for_batch(params, "batch.atom_pair_sparse_count_fingerprint")?;
        self.collect_optional_values_with_options(
            "batch.atom_pair_sparse_count_fingerprint",
            |molecule| {
                generator
                    .sparse_count_fingerprint(
                        molecule,
                        &mut crate::fingerprint::atom_pair_function_arguments(params),
                    )
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn atom_pair_count_fingerprint_list_with_progress(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.atom_pair_count_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn atom_pair_count_fingerprint_list_with_options(
        &self,
        params: &crate::AtomPairFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing count AtomPair fingerprints",
            |progress| self.atom_pair_count_fingerprint_list_with_runtime(params, progress, n_jobs),
        )
    }

    fn atom_pair_count_fingerprint_list_with_runtime(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        let generator = atom_pair_generator_for_batch(params, "batch.atom_pair_count_fingerprint")?;
        self.collect_optional_values_with_options(
            "batch.atom_pair_count_fingerprint",
            |molecule| {
                generator
                    .count_fingerprint(
                        molecule,
                        &mut crate::fingerprint::atom_pair_function_arguments(params),
                    )
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn atom_pair_sparse_bit_fingerprint_list_with_progress(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.atom_pair_sparse_bit_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn atom_pair_sparse_bit_fingerprint_list_with_options(
        &self,
        params: &crate::AtomPairFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing sparse-bit AtomPair fingerprints",
            |progress| {
                self.atom_pair_sparse_bit_fingerprint_list_with_runtime(params, progress, n_jobs)
            },
        )
    }

    fn atom_pair_sparse_bit_fingerprint_list_with_runtime(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        let generator =
            atom_pair_generator_for_batch(params, "batch.atom_pair_sparse_bit_fingerprint")?;
        self.collect_optional_values_with_options(
            "batch.atom_pair_sparse_bit_fingerprint",
            |molecule| {
                generator
                    .sparse_bit_fingerprint(
                        molecule,
                        &mut crate::fingerprint::atom_pair_function_arguments(params),
                    )
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn atom_pair_fingerprint_list_with_progress(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.atom_pair_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn atom_pair_fingerprint_list_with_options(
        &self,
        params: &crate::AtomPairFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing AtomPair fingerprints",
            |progress| self.atom_pair_fingerprint_list_with_runtime(params, progress, n_jobs),
        )
    }

    fn atom_pair_fingerprint_list_with_runtime(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        let generator = atom_pair_generator_for_batch(params, "batch.atom_pair_fingerprint")?;
        self.collect_optional_values_with_options(
            "batch.atom_pair_fingerprint",
            |molecule| {
                generator
                    .fingerprint(
                        molecule,
                        &mut crate::fingerprint::atom_pair_function_arguments(params),
                    )
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn atom_pair_fingerprint_with_output_list_with_progress(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::AtomPairFingerprintOutput>>, BatchValidationError> {
        self.atom_pair_fingerprint_with_output_list_with_runtime(params, progress, None)
    }

    pub fn atom_pair_fingerprint_with_output_list_with_options(
        &self,
        params: &crate::AtomPairFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::AtomPairFingerprintOutput>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing AtomPair fingerprints",
            |progress| {
                self.atom_pair_fingerprint_with_output_list_with_runtime(params, progress, n_jobs)
            },
        )
    }

    fn atom_pair_fingerprint_with_output_list_with_runtime(
        &self,
        params: &crate::AtomPairFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::AtomPairFingerprintOutput>>, BatchValidationError> {
        let generator =
            atom_pair_generator_for_batch(params, "batch.atom_pair_fingerprint_with_output")?;
        self.collect_optional_values_with_options(
            "batch.atom_pair_fingerprint_with_output",
            |molecule| {
                let mut arguments = crate::fingerprint::atom_pair_function_arguments(params);
                let fingerprint = generator
                    .fingerprint(molecule, &mut arguments)
                    .map_err(|error| error.to_string())?;
                Ok(crate::AtomPairFingerprintOutput {
                    fingerprint,
                    additional_output: arguments.additional_output,
                })
            },
            progress,
            n_jobs,
        )
    }

    pub fn morgan_fingerprint_list_with_progress(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.morgan_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn morgan_fingerprint_list_with_options(
        &self,
        params: &crate::MorganFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Morgan fingerprints",
            |progress| self.morgan_fingerprint_list_with_runtime(params, progress, n_jobs),
        )
    }

    fn morgan_fingerprint_list_with_runtime(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.morgan_fingerprint",
            |molecule| {
                molecule
                    .morgan_fingerprint(params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn morgan_fingerprint_with_output_list_with_progress(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        self.morgan_fingerprint_with_output_list_with_runtime(params, progress, None)
    }

    pub fn morgan_fingerprint_with_output_list_with_options(
        &self,
        params: &crate::MorganFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Morgan fingerprints",
            |progress| {
                self.morgan_fingerprint_with_output_list_with_runtime(params, progress, n_jobs)
            },
        )
    }

    fn morgan_fingerprint_with_output_list_with_runtime(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.morgan_fingerprint_with_output",
            |molecule| {
                molecule
                    .morgan_fingerprint_with_output(params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn topological_torsion_sparse_count_fingerprint_list_with_progress(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.topological_torsion_sparse_count_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn topological_torsion_sparse_count_fingerprint_list_with_options(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Topological Torsion sparse-count fingerprints",
            |progress| {
                self.topological_torsion_sparse_count_fingerprint_list_with_runtime(
                    params, progress, n_jobs,
                )
            },
        )
    }

    fn topological_torsion_sparse_count_fingerprint_list_with_runtime(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.topological_torsion_sparse_count_fingerprint",
            |molecule| {
                crate::topological_torsion_sparse_count_fingerprint(molecule, params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn topological_torsion_sparse_fingerprint_list_with_progress(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.topological_torsion_sparse_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn topological_torsion_sparse_fingerprint_list_with_options(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Topological Torsion sparse-bit fingerprints",
            |progress| {
                self.topological_torsion_sparse_fingerprint_list_with_runtime(
                    params, progress, n_jobs,
                )
            },
        )
    }

    fn topological_torsion_sparse_fingerprint_list_with_runtime(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseBitFingerprint>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.topological_torsion_sparse_fingerprint",
            |molecule| {
                crate::topological_torsion_sparse_fingerprint(molecule, params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn topological_torsion_count_fingerprint_list_with_progress(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.topological_torsion_count_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn topological_torsion_count_fingerprint_list_with_options(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Topological Torsion count fingerprints",
            |progress| {
                self.topological_torsion_count_fingerprint_list_with_runtime(
                    params, progress, n_jobs,
                )
            },
        )
    }

    fn topological_torsion_count_fingerprint_list_with_runtime(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::SparseCountFingerprint>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.topological_torsion_count_fingerprint",
            |molecule| {
                crate::topological_torsion_count_fingerprint(molecule, params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn topological_torsion_fingerprint_list_with_progress(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.topological_torsion_fingerprint_list_with_runtime(params, progress, None)
    }

    pub fn topological_torsion_fingerprint_list_with_options(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Topological Torsion fingerprints",
            |progress| {
                self.topological_torsion_fingerprint_list_with_runtime(params, progress, n_jobs)
            },
        )
    }

    fn topological_torsion_fingerprint_list_with_runtime(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.topological_torsion_fingerprint",
            |molecule| {
                crate::topological_torsion_fingerprint(molecule, params)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn topological_torsion_fingerprint_with_output_list_with_progress(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        request: crate::TopologicalTorsionFingerprintOutputRequest,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::TopologicalTorsionFingerprintResult>>, BatchValidationError> {
        self.topological_torsion_fingerprint_with_output_list_with_runtime(
            params, request, progress, None,
        )
    }

    pub fn topological_torsion_fingerprint_with_output_list_with_options(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        request: crate::TopologicalTorsionFingerprintOutputRequest,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<crate::TopologicalTorsionFingerprintResult>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Computing Topological Torsion fingerprints with provenance",
            |progress| {
                self.topological_torsion_fingerprint_with_output_list_with_runtime(
                    params, request, progress, n_jobs,
                )
            },
        )
    }

    fn topological_torsion_fingerprint_with_output_list_with_runtime(
        &self,
        params: &crate::TopologicalTorsionFingerprintParams,
        request: crate::TopologicalTorsionFingerprintOutputRequest,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<crate::TopologicalTorsionFingerprintResult>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.topological_torsion_fingerprint_with_output",
            |molecule| {
                crate::topological_torsion_fingerprint_with_output(molecule, params, request)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    pub fn to_svg_list_with_options(
        &self,
        width: u32,
        height: u32,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Drawing SVG molecules",
            |progress| self.to_svg_list_with_runtime(width, height, progress, n_jobs),
        )
    }

    pub fn to_svg_list_with_progress(
        &self,
        width: u32,
        height: u32,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.to_svg_list_with_runtime(width, height, progress, None)
    }

    fn to_svg_list_with_runtime(
        &self,
        width: u32,
        height: u32,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.to_svg",
            |molecule| {
                molecule
                    .to_svg(width, height)
                    .map_err(|error| error.to_string())
            },
            progress,
            n_jobs,
        )
    }

    #[cfg(test)]
    pub(crate) fn prepare_for_drawing_parity_list(
        &self,
    ) -> Result<Vec<Option<PreparedDrawMolecule>>, BatchValidationError> {
        self.collect_optional_values_with_options(
            "batch.prepare_for_drawing_parity",
            |molecule| {
                molecule
                    .prepared_for_drawing_parity()
                    .map_err(|error| error.to_string())
            },
            None,
            None,
        )
    }

    fn collect_optional_values_with_options<T, F>(
        &self,
        operation: &'static str,
        collect: F,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<Vec<Option<T>>, BatchValidationError>
    where
        T: Send,
        F: Fn(&Molecule) -> Result<T, String> + Sync + Send,
    {
        let pairs: Vec<(Option<T>, Option<BatchRecordError>)> =
            self.run_with_parallel_jobs(n_jobs, || {
                self.records
                    .par_iter()
                    .enumerate()
                    .map(|(index, record)| {
                        let out = match record {
                            BatchRecord::Molecule(molecule) => match collect(molecule) {
                                Ok(value) => (Some(value), None),
                                Err(error) => {
                                    (None, Some(BatchRecordError::new(index, operation, error)))
                                }
                            },
                            BatchRecord::Error(_) => (None, None),
                        };
                        tick_progress(progress);
                        out
                    })
                    .collect()
            });
        let mut values = Vec::with_capacity(pairs.len());
        let mut record_errors = Vec::new();
        for (value, error) in pairs {
            values.push(value);
            if let Some(error) = error {
                record_errors.push(error);
            }
        }
        if !record_errors.is_empty() {
            return Err(BatchValidationError::from_record_errors(record_errors));
        }
        Ok(values)
    }

    pub fn write_images(
        &self,
        output_dir: impl AsRef<Path>,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_images_with_options(
            output_dir.as_ref(),
            "png",
            500,
            500,
            errors,
            filenames,
            None,
            None,
        )
    }

    pub fn write_images_with_options(
        &self,
        output_dir: &Path,
        format: &str,
        width: u32,
        height: u32,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Writing molecule images",
            |progress| {
                self.write_images_with_runtime(
                    output_dir, format, width, height, errors, filenames, progress, n_jobs,
                )
            },
        )
    }

    pub fn write_images_with_progress(
        &self,
        output_dir: &Path,
        format: &str,
        width: u32,
        height: u32,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_images_with_runtime(
            output_dir, format, width, height, errors, filenames, progress, None,
        )
    }

    fn write_images_with_runtime(
        &self,
        output_dir: &Path,
        format: &str,
        width: u32,
        height: u32,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        fs::create_dir_all(output_dir).map_err(|_| {
            BatchValidationError::unsupported("batch.write_images", "directory creation failed")
        })?;
        let paths = output_paths(output_dir, self.records.len(), format, filenames)?;
        let outcomes: Vec<Result<bool, BatchRecordError>> =
            self.run_with_parallel_jobs(n_jobs, || {
                self.records
                    .par_iter()
                    .enumerate()
                    .map(|(index, record)| {
                        let out = match record {
                            BatchRecord::Molecule(molecule) => {
                                write_one_image(molecule, &paths[index], format, width, height)
                                    .map(|()| true)
                                    .map_err(|error| {
                                        BatchRecordError::new(index, "batch.write_images", error)
                                    })
                            }
                            BatchRecord::Error(_) => Ok(false),
                        };
                        tick_progress(progress);
                        out
                    })
                    .collect()
            });
        let mut written = 0usize;
        let mut skipped = 0usize;
        let mut record_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(true) => written += 1,
                Ok(false) => skipped += 1,
                Err(error) => record_errors.push(error),
            }
        }
        if errors.raise_on_errors() && (skipped != 0 || !record_errors.is_empty()) {
            let mut all_errors = record_errors.clone();
            all_errors.extend(record_errors_from_records(&self.records));
            return Err(BatchValidationError::from_record_errors(all_errors));
        }
        Ok(BatchExportReport {
            written,
            skipped,
            failed: record_errors.len(),
            errors: record_errors,
        })
    }

    pub fn write_sdf(
        &self,
        output_path: impl AsRef<Path>,
        errors: BatchErrorMode,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_sdf_with_options(output_path.as_ref(), SdfFormat::V2000, errors, None, None)
    }

    pub fn write_sdf_with_options(
        &self,
        output_path: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Writing SDF records",
            |progress| self.write_sdf_with_runtime(output_path, format, errors, progress, n_jobs),
        )
    }

    pub fn write_sdf_with_progress(
        &self,
        output_path: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_sdf_with_runtime(output_path, format, errors, progress, None)
    }

    fn write_sdf_with_runtime(
        &self,
        output_path: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let outcomes: Vec<Result<String, BatchRecordError>> =
            self.run_with_parallel_jobs(n_jobs, || {
                self.records
                    .par_iter()
                    .enumerate()
                    .map(|(index, record)| {
                        let out = match record {
                            BatchRecord::Molecule(molecule) => {
                                molecule_to_sdf_record_string(molecule, format).map_err(|error| {
                                    BatchRecordError::new(index, "batch.write_sdf", error)
                                })
                            }
                            BatchRecord::Error(error) => Err(error.clone()),
                        };
                        tick_progress(progress);
                        out
                    })
                    .collect()
            });
        let mut blocks = Vec::new();
        let mut record_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(block) => blocks.push(block),
                Err(error) => record_errors.push(error),
            }
        }
        if errors.raise_on_errors() && !record_errors.is_empty() {
            return Err(BatchValidationError::from_record_errors(record_errors));
        }
        if let Some(parent) = output_path.parent() {
            fs::create_dir_all(parent).map_err(|_| {
                BatchValidationError::unsupported("batch.write_sdf", "directory creation failed")
            })?;
        }
        let mut file = File::create(output_path).map_err(|_| {
            BatchValidationError::unsupported("batch.write_sdf", "file create failed")
        })?;
        for block in &blocks {
            file.write_all(block.as_bytes()).map_err(|_| {
                BatchValidationError::unsupported("batch.write_sdf", "file write failed")
            })?;
        }
        Ok(BatchExportReport {
            written: blocks.len(),
            skipped: record_errors.len(),
            failed: 0,
            errors: Vec::new(),
        })
    }

    pub fn write_sdf_files_with_progress(
        &self,
        output_dir: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_sdf_files_with_runtime(output_dir, format, errors, filenames, progress, None)
    }

    pub fn write_sdf_files_with_options(
        &self,
        output_dir: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.with_progress_bar_for(
            progress_bar,
            self.records.len(),
            "Writing SDF files",
            |progress| {
                self.write_sdf_files_with_runtime(
                    output_dir, format, errors, filenames, progress, n_jobs,
                )
            },
        )
    }

    fn write_sdf_files_with_runtime(
        &self,
        output_dir: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
        progress: BatchProgress<'_>,
        n_jobs: Option<usize>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        fs::create_dir_all(output_dir).map_err(|_| {
            BatchValidationError::unsupported("batch.write_sdf_files", "directory creation failed")
        })?;
        let paths = output_paths(output_dir, self.records.len(), "sdf", filenames)?;
        let outcomes: Vec<Result<bool, BatchRecordError>> =
            self.run_with_parallel_jobs(n_jobs, || {
                self.records
                    .par_iter()
                    .enumerate()
                    .map(|(index, record)| {
                        let out = match record {
                            BatchRecord::Molecule(molecule) => {
                                write_one_sdf_file(molecule, &paths[index], format)
                                    .map(|()| true)
                                    .map_err(|error| {
                                        BatchRecordError::new(index, "batch.write_sdf_files", error)
                                    })
                            }
                            BatchRecord::Error(_) => Ok(false),
                        };
                        tick_progress(progress);
                        out
                    })
                    .collect()
            });
        let mut written = 0usize;
        let mut skipped = 0usize;
        let mut record_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(true) => written += 1,
                Ok(false) => skipped += 1,
                Err(error) => record_errors.push(error),
            }
        }
        if errors.raise_on_errors() && (skipped != 0 || !record_errors.is_empty()) {
            let mut all_errors = record_errors.clone();
            all_errors.extend(record_errors_from_records(&self.records));
            return Err(BatchValidationError::from_record_errors(all_errors));
        }
        Ok(BatchExportReport {
            written,
            skipped,
            failed: record_errors.len(),
            errors: record_errors,
        })
    }
}

fn write_one_image(
    molecule: &Molecule,
    path: &Path,
    format: &str,
    width: u32,
    height: u32,
) -> Result<(), String> {
    match format {
        "svg" => {
            let svg = molecule
                .to_svg(width, height)
                .map_err(|error| error.to_string())?;
            fs::write(path, svg).map_err(|error| error.to_string())
        }
        "png" => {
            let png = molecule
                .to_png(width, height)
                .map_err(|error| error.to_string())?;
            fs::write(path, png).map_err(|error| error.to_string())
        }
        other => Err(format!(
            "unsupported image format '{other}', expected 'png' or 'svg'"
        )),
    }
}

fn write_one_sdf_file(molecule: &Molecule, path: &Path, format: SdfFormat) -> Result<(), String> {
    let block = molecule_to_sdf_record_string(molecule, format)?;
    fs::write(path, block).map_err(|error| error.to_string())
}

fn molecule_to_sdf_record_string(molecule: &Molecule, format: SdfFormat) -> Result<String, String> {
    let params = molblock::MolBlockWriteParams {
        format,
        force_2d: molecule.coordinates_2d().is_some(),
        ..Default::default()
    };
    molblock::mol_to_sdf_record_with_params(molecule, &params).map_err(|error| error.to_string())
}

fn output_paths(
    out_dir: &Path,
    total: usize,
    extension: &str,
    filenames: Option<&[String]>,
) -> Result<Vec<PathBuf>, BatchValidationError> {
    if let Some(filenames) = filenames
        && filenames.len() != total
    {
        return Err(BatchValidationError::unsupported(
            "batch.output_paths",
            "filenames length must match batch length",
        ));
    }
    let mut seen = HashSet::new();
    let mut paths = Vec::with_capacity(total);
    for index in 0..total {
        let filename = match filenames.and_then(|names| names.get(index)) {
            Some(raw) => normalize_output_filename(raw, extension).map_err(|_| {
                BatchValidationError::unsupported("batch.output_paths", "invalid filename")
            })?,
            None => format!("mol_{index}.{extension}"),
        };
        if !seen.insert(filename.clone()) {
            return Err(BatchValidationError::unsupported(
                "batch.output_paths",
                "duplicate output filename",
            ));
        }
        paths.push(out_dir.join(filename));
    }
    Ok(paths)
}

fn normalize_output_filename(raw: &str, extension: &str) -> Result<String, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("filename must not be empty".to_string());
    }
    let path = Path::new(trimmed);
    if path.is_absolute() {
        return Err("filename must be relative to the output directory".to_string());
    }
    let components = path.components().collect::<Vec<_>>();
    if components.len() != 1 || !matches!(components[0], Component::Normal(_)) {
        return Err("filename must not include path separators or '..'".to_string());
    }
    let file_name = path
        .file_name()
        .and_then(|value| value.to_str())
        .ok_or_else(|| "filename must be valid UTF-8".to_string())?;
    match path.extension().and_then(|value| value.to_str()) {
        Some(actual) if actual.eq_ignore_ascii_case(extension) => Ok(file_name.to_string()),
        Some(actual) => Err(format!(
            "filename extension '.{actual}' does not match expected '.{extension}'"
        )),
        None => Ok(format!("{file_name}.{extension}")),
    }
}

fn split_sdf_record_strings(sdf_text: &str) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    for line in sdf_text.lines() {
        current.push_str(line);
        current.push('\n');
        if line.trim_end() == "$$$$" {
            records.push(std::mem::take(&mut current));
        }
    }
    if !current.trim().is_empty() {
        records.push(current);
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_smiles_list_preserves_order_and_keeps_record_errors() {
        let smiles = vec!["CCO".to_string(), "C1".to_string(), "N".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles);
        assert_eq!(batch.len(), 3);
        assert_eq!(batch.valid_mask(), vec![true, false, true]);
        assert_eq!(batch.valid_count(), 2);
        assert_eq!(batch.invalid_count(), 1);
        let errors = batch.errors();
        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].index, 1);
        assert_eq!(errors[0].operation, "batch.from_smiles_list");
    }

    #[test]
    fn filter_valid_preserves_surviving_input_order() {
        let smiles = vec!["C".to_string(), "C1".to_string(), "O".to_string()];
        let filtered = MoleculeBatch::from_smiles_list(&smiles).filter_valid();
        assert_eq!(filtered.len(), 2);
        assert_eq!(filtered.valid_mask(), vec![true, true]);
    }

    #[test]
    fn read_sdf_records_from_reader_keeps_recovered_supplier_none_records() {
        let sdf = "\
ok
  COSMolKit      2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
bad


  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
$$$$
after
  COSMolKit      2D

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
";

        let kept = MoleculeBatch::read_sdf_records_from_reader(
            std::io::Cursor::new(sdf.as_bytes()),
            SdfCoordinateMode::Preserve,
            BatchErrorMode::KeepErrors,
        )
        .unwrap();

        assert_eq!(kept.len(), 3);
        assert_eq!(kept.valid_mask(), vec![true, false, true]);
        let errors = kept.errors();
        assert_eq!(errors.len(), 1);
        assert_eq!(errors[0].index, 1);
        assert_eq!(errors[0].operation, "read_sdf");

        let strict = MoleculeBatch::read_sdf_records_from_reader(
            std::io::Cursor::new(sdf.as_bytes()),
            SdfCoordinateMode::Preserve,
            BatchErrorMode::Strict,
        )
        .unwrap_err();
        assert_eq!(strict.errors, 1);
        assert_eq!(strict.record_errors[0].index, 1);
    }

    #[test]
    fn batch_configuration_defaults_are_none() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]);
        assert_eq!(batch.parallel_jobs(), None);
        assert_eq!(batch.progress_bar(), None);
    }

    #[test]
    fn batch_configuration_can_be_set_and_cleared() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()])
            .with_parallel_jobs(Some(4))
            .with_progress_bar(Some(true));
        assert_eq!(batch.parallel_jobs(), Some(4));
        assert_eq!(batch.progress_bar(), Some(true));

        let cleared = batch.with_parallel_jobs(None).with_progress_bar(None);
        assert_eq!(cleared.parallel_jobs(), None);
        assert_eq!(cleared.progress_bar(), None);
    }

    #[test]
    fn batch_configuration_is_preserved_across_transforms() {
        let batch = MoleculeBatch::from_smiles_list(&["CC".to_string()])
            .with_parallel_jobs(Some(2))
            .with_progress_bar(Some(false));

        let transformed = batch
            .with_2d_coordinates(BatchErrorMode::Strict)
            .expect("2d coordinates should succeed");

        assert_eq!(transformed.parallel_jobs(), Some(2));
        assert_eq!(transformed.progress_bar(), Some(false));
    }

    #[test]
    fn batch_configuration_is_preserved_across_filter_valid() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "C1".to_string()])
            .with_parallel_jobs(Some(8))
            .with_progress_bar(Some(true));

        let filtered = batch.filter_valid();

        assert_eq!(filtered.parallel_jobs(), Some(8));
        assert_eq!(filtered.progress_bar(), Some(true));
    }

    #[test]
    fn transform_options_preserve_batch_configuration() {
        let batch = MoleculeBatch::from_smiles_list(&["CC".to_string()])
            .with_parallel_jobs(Some(4))
            .with_progress_bar(Some(true));

        let transformed = batch
            .with_2d_coordinates_with_options(BatchErrorMode::Strict, Some(1), Some(false))
            .expect("2d coordinates should succeed");

        assert_eq!(transformed.parallel_jobs(), Some(4));
        assert_eq!(transformed.progress_bar(), Some(true));
    }

    #[test]
    fn smiles_writer_options_use_batch_runtime_overrides_without_mutating_batch_defaults() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "CC".to_string()])
            .with_parallel_jobs(Some(2))
            .with_progress_bar(Some(true));

        let smiles = batch
            .to_smiles_list_with_options(BatchErrorMode::Strict, Some(1), Some(false))
            .expect("smiles writing should succeed");

        assert_eq!(smiles, vec!["CCO".to_string(), "CC".to_string()]);
        assert_eq!(batch.parallel_jobs(), Some(2));
        assert_eq!(batch.progress_bar(), Some(true));
    }

    #[test]
    fn svg_collection_options_use_batch_runtime_overrides_without_mutating_batch_defaults() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()])
            .with_progress_bar(Some(true))
            .with_parallel_jobs(Some(2));

        let prepared = batch
            .with_2d_coordinates_with_options(BatchErrorMode::Strict, Some(1), Some(false))
            .expect("2d coordinates should succeed");
        let svgs = prepared
            .to_svg_list_with_options(320, 240, Some(1), Some(false))
            .expect("svg export should succeed");

        assert_eq!(svgs.len(), 1);
        assert!(svgs[0].as_ref().is_some_and(|svg| svg.contains("<svg")));
        assert_eq!(prepared.parallel_jobs(), Some(2));
        assert_eq!(prepared.progress_bar(), Some(true));
    }
}
