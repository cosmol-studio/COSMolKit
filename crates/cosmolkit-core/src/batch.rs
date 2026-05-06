use crate::io::molblock::{self, SdfFormat};
use crate::io::sdf::{SdfCoordinateMode, read_sdf_from_str_with_coordinate_mode};
use crate::{Molecule, PreparedDrawMolecule, SmilesWriteParams};
use indicatif::{ProgressBar, ProgressDrawTarget, ProgressStyle};
use rayon::prelude::*;
use std::collections::HashSet;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Component, Path, PathBuf};
use thiserror::Error;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub enum BatchErrorMode {
    Raise,
    Keep,
    Skip,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchRecordError {
    pub index: usize,
    pub input: Option<String>,
    pub stage: String,
    pub error_type: String,
    pub message: String,
}

impl BatchRecordError {
    #[must_use]
    pub fn new(
        index: usize,
        input: Option<String>,
        stage: impl Into<String>,
        error_type: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            index,
            input,
            stage: stage.into(),
            error_type: error_type.into(),
            message: message.into(),
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum BatchRecord {
    Valid(Molecule),
    Invalid(BatchRecordError),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MoleculeBatch {
    pub records: Vec<BatchRecord>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BatchExportReport {
    pub total: usize,
    pub success: usize,
    pub failed: usize,
    pub errors: Vec<BatchRecordError>,
}

#[derive(Debug, Clone, PartialEq, Eq, Error)]
#[error("batch validation failed with {} error(s)", .errors.len())]
pub struct BatchValidationError {
    pub errors: Vec<BatchRecordError>,
}

/// Optional per-record progress callback used by parallel batch operations.
pub type BatchProgress<'a> = Option<&'a (dyn Fn() + Sync)>;

/// Terminal progress bar for Rust batch workflows.
pub struct BatchProgressBar {
    inner: ProgressBar,
}

impl BatchProgressBar {
    /// Create a stderr progress bar that adapts to terminal width.
    #[must_use]
    pub fn new(total: usize, message: impl Into<String>) -> Self {
        let progress_bar = ProgressBar::with_draw_target(
            Some(total as u64),
            ProgressDrawTarget::term_like_with_hz(Box::new(console::Term::buffered_stderr()), 20),
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

    /// Return a callback suitable for `*_with_progress()` batch methods.
    #[must_use]
    pub fn callback(&self) -> Box<dyn Fn() + Sync + '_> {
        Box::new(|| self.inner.inc(1))
    }

    /// Finish rendering the progress bar and emit a newline.
    pub fn finish(&self) {
        self.inner.finish();
        eprintln!();
    }
}

/// Create a progress bar when `enabled` is true and `total` is non-zero.
#[must_use]
pub fn batch_progress_bar(
    enabled: bool,
    total: usize,
    message: impl Into<String>,
) -> Option<BatchProgressBar> {
    if !enabled || total == 0 {
        return None;
    }
    Some(BatchProgressBar::new(total, message))
}

fn tick_progress(progress: BatchProgress<'_>) {
    if let Some(progress) = progress {
        progress();
    }
}

impl MoleculeBatch {
    #[must_use]
    pub fn new(records: Vec<BatchRecord>) -> Self {
        Self { records }
    }

    pub fn from_smiles_list(
        smiles: &[String],
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        Self::from_smiles_list_with_sanitize(smiles, true, errors)
    }

    pub fn from_smiles_list_with_sanitize(
        smiles: &[String],
        sanitize: bool,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let records: Vec<BatchRecord> = smiles
            .par_iter()
            .enumerate()
            .filter_map(|(index, smiles)| {
                match Molecule::from_smiles_with_sanitize(smiles, sanitize) {
                    Ok(molecule) => Some(BatchRecord::Valid(molecule)),
                    Err(error) => {
                        let record_error = BatchRecordError::new(
                            index,
                            Some(smiles.clone()),
                            "parse_smiles",
                            "SmilesParseError",
                            error.to_string(),
                        );
                        match errors {
                            BatchErrorMode::Raise | BatchErrorMode::Keep => {
                                Some(BatchRecord::Invalid(record_error))
                            }
                            BatchErrorMode::Skip => None,
                        }
                    }
                }
            })
            .collect();
        Self::from_records_with_mode(records, errors)
    }

    fn from_sdf_record_strings(
        records: &[String],
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let batch_records: Vec<BatchRecord> = records
            .par_iter()
            .enumerate()
            .filter_map(|(index, sdf)| {
                match read_sdf_from_str_with_coordinate_mode(sdf, coordinate_mode) {
                    Ok(record) => Some(BatchRecord::Valid(record.molecule)),
                    Err(error) => {
                        let record_error = BatchRecordError::new(
                            index,
                            None,
                            "read_sdf",
                            "SdfReadError",
                            error.to_string(),
                        );
                        match errors {
                            BatchErrorMode::Raise | BatchErrorMode::Keep => {
                                Some(BatchRecord::Invalid(record_error))
                            }
                            BatchErrorMode::Skip => None,
                        }
                    }
                }
            })
            .collect();
        Self::from_records_with_mode(batch_records, errors)
    }

    pub fn read_sdf_records_from_str(
        sdf_text: &str,
        coordinate_mode: SdfCoordinateMode,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        let records = split_sdf_record_strings(sdf_text);
        Self::from_sdf_record_strings(&records, coordinate_mode, errors)
    }

    fn from_records_with_mode(
        records: Vec<BatchRecord>,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        if matches!(errors, BatchErrorMode::Raise) {
            let collected = records
                .iter()
                .filter_map(|record| match record {
                    BatchRecord::Valid(_) => None,
                    BatchRecord::Invalid(error) => Some(error.clone()),
                })
                .collect::<Vec<_>>();
            if !collected.is_empty() {
                return Err(BatchValidationError { errors: collected });
            }
        }
        Ok(Self { records })
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
            .map(|record| matches!(record, BatchRecord::Valid(_)))
            .collect()
    }

    #[must_use]
    pub fn errors(&self) -> Vec<BatchRecordError> {
        self.records
            .iter()
            .filter_map(|record| match record {
                BatchRecord::Valid(_) => None,
                BatchRecord::Invalid(error) => Some(error.clone()),
            })
            .collect()
    }

    #[must_use]
    pub fn valid_count(&self) -> usize {
        self.records
            .iter()
            .filter(|record| matches!(record, BatchRecord::Valid(_)))
            .count()
    }

    #[must_use]
    pub fn invalid_count(&self) -> usize {
        self.records.len() - self.valid_count()
    }

    #[must_use]
    pub fn filter_valid(&self) -> Self {
        let records = self
            .records
            .iter()
            .filter_map(|record| match record {
                BatchRecord::Valid(molecule) => Some(BatchRecord::Valid(molecule.clone())),
                BatchRecord::Invalid(_) => None,
            })
            .collect();
        Self { records }
    }

    pub fn add_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.add_hydrogens_with_progress(errors, None)
    }

    pub fn add_hydrogens_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_valid(
            "add_hydrogens",
            "AddHydrogensError",
            errors,
            |molecule| {
                molecule
                    .with_hydrogens()
                    .map_err(|error| format!("{error:?}"))
            },
            progress,
        )
    }

    pub fn remove_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.remove_hydrogens_with_progress(errors, None)
    }

    pub fn remove_hydrogens_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_valid(
            "remove_hydrogens",
            "RemoveHydrogensError",
            errors,
            |molecule| {
                molecule
                    .without_hydrogens()
                    .map_err(|error| format!("{error:?}"))
            },
            progress,
        )
    }

    pub fn sanitize(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.sanitize_with_progress(errors, None)
    }

    pub fn sanitize_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_valid(
            "sanitize",
            "SanitizeError",
            errors,
            |molecule| molecule.sanitize().map_err(|error| error.to_string()),
            progress,
        )
    }

    pub fn kekulize(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.kekulize_with_sanitize(false, errors)
    }

    pub fn kekulize_with_sanitize(
        &self,
        sanitize: bool,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.kekulize_with_sanitize_and_progress(sanitize, errors, None)
    }

    pub fn kekulize_with_sanitize_and_progress(
        &self,
        sanitize: bool,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_valid(
            "kekulize",
            "KekulizeError",
            errors,
            |molecule| {
                molecule
                    .with_kekulized_bonds(sanitize)
                    .map_err(|error| format!("{error:?}"))
            },
            progress,
        )
    }

    pub fn compute_2d_coords(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.compute_2d_coords_with_progress(errors, None)
    }

    pub fn compute_2d_coords_with_progress(
        &self,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError> {
        self.transform_valid(
            "compute_2d_coords",
            "CoordinateGenerationError",
            errors,
            |molecule| molecule.with_2d_coords().map_err(|error| error.to_string()),
            progress,
        )
    }

    pub fn to_smiles_list(
        &self,
        isomeric_smiles: bool,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.collect_optional_values(
            "to_smiles",
            "SmilesWriteError",
            |molecule| {
                molecule
                    .to_smiles(isomeric_smiles)
                    .map_err(|error| error.to_string())
            },
            None,
        )
    }

    pub fn to_smiles_list_with_params(
        &self,
        params: &SmilesWriteParams,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.to_smiles_list_with_params_and_progress(params, None)
    }

    pub fn to_smiles_list_with_params_and_progress(
        &self,
        params: &SmilesWriteParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.collect_optional_values(
            "to_smiles",
            "SmilesWriteError",
            |molecule| {
                molecule
                    .to_smiles_with_params(params)
                    .map_err(|error| error.to_string())
            },
            progress,
        )
    }

    pub fn dg_bounds_matrix_list(
        &self,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.dg_bounds_matrix_list_with_progress(None)
    }

    pub fn dg_bounds_matrix_list_with_progress(
        &self,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<Vec<Vec<f64>>>>, BatchValidationError> {
        self.collect_optional_values(
            "dg_bounds_matrix",
            "DistanceGeometryError",
            |molecule| {
                molecule
                    .dg_bounds_matrix()
                    .map_err(|error| error.to_string())
            },
            progress,
        )
    }

    pub fn morgan_fingerprint_list(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.morgan_fingerprint_list_with_progress(params, None)
    }

    pub fn morgan_fingerprint_list_with_progress(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::Fingerprint>>, BatchValidationError> {
        self.collect_optional_values(
            "morgan_fingerprint",
            "FingerprintError",
            |molecule| {
                molecule
                    .morgan_fingerprint(params)
                    .map_err(|error| error.to_string())
            },
            progress,
        )
    }

    pub fn morgan_fingerprint_with_output_list(
        &self,
        params: &crate::MorganFingerprintParams,
    ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        self.morgan_fingerprint_with_output_list_with_progress(params, None)
    }

    pub fn morgan_fingerprint_with_output_list_with_progress(
        &self,
        params: &crate::MorganFingerprintParams,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<crate::MorganFingerprintOutput>>, BatchValidationError> {
        self.collect_optional_values(
            "morgan_fingerprint",
            "FingerprintError",
            |molecule| {
                molecule
                    .morgan_fingerprint_with_output(params)
                    .map_err(|error| error.to_string())
            },
            progress,
        )
    }

    pub fn to_svg_list(
        &self,
        width: u32,
        height: u32,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.to_svg_list_with_progress(width, height, None)
    }

    pub fn to_svg_list_with_progress(
        &self,
        width: u32,
        height: u32,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<String>>, BatchValidationError> {
        self.collect_optional_values(
            "to_svg",
            "SvgDrawError",
            |molecule| {
                molecule
                    .to_svg(width, height)
                    .map_err(|error| error.to_string())
            },
            progress,
        )
    }

    pub fn prepare_for_drawing_parity_list(
        &self,
    ) -> Result<Vec<Option<PreparedDrawMolecule>>, BatchValidationError> {
        self.collect_optional_values(
            "prepare_for_drawing_parity",
            "PreparedDrawError",
            |molecule| {
                molecule
                    .prepare_for_drawing_parity()
                    .map_err(|error| error.to_string())
            },
            None,
        )
    }

    fn transform_valid<F>(
        &self,
        stage: &'static str,
        error_type: &'static str,
        errors: BatchErrorMode,
        transform: F,
        progress: BatchProgress<'_>,
    ) -> Result<Self, BatchValidationError>
    where
        F: Fn(&Molecule) -> Result<Molecule, String> + Sync,
    {
        let records: Vec<BatchRecord> = self
            .records
            .par_iter()
            .enumerate()
            .filter_map(|(index, record)| {
                let out = match record {
                    BatchRecord::Valid(molecule) => match transform(molecule) {
                        Ok(molecule) => Some(BatchRecord::Valid(molecule)),
                        Err(message) => {
                            let error =
                                BatchRecordError::new(index, None, stage, error_type, message);
                            match errors {
                                BatchErrorMode::Raise | BatchErrorMode::Keep => {
                                    Some(BatchRecord::Invalid(error))
                                }
                                BatchErrorMode::Skip => None,
                            }
                        }
                    },
                    BatchRecord::Invalid(error) => match errors {
                        BatchErrorMode::Raise | BatchErrorMode::Keep => {
                            Some(BatchRecord::Invalid(error.clone()))
                        }
                        BatchErrorMode::Skip => None,
                    },
                };
                tick_progress(progress);
                out
            })
            .collect();
        Self::from_records_with_mode(records, errors)
    }

    fn collect_optional_values<T, F>(
        &self,
        stage: &'static str,
        error_type: &'static str,
        collect: F,
        progress: BatchProgress<'_>,
    ) -> Result<Vec<Option<T>>, BatchValidationError>
    where
        T: Send,
        F: Fn(&Molecule) -> Result<T, String> + Sync,
    {
        let pairs: Vec<(Option<T>, Option<BatchRecordError>)> = self
            .records
            .par_iter()
            .enumerate()
            .map(|(index, record)| {
                let out = match record {
                    BatchRecord::Valid(molecule) => match collect(molecule) {
                        Ok(value) => (Some(value), None),
                        Err(message) => (
                            None,
                            Some(BatchRecordError::new(
                                index, None, stage, error_type, message,
                            )),
                        ),
                    },
                    BatchRecord::Invalid(_) => (None, None),
                };
                tick_progress(progress);
                out
            })
            .collect();
        let mut values = Vec::with_capacity(pairs.len());
        let mut errors = Vec::new();
        for (value, error) in pairs {
            values.push(value);
            if let Some(error) = error {
                errors.push(error);
            }
        }
        if errors.is_empty() {
            Ok(values)
        } else {
            Err(BatchValidationError { errors })
        }
    }

    pub fn write_images(
        &self,
        out_dir: &Path,
        format: &str,
        width: u32,
        height: u32,
        errors: BatchErrorMode,
        filenames: Option<&[Option<String>]>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_images_with_progress(out_dir, format, width, height, errors, filenames, None)
    }

    pub fn write_images_with_progress(
        &self,
        out_dir: &Path,
        format: &str,
        width: u32,
        height: u32,
        errors: BatchErrorMode,
        filenames: Option<&[Option<String>]>,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let format = format.to_ascii_lowercase();
        let paths = output_paths(out_dir, self.records.len(), &format, filenames, "to_images")?;
        fs::create_dir_all(out_dir).map_err(|error| BatchValidationError {
            errors: vec![BatchRecordError::new(
                0,
                Some(out_dir.display().to_string()),
                "to_images",
                "IoError",
                format!("create output directory failed: {error}"),
            )],
        })?;
        let outcomes: Vec<Result<(), BatchRecordError>> = self
            .records
            .par_iter()
            .enumerate()
            .filter_map(|(index, record)| {
                let out = match record {
                    BatchRecord::Valid(molecule) => Some(
                        write_one_image(molecule, &paths[index], &format, width, height).map_err(
                            |message| {
                                BatchRecordError::new(
                                    index,
                                    Some(paths[index].display().to_string()),
                                    "to_images",
                                    "ImageExportError",
                                    message,
                                )
                            },
                        ),
                    ),
                    BatchRecord::Invalid(error) => match errors {
                        BatchErrorMode::Raise | BatchErrorMode::Keep => Some(Err(error.clone())),
                        BatchErrorMode::Skip => None,
                    },
                };
                tick_progress(progress);
                out
            })
            .collect();
        let mut success = 0usize;
        let mut report_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(()) => success += 1,
                Err(error) => report_errors.push(error),
            }
        }
        if matches!(errors, BatchErrorMode::Raise) && !report_errors.is_empty() {
            return Err(BatchValidationError {
                errors: report_errors,
            });
        }
        Ok(BatchExportReport {
            total: self.records.len(),
            success,
            failed: report_errors.len(),
            errors: report_errors,
        })
    }

    pub fn write_sdf(
        &self,
        path: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_sdf_with_progress(path, format, errors, None)
    }

    pub fn write_sdf_with_progress(
        &self,
        path: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let outcomes: Vec<Result<String, BatchRecordError>> = self
            .records
            .par_iter()
            .enumerate()
            .filter_map(|(index, record)| {
                let out = match record {
                    BatchRecord::Valid(molecule) => Some(
                        molecule_to_sdf_record_string(molecule, format).map_err(|error| {
                            BatchRecordError::new(index, None, "to_sdf", "SdfWriteError", error)
                        }),
                    ),
                    BatchRecord::Invalid(error) => match errors {
                        BatchErrorMode::Raise | BatchErrorMode::Keep => Some(Err(error.clone())),
                        BatchErrorMode::Skip => None,
                    },
                };
                tick_progress(progress);
                out
            })
            .collect();
        let mut blocks = Vec::new();
        let mut report_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(block) => blocks.push(block),
                Err(error) => report_errors.push(error),
            }
        }
        if matches!(errors, BatchErrorMode::Raise) && !report_errors.is_empty() {
            return Err(BatchValidationError {
                errors: report_errors,
            });
        }
        let mut file = File::create(path).map_err(|error| BatchValidationError {
            errors: vec![BatchRecordError::new(
                0,
                Some(path.display().to_string()),
                "to_sdf",
                "IoError",
                format!("create SDF failed: {error}"),
            )],
        })?;
        let success = blocks.len();
        for block in blocks {
            file.write_all(block.as_bytes())
                .map_err(|error| BatchValidationError {
                    errors: vec![BatchRecordError::new(
                        0,
                        Some(path.display().to_string()),
                        "to_sdf",
                        "IoError",
                        format!("write SDF failed: {error}"),
                    )],
                })?;
        }
        Ok(BatchExportReport {
            total: self.records.len(),
            success,
            failed: report_errors.len(),
            errors: report_errors,
        })
    }

    pub fn write_sdf_files(
        &self,
        out_dir: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        filenames: Option<&[Option<String>]>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        self.write_sdf_files_with_progress(out_dir, format, errors, filenames, None)
    }

    pub fn write_sdf_files_with_progress(
        &self,
        out_dir: &Path,
        format: SdfFormat,
        errors: BatchErrorMode,
        filenames: Option<&[Option<String>]>,
        progress: BatchProgress<'_>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let paths = output_paths(
            out_dir,
            self.records.len(),
            "sdf",
            filenames,
            "to_sdf_files",
        )?;
        fs::create_dir_all(out_dir).map_err(|error| BatchValidationError {
            errors: vec![BatchRecordError::new(
                0,
                Some(out_dir.display().to_string()),
                "to_sdf_files",
                "IoError",
                format!("create output directory failed: {error}"),
            )],
        })?;
        let outcomes: Vec<Result<(), BatchRecordError>> = self
            .records
            .par_iter()
            .enumerate()
            .filter_map(|(index, record)| {
                let out = match record {
                    BatchRecord::Valid(molecule) => Some(
                        write_one_sdf_file(molecule, &paths[index], format).map_err(|message| {
                            BatchRecordError::new(
                                index,
                                Some(paths[index].display().to_string()),
                                "to_sdf_files",
                                "SdfWriteError",
                                message,
                            )
                        }),
                    ),
                    BatchRecord::Invalid(error) => match errors {
                        BatchErrorMode::Raise | BatchErrorMode::Keep => Some(Err(error.clone())),
                        BatchErrorMode::Skip => None,
                    },
                };
                tick_progress(progress);
                out
            })
            .collect();
        let mut success = 0usize;
        let mut report_errors = Vec::new();
        for outcome in outcomes {
            match outcome {
                Ok(()) => success += 1,
                Err(error) => report_errors.push(error),
            }
        }
        if matches!(errors, BatchErrorMode::Raise) && !report_errors.is_empty() {
            return Err(BatchValidationError {
                errors: report_errors,
            });
        }
        Ok(BatchExportReport {
            total: self.records.len(),
            success,
            failed: report_errors.len(),
            errors: report_errors,
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
    if molecule.coords_2d().is_some() {
        molblock::mol_to_2d_sdf_record(molecule, format).map_err(|error| error.to_string())
    } else if molecule.coords_3d().is_some() {
        molblock::mol_to_3d_sdf_record(molecule, format).map_err(|error| error.to_string())
    } else {
        Err(
            "SDF writing requires coordinates; call with_2d_coords() or read a molecule with 2D/3D coordinates before writing SDF"
                .to_owned(),
        )
    }
}

fn output_paths(
    out_dir: &Path,
    total: usize,
    extension: &str,
    filenames: Option<&[Option<String>]>,
    stage: &'static str,
) -> Result<Vec<PathBuf>, BatchValidationError> {
    if let Some(filenames) = filenames
        && filenames.len() != total
    {
        return Err(BatchValidationError {
            errors: vec![BatchRecordError::new(
                0,
                None,
                stage,
                "FilenameError",
                format!(
                    "filenames length {} must match batch length {total}",
                    filenames.len()
                ),
            )],
        });
    }

    let mut seen = HashSet::new();
    let mut paths = Vec::with_capacity(total);
    for index in 0..total {
        let filename = match filenames.and_then(|items| items[index].as_deref()) {
            Some(raw) => normalize_output_filename(raw, extension).map_err(|message| {
                BatchValidationError {
                    errors: vec![BatchRecordError::new(
                        index,
                        Some(raw.to_string()),
                        stage,
                        "FilenameError",
                        message,
                    )],
                }
            })?,
            None => format!("{index:06}.{extension}"),
        };
        if !seen.insert(filename.clone()) {
            return Err(BatchValidationError {
                errors: vec![BatchRecordError::new(
                    index,
                    Some(filename),
                    stage,
                    "FilenameError",
                    "duplicate output filename",
                )],
            });
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
    let Some(file_name) = path.file_name().and_then(|value| value.to_str()) else {
        return Err("filename must be valid UTF-8".to_string());
    };
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
    let mut after_record = false;
    let lines = sdf_text.lines().collect::<Vec<_>>();

    for (line_idx, line) in lines.iter().enumerate() {
        if current.is_empty()
            && after_record
            && line.trim().is_empty()
            && lines
                .get(line_idx + 1)
                .is_some_and(|next| next.trim().is_empty())
        {
            after_record = false;
            continue;
        }
        after_record = false;
        current.push_str(line);
        current.push('\n');
        if line.trim_end() == "$$$$" {
            records.push(std::mem::take(&mut current));
            after_record = true;
        }
    }

    if !current.trim().is_empty() {
        records.push(current);
    }

    records
}

#[cfg(test)]
mod tests {
    use super::{BatchErrorMode, BatchRecord, MoleculeBatch};
    use crate::Molecule;
    use crate::io::molblock::SdfFormat;
    use std::fs;
    use std::path::PathBuf;

    fn unique_temp_dir(name: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "cosmolkit-batch-{name}-{}-{:?}",
            std::process::id(),
            std::thread::current().id()
        ))
    }

    #[test]
    fn from_smiles_list_keeps_invalid_records_in_order() {
        let smiles = vec!["CCO".to_string(), "C1CC".to_string(), "CC".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
            .expect("keep mode should not raise");

        assert_eq!(batch.len(), 3);
        assert_eq!(batch.valid_mask(), vec![true, false, true]);
        assert_eq!(batch.errors()[0].index, 1);
        assert_eq!(
            batch
                .filter_valid()
                .to_smiles_list(true)
                .expect("valid records should serialize to SMILES"),
            vec![Some("CCO".to_string()), Some("CC".to_string())]
        );
    }

    #[test]
    fn from_smiles_list_raise_aggregates_errors() {
        let smiles = vec!["C1CC".to_string(), "N1".to_string()];
        let error = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Raise)
            .expect_err("raise mode should aggregate invalid inputs");

        assert_eq!(error.errors.len(), 2);
        assert_eq!(error.errors[0].stage, "parse_smiles");
        assert_eq!(error.errors[1].index, 1);
    }

    #[test]
    fn transforms_skip_invalid_records_when_requested() {
        let smiles = vec!["CCO".to_string(), "C1CC".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
            .expect("keep mode should not raise");
        let prepared = batch
            .add_hydrogens(BatchErrorMode::Skip)
            .expect("skip mode should drop invalid records")
            .compute_2d_coords(BatchErrorMode::Skip)
            .expect("2D coords should compute for valid record");

        assert_eq!(prepared.len(), 1);
        assert_eq!(prepared.valid_count(), 1);
        assert_eq!(prepared.invalid_count(), 0);
    }

    #[test]
    fn sanitize_transform_applies_pipeline_to_valid_records() {
        let raw = Molecule::from_smiles_with_sanitize("CN(=O)=O", false)
            .expect("unsanitized nitro SMILES should parse");
        let batch = MoleculeBatch::new(vec![BatchRecord::Valid(raw)]);

        let sanitized = batch
            .sanitize(BatchErrorMode::Raise)
            .expect("sanitize should transform valid records");

        let BatchRecord::Valid(mol) = &sanitized.records[0] else {
            panic!("record should remain valid");
        };
        assert_eq!(
            mol.atoms()
                .iter()
                .map(|atom| atom.formal_charge)
                .collect::<Vec<_>>(),
            vec![0, 1, -1, 0]
        );
    }

    #[test]
    fn write_images_uses_custom_filenames_and_rejects_unsafe_names() {
        let smiles = vec!["CCO".to_string(), "C1CC".to_string(), "CC".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Keep)
            .expect("keep mode should not raise");
        let out_dir = unique_temp_dir("images");
        let filenames = vec![
            Some("ethanol".to_string()),
            Some("invalid.svg".to_string()),
            None,
        ];

        let report = batch
            .write_images(
                &out_dir,
                "svg",
                120,
                100,
                BatchErrorMode::Skip,
                Some(&filenames),
            )
            .expect("custom image filenames should write");

        assert_eq!(report.success, 2);
        assert!(out_dir.join("ethanol.svg").exists());
        assert!(out_dir.join("000002.svg").exists());
        let _ = fs::remove_dir_all(&out_dir);

        let bad = vec![Some("../escape".to_string()), None, None];
        let error = batch
            .write_images(
                &unique_temp_dir("bad-images"),
                "svg",
                120,
                100,
                BatchErrorMode::Skip,
                Some(&bad),
            )
            .expect_err("unsafe filename should be rejected");
        assert_eq!(error.errors[0].error_type, "FilenameError");
    }

    #[test]
    fn write_sdf_files_uses_custom_filenames() {
        let smiles = vec!["CCO".to_string(), "CC".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles, BatchErrorMode::Raise)
            .expect("SMILES should parse")
            .compute_2d_coords(BatchErrorMode::Raise)
            .expect("2D coords should compute");
        let out_dir = unique_temp_dir("sdf-files");
        let filenames = vec![Some("ethanol".to_string()), Some("ethane.sdf".to_string())];

        let report = batch
            .write_sdf_files(
                &out_dir,
                SdfFormat::V2000,
                BatchErrorMode::Raise,
                Some(&filenames),
            )
            .expect("custom SDF filenames should write");

        assert_eq!(report.success, 2);
        assert!(out_dir.join("ethanol.sdf").exists());
        assert!(out_dir.join("ethane.sdf").exists());
        let _ = fs::remove_dir_all(&out_dir);
    }
}
