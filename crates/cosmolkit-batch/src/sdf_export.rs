//! Ordered source SDF export over borrowed detached values.
use crate::{BatchErrorMode, BatchProgressBar, BatchRecordError, BatchValidationError};
use cosmolkit_io::{MolBlockWriteParams, SdfFormat};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::Write;
use std::path::{Component, Path, PathBuf};

pub enum SdfExportRecord<'a> {
    Molecule(cosmolkit_io::MolWriteInput<'a>),
    Error(&'a BatchRecordError),
}
#[derive(Clone, Debug)]
pub struct BatchExportReport {
    pub written: usize,
    /// All unsuccessful records, including errors already present in the input.
    pub failed: usize,
    pub errors: Vec<BatchRecordError>,
}
impl BatchExportReport {
    /// Write the source-defined count-only JSON/CSV report to a literal path.
    /// Path expansion is the caller's language-boundary responsibility.
    pub fn write_report(&self, path: &Path) -> Result<(), BatchValidationError> {
        // COSMolKit❗✔️: pinned d892ec3507c5b568c5ed5d86ae44e466f7d03855
        // python/src/lib.rs:4882, report writing after image export returns:
        //             write_batch_report(path, &report)?;
        // Reuse the sole report serializer/writer without cloning or buffering.
        write_export_report(path, self)
    }
    pub fn total(&self) -> usize {
        self.written + self.failed
    }
    pub fn success(&self) -> usize {
        self.written
    }
    pub fn failed(&self) -> usize {
        self.failed
    }
    pub fn errors(&self) -> &[BatchRecordError] {
        &self.errors
    }
}
#[derive(Clone, Copy, Debug)]
pub struct BatchExportParams {
    pub format: SdfFormat,
    /// The live batch resolves inheritance; detached exports default to Strict.
    pub errors: Option<BatchErrorMode>,
    pub n_jobs: Option<usize>,
    pub progress_bar: Option<bool>,
}
impl Default for BatchExportParams {
    fn default() -> Self {
        Self {
            format: SdfFormat::V2000,
            errors: None,
            n_jobs: None,
            progress_bar: None,
        }
    }
}
fn failure(
    operation: &'static str,
    error: impl std::error::Error + Send + Sync + 'static,
) -> BatchValidationError {
    BatchValidationError::from_record_errors(vec![BatchRecordError::with_source(
        0, operation, error,
    )])
}
fn serialize(
    record: &SdfExportRecord<'_>,
    format: SdfFormat,
    index: usize,
    operation: &'static str,
) -> Result<cosmolkit_model::PropertyText, BatchRecordError> {
    match record {
        SdfExportRecord::Error(error) => Err((*error).clone()),
        SdfExportRecord::Molecule(data) => {
            let params = MolBlockWriteParams {
                format,
                force_2d: !data.coordinates.conformers_2d.is_empty(),
                ..Default::default()
            };
            cosmolkit_io::write_sdf_with_params(*data, &params)
                .map_err(|e| BatchRecordError::with_source(index, operation, e))
        }
    }
}
pub fn export_sdf(
    records: &[SdfExportRecord<'_>],
    path: &Path,
    params: &BatchExportParams,
) -> Result<BatchExportReport, BatchValidationError> {
    // d892 properties/batch.rs::write_sdf_with_runtime: prepare every indexed
    // result, validate strict errors before opening, then append valid blocks
    // in input order. The detached owner never receives live runtime values.
    let progress = params
        .progress_bar
        .unwrap_or(false)
        .then(|| BatchProgressBar::new(records.len(), "Writing SDF records"));
    let result = crate::sdf::parallel(params.n_jobs, || {
        records
            .par_iter()
            .enumerate()
            .map(|(index, record)| {
                let out = serialize(record, params.format, index, "batch.write_sdf");
                if let Some(progress) = &progress {
                    progress.inc(1);
                }
                out
            })
            .collect::<Vec<_>>()
    });
    if let Some(progress) = progress {
        progress.finish();
    }
    let mut blocks = Vec::new();
    let mut errors = Vec::new();
    for result in result? {
        match result {
            Ok(block) => blocks.push(block),
            Err(error) => errors.push(error),
        }
    }
    // User-approved export contract: keep original input errors alongside new
    // serialization errors; KEEP continues other rows without erasing failures.
    if params.errors.unwrap_or_default() == BatchErrorMode::Strict && !errors.is_empty() {
        return Err(BatchValidationError::from_record_errors(errors));
    }
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| failure("batch.write_sdf", e))?;
    }
    let mut file = std::fs::File::create(path).map_err(|e| failure("batch.write_sdf", e))?;
    for block in &blocks {
        file.write_all(block.as_bytes())
            .map_err(|e| failure("batch.write_sdf", e))?;
    }
    Ok(BatchExportReport {
        written: blocks.len(),
        failed: errors.len(),
        errors,
    })
}
pub fn export_sdf_files(
    records: &[SdfExportRecord<'_>],
    directory: &Path,
    params: &BatchExportParams,
    filenames: Option<&[Option<String>]>,
) -> Result<BatchExportReport, BatchValidationError> {
    // The old Python projection completed missing filename entries before
    // the owner created a directory. Retain that ordering at this Rust edge.
    if let Some(names) = filenames
        && names.len() != records.len()
    {
        return Err(BatchValidationError::from_record_errors(vec![
            BatchRecordError::new(
                0,
                "filenames",
                format!(
                    "filenames length must match batch length: expected {}, got {}",
                    records.len(),
                    names.len()
                ),
            ),
        ]));
    }
    std::fs::create_dir_all(directory).map_err(|e| failure("batch.write_sdf_files", e))?;
    let paths = output_paths(directory, records.len(), filenames)?;
    let progress = params
        .progress_bar
        .unwrap_or(false)
        .then(|| BatchProgressBar::new(records.len(), "Writing SDF files"));
    let outcomes = crate::sdf::parallel(params.n_jobs, || {
        records
            .par_iter()
            .enumerate()
            .map(|(index, record)| {
                let result = match record {
                    SdfExportRecord::Error(error) => Err((*error).clone()),
                    SdfExportRecord::Molecule(_) => {
                        serialize(record, params.format, index, "batch.write_sdf_files").and_then(
                            |block| {
                                std::fs::write(&paths[index], block).map_err(|e| {
                                    BatchRecordError::with_source(index, "batch.write_sdf_files", e)
                                })
                            },
                        )
                    }
                };
                if let Some(progress) = &progress {
                    progress.inc(1);
                }
                result
            })
            .collect::<Vec<_>>()
    });
    if let Some(progress) = progress {
        progress.finish();
    }
    let mut written = 0;
    let mut errors = Vec::new();
    for result in outcomes? {
        match result {
            Ok(()) => written += 1,
            Err(error) => errors.push(error),
        }
    }
    if params.errors.unwrap_or_default() == BatchErrorMode::Strict && !errors.is_empty() {
        return Err(BatchValidationError::from_record_errors(errors));
    }
    Ok(BatchExportReport {
        written,
        failed: errors.len(),
        errors,
    })
}
fn output_paths(
    out_dir: &Path,
    total: usize,
    filenames: Option<&[Option<String>]>,
) -> Result<Vec<PathBuf>, BatchValidationError> {
    output_paths_with_extension(out_dir, total, "sdf", filenames)
}
pub fn output_paths_with_extension(
    out_dir: &Path,
    total: usize,
    extension: &str,
    filenames: Option<&[Option<String>]>,
) -> Result<Vec<PathBuf>, BatchValidationError> {
    if let Some(names) = filenames
        && names.len() != total
    {
        return Err(BatchValidationError::parameter(
            "batch.output_paths",
            "filenames length must match batch length",
        ));
    }
    // COSMolKit❗✔️: pinned d892ec3 properties/batch.rs::output_paths; same indexed order/all-result collection and scalar owner calls.
    // fn output_paths(
    //     out_dir: &Path,
    //     total: usize,
    //     extension: &str,
    //     filenames: Option<&[String]>,
    // ) -> Result<Vec<PathBuf>, BatchValidationError> {
    //     if let Some(filenames) = filenames
    //         && filenames.len() != total
    //     {
    //         return Err(BatchValidationError::unsupported(
    //             "batch.output_paths",
    //             "filenames length must match batch length",
    //         ));
    //     }
    //     let mut seen = HashSet::new();
    //     let mut paths = Vec::with_capacity(total);
    //     for index in 0..total {
    //         let filename = match filenames.and_then(|names| names.get(index)) {
    //             Some(raw) => normalize_output_filename(raw, extension).map_err(|_| {
    //                 BatchValidationError::unsupported("batch.output_paths", "invalid filename")
    //             })?,
    //             None => format!("mol_{index}.{extension}"),
    //         };
    //         if !seen.insert(filename.clone()) {
    //             return Err(BatchValidationError::unsupported(
    //                 "batch.output_paths",
    //                 "duplicate output filename",
    //             ));
    //         }
    //         paths.push(out_dir.join(filename));
    //     }
    //     Ok(paths)
    // }

    let mut seen = HashSet::new();
    let mut paths = Vec::with_capacity(total);
    for index in 0..total {
        let name = match filenames.and_then(|names| names[index].as_deref()) {
            Some(raw) => normalize_output_filename(raw, extension).map_err(|reason| {
                BatchValidationError::from_record_errors(vec![BatchRecordError::new(
                    index,
                    "batch.output_paths",
                    reason,
                )])
            })?,
            None => format!("mol_{index}.{extension}"),
        };
        if !seen.insert(name.clone()) {
            return Err(BatchValidationError::parameter(
                "batch.output_paths",
                "duplicate output filename",
            ));
        }
        paths.push(out_dir.join(name));
    }
    Ok(paths)
}
fn normalize_output_filename(raw: &str, extension: &str) -> Result<String, String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return Err("filename must not be empty".into());
    }
    let path = Path::new(trimmed);
    if path.is_absolute() {
        return Err("filename must be relative to the output directory".into());
    }
    let components = path.components().collect::<Vec<_>>();
    if components.len() != 1 || !matches!(components[0], Component::Normal(_)) {
        return Err("filename must not include path separators or '..'".into());
    }
    let name = path
        .file_name()
        .and_then(|v| v.to_str())
        .ok_or_else(|| "filename must be valid UTF-8".to_string())?;
    match path.extension().and_then(|v| v.to_str()) {
        Some(actual) if actual.eq_ignore_ascii_case(extension) => Ok(name.to_string()),
        Some(actual) => Err(format!(
            "filename extension '.{actual}' does not match expected '.{extension}'"
        )),
        None => Ok(format!("{name}.{extension}")),
    }
}
pub fn write_export_report(
    path: &Path,
    report: &BatchExportReport,
) -> Result<(), BatchValidationError> {
    // COSMolKit❗✔️: pinned d892ec3 Python write_batch_report:530-553.
    // fn write_batch_report(path: &str, report: &cosmolkit_core::BatchExportReport) -> PyResult<()> {
    //     let expanded_path = expand_user_path(path)?;
    //     let ext = expanded_path
    //         .extension()
    //         .and_then(|s| s.to_str())
    //         .unwrap_or("json")
    //         .to_ascii_lowercase();
    //     let content = if ext == "csv" {
    //         format!(
    //             "written,skipped,failed\n{},{},{}\n",
    //             report.written, report.skipped, report.failed
    //         )
    //     } else {
    //         format!(
    //             "{{\n  \"written\": {},\n  \"skipped\": {},\n  \"failed\": {}\n}}\n",
    //             report.written, report.skipped, report.failed
    //         )
    //     };
    //     fs::write(&expanded_path, content)
    //         .map_err(|err| PyValueError::new_err(format!("write error report failed: {err}")))
    // }
    //
    // fn complete_batch_filenames(
    //     filenames: Option<Vec<Option<String>>>,
    // Language boundary supplies the expanded path; Rust Path remains literal.
    // User-approved report contract removes skipped: every unsuccessful input
    // contributes once to failed. Structured errors remain in the report value.
    let ext = path
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("json")
        .to_ascii_lowercase();
    let content = if ext == "csv" {
        format!("written,failed\n{},{}\n", report.written, report.failed)
    } else {
        format!(
            "{{\n  \"written\": {},\n  \"failed\": {}\n}}\n",
            report.written, report.failed
        )
    };
    std::fs::write(path, content).map_err(|e| failure("write error report", e))
}
