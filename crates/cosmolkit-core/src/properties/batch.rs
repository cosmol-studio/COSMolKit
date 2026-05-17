use std::path::Path;

use crate::{Molecule, SmilesWriteParams};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchErrorMode {
    Strict,
    KeepErrors,
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
    n_jobs: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct BatchExportReport {
    pub written: usize,
    pub skipped: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("batch validation failed with {errors} errors")]
pub struct BatchValidationError {
    pub errors: usize,
    pub reason: Option<crate::UnsupportedFeatureError>,
}

pub type BatchProgress<'a> = Option<&'a (dyn Fn() + Sync)>;

pub struct BatchProgressBar;

impl BatchProgressBar {
    #[must_use]
    pub fn new(_total: usize, _message: impl Into<String>) -> Self {
        Self
    }

    pub fn callback(&self) -> Box<dyn Fn() + Sync + '_> {
        Box::new(|| {})
    }

    pub fn inc(&self, _delta: u64) {}

    pub fn finish(&self) {}
}

pub fn batch_progress_bar(total: usize, message: impl Into<String>) -> BatchProgressBar {
    BatchProgressBar::new(total, message)
}

impl MoleculeBatch {
    #[must_use]
    pub fn new(records: Vec<BatchRecord>) -> Self {
        Self { records, n_jobs: 1 }
    }

    /// Set the number of parallel jobs for batch processing.
    ///
    /// This is a design placeholder for future parallel execution. The actual
    /// threading/rayon integration is deferred; this sets the target degree of
    /// parallelism that transforms and exports will use when threading support
    /// is wired in.
    #[must_use]
    pub fn with_parallel_jobs(mut self, n_jobs: usize) -> Self {
        self.n_jobs = n_jobs.max(1);
        self
    }

    /// Return the configured number of parallel jobs.
    #[must_use]
    pub fn parallel_jobs(&self) -> usize {
        self.n_jobs
    }

    /// Return a reference to the record at the given index, or `None` if out of bounds.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<&BatchRecord> {
        self.records.get(index)
    }

    /// Iterate over all records in order.
    pub fn iter(&self) -> impl Iterator<Item = &BatchRecord> {
        self.records.iter()
    }

    /// Return a human-readable summary of all errors in the batch.
    ///
    /// Returns an empty string when there are no errors.
    #[must_use]
    pub fn error_summary(&self) -> String {
        let errors: Vec<_> = self
            .records
            .iter()
            .filter_map(|record| match record {
                BatchRecord::Error(error) => Some(error.clone()),
                BatchRecord::Molecule(_) => None,
            })
            .collect();
        if errors.is_empty() {
            return String::new();
        }
        let total = self.records.len();
        let valid = total - errors.len();
        let mut lines = Vec::with_capacity(errors.len() + 2);
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

    /// Return a structured error report.
    ///
    /// Each entry pairs the index with the error details: `(index, operation, message)`.
    #[must_use]
    pub fn error_report(&self) -> Vec<(usize, &'static str, String)> {
        self.records
            .iter()
            .filter_map(|record| match record {
                BatchRecord::Error(error) => {
                    Some((error.index, error.operation, error.message.clone()))
                }
                BatchRecord::Molecule(_) => None,
            })
            .collect()
    }

    pub fn from_smiles_list(smiles: &[String]) -> Self {
        // BEGIN COSMolKit SOURCE-PORT FRAME README checked batch.from_smiles_list
        // preserve input order
        // keep per-record errors as typed batch records
        // do not add batch-specific chemistry; delegate parsing to Molecule::from_smiles
        // END COSMolKit SOURCE-PORT FRAME README checked batch.from_smiles_list
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
            n_jobs: 1,
        }
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
                BatchRecord::Error(error) => Some(error.clone()),
                BatchRecord::Molecule(_) => None,
            })
            .collect()
    }

    #[must_use]
    pub fn valid_count(&self) -> usize {
        self.valid_mask().into_iter().filter(|valid| *valid).count()
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
        }
    }

    pub fn sanitized(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.transform("batch.sanitized", errors, Molecule::sanitized)
    }

    pub fn with_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.transform("batch.with_hydrogens", errors, Molecule::with_hydrogens)
    }

    pub fn without_hydrogens(&self, errors: BatchErrorMode) -> Result<Self, BatchValidationError> {
        self.transform(
            "batch.without_hydrogens",
            errors,
            Molecule::without_hydrogens,
        )
    }

    pub fn with_kekulized_bonds(
        &self,
        clear_aromatic_flags: bool,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.transform("batch.with_kekulized_bonds", errors, |molecule| {
            molecule.with_kekulized_bonds(clear_aromatic_flags)
        })
    }

    pub fn with_2d_coordinates(
        &self,
        errors: BatchErrorMode,
    ) -> Result<Self, BatchValidationError> {
        self.transform(
            "batch.with_2d_coordinates",
            errors,
            Molecule::with_2d_coordinates,
        )
    }

    fn transform(
        &self,
        operation: &'static str,
        errors: BatchErrorMode,
        mut transform: impl FnMut(&Molecule) -> Result<Molecule, crate::OperationError>,
    ) -> Result<Self, BatchValidationError> {
        // BEGIN COSMolKit SOURCE-PORT FRAME README checked batch transformations
        // preserve input order
        // keep batch as orchestration over registered molecule operations
        // never mutate source batch records in place
        // keep or raise per-record errors according to BatchErrorMode
        // END COSMolKit SOURCE-PORT FRAME README checked batch transformations
        let mut records = Vec::with_capacity(self.records.len());
        let mut error_count = 0_usize;
        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => match transform(molecule) {
                    Ok(molecule) => records.push(BatchRecord::Molecule(molecule)),
                    Err(error) => {
                        error_count += 1;
                        records.push(BatchRecord::Error(BatchRecordError::new(
                            index,
                            operation,
                            error.to_string(),
                        )));
                    }
                },
                BatchRecord::Error(error) => {
                    error_count += 1;
                    records.push(BatchRecord::Error(error.clone()));
                }
            }
        }
        if errors == BatchErrorMode::Strict && error_count != 0 {
            return Err(BatchValidationError {
                errors: error_count,
                reason: None,
            });
        }
        Ok(Self {
            records,
            n_jobs: self.n_jobs,
        })
    }

    pub fn to_smiles_list(
        &self,
        errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        let mut results = Vec::with_capacity(self.records.len());
        let mut error_count = 0_usize;
        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => {
                    results.push(molecule.to_smiles(true).unwrap_or_else(|_| "?".to_string()));
                }
                BatchRecord::Error(_) => {
                    error_count += 1;
                    results.push("?".to_string());
                }
            }
        }
        if errors == BatchErrorMode::Strict && error_count != 0 {
            return Err(BatchValidationError {
                errors: error_count,
                reason: None,
            });
        }
        Ok(results)
    }

    pub fn to_smiles_list_with_params(
        &self,
        params: &SmilesWriteParams,
        errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        let mut results = Vec::with_capacity(self.records.len());
        let mut error_count = 0_usize;
        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => match molecule.to_smiles_with_params(params) {
                    Ok(s) => results.push(s),
                    Err(e) => {
                        error_count += 1;
                        if errors == BatchErrorMode::Strict {
                            return Err(BatchValidationError {
                                errors: error_count,
                                reason: None,
                            });
                        }
                        results.push("?".to_string());
                    }
                },
                BatchRecord::Error(_) => {
                    error_count += 1;
                    if errors == BatchErrorMode::Strict {
                        return Err(BatchValidationError {
                            errors: error_count,
                            reason: None,
                        });
                    }
                    results.push("?".to_string());
                }
            }
        }
        if errors == BatchErrorMode::Strict && error_count != 0 {
            return Err(BatchValidationError {
                errors: error_count,
                reason: None,
            });
        }
        Ok(results)
    }

    /// Write molecular images (PNG) for valid records.
    ///
    /// `filenames` is an optional list of filenames (without directory prefix).
    /// When provided, its length must equal the batch length; `BatchRecord::Error`
    /// entries still consume a slot in the filenames list (the corresponding file
    /// is skipped). When `None`, filenames are generated as `mol_{index}.png`.
    pub fn write_images(
        &self,
        output_dir: impl AsRef<Path>,
        errors: BatchErrorMode,
        filenames: Option<&[String]>,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let dir = output_dir.as_ref();
        std::fs::create_dir_all(dir).map_err(|e| BatchValidationError {
            errors: 1,
            reason: Some(crate::UnsupportedFeatureError {
                feature: "batch.write_images",
                reason: "directory creation failed",
            }),
        })?;

        if let Some(names) = filenames {
            if names.len() != self.records.len() {
                return Err(BatchValidationError {
                    errors: 1,
                    reason: Some(crate::UnsupportedFeatureError {
                        feature: "batch.write_images",
                        reason: "filenames length must match batch length",
                    }),
                });
            }
        }

        let mut written = 0_usize;
        let mut skipped = 0_usize;
        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => {
                    let filename = match filenames {
                        Some(names) => dir.join(&names[index]),
                        None => dir.join(format!("mol_{}.png", index)),
                    };
                    match crate::draw::mol_to_png(molecule, 500, 500) {
                        Ok(png_bytes) => {
                            if std::fs::write(&filename, &png_bytes).is_ok() {
                                written += 1;
                            } else {
                                skipped += 1;
                            }
                        }
                        Err(_) => {
                            skipped += 1;
                        }
                    }
                }
                BatchRecord::Error(_) => {
                    skipped += 1;
                }
            }
        }
        if errors == BatchErrorMode::Strict && skipped != 0 {
            return Err(BatchValidationError {
                errors: skipped,
                reason: None,
            });
        }
        Ok(BatchExportReport { written, skipped })
    }

    /// Write all valid molecules as a multi-record SDF file.
    ///
    /// Each valid `BatchRecord::Molecule` is written as a separate SDF record
    /// using the V2000 format with 2D coordinates. Error records are skipped.
    pub fn write_sdf(
        &self,
        output_path: impl AsRef<Path>,
        errors: BatchErrorMode,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let path = output_path.as_ref();

        let mut sdf_content = String::new();
        let mut written = 0_usize;
        let mut skipped = 0_usize;

        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => {
                    match crate::io::molblock::mol_to_sdf_record_with_params(
                        molecule,
                        &crate::io::molblock::MolBlockWriteParams::default(),
                    ) {
                        Ok(record_str) => {
                            sdf_content.push_str(&record_str);
                            written += 1;
                        }
                        Err(_) => {
                            skipped += 1;
                        }
                    }
                }
                BatchRecord::Error(_) => {
                    skipped += 1;
                }
            }
        }

        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent).map_err(|e| BatchValidationError {
                errors: 1,
                reason: Some(crate::UnsupportedFeatureError {
                    feature: "batch.write_sdf",
                    reason: "directory creation failed",
                }),
            })?;
        }

        std::fs::write(path, &sdf_content).map_err(|e| BatchValidationError {
            errors: 1,
            reason: Some(crate::UnsupportedFeatureError {
                feature: "batch.write_sdf",
                reason: &"file write failed",
            }),
        })?;

        if errors == BatchErrorMode::Strict && skipped != 0 {
            return Err(BatchValidationError {
                errors: skipped,
                reason: None,
            });
        }
        Ok(BatchExportReport { written, skipped })
    }
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
    fn batch_transforms_delegate_to_registered_ops_and_keep_errors_in_order() {
        let smiles = vec!["CCO".to_string(), "C1".to_string(), "N".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles);

        let transformed = batch
            .with_kekulized_bonds(true, BatchErrorMode::KeepErrors)
            .unwrap();

        assert_eq!(transformed.len(), 3);
        assert_eq!(transformed.valid_mask(), vec![true, false, true]);
        assert_eq!(transformed.errors()[0].index, 1);
    }

    #[test]
    fn batch_transform_strict_reports_existing_or_operation_errors() {
        let smiles = vec!["C1CC1".to_string(), "C1".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles);

        let existing_error = batch
            .sanitized(BatchErrorMode::Strict)
            .expect_err("strict mode should reject existing record errors");

        // "C1CC1" (cyclopropane) — coordinate computation now handles rings
        // "C1" was filtered out by the sanitization check above
        let operation_ok = batch
            .filter_valid()
            .with_2d_coordinates(BatchErrorMode::Strict);

        assert_eq!(existing_error.errors, 1);
        assert!(
            operation_ok.is_ok(),
            "coord computation now succeeds for cyclopropane"
        );
    }

    #[test]
    fn batch_transform_keep_errors_records_coordinate_errors() {
        let smiles = vec!["C1CC1".to_string(), "N".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles);

        let transformed = batch
            .with_2d_coordinates(BatchErrorMode::KeepErrors)
            .unwrap();

        // "C1CC1" (cyclopropane) — rings are now handled by the ported SSSR algorithm
        // "N" (1 atom) — strict subset handles single atoms, succeeds
        assert_eq!(transformed.valid_mask(), vec![true, true]);
        assert_eq!(transformed.errors().len(), 0);
    }

    #[test]
    fn with_parallel_jobs_default_is_one() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]);
        assert_eq!(batch.parallel_jobs(), 1);
    }

    #[test]
    fn with_parallel_jobs_sets_and_gets_value() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]).with_parallel_jobs(4);
        assert_eq!(batch.parallel_jobs(), 4);
    }

    #[test]
    fn with_parallel_jobs_clamps_to_at_least_one() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]).with_parallel_jobs(0);
        assert_eq!(batch.parallel_jobs(), 1);
    }

    #[test]
    fn with_parallel_jobs_preserved_across_filter_valid() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "C1".to_string()])
            .with_parallel_jobs(8);
        let filtered = batch.filter_valid();
        assert_eq!(filtered.parallel_jobs(), 8);
    }

    #[test]
    fn get_returns_record_at_index() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let rec0 = batch.get(0);
        assert!(rec0.is_some());
        assert!(matches!(rec0.unwrap(), BatchRecord::Molecule(_)));
        let rec1 = batch.get(1);
        assert!(rec1.is_some());
        assert!(matches!(rec1.unwrap(), BatchRecord::Molecule(_)));
    }

    #[test]
    fn get_returns_none_for_out_of_bounds() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]);
        assert!(batch.get(5).is_none());
        assert!(batch.get(usize::MAX).is_none());
    }

    #[test]
    fn iter_yields_all_records_in_order() {
        let smiles = vec!["CCO".to_string(), "C1".to_string(), "N".to_string()];
        let batch = MoleculeBatch::from_smiles_list(&smiles);
        let records: Vec<&BatchRecord> = batch.iter().collect();
        assert_eq!(records.len(), 3);
        assert!(matches!(records[0], BatchRecord::Molecule(_)));
        assert!(matches!(records[1], BatchRecord::Error(_)));
        assert!(matches!(records[2], BatchRecord::Molecule(_)));
    }

    #[test]
    fn error_summary_empty_when_no_errors() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let summary = batch.error_summary();
        assert!(summary.is_empty());
    }

    #[test]
    fn error_summary_includes_counts_and_details() {
        let batch = MoleculeBatch::from_smiles_list(&[
            "CCO".to_string(),
            "C1".to_string(),
            "N".to_string(),
        ]);
        let summary = batch.error_summary();
        assert!(summary.contains("1 errors out of 3 records"));
        assert!(summary.contains("2 valid"));
        assert!(summary.contains("[1]"));
        assert!(summary.contains("batch.from_smiles_list"));
    }

    #[test]
    fn error_report_returns_structured_tuples() {
        let batch = MoleculeBatch::from_smiles_list(&[
            "CCO".to_string(),
            "C1".to_string(),
            "N".to_string(),
        ]);
        let report = batch.error_report();
        assert_eq!(report.len(), 1);
        assert_eq!(report[0].0, 1);
        assert_eq!(report[0].1, "batch.from_smiles_list");
    }

    #[test]
    fn write_images_custom_filenames() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let filenames = vec!["first.png".to_string(), "second.png".to_string()];
        let report = batch
            .write_images(dir.path(), BatchErrorMode::KeepErrors, Some(&filenames))
            .unwrap();
        assert_eq!(report.written, 2);
        assert!(dir.path().join("first.png").exists());
        assert!(dir.path().join("second.png").exists());
    }

    #[test]
    fn write_images_custom_filenames_error_records_slot() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "C1".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let filenames = vec!["ok.png".to_string(), "bad.png".to_string()];
        let report = batch
            .write_images(dir.path(), BatchErrorMode::KeepErrors, Some(&filenames))
            .unwrap();
        // Only one valid molecule, but bad.png slot also reserved
        assert_eq!(report.written, 1);
        assert_eq!(report.skipped, 1);
        assert!(dir.path().join("ok.png").exists());
        assert!(!dir.path().join("bad.png").exists());
    }

    #[test]
    fn write_images_custom_filenames_requires_correct_length() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let filenames = vec!["only.png".to_string()];
        let result = batch.write_images(dir.path(), BatchErrorMode::KeepErrors, Some(&filenames));
        assert!(result.is_err());
    }

    #[test]
    fn write_images_no_filenames_uses_default_index_names() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let report = batch
            .write_images(dir.path(), BatchErrorMode::KeepErrors, None)
            .unwrap();
        assert_eq!(report.written, 2);
        assert!(dir.path().join("mol_0.png").exists());
        assert!(dir.path().join("mol_1.png").exists());
    }

    #[test]
    fn write_sdf_writes_multi_record_sdf() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "N".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("output.sdf");
        let report = batch.write_sdf(&path, BatchErrorMode::KeepErrors).unwrap();
        assert_eq!(report.written, 2);
        assert_eq!(report.skipped, 0);
        let content = std::fs::read_to_string(&path).unwrap();
        // Each SDF record ends with $$$$
        assert_eq!(content.matches("$$$$").count(), 2);
        // Both molecules appear (count V2000 headers)
        assert_eq!(content.matches("V2000").count(), 2);
    }

    #[test]
    fn write_sdf_skips_error_records() {
        let batch = MoleculeBatch::from_smiles_list(&[
            "CCO".to_string(),
            "C1".to_string(),
            "N".to_string(),
        ]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("output.sdf");
        let report = batch.write_sdf(&path, BatchErrorMode::KeepErrors).unwrap();
        // C1 is invalid, so only 2 written
        assert_eq!(report.written, 2);
        assert_eq!(report.skipped, 1);
        let content = std::fs::read_to_string(&path).unwrap();
        assert_eq!(content.matches("$$$$").count(), 2);
    }

    #[test]
    fn write_sdf_strict_reports_skipped_records() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "C1".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("output.sdf");
        let result = batch.write_sdf(&path, BatchErrorMode::Strict);
        assert!(result.is_err());
        assert_eq!(result.unwrap_err().errors, 1);
    }

    #[test]
    fn write_sdf_empty_batch() {
        let batch = MoleculeBatch::new(vec![]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("output.sdf");
        let report = batch.write_sdf(&path, BatchErrorMode::KeepErrors).unwrap();
        assert_eq!(report.written, 0);
        assert_eq!(report.skipped, 0);
        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.is_empty());
    }

    #[test]
    fn error_summary_escape_does_not_panic() {
        let batch = MoleculeBatch::new(vec![BatchRecord::Error(BatchRecordError::new(
            42,
            "test.op",
            "something went wrong: Ω≈ç√∫",
        ))]);
        let summary = batch.error_summary();
        assert!(summary.contains("42"));
        assert!(summary.contains("test.op"));
    }

    #[test]
    fn write_sdf_creates_parent_directory() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string()]);
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("subdir").join("nested").join("output.sdf");
        let report = batch.write_sdf(&path, BatchErrorMode::KeepErrors).unwrap();
        assert_eq!(report.written, 1);
        assert!(path.exists());
    }

    #[test]
    fn to_smiles_list_with_params_keep_errors_preserves_record_slots() {
        let batch = MoleculeBatch::from_smiles_list(&[
            "CCO".to_string(),
            "C1".to_string(),
            "N".to_string(),
        ]);

        let smiles = batch
            .to_smiles_list_with_params(&SmilesWriteParams::default(), BatchErrorMode::KeepErrors)
            .unwrap();

        assert_eq!(smiles.len(), 3);
        assert_eq!(smiles[1], "?");
    }

    #[test]
    fn to_smiles_list_with_params_strict_rejects_existing_error_slots() {
        let batch = MoleculeBatch::from_smiles_list(&["CCO".to_string(), "C1".to_string()]);

        let error = batch
            .to_smiles_list_with_params(&SmilesWriteParams::default(), BatchErrorMode::Strict)
            .unwrap_err();

        assert_eq!(error.errors, 1);
    }
}
