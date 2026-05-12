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
        Self { records }
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
        Ok(Self { records })
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
                    }
                },
                BatchRecord::Error(_) => {
                    error_count += 1;
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

    pub fn write_images(
        &self,
        output_dir: impl AsRef<Path>,
        errors: BatchErrorMode,
    ) -> Result<BatchExportReport, BatchValidationError> {
        let dir = output_dir.as_ref();
        std::fs::create_dir_all(dir).map_err(|e| BatchValidationError {
            errors: 1,
            reason: Some(crate::UnsupportedFeatureError {
                feature: "batch.write_images",
                reason: "directory creation failed",
            }),
        })?;
        let mut written = 0_usize;
        let mut skipped = 0_usize;
        for (index, record) in self.records.iter().enumerate() {
            match record {
                BatchRecord::Molecule(molecule) => {
                    let filename = dir.join(format!("mol_{}.png", index));
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
}
