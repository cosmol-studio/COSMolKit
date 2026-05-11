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
        let unsupported = crate::UnsupportedFeatureError::from_spec(&crate::BATCH_FEATURE);
        Self {
            records: smiles
                .iter()
                .enumerate()
                .map(|(index, _)| {
                    BatchRecord::Error(BatchRecordError::new(
                        index,
                        "batch.from_smiles_list",
                        unsupported.to_string(),
                    ))
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

    pub fn to_smiles_list(
        &self,
        _errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        Err(BatchValidationError {
            errors: self.len(),
            reason: Some(crate::UnsupportedFeatureError::from_spec(
                &crate::BATCH_FEATURE,
            )),
        })
    }

    pub fn to_smiles_list_with_params(
        &self,
        _params: &SmilesWriteParams,
        _errors: BatchErrorMode,
    ) -> Result<Vec<String>, BatchValidationError> {
        Err(BatchValidationError {
            errors: self.len(),
            reason: Some(crate::UnsupportedFeatureError::from_spec(
                &crate::BATCH_FEATURE,
            )),
        })
    }

    pub fn write_images(
        &self,
        _output_dir: impl AsRef<Path>,
        _errors: BatchErrorMode,
    ) -> Result<BatchExportReport, BatchValidationError> {
        Err(BatchValidationError {
            errors: self.len(),
            reason: Some(crate::UnsupportedFeatureError::from_spec(
                &crate::BATCH_FEATURE,
            )),
        })
    }
}
