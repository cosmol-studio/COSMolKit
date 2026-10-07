//! Detached batch records and processing boundaries.

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};
mod sdf;
mod sdf_export;
pub use sdf::{
    BatchErrorMode, BatchProgressBar, BatchReadParams, BatchRecordError, BatchValidationError,
    finalize_sdf_record, read_sdf_dataset, read_sdf_reader, read_sdf_text, validate_record_errors,
};
pub use sdf_export::{
    BatchExportParams, BatchExportReport, SdfExportRecord, export_sdf, export_sdf_files,
    output_paths_with_extension, write_export_report,
};

#[derive(Debug, Clone, PartialEq)]
pub struct BatchRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    /// Moved final assignments from the unique IO finalizer, without live commit authority.
    pub post_state: cosmolkit_io::MolPostDerivedState,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchError {
    Unsupported,
}

pub fn process(records: Vec<BatchRecord>) -> Result<Vec<BatchRecord>, BatchError> {
    let _ = records;
    Err(BatchError::Unsupported)
}

mod scheduler;
pub use scheduler::run_indexed;
