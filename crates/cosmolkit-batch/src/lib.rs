//! Detached batch records and processing boundaries.

use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

#[derive(Debug, Clone, PartialEq)]
pub struct BatchRecord {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatchError {
    Unsupported,
}

pub fn process(records: Vec<BatchRecord>) -> Result<Vec<BatchRecord>, BatchError> {
    let _ = records;
    Err(BatchError::Unsupported)
}
