//! Detached tautomer transformation and enumeration boundaries.

use cosmolkit_model::{MoleculeProperties, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TautomerError {
    Unsupported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct TautomerOptions {
    pub max_results: Option<usize>,
}

pub fn enumerate(
    topology: &TopologyBlock,
    properties: &MoleculeProperties,
    options: &TautomerOptions,
) -> Result<Vec<(TopologyBlock, MoleculeProperties)>, TautomerError> {
    let _ = (topology, properties, options);
    Err(TautomerError::Unsupported)
}
