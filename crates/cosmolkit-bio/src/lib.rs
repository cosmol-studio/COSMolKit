//! Detached structural-biology value and algorithm boundaries.

use cosmolkit_model::TopologyBlock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BioError {
    Unsupported,
}

pub fn select_residues(topology: &TopologyBlock, query: &str) -> Result<TopologyBlock, BioError> {
    let _ = (topology, query);
    Err(BioError::Unsupported)
}
