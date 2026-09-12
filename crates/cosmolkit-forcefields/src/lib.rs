//! Detached UFF/MMFF force-field boundaries.

use cosmolkit_model::{CoordinateBlock, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForceFieldError {
    Unsupported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ForceFieldOptions {
    pub max_iterations: usize,
}

pub fn mmff_has_all_molecule_params(topology: &TopologyBlock) -> Result<bool, ForceFieldError> {
    let _ = topology;
    Err(ForceFieldError::Unsupported)
}

pub fn uff_has_all_molecule_params(topology: &TopologyBlock) -> Result<bool, ForceFieldError> {
    let _ = topology;
    Err(ForceFieldError::Unsupported)
}

pub fn mmff_optimize(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    options: &ForceFieldOptions,
) -> Result<CoordinateBlock, ForceFieldError> {
    let _ = (topology, coordinates, options);
    Err(ForceFieldError::Unsupported)
}
