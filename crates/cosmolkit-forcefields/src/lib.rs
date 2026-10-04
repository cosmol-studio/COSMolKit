//! Detached UFF/MMFF force-field boundaries.

mod geometry;
mod kernel;
mod mmff;
mod optimizer;
mod uff;

use cosmolkit_model::{CoordinateBlock, TopologyBlock};

pub use uff::{
    UffConformerError, UffConformerOptions, UffConformerOutcome, UffSingleError, UffSingleOptions,
    UffSingleOutcome, optimize_uff_conformers_prepared, optimize_uff_single_prepared,
};
pub use uff::{UffParameterError, UffParameterErrorKind, uff_has_all_molecule_params};

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

pub fn mmff_optimize(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    options: &ForceFieldOptions,
) -> Result<CoordinateBlock, ForceFieldError> {
    let _ = (topology, coordinates, options);
    Err(ForceFieldError::Unsupported)
}
