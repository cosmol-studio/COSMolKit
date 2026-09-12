//! Detached conformer-generation boundaries.

use cosmolkit_model::{Conformer3D, CoordinateBlock, TopologyBlock};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConformerError {
    Unsupported,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct ConformerOptions {
    pub max_conformers: usize,
}

pub fn embed(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    options: &ConformerOptions,
) -> Result<Conformer3D, ConformerError> {
    let _ = (topology, coordinates, options);
    Err(ConformerError::Unsupported)
}

pub fn embed_multiple(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
    options: &ConformerOptions,
) -> Result<Vec<Conformer3D>, ConformerError> {
    let _ = (topology, coordinates, options);
    Err(ConformerError::Unsupported)
}
