//! Detached stereochemistry and stereoisomer boundaries.

mod bond_dirs;

pub use bond_dirs::assign_chiral_types_from_bond_dirs;
use cosmolkit_core::ValenceError;
use cosmolkit_model::{Conformer3D, TopologyBlock};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum StereoError {
    #[error("unsupported stereochemistry branch: {reason}")]
    Unsupported { reason: &'static str },
    #[error("invalid stereochemistry state: {0}")]
    InvalidState(String),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct StereoOptions {
    pub max_isomers: Option<usize>,
}

pub fn perceive(topology: &TopologyBlock) -> Result<TopologyBlock, StereoError> {
    let _ = topology;
    Err(StereoError::Unsupported {
        reason: "high-level stereochemistry perception is not detached yet",
    })
}

pub fn enumerate(
    topology: &TopologyBlock,
    conformer: Option<&Conformer3D>,
    options: &StereoOptions,
) -> Result<Vec<TopologyBlock>, StereoError> {
    let _ = (topology, conformer, options);
    Err(StereoError::Unsupported {
        reason: "stereoisomer enumeration is not detached yet",
    })
}
