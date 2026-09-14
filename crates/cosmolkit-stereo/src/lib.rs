//! Detached stereochemistry and stereoisomer boundaries.

mod cip_graph;
mod cip_labels;

pub use cip_graph::CipLabelerError;
pub use cip_labels::{CipLabelAssignment, CipLabelOptions, assign_cip_labels};
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

pub fn assign_chiral_types_from_bond_dirs(
    topology: &mut TopologyBlock,
    conformer: &Conformer3D,
    replace_existing_tags: bool,
) -> Result<(), StereoError> {
    cosmolkit_core::assign_chiral_types_from_bond_dirs(topology, conformer, replace_existing_tags)
        .map_err(|error| match error {
            cosmolkit_core::BondDirectionStereoError::InvalidState(message) => {
                StereoError::InvalidState(message)
            }
            cosmolkit_core::BondDirectionStereoError::Valence(error) => StereoError::Valence(error),
        })
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
