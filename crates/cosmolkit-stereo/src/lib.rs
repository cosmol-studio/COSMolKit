//! Detached stereochemistry and stereoisomer boundaries.

mod chiral_centers;
mod cip_graph;
mod cip_labels;
mod tetrahedral;
pub use chiral_centers::find_chiral_centers;
pub use tetrahedral::{StereoReadError, perceive_stereochemistry, tetrahedral_stereo};

pub use cip_graph::CipLabelerError;
#[doc(hidden)]
pub use cip_labels::assign_cip_labels_cow;
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

#[cfg(feature = "enumeration")]
mod stereo_enumerate;
#[cfg(feature = "enumeration")]
pub use stereo_enumerate::{
    EnumerationError, StereoisomerIterator, StereoisomerOptions, StereoisomerRandomSource,
    enumerate_stereoisomers, enumerate_stereoisomers_with_random_bits, stereoisomer_count,
};
