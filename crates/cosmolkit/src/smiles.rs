//! Thin public SMILES construction over the detached parser and chemistry owners.

use std::fmt;

use crate::{Molecule, OperationError};

/// Structured failure from the public SMILES construction pipeline.
#[derive(Debug)]
pub enum SmilesError {
    /// SMILES grammar, CXSMILES, or detached parser validation failed.
    Parse(cosmolkit_smiles::SmilesParseError),
    /// Source-defined hydrogen removal failed.
    Hydrogen(cosmolkit_core::HydrogenError),
    /// Source-defined sanitization failed.
    Sanitize(cosmolkit_core::SanitizeError),
    /// Source-defined post-chemistry stereo completion failed.
    Stereo(cosmolkit_smiles::SmilesStereoError),
    /// Final live-state validation failed.
    Construction(OperationError),
}

impl fmt::Display for SmilesError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Parse(error) => write!(formatter, "SMILES parsing failed: {error}"),
            Self::Hydrogen(error) => write!(formatter, "SMILES hydrogen removal failed: {error}"),
            Self::Sanitize(error) => write!(formatter, "SMILES sanitization failed: {error}"),
            Self::Stereo(error) => write!(formatter, "SMILES stereo completion failed: {error}"),
            Self::Construction(error) => {
                write!(formatter, "SMILES molecule construction failed: {error}")
            }
        }
    }
}

impl std::error::Error for SmilesError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Parse(error) => Some(error),
            Self::Hydrogen(error) => Some(error),
            Self::Sanitize(error) => Some(error),
            Self::Stereo(error) => Some(error),
            Self::Construction(error) => Some(error),
        }
    }
}

impl From<cosmolkit_smiles::SmilesParseError> for SmilesError {
    fn from(error: cosmolkit_smiles::SmilesParseError) -> Self {
        Self::Parse(error)
    }
}

impl From<cosmolkit_core::HydrogenError> for SmilesError {
    fn from(error: cosmolkit_core::HydrogenError) -> Self {
        Self::Hydrogen(error)
    }
}

impl From<cosmolkit_core::SanitizeError> for SmilesError {
    fn from(error: cosmolkit_core::SanitizeError) -> Self {
        Self::Sanitize(error)
    }
}

impl From<OperationError> for SmilesError {
    fn from(error: OperationError) -> Self {
        Self::Construction(error)
    }
}

impl From<cosmolkit_smiles::SmilesStereoError> for SmilesError {
    fn from(error: cosmolkit_smiles::SmilesStereoError) -> Self {
        Self::Stereo(error)
    }
}

impl Molecule {
    /// Constructs a molecule from SMILES with the pinned source defaults.
    pub fn from_smiles(input: &str) -> Result<Self, SmilesError> {
        Self::from_smiles_with_params(input, &cosmolkit_smiles::SmilesParseParams::default())
    }

    /// Constructs a molecule from SMILES using explicit parser parameters.
    ///
    /// Parsing, sanitization, and hydrogen removal remain in their detached
    /// algorithm owners. A live molecule is installed only after every
    /// requested stage and final structural validation succeeds.
    pub fn from_smiles_with_params(
        input: &str,
        params: &cosmolkit_smiles::SmilesParseParams,
    ) -> Result<Self, SmilesError> {
        let cosmolkit_smiles::SmilesRecord {
            mut topology,
            mut coordinates,
            mut properties,
        } = cosmolkit_smiles::parse_smiles(input, params)?;

        if params.remove_hydrogens {
            let remove_params = cosmolkit_core::RemoveHsParams {
                update_explicit_count: true,
                sanitize: params.sanitize,
                ..cosmolkit_core::RemoveHsParams::default()
            };
            let result = cosmolkit_core::remove_hydrogens_with_params(
                topology,
                coordinates,
                properties,
                &remove_params,
            )?;
            topology = result.topology;
            coordinates = result.coordinates;
            properties = result.properties;
        } else if params.sanitize {
            topology = cosmolkit_core::sanitize_topology(
                &topology,
                &cosmolkit_core::SanitizeParams::default(),
            )?
            .topology;
        }

        let record = cosmolkit_smiles::finalize_smiles_stereo(
            cosmolkit_smiles::SmilesRecord {
                topology,
                coordinates,
                properties,
            },
            params,
        )?;
        Self::from_parts(record.topology, record.coordinates, record.properties)
            .map_err(SmilesError::Construction)
    }
}
