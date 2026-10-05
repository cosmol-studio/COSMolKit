//! Shared failures while preparing molecule state for fingerprint generation.

use cosmolkit_core::RingFindingError;
use std::fmt;

/// A preparation failure shared by all public fingerprint families.
#[derive(Debug)]
pub enum FingerprintPreparationError {
    /// The molecule has no valid prepared valence assignment.
    MissingPreparedValence,
    /// Detached ring preparation failed with its original core error.
    RingPreparation(RingFindingError),
}

impl fmt::Display for FingerprintPreparationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingPreparedValence => formatter
                .write_str("Fingerprint preparation requires a valid prepared valence assignment"),
            Self::RingPreparation(error) => {
                write!(formatter, "Fingerprint ring preparation failed: {error}")
            }
        }
    }
}

impl std::error::Error for FingerprintPreparationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::MissingPreparedValence => None,
            Self::RingPreparation(error) => Some(error),
        }
    }
}
