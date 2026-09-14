//! Structured operation errors.

use std::fmt;

use cosmolkit_model::{
    CoordinateValidationError, MappingValidationError, MoleculePropertyError, TopologyEditError,
    TopologyValidationError,
};

use super::{
    CipStatePolicy, DerivedState, MappingRequirement, MoleculeOpOutput, MoleculeOpSpec,
    TopologyEditKind, UnsupportedFeatureError,
};

/// Structured failure used while a capability is not yet implemented.
#[derive(Clone, Debug, PartialEq)]
pub enum OperationError {
    UnsupportedFeature {
        operation: &'static MoleculeOpSpec,
        source: UnsupportedFeatureError,
    },
    Unsupported {
        operation: &'static str,
    },
    OutputMismatch {
        operation: &'static str,
        expected: MoleculeOpOutput,
        actual: MoleculeOpOutput,
    },
    AccessDenied {
        operation: &'static str,
        block: &'static str,
    },
    BlockCheckedOut {
        operation: &'static str,
        block: &'static str,
    },
    BlockNotCheckedOut {
        operation: &'static str,
        block: &'static str,
    },
    IncompleteCommit {
        operation: &'static str,
        block: &'static str,
    },
    TopologyEditContract {
        operation: &'static str,
        issue: &'static str,
        expected: TopologyEditKind,
        actual: Option<TopologyEditKind>,
    },
    MappingContract {
        operation: &'static str,
        issue: &'static str,
        requirement: MappingRequirement,
    },
    InvalidTopologyMapping {
        operation: &'static str,
        source: MappingValidationError,
    },
    AutoRemapContract {
        operation: &'static str,
        block: &'static str,
        issue: &'static str,
    },
    OperationContract {
        operation: &'static str,
        field: &'static str,
        issue: &'static str,
        expected: u8,
        actual: u8,
    },
    SemanticPreconditionContract {
        operation: &'static str,
        missing: super::SemanticPreconditionSet,
        issue: &'static str,
    },
    CoordinateAppendRequiresValues {
        operation: &'static str,
    },
    DerivedEffectContract {
        operation: &'static str,
        action: &'static str,
        states: DerivedState,
        issue: &'static str,
    },
    CipStateContract {
        operation: &'static str,
        policy: CipStatePolicy,
        issue: &'static str,
    },
    InvalidTopology(TopologyValidationError),
    InvalidTopologyEdit(TopologyEditError),
    InvalidCoordinates(CoordinateValidationError),
    InvalidProperty(MoleculePropertyError),
    InvalidPropertyList {
        target: &'static str,
        name: String,
        values: usize,
        expected: usize,
    },
    InvalidDerivedCache {
        state: &'static str,
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    #[cfg(any(feature = "valence", feature = "stereo"))]
    Valence(cosmolkit_core::ValenceError),
    #[cfg(feature = "radicals")]
    Radical(cosmolkit_core::RadicalError),
    #[cfg(any(feature = "rings", feature = "stereo", feature = "aromaticity"))]
    Rings(cosmolkit_core::RingFindingError),
    #[cfg(feature = "stereo")]
    PotentialStereo(cosmolkit_core::PotentialStereoError),
    #[cfg(feature = "stereo")]
    Stereo(cosmolkit_core::StereoError),
    #[cfg(feature = "stereo")]
    CipLabeler(cosmolkit_stereo::CipLabelerError),
    #[cfg(feature = "transforms")]
    Transform(cosmolkit_core::TransformError),
    #[cfg(feature = "kekulize")]
    Kekulize(cosmolkit_core::KekulizeError),
    #[cfg(feature = "aromaticity")]
    Aromaticity(cosmolkit_core::AromaticityError),
    #[cfg(feature = "sanitize")]
    Sanitize(cosmolkit_core::SanitizeError),
    #[cfg(feature = "hydrogens")]
    Hydrogen(cosmolkit_core::HydrogenError),
    InvalidAlgorithmResult {
        operation: &'static str,
        field: &'static str,
        actual: usize,
        expected: usize,
    },
    Algorithm {
        operation: &'static str,
        detail: String,
    },
}

impl fmt::Display for OperationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedFeature { operation, source } => write!(
                formatter,
                "operation `{}` cannot run because feature `{}` is unsupported: {}",
                operation.method, source.feature, source.reason
            ),
            Self::Unsupported { operation } => {
                write!(formatter, "operation `{operation}` is not implemented")
            }
            Self::OutputMismatch {
                operation,
                expected,
                actual,
            } => write!(
                formatter,
                "operation `{operation}` has output {actual:?}, expected {expected:?}"
            ),
            Self::AccessDenied { operation, block } => {
                write!(formatter, "operation `{operation}` cannot access {block}")
            }
            Self::BlockCheckedOut { operation, block } => {
                write!(formatter, "operation `{operation}` has {block} checked out")
            }
            Self::BlockNotCheckedOut { operation, block } => write!(
                formatter,
                "operation `{operation}` cannot install {block} without a checkout"
            ),
            Self::IncompleteCommit { operation, block } => {
                write!(formatter, "operation `{operation}` did not commit {block}")
            }
            Self::TopologyEditContract {
                operation,
                issue,
                expected,
                actual,
            } => write!(
                formatter,
                "operation `{operation}` has invalid topology-edit evidence ({issue}): expected {expected:?}, got {actual:?}"
            ),
            Self::MappingContract {
                operation,
                issue,
                requirement,
            } => write!(
                formatter,
                "operation `{operation}` violates mapping requirement {requirement:?}: {issue}"
            ),
            Self::InvalidTopologyMapping { operation, source } => write!(
                formatter,
                "operation `{operation}` supplied an invalid topology mapping: {source}"
            ),
            Self::AutoRemapContract {
                operation,
                block,
                issue,
            } => write!(
                formatter,
                "operation `{operation}` cannot auto-remap {block}: {issue}"
            ),
            Self::OperationContract {
                operation,
                field,
                issue,
                expected,
                actual,
            } => write!(
                formatter,
                "operation `{operation}` has invalid {field} contract ({issue}): expected mask {expected:#x}, got {actual:#x}"
            ),
            Self::SemanticPreconditionContract {
                operation,
                missing,
                issue,
            } => write!(
                formatter,
                "operation `{operation}` lacks semantic-precondition evidence for mask {:#x}: {issue}",
                missing.bits()
            ),
            Self::CoordinateAppendRequiresValues { operation } => write!(
                formatter,
                "operation `{operation}` appends atoms to populated conformers without complete coordinates"
            ),
            Self::DerivedEffectContract {
                operation,
                action,
                states,
                issue,
            } => write!(
                formatter,
                "operation `{operation}` has invalid derived-effect {action} for mask {:#x}: {issue}",
                states.bits()
            ),
            Self::CipStateContract {
                operation,
                policy,
                issue,
            } => write!(
                formatter,
                "operation `{operation}` has invalid CIP transition {policy:?}: {issue}"
            ),
            Self::InvalidTopology(error) => error.fmt(formatter),
            Self::InvalidTopologyEdit(error) => error.fmt(formatter),
            Self::InvalidCoordinates(error) => error.fmt(formatter),
            Self::InvalidProperty(error) => error.fmt(formatter),
            Self::InvalidPropertyList {
                target,
                name,
                values,
                expected,
            } => write!(
                formatter,
                "{target} SDF property list `{name}` has {values} rows, expected {expected}"
            ),
            Self::InvalidDerivedCache {
                state,
                field,
                actual,
                expected,
            } => write!(
                formatter,
                "derived cache state `{state}` field `{field}` has {actual} entries, expected {expected}"
            ),
            #[cfg(any(feature = "valence", feature = "stereo"))]
            Self::Valence(error) => write!(formatter, "valence assignment failed: {error}"),
            #[cfg(feature = "radicals")]
            Self::Radical(error) => write!(formatter, "radical assignment failed: {error}"),
            #[cfg(any(feature = "rings", feature = "stereo", feature = "aromaticity"))]
            Self::Rings(error) => write!(formatter, "ring assignment failed: {error}"),
            #[cfg(feature = "stereo")]
            Self::PotentialStereo(error) => {
                write!(formatter, "potential-stereo perception failed: {error}")
            }
            #[cfg(feature = "stereo")]
            Self::Stereo(error) => write!(formatter, "structure-tag assignment failed: {error}"),
            #[cfg(feature = "stereo")]
            Self::CipLabeler(error) => write!(formatter, "CIP label assignment failed: {error}"),
            #[cfg(feature = "transforms")]
            Self::Transform(error) => write!(formatter, "coordinate transform failed: {error}"),
            #[cfg(feature = "kekulize")]
            Self::Kekulize(error) => write!(formatter, "kekulization failed: {error}"),
            #[cfg(feature = "aromaticity")]
            Self::Aromaticity(error) => write!(formatter, "aromaticity assignment failed: {error}"),
            #[cfg(feature = "sanitize")]
            Self::Sanitize(error) => write!(formatter, "sanitization failed: {error}"),
            #[cfg(feature = "hydrogens")]
            Self::Hydrogen(error) => write!(formatter, "hydrogen transformation failed: {error}"),
            Self::InvalidAlgorithmResult {
                operation,
                field,
                actual,
                expected,
            } => write!(
                formatter,
                "operation `{operation}` returned {actual} {field} rows, expected {expected}"
            ),
            Self::Algorithm { operation, detail } => {
                write!(formatter, "operation `{operation}` failed: {detail}")
            }
        }
    }
}

impl std::error::Error for OperationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            #[cfg(any(feature = "valence", feature = "stereo"))]
            Self::Valence(error) => Some(error),
            #[cfg(feature = "radicals")]
            Self::Radical(error) => Some(error),
            #[cfg(any(feature = "rings", feature = "stereo", feature = "aromaticity"))]
            Self::Rings(error) => Some(error),
            #[cfg(feature = "stereo")]
            Self::PotentialStereo(error) => Some(error),
            #[cfg(feature = "stereo")]
            Self::Stereo(error) => Some(error),
            #[cfg(feature = "stereo")]
            Self::CipLabeler(error) => Some(error),
            #[cfg(feature = "transforms")]
            Self::Transform(error) => Some(error),
            #[cfg(feature = "kekulize")]
            Self::Kekulize(error) => Some(error),
            #[cfg(feature = "aromaticity")]
            Self::Aromaticity(error) => Some(error),
            #[cfg(feature = "sanitize")]
            Self::Sanitize(error) => Some(error),
            #[cfg(feature = "hydrogens")]
            Self::Hydrogen(error) => Some(error),
            _ => None,
        }
    }
}
