//! Structured operation errors.

use std::fmt;

use cosmolkit_model::{CoordinateValidationError, TopologyValidationError};

/// Structured failure used while a capability is not yet implemented.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum OperationError {
    Unsupported {
        operation: &'static str,
    },
    AccessDenied {
        operation: &'static str,
        block: &'static str,
    },
    BlockCheckedOut {
        block: &'static str,
    },
    IncompleteCommit {
        block: &'static str,
    },
    InvalidTopology(TopologyValidationError),
    InvalidCoordinates(CoordinateValidationError),
    Algorithm {
        operation: &'static str,
        detail: String,
    },
}

impl fmt::Display for OperationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Unsupported { operation } => {
                write!(formatter, "operation `{operation}` is not implemented")
            }
            Self::AccessDenied { operation, block } => {
                write!(formatter, "operation `{operation}` cannot access {block}")
            }
            Self::BlockCheckedOut { block } => write!(formatter, "{block} is checked out"),
            Self::IncompleteCommit { block } => write!(formatter, "{block} was not committed"),
            Self::InvalidTopology(error) => error.fmt(formatter),
            Self::InvalidCoordinates(error) => error.fmt(formatter),
            Self::Algorithm { operation, detail } => {
                write!(formatter, "operation `{operation}` failed: {detail}")
            }
        }
    }
}

impl std::error::Error for OperationError {}
