use crate::{ReactionRole, ReactionValidationReport};
use cosmolkit_model::AtomId;

#[derive(Debug, thiserror::Error)]
pub enum ReactionValidationError {
    #[error("{role:?} template {template} atom {atom} property {property}: {source}")]
    Property {
        role: ReactionRole,
        template: usize,
        atom: AtomId,
        property: &'static str,
        #[source]
        source: cosmolkit_core::PropertyIntReadError,
    },
    #[error(
        "{role:?} template {template} atom {atom} property {property} unsigned map {value} exceeds signed source range"
    )]
    MapOverflow {
        role: ReactionRole,
        template: usize,
        atom: AtomId,
        property: &'static str,
        value: u32,
    },
    #[error("source invariant missing atom for product {template}, atom {atom}, map {map}")]
    MissingReactingAtom {
        template: usize,
        atom: AtomId,
        map: i32,
    },
    #[error("product {template} atom {atom} signed query arithmetic overflow")]
    QueryArithmetic { template: usize, atom: AtomId },
    #[error("product {template} atom {atom} annotation {property}: {source}")]
    Annotation {
        template: usize,
        atom: AtomId,
        property: &'static str,
        #[source]
        source: cosmolkit_model::AtomPropertyError,
    },
}

#[derive(Debug, thiserror::Error)]
pub enum ReactionInitializationError {
    #[error("reaction validation execution failed: {source}")]
    Validation {
        #[source]
        source: ReactionValidationError,
    },
    #[error("initialization failed with {num_errors} validation errors", num_errors = .report.num_errors())]
    Invalid { report: ReactionValidationReport },
}
