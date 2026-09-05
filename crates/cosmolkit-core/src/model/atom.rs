//! Core aliases for the shared concrete atom model.
//!
//! Storage and local value operations live in `cosmolkit-model`; this module
//! intentionally contains no second atom implementation or query payload.

pub use cosmolkit_model::AtomId;
pub use cosmolkit_model::{
    AtomPdbResidueInfo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element, ElementInfo,
    ElementParseError, Hybridization,
};

pub use cosmolkit_model::Atom;
pub use cosmolkit_model::AtomSpec;
