//! Core aliases for the shared concrete atom model.
//!
//! Storage and local value operations live in `cosmolkit-model`. The core
//! specializes the opaque query payload with its search-layer representation;
//! this module intentionally contains no second atom implementation.

pub use cosmolkit_model::AtomId;
pub use cosmolkit_model::{
    AtomPdbResidueInfo, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element, ElementInfo,
    ElementParseError, Hybridization,
};

pub type AtomSpec = cosmolkit_model::AtomSpec<crate::QueryNode<crate::AtomQueryPredicate>>;
pub type Atom = cosmolkit_model::Atom<crate::QueryNode<crate::AtomQueryPredicate>>;
