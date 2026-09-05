//! Compatibility re-exports for query values now owned by `cosmolkit-model`.
//!
//! Search behavior remains in `crate::search`; this module is not a second
//! query model.

pub use cosmolkit_model::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryAtom, QueryBond, QueryGraph, QueryGraphError, QueryNode, RecursiveStructureQuery,
};
