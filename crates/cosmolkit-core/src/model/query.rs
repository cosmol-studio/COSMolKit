//! Query vocabulary exposed alongside the molecule data model.
//!
//! Query predicates and graph values are model-level data. Parsing, matching,
//! and compilation remain in `search`, but consumers can depend on this
//! vocabulary without importing parser implementation modules.

pub use crate::search::query::{AtomQueryPredicate, BondQueryPredicate, QueryNode};
pub use crate::search::query_graph::{
    CompiledQuery, MatchResult, McsResult, QueryAtom, QueryBond, QueryGraph, QueryGraphError,
};
