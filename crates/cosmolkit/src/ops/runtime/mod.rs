//! Private molecule-operation runtime boundary.
//!
//! Operation bodies are siblings of this module. Only generated, marker-
//! specific capabilities cross this boundary; unrestricted block access and
//! wrapper-owned transaction lifecycle methods remain visible solely within
//! this subtree.

pub(super) use super::{
    BlockSet, CipStatePolicy, DerivedState, FeatureSpec, MappingRequirement, MoleculeOpOutput,
    MoleculeOpSpec, OperationError, SupportStatus, TopologyEditKind,
};

#[path = "../context.rs"]
pub(super) mod context;
#[path = "../multiple.rs"]
pub(super) mod multiple;
#[path = "../registry.rs"]
pub(super) mod registry;
