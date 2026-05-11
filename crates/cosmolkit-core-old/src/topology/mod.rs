//! Lightweight internal helpers for topology-related molecule operations.
//!
//! Strong operations carry index mappings so dependent data can be remapped
//! from one factual source of truth. Weak operations keep indices stable but
//! still declare which derived state must be invalidated or recomputed.

#![allow(dead_code)]

mod mapping;
mod ops;
mod spec;

pub(crate) use mapping::{AtomMapping, BondMapping, IndexMapping};
pub(crate) use ops::{
    TopologyEditResult, apply_strong_topology_op, apply_topology_op, identity_then_created_mapping,
};
pub(crate) use spec::{
    KEKULIZE_SPEC, SANITIZE_SPEC, TopologyOpKind, TopologyOpSpec, WITH_HYDROGENS_SPEC,
    WITHOUT_HYDROGENS_SPEC,
};
