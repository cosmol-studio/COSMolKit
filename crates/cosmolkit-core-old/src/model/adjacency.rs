//! Core re-exports for the shared adjacency value.
//!
//! The CSR representation and its local endpoint validation are owned by
//! `cosmolkit-model`; core keeps only molecule/runtime lifecycle checks.

pub use cosmolkit_model::{AdjacencyError, AdjacencyList, NeighborRef};
