//! Core molecule state model.
//!
//! This directory owns storage, identifiers, builders, read-only operation
//! views, and invariants for the value-style `Molecule` core.

pub mod adjacency;
pub mod atom;
pub mod bond;
pub mod builder;
pub mod derived;
pub mod error;
pub(crate) mod invariants;
pub mod molecule;
pub mod query;
pub mod read_parts;
pub mod sgroup;
