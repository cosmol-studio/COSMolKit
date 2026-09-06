//! Detached chemistry algorithms for the crate migration.
//!
//! The crate is intentionally below the public runtime boundary. It accepts
//! owned values from `cosmolkit-model` and never accepts a live `Molecule`, an
//! operation context, or runtime cache state.

mod hydrogens;

pub use hydrogens::{
    AddHsParams, CoreOperationError, DetachedBlocks, RemoveHsParams, add_hydrogens_impl,
    add_hydrogens_with_params, remove_hydrogens_impl, remove_hydrogens_with_params,
};
