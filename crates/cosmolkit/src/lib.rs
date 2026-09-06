//! Public COSMolKit runtime.
//!
//! This crate owns the live [`Molecule`] value and the operation lifecycle.
//! Chemistry algorithms are intentionally not implemented in this migration
//! step. They will be moved behind this boundary one operation at a time and
//! will receive detached [`cosmolkit_model`] blocks rather than a live
//! molecule.

pub mod binding_contract;
#[cfg(feature = "descriptors")]
mod descriptors;
mod molecule;
pub mod ops;

pub use binding_contract::{
    BINDING_CONTRACT, BindingContractEntry, BindingKind, BindingOwner, ReturnKind, StateModel,
};
pub use cosmolkit_model as model;
pub use cosmolkit_model::*;
pub use molecule::Molecule;
#[cfg(feature = "hydrogens")]
pub(crate) use ops::OpParts;
pub use ops::{
    AccessMode, BlockAccess, MOLECULE_OPS, MoleculeOpSpec, OperationError, OperationOutput,
    TopologyEditKind,
};
#[cfg(feature = "hydrogens")]
pub(crate) use ops::{AddHydrogensAccess, RemoveHydrogensAccess};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
