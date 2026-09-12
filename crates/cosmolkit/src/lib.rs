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
mod molecule_builder;
pub mod ops;
mod strict;

pub use binding_contract::{
    BINDING_CONTRACT, BindingCallableContract, BindingContractEntry, BindingDefault,
    BindingExposure, BindingItem, BindingKind, BindingOwner, BindingParameterContract,
    BindingParity, BindingSupport, BindingTypeRole, StateModel,
};
pub use cosmolkit_model as model;
pub use cosmolkit_model::*;
pub use molecule::Molecule;
pub use molecule_builder::MoleculeBuilder;
pub(crate) use ops::DerivedState;
pub use ops::{
    BlockAccess, BlockSet, FeatureSpec, FeatureSpecIter, MOLECULE_OPS, MoleculeOpKind,
    MoleculeOpOutput, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OperationDomain, OperationError,
    OperationInvariantEntry, PARITY_MATRIX, ParityMatrixEntry, ParityPolicy, SUPPORT_MATRIX,
    SupportMatrixEntry, SupportStatus, TopologyEditKind, UnsupportedFeatureError, feature_spec,
    feature_specs, operation_invariant, operation_invariant_matrix, operation_parity,
    operation_spec, operation_specs, parity_matrix, support_matrix,
};
#[cfg(test)]
pub(crate) use ops::{CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess};
pub(crate) use ops::{MultiOutputOpParts, OpParts, PreservationProof};
#[cfg(feature = "hydrogens")]
pub(crate) use ops::{WithHydrogensAccess, WithoutHydrogensAccess};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
