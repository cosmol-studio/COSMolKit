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
#[cfg(feature = "matrices")]
mod matrices;
mod molecule;
mod molecule_builder;
pub mod ops;
mod strict;

pub use binding_contract::{
    BINDING_CONTRACT, BindingCallableContract, BindingContractEntry, BindingDefault,
    BindingExposure, BindingItem, BindingKind, BindingOwner, BindingParameterContract,
    BindingParity, BindingSupport, BindingTypeRole, StateModel,
};
#[cfg(feature = "rings")]
pub use cosmolkit_core::RingSearchParams;
#[cfg(feature = "hydrogens")]
pub use cosmolkit_core::{AddHsParams, HydrogenError, RemoveHsParams};
#[cfg(feature = "aromaticity")]
pub use cosmolkit_core::{AromaticityError, AromaticityModel, AromaticityParams};
#[cfg(feature = "transforms")]
pub use cosmolkit_core::{AtomPositionParams, TransformError};
#[cfg(feature = "sanitize")]
pub use cosmolkit_core::{
    ChemistryProblem, ChemistryProblemError, ChemistryProblemReport, SanitizeError,
    SanitizeOperations, SanitizeParams, SanitizeStage,
};
#[cfg(feature = "matrices")]
pub use cosmolkit_core::{DenseMatrix, DistanceMatrix3dParams, MatrixError};
#[cfg(feature = "kekulize")]
pub use cosmolkit_core::{KekulizeError, KekulizeParams};
#[cfg(feature = "stereo")]
pub use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoError, PotentialStereoInfo,
    PotentialStereoParams, PotentialStereoSpecified, PotentialStereoType, RingStereoRelation,
    StereoError, StructureTagParams,
};
#[cfg(feature = "valence")]
pub use cosmolkit_core::{ValenceError, ValenceModel, ValenceParams};
pub use cosmolkit_model as model;
pub use cosmolkit_model::*;
#[cfg(feature = "matrices")]
pub use matrices::DistanceMatrixParams;
pub use molecule::Molecule;
pub use molecule_builder::MoleculeBuilder;
pub(crate) use ops::DerivedState;
#[cfg(feature = "stereo")]
pub(crate) use ops::PotentialStereoAccess;
#[cfg(feature = "stereo")]
pub use ops::PotentialStereoResult;
#[cfg(feature = "sanitize")]
pub(crate) use ops::SanitizeAccess;
#[cfg(feature = "aromaticity")]
pub(crate) use ops::WithAssignedAromaticityAccess;
#[cfg(feature = "radicals")]
pub(crate) use ops::WithAssignedRadicalsAccess;
#[cfg(feature = "valence")]
pub(crate) use ops::WithAssignedValenceAccess;
#[cfg(feature = "transforms")]
pub(crate) use ops::WithAtomPositionAccess;
#[cfg(feature = "stereo")]
pub(crate) use ops::WithChiralTagsFromStructureAccess;
#[cfg(feature = "kekulize")]
pub(crate) use ops::WithKekulizedBondsAccess;
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
#[cfg(feature = "rings")]
pub(crate) use ops::{WithAssignedRingFamiliesAccess, WithAssignedRingsAccess};
#[cfg(feature = "hydrogens")]
pub(crate) use ops::{WithHydrogensAccess, WithoutHydrogensAccess};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
