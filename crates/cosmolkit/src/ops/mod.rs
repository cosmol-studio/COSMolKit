//! Runtime operation system.

#[cfg(feature = "aromaticity")]
mod aromaticity;
#[cfg(feature = "stereo")]
mod cip_labels;
#[cfg(test)]
mod cow_tests;
mod error;
#[cfg(feature = "hydrogens")]
mod hydrogens;
#[cfg(feature = "kekulize")]
mod kekulize;
mod metadata;
#[cfg(feature = "stereo")]
mod potential_stereo;
#[cfg(feature = "radicals")]
mod radicals;
#[cfg(feature = "rings")]
mod rings;
mod runtime;
#[allow(unexpected_cfgs)]
#[cfg(cosmolkit_runtime_privacy_probe)]
mod runtime_privacy_probe;
#[cfg(feature = "sanitize")]
mod sanitize;
#[cfg(feature = "stereo")]
mod structure_tags;
#[cfg(feature = "transforms")]
mod transforms;
#[cfg(feature = "valence")]
mod valence;

pub use error::OperationError;
pub use metadata::{
    BlockAccess, BlockSet, CipStatePolicy, DerivedEffects, DerivedState, FeatureSpec,
    FeatureSpecIter, MappingRequirement, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec,
    OperationDomain, OperationInvariantEntry, ParityMatrixEntry, ParityPolicy,
    SemanticPreconditionSet, SupportMatrixEntry, SupportStatus, TopologyEditKind,
    UnsupportedFeatureError, feature_spec, feature_specs, operation_invariant,
    operation_invariant_matrix, operation_parity, operation_spec, operation_specs, parity_matrix,
    support_matrix,
};
#[cfg(feature = "stereo")]
pub use potential_stereo::PotentialStereoResult;
pub use runtime::registry::{
    MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX, SUPPORT_MATRIX,
};

pub(crate) use runtime::context::{OpParts, PreservationProof};
pub(crate) use runtime::context::{PendingMolecule, PendingResult, ResultFinalizer};
pub(crate) use runtime::multiple::MultiOutputOpParts;
#[cfg(feature = "stereo")]
pub(crate) use runtime::registry::PotentialStereoAccess;
#[cfg(feature = "sanitize")]
pub(crate) use runtime::registry::SanitizeAccess;
#[cfg(feature = "aromaticity")]
pub(crate) use runtime::registry::WithAssignedAromaticityAccess;
#[cfg(feature = "radicals")]
pub(crate) use runtime::registry::WithAssignedRadicalsAccess;
#[cfg(feature = "valence")]
pub(crate) use runtime::registry::WithAssignedValenceAccess;
#[cfg(feature = "transforms")]
pub(crate) use runtime::registry::WithAtomPositionAccess;
#[cfg(feature = "stereo")]
pub(crate) use runtime::registry::WithChiralTagsFromStructureAccess;
#[cfg(feature = "stereo")]
pub(crate) use runtime::registry::WithCipLabelsAccess;
#[cfg(feature = "kekulize")]
pub(crate) use runtime::registry::WithKekulizedBondsAccess;
#[cfg(test)]
pub(crate) use runtime::registry::{
    CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess,
};
#[cfg(feature = "rings")]
pub(crate) use runtime::registry::{WithAssignedRingFamiliesAccess, WithAssignedRingsAccess};
#[cfg(feature = "hydrogens")]
pub(crate) use runtime::registry::{WithHydrogensAccess, WithoutHydrogensAccess};
