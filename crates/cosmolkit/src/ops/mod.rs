//! Runtime operation system.

#[cfg(feature = "cap-aromaticity")]
mod aromaticity;
#[cfg(feature = "cap-stereo")]
mod cip_labels;
#[cfg(test)]
mod cow_tests;
#[cfg(feature = "cap-depict")]
mod depict;
mod error;
#[cfg(feature = "cap-hydrogens")]
mod hydrogens;
#[cfg(feature = "cap-kekulize")]
mod kekulize;
mod metadata;
#[cfg(feature = "cap-stereo")]
mod potential_stereo;
#[cfg(feature = "cap-radicals")]
mod radicals;
#[cfg(feature = "cap-rings")]
mod rings;
mod runtime;
#[allow(unexpected_cfgs)]
#[cfg(cosmolkit_runtime_privacy_probe)]
mod runtime_privacy_probe;
#[cfg(feature = "cap-sanitize")]
mod sanitize;
#[cfg(feature = "cap-stereo")]
mod structure_tags;
#[cfg(feature = "cap-transforms")]
mod transforms;
#[cfg(feature = "cap-valence")]
mod valence;

pub use error::OperationError;
pub use metadata::{
    BlockAccess, BlockSet, CipStatePolicy, DerivedEffects, DerivedState, FeatureSpec,
    FeatureSpecIter, MappingRequirement, MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec,
    OperationDomain, OperationInvariantEntry, ParityMatrixEntry, ParityPolicy,
    SemanticPreconditionSet, SupportMatrixEntry, TopologyEditKind, UnsupportedFeatureError,
    feature_spec, feature_specs, operation_invariant, operation_invariant_matrix, operation_parity,
    operation_spec, operation_specs, parity_matrix, support_matrix,
};
#[cfg(feature = "cap-stereo")]
pub use potential_stereo::PotentialStereoResult;
pub use runtime::registry::{
    MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX, SUPPORT_MATRIX,
};

pub use crate::FunctionStatus;
pub(crate) use runtime::context::{OpParts, PreservationProof};
pub(crate) use runtime::context::{PendingMolecule, PendingResult, ResultFinalizer};
pub(crate) use runtime::multiple::MultiOutputOpParts;
#[cfg(feature = "cap-stereo")]
pub(crate) use runtime::registry::PotentialStereoAccess;
#[cfg(feature = "cap-sanitize")]
pub(crate) use runtime::registry::SanitizeAccess;
#[cfg(feature = "cap-depict")]
pub(crate) use runtime::registry::With2dCoordinatesAccess;
#[cfg(feature = "cap-aromaticity")]
pub(crate) use runtime::registry::WithAssignedAromaticityAccess;
#[cfg(feature = "cap-radicals")]
pub(crate) use runtime::registry::WithAssignedRadicalsAccess;
#[cfg(feature = "cap-valence")]
pub(crate) use runtime::registry::WithAssignedValenceAccess;
#[cfg(feature = "cap-transforms")]
pub(crate) use runtime::registry::WithAtomPositionAccess;
#[cfg(feature = "cap-stereo")]
pub(crate) use runtime::registry::WithChiralTagsFromStructureAccess;
#[cfg(feature = "cap-stereo")]
pub(crate) use runtime::registry::WithCipLabelsAccess;
#[cfg(feature = "cap-kekulize")]
pub(crate) use runtime::registry::WithKekulizedBondsAccess;
#[cfg(test)]
pub(crate) use runtime::registry::{
    CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess,
};
#[cfg(feature = "cap-rings")]
pub(crate) use runtime::registry::{WithAssignedRingFamiliesAccess, WithAssignedRingsAccess};
#[cfg(feature = "cap-hydrogens")]
pub(crate) use runtime::registry::{WithHydrogensAccess, WithoutHydrogensAccess};
