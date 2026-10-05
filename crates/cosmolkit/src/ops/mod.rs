//! Runtime operation system.

#[cfg(feature = "cap-aromaticity")]
mod aromaticity;
#[cfg(feature = "cap-stereo")]
mod cip_labels;
#[cfg(test)]
mod cow_tests;
#[cfg(all(test, feature = "cap-aromaticity"))]
pub(crate) use cow_tests::ring_aromaticity_probe;
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
#[cfg(feature = "cap-forcefields")]
mod uff_optimization;
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
#[cfg(feature = "cap-forcefields")]
pub(crate) use runtime::registry::WithUffOptimizedConformersAccess;
#[cfg(feature = "cap-forcefields")]
pub(crate) use runtime::registry::WithUffOptimizedCoordinatesAccess;
pub use runtime::registry::{
    MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX, SUPPORT_MATRIX,
};
#[cfg(feature = "cap-forcefields")]
pub use uff_optimization::{
    UffConformerOptimizationParams, UffConformerOptimizationResult, UffConformerResult,
    UffOptimizationError, UffOptimizationErrorKind, UffOptimizationParams, UffOptimizationResult,
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
    RingLiveCowCheckoutConflictForTestAccess,
};
#[cfg(feature = "cap-rings")]
pub(crate) use runtime::registry::{WithAssignedRingFamiliesAccess, WithAssignedRingsAccess};
#[cfg(feature = "cap-hydrogens")]
pub(crate) use runtime::registry::{WithHydrogensAccess, WithoutHydrogensAccess};

#[cfg(feature = "cap-forcefields")]
pub(crate) mod mmff_optimization;
#[cfg(feature = "cap-forcefields")]
pub use mmff_optimization::{
    MmffConformerOptimizationParams, MmffOptimizationError, MmffOptimizationParams,
    MmffOptimizeMoleculeConfResult, MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult,
};
#[cfg(feature = "cap-forcefields")]
pub(crate) use runtime::registry::{WithMmffOptimizedAccess, WithMmffOptimizedConfsAccess};
