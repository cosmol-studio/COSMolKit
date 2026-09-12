//! Runtime operation system.

mod context;
#[cfg(test)]
mod cow_tests;
mod error;
#[cfg(feature = "hydrogens")]
mod hydrogens;
mod metadata;
mod multiple;
mod registry;

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
pub use registry::{MOLECULE_OPS, OPERATION_INVARIANT_MATRIX, PARITY_MATRIX, SUPPORT_MATRIX};

pub(crate) use context::{OpParts, PreservationProof};
#[cfg(test)]
pub(crate) use cow_tests::{CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess};
pub(crate) use multiple::MultiOutputOpParts;
#[cfg(feature = "hydrogens")]
pub(crate) use registry::{WithHydrogensAccess, WithoutHydrogensAccess};
