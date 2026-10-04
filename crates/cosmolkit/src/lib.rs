//! Public COSMolKit runtime.
//!
//! This crate owns the live [`Molecule`] value and the operation lifecycle.
//! Chemistry algorithms are intentionally not implemented in this migration
//! step. They will be moved behind this boundary one operation at a time and
//! will receive detached [`cosmolkit_model`] blocks rather than a live
//! molecule.
//!
//! # Cargo features
//! Default features enable `full`. Plain names such as `core`, `bio`, and
//! `fingerprints` select bundles; `cap-*` names select individual capabilities,
//! such as `cap-io` or `cap-kekulize`. With defaults disabled, `core` is not
//! implicit. Features compose additively and do not change operation behavior.
//! See the crate README for bundle membership and advanced selection examples.

pub mod binding_contract;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    FingerprintError, SparseCountFingerprint, SparseCountFingerprint32,
};
#[cfg(feature = "cap-descriptors")]
mod descriptors;
#[cfg(feature = "cap-depict")]
pub use cosmolkit_depict::{
    Compute2DCoordinatesParams as Coordinate2DParams, Coordinate2DLayoutError,
    Coordinate2DTemplateError, DepictError as Coordinate2DError,
};
#[cfg(feature = "cap-descriptors")]
pub use cosmolkit_descriptors::DescriptorError;
#[cfg(feature = "cap-descriptors")]
pub use descriptors::DescriptorReadError;
#[cfg(feature = "cap-matrices")]
mod matrices;
mod molecule;
mod molecule_builder;
pub mod ops;
#[cfg(feature = "cap-io")]
mod sdf;
#[cfg(feature = "cap-smiles")]
mod smiles;
mod strict;

pub use binding_contract::{
    BINDING_CONTRACT, BindingCallableContract, BindingContractEntry, BindingDefault, BindingItem,
    BindingKind, BindingOwner, BindingParameterContract, BindingReceiver, BindingTypeRole,
    FunctionStatus, StateModel,
};
#[cfg(feature = "cap-bio")]
mod bio;
#[cfg(feature = "cap-bio")]
pub use bio::{BioOperationError, BioSelection, BioStructure, Protein, ProteinReadError};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_bio::{
    AltLocLabel, AltLocRequest, AtomName, AtomSourceIds, BioAltLocGroupId, BioAssembly,
    BioAssemblyGenerator, BioAssemblyId, BioAssemblyOperator, BioAssemblySpecialKind, BioAtomId,
    BioAtomRow, BioCalcFlag, BioChainId, BioChainRow, BioCisPep, BioConnection, BioCoordinateBlock,
    BioCoordinateFormat, BioCrystalCell, BioCrystalInfo, BioEntityDbRef, BioEntityId, BioEntityRow,
    BioHelix, BioMetadata, BioModRes, BioModelId, BioModelRow, BioNcsOperator, BioResidueId,
    BioResidueRow, BioRowChainError, BioRowModelError, BioRowSpan, BioRowTraverseError,
    BioSelectionCopyCause, BioSelectionCopyError, BioSelectionMatchError, BioSheet,
    BioSiftsUnpResidue, BioStructureError, BioStructureParts, BioStructureSourceState,
    BioTransform, ChainKind, ChainSourceIds, EntityKind, EntitySourceIds, PdbAtomSerial,
    PdbChainId, PdbSeqId, PolymerKind, ProteinAtomIter, ProteinAtomRef, ProteinChainIter,
    ProteinChainRef, ProteinProjectionError, ProteinResidueIter, ProteinResidueRef,
    ProteinSelectionSummary, ResidueCode, ResidueInfo, ResidueInfoKind, ResidueKind, ResidueName,
    ResidueSequenceError, ResidueSourceIds, UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter,
    expand_one_letter_sequence, find_residue_info, find_residue_info_index, residue_code,
    residue_info, residue_info_checked,
};
#[cfg(feature = "cap-rings")]
pub use cosmolkit_core::RingSearchParams;
#[cfg(feature = "cap-hydrogens")]
pub use cosmolkit_core::{AddHsParams, HydrogenError, RemoveHsParams};
#[cfg(feature = "cap-aromaticity")]
pub use cosmolkit_core::{AromaticityError, AromaticityModel, AromaticityParams};
#[cfg(feature = "cap-transforms")]
pub use cosmolkit_core::{AtomPositionParams, TransformError};
#[cfg(feature = "cap-sanitize")]
pub use cosmolkit_core::{
    ChemistryProblem, ChemistryProblemError, ChemistryProblemReport, SanitizeError,
    SanitizeOperations, SanitizeParams, SanitizeStage,
};
#[cfg(feature = "cap-matrices")]
pub use cosmolkit_core::{DenseMatrix, DistanceMatrix3dParams, MatrixError};
#[cfg(feature = "cap-kekulize")]
pub use cosmolkit_core::{KekulizeError, KekulizeParams};
#[cfg(feature = "cap-stereo")]
pub use cosmolkit_core::{
    PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoError, PotentialStereoInfo,
    PotentialStereoParams, PotentialStereoSpecified, PotentialStereoType, RingStereoRelation,
    StereoError, StructureTagParams,
};
#[cfg(feature = "cap-valence")]
pub use cosmolkit_core::{ValenceError, ValenceModel, ValenceParams};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{
    BioMmcifReadError, BioMmcifReadStage, BioMmcifWriteError, BioMmcifWriteParams, BioPdbReadError,
    BioPdbReadParams, BioPdbReadStage, BioReadError, BioReadParams,
};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{BioSelectionParseError, SelectionSeqidRangeError, SelectionSyntaxError};
pub use cosmolkit_model as model;
pub use cosmolkit_model::*;
#[cfg(feature = "cap-smiles")]
pub use cosmolkit_smiles::{
    CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, RandomSmilesWriteParams,
    SmilesParseParams, SmilesStereoError, SmilesWriteParams,
};
#[cfg(feature = "cap-stereo")]
pub use cosmolkit_stereo::{CipLabelOptions, CipLabelerError};
#[cfg(feature = "cap-matrices")]
pub use matrices::DistanceMatrixParams;
pub use molecule::Molecule;
pub use molecule_builder::MoleculeBuilder;
pub(crate) use ops::DerivedState;
#[cfg(feature = "cap-stereo")]
pub(crate) use ops::PotentialStereoAccess;
#[cfg(feature = "cap-stereo")]
pub use ops::PotentialStereoResult;
#[cfg(feature = "cap-sanitize")]
pub(crate) use ops::SanitizeAccess;
#[cfg(feature = "cap-depict")]
pub(crate) use ops::With2dCoordinatesAccess;
#[cfg(feature = "cap-aromaticity")]
pub(crate) use ops::WithAssignedAromaticityAccess;
#[cfg(feature = "cap-radicals")]
pub(crate) use ops::WithAssignedRadicalsAccess;
#[cfg(feature = "cap-valence")]
pub(crate) use ops::WithAssignedValenceAccess;
#[cfg(feature = "cap-transforms")]
pub(crate) use ops::WithAtomPositionAccess;
#[cfg(feature = "cap-stereo")]
pub(crate) use ops::WithChiralTagsFromStructureAccess;
#[cfg(feature = "cap-stereo")]
pub(crate) use ops::WithCipLabelsAccess;
#[cfg(feature = "cap-kekulize")]
pub(crate) use ops::WithKekulizedBondsAccess;
pub use ops::{
    BlockAccess, BlockSet, FeatureSpec, FeatureSpecIter, MOLECULE_OPS, MoleculeOpKind,
    MoleculeOpOutput, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OperationDomain, OperationError,
    OperationInvariantEntry, PARITY_MATRIX, ParityMatrixEntry, ParityPolicy, SUPPORT_MATRIX,
    SupportMatrixEntry, TopologyEditKind, UnsupportedFeatureError, feature_spec, feature_specs,
    operation_invariant, operation_invariant_matrix, operation_parity, operation_spec,
    operation_specs, parity_matrix, support_matrix,
};
#[cfg(test)]
pub(crate) use ops::{
    CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess,
    RingLiveCowCheckoutConflictForTestAccess,
};
pub(crate) use ops::{MultiOutputOpParts, OpParts, PreservationProof};
pub(crate) use ops::{PendingMolecule, PendingResult, ResultFinalizer};
#[cfg(feature = "cap-rings")]
pub(crate) use ops::{WithAssignedRingFamiliesAccess, WithAssignedRingsAccess};
#[cfg(feature = "cap-hydrogens")]
pub(crate) use ops::{WithHydrogensAccess, WithoutHydrogensAccess};
#[cfg(feature = "cap-io")]
pub use sdf::{SdfCoordinateMode, SdfError, SdfGraph, SdfReadParams, SdfRecord};
#[cfg(feature = "cap-smiles")]
pub use smiles::{
    FragmentCxSmilesWriteParams, FragmentSmilesWriteParams, SmilesError, SmilesWriteError,
};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
