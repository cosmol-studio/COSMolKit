//! Public COSMolKit runtime.
//!
//! This crate owns the live [`Molecule`] value and the operation lifecycle.
//! Domain crates own chemistry algorithms and receive detached
//! [`cosmolkit_model`] blocks rather than a live molecule. Public adapters
//! delegate algorithms through declared operation capabilities.
//!
//! # Cargo features
//! Default features enable `full`. Plain names such as `core`, `bio`, and
//! `fingerprints` select bundles; `cap-*` names select individual capabilities,
//! such as `cap-io` or `cap-kekulize`. With defaults disabled, `core` is not
//! implicit. Features compose additively and do not change operation behavior.
//! See the crate README for bundle membership and advanced selection examples.

#[doc(hidden)]
pub mod binding_contract;
#[cfg(feature = "cap-forcefields")]
mod forcefields;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    AtomPairAtomInvariantsGenerator, AtomPairParams, Fingerprint, FingerprintAdditionalOutput,
    FingerprintError, FingerprintJsonError, MorganParams, SparseBitFingerprint,
    SparseCountFingerprint, SparseCountFingerprint32, TopologicalTorsionParams,
};
#[cfg(feature = "cap-forcefields")]
pub use forcefields::{
    MmffAtomProperties, MmffMolPropertiesError, MmffProperties, MmffPropertiesParams, MmffVariant,
    UffParameterError, UffParameterErrorKind, UffParameterQueryError,
};
#[cfg(feature = "cap-forcefields")]
pub(crate) use ops::WithUffOptimizedAccess;
#[cfg(feature = "cap-forcefields")]
pub(crate) use ops::WithUffOptimizedConfsAccess;
#[cfg(feature = "cap-forcefields")]
pub use ops::{
    UffConformerOptimizationParams, UffConformerOptimizationResult, UffConformerResult,
    UffOptimizationError, UffOptimizationErrorKind, UffOptimizationParams, UffOptimizationResult,
};
#[cfg(feature = "cap-descriptors")]
mod descriptors;
#[cfg(feature = "cap-depict")]
mod drawing;
#[cfg(feature = "cap-depict")]
pub use cosmolkit_depict::{
    Compute2DCoordinatesParams as Coordinate2DParams, Coordinate2DLayoutError,
    Coordinate2DTemplateError, DepictError as Coordinate2DError, DrawingError,
};
#[cfg(feature = "cap-descriptors")]
pub use cosmolkit_descriptors::{
    CrippenTotals, DescriptorError, LabuteAsaContributions, RotatableBondsOptions,
};
#[cfg(feature = "cap-descriptors")]
pub use descriptors::DescriptorReadError;
#[cfg(feature = "cap-depict")]
pub use drawing::DrawingWriteError;
#[cfg(feature = "cap-alignment")]
mod alignment;
#[cfg(feature = "cap-matrices")]
mod matrices;
mod molecule;
mod molecule_builder;
#[cfg(feature = "cap-fingerprints")]
mod morgan;
#[cfg(feature = "cap-fingerprints")]
mod morgan_reusable;
#[cfg(feature = "cap-fingerprints")]
pub use morgan_reusable::{
    MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganCallParams,
    MorganFingerprintGenerator, MorganSettings,
};
#[doc(hidden)]
pub mod ops;
#[cfg(feature = "cap-io")]
mod sdf;
#[cfg(feature = "cap-search")]
pub mod search;
#[cfg(feature = "cap-search")]
pub use search::{
    compile_query, parse_smarts, parse_smarts_with_params, write_cx_smarts, write_smarts,
};
#[cfg(feature = "cap-smiles")]
mod smiles;
mod strict;
#[cfg(feature = "cap-search")]
pub use cosmolkit_search::{
    CompiledQuery, MatchError, MatchResult, QueryCompileError, SmartsParseError, SmartsParseParams,
    SmartsWriteError, SmartsWriteParams, SubstructMatchError, SubstructMatchParams,
};

#[doc(hidden)]
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
    ProteinSelectionSummary, ResidueCode, ResidueCodeParseError, ResidueIdentity, ResidueInfo,
    ResidueInfoKind, ResidueKind, ResidueName, ResidueSequenceError, ResidueSourceIds,
    UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter, expand_one_letter_sequence,
    find_residue_info, find_residue_info_index, residue_code, residue_info, residue_info_checked,
};
#[cfg(feature = "cap-rings")]
pub use cosmolkit_core::RingSearchParams;
#[cfg(feature = "cap-hydrogens")]
pub use cosmolkit_core::{AddHsParams, HydrogenError, RemoveHsParams};
#[cfg(feature = "cap-aromaticity")]
pub use cosmolkit_core::{AromaticityError, AromaticityModel, AromaticityParams};
#[cfg(feature = "cap-valence")]
pub use cosmolkit_core::{ValenceModel, ValenceParams};
// Detached metadata and its error are always-present model vocabulary.
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
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{
    BioMmcifReadError, BioMmcifReadStage, BioMmcifWriteError, BioMmcifWriteParams, BioPdbReadError,
    BioPdbReadParams, BioPdbReadStage, BioPdbWriteError, BioPdbWriteParams, BioReadError,
    BioReadParams,
};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{BioSelectionParseError, SelectionSeqidRangeError, SelectionSyntaxError};
pub use cosmolkit_model as model;
pub use cosmolkit_model::*;
pub use cosmolkit_model::{AtomMetadata, ValenceError};
#[cfg(feature = "cap-smiles")]
pub use cosmolkit_smiles::{
    CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, RandomSmilesWriteParams,
    SmilesParseParams, SmilesStereoError, SmilesWriteParams,
};
#[cfg(feature = "cap-stereo")]
pub use cosmolkit_stereo::{CipLabelOptions, CipLabelerError, StereoReadError};
#[cfg(feature = "cap-matrices")]
pub use matrices::DistanceMatrixParams;
pub use molecule::Molecule;
pub use molecule_builder::MoleculeBuilder;
#[cfg(feature = "cap-fingerprints")]
pub use morgan::{MorganFingerprintParams, MorganInvariants, MorganReadError};
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_preparation;
#[cfg(feature = "cap-fingerprints")]
pub use fingerprint_preparation::FingerprintPreparationError;
#[cfg(feature = "cap-fingerprints")]
mod atom_pair;
#[cfg(feature = "cap-fingerprints")]
pub use atom_pair::{AtomPairFingerprintParams, AtomPairReadError};
#[cfg(feature = "cap-fingerprints")]
mod legacy_topological_torsion;
#[cfg(feature = "cap-fingerprints")]
pub use legacy_topological_torsion::LegacyTopologicalTorsionParams;
#[cfg(feature = "cap-fingerprints")]
mod topological_torsion;
pub(crate) use ops::DerivedState;
#[cfg(feature = "cap-stereo")]
pub(crate) use ops::PotentialStereoAccess;
#[cfg(feature = "cap-stereo")]
#[doc(inline)]
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
#[doc(hidden)]
pub use ops::{
    BlockAccess, BlockSet, FeatureSpec, FeatureSpecIter, MOLECULE_OPS, MoleculeOpKind,
    MoleculeOpOutput, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OperationDomain,
    OperationInvariantEntry, PARITY_MATRIX, ParityMatrixEntry, ParityPolicy, SUPPORT_MATRIX,
    SupportMatrixEntry, TopologyEditKind, feature_spec, feature_specs, operation_invariant,
    operation_invariant_matrix, operation_parity, operation_spec, operation_specs, parity_matrix,
    support_matrix,
};
#[cfg(test)]
pub(crate) use ops::{
    CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess,
    RingLiveCowCheckoutConflictForTestAccess,
};
pub(crate) use ops::{MultiOutputOpParts, OpParts, PreservationProof};
#[doc(inline)]
pub use ops::{OperationError, UnsupportedFeatureError};
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
#[cfg(feature = "cap-fingerprints")]
pub use topological_torsion::{
    TopologicalTorsionCallParams, TopologicalTorsionFingerprintGenerator,
    TopologicalTorsionFingerprintParams, TopologicalTorsionReadError, TopologicalTorsionSettings,
};

/// Returns RDKit periodic-table metadata for an element, including the dummy (`*`).
///
/// All elements with atomic numbers `0..=118` are supported. The record contains
/// the canonical symbol, period, outer electron count, complete ordered valence
/// list, `Rb0` bond radius in angstroms, and atomic weight. Source-defined zeros
/// and the unrestricted-valence sentinel `-1` are preserved.
///
/// The returned symbol and valence slice borrow immutable shared table data.
/// This query has no options or errors and does not change any molecule.
#[cfg(feature = "cap-valence")]
#[must_use]
pub fn element_info(element: Element) -> ElementInfo {
    // Source: RDKit 2026.03.1, Code/GraphMol/{PeriodicTable,atomic_data}.h.
    // These field-access anchors describe the delegated core table lookup;
    // Element's checked identity satisfies the numeric source preconditions.
    // RDKit✔️✔️: std::string getElementSymbol(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Symbol();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getAtomicWeight(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   double mass = byanum[atomicNumber].Mass();
    // RDKit✔️✔️:   return mass;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: double getRb0(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].Rb0();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const INT_VECT &getValenceList(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].ValenceList();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: int getNouterElecs(UINT atomicNumber) const {
    // RDKit✔️✔️:   PRECONDITION(atomicNumber < byanum.size(), "Atomic number not found");
    // RDKit✔️✔️:   return byanum[atomicNumber].NumOuterShellElec();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: int AtomicNum() const { return anum; }
    // RDKit✔️✔️: unsigned int Row() const { return row; }
    // Warm lookup is O(1), returns a Copy record and borrows the existing
    // static valence slice without allocation. Cold initialization remains
    // the single core-owned table initialization.
    cosmolkit_core::element_info(element)
}

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

#[cfg(feature = "cap-forcefields")]
pub use ops::{
    MmffConformerOptimizationParams, MmffOptimizationError, MmffOptimizationParams,
    MmffOptimizeMoleculeConfResult, MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult,
};

#[cfg(feature = "cap-forcefields")]
pub(crate) use ops::{WithMmffOptimizedAccess, WithMmffOptimizedConfsAccess};

#[cfg(feature = "cap-forcefields")]
pub use forcefields::{MmffEnergyGradient, MmffEvaluationParams};

#[cfg(feature = "cap-tautomer")]
mod tautomer;
#[cfg(feature = "cap-tautomer")]
pub(crate) use ops::{CanonicalTautomerWithParamsAccess, EnumerateTautomersWithParamsAccess};
#[cfg(feature = "cap-tautomer")]
pub use tautomer::{
    TautomerCatalogError, TautomerEnumeration, TautomerEnumerationCallback,
    TautomerEnumerationStatus, TautomerMoleculeView, TautomerParams, TautomerProgress,
    TautomerRunError, TautomerScore, TautomerScoreParams, TautomerScoreTerm, TautomerScorer,
    canonical_tautomer_from_molecules, canonical_tautomer_from_molecules_with_params,
    default_tautomer_score_terms,
};

#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
pub use bio::BioMoleculeError;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    AtomCodeExplanation, AtomCodeExplanationError, AtomPairsParameters,
};
#[cfg(feature = "cap-serialization")]
pub use cosmolkit_io::PickleError;
#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
pub use cosmolkit_io::{BioMoleculeConversionError, BioMoleculeParams};
#[cfg(feature = "cap-fingerprints")]
pub use ops::AtomPairAtomCodeResult;
#[cfg(feature = "cap-fingerprints")]
pub(crate) use ops::WithAtomPairAtomCodeAccess;
#[cfg(feature = "cap-forcefields")]
pub use ops::{UffEnergyGradient, UffEvaluationParams};

#[cfg(feature = "cap-hashing")]
mod molecular_hash;
#[cfg(feature = "cap-hashing")]
pub use cosmolkit_core::CipRankError;
#[cfg(feature = "cap-hashing")]
pub use cosmolkit_fingerprints::MoleculeHashError;

#[cfg(feature = "cap-fingerprints")]
mod topological_torsion_path;
#[cfg(feature = "cap-alignment")]
pub use cosmolkit_alignment::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AlignmentResult, AlignmentTransform,
    AllConformerRmsdParameters, BestAlignmentParameters, ConformerAlignmentParameters,
    ConformerAlignmentReport, ConformerRmsd, CoordinateRmsdParameters,
};
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{TopologicalTorsionPathScoreError, explain_path_score};

#[cfg(feature = "cap-alignment")]
pub(crate) use ops::{WithAlignedConformersAccess, WithAlignmentToAccess};
