//! Public COSMolKit runtime.
//!
//! This crate owns the live [`Molecule`] value and the operation lifecycle.
//! Domain crates own chemistry algorithms and receive detached
//! [`cosmolkit_model`] blocks rather than a live molecule. Public adapters
//! delegate algorithms through declared operation capabilities.
//!
//! # Cargo features
//! Default features enable `full`. Plain names such as `core`, `bio`, and
//! `fingerprints` select functionality and include its public prerequisites.
//! For example, `reaction` includes search and SMILES APIs. With defaults
//! disabled, `core` is not implicit. Features compose additively and do not
//! change operation behavior.
//! `core` includes molecular format parsing/writing for text and files;
//! native binary archives and ordered batch processing are included in `core`.
//! Binary archives are not compiled on wasm32.
//! See the dependency tree at the top of the crate README.

#[cfg(doctest)]
#[path = "../tests/support/feature_gate_docs.rs"]
mod feature_gate_docs;

#[doc(hidden)]
pub mod binding_contract;
#[cfg(feature = "cap-forcefields")]
mod forcefields;
#[cfg(feature = "cap-forcefields")]
pub use cosmolkit_core::FragmentCoordinateViewError;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    AtomPairAtomInvariantsGenerator, AtomPairError, AtomPairGeneratorParams, Fingerprint,
    FingerprintAdditionalOutput, FingerprintError, FingerprintJsonError, FingerprintWorkerError,
    MorganError, MorganGeneratorParams, SparseBitFingerprint, SparseCountFingerprint,
    SparseCountFingerprint32, TopologicalTorsionError, TopologicalTorsionGeneratorParams,
};
#[cfg(all(feature = "cap-forcefields", not(target_family = "wasm")))]
pub use cosmolkit_forcefields::DispatchedUffOptimizationError;
#[cfg(feature = "cap-forcefields")]
pub use cosmolkit_forcefields::{
    AutomaticForceFieldConstructionError, ForceFieldConstructionError, MmffBuilderError,
    MmffOptimizationError as MmffOptimizationSourceError, MmffOptimizationFailure,
    NonbondedAssemblyError, NonbondedFragmentMappingError, OptimizationStageError,
    SerialConformerOptimizationError, SerialCoordinateCountError, SerialResultCapacityError,
    SerialUffOptimizationError, SingleConformerOptimizationError, UffConformerError,
    UffEvaluationError, UffEvaluationFailure, UffPreparedOptimizationError, UffSingleError,
    UffThreadCountError,
};
#[cfg(feature = "cap-forcefields")]
pub use cosmolkit_forcefields::{
    ForceFieldAngleIndexArgument, ForceFieldAngleRangeBound, ForceFieldBondIndexArgument,
    ForceFieldIndexArgument, ForceFieldKernelError, ForceFieldTorsionIndexArgument,
    MmffAngleLookupError, MmffOopLookupError, MmffParamError, MmffParamParseCause,
    MmffParamParseError, MmffParamTable, MmffTorLookupError, MmffVdwParseError, UffBondMathError,
    UffBuilderError, UffDefaultTorsionQueryError, UffInversionContributionError,
    UffInversionIndexArgument, UffParamError, UffPreparedValenceField, UffTorsionBondQueryError,
    UffTypingError, UffTypingInput,
};
#[cfg(feature = "cap-forcefields")]
pub use forcefields::{
    MmffAtomProperties, MmffMolPropertiesError, MmffProperties, MmffPropertiesParams, MmffVariant,
    UffParameterError, UffParameterErrorKind, UffParameterQueryError,
};
#[cfg(feature = "cap-forcefields")]
pub(crate) use ops::WithUffOptimizedAccess;
#[cfg(feature = "cap-forcefields")]
pub(crate) use ops::WithUffOptimizedConformersAccess;
#[cfg(feature = "cap-forcefields")]
pub use ops::{
    UffConformerOptimizationParams, UffConformerOptimizationResult, UffConformerResult,
    UffOptimizationError, UffOptimizationErrorKind, UffOptimizationParams, UffOptimizationResult,
};
#[cfg(feature = "cap-descriptors")]
mod descriptors;
#[cfg(feature = "cap-depict")]
mod drawing;
#[cfg(feature = "cap-io")]
pub use cosmolkit_core::{AttachmentExpansionError, BondDirectionStereoError};
#[cfg(feature = "cap-depict")]
pub use cosmolkit_core::{
    DoubleLexicalReadError, DoubleLexicalReadErrorKind, SourceNumericStreamState,
    UnsignedStreamArrayError,
};
#[cfg(feature = "cap-depict")]
pub use cosmolkit_depict::{
    Compute2DCoordinatesParams as Coordinate2DParams, Coordinate2DLayoutError,
    Coordinate2DTemplateError, DepictError as Coordinate2DError, DrawingError, SvgParseError,
    XmlParseError, XmlTextPos,
};
#[cfg(feature = "cap-descriptors")]
pub use cosmolkit_descriptors::{
    CrippenTotals, DescriptorError, DescriptorSearchCause, LabuteAsaContributions,
    RotatableBondsOptions,
};
#[cfg(feature = "cap-io")]
pub use cosmolkit_io::{MolPostError, MolProcessingError};
#[cfg(feature = "cap-smiles")]
pub use cosmolkit_smiles::FragmentWriteInputError;
#[cfg(feature = "cap-descriptors")]
pub use descriptors::DescriptorReadError;
#[cfg(feature = "cap-depict")]
pub use drawing::DrawingWriteError;
#[cfg(feature = "cap-alignment")]
mod alignment;
#[cfg(feature = "cap-inchi")]
mod inchi;
#[cfg(feature = "cap-matrices")]
mod matrices;
mod molecule;
#[doc(hidden)]
pub use cosmolkit_macros::{Configuration, project_configuration};
mod molecule_builder;
#[cfg(feature = "cap-inchi")]
pub use cosmolkit_inchi::{InchiError, InchiErrorKind};
#[cfg(feature = "cap-inchi")]
pub use inchi::{InchiReadParams, InchiWriteParams, inchi_to_key};
#[cfg(feature = "cap-reaction")]
mod reaction;
#[cfg(all(feature = "cap-reaction", not(feature = "cap-smiles")))]
pub use cosmolkit_reaction::CxSmilesFields;
#[cfg(all(feature = "cap-reaction", not(feature = "cap-search")))]
pub use cosmolkit_reaction::SubstructMatchParams;
#[cfg(feature = "cap-reaction")]
pub use cosmolkit_reaction::{
    CxParseError, ReactionApplyError, ReactionApplyParams, ReactionCoordinateSelection,
    ReactionInitializationError, ReactionModelError, ReactionParseError, ReactionParseParams,
    ReactionProductError, ReactionRole, ReactionRunError, ReactionRunParams,
    ReactionSingleRunParams, ReactionTemplateRemovalParams, ReactionValidationError,
    ReactionValidationIssue, ReactionValidationIssueKind, ReactionValidationParams,
    ReactionValidationReport, ReactionValidationSeverity, ReactionWriteError, ReactionWriteParams,
};
#[cfg(feature = "cap-reaction")]
pub(crate) use ops::{
    ApplyReactionAccess, ReactionProductsAccess, ReactionProductsFromInputsAccess, RunAccess,
};
#[cfg(feature = "cap-reaction")]
pub use reaction::{
    Reaction, ReactionApplyResult, ReactionTemplateRemoval, parse_smirks, parse_smirks_with_params,
};
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_batch;
#[cfg(feature = "cap-fingerprints")]
mod fingerprint_params;
#[cfg(feature = "cap-fingerprints")]
pub use fingerprint_batch::{FingerprintBatch, FingerprintBatchError};
#[cfg(feature = "cap-fingerprints")]
mod morgan;
#[cfg(feature = "cap-fingerprints")]
mod morgan_reusable;
#[cfg(feature = "cap-fingerprints")]
pub use morgan_reusable::{
    MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganCallParams,
    MorganFingerprintGenerator, MorganSettings,
};
#[cfg(feature = "cap-batch")]
mod batch;
#[doc(hidden)]
pub mod ops;
#[cfg(feature = "cap-io")]
pub use cosmolkit_io::SdfReadError;
#[cfg(feature = "cap-tautomer")]
pub(crate) use ops::{WithAssignedSymmSssrAccess, WithInstalledTautomerScoreCacheAccess};
#[cfg(feature = "cap-io")]
mod sdf;
#[cfg(feature = "cap-batch")]
pub use batch::{
    BatchError, BatchErrorMode, BatchParams, BatchQueryParams, BatchRecord, BatchValidationError,
    MoleculeBatch,
};
#[cfg(all(feature = "cap-batch", feature = "cap-depict"))]
pub use batch::{BatchImageError, BatchImageParams};
#[cfg(feature = "cap-batch")]
pub use cosmolkit_batch::BatchExportReport;
#[cfg(all(feature = "cap-io", feature = "cap-batch"))]
mod sdf_batch;
#[cfg(all(feature = "cap-io", feature = "cap-batch"))]
pub use sdf_batch::{BatchExportParams, SdfBatchIterator, SdfReaderBatchIterator};
#[cfg(feature = "cap-search")]
pub mod search;
#[cfg(feature = "cap-search")]
pub use search::{
    compile_query, maximum_common_substructure, maximum_common_substructure_with_params,
    parse_smarts, parse_smarts_with_params, write_cx_smarts, write_smarts,
};
#[cfg(feature = "cap-smiles")]
mod smiles;
mod strict;
#[cfg(feature = "cap-search")]
pub use cosmolkit_core::WedgeError;
#[cfg(feature = "cap-search")]
pub use cosmolkit_core::parser_stereo_order::CxStereoGroupMergeError;
#[cfg(feature = "cap-search")]
pub use cosmolkit_search::{
    AtomComparator as McsAtomComparator, BondComparator as McsBondComparator, CompiledQuery,
    CxQueryLoweringError, MatchError, MatchResult, McsAtomCompareParameters,
    McsBondCompareParameters, McsCandidateMatchError as McsError, McsError as SearchMcsError,
    McsParameters, McsProgressError, McsResult, QueryCompileError, QueryMatchContextError,
    SmartsParseError, SmartsParseParams, SmartsWriteError, SmartsWriteParams, SubstructMatchError,
    SubstructMatchParams,
};
#[cfg(feature = "cap-search")]
pub use cosmolkit_smiles::SmartsTraversalError;

pub use binding_contract::conveniences::*;
#[doc(hidden)]
pub use binding_contract::{
    BINDING_CONTRACT, BINDING_CONTRACT_CONVENIENCES, BINDING_CONTRACT_KEYWORDS,
    BINDING_CONTRACT_PROPERTIES, BINDING_CONTRACT_PYTHON_ADAPTERS,
    BINDING_CONTRACT_PYTHON_COLLECTIONS, BindingCallableContract, BindingConfigurationField,
    BindingContractEntry, BindingConvenienceContract, BindingDefault, BindingItem,
    BindingKeywordContract, BindingKind, BindingOwner, BindingParameterContract,
    BindingPropertyAccess, BindingPropertyContract, BindingPythonAdapterContract, BindingReceiver,
    BindingTypeRole, FunctionStatus, StateModel,
};
#[cfg(feature = "cap-bio")]
mod bio;
#[cfg(feature = "cap-bio")]
pub use bio::{BioOperationError, BioSelection, BioStructure, Protein, ProteinReadError};
#[cfg(feature = "cap-bio")]
#[doc(hidden)]
pub use cosmolkit_bio::__ck_configuration_BioCrystalCell;
#[cfg(feature = "cap-bio")]
pub use cosmolkit_bio::{
    AltLocLabel, AltLocRequest, AtomAddress, AtomName, AtomSourceIds, BioAltLocGroupId,
    BioAssembly, BioAssemblyGenerator, BioAssemblyId, BioAssemblyOperator, BioAssemblySpecialKind,
    BioAsu, BioAtomId, BioAtomRow, BioBasicRefinementInfo, BioCalcFlag, BioChainId, BioChainRow,
    BioCisPep, BioConnection, BioConnectionKind, BioCoordinateBlock, BioCoordinateFormat,
    BioCrystalCell, BioCrystalInfo, BioDiffractionInfo, BioEntityDbRef, BioEntityId, BioEntityRow,
    BioExperimentInfo, BioExperimentalCrystalInfo, BioHelix, BioHelixClass, BioMetadata, BioModRes,
    BioModelId, BioModelRow, BioNcsOperator, BioRefinementInfo, BioRefinementRestraint,
    BioReflectionsInfo, BioResidueId, BioResidueRow, BioRowChainError, BioRowModelError,
    BioRowSpan, BioRowTraverseError, BioSelectionCopyCause, BioSelectionCopyError,
    BioSelectionMatchError, BioSheet, BioSiftsUnpResidue, BioSoftwareClassification,
    BioSoftwareItem, BioStrand, BioStructureError, BioStructureParts, BioStructureSourceState,
    BioTlsGroup, BioTlsSelection, BioTransform, ChainKind, ChainSourceIds, EntityKind,
    EntitySourceIds, PdbAtomSerial, PdbChainId, PdbSeqId, PolymerKind, ProteinAtomIter,
    ProteinAtomRef, ProteinChainIter, ProteinChainRef, ProteinProjectionError, ProteinResidueIter,
    ProteinResidueRef, ProteinSelectionSummary, ResidueAddress, ResidueCode, ResidueCodeParseError,
    ResidueIdentity, ResidueInfo, ResidueInfoKind, ResidueKind, ResidueName, ResidueSequenceError,
    ResidueSourceIds, UNKNOWN_TABULATED_RESIDUE_INDEX, expand_one_letter,
    expand_one_letter_sequence, find_residue_info, find_residue_info_index, residue_code,
    residue_info, residue_info_checked,
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
#[cfg(feature = "cap-smiles")]
pub use cosmolkit_core::{
    AdjustHsError, AtropisomerError, CanonicalRankError, CleanupError, ConjugationError,
    DoubleBondStereoError, HybridizationError, LegacyStereoError, PathError, PeriodicTableError,
    PropertyCacheError, PropertyIntReadError, PropertyUIntReadError, PropertyULongReadError,
    RadicalError, RequiredPropertyStringError, RingDecomposerError, RingFindingError,
    StereoOrderError, UIntLexicalReadError,
};
#[cfg(feature = "cap-transforms")]
pub use cosmolkit_core::{
    AtomPathSubsetCopyError, AtomPositionParams, FullCopyComponentError, MoleculeFragmentsError,
    MoleculeFragmentsFailure, OrderedFragmentBuildError, SelectedBondCopyError,
    SelectedSubstanceGroupCopyError, SubstanceGroupAtomList, TransformError,
};
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
#[cfg(feature = "cap-smiles")]
pub use cosmolkit_core::{
    parser_helpers::ParserCarrierError, parser_stereo_order::ParserStereoOrderError,
};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::cif::{CifReadError, CifReadErrorKind};
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
    SmilesParseError, SmilesParseParams, SmilesStereoError, SmilesWriteParams,
};
#[cfg(feature = "cap-stereo")]
pub use cosmolkit_stereo::{CipLabelOptions, CipLabelerError, StereoReadError};
#[cfg(feature = "cap-matrices")]
pub use matrices::DistanceMatrixParams;
pub use molecule::{Molecule, MoleculeAtom};
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
#[cfg(cosmolkit_runtime_privacy_probe)]
pub(crate) use ops::LazyOutputProbeAccess;
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
#[cfg(feature = "cap-transforms")]
pub(crate) use ops::WithAtomPropertyAccess;
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
    SupportMatrixEntry, feature_spec, feature_specs, operation_invariant,
    operation_invariant_matrix, operation_parity, operation_spec, operation_specs, parity_matrix,
    support_matrix,
};
#[cfg(test)]
pub(crate) use ops::{
    CowCoordinatesFailureForTestAccess, CowCoordinatesForTestAccess,
    ReactionApplyReportForTestAccess, RingLiveCowCheckoutConflictForTestAccess,
};
#[cfg(feature = "cap-transforms")]
pub(crate) use ops::{FragmentsAccess, LargestFragmentAccess};
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
pub(crate) use ops::{WithMmffOptimizedAccess, WithMmffOptimizedConformersAccess};

#[cfg(feature = "cap-forcefields")]
pub use forcefields::{MmffEnergyGradient, MmffEvaluationParams};

#[cfg(feature = "cap-tautomer")]
mod tautomer;
#[cfg(feature = "cap-tautomer")]
pub use cosmolkit_tautomer::TautomerTransformError;
#[cfg(feature = "cap-tautomer")]
pub(crate) use ops::{CanonicalTautomerWithParamsAccess, EnumerateTautomersWithParamsAccess};
#[cfg(feature = "cap-tautomer")]
#[doc(hidden)]
pub use tautomer::canonical_tautomer_from_molecule_hosts_with_params;
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
#[cfg(all(feature = "cap-serialization", not(target_arch = "wasm32")))]
pub use cosmolkit_io::PickleError;
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{
    AtomSiteIdentityField, CoordinateBlockSelectionError, MissingAtomSiteIdentityField,
    MmcifAssemblyError, MmcifAtomSiteGroupingError, MmcifAtomSiteRowOutOfBounds,
    MmcifAtomSiteScalarError, MmcifCcdRestoreError, MmcifCisPepError, MmcifConnectionError,
    MmcifCrystalInfoError, MmcifEntityInfoError, MmcifHelixError, MmcifLabelAddressError,
    MmcifModResError, MmcifOperationExprError, MmcifSequenceIdError, MmcifSheetError,
    MmcifSiftsUnpError, MmcifTlsInfoError, PdbResidueKeyError,
};
#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
pub use cosmolkit_io::{BioMoleculeConversionError, BioMoleculeParams, PdbPostprocessError};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{
    ChemCompAtomError, ChemCompModelError, JsonError, JsonErrorKind, MmjsonReadError,
    MmjsonStructureError, MmjsonUnexpectedType, MmjsonValueError,
};
#[cfg(feature = "cap-bio")]
pub use cosmolkit_io::{
    PdbAtomFieldError, PdbChargeError, PdbCisPepRecordError, PdbConectError,
    PdbConnectionFinalizationError, PdbCrystalError, PdbDbRefError, PdbHeaderError, PdbHelixError,
    PdbInputFormatError, PdbLinkNameCounter, PdbLinkRecordError, PdbModResError,
    PdbModelTransitionError, PdbMtrixError, PdbReaderCause, PdbRemark3Error, PdbRemark200Error,
    PdbRemark350Error, PdbRemarkError, PdbRemarkFinalizationError, PdbRemarkMetadataError,
    PdbSheetError, PdbSsbondRecordError, PdbSsbondResolutionError, PdbStructureFinalizationError,
};
#[cfg(feature = "cap-fingerprints")]
pub use ops::AtomPairAtomCodeResult;
#[cfg(feature = "cap-fingerprints")]
pub(crate) use ops::WithAtomPairAtomCodeAccess;
#[cfg(feature = "cap-hashing")]
pub(crate) use ops::{MurckoDecomposeAccess, MurckoScaffoldAccess, NetScaffoldAccess};
#[cfg(feature = "cap-forcefields")]
pub use ops::{UffEnergyGradient, UffEvaluationParams};

#[cfg(feature = "cap-hashing")]
mod molecular_hash;
#[cfg(feature = "cap-hashing")]
pub use cosmolkit_core::{CipRankError, ScaffoldError};
#[cfg(feature = "cap-hashing")]
pub use cosmolkit_fingerprints::MoleculeHashError;

#[cfg(feature = "cap-transforms")]
pub use cosmolkit_core::{
    Coordinate2DInputParams, Coordinate3DInputParams, Coordinate3DReadError, CoordinateInputError,
    CoordinateZPolicy, Replace3DCoordinatesParams,
};

#[cfg(feature = "cap-transforms")]
pub(crate) use ops::{
    With2dCoordinateBlockAccess, With3dCoordinatesAccess, WithAdded3dConformerAccess,
    WithCleared3dConformersAccess, WithOnly3dConformerAccess,
};

#[cfg(feature = "cap-fingerprints")]
mod avalon_fingerprint;
#[cfg(feature = "cap-fingerprints")]
mod layered_fingerprint;
#[cfg(feature = "cap-fingerprints")]
pub use avalon_fingerprint::AvalonFingerprintError;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    AvalonError as AvalonEngineError, AvalonFingerprintFlags, AvalonFingerprintParams,
};
#[cfg(feature = "cap-fingerprints")]
mod topological_torsion_path;
#[cfg(feature = "cap-alignment")]
pub use cosmolkit_alignment::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AlignmentResult, AlignmentTransform,
    AllConformerRmsdParameters, BestAlignmentParameters, ConformerAlignmentParameters,
    ConformerAlignmentReport, ConformerRmsd, CoordinateRmsdParameters, ThreadSelectionError,
};
#[cfg(feature = "cap-alignment")]
pub use cosmolkit_core::ThreadCountError;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    LAYERED_FINGERPRINT_MAX_LAYERS, LAYERED_FINGERPRINT_SUBSTRUCTURE_LAYERS,
    LAYERED_FINGERPRINT_VERSION, LayeredFingerprintError, LayeredFingerprintLayers,
    LayeredFingerprintParams, LayeredFingerprintResult,
};
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{TopologicalTorsionPathScoreError, explain_path_score};
#[cfg(feature = "cap-fingerprints")]
pub use layered_fingerprint::{
    fingerprint_layered_query_with_output_with_params, fingerprint_layered_query_with_params,
};

#[cfg(feature = "cap-alignment")]
pub(crate) use ops::{WithAlignedConformersAccess, WithAlignmentToAccess};

#[cfg(feature = "cap-fingerprints")]
mod maccs_fingerprint;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_MaccsFingerprintParams;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{AtomCodeError, MaccsFingerprintError, MaccsFingerprintParams};

#[cfg(any(feature = "cap-conformer", feature = "cap-stereoisomers"))]
mod conformer_projection;

#[cfg(feature = "cap-conformer")]
mod conformer;
#[cfg(feature = "cap-conformer")]
pub use conformer::{
    ConformerError, ConformerRunError, EmbedFailureCause, EmbedFailureCountError, EmbedParams,
};
#[cfg(feature = "cap-conformer")]
pub use cosmolkit_conformer::{
    BoundKind, BoundsMatrixError, BoundsQueryError, EmbedSeedError, GenerationError,
    GenerationFailure, GraphBoundsError, MatrixAxis, PruningError,
};
#[cfg(feature = "cap-forcefields")]
pub use cosmolkit_forcefields::{
    ConformerOptimizerCause, ConformerOptimizerError, CrystalffTorsionPreferencesError,
    DistGeomForceFieldError, MissingExplicitHydrogensError, UffBoundsCause, UffBoundsError,
};

#[cfg(feature = "cap-conformer")]
pub use conformer::{EmbedMoleculeResult, EmbedMultipleConfsResult};

#[cfg(feature = "cap-conformer")]
pub(crate) use ops::{
    With3dConformerAccess, With3dConformerResultAccess, With3dConformersAccess,
    With3dConformersResultAccess,
};
#[cfg(feature = "cap-fingerprints")]
mod pattern_fingerprint;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    PATTERN_FINGERPRINT_VERSION, PatternFingerprintError, PatternFingerprintParams,
};
#[cfg(feature = "cap-fingerprints")]
pub use pattern_fingerprint::{fingerprint_pattern_query, fingerprint_pattern_query_with_params};

#[cfg(feature = "cap-fingerprints")]
mod topological_fingerprint;
#[cfg(feature = "cap-fingerprints")]
pub use cosmolkit_fingerprints::{
    TopologicalFingerprintError, TopologicalFingerprintOutput, TopologicalFingerprintOutputRequest,
    TopologicalFingerprintParams, TopologicalFingerprintResult,
};
#[cfg(feature = "cap-fingerprints")]
pub use topological_fingerprint::{
    fingerprint_topological_query_with_output_with_params,
    fingerprint_topological_query_with_params,
};

#[cfg(all(feature = "cap-batch", feature = "cap-fingerprints"))]
pub use cosmolkit_fingerprints::{
    BatchFingerprintAdditionalOutput, BatchFingerprintOutput, BatchFingerprintOutputError,
};

#[cfg(feature = "cap-io")]
mod molecular_io;
#[cfg(feature = "cap-io")]
pub use cosmolkit_io::SdfWriteError;
#[cfg(feature = "cap-io")]
mod property_text;
#[cfg(feature = "cap-io")]
mod sdf_supplier;
#[cfg(feature = "cap-io")]
pub use molecular_io::{
    Mol2PostError, Mol2ReadError, Mol2ReadParams, Mol2Type, MolBlockWriteParams,
    MolCoordinateSelection, MolWriteError, MolecularIoError, PdbWriteError, PdbWriteParams,
    SdfFormat, XyzReadError, XyzWriteError, XyzWriteParams,
};
#[cfg(feature = "cap-stereoisomers")]
pub use ops::StereoisomerIterator;
#[cfg(feature = "cap-io")]
pub use property_text::{PropertyStringError, property_value_to_text};
#[cfg(feature = "cap-io")]
pub use sdf_supplier::{
    SdfDataset, SdfDatasetIterator, SdfReader, SdfRecordMetadata, SdfRecordStream,
};

#[cfg(feature = "cap-stereoisomers")]
mod stereoisomers;
#[cfg(feature = "cap-stereoisomers")]
pub(crate) use ops::{EnumerateStereoisomersAccess, EnumerateStereoisomersWithRandomBitsAccess};
#[cfg(feature = "cap-stereoisomers")]
pub use stereoisomers::{
    EnumerationError, EnumerationRunError, StereoisomerOptions, StereoisomerRandomSource,
};

#[cfg(feature = "cap-forcefields")]
mod persistent_forcefields;
#[cfg(feature = "cap-forcefields")]
pub use persistent_forcefields::{
    ForceFieldAngleConstraintParams, ForceFieldDistanceConstraintParams, ForceFieldEnergyGradient,
    ForceFieldError, ForceFieldMinimizeOutcome, ForceFieldMinimizeParams,
    ForceFieldPositionConstraintParams, ForceFieldTorsionConstraintParams, MmffForceFieldError,
    MmffForceFieldParams, MolecularForceField, MolecularForceFieldErrorKind, UffForceFieldError,
    UffForceFieldParams,
};

#[cfg(feature = "cap-transforms")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_AtomPositionParams;
#[cfg(feature = "cap-transforms")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_Coordinate3DInputParams;
#[cfg(feature = "cap-matrices")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_DistanceMatrix3dParams;
#[cfg(feature = "cap-kekulize")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_KekulizeParams;
#[cfg(feature = "cap-stereo")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_PotentialStereoParams;
#[cfg(feature = "cap-hydrogens")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_RemoveHsParams;
#[cfg(feature = "cap-transforms")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_Replace3DCoordinatesParams;
#[cfg(feature = "cap-rings")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_RingSearchParams;
#[cfg(feature = "cap-stereo")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_StructureTagParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_PatternFingerprintParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_TopologicalFingerprintOutputRequest;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use cosmolkit_forcefields::__ck_configuration_MmffConformerOptimizationParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use cosmolkit_forcefields::__ck_configuration_MmffEvaluationParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use cosmolkit_forcefields::__ck_configuration_MmffOptimizationParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use cosmolkit_forcefields::__ck_configuration_MmffPropertiesParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use cosmolkit_forcefields::__ck_configuration_UffEvaluationParams;
#[cfg(all(feature = "cap-bio", feature = "cap-io"))]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_BioMoleculeParams;
#[cfg(feature = "cap-bio")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_BioPdbReadParams;
#[cfg(feature = "cap-bio")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_BioPdbWriteParams;
#[cfg(feature = "cap-io")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_XyzWriteParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionApplyParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionTemplateRemovalParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionValidationParams;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_McsAtomCompareParameters;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_McsBondCompareParameters;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_SmartsWriteParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use cosmolkit_smiles::__ck_configuration_RandomSmilesWriteParams;
#[cfg(feature = "cap-inchi")]
#[doc(hidden)]
pub use inchi::__configuration_InchiReadParams;
#[cfg(feature = "cap-inchi")]
#[doc(hidden)]
pub use inchi::__configuration_InchiWriteParams;
#[cfg(feature = "cap-matrices")]
#[doc(hidden)]
pub use matrices::__configuration_DistanceMatrixParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use ops::__configuration_UffConformerOptimizationParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use ops::__configuration_UffOptimizationParams;

#[cfg(all(feature = "cap-batch", feature = "cap-depict"))]
#[doc(hidden)]
pub use batch::__configuration_BatchImageParams;
#[cfg(feature = "cap-alignment")]
#[doc(hidden)]
pub use cosmolkit_alignment::__ck_configuration_AlignmentParameters;
#[cfg(feature = "cap-alignment")]
#[doc(hidden)]
pub use cosmolkit_alignment::__ck_configuration_AllConformerRmsdParameters;
#[cfg(feature = "cap-alignment")]
#[doc(hidden)]
pub use cosmolkit_alignment::__ck_configuration_BestAlignmentParameters;
#[cfg(feature = "cap-alignment")]
#[doc(hidden)]
pub use cosmolkit_alignment::__ck_configuration_ConformerAlignmentParameters;
#[cfg(feature = "cap-alignment")]
#[doc(hidden)]
pub use cosmolkit_alignment::__ck_configuration_CoordinateRmsdParameters;
#[cfg(feature = "cap-hydrogens")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_AddHsParams;
#[cfg(feature = "cap-depict")]
#[doc(hidden)]
pub use cosmolkit_depict::__ck_configuration_Compute2DCoordinatesParams as __ck_configuration_Coordinate2DParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_AtomPairGeneratorParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_LayeredFingerprintParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_MorganGeneratorParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_TopologicalFingerprintParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_TopologicalTorsionGeneratorParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionParseParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionRunParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionWriteParams;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_McsParameters;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_SmartsParseParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use cosmolkit_smiles::__ck_configuration_CxSmilesWriteParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use cosmolkit_smiles::__ck_configuration_SmilesParseParams;
#[cfg(feature = "cap-stereo")]
#[doc(hidden)]
pub use cosmolkit_stereo::__ck_configuration_CipLabelOptions;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use legacy_topological_torsion::__configuration_LegacyTopologicalTorsionParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use morgan_reusable::__configuration_MorganCallParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use smiles::__configuration_FragmentCxSmilesWriteParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use smiles::__configuration_FragmentSmilesWriteParams;
#[cfg(feature = "cap-tautomer")]
#[doc(hidden)]
pub use tautomer::__configuration_TautomerScoreParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use topological_torsion::__configuration_TopologicalTorsionCallParams;

#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use atom_pair::__configuration_AtomPairFingerprintParams;
#[cfg(feature = "cap-batch")]
#[doc(hidden)]
pub use batch::__configuration_BatchParams;
#[cfg(feature = "cap-batch")]
#[doc(hidden)]
pub use batch::__configuration_BatchQueryParams;
#[cfg(feature = "cap-batch")]
#[doc(hidden)]
pub use cosmolkit_batch::__ck_configuration_BatchExportParams;
#[cfg(feature = "cap-conformer")]
#[doc(hidden)]
pub use cosmolkit_conformer::__ck_configuration_EmbedParams;
#[cfg(feature = "cap-aromaticity")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_AromaticityParams;
#[cfg(feature = "cap-transforms")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_Coordinate2DInputParams;
#[cfg(feature = "cap-sanitize")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_SanitizeParams;
#[cfg(feature = "cap-valence")]
#[doc(hidden)]
pub use cosmolkit_core::__ck_configuration_ValenceParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_AtomPairAtomInvariantsGenerator;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use cosmolkit_fingerprints::__ck_configuration_AvalonFingerprintParams;
#[cfg(feature = "cap-bio")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_BioMmcifWriteParams;
#[cfg(feature = "cap-bio")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_BioReadParams;
#[cfg(feature = "cap-io")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_MolBlockWriteParams;
#[cfg(feature = "cap-io")]
#[doc(hidden)]
pub use cosmolkit_io::__ck_configuration_PdbExportParams as __ck_configuration_PdbWriteParams;
#[cfg(feature = "cap-reaction")]
#[doc(hidden)]
pub use cosmolkit_reaction::__ck_configuration_ReactionSingleRunParams;
#[cfg(feature = "cap-search")]
#[doc(hidden)]
pub use cosmolkit_search::__ck_configuration_SubstructMatchParams;
#[cfg(feature = "cap-smiles")]
#[doc(hidden)]
pub use cosmolkit_smiles::__ck_configuration_SmilesWriteParams;
#[cfg(feature = "cap-stereoisomers")]
#[doc(hidden)]
pub use cosmolkit_stereo::__ck_configuration_StereoisomerOptions;
#[cfg(feature = "cap-io")]
#[doc(hidden)]
pub use molecular_io::__configuration_Mol2ReadParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use morgan::__configuration_MorganFingerprintParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use morgan_reusable::__configuration_MorganBondInvariantsGenerator;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_ForceFieldAngleConstraintParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_ForceFieldDistanceConstraintParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_ForceFieldMinimizeParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_ForceFieldPositionConstraintParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_ForceFieldTorsionConstraintParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_MmffForceFieldParams;
#[cfg(feature = "cap-forcefields")]
#[doc(hidden)]
pub use persistent_forcefields::__configuration_UffForceFieldParams;
#[cfg(feature = "cap-io")]
#[doc(hidden)]
pub use sdf::__configuration_SdfReadParams;
#[cfg(feature = "cap-tautomer")]
#[doc(hidden)]
pub use tautomer::__configuration_TautomerParams;
#[cfg(feature = "cap-fingerprints")]
#[doc(hidden)]
pub use topological_torsion::__configuration_TopologicalTorsionFingerprintParams;
