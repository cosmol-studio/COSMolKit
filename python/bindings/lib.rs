mod alignment_calls;
mod alignment_errors;
mod alignment_queries;
mod alignment_values;
mod allocation_errors;
#[path = "../framework/arguments.rs"]
mod arguments;
mod atom_bond_values;
#[cfg(test)]
mod atom_property_calls;
mod attachment_values;
mod batch_calls;
mod batch_errors;
mod batch_fingerprint_values;
mod batch_query_params;
mod binary_archive;
mod bio_errors;
mod bio_geometry_values;
mod bio_hierarchy;
mod bio_metadata_values;
mod bio_mmcif_reader;
mod bio_mmcif_writer;
mod bio_molecule_conversion;
mod bio_pdb_reader;
mod bio_pdb_writer;
mod bio_read_errors;
mod bio_reader_calls;
mod bio_relationships;
mod bio_rows;
mod bio_scalar_metadata;
mod bio_selection_errors;
mod bio_selection_values;
mod bio_source_values;
mod bio_structure;
mod bio_values;
mod bio_write_configuration;
mod builder;
#[path = "../framework/callbacks.rs"]
mod callbacks;
mod canonical_errors;
mod chemical_params;
mod chemical_values;
mod chemistry_calls;
mod chemistry_problems;
mod cip_atom_code_errors;
mod collection_params;
#[path = "../framework/configuration.rs"]
mod configuration;
#[cfg(test)]
pub use configuration::fixture::{
    __configuration_ProjectionChild, __configuration_ProjectionExample,
};
mod conformer_execution;
mod conformer_queries;
mod construction_values;
mod coordinate_2d_errors;
mod coordinate_calls;
mod coordinate_input;
mod coordinate_input_errors;
mod coordinates;
mod descriptor_errors;
mod detached_values;
mod double_stereo_errors;
mod drawing_errors;
mod elements;
mod embedding_calls;
mod embedding_values;
#[path = "../framework/enums.rs"]
mod enums;
#[path = "../framework/errors.rs"]
mod error_framework;
mod errors;
mod feature_metadata;
mod fingerprint_batch;
mod fingerprint_calls;
mod fingerprint_configuration;
mod fingerprint_explanations;
mod fingerprint_json;
mod fingerprint_providers;
mod fingerprint_values;
mod forcefield_calls;
mod forcefield_execution_errors;
mod fragment_errors;
mod group_values;
mod hydrogen_errors;
mod inchi_values;
#[path = "../framework/indexing.rs"]
mod indexing;
mod io_batch_params;
mod io_source_errors;
mod io_write_errors;
mod kekulize_errors;
mod layered_fingerprint_errors;
mod layered_producers;
mod leaf_values;
mod legacy_torsion_calls;
#[cfg(test)]
mod maccs_entrypoints;
mod maccs_errors;
mod match_errors;
mod matrix_errors;
mod mcs;
mod mcs_errors;
mod mcs_producers;
mod metadata_values;
mod mmff_properties;
mod model_edit_errors;
mod model_errors;
mod mol2_post_errors;
mod mol2_read_errors;
mod mol_post_errors;
mod mol_post_support_errors;
mod molecular_hash;
mod molecular_io_calls;
mod molecular_io_errors;
mod morgan_configuration;
mod morgan_errors;
mod morgan_reusable;
mod native_numerical_values;
mod notation_params;
mod numerical_params;
mod numerical_results;
mod operation_errors;
mod operation_metadata;
mod optional_configuration;
mod ordered_results;
#[path = "../framework/parameters.rs"]
mod parameters;
mod params;
mod parser_errors;
mod path_errors;
mod path_score_errors;
mod pattern_fingerprint_errors;
mod pattern_producers;
mod pdb_flags;
mod persistent_forcefields;
mod potential_stereo_errors;
mod property_numeric_errors;
#[cfg(test)]
mod property_string_queries;
mod property_values;
mod protein_values;
mod public_selectors;
mod query_compiler;
mod query_model_errors;
mod query_payloads;
mod rdkit_adapter;
mod reaction;
mod reaction_calls;
mod reaction_execution_errors;
mod reaction_parse;
mod reaction_validation;
mod reaction_values;
mod reaction_writer;
mod residue_values;
#[path = "../framework/retained_error.rs"]
mod retained_error;
mod ring_errors;
mod ring_results;
mod sanitize_errors;
mod sanitize_leaf_errors;
mod sanitize_support_errors;
mod sanitize_values;
mod scaffold_errors;
mod sdf_read_errors;
mod sdf_records;
mod sdf_sources;
mod sdf_values;
mod search_callbacks;
mod search_reaction_params;
mod selectors;
#[path = "../framework/sequence.rs"]
mod sequence;
mod smarts_errors;
mod smarts_factory;
mod smarts_writer;
mod smarts_writer_support_errors;
mod smiles_parse_errors;
mod smiles_stereo_errors;
mod smiles_writer;
mod stereo_errors;
mod stereo_groups;
mod stereo_order_errors;
mod stereo_read_errors;
mod stereo_values;
mod stereoisomer_iterator;
mod stereoisomer_options;
mod string_errors;
mod tautomer_callbacks;
mod tautomer_catalog;
mod tautomer_configuration;
mod tautomer_errors;
mod tautomer_producers;
mod tautomer_transform_errors;
mod tautomer_values;
mod topological_fingerprint_errors;
mod topological_producers;
mod topology_edit_errors;
mod topology_values;
mod torsion_reusable;
mod transform_errors;
mod uff_parameter_errors;
mod utf8_errors;
mod valence_errors;
mod xyz_errors;

#[pyo3::pymodule(name = "_native")]
pub mod native {
    #[pymodule_init]
    fn register_errors(module: &pyo3::Bound<'_, pyo3::types::PyModule>) -> pyo3::PyResult<()> {
        super::error_framework::register(module)
    }
    #[pymodule_export]
    use super::batch_fingerprint_values::{
        BatchFingerprintAdditionalOutput, BatchFingerprintOutput, BatchFingerprintOutputError,
        BatchFingerprintOutputErrorKind, TopologicalFingerprintParams,
    };
    #[pymodule_export]
    use super::batch_query_params::BatchQueryParams;
    #[pymodule_export]
    use super::binary_archive::{
        _molecule_from_binary, PickleError, PickleErrorKind, mol_from_binary,
    };
    #[pymodule_export]
    use super::ring_results::RingSearchResult;

    #[pymodule_export]
    use super::bio_selection_errors::{
        BioOperationError, BioOperationErrorKind, BioRowChainError, BioRowChainErrorKind,
        BioRowModelError, BioRowModelErrorKind, BioRowTraverseError, BioRowTraverseErrorKind,
        BioSelectionCopyCause, BioSelectionCopyCauseKind, BioSelectionCopyError,
        BioSelectionMatchError,
    };
    #[pymodule_export]
    use super::bio_selection_values::{
        BioSelection, BioSelectionParseError, BioSelectionParseErrorKind,
    };
    #[pymodule_export]
    use super::builder::MoleculeBuilder;
    #[pymodule_export]
    use super::canonical_errors::{
        CanonicalRankError, CanonicalRankErrorKind, PropertyStringError, PropertyStringErrorKind,
        property_value_to_text,
    };
    #[pymodule_export]
    use super::canonical_errors::{RequiredPropertyStringError, RequiredPropertyStringErrorKind};
    #[pymodule_export]
    use super::chemistry_problems::{
        ChemistryProblem, ChemistryProblemError, ChemistryProblemErrorKind, ChemistryProblemReport,
    };
    #[pymodule_export]
    use super::cip_atom_code_errors::{AtomCodeError, AtomCodeErrorKind};
    #[pymodule_export]
    use super::cip_atom_code_errors::{CipLabelerError, CipLabelerErrorKind};
    #[pymodule_export]
    use super::conformer_execution::{
        BoundKind, BoundsMatrixError, BoundsMatrixErrorKind, BoundsQueryError,
        BoundsQueryErrorKind, ConformerOptimizerCause, ConformerOptimizerCauseKind,
        ConformerOptimizerError, ConformerRunError, CrystalffTorsionPreferencesError,
        CrystalffTorsionPreferencesErrorKind, DistGeomForceFieldError, DistGeomForceFieldErrorKind,
        EmbedSeedError, EmbedSeedErrorKind, EnumerationError, EnumerationErrorKind,
        EnumerationRunError, GenerationError, GenerationErrorKind, GenerationFailure,
        GenerationFailureKind, GraphBoundsError, GraphBoundsErrorKind, MatrixAxis,
        MissingExplicitHydrogensError, PoisonError, PruningError, PruningErrorKind, UffBoundsCause,
        UffBoundsCauseKind, UffBoundsError,
    };
    #[pymodule_export]
    use super::coordinate_input_errors::{CoordinateInputError, CoordinateInputErrorKind};
    #[pymodule_export]
    use super::descriptor_errors::{
        DescriptorError, DescriptorErrorKind, DescriptorReadError, DescriptorReadErrorKind,
        DescriptorSearchCause, DescriptorSearchCauseKind,
    };
    #[pymodule_export]
    use super::detached_values::{CoordinateBlock, SdfRecordMetadata};
    #[pymodule_export]
    use super::feature_metadata::{FeatureSpec, FeatureSpecIter, feature_spec, feature_specs};
    #[pymodule_export]
    use super::forcefield_calls::{UffParameterQueryError, UffParameterQueryErrorKind};
    #[pymodule_export]
    use super::forcefield_execution_errors::{
        AutomaticForceFieldConstructionError, AutomaticForceFieldConstructionErrorKind,
        DispatchedUffOptimizationError, DispatchedUffOptimizationErrorKind,
        ForceFieldConstructionError, ForceFieldConstructionErrorKind, FragmentCoordinateViewError,
        FragmentCoordinateViewErrorKind, MmffBuilderError, MmffBuilderErrorKind,
        MmffOptimizationError, MmffOptimizationFailure, MmffOptimizationFailureKind,
        MmffOptimizationSourceError, NonbondedAssemblyError, NonbondedAssemblyErrorKind,
        NonbondedFragmentMappingError, NonbondedFragmentMappingErrorKind, OptimizationStageError,
        OptimizationStageErrorKind, SerialConformerOptimizationError,
        SerialConformerOptimizationErrorKind, SerialCoordinateCountError,
        SerialResultCapacityError, SerialUffOptimizationError, SerialUffOptimizationErrorKind,
        SingleConformerOptimizationError, SingleConformerOptimizationErrorKind, UffConformerError,
        UffEvaluationError, UffEvaluationFailure, UffEvaluationFailureKind, UffOptimizationError,
        UffPreparedOptimizationError, UffPreparedOptimizationErrorKind, UffSingleError,
        UffThreadCountError, UffThreadCountErrorKind,
    };
    #[pymodule_export]
    use super::fragment_errors::{
        AtomPathSubsetCopyError, AtomPathSubsetCopyErrorKind, FullCopyComponentError,
        FullCopyComponentErrorKind, MoleculeFragmentsError, MoleculeFragmentsFailure,
        MoleculeFragmentsFailureKind, OrderedFragmentBuildError, OrderedFragmentBuildErrorKind,
        SelectedBondCopyError, SelectedBondCopyErrorKind, SelectedSubstanceGroupCopyError,
        SelectedSubstanceGroupCopyErrorKind, SubstanceGroupAtomList,
    };
    #[pymodule_export]
    use super::hydrogen_errors::{HydrogenError, HydrogenErrorKind};
    #[pymodule_export]
    use super::inchi_values::{
        _mol_from_inchi_keywords, _mol_from_inchi_params, InchiDiagnostic, InchiDiagnosticLevel,
        InchiError, InchiErrorKind, InchiResult, inchi_to_key, mol_from_inchi,
    };
    #[pymodule_export]
    use super::kekulize_errors::{KekulizeError, KekulizeErrorKind};
    #[pymodule_export]
    use super::layered_fingerprint_errors::{LayeredFingerprintError, LayeredFingerprintErrorKind};
    #[pymodule_export]
    use super::layered_producers::{
        fingerprint_layered_query_with_output_with_params, fingerprint_layered_query_with_params,
    };
    #[pymodule_export]
    use super::leaf_values::{
        MmffAtomProperties, MmffVariant, OptimizationStatus, TautomerEnumerationStatus,
        TopologicalFingerprintOutput, UffOptimizationErrorKind, UffParameterErrorKind,
    };
    #[pymodule_export]
    use super::match_errors::{
        MatchError, MatchErrorKind, QueryMatchContextError, QueryMatchContextErrorKind,
        SubstructMatchError, SubstructMatchErrorKind,
    };
    #[pymodule_export]
    use super::matrix_errors::{MatrixError, MatrixErrorKind};
    #[pymodule_export]
    use super::mcs::McsResult;
    #[pymodule_export]
    use super::mcs_errors::{
        McsError, McsErrorKind, McsProgressError, McsProgressErrorKind, SearchMcsError,
        SearchMcsErrorKind,
    };
    #[pymodule_export]
    use super::mcs_producers::{
        _maximum_common_substructure_keywords, _maximum_common_substructure_params,
        maximum_common_substructure, maximum_common_substructure_with_params,
    };
    #[pymodule_export]
    use super::mmff_properties::{
        MmffAngleLookupError, MmffAngleLookupErrorKind, MmffMolPropertiesError,
        MmffMolPropertiesErrorKind, MmffOopLookupError, MmffOopLookupErrorKind, MmffParamError,
        MmffParamErrorKind, MmffParamParseCause, MmffParamParseCauseKind, MmffParamParseError,
        MmffParamTable, MmffProperties, MmffTorLookupError, MmffTorLookupErrorKind,
        MmffVdwParseError, MmffVdwParseErrorKind,
    };
    #[pymodule_export]
    use super::model_edit_errors::{
        AdjacencyError, AdjacencyErrorKind, MappingValidationError, MappingValidationErrorKind,
        QueryAtomConversionError, QueryAtomConversionErrorKind, QueryStateError,
        QueryStateErrorKind,
    };
    #[pymodule_export]
    use super::model_edit_errors::{MoleculePropertyError, MoleculePropertyErrorKind};
    #[pymodule_export]
    use super::model_errors::{
        CoordinateValidationError, CoordinateValidationErrorKind, StereoGroupError,
        StereoGroupErrorKind, TopologyValidationError, TopologyValidationErrorKind,
    };
    #[pymodule_export]
    use super::molecular_hash::{
        CipRankError, CipRankErrorKind, MoleculeHashError, MoleculeHashErrorKind,
    };
    #[pymodule_export]
    use super::morgan_configuration::{
        MorganAtomInvariantsGenerator, MorganFingerprintParams, MorganInvariants,
    };
    #[pymodule_export]
    use super::native_numerical_values::{
        DenseMatrix, ForceFieldEnergyGradient, ForceFieldMinimizeOutcome,
    };
    #[pymodule_export]
    use super::numerical_results::{
        CrippenTotals, LabuteAsaContributions, MmffEnergyGradient, MmffOptimizeMoleculeConfResult,
        MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult, OptimizationOutcome,
        UffConformerOptimizationResult, UffConformerResult, UffEnergyGradient,
        UffOptimizationResult,
    };
    #[pymodule_export]
    use super::operation_errors::{
        MoleculeOpOutput, OperationCipStatePolicy, OperationDerivedState, OperationError,
        OperationErrorKind, OperationMappingRequirement, OperationSemanticPreconditions,
        UnsupportedFeatureError,
    };
    #[pymodule_export]
    use super::operation_metadata::{
        FunctionStatus, MoleculeOpSpec, OperationInvariantEntry, ParityMatrixEntry,
        SupportMatrixEntry, operation_invariant, operation_invariant_matrix, operation_parity,
        operation_spec, operation_specs, parity_matrix, support_matrix,
    };
    #[pymodule_export]
    use super::ordered_results::{
        AtomPairAtomCodeResult, LayeredFingerprintResult, LigandRef, MatchResult, TetrahedralStereo,
    };
    #[pymodule_export]
    use super::parser_errors::{ParserCarrierError, ParserCarrierErrorKind};
    #[pymodule_export]
    use super::parser_errors::{ParserStereoOrderError, ParserStereoOrderErrorKind};
    #[pymodule_export]
    use super::path_errors::{
        PathError, PathErrorKind, PeriodicTableError, PeriodicTableErrorKind,
    };
    #[pymodule_export]
    use super::path_score_errors::{
        TopologicalTorsionPathScoreError, TopologicalTorsionPathScoreErrorKind,
    };
    #[pymodule_export]
    use super::pattern_fingerprint_errors::{PatternFingerprintError, PatternFingerprintErrorKind};
    #[pymodule_export]
    use super::pattern_producers::{
        _fingerprint_pattern_query_keywords, _fingerprint_pattern_query_params,
        fingerprint_pattern_query, fingerprint_pattern_query_with_params,
    };
    #[pymodule_export]
    use super::persistent_forcefields::{
        ForceFieldError, MmffForceFieldError, MolecularForceField, UffForceFieldError,
    };
    #[pymodule_export]
    use super::property_numeric_errors::{PropertyULongReadError, PropertyULongReadErrorKind};
    #[pymodule_export]
    use super::public_selectors::{
        AvalonFingerprintFlags, BioMmcifReadStage, BioPdbReadStage, LayeredFingerprintLayers,
        MolecularForceFieldErrorKind, ParityPolicy, RotatableBondsOptions,
    };
    #[pymodule_export]
    use super::query_compiler::{CompiledQuery, QueryGraph, compile_query};
    #[pymodule_export]
    use super::query_model_errors::{QueryGraphError, QueryGraphErrorKind};
    #[pymodule_export]
    use super::reaction::{
        Reaction, ReactionModelError, ReactionModelErrorKind, ReactionTemplateRemoval,
    };
    #[pymodule_export]
    use super::reaction_validation::{
        ReactionInitializationError, ReactionInitializationErrorKind, ReactionValidationError,
        ReactionValidationErrorKind,
    };
    #[pymodule_export]
    use super::reaction_writer::{ReactionWriteError, ReactionWriteErrorKind};
    #[pymodule_export]
    use super::scaffold_errors::{ScaffoldError, ScaffoldErrorKind};
    #[pymodule_export]
    use super::sdf_records::{
        _mol_from_mol_keywords, _mol_from_mol_params, _mol_from_sdf_keywords, _mol_from_sdf_params,
        SdfError, SdfErrorKind, SdfRecord, mol_from_mol, mol_from_sdf,
    };
    #[pymodule_export]
    use super::stereo_groups::{StereoGroup, StereoGroupKind};
    #[pymodule_export]
    use super::stereoisomer_iterator::StereoisomerIterator;
    #[pymodule_export]
    use super::tautomer_producers::{
        canonical_tautomer_from_molecules, canonical_tautomer_from_molecules_with_params,
    };
    #[pymodule_export]
    use super::topological_fingerprint_errors::{
        TopologicalFingerprintError, TopologicalFingerprintErrorKind,
    };
    #[pymodule_export]
    use super::topological_producers::{
        TopologicalFingerprintResult, fingerprint_topological_query_with_output_with_params,
        fingerprint_topological_query_with_params,
    };
    #[pymodule_export]
    use super::uff_parameter_errors::{
        ForceFieldAngleIndexArgument, ForceFieldAngleRangeBound, ForceFieldBondIndexArgument,
        ForceFieldIndexArgument, ForceFieldKernelError, ForceFieldKernelErrorKind,
        ForceFieldTorsionIndexArgument, UffBondMathError, UffBondMathErrorKind, UffBuilderError,
        UffBuilderErrorKind, UffDefaultTorsionQueryError, UffDefaultTorsionQueryErrorKind,
        UffInversionContributionError, UffInversionContributionErrorKind,
        UffInversionIndexArgument, UffParamError, UffParamErrorKind, UffParameterError,
        UffPreparedValenceField, UffTorsionBondQueryError, UffTorsionBondQueryErrorKind,
        UffTypingError, UffTypingErrorKind, UffTypingInput,
    };
    use crate::parameters::{ConfigurationField, FieldInput, Scalar};

    #[pymodule_export]
    use super::allocation_errors::{
        TryReserveError, UnsignedStreamArrayError, UnsignedStreamArrayErrorKind,
    };
    #[pymodule_export]
    use super::io_write_errors::{
        PdbWriteError, PdbWriteErrorKind, XyzWriteError, XyzWriteErrorKind,
    };
    #[pymodule_export]
    use super::mol_post_errors::{
        MolPostError, MolPostErrorKind, MolProcessingError, MolProcessingErrorKind,
    };
    #[pymodule_export]
    use super::mol_post_support_errors::{
        AttachmentExpansionError, AttachmentExpansionErrorKind, BondDirectionStereoError,
        BondDirectionStereoErrorKind,
    };
    #[pymodule_export]
    use super::mol2_read_errors::{Mol2ReadError, Mol2ReadErrorKind};
    #[pymodule_export]
    use super::morgan_errors::{
        AtomPairError, AtomPairErrorKind, AtomPairReadError, AtomPairReadErrorKind,
        FingerprintPreparationError, FingerprintPreparationErrorKind, FingerprintWorkerError,
        FingerprintWorkerErrorKind, MorganError, MorganErrorKind, MorganReadError,
        MorganReadErrorKind, TopologicalTorsionError, TopologicalTorsionErrorKind,
        TopologicalTorsionReadError, TopologicalTorsionReadErrorKind,
    };
    #[pymodule_export]
    use super::morgan_reusable::{MorganFingerprintGenerator, MorganSettings};
    #[pymodule_export]
    use super::query_payloads::{
        AtomQueryPredicate, AtomQueryPredicateKind, AtomRangeBounds, AtomRangeBoundsKind,
        AtomRangeDataFunction, AtomRangeDataFunctionKind, AtomRangeQuery, BondQueryPredicate,
        BondQueryPredicateKind, RecursiveStructureQuery,
    };
    #[pymodule_export]
    use super::reaction_execution_errors::{
        ReactionApplyError, ReactionApplyErrorKind, ReactionProductError, ReactionProductErrorKind,
        ReactionRunError, ReactionRunErrorKind,
    };
    #[pymodule_export]
    use super::reaction_parse::{
        _parse_smirks_keywords, _parse_smirks_params, _reaction_from_smirks_keywords,
        _reaction_from_smirks_params,
    };
    #[pymodule_export]
    use super::reaction_parse::{
        CxParseError, ReactionParseError, ReactionParseErrorKind, parse_smirks,
        parse_smirks_with_params, reaction_from_smirks,
    };
    #[pymodule_export]
    use super::reaction_values::{
        ReactionApplyResult, ReactionRole, ReactionValidationIssue, ReactionValidationIssueKind,
        ReactionValidationReport, ReactionValidationSeverity,
    };
    #[pymodule_export]
    use super::ring_errors::{
        RingDecomposerError, RingDecomposerErrorKind, RingFindingError, RingFindingErrorKind,
    };
    #[pymodule_export]
    use super::sdf_read_errors::{SdfReadError, SdfReadErrorKind};
    #[pymodule_export]
    use super::sdf_values::SdfGraph;
    use super::search_callbacks::{AtomMatch, BondMatch, FinalMatch};
    #[pymodule_export]
    use super::search_callbacks::{QueryAtom, SubstructMatchParams};
    #[pymodule_export]
    use super::smarts_errors::{
        CxQueryLoweringError, CxQueryLoweringErrorKind, CxStereoGroupMergeError,
        CxStereoGroupMergeErrorKind, SmartsParseError, SmartsParseErrorKind,
    };
    #[pymodule_export]
    use super::smarts_factory::{
        _parse_smarts_keywords, _parse_smarts_params, parse_smarts, parse_smarts_with_params,
    };
    #[pymodule_export]
    use super::smarts_writer::{
        _write_cx_smarts_keywords, _write_cx_smarts_params, _write_smarts_keywords,
        _write_smarts_params, SmartsWriteError, SmartsWriteErrorKind, write_cx_smarts,
        write_smarts,
    };
    #[pymodule_export]
    use super::smarts_writer_support_errors::{
        SmartsTraversalError, SmartsTraversalErrorKind, WedgeError, WedgeErrorKind,
    };
    #[pymodule_export]
    use super::smiles_parse_errors::{SmilesParseError, SmilesParseErrorKind};
    #[pymodule_export]
    use super::smiles_stereo_errors::{SmilesStereoError, SmilesStereoErrorKind};
    #[pymodule_export]
    use super::smiles_writer::{
        FragmentWriteInputError, FragmentWriteInputErrorKind, SmilesWriteError,
        SmilesWriteErrorKind,
    };
    #[pymodule_export]
    use super::stereo_errors::{LegacyStereoError, LegacyStereoErrorKind};
    #[pymodule_export]
    use super::string_errors::{
        AvalonEngineError, AvalonEngineErrorKind, QueryCompileError, QueryCompileErrorKind,
    };
    #[pymodule_export]
    use super::tautomer_errors::{TautomerRunError, TautomerRunErrorKind};
    #[pymodule_export]
    use super::topology_edit_errors::{
        BondEndPointsParseErrorKind, TopologyEditError, TopologyEditErrorKind,
    };
    #[pymodule_export]
    use super::topology_values::TopologyBlock;
    #[pymodule_export]
    use super::torsion_reusable::{
        TopologicalTorsionFingerprintGenerator, TopologicalTorsionSettings,
    };
    #[pymodule_export]
    use super::transform_errors::{TransformError, TransformErrorKind};
    #[pymodule_export]
    use super::valence_errors::{ValenceError, ValenceErrorKind, ValencePhase};
    #[pymodule_export]
    use super::xyz_errors::{ParseFloatError, XyzReadError, XyzReadErrorKind};
    use super::{arguments::Argument, coordinates, errors};
    use cosmolkit as ck;
    use numpy::IntoPyArray;
    use pyo3::prelude::*;

    #[pymodule_export]
    use super::double_stereo_errors::{DoubleBondStereoError, DoubleBondStereoErrorKind};
    #[pymodule_export]
    use super::embedding_values::{
        ConformerError, ConformerErrorKind, CoordinateDimension, EmbedFailureCause,
        EmbedFailureCountError, EmbedFailureCountErrorKind, EmbedMoleculeResult,
        EmbedMultipleConfsResult, EmbedParams, EmbedTermination,
    };
    #[pymodule_export]
    use super::io_source_errors::{ExternalError, IoError, IoErrorKind};
    #[pymodule_export]
    use super::model_edit_errors::{
        AtomPropertyError, AtomPropertyErrorKind, BondValueError, BondValueErrorKind,
    };
    #[pymodule_export]
    use super::pdb_flags::{PdbReadFlags, PdbWriteFlags};
    #[pymodule_export]
    use super::potential_stereo_errors::{PotentialStereoError, PotentialStereoErrorKind};
    #[pymodule_export]
    use super::property_numeric_errors::{
        PropertyIntReadError, PropertyIntReadErrorKind, PropertyUIntReadError,
        PropertyUIntReadErrorKind, UIntLexicalReadError, UIntLexicalReadErrorKind,
    };
    #[pymodule_export]
    use super::sanitize_errors::{SanitizeError, SanitizeErrorKind};
    #[pymodule_export]
    use super::sanitize_leaf_errors::{
        AdjustHsError, AdjustHsErrorKind, HybridizationError, HybridizationErrorKind, RadicalError,
        RadicalErrorKind,
    };
    #[pymodule_export]
    use super::sanitize_support_errors::{
        AromaticityError, AromaticityErrorKind, CleanupError, CleanupErrorKind, ConjugationError,
        ConjugationErrorKind, PropertyCacheError, PropertyCacheErrorKind,
    };
    #[pymodule_export]
    use super::sanitize_values::{SanitizeOperations, SanitizeParams, SanitizeStage};
    #[pymodule_export]
    use super::stereo_errors::{
        AtropisomerError, AtropisomerErrorKind, StereoError, StereoErrorKind,
    };
    #[pymodule_export]
    use super::stereo_order_errors::{StereoOrderError, StereoOrderErrorKind};
    #[pymodule_export]
    use super::stereo_read_errors::{StereoReadError, StereoReadErrorKind};
    #[pymodule_export]
    use super::stereo_values::{
        PotentialStereoCenter, PotentialStereoDescriptor, PotentialStereoInfo,
        PotentialStereoResult, PotentialStereoSpecified, PotentialStereoType, RingStereoRelation,
    };
    #[pymodule_export]
    use super::stereoisomer_options::{StereoisomerOptions, StereoisomerRandomSource};
    #[pymodule_export]
    use super::tautomer_catalog::{TautomerCatalogError, TautomerCatalogErrorKind};
    #[pymodule_export]
    use super::tautomer_configuration::TautomerParams;
    #[pymodule_export]
    use super::tautomer_transform_errors::{TautomerTransformError, TautomerTransformErrorKind};
    #[pymodule_export]
    use super::utf8_errors::{FromUtf8Error, Utf8Error};

    #[pymodule_export]
    use super::io_batch_params::{
        BatchErrorMode, BatchExportParams, BatchImageParams, BatchParams, Mol2ReadParams, Mol2Type,
        MolBlockWriteParams, MolCoordinateSelection, PdbWriteParams, SdfFormat,
    };

    #[pymodule_export]
    use super::search_reaction_params::{
        McsParameters, ReactionSingleRunParams, ReactionWriteParams,
    };

    #[pymodule_export]
    use super::notation_params::{
        CxCoordinateSelection, CxSmilesFields, CxSmilesWriteParams, FragmentCxSmilesWriteParams,
        FragmentSmilesWriteParams, SmartsParseParams,
    };

    #[pymodule_export]
    use super::attachment_values::{
        AtomMetadata, TemplateAttachment, TemplateAttachmentOrder, TemplateAttachmentOrderError,
        TemplateAttachmentOrderErrorKind,
    };

    #[pymodule_export]
    use super::group_values::{
        SGroupBracket, SGroupCState, SGroupDisplay, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind,
    };

    #[pymodule_export]
    use super::atom_bond_values::{
        Atom, Bond, CipDescriptor, CipDescriptorError, CipDescriptorErrorKind, MoleculeAtom,
    };

    #[pymodule_export]
    use super::alignment_errors::{
        AlignmentError, AlignmentErrorKind, ThreadCountError, ThreadCountErrorKind,
        ThreadSelectionError,
    };
    #[pymodule_export]
    use super::alignment_values::{
        AlignmentResult, AlignmentTransform, AllConformerRmsdParameters, BestAlignmentParameters,
        ConformerAlignmentReport, ConformerRmsd, CoordinateRmsdParameters,
    };

    #[pymodule_export]
    use super::bio_hierarchy::{BioAtom, BioChain, BioModel, BioResidue};
    #[pymodule_export]
    use super::bio_structure::{BioStructure, BioStructureParts};

    #[pymodule_export]
    use super::protein_values::{
        Protein, ProteinAtomRef, ProteinChainRef, ProteinProjectionError,
        ProteinProjectionErrorKind, ProteinResidueRef, ProteinSelectionSummary,
    };

    #[pymodule_export]
    use super::bio_relationships::{
        AtomAddress, BioAsu, BioCisPep, BioConnection, BioConnectionKind, BioHelix, BioHelixClass,
        BioModRes, BioSheet, BioStrand, ResidueAddress,
    };

    #[pymodule_export]
    use super::bio_metadata_values::{
        BioExperimentInfo, BioExperimentalCrystalInfo, BioMetadata, BioRefinementInfo,
        BioSoftwareClassification, BioSoftwareItem, BioTlsGroup, BioTlsSelection,
    };

    #[pymodule_export]
    use super::bio_scalar_metadata::{
        BioBasicRefinementInfo, BioDiffractionInfo, BioRefinementRestraint, BioReflectionsInfo,
        BioStructureSourceState,
    };

    #[pymodule_export]
    use super::bio_geometry_values::{
        BioAssembly, BioAssemblyGenerator, BioAssemblyOperator, BioAssemblySpecialKind,
        BioCoordinateBlock, BioCrystalInfo, BioNcsOperator,
    };

    #[pymodule_export]
    use super::bio_rows::{
        BioAtomRow, BioCalcFlag, BioEntityDbRef, BioEntityRow, BioResidueRow, BioSiftsUnpResidue,
    };

    #[pymodule_export]
    use super::bio_source_values::{
        AltLocLabel, AltLocRequest, AtomSourceIds, EntityKind, EntitySourceIds, PdbSeqId,
        PolymerKind, ResidueKind, ResidueSourceIds,
    };

    #[pymodule_export]
    use super::residue_values::{
        ResidueCode, ResidueCodeParseError, ResidueIdentity, ResidueInfo, ResidueInfoKind,
        ResidueSequenceError, ResidueSequenceErrorKind, expand_one_letter,
        expand_one_letter_sequence, find_residue_info, find_residue_info_index, residue_code,
        residue_info, residue_info_checked,
    };

    #[pymodule_export]
    use super::metadata_values::{MoleculeProperties, SdfPropertyList, SdfPropertyListTarget};

    #[pymodule_export]
    use super::bio_mmcif_reader::{
        AtomSiteIdentityField, BioMmcifReadError, CoordinateBlockSelectionError,
        CoordinateBlockSelectionErrorKind, MissingAtomSiteIdentityField, MmcifAssemblyError,
        MmcifAssemblyErrorKind, MmcifAtomSiteGroupingError, MmcifAtomSiteGroupingErrorKind,
        MmcifAtomSiteRowOutOfBounds, MmcifAtomSiteScalarError, MmcifAtomSiteScalarErrorKind,
        MmcifCcdRestoreError, MmcifCcdRestoreErrorKind, MmcifCisPepError, MmcifCisPepErrorKind,
        MmcifConnectionError, MmcifConnectionErrorKind, MmcifCrystalInfoError,
        MmcifCrystalInfoErrorKind, MmcifEntityInfoError, MmcifEntityInfoErrorKind, MmcifHelixError,
        MmcifHelixErrorKind, MmcifLabelAddressError, MmcifLabelAddressErrorKind, MmcifModResError,
        MmcifModResErrorKind, MmcifOperationExprError, MmcifOperationExprErrorKind,
        MmcifSequenceIdError, MmcifSequenceIdErrorKind, MmcifSheetError, MmcifSheetErrorKind,
        MmcifSiftsUnpError, MmcifSiftsUnpErrorKind, MmcifTlsInfoError, MmcifTlsInfoErrorKind,
        PdbResidueKeyError, PdbResidueKeyErrorKind, bio_from_mmcif,
    };
    #[pymodule_export]
    use super::bio_mmcif_writer::{
        BioMmcifWriteError, BioMmcifWriteErrorKind, CifReadError, CifReadErrorKind,
    };
    #[pymodule_export]
    use super::bio_molecule_conversion::{
        BioMoleculeConversionError, BioMoleculeConversionErrorKind, BioMoleculeError,
        BioMoleculeErrorKind, PdbPostprocessError, PdbPostprocessErrorKind,
    };
    #[pymodule_export]
    use super::bio_pdb_reader::{
        _bio_from_pdb_keywords, _bio_from_pdb_params, BioPdbReadError, PdbAtomFieldError,
        PdbAtomFieldErrorKind, PdbChargeError, PdbCisPepRecordError, PdbCisPepRecordErrorKind,
        PdbConectError, PdbConectErrorKind, PdbConnectionFinalizationError,
        PdbConnectionFinalizationErrorKind, PdbCrystalError, PdbCrystalErrorKind, PdbDbRefError,
        PdbDbRefErrorKind, PdbHeaderError, PdbHeaderErrorKind, PdbHelixError, PdbHelixErrorKind,
        PdbInputFormatError, PdbInputFormatErrorKind, PdbLinkNameCounter, PdbLinkRecordError,
        PdbLinkRecordErrorKind, PdbModResError, PdbModResErrorKind, PdbModelTransitionError,
        PdbModelTransitionErrorKind, PdbMtrixError, PdbMtrixErrorKind, PdbReaderCause,
        PdbReaderCauseKind, PdbRemark3Error, PdbRemark3ErrorKind, PdbRemark200Error,
        PdbRemark200ErrorKind, PdbRemark350Error, PdbRemark350ErrorKind, PdbRemarkError,
        PdbRemarkErrorKind, PdbRemarkFinalizationError, PdbRemarkFinalizationErrorKind,
        PdbRemarkMetadataError, PdbRemarkMetadataErrorKind, PdbSheetError, PdbSheetErrorKind,
        PdbSsbondRecordError, PdbSsbondRecordErrorKind, PdbSsbondResolutionError,
        PdbStructureFinalizationError, PdbStructureFinalizationErrorKind, bio_from_pdb,
    };
    #[pymodule_export]
    use super::bio_pdb_writer::{BioPdbWriteError, BioPdbWriteErrorKind};
    #[pymodule_export]
    use super::bio_read_errors::{
        BioReadError, BioReadErrorKind, ChemCompAtomError, ChemCompAtomErrorKind,
        ChemCompModelError, ChemCompModelErrorKind, JsonError, JsonErrorKind, MmjsonReadError,
        MmjsonReadErrorKind, MmjsonStructureError, MmjsonStructureErrorKind, MmjsonUnexpectedType,
        MmjsonValueError, MmjsonValueErrorKind, ProteinReadError, ProteinReadErrorKind,
    };
    #[pymodule_export]
    use super::bio_reader_calls::{
        _protein_from_pdb_keywords, _protein_from_pdb_params, protein_from_mmcif, protein_from_pdb,
    };
    #[pymodule_export]
    use super::bio_write_configuration::BioMmcifWriteParams;
    #[pymodule_export]
    use super::mol2_post_errors::Mol2PostError;
    #[pymodule_export]
    use super::rdkit_adapter::mol_from_rdkit;

    #[pymodule_export]
    use super::fingerprint_configuration::{
        AtomPairFingerprintParams, AtomPairGeneratorParams, MorganGeneratorParams,
        TopologicalTorsionFingerprintParams, TopologicalTorsionGeneratorParams,
    };
    #[pymodule_export]
    use super::fingerprint_json::{
        FingerprintJsonError, FingerprintJsonErrorKind, SerdeJsonCategory, SerdeJsonError,
    };

    #[pymodule_export]
    use super::optional_configuration::{Coordinate2DParams, ReactionRunParams};

    #[pymodule_export]
    use super::property_values::{
        AtomPropertyReadError, AtomPropertyReadErrorKind, MissingPropertyError, PropertyValue,
        PropertyValueError, PropertyValueKind,
    };

    #[pymodule_export]
    use super::bio_values::{
        AtomName, BioChainRow, BioCoordinateFormat, BioCrystalCell, BioModelRow, BioReadParams,
        BioRowSpan, BioTransform, ChainKind, ChainSourceIds,
    };

    #[pymodule_export]
    use super::bio_errors::{BioStructureError, BioStructureErrorKind};

    #[pymodule_export]
    use super::collection_params::{
        AlignmentAtomMap, AlignmentParameters, LayeredFingerprintParams,
    };

    #[pymodule_export]
    use super::numerical_params::{
        AvalonFingerprintParams, CipLabelOptions, ConformerAlignmentParameters,
        LegacyTopologicalTorsionParams,
    };

    #[pymodule_export]
    use super::construction_values::{AtomSpec, BondSpec, ReactionCoordinateSelection};

    #[pymodule_export]
    use super::chemical_values::{
        BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization, McsAtomComparator,
        McsBondComparator,
    };

    #[pymodule_export]
    use super::fingerprint_explanations::{AtomCodeExplanation, explain_path_score};

    #[pymodule_export]
    use super::tautomer_callbacks::{TautomerMoleculeView, TautomerProgress};
    #[pymodule_export]
    use super::tautomer_values::{
        TautomerEnumeration, TautomerScore, TautomerScoreParams, TautomerScoreTerm,
        default_tautomer_score_terms,
    };

    #[pymodule_export]
    use super::chemical_params::{
        AddHsParams, ForceFieldMinimizeParams, MmffForceFieldParams, SmilesWriteParams,
        UffForceFieldParams,
    };

    #[pymodule_export]
    use super::selectors::{
        AromaticityModel, AromaticityParams, Coordinate2DInputParams, CoordinateZPolicy,
        SdfCoordinateMode, SdfReadParams, ValenceModel, ValenceParams,
    };

    #[pymodule_export]
    use super::coordinate_2d_errors::{
        Coordinate2DError, Coordinate2DErrorKind, Coordinate2DLayoutError,
        Coordinate2DLayoutErrorKind, Coordinate2DTemplateError, Coordinate2DTemplateErrorKind,
    };
    #[pymodule_export]
    use super::coordinates::{Conformer2D, Conformer3D};

    #[pymodule_export]
    use super::elements::{Element, ElementInfo, element_info};

    #[pymodule_export]
    use super::errors::{
        AtomCodeExplanationError, AtomCodeExplanationErrorKind, Coordinate3DReadError,
        Coordinate3DReadErrorKind, FingerprintError, FingerprintErrorKind, MaccsFingerprintError,
        MaccsFingerprintErrorKind, SmilesError, SmilesErrorKind,
    };

    #[pymodule_export]
    use super::fingerprint_providers::{
        AtomPairAtomInvariantsGenerator, FingerprintAdditionalOutput, MorganBondInvariantsGenerator,
    };

    #[pymodule_export]
    use super::fingerprint_values::{
        AtomPairsParameters, SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
    };

    #[pymodule_export]
    use super::params::{
        AtomPositionParams, BioMoleculeParams, BioPdbReadParams, BioPdbWriteParams,
        Coordinate3DInputParams, DistanceMatrix3dParams, DistanceMatrixParams,
        ForceFieldAngleConstraintParams, ForceFieldDistanceConstraintParams,
        ForceFieldPositionConstraintParams, ForceFieldTorsionConstraintParams, InchiReadParams,
        InchiWriteParams, KekulizeParams, MaccsFingerprintParams, McsAtomCompareParameters,
        McsBondCompareParameters, MmffConformerOptimizationParams, MmffEvaluationParams,
        MmffOptimizationParams, MmffPropertiesParams, MorganCallParams, PatternFingerprintParams,
        PotentialStereoParams, RandomSmilesWriteParams, ReactionApplyParams, ReactionParseParams,
        ReactionTemplateRemovalParams, ReactionValidationParams, RemoveHsParams,
        Replace3DCoordinatesParams, RingSearchParams, SmartsWriteParams, StructureTagParams,
        TopologicalFingerprintOutputRequest, TopologicalTorsionCallParams,
        UffConformerOptimizationParams, UffEvaluationParams, UffOptimizationParams, XyzWriteParams,
    };

    #[pyfunction]
    /// Return the COSMolKit package version compiled into this extension.
    fn version() -> &'static str {
        ck::version()
    }

    #[pyclass(module = "cosmolkit", from_py_object)]
    #[derive(Clone)]
    /// Molecular graph, properties and independently selected 2D/3D coordinates.
    /// Value-returning operations leave their source unchanged; methods ending
    /// in an underscore modify this value in place. copy() initially shares
    /// molecular storage, with changes isolated by copy-on-write.
    pub struct Molecule {
        pub(crate) inner: ck::Molecule,
    }

    impl Molecule {
        fn _substruct_match_params_native(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<super::search_callbacks::MatchInvocation>,
        ) -> PyResult<Option<MatchResult>> {
            let py = slf.py();
            let target = slf.clone().unbind();
            let holder = slf.borrow();
            let native = &holder.inner;
            let query = query.inner();
            let result = match params {
                Some(params) => params.with_callbacks(py, &target, |params| {
                    py.detach(|| native.substruct_match_with_params(query, params))
                })?,
                None => py.detach(|| native.substruct_match(query)),
            };
            result
                .map(|value| value.map(|inner| MatchResult { inner }))
                .map_err(SubstructMatchError::from_rust)
        }
        fn _substruct_matches_params_native(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<super::search_callbacks::MatchInvocation>,
        ) -> PyResult<Vec<MatchResult>> {
            let py = slf.py();
            let target = slf.clone().unbind();
            let holder = slf.borrow();
            let native = &holder.inner;
            let query = query.inner();
            let result = match params {
                Some(params) => params.with_callbacks(py, &target, |params| {
                    py.detach(|| native.substruct_matches_with_params(query, params))
                })?,
                None => py.detach(|| native.substruct_matches(query)),
            };
            result
                .map(|value| {
                    value
                        .into_iter()
                        .map(|inner| MatchResult { inner })
                        .collect()
                })
                .map_err(SubstructMatchError::from_rust)
        }
        fn _has_substruct_match_params_native(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<super::search_callbacks::MatchInvocation>,
        ) -> PyResult<bool> {
            let py = slf.py();
            let target = slf.clone().unbind();
            let holder = slf.borrow();
            let native = &holder.inner;
            let query = query.inner();
            let result = match params {
                Some(params) => params.with_callbacks(py, &target, |params| {
                    py.detach(|| native.has_substruct_match_with_params(query, params))
                })?,
                None => py.detach(|| native.has_substruct_match(query)),
            };
            result
                .map(|value| value)
                .map_err(SubstructMatchError::from_rust)
        }
        fn _chem_tautomer_score_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::TautomerScoreParams>,
        ) -> PyResult<super::tautomer_values::TautomerScore> {
            match params {
                Some(value) => py.detach(|| self.inner.tautomer_score_with_params(&value)),
                None => py.detach(|| self.inner.tautomer_score()),
            }
            .map(|inner| super::tautomer_values::TautomerScore { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        fn _with_3d_conformer_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<Molecule> {
            match params {
                Some(params) => py.detach(|| self.inner.with_3d_conformer_with_params(&params)),
                None => py.detach(|| self.inner.with_3d_conformer()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _embed_3d_conformer_params_native(
            &mut self,
            py: Python<'_>,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<()> {
            match params {
                Some(params) => py.detach(|| self.inner.embed_3d_conformer_with_params_(&params)),
                None => py.detach(|| self.inner.embed_3d_conformer_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _with_3d_conformer_result_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            match params {
                Some(params) => {
                    py.detach(|| self.inner.with_3d_conformer_result_with_params(&params))
                }
                None => py.detach(|| self.inner.with_3d_conformer_result()),
            }
            .map(|inner| super::embedding_values::EmbedMoleculeResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _embed_3d_conformer_result_params_native(
            &mut self,
            py: Python<'_>,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            match params {
                Some(params) => {
                    py.detach(|| self.inner.embed_3d_conformer_result_with_params_(&params))
                }
                None => py.detach(|| self.inner.embed_3d_conformer_result_()),
            }
            .map(|inner| super::embedding_values::EmbedMoleculeResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _with_3d_conformers_params_native(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<Molecule> {
            match params {
                Some(params) => py.detach(|| {
                    self.inner
                        .with_3d_conformers_with_params(num_confs, &params)
                }),
                None => py.detach(|| self.inner.with_3d_conformers(num_confs)),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _embed_3d_conformers_params_native(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<()> {
            match params {
                Some(params) => py.detach(|| {
                    self.inner
                        .embed_3d_conformers_with_params_(num_confs, &params)
                }),
                None => py.detach(|| self.inner.embed_3d_conformers_(num_confs)),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _with_3d_conformers_result_params_native(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            match params {
                Some(params) => py.detach(|| {
                    self.inner
                        .with_3d_conformers_result_with_params(num_confs, &params)
                }),
                None => py.detach(|| self.inner.with_3d_conformers_result(num_confs)),
            }
            .map(|inner| super::embedding_values::EmbedMultipleConfsResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _embed_3d_conformers_result_params_native(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&ck::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            match params {
                Some(params) => py.detach(|| {
                    self.inner
                        .embed_3d_conformers_result_with_params_(num_confs, &params)
                }),
                None => py.detach(|| self.inner.embed_3d_conformers_result_(num_confs)),
            }
            .map(|inner| super::embedding_values::EmbedMultipleConfsResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        fn _fingerprint_topological_torsion_legacy_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::LegacyTopologicalTorsionParams>,
        ) -> PyResult<Fingerprint> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .fingerprint_topological_torsion_legacy_with_params(&params),
                None => self.inner.fingerprint_topological_torsion_legacy(),
            })
            .map(|inner| Fingerprint { inner })
            .map_err(crate::morgan_errors::TopologicalTorsionReadError::from_rust)
        }

        fn _fingerprint_topological_torsion_count_legacy_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::LegacyTopologicalTorsionParams>,
        ) -> PyResult<SparseCountFingerprint> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .fingerprint_topological_torsion_count_legacy_with_params(&params),
                None => self.inner.fingerprint_topological_torsion_count_legacy(),
            })
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(crate::morgan_errors::TopologicalTorsionReadError::from_rust)
        }

        fn _fingerprint_topological_torsion_sparse_count_legacy_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::LegacyTopologicalTorsionParams>,
        ) -> PyResult<SparseCountFingerprint> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .fingerprint_topological_torsion_sparse_count_legacy_with_params(&params),
                None => self
                    .inner
                    .fingerprint_topological_torsion_sparse_count_legacy(),
            })
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(crate::morgan_errors::TopologicalTorsionReadError::from_rust)
        }

        fn _alignment_transform_to_params_native(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&ck::AlignmentParameters>,
        ) -> PyResult<AlignmentResult> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .alignment_transform_to_with_params(&reference.inner, &params),
                None => self.inner.alignment_transform_to(&reference.inner),
            })
            .map(|inner| AlignmentResult { inner })
            .map_err(crate::alignment_errors::AlignmentError::from_rust)
        }

        fn _best_alignment_to_params_native(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&ck::BestAlignmentParameters>,
        ) -> PyResult<AlignmentResult> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .best_alignment_to_with_params(&reference.inner, &params),
                None => self.inner.best_alignment_to(&reference.inner),
            })
            .map(|inner| AlignmentResult { inner })
            .map_err(crate::alignment_errors::AlignmentError::from_rust)
        }

        fn _best_rmsd_to_params_native(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&ck::BestAlignmentParameters>,
        ) -> PyResult<f64> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .best_rmsd_to_with_params(&reference.inner, &params),
                None => self.inner.best_rmsd_to(&reference.inner),
            })
            .map_err(crate::alignment_errors::AlignmentError::from_rust)
        }

        fn _coordinate_rmsd_to_params_native(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&ck::CoordinateRmsdParameters>,
        ) -> PyResult<f64> {
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .coordinate_rmsd_to_with_params(&reference.inner, &params),
                None => self.inner.coordinate_rmsd_to(&reference.inner),
            })
            .map_err(crate::alignment_errors::AlignmentError::from_rust)
        }

        fn _all_conformer_best_rmsds_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::AllConformerRmsdParameters>,
        ) -> PyResult<Vec<ConformerRmsd>> {
            py.detach(|| match params {
                Some(params) => self.inner.all_conformer_best_rmsds_with_params(&params),
                None => self.inner.all_conformer_best_rmsds(),
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|inner| ConformerRmsd { inner })
                    .collect()
            })
            .map_err(crate::alignment_errors::AlignmentError::from_rust)
        }

        fn _fingerprint_topological_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::TopologicalFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            match params {
                Some(params) => self.inner.fingerprint_topological_with_params(&params),
                None => self.inner.fingerprint_topological(),
            }
            .map(|inner| Fingerprint { inner })
            .map_err(TopologicalFingerprintError::from_rust)
        }

        fn _fingerprint_topological_with_output_params_native(
            &self,
            _py: Python<'_>,
            params: Option<&ck::TopologicalFingerprintParams>,
            request: Option<ck::TopologicalFingerprintOutputRequest>,
        ) -> PyResult<TopologicalFingerprintResult> {
            match (params, request) {
                (None, None) => self.inner.fingerprint_topological_with_output(),
                (Some(params), request) => {
                    self.inner.fingerprint_topological_with_output_with_params(
                        &params,
                        request.unwrap_or_default(),
                    )
                }
                (None, Some(request)) => {
                    self.inner.fingerprint_topological_with_output_with_params(
                        &ck::TopologicalFingerprintParams::default(),
                        request,
                    )
                }
            }
            .map(TopologicalFingerprintResult::from_rust)
            .map_err(TopologicalFingerprintError::from_rust)
        }

        fn _fingerprint_layered_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::LayeredFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            match params {
                Some(params) => self.inner.fingerprint_layered_with_params(&params),
                None => self.inner.fingerprint_layered(),
            }
            .map(|inner| Fingerprint { inner })
            .map_err(LayeredFingerprintError::from_rust)
        }

        fn _fingerprint_layered_with_output_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::LayeredFingerprintParams>,
        ) -> PyResult<LayeredFingerprintResult> {
            match params {
                Some(params) => self
                    .inner
                    .fingerprint_layered_with_output_with_params(&params),
                None => self.inner.fingerprint_layered_with_output(),
            }
            .map(|inner| LayeredFingerprintResult { inner })
            .map_err(LayeredFingerprintError::from_rust)
        }

        fn _chem_with_cip_labels_params_native(
            &self,
            py: Python<'_>,
            options: Option<&ck::CipLabelOptions>,
        ) -> PyResult<Molecule> {
            match options {
                Some(value) => py.detach(|| self.inner.with_cip_labels_with_options(&value)),
                None => py.detach(|| self.inner.with_cip_labels()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _chem_assign_cip_labels_params_native(
            &mut self,
            py: Python<'_>,
            options: Option<&ck::CipLabelOptions>,
        ) -> PyResult<()> {
            match options {
                Some(value) => py.detach(|| self.inner.assign_cip_labels_with_options_(&value)),
                None => py.detach(|| self.inner.assign_cip_labels_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _chem_with_hydrogens_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::AddHsParams>,
        ) -> PyResult<Molecule> {
            match params {
                Some(value) => py.detach(|| self.inner.with_hydrogens_with_params(&value)),
                None => py.detach(|| self.inner.with_hydrogens()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _chem_add_hydrogens_params_native(
            &mut self,
            py: Python<'_>,
            params: Option<&ck::AddHsParams>,
        ) -> PyResult<()> {
            match params {
                Some(value) => py.detach(|| self.inner.add_hydrogens_with_params_(&value)),
                None => py.detach(|| self.inner.add_hydrogens_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _alignment_with_alignment_to_params_native(
            &self,
            py: Python<'_>,
            reference: &Self,
            params: Option<&ck::AlignmentParameters>,
        ) -> PyResult<(Self, super::alignment_values::AlignmentResult)> {
            match params {
                Some(params) => py.detach(|| {
                    self.inner
                        .with_alignment_to_with_params(&reference.inner, &params)
                }),
                None => py.detach(|| self.inner.with_alignment_to(&reference.inner)),
            }
            .map(|(inner, report)| {
                (
                    Self { inner },
                    super::alignment_values::AlignmentResult { inner: report },
                )
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _alignment_with_alignment_to_inplace_params_native(
            &mut self,
            py: Python<'_>,
            reference: &Self,
            params: Option<&ck::AlignmentParameters>,
        ) -> PyResult<super::alignment_values::AlignmentResult> {
            match params {
                Some(params) => {
                    py.detach(|| self.inner.align_to_with_params_(&reference.inner, &params))
                }
                None => py.detach(|| self.inner.align_to_(&reference.inner)),
            }
            .map(|inner| super::alignment_values::AlignmentResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _alignment_with_aligned_conformers_params_native(
            &self,
            py: Python<'_>,
            params: Option<&ck::ConformerAlignmentParameters>,
        ) -> PyResult<(Self, super::alignment_values::ConformerAlignmentReport)> {
            match params {
                Some(params) => {
                    py.detach(|| self.inner.with_aligned_conformers_with_params(&params))
                }
                None => py.detach(|| self.inner.with_aligned_conformers()),
            }
            .map(|(inner, report)| {
                (
                    Self { inner },
                    super::alignment_values::ConformerAlignmentReport { inner: report },
                )
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        fn _alignment_with_aligned_conformers_inplace_params_native(
            &mut self,
            py: Python<'_>,
            params: Option<&ck::ConformerAlignmentParameters>,
        ) -> PyResult<super::alignment_values::ConformerAlignmentReport> {
            match params {
                Some(params) => py.detach(|| self.inner.align_conformers_with_params_(&params)),
                None => py.detach(|| self.inner.align_conformers_()),
            }
            .map(|inner| super::alignment_values::ConformerAlignmentReport { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
    }

    #[pymethods]
    impl Molecule {
        /// Return an independent molecule value using copy-on-write.
        /// Molecular storage is initially shared, not eagerly deep-copied.
        /// Subsequent modifications to either molecule do not affect the other.
        fn copy(&self) -> Self {
            Self {
                inner: self.inner.copy(),
            }
        }
        fn __copy__(&self) -> Self {
            self.copy()
        }
        #[pyo3(signature=(reaction,reactant_template,params=crate::arguments::Argument::Missing,*,coordinate_selection=crate::arguments::Argument::Missing))]
        /// Apply this molecule to the selected reactant template and return product sets. Matching callbacks and coordinate selection follow the supplied reaction policy.
        fn reaction_products(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            reactant_template: usize,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::search_reaction_params::ReactionSingleRunParams>>,
            >,
            coordinate_selection: crate::arguments::Argument<
                Option<PyRef<'_, super::construction_values::ReactionCoordinateSelection>>,
            >,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let keyword_form = coordinate_selection.is_supplied();
            crate::arguments::configuration_arguments!(params; coordinate_selection);
            let params = params.as_deref();
            let coordinate_selection = coordinate_selection.as_deref();
            if keyword_form {
                Self::_reaction_products_keywords(
                    slf,
                    reaction,
                    reactant_template,
                    coordinate_selection,
                )
            } else {
                Self::_reaction_products_params(slf, reaction, reactant_template, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.reaction_products"}"#]
        #[pyo3(signature=(reaction,reactant_template,params=None))]
        fn _reaction_products_params(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            reactant_template: usize,
            params: Option<&super::search_reaction_params::ReactionSingleRunParams>,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let py = slf.py();
            let owner = slf.borrow();
            let native = &owner.inner;
            reaction
                .invoke(py, &[native], |reaction| match params {
                    Some(params) => {
                        native.reaction_products_with_params(reaction, reactant_template, &params)
                    }
                    None => native.reaction_products(reaction, reactant_template),
                })
                .map(super::reaction_calls::products)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.reaction_products"}"#]
        #[pyo3(signature=(reaction,reactant_template,*,coordinate_selection=None))]
        fn _reaction_products_keywords(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            reactant_template: usize,
            coordinate_selection: Option<&super::construction_values::ReactionCoordinateSelection>,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let params = super::search_reaction_params::ReactionSingleRunParams {
                inner: crate::parameters::NativeField::Owned(ck::ReactionSingleRunParams {
                    coordinate_selection: coordinate_selection
                        .map(|value| value.inner)
                        .unwrap_or_else(|| {
                            ck::ReactionSingleRunParams::default().coordinate_selection
                        }),
                }),
            };
            Self::_reaction_products_params(slf, reaction, reactant_template, Some(&params))
        }
        /// Apply this molecule to the selected reactant template and return product sets. Matching callbacks and coordinate selection follow the supplied reaction policy. Uses the explicit configuration argument.
        fn reaction_products_with_params(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            reactant_template: usize,
            params: &super::search_reaction_params::ReactionSingleRunParams,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            Self::_reaction_products_params(slf, reaction, reactant_template, Some(params))
        }

        #[pyo3(signature=(reaction,reactants,params=crate::arguments::Argument::Missing,*,max_products=crate::arguments::Argument::Missing,coordinate_selections=crate::arguments::Argument::Missing,copy_atom_properties=crate::arguments::Argument::Missing))]
        /// Run the reaction with this molecule and the supplied additional reactant inputs, returning product sets in native order.
        fn reaction_products_from_inputs(
            &self,
            py: Python<'_>,
            reaction: &mut super::reaction::Reaction,
            reactants: Vec<Py<Molecule>>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::optional_configuration::ReactionRunParams>>,
            >,
            max_products: crate::arguments::Argument<Option<u32>>,
            coordinate_selections: crate::arguments::Argument<
                Option<Vec<super::construction_values::ReactionCoordinateSelection>>,
            >,
            copy_atom_properties: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let keyword_form = max_products.is_supplied()
                || coordinate_selections.is_supplied()
                || copy_atom_properties.is_supplied();
            crate::arguments::configuration_arguments!(params; max_products, coordinate_selections, copy_atom_properties);
            let params = params.as_deref();
            if keyword_form {
                self._reaction_inputs_keywords(
                    py,
                    reaction,
                    reactants,
                    max_products,
                    coordinate_selections,
                    copy_atom_properties,
                )
            } else {
                self._reaction_inputs_params(py, reaction, reactants, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.reaction_products_from_inputs"}"#]
        #[pyo3(signature=(reaction,reactants,params=None))]
        fn _reaction_inputs_params(
            &self,
            py: Python<'_>,
            reaction: &mut super::reaction::Reaction,
            reactants: Vec<Py<Molecule>>,
            params: Option<&super::optional_configuration::ReactionRunParams>,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let defaults;
            let params = match params {
                Some(params) => {
                    defaults = params.to_rust(py)?;
                    &defaults
                }
                None => {
                    defaults = ck::ReactionRunParams::default();
                    &defaults
                }
            };
            super::reaction_calls::run_from_inputs(self, reaction, py, &reactants, params)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.reaction_products_from_inputs"}"#]
        #[pyo3(signature=(reaction,reactants,*,max_products=None,coordinate_selections=None,copy_atom_properties=None))]
        fn _reaction_inputs_keywords(
            &self,
            py: Python<'_>,
            reaction: &mut super::reaction::Reaction,
            reactants: Vec<Py<Molecule>>,
            max_products: Option<u32>,
            coordinate_selections: Option<
                Vec<super::construction_values::ReactionCoordinateSelection>,
            >,
            copy_atom_properties: Option<bool>,
        ) -> PyResult<Vec<Vec<Molecule>>> {
            let params = super::reaction_calls::keywords(
                max_products,
                coordinate_selections,
                copy_atom_properties,
            );
            super::reaction_calls::run_from_inputs(self, reaction, py, &reactants, &params)
        }

        #[pyo3(signature=(reaction,params=crate::arguments::Argument::Missing,*,remove_unmatched_atoms=crate::arguments::Argument::Missing))]
        /// Return the native single-reactant application result as a new molecular value, leaving this molecule unchanged.
        fn apply_reaction(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::ReactionApplyParams>>,
            >,
            remove_unmatched_atoms: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::reaction_values::ReactionApplyResult> {
            let keyword_form = remove_unmatched_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; remove_unmatched_atoms);
            let params = params.as_deref();
            if keyword_form {
                Self::_apply_reaction_keywords(slf, reaction, remove_unmatched_atoms)
            } else {
                Self::_apply_reaction_params(slf, reaction, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.apply_reaction"}"#]
        #[pyo3(signature=(reaction,params=None))]
        fn _apply_reaction_params(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: Option<&super::params::ReactionApplyParams>,
        ) -> PyResult<super::reaction_values::ReactionApplyResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let py = slf.py();
            let owner = slf.borrow();
            let native = &owner.inner;
            reaction
                .invoke(py, &[native], |reaction| match params {
                    Some(params) => native.apply_reaction_with_params(reaction, &params),
                    None => native.apply_reaction(reaction),
                })
                .map(|inner| super::reaction_values::ReactionApplyResult { inner })
        }
        #[doc = r#"{"role":"overload","target":"Molecule.apply_reaction"}"#]
        #[pyo3(signature=(reaction,*,remove_unmatched_atoms=None))]
        fn _apply_reaction_keywords(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            remove_unmatched_atoms: Option<bool>,
        ) -> PyResult<super::reaction_values::ReactionApplyResult> {
            let params = super::params::ReactionApplyParams {
                inner: crate::parameters::NativeField::Owned(ck::ReactionApplyParams {
                    remove_unmatched_atoms: remove_unmatched_atoms.unwrap_or_else(|| {
                        ck::ReactionApplyParams::default().remove_unmatched_atoms
                    }),
                }),
            };
            Self::_apply_reaction_params(slf, reaction, Some(&params))
        }
        /// Return the native single-reactant application result as a new molecular value, leaving this molecule unchanged. Uses the explicit configuration argument.
        fn apply_reaction_with_params(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: &super::params::ReactionApplyParams,
        ) -> PyResult<super::reaction_values::ReactionApplyResult> {
            Self::_apply_reaction_params(slf, reaction, Some(params))
        }

        #[pyo3(signature=(reaction,params=crate::arguments::Argument::Missing,*,remove_unmatched_atoms=crate::arguments::Argument::Missing))]
        /// Apply a single-reactant reaction in place and return whether it changed the molecule. COW isolates other molecular values.
        fn apply_reaction_(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::ReactionApplyParams>>,
            >,
            remove_unmatched_atoms: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<bool> {
            let keyword_form = remove_unmatched_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; remove_unmatched_atoms);
            let params = params.as_deref();
            if keyword_form {
                Self::_apply_reaction_inplace_keywords(slf, reaction, remove_unmatched_atoms)
            } else {
                Self::_apply_reaction_inplace_params(slf, reaction, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.apply_reaction_"}"#]
        #[pyo3(signature=(reaction,params=None))]
        fn _apply_reaction_inplace_params(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: Option<&super::params::ReactionApplyParams>,
        ) -> PyResult<bool> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let py = slf.py();
            // Only callbacks requesting a target need a retained input value.
            // This is callback ownership, not a rollback snapshot or writeback.
            let snapshot = if reaction.needs_callback_target(py) {
                Some(slf.try_borrow()?.inner.clone())
            } else {
                None
            };
            let targets: Vec<_> = snapshot.iter().collect();
            let mut receiver = slf.try_borrow_mut()?;
            let native = &mut receiver.inner;
            reaction.invoke(py, &targets, |reaction| match params {
                Some(params) => native.apply_reaction_with_params_(reaction, &params),
                None => native.apply_reaction_(reaction),
            })
        }
        #[doc = r#"{"role":"overload","target":"Molecule.apply_reaction_"}"#]
        #[pyo3(signature=(reaction,*,remove_unmatched_atoms=None))]
        fn _apply_reaction_inplace_keywords(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            remove_unmatched_atoms: Option<bool>,
        ) -> PyResult<bool> {
            let params = super::params::ReactionApplyParams {
                inner: crate::parameters::NativeField::Owned(ck::ReactionApplyParams {
                    remove_unmatched_atoms: remove_unmatched_atoms.unwrap_or_else(|| {
                        ck::ReactionApplyParams::default().remove_unmatched_atoms
                    }),
                }),
            };
            Self::_apply_reaction_inplace_params(slf, reaction, Some(&params))
        }
        /// Apply a single-reactant reaction in place and return whether it changed the molecule. COW isolates other molecular values. Uses the explicit configuration argument.
        fn apply_reaction_with_params_(
            slf: Bound<'_, Self>,
            reaction: &mut super::reaction::Reaction,
            params: &super::params::ReactionApplyParams,
        ) -> PyResult<bool> {
            Self::_apply_reaction_inplace_params(slf, reaction, Some(params))
        }

        #[staticmethod]
        #[pyo3(signature=(rdmol, sanitize=None))]
        /// Import an RDKit molecule through the canonical adapter; sanitize selects whether to sanitize the imported molecular value.
        fn from_rdkit(
            py: Python<'_>,
            rdmol: super::rdkit_adapter::RdkitMol,
            sanitize: Option<bool>,
        ) -> PyResult<Self> {
            super::rdkit_adapter::import(py, &rdmol, sanitize)
        }
        /// Read the atom's signed 32-bit integer property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_int(&self, atom: usize, key: &str) -> PyResult<Option<i32>> {
            self.inner
                .atom_property_int(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's signed 32-bit integer property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_int(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: i32,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_int(ck::AtomId::new(atom), key, value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's signed 32-bit integer property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_int_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: i32,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_int_(ck::AtomId::new(atom), key, value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's unsigned 32-bit integer property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_uint(&self, atom: usize, key: &str) -> PyResult<Option<u32>> {
            self.inner
                .atom_property_uint(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's unsigned 32-bit integer property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_uint(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: u32,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_uint(ck::AtomId::new(atom), key, value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's unsigned 32-bit integer property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_uint_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: u32,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_uint_(ck::AtomId::new(atom), key, value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's 64-bit floating-point property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_double(&self, atom: usize, key: &str) -> PyResult<Option<f64>> {
            self.inner
                .atom_property_double(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's 64-bit floating-point property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_double(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: f64,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_double(ck::AtomId::new(atom), key, value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's 64-bit floating-point property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_double_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: f64,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_double_(ck::AtomId::new(atom), key, value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's boolean property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_bool(&self, atom: usize, key: &str) -> PyResult<Option<bool>> {
            self.inner
                .atom_property_bool(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's boolean property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_bool(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: bool,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_bool(ck::AtomId::new(atom), key, value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's boolean property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_bool_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: bool,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_bool_(ck::AtomId::new(atom), key, value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's text property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_string(&self, atom: usize, key: &str) -> PyResult<Option<String>> {
            self.inner
                .atom_property_string(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's text property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_string(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: &str,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_string(ck::AtomId::new(atom), key, value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's text property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_string_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: &str,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_string_(ck::AtomId::new(atom), key, value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's signed-integer vector property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_int_vector(&self, atom: usize, key: &str) -> PyResult<Option<Vec<i32>>> {
            self.inner
                .atom_property_int_vector(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's signed-integer vector property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_int_vector(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: Vec<i32>,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_int_vector(ck::AtomId::new(atom), key, &value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's signed-integer vector property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_int_vector_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: Vec<i32>,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_int_vector_(ck::AtomId::new(atom), key, &value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Read the atom's text-vector property. A missing key returns None;
        /// an incompatible stored type or invalid atom index raises a structured error.
        fn atom_property_string_vector(
            &self,
            atom: usize,
            key: &str,
        ) -> PyResult<Option<Vec<String>>> {
            self.inner
                .atom_property_string_vector(ck::AtomId::new(atom), key)
                .map_err(super::property_values::AtomPropertyReadError::from_rust)
        }
        /// Return a new molecule with the atom's text-vector property set to value.
        /// The source molecule is unchanged; atom is a zero-based index.
        fn with_atom_property_string_vector(
            &self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: Vec<String>,
        ) -> PyResult<Self> {
            py.detach(|| {
                self.inner
                    .with_atom_property_string_vector(ck::AtomId::new(atom), key, &value)
            })
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Set the atom's text-vector property in place. atom is a zero-based
        /// index; COW isolates other molecular values sharing the original storage.
        fn set_atom_property_string_vector_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            key: &str,
            value: Vec<String>,
        ) -> PyResult<()> {
            py.detach(|| {
                self.inner
                    .set_atom_property_string_vector_(ck::AtomId::new(atom), key, &value)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[pyo3(signature=(atom_id,branch_subtract=0,include_chirality=false,use_legacy_stereo_perception=true))]
        /// Return the native atom-pair code result for the selected atom and encoding options.
        fn with_atom_pair_atom_code(
            &self,
            py: Python<'_>,
            atom_id: usize,
            branch_subtract: u32,
            include_chirality: bool,
            use_legacy_stereo_perception: bool,
        ) -> PyResult<super::ordered_results::AtomPairAtomCodeResult> {
            py.detach(|| {
                self.inner.with_atom_pair_atom_code(
                    ck::AtomId::new(atom_id),
                    branch_subtract,
                    include_chirality,
                    use_legacy_stereo_perception,
                )
            })
            .map(|inner| super::ordered_results::AtomPairAtomCodeResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[staticmethod]
        #[pyo3(signature=(input,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
        /// Parse SDF text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_sdf(
            py: Python<'_>,
            input: &str,
            params: crate::arguments::Argument<Option<PyRef<'_, super::selectors::SdfReadParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            strict_parsing: crate::arguments::Argument<Option<bool>>,
            expand_attachment_points: crate::arguments::Argument<Option<bool>>,
            process_property_lists: crate::arguments::Argument<Option<bool>>,
            coordinate_mode: crate::arguments::Argument<
                Option<crate::selectors::SdfCoordinateModeInput>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || strict_parsing.is_supplied()
                || expand_attachment_points.is_supplied()
                || process_property_lists.is_supplied()
                || coordinate_mode.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
            let params = params.as_deref();
            if keyword_form {
                Self::_from_sdf_keywords(
                    py,
                    input,
                    sanitize,
                    remove_hs,
                    strict_parsing,
                    expand_attachment_points,
                    process_property_lists,
                    coordinate_mode,
                )
            } else {
                Self::_from_sdf_params(py, input, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_sdf"}"#]
        #[pyo3(signature=(input,params=None))]
        pub(crate) fn _from_sdf_params(
            py: Python<'_>,
            input: &str,
            params: Option<&super::selectors::SdfReadParams>,
        ) -> PyResult<Molecule> {
            super::sdf_records::molecule_sdf(py, input, params)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_sdf"}"#]
        #[pyo3(signature=(input,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
        pub(crate) fn _from_sdf_keywords(
            py: Python<'_>,
            input: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            strict_parsing: Option<bool>,
            expand_attachment_points: Option<bool>,
            process_property_lists: Option<bool>,
            coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
        ) -> PyResult<Molecule> {
            let params = super::sdf_records::configuration(
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )?;
            Self::_from_sdf_params(py, input, Some(&params))
        }
        #[staticmethod]
        /// Parse SDF text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path. Uses the explicit reader parameters.
        pub(crate) fn from_sdf_with_params(
            py: Python<'_>,
            input: &str,
            params: &super::selectors::SdfReadParams,
        ) -> PyResult<Molecule> {
            Self::_from_sdf_params(py, input, Some(params))
        }

        #[staticmethod]
        #[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
        /// Parse MDL MolBlock text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_mol(
            py: Python<'_>,
            text: &str,
            params: crate::arguments::Argument<Option<PyRef<'_, super::selectors::SdfReadParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            strict_parsing: crate::arguments::Argument<Option<bool>>,
            expand_attachment_points: crate::arguments::Argument<Option<bool>>,
            process_property_lists: crate::arguments::Argument<Option<bool>>,
            coordinate_mode: crate::arguments::Argument<
                Option<crate::selectors::SdfCoordinateModeInput>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || strict_parsing.is_supplied()
                || expand_attachment_points.is_supplied()
                || process_property_lists.is_supplied()
                || coordinate_mode.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
            let params = params.as_deref();
            if keyword_form {
                Self::_from_mol_keywords(
                    py,
                    text,
                    sanitize,
                    remove_hs,
                    strict_parsing,
                    expand_attachment_points,
                    process_property_lists,
                    coordinate_mode,
                )
            } else {
                Self::_from_mol_params(py, text, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mol"}"#]
        #[pyo3(signature=(text,params=None))]
        pub(crate) fn _from_mol_params(
            py: Python<'_>,
            text: &str,
            params: Option<&super::selectors::SdfReadParams>,
        ) -> PyResult<Molecule> {
            super::sdf_records::molecule_mol(py, text, params)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mol"}"#]
        #[pyo3(signature=(text,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
        pub(crate) fn _from_mol_keywords(
            py: Python<'_>,
            text: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            strict_parsing: Option<bool>,
            expand_attachment_points: Option<bool>,
            process_property_lists: Option<bool>,
            coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
        ) -> PyResult<Molecule> {
            let params = super::sdf_records::configuration(
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )?;
            Self::_from_mol_params(py, text, Some(&params))
        }
        #[staticmethod]
        /// Parse MDL MolBlock text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path. Uses the explicit reader parameters.
        pub(crate) fn from_mol_with_params(
            py: Python<'_>,
            text: &str,
            params: &super::selectors::SdfReadParams,
        ) -> PyResult<Molecule> {
            Self::_from_mol_params(py, text, Some(params))
        }

        #[pyo3(signature=(atom, position, params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing))]
        /// Return a new molecule with one atom position replaced under the selected coordinate and Z policy. The source is unchanged.
        fn with_atom_position(
            &self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: crate::arguments::Argument<Option<PyRef<'_, AtomPositionParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<Molecule> {
            let keyword_form = conformer_id.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id);
            let params = params.as_deref();
            if keyword_form {
                self._with_atom_position_keywords(py, atom, position, conformer_id)
            } else {
                self._with_atom_position_params(py, atom, position, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_atom_position"}"#]
        #[pyo3(signature=(atom, position, params=None))]
        fn _with_atom_position_params(
            &self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: Option<&AtomPositionParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.with_atom_position_with_params(
                    ck::AtomId::new(atom),
                    position,
                    &params,
                ),
                None => self
                    .inner
                    .with_atom_position(ck::AtomId::new(atom), position),
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_atom_position"}"#]
        #[pyo3(signature=(atom, position, *, conformer_id=None))]
        fn _with_atom_position_keywords(
            &self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            conformer_id: Option<usize>,
        ) -> PyResult<Molecule> {
            let params = AtomPositionParams {
                inner: crate::parameters::NativeField::Owned(ck::AtomPositionParams {
                    conformer_id,
                }),
            };
            self._with_atom_position_params(py, atom, position, Some(&params))
        }
        /// Return a new molecule with one atom position replaced under the selected coordinate and Z policy. The source is unchanged. Uses the explicit configuration argument.
        fn with_atom_position_with_params(
            &self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: &AtomPositionParams,
        ) -> PyResult<Molecule> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .with_atom_position_with_params(ck::AtomId::new(atom), position, &params)
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(atom, position, params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing))]
        /// Replace one atom position in place under the selected coordinate and Z policy.
        fn set_atom_position_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: crate::arguments::Argument<Option<PyRef<'_, AtomPositionParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<()> {
            let keyword_form = conformer_id.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id);
            let params = params.as_deref();
            if keyword_form {
                self._set_atom_position_keywords(py, atom, position, conformer_id)
            } else {
                self._set_atom_position_params(py, atom, position, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_atom_position_"}"#]
        #[pyo3(signature=(atom, position, params=None))]
        fn _set_atom_position_params(
            &mut self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: Option<&AtomPositionParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.set_atom_position_with_params_(
                    ck::AtomId::new(atom),
                    position,
                    &params,
                ),
                None => self
                    .inner
                    .set_atom_position_(ck::AtomId::new(atom), position),
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_atom_position_"}"#]
        #[pyo3(signature=(atom, position, *, conformer_id=None))]
        fn _set_atom_position_keywords(
            &mut self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            conformer_id: Option<usize>,
        ) -> PyResult<()> {
            let params = AtomPositionParams {
                inner: crate::parameters::NativeField::Owned(ck::AtomPositionParams {
                    conformer_id,
                }),
            };
            self._set_atom_position_params(py, atom, position, Some(&params))
        }
        /// Replace one atom position in place under the selected coordinate and Z policy. Uses the explicit configuration argument.
        fn set_atom_position_with_params_(
            &mut self,
            py: Python<'_>,
            atom: usize,
            position: [f64; 3],
            params: &AtomPositionParams,
        ) -> PyResult<()> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .set_atom_position_with_params_(ck::AtomId::new(atom), position, &params)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, coordinate_map=crate::arguments::Argument::Missing,canonical_orientation=crate::arguments::Argument::Missing,clear_existing_2d=crate::arguments::Argument::Missing,flips_per_sample=crate::arguments::Argument::Missing,samples=crate::arguments::Argument::Missing,sample_seed=crate::arguments::Argument::Missing,permute_degree_four=crate::arguments::Argument::Missing,force_rdkit=crate::arguments::Argument::Missing,use_ring_templates=crate::arguments::Argument::Missing,use_path_angle_expansion=crate::arguments::Argument::Missing))]
        /// Generate a 2D depiction layout in a new molecule. Supply Coordinate2DParams or keyword layout options, not both; the source is unchanged.
        fn with_2d_coordinates(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate2DParams>>>,
            coordinate_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<usize, [f64; 2]>>,
            >,
            canonical_orientation: crate::arguments::Argument<Option<bool>>,
            clear_existing_2d: crate::arguments::Argument<Option<bool>>,
            flips_per_sample: crate::arguments::Argument<Option<u32>>,
            samples: crate::arguments::Argument<Option<u32>>,
            sample_seed: crate::arguments::Argument<Option<i32>>,
            permute_degree_four: crate::arguments::Argument<Option<bool>>,
            force_rdkit: crate::arguments::Argument<Option<bool>>,
            use_ring_templates: crate::arguments::Argument<Option<bool>>,
            use_path_angle_expansion: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = coordinate_map.is_supplied()
                || canonical_orientation.is_supplied()
                || clear_existing_2d.is_supplied()
                || flips_per_sample.is_supplied()
                || samples.is_supplied()
                || sample_seed.is_supplied()
                || permute_degree_four.is_supplied()
                || force_rdkit.is_supplied()
                || use_ring_templates.is_supplied()
                || use_path_angle_expansion.is_supplied();
            crate::arguments::configuration_arguments!(params; coordinate_map, canonical_orientation, clear_existing_2d, flips_per_sample, samples, sample_seed, permute_degree_four, force_rdkit, use_ring_templates, use_path_angle_expansion);
            let params = params.as_deref();
            if keyword_form {
                self._with_2d_coordinates_keywords(
                    py,
                    coordinate_map,
                    canonical_orientation,
                    clear_existing_2d,
                    flips_per_sample,
                    samples,
                    sample_seed,
                    permute_degree_four,
                    force_rdkit,
                    use_ring_templates,
                    use_path_angle_expansion,
                )
            } else {
                self._with_2d_coordinates_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_2d_coordinates"}"#]
        #[pyo3(signature=(params=None))]
        fn _with_2d_coordinates_params(
            &self,
            py: Python<'_>,
            params: Option<&Coordinate2DParams>,
        ) -> PyResult<Molecule> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            super::coordinate_calls::with_layout(self, py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_2d_coordinates"}"#]
        #[pyo3(signature=(*, coordinate_map=None,canonical_orientation=None,clear_existing_2d=None,flips_per_sample=None,samples=None,sample_seed=None,permute_degree_four=None,force_rdkit=None,use_ring_templates=None,use_path_angle_expansion=None))]
        fn _with_2d_coordinates_keywords(
            &self,
            py: Python<'_>,
            coordinate_map: Option<std::collections::BTreeMap<usize, [f64; 2]>>,
            canonical_orientation: Option<bool>,
            clear_existing_2d: Option<bool>,
            flips_per_sample: Option<u32>,
            samples: Option<u32>,
            sample_seed: Option<i32>,
            permute_degree_four: Option<bool>,
            force_rdkit: Option<bool>,
            use_ring_templates: Option<bool>,
            use_path_angle_expansion: Option<bool>,
        ) -> PyResult<Molecule> {
            let params = super::coordinate_calls::layout(
                coordinate_map,
                canonical_orientation,
                clear_existing_2d,
                flips_per_sample,
                samples,
                sample_seed,
                permute_degree_four,
                force_rdkit,
                use_ring_templates,
                use_path_angle_expansion,
            );
            super::coordinate_calls::with_layout(self, py, Some(&params))
        }
        /// Generate a 2D depiction layout in a new molecule using Coordinate2DParams; the source is unchanged.
        fn with_2d_coordinates_with_params(
            &self,
            py: Python<'_>,
            params: &Coordinate2DParams,
        ) -> PyResult<Molecule> {
            self._with_2d_coordinates_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, coordinate_map=crate::arguments::Argument::Missing,canonical_orientation=crate::arguments::Argument::Missing,clear_existing_2d=crate::arguments::Argument::Missing,flips_per_sample=crate::arguments::Argument::Missing,samples=crate::arguments::Argument::Missing,sample_seed=crate::arguments::Argument::Missing,permute_degree_four=crate::arguments::Argument::Missing,force_rdkit=crate::arguments::Argument::Missing,use_ring_templates=crate::arguments::Argument::Missing,use_path_angle_expansion=crate::arguments::Argument::Missing))]
        /// Generate and install a 2D depiction layout in place. Supply Coordinate2DParams or keyword layout options, not both.
        fn add_2d_coordinates_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate2DParams>>>,
            coordinate_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<usize, [f64; 2]>>,
            >,
            canonical_orientation: crate::arguments::Argument<Option<bool>>,
            clear_existing_2d: crate::arguments::Argument<Option<bool>>,
            flips_per_sample: crate::arguments::Argument<Option<u32>>,
            samples: crate::arguments::Argument<Option<u32>>,
            sample_seed: crate::arguments::Argument<Option<i32>>,
            permute_degree_four: crate::arguments::Argument<Option<bool>>,
            force_rdkit: crate::arguments::Argument<Option<bool>>,
            use_ring_templates: crate::arguments::Argument<Option<bool>>,
            use_path_angle_expansion: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = coordinate_map.is_supplied()
                || canonical_orientation.is_supplied()
                || clear_existing_2d.is_supplied()
                || flips_per_sample.is_supplied()
                || samples.is_supplied()
                || sample_seed.is_supplied()
                || permute_degree_four.is_supplied()
                || force_rdkit.is_supplied()
                || use_ring_templates.is_supplied()
                || use_path_angle_expansion.is_supplied();
            crate::arguments::configuration_arguments!(params; coordinate_map, canonical_orientation, clear_existing_2d, flips_per_sample, samples, sample_seed, permute_degree_four, force_rdkit, use_ring_templates, use_path_angle_expansion);
            let params = params.as_deref();
            if keyword_form {
                self._add_2d_coordinates_keywords(
                    py,
                    coordinate_map,
                    canonical_orientation,
                    clear_existing_2d,
                    flips_per_sample,
                    samples,
                    sample_seed,
                    permute_degree_four,
                    force_rdkit,
                    use_ring_templates,
                    use_path_angle_expansion,
                )
            } else {
                self._add_2d_coordinates_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_2d_coordinates_"}"#]
        #[pyo3(signature=(params=None))]
        fn _add_2d_coordinates_params(
            &mut self,
            py: Python<'_>,
            params: Option<&Coordinate2DParams>,
        ) -> PyResult<()> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            super::coordinate_calls::add_layout(self, py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_2d_coordinates_"}"#]
        #[pyo3(signature=(*, coordinate_map=None,canonical_orientation=None,clear_existing_2d=None,flips_per_sample=None,samples=None,sample_seed=None,permute_degree_four=None,force_rdkit=None,use_ring_templates=None,use_path_angle_expansion=None))]
        fn _add_2d_coordinates_keywords(
            &mut self,
            py: Python<'_>,
            coordinate_map: Option<std::collections::BTreeMap<usize, [f64; 2]>>,
            canonical_orientation: Option<bool>,
            clear_existing_2d: Option<bool>,
            flips_per_sample: Option<u32>,
            samples: Option<u32>,
            sample_seed: Option<i32>,
            permute_degree_four: Option<bool>,
            force_rdkit: Option<bool>,
            use_ring_templates: Option<bool>,
            use_path_angle_expansion: Option<bool>,
        ) -> PyResult<()> {
            let params = super::coordinate_calls::layout(
                coordinate_map,
                canonical_orientation,
                clear_existing_2d,
                flips_per_sample,
                samples,
                sample_seed,
                permute_degree_four,
                force_rdkit,
                use_ring_templates,
                use_path_angle_expansion,
            );
            super::coordinate_calls::add_layout(self, py, Some(&params))
        }
        /// Generate and install a 2D depiction layout in place using Coordinate2DParams.
        fn add_2d_coordinates_with_params_(
            &mut self,
            py: Python<'_>,
            params: &Coordinate2DParams,
        ) -> PyResult<()> {
            self._add_2d_coordinates_params(py, Some(params))
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, z_policy=crate::arguments::Argument::Missing))]
        /// Return a new molecule with the supplied 2D coordinate block installed. This installs coordinates rather than generating a layout.
        fn with_2d_coordinate_block(
            &self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate2DInputParams>>>,
            z_policy: crate::arguments::Argument<Option<super::selectors::CoordinateZPolicyInput>>,
        ) -> PyResult<Molecule> {
            let keyword_form = z_policy.is_supplied();
            crate::arguments::configuration_arguments!(params; z_policy);
            let params = params.as_deref();
            if keyword_form {
                self._with_2d_coordinate_block_input_keywords(py, coordinates, z_policy)
            } else {
                self._with_2d_coordinate_block_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_2d_coordinate_block"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _with_2d_coordinate_block_input_params(
            &self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: Option<&Coordinate2DInputParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .with_2d_coordinate_block_with_params(coordinates, &params),
                None => self.inner.with_2d_coordinate_block(coordinates),
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_2d_coordinate_block"}"#]
        #[pyo3(signature=(coordinates, *, z_policy=None))]
        fn _with_2d_coordinate_block_input_keywords(
            &self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            z_policy: Option<super::selectors::CoordinateZPolicyInput>,
        ) -> PyResult<Molecule> {
            let params = super::coordinate_calls::input_2d(z_policy)?;
            self._with_2d_coordinate_block_input_params(py, coordinates, Some(&params))
        }
        /// Return a new molecule with the supplied 2D coordinate block installed. This installs coordinates rather than generating a layout. Uses the explicit configuration argument.
        fn with_2d_coordinate_block_with_params(
            &self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: &Coordinate2DInputParams,
        ) -> PyResult<Molecule> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .with_2d_coordinate_block_with_params(coordinates, &params)
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, z_policy=crate::arguments::Argument::Missing))]
        /// Install the supplied 2D coordinates in place. Coordinate rows must match the atom count and the declared coordinate-input policy.
        fn set_2d_coordinates_(
            &mut self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate2DInputParams>>>,
            z_policy: crate::arguments::Argument<Option<super::selectors::CoordinateZPolicyInput>>,
        ) -> PyResult<()> {
            let keyword_form = z_policy.is_supplied();
            crate::arguments::configuration_arguments!(params; z_policy);
            let params = params.as_deref();
            if keyword_form {
                self._set_2d_coordinates_input_keywords(py, coordinates, z_policy)
            } else {
                self._set_2d_coordinates_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_2d_coordinates_"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _set_2d_coordinates_input_params(
            &mut self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: Option<&Coordinate2DInputParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .set_2d_coordinates_with_params_(coordinates, &params),
                None => self.inner.set_2d_coordinates_(coordinates),
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_2d_coordinates_"}"#]
        #[pyo3(signature=(coordinates, *, z_policy=None))]
        fn _set_2d_coordinates_input_keywords(
            &mut self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            z_policy: Option<super::selectors::CoordinateZPolicyInput>,
        ) -> PyResult<()> {
            let params = super::coordinate_calls::input_2d(z_policy)?;
            self._set_2d_coordinates_input_params(py, coordinates, Some(&params))
        }
        /// Install the supplied 2D coordinates in place. Coordinate rows must match the atom count and the declared coordinate-input policy. Uses the explicit configuration argument.
        fn set_2d_coordinates_with_params_(
            &mut self,
            py: Python<'_>,
            coordinates: Vec<Vec<f64>>,
            params: &Coordinate2DInputParams,
        ) -> PyResult<()> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .set_2d_coordinates_with_params_(coordinates, &params)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing))]
        /// Return a new molecule with supplied 3D coordinate rows installed under the selected input policy. Does not run conformer embedding.
        fn with_3d_coordinates(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Replace3DCoordinatesParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<Molecule> {
            let keyword_form = conformer_id.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id);
            let params = params.as_deref();
            if keyword_form {
                self._with_3d_coordinates_input_keywords(py, coordinates, conformer_id)
            } else {
                self._with_3d_coordinates_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_coordinates"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _with_3d_coordinates_input_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Replace3DCoordinatesParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .with_3d_coordinates_with_params(coordinates, &params),
                None => self.inner.with_3d_coordinates(coordinates),
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_coordinates"}"#]
        #[pyo3(signature=(coordinates, *, conformer_id=None))]
        fn _with_3d_coordinates_input_keywords(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            conformer_id: Option<usize>,
        ) -> PyResult<Molecule> {
            let params = Replace3DCoordinatesParams {
                inner: crate::parameters::NativeField::Owned(ck::Replace3DCoordinatesParams {
                    conformer_id: conformer_id
                        .unwrap_or(ck::Replace3DCoordinatesParams::default().conformer_id),
                }),
            };
            self._with_3d_coordinates_input_params(py, coordinates, Some(&params))
        }
        /// Return a new molecule with supplied 3D coordinate rows installed under the selected input policy. Does not run conformer embedding. Uses the explicit configuration argument.
        fn with_3d_coordinates_with_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Replace3DCoordinatesParams,
        ) -> PyResult<Molecule> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .with_3d_coordinates_with_params(coordinates, &params)
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing))]
        /// Install supplied 3D coordinate rows in place under the selected input policy. Does not run conformer embedding.
        fn set_3d_coordinates_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Replace3DCoordinatesParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<()> {
            let keyword_form = conformer_id.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id);
            let params = params.as_deref();
            if keyword_form {
                self._set_3d_coordinates_input_keywords(py, coordinates, conformer_id)
            } else {
                self._set_3d_coordinates_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_3d_coordinates_"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _set_3d_coordinates_input_params(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Replace3DCoordinatesParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .set_3d_coordinates_with_params_(coordinates, &params),
                None => self.inner.set_3d_coordinates_(coordinates),
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_3d_coordinates_"}"#]
        #[pyo3(signature=(coordinates, *, conformer_id=None))]
        fn _set_3d_coordinates_input_keywords(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            conformer_id: Option<usize>,
        ) -> PyResult<()> {
            let params = Replace3DCoordinatesParams {
                inner: crate::parameters::NativeField::Owned(ck::Replace3DCoordinatesParams {
                    conformer_id: conformer_id
                        .unwrap_or(ck::Replace3DCoordinatesParams::default().conformer_id),
                }),
            };
            self._set_3d_coordinates_input_params(py, coordinates, Some(&params))
        }
        /// Install supplied 3D coordinate rows in place under the selected input policy. Does not run conformer embedding. Uses the explicit configuration argument.
        fn set_3d_coordinates_with_params_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Replace3DCoordinatesParams,
        ) -> PyResult<()> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .set_3d_coordinates_with_params_(coordinates, &params)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, is_3d=crate::arguments::Argument::Missing))]
        /// Return a new molecule with an additional 3D conformer containing the supplied coordinates.
        fn with_added_3d_conformer(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate3DInputParams>>>,
            is_3d: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = is_3d.is_supplied();
            crate::arguments::configuration_arguments!(params; is_3d);
            let params = params.as_deref();
            if keyword_form {
                self._with_added_3d_conformer_input_keywords(py, coordinates, is_3d)
            } else {
                self._with_added_3d_conformer_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_added_3d_conformer"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _with_added_3d_conformer_input_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Coordinate3DInputParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .with_added_3d_conformer_with_params(coordinates, &params),
                None => self.inner.with_added_3d_conformer(coordinates),
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_added_3d_conformer"}"#]
        #[pyo3(signature=(coordinates, *, is_3d=None))]
        fn _with_added_3d_conformer_input_keywords(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            is_3d: Option<bool>,
        ) -> PyResult<Molecule> {
            let params = Coordinate3DInputParams {
                inner: crate::parameters::NativeField::Owned(ck::Coordinate3DInputParams {
                    is_3d: is_3d.unwrap_or(ck::Coordinate3DInputParams::default().is_3d),
                }),
            };
            self._with_added_3d_conformer_input_params(py, coordinates, Some(&params))
        }
        /// Return a new molecule with an additional 3D conformer containing the supplied coordinates. Uses the explicit configuration argument.
        fn with_added_3d_conformer_with_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Coordinate3DInputParams,
        ) -> PyResult<Molecule> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .with_added_3d_conformer_with_params(coordinates, &params)
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, is_3d=crate::arguments::Argument::Missing))]
        /// Add a 3D conformer from supplied coordinates in place; existing conformers are retained.
        fn add_3d_conformer_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate3DInputParams>>>,
            is_3d: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<usize> {
            let keyword_form = is_3d.is_supplied();
            crate::arguments::configuration_arguments!(params; is_3d);
            let params = params.as_deref();
            if keyword_form {
                self._add_3d_conformer_input_keywords(py, coordinates, is_3d)
            } else {
                self._add_3d_conformer_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_3d_conformer_"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _add_3d_conformer_input_params(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Coordinate3DInputParams>,
        ) -> PyResult<usize> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .add_3d_conformer_with_params_(coordinates, &params),
                None => self.inner.add_3d_conformer_(coordinates),
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_3d_conformer_"}"#]
        #[pyo3(signature=(coordinates, *, is_3d=None))]
        fn _add_3d_conformer_input_keywords(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            is_3d: Option<bool>,
        ) -> PyResult<usize> {
            let params = Coordinate3DInputParams {
                inner: crate::parameters::NativeField::Owned(ck::Coordinate3DInputParams {
                    is_3d: is_3d.unwrap_or(ck::Coordinate3DInputParams::default().is_3d),
                }),
            };
            self._add_3d_conformer_input_params(py, coordinates, Some(&params))
        }
        /// Add a 3D conformer from supplied coordinates in place; existing conformers are retained. Uses the explicit configuration argument.
        fn add_3d_conformer_with_params_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Coordinate3DInputParams,
        ) -> PyResult<usize> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .add_3d_conformer_with_params_(coordinates, &params)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, is_3d=crate::arguments::Argument::Missing))]
        /// Return a new molecule whose only 3D conformer contains the supplied coordinates. Rows contain exactly three numeric components.
        fn with_only_3d_conformer(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate3DInputParams>>>,
            is_3d: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = is_3d.is_supplied();
            crate::arguments::configuration_arguments!(params; is_3d);
            let params = params.as_deref();
            if keyword_form {
                self._with_only_3d_conformer_input_keywords(py, coordinates, is_3d)
            } else {
                self._with_only_3d_conformer_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_only_3d_conformer"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _with_only_3d_conformer_input_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Coordinate3DInputParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .with_only_3d_conformer_with_params(coordinates, &params),
                None => self.inner.with_only_3d_conformer(coordinates),
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_only_3d_conformer"}"#]
        #[pyo3(signature=(coordinates, *, is_3d=None))]
        fn _with_only_3d_conformer_input_keywords(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            is_3d: Option<bool>,
        ) -> PyResult<Molecule> {
            let params = Coordinate3DInputParams {
                inner: crate::parameters::NativeField::Owned(ck::Coordinate3DInputParams {
                    is_3d: is_3d.unwrap_or(ck::Coordinate3DInputParams::default().is_3d),
                }),
            };
            self._with_only_3d_conformer_input_params(py, coordinates, Some(&params))
        }
        /// Return a new molecule whose only 3D conformer contains the supplied coordinates. Rows contain exactly three numeric components. Uses the explicit configuration argument.
        fn with_only_3d_conformer_with_params(
            &self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Coordinate3DInputParams,
        ) -> PyResult<Molecule> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .with_only_3d_conformer_with_params(coordinates, &params)
            })
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(coordinates, params=crate::arguments::Argument::Missing, *, is_3d=crate::arguments::Argument::Missing))]
        /// Replace all 3D conformers with one containing the supplied coordinates. Each row contains exactly three numeric components.
        fn set_only_3d_conformer_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, Coordinate3DInputParams>>>,
            is_3d: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<usize> {
            let keyword_form = is_3d.is_supplied();
            crate::arguments::configuration_arguments!(params; is_3d);
            let params = params.as_deref();
            if keyword_form {
                self._set_only_3d_conformer_input_keywords(py, coordinates, is_3d)
            } else {
                self._set_only_3d_conformer_input_params(py, coordinates, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_only_3d_conformer_"}"#]
        #[pyo3(signature=(coordinates, params=None))]
        fn _set_only_3d_conformer_input_params(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: Option<&Coordinate3DInputParams>,
        ) -> PyResult<usize> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .set_only_3d_conformer_with_params_(coordinates, &params),
                None => self.inner.set_only_3d_conformer_(coordinates),
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.set_only_3d_conformer_"}"#]
        #[pyo3(signature=(coordinates, *, is_3d=None))]
        fn _set_only_3d_conformer_input_keywords(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            is_3d: Option<bool>,
        ) -> PyResult<usize> {
            let params = Coordinate3DInputParams {
                inner: crate::parameters::NativeField::Owned(ck::Coordinate3DInputParams {
                    is_3d: is_3d.unwrap_or(ck::Coordinate3DInputParams::default().is_3d),
                }),
            };
            self._set_only_3d_conformer_input_params(py, coordinates, Some(&params))
        }
        /// Replace all 3D conformers with one containing the supplied coordinates. Each row contains exactly three numeric components. Uses the explicit configuration argument.
        fn set_only_3d_conformer_with_params_(
            &mut self,
            py: Python<'_>,
            coordinates: crate::coordinate_input::Coordinates3D<'_>,
            params: &Coordinate3DInputParams,
        ) -> PyResult<usize> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            let coordinates = coordinates.into_rows(self.inner.num_atoms())?;

            py.detach(|| {
                self.inner
                    .set_only_3d_conformer_with_params_(coordinates, &params)
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return a new molecule with all 3D conformers removed. The source molecule and its 2D coordinate block are unchanged.
        fn with_cleared_3d_conformers(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.with_cleared_3d_conformers())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Remove all 3D conformers in place. The 2D coordinate block is unchanged; returns None.
        fn clear_3d_conformers_(&mut self, py: Python<'_>) -> PyResult<()> {
            py.detach(|| self.inner.clear_3d_conformers_())
                .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing, vdw_threshold=crate::arguments::Argument::Missing, ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Construct a persistent UFF force field for the selected conformer, for evaluation, constraints and minimization without rebuilding the field.
        fn uff_force_field(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, UffForceFieldParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            vdw_threshold: crate::arguments::Argument<Option<f64>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<MolecularForceField> {
            let keyword_form = conformer_id.is_supplied()
                || vdw_threshold.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, vdw_threshold, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._uff_field_keywords(
                    py,
                    conformer_id,
                    vdw_threshold,
                    ignore_interfragment_interactions,
                )
            } else {
                self._uff_field_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.uff_force_field"}"#]
        #[pyo3(signature=(params=None))]
        fn _uff_field_params(
            &self,
            py: Python<'_>,
            params: Option<&UffForceFieldParams>,
        ) -> PyResult<MolecularForceField> {
            self.uff_field_params(py, params)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.uff_force_field"}"#]
        #[pyo3(signature=(*, conformer_id=None, vdw_threshold=None, ignore_interfragment_interactions=None))]
        fn _uff_field_keywords(
            &self,
            py: Python<'_>,
            conformer_id: Option<usize>,
            vdw_threshold: Option<f64>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<MolecularForceField> {
            self.uff_field_keywords(
                py,
                conformer_id,
                vdw_threshold,
                ignore_interfragment_interactions,
            )
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing, mmff_variant=crate::arguments::Argument::Missing, non_bonded_threshold=crate::arguments::Argument::Missing, ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Construct a persistent MMFF force field for the selected conformer and MMFF variant, for evaluation, constraints and minimization without rebuilding the field.
        fn mmff_force_field(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, MmffForceFieldParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            mmff_variant: crate::arguments::Argument<Option<crate::leaf_values::MmffVariantInput>>,
            non_bonded_threshold: crate::arguments::Argument<Option<f64>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<MolecularForceField> {
            let keyword_form = conformer_id.is_supplied()
                || mmff_variant.is_supplied()
                || non_bonded_threshold.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, mmff_variant, non_bonded_threshold, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._mmff_field_keywords(
                    py,
                    conformer_id,
                    mmff_variant,
                    non_bonded_threshold,
                    ignore_interfragment_interactions,
                )
            } else {
                self._mmff_field_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_force_field"}"#]
        #[pyo3(signature=(params=None))]
        fn _mmff_field_params(
            &self,
            py: Python<'_>,
            params: Option<&MmffForceFieldParams>,
        ) -> PyResult<MolecularForceField> {
            self.mmff_field_params(py, params)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_force_field"}"#]
        #[pyo3(signature=(*, conformer_id=None, mmff_variant=None, non_bonded_threshold=None, ignore_interfragment_interactions=None))]
        fn _mmff_field_keywords(
            &self,
            py: Python<'_>,
            conformer_id: Option<usize>,
            mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
            non_bonded_threshold: Option<f64>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<MolecularForceField> {
            self.mmff_field_keywords(
                py,
                conformer_id,
                mmff_variant,
                non_bonded_threshold,
                ignore_interfragment_interactions,
            )
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,torsion_atom_count=crate::arguments::Argument::Missing,include_chirality=crate::arguments::Argument::Missing,fp_size=crate::arguments::Argument::Missing,bits_per_entry=crate::arguments::Argument::Missing,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing))]
        /// Calculate the legacy topological-torsion bit fingerprint with the selected torsion length and chirality policy.
        fn fingerprint_topological_torsion_legacy(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, LegacyTopologicalTorsionParams>>>,
            torsion_atom_count: crate::arguments::Argument<Option<u32>>,
            include_chirality: crate::arguments::Argument<Option<bool>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            bits_per_entry: crate::arguments::Argument<Option<u32>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = torsion_atom_count.is_supplied()
                || include_chirality.is_supplied()
                || fp_size.is_supplied()
                || bits_per_entry.is_supplied()
                || from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied();
            crate::arguments::configuration_arguments!(params; torsion_atom_count, include_chirality, fp_size, bits_per_entry, from_atoms, ignore_atoms, custom_atom_invariants);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_legacy_keywords(
                    py,
                    torsion_atom_count,
                    include_chirality,
                    fp_size,
                    bits_per_entry,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                )
            } else {
                self._fingerprint_topological_torsion_legacy_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_legacy"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_topological_torsion_legacy_params(
            &self,
            py: Python<'_>,
            params: Option<&LegacyTopologicalTorsionParams>,
        ) -> PyResult<Fingerprint> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_torsion_legacy_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_legacy"}"#]
        #[pyo3(signature=(*,torsion_atom_count=None,include_chirality=None,fp_size=None,bits_per_entry=None,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None))]
        fn _fingerprint_topological_torsion_legacy_keywords(
            &self,
            py: Python<'_>,
            torsion_atom_count: Option<u32>,
            include_chirality: Option<bool>,
            fp_size: Option<u32>,
            bits_per_entry: Option<u32>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
        ) -> PyResult<Fingerprint> {
            let params = crate::legacy_torsion_calls::keywords(
                torsion_atom_count,
                include_chirality,
                fp_size,
                bits_per_entry,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
            );
            self._fingerprint_topological_torsion_legacy_params_native(py, Some(&params))
        }
        /// Calculate the legacy topological-torsion bit fingerprint with the selected torsion length and chirality policy. Uses the explicit configuration argument.
        fn fingerprint_topological_torsion_legacy_with_params(
            &self,
            py: Python<'_>,
            params: &LegacyTopologicalTorsionParams,
        ) -> PyResult<Fingerprint> {
            self._fingerprint_topological_torsion_legacy_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,torsion_atom_count=crate::arguments::Argument::Missing,include_chirality=crate::arguments::Argument::Missing,fp_size=crate::arguments::Argument::Missing,bits_per_entry=crate::arguments::Argument::Missing,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing))]
        /// Calculate the legacy topological-torsion count fingerprint with the selected encoding policy.
        fn fingerprint_topological_torsion_count_legacy(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, LegacyTopologicalTorsionParams>>>,
            torsion_atom_count: crate::arguments::Argument<Option<u32>>,
            include_chirality: crate::arguments::Argument<Option<bool>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            bits_per_entry: crate::arguments::Argument<Option<u32>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<SparseCountFingerprint> {
            let keyword_form = torsion_atom_count.is_supplied()
                || include_chirality.is_supplied()
                || fp_size.is_supplied()
                || bits_per_entry.is_supplied()
                || from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied();
            crate::arguments::configuration_arguments!(params; torsion_atom_count, include_chirality, fp_size, bits_per_entry, from_atoms, ignore_atoms, custom_atom_invariants);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_count_legacy_keywords(
                    py,
                    torsion_atom_count,
                    include_chirality,
                    fp_size,
                    bits_per_entry,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                )
            } else {
                self._fingerprint_topological_torsion_count_legacy_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count_legacy"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_topological_torsion_count_legacy_params(
            &self,
            py: Python<'_>,
            params: Option<&LegacyTopologicalTorsionParams>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_torsion_count_legacy_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count_legacy"}"#]
        #[pyo3(signature=(*,torsion_atom_count=None,include_chirality=None,fp_size=None,bits_per_entry=None,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None))]
        fn _fingerprint_topological_torsion_count_legacy_keywords(
            &self,
            py: Python<'_>,
            torsion_atom_count: Option<u32>,
            include_chirality: Option<bool>,
            fp_size: Option<u32>,
            bits_per_entry: Option<u32>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::legacy_torsion_calls::keywords(
                torsion_atom_count,
                include_chirality,
                fp_size,
                bits_per_entry,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
            );
            self._fingerprint_topological_torsion_count_legacy_params_native(py, Some(&params))
        }
        /// Calculate the legacy topological-torsion count fingerprint with the selected encoding policy. Uses the explicit configuration argument.
        fn fingerprint_topological_torsion_count_legacy_with_params(
            &self,
            py: Python<'_>,
            params: &LegacyTopologicalTorsionParams,
        ) -> PyResult<SparseCountFingerprint> {
            self._fingerprint_topological_torsion_count_legacy_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,torsion_atom_count=crate::arguments::Argument::Missing,include_chirality=crate::arguments::Argument::Missing,fp_size=crate::arguments::Argument::Missing,bits_per_entry=crate::arguments::Argument::Missing,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing))]
        /// Calculate sparse legacy topological-torsion counts keyed by torsion identifiers.
        fn fingerprint_topological_torsion_sparse_count_legacy(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, LegacyTopologicalTorsionParams>>>,
            torsion_atom_count: crate::arguments::Argument<Option<u32>>,
            include_chirality: crate::arguments::Argument<Option<bool>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            bits_per_entry: crate::arguments::Argument<Option<u32>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<SparseCountFingerprint> {
            let keyword_form = torsion_atom_count.is_supplied()
                || include_chirality.is_supplied()
                || fp_size.is_supplied()
                || bits_per_entry.is_supplied()
                || from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied();
            crate::arguments::configuration_arguments!(params; torsion_atom_count, include_chirality, fp_size, bits_per_entry, from_atoms, ignore_atoms, custom_atom_invariants);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_sparse_count_legacy_keywords(
                    py,
                    torsion_atom_count,
                    include_chirality,
                    fp_size,
                    bits_per_entry,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                )
            } else {
                self._fingerprint_topological_torsion_sparse_count_legacy_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count_legacy"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_topological_torsion_sparse_count_legacy_params(
            &self,
            py: Python<'_>,
            params: Option<&LegacyTopologicalTorsionParams>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_torsion_sparse_count_legacy_params_native(
                py,
                params.as_ref(),
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count_legacy"}"#]
        #[pyo3(signature=(*,torsion_atom_count=None,include_chirality=None,fp_size=None,bits_per_entry=None,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None))]
        fn _fingerprint_topological_torsion_sparse_count_legacy_keywords(
            &self,
            py: Python<'_>,
            torsion_atom_count: Option<u32>,
            include_chirality: Option<bool>,
            fp_size: Option<u32>,
            bits_per_entry: Option<u32>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::legacy_torsion_calls::keywords(
                torsion_atom_count,
                include_chirality,
                fp_size,
                bits_per_entry,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
            );
            self._fingerprint_topological_torsion_sparse_count_legacy_params_native(
                py,
                Some(&params),
            )
        }
        /// Calculate sparse legacy topological-torsion counts keyed by torsion identifiers. Uses the explicit configuration argument.
        fn fingerprint_topological_torsion_sparse_count_legacy_with_params(
            &self,
            py: Python<'_>,
            params: &LegacyTopologicalTorsionParams,
        ) -> PyResult<SparseCountFingerprint> {
            self._fingerprint_topological_torsion_sparse_count_legacy_params(py, Some(params))
        }

        #[pyo3(signature=(torsion_atom_count=None))]
        /// Return legacy topological torsion identifiers with the requested torsion and chirality settings.
        fn topological_torsion_ids(
            &self,
            py: Python<'_>,
            torsion_atom_count: Option<u32>,
        ) -> PyResult<Vec<u64>> {
            py.detach(|| match torsion_atom_count {
                Some(count) => self.inner.topological_torsion_ids_with_params(count),
                None => self.inner.topological_torsion_ids(),
            })
            .map_err(crate::morgan_errors::TopologicalTorsionReadError::from_rust)
        }
        /// Return legacy topological torsion identifiers with the requested torsion and chirality settings. This entry accepts the explicit native parameter form.
        fn topological_torsion_ids_with_params(
            &self,
            py: Python<'_>,
            torsion_atom_count: u32,
        ) -> PyResult<Vec<u64>> {
            py.detach(|| {
                self.inner
                    .topological_torsion_ids_with_params(torsion_atom_count)
            })
            .map_err(crate::morgan_errors::TopologicalTorsionReadError::from_rust)
        }
        #[pyo3(signature=(generator,params=crate::arguments::Argument::Missing,additional_output=None,*,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing,custom_bond_invariants=crate::arguments::Argument::Missing,use_legacy_stereo_perception=crate::arguments::Argument::Missing))]
        /// Calculate a fixed-width bit fingerprint from topological torsion paths with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_topological_torsion_with_generator(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalTorsionCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_bond_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            use_legacy_stereo_perception: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied()
                || custom_bond_invariants.is_supplied()
                || use_legacy_stereo_perception.is_supplied();
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, use_legacy_stereo_perception);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    use_legacy_stereo_perception,
                )
            } else {
                self._fingerprint_topological_torsion_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_with_generator"}"#]
        #[pyo3(signature=(generator,params=None,additional_output=None))]
        fn _fingerprint_topological_torsion_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: Option<&TopologicalTorsionCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            crate::torsion_reusable::fingerprint_topological_torsion_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_with_generator"}"#]
        #[pyo3(signature=(generator,additional_output=None,*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,use_legacy_stereo_perception=None))]
        fn _fingerprint_topological_torsion_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
            use_legacy_stereo_perception: Option<bool>,
        ) -> PyResult<Fingerprint> {
            let params = crate::torsion_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                use_legacy_stereo_perception,
            )?;
            self._fingerprint_topological_torsion_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator,params=crate::arguments::Argument::Missing,additional_output=None,*,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing,custom_bond_invariants=crate::arguments::Argument::Missing,use_legacy_stereo_perception=crate::arguments::Argument::Missing))]
        /// Calculate a fixed-width count fingerprint from topological torsion paths with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_topological_torsion_count_with_generator(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalTorsionCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_bond_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            use_legacy_stereo_perception: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<SparseCountFingerprint32> {
            let keyword_form = from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied()
                || custom_bond_invariants.is_supplied()
                || use_legacy_stereo_perception.is_supplied();
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, use_legacy_stereo_perception);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_count_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    use_legacy_stereo_perception,
                )
            } else {
                self._fingerprint_topological_torsion_count_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count_with_generator"}"#]
        #[pyo3(signature=(generator,params=None,additional_output=None))]
        fn _fingerprint_topological_torsion_count_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: Option<&TopologicalTorsionCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::torsion_reusable::fingerprint_topological_torsion_count_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count_with_generator"}"#]
        #[pyo3(signature=(generator,additional_output=None,*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,use_legacy_stereo_perception=None))]
        fn _fingerprint_topological_torsion_count_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
            use_legacy_stereo_perception: Option<bool>,
        ) -> PyResult<SparseCountFingerprint32> {
            let params = crate::torsion_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                use_legacy_stereo_perception,
            )?;
            self._fingerprint_topological_torsion_count_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator,params=crate::arguments::Argument::Missing,additional_output=None,*,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing,custom_bond_invariants=crate::arguments::Argument::Missing,use_legacy_stereo_perception=crate::arguments::Argument::Missing))]
        /// Calculate a sparse bit fingerprint from topological torsion paths with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_topological_torsion_sparse_with_generator(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalTorsionCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_bond_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            use_legacy_stereo_perception: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<SparseBitFingerprint> {
            let keyword_form = from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied()
                || custom_bond_invariants.is_supplied()
                || use_legacy_stereo_perception.is_supplied();
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, use_legacy_stereo_perception);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_sparse_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    use_legacy_stereo_perception,
                )
            } else {
                self._fingerprint_topological_torsion_sparse_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_with_generator"}"#]
        #[pyo3(signature=(generator,params=None,additional_output=None))]
        fn _fingerprint_topological_torsion_sparse_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: Option<&TopologicalTorsionCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::torsion_reusable::fingerprint_topological_torsion_sparse_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_with_generator"}"#]
        #[pyo3(signature=(generator,additional_output=None,*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,use_legacy_stereo_perception=None))]
        fn _fingerprint_topological_torsion_sparse_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
            use_legacy_stereo_perception: Option<bool>,
        ) -> PyResult<SparseBitFingerprint> {
            let params = crate::torsion_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                use_legacy_stereo_perception,
            )?;
            self._fingerprint_topological_torsion_sparse_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator,params=crate::arguments::Argument::Missing,additional_output=None,*,from_atoms=crate::arguments::Argument::Missing,ignore_atoms=crate::arguments::Argument::Missing,custom_atom_invariants=crate::arguments::Argument::Missing,custom_bond_invariants=crate::arguments::Argument::Missing,use_legacy_stereo_perception=crate::arguments::Argument::Missing))]
        /// Calculate a sparse count fingerprint from topological torsion paths with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_topological_torsion_sparse_count_with_generator(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalTorsionCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            custom_bond_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            use_legacy_stereo_perception: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<SparseCountFingerprint> {
            let keyword_form = from_atoms.is_supplied()
                || ignore_atoms.is_supplied()
                || custom_atom_invariants.is_supplied()
                || custom_bond_invariants.is_supplied()
                || use_legacy_stereo_perception.is_supplied();
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, use_legacy_stereo_perception);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_torsion_sparse_count_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    use_legacy_stereo_perception,
                )
            } else {
                self._fingerprint_topological_torsion_sparse_count_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count_with_generator"}"#]
        #[pyo3(signature=(generator,params=None,additional_output=None))]
        fn _fingerprint_topological_torsion_sparse_count_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            params: Option<&TopologicalTorsionCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::torsion_reusable::fingerprint_topological_torsion_sparse_count_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count_with_generator"}"#]
        #[pyo3(signature=(generator,additional_output=None,*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,use_legacy_stereo_perception=None))]
        fn _fingerprint_topological_torsion_sparse_count_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &TopologicalTorsionFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
            use_legacy_stereo_perception: Option<bool>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::torsion_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                use_legacy_stereo_perception,
            )?;
            self._fingerprint_topological_torsion_sparse_count_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }
        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, radius=Argument::Missing, include_chirality=Argument::Missing, use_bond_types=Argument::Missing, include_ring_membership=Argument::Missing, only_nonzero_invariants=Argument::Missing, include_redundant_environments=Argument::Missing, fp_size=Argument::Missing, count_simulation=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, invariants=Argument::Missing))]
        /// Calculate a fixed-width bit fingerprint from Morgan circular atom environments.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_morgan(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, MorganFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_bond_types: Argument<FieldInput<Scalar<bool>>>,
            include_ring_membership: Argument<FieldInput<Scalar<bool>>>,
            only_nonzero_invariants: Argument<FieldInput<Scalar<bool>>>,
            include_redundant_environments: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            invariants: Argument<FieldInput<ConfigurationField<MorganInvariants>>>,
        ) -> PyResult<Fingerprint> {
            crate::arguments::configuration_arguments!(params; radius, include_chirality, use_bond_types, include_ring_membership, only_nonzero_invariants, include_redundant_environments, fp_size, count_simulation, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, invariants);
            if let Some(params) = params.as_deref() {
                self._fingerprint_morgan_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_morgan_keywords(
                    py,
                    additional_output,
                    radius,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fp_size,
                    count_simulation,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    invariants,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_morgan_params(
            &self,
            py: Python<'_>,
            params: Option<&MorganFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_morgan(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan"}"#]
        #[pyo3(signature=(additional_output=None, *, radius=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_bond_types=FieldInput::<Scalar<bool>>::Default, include_ring_membership=FieldInput::<Scalar<bool>>::Default, only_nonzero_invariants=FieldInput::<Scalar<bool>>::Default, include_redundant_environments=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, invariants=FieldInput::<ConfigurationField<MorganInvariants>>::Default))]
        fn _fingerprint_morgan_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_bond_types: FieldInput<Scalar<bool>>,
            include_ring_membership: FieldInput<Scalar<bool>>,
            only_nonzero_invariants: FieldInput<Scalar<bool>>,
            include_redundant_environments: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            count_simulation: FieldInput<Scalar<bool>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            invariants: FieldInput<ConfigurationField<MorganInvariants>>,
        ) -> PyResult<Fingerprint> {
            let params = crate::fingerprint_calls::morgan_keywords(
                radius,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
                fp_size,
                count_simulation,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                invariants,
            );
            crate::fingerprint_calls::run_morgan(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width bit fingerprint from Morgan circular atom environments using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_morgan_with_params(
            &self,
            py: Python<'_>,
            params: &MorganFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            self._fingerprint_morgan_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, radius=Argument::Missing, include_chirality=Argument::Missing, use_bond_types=Argument::Missing, include_ring_membership=Argument::Missing, only_nonzero_invariants=Argument::Missing, include_redundant_environments=Argument::Missing, fp_size=Argument::Missing, count_simulation=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, invariants=Argument::Missing))]
        /// Calculate a fixed-width count fingerprint from Morgan circular atom environments.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_morgan_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, MorganFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_bond_types: Argument<FieldInput<Scalar<bool>>>,
            include_ring_membership: Argument<FieldInput<Scalar<bool>>>,
            only_nonzero_invariants: Argument<FieldInput<Scalar<bool>>>,
            include_redundant_environments: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            invariants: Argument<FieldInput<ConfigurationField<MorganInvariants>>>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::arguments::configuration_arguments!(params; radius, include_chirality, use_bond_types, include_ring_membership, only_nonzero_invariants, include_redundant_environments, fp_size, count_simulation, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, invariants);
            if let Some(params) = params.as_deref() {
                self._fingerprint_morgan_count_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_morgan_count_keywords(
                    py,
                    additional_output,
                    radius,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fp_size,
                    count_simulation,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    invariants,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_morgan_count_params(
            &self,
            py: Python<'_>,
            params: Option<&MorganFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_morgan_count(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_count"}"#]
        #[pyo3(signature=(additional_output=None, *, radius=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_bond_types=FieldInput::<Scalar<bool>>::Default, include_ring_membership=FieldInput::<Scalar<bool>>::Default, only_nonzero_invariants=FieldInput::<Scalar<bool>>::Default, include_redundant_environments=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, invariants=FieldInput::<ConfigurationField<MorganInvariants>>::Default))]
        fn _fingerprint_morgan_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_bond_types: FieldInput<Scalar<bool>>,
            include_ring_membership: FieldInput<Scalar<bool>>,
            only_nonzero_invariants: FieldInput<Scalar<bool>>,
            include_redundant_environments: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            count_simulation: FieldInput<Scalar<bool>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            invariants: FieldInput<ConfigurationField<MorganInvariants>>,
        ) -> PyResult<SparseCountFingerprint32> {
            let params = crate::fingerprint_calls::morgan_keywords(
                radius,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
                fp_size,
                count_simulation,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                invariants,
            );
            crate::fingerprint_calls::run_morgan_count(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width count fingerprint from Morgan circular atom environments using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_morgan_count_with_params(
            &self,
            py: Python<'_>,
            params: &MorganFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            self._fingerprint_morgan_count_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, radius=Argument::Missing, include_chirality=Argument::Missing, use_bond_types=Argument::Missing, include_ring_membership=Argument::Missing, only_nonzero_invariants=Argument::Missing, include_redundant_environments=Argument::Missing, fp_size=Argument::Missing, count_simulation=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, invariants=Argument::Missing))]
        /// Calculate a sparse bit fingerprint from Morgan circular atom environments.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_morgan_sparse(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, MorganFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_bond_types: Argument<FieldInput<Scalar<bool>>>,
            include_ring_membership: Argument<FieldInput<Scalar<bool>>>,
            only_nonzero_invariants: Argument<FieldInput<Scalar<bool>>>,
            include_redundant_environments: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            invariants: Argument<FieldInput<ConfigurationField<MorganInvariants>>>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::arguments::configuration_arguments!(params; radius, include_chirality, use_bond_types, include_ring_membership, only_nonzero_invariants, include_redundant_environments, fp_size, count_simulation, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, invariants);
            if let Some(params) = params.as_deref() {
                self._fingerprint_morgan_sparse_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_morgan_sparse_keywords(
                    py,
                    additional_output,
                    radius,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fp_size,
                    count_simulation,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    invariants,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_morgan_sparse_params(
            &self,
            py: Python<'_>,
            params: Option<&MorganFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_morgan_sparse(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse"}"#]
        #[pyo3(signature=(additional_output=None, *, radius=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_bond_types=FieldInput::<Scalar<bool>>::Default, include_ring_membership=FieldInput::<Scalar<bool>>::Default, only_nonzero_invariants=FieldInput::<Scalar<bool>>::Default, include_redundant_environments=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, invariants=FieldInput::<ConfigurationField<MorganInvariants>>::Default))]
        fn _fingerprint_morgan_sparse_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_bond_types: FieldInput<Scalar<bool>>,
            include_ring_membership: FieldInput<Scalar<bool>>,
            only_nonzero_invariants: FieldInput<Scalar<bool>>,
            include_redundant_environments: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            count_simulation: FieldInput<Scalar<bool>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            invariants: FieldInput<ConfigurationField<MorganInvariants>>,
        ) -> PyResult<SparseBitFingerprint> {
            let params = crate::fingerprint_calls::morgan_keywords(
                radius,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
                fp_size,
                count_simulation,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                invariants,
            );
            crate::fingerprint_calls::run_morgan_sparse(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse bit fingerprint from Morgan circular atom environments using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_morgan_sparse_with_params(
            &self,
            py: Python<'_>,
            params: &MorganFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            self._fingerprint_morgan_sparse_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, radius=Argument::Missing, include_chirality=Argument::Missing, use_bond_types=Argument::Missing, include_ring_membership=Argument::Missing, only_nonzero_invariants=Argument::Missing, include_redundant_environments=Argument::Missing, fp_size=Argument::Missing, count_simulation=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, invariants=Argument::Missing))]
        /// Calculate a sparse count fingerprint from Morgan circular atom environments.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_morgan_sparse_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, MorganFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_bond_types: Argument<FieldInput<Scalar<bool>>>,
            include_ring_membership: Argument<FieldInput<Scalar<bool>>>,
            only_nonzero_invariants: Argument<FieldInput<Scalar<bool>>>,
            include_redundant_environments: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            invariants: Argument<FieldInput<ConfigurationField<MorganInvariants>>>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::arguments::configuration_arguments!(params; radius, include_chirality, use_bond_types, include_ring_membership, only_nonzero_invariants, include_redundant_environments, fp_size, count_simulation, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, invariants);
            if let Some(params) = params.as_deref() {
                self._fingerprint_morgan_sparse_count_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_morgan_sparse_count_keywords(
                    py,
                    additional_output,
                    radius,
                    include_chirality,
                    use_bond_types,
                    include_ring_membership,
                    only_nonzero_invariants,
                    include_redundant_environments,
                    fp_size,
                    count_simulation,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    invariants,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_morgan_sparse_count_params(
            &self,
            py: Python<'_>,
            params: Option<&MorganFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_morgan_sparse_count(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_count"}"#]
        #[pyo3(signature=(additional_output=None, *, radius=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_bond_types=FieldInput::<Scalar<bool>>::Default, include_ring_membership=FieldInput::<Scalar<bool>>::Default, only_nonzero_invariants=FieldInput::<Scalar<bool>>::Default, include_redundant_environments=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, invariants=FieldInput::<ConfigurationField<MorganInvariants>>::Default))]
        fn _fingerprint_morgan_sparse_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            radius: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_bond_types: FieldInput<Scalar<bool>>,
            include_ring_membership: FieldInput<Scalar<bool>>,
            only_nonzero_invariants: FieldInput<Scalar<bool>>,
            include_redundant_environments: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            count_simulation: FieldInput<Scalar<bool>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            invariants: FieldInput<ConfigurationField<MorganInvariants>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::fingerprint_calls::morgan_keywords(
                radius,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
                fp_size,
                count_simulation,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                invariants,
            );
            crate::fingerprint_calls::run_morgan_sparse_count(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse count fingerprint from Morgan circular atom environments using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_morgan_sparse_count_with_params(
            &self,
            py: Python<'_>,
            params: &MorganFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            self._fingerprint_morgan_sparse_count_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, min_distance=Argument::Missing, max_distance=Argument::Missing, include_chirality=Argument::Missing, use_2d=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, conformer_id=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a fixed-width bit fingerprint from atom-pair atom types and distances.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_atom_pair(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, AtomPairFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: Argument<FieldInput<Scalar<u32>>>,
            max_distance: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_2d: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            conformer_id: Argument<FieldInput<Option<Scalar<usize>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<Fingerprint> {
            crate::arguments::configuration_arguments!(params; min_distance, max_distance, include_chirality, use_2d, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, conformer_id, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_atom_pair_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_atom_pair_keywords(
                    py,
                    additional_output,
                    min_distance,
                    max_distance,
                    include_chirality,
                    use_2d,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    conformer_id,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_atom_pair_params(
            &self,
            py: Python<'_>,
            params: Option<&AtomPairFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_atom_pair(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair"}"#]
        #[pyo3(signature=(additional_output=None, *, min_distance=FieldInput::<Scalar<u32>>::Default, max_distance=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_2d=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, conformer_id=FieldInput::<Option<Scalar<usize>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_atom_pair_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: FieldInput<Scalar<u32>>,
            max_distance: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_2d: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            conformer_id: FieldInput<Option<Scalar<usize>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<Fingerprint> {
            let params = crate::fingerprint_calls::atom_pair_keywords(
                min_distance,
                max_distance,
                include_chirality,
                use_2d,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_atom_pair(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width bit fingerprint from atom-pair atom types and distances using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_atom_pair_with_params(
            &self,
            py: Python<'_>,
            params: &AtomPairFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            self._fingerprint_atom_pair_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, min_distance=Argument::Missing, max_distance=Argument::Missing, include_chirality=Argument::Missing, use_2d=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, conformer_id=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a fixed-width count fingerprint from atom-pair atom types and distances.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_atom_pair_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, AtomPairFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: Argument<FieldInput<Scalar<u32>>>,
            max_distance: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_2d: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            conformer_id: Argument<FieldInput<Option<Scalar<usize>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::arguments::configuration_arguments!(params; min_distance, max_distance, include_chirality, use_2d, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, conformer_id, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_atom_pair_count_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_atom_pair_count_keywords(
                    py,
                    additional_output,
                    min_distance,
                    max_distance,
                    include_chirality,
                    use_2d,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    conformer_id,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_atom_pair_count_params(
            &self,
            py: Python<'_>,
            params: Option<&AtomPairFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_atom_pair_count(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_count"}"#]
        #[pyo3(signature=(additional_output=None, *, min_distance=FieldInput::<Scalar<u32>>::Default, max_distance=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_2d=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, conformer_id=FieldInput::<Option<Scalar<usize>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_atom_pair_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: FieldInput<Scalar<u32>>,
            max_distance: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_2d: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            conformer_id: FieldInput<Option<Scalar<usize>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseCountFingerprint32> {
            let params = crate::fingerprint_calls::atom_pair_keywords(
                min_distance,
                max_distance,
                include_chirality,
                use_2d,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_atom_pair_count(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width count fingerprint from atom-pair atom types and distances using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_atom_pair_count_with_params(
            &self,
            py: Python<'_>,
            params: &AtomPairFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            self._fingerprint_atom_pair_count_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, min_distance=Argument::Missing, max_distance=Argument::Missing, include_chirality=Argument::Missing, use_2d=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, conformer_id=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a sparse bit fingerprint from atom-pair atom types and distances.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_atom_pair_sparse(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, AtomPairFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: Argument<FieldInput<Scalar<u32>>>,
            max_distance: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_2d: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            conformer_id: Argument<FieldInput<Option<Scalar<usize>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::arguments::configuration_arguments!(params; min_distance, max_distance, include_chirality, use_2d, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, conformer_id, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_atom_pair_sparse_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_atom_pair_sparse_keywords(
                    py,
                    additional_output,
                    min_distance,
                    max_distance,
                    include_chirality,
                    use_2d,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    conformer_id,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_sparse"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_atom_pair_sparse_params(
            &self,
            py: Python<'_>,
            params: Option<&AtomPairFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_atom_pair_sparse(self, py, &effective, additional_output)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_sparse"}"#]
        #[pyo3(signature=(additional_output=None, *, min_distance=FieldInput::<Scalar<u32>>::Default, max_distance=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_2d=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, conformer_id=FieldInput::<Option<Scalar<usize>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_atom_pair_sparse_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: FieldInput<Scalar<u32>>,
            max_distance: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_2d: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            conformer_id: FieldInput<Option<Scalar<usize>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseBitFingerprint> {
            let params = crate::fingerprint_calls::atom_pair_keywords(
                min_distance,
                max_distance,
                include_chirality,
                use_2d,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_atom_pair_sparse(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse bit fingerprint from atom-pair atom types and distances using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_atom_pair_sparse_with_params(
            &self,
            py: Python<'_>,
            params: &AtomPairFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            self._fingerprint_atom_pair_sparse_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, min_distance=Argument::Missing, max_distance=Argument::Missing, include_chirality=Argument::Missing, use_2d=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, conformer_id=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a sparse count fingerprint from atom-pair atom types and distances.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_atom_pair_sparse_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, AtomPairFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: Argument<FieldInput<Scalar<u32>>>,
            max_distance: Argument<FieldInput<Scalar<u32>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            use_2d: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            conformer_id: Argument<FieldInput<Option<Scalar<usize>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::arguments::configuration_arguments!(params; min_distance, max_distance, include_chirality, use_2d, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, conformer_id, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_atom_pair_sparse_count_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_atom_pair_sparse_count_keywords(
                    py,
                    additional_output,
                    min_distance,
                    max_distance,
                    include_chirality,
                    use_2d,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    conformer_id,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_sparse_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_atom_pair_sparse_count_params(
            &self,
            py: Python<'_>,
            params: Option<&AtomPairFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_atom_pair_sparse_count(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_atom_pair_sparse_count"}"#]
        #[pyo3(signature=(additional_output=None, *, min_distance=FieldInput::<Scalar<u32>>::Default, max_distance=FieldInput::<Scalar<u32>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, use_2d=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, conformer_id=FieldInput::<Option<Scalar<usize>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_atom_pair_sparse_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            min_distance: FieldInput<Scalar<u32>>,
            max_distance: FieldInput<Scalar<u32>>,
            include_chirality: FieldInput<Scalar<bool>>,
            use_2d: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            conformer_id: FieldInput<Option<Scalar<usize>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::fingerprint_calls::atom_pair_keywords(
                min_distance,
                max_distance,
                include_chirality,
                use_2d,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_atom_pair_sparse_count(
                self,
                py,
                &params,
                additional_output,
            )
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse count fingerprint from atom-pair atom types and distances using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_atom_pair_sparse_count_with_params(
            &self,
            py: Python<'_>,
            params: &AtomPairFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            self._fingerprint_atom_pair_sparse_count_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, torsion_atom_count=Argument::Missing, only_shortest_paths=Argument::Missing, include_chirality=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a fixed-width bit fingerprint from topological torsion paths.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_topological_torsion(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, TopologicalTorsionFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: Argument<FieldInput<Scalar<u32>>>,
            only_shortest_paths: Argument<FieldInput<Scalar<bool>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<Fingerprint> {
            crate::arguments::configuration_arguments!(params; torsion_atom_count, only_shortest_paths, include_chirality, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_topological_torsion_params(py, Some(params), additional_output)
            } else {
                self._fingerprint_topological_torsion_keywords(
                    py,
                    additional_output,
                    torsion_atom_count,
                    only_shortest_paths,
                    include_chirality,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_topological_torsion_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalTorsionFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_topological_torsion(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion"}"#]
        #[pyo3(signature=(additional_output=None, *, torsion_atom_count=FieldInput::<Scalar<u32>>::Default, only_shortest_paths=FieldInput::<Scalar<bool>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_topological_torsion_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: FieldInput<Scalar<u32>>,
            only_shortest_paths: FieldInput<Scalar<bool>>,
            include_chirality: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<Fingerprint> {
            let params = crate::fingerprint_calls::topological_torsion_keywords(
                torsion_atom_count,
                only_shortest_paths,
                include_chirality,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_topological_torsion(self, py, &params, additional_output)
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width bit fingerprint from topological torsion paths using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_topological_torsion_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalTorsionFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            self._fingerprint_topological_torsion_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, torsion_atom_count=Argument::Missing, only_shortest_paths=Argument::Missing, include_chirality=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a fixed-width count fingerprint from topological torsion paths.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_topological_torsion_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, TopologicalTorsionFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: Argument<FieldInput<Scalar<u32>>>,
            only_shortest_paths: Argument<FieldInput<Scalar<bool>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::arguments::configuration_arguments!(params; torsion_atom_count, only_shortest_paths, include_chirality, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_topological_torsion_count_params(
                    py,
                    Some(params),
                    additional_output,
                )
            } else {
                self._fingerprint_topological_torsion_count_keywords(
                    py,
                    additional_output,
                    torsion_atom_count,
                    only_shortest_paths,
                    include_chirality,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_topological_torsion_count_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalTorsionFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_topological_torsion_count(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_count"}"#]
        #[pyo3(signature=(additional_output=None, *, torsion_atom_count=FieldInput::<Scalar<u32>>::Default, only_shortest_paths=FieldInput::<Scalar<bool>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_topological_torsion_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: FieldInput<Scalar<u32>>,
            only_shortest_paths: FieldInput<Scalar<bool>>,
            include_chirality: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseCountFingerprint32> {
            let params = crate::fingerprint_calls::topological_torsion_keywords(
                torsion_atom_count,
                only_shortest_paths,
                include_chirality,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_topological_torsion_count(
                self,
                py,
                &params,
                additional_output,
            )
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a fixed-width count fingerprint from topological torsion paths using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_topological_torsion_count_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalTorsionFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            self._fingerprint_topological_torsion_count_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, torsion_atom_count=Argument::Missing, only_shortest_paths=Argument::Missing, include_chirality=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a sparse bit fingerprint from topological torsion paths.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_topological_torsion_sparse(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, TopologicalTorsionFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: Argument<FieldInput<Scalar<u32>>>,
            only_shortest_paths: Argument<FieldInput<Scalar<bool>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::arguments::configuration_arguments!(params; torsion_atom_count, only_shortest_paths, include_chirality, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_topological_torsion_sparse_params(
                    py,
                    Some(params),
                    additional_output,
                )
            } else {
                self._fingerprint_topological_torsion_sparse_keywords(
                    py,
                    additional_output,
                    torsion_atom_count,
                    only_shortest_paths,
                    include_chirality,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_topological_torsion_sparse_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalTorsionFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_topological_torsion_sparse(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse"}"#]
        #[pyo3(signature=(additional_output=None, *, torsion_atom_count=FieldInput::<Scalar<u32>>::Default, only_shortest_paths=FieldInput::<Scalar<bool>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_topological_torsion_sparse_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: FieldInput<Scalar<u32>>,
            only_shortest_paths: FieldInput<Scalar<bool>>,
            include_chirality: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseBitFingerprint> {
            let params = crate::fingerprint_calls::topological_torsion_keywords(
                torsion_atom_count,
                only_shortest_paths,
                include_chirality,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_topological_torsion_sparse(
                self,
                py,
                &params,
                additional_output,
            )
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse bit fingerprint from topological torsion paths using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_topological_torsion_sparse_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalTorsionFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            self._fingerprint_topological_torsion_sparse_params(py, Some(params), additional_output)
        }

        #[pyo3(signature=(params=Argument::Missing, additional_output=None, *, torsion_atom_count=Argument::Missing, only_shortest_paths=Argument::Missing, include_chirality=Argument::Missing, count_simulation=Argument::Missing, fp_size=Argument::Missing, bits_per_feature=Argument::Missing, count_bounds=Argument::Missing, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing, atom_invariants_generator=Argument::Missing, use_legacy_stereo_perception=Argument::Missing))]
        /// Calculate a sparse count fingerprint from topological torsion paths.
        /// Supply a parameter object or keyword options, not both. additional_output,
        /// when supplied, collects the requested provenance in that output object.
        fn fingerprint_topological_torsion_sparse_count(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, TopologicalTorsionFingerprintParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: Argument<FieldInput<Scalar<u32>>>,
            only_shortest_paths: Argument<FieldInput<Scalar<bool>>>,
            include_chirality: Argument<FieldInput<Scalar<bool>>>,
            count_simulation: Argument<FieldInput<Scalar<bool>>>,
            fp_size: Argument<FieldInput<Scalar<u32>>>,
            bits_per_feature: Argument<FieldInput<Scalar<u32>>>,
            count_bounds: Argument<FieldInput<Vec<u32>>>,
            from_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            ignore_atoms: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_atom_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            custom_bond_invariants: Argument<FieldInput<Option<Vec<u32>>>>,
            atom_invariants_generator: Argument<
                FieldInput<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>,
            >,
            use_legacy_stereo_perception: Argument<FieldInput<Scalar<bool>>>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::arguments::configuration_arguments!(params; torsion_atom_count, only_shortest_paths, include_chirality, count_simulation, fp_size, bits_per_feature, count_bounds, from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants, atom_invariants_generator, use_legacy_stereo_perception);
            if let Some(params) = params.as_deref() {
                self._fingerprint_topological_torsion_sparse_count_params(
                    py,
                    Some(params),
                    additional_output,
                )
            } else {
                self._fingerprint_topological_torsion_sparse_count_keywords(
                    py,
                    additional_output,
                    torsion_atom_count,
                    only_shortest_paths,
                    include_chirality,
                    count_simulation,
                    fp_size,
                    bits_per_feature,
                    count_bounds,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                    atom_invariants_generator,
                    use_legacy_stereo_perception,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count"}"#]
        #[pyo3(signature=(params=None, additional_output=None))]
        fn _fingerprint_topological_torsion_sparse_count_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalTorsionFingerprintParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            let effective = params
                .map(|value| value.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            crate::fingerprint_calls::run_topological_torsion_sparse_count(
                self,
                py,
                &effective,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_torsion_sparse_count"}"#]
        #[pyo3(signature=(additional_output=None, *, torsion_atom_count=FieldInput::<Scalar<u32>>::Default, only_shortest_paths=FieldInput::<Scalar<bool>>::Default, include_chirality=FieldInput::<Scalar<bool>>::Default, count_simulation=FieldInput::<Scalar<bool>>::Default, fp_size=FieldInput::<Scalar<u32>>::Default, bits_per_feature=FieldInput::<Scalar<u32>>::Default, count_bounds=FieldInput::<Vec<u32>>::Default, from_atoms=FieldInput::<Option<Vec<u32>>>::Default, ignore_atoms=FieldInput::<Option<Vec<u32>>>::Default, custom_atom_invariants=FieldInput::<Option<Vec<u32>>>::Default, custom_bond_invariants=FieldInput::<Option<Vec<u32>>>::Default, atom_invariants_generator=FieldInput::<Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>>::Default, use_legacy_stereo_perception=FieldInput::<Scalar<bool>>::Default))]
        fn _fingerprint_topological_torsion_sparse_count_keywords(
            &self,
            py: Python<'_>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            torsion_atom_count: FieldInput<Scalar<u32>>,
            only_shortest_paths: FieldInput<Scalar<bool>>,
            include_chirality: FieldInput<Scalar<bool>>,
            count_simulation: FieldInput<Scalar<bool>>,
            fp_size: FieldInput<Scalar<u32>>,
            bits_per_feature: FieldInput<Scalar<u32>>,
            count_bounds: FieldInput<Vec<u32>>,
            from_atoms: FieldInput<Option<Vec<u32>>>,
            ignore_atoms: FieldInput<Option<Vec<u32>>>,
            custom_atom_invariants: FieldInput<Option<Vec<u32>>>,
            custom_bond_invariants: FieldInput<Option<Vec<u32>>>,
            atom_invariants_generator: FieldInput<
                Option<ConfigurationField<AtomPairAtomInvariantsGenerator>>,
            >,
            use_legacy_stereo_perception: FieldInput<Scalar<bool>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::fingerprint_calls::topological_torsion_keywords(
                torsion_atom_count,
                only_shortest_paths,
                include_chirality,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                atom_invariants_generator,
                use_legacy_stereo_perception,
            );
            crate::fingerprint_calls::run_topological_torsion_sparse_count(
                self,
                py,
                &params,
                additional_output,
            )
        }
        #[pyo3(signature=(params, additional_output=None))]
        /// Calculate a sparse count fingerprint from topological torsion paths using the supplied parameter object.
        /// The input molecule is unchanged; additional_output collects requested provenance.
        fn fingerprint_topological_torsion_sparse_count_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalTorsionFingerprintParams,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            self._fingerprint_topological_torsion_sparse_count_params(
                py,
                Some(params),
                additional_output,
            )
        }
        #[pyo3(signature=(generator, params=Argument::Missing, additional_output=None, *, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing))]
        /// Calculate a fixed-width bit fingerprint from Morgan circular atom environments with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_morgan_with_generator(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Argument<Option<PyRef<'_, MorganCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Argument<Option<Vec<u32>>>,
            ignore_atoms: Argument<Option<Vec<u32>>>,
            custom_atom_invariants: Argument<Option<Vec<u32>>>,
            custom_bond_invariants: Argument<Option<Vec<u32>>>,
        ) -> PyResult<Fingerprint> {
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants);
            let params = params.as_deref();
            let configured = from_atoms.is_some()
                || ignore_atoms.is_some()
                || custom_atom_invariants.is_some()
                || custom_bond_invariants.is_some();
            if configured {
                self._fingerprint_morgan_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                )
            } else {
                self._fingerprint_morgan_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_with_generator"}"#]
        #[pyo3(signature=(generator, params=None, additional_output=None))]
        fn _fingerprint_morgan_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Option<&MorganCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<Fingerprint> {
            crate::morgan_reusable::fingerprint_morgan_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_with_generator"}"#]
        #[pyo3(signature=(generator, additional_output=None, *, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None))]
        fn _fingerprint_morgan_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
        ) -> PyResult<Fingerprint> {
            let params = crate::morgan_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
            )?;
            self._fingerprint_morgan_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator, params=Argument::Missing, additional_output=None, *, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing))]
        /// Calculate a fixed-width count fingerprint from Morgan circular atom environments with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_morgan_count_with_generator(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Argument<Option<PyRef<'_, MorganCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Argument<Option<Vec<u32>>>,
            ignore_atoms: Argument<Option<Vec<u32>>>,
            custom_atom_invariants: Argument<Option<Vec<u32>>>,
            custom_bond_invariants: Argument<Option<Vec<u32>>>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants);
            let params = params.as_deref();
            let configured = from_atoms.is_some()
                || ignore_atoms.is_some()
                || custom_atom_invariants.is_some()
                || custom_bond_invariants.is_some();
            if configured {
                self._fingerprint_morgan_count_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                )
            } else {
                self._fingerprint_morgan_count_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_count_with_generator"}"#]
        #[pyo3(signature=(generator, params=None, additional_output=None))]
        fn _fingerprint_morgan_count_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Option<&MorganCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint32> {
            crate::morgan_reusable::fingerprint_morgan_count_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_count_with_generator"}"#]
        #[pyo3(signature=(generator, additional_output=None, *, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None))]
        fn _fingerprint_morgan_count_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
        ) -> PyResult<SparseCountFingerprint32> {
            let params = crate::morgan_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
            )?;
            self._fingerprint_morgan_count_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator, params=Argument::Missing, additional_output=None, *, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing))]
        /// Calculate a sparse bit fingerprint from Morgan circular atom environments with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_morgan_sparse_with_generator(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Argument<Option<PyRef<'_, MorganCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Argument<Option<Vec<u32>>>,
            ignore_atoms: Argument<Option<Vec<u32>>>,
            custom_atom_invariants: Argument<Option<Vec<u32>>>,
            custom_bond_invariants: Argument<Option<Vec<u32>>>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants);
            let params = params.as_deref();
            let configured = from_atoms.is_some()
                || ignore_atoms.is_some()
                || custom_atom_invariants.is_some()
                || custom_bond_invariants.is_some();
            if configured {
                self._fingerprint_morgan_sparse_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                )
            } else {
                self._fingerprint_morgan_sparse_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_with_generator"}"#]
        #[pyo3(signature=(generator, params=None, additional_output=None))]
        fn _fingerprint_morgan_sparse_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Option<&MorganCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseBitFingerprint> {
            crate::morgan_reusable::fingerprint_morgan_sparse_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_with_generator"}"#]
        #[pyo3(signature=(generator, additional_output=None, *, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None))]
        fn _fingerprint_morgan_sparse_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
        ) -> PyResult<SparseBitFingerprint> {
            let params = crate::morgan_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
            )?;
            self._fingerprint_morgan_sparse_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(generator, params=Argument::Missing, additional_output=None, *, from_atoms=Argument::Missing, ignore_atoms=Argument::Missing, custom_atom_invariants=Argument::Missing, custom_bond_invariants=Argument::Missing))]
        /// Calculate a sparse count fingerprint from Morgan circular atom environments with a reusable generator.
        /// The generator supplies algorithm settings; per-call parameters select atoms and
        /// custom invariants. Supply per-call Params or keywords, not both.
        fn fingerprint_morgan_sparse_count_with_generator(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Argument<Option<PyRef<'_, MorganCallParams>>>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Argument<Option<Vec<u32>>>,
            ignore_atoms: Argument<Option<Vec<u32>>>,
            custom_atom_invariants: Argument<Option<Vec<u32>>>,
            custom_bond_invariants: Argument<Option<Vec<u32>>>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::arguments::configuration_arguments!(params; from_atoms, ignore_atoms, custom_atom_invariants, custom_bond_invariants);
            let params = params.as_deref();
            let configured = from_atoms.is_some()
                || ignore_atoms.is_some()
                || custom_atom_invariants.is_some()
                || custom_bond_invariants.is_some();
            if configured {
                self._fingerprint_morgan_sparse_count_with_generator_keywords(
                    py,
                    generator,
                    additional_output,
                    from_atoms,
                    ignore_atoms,
                    custom_atom_invariants,
                    custom_bond_invariants,
                )
            } else {
                self._fingerprint_morgan_sparse_count_with_generator_params(
                    py,
                    generator,
                    params,
                    additional_output,
                )
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_count_with_generator"}"#]
        #[pyo3(signature=(generator, params=None, additional_output=None))]
        fn _fingerprint_morgan_sparse_count_with_generator_params(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            params: Option<&MorganCallParams>,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
        ) -> PyResult<SparseCountFingerprint> {
            crate::morgan_reusable::fingerprint_morgan_sparse_count_with_generator(
                self,
                py,
                generator,
                params,
                additional_output,
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_morgan_sparse_count_with_generator"}"#]
        #[pyo3(signature=(generator, additional_output=None, *, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None))]
        fn _fingerprint_morgan_sparse_count_with_generator_keywords(
            &self,
            py: Python<'_>,
            generator: &MorganFingerprintGenerator,
            additional_output: Option<&mut FingerprintAdditionalOutput>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
            custom_atom_invariants: Option<Vec<u32>>,
            custom_bond_invariants: Option<Vec<u32>>,
        ) -> PyResult<SparseCountFingerprint> {
            let params = crate::morgan_reusable::call_keywords(
                py,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
            )?;
            self._fingerprint_morgan_sparse_count_with_generator_params(
                py,
                generator,
                Some(&params),
                additional_output,
            )
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, include_atom_maps=crate::arguments::Argument::Missing, isomeric_smiles=crate::arguments::Argument::Missing, include_dative_bonds=crate::arguments::Argument::Missing, rooted_at_atom=crate::arguments::Argument::Missing))]
        /// Serialize the molecule as SMARTS without changing it. Supply SmartsWriteParams or keyword options, not both.
        fn to_smarts(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, SmartsWriteParams>>>,
            include_atom_maps: crate::arguments::Argument<Option<bool>>,
            isomeric_smiles: crate::arguments::Argument<Option<bool>>,
            include_dative_bonds: crate::arguments::Argument<Option<bool>>,
            rooted_at_atom: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<String> {
            let keyword_form = include_atom_maps.is_supplied()
                || isomeric_smiles.is_supplied()
                || include_dative_bonds.is_supplied()
                || rooted_at_atom.is_supplied();
            crate::arguments::configuration_arguments!(params; include_atom_maps, isomeric_smiles, include_dative_bonds, rooted_at_atom);
            let params = params.as_deref();
            if keyword_form {
                self._to_smarts_keywords(
                    include_atom_maps,
                    isomeric_smiles,
                    include_dative_bonds,
                    rooted_at_atom,
                )
            } else {
                self._to_smarts_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_smarts"}"#]
        #[pyo3(signature=(params=None))]
        fn _to_smarts_params(&self, params: Option<&SmartsWriteParams>) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
            let defaults = ck::SmartsWriteParams::default();
            self.inner
                .to_smarts_with_params(params.as_ref().unwrap_or(&defaults))
                .map_err(SmartsWriteError::from_rust)
                .and_then(crate::smarts_writer::text)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_smarts"}"#]
        #[pyo3(signature=(*, include_atom_maps=None, isomeric_smiles=None, include_dative_bonds=None, rooted_at_atom=None))]
        fn _to_smarts_keywords(
            &self,
            include_atom_maps: Option<bool>,
            isomeric_smiles: Option<bool>,
            include_dative_bonds: Option<bool>,
            rooted_at_atom: Option<usize>,
        ) -> PyResult<String> {
            self._to_smarts_params(Some(&crate::smarts_writer::keywords(
                include_atom_maps,
                isomeric_smiles,
                include_dative_bonds,
                rooted_at_atom,
            )))
        }
        /// Serialize the molecule as SMARTS without changing it. Supply SmartsWriteParams or keyword options, not both. This entry accepts the explicit native parameter form.
        fn to_smarts_with_params(&self, params: &SmartsWriteParams) -> PyResult<String> {
            self._to_smarts_params(Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, include_atom_maps=crate::arguments::Argument::Missing, isomeric_smiles=crate::arguments::Argument::Missing, include_dative_bonds=crate::arguments::Argument::Missing, rooted_at_atom=crate::arguments::Argument::Missing))]
        /// Serialize SMARTS with CX extensions without changing the molecule. Supply SmartsWriteParams or keyword options, not both.
        fn to_cx_smarts(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, SmartsWriteParams>>>,
            include_atom_maps: crate::arguments::Argument<Option<bool>>,
            isomeric_smiles: crate::arguments::Argument<Option<bool>>,
            include_dative_bonds: crate::arguments::Argument<Option<bool>>,
            rooted_at_atom: crate::arguments::Argument<Option<usize>>,
        ) -> PyResult<String> {
            let keyword_form = include_atom_maps.is_supplied()
                || isomeric_smiles.is_supplied()
                || include_dative_bonds.is_supplied()
                || rooted_at_atom.is_supplied();
            crate::arguments::configuration_arguments!(params; include_atom_maps, isomeric_smiles, include_dative_bonds, rooted_at_atom);
            let params = params.as_deref();
            if keyword_form {
                self._to_cx_smarts_keywords(
                    include_atom_maps,
                    isomeric_smiles,
                    include_dative_bonds,
                    rooted_at_atom,
                )
            } else {
                self._to_cx_smarts_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_cx_smarts"}"#]
        #[pyo3(signature=(params=None))]
        fn _to_cx_smarts_params(&self, params: Option<&SmartsWriteParams>) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
            let defaults = ck::SmartsWriteParams::default();
            self.inner
                .to_cx_smarts_with_params(params.as_ref().unwrap_or(&defaults))
                .map_err(SmartsWriteError::from_rust)
                .and_then(crate::smarts_writer::text)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_cx_smarts"}"#]
        #[pyo3(signature=(*, include_atom_maps=None, isomeric_smiles=None, include_dative_bonds=None, rooted_at_atom=None))]
        fn _to_cx_smarts_keywords(
            &self,
            include_atom_maps: Option<bool>,
            isomeric_smiles: Option<bool>,
            include_dative_bonds: Option<bool>,
            rooted_at_atom: Option<usize>,
        ) -> PyResult<String> {
            self._to_cx_smarts_params(Some(&crate::smarts_writer::keywords(
                include_atom_maps,
                isomeric_smiles,
                include_dative_bonds,
                rooted_at_atom,
            )))
        }
        /// Serialize SMARTS with CX extensions without changing the molecule. Supply SmartsWriteParams or keyword options, not both. This entry accepts the explicit native parameter form.
        fn to_cx_smarts_with_params(&self, params: &SmartsWriteParams) -> PyResult<String> {
            self._to_cx_smarts_params(Some(params))
        }

        #[pyo3(signature=(reference, params=crate::arguments::Argument::Missing, *, probe_conformer_id=crate::arguments::Argument::Missing, reference_conformer_id=crate::arguments::Argument::Missing, atom_map=crate::arguments::Argument::Missing, weights=crate::arguments::Argument::Missing, reflect=crate::arguments::Argument::Missing, max_iterations=crate::arguments::Argument::Missing))]
        /// Calculate a coordinate alignment transform to a reference without applying it. Parameters select conformers, atom correspondences, weights and reflection.
        fn alignment_transform_to(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: crate::arguments::Argument<Option<PyRef<'_, AlignmentParameters>>>,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_map: crate::arguments::Argument<Option<Vec<AlignmentAtomMap>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<AlignmentResult> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_map.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_map, weights, reflect, max_iterations);
            let params = params.as_deref();
            if keyword_form {
                self._alignment_transform_to_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_map,
                    weights,
                    reflect,
                    max_iterations,
                )
            } else {
                self._alignment_transform_to_params(py, reference, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.alignment_transform_to"}"#]
        #[pyo3(signature=(reference, params=None))]
        fn _alignment_transform_to_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&AlignmentParameters>,
        ) -> PyResult<AlignmentResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._alignment_transform_to_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.alignment_transform_to"}"#]
        #[pyo3(signature=(reference, *, probe_conformer_id=None, reference_conformer_id=None, atom_map=None, weights=None, reflect=None, max_iterations=None))]
        fn _alignment_transform_to_keywords(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_map: Option<Vec<AlignmentAtomMap>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
        ) -> PyResult<AlignmentResult> {
            let params = crate::alignment_queries::transform_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_map,
                weights,
                reflect,
                max_iterations,
            );
            self._alignment_transform_to_params_native(py, reference, Some(&params))
        }
        /// Calculate a coordinate alignment transform to a reference without applying it. Parameters select conformers, atom correspondences, weights and reflection. This entry accepts the explicit native parameter form.
        fn alignment_transform_to_with_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: &AlignmentParameters,
        ) -> PyResult<AlignmentResult> {
            self._alignment_transform_to_params(py, reference, Some(params))
        }

        #[pyo3(signature=(reference, params=crate::arguments::Argument::Missing, *, probe_conformer_id=crate::arguments::Argument::Missing, reference_conformer_id=crate::arguments::Argument::Missing, atom_maps=crate::arguments::Argument::Missing, weights=crate::arguments::Argument::Missing, reflect=crate::arguments::Argument::Missing, max_iterations=crate::arguments::Argument::Missing, max_matches=crate::arguments::Argument::Missing, symmetrize_conjugated_terminal_groups=crate::arguments::Argument::Missing, ignore_hydrogens=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing))]
        /// Find the best alignment to a reference across allowed atom mappings, returning its transform and RMSD without changing either molecule.
        fn best_alignment_to(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: crate::arguments::Argument<Option<PyRef<'_, BestAlignmentParameters>>>,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_maps: crate::arguments::Argument<Option<Vec<Vec<AlignmentAtomMap>>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            max_matches: crate::arguments::Argument<Option<i32>>,
            symmetrize_conjugated_terminal_groups: crate::arguments::Argument<Option<bool>>,
            ignore_hydrogens: crate::arguments::Argument<Option<bool>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
        ) -> PyResult<AlignmentResult> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_maps.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied()
                || max_matches.is_supplied()
                || symmetrize_conjugated_terminal_groups.is_supplied()
                || ignore_hydrogens.is_supplied()
                || num_threads.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_maps, weights, reflect, max_iterations, max_matches, symmetrize_conjugated_terminal_groups, ignore_hydrogens, num_threads);
            let params = params.as_deref();
            if keyword_form {
                self._best_alignment_to_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_maps,
                    weights,
                    reflect,
                    max_iterations,
                    max_matches,
                    symmetrize_conjugated_terminal_groups,
                    ignore_hydrogens,
                    num_threads,
                )
            } else {
                self._best_alignment_to_params(py, reference, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.best_alignment_to"}"#]
        #[pyo3(signature=(reference, params=None))]
        fn _best_alignment_to_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&BestAlignmentParameters>,
        ) -> PyResult<AlignmentResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._best_alignment_to_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.best_alignment_to"}"#]
        #[pyo3(signature=(reference, *, probe_conformer_id=None, reference_conformer_id=None, atom_maps=None, weights=None, reflect=None, max_iterations=None, max_matches=None, symmetrize_conjugated_terminal_groups=None, ignore_hydrogens=None, num_threads=None))]
        fn _best_alignment_to_keywords(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
            max_matches: Option<i32>,
            symmetrize_conjugated_terminal_groups: Option<bool>,
            ignore_hydrogens: Option<bool>,
            num_threads: Option<i32>,
        ) -> PyResult<AlignmentResult> {
            let params = crate::alignment_queries::best_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_maps,
                weights,
                reflect,
                max_iterations,
                max_matches,
                symmetrize_conjugated_terminal_groups,
                ignore_hydrogens,
                num_threads,
            );
            self._best_alignment_to_params_native(py, reference, Some(&params))
        }
        /// Find the best alignment to a reference across allowed atom mappings, returning its transform and RMSD without changing either molecule. This entry accepts the explicit native parameter form.
        fn best_alignment_to_with_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: &BestAlignmentParameters,
        ) -> PyResult<AlignmentResult> {
            self._best_alignment_to_params(py, reference, Some(params))
        }

        #[pyo3(signature=(reference, params=crate::arguments::Argument::Missing, *, probe_conformer_id=crate::arguments::Argument::Missing, reference_conformer_id=crate::arguments::Argument::Missing, atom_maps=crate::arguments::Argument::Missing, weights=crate::arguments::Argument::Missing, reflect=crate::arguments::Argument::Missing, max_iterations=crate::arguments::Argument::Missing, max_matches=crate::arguments::Argument::Missing, symmetrize_conjugated_terminal_groups=crate::arguments::Argument::Missing, ignore_hydrogens=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing))]
        /// Calculate the best alignment RMSD to a reference across allowed atom mappings without changing either molecule.
        fn best_rmsd_to(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: crate::arguments::Argument<Option<PyRef<'_, BestAlignmentParameters>>>,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_maps: crate::arguments::Argument<Option<Vec<Vec<AlignmentAtomMap>>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            max_matches: crate::arguments::Argument<Option<i32>>,
            symmetrize_conjugated_terminal_groups: crate::arguments::Argument<Option<bool>>,
            ignore_hydrogens: crate::arguments::Argument<Option<bool>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
        ) -> PyResult<f64> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_maps.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied()
                || max_matches.is_supplied()
                || symmetrize_conjugated_terminal_groups.is_supplied()
                || ignore_hydrogens.is_supplied()
                || num_threads.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_maps, weights, reflect, max_iterations, max_matches, symmetrize_conjugated_terminal_groups, ignore_hydrogens, num_threads);
            let params = params.as_deref();
            if keyword_form {
                self._best_rmsd_to_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_maps,
                    weights,
                    reflect,
                    max_iterations,
                    max_matches,
                    symmetrize_conjugated_terminal_groups,
                    ignore_hydrogens,
                    num_threads,
                )
            } else {
                self._best_rmsd_to_params(py, reference, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.best_rmsd_to"}"#]
        #[pyo3(signature=(reference, params=None))]
        fn _best_rmsd_to_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&BestAlignmentParameters>,
        ) -> PyResult<f64> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._best_rmsd_to_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.best_rmsd_to"}"#]
        #[pyo3(signature=(reference, *, probe_conformer_id=None, reference_conformer_id=None, atom_maps=None, weights=None, reflect=None, max_iterations=None, max_matches=None, symmetrize_conjugated_terminal_groups=None, ignore_hydrogens=None, num_threads=None))]
        fn _best_rmsd_to_keywords(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
            max_matches: Option<i32>,
            symmetrize_conjugated_terminal_groups: Option<bool>,
            ignore_hydrogens: Option<bool>,
            num_threads: Option<i32>,
        ) -> PyResult<f64> {
            let params = crate::alignment_queries::best_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_maps,
                weights,
                reflect,
                max_iterations,
                max_matches,
                symmetrize_conjugated_terminal_groups,
                ignore_hydrogens,
                num_threads,
            );
            self._best_rmsd_to_params_native(py, reference, Some(&params))
        }
        /// Calculate the best alignment RMSD to a reference across allowed atom mappings without changing either molecule. This entry accepts the explicit native parameter form.
        fn best_rmsd_to_with_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: &BestAlignmentParameters,
        ) -> PyResult<f64> {
            self._best_rmsd_to_params(py, reference, Some(params))
        }

        #[pyo3(signature=(reference, params=crate::arguments::Argument::Missing, *, probe_conformer_id=crate::arguments::Argument::Missing, reference_conformer_id=crate::arguments::Argument::Missing, atom_maps=crate::arguments::Argument::Missing, weights=crate::arguments::Argument::Missing, max_matches=crate::arguments::Argument::Missing, symmetrize_conjugated_terminal_groups=crate::arguments::Argument::Missing))]
        /// Calculate coordinate RMSD to the reference without first aligning the coordinates.
        fn coordinate_rmsd_to(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: crate::arguments::Argument<Option<PyRef<'_, CoordinateRmsdParameters>>>,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_maps: crate::arguments::Argument<Option<Vec<Vec<AlignmentAtomMap>>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            max_matches: crate::arguments::Argument<Option<i32>>,
            symmetrize_conjugated_terminal_groups: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<f64> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_maps.is_supplied()
                || weights.is_supplied()
                || max_matches.is_supplied()
                || symmetrize_conjugated_terminal_groups.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_maps, weights, max_matches, symmetrize_conjugated_terminal_groups);
            let params = params.as_deref();
            if keyword_form {
                self._coordinate_rmsd_to_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_maps,
                    weights,
                    max_matches,
                    symmetrize_conjugated_terminal_groups,
                )
            } else {
                self._coordinate_rmsd_to_params(py, reference, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.coordinate_rmsd_to"}"#]
        #[pyo3(signature=(reference, params=None))]
        fn _coordinate_rmsd_to_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: Option<&CoordinateRmsdParameters>,
        ) -> PyResult<f64> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._coordinate_rmsd_to_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.coordinate_rmsd_to"}"#]
        #[pyo3(signature=(reference, *, probe_conformer_id=None, reference_conformer_id=None, atom_maps=None, weights=None, max_matches=None, symmetrize_conjugated_terminal_groups=None))]
        fn _coordinate_rmsd_to_keywords(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
            weights: Option<Vec<f64>>,
            max_matches: Option<i32>,
            symmetrize_conjugated_terminal_groups: Option<bool>,
        ) -> PyResult<f64> {
            let params = crate::alignment_queries::coordinate_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_maps,
                weights,
                max_matches,
                symmetrize_conjugated_terminal_groups,
            );
            self._coordinate_rmsd_to_params_native(py, reference, Some(&params))
        }
        /// Calculate coordinate RMSD to the reference without first aligning the coordinates. This entry accepts the explicit native parameter form.
        fn coordinate_rmsd_to_with_params(
            &self,
            py: Python<'_>,
            reference: &Molecule,
            params: &CoordinateRmsdParameters,
        ) -> PyResult<f64> {
            self._coordinate_rmsd_to_params(py, reference, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, atom_maps=crate::arguments::Argument::Missing, weights=crate::arguments::Argument::Missing, max_matches=crate::arguments::Argument::Missing, symmetrize_conjugated_terminal_groups=crate::arguments::Argument::Missing, ignore_hydrogens=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing))]
        /// Calculate symmetry-aware best RMSDs between this molecule's conformers using the selected matching and alignment parameters.
        fn all_conformer_best_rmsds(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, AllConformerRmsdParameters>>>,
            atom_maps: crate::arguments::Argument<Option<Vec<Vec<AlignmentAtomMap>>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            max_matches: crate::arguments::Argument<Option<i32>>,
            symmetrize_conjugated_terminal_groups: crate::arguments::Argument<Option<bool>>,
            ignore_hydrogens: crate::arguments::Argument<Option<bool>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
        ) -> PyResult<Vec<ConformerRmsd>> {
            let keyword_form = atom_maps.is_supplied()
                || weights.is_supplied()
                || max_matches.is_supplied()
                || symmetrize_conjugated_terminal_groups.is_supplied()
                || ignore_hydrogens.is_supplied()
                || num_threads.is_supplied();
            crate::arguments::configuration_arguments!(params; atom_maps, weights, max_matches, symmetrize_conjugated_terminal_groups, ignore_hydrogens, num_threads);
            let params = params.as_deref();
            if keyword_form {
                self._all_conformer_best_rmsds_keywords(
                    py,
                    atom_maps,
                    weights,
                    max_matches,
                    symmetrize_conjugated_terminal_groups,
                    ignore_hydrogens,
                    num_threads,
                )
            } else {
                self._all_conformer_best_rmsds_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.all_conformer_best_rmsds"}"#]
        #[pyo3(signature=(params=None))]
        fn _all_conformer_best_rmsds_params(
            &self,
            py: Python<'_>,
            params: Option<&AllConformerRmsdParameters>,
        ) -> PyResult<Vec<ConformerRmsd>> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._all_conformer_best_rmsds_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.all_conformer_best_rmsds"}"#]
        #[pyo3(signature=(*, atom_maps=None, weights=None, max_matches=None, symmetrize_conjugated_terminal_groups=None, ignore_hydrogens=None, num_threads=None))]
        fn _all_conformer_best_rmsds_keywords(
            &self,
            py: Python<'_>,
            atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>,
            weights: Option<Vec<f64>>,
            max_matches: Option<i32>,
            symmetrize_conjugated_terminal_groups: Option<bool>,
            ignore_hydrogens: Option<bool>,
            num_threads: Option<i32>,
        ) -> PyResult<Vec<ConformerRmsd>> {
            let params = crate::alignment_queries::all_keywords(
                atom_maps,
                weights,
                max_matches,
                symmetrize_conjugated_terminal_groups,
                ignore_hydrogens,
                num_threads,
            );
            self._all_conformer_best_rmsds_params_native(py, Some(&params))
        }
        /// Calculate symmetry-aware best RMSDs between this molecule's conformers using the selected matching and alignment parameters. This entry accepts the explicit native parameter form.
        fn all_conformer_best_rmsds_with_params(
            &self,
            py: Python<'_>,
            params: &AllConformerRmsdParameters,
        ) -> PyResult<Vec<ConformerRmsd>> {
            self._all_conformer_best_rmsds_params(py, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, min_path=crate::arguments::Argument::Missing, max_path=crate::arguments::Argument::Missing, fp_size=crate::arguments::Argument::Missing, num_bits_per_feature=crate::arguments::Argument::Missing, use_hs=crate::arguments::Argument::Missing, target_density=crate::arguments::Argument::Missing, min_size=crate::arguments::Argument::Missing, branched_paths=crate::arguments::Argument::Missing, use_bond_order=crate::arguments::Argument::Missing, atom_invariants=crate::arguments::Argument::Missing, from_atoms=crate::arguments::Argument::Missing, ignore_atoms=crate::arguments::Argument::Missing))]
        /// Calculate the RDKit-style path-based topological bit fingerprint with the supplied path and folding settings.
        fn fingerprint_topological(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalFingerprintParams>>>,
            min_path: crate::arguments::Argument<Option<u32>>,
            max_path: crate::arguments::Argument<Option<u32>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            num_bits_per_feature: crate::arguments::Argument<Option<u32>>,
            use_hs: crate::arguments::Argument<Option<bool>>,
            target_density: crate::arguments::Argument<Option<f64>>,
            min_size: crate::arguments::Argument<Option<u32>>,
            branched_paths: crate::arguments::Argument<Option<bool>>,
            use_bond_order: crate::arguments::Argument<Option<bool>>,
            atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = min_path.is_supplied()
                || max_path.is_supplied()
                || fp_size.is_supplied()
                || num_bits_per_feature.is_supplied()
                || use_hs.is_supplied()
                || target_density.is_supplied()
                || min_size.is_supplied()
                || branched_paths.is_supplied()
                || use_bond_order.is_supplied()
                || atom_invariants.is_supplied()
                || from_atoms.is_supplied()
                || ignore_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; min_path, max_path, fp_size, num_bits_per_feature, use_hs, target_density, min_size, branched_paths, use_bond_order, atom_invariants, from_atoms, ignore_atoms);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_topological_keywords(
                    py,
                    min_path,
                    max_path,
                    fp_size,
                    num_bits_per_feature,
                    use_hs,
                    target_density,
                    min_size,
                    branched_paths,
                    use_bond_order,
                    atom_invariants,
                    from_atoms,
                    ignore_atoms,
                )
            } else {
                self._fingerprint_topological_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_topological_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological"}"#]
        #[pyo3(signature=(*, min_path=None, max_path=None, fp_size=None, num_bits_per_feature=None, use_hs=None, target_density=None, min_size=None, branched_paths=None, use_bond_order=None, atom_invariants=None, from_atoms=None, ignore_atoms=None))]
        fn _fingerprint_topological_keywords(
            &self,
            py: Python<'_>,
            min_path: Option<u32>,
            max_path: Option<u32>,
            fp_size: Option<u32>,
            num_bits_per_feature: Option<u32>,
            use_hs: Option<bool>,
            target_density: Option<f64>,
            min_size: Option<u32>,
            branched_paths: Option<bool>,
            use_bond_order: Option<bool>,
            atom_invariants: Option<Vec<u32>>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
        ) -> PyResult<Fingerprint> {
            let params = super::topological_producers::keyword_params(
                min_path,
                max_path,
                fp_size,
                num_bits_per_feature,
                use_hs,
                target_density,
                min_size,
                branched_paths,
                use_bond_order,
                atom_invariants,
                from_atoms,
                ignore_atoms,
            );
            self._fingerprint_topological_params_native(py, Some(&params))
        }
        /// Calculate the RDKit-style path-based topological bit fingerprint with the supplied path and folding settings. Uses the explicit configuration argument.
        fn fingerprint_topological_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalFingerprintParams,
        ) -> PyResult<Fingerprint> {
            let params = params.to_rust(py)?;
            self.inner
                .fingerprint_topological_with_params(&params)
                .map(|inner| Fingerprint { inner })
                .map_err(TopologicalFingerprintError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing, request=crate::arguments::Argument::Missing, *, atom_bits=crate::arguments::Argument::Missing, bit_info=crate::arguments::Argument::Missing, min_path=crate::arguments::Argument::Missing, max_path=crate::arguments::Argument::Missing, fp_size=crate::arguments::Argument::Missing, num_bits_per_feature=crate::arguments::Argument::Missing, use_hs=crate::arguments::Argument::Missing, target_density=crate::arguments::Argument::Missing, min_size=crate::arguments::Argument::Missing, branched_paths=crate::arguments::Argument::Missing, use_bond_order=crate::arguments::Argument::Missing, atom_invariants=crate::arguments::Argument::Missing, from_atoms=crate::arguments::Argument::Missing, ignore_atoms=crate::arguments::Argument::Missing))]
        /// Calculate the path-based topological fingerprint and requested atom/bit provenance.
        fn fingerprint_topological_with_output(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, TopologicalFingerprintParams>>>,
            request: crate::arguments::Argument<
                Option<PyRef<'_, TopologicalFingerprintOutputRequest>>,
            >,
            atom_bits: crate::arguments::Argument<Option<bool>>,
            bit_info: crate::arguments::Argument<Option<bool>>,
            min_path: crate::arguments::Argument<Option<u32>>,
            max_path: crate::arguments::Argument<Option<u32>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            num_bits_per_feature: crate::arguments::Argument<Option<u32>>,
            use_hs: crate::arguments::Argument<Option<bool>>,
            target_density: crate::arguments::Argument<Option<f64>>,
            min_size: crate::arguments::Argument<Option<u32>>,
            branched_paths: crate::arguments::Argument<Option<bool>>,
            use_bond_order: crate::arguments::Argument<Option<bool>>,
            atom_invariants: crate::arguments::Argument<Option<Vec<u32>>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
            ignore_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<TopologicalFingerprintResult> {
            let request_keyword_form = atom_bits.is_supplied() || bit_info.is_supplied();
            crate::arguments::configuration_arguments!(request; atom_bits, bit_info);
            let request = request.as_deref();
            let keyword_form = min_path.is_supplied()
                || max_path.is_supplied()
                || fp_size.is_supplied()
                || num_bits_per_feature.is_supplied()
                || use_hs.is_supplied()
                || target_density.is_supplied()
                || min_size.is_supplied()
                || branched_paths.is_supplied()
                || use_bond_order.is_supplied()
                || atom_invariants.is_supplied()
                || from_atoms.is_supplied()
                || ignore_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; min_path, max_path, fp_size, num_bits_per_feature, use_hs, target_density, min_size, branched_paths, use_bond_order, atom_invariants, from_atoms, ignore_atoms);
            let params = params.as_deref();
            if request_keyword_form {
                return if keyword_form {
                    self._fingerprint_topological_with_output_keywords_request_keywords(
                        py,
                        atom_bits,
                        bit_info,
                        min_path,
                        max_path,
                        fp_size,
                        num_bits_per_feature,
                        use_hs,
                        target_density,
                        min_size,
                        branched_paths,
                        use_bond_order,
                        atom_invariants,
                        from_atoms,
                        ignore_atoms,
                    )
                } else {
                    self._fingerprint_topological_with_output_params_request_keywords(
                        py, params, atom_bits, bit_info,
                    )
                };
            }
            if keyword_form {
                self._fingerprint_topological_with_output_keywords(
                    py,
                    request,
                    min_path,
                    max_path,
                    fp_size,
                    num_bits_per_feature,
                    use_hs,
                    target_density,
                    min_size,
                    branched_paths,
                    use_bond_order,
                    atom_invariants,
                    from_atoms,
                    ignore_atoms,
                )
            } else {
                self._fingerprint_topological_with_output_params(py, params, request)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_with_output"}"#]
        #[pyo3(signature=(params=None, request=None))]
        fn _fingerprint_topological_with_output_params(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalFingerprintParams>,
            request: Option<&TopologicalFingerprintOutputRequest>,
        ) -> PyResult<TopologicalFingerprintResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            let request = request.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_with_output_params_native(py, params.as_ref(), request)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_with_output"}"#]
        #[pyo3(signature=(*, request=None, min_path=None, max_path=None, fp_size=None, num_bits_per_feature=None, use_hs=None, target_density=None, min_size=None, branched_paths=None, use_bond_order=None, atom_invariants=None, from_atoms=None, ignore_atoms=None))]
        fn _fingerprint_topological_with_output_keywords(
            &self,
            py: Python<'_>,
            request: Option<&TopologicalFingerprintOutputRequest>,
            min_path: Option<u32>,
            max_path: Option<u32>,
            fp_size: Option<u32>,
            num_bits_per_feature: Option<u32>,
            use_hs: Option<bool>,
            target_density: Option<f64>,
            min_size: Option<u32>,
            branched_paths: Option<bool>,
            use_bond_order: Option<bool>,
            atom_invariants: Option<Vec<u32>>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
        ) -> PyResult<TopologicalFingerprintResult> {
            let params = super::topological_producers::keyword_params(
                min_path,
                max_path,
                fp_size,
                num_bits_per_feature,
                use_hs,
                target_density,
                min_size,
                branched_paths,
                use_bond_order,
                atom_invariants,
                from_atoms,
                ignore_atoms,
            );
            let request = request.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_topological_with_output_params_native(py, Some(&params), request)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_with_output"}"#]
        #[pyo3(signature=(params=None, *, atom_bits=None, bit_info=None))]
        fn _fingerprint_topological_with_output_params_request_keywords(
            &self,
            py: Python<'_>,
            params: Option<&TopologicalFingerprintParams>,
            atom_bits: Option<bool>,
            bit_info: Option<bool>,
        ) -> PyResult<TopologicalFingerprintResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            let request = super::topological_producers::request_keywords(atom_bits, bit_info);
            self._fingerprint_topological_with_output_params_native(
                py,
                params.as_ref(),
                Some(request),
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_topological_with_output"}"#]
        #[pyo3(signature=(*, atom_bits=None, bit_info=None, min_path=None, max_path=None, fp_size=None, num_bits_per_feature=None, use_hs=None, target_density=None, min_size=None, branched_paths=None, use_bond_order=None, atom_invariants=None, from_atoms=None, ignore_atoms=None))]
        fn _fingerprint_topological_with_output_keywords_request_keywords(
            &self,
            py: Python<'_>,
            atom_bits: Option<bool>,
            bit_info: Option<bool>,
            min_path: Option<u32>,
            max_path: Option<u32>,
            fp_size: Option<u32>,
            num_bits_per_feature: Option<u32>,
            use_hs: Option<bool>,
            target_density: Option<f64>,
            min_size: Option<u32>,
            branched_paths: Option<bool>,
            use_bond_order: Option<bool>,
            atom_invariants: Option<Vec<u32>>,
            from_atoms: Option<Vec<u32>>,
            ignore_atoms: Option<Vec<u32>>,
        ) -> PyResult<TopologicalFingerprintResult> {
            let params = super::topological_producers::keyword_params(
                min_path,
                max_path,
                fp_size,
                num_bits_per_feature,
                use_hs,
                target_density,
                min_size,
                branched_paths,
                use_bond_order,
                atom_invariants,
                from_atoms,
                ignore_atoms,
            );
            let request = super::topological_producers::request_keywords(atom_bits, bit_info);
            self._fingerprint_topological_with_output_params_native(
                py,
                Some(&params),
                Some(request),
            )
        }
        /// Calculate the path-based topological fingerprint and requested atom/bit provenance. Uses the explicit configuration argument.
        fn fingerprint_topological_with_output_with_params(
            &self,
            py: Python<'_>,
            params: &TopologicalFingerprintParams,
            request: &TopologicalFingerprintOutputRequest,
        ) -> PyResult<TopologicalFingerprintResult> {
            let request = pyo3::Python::attach(|py| request.to_rust(py))?;

            let params = params.to_rust(py)?;
            self.inner
                .fingerprint_topological_with_output_with_params(&params, request)
                .map(TopologicalFingerprintResult::from_rust)
                .map_err(TopologicalFingerprintError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, layers=crate::arguments::Argument::Missing, min_path=crate::arguments::Argument::Missing, max_path=crate::arguments::Argument::Missing, fp_size=crate::arguments::Argument::Missing, atom_counts=crate::arguments::Argument::Missing, set_only_bits=crate::arguments::Argument::Missing, branched_paths=crate::arguments::Argument::Missing, from_atoms=crate::arguments::Argument::Missing))]
        /// Calculate a layered path fingerprint using the selected layer flags, path bounds and bit width.
        fn fingerprint_layered(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, LayeredFingerprintParams>>>,
            layers: crate::arguments::Argument<Option<u32>>,
            min_path: crate::arguments::Argument<Option<u32>>,
            max_path: crate::arguments::Argument<Option<u32>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            atom_counts: crate::arguments::Argument<Option<Vec<u32>>>,
            set_only_bits: crate::arguments::Argument<Option<PyRef<'_, Fingerprint>>>,
            branched_paths: crate::arguments::Argument<Option<bool>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = layers.is_supplied()
                || min_path.is_supplied()
                || max_path.is_supplied()
                || fp_size.is_supplied()
                || atom_counts.is_supplied()
                || set_only_bits.is_supplied()
                || branched_paths.is_supplied()
                || from_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; layers, min_path, max_path, fp_size, atom_counts, set_only_bits, branched_paths, from_atoms);
            let params = params.as_deref();
            let set_only_bits = set_only_bits.as_deref();
            if keyword_form {
                self._fingerprint_layered_keywords(
                    py,
                    layers,
                    min_path,
                    max_path,
                    fp_size,
                    atom_counts,
                    set_only_bits,
                    branched_paths,
                    from_atoms,
                )
            } else {
                self._fingerprint_layered_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_layered"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_layered_params(
            &self,
            py: Python<'_>,
            params: Option<&LayeredFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_layered_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_layered"}"#]
        #[pyo3(signature=(*, layers=None, min_path=None, max_path=None, fp_size=None, atom_counts=None, set_only_bits=None, branched_paths=None, from_atoms=None))]
        fn _fingerprint_layered_keywords(
            &self,
            py: Python<'_>,
            layers: Option<u32>,
            min_path: Option<u32>,
            max_path: Option<u32>,
            fp_size: Option<u32>,
            atom_counts: Option<Vec<u32>>,
            set_only_bits: Option<&Fingerprint>,
            branched_paths: Option<bool>,
            from_atoms: Option<Vec<u32>>,
        ) -> PyResult<Fingerprint> {
            let params = super::layered_producers::keyword_params(
                layers,
                min_path,
                max_path,
                fp_size,
                atom_counts,
                set_only_bits,
                branched_paths,
                from_atoms,
            );
            self._fingerprint_layered_params_native(py, Some(&params))
        }
        /// Calculate a layered path fingerprint using the selected layer flags, path bounds and bit width. Uses the explicit configuration argument.
        fn fingerprint_layered_with_params(
            &self,
            py: Python<'_>,
            params: &LayeredFingerprintParams,
        ) -> PyResult<Fingerprint> {
            let params = params.to_rust(py)?;
            self.inner
                .fingerprint_layered_with_params(&params)
                .map(|inner| Fingerprint { inner })
                .map_err(LayeredFingerprintError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing, *, layers=crate::arguments::Argument::Missing, min_path=crate::arguments::Argument::Missing, max_path=crate::arguments::Argument::Missing, fp_size=crate::arguments::Argument::Missing, atom_counts=crate::arguments::Argument::Missing, set_only_bits=crate::arguments::Argument::Missing, branched_paths=crate::arguments::Argument::Missing, from_atoms=crate::arguments::Argument::Missing))]
        /// Calculate a layered path fingerprint together with its native additional output.
        fn fingerprint_layered_with_output(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, LayeredFingerprintParams>>>,
            layers: crate::arguments::Argument<Option<u32>>,
            min_path: crate::arguments::Argument<Option<u32>>,
            max_path: crate::arguments::Argument<Option<u32>>,
            fp_size: crate::arguments::Argument<Option<u32>>,
            atom_counts: crate::arguments::Argument<Option<Vec<u32>>>,
            set_only_bits: crate::arguments::Argument<Option<PyRef<'_, Fingerprint>>>,
            branched_paths: crate::arguments::Argument<Option<bool>>,
            from_atoms: crate::arguments::Argument<Option<Vec<u32>>>,
        ) -> PyResult<LayeredFingerprintResult> {
            let keyword_form = layers.is_supplied()
                || min_path.is_supplied()
                || max_path.is_supplied()
                || fp_size.is_supplied()
                || atom_counts.is_supplied()
                || set_only_bits.is_supplied()
                || branched_paths.is_supplied()
                || from_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; layers, min_path, max_path, fp_size, atom_counts, set_only_bits, branched_paths, from_atoms);
            let params = params.as_deref();
            let set_only_bits = set_only_bits.as_deref();
            if keyword_form {
                self._fingerprint_layered_with_output_keywords(
                    py,
                    layers,
                    min_path,
                    max_path,
                    fp_size,
                    atom_counts,
                    set_only_bits,
                    branched_paths,
                    from_atoms,
                )
            } else {
                self._fingerprint_layered_with_output_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_layered_with_output"}"#]
        #[pyo3(signature=(params=None))]
        fn _fingerprint_layered_with_output_params(
            &self,
            py: Python<'_>,
            params: Option<&LayeredFingerprintParams>,
        ) -> PyResult<LayeredFingerprintResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._fingerprint_layered_with_output_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_layered_with_output"}"#]
        #[pyo3(signature=(*, layers=None, min_path=None, max_path=None, fp_size=None, atom_counts=None, set_only_bits=None, branched_paths=None, from_atoms=None))]
        fn _fingerprint_layered_with_output_keywords(
            &self,
            py: Python<'_>,
            layers: Option<u32>,
            min_path: Option<u32>,
            max_path: Option<u32>,
            fp_size: Option<u32>,
            atom_counts: Option<Vec<u32>>,
            set_only_bits: Option<&Fingerprint>,
            branched_paths: Option<bool>,
            from_atoms: Option<Vec<u32>>,
        ) -> PyResult<LayeredFingerprintResult> {
            let params = super::layered_producers::keyword_params(
                layers,
                min_path,
                max_path,
                fp_size,
                atom_counts,
                set_only_bits,
                branched_paths,
                from_atoms,
            );
            self._fingerprint_layered_with_output_params_native(py, Some(&params))
        }
        /// Calculate a layered path fingerprint together with its native additional output. Uses the explicit configuration argument.
        fn fingerprint_layered_with_output_with_params(
            &self,
            py: Python<'_>,
            params: &LayeredFingerprintParams,
        ) -> PyResult<LayeredFingerprintResult> {
            let params = params.to_rust(py)?;
            self.inner
                .fingerprint_layered_with_output_with_params(&params)
                .map(|inner| LayeredFingerprintResult { inner })
                .map_err(LayeredFingerprintError::from_rust)
        }
        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, n_bits=crate::arguments::Argument::Missing, tautomeric=crate::arguments::Argument::Missing))]
        /// Calculate the canonical pattern fingerprint using the requested bit width and pattern options.
        fn fingerprint_pattern(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, PatternFingerprintParams>>>,
            n_bits: crate::arguments::Argument<Option<usize>>,
            tautomeric: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = n_bits.is_supplied() || tautomeric.is_supplied();
            crate::arguments::configuration_arguments!(params; n_bits, tautomeric);
            let params = params.as_deref();
            if keyword_form {
                self._fingerprint_pattern_keywords(n_bits, tautomeric)
            } else {
                self._fingerprint_pattern_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_pattern"}"#]
        #[pyo3(signature = (params=None))]
        fn _fingerprint_pattern_params(
            &self,
            params: Option<&PatternFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(params) => self.inner.fingerprint_pattern_with_params(&params),
                None => self.inner.fingerprint_pattern(),
            }
            .map(|inner| Fingerprint { inner })
            .map_err(PatternFingerprintError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_pattern"}"#]
        #[pyo3(signature = (*, n_bits=None, tautomeric=None))]
        fn _fingerprint_pattern_keywords(
            &self,
            n_bits: Option<usize>,
            tautomeric: Option<bool>,
        ) -> PyResult<Fingerprint> {
            let params = super::pattern_producers::keyword_params(n_bits, tautomeric);
            self._fingerprint_pattern_params(Some(&params))
        }
        /// Calculate the canonical pattern fingerprint using the requested bit width and pattern options. Uses the explicit configuration argument.
        fn fingerprint_pattern_with_params(
            &self,
            params: &PatternFingerprintParams,
        ) -> PyResult<Fingerprint> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            self.inner
                .fingerprint_pattern_with_params(&params)
                .map(|inner| Fingerprint { inner })
                .map_err(PatternFingerprintError::from_rust)
        }
        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, operations=crate::arguments::Argument::Missing))]
        /// Inspect the molecule under the selected sanitization checks and return structured chemical-problem diagnostics.
        fn detect_chemistry_problems(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, SanitizeParams>>>,
            operations: crate::arguments::Argument<
                Option<super::sanitize_values::SanitizeOperationsInput>,
            >,
        ) -> PyResult<ChemistryProblemReport> {
            let keyword_form = operations.is_supplied();
            crate::arguments::configuration_arguments!(params; operations);
            let params = params.as_deref();
            if keyword_form {
                self._detect_chemistry_problems_keywords(operations)
            } else {
                self._detect_chemistry_problems_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.detect_chemistry_problems"}"#]
        #[pyo3(signature = (params=None))]
        fn _detect_chemistry_problems_params(
            &self,
            params: Option<&SanitizeParams>,
        ) -> PyResult<ChemistryProblemReport> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(params) => self.inner.detect_chemistry_problems_with_params(&params),
                None => self.inner.detect_chemistry_problems(),
            }
            .map(|inner| ChemistryProblemReport { inner })
            .map_err(SanitizeError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.detect_chemistry_problems"}"#]
        #[pyo3(signature = (*, operations=None))]
        fn _detect_chemistry_problems_keywords(
            &self,
            operations: Option<super::sanitize_values::SanitizeOperationsInput>,
        ) -> PyResult<ChemistryProblemReport> {
            let params = SanitizeParams {
                inner: crate::parameters::NativeField::Owned(ck::SanitizeParams {
                    operations: operations
                        .map(super::sanitize_values::SanitizeOperationsInput::into_rust)
                        .transpose()?
                        .unwrap_or_else(|| ck::SanitizeParams::default().operations),
                }),
            };
            self._detect_chemistry_problems_params(Some(&params))
        }
        /// Inspect the molecule under the selected sanitization checks and return structured chemical-problem diagnostics. This entry accepts the explicit native parameter form.
        fn detect_chemistry_problems_with_params(
            &self,
            params: &SanitizeParams,
        ) -> PyResult<ChemistryProblemReport> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            self.inner
                .detect_chemistry_problems_with_params(&params)
                .map(|inner| ChemistryProblemReport { inner })
                .map_err(SanitizeError::from_rust)
        }
        #[pyo3(signature = (path, size, atom_codes=None))]
        /// Calculate the canonical topological torsion score for the supplied atom path.
        fn topological_torsion_path_score(
            &self,
            path: Vec<usize>,
            size: usize,
            atom_codes: Option<Vec<u32>>,
        ) -> PyResult<u64> {
            self.inner
                .topological_torsion_path_score(&path, size, atom_codes.as_deref())
                .map_err(TopologicalTorsionPathScoreError::from_rust)
        }
        /// Count atoms other than hydrogen.
        fn num_heavy_atoms(&self) -> PyResult<u32> {
            self.inner
                .num_heavy_atoms()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count atoms including explicit and implicit hydrogens, rather than only stored atom rows.
        fn total_atom_count(&self) -> PyResult<u32> {
            self.inner
                .total_atom_count()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return the Lipinski nitrogen/oxygen acceptor count, not the SMARTS-based num_hba definition.
        fn lipinski_hba(&self) -> PyResult<u32> {
            self.inner
                .lipinski_hba()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return the Lipinski donor-hydrogen count, not the SMARTS-based donor-atom count num_hbd.
        fn lipinski_hbd(&self) -> PyResult<u32> {
            self.inner
                .lipinski_hbd()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count atoms other than carbon and hydrogen using the canonical heteroatom query.
        fn num_heteroatoms(&self) -> PyResult<u32> {
            self.inner
                .num_heteroatoms()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count hydrogen-bond acceptor atoms using the canonical SMARTS definitions. Requires prepared valence and ring information.
        fn num_hba(&self) -> PyResult<u32> {
            self.inner.num_hba().map_err(DescriptorReadError::from_rust)
        }
        /// Count hydrogen-bond donor atoms using the canonical SMARTS definitions. Requires prepared valence information.
        fn num_hbd(&self) -> PyResult<u32> {
            self.inner.num_hbd().map_err(DescriptorReadError::from_rust)
        }
        /// Count amide bonds using the canonical amide pattern.
        fn num_amide_bonds(&self) -> PyResult<u32> {
            self.inner
                .num_amide_bonds()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count assigned and unassigned atom stereocenters recognized by the canonical stereo analysis.
        fn num_atom_stereo_centers(&self) -> PyResult<u32> {
            self.inner
                .num_atom_stereo_centers()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count recognized atom stereocenters without a specified stereochemical assignment.
        fn num_unspecified_atom_stereo_centers(&self) -> PyResult<u32> {
            self.inner
                .num_unspecified_atom_stereo_centers()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count spiro atoms from the molecule's ring information.
        fn num_spiro_atoms(&self) -> PyResult<u32> {
            self.inner
                .num_spiro_atoms()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count bridgehead atoms from the molecule's ring information.
        fn num_bridgehead_atoms(&self) -> PyResult<u32> {
            self.inner
                .num_bridgehead_atoms()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return the fraction of carbon atoms with sp3 hybridization.
        fn fraction_csp3(&self) -> PyResult<f64> {
            self.inner
                .fraction_csp3()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the quantitative estimate of drug-likeness with the canonical default weights.
        fn qed(&self) -> PyResult<f64> {
            self.inner.qed().map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 1 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_1(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_1()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 2 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_2(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_2()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 3 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_3(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_3()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 4 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_4(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_4()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 5 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_5(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_5()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 6 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_6(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_6()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 7 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_7(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_7()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 8 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_8(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_8()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 9 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_9(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_9()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 10 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_10(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_10()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 11 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_11(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_11()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 12 of the default surface-area descriptor classified by atomic Crippen logP.
        fn slogp_vsa_12(&self) -> PyResult<f64> {
            self.inner
                .slogp_vsa_12()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 1 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_1(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_1()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 2 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_2(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_2()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 3 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_3(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_3()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 4 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_4(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_4()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 5 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_5(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_5()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 6 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_6(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_6()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 7 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_7(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_7()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 8 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_8(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_8()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 9 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_9(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_9()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return bin 10 of the default surface-area descriptor classified by atomic Crippen molar refractivity.
        fn smr_vsa_10(&self) -> PyResult<f64> {
            self.inner
                .smr_vsa_10()
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (params=None))]
        /// Count rotatable bonds under the selected strictness definition; None uses the native default.
        fn num_rotatable_bonds(&self, params: Option<RotatableBondsOptions>) -> PyResult<u32> {
            match params {
                Some(params) => self.inner.num_rotatable_bonds_with_params(&params.into()),
                None => self.inner.num_rotatable_bonds(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Count rotatable bonds under the explicitly selected strictness definition.
        fn num_rotatable_bonds_with_params(&self, params: RotatableBondsOptions) -> PyResult<u32> {
            self.inner
                .num_rotatable_bonds_with_params(&params.into())
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, only_heavy=None))]
        /// Calculate molecular weight using average atomic masses. only_heavy excludes hydrogens.
        fn molecular_weight(&self, only_heavy: Option<bool>) -> PyResult<f64> {
            match only_heavy {
                Some(only_heavy) => self.inner.molecular_weight_with_params(only_heavy),
                None => self.inner.molecular_weight(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate molecular weight using average atomic masses. only_heavy excludes hydrogens. This entry accepts the explicit native parameter form.
        fn molecular_weight_with_params(&self, only_heavy: bool) -> PyResult<f64> {
            self.inner
                .molecular_weight_with_params(only_heavy)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, only_heavy=None))]
        /// Calculate exact molecular weight using isotopic masses. only_heavy excludes hydrogens.
        fn exact_molecular_weight(&self, only_heavy: Option<bool>) -> PyResult<f64> {
            match only_heavy {
                Some(only_heavy) => self.inner.exact_molecular_weight_with_params(only_heavy),
                None => self.inner.exact_molecular_weight(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate exact molecular weight using isotopic masses. only_heavy excludes hydrogens. This entry accepts the explicit native parameter form.
        fn exact_molecular_weight_with_params(&self, only_heavy: bool) -> PyResult<f64> {
            self.inner
                .exact_molecular_weight_with_params(only_heavy)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, separate_isotopes=None, abbreviate_h_isotopes=None))]
        /// Return the molecular formula, with optional isotope separation and hydrogen-isotope abbreviations.
        fn molecular_formula(
            &self,
            separate_isotopes: Option<bool>,
            abbreviate_h_isotopes: Option<bool>,
        ) -> PyResult<String> {
            match (separate_isotopes, abbreviate_h_isotopes) {
                (None, None) => self.inner.molecular_formula(),
                (separate, abbreviate) => self.inner.molecular_formula_with_params(
                    separate.unwrap_or(false),
                    abbreviate.unwrap_or(false),
                ),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Return the molecular formula, with optional isotope separation and hydrogen-isotope abbreviations. This entry accepts the explicit native parameter form.
        fn molecular_formula_with_params(
            &self,
            separate_isotopes: bool,
            abbreviate_h_isotopes: bool,
        ) -> PyResult<String> {
            self.inner
                .molecular_formula_with_params(separate_isotopes, abbreviate_h_isotopes)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, include_hydrogens=None, force=None))]
        /// Calculate the Wildman-Crippen logP and molar-refractivity descriptors. include_hydrogens controls hydrogen contributions; force requests recalculation.
        fn crippen_descriptors(
            &self,
            include_hydrogens: Option<bool>,
            force: Option<bool>,
        ) -> PyResult<CrippenTotals> {
            match (include_hydrogens, force) {
                (None, None) => self.inner.crippen_descriptors(),
                (option, force) => self.inner.crippen_descriptors_with_params(
                    option.unwrap_or(true),
                    force.unwrap_or(false),
                ),
            }
            .map(|inner| CrippenTotals { inner })
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the Wildman-Crippen logP and molar-refractivity descriptors. include_hydrogens controls hydrogen contributions; force requests recalculation. This entry accepts the explicit native parameter form.
        fn crippen_descriptors_with_params(
            &self,
            include_hydrogens: bool,
            force: bool,
        ) -> PyResult<CrippenTotals> {
            self.inner
                .crippen_descriptors_with_params(include_hydrogens, force)
                .map(|inner| CrippenTotals { inner })
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, include_hydrogens=None, force=None))]
        /// Calculate Labute's approximate molecular surface area. include_hydrogens controls hydrogen contributions; force requests recalculation.
        fn labute_asa(
            &self,
            include_hydrogens: Option<bool>,
            force: Option<bool>,
        ) -> PyResult<f64> {
            match (include_hydrogens, force) {
                (None, None) => self.inner.labute_asa(),
                (option, force) => self
                    .inner
                    .labute_asa_with_params(option.unwrap_or(true), force.unwrap_or(false)),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate Labute's approximate molecular surface area. include_hydrogens controls hydrogen contributions; force requests recalculation. This entry accepts the explicit native parameter form.
        fn labute_asa_with_params(&self, include_hydrogens: bool, force: bool) -> PyResult<f64> {
            self.inner
                .labute_asa_with_params(include_hydrogens, force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, include_hydrogens=None, force=None))]
        /// Return per-atom Labute surface contributions and the hydrogen contribution, preserving atom-index order.
        fn labute_asa_contributions(
            &self,
            include_hydrogens: Option<bool>,
            force: Option<bool>,
        ) -> PyResult<LabuteAsaContributions> {
            match (include_hydrogens, force) {
                (None, None) => self.inner.labute_asa_contributions(),
                (option, force) => self.inner.labute_asa_contributions_with_params(
                    option.unwrap_or(true),
                    force.unwrap_or(false),
                ),
            }
            .map(|inner| LabuteAsaContributions { inner })
            .map_err(DescriptorReadError::from_rust)
        }
        /// Return per-atom Labute surface contributions and the hydrogen contribution, preserving atom-index order. This entry accepts the explicit native parameter form.
        fn labute_asa_contributions_with_params(
            &self,
            include_hydrogens: bool,
            force: bool,
        ) -> PyResult<LabuteAsaContributions> {
            self.inner
                .labute_asa_contributions_with_params(include_hydrogens, force)
                .map(|inner| LabuteAsaContributions { inner })
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, include_sulfur_phosphorus=None, force=None))]
        /// Calculate topological polar surface area. include_sulfur_phosphorus includes sulfur and phosphorus fragments; force requests recalculation.
        fn tpsa(
            &self,
            include_sulfur_phosphorus: Option<bool>,
            force: Option<bool>,
        ) -> PyResult<f64> {
            match (include_sulfur_phosphorus, force) {
                (None, None) => self.inner.tpsa(),
                (option, force) => self
                    .inner
                    .tpsa_with_params(option.unwrap_or(false), force.unwrap_or(false)),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate topological polar surface area. include_sulfur_phosphorus includes sulfur and phosphorus fragments; force requests recalculation. This entry accepts the explicit native parameter form.
        fn tpsa_with_params(&self, include_sulfur_phosphorus: bool, force: bool) -> PyResult<f64> {
            self.inner
                .tpsa_with_params(include_sulfur_phosphorus, force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, bins=None, force=None))]
        /// Return surface-area bins classified by atomic Crippen logP contributions. bins supplies ordered cut points; None uses native default bins.
        fn slogp_vsa(&self, bins: Option<Vec<f64>>, force: Option<bool>) -> PyResult<Vec<f64>> {
            match (bins, force) {
                (None, None) => self.inner.slogp_vsa(),
                (bins, force) => self
                    .inner
                    .slogp_vsa_with_params(bins.as_deref(), force.unwrap_or(false)),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Return surface-area bins classified by atomic Crippen logP contributions. bins supplies ordered cut points; None uses native default bins. This entry accepts the explicit native parameter form.
        fn slogp_vsa_with_params(&self, bins: Option<Vec<f64>>, force: bool) -> PyResult<Vec<f64>> {
            self.inner
                .slogp_vsa_with_params(bins.as_deref(), force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, bins=None, force=None))]
        /// Return surface-area bins classified by atomic Crippen molar-refractivity contributions. bins supplies ordered cut points; None uses native defaults.
        fn smr_vsa(&self, bins: Option<Vec<f64>>, force: Option<bool>) -> PyResult<Vec<f64>> {
            match (bins, force) {
                (None, None) => self.inner.smr_vsa(),
                (bins, force) => self
                    .inner
                    .smr_vsa_with_params(bins.as_deref(), force.unwrap_or(false)),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Return surface-area bins classified by atomic Crippen molar-refractivity contributions. bins supplies ordered cut points; None uses native defaults. This entry accepts the explicit native parameter form.
        fn smr_vsa_with_params(&self, bins: Option<Vec<f64>>, force: bool) -> PyResult<Vec<f64>> {
            self.inner
                .smr_vsa_with_params(bins.as_deref(), force)
                .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the Hall-Kier alpha correction used by molecular shape descriptors.
        fn hall_kier_alpha(&self) -> PyResult<f64> {
            self.inner
                .hall_kier_alpha()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return Hall-Kier alpha and its per-atom contributions in atom-index order.
        fn hall_kier_alpha_with_contributions(&self) -> PyResult<(f64, Vec<f64>)> {
            self.inner
                .hall_kier_alpha_with_contributions()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the first-order Kier kappa molecular shape index.
        fn kappa_1(&self) -> PyResult<f64> {
            self.inner.kappa_1().map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the second-order Kier kappa molecular shape index.
        fn kappa_2(&self) -> PyResult<f64> {
            self.inner.kappa_2().map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the third-order Kier kappa molecular shape index.
        fn kappa_3(&self) -> PyResult<f64> {
            self.inner.kappa_3().map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the molecular flexibility index from the canonical kappa descriptors.
        fn phi(&self) -> PyResult<f64> {
            self.inner.phi().map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (force=false))]
        /// Return the 42 molecular quantum numbers in their canonical feature order. force requests recalculation.
        fn mqns(&self, force: bool) -> PyResult<Vec<u32>> {
            self.inner
                .mqns(force)
                .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the zeroth-order connectivity index from atomic degrees.
        fn chi_0(&self) -> PyResult<f64> {
            self.inner.chi_0().map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the first-order connectivity index from bonded atomic degrees.
        fn chi_1(&self) -> PyResult<f64> {
            self.inner.chi_1().map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-0 valence connectivity index. force requests recalculation.
        fn chi_0_v(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_0_v_with_params(force),
                None => self.inner.chi_0_v(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-0 valence connectivity index. force requests recalculation.
        fn chi_0_v_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_0_v_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-1 valence connectivity index. force requests recalculation.
        fn chi_1_v(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_1_v_with_params(force),
                None => self.inner.chi_1_v(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-1 valence connectivity index. force requests recalculation.
        fn chi_1_v_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_1_v_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-2 valence connectivity index. force requests recalculation.
        fn chi_2_v(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_2_v_with_params(force),
                None => self.inner.chi_2_v(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-2 valence connectivity index. force requests recalculation.
        fn chi_2_v_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_2_v_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-3 valence connectivity index. force requests recalculation.
        fn chi_3_v(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_3_v_with_params(force),
                None => self.inner.chi_3_v(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-3 valence connectivity index. force requests recalculation.
        fn chi_3_v_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_3_v_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-4 valence connectivity index. force requests recalculation.
        fn chi_4_v(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_4_v_with_params(force),
                None => self.inner.chi_4_v(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-4 valence connectivity index. force requests recalculation.
        fn chi_4_v_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_4_v_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-0 normalized connectivity index. force requests recalculation.
        fn chi_0_n(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_0_n_with_params(force),
                None => self.inner.chi_0_n(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-0 normalized connectivity index. force requests recalculation.
        fn chi_0_n_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_0_n_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-1 normalized connectivity index. force requests recalculation.
        fn chi_1_n(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_1_n_with_params(force),
                None => self.inner.chi_1_n(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-1 normalized connectivity index. force requests recalculation.
        fn chi_1_n_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_1_n_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-2 normalized connectivity index. force requests recalculation.
        fn chi_2_n(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_2_n_with_params(force),
                None => self.inner.chi_2_n(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-2 normalized connectivity index. force requests recalculation.
        fn chi_2_n_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_2_n_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-3 normalized connectivity index. force requests recalculation.
        fn chi_3_n(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_3_n_with_params(force),
                None => self.inner.chi_3_n(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-3 normalized connectivity index. force requests recalculation.
        fn chi_3_n_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_3_n_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (*, force=None))]
        /// Calculate the order-4 normalized connectivity index. force requests recalculation.
        fn chi_4_n(&self, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_4_n_with_params(force),
                None => self.inner.chi_4_n(),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the order-4 normalized connectivity index. force requests recalculation.
        fn chi_4_n_with_params(&self, force: bool) -> PyResult<f64> {
            self.inner
                .chi_4_n_with_params(force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (order, *, force=None))]
        /// Calculate the valence connectivity index at the requested path order. force requests recalculation.
        fn chi_n_v(&self, order: u32, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_n_v_with_params(order, force),
                None => self.inner.chi_n_v(order),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the valence connectivity index at the requested path order. force requests recalculation. Uses the explicitly supplied parameter form.
        fn chi_n_v_with_params(&self, order: u32, force: bool) -> PyResult<f64> {
            self.inner
                .chi_n_v_with_params(order, force)
                .map_err(DescriptorReadError::from_rust)
        }
        #[pyo3(signature = (order, *, force=None))]
        /// Calculate the normalized connectivity index at the requested path order. force requests recalculation.
        fn chi_n_n(&self, order: u32, force: Option<bool>) -> PyResult<f64> {
            match force {
                Some(force) => self.inner.chi_n_n_with_params(order, force),
                None => self.inner.chi_n_n(order),
            }
            .map_err(DescriptorReadError::from_rust)
        }
        /// Calculate the normalized connectivity index at the requested path order. force requests recalculation. Uses the explicitly supplied parameter form.
        fn chi_n_n_with_params(&self, order: u32, force: bool) -> PyResult<f64> {
            self.inner
                .chi_n_n_with_params(order, force)
                .map_err(DescriptorReadError::from_rust)
        }

        /// Read the initialized ring count. Does not calculate rings; missing initialized ring information raises DescriptorReadError.
        fn num_rings(&self) -> PyResult<u32> {
            self.inner
                .num_rings()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count rings containing a non-carbon atom in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_heterocycles(&self) -> PyResult<u32> {
            self.inner
                .num_heterocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count aromatic rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aromatic_rings(&self) -> PyResult<u32> {
            self.inner
                .num_aromatic_rings()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count saturated rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_saturated_rings(&self) -> PyResult<u32> {
            self.inner
                .num_saturated_rings()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count non-aromatic rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aliphatic_rings(&self) -> PyResult<u32> {
            self.inner
                .num_aliphatic_rings()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count aromatic rings containing a non-carbon atom in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aromatic_heterocycles(&self) -> PyResult<u32> {
            self.inner
                .num_aromatic_heterocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count aromatic all-carbon rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aromatic_carbocycles(&self) -> PyResult<u32> {
            self.inner
                .num_aromatic_carbocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count non-aromatic rings containing a non-carbon atom in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aliphatic_heterocycles(&self) -> PyResult<u32> {
            self.inner
                .num_aliphatic_heterocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count non-aromatic all-carbon rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_aliphatic_carbocycles(&self) -> PyResult<u32> {
            self.inner
                .num_aliphatic_carbocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count saturated rings containing a non-carbon atom in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_saturated_heterocycles(&self) -> PyResult<u32> {
            self.inner
                .num_saturated_heterocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Count saturated all-carbon rings in the stored ring information.
        /// Does not calculate rings; an empty or absent ring assignment yields zero.
        fn num_saturated_carbocycles(&self) -> PyResult<u32> {
            self.inner
                .num_saturated_carbocycles()
                .map_err(DescriptorReadError::from_rust)
        }
        /// Return (center, configuration) pairs. An undetermined configuration is None;
        /// within an ordered configuration, None denotes hydrogen or a lone pair.
        #[pyo3(signature = (include_unassigned=true))]
        fn tetrahedral_stereo(
            &self,
            include_unassigned: bool,
        ) -> PyResult<Vec<(usize, Option<Vec<Option<usize>>>)>> {
            self.inner
                .tetrahedral_stereo(include_unassigned)
                .map(|values| {
                    values
                        .into_iter()
                        .map(|value| {
                            (
                                value.center.index(),
                                value.ligands.map(|ligands| {
                                    ligands
                                        .into_iter()
                                        .map(|ligand| match ligand {
                                            ck::LigandRef::Atom(atom) => Some(atom.index()),
                                            ck::LigandRef::ImplicitHydrogen
                                            | ck::LigandRef::LonePair => None,
                                        })
                                        .collect()
                                }),
                            )
                        })
                        .collect()
                })
                .map_err(StereoReadError::from_rust)
        }
        /// Run the native stereochemistry perception query and propagate any error.
        fn perceive_stereochemistry(&self) -> PyResult<()> {
            self.inner
                .perceive_stereochemistry()
                .map_err(StereoReadError::from_rust)
        }
        #[pyo3(signature = (include_unassigned=false))]
        /// Return pairs of atom indices and chiral labels; include_unassigned controls inclusion of unassigned centers.
        fn find_chiral_centers(&self, include_unassigned: bool) -> PyResult<Vec<(usize, String)>> {
            self.inner
                .find_chiral_centers(include_unassigned)
                .map_err(StereoReadError::from_rust)
        }
        #[pyo3(signature = (query, params=crate::arguments::Argument::Missing, *, max_matches=crate::arguments::Argument::Missing, uniquify=crate::arguments::Argument::Missing, use_chirality=crate::arguments::Argument::Missing, use_enhanced_stereo=crate::arguments::Argument::Missing, specified_stereo_query_matches_unspecified=crate::arguments::Argument::Missing, use_query_query_matches=crate::arguments::Argument::Missing, recursion_possible=crate::arguments::Argument::Missing, max_recursive_matches=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing, aromatic_matches_conjugated=crate::arguments::Argument::Missing, aromatic_matches_single_or_double=crate::arguments::Argument::Missing, extra_atom_check_overrides_default_check=crate::arguments::Argument::Missing, extra_bond_check_overrides_default_check=crate::arguments::Argument::Missing, use_generic_matchers=crate::arguments::Argument::Missing, atom_properties=crate::arguments::Argument::Missing, bond_properties=crate::arguments::Argument::Missing, final_match=crate::arguments::Argument::Missing, atom_match=crate::arguments::Argument::Missing, bond_match=crate::arguments::Argument::Missing, extra_final_check=crate::arguments::Argument::Missing, extra_atom_check=crate::arguments::Argument::Missing, extra_bond_check=crate::arguments::Argument::Missing))]
        /// Find the first query match, returning its atom mapping or None. Supply SubstructMatchParams or keyword options, not both.
        fn substruct_match(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: crate::arguments::Argument<Option<PyRef<'_, SubstructMatchParams>>>,
            max_matches: crate::arguments::Argument<Option<usize>>,
            uniquify: crate::arguments::Argument<Option<bool>>,
            use_chirality: crate::arguments::Argument<Option<bool>>,
            use_enhanced_stereo: crate::arguments::Argument<Option<bool>>,
            specified_stereo_query_matches_unspecified: crate::arguments::Argument<Option<bool>>,
            use_query_query_matches: crate::arguments::Argument<Option<bool>>,
            recursion_possible: crate::arguments::Argument<Option<bool>>,
            max_recursive_matches: crate::arguments::Argument<Option<usize>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            aromatic_matches_conjugated: crate::arguments::Argument<Option<bool>>,
            aromatic_matches_single_or_double: crate::arguments::Argument<Option<bool>>,
            extra_atom_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            extra_bond_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            use_generic_matchers: crate::arguments::Argument<Option<bool>>,
            atom_properties: crate::arguments::Argument<Option<Vec<String>>>,
            bond_properties: crate::arguments::Argument<Option<Vec<String>>>,
            final_match: crate::arguments::Argument<Option<FinalMatch>>,
            atom_match: crate::arguments::Argument<Option<AtomMatch>>,
            bond_match: crate::arguments::Argument<Option<BondMatch>>,
            extra_final_check: crate::arguments::Argument<Option<FinalMatch>>,
            extra_atom_check: crate::arguments::Argument<Option<AtomMatch>>,
            extra_bond_check: crate::arguments::Argument<Option<BondMatch>>,
        ) -> PyResult<Option<MatchResult>> {
            let keyword_form = max_matches.is_supplied()
                || uniquify.is_supplied()
                || use_chirality.is_supplied()
                || use_enhanced_stereo.is_supplied()
                || specified_stereo_query_matches_unspecified.is_supplied()
                || use_query_query_matches.is_supplied()
                || recursion_possible.is_supplied()
                || max_recursive_matches.is_supplied()
                || num_threads.is_supplied()
                || aromatic_matches_conjugated.is_supplied()
                || aromatic_matches_single_or_double.is_supplied()
                || extra_atom_check_overrides_default_check.is_supplied()
                || extra_bond_check_overrides_default_check.is_supplied()
                || use_generic_matchers.is_supplied()
                || atom_properties.is_supplied()
                || bond_properties.is_supplied()
                || final_match.is_supplied()
                || atom_match.is_supplied()
                || bond_match.is_supplied()
                || extra_final_check.is_supplied()
                || extra_atom_check.is_supplied()
                || extra_bond_check.is_supplied();
            crate::arguments::configuration_arguments!(params; max_matches, uniquify, use_chirality, use_enhanced_stereo, specified_stereo_query_matches_unspecified, use_query_query_matches, recursion_possible, max_recursive_matches, num_threads, aromatic_matches_conjugated, aromatic_matches_single_or_double, extra_atom_check_overrides_default_check, extra_bond_check_overrides_default_check, use_generic_matchers, atom_properties, bond_properties, final_match, atom_match, bond_match, extra_final_check, extra_atom_check, extra_bond_check);
            let params = params.as_deref();
            if keyword_form {
                Self::_substruct_match_keywords(
                    slf,
                    query,
                    max_matches,
                    uniquify,
                    use_chirality,
                    use_enhanced_stereo,
                    specified_stereo_query_matches_unspecified,
                    use_query_query_matches,
                    recursion_possible,
                    max_recursive_matches,
                    num_threads,
                    aromatic_matches_conjugated,
                    aromatic_matches_single_or_double,
                    extra_atom_check_overrides_default_check,
                    extra_bond_check_overrides_default_check,
                    use_generic_matchers,
                    atom_properties,
                    bond_properties,
                    final_match,
                    atom_match,
                    bond_match,
                    extra_final_check,
                    extra_atom_check,
                    extra_bond_check,
                )
            } else {
                Self::_substruct_match_params(slf, query, params)
            }
        }

        #[doc = r#"{"role":"overload","target":"Molecule.substruct_match"}"#]
        #[pyo3(signature = (query, params=None))]
        fn _substruct_match_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<&SubstructMatchParams>,
        ) -> PyResult<Option<MatchResult>> {
            let py = slf.py();
            let params = params.map(|value| value.capture(py)).transpose()?;
            Self::_substruct_match_params_native(slf, query, params)
        }

        #[doc = r#"{"role":"overload","target":"Molecule.substruct_match"}"#]
        #[pyo3(signature = (query, *, max_matches=None, uniquify=None, use_chirality=None, use_enhanced_stereo=None, specified_stereo_query_matches_unspecified=None, use_query_query_matches=None, recursion_possible=None, max_recursive_matches=None, num_threads=None, aromatic_matches_conjugated=None, aromatic_matches_single_or_double=None, extra_atom_check_overrides_default_check=None, extra_bond_check_overrides_default_check=None, use_generic_matchers=None, atom_properties=None, bond_properties=None, final_match=None, atom_match=None, bond_match=None, extra_final_check=None, extra_atom_check=None, extra_bond_check=None))]
        fn _substruct_match_keywords(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            max_matches: Option<usize>,
            uniquify: Option<bool>,
            use_chirality: Option<bool>,
            use_enhanced_stereo: Option<bool>,
            specified_stereo_query_matches_unspecified: Option<bool>,
            use_query_query_matches: Option<bool>,
            recursion_possible: Option<bool>,
            max_recursive_matches: Option<usize>,
            num_threads: Option<i32>,
            aromatic_matches_conjugated: Option<bool>,
            aromatic_matches_single_or_double: Option<bool>,
            extra_atom_check_overrides_default_check: Option<bool>,
            extra_bond_check_overrides_default_check: Option<bool>,
            use_generic_matchers: Option<bool>,
            atom_properties: Option<Vec<String>>,
            bond_properties: Option<Vec<String>>,
            final_match: Option<FinalMatch>,
            atom_match: Option<AtomMatch>,
            bond_match: Option<BondMatch>,
            extra_final_check: Option<FinalMatch>,
            extra_atom_check: Option<AtomMatch>,
            extra_bond_check: Option<BondMatch>,
        ) -> PyResult<Option<MatchResult>> {
            let params = SubstructMatchParams::from_keywords(
                max_matches,
                uniquify,
                use_chirality,
                use_enhanced_stereo,
                specified_stereo_query_matches_unspecified,
                use_query_query_matches,
                recursion_possible,
                max_recursive_matches,
                num_threads,
                aromatic_matches_conjugated,
                aromatic_matches_single_or_double,
                extra_atom_check_overrides_default_check,
                extra_bond_check_overrides_default_check,
                use_generic_matchers,
                atom_properties,
                bond_properties,
                final_match,
                atom_match,
                bond_match,
                extra_final_check,
                extra_atom_check,
                extra_bond_check,
            )?;
            Self::_substruct_match_params_native(slf, query, Some(params))
        }

        /// Find the first query match, returning its atom mapping or None. Supply SubstructMatchParams or keyword options, not both. Uses the explicitly supplied parameter form.
        fn substruct_match_with_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: &SubstructMatchParams,
        ) -> PyResult<Option<MatchResult>> {
            Self::_substruct_match_params(slf, query, Some(params))
        }

        #[pyo3(signature = (query, params=crate::arguments::Argument::Missing, *, max_matches=crate::arguments::Argument::Missing, uniquify=crate::arguments::Argument::Missing, use_chirality=crate::arguments::Argument::Missing, use_enhanced_stereo=crate::arguments::Argument::Missing, specified_stereo_query_matches_unspecified=crate::arguments::Argument::Missing, use_query_query_matches=crate::arguments::Argument::Missing, recursion_possible=crate::arguments::Argument::Missing, max_recursive_matches=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing, aromatic_matches_conjugated=crate::arguments::Argument::Missing, aromatic_matches_single_or_double=crate::arguments::Argument::Missing, extra_atom_check_overrides_default_check=crate::arguments::Argument::Missing, extra_bond_check_overrides_default_check=crate::arguments::Argument::Missing, use_generic_matchers=crate::arguments::Argument::Missing, atom_properties=crate::arguments::Argument::Missing, bond_properties=crate::arguments::Argument::Missing, final_match=crate::arguments::Argument::Missing, atom_match=crate::arguments::Argument::Missing, bond_match=crate::arguments::Argument::Missing, extra_final_check=crate::arguments::Argument::Missing, extra_atom_check=crate::arguments::Argument::Missing, extra_bond_check=crate::arguments::Argument::Missing))]
        /// Return query matches in native order under the selected uniqueness, chirality and match-limit settings. Supply Params or keywords, not both.
        fn substruct_matches(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: crate::arguments::Argument<Option<PyRef<'_, SubstructMatchParams>>>,
            max_matches: crate::arguments::Argument<Option<usize>>,
            uniquify: crate::arguments::Argument<Option<bool>>,
            use_chirality: crate::arguments::Argument<Option<bool>>,
            use_enhanced_stereo: crate::arguments::Argument<Option<bool>>,
            specified_stereo_query_matches_unspecified: crate::arguments::Argument<Option<bool>>,
            use_query_query_matches: crate::arguments::Argument<Option<bool>>,
            recursion_possible: crate::arguments::Argument<Option<bool>>,
            max_recursive_matches: crate::arguments::Argument<Option<usize>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            aromatic_matches_conjugated: crate::arguments::Argument<Option<bool>>,
            aromatic_matches_single_or_double: crate::arguments::Argument<Option<bool>>,
            extra_atom_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            extra_bond_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            use_generic_matchers: crate::arguments::Argument<Option<bool>>,
            atom_properties: crate::arguments::Argument<Option<Vec<String>>>,
            bond_properties: crate::arguments::Argument<Option<Vec<String>>>,
            final_match: crate::arguments::Argument<Option<FinalMatch>>,
            atom_match: crate::arguments::Argument<Option<AtomMatch>>,
            bond_match: crate::arguments::Argument<Option<BondMatch>>,
            extra_final_check: crate::arguments::Argument<Option<FinalMatch>>,
            extra_atom_check: crate::arguments::Argument<Option<AtomMatch>>,
            extra_bond_check: crate::arguments::Argument<Option<BondMatch>>,
        ) -> PyResult<Vec<MatchResult>> {
            let keyword_form = max_matches.is_supplied()
                || uniquify.is_supplied()
                || use_chirality.is_supplied()
                || use_enhanced_stereo.is_supplied()
                || specified_stereo_query_matches_unspecified.is_supplied()
                || use_query_query_matches.is_supplied()
                || recursion_possible.is_supplied()
                || max_recursive_matches.is_supplied()
                || num_threads.is_supplied()
                || aromatic_matches_conjugated.is_supplied()
                || aromatic_matches_single_or_double.is_supplied()
                || extra_atom_check_overrides_default_check.is_supplied()
                || extra_bond_check_overrides_default_check.is_supplied()
                || use_generic_matchers.is_supplied()
                || atom_properties.is_supplied()
                || bond_properties.is_supplied()
                || final_match.is_supplied()
                || atom_match.is_supplied()
                || bond_match.is_supplied()
                || extra_final_check.is_supplied()
                || extra_atom_check.is_supplied()
                || extra_bond_check.is_supplied();
            crate::arguments::configuration_arguments!(params; max_matches, uniquify, use_chirality, use_enhanced_stereo, specified_stereo_query_matches_unspecified, use_query_query_matches, recursion_possible, max_recursive_matches, num_threads, aromatic_matches_conjugated, aromatic_matches_single_or_double, extra_atom_check_overrides_default_check, extra_bond_check_overrides_default_check, use_generic_matchers, atom_properties, bond_properties, final_match, atom_match, bond_match, extra_final_check, extra_atom_check, extra_bond_check);
            let params = params.as_deref();
            if keyword_form {
                Self::_substruct_matches_keywords(
                    slf,
                    query,
                    max_matches,
                    uniquify,
                    use_chirality,
                    use_enhanced_stereo,
                    specified_stereo_query_matches_unspecified,
                    use_query_query_matches,
                    recursion_possible,
                    max_recursive_matches,
                    num_threads,
                    aromatic_matches_conjugated,
                    aromatic_matches_single_or_double,
                    extra_atom_check_overrides_default_check,
                    extra_bond_check_overrides_default_check,
                    use_generic_matchers,
                    atom_properties,
                    bond_properties,
                    final_match,
                    atom_match,
                    bond_match,
                    extra_final_check,
                    extra_atom_check,
                    extra_bond_check,
                )
            } else {
                Self::_substruct_matches_params(slf, query, params)
            }
        }

        #[doc = r#"{"role":"overload","target":"Molecule.substruct_matches"}"#]
        #[pyo3(signature = (query, params=None))]
        fn _substruct_matches_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<&SubstructMatchParams>,
        ) -> PyResult<Vec<MatchResult>> {
            let py = slf.py();
            let params = params.map(|value| value.capture(py)).transpose()?;
            Self::_substruct_matches_params_native(slf, query, params)
        }

        #[doc = r#"{"role":"overload","target":"Molecule.substruct_matches"}"#]
        #[pyo3(signature = (query, *, max_matches=None, uniquify=None, use_chirality=None, use_enhanced_stereo=None, specified_stereo_query_matches_unspecified=None, use_query_query_matches=None, recursion_possible=None, max_recursive_matches=None, num_threads=None, aromatic_matches_conjugated=None, aromatic_matches_single_or_double=None, extra_atom_check_overrides_default_check=None, extra_bond_check_overrides_default_check=None, use_generic_matchers=None, atom_properties=None, bond_properties=None, final_match=None, atom_match=None, bond_match=None, extra_final_check=None, extra_atom_check=None, extra_bond_check=None))]
        fn _substruct_matches_keywords(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            max_matches: Option<usize>,
            uniquify: Option<bool>,
            use_chirality: Option<bool>,
            use_enhanced_stereo: Option<bool>,
            specified_stereo_query_matches_unspecified: Option<bool>,
            use_query_query_matches: Option<bool>,
            recursion_possible: Option<bool>,
            max_recursive_matches: Option<usize>,
            num_threads: Option<i32>,
            aromatic_matches_conjugated: Option<bool>,
            aromatic_matches_single_or_double: Option<bool>,
            extra_atom_check_overrides_default_check: Option<bool>,
            extra_bond_check_overrides_default_check: Option<bool>,
            use_generic_matchers: Option<bool>,
            atom_properties: Option<Vec<String>>,
            bond_properties: Option<Vec<String>>,
            final_match: Option<FinalMatch>,
            atom_match: Option<AtomMatch>,
            bond_match: Option<BondMatch>,
            extra_final_check: Option<FinalMatch>,
            extra_atom_check: Option<AtomMatch>,
            extra_bond_check: Option<BondMatch>,
        ) -> PyResult<Vec<MatchResult>> {
            let params = SubstructMatchParams::from_keywords(
                max_matches,
                uniquify,
                use_chirality,
                use_enhanced_stereo,
                specified_stereo_query_matches_unspecified,
                use_query_query_matches,
                recursion_possible,
                max_recursive_matches,
                num_threads,
                aromatic_matches_conjugated,
                aromatic_matches_single_or_double,
                extra_atom_check_overrides_default_check,
                extra_bond_check_overrides_default_check,
                use_generic_matchers,
                atom_properties,
                bond_properties,
                final_match,
                atom_match,
                bond_match,
                extra_final_check,
                extra_atom_check,
                extra_bond_check,
            )?;
            Self::_substruct_matches_params_native(slf, query, Some(params))
        }

        /// Return query matches in native order under the selected uniqueness, chirality and match-limit settings. Supply Params or keywords, not both. Uses the explicitly supplied parameter form.
        fn substruct_matches_with_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: &SubstructMatchParams,
        ) -> PyResult<Vec<MatchResult>> {
            Self::_substruct_matches_params(slf, query, Some(params))
        }

        #[pyo3(signature = (query, params=crate::arguments::Argument::Missing, *, max_matches=crate::arguments::Argument::Missing, uniquify=crate::arguments::Argument::Missing, use_chirality=crate::arguments::Argument::Missing, use_enhanced_stereo=crate::arguments::Argument::Missing, specified_stereo_query_matches_unspecified=crate::arguments::Argument::Missing, use_query_query_matches=crate::arguments::Argument::Missing, recursion_possible=crate::arguments::Argument::Missing, max_recursive_matches=crate::arguments::Argument::Missing, num_threads=crate::arguments::Argument::Missing, aromatic_matches_conjugated=crate::arguments::Argument::Missing, aromatic_matches_single_or_double=crate::arguments::Argument::Missing, extra_atom_check_overrides_default_check=crate::arguments::Argument::Missing, extra_bond_check_overrides_default_check=crate::arguments::Argument::Missing, use_generic_matchers=crate::arguments::Argument::Missing, atom_properties=crate::arguments::Argument::Missing, bond_properties=crate::arguments::Argument::Missing, final_match=crate::arguments::Argument::Missing, atom_match=crate::arguments::Argument::Missing, bond_match=crate::arguments::Argument::Missing, extra_final_check=crate::arguments::Argument::Missing, extra_atom_check=crate::arguments::Argument::Missing, extra_bond_check=crate::arguments::Argument::Missing))]
        /// Return whether the query has an accepted match. Supply SubstructMatchParams or keyword options, not both.
        fn has_substruct_match(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: crate::arguments::Argument<Option<PyRef<'_, SubstructMatchParams>>>,
            max_matches: crate::arguments::Argument<Option<usize>>,
            uniquify: crate::arguments::Argument<Option<bool>>,
            use_chirality: crate::arguments::Argument<Option<bool>>,
            use_enhanced_stereo: crate::arguments::Argument<Option<bool>>,
            specified_stereo_query_matches_unspecified: crate::arguments::Argument<Option<bool>>,
            use_query_query_matches: crate::arguments::Argument<Option<bool>>,
            recursion_possible: crate::arguments::Argument<Option<bool>>,
            max_recursive_matches: crate::arguments::Argument<Option<usize>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            aromatic_matches_conjugated: crate::arguments::Argument<Option<bool>>,
            aromatic_matches_single_or_double: crate::arguments::Argument<Option<bool>>,
            extra_atom_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            extra_bond_check_overrides_default_check: crate::arguments::Argument<Option<bool>>,
            use_generic_matchers: crate::arguments::Argument<Option<bool>>,
            atom_properties: crate::arguments::Argument<Option<Vec<String>>>,
            bond_properties: crate::arguments::Argument<Option<Vec<String>>>,
            final_match: crate::arguments::Argument<Option<FinalMatch>>,
            atom_match: crate::arguments::Argument<Option<AtomMatch>>,
            bond_match: crate::arguments::Argument<Option<BondMatch>>,
            extra_final_check: crate::arguments::Argument<Option<FinalMatch>>,
            extra_atom_check: crate::arguments::Argument<Option<AtomMatch>>,
            extra_bond_check: crate::arguments::Argument<Option<BondMatch>>,
        ) -> PyResult<bool> {
            let keyword_form = max_matches.is_supplied()
                || uniquify.is_supplied()
                || use_chirality.is_supplied()
                || use_enhanced_stereo.is_supplied()
                || specified_stereo_query_matches_unspecified.is_supplied()
                || use_query_query_matches.is_supplied()
                || recursion_possible.is_supplied()
                || max_recursive_matches.is_supplied()
                || num_threads.is_supplied()
                || aromatic_matches_conjugated.is_supplied()
                || aromatic_matches_single_or_double.is_supplied()
                || extra_atom_check_overrides_default_check.is_supplied()
                || extra_bond_check_overrides_default_check.is_supplied()
                || use_generic_matchers.is_supplied()
                || atom_properties.is_supplied()
                || bond_properties.is_supplied()
                || final_match.is_supplied()
                || atom_match.is_supplied()
                || bond_match.is_supplied()
                || extra_final_check.is_supplied()
                || extra_atom_check.is_supplied()
                || extra_bond_check.is_supplied();
            crate::arguments::configuration_arguments!(params; max_matches, uniquify, use_chirality, use_enhanced_stereo, specified_stereo_query_matches_unspecified, use_query_query_matches, recursion_possible, max_recursive_matches, num_threads, aromatic_matches_conjugated, aromatic_matches_single_or_double, extra_atom_check_overrides_default_check, extra_bond_check_overrides_default_check, use_generic_matchers, atom_properties, bond_properties, final_match, atom_match, bond_match, extra_final_check, extra_atom_check, extra_bond_check);
            let params = params.as_deref();
            if keyword_form {
                Self::_has_substruct_match_keywords(
                    slf,
                    query,
                    max_matches,
                    uniquify,
                    use_chirality,
                    use_enhanced_stereo,
                    specified_stereo_query_matches_unspecified,
                    use_query_query_matches,
                    recursion_possible,
                    max_recursive_matches,
                    num_threads,
                    aromatic_matches_conjugated,
                    aromatic_matches_single_or_double,
                    extra_atom_check_overrides_default_check,
                    extra_bond_check_overrides_default_check,
                    use_generic_matchers,
                    atom_properties,
                    bond_properties,
                    final_match,
                    atom_match,
                    bond_match,
                    extra_final_check,
                    extra_atom_check,
                    extra_bond_check,
                )
            } else {
                Self::_has_substruct_match_params(slf, query, params)
            }
        }

        #[doc = r#"{"role":"overload","target":"Molecule.has_substruct_match"}"#]
        #[pyo3(signature = (query, params=None))]
        fn _has_substruct_match_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: Option<&SubstructMatchParams>,
        ) -> PyResult<bool> {
            let py = slf.py();
            let params = params.map(|value| value.capture(py)).transpose()?;
            Self::_has_substruct_match_params_native(slf, query, params)
        }

        #[doc = r#"{"role":"overload","target":"Molecule.has_substruct_match"}"#]
        #[pyo3(signature = (query, *, max_matches=None, uniquify=None, use_chirality=None, use_enhanced_stereo=None, specified_stereo_query_matches_unspecified=None, use_query_query_matches=None, recursion_possible=None, max_recursive_matches=None, num_threads=None, aromatic_matches_conjugated=None, aromatic_matches_single_or_double=None, extra_atom_check_overrides_default_check=None, extra_bond_check_overrides_default_check=None, use_generic_matchers=None, atom_properties=None, bond_properties=None, final_match=None, atom_match=None, bond_match=None, extra_final_check=None, extra_atom_check=None, extra_bond_check=None))]
        fn _has_substruct_match_keywords(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            max_matches: Option<usize>,
            uniquify: Option<bool>,
            use_chirality: Option<bool>,
            use_enhanced_stereo: Option<bool>,
            specified_stereo_query_matches_unspecified: Option<bool>,
            use_query_query_matches: Option<bool>,
            recursion_possible: Option<bool>,
            max_recursive_matches: Option<usize>,
            num_threads: Option<i32>,
            aromatic_matches_conjugated: Option<bool>,
            aromatic_matches_single_or_double: Option<bool>,
            extra_atom_check_overrides_default_check: Option<bool>,
            extra_bond_check_overrides_default_check: Option<bool>,
            use_generic_matchers: Option<bool>,
            atom_properties: Option<Vec<String>>,
            bond_properties: Option<Vec<String>>,
            final_match: Option<FinalMatch>,
            atom_match: Option<AtomMatch>,
            bond_match: Option<BondMatch>,
            extra_final_check: Option<FinalMatch>,
            extra_atom_check: Option<AtomMatch>,
            extra_bond_check: Option<BondMatch>,
        ) -> PyResult<bool> {
            let params = SubstructMatchParams::from_keywords(
                max_matches,
                uniquify,
                use_chirality,
                use_enhanced_stereo,
                specified_stereo_query_matches_unspecified,
                use_query_query_matches,
                recursion_possible,
                max_recursive_matches,
                num_threads,
                aromatic_matches_conjugated,
                aromatic_matches_single_or_double,
                extra_atom_check_overrides_default_check,
                extra_bond_check_overrides_default_check,
                use_generic_matchers,
                atom_properties,
                bond_properties,
                final_match,
                atom_match,
                bond_match,
                extra_final_check,
                extra_atom_check,
                extra_bond_check,
            )?;
            Self::_has_substruct_match_params_native(slf, query, Some(params))
        }

        /// Return whether the query has an accepted match. Supply SubstructMatchParams or keyword options, not both. Uses the explicitly supplied parameter form.
        fn has_substruct_match_with_params(
            slf: Bound<'_, Self>,
            query: &QueryGraph,
            params: &SubstructMatchParams,
        ) -> PyResult<bool> {
            Self::_has_substruct_match_params(slf, query, Some(params))
        }

        /// Match an already compiled query without recompiling its query representation.
        fn substruct_matches_compiled(
            &self,
            py: Python<'_>,
            query: &CompiledQuery,
        ) -> PyResult<Vec<MatchResult>> {
            let native = &self.inner;
            let query = query.inner();
            py.detach(|| native.substruct_matches_compiled(query))
                .map(|values| {
                    values
                        .into_iter()
                        .map(|inner| MatchResult { inner })
                        .collect()
                })
                .map_err(MatchError::from_rust)
        }

        #[staticmethod]
        /// Decode a COSMolKit molecular archive, preserving the serialized molecular state. Malformed or unsupported data raises its native archive error.
        pub(crate) fn from_binary(data: &[u8]) -> PyResult<Self> {
            ck::Molecule::from_binary(data)
                .map(|inner| Self { inner })
                .map_err(PickleError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Serialize the molecule as MDL MolBlock text with the selected writer and coordinate options; does not write a file.
        pub(crate) fn to_mol(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<String> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_mol_keywords(
                    py,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_to_mol_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_mol"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_mol_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_mol_with_params(&params),
                None => self.inner.to_mol(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_mol"}"#]
        #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_to_mol_keywords(
            &self,
            py: Python<'_>,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_to_mol_params(py, Some(&params))
        }
        /// Serialize the molecule as MDL MolBlock text with the selected writer and coordinate options; does not write a file. Uses explicit writer parameters.
        pub(crate) fn to_mol_with_params(
            &self,
            py: Python<'_>,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<String> {
            self._io_to_mol_params(py, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Serialize the molecule as SDF text with the selected writer and coordinate options; does not write a file.
        pub(crate) fn to_sdf(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<String> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_sdf_keywords(
                    py,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_to_sdf_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_sdf_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_sdf_with_params(&params),
                None => self.inner.to_sdf(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf"}"#]
        #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_to_sdf_keywords(
            &self,
            py: Python<'_>,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_to_sdf_params(py, Some(&params))
        }
        /// Serialize the molecule as SDF text with the selected writer and coordinate options; does not write a file. Uses explicit writer parameters.
        pub(crate) fn to_sdf_with_params(
            &self,
            py: Python<'_>,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<String> {
            self._io_to_sdf_params(py, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Serialize SDF using the molecule's 2D coordinate selection and the supplied writer options. Returns text without writing a file.
        pub(crate) fn to_sdf_2d(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<String> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_sdf_2d_keywords(
                    py,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_to_sdf_2d_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf_2d"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_sdf_2d_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_sdf_2d_with_params(&params),
                None => self.inner.to_sdf_2d(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf_2d"}"#]
        #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_to_sdf_2d_keywords(
            &self,
            py: Python<'_>,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_to_sdf_2d_params(py, Some(&params))
        }
        /// Serialize SDF using the molecule's 2D coordinate selection and the supplied writer options. Returns text without writing a file.
        pub(crate) fn to_sdf_2d_with_params(
            &self,
            py: Python<'_>,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<String> {
            self._io_to_sdf_2d_params(py, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Serialize SDF using the molecule's 3D coordinate selection and the supplied writer options. Returns text without writing a file.
        pub(crate) fn to_sdf_3d(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<String> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_sdf_3d_keywords(
                    py,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_to_sdf_3d_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf_3d"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_sdf_3d_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_sdf_3d_with_params(&params),
                None => self.inner.to_sdf_3d(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_sdf_3d"}"#]
        #[pyo3(signature=(*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_to_sdf_3d_keywords(
            &self,
            py: Python<'_>,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_to_sdf_3d_params(py, Some(&params))
        }
        /// Serialize SDF using the molecule's 3D coordinate selection and the supplied writer options. Returns text without writing a file.
        pub(crate) fn to_sdf_3d_with_params(
            &self,
            py: Python<'_>,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<String> {
            self._io_to_sdf_3d_params(py, Some(params))
        }
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Write the molecule as MDL MolBlock to the supplied filesystem path with the selected writer options. IO errors propagate.
        pub(crate) fn write_mol(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_write_mol_keywords(
                    py,
                    path,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_write_mol_params(py, path, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_mol"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_write_mol_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => self.inner.write_mol_with_params(path_text, &params),
                None => self.inner.write_mol(path_text),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_mol"}"#]
        #[pyo3(signature=(path,*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_write_mol_keywords(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<()> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_write_mol_params(py, path, Some(&params))
        }
        /// Write the molecule as MDL MolBlock to the supplied filesystem path with the selected writer options. IO errors propagate. Uses explicit writer parameters.
        pub(crate) fn write_mol_with_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<()> {
            self._io_write_mol_params(py, path, Some(params))
        }
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Write the molecule as SDF to the supplied filesystem path with the selected writer options. IO errors propagate.
        pub(crate) fn write_sdf(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_write_sdf_keywords(
                    py,
                    path,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_write_sdf_params(py, path, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_sdf"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_write_sdf_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => self.inner.write_sdf_with_params(path_text, &params),
                None => self.inner.write_sdf(path_text),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_sdf"}"#]
        #[pyo3(signature=(path,*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_write_sdf_keywords(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<()> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_write_sdf_params(py, path, Some(&params))
        }
        /// Write the molecule as SDF to the supplied filesystem path with the selected writer options. IO errors propagate. Uses explicit writer parameters.
        pub(crate) fn write_sdf_with_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<()> {
            self._io_write_sdf_params(py, path, Some(params))
        }
        #[pyo3(signature=(directory,file_name=None,params=crate::arguments::Argument::Missing,*,format=crate::arguments::Argument::Missing,force_2d=crate::arguments::Argument::Missing,include_stereo=crate::arguments::Argument::Missing,kekulize=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing,coordinate_selection=crate::arguments::Argument::Missing,include_coordinates=crate::arguments::Argument::Missing))]
        /// Write native SDF file outputs under the requested coordinate and writer policy, preserving the canonical file naming and error behavior.
        pub(crate) fn write_sdf_files(
            &self,
            py: Python<'_>,
            directory: std::path::PathBuf,
            file_name: Option<&str>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::MolBlockWriteParams>>,
            >,
            format: crate::arguments::Argument<Option<crate::io_batch_params::SdfFormatInput>>,
            force_2d: crate::arguments::Argument<Option<bool>>,
            include_stereo: crate::arguments::Argument<Option<bool>>,
            kekulize: crate::arguments::Argument<Option<bool>>,
            precision: crate::arguments::Argument<Option<usize>>,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            include_coordinates: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<std::path::PathBuf> {
            let keyword_form = format.is_supplied()
                || force_2d.is_supplied()
                || include_stereo.is_supplied()
                || kekulize.is_supplied()
                || precision.is_supplied()
                || coordinate_selection.is_supplied()
                || include_coordinates.is_supplied();
            crate::arguments::configuration_arguments!(params; format, force_2d, include_stereo, kekulize, precision, coordinate_selection, include_coordinates);
            let params = params.as_deref();
            if keyword_form {
                self._io_write_sdf_files_keywords(
                    py,
                    directory,
                    file_name,
                    format,
                    force_2d,
                    include_stereo,
                    kekulize,
                    precision,
                    coordinate_selection,
                    include_coordinates,
                )
            } else {
                self._io_write_sdf_files_params(py, directory, file_name, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_sdf_files"}"#]
        #[pyo3(signature=(directory,file_name=None,params=None))]
        fn _io_write_sdf_files_params(
            &self,
            py: Python<'_>,
            directory: std::path::PathBuf,
            file_name: Option<&str>,
            params: Option<&crate::io_batch_params::MolBlockWriteParams>,
        ) -> PyResult<std::path::PathBuf> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&directory)?;
            py.detach(|| match params {
                Some(params) => self
                    .inner
                    .write_sdf_files_with_params(path_text, file_name, &params),
                None => self.inner.write_sdf_files(path_text, file_name),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_sdf_files"}"#]
        #[pyo3(signature=(directory,file_name=None,*,format=None,force_2d=None,include_stereo=None,kekulize=None,precision=None,coordinate_selection=None,include_coordinates=None))]
        fn _io_write_sdf_files_keywords(
            &self,
            py: Python<'_>,
            directory: std::path::PathBuf,
            file_name: Option<&str>,
            format: Option<crate::io_batch_params::SdfFormatInput>,
            force_2d: Option<bool>,
            include_stereo: Option<bool>,
            kekulize: Option<bool>,
            precision: Option<usize>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            include_coordinates: Option<bool>,
        ) -> PyResult<std::path::PathBuf> {
            let params = crate::molecular_io_calls::mol_write_keywords(
                format,
                force_2d,
                include_stereo,
                kekulize,
                precision,
                coordinate_selection,
                include_coordinates,
            )?;
            self._io_write_sdf_files_params(py, directory, file_name, Some(&params))
        }
        /// Write native SDF file outputs under the requested coordinate and writer policy, preserving the canonical file naming and error behavior.
        pub(crate) fn write_sdf_files_with_params(
            &self,
            py: Python<'_>,
            directory: std::path::PathBuf,
            file_name: Option<&str>,
            params: &crate::io_batch_params::MolBlockWriteParams,
        ) -> PyResult<std::path::PathBuf> {
            self._io_write_sdf_files_params(py, directory, file_name, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,coordinate_selection=crate::arguments::Argument::Missing,all_conformers=crate::arguments::Argument::Missing,flavor=crate::arguments::Argument::Missing))]
        /// Serialize the molecule as PDB text with the selected writer and coordinate options; does not write a file.
        pub(crate) fn to_pdb(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::PdbWriteParams>>,
            >,
            coordinate_selection: crate::arguments::Argument<
                Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            >,
            all_conformers: crate::arguments::Argument<Option<bool>>,
            flavor: crate::arguments::Argument<Option<crate::pdb_flags::PdbWriteFlags>>,
        ) -> PyResult<String> {
            let keyword_form = coordinate_selection.is_supplied()
                || all_conformers.is_supplied()
                || flavor.is_supplied();
            crate::arguments::configuration_arguments!(params; coordinate_selection, all_conformers, flavor);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_pdb_keywords(py, coordinate_selection, all_conformers, flavor)
            } else {
                self._io_to_pdb_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_pdb"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_pdb_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::io_batch_params::PdbWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_pdb_with_params(&params),
                None => self.inner.to_pdb(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_pdb"}"#]
        #[pyo3(signature=(*,coordinate_selection=None,all_conformers=None,flavor=None))]
        fn _io_to_pdb_keywords(
            &self,
            py: Python<'_>,
            coordinate_selection: Option<crate::io_batch_params::MolCoordinateSelectionInput>,
            all_conformers: Option<bool>,
            flavor: Option<crate::pdb_flags::PdbWriteFlags>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::pdb_write_keywords(
                coordinate_selection,
                all_conformers,
                flavor,
            )?;
            self._io_to_pdb_params(py, Some(&params))
        }
        /// Serialize the molecule as PDB text with the selected writer and coordinate options; does not write a file. Uses explicit writer parameters.
        pub(crate) fn to_pdb_with_params(
            &self,
            py: Python<'_>,
            params: &crate::io_batch_params::PdbWriteParams,
        ) -> PyResult<String> {
            self._io_to_pdb_params(py, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,conformer_id=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing))]
        /// Serialize the molecule as XYZ text with the selected writer and coordinate options; does not write a file.
        pub(crate) fn to_xyz(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::params::XyzWriteParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            precision: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<String> {
            let keyword_form = conformer_id.is_supplied() || precision.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, precision);
            let params = params.as_deref();
            if keyword_form {
                self._io_to_xyz_keywords(py, conformer_id, precision)
            } else {
                self._io_to_xyz_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_xyz"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_to_xyz_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::params::XyzWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.to_xyz_with_params(&params),
                None => self.inner.to_xyz(),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_xyz"}"#]
        #[pyo3(signature=(*,conformer_id=None,precision=None))]
        fn _io_to_xyz_keywords(
            &self,
            py: Python<'_>,
            conformer_id: Option<usize>,
            precision: Option<u32>,
        ) -> PyResult<String> {
            let params = crate::molecular_io_calls::xyz_write_keywords(conformer_id, precision)?;
            self._io_to_xyz_params(py, Some(&params))
        }
        /// Serialize the molecule as XYZ text with the selected writer and coordinate options; does not write a file. Uses explicit writer parameters.
        pub(crate) fn to_xyz_with_params(
            &self,
            py: Python<'_>,
            params: &crate::params::XyzWriteParams,
        ) -> PyResult<String> {
            self._io_to_xyz_params(py, Some(params))
        }
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,conformer_id=crate::arguments::Argument::Missing,precision=crate::arguments::Argument::Missing))]
        /// Write the molecule as XYZ to the supplied filesystem path with the selected writer options. IO errors propagate.
        pub(crate) fn write_xyz(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::params::XyzWriteParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            precision: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<()> {
            let keyword_form = conformer_id.is_supplied() || precision.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, precision);
            let params = params.as_deref();
            if keyword_form {
                self._io_write_xyz_keywords(py, path, conformer_id, precision)
            } else {
                self._io_write_xyz_params(py, path, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_xyz"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_write_xyz_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::params::XyzWriteParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => self.inner.write_xyz_with_params(path_text, &params),
                None => self.inner.write_xyz(path_text),
            })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.write_xyz"}"#]
        #[pyo3(signature=(path,*,conformer_id=None,precision=None))]
        fn _io_write_xyz_keywords(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            conformer_id: Option<usize>,
            precision: Option<u32>,
        ) -> PyResult<()> {
            let params = crate::molecular_io_calls::xyz_write_keywords(conformer_id, precision)?;
            self._io_write_xyz_params(py, path, Some(&params))
        }
        /// Write the molecule as XYZ to the supplied filesystem path with the selected writer options. IO errors propagate. Uses explicit writer parameters.
        pub(crate) fn write_xyz_with_params(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::params::XyzWriteParams,
        ) -> PyResult<()> {
            self._io_write_xyz_params(py, path, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,flavor=crate::arguments::Argument::Missing,proximity_bonding=crate::arguments::Argument::Missing))]
        /// Parse PDB text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_pdb(
            py: Python<'_>,
            text: &str,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioMoleculeParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            flavor: crate::arguments::Argument<Option<crate::pdb_flags::PdbReadFlags>>,
            proximity_bonding: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || flavor.is_supplied()
                || proximity_bonding.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, flavor, proximity_bonding);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_from_pdb_keywords(
                    py,
                    text,
                    sanitize,
                    remove_hs,
                    flavor,
                    proximity_bonding,
                )
            } else {
                Self::_io_from_pdb_params(py, text, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_pdb"}"#]
        #[pyo3(signature=(text,params=None))]
        pub(super) fn _io_from_pdb_params(
            py: Python<'_>,
            text: &str,
            params: Option<&crate::params::BioMoleculeParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => ck::Molecule::from_pdb_with_params(text, &params),
                None => ck::Molecule::from_pdb(text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_pdb"}"#]
        #[pyo3(signature=(text,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
        pub(super) fn _io_from_pdb_keywords(
            py: Python<'_>,
            text: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            flavor: Option<crate::pdb_flags::PdbReadFlags>,
            proximity_bonding: Option<bool>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::bio_keywords(
                sanitize,
                remove_hs,
                flavor,
                proximity_bonding,
            )?;
            Self::_io_from_pdb_params(py, text, Some(&params))
        }
        #[staticmethod]
        /// Parse PDB text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path. Uses the explicit reader parameters.
        pub(crate) fn from_pdb_with_params(
            py: Python<'_>,
            text: &str,
            params: &crate::params::BioMoleculeParams,
        ) -> PyResult<Self> {
            Self::_io_from_pdb_params(py, text, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,flavor=crate::arguments::Argument::Missing,proximity_bonding=crate::arguments::Argument::Missing))]
        /// Parse mmCIF text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_mmcif(
            py: Python<'_>,
            text: &str,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::params::BioMoleculeParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            flavor: crate::arguments::Argument<Option<crate::pdb_flags::PdbReadFlags>>,
            proximity_bonding: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || flavor.is_supplied()
                || proximity_bonding.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, flavor, proximity_bonding);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_from_mmcif_keywords(
                    py,
                    text,
                    sanitize,
                    remove_hs,
                    flavor,
                    proximity_bonding,
                )
            } else {
                Self::_io_from_mmcif_params(py, text, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mmcif"}"#]
        #[pyo3(signature=(text,params=None))]
        pub(super) fn _io_from_mmcif_params(
            py: Python<'_>,
            text: &str,
            params: Option<&crate::params::BioMoleculeParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => ck::Molecule::from_mmcif_with_params(text, &params),
                None => ck::Molecule::from_mmcif(text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mmcif"}"#]
        #[pyo3(signature=(text,*,sanitize=None,remove_hs=None,flavor=None,proximity_bonding=None))]
        pub(super) fn _io_from_mmcif_keywords(
            py: Python<'_>,
            text: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            flavor: Option<crate::pdb_flags::PdbReadFlags>,
            proximity_bonding: Option<bool>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::bio_keywords(
                sanitize,
                remove_hs,
                flavor,
                proximity_bonding,
            )?;
            Self::_io_from_mmcif_params(py, text, Some(&params))
        }
        #[staticmethod]
        /// Parse mmCIF text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path. Uses the explicit reader parameters.
        pub(crate) fn from_mmcif_with_params(
            py: Python<'_>,
            text: &str,
            params: &crate::params::BioMoleculeParams,
        ) -> PyResult<Self> {
            Self::_io_from_mmcif_params(py, text, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(text,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,variant=crate::arguments::Argument::Missing,cleanup_substructures=crate::arguments::Argument::Missing))]
        /// Parse Tripos MOL2 text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_mol2(
            py: Python<'_>,
            text: &str,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::Mol2ReadParams>>,
            >,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            variant: crate::arguments::Argument<Option<crate::io_batch_params::Mol2TypeInput>>,
            cleanup_substructures: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || variant.is_supplied()
                || cleanup_substructures.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, variant, cleanup_substructures);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_from_mol2_keywords(
                    py,
                    text,
                    sanitize,
                    remove_hs,
                    variant,
                    cleanup_substructures,
                )
            } else {
                Self::_io_from_mol2_params(py, text, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mol2"}"#]
        #[pyo3(signature=(text,params=None))]
        pub(super) fn _io_from_mol2_params(
            py: Python<'_>,
            text: &str,
            params: Option<&crate::io_batch_params::Mol2ReadParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => ck::Molecule::from_mol2_with_params(text, &params),
                None => ck::Molecule::from_mol2(text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_mol2"}"#]
        #[pyo3(signature=(text,*,sanitize=None,remove_hs=None,variant=None,cleanup_substructures=None))]
        pub(super) fn _io_from_mol2_keywords(
            py: Python<'_>,
            text: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            variant: Option<crate::io_batch_params::Mol2TypeInput>,
            cleanup_substructures: Option<bool>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::mol2_read_keywords(
                sanitize,
                remove_hs,
                variant,
                cleanup_substructures,
            )?;
            Self::_io_from_mol2_params(py, text, Some(&params))
        }
        #[staticmethod]
        /// Parse Tripos MOL2 text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path. Uses the explicit reader parameters.
        pub(crate) fn from_mol2_with_params(
            py: Python<'_>,
            text: &str,
            params: &crate::io_batch_params::Mol2ReadParams,
        ) -> PyResult<Self> {
            Self::_io_from_mol2_params(py, text, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,variant=crate::arguments::Argument::Missing,cleanup_substructures=crate::arguments::Argument::Missing))]
        /// Read a molecule from a Tripos MOL2 file at the supplied path using the canonical reader. IO and parsing errors propagate.
        pub(crate) fn read_mol2(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::io_batch_params::Mol2ReadParams>>,
            >,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            variant: crate::arguments::Argument<Option<crate::io_batch_params::Mol2TypeInput>>,
            cleanup_substructures: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || variant.is_supplied()
                || cleanup_substructures.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, variant, cleanup_substructures);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_read_mol2_keywords(
                    py,
                    path,
                    sanitize,
                    remove_hs,
                    variant,
                    cleanup_substructures,
                )
            } else {
                Self::_io_read_mol2_params(py, path, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_mol2"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_read_mol2_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::io_batch_params::Mol2ReadParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => ck::Molecule::read_mol2_with_params(path_text, &params),
                None => ck::Molecule::read_mol2(path_text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_mol2"}"#]
        #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,variant=None,cleanup_substructures=None))]
        fn _io_read_mol2_keywords(
            py: Python<'_>,
            path: std::path::PathBuf,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            variant: Option<crate::io_batch_params::Mol2TypeInput>,
            cleanup_substructures: Option<bool>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::mol2_read_keywords(
                sanitize,
                remove_hs,
                variant,
                cleanup_substructures,
            )?;
            Self::_io_read_mol2_params(py, path, Some(&params))
        }
        #[staticmethod]
        /// Read a molecule from a Tripos MOL2 file at the supplied path using the canonical reader. IO and parsing errors propagate. Uses explicit reader parameters.
        pub(crate) fn read_mol2_with_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::io_batch_params::Mol2ReadParams,
        ) -> PyResult<Self> {
            Self::_io_read_mol2_params(py, path, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
        /// Read a molecule from a MDL MolBlock file at the supplied path using the canonical reader. IO and parsing errors propagate.
        pub(crate) fn read_mol(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::selectors::SdfReadParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            strict_parsing: crate::arguments::Argument<Option<bool>>,
            expand_attachment_points: crate::arguments::Argument<Option<bool>>,
            process_property_lists: crate::arguments::Argument<Option<bool>>,
            coordinate_mode: crate::arguments::Argument<
                Option<crate::selectors::SdfCoordinateModeInput>,
            >,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || strict_parsing.is_supplied()
                || expand_attachment_points.is_supplied()
                || process_property_lists.is_supplied()
                || coordinate_mode.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_read_mol_keywords(
                    py,
                    path,
                    sanitize,
                    remove_hs,
                    strict_parsing,
                    expand_attachment_points,
                    process_property_lists,
                    coordinate_mode,
                )
            } else {
                Self::_io_read_mol_params(py, path, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_mol"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_read_mol_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::selectors::SdfReadParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => ck::Molecule::read_mol_with_params(path_text, &params),
                None => ck::Molecule::read_mol(path_text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_mol"}"#]
        #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
        fn _io_read_mol_keywords(
            py: Python<'_>,
            path: std::path::PathBuf,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            strict_parsing: Option<bool>,
            expand_attachment_points: Option<bool>,
            process_property_lists: Option<bool>,
            coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::sdf_keywords(
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )?;
            Self::_io_read_mol_params(py, path, Some(&params))
        }
        #[staticmethod]
        /// Read a molecule from a MDL MolBlock file at the supplied path using the canonical reader. IO and parsing errors propagate. Uses explicit reader parameters.
        pub(crate) fn read_mol_with_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::selectors::SdfReadParams,
        ) -> PyResult<Self> {
            Self::_io_read_mol_params(py, path, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
        /// Read a molecule from a SDF file at the supplied path using the canonical reader. IO and parsing errors propagate.
        pub(crate) fn read_sdf(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: crate::arguments::Argument<Option<PyRef<'_, crate::selectors::SdfReadParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
            strict_parsing: crate::arguments::Argument<Option<bool>>,
            expand_attachment_points: crate::arguments::Argument<Option<bool>>,
            process_property_lists: crate::arguments::Argument<Option<bool>>,
            coordinate_mode: crate::arguments::Argument<
                Option<crate::selectors::SdfCoordinateModeInput>,
            >,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied()
                || remove_hs.is_supplied()
                || strict_parsing.is_supplied()
                || expand_attachment_points.is_supplied()
                || process_property_lists.is_supplied()
                || coordinate_mode.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
            let params = params.as_deref();
            if keyword_form {
                Self::_io_read_sdf_keywords(
                    py,
                    path,
                    sanitize,
                    remove_hs,
                    strict_parsing,
                    expand_attachment_points,
                    process_property_lists,
                    coordinate_mode,
                )
            } else {
                Self::_io_read_sdf_params(py, path, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_sdf"}"#]
        #[pyo3(signature=(path,params=None))]
        fn _io_read_sdf_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: Option<&crate::selectors::SdfReadParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let path_text = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| match params {
                Some(params) => ck::Molecule::read_sdf_with_params(path_text, &params),
                None => ck::Molecule::read_sdf(path_text),
            })
            .map(|inner| Self { inner })
            .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.read_sdf"}"#]
        #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
        fn _io_read_sdf_keywords(
            py: Python<'_>,
            path: std::path::PathBuf,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
            strict_parsing: Option<bool>,
            expand_attachment_points: Option<bool>,
            process_property_lists: Option<bool>,
            coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
        ) -> PyResult<Self> {
            let params = crate::molecular_io_calls::sdf_keywords(
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )?;
            Self::_io_read_sdf_params(py, path, Some(&params))
        }
        #[staticmethod]
        /// Read a molecule from a SDF file at the supplied path using the canonical reader. IO and parsing errors propagate. Uses explicit reader parameters.
        pub(crate) fn read_sdf_with_params(
            py: Python<'_>,
            path: std::path::PathBuf,
            params: &crate::selectors::SdfReadParams,
        ) -> PyResult<Self> {
            Self::_io_read_sdf_params(py, path, Some(params))
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,n_bits=crate::arguments::Argument::Missing,is_query=crate::arguments::Argument::Missing,bit_flags=crate::arguments::Argument::Missing))]
        /// Calculate the Avalon fingerprint using the selected native Avalon options.
        pub(crate) fn fingerprint_avalon(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, crate::numerical_params::AvalonFingerprintParams>>,
            >,
            n_bits: crate::arguments::Argument<Option<u32>>,
            is_query: crate::arguments::Argument<Option<bool>>,
            bit_flags: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<crate::native::Fingerprint> {
            let keyword_form =
                n_bits.is_supplied() || is_query.is_supplied() || bit_flags.is_supplied();
            crate::arguments::configuration_arguments!(params; n_bits, is_query, bit_flags);
            let params = params.as_deref();
            if keyword_form {
                self._io_fingerprint_avalon_keywords(py, n_bits, is_query, bit_flags)
            } else {
                self._io_fingerprint_avalon_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_avalon"}"#]
        #[pyo3(signature=(params=None))]
        fn _io_fingerprint_avalon_params(
            &self,
            py: Python<'_>,
            params: Option<&crate::numerical_params::AvalonFingerprintParams>,
        ) -> PyResult<crate::native::Fingerprint> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            py.detach(|| match params {
                Some(params) => self.inner.fingerprint_avalon_with_params(&params),
                None => self.inner.fingerprint_avalon(),
            })
            .map(|inner| crate::native::Fingerprint { inner })
            .map_err(crate::molecular_io_errors::AvalonFingerprintError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.fingerprint_avalon"}"#]
        #[pyo3(signature=(*,n_bits=None,is_query=None,bit_flags=None))]
        fn _io_fingerprint_avalon_keywords(
            &self,
            py: Python<'_>,
            n_bits: Option<u32>,
            is_query: Option<bool>,
            bit_flags: Option<u32>,
        ) -> PyResult<crate::native::Fingerprint> {
            let params = crate::molecular_io_calls::avalon_keywords(n_bits, is_query, bit_flags)?;
            self._io_fingerprint_avalon_params(py, Some(&params))
        }
        /// Calculate the Avalon fingerprint using the selected native Avalon options. Uses the explicit configuration argument.
        pub(crate) fn fingerprint_avalon_with_params(
            &self,
            py: Python<'_>,
            params: &crate::numerical_params::AvalonFingerprintParams,
        ) -> PyResult<crate::native::Fingerprint> {
            self._io_fingerprint_avalon_params(py, Some(params))
        }
        #[staticmethod]
        /// Parse XYZ text into a new molecule using the canonical reader and selected input policy. This accepts content, not a filesystem path.
        pub(crate) fn from_xyz_block(py: Python<'_>, text: &str) -> PyResult<Self> {
            py.detach(|| ck::Molecule::from_xyz_block(text))
                .map(|inner| Self { inner })
                .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }
        #[staticmethod]
        /// Read a molecule from a XYZ file at the supplied path using the canonical reader. IO and parsing errors propagate.
        pub(crate) fn read_xyz(py: Python<'_>, path: std::path::PathBuf) -> PyResult<Self> {
            let value = crate::molecular_io_calls::path_text(&path)?;
            py.detach(|| ck::Molecule::read_xyz(value))
                .map(|inner| Self { inner })
                .map_err(crate::molecular_io_errors::MolecularIoError::from_rust)
        }

        /// Render a molecular SVG image at the supplied width and height, returning its XML text.
        fn to_svg(&self, py: Python<'_>, width: u32, height: u32) -> PyResult<String> {
            py.detach(|| self.inner.to_svg(width, height))
                .map_err(crate::drawing_errors::DrawingError::from_rust)
        }
        /// Render a molecular PNG image at the supplied width and height, returning PNG bytes.
        fn to_png<'py>(
            &self,
            py: Python<'py>,
            width: u32,
            height: u32,
        ) -> PyResult<Bound<'py, pyo3::types::PyBytes>> {
            let bytes = py
                .detach(|| self.inner.to_png(width, height))
                .map_err(crate::drawing_errors::DrawingError::from_rust)?;
            Ok(pyo3::types::PyBytes::new(py, &bytes))
        }
        /// Render a molecular SVG image and write it to the supplied filesystem path.
        fn write_svg(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            width: u32,
            height: u32,
        ) -> PyResult<()> {
            py.detach(|| self.inner.write_svg(&path, width, height))
                .map_err(crate::drawing_errors::DrawingWriteError::from_rust)
        }
        /// Render a molecular PNG image and write it to the supplied filesystem path.
        fn write_png(
            &self,
            py: Python<'_>,
            path: std::path::PathBuf,
            width: u32,
            height: u32,
        ) -> PyResult<()> {
            py.detach(|| self.inner.write_png(&path, width, height))
                .map_err(crate::drawing_errors::DrawingWriteError::from_rust)
        }

        /// Serialize the molecular value to COSMolKit archive bytes for later from_binary or pickle reconstruction.
        fn to_binary<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, pyo3::types::PyBytes>> {
            self.inner
                .to_binary()
                .map(|data| pyo3::types::PyBytes::new(py, &data))
                .map_err(PickleError::from_rust)
        }

        fn __reduce__<'py>(
            &self,
            py: Python<'py>,
        ) -> PyResult<(
            crate::binary_archive::MoleculeRebuild,
            (Bound<'py, pyo3::types::PyBytes>,),
        )> {
            let rebuild = py
                .import("cosmolkit")?
                .getattr("_molecule_from_binary")?
                .extract()?;
            Ok((rebuild, (self.to_binary(py)?,)))
        }

        fn __reduce_ex__<'py>(
            &self,
            py: Python<'py>,
            protocol: crate::binary_archive::PickleProtocol,
        ) -> PyResult<(
            crate::binary_archive::MoleculeRebuild,
            (Bound<'py, pyo3::types::PyBytes>,),
        )> {
            let _ = protocol;
            self.__reduce__(py)
        }
        #[getter]
        /// Atoms in index order, each retaining this molecular version for contextual reads.
        fn atoms(&self) -> Vec<MoleculeAtom> {
            self.inner
                .atoms()
                .into_iter()
                .map(|inner| MoleculeAtom { inner })
                .collect()
        }
        #[getter]
        /// Detached bond values in bond-index order.
        fn bonds(&self) -> Vec<Bond> {
            self.inner
                .bonds()
                .iter()
                .cloned()
                .map(|inner| Bond { inner })
                .collect()
        }
        /// Return the atom at a zero-based index, or None when the index is absent.
        /// The returned atom retains this molecular version; later edits to the
        /// source molecule do not change its topology or contextual properties.
        fn atom(&self, atom_id: usize) -> Option<MoleculeAtom> {
            self.inner
                .atom(ck::AtomId::new(atom_id))
                .map(|inner| MoleculeAtom { inner })
        }
        /// Return a detached bond value at a zero-based index, or None when the index is absent.
        fn bond(&self, bond_id: usize) -> Option<Bond> {
            self.inner
                .bond(ck::BondId::new(bond_id))
                .cloned()
                .map(|inner| Bond { inner })
        }
        #[getter]
        /// Whether a complete molecular CIP assignment is recorded. A partial assignment does not imply this flag.
        fn cip_computed(&self) -> PyResult<bool> {
            self.inner
                .cip_computed()
                .map_err(super::property_values::PropertyValueError::from_rust)
        }
        #[getter]
        /// Return the molecule-level property value collection.
        fn properties(&self) -> MoleculeProperties {
            MoleculeProperties {
                inner: self.inner.properties().clone(),
            }
        }
        #[getter]
        /// Return the detached topology value; modifying another molecule cannot invalidate this read.
        fn topology(&self) -> TopologyBlock {
            TopologyBlock::from_molecule(&self.inner)
        }
        /// Create a molecule builder initialized from this molecular value.
        fn to_builder(&self) -> MoleculeBuilder {
            MoleculeBuilder {
                inner: self.inner.to_builder(),
            }
        }
        /// Return whether MMFF parameters are available for every required atom in the molecule.
        fn mmff_has_all_molecule_params(&self) -> PyResult<bool> {
            self.inner
                .mmff_has_all_molecule_params()
                .map_err(super::mmff_properties::MmffMolPropertiesError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,mmff_variant=crate::arguments::Argument::Missing))]
        /// Prepare MMFF atom typing and molecular properties with the selected MMFF variant.
        fn mmff_properties(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, MmffPropertiesParams>>>,
            mmff_variant: crate::arguments::Argument<Option<crate::leaf_values::MmffVariantInput>>,
        ) -> PyResult<super::mmff_properties::MmffProperties> {
            let keyword_form = mmff_variant.is_supplied();
            crate::arguments::configuration_arguments!(params; mmff_variant);
            let params = params.as_deref();
            if keyword_form {
                self._mmff_properties_keywords(mmff_variant)
            } else {
                self._mmff_properties_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_properties"}"#]
        #[pyo3(signature=(params=None))]
        fn _mmff_properties_params(
            &self,
            params: Option<&MmffPropertiesParams>,
        ) -> PyResult<super::mmff_properties::MmffProperties> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let result = match params {
                Some(params) => self.inner.mmff_properties_with_params(&params),
                None => self.inner.mmff_properties(),
            };
            result
                .map(|inner| super::mmff_properties::MmffProperties { inner })
                .map_err(super::mmff_properties::MmffMolPropertiesError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_properties"}"#]
        #[pyo3(signature=(*,mmff_variant=None))]
        fn _mmff_properties_keywords(
            &self,
            mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
        ) -> PyResult<super::mmff_properties::MmffProperties> {
            if let Some(mmff_variant) = mmff_variant {
                self._mmff_properties_params(Some(&MmffPropertiesParams {
                    inner: crate::parameters::NativeField::Owned(ck::MmffPropertiesParams {
                        mmff_variant: mmff_variant.0,
                    }),
                }))
            } else {
                self._mmff_properties_params(None)
            }
        }
        /// Prepare MMFF atom typing and molecular properties with the selected MMFF variant. Uses the explicitly supplied parameter form.
        fn mmff_properties_with_params(
            &self,
            params: &MmffPropertiesParams,
        ) -> PyResult<super::mmff_properties::MmffProperties> {
            self._mmff_properties_params(Some(params))
        }
        #[pyo3(signature = (recalculate=true))]
        /// Return per-atom degree, valence and hydrogen information in atom-index order. recalculate explicitly requests rebuilding the required information.
        fn atom_metadata(&self, recalculate: bool) -> PyResult<Vec<AtomMetadata>> {
            self.inner
                .atom_metadata(recalculate)
                .map(|rows| {
                    rows.into_iter()
                        .map(|inner| AtomMetadata { inner })
                        .collect()
                })
                .map_err(ValenceError::from_rust)
        }
        /// Check whether the selected atom violates the canonical valence rules.
        fn has_valence_violation(&self, atom_id: usize) -> PyResult<bool> {
            self.inner
                .has_valence_violation(ck::AtomId::new(atom_id))
                .map_err(ValenceError::from_rust)
        }

        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, use_bond_order=crate::arguments::Argument::Missing, use_atom_weights=crate::arguments::Argument::Missing))]
        /// Return the topological atom-distance matrix, with rows and columns in atom-index order.
        fn distance_matrix(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, DistanceMatrixParams>>>,
            use_bond_order: crate::arguments::Argument<Option<bool>>,
            use_atom_weights: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<DenseMatrix> {
            let keyword_form = use_bond_order.is_supplied() || use_atom_weights.is_supplied();
            crate::arguments::configuration_arguments!(params; use_bond_order, use_atom_weights);
            let params = params.as_deref();
            if keyword_form {
                self._distance_matrix_keywords(use_bond_order, use_atom_weights)
            } else {
                self._distance_matrix_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.distance_matrix"}"#]
        #[pyo3(signature = (params=None))]
        fn _distance_matrix_params(
            &self,
            params: Option<&DistanceMatrixParams>,
        ) -> PyResult<DenseMatrix> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults;
            let params = if let Some(params) = params.as_ref() {
                params
            } else {
                defaults = ck::DistanceMatrixParams::default();
                &defaults
            };
            self.inner
                .distance_matrix_with_params(params)
                .map(|inner| DenseMatrix { inner })
                .map_err(MatrixError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.distance_matrix"}"#]
        #[pyo3(signature = (*, use_bond_order=None, use_atom_weights=None))]
        fn _distance_matrix_keywords(
            &self,
            use_bond_order: Option<bool>,
            use_atom_weights: Option<bool>,
        ) -> PyResult<DenseMatrix> {
            let defaults = ck::DistanceMatrixParams::default();
            let params = DistanceMatrixParams {
                inner: crate::parameters::NativeField::Owned(ck::DistanceMatrixParams {
                    use_bond_order: use_bond_order.unwrap_or(defaults.use_bond_order),
                    use_atom_weights: use_atom_weights.unwrap_or(defaults.use_atom_weights),
                }),
            };
            self._distance_matrix_params(Some(&params))
        }
        /// Return the topological atom-distance matrix, with rows and columns in atom-index order. Uses the explicitly supplied parameter form.
        fn distance_matrix_with_params(
            &self,
            params: &DistanceMatrixParams,
        ) -> PyResult<DenseMatrix> {
            self._distance_matrix_params(Some(params))
        }

        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, conformer_id=crate::arguments::Argument::Missing, use_atom_weights=crate::arguments::Argument::Missing))]
        /// Return the coordinate-distance matrix for the selected 3D conformer, in atom-index order.
        fn distance_matrix_3d(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, DistanceMatrix3dParams>>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            use_atom_weights: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<DenseMatrix> {
            let keyword_form = conformer_id.is_supplied() || use_atom_weights.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, use_atom_weights);
            let params = params.as_deref();
            if keyword_form {
                self._distance_matrix_3d_keywords(conformer_id, use_atom_weights)
            } else {
                self._distance_matrix_3d_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.distance_matrix_3d"}"#]
        #[pyo3(signature = (params=None))]
        fn _distance_matrix_3d_params(
            &self,
            params: Option<&DistanceMatrix3dParams>,
        ) -> PyResult<DenseMatrix> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults;
            let params = if let Some(params) = params.as_ref() {
                params
            } else {
                defaults = ck::DistanceMatrix3dParams::default();
                &defaults
            };
            self.inner
                .distance_matrix_3d_with_params(params)
                .map(|inner| DenseMatrix { inner })
                .map_err(MatrixError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.distance_matrix_3d"}"#]
        #[pyo3(signature = (*, conformer_id=None, use_atom_weights=None))]
        fn _distance_matrix_3d_keywords(
            &self,
            conformer_id: Option<usize>,
            use_atom_weights: Option<bool>,
        ) -> PyResult<DenseMatrix> {
            let defaults = ck::DistanceMatrix3dParams::default();
            let params = DistanceMatrix3dParams {
                inner: crate::parameters::NativeField::Owned(ck::DistanceMatrix3dParams {
                    conformer_id,
                    use_atom_weights: use_atom_weights.unwrap_or(defaults.use_atom_weights),
                }),
            };
            self._distance_matrix_3d_params(Some(&params))
        }
        /// Return the coordinate-distance matrix for the selected 3D conformer, in atom-index order. Uses the explicitly supplied parameter form.
        fn distance_matrix_3d_with_params(
            &self,
            params: &DistanceMatrix3dParams,
        ) -> PyResult<DenseMatrix> {
            self._distance_matrix_3d_params(Some(params))
        }
        /// Read a molecule-level property as a typed PropertyValue, or None when the key is absent.
        fn property(&self, key: &str) -> Option<super::property_values::PropertyValue> {
            self.inner
                .property(key)
                .cloned()
                .map(|inner| super::property_values::PropertyValue { inner })
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_tautomers=crate::arguments::Argument::Missing,max_transforms=crate::arguments::Argument::Missing,remove_sp3_stereo=crate::arguments::Argument::Missing,remove_bond_stereo=crate::arguments::Argument::Missing,remove_isotopic_hydrogens=crate::arguments::Argument::Missing,reassign_stereo=crate::arguments::Argument::Missing,callback=crate::arguments::Argument::Missing,scorer=crate::arguments::Argument::Missing,score_params=crate::arguments::Argument::Missing))]
        /// Enumerate tautomers using the supplied transformation, limit and callback policy; return the native result and enumeration status.
        fn enumerate_tautomers(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<Py<crate::tautomer_configuration::TautomerParams>>,
            >,
            max_tautomers: crate::arguments::Argument<Option<u32>>,
            max_transforms: crate::arguments::Argument<Option<u32>>,
            remove_sp3_stereo: crate::arguments::Argument<Option<bool>>,
            remove_bond_stereo: crate::arguments::Argument<Option<bool>>,
            remove_isotopic_hydrogens: crate::arguments::Argument<Option<bool>>,
            reassign_stereo: crate::arguments::Argument<Option<bool>>,
            callback: crate::arguments::Argument<
                Option<crate::tautomer_callbacks::EnumerationCallback>,
            >,
            scorer: crate::arguments::Argument<Option<crate::tautomer_callbacks::Scorer>>,
            score_params: crate::arguments::Argument<
                Option<Py<crate::tautomer_values::TautomerScoreParams>>,
            >,
        ) -> PyResult<crate::tautomer_values::TautomerEnumeration> {
            let keyword_form = max_tautomers.is_supplied()
                || max_transforms.is_supplied()
                || remove_sp3_stereo.is_supplied()
                || remove_bond_stereo.is_supplied()
                || remove_isotopic_hydrogens.is_supplied()
                || reassign_stereo.is_supplied()
                || callback.is_supplied()
                || scorer.is_supplied()
                || score_params.is_supplied();
            crate::arguments::configuration_arguments!(params; max_tautomers, max_transforms, remove_sp3_stereo, remove_bond_stereo, remove_isotopic_hydrogens, reassign_stereo, callback, scorer, score_params);
            if keyword_form {
                self._enumerate_tautomers_keywords(
                    py,
                    max_tautomers,
                    max_transforms,
                    remove_sp3_stereo,
                    remove_bond_stereo,
                    remove_isotopic_hydrogens,
                    reassign_stereo,
                    callback,
                    scorer,
                    score_params,
                )
            } else {
                self._enumerate_tautomers_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.enumerate_tautomers"}"#]
        #[pyo3(signature=(params=None))]
        fn _enumerate_tautomers_params(
            &self,
            py: Python<'_>,
            params: Option<Py<crate::tautomer_configuration::TautomerParams>>,
        ) -> PyResult<crate::tautomer_values::TautomerEnumeration> {
            let params = params
                .as_ref()
                .map(|value| value.try_borrow(py)?.capture(py))
                .transpose()?;
            crate::tautomer_producers::enumerate(py, self, params)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.enumerate_tautomers"}"#]
        #[pyo3(signature=(*,max_tautomers=None,max_transforms=None,remove_sp3_stereo=None,remove_bond_stereo=None,remove_isotopic_hydrogens=None,reassign_stereo=None,callback=None,scorer=None,score_params=None))]
        fn _enumerate_tautomers_keywords(
            &self,
            py: Python<'_>,
            max_tautomers: Option<u32>,
            max_transforms: Option<u32>,
            remove_sp3_stereo: Option<bool>,
            remove_bond_stereo: Option<bool>,
            remove_isotopic_hydrogens: Option<bool>,
            reassign_stereo: Option<bool>,
            callback: Option<crate::tautomer_callbacks::EnumerationCallback>,
            scorer: Option<crate::tautomer_callbacks::Scorer>,
            score_params: Option<Py<crate::tautomer_values::TautomerScoreParams>>,
        ) -> PyResult<crate::tautomer_values::TautomerEnumeration> {
            if max_tautomers.is_none()
                && max_transforms.is_none()
                && remove_sp3_stereo.is_none()
                && remove_bond_stereo.is_none()
                && remove_isotopic_hydrogens.is_none()
                && reassign_stereo.is_none()
                && callback.is_none()
                && scorer.is_none()
                && score_params.is_none()
            {
                return self._enumerate_tautomers_params(py, None);
            }
            let params = crate::tautomer_producers::keywords(
                py,
                max_tautomers,
                max_transforms,
                remove_sp3_stereo,
                remove_bond_stereo,
                remove_isotopic_hydrogens,
                reassign_stereo,
                callback,
                scorer,
                score_params,
            )?;
            crate::tautomer_producers::enumerate(py, self, Some(params))
        }
        /// Enumerate tautomers using the supplied transformation, limit and callback policy; return the native result and enumeration status. Uses the explicitly supplied parameter form.
        fn enumerate_tautomers_with_params(
            &self,
            py: Python<'_>,
            params: Py<crate::tautomer_configuration::TautomerParams>,
        ) -> PyResult<crate::tautomer_values::TautomerEnumeration> {
            self._enumerate_tautomers_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_tautomers=crate::arguments::Argument::Missing,max_transforms=crate::arguments::Argument::Missing,remove_sp3_stereo=crate::arguments::Argument::Missing,remove_bond_stereo=crate::arguments::Argument::Missing,remove_isotopic_hydrogens=crate::arguments::Argument::Missing,reassign_stereo=crate::arguments::Argument::Missing,callback=crate::arguments::Argument::Missing,scorer=crate::arguments::Argument::Missing,score_params=crate::arguments::Argument::Missing))]
        /// Select the canonical tautomer using the supplied enumeration and scoring policy.
        fn canonical_tautomer(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<Py<crate::tautomer_configuration::TautomerParams>>,
            >,
            max_tautomers: crate::arguments::Argument<Option<u32>>,
            max_transforms: crate::arguments::Argument<Option<u32>>,
            remove_sp3_stereo: crate::arguments::Argument<Option<bool>>,
            remove_bond_stereo: crate::arguments::Argument<Option<bool>>,
            remove_isotopic_hydrogens: crate::arguments::Argument<Option<bool>>,
            reassign_stereo: crate::arguments::Argument<Option<bool>>,
            callback: crate::arguments::Argument<
                Option<crate::tautomer_callbacks::EnumerationCallback>,
            >,
            scorer: crate::arguments::Argument<Option<crate::tautomer_callbacks::Scorer>>,
            score_params: crate::arguments::Argument<
                Option<Py<crate::tautomer_values::TautomerScoreParams>>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = max_tautomers.is_supplied()
                || max_transforms.is_supplied()
                || remove_sp3_stereo.is_supplied()
                || remove_bond_stereo.is_supplied()
                || remove_isotopic_hydrogens.is_supplied()
                || reassign_stereo.is_supplied()
                || callback.is_supplied()
                || scorer.is_supplied()
                || score_params.is_supplied();
            crate::arguments::configuration_arguments!(params; max_tautomers, max_transforms, remove_sp3_stereo, remove_bond_stereo, remove_isotopic_hydrogens, reassign_stereo, callback, scorer, score_params);
            if keyword_form {
                self._canonical_tautomer_keywords(
                    py,
                    max_tautomers,
                    max_transforms,
                    remove_sp3_stereo,
                    remove_bond_stereo,
                    remove_isotopic_hydrogens,
                    reassign_stereo,
                    callback,
                    scorer,
                    score_params,
                )
            } else {
                self._canonical_tautomer_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.canonical_tautomer"}"#]
        #[pyo3(signature=(params=None))]
        fn _canonical_tautomer_params(
            &self,
            py: Python<'_>,
            params: Option<Py<crate::tautomer_configuration::TautomerParams>>,
        ) -> PyResult<Molecule> {
            let params = params
                .as_ref()
                .map(|value| value.try_borrow(py)?.capture(py))
                .transpose()?;
            crate::tautomer_producers::canonical(py, self, params)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.canonical_tautomer"}"#]
        #[pyo3(signature=(*,max_tautomers=None,max_transforms=None,remove_sp3_stereo=None,remove_bond_stereo=None,remove_isotopic_hydrogens=None,reassign_stereo=None,callback=None,scorer=None,score_params=None))]
        fn _canonical_tautomer_keywords(
            &self,
            py: Python<'_>,
            max_tautomers: Option<u32>,
            max_transforms: Option<u32>,
            remove_sp3_stereo: Option<bool>,
            remove_bond_stereo: Option<bool>,
            remove_isotopic_hydrogens: Option<bool>,
            reassign_stereo: Option<bool>,
            callback: Option<crate::tautomer_callbacks::EnumerationCallback>,
            scorer: Option<crate::tautomer_callbacks::Scorer>,
            score_params: Option<Py<crate::tautomer_values::TautomerScoreParams>>,
        ) -> PyResult<Molecule> {
            if max_tautomers.is_none()
                && max_transforms.is_none()
                && remove_sp3_stereo.is_none()
                && remove_bond_stereo.is_none()
                && remove_isotopic_hydrogens.is_none()
                && reassign_stereo.is_none()
                && callback.is_none()
                && scorer.is_none()
                && score_params.is_none()
            {
                return self._canonical_tautomer_params(py, None);
            }
            let params = crate::tautomer_producers::keywords(
                py,
                max_tautomers,
                max_transforms,
                remove_sp3_stereo,
                remove_bond_stereo,
                remove_isotopic_hydrogens,
                reassign_stereo,
                callback,
                scorer,
                score_params,
            )?;
            crate::tautomer_producers::canonical(py, self, Some(params))
        }
        /// Select the canonical tautomer using the supplied enumeration and scoring policy. Uses the explicitly supplied parameter form.
        fn canonical_tautomer_with_params(
            &self,
            py: Python<'_>,
            params: Py<crate::tautomer_configuration::TautomerParams>,
        ) -> PyResult<Molecule> {
            self._canonical_tautomer_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed a 3D conformer in a new molecular value using the selected distance-geometry parameters.
        fn with_3d_conformer(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._with_3d_conformer_keywords(
                    py,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._with_3d_conformer_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformer"}"#]
        #[pyo3(signature=(params=None))]
        fn _with_3d_conformer_params(
            &self,
            py: Python<'_>,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<Molecule> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._with_3d_conformer_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformer"}"#]
        #[pyo3(signature=(*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _with_3d_conformer_keywords(
            &self,
            py: Python<'_>,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<Molecule> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._with_3d_conformer_params(py, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._with_3d_conformer_params_native(py, Some(&params))
        }
        /// Embed a 3D conformer in a new molecular value using the selected distance-geometry parameters. Uses the explicit configuration argument.
        fn with_3d_conformer_with_params(
            &self,
            py: Python<'_>,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<Molecule> {
            self._with_3d_conformer_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed a 3D conformer in place using the selected distance-geometry parameters.
        fn embed_3d_conformer_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<()> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._embed_3d_conformer_keywords(
                    py,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._embed_3d_conformer_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformer_"}"#]
        #[pyo3(signature=(params=None))]
        fn _embed_3d_conformer_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<()> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._embed_3d_conformer_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformer_"}"#]
        #[pyo3(signature=(*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _embed_3d_conformer_keywords(
            &mut self,
            py: Python<'_>,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<()> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._embed_3d_conformer_params(py, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._embed_3d_conformer_params_native(py, Some(&params))
        }
        /// Embed a 3D conformer in place using the selected distance-geometry parameters. Uses the explicit configuration argument.
        fn embed_3d_conformer_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<()> {
            self._embed_3d_conformer_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed a 3D conformer and return a new molecule with its generation outcome and diagnostics. A failed attempt is not a valid conformer identifier.
        fn with_3d_conformer_result(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._with_3d_conformer_result_keywords(
                    py,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._with_3d_conformer_result_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformer_result"}"#]
        #[pyo3(signature=(params=None))]
        fn _with_3d_conformer_result_params(
            &self,
            py: Python<'_>,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._with_3d_conformer_result_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformer_result"}"#]
        #[pyo3(signature=(*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _with_3d_conformer_result_keywords(
            &self,
            py: Python<'_>,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._with_3d_conformer_result_params(py, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._with_3d_conformer_result_params_native(py, Some(&params))
        }
        /// Embed a 3D conformer and return a new molecule with its generation outcome and diagnostics. A failed attempt is not a valid conformer identifier. Uses the explicit configuration argument.
        fn with_3d_conformer_result_with_params(
            &self,
            py: Python<'_>,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            self._with_3d_conformer_result_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed a 3D conformer in place and return the native generation outcome and diagnostics.
        fn embed_3d_conformer_result_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._embed_3d_conformer_result_keywords(
                    py,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._embed_3d_conformer_result_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformer_result_"}"#]
        #[pyo3(signature=(params=None))]
        fn _embed_3d_conformer_result_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._embed_3d_conformer_result_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformer_result_"}"#]
        #[pyo3(signature=(*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _embed_3d_conformer_result_keywords(
            &mut self,
            py: Python<'_>,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._embed_3d_conformer_result_params(py, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._embed_3d_conformer_result_params_native(py, Some(&params))
        }
        /// Embed a 3D conformer in place and return the native generation outcome and diagnostics. Uses the explicit configuration argument.
        fn embed_3d_conformer_result_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<super::embedding_values::EmbedMoleculeResult> {
            self._embed_3d_conformer_result_params(py, Some(params))
        }

        #[pyo3(signature=(num_confs,params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed the requested number of 3D conformers in a new molecular value under the selected sampling and pruning policy.
        fn with_3d_conformers(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._with_3d_conformers_keywords(
                    py,
                    num_confs,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._with_3d_conformers_params(py, num_confs, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformers"}"#]
        #[pyo3(signature=(num_confs,params=None))]
        fn _with_3d_conformers_params(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<Molecule> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._with_3d_conformers_params_native(py, num_confs, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformers"}"#]
        #[pyo3(signature=(num_confs,*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _with_3d_conformers_keywords(
            &self,
            py: Python<'_>,
            num_confs: u32,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<Molecule> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._with_3d_conformers_params(py, num_confs, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._with_3d_conformers_params_native(py, num_confs, Some(&params))
        }
        /// Embed the requested number of 3D conformers in a new molecular value under the selected sampling and pruning policy. Uses the explicit configuration argument.
        fn with_3d_conformers_with_params(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<Molecule> {
            self._with_3d_conformers_params(py, num_confs, Some(params))
        }

        #[pyo3(signature=(num_confs,params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed the requested number of 3D conformers in place under the selected sampling and pruning policy.
        fn embed_3d_conformers_(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<()> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._embed_3d_conformers_keywords(
                    py,
                    num_confs,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._embed_3d_conformers_params(py, num_confs, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformers_"}"#]
        #[pyo3(signature=(num_confs,params=None))]
        fn _embed_3d_conformers_params(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<()> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._embed_3d_conformers_params_native(py, num_confs, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformers_"}"#]
        #[pyo3(signature=(num_confs,*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _embed_3d_conformers_keywords(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<()> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._embed_3d_conformers_params(py, num_confs, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._embed_3d_conformers_params_native(py, num_confs, Some(&params))
        }
        /// Embed the requested number of 3D conformers in place under the selected sampling and pruning policy. Uses the explicit configuration argument.
        fn embed_3d_conformers_with_params_(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<()> {
            self._embed_3d_conformers_params(py, num_confs, Some(params))
        }

        #[pyo3(signature=(num_confs,params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Return a new molecule with multiple-conformer generation results and diagnostics. Requested attempts and successful conformers are distinct counts.
        fn with_3d_conformers_result(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._with_3d_conformers_result_keywords(
                    py,
                    num_confs,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._with_3d_conformers_result_params(py, num_confs, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformers_result"}"#]
        #[pyo3(signature=(num_confs,params=None))]
        fn _with_3d_conformers_result_params(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._with_3d_conformers_result_params_native(py, num_confs, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_3d_conformers_result"}"#]
        #[pyo3(signature=(num_confs,*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _with_3d_conformers_result_keywords(
            &self,
            py: Python<'_>,
            num_confs: u32,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._with_3d_conformers_result_params(py, num_confs, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._with_3d_conformers_result_params_native(py, num_confs, Some(&params))
        }
        /// Return a new molecule with multiple-conformer generation results and diagnostics. Requested attempts and successful conformers are distinct counts. Uses the explicit configuration argument.
        fn with_3d_conformers_result_with_params(
            &self,
            py: Python<'_>,
            num_confs: u32,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            self._with_3d_conformers_result_params(py, num_confs, Some(params))
        }

        #[pyo3(signature=(num_confs,params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,num_threads=crate::arguments::Argument::Missing,random_seed=crate::arguments::Argument::Missing,clear_conformers=crate::arguments::Argument::Missing,use_random_coords=crate::arguments::Argument::Missing,box_size_mult=crate::arguments::Argument::Missing,rand_neg_eig=crate::arguments::Argument::Missing,num_zero_fail=crate::arguments::Argument::Missing,coord_map=crate::arguments::Argument::Missing,optimizer_force_tol=crate::arguments::Argument::Missing,ignore_smoothing_failures=crate::arguments::Argument::Missing,enforce_chirality=crate::arguments::Argument::Missing,use_exp_torsion_angle_prefs=crate::arguments::Argument::Missing,use_basic_knowledge=crate::arguments::Argument::Missing,verbose=crate::arguments::Argument::Missing,basin_thresh=crate::arguments::Argument::Missing,prune_rms_thresh=crate::arguments::Argument::Missing,only_heavy_atoms_for_rms=crate::arguments::Argument::Missing,et_version=crate::arguments::Argument::Missing,embed_fragments_separately=crate::arguments::Argument::Missing,use_small_ring_torsions=crate::arguments::Argument::Missing,use_macrocycle_torsions=crate::arguments::Argument::Missing,use_macrocycle14config=crate::arguments::Argument::Missing,timeout=crate::arguments::Argument::Missing,cpci=crate::arguments::Argument::Missing,force_trans_amides=crate::arguments::Argument::Missing,use_symmetry_for_pruning=crate::arguments::Argument::Missing,bounds_mat_force_scaling=crate::arguments::Argument::Missing,track_failures=crate::arguments::Argument::Missing,enable_sequential_random_seeds=crate::arguments::Argument::Missing,symmetrize_conjugated_terminal_groups_for_pruning=crate::arguments::Argument::Missing))]
        /// Embed multiple conformers in place and return native generation results and diagnostics.
        fn embed_3d_conformers_result_(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::embedding_values::EmbedParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<u32>>,
            num_threads: crate::arguments::Argument<Option<i32>>,
            random_seed: crate::arguments::Argument<Option<i32>>,
            clear_conformers: crate::arguments::Argument<Option<bool>>,
            use_random_coords: crate::arguments::Argument<Option<bool>>,
            box_size_mult: crate::arguments::Argument<Option<f64>>,
            rand_neg_eig: crate::arguments::Argument<Option<bool>>,
            num_zero_fail: crate::arguments::Argument<Option<u32>>,
            coord_map: crate::arguments::Argument<
                Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            >,
            optimizer_force_tol: crate::arguments::Argument<Option<f64>>,
            ignore_smoothing_failures: crate::arguments::Argument<Option<bool>>,
            enforce_chirality: crate::arguments::Argument<Option<bool>>,
            use_exp_torsion_angle_prefs: crate::arguments::Argument<Option<bool>>,
            use_basic_knowledge: crate::arguments::Argument<Option<bool>>,
            verbose: crate::arguments::Argument<Option<bool>>,
            basin_thresh: crate::arguments::Argument<Option<f64>>,
            prune_rms_thresh: crate::arguments::Argument<Option<f64>>,
            only_heavy_atoms_for_rms: crate::arguments::Argument<Option<bool>>,
            et_version: crate::arguments::Argument<Option<u32>>,
            embed_fragments_separately: crate::arguments::Argument<Option<bool>>,
            use_small_ring_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle_torsions: crate::arguments::Argument<Option<bool>>,
            use_macrocycle14config: crate::arguments::Argument<Option<bool>>,
            timeout: crate::arguments::Argument<Option<u32>>,
            cpci: crate::arguments::Argument<Option<std::collections::BTreeMap<(u32, u32), f64>>>,
            force_trans_amides: crate::arguments::Argument<Option<bool>>,
            use_symmetry_for_pruning: crate::arguments::Argument<Option<bool>>,
            bounds_mat_force_scaling: crate::arguments::Argument<Option<f64>>,
            track_failures: crate::arguments::Argument<Option<bool>>,
            enable_sequential_random_seeds: crate::arguments::Argument<Option<bool>>,
            symmetrize_conjugated_terminal_groups_for_pruning: crate::arguments::Argument<
                Option<bool>,
            >,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            let keyword_form = max_iterations.is_supplied()
                || num_threads.is_supplied()
                || random_seed.is_supplied()
                || clear_conformers.is_supplied()
                || use_random_coords.is_supplied()
                || box_size_mult.is_supplied()
                || rand_neg_eig.is_supplied()
                || num_zero_fail.is_supplied()
                || coord_map.is_supplied()
                || optimizer_force_tol.is_supplied()
                || ignore_smoothing_failures.is_supplied()
                || enforce_chirality.is_supplied()
                || use_exp_torsion_angle_prefs.is_supplied()
                || use_basic_knowledge.is_supplied()
                || verbose.is_supplied()
                || basin_thresh.is_supplied()
                || prune_rms_thresh.is_supplied()
                || only_heavy_atoms_for_rms.is_supplied()
                || et_version.is_supplied()
                || embed_fragments_separately.is_supplied()
                || use_small_ring_torsions.is_supplied()
                || use_macrocycle_torsions.is_supplied()
                || use_macrocycle14config.is_supplied()
                || timeout.is_supplied()
                || cpci.is_supplied()
                || force_trans_amides.is_supplied()
                || use_symmetry_for_pruning.is_supplied()
                || bounds_mat_force_scaling.is_supplied()
                || track_failures.is_supplied()
                || enable_sequential_random_seeds.is_supplied()
                || symmetrize_conjugated_terminal_groups_for_pruning.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, num_threads, random_seed, clear_conformers, use_random_coords, box_size_mult, rand_neg_eig, num_zero_fail, coord_map, optimizer_force_tol, ignore_smoothing_failures, enforce_chirality, use_exp_torsion_angle_prefs, use_basic_knowledge, verbose, basin_thresh, prune_rms_thresh, only_heavy_atoms_for_rms, et_version, embed_fragments_separately, use_small_ring_torsions, use_macrocycle_torsions, use_macrocycle14config, timeout, cpci, force_trans_amides, use_symmetry_for_pruning, bounds_mat_force_scaling, track_failures, enable_sequential_random_seeds, symmetrize_conjugated_terminal_groups_for_pruning);
            let params = params.as_deref();
            if keyword_form {
                self._embed_3d_conformers_result_keywords(
                    py,
                    num_confs,
                    max_iterations,
                    num_threads,
                    random_seed,
                    clear_conformers,
                    use_random_coords,
                    box_size_mult,
                    rand_neg_eig,
                    num_zero_fail,
                    coord_map,
                    optimizer_force_tol,
                    ignore_smoothing_failures,
                    enforce_chirality,
                    use_exp_torsion_angle_prefs,
                    use_basic_knowledge,
                    verbose,
                    basin_thresh,
                    prune_rms_thresh,
                    only_heavy_atoms_for_rms,
                    et_version,
                    embed_fragments_separately,
                    use_small_ring_torsions,
                    use_macrocycle_torsions,
                    use_macrocycle14config,
                    timeout,
                    cpci,
                    force_trans_amides,
                    use_symmetry_for_pruning,
                    bounds_mat_force_scaling,
                    track_failures,
                    enable_sequential_random_seeds,
                    symmetrize_conjugated_terminal_groups_for_pruning,
                )
            } else {
                self._embed_3d_conformers_result_params(py, num_confs, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformers_result_"}"#]
        #[pyo3(signature=(num_confs,params=None))]
        fn _embed_3d_conformers_result_params(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: Option<&super::embedding_values::EmbedParams>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._embed_3d_conformers_result_params_native(py, num_confs, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.embed_3d_conformers_result_"}"#]
        #[pyo3(signature=(num_confs,*,max_iterations=None,num_threads=None,random_seed=None,clear_conformers=None,use_random_coords=None,box_size_mult=None,rand_neg_eig=None,num_zero_fail=None,coord_map=None,optimizer_force_tol=None,ignore_smoothing_failures=None,enforce_chirality=None,use_exp_torsion_angle_prefs=None,use_basic_knowledge=None,verbose=None,basin_thresh=None,prune_rms_thresh=None,only_heavy_atoms_for_rms=None,et_version=None,embed_fragments_separately=None,use_small_ring_torsions=None,use_macrocycle_torsions=None,use_macrocycle14config=None,timeout=None,cpci=None,force_trans_amides=None,use_symmetry_for_pruning=None,bounds_mat_force_scaling=None,track_failures=None,enable_sequential_random_seeds=None,symmetrize_conjugated_terminal_groups_for_pruning=None))]
        fn _embed_3d_conformers_result_keywords(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            max_iterations: Option<u32>,
            num_threads: Option<i32>,
            random_seed: Option<i32>,
            clear_conformers: Option<bool>,
            use_random_coords: Option<bool>,
            box_size_mult: Option<f64>,
            rand_neg_eig: Option<bool>,
            num_zero_fail: Option<u32>,
            coord_map: Option<std::collections::BTreeMap<i32, [f64; 3]>>,
            optimizer_force_tol: Option<f64>,
            ignore_smoothing_failures: Option<bool>,
            enforce_chirality: Option<bool>,
            use_exp_torsion_angle_prefs: Option<bool>,
            use_basic_knowledge: Option<bool>,
            verbose: Option<bool>,
            basin_thresh: Option<f64>,
            prune_rms_thresh: Option<f64>,
            only_heavy_atoms_for_rms: Option<bool>,
            et_version: Option<u32>,
            embed_fragments_separately: Option<bool>,
            use_small_ring_torsions: Option<bool>,
            use_macrocycle_torsions: Option<bool>,
            use_macrocycle14config: Option<bool>,
            timeout: Option<u32>,
            cpci: Option<std::collections::BTreeMap<(u32, u32), f64>>,
            force_trans_amides: Option<bool>,
            use_symmetry_for_pruning: Option<bool>,
            bounds_mat_force_scaling: Option<f64>,
            track_failures: Option<bool>,
            enable_sequential_random_seeds: Option<bool>,
            symmetrize_conjugated_terminal_groups_for_pruning: Option<bool>,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            if max_iterations.is_none()
                && num_threads.is_none()
                && random_seed.is_none()
                && clear_conformers.is_none()
                && use_random_coords.is_none()
                && box_size_mult.is_none()
                && rand_neg_eig.is_none()
                && num_zero_fail.is_none()
                && coord_map.is_none()
                && optimizer_force_tol.is_none()
                && ignore_smoothing_failures.is_none()
                && enforce_chirality.is_none()
                && use_exp_torsion_angle_prefs.is_none()
                && use_basic_knowledge.is_none()
                && verbose.is_none()
                && basin_thresh.is_none()
                && prune_rms_thresh.is_none()
                && only_heavy_atoms_for_rms.is_none()
                && et_version.is_none()
                && embed_fragments_separately.is_none()
                && use_small_ring_torsions.is_none()
                && use_macrocycle_torsions.is_none()
                && use_macrocycle14config.is_none()
                && timeout.is_none()
                && cpci.is_none()
                && force_trans_amides.is_none()
                && use_symmetry_for_pruning.is_none()
                && bounds_mat_force_scaling.is_none()
                && track_failures.is_none()
                && enable_sequential_random_seeds.is_none()
                && symmetrize_conjugated_terminal_groups_for_pruning.is_none()
            {
                return self._embed_3d_conformers_result_params(py, num_confs, None);
            }
            let params = super::embedding_calls::keywords(
                max_iterations,
                num_threads,
                random_seed,
                clear_conformers,
                use_random_coords,
                box_size_mult,
                rand_neg_eig,
                num_zero_fail,
                coord_map,
                optimizer_force_tol,
                ignore_smoothing_failures,
                enforce_chirality,
                use_exp_torsion_angle_prefs,
                use_basic_knowledge,
                verbose,
                basin_thresh,
                prune_rms_thresh,
                only_heavy_atoms_for_rms,
                et_version,
                embed_fragments_separately,
                use_small_ring_torsions,
                use_macrocycle_torsions,
                use_macrocycle14config,
                timeout,
                cpci,
                force_trans_amides,
                use_symmetry_for_pruning,
                bounds_mat_force_scaling,
                track_failures,
                enable_sequential_random_seeds,
                symmetrize_conjugated_terminal_groups_for_pruning,
            );
            self._embed_3d_conformers_result_params_native(py, num_confs, Some(&params))
        }
        /// Embed multiple conformers in place and return native generation results and diagnostics. Uses the explicit configuration argument.
        fn embed_3d_conformers_result_with_params_(
            &mut self,
            py: Python<'_>,
            num_confs: u32,
            params: &super::embedding_values::EmbedParams,
        ) -> PyResult<super::embedding_values::EmbedMultipleConfsResult> {
            self._embed_3d_conformers_result_params(py, num_confs, Some(params))
        }

        #[pyo3(signature=(options=crate::arguments::Argument::Missing,*,try_embedding=crate::arguments::Argument::Missing,only_unassigned=crate::arguments::Argument::Missing,max_isomers=crate::arguments::Argument::Missing,random_source=crate::arguments::Argument::Missing,unique=crate::arguments::Argument::Missing,only_stereo_groups=crate::arguments::Argument::Missing))]
        /// Enumerate stereoisomer molecules under the configured assignment, uniqueness and sampling policy.
        fn enumerate_stereoisomers(
            &self,
            py: Python<'_>,
            options: crate::arguments::Argument<
                Option<Py<super::stereoisomer_options::StereoisomerOptions>>,
            >,
            try_embedding: crate::arguments::Argument<Option<bool>>,
            only_unassigned: crate::arguments::Argument<Option<bool>>,
            max_isomers: crate::arguments::Argument<Option<usize>>,
            random_source: crate::arguments::Argument<
                Option<Py<super::stereoisomer_options::StereoisomerRandomSource>>,
            >,
            unique: crate::arguments::Argument<Option<bool>>,
            only_stereo_groups: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::stereoisomer_iterator::StereoisomerIterator> {
            let keyword_form = try_embedding.is_supplied()
                || only_unassigned.is_supplied()
                || max_isomers.is_supplied()
                || random_source.is_supplied()
                || unique.is_supplied()
                || only_stereo_groups.is_supplied();
            crate::arguments::configuration_arguments!(options; try_embedding, only_unassigned, max_isomers, random_source, unique, only_stereo_groups);
            if keyword_form {
                self._enumerate_stereoisomers_keywords(
                    py,
                    try_embedding,
                    only_unassigned,
                    max_isomers,
                    random_source,
                    unique,
                    only_stereo_groups,
                )
            } else {
                self._enumerate_stereoisomers_options(py, options)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.enumerate_stereoisomers"}"#]
        #[pyo3(signature=(options=None))]
        fn _enumerate_stereoisomers_options(
            &self,
            py: Python<'_>,
            options: Option<Py<super::stereoisomer_options::StereoisomerOptions>>,
        ) -> PyResult<super::stereoisomer_iterator::StereoisomerIterator> {
            super::stereoisomer_iterator::enumerate(py, self, options)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.enumerate_stereoisomers"}"#]
        #[pyo3(signature=(*,try_embedding=None,only_unassigned=None,max_isomers=None,random_source=None,unique=None,only_stereo_groups=None))]
        fn _enumerate_stereoisomers_keywords(
            &self,
            py: Python<'_>,
            try_embedding: Option<bool>,
            only_unassigned: Option<bool>,
            max_isomers: Option<usize>,
            random_source: Option<Py<super::stereoisomer_options::StereoisomerRandomSource>>,
            unique: Option<bool>,
            only_stereo_groups: Option<bool>,
        ) -> PyResult<super::stereoisomer_iterator::StereoisomerIterator> {
            if try_embedding.is_none()
                && only_unassigned.is_none()
                && max_isomers.is_none()
                && random_source.is_none()
                && unique.is_none()
                && only_stereo_groups.is_none()
            {
                return self._enumerate_stereoisomers_options(py, None);
            }
            let defaults = ck::StereoisomerOptions::default();
            let options = super::stereoisomer_options::StereoisomerOptions::new(
                py,
                try_embedding.unwrap_or(defaults.try_embedding),
                only_unassigned.unwrap_or(defaults.only_unassigned),
                only_stereo_groups.unwrap_or(defaults.only_stereo_groups),
                max_isomers.unwrap_or(defaults.max_isomers),
                unique.unwrap_or(defaults.unique),
                random_source,
            );
            self._enumerate_stereoisomers_options(py, Some(Py::new(py, options)?))
        }
        /// Enumerate stereoisomer molecules under the configured assignment, uniqueness and sampling policy. Uses the explicitly supplied options.
        fn enumerate_stereoisomers_with_options(
            &self,
            py: Python<'_>,
            options: Py<super::stereoisomer_options::StereoisomerOptions>,
        ) -> PyResult<super::stereoisomer_iterator::StereoisomerIterator> {
            self._enumerate_stereoisomers_options(py, Some(options))
        }
        /// Enumerate stereoisomers with the explicitly supplied random-bit source, preserving its native state and sampling behavior.
        fn enumerate_stereoisomers_with_random_bits(
            &self,
            py: Python<'_>,
            options: Py<super::stereoisomer_options::StereoisomerOptions>,
            callback: super::stereoisomer_options::RandomBitsCallback,
        ) -> PyResult<super::stereoisomer_iterator::StereoisomerIterator> {
            super::stereoisomer_iterator::with_random_bits(py, self, options, callback)
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,operations=crate::arguments::Argument::Missing))]
        /// Return a sanitized molecular value using the requested sanitization stages. The source is unchanged; errors retain their native type and cause.
        fn sanitize(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::sanitize_values::SanitizeParams>>,
            >,
            operations: crate::arguments::Argument<
                Option<super::sanitize_values::SanitizeOperationsInput>,
            >,
        ) -> PyResult<Molecule> {
            let keyword_form = operations.is_supplied();
            crate::arguments::configuration_arguments!(params; operations);
            let params = params.as_deref();
            if keyword_form {
                self._chem_sanitize_keywords(py, operations)
            } else {
                self._chem_sanitize_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.sanitize"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_sanitize_params(
            &self,
            py: Python<'_>,
            params: Option<&super::sanitize_values::SanitizeParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.sanitize_with_params(&value)),
                None => py.detach(|| self.inner.sanitize()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.sanitize"}"#]
        #[pyo3(signature=(*,operations=None))]
        fn _chem_sanitize_keywords(
            &self,
            py: Python<'_>,
            operations: Option<super::sanitize_values::SanitizeOperationsInput>,
        ) -> PyResult<Molecule> {
            if operations.is_none() {
                return self._chem_sanitize_params(py, None);
            }
            let params = super::chemistry_calls::sanitize(operations)?;
            self._chem_sanitize_params(py, Some(&params))
        }
        /// Return a sanitized molecular value using the requested sanitization stages. The source is unchanged; errors retain their native type and cause. Uses the explicit configuration argument.
        fn sanitize_with_params(
            &self,
            py: Python<'_>,
            params: &super::sanitize_values::SanitizeParams,
        ) -> PyResult<Molecule> {
            self._chem_sanitize_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,operations=crate::arguments::Argument::Missing))]
        /// Run the requested sanitization stages in place. Errors retain their native type and cause.
        fn sanitize_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::sanitize_values::SanitizeParams>>,
            >,
            operations: crate::arguments::Argument<
                Option<super::sanitize_values::SanitizeOperationsInput>,
            >,
        ) -> PyResult<()> {
            let keyword_form = operations.is_supplied();
            crate::arguments::configuration_arguments!(params; operations);
            let params = params.as_deref();
            if keyword_form {
                self._chem_sanitize_inplace_keywords(py, operations)
            } else {
                self._chem_sanitize_inplace_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.sanitize_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_sanitize_inplace_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::sanitize_values::SanitizeParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.sanitize_with_params_(&value)),
                None => py.detach(|| self.inner.sanitize_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.sanitize_"}"#]
        #[pyo3(signature=(*,operations=None))]
        fn _chem_sanitize_inplace_keywords(
            &mut self,
            py: Python<'_>,
            operations: Option<super::sanitize_values::SanitizeOperationsInput>,
        ) -> PyResult<()> {
            if operations.is_none() {
                return self._chem_sanitize_inplace_params(py, None);
            }
            let params = super::chemistry_calls::sanitize(operations)?;
            self._chem_sanitize_inplace_params(py, Some(&params))
        }
        /// Run the requested sanitization stages in place. Errors retain their native type and cause. Uses the explicit configuration argument.
        fn sanitize_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::sanitize_values::SanitizeParams,
        ) -> PyResult<()> {
            self._chem_sanitize_inplace_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,mark_atoms_bonds=crate::arguments::Argument::Missing,canonical=crate::arguments::Argument::Missing,max_backtracks=crate::arguments::Argument::Missing))]
        /// Return a new molecule with Kekule bond orders assigned under the selected aromatic-flag policy.
        fn with_kekulized_bonds(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::KekulizeParams>>>,
            mark_atoms_bonds: crate::arguments::Argument<Option<bool>>,
            canonical: crate::arguments::Argument<Option<bool>>,
            max_backtracks: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<Molecule> {
            let keyword_form = mark_atoms_bonds.is_supplied()
                || canonical.is_supplied()
                || max_backtracks.is_supplied();
            crate::arguments::configuration_arguments!(params; mark_atoms_bonds, canonical, max_backtracks);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_kekulized_bonds_keywords(
                    py,
                    mark_atoms_bonds,
                    canonical,
                    max_backtracks,
                )
            } else {
                self._chem_with_kekulized_bonds_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_kekulized_bonds"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_kekulized_bonds_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::KekulizeParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.with_kekulized_bonds_with_params(&value)),
                None => py.detach(|| self.inner.with_kekulized_bonds()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_kekulized_bonds"}"#]
        #[pyo3(signature=(*,mark_atoms_bonds=None,canonical=None,max_backtracks=None))]
        fn _chem_with_kekulized_bonds_keywords(
            &self,
            py: Python<'_>,
            mark_atoms_bonds: Option<bool>,
            canonical: Option<bool>,
            max_backtracks: Option<u32>,
        ) -> PyResult<Molecule> {
            if mark_atoms_bonds.is_none() && canonical.is_none() && max_backtracks.is_none() {
                return self._chem_with_kekulized_bonds_params(py, None);
            }
            let params =
                super::chemistry_calls::kekulize(mark_atoms_bonds, canonical, max_backtracks)?;
            self._chem_with_kekulized_bonds_params(py, Some(&params))
        }
        /// Return a new molecule with Kekule bond orders assigned under the selected aromatic-flag policy. Uses the explicit configuration argument.
        fn with_kekulized_bonds_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::KekulizeParams,
        ) -> PyResult<Molecule> {
            self._chem_with_kekulized_bonds_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,mark_atoms_bonds=crate::arguments::Argument::Missing,canonical=crate::arguments::Argument::Missing,max_backtracks=crate::arguments::Argument::Missing))]
        /// Assign Kekule bond orders in place under the selected aromatic-flag policy.
        fn kekulize_bonds_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::KekulizeParams>>>,
            mark_atoms_bonds: crate::arguments::Argument<Option<bool>>,
            canonical: crate::arguments::Argument<Option<bool>>,
            max_backtracks: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<()> {
            let keyword_form = mark_atoms_bonds.is_supplied()
                || canonical.is_supplied()
                || max_backtracks.is_supplied();
            crate::arguments::configuration_arguments!(params; mark_atoms_bonds, canonical, max_backtracks);
            let params = params.as_deref();
            if keyword_form {
                self._chem_kekulize_bonds_keywords(py, mark_atoms_bonds, canonical, max_backtracks)
            } else {
                self._chem_kekulize_bonds_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.kekulize_bonds_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_kekulize_bonds_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::params::KekulizeParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.kekulize_bonds_with_params_(&value)),
                None => py.detach(|| self.inner.kekulize_bonds_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.kekulize_bonds_"}"#]
        #[pyo3(signature=(*,mark_atoms_bonds=None,canonical=None,max_backtracks=None))]
        fn _chem_kekulize_bonds_keywords(
            &mut self,
            py: Python<'_>,
            mark_atoms_bonds: Option<bool>,
            canonical: Option<bool>,
            max_backtracks: Option<u32>,
        ) -> PyResult<()> {
            if mark_atoms_bonds.is_none() && canonical.is_none() && max_backtracks.is_none() {
                return self._chem_kekulize_bonds_params(py, None);
            }
            let params =
                super::chemistry_calls::kekulize(mark_atoms_bonds, canonical, max_backtracks)?;
            self._chem_kekulize_bonds_params(py, Some(&params))
        }
        /// Assign Kekule bond orders in place under the selected aromatic-flag policy. Uses the explicit configuration argument.
        fn kekulize_bonds_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::params::KekulizeParams,
        ) -> PyResult<()> {
            self._chem_kekulize_bonds_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,model=crate::arguments::Argument::Missing))]
        /// Return a new molecule with aromaticity assigned by the selected aromaticity model.
        fn with_assigned_aromaticity(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::selectors::AromaticityParams>>,
            >,
            model: crate::arguments::Argument<Option<super::selectors::AromaticityModelInput>>,
        ) -> PyResult<Molecule> {
            let keyword_form = model.is_supplied();
            crate::arguments::configuration_arguments!(params; model);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_assigned_aromaticity_keywords(py, model)
            } else {
                self._chem_with_assigned_aromaticity_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_aromaticity"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_assigned_aromaticity_params(
            &self,
            py: Python<'_>,
            params: Option<&super::selectors::AromaticityParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => {
                    py.detach(|| self.inner.with_assigned_aromaticity_with_params(&value))
                }
                None => py.detach(|| self.inner.with_assigned_aromaticity()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_aromaticity"}"#]
        #[pyo3(signature=(*,model=None))]
        fn _chem_with_assigned_aromaticity_keywords(
            &self,
            py: Python<'_>,
            model: Option<super::selectors::AromaticityModelInput>,
        ) -> PyResult<Molecule> {
            if model.is_none() {
                return self._chem_with_assigned_aromaticity_params(py, None);
            }
            let params = super::chemistry_calls::aromaticity(model)?;
            self._chem_with_assigned_aromaticity_params(py, Some(&params))
        }
        /// Return a new molecule with aromaticity assigned by the selected aromaticity model. Uses the explicit configuration argument.
        fn with_assigned_aromaticity_with_params(
            &self,
            py: Python<'_>,
            params: &super::selectors::AromaticityParams,
        ) -> PyResult<Molecule> {
            self._chem_with_assigned_aromaticity_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,model=crate::arguments::Argument::Missing))]
        /// Assign aromaticity in place using the selected aromaticity model.
        fn assign_aromaticity_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::selectors::AromaticityParams>>,
            >,
            model: crate::arguments::Argument<Option<super::selectors::AromaticityModelInput>>,
        ) -> PyResult<()> {
            let keyword_form = model.is_supplied();
            crate::arguments::configuration_arguments!(params; model);
            let params = params.as_deref();
            if keyword_form {
                self._chem_assign_aromaticity_keywords(py, model)
            } else {
                self._chem_assign_aromaticity_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_aromaticity_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_assign_aromaticity_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::selectors::AromaticityParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.assign_aromaticity_with_params_(&value)),
                None => py.detach(|| self.inner.assign_aromaticity_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_aromaticity_"}"#]
        #[pyo3(signature=(*,model=None))]
        fn _chem_assign_aromaticity_keywords(
            &mut self,
            py: Python<'_>,
            model: Option<super::selectors::AromaticityModelInput>,
        ) -> PyResult<()> {
            if model.is_none() {
                return self._chem_assign_aromaticity_params(py, None);
            }
            let params = super::chemistry_calls::aromaticity(model)?;
            self._chem_assign_aromaticity_params(py, Some(&params))
        }
        /// Assign aromaticity in place using the selected aromaticity model. Uses the explicit configuration argument.
        fn assign_aromaticity_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::selectors::AromaticityParams,
        ) -> PyResult<()> {
            self._chem_assign_aromaticity_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,model=crate::arguments::Argument::Missing,strict=crate::arguments::Argument::Missing))]
        /// Return a new molecule with valence information assigned under the selected valence model and strictness.
        fn with_assigned_valence(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::selectors::ValenceParams>>>,
            model: crate::arguments::Argument<Option<super::selectors::ValenceModelInput>>,
            strict: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = model.is_supplied() || strict.is_supplied();
            crate::arguments::configuration_arguments!(params; model, strict);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_assigned_valence_keywords(py, model, strict)
            } else {
                self._chem_with_assigned_valence_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_valence"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_assigned_valence_params(
            &self,
            py: Python<'_>,
            params: Option<&super::selectors::ValenceParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.with_assigned_valence_with_params(&value)),
                None => py.detach(|| self.inner.with_assigned_valence()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_valence"}"#]
        #[pyo3(signature=(*,model=None,strict=None))]
        fn _chem_with_assigned_valence_keywords(
            &self,
            py: Python<'_>,
            model: Option<super::selectors::ValenceModelInput>,
            strict: Option<bool>,
        ) -> PyResult<Molecule> {
            if model.is_none() && strict.is_none() {
                return self._chem_with_assigned_valence_params(py, None);
            }
            let params = super::chemistry_calls::valence(model, strict)?;
            self._chem_with_assigned_valence_params(py, Some(&params))
        }
        /// Return a new molecule with valence information assigned under the selected valence model and strictness. Uses the explicit configuration argument.
        fn with_assigned_valence_with_params(
            &self,
            py: Python<'_>,
            params: &super::selectors::ValenceParams,
        ) -> PyResult<Molecule> {
            self._chem_with_assigned_valence_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,model=crate::arguments::Argument::Missing,strict=crate::arguments::Argument::Missing))]
        /// Assign valence information in place under the selected valence model and strictness.
        fn assign_valence_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::selectors::ValenceParams>>>,
            model: crate::arguments::Argument<Option<super::selectors::ValenceModelInput>>,
            strict: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = model.is_supplied() || strict.is_supplied();
            crate::arguments::configuration_arguments!(params; model, strict);
            let params = params.as_deref();
            if keyword_form {
                self._chem_assign_valence_keywords(py, model, strict)
            } else {
                self._chem_assign_valence_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_valence_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_assign_valence_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::selectors::ValenceParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.assign_valence_with_params_(&value)),
                None => py.detach(|| self.inner.assign_valence_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_valence_"}"#]
        #[pyo3(signature=(*,model=None,strict=None))]
        fn _chem_assign_valence_keywords(
            &mut self,
            py: Python<'_>,
            model: Option<super::selectors::ValenceModelInput>,
            strict: Option<bool>,
        ) -> PyResult<()> {
            if model.is_none() && strict.is_none() {
                return self._chem_assign_valence_params(py, None);
            }
            let params = super::chemistry_calls::valence(model, strict)?;
            self._chem_assign_valence_params(py, Some(&params))
        }
        /// Assign valence information in place under the selected valence model and strictness. Uses the explicit configuration argument.
        fn assign_valence_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::selectors::ValenceParams,
        ) -> PyResult<()> {
            self._chem_assign_valence_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,include_dative_bonds=crate::arguments::Argument::Missing,include_hydrogen_bonds=crate::arguments::Argument::Missing))]
        /// Return a new molecule with ring families calculated using the selected ring-search policy.
        fn with_assigned_ring_families(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::RingSearchParams>>>,
            include_dative_bonds: crate::arguments::Argument<Option<bool>>,
            include_hydrogen_bonds: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form =
                include_dative_bonds.is_supplied() || include_hydrogen_bonds.is_supplied();
            crate::arguments::configuration_arguments!(params; include_dative_bonds, include_hydrogen_bonds);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_assigned_ring_families_keywords(
                    py,
                    include_dative_bonds,
                    include_hydrogen_bonds,
                )
            } else {
                self._chem_with_assigned_ring_families_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_ring_families"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_assigned_ring_families_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::RingSearchParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => {
                    py.detach(|| self.inner.with_assigned_ring_families_with_params(&value))
                }
                None => py.detach(|| self.inner.with_assigned_ring_families()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_assigned_ring_families"}"#]
        #[pyo3(signature=(*,include_dative_bonds=None,include_hydrogen_bonds=None))]
        fn _chem_with_assigned_ring_families_keywords(
            &self,
            py: Python<'_>,
            include_dative_bonds: Option<bool>,
            include_hydrogen_bonds: Option<bool>,
        ) -> PyResult<Molecule> {
            if include_dative_bonds.is_none() && include_hydrogen_bonds.is_none() {
                return self._chem_with_assigned_ring_families_params(py, None);
            }
            let params =
                super::chemistry_calls::rings(include_dative_bonds, include_hydrogen_bonds)?;
            self._chem_with_assigned_ring_families_params(py, Some(&params))
        }
        /// Return a new molecule with ring families calculated using the selected ring-search policy. Uses the explicit configuration argument.
        fn with_assigned_ring_families_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::RingSearchParams,
        ) -> PyResult<Molecule> {
            self._chem_with_assigned_ring_families_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,include_dative_bonds=crate::arguments::Argument::Missing,include_hydrogen_bonds=crate::arguments::Argument::Missing))]
        /// Calculate and install ring families in place using the selected ring-search policy.
        fn assign_ring_families_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::RingSearchParams>>>,
            include_dative_bonds: crate::arguments::Argument<Option<bool>>,
            include_hydrogen_bonds: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form =
                include_dative_bonds.is_supplied() || include_hydrogen_bonds.is_supplied();
            crate::arguments::configuration_arguments!(params; include_dative_bonds, include_hydrogen_bonds);
            let params = params.as_deref();
            if keyword_form {
                self._chem_assign_ring_families_keywords(
                    py,
                    include_dative_bonds,
                    include_hydrogen_bonds,
                )
            } else {
                self._chem_assign_ring_families_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_ring_families_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_assign_ring_families_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::params::RingSearchParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.assign_ring_families_with_params_(&value)),
                None => py.detach(|| self.inner.assign_ring_families_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_ring_families_"}"#]
        #[pyo3(signature=(*,include_dative_bonds=None,include_hydrogen_bonds=None))]
        fn _chem_assign_ring_families_keywords(
            &mut self,
            py: Python<'_>,
            include_dative_bonds: Option<bool>,
            include_hydrogen_bonds: Option<bool>,
        ) -> PyResult<()> {
            if include_dative_bonds.is_none() && include_hydrogen_bonds.is_none() {
                return self._chem_assign_ring_families_params(py, None);
            }
            let params =
                super::chemistry_calls::rings(include_dative_bonds, include_hydrogen_bonds)?;
            self._chem_assign_ring_families_params(py, Some(&params))
        }
        /// Calculate and install ring families in place using the selected ring-search policy. Uses the explicit configuration argument.
        fn assign_ring_families_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::params::RingSearchParams,
        ) -> PyResult<()> {
            self._chem_assign_ring_families_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,conformer_id=crate::arguments::Argument::Missing,replace_existing_tags=crate::arguments::Argument::Missing))]
        /// Return a new molecule with chiral tags assigned from the selected coordinate structure.
        fn with_chiral_tags_from_structure(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::StructureTagParams>>,
            >,
            conformer_id: crate::arguments::Argument<Option<i32>>,
            replace_existing_tags: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = conformer_id.is_supplied() || replace_existing_tags.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, replace_existing_tags);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_chiral_tags_from_structure_keywords(
                    py,
                    conformer_id,
                    replace_existing_tags,
                )
            } else {
                self._chem_with_chiral_tags_from_structure_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_chiral_tags_from_structure"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_chiral_tags_from_structure_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::StructureTagParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| {
                    self.inner
                        .with_chiral_tags_from_structure_with_params(&value)
                }),
                None => py.detach(|| self.inner.with_chiral_tags_from_structure()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_chiral_tags_from_structure"}"#]
        #[pyo3(signature=(*,conformer_id=None,replace_existing_tags=None))]
        fn _chem_with_chiral_tags_from_structure_keywords(
            &self,
            py: Python<'_>,
            conformer_id: Option<i32>,
            replace_existing_tags: Option<bool>,
        ) -> PyResult<Molecule> {
            if conformer_id.is_none() && replace_existing_tags.is_none() {
                return self._chem_with_chiral_tags_from_structure_params(py, None);
            }
            let params = super::chemistry_calls::structure(conformer_id, replace_existing_tags)?;
            self._chem_with_chiral_tags_from_structure_params(py, Some(&params))
        }
        /// Return a new molecule with chiral tags assigned from the selected coordinate structure. Uses the explicit configuration argument.
        fn with_chiral_tags_from_structure_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::StructureTagParams,
        ) -> PyResult<Molecule> {
            self._chem_with_chiral_tags_from_structure_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,conformer_id=crate::arguments::Argument::Missing,replace_existing_tags=crate::arguments::Argument::Missing))]
        /// Assign chiral tags in place from the selected coordinate structure.
        fn assign_chiral_tags_from_structure_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::StructureTagParams>>,
            >,
            conformer_id: crate::arguments::Argument<Option<i32>>,
            replace_existing_tags: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = conformer_id.is_supplied() || replace_existing_tags.is_supplied();
            crate::arguments::configuration_arguments!(params; conformer_id, replace_existing_tags);
            let params = params.as_deref();
            if keyword_form {
                self._chem_assign_chiral_tags_from_structure_keywords(
                    py,
                    conformer_id,
                    replace_existing_tags,
                )
            } else {
                self._chem_assign_chiral_tags_from_structure_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_chiral_tags_from_structure_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_assign_chiral_tags_from_structure_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::params::StructureTagParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| {
                    self.inner
                        .assign_chiral_tags_from_structure_with_params_(&value)
                }),
                None => py.detach(|| self.inner.assign_chiral_tags_from_structure_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_chiral_tags_from_structure_"}"#]
        #[pyo3(signature=(*,conformer_id=None,replace_existing_tags=None))]
        fn _chem_assign_chiral_tags_from_structure_keywords(
            &mut self,
            py: Python<'_>,
            conformer_id: Option<i32>,
            replace_existing_tags: Option<bool>,
        ) -> PyResult<()> {
            if conformer_id.is_none() && replace_existing_tags.is_none() {
                return self._chem_assign_chiral_tags_from_structure_params(py, None);
            }
            let params = super::chemistry_calls::structure(conformer_id, replace_existing_tags)?;
            self._chem_assign_chiral_tags_from_structure_params(py, Some(&params))
        }
        /// Assign chiral tags in place from the selected coordinate structure. Uses the explicit configuration argument.
        fn assign_chiral_tags_from_structure_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::params::StructureTagParams,
        ) -> PyResult<()> {
            self._chem_assign_chiral_tags_from_structure_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,clean=crate::arguments::Argument::Missing,flag_possible=crate::arguments::Argument::Missing,allow_nontetrahedral=crate::arguments::Argument::Missing))]
        /// Return the native potential-stereo information for atoms and bonds under the selected options.
        fn potential_stereo(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::PotentialStereoParams>>,
            >,
            clean: crate::arguments::Argument<Option<bool>>,
            flag_possible: crate::arguments::Argument<Option<bool>>,
            allow_nontetrahedral: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::stereo_values::PotentialStereoResult> {
            let keyword_form = clean.is_supplied()
                || flag_possible.is_supplied()
                || allow_nontetrahedral.is_supplied();
            crate::arguments::configuration_arguments!(params; clean, flag_possible, allow_nontetrahedral);
            let params = params.as_deref();
            if keyword_form {
                self._chem_potential_stereo_keywords(py, clean, flag_possible, allow_nontetrahedral)
            } else {
                self._chem_potential_stereo_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.potential_stereo"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_potential_stereo_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::PotentialStereoParams>,
        ) -> PyResult<super::stereo_values::PotentialStereoResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.potential_stereo_with_params(&value)),
                None => py.detach(|| self.inner.potential_stereo()),
            }
            .map(|inner| super::stereo_values::PotentialStereoResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.potential_stereo"}"#]
        #[pyo3(signature=(*,clean=None,flag_possible=None,allow_nontetrahedral=None))]
        fn _chem_potential_stereo_keywords(
            &self,
            py: Python<'_>,
            clean: Option<bool>,
            flag_possible: Option<bool>,
            allow_nontetrahedral: Option<bool>,
        ) -> PyResult<super::stereo_values::PotentialStereoResult> {
            if clean.is_none() && flag_possible.is_none() && allow_nontetrahedral.is_none() {
                return self._chem_potential_stereo_params(py, None);
            }
            let params =
                super::chemistry_calls::potential(clean, flag_possible, allow_nontetrahedral)?;
            self._chem_potential_stereo_params(py, Some(&params))
        }
        /// Return the native potential-stereo information for atoms and bonds under the selected options. Uses the explicitly supplied parameter form.
        fn potential_stereo_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::PotentialStereoParams,
        ) -> PyResult<super::stereo_values::PotentialStereoResult> {
            self._chem_potential_stereo_params(py, Some(params))
        }

        #[pyo3(signature=(options=crate::arguments::Argument::Missing,*,atoms=crate::arguments::Argument::Missing,bonds=crate::arguments::Argument::Missing,max_recursive_iterations=crate::arguments::Argument::Missing))]
        /// Return a new molecule with CIP labels assigned under the selected complete or partial assignment options.
        fn with_cip_labels(
            &self,
            py: Python<'_>,
            options: crate::arguments::Argument<
                Option<PyRef<'_, super::numerical_params::CipLabelOptions>>,
            >,
            atoms: crate::arguments::Argument<Option<Vec<usize>>>,
            bonds: crate::arguments::Argument<Option<Vec<usize>>>,
            max_recursive_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<Molecule> {
            let keyword_form = atoms.is_supplied()
                || bonds.is_supplied()
                || max_recursive_iterations.is_supplied();
            crate::arguments::configuration_arguments!(options; atoms, bonds, max_recursive_iterations);
            let options = options.as_deref();
            if keyword_form {
                self._chem_with_cip_labels_keywords(py, atoms, bonds, max_recursive_iterations)
            } else {
                self._chem_with_cip_labels_params(py, options)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_cip_labels"}"#]
        #[pyo3(signature=(options=None))]
        fn _chem_with_cip_labels_params(
            &self,
            py: Python<'_>,
            options: Option<&super::numerical_params::CipLabelOptions>,
        ) -> PyResult<Molecule> {
            let options = options.map(|value| value.to_rust(py)).transpose()?;
            self._chem_with_cip_labels_params_native(py, options.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_cip_labels"}"#]
        #[pyo3(signature=(*,atoms=None,bonds=None,max_recursive_iterations=None))]
        fn _chem_with_cip_labels_keywords(
            &self,
            py: Python<'_>,
            atoms: Option<Vec<usize>>,
            bonds: Option<Vec<usize>>,
            max_recursive_iterations: Option<u32>,
        ) -> PyResult<Molecule> {
            if atoms.is_none() && bonds.is_none() && max_recursive_iterations.is_none() {
                return self._chem_with_cip_labels_params(py, None);
            }
            let params = super::chemistry_calls::cip(atoms, bonds, max_recursive_iterations)?;
            self._chem_with_cip_labels_params_native(py, Some(&params))
        }
        /// Return a new molecule with CIP labels assigned under the selected complete or partial assignment options. Uses the explicit configuration argument.
        fn with_cip_labels_with_options(
            &self,
            py: Python<'_>,
            options: &super::numerical_params::CipLabelOptions,
        ) -> PyResult<Molecule> {
            self._chem_with_cip_labels_params(py, Some(options))
        }

        #[pyo3(signature=(options=crate::arguments::Argument::Missing,*,atoms=crate::arguments::Argument::Missing,bonds=crate::arguments::Argument::Missing,max_recursive_iterations=crate::arguments::Argument::Missing))]
        /// Assign CIP labels in place under the selected complete or partial assignment options.
        fn assign_cip_labels_(
            &mut self,
            py: Python<'_>,
            options: crate::arguments::Argument<
                Option<PyRef<'_, super::numerical_params::CipLabelOptions>>,
            >,
            atoms: crate::arguments::Argument<Option<Vec<usize>>>,
            bonds: crate::arguments::Argument<Option<Vec<usize>>>,
            max_recursive_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<()> {
            let keyword_form = atoms.is_supplied()
                || bonds.is_supplied()
                || max_recursive_iterations.is_supplied();
            crate::arguments::configuration_arguments!(options; atoms, bonds, max_recursive_iterations);
            let options = options.as_deref();
            if keyword_form {
                self._chem_assign_cip_labels_keywords(py, atoms, bonds, max_recursive_iterations)
            } else {
                self._chem_assign_cip_labels_params(py, options)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_cip_labels_"}"#]
        #[pyo3(signature=(options=None))]
        fn _chem_assign_cip_labels_params(
            &mut self,
            py: Python<'_>,
            options: Option<&super::numerical_params::CipLabelOptions>,
        ) -> PyResult<()> {
            let options = options.map(|value| value.to_rust(py)).transpose()?;
            self._chem_assign_cip_labels_params_native(py, options.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.assign_cip_labels_"}"#]
        #[pyo3(signature=(*,atoms=None,bonds=None,max_recursive_iterations=None))]
        fn _chem_assign_cip_labels_keywords(
            &mut self,
            py: Python<'_>,
            atoms: Option<Vec<usize>>,
            bonds: Option<Vec<usize>>,
            max_recursive_iterations: Option<u32>,
        ) -> PyResult<()> {
            if atoms.is_none() && bonds.is_none() && max_recursive_iterations.is_none() {
                return self._chem_assign_cip_labels_params(py, None);
            }
            let params = super::chemistry_calls::cip(atoms, bonds, max_recursive_iterations)?;
            self._chem_assign_cip_labels_params_native(py, Some(&params))
        }
        /// Assign CIP labels in place under the selected complete or partial assignment options. Uses the explicit configuration argument.
        fn assign_cip_labels_with_options_(
            &mut self,
            py: Python<'_>,
            options: &super::numerical_params::CipLabelOptions,
        ) -> PyResult<()> {
            self._chem_assign_cip_labels_params(py, Some(options))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,explicit_only=crate::arguments::Argument::Missing,add_coords=crate::arguments::Argument::Missing,add_residue_info=crate::arguments::Argument::Missing,skip_queries=crate::arguments::Argument::Missing,only_on_atoms=crate::arguments::Argument::Missing))]
        /// Return a molecule with explicit hydrogen atoms added; the source molecule is unchanged.
        /// Supply either AddHsParams or keyword options, never both. With no
        /// options, use the Rust defaults. only_on_atoms limits addition to the
        /// selected zero-based atom indices; add_coords requests hydrogen coordinates.
        fn with_hydrogens(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::chemical_params::AddHsParams>>,
            >,
            explicit_only: crate::arguments::Argument<Option<bool>>,
            add_coords: crate::arguments::Argument<Option<bool>>,
            add_residue_info: crate::arguments::Argument<Option<bool>>,
            skip_queries: crate::arguments::Argument<Option<bool>>,
            only_on_atoms: crate::arguments::Argument<Option<Vec<usize>>>,
        ) -> PyResult<Molecule> {
            let keyword_form = explicit_only.is_supplied()
                || add_coords.is_supplied()
                || add_residue_info.is_supplied()
                || skip_queries.is_supplied()
                || only_on_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; explicit_only, add_coords, add_residue_info, skip_queries, only_on_atoms);
            let params = params.as_deref();
            if keyword_form {
                self._chem_with_hydrogens_keywords(
                    py,
                    explicit_only,
                    add_coords,
                    add_residue_info,
                    skip_queries,
                    only_on_atoms,
                )
            } else {
                self._chem_with_hydrogens_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_hydrogens"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_with_hydrogens_params(
            &self,
            py: Python<'_>,
            params: Option<&super::chemical_params::AddHsParams>,
        ) -> PyResult<Molecule> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._chem_with_hydrogens_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_hydrogens"}"#]
        #[pyo3(signature=(*,explicit_only=None,add_coords=None,add_residue_info=None,skip_queries=None,only_on_atoms=None))]
        fn _chem_with_hydrogens_keywords(
            &self,
            py: Python<'_>,
            explicit_only: Option<bool>,
            add_coords: Option<bool>,
            add_residue_info: Option<bool>,
            skip_queries: Option<bool>,
            only_on_atoms: Option<Vec<usize>>,
        ) -> PyResult<Molecule> {
            if explicit_only.is_none()
                && add_coords.is_none()
                && add_residue_info.is_none()
                && skip_queries.is_none()
                && only_on_atoms.is_none()
            {
                return self._chem_with_hydrogens_params(py, None);
            }
            let params = super::chemistry_calls::hydrogens(
                explicit_only,
                add_coords,
                add_residue_info,
                skip_queries,
                only_on_atoms,
            )?;
            self._chem_with_hydrogens_params_native(py, Some(&params))
        }
        /// Return a new molecule with explicit hydrogen atoms added using AddHsParams; the source molecule is unchanged.
        fn with_hydrogens_with_params(
            &self,
            py: Python<'_>,
            params: &super::chemical_params::AddHsParams,
        ) -> PyResult<Molecule> {
            self._chem_with_hydrogens_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,explicit_only=crate::arguments::Argument::Missing,add_coords=crate::arguments::Argument::Missing,add_residue_info=crate::arguments::Argument::Missing,skip_queries=crate::arguments::Argument::Missing,only_on_atoms=crate::arguments::Argument::Missing))]
        /// In place, add explicit hydrogen atoms using AddHsParams or keyword options.
        /// Copy-on-write isolates molecular values sharing the original storage.
        /// Returns None; use with_hydrogens to obtain a new molecule instead.
        fn add_hydrogens_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::chemical_params::AddHsParams>>,
            >,
            explicit_only: crate::arguments::Argument<Option<bool>>,
            add_coords: crate::arguments::Argument<Option<bool>>,
            add_residue_info: crate::arguments::Argument<Option<bool>>,
            skip_queries: crate::arguments::Argument<Option<bool>>,
            only_on_atoms: crate::arguments::Argument<Option<Vec<usize>>>,
        ) -> PyResult<()> {
            let keyword_form = explicit_only.is_supplied()
                || add_coords.is_supplied()
                || add_residue_info.is_supplied()
                || skip_queries.is_supplied()
                || only_on_atoms.is_supplied();
            crate::arguments::configuration_arguments!(params; explicit_only, add_coords, add_residue_info, skip_queries, only_on_atoms);
            let params = params.as_deref();
            if keyword_form {
                self._chem_add_hydrogens_keywords(
                    py,
                    explicit_only,
                    add_coords,
                    add_residue_info,
                    skip_queries,
                    only_on_atoms,
                )
            } else {
                self._chem_add_hydrogens_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_hydrogens_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_add_hydrogens_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::chemical_params::AddHsParams>,
        ) -> PyResult<()> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._chem_add_hydrogens_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.add_hydrogens_"}"#]
        #[pyo3(signature=(*,explicit_only=None,add_coords=None,add_residue_info=None,skip_queries=None,only_on_atoms=None))]
        fn _chem_add_hydrogens_keywords(
            &mut self,
            py: Python<'_>,
            explicit_only: Option<bool>,
            add_coords: Option<bool>,
            add_residue_info: Option<bool>,
            skip_queries: Option<bool>,
            only_on_atoms: Option<Vec<usize>>,
        ) -> PyResult<()> {
            if explicit_only.is_none()
                && add_coords.is_none()
                && add_residue_info.is_none()
                && skip_queries.is_none()
                && only_on_atoms.is_none()
            {
                return self._chem_add_hydrogens_params(py, None);
            }
            let params = super::chemistry_calls::hydrogens(
                explicit_only,
                add_coords,
                add_residue_info,
                skip_queries,
                only_on_atoms,
            )?;
            self._chem_add_hydrogens_params_native(py, Some(&params))
        }
        /// Add hydrogen atoms in place using AddHsParams and copy-on-write isolation from other molecule values.
        fn add_hydrogens_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::chemical_params::AddHsParams,
        ) -> PyResult<()> {
            self._chem_add_hydrogens_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,remove_degree_zero=crate::arguments::Argument::Missing,remove_higher_degrees=crate::arguments::Argument::Missing,remove_only_h_neighbors=crate::arguments::Argument::Missing,remove_isotopes=crate::arguments::Argument::Missing,remove_and_track_isotopes=crate::arguments::Argument::Missing,remove_dummy_neighbors=crate::arguments::Argument::Missing,remove_defining_bond_stereo=crate::arguments::Argument::Missing,remove_with_wedged_bond=crate::arguments::Argument::Missing,remove_with_query=crate::arguments::Argument::Missing,remove_mapped=crate::arguments::Argument::Missing,remove_in_sgroups=crate::arguments::Argument::Missing,show_warnings=crate::arguments::Argument::Missing,remove_nonimplicit=crate::arguments::Argument::Missing,update_explicit_count=crate::arguments::Argument::Missing,remove_hydrides=crate::arguments::Argument::Missing,remove_nontetrahedral_neighbors=crate::arguments::Argument::Missing,sanitize=crate::arguments::Argument::Missing))]
        /// Return a new molecule with eligible hydrogen atoms removed under RemoveHsParams. The removal policy determines which isotopic, mapped and stereochemistry-defining hydrogens remain.
        fn without_hydrogens(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::RemoveHsParams>>>,
            remove_degree_zero: crate::arguments::Argument<Option<bool>>,
            remove_higher_degrees: crate::arguments::Argument<Option<bool>>,
            remove_only_h_neighbors: crate::arguments::Argument<Option<bool>>,
            remove_isotopes: crate::arguments::Argument<Option<bool>>,
            remove_and_track_isotopes: crate::arguments::Argument<Option<bool>>,
            remove_dummy_neighbors: crate::arguments::Argument<Option<bool>>,
            remove_defining_bond_stereo: crate::arguments::Argument<Option<bool>>,
            remove_with_wedged_bond: crate::arguments::Argument<Option<bool>>,
            remove_with_query: crate::arguments::Argument<Option<bool>>,
            remove_mapped: crate::arguments::Argument<Option<bool>>,
            remove_in_sgroups: crate::arguments::Argument<Option<bool>>,
            show_warnings: crate::arguments::Argument<Option<bool>>,
            remove_nonimplicit: crate::arguments::Argument<Option<bool>>,
            update_explicit_count: crate::arguments::Argument<Option<bool>>,
            remove_hydrides: crate::arguments::Argument<Option<bool>>,
            remove_nontetrahedral_neighbors: crate::arguments::Argument<Option<bool>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Molecule> {
            let keyword_form = remove_degree_zero.is_supplied()
                || remove_higher_degrees.is_supplied()
                || remove_only_h_neighbors.is_supplied()
                || remove_isotopes.is_supplied()
                || remove_and_track_isotopes.is_supplied()
                || remove_dummy_neighbors.is_supplied()
                || remove_defining_bond_stereo.is_supplied()
                || remove_with_wedged_bond.is_supplied()
                || remove_with_query.is_supplied()
                || remove_mapped.is_supplied()
                || remove_in_sgroups.is_supplied()
                || show_warnings.is_supplied()
                || remove_nonimplicit.is_supplied()
                || update_explicit_count.is_supplied()
                || remove_hydrides.is_supplied()
                || remove_nontetrahedral_neighbors.is_supplied()
                || sanitize.is_supplied();
            crate::arguments::configuration_arguments!(params; remove_degree_zero, remove_higher_degrees, remove_only_h_neighbors, remove_isotopes, remove_and_track_isotopes, remove_dummy_neighbors, remove_defining_bond_stereo, remove_with_wedged_bond, remove_with_query, remove_mapped, remove_in_sgroups, show_warnings, remove_nonimplicit, update_explicit_count, remove_hydrides, remove_nontetrahedral_neighbors, sanitize);
            let params = params.as_deref();
            if keyword_form {
                self._chem_without_hydrogens_keywords(
                    py,
                    remove_degree_zero,
                    remove_higher_degrees,
                    remove_only_h_neighbors,
                    remove_isotopes,
                    remove_and_track_isotopes,
                    remove_dummy_neighbors,
                    remove_defining_bond_stereo,
                    remove_with_wedged_bond,
                    remove_with_query,
                    remove_mapped,
                    remove_in_sgroups,
                    show_warnings,
                    remove_nonimplicit,
                    update_explicit_count,
                    remove_hydrides,
                    remove_nontetrahedral_neighbors,
                    sanitize,
                )
            } else {
                self._chem_without_hydrogens_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.without_hydrogens"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_without_hydrogens_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::RemoveHsParams>,
        ) -> PyResult<Molecule> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.without_hydrogens_with_params(&value)),
                None => py.detach(|| self.inner.without_hydrogens()),
            }
            .map(|inner| Molecule { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.without_hydrogens"}"#]
        #[pyo3(signature=(*,remove_degree_zero=None,remove_higher_degrees=None,remove_only_h_neighbors=None,remove_isotopes=None,remove_and_track_isotopes=None,remove_dummy_neighbors=None,remove_defining_bond_stereo=None,remove_with_wedged_bond=None,remove_with_query=None,remove_mapped=None,remove_in_sgroups=None,show_warnings=None,remove_nonimplicit=None,update_explicit_count=None,remove_hydrides=None,remove_nontetrahedral_neighbors=None,sanitize=None))]
        fn _chem_without_hydrogens_keywords(
            &self,
            py: Python<'_>,
            remove_degree_zero: Option<bool>,
            remove_higher_degrees: Option<bool>,
            remove_only_h_neighbors: Option<bool>,
            remove_isotopes: Option<bool>,
            remove_and_track_isotopes: Option<bool>,
            remove_dummy_neighbors: Option<bool>,
            remove_defining_bond_stereo: Option<bool>,
            remove_with_wedged_bond: Option<bool>,
            remove_with_query: Option<bool>,
            remove_mapped: Option<bool>,
            remove_in_sgroups: Option<bool>,
            show_warnings: Option<bool>,
            remove_nonimplicit: Option<bool>,
            update_explicit_count: Option<bool>,
            remove_hydrides: Option<bool>,
            remove_nontetrahedral_neighbors: Option<bool>,
            sanitize: Option<bool>,
        ) -> PyResult<Molecule> {
            if remove_degree_zero.is_none()
                && remove_higher_degrees.is_none()
                && remove_only_h_neighbors.is_none()
                && remove_isotopes.is_none()
                && remove_and_track_isotopes.is_none()
                && remove_dummy_neighbors.is_none()
                && remove_defining_bond_stereo.is_none()
                && remove_with_wedged_bond.is_none()
                && remove_with_query.is_none()
                && remove_mapped.is_none()
                && remove_in_sgroups.is_none()
                && show_warnings.is_none()
                && remove_nonimplicit.is_none()
                && update_explicit_count.is_none()
                && remove_hydrides.is_none()
                && remove_nontetrahedral_neighbors.is_none()
                && sanitize.is_none()
            {
                return self._chem_without_hydrogens_params(py, None);
            }
            let params = super::chemistry_calls::remove_hydrogens(
                remove_degree_zero,
                remove_higher_degrees,
                remove_only_h_neighbors,
                remove_isotopes,
                remove_and_track_isotopes,
                remove_dummy_neighbors,
                remove_defining_bond_stereo,
                remove_with_wedged_bond,
                remove_with_query,
                remove_mapped,
                remove_in_sgroups,
                show_warnings,
                remove_nonimplicit,
                update_explicit_count,
                remove_hydrides,
                remove_nontetrahedral_neighbors,
                sanitize,
            )?;
            self._chem_without_hydrogens_params(py, Some(&params))
        }
        /// Return a new molecule with eligible hydrogen atoms removed under RemoveHsParams. The removal policy determines which isotopic, mapped and stereochemistry-defining hydrogens remain. Uses the explicit configuration argument.
        fn without_hydrogens_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::RemoveHsParams,
        ) -> PyResult<Molecule> {
            self._chem_without_hydrogens_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,remove_degree_zero=crate::arguments::Argument::Missing,remove_higher_degrees=crate::arguments::Argument::Missing,remove_only_h_neighbors=crate::arguments::Argument::Missing,remove_isotopes=crate::arguments::Argument::Missing,remove_and_track_isotopes=crate::arguments::Argument::Missing,remove_dummy_neighbors=crate::arguments::Argument::Missing,remove_defining_bond_stereo=crate::arguments::Argument::Missing,remove_with_wedged_bond=crate::arguments::Argument::Missing,remove_with_query=crate::arguments::Argument::Missing,remove_mapped=crate::arguments::Argument::Missing,remove_in_sgroups=crate::arguments::Argument::Missing,show_warnings=crate::arguments::Argument::Missing,remove_nonimplicit=crate::arguments::Argument::Missing,update_explicit_count=crate::arguments::Argument::Missing,remove_hydrides=crate::arguments::Argument::Missing,remove_nontetrahedral_neighbors=crate::arguments::Argument::Missing,sanitize=crate::arguments::Argument::Missing))]
        /// Remove eligible hydrogen atoms in place under RemoveHsParams. The removal policy determines which isotopic, mapped and stereochemistry-defining hydrogens remain.
        fn remove_hydrogens_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::RemoveHsParams>>>,
            remove_degree_zero: crate::arguments::Argument<Option<bool>>,
            remove_higher_degrees: crate::arguments::Argument<Option<bool>>,
            remove_only_h_neighbors: crate::arguments::Argument<Option<bool>>,
            remove_isotopes: crate::arguments::Argument<Option<bool>>,
            remove_and_track_isotopes: crate::arguments::Argument<Option<bool>>,
            remove_dummy_neighbors: crate::arguments::Argument<Option<bool>>,
            remove_defining_bond_stereo: crate::arguments::Argument<Option<bool>>,
            remove_with_wedged_bond: crate::arguments::Argument<Option<bool>>,
            remove_with_query: crate::arguments::Argument<Option<bool>>,
            remove_mapped: crate::arguments::Argument<Option<bool>>,
            remove_in_sgroups: crate::arguments::Argument<Option<bool>>,
            show_warnings: crate::arguments::Argument<Option<bool>>,
            remove_nonimplicit: crate::arguments::Argument<Option<bool>>,
            update_explicit_count: crate::arguments::Argument<Option<bool>>,
            remove_hydrides: crate::arguments::Argument<Option<bool>>,
            remove_nontetrahedral_neighbors: crate::arguments::Argument<Option<bool>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<()> {
            let keyword_form = remove_degree_zero.is_supplied()
                || remove_higher_degrees.is_supplied()
                || remove_only_h_neighbors.is_supplied()
                || remove_isotopes.is_supplied()
                || remove_and_track_isotopes.is_supplied()
                || remove_dummy_neighbors.is_supplied()
                || remove_defining_bond_stereo.is_supplied()
                || remove_with_wedged_bond.is_supplied()
                || remove_with_query.is_supplied()
                || remove_mapped.is_supplied()
                || remove_in_sgroups.is_supplied()
                || show_warnings.is_supplied()
                || remove_nonimplicit.is_supplied()
                || update_explicit_count.is_supplied()
                || remove_hydrides.is_supplied()
                || remove_nontetrahedral_neighbors.is_supplied()
                || sanitize.is_supplied();
            crate::arguments::configuration_arguments!(params; remove_degree_zero, remove_higher_degrees, remove_only_h_neighbors, remove_isotopes, remove_and_track_isotopes, remove_dummy_neighbors, remove_defining_bond_stereo, remove_with_wedged_bond, remove_with_query, remove_mapped, remove_in_sgroups, show_warnings, remove_nonimplicit, update_explicit_count, remove_hydrides, remove_nontetrahedral_neighbors, sanitize);
            let params = params.as_deref();
            if keyword_form {
                self._chem_remove_hydrogens_keywords(
                    py,
                    remove_degree_zero,
                    remove_higher_degrees,
                    remove_only_h_neighbors,
                    remove_isotopes,
                    remove_and_track_isotopes,
                    remove_dummy_neighbors,
                    remove_defining_bond_stereo,
                    remove_with_wedged_bond,
                    remove_with_query,
                    remove_mapped,
                    remove_in_sgroups,
                    show_warnings,
                    remove_nonimplicit,
                    update_explicit_count,
                    remove_hydrides,
                    remove_nontetrahedral_neighbors,
                    sanitize,
                )
            } else {
                self._chem_remove_hydrogens_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.remove_hydrogens_"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_remove_hydrogens_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::params::RemoveHsParams>,
        ) -> PyResult<()> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.remove_hydrogens_with_params_(&value)),
                None => py.detach(|| self.inner.remove_hydrogens_()),
            }
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.remove_hydrogens_"}"#]
        #[pyo3(signature=(*,remove_degree_zero=None,remove_higher_degrees=None,remove_only_h_neighbors=None,remove_isotopes=None,remove_and_track_isotopes=None,remove_dummy_neighbors=None,remove_defining_bond_stereo=None,remove_with_wedged_bond=None,remove_with_query=None,remove_mapped=None,remove_in_sgroups=None,show_warnings=None,remove_nonimplicit=None,update_explicit_count=None,remove_hydrides=None,remove_nontetrahedral_neighbors=None,sanitize=None))]
        fn _chem_remove_hydrogens_keywords(
            &mut self,
            py: Python<'_>,
            remove_degree_zero: Option<bool>,
            remove_higher_degrees: Option<bool>,
            remove_only_h_neighbors: Option<bool>,
            remove_isotopes: Option<bool>,
            remove_and_track_isotopes: Option<bool>,
            remove_dummy_neighbors: Option<bool>,
            remove_defining_bond_stereo: Option<bool>,
            remove_with_wedged_bond: Option<bool>,
            remove_with_query: Option<bool>,
            remove_mapped: Option<bool>,
            remove_in_sgroups: Option<bool>,
            show_warnings: Option<bool>,
            remove_nonimplicit: Option<bool>,
            update_explicit_count: Option<bool>,
            remove_hydrides: Option<bool>,
            remove_nontetrahedral_neighbors: Option<bool>,
            sanitize: Option<bool>,
        ) -> PyResult<()> {
            if remove_degree_zero.is_none()
                && remove_higher_degrees.is_none()
                && remove_only_h_neighbors.is_none()
                && remove_isotopes.is_none()
                && remove_and_track_isotopes.is_none()
                && remove_dummy_neighbors.is_none()
                && remove_defining_bond_stereo.is_none()
                && remove_with_wedged_bond.is_none()
                && remove_with_query.is_none()
                && remove_mapped.is_none()
                && remove_in_sgroups.is_none()
                && show_warnings.is_none()
                && remove_nonimplicit.is_none()
                && update_explicit_count.is_none()
                && remove_hydrides.is_none()
                && remove_nontetrahedral_neighbors.is_none()
                && sanitize.is_none()
            {
                return self._chem_remove_hydrogens_params(py, None);
            }
            let params = super::chemistry_calls::remove_hydrogens(
                remove_degree_zero,
                remove_higher_degrees,
                remove_only_h_neighbors,
                remove_isotopes,
                remove_and_track_isotopes,
                remove_dummy_neighbors,
                remove_defining_bond_stereo,
                remove_with_wedged_bond,
                remove_with_query,
                remove_mapped,
                remove_in_sgroups,
                show_warnings,
                remove_nonimplicit,
                update_explicit_count,
                remove_hydrides,
                remove_nontetrahedral_neighbors,
                sanitize,
            )?;
            self._chem_remove_hydrogens_params(py, Some(&params))
        }
        /// Remove eligible hydrogen atoms in place under RemoveHsParams. The removal policy determines which isotopic, mapped and stereochemistry-defining hydrogens remain. Uses the explicit configuration argument.
        fn remove_hydrogens_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::params::RemoveHsParams,
        ) -> PyResult<()> {
            self._chem_remove_hydrogens_params(py, Some(params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,terms=crate::arguments::Argument::Missing))]
        /// Calculate the tautomer score using the supplied scoring terms and options.
        fn tautomer_score(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::tautomer_values::TautomerScoreParams>>,
            >,
            terms: crate::arguments::Argument<
                Option<Vec<Py<super::tautomer_values::TautomerScoreTerm>>>,
            >,
        ) -> PyResult<super::tautomer_values::TautomerScore> {
            let keyword_form = terms.is_supplied();
            crate::arguments::configuration_arguments!(params; terms);
            let params = params.as_deref();
            if keyword_form {
                self._chem_tautomer_score_keywords(py, terms)
            } else {
                self._chem_tautomer_score_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.tautomer_score"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_tautomer_score_params(
            &self,
            py: Python<'_>,
            params: Option<&super::tautomer_values::TautomerScoreParams>,
        ) -> PyResult<super::tautomer_values::TautomerScore> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._chem_tautomer_score_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.tautomer_score"}"#]
        #[pyo3(signature=(*,terms=None))]
        fn _chem_tautomer_score_keywords(
            &self,
            py: Python<'_>,
            terms: Option<Vec<Py<super::tautomer_values::TautomerScoreTerm>>>,
        ) -> PyResult<super::tautomer_values::TautomerScore> {
            if terms.is_none() {
                return self._chem_tautomer_score_params(py, None);
            }
            let params = super::chemistry_calls::tautomer_score(py, terms)?;
            self._chem_tautomer_score_params_native(py, Some(&params))
        }
        /// Calculate the tautomer score using the supplied scoring terms and options. Uses the explicitly supplied parameter form.
        fn tautomer_score_with_params(
            &self,
            py: Python<'_>,
            params: &super::tautomer_values::TautomerScoreParams,
        ) -> PyResult<super::tautomer_values::TautomerScore> {
            self._chem_tautomer_score_params(py, Some(params))
        }
        /// Return a new molecule with radical-electron counts assigned by the native radical assignment.
        fn with_assigned_radicals(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.with_assigned_radicals())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Assign radical-electron counts in place using the native radical assignment.
        fn assign_radicals_(&mut self, py: Python<'_>) -> PyResult<()> {
            py.detach(|| self.inner.assign_radicals_())
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,include_dative_bonds=crate::arguments::Argument::Missing,include_hydrogen_bonds=crate::arguments::Argument::Missing))]
        /// Calculate the complete symmetrized SSSR ring collection without modifying the source molecule.
        fn symmetrized_sssr(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::RingSearchParams>>>,
            include_dative_bonds: crate::arguments::Argument<Option<bool>>,
            include_hydrogen_bonds: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<RingSearchResult> {
            let keyword_form =
                include_dative_bonds.is_supplied() || include_hydrogen_bonds.is_supplied();
            crate::arguments::configuration_arguments!(params; include_dative_bonds, include_hydrogen_bonds);
            let params = params.as_deref();
            if keyword_form {
                self._chem_symmetrized_sssr_keywords(
                    py,
                    include_dative_bonds,
                    include_hydrogen_bonds,
                )
            } else {
                self._chem_symmetrized_sssr_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.symmetrized_sssr"}"#]
        #[pyo3(signature=(params=None))]
        fn _chem_symmetrized_sssr_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::RingSearchParams>,
        ) -> PyResult<RingSearchResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.symmetrized_sssr_with_params(&value)),
                None => py.detach(|| self.inner.symmetrized_sssr()),
            }
            .map(|inner| RingSearchResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.symmetrized_sssr"}"#]
        #[pyo3(signature=(*,include_dative_bonds=None,include_hydrogen_bonds=None))]
        fn _chem_symmetrized_sssr_keywords(
            &self,
            py: Python<'_>,
            include_dative_bonds: Option<bool>,
            include_hydrogen_bonds: Option<bool>,
        ) -> PyResult<RingSearchResult> {
            if include_dative_bonds.is_none() && include_hydrogen_bonds.is_none() {
                return self._chem_symmetrized_sssr_params(py, None);
            }
            let params =
                super::chemistry_calls::rings(include_dative_bonds, include_hydrogen_bonds)?;
            self._chem_symmetrized_sssr_params(py, Some(&params))
        }
        /// Calculate the complete symmetrized SSSR ring collection without modifying the source molecule. Uses the explicitly supplied parameter form.
        fn symmetrized_sssr_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::RingSearchParams,
        ) -> PyResult<RingSearchResult> {
            self._chem_symmetrized_sssr_params(py, Some(params))
        }
        /// Return a new molecule with ordinary ring information assigned.
        fn with_assigned_rings(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.with_assigned_rings())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Calculate and install ordinary ring information in place.
        fn assign_rings_(&mut self, py: Python<'_>) -> PyResult<()> {
            py.detach(|| self.inner.assign_rings_())
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return the Bemis-Murcko ring-and-linker scaffold as a new molecule.
        fn murcko_scaffold(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.murcko_scaffold())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return the native net scaffold as a new molecule, leaving the source unchanged.
        fn net_scaffold(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.net_scaffold())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return the native Murcko decomposition as a new molecular value.
        fn murcko_decompose(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.murcko_decompose())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return the fragment selected by the native largest-fragment policy; the source is unchanged.
        fn largest_fragment(&self, py: Python<'_>) -> PyResult<Molecule> {
            py.detach(|| self.inner.largest_fragment())
                .map(|inner| Molecule { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Split disconnected components into new molecular values in native component order.
        fn fragments(&self, py: Python<'_>) -> PyResult<Vec<Molecule>> {
            py.detach(|| self.inner.fragments())
                .map(|values| values.into_iter().map(|inner| Molecule { inner }).collect())
                .map_err(super::operation_errors::OperationError::from_rust)
        }

        #[staticmethod]
        /// Construct a molecule from explicit topology, coordinates and molecular properties using the canonical validation boundary.
        fn from_parts(
            topology: &super::topology_values::TopologyBlock,
            coordinates: &super::detached_values::CoordinateBlock,
            properties: &super::metadata_values::MoleculeProperties,
        ) -> PyResult<Self> {
            ck::Molecule::from_parts(
                topology.inner().clone(),
                coordinates.inner.clone(),
                properties.inner.clone(),
            )
            .map(|inner| Self { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        /// Return whether UFF parameters are available for every required atom in the molecule.
        fn uff_has_all_molecule_params(&self, py: Python<'_>) -> PyResult<bool> {
            py.detach(|| self.inner.uff_has_all_molecule_params())
                .map_err(super::forcefield_calls::UffParameterQueryError::from_rust)
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,vdw_threshold=crate::arguments::Argument::Missing,conformer_id=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Evaluate UFF energy and coordinate gradients for the selected conformer without minimizing it.
        fn uff_energy_gradient(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::UffEvaluationParams>>,
            >,
            vdw_threshold: crate::arguments::Argument<Option<f64>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::numerical_results::UffEnergyGradient> {
            let keyword_form = vdw_threshold.is_supplied()
                || conformer_id.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; vdw_threshold, conformer_id, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_uff_energy_gradient_keywords(
                    py,
                    vdw_threshold,
                    conformer_id,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_uff_energy_gradient_params(py, params)
            }
        }
        /// Evaluate UFF energy and coordinate gradients for the selected conformer without minimizing it. Uses the explicitly supplied parameter form.
        fn uff_energy_gradient_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::UffEvaluationParams,
        ) -> PyResult<super::numerical_results::UffEnergyGradient> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| self.inner.uff_energy_gradient_with_params(&params))
                .map(|inner| super::numerical_results::UffEnergyGradient { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.uff_energy_gradient"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_uff_energy_gradient_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::UffEvaluationParams>,
        ) -> PyResult<super::numerical_results::UffEnergyGradient> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.uff_energy_gradient_with_params(&value)),
                None => py.detach(|| self.inner.uff_energy_gradient()),
            }
            .map(|inner| super::numerical_results::UffEnergyGradient { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.uff_energy_gradient"}"#]
        #[pyo3(signature=(*,vdw_threshold=None,conformer_id=None,ignore_interfragment_interactions=None))]
        fn _ff_uff_energy_gradient_keywords(
            &self,
            py: Python<'_>,
            vdw_threshold: Option<f64>,
            conformer_id: Option<usize>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<super::numerical_results::UffEnergyGradient> {
            let params = super::forcefield_calls::uff_energy_gradient_params(
                vdw_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            );
            self._ff_uff_energy_gradient_params(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,mmff_variant=crate::arguments::Argument::Missing,non_bonded_threshold=crate::arguments::Argument::Missing,conformer_id=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Evaluate MMFF energy and coordinate gradients for the selected conformer and MMFF variant without minimizing it.
        fn mmff_energy_gradient(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::MmffEvaluationParams>>,
            >,
            mmff_variant: crate::arguments::Argument<Option<crate::leaf_values::MmffVariantInput>>,
            non_bonded_threshold: crate::arguments::Argument<Option<f64>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Option<super::numerical_results::MmffEnergyGradient>> {
            let keyword_form = mmff_variant.is_supplied()
                || non_bonded_threshold.is_supplied()
                || conformer_id.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; mmff_variant, non_bonded_threshold, conformer_id, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_mmff_energy_gradient_keywords(
                    py,
                    mmff_variant,
                    non_bonded_threshold,
                    conformer_id,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_mmff_energy_gradient_params(py, params)
            }
        }
        /// Evaluate MMFF energy and coordinate gradients for the selected conformer and MMFF variant without minimizing it. Uses the explicitly supplied parameter form.
        fn mmff_energy_gradient_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::MmffEvaluationParams,
        ) -> PyResult<Option<super::numerical_results::MmffEnergyGradient>> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| self.inner.mmff_energy_gradient_with_params(&params))
                .map(|value| {
                    value.map(|inner| super::numerical_results::MmffEnergyGradient { inner })
                })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_energy_gradient"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_mmff_energy_gradient_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::MmffEvaluationParams>,
        ) -> PyResult<Option<super::numerical_results::MmffEnergyGradient>> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.mmff_energy_gradient_with_params(&value)),
                None => py.detach(|| self.inner.mmff_energy_gradient()),
            }
            .map(|value| value.map(|inner| super::numerical_results::MmffEnergyGradient { inner }))
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.mmff_energy_gradient"}"#]
        #[pyo3(signature=(*,mmff_variant=None,non_bonded_threshold=None,conformer_id=None,ignore_interfragment_interactions=None))]
        fn _ff_mmff_energy_gradient_keywords(
            &self,
            py: Python<'_>,
            mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
            non_bonded_threshold: Option<f64>,
            conformer_id: Option<usize>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<Option<super::numerical_results::MmffEnergyGradient>> {
            let params = super::forcefield_calls::mmff_energy_gradient_params(
                mmff_variant,
                non_bonded_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            );
            self._ff_mmff_energy_gradient_params(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,max_iterations=crate::arguments::Argument::Missing,vdw_threshold=crate::arguments::Argument::Missing,conformer_id=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Return a new molecule and UFF minimization outcome for the selected conformer, leaving the source unchanged.
        fn with_uff_optimized(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::UffOptimizationParams>>,
            >,
            max_iterations: crate::arguments::Argument<Option<i32>>,
            vdw_threshold: crate::arguments::Argument<Option<f64>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::numerical_results::UffOptimizationResult> {
            let keyword_form = max_iterations.is_supplied()
                || vdw_threshold.is_supplied()
                || conformer_id.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; max_iterations, vdw_threshold, conformer_id, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_with_uff_optimized_keywords(
                    py,
                    max_iterations,
                    vdw_threshold,
                    conformer_id,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_with_uff_optimized_params(py, params)
            }
        }
        /// Return a new molecule and UFF minimization outcome for the selected conformer, leaving the source unchanged. Uses the explicitly supplied parameter form.
        fn with_uff_optimized_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::UffOptimizationParams,
        ) -> PyResult<super::numerical_results::UffOptimizationResult> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| self.inner.with_uff_optimized_with_params(&params))
                .map(|inner| super::numerical_results::UffOptimizationResult { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_uff_optimized"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_with_uff_optimized_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::UffOptimizationParams>,
        ) -> PyResult<super::numerical_results::UffOptimizationResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.with_uff_optimized_with_params(&value)),
                None => py.detach(|| self.inner.with_uff_optimized()),
            }
            .map(|inner| super::numerical_results::UffOptimizationResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_uff_optimized"}"#]
        #[pyo3(signature=(*,max_iterations=None,vdw_threshold=None,conformer_id=None,ignore_interfragment_interactions=None))]
        fn _ff_with_uff_optimized_keywords(
            &self,
            py: Python<'_>,
            max_iterations: Option<i32>,
            vdw_threshold: Option<f64>,
            conformer_id: Option<usize>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<super::numerical_results::UffOptimizationResult> {
            let params = super::forcefield_calls::with_uff_optimized_params(
                max_iterations,
                vdw_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            );
            self._ff_with_uff_optimized_params(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,num_threads=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing,vdw_threshold=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Optimize the molecule's conformers with UFF and return the new molecule with per-conformer outcomes.
        fn with_uff_optimized_conformers(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::UffConformerOptimizationParams>>,
            >,
            num_threads: crate::arguments::Argument<Option<i32>>,
            max_iterations: crate::arguments::Argument<Option<i32>>,
            vdw_threshold: crate::arguments::Argument<Option<f64>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::numerical_results::UffConformerOptimizationResult> {
            let keyword_form = num_threads.is_supplied()
                || max_iterations.is_supplied()
                || vdw_threshold.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; num_threads, max_iterations, vdw_threshold, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_with_uff_optimized_conformers_keywords(
                    py,
                    num_threads,
                    max_iterations,
                    vdw_threshold,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_with_uff_optimized_conformers_params(py, params)
            }
        }
        /// Optimize the molecule's conformers with UFF and return the new molecule with per-conformer outcomes. Uses the explicitly supplied parameter form.
        fn with_uff_optimized_conformers_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::UffConformerOptimizationParams,
        ) -> PyResult<super::numerical_results::UffConformerOptimizationResult> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .with_uff_optimized_conformers_with_params(&params)
            })
            .map(|inner| super::numerical_results::UffConformerOptimizationResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_uff_optimized_conformers"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_with_uff_optimized_conformers_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::UffConformerOptimizationParams>,
        ) -> PyResult<super::numerical_results::UffConformerOptimizationResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => {
                    py.detach(|| self.inner.with_uff_optimized_conformers_with_params(&value))
                }
                None => py.detach(|| self.inner.with_uff_optimized_conformers()),
            }
            .map(|inner| super::numerical_results::UffConformerOptimizationResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_uff_optimized_conformers"}"#]
        #[pyo3(signature=(*,num_threads=None,max_iterations=None,vdw_threshold=None,ignore_interfragment_interactions=None))]
        fn _ff_with_uff_optimized_conformers_keywords(
            &self,
            py: Python<'_>,
            num_threads: Option<i32>,
            max_iterations: Option<i32>,
            vdw_threshold: Option<f64>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<super::numerical_results::UffConformerOptimizationResult> {
            let params = super::forcefield_calls::with_uff_optimized_conformers_params(
                num_threads,
                max_iterations,
                vdw_threshold,
                ignore_interfragment_interactions,
            );
            self._ff_with_uff_optimized_conformers_params(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,mmff_variant=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing,non_bonded_threshold=crate::arguments::Argument::Missing,conformer_id=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Return a new molecule and MMFF minimization outcome for the selected conformer, leaving the source unchanged.
        fn with_mmff_optimized(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::MmffOptimizationParams>>,
            >,
            mmff_variant: crate::arguments::Argument<Option<crate::leaf_values::MmffVariantInput>>,
            max_iterations: crate::arguments::Argument<Option<i32>>,
            non_bonded_threshold: crate::arguments::Argument<Option<f64>>,
            conformer_id: crate::arguments::Argument<Option<usize>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeResult> {
            let keyword_form = mmff_variant.is_supplied()
                || max_iterations.is_supplied()
                || non_bonded_threshold.is_supplied()
                || conformer_id.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; mmff_variant, max_iterations, non_bonded_threshold, conformer_id, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_with_mmff_optimized_keywords(
                    py,
                    mmff_variant,
                    max_iterations,
                    non_bonded_threshold,
                    conformer_id,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_with_mmff_optimized_params(py, params)
            }
        }
        /// Return a new molecule and MMFF minimization outcome for the selected conformer, leaving the source unchanged. Uses the explicitly supplied parameter form.
        fn with_mmff_optimized_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::MmffOptimizationParams,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeResult> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| self.inner.with_mmff_optimized_with_params(&params))
                .map(|inner| super::numerical_results::MmffOptimizeMoleculeResult { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_mmff_optimized"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_with_mmff_optimized_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::MmffOptimizationParams>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| self.inner.with_mmff_optimized_with_params(&value)),
                None => py.detach(|| self.inner.with_mmff_optimized()),
            }
            .map(|inner| super::numerical_results::MmffOptimizeMoleculeResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_mmff_optimized"}"#]
        #[pyo3(signature=(*,mmff_variant=None,max_iterations=None,non_bonded_threshold=None,conformer_id=None,ignore_interfragment_interactions=None))]
        fn _ff_with_mmff_optimized_keywords(
            &self,
            py: Python<'_>,
            mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
            max_iterations: Option<i32>,
            non_bonded_threshold: Option<f64>,
            conformer_id: Option<usize>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeResult> {
            let params = super::forcefield_calls::with_mmff_optimized_params(
                mmff_variant,
                max_iterations,
                non_bonded_threshold,
                conformer_id,
                ignore_interfragment_interactions,
            );
            self._ff_with_mmff_optimized_params(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,num_threads=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing,mmff_variant=crate::arguments::Argument::Missing,non_bonded_threshold=crate::arguments::Argument::Missing,ignore_interfragment_interactions=crate::arguments::Argument::Missing))]
        /// Optimize the molecule's conformers with MMFF and return the new molecule with per-conformer outcomes.
        fn with_mmff_optimized_conformers(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::params::MmffConformerOptimizationParams>>,
            >,
            num_threads: crate::arguments::Argument<Option<i32>>,
            max_iterations: crate::arguments::Argument<Option<i32>>,
            mmff_variant: crate::arguments::Argument<Option<crate::leaf_values::MmffVariantInput>>,
            non_bonded_threshold: crate::arguments::Argument<Option<f64>>,
            ignore_interfragment_interactions: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeConfsResult> {
            let keyword_form = num_threads.is_supplied()
                || max_iterations.is_supplied()
                || mmff_variant.is_supplied()
                || non_bonded_threshold.is_supplied()
                || ignore_interfragment_interactions.is_supplied();
            crate::arguments::configuration_arguments!(params; num_threads, max_iterations, mmff_variant, non_bonded_threshold, ignore_interfragment_interactions);
            let params = params.as_deref();
            if keyword_form {
                self._ff_with_mmff_optimized_conformers_keywords(
                    py,
                    num_threads,
                    max_iterations,
                    mmff_variant,
                    non_bonded_threshold,
                    ignore_interfragment_interactions,
                )
            } else {
                self._ff_with_mmff_optimized_conformers_params(py, params)
            }
        }
        /// Optimize the molecule's conformers with MMFF and return the new molecule with per-conformer outcomes. Uses the explicitly supplied parameter form.
        fn with_mmff_optimized_conformers_with_params(
            &self,
            py: Python<'_>,
            params: &super::params::MmffConformerOptimizationParams,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeConfsResult> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            py.detach(|| {
                self.inner
                    .with_mmff_optimized_conformers_with_params(&params)
            })
            .map(|inner| super::numerical_results::MmffOptimizeMoleculeConfsResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_mmff_optimized_conformers"}"#]
        #[pyo3(signature=(params=None))]
        fn _ff_with_mmff_optimized_conformers_params(
            &self,
            py: Python<'_>,
            params: Option<&super::params::MmffConformerOptimizationParams>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeConfsResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            match params {
                Some(value) => py.detach(|| {
                    self.inner
                        .with_mmff_optimized_conformers_with_params(&value)
                }),
                None => py.detach(|| self.inner.with_mmff_optimized_conformers()),
            }
            .map(|inner| super::numerical_results::MmffOptimizeMoleculeConfsResult { inner })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_mmff_optimized_conformers"}"#]
        #[pyo3(signature=(*,num_threads=None,max_iterations=None,mmff_variant=None,non_bonded_threshold=None,ignore_interfragment_interactions=None))]
        fn _ff_with_mmff_optimized_conformers_keywords(
            &self,
            py: Python<'_>,
            num_threads: Option<i32>,
            max_iterations: Option<i32>,
            mmff_variant: Option<crate::leaf_values::MmffVariantInput>,
            non_bonded_threshold: Option<f64>,
            ignore_interfragment_interactions: Option<bool>,
        ) -> PyResult<super::numerical_results::MmffOptimizeMoleculeConfsResult> {
            let params = super::forcefield_calls::with_mmff_optimized_conformers_params(
                num_threads,
                max_iterations,
                mmff_variant,
                non_bonded_threshold,
                ignore_interfragment_interactions,
            );
            self._ff_with_mmff_optimized_conformers_params(py, Some(&params))
        }
        #[pyo3(signature=(reference,params=crate::arguments::Argument::Missing,*,probe_conformer_id=crate::arguments::Argument::Missing,reference_conformer_id=crate::arguments::Argument::Missing,atom_map=crate::arguments::Argument::Missing,weights=crate::arguments::Argument::Missing,reflect=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing))]
        /// Return a new molecule aligned to the reference coordinates, preserving the source molecule.
        fn with_alignment_to(
            &self,
            py: Python<'_>,
            reference: &Self,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::collection_params::AlignmentParameters>>,
            >,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_map: crate::arguments::Argument<
                Option<Vec<super::collection_params::AlignmentAtomMap>>,
            >,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<(Self, super::alignment_values::AlignmentResult)> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_map.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_map, weights, reflect, max_iterations);
            let params = params.as_deref();
            if keyword_form {
                self._alignment_with_alignment_to_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_map,
                    weights,
                    reflect,
                    max_iterations,
                )
            } else {
                self._alignment_with_alignment_to_params(py, reference, params)
            }
        }
        /// Return a new molecule aligned to the reference coordinates, preserving the source molecule. Uses the explicitly supplied parameter form.
        fn with_alignment_to_with_params(
            &self,
            py: Python<'_>,
            reference: &Self,
            params: &super::collection_params::AlignmentParameters,
        ) -> PyResult<(Self, super::alignment_values::AlignmentResult)> {
            let params = params.to_rust(py)?;
            py.detach(|| {
                self.inner
                    .with_alignment_to_with_params(&reference.inner, &params)
            })
            .map(|(inner, report)| {
                (
                    Self { inner },
                    super::alignment_values::AlignmentResult { inner: report },
                )
            })
            .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_alignment_to"}"#]
        #[pyo3(signature=(reference,params=None))]
        fn _alignment_with_alignment_to_params(
            &self,
            py: Python<'_>,
            reference: &Self,
            params: Option<&super::collection_params::AlignmentParameters>,
        ) -> PyResult<(Self, super::alignment_values::AlignmentResult)> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._alignment_with_alignment_to_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_alignment_to"}"#]
        #[pyo3(signature=(reference,*,probe_conformer_id=None,reference_conformer_id=None,atom_map=None,weights=None,reflect=None,max_iterations=None))]
        fn _alignment_with_alignment_to_keywords(
            &self,
            py: Python<'_>,
            reference: &Self,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_map: Option<Vec<super::collection_params::AlignmentAtomMap>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
        ) -> PyResult<(Self, super::alignment_values::AlignmentResult)> {
            let params = super::alignment_queries::transform_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_map,
                weights,
                reflect,
                max_iterations,
            );
            self._alignment_with_alignment_to_params_native(py, reference, Some(&params))
        }

        #[pyo3(signature=(reference,params=crate::arguments::Argument::Missing,*,probe_conformer_id=crate::arguments::Argument::Missing,reference_conformer_id=crate::arguments::Argument::Missing,atom_map=crate::arguments::Argument::Missing,weights=crate::arguments::Argument::Missing,reflect=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing))]
        /// Align this molecule to the reference coordinates in place; COW isolates other molecular values.
        fn align_to_(
            &mut self,
            py: Python<'_>,
            reference: &Self,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::collection_params::AlignmentParameters>>,
            >,
            probe_conformer_id: crate::arguments::Argument<Option<i32>>,
            reference_conformer_id: crate::arguments::Argument<Option<i32>>,
            atom_map: crate::arguments::Argument<
                Option<Vec<super::collection_params::AlignmentAtomMap>>,
            >,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<super::alignment_values::AlignmentResult> {
            let keyword_form = probe_conformer_id.is_supplied()
                || reference_conformer_id.is_supplied()
                || atom_map.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied();
            crate::arguments::configuration_arguments!(params; probe_conformer_id, reference_conformer_id, atom_map, weights, reflect, max_iterations);
            let params = params.as_deref();
            if keyword_form {
                self._alignment_with_alignment_to_inplace_keywords(
                    py,
                    reference,
                    probe_conformer_id,
                    reference_conformer_id,
                    atom_map,
                    weights,
                    reflect,
                    max_iterations,
                )
            } else {
                self._alignment_with_alignment_to_inplace_params(py, reference, params)
            }
        }
        /// Align this molecule's selected conformer to a reference in place using explicit alignment parameters.
        fn align_to_with_params_(
            &mut self,
            py: Python<'_>,
            reference: &Self,
            params: &super::collection_params::AlignmentParameters,
        ) -> PyResult<super::alignment_values::AlignmentResult> {
            let params = params.to_rust(py)?;
            py.detach(|| self.inner.align_to_with_params_(&reference.inner, &params))
                .map(|inner| super::alignment_values::AlignmentResult { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.align_to_"}"#]
        #[pyo3(signature=(reference,params=None))]
        fn _alignment_with_alignment_to_inplace_params(
            &mut self,
            py: Python<'_>,
            reference: &Self,
            params: Option<&super::collection_params::AlignmentParameters>,
        ) -> PyResult<super::alignment_values::AlignmentResult> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._alignment_with_alignment_to_inplace_params_native(py, reference, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.align_to_"}"#]
        #[pyo3(signature=(reference,*,probe_conformer_id=None,reference_conformer_id=None,atom_map=None,weights=None,reflect=None,max_iterations=None))]
        fn _alignment_with_alignment_to_inplace_keywords(
            &mut self,
            py: Python<'_>,
            reference: &Self,
            probe_conformer_id: Option<i32>,
            reference_conformer_id: Option<i32>,
            atom_map: Option<Vec<super::collection_params::AlignmentAtomMap>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
        ) -> PyResult<super::alignment_values::AlignmentResult> {
            let params = super::alignment_queries::transform_keywords(
                probe_conformer_id,
                reference_conformer_id,
                atom_map,
                weights,
                reflect,
                max_iterations,
            );
            self._alignment_with_alignment_to_inplace_params_native(py, reference, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,atom_indices=crate::arguments::Argument::Missing,conformer_ids=crate::arguments::Argument::Missing,weights=crate::arguments::Argument::Missing,reflect=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing))]
        /// Return a new molecule with its conformers aligned according to the supplied parameters.
        fn with_aligned_conformers(
            &self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::numerical_params::ConformerAlignmentParameters>>,
            >,
            atom_indices: crate::arguments::Argument<Option<Vec<usize>>>,
            conformer_ids: crate::arguments::Argument<Option<Vec<usize>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<(Self, super::alignment_values::ConformerAlignmentReport)> {
            let keyword_form = atom_indices.is_supplied()
                || conformer_ids.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied();
            crate::arguments::configuration_arguments!(params; atom_indices, conformer_ids, weights, reflect, max_iterations);
            let params = params.as_deref();
            if keyword_form {
                self._alignment_with_aligned_conformers_keywords(
                    py,
                    atom_indices,
                    conformer_ids,
                    weights,
                    reflect,
                    max_iterations,
                )
            } else {
                self._alignment_with_aligned_conformers_params(py, params)
            }
        }
        /// Return a new molecule with its conformers aligned according to the supplied parameters. Uses the explicitly supplied parameter form.
        fn with_aligned_conformers_with_params(
            &self,
            py: Python<'_>,
            params: &super::numerical_params::ConformerAlignmentParameters,
        ) -> PyResult<(Self, super::alignment_values::ConformerAlignmentReport)> {
            let params = params.to_rust(py)?;
            py.detach(|| self.inner.with_aligned_conformers_with_params(&params))
                .map(|(inner, report)| {
                    (
                        Self { inner },
                        super::alignment_values::ConformerAlignmentReport { inner: report },
                    )
                })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_aligned_conformers"}"#]
        #[pyo3(signature=(params=None))]
        fn _alignment_with_aligned_conformers_params(
            &self,
            py: Python<'_>,
            params: Option<&super::numerical_params::ConformerAlignmentParameters>,
        ) -> PyResult<(Self, super::alignment_values::ConformerAlignmentReport)> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._alignment_with_aligned_conformers_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.with_aligned_conformers"}"#]
        #[pyo3(signature=(*,atom_indices=None,conformer_ids=None,weights=None,reflect=None,max_iterations=None))]
        fn _alignment_with_aligned_conformers_keywords(
            &self,
            py: Python<'_>,
            atom_indices: Option<Vec<usize>>,
            conformer_ids: Option<Vec<usize>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
        ) -> PyResult<(Self, super::alignment_values::ConformerAlignmentReport)> {
            let params = super::alignment_calls::conformer_keywords(
                atom_indices,
                conformer_ids,
                weights,
                reflect,
                max_iterations,
            );
            self._alignment_with_aligned_conformers_params_native(py, Some(&params))
        }

        #[pyo3(signature=(params=crate::arguments::Argument::Missing,*,atom_indices=crate::arguments::Argument::Missing,conformer_ids=crate::arguments::Argument::Missing,weights=crate::arguments::Argument::Missing,reflect=crate::arguments::Argument::Missing,max_iterations=crate::arguments::Argument::Missing))]
        /// Align the stored conformers in place according to the supplied parameters.
        fn align_conformers_(
            &mut self,
            py: Python<'_>,
            params: crate::arguments::Argument<
                Option<PyRef<'_, super::numerical_params::ConformerAlignmentParameters>>,
            >,
            atom_indices: crate::arguments::Argument<Option<Vec<usize>>>,
            conformer_ids: crate::arguments::Argument<Option<Vec<usize>>>,
            weights: crate::arguments::Argument<Option<Vec<f64>>>,
            reflect: crate::arguments::Argument<Option<bool>>,
            max_iterations: crate::arguments::Argument<Option<u32>>,
        ) -> PyResult<super::alignment_values::ConformerAlignmentReport> {
            let keyword_form = atom_indices.is_supplied()
                || conformer_ids.is_supplied()
                || weights.is_supplied()
                || reflect.is_supplied()
                || max_iterations.is_supplied();
            crate::arguments::configuration_arguments!(params; atom_indices, conformer_ids, weights, reflect, max_iterations);
            let params = params.as_deref();
            if keyword_form {
                self._alignment_with_aligned_conformers_inplace_keywords(
                    py,
                    atom_indices,
                    conformer_ids,
                    weights,
                    reflect,
                    max_iterations,
                )
            } else {
                self._alignment_with_aligned_conformers_inplace_params(py, params)
            }
        }
        /// Align this molecule's conformers in place using explicit alignment parameters.
        fn align_conformers_with_params_(
            &mut self,
            py: Python<'_>,
            params: &super::numerical_params::ConformerAlignmentParameters,
        ) -> PyResult<super::alignment_values::ConformerAlignmentReport> {
            let params = params.to_rust(py)?;
            py.detach(|| self.inner.align_conformers_with_params_(&params))
                .map(|inner| super::alignment_values::ConformerAlignmentReport { inner })
                .map_err(super::operation_errors::OperationError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.align_conformers_"}"#]
        #[pyo3(signature=(params=None))]
        fn _alignment_with_aligned_conformers_inplace_params(
            &mut self,
            py: Python<'_>,
            params: Option<&super::numerical_params::ConformerAlignmentParameters>,
        ) -> PyResult<super::alignment_values::ConformerAlignmentReport> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            self._alignment_with_aligned_conformers_inplace_params_native(py, params.as_ref())
        }
        #[doc = r#"{"role":"overload","target":"Molecule.align_conformers_"}"#]
        #[pyo3(signature=(*,atom_indices=None,conformer_ids=None,weights=None,reflect=None,max_iterations=None))]
        fn _alignment_with_aligned_conformers_inplace_keywords(
            &mut self,
            py: Python<'_>,
            atom_indices: Option<Vec<usize>>,
            conformer_ids: Option<Vec<usize>>,
            weights: Option<Vec<f64>>,
            reflect: Option<bool>,
            max_iterations: Option<u32>,
        ) -> PyResult<super::alignment_values::ConformerAlignmentReport> {
            let params = super::alignment_calls::conformer_keywords(
                atom_indices,
                conformer_ids,
                weights,
                reflect,
                max_iterations,
            );
            self._alignment_with_aligned_conformers_inplace_params_native(py, Some(&params))
        }
        #[getter]
        /// Number of stored 3D conformers, not the count of attempted embeddings.
        fn num_3d_conformers(&self) -> usize {
            self.inner.num_3d_conformers()
        }
        #[doc = r#"{"role":"return_type","imports":["numpy"],"type":"numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
        /// Return an independent NumPy distance-geometry bounds matrix in atom-index order.
        fn dg_bounds_matrix<'py>(
            &self,
            py: Python<'py>,
        ) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
            super::conformer_queries::bounds(self, py)
        }
        #[pyo3(signature=(options=crate::arguments::Argument::Missing,*,try_embedding=crate::arguments::Argument::Missing,only_unassigned=crate::arguments::Argument::Missing,max_isomers=crate::arguments::Argument::Missing,random_source=crate::arguments::Argument::Missing,unique=crate::arguments::Argument::Missing,only_stereo_groups=crate::arguments::Argument::Missing))]
        /// Return the arbitrary-precision count from the native stereoisomer enumeration policy.
        fn stereoisomer_count(
            &self,
            py: Python<'_>,
            options: crate::arguments::Argument<
                Option<Py<super::stereoisomer_options::StereoisomerOptions>>,
            >,
            try_embedding: crate::arguments::Argument<Option<bool>>,
            only_unassigned: crate::arguments::Argument<Option<bool>>,
            max_isomers: crate::arguments::Argument<Option<usize>>,
            random_source: crate::arguments::Argument<
                Option<Py<super::stereoisomer_options::StereoisomerRandomSource>>,
            >,
            unique: crate::arguments::Argument<Option<bool>>,
            only_stereo_groups: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<num_bigint::BigUint> {
            let keyword_form = try_embedding.is_supplied()
                || only_unassigned.is_supplied()
                || max_isomers.is_supplied()
                || random_source.is_supplied()
                || unique.is_supplied()
                || only_stereo_groups.is_supplied();
            crate::arguments::configuration_arguments!(options; try_embedding, only_unassigned, max_isomers, random_source, unique, only_stereo_groups);
            if keyword_form {
                self._stereoisomer_count_keywords(
                    py,
                    try_embedding,
                    only_unassigned,
                    max_isomers,
                    random_source,
                    unique,
                    only_stereo_groups,
                )
            } else {
                self._stereoisomer_count_options(py, options)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.stereoisomer_count"}"#]
        #[pyo3(signature=(options=None))]
        fn _stereoisomer_count_options(
            &self,
            py: Python<'_>,
            options: Option<Py<super::stereoisomer_options::StereoisomerOptions>>,
        ) -> PyResult<num_bigint::BigUint> {
            super::conformer_queries::count(py, self, options)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.stereoisomer_count"}"#]
        #[pyo3(signature=(*,try_embedding=None,only_unassigned=None,max_isomers=None,random_source=None,unique=None,only_stereo_groups=None))]
        fn _stereoisomer_count_keywords(
            &self,
            py: Python<'_>,
            try_embedding: Option<bool>,
            only_unassigned: Option<bool>,
            max_isomers: Option<usize>,
            random_source: Option<Py<super::stereoisomer_options::StereoisomerRandomSource>>,
            unique: Option<bool>,
            only_stereo_groups: Option<bool>,
        ) -> PyResult<num_bigint::BigUint> {
            super::conformer_queries::count_keywords(
                py,
                self,
                try_embedding,
                only_unassigned,
                max_isomers,
                random_source,
                unique,
                only_stereo_groups,
            )
        }
        /// Return the arbitrary-precision count from the native stereoisomer enumeration policy. Uses the explicitly supplied options.
        fn stereoisomer_count_with_options(
            &self,
            py: Python<'_>,
            options: Py<super::stereoisomer_options::StereoisomerOptions>,
        ) -> PyResult<num_bigint::BigUint> {
            self._stereoisomer_count_options(py, Some(options))
        }
        #[new]
        /// Construct an empty molecular value.
        fn py_new() -> Self {
            Self::new()
        }
        #[staticmethod]
        /// Construct an empty molecular value.
        fn new() -> Self {
            Self {
                inner: ck::Molecule::new(),
            }
        }
        #[staticmethod]
        #[pyo3(signature = (text, params=Argument::Missing, *, sanitize=Argument::Missing, allow_cxsmiles=Argument::Missing, strict_cxsmiles=Argument::Missing, parse_name=Argument::Missing, remove_hs=Argument::Missing, skip_cleanup=Argument::Missing, debug_parse=Argument::Missing, replacements=Argument::Missing))]
        /// Parse SMILES into a new molecule. Supply SmilesParseParams or keyword options, not both; sanitization follows the selected parsing policy.
        fn from_smiles(
            text: &str,
            params: Argument<Option<PyRef<'_, super::extension::SmilesParseParams>>>,
            sanitize: Argument<Option<bool>>,
            allow_cxsmiles: Argument<Option<bool>>,
            strict_cxsmiles: Argument<Option<bool>>,
            parse_name: Argument<Option<bool>>,
            remove_hs: Argument<Option<bool>>,
            skip_cleanup: Argument<Option<bool>>,
            debug_parse: Argument<Option<bool>>,
            replacements: Argument<Option<std::collections::BTreeMap<String, String>>>,
        ) -> PyResult<Self> {
            crate::arguments::configuration_arguments!(params; sanitize, allow_cxsmiles, strict_cxsmiles, parse_name, remove_hs, skip_cleanup, debug_parse, replacements);
            let params = params.as_deref();
            if let Some(params) = params {
                Self::_from_smiles_params(text, Some(params))
            } else {
                Self::_from_smiles_keywords(
                    text,
                    sanitize,
                    allow_cxsmiles,
                    strict_cxsmiles,
                    parse_name,
                    remove_hs,
                    skip_cleanup,
                    debug_parse,
                    replacements,
                )
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_smiles"}"#]
        #[pyo3(signature = (text, params=None))]
        fn _from_smiles_params(
            text: &str,
            params: Option<&super::extension::SmilesParseParams>,
        ) -> PyResult<Self> {
            super::extension::_mol_from_smiles_params(text, params)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_smiles"}"#]
        #[pyo3(signature = (text, *, sanitize=None, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, remove_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
        fn _from_smiles_keywords(
            text: &str,
            sanitize: Option<bool>,
            allow_cxsmiles: Option<bool>,
            strict_cxsmiles: Option<bool>,
            parse_name: Option<bool>,
            remove_hs: Option<bool>,
            skip_cleanup: Option<bool>,
            debug_parse: Option<bool>,
            replacements: Option<std::collections::BTreeMap<String, String>>,
        ) -> PyResult<Self> {
            super::extension::_mol_from_smiles_keywords(
                text,
                sanitize,
                allow_cxsmiles,
                strict_cxsmiles,
                parse_name,
                remove_hs,
                skip_cleanup,
                debug_parse,
                replacements,
            )
        }
        #[staticmethod]
        /// Parse SMILES into a new molecule. Supply SmilesParseParams or keyword options, not both; sanitization follows the selected parsing policy. Uses the explicitly supplied parameter form.
        fn from_smiles_with_params(
            text: &str,
            params: &super::extension::SmilesParseParams,
        ) -> PyResult<Self> {
            Self::_from_smiles_params(text, Some(params))
        }
        #[staticmethod]
        #[pyo3(signature = (text, params=crate::arguments::Argument::Missing, *, sanitize=crate::arguments::Argument::Missing, remove_hs=crate::arguments::Argument::Missing))]
        /// Parse InChI into a new molecular value using the selected sanitization and hydrogen-removal policy.
        fn from_inchi(
            text: &str,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::InchiReadParams>>>,
            sanitize: crate::arguments::Argument<Option<bool>>,
            remove_hs: crate::arguments::Argument<Option<bool>>,
        ) -> PyResult<Self> {
            let keyword_form = sanitize.is_supplied() || remove_hs.is_supplied();
            crate::arguments::configuration_arguments!(params; sanitize, remove_hs);
            let params = params.as_deref();
            if keyword_form {
                Self::_from_inchi_keywords(text, sanitize, remove_hs)
            } else {
                Self::_from_inchi_params(text, params)
            }
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_inchi"}"#]
        #[pyo3(signature = (text, params=None))]
        pub(crate) fn _from_inchi_params(
            text: &str,
            params: Option<&super::params::InchiReadParams>,
        ) -> PyResult<Self> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults = ck::InchiReadParams::default();
            ck::Molecule::from_inchi_with_params(text, params.as_ref().unwrap_or(&defaults))
                .map(|inner| Self { inner })
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[staticmethod]
        #[doc = r#"{"role":"overload","target":"Molecule.from_inchi"}"#]
        #[pyo3(signature = (text, *, sanitize=None, remove_hs=None))]
        pub(crate) fn _from_inchi_keywords(
            text: &str,
            sanitize: Option<bool>,
            remove_hs: Option<bool>,
        ) -> PyResult<Self> {
            let defaults = ck::InchiReadParams::default();
            let params = super::params::InchiReadParams {
                inner: crate::parameters::NativeField::Owned(ck::InchiReadParams {
                    sanitize: sanitize.unwrap_or(defaults.sanitize),
                    remove_hs: remove_hs.unwrap_or(defaults.remove_hs),
                }),
            };
            Self::_from_inchi_params(text, Some(&params))
        }
        #[staticmethod]
        /// Parse InChI into a new molecular value using the selected sanitization and hydrogen-removal policy. Uses the explicitly supplied parameter form.
        fn from_inchi_with_params(
            text: &str,
            params: &super::params::InchiReadParams,
        ) -> PyResult<Self> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            ck::Molecule::from_inchi_with_params(text, &params)
                .map(|inner| Self { inner })
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, options=crate::arguments::Argument::Missing))]
        /// Generate an InChI string using the selected native InChI options.
        fn to_inchi(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::InchiWriteParams>>>,
            options: crate::arguments::Argument<Option<String>>,
        ) -> PyResult<String> {
            let keyword_form = options.is_supplied();
            crate::arguments::configuration_arguments!(params; options);
            let params = params.as_deref();
            if keyword_form {
                self._to_inchi_keywords(options)
            } else {
                self._to_inchi_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi"}"#]
        #[pyo3(signature = (params=None))]
        fn _to_inchi_params(
            &self,
            params: Option<&super::params::InchiWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults = ck::InchiWriteParams::default();
            self.inner
                .to_inchi_with_params(params.as_ref().unwrap_or(&defaults))
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi"}"#]
        #[pyo3(signature = (*, options=None))]
        fn _to_inchi_keywords(&self, options: Option<String>) -> PyResult<String> {
            let params = super::params::InchiWriteParams {
                inner: crate::parameters::NativeField::Owned(ck::InchiWriteParams {
                    options: options.unwrap_or_else(|| ck::InchiWriteParams::default().options),
                }),
            };
            self._to_inchi_params(Some(&params))
        }
        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, options=crate::arguments::Argument::Missing))]
        /// Generate InChI with the native result metadata and diagnostic information.
        fn to_inchi_result(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::InchiWriteParams>>>,
            options: crate::arguments::Argument<Option<String>>,
        ) -> PyResult<super::inchi_values::InchiResult> {
            let keyword_form = options.is_supplied();
            crate::arguments::configuration_arguments!(params; options);
            let params = params.as_deref();
            if keyword_form {
                self._to_inchi_result_keywords(options)
            } else {
                self._to_inchi_result_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi_result"}"#]
        #[pyo3(signature = (params=None))]
        fn _to_inchi_result_params(
            &self,
            params: Option<&super::params::InchiWriteParams>,
        ) -> PyResult<super::inchi_values::InchiResult> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults = ck::InchiWriteParams::default();
            self.inner
                .to_inchi_result_with_params(params.as_ref().unwrap_or(&defaults))
                .map(|inner| super::inchi_values::InchiResult { inner })
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi_result"}"#]
        #[pyo3(signature = (*, options=None))]
        fn _to_inchi_result_keywords(
            &self,
            options: Option<String>,
        ) -> PyResult<super::inchi_values::InchiResult> {
            let params = super::params::InchiWriteParams {
                inner: crate::parameters::NativeField::Owned(ck::InchiWriteParams {
                    options: options.unwrap_or_else(|| ck::InchiWriteParams::default().options),
                }),
            };
            self._to_inchi_result_params(Some(&params))
        }
        /// Generate InChI with the native result metadata and diagnostic information. Uses the explicitly supplied parameter form.
        fn to_inchi_result_with_params(
            &self,
            params: &super::params::InchiWriteParams,
        ) -> PyResult<super::inchi_values::InchiResult> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            self.inner
                .to_inchi_result_with_params(&params)
                .map(|inner| super::inchi_values::InchiResult { inner })
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        /// Generate an InChI string using the selected native InChI options. Uses the explicitly supplied parameter form.
        fn to_inchi_with_params(
            &self,
            params: &super::params::InchiWriteParams,
        ) -> PyResult<String> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            self.inner
                .to_inchi_with_params(&params)
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[pyo3(signature = (params=crate::arguments::Argument::Missing, *, options=crate::arguments::Argument::Missing))]
        /// Generate the InChIKey for this molecule with the selected InChI options.
        fn to_inchi_key(
            &self,
            params: crate::arguments::Argument<Option<PyRef<'_, super::params::InchiWriteParams>>>,
            options: crate::arguments::Argument<Option<String>>,
        ) -> PyResult<String> {
            let keyword_form = options.is_supplied();
            crate::arguments::configuration_arguments!(params; options);
            let params = params.as_deref();
            if keyword_form {
                self._to_inchi_key_keywords(options)
            } else {
                self._to_inchi_key_params(params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi_key"}"#]
        #[pyo3(signature = (params=None))]
        fn _to_inchi_key_params(
            &self,
            params: Option<&super::params::InchiWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let defaults = ck::InchiWriteParams::default();
            self.inner
                .to_inchi_key_with_params(params.as_ref().unwrap_or(&defaults))
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_inchi_key"}"#]
        #[pyo3(signature = (*, options=None))]
        fn _to_inchi_key_keywords(&self, options: Option<String>) -> PyResult<String> {
            let params = super::params::InchiWriteParams {
                inner: crate::parameters::NativeField::Owned(ck::InchiWriteParams {
                    options: options.unwrap_or_else(|| ck::InchiWriteParams::default().options),
                }),
            };
            self._to_inchi_key_params(Some(&params))
        }
        /// Generate the InChIKey for this molecule with the selected InChI options. Uses the explicitly supplied parameter form.
        fn to_inchi_key_with_params(
            &self,
            params: &super::params::InchiWriteParams,
        ) -> PyResult<String> {
            let params = pyo3::Python::attach(|py| params.to_rust(py))?;

            self.inner
                .to_inchi_key_with_params(&params)
                .map_err(super::inchi_values::InchiError::from_rust)
        }
        #[getter]
        /// Number of stored atom rows. Implicit hydrogens are not separate rows.
        fn num_atoms(&self) -> usize {
            self.inner.num_atoms()
        }
        #[getter]
        /// Number of stored bond rows.
        fn num_bonds(&self) -> usize {
            self.inner.num_bonds()
        }
        /// Calculate the canonical molecular graph hash.
        fn molecular_hash(&self) -> PyResult<u64> {
            self.inner
                .molecular_hash()
                .map_err(MoleculeHashError::from_rust)
        }
        /// Calculate the molecular graph hash using the supplied per-atom ranks.
        fn molecular_hash_with_ranks(&self, ranks: Vec<u32>) -> PyResult<u64> {
            self.inner
                .molecular_hash_with_ranks(&ranks)
                .map_err(MoleculeHashError::from_rust)
        }
        #[getter]
        /// Whether a 2D coordinate block is stored; does not generate coordinates.
        fn has_2d_coordinates(&self) -> bool {
            self.inner.has_2d_coordinates()
        }
        #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]] | None"}"#]
        /// Independent float64 array of shape (num_atoms, 2), in atom-index order,
        /// or None when no 2D coordinate block is stored. Does not generate a layout.
        #[getter]
        fn coordinates_2d<'py>(&self, py: Python<'py>) -> Option<Bound<'py, numpy::PyArray2<f64>>> {
            self.inner
                .coordinates_2d()
                .map(|rows| coordinates::array(py, rows))
        }
        #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int, int], numpy.dtype[numpy.float64]]"}"#]
        #[pyo3(signature = (conformer_id=0))]
        /// Return an independent float64 array of shape (num_atoms, 3) for conformer_id. An absent conformer raises Coordinate3DReadError.
        fn coordinates_3d<'py>(
            &self,
            py: Python<'py>,
            conformer_id: usize,
        ) -> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
            self.inner
                .coordinates_3d(conformer_id)
                .map(|rows| coordinates::array(py, rows))
                .map_err(errors::Coordinate3DReadError::from_rust)
        }
        #[getter]
        /// Return stored 3D conformer values in native storage order.
        fn conformers_3d(&self) -> Vec<Conformer3D> {
            self.inner
                .conformers_3d()
                .iter()
                .cloned()
                .map(|inner| Conformer3D { inner })
                .collect()
        }
        #[pyo3(signature = (params=Argument::Missing, *, n_bits=Argument::Missing))]
        /// Calculate MACCS keys using MaccsFingerprintParams or the n_bits keyword, not both.
        fn fingerprint_maccs(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, MaccsFingerprintParams>>>,
            n_bits: Argument<Option<usize>>,
        ) -> PyResult<Fingerprint> {
            let keyword_form = n_bits.is_supplied();
            crate::arguments::configuration_arguments!(params; n_bits);
            if keyword_form {
                self._fingerprint_maccs_keywords(n_bits)
            } else {
                self._fingerprint_maccs_params(py, params.as_deref())
            }
        }

        #[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
        #[pyo3(signature = (params=None))]
        fn _fingerprint_maccs_params(
            &self,
            py: Python<'_>,
            params: Option<&MaccsFingerprintParams>,
        ) -> PyResult<Fingerprint> {
            let effective = params
                .map(|params| params.to_rust(py))
                .transpose()?
                .unwrap_or_default();
            self.inner
                .fingerprint_maccs_with_params(&effective)
                .map(|inner| Fingerprint { inner })
                .map_err(errors::MaccsFingerprintError::from_rust)
        }

        #[doc = r#"{"role": "overload", "target": "Molecule.fingerprint_maccs"}"#]
        #[pyo3(signature = (*, n_bits=None))]
        fn _fingerprint_maccs_keywords(&self, n_bits: Option<usize>) -> PyResult<Fingerprint> {
            let params = ck::MaccsFingerprintParams {
                n_bits: n_bits.unwrap_or_else(|| ck::MaccsFingerprintParams::default().n_bits),
            };
            self.inner
                .fingerprint_maccs_with_params(&params)
                .map(|inner| Fingerprint { inner })
                .map_err(errors::MaccsFingerprintError::from_rust)
        }
        /// Calculate MACCS keys using MaccsFingerprintParams or the n_bits keyword, not both. Uses the explicitly supplied parameter form.
        fn fingerprint_maccs_with_params(
            &self,
            py: Python<'_>,
            params: &MaccsFingerprintParams,
        ) -> PyResult<Fingerprint> {
            self.inner
                .fingerprint_maccs_with_params(&params.to_rust(py)?)
                .map(|inner| Fingerprint { inner })
                .map_err(errors::MaccsFingerprintError::from_rust)
        }
        /// Return the raw 167-bit MACCS fingerprint with bit zero unused.
        /// The ordinary fingerprint_maccs entry projects keys 1..166 to 166 bits.
        fn fingerprint_maccs_raw(&self) -> PyResult<Fingerprint> {
            self.inner
                .fingerprint_maccs_raw()
                .map(|inner| Fingerprint { inner })
                .map_err(errors::MaccsFingerprintError::from_rust)
        }
        #[pyo3(signature=(params=Argument::Missing, *, isomeric_smiles=Argument::Missing, kekule=Argument::Missing, canonical=Argument::Missing, clean_stereo=Argument::Missing, rooted_at_atom=Argument::Missing, all_bonds_explicit=Argument::Missing, all_hydrogens_explicit=Argument::Missing, include_dative_bonds=Argument::Missing, ignore_atom_map_numbers=Argument::Missing))]
        /// Serialize SMILES with the selected canonicalization, stereochemistry and atom-order options. Supply SmilesWriteParams or keywords, not both.
        fn to_smiles(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, SmilesWriteParams>>>,
            isomeric_smiles: Argument<Option<bool>>,
            kekule: Argument<Option<bool>>,
            canonical: Argument<Option<bool>>,
            clean_stereo: Argument<Option<bool>>,
            rooted_at_atom: Argument<Option<usize>>,
            all_bonds_explicit: Argument<Option<bool>>,
            all_hydrogens_explicit: Argument<Option<bool>>,
            include_dative_bonds: Argument<Option<bool>>,
            ignore_atom_map_numbers: Argument<Option<bool>>,
        ) -> PyResult<String> {
            crate::arguments::configuration_arguments!(params; isomeric_smiles, kekule, canonical, clean_stereo, rooted_at_atom, all_bonds_explicit, all_hydrogens_explicit, include_dative_bonds, ignore_atom_map_numbers);
            let params = params.as_deref();
            let configured = isomeric_smiles.is_some()
                || kekule.is_some()
                || canonical.is_some()
                || clean_stereo.is_some()
                || rooted_at_atom.is_some()
                || all_bonds_explicit.is_some()
                || all_hydrogens_explicit.is_some()
                || include_dative_bonds.is_some()
                || ignore_atom_map_numbers.is_some();
            if configured {
                self._to_smiles_keywords(
                    py,
                    isomeric_smiles,
                    kekule,
                    canonical,
                    clean_stereo,
                    rooted_at_atom,
                    all_bonds_explicit,
                    all_hydrogens_explicit,
                    include_dative_bonds,
                    ignore_atom_map_numbers,
                )
            } else {
                self._to_smiles_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_smiles"}"#]
        #[pyo3(signature=(params=None))]
        fn _to_smiles_params(
            &self,
            py: Python<'_>,
            params: Option<&SmilesWriteParams>,
        ) -> PyResult<String> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let value = py.detach(|| match params {
                Some(params) => self.inner.to_smiles_with_params(&params),
                None => self.inner.to_smiles(),
            });
            super::smiles_writer::text(py, value)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_smiles"}"#]
        #[pyo3(signature=(*, isomeric_smiles=None, kekule=None, canonical=None, clean_stereo=None, rooted_at_atom=None, all_bonds_explicit=None, all_hydrogens_explicit=None, include_dative_bonds=None, ignore_atom_map_numbers=None))]
        fn _to_smiles_keywords(
            &self,
            py: Python<'_>,
            isomeric_smiles: Option<bool>,
            kekule: Option<bool>,
            canonical: Option<bool>,
            clean_stereo: Option<bool>,
            rooted_at_atom: Option<usize>,
            all_bonds_explicit: Option<bool>,
            all_hydrogens_explicit: Option<bool>,
            include_dative_bonds: Option<bool>,
            ignore_atom_map_numbers: Option<bool>,
        ) -> PyResult<String> {
            let params = super::smiles_writer::keywords(
                isomeric_smiles,
                kekule,
                canonical,
                clean_stereo,
                rooted_at_atom,
                all_bonds_explicit,
                all_hydrogens_explicit,
                include_dative_bonds,
                ignore_atom_map_numbers,
            );
            self._to_smiles_params(py, Some(&params))
        }
        /// Serialize SMILES with the selected canonicalization, stereochemistry and atom-order options. Supply SmilesWriteParams or keywords, not both. Uses the explicitly supplied parameter form.
        fn to_smiles_with_params(
            &self,
            py: Python<'_>,
            params: &SmilesWriteParams,
        ) -> PyResult<String> {
            self._to_smiles_params(py, Some(params))
        }

        #[pyo3(signature=(params=Argument::Missing, *, smiles=Argument::Missing, fields=Argument::Missing, coordinate_selection=Argument::Missing))]
        /// Serialize SMILES with the selected CX extensions. Supply CxSmilesWriteParams or keyword options, not both.
        fn to_cx_smiles(
            &self,
            py: Python<'_>,
            params: Argument<Option<PyRef<'_, CxSmilesWriteParams>>>,
            smiles: Argument<Option<PyRef<'_, SmilesWriteParams>>>,
            fields: Argument<Option<super::notation_params::CxSmilesFieldsInput>>,
            coordinate_selection: Argument<
                Option<super::notation_params::CxCoordinateSelectionInput>,
            >,
        ) -> PyResult<String> {
            crate::arguments::configuration_arguments!(params; smiles, fields, coordinate_selection);
            let params = params.as_deref();
            let smiles = smiles.as_deref();
            let configured = smiles.is_some() || fields.is_some() || coordinate_selection.is_some();
            if configured {
                self._to_cx_smiles_keywords(py, smiles, fields, coordinate_selection)
            } else {
                self._to_cx_smiles_params(py, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_cx_smiles"}"#]
        #[pyo3(signature=(params=None))]
        fn _to_cx_smiles_params(
            &self,
            py: Python<'_>,
            params: Option<&CxSmilesWriteParams>,
        ) -> PyResult<String> {
            let params = params.map(|value| value.to_rust(py)).transpose()?;
            let value = py.detach(|| match params {
                Some(ref params) => self.inner.to_cx_smiles_with_params(params),
                None => self.inner.to_cx_smiles(),
            });
            super::smiles_writer::text(py, value)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_cx_smiles"}"#]
        #[pyo3(signature=(*, smiles=None, fields=None, coordinate_selection=None))]
        fn _to_cx_smiles_keywords(
            &self,
            py: Python<'_>,
            smiles: Option<&SmilesWriteParams>,
            fields: Option<super::notation_params::CxSmilesFieldsInput>,
            coordinate_selection: Option<super::notation_params::CxCoordinateSelectionInput>,
        ) -> PyResult<String> {
            let smiles =
                pyo3::Python::attach(|py| smiles.map(|value| value.to_rust(py)).transpose())?;

            let defaults = ck::CxSmilesWriteParams::default();
            let params = ck::CxSmilesWriteParams {
                smiles: smiles.unwrap_or(defaults.smiles),
                fields: fields
                    .map(super::notation_params::CxSmilesFieldsInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.fields),
                coordinate_selection: coordinate_selection
                    .map(super::notation_params::CxCoordinateSelectionInput::into_rust)
                    .transpose()?
                    .unwrap_or(defaults.coordinate_selection),
            };
            super::smiles_writer::text(
                py,
                py.detach(|| self.inner.to_cx_smiles_with_params(&params)),
            )
        }
        /// Serialize SMILES with the selected CX extensions. Supply CxSmilesWriteParams or keyword options, not both. Uses the explicitly supplied parameter form.
        fn to_cx_smiles_with_params(
            &self,
            py: Python<'_>,
            params: &CxSmilesWriteParams,
        ) -> PyResult<String> {
            self._to_cx_smiles_params(py, Some(params))
        }

        #[pyo3(signature=(atoms=Argument::Missing, *, params=Argument::Missing, smiles=Argument::Missing, bonds=Argument::Missing, atom_symbols=Argument::Missing, bond_symbols=Argument::Missing))]
        /// Serialize only the selected atoms and bonds as fragment SMILES, without changing the molecule.
        fn to_fragment_smiles(
            &self,
            py: Python<'_>,
            atoms: Argument<Option<super::smiles_writer::FragmentSmilesInput>>,
            params: Argument<Option<PyRef<'_, FragmentSmilesWriteParams>>>,
            smiles: Argument<Option<PyRef<'_, SmilesWriteParams>>>,
            bonds: Argument<Option<Vec<usize>>>,
            atom_symbols: Argument<Option<Vec<String>>>,
            bond_symbols: Argument<Option<Vec<String>>>,
        ) -> PyResult<String> {
            if matches!(
                &atoms,
                Argument::Value(Some(super::smiles_writer::FragmentSmilesInput::Params(_)))
            ) && (params.is_supplied()
                || smiles.is_supplied()
                || bonds.is_supplied()
                || atom_symbols.is_supplied()
                || bond_symbols.is_supplied())
            {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "parameter object and configuration keywords are mutually exclusive",
                ));
            }
            crate::arguments::configuration_arguments!(params; atoms, smiles, bonds, atom_symbols, bond_symbols);
            let params = params.as_deref();
            let smiles = smiles.as_deref();
            use super::smiles_writer::FragmentSmilesInput;
            let configured = smiles.is_some()
                || bonds.is_some()
                || atom_symbols.is_some()
                || bond_symbols.is_some();
            match (atoms, params) {
                (Some(FragmentSmilesInput::Params(value)), None) if !configured => {
                    self._to_fragment_smiles_params(py, &*value.try_borrow(py)?)
                }
                (None, Some(params)) if !configured => self._to_fragment_smiles_params(py, params),
                (Some(FragmentSmilesInput::Atoms(atoms)), None) => self
                    ._to_fragment_smiles_keywords(
                        py,
                        atoms,
                        smiles,
                        bonds,
                        atom_symbols,
                        bond_symbols,
                    ),
                (None, None) => Err(pyo3::exceptions::PyTypeError::new_err(
                    "atoms or Params are required",
                )),
                _ => Err(pyo3::exceptions::PyTypeError::new_err(
                    "Params and atom/keyword configuration are mutually exclusive",
                )),
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_fragment_smiles"}"#]
        fn _to_fragment_smiles_params(
            &self,
            py: Python<'_>,
            params: &FragmentSmilesWriteParams,
        ) -> PyResult<String> {
            let params = params.to_rust(py)?;
            super::smiles_writer::text(
                py,
                py.detach(|| self.inner.to_fragment_smiles_with_params(&params)),
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_fragment_smiles"}"#]
        #[pyo3(signature=(atoms, *, smiles=None, bonds=None, atom_symbols=None, bond_symbols=None))]
        fn _to_fragment_smiles_keywords(
            &self,
            py: Python<'_>,
            atoms: Vec<usize>,
            smiles: Option<&SmilesWriteParams>,
            bonds: Option<Vec<usize>>,
            atom_symbols: Option<Vec<String>>,
            bond_symbols: Option<Vec<String>>,
        ) -> PyResult<String> {
            let smiles =
                pyo3::Python::attach(|py| smiles.map(|value| value.to_rust(py)).transpose())?;

            let atoms = atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>();
            if smiles.is_none()
                && bonds.is_none()
                && atom_symbols.is_none()
                && bond_symbols.is_none()
            {
                return super::smiles_writer::text(
                    py,
                    py.detach(|| self.inner.to_fragment_smiles(&atoms)),
                );
            }
            let params = ck::FragmentSmilesWriteParams {
                atoms,
                smiles: smiles.unwrap_or_else(|| ck::FragmentSmilesWriteParams::default().smiles),
                bonds: bonds.map(|values| values.into_iter().map(ck::BondId::new).collect()),
                atom_symbols,
                bond_symbols,
            };
            super::smiles_writer::text(
                py,
                py.detach(|| self.inner.to_fragment_smiles_with_params(&params)),
            )
        }
        /// Serialize only the selected atoms and bonds as fragment SMILES, without changing the molecule. Uses the explicitly supplied parameter form.
        fn to_fragment_smiles_with_params(
            &self,
            py: Python<'_>,
            params: &FragmentSmilesWriteParams,
        ) -> PyResult<String> {
            self._to_fragment_smiles_params(py, params)
        }

        #[pyo3(signature=(atoms=Argument::Missing, *, params=Argument::Missing, cx=Argument::Missing, bonds=Argument::Missing, atom_symbols=Argument::Missing, bond_symbols=Argument::Missing))]
        /// Serialize selected atoms and bonds as fragment CXSMILES with the requested extensions.
        fn to_fragment_cx_smiles(
            &self,
            py: Python<'_>,
            atoms: Argument<Option<super::smiles_writer::FragmentCxSmilesInput>>,
            params: Argument<Option<PyRef<'_, FragmentCxSmilesWriteParams>>>,
            cx: Argument<Option<PyRef<'_, CxSmilesWriteParams>>>,
            bonds: Argument<Option<Vec<usize>>>,
            atom_symbols: Argument<Option<Vec<String>>>,
            bond_symbols: Argument<Option<Vec<String>>>,
        ) -> PyResult<String> {
            if matches!(
                &atoms,
                Argument::Value(Some(super::smiles_writer::FragmentCxSmilesInput::Params(_)))
            ) && (params.is_supplied()
                || cx.is_supplied()
                || bonds.is_supplied()
                || atom_symbols.is_supplied()
                || bond_symbols.is_supplied())
            {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "parameter object and configuration keywords are mutually exclusive",
                ));
            }
            crate::arguments::configuration_arguments!(params; atoms, cx, bonds, atom_symbols, bond_symbols);
            let params = params.as_deref();
            let cx = cx.as_deref();
            use super::smiles_writer::FragmentCxSmilesInput;
            let configured =
                cx.is_some() || bonds.is_some() || atom_symbols.is_some() || bond_symbols.is_some();
            match (atoms, params) {
                (Some(FragmentCxSmilesInput::Params(value)), None) if !configured => {
                    self._to_fragment_cx_smiles_params(py, &*value.try_borrow(py)?)
                }
                (None, Some(params)) if !configured => {
                    self._to_fragment_cx_smiles_params(py, params)
                }
                (Some(FragmentCxSmilesInput::Atoms(atoms)), None) => self
                    ._to_fragment_cx_smiles_keywords(
                        py,
                        atoms,
                        cx,
                        bonds,
                        atom_symbols,
                        bond_symbols,
                    ),
                (None, None) => Err(pyo3::exceptions::PyTypeError::new_err(
                    "atoms or Params are required",
                )),
                _ => Err(pyo3::exceptions::PyTypeError::new_err(
                    "Params and atom/keyword configuration are mutually exclusive",
                )),
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_fragment_cx_smiles"}"#]
        fn _to_fragment_cx_smiles_params(
            &self,
            py: Python<'_>,
            params: &FragmentCxSmilesWriteParams,
        ) -> PyResult<String> {
            let params = params.to_rust(py)?;
            super::smiles_writer::text(
                py,
                py.detach(|| self.inner.to_fragment_cx_smiles_with_params(&params)),
            )
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_fragment_cx_smiles"}"#]
        #[pyo3(signature=(atoms, *, cx=None, bonds=None, atom_symbols=None, bond_symbols=None))]
        fn _to_fragment_cx_smiles_keywords(
            &self,
            py: Python<'_>,
            atoms: Vec<usize>,
            cx: Option<&CxSmilesWriteParams>,
            bonds: Option<Vec<usize>>,
            atom_symbols: Option<Vec<String>>,
            bond_symbols: Option<Vec<String>>,
        ) -> PyResult<String> {
            let atoms = atoms.into_iter().map(ck::AtomId::new).collect::<Vec<_>>();
            if cx.is_none() && bonds.is_none() && atom_symbols.is_none() && bond_symbols.is_none() {
                return super::smiles_writer::text(
                    py,
                    py.detach(|| self.inner.to_fragment_cx_smiles(&atoms)),
                );
            }
            let params = ck::FragmentCxSmilesWriteParams {
                atoms,
                cx: cx
                    .map(|value| value.to_rust(py))
                    .transpose()?
                    .unwrap_or_else(|| ck::FragmentCxSmilesWriteParams::default().cx),
                bonds: bonds.map(|values| values.into_iter().map(ck::BondId::new).collect()),
                atom_symbols,
                bond_symbols,
            };
            super::smiles_writer::text(
                py,
                py.detach(|| self.inner.to_fragment_cx_smiles_with_params(&params)),
            )
        }
        /// Serialize selected atoms and bonds as fragment CXSMILES with the requested extensions. Uses the explicitly supplied parameter form.
        fn to_fragment_cx_smiles_with_params(
            &self,
            py: Python<'_>,
            params: &FragmentCxSmilesWriteParams,
        ) -> PyResult<String> {
            self._to_fragment_cx_smiles_params(py, params)
        }

        #[pyo3(signature=(count, seed, params=Argument::Missing, *, isomeric_smiles=Argument::Missing, kekule=Argument::Missing, all_bonds_explicit=Argument::Missing, all_hydrogens_explicit=Argument::Missing))]
        /// Generate randomized SMILES using the requested count and random seed; molecular chemistry is unchanged.
        fn to_random_smiles(
            &self,
            py: Python<'_>,
            count: u32,
            seed: u32,
            params: Argument<Option<PyRef<'_, RandomSmilesWriteParams>>>,
            isomeric_smiles: Argument<Option<bool>>,
            kekule: Argument<Option<bool>>,
            all_bonds_explicit: Argument<Option<bool>>,
            all_hydrogens_explicit: Argument<Option<bool>>,
        ) -> PyResult<Vec<String>> {
            crate::arguments::configuration_arguments!(params; isomeric_smiles, kekule, all_bonds_explicit, all_hydrogens_explicit);
            let params = params.as_deref();
            let configured = isomeric_smiles.is_some()
                || kekule.is_some()
                || all_bonds_explicit.is_some()
                || all_hydrogens_explicit.is_some();
            if configured {
                self._to_random_smiles_keywords(
                    py,
                    count,
                    seed,
                    isomeric_smiles,
                    kekule,
                    all_bonds_explicit,
                    all_hydrogens_explicit,
                )
            } else {
                self._to_random_smiles_params(py, count, seed, params)
            }
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_random_smiles"}"#]
        #[pyo3(signature=(count, seed, params=None))]
        fn _to_random_smiles_params(
            &self,
            py: Python<'_>,
            count: u32,
            seed: u32,
            params: Option<&RandomSmilesWriteParams>,
        ) -> PyResult<Vec<String>> {
            let params =
                pyo3::Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;

            let value = py.detach(|| match params {
                Some(params) => self
                    .inner
                    .to_random_smiles_with_params(count, seed, &params),
                None => self.inner.to_random_smiles(count, seed),
            });
            super::smiles_writer::texts(py, value)
        }
        #[doc = r#"{"role":"overload","target":"Molecule.to_random_smiles"}"#]
        #[pyo3(signature=(count, seed, *, isomeric_smiles=None, kekule=None, all_bonds_explicit=None, all_hydrogens_explicit=None))]
        fn _to_random_smiles_keywords(
            &self,
            py: Python<'_>,
            count: u32,
            seed: u32,
            isomeric_smiles: Option<bool>,
            kekule: Option<bool>,
            all_bonds_explicit: Option<bool>,
            all_hydrogens_explicit: Option<bool>,
        ) -> PyResult<Vec<String>> {
            let defaults = ck::RandomSmilesWriteParams::default();
            let params = RandomSmilesWriteParams {
                inner: crate::parameters::NativeField::Owned(ck::RandomSmilesWriteParams {
                    isomeric_smiles: isomeric_smiles.unwrap_or(defaults.isomeric_smiles),
                    kekule: kekule.unwrap_or(defaults.kekule),
                    all_bonds_explicit: all_bonds_explicit.unwrap_or(defaults.all_bonds_explicit),
                    all_hydrogens_explicit: all_hydrogens_explicit
                        .unwrap_or(defaults.all_hydrogens_explicit),
                }),
            };
            self._to_random_smiles_params(py, count, seed, Some(&params))
        }
        /// Generate randomized SMILES using the requested count and random seed; molecular chemistry is unchanged. Uses the explicitly supplied parameter form.
        fn to_random_smiles_with_params(
            &self,
            py: Python<'_>,
            count: u32,
            seed: u32,
            params: &RandomSmilesWriteParams,
        ) -> PyResult<Vec<String>> {
            self._to_random_smiles_params(py, count, seed, Some(params))
        }

        /// Read the selected bond's text property, returning None for an absent key and a structured error for an incompatible value.
        fn bond_property_string(
            &self,
            py: Python<'_>,
            id: usize,
            key: &str,
        ) -> PyResult<Option<String>> {
            self.inner
                .bond_property_string(ck::BondId::new(id), key)
                .map_err(crate::canonical_errors::PropertyStringError::from_rust)?
                .map(|text| crate::property_values::decode_source_text(py, &text))
                .transpose()
        }
    }

    #[pymodule_export]
    pub use crate::fingerprint_values::Fingerprint;

    #[pymethods]
    impl Fingerprint {
        #[staticmethod]
        /// Construct a bit fingerprint of n_bits positions from the supplied set-bit indices; invalid indices raise FingerprintError.
        fn from_on_bits(n_bits: u32, on_bits: Vec<u32>) -> PyResult<Self> {
            ck::Fingerprint::from_on_bits(n_bits, on_bits)
                .map(|inner| Self { inner })
                .map_err(errors::FingerprintError::from_rust)
        }
        #[getter]
        /// Size of the fingerprint's bit index space, not the number of set bits.
        fn n_bits(&self) -> u32 {
            self.inner.n_bits()
        }
        /// Calculate Tanimoto similarity to another fingerprint using the native size and empty-fingerprint rules.
        fn tanimoto(&self, other: &Self) -> PyResult<f64> {
            self.inner
                .tanimoto(&other.inner)
                .map_err(errors::FingerprintError::from_rust)
        }
        #[getter]
        /// Return set-bit indices in increasing order.
        fn on_bits(&self) -> Vec<u32> {
            self.inner.on_bits()
        }
        /// Read the bit at a zero-based index; out-of-range indices raise FingerprintError.
        fn get_bit(&self, bit: u32) -> PyResult<bool> {
            self.inner
                .get_bit(bit)
                .map_err(errors::FingerprintError::from_rust)
        }
        #[doc = r#"{"role": "return_type", "imports": ["numpy"], "type": "numpy.ndarray[tuple[int], numpy.dtype[numpy.uint8]]"}"#]
        /// Return an independent uint8 array with one entry per fingerprint bit.
        /// Shape is (n_bits,); entries are 0 or 1 in increasing bit-index order.
        fn to_numpy<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, numpy::PyArray1<u8>>> {
            let values = (0..self.inner.n_bits())
                .map(|bit| self.inner.get_bit(bit).map(u8::from))
                .collect::<Result<Vec<_>, _>>()
                .map_err(errors::FingerprintError::from_rust)?;
            Ok(values.into_pyarray(py))
        }
    }

    #[pymodule_export]
    use super::fingerprint_batch::{
        FingerprintBatch, FingerprintBatchError, FingerprintBatchErrorKind,
    };

    // Register SMILES configuration calls and the molecule collection.
    #[pymodule_export]
    use super::molecular_io_errors::{
        AvalonFingerprintError, AvalonFingerprintErrorKind, MolWriteError, MolWriteErrorKind,
        MolecularIoError, MolecularIoErrorKind, SdfWriteError, SdfWriteErrorKind,
    };

    #[pymodule_export]
    use super::drawing_errors::{
        DoubleLexicalReadError, DoubleLexicalReadErrorKind, DrawingError, DrawingErrorKind,
        DrawingWriteError, DrawingWriteErrorKind, SourceNumericStreamState, SvgParseError,
        SvgParseErrorKind, XmlParseError, XmlParseErrorKind, XmlTextPos,
    };

    #[pymodule_export]
    use super::batch_errors::{
        BatchError, BatchExportReport, BatchImageError, BatchImageErrorKind, BatchValidationError,
    };

    #[pymodule_export]
    use super::sdf_sources::{
        SdfBatchIterator, SdfDataset, SdfDatasetIterator, SdfReader, SdfReaderBatchIterator,
        SdfRecordStream,
    };

    #[pymodule_export]
    use super::batch_calls::{
        _mols_from_sdf_records_keywords, _mols_from_sdf_records_params,
        _mols_from_smiles_list_keywords, _mols_from_smiles_list_params, mols_from_records,
        mols_from_sdf_records, mols_from_smiles_list,
    };

    #[pymodule_export]
    use super::molecular_io_calls::{
        _mol_from_mmcif_keywords, _mol_from_mmcif_params, _mol_from_mol2_keywords,
        _mol_from_mol2_params, _mol_from_pdb_keywords, _mol_from_pdb_params, mol_from_mmcif,
        mol_from_mol2, mol_from_pdb, mol_from_xyz_block,
    };

    #[pymodule_export]
    use super::extension::{MoleculeBatch, SmilesParseParams, mol_from_smiles};

    #[pymodule_export]
    use super::extension::{_mol_from_smiles_keywords, _mol_from_smiles_params};
}

mod extension;
