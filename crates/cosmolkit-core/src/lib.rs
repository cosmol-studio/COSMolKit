//! COSMolKit core built around value-style molecule state.
//!
//! This crate is the active molecular graph, state, operation, IO, and chemistry
//! core. Reference-only crates must not be used as dependencies for this crate.
//!
//! # Non-negotiable architecture rules
//!
//! - `Molecule` is a value object. Public transforms return a new `Molecule`.
//! - Public APIs must not expose `atoms_mut`, `bonds_mut`, `topology_mut`, or
//!   direct mutable access to internal storage.
//! - Construction happens through `MoleculeBuilder`.
//! - Mutation of existing molecules happens through registered operations.
//! - Derived state must be invalidated or recomputed through an operation
//!   contract, not by ad hoc field edits.
//! - RDKit compatibility belongs in a future `compat::rdkit` layer. It must not
//!   define the canonical shape of core molecule state.
//!
//! Agent guardrail: if a future change requires bypassing any rule above, the
//! agent is not allowed to "just make it work". The agent must stop and confirm
//! the design exception with the human author before editing code that violates
//! these standards.
//!
//! # Examples
//!
//! Value-style operations return a new molecule and leave the receiver
//! unchanged:
//!
//! ```
//! use cosmolkit_core::Molecule;
//!
//! let mol = Molecule::from_smiles("CCO").unwrap();
//! let with_hydrogens = mol.with_hydrogens().unwrap();
//!
//! assert_eq!(mol.num_atoms(), 3);
//! assert!(with_hydrogens.num_atoms() > mol.num_atoms());
//! ```
//!
//! Public in-place molecule operations end with `_` and mutate the receiver:
//!
//! ```
//! use cosmolkit_core::Molecule;
//!
//! let mut mol = Molecule::from_smiles("CCO").unwrap();
//! let expected_atoms = mol.with_hydrogens().unwrap().num_atoms();
//!
//! mol.add_hydrogens_().unwrap();
//!
//! assert_eq!(mol.num_atoms(), expected_atoms);
//! ```
//!
//! Atom chiral tags can be assigned from a selected 3D conformer through a
//! registered, transactional operation with pinned-RDKit parity:
//!
//! ```
//! use cosmolkit_core::{ChiralTag, Molecule};
//!
//! let mol = Molecule::from_smiles("C(F)(Cl)Br")
//!     .unwrap()
//!     .with_only_3d_conformer(
//!         vec![
//!             [0.0, 0.0, 0.0],
//!             [1.0, 0.0, 0.0],
//!             [0.0, 1.0, 0.0],
//!             [0.0, 0.0, 1.0],
//!         ],
//!         true,
//!     )
//!     .unwrap();
//! let assigned = mol.with_chiral_tags_from_structure(-1, true).unwrap();
//!
//! assert_eq!(mol.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
//! assert_ne!(assigned.atoms()[0].chiral_tag(), ChiralTag::Unspecified);
//! ```
//!
//! Coordinate operations follow the same value-style and in-place split:
//!
//! ```
//! use cosmolkit_core::Molecule;
//!
//! let mol = Molecule::from_smiles("c1ccccc1").unwrap();
//! let with_coords = mol.with_2d_coordinates().unwrap();
//! assert!(mol.coordinates_2d().is_none());
//! assert_eq!(with_coords.coordinates_2d().unwrap().len(), with_coords.num_atoms());
//!
//! let mut in_place = Molecule::from_smiles("CCO").unwrap();
//! in_place.compute_2d_coordinates_().unwrap();
//! assert_eq!(in_place.coordinates_2d().unwrap().len(), in_place.num_atoms());
//! ```

pub mod bio;
pub mod chemistry;
pub mod confseq;
pub mod inchi;
pub mod io;
pub mod model;
pub mod notation;
pub mod operations;
pub mod properties;
pub mod search;
mod source_port_helpers;
pub mod support;
#[cfg(test)]
pub(crate) use cosmolkit_test_support as test_data;

pub use bio::invariants as bio_invariants;
pub use bio::ops as bio_ops;
pub use chemistry::{
    aromaticity, atropisomer, coordinates, distgeom, hydrogens, kekulize, mol_align,
    mol_transforms, rings, stereo, stereo_enumerate, tautomer, valence,
};
pub use io::pdb_writer;
pub(crate) use model::invariants;
pub use model::{adjacency, atom, bond, builder, derived, error, molecule, read_parts, sgroup};
pub use mol_align::{
    AlignmentAtomMap, AlignmentError, AlignmentParameters, AlignmentResult, AlignmentTransform,
    AllConformerRmsdParameters, BestAlignmentParameters, ConformerAlignmentParameters,
    ConformerAlignmentReport, ConformerRmsd, CoordinateRmsdParameters,
};
pub(crate) use notation::{canon_rank, smiles};
pub use notation::{canon_smiles, fragment, sequence, smiles_write};
pub use operations::{ops, sanitize};
pub use properties::avalon_fingerprint::{AvalonFingerprintFlags, AvalonFingerprintParams};
pub use properties::{avalon_fingerprint, batch, draw, fingerprint, mol_hash, mol_pickler};
pub use search::{
    SmartsParseParams, SmartsWriteError, get_atom_smarts, get_bond_smarts,
    mol_fragment_to_cx_smarts, mol_fragment_to_smarts, mol_from_smarts, mol_to_cx_smarts,
    mol_to_smarts,
};
pub use search::{query, smarts_parse, substruct};

pub use adjacency::{AdjacencyError, AdjacencyList, NeighborRef};
pub use aromaticity::{AromaticityAssignment, AromaticityError, AromaticityModel, set_aromaticity};
pub use atom::{
    Atom, AtomId, AtomPdbResidueInfo, AtomSpec, ChiralTag, ELEMENTS, ELEMENTS_WITH_DUMMY, Element,
    ElementInfo, ElementParseError, Hybridization,
};
pub use batch::{
    BatchErrorMode, BatchExportReport, BatchProgress, BatchProgressBar, BatchRecord,
    BatchRecordError, BatchValidationError, MoleculeBatch, batch_progress_bar,
};
pub use bio::protein::{
    Protein, ProteinAtomRef, ProteinChainRef, ProteinResidueRef, ProteinSelectionSummary,
};
pub use bio::resinfo::{
    RESIDUE_INFO_TABLE, ResidueCode, ResidueCodeParseError, ResidueIdentity, ResidueInfo,
    ResidueInfoKind, ResidueInfoSequenceError, UNKNOWN_TABULATED_RESIDUE_IDX, expand_one_letter,
    expand_one_letter_sequence, expand_protein_one_letter, expand_protein_one_letter_string,
    find_tabulated_residue, find_tabulated_residue_idx, get_residue_info, get_residue_info_checked,
    residue_code_from_name,
};
pub use bio::{
    AltLocLabel, AtomName, AtomRow, AtomSourceIds, BioAssembly, BioCisPep, BioConnection,
    BioConnectionType, BioCoorFormat, BioMetadata, BioModRes, BioNcsOperator, BioStructure,
    BioTransform, ChainKind, ChainRow, ChainSourceIds, CoordinateBlock, CrystalCell, CrystalInfo,
    EntityKind, EntityRow, EntitySourceIds, ModelRow, PolymerKind, ResidueKind, ResidueName,
    ResidueRow, RowSpan, classify_residue_name,
};
pub use bio::{
    AtomId as BioAtomId, ChainId as BioChainId, EntityId as BioEntityId, ModelId as BioModelId,
    ResidueId as BioResidueId,
};
pub use bio_ops::{
    BioBlockSet, BioDerivedState, BioEditKind, BioOpDomain, BioOpKind, BioOperationError,
    BioParityPolicy, BioRowMapping, BioStateSet, BioStructureMapping, BioStructureOpSpec,
};
pub use bond::{Bond, BondDirection, BondId, BondOrder, BondSpec, BondStereo};
pub use builder::MoleculeBuilder;
pub use chemistry::cip::{CipDescriptor, CipDescriptorError, CipLabelOptions, CipLabelerError};
pub use chemistry::forcefield::mmff::{
    MMFF_MOL_PROPERTIES_FEATURE, MmffAngle, MmffAngleCollection, MmffAtomProperties, MmffBond,
    MmffBondCollection, MmffChg, MmffChgCollection, MmffDef, MmffDefCollection, MmffDfsbCollection,
    MmffMolProperties, MmffMolPropertiesError, MmffOop, MmffOopCollection,
    MmffOptimizeMoleculeConfResult, MmffOptimizeMoleculeConfsResult, MmffOptimizeMoleculeResult,
    MmffParamError, MmffPbci, MmffPbciCollection, MmffProp, MmffPropCollection, MmffPublicApiError,
    MmffStbn, MmffStbnCollection, MmffTor, MmffTorCollection, MmffVariant, MmffVdw,
    MmffVdwCollection, MmffVdwRijstarEps, mmff_has_all_molecule_params,
    mmff_initial_gradient_for_parity, mmff_optimize_molecule, mmff_optimize_molecule_confs,
    mmff_sanitize_ops, sanitize_mmff_mol,
};
pub use chemistry::forcefield::uff::{
    UffAngle, UffBond, UffInv, UffOptimizeMoleculeConfResult, UffOptimizeMoleculeConfsResult,
    UffOptimizeMoleculeResult, UffPublicApiError, UffTor, UffVdw, get_uff_angle_bend_params,
    get_uff_bond_stretch_params, get_uff_inversion_params, get_uff_torsion_params,
    get_uff_vdw_params, uff_has_all_molecule_params, uff_initial_gradient_for_parity,
    uff_optimize_molecule, uff_optimize_molecule_confs,
};
pub use chemistry::forcefield::{
    AngleConstraintContrib, AngleConstraintContribs, AngleConstraintContribsParams, DihedralOutput,
    DistanceConstraintContrib, DistanceConstraintContribs, DistanceConstraintContribsParams,
    ForceField, ForceFieldContrib, ForceFieldSnapshot, ForceFieldVec3, PositionConstraintContrib,
    TorsionAngleContribM6, TorsionAngleContribs, TorsionAngleContribsParams,
    TorsionConstraintContrib, calc_torsion_energy, calc_torsion_energy_m6,
    compute_dihedral_from_flat, compute_dihedral_from_points, compute_dihedral_from_position_vec,
    normalize_angle_deg,
};
pub use chemistry::tautomer::{
    TautomerCanonicalizationError, TautomerCatalog, TautomerCatalogError, TautomerEnumeration,
    TautomerEnumerationCallback, TautomerEnumerationError, TautomerEnumerationStatus,
    TautomerEnumerator, TautomerOptions, TautomerRunError, TautomerScore, TautomerScoreError,
    TautomerScoreTerm, TautomerTransform, TautomerTransformError, default_tautomer_score_terms,
    score_tautomer, score_tautomer_hetero_hydrogens, score_tautomer_rings,
    score_tautomer_substructures, score_tautomer_with_terms,
};
pub use confseq::{
    ConfSeqBatchDecodeResult, ConfSeqDecodeError, ConfSeqDecodeOptions, ConfSeqFastGeometryError,
    ConfSeqTemplateBackend, decode_confseq, decode_confseq_batch,
    decode_confseq_batch_with_options, decode_confseq_record, decode_confseq_record_batch,
    decode_confseq_record_batch_with_options, decode_confseq_record_with_options,
    decode_confseq_with_options,
};
pub use coordinates::With2DCoordinatesParams;
pub use derived::DerivedState;
pub use distgeom::{
    DgBoundsError, EmbedFailureCause, EmbedMoleculeResult, EmbedMultipleConfsResult,
    EmbedParameters, embed_molecule, embed_molecule_result, embed_multiple_confs,
    embed_multiple_confs_result, embed_multiple_confs_return_vector,
};
pub use draw::SvgDrawError;
pub use error::{InvariantError, MoleculeBuildError};
pub use fingerprint::{
    ATOM_PAIR_ATOM_NUMBER_TYPES, ATOM_PAIR_CODE_SIZE, ATOM_PAIR_MAX_NUM_BRANCHES,
    ATOM_PAIR_MAX_NUM_PI, ATOM_PAIR_MAX_PATH_LENGTH, ATOM_PAIR_NUM_BRANCH_BITS,
    ATOM_PAIR_NUM_CHIRAL_BITS, ATOM_PAIR_NUM_FINGERPRINT_BITS, ATOM_PAIR_NUM_PATH_BITS,
    ATOM_PAIR_NUM_PI_BITS, ATOM_PAIR_NUM_TYPE_BITS, ATOM_PAIRS_VERSION, AdditionalOutput,
    AtomPairFingerprintGenerator, AtomPairFingerprintOutput, AtomPairFingerprintParams,
    Fingerprint, FingerprintArguments, FingerprintError, FingerprintFuncArguments,
    LAYERED_FINGERPRINT_MAX_LAYERS, LAYERED_FINGERPRINT_SUBSTRUCTURE_LAYERS,
    LAYERED_FINGERPRINT_VERSION, LayeredFingerprintLayers, LayeredFingerprintParams,
    LayeredFingerprintResult, MaccsFingerprintParams, MorganAdditionalOutput,
    MorganAtomInvariantsGenerator, MorganBitFingerprintOutput, MorganBondInvariantsGenerator,
    MorganFingerprintOutput, MorganFingerprintParams, MorganSparseFingerprintOutput,
    PATTERN_FINGERPRINT_VERSION, PatternFingerprintParams, SparseBitFingerprint,
    SparseCountFingerprint, TopologicalFingerprintOutput, TopologicalFingerprintOutputRequest,
    TopologicalFingerprintParams, TopologicalFingerprintResult,
    TopologicalTorsionFingerprintGenerator, TopologicalTorsionFingerprintOutputRequest,
    TopologicalTorsionFingerprintParams, TopologicalTorsionFingerprintResult,
    TopologicalTorsionFingerprintValue, TopologicalTorsionFingerprintVector,
    TopologicalTorsionLegacyKind, TopologicalTorsionLegacyParams, TopologicalTorsionLegacyResult,
    atom_pair_fingerprint, atom_pair_fingerprint_with_output, get_atom_code,
    get_topological_torsion_code, get_topological_torsion_hash, layered_fingerprint,
    layered_fingerprint_with_output, maccs_fingerprint, maccs_get_fingerprint_as_bit_vect,
    morgan_get_fingerprint, morgan_get_fingerprint_as_bit_vect, morgan_get_hashed_fingerprint,
    pattern_fingerprint, topological_fingerprint_with_output,
    topological_torsion_count_fingerprint, topological_torsion_fingerprint,
    topological_torsion_fingerprint_with_output, topological_torsion_generator,
    topological_torsion_legacy_fingerprint, topological_torsion_sparse_count_fingerprint,
    topological_torsion_sparse_fingerprint,
};
pub use hydrogens::{AddHsParams, AddHydrogensError, RemoveHsParams, RemoveHydrogensError};
pub use inchi::{
    INCHI_API_PARITY_MATRIX, InchiApiParitySpec, InchiDiagnostic, InchiDiagnosticLevel, InchiError,
    InchiErrorKind, InchiReturnValues, InchiToInchiKeyOutput, MolFromInchiOutput,
    MolToInchiKeyOutput, MolToInchiOutput, inchi_to_inchi_key, mol_from_inchi, mol_to_inchi,
    mol_to_inchi_key,
};
#[allow(deprecated)]
pub use io::bio::read_mmcif_atom_site_subset_from_str;
pub use io::bio::{
    BioPdbReadParams, BioReadError, BioWriteError, MmcifOutputGroups, MmcifWriteOptions,
};
pub use io::mol2::{
    Mol2ReadError, Mol2ReadParams, Mol2Record, Mol2Type, mol_from_mol2_block_like_rdkit,
    mol_from_mol2_data_stream_like_rdkit, mol_from_mol2_file_like_rdkit, read_mol2_file,
    read_mol2_file_with_params, read_mol2_from_str, read_mol2_from_str_with_params,
};
#[allow(deprecated)]
pub use io::pdb_molecule::{
    PdbMoleculeConversionError, RdkitPdbMolProfile, bio_structure_to_rdkit_pdb_molecule,
    molecule_from_mmcif_block_with_options, molecule_from_pdb_block_with_options,
};
pub use io::pdb_molecule::{StructureMoleculeConversionError, StructureMoleculeOptions};
pub use io::sdf::{SdfCoordinateMode, SdfDataset, SdfReadParams, SdfRecordMetadata};
pub use io::xyz::{XyzReadError, read_xyz_from_str};
pub use kekulize::KekulizeError;
pub use mol_pickler::{PickleError, mol_from_binary, mol_to_binary};
pub use molecule::{
    AtomMapping, BondMapping, Conformer2D, Conformer3D, ConformerStore, CoordinateDimension,
    Molecule, MoleculeCapabilities, MoleculeProperties, PropertyStore, SdfPropertyList,
    SdfPropertyListTarget, SmilesParseError, SmilesWriteError, TopologyMapping, TopologyTrust,
};
pub use notation::smiles::SmilesParseParams;
pub use ops::{
    ASSIGNED_AROMATICITY_SPEC, ASSIGNED_RING_FAMILIES_SPEC, ASSIGNED_RINGS_SPEC,
    ASSIGNED_VALENCE_SPEC, BlockAccess, BlockSet, CipStatePolicy,
    ENUMERATE_TAUTOMERS_WITH_OPTIONS_SPEC, InvariantCheckSet, MOLECULE_OPS, MappingRequirement,
    MoleculeOpKind, MoleculeOpOutput, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OperationDomain,
    OperationError, OperationInvariantEntry, OperationTrace, PARITY_MATRIX, ParityMatrixEntry,
    ParityPolicy, SANITIZE_SPEC, SUPPORT_MATRIX, SemanticPrecondition, SemanticPreconditionSet,
    SupportMatrixEntry, TopologyEditKind, WITH_2D_COORDINATES_SPEC, WITH_3D_CONFORMER_SPEC,
    WITH_3D_CONFORMERS_SPEC, WITH_ADDED_3D_CONFORMER_SPEC, WITH_CLEARED_3D_CONFORMERS_SPEC,
    WITH_HYDROGENS_SPEC, WITH_KEKULIZED_BONDS_SPEC, WITH_ONLY_3D_CONFORMER_SPEC,
    WITHOUT_HYDROGENS_SPEC, WITHOUT_HYDROGENS_WITH_PARAMS_SPEC,
};
pub use properties::descriptors::{
    CrippenDescriptorValues, DescriptorError, DescriptorResult, NumRotatableBondsOptions,
    calc_crippen_descriptors, calc_exact_mol_wt, calc_fraction_csp3, calc_mol_formula, calc_mol_wt,
    calc_num_aromatic_rings, calc_num_hba, calc_num_hbd, calc_num_rotatable_bonds, calc_qed,
    calc_tpsa,
};
pub use query::{AtomQueryPredicate, BondQueryPredicate, QueryNode, SmartsParseError};
pub use rings::{
    RingFindType, RingFindingError, RingInfo, fast_find_rings, find_ring_families, find_sssr,
    symmetrize_sssr,
};
pub use sanitize::{SanitizeError, SanitizeOps, SanitizeStep, detect_chemistry_problems};
pub use sgroup::{
    SGroupAttachPoint, SGroupBondRole, SGroupBracket, SGroupBracketStyle, SGroupCState,
    SGroupConnection, SGroupData, SGroupDisplay, SubstanceGroup, SubstanceGroupId,
    SubstanceGroupKind,
};
pub use smiles_write::{
    CxSmilesFields, RestoreBondDirOption, SmilesWriteParams, mol_to_random_smiles_vect,
};
pub use stereo::{
    DoubleBondStereo, LigandRef, StereoError, StereoGroup, StereoGroupKind, TetrahedralStereo,
    assign_stereochemistry, perceive_stereochemistry,
};
pub use substruct::{
    AtomCoordsMatchFunctor, ExtraAtomCheck, ExtraBondCheck, ExtraFinalCheck, SubstructMatchError,
    SubstructMatchOverload, SubstructMatchParams, SubstructMatchParamsJsonError,
    SubstructMatchResult, check_substruct_match_overload_support, get_substruct_match,
    get_substruct_matches, get_substruct_matches_with_params, has_substruct_match,
    substruct_match_params_to_json, try_get_substruct_matches_with_params,
    update_substruct_match_params_from_json,
};
pub use support::{
    AROMATICITY_FEATURE, ATOM_PAIR_FINGERPRINT_FEATURE, AVALON_FINGERPRINT_FEATURE, BATCH_FEATURE,
    BIO_MMCIF_READ_FEATURE, BIO_MMCIF_WRITE_FEATURE, BIO_PDB_READ_FEATURE, BIO_SELECTION_FEATURE,
    BIO_STRUCTURE_FEATURE, CIP_LABELER_FEATURE, CONFORMER_GENERATION_FEATURE,
    COORDINATE_2D_FEATURE, COORDINATE_EDIT_FEATURE, DESCRIPTORS_FEATURE, DG_BOUNDS_FEATURE,
    DRAWING_FEATURE, FINGERPRINT_FEATURE, FeatureCategory, FeatureSpec, HYDROGENS_FEATURE,
    INCHI_FEATURE, KEKULIZE_FEATURE, LAYERED_FINGERPRINT_FEATURE, MOLALIGN_FEATURE,
    MOLBLOCK_IO_FEATURE, PATTERN_FINGERPRINT_FEATURE, PUBLIC_FEATURES, RINGS_FEATURE,
    SANITIZE_FEATURE, SMILES_PARSE_FEATURE, SMILES_WRITE_FEATURE, STEREO_FEATURE,
    SUBSTRUCTURE_FEATURE, SupportStatus, TAUTOMER_ENUMERATION_FEATURE,
    TOPOLOGICAL_TORSION_FINGERPRINT_FEATURE, UnsupportedFeatureError, VALENCE_FEATURE,
};
#[allow(deprecated)]
pub use support::{
    BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE, BIO_PDB_COORDINATE_SUBSET_READ_FEATURE,
};
pub use valence::{
    ValenceAssignment, ValenceError, ValenceModel, assign_radicals, assign_valence,
    assign_valence_with_options, atom_has_valence_violation, cached_valence_assignment,
    rdkit_atomic_number_from_symbol, rdkit_element_symbol, rdkit_valence_list,
};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
