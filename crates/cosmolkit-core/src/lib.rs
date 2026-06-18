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
pub mod io;
pub mod model;
pub mod notation;
pub mod operations;
pub mod properties;
pub mod search;
pub mod support;

pub use bio::invariants as bio_invariants;
pub use bio::ops as bio_ops;
pub use chemistry::{
    aromaticity, atropisomer, coordinates, distgeom, hydrogens, kekulize, mol_transforms, rings,
    stereo, stereo_enumerate, valence,
};
pub use io::pdb_writer;
pub(crate) use model::invariants;
pub use model::{adjacency, atom, bond, builder, derived, error, molecule, read_parts, sgroup};
pub(crate) use notation::{canon_rank, smiles};
pub use notation::{canon_smiles, fragment, sequence, smiles_write};
pub use operations::{ops, sanitize};
pub use properties::{avalon_fingerprint, batch, draw, fingerprint, mol_hash, mol_pickler};
pub use search::{query, smarts_parse, substruct};

pub use adjacency::{AdjacencyError, AdjacencyList, NeighborRef};
pub use aromaticity::{AromaticityAssignment, AromaticityError, AromaticityModel, set_aromaticity};
pub use atom::{Atom, AtomId, AtomPdbResidueInfo, AtomSpec, ChiralTag, Element, Hybridization};
pub use batch::{
    BatchErrorMode, BatchExportReport, BatchProgress, BatchProgressBar, BatchRecord,
    BatchRecordError, BatchValidationError, MoleculeBatch, batch_progress_bar,
};
pub use bio::protein::{
    Protein, ProteinAtomRef, ProteinChainRef, ProteinResidueRef, ProteinSelectionSummary,
};
pub use bio::resinfo::{
    RESIDUE_INFO_TABLE, ResidueCode, ResidueInfo, ResidueInfoKind, ResidueInfoSequenceError,
    UNKNOWN_TABULATED_RESIDUE_IDX, expand_one_letter, expand_one_letter_sequence,
    expand_protein_one_letter, expand_protein_one_letter_string, find_tabulated_residue,
    find_tabulated_residue_idx, get_residue_info, get_residue_info_checked, residue_code_from_name,
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
    BioBlockSet, BioDerivedState, BioEditKind, BioOpDomain, BioOpKind, BioOpOutcome,
    BioOperationError, BioParityPolicy, BioRowMapping, BioStateSet, BioStructureMapping,
    BioStructureOpSpec,
};
pub use bond::{Bond, BondDirection, BondId, BondOrder, BondSpec, BondStereo};
pub use builder::MoleculeBuilder;
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
pub use confseq::{
    ConfSeqBaseConformerError, ConfSeqBatchDecodeResult, ConfSeqDecodeError, ConfSeqDecodeOptions,
    ConfSeqTemplateBackend, decode_confseq, decode_confseq_batch,
    decode_confseq_batch_with_options, decode_confseq_with_options,
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
    Fingerprint, FingerprintError, MorganAdditionalOutput, MorganAtomInvariantsGenerator,
    MorganBondInvariantsGenerator, MorganFingerprintOutput, MorganFingerprintParams,
};
pub use hydrogens::{AddHsParams, AddHydrogensError, RemoveHsParams, RemoveHydrogensError};
pub use io::bio::{BioPdbReadParams, BioReadError, read_mmcif_atom_site_subset_from_str};
pub use io::mol2::{
    Mol2ReadError, Mol2ReadParams, Mol2Record, Mol2Type, mol_from_mol2_block_like_rdkit,
    mol_from_mol2_data_stream_like_rdkit, mol_from_mol2_file_like_rdkit, read_mol2_file,
    read_mol2_file_with_params, read_mol2_from_str, read_mol2_from_str_with_params,
};
pub use io::pdb_molecule::{
    PdbMoleculeConversionError, RdkitPdbMolProfile, bio_structure_to_rdkit_pdb_molecule,
    molecule_from_mmcif_block_with_options, molecule_from_pdb_block_with_options,
};
pub use io::sdf::{SdfCoordinateMode, SdfDataset, SdfReadParams, SdfRecordMetadata};
pub use io::xyz::{XyzReadError, read_xyz_from_str};
pub use kekulize::KekulizeError;
pub use mol_pickler::{PickleError, mol_from_binary, mol_to_binary};
pub use molecule::{
    AtomMapping, BondMapping, Conformer2D, Conformer3D, ConformerStore, CoordinateDimension,
    Molecule, MoleculeCapabilities, MoleculeProperties, PropertyStore, SdfPropertyList,
    SdfPropertyListTarget, SmilesParseError, SmilesWriteError, TopologyMapping, TopologyTrust,
};
pub use ops::{
    ASSIGNED_AROMATICITY_SPEC, ASSIGNED_RING_FAMILIES_SPEC, ASSIGNED_RINGS_SPEC,
    ASSIGNED_VALENCE_SPEC, BlockAccess, BlockSet, InvariantCheckSet, MOLECULE_OPS,
    MappingRequirement, MoleculeOpKind, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OpOutcome,
    OperationDomain, OperationError, OperationInvariantEntry, OperationTrace, PARITY_MATRIX,
    ParityMatrixEntry, ParityPolicy, SANITIZE_SPEC, SUPPORT_MATRIX, SemanticPrecondition,
    SemanticPreconditionSet, SupportMatrixEntry, TopologyEditKind, WITH_2D_COORDINATES_SPEC,
    WITH_3D_CONFORMER_SPEC, WITH_3D_CONFORMERS_SPEC, WITH_HYDROGENS_SPEC,
    WITH_KEKULIZED_BONDS_SPEC, WITHOUT_HYDROGENS_SPEC, WITHOUT_HYDROGENS_WITH_PARAMS_SPEC,
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
pub use smiles::assign_double_bond_stereo_from_directions;
pub use smiles_write::{
    CxSmilesFields, RestoreBondDirOption, SmilesWriteParams, mol_to_random_smiles_vect,
};
pub use stereo::{
    DoubleBondStereo, LigandRef, StereoError, StereoGroup, StereoGroupKind, TetrahedralStereo,
    assign_stereochemistry, perceive_stereochemistry,
};
pub use substruct::{
    SubstructMatchParams, SubstructMatchResult, get_substruct_match, get_substruct_matches,
    get_substruct_matches_with_params, has_substruct_match,
};
pub use support::{
    AROMATICITY_FEATURE, BATCH_FEATURE, BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
    BIO_PDB_COORDINATE_SUBSET_READ_FEATURE, BIO_SELECTION_FEATURE, BIO_STRUCTURE_FEATURE,
    CONFORMER_GENERATION_FEATURE, COORDINATE_2D_FEATURE, DG_BOUNDS_FEATURE, DRAWING_FEATURE,
    FINGERPRINT_FEATURE, FeatureCategory, FeatureSpec, HYDROGENS_FEATURE, KEKULIZE_FEATURE,
    MOLBLOCK_IO_FEATURE, PUBLIC_FEATURES, RINGS_FEATURE, SANITIZE_FEATURE, SMILES_PARSE_FEATURE,
    SMILES_WRITE_FEATURE, STEREO_FEATURE, SupportStatus, UnsupportedFeatureError, VALENCE_FEATURE,
};
pub use valence::{
    ValenceAssignment, ValenceError, ValenceModel, assign_radicals, assign_valence,
    assign_valence_with_options, atom_has_valence_violation, rdkit_valence_list,
};

/// Returns the crate version at compile time.
#[must_use]
pub fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}
