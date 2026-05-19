//! COSMolKit core, redesigned around value-style molecule state.
//!
//! This crate intentionally starts small. The previous implementation has been
//! moved to `crates/cosmolkit-core-old` and is reference material only; it must
//! not be used as a dependency for this crate.
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

pub mod bio;
pub mod chemistry;
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
pub use bio::{
    AltLocLabel, AtomName, AtomRow, AtomSourceIds, BioMetadata, BioStructure, ChainKind, ChainRow,
    ChainSourceIds, CoordinateBlock, CrystalCell, CrystalInfo, EntityKind, EntityRow,
    EntitySourceIds, ModelRow, PolymerKind, ResidueKind, ResidueName, ResidueRow, RowSpan,
    classify_residue_name,
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
pub use coordinates::With2DCoordinatesParams;
pub use derived::DerivedState;
pub use distgeom::DgBoundsError;
pub use draw::{PreparedDrawAtom, PreparedDrawBond, PreparedDrawMolecule, SvgDrawError};
pub use error::{InvariantError, MoleculeBuildError};
pub use fingerprint::{
    Fingerprint, FingerprintError, MorganAdditionalOutput, MorganAtomInvariantsGenerator,
    MorganBondInvariantsGenerator, MorganFingerprintOutput, MorganFingerprintParams,
};
pub use hydrogens::{AddHsParams, AddHydrogensError, RemoveHsParams, RemoveHydrogensError};
pub use io::bio::{
    BioPdbReadParams, BioReadError, read_mmcif_atom_site_subset_from_str,
    read_pdb_coordinate_subset_from_str, read_pdb_coordinate_subset_from_str_with_params,
};
pub use io::sdf::{SdfDataset, SdfRecordMetadata};
pub use kekulize::KekulizeError;
pub use mol_pickler::{PickleError, mol_from_binary, mol_to_binary};
pub use molecule::{
    AtomMapping, BondMapping, Conformer2D, Conformer3D, ConformerStore, CoordinateDimension,
    Molecule, MoleculeProperties, PropertyStore, SdfPropertyList, SdfPropertyListTarget,
    SmilesParseError, SmilesWriteError, TopologyMapping,
};
pub use ops::{
    ASSIGNED_AROMATICITY_SPEC, ASSIGNED_RING_FAMILIES_SPEC, ASSIGNED_RINGS_SPEC,
    ASSIGNED_VALENCE_SPEC, BlockAccess, BlockSet, InvariantCheckSet, MOLECULE_OPS,
    MappingRequirement, MoleculeOpKind, MoleculeOpSpec, OPERATION_INVARIANT_MATRIX, OpOutcome,
    OperationDomain, OperationError, OperationInvariantEntry, OperationTrace, PARITY_MATRIX,
    ParityMatrixEntry, ParityPolicy, SANITIZED_SPEC, SUPPORT_MATRIX, SupportMatrixEntry,
    TopologyEditKind, WITH_2D_COORDINATES_SPEC, WITH_HYDROGENS_SPEC, WITH_KEKULIZED_BONDS_SPEC,
    WITHOUT_HYDROGENS_SPEC, WITHOUT_HYDROGENS_WITH_PARAMS_SPEC,
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
    assign_stereochemistry,
};
pub use substruct::{
    SubstructMatchParams, SubstructMatchResult, get_substruct_match, get_substruct_matches,
    get_substruct_matches_with_params, has_substruct_match,
};
pub use support::{
    AROMATICITY_FEATURE, BATCH_FEATURE, BIO_MMCIF_ATOM_SITE_SUBSET_READ_FEATURE,
    BIO_PDB_COORDINATE_SUBSET_READ_FEATURE, BIO_SELECTION_FEATURE, BIO_STRUCTURE_FEATURE,
    COORDINATE_2D_FEATURE, DG_BOUNDS_FEATURE, DRAWING_FEATURE, FINGERPRINT_FEATURE,
    FeatureCategory, FeatureSpec, HYDROGENS_FEATURE, KEKULIZE_FEATURE, MOLBLOCK_IO_FEATURE,
    PUBLIC_FEATURES, RINGS_FEATURE, SANITIZE_FEATURE, SMILES_PARSE_FEATURE, SMILES_WRITE_FEATURE,
    STEREO_FEATURE, SupportStatus, UnsupportedFeatureError, VALENCE_FEATURE,
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
