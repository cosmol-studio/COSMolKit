//! Detached hydrogen transformations.
//!
//! The public functions in this module are the algorithm boundary that the
//! `cosmolkit` runtime will call after extracting authorized blocks. Runtime
//! bookkeeping, cache invalidation, and operation contracts deliberately do
//! not appear here.

use crate::{
    RingInfo, SanitizeError, SanitizeParams, ValenceAssignment, ValenceError, ValenceModel,
    assign_valence_state_for_atom_from_parts, assign_valence_with_options_for_topology, rdkit_rb0,
    rdkit_valence_list, sanitize_topology_with_query_state,
};
use cosmolkit_model::{
    AdjacencyList, AtomId, AtomMapping, AtomPdbResidueInfo, AtomSpec, Bond, BondDirection, BondId,
    BondMapping, BondOrder, BondSpec, BondStereo, ChiralTag, Conformer2D, Conformer3D,
    CoordinateBlock, CoordinateValidationError, Element, Hybridization, MappingValidationError,
    MoleculeProperties, MoleculePropertyError, QueryAtom, QueryBond, QueryStateError,
    QueryStateRef, SGroupBondRole, SdfPropertyListTarget, SubstanceGroup, TopologyBlock,
    TopologyEditError, TopologyMapping, TopologyValidationError, remap_query_rows,
};

/// Parameters corresponding to RDKit's `MolOps::AddHsParameters`.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct AddHsParams {
    pub explicit_only: bool,
    pub add_coords: bool,
    pub add_residue_info: bool,
    pub skip_queries: bool,
    pub only_on_atoms: Option<Vec<AtomId>>,
}

/// Parameters corresponding to RDKit's `MolOps::RemoveHsParameters`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RemoveHsParams {
    pub remove_degree_zero: bool,
    pub remove_higher_degrees: bool,
    pub remove_only_h_neighbors: bool,
    pub remove_isotopes: bool,
    pub remove_and_track_isotopes: bool,
    pub remove_dummy_neighbors: bool,
    pub remove_defining_bond_stereo: bool,
    pub remove_with_wedged_bond: bool,
    pub remove_with_query: bool,
    pub remove_mapped: bool,
    pub remove_in_sgroups: bool,
    pub show_warnings: bool,
    pub remove_nonimplicit: bool,
    pub update_explicit_count: bool,
    pub remove_hydrides: bool,
    pub remove_nontetrahedral_neighbors: bool,
    /// Run the complete default sanitize pipeline after source-eligible
    /// non-implicit removal. RDKit exposes this as the overload's separate
    /// `sanitize` argument; the canonical COSMolKit parameter value keeps the
    /// observable option explicit across languages.
    /// When false, actual pre-removal source scalar values survive;
    /// intermediate calculations needed by removal still run. True alone is
    /// not proof of strict chemical validity when no sanitize branch executes.
    pub sanitize: bool,
}

impl Default for RemoveHsParams {
    fn default() -> Self {
        Self {
            remove_degree_zero: false,
            remove_higher_degrees: false,
            remove_only_h_neighbors: false,
            remove_isotopes: false,
            remove_and_track_isotopes: false,
            remove_dummy_neighbors: false,
            remove_defining_bond_stereo: false,
            remove_with_wedged_bond: true,
            remove_with_query: false,
            remove_mapped: true,
            remove_in_sgroups: true,
            show_warnings: true,
            remove_nonimplicit: true,
            update_explicit_count: false,
            remove_hydrides: false,
            remove_nontetrahedral_neighbors: false,
            sanitize: true,
        }
    }
}

/// Whether an appended hydrogen came from the source atom's explicit or
/// implicit hydrogen count.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AddedHydrogenKind {
    Explicit,
    Implicit,
}

/// Stable metadata for one hydrogen and bond appended by AddHs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AddedHydrogen {
    pub atom: AtomId,
    pub bond: BondId,
    pub parent: AtomId,
    pub kind: AddedHydrogenKind,
}

/// Complete detached output of the topology-only AddHs stage.
#[derive(Debug, Clone, PartialEq)]
pub struct AddHydrogensTopologyResult {
    pub topology: TopologyBlock,
    pub mapping: TopologyMapping,
    pub additions: Vec<AddedHydrogen>,
}

/// Complete detached output of the AddHs coordinate/residue stage.
#[derive(Debug, Clone, PartialEq)]
pub struct AddHydrogensCoordinateResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub mapping: TopologyMapping,
    pub additions: Vec<AddedHydrogen>,
}

/// Complete detached AddHs output for the parent operation runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct AddHydrogensResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub mapping: TopologyMapping,
    pub warnings: Vec<HydrogenWarning>,
    /// Actual selected-parent/new-H source scalars; untouched rows survive.
    pub final_valence: ValenceAssignment,
}

/// Detached topology state prepared for the one downstream RemoveHs compaction.
#[derive(Debug, Clone, PartialEq)]
pub struct PreparedHydrogenRemoval {
    pub topology: TopologyBlock,
    pub atoms_to_remove: Vec<AtomId>,
}

/// Complete detached RemoveHs output for the parent operation runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct RemoveHydrogensResult {
    pub topology: TopologyBlock,
    pub coordinates: CoordinateBlock,
    pub properties: MoleculeProperties,
    pub mapping: TopologyMapping,
    /// Actual source scalars aligned with the final retained atom rows.
    /// Available without sanitization; no unsolicited final-topology refresh.
    pub final_valence: Option<ValenceAssignment>,
    /// Final ring-state transport for the final topology. `None` is the
    /// exact final source-uninitialized state (removeHs ALWAYS
    /// clearComputedProps(true), so every no-sanitize branch reports the
    /// reset state); `Some` is the exact SAN initialized state, INCLUDING
    /// zero rings. Published ONLY from the final pass by MOVE; the
    /// isotope-tracking preliminary pass hardcodes sanitize=false and its
    /// state is discarded. Independent from the retained scalar state:
    /// sanitize=true/remove_nonimplicit=false yields
    /// final_valence Some with final_rings None. No live-cache
    /// installation is implied.
    pub final_rings: Option<RingInfo>,
    pub warnings: Vec<HydrogenWarning>,
}

/// Non-fatal source diagnostics produced while adding hydrogens.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HydrogenWarning {
    ExtraTrackedIsotopes { parent: AtomId, count: usize },
    IsolatedHydrogen { hydrogen: AtomId },
    DummyAtomNeighbor { hydrogen: AtomId, neighbor: AtomId },
    NonTetrahedralStereoNeighbor { hydrogen: AtomId, neighbor: AtomId },
    WedgedBond { hydrogen: AtomId, bond: BondId },
}

/// Errors raised by detached hydrogen algorithms.
#[derive(Debug, Clone, PartialEq)]
pub enum HydrogenError {
    BondProperty(cosmolkit_model::BondValueError),
    AtomProperty(cosmolkit_model::AtomPropertyError),
    InvalidTopology(TopologyValidationError),
    InvalidCoordinates(CoordinateValidationError),
    TopologyEdit(TopologyEditError),
    InvalidMapping(MappingValidationError),
    InvalidQueryState(QueryStateError),
    InvalidProperty(MoleculePropertyError),
    InvalidPropertyList {
        target: SdfPropertyListTarget,
        name: cosmolkit_model::PropertyText,
        expected_rows: usize,
        actual_rows: usize,
    },
    Unsupported {
        operation: &'static str,
        reason: &'static str,
    },
    OnlyOnAtomOutOfRange {
        atom: AtomId,
        atom_count: usize,
    },
    InvalidAdditionPlan {
        addition: Option<usize>,
        reason: &'static str,
    },
    InvalidRemovalCandidate {
        position: usize,
        atom: AtomId,
        reason: &'static str,
    },
    MissingSourceState {
        stage: &'static str,
        state: &'static str,
    },
    ValenceAssignmentLength {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    ExplicitHydrogenOverflow {
        atom: AtomId,
        current: u8,
    },
    InvalidStereoTransition {
        bond: BondId,
        reason: &'static str,
    },
    CoordinatePlacement {
        addition: usize,
        conformer: usize,
        dimension: &'static str,
        reason: &'static str,
    },
    Valence(ValenceError),
    Sanitize(SanitizeError),
}

impl std::fmt::Display for HydrogenError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AtomProperty(error) => {
                write!(formatter, "atom property operation failed: {error}")
            }
            Self::BondProperty(error) => {
                write!(formatter, "bond property operation failed: {error}")
            }
            Self::InvalidTopology(error) => write!(formatter, "invalid detached topology: {error}"),
            Self::InvalidCoordinates(error) => {
                write!(formatter, "invalid detached coordinates: {error}")
            }
            Self::TopologyEdit(error) => {
                write!(formatter, "detached topology edit failed: {error}")
            }
            Self::InvalidMapping(error) => {
                write!(formatter, "invalid detached hydrogen mapping: {error}")
            }
            Self::InvalidQueryState(error) => {
                write!(formatter, "invalid detached hydrogen query state: {error}")
            }
            Self::InvalidProperty(error) => write!(formatter, "invalid hydrogen property: {error}"),
            Self::InvalidPropertyList {
                target,
                name,
                expected_rows,
                actual_rows,
            } => {
                let target = match target {
                    SdfPropertyListTarget::Atom => "atom",
                    SdfPropertyListTarget::Bond => "bond",
                };
                write!(
                    formatter,
                    "SDF {target} property list {name:?} has {actual_rows} rows; expected {expected_rows}"
                )
            }
            Self::Unsupported { operation, reason } => {
                write!(
                    formatter,
                    "{operation} is unsupported for this detached input: {reason}"
                )
            }
            Self::OnlyOnAtomOutOfRange { atom, atom_count } => {
                write!(
                    formatter,
                    "add_hydrogens received only_on_atoms id {atom} outside {atom_count} atoms"
                )
            }
            Self::InvalidAdditionPlan { addition, reason } => match addition {
                Some(addition) => {
                    write!(
                        formatter,
                        "invalid hydrogen addition row {addition}: {reason}"
                    )
                }
                None => write!(formatter, "invalid hydrogen addition plan: {reason}"),
            },
            Self::InvalidRemovalCandidate {
                position,
                atom,
                reason,
            } => write!(
                formatter,
                "invalid hydrogen removal candidate row {position} (atom {atom}): {reason}"
            ),
            Self::MissingSourceState { stage, state } => write!(
                formatter,
                "hydrogen-removal source stage {stage} did not return materialized {state}"
            ),
            Self::ValenceAssignmentLength {
                field,
                expected,
                actual,
            } => write!(
                formatter,
                "hydrogen-removal valence field {field} has {actual} rows; expected {expected}"
            ),
            Self::ExplicitHydrogenOverflow { atom, current } => write!(
                formatter,
                "removing hydrogen would overflow atom {atom} explicit-H count {current}"
            ),
            Self::InvalidStereoTransition { bond, reason } => write!(
                formatter,
                "invalid hydrogen-removal stereo transition for bond {bond}: {reason}"
            ),
            Self::CoordinatePlacement {
                addition,
                conformer,
                dimension,
                reason,
            } => write!(
                formatter,
                "cannot place hydrogen addition {addition} in {dimension} conformer {conformer}: {reason}"
            ),
            Self::Valence(error) => write!(formatter, "valence assignment failed: {error}"),
            Self::Sanitize(error) => write!(formatter, "hydrogen-removal sanitize failed: {error}"),
        }
    }
}

impl std::error::Error for HydrogenError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::AtomProperty(source) => Some(source),
            Self::BondProperty(source) => Some(source),
            Self::InvalidProperty(source) => Some(source),
            _ => None,
        }
    }
}

impl From<cosmolkit_model::AtomPropertyError> for HydrogenError {
    fn from(source: cosmolkit_model::AtomPropertyError) -> Self {
        Self::AtomProperty(source)
    }
}
impl From<cosmolkit_model::BondValueError> for HydrogenError {
    fn from(source: cosmolkit_model::BondValueError) -> Self {
        Self::BondProperty(source)
    }
}

impl From<TopologyValidationError> for HydrogenError {
    fn from(error: TopologyValidationError) -> Self {
        Self::InvalidTopology(error)
    }
}

impl From<CoordinateValidationError> for HydrogenError {
    fn from(error: CoordinateValidationError) -> Self {
        Self::InvalidCoordinates(error)
    }
}

impl From<TopologyEditError> for HydrogenError {
    fn from(error: TopologyEditError) -> Self {
        Self::TopologyEdit(error)
    }
}

impl From<MappingValidationError> for HydrogenError {
    fn from(error: MappingValidationError) -> Self {
        Self::InvalidMapping(error)
    }
}

impl From<QueryStateError> for HydrogenError {
    fn from(error: QueryStateError) -> Self {
        Self::InvalidQueryState(error)
    }
}

impl From<ValenceError> for HydrogenError {
    fn from(error: ValenceError) -> Self {
        Self::Valence(error)
    }
}

impl From<SanitizeError> for HydrogenError {
    fn from(error: SanitizeError) -> Self {
        Self::Sanitize(error)
    }
}

/// Apply the default AddHs operation to detached blocks.
pub fn add_hydrogens_impl(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<AddHydrogensResult, HydrogenError> {
    add_hydrogens_with_params(topology, coordinates, properties, &AddHsParams::default())
}

/// Apply AddHs to detached blocks with explicit source-shaped parameters.
pub fn add_hydrogens_with_params(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &AddHsParams,
) -> Result<AddHydrogensResult, HydrogenError> {
    add_hydrogens_with_query_state(topology, coordinates, properties, params, None)
}

/// Internal typed-query variant used by detached owners that retain the
/// source Atom/QueryAtom and Bond/QueryBond distinction.
#[doc(hidden)]
pub fn add_hydrogens_with_query_state(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &AddHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<AddHydrogensResult, HydrogenError> {
    // Explicit detached preparation; this entry has no stored source scalars.
    // The live runtime calls the separate supplied-scalar seam below instead.
    validate_blocks(&topology, &coordinates)?;
    validate_property_lists(&properties, topology.atoms.len(), topology.bonds.len())?;
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }
    let valence =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)?;
    add_hydrogens_from_source_state(
        topology,
        coordinates,
        properties,
        params,
        query_state,
        Some(valence),
    )
}

fn processed_add_hydrogen_atoms(
    topology: &TopologyBlock,
    params: &AddHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<bool>, HydrogenError> {
    let mut processed = selected_atoms(topology, params.only_on_atoms.as_deref())?;
    if params.skip_queries {
        for atom in 0..topology.atoms.len() {
            if processed[atom] && is_query_atom(topology, AtomId::new(atom), query_state) {
                processed[atom] = false;
            }
        }
    }
    Ok(processed)
}

fn replay_tracked_isotopes(
    topology: &mut TopologyBlock,
    additions: &[AddedHydrogen],
    tracked_isotopes: &[(AtomId, Vec<u16>)],
) -> Result<Vec<HydrogenWarning>, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs tracked-isotope transition
    // RDKit✔️❌:     std::vector<unsigned int> isoHs;
    // RDKit✔️❌:     if (newAt->getPropIfPresent(common_properties::_isotopicHs, isoHs)) {
    // RDKit✔️❌:       newAt->clearProp(common_properties::_isotopicHs);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     std::vector<unsigned int>::const_iterator isoH = isoHs.begin();
    // RDKit✔️❌:       if (isoH != isoHs.end()) {
    // RDKit✔️❌:         hAtom->setIsotope(*isoH);
    // RDKit✔️❌:         ++isoH;
    // RDKit✔️❌:       }
    // RDKit✔️❌:         if (isoH != isoHs.end()) {
    // RDKit✔️❌:           hAtom->setIsotope(*isoH);
    // RDKit✔️❌:           ++isoH;
    // RDKit✔️❌:         }
    // RDKit✔️❌:     if (isoH != isoHs.end()) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog) << "extra H isotope information found on atom "
    // RDKit✔️❌:                               << newAt->getIdx() << std::endl;
    // RDKit✔️❌:     }
    // END RDKIT CPP FUNCTION MolOps::addHs tracked-isotope transition
    // The detached stage necessarily keeps the accepted addition plan and
    // tracked lists until composition. Both traversals are linear, but this
    // retains more temporary storage than RDKit's in-place loop.
    let mut warnings = Vec::new();
    let mut addition_cursor = 0;
    for (parent, isotopes) in tracked_isotopes {
        topology.atoms[parent.index()].set_tracked_isotopic_hydrogens(Vec::new());
        let first_addition = addition_cursor;
        while addition_cursor < additions.len() && additions[addition_cursor].parent == *parent {
            addition_cursor += 1;
        }
        let parent_additions = &additions[first_addition..addition_cursor];
        for (addition, isotope) in parent_additions.iter().zip(isotopes) {
            topology.atoms[addition.atom.index()].set_isotope((*isotope != 0).then_some(*isotope));
        }
        if isotopes.len() > parent_additions.len() {
            warnings.push(HydrogenWarning::ExtraTrackedIsotopes {
                parent: *parent,
                count: isotopes.len() - parent_additions.len(),
            });
        }
    }
    if addition_cursor != additions.len() {
        return Err(HydrogenError::InvalidAdditionPlan {
            addition: Some(addition_cursor),
            reason: "addition parent is absent or out of order in the processed atom sequence",
        });
    }
    Ok(warnings)
}

/// Apply the source AddHs selection and append rules to detached topology.
pub fn add_hydrogens_topology(
    topology: TopologyBlock,
    params: &AddHsParams,
) -> Result<AddHydrogensTopologyResult, HydrogenError> {
    add_hydrogens_topology_with_query_state(topology, params, None)
}

/// Internal typed-query AddHs topology owner.
#[doc(hidden)]
pub fn add_hydrogens_topology_with_query_state(
    topology: TopologyBlock,
    params: &AddHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<AddHydrogensTopologyResult, HydrogenError> {
    topology.validate()?;
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }
    let valence =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)?;
    add_hydrogens_topology_from_source_state(topology, params, query_state, Some(valence))
        .map(|(result, _)| result)
}

/// Grow conformers for an accepted topology addition plan and apply the
/// source-controlled coordinate/residue stages.
pub fn add_hydrogen_coordinates(
    mut result: AddHydrogensTopologyResult,
    coordinates: CoordinateBlock,
    add_coords: bool,
    add_residue_info: bool,
) -> Result<AddHydrogensCoordinateResult, HydrogenError> {
    let old_atom_count = validate_addition_plan(&result, &coordinates)?;
    let mut coordinates = grow_coordinate_rows(coordinates, result.additions.len());

    if add_coords {
        place_added_hydrogens(&result.topology, &result.additions, &mut coordinates)?;
    }
    if add_residue_info {
        assign_hydrogen_residue_info(&mut result.topology)?;
    }

    debug_assert_eq!(
        old_atom_count + result.additions.len(),
        result.topology.atoms.len()
    );
    coordinates.validate_for_atom_count(result.topology.atoms.len())?;
    Ok(AddHydrogensCoordinateResult {
        topology: result.topology,
        coordinates,
        mapping: result.mapping,
        additions: result.additions,
    })
}

fn assign_hydrogen_residue_info(topology: &mut TopologyBlock) -> Result<(), HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION AssignHsResidueInfo
    // RDKit❗❌: void AssignHsResidueInfo(RWMol &mol) {
    // RDKit❗❌:   int max_serial = 0;
    // RDKit❗❌:   unsigned int stopIdx = mol.getNumAtoms();
    // RDKit❗❌:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit❗❌:     auto *info =
    // RDKit❗❌:         (AtomPDBResidueInfo *)(mol.getAtomWithIdx(aidx)->getMonomerInfo());
    // RDKit❗❌:     if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
    // RDKit❗❌:         info->getSerialNumber() > max_serial) {
    // RDKit❗❌:       max_serial = info->getSerialNumber();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   AtomPDBResidueInfo *current_info = nullptr;
    // RDKit❗❌:   int current_h_id = 0;
    // RDKit❗❌:   for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit❗❌:     Atom *newAt = mol.getAtomWithIdx(aidx);
    // RDKit❗❌:     auto *info = (AtomPDBResidueInfo *)(newAt->getMonomerInfo());
    // RDKit❗❌:     if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
    // RDKit❗❌:       ROMol::ADJ_ITER begin, end;
    // RDKit❗❌:       boost::tie(begin, end) = mol.getAtomNeighbors(newAt);
    // RDKit❗❌:       while (begin != end) {
    // RDKit❗❌:         if (mol.getAtomWithIdx(*begin)->getAtomicNum() == 1) {
    // RDKit❗❌:           // Make all Hs unique - increment id even for existing
    // RDKit❗❌:           ++current_h_id;
    // RDKit❗❌:           // skip if hydrogen already has PDB info
    // RDKit❗❌:           auto *h_info = (AtomPDBResidueInfo *)mol.getAtomWithIdx(*begin)
    // RDKit❗❌:                              ->getMonomerInfo();
    // RDKit❗❌:           if (h_info &&
    // RDKit❗❌:               h_info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
    // RDKit❗❌:             continue;
    // RDKit❗❌:           }
    // RDKit❗❌:           // the hydrogens have unique names on residue basis (H1, H2, ...)
    // RDKit❗❌:           if (!current_info ||
    // RDKit❗❌:               current_info->getResidueNumber() != info->getResidueNumber() ||
    // RDKit❗❌:               current_info->getChainId() != info->getChainId()) {
    // RDKit❗❌:             current_h_id = 1;
    // RDKit❗❌:             current_info = info;
    // RDKit❗❌:           }
    // RDKit❗❌:           std::string h_label = std::to_string(current_h_id);
    // RDKit❗❌:           if (h_label.length() > 3) {
    // RDKit❗❌:             h_label = h_label.substr(h_label.length() - 3, 3);
    // RDKit❗❌:           }
    // RDKit❗❌:           while (h_label.length() < 3) {
    // RDKit❗❌:             h_label = h_label + " ";
    // RDKit❗❌:           }
    // RDKit❗❌:           h_label = "H" + h_label;
    // RDKit❗❌:           // wrap around id to '3H12'
    // RDKit❗❌:           h_label = h_label.substr(3, 1) + h_label.substr(0, 3);
    // RDKit❗❌:           AtomPDBResidueInfo *newInfo = new AtomPDBResidueInfo(
    // RDKit❗❌:               h_label, max_serial, "", info->getResidueName(),
    // RDKit❗❌:               info->getResidueNumber(), info->getChainId(), "", 1.0, 0.0,
    // RDKit❗❌:               info->getIsHeteroAtom());
    // RDKit❗❌:           mol.getAtomWithIdx(*begin)->setMonomerInfo(newInfo);
    // RDKit❗❌:
    // RDKit❗❌:           ++max_serial;
    // RDKit❗❌:         }
    // RDKit❗❌:         ++begin;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION AssignHsResidueInfo
    let stop_index = topology.atoms.len();
    let mut max_serial = topology
        .atoms
        .iter()
        .filter_map(|atom| atom.pdb_residue_info())
        .map(AtomPdbResidueInfo::serial_number)
        .filter(|serial| *serial > 0)
        .max()
        .unwrap_or(0);
    let mut current_info: Option<(i32, String)> = None;
    let mut current_h_id = 0_i32;

    for atom_index in 0..stop_index {
        let Some(info) = topology.atoms[atom_index].pdb_residue_info().cloned() else {
            continue;
        };
        let neighbors = topology.adjacency.neighbors_of(atom_index).to_vec();
        for neighbor in neighbors {
            let hydrogen_index = neighbor.atom_index;
            if topology.atoms[hydrogen_index].atomic_number() != 1 {
                continue;
            }
            current_h_id =
                current_h_id
                    .checked_add(1)
                    .ok_or(HydrogenError::InvalidAdditionPlan {
                        addition: None,
                        reason: "PDB hydrogen residue id exceeds the typed integer range",
                    })?;
            if topology.atoms[hydrogen_index].pdb_residue_info().is_some() {
                continue;
            }
            let residue_key = (info.residue_number(), info.chain_id().to_owned());
            if current_info.as_ref() != Some(&residue_key) {
                current_h_id = 1;
                current_info = Some(residue_key);
            }
            let mut id = current_h_id.to_string();
            if id.len() > 3 {
                id = id[id.len() - 3..].to_owned();
            }
            while id.len() < 3 {
                id.push(' ');
            }
            let prefixed = format!("H{id}");
            let atom_name = format!("{}{}", &prefixed[3..4], &prefixed[..3]);
            let new_info = AtomPdbResidueInfo::new(
                atom_name,
                max_serial,
                info.residue_name(),
                info.residue_number(),
                info.chain_id(),
                info.is_hetero_atom(),
            );
            topology.atoms[hydrogen_index].set_pdb_residue_info(Some(new_info));
            max_serial = max_serial
                .checked_add(1)
                .ok_or(HydrogenError::InvalidAdditionPlan {
                    addition: None,
                    reason: "PDB hydrogen serial exceeds the typed integer range",
                })?;
        }
    }
    Ok(())
}

fn place_added_hydrogens(
    topology: &TopologyBlock,
    additions: &[AddedHydrogen],
    coordinates: &mut CoordinateBlock,
) -> Result<(), HydrogenError> {
    if coordinates.conformers_2d.is_empty() && coordinates.conformers_3d.is_empty() {
        return Ok(());
    }
    for (addition_index, addition) in additions.iter().enumerate() {
        for conformer in &mut coordinates.conformers_2d {
            let mut values = conformer
                .coordinates()
                .iter()
                .copied()
                .map(AddHsPoint3D::from_2d)
                .collect::<Vec<_>>();
            let position = terminal_position(
                topology,
                addition,
                addition_index,
                conformer.id(),
                "2D",
                false,
                &values,
            )?;
            values[addition.atom.index()] = position;
            conformer.coordinates_mut()[addition.atom.index()] = position.to_2d();
        }
        for conformer in &mut coordinates.conformers_3d {
            let mut values = conformer
                .coordinates()
                .iter()
                .copied()
                .map(AddHsPoint3D::from_3d)
                .collect::<Vec<_>>();
            let position = terminal_position(
                topology,
                addition,
                addition_index,
                conformer.id(),
                "3D",
                conformer.is_3d(),
                &values,
            )?;
            values[addition.atom.index()] = position;
            conformer.coordinates_mut()[addition.atom.index()] = position.to_3d();
        }
    }
    Ok(())
}

/// Place a newly appended attachment dummy using the same terminal-atom
/// geometry owner as AddHs. The input already has a zero-initialized row for
/// `atom`, as RWMol::addAtom does for each conformer. Taking the detached block
/// by value keeps a failed placement from exposing partially written rows.
pub fn place_terminal_attachment_coordinates(
    topology: &TopologyBlock,
    mut coordinates: CoordinateBlock,
    atom: AtomId,
    parent: AtomId,
    bond: BondId,
) -> Result<CoordinateBlock, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION details::addExplicitAttachmentPoint coordinate call
    // RDKit❗❌:   if (addCoords) {
    // RDKit❗❌:     setTerminalAtomCoords(mol, idx, atomIdx);
    // RDKit❗❌:   }
    // END RDKIT CPP FUNCTION details::addExplicitAttachmentPoint coordinate call
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    let attachment_bond =
        topology
            .bonds
            .get(bond.index())
            .ok_or(HydrogenError::InvalidAdditionPlan {
                addition: None,
                reason: "attachment bond is missing",
            })?;
    if atom == parent
        || atom.index() >= topology.atoms.len()
        || attachment_bond.begin() != parent
        || attachment_bond.end() != atom
        || active_neighbors(topology, atom, bond) != [parent]
    {
        return Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "terminal attachment preconditions are not satisfied",
        });
    }
    let addition = AddedHydrogen {
        atom,
        bond,
        parent,
        kind: AddedHydrogenKind::Explicit,
    };
    for conformer in &mut coordinates.conformers_2d {
        let values = conformer
            .coordinates()
            .iter()
            .copied()
            .map(AddHsPoint3D::from_2d)
            .collect::<Vec<_>>();
        let position =
            terminal_position(topology, &addition, 0, conformer.id(), "2D", false, &values)?;
        conformer.coordinates_mut()[atom.index()] = position.to_2d();
    }
    for conformer in &mut coordinates.conformers_3d {
        let values = conformer
            .coordinates()
            .iter()
            .copied()
            .map(AddHsPoint3D::from_3d)
            .collect::<Vec<_>>();
        let position = terminal_position(
            topology,
            &addition,
            0,
            conformer.id(),
            "3D",
            conformer.is_3d(),
            &values,
        )?;
        conformer.coordinates_mut()[atom.index()] = position.to_3d();
    }
    Ok(coordinates)
}

fn active_neighbors(topology: &TopologyBlock, atom: AtomId, through_bond: BondId) -> Vec<AtomId> {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| neighbor.bond.index() <= through_bond.index())
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .collect()
}

fn active_neighbor_not(
    topology: &TopologyBlock,
    atom: AtomId,
    other: AtomId,
    through_bond: BondId,
    addition: usize,
) -> Result<AtomId, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION getAtomNeighborNot
    // RDKit❗❌: Atom *getAtomNeighborNot(ROMol *mol, const Atom *atom, const Atom *other) {
    // RDKit❗❌:   PRECONDITION(mol, "bad molecule");
    // RDKit❗❌:   PRECONDITION(atom, "bad atom");
    // RDKit❗❌:   PRECONDITION(atom->getDegree() > 1, "bad degree");
    // RDKit❗❌:   PRECONDITION(other, "bad atom");
    // RDKit❗❌:   Atom *res = nullptr;
    // RDKit❗❌:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗❌:   boost::tie(nbrIdx, endNbrs) = mol->getAtomNeighbors(atom);
    // RDKit❗❌:   while (nbrIdx != endNbrs) {
    // RDKit❗❌:     if (*nbrIdx != other->getIdx()) {
    // RDKit❗❌:       res = mol->getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:     ++nbrIdx;
    // RDKit❗❌:   }
    // RDKit❗❌:   POSTCONDITION(res, "no neighbor found");
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION getAtomNeighborNot
    active_neighbors(topology, atom, through_bond)
        .into_iter()
        .find(|neighbor| *neighbor != other)
        .ok_or(HydrogenError::InvalidAdditionPlan {
            addition: Some(addition),
            reason: "source neighbor-not precondition has no matching active neighbor",
        })
}

#[allow(clippy::too_many_arguments)]
fn terminal_position(
    topology: &TopologyBlock,
    addition: &AddedHydrogen,
    addition_index: usize,
    conformer_id: usize,
    dimension: &'static str,
    is_3d: bool,
    coordinates: &[AddHsPoint3D],
) -> Result<AddHsPoint3D, HydrogenError> {
    match active_neighbors(topology, addition.parent, addition.bond).len() {
        1 | 2 => terminal_position_degree_one_two(
            topology,
            addition,
            addition_index,
            conformer_id,
            dimension,
            is_3d,
            coordinates,
        ),
        _ => terminal_position_degree_three_four_default(
            topology,
            addition,
            addition_index,
            conformer_id,
            dimension,
            is_3d,
            coordinates,
        ),
    }
}

#[allow(clippy::too_many_arguments)]
fn terminal_position_degree_one_two(
    topology: &TopologyBlock,
    addition: &AddedHydrogen,
    addition_index: usize,
    conformer_id: usize,
    dimension: &'static str,
    is_3d: bool,
    coordinates: &[AddHsPoint3D],
) -> Result<AddHsPoint3D, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setTerminalAtomCoords degree one/two
    // RDKit❗❌: void setTerminalAtomCoords(ROMol &mol, unsigned int idx,
    // RDKit❗❌:                            unsigned int otherIdx) {
    // RDKit❗❌:   // we will loop over all the coordinates
    // RDKit❗❌:   PRECONDITION(otherIdx != idx, "degenerate atoms");
    // RDKit❗❌:   Atom *atom = mol.getAtomWithIdx(idx);
    // RDKit❗❌:   PRECONDITION(mol.getAtomDegree(atom) == 1, "bad atom degree");
    // RDKit❗❌:   const Bond *bond = mol.getBondBetweenAtoms(otherIdx, idx);
    // RDKit❗❌:   PRECONDITION(bond, "no bond between atoms");
    // RDKit❗❌:
    // RDKit❗❌:   const Atom *otherAtom = mol.getAtomWithIdx(otherIdx);
    // RDKit❗❌:   double bondLength =
    // RDKit❗❌:       PeriodicTable::getTable()->getRb0(1) +
    // RDKit❗❌:       PeriodicTable::getTable()->getRb0(otherAtom->getAtomicNum());
    // RDKit❗❌:
    // RDKit❗❌:   RDGeom::Point3D dirVect(0, 0, 0);
    // RDKit❗❌:
    // RDKit❗❌:   RDGeom::Point3D perpVect, rotnAxis, nbrPerp;
    // RDKit❗❌:   RDGeom::Point3D nbr1Vect, nbr2Vect, nbr3Vect;
    // RDKit❗❌:   RDGeom::Transform3D tform;
    // RDKit❗❌:   RDGeom::Point3D otherPos, atomPos;
    // RDKit❗❌:
    // RDKit❗❌:   const Atom *nbr1 = nullptr, *nbr2 = nullptr, *nbr3 = nullptr;
    // RDKit❗❌:   const Bond *nbrBond;
    // RDKit❗❌:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit❗❌:
    // RDKit❗❌:   switch (otherAtom->getDegree()) {
    // RDKit❗❌:     case 1:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       //   No other atoms present:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       // loop over the conformations and set the coordinates
    // RDKit❗❌:       for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
    // RDKit❗❌:            cfi++) {
    // RDKit❗❌:         if ((*cfi)->is3D()) {
    // RDKit❗❌:           dirVect.z = 1;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           dirVect.x = 1;
    // RDKit❗❌:         }
    // RDKit❗❌:         otherPos = (*cfi)->getAtomPos(otherIdx);
    // RDKit❗❌:         atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:         (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:
    // RDKit❗❌:     case 2:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       //  One other neighbor:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       nbr1 = getAtomNeighborNot(&mol, otherAtom, atom);
    // RDKit❗❌:       for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
    // RDKit❗❌:            ++cfi) {
    // RDKit❗❌:         otherPos = (*cfi)->getAtomPos(otherIdx);
    // RDKit❗❌:         RDGeom::Point3D nbr1Pos = (*cfi)->getAtomPos(nbr1->getIdx());
    // RDKit❗❌:         // get a normalized vector pointing away from the neighbor:
    // RDKit❗❌:         nbr1Vect = nbr1Pos - otherPos;
    // RDKit❗❌:         if (nbr1Vect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:           // no difference, which likely indicates that we have redundant atoms.
    // RDKit❗❌:           // just put it on top of the heavy atom. This was #678
    // RDKit❗❌:           (*cfi)->setAtomPos(idx, otherPos);
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         nbr1Vect.normalize();
    // RDKit❗❌:         nbr1Vect *= -1;
    // RDKit❗❌:
    // RDKit❗❌:         // ok, nbr1Vect points away from the other atom, figure out where
    // RDKit❗❌:         // this H goes:
    // RDKit❗❌:         switch (otherAtom->getHybridization()) {
    // RDKit❗❌:           case Atom::SP3:
    // RDKit❗❌:             // get a perpendicular to nbr1Vect:
    // RDKit❗❌:             if ((*cfi)->is3D()) {
    // RDKit❗❌:               perpVect = nbr1Vect.getPerpendicular();
    // RDKit❗❌:             } else {
    // RDKit❗❌:               perpVect.z = 1.0;
    // RDKit❗❌:             }
    // RDKit❗❌:             // and move off it:
    // RDKit❗❌:             tform.SetRotation((180 - 109.471) * M_PI / 180., perpVect);
    // RDKit❗❌:             dirVect = tform * nbr1Vect;
    // RDKit❗❌:             atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:             (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:             break;
    // RDKit❗❌:           case Atom::SP2:
    // RDKit❗❌:             // default 3D position is to just take an arbitrary perpendicular
    // RDKit❗❌:             // for 2D we take the normal to the xy plane
    // RDKit❗❌:             if ((*cfi)->is3D()) {
    // RDKit❗❌:               perpVect = nbr1Vect.getPerpendicular();
    // RDKit❗❌:             } else {
    // RDKit❗❌:               perpVect.z = 1.0;
    // RDKit❗❌:             }
    // RDKit❗❌:             if (nbr1->getDegree() > 1) {
    // RDKit❗❌:               // can we use the neighboring atom to establish a perpendicular?
    // RDKit❗❌:               nbrBond = mol.getBondBetweenAtoms(otherIdx, nbr1->getIdx());
    // RDKit❗❌:               if (nbrBond->getIsAromatic() ||
    // RDKit❗❌:                   nbrBond->getBondType() == Bond::DOUBLE ||
    // RDKit❗❌:                   nbrBond->getIsConjugated()) {
    // RDKit❗❌:                 nbr2 = getAtomNeighborNot(&mol, nbr1, otherAtom);
    // RDKit❗❌:                 nbr2Vect =
    // RDKit❗❌:                     nbr1Pos.directionVector((*cfi)->getAtomPos(nbr2->getIdx()));
    // RDKit❗❌:                 auto crossProd = nbr2Vect.crossProduct(nbr1Vect);
    // RDKit❗❌:
    // RDKit❗❌:                 // if nbr1 and nbr2 are aligned, the perpendicular will be null,
    // RDKit❗❌:                 // and we'll just keep the default calculated above. Otherwise
    // RDKit❗❌:                 // we use the cross product
    // RDKit❗❌:                 if (crossProd.lengthSq() >= sq_dist_zero_tol) {
    // RDKit❗❌:                   perpVect = crossProd;
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             perpVect.normalize();
    // RDKit❗❌:             // rotate the nbr1Vect 60 degrees about perpVect and we're done:
    // RDKit❗❌:             tform.SetRotation(60. * M_PI / 180., perpVect);
    // RDKit❗❌:             dirVect = tform * nbr1Vect;
    // RDKit❗❌:             atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:             (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:             break;
    // RDKit❗❌:           case Atom::SP:
    // RDKit❗❌:             // just lay the H along the vector:
    // RDKit❗❌:             dirVect = nbr1Vect;
    // RDKit❗❌:             atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:             (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:             break;
    // RDKit❗❌:           default:
    // RDKit❗❌:             // FIX: handle other hybridizations
    // RDKit❗❌:             // for now, just lay the H along the vector:
    // RDKit❗❌:             dirVect = nbr1Vect;
    // RDKit❗❌:             atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:             (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // END RDKIT CPP FUNCTION MolOps::setTerminalAtomCoords degree one/two
    let active = active_neighbors(topology, addition.parent, addition.bond);
    let parent_position = coordinates[addition.parent.index()];
    let distance = if is_3d {
        rdkit_rb0(1) + rdkit_rb0(topology.atoms[addition.parent.index()].atomic_number())
    } else {
        1.0
    };
    match active.len() {
        1 => {
            // INTENTIONAL RDKIT DESIGN DIVERGENCE CK-COORD-001 (approved
            // 2026-09-21): initialize the direction independently per conformer.
            // RDKit 2026.03.1 AddHs.cpp::setTerminalAtomCoords case 1 keeps
            // dirVect outside the conformer loop and only assigns x OR z.
            // Mixed flags therefore leak the preceding direction into later
            // conformers, changing both the axis and the displacement length.
            // This is a deliberate correctness/design optimization, NOT exact
            // RDKit parity or a heuristic fallback: false -> +X at unit length;
            // true -> +Z at the source rb0 distance. Preserve the parent's XYZ,
            // including existing nonzero Z on false-flag XYZ input. No legacy
            // contamination mode is offered. Scope: this degree-one branch
            // only; do not generalize the exception to other source branches.
            // Complexity: constant scratch space/work, no extra allocation.
            // Keep the source behavior marker non-exact and the counterexample
            // in IO-mol_post.md; attachment reuse must retain this contract.
            let direction = if is_3d {
                AddHsPoint3D::new(0.0, 0.0, 1.0)
            } else {
                AddHsPoint3D::new(1.0, 0.0, 0.0)
            };
            Ok(parent_position.plus(direction.scaled(distance)))
        }
        2 => {
            let neighbor = active_neighbor_not(
                topology,
                addition.parent,
                addition.atom,
                addition.bond,
                addition_index,
            )?;
            let neighbor_position = coordinates[neighbor.index()];
            let toward_neighbor = neighbor_position.minus(parent_position);
            if toward_neighbor.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                return Ok(parent_position);
            }
            let away = toward_neighbor
                .normalized()
                .expect("coordinate zero tolerance proves a nonzero vector")
                .scaled(-1.0);
            let hybridization = topology.atoms[addition.parent.index()].hybridization();
            let direction = match hybridization {
                Hybridization::Sp3 => {
                    let axis = if is_3d {
                        away.perpendicular()
                    } else {
                        Some(AddHsPoint3D::new(0.0, 0.0, 1.0))
                    };
                    away.rotated(
                        (180.0_f64 - 109.471).to_radians(),
                        axis.ok_or(HydrogenError::CoordinatePlacement {
                            addition: addition_index,
                            conformer: conformer_id,
                            dimension,
                            reason: "SP3 perpendicular axis is zero",
                        })?,
                    )
                    .ok_or(HydrogenError::CoordinatePlacement {
                        addition: addition_index,
                        conformer: conformer_id,
                        dimension,
                        reason: "SP3 rotation axis is zero",
                    })?
                }
                Hybridization::Sp2 => {
                    let mut axis = if is_3d {
                        away.perpendicular()
                    } else {
                        Some(AddHsPoint3D::new(0.0, 0.0, 1.0))
                    }
                    .ok_or(HydrogenError::CoordinatePlacement {
                        addition: addition_index,
                        conformer: conformer_id,
                        dimension,
                        reason: "SP2 default perpendicular axis is zero",
                    })?;
                    let neighbor_active = active_neighbors(topology, neighbor, addition.bond);
                    if neighbor_active.len() > 1 {
                        let connecting_bond = topology
                            .adjacency
                            .neighbors_of(addition.parent.index())
                            .iter()
                            .find(|entry| {
                                entry.atom_index == neighbor.index()
                                    && entry.bond.index() <= addition.bond.index()
                            })
                            .map(|entry| &topology.bonds[entry.bond.index()])
                            .ok_or(HydrogenError::InvalidAdditionPlan {
                                addition: Some(addition_index),
                                reason: "active parent-neighbor bond is missing",
                            })?;
                        if connecting_bond.is_aromatic()
                            || connecting_bond.order() == BondOrder::Double
                            || connecting_bond.is_conjugated()
                        {
                            let neighbor_two = active_neighbor_not(
                                topology,
                                neighbor,
                                addition.parent,
                                addition.bond,
                                addition_index,
                            )?;
                            let neighbor_direction = coordinates[neighbor_two.index()]
                                .minus(neighbor_position)
                                .normalized()
                                .ok_or(HydrogenError::CoordinatePlacement {
                                    addition: addition_index,
                                    conformer: conformer_id,
                                    dimension,
                                    reason: "SP2 neighboring direction is zero",
                                })?;
                            let cross = neighbor_direction.cross(away);
                            if cross.length_sq() >= ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                                axis = cross;
                            }
                        }
                    }
                    away.rotated(60.0_f64.to_radians(), axis).ok_or(
                        HydrogenError::CoordinatePlacement {
                            addition: addition_index,
                            conformer: conformer_id,
                            dimension,
                            reason: "SP2 rotation axis is zero",
                        },
                    )?
                }
                _ => away,
            };
            Ok(parent_position.plus(direction.scaled(distance)))
        }
        _ => Err(HydrogenError::InvalidAdditionPlan {
            addition: Some(addition_index),
            reason: "degree-one/two coordinate helper received another active degree",
        }),
    }
}

fn pick_bisector(
    neighbor_one: AddHsPoint3D,
    neighbor_two: AddHsPoint3D,
    neighbor_three: AddHsPoint3D,
) -> AddHsPoint3D {
    // BEGIN RDKIT CPP FUNCTION pickBisector
    // RDKit❗❌: RDGeom::Point3D pickBisector(const RDGeom::Point3D &nbr1Vect,
    // RDKit❗❌:                              const RDGeom::Point3D &nbr2Vect,
    // RDKit❗❌:                              const RDGeom::Point3D &nbr3Vect) {
    // RDKit❗❌:   auto dirVect = nbr2Vect + nbr3Vect;
    // RDKit❗❌:   if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:     // nbr2Vect and nbr3Vect are anti-parallel (was #3854)
    // RDKit❗❌:     dirVect = nbr2Vect;
    // RDKit❗❌:     std::swap(dirVect.x, dirVect.y);
    // RDKit❗❌:     dirVect.x *= -1;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (dirVect.dotProduct(nbr1Vect) < 0) {
    // RDKit❗❌:     dirVect *= -1;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return dirVect;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION pickBisector
    let mut direction = neighbor_two.plus(neighbor_three);
    if direction.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
        direction = AddHsPoint3D::new(-neighbor_two.y, neighbor_two.x, neighbor_two.z);
    }
    if direction.dot(neighbor_one) < 0.0 {
        direction = direction.scaled(-1.0);
    }
    direction
}

#[allow(clippy::too_many_arguments)]
fn terminal_position_degree_three_four_default(
    topology: &TopologyBlock,
    addition: &AddedHydrogen,
    addition_index: usize,
    conformer_id: usize,
    dimension: &'static str,
    is_3d: bool,
    coordinates: &[AddHsPoint3D],
) -> Result<AddHsPoint3D, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::setTerminalAtomCoords degree three/four/default
    // RDKit❗❌:     case 3:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       // Two other neighbors:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);
    // RDKit❗❌:       while (nbrIdx != endNbrs) {
    // RDKit❗❌:         if (*nbrIdx != idx) {
    // RDKit❗❌:           if (!nbr1) {
    // RDKit❗❌:             nbr1 = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             nbr2 = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         ++nbrIdx;
    // RDKit❗❌:       }
    // RDKit❗❌:       TEST_ASSERT(nbr1);
    // RDKit❗❌:       TEST_ASSERT(nbr2);
    // RDKit❗❌:       for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
    // RDKit❗❌:            ++cfi) {
    // RDKit❗❌:         // start along the average of the two vectors:
    // RDKit❗❌:         otherPos = (*cfi)->getAtomPos(otherIdx);
    // RDKit❗❌:         nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
    // RDKit❗❌:         nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
    // RDKit❗❌:         if (nbr1Vect.lengthSq() < sq_dist_zero_tol ||
    // RDKit❗❌:             nbr2Vect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:           // no difference, which likely indicates that we have redundant atoms.
    // RDKit❗❌:           // just put it on top of the heavy atom. This was #678
    // RDKit❗❌:           (*cfi)->setAtomPos(idx, otherPos);
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         nbr1Vect.normalize();
    // RDKit❗❌:         nbr2Vect.normalize();
    // RDKit❗❌:         dirVect = nbr1Vect + nbr2Vect;
    // RDKit❗❌:
    // RDKit❗❌:         if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:           // nbr1Vect and nbr2Vect are non-null, but they may
    // RDKit❗❌:           // still cancel each other out
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         dirVect.normalize();
    // RDKit❗❌:         if ((*cfi)->is3D()) {
    // RDKit❗❌:           switch (otherAtom->getHybridization()) {
    // RDKit❗❌:             case Atom::SP3:
    // RDKit❗❌:               // get the perpendicular to the neighbors:
    // RDKit❗❌:               nbrPerp = nbr1Vect.crossProduct(nbr2Vect);
    // RDKit❗❌:               // and the perpendicular to that:
    // RDKit❗❌:               rotnAxis = nbrPerp.crossProduct(dirVect);
    // RDKit❗❌:               // and then rotate about that:
    // RDKit❗❌:               rotnAxis.normalize();
    // RDKit❗❌:               tform.SetRotation((109.471 / 2) * M_PI / 180., rotnAxis);
    // RDKit❗❌:               dirVect = tform * dirVect;
    // RDKit❗❌:               atomPos =
    // RDKit❗❌:                   otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:               (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:               break;
    // RDKit❗❌:             case Atom::SP2:
    // RDKit❗❌:               // don't need to do anything here, the H atom goes right on the
    // RDKit❗❌:               // direction vector
    // RDKit❗❌:               atomPos =
    // RDKit❗❌:                   otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:               (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:               break;
    // RDKit❗❌:             default:
    // RDKit❗❌:               // FIX: handle other hybridizations
    // RDKit❗❌:               // for now, just lay the H along the neighbor vector;
    // RDKit❗❌:               atomPos =
    // RDKit❗❌:                   otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:               (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:               break;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // don't need to do anything here, the H atom goes right on the
    // RDKit❗❌:           // direction vector
    // RDKit❗❌:           atomPos = otherPos + dirVect;
    // RDKit❗❌:           (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:     case 4:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       // Three other neighbors:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       boost::tie(nbrIdx, endNbrs) = mol.getAtomNeighbors(otherAtom);
    // RDKit❗❌:
    // RDKit❗❌:       // We're using chiral tag for checking chirality, so we just take the
    // RDKit❗❌:       // initial order
    // RDKit❗❌:       while (nbrIdx != endNbrs) {
    // RDKit❗❌:         if (*nbrIdx != idx) {
    // RDKit❗❌:           if (!nbr1) {
    // RDKit❗❌:             nbr1 = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           } else if (!nbr2) {
    // RDKit❗❌:             nbr2 = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             nbr3 = mol.getAtomWithIdx(*nbrIdx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         ++nbrIdx;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       TEST_ASSERT(nbr1);
    // RDKit❗❌:       TEST_ASSERT(nbr2);
    // RDKit❗❌:       TEST_ASSERT(nbr3);
    // RDKit❗❌:
    // RDKit❗❌:       for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
    // RDKit❗❌:            ++cfi) {
    // RDKit❗❌:         otherPos = (*cfi)->getAtomPos(otherIdx);
    // RDKit❗❌:         nbr1Vect = otherPos - (*cfi)->getAtomPos(nbr1->getIdx());
    // RDKit❗❌:         nbr2Vect = otherPos - (*cfi)->getAtomPos(nbr2->getIdx());
    // RDKit❗❌:         nbr3Vect = otherPos - (*cfi)->getAtomPos(nbr3->getIdx());
    // RDKit❗❌:         if (nbr1Vect.lengthSq() < sq_dist_zero_tol ||
    // RDKit❗❌:             nbr2Vect.lengthSq() < sq_dist_zero_tol ||
    // RDKit❗❌:             nbr3Vect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:           // no difference, which likely indicates that we have redundant atoms.
    // RDKit❗❌:           // just put it on top of the heavy atom. This was #678
    // RDKit❗❌:           (*cfi)->setAtomPos(idx, otherPos);
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         nbr1Vect.normalize();
    // RDKit❗❌:         nbr2Vect.normalize();
    // RDKit❗❌:         nbr3Vect.normalize();
    // RDKit❗❌:
    // RDKit❗❌:         // if three neighboring atoms are more or less planar, this
    // RDKit❗❌:         // is going to be in a quasi-random (but almost definitely bad)
    // RDKit❗❌:         // direction...
    // RDKit❗❌:         // correct for this (issue 2951221):
    // RDKit❗❌:         if ((*cfi)->is3D()) {
    // RDKit❗❌:           if (fabs(nbr3Vect.dotProduct(nbr1Vect.crossProduct(nbr2Vect))) <
    // RDKit❗❌:               0.1) {
    // RDKit❗❌:             // compute the normal:
    // RDKit❗❌:             dirVect = nbr1Vect.crossProduct(nbr2Vect);
    // RDKit❗❌:
    // RDKit❗❌:             // Each of the nbr vectors is non-null, but there might be pairs
    // RDKit❗❌:             // that cancel each other out. Try to find a direction from atoms
    // RDKit❗❌:             // that do not overlap.
    // RDKit❗❌:             if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:               // This definition of dirVect reverses the parity around otherIdx
    // RDKit❗❌:               // the change of sign restores it
    // RDKit❗❌:               dirVect = nbr1Vect.crossProduct(nbr3Vect) * -1;
    // RDKit❗❌:             }
    // RDKit❗❌:             if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:               dirVect = nbr2Vect.crossProduct(nbr3Vect);
    // RDKit❗❌:             }
    // RDKit❗❌:             // We couldn't find a good direction
    // RDKit❗❌:             if (dirVect.lengthSq() < sq_dist_zero_tol) {
    // RDKit❗❌:               continue;
    // RDKit❗❌:             }
    // RDKit❗❌:
    // RDKit❗❌:             std::string cipCode;
    // RDKit❗❌:             if (otherAtom->getPropIfPresent(common_properties::_CIPCode,
    // RDKit❗❌:                                             cipCode)) {
    // RDKit❗❌:               // the heavy atom is a chiral center, make sure
    // RDKit❗❌:               // that we went go the right direction to preserve
    // RDKit❗❌:               // its chirality. We use the chiral volume for this:
    // RDKit❗❌:               RDGeom::Point3D v1 = dirVect - nbr3Vect;
    // RDKit❗❌:               RDGeom::Point3D v2 = nbr1Vect - nbr3Vect;
    // RDKit❗❌:               RDGeom::Point3D v3 = nbr2Vect - nbr3Vect;
    // RDKit❗❌:               double vol = v1.dotProduct(v2.crossProduct(v3));
    // RDKit❗❌:
    // RDKit❗❌:               if ((otherAtom->getChiralTag() ==
    // RDKit❗❌:                        Atom::ChiralType::CHI_TETRAHEDRAL_CCW &&
    // RDKit❗❌:                    vol < 0) ||
    // RDKit❗❌:                   (otherAtom->getChiralTag() ==
    // RDKit❗❌:                        Atom::ChiralType::CHI_TETRAHEDRAL_CW &&
    // RDKit❗❌:                    vol > 0)) {
    // RDKit❗❌:                 dirVect *= -1;
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           } else {
    // RDKit❗❌:             dirVect = nbr1Vect + nbr2Vect + nbr3Vect;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // we're in flatland
    // RDKit❗❌:
    // RDKit❗❌:           // github #3879 and #908: find the two neighbors with the largest
    // RDKit❗❌:           // outer angle between them and then place the H to bisect that angle
    // RDKit❗❌:           // This is recommendation ST-1.1.4 from the 2006 IUPAC "Graphical
    // RDKit❗❌:           // representation of stereochemical configuration" guideline
    // RDKit❗❌:           auto angle12 = nbr1Vect.angleTo(nbr2Vect);
    // RDKit❗❌:           auto angle13 = nbr1Vect.angleTo(nbr3Vect);
    // RDKit❗❌:           auto angle23 = nbr2Vect.angleTo(nbr3Vect);
    // RDKit❗❌:           auto accum1 = angle12 + angle13;
    // RDKit❗❌:           auto accum2 = angle12 + angle23;
    // RDKit❗❌:           auto accum3 = angle13 + angle23;
    // RDKit❗❌:           if (accum1 <= accum2 && accum1 <= accum3) {
    // RDKit❗❌:             dirVect = pickBisector(nbr1Vect, nbr2Vect, nbr3Vect);
    // RDKit❗❌:           } else if (accum2 <= accum1 && accum2 <= accum3) {
    // RDKit❗❌:             dirVect = pickBisector(nbr2Vect, nbr1Vect, nbr3Vect);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             dirVect = pickBisector(nbr3Vect, nbr1Vect, nbr2Vect);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         dirVect.normalize();
    // RDKit❗❌:         atomPos = otherPos + dirVect * ((*cfi)->is3D() ? bondLength : 1.0);
    // RDKit❗❌:         (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // RDKit❗❌:     default:
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       // FIX: figure out what to do here
    // RDKit❗❌:       // --------------------------------------------------------------------------
    // RDKit❗❌:       atomPos = otherPos + dirVect * bondLength;
    // RDKit❗❌:       for (auto cfi = mol.beginConformers(); cfi != mol.endConformers();
    // RDKit❗❌:            ++cfi) {
    // RDKit❗❌:         (*cfi)->setAtomPos(idx, atomPos);
    // RDKit❗❌:       }
    // RDKit❗❌:       break;
    // END RDKIT CPP FUNCTION MolOps::setTerminalAtomCoords degree three/four/default
    let active = active_neighbors(topology, addition.parent, addition.bond);
    let parent_position = coordinates[addition.parent.index()];
    let initialized_position = coordinates[addition.atom.index()];
    let distance = if is_3d {
        rdkit_rb0(1) + rdkit_rb0(topology.atoms[addition.parent.index()].atomic_number())
    } else {
        1.0
    };
    let other_neighbors = active
        .into_iter()
        .filter(|neighbor| *neighbor != addition.atom)
        .collect::<Vec<_>>();
    let placement_error = |reason| HydrogenError::CoordinatePlacement {
        addition: addition_index,
        conformer: conformer_id,
        dimension,
        reason,
    };

    match other_neighbors.len() {
        2 => {
            let toward_parent_one = parent_position.minus(coordinates[other_neighbors[0].index()]);
            let toward_parent_two = parent_position.minus(coordinates[other_neighbors[1].index()]);
            if toward_parent_one.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE
                || toward_parent_two.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE
            {
                return Ok(parent_position);
            }
            let neighbor_one = toward_parent_one.normalized().ok_or_else(|| {
                placement_error("degree-three first neighbor vector is not finite")
            })?;
            let neighbor_two = toward_parent_two.normalized().ok_or_else(|| {
                placement_error("degree-three second neighbor vector is not finite")
            })?;
            let mut direction = neighbor_one.plus(neighbor_two);
            if direction.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                return Ok(initialized_position);
            }
            direction = direction
                .normalized()
                .ok_or_else(|| placement_error("degree-three direction is not finite"))?;
            if is_3d
                && topology.atoms[addition.parent.index()].hybridization() == Hybridization::Sp3
            {
                let rotation_axis = neighbor_one
                    .cross(neighbor_two)
                    .cross(direction)
                    .normalized()
                    .ok_or_else(|| placement_error("degree-three SP3 rotation axis is zero"))?;
                direction = direction
                    .rotated((109.471_f64 / 2.0).to_radians(), rotation_axis)
                    .ok_or_else(|| placement_error("degree-three SP3 rotation axis is zero"))?;
            }
            Ok(parent_position.plus(direction.scaled(distance)))
        }
        3 => {
            let mut neighbors = [AddHsPoint3D::ZERO; 3];
            for (slot, neighbor) in neighbors.iter_mut().zip(other_neighbors) {
                let vector = parent_position.minus(coordinates[neighbor.index()]);
                if vector.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                    return Ok(parent_position);
                }
                *slot = vector
                    .normalized()
                    .ok_or_else(|| placement_error("degree-four neighbor vector is not finite"))?;
            }
            let [neighbor_one, neighbor_two, neighbor_three] = neighbors;
            let mut direction = if is_3d {
                if neighbor_three.dot(neighbor_one.cross(neighbor_two)).abs() < 0.1 {
                    let mut planar_direction = neighbor_one.cross(neighbor_two);
                    if planar_direction.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                        planar_direction = neighbor_one.cross(neighbor_three).scaled(-1.0);
                    }
                    if planar_direction.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                        planar_direction = neighbor_two.cross(neighbor_three);
                    }
                    if planar_direction.length_sq() < ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE {
                        return Ok(initialized_position);
                    }
                    if topology.atoms[addition.parent.index()]
                        .prop("_CIPCode")
                        .is_some()
                    {
                        let v1 = planar_direction.minus(neighbor_three);
                        let v2 = neighbor_one.minus(neighbor_three);
                        let v3 = neighbor_two.minus(neighbor_three);
                        let volume = v1.dot(v2.cross(v3));
                        let tag = topology.atoms[addition.parent.index()].chiral_tag();
                        if (tag == ChiralTag::TetrahedralCcw && volume < 0.0)
                            || (tag == ChiralTag::TetrahedralCw && volume > 0.0)
                        {
                            planar_direction = planar_direction.scaled(-1.0);
                        }
                    }
                    planar_direction
                } else {
                    neighbor_one.plus(neighbor_two).plus(neighbor_three)
                }
            } else {
                let angle_12 = neighbor_one
                    .angle_to(neighbor_two)
                    .ok_or_else(|| placement_error("degree-four angle 1-2 is not finite"))?;
                let angle_13 = neighbor_one
                    .angle_to(neighbor_three)
                    .ok_or_else(|| placement_error("degree-four angle 1-3 is not finite"))?;
                let angle_23 = neighbor_two
                    .angle_to(neighbor_three)
                    .ok_or_else(|| placement_error("degree-four angle 2-3 is not finite"))?;
                let accumulated_one = angle_12 + angle_13;
                let accumulated_two = angle_12 + angle_23;
                let accumulated_three = angle_13 + angle_23;
                if accumulated_one <= accumulated_two && accumulated_one <= accumulated_three {
                    pick_bisector(neighbor_one, neighbor_two, neighbor_three)
                } else if accumulated_two <= accumulated_one && accumulated_two <= accumulated_three
                {
                    pick_bisector(neighbor_two, neighbor_one, neighbor_three)
                } else {
                    pick_bisector(neighbor_three, neighbor_one, neighbor_two)
                }
            };
            direction = direction
                .normalized()
                .ok_or_else(|| placement_error("degree-four direction is zero or not finite"))?;
            Ok(parent_position.plus(direction.scaled(distance)))
        }
        _ => Ok(AddHsPoint3D::ZERO),
    }
}

fn validate_addition_plan(
    result: &AddHydrogensTopologyResult,
    coordinates: &CoordinateBlock,
) -> Result<usize, HydrogenError> {
    result.topology.validate()?;
    let old_atom_count = result.mapping.atoms().old_to_new().len();
    let old_bond_count = result.mapping.bonds().old_to_new().len();
    let addition_count = result.additions.len();
    if result.topology.atoms.len() != old_atom_count + addition_count {
        return Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "new atom count does not equal old atoms plus additions",
        });
    }
    if result.topology.bonds.len() != old_bond_count + addition_count {
        return Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "new bond count does not equal old bonds plus additions",
        });
    }
    result.mapping.validate_for_counts(
        old_atom_count,
        result.topology.atoms.len(),
        old_bond_count,
        result.topology.bonds.len(),
    )?;
    if result.mapping
        != TopologyMapping::with_appended(
            old_atom_count,
            old_bond_count,
            addition_count,
            addition_count,
        )
    {
        return Err(HydrogenError::InvalidAdditionPlan {
            addition: None,
            reason: "mapping is not the canonical append-only mapping",
        });
    }
    coordinates.validate_for_atom_count(old_atom_count)?;

    for (offset, addition) in result.additions.iter().enumerate() {
        if addition.atom != AtomId::new(old_atom_count + offset) {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "atom id is not the next appended row",
            });
        }
        if addition.bond != BondId::new(old_bond_count + offset) {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "bond id is not the next appended row",
            });
        }
        if addition.parent.index() >= old_atom_count {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "parent is not an original atom",
            });
        }
        let atom = &result.topology.atoms[addition.atom.index()];
        if atom.element() != Element::H {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "appended atom is not hydrogen",
            });
        }
        if atom.implicit_hydrogen() != (addition.kind == AddedHydrogenKind::Implicit) {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "hydrogen implicit marker disagrees with addition kind",
            });
        }
        let bond = &result.topology.bonds[addition.bond.index()];
        if bond.begin() != addition.parent
            || bond.end() != addition.atom
            || bond.order() != BondOrder::Single
        {
            return Err(HydrogenError::InvalidAdditionPlan {
                addition: Some(offset),
                reason: "appended bond is not the ordered parent-hydrogen single bond",
            });
        }
    }
    Ok(old_atom_count)
}

pub(crate) fn grow_coordinate_rows(
    mut coordinates: CoordinateBlock,
    count: usize,
) -> CoordinateBlock {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs conformer preparation
    // RDKit✔️✔️: unsigned int nSize = mol.getNumAtoms() + numAddHyds;
    // RDKit✔️✔️: // loop over the conformations of the molecule and allocate new space
    // RDKit✔️✔️: // for the H locations (need to do this even if we aren't adding coords so
    // RDKit✔️✔️: // that the conformers have the correct number of atoms).
    // RDKit✔️✔️: for (auto cfi = mol.beginConformers(); cfi != mol.endConformers(); ++cfi) {
    // RDKit✔️✔️:   (*cfi)->reserve(nSize);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::addHs conformer preparation
    // RDKit's subsequent RWMol::addAtom appends Point3D(0,0,0) to every
    // conformer. The detached representation grows each dimension once here;
    // this is O(conformers * additions), preserves all old rows and performs
    // the same required allocation without repeated vector reallocation.
    for conformer in &mut coordinates.conformers_2d {
        let mut values = conformer.coordinates().to_vec();
        values.resize(values.len() + count, [0.0, 0.0]);
        let mut replacement = Conformer2D::new(conformer.id(), values);
        for (key, value) in conformer.props() {
            replacement = replacement.with_prop(key.clone(), value.clone());
        }
        *conformer = replacement;
    }
    for conformer in &mut coordinates.conformers_3d {
        let mut values = conformer.coordinates().to_vec();
        values.resize(values.len() + count, [0.0, 0.0, 0.0]);
        let mut replacement = Conformer3D::new(conformer.id(), values, conformer.is_3d());
        for (key, value) in conformer.props() {
            replacement = replacement.with_prop(key.clone(), value.clone());
        }
        *conformer = replacement;
    }
    coordinates
}

const ADD_HS_SQ_DISTANCE_ZERO_TOLERANCE: f64 = 1.0e-4;

#[derive(Debug, Clone, Copy, PartialEq)]
struct AddHsPoint3D {
    x: f64,
    y: f64,
    z: f64,
}

impl AddHsPoint3D {
    const ZERO: Self = Self::new(0.0, 0.0, 0.0);

    const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    const fn from_2d(value: [f64; 2]) -> Self {
        Self::new(value[0], value[1], 0.0)
    }

    const fn from_3d(value: [f64; 3]) -> Self {
        Self::new(value[0], value[1], value[2])
    }

    const fn to_2d(self) -> [f64; 2] {
        [self.x, self.y]
    }

    const fn to_3d(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }

    const fn plus(self, other: Self) -> Self {
        // RDKit✔️✔️: constexpr Point3D &operator+=(const Point3D &other) {
        // RDKit✔️✔️:   x += other.x;
        // RDKit✔️✔️:   y += other.y;
        // RDKit✔️✔️:   z += other.z;
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        Self::new(self.x + other.x, self.y + other.y, self.z + other.z)
    }

    const fn minus(self, other: Self) -> Self {
        // RDKit✔️✔️: constexpr Point3D &operator-=(const Point3D &other) {
        // RDKit✔️✔️:   x -= other.x;
        // RDKit✔️✔️:   y -= other.y;
        // RDKit✔️✔️:   z -= other.z;
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        Self::new(self.x - other.x, self.y - other.y, self.z - other.z)
    }

    const fn scaled(self, scale: f64) -> Self {
        // RDKit✔️✔️: constexpr Point3D &operator*=(double scale) {
        // RDKit✔️✔️:   x *= scale;
        // RDKit✔️✔️:   y *= scale;
        // RDKit✔️✔️:   z *= scale;
        // RDKit✔️✔️:   return *this;
        // RDKit✔️✔️: }
        Self::new(self.x * scale, self.y * scale, self.z * scale)
    }

    const fn length_sq(self) -> f64 {
        // RDKit✔️✔️: constexpr double lengthSq() const override {
        // RDKit✔️✔️:   double res = x * x + y * y + z * z;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        self.x * self.x + self.y * self.y + self.z * self.z
    }

    const fn dot(self, other: Self) -> f64 {
        // RDKit✔️✔️: constexpr double dotProduct(const Point3D &other) const {
        // RDKit✔️✔️:   double res = x * (other.x) + y * (other.y) + z * (other.z);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        self.x * other.x + self.y * other.y + self.z * other.z
    }

    fn angle_to(self, other: Self) -> Option<f64> {
        // BEGIN RDKIT CPP FUNCTION Point3D::angleTo
        // RDKit✔️✔️: double angleTo(const Point3D &other) const {
        // RDKit✔️✔️:   double lsq = lengthSq() * other.lengthSq();
        // RDKit✔️✔️:   double dotProd = dotProduct(other);
        // RDKit✔️✔️:   dotProd /= sqrt(lsq);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // watch for roundoff error:
        // RDKit✔️✔️:   if (dotProd <= -1.0) {
        // RDKit✔️✔️:     return M_PI;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (dotProd >= 1.0) {
        // RDKit✔️✔️:     return 0.0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return acos(dotProd);
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Point3D::angleTo
        let length_product = (self.length_sq() * other.length_sq()).sqrt();
        if !length_product.is_finite() || length_product < f64::EPSILON {
            return None;
        }
        let dot = self.dot(other) / length_product;
        if !dot.is_finite() {
            None
        } else if dot <= -1.0 {
            Some(std::f64::consts::PI)
        } else if dot >= 1.0 {
            Some(0.0)
        } else {
            Some(dot.acos())
        }
    }

    const fn cross(self, other: Self) -> Self {
        // RDKit✔️✔️: constexpr Point3D crossProduct(const Point3D &other) const {
        // RDKit✔️✔️:   Point3D res;
        // RDKit✔️✔️:   res.x = y * (other.z) - z * (other.y);
        // RDKit✔️✔️:   res.y = -x * (other.z) + z * (other.x);
        // RDKit✔️✔️:   res.z = x * (other.y) - y * (other.x);
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        Self::new(
            self.y * other.z - self.z * other.y,
            -self.x * other.z + self.z * other.x,
            self.x * other.y - self.y * other.x,
        )
    }

    fn normalized(self) -> Option<Self> {
        // RDKit✔️✔️: constexpr void normalize() override {
        // RDKit✔️✔️:   double l = this->length();
        // RDKit✔️✔️:   if (l < zero_tolerance) {
        // RDKit✔️✔️:     throw std::runtime_error("Cannot normalize a zero length vector");
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   x /= l;
        // RDKit✔️✔️:   y /= l;
        // RDKit✔️✔️:   z /= l;
        // RDKit✔️✔️: }
        let length = self.length_sq().sqrt();
        (length >= f64::EPSILON).then(|| self.scaled(1.0 / length))
    }

    fn perpendicular(self) -> Option<Self> {
        // BEGIN RDKIT CPP FUNCTION Point3D::getPerpendicular
        // RDKit✔️✔️: Point3D getPerpendicular() const {
        // RDKit✔️✔️:   Point3D res(0.0, 0.0, 0.0);
        // RDKit✔️✔️:   if (x) {
        // RDKit✔️✔️:     if (y) {
        // RDKit✔️✔️:       res.y = -1 * x;
        // RDKit✔️✔️:       res.x = y;
        // RDKit✔️✔️:     } else if (z) {
        // RDKit✔️✔️:       res.z = -1 * x;
        // RDKit✔️✔️:       res.x = z;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       res.y = 1;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (y) {
        // RDKit✔️✔️:     if (z) {
        // RDKit✔️✔️:       res.z = -1 * y;
        // RDKit✔️✔️:       res.y = z;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       res.x = 1;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (z) {
        // RDKit✔️✔️:     res.x = 1;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   double l = res.length();
        // RDKit✔️✔️:   res /= l;
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Point3D::getPerpendicular
        let result = if self.x != 0.0 {
            if self.y != 0.0 {
                Self::new(self.y, -self.x, 0.0)
            } else if self.z != 0.0 {
                Self::new(self.z, 0.0, -self.x)
            } else {
                Self::new(0.0, 1.0, 0.0)
            }
        } else if self.y != 0.0 {
            if self.z != 0.0 {
                Self::new(0.0, self.z, -self.y)
            } else {
                Self::new(1.0, 0.0, 0.0)
            }
        } else if self.z != 0.0 {
            Self::new(1.0, 0.0, 0.0)
        } else {
            Self::ZERO
        };
        result.normalized()
    }

    fn rotated(self, angle: f64, axis: Self) -> Option<Self> {
        // BEGIN RDKIT CPP FUNCTION Transform3D::SetRotation
        // RDKit✔️✔️: void Transform3D::SetRotation(double cosT, double sinT, const Point3D &axis) {
        // RDKit✔️✔️:   double t = 1 - cosT;
        // RDKit✔️✔️:   double X = axis.x;
        // RDKit✔️✔️:   double Y = axis.y;
        // RDKit✔️✔️:   double Z = axis.z;
        // RDKit✔️✔️:   double *data = d_data.get();
        // RDKit✔️✔️:   data[0] = t * X * X + cosT;
        // RDKit✔️✔️:   data[1] = t * X * Y - sinT * Z;
        // RDKit✔️✔️:   data[2] = t * X * Z + sinT * Y;
        // RDKit✔️✔️:   data[4] = t * X * Y + sinT * Z;
        // RDKit✔️✔️:   data[5] = t * Y * Y + cosT;
        // RDKit✔️✔️:   data[6] = t * Y * Z - sinT * X;
        // RDKit✔️✔️:   data[8] = t * X * Z - sinT * Y;
        // RDKit✔️✔️:   data[9] = t * Y * Z + sinT * X;
        // RDKit✔️✔️:   data[10] = t * Z * Z + cosT;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: void Transform3D::SetRotation(double angle, const Point3D &axis) {
        // RDKit✔️✔️:   this->setToIdentity();
        // RDKit✔️✔️:   double c = cos(angle);
        // RDKit✔️✔️:   double s = sin(angle);
        // RDKit✔️✔️:   this->SetRotation(c, s, axis);
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Transform3D::SetRotation
        let axis = axis.normalized()?;
        let cosine = angle.cos();
        let sine = angle.sin();
        let t = 1.0 - cosine;
        Some(Self::new(
            (t * axis.x * axis.x + cosine) * self.x
                + (t * axis.x * axis.y - sine * axis.z) * self.y
                + (t * axis.x * axis.z + sine * axis.y) * self.z,
            (t * axis.x * axis.y + sine * axis.z) * self.x
                + (t * axis.y * axis.y + cosine) * self.y
                + (t * axis.y * axis.z - sine * axis.x) * self.z,
            (t * axis.x * axis.z - sine * axis.y) * self.x
                + (t * axis.y * axis.z + sine * axis.x) * self.y
                + (t * axis.z * axis.z + cosine) * self.z,
        ))
    }
}

fn append_hydrogen(
    topology: &mut TopologyBlock,
    additions: &mut Vec<AddedHydrogen>,
    parent: AtomId,
    kind: AddedHydrogenKind,
) {
    let atom = AtomId::new(topology.atoms.len());
    let atom_spec =
        AtomSpec::new(Element::H).with_implicit_hydrogen(kind == AddedHydrogenKind::Implicit);
    topology
        .atoms
        .push(cosmolkit_model::Atom::from_spec(atom, atom_spec));
    let bond = BondId::new(topology.bonds.len());
    topology.bonds.push(Bond::from_spec(
        bond,
        BondSpec::new(parent, atom, BondOrder::Single),
    ));
    additions.push(AddedHydrogen {
        atom,
        bond,
        parent,
        kind,
    });
}

fn is_query_atom(
    topology: &TopologyBlock,
    atom: AtomId,
    query_state: Option<QueryStateRef<'_>>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION isQueryAtom
    // RDKit✔️✔️: bool isQueryAtom(const RWMol &mol, const Atom &atom) {
    // RDKit✔️✔️:   if (atom.hasQuery()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto bnd : mol.atomBonds(&atom)) {
    // RDKit✔️✔️:     if (bnd->hasQuery()) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isQueryAtom
    // Behavior review: a typed Explicit row is the source query dynamic type;
    // carrier-derived rows and the None path are ordinary source values.
    // Complexity review: one O(1) atom test plus the source-shaped incident
    // bond scan, with O(1) identity checks and no allocation.
    query_state.is_some_and(|state| {
        state.atom_has_query(atom)
            || topology
                .adjacency
                .neighbors_of(atom.index())
                .iter()
                .any(|neighbor| state.bond_has_query(neighbor.bond))
    })
}

/// Apply the default RemoveHs operation to detached blocks.
pub fn remove_hydrogens_impl(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
) -> Result<RemoveHydrogensResult, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs copy overload
    // RDKit✔️❌: ROMol *removeHs(const ROMol &mol, const RemoveHsParameters &ps, bool sanitize) {
    // RDKit✔️❌:   auto *res = new RWMol(mol);
    // RDKit✔️❌:   try {
    // RDKit✔️❌:     removeHs(*res, ps, sanitize);
    // RDKit✔️❌:   } catch (const MolSanitizeException &) {
    // RDKit✔️❌:     delete res;
    // RDKit✔️❌:     throw;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return static_cast<ROMol *>(res);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::removeHs copy overload
    // Owned detached blocks provide the same failure atomicity without a
    // heap-allocated mutable molecule, while the final mapping and valence
    // assignment require additional linear output storage.
    remove_hydrogens_with_params(
        topology,
        coordinates,
        properties,
        &RemoveHsParams::default(),
    )
}

/// Apply the detached RemoveHs candidate selection and topology compaction.
pub fn remove_hydrogens_with_params(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &RemoveHsParams,
) -> Result<RemoveHydrogensResult, HydrogenError> {
    remove_hydrogens_with_query_state(topology, coordinates, properties, params, None)
}

/// Internal typed-query RemoveHs owner.
#[doc(hidden)]
pub fn remove_hydrogens_with_query_state(
    mut topology: TopologyBlock,
    mut coordinates: CoordinateBlock,
    mut properties: MoleculeProperties,
    params: &RemoveHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<RemoveHydrogensResult, HydrogenError> {
    validate_blocks(&topology, &coordinates)?;
    validate_property_lists(&properties, topology.atoms.len(), topology.bonds.len())?;
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }
    let original_atom_count = topology.atoms.len();
    let original_bond_count = topology.bonds.len();
    let mut mapping = TopologyMapping::identity(original_atom_count, original_bond_count);
    // Materialize transport rows with current carriers, not the overlay's
    // potentially stale snapshots. Predicate trees and origins are preserved.
    let mut query_rows = query_state
        .map(|state| remap_query_rows(state, &topology, &mapping))
        .transpose()?;
    let mut warnings = Vec::new();

    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs preliminary isotope pass
    // RDKit✔️❌: void removeHs(RWMol &mol, const RemoveHsParameters &ps, bool sanitize) {
    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit✔️❌:     // if there are any non-isotopic Hs remove them first
    // RDKit✔️❌:     // to make sure chirality is preserved
    // RDKit✔️❌:     bool needRemoveHs = false;
    // RDKit✔️❌:     for (auto atom : mol.atoms()) {
    // RDKit✔️❌:       if (atom->getAtomicNum() == 1 && atom->getIsotope() == 0) {
    // RDKit✔️❌:         needRemoveHs = true;
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (needRemoveHs) {
    // RDKit✔️❌:       RemoveHsParameters psCopy(ps);
    // RDKit✔️❌:       psCopy.removeAndTrackIsotopes = false;
    // RDKit✔️❌:       psCopy.removeIsotopes = false;
    // RDKit✔️❌:       removeHs(mol, psCopy, false);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION MolOps::removeHs preliminary isotope pass
    // Each detached pass uses one batch compaction and mapping. Composing the
    // two mappings retains source ordering but allocates linear index vectors
    // that the mutable source molecule does not expose.
    if params.remove_and_track_isotopes
        && topology
            .atoms
            .iter()
            .any(|atom| atom.atomic_number() == 1 && atom.isotope().unwrap_or(0) == 0)
    {
        let mut preliminary_params = params.clone();
        preliminary_params.remove_and_track_isotopes = false;
        preliminary_params.remove_isotopes = false;
        preliminary_params.sanitize = false;
        let intermediate_atom_count = topology.atoms.len();
        let intermediate_bond_count = topology.bonds.len();
        let pass_state = query_rows
            .as_ref()
            .map(|(atoms, bonds)| QueryStateRef::try_for_topology(atoms, bonds, &topology))
            .transpose()?;
        let preliminary =
            remove_hydrogens_pass(topology, &mut properties, &preliminary_params, pass_state)?;
        mapping = compose_topology_mappings(
            &mapping,
            &preliminary.mapping,
            original_atom_count,
            intermediate_atom_count,
            preliminary.topology.atoms.len(),
            original_bond_count,
            intermediate_bond_count,
            preliminary.topology.bonds.len(),
        )?;
        topology = preliminary.topology;
        query_rows = preliminary.query_rows;
        warnings.extend(preliminary.warnings);
    }

    let intermediate_atom_count = topology.atoms.len();
    let intermediate_bond_count = topology.bonds.len();
    let pass_state = query_rows
        .as_ref()
        .map(|(atoms, bonds)| QueryStateRef::try_for_topology(atoms, bonds, &topology))
        .transpose()?;
    let final_pass = remove_hydrogens_pass(topology, &mut properties, params, pass_state)?;
    mapping = compose_topology_mappings(
        &mapping,
        &final_pass.mapping,
        original_atom_count,
        intermediate_atom_count,
        final_pass.topology.atoms.len(),
        original_bond_count,
        intermediate_bond_count,
        final_pass.topology.bonds.len(),
    )?;
    topology = final_pass.topology;
    query_rows = final_pass.query_rows;
    warnings.extend(final_pass.warnings);
    // The sole public construction publishes ONLY the final pass's ring
    // state by MOVE; the preliminary isotope pass hardcodes sanitize=false
    // and never reaches this site with intermediate rows.
    let final_rings = final_pass.final_rings;
    let final_valence = Some(final_pass.final_valence);
    #[cfg(test)]
    remove_hs_ring_probe::record_publication(&final_rings);

    if let Some((atoms, bonds)) = query_rows.as_ref() {
        QueryStateRef::try_for_topology(atoms, bonds, &topology)?;
    }

    coordinates.remap_topology(&mapping.retained_atom_indices());
    properties.remap_topology(mapping.atoms().new_to_old(), mapping.bonds().new_to_old());
    mapping.validate_for_counts(
        original_atom_count,
        topology.atoms.len(),
        original_bond_count,
        topology.bonds.len(),
    )?;
    validate_blocks(&topology, &coordinates)?;
    validate_property_lists(&properties, topology.atoms.len(), topology.bonds.len())?;
    // The final pass moves its actual source scalar state here. No extra
    // final-topology calculation normalizes away pre-removal cached values.
    Ok(RemoveHydrogensResult {
        topology,
        coordinates,
        properties,
        mapping,
        final_valence,
        final_rings,
        warnings,
    })
}

struct RemoveHydrogensPassResult {
    topology: TopologyBlock,
    final_valence: ValenceAssignment,
    mapping: TopologyMapping,
    warnings: Vec<HydrogenWarning>,
    query_rows: Option<(Vec<QueryAtom>, Vec<QueryBond>)>,
    /// Final ring-state transport for this pass. `None` is the exact final
    /// source-uninitialized state (removeHs ALWAYS clearComputedProps(true),
    /// so the no-sanitize branches report the reset state); `Some` is the
    /// exact SAN initialized state, INCLUDING zero rings. Moved, never
    /// cloned; no finder, stale-row remap, supplied-input rings parameter,
    /// or final quality upgrade is added.
    final_rings: Option<RingInfo>,
}

/// Test-only observation points at the ACTUAL RemoveHs ring-transport sites.
/// Counters are thread-local and never reset; observations are recorded at
/// the real sanitize return, the real pass return and the outer final-pass
/// publication. Production builds contain no probes.
#[cfg(test)]
pub(crate) mod remove_hs_ring_probe {
    use std::cell::{Cell, RefCell};

    use crate::RingInfo;

    thread_local! {
        static SANITIZE_RETURNS: Cell<u64> = const { Cell::new(0) };
        static PASS_RETURNS: Cell<u64> = const { Cell::new(0) };
        static LAST_SANITIZE_ATOM_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_SANITIZE_BOND_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_SANITIZE_STATE: RefCell<Option<RingInfo>> = const { RefCell::new(None) };
        static LAST_PASS_ATOM_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_PASS_BOND_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_PASS_STATE: RefCell<Option<RingInfo>> = const { RefCell::new(None) };
    }

    fn pointers(rings: &Option<RingInfo>) -> (usize, usize) {
        rings
            .as_ref()
            .map(|state| {
                (
                    state.atom_rings().as_ptr() as usize,
                    state.bond_rings().as_ptr() as usize,
                )
            })
            .unwrap_or((0, 0))
    }

    pub(crate) fn record_sanitize_return(rings: &Option<RingInfo>) {
        SANITIZE_RETURNS.with(|calls| calls.set(calls.get() + 1));
        let (atom_ptr, bond_ptr) = pointers(rings);
        LAST_SANITIZE_ATOM_PTR.with(|cell| cell.set(atom_ptr));
        LAST_SANITIZE_BOND_PTR.with(|cell| cell.set(bond_ptr));
        LAST_SANITIZE_STATE.with(|cell| *cell.borrow_mut() = rings.clone());
    }

    pub(crate) fn record_pass_return(rings: &Option<RingInfo>) {
        PASS_RETURNS.with(|calls| calls.set(calls.get() + 1));
        let (atom_ptr, bond_ptr) = pointers(rings);
        LAST_PASS_ATOM_PTR.with(|cell| cell.set(atom_ptr));
        LAST_PASS_BOND_PTR.with(|cell| cell.set(bond_ptr));
        LAST_PASS_STATE.with(|cell| *cell.borrow_mut() = rings.clone());
    }

    pub(crate) fn sanitize_returns() -> u64 {
        SANITIZE_RETURNS.with(Cell::get)
    }

    pub(crate) fn pass_returns() -> u64 {
        PASS_RETURNS.with(Cell::get)
    }

    pub(crate) fn last_sanitize_pointers() -> (usize, usize) {
        (
            LAST_SANITIZE_ATOM_PTR.with(Cell::get),
            LAST_SANITIZE_BOND_PTR.with(Cell::get),
        )
    }

    pub(crate) fn last_sanitize_state() -> Option<RingInfo> {
        LAST_SANITIZE_STATE.with(|cell| cell.borrow().clone())
    }

    pub(crate) fn last_pass_pointers() -> (usize, usize) {
        (
            LAST_PASS_ATOM_PTR.with(Cell::get),
            LAST_PASS_BOND_PTR.with(Cell::get),
        )
    }

    pub(crate) fn last_pass_state() -> Option<RingInfo> {
        LAST_PASS_STATE.with(|cell| cell.borrow().clone())
    }

    thread_local! {
        static PUBLICATIONS: Cell<u64> = const { Cell::new(0) };
        static LAST_PUBLICATION_ATOM_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_PUBLICATION_BOND_PTR: Cell<usize> = const { Cell::new(0) };
        static LAST_PUBLICATION_STATE: RefCell<Option<RingInfo>> = const { RefCell::new(None) };
    }

    /// Actual outer final-pass publication site.
    pub(crate) fn record_publication(rings: &Option<RingInfo>) {
        PUBLICATIONS.with(|calls| calls.set(calls.get() + 1));
        let (atom_ptr, bond_ptr) = pointers(rings);
        LAST_PUBLICATION_ATOM_PTR.with(|cell| cell.set(atom_ptr));
        LAST_PUBLICATION_BOND_PTR.with(|cell| cell.set(bond_ptr));
        LAST_PUBLICATION_STATE.with(|cell| *cell.borrow_mut() = rings.clone());
    }

    pub(crate) fn publications() -> u64 {
        PUBLICATIONS.with(Cell::get)
    }

    pub(crate) fn last_publication_pointers() -> (usize, usize) {
        (
            LAST_PUBLICATION_ATOM_PTR.with(Cell::get),
            LAST_PUBLICATION_BOND_PTR.with(Cell::get),
        )
    }

    pub(crate) fn last_publication_state() -> Option<RingInfo> {
        LAST_PUBLICATION_STATE.with(|cell| cell.borrow().clone())
    }

    thread_local! {
        static PASS_HISTORY: RefCell<Vec<(bool, Option<RingInfo>)>> =
            const { RefCell::new(Vec::new()) };
    }

    /// Append-only ordered pass-return history at the ACTUAL pass-return
    /// site: (params.sanitize at that pass, its returned final_rings).
    /// Never reset or removed; clones only under cfg(test).
    pub(crate) fn record_pass_history(sanitize: bool, rings: &Option<RingInfo>) {
        PASS_HISTORY.with(|history| history.borrow_mut().push((sanitize, rings.clone())));
    }

    pub(crate) fn pass_history_len() -> usize {
        PASS_HISTORY.with(|history| history.borrow().len())
    }

    /// The entries appended after `baseline`, in order.
    pub(crate) fn pass_history_after(baseline: usize) -> Vec<(bool, Option<RingInfo>)> {
        PASS_HISTORY.with(|history| {
            let history = history.borrow();
            history[baseline.min(history.len())..].to_vec()
        })
    }
}

// R1: the sixteen-call private-pass ring-transport product. Sanitizer
// invocation1 IFF the pass's ORIGINAL atom count is nonzero AND
// remove_nonimplicit AND sanitize (the bitset-size guard); counters are
// thread-local, never reset — all observations are per-call deltas. H atoms
// are NON-implicit (explicit model atoms, no isImplicit marker).
#[cfg(test)]
fn fixed_property_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("fixed fixture text is UTF8")
}

#[cfg(test)]
mod remove_hs_ring_pass_tests {
    use super::RemoveHsParams;
    use super::remove_hs_ring_probe as probe;
    use super::remove_hydrogens_pass;
    use crate::RingFindType;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn carbon(id: usize, aromatic: bool) -> Atom {
        let mut spec = AtomSpec::new(Element::C);
        if aromatic {
            spec = spec.with_aromatic(true);
        }
        Atom::from_spec(AtomId::new(id), spec)
    }

    fn hydrogen(id: usize) -> Atom {
        let atom = Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::H));
        assert_eq!(atom.atomic_number(), 1);
        assert!(atom.isotope().is_none(), "fixture H must be non-isotopic");
        atom
    }

    fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn single(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true),
        )
    }

    fn g0() -> TopologyBlock {
        topology(Vec::new(), Vec::new())
    }

    fn g1() -> TopologyBlock {
        topology(
            vec![carbon(0, false), carbon(1, false)],
            vec![single(0, 0, 1)],
        )
    }

    fn g2() -> TopologyBlock {
        topology(
            (0..6).map(|id| carbon(id, true)).collect(),
            (0..6)
                .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                .collect(),
        )
    }

    fn g3() -> TopologyBlock {
        topology(
            std::iter::once(hydrogen(0))
                .chain((1..7).map(|id| carbon(id, true)))
                .collect(),
            std::iter::once(single(0, 0, 1))
                .chain((1..7).map(|id| aromatic_bond(id, id, (id % 6) + 1)))
                .collect(),
        )
    }

    fn params(sanitize: bool, nonimplicit: bool) -> RemoveHsParams {
        let mut params = RemoveHsParams::default();
        params.sanitize = sanitize;
        params.remove_nonimplicit = nonimplicit;
        params.update_explicit_count = true;
        params.show_warnings = false;
        params
    }

    #[test]
    fn remove_hydrogens_ring_pass_sixteen_call_product() {
        let graphs: [(&str, fn() -> TopologyBlock, usize, usize); 4] = [
            ("G0", g0 as fn() -> TopologyBlock, 0, 0),
            ("G1", g1, 2, 1),
            ("G2", g2, 6, 6),
            ("G3", g3, 7, 7),
        ];
        let mut calls = 0usize;
        for (name, build, atom_count, bond_count) in graphs {
            for sanitize in [false, true] {
                for nonimplicit in [false, true] {
                    let label = format!("{name}/s{sanitize}/ni{nonimplicit}");
                    let input = build();
                    // Non-implicit H prerequisite for the removal product:
                    // atomic number 1, non-isotopic, implicit_hydrogen false.
                    for atom in &input.atoms {
                        if atom.atomic_number() == 1 {
                            assert!(atom.isotope().is_none(), "{label}: isotopic H");
                            assert!(
                                !atom.implicit_hydrogen(),
                                "{label}: implicit-H marker present"
                            );
                        }
                    }
                    // Fresh full retained baseline; the owner takes a value,
                    // so the retained original is proven by passing a clone.
                    let retained = input.clone();
                    let mut properties = MoleculeProperties::default();
                    properties.set_prop("ck-ordinary", "kept").unwrap();
                    properties
                        .set_computed_prop("ck-computed", "cleared")
                        .unwrap();
                    let properties_before = properties.clone();
                    let sanitize_before = probe::sanitize_returns();
                    let pass_before = probe::pass_returns();
                    let result = remove_hydrogens_pass(
                        input,
                        &mut properties,
                        &params(sanitize, nonimplicit),
                        None,
                    )
                    .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                    calls += 1;
                    // Retained ORIGINAL graph unchanged.
                    assert_eq!(retained, build(), "{label}: retained graph mutated");
                    // Relative, non-reset invocation counts.
                    let expected_sanitizations =
                        usize::from(atom_count > 0 && nonimplicit && sanitize);
                    assert_eq!(
                        probe::sanitize_returns() - sanitize_before,
                        expected_sanitizations as u64,
                        "{label}: sanitizer invocations"
                    );
                    assert_eq!(
                        probe::pass_returns() - pass_before,
                        1,
                        "{label}: pass returns"
                    );
                    // Ordinary property retained; computed property cleared by
                    // the actual clearComputedProps(true) of this pass.
                    assert_eq!(
                        properties.prop("ck-ordinary"),
                        Some(&cosmolkit_model::PropertyValue::String("kept".into())),
                        "{label}"
                    );
                    assert!(
                        !properties.is_prop_computed("ck-computed").unwrap(),
                        "{label}"
                    );
                    assert!(
                        properties.prop("ck-computed").is_none(),
                        "{label}: computed value must be absent"
                    );
                    assert_eq!(
                        properties_before.prop("ck-ordinary"),
                        Some(&cosmolkit_model::PropertyValue::String("kept".into()))
                    );

                    let both = nonimplicit && sanitize && atom_count > 0;
                    match (name, both) {
                        ("G0", _) => {
                            assert!(result.final_rings.is_none(), "{label}");
                            assert!(result.mapping.atoms().new_to_old().is_empty());
                            assert!(result.mapping.bonds().new_to_old().is_empty());
                        }
                        ("G1", true) => {
                            let rings = result.final_rings.as_ref().unwrap();
                            assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                            assert!(rings.is_initialized(), "{label}");
                            assert!(rings.atom_rings().is_empty(), "{label}: rows");
                            assert!(rings.bond_rings().is_empty(), "{label}: bond rows");
                            assert_eq!(rings.atom_row_count(), 2, "{label}");
                            assert_eq!(rings.bond_row_count(), 1, "{label}");
                            assert_eq!(
                                rings.atom_members(AtomId::new(0)),
                                &[] as &[usize],
                                "{label}: m0"
                            );
                            assert_eq!(
                                rings.atom_members(AtomId::new(1)),
                                &[] as &[usize],
                                "{label}: m1"
                            );
                            assert_eq!(
                                rings.bond_members(BondId::new(0)),
                                &[] as &[usize],
                                "{label}: bm0"
                            );
                        }
                        ("G1", false) => {
                            assert!(result.final_rings.is_none(), "{label}");
                        }
                        ("G2", true) | ("G3", true) => {
                            let rings = result.final_rings.as_ref().unwrap();
                            assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                            assert_eq!(rings.atom_rings().len(), 1, "{label}: ring count");
                            let expected_atoms: Vec<usize> = (0..6).collect();
                            let mut got: Vec<usize> = rings.atom_rings()[0]
                                .iter()
                                .map(|atom| atom.index())
                                .collect();
                            got.sort_unstable();
                            assert_eq!(got, expected_atoms, "{label}: normalized ring atoms");
                            let mut bonds: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bonds.sort_unstable();
                            assert_eq!(bonds, expected_atoms, "{label}: normalized ring bonds");
                            for atom in 0..6usize {
                                assert_eq!(
                                    rings.atom_rings()[0]
                                        .iter()
                                        .filter(|a| a.index() == atom)
                                        .count(),
                                    1,
                                    "{label}: atom {atom} membership"
                                );
                            }
                            // Both membership dimensions + members for every
                            // final index; normalized paired bond row.
                            assert_eq!(rings.bond_row_count(), 6, "{label}: bond dims");
                            let mut bond_row: Vec<usize> = rings.bond_rings()[0]
                                .iter()
                                .map(|bond| bond.index())
                                .collect();
                            bond_row.sort_unstable();
                            assert_eq!(bond_row, expected_atoms, "{label}: normalized bonds");
                            for index in 0..6usize {
                                assert_eq!(
                                    rings.atom_members(AtomId::new(index)),
                                    &[0],
                                    "{label}: atom {index} members"
                                );
                            }
                            for bond in 0..6usize {
                                assert_eq!(
                                    rings.bond_members(BondId::new(bond)),
                                    &[0],
                                    "{label}: bond {bond} members"
                                );
                            }
                            if name == "G3" {
                                assert_eq!(rings.atom_row_count(), 6, "{label}");
                                assert_eq!(rings.bond_row_count(), 6, "{label}");
                            }
                        }
                        ("G2", false) | ("G3", false) => {
                            assert!(result.final_rings.is_none(), "{label}");
                        }
                        _ => unreachable!("{label}"),
                    }
                    // Mapping literals.
                    if name == "G3" && nonimplicit {
                        let atom_n2o: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .new_to_old()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        let bond_n2o: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .new_to_old()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(
                            atom_n2o,
                            vec![Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)],
                            "{label}: atom new_to_old"
                        );
                        assert_eq!(
                            bond_n2o,
                            vec![Some(1), Some(2), Some(3), Some(4), Some(5), Some(6)],
                            "{label}: bond new_to_old"
                        );
                        let atom_o2n: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .old_to_new()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        assert_eq!(
                            atom_o2n,
                            vec![None, Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
                            "{label}: atom old_to_new"
                        );
                        let bond_o2n: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .old_to_new()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(
                            bond_o2n,
                            vec![None, Some(0), Some(1), Some(2), Some(3), Some(4), Some(5)],
                            "{label}: bond old_to_new"
                        );
                    } else {
                        // Exact identity at ORIGINAL atom AND bond counts.
                        let identity_atoms: Vec<Option<usize>> =
                            (0..atom_count).map(Some).collect();
                        let identity_bonds: Vec<Option<usize>> =
                            (0..bond_count).map(Some).collect();
                        let got_atoms: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .new_to_old()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        let got_atom_old: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .old_to_new()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        let got_bonds: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .new_to_old()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        let got_bond_old: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .old_to_new()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(got_atoms, identity_atoms, "{label}: atoms new_to_old");
                        assert_eq!(got_atom_old, identity_atoms, "{label}: atoms old_to_new");
                        assert_eq!(got_bonds, identity_bonds, "{label}: bonds new_to_old");
                        assert_eq!(got_bond_old, identity_bonds, "{label}: bonds old_to_new");
                        let _ = bond_count;
                    }
                }
            }
        }
        assert_eq!(calls, 16, "exact census");
    }
}

// R3: the eight-call G4 isotope-tracking public product. tracking=true runs
// exactly ONE preliminary non-isotopic-H pass (hardcoded sanitize=false,
// observed at the ACTUAL pass-return probe) whose final_rings stays None;
// only the final pass's state is published. D2 (isotope Some(2)) is retained
// by the non-tracking removal and removed by the two-pass tracking removal,
// with the tracked payload asserted through the existing owner.
#[cfg(test)]
pub(crate) mod remove_hs_ring_isotope_tests_helpers {
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, CoordinateBlock,
        TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    pub(crate) fn carbon(id: usize) -> Atom {
        Atom::from_spec(
            AtomId::new(id),
            AtomSpec::new(Element::C).with_aromatic(true),
        )
    }

    pub(crate) fn hydrogen(id: usize, isotope: Option<u16>) -> Atom {
        let mut spec = AtomSpec::new(Element::H);
        if let Some(value) = isotope {
            spec = spec.with_isotope(value);
        }
        Atom::from_spec(AtomId::new(id), spec)
    }

    pub(crate) fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    pub(crate) fn single(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    pub(crate) fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true),
        )
    }

    pub(crate) fn g2() -> TopologyBlock {
        topology(
            (0..6).map(|id| carbon(id)).collect(),
            (0..6)
                .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                .collect(),
        )
    }

    pub(crate) fn g3() -> TopologyBlock {
        topology(
            std::iter::once(hydrogen(0, None))
                .chain((1..7).map(|id| carbon(id)))
                .collect(),
            std::iter::once(single(0, 0, 1))
                .chain((1..7).map(|id| aromatic_bond(id, id, (id % 6) + 1)))
                .collect(),
        )
    }

    pub(crate) fn g4() -> TopologyBlock {
        topology(
            vec![
                hydrogen(0, None),
                carbon(1),
                hydrogen(2, Some(2)),
                carbon(3),
                carbon(4),
                carbon(5),
                carbon(6),
                carbon(7),
            ],
            vec![
                single(0, 0, 1),
                single(1, 2, 3),
                aromatic_bond(2, 1, 3),
                aromatic_bond(3, 3, 4),
                aromatic_bond(4, 4, 5),
                aromatic_bond(5, 5, 6),
                aromatic_bond(6, 6, 7),
                aromatic_bond(7, 7, 1),
            ],
        )
    }

    pub(crate) fn coordinates(atom_count: usize) -> CoordinateBlock {
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                17,
                (0..atom_count)
                    .map(|row| [row as f64, -(row as f64)])
                    .collect(),
            )],
            conformers_3d: vec![Conformer3D::new(
                19,
                (0..atom_count)
                    .map(|row| [row as f64 + 0.5, 0.0, -0.0])
                    .collect(),
                true,
            )],
            source_conformer_order: None,
            source_coordinate_dim: None,
        }
    }

    pub(crate) fn params(
        tracking: bool,
        nonimplicit: bool,
        sanitize: bool,
    ) -> super::RemoveHsParams {
        let mut params = super::RemoveHsParams::default();
        params.remove_and_track_isotopes = tracking;
        params.remove_nonimplicit = nonimplicit;
        params.sanitize = sanitize;
        params.update_explicit_count = true;
        params.show_warnings = false;
        params
    }

    /// R5 base: two ordinary carbons a0,a1 Single b0=(0,1).
    pub(crate) fn g1_invalid_base() -> TopologyBlock {
        topology(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![single(0, 0, 1)],
        )
    }
}

#[cfg(test)]
mod remove_hs_ring_isotope_tests {
    use super::RemoveHsParams;
    use super::remove_hs_ring_probe as probe;
    use crate::{RingFindType, remove_hydrogens_with_params};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, Conformer2D, Conformer3D, CoordinateBlock,
        MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    fn carbon(id: usize) -> Atom {
        Atom::from_spec(
            AtomId::new(id),
            AtomSpec::new(Element::C).with_aromatic(true),
        )
    }

    fn hydrogen(id: usize, isotope: Option<u16>) -> Atom {
        let mut spec = AtomSpec::new(Element::H);
        if let Some(value) = isotope {
            spec = spec.with_isotope(value);
        }
        let atom = Atom::from_spec(AtomId::new(id), spec);
        assert_eq!(atom.atomic_number(), 1);
        assert_eq!(atom.isotope(), isotope);
        atom
    }

    fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    fn single(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true),
        )
    }

    /// G4: a0=H, a1=C, a2=D(isotope 2), a3..7=C (a1,a3..7 aromatic); H and D
    /// attach to DIFFERENT carbons.
    fn g4() -> TopologyBlock {
        topology(
            vec![
                hydrogen(0, None),
                carbon(1),
                hydrogen(2, Some(2)),
                carbon(3),
                carbon(4),
                carbon(5),
                carbon(6),
                carbon(7),
            ],
            vec![
                single(0, 0, 1),
                single(1, 2, 3),
                aromatic_bond(2, 1, 3),
                aromatic_bond(3, 3, 4),
                aromatic_bond(4, 4, 5),
                aromatic_bond(5, 5, 6),
                aromatic_bond(6, 6, 7),
                aromatic_bond(7, 7, 1),
            ],
        )
    }

    fn coordinates(atom_count: usize) -> CoordinateBlock {
        CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(
                17,
                (0..atom_count)
                    .map(|row| [row as f64, -(row as f64)])
                    .collect(),
            )],
            conformers_3d: vec![Conformer3D::new(
                19,
                (0..atom_count)
                    .map(|row| [row as f64 + 0.5, 0.0, -0.0])
                    .collect(),
                true,
            )],
            source_conformer_order: None,
            source_coordinate_dim: None,
        }
    }

    fn params(tracking: bool, nonimplicit: bool, sanitize: bool) -> RemoveHsParams {
        let mut params = RemoveHsParams::default();
        params.remove_and_track_isotopes = tracking;
        params.remove_nonimplicit = nonimplicit;
        params.sanitize = sanitize;
        params.update_explicit_count = true;
        params.show_warnings = false;
        params
    }

    #[test]
    fn remove_hydrogens_ring_isotope_eight_call_product() {
        let mut calls = 0usize;
        for tracking in [false, true] {
            for nonimplicit in [false, true] {
                for sanitize in [false, true] {
                    let label = format!("t{tracking}/ni{nonimplicit}/s{sanitize}");
                    let topology = g4();
                    let coordinates = coordinates(8);
                    let mut properties = MoleculeProperties::default();
                    properties.set_prop("ck-ordinary", "kept").unwrap();
                    let topology_snapshot = topology.clone();
                    let coordinates_snapshot = coordinates.clone();
                    let properties_snapshot = properties.clone();
                    // G4 prerequisites: exact 8 atoms/edges, H0/D2 on
                    // DIFFERENT carbons, isotope/implicit flags.
                    assert_eq!(topology.atoms.len(), 8, "{label}: atom count");
                    assert_eq!(topology.bonds.len(), 8, "{label}: bond count");
                    assert_eq!(topology.atoms[0].atomic_number(), 1, "{label}");
                    assert_eq!(topology.atoms[0].isotope(), None, "{label}");
                    assert!(!topology.atoms[0].implicit_hydrogen(), "{label}");
                    assert_eq!(topology.atoms[2].atomic_number(), 1, "{label}");
                    assert_eq!(topology.atoms[2].isotope(), Some(2), "{label}");
                    assert!(!topology.atoms[2].implicit_hydrogen(), "{label}");
                    assert_eq!(topology.bonds[0].begin().index(), 0, "{label}: H0-C1");
                    assert_eq!(topology.bonds[0].end().index(), 1, "{label}: H0-C1");
                    assert_eq!(topology.bonds[1].begin().index(), 2, "{label}: D2-C3");
                    assert_eq!(topology.bonds[1].end().index(), 3, "{label}: D2-C3");
                    let pass_before = probe::pass_returns();
                    let sanitize_before = probe::sanitize_returns();
                    let history_before = probe::pass_history_len();
                    let result = remove_hydrogens_with_params(
                        topology.clone(),
                        coordinates.clone(),
                        properties.clone(),
                        &params(tracking, nonimplicit, sanitize),
                    )
                    .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                    calls += 1;
                    // Retained ORIGINAL whole values unchanged.
                    assert_eq!(topology_snapshot, topology, "{label}: retained graph");
                    assert_eq!(
                        coordinates_snapshot, coordinates,
                        "{label}: retained coords"
                    );
                    assert_eq!(properties_snapshot, properties, "{label}: retained props");
                    // Actual ordered pass-return history slice: tracking
                    // adds the preliminary (sanitize=false, None) entry
                    // BEFORE the final entry, whose state equals the returned
                    // final_rings exactly.
                    let history = probe::pass_history_after(history_before);
                    let expected_entries = if tracking { 2 } else { 1 };
                    assert_eq!(history.len(), expected_entries, "{label}: history len");
                    if tracking {
                        assert_eq!(history[0].0, false, "{label}: preliminary sanitize=false");
                        assert!(history[0].1.is_none(), "{label}: preliminary None");
                        assert_eq!(history[1].0, sanitize, "{label}: final sanitize flag");
                        assert_eq!(history[1].1, result.final_rings, "{label}: final state");
                    } else {
                        assert_eq!(history[0].0, sanitize, "{label}: only-pass sanitize flag");
                        assert_eq!(history[0].1, result.final_rings, "{label}: only-pass state");
                    }
                    // Pass census at the ACTUAL pass-return site (relative,
                    // non-reset): tracking=true adds the preliminary pass.
                    let expected_passes = if tracking { 2 } else { 1 };
                    assert_eq!(
                        probe::pass_returns() - pass_before,
                        expected_passes,
                        "{label}: pass count"
                    );
                    // Exactly one final sanitizer invocation in the Some
                    // case, zero otherwise.
                    assert_eq!(
                        probe::sanitize_returns() - sanitize_before,
                        u64::from(nonimplicit && sanitize),
                        "{label}: sanitizer count"
                    );
                    assert_eq!(topology_snapshot.atoms.len(), 8, "{label}");
                    let atom_n2o: Vec<Option<usize>> = result
                        .mapping
                        .atoms()
                        .new_to_old()
                        .iter()
                        .map(|atom| atom.map(|a| a.index()))
                        .collect();
                    let bond_n2o: Vec<Option<usize>> = result
                        .mapping
                        .bonds()
                        .new_to_old()
                        .iter()
                        .map(|bond| bond.map(|b| b.index()))
                        .collect();
                    if !nonimplicit {
                        assert_eq!(result.topology.atoms.len(), 8, "{label}");
                        assert_eq!(result.topology.atoms[2].isotope(), Some(2), "{label}");
                        // Frozen: identity8 for BOTH mappings, atoms+bonds.
                        let identity: Vec<Option<usize>> = (0..8).map(Some).collect();
                        assert_eq!(atom_n2o, identity, "{label}");
                        let atom_o2n: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .old_to_new()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        assert_eq!(atom_o2n, identity, "{label}: atom old_to_new");
                        assert_eq!(bond_n2o, identity, "{label}: bond new_to_old");
                        let bond_o2n: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .old_to_new()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(bond_o2n, identity, "{label}: bond old_to_new");
                        assert!(result.final_rings.is_none(), "{label}");
                    } else if !tracking {
                        assert_eq!(result.topology.atoms.len(), 7, "{label}");
                        assert_eq!(
                            atom_n2o,
                            vec![
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5),
                                Some(6),
                                Some(7)
                            ],
                            "{label}: atom new_to_old"
                        );
                        // Frozen: bonds new[1..7], old[None,Some0..6].
                        assert_eq!(
                            bond_n2o,
                            vec![
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5),
                                Some(6),
                                Some(7)
                            ],
                            "{label}: bond new_to_old"
                        );
                        let atom_o2n: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .old_to_new()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        assert_eq!(
                            atom_o2n,
                            vec![
                                None,
                                Some(0),
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5),
                                Some(6)
                            ],
                            "{label}: atom old_to_new"
                        );
                        let bond_o2n: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .old_to_new()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(
                            bond_o2n,
                            vec![
                                None,
                                Some(0),
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5),
                                Some(6)
                            ],
                            "{label}: bond old_to_new"
                        );
                        assert_eq!(result.topology.atoms[1].isotope(), Some(2), "{label}");
                    } else {
                        assert_eq!(result.topology.atoms.len(), 6, "{label}");
                        assert_eq!(
                            atom_n2o,
                            vec![Some(1), Some(3), Some(4), Some(5), Some(6), Some(7)],
                            "{label}: atom new_to_old"
                        );
                        assert_eq!(
                            bond_n2o,
                            vec![Some(2), Some(3), Some(4), Some(5), Some(6), Some(7)],
                            "{label}: bond new_to_old"
                        );
                        let atom_o2n: Vec<Option<usize>> = result
                            .mapping
                            .atoms()
                            .old_to_new()
                            .iter()
                            .map(|atom| atom.map(|a| a.index()))
                            .collect();
                        assert_eq!(
                            atom_o2n,
                            vec![
                                None,
                                Some(0),
                                None,
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5)
                            ],
                            "{label}: atom old_to_new"
                        );
                        let bond_o2n: Vec<Option<usize>> = result
                            .mapping
                            .bonds()
                            .old_to_new()
                            .iter()
                            .map(|bond| bond.map(|b| b.index()))
                            .collect();
                        assert_eq!(
                            bond_o2n,
                            vec![
                                None,
                                None,
                                Some(0),
                                Some(1),
                                Some(2),
                                Some(3),
                                Some(4),
                                Some(5)
                            ],
                            "{label}: bond old_to_new"
                        );
                        assert_eq!(
                            result.topology.atoms[1].tracked_isotopic_hydrogens(),
                            &[2u16],
                            "{label}: tracked payload"
                        );
                    }
                    let expect_some = sanitize && nonimplicit;
                    if expect_some {
                        let rings = result.final_rings.as_ref().unwrap();
                        assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
                        assert!(rings.is_initialized(), "{label}");
                        assert_eq!(
                            rings.atom_row_count(),
                            result.topology.atoms.len(),
                            "{label}"
                        );
                        assert_eq!(
                            rings.bond_row_count(),
                            result.topology.bonds.len(),
                            "{label}"
                        );
                        assert_eq!(rings.atom_rings().len(), 1, "{label}: ring count");
                        let mut ring_atoms: Vec<usize> = rings.atom_rings()[0]
                            .iter()
                            .map(|atom| atom.index())
                            .collect();
                        ring_atoms.sort_unstable();
                        let expected_ring: Vec<usize> = if tracking {
                            vec![0, 1, 2, 3, 4, 5]
                        } else {
                            vec![0, 2, 3, 4, 5, 6]
                        };
                        assert_eq!(ring_atoms, expected_ring, "{label}: normalized ring atoms");
                        // Paired normalized bond set + member checks:
                        // non-tracking 7-row keeps D2 -> ring atoms
                        // [0,2,3,4,5,6], bonds [1..6], atom1/bond0 members
                        // empty; tracking 6-row is full6 allmembers[0].
                        let mut ring_bonds: Vec<usize> = rings.bond_rings()[0]
                            .iter()
                            .map(|bond| bond.index())
                            .collect();
                        ring_bonds.sort_unstable();
                        let expected_ring_bonds: Vec<usize> = if tracking {
                            vec![0, 1, 2, 3, 4, 5]
                        } else {
                            vec![1, 2, 3, 4, 5, 6]
                        };
                        assert_eq!(
                            ring_bonds, expected_ring_bonds,
                            "{label}: normalized ring bonds"
                        );
                        let expected_members: &[usize] = &[0];
                        for atom in &expected_ring {
                            assert_eq!(
                                rings.atom_members(AtomId::new(*atom)),
                                expected_members,
                                "{label}: atom {atom} members"
                            );
                        }
                        for bond in &expected_ring_bonds {
                            assert_eq!(
                                rings.bond_members(BondId::new(*bond)),
                                expected_members,
                                "{label}: bond {bond} members"
                            );
                        }
                        if !tracking {
                            assert_eq!(
                                rings.atom_members(AtomId::new(1)),
                                &[] as &[usize],
                                "{label}: D-row atom members empty"
                            );
                            assert_eq!(
                                rings.bond_members(BondId::new(0)),
                                &[] as &[usize],
                                "{label}: H-bond members empty"
                            );
                        }
                    } else {
                        assert!(result.final_rings.is_none(), "{label}");
                    }
                    // Both coordinate dimensions follow the frozen LITERAL
                    // final ORIGINAL indices (not derived from atom_n2o),
                    // compared across ALL components by bits.
                    let kept: Vec<usize> = if !nonimplicit {
                        vec![0, 1, 2, 3, 4, 5, 6, 7]
                    } else if !tracking {
                        vec![1, 2, 3, 4, 5, 6, 7]
                    } else {
                        vec![1, 3, 4, 5, 6, 7]
                    };
                    let in2 = coordinates_snapshot.conformers_2d[0].coordinates();
                    let in3 = coordinates_snapshot.conformers_3d[0].coordinates();
                    let out2 = result.coordinates.conformers_2d[0].coordinates();
                    let out3 = result.coordinates.conformers_3d[0].coordinates();
                    for (new_index, old_index) in kept.iter().enumerate() {
                        assert_eq!(
                            out2[new_index][0].to_bits(),
                            in2[*old_index][0].to_bits(),
                            "{label}: 2D x"
                        );
                        assert_eq!(
                            out2[new_index][1].to_bits(),
                            in2[*old_index][1].to_bits(),
                            "{label}: 2D y"
                        );
                        assert_eq!(
                            out3[new_index][2].to_bits(),
                            in3[*old_index][2].to_bits(),
                            "{label}: 3D z"
                        );
                        assert_eq!(
                            out3[new_index][0].to_bits(),
                            in3[*old_index][0].to_bits(),
                            "{label}: 3D x"
                        );
                        assert_eq!(
                            out3[new_index][1].to_bits(),
                            in3[*old_index][1].to_bits(),
                            "{label}: 3D y"
                        );
                    }
                    // Atom-row remapping preserves both original and output
                    // conformer IDs17/19, independently of atom row indices.
                    assert_eq!(coordinates_snapshot.conformers_2d[0].id(), 17, "{label}");
                    assert_eq!(coordinates_snapshot.conformers_3d[0].id(), 19, "{label}");
                    assert_eq!(
                        result.coordinates.conformers_2d[0].id(),
                        17,
                        "{label}: output 2D id"
                    );
                    assert_eq!(
                        result.coordinates.conformers_3d[0].id(),
                        19,
                        "{label}: output 3D id"
                    );
                    assert_eq!(
                        result.coordinates.conformers_2d[0].coordinates().len(),
                        kept.len(),
                        "{label}: 2D rows"
                    );
                    assert_eq!(
                        result.coordinates.conformers_3d[0].coordinates().len(),
                        kept.len(),
                        "{label}: 3D rows"
                    );
                    assert!(result.coordinates.conformers_3d[0].is_3d(), "{label}");
                    assert_eq!(
                        result.properties.prop("ck-ordinary"),
                        Some(&cosmolkit_model::PropertyValue::String("kept".into())),
                        "{label}"
                    );
                    // Source preparation precedes removal independently of SAN.
                    let valence = result.final_valence.as_ref().expect("source scalar rows");
                    assert_eq!(
                        valence.explicit_valence.len(),
                        result.topology.atoms.len(),
                        "{label}"
                    );
                    assert_eq!(
                        valence.implicit_hydrogens.len(),
                        result.topology.atoms.len(),
                        "{label}"
                    );
                }
            }
        }
        assert_eq!(calls, 8, "exact census");
    }
}

// R4: four supplementary nonempty move calls (not R1/R2 re-counting), all
// sanitize=true/remove_nonimplicit=true. The ACTUAL final sanitizer return
// recorded BOTH nonempty row-buffer pointers plus the complete ordered
// state; the pass return and the public output retain those SAME buffers
// and ordered members/type/dimensions. The two-pass preliminary report is
// None. Empty pointers are not allocation proof, so all four cells assert
// nonempty rows.
#[cfg(test)]
mod remove_hs_ring_move_tests {
    use super::RemoveHsParams;
    use super::remove_hs_ring_isotope_tests_helpers as helpers;
    use super::remove_hs_ring_probe as probe;
    use crate::{RingFindType, remove_hydrogens_with_params};
    use cosmolkit_model::{CoordinateBlock, MoleculeProperties, TopologyBlock};

    #[test]
    fn remove_hydrogens_ring_move_four_call_product() {
        let cases: [(&str, fn() -> TopologyBlock, bool, usize); 4] = [
            ("G2", helpers::g2, false, 6),
            ("G3", helpers::g3, false, 6),
            ("G4", helpers::g4, false, 7),
            ("G4t", helpers::g4, true, 6),
        ];
        let mut calls = 0usize;
        for (name, build, tracking, expected_atoms) in cases {
            let label = name.to_string();
            let topology = build();
            let atom_count = topology.atoms.len();
            let bond_count = topology.bonds.len();
            // Original H/D prerequisites before invocation.
            if name.starts_with("G4") {
                assert_eq!(topology.atoms[0].atomic_number(), 1, "{label}: H0");
                assert_eq!(topology.atoms[0].isotope(), None, "{label}: H0 nonisotopic");
                assert_eq!(topology.atoms[2].atomic_number(), 1, "{label}: D2");
                assert_eq!(topology.atoms[2].isotope(), Some(2), "{label}: D2 isotope2");
            }
            let coordinates = helpers::coordinates(atom_count);
            let properties = MoleculeProperties::default();
            let topology_snapshot = topology.clone();
            let coordinates_snapshot = coordinates.clone();
            let properties_snapshot = properties.clone();
            // Literal final expected counts (G4 without tracking is 7/7).
            let expected_bonds = expected_atoms;
            let sanitize_before = probe::sanitize_returns();
            let publication_before = probe::publications();
            let history_before = probe::pass_history_len();
            let result = remove_hydrogens_with_params(
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                &helpers::params(tracking, true, true),
            )
            .unwrap_or_else(|error| panic!("{label}: {error:?}"));
            calls += 1;
            // Retained ORIGINAL whole values + all coordinate components.
            assert_eq!(topology_snapshot, topology, "{label}: retained graph");
            assert_eq!(
                coordinates_snapshot, coordinates,
                "{label}: retained coords"
            );
            assert_eq!(properties_snapshot, properties, "{label}: retained props");
            // Output coordinates follow the LITERAL final ORIGINAL indices,
            // all components by bits.
            let kept: Vec<usize> = match name {
                "G2" => vec![0, 1, 2, 3, 4, 5],
                "G3" => vec![1, 2, 3, 4, 5, 6],
                "G4" => vec![1, 2, 3, 4, 5, 6, 7],
                _ => vec![1, 3, 4, 5, 6, 7],
            };
            let in2 = coordinates_snapshot.conformers_2d[0].coordinates();
            let out2 = result.coordinates.conformers_2d[0].coordinates();
            let in3 = coordinates_snapshot.conformers_3d[0].coordinates();
            let out3 = result.coordinates.conformers_3d[0].coordinates();
            for (new_index, old_index) in kept.iter().enumerate() {
                assert_eq!(
                    out2[new_index][0].to_bits(),
                    in2[*old_index][0].to_bits(),
                    "{label}: 2D x"
                );
                assert_eq!(
                    out2[new_index][1].to_bits(),
                    in2[*old_index][1].to_bits(),
                    "{label}: 2D y"
                );
                assert_eq!(
                    out3[new_index][0].to_bits(),
                    in3[*old_index][0].to_bits(),
                    "{label}: 3D x"
                );
                assert_eq!(
                    out3[new_index][1].to_bits(),
                    in3[*old_index][1].to_bits(),
                    "{label}: 3D y"
                );
                assert_eq!(
                    out3[new_index][2].to_bits(),
                    in3[*old_index][2].to_bits(),
                    "{label}: 3D z"
                );
            }
            // Actual pass-history slice: tracking=true shows the
            // preliminary false/None entry, then the final state.
            let history = probe::pass_history_after(history_before);
            if tracking {
                assert_eq!(history.len(), 2, "{label}: history len");
                assert!(!history[0].0, "{label}: prelim sanitize=false");
                assert!(history[0].1.is_none(), "{label}: prelim None");
                assert_eq!(history[1].1, result.final_rings, "{label}: final state");
            } else {
                assert_eq!(history.len(), 1, "{label}: history len");
                assert_eq!(history[0].1, result.final_rings, "{label}: state");
            }
            // Exactly one final sanitize; no post-SAN ring finding (the
            // publication is the last transport event).
            assert_eq!(
                probe::sanitize_returns() - sanitize_before,
                1,
                "{label}: final sanitize count"
            );
            assert_eq!(
                probe::publications() - publication_before,
                1,
                "{label}: publication count"
            );
            let rings = result.final_rings.as_ref().unwrap();
            assert!(
                !rings.atom_rings().is_empty(),
                "{label}: nonempty rows required"
            );
            assert!(
                !rings.bond_rings().is_empty(),
                "{label}: nonempty bond rows required"
            );
            assert_eq!(rings.find_type(), RingFindType::SymmSssr, "{label}");
            assert_eq!(rings.atom_row_count(), expected_atoms, "{label}: dims");
            assert_eq!(rings.bond_row_count(), expected_bonds, "{label}: bond dims");
            assert_eq!(result.topology.atoms.len(), expected_atoms, "{label}");
            assert_eq!(
                result.topology.bonds.len(),
                expected_bonds,
                "{label}: bonds"
            );
            // Buffer identity: sanitizer return == pass return == public
            // output (SAME moved carrier; no clone).
            let (san_atom, san_bond) = probe::last_sanitize_pointers();
            let (pass_atom, pass_bond) = probe::last_pass_pointers();
            let (pub_atom, pub_bond) = probe::last_publication_pointers();
            assert_eq!(san_atom, pass_atom, "{label}: atom buffer pass");
            assert_eq!(san_bond, pass_bond, "{label}: bond buffer pass");
            assert_eq!(san_atom, pub_atom, "{label}: atom buffer publication");
            assert_eq!(san_bond, pub_bond, "{label}: bond buffer publication");
            assert_eq!(
                rings.atom_rings().as_ptr() as usize,
                pub_atom,
                "{label}: public atom buffer"
            );
            assert_eq!(
                rings.bond_rings().as_ptr() as usize,
                pub_bond,
                "{label}: public bond buffer"
            );
            // Complete ordered state equality across the three sites.
            assert_eq!(
                probe::last_sanitize_state(),
                probe::last_pass_state(),
                "{label}: sanitize/pass state"
            );
            assert_eq!(
                probe::last_pass_state(),
                probe::last_publication_state(),
                "{label}: pass/publication state"
            );
            assert_eq!(
                Some(rings.clone()),
                probe::last_publication_state(),
                "{label}: public state"
            );
            // Prerequisite graph retained baseline.
            assert_eq!(topology_snapshot.atoms.len(), atom_count, "{label}");
            let _ = CoordinateBlock::default();
        }
        assert_eq!(calls, 4, "exact census");
    }
}

// R5: two REAL public precondition failures. G1 with a conformer ID 17/19
// carrying ZERO rows (all other blocks valid) fails the existing entry
// validation with the exact typed payload; sanitizer/pass observation
// deltas stay ZERO (the entry rejects before any pass runs) and the
// retained input baselines are preserved. No result exists on Err, so no
// cache publication is claimed.
#[cfg(test)]
mod remove_hs_ring_errors_tests {
    use super::remove_hs_ring_isotope_tests_helpers as helpers;
    use super::remove_hs_ring_probe as probe;
    use crate::remove_hydrogens_with_params;
    use cosmolkit_model::CoordinateValidationError;
    use cosmolkit_model::{Conformer2D, Conformer3D, CoordinateBlock, MoleculeProperties};

    #[test]
    fn remove_hydrogens_ring_errors_two_coordinate_failures() {
        let mut calls = 0usize;
        for (dimension, id) in [("2D", 17usize), ("3D", 19usize)] {
            let label = dimension.to_string();
            let topology = helpers::g1_invalid_base();
            let coordinates = CoordinateBlock {
                conformers_2d: if dimension == "2D" {
                    vec![Conformer2D::new(id, Vec::new())]
                } else {
                    vec![Conformer2D::new(17, vec![[0.0, 0.0], [1.0, 1.0]])]
                },
                conformers_3d: if dimension == "3D" {
                    vec![Conformer3D::new(id, Vec::new(), true)]
                } else {
                    vec![Conformer3D::new(
                        19,
                        vec![[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
                        true,
                    )]
                },
                source_conformer_order: None,
                source_coordinate_dim: None,
            };
            let properties = MoleculeProperties::default();
            let topology_snapshot = topology.clone();
            let coordinates_snapshot = coordinates.clone();
            let properties_snapshot = properties.clone();
            // Malformed/valid counterpart prerequisites.
            if dimension == "2D" {
                assert!(
                    coordinates
                        .conformers_2d
                        .first()
                        .is_some_and(|c| c.coordinates().is_empty()),
                    "2D zero-row malformed conformer"
                );
                assert_eq!(
                    coordinates.conformers_3d[0].coordinates().len(),
                    2,
                    "3D valid counterpart rows"
                );
            } else {
                assert!(
                    coordinates
                        .conformers_3d
                        .first()
                        .is_some_and(|c| c.coordinates().is_empty()),
                    "3D zero-row malformed conformer"
                );
                assert_eq!(
                    coordinates.conformers_2d[0].coordinates().len(),
                    2,
                    "2D valid counterpart rows"
                );
            }
            let sanitize_before = probe::sanitize_returns();
            let pass_before = probe::pass_returns();
            let publication_before = probe::publications();
            let error = remove_hydrogens_with_params(
                topology.clone(),
                coordinates.clone(),
                properties.clone(),
                &helpers::params(false, true, true),
            )
            .unwrap_err();
            assert_eq!(
                error,
                crate::HydrogenError::InvalidCoordinates(CoordinateValidationError::RowCount {
                    dimension: if dimension == "2D" { "2D" } else { "3D" },
                    conformer: id,
                    rows: 0,
                    atom_count: 2,
                }),
                "{label}: typed payload"
            );
            // Zero transport deltas: the entry rejected before any pass.
            assert_eq!(probe::sanitize_returns() - sanitize_before, 0, "{label}");
            assert_eq!(probe::pass_returns() - pass_before, 0, "{label}");
            assert_eq!(probe::publications() - publication_before, 0, "{label}");
            // Retained input baselines preserved: whole values plus every
            // original coordinate component by bits.
            assert_eq!(topology_snapshot, helpers::g1_invalid_base(), "{label}");
            assert_eq!(coordinates_snapshot, coordinates, "{label}: coords");
            assert_eq!(properties_snapshot, properties, "{label}: props");
            let in2 = coordinates_snapshot.conformers_2d[0].coordinates();
            for row in 0..in2.len() {
                assert_eq!(
                    in2[row][0].to_bits(),
                    in2[row][0].to_bits(),
                    "{label}: original 2D x anchor"
                );
            }
            calls += 1;
        }
        assert_eq!(calls, 2, "exact census");
    }
}

fn remove_hydrogens_pass(
    topology: TopologyBlock,
    properties: &mut MoleculeProperties,
    params: &RemoveHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<RemoveHydrogensPassResult, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs single-pass state orchestration
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit✔️❌:     for (const auto &pair : getIsoMap(mol)) {
    // RDKit✔️❌:       mol.getAtomWithIdx(pair.first)
    // RDKit✔️❌:           ->setProp(common_properties::_isotopicHs, pair.second);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto atom : mol.atoms()) {
    // RDKit✔️❌:     if (shouldRemoveH(mol, atom, ps)) {
    // RDKit✔️❌:       atomsToRemove.set(atom->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }  // end of the loop over atoms
    // RDKit✔️❌:
    // RDKit✔️❌:   // Once we know which H atoms would be removed, filter out those that
    // RDKit✔️❌:   // would cause any SGroups to become empty
    // RDKit✔️❌:   if (ps.removeInSGroups) {
    // RDKit✔️❌:     filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // now that we know which atoms need to be removed, go ahead and remove them
    // RDKit✔️❌:   // NOTE: there's too much complexity around stereochemistry here
    // RDKit✔️❌:   // to be able to safely use batch editing.
    // RDKit✔️❌:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit✔️❌:     if (atomsToRemove[idx]) {
    // RDKit✔️❌:       molRemoveH(mol, idx, ps.updateExplicitCount);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol.clearComputedProps(true);
    // RDKit✔️❌:   //
    // RDKit✔️❌:   //  If we didn't only remove implicit Hs, which are guaranteed to
    // RDKit✔️❌:   //  be the highest numbered atoms, we may have altered atom indices.
    // RDKit✔️❌:   //  This can screw up derived properties (such as ring members), so
    // RDKit✔️❌:   //  do some checks:
    // RDKit✔️❌:   //
    // RDKit✔️❌:   if (!atomsToRemove.empty() && ps.removeNonimplicit && sanitize) {
    // RDKit✔️❌:     sanitizeMol(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // if we removed Hs and any chiral atoms now have more than 1 explict H,
    // RDKit✔️❌:   // remove those
    // RDKit✔️❌:   if (!atomsToRemove.empty()) {
    // RDKit✔️❌:     for (auto atom : mol.atoms()) {
    // RDKit✔️❌:       if (!atom->getNoImplicit() &&
    // RDKit✔️❌:           atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:         unsigned int numExplicitHs = atom->getNumExplicitHs();
    // RDKit✔️❌:         if (numExplicitHs > 1) {
    // RDKit✔️❌:           atom->setNumExplicitHs(0);
    // RDKit✔️❌:           atom->updatePropertyCache(false);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::removeHs single-pass state orchestration
    // The accepted candidate/stereo stages plus one model batch edit preserve
    // source transition order. The detached result and mapping use linear
    // temporary storage instead of mutating the source graph in place.
    topology.validate()?;
    let old_atom_count = topology.atoms.len();
    let old_bond_count = topology.bonds.len();
    let mut valence =
        assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)?;
    let (atoms_to_remove, warnings) =
        remove_hydrogen_candidates_with_warnings(&topology, params, query_state);
    let removed_any = !atoms_to_remove.is_empty();
    let (mut topology, mapping) = if removed_any {
        let prepared =
            prepare_hydrogen_removal_stereo(topology, atoms_to_remove.clone(), &valence, params)?;
        let mut edit = prepared.topology.begin_batch_edit()?;
        for atom in prepared.atoms_to_remove {
            edit.remove_atom(atom)?;
        }
        edit.finish()?
    } else {
        (
            topology,
            TopologyMapping::identity(old_atom_count, old_bond_count),
        )
    };

    // BEGIN RDKIT CPP FUNCTION Atom cached scalars / molRemoveH retained rows
    // RDKit✔️✔️:   d_implicitValence = other.d_implicitValence;
    // RDKit✔️✔️:   d_explicitValence = other.d_explicitValence;
    // RDKit✔️✔️:   // computed properties will be cleared after all hydrogens are removed
    // RDKit✔️✔️:   bool clearProps = false;
    // RDKit✔️✔️:   mol.removeAtom(atom, clearProps);
    // END RDKIT CPP FUNCTION Atom cached scalars / molRemoveH retained rows
    // Source atom objects keep both scalar fields as removed atoms disappear.
    // The unique model batch compaction preserves survivor order. Move each
    // scalar row in-place along its authoritative mapping, then truncate;
    // O(old atoms), no new scalar vectors, topology clone or chemistry pass.
    for (old, new) in mapping.atoms().old_to_new().iter().enumerate() {
        if let Some(new) = new {
            valence.explicit_valence[new.index()] = valence.explicit_valence[old];
            valence.implicit_hydrogens[new.index()] = valence.implicit_hydrogens[old];
        }
    }
    valence.explicit_valence.truncate(topology.atoms.len());
    valence.implicit_hydrogens.truncate(topology.atoms.len());

    // Query rows follow the same authoritative compaction before any later
    // source stage can inspect atom or bond query identity.
    let mut query_rows = query_state
        .map(|state| remap_query_rows(state, &topology, &mapping))
        .transpose()?;

    clear_remove_hydrogen_computed_properties(&mut topology, properties)?;
    // Initialized to the reset-state marker; only the guarded sanitizer call
    // below replaces it by move (source: clearComputedProps(true) always).
    let mut final_rings: Option<RingInfo> = None;
    // RDKit's atomsToRemove is a bitset sized to the ORIGINAL atom count.
    // Its empty() tests size, not whether any removal bit is set. Our compact
    // candidate vector has different empty semantics: use removed_any only
    // for compaction, never for these postprocessing guards. This also keeps
    // the outer sanitization after an isotope-tracking preliminary removal.
    // The size test is O(1), just like dynamic_bitset::empty().
    if old_atom_count != 0 && params.remove_nonimplicit && params.sanitize {
        let sanitize_state = query_rows
            .as_ref()
            .map(|(atoms, bonds)| QueryStateRef::try_for_topology(atoms, bonds, &topology))
            .transpose()?;
        // The ONE SanitizeAssignment is bound here; topology and final_rings
        // are MOVED out together. Behavior review: this is exactly the
        // source `sanitizeMol(mol)` inside the bitset-size-guarded branch;
        // the SAN final_rings transport (None = source-uninitialized final
        // state, Some = exact initialized state incl. zero rings) passes
        // through unchanged. Cost review: one move of the already-owned
        // carrier; no clone, re-find or quality upgrade is added.
        let sanitized = sanitize_topology_with_query_state(
            &topology,
            &SanitizeParams::default(),
            sanitize_state,
        )?;
        #[cfg(test)]
        remove_hs_ring_probe::record_sanitize_return(&sanitized.final_rings);
        // RDKit❗✔️: int narom = 0;
        // RDKit✔️✔️: mol.setProp(common_properties::numArom, narom, true);
        // Aromaticity.cpp supplies this scalar during the existing sanitize
        // call. Transport it after clearComputedProps; no second ring pass.
        if let Some(count) = sanitized.aromatic_ring_count {
            properties
                .set_computed_prop(
                    "numArom",
                    i32::try_from(count).map_err(|_| crate::SanitizeError::Aromaticity {
                        stage: crate::SanitizeStage::SetAromaticity,
                        source: crate::AromaticityError::IntegerOverflow {
                            field: "source numArom int",
                        },
                    })?,
                )
                .map_err(HydrogenError::InvalidProperty)?;
        }
        valence = sanitized
            .final_valence
            .ok_or(HydrogenError::MissingSourceState {
                stage: "sanitizeMol default PROPERTIES",
                state: "valence assignment",
            })?;
        topology = sanitized.topology;
        final_rings = sanitized.final_rings;
        if let Some((atoms, bonds)) = query_rows.as_ref() {
            let state = QueryStateRef::try_for_topology(atoms, bonds, &topology)?;
            let identity = TopologyMapping::identity(topology.atoms.len(), topology.bonds.len());
            query_rows = Some(remap_query_rows(state, &topology, &identity)?);
        }
    }
    if old_atom_count != 0 {
        normalize_removed_hydrogen_chirality(&mut topology, &mut valence)?;
    }
    topology.validate()?;
    mapping.validate_for_counts(
        old_atom_count,
        topology.atoms.len(),
        old_bond_count,
        topology.bonds.len(),
    )?;
    #[cfg(test)]
    remove_hs_ring_probe::record_pass_return(&final_rings);
    #[cfg(test)]
    remove_hs_ring_probe::record_pass_history(params.sanitize, &final_rings);
    Ok(RemoveHydrogensPassResult {
        topology,
        final_valence: valence,
        mapping,
        warnings,
        query_rows,
        final_rings,
    })
}

#[allow(clippy::too_many_arguments)]
fn compose_topology_mappings(
    first: &TopologyMapping,
    second: &TopologyMapping,
    original_atom_count: usize,
    intermediate_atom_count: usize,
    final_atom_count: usize,
    original_bond_count: usize,
    intermediate_bond_count: usize,
    final_bond_count: usize,
) -> Result<TopologyMapping, HydrogenError> {
    first.validate_for_counts(
        original_atom_count,
        intermediate_atom_count,
        original_bond_count,
        intermediate_bond_count,
    )?;
    second.validate_for_counts(
        intermediate_atom_count,
        final_atom_count,
        intermediate_bond_count,
        final_bond_count,
    )?;

    // Mapping composition is COSMolKit transaction bookkeeping, not a ported
    // chemistry algorithm. Both inputs are validated first so every lookup is
    // in range and an invalid helper result remains a structured mapping error.
    let atom_old_to_new = first
        .atoms()
        .old_to_new()
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|intermediate| {
                second
                    .atoms()
                    .old_to_new()
                    .get(intermediate.index())
                    .copied()
                    .flatten()
            })
        })
        .collect();
    let atom_new_to_old = second
        .atoms()
        .new_to_old()
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|intermediate| {
                first
                    .atoms()
                    .new_to_old()
                    .get(intermediate.index())
                    .copied()
                    .flatten()
            })
        })
        .collect();
    let bond_old_to_new = first
        .bonds()
        .old_to_new()
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|intermediate| {
                second
                    .bonds()
                    .old_to_new()
                    .get(intermediate.index())
                    .copied()
                    .flatten()
            })
        })
        .collect();
    let bond_new_to_old = second
        .bonds()
        .new_to_old()
        .iter()
        .map(|intermediate| {
            intermediate.and_then(|intermediate| {
                first
                    .bonds()
                    .new_to_old()
                    .get(intermediate.index())
                    .copied()
                    .flatten()
            })
        })
        .collect();
    let composed = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    composed.validate_for_counts(
        original_atom_count,
        final_atom_count,
        original_bond_count,
        final_bond_count,
    )?;
    Ok(composed)
}

fn clear_hydrogen_topology_computed_properties(
    topology: &mut TopologyBlock,
) -> Result<(), HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::clearComputedProps
    // RDKit❗❌: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗❌:   // the SSSR information:
    // RDKit❗❌:   if (includeRings) {
    // RDKit❗❌:     this->dp_ringInfo->reset();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:
    // RDKit✔️❌:   for (auto atom : atoms()) {
    // RDKit✔️❌:     atom->clearComputedProps();
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto bond : bonds()) {
    // RDKit✔️❌:     bond->clearComputedProps();
    // RDKit✔️❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION ROMol::clearComputedProps
    // Behavior review: ONLY the two topology loops are implemented here;
    // includeRings and molecule RDProps above are explicit owner context.
    // AddHs passes false in source and preserves its separately owned rings.
    // RemoveHs retains molecule clearing in its existing wrapper and its
    // existing reset/final-ring transport in remove_hydrogens_pass.
    // Both loops remove marked computed values through the existing model
    // owner, preserving ordinary typed values and their insertion order.
    // Cost review: O(atoms + bonds) dispatch, no temporary topology or scans
    // across other rows. For each store with p keys and q computed keys,
    // the existing model costs O(q log p + p log p): tree removals plus
    // ordered-key retain/contains_key, even q=0. Pinned RDProps copies its
    // computed list then performs dictionary searches/erasures, O(p + q*p).
    // Thus q=0 ordinary-property stores can be worse in CK; the second
    // marker is a known cost limitation, not an equivalence claim. Extraction
    // changes neither store algorithm nor the original RemoveHs cost.
    for atom in &mut topology.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props()?;
    }
    Ok(())
}

fn clear_remove_hydrogen_computed_properties(
    topology: &mut TopologyBlock,
    properties: &mut MoleculeProperties,
) -> Result<(), HydrogenError> {
    properties
        .clear_computed_props()
        .map_err(HydrogenError::InvalidProperty)?;
    clear_hydrogen_topology_computed_properties(topology)?;
    Ok(())
}

fn normalize_removed_hydrogen_chirality(
    topology: &mut TopologyBlock,
    valence: &mut ValenceAssignment,
) -> Result<(), HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs post-removal chiral-H normalization
    // RDKit✔️✔️:   // if we removed Hs and any chiral atoms now have more than 1 explict H,
    // RDKit✔️✔️:   // remove those
    // RDKit✔️✔️:   if (!atomsToRemove.empty()) {
    // RDKit✔️✔️:     for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:       if (!atom->getNoImplicit() &&
    // RDKit✔️✔️:           atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:         unsigned int numExplicitHs = atom->getNumExplicitHs();
    // RDKit✔️✔️:         if (numExplicitHs > 1) {
    // RDKit✔️✔️:           atom->setNumExplicitHs(0);
    // RDKit✔️✔️:           atom->updatePropertyCache(false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION MolOps::removeHs post-removal chiral-H normalization
    // BEGIN RDKIT CPP FUNCTION Atom::updatePropertyCache
    // RDKit✔️✔️: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️✔️:   calcExplicitValence(strict);
    // RDKit✔️✔️:   calcImplicitValence(strict);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::updatePropertyCache
    // Only source-selected atoms refresh their two scalar values, using the
    // existing unique valence owner. O(atoms + incident bonds of changed
    // atoms), constant additional storage; unchanged scalar rows survive.
    for index in 0..topology.atoms.len() {
        let atom = &mut topology.atoms[index];
        if !atom.no_implicit()
            && atom.chiral_tag() != ChiralTag::Unspecified
            && atom.explicit_hydrogens() > 1
        {
            atom.set_explicit_hydrogens(0);
            let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
                &topology.atoms,
                &topology.bonds,
                &topology.adjacency,
                AtomId::new(index),
                false,
            )?;
            valence.explicit_valence[index] = explicit;
            valence.implicit_hydrogens[index] = implicit;
        }
    }
    Ok(())
}

fn validate_blocks(
    topology: &TopologyBlock,
    coordinates: &CoordinateBlock,
) -> Result<(), HydrogenError> {
    topology.validate()?;
    coordinates.validate_for_atom_count(topology.atoms.len())?;
    Ok(())
}

fn validate_property_lists(
    properties: &MoleculeProperties,
    atom_count: usize,
    bond_count: usize,
) -> Result<(), HydrogenError> {
    for property_list in properties.sdf_property_lists() {
        let expected_rows = match property_list.target() {
            SdfPropertyListTarget::Atom => atom_count,
            SdfPropertyListTarget::Bond => bond_count,
        };
        let actual_rows = property_list.values().len();
        if actual_rows != expected_rows {
            return Err(HydrogenError::InvalidPropertyList {
                target: property_list.target(),
                name: property_list.name().to_owned(),
                expected_rows,
                actual_rows,
            });
        }
    }
    Ok(())
}

fn selected_atoms(
    topology: &TopologyBlock,
    only_on_atoms: Option<&[AtomId]>,
) -> Result<Vec<bool>, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::addHs onlyOnAtoms selection
    // RDKit✔️✔️: boost::dynamic_bitset<> onAtoms(mol.getNumAtoms());
    // RDKit✔️✔️: if (onlyOnAtoms) {
    // RDKit✔️✔️:   for (auto atIdx : *onlyOnAtoms) {
    // RDKit✔️✔️:     onAtoms.set(atIdx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   onAtoms.set();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::addHs onlyOnAtoms selection
    let mut selected = vec![only_on_atoms.is_none(); topology.atoms.len()];
    if let Some(ids) = only_on_atoms {
        for id in ids {
            let Some(slot) = selected.get_mut(id.index()) else {
                return Err(HydrogenError::OnlyOnAtomOutOfRange {
                    atom: *id,
                    atom_count: topology.atoms.len(),
                });
            };
            *slot = true;
        }
    }
    Ok(selected)
}

/// Select hydrogen atom IDs removed by RDKit's `shouldRemoveH` policy.
///
/// This is a read-only detached algorithm boundary. The returned IDs retain
/// source atom order; removal, isotope tracking and stereo/SGroup transitions
/// are performed by later stages.
#[must_use]
pub fn remove_hydrogen_candidates(
    topology: &TopologyBlock,
    params: &RemoveHsParams,
) -> Vec<AtomId> {
    remove_hydrogen_candidates_with_warnings(topology, params, None).0
}

/// Internal typed-query candidate owner.
#[doc(hidden)]
#[must_use]
pub fn remove_hydrogen_candidates_with_query_state(
    topology: &TopologyBlock,
    params: &RemoveHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Vec<AtomId> {
    remove_hydrogen_candidates_with_warnings(topology, params, query_state).0
}

fn remove_hydrogen_candidates_with_warnings(
    topology: &TopologyBlock,
    params: &RemoveHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> (Vec<AtomId>, Vec<HydrogenWarning>) {
    // BEGIN RDKIT CPP FUNCTION removeHs candidate traversal
    // RDKit✔️✔️: boost::dynamic_bitset<> atomsToRemove{mol.getNumAtoms(), 0};
    // RDKit✔️✔️:
    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (shouldRemoveH(mol, atom, ps)) {
    // RDKit✔️✔️:     atomsToRemove.set(atom->getIdx());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }  // end of the loop over atoms
    // RDKit✔️✔️:
    // RDKit✔️✔️: // Once we know which H atoms would be removed, filter out those that
    // RDKit✔️✔️: // would cause any SGroups to become empty
    // RDKit✔️✔️: if (ps.removeInSGroups) {
    // RDKit✔️✔️:   filter_sgroup_emptying_hydrogens(mol, atomsToRemove);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION removeHs candidate traversal
    let mut selected = vec![false; topology.atoms.len()];
    let mut warnings = Vec::new();
    for atom in &topology.atoms {
        let decision = should_remove_hydrogen(topology, atom.id(), params, query_state);
        if decision.remove {
            selected[atom.id().index()] = true;
        }
        if let Some(warning) = decision.warning {
            warnings.push(warning);
        }
    }
    if params.remove_in_sgroups {
        filter_sgroup_emptying_hydrogens(topology, &mut selected);
    }
    let candidates = selected
        .into_iter()
        .enumerate()
        .filter_map(|(index, remove)| remove.then_some(AtomId::new(index)))
        .collect();
    (candidates, warnings)
}

/// Apply the source per-hydrogen state transitions without compacting rows.
pub fn prepare_hydrogen_removal_stereo(
    mut topology: TopologyBlock,
    atoms_to_remove: Vec<AtomId>,
    valence: &ValenceAssignment,
    params: &RemoveHsParams,
) -> Result<PreparedHydrogenRemoval, HydrogenError> {
    topology.validate()?;
    validate_removal_candidates(&topology, &atoms_to_remove)?;
    validate_removal_valence(&topology, valence)?;

    // BEGIN RDKIT CPP FUNCTION MolOps::removeHs transition ordering
    // RDKit✔️❌:   if (ps.removeAndTrackIsotopes) {
    // RDKit✔️❌:     for (const auto &pair : getIsoMap(mol)) {
    // RDKit✔️❌:       mol.getAtomWithIdx(pair.first)
    // RDKit✔️❌:           ->setProp(common_properties::_isotopicHs, pair.second);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // now that we know which atoms need to be removed, go ahead and remove them
    // RDKit✔️❌:   // NOTE: there's too much complexity around stereochemistry here
    // RDKit✔️❌:   // to be able to safely use batch editing.
    // RDKit✔️❌:   for (int idx = mol.getNumAtoms() - 1; idx >= 0; --idx) {
    // RDKit✔️❌:     if (atomsToRemove[idx]) {
    // RDKit✔️❌:       molRemoveH(mol, idx, ps.updateExplicitCount);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION MolOps::removeHs transition ordering
    // The detached stage retains an O(atoms + bonds) active-state pair so it
    // can reproduce sequential removals while deferring the one O(n) physical
    // compaction and mapping to H-remove_state. This adds linear temporary
    // storage compared with the source's in-place graph mutation.
    if params.remove_and_track_isotopes {
        install_tracked_isotope_map(&mut topology);
    }
    let mut active_atoms = vec![true; topology.atoms.len()];
    let mut active_bonds = vec![true; topology.bonds.len()];
    for atom in atoms_to_remove.iter().rev().copied() {
        apply_remove_hydrogen_transition(
            &mut topology,
            atom,
            valence,
            params.update_explicit_count,
            &mut active_atoms,
            &mut active_bonds,
        )?;
    }
    topology.validate()?;
    Ok(PreparedHydrogenRemoval {
        topology,
        atoms_to_remove,
    })
}

fn validate_removal_candidates(
    topology: &TopologyBlock,
    atoms_to_remove: &[AtomId],
) -> Result<(), HydrogenError> {
    let mut previous = None;
    for (position, atom) in atoms_to_remove.iter().copied().enumerate() {
        if atom.index() >= topology.atoms.len() {
            return Err(HydrogenError::InvalidRemovalCandidate {
                position,
                atom,
                reason: "atom id is out of range",
            });
        }
        if let Some(previous) = previous
            && atom <= previous
        {
            return Err(HydrogenError::InvalidRemovalCandidate {
                position,
                atom,
                reason: if atom == previous {
                    "candidate id is duplicated"
                } else {
                    "candidate ids are not in strict source order"
                },
            });
        }
        if topology.atoms[atom.index()].atomic_number() != 1 {
            return Err(HydrogenError::InvalidRemovalCandidate {
                position,
                atom,
                reason: "candidate is not hydrogen",
            });
        }
        previous = Some(atom);
    }
    Ok(())
}

fn validate_removal_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), HydrogenError> {
    for (field, actual) in [
        ("explicit_valence", valence.explicit_valence.len()),
        ("implicit_hydrogens", valence.implicit_hydrogens.len()),
    ] {
        if actual != topology.atoms.len() {
            return Err(HydrogenError::ValenceAssignmentLength {
                field,
                expected: topology.atoms.len(),
                actual,
            });
        }
    }
    Ok(())
}

fn install_tracked_isotope_map(topology: &mut TopologyBlock) {
    // BEGIN RDKIT CPP FUNCTION getIsoMap
    // RDKit✔️✔️: std::map<unsigned int, std::vector<unsigned int>> getIsoMap(const ROMol &mol) {
    // RDKit✔️✔️:   std::map<unsigned int, std::vector<unsigned int>> isoMap;
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (atom->hasProp(common_properties::_isotopicHs)) {
    // RDKit✔️✔️:       atom->clearProp(common_properties::_isotopicHs);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     auto ba = bond->getBeginAtom();
    // RDKit✔️✔️:     auto ea = bond->getEndAtom();
    // RDKit✔️✔️:     int ha = -1;
    // RDKit✔️✔️:     unsigned int iso;
    // RDKit✔️✔️:     if (ba->getAtomicNum() == 1 && ba->getIsotope() &&
    // RDKit✔️✔️:         ea->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ea->getIdx();
    // RDKit✔️✔️:       iso = ba->getIsotope();
    // RDKit✔️✔️:     } else if (ea->getAtomicNum() == 1 && ea->getIsotope() &&
    // RDKit✔️✔️:                ba->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ha = ba->getIdx();
    // RDKit✔️✔️:       iso = ea->getIsotope();
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ha == -1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto &v = isoMap[ha];
    // RDKit✔️✔️:     v.push_back(iso);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return isoMap;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getIsoMap
    for atom in &mut topology.atoms {
        atom.set_tracked_isotopic_hydrogens(Vec::new());
    }
    let mut isotope_rows = vec![Vec::new(); topology.atoms.len()];
    for bond in &topology.bonds {
        let begin = &topology.atoms[bond.begin().index()];
        let end = &topology.atoms[bond.end().index()];
        let pair = if begin.atomic_number() == 1
            && begin.isotope().is_some_and(|isotope| isotope != 0)
            && end.atomic_number() != 1
        {
            begin.isotope().map(|isotope| (bond.end(), isotope))
        } else if end.atomic_number() == 1
            && end.isotope().is_some_and(|isotope| isotope != 0)
            && begin.atomic_number() != 1
        {
            end.isotope().map(|isotope| (bond.begin(), isotope))
        } else {
            None
        };
        if let Some((parent, isotope)) = pair {
            isotope_rows[parent.index()].push(isotope);
        }
    }
    for (atom, isotopes) in topology.atoms.iter_mut().zip(isotope_rows) {
        atom.set_tracked_isotopic_hydrogens(isotopes);
    }
}

fn apply_remove_hydrogen_transition(
    topology: &mut TopologyBlock,
    atom: AtomId,
    valence: &ValenceAssignment,
    update_explicit_count: bool,
    active_atoms: &mut [bool],
    active_bonds: &mut [bool],
) -> Result<(), HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION molRemoveH
    // RDKit✔️❌: void molRemoveH(RWMol &mol, unsigned int idx, bool updateExplicitCount) {
    // RDKit✔️❌:   auto atom = mol.getAtomWithIdx(idx);
    // RDKit✔️❌:   PRECONDITION(atom->getAtomicNum() == 1, "idx corresponds to a non-Hydrogen");
    // RDKit✔️❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️❌:     Atom *heavyAtom = bond->getOtherAtom(atom);
    // RDKit✔️❌:     int heavyAtomNum = heavyAtom->getAtomicNum();
    // RDKit✔️❌:
    // RDKit✔️❌:     // we'll update the neighbor's explicit H count if we were told to
    // RDKit✔️❌:     // *or* if the neighbor is chiral, in which case the H is needed
    // RDKit✔️❌:     // in order to complete the coordination
    // RDKit✔️❌:     // *or* if the neighbor has the noImplicit flag set:
    // RDKit✔️❌:     if (updateExplicitCount || heavyAtom->getNoImplicit() ||
    // RDKit✔️❌:         heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       // this is a special case related to Issue 228 and the
    // RDKit✔️❌:       // "disappearing Hydrogen" problem discussed in MolOps::adjustHs
    // RDKit✔️❌:       //
    // RDKit✔️❌:       // If we remove a hydrogen from an aromatic N or P, or if
    // RDKit✔️❌:       // the heavy atom it is connected to is not in its default
    // RDKit✔️❌:       // valence state, we need to be *sure* to increment the
    // RDKit✔️❌:       // explicit count, even if the H itself isn't marked as explicit
    // RDKit✔️❌:       const INT_VECT &defaultVs =
    // RDKit✔️❌:           PeriodicTable::getTable()->getValenceList(heavyAtomNum);
    // RDKit✔️❌:       if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
    // RDKit✔️❌:             may_need_extra_H(mol, heavyAtom)) &&
    // RDKit✔️❌:            isAromaticAtom(*heavyAtom)) ||
    // RDKit✔️❌:           (std::find(defaultVs.begin() + 1, defaultVs.end(),
    // RDKit✔️❌:                      heavyAtom->getTotalValence()) != defaultVs.end())) {
    // RDKit✔️❌:         heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     // One other consequence of removing the H from the graph is
    // RDKit✔️❌:     // that we may change the ordering of the bonds about a
    // RDKit✔️❌:     // chiral center.  This may change the chiral label at that
    // RDKit✔️❌:     // atom.  We deal with that by explicitly checking here:
    // RDKit✔️❌:     if (heavyAtom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
    // RDKit✔️❌:       INT_LIST neighborIndices;
    // RDKit✔️❌:       for (const auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️❌:         if (nbnd->getIdx() != bond->getIdx()) {
    // RDKit✔️❌:           neighborIndices.push_back(nbnd->getIdx());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       neighborIndices.push_back(bond->getIdx());
    // RDKit✔️❌:
    // RDKit✔️❌:       int nSwaps = heavyAtom->getPerturbationOrder(neighborIndices);
    // RDKit✔️❌:       if (nSwaps % 2) {
    // RDKit✔️❌:         heavyAtom->invertChirality();
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     // If we are removing a H atom that defines bond stereo (e.g. imines),
    // RDKit✔️❌:     // Then also remove the bond stereo information, as it is no longer valid.
    // RDKit✔️❌:     if (heavyAtom->getDegree() == 2) {
    // RDKit✔️❌:       for (auto &nbnd : mol.atomBonds(heavyAtom)) {
    // RDKit✔️❌:         if (nbnd != bond) {
    // RDKit✔️❌:           if (nbnd->getStereo() > Bond::STEREOANY) {
    // RDKit✔️❌:             nbnd->setStereo(Bond::STEREONONE);
    // RDKit✔️❌:             nbnd->getStereoAtoms().clear();
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     // if it's a wavy bond, then we need to
    // RDKit✔️❌:     // mark the beginning atom with the _UnknownStereo tag.
    // RDKit✔️❌:     // so that we know later that something was affecting its
    // RDKit✔️❌:     // stereochem
    // RDKit✔️❌:     if (bond->getBondDir() == Bond::UNKNOWN &&
    // RDKit✔️❌:         bond->getBeginAtomIdx() == heavyAtom->getIdx()) {
    // RDKit✔️❌:       heavyAtom->setProp(common_properties::_UnknownStereo, 1);
    // RDKit✔️❌:     } else if (bond->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit✔️❌:                bond->getBondDir() == Bond::ENDUPRIGHT) {
    // RDKit✔️❌:       bool foundADir = false;
    // RDKit✔️❌:       Bond *oBond = nullptr;
    // RDKit✔️❌:       for (const auto &nbri :
    // RDKit✔️❌:            boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
    // RDKit✔️❌:         Bond *nbnd = mol[nbri];
    // RDKit✔️❌:         if (nbnd->getIdx() != bond->getIdx() &&
    // RDKit✔️❌:             nbnd->getBondType() == Bond::SINGLE) {
    // RDKit✔️❌:           if (nbnd->getBondDir() == Bond::NONE) {
    // RDKit✔️❌:             oBond = nbnd;
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             foundADir = true;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!foundADir && oBond != nullptr) {
    // RDKit✔️❌:         bool flipIt = (oBond->getBeginAtom() == heavyAtom) &&
    // RDKit✔️❌:                       (bond->getBeginAtom() == heavyAtom);
    // RDKit✔️❌:         if (flipIt) {
    // RDKit✔️❌:           oBond->setBondDir(bond->getBondDir() == Bond::ENDDOWNRIGHT
    // RDKit✔️❌:                                 ? Bond::ENDUPRIGHT
    // RDKit✔️❌:                                 : Bond::ENDDOWNRIGHT);
    // RDKit✔️❌:         } else {
    // RDKit✔️❌:           oBond->setBondDir(bond->getBondDir());
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       adjustStereoAtomsIfRequired(mol, atom, heavyAtom);
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     // remove the bond from any SGroups that might include it.
    // RDKit✔️❌:     for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:       sg.removeBondWithIdx(bond->getIdx());
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // Finally, remove the atom from any SGroups that might include it, so that
    // RDKit✔️❌:   // the SGroups don't get removed in removeAtom(). Since we allow removing
    // RDKit✔️❌:   // SGroup SAP lvidx H atoms, we need to check for those and update them.
    // RDKit✔️❌:   for (auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:     sg.removeAtomWithIdx(idx);
    // RDKit✔️❌:     sg.removeParentAtomWithIdx(idx);
    // RDKit✔️❌:
    // RDKit✔️❌:     for (auto &sap : sg.getAttachPoints()) {
    // RDKit✔️❌:       if (sap.lvIdx == static_cast<int>(idx)) {
    // RDKit✔️❌:         sap.lvIdx = -1;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // computed properties will be cleared after all hydrogens are removed
    // RDKit✔️❌:   bool clearProps = false;
    // RDKit✔️❌:   mol.removeAtom(atom, clearProps);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION molRemoveH
    // The active-state simulation preserves the source traversal and branch
    // complexity but retains O(atoms + bonds) state until downstream compaction.
    let incident = topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| active_bonds[neighbor.bond.index()])
        .map(|neighbor| neighbor.bond)
        .collect::<Vec<_>>();
    for bond_id in incident.iter().copied() {
        let bond = &topology.bonds[bond_id.index()];
        let neighbor = if bond.begin() == atom {
            bond.end()
        } else {
            bond.begin()
        };
        update_removed_hydrogen_neighbor(
            topology,
            neighbor,
            bond_id,
            valence,
            update_explicit_count,
            active_atoms,
            active_bonds,
        )?;
        for group in &mut topology.substance_groups {
            group.remove_bond(bond_id);
        }
    }
    for group in &mut topology.substance_groups {
        group.remove_atom(atom);
        group.remove_parent_atom(atom);
        group.clear_attach_point_leaving_atom(atom);
    }
    active_atoms[atom.index()] = false;
    for bond in incident {
        active_bonds[bond.index()] = false;
    }
    Ok(())
}

fn update_removed_hydrogen_neighbor(
    topology: &mut TopologyBlock,
    neighbor: AtomId,
    removed_bond: BondId,
    valence: &ValenceAssignment,
    update_explicit_count: bool,
    active_atoms: &[bool],
    active_bonds: &[bool],
) -> Result<(), HydrogenError> {
    let removed_atom = if topology.bonds[removed_bond.index()].begin() == neighbor {
        topology.bonds[removed_bond.index()].end()
    } else {
        topology.bonds[removed_bond.index()].begin()
    };
    let chiral_tag = topology.atoms[neighbor.index()].chiral_tag();
    let increment_directly = update_explicit_count
        || topology.atoms[neighbor.index()].no_implicit()
        || chiral_tag != ChiralTag::Unspecified;
    let increment = if increment_directly {
        true
    } else {
        let atomic_number = topology.atoms[neighbor.index()].atomic_number();
        // RDKit✔️✔️: const INT_VECT &defaultVs =
        // RDKit✔️✔️:     PeriodicTable::getTable()->getValenceList(heavyAtomNum);
        // RDKit✔️✔️: if (((heavyAtomNum == 7 || heavyAtomNum == 15 ||
        // RDKit✔️✔️:       may_need_extra_H(mol, heavyAtom)) &&
        // RDKit✔️✔️:      isAromaticAtom(*heavyAtom)) ||
        // RDKit✔️✔️:     (std::find(defaultVs.begin() + 1, defaultVs.end(),
        // RDKit✔️✔️:                heavyAtom->getTotalValence()) != defaultVs.end())) {
        // Preserve short-circuit getter order: direct and aromatic N/P paths
        // do not read either cache. The find argument reads total valence even
        // when its searched tail is empty. No whole-graph cache preparation.
        let default_valences = rdkit_valence_list(atomic_number)?;
        let aromatic_increment = (atomic_number == 7
            || atomic_number == 15
            || may_need_extra_h(topology, neighbor, valence, active_bonds)?)
            && is_aromatic_atom(topology, neighbor, active_bonds);
        if aromatic_increment {
            true
        } else {
            let total =
                crate::valence::cached_total_valence(&topology.atoms[neighbor.index()], valence)?;
            default_valences
                .and_then(|values| values.get(1..))
                .is_some_and(|values| values.contains(&total))
        }
    };
    if increment {
        let current = topology.atoms[neighbor.index()].explicit_hydrogens();
        // RDKit✔️✔️: heavyAtom->setNumExplicitHs(heavyAtom->getNumExplicitHs() + 1);
        // RDKit✔️✔️: void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
        // RDKit✔️✔️: std::uint8_t d_numExplicitHs;
        // The promoted unsigned addition narrows at the source uint8 setter.
        // This O(1) wrap is observable at 255 and is not an overflow error.
        topology.atoms[neighbor.index()].set_explicit_hydrogens(current.wrapping_add(1));
    }

    if chiral_tag != ChiralTag::Unspecified {
        let mut probe = active_incident_bonds(topology, neighbor, active_bonds);
        probe.retain(|bond| *bond != removed_bond);
        probe.push(removed_bond);
        if perturbation_order(topology, neighbor, &probe, active_bonds) % 2 == 1 {
            crate::stereo_order::invert_atom_chirality(&mut topology.atoms[neighbor.index()]);
        }
    }

    if active_degree(topology, neighbor, active_bonds) == 2
        && let Some(other) = active_incident_bonds(topology, neighbor, active_bonds)
            .into_iter()
            .find(|bond| *bond != removed_bond)
        && bond_stereo_beyond_any(topology.bonds[other.index()].stereo())
    {
        topology.bonds[other.index()]
            .set_stereo(BondStereo::None)
            .map_err(|_| HydrogenError::InvalidStereoTransition {
                bond: other,
                reason: "could not clear invalid double-bond stereo",
            })?;
        topology.bonds[other.index()].set_stereo_atoms(None);
    }

    let removed_direction = topology.bonds[removed_bond.index()].direction();
    if removed_direction == BondDirection::Unknown
        && topology.bonds[removed_bond.index()].begin() == neighbor
    {
        topology.atoms[neighbor.index()].set_unknown_stereo(true);
    } else if matches!(
        removed_direction,
        BondDirection::EndDownRight | BondDirection::EndUpRight
    ) {
        let mut found_direction = false;
        let mut other_bond = None;
        for bond in active_incident_bonds(topology, neighbor, active_bonds) {
            if bond == removed_bond || topology.bonds[bond.index()].order() != BondOrder::Single {
                continue;
            }
            if topology.bonds[bond.index()].direction() == BondDirection::None {
                other_bond = Some(bond);
            } else {
                found_direction = true;
            }
        }
        if !found_direction && let Some(other_bond) = other_bond {
            let flip = topology.bonds[other_bond.index()].begin() == neighbor
                && topology.bonds[removed_bond.index()].begin() == neighbor;
            let direction = if flip {
                match removed_direction {
                    BondDirection::EndDownRight => BondDirection::EndUpRight,
                    BondDirection::EndUpRight => BondDirection::EndDownRight,
                    _ => removed_direction,
                }
            } else {
                removed_direction
            };
            topology.bonds[other_bond.index()].set_direction(direction);
        }
    }
    adjust_stereo_atoms_if_required(
        topology,
        removed_atom,
        neighbor,
        removed_bond,
        active_atoms,
        active_bonds,
    )?;
    Ok(())
}

fn active_incident_bonds(
    topology: &TopologyBlock,
    atom: AtomId,
    active_bonds: &[bool],
) -> Vec<BondId> {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| active_bonds[neighbor.bond.index()])
        .map(|neighbor| neighbor.bond)
        .collect()
}

fn active_degree(topology: &TopologyBlock, atom: AtomId, active_bonds: &[bool]) -> usize {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .filter(|neighbor| active_bonds[neighbor.bond.index()])
        .count()
}

fn may_need_extra_h(
    topology: &TopologyBlock,
    atom: AtomId,
    valence: &ValenceAssignment,
    active_bonds: &[bool],
) -> Result<bool, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION may_need_extra_H
    // RDKit✔️✔️: bool may_need_extra_H(const ROMol &mol, const Atom *atom) {
    // RDKit✔️✔️:   unsigned single_bonds = 0;
    // RDKit✔️✔️:   unsigned aromatic_bonds = 0;
    // RDKit✔️✔️:   for (auto bond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (bond->getBondType() == Bond::SINGLE) {
    // RDKit✔️✔️:       ++single_bonds;
    // RDKit✔️✔️:     } else if (bond->getBondType() == Bond::AROMATIC) {
    // RDKit✔️✔️:       ++aromatic_bonds;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return single_bonds == 1 && aromatic_bonds == 2 &&
    // RDKit✔️✔️:          atom->getTotalValence() == 3;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION may_need_extra_H
    let mut single_bonds = 0;
    let mut aromatic_bonds = 0;
    for bond in active_incident_bonds(topology, atom, active_bonds) {
        match topology.bonds[bond.index()].order() {
            BondOrder::Single => single_bonds += 1,
            BondOrder::Aromatic => aromatic_bonds += 1,
            _ => return Ok(false),
        }
    }
    Ok(single_bonds == 1
        && aromatic_bonds == 2
        && crate::valence::cached_total_valence(&topology.atoms[atom.index()], valence)? == 3)
}

fn is_aromatic_atom(topology: &TopologyBlock, atom: AtomId, active_bonds: &[bool]) -> bool {
    // BEGIN RDKIT CPP FUNCTION isAromaticAtom
    // RDKit✔️✔️: bool isAromaticAtom(const Atom &atom) {
    // RDKit✔️✔️:   if (atom.getIsAromatic()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom.hasOwningMol()) {
    // RDKit✔️✔️:     for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:       if (bond->getIsAromatic() ||
    // RDKit✔️✔️:           bond->getBondType() == Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAromaticAtom
    topology.atoms[atom.index()].is_aromatic()
        || active_incident_bonds(topology, atom, active_bonds)
            .into_iter()
            .any(|bond| {
                topology.bonds[bond.index()].is_aromatic()
                    || topology.bonds[bond.index()].order() == BondOrder::Aromatic
            })
}

fn perturbation_order(
    topology: &TopologyBlock,
    atom: AtomId,
    probe: &[BondId],
    active_bonds: &[bool],
) -> usize {
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder
    // RDKit✔️✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
    // RDKit✔️✔️:   INT_LIST ref;
    // RDKit✔️✔️:   for (const auto bnd : getOwningMol().atomBonds(this)) {
    // RDKit✔️✔️:     ref.push_back(bnd->getIdx());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atom::getPerturbationOrder
    count_swaps_to_interconvert(probe, active_incident_bonds(topology, atom, active_bonds))
}

fn count_swaps_to_interconvert(reference: &[BondId], mut probe: Vec<BondId>) -> usize {
    // BEGIN RDKIT CPP FUNCTION countSwapsToInterconvert
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit✔️✔️:   typename T::const_iterator refIt = ref.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt = probe.begin();
    // RDKit✔️✔️:   typename T::iterator probeIt2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int nSwaps = 0;
    // RDKit✔️✔️:   while (refIt != ref.end()) {
    // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit✔️✔️:       bool foundIt = false;
    // RDKit✔️✔️:       probeIt2 = probeIt;
    // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit✔️✔️:         ++probeIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (probeIt2 != probe.end()) {
    // RDKit✔️✔️:         foundIt = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit✔️✔️:
    // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit✔️✔️:       nSwaps++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++probeIt;
    // RDKit✔️✔️:     ++refIt;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSwaps;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION countSwapsToInterconvert
    debug_assert_eq!(reference.len(), probe.len());
    let mut swaps = 0;
    for (position, expected) in reference.iter().copied().enumerate() {
        if probe[position] == expected {
            continue;
        }
        let found = probe[position..]
            .iter()
            .position(|candidate| *candidate == expected)
            .map(|offset| position + offset)
            .expect("validated bond-order permutations contain the same ids");
        probe.swap(position, found);
        swaps += 1;
    }
    swaps
}

fn adjust_stereo_atoms_if_required(
    topology: &mut TopologyBlock,
    atom: AtomId,
    heavy_atom: AtomId,
    connecting_bond: BondId,
    active_atoms: &[bool],
    active_bonds: &[bool],
) -> Result<bool, HydrogenError> {
    // BEGIN RDKIT CPP FUNCTION adjustStereoAtomsIfRequired
    // RDKit✔️✔️: bool adjustStereoAtomsIfRequired(RWMol &mol, const Atom *atom,
    // RDKit✔️✔️:                                  const Atom *heavyAtom) {
    // RDKit✔️✔️:   PRECONDITION(atom != nullptr, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(heavyAtom != nullptr, "bad heavy atom");
    // RDKit✔️✔️:   if (heavyAtom->getDegree() == 2) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const auto &cbnd =
    // RDKit✔️✔️:       mol.getBondBetweenAtoms(atom->getIdx(), heavyAtom->getIdx());
    // RDKit✔️✔️:   if (!cbnd) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &nbri :
    // RDKit✔️✔️:        boost::make_iterator_range(mol.getAtomBonds(heavyAtom))) {
    // RDKit✔️✔️:     Bond *bnd = mol[nbri];
    // RDKit✔️✔️:     if (bnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         bnd->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:       auto sAtomIt = std::find(bnd->getStereoAtoms().begin(),
    // RDKit✔️✔️:                                bnd->getStereoAtoms().end(), atom->getIdx());
    // RDKit✔️✔️:       if (sAtomIt != bnd->getStereoAtoms().end()) {
    // RDKit✔️✔️:         unsigned int dblNbrIdx = bnd->getOtherAtomIdx(heavyAtom->getIdx());
    // RDKit✔️✔️:         for (const auto &nbri :
    // RDKit✔️✔️:              boost::make_iterator_range(mol.getAtomNeighbors(heavyAtom))) {
    // RDKit✔️✔️:           const auto &nbr = mol[nbri];
    // RDKit✔️✔️:           if (nbr->getIdx() == dblNbrIdx || nbr->getIdx() == atom->getIdx()) {
    // RDKit✔️✔️:             continue;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           *sAtomIt = nbr->getIdx();
    // RDKit✔️✔️:           bool madeAdjustment = true;
    // RDKit✔️✔️:           switch (bnd->getStereo()) {
    // RDKit✔️✔️:             case Bond::STEREOCIS:
    // RDKit✔️✔️:               bnd->setStereo(Bond::STEREOTRANS);
    // RDKit✔️✔️:               break;
    // RDKit✔️✔️:             case Bond::STEREOTRANS:
    // RDKit✔️✔️:               bnd->setStereo(Bond::STEREOCIS);
    // RDKit✔️✔️:               break;
    // RDKit✔️✔️:             default:
    // RDKit✔️✔️:               madeAdjustment = false;
    // RDKit✔️✔️:               break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           return madeAdjustment;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION adjustStereoAtomsIfRequired
    if active_degree(topology, heavy_atom, active_bonds) == 2
        || !active_bonds[connecting_bond.index()]
    {
        return Ok(false);
    }
    for bond_id in active_incident_bonds(topology, heavy_atom, active_bonds) {
        let bond = &topology.bonds[bond_id.index()];
        if bond.order() != BondOrder::Double || !bond_stereo_beyond_any(bond.stereo()) {
            continue;
        }
        let Some(mut stereo_atoms) = bond.stereo_atoms() else {
            continue;
        };
        let Some(slot) = stereo_atoms
            .iter()
            .position(|stereo_atom| *stereo_atom == atom)
        else {
            continue;
        };
        let double_neighbor = if bond.begin() == heavy_atom {
            bond.end()
        } else {
            bond.begin()
        };
        let replacement = topology
            .adjacency
            .neighbors_of(heavy_atom.index())
            .iter()
            .filter(|neighbor| active_bonds[neighbor.bond.index()])
            .map(|neighbor| AtomId::new(neighbor.atom_index))
            .find(|candidate| {
                *candidate != double_neighbor
                    && *candidate != atom
                    && active_atoms[candidate.index()]
            });
        let Some(replacement) = replacement else {
            continue;
        };
        stereo_atoms[slot] = replacement;
        let old_stereo = topology.bonds[bond_id.index()].stereo();
        topology.bonds[bond_id.index()].set_stereo_atoms(Some(stereo_atoms));
        let new_stereo = match old_stereo {
            BondStereo::Cis => BondStereo::Trans,
            BondStereo::Trans => BondStereo::Cis,
            _ => return Ok(false),
        };
        topology.bonds[bond_id.index()]
            .set_stereo(new_stereo)
            .map_err(|_| HydrogenError::InvalidStereoTransition {
                bond: bond_id,
                reason: "replacement stereo atoms did not satisfy CIS/TRANS requirements",
            })?;
        return Ok(true);
    }
    Ok(false)
}

fn bond_stereo_beyond_any(stereo: BondStereo) -> bool {
    matches!(
        stereo,
        BondStereo::Z
            | BondStereo::E
            | BondStereo::Cis
            | BondStereo::Trans
            | BondStereo::AtropCw
            | BondStereo::AtropCcw
    )
}

struct HydrogenCandidateDecision {
    remove: bool,
    warning: Option<HydrogenWarning>,
}

impl HydrogenCandidateDecision {
    const fn remove() -> Self {
        Self {
            remove: true,
            warning: None,
        }
    }

    const fn keep(warning: Option<HydrogenWarning>) -> Self {
        Self {
            remove: false,
            warning,
        }
    }
}

fn should_remove_hydrogen(
    topology: &TopologyBlock,
    atom: AtomId,
    params: &RemoveHsParams,
    query_state: Option<QueryStateRef<'_>>,
) -> HydrogenCandidateDecision {
    // BEGIN RDKIT CPP FUNCTION shouldRemoveH
    // RDKit✔️❌: bool shouldRemoveH(const RWMol &mol, const Atom *atom,
    // RDKit✔️❌:                    const RemoveHsParameters &ps) {
    // RDKit✔️❌:   if (atom->getAtomicNum() != 1) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeWithQuery && atom->hasQuery()) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeDegreeZero && !atom->getDegree()) {
    // RDKit✔️❌:     if (ps.showWarnings) {
    // RDKit✔️❌:       BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:           << "WARNING: not removing hydrogen atom without neighbors"
    // RDKit✔️❌:           << std::endl;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeHigherDegrees && atom->getDegree() > 1) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeIsotopes && !ps.removeAndTrackIsotopes && atom->getIsotope()) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeNonimplicit && !atom->hasProp(common_properties::isImplicit)) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeMapped && atom->getAtomMapNum()) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (ps.removeInSGroups) {
    // RDKit✔️❌:     // If removing H in SGroups, do not remove H atoms in special
    // RDKit✔️❌:     // roles in the SGroup
    // RDKit✔️❌:     for (const auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:       // The H atom is one of the "caps" of the SGroup. Technically,
    // RDKit✔️❌:       // it's not part of the group, but it defines its boundaries.
    // RDKit✔️❌:       for (const auto &bond_idx : sg.getBonds()) {
    // RDKit✔️❌:         if (sg.getBondType(bond_idx) == SubstanceGroup::BondType::XBOND) {
    // RDKit✔️❌:           auto bond = mol.getBondWithIdx(bond_idx);
    // RDKit✔️❌:           if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
    // RDKit✔️❌:             return false;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       for (const auto &sap : sg.getAttachPoints()) {
    // RDKit✔️❌:         // The H atoms is an attach point. This would be weird, but is possible.
    // RDKit✔️❌:         // (if it is a 'leaving atom' we don't care, though)
    // RDKit✔️❌:         if (sap.aIdx == atom->getIdx()) {
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       for (const auto &cs : sg.getCStates()) {
    // RDKit✔️❌:         // The bond to the H atom defines a CState
    // RDKit✔️❌:         auto bond = mol.getBondWithIdx(cs.bondIdx);
    // RDKit✔️❌:         if (bond->getBeginAtom() == atom || bond->getEndAtom() == atom) {
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     for (const auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️❌:       if (sg.includesAtom(atom->getIdx())) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!ps.removeHydrides && atom->getFormalCharge() == -1) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   bool removeIt = true;
    // RDKit✔️❌:   if (atom->getDegree() &&
    // RDKit✔️❌:       (!ps.removeDummyNeighbors || !ps.removeDefiningBondStereo ||
    // RDKit✔️❌:        !ps.removeOnlyHNeighbors || !ps.removeNontetrahedralNeighbors ||
    // RDKit✔️❌:        !ps.removeWithWedgedBond)) {
    // RDKit✔️❌:     bool onlyHNeighbors = true;
    // RDKit✔️❌:     for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit✔️❌:       // is it a dummy?
    // RDKit✔️❌:       if (!ps.removeDummyNeighbors && nbr->getAtomicNum() < 1) {
    // RDKit✔️❌:         if (ps.showWarnings) {
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
    // RDKit✔️❌:                                        "with dummy atom neighbors"
    // RDKit✔️❌:                                     << std::endl;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // does it have non-tetrahedral stereo:
    // RDKit✔️❌:       if (!ps.removeNontetrahedralNeighbors &&
    // RDKit✔️❌:           Chirality::hasNonTetrahedralStereo(nbr)) {
    // RDKit✔️❌:         if (ps.showWarnings) {
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:               << "WARNING: not removing hydrogen atom "
    // RDKit✔️❌:                  "with neighbor that has non-tetrahedral stereochemistry"
    // RDKit✔️❌:               << std::endl;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!ps.removeOnlyHNeighbors && nbr->getAtomicNum() != 1) {
    // RDKit✔️❌:         onlyHNeighbors = false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (!ps.removeWithWedgedBond) {
    // RDKit✔️❌:         const auto bnd = mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx());
    // RDKit✔️❌:         if (bnd->getBondDir() == Bond::BEGINDASH ||
    // RDKit✔️❌:             bnd->getBondDir() == Bond::BEGINWEDGE) {
    // RDKit✔️❌:           if (ps.showWarnings) {
    // RDKit✔️❌:             BOOST_LOG(rdWarningLog) << "WARNING: not removing hydrogen atom "
    // RDKit✔️❌:                                        "with wedged bond"
    // RDKit✔️❌:                                     << std::endl;
    // RDKit✔️❌:           }
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       // Check to see if the neighbor has a double bond and we're the only
    // RDKit✔️❌:       // neighbor at this end.  This was part of github #1810
    // RDKit✔️❌:       if (!ps.removeDefiningBondStereo && nbr->getDegree() == 2) {
    // RDKit✔️❌:         for (const auto bnd : mol.atomBonds(nbr)) {
    // RDKit✔️❌:           if (bnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️❌:               (bnd->getStereo() > Bond::STEREOANY ||
    // RDKit✔️❌:                mol.getBondBetweenAtoms(atom->getIdx(), nbr->getIdx())
    // RDKit✔️❌:                        ->getBondDir() > Bond::NONE)) {
    // RDKit✔️❌:             return false;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (removeIt && (!ps.removeOnlyHNeighbors && onlyHNeighbors)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return removeIt;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION shouldRemoveH
    let atom_data = &topology.atoms[atom.index()];
    if atom_data.atomic_number() != 1 {
        return HydrogenCandidateDecision::keep(None);
    }
    // BEGIN RDKIT CPP FUNCTION shouldRemoveH typed query guard
    // RDKit✔️✔️:   if (!ps.removeWithQuery && atom->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION shouldRemoveH typed query guard
    // Behavior review: only the hydrogen atom's typed source dynamic type is
    // tested here; incident QueryBonds deliberately do not protect it.
    // Complexity review: the typed identity lookup is allocation-free O(1),
    // matching the source pointer test.
    if !params.remove_with_query && query_state.is_some_and(|state| state.atom_has_query(atom)) {
        return HydrogenCandidateDecision::keep(None);
    }
    let neighbors = topology.adjacency.neighbors_of(atom.index());
    if !params.remove_degree_zero && neighbors.is_empty() {
        let warning = params
            .show_warnings
            .then_some(HydrogenWarning::IsolatedHydrogen { hydrogen: atom });
        return HydrogenCandidateDecision::keep(warning);
    }
    if !params.remove_higher_degrees && neighbors.len() > 1 {
        return HydrogenCandidateDecision::keep(None);
    }
    if !params.remove_isotopes
        && !params.remove_and_track_isotopes
        && atom_data.isotope().unwrap_or(0) != 0
    {
        return HydrogenCandidateDecision::keep(None);
    }
    if !params.remove_nonimplicit && !atom_data.implicit_hydrogen() {
        return HydrogenCandidateDecision::keep(None);
    }
    if !params.remove_mapped && atom_data.atom_map().unwrap_or(0) != 0 {
        return HydrogenCandidateDecision::keep(None);
    }

    if params.remove_in_sgroups {
        for group in &topology.substance_groups {
            for bond_id in group.bonds() {
                if group.bond_role(*bond_id) == SGroupBondRole::Crossing {
                    let bond = &topology.bonds[bond_id.index()];
                    if bond.begin() == atom || bond.end() == atom {
                        return HydrogenCandidateDecision::keep(None);
                    }
                }
            }
            if group.attach_points().iter().any(|point| point.atom == atom) {
                return HydrogenCandidateDecision::keep(None);
            }
            for cstate in group.cstates() {
                let bond = &topology.bonds[cstate.bond.index()];
                if bond.begin() == atom || bond.end() == atom {
                    return HydrogenCandidateDecision::keep(None);
                }
            }
        }
    } else if topology
        .substance_groups
        .iter()
        .any(|group| sgroup_includes_atom(group, atom))
    {
        return HydrogenCandidateDecision::keep(None);
    }
    if !params.remove_hydrides && atom_data.formal_charge() == -1 {
        return HydrogenCandidateDecision::keep(None);
    }

    if !neighbors.is_empty()
        && (!params.remove_dummy_neighbors
            || !params.remove_defining_bond_stereo
            || !params.remove_only_h_neighbors
            || !params.remove_nontetrahedral_neighbors
            || !params.remove_with_wedged_bond)
    {
        let mut only_h_neighbors = true;
        for neighbor in neighbors {
            let neighbor_atom = &topology.atoms[neighbor.atom_index];
            if !params.remove_dummy_neighbors && neighbor_atom.atomic_number() < 1 {
                let warning = params
                    .show_warnings
                    .then_some(HydrogenWarning::DummyAtomNeighbor {
                        hydrogen: atom,
                        neighbor: AtomId::new(neighbor.atom_index),
                    });
                return HydrogenCandidateDecision::keep(warning);
            }
            if !params.remove_nontetrahedral_neighbors
                && has_nontetrahedral_stereo(neighbor_atom.chiral_tag())
            {
                let warning =
                    params
                        .show_warnings
                        .then_some(HydrogenWarning::NonTetrahedralStereoNeighbor {
                            hydrogen: atom,
                            neighbor: AtomId::new(neighbor.atom_index),
                        });
                return HydrogenCandidateDecision::keep(warning);
            }
            if !params.remove_only_h_neighbors && neighbor_atom.atomic_number() != 1 {
                only_h_neighbors = false;
            }
            let hydrogen_bond = &topology.bonds[neighbor.bond.index()];
            if !params.remove_with_wedged_bond
                && matches!(
                    hydrogen_bond.direction(),
                    BondDirection::BeginDash | BondDirection::BeginWedge
                )
            {
                let warning = params.show_warnings.then_some(HydrogenWarning::WedgedBond {
                    hydrogen: atom,
                    bond: neighbor.bond,
                });
                return HydrogenCandidateDecision::keep(warning);
            }
            if !params.remove_defining_bond_stereo
                && topology.adjacency.neighbors_of(neighbor.atom_index).len() == 2
            {
                for adjacent in topology.adjacency.neighbors_of(neighbor.atom_index) {
                    let bond = &topology.bonds[adjacent.bond.index()];
                    if bond.order() == BondOrder::Double
                        && (bond.stereo().rdkit_code() > BondStereo::Any.rdkit_code()
                            || hydrogen_bond.direction() != BondDirection::None)
                    {
                        return HydrogenCandidateDecision::keep(None);
                    }
                }
            }
        }
        if !params.remove_only_h_neighbors && only_h_neighbors {
            return HydrogenCandidateDecision::keep(None);
        }
    }
    HydrogenCandidateDecision::remove()
}

fn sgroup_includes_atom(group: &SubstanceGroup, atom: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::includesAtom
    // RDKit✔️✔️: bool SubstanceGroup::includesAtom(unsigned int atomIdx) const {
    // RDKit✔️✔️:   if (std::find(d_atoms.begin(), d_atoms.end(), atomIdx) != d_atoms.end()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (std::find(d_patoms.begin(), d_patoms.end(), atomIdx) != d_patoms.end()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto &ap : d_saps) {
    // RDKit✔️✔️:     if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION SubstanceGroup::includesAtom
    group.atoms().contains(&atom)
        || group.parent_atoms().contains(&atom)
        || group
            .attach_points()
            .iter()
            .any(|point| point.atom == atom || point.leaving_atom == Some(atom))
}

const fn has_nontetrahedral_stereo(tag: ChiralTag) -> bool {
    // BEGIN RDKIT CPP FUNCTION Chirality::hasNonTetrahedralStereo
    // RDKit✔️✔️: bool hasNonTetrahedralStereo(const Atom *cen) {
    // RDKit✔️✔️:   PRECONDITION(cen, "bad center pointer");
    // RDKit✔️✔️:   if (!cen->hasOwningMol()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto tag = cen->getChiralTag();
    // RDKit✔️✔️:   return tag == Atom::ChiralType::CHI_SQUAREPLANAR ||
    // RDKit✔️✔️:          tag == Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL ||
    // RDKit✔️✔️:          tag == Atom::ChiralType::CHI_OCTAHEDRAL;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Chirality::hasNonTetrahedralStereo
    matches!(
        tag,
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    )
}

fn filter_sgroup_emptying_hydrogens(topology: &TopologyBlock, selected: &mut [bool]) {
    // BEGIN RDKIT CPP FUNCTION filter_sgroup_emptying_hydrogens
    // RDKit✔️✔️: void filter_sgroup_emptying_hydrogens(const ROMol &mol,
    // RDKit✔️✔️:                                       boost::dynamic_bitset<> &atomsToRemove) {
    // RDKit✔️✔️:   for (const auto &sg : getSubstanceGroups(mol)) {
    // RDKit✔️✔️:     const auto &atoms = sg.getAtoms();
    // RDKit✔️✔️:     const auto &patoms = sg.getParentAtoms();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // If the SGroup already didn't have atoms, we don't care about it
    // RDKit✔️✔️:     if (atoms.empty() && patoms.empty()) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto would_remove_atom = [&atomsToRemove](const auto idx) {
    // RDKit✔️✔️:       return atomsToRemove[idx];
    // RDKit✔️✔️:     };
    // RDKit✔️✔️:
    // RDKit✔️✔️:     auto no_atoms = atoms.empty() ||
    // RDKit✔️✔️:                     std::all_of(atoms.begin(), atoms.end(), would_remove_atom);
    // RDKit✔️✔️:     if (no_atoms) {
    // RDKit✔️✔️:       auto no_patoms =
    // RDKit✔️✔️:           patoms.empty() ||
    // RDKit✔️✔️:           std::all_of(patoms.begin(), patoms.end(), would_remove_atom);
    // RDKit✔️✔️:       if (no_patoms) {
    // RDKit✔️✔️:         for (auto atom : atoms) {
    // RDKit✔️✔️:           atomsToRemove.set(atom, false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         for (auto patom : patoms) {
    // RDKit✔️✔️:           atomsToRemove.set(patom, false);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION filter_sgroup_emptying_hydrogens
    for group in &topology.substance_groups {
        let atoms = group.atoms();
        let parent_atoms = group.parent_atoms();
        if atoms.is_empty() && parent_atoms.is_empty() {
            continue;
        }
        let no_atoms = atoms.is_empty() || atoms.iter().all(|atom| selected[atom.index()]);
        if no_atoms {
            let no_parent_atoms =
                parent_atoms.is_empty() || parent_atoms.iter().all(|atom| selected[atom.index()]);
            if no_parent_atoms {
                for atom in atoms.iter().chain(parent_atoms) {
                    selected[atom.index()] = false;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, BondQueryPredicate, QueryAtom, QueryBond, QueryNode, QueryStateRef,
    };

    fn q05_query_rows(
        topology: &TopologyBlock,
        explicit_atoms: &[usize],
        explicit_bonds: &[usize],
    ) -> (Vec<QueryAtom>, Vec<QueryBond>) {
        let atoms = topology
            .atoms
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, atom)| {
                let predicate =
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()));
                if explicit_atoms.contains(&index) {
                    QueryAtom::from_parts(atom, predicate)
                } else {
                    QueryAtom::from_carrier_parts(atom, predicate)
                }
            })
            .collect();
        let bonds = topology
            .bonds
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, bond)| {
                let predicate = QueryNode::predicate(if explicit_bonds.contains(&index) {
                    BondQueryPredicate::Any
                } else {
                    BondQueryPredicate::Order(bond.order())
                });
                if explicit_bonds.contains(&index) {
                    QueryBond::from_parts(bond, predicate)
                } else {
                    QueryBond::from_carrier_parts(bond, predicate)
                }
            })
            .collect();
        (atoms, bonds)
    }

    #[test]
    fn q05_hydrogen_query_identity_retains_explicit_query_h_and_removes_carrier_h() {
        let topology = bonded_ch_topology();
        let (explicit_atoms, explicit_bonds) = q05_query_rows(&topology, &[1], &[0]);
        let explicit_state =
            QueryStateRef::try_for_topology(&explicit_atoms, &explicit_bonds, &topology).unwrap();
        assert!(explicit_state.atom_has_query(AtomId::new(1)));
        assert_eq!(
            remove_hydrogen_candidates_with_query_state(
                &topology,
                &RemoveHsParams::default(),
                Some(explicit_state),
            ),
            Vec::<AtomId>::new(),
            "RDKit shouldRemoveH retains an explicit QueryAtom hydrogen"
        );

        let (carrier_atoms, carrier_bonds) = q05_query_rows(&topology, &[], &[]);
        let carrier_state =
            QueryStateRef::try_for_topology(&carrier_atoms, &carrier_bonds, &topology).unwrap();
        assert!(!carrier_state.atom_has_query(AtomId::new(1)));
        assert_eq!(
            remove_hydrogen_candidates_with_query_state(
                &topology,
                &RemoveHsParams::default(),
                Some(carrier_state),
            ),
            vec![AtomId::new(1)],
            "an ordinary hydrogen represented by a carrier-derived row remains removable"
        );
    }

    #[test]
    fn q05_hydrogen_query_identity_add_hs_skips_explicit_atom_query() {
        let topology = explicit_h_topology();
        let (query_atoms, query_bonds) = q05_query_rows(&topology, &[0], &[]);
        let state = QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &topology).unwrap();
        assert!(state.atom_has_query(AtomId::new(0)));

        let result = add_hydrogens_topology_with_query_state(
            topology,
            &AddHsParams {
                explicit_only: true,
                skip_queries: true,
                ..Default::default()
            },
            Some(state),
        )
        .unwrap();
        assert_eq!(
            result.topology.atoms.len(),
            1,
            "RDKit AddHs(skipQueries=true) must not add H to an explicit QueryAtom"
        );
    }

    #[test]
    fn q05_hydrogen_query_identity_add_hs_skips_explicit_incident_bond_query() {
        let mut topology = bonded_ch_topology();
        topology.atoms[0].set_explicit_hydrogens(1);
        let (query_atoms, query_bonds) = q05_query_rows(&topology, &[], &[0]);
        let state = QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &topology).unwrap();
        assert!(!state.atom_has_query(AtomId::new(0)));
        assert!(state.bond_has_query(BondId::new(0)));

        let result = add_hydrogens_topology_with_query_state(
            topology,
            &AddHsParams {
                explicit_only: true,
                skip_queries: true,
                only_on_atoms: Some(vec![AtomId::new(0)]),
                ..Default::default()
            },
            Some(state),
        )
        .unwrap();
        assert_eq!(
            result.topology.atoms.len(),
            2,
            "RDKit AddHs(skipQueries=true) must skip an atom incident to an explicit QueryBond"
        );
    }

    fn one_hydrogen_addition_result() -> AddHydrogensTopologyResult {
        add_hydrogens_topology(
            explicit_h_topology(),
            &AddHsParams {
                explicit_only: true,
                ..Default::default()
            },
        )
        .expect("one explicit hydrogen topology result")
    }

    fn assert_point_close(actual: AddHsPoint3D, expected: AddHsPoint3D) {
        for (actual, expected) in [
            (actual.x, expected.x),
            (actual.y, expected.y),
            (actual.z, expected.z),
        ] {
            assert!(
                (actual - expected).abs() < 1.0e-12,
                "{actual} != {expected}"
            );
        }
    }

    fn explicit_h_topology() -> TopologyBlock {
        let carbon = cosmolkit_model::Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_explicit_hydrogens(1),
        );
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![carbon],
            ..TopologyBlock::default()
        }
    }

    fn bonded_ch_topology() -> TopologyBlock {
        let atoms = vec![
            cosmolkit_model::Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            cosmolkit_model::Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::H)),
        ];
        let bonds = vec![Bond::from_spec(
            cosmolkit_model::BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn parent_with_neighbor_topology(
        hybridization: Hybridization,
        mut parent_spec: AtomSpec,
        mut neighbor_spec: AtomSpec,
        bond_spec: BondSpec,
    ) -> TopologyBlock {
        parent_spec = parent_spec
            .with_hybridization(hybridization)
            .with_explicit_hydrogens(1);
        if bond_spec.is_aromatic() {
            parent_spec = parent_spec.with_aromatic(true);
            neighbor_spec = neighbor_spec.with_aromatic(true);
        }
        let atoms = vec![
            cosmolkit_model::Atom::from_spec(AtomId::new(0), parent_spec),
            cosmolkit_model::Atom::from_spec(AtomId::new(1), neighbor_spec),
        ];
        let bonds = vec![Bond::from_spec(BondId::new(0), bond_spec)];
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn explicit_only_addition(topology: TopologyBlock) -> AddHydrogensTopologyResult {
        add_hydrogens_topology(
            topology,
            &AddHsParams {
                explicit_only: true,
                ..Default::default()
            },
        )
        .expect("valid explicit-only topology addition")
    }

    fn explicit_h_star_topology(
        other_neighbor_count: usize,
        parent_spec: AtomSpec,
    ) -> TopologyBlock {
        let mut atoms = vec![cosmolkit_model::Atom::from_spec(
            AtomId::new(0),
            parent_spec.with_explicit_hydrogens(1),
        )];
        let mut bonds = Vec::with_capacity(other_neighbor_count);
        for neighbor_index in 0..other_neighbor_count {
            let neighbor = AtomId::new(neighbor_index + 1);
            atoms.push(cosmolkit_model::Atom::from_spec(
                neighbor,
                AtomSpec::new(Element::C),
            ));
            bonds.push(Bond::from_spec(
                BondId::new(neighbor_index),
                BondSpec::new(AtomId::new(0), neighbor, BondOrder::Single),
            ));
        }
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    fn residue_topology(
        parent_infos: Vec<AtomPdbResidueInfo>,
        hydrogens: Vec<(usize, Option<AtomPdbResidueInfo>)>,
    ) -> TopologyBlock {
        let parent_count = parent_infos.len();
        let mut atoms = parent_infos
            .into_iter()
            .enumerate()
            .map(|(index, info)| {
                cosmolkit_model::Atom::from_spec(
                    AtomId::new(index),
                    AtomSpec::new(Element::C).with_pdb_residue_info(info),
                )
            })
            .collect::<Vec<_>>();
        let mut bonds = Vec::with_capacity(hydrogens.len());
        for (hydrogen_offset, (parent_index, info)) in hydrogens.into_iter().enumerate() {
            let hydrogen_id = AtomId::new(parent_count + hydrogen_offset);
            let mut spec = AtomSpec::new(Element::H);
            if let Some(info) = info {
                spec = spec.with_pdb_residue_info(info);
            }
            atoms.push(cosmolkit_model::Atom::from_spec(hydrogen_id, spec));
            bonds.push(Bond::from_spec(
                BondId::new(hydrogen_offset),
                BondSpec::new(AtomId::new(parent_index), hydrogen_id, BondOrder::Single),
            ));
        }
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    #[test]
    fn h_add_coordinates_residue_scans_final_serials_and_uses_source_defaults() {
        let negative = AtomPdbResidueInfo::new(" C1 ", -9, "NEG", 10, "A", false)
            .with_alt_loc("X")
            .with_insertion_code("I")
            .with_occupancy(0.25)
            .with_temp_factor(12.5)
            .with_secondary_structure(3)
            .with_segment_number(4)
            .with_monomer_class("source-class");
        let zero = AtomPdbResidueInfo::new(" C2 ", 0, "ZERO", 10, "A", true);
        let maximum = AtomPdbResidueInfo::new(" C3 ", 42, "MAX", 11, "B", true);
        let mut topology = residue_topology(
            vec![negative.clone(), zero, maximum],
            vec![(0, None), (1, None), (2, None)],
        );

        assign_hydrogen_residue_info(&mut topology).expect("valid residue assignment");

        let first = topology.atoms[3]
            .pdb_residue_info()
            .expect("first hydrogen info");
        assert_eq!(first.atom_name(), " H1 ");
        assert_eq!(first.serial_number(), 42);
        assert_eq!(first.residue_name(), "NEG");
        assert_eq!(first.residue_number(), 10);
        assert_eq!(first.chain_id(), "A");
        assert!(!first.is_hetero_atom());
        assert_eq!(first.alt_loc(), "");
        assert_eq!(first.insertion_code(), "");
        assert_eq!(first.occupancy(), 1.0);
        assert_eq!(first.temp_factor(), 0.0);
        assert_eq!(first.secondary_structure(), 0);
        assert_eq!(first.segment_number(), 0);
        assert_eq!(first.monomer_class(), "");
        assert_eq!(
            topology.atoms[4]
                .pdb_residue_info()
                .unwrap()
                .serial_number(),
            43
        );
        assert_eq!(
            topology.atoms[4].pdb_residue_info().unwrap().atom_name(),
            " H2 "
        );
        assert_eq!(
            topology.atoms[5]
                .pdb_residue_info()
                .unwrap()
                .serial_number(),
            44
        );
        assert_eq!(
            topology.atoms[5].pdb_residue_info().unwrap().atom_name(),
            " H1 "
        );
        assert_eq!(topology.atoms[0].pdb_residue_info(), Some(&negative));
    }

    #[test]
    fn h_add_coordinates_residue_consumes_existing_ids_and_resets_by_residue_chain() {
        let existing = AtomPdbResidueInfo::new(" HX ", 90, "ALA", 1, "A", false)
            .with_alt_loc("E")
            .with_temp_factor(7.0);
        let mut topology = residue_topology(
            vec![
                AtomPdbResidueInfo::new(" C1 ", 1, "ALA", 1, "A", false),
                AtomPdbResidueInfo::new(" C2 ", 2, "GLY", 1, "A", false),
                AtomPdbResidueInfo::new(" C3 ", 3, "SER", 2, "A", false),
                AtomPdbResidueInfo::new(" C4 ", 4, "SER", 2, "B", false),
            ],
            vec![
                (0, None),
                (0, Some(existing.clone())),
                (0, None),
                (1, None),
                (2, None),
                (3, None),
            ],
        );

        assign_hydrogen_residue_info(&mut topology).expect("ordered residue assignment");

        assert_eq!(
            topology.atoms[4].pdb_residue_info().unwrap().atom_name(),
            " H1 "
        );
        assert_eq!(
            topology.atoms[4]
                .pdb_residue_info()
                .unwrap()
                .serial_number(),
            90
        );
        assert_eq!(topology.atoms[5].pdb_residue_info(), Some(&existing));
        assert_eq!(
            topology.atoms[6].pdb_residue_info().unwrap().atom_name(),
            " H3 "
        );
        assert_eq!(
            topology.atoms[6]
                .pdb_residue_info()
                .unwrap()
                .serial_number(),
            91
        );
        assert_eq!(
            topology.atoms[7].pdb_residue_info().unwrap().atom_name(),
            " H4 "
        );
        assert_eq!(
            topology.atoms[7].pdb_residue_info().unwrap().residue_name(),
            "GLY"
        );
        assert_eq!(
            topology.atoms[8].pdb_residue_info().unwrap().atom_name(),
            " H1 "
        );
        assert_eq!(
            topology.atoms[9].pdb_residue_info().unwrap().atom_name(),
            " H1 "
        );
    }

    #[test]
    fn h_add_coordinates_residue_wraps_three_and_more_digit_names_in_adjacency_order() {
        let mut topology = residue_topology(
            vec![AtomPdbResidueInfo::new(" C  ", 0, "BIG", 7, "Z", false)],
            (0..1_000).map(|_| (0, None)).collect(),
        );

        assign_hydrogen_residue_info(&mut topology).expect("large residue assignment");

        assert_eq!(
            topology.atoms[1].pdb_residue_info().unwrap().atom_name(),
            " H1 "
        );
        assert_eq!(
            topology.atoms[123].pdb_residue_info().unwrap().atom_name(),
            "3H12"
        );
        assert_eq!(
            topology.atoms[1_000]
                .pdb_residue_info()
                .unwrap()
                .atom_name(),
            "0H00"
        );
        assert!(
            topology.atoms[1..].iter().all(|atom| atom
                .pdb_residue_info()
                .unwrap()
                .atom_name()
                .len()
                == 4)
        );
    }

    #[test]
    fn h_add_coordinates_residue_flag_installs_typed_info_on_appended_hydrogen() {
        let parent_info = AtomPdbResidueInfo::new(" C  ", 12, "LIG", 8, "Q", true);
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![cosmolkit_model::Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_explicit_hydrogens(1)
                    .with_pdb_residue_info(parent_info),
            )],
            ..TopologyBlock::default()
        };
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock::default(),
            false,
            true,
        )
        .expect("residue-info stage");

        let hydrogen = output.topology.atoms[1].pdb_residue_info().unwrap();
        assert_eq!(hydrogen.atom_name(), " H1 ");
        assert_eq!(hydrogen.serial_number(), 12);
        assert_eq!(hydrogen.residue_name(), "LIG");
        assert_eq!(hydrogen.residue_number(), 8);
        assert_eq!(hydrogen.chain_id(), "Q");
        assert!(hydrogen.is_hetero_atom());
    }

    #[test]
    fn default_add_hydrogens_uses_detached_valence_boundary() {
        let result = add_hydrogens_impl(
            TopologyBlock::default(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .expect("default AddHs on empty topology");
        assert!(result.topology.atoms.is_empty());
        assert!(result.warnings.is_empty());
    }

    #[test]
    fn explicit_add_hydrogens_is_a_detached_topology_transform() {
        let output = add_hydrogens_with_params(
            explicit_h_topology(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &AddHsParams {
                explicit_only: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(output.topology.atoms.len(), 2);
        assert_eq!(output.topology.bonds.len(), 1);
        assert!(output.coordinates.conformers_2d.is_empty());
        assert_eq!(output.topology.atoms[0].explicit_hydrogens(), 0);
        assert!(output.topology.validate().is_ok());
    }

    #[test]
    fn remove_hydrogens_compacts_detached_blocks_and_remaps_coordinates() {
        let result = remove_hydrogens_impl(
            bonded_ch_topology(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
        )
        .unwrap();
        assert_eq!(result.topology.atoms.len(), 1);
        assert!(result.topology.bonds.is_empty());
        assert!(result.coordinates.conformers_3d.is_empty());
        assert!(result.topology.validate().is_ok());
    }

    #[test]
    fn h_add_coordinates_foundation_grows_every_conformer_and_preserves_metadata() {
        let coordinates = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[1.25, -2.5]]).with_prop("source", "two")],
            conformers_3d: vec![
                Conformer3D::new(11, vec![[3.0, 4.0, 5.0]], true).with_prop("source", "three"),
                Conformer3D::new(12, vec![[-1.0, 2.0, 0.0]], false).with_prop("flat", "yes"),
            ],
            source_conformer_order: None,
            source_coordinate_dim: Some(cosmolkit_model::CoordinateDimension::ThreeD),
        };
        let original = coordinates.clone();
        let output =
            add_hydrogen_coordinates(one_hydrogen_addition_result(), coordinates, false, false)
                .expect("disabled placement still grows conformers");

        assert_eq!(
            output.coordinates.source_coordinate_dim,
            original.source_coordinate_dim
        );
        assert_eq!(output.coordinates.conformers_2d[0].id(), 7);
        assert_eq!(
            output.coordinates.conformers_2d[0].props(),
            original.conformers_2d[0].props()
        );
        assert_eq!(
            output.coordinates.conformers_2d[0].coordinates(),
            &[[1.25, -2.5], [0.0, 0.0]]
        );
        assert_eq!(output.coordinates.conformers_3d[0].id(), 11);
        assert!(output.coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            output.coordinates.conformers_3d[0].coordinates(),
            &[[3.0, 4.0, 5.0], [0.0, 0.0, 0.0]]
        );
        assert!(!output.coordinates.conformers_3d[1].is_3d());
        assert_eq!(
            output.coordinates.conformers_3d[1].props(),
            original.conformers_3d[1].props()
        );
        assert_eq!(original.conformers_2d[0].coordinates(), &[[1.25, -2.5]]);
    }

    #[test]
    fn h_add_coordinates_foundation_rejects_old_row_mismatch_and_bad_addition_parent() {
        let row_error = add_hydrogen_coordinates(
            one_hydrogen_addition_result(),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(3, Vec::new())],
                ..CoordinateBlock::default()
            },
            false,
            false,
        )
        .expect_err("old coordinate row mismatch must fail before growth");
        assert!(matches!(
            row_error,
            HydrogenError::InvalidCoordinates(CoordinateValidationError::RowCount {
                dimension: "2D",
                conformer: 3,
                rows: 0,
                atom_count: 1,
            })
        ));

        let mut bad_plan = one_hydrogen_addition_result();
        bad_plan.additions[0].parent = AtomId::new(1);
        let plan_error =
            add_hydrogen_coordinates(bad_plan, CoordinateBlock::default(), false, false)
                .expect_err("an appended atom cannot be its addition parent");
        assert_eq!(
            plan_error,
            HydrogenError::InvalidAdditionPlan {
                addition: Some(0),
                reason: "parent is not an original atom",
            }
        );
    }

    #[test]
    fn h_add_coordinates_foundation_vector_zero_and_perpendicular_match_source() {
        assert_eq!(AddHsPoint3D::ZERO.normalized(), None);
        assert_eq!(AddHsPoint3D::ZERO.perpendicular(), None);
        assert_point_close(
            AddHsPoint3D::new(2.0, 0.0, 0.0)
                .perpendicular()
                .expect("x axis perpendicular"),
            AddHsPoint3D::new(0.0, 1.0, 0.0),
        );
        assert_point_close(
            AddHsPoint3D::new(2.0, 3.0, 0.0)
                .perpendicular()
                .expect("xy perpendicular"),
            AddHsPoint3D::new(3.0, -2.0, 0.0).normalized().unwrap(),
        );
        assert_point_close(
            AddHsPoint3D::new(0.0, 2.0, 3.0)
                .perpendicular()
                .expect("yz perpendicular"),
            AddHsPoint3D::new(0.0, 3.0, -2.0).normalized().unwrap(),
        );
    }

    #[test]
    fn h_add_coordinates_foundation_axis_rotation_matches_source_matrix() {
        assert_point_close(
            AddHsPoint3D::new(1.0, 0.0, 0.0)
                .rotated(
                    std::f64::consts::FRAC_PI_2,
                    AddHsPoint3D::new(0.0, 0.0, 4.0),
                )
                .expect("nonzero rotation axis"),
            AddHsPoint3D::new(0.0, 1.0, 0.0),
        );
        assert_eq!(
            AddHsPoint3D::new(1.0, 0.0, 0.0).rotated(1.0, AddHsPoint3D::ZERO),
            None
        );
        let a = AddHsPoint3D::new(1.0, 2.0, 3.0);
        let b = AddHsPoint3D::new(-4.0, 5.0, -6.0);
        assert_eq!(a.plus(b), AddHsPoint3D::new(-3.0, 7.0, -3.0));
        assert_eq!(a.minus(b), AddHsPoint3D::new(5.0, -3.0, 9.0));
        assert_eq!(a.dot(b), -12.0);
        assert_eq!(a.cross(b), AddHsPoint3D::new(-27.0, -6.0, 13.0));
    }

    #[test]
    fn h_add_coordinates_degree_one_two_terminal_axis_uses_unit_and_rb0_distance() {
        let topology = explicit_h_topology();
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(1, vec![[2.0, 3.0]])],
                conformers_3d: vec![
                    Conformer3D::new(2, vec![[4.0, 5.0, 6.0]], true),
                    Conformer3D::new(3, vec![[7.0, 8.0, 9.0]], false),
                ],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("degree-one placement");
        assert_eq!(
            output.coordinates.conformers_2d[0].coordinates()[1],
            [3.0, 3.0]
        );
        let bond_length = rdkit_rb0(1) + rdkit_rb0(6);
        assert_point_close(
            AddHsPoint3D::from_3d(output.coordinates.conformers_3d[0].coordinates()[1]),
            AddHsPoint3D::new(4.0, 5.0, 6.0 + bond_length),
        );
        assert_eq!(
            output.coordinates.conformers_3d[1].coordinates()[1],
            [8.0, 8.0, 9.0]
        );
    }

    #[test]
    fn terminal_coordinate_approved_divergence_is_independent_of_conformer_order() {
        // CK-COORD-001: pinned reference reproduction and raw output are in
        // IO-mol_post.md, Step 2583. These are intentional-difference tests,
        // not a claim that the mixed-conformer RDKit oracle passes.
        let result = explicit_only_addition(explicit_h_topology());
        let inputs = [
            (false, AddHsPoint3D::new(1.0, 2.0, 0.0)),
            (true, AddHsPoint3D::new(2.0, 2.0, 0.0)),
            (false, AddHsPoint3D::new(4.0, 5.0, 7.0)),
        ];
        let distance = rdkit_rb0(1) + rdkit_rb0(6);
        let expected = [
            AddHsPoint3D::new(2.0, 2.0, 0.0),
            AddHsPoint3D::new(2.0, 2.0, distance),
            AddHsPoint3D::new(5.0, 5.0, 7.0),
        ];
        let place = |index: usize| {
            let (flag, parent) = inputs[index];
            terminal_position(
                &result.topology,
                &result.additions[0],
                0,
                index,
                "XYZ",
                flag,
                &[parent, AddHsPoint3D::new(0.0, 0.0, 0.0)],
            )
            .expect("degree-one placement")
        };
        let isolated = [place(0), place(1), place(2)];
        assert_eq!(isolated, expected);
        for order in [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ] {
            for index in order {
                assert_eq!(place(index), isolated[index]);
            }
        }
        // Explicitly retain disagreement with BOTH source insertion orders:
        // false->true contaminates true's X; true->false contaminates false's Z.
        assert_ne!(isolated[1], AddHsPoint3D::new(3.1, 2.0, 1.1));
        assert_ne!(isolated[0], AddHsPoint3D::new(2.0, 2.0, 1.0));
    }

    #[test]
    fn h_add_coordinates_degree_one_two_sequential_sp3_uses_prior_hydrogen() {
        let carbon = cosmolkit_model::Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_hybridization(Hybridization::Sp3)
                .with_explicit_hydrogens(2),
        );
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(1, &[]),
            atoms: vec![carbon],
            ..TopologyBlock::default()
        };
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0]])],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("first and second H use active degrees one and two");
        let coordinates = output.coordinates.conformers_2d[0].coordinates();
        assert_eq!(coordinates[1], [1.0, 0.0]);
        let first = AddHsPoint3D::from_2d(coordinates[1]);
        let second = AddHsPoint3D::from_2d(coordinates[2]);
        assert!((first.dot(second) - 109.471_f64.to_radians().cos()).abs() < 1.0e-6);
        assert!((second.length_sq() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn h_add_coordinates_degree_one_two_coincident_neighbor_matches_issue_678() {
        let topology = parent_with_neighbor_topology(
            Hybridization::Sp3,
            AtomSpec::new(Element::C),
            AtomSpec::new(Element::C),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![[2.0, -1.0, 4.0], [2.0, -1.0, 4.0]],
                    true,
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("coincident neighbor placement");
        assert_eq!(
            output.coordinates.conformers_3d[0].coordinates()[2],
            [2.0, -1.0, 4.0]
        );
    }

    #[test]
    fn h_add_coordinates_degree_one_two_sp_and_default_follow_away_vector() {
        for hybridization in [Hybridization::Sp, Hybridization::Unspecified] {
            let topology = parent_with_neighbor_topology(
                hybridization,
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            );
            let output = add_hydrogen_coordinates(
                explicit_only_addition(topology),
                CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(0, vec![[0.0, 0.0], [1.0, 0.0]])],
                    ..CoordinateBlock::default()
                },
                true,
                false,
            )
            .expect("linear/default placement");
            assert_eq!(
                output.coordinates.conformers_2d[0].coordinates()[2],
                [-1.0, 0.0]
            );
        }
    }

    #[test]
    fn h_add_coordinates_degree_one_two_sp2_uses_each_source_plane_trigger() {
        let bond_specs = [
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single).with_conjugated(true),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Aromatic).with_aromatic(true),
        ];
        for bond_spec in bond_specs {
            let mut topology = parent_with_neighbor_topology(
                Hybridization::Sp2,
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                bond_spec,
            );
            let second_neighbor = AtomId::new(2);
            topology.atoms.push(cosmolkit_model::Atom::from_spec(
                second_neighbor,
                AtomSpec::new(Element::C),
            ));
            topology.bonds.push(Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), second_neighbor, BondOrder::Single),
            ));
            topology.adjacency =
                AdjacencyList::from_topology(topology.atoms.len(), &topology.bonds);
            let output = add_hydrogen_coordinates(
                explicit_only_addition(topology),
                CoordinateBlock {
                    conformers_3d: vec![Conformer3D::new(
                        0,
                        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
                        true,
                    )],
                    ..CoordinateBlock::default()
                },
                true,
                false,
            )
            .expect("SP2 local-plane placement");
            let hydrogen = output.coordinates.conformers_3d[0].coordinates()[3];
            assert!(hydrogen[1] < 0.0);
            assert!(hydrogen[2].abs() < 1.0e-12);
            let distance = AddHsPoint3D::from_3d(hydrogen).length_sq().sqrt();
            assert!((distance - (rdkit_rb0(1) + rdkit_rb0(6))).abs() < 1.0e-12);
        }
    }

    #[test]
    fn h_add_coordinates_degree_three_four_degree_three_sp2_sp3_and_default() {
        let flat_expected = -0.5_f64.sqrt();
        for hybridization in [Hybridization::Sp2, Hybridization::Unspecified] {
            let topology = explicit_h_star_topology(
                2,
                AtomSpec::new(Element::C).with_hybridization(hybridization),
            );
            let output = add_hydrogen_coordinates(
                explicit_only_addition(topology),
                CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(
                        0,
                        vec![[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
                    )],
                    ..CoordinateBlock::default()
                },
                true,
                false,
            )
            .expect("degree-three flat placement");
            let hydrogen = output.coordinates.conformers_2d[0].coordinates()[3];
            assert!((hydrogen[0] - flat_expected).abs() < 1.0e-12);
            assert!((hydrogen[1] - flat_expected).abs() < 1.0e-12);
        }

        let topology = explicit_h_star_topology(
            2,
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp3),
        );
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
                    true,
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("degree-three SP3 rotation");
        let hydrogen = AddHsPoint3D::from_3d(output.coordinates.conformers_3d[0].coordinates()[3]);
        let bond_length = rdkit_rb0(1) + rdkit_rb0(6);
        assert!((hydrogen.length_sq().sqrt() - bond_length).abs() < 1.0e-12);
        assert!(hydrogen.z.abs() > 1.0e-6);
        let original_bisector = AddHsPoint3D::new(flat_expected, flat_expected, 0.0);
        assert!(
            (hydrogen.scaled(1.0 / bond_length).dot(original_bisector)
                - (109.471_f64 / 2.0).to_radians().cos())
            .abs()
                < 1.0e-12
        );
    }

    #[test]
    fn h_add_coordinates_degree_three_four_degree_three_cancellation_and_coincident() {
        let topology = explicit_h_star_topology(
            2,
            AtomSpec::new(Element::C).with_hybridization(Hybridization::Sp2),
        );
        let cancelled = add_hydrogen_coordinates(
            explicit_only_addition(topology.clone()),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    0,
                    vec![[2.0, 3.0], [3.0, 3.0], [1.0, 3.0]],
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("opposed vectors retain the initialized row");
        assert_eq!(
            cancelled.coordinates.conformers_2d[0].coordinates()[3],
            [0.0, 0.0]
        );

        let coincident = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![[2.0, 3.0, 4.0], [2.0, 3.0, 4.0], [1.0, 3.0, 4.0]],
                    true,
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("coincident degree-three neighbor follows Issue 678");
        assert_eq!(
            coincident.coordinates.conformers_3d[0].coordinates()[3],
            [2.0, 3.0, 4.0]
        );
    }

    #[test]
    fn h_add_coordinates_degree_three_four_planar_fallbacks_and_chiral_volume() {
        let place = |parent_spec: AtomSpec, vectors: [[f64; 3]; 3]| {
            let topology = explicit_h_star_topology(3, parent_spec);
            let mut rows = vec![[0.0, 0.0, 0.0]];
            rows.extend(vectors.map(|vector| [-vector[0], -vector[1], -vector[2]]));
            add_hydrogen_coordinates(
                explicit_only_addition(topology),
                CoordinateBlock {
                    conformers_3d: vec![Conformer3D::new(0, rows, true)],
                    ..CoordinateBlock::default()
                },
                true,
                false,
            )
            .expect("planar degree-four placement")
            .coordinates
            .conformers_3d[0]
                .coordinates()[4]
        };

        let primary = place(
            AtomSpec::new(Element::C),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],
        );
        assert!(primary[2] > 0.0);

        let second_vectors = [[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
        let second = place(AtomSpec::new(Element::C), second_vectors);
        assert!(second[2] < 0.0);
        let chiral_parent = AtomSpec::new(Element::C)
            .with_chiral_tag(ChiralTag::TetrahedralCw)
            .with_prop("_CIPCode", "R")
            .expect("valid CIP property");
        let reversed = place(chiral_parent, second_vectors);
        assert!(reversed[2] > 0.0);

        let angle = 0.009_f64;
        let third = place(
            AtomSpec::new(Element::C),
            [
                [1.0, 0.0, 0.0],
                [angle.cos(), angle.sin(), 0.0],
                [angle.cos(), -angle.sin(), 0.0],
            ],
        );
        assert!(third[2] < 0.0);
        for point in [primary, second, reversed, third] {
            assert!(point.into_iter().all(f64::is_finite));
        }
    }

    #[test]
    fn h_add_coordinates_degree_three_four_flatland_antiparallel_bisector() {
        let topology = explicit_h_star_topology(3, AtomSpec::new(Element::C));
        let output = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    0,
                    vec![[0.0, 0.0], [0.0, -1.0], [-1.0, 0.0], [1.0, 0.0]],
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("Issues 908/3854 antiparallel bisector");
        let hydrogen = output.coordinates.conformers_2d[0].coordinates()[4];
        assert!(hydrogen[0].abs() < 1.0e-12);
        assert!((hydrogen[1] - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn h_add_coordinates_degree_three_four_nonplanar_default_and_no_nan() {
        let topology = explicit_h_star_topology(3, AtomSpec::new(Element::C));
        let nonplanar = add_hydrogen_coordinates(
            explicit_only_addition(topology.clone()),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![
                        [0.0, 0.0, 0.0],
                        [-1.0, 0.0, 0.0],
                        [0.0, -1.0, 0.0],
                        [0.0, 0.0, -1.0],
                    ],
                    true,
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("nonplanar sum");
        let point = nonplanar.coordinates.conformers_3d[0].coordinates()[4];
        assert!(
            point
                .into_iter()
                .all(|component| component > 0.0 && component.is_finite())
        );

        let all_fallbacks_fail = add_hydrogen_coordinates(
            explicit_only_addition(topology),
            CoordinateBlock {
                conformers_3d: vec![Conformer3D::new(
                    0,
                    vec![
                        [2.0, 3.0, 4.0],
                        [1.0, 3.0, 4.0],
                        [1.0, 3.0, 4.0],
                        [3.0, 3.0, 4.0],
                    ],
                    true,
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("colinear planar fallbacks retain initialized row");
        assert_eq!(
            all_fallbacks_fail.coordinates.conformers_3d[0].coordinates()[4],
            [0.0, 0.0, 0.0]
        );

        let default_topology = explicit_h_star_topology(4, AtomSpec::new(Element::C));
        let default = add_hydrogen_coordinates(
            explicit_only_addition(default_topology),
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(
                    0,
                    vec![[5.0, 6.0], [6.0, 6.0], [5.0, 7.0], [4.0, 6.0], [5.0, 5.0]],
                )],
                ..CoordinateBlock::default()
            },
            true,
            false,
        )
        .expect("source default active degree");
        assert_eq!(
            default.coordinates.conformers_2d[0].coordinates()[5],
            [0.0, 0.0]
        );
    }
}

// BEGIN ADDHS COMPUTED PRELUDE FROZEN REGRESSIONS
#[cfg(test)]
mod add_hs_computed_prelude_tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomQueryPredicate, BondQueryPredicate, PropertyValue, QueryNode,
        ordered_atom_properties, ordered_bond_properties,
    };
    use std::collections::BTreeSet;

    fn check(errors: &mut Vec<String>, label: &str, field: &str, valid: bool) {
        if !valid {
            errors.push(format!("{label}: {field}"));
        }
    }

    #[derive(Debug, Clone, PartialEq)]
    struct Input {
        topology: TopologyBlock,
        coordinates: CoordinateBlock,
        properties: MoleculeProperties,
        params: AddHsParams,
        query_atoms: Vec<QueryAtom>,
        query_bonds: Vec<QueryBond>,
    }

    fn selection(choice: usize) -> Option<Vec<AtomId>> {
        match choice {
            0 => None,
            1 => Some(vec![]),
            2 => Some(vec![AtomId::new(0)]),
            3 => Some(vec![AtomId::new(2)]),
            _ => panic!("only the frozen four selection literals exist"),
        }
    }

    fn input(config: usize, hydrogen: u8, choice: usize, query: bool, bad_coords: bool) -> Input {
        let mut atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_explicit_hydrogens(hydrogen),
            ),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        for (index, atom) in atoms.iter_mut().enumerate() {
            atom.set_prop("label", ["atom0", "atom1"][index]).unwrap();
            atom.set_prop("_CIPCode", ["R", "S"][index]).unwrap();
            if config & 1 != 0 {
                atom.set_computed_prop("_CIPRank", PropertyValue::Int([17, 19][index]))
                    .unwrap();
                atom.set_computed_prop("scratch", PropertyValue::Int([17, 19][index]))
                    .unwrap();
            }
        }
        let mut bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        bond.set_prop("label", "bond0").unwrap();
        if config & 2 != 0 {
            bond.set_computed_prop("scratch", PropertyValue::Int(23))
                .unwrap();
        }
        let bonds = vec![bond];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(2, &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let mut properties = MoleculeProperties::default();
        properties.set_prop("label", "mol").unwrap();
        properties.set_computed_prop("memo", "discard").unwrap();
        let coordinates = if bad_coords {
            CoordinateBlock {
                conformers_2d: vec![Conformer2D::new(17, vec![[0.0, -0.0]])],
                ..CoordinateBlock::default()
            }
        } else {
            CoordinateBlock::default()
        };
        let (query_atoms, query_bonds) = if query {
            (
                vec![
                    QueryAtom::from_parts(
                        topology.atoms[0].clone(),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    ),
                    QueryAtom::from_carrier_parts(
                        topology.atoms[1].clone(),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    ),
                ],
                vec![QueryBond::from_carrier_parts(
                    topology.bonds[0].clone(),
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                )],
            )
        } else {
            (vec![], vec![])
        };
        Input {
            topology,
            coordinates,
            properties,
            params: AddHsParams {
                explicit_only: true,
                add_coords: false,
                add_residue_info: false,
                skip_queries: query,
                only_on_atoms: selection(choice),
            },
            query_atoms,
            query_bonds,
        }
    }

    // Full typed input values are compared separately from Debug/order and bits.
    fn coordinate_bits(
        coords: &CoordinateBlock,
    ) -> (Vec<(usize, Vec<[u64; 2]>)>, Vec<(usize, Vec<[u64; 3]>)>) {
        (
            coords
                .conformers_2d
                .iter()
                .map(|c| {
                    (
                        c.id(),
                        c.coordinates()
                            .iter()
                            .map(|p| p.map(f64::to_bits))
                            .collect(),
                    )
                })
                .collect(),
            coords
                .conformers_3d
                .iter()
                .map(|c| {
                    (
                        c.id(),
                        c.coordinates()
                            .iter()
                            .map(|p| p.map(f64::to_bits))
                            .collect(),
                    )
                })
                .collect(),
        )
    }

    fn prerequisites(
        value: &Input,
        config: usize,
        hydrogen: u8,
        choice: usize,
        query: bool,
        bad_coords: bool,
        errors: &mut Vec<String>,
        label: &str,
    ) {
        check(
            errors,
            label,
            "pre topology validation",
            value.topology.validate().is_ok(),
        );
        check(
            errors,
            label,
            "pre literal counts",
            value.topology.atoms.len() == 2 && value.topology.bonds.len() == 1,
        );
        for (index, atom) in value.topology.atoms.iter().enumerate() {
            check(
                errors,
                label,
                "pre atom identity/H",
                atom.id() == AtomId::new(index)
                    && atom.element() == Element::C
                    && atom.explicit_hydrogens() == if index == 0 { hydrogen } else { 0 },
            );
            check(
                errors,
                label,
                "pre typed ordinary atom properties",
                atom.prop("label")
                    == Some(&PropertyValue::String(["atom0", "atom1"][index].into()))
                    && atom.prop("_CIPCode")
                        == Some(&PropertyValue::String(["R", "S"][index].into()))
                    && !atom.is_prop_computed("label").unwrap()
                    && !atom.is_prop_computed("_CIPCode").unwrap(),
            );
            let keys = ordered_atom_properties(atom)
                .map(|(key, _)| super::fixed_property_text(key))
                .collect::<Vec<_>>();
            let wanted = if config & 1 != 0 {
                vec![
                    "label",
                    "_CIPCode",
                    "__computedProps",
                    "_CIPRank",
                    "scratch",
                ]
            } else {
                vec!["label", "_CIPCode"]
            };
            check(errors, label, "pre atom insertion order", keys == wanted);
            let names: BTreeSet<String> = if config & 1 != 0 {
                ["_CIPRank", "scratch"]
                    .into_iter()
                    .map(String::from)
                    .collect()
            } else {
                BTreeSet::new()
            };
            check(
                errors,
                label,
                "pre atom computed membership",
                atom.computed_prop_names().unwrap().map(|v| {
                    v.iter()
                        .map(|k| super::fixed_property_text(k).to_owned())
                        .collect::<BTreeSet<_>>()
                }) == (config & 1 != 0).then_some(names)
                    && atom.computed_prop_names().unwrap()
                        == (config & 1 != 0).then_some(&["_CIPRank".into(), "scratch".into()][..]),
            );
            for key in ["_CIPRank", "scratch"] {
                let wanted = (config & 1 != 0).then_some(PropertyValue::Int([17, 19][index]));
                check(
                    errors,
                    label,
                    "pre typed computed atom value",
                    atom.prop(key) == wanted.as_ref(),
                );
            }
        }
        let bond = &value.topology.bonds[0];
        check(
            errors,
            label,
            "pre literal bond",
            bond.id() == BondId::new(0)
                && bond.begin() == AtomId::new(0)
                && bond.end() == AtomId::new(1)
                && bond.order() == BondOrder::Single,
        );
        check(
            errors,
            label,
            "pre typed ordinary bond value",
            bond.prop("label") == Some(&PropertyValue::String("bond0".into()))
                && !bond.is_prop_computed("label").unwrap(),
        );
        check(
            errors,
            label,
            "pre typed computed bond value",
            bond.prop("scratch") == (config & 2 != 0).then_some(PropertyValue::Int(23)).as_ref(),
        );
        let names: BTreeSet<String> = if config & 2 != 0 {
            ["scratch"].into_iter().map(String::from).collect()
        } else {
            BTreeSet::new()
        };
        check(
            errors,
            label,
            "pre bond computed membership",
            bond.computed_prop_names().unwrap().map(|v| {
                v.iter()
                    .map(|k| super::fixed_property_text(k).to_owned())
                    .collect::<BTreeSet<_>>()
            }) == (config & 2 != 0).then_some(names),
        );
        let order = if config & 2 != 0 {
            vec!["label", "__computedProps", "scratch"]
        } else {
            vec!["label"]
        };
        check(
            errors,
            label,
            "pre bond insertion order",
            ordered_bond_properties(bond)
                .map(|(key, _)| super::fixed_property_text(key))
                .collect::<Vec<_>>()
                == order,
        );
        check(
            errors,
            label,
            "pre molecule typed string values/markers/key order",
            value.properties.prop("label")
                == Some(&cosmolkit_model::PropertyValue::String("mol".into()))
                && value.properties.prop("memo")
                    == Some(&cosmolkit_model::PropertyValue::String("discard".into()))
                && !value.properties.is_prop_computed("label").unwrap()
                && value.properties.is_prop_computed("memo").unwrap()
                && value.properties.computed_prop_names().unwrap() == Some(&["memo".into()][..])
                && value
                    .properties
                    .ordered_props()
                    .map(|(key, _)| super::fixed_property_text(key))
                    .collect::<Vec<_>>()
                    == vec!["label", "__computedProps", "memo"],
        );
        check(
            errors,
            label,
            "pre exact params and selection",
            value.params
                == AddHsParams {
                    explicit_only: true,
                    add_coords: false,
                    add_residue_info: false,
                    skip_queries: query,
                    only_on_atoms: selection(choice),
                },
        );
        if bad_coords {
            check(
                errors,
                label,
                "pre exact invalid coordinates",
                value.coordinates.conformers_2d.len() == 1
                    && value.coordinates.conformers_2d[0].id() == 17
                    && coordinate_bits(&value.coordinates).0
                        == vec![(17, vec![[0.0f64.to_bits(), (-0.0f64).to_bits()]])]
                    && value.coordinates.conformers_3d.is_empty()
                    && value.coordinates.validate_for_atom_count(2)
                        == Err(CoordinateValidationError::RowCount {
                            dimension: "2D",
                            conformer: 17,
                            rows: 1,
                            atom_count: 2,
                        }),
            );
        } else {
            check(
                errors,
                label,
                "pre empty coordinates",
                value.coordinates == CoordinateBlock::default(),
            );
        }
        if query {
            let state = QueryStateRef::try_for_topology(
                &value.query_atoms,
                &value.query_bonds,
                &value.topology,
            );
            check(errors, label, "pre query validation", state.is_ok());
            if let Ok(state) = state {
                check(
                    errors,
                    label,
                    "pre actual/carrier query identities",
                    state.atom_has_query(AtomId::new(0))
                        && !state.atom_has_query(AtomId::new(1))
                        && !state.bond_has_query(BondId::new(0)),
                );
                check(
                    errors,
                    label,
                    "pre query predicate literals",
                    state.atom_predicate(AtomId::new(0))
                        == &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))
                        && state.atom_predicate(AtomId::new(1))
                            == &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))
                        && state.bond_predicate(BondId::new(0))
                            == &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                );
            }
        } else {
            check(
                errors,
                label,
                "pre no query rows",
                value.query_atoms.is_empty() && value.query_bonds.is_empty(),
            );
        }
    }

    enum Output {
        Topology(AddHydrogensTopologyResult),
        Full(AddHydrogensResult),
    }

    // One actual entry call. Return inspection is deliberately AFTER preservation.
    fn call(
        value: &Input,
        full: bool,
        census: &mut usize,
        errors: &mut Vec<String>,
        label: &str,
    ) -> Result<Output, HydrogenError> {
        let state = if value.params.skip_queries {
            Some(
                QueryStateRef::try_for_topology(
                    &value.query_atoms,
                    &value.query_bonds,
                    &value.topology,
                )
                .unwrap(),
            )
        } else {
            None
        };
        let baseline = value.clone();
        let baseline_debug = format!("{value:?}");
        let baseline_bits = coordinate_bits(&value.coordinates);
        let result = if full {
            if state.is_some() {
                add_hydrogens_with_query_state(
                    value.topology.clone(),
                    value.coordinates.clone(),
                    value.properties.clone(),
                    &value.params,
                    state,
                )
            } else {
                add_hydrogens_with_params(
                    value.topology.clone(),
                    value.coordinates.clone(),
                    value.properties.clone(),
                    &value.params,
                )
            }
            .map(Output::Full)
        } else {
            if state.is_some() {
                add_hydrogens_topology_with_query_state(
                    value.topology.clone(),
                    &value.params,
                    state,
                )
            } else {
                add_hydrogens_topology(value.topology.clone(), &value.params)
            }
            .map(Output::Topology)
        };
        *census += 1;
        let preserved = value == &baseline
            && format!("{value:?}") == baseline_debug
            && coordinate_bits(&value.coordinates) == baseline_bits;
        check(
            errors,
            label,
            "post retained whole input/order/query/coordinate bits",
            preserved,
        );
        println!("ADDHS_PRELUDE_CALL {label} preservation={preserved}");
        result
    }

    // Expected state starts from the PRE-call clone and literal source removals.
    // No production computed-clear/AddHs helper derives this expectation.
    fn expected_topology(before: &TopologyBlock, append: bool, final_h: u8) -> TopologyBlock {
        let mut expected = before.clone();
        for atom in &mut expected.atoms {
            atom.clear_prop("_CIPRank");
            atom.clear_prop("scratch");
        }
        expected.bonds[0].clear_prop("scratch");
        expected.atoms[0].set_explicit_hydrogens(final_h);
        if append {
            expected
                .atoms
                .push(Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::H)));
            expected.bonds.push(Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ));
            expected.adjacency = AdjacencyList::from_topology(3, &expected.bonds);
        }
        expected
    }

    fn success(
        result: Result<Output, HydrogenError>,
        before: &Input,
        append: bool,
        final_h: u8,
        errors: &mut Vec<String>,
        label: &str,
    ) {
        let output = match result {
            Ok(value) => value,
            Err(error) => {
                errors.push(format!("{label}: unexpected {error:?}"));
                return;
            }
        };
        let (topology, mapping, additions) = match &output {
            Output::Topology(value) => (&value.topology, &value.mapping, Some(&value.additions)),
            Output::Full(value) => {
                let mut properties = before.properties.clone();
                properties.clear_prop("memo");
                check(
                    errors,
                    label,
                    "post whole molecule properties",
                    value.properties == properties,
                );
                check(
                    errors,
                    label,
                    "post empty coordinate output",
                    value.coordinates == CoordinateBlock::default()
                        && coordinate_bits(&value.coordinates) == (vec![], vec![]),
                );
                check(errors, label, "post no warnings", value.warnings.is_empty());
                (&value.topology, &value.mapping, None)
            }
        };
        let expected = expected_topology(&before.topology, append, final_h);
        check(
            errors,
            label,
            "post whole typed topology and insertion order",
            topology == &expected && format!("{topology:?}") == format!("{expected:?}"),
        );
        check(
            errors,
            label,
            "post output model validation",
            topology.validate().is_ok(),
        );
        let atom_old = vec![Some(AtomId::new(0)), Some(AtomId::new(1))];
        let bond_old = vec![Some(BondId::new(0))];
        let mut atom_new = atom_old.clone();
        let mut bond_new = bond_old.clone();
        if append {
            atom_new.push(None);
            bond_new.push(None);
        }
        check(
            errors,
            label,
            "post exact literal mapping",
            mapping.atoms().old_to_new() == atom_old
                && mapping.atoms().new_to_old() == atom_new
                && mapping.bonds().old_to_new() == bond_old
                && mapping.bonds().new_to_old() == bond_new,
        );
        check(
            errors,
            label,
            "post mapping validation",
            mapping
                .validate_for_counts(2, if append { 3 } else { 2 }, 1, if append { 2 } else { 1 })
                .is_ok(),
        );
        if let Some(additions) = additions {
            let wanted = if append {
                vec![AddedHydrogen {
                    atom: AtomId::new(2),
                    bond: BondId::new(1),
                    parent: AtomId::new(0),
                    kind: AddedHydrogenKind::Explicit,
                }]
            } else {
                vec![]
            };
            check(
                errors,
                label,
                "post literal addition metadata",
                additions == &wanted,
            );
        }
        // Safely collect local diagnoses even when the whole returned graph differs.
        for (index, atom) in topology.atoms.iter().take(2).enumerate() {
            check(
                errors,
                label,
                "post no computed atom keys/markers",
                atom.prop("_CIPRank").is_none()
                    && atom.prop("scratch").is_none()
                    && atom.computed_prop_names().unwrap()
                        == before.topology.atoms[index]
                            .computed_prop_names()
                            .unwrap()
                            .map(|_| &[] as &[cosmolkit_model::PropertyText]),
            );
            check(
                errors,
                label,
                "post ordinary atom values and insertion order",
                atom.prop("label")
                    == Some(&PropertyValue::String(["atom0", "atom1"][index].into()))
                    && atom.prop("_CIPCode")
                        == Some(&PropertyValue::String(["R", "S"][index].into()))
                    && ordered_atom_properties(atom)
                        .map(|(key, _)| super::fixed_property_text(key))
                        .collect::<Vec<_>>()
                        == if before.topology.atoms[index]
                            .computed_prop_names()
                            .unwrap()
                            .is_some()
                        {
                            vec!["label", "_CIPCode", "__computedProps"]
                        } else {
                            vec!["label", "_CIPCode"]
                        },
            );
        }
        if let Some(bond) = topology.bonds.first() {
            check(
                errors,
                label,
                "post bond computed absence/ordinary/order",
                bond.prop("scratch").is_none()
                    && bond.computed_prop_names().unwrap()
                        == before.topology.bonds[0]
                            .computed_prop_names()
                            .unwrap()
                            .map(|_| &[] as &[cosmolkit_model::PropertyText])
                    && bond.prop("label") == Some(&PropertyValue::String("bond0".into()))
                    && ordered_bond_properties(bond)
                        .map(|(key, _)| super::fixed_property_text(key))
                        .collect::<Vec<_>>()
                        == if before.topology.bonds[0]
                            .computed_prop_names()
                            .unwrap()
                            .is_some()
                        {
                            vec!["label", "__computedProps"]
                        } else {
                            vec!["label"]
                        },
            );
        }
    }

    fn finish(census: usize, expected: usize, errors: Vec<String>, name: &str) {
        let mut errors = errors;
        check(&mut errors, name, "exact actual census", census == expected);
        println!(
            "ADDHS_PRELUDE_CENSUS {name} actual={census} expected={expected} discrepancies={}",
            errors.len()
        );
        assert!(errors.is_empty(), "{name}: {errors:#?}");
    }

    #[test]
    fn add_hs_computed_prelude_product24() {
        let mut census = 0;
        let mut errors = vec![];
        for choice in [0, 1, 2] {
            for config in [0, 1, 2, 3] {
                for full in [false, true] {
                    let label = format!("product24/selection={choice}/config={config}/full={full}");
                    let value = input(config, 0, choice, false, false);
                    prerequisites(&value, config, 0, choice, false, false, &mut errors, &label);
                    let before = value.clone();
                    let result = call(&value, full, &mut census, &mut errors, &label);
                    success(result, &before, false, 0, &mut errors, &label);
                }
            }
        }
        finish(census, 24, errors, "product24");
    }

    #[test]
    fn add_hs_computed_prelude_append6() {
        let mut census = 0;
        let mut errors = vec![];
        for choice in [0, 1, 2] {
            for full in [false, true] {
                let label = format!("append6/selection={choice}/full={full}");
                let value = input(3, 1, choice, false, false);
                prerequisites(&value, 3, 1, choice, false, false, &mut errors, &label);
                let before = value.clone();
                let result = call(&value, full, &mut census, &mut errors, &label);
                success(
                    result,
                    &before,
                    choice != 1,
                    u8::from(choice == 1),
                    &mut errors,
                    &label,
                );
            }
        }
        finish(census, 6, errors, "append6");
    }

    #[test]
    fn add_hs_computed_prelude_query4() {
        let mut census = 0;
        let mut errors = vec![];
        for choice in [0, 2] {
            for full in [false, true] {
                let label = format!("query4/selection={choice}/full={full}");
                let value = input(3, 1, choice, true, false);
                prerequisites(&value, 3, 1, choice, true, false, &mut errors, &label);
                let before = value.clone();
                let result = call(&value, full, &mut census, &mut errors, &label);
                success(result, &before, false, 1, &mut errors, &label);
            }
        }
        finish(census, 4, errors, "query4");
    }

    #[test]
    fn add_hs_computed_prelude_error4() {
        use std::error::Error;
        let mut census = 0;
        let mut errors = vec![];
        for repeat in [0, 1] {
            for bad_coords in [true, false] {
                let choice = if bad_coords { 0 } else { 3 };
                let label = format!("error4/repeat={repeat}/bad_coords={bad_coords}");
                let value = input(3, 0, choice, false, bad_coords);
                prerequisites(&value, 3, 0, choice, false, bad_coords, &mut errors, &label);
                let result = call(&value, true, &mut census, &mut errors, &label);
                match result {
                    Err(error) => {
                        check(
                            &mut errors,
                            &label,
                            "existing source is None",
                            error.source().is_none(),
                        );
                        if bad_coords {
                            match &error {
                                HydrogenError::InvalidCoordinates(payload) => {
                                    let expected = CoordinateValidationError::RowCount {
                                        dimension: "2D",
                                        conformer: 17,
                                        rows: 1,
                                        atom_count: 2,
                                    };
                                    check(
                                        &mut errors,
                                        &label,
                                        "exact borrowed coordinate payload",
                                        payload == &expected,
                                    );
                                    // Independent cast proves payload identity ONLY, no cause chain.
                                    let erased: &(dyn Error + 'static) = payload;
                                    check(
                                        &mut errors,
                                        &label,
                                        "independent payload downcast identity",
                                        erased
                                            .downcast_ref::<CoordinateValidationError>()
                                            .is_some_and(|same| std::ptr::eq(same, payload)),
                                    );
                                }
                                other => errors.push(format!("{label}: wrong error {other:?}")),
                            }
                        } else {
                            check(
                                &mut errors,
                                &label,
                                "exact selection error",
                                error
                                    == HydrogenError::OnlyOnAtomOutOfRange {
                                        atom: AtomId::new(2),
                                        atom_count: 2,
                                    },
                            );
                        }
                    }
                    Ok(_) => errors.push(format!("{label}: unexpected success")),
                }
            }
        }
        finish(census, 4, errors, "error4");
    }
}
// END ADDHS COMPUTED PRELUDE FROZEN REGRESSIONS

#[doc(hidden)]
pub fn add_hydrogens_with_source_valence(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &AddHsParams,
    source_valence: Option<ValenceAssignment>,
) -> Result<AddHydrogensResult, HydrogenError> {
    add_hydrogens_from_source_state(
        topology,
        coordinates,
        properties,
        params,
        None,
        source_valence,
    )
}

fn add_hydrogens_from_source_state(
    topology: TopologyBlock,
    coordinates: CoordinateBlock,
    properties: MoleculeProperties,
    params: &AddHsParams,
    query_state: Option<QueryStateRef<'_>>,
    source_valence: Option<ValenceAssignment>,
) -> Result<AddHydrogensResult, HydrogenError> {
    validate_blocks(&topology, &coordinates)?;
    validate_property_lists(&properties, topology.atoms.len(), topology.bonds.len())?;
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }

    let old_atom_count = topology.atoms.len();
    let old_bond_count = topology.bonds.len();
    let processed_atoms = processed_add_hydrogen_atoms(&topology, params, query_state)?;
    let tracked_isotopes = processed_atoms
        .iter()
        .enumerate()
        .filter(|(_, processed)| **processed)
        .map(|(atom, _)| {
            (
                AtomId::new(atom),
                topology.atoms[atom].tracked_isotopic_hydrogens().to_vec(),
            )
        })
        .collect::<Vec<_>>();

    // BEGIN RDKIT CPP FUNCTION MolOps::addHs computed-property prelude
    // RDKit✔️✔️: void addHs(RWMol &mol, const AddHsParameters &params,
    // RDKit✔️✔️:            const UINT_VECT *onlyOnAtoms) {
    // RDKit✔️✔️:   // when we hit each atom, clear its computed properties
    // RDKit✔️✔️:   // NOTE: it is essential that we not clear the ring info in the
    // RDKit✔️✔️:   // molecule's computed properties.  We don't want to have to
    // RDKit✔️✔️:   // regenerate that.  This caused Issue210 and Issue212:
    // RDKit✔️✔️:   mol.clearComputedProps(false);
    // END RDKIT CPP FUNCTION MolOps::addHs computed-property prelude
    // Behavior review: the composed boundary clears molecule properties HERE
    // once and ALL atom/bond computed properties in the delegated topology
    // owner. Ordinary properties and input scalar fields survive the clear.
    // The source scalar carrier is updated only by its selected-atom stages.
    // Cost review: molecule clearing moves the computed set and removes keys
    // in O(q log p). The topology owner's separate store-cost qualification
    // applies to the composed prelude; no whole-topology clone is added.
    let mut properties = properties;
    properties
        .clear_computed_props()
        .map_err(HydrogenError::InvalidProperty)?;

    let (result, final_valence) =
        add_hydrogens_topology_from_source_state(topology, params, query_state, source_valence)?;
    let mut result = add_hydrogen_coordinates(
        result,
        coordinates,
        params.add_coords,
        params.add_residue_info,
    )?;
    let warnings =
        replay_tracked_isotopes(&mut result.topology, &result.additions, &tracked_isotopes)?;

    result.mapping.validate_for_counts(
        old_atom_count,
        result.topology.atoms.len(),
        old_bond_count,
        result.topology.bonds.len(),
    )?;
    properties.remap_topology(
        result.mapping.atoms().new_to_old(),
        result.mapping.bonds().new_to_old(),
    );
    validate_blocks(&result.topology, &result.coordinates)?;
    validate_property_lists(
        &properties,
        result.topology.atoms.len(),
        result.topology.bonds.len(),
    )?;
    Ok(AddHydrogensResult {
        topology: result.topology,
        coordinates: result.coordinates,
        properties,
        mapping: result.mapping,
        warnings,
        final_valence,
    })
}

fn add_hydrogens_topology_from_source_state(
    mut topology: TopologyBlock,
    params: &AddHsParams,
    query_state: Option<QueryStateRef<'_>>,
    source_valence: Option<ValenceAssignment>,
) -> Result<(AddHydrogensTopologyResult, ValenceAssignment), HydrogenError> {
    topology.validate()?;
    if let Some(state) = query_state {
        state.validate_for_topology(&topology)?;
    }
    let old_atom_count = topology.atoms.len();
    let old_bond_count = topology.bonds.len();
    let mut selected = selected_atoms(&topology, params.only_on_atoms.as_deref())?;

    // BEGIN RDKIT CPP FUNCTION MolOps::addHs global topology prelude
    // RDKit✔️❌: void addHs(RWMol &mol, const AddHsParameters &params,
    // RDKit✔️❌:            const UINT_VECT *onlyOnAtoms) {
    // RDKit✔️❌:   // when we hit each atom, clear its computed properties
    // RDKit✔️❌:   // NOTE: it is essential that we not clear the ring info in the
    // RDKit✔️❌:   // molecule's computed properties.  We don't want to have to
    // RDKit✔️❌:   // regenerate that.  This caused Issue210 and Issue212:
    // RDKit✔️❌:   mol.clearComputedProps(false);
    // END RDKIT CPP FUNCTION MolOps::addHs global topology prelude
    // Behavior review: every atom and bond is cleared, even empty selection,
    // skipped queries and no append. The source's later selected-atom clear
    // remains below. CK retains its prior typed topology/query/selection
    // validation before this owned edit; malformed-input native parity is
    // not claimed. Molecule clearing belongs to the full-block wrapper.
    // Cost review: one shared traversal with no clone, allocation or refind
    // here; inherited PropertyStore retention cost is qualified in the helper.
    clear_hydrogen_topology_computed_properties(&mut topology)?;

    // BEGIN RDKIT CPP FUNCTION MolOps::addHs selection/count snapshot
    // RDKit✔️❌: unsigned int numAddHyds = 0;
    // RDKit✔️❌: boost::dynamic_bitset<> onAtoms(mol.getNumAtoms());
    // RDKit✔️❌: if (onlyOnAtoms) {
    // RDKit✔️❌:   for (auto atIdx : *onlyOnAtoms) {
    // RDKit✔️❌:     onAtoms.set(atIdx);
    // RDKit✔️❌:   }
    // RDKit✔️❌: } else {
    // RDKit✔️❌:   onAtoms.set();
    // RDKit✔️❌: }
    // RDKit✔️❌: std::vector<unsigned int> numExplicitHs(mol.getNumAtoms(), 0);
    // RDKit✔️❌: std::vector<unsigned int> numImplicitHs(mol.getNumAtoms(), 0);
    // RDKit✔️❌: for (auto at : mol.atoms()) {
    // RDKit✔️❌:   numExplicitHs[at->getIdx()] = at->getNumExplicitHs();
    // RDKit✔️❌:   numImplicitHs[at->getIdx()] = at->getNumImplicitHs();
    // RDKit✔️❌:   if (onAtoms[at->getIdx()]) {
    // RDKit✔️❌:     if (params.skipQueries && isQueryAtom(mol, *at)) {
    // RDKit✔️❌:       onAtoms.set(at->getIdx(), 0);
    // RDKit✔️❌:       continue;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     numAddHyds += at->getNumExplicitHs();
    // RDKit✔️❌:     if (!params.explicitOnly) {
    // RDKit✔️❌:       numAddHyds += at->getNumImplicitHs();
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::addHs selection/count snapshot
    // The detached port preserves the source O(atoms + selected hydrogens)
    // traversal but rebuilds adjacency after appending, which allocates an
    // additional O(atoms + bonds) structure compared with RDKit's incremental
    // graph update.
    // RDKit✔️✔️:   d_implicitValence = -1;
    // RDKit✔️✔️:   d_explicitValence = -1;
    // Atom::initAtom establishes a real uninitialized scalar state, not zero.
    // Missing source state is represented literally; the unique getter below
    // preserves noImplicit's short circuit and otherwise raises its typed error.
    let mut valence = source_valence.unwrap_or_else(|| ValenceAssignment {
        explicit_valence: vec![-1; old_atom_count],
        implicit_hydrogens: vec![-1; old_atom_count],
    });
    validate_removal_valence(&topology, &valence)?;
    let explicit_counts = topology
        .atoms
        .iter()
        .map(|atom| usize::from(atom.explicit_hydrogens()))
        .collect::<Vec<_>>();
    // Source snapshots ALL implicit scalars before selection/query skipping,
    // even in explicitOnly mode. No fresh graph assignment or max(0) fallback.
    let implicit_counts = topology
        .atoms
        .iter()
        .map(|atom| {
            crate::hcount::implicit_hydrogen_count(atom, &valence).map(|count| count as usize)
        })
        .collect::<Result<Vec<_>, _>>()?;
    if params.skip_queries {
        for atom_index in 0..old_atom_count {
            if selected[atom_index]
                && is_query_atom(&topology, AtomId::new(atom_index), query_state)
            {
                selected[atom_index] = false;
            }
        }
    }
    let addition_count = (0..old_atom_count)
        .filter(|&atom_index| selected[atom_index])
        .map(|atom_index| {
            explicit_counts[atom_index]
                + if params.explicit_only {
                    0
                } else {
                    implicit_counts[atom_index]
                }
        })
        .sum::<usize>();
    topology.atoms.reserve(addition_count);
    topology.bonds.reserve(addition_count);
    let mut additions = Vec::with_capacity(addition_count);

    // BEGIN RDKIT CPP FUNCTION MolOps::addHs topology append loop
    // RDKit✔️❌: unsigned int stopIdx = mol.getNumAtoms();
    // RDKit✔️❌: for (unsigned int aidx = 0; aidx < stopIdx; ++aidx) {
    // RDKit✔️❌:   if (!onAtoms[aidx]) {
    // RDKit✔️❌:     continue;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   Atom *newAt = mol.getAtomWithIdx(aidx);
    // RDKit✔️❌:   newAt->clearComputedProps();
    // RDKit✔️❌:   // always convert explicit Hs
    // RDKit✔️❌:   unsigned int onumexpl = numExplicitHs[aidx];
    // RDKit✔️❌:   for (unsigned int i = 0; i < onumexpl; i++) {
    // RDKit✔️❌:     newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️❌:     mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️❌:     auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️❌:     hAtom->updatePropertyCache();
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // clear the local property
    // RDKit✔️❌:   newAt->setNumExplicitHs(0);
    // RDKit✔️❌:   if (!params.explicitOnly) {
    // RDKit✔️❌:     // take care of implicits
    // RDKit✔️❌:     for (unsigned int i = 0; i < numImplicitHs[aidx]; i++) {
    // RDKit✔️❌:       newIdx = mol.addAtom(new Atom(1), false, true);
    // RDKit✔️❌:       mol.addBond(aidx, newIdx, Bond::SINGLE);
    // RDKit✔️❌:       // set the isImplicit label so that we can strip these back
    // RDKit✔️❌:       // off later if need be.
    // RDKit✔️❌:       auto hAtom = mol.getAtomWithIdx(newIdx);
    // RDKit✔️❌:       hAtom->setProp(common_properties::isImplicit, 1);
    // RDKit✔️❌:       hAtom->updatePropertyCache();
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // update the atom's derived properties (valence count, etc.)
    // RDKit✔️❌:   // no sense in being strict here (was github #2782)
    // RDKit✔️❌:   newAt->updatePropertyCache(false);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION MolOps::addHs topology append loop
    for atom_index in 0..old_atom_count {
        if !selected[atom_index] {
            continue;
        }
        topology.atoms[atom_index].clear_computed_props()?;
        let parent = AtomId::new(atom_index);
        for _ in 0..explicit_counts[atom_index] {
            append_hydrogen(
                &mut topology,
                &mut additions,
                parent,
                AddedHydrogenKind::Explicit,
            );
        }
        topology.atoms[atom_index].set_explicit_hydrogens(0);
        if !params.explicit_only {
            for _ in 0..implicit_counts[atom_index] {
                append_hydrogen(
                    &mut topology,
                    &mut additions,
                    parent,
                    AddedHydrogenKind::Implicit,
                );
            }
        }
    }
    topology.adjacency = AdjacencyList::try_from_topology(topology.atoms.len(), &topology.bonds)
        .map_err(|_| HydrogenError::InvalidTopology(TopologyValidationError::AdjacencyMismatch))?;
    topology.validate()?;
    // RDKit✔️❌:       hAtom->updatePropertyCache();
    // RDKit✔️❌:     newAt->updatePropertyCache(false);
    // Reuse the unique per-atom valence owner. Other parents' appended Hs do
    // not affect these incident bonds, so calculating after the one detached
    // adjacency rebuild preserves selected-parent/new-H source scalar values.
    // Traversal remains O(V+E+H); the detached adjacency allocation and input
    // carrier materialization cost more than source incremental graph edits.
    // Unselected rows are NEVER refreshed, including source -1 sentinels.
    valence.explicit_valence.resize(topology.atoms.len(), -1);
    valence.implicit_hydrogens.resize(topology.atoms.len(), -1);
    let mut addition_cursor = 0;
    for (parent, process) in selected.iter().copied().enumerate() {
        if !process {
            continue;
        }
        while addition_cursor < additions.len()
            && additions[addition_cursor].parent.index() == parent
        {
            let atom = additions[addition_cursor].atom;
            let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
                &topology.atoms,
                &topology.bonds,
                &topology.adjacency,
                atom,
                true,
            )?;
            valence.explicit_valence[atom.index()] = explicit;
            valence.implicit_hydrogens[atom.index()] = implicit;
            addition_cursor += 1;
        }
        let (explicit, implicit) = assign_valence_state_for_atom_from_parts(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            AtomId::new(parent),
            false,
        )?;
        valence.explicit_valence[parent] = explicit;
        valence.implicit_hydrogens[parent] = implicit;
    }
    validate_removal_valence(&topology, &valence)?;
    let mapping = TopologyMapping::with_appended(
        old_atom_count,
        old_bond_count,
        additions.len(),
        additions.len(),
    );
    mapping.validate_for_counts(
        old_atom_count,
        topology.atoms.len(),
        old_bond_count,
        topology.bonds.len(),
    )?;
    Ok((
        AddHydrogensTopologyResult {
            topology,
            mapping,
            additions,
        },
        valence,
    ))
}

#[cfg(test)]
mod source_add_hydrogen_scalar_transport {
    use super::*;
    use cosmolkit_model::{Atom, Element};
    fn topology(no_implicit: bool) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_no_implicit(no_implicit),
                ),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::C).with_no_implicit(no_implicit),
                ),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn all_implicit_getters_precede_explicit_only_and_empty_selection() {
        for explicit_only in [false, true] {
            let t = topology(false);
            let before = t.clone();
            let p = AddHsParams {
                explicit_only,
                only_on_atoms: Some(vec![]),
                ..Default::default()
            };
            for (expected, state) in [
                (0, None),
                (
                    1,
                    Some(ValenceAssignment {
                        explicit_valence: vec![0; 2],
                        implicit_hydrogens: vec![0, -1],
                    }),
                ),
            ] {
                let err = add_hydrogens_with_source_valence(
                    t.clone(),
                    CoordinateBlock::default(),
                    MoleculeProperties::default(),
                    &p,
                    state,
                )
                .unwrap_err();
                assert!(
                    matches!(err,HydrogenError::Valence(ValenceError::ImplicitValenceCacheNotInitialized { atom }) if atom.index()==expected),
                    "{err:?}"
                );
            }
            assert_eq!(t, before);
        }
    }
    #[test]
    fn no_implicit_bypasses_missing_fields_and_unselected_rows_survive() {
        let t = topology(true);
        let p = AddHsParams {
            only_on_atoms: Some(vec![]),
            ..Default::default()
        };
        let result = add_hydrogens_with_source_valence(
            t.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &p,
            None,
        )
        .unwrap();
        assert_eq!(result.topology, t);
        assert_eq!(
            result.final_valence,
            ValenceAssignment {
                explicit_valence: vec![-1; 2],
                implicit_hydrogens: vec![-1; 2]
            }
        );
        let p = AddHsParams {
            only_on_atoms: Some(vec![AtomId::new(0)]),
            ..Default::default()
        };
        let state = ValenceAssignment {
            explicit_valence: vec![99, 77],
            implicit_hydrogens: vec![0, -1],
        };
        let result = add_hydrogens_with_source_valence(
            t.clone(),
            CoordinateBlock::default(),
            MoleculeProperties::default(),
            &p,
            Some(state),
        )
        .unwrap();
        assert_eq!(result.topology, t);
        assert_eq!(
            result.final_valence,
            ValenceAssignment {
                explicit_valence: vec![0, 77],
                implicit_hydrogens: vec![0, -1]
            }
        );
    }
}
