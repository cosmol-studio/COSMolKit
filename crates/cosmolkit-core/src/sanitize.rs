//! RDKit-aligned detached sanitization orchestration and assignments.
//!
//! This module sequences source-backed algorithms over detached topology
//! values. It deliberately stops before live cache installation: cache
//! validity, invalidation, operation capabilities, and commit authority belong
//! exclusively to the parent `cosmolkit` runtime.

use std::ops::{BitAnd, BitOr, BitOrAssign};

use cosmolkit_model::{
    AtomId, QueryStateError, QueryStateRef, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::Hybridization;

use crate::aromaticity::assign_default_aromaticity_with_cached_valence;
use crate::cleanup::cleanup_in_place;
use crate::{
    AromaticityError, AtropisomerError, CleanupError, CleanupParams, ConjugationError,
    HybridizationAssignment, HybridizationError, KekulizeError, KekulizeParams, RadicalError,
    RingFindingError, RingInfo, RingSearchParams, StereoError, ValenceAssignment, ValenceError,
    ValenceModel, ValenceParams, assign_conjugation, assign_hybridization, assign_radicals,
    assign_valence, assign_valence_state_for_atom_from_parts, cleanup, find_sssr, kekulize,
    symmetrized_sssr,
};

use crate::hcount::{AdjustHsError, adjust_hs};

const NAMED_SANITIZE_BITS: u32 = 0x0000_0fff;
const ALL_SANITIZE_BITS: u32 = 0x0fff_ffff;

/// Independently selectable RDKit sanitization stages.
///
/// `ALL` intentionally retains RDKit's reserved high bits. Raw construction
/// accepts either the exact `ALL` sentinel or a combination of named bits;
/// arbitrary unknown-bit patterns are rejected instead of silently doing
/// nothing.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SanitizeOperations(u32);

impl SanitizeOperations {
    pub const NONE: Self = Self(0x000);
    pub const CLEANUP: Self = Self(0x001);
    pub const PROPERTIES: Self = Self(0x002);
    pub const SYMM_RINGS: Self = Self(0x004);
    pub const KEKULIZE: Self = Self(0x008);
    pub const FIND_RADICALS: Self = Self(0x010);
    pub const SET_AROMATICITY: Self = Self(0x020);
    pub const SET_CONJUGATION: Self = Self(0x040);
    pub const SET_HYBRIDIZATION: Self = Self(0x080);
    pub const CLEANUP_CHIRALITY: Self = Self(0x100);
    pub const ADJUST_HS: Self = Self(0x200);
    pub const CLEANUP_ORGANOMETALLICS: Self = Self(0x400);
    pub const CLEANUP_ATROPISOMERS: Self = Self(0x800);
    pub const ALL: Self = Self(ALL_SANITIZE_BITS);

    pub fn from_bits(bits: u32) -> Result<Self, SanitizeError> {
        let unknown_bits = bits & !NAMED_SANITIZE_BITS;
        if unknown_bits != 0 && bits != ALL_SANITIZE_BITS {
            return Err(SanitizeError::InvalidOperations { bits, unknown_bits });
        }
        Ok(Self(bits))
    }

    #[must_use]
    pub const fn bits(self) -> u32 {
        self.0
    }

    #[must_use]
    pub const fn contains(self, operation: Self) -> bool {
        operation.0 != 0 && self.0 & operation.0 == operation.0
    }

    #[must_use]
    pub const fn is_empty(self) -> bool {
        self.0 == 0
    }
}

impl Default for SanitizeOperations {
    fn default() -> Self {
        Self::ALL
    }
}

impl TryFrom<u32> for SanitizeOperations {
    type Error = SanitizeError;

    fn try_from(bits: u32) -> Result<Self, Self::Error> {
        Self::from_bits(bits)
    }
}

impl From<SanitizeOperations> for u32 {
    fn from(operations: SanitizeOperations) -> Self {
        operations.bits()
    }
}

impl BitOr for SanitizeOperations {
    type Output = Self;

    fn bitor(self, rhs: Self) -> Self::Output {
        Self(self.0 | rhs.0)
    }
}

impl BitOrAssign for SanitizeOperations {
    fn bitor_assign(&mut self, rhs: Self) {
        self.0 |= rhs.0;
    }
}

impl BitAnd for SanitizeOperations {
    type Output = Self;

    fn bitand(self, rhs: Self) -> Self::Output {
        Self(self.0 & rhs.0)
    }
}

/// Exact sanitization stage used for source-compatible failure reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum SanitizeStage {
    None = 0x000,
    Cleanup = 0x001,
    Properties = 0x002,
    SymmRings = 0x004,
    Kekulize = 0x008,
    FindRadicals = 0x010,
    SetAromaticity = 0x020,
    SetConjugation = 0x040,
    SetHybridization = 0x080,
    CleanupChirality = 0x100,
    AdjustHs = 0x200,
    CleanupOrganometallics = 0x400,
    CleanupAtropisomers = 0x800,
}

/// Parameters for the detached sanitization pipeline.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SanitizeParams {
    pub operations: SanitizeOperations,
}

impl Default for SanitizeParams {
    fn default() -> Self {
        Self {
            operations: SanitizeOperations::ALL,
        }
    }
}

/// The authoritative detached output of a successful sanitize operation.
#[derive(Debug, Clone, PartialEq)]
pub struct SanitizeAssignment {
    pub topology: TopologyBlock,
    /// The existing final strict property-cache result for the returned topology.
    /// Absent when PROPERTIES is disabled; earlier intermediate assignments do
    /// not establish a valid final-topology cache.
    pub final_valence: Option<ValenceAssignment>,
    /// The same source-stage cache when PROPERTIES is disabled. This is a
    /// detached algorithm value; it does not certify strict runtime state.
    pub non_strict_valence: Option<ValenceAssignment>,
    /// Source `numArom` value from the executed aromaticity assignment, absent
    /// when SET_AROMATICITY is disabled. No ring-count recomputation occurs.
    pub aromatic_ring_count: Option<usize>,
    /// Final ring-state transport for the source ring lifecycle.
    ///
    /// `None` means the FINAL SOURCE-UNINITIALIZED ring state — sanitizeMol's
    /// entry `mol.clearComputedProps()` resets ring info for EVERY operation
    /// mask (including `NONE`), so absence here reports the source's final
    /// uninitialized state, never "caller input preserved" (the Kekulize
    /// ring-update convention uses `None` for that different meaning).
    /// `Some` is the initialized exact final stage state — including
    /// initialized-empty rows — MOVED here after successful final validation.
    /// No final re-find, quality upgrade, or runtime cache authority is
    /// implied.
    pub final_rings: Option<RingInfo>,
}

/// A sanitization failure with the exact active source stage and typed cause.
#[derive(Clone, Debug, PartialEq, thiserror::Error)]
pub enum SanitizeError {
    #[error("molecule property operation failed: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("bond property operation failed: {0}")]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error("invalid sanitize operation bits 0x{bits:08x}; unknown bits 0x{unknown_bits:08x}")]
    InvalidOperations { bits: u32, unknown_bits: u32 },
    #[error("sanitization failed at {stage:?}: invalid topology: {source}")]
    InvalidTopology {
        stage: SanitizeStage,
        source: TopologyValidationError,
    },
    #[error("sanitization failed at {stage:?}: invalid query state: {source}")]
    InvalidQueryState {
        stage: SanitizeStage,
        source: QueryStateError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Cleanup {
        stage: SanitizeStage,
        source: CleanupError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Properties {
        stage: SanitizeStage,
        source: PropertyCacheError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Rings {
        stage: SanitizeStage,
        source: RingFindingError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Kekulize {
        stage: SanitizeStage,
        source: KekulizeError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Radicals {
        stage: SanitizeStage,
        source: RadicalError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Aromaticity {
        stage: SanitizeStage,
        source: AromaticityError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Conjugation {
        stage: SanitizeStage,
        source: ConjugationError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Hybridization {
        stage: SanitizeStage,
        source: HybridizationError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Atropisomers {
        stage: SanitizeStage,
        source: AtropisomerError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    Chirality {
        stage: SanitizeStage,
        source: StereoError,
    },
    #[error("sanitization failed at {stage:?}: {source}")]
    AdjustHs {
        stage: SanitizeStage,
        source: AdjustHsError,
    },
}

/// One source-equivalent chemistry problem found without stopping detection.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ChemistryProblem {
    pub operation: SanitizeStage,
    pub error: ChemistryProblemError,
}

/// Problems that RDKit's detector catches and returns to the caller.
#[derive(Clone, Debug, PartialEq, Eq, thiserror::Error)]
pub enum ChemistryProblemError {
    #[error(transparent)]
    Valence(ValenceError),
    #[error(transparent)]
    Kekulize(KekulizeError),
}

/// Ordered chemistry problems found while inspecting a detached topology.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ChemistryProblemReport {
    pub problems: Vec<ChemistryProblem>,
}

/// Configuration for detached property-cache calculation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PropertyCacheParams {
    pub strict: bool,
}

impl Default for PropertyCacheParams {
    fn default() -> Self {
        Self { strict: true }
    }
}

/// Errors produced while calculating detached property-cache values.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PropertyCacheError {
    #[error("property-cache valence assignment failed: {0}")]
    Valence(#[from] ValenceError),
}

/// The complete value-level result of the RDKit property-cache stage.
///
/// This is intentionally an assignment rather than a runtime cache object.
/// The parent runtime decides whether and when these values become
/// authoritative derived state on a live molecule.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PropertyCacheAssignment {
    pub valence: ValenceAssignment,
}

impl PropertyCacheAssignment {
    #[must_use]
    pub fn valence(&self) -> &ValenceAssignment {
        &self.valence
    }

    #[must_use]
    pub fn into_valence(self) -> ValenceAssignment {
        self.valence
    }
}

/// Calculate RDKit-like explicit valence and implicit-hydrogen facts for a
/// detached topology.
///
/// The source operation writes existing atom cache fields in place. This
/// boundary instead returns two newly allocated vectors so that no algorithm
/// crate can acquire live-cache authority. That preserves source behavior but
/// is a material allocation difference from the source implementation.
pub fn assign_property_cache(
    topology: &TopologyBlock,
    params: &PropertyCacheParams,
) -> Result<PropertyCacheAssignment, PropertyCacheError> {
    // BEGIN RDKIT CPP FUNCTION ROMol::updatePropertyCache
    // RDKit✔️❌: void ROMol::updatePropertyCache(bool strict) {
    // RDKit✔️❌:   for (auto atom : atoms()) {
    // RDKit✔️❌:     atom->updatePropertyCache(strict);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for (auto bond : bonds()) {
    // RDKit✔️❌:     bond->updatePropertyCache(strict);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION ROMol::updatePropertyCache

    // BEGIN RDKIT CPP FUNCTION Atom::updatePropertyCache
    // RDKit✔️❌: void Atom::updatePropertyCache(bool strict) {
    // RDKit✔️❌:   calcExplicitValence(strict);
    // RDKit✔️❌:   calcImplicitValence(strict);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION Atom::updatePropertyCache

    // BEGIN RDKIT CPP FUNCTION Bond::updatePropertyCache
    // RDKit✔️❌: void updatePropertyCache(bool strict = true) { (void)strict; }
    // END RDKIT CPP FUNCTION Bond::updatePropertyCache

    let valence = assign_valence(
        topology,
        &ValenceParams {
            model: ValenceModel::RdkitLike,
            strict: params.strict,
        },
    )?;
    Ok(PropertyCacheAssignment { valence })
}

fn clear_topology_computed_properties(topology: &mut TopologyBlock) -> Result<(), SanitizeError> {
    for atom in &mut topology.atoms {
        atom.clear_computed_props()?;
    }
    for bond in &mut topology.bonds {
        bond.clear_computed_props()?;
    }
    Ok(())
}

fn materialize_radicals(topology: &mut TopologyBlock, values: &[u8]) {
    for (atom, &value) in topology.atoms.iter_mut().zip(values) {
        atom.set_radical_electrons(value);
    }
}

fn materialize_hybridization(topology: &mut TopologyBlock, assignment: &HybridizationAssignment) {
    for (atom, &value) in topology.atoms.iter_mut().zip(&assignment.values) {
        atom.set_hybridization(value);
    }
}

fn current_hybridization(topology: &TopologyBlock) -> HybridizationAssignment {
    HybridizationAssignment {
        values: topology
            .atoms
            .iter()
            .map(cosmolkit_model::Atom::hybridization)
            .collect(),
    }
}

/// One private sanitizer helper owning the source `MolOps::Hybridizations`
/// acquisition for the CLEANUP_ATROPISOMERS stage.
///
/// Branches ALWAYS on the ACTUAL first atom of the CURRENT topology: an
/// atomic-number-zero dummy stays Unspecified
/// (ConjugHybrid.cpp setHybridization), so a dummy-first topology takes the
/// copied branch even when an earlier SET_HYBRIDIZATION assignment exists.
/// The carried assignment is reused ONLY when the source first-atom guard
/// holds AND its values equal the current topology's atom fields. The copied
/// branch runs the existing selected sanitize SET_CONJUGATION|
/// SET_HYBRIDIZATION path (nonstrict property cache, conjugation,
/// hybridization) on an OWNED copy; the copy's rings/props/hybs never leak
/// into the outer topology — only the materialized hybridization values are
/// read out of the copied output.
fn cleanup_stage_hybridizations(
    topology: &TopologyBlock,
    carried: Option<HybridizationAssignment>,
) -> Result<HybridizationAssignment, SanitizeError> {
    #[cfg(test)]
    hybridizations_probe::record_helper();
    // Complete pinned source: MolOps.cpp MolOps::Hybridizations::Hybridizations.
    // RDKit✔️✔️: MolOps::Hybridizations::Hybridizations(const ROMol &mol) {
    // RDKit✔️✔️:   d_hybridizations.clear();
    // RDKit✔️✔️:   // see if the mol already has computed hybridizations:
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (mol.getNumAtoms() == 0) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if topology.atoms.is_empty() {
        return Ok(HybridizationAssignment { values: Vec::new() });
    }
    // RDKit✔️✔️:   if ((*mol.atoms().begin())->getHybridization() !=
    // RDKit✔️✔️:       Atom::HybridizationType::UNSPECIFIED) {
    // RDKit✔️✔️:     for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:       d_hybridizations.push_back((int)atom->getHybridization());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    let first_specified = topology
        .atoms
        .first()
        .is_some_and(|atom| atom.hybridization() != Hybridization::Unspecified);
    if first_specified {
        if let Some(assignment) = carried {
            if assignment.values.len() == topology.atoms.len()
                && assignment
                    .values
                    .iter()
                    .zip(&topology.atoms)
                    .all(|(value, atom)| *value == atom.hybridization())
            {
                return Ok(assignment);
            }
        }
        return Ok(current_hybridization(topology));
    }
    // RDKit✔️✔️:   // compute them in a copy of the mol, so as not to change the mol passed in
    // RDKit✔️✔️:
    // RDKit✔️✔️:   RWMol molCopy(mol);
    // RDKit✔️✔️:   unsigned int operationThatFailed;
    // RDKit✔️✔️:   unsigned int santitizeOps =
    // RDKit✔️✔️:       MolOps::SANITIZE_SETCONJUGATION | MolOps::SANITIZE_SETHYBRIDIZATION;
    // RDKit✔️✔️:   MolOps::sanitizeMol(molCopy, operationThatFailed, santitizeOps);
    // RDKit✔️✔️:   for (auto atom : molCopy.atoms()) {
    // RDKit✔️✔️:     // determine hybridization and remove chiral atoms that are not sp3
    // RDKit✔️✔️:     d_hybridizations.push_back((int)atom->getHybridization());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    // The ONE existing selected-sanitizer owner performs the entry
    // clearComputedProps, the nonstrict property cache, the selected
    // SET_CONJUGATION|SET_HYBRIDIZATION stages and final validation on its
    // OWN clone of this borrowed input (source: RWMol molCopy(mol);
    // sanitizeMol(molCopy, ...)). This helper preclones NOTHING and chains
    // no stages manually. Actual copy costs: the initial O(atoms+bonds)
    // sanitizer working clone plus its existing owner-stage clones (including
    // assign_conjugation), and one O(atoms) field read below. No helper
    // preclone or carried-assignment Vec clone precedes reuse/drop; this does
    // not establish whole-operation allocation parity with source molCopy.
    #[cfg(test)]
    hybridizations_probe::record_copy();
    let selected = sanitize_topology(
        topology,
        &SanitizeParams {
            operations: SanitizeOperations::SET_CONJUGATION | SanitizeOperations::SET_HYBRIDIZATION,
            ..SanitizeParams::default()
        },
    )?;
    #[cfg(test)]
    hybridizations_probe::record_selected_return(&selected);
    // Only the copied output's materialized hybridization fields are read;
    // the copy (rings/props/hybs) never leaks into the outer topology.
    Ok(current_hybridization(&selected.topology))
}

/// Test-only observation points at the actual Hybridizations helper sites.
#[cfg(test)]
pub(crate) mod hybridizations_probe {
    use std::cell::{Cell, RefCell};

    use cosmolkit_types::Hybridization;

    thread_local! {
        static HELPER_CALLS: Cell<u64> = const { Cell::new(0) };
        static COPY_CALLS: Cell<u64> = const { Cell::new(0) };
        static SELECTED_ENTRY_CALLS: Cell<u64> = const { Cell::new(0) };
        static SELECTED_RETURN_HYBS: RefCell<Option<Vec<Hybridization>>> =
            const { RefCell::new(None) };
        static SELECTED_SENTINELS_CLEARED: Cell<Option<(bool, bool)>> = const { Cell::new(None) };
        static SELECTED_FINAL_RINGS_NONE: Cell<Option<bool>> = const { Cell::new(None) };
    }

    /// Computed-property sentinel key shared with the L06 two-call proof.
    pub(crate) const SENTINEL_KEY: &str = "ck-l06-computed-sentinel";

    pub(crate) fn record_helper() {
        HELPER_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn record_copy() {
        COPY_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn helper_calls() -> u64 {
        HELPER_CALLS.with(Cell::get)
    }

    pub(crate) fn copy_calls() -> u64 {
        COPY_CALLS.with(Cell::get)
    }

    pub(crate) fn record_selected_entry() {
        SELECTED_ENTRY_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn selected_entry_calls() -> u64 {
        SELECTED_ENTRY_CALLS.with(Cell::get)
    }

    /// Observation immediately AFTER the nested selected sanitizer returns:
    /// materialized hybs, whether the atom/bond computed sentinels are
    /// absent in the returned copy, and final_rings None.
    pub(crate) fn record_selected_return(assignment: &super::SanitizeAssignment) {
        SELECTED_RETURN_HYBS.with(|cell| {
            *cell.borrow_mut() = Some(
                assignment
                    .topology
                    .atoms
                    .iter()
                    .map(|atom| atom.hybridization())
                    .collect(),
            );
        });
        let atom_cleared = assignment
            .topology
            .atoms
            .first()
            .is_none_or(|atom| !atom.props().contains_key(SENTINEL_KEY.as_bytes()));
        let bond_cleared = assignment
            .topology
            .bonds
            .first()
            .is_none_or(|bond| !bond.props().contains_key(SENTINEL_KEY.as_bytes()));
        SELECTED_SENTINELS_CLEARED.with(|cell| cell.set(Some((atom_cleared, bond_cleared))));
        SELECTED_FINAL_RINGS_NONE.with(|cell| cell.set(Some(assignment.final_rings.is_none())));
    }

    pub(crate) fn selected_return_hybs() -> Option<Vec<Hybridization>> {
        SELECTED_RETURN_HYBS.with(|cell| cell.borrow().clone())
    }

    pub(crate) fn selected_sentinels_cleared() -> Option<(bool, bool)> {
        SELECTED_SENTINELS_CLEARED.with(Cell::take)
    }

    pub(crate) fn selected_final_rings_none() -> Option<bool> {
        SELECTED_FINAL_RINGS_NONE.with(Cell::take)
    }
}

/// Run RDKit's selected sanitization stages in exact source order over a
/// detached topology value.
///
/// The borrowed input is never mutated and a failed stage returns no partial
/// topology. Intermediate assignments reproduce stage dependencies but do not
/// become a second live-cache authority.
pub fn sanitize_topology(
    topology: &TopologyBlock,
    params: &SanitizeParams,
) -> Result<SanitizeAssignment, SanitizeError> {
    sanitize_topology_with_query_state(topology, params, None)
}

#[doc(hidden)]
pub fn sanitize_topology_with_query_state(
    topology: &TopologyBlock,
    params: &SanitizeParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<SanitizeAssignment, SanitizeError> {
    #[cfg(test)]
    if params.operations
        == (SanitizeOperations::SET_CONJUGATION | SanitizeOperations::SET_HYBRIDIZATION)
    {
        hybridizations_probe::record_selected_entry();
    }
    topology
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;
    if let Some(state) = query_state {
        state.validate_for_topology(topology).map_err(|source| {
            SanitizeError::InvalidQueryState {
                stage: SanitizeStage::None,
                source,
            }
        })?;
    }

    // Complete pinned source: MolOps.cpp::sanitizeMol(RWMol &, unsigned int &, unsigned int).
    // RDKit✔️❌: void sanitizeMol(RWMol &mol, unsigned int &operationThatFailed,
    // RDKit✔️❌:                  unsigned int sanitizeOps) {
    // RDKit✔️❌:   // clear out any cached properties
    // RDKit✔️❌:   mol.clearComputedProps();
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUP;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     // clean up things like nitro groups
    // RDKit✔️❌:     cleanUp(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // fix things like non-metal to metal bonds that should be dative.
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUP_ORGANOMETALLICS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanUpOrganometallics(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // update computed properties on atoms and bonds:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     mol.updatePropertyCache(true);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SYMMRINGS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     VECT_INT_VECT arings;
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(mol, arings);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // kekulizations
    // RDKit✔️❌:   operationThatFailed = SANITIZE_KEKULIZE;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     Kekulize(mol, true, false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // look for radicals:
    // RDKit✔️❌:   // We do this now because we need to know
    // RDKit✔️❌:   // that the N in [N]1C=CC=C1 has a radical
    // RDKit✔️❌:   // before we move into setAromaticity().
    // RDKit✔️❌:   // It's important that this happen post-Kekulization
    // RDKit✔️❌:   // because there's no way of telling what to do
    // RDKit✔️❌:   // with the same molecule if it's in the form
    // RDKit✔️❌:   // [n]1cccc1
    // RDKit✔️❌:   operationThatFailed = SANITIZE_FINDRADICALS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     assignRadicals(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // then do aromaticity perception
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETAROMATICITY;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setAromaticity(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // set conjugation
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETCONJUGATION;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setConjugation(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // set hybridization
    // RDKit✔️❌:   operationThatFailed = SANITIZE_SETHYBRIDIZATION;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     setHybridization(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUPATROPISOMERS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanupAtropisomers(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // remove bogus chirality specs:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_CLEANUPCHIRALITY;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     cleanupChirality(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // adjust Hydrogen counts:
    // RDKit✔️❌:   operationThatFailed = SANITIZE_ADJUSTHS;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     adjustHs(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // now that everything has been cleaned up, go through and check/update the
    // RDKit✔️❌:   // computed valences on atoms and bonds one more time
    // RDKit✔️❌:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️❌:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️❌:     mol.updatePropertyCache(true);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   operationThatFailed = 0;
    // RDKit✔️❌: }
    // The source mutates one `RWMol`. This detached boundary clones the full
    // topology to isolate the borrowed input. Cleanup stages now mutate that
    // one detached attempt, retaining their input/output validation and stage
    // order without cloning the topology again. Other owner-stage overhead
    // remains; this does not promote the full sanitizer's complexity markers.

    let operations = params.operations;
    let mut working = topology.clone();
    clear_topology_computed_properties(&mut working)?;

    if operations.contains(SanitizeOperations::CLEANUP) {
        cleanup_in_place(
            &mut working,
            &CleanupParams {
                charge_normalization: true,
                organometallics: false,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::Cleanup,
            source,
        })?;
    }

    if operations.contains(SanitizeOperations::CLEANUP_ORGANOMETALLICS) {
        cleanup_in_place(
            &mut working,
            &CleanupParams {
                charge_normalization: false,
                organometallics: true,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::CleanupOrganometallics,
            source,
        })?;
    }

    let mut valence = assign_property_cache(
        &working,
        &PropertyCacheParams {
            strict: operations.contains(SanitizeOperations::PROPERTIES),
        },
    )
    .map_err(|source| SanitizeError::Properties {
        stage: SanitizeStage::Properties,
        source,
    })?
    .into_valence();

    let mut rings: Option<RingInfo> = None;
    let mut aromatic_ring_count = None;
    if operations.contains(SanitizeOperations::SYMM_RINGS) {
        #[cfg(test)]
        final_rings_probe::record_symm_acquisition();
        rings = Some(
            symmetrized_sssr(&working, &RingSearchParams::default()).map_err(|source| {
                SanitizeError::Rings {
                    stage: SanitizeStage::SymmRings,
                    source,
                }
            })?,
        );
        // Actual SYMM-stage completion observation: both nonempty outer
        // buffer pointers and the complete ordered state, immediately
        // after the local carrier assignment and BEFORE any K/A/T stage
        // (the nested selected sanitize for T carries no S bit and never
        // reaches this site).
        #[cfg(test)]
        final_rings_probe::record_early_s(rings.as_ref().expect("SYMM state"));
    }

    if operations.contains(SanitizeOperations::KEKULIZE) {
        // Source retry controller shares the one mutable detached attempt;
        // only true MolSanitizeException equivalents trigger canonical retry.
        // Borrow current intermediate valence/rings, without topology/ring clone.
        crate::kekulize::source_kekulize_for_sanitize(
            &mut working,
            &mut valence,
            &mut rings,
            query_state,
        )
        .map_err(|source| SanitizeError::Kekulize {
            stage: SanitizeStage::Kekulize,
            source,
        })?;
    }

    if operations.contains(SanitizeOperations::FIND_RADICALS) {
        let assignment = assign_radicals(&working).map_err(|source| SanitizeError::Radicals {
            stage: SanitizeStage::FindRadicals,
            source,
        })?;
        materialize_radicals(&mut working, &assignment.radical_electrons);
    }

    if operations.contains(SanitizeOperations::SET_AROMATICITY) {
        // BEGIN RDKIT CPP FUNCTION: Aromaticity.cpp :: setAromaticity ring guard
        // RDKit✔️✔️: void setAromaticity(RWMol &mol, AromaticityModel model, int (*func)(RWMol &)) {
        // RDKit✔️✔️:   // This function used to check if the input molecule came
        // RDKit✔️✔️:   // with aromaticity information, assumed it is correct and
        // RDKit✔️✔️:   // did not touch it. Now it ignores that information entirely.
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // first find the all the simple rings in the molecule
        // RDKit✔️✔️:   VECT_INT_VECT srings;
        // RDKit✔️✔️:   if (mol.getRingInfo()->isInitialized()) {
        // RDKit✔️✔️:     srings = mol.getRingInfo()->atomRings();
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     MolOps::symmetrizeSSSR(mol, srings);
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION: Aromaticity.cpp :: setAromaticity ring guard
        // Behavior review: the source borrows the mol's initialized rows
        // (whatever finder produced them — SSSR stays SSSR, no quality
        // upgrade) and only symmetrizes when uninitialized. The carrier
        // mirrors that: initialized ⇒ BORROW the one owned state; absent ⇒
        // one fresh symmetrized acquisition owned by this stage.
        // Complexity review: the borrow replaces the previous full-RingInfo
        // clone on both arms; the only allocation is the single fresh
        // symmetrized result when the carrier is absent.
        if rings.is_none() {
            rings = Some(
                symmetrized_sssr(&working, &RingSearchParams::default()).map_err(|source| {
                    SanitizeError::Rings {
                        stage: SanitizeStage::SetAromaticity,
                        source,
                    }
                })?,
            );
        }
        let ring_assignment = rings.as_ref().expect("aromaticity state present");
        // Behavior: RDKit SetAromaticity consumes the existing atom cache;
        // `valence` is the same strict Properties assignment with only the
        // source-refreshed Kekulize N/P rows applied above.
        // Complexity: borrow those rows into the default aromaticity owner,
        // avoiding a second O(V+E) assignment while keeping its validation.
        let aromaticity = assign_default_aromaticity_with_cached_valence(
            &working,
            &ring_assignment,
            &valence,
            query_state,
        )
        .map_err(|source| SanitizeError::Aromaticity {
            stage: SanitizeStage::SetAromaticity,
            source,
        })?;
        // RDKit✔️✔️:   mol.setProp(common_properties::numArom, narom, true);
        // Transport the scalar computed by the aromaticity owner, O(1).
        aromatic_ring_count = Some(aromaticity.aromatic_ring_count);
        working = aromaticity.topology;
    }

    if operations.contains(SanitizeOperations::SET_CONJUGATION) {
        working = assign_conjugation(&working, &valence).map_err(|source| {
            SanitizeError::Conjugation {
                stage: SanitizeStage::SetConjugation,
                source,
            }
        })?;
    }

    let mut hybridization = None;
    if operations.contains(SanitizeOperations::SET_HYBRIDIZATION) {
        let assignment = assign_hybridization(&working, &valence).map_err(|source| {
            SanitizeError::Hybridization {
                stage: SanitizeStage::SetHybridization,
                source,
            }
        })?;
        materialize_hybridization(&mut working, &assignment);
        hybridization = Some(assignment);
    }

    if operations.contains(SanitizeOperations::CLEANUP_ATROPISOMERS) {
        // Complete pinned source: MolOps.cpp cleanupAtropisomers wrapper.
        // RDKit✔️✔️: void cleanupAtropisomers(RWMol &mol) {
        // RDKit✔️✔️:   auto hybs = MolOps::Hybridizations(mol);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   MolOps::cleanupAtropisomers(mol, hybs);
        // RDKit✔️✔️: }
        // Hybridizations are ALWAYS constructed first (even without tags);
        // the owned-ring-state wrapper then prevalidates topology/hybridization
        // before its tag loop and maps its two private causes structurally.
        let assignment = cleanup_stage_hybridizations(&working, hybridization.take())?;
        let (output, ring_state) =
            crate::atropisomer::cleanup_invalid_atropisomers_with_ring_state(
                &working,
                &assignment,
                rings.take(),
            )
            .map_err(|source| match source {
                crate::atropisomer::AtropisomerCleanupError::Rings(inner) => SanitizeError::Rings {
                    stage: SanitizeStage::CleanupAtropisomers,
                    source: inner,
                },
                crate::atropisomer::AtropisomerCleanupError::Algorithm(inner) => {
                    SanitizeError::Atropisomers {
                        stage: SanitizeStage::CleanupAtropisomers,
                        source: inner,
                    }
                }
            })?;
        working = output;
        rings = ring_state;
    }

    if operations.contains(SanitizeOperations::CLEANUP_CHIRALITY) {
        working =
            crate::structure_tags::cleanup_chirality(&working, &valence).map_err(|source| {
                SanitizeError::Chirality {
                    stage: SanitizeStage::CleanupChirality,
                    source,
                }
            })?;
    }

    if operations.contains(SanitizeOperations::ADJUST_HS) {
        let assignment =
            adjust_hs(&working, &valence).map_err(|source| SanitizeError::AdjustHs {
                stage: SanitizeStage::AdjustHs,
                source,
            })?;
        working = assignment.topology;
        valence = assignment.valence;
    }

    // RDKit✔️✔️:   operationThatFailed = SANITIZE_PROPERTIES;
    // RDKit✔️✔️:   if (sanitizeOps & operationThatFailed) {
    // RDKit✔️✔️:     mol.updatePropertyCache(true);
    // RDKit✔️✔️:   }
    // Behavior review: transport the final source-stage assignment already
    // computed here, never the earlier intermediate property-cache state.
    // PROPERTIES disabled leaves no final assignment to certify or install.
    // Complexity review: move the existing assignment into the return value;
    // no new property-cache evaluation, traversal, or vector clone is added.
    if operations.contains(SanitizeOperations::PROPERTIES) {
        valence = assign_property_cache(&working, &PropertyCacheParams { strict: true })
            .map_err(|source| SanitizeError::Properties {
                stage: SanitizeStage::Properties,
                source,
            })?
            .into_valence();
    }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // Move the final source-stage rows to exactly one typed return field.
    // This preserves the strict-certification distinction without recalculating
    // or cloning the rows when PROPERTIES is disabled.
    let (final_valence, non_strict_valence) = if operations.contains(SanitizeOperations::PROPERTIES)
    {
        (Some(valence), None)
    } else {
        (None, Some(valence))
    };

    working
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;
    // RDKit✔️✔️:   operationThatFailed = 0;
    // Behavior review: the source leaves the mol's RingInfo in whatever state
    // the executed stages produced (entry clearComputedProps reset it first);
    // the local carrier is MOVED into the return only after this final
    // validation succeeds, so failures publish no partial ring result.
    // Complexity review: one move of the owned carrier; no clone, re-find,
    // or quality upgrade is added at the transport point.
    #[cfg(test)]
    final_rings_probe::record(&rings);
    Ok(SanitizeAssignment {
        topology: working,
        final_valence,
        non_strict_valence,
        aromatic_ring_count,
        final_rings: rings,
    })
}

/// Test-only observation at the ACTUAL final_rings move point: captures the
/// carrier's buffer pointer and a clone of its complete state immediately
/// BEFORE the move. Production builds contain no probe.
#[cfg(test)]
pub(crate) mod final_rings_probe {
    use std::cell::{Cell, RefCell};

    use crate::RingInfo;

    thread_local! {
        static LAST_POINTER: Cell<usize> = const { Cell::new(0) };
        static LAST_STATE: RefCell<Option<RingInfo>> = const { RefCell::new(None) };
        static SYMM_ACQUISITIONS: Cell<u64> = const { Cell::new(0) };
        static EARLY_S_ATOM_PTR: Cell<usize> = const { Cell::new(0) };
        static EARLY_S_BOND_PTR: Cell<usize> = const { Cell::new(0) };
        static EARLY_S_STATE: RefCell<Option<RingInfo>> = const { RefCell::new(None) };
    }

    pub(crate) fn record(rings: &Option<RingInfo>) {
        LAST_POINTER.with(|cell| {
            cell.set(
                rings
                    .as_ref()
                    .map(|state| state.atom_rings().as_ptr() as usize)
                    .unwrap_or(0),
            )
        });
        LAST_STATE.with(|cell| {
            *cell.borrow_mut() = rings.clone();
        });
    }

    pub(crate) fn last_pointer() -> usize {
        LAST_POINTER.with(Cell::get)
    }

    pub(crate) fn last_state() -> Option<RingInfo> {
        LAST_STATE.with(|cell| cell.borrow().clone())
    }

    pub(crate) fn record_symm_acquisition() {
        SYMM_ACQUISITIONS.with(|calls| calls.set(calls.get() + 1));
    }

    pub(crate) fn symm_acquisitions() -> u64 {
        SYMM_ACQUISITIONS.with(Cell::get)
    }

    /// Both nonempty outer-buffer pointers plus the complete ordered state
    /// captured at ACTUAL SYMM-stage completion.
    pub(crate) fn record_early_s(rings: &RingInfo) {
        EARLY_S_ATOM_PTR.with(|cell| cell.set(rings.atom_rings().as_ptr() as usize));
        EARLY_S_BOND_PTR.with(|cell| cell.set(rings.bond_rings().as_ptr() as usize));
        EARLY_S_STATE.with(|cell| *cell.borrow_mut() = Some(rings.clone()));
    }

    pub(crate) fn early_s_pointers() -> (usize, usize) {
        (
            EARLY_S_ATOM_PTR.with(Cell::get),
            EARLY_S_BOND_PTR.with(Cell::get),
        )
    }

    pub(crate) fn early_s_state() -> Option<RingInfo> {
        EARLY_S_STATE.with(|cell| cell.borrow().clone())
    }
}

fn detection_property_cache(
    topology: &mut TopologyBlock,
    strict: bool,
    report: &mut ChemistryProblemReport,
) -> Result<ValenceAssignment, SanitizeError> {
    // The changed detectChemistryProblems owner needs the actual cache input
    // for both attempts. Reuse the existing Atom updatePropertyCache kernel;
    // each completed explicit store survives an implicit failure before next atom.
    for index in 0..topology.atoms.len() {
        match crate::valence::update_source_atom_cache(topology, AtomId::new(index), strict) {
            Ok(()) => {}
            Err(source @ ValenceError::InvalidValence { .. }) if strict => {
                report.problems.push(ChemistryProblem {
                    operation: SanitizeStage::Properties,
                    error: ChemistryProblemError::Valence(source),
                });
            }
            Err(source) => {
                return Err(SanitizeError::Properties {
                    stage: SanitizeStage::Properties,
                    source: PropertyCacheError::Valence(source),
                });
            }
        }
    }
    Ok(ValenceAssignment {
        explicit_valence: topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().explicit_valence))
            .collect(),
        implicit_hydrogens: topology
            .atoms
            .iter()
            .map(|a| i32::from(a.source_valence_facts().implicit_valence))
            .collect(),
    })
}

fn is_source_kekulize_problem(error: &KekulizeError) -> bool {
    match error {
        KekulizeError::AromaticAtomOutsideRing { .. }
        | KekulizeError::NotKekulizable { .. }
        | KekulizeError::PostconditionValenceMismatch { .. } => true,
        KekulizeError::Valence(source) => {
            matches!(source, ValenceError::InvalidValence { .. })
        }
        _ => false,
    }
}

/// Detect all source-equivalent atom-valence problems followed by at most one
/// source-equivalent kekulization problem.
///
/// The borrowed input is never mutated. Errors outside RDKit's caught
/// `MolSanitizeException` boundary remain top-level [`SanitizeError`] values.
pub fn detect_chemistry_problems(
    topology: &TopologyBlock,
    params: &SanitizeParams,
) -> Result<ChemistryProblemReport, SanitizeError> {
    topology
        .validate()
        .map_err(|source| SanitizeError::InvalidTopology {
            stage: SanitizeStage::None,
            source,
        })?;

    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: detectChemistryProblems (2026.03.6 complete)
    // RDKit❗❌: std::vector<std::unique_ptr<MolSanitizeException>> detectChemistryProblems(
    // RDKit❗❌:     const ROMol &imol, unsigned int sanitizeOps) {
    // RDKit❗❌:   RWMol mol(imol);
    // RDKit❗❌:   std::vector<std::unique_ptr<MolSanitizeException>> res;
    // RDKit❗❌:
    // RDKit❗❌:   // clear out any cached properties
    // RDKit❗❌:   mol.clearComputedProps();
    // RDKit❗❌:
    // RDKit❗❌:   int operation;
    // RDKit❗❌:   operation = SANITIZE_CLEANUP;
    // RDKit❗❌:   if (sanitizeOps & operation) {
    // RDKit❗❌:     // clean up things like nitro groups
    // RDKit❗❌:     cleanUp(mol);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // update computed properties on atoms and bonds:
    // RDKit❗❌:   operation = SANITIZE_PROPERTIES;
    // RDKit❗❌:   if (sanitizeOps & operation) {
    // RDKit❗❌:     for (auto &atom : mol.atoms()) {
    // RDKit❗❌:       try {
    // RDKit❗❌:         bool strict = true;
    // RDKit❗❌:         atom->updatePropertyCache(strict);
    // RDKit❗❌:       } catch (const MolSanitizeException &e) {
    // RDKit❗❌:         res.emplace_back(e.copy());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   } else {
    // RDKit❗❌:     mol.updatePropertyCache(false);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // kekulizations
    // RDKit❗❌:   operation = SANITIZE_KEKULIZE;
    // RDKit❗❌:   if (sanitizeOps & operation) {
    // RDKit❗❌:     try {
    // RDKit❗❌:       kekulizeForSanitize(mol);
    // RDKit❗❌:     } catch (const MolSanitizeException &e) {
    // RDKit❗❌:       res.emplace_back(e.copy());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: detectChemistryProblems (2026.03.6 complete)
    let operations = params.operations;
    let mut working = topology.clone();
    clear_topology_computed_properties(&mut working)?;
    let mut report = ChemistryProblemReport::default();

    if operations.contains(SanitizeOperations::CLEANUP) {
        working = cleanup(
            &working,
            &CleanupParams {
                charge_normalization: true,
                organometallics: false,
            },
        )
        .map_err(|source| SanitizeError::Cleanup {
            stage: SanitizeStage::Cleanup,
            source,
        })?;
    }

    let mut valence = detection_property_cache(
        &mut working,
        operations.contains(SanitizeOperations::PROPERTIES),
        &mut report,
    )?;
    if operations.contains(SanitizeOperations::KEKULIZE) {
        let mut rings = None;
        if let Err(source) = crate::kekulize::source_kekulize_for_sanitize(
            &mut working,
            &mut valence,
            &mut rings,
            None,
        ) {
            if is_source_kekulize_problem(&source) {
                report.problems.push(ChemistryProblem {
                    operation: SanitizeStage::Kekulize,
                    error: ChemistryProblemError::Kekulize(source),
                });
            } else {
                return Err(SanitizeError::Kekulize {
                    stage: SanitizeStage::Kekulize,
                    source,
                });
            }
        }
    }

    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondOrder, Element};

    fn ethanol_topology() -> TopologyBlock {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        let adjacency = AdjacencyList::from_topology(3, &bonds);
        TopologyBlock {
            atoms,
            bonds,
            adjacency,
            stereo_groups: Vec::new(),
            substance_groups: Vec::new(),
        }
    }

    #[test]
    fn detached_assignment_matches_formal_valence_boundary() {
        let topology = ethanol_topology();
        let assignment =
            assign_property_cache(&topology, &PropertyCacheParams { strict: false }).unwrap();
        assert_eq!(assignment.valence.explicit_valence, vec![1, 2, 1]);
        assert_eq!(assignment.valence.implicit_hydrogens, vec![3, 2, 1]);
    }
}

// L06: the frozen eleven-call Hybridizations helper/composition table.
// Six tagged CC cells (three initial vectors x T-only/HYB|T), four no-tag
// controls (empty + CC x the three vectors) and the dummy-first eleventh
// call proving the copied branch still runs despite a carried HYB
// assignment. Counters are thread-local, never reset; all observations are
// per-call deltas. Full-topology snapshot equality is the computed-prop
// sentinel (properties participate in equality).
#[cfg(test)]
mod cleanup_hybridizations_tests {
    use super::cleanup_stage_hybridizations;
    use super::hybridizations_probe as hyb_probe;
    use crate::atropisomer::cleanup_invalid_atropisomers_with_ring_state;
    use crate::atropisomer::cleanup_ring_state_probe as atrop_probe;
    use crate::sanitize_topology;
    use crate::{HybridizationAssignment, SanitizeOperations, SanitizeParams};
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
    use cosmolkit_types::{BondOrder, BondStereo, Element, Hybridization};

    fn cc(initial: [Hybridization; 2], stereo: BondStereo) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C).with_hybridization(initial[0]),
                ),
                Atom::from_spec(
                    AtomId::new(1),
                    AtomSpec::new(Element::C).with_hybridization(initial[1]),
                ),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_stereo(stereo),
            )],
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn hybs(values: &[Hybridization]) -> HybridizationAssignment {
        HybridizationAssignment {
            values: values.to_vec(),
        }
    }

    // HYB materialization for the HYB|T cells, exactly as the
    // SET_HYBRIDIZATION stage would leave the outer topology.
    fn materialized_cc(initial: [Hybridization; 2]) -> (TopologyBlock, HybridizationAssignment) {
        let mut working = cc(initial, BondStereo::AtropCw);
        let valence =
            super::assign_property_cache(&working, &super::PropertyCacheParams { strict: false })
                .unwrap()
                .into_valence();
        let assignment = crate::assign_hybridization(&working, &valence).unwrap();
        super::materialize_hybridization(&mut working, &assignment);
        (working, assignment)
    }

    #[test]
    fn sanitize_ring_l06_hybridizations_eleven_call_table() {
        let initials = [
            [Hybridization::Sp2, Hybridization::Unspecified],
            [Hybridization::Unspecified, Hybridization::Sp2],
            [Hybridization::Sp2, Hybridization::Sp2],
        ];
        let mut calls = 0usize;

        // --- Six tagged CC cells -----------------------------------------
        for (index, initial) in initials.iter().enumerate() {
            // T-only: carrier absent, input fields are the initial vector.
            let topology = cc(*initial, BondStereo::AtropCw);
            let snapshot = topology.clone();
            let helper_before = hyb_probe::helper_calls();
            let copy_before = hyb_probe::copy_calls();
            let find_before = atrop_probe::find_calls();
            let group_before = atrop_probe::group_cleanup_calls();
            let assignment = cleanup_stage_hybridizations(&topology, None).unwrap();
            calls += 1;
            assert_eq!(hyb_probe::helper_calls() - helper_before, 1, "t{index}");
            let copy_delta = hyb_probe::copy_calls() - copy_before;
            let expected_copy = if index == 1 { 1 } else { 0 };
            assert_eq!(copy_delta, expected_copy, "t{index}: copied branch");
            let expected_values: &[Hybridization] = if index == 1 {
                // Copied branch computes CH3-CH3-like carbons as Sp3/Sp3.
                &[Hybridization::Sp3, Hybridization::Sp3]
            } else {
                initial.as_slice()
            };
            assert_eq!(assignment.values, expected_values, "t{index}: values");
            assert_eq!(topology, snapshot, "t{index}: input mutated");
            // Composition: the wrapper consumes the helper output with an
            // absent carrier (outer entry reset semantics) and acquires the
            // initialized-empty SSSR exactly once for the tagged bond.
            let (output, state) =
                cleanup_invalid_atropisomers_with_ring_state(&topology, &assignment, None).unwrap();
            assert_eq!(
                atrop_probe::find_calls() - find_before,
                1,
                "t{index}: acquire"
            );
            let clears = index != 2;
            if clears {
                assert_eq!(output.bonds[0].stereo(), BondStereo::None, "t{index}");
                assert_eq!(
                    atrop_probe::group_cleanup_calls() - group_before,
                    1,
                    "t{index}: group cleanup"
                );
            } else {
                assert_eq!(output.bonds[0].stereo(), BondStereo::AtropCw, "t{index}");
                assert_eq!(
                    atrop_probe::group_cleanup_calls() - group_before,
                    0,
                    "t{index}: group cleanup"
                );
            }
            let acquired = state.unwrap();
            assert_eq!(acquired.find_type(), crate::RingFindType::Sssr, "t{index}");
            assert!(acquired.atom_rings().is_empty(), "t{index}: SSSR empty");

            // HYB|T: the outer topology carries materialized [Sp3,Sp3] and
            // the carried assignment matches, so the guard-true reuse path
            // applies (copied branch 0) and every tag clears.
            let (working, carried) = materialized_cc(*initial);
            let working_snapshot = working.clone();
            let copy_before = hyb_probe::copy_calls();
            let find_before = atrop_probe::find_calls();
            let assignment = cleanup_stage_hybridizations(&working, Some(carried)).unwrap();
            calls += 1;
            assert_eq!(
                hyb_probe::copy_calls() - copy_before,
                0,
                "h{index}: copied branch"
            );
            assert_eq!(
                assignment.values,
                &[Hybridization::Sp3, Hybridization::Sp3],
                "h{index}: values"
            );
            assert_eq!(working, working_snapshot, "h{index}: input mutated");
            let (output, _) =
                cleanup_invalid_atropisomers_with_ring_state(&working, &assignment, None).unwrap();
            assert_eq!(
                atrop_probe::find_calls() - find_before,
                1,
                "h{index}: acquire"
            );
            assert_eq!(output.bonds[0].stereo(), BondStereo::None, "h{index}");
        }

        // --- Four no-tag controls ----------------------------------------
        let empty =
            TopologyBlock::try_from_parts(Vec::new(), Vec::new(), Vec::new(), Vec::new()).unwrap();
        let copy_before = hyb_probe::copy_calls();
        let assignment = cleanup_stage_hybridizations(&empty, None).unwrap();
        calls += 1;
        assert!(assignment.values.is_empty(), "empty: values");
        assert_eq!(hyb_probe::copy_calls() - copy_before, 0, "empty: copied");

        for (index, initial) in initials.iter().enumerate() {
            let topology = cc(*initial, BondStereo::None);
            let snapshot = topology.clone();
            let copy_before = hyb_probe::copy_calls();
            let find_before = atrop_probe::find_calls();
            let assignment = cleanup_stage_hybridizations(&topology, None).unwrap();
            calls += 1;
            // No tags: fields preserved; acquire0; copied branch ONLY for
            // the middle CC (first atom Unspecified).
            let expected_copy = if index == 1 { 1 } else { 0 };
            assert_eq!(
                hyb_probe::copy_calls() - copy_before,
                expected_copy,
                "n{index}"
            );
            assert_eq!(
                atrop_probe::find_calls() - find_before,
                0,
                "n{index}: acquire"
            );
            assert_eq!(topology, snapshot, "n{index}: fields mutated");
            if index == 1 {
                assert_eq!(
                    assignment.values,
                    &[Hybridization::Sp3, Hybridization::Sp3],
                    "n{index}: values"
                );
            } else {
                assert_eq!(assignment.values, initial.as_slice(), "n{index}: values");
            }
        }

        // --- Eleventh call: dummy-first [*,C], no tags, carried HYB ------
        let dummy_first = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::DUMMY)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let snapshot = dummy_first.clone();
        // After HYB the dummy stays Unspecified (ConjugHybrid dummy rule);
        // the carried assignment therefore matches the outer fields, but the
        // ACTUAL first atom is Unspecified, so the source takes the COPIED
        // branch anyway. This is the guard the previous unconditional
        // prior-assignment reuse skipped.
        let carried = hybs(&[Hybridization::Unspecified, Hybridization::Sp3]);
        let helper_before = hyb_probe::helper_calls();
        let copy_before = hyb_probe::copy_calls();
        let find_before = atrop_probe::find_calls();
        let assignment = cleanup_stage_hybridizations(&dummy_first, Some(carried)).unwrap();
        calls += 1;
        assert_eq!(hyb_probe::helper_calls() - helper_before, 1, "11th: helper");
        assert_eq!(
            hyb_probe::copy_calls() - copy_before,
            1,
            "11th: copied branch"
        );
        assert_eq!(atrop_probe::find_calls() - find_before, 0, "11th: acquire");
        assert_eq!(
            assignment.values,
            &[Hybridization::Unspecified, Hybridization::Sp3],
            "11th: values"
        );
        assert_eq!(dummy_first, snapshot, "11th: input mutated");

        assert_eq!(calls, 11, "exact census");

        // Composition sanity for the eleventh graph: the PUBLIC T-only and
        // HYB|T sanitize paths still succeed on this graph.
        sanitize_topology(
            &dummy_first,
            &SanitizeParams {
                operations: SanitizeOperations::CLEANUP_ATROPISOMERS,
                ..SanitizeParams::default()
            },
        )
        .unwrap();
        sanitize_topology(
            &dummy_first,
            &SanitizeParams {
                operations: SanitizeOperations::SET_HYBRIDIZATION
                    | SanitizeOperations::CLEANUP_ATROPISOMERS,
                ..SanitizeParams::default()
            },
        )
        .unwrap();
    }

    // Fixed closure proof: the copied branch goes through the ACTUAL
    // selected sanitize entry (exact SET_CONJUGATION|SET_HYBRIDIZATION
    // mask, counted at that entry, not at the helper), which owns the
    // cloned copy; the returned copy has materialized [Sp3,Sp3], cleared
    // atom/bond computed sentinels and final_rings None. The ORIGINAL
    // topology/props and the supplied carried assignment stay unchanged.
    #[test]
    fn sanitize_ring_l06_selected_sanitizer_entry_two_calls() {
        use super::hybridizations_probe::{self, SENTINEL_KEY};
        let mut graph = cc(
            [Hybridization::Unspecified, Hybridization::Sp2],
            BondStereo::None,
        );
        // Nondefault computed-property sentinels through the existing
        // owners; the selected entry clearComputedProps must clear them
        // INSIDE the copy only.
        graph.atoms[0]
            .set_computed_prop(SENTINEL_KEY, "atom-sentinel")
            .unwrap();
        graph.bonds[0]
            .set_computed_prop(SENTINEL_KEY, "bond-sentinel")
            .unwrap();
        let snapshot = graph.clone();
        let carried_values = [
            None,
            Some(HybridizationAssignment {
                values: vec![Hybridization::Sp3, Hybridization::Sp3],
            }),
        ];
        let mut calls = 0usize;
        for carried in carried_values {
            let carried_snapshot = carried.clone();
            let entry_before = hybridizations_probe::selected_entry_calls();
            let helper_before = hybridizations_probe::helper_calls();
            let copy_before = hybridizations_probe::copy_calls();
            let assignment = cleanup_stage_hybridizations(&graph, carried.clone()).unwrap();
            calls += 1;
            assert_eq!(
                hybridizations_probe::helper_calls() - helper_before,
                1,
                "helper entry"
            );
            assert_eq!(
                hybridizations_probe::copy_calls() - copy_before,
                1,
                "copied branch"
            );
            // The ACTUAL selected sanitize entry ran exactly once with the
            // exact CONJ|HYB mask for this helper call.
            assert_eq!(
                hybridizations_probe::selected_entry_calls() - entry_before,
                1,
                "selected sanitize entry"
            );
            assert_eq!(
                assignment.values,
                &[Hybridization::Sp3, Hybridization::Sp3],
                "materialized hybs"
            );
            // Observation recorded immediately after the nested sanitizer
            // returned: materialized hybs, cleared sentinels, no rings.
            assert_eq!(
                hybridizations_probe::selected_return_hybs(),
                Some(vec![Hybridization::Sp3, Hybridization::Sp3]),
                "returned-copy hybs"
            );
            assert_eq!(
                hybridizations_probe::selected_sentinels_cleared(),
                Some((true, true)),
                "returned-copy sentinels cleared"
            );
            assert_eq!(
                hybridizations_probe::selected_final_rings_none(),
                Some(true),
                "returned-copy final_rings None"
            );
            // Original topology/props and supplied assignment unchanged.
            assert_eq!(graph, snapshot, "original topology/props mutated");
            assert_eq!(carried_snapshot, carried, "supplied assignment consumed");
        }
        assert_eq!(calls, 2, "exact census");
    }
}

// L07: composed CLEANUP_ATROPISOMERS regressions through the public
// sanitize entry. All success/topology assertions are unconditional.
// Probe counters are thread-local, never reset; observations are deltas.
#[cfg(test)]
mod cleanup_composed_tests {
    use crate::atropisomer::cleanup_ring_state_probe as atrop_probe;
    use crate::{RingFindType, SanitizeOperations, SanitizeParams, sanitize_topology};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, StereoGroup, StereoGroupKind, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, BondStereo, Element, Hybridization};

    fn cycle(n: usize, stereo: BondStereo, aromatic: bool) -> TopologyBlock {
        debug_assert!(n >= 3, "cycle fixture requires n >= 3");
        let atoms = (0..n)
            .map(|id| {
                let mut spec = AtomSpec::new(Element::C);
                if aromatic {
                    spec = spec.with_aromatic(true);
                }
                Atom::from_spec(AtomId::new(id), spec)
            })
            .collect::<Vec<_>>();
        let bonds = (0..n)
            .map(|id| {
                // n == 8 uses alternating Single/Double orders
                // (cyclooctatetraene): the copied branch computes Sp2 at
                // every atom, which is the composed macrocycle-retention
                // case; plain sp3 cycles clear regardless of ring size.
                let order = if n == 8 && id % 2 == 1 {
                    BondOrder::Double
                } else if aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                };
                let mut spec = BondSpec::new(AtomId::new(id), AtomId::new((id + 1) % n), order);
                if aromatic {
                    spec = spec.with_aromatic(true);
                }
                if id == 0 {
                    spec = spec.with_stereo(stereo);
                }
                Bond::from_spec(BondId::new(id), spec)
            })
            .collect::<Vec<_>>();
        let sentinel = StereoGroup::new(
            StereoGroupKind::Or,
            (0..n).map(AtomId::new).collect(),
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(7);
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), vec![sentinel]).unwrap()
    }

    fn cc(stereo: BondStereo) -> TopologyBlock {
        // One SINGLE bond between two carbons (no duplicate parallel edge).
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single).with_stereo(stereo),
        )];
        let sentinel = StereoGroup::new(
            StereoGroupKind::Or,
            vec![AtomId::new(0), AtomId::new(1)],
            Vec::new(),
        )
        .expect("valid distinct stereo members")
        .with_id(7);
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), vec![sentinel]).unwrap()
    }

    fn run(topology: &TopologyBlock, operations: SanitizeOperations) -> crate::SanitizeAssignment {
        sanitize_topology(
            topology,
            &SanitizeParams {
                operations,
                ..SanitizeParams::default()
            },
        )
        .unwrap()
    }

    #[test]
    fn tautomer_partial_sanitize_transports_final_non_strict_cache_and_num_arom() {
        let topology = cycle(6, BondStereo::None, true);
        let before = topology.clone();
        let assignment = run(
            &topology,
            SanitizeOperations::KEKULIZE
                | SanitizeOperations::SET_AROMATICITY
                | SanitizeOperations::SET_CONJUGATION
                | SanitizeOperations::SET_HYBRIDIZATION
                | SanitizeOperations::ADJUST_HS,
        );
        assert!(assignment.final_valence.is_none());
        let valence = assignment
            .non_strict_valence
            .as_ref()
            .expect("source partial cache");
        assert_eq!(valence.implicit_hydrogens, vec![1; 6]);
        assert_eq!(assignment.aromatic_ring_count, Some(1));
        assert!(
            assignment
                .topology
                .bonds
                .iter()
                .all(|bond| bond.is_aromatic())
        );
        assert_eq!(topology, before);
        let strict = run(&topology, SanitizeOperations::ALL);
        assert!(strict.final_valence.is_some());
        assert!(strict.non_strict_valence.is_none());
        assert_eq!(strict.aromatic_ring_count, Some(1));
        let none = run(&topology, SanitizeOperations::NONE);
        assert!(none.final_valence.is_none());
        assert!(none.non_strict_valence.is_some());
        assert_eq!(none.aromatic_ring_count, None);
    }

    #[test]
    fn sanitize_ring_l07_composed_no_tag_cc_returns_none_despite_cold_hybridizations() {
        let cc = cc(BondStereo::None);
        let snapshot = cc.clone();
        let helper_before = super::hybridizations_probe::helper_calls();
        let find_before = atrop_probe::find_calls();
        let assignment = run(&cc, SanitizeOperations::CLEANUP_ATROPISOMERS);
        // Hybridizations are ALWAYS constructed; NO ring state is acquired.
        assert!(
            super::hybridizations_probe::helper_calls() - helper_before >= 1,
            "hybridizations not constructed"
        );
        assert_eq!(atrop_probe::find_calls() - find_before, 0, "outer find");
        assert!(assignment.final_rings.is_none(), "final rings not None");
        assert_eq!(assignment.topology, snapshot, "topology changed");
    }

    #[test]
    fn sanitize_ring_l07_composed_tagged_cc_acquires_initialized_empty_sssr() {
        let cc = cc(BondStereo::AtropCw);
        let find_before = atrop_probe::find_calls();
        let assignment = run(&cc, SanitizeOperations::CLEANUP_ATROPISOMERS);
        assert_eq!(atrop_probe::find_calls() - find_before, 1, "acquire");
        assert_eq!(assignment.topology.bonds[0].stereo(), BondStereo::None);
        let rings = assignment.final_rings.as_ref().unwrap();
        assert_eq!(rings.find_type(), RingFindType::Sssr);
        assert!(rings.is_initialized());
        assert!(rings.atom_rings().is_empty(), "SSSR not empty on CC");
    }

    #[test]
    fn sanitize_ring_l07_composed_tagged_cycles_with_group_sentinels() {
        for (n, clears) in [(6usize, true), (8usize, false)] {
            let topology = cycle(n, BondStereo::AtropCcw, false);
            let groups_snapshot = topology.stereo_groups.clone();
            let find_before = atrop_probe::find_calls();
            let group_before = atrop_probe::group_cleanup_calls();
            let assignment = run(&topology, SanitizeOperations::CLEANUP_ATROPISOMERS);
            assert_eq!(atrop_probe::find_calls() - find_before, 1, "c{n}: acquire");
            let rings = assignment.final_rings.as_ref().unwrap();
            assert_eq!(rings.find_type(), RingFindType::Sssr, "c{n}");
            let atom_row: Vec<usize> = rings.atom_rings()[0]
                .iter()
                .map(|atom| atom.index())
                .collect();
            // Frozen row convention: atoms start at 0 and run backwards.
            let expected_atoms: Vec<usize> = std::iter::once(0).chain((1..n).rev()).collect();
            assert_eq!(atom_row, expected_atoms, "c{n}: row");
            if clears {
                assert_eq!(
                    assignment.topology.bonds[0].stereo(),
                    BondStereo::None,
                    "c{n}"
                );
                assert_eq!(
                    atrop_probe::group_cleanup_calls() - group_before,
                    1,
                    "c{n}: group cleanup"
                );
            } else {
                assert_eq!(
                    assignment.topology.bonds[0].stereo(),
                    BondStereo::AtropCcw,
                    "c{n}"
                );
                assert_eq!(
                    atrop_probe::group_cleanup_calls() - group_before,
                    0,
                    "c{n}: group cleanup"
                );
            }
            assert_eq!(
                assignment.topology.stereo_groups, groups_snapshot,
                "c{n}: ordered group output"
            );
        }
    }

    #[test]
    fn sanitize_ring_l07_composed_symm_rings_preserved_through_cleanup() {
        let benzene = cycle(6, BondStereo::None, true);
        let find_before = atrop_probe::find_calls();
        let assignment = run(
            &benzene,
            SanitizeOperations::SYMM_RINGS | SanitizeOperations::CLEANUP_ATROPISOMERS,
        );
        // The T stage borrows the SYMM-owned state: no wrapper acquisition.
        assert_eq!(atrop_probe::find_calls() - find_before, 0, "T find");
        let rings = assignment.final_rings.as_ref().unwrap();
        assert_eq!(rings.find_type(), RingFindType::SymmSssr);
        assert_eq!(rings.atom_rings().len(), 1, "benzene Symm row");
    }

    #[test]
    fn sanitize_ring_l07_composed_error_chain_reachability_record() {
        // Recording a structural fact, not a weakened assertion: through the
        // composed public T stage, no VALID topology reaches the private
        // error paths — the copied-branch property cache runs nonstrict
        // (overvalence returns the source -1 sentinel without failing),
        // find_sssr succeeds on valid topologies, and the wrapper's
        // hybridization-length prevalidation cannot trigger because the
        // helper always sizes the assignment to the topology. An invalid
        // topology is rejected at the sanitize entry (stage None), before
        // the T stage. The Rings/Algorithm mapping itself is exercised by
        // the l05 vocabulary controls and the public-entry l04 errors; both
        // remain unconditional. The overvalent first-unspecified input
        // below therefore SUCCEEDS through the copied branch (fields Sp3):
        let mut overvalent_atoms = Vec::new();
        overvalent_atoms.push(Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)));
        overvalent_atoms.push(Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)));
        for id in 2..7usize {
            overvalent_atoms.push(Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::H)));
        }
        let overvalent_bonds = (1..7usize)
            .map(|id| {
                Bond::from_spec(
                    BondId::new(id - 1),
                    BondSpec::new(AtomId::new(0), AtomId::new(id), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        let overvalent = TopologyBlock::try_from_parts(
            overvalent_atoms,
            overvalent_bonds,
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let copy_before = super::hybridizations_probe::copy_calls();
        let overvalent_snapshot = overvalent.clone();
        let assignment = run(&overvalent, SanitizeOperations::CLEANUP_ATROPISOMERS);
        assert_eq!(
            super::hybridizations_probe::copy_calls() - copy_before,
            1,
            "copied branch ran on the overvalent sentinel input"
        );
        // No tags: the outer topology is preserved exactly; the copied
        // branch's computed fields never leak.
        assert_eq!(assignment.topology, overvalent_snapshot, "outer mutated");
        assert!(assignment.final_rings.is_none(), "no outer acquisition");
    }
}

// L09: the eight-call nonempty buffer-transport product. The cfg(test)
// probe above captures the owning carrier's buffer pointer and complete
// state at the ACTUAL final_rings move point, before any consumer; the
// returned atom/bond outer pointers and complete ordered state must equal
// that captured state with no extra acquisition. This is move/preservation
// evidence, not independent finder-order correctness; empty pointers never
// prove allocation, so all eight cells use NONEMPTY states.
#[cfg(test)]
mod final_rings_transport_tests {
    use super::final_rings_probe;
    use crate::atropisomer::cleanup_ring_state_probe as atrop_probe;
    use crate::{RingFindType, SanitizeOperations, SanitizeParams, sanitize_topology};
    use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock};
    use cosmolkit_types::{BondOrder, Element};

    fn benzene() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..6)
                .map(|id| {
                    Atom::from_spec(
                        AtomId::new(id),
                        AtomSpec::new(Element::C).with_aromatic(true),
                    )
                })
                .collect(),
            (0..6)
                .map(|id| {
                    Bond::from_spec(
                        BondId::new(id),
                        BondSpec::new(
                            AtomId::new(id),
                            AtomId::new((id + 1) % 6),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true),
                    )
                })
                .collect(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn cube() -> TopologyBlock {
        let atoms = (0..8)
            .map(|id| Atom::from_spec(AtomId::new(id), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
        ];
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(id, (begin, end))| {
                Bond::from_spec(
                    BondId::new(id),
                    BondSpec::new(AtomId::new(*begin), AtomId::new(*end), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap()
    }

    #[test]
    fn sanitize_ring_l09_final_rings_buffer_transport_eight_calls() {
        let masks = [
            ("S", SanitizeOperations::SYMM_RINGS),
            (
                "S|K",
                SanitizeOperations::SYMM_RINGS | SanitizeOperations::KEKULIZE,
            ),
            (
                "S|A",
                SanitizeOperations::SYMM_RINGS | SanitizeOperations::SET_AROMATICITY,
            ),
            (
                "S|K|A|T",
                SanitizeOperations::SYMM_RINGS
                    | SanitizeOperations::KEKULIZE
                    | SanitizeOperations::SET_AROMATICITY
                    | SanitizeOperations::CLEANUP_ATROPISOMERS,
            ),
        ];
        let mut calls = 0usize;
        for (graph_name, graph, expected_rows) in
            [("benzene", benzene(), 1usize), ("cube", cube(), 6usize)]
        {
            for (mask_name, operations) in masks {
                let graph_snapshot = graph.clone();
                let symm_before = final_rings_probe::symm_acquisitions();
                let find_before = atrop_probe::find_calls();
                let assignment = sanitize_topology(
                    &graph,
                    &SanitizeParams {
                        operations,
                        ..SanitizeParams::default()
                    },
                )
                .unwrap();
                calls += 1;
                // Input preserved immediately after this call.
                assert_eq!(
                    graph, graph_snapshot,
                    "{graph_name}/{mask_name}: input mutated"
                );
                // The ACTUAL S acquisition entry ran exactly once for this
                // call (non-reset counter, relative baseline). The nested
                // selected sanitize for T carries no S bit, so it never
                // bumps this counter or the early-S observation.
                assert_eq!(
                    final_rings_probe::symm_acquisitions() - symm_before,
                    1,
                    "{graph_name}/{mask_name}: S acquisition count"
                );
                // No tags on these graphs: the T stage never acquires, so
                // the SYMM-owned state is the transported state with no
                // extra acquisition.
                assert_eq!(
                    atrop_probe::find_calls() - find_before,
                    0,
                    "{graph_name}/{mask_name}: extra acquisition"
                );
                let rings = assignment
                    .final_rings
                    .as_ref()
                    .unwrap_or_else(|| panic!("{graph_name}/{mask_name}: no rings"));
                assert!(
                    !rings.atom_rings().is_empty(),
                    "{graph_name}/{mask_name}: empty state cannot prove transport"
                );
                assert!(
                    !rings.bond_rings().is_empty(),
                    "{graph_name}/{mask_name}: empty bonds cannot prove transport"
                );
                assert_eq!(
                    rings.find_type(),
                    RingFindType::SymmSssr,
                    "{graph_name}/{mask_name}"
                );
                assert_eq!(
                    rings.atom_rings().len(),
                    expected_rows,
                    "{graph_name}/{mask_name}"
                );
                // BOTH returned outer pointers and the complete ordered
                // state equal the EARLY S-stage observation (survival across
                // K/A/T) AND the final move-point observation.
                let (early_atom_ptr, early_bond_ptr) = final_rings_probe::early_s_pointers();
                assert_eq!(
                    rings.atom_rings().as_ptr() as usize,
                    early_atom_ptr,
                    "{graph_name}/{mask_name}: atom buffer differs from early S stage"
                );
                assert_eq!(
                    rings.bond_rings().as_ptr() as usize,
                    early_bond_ptr,
                    "{graph_name}/{mask_name}: bond buffer differs from early S stage"
                );
                assert_eq!(
                    Some(rings.clone()),
                    final_rings_probe::early_s_state(),
                    "{graph_name}/{mask_name}: ordered state differs from early S stage"
                );
                // The final move-point observation agrees.
                assert_eq!(
                    rings.atom_rings().as_ptr() as usize,
                    final_rings_probe::last_pointer(),
                    "{graph_name}/{mask_name}: buffer not the moved carrier"
                );
                assert_eq!(
                    Some(rings.clone()),
                    final_rings_probe::last_state(),
                    "{graph_name}/{mask_name}: state differs from move-point capture"
                );
                assert!(rings.is_symm_sssr(), "{graph_name}/{mask_name}: quality");
            }
        }
        assert_eq!(calls, 8, "exact census");
    }
}

#[cfg(test)]
#[path = "tests/property_cache.rs"]
mod property_cache_tests;
#[cfg(test)]
mod recovery_chem17_detection {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondSpec, SourceAtomValenceFacts};
    use cosmolkit_types::{BondOrder, Element};
    #[test]
    fn properties_prefix_caches_explicit_even_when_implicit_fails_then_continues() {
        let mut t = TopologyBlock::try_from_parts(
            vec![
                Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C)
                        .with_explicit_hydrogens(4)
                        .with_radical_electrons(1),
                ),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            ],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        for a in &mut t.atoms {
            a.set_source_valence_facts(SourceAtomValenceFacts::UNINITIALIZED);
        }
        let mut report = ChemistryProblemReport::default();
        let v = detection_property_cache(&mut t, true, &mut report).unwrap();
        assert_eq!(report.problems.len(), 1);
        assert!(matches!(
            report.problems[0].error,
            ChemistryProblemError::Valence(ValenceError::InvalidValence {
                phase: crate::ValencePhase::Implicit,
                ..
            })
        ));
        assert_eq!(v.explicit_valence, vec![4, 0]);
        assert_eq!(v.implicit_hydrogens, vec![-1, 4]);
        assert_eq!(t.atoms[0].source_valence_facts().explicit_valence, 4);
        assert_eq!(t.atoms[0].source_valence_facts().implicit_valence, -1);
    }
    #[test]
    fn detection_prepares_actual_cache_before_both_failed_kekulize_attempts() {
        let mut t = TopologyBlock::try_from_parts(
            (0..5)
                .map(|i| {
                    Atom::from_spec(
                        AtomId::new(i),
                        AtomSpec::new(Element::C).with_aromatic(true),
                    )
                })
                .collect(),
            (0..5)
                .map(|i| {
                    Bond::from_spec(
                        cosmolkit_model::BondId::new(i),
                        BondSpec::new(
                            AtomId::new(i),
                            AtomId::new((i + 1) % 5),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap();
        let before = t.clone();
        let mut report = ChemistryProblemReport::default();
        let mut v = detection_property_cache(&mut t, true, &mut report).unwrap();
        assert!(report.problems.is_empty());
        assert_eq!(v.explicit_valence, vec![3; 5]);
        assert_eq!(v.implicit_hydrogens, vec![1; 5]);
        let mut rings = None;
        assert!(matches!(
            crate::kekulize::source_kekulize_for_sanitize(&mut t, &mut v, &mut rings, None),
            Err(KekulizeError::NotKekulizable { .. })
        ));
        assert_eq!(v.explicit_valence, vec![3; 5]);
        assert_eq!(v.implicit_hydrogens, vec![1; 5]);
        assert!(t.bonds.iter().all(|b| b.order() == BondOrder::Single));
        let report = detect_chemistry_problems(&before, &SanitizeParams::default()).unwrap();
        assert_eq!(report.problems.len(), 1);
        assert_eq!(report.problems[0].operation, SanitizeStage::Kekulize);
    }
}

/// The source tautomer-product sequence over borrowed detached state.
/// Ring topology is unchanged; this is neither a live-cache capability nor
/// sanitizeMol's retry controller. Completed stores remain observable on Err.
#[doc(hidden)]
pub fn source_sanitize_tautomer_product(
    topology: &mut TopologyBlock,
    properties: &mut cosmolkit_model::MoleculeProperties,
    valence: &mut ValenceAssignment,
    rings: &mut RingInfo,
) -> Result<(), SanitizeError> {
    // BEGIN RECOVERY SEARCH04 SOURCE TautomerProductStages
    // RDKit❗❌:             // We only change bond orders/H counts/charges; the molecular graph
    // RDKit❗❌:             // (and therefore ring topology) is unchanged.
    // RDKit❗❌:             // `sanitizeMol()` always calls `clearComputedProps()` which resets
    // RDKit❗❌:             // ring info and forces ring-finding for each generated tautomer.
    // RDKit❗❌:             // Avoid that by clearing computed props without touching rings,
    // RDKit❗❌:             // then running the specific sanitize steps we need.
    // RDKit❗❌:             product->clearComputedProps(false);
    // RDKit❗❌:             product->updatePropertyCache(false);
    // RDKit❗❌:             MolOps::Kekulize(*product);
    // RDKit❗❌:             MolOps::setAromaticity(*product);
    // RDKit❗❌:             MolOps::setConjugation(*product);
    // RDKit❗❌:             MolOps::setHybridization(*product);
    // RDKit❗❌:             MolOps::adjustHs(*product);
    // END RECOVERY SEARCH04 SOURCE TautomerProductStages
    // BEGIN RECOVERY SEARCH04 SOURCE ROMolClearComputedProps
    // RDKit❗❌: void ROMol::clearComputedProps(bool includeRings) const {
    // RDKit❗❌:   // the SSSR information:
    // RDKit❗❌:   if (includeRings) {
    // RDKit❗❌:     this->dp_ringInfo->reset();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   RDProps::clearComputedProps();
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atom : atoms()) {
    // RDKit❗❌:     atom->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   for (auto bond : bonds()) {
    // RDKit❗❌:     bond->clearComputedProps();
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RECOVERY SEARCH04 SOURCE ROMolClearComputedProps
    // Source includeRings=false deliberately preserves the same carrier.
    // Existing detached aromaticity/conjugation/adjustHs owners allocate owned
    // topology results; their additional cloning/validation costs remain ❌.
    // This narrow orchestrator adds no second chemistry algorithm or retry.
    properties.clear_computed_props()?;
    clear_topology_computed_properties(topology)?;
    *valence = assign_property_cache(topology, &PropertyCacheParams { strict: false })
        .map_err(|source| SanitizeError::Properties {
            stage: SanitizeStage::Properties,
            source,
        })?
        .into_valence();
    crate::source_kekulize_attempt(
        topology,
        valence,
        rings,
        &KekulizeParams {
            mark_atoms_bonds: true,
            canonical: true,
            max_backtracks: 100,
        },
    )
    .map_err(|source| SanitizeError::Kekulize {
        stage: SanitizeStage::Kekulize,
        source,
    })?;
    if !rings.is_initialized() {
        *rings = symmetrized_sssr(topology, &RingSearchParams::default()).map_err(|source| {
            SanitizeError::Rings {
                stage: SanitizeStage::SetAromaticity,
                source,
            }
        })?;
    }
    let aromaticity = assign_default_aromaticity_with_cached_valence(
        topology, rings, valence, None,
    )
    .map_err(|source| SanitizeError::Aromaticity {
        stage: SanitizeStage::SetAromaticity,
        source,
    })?;
    *topology = aromaticity.topology;
    properties.set_computed_prop(
        "numArom",
        i32::try_from(aromaticity.aromatic_ring_count).map_err(|_| SanitizeError::Aromaticity {
            stage: SanitizeStage::SetAromaticity,
            source: AromaticityError::IntegerOverflow {
                field: "source numArom int",
            },
        })?,
    )?;
    *topology =
        assign_conjugation(topology, valence).map_err(|source| SanitizeError::Conjugation {
            stage: SanitizeStage::SetConjugation,
            source,
        })?;
    let hybridization =
        assign_hybridization(topology, valence).map_err(|source| SanitizeError::Hybridization {
            stage: SanitizeStage::SetHybridization,
            source,
        })?;
    materialize_hybridization(topology, &hybridization);
    let adjusted = adjust_hs(topology, valence).map_err(|source| SanitizeError::AdjustHs {
        stage: SanitizeStage::AdjustHs,
        source,
    })?;
    *topology = adjusted.topology;
    *valence = adjusted.valence;
    Ok(())
}

#[cfg(test)]
mod search04_product_stage_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondId, BondSpec, MoleculeProperties};
    use cosmolkit_types::{BondOrder, Element};
    fn cycle(n: usize) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..n)
                .map(|i| {
                    Atom::from_spec(
                        AtomId::new(i),
                        AtomSpec::new(Element::C).with_aromatic(true),
                    )
                })
                .collect(),
            (0..n)
                .map(|i| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(
                            AtomId::new(i),
                            AtomId::new((i + 1) % n),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn search04_product_stages_preserve_initialized_ring_carrier_and_clear_only_computed_props() {
        let mut topology = cycle(6);
        topology.atoms[0]
            .set_computed_prop("drop-atom", 1_i32)
            .unwrap();
        topology.atoms[0].set_prop("keep-atom", 2_i32).unwrap();
        topology.bonds[0]
            .set_computed_prop("drop-bond", 3_i32)
            .unwrap();
        topology.bonds[0].set_prop("keep-bond", 4_i32).unwrap();
        let mut properties = MoleculeProperties::default();
        properties
            .set_computed_prop("drop-molecule", 5_i32)
            .unwrap();
        properties.set_prop("keep-molecule", 6_i32).unwrap();
        let mut rings = crate::fast_find_rings(&topology).unwrap();
        let rings_before = rings.clone();
        let pointer = rings.atom_rings().as_ptr();
        let mut valence = ValenceAssignment {
            explicit_valence: vec![],
            implicit_hydrogens: vec![],
        };
        source_sanitize_tautomer_product(&mut topology, &mut properties, &mut valence, &mut rings)
            .unwrap();
        assert_eq!(rings, rings_before);
        assert_eq!(rings.atom_rings().as_ptr(), pointer);
        assert!(rings.is_find_fast_or_better());
        assert!(!rings.is_symm_sssr());
        assert!(topology.atoms.iter().all(|a| a.is_aromatic()));
        assert!(
            topology
                .bonds
                .iter()
                .all(|b| b.order() == BondOrder::Aromatic && b.is_conjugated())
        );
        assert!(
            topology
                .atoms
                .iter()
                .all(|a| a.hybridization() == Hybridization::Sp2)
        );
        assert_eq!(valence.explicit_valence, vec![3; 6]);
        assert_eq!(valence.implicit_hydrogens, vec![1; 6]);
        assert_eq!(topology.atoms[0].prop("drop-atom"), None);
        assert!(topology.atoms[0].prop("keep-atom").is_some());
        assert_eq!(topology.bonds[0].prop("drop-bond"), None);
        assert!(topology.bonds[0].prop("keep-bond").is_some());
        assert_eq!(properties.prop("drop-molecule"), None);
        assert!(properties.prop("keep-molecule").is_some());
        assert_eq!(
            properties.prop("numArom"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
    }
    #[test]
    fn search04_product_stages_use_one_attempt_and_retain_completed_failure_state() {
        let mut topology = cycle(5);
        let mut rings = crate::symmetrized_sssr(&topology, &Default::default()).unwrap();
        let before = rings.clone();
        let mut valence = ValenceAssignment {
            explicit_valence: vec![],
            implicit_hydrogens: vec![],
        };
        let mut properties = MoleculeProperties::default();
        properties.set_computed_prop("drop", 1_i32).unwrap();
        let error = source_sanitize_tautomer_product(
            &mut topology,
            &mut properties,
            &mut valence,
            &mut rings,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            SanitizeError::Kekulize {
                stage: SanitizeStage::Kekulize,
                source: KekulizeError::NotKekulizable { .. }
            }
        ));
        assert_eq!(rings, before);
        assert_eq!(properties.prop("drop"), None);
        assert_eq!(properties.prop("numArom"), None);
        assert_eq!(valence.explicit_valence.len(), 5);
        assert_eq!(valence.implicit_hydrogens.len(), 5);
        assert!(topology.atoms.iter().all(|a| a.is_aromatic()));
        assert!(
            topology
                .bonds
                .iter()
                .all(|b| b.order() == BondOrder::Single && b.is_aromatic())
        );
    }
}
