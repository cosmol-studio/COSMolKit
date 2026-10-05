//! Detached valence result and error vocabulary shared by source-backed owners.
//! Calculation and cache preparation belong to algorithm crates.

use crate::{AtomId, BondId, TopologyValidationError};
use cosmolkit_types::BondOrder;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ValencePhase {
    EffectiveAtomicNumber,
    Explicit,
    Implicit,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum ValenceError {
    /// Source explicit getter PRECONDITION used by numPi; preserve literal text.
    #[error("getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()")]
    PiElectronExplicitValenceCacheNotInitialized { atom: AtomId },
    /// Source CHECK_INVARIANT in numPiElectrons; retain its literal message.
    #[error("explicit valence exceeds atom degree")]
    PiElectronInvariant {
        atom: AtomId,
        explicit_valence: u32,
        physical_bonds: u32,
    },
    #[error("{message}")]
    InvalidValence {
        atom: AtomId,
        atomic_number: u8,
        formal_charge: i8,
        phase: ValencePhase,
        calculated: Option<i32>,
        reason: &'static str,
        message: String,
    },
    #[error("invalid topology: {source}")]
    InvalidTopology { source: TopologyValidationError },
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error(
        "adjacency for atom {atom} references neighbor atom row {neighbor_atom}, out of range for {atom_count} atoms"
    )]
    AdjacencyAtomOutOfRange {
        atom: AtomId,
        neighbor_atom: usize,
        atom_count: usize,
    },
    #[error(
        "adjacency for atom {atom} references bond {bond}, out of range for {bond_count} bonds"
    )]
    AdjacencyBondOutOfRange {
        atom: AtomId,
        bond: BondId,
        bond_count: usize,
    },
    #[error(
        "adjacency for atom {atom} references neighbor row {neighbor_atom} through bond {bond}, but that bond has endpoints {begin}-{end}"
    )]
    AdjacencyEndpointMismatch {
        atom: AtomId,
        neighbor_atom: usize,
        bond: BondId,
        begin: AtomId,
        end: AtomId,
    },
    #[error("explicit valence input for atom {atom} must be nonnegative, got {value}")]
    InvalidExplicitValenceInput { atom: AtomId, value: i32 },
    #[error("periodic-table field {field} is unavailable for atomic number {atomic_number}")]
    PeriodicTableLookup {
        atomic_number: u8,
        field: &'static str,
    },
    #[error("explicit valence is not available for atom {atom}")]
    ExplicitValenceCacheNotInitialized { atom: AtomId },
    #[error("implicit valence is not available for atom {atom}")]
    ImplicitValenceCacheNotInitialized { atom: AtomId },
    #[error(
        "hydrogen count overflow at atom {atom}: explicit={explicit}, implicit={implicit}, neighbor_hydrogens={neighbor_hydrogens}"
    )]
    HydrogenCountOverflow {
        atom: AtomId,
        explicit: u32,
        implicit: u32,
        neighbor_hydrogens: usize,
    },
    #[error("Bad bond type")]
    BadBondType {
        bond: Option<BondId>,
        order: BondOrder,
    },
}

/// Owned context-dependent atom read results. Canonical Atom vocabulary stays
/// in the model; these rows hold only degree and calculated valence metadata.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomMetadata {
    pub degree: usize,
    pub explicit_valence: i32,
    pub implicit_hydrogens: i32,
    pub total_hydrogens: i32,
    pub total_valence: i32,
}
