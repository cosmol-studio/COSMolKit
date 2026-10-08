// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::collections::{BTreeMap, BTreeSet, VecDeque};

use cosmolkit_model::{
    Atom, AtomId, Bond, BondId, BondValueError, PropertyValue, PropertyValueKind, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

use crate::hcount::total_hydrogen_count_from_validated;
use crate::{
    DoubleBondControl, DoubleBondStereoDescriptor, DoubleBondStereoError,
    DoubleBondStereoSpecified, RingInfo, StereoOrderError, ValenceAssignment, ValenceError,
    bond_affects_atom_chirality, count_swaps_to_interconvert, double_bond_stereo_info,
    is_atom_bridgehead_from_topology,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PotentialStereoParams {
    pub clean: bool,
    pub flag_possible: bool,
    pub allow_nontetrahedral: bool,
}

impl Default for PotentialStereoParams {
    fn default() -> Self {
        Self {
            clean: false,
            flag_possible: true,
            allow_nontetrahedral: true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoType {
    AtomTetrahedral,
    AtomSquarePlanar,
    AtomTrigonalBipyramidal,
    AtomOctahedral,
    BondDouble,
    BondCumuleneEven,
    BondAtropisomer,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoSpecified {
    Unspecified,
    Specified,
    Unknown,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoDescriptor {
    None,
    TetrahedralClockwise,
    TetrahedralCounterclockwise,
    BondCis,
    BondTrans,
    BondAtropCw,
    BondAtropCcw,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PotentialStereoCenter {
    Atom(AtomId),
    Bond(BondId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PotentialStereoInfo {
    pub stereo_type: PotentialStereoType,
    pub specified: PotentialStereoSpecified,
    pub centered_on: PotentialStereoCenter,
    pub descriptor: PotentialStereoDescriptor,
    pub permutation: u32,
    pub controlling_atoms: Vec<Option<AtomId>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RingStereoRelation {
    pub atom: AtomId,
    pub other: AtomId,
    pub same_orientation: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PotentialStereoAssignment {
    pub stereo: Vec<PotentialStereoInfo>,
    pub atom_ranks: Vec<u32>,
    pub ring_relations: Vec<RingStereoRelation>,
    pub cleaned_topology: Option<TopologyBlock>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PotentialStereoError {
    #[error("atom {atom} property {property} has invalid kind {kind:?}")]
    InvalidPropertyKind {
        atom: AtomId,
        property: &'static str,
        kind: PropertyValueKind,
    },
    #[error("atom {atom} has invalid ring stereo reference {value}")]
    InvalidRingStereoReference { atom: AtomId, value: i32 },
    #[error("atom {atom} has no other ring atoms")]
    EmptyRingStereoReferences { atom: AtomId },
    #[error("ring preparation failed: {0}")]
    Ring(#[from] crate::RingFindingError),
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("valence field {field} has {actual} rows, expected {atom_count}")]
    InvalidValence {
        field: &'static str,
        actual: usize,
        atom_count: usize,
    },
    #[error("valence field {field} has invalid value {value} for atom {atom}")]
    InvalidValenceValue {
        field: &'static str,
        atom: AtomId,
        value: i32,
    },
    #[error("atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("invalid ring information: {reason} at row {row} (value {value}, limit {limit})")]
    InvalidRingInfo {
        reason: &'static str,
        row: usize,
        value: usize,
        limit: usize,
    },
    #[error("atom {atom} has invalid nonzero degree {degree}")]
    InvalidAtomDegree { atom: AtomId, degree: usize },
    #[error("bond {bond} {endpoint} endpoint has invalid degree {degree}")]
    InvalidBondDegree {
        bond: BondId,
        endpoint: &'static str,
        degree: usize,
    },
    #[error("bond {bond} has invalid stereo references: {reason}")]
    InvalidStereoReferences { bond: BondId, reason: &'static str },
    #[error("bond {bond} has unsupported order {order:?}")]
    UnsupportedBondOrder { bond: BondId, order: BondOrder },
    #[error("atom {atom} has invalid chiral permutation {permutation}")]
    InvalidChiralPermutation { atom: AtomId, permutation: u32 },
    #[error("atropisomer support is unavailable for bond {bond}")]
    AtropisomerDependencyUnavailable { bond: BondId },
    #[error("potential-stereo refinement did not converge after {iterations} iterations")]
    RefinementDidNotConverge { iterations: usize },
    #[error("valence calculation failed: {0}")]
    Valence(#[from] ValenceError),
    #[error("tetrahedral order calculation failed: {0}")]
    StereoOrder(#[from] StereoOrderError),
    #[error("double-bond stereo calculation failed: {0}")]
    DoubleStereo(#[from] DoubleBondStereoError),
    #[error("invalid detached bond value: {0}")]
    BondValue(#[from] BondValueError),
}

fn validate_topology_and_valence(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
) -> Result<(), PotentialStereoError> {
    topology.validate()?;
    let atom_count = topology.atoms.len();
    for (field, rows) in [
        ("explicit_valence", &valence.explicit_valence),
        ("implicit_hydrogens", &valence.implicit_hydrogens),
    ] {
        if rows.len() != atom_count {
            return Err(PotentialStereoError::InvalidValence {
                field,
                actual: rows.len(),
                atom_count,
            });
        }
        if let Some((index, value)) = rows.iter().copied().enumerate().find(|(index, _)| {
            // RDKit✔️✔️:   return df_noImplicit ? 0 : d_implicitValence;
            // Source initialization/getter semantics belong to the unique
            // valence and hcount owners. Unused NoImplicit I is retained.
            // Two O(1) indexed getters preserve the existing O(V) validation.
            let atom = &topology.atoms[*index];
            if field == "explicit_valence" {
                crate::valence::cached_explicit_valence(atom, Some(valence)).is_err()
            } else {
                crate::hcount::implicit_hydrogen_count(atom, valence).is_err()
            }
        }) {
            return Err(PotentialStereoError::InvalidValenceValue {
                field,
                atom: AtomId::new(index),
                value,
            });
        }
    }
    Ok(())
}

fn validate_ring_rows(
    topology: &TopologyBlock,
    rings: &RingInfo,
) -> Result<(), PotentialStereoError> {
    if !rings.is_initialized() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "ring information is not initialized",
            row: 0,
            value: 0,
            limit: 0,
        });
    }
    let atom_count = topology.atoms.len();
    let preserved_prefix = (rings.atom_row_count() != atom_count
        || rings.bond_row_count() != topology.bonds.len())
        && crate::rings::preserves_appended_terminal_hydrogen_ring_prefix(topology, rings);
    if rings.atom_row_count() != atom_count && !preserved_prefix {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "atom membership row count mismatch",
            row: 0,
            value: rings.atom_row_count(),
            limit: atom_count,
        });
    }
    if rings.bond_row_count() != topology.bonds.len() && !preserved_prefix {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "bond membership row count mismatch",
            row: 0,
            value: rings.bond_row_count(),
            limit: topology.bonds.len(),
        });
    }
    if rings.atom_rings().len() != rings.bond_rings().len() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "atom/bond ring table length mismatch",
            row: 0,
            value: rings.atom_rings().len(),
            limit: rings.bond_rings().len(),
        });
    }
    for (row, (atoms, bonds)) in rings
        .atom_rings()
        .iter()
        .zip(rings.bond_rings())
        .enumerate()
    {
        if atoms.len() != bonds.len() {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring atom/bond size mismatch",
                row,
                value: atoms.len(),
                limit: bonds.len(),
            });
        }
        if let Some(atom) = atoms.iter().find(|atom| atom.index() >= atom_count) {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring atom out of range",
                row,
                value: atom.index(),
                limit: atom_count,
            });
        }
        if let Some(bond) = bonds
            .iter()
            .find(|bond| bond.index() >= topology.bonds.len())
        {
            return Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring bond out of range",
                row,
                value: bond.index(),
                limit: topology.bonds.len(),
            });
        }
    }
    Ok(())
}

fn validate_inputs(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) -> Result<(), PotentialStereoError> {
    validate_topology_and_valence(topology, valence)?;
    if !rings.is_symm_sssr() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "ring information is not symmetric SSSR",
            row: 0,
            value: rings.find_type() as usize,
            limit: 0,
        });
    }
    validate_ring_rows(topology, rings)?;
    Ok(())
}

fn graph_degree(topology: &TopologyBlock, atom: AtomId) -> usize {
    topology.adjacency.neighbors_of(atom.index()).len()
}

fn nonzero_degree(topology: &TopologyBlock, atom: AtomId) -> Result<usize, PotentialStereoError> {
    let mut degree = 0;
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let bond = &topology.bonds[neighbor.bond.index()];
        degree += usize::from(bond_affects_atom_chirality(bond, atom)?);
    }
    Ok(degree)
}

fn total_hydrogens(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
    include_neighbors: bool,
) -> Result<usize, PotentialStereoError> {
    // Every call is below shared topology/valence validation in either the
    // whole-analysis path or the selected-atom writer batch.
    Ok(total_hydrogen_count_from_validated(topology, valence, atom, include_neighbors)? as usize)
}

pub(crate) fn total_degree(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
) -> Result<usize, PotentialStereoError> {
    // RDKit❗✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit❗✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION Atom::getTotalDegree
    // The source GetTotalNumHs read includes NoImplicit and signed field width.
    // Reuse its unique owner rather than interpreting detached rows directly.
    // O(1): includeNeighbors=false and adjacency length do not scan the graph.
    Ok(graph_degree(topology, atom)
        + total_hydrogen_count_from_validated(topology, valence, atom, false)? as usize)
}

fn has_protium_neighbor(topology: &TopologyBlock, atom: AtomId) -> bool {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| {
            let value = &topology.atoms[neighbor.atom_index];
            value.atomic_number() == 1 && value.isotope().is_none()
        })
}

fn has_conjugated_bond(topology: &TopologyBlock, atom: AtomId) -> bool {
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| topology.bonds[neighbor.bond.index()].is_conjugated())
}

fn is_potential_nontetrahedral(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION isAtomPotentialNontetrahedralCenter
    // RDKit✔️❌: auto nzdegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️❌: auto impHDegree = atom->getTotalNumHs();
    // RDKit✔️❌: auto tnzdegree = nzdegree + impHDegree;
    // RDKit✔️❌: if (tnzdegree > 6 || tnzdegree < 2 || (anum < 12 && anum != 4)) return false;
    // RDKit✔️❌: if (chiralType >= Atom::CHI_SQUAREPLANAR &&
    // RDKit✔️❌:     chiralType <= Atom::CHI_OCTAHEDRAL) return true;
    // RDKit✔️❌: if (chiralType == Atom::CHI_UNSPECIFIED && tnzdegree >= 4) return true;
    // END RDKIT CPP FUNCTION isAtomPotentialNontetrahedralCenter
    // Hydrogen composition reuses the source-backed helper for already
    // validated topology; it does not repeat full-topology validation per atom.
    let value = &topology.atoms[atom.index()];
    let degree = nonzero_degree(topology, atom)? + total_hydrogens(topology, valence, atom, false)?;
    if degree > 6 || degree < 2 || (value.atomic_number() < 12 && value.atomic_number() != 4) {
        return Ok(false);
    }
    Ok(matches!(
        value.chiral_tag(),
        ChiralTag::SquarePlanar | ChiralTag::TrigonalBipyramidal | ChiralTag::Octahedral
    ) || (value.chiral_tag() == ChiralTag::Unspecified && degree >= 4))
}

fn is_potential_tetrahedral(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&RingInfo>,
    ring_rows_validated: &mut bool,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // RDKit❗✔️: bool isAtomPotentialTetrahedralCenter(const Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(atom, "atom is null");
    // RDKit❗✔️:   auto nzDegree = getAtomNonzeroDegree(atom);
    // RDKit❗✔️:   auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit❗✔️:   if (tnzDegree > 4) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     const auto &mol = atom->getOwningMol();
    // RDKit❗✔️:     if (nzDegree == 4) {
    // RDKit❗✔️:       // chirality is always possible with 4 nbrs
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     } else if (nzDegree <= 1) {
    // RDKit❗✔️:       // chirality is never possible with 0 or 1 nbr
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     } else if (nzDegree < 3 &&
    // RDKit❗✔️:                (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
    // RDKit❗✔️:       // less than three neighbors is never stereogenic
    // RDKit❗✔️:       // unless it is a phosphine/arsine with implicit H
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
    // RDKit❗✔️:       // from logical flow: degree is 2 or 3 (implicit H)
    // RDKit❗✔️:       // Since InChI Software v. 1.02-standard (2009), phosphines and arsines
    // RDKit❗✔️:       // are always treated as stereogenic even with H atom neighbors.
    // RDKit❗✔️:       // Accept automatically.
    // RDKit❗✔️:       return true;
    // RDKit❗✔️:     } else if (nzDegree == 3) {
    // RDKit❗✔️:       // three-coordinate with a single H we'll accept automatically:
    // RDKit❗✔️:       if (atom->getTotalNumHs() == 1) {
    // RDKit❗✔️:         if (detail::has_protium_neighbor(mol, atom)) {
    // RDKit❗✔️:           // more than one H is never stereogenic
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         return true;
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // otherwise we default to not being a legal center
    // RDKit❗✔️:         bool legalCenter = false;
    // RDKit❗✔️:         // but there are a few special cases we'll accept
    // RDKit❗✔️:         // sulfur or selenium with either a positive charge or a double
    // RDKit❗✔️:         // bond:
    // RDKit❗✔️:         if ((atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) &&
    // RDKit❗✔️:             (atom->getValence(Atom::ValenceType::EXPLICIT) == 4 ||
    // RDKit❗✔️:              (atom->getValence(Atom::ValenceType::EXPLICIT) == 3 &&
    // RDKit❗✔️:               atom->getFormalCharge() == 1))) {
    // RDKit❗✔️:           legalCenter = true;
    // RDKit❗✔️:         } else if (atom->getAtomicNum() == 7) {
    // RDKit❗✔️:           // three-coordinate N additional requirements:
    // RDKit❗✔️:           //   in a ring of size 3  (from InChI)
    // RDKit❗✔️:           // OR
    // RDKit❗✔️:           /// is a bridgehead atom (RDKit extension)
    // RDKit❗✔️:           // Also: cannot be SP2 hybridized or have a conjugated bond
    // RDKit❗✔️:           //   (this was Github #7434)
    // RDKit❗✔️:           if (atom->getHybridization() == Atom::HybridizationType::SP3 &&
    // RDKit❗✔️:               !MolOps::atomHasConjugatedBond(atom) &&
    // RDKit❗✔️:               (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
    // RDKit❗✔️:                queryIsAtomBridgehead(atom))) {
    // RDKit❗✔️:             legalCenter = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         return legalCenter;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       return false;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Source scalar and checked batch share this predicate. Scalar callers
    // carry actual source cache membership vectors, which may be sparse after
    // reset/addRing; the batch interface retains its full detached validation.
    let value = &topology.atoms[atom.index()];
    let degree = nonzero_degree(topology, atom)?;
    let hydrogens = total_hydrogens(topology, valence, atom, false)?;
    if degree + hydrogens > 4 {
        return Ok(false);
    }
    if degree == 4 {
        return Ok(true);
    }
    if degree <= 1 {
        return Ok(false);
    }
    if degree < 3 && !matches!(value.atomic_number(), 15 | 33) {
        return Ok(false);
    }
    if matches!(value.atomic_number(), 15 | 33) {
        return Ok(true);
    }
    if degree == 3 {
        if hydrogens == 1 {
            return Ok(!has_protium_neighbor(topology, atom));
        }
        if matches!(value.atomic_number(), 16 | 34) {
            // RDKit✔️✔️:         if ((atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) &&
            // RDKit✔️✔️:             (atom->getValence(Atom::ValenceType::EXPLICIT) == 4 ||
            // RDKit✔️✔️:              (atom->getValence(Atom::ValenceType::EXPLICIT) == 3 &&
            // RDKit✔️✔️:               atom->getFormalCharge() == 1))) {
            // Reuse the actual source getter, with O(1) signed-width reads.
            let explicit = crate::valence::cached_explicit_valence(value, Some(valence))?;
            if explicit == 4 || (explicit == 3 && value.formal_charge() == 1) {
                return Ok(true);
            }
        }
        if value.atomic_number() == 7
            && value.hybridization() == Hybridization::Sp3
            && !has_conjugated_bond(topology, atom)
        {
            let rings = rings.ok_or(PotentialStereoError::InvalidRingInfo {
                reason: "ring information is unavailable for nitrogen potential-stereo branch",
                row: 0,
                value: 0,
                limit: 0,
            })?;
            if !rings.is_initialized() {
                return Err(PotentialStereoError::InvalidRingInfo {
                    reason: "ring information is not initialized",
                    row: 0,
                    value: 0,
                    limit: 0,
                });
            }
            if !*ring_rows_validated {
                validate_ring_rows(topology, rings)?;
                *ring_rows_validated = true;
            }
            if rings.is_atom_in_ring_of_size(atom, 3)
                || is_atom_bridgehead_from_topology(topology, atom.index(), rings) != 0
            {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

fn is_potential_atom(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    atom: AtomId,
    allow_nontetrahedral: bool,
) -> Result<bool, PotentialStereoError> {
    // `potential_stereo` validates the complete symmetric ring assignment
    // before this whole-topology analysis begins.
    let mut ring_rows_validated = true;
    Ok(is_potential_tetrahedral(
        topology,
        valence,
        Some(rings),
        &mut ring_rows_validated,
        atom,
    )? || (allow_nontetrahedral && is_potential_nontetrahedral(topology, valence, atom)?))
}

/// Evaluate the source tetrahedral-potential predicate for the requested atom
/// IDs in order, sharing detached-state validation across the batch.
///
/// Ring data is consulted and validated only when an atom reaches the source's
/// ring-dependent nitrogen branch; the supplied ring finding type is retained.
/// Scalar getter over actual source topology, valence cache and ring cache.
/// Source reset/addRing caches permit sparse membership rows. Initialization
/// is checked at the nitrogen ring branch, as in RingInfo's source getters;
/// the separate checked-batch interface retains its complete shape checks.
#[doc(hidden)]
pub fn potential_tetrahedral_center_from_source(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&RingInfo>,
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // RDKit❗✔️: Chirality::detail::isAtomPotentialTetrahedralCenter(atom)
    // The shared predicate below contains the complete pinned helper body.
    // O(degree) source getters, no complete graph scan or per-center cache copy.
    if atom.index() >= topology.atoms.len() {
        return Err(PotentialStereoError::AtomOutOfRange {
            atom,
            atom_count: topology.atoms.len(),
        });
    }
    let mut source_cache_rows = true;
    is_potential_tetrahedral(topology, valence, rings, &mut source_cache_rows, atom)
}

pub fn potential_tetrahedral_centers_for_atoms(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&RingInfo>,
    atoms: &[AtomId],
) -> Result<Vec<bool>, PotentialStereoError> {
    validate_topology_and_valence(topology, valence)?;
    let atom_count = topology.atoms.len();
    let mut ring_rows_validated = false;
    let mut centers = Vec::with_capacity(atoms.len());
    for &atom in atoms {
        if atom.index() >= atom_count {
            return Err(PotentialStereoError::AtomOutOfRange { atom, atom_count });
        }
        centers.push(is_potential_tetrahedral(
            topology,
            valence,
            rings,
            &mut ring_rows_validated,
            atom,
        )?);
    }
    Ok(centers)
}

fn atom_info(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atom: AtomId,
    allow_nontetrahedral: bool,
) -> Result<PotentialStereoInfo, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION getStereoInfo(const Atom *)
    // RDKit✔️✔️: for (const auto &nbri : mol.getAtomBonds(atom)) {
    // RDKit✔️✔️:   if (bnd->getBondDir() == Bond::UNKNOWN) explicitUnknownStereo = 1;
    // RDKit✔️✔️:   sinfo.controllingAtoms.push_back(bnd->getOtherAtomIdx(atom->getIdx()));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: std::vector<unsigned> origNbrOrder = sinfo.controllingAtoms;
    // RDKit✔️✔️: std::sort(sinfo.controllingAtoms.begin(), sinfo.controllingAtoms.end());
    // RDKit✔️✔️: if (explicitUnknownStereo) sinfo.specified = StereoSpecified::Unknown;
    // RDKit✔️✔️: else if (stereo == CHI_TETRAHEDRAL_CCW || stereo == CHI_TETRAHEDRAL_CW) {
    // RDKit✔️✔️:   unsigned nSwaps = countSwapsToInterconvert(origNbrOrder, controllingAtoms);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getStereoInfo(const Atom *)
    let value = &topology.atoms[atom.index()];
    let original = topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .map(|neighbor| AtomId::new(neighbor.atom_index))
        .collect::<Vec<_>>();
    let mut controls = original.clone();
    controls.sort_unstable();
    let explicit_unknown = value.unknown_stereo()
        || topology
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .any(|neighbor| {
                let bond = &topology.bonds[neighbor.bond.index()];
                bond.direction() == BondDirection::Unknown || bond.unknown_stereo()
            });
    let mut info = PotentialStereoInfo {
        stereo_type: PotentialStereoType::AtomTetrahedral,
        specified: PotentialStereoSpecified::Unspecified,
        centered_on: PotentialStereoCenter::Atom(atom),
        descriptor: PotentialStereoDescriptor::None,
        permutation: 0,
        controlling_atoms: controls.iter().copied().map(Some).collect(),
    };
    if explicit_unknown {
        info.specified = PotentialStereoSpecified::Unknown;
        return Ok(info);
    }
    if matches!(
        value.chiral_tag(),
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    ) {
        info.specified = PotentialStereoSpecified::Specified;
        let odd = count_swaps_to_interconvert(&original, &controls)? % 2 == 1;
        let clockwise = (value.chiral_tag() == ChiralTag::TetrahedralCw) ^ odd;
        info.descriptor = if clockwise {
            PotentialStereoDescriptor::TetrahedralClockwise
        } else {
            PotentialStereoDescriptor::TetrahedralCounterclockwise
        };
        return Ok(info);
    }
    if allow_nontetrahedral && is_potential_nontetrahedral(topology, valence, atom)? {
        let degree = total_degree(topology, valence, atom)?;
        info.stereo_type = match value.chiral_tag() {
            ChiralTag::SquarePlanar => PotentialStereoType::AtomSquarePlanar,
            ChiralTag::TrigonalBipyramidal => PotentialStereoType::AtomTrigonalBipyramidal,
            ChiralTag::Octahedral => PotentialStereoType::AtomOctahedral,
            ChiralTag::Unspecified if degree == 5 => PotentialStereoType::AtomTrigonalBipyramidal,
            ChiralTag::Unspecified if degree == 6 => PotentialStereoType::AtomOctahedral,
            _ => PotentialStereoType::AtomTetrahedral,
        };
        if let Some(permutation) = value.chiral_permutation() {
            info.permutation = permutation;
            info.specified = if permutation == 0 {
                PotentialStereoSpecified::Unknown
            } else {
                PotentialStereoSpecified::Specified
            };
        }
    }
    Ok(info)
}

pub(crate) fn is_potential_bond(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    bond: &Bond,
) -> Result<bool, PotentialStereoError> {
    // A standalone detached candidate check owns a fresh structural proof.
    let mut preserved_prefix = None;
    is_potential_bond_with_prefix_check(topology, valence, rings, bond, &mut preserved_prefix)
}

/// Native FileParsers consumes bondRings directly, not checked membership rows.
pub(crate) fn is_potential_bond_source(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    bond: &Bond,
) -> Result<bool, PotentialStereoError> {
    is_potential_bond_impl(topology, valence, rings, bond, None)
}

fn bond_info(
    topology: &TopologyBlock,
    bond: BondId,
) -> Result<PotentialStereoInfo, PotentialStereoError> {
    let state = &topology.bonds[bond.index()];
    if state.order() == BondOrder::Single
        && matches!(state.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw)
    {
        return represented_atropisomer_info(topology, state);
    }
    let source = double_bond_stereo_info(topology, bond)?;
    Ok(PotentialStereoInfo {
        stereo_type: PotentialStereoType::BondDouble,
        specified: match source.specified {
            DoubleBondStereoSpecified::Unspecified => PotentialStereoSpecified::Unspecified,
            DoubleBondStereoSpecified::Specified => PotentialStereoSpecified::Specified,
            DoubleBondStereoSpecified::Unknown => PotentialStereoSpecified::Unknown,
        },
        centered_on: PotentialStereoCenter::Bond(bond),
        descriptor: match source.descriptor {
            Some(DoubleBondStereoDescriptor::Cis) => PotentialStereoDescriptor::BondCis,
            Some(DoubleBondStereoDescriptor::Trans) => PotentialStereoDescriptor::BondTrans,
            None => PotentialStereoDescriptor::None,
        },
        permutation: 0,
        controlling_atoms: source
            .controlling_atoms
            .into_iter()
            .map(|control| match control {
                DoubleBondControl::Atom(atom) => Some(atom),
                DoubleBondControl::Implicit => None,
            })
            .collect(),
    })
}

fn atom_symbol(atom: &Atom) -> String {
    format!(
        "{}{}{}",
        atom.isotope().unwrap_or(0),
        atom.element().symbol(),
        atom.formal_charge()
    )
}

fn bond_symbol(bond: &Bond) -> &'static str {
    // BEGIN RDKIT CPP FUNCTION getBondSymbol
    // RDKit✔️✔️: if (bond->getIsAromatic()) res = ":";
    // RDKit✔️✔️: else switch (bond->getBondType()) {
    // RDKit✔️✔️: case SINGLE: res = "-"; break; case DOUBLE: res = "="; break;
    // RDKit✔️✔️: case TRIPLE: res = "#"; break; case AROMATIC: res = ":"; break;
    // RDKit✔️✔️: default: res = "?"; break; }
    // END RDKIT CPP FUNCTION getBondSymbol
    if bond.is_aromatic() {
        ":"
    } else {
        match bond.order() {
            BondOrder::Single => "-",
            BondOrder::Double => "=",
            BondOrder::Triple => "#",
            BondOrder::Aromatic => ":",
            _ => "?",
        }
    }
}

fn connectivity_ranks(
    topology: &TopologyBlock,
    atom_symbols: &[String],
    bond_symbols: &[String],
) -> Result<Vec<u32>, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION rankFragmentAtoms restricted call from runCleanup
    // RDKit✔️✔️: Canon::rankFragmentAtoms(mol, aranks, atomsInPlay, bondsInPlay,
    // RDKit✔️✔️:   &atomSymbols, &bondSymbols, false, false, false, false, false, false);
    // END RDKIT CPP FUNCTION rankFragmentAtoms restricted call from runCleanup
    // With every atom and bond in play and all six optional dimensions false,
    // the source comparison is the stable equitable refinement of the supplied
    // atom label and sorted incident (bond label, neighbor class) multiset.
    let count = topology.atoms.len();
    if count == 0 {
        return Ok(Vec::new());
    }
    let mut ranks = vec![0; count];
    for iteration in 0..=count {
        let keys = (0..count)
            .map(|index| {
                let mut neighbors = topology
                    .adjacency
                    .neighbors_of(index)
                    .iter()
                    .map(|neighbor| {
                        (
                            bond_symbols[neighbor.bond.index()].clone(),
                            ranks[neighbor.atom_index],
                        )
                    })
                    .collect::<Vec<_>>();
                neighbors.sort_unstable();
                (ranks[index], atom_symbols[index].clone(), neighbors)
            })
            .collect::<Vec<_>>();
        let unique = keys.iter().cloned().collect::<BTreeSet<_>>();
        let classes = unique
            .into_iter()
            .enumerate()
            .map(|(rank, key)| (key, rank as u32))
            .collect::<BTreeMap<_, _>>();
        let next = keys.iter().map(|key| classes[key]).collect::<Vec<_>>();
        if iteration > 0 && next == ranks {
            return Ok(next);
        }
        ranks = next;
    }
    Err(PotentialStereoError::RefinementDidNotConverge {
        iterations: count + 1,
    })
}

fn initialize_atoms(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
    known: &mut [bool],
    possible: &mut [bool],
    symbols: &mut [String],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION initAtomInfo
    // RDKit✔️❌: atomSymbols[aidx] = getAtomCompareSymbol(*atom);
    // RDKit✔️❌: if (detail::isAtomPotentialStereoAtom(atom, allowNontetrahedralStereo)) {
    // RDKit✔️❌:   auto sinfo = detail::getStereoInfo(atom);
    // RDKit✔️❌:   switch (sinfo.specified) {
    // RDKit✔️❌:   case Unknown: knownAtoms.set(aidx); atomSymbols[aidx] += std::to_string(aidx); break;
    // RDKit✔️❌:   case Chirality::StereoSpecified::Specified:
    // RDKit✔️❌:     knownAtoms.set(aidx);
    // RDKit✔️❌:     if (sinfo.descriptor == StereoDescriptor::Tet_CCW) {
    // RDKit✔️❌:       atomSymbols[aidx] += "_CCW";
    // RDKit✔️❌:     } else if (sinfo.descriptor == StereoDescriptor::Tet_CW) {
    // RDKit✔️❌:       atomSymbols[aidx] += "_CW";
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       atomSymbols[aidx] += "_STEREO";
    // RDKit✔️❌:     }
    // RDKit✔️❌:     break;
    // RDKit✔️❌:   case Unspecified: if (flagPossible) possibleAtoms.set(aidx); break;
    // RDKit✔️❌:   }
    // RDKit✔️❌: } else if (cleanIt) atom->setChiralTag(CHI_UNSPECIFIED);
    // END RDKIT CPP FUNCTION initAtomInfo
    for index in 0..topology.atoms.len() {
        let atom = AtomId::new(index);
        symbols[index] = atom_symbol(&topology.atoms[index]);
        if is_potential_atom(topology, valence, rings, atom, params.allow_nontetrahedral)? {
            let info = atom_info(topology, valence, atom, params.allow_nontetrahedral)?;
            match info.specified {
                PotentialStereoSpecified::Unknown => {
                    known[index] = true;
                    symbols[index].push_str(&index.to_string());
                }
                PotentialStereoSpecified::Specified => {
                    known[index] = true;
                    symbols[index].push_str(match info.descriptor {
                        PotentialStereoDescriptor::TetrahedralCounterclockwise => "_CCW",
                        PotentialStereoDescriptor::TetrahedralClockwise => "_CW",
                        _ => "_STEREO",
                    });
                }
                PotentialStereoSpecified::Unspecified if params.flag_possible => {
                    possible[index] = true;
                    if !params.clean {
                        symbols[index].push('_');
                        symbols[index].push_str(&index.to_string());
                    }
                }
                PotentialStereoSpecified::Unspecified => {}
            }
        } else if params.clean {
            topology.atoms[index].set_chiral_tag(ChiralTag::Unspecified);
        }
    }
    Ok(())
}

fn clear_bond_stereo_value(bond: &mut Bond) -> Result<(), PotentialStereoError> {
    bond.set_stereo(BondStereo::None)?;
    bond.set_stereo_atoms(None);
    Ok(())
}

fn initialize_bonds(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
    known: &mut [bool],
    possible: &mut [bool],
    symbols: &mut [String],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION initBondInfo
    // RDKit✔️❌: void initBondInfo(ROMol &mol, bool flagPossible, bool cleanIt,
    // RDKit✔️❌:                   boost::dynamic_bitset<> &knownBonds,
    // RDKit✔️❌:                   std::vector<std::string> &bondSymbols,
    // RDKit✔️❌:                   boost::dynamic_bitset<> &possibleBonds) {
    // RDKit✔️❌:   for (const auto bond : mol.bonds()) {
    // RDKit✔️❌:     auto bidx = bond->getIdx();
    // RDKit✔️❌:     bondSymbols[bidx] = getBondSymbol(bond);
    // RDKit✔️❌:     if (detail::isBondPotentialStereoBond(bond)) {
    // RDKit✔️❌:       auto sinfo = detail::getStereoInfo(bond);
    // RDKit✔️❌:       switch (sinfo.specified) {
    // RDKit✔️❌:         case Chirality::StereoSpecified::Unknown:
    // RDKit✔️❌:           knownBonds.set(bidx);
    // RDKit✔️❌:           bondSymbols[bidx] += "_" + std::to_string(bidx);
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Chirality::StereoSpecified::Specified:
    // RDKit✔️❌:           knownBonds.set(bidx);
    // RDKit✔️❌:           if (sinfo.descriptor == StereoDescriptor::Bond_Cis) {
    // RDKit✔️❌:             bondSymbols[bidx] += "_cis";
    // RDKit✔️❌:           } else if (sinfo.descriptor == StereoDescriptor::Bond_Trans) {
    // RDKit✔️❌:             bondSymbols[bidx] += "_trans";
    // RDKit✔️❌:           } else {
    // RDKit✔️❌:             bondSymbols[bidx] += "_STEREO";
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         case Chirality::StereoSpecified::Unspecified:
    // RDKit✔️❌:           if (flagPossible) {
    // RDKit✔️❌:             possibleBonds.set(bidx);
    // RDKit✔️❌:             if (!cleanIt) {
    // RDKit✔️❌:               bondSymbols[bidx] += "_" + std::to_string(bidx);
    // RDKit✔️❌:             }
    // RDKit✔️❌:           }
    // RDKit✔️❌:           break;
    // RDKit✔️❌:         default:
    // RDKit✔️❌:           throw ValueErrorException("bad StereoInfo.specified type");
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       auto currentStereo = bond->getStereo();
    // RDKit✔️❌:       if (currentStereo != Bond::BondStereo::STEREOATROPCW &&
    // RDKit✔️❌:           currentStereo != Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️❌:         if (cleanIt) {
    // RDKit✔️❌:           bond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         knownBonds.set(bidx);
    // RDKit✔️❌:         if (currentStereo == Bond::BondStereo::STEREOATROPCW) {
    // RDKit✔️❌:           bondSymbols[bidx] += "_atropcw";
    // RDKit✔️❌:         } else if (currentStereo == Bond::BondStereo::STEREOATROPCCW) {
    // RDKit✔️❌:           bondSymbols[bidx] += "_atropccw";
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION initBondInfo
    // The source loop changes stereo values only: atomic numbers, bond order,
    // endpoints, adjacency and RingInfo are fixed for this invocation. Share
    // its one lazily checked structural fact; do not rescan A+B+membership
    // for each eligible bond. Standalone calls never reuse this local scratch.
    // Constant space, one linear proof at most; fully sized/early-return paths
    // retain no proof. All source degree/H/init/error ordering stays below.
    let mut preserved_prefix = None;
    for index in 0..topology.bonds.len() {
        symbols[index] = bond_symbol(&topology.bonds[index]).to_owned();
        let candidate = {
            let bond = &topology.bonds[index];
            is_potential_bond_with_prefix_check(
                topology,
                valence,
                rings,
                bond,
                &mut preserved_prefix,
            )?
        };
        if candidate {
            let info = bond_info(topology, BondId::new(index))?;
            match info.specified {
                PotentialStereoSpecified::Unknown => {
                    known[index] = true;
                    symbols[index].push('_');
                    symbols[index].push_str(&index.to_string());
                }
                PotentialStereoSpecified::Specified => {
                    known[index] = true;
                    symbols[index].push_str(match info.descriptor {
                        PotentialStereoDescriptor::BondCis => "_cis",
                        PotentialStereoDescriptor::BondTrans => "_trans",
                        _ => "_STEREO",
                    });
                }
                PotentialStereoSpecified::Unspecified if params.flag_possible => {
                    possible[index] = true;
                    if !params.clean {
                        symbols[index].push('_');
                        symbols[index].push_str(&index.to_string());
                    }
                }
                PotentialStereoSpecified::Unspecified => {}
            }
        } else {
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       auto currentStereo = bond->getStereo();
            // RDKit✔️✔️:       if (currentStereo != Bond::BondStereo::STEREOATROPCW &&
            // RDKit✔️✔️:           currentStereo != Bond::BondStereo::STEREOATROPCCW) {
            // RDKit✔️✔️:         if (cleanIt) {
            // RDKit✔️✔️:           bond->setStereo(Bond::BondStereo::STEREONONE);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       } else {
            // RDKit✔️✔️:         knownBonds.set(bidx);
            // RDKit✔️✔️:         if (currentStereo == Bond::BondStereo::STEREOATROPCW) {
            // RDKit✔️✔️:           bondSymbols[bidx] += "_atropcw";
            // RDKit✔️✔️:         } else if (currentStereo == Bond::BondStereo::STEREOATROPCCW) {
            // RDKit✔️✔️:           bondSymbols[bidx] += "_atropccw";
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
            match topology.bonds[index].stereo() {
                BondStereo::AtropCw => {
                    known[index] = true;
                    symbols[index].push_str("_atropcw");
                }
                BondStereo::AtropCcw => {
                    known[index] = true;
                    symbols[index].push_str("_atropccw");
                }
                _ if params.clean => clear_bond_stereo_value(&mut topology.bonds[index])?,
                _ => {}
            }
        }
    }
    Ok(())
}

fn bond_between(topology: &TopologyBlock, left: AtomId, right: AtomId) -> Option<BondId> {
    topology
        .adjacency
        .neighbors_of(left.index())
        .iter()
        .find_map(|neighbor| (neighbor.atom_index == right.index()).then_some(neighbor.bond))
}

fn flag_ring_stereo(
    topology: &TopologyBlock,
    rings: &RingInfo,
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
    known_atoms: &[bool],
    possible_atoms: Option<&[bool]>,
    known_bonds: &[bool],
    possible_bonds: Option<&[bool]>,
) {
    // BEGIN RDKIT CPP FUNCTION flagRingStereo
    // RDKit✔️✔️: for (unsigned int ridx = 0; ridx < ringInfo->atomRings().size(); ++ridx) {
    // RDKit✔️✔️:   for (unsigned int ai = 0; ai < sz; ++ai) {
    // RDKit✔️✔️:     for (unsigned int ringDivisor : {2, 3}) {
    // RDKit✔️✔️:       bool ringIsMultipleOfDivisor = ((sz % ringDivisor) == 0);
    // RDKit✔️✔️:       auto incrementSize = sz / ringDivisor;
    // RDKit✔️✔️:       if (ringIsMultipleOfDivisor) {
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (ringInfo->numAtomRings(aidx) > 1) {
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (nHere > 1) {
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION flagRingStereo
    let active_atom =
        |index: usize| known_atoms[index] || possible_atoms.is_some_and(|rows| rows[index]);
    let active_bond =
        |index: usize| known_bonds[index] || possible_bonds.is_some_and(|rows| rows[index]);
    for (ring_index, ring) in rings.atom_rings().iter().enumerate() {
        let ring_bonds = &rings.bond_rings()[ring_index];
        let size = ring.len();
        if size == 0 {
            continue;
        }
        let half_size = size / 2 + usize::from(size % 2 != 0);
        let mut marked = vec![false; topology.atoms.len()];
        let mut count_here = 0usize;
        for position in 0..size {
            let atom = ring[position];
            if !active_atom(atom.index()) {
                continue;
            }
            for divisor in [2usize, 3] {
                if size % divisor != 0 {
                    continue;
                }
                let increment = size / divisor;
                let mut by_bond = 0usize;
                let mut by_atom = 0usize;
                let mut offset = increment;
                while offset < size {
                    let other = ring[(position + offset) % size];
                    let has_external =
                        topology
                            .adjacency
                            .neighbors_of(other.index())
                            .iter()
                            .any(|neighbor| {
                                active_bond(neighbor.bond.index())
                                    && !ring_bonds.contains(&neighbor.bond)
                            });
                    if has_external {
                        by_bond += 1;
                    } else if by_bond == 0 && active_atom(other.index()) {
                        by_atom += 1;
                    }
                    offset += increment;
                }
                if by_bond == divisor - 1 || by_atom == divisor - 1 {
                    count_here += 1 + by_bond;
                    let mut offset = 0;
                    while offset < size {
                        marked[ring[(position + offset) % size].index()] = true;
                        offset += increment;
                    }
                }
            }
            if rings.num_atom_rings(atom) > 1 {
                let mut previous = atom;
                for step in 1..=half_size {
                    let other = ring[(position + step) % size];
                    let Some(edge) = bond_between(topology, previous, other) else {
                        break;
                    };
                    if rings.num_bond_rings(edge) < 2 {
                        break;
                    }
                    if active_atom(other.index()) {
                        count_here += 2;
                        marked[atom.index()] = true;
                        marked[other.index()] = true;
                        break;
                    }
                    previous = other;
                }
            }
        }
        if count_here > 1 {
            for atom in ring {
                if marked[atom.index()] {
                    possible_ring_atoms[atom.index()] += 1;
                }
            }
            for bond in ring_bonds {
                possible_ring_bonds[bond.index()] += 1;
            }
        }
    }
}

fn controlling_atoms_are_duplicates(
    topology: &TopologyBlock,
    rings: &RingInfo,
    bond: BondId,
    left: AtomId,
    right: AtomId,
    ranks: &[u32],
    possible_atoms: &[bool],
    known_atoms: &[bool],
    possible_bonds: &[bool],
    known_bonds: &[bool],
) -> bool {
    // BEGIN RDKIT CPP FUNCTION areStereobondControllingAtomsDupes
    // RDKit✔️✔️: if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) return false;
    // RDKit✔️✔️: if (atomRanks[controllingAtom1] != atomRanks[controllingAtom2]) {
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION areStereobondControllingAtomsDupes
    if ranks[left.index()] != ranks[right.index()] {
        return false;
    }
    let left_members = rings.atom_members(left);
    let right_members = rings.atom_members(right);
    let mut li = 0;
    let mut ri = 0;
    while li < left_members.len() && ri < right_members.len() {
        match left_members[li].cmp(&right_members[ri]) {
            std::cmp::Ordering::Less => li += 1,
            std::cmp::Ordering::Greater => ri += 1,
            std::cmp::Ordering::Equal => {
                let ring = &rings.atom_rings()[left_members[li]];
                li += 1;
                ri += 1;
                if ring.len() % 2 != 0 {
                    continue;
                }
                for endpoint in [
                    topology.bonds[bond.index()].begin(),
                    topology.bonds[bond.index()].end(),
                ] {
                    let Some(position) = ring.iter().position(|atom| *atom == endpoint) else {
                        continue;
                    };
                    let opposite = ring[(position + ring.len() / 2) % ring.len()];
                    if possible_atoms[opposite.index()] || known_atoms[opposite.index()] {
                        return false;
                    }
                    if graph_degree(topology, opposite) == 3 {
                        for neighbor in topology.adjacency.neighbors_of(opposite.index()) {
                            if !ring.contains(&AtomId::new(neighbor.atom_index))
                                && (possible_bonds[neighbor.bond.index()]
                                    || known_bonds[neighbor.bond.index()])
                            {
                                return false;
                            }
                        }
                    }
                }
            }
        }
    }
    true
}

#[allow(clippy::too_many_arguments)]
fn update_atoms(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    ranks: &[u32],
    symbols: &mut [String],
    possible: &mut [bool],
    known: &[bool],
    fixed: &mut [bool],
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
    rings: &RingInfo,
    allow_nontetrahedral: bool,
    output: &mut Vec<PotentialStereoInfo>,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION updateAtoms
    // RDKit✔️✔️: if (knownAtoms[aidx] || possibleAtoms[aidx]) {
    // RDKit✔️✔️:   auto sinfo = detail::getStereoInfo(atom);
    // RDKit✔️✔️:   if (fixedAtoms[aidx]) sinfos.push_back(std::move(sinfo));
    // RDKit✔️✔️:   if (fixedAtoms[aidx]) {
    // RDKit✔️✔️:     sinfos.push_back(std::move(sinfo));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::vector<unsigned int> nbrs;
    // RDKit✔️✔️:     nbrs.reserve(sinfo.controllingAtoms.size());
    // RDKit✔️✔️:     bool haveADupe = false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION updateAtoms
    let mut another = false;
    for index in 0..topology.atoms.len() {
        if !known[index] && !possible[index] {
            continue;
        }
        let atom = AtomId::new(index);
        let mut info = atom_info(topology, valence, atom, allow_nontetrahedral)?;
        if fixed[index] {
            output.push(info);
            continue;
        }
        let mut neighbor_ranks = Vec::new();
        let mut duplicate = false;
        if info.stereo_type == PotentialStereoType::AtomTetrahedral {
            for neighbor in info.controlling_atoms.iter().flatten().copied() {
                let rank = ranks[neighbor.index()];
                if neighbor_ranks.contains(&rank) {
                    if possible_ring_atoms[index] > 0 {
                        let transmitting = bond_between(topology, atom, neighbor)
                            .is_some_and(|bond| possible_ring_bonds[bond.index()] > 0);
                        if !transmitting {
                            duplicate = true;
                            break;
                        }
                    } else {
                        duplicate = true;
                        break;
                    }
                } else {
                    neighbor_ranks.push(rank);
                }
            }
        }
        if !duplicate {
            let mut next_symbol = symbols[index].clone();
            if !possible[index] {
                let mut sorted = neighbor_ranks.clone();
                sorted.sort_unstable();
                if info.stereo_type == PotentialStereoType::AtomTetrahedral
                    && count_swaps_to_interconvert(&neighbor_ranks, &sorted)? % 2 == 1
                {
                    info.descriptor = match info.descriptor {
                        PotentialStereoDescriptor::TetrahedralClockwise => {
                            PotentialStereoDescriptor::TetrahedralCounterclockwise
                        }
                        PotentialStereoDescriptor::TetrahedralCounterclockwise => {
                            PotentialStereoDescriptor::TetrahedralClockwise
                        }
                        descriptor => descriptor,
                    };
                }
                next_symbol = atom_symbol(&topology.atoms[index]);
                next_symbol.push_str(match info.descriptor {
                    PotentialStereoDescriptor::TetrahedralClockwise => "_CW",
                    PotentialStereoDescriptor::TetrahedralCounterclockwise => "_CCW",
                    _ => "",
                });
                fixed[index] = true;
            }
            if symbols[index] != next_symbol {
                symbols[index] = next_symbol;
                another = true;
            }
            output.push(info);
        } else {
            another |= possible[index];
            possible[index] = false;
            symbols[index] = atom_symbol(&topology.atoms[index]);
            if possible_ring_atoms[index] > 0 {
                possible_ring_atoms[index] = 0;
                another = true;
                for (ring_index, ring) in rings.atom_rings().iter().enumerate() {
                    let mut remaining = 0usize;
                    for ring_atom in ring {
                        fixed[ring_atom.index()] = false;
                        remaining += usize::from(possible_ring_atoms[ring_atom.index()] > 0);
                    }
                    if remaining <= 1 {
                        if remaining == 1 {
                            if let Some(last) = ring
                                .iter()
                                .find(|ring_atom| possible_ring_atoms[ring_atom.index()] > 0)
                            {
                                possible_ring_atoms[last.index()] -= 1;
                            }
                        }
                        for ring_bond in &rings.bond_rings()[ring_index] {
                            if possible_ring_bonds[ring_bond.index()] > 0 {
                                possible_ring_bonds[ring_bond.index()] -= 1;
                            }
                        }
                    }
                }
            }
        }
    }
    Ok(another)
}

fn swap_bond_descriptor(info: &mut PotentialStereoInfo) {
    info.descriptor = match info.descriptor {
        PotentialStereoDescriptor::BondCis => PotentialStereoDescriptor::BondTrans,
        PotentialStereoDescriptor::BondTrans => PotentialStereoDescriptor::BondCis,
        descriptor => descriptor,
    };
}

#[allow(clippy::too_many_arguments)]
fn update_bonds(
    topology: &TopologyBlock,
    ranks: &[u32],
    symbols: &mut [String],
    possible_atoms: &[bool],
    possible_bonds: &mut [bool],
    known_atoms: &[bool],
    known_bonds: &[bool],
    fixed_bonds: &mut [bool],
    rings: &RingInfo,
    output: &mut Vec<PotentialStereoInfo>,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION updateBonds
    // RDKit✔️✔️: if (knownBonds[bidx] || possibleBonds[bidx]) {
    // RDKit✔️✔️:   auto sinfo = detail::getStereoInfo(bond);
    // RDKit✔️✔️:   if (both controlling slots missing on either end) fixedBonds.set(bidx);
    // RDKit✔️✔️:   if (!fixedBonds[bidx]) {
    // RDKit✔️✔️:     if (sinfo.controllingAtoms[0] != Atom::NOATOM &&
    // RDKit✔️✔️:         sinfo.controllingAtoms[1] != Atom::NOATOM) {
    // RDKit✔️✔️:     if (haveADupe && possibleBonds[bidx]) possibleBonds[bidx] = 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION updateBonds
    let mut another = false;
    for index in 0..topology.bonds.len() {
        if !known_bonds[index] && !possible_bonds[index] {
            continue;
        }
        let bond_id = BondId::new(index);
        // RDKit✔️✔️:       if (sinfo.type == Chirality::StereoType::Unspecified) {
        // RDKit✔️✔️:         continue;  // not a double bond nor an atropisomer bond
        // RDKit✔️✔️:       }
        let bond = &topology.bonds[index];
        if bond.order() != BondOrder::Double
            && !(bond.order() == BondOrder::Single
                && matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
        {
            continue;
        }
        let mut info = bond_info(topology, bond_id)?;
        if info.controlling_atoms.len() != 4 {
            return Err(PotentialStereoError::InvalidStereoReferences {
                bond: bond_id,
                reason: "potential stereo bond must have four controlling slots",
            });
        }
        if (info.controlling_atoms[0].is_none() && info.controlling_atoms[1].is_none())
            || (info.controlling_atoms[2].is_none() && info.controlling_atoms[3].is_none())
        {
            if info.specified == PotentialStereoSpecified::Specified {
                return Err(PotentialStereoError::InvalidStereoReferences {
                    bond: bond_id,
                    reason: "specified bond has no controlling atom on one endpoint",
                });
            }
            fixed_bonds[index] = true;
        }
        if fixed_bonds[index] {
            output.push(info);
            continue;
        }
        let mut duplicate = false;
        let mut needs_swap = false;
        for offset in [0usize, 2] {
            if let (Some(left), Some(right)) = (
                info.controlling_atoms[offset],
                info.controlling_atoms[offset + 1],
            ) {
                if controlling_atoms_are_duplicates(
                    topology,
                    rings,
                    bond_id,
                    left,
                    right,
                    ranks,
                    possible_atoms,
                    known_atoms,
                    possible_bonds,
                    known_bonds,
                ) {
                    duplicate = true;
                } else if ranks[left.index()] < ranks[right.index()] {
                    info.controlling_atoms.swap(offset, offset + 1);
                    needs_swap = !needs_swap;
                }
            }
        }
        if !duplicate {
            if needs_swap {
                swap_bond_descriptor(&mut info);
            }
            let mut next_symbol = symbols[index].clone();
            match (info.specified, info.descriptor) {
                (PotentialStereoSpecified::Specified, PotentialStereoDescriptor::BondCis) => {
                    next_symbol.push_str("_cis")
                }
                (PotentialStereoSpecified::Specified, PotentialStereoDescriptor::BondTrans) => {
                    next_symbol.push_str("_trans")
                }
                (PotentialStereoSpecified::Unknown, _) => next_symbol.push_str("_unk"),
                _ => {}
            }
            if symbols[index] != next_symbol {
                symbols[index] = next_symbol;
                another = true;
            }
            if !possible_bonds[index] {
                fixed_bonds[index] = true;
            }
            output.push(info);
        } else if possible_bonds[index] {
            possible_bonds[index] = false;
            symbols[index] = bond_symbol(&topology.bonds[index]).to_owned();
            another = true;
        }
    }
    Ok(another)
}

fn clean_invalid_stereo(
    topology: &mut TopologyBlock,
    fixed_atoms: &[bool],
    known_atoms: &[bool],
    fixed_bonds: &[bool],
    known_bonds: &[bool],
) -> Result<(), PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION cleanMolStereo
    // RDKit✔️✔️: if (!fixedAtoms[i] && knownAtoms[i]) {
    // RDKit✔️✔️:   if (atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:       atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (!fixedBonds[i] && knownBonds[i]) {
    // RDKit✔️✔️:   bond->setStereo(STEREONONE); bond->setBondDir(NONE);
    // RDKit✔️✔️:   bond->getStereoAtoms().clear();
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (removedStereo) {
    // END RDKIT CPP FUNCTION cleanMolStereo
    let mut wedge_centers = Vec::new();
    for index in 0..topology.atoms.len() {
        if fixed_atoms[index] || !known_atoms[index] {
            continue;
        }
        match topology.atoms[index].chiral_tag() {
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw => {
                topology.atoms[index].set_chiral_tag(ChiralTag::Unspecified);
                wedge_centers.push(AtomId::new(index));
            }
            ChiralTag::Tetrahedral
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral => {
                topology.atoms[index].set_chiral_permutation(Some(0));
            }
            _ => {}
        }
    }
    for atom in wedge_centers {
        let bonds = topology
            .adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| neighbor.bond)
            .collect::<Vec<_>>();
        for bond in bonds {
            if matches!(
                topology.bonds[bond.index()].direction(),
                BondDirection::BeginDash | BondDirection::BeginWedge
            ) {
                topology.bonds[bond.index()].set_direction(BondDirection::None);
            }
        }
    }
    let mut removed_bond_stereo = false;
    for index in 0..topology.bonds.len() {
        if !fixed_bonds[index] && known_bonds[index] {
            clear_bond_stereo_value(&mut topology.bonds[index])?;
            topology.bonds[index].set_direction(BondDirection::None);
            removed_bond_stereo = true;
        }
    }
    if removed_bond_stereo {
        for index in 0..topology.bonds.len() {
            if !matches!(
                topology.bonds[index].direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ) {
                continue;
            }
            let endpoints = [topology.bonds[index].begin(), topology.bonds[index].end()];
            let mut direction_ok = false;
            for endpoint in endpoints {
                direction_ok = topology
                    .adjacency
                    .neighbors_of(endpoint.index())
                    .iter()
                    .any(|neighbor| {
                        neighbor.bond.index() != index
                            && topology.bonds[neighbor.bond.index()].stereo() != BondStereo::None
                    });
                if !direction_ok {
                    topology.bonds[index].set_direction(BondDirection::None);
                }
            }
        }
    }
    Ok(())
}

fn atom_is_candidate_for_ring_stereochemistry(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
    atom: AtomId,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION atomIsCandidateForRingStereochem
    // RDKit✔️❌: if (ringInfo->isInitialized() && ringInfo->numAtomRings(atom->getIdx())) {
    // RDKit✔️❌:   if (atom->getAtomicNum() == 7 && atom->getTotalDegree() == 3 &&
    // RDKit✔️❌:       !ringInfo->isAtomInRingOfSize(atom->getIdx(), 3) &&
    // RDKit✔️❌:       !queryIsAtomBridgehead(atom)) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<const Atom *> nonRingNbrs;
    // RDKit✔️❌:   std::vector<const Atom *> ringNbrs;
    // RDKit✔️❌:   std::set<unsigned int> ringNbrRanks;
    // RDKit✔️❌:   for (const auto bond : mol.atomBonds(atom)) {
    // RDKit✔️❌:     if (!ringInfo->numBondRings(bond->getIdx())) {
    // RDKit✔️❌:       nonRingNbrs.push_back(bond->getOtherAtom(atom));
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       const Atom *nbr = bond->getOtherAtom(atom);
    // RDKit✔️❌:       ringNbrs.push_back(nbr);
    // RDKit✔️❌:       ringNbrRanks.insert(atomRanks[nbr->getIdx()]);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION atomIsCandidateForRingStereochem
    // The detached helper repeats validated total-H lookup for each candidate,
    // so behavior is reproduced while the performance axis remains lower.
    if !rings.is_initialized() || rings.num_atom_rings(atom) == 0 {
        return Ok(false);
    }
    let value = &topology.atoms[atom.index()];
    if value.atomic_number() == 7
        && total_degree(topology, valence, atom)? == 3
        && !rings.is_atom_in_ring_of_size(atom, 3)
        && is_atom_bridgehead_from_topology(topology, atom.index(), rings) == 0
    {
        return Ok(false);
    }
    let mut non_ring_neighbors = Vec::new();
    let mut ring_neighbor_count = 0usize;
    let mut ring_neighbor_ranks = BTreeSet::new();
    for neighbor in topology.adjacency.neighbors_of(atom.index()) {
        let neighbor_atom = AtomId::new(neighbor.atom_index);
        if rings.num_bond_rings(neighbor.bond) == 0 {
            non_ring_neighbors.push(neighbor_atom);
        } else {
            ring_neighbor_count += 1;
            ring_neighbor_ranks.insert(ranks[neighbor.atom_index]);
        }
    }
    Ok(match non_ring_neighbors.as_slice() {
        [left, right] => {
            ranks[left.index()] != ranks[right.index()]
                && ring_neighbor_count != ring_neighbor_ranks.len()
        }
        [_] => ring_neighbor_count > ring_neighbor_ranks.len(),
        [] => {
            (ring_neighbor_count == 4 && ring_neighbor_ranks.len() == 3)
                || (ring_neighbor_count == 3 && ring_neighbor_ranks.len() == 2)
        }
        _ => false,
    })
}

pub(crate) struct RingPropertyUpdate {
    pub atom: AtomId,
    pub key: &'static str,
    pub value: PropertyValue,
    pub computed: bool,
}

pub(crate) struct RingSpecialCases {
    pub relations: Vec<RingStereoRelation>,
    pub flags: Vec<bool>,
    pub updates: Vec<RingPropertyUpdate>,
}

fn signed_ring_reference(atom: AtomId, same: bool) -> Result<i32, PotentialStereoError> {
    let value = atom
        .index()
        .checked_add(1)
        .and_then(|n| i32::try_from(n).ok())
        .ok_or(PotentialStereoError::InvalidRingStereoReference {
            atom,
            value: i32::MAX,
        })?;
    Ok(if same { value } else { -value })
}

fn ring_reference_index(
    owner: AtomId,
    value: i32,
    count: usize,
) -> Result<usize, PotentialStereoError> {
    let index = value
        .checked_abs()
        .and_then(|n| n.checked_sub(1))
        .and_then(|n| usize::try_from(n).ok())
        .filter(|n| *n < count)
        .ok_or(PotentialStereoError::InvalidRingStereoReference { atom: owner, value })?;
    Ok(index)
}

fn ring_vector(
    topology: &TopologyBlock,
    overlay: &[Option<Vec<i32>>],
    atom: AtomId,
) -> Result<Vec<i32>, PotentialStereoError> {
    // RDKit❗✔️: ratom->getPropIfPresent(common_properties::_ringStereoAtoms,
    // RDKit❗✔️:                         oringatoms);
    // Exact tag read when reached; copy only the requested detached vector.
    if let Some(value) = &overlay[atom.index()] {
        return Ok(value.clone());
    }
    match topology.atoms[atom.index()].prop("_ringStereoAtoms") {
        None => Ok(Vec::new()),
        Some(PropertyValue::IntVector(value)) => Ok(value.clone()),
        Some(value) => Err(PotentialStereoError::InvalidPropertyKind {
            atom,
            property: "_ringStereoAtoms",
            kind: value.kind(),
        }),
    }
}

pub(crate) fn special_ring_cases(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    input_rings: &RingInfo,
    ranks: &[u32],
) -> Result<RingSpecialCases, PotentialStereoError> {
    // RDKit❗❌:   boost::dynamic_bitset<> atomsSeen(mol.getNumAtoms());
    // RDKit❗❌:   boost::dynamic_bitset<> atomsUsed(mol.getNumAtoms());
    // RDKit❗❌:   boost::dynamic_bitset<> bondsSeen(mol.getNumBonds());
    // RDKit❗❌:
    // RDKit❗❌:   for (const auto atom : mol.atoms()) {
    // RDKit❗❌:     if (atomsSeen[atom->getIdx()]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (atom->getChiralTag() == Atom::CHI_UNSPECIFIED ||
    // RDKit❗❌:         atom->hasProp(common_properties::_CIPCode) ||
    // RDKit❗❌:         !mol.getRingInfo()->numAtomRings(atom->getIdx()) ||
    // RDKit❗❌:         !atomIsCandidateForRingStereochem(mol, atom, atomRanks)) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     // do a BFS from this ring atom along ring bonds and find other
    // RDKit❗❌:     // stereochemistry candidates.
    // RDKit❗❌:     std::list<const Atom *> nextAtoms;
    // RDKit❗❌:     // start with finding viable neighbors
    // RDKit❗❌:     for (const auto bond : mol.atomBonds(atom)) {
    // RDKit❗❌:       unsigned int bidx = bond->getIdx();
    // RDKit❗❌:       if (!bondsSeen[bidx]) {
    // RDKit❗❌:         bondsSeen.set(bidx);
    // RDKit❗❌:         if (mol.getRingInfo()->numBondRings(bidx)) {
    // RDKit❗❌:           const Atom *oatom = bond->getOtherAtom(atom);
    // RDKit❗❌:           if (!atomsSeen[oatom->getIdx()]) {
    // RDKit❗❌:             nextAtoms.push_back(oatom);
    // RDKit❗❌:             atomsUsed.set(oatom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     INT_VECT ringStereoAtoms(0);
    // RDKit❗❌:     if (!nextAtoms.empty()) {
    // RDKit❗❌:       atom->getPropIfPresent(common_properties::_ringStereoAtoms,
    // RDKit❗❌:                              ringStereoAtoms);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     while (!nextAtoms.empty()) {
    // RDKit❗❌:       const Atom *ratom = nextAtoms.front();
    // RDKit❗❌:       nextAtoms.pop_front();
    // RDKit❗❌:       atomsSeen.set(ratom->getIdx());
    // RDKit❗❌:       if (ratom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗❌:           !ratom->hasProp(common_properties::_CIPCode) &&
    // RDKit❗❌:           atomIsCandidateForRingStereochem(mol, ratom, atomRanks)) {
    // RDKit❗❌:         int same = (ratom->getChiralTag() == atom->getChiralTag()) ? 1 : -1;
    // RDKit❗❌:         ringStereoAtoms.push_back(same * (ratom->getIdx() + 1));
    // RDKit❗❌:         INT_VECT oringatoms(0);
    // RDKit❗❌:         ratom->getPropIfPresent(common_properties::_ringStereoAtoms,
    // RDKit❗❌:                                 oringatoms);
    // RDKit❗❌:         oringatoms.push_back(same * (atom->getIdx() + 1));
    // RDKit❗❌:         ratom->setProp(common_properties::_ringStereoAtoms, oringatoms, true);
    // RDKit❗❌:         possibleSpecialCases.set(ratom->getIdx());
    // RDKit❗❌:         possibleSpecialCases.set(atom->getIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:       // now push this atom's neighbors
    // RDKit❗❌:       for (const auto bond : mol.atomBonds(ratom)) {
    // RDKit❗❌:         unsigned int bidx = bond->getIdx();
    // RDKit❗❌:         if (!bondsSeen[bidx]) {
    // RDKit❗❌:           bondsSeen.set(bidx);
    // RDKit❗❌:           if (mol.getRingInfo()->numBondRings(bidx)) {
    // RDKit❗❌:             const Atom *oatom = bond->getOtherAtom(ratom);
    // RDKit❗❌:             if (!atomsSeen[oatom->getIdx()] && !atomsUsed[oatom->getIdx()]) {
    // RDKit❗❌:               nextAtoms.push_back(oatom);
    // RDKit❗❌:               atomsUsed.set(oatom->getIdx());
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }  // end of BFS
    // RDKit❗❌:     if (ringStereoAtoms.size() != 0) {
    // RDKit❗❌:       atom->setProp(common_properties::_ringStereoAtoms, ringStereoAtoms, true);
    // RDKit❗❌:       // because we're only going to hit each ring atom once, the first atom we
    // RDKit❗❌:       // encounter in a ring is going to end up with all the other atoms set as
    // RDKit❗❌:       // stereoAtoms, but each of them will only have the first atom present. We
    // RDKit❗❌:       // need to fix that. because the traverse from the first atom only
    // RDKit❗❌:       // followed ring bonds, these things are all by definition in one ring
    // RDKit❗❌:       // system. (Q: is this true if there's a spiro center in there?)
    // RDKit❗❌:       INT_VECT same(mol.getNumAtoms(), 0);
    // RDKit❗❌:       for (auto ringAtomEntry : ringStereoAtoms) {
    // RDKit❗❌:         int ringAtomIdx =
    // RDKit❗❌:             ringAtomEntry < 0 ? -ringAtomEntry - 1 : ringAtomEntry - 1;
    // RDKit❗❌:         same[ringAtomIdx] = ringAtomEntry;
    // RDKit❗❌:       }
    // RDKit❗❌:       for (INT_VECT_CI rae = ringStereoAtoms.begin();
    // RDKit❗❌:            rae != ringStereoAtoms.end(); ++rae) {
    // RDKit❗❌:         int ringAtomEntry = *rae;
    // RDKit❗❌:         int ringAtomIdx =
    // RDKit❗❌:             ringAtomEntry < 0 ? -ringAtomEntry - 1 : ringAtomEntry - 1;
    // RDKit❗❌:         INT_VECT lringatoms(0);
    // RDKit❗❌:         mol.getAtomWithIdx(ringAtomIdx)
    // RDKit❗❌:             ->getPropIfPresent(common_properties::_ringStereoAtoms, lringatoms);
    // RDKit❗❌:         CHECK_INVARIANT(lringatoms.size() > 0, "no other ring atoms found.");
    // RDKit❗❌:         for (auto orae = rae + 1; orae != ringStereoAtoms.end(); ++orae) {
    // RDKit❗❌:           int oringAtomEntry = *orae;
    // RDKit❗❌:           int oringAtomIdx =
    // RDKit❗❌:               oringAtomEntry < 0 ? -oringAtomEntry - 1 : oringAtomEntry - 1;
    // RDKit❗❌:           int theseDifferent = (ringAtomEntry < 0) ^ (oringAtomEntry < 0);
    // RDKit❗❌:           lringatoms.push_back(theseDifferent ? -(oringAtomIdx + 1)
    // RDKit❗❌:                                               : (oringAtomIdx + 1));
    // RDKit❗❌:           INT_VECT olringatoms(0);
    // RDKit❗❌:           mol.getAtomWithIdx(oringAtomIdx)
    // RDKit❗❌:               ->getPropIfPresent(common_properties::_ringStereoAtoms,
    // RDKit❗❌:                                  olringatoms);
    // RDKit❗❌:           CHECK_INVARIANT(olringatoms.size() > 0, "no other ring atoms found.");
    // RDKit❗❌:           olringatoms.push_back(theseDifferent ? -(ringAtomIdx + 1)
    // RDKit❗❌:                                                : (ringAtomIdx + 1));
    // RDKit❗❌:           mol.getAtomWithIdx(oringAtomIdx)
    // RDKit❗❌:               ->setProp(common_properties::_ringStereoAtoms, olringatoms);
    // RDKit❗❌:         }
    // RDKit❗❌:         mol.getAtomWithIdx(ringAtomIdx)
    // RDKit❗❌:             ->setProp(common_properties::_ringStereoAtoms, lringatoms);
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:     } else {
    // RDKit❗❌:       possibleSpecialCases.reset(atom->getIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:     atomsSeen.set(atom->getIdx());
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // RDKit❗❌:
    // RDKit❗❌: std::pair<bool, bool> isAtomPotentialChiralCenter(
    // RDKit❗❌:     const Atom *atom, const ROMol &mol, const UINT_VECT &ranks,
    // RDKit❗❌:     Chirality::INT_PAIR_VECT &nbrs) {
    // Behavior: preserve source BFS, visited-edge guards, lazy property reads,
    // duplicate signed references and computed membership. Detached changes are
    // applied only by the legacy adapter, never by relation-only perception.
    // Complexity: indexed bit vectors and ordered queue match source traversal;
    // copy only touched property vectors. Final relation presentation uses a
    // tree set (existing public order), so performance marker remains worse.
    let prepared_rings;
    let rings = if input_rings.is_symm_sssr() {
        input_rings
    } else {
        prepared_rings = crate::symmetrized_sssr(topology, &crate::RingSearchParams::default())?;
        &prepared_rings
    };
    let count = topology.atoms.len();
    let mut seen = vec![false; count];
    let mut used = vec![false; count];
    let mut bonds_seen = vec![false; topology.bonds.len()];
    let mut flags = vec![false; count];
    let mut vectors: Vec<Option<Vec<i32>>> = vec![None; count];
    let mut caches = vec![None; count];
    let mut computed = vec![false; count];
    let mut order = Vec::<(usize, bool)>::new();
    for start in 0..count {
        if seen[start] {
            continue;
        }
        let atom = AtomId::new(start);
        let value = &topology.atoms[start];
        if value.chiral_tag() == ChiralTag::Unspecified
            || value.prop("_CIPCode").is_some()
            || rings.num_atom_rings(atom) == 0
            || !ring_candidate_cached(
                topology,
                valence,
                rings,
                ranks,
                atom,
                &mut caches,
                &mut order,
            )?
        {
            continue;
        }
        let mut queue = VecDeque::new();
        for neighbor in topology.adjacency.neighbors_of(start) {
            let b = neighbor.bond.index();
            if !bonds_seen[b] {
                bonds_seen[b] = true;
                if rings.num_bond_rings(neighbor.bond) != 0 && !seen[neighbor.atom_index] {
                    queue.push_back(neighbor.atom_index);
                    used[neighbor.atom_index] = true;
                }
            }
        }
        let mut start_vector = if queue.is_empty() {
            Vec::new()
        } else {
            ring_vector(topology, &vectors, atom)?
        };
        while let Some(other) = queue.pop_front() {
            seen[other] = true;
            let other_atom = AtomId::new(other);
            let other_value = &topology.atoms[other];
            if other_value.chiral_tag() != ChiralTag::Unspecified
                && other_value.prop("_CIPCode").is_none()
                && ring_candidate_cached(
                    topology,
                    valence,
                    rings,
                    ranks,
                    other_atom,
                    &mut caches,
                    &mut order,
                )?
            {
                let same = value.chiral_tag() == other_value.chiral_tag();
                start_vector.push(signed_ring_reference(other_atom, same)?);
                let mut other_vector = ring_vector(topology, &vectors, other_atom)?;
                other_vector.push(signed_ring_reference(atom, same)?);
                if vectors[other].is_none() {
                    order.push((other, false));
                }
                vectors[other] = Some(other_vector);
                computed[other] = true;
                flags[other] = true;
                flags[start] = true;
            }
            for neighbor in topology.adjacency.neighbors_of(other) {
                let b = neighbor.bond.index();
                if !bonds_seen[b] {
                    bonds_seen[b] = true;
                    if rings.num_bond_rings(neighbor.bond) != 0
                        && !seen[neighbor.atom_index]
                        && !used[neighbor.atom_index]
                    {
                        queue.push_back(neighbor.atom_index);
                        used[neighbor.atom_index] = true;
                    }
                }
            }
        }
        if !start_vector.is_empty() {
            if vectors[start].is_none() {
                order.push((start, false));
            }
            vectors[start] = Some(start_vector.clone());
            computed[start] = true;
            // Source first constructs the complete signed index table, before
            // reading any reciprocal array. Preserve invalid-index error order.
            let mut same = vec![0; count];
            for &entry in &start_vector {
                same[ring_reference_index(atom, entry, count)?] = entry;
            }
            for (position, &entry) in start_vector.iter().enumerate() {
                let index = ring_reference_index(atom, entry, count)?;
                let owner = AtomId::new(index);
                let mut local = ring_vector(topology, &vectors, owner)?;
                if local.is_empty() {
                    return Err(PotentialStereoError::EmptyRingStereoReferences { atom: owner });
                }
                for &other_entry in &start_vector[position + 1..] {
                    let other = ring_reference_index(atom, other_entry, count)?;
                    let same = (entry < 0) == (other_entry < 0);
                    local.push(signed_ring_reference(AtomId::new(other), same)?);
                    let mut other_vector = ring_vector(topology, &vectors, AtomId::new(other))?;
                    if other_vector.is_empty() {
                        return Err(PotentialStereoError::EmptyRingStereoReferences {
                            atom: AtomId::new(other),
                        });
                    }
                    other_vector.push(signed_ring_reference(owner, same)?);
                    if vectors[other].is_none() {
                        order.push((other, false));
                    }
                    vectors[other] = Some(other_vector);
                }
                if vectors[index].is_none() {
                    order.push((index, false));
                }
                vectors[index] = Some(local);
            }
        } else {
            flags[start] = false;
        }
        seen[start] = true;
    }
    let mut relations = BTreeSet::new();
    for (index, row) in vectors.iter().enumerate() {
        if let Some(row) = row {
            for &entry in row {
                relations.insert(RingStereoRelation {
                    atom: AtomId::new(index),
                    other: AtomId::new(ring_reference_index(AtomId::new(index), entry, count)?),
                    same_orientation: entry > 0,
                });
            }
        }
    }
    let mut updates = Vec::with_capacity(order.len());
    for (index, cache) in order {
        updates.push(if cache {
            RingPropertyUpdate {
                atom: AtomId::new(index),
                key: "_ringStereochemCand",
                value: PropertyValue::Bool(caches[index].take().expect("registered cache update")),
                computed: true,
            }
        } else {
            RingPropertyUpdate {
                atom: AtomId::new(index),
                key: "_ringStereoAtoms",
                value: PropertyValue::IntVector(
                    vectors[index].take().expect("registered vector update"),
                ),
                computed: computed[index],
            }
        });
    }
    Ok(RingSpecialCases {
        relations: relations.into_iter().collect(),
        flags,
        updates,
    })
}

#[allow(clippy::too_many_arguments)]
fn ring_candidate_cached(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
    atom: AtomId,
    caches: &mut [Option<bool>],
    order: &mut Vec<(usize, bool)>,
) -> Result<bool, PotentialStereoError> {
    // RDKit❗✔️: bool res = false;
    // RDKit❗✔️: if (!atom->getPropIfPresent(common_properties::_ringStereochemCand, res)) {
    // RDKit❗✔️: atom->setProp(common_properties::_ringStereochemCand, res, 1);
    // Cached values are read only when the source outer/BFS guard reaches this
    // helper. Missing N exclusion returns without a computed cache write.
    if let Some(value) = caches[atom.index()] {
        return Ok(value);
    }
    if let Some(value) = topology.atoms[atom.index()].prop("_ringStereochemCand") {
        return value
            .as_bool()
            .map_err(|_| PotentialStereoError::InvalidPropertyKind {
                atom,
                property: "_ringStereochemCand",
                kind: value.kind(),
            });
    }
    let source_atom = &topology.atoms[atom.index()];
    if rings.is_initialized()
        && rings.num_atom_rings(atom) != 0
        && source_atom.atomic_number() == 7
        && total_degree(topology, valence, atom)? == 3
        && !rings.is_atom_in_ring_of_size(atom, 3)
        && is_atom_bridgehead_from_topology(topology, atom.index(), rings) == 0
    {
        return Ok(false);
    }
    let value = atom_is_candidate_for_ring_stereochemistry(topology, valence, rings, ranks, atom)?;
    caches[atom.index()] = Some(value);
    order.push((atom.index(), true));
    Ok(value)
}

pub(crate) fn special_ring_relations(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
) -> Result<Vec<RingStereoRelation>, PotentialStereoError> {
    Ok(special_ring_cases(topology, valence, rings, ranks)?.relations)
}

#[allow(clippy::too_many_arguments)]
fn refinement_loop(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    allow_nontetrahedral: bool,
    atom_symbols: &mut [String],
    bond_symbols: &mut [String],
    possible_atoms: &mut [bool],
    known_atoms: &[bool],
    fixed_atoms: &mut [bool],
    possible_bonds: &mut [bool],
    known_bonds: &[bool],
    fixed_bonds: &mut [bool],
    possible_ring_atoms: &mut [usize],
    possible_ring_bonds: &mut [usize],
) -> Result<(Vec<PotentialStereoInfo>, Vec<u32>), PotentialStereoError> {
    let limit = topology.atoms.len() + topology.bonds.len() + 2;
    let mut ranks = vec![0; topology.atoms.len()];
    let mut result = Vec::new();
    for iteration in 0..limit {
        result.clear();
        ranks = connectivity_ranks(topology, atom_symbols, bond_symbols)?;
        let mut another = update_atoms(
            topology,
            valence,
            &ranks,
            atom_symbols,
            possible_atoms,
            known_atoms,
            fixed_atoms,
            possible_ring_atoms,
            possible_ring_bonds,
            rings,
            allow_nontetrahedral,
            &mut result,
        )?;
        another |= update_bonds(
            topology,
            &ranks,
            bond_symbols,
            possible_atoms,
            possible_bonds,
            known_atoms,
            known_bonds,
            fixed_bonds,
            rings,
            &mut result,
        )?;
        if !another {
            return Ok((result, ranks));
        }
        if iteration + 1 == limit {
            return Err(PotentialStereoError::RefinementDidNotConverge { iterations: limit });
        }
    }
    unreachable!("nonzero refinement limit always returns from the loop")
}

pub fn potential_stereo(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    params: &PotentialStereoParams,
) -> Result<PotentialStereoAssignment, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION findPotentialStereo
    // RDKit✔️❌: std::vector<StereoInfo> findPotentialStereo(ROMol &mol, bool cleanIt,
    // RDKit✔️❌:                                             bool findPossible) {
    // RDKit✔️❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(mol);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️❌:     mol.updatePropertyCache(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::vector<StereoInfo> res = runCleanup(mol, findPossible, cleanIt);
    // RDKit✔️❌:   mol.setProp("_potentialStereo", res, true);
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION findPotentialStereo
    // The caller supplies the already materialized valence and symmetric-ring
    // assignments. The validated total-H helper avoids rescanning topology for
    // each atom and bond; the private rank implementation remains slower than
    // the C++ path.
    validate_inputs(topology, valence, rings)?;
    let mut working = topology.clone();
    let atom_count = working.atoms.len();
    let bond_count = working.bonds.len();
    let mut known_atoms = vec![false; atom_count];
    let mut possible_atoms = vec![false; atom_count];
    let mut atom_symbols = vec![String::new(); atom_count];
    initialize_atoms(
        &mut working,
        valence,
        rings,
        params,
        &mut known_atoms,
        &mut possible_atoms,
        &mut atom_symbols,
    )?;
    let mut known_bonds = vec![false; bond_count];
    let mut possible_bonds = vec![false; bond_count];
    let mut bond_symbols = vec![String::new(); bond_count];
    initialize_bonds(
        &mut working,
        valence,
        rings,
        params,
        &mut known_bonds,
        &mut possible_bonds,
        &mut bond_symbols,
    )?;
    let original_possible_atoms = possible_atoms.clone();
    let original_possible_bonds = possible_bonds.clone();
    let mut possible_ring_atoms = vec![0usize; atom_count];
    let mut possible_ring_bonds = vec![0usize; bond_count];
    let possible_atom_view = if params.clean {
        None
    } else {
        Some(possible_atoms.as_slice())
    };
    let possible_bond_view = if params.clean {
        None
    } else {
        Some(possible_bonds.as_slice())
    };
    flag_ring_stereo(
        &working,
        rings,
        &mut possible_ring_atoms,
        &mut possible_ring_bonds,
        &known_atoms,
        possible_atom_view,
        &known_bonds,
        possible_bond_view,
    );
    let mut fixed_atoms = vec![false; atom_count];
    let mut fixed_bonds = vec![false; bond_count];
    let (mut stereo, mut ranks) = refinement_loop(
        &working,
        valence,
        rings,
        params.allow_nontetrahedral,
        &mut atom_symbols,
        &mut bond_symbols,
        &mut possible_atoms,
        &known_atoms,
        &mut fixed_atoms,
        &mut possible_bonds,
        &known_bonds,
        &mut fixed_bonds,
        &mut possible_ring_atoms,
        &mut possible_ring_bonds,
    )?;
    if params.clean {
        clean_invalid_stereo(
            &mut working,
            &fixed_atoms,
            &known_atoms,
            &fixed_bonds,
            &known_bonds,
        )?;
    }
    if params.flag_possible
        && (possible_atoms != original_possible_atoms || possible_bonds != original_possible_bonds)
    {
        possible_atoms = original_possible_atoms;
        for index in 0..atom_count {
            if !fixed_atoms[index] && known_atoms[index] {
                possible_atoms[index] = true;
                known_atoms[index] = false;
            }
            if possible_atoms[index] {
                atom_symbols[index].push('_');
                atom_symbols[index].push_str(&index.to_string());
            }
        }
        possible_bonds = original_possible_bonds;
        for index in 0..bond_count {
            if !fixed_bonds[index] && known_bonds[index] {
                possible_bonds[index] = true;
                known_bonds[index] = false;
            }
            if possible_bonds[index] {
                bond_symbols[index].push('_');
                bond_symbols[index].push_str(&index.to_string());
            }
        }
        flag_ring_stereo(
            &working,
            rings,
            &mut possible_ring_atoms,
            &mut possible_ring_bonds,
            &known_atoms,
            Some(&possible_atoms),
            &known_bonds,
            Some(&possible_bonds),
        );
        (stereo, ranks) = refinement_loop(
            &working,
            valence,
            rings,
            params.allow_nontetrahedral,
            &mut atom_symbols,
            &mut bond_symbols,
            &mut possible_atoms,
            &known_atoms,
            &mut fixed_atoms,
            &mut possible_bonds,
            &known_bonds,
            &mut fixed_bonds,
            &mut possible_ring_atoms,
            &mut possible_ring_bonds,
        )?;
    }
    let ring_relations = special_ring_relations(&working, valence, rings, &ranks)?;
    working.validate()?;
    Ok(PotentialStereoAssignment {
        stereo,
        atom_ranks: ranks,
        ring_relations,
        cleaned_topology: params.clean.then_some(working),
    })
}

#[cfg(test)]
mod cf_smi_tetra_tests {
    use super::{
        PotentialStereoError, is_potential_tetrahedral, potential_tetrahedral_centers_for_atoms,
        validate_topology_and_valence,
    };
    use crate::{RingFindType, RingInfo, ValenceAssignment};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, TopologyBlock, TopologyValidationError,
    };
    use cosmolkit_types::{BondOrder, Element, Hybridization};

    fn topology(atom_specs: Vec<AtomSpec>, bond_specs: Vec<BondSpec>) -> TopologyBlock {
        let atoms = atom_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Atom::from_spec(AtomId::new(index), spec))
            .collect();
        let bonds = bond_specs
            .into_iter()
            .enumerate()
            .map(|(index, spec)| Bond::from_spec(BondId::new(index), spec))
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    fn edge(begin: usize, end: usize, order: BondOrder) -> BondSpec {
        BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
    }

    fn assignment(
        topology: &TopologyBlock,
        center_explicit_valence: i32,
        center_implicit_hydrogens: i32,
    ) -> ValenceAssignment {
        let mut explicit_valence = vec![0; topology.atoms.len()];
        let mut implicit_hydrogens = vec![0; topology.atoms.len()];
        if !topology.atoms.is_empty() {
            explicit_valence[0] = center_explicit_valence;
            implicit_hydrogens[0] = center_implicit_hydrogens;
        }
        ValenceAssignment {
            explicit_valence,
            implicit_hydrogens,
        }
    }

    fn star(center: AtomSpec, neighbors: Vec<AtomSpec>, bonds: Vec<BondSpec>) -> TopologyBlock {
        let mut atoms = Vec::with_capacity(neighbors.len() + 1);
        atoms.push(center);
        atoms.extend(neighbors);
        topology(atoms, bonds)
    }

    fn retained_rings(
        atom_count: usize,
        bond_count: usize,
        atom_rows: Vec<Vec<usize>>,
        bond_rows: Vec<Vec<usize>>,
    ) -> RingInfo {
        let atom_rings = atom_rows
            .into_iter()
            .map(|row| row.into_iter().map(AtomId::new).collect())
            .collect();
        let bond_rings = bond_rows
            .into_iter()
            .map(|row| row.into_iter().map(BondId::new).collect())
            .collect();
        RingInfo::from_persisted_components(
            true,
            RingFindType::Fast,
            atom_count,
            bond_count,
            atom_rings,
            bond_rings,
            vec![],
            vec![],
            None,
            vec![],
            vec![],
        )
        .unwrap()
    }

    fn scalar(
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        rings: Option<&RingInfo>,
        atom: AtomId,
    ) -> Result<bool, PotentialStereoError> {
        validate_topology_and_valence(topology, valence)?;
        let mut ring_rows_validated = false;
        is_potential_tetrahedral(topology, valence, rings, &mut ring_rows_validated, atom)
    }

    fn assert_source_center(
        topology: &TopologyBlock,
        valence: &ValenceAssignment,
        rings: Option<&RingInfo>,
        expected: bool,
    ) {
        // Fixed outcomes are derived from the pinned FindStereo.cpp branch
        // order. The scalar comparison checks batch routing, not expectations.
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(topology, valence, rings, &[AtomId::new(0)],),
            Ok(vec![expected])
        );
        assert_eq!(
            scalar(topology, valence, rings, AtomId::new(0)),
            Ok(expected)
        );
    }

    #[test]
    fn cf_smi_tetra_degree_boundaries_preserve_requested_order_and_duplicates() {
        let mut atom_specs = Vec::new();
        let mut bond_specs = Vec::new();
        let mut centers = Vec::new();
        for degree in 0..=5 {
            let center = atom_specs.len();
            centers.push(AtomId::new(center));
            atom_specs.push(AtomSpec::new(Element::C));
            for _ in 0..degree {
                let neighbor = atom_specs.len();
                atom_specs.push(AtomSpec::new(Element::F));
                bond_specs.push(edge(center, neighbor, BondOrder::Single));
            }
        }
        let topology = topology(atom_specs, bond_specs);
        let valence = assignment(&topology, 0, 0);
        let requested = [
            centers[5], centers[4], centers[0], centers[1], centers[2], centers[3], centers[4],
        ];
        let expected = [false, true, false, false, false, false, true];
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(&topology, &valence, None, &requested),
            Ok(expected.to_vec())
        );
        for (&atom, &is_potential) in requested.iter().zip(&expected) {
            assert_eq!(scalar(&topology, &valence, None, atom), Ok(is_potential));
        }
    }

    #[test]
    fn cf_smi_tetra_zero_and_directional_dative_bonds_use_source_nonzero_degree() {
        let ligands = || {
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
                AtomSpec::new(Element::I),
            ]
        };
        let ordinary = |last_order| {
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
                edge(0, 4, last_order),
            ]
        };

        let zero = star(
            AtomSpec::new(Element::C),
            ligands(),
            ordinary(BondOrder::Zero),
        );
        assert_source_center(&zero, &assignment(&zero, 0, 0), None, false);

        let dative_out = star(
            AtomSpec::new(Element::C),
            ligands(),
            ordinary(BondOrder::Dative),
        );
        assert_source_center(&dative_out, &assignment(&dative_out, 0, 0), None, false);

        let mut incoming = ordinary(BondOrder::Single);
        incoming[3] = edge(4, 0, BondOrder::Dative);
        let dative_in = star(AtomSpec::new(Element::C), ligands(), incoming);
        assert_source_center(&dative_in, &assignment(&dative_in, 0, 0), None, true);
    }

    #[test]
    fn cf_smi_tetra_phosphorus_and_arsenic_exceptions_precede_hydrogen_branch() {
        let phosphorus = star(
            AtomSpec::new(Element::P),
            vec![AtomSpec::new(Element::F), AtomSpec::new(Element::CL)],
            vec![edge(0, 1, BondOrder::Single), edge(0, 2, BondOrder::Single)],
        );
        assert_source_center(&phosphorus, &assignment(&phosphorus, 2, 0), None, true);

        let arsenic = star(
            AtomSpec::new(Element::AS),
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
                AtomSpec::new(Element::H),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
                edge(0, 4, BondOrder::Zero),
            ],
        );
        assert_source_center(&arsenic, &assignment(&arsenic, 3, 1), None, true);
    }

    #[test]
    fn cf_smi_tetra_explicit_implicit_protium_and_isotopic_hydrogens_follow_source() {
        let heavy_bonds = || {
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ]
        };
        let heavy_neighbors = || {
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
            ]
        };

        let explicit = star(
            AtomSpec::new(Element::C).with_explicit_hydrogens(1),
            heavy_neighbors(),
            heavy_bonds(),
        );
        assert_source_center(&explicit, &assignment(&explicit, 0, 0), None, true);

        let implicit = star(AtomSpec::new(Element::C), heavy_neighbors(), heavy_bonds());
        assert_source_center(&implicit, &assignment(&implicit, 0, 1), None, true);

        let mut protium_bonds = heavy_bonds();
        protium_bonds.push(edge(0, 4, BondOrder::Zero));
        let protium = star(
            AtomSpec::new(Element::C),
            [heavy_neighbors(), vec![AtomSpec::new(Element::H)]].concat(),
            protium_bonds,
        );
        assert_source_center(&protium, &assignment(&protium, 0, 1), None, false);

        let mut isotope_bonds = heavy_bonds();
        isotope_bonds.push(edge(0, 4, BondOrder::Zero));
        let isotopic_h = star(
            AtomSpec::new(Element::C),
            [
                heavy_neighbors(),
                vec![AtomSpec::new(Element::H).with_isotope(2)],
            ]
            .concat(),
            isotope_bonds,
        );
        assert_source_center(&isotopic_h, &assignment(&isotopic_h, 0, 1), None, true);
    }

    #[test]
    fn cf_smi_tetra_sulfur_and_selenium_explicit_valence_and_charge_cases() {
        let neighbors = || {
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
            ]
        };
        let bonds = || {
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ]
        };

        let sulfur_valence_four = star(AtomSpec::new(Element::S), neighbors(), bonds());
        assert_source_center(
            &sulfur_valence_four,
            &assignment(&sulfur_valence_four, 4, 0),
            None,
            true,
        );

        let sulfur_cation = star(
            AtomSpec::new(Element::S).with_formal_charge(1),
            neighbors(),
            bonds(),
        );
        assert_source_center(
            &sulfur_cation,
            &assignment(&sulfur_cation, 3, 0),
            None,
            true,
        );

        let sulfur_neutral = star(AtomSpec::new(Element::S), neighbors(), bonds());
        assert_source_center(
            &sulfur_neutral,
            &assignment(&sulfur_neutral, 3, 0),
            None,
            false,
        );

        let selenium_valence_four = star(AtomSpec::new(Element::SE), neighbors(), bonds());
        assert_source_center(
            &selenium_valence_four,
            &assignment(&selenium_valence_four, 4, 0),
            None,
            true,
        );
    }

    #[test]
    fn cf_smi_tetra_nitrogen_ring_bridgehead_hybridization_and_conjugation_routes() {
        let triangle = |hybridization, conjugated| {
            let center = AtomSpec::new(Element::N).with_hybridization(hybridization);
            let mut first = edge(0, 1, BondOrder::Single);
            first = first.with_conjugated(conjugated);
            let topology = star(
                center,
                vec![
                    AtomSpec::new(Element::C),
                    AtomSpec::new(Element::C),
                    AtomSpec::new(Element::C),
                ],
                vec![
                    first,
                    edge(1, 2, BondOrder::Single),
                    edge(2, 0, BondOrder::Single),
                    edge(0, 3, BondOrder::Single),
                ],
            );
            let rings = retained_rings(4, 4, vec![vec![0, 1, 2]], vec![vec![0, 1, 2]]);
            (topology, rings)
        };

        let (three_ring_nitrogen, three_ring_info) = triangle(Hybridization::Sp3, false);
        assert_source_center(
            &three_ring_nitrogen,
            &assignment(&three_ring_nitrogen, 3, 0),
            Some(&three_ring_info),
            true,
        );

        let (sp2_nitrogen, _) = triangle(Hybridization::Sp2, false);
        assert_source_center(&sp2_nitrogen, &assignment(&sp2_nitrogen, 3, 0), None, false);

        let (conjugated_nitrogen, _) = triangle(Hybridization::Sp3, true);
        assert_source_center(
            &conjugated_nitrogen,
            &assignment(&conjugated_nitrogen, 3, 0),
            None,
            false,
        );

        let acyclic = star(
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp3),
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ],
        );
        let acyclic_rings = RingInfo::new(RingFindType::Fast, 4, 3);
        assert_source_center(
            &acyclic,
            &assignment(&acyclic, 3, 0),
            Some(&acyclic_rings),
            false,
        );

        let bridgehead = star(
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp3),
            vec![
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
                AtomSpec::new(Element::C),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(1, 2, BondOrder::Single),
                edge(2, 0, BondOrder::Single),
                edge(2, 3, BondOrder::Single),
                edge(3, 0, BondOrder::Single),
            ],
        );
        let bridgehead_rings = retained_rings(
            4,
            5,
            vec![vec![0, 1, 2], vec![0, 1, 2, 3]],
            vec![vec![0, 1, 2], vec![0, 1, 3, 4]],
        );
        assert_source_center(
            &bridgehead,
            &assignment(&bridgehead, 3, 0),
            Some(&bridgehead_rings),
            true,
        );
    }

    #[test]
    fn cf_smi_tetra_required_ring_state_errors_only_when_source_reaches_nitrogen_query() {
        let nitrogen = star(
            AtomSpec::new(Element::N).with_hybridization(Hybridization::Sp3),
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
            ],
        );
        let valence = assignment(&nitrogen, 3, 0);
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(&nitrogen, &valence, None, &[AtomId::new(0)]),
            Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring information is unavailable for nitrogen potential-stereo branch",
                row: 0,
                value: 0,
                limit: 0,
            })
        );
        assert_eq!(
            scalar(&nitrogen, &valence, None, AtomId::new(0)),
            Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring information is unavailable for nitrogen potential-stereo branch",
                row: 0,
                value: 0,
                limit: 0,
            })
        );

        let mismatched_rings = RingInfo::new(RingFindType::Fast, 3, 3);
        assert!(matches!(
            potential_tetrahedral_centers_for_atoms(
                &nitrogen,
                &valence,
                Some(&mismatched_rings),
                &[AtomId::new(0)],
            ),
            Err(PotentialStereoError::InvalidRingInfo {
                reason: "atom membership row count mismatch",
                value: 3,
                limit: 4,
                ..
            })
        ));

        let carbon = star(
            AtomSpec::new(Element::C),
            vec![AtomSpec::new(Element::F)],
            vec![edge(0, 1, BondOrder::Single)],
        );
        assert_source_center(&carbon, &assignment(&carbon, 0, 0), None, false);
    }

    #[test]
    fn cf_smi_tetra_rejects_malformed_atom_topology_and_valence_state() {
        let source = star(
            AtomSpec::new(Element::C),
            vec![
                AtomSpec::new(Element::F),
                AtomSpec::new(Element::CL),
                AtomSpec::new(Element::BR),
                AtomSpec::new(Element::I),
            ],
            vec![
                edge(0, 1, BondOrder::Single),
                edge(0, 2, BondOrder::Single),
                edge(0, 3, BondOrder::Single),
                edge(0, 4, BondOrder::Single),
            ],
        );
        let valence = assignment(&source, 4, 0);
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(&source, &valence, None, &[AtomId::new(99)],),
            Err(PotentialStereoError::AtomOutOfRange {
                atom: AtomId::new(99),
                atom_count: 5,
            })
        );
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(
                &source,
                &ValenceAssignment {
                    explicit_valence: vec![4; 4],
                    implicit_hydrogens: vec![0; 5],
                },
                None,
                &[AtomId::new(0)],
            ),
            Err(PotentialStereoError::InvalidValence {
                field: "explicit_valence",
                actual: 4,
                atom_count: 5,
            })
        );
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(
                &source,
                &ValenceAssignment {
                    explicit_valence: vec![4, -1, 1, 1, 1],
                    implicit_hydrogens: vec![0; 5],
                },
                None,
                &[AtomId::new(0)],
            ),
            Err(PotentialStereoError::InvalidValenceValue {
                field: "explicit_valence",
                atom: AtomId::new(1),
                value: -1,
            })
        );

        let mut malformed = source.clone();
        malformed.atoms[0] = malformed.atoms[0].clone().with_id(AtomId::new(1));
        assert_eq!(
            potential_tetrahedral_centers_for_atoms(&malformed, &valence, None, &[AtomId::new(0)],),
            Err(PotentialStereoError::InvalidTopology(
                TopologyValidationError::AtomIdMismatch {
                    position: 0,
                    id: AtomId::new(1),
                }
            ))
        );
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, PropertyValue, PropertyValueKind};
    use cosmolkit_types::Element;
    // FROZEN UINT CONDITION: STRICT_BOOL_VECTOR
    #[test]
    fn uint_cell_strict_bool_vector_potential_stereo() {
        for n in [0_u32, 1, 4294967295] {
            let atom = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_prop("_ringStereochemCand", PropertyValue::UInt(n))
                    .unwrap(),
            );
            let g = TopologyBlock::try_from_parts(vec![atom], vec![], vec![], vec![]).unwrap();
            let before = g.clone();
            let v = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
            let mut cache = [None];
            let mut order = vec![];
            assert_eq!(
                ring_candidate_cached(
                    &g,
                    &v,
                    &RingInfo::new(crate::RingFindType::Fast, 1, 0),
                    &[0],
                    AtomId::new(0),
                    &mut cache,
                    &mut order
                ),
                Err(PotentialStereoError::InvalidPropertyKind {
                    atom: AtomId::new(0),
                    property: "_ringStereochemCand",
                    kind: PropertyValueKind::UInt
                })
            );
            assert_eq!(cache, [None]);
            assert!(order.is_empty());
            assert_eq!(g, before);
        }
    }
}

fn is_potential_bond_with_prefix_check(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    bond: &Bond,
    preserved_prefix: &mut Option<bool>,
) -> Result<bool, PotentialStereoError> {
    is_potential_bond_impl(topology, valence, rings, bond, Some(preserved_prefix))
}

fn is_potential_bond_impl(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    bond: &Bond,
    preserved_prefix: Option<&mut Option<bool>>,
) -> Result<bool, PotentialStereoError> {
    // BEGIN RDKIT CPP FUNCTION isBondPotentialStereoBond
    // RDKit✔️❌: bool isBondPotentialStereoBond(const Bond *bond) {
    // RDKit✔️❌:   PRECONDITION(bond, "bond is null");
    // RDKit✔️❌:   if (bond->getBondType() != Bond::BondType::DOUBLE) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   // at the moment the condition for being a potential stereo bond is that
    // RDKit✔️❌:   // each of the beginning and end neighbors must have at least 2 explicit
    // RDKit✔️❌:   // neighbors but no more than 3 total neighbors.
    // RDKit✔️❌:   // if it's a ring bond, the smallest ring it's in must have at least 8
    // RDKit✔️❌:   // members
    // RDKit✔️❌:   //  (this is common with InChI)
    // RDKit✔️❌:   const auto beginAtom = bond->getBeginAtom();
    // RDKit✔️❌:   auto begDegree = beginAtom->getTotalDegree();
    // RDKit✔️❌:   const auto endAtom = bond->getEndAtom();
    // RDKit✔️❌:   auto endDegree = endAtom->getTotalDegree();
    // RDKit✔️❌:   if (begDegree > 1 && begDegree < 4 && endDegree > 1 && endDegree < 4 &&
    // RDKit✔️❌:       beginAtom->getTotalNumHs(true) < 2 && endAtom->getTotalNumHs(true) < 2) {
    // RDKit✔️❌:     // check rings
    // RDKit✔️❌:     const auto ri = bond->getOwningMol().getRingInfo();
    // RDKit✔️❌:     for (const auto &bring : ri->bondRings()) {
    // RDKit✔️❌:       if (bring.size() < minRingSizeForDoubleBondStereo &&
    // RDKit✔️❌:           std::find(bring.begin(), bring.end(), bond->getIdx()) !=
    // RDKit✔️❌:               bring.end()) {
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     return true;
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION isBondPotentialStereoBond
    if bond.order() != BondOrder::Double {
        return Ok(false);
    }
    let begin_degree = total_degree(topology, valence, bond.begin())?;
    let end_degree = total_degree(topology, valence, bond.end())?;
    if !(2..4).contains(&begin_degree) || !(2..4).contains(&end_degree) {
        return Ok(false);
    }
    if total_hydrogens(topology, valence, bond.begin(), true)? >= 2
        || total_hydrogens(topology, valence, bond.end(), true)? >= 2
    {
        return Ok(false);
    }
    let Some(preserved_prefix) = preserved_prefix else {
        // Native bondRings() is the stored ring vector, including an empty
        // vector before initialization. Membership tables are not consulted.
        return Ok(!rings
            .bond_rings()
            .iter()
            .any(|ring| ring.len() < 8 && ring.contains(&bond.id())));
    };
    if !rings.is_initialized() {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "ring information is not initialized",
            row: 0,
            value: 0,
            limit: 0,
        });
    }
    if rings.bond_row_count() != topology.bonds.len()
        && !*preserved_prefix.get_or_insert_with(|| {
            crate::rings::preserves_appended_terminal_hydrogen_ring_prefix(topology, rings)
        })
    {
        return Err(PotentialStereoError::InvalidRingInfo {
            reason: "bond membership row count mismatch",
            row: 0,
            value: rings.bond_row_count(),
            limit: topology.bonds.len(),
        });
    }
    Ok(!rings
        .bond_ring_sizes(bond.id())
        .into_iter()
        .any(|size| size < 8))
}

fn represented_atropisomer_info(
    topology: &TopologyBlock,
    bond: &Bond,
) -> Result<PotentialStereoInfo, PotentialStereoError> {
    // RDKit✔️✔️:   } else if (bond->getBondType() == Bond::BondType::SINGLE &&
    // RDKit✔️✔️:              (bond->getStereo() == Bond::BondStereo::STEREOATROPCCW ||
    // RDKit✔️✔️:               bond->getStereo() == Bond::BondStereo::STEREOATROPCW)) {
    // RDKit✔️✔️:     if (beginAtom->getDegree() < 2 || endAtom->getDegree() < 2 ||
    // RDKit✔️✔️:         beginAtom->getDegree() > 3 || endAtom->getDegree() > 3) {
    // RDKit✔️✔️:       throw ValueErrorException("invalid atom degree in getStereoInfo(bond)");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     sinfo.type = StereoType::Bond_Atropisomer;
    // RDKit✔️✔️:     sinfo.centeredOn = bond->getIdx();
    // RDKit✔️✔️:     sinfo.controllingAtoms.reserve(4);
    // RDKit✔️✔️:
    // RDKit✔️✔️:     const auto &mol = bond->getOwningMol();
    // RDKit✔️✔️:     for (const auto nbr : mol.atomBonds(beginAtom)) {
    // RDKit✔️✔️:       if (nbr->getIdx() != bond->getIdx()) {
    // RDKit✔️✔️:         sinfo.controllingAtoms.push_back(
    // RDKit✔️✔️:             nbr->getOtherAtomIdx(beginAtom->getIdx()));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (beginAtom->getDegree() == 2) {
    // RDKit✔️✔️:       sinfo.controllingAtoms.push_back(Atom::NOATOM);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto nbr : mol.atomBonds(endAtom)) {
    // RDKit✔️✔️:       if (nbr->getIdx() != bond->getIdx()) {
    // RDKit✔️✔️:         sinfo.controllingAtoms.push_back(
    // RDKit✔️✔️:             nbr->getOtherAtomIdx(endAtom->getIdx()));
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (endAtom->getDegree() == 2) {
    // RDKit✔️✔️:       sinfo.controllingAtoms.push_back(Atom::NOATOM);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:
    // RDKit✔️✔️:     Bond::BondStereo stereo = bond->getStereo();
    // RDKit✔️✔️:     sinfo.specified = Chirality::StereoSpecified::Specified;
    // RDKit✔️✔️:     switch (stereo) {
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCW:
    // RDKit✔️✔️:         sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCW;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       case Bond::BondStereo::STEREOATROPCCW:
    // RDKit✔️✔️:         sinfo.descriptor = Chirality::StereoDescriptor::Bond_AtropCCW;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:
    // RDKit✔️✔️:         UNDER_CONSTRUCTION("unrecognized bond stereo type");
    // RDKit✔️✔️:     }
    // Degree lookup is O(1); both source and detached carriers visit at most
    // three bonds at each end and allocate the same four controlling slots.
    for (endpoint, atom) in [("begin", bond.begin()), ("end", bond.end())] {
        let degree = graph_degree(topology, atom);
        if !(2..=3).contains(&degree) {
            return Err(PotentialStereoError::InvalidBondDegree {
                bond: bond.id(),
                endpoint,
                degree,
            });
        }
    }
    let mut controlling_atoms = Vec::with_capacity(4);
    for atom in [bond.begin(), bond.end()] {
        for neighbor in topology.adjacency.neighbors_of(atom.index()) {
            if neighbor.bond != bond.id() {
                controlling_atoms.push(Some(AtomId::new(neighbor.atom_index)));
            }
        }
        if graph_degree(topology, atom) == 2 {
            controlling_atoms.push(None);
        }
    }
    Ok(PotentialStereoInfo {
        stereo_type: PotentialStereoType::BondAtropisomer,
        specified: PotentialStereoSpecified::Specified,
        centered_on: PotentialStereoCenter::Bond(bond.id()),
        descriptor: match bond.stereo() {
            BondStereo::AtropCw => PotentialStereoDescriptor::BondAtropCw,
            BondStereo::AtropCcw => PotentialStereoDescriptor::BondAtropCcw,
            _ => {
                return Err(PotentialStereoError::InvalidStereoReferences {
                    bond: bond.id(),
                    reason: "represented atropisomer requires CW or CCW stereo",
                });
            }
        },
        permutation: 0,
        controlling_atoms,
    })
}

#[cfg(test)]
mod source_cached_getter_conditions {
    use super::*;
    use cosmolkit_model::{AtomSpec, Element};
    #[test]
    fn no_implicit_and_signed_width_use_source_counts_without_changing_rows() {
        for (no_implicit, implicit, expected) in [
            (true, -1, 2),
            (true, 128, 2),
            (false, 256, 2),
            (false, -256, 2),
            (false, 1, 3),
        ] {
            let topology = TopologyBlock::try_from_parts(
                vec![Atom::from_spec(
                    AtomId::new(0),
                    AtomSpec::new(Element::C)
                        .with_no_implicit(no_implicit)
                        .with_explicit_hydrogens(2),
                )],
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            let valence = ValenceAssignment {
                explicit_valence: vec![2],
                implicit_hydrogens: vec![implicit],
            };
            let before = valence.clone();
            validate_topology_and_valence(&topology, &valence).unwrap();
            assert_eq!(
                total_degree(&topology, &valence, AtomId::new(0)).unwrap(),
                expected
            );
            assert_eq!(valence, before);
        }
        let topology = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        assert_eq!(
            validate_topology_and_valence(
                &topology,
                &ValenceAssignment {
                    explicit_valence: vec![0],
                    implicit_hydrogens: vec![-1]
                }
            ),
            Err(PotentialStereoError::InvalidValenceValue {
                field: "implicit_hydrogens",
                atom: AtomId::new(0),
                value: -1
            })
        );
    }
}

#[cfg(test)]
mod source590_scalar_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, Bond, BondSpec, Element};
    #[test]
    fn source590_scalar_nitrogen_accepts_initialized_sparse_source_cache() {
        let mut atoms = (0..4)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    AtomSpec::new(if i == 0 { Element::N } else { Element::C }),
                )
            })
            .collect::<Vec<_>>();
        atoms[0].set_hybridization(Hybridization::Sp3);
        let topology = TopologyBlock::try_from_parts(
            atoms,
            (1..4)
                .map(|i| {
                    Bond::from_spec(
                        BondId::new(i - 1),
                        BondSpec::new(AtomId::new(0), AtomId::new(i), BondOrder::Single),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap();
        let valence = crate::assign_valence_with_options_for_topology(
            &topology,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let mut rings = crate::RingInfo::new(crate::RingFindType::Sssr, 0, 0);
        rings.reset();
        crate::find_sssr_with_source_outputs_from_parts(
            4,
            &topology.bonds,
            &topology.adjacency,
            &mut rings,
            None,
            None,
            false,
            false,
        )
        .unwrap();
        assert!(rings.is_initialized());
        assert_eq!(rings.atom_row_count(), 0);
        assert!(
            !potential_tetrahedral_center_from_source(
                &topology,
                &valence,
                Some(&rings),
                AtomId::new(0)
            )
            .unwrap()
        );
        assert!(matches!(
            potential_tetrahedral_centers_for_atoms(
                &topology,
                &valence,
                Some(&rings),
                &[AtomId::new(0)]
            ),
            Err(PotentialStereoError::InvalidRingInfo {
                reason: "atom membership row count mismatch",
                ..
            })
        ));
        assert!(matches!(
            potential_tetrahedral_center_from_source(
                &topology,
                &valence,
                Some(&{
                    let mut empty = crate::RingInfo::new(crate::RingFindType::Sssr, 0, 0);
                    empty.reset();
                    empty
                }),
                AtomId::new(0)
            ),
            Err(PotentialStereoError::InvalidRingInfo {
                reason: "ring information is not initialized",
                ..
            })
        ));
    }
}
