// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    borrow::Cow,
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, VecDeque},
};

use crate::{
    RingFindType, RingFindingError, RingInfo, ValenceAssignment, ValenceError, ValenceModel,
    fast_find_rings_from_parts, find_sssr_from_parts,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, Bond, BondId, QueryStateError, QueryStateRef, StereoGroupKind,
    TemplateAttachmentOrderError, TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

#[derive(Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum CanonicalRankError {
    #[error("hanoi scratch is too small")]
    HanoiScratchTooSmall,
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[error(
        "atom {atom_index} property {property} unsigned value {value} causes positive_overflow converting UInt to signed int"
    )]
    UnsignedRankOverflow {
        atom_index: usize,
        property: &'static str,
        value: u32,
    },
    #[error("atom {atom_index} property {property} has invalid kind {kind:?}")]
    InvalidPropertyKind {
        atom_index: usize,
        property: &'static str,
        kind: cosmolkit_model::PropertyValueKind,
    },
    #[error("component atom index {atom_index} is outside topology atom count {atom_count}")]
    ComponentAtomOutOfRange {
        atom_index: usize,
        atom_count: usize,
    },
    #[error("component atom index {atom_index} occurs more than once")]
    DuplicateComponentAtom { atom_index: usize },
    #[error("atom mask length {actual} does not match topology atom count {expected}")]
    AtomMaskLength { expected: usize, actual: usize },
    #[error("bond mask length {actual} does not match topology bond count {expected}")]
    BondMaskLength { expected: usize, actual: usize },
    #[error("atom symbol length {actual} does not match topology atom count {expected}")]
    AtomSymbolLength { expected: usize, actual: usize },
    #[error("bond symbol length {actual} does not match topology bond count {expected}")]
    BondSymbolLength { expected: usize, actual: usize },
    #[error(
        "prepared valence dimensions explicit={explicit_len}, implicit={implicit_len} do not match {atom_count} atoms"
    )]
    PreparedValenceLength {
        atom_count: usize,
        explicit_len: usize,
        implicit_len: usize,
    },
    #[error("prepared valence is invalid for atom {atom_index}")]
    PreparedValenceInvalid { atom_index: usize },
    #[error(
        "prepared ring dimensions atoms={actual_atoms}, bonds={actual_bonds} do not match topology atoms={expected_atoms}, bonds={expected_bonds}"
    )]
    PreparedRingLength {
        expected_atoms: usize,
        actual_atoms: usize,
        expected_bonds: usize,
        actual_bonds: usize,
    },
    #[error("incomplete RDKit canonical-rank port for {branch}: {reason}")]
    ProtocolDebt {
        branch: &'static str,
        reason: &'static str,
    },
    #[error("template attachment remap failed for carrier atom {carrier}: {source}")]
    TemplateAttachmentRemap {
        carrier: AtomId,
        source: TemplateAttachmentOrderError,
    },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

#[derive(Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum KekulizeError {
    #[error("invalid topology: {0}")]
    InvalidTopology(#[from] TopologyValidationError),
    #[error("atom selection length {actual} does not match topology atom count {expected}")]
    AtomSelectionLength { expected: usize, actual: usize },
    #[error("bond selection length {actual} does not match topology bond count {expected}")]
    BondSelectionLength { expected: usize, actual: usize },
    #[error("candidate atom {atom} is outside {atom_count} topology atoms")]
    CandidateAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("candidate atom {atom} occurs more than once")]
    DuplicateCandidateAtom { atom: AtomId },
    #[error("{field} length {actual} does not match expected length {expected}")]
    MatchingStateLength {
        field: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error("matching state contains duplicate done atom {atom}")]
    DuplicateDoneAtom { atom: AtomId },
    #[error("matching state references missing bond between atoms {begin} and {end}")]
    MissingCandidateBond { begin: AtomId, end: AtomId },
    #[error("matching backtrack anchor {atom} is absent from the completed-atom list")]
    MissingBacktrackAnchor { atom: AtomId },
    #[error(
        "cannot enumerate {questions} dummy-atom questions with the source {bit_width}-bit subset counter"
    )]
    QuestionSubsetOverflow { questions: usize, bit_width: u32 },
    #[error("aromatic atom {atom} is not in a ring")]
    AromaticAtomOutsideRing { atom: AtomId },
    #[error("could not kekulize molecule; remaining atoms: {problem_atoms:?}")]
    NotKekulizable { problem_atoms: Vec<AtomId> },
    #[error("kekulization changed total valence for atom {atom} from {before} to {after}")]
    PostconditionValenceMismatch {
        atom: AtomId,
        before: i32,
        after: i32,
    },
    #[error("unsupported query state on bond {bond}: {detail}")]
    UnsupportedQueryState { bond: BondId, detail: &'static str },
    #[error(transparent)]
    InvalidQueryState(#[from] QueryStateError),
    #[error("integer overflow while computing {field} for atom {atom}")]
    IntegerOverflow { atom: AtomId, field: &'static str },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    CanonicalRank(#[from] CanonicalRankError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KekulizeParams {
    pub mark_atoms_bonds: bool,
    pub canonical: bool,
    pub max_backtracks: u32,
}

impl Default for KekulizeParams {
    fn default() -> Self {
        Self {
            mark_atoms_bonds: true,
            canonical: true,
            max_backtracks: 100,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct KekulizeAssignment {
    pub topology: TopologyBlock,
    /// Actual source scalar state after selected implicit updates and neutral
    /// aromatic N/P normalization. Unselected rows retain their input state.
    /// Nonaromatic returns retain selected writes; atoms-none retains input.
    /// `refreshed_valence_atoms` identifies explicit N/P normalization only.
    pub final_valence: Option<ValenceAssignment>,
    pub refreshed_valence_atoms: Vec<AtomId>,
    /// Final ring-state transport for the Kekulize ring-update calling
    /// convention: `None` means the supplied caller state is unchanged
    /// (absent, already reset, or an initialized assignment that was
    /// borrowed without modification); `Some(initialized)` moves the newly
    /// acquired final state; `Some(uninitialized)` replaces a previously
    /// initialized Other assignment that canonical ranking reset.
    pub ring_update: Option<RingInfo>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum KekulizeAttempt {
    Applied(KekulizeAssignment),
    NotKekulizable {
        topology: TopologyBlock,
        problem_atoms: Vec<AtomId>,
    },
}

#[derive(Debug)]
struct PreparedKekulizeSelection {
    atoms_in_play: Vec<bool>,
    bonds_in_play: Vec<bool>,
    original_total_valences: Vec<i32>,
    dummy_atoms: Vec<bool>,
    rings: RingInfo,
    candidate_atom_rings: Vec<Vec<AtomId>>,
    candidate_bond_rings: Vec<Vec<BondId>>,
    found_aromatic: bool,
    valence: ValenceAssignment,
}

/// Selection-stage results that do not depend on ring state: dimension
/// validation, query-type filtering, valence work, dummy detection and the
/// wedge scan, in source order and before any ring acquisition.
#[derive(Debug)]
struct PreparedKekulizeCore {
    atoms_in_play: Vec<bool>,
    bonds_in_play: Vec<bool>,
    original_total_valences: Vec<i32>,
    dummy_atoms: Vec<bool>,
    wedged_atoms: Vec<bool>,
    found_aromatic: bool,
    valence: ValenceAssignment,
}

#[derive(Debug)]
struct CandidateState {
    topology: TopologyBlock,
    double_bond_candidates: Vec<bool>,
    questions: Vec<AtomId>,
    done: Vec<AtomId>,
}

#[derive(Debug)]
struct MatchingState {
    topology: TopologyBlock,
    succeeded: bool,
    double_bond_candidates: Vec<bool>,
    double_bonds_added: Vec<bool>,
    done: Vec<AtomId>,
    problem_atoms: Vec<AtomId>,
    backtracks: u32,
}

#[derive(Debug)]
struct FusedKekulizeState {
    topology: TopologyBlock,
    succeeded: bool,
    problem_atoms: Vec<AtomId>,
}

#[derive(Debug)]
struct CandidateAttempt {
    double_bond_candidates: Vec<bool>,
    questions: Vec<AtomId>,
    done: Vec<AtomId>,
}
#[derive(Debug)]
struct MatchingAttempt {
    succeeded: bool,
    double_bond_candidates: Vec<bool>,
    double_bonds_added: Vec<bool>,
    done: Vec<AtomId>,
    problem_atoms: Vec<AtomId>,
    backtracks: u32,
}
#[derive(Debug)]
struct FusedAttempt {
    succeeded: bool,
    problem_atoms: Vec<AtomId>,
}
#[derive(Debug)]
struct PreparedKekulizeInputs {
    atoms_in_play: Vec<bool>,
    bonds_in_play: Vec<bool>,
    original_total_valences: Vec<i32>,
    dummy_atoms: Vec<bool>,
    wedged_atoms: Vec<bool>,
    found_aromatic: bool,
}

#[derive(Debug)]
struct QuestionEnumerator {
    questions: Vec<AtomId>,
    // Packed little-endian binary counter; shifts are confined to 0..64.
    state: Vec<u64>,
    done: bool,
}

fn atom_is_aromatic_for_kekulize(topology: &TopologyBlock, atom: AtomId) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: isAromaticAtom
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
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: isAromaticAtom
    // Behavior: a validated detached topology supplies the owning neighborhood;
    // preserve the atom flag OR incident flag OR incident aromatic-order test.
    // Cost: O(degree) indexed adjacency traversal with early return, no allocation
    // or clone, matching the source neighborhood scan and flag short circuit.
    if topology.atoms[atom.index()].is_aromatic() {
        return true;
    }
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| {
            let bond = &topology.bonds[neighbor.bond.index()];
            bond.is_aromatic() || bond.order() == BondOrder::Aromatic
        })
}

fn selected_bond_has_type_query(
    bond: &Bond,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<bool, KekulizeError> {
    // RDKit✔️✔️: inline bool hasBondTypeQuery(const Bond &bond) {
    // RDKit✔️✔️:   if (!bond.hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return hasBondTypeQuery(*bond.getQuery());
    // RDKit✔️✔️: }
    // Behavior review: typed Explicit origin is the source dynamic query-bond
    // test; CarrierDerived and absent state are ordinary bonds. The shared
    // QueryOps traversal owns complete recursive bond-type classification.
    // Complexity review: two O(1) indexed origin/predicate accesses followed
    // by the source-shaped O(query nodes) traversal, without allocation.
    let Some(state) = query_state else {
        return Ok(false);
    };
    if !state.bond_has_query(bond.id()) {
        return Ok(false);
    }
    Ok(crate::query_ops::bond_predicate_has_type_query(
        state.bond_predicate(bond.id()),
    ))
}

fn checked_total_valence(
    valence: &ValenceAssignment,
    atom: &cosmolkit_model::Atom,
) -> Result<i32, KekulizeError> {
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Unique O(1) cached getters retain signed storage, initialization errors
    // and the source NoImplicit branch; no topology calculation occurs here.
    Ok(crate::valence::cached_total_valence(atom, valence)?)
}

/// Cold selection assembly retained by the owning selection tests: the same
/// core stages plus the fresh-SSSR acquisition the rings=None engine path
/// performs through the transition owner. Production callers use
/// prepare_kekulize_core plus the transition owner directly.
#[cfg(test)]
fn prepare_kekulize_selection(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    query_state: Option<QueryStateRef<'_>>,
) -> Result<PreparedKekulizeSelection, KekulizeError> {
    let core = prepare_kekulize_core(topology, atoms_in_play, bonds_in_play, query_state, None)?;
    let rings = if core.found_aromatic {
        // BEGIN RDKIT CPP FUNCTION KekulizeFragment ring selection
        // RDKit✔️✔️: VECT_INT_VECT allringsSSSR;
        // RDKit✔️✔️: if (!mol.getRingInfo()->isInitialized()) {
        // RDKit✔️✔️:   MolOps::findSSSR(mol, allringsSSSR);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: const VECT_INT_VECT &allrings =
        // RDKit✔️✔️:     allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
        // Detached TopologyBlock has no RingInfo input, so this owner follows
        // the source's uninitialized-state SSSR branch after canonical rank's
        // temporary fast-ring preparation has been reset.
        find_sssr_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)?
    } else {
        RingInfo::new(
            RingFindType::Fast,
            topology.atoms.len(),
            topology.bonds.len(),
        )
    };
    let (candidate_atom_rings, candidate_bond_rings) = if core.found_aromatic {
        collect_kekulize_candidate_rings(
            topology,
            &rings,
            &core.atoms_in_play,
            &core.bonds_in_play,
            &core.dummy_atoms,
            &core.wedged_atoms,
        )?
    } else {
        (Vec::new(), Vec::new())
    };
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment selection
    Ok(PreparedKekulizeSelection {
        atoms_in_play: core.atoms_in_play,
        bonds_in_play: core.bonds_in_play,
        original_total_valences: core.original_total_valences,
        dummy_atoms: core.dummy_atoms,
        rings,
        candidate_atom_rings,
        candidate_bond_rings,
        found_aromatic: core.found_aromatic,
        valence: core.valence,
    })
}

fn prepare_kekulize_core(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    query_state: Option<QueryStateRef<'_>>,
    source_valence: Option<&ValenceAssignment>,
) -> Result<PreparedKekulizeCore, KekulizeError> {
    let mut valence = source_valence
        .cloned()
        .unwrap_or_else(|| ValenceAssignment {
            explicit_valence: vec![-1; topology.atoms.len()],
            implicit_hydrogens: vec![-1; topology.atoms.len()],
        });
    let prepared = prepare_kekulize_inputs(
        topology,
        atoms_in_play,
        bonds_in_play,
        query_state,
        &mut valence,
    )?;
    // Retain the old private atoms-none scratch result for existing callers.
    if !atoms_in_play.iter().any(|selected| *selected) {
        valence = ValenceAssignment {
            explicit_valence: vec![0; topology.atoms.len()],
            implicit_hydrogens: vec![0; topology.atoms.len()],
        };
    }
    Ok(PreparedKekulizeCore {
        atoms_in_play: prepared.atoms_in_play,
        bonds_in_play: prepared.bonds_in_play,
        original_total_valences: prepared.original_total_valences,
        dummy_atoms: prepared.dummy_atoms,
        wedged_atoms: prepared.wedged_atoms,
        found_aromatic: prepared.found_aromatic,
        valence,
    })
}

fn prepare_kekulize_inputs(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    query_state: Option<QueryStateRef<'_>>,
    valence: &mut ValenceAssignment,
) -> Result<PreparedKekulizeInputs, KekulizeError> {
    // Existing detached graph/query/width validation order is preserved.
    topology.validate()?;
    if atoms_in_play.len() != topology.atoms.len() {
        return Err(KekulizeError::AtomSelectionLength {
            expected: topology.atoms.len(),
            actual: atoms_in_play.len(),
        });
    }
    if bonds_in_play.len() != topology.bonds.len() {
        return Err(KekulizeError::BondSelectionLength {
            expected: topology.bonds.len(),
            actual: bonds_in_play.len(),
        });
    }
    if !atoms_in_play.iter().any(|selected| *selected) {
        return Ok(PreparedKekulizeInputs {
            atoms_in_play: atoms_in_play.to_vec(),
            bonds_in_play: bonds_in_play.to_vec(),
            original_total_valences: vec![0; topology.atoms.len()],
            dummy_atoms: vec![false; topology.atoms.len()],
            wedged_atoms: vec![false; topology.atoms.len()],
            found_aromatic: false,
        });
    }

    let mut selected_bonds = bonds_in_play.to_vec();
    let mut found_aromatic = false;
    for bond in &topology.bonds {
        if !selected_bonds[bond.id().index()] {
            continue;
        }
        if selected_bond_has_type_query(bond, query_state)? {
            selected_bonds[bond.id().index()] = false;
            continue;
        }
        if bond.is_aromatic() {
            found_aromatic = true;
        }
    }

    for (field, actual) in [
        ("explicit valence", valence.explicit_valence.len()),
        ("implicit valence", valence.implicit_hydrogens.len()),
    ] {
        if actual != topology.atoms.len() {
            return Err(KekulizeError::MatchingStateLength {
                field,
                expected: topology.atoms.len(),
                actual,
            });
        }
    }
    let mut original_total_valences = vec![0; topology.atoms.len()];
    let mut dummy_atoms = vec![false; topology.atoms.len()];
    for (atom_idx, selected) in atoms_in_play.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom_id = AtomId::new(atom_idx);
        crate::valence::source_calc_implicit_cache_row(
            &topology.atoms,
            &topology.bonds,
            &topology.adjacency,
            atom_id,
            &mut valence.explicit_valence[atom_idx],
            &mut valence.implicit_hydrogens[atom_idx],
            false,
        )?;
        original_total_valences[atom_idx] =
            checked_total_valence(valence, &topology.atoms[atom_idx])?;
        found_aromatic |= atom_is_aromatic_for_kekulize(topology, atom_id);
        dummy_atoms[atom_idx] = topology.atoms[atom_idx].atomic_number() == 0;
    }

    let mut wedged_atoms = vec![false; topology.atoms.len()];
    for bond in &topology.bonds {
        if selected_bonds[bond.id().index()]
            && matches!(
                bond.direction(),
                BondDirection::BeginWedge | BondDirection::BeginDash
            )
        {
            wedged_atoms[bond.begin().index()] = true;
        }
    }
    let _ = &selected_bonds;
    Ok(PreparedKekulizeInputs {
        atoms_in_play: atoms_in_play.to_vec(),
        bonds_in_play: selected_bonds,
        original_total_valences,
        dummy_atoms,
        wedged_atoms,
        found_aromatic,
    })
}

/// Copy the source-defined candidate ring rows from one borrowed assignment.
///
/// Owner of the KekulizeFragment candidate rotation/filter/conversion over a
/// borrowed ring state: the caller supplies the rows it wants read (freshly
/// acquired SSSR in the cold production caller, or an initialized caller
/// assignment); only the source-defined candidate row copies are produced.
fn collect_kekulize_candidate_rings(
    topology: &TopologyBlock,
    rings: &RingInfo,
    atoms_in_play: &[bool],
    selected_bonds: &[bool],
    dummy_atoms: &[bool],
    wedged_atoms: &[bool],
) -> Result<(Vec<Vec<AtomId>>, Vec<Vec<BondId>>), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: Kekulize.cpp:658-717 candidate rotation, all-row conversion, paired filter (complete)
    // RDKit✔️✔️:     std::deque<INT_VECT> tmpRings;
    // RDKit✔️✔️:     auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
    // RDKit✔️✔️:     bool ringOk = false;
    // RDKit✔️✔️:     for (auto ai : ring) {
    // RDKit✔️✔️:       if (!atomsToUse[ai]) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!dummyAts[ai]) {
    // RDKit✔️✔️:         ringOk = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return ringOk;
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   // we can't just copy the rings over: we're going to rearrange them so that
    // RDKit✔️✔️:   // we try to favor starting the traversal of any ring from an atom that is
    // RDKit✔️✔️:   // at the end of a wedged ring bond. This is part of our attempt to avoid
    // RDKit✔️✔️:   // assigning double bonds to bonds with wedging
    // RDKit✔️✔️:   for (const auto &ring : allrings) {
    // RDKit✔️✔️:     if (containsNonDummy(ring)) {
    // RDKit✔️✔️:       unsigned int startPos = 0;
    // RDKit✔️✔️:       bool hasWedge = false;
    // RDKit✔️✔️:       for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit✔️✔️:         if (wedgedAtoms[ring[ri]]) {
    // RDKit✔️✔️:           startPos = ri;
    // RDKit✔️✔️:           hasWedge = true;
    // RDKit✔️✔️:           break;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       INT_VECT nring(ring.size());
    // RDKit✔️✔️:       for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit✔️✔️:         nring[ri] = ring.at((ri + startPos) % ring.size());
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!hasWedge) {
    // RDKit✔️✔️:         tmpRings.push_back(nring);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         tmpRings.push_front(nring);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   VECT_INT_VECT arings;
    // RDKit✔️✔️:   arings.reserve(allrings.size());
    // RDKit✔️✔️:   arings.insert(arings.end(), tmpRings.begin(), tmpRings.end());
    // RDKit✔️✔️:   VECT_INT_VECT allbrings;
    // RDKit✔️✔️:   RingUtils::convertToBonds(arings, allbrings, mol);
    // RDKit✔️✔️:   VECT_INT_VECT brings;
    // RDKit✔️✔️:   brings.reserve(allbrings.size());
    // RDKit✔️✔️:   auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
    // RDKit✔️✔️:     return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
    // RDKit✔️✔️:       return bondsToUse[bi];
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   VECT_INT_VECT aringsRemaining;
    // RDKit✔️✔️:   aringsRemaining.reserve(arings.size());
    // RDKit✔️✔️:   for (unsigned i = 0; i < allbrings.size(); ++i) {
    // RDKit✔️✔️:     if (copyBondRingsWithinFragment(allbrings[i])) {
    // RDKit✔️✔️:         brings.push_back(allbrings[i]);
    // RDKit✔️✔️:         aringsRemaining.push_back(arings[i]);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     arings = std::move(aringsRemaining);
    // END RDKIT CPP FUNCTION: Kekulize.cpp:658-717 candidate rotation, all-row conversion, paired filter (complete)
    // Behavior review: the source queues ONLY rotated ATOM rows selected by
    // containsNonDummy (all ring atoms selected + at least one non-dummy),
    // with wedged rows at the deque FRONT; it then derives EVERY candidate
    // bond row from consecutive graph atom pairs plus the closing edge
    // (RingUtils::convertToBonds) and only afterwards filters the DERIVED
    // rows by bondsToUse. Stored bond rows never supply conversion output
    // or determine the filter; a missing graph edge is detected BEFORE the
    // bond filter via the existing typed RingFinding cause.
    // Complexity review: one clone per selected atom row (the source copies
    // nring identically), ONE reserved output per derived row, no whole
    // RingInfo clone, no index-staging vector, no duplicate traversal.
    let mut candidate_rings = VecDeque::new();
    // RDKit✔️✔️: auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
    // RDKit✔️✔️:   bool ringOk = false;
    // RDKit✔️✔️:   for (auto ai : ring) {
    // RDKit✔️✔️:     if (!atomsToUse[ai]) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dummyAts[ai]) {
    // RDKit✔️✔️:       ringOk = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return ringOk;
    // RDKit✔️✔️: };
    // RDKit✔️✔️: auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
    // RDKit✔️✔️:   return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
    // RDKit✔️✔️:     return bondsToUse[bi];
    // RDKit✔️✔️:   });
    // RDKit✔️✔️: };
    // RDKit✔️✔️: // we can't just copy the rings over: we're going to rearrange them so that
    // RDKit✔️✔️: // we try to favor starting the traversal of any ring from an atom that is
    // RDKit✔️✔️: // at the end of a wedged ring bond. This is part of our attempt to avoid
    // RDKit✔️✔️: // assigning double bonds to bonds with wedging
    for atom_ring in rings.atom_rings() {
        let atoms_are_selected = atom_ring.iter().all(|atom| atoms_in_play[atom.index()]);
        let contains_non_dummy = atom_ring.iter().any(|atom| !dummy_atoms[atom.index()]);
        if atoms_are_selected && contains_non_dummy {
            // RDKit✔️✔️:         unsigned int startPos = 0;
            // RDKit✔️✔️:         bool hasWedge = false;
            // RDKit✔️✔️:         for (auto ri = 0u; ri < ring.size(); ++ri) {
            // RDKit✔️✔️:           if (wedgedAtoms[ring[ri]]) {
            // RDKit✔️✔️:             startPos = ri;
            // RDKit✔️✔️:             hasWedge = true;
            // RDKit✔️✔️:             break;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            let wedge_start = atom_ring.iter().position(|atom| wedged_atoms[atom.index()]);
            // RDKit✔️✔️:         INT_VECT nring(ring.size());
            // RDKit✔️✔️:         for (auto ri = 0u; ri < ring.size(); ++ri) {
            // RDKit✔️✔️:           nring[ri] = ring.at((ri + startPos) % ring.size());
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         if (!hasWedge) {
            // RDKit✔️✔️:           tmpRings.push_back(nring);
            // RDKit✔️✔️:         } else {
            // RDKit✔️✔️:           tmpRings.push_front(nring);
            // RDKit✔️✔️:         }
            let start = wedge_start.unwrap_or(0);
            let mut rotated_atoms = atom_ring.clone();
            rotated_atoms.rotate_left(start);
            if wedge_start.is_some() {
                candidate_rings.push_front(rotated_atoms);
            } else {
                candidate_rings.push_back(rotated_atoms);
            }
        }
    }
    // Behavior review (actual three phases): phase 1 rotates each
    // containsNonDummy-surviving atom row by clone+rotate_left (the source
    // copy-constructs nring) and places wedged rows at the deque FRONT;
    // phase 2 materializes ALL rotated rows into ONE owned ordered Vec
    // (rows MOVED out of the deque, matching arings.insert of tmpRings)
    // and derives EVERY bond row from the graph BEFORE any bond filtering;
    // phase 3 filters the completed derived pairs by bondsToUse exactly as
    // copyBondRingsWithinFragment over allbrings and MOVES surviving paired
    // rows into the result (the source copies them into brings/
    // aringsRemaining). No partial result is published; existing
    // RingFinding errors propagate by ?.
    // Allocation/complexity review (actual): phase 1 pays one row clone +
    // one rotate per surviving row (O(row) each); phase 2 pays one reserved
    // arings Vec + one reserved allbrings Vec, with per-edge neighbor-scan
    // conversion cost O(sum of degrees over consecutive pairs, including
    // the closing edge) — NOT O(1) per edge; phase 3 consumes both vectors
    // by value with an O(members) mask scan. These are actual costs of
    // this owner; no whole-Kekulize performance acceptance is claimed.
    let count = candidate_rings.len();
    let mut arings: Vec<Vec<AtomId>> = Vec::with_capacity(count);
    while let Some(atom_row) = candidate_rings.pop_front() {
        arings.push(atom_row);
    }
    let mut allbrings: Vec<Vec<BondId>> = Vec::with_capacity(arings.len());
    for atom_row in &arings {
        allbrings.push(crate::rings::ring_atom_ids_to_bond_ids(topology, atom_row)?);
    }
    // Phase 3 (source copyBondRingsWithinFragment over allbrings): the
    // paired derived-bond filter consumes both completed vectors by value;
    // surviving paired rows are MOVED into the result.
    let mut candidate_atom_rings = Vec::new();
    let mut candidate_bond_rings = Vec::new();
    for (atom_row, bond_row) in arings.into_iter().zip(allbrings) {
        if bond_row.iter().all(|bond| selected_bonds[bond.index()]) {
            candidate_atom_rings.push(atom_row);
            candidate_bond_rings.push(bond_row);
        }
    }
    Ok((candidate_atom_rings, candidate_bond_rings))
}

/// Test-only observation points at the actual rank/SSSR/borrow sites.
/// Baselines are captured at these sites (the counters are never reset);
/// production builds contain no counters, probes or exported hooks.
#[cfg(test)]
mod ring_transport_probe {
    use std::cell::Cell;

    thread_local! {
        static RANK_CALLS: Cell<u64> = const { Cell::new(0) };
        static SSSR_CALLS: Cell<u64> = const { Cell::new(0) };
        static BORROW_SITES: Cell<Vec<usize>> = const { Cell::new(Vec::new()) };
    }

    pub fn record_rank() {
        RANK_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub fn record_sssr() {
        SSSR_CALLS.with(|calls| calls.set(calls.get() + 1));
    }

    pub fn record_borrow(rings: &crate::rings::RingInfo) {
        BORROW_SITES.with(|sites| {
            let mut recorded = sites.take();
            recorded.push(rings as *const _ as usize);
            sites.set(recorded);
        });
    }

    pub fn rank_calls() -> u64 {
        RANK_CALLS.with(Cell::get)
    }

    pub fn sssr_calls() -> u64 {
        SSSR_CALLS.with(Cell::get)
    }

    pub fn borrow_sites() -> Vec<usize> {
        BORROW_SITES.with(|sites| {
            let recorded = sites.take();
            sites.set(recorded.clone());
            recorded
        })
    }

    /// Observation at the actual permanent findSSSR completion: the outer
    /// atom/bond row slice pointers and lengths of the one acquired buffer.
    /// The ranking helper's temporary fast discovery never reaches here.
    thread_local! {
        static ACQUIRED_BUFFERS: Cell<Vec<(usize, usize, usize, usize)>> =
            const { Cell::new(Vec::new()) };
    }

    pub fn record_acquired(rings: &crate::rings::RingInfo) {
        ACQUIRED_BUFFERS.with(|buffers| {
            let mut recorded = buffers.take();
            recorded.push((
                rings.atom_rings().as_ptr() as usize,
                rings.atom_rings().len(),
                rings.bond_rings().as_ptr() as usize,
                rings.bond_rings().len(),
            ));
            buffers.set(recorded);
        });
    }

    pub fn acquired_buffers() -> Vec<(usize, usize, usize, usize)> {
        ACQUIRED_BUFFERS.with(|buffers| {
            let recorded = buffers.take();
            buffers.set(recorded.clone());
            recorded
        })
    }
}

/// The single permanent SSSR acquisition entry used by both source guards.
/// Under cfg(test) the completion of the real finder is observed once here,
/// capturing the acquired buffer's outer row slice pointers and lengths for
/// the moved-output proof; production builds contain no probe.
fn acquire_sssr_rows(topology: &TopologyBlock) -> Result<RingInfo, KekulizeError> {
    let acquired =
        find_sssr_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)?;
    #[cfg(test)]
    ring_transport_probe::record_acquired(&acquired);
    Ok(acquired)
}

/// Effective ring rows available to the KekulizeFragment dispatch, tracked
/// without cloning a supplied assignment: either the caller's initialized
/// state is borrowed unchanged, a fresh acquisition is owned, or no
/// initialized rows exist for the dispatch.
enum KekulizeRingRows<'a> {
    /// Initialized caller assignment, borrowed and never modified.
    Borrowed(&'a RingInfo),
    Live(&'a mut RingInfo),
    /// Fresh SSSR rows moved from the finder (source findSSSR install).
    Acquired(RingInfo),
    /// No initialized rows exist (absent/reset caller state with no
    /// acquisition); the independent marking guard may still acquire.
    Uninitialized,
}

impl KekulizeRingRows<'_> {
    fn as_ring_info(&self) -> Option<&RingInfo> {
        match self {
            KekulizeRingRows::Borrowed(rings) => Some(rings),
            KekulizeRingRows::Live(rings) => rings.is_initialized().then_some(&**rings),
            KekulizeRingRows::Acquired(rings) => Some(rings),
            KekulizeRingRows::Uninitialized => None,
        }
    }
}

/// Validate the membership dimensions of an initialized borrowed assignment
/// at a site that actually consumes its rows or membership tables. The
/// category is the existing PreparedRingError vocabulary surfaced through
/// KekulizeError::CanonicalRank; reset storage (zero dimensions) is never
/// validated here — it is reacquired instead.
fn validate_consumed_ring_dimensions(
    topology: &TopologyBlock,
    rings: &RingInfo,
) -> Result<(), KekulizeError> {
    if rings.atom_row_count() != topology.atoms.len()
        || rings.bond_row_count() != topology.bonds.len()
    {
        return Err(KekulizeError::CanonicalRank(
            CanonicalRankError::PreparedRingLength {
                expected_atoms: topology.atoms.len(),
                actual_atoms: rings.atom_row_count(),
                expected_bonds: topology.bonds.len(),
                actual_bonds: rings.bond_row_count(),
            },
        ));
    }
    Ok(())
}

/// Run the KekulizeFragment ring-state transitions in source order over one
/// borrowed-or-owned state: optional canonical ranking (new_canon.cpp
/// rankFragmentAtoms guard and post-ranking reset) followed by the candidate
/// acquisition restricted to effective bondsToUse.any(). Returns the rows the
/// fused dispatch and marking read, plus the caller-visible ring update.
impl KekulizeRingRows<'_> {
    fn install(&mut self, next: RingInfo) {
        match self {
            Self::Live(rings) => **rings = next,
            _ => *self = Self::Acquired(next),
        }
    }
    fn reset(&mut self) {
        match self {
            Self::Live(rings) => rings.reset(),
            _ => *self = Self::Uninitialized,
        }
    }
    fn into_update(self, replaced_by_reset: bool) -> Option<RingInfo> {
        match self {
            Self::Acquired(rings) => Some(rings),
            Self::Live(_) => None,
            _ if replaced_by_reset => {
                let mut reset = RingInfo::new(RingFindType::OtherOrUnknown, 0, 0);
                reset.reset();
                Some(reset)
            }
            _ => None,
        }
    }
}

fn kekulize_ring_state_transition<'a>(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    canonical: bool,
    acquisition_required: bool,
    rings: Option<&'a RingInfo>,
) -> Result<(KekulizeRingRows<'a>, bool, Vec<usize>), KekulizeError> {
    let mut rows = match rings {
        Some(rings) if rings.is_initialized() => KekulizeRingRows::Borrowed(rings),
        _ => KekulizeRingRows::Uninitialized,
    };
    let (reset, ranks) = kekulize_ring_state_transition_mut(
        topology,
        valence,
        atoms_in_play,
        bonds_in_play,
        canonical,
        acquisition_required,
        &mut rows,
    )?;
    Ok((rows, reset, ranks))
}

fn kekulize_ring_state_transition_mut(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    canonical: bool,
    acquisition_required: bool,
    rows: &mut KekulizeRingRows<'_>,
) -> Result<(bool, Vec<usize>), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankFragmentAtoms (2026.03.6 complete)
    // RDKit❗❌: void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
    // RDKit❗❌:                        const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit❗❌:                        const boost::dynamic_bitset<> &bondsInPlay,
    // RDKit❗❌:                        const std::vector<std::string> *atomSymbols,
    // RDKit❗❌:                        const std::vector<std::string> *bondSymbols,
    // RDKit❗❌:                        bool breakTies, bool includeChirality,
    // RDKit❗❌:                        bool includeIsotopes, bool includeAtomMaps,
    // RDKit❗❌:                        bool includeChiralPresence, bool includeRingStereo) {
    // RDKit❗❌:   PRECONDITION(atomsInPlay.size() == mol.getNumAtoms(), "bad atomsInPlay size");
    // RDKit❗❌:   PRECONDITION(bondsInPlay.size() == mol.getNumBonds(), "bad bondsInPlay size");
    // RDKit❗❌:   PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
    // RDKit❗❌:                "bad atomSymbols size");
    // RDKit❗❌:   PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
    // RDKit❗❌:                "bad bondSymbols size");
    // RDKit❗❌:
    // RDKit❗❌:   if (!mol.getNumAtoms()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   bool clearRings = false;
    // RDKit❗❌:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit❗❌:     MolOps::fastFindRings(mol);
    // RDKit❗❌:     clearRings = true;
    // RDKit❗❌:   }
    // RDKit❗❌:   res.resize(mol.getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
    // RDKit❗❌:   detail::initFragmentCanonAtoms(mol, atoms, includeChirality, atomSymbols,
    // RDKit❗❌:                                  bondSymbols, atomsInPlay, bondsInPlay, true);
    // RDKit❗❌:
    // RDKit❗❌:   AtomCompareFunctor ftor(&atoms.front(), mol, &atomsInPlay, &bondsInPlay);
    // RDKit❗❌:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit❗❌:   ftor.df_useChirality = includeChirality;
    // RDKit❗❌:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit❗❌:   ftor.df_useChiralityRings = includeChirality;
    // RDKit❗❌:   ftor.df_useChiralPresence = includeChiralPresence;
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> order(mol.getNumAtoms());
    // RDKit❗❌:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit❗❌:                           includeRingStereo, &atomsInPlay, &bondsInPlay);
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     res[order[i]] = atoms[order[i]].index;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (clearRings) {
    // RDKit❗❌:     mol.getRingInfo()->reset();
    // RDKit❗❌:   }
    // RDKit❗❌: }  // end of rankFragmentAtoms()
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankFragmentAtoms (2026.03.6 complete)
    // Source publishes Fast BEFORE rank; rank errors retain it. The existing
    // prepared rank algorithm consumes these very rows without reacquiring.
    let mut replaced_by_reset = false;
    let atom_ranks;
    if canonical && !topology.atoms.is_empty() {
        #[cfg(test)]
        ring_transport_probe::record_rank();
        // A detached borrowed carrier keeps the existing prepared-rank input
        // contract. Reject its malformed dimensions through the same rank
        // owner before replacing scratch rows, preserving valence-before-ring
        // error priority without duplicating validation or running a second
        // rank on valid inputs. Live source state still publishes Fast before
        // any later rank error, as required by the upstream mutation order.
        if let KekulizeRingRows::Borrowed(supplied) = rows
            && (supplied.atom_row_count() != topology.atoms.len()
                || supplied.bond_row_count() != topology.bonds.len())
        {
            rank_fragment_atoms_with_prepared_state(
                topology,
                valence,
                Some(supplied),
                atoms_in_play,
                bonds_in_play,
                None,
                None,
                &CanonicalRankParams::kekulize_fragment_default(),
            )?;
        }
        let clear_rings = !rows
            .as_ring_info()
            .is_some_and(RingInfo::is_find_fast_or_better);
        replaced_by_reset = rows
            .as_ring_info()
            .is_some_and(|rings| !rings.is_find_fast_or_better());
        if clear_rings {
            rows.install(fast_find_rings_from_parts(
                topology.atoms.len(),
                &topology.bonds,
                &topology.adjacency,
            )?);
        }
        atom_ranks = rank_fragment_atoms_with_prepared_state(
            topology,
            valence,
            rows.as_ring_info(),
            atoms_in_play,
            bonds_in_play,
            None,
            None,
            &CanonicalRankParams::kekulize_fragment_default(),
        )?;
        if clear_rings {
            rows.reset();
        }
    } else {
        atom_ranks = (0..topology.atoms.len()).collect();
    }
    if acquisition_required && rows.as_ring_info().is_none() {
        #[cfg(test)]
        ring_transport_probe::record_sssr();
        rows.install(acquire_sssr_rows(topology)?);
    }
    Ok((replaced_by_reset, atom_ranks))
}

fn is_early_atom_for_kekulize(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: isEarlyAtom
    // RDKit✔️✔️: bool isEarlyAtom(int atomicNum) {
    // RDKit✔️✔️:   static const bool table[119] = {
    // RDKit✔️✔️:       false, false, false, true, true, true, false, false, false, false,
    // RDKit✔️✔️:       false, true, true, true, false, false, false, false, false, true,
    // RDKit✔️✔️:       true, true, true, false, false, false, false, false, false, false,
    // RDKit✔️✔️:       true, true, true, false, false, false, false, true, true, true,
    // RDKit✔️✔️:       true, true, false, false, false, false, false, false, true, true,
    // RDKit✔️✔️:       true, true, false, false, false, true, true, true, true, true,
    // RDKit✔️✔️:       true, true, false, false, false, false, false, false, false, false,
    // RDKit✔️✔️:       false, false, true, true, false, false, false, false, false, false,
    // RDKit✔️✔️:       true, true, true, true, false, false, false, true, true, true,
    // RDKit✔️✔️:       true, true, true, true, false, false, false, false, false, false,
    // RDKit✔️✔️:       false, false, false, false, true, true, true, true, true, true,
    // RDKit✔️✔️:       true, true, true, true, true, true, true, true, true,
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return ((unsigned int)atomicNum < 119) && table[atomicNum];
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: isEarlyAtom
    matches!(
        atomic_number,
        3..=5
            | 11..=13
            | 19..=22
            | 30..=32
            | 37..=41
            | 48..=51
            | 55..=61
            | 72..=73
            | 80..=83
            | 87..=93
            | 104..=118
    )
}

fn mark_double_bond_candidates(
    topology: &TopologyBlock,
    all_atoms: &[AtomId],
    rings: &RingInfo,
    valence: &ValenceAssignment,
) -> Result<CandidateState, KekulizeError> {
    let mut working = topology.clone();
    let state = mark_double_bond_candidates_mut(&mut working, all_atoms, rings, valence)?;
    Ok(CandidateState {
        topology: working,
        double_bond_candidates: state.double_bond_candidates,
        questions: state.questions,
        done: state.done,
    })
}

fn mark_double_bond_candidates_mut(
    topology: &mut TopologyBlock,
    all_atoms: &[AtomId],
    rings: &RingInfo,
    valence: &ValenceAssignment,
) -> Result<CandidateAttempt, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: markDbondCands (2026.03.6 complete)
    // RDKit❗❌: void markDbondCands(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗❌:                     boost::dynamic_bitset<> &dBndCands, INT_VECT &questions,
    // RDKit❗❌:                     INT_VECT &done) {
    // RDKit❗❌:   // ok this function does more than mark atoms that are candidates for
    // RDKit❗❌:   // double bonds during kekulization
    // RDKit❗❌:   // - check that a non-aromatic atom does not have any aromatic bonds
    // RDKit❗❌:   // - marks all aromatic bonds to single bonds
    // RDKit❗❌:   // - marks atoms that can take a double bond
    // RDKit❗❌:
    // RDKit❗❌:   bool hasAromaticOrDummyAtom =
    // RDKit❗❌:       std::any_of(allAtms.begin(), allAtms.end(), [&mol](int allAtm) {
    // RDKit❗❌:         return (!mol.getAtomWithIdx(allAtm)->getAtomicNum() ||
    // RDKit❗❌:                 isAromaticAtom(*mol.getAtomWithIdx(allAtm)));
    // RDKit❗❌:       });
    // RDKit❗❌:   // if there's not at least one atom in the ring that's
    // RDKit❗❌:   // marked as being aromatic or a dummy,
    // RDKit❗❌:   // there's no point in continuing:
    // RDKit❗❌:   if (!hasAromaticOrDummyAtom) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   // mark rings which are not candidates for double bonds
    // RDKit❗❌:   // i.e. that have at least one atom which is in a single ring
    // RDKit❗❌:   // and is not aromatic
    // RDKit❗❌:   boost::dynamic_bitset<> isRingNotCand(mol.getRingInfo()->numRings());
    // RDKit❗❌:   unsigned int ri = 0;
    // RDKit❗❌:   for (const auto &aring : mol.getRingInfo()->atomRings()) {
    // RDKit❗❌:     isRingNotCand.set(ri);
    // RDKit❗❌:     for (auto ai : aring) {
    // RDKit❗❌:       const auto at = mol.getAtomWithIdx(ai);
    // RDKit❗❌:       if (isAromaticAtom(*at) && mol.getRingInfo()->numAtomRings(ai) == 1) {
    // RDKit❗❌:         isRingNotCand.reset(ri);
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     ++ri;
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<Bond *> makeSingle;
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
    // RDKit❗❌:   for (int allAtm : allAtms) {
    // RDKit❗❌:     inAllAtms.set(allAtm);
    // RDKit❗❌:     Atom *at = mol.getAtomWithIdx(allAtm);
    // RDKit❗❌:
    // RDKit❗❌:     if (at->getAtomicNum() && !isAromaticAtom(*at)) {
    // RDKit❗❌:       done.push_back(allAtm);
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // count the number of neighbors connected with single,
    // RDKit❗❌:     // double, or aromatic bonds. Along the way, mark
    // RDKit❗❌:     // bonds that we will later mark as being single:
    // RDKit❗❌:     int sbo = 0;
    // RDKit❗❌:     unsigned nToIgnore = 0;
    // RDKit❗❌:     unsigned int nonArNonDummyNbr = 0;
    // RDKit❗❌:     for (const auto bond : mol.atomBonds(at)) {
    // RDKit❗❌:       auto otherAt = bond->getOtherAtom(at);
    // RDKit❗❌:       if (otherAt->getAtomicNum() && !otherAt->getIsAromatic() &&
    // RDKit❗❌:           inAllAtms.test(otherAt->getIdx())) {
    // RDKit❗❌:         ++nonArNonDummyNbr;
    // RDKit❗❌:       }
    // RDKit❗❌:       if (bond->getIsAromatic() && (bond->getBondType() == Bond::SINGLE ||
    // RDKit❗❌:                                     bond->getBondType() == Bond::DOUBLE ||
    // RDKit❗❌:                                     bond->getBondType() == Bond::AROMATIC)) {
    // RDKit❗❌:         ++sbo;
    // RDKit❗❌:         // mark this bond to be marked single later
    // RDKit❗❌:         // we don't want to do right now because it can screw-up the
    // RDKit❗❌:         // valence calculation to determine the number of hydrogens below
    // RDKit❗❌:         makeSingle.push_back(bond);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         int bondContrib = std::lround(bond->getValenceContrib(at));
    // RDKit❗❌:         sbo += bondContrib;
    // RDKit❗❌:         if (!bondContrib) {
    // RDKit❗❌:           ++nToIgnore;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto numAtomRings = mol.getRingInfo()->numAtomRings(at->getIdx());
    // RDKit❗❌:     const auto &riVect = mol.getRingInfo()->atomMembers(at->getIdx());
    // RDKit❗❌:     size_t numNonCandRings = std::count_if(
    // RDKit❗❌:         riVect.begin(), riVect.end(),
    // RDKit❗❌:         [&isRingNotCand](int ri) { return isRingNotCand.test(ri); });
    // RDKit❗❌:     if (!at->getAtomicNum() && nonArNonDummyNbr < numAtomRings &&
    // RDKit❗❌:         numNonCandRings < numAtomRings) {
    // RDKit❗❌:       // dummies always start as candidates to have a double bond:
    // RDKit❗❌:       dBndCands[allAtm] = 1;
    // RDKit❗❌:       // but they don't have to have one, so mark them as questionable:
    // RDKit❗❌:       questions.push_back(allAtm);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // for non dummies, it's a bit more work to figure out if they
    // RDKit❗❌:       // can take a double bond:
    // RDKit❗❌:
    // RDKit❗❌:       sbo += at->getTotalNumHs();
    // RDKit❗❌:       auto dv =
    // RDKit❗❌:           PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
    // RDKit❗❌:       auto chrg = at->getFormalCharge();
    // RDKit❗❌:       if (isEarlyAtom(at->getAtomicNum())) {
    // RDKit❗❌:         chrg = -chrg;  // fix for GitHub #65
    // RDKit❗❌:       }
    // RDKit❗❌:       // special case for carbon - see GitHub #539
    // RDKit❗❌:       if (at->getAtomicNum() == 6 && chrg > 0) {
    // RDKit❗❌:         chrg = -chrg;
    // RDKit❗❌:       }
    // RDKit❗❌:       dv += chrg;
    // RDKit❗❌:       int tbo = at->getTotalValence();
    // RDKit❗❌:       int nRadicals = at->getNumRadicalElectrons();
    // RDKit❗❌:       int totalDegree = at->getDegree() +
    // RDKit❗❌:                         at->getValence(Atom::ValenceType::IMPLICIT) - nToIgnore;
    // RDKit❗❌:
    // RDKit❗❌:       const auto &valList =
    // RDKit❗❌:           PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
    // RDKit❗❌:       unsigned int vi = 1;
    // RDKit❗❌:
    // RDKit❗❌:       while (tbo > dv && vi < valList.size() && valList[vi] > 0) {
    // RDKit❗❌:         dv = valList[vi] + chrg;
    // RDKit❗❌:         ++vi;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // Kekulize aromatic N-oxides, such as O=n1ccccc1
    // RDKit❗❌:       // These only reach here if SANITIZE_CLEANUP is disabled.
    // RDKit❗❌:       if (tbo == 5 && sbo == 4 && dv == 3 && totalDegree == 3 &&
    // RDKit❗❌:           nRadicals == 0 && chrg == 0 && at->getTotalNumHs() == 0) {
    // RDKit❗❌:         switch (at->getAtomicNum()) {
    // RDKit❗❌:           case 7:   // N
    // RDKit❗❌:           case 15:  // P
    // RDKit❗❌:           case 33:  // As
    // RDKit❗❌:             dv = 5;
    // RDKit❗❌:             break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // std::cerr << "  kek: " << at->getIdx() << " tbo:" << tbo << " sbo:" <<
    // RDKit❗❌:       // sbo
    // RDKit❗❌:       //           << "  dv : " << dv << " totalDegree : " << totalDegree
    // RDKit❗❌:       //           << " nRadicals: " << nRadicals << std::endl;
    // RDKit❗❌:       if (totalDegree + nRadicals >= dv) {
    // RDKit❗❌:         // if our degree + nRadicals exceeds the default valence,
    // RDKit❗❌:         // there's no way we can take a double bond, just continue.
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       // we're a candidate if our total current bond order + nRadicals + 1
    // RDKit❗❌:       // matches the valence state
    // RDKit❗❌:       // (including nRadicals here was SF.net issue 3349243)
    // RDKit❗❌:       if (dv == (sbo + 1 + nRadicals)) {
    // RDKit❗❌:         dBndCands[allAtm] = 1;
    // RDKit❗❌:       } else if (!nRadicals && at->getNoImplicit() && dv == (sbo + 2)) {
    // RDKit❗❌:         // special case: there is currently no radical on the atom, but if
    // RDKit❗❌:         // if we allow one then this is a candidate:
    // RDKit❗❌:         dBndCands[allAtm] = 1;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }  // loop over all atoms in the fused system
    // RDKit❗❌:
    // RDKit❗❌:   // now turn all the aromatic bond in this fused system to single
    // RDKit❗❌:   for (auto &bi : makeSingle) {
    // RDKit❗❌:     bi->setBondType(Bond::SINGLE);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: markDbondCands (2026.03.6 complete)

    let mut seen = vec![false; topology.atoms.len()];
    for &atom in all_atoms {
        if atom.index() >= topology.atoms.len() {
            return Err(KekulizeError::CandidateAtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
        if std::mem::replace(&mut seen[atom.index()], true) {
            return Err(KekulizeError::DuplicateCandidateAtom { atom });
        }
    }

    let mut double_bond_candidates = vec![false; topology.atoms.len()];
    let mut questions = Vec::new();
    let mut done = Vec::new();
    let has_aromatic_or_dummy_atom = all_atoms.iter().any(|&atom| {
        topology.atoms[atom.index()].atomic_number() == 0
            || atom_is_aromatic_for_kekulize(topology, atom)
    });
    if !has_aromatic_or_dummy_atom {
        return Ok(CandidateAttempt {
            double_bond_candidates,
            questions,
            done,
        });
    }

    let is_ring_not_candidate = rings
        .atom_rings()
        .iter()
        .map(|ring| {
            !ring.iter().any(|&atom| {
                atom_is_aromatic_for_kekulize(topology, atom) && rings.num_atom_rings(atom) == 1
            })
        })
        .collect::<Vec<_>>();
    let mut make_single = vec![false; topology.bonds.len()];
    let mut in_all_atoms = vec![false; topology.atoms.len()];

    for &atom_id in all_atoms {
        in_all_atoms[atom_id.index()] = true;
        let atom = &topology.atoms[atom_id.index()];
        if atom.atomic_number() != 0 && !atom_is_aromatic_for_kekulize(topology, atom_id) {
            done.push(atom_id);
            continue;
        }

        let mut single_bond_order = 0i32;
        let mut neighbors_to_ignore = 0usize;
        let mut non_aromatic_non_dummy_neighbors = 0usize;
        for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
            let bond = &topology.bonds[neighbor.bond.index()];
            let other = &topology.atoms[neighbor.atom_index];
            if other.atomic_number() != 0
                && !other.is_aromatic()
                && in_all_atoms[neighbor.atom_index]
            {
                non_aromatic_non_dummy_neighbors += 1;
            }
            if bond.is_aromatic()
                && matches!(
                    bond.order(),
                    BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                )
            {
                single_bond_order =
                    single_bond_order
                        .checked_add(1)
                        .ok_or(KekulizeError::IntegerOverflow {
                            atom: atom_id,
                            field: "candidate bond-order sum",
                        })?;
                make_single[bond.id().index()] = true;
            } else {
                let contribution = crate::bond_valence_contrib(bond, atom_id)?.round() as i32;
                single_bond_order = single_bond_order.checked_add(contribution).ok_or(
                    KekulizeError::IntegerOverflow {
                        atom: atom_id,
                        field: "candidate bond-order sum",
                    },
                )?;
                if contribution == 0 {
                    neighbors_to_ignore += 1;
                }
            }
        }

        let atom_ring_count = rings.num_atom_rings(atom_id);
        let non_candidate_ring_count = rings
            .atom_members(atom_id)
            .iter()
            .filter(|&&ring| is_ring_not_candidate[ring])
            .count();
        if atom.atomic_number() == 0
            && non_aromatic_non_dummy_neighbors < atom_ring_count
            && non_candidate_ring_count < atom_ring_count
        {
            double_bond_candidates[atom_id.index()] = true;
            questions.push(atom_id);
        } else {
            let total_hydrogens = i32::from(atom.explicit_hydrogens())
                .checked_add(valence.implicit_hydrogens[atom_id.index()])
                .ok_or(KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "total hydrogen count",
                })?;
            single_bond_order = single_bond_order.checked_add(total_hydrogens).ok_or(
                KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "candidate bond-order and hydrogen sum",
                },
            )?;
            let mut charge = i32::from(atom.formal_charge());
            if is_early_atom_for_kekulize(atom.atomic_number()) {
                charge = -charge;
            }
            if atom.atomic_number() == 6 && charge > 0 {
                charge = -charge;
            }
            let mut default_valence = crate::rdkit_default_valence(atom.atomic_number())?
                .checked_add(charge)
                .ok_or(KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "charged default valence",
                })?;
            let total_bond_order = checked_total_valence(valence, atom)?;
            let radical_electrons = i32::from(atom.radical_electrons());
            let degree = i32::try_from(topology.adjacency.neighbors_of(atom_id.index()).len())
                .map_err(|_| KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "atom degree",
                })?;
            let ignored =
                i32::try_from(neighbors_to_ignore).map_err(|_| KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "ignored-neighbor count",
                })?;
            let total_degree = degree
                .checked_add(valence.implicit_hydrogens[atom_id.index()])
                .and_then(|value| value.checked_sub(ignored))
                .ok_or(KekulizeError::IntegerOverflow {
                    atom: atom_id,
                    field: "total degree",
                })?;
            let valence_list = crate::required_valence_list(atom.atomic_number())?;
            let mut valence_index = 1usize;
            while total_bond_order > default_valence
                && valence_index < valence_list.len()
                && valence_list[valence_index] > 0
            {
                default_valence = valence_list[valence_index].checked_add(charge).ok_or(
                    KekulizeError::IntegerOverflow {
                        atom: atom_id,
                        field: "charged alternate valence",
                    },
                )?;
                valence_index += 1;
            }
            if total_bond_order == 5
                && single_bond_order == 4
                && default_valence == 3
                && total_degree == 3
                && radical_electrons == 0
                && charge == 0
                && total_hydrogens == 0
                && matches!(atom.atomic_number(), 7 | 15 | 33)
            {
                default_valence = 5;
            }
            if total_degree + radical_electrons >= default_valence {
                continue;
            }
            if default_valence == single_bond_order + 1 + radical_electrons
                || (radical_electrons == 0
                    && atom.no_implicit()
                    && default_valence == single_bond_order + 2)
            {
                double_bond_candidates[atom_id.index()] = true;
            }
        }
    }
    for (bond_idx, should_make_single) in make_single.into_iter().enumerate() {
        if should_make_single {
            topology.bonds[bond_idx].set_order(BondOrder::Single);
        }
    }
    Ok(CandidateAttempt {
        double_bond_candidates,
        questions,
        done,
    })
}

fn bond_between_atoms(
    topology: &TopologyBlock,
    begin: AtomId,
    end: AtomId,
) -> Result<BondId, KekulizeError> {
    topology
        .adjacency
        .neighbors_of(begin.index())
        .iter()
        .find(|neighbor| neighbor.atom_index == end.index())
        .map(|neighbor| neighbor.bond)
        .ok_or(KekulizeError::MissingCandidateBond { begin, end })
}

fn backtrack_kekulize(
    topology: &mut TopologyBlock,
    last_option: AtomId,
    done: &mut Vec<AtomId>,
    atom_queue: &mut VecDeque<AtomId>,
    double_bond_candidates: &mut [bool],
    double_bonds_added: &mut [bool],
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: backTrack
    // RDKit✔️✔️: void backTrack(RWMol &mol, INT_INT_DEQ_MAP &, int lastOpt, INT_VECT &done,
    // RDKit✔️✔️:                INT_DEQUE &aqueue, boost::dynamic_bitset<> &dBndCands,
    // RDKit✔️✔️:                boost::dynamic_bitset<> &dBndAdds) {
    // RDKit✔️✔️:   auto ei = std::find(done.begin(), done.end(), lastOpt);
    // RDKit✔️✔️:   INT_VECT tdone;
    // RDKit✔️✔️:   tdone.insert(tdone.end(), done.begin(), ei);
    let anchor = done
        .iter()
        .position(|atom| *atom == last_option)
        .ok_or(KekulizeError::MissingBacktrackAnchor { atom: last_option })?;
    let retained_done = done[..anchor].to_vec();

    // RDKit✔️✔️:   INT_VECT_CRI eri = std::find(done.rbegin(), done.rend(), lastOpt);
    // RDKit✔️✔️:   ++eri;
    // RDKit✔️✔️:   for (INT_VECT_CRI ri = done.rbegin(); ri != eri; ++ri) {
    // RDKit✔️✔️:     aqueue.push_front(*ri);
    // RDKit✔️✔️:   }
    for &atom in done[anchor..].iter().rev() {
        atom_queue.push_front(atom);
    }

    // RDKit✔️✔️:   Bond *bnd;
    // RDKit✔️✔️:   unsigned int nbnds = mol.getNumBonds();
    // RDKit✔️✔️:   for (unsigned int bi = 0; bi < nbnds; ++bi) {
    // RDKit✔️✔️:     if (dBndAdds[bi]) {
    // RDKit✔️✔️:       bnd = mol.getBondWithIdx(bi);
    // RDKit✔️✔️:       int aid1 = bnd->getBeginAtomIdx();
    // RDKit✔️✔️:       int aid2 = bnd->getEndAtomIdx();
    // RDKit✔️✔️:       if ((std::find(tdone.begin(), tdone.end(), aid1) == tdone.end()) &&
    // RDKit✔️✔️:           (std::find(tdone.begin(), tdone.end(), aid2) == tdone.end())) {
    // RDKit✔️✔️:         dBndAdds[bi] = 0;
    // RDKit✔️✔️:         bnd->setBondType(Bond::SINGLE);
    // RDKit✔️✔️:         dBndCands[aid1] = 1;
    // RDKit✔️✔️:         dBndCands[aid2] = 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond_idx in 0..topology.bonds.len() {
        if !double_bonds_added[bond_idx] {
            continue;
        }
        let begin = topology.bonds[bond_idx].begin();
        let end = topology.bonds[bond_idx].end();
        if !retained_done.contains(&begin) && !retained_done.contains(&end) {
            double_bonds_added[bond_idx] = false;
            topology.bonds[bond_idx].set_order(BondOrder::Single);
            double_bond_candidates[begin.index()] = true;
            double_bond_candidates[end.index()] = true;
        }
    }
    // RDKit✔️✔️:   done = tdone;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: backTrack
    *done = retained_done;
    Ok(())
}

fn kekulize_matching_worker(
    mut topology: TopologyBlock,
    all_atoms: &[AtomId],
    initial_candidates: &[bool],
    initial_bonds_added: &[bool],
    initial_done: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<MatchingState, KekulizeError> {
    let state = kekulize_matching_worker_mut(
        &mut topology,
        all_atoms,
        initial_candidates,
        initial_bonds_added,
        initial_done,
        atom_ranks,
        max_backtracks,
    )?;
    Ok(MatchingState {
        topology,
        succeeded: state.succeeded,
        double_bond_candidates: state.double_bond_candidates,
        double_bonds_added: state.double_bonds_added,
        done: state.done,
        problem_atoms: state.problem_atoms,
        backtracks: state.backtracks,
    })
}

fn kekulize_matching_worker_mut(
    topology: &mut TopologyBlock,
    all_atoms: &[AtomId],
    initial_candidates: &[bool],
    initial_bonds_added: &[bool],
    initial_done: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<MatchingAttempt, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeWorker (2026.03.6 complete)
    // RDKit❗❌: bool kekulizeWorker(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗❌:                     boost::dynamic_bitset<> dBndCands,
    // RDKit❗❌:                     boost::dynamic_bitset<> dBndAdds, INT_VECT done,
    // RDKit❗❌:                     const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit❗❌:   INT_DEQUE astack;
    // RDKit❗❌:   INT_INT_DEQ_MAP options;
    // RDKit❗❌:   int lastOpt = -1;
    // RDKit❗❌:   boost::dynamic_bitset<> localBondsAdded(mol.getNumBonds());
    // RDKit❗❌:   boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
    // RDKit❗❌:   for (int allAtm : allAtms) {
    // RDKit❗❌:     inAllAtms.set(allAtm);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto lessByRank = [&atomRanks](int a, int b) {
    // RDKit❗❌:     const auto ra = atomRanks.at(static_cast<unsigned int>(a));
    // RDKit❗❌:     const auto rb = atomRanks.at(static_cast<unsigned int>(b));
    // RDKit❗❌:     return (ra < rb) || (ra == rb && a < b);
    // RDKit❗❌:   };
    // RDKit❗❌:
    // RDKit❗❌:   // Prefer starting traversal at atoms which are the *end* of wedged/dashed
    // RDKit❗❌:   // bonds. Wedged bonds encode stereo and must remain single bonds; by starting
    // RDKit❗❌:   // the kekulization walk at wedge-end atoms we assign their double bond to a
    // RDKit❗❌:   // *different* neighbor first, giving the algorithm more freedom to keep the
    // RDKit❗❌:   // wedged bond single.
    // RDKit❗❌:   boost::dynamic_bitset<> wedgeEndAtoms(mol.getNumAtoms());
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     if (bond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit❗❌:         bond->getBondDir() == Bond::BondDir::BEGINDASH) {
    // RDKit❗❌:       const auto endIdx = bond->getEndAtomIdx();
    // RDKit❗❌:       if (inAllAtms.test(endIdx)) {
    // RDKit❗❌:         wedgeEndAtoms.set(endIdx);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // Pre-sort allAtms: wedge-end atoms first, then by canonical rank.
    // RDKit❗❌:   // This way the first not-yet-done atom is always the best starting point.
    // RDKit❗❌:   INT_VECT sortedAtms(allAtms);
    // RDKit❗❌:   std::sort(sortedAtms.begin(), sortedAtms.end(),
    // RDKit❗❌:             [&wedgeEndAtoms, &lessByRank](int a, int b) {
    // RDKit❗❌:               const bool wa = wedgeEndAtoms.test(a);
    // RDKit❗❌:               const bool wb = wedgeEndAtoms.test(b);
    // RDKit❗❌:               if (wa != wb) {
    // RDKit❗❌:                 return wa;  // wedge-end atoms come first
    // RDKit❗❌:               }
    // RDKit❗❌:               return lessByRank(a, b);
    // RDKit❗❌:             });
    // RDKit❗❌:
    // RDKit❗❌:   // ok the algorithm goes something like this
    // RDKit❗❌:   // - start with an atom that has been marked aromatic before
    // RDKit❗❌:   // - check if it can have a double bond
    // RDKit❗❌:   // - add its neighbors to the stack
    // RDKit❗❌:   // - check if one of its neighbors can also have a double bond
    // RDKit❗❌:   // - if yes add a double bond.
    // RDKit❗❌:   // - if multiple neighbors can have double bonds - add them to a
    // RDKit❗❌:   //   options stack we may have to retrace out path if we chose the
    // RDKit❗❌:   //   wrong neighbor to add the double bond
    // RDKit❗❌:   // - if double bond added update the candidates for double bond
    // RDKit❗❌:   // - move to the next atom on the stack and repeat the process
    // RDKit❗❌:   // - if an atom that can have multiple a double bond has no
    // RDKit❗❌:   //   neighbors that can take double bond - we made a mistake
    // RDKit❗❌:   //   earlier by picking a wrong candidate for double bond
    // RDKit❗❌:   // - in this case back track to where we made the mistake
    // RDKit❗❌:
    // RDKit❗❌:   int curr = -1;
    // RDKit❗❌:   INT_DEQUE btmoves;
    // RDKit❗❌:   unsigned int numBT = 0;  // number of back tracks so far
    // RDKit❗❌:   while ((done.size() < sortedAtms.size()) || !astack.empty()) {
    // RDKit❗❌:     // pick a curr atom to work with
    // RDKit❗❌:     if (astack.size() > 0) {
    // RDKit❗❌:       curr = astack.front();
    // RDKit❗❌:       astack.pop_front();
    // RDKit❗❌:     } else {
    // RDKit❗❌:       curr = -1;
    // RDKit❗❌:       for (int allAtm : sortedAtms) {
    // RDKit❗❌:         if (std::find(done.begin(), done.end(), allAtm) == done.end()) {
    // RDKit❗❌:           curr = allAtm;
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     CHECK_INVARIANT(curr >= 0, "starting point not found");
    // RDKit❗❌:     done.push_back(curr);
    // RDKit❗❌:
    // RDKit❗❌:     // loop over the neighbors if we can add double bonds or
    // RDKit❗❌:     // simply push them onto the stack
    // RDKit❗❌:     INT_DEQUE opts;
    // RDKit❗❌:     bool cCand = false;
    // RDKit❗❌:     if (dBndCands[curr]) {
    // RDKit❗❌:       cCand = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     int ncnd;
    // RDKit❗❌:     // if we are here because of backtracking
    // RDKit❗❌:     if (options.find(curr) != options.end()) {
    // RDKit❗❌:       opts = options[curr];
    // RDKit❗❌:       CHECK_INVARIANT(opts.size() > 0, "");
    // RDKit❗❌:     } else {
    // RDKit❗❌:       INT_DEQUE lstack;
    // RDKit❗❌:       std::vector<int> optsV;
    // RDKit❗❌:       std::vector<int> wedgedOptsV;
    // RDKit❗❌:       std::vector<int> nbrs;
    // RDKit❗❌:       for (auto nbrAtom : mol.atomNeighbors(mol.getAtomWithIdx(curr))) {
    // RDKit❗❌:         const auto nbrIdx = static_cast<int>(nbrAtom->getIdx());
    // RDKit❗❌:         // ignore if the neighbor is not part of the fused system
    // RDKit❗❌:         if (!inAllAtms.test(nbrIdx)) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         // ignore if the neighbor has already been dealt with before
    // RDKit❗❌:         if (std::find(done.begin(), done.end(), nbrIdx) != done.end()) {
    // RDKit❗❌:           continue;
    // RDKit❗❌:         }
    // RDKit❗❌:         nbrs.push_back(nbrIdx);
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       std::sort(nbrs.begin(), nbrs.end(), lessByRank);
    // RDKit❗❌:
    // RDKit❗❌:       for (int nbrIdx : nbrs) {
    // RDKit❗❌:         auto nbrBond = mol.getBondBetweenAtoms(curr, nbrIdx);
    // RDKit❗❌:
    // RDKit❗❌:         // if the neighbor is not on the stack add it
    // RDKit❗❌:         if (std::find(astack.begin(), astack.end(), nbrIdx) == astack.end()) {
    // RDKit❗❌:           lstack.push_back(nbrIdx);
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // check if the neighbor is also a candidate for a double bond
    // RDKit❗❌:         // the refinement that we'll make to the candidate check we've already
    // RDKit❗❌:         // done is to make sure that the bond is either flagged as aromatic
    // RDKit❗❌:         // or involves a dummy atom. This was Issue 3525076.
    // RDKit❗❌:         // This fix is not really 100% of the way there: a situation like
    // RDKit❗❌:         // that for Issue 3525076 but involving a dummy atom in the cage
    // RDKit❗❌:         // could lead to the same failure. The full fix would require
    // RDKit❗❌:         // a fairly detailed analysis of all bonds in the molecule to determine
    // RDKit❗❌:         // which of them is eligible to be converted.
    // RDKit❗❌:         if (cCand && dBndCands[nbrIdx] &&
    // RDKit❗❌:             (nbrBond->getIsAromatic() ||
    // RDKit❗❌:              mol.getAtomWithIdx(curr)->getAtomicNum() == 0 ||
    // RDKit❗❌:              mol.getAtomWithIdx(nbrIdx)->getAtomicNum() == 0)) {
    // RDKit❗❌:           // in order to try and avoid making wedged bonds double, we will add
    // RDKit❗❌:           // this neighbor at the back of the options after this loop if the
    // RDKit❗❌:           // bond is wedged. otherwise we append it to the options directly
    // RDKit❗❌:           if (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit❗❌:               nbrBond->getBondDir() == Bond::BondDir::BEGINDASH) {
    // RDKit❗❌:             wedgedOptsV.push_back(nbrIdx);
    // RDKit❗❌:           } else {
    // RDKit❗❌:             optsV.push_back(nbrIdx);
    // RDKit❗❌:           }
    // RDKit❗❌:         }  // end of curr atoms can have a double bond
    // RDKit❗❌:       }  // end of looping over neighbors
    // RDKit❗❌:
    // RDKit❗❌:       // Non-wedged options first, then wedged — both already in rank order
    // RDKit❗❌:       // because nbrs was pre-sorted by lessByRank above.
    // RDKit❗❌:       for (int v : optsV) {
    // RDKit❗❌:         opts.push_back(v);
    // RDKit❗❌:       }
    // RDKit❗❌:       for (int v : wedgedOptsV) {
    // RDKit❗❌:         opts.push_back(v);
    // RDKit❗❌:       }
    // RDKit❗❌:       astack.insert(astack.end(), lstack.begin(), lstack.end());
    // RDKit❗❌:     }
    // RDKit❗❌:     // now add a double bond from current to one of the neighbors if we can
    // RDKit❗❌:     if (cCand) {
    // RDKit❗❌:       if (!opts.empty()) {
    // RDKit❗❌:         ncnd = opts.front();
    // RDKit❗❌:         opts.pop_front();
    // RDKit❗❌:         auto bnd = mol.getBondBetweenAtoms(curr, ncnd);
    // RDKit❗❌:         bnd->setBondType(Bond::DOUBLE);
    // RDKit❗❌:         if (bnd->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗❌:           bnd->setBondDir(Bond::BondDir::NONE);
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // remove current and the neighbor from the dBndCands list
    // RDKit❗❌:         dBndCands[curr] = 0;
    // RDKit❗❌:         dBndCands[ncnd] = 0;
    // RDKit❗❌:
    // RDKit❗❌:         // add them to the list of bonds to which have been made double
    // RDKit❗❌:         dBndAdds[bnd->getIdx()] = 1;
    // RDKit❗❌:         localBondsAdded[bnd->getIdx()] = 1;
    // RDKit❗❌:
    // RDKit❗❌:         // if this is an atom we previously visted and picked we
    // RDKit❗❌:         // simply tried a different option now, overwrite the options
    // RDKit❗❌:         // stored for this atoms
    // RDKit❗❌:         if (options.find(curr) != options.end()) {
    // RDKit❗❌:           if (opts.size() == 0) {
    // RDKit❗❌:             options.erase(curr);
    // RDKit❗❌:             btmoves.pop_back();
    // RDKit❗❌:             if (btmoves.size() > 0) {
    // RDKit❗❌:               lastOpt = btmoves.back();
    // RDKit❗❌:             } else {
    // RDKit❗❌:               lastOpt = -1;
    // RDKit❗❌:             }
    // RDKit❗❌:           } else {
    // RDKit❗❌:             options[curr] = opts;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // this is new atoms we are trying and have other
    // RDKit❗❌:           // neighbors as options to add double bond store this to
    // RDKit❗❌:           // the options stack, we may have made a mistake in
    // RDKit❗❌:           // which one we chose and have to return here
    // RDKit❗❌:           if (opts.size() > 0) {
    // RDKit❗❌:             lastOpt = curr;
    // RDKit❗❌:             btmoves.push_back(lastOpt);
    // RDKit❗❌:             options[curr] = opts;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:       }  // end of adding a double bond
    // RDKit❗❌:       else if (mol.getAtomWithIdx(curr)->getAtomicNum()) {
    // RDKit❗❌:         // we have a non-dummy atom that should be getting a double
    // RDKit❗❌:         // bond but none of the neighbors can take one. Most likely
    // RDKit❗❌:         // because of a wrong choice earlier so back track
    // RDKit❗❌:         if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
    // RDKit❗❌:           // std::cerr << "PRE BACKTRACK" << std::endl;
    // RDKit❗❌:           // mol.debugMol(std::cerr);
    // RDKit❗❌:           backTrack(mol, options, lastOpt, done, astack, dBndCands, dBndAdds);
    // RDKit❗❌:           // std::cerr << "POST BACKTRACK" << std::endl;
    // RDKit❗❌:           // mol.debugMol(std::cerr);
    // RDKit❗❌:           ++numBT;
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // undo any remaining changes we made while here
    // RDKit❗❌:           // this was github #962
    // RDKit❗❌:           for (unsigned int bidx = 0; bidx < mol.getNumBonds(); ++bidx) {
    // RDKit❗❌:             if (localBondsAdded[bidx]) {
    // RDKit❗❌:               mol.getBondWithIdx(bidx)->setBondType(Bond::SINGLE);
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:       }  // end of else try to backtrack
    // RDKit❗❌:     }  // end of curr atom atom being a cand for double bond
    // RDKit❗❌:   }  // end of while we are not done with all atoms
    // RDKit❗❌:   return true;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeWorker (2026.03.6 complete)

    topology.validate()?;
    for (field, expected, actual) in [
        (
            "double-bond candidates",
            topology.atoms.len(),
            initial_candidates.len(),
        ),
        (
            "double bonds added",
            topology.bonds.len(),
            initial_bonds_added.len(),
        ),
        ("atom ranks", topology.atoms.len(), atom_ranks.len()),
    ] {
        if actual != expected {
            return Err(KekulizeError::MatchingStateLength {
                field,
                expected,
                actual,
            });
        }
    }
    let mut in_all_atoms = vec![false; topology.atoms.len()];
    for &atom in all_atoms {
        if atom.index() >= topology.atoms.len() {
            return Err(KekulizeError::CandidateAtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
        if std::mem::replace(&mut in_all_atoms[atom.index()], true) {
            return Err(KekulizeError::DuplicateCandidateAtom { atom });
        }
    }
    let mut seen_done = vec![false; topology.atoms.len()];
    for &atom in initial_done {
        if atom.index() >= topology.atoms.len() {
            return Err(KekulizeError::CandidateAtomOutOfRange {
                atom,
                atom_count: topology.atoms.len(),
            });
        }
        if std::mem::replace(&mut seen_done[atom.index()], true) {
            return Err(KekulizeError::DuplicateDoneAtom { atom });
        }
    }

    let mut atom_stack = VecDeque::new();
    let mut options = BTreeMap::<AtomId, VecDeque<AtomId>>::new();
    let mut last_option = None;
    let mut local_bonds_added = vec![false; topology.bonds.len()];
    let mut double_bond_candidates = initial_candidates.to_vec();
    let mut double_bonds_added = initial_bonds_added.to_vec();
    let mut done = initial_done.to_vec();

    let mut wedge_end_atoms = vec![false; topology.atoms.len()];
    for bond in &topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) && in_all_atoms[bond.end().index()]
        {
            wedge_end_atoms[bond.end().index()] = true;
        }
    }

    let mut sorted_atoms = all_atoms.to_vec();
    sorted_atoms.sort_by_key(|atom| {
        (
            !wedge_end_atoms[atom.index()],
            atom_ranks[atom.index()],
            atom.index(),
        )
    });

    let mut backtrack_moves = Vec::new();
    let mut backtracks = 0u32;
    while done.len() < sorted_atoms.len() || !atom_stack.is_empty() {
        let current = atom_stack
            .pop_front()
            .or_else(|| {
                sorted_atoms
                    .iter()
                    .copied()
                    .find(|atom| !done.contains(atom))
            })
            .ok_or(KekulizeError::MissingBacktrackAnchor {
                atom: AtomId::new(topology.atoms.len()),
            })?;
        done.push(current);

        let current_is_candidate = double_bond_candidates[current.index()];
        let mut current_options = if let Some(stored) = options.get(&current) {
            stored.clone()
        } else {
            let mut local_stack = VecDeque::new();
            let mut ordinary_options = VecDeque::new();
            let mut wedged_options = VecDeque::new();
            let mut neighbors = topology
                .adjacency
                .neighbors_of(current.index())
                .iter()
                .map(|neighbor| AtomId::new(neighbor.atom_index))
                .filter(|neighbor| in_all_atoms[neighbor.index()] && !done.contains(neighbor))
                .collect::<Vec<_>>();
            neighbors.sort_by_key(|atom| (atom_ranks[atom.index()], atom.index()));

            for neighbor in neighbors {
                let bond_id = bond_between_atoms(topology, current, neighbor)?;
                let bond = &topology.bonds[bond_id.index()];
                if !atom_stack.contains(&neighbor) {
                    local_stack.push_back(neighbor);
                }
                if current_is_candidate
                    && double_bond_candidates[neighbor.index()]
                    && (bond.is_aromatic()
                        || topology.atoms[current.index()].atomic_number() == 0
                        || topology.atoms[neighbor.index()].atomic_number() == 0)
                {
                    if matches!(
                        bond.direction(),
                        BondDirection::BeginWedge | BondDirection::BeginDash
                    ) {
                        wedged_options.push_back(neighbor);
                    } else {
                        ordinary_options.push_back(neighbor);
                    }
                }
            }
            ordinary_options.append(&mut wedged_options);
            atom_stack.append(&mut local_stack);
            ordinary_options
        };

        if current_is_candidate {
            if let Some(neighbor) = current_options.pop_front() {
                let bond_id = bond_between_atoms(topology, current, neighbor)?;
                topology.bonds[bond_id.index()].set_order(BondOrder::Double);
                if topology.bonds[bond_id.index()].direction() != BondDirection::None {
                    topology.bonds[bond_id.index()].set_direction(BondDirection::None);
                }
                double_bond_candidates[current.index()] = false;
                double_bond_candidates[neighbor.index()] = false;
                double_bonds_added[bond_id.index()] = true;
                local_bonds_added[bond_id.index()] = true;

                if options.contains_key(&current) {
                    if current_options.is_empty() {
                        options.remove(&current);
                        backtrack_moves.pop();
                        last_option = backtrack_moves.last().copied();
                    } else {
                        options.insert(current, current_options);
                    }
                } else {
                    if !current_options.is_empty() {
                        last_option = Some(current);
                        backtrack_moves.push(current);
                        options.insert(current, current_options);
                    }
                }
            } else if topology.atoms[current.index()].atomic_number() != 0 {
                if let Some(anchor) = last_option.filter(|_| backtracks < max_backtracks) {
                    backtrack_kekulize(
                        topology,
                        anchor,
                        &mut done,
                        &mut atom_stack,
                        &mut double_bond_candidates,
                        &mut double_bonds_added,
                    )?;
                    backtracks += 1;
                } else {
                    for (bond_idx, was_added) in local_bonds_added.iter().copied().enumerate() {
                        if was_added {
                            topology.bonds[bond_idx].set_order(BondOrder::Single);
                        }
                    }
                    let problem_atoms = all_atoms
                        .iter()
                        .copied()
                        .filter(|atom| double_bond_candidates[atom.index()])
                        .collect();
                    return Ok(MatchingAttempt {
                        succeeded: false,
                        double_bond_candidates,
                        double_bonds_added,
                        done,
                        problem_atoms,
                        backtracks,
                    });
                }
            }
        }
    }
    Ok(MatchingAttempt {
        succeeded: true,
        double_bond_candidates,
        double_bonds_added,
        done,
        problem_atoms: Vec::new(),
        backtracks,
    })
}

impl QuestionEnumerator {
    fn new(questions: Vec<AtomId>) -> Result<Self, KekulizeError> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::QuestionEnumerator
        // RDKit✔️✔️:   QuestionEnumerator(INT_VECT questions)
        // RDKit✔️✔️:       : d_questions(std::move(questions)), d_state(d_questions.size()) {
        // RDKit✔️✔️:     if (!d_state.empty()) {
        // RDKit✔️✔️:       // Start at one because the empty subset has already been attempted.
        // RDKit✔️✔️:       d_state.set(0);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::QuestionEnumerator
        // RDKit✔️✔️: dynamic_bitset-equivalent packed words; preserve the
        // existing private Result signature while removing its fixed-width gate.
        // Space is ceil(Q/64) words, not the one-byte-per-bool Vec<bool> layout.
        let mut state = vec![0; questions.len().div_ceil(64)];
        if !questions.is_empty() {
            state[0] = 1;
        }
        Ok(Self {
            questions,
            state,
            done: false,
        })
    }

    fn next(&mut self) -> Vec<AtomId> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::next
        // RDKit✔️✔️:   INT_VECT next() {
        // RDKit✔️✔️:     INT_VECT res;
        // RDKit✔️✔️:     if (d_done) {
        // RDKit✔️✔️:       return res;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     for (size_t i = 0; i < d_questions.size(); ++i) {
        // RDKit✔️✔️:       if (d_state.test(i)) {
        // RDKit✔️✔️:         res.push_back(d_questions[i]);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     size_t pos = 0;
        // RDKit✔️✔️:     while (pos < d_state.size() && d_state.test(pos)) {
        // RDKit✔️✔️:       d_state.reset(pos++);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (pos == d_state.size()) {
        // RDKit✔️✔️:       d_done = true;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       d_state.set(pos);
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::next
        // RDKit✔️✔️: scan Q positions and allocate only this subset, as source.
        // The packed carry performs the same ordered bit reset/set operations;
        // bits outside Q in the last word are never read, set or shifted into.
        if self.done {
            return Vec::new();
        }
        let mut selected = Vec::new();
        for (index, &question) in self.questions.iter().enumerate() {
            if self.state[index / 64] & (1u64 << (index % 64)) != 0 {
                selected.push(question);
            }
        }
        let mut pos = 0;
        while pos < self.questions.len() && self.state[pos / 64] & (1u64 << (pos % 64)) != 0 {
            self.state[pos / 64] &= !(1u64 << (pos % 64));
            pos += 1;
        }
        if pos == self.questions.len() {
            self.done = true;
        } else {
            self.state[pos / 64] |= 1u64 << (pos % 64);
        }
        selected
    }
}

fn permute_dummies_and_kekulize(
    mut topology: TopologyBlock,
    all_atoms: &[AtomId],
    initial_candidates: &[bool],
    questions: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedKekulizeState, KekulizeError> {
    let state = permute_dummies_and_kekulize_mut(
        &mut topology,
        all_atoms,
        initial_candidates,
        questions,
        atom_ranks,
        max_backtracks,
    )?;
    Ok(FusedKekulizeState {
        topology,
        succeeded: state.succeeded,
        problem_atoms: state.problem_atoms,
    })
}

fn permute_dummies_and_kekulize_mut(
    topology: &mut TopologyBlock,
    all_atoms: &[AtomId],
    initial_candidates: &[bool],
    questions: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedAttempt, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: permuteDummiesAndKekulize (2026.03.6 complete)
    // RDKit❗❌: bool permuteDummiesAndKekulize(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗❌:                                boost::dynamic_bitset<> dBndCands,
    // RDKit❗❌:                                INT_VECT &questions, const UINT_VECT &atomRanks,
    // RDKit❗❌:                                unsigned int maxBackTracks) {
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    // RDKit❗❌:   for (int allAtm : allAtms) {
    // RDKit❗❌:     atomsInPlay[allAtm] = 1;
    // RDKit❗❌:   }
    // RDKit❗❌:   bool kekulized = false;
    // RDKit❗❌:   QuestionEnumerator qEnum(questions);
    // RDKit❗❌:   while (!kekulized && questions.size()) {
    // RDKit❗❌:     boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
    // RDKit❗❌:     INT_VECT done;
    // RDKit❗❌:     // reset the state: all aromatic bonds are remarked to single:
    // RDKit❗❌:     for (const auto bond : mol.bonds()) {
    // RDKit❗❌:       if (bond->getIsAromatic() && bond->getBondType() != Bond::SINGLE &&
    // RDKit❗❌:           atomsInPlay[bond->getBeginAtomIdx()] &&
    // RDKit❗❌:           atomsInPlay[bond->getEndAtomIdx()]) {
    // RDKit❗❌:         bond->setBondType(Bond::SINGLE);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // pick a new permutation of the questionable atoms:
    // RDKit❗❌:     const auto &switchOff = qEnum.next();
    // RDKit❗❌:     if (!switchOff.size()) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto tCands = dBndCands;
    // RDKit❗❌:     for (int it : switchOff) {
    // RDKit❗❌:       tCands[it] = 0;
    // RDKit❗❌:     }
    // RDKit❗❌:     // try kekulizing again:
    // RDKit❗❌:     kekulized = kekulizeWorker(mol, allAtms, tCands, dBndAdds, done, atomRanks,
    // RDKit❗❌:                                maxBackTracks);
    // RDKit❗❌:   }
    // RDKit❗❌:   return kekulized;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: permuteDummiesAndKekulize (2026.03.6 complete)

    let mut atoms_in_play = vec![false; topology.atoms.len()];
    for &atom in all_atoms {
        atoms_in_play[atom.index()] = true;
    }
    let mut question_enumerator = QuestionEnumerator::new(questions.to_vec())?;
    while !questions.is_empty() {
        let double_bonds_added = vec![false; topology.bonds.len()];
        for bond in &mut topology.bonds {
            if bond.is_aromatic()
                && bond.order() != BondOrder::Single
                && atoms_in_play[bond.begin().index()]
                && atoms_in_play[bond.end().index()]
            {
                bond.set_order(BondOrder::Single);
            }
        }
        let switch_off = question_enumerator.next();
        if switch_off.is_empty() {
            break;
        }
        let mut trial_candidates = initial_candidates.to_vec();
        for atom in switch_off {
            trial_candidates[atom.index()] = false;
        }
        let trial = kekulize_matching_worker_mut(
            topology,
            all_atoms,
            &trial_candidates,
            &double_bonds_added,
            &[],
            atom_ranks,
            max_backtracks,
        )?;
        if trial.succeeded {
            return Ok(FusedAttempt {
                succeeded: true,
                problem_atoms: Vec::new(),
            });
        }
    }
    Ok(FusedAttempt {
        succeeded: false,
        problem_atoms: initial_candidates
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(index, candidate)| candidate.then(|| AtomId::new(index)))
            .collect(),
    })
}

fn make_ring_neighbor_map(bond_rings: &[Vec<BondId>]) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::makeRingNeighborMap
    // RDKit✔️✔️: void makeRingNeighborMap(const VECT_INT_VECT &brings,
    // RDKit✔️✔️:                          INT_INT_VECT_MAP &neighMap, unsigned int maxSize,
    // RDKit✔️✔️:                          unsigned int maxOverlapSize) {
    // RDKit✔️✔️:   auto nrings = rdcast<int>(brings.size());
    // RDKit✔️✔️:   int i, j;
    // RDKit✔️✔️:   INT_VECT ring1;
    let mut neighbors = vec![Vec::new(); bond_rings.len()];
    // RDKit✔️✔️:   for (i = 0; i < nrings; ++i) {
    for first in 0..bond_rings.len() {
        // RDKit✔️✔️:     // create an empty INT_VECT at neighMap[i] if it does not yet exist
        // RDKit✔️✔️:     neighMap[i];
        // RDKit✔️✔️:     if (maxSize && brings[i].size() > maxSize) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ring1 = brings[i];
        // The kekulization call uses both source defaults (`maxSize == 0` and
        // `maxOverlapSize == 0`), so no ring is filtered by either option.
        // RDKit✔️✔️:     for (j = i + 1; j < nrings; ++j) {
        for second in first + 1..bond_rings.len() {
            // RDKit✔️✔️:       if (maxSize && brings[j].size() > maxSize) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       INT_VECT inter;
            // RDKit✔️✔️:       Intersect(ring1, brings[j], inter);
            // RDKit✔️✔️:       if (inter.size() > 0 &&
            // RDKit✔️✔️:           (!maxOverlapSize || inter.size() <= maxOverlapSize)) {
            if bond_rings[first]
                .iter()
                .any(|bond| bond_rings[second].contains(bond))
            {
                // RDKit✔️✔️:         neighMap[i].push_back(j);
                // RDKit✔️✔️:         neighMap[j].push_back(i);
                neighbors[first].push(second);
                neighbors[second].push(first);
                // RDKit✔️✔️:       }
            }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::makeRingNeighborMap
    neighbors
}

fn pick_fused_rings(current: usize, neighbor_map: &[Vec<usize>], done: &mut [bool]) -> Vec<usize> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::pickFusedRings
    // RDKit✔️✔️: void pickFusedRings(int curr, const INT_INT_VECT_MAP &neighMap, INT_VECT &res,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> &done, int depth) {
    // RDKit✔️✔️:   auto pos = neighMap.find(curr);
    // RDKit✔️✔️:   PRECONDITION(pos != neighMap.end(), "bad argument");
    // RDKit✔️✔️:   done[curr] = 1;
    // RDKit✔️✔️:   res.push_back(curr);
    done[current] = true;
    let mut result = vec![current];
    // An explicit frame stack preserves the recursive source's exact DFS
    // order while avoiding call-stack exhaustion on a large ring graph.
    let mut frames = vec![(current, 0usize)];
    // RDKit✔️✔️:   const auto &neighs = pos->second;
    // RDKit✔️✔️:   for (int neigh : neighs) {
    // RDKit✔️✔️:     if (!done[neigh]) {
    // RDKit✔️✔️:       pickFusedRings(neigh, neighMap, res, done, depth + 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    while let Some((ring, next_neighbor)) = frames.last_mut() {
        let Some(&neighbor) = neighbor_map[*ring].get(*next_neighbor) else {
            frames.pop();
            continue;
        };
        *next_neighbor += 1;
        if !done[neighbor] {
            done[neighbor] = true;
            result.push(neighbor);
            frames.push((neighbor, 0));
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Aromaticity.cpp :: RingUtils::pickFusedRings
    result
}

fn kekulize_fused_system(
    mut topology: TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    all_rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedKekulizeState, KekulizeError> {
    let state = kekulize_fused_system_mut(
        &mut topology,
        atom_rings,
        all_rings,
        valence,
        atom_ranks,
        max_backtracks,
    )?;
    Ok(FusedKekulizeState {
        topology,
        succeeded: state.succeeded,
        problem_atoms: state.problem_atoms,
    })
}

fn kekulize_fused_system_mut(
    topology: &mut TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    all_rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedAttempt, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeFused (2026.03.6 complete)
    // RDKit❗❌: void kekulizeFused(RWMol &mol, const VECT_INT_VECT &arings,
    // RDKit❗❌:                    const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit❗❌:   // get all the atoms in the ring system
    // RDKit❗❌:   INT_VECT allAtms;
    // RDKit❗❌:   Union(arings, allAtms);
    // RDKit❗❌:   // get all the atoms that are candidates to receive a double bond
    // RDKit❗❌:   // also mark atoms in the fused system that are not aromatic to begin with
    // RDKit❗❌:   // as done. Mark all the bonds that are part of the aromatic system
    // RDKit❗❌:   // to be single bonds
    // RDKit❗❌:   INT_VECT done;
    // RDKit❗❌:   INT_VECT questions;
    // RDKit❗❌:   auto nats = mol.getNumAtoms();
    // RDKit❗❌:   auto nbnds = mol.getNumBonds();
    // RDKit❗❌:   boost::dynamic_bitset<> dBndCands(nats);
    // RDKit❗❌:   boost::dynamic_bitset<> dBndAdds(nbnds);
    // RDKit❗❌:   markDbondCands(mol, allAtms, dBndCands, questions, done);
    // RDKit❗❌:
    // RDKit❗❌:   auto kekulized = kekulizeWorker(mol, allAtms, dBndCands, dBndAdds, done,
    // RDKit❗❌:                                   atomRanks, maxBackTracks);
    // RDKit❗❌:   if (!kekulized && questions.size()) {
    // RDKit❗❌:     // we failed, but there are some dummy atoms we can try permuting.
    // RDKit❗❌:     kekulized = permuteDummiesAndKekulize(mol, allAtms, dBndCands, questions,
    // RDKit❗❌:                                           atomRanks, maxBackTracks);
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!kekulized) {
    // RDKit❗❌:     // we exhausted all option (or crossed the allowed
    // RDKit❗❌:     // number of backTracks) and we still need to backtrack
    // RDKit❗❌:     // can't kekulize this thing
    // RDKit❗❌:     std::vector<unsigned int> problemAtoms;
    // RDKit❗❌:     std::ostringstream errout;
    // RDKit❗❌:     errout << "Can't kekulize mol.";
    // RDKit❗❌:     errout << "  Unkekulized atoms:";
    // RDKit❗❌:     for (unsigned int i = 0; i < nats; ++i) {
    // RDKit❗❌:       if (dBndCands[i]) {
    // RDKit❗❌:         errout << " " << i;
    // RDKit❗❌:         problemAtoms.push_back(i);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     std::string msg = errout.str();
    // RDKit❗❌:     BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit❗❌:     throw KekulizeException(msg, problemAtoms);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeFused (2026.03.6 complete)

    let mut all_atoms = Vec::new();
    for ring in atom_rings {
        for &atom in ring {
            if !all_atoms.contains(&atom) {
                all_atoms.push(atom);
            }
        }
    }
    let candidate_state =
        mark_double_bond_candidates_mut(topology, &all_atoms, all_rings, valence)?;
    let initial_candidates = candidate_state.double_bond_candidates.clone();
    let questions = candidate_state.questions.clone();
    let initial_bonds_added = vec![false; topology.bonds.len()];
    let first_attempt = kekulize_matching_worker_mut(
        topology,
        &all_atoms,
        &initial_candidates,
        &initial_bonds_added,
        &candidate_state.done,
        atom_ranks,
        max_backtracks,
    )?;
    if first_attempt.succeeded {
        return Ok(FusedAttempt {
            succeeded: true,
            problem_atoms: Vec::new(),
        });
    }
    if !questions.is_empty() {
        let permuted = permute_dummies_and_kekulize_mut(
            topology,
            &all_atoms,
            &initial_candidates,
            &questions,
            atom_ranks,
            max_backtracks,
        )?;
        if permuted.succeeded {
            return Ok(permuted);
        }
        return Ok(permuted);
    }
    Ok(FusedAttempt {
        succeeded: false,
        problem_atoms: initial_candidates
            .iter()
            .copied()
            .enumerate()
            .filter_map(|(index, candidate)| candidate.then(|| AtomId::new(index)))
            .collect(),
    })
}

fn kekulize_fused_components(
    mut topology: TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    bond_rings: &[Vec<BondId>],
    all_rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedKekulizeState, KekulizeError> {
    let state = kekulize_fused_components_mut(
        &mut topology,
        atom_rings,
        bond_rings,
        all_rings,
        valence,
        atom_ranks,
        max_backtracks,
    )?;
    Ok(FusedKekulizeState {
        topology,
        succeeded: state.succeeded,
        problem_atoms: state.problem_atoms,
    })
}

fn kekulize_fused_components_mut(
    topology: &mut TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    bond_rings: &[Vec<BondId>],
    all_rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedAttempt, KekulizeError> {
    let neighbor_map = make_ring_neighbor_map(bond_rings);
    let mut done = vec![false; atom_rings.len()];
    for current in 0..atom_rings.len() {
        if done[current] {
            continue;
        }
        let fused = pick_fused_rings(current, &neighbor_map, &mut done);
        let fused_atom_rings = fused
            .into_iter()
            .map(|ring| atom_rings[ring].clone())
            .collect::<Vec<_>>();
        let state = kekulize_fused_system_mut(
            topology,
            &fused_atom_rings,
            all_rings,
            valence,
            atom_ranks,
            max_backtracks,
        )?;
        if !state.succeeded {
            return Ok(FusedAttempt {
                succeeded: false,
                problem_atoms: state.problem_atoms,
            });
        }
    }
    Ok(FusedAttempt {
        succeeded: true,
        problem_atoms: Vec::new(),
    })
}

fn kekulize_fragment(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
    rings: Option<&RingInfo>,
    source_valence: Option<&ValenceAssignment>,
) -> Result<KekulizeAssignment, KekulizeError> {
    let mut working = topology.clone();
    let mut valence = source_valence
        .cloned()
        .unwrap_or_else(|| ValenceAssignment {
            explicit_valence: vec![-1; topology.atoms.len()],
            implicit_hydrogens: vec![-1; topology.atoms.len()],
        });
    let mut rows = match rings {
        Some(rings) if rings.is_initialized() => KekulizeRingRows::Borrowed(rings),
        _ => KekulizeRingRows::Uninitialized,
    };
    let (refreshed_valence_atoms, reset) = kekulize_fragment_attempt(
        &mut working,
        &mut valence,
        &mut rows,
        atoms_in_play,
        bonds_in_play,
        params,
        query_state,
    )?;
    Ok(KekulizeAssignment {
        topology: working,
        final_valence: if atoms_in_play.iter().any(|selected| *selected) {
            Some(valence)
        } else {
            source_valence.cloned()
        },
        refreshed_valence_atoms,
        ring_update: rows.into_update(reset),
    })
}

fn kekulize_fragment_attempt(
    topology: &mut TopologyBlock,
    valence: &mut ValenceAssignment,
    rows: &mut KekulizeRingRows<'_>,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<(Vec<AtomId>, bool), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment (2026.03.6 complete)
    // RDKit❗❌: void KekulizeFragment(RWMol &mol, const boost::dynamic_bitset<> &atomsToUse,
    // RDKit❗❌:                       boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds,
    // RDKit❗❌:                       bool canonical, unsigned int maxBackTracks) {
    // RDKit❗❌:   PRECONDITION(atomsToUse.size() == mol.getNumAtoms(),
    // RDKit❗❌:                "atomsToUse is wrong size");
    // RDKit❗❌:   PRECONDITION(bondsToUse.size() == mol.getNumBonds(),
    // RDKit❗❌:                "bondsToUse is wrong size");
    // RDKit❗❌:   // if there are no atoms to use we can directly return
    // RDKit❗❌:   if (atomsToUse.none()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // there's no point doing kekulization if there are no aromatic bonds
    // RDKit❗❌:   // without queries:
    // RDKit❗❌:   bool foundAromatic = false;
    // RDKit❗❌:   for (const auto bond : mol.bonds()) {
    // RDKit❗❌:     if (bondsToUse[bond->getIdx()]) {
    // RDKit❗❌:       if (QueryOps::hasBondTypeQuery(*bond)) {
    // RDKit❗❌:         // we don't kekulize bonds with bond type queries
    // RDKit❗❌:         bondsToUse[bond->getIdx()] = 0;
    // RDKit❗❌:       } else if (bond->getIsAromatic()) {
    // RDKit❗❌:         foundAromatic = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // before everything do implicit valence calculation and store them
    // RDKit❗❌:   // we will repeat after kekulization and compare for the sake of error
    // RDKit❗❌:   // checking
    // RDKit❗❌:   auto numAtoms = mol.getNumAtoms();
    // RDKit❗❌:   INT_VECT valences(numAtoms);
    // RDKit❗❌:   boost::dynamic_bitset<> dummyAts(numAtoms);
    // RDKit❗❌:
    // RDKit❗❌:   for (auto atom : mol.atoms()) {
    // RDKit❗❌:     if (!atomsToUse[atom->getIdx()]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     atom->calcImplicitValence(false);
    // RDKit❗❌:     valences[atom->getIdx()] = atom->getTotalValence();
    // RDKit❗❌:     if (isAromaticAtom(*atom)) {
    // RDKit❗❌:       foundAromatic = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!atom->getAtomicNum()) {
    // RDKit❗❌:       dummyAts[atom->getIdx()] = 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!foundAromatic) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   UINT_VECT atomRanks(mol.getNumAtoms());
    // RDKit❗❌:   if (canonical) {
    // RDKit❗❌:     Canon::rankFragmentAtoms(mol, atomRanks, atomsToUse, bondsToUse);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     // When canonical=false (e.g. during sanitization), we skip the
    // RDKit❗❌:     // expensive ranking step and use atom indices directly.  This is
    // RDKit❗❌:     // appropriate because sanitization runs *before* stereo perception:
    // RDKit❗❌:     // canonical ranking would be based on incomplete chemistry and the
    // RDKit❗❌:     // "deterministic" result would be meaningless.  Callers who need a
    // RDKit❗❌:     // canonical Kekulé form should call Kekulize() with canonical=true
    // RDKit❗❌:     // after the molecule is fully sanitized and stereo has been assigned.
    // RDKit❗❌:     std::iota(atomRanks.begin(), atomRanks.end(), 0u);
    // RDKit❗❌:   }
    // RDKit❗❌:   // if any bonds to kekulize then give it a try:
    // RDKit❗❌:   if (bondsToUse.any()) {
    // RDKit❗❌:     // mark atoms at the beginning of wedged bonds
    // RDKit❗❌:     boost::dynamic_bitset<> wedgedAtoms(numAtoms);
    // RDKit❗❌:     for (const auto bond : mol.bonds()) {
    // RDKit❗❌:       if (bondsToUse[bond->getIdx()] &&
    // RDKit❗❌:           (bond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗❌:            bond->getBondDir() == Bond::BEGINDASH)) {
    // RDKit❗❌:         wedgedAtoms.set(bond->getBeginAtomIdx());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // A bit on the state of the molecule at this point
    // RDKit❗❌:     // - aromatic and non aromatic atoms and bonds may be mixed up
    // RDKit❗❌:
    // RDKit❗❌:     // - for all aromatic bonds it is assumed that that both the following
    // RDKit❗❌:     //   are true:
    // RDKit❗❌:     //       - getIsAromatic returns true
    // RDKit❗❌:     //       - getBondType return aromatic
    // RDKit❗❌:     // - all aromatic atoms return true for "getIsAromatic"
    // RDKit❗❌:
    // RDKit❗❌:     // first find all the simple rings in the molecule that are not
    // RDKit❗❌:     // completely composed of dummy atoms
    // RDKit❗❌:     VECT_INT_VECT allringsSSSR;
    // RDKit❗❌:     if (!mol.getRingInfo()->isInitialized()) {
    // RDKit❗❌:       MolOps::findSSSR(mol, allringsSSSR);
    // RDKit❗❌:     }
    // RDKit❗❌:     const VECT_INT_VECT &allrings =
    // RDKit❗❌:         allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
    // RDKit❗❌:     std::deque<INT_VECT> tmpRings;
    // RDKit❗❌:     auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
    // RDKit❗❌:       bool ringOk = false;
    // RDKit❗❌:       for (auto ai : ring) {
    // RDKit❗❌:         if (!atomsToUse[ai]) {
    // RDKit❗❌:           return false;
    // RDKit❗❌:         }
    // RDKit❗❌:         if (!dummyAts[ai]) {
    // RDKit❗❌:           ringOk = true;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       return ringOk;
    // RDKit❗❌:     };
    // RDKit❗❌:     // we can't just copy the rings over: we're going to rearrange them so that
    // RDKit❗❌:     // we try to favor starting the traversal of any ring from an atom that is
    // RDKit❗❌:     // at the end of a wedged ring bond. This is part of our attempt to avoid
    // RDKit❗❌:     // assigning double bonds to bonds with wedging
    // RDKit❗❌:     for (const auto &ring : allrings) {
    // RDKit❗❌:       if (containsNonDummy(ring)) {
    // RDKit❗❌:         unsigned int startPos = 0;
    // RDKit❗❌:         bool hasWedge = false;
    // RDKit❗❌:         for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit❗❌:           if (wedgedAtoms[ring[ri]]) {
    // RDKit❗❌:             startPos = ri;
    // RDKit❗❌:             hasWedge = true;
    // RDKit❗❌:             break;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         INT_VECT nring(ring.size());
    // RDKit❗❌:         for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit❗❌:           nring[ri] = ring.at((ri + startPos) % ring.size());
    // RDKit❗❌:         }
    // RDKit❗❌:         if (!hasWedge) {
    // RDKit❗❌:           tmpRings.push_back(nring);
    // RDKit❗❌:         } else {
    // RDKit❗❌:           tmpRings.push_front(nring);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     VECT_INT_VECT arings;
    // RDKit❗❌:     arings.reserve(allrings.size());
    // RDKit❗❌:     arings.insert(arings.end(), tmpRings.begin(), tmpRings.end());
    // RDKit❗❌:     VECT_INT_VECT allbrings;
    // RDKit❗❌:     RingUtils::convertToBonds(arings, allbrings, mol);
    // RDKit❗❌:     VECT_INT_VECT brings;
    // RDKit❗❌:     brings.reserve(allbrings.size());
    // RDKit❗❌:     auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
    // RDKit❗❌:       return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
    // RDKit❗❌:         return bondsToUse[bi];
    // RDKit❗❌:       });
    // RDKit❗❌:     };
    // RDKit❗❌:     VECT_INT_VECT aringsRemaining;
    // RDKit❗❌:     aringsRemaining.reserve(arings.size());
    // RDKit❗❌:     for (unsigned i = 0; i < allbrings.size(); ++i) {
    // RDKit❗❌:       if (copyBondRingsWithinFragment(allbrings[i])) {
    // RDKit❗❌:         brings.push_back(allbrings[i]);
    // RDKit❗❌:         aringsRemaining.push_back(arings[i]);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     arings = std::move(aringsRemaining);
    // RDKit❗❌:
    // RDKit❗❌:     // make a neighbor map for the rings i.e. a ring is a
    // RDKit❗❌:     // neighbor to another candidate ring if it shares at least
    // RDKit❗❌:     // one bond
    // RDKit❗❌:     // useful to figure out fused systems
    // RDKit❗❌:     INT_INT_VECT_MAP neighMap;
    // RDKit❗❌:     RingUtils::makeRingNeighborMap(brings, neighMap);
    // RDKit❗❌:
    // RDKit❗❌:     int curr = 0;
    // RDKit❗❌:     int cnrs = rdcast<int>(arings.size());
    // RDKit❗❌:     boost::dynamic_bitset<> fusDone(cnrs);
    // RDKit❗❌:     while (curr < cnrs) {
    // RDKit❗❌:       INT_VECT fused;
    // RDKit❗❌:       RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
    // RDKit❗❌:       VECT_INT_VECT frings(fused.size());
    // RDKit❗❌:       std::transform(fused.begin(), fused.end(), frings.begin(),
    // RDKit❗❌:                      [&arings](const int ri) { return arings[ri]; });
    // RDKit❗❌:       kekulizeFused(mol, frings, atomRanks, maxBackTracks);
    // RDKit❗❌:       int rix;
    // RDKit❗❌:       for (rix = 0; rix < cnrs; ++rix) {
    // RDKit❗❌:         if (!fusDone[rix]) {
    // RDKit❗❌:           curr = rix;
    // RDKit❗❌:           break;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (rix == cnrs) {
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (markAtomsBonds) {
    // RDKit❗❌:     // if we want the atoms and bonds to be marked non-aromatic do
    // RDKit❗❌:     // that here.
    // RDKit❗❌:     if (!mol.getRingInfo()->isInitialized()) {
    // RDKit❗❌:       MolOps::findSSSR(mol);
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto bond : mol.bonds()) {
    // RDKit❗❌:       if (bondsToUse[bond->getIdx()]) {
    // RDKit❗❌:         bond->setIsAromatic(false);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     for (auto atom : mol.atoms()) {
    // RDKit❗❌:       if (atomsToUse[atom->getIdx()] && atom->getIsAromatic()) {
    // RDKit❗❌:         // if we're doing the full molecule and there are aromatic atoms not in
    // RDKit❗❌:         // a ring, throw an exception
    // RDKit❗❌:         if (atomsToUse.all() && bondsToUse.all() &&
    // RDKit❗❌:             !mol.getRingInfo()->numAtomRings(atom->getIdx())) {
    // RDKit❗❌:           std::ostringstream errout;
    // RDKit❗❌:           errout << "non-ring atom " << atom->getIdx() << " marked aromatic";
    // RDKit❗❌:           auto msg = errout.str();
    // RDKit❗❌:           BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit❗❌:           throw AtomKekulizeException(msg, atom->getIdx());
    // RDKit❗❌:         }
    // RDKit❗❌:         atom->setIsAromatic(false);
    // RDKit❗❌:         // make sure "explicit" Hs on things like pyrroles don't hang around
    // RDKit❗❌:         // this was Github Issue 141
    // RDKit❗❌:         if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
    // RDKit❗❌:             atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
    // RDKit❗❌:           atom->setNoImplicit(false);
    // RDKit❗❌:           atom->setNumExplicitHs(0);
    // RDKit❗❌:           atom->updatePropertyCache(false);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // ok some error checking here force a implicit valence
    // RDKit❗❌:   // calculation that should do some error checking by itself. In
    // RDKit❗❌:   // addition compare them to what they were before kekulizing
    // RDKit❗❌:   for (auto atom : mol.atoms()) {
    // RDKit❗❌:     if (!atomsToUse[atom->getIdx()]) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     int val = atom->getTotalValence();
    // RDKit❗❌:     if (val != valences[atom->getIdx()]) {
    // RDKit❗❌:       std::ostringstream errout;
    // RDKit❗❌:       errout << "Kekulization somehow screwed up valence on " << atom->getIdx()
    // RDKit❗❌:              << ": " << val << "!=" << valences[atom->getIdx()] << std::endl;
    // RDKit❗❌:       auto msg = errout.str();
    // RDKit❗❌:       BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit❗❌:       throw AtomKekulizeException(msg, atom->getIdx());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment (2026.03.6 complete)
    // ❗❌ qualification applies to unchanged detached validation, Vec<bool>
    // masks and candidate copying vs source packed bitsets, not silent fallback.
    // Graph, valence rows and live ring state are mutated at each source store.
    let prepared =
        prepare_kekulize_inputs(topology, atoms_in_play, bonds_in_play, query_state, valence)?;
    if !prepared.found_aromatic {
        return Ok((Vec::new(), false));
    }
    let effective_bonds_any = prepared.bonds_in_play.iter().any(|selected| *selected);
    let (replaced_by_reset, atom_ranks) = kekulize_ring_state_transition_mut(
        topology,
        valence,
        &prepared.atoms_in_play,
        &prepared.bonds_in_play,
        params.canonical,
        effective_bonds_any,
        rows,
    )?;

    let (candidate_atom_rings, candidate_bond_rings) = if effective_bonds_any {
        match rows.as_ring_info() {
            Some(effective) => {
                validate_consumed_ring_dimensions(topology, effective)?;
                #[cfg(test)]
                ring_transport_probe::record_borrow(effective);
                collect_kekulize_candidate_rings(
                    topology,
                    effective,
                    &prepared.atoms_in_play,
                    &prepared.bonds_in_play,
                    &prepared.dummy_atoms,
                    &prepared.wedged_atoms,
                )?
            }
            None => (Vec::new(), Vec::new()),
        }
    } else {
        (Vec::new(), Vec::new())
    };
    if effective_bonds_any && !candidate_atom_rings.is_empty() {
        let fused = kekulize_fused_components_mut(
            topology,
            &candidate_atom_rings,
            &candidate_bond_rings,
            rows.as_ring_info()
                .expect("candidate rows imply initialized ring state"),
            valence,
            &atom_ranks,
            params.max_backtracks,
        )?;
        if !fused.succeeded {
            return Err(KekulizeError::NotKekulizable {
                problem_atoms: fused.problem_atoms,
            });
        }
    }

    let mut refreshed_valence_atoms = Vec::new();
    if params.mark_atoms_bonds {
        if rows.as_ring_info().is_none() {
            #[cfg(test)]
            ring_transport_probe::record_sssr();
            rows.install(acquire_sssr_rows(topology)?);
        }
        let marking_rows = rows.as_ring_info().expect("marking rows installed");
        validate_consumed_ring_dimensions(topology, marking_rows)?;
        for (bond_idx, selected) in prepared.bonds_in_play.iter().copied().enumerate() {
            if selected {
                topology.bonds[bond_idx].set_aromatic(false);
            }
        }
        for (atom_idx, selected) in prepared.atoms_in_play.iter().copied().enumerate() {
            if !selected || !topology.atoms[atom_idx].is_aromatic() {
                continue;
            }
            let atom_id = AtomId::new(atom_idx);
            if prepared.atoms_in_play.iter().all(|selected| *selected)
                && prepared.bonds_in_play.iter().all(|selected| *selected)
                && marking_rows.num_atom_rings(atom_id) == 0
            {
                return Err(KekulizeError::AromaticAtomOutsideRing { atom: atom_id });
            }
            topology.atoms[atom_idx].set_aromatic(false);
            if matches!(topology.atoms[atom_idx].atomic_number(), 7 | 15)
                && topology.atoms[atom_idx].formal_charge() == 0
                && topology.atoms[atom_idx].explicit_hydrogens() == 1
            {
                topology.atoms[atom_idx].set_no_implicit(false);
                topology.atoms[atom_idx].set_explicit_hydrogens(0);
                valence.explicit_valence[atom_idx] =
                    crate::assign_explicit_valence_for_atom_from_parts(
                        &topology.atoms,
                        &topology.bonds,
                        &topology.adjacency,
                        atom_id,
                        false,
                    )?;
                crate::valence::source_calc_implicit_cache_row(
                    &topology.atoms,
                    &topology.bonds,
                    &topology.adjacency,
                    atom_id,
                    &mut valence.explicit_valence[atom_idx],
                    &mut valence.implicit_hydrogens[atom_idx],
                    false,
                )?;
                refreshed_valence_atoms.push(atom_id);
            }
        }
    }

    for (atom_idx, selected) in prepared.atoms_in_play.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom_id = AtomId::new(atom_idx);
        let after = checked_total_valence(valence, &topology.atoms[atom_idx])?;
        let before = prepared.original_total_valences[atom_idx];
        if after != before {
            return Err(KekulizeError::PostconditionValenceMismatch {
                atom: atom_id,
                before,
                after,
            });
        }
    }
    topology.validate()?;
    Ok((refreshed_valence_atoms, replaced_by_reset))
}

pub fn kekulize(
    topology: &TopologyBlock,
    params: &KekulizeParams,
) -> Result<KekulizeAssignment, KekulizeError> {
    kekulize_with_query_state(topology, params, None)
}

/// Kekulizes a whole topology while transporting the caller's ring state.
///
/// `rings` supplies the caller's current ring assignment when one exists.
/// The returned `ring_update` field reports the final state: `None` means
/// the supplied state is unchanged (absent, already reset, or an initialized
/// Fast/SSSR/Symm assignment that was borrowed without modification);
/// `Some(initialized)` moves the newly acquired final state; and
/// `Some(uninitialized)` replaces a previously initialized Other assignment
/// that canonical ranking reset. The supplied assignment is immutable on
/// success and error.
pub fn kekulize_with_query_state_and_ring_info(
    topology: &TopologyBlock,
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
    rings: Option<&RingInfo>,
    source_valence: Option<&ValenceAssignment>,
) -> Result<KekulizeAssignment, KekulizeError> {
    if let Some(state) = query_state {
        state.validate_for_topology(topology)?;
    }
    let atoms_in_play = vec![true; topology.atoms.len()];
    let bonds_in_play = vec![true; topology.bonds.len()];
    kekulize_fragment(
        topology,
        &atoms_in_play,
        &bonds_in_play,
        params,
        query_state,
        rings,
        source_valence,
    )
}

/// One source attempt over borrowed detached values; completed stores survive Err.
/// This carries no live Molecule/runtime authority and performs no recovery.
#[doc(hidden)]
pub fn source_kekulize_attempt(
    topology: &mut TopologyBlock,
    valence: &mut ValenceAssignment,
    rings: &mut RingInfo,
    params: &KekulizeParams,
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::Kekulize (2026.03.6 complete)
    // RDKit✔️✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit✔️✔️:               unsigned int maxBackTracks) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
    // RDKit✔️✔️:   atomsToUse.set();
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
    // RDKit✔️✔️:   bondsToUse.set();
    // RDKit✔️✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
    // RDKit✔️✔️:                             canonical, maxBackTracks);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::Kekulize (2026.03.6 complete)
    let atoms = vec![true; topology.atoms.len()];
    let bonds = vec![true; topology.bonds.len()];
    let mut rows = KekulizeRingRows::Live(rings);
    kekulize_fragment_attempt(topology, valence, &mut rows, &atoms, &bonds, params, None)
        .map(|_| ())
}

fn source_sanitize_kekulize_error(error: &KekulizeError) -> bool {
    matches!(
        error,
        KekulizeError::NotKekulizable { .. }
            | KekulizeError::AromaticAtomOutsideRing { .. }
            | KekulizeError::PostconditionValenceMismatch { .. }
            | KekulizeError::Valence(ValenceError::InvalidValence { .. })
    )
}

fn restore_source_aromatic_membership(
    topology: &mut TopologyBlock,
    aromatic_atoms: &[bool],
    aromatic_bonds: &[bool],
) {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeIfPossible aromatic recovery block (2026.03.6 complete)
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit✔️✔️:       if (aromaticBonds[i]) {
    // RDKit✔️✔️:         auto bond = mol.getBondWithIdx(i);
    // RDKit✔️✔️:         bond->setIsAromatic(true);
    // RDKit✔️✔️:         bond->setBondType(Bond::BondType::AROMATIC);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       if (aromaticAtoms[i]) {
    // RDKit✔️✔️:         mol.getAtomWithIdx(i)->setIsAromatic(true);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeIfPossible aromatic recovery block (2026.03.6 complete)
    // The complete source recovery owner is source_kekulize_for_sanitize below.
    // Deliberately restore AROMATIC order, not the original single/double order.
    for (bond, &was_aromatic) in topology.bonds.iter_mut().zip(aromatic_bonds) {
        if was_aromatic {
            bond.set_aromatic(true);
            bond.set_order(BondOrder::Aromatic);
        }
    }
    for (atom, &was_aromatic) in topology.atoms.iter_mut().zip(aromatic_atoms) {
        if was_aromatic {
            atom.set_aromatic(true);
        }
    }
}

pub(crate) fn source_kekulize_for_sanitize(
    topology: &mut TopologyBlock,
    valence: &mut ValenceAssignment,
    rings: &mut Option<RingInfo>,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: kekulizeForSanitize (2026.03.6 complete)
    // RDKit✔️✔️: void kekulizeForSanitize(RWMol &mol) {
    // RDKit✔️✔️:   if (!MolOps::KekulizeIfPossible(mol, true, false)) {
    // RDKit✔️✔️:     MolOps::Kekulize(mol, true, true);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/MolOps.cpp :: kekulizeForSanitize (2026.03.6 complete)
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::KekulizeIfPossible (2026.03.6 complete)
    // RDKit✔️✔️: bool KekulizeIfPossible(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit✔️✔️:                         unsigned int maxBackTracks) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromaticBonds(mol.getNumBonds());
    // RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getIsAromatic()) {
    // RDKit✔️✔️:       aromaticBonds.set(bond->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromaticAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (isAromaticAtom(*atom)) {
    // RDKit✔️✔️:       aromaticAtoms.set(atom->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool res = true;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     Kekulize(mol, markAtomsBonds, canonical, maxBackTracks);
    // RDKit✔️✔️:   } catch (const MolSanitizeException &) {
    // RDKit✔️✔️:     res = false;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit✔️✔️:       if (aromaticBonds[i]) {
    // RDKit✔️✔️:         auto bond = mol.getBondWithIdx(i);
    // RDKit✔️✔️:         bond->setIsAromatic(true);
    // RDKit✔️✔️:         bond->setBondType(Bond::BondType::AROMATIC);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       if (aromaticAtoms[i]) {
    // RDKit✔️✔️:         mol.getAtomWithIdx(i)->setIsAromatic(true);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::KekulizeIfPossible (2026.03.6 complete)
    let aromatic_atoms = topology
        .atoms
        .iter()
        .map(|a| atom_is_aromatic_for_kekulize(topology, a.id()))
        .collect::<Vec<_>>();
    let aromatic_bonds = topology
        .bonds
        .iter()
        .map(Bond::is_aromatic)
        .collect::<Vec<_>>();
    let atoms = vec![true; topology.atoms.len()];
    let bonds = vec![true; topology.bonds.len()];
    let was_absent = rings.is_none();
    let live = rings.get_or_insert_with(|| {
        let mut r = RingInfo::new(RingFindType::OtherOrUnknown, 0, 0);
        r.reset();
        r
    });
    let mut rows = KekulizeRingRows::Live(live);
    let first = kekulize_fragment_attempt(
        topology,
        valence,
        &mut rows,
        &atoms,
        &bonds,
        &KekulizeParams {
            mark_atoms_bonds: true,
            canonical: false,
            max_backtracks: 100,
        },
        query_state,
    );
    let result = match first {
        Ok(_) => Ok(()),
        Err(error) if source_sanitize_kekulize_error(&error) => {
            restore_source_aromatic_membership(topology, &aromatic_atoms, &aromatic_bonds);
            kekulize_fragment_attempt(
                topology,
                valence,
                &mut rows,
                &atoms,
                &bonds,
                &KekulizeParams {
                    mark_atoms_bonds: true,
                    canonical: true,
                    max_backtracks: 100,
                },
                query_state,
            )
            .map(|_| ())
        }
        Err(error) => Err(error),
    };
    // Option is detached transport only: do not manufacture an initialized
    // carrier for a source-uninitialized early return; actual Fast/SSSR survives.
    let remains_uninitialized = rows.as_ring_info().is_none();
    if was_absent && remains_uninitialized {
        *rings = None;
    }
    result
}

/// Kekulizes selected original-index atoms and bonds in a detached topology.
pub fn kekulize_selected_fragment(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    params: &KekulizeParams,
) -> Result<KekulizeAssignment, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct selected kekulization
    // RDKit✔️❌:     if (atomsInPlay && bondsInPlay) {
    // RDKit✔️❌:       MolOps::details::KekulizeFragment(static_cast<RWMol &>(mol), *atomsInPlay,
    // RDKit✔️❌:                                         *bondsInPlay);
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       MolOps::Kekulize(static_cast<RWMol &>(mol));
    // RDKit✔️❌:     }
    // END RDKIT CPP FUNCTION SmilesWrite::FragmentSmilesConstruct selected kekulization
    // Behavior: the caller supplies both original-index masks, so the single
    // core masked owner performs source validation and all selected chemistry.
    // Cost: detached return ownership clones the topology, including the
    // source early-return path; the upstream caller mutates its private copy.
    kekulize_fragment(
        topology,
        atoms_in_play,
        bonds_in_play,
        params,
        None,
        None,
        None,
    )
}

#[doc(hidden)]
pub fn kekulize_with_query_state(
    topology: &TopologyBlock,
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<KekulizeAssignment, KekulizeError> {
    if let Some(state) = query_state {
        state.validate_for_topology(topology)?;
    }
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::Kekulize
    // RDKit✔️✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit✔️✔️:               unsigned int maxBackTracks) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
    // RDKit✔️✔️:   atomsToUse.set();
    // RDKit✔️✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
    // RDKit✔️✔️:   bondsToUse.set();
    let atoms_in_play = vec![true; topology.atoms.len()];
    let bonds_in_play = vec![true; topology.bonds.len()];
    // RDKit✔️✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
    // RDKit✔️✔️:                             canonical, maxBackTracks);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::Kekulize
    kekulize_fragment(
        topology,
        &atoms_in_play,
        &bonds_in_play,
        params,
        query_state,
        None,
        None,
    )
}

pub fn kekulize_if_possible(
    topology: &TopologyBlock,
    params: &KekulizeParams,
) -> Result<KekulizeAttempt, KekulizeError> {
    kekulize_if_possible_with_query_state(topology, params, None)
}

#[doc(hidden)]
pub fn kekulize_if_possible_with_query_state(
    topology: &TopologyBlock,
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<KekulizeAttempt, KekulizeError> {
    kekulize_if_possible_with_query_state_and_ring_info(topology, params, query_state, None)
}

/// Source IfPossible fallback over the existing ring-aware Kekulize owner.
/// Supplied rings are immutable; Applied transports the owner's ring_update.
/// The existing NotKekulizable attempt does not transport exception-time ring
/// updates: failure-side upstream RingInfo parity remains qualified.
pub fn kekulize_if_possible_with_query_state_and_ring_info(
    topology: &TopologyBlock,
    params: &KekulizeParams,
    query_state: Option<QueryStateRef<'_>>,
    rings: Option<&RingInfo>,
) -> Result<KekulizeAttempt, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::KekulizeIfPossible
    // RDKit✔️✔️: bool KekulizeIfPossible(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit✔️✔️:                         unsigned int maxBackTracks) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromaticBonds(mol.getNumBonds());
    // RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getIsAromatic()) {
    // RDKit✔️✔️:       aromaticBonds.set(bond->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   boost::dynamic_bitset<> aromaticAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (isAromaticAtom(*atom)) {
    // RDKit✔️✔️:       aromaticAtoms.set(atom->getIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool res = true;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     Kekulize(mol, markAtomsBonds, canonical, maxBackTracks);
    // RDKit✔️✔️:   } catch (const MolSanitizeException &) {
    // RDKit✔️✔️:     res = false;
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit✔️✔️:       if (aromaticBonds[i]) {
    // RDKit✔️✔️:         auto bond = mol.getBondWithIdx(i);
    // RDKit✔️✔️:         bond->setIsAromatic(true);
    // RDKit✔️✔️:         bond->setBondType(Bond::BondType::AROMATIC);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       if (aromaticAtoms[i]) {
    // RDKit✔️✔️:         mol.getAtomWithIdx(i)->setIsAromatic(true);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: MolOps::KekulizeIfPossible
    // The detached input is immutable, so its complete state is the source
    // snapshot and is returned unchanged for the audited sanitize failures.
    // Input review: forward all four inputs unchanged to ONE existing owner.
    // Behavior review: retain all modeled fallback categories and problem order;
    // all other typed errors propagate. Failure ring-update parity is qualified.
    // Cost review: no extra ranking, acquisition, or success topology/ring clone;
    // Applied moves its complete assignment, fallback retains its original clone.
    #[cfg(test)]
    drawing_ring_if_possible_probe::forward(topology, params, query_state.is_some(), rings);
    match kekulize_with_query_state_and_ring_info(topology, params, query_state, rings, None) {
        Ok(assignment) => Ok(KekulizeAttempt::Applied(assignment)),
        Err(KekulizeError::NotKekulizable { problem_atoms }) => {
            Ok(KekulizeAttempt::NotKekulizable {
                topology: topology.clone(),
                problem_atoms,
            })
        }
        Err(KekulizeError::AromaticAtomOutsideRing { atom })
        | Err(KekulizeError::PostconditionValenceMismatch { atom, .. }) => {
            Ok(KekulizeAttempt::NotKekulizable {
                topology: topology.clone(),
                problem_atoms: vec![atom],
            })
        }
        Err(error) => Err(error),
    }
}

#[cfg(test)]
mod drawing_ring_if_possible_probe {
    use super::*;
    use std::cell::RefCell;
    #[derive(Clone, Debug, PartialEq, Eq)]
    pub(super) struct Forward {
        pub topology: usize,
        pub params: KekulizeParams,
        pub query_present: bool,
        pub rings: Option<usize>,
    }
    thread_local! {
        static CALLS: RefCell<Vec<Forward>> = const { RefCell::new(Vec::new()) };
    }
    pub(super) fn forward(
        topology: &TopologyBlock,
        params: &KekulizeParams,
        query_present: bool,
        rings: Option<&RingInfo>,
    ) {
        CALLS.with(|calls| {
            calls.borrow_mut().push(Forward {
                topology: topology as *const _ as usize,
                params: params.clone(),
                query_present,
                rings: rings.map(|rings| rings as *const _ as usize),
            })
        });
    }
    pub(super) fn calls() -> Vec<Forward> {
        CALLS.with(|calls| calls.borrow().clone())
    }
}

struct CanonRankReadView<'a> {
    atoms: &'a [Atom],
    bonds: &'a [Bond],
    adjacency: &'a AdjacencyList,
    rings: Cow<'a, RingInfo>,
    valence: Cow<'a, ValenceAssignment>,
}

impl<'a> CanonRankReadView<'a> {
    fn from_topology(topology: &'a TopologyBlock) -> Result<Self, CanonicalRankError> {
        Self::from_parts(&topology.atoms, &topology.bonds, &topology.adjacency)
    }

    fn from_parts(
        atoms: &'a [Atom],
        bonds: &'a [Bond],
        adjacency: &'a AdjacencyList,
    ) -> Result<Self, CanonicalRankError> {
        // RDKit✔️✔️:   bool clearRings = false;
        // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
        // RDKit✔️✔️:     MolOps::fastFindRings(mol);
        // RDKit✔️✔️:     clearRings = true;
        // RDKit✔️✔️:   }
        let rings = fast_find_rings_from_parts(atoms.len(), bonds, adjacency)?;
        let valence = crate::assign_valence_with_options_from_parts(
            atoms,
            bonds,
            adjacency,
            ValenceModel::RdkitLike,
            false,
        )?;
        Ok(Self {
            atoms,
            bonds,
            adjacency,
            rings: Cow::Owned(rings),
            valence: Cow::Owned(valence),
        })
    }

    fn from_prepared_state(
        topology: &'a TopologyBlock,
        valence: &'a ValenceAssignment,
        rings: Option<&'a RingInfo>,
        atoms_in_play: Option<&[bool]>,
    ) -> Result<Self, CanonicalRankError> {
        // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
        // RDKit✔️✔️:     MolOps::fastFindRings(mol);
        // RDKit✔️✔️:     clearRings = true;
        // RDKit✔️✔️:   }
        // A caller supplies already prepared property-cache values for this
        // exact topology. Check dimensions and source per-atom validity before
        // borrowing; a missing/unknown ring set uses the one existing fast
        // owner. Borrowing is O(1) and avoids duplicate O(V+E) preparation.
        if valence.explicit_valence.len() != topology.atoms.len()
            || valence.implicit_hydrogens.len() != topology.atoms.len()
        {
            return Err(CanonicalRankError::PreparedValenceLength {
                atom_count: topology.atoms.len(),
                explicit_len: valence.explicit_valence.len(),
                implicit_len: valence.implicit_hydrogens.len(),
            });
        }
        // RDKit✔️✔️:     if (atomsInPlay[i]) {
        // RDKit✔️✔️:         advancedInitCanonAtom(mol, atomsi, i);
        // Only selected fragment rows consume scalar fields. Unselected source
        // -1 sentinels stay untouched; ordinary whole views still check all rows.
        for (atom_index, atom) in topology.atoms.iter().enumerate() {
            if atoms_in_play.is_some_and(|mask| !mask[atom_index]) {
                continue;
            }
            if crate::valence::cached_explicit_valence(atom, Some(valence)).is_err()
                || crate::hcount::implicit_hydrogen_count(atom, valence).is_err()
            {
                return Err(CanonicalRankError::PreparedValenceInvalid { atom_index });
            }
        }
        if let Some(rings) = rings
            && (rings.atom_row_count() != topology.atoms.len()
                || rings.bond_row_count() != topology.bonds.len())
        {
            return Err(CanonicalRankError::PreparedRingLength {
                expected_atoms: topology.atoms.len(),
                actual_atoms: rings.atom_row_count(),
                expected_bonds: topology.bonds.len(),
                actual_bonds: rings.bond_row_count(),
            });
        }
        let rings = if let Some(rings) = rings.filter(|rings| rings.is_find_fast_or_better()) {
            Cow::Borrowed(rings)
        } else {
            Cow::Owned(fast_find_rings_from_parts(
                topology.atoms.len(),
                &topology.bonds,
                &topology.adjacency,
            )?)
        };
        Ok(Self {
            atoms: &topology.atoms,
            bonds: &topology.bonds,
            adjacency: &topology.adjacency,
            rings,
            valence: Cow::Borrowed(valence),
        })
    }

    fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    fn atom_degree(&self, atom: AtomId) -> usize {
        self.adjacency.neighbors_of(atom.index()).len()
    }

    fn atom_neighbors(&self, atom: AtomId) -> Vec<usize> {
        self.adjacency
            .neighbors_of(atom.index())
            .iter()
            .map(|neighbor| neighbor.atom_index)
            .collect()
    }

    fn bond_other_atom_index(&self, bond_id: BondId, atom_id: AtomId) -> Option<usize> {
        let bond = self.bonds.get(bond_id.index())?;
        if bond.begin() == atom_id {
            Some(bond.end().index())
        } else if bond.end() == atom_id {
            Some(bond.begin().index())
        } else {
            None
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct CanonicalRankParams {
    pub break_ties: bool,
    pub include_chirality: bool,
    pub include_isotopes: bool,
    pub include_atom_maps: bool,
    pub include_chiral_presence: bool,
    pub include_stereo_groups: bool,
    pub use_non_stereo_ranks: bool,
    pub include_ring_stereo: bool,
    chirality_rings_use_ring_stereo: bool,
}

impl Default for CanonicalRankParams {
    fn default() -> Self {
        Self {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: true,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
            chirality_rings_use_ring_stereo: true,
        }
    }
}

impl CanonicalRankParams {
    const fn kekulize_fragment_default() -> Self {
        Self {
            break_ties: true,
            include_chirality: true,
            include_isotopes: true,
            include_atom_maps: true,
            include_chiral_presence: false,
            include_stereo_groups: false,
            use_non_stereo_ranks: false,
            include_ring_stereo: true,
            // rankFragmentAtoms sets df_useChiralityRings directly from
            // includeChirality instead of gating it on includeRingStereo.
            chirality_rings_use_ring_stereo: false,
        }
    }
}

pub(crate) fn rank_mol_atoms(topology: &TopologyBlock) -> Result<Vec<usize>, CanonicalRankError> {
    rank_mol_atoms_with_params(topology, &CanonicalRankParams::default())
}

/// Returns RDKit-compatible canonical ranks for all atoms with explicit source options.
pub fn rank_mol_atoms_with_params(
    topology: &TopologyBlock,
    params: &CanonicalRankParams,
) -> Result<Vec<usize>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankMolAtoms
    // RDKit✔️✔️: void rankMolAtoms(const ROMol &mol, std::vector<unsigned int> &res,
    // RDKit✔️✔️:                   bool breakTies, bool includeChirality, bool includeIsotopes,
    // RDKit✔️✔️:                   bool includeAtomMaps, bool includeChiralPresence,
    // RDKit✔️✔️:                   bool includeStereoGroups, bool useNonStereoRanks,
    // RDKit✔️✔️:                   bool includeRingStereo) {
    // RDKit✔️✔️:   if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    if topology.atoms.is_empty() {
        return Ok(Vec::new());
    }
    // RDKit✔️✔️:   bool clearRings = false;
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:     clearRings = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res.resize(mol.getNumAtoms());
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
    // RDKit✔️✔️:   initCanonAtoms(mol, atoms, includeChirality, includeStereoGroups);
    let view = CanonRankReadView::from_topology(topology)?;
    let mut atoms = init_canon_atoms(
        &view,
        topology,
        params.include_chirality,
        params.include_stereo_groups,
    )?;
    // RDKit✔️✔️:   AtomCompareFunctor ftor(&atoms.front(), mol);
    // RDKit✔️✔️:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit✔️✔️:   ftor.df_useChirality = includeChirality;
    // RDKit✔️✔️:   ftor.df_useChiralityRings = includeChirality && includeRingStereo;
    // RDKit✔️✔️:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit✔️✔️:   ftor.df_useNonStereoRanks = useNonStereoRanks;
    // RDKit✔️✔️:   ftor.df_useChiralPresence = includeChiralPresence;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::vector<int> order(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit✔️✔️:                           includeRingStereo);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     res[order[i]] = atoms[order[i]].index;
    // RDKit✔️✔️:   }
    let result = rank_initialized_atoms(&view, &mut atoms, *params)?;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (clearRings) {
    // RDKit✔️✔️:     mol.getRingInfo()->reset();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }  // end of rankMolAtoms()
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankMolAtoms
    Ok(result)
}

/// Returns RDKit-compatible stable canonical ranks for a selected graph fragment.
///
/// The result has one entry per topology atom, matching `rankFragmentAtoms`.
/// `atoms_in_play` and `bonds_in_play` are independent source masks; a bond is
/// included in the fragment only when its bit and both endpoint bits are set.
pub fn rank_fragment_atoms(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
) -> Result<Vec<usize>, CanonicalRankError> {
    rank_fragment_atoms_with_params(
        topology,
        atoms_in_play,
        bonds_in_play,
        None,
        None,
        &CanonicalRankParams::kekulize_fragment_default(),
    )
}

/// Ranks an original-index fragment with source symbol tables and options.
pub fn rank_fragment_atoms_with_params(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    params: &CanonicalRankParams,
) -> Result<Vec<usize>, CanonicalRankError> {
    rank_fragment_atoms_engine(
        topology,
        atoms_in_play,
        bonds_in_play,
        atom_symbols,
        bond_symbols,
        params,
        None,
    )
}

/// Ranks a fragment while borrowing exact-topology prepared valence and rings.
pub fn rank_fragment_atoms_with_prepared_state(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: Option<&RingInfo>,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    params: &CanonicalRankParams,
) -> Result<Vec<usize>, CanonicalRankError> {
    rank_fragment_atoms_engine(
        topology,
        atoms_in_play,
        bonds_in_play,
        atom_symbols,
        bond_symbols,
        params,
        Some((valence, rings)),
    )
}

fn rank_fragment_atoms_engine(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    atom_symbols: Option<&[String]>,
    bond_symbols: Option<&[String]>,
    params: &CanonicalRankParams,
    prepared: Option<(&ValenceAssignment, Option<&RingInfo>)>,
) -> Result<Vec<usize>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankFragmentAtoms
    // RDKit✔️✔️: void rankFragmentAtoms(const ROMol &mol, std::vector<unsigned int> &res,
    // RDKit✔️✔️:                        const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️✔️:                        const boost::dynamic_bitset<> &bondsInPlay,
    // RDKit✔️✔️:                        const std::vector<std::string> *atomSymbols,
    // RDKit✔️✔️:                        const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                        bool breakTies, bool includeChirality,
    // RDKit✔️✔️:                        bool includeIsotopes, bool includeAtomMaps,
    // RDKit✔️✔️:                        bool includeChiralPresence, bool includeRingStereo) {
    // RDKit✔️✔️: PRECONDITION(atomsInPlay.size() == mol.getNumAtoms(), "bad atomsInPlay size");
    // RDKit✔️✔️: PRECONDITION(bondsInPlay.size() == mol.getNumBonds(), "bad bondsInPlay size");
    // Symbol and mask validation is O(1) per table; the existing ranking
    // engine retains its source-shaped partitions and original-index output.
    if atoms_in_play.len() != topology.atoms.len() {
        return Err(CanonicalRankError::AtomMaskLength {
            expected: topology.atoms.len(),
            actual: atoms_in_play.len(),
        });
    }
    if bonds_in_play.len() != topology.bonds.len() {
        return Err(CanonicalRankError::BondMaskLength {
            expected: topology.bonds.len(),
            actual: bonds_in_play.len(),
        });
    }
    // RDKit✔️✔️: PRECONDITION(!atomSymbols || atomSymbols->size() == mol.getNumAtoms(),
    // RDKit✔️✔️:              "bad atomSymbols size");
    // RDKit✔️✔️: PRECONDITION(!bondSymbols || bondSymbols->size() == mol.getNumBonds(),
    // RDKit✔️✔️:              "bad bondSymbols size");
    if let Some(symbols) = atom_symbols
        && symbols.len() != topology.atoms.len()
    {
        return Err(CanonicalRankError::AtomSymbolLength {
            expected: topology.atoms.len(),
            actual: symbols.len(),
        });
    }
    if let Some(symbols) = bond_symbols
        && symbols.len() != topology.bonds.len()
    {
        return Err(CanonicalRankError::BondSymbolLength {
            expected: topology.bonds.len(),
            actual: symbols.len(),
        });
    }
    // RDKit✔️✔️: if (!mol.getNumAtoms()) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    if topology.atoms.is_empty() {
        return Ok(Vec::new());
    }
    // RDKit✔️✔️:   bool clearRings = false;
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:     clearRings = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   std::vector<Canon::canon_atom> atoms(mol.getNumAtoms());
    // RDKit✔️✔️:   detail::initFragmentCanonAtoms(mol, atoms, includeChirality, atomSymbols,
    // RDKit✔️✔️:                                  bondSymbols, atomsInPlay, bondsInPlay, true);
    let view = if let Some((valence, rings)) = prepared {
        CanonRankReadView::from_prepared_state(topology, valence, rings, Some(atoms_in_play))?
    } else {
        CanonRankReadView::from_topology(topology)?
    };
    let mut atoms = init_fragment_canon_atoms(
        &view,
        atoms_in_play,
        bonds_in_play,
        params.include_chirality,
        atom_symbols,
        bond_symbols,
    )?;
    // RDKit✔️✔️:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit✔️✔️:   ftor.df_useChirality = includeChirality;
    // RDKit✔️✔️:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit✔️✔️:   ftor.df_useChiralityRings = includeChirality;
    // RDKit✔️✔️:   ftor.df_useChiralPresence = includeChiralPresence;
    // RDKit✔️✔️:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit✔️✔️:                           includeRingStereo, &atomsInPlay, &bondsInPlay);
    let mut fragment_params = *params;
    fragment_params.include_stereo_groups = false;
    fragment_params.use_non_stereo_ranks = false;
    // The source fragment path always enables the chirality-ring comparison
    // when chirality is enabled; the whole-molecule path also gates on ring
    // stereo. These are independent policy branches within the same engine.
    fragment_params.chirality_rings_use_ring_stereo = false;
    rank_initialized_atoms(&view, &mut atoms, fragment_params)
    // RDKit✔️✔️:   if (clearRings) {
    // RDKit✔️✔️:     mol.getRingInfo()->reset();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: rankFragmentAtoms
}

#[derive(Debug, Clone, Copy)]
struct CanonRankFlags {
    use_isotopes: bool,
    use_chirality: bool,
    use_chirality_rings: bool,
    use_atom_maps: bool,
    use_non_stereo_ranks: bool,
    use_chiral_presence: bool,
    use_atom_maps_on_dummies: bool,
}

impl CanonRankFlags {
    const fn from_fragment_options(options: CanonicalRankParams) -> Self {
        Self {
            use_isotopes: options.include_isotopes,
            use_chirality: options.include_chirality,
            use_chirality_rings: if options.chirality_rings_use_ring_stereo {
                options.include_chirality && options.include_ring_stereo
            } else {
                options.include_chirality
            },
            use_atom_maps: options.include_atom_maps,
            use_non_stereo_ranks: options.use_non_stereo_ranks,
            use_chiral_presence: options.include_chiral_presence,
            use_atom_maps_on_dummies: true,
        }
    }
}

fn rank_initialized_atoms(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    options: CanonicalRankParams,
) -> Result<Vec<usize>, CanonicalRankError> {
    // RDKit✔️✔️:   AtomCompareFunctor ftor(&atoms.front(), mol, &atomsInPlay, &bondsInPlay);
    // RDKit✔️✔️:   ftor.df_useIsotopes = includeIsotopes;
    // RDKit✔️✔️:   ftor.df_useChirality = includeChirality;
    // RDKit✔️✔️:   ftor.df_useChiralityRings = includeChirality && includeRingStereo;
    // RDKit✔️✔️:   ftor.df_useAtomMaps = includeAtomMaps;
    // RDKit✔️✔️:   ftor.df_useNonStereoRanks = useNonStereoRanks;
    // RDKit✔️✔️:   ftor.df_useChiralPresence = includeChiralPresence;
    let mut order = vec![0usize; view.num_atoms()];
    let flags = CanonRankFlags::from_fragment_options(options);
    // RDKit✔️✔️:   detail::rankWithFunctor(ftor, breakTies, order, true, includeChirality,
    // RDKit✔️✔️:                           includeRingStereo, &atomsInPlay, &bondsInPlay);
    rank_with_atom_compare_functor_for_kekulize(
        view,
        atoms,
        options.break_ties,
        options.include_ring_stereo,
        flags,
        &mut order,
    )?;

    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     res[order[i]] = atoms[order[i]].index;
    // RDKit✔️✔️:   }
    let mut res = vec![0usize; view.num_atoms()];
    for idx in 0..view.num_atoms() {
        res[order[idx]] = usize::try_from(atoms[order[idx]].index).unwrap_or(usize::MAX);
    }
    Ok(res)
}

#[derive(Debug, Clone, Copy)]
struct CanonBondHolder<'a> {
    bond_type: BondOrder,
    bond_stereo: u8,
    stype: BondStereo,
    controlling_atoms: [Option<usize>; 4],
    nbr_sym_class: usize,
    nbr_idx: usize,
    p_symbol: Option<&'a str>,
    // Source-aligned storage for RDKit's needsInit=false symbol update branch.
    // rankFragmentAtoms forces needsInit=true today, so this remains dormant.
    #[allow(dead_code)]
    bond_idx: usize,
}

#[derive(Debug, Clone)]
struct CanonAtom<'a> {
    index: i32,
    is_in_play: bool,
    degree: usize,
    atomic_number: u8,
    isotope: u16,
    atom_map: u32,
    source_atom: &'a Atom,
    formal_charge: i8,
    chiral_tag: ChiralTag,
    total_num_hs: usize,
    is_ring_atom: bool,
    has_ring_nbr: bool,
    is_ring_stereo_atom: bool,
    which_stereo_group: usize,
    type_of_stereo_group: CanonStereoGroupType,
    neighbor_num: Vec<i32>,
    revisted_neighbors: Vec<i32>,
    all_nbr_ids: Vec<usize>,
    nbr_ids: Vec<usize>,
    p_symbol: Option<&'a str>,
    bonds: Vec<CanonBondHolder<'a>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum CanonStereoGroupType {
    Absolute = 0,
    Or = 1,
    And = 2,
}

#[derive(Debug, Clone, Copy)]
enum CanonCompareMode {
    Atom,
    SpecialChirality,
    SpecialSymmetry,
}

fn atropisomer_atoms_and_bonds_for_canonical_rank(
    view: &CanonRankReadView<'_>,
    bond_id: BondId,
) -> Option<[(AtomId, Vec<BondId>); 2]> {
    // BEGIN RDKIT CPP FUNCTION Atropisomers::getAtropisomerAtomsAndBonds
    // RDKit✔️✔️: bool getAtropisomerAtomsAndBonds(const Bond *bond,
    // RDKit✔️✔️:                                  AtropAtomAndBondVec atomsAndBondVects[2],
    // RDKit✔️✔️:                                  const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(bond, "no bond");
    let bond = view.bonds.get(bond_id.index())?;
    // RDKit✔️✔️:   atomsAndBondVects[0].first = bond->getBeginAtom();
    // RDKit✔️✔️:   atomsAndBondVects[1].first = bond->getEndAtom();
    let atoms = [bond.begin(), bond.end()];
    let mut result = [(atoms[0], Vec::new()), (atoms[1], Vec::new())];
    // RDKit✔️✔️:   for (int bondAtomIndex = 0; bondAtomIndex < 2; ++bondAtomIndex) {
    for bond_atom_index in 0..2 {
        // RDKit✔️✔️:     for (const auto nbrBond :
        // RDKit✔️✔️:          mol.atomBonds(atomsAndBondVects[bondAtomIndex].first)) {
        for neighbor_bond in view.bonds {
            if neighbor_bond.begin() != atoms[bond_atom_index]
                && neighbor_bond.end() != atoms[bond_atom_index]
            {
                continue;
            }
            // RDKit✔️✔️:       if (nbrBond == bond) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if neighbor_bond.id() == bond_id {
                continue;
            }
            // RDKit✔️✔️:       atomsAndBondVects[bondAtomIndex].second.push_back(nbrBond);
            result[bond_atom_index].1.push(neighbor_bond.id());
        }
        // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 0) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if result[bond_atom_index].1.is_empty() {
            return None;
        }
        // RDKit✔️✔️:     if (atomsAndBondVects[bondAtomIndex].second.size() == 2 &&
        // RDKit✔️✔️:         atomsAndBondVects[bondAtomIndex]
        // RDKit✔️✔️:                 .second[1]
        // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
        // RDKit✔️✔️:                 ->getIdx() <
        // RDKit✔️✔️:             atomsAndBondVects[bondAtomIndex]
        // RDKit✔️✔️:                 .second[0]
        // RDKit✔️✔️:                 ->getOtherAtom(atomsAndBondVects[bondAtomIndex].first)
        // RDKit✔️✔️:                 ->getIdx()) {
        // RDKit✔️✔️:       std::swap(atomsAndBondVects[bondAtomIndex].second[0],
        // RDKit✔️✔️:                 atomsAndBondVects[bondAtomIndex].second[1]);
        // RDKit✔️✔️:     }
        if result[bond_atom_index].1.len() == 2 {
            let first_other =
                bond_other_atom_index(view, result[bond_atom_index].1[0], atoms[bond_atom_index])?;
            let second_other =
                bond_other_atom_index(view, result[bond_atom_index].1[1], atoms[bond_atom_index])?;
            if second_other < first_other {
                result[bond_atom_index].1.swap(0, 1);
            }
        }
    }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION Atropisomers::getAtropisomerAtomsAndBonds
    Some(result)
}

fn bond_other_atom_index(
    view: &CanonRankReadView<'_>,
    bond_id: BondId,
    atom_id: AtomId,
) -> Option<usize> {
    view.bond_other_atom_index(bond_id, atom_id)
}

fn count_swaps_to_interconvert<T: Copy + Eq>(reference: &[T], mut probe: Vec<T>) -> usize {
    // BEGIN RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/utils.h :: countSwapsToInterconvert
    // RDKit✔️✔️: template <class T>
    // RDKit✔️✔️: unsigned int countSwapsToInterconvert(const T &ref, T probe) {
    // RDKit✔️✔️:   PRECONDITION(ref.size() == probe.size(), "size mismatch");
    // RDKit✔️✔️:   unsigned int nSwaps = 0;
    // RDKit✔️✔️:   while (refIt != ref.end()) {
    // RDKit✔️✔️:     if ((*probeIt) != (*refIt)) {
    // RDKit✔️✔️:       while ((*probeIt2) != (*refIt) && probeIt2 != probe.end()) {
    // RDKit✔️✔️:         ++probeIt2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       CHECK_INVARIANT(foundIt, "could not find probe element");
    // RDKit✔️✔️:       std::swap(*probeIt, *probeIt2);
    // RDKit✔️✔️:       nSwaps++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return nSwaps;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION third_party/rdkit/Code/RDGeneral/utils.h :: countSwapsToInterconvert
    assert_eq!(reference.len(), probe.len(), "size mismatch");
    let mut swaps = 0;
    for (index, expected) in reference.iter().enumerate() {
        if probe[index] == *expected {
            continue;
        }
        let found = probe[index..]
            .iter()
            .position(|candidate| candidate == expected)
            .map(|offset| index + offset)
            .expect("could not find probe element");
        probe.swap(index, found);
        swaps += 1;
    }
    swaps
}

fn empty_canon_atom_from_source_atom<'a>(atom: &'a Atom) -> CanonAtom<'a> {
    CanonAtom {
        index: i32::try_from(atom.id().index()).unwrap_or(i32::MAX),
        is_in_play: true,
        degree: 0,
        atomic_number: atom.atomic_number(),
        isotope: atom.isotope().unwrap_or(0),
        atom_map: atom.atom_map().unwrap_or(0),
        // Like source canon_atom::atom, borrow without reading rank metadata.
        // Property lookup and conversion belong to the guarded comparator.
        source_atom: atom,
        formal_charge: atom.formal_charge(),
        chiral_tag: atom.chiral_tag(),
        total_num_hs: 0,
        is_ring_atom: false,
        has_ring_nbr: false,
        is_ring_stereo_atom: false,
        which_stereo_group: 0,
        type_of_stereo_group: CanonStereoGroupType::Absolute,
        neighbor_num: Vec::new(),
        revisted_neighbors: Vec::new(),
        all_nbr_ids: Vec::new(),
        nbr_ids: Vec::new(),
        p_symbol: None,
        bonds: Vec::new(),
    }
}

fn init_canon_atoms<'a>(
    view: &CanonRankReadView<'a>,
    topology: &TopologyBlock,
    include_chirality: bool,
    include_stereo_groups: bool,
) -> Result<Vec<CanonAtom<'a>>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: initCanonAtoms
    // RDKit✔️✔️: void initCanonAtoms(const ROMol &mol, std::vector<Canon::canon_atom> &atoms,
    // RDKit✔️✔️:                     bool includeChirality, bool includeStereoGroups) {
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     basicInitCanonAtom(mol, atoms[i], i);
    // RDKit✔️✔️:     advancedInitCanonAtom(mol, atoms[i], i);
    // RDKit✔️✔️:     atoms[i].bonds.reserve(atoms[i].degree);
    // RDKit✔️✔️:     getBonds(mol, atoms[i].atom, atoms[i].bonds, includeChirality, atoms);
    // RDKit✔️✔️:   }
    let atoms_in_play = vec![true; view.num_atoms()];
    let bonds_in_play = vec![true; view.bonds.len()];
    let mut atoms = init_fragment_canon_atoms(
        view,
        &atoms_in_play,
        &bonds_in_play,
        include_chirality,
        None,
        None,
    )?;
    // RDKit✔️✔️:   if (includeChirality && includeStereoGroups) {
    // RDKit✔️✔️:     unsigned int sgidx = 1;
    // RDKit✔️✔️:     for (auto &sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:       for (auto atom : sg.getAtoms()) {
    // RDKit✔️✔️:         atoms[atom->getIdx()].whichStereoGroup = sgidx;
    // RDKit✔️✔️:         atoms[atom->getIdx()].typeOfStereoGroup = sg.getGroupType();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++sgidx;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if include_chirality && include_stereo_groups {
        for (group_index, group) in topology.stereo_groups.iter().enumerate() {
            let kind = match group.kind() {
                StereoGroupKind::Absolute => CanonStereoGroupType::Absolute,
                StereoGroupKind::Or => CanonStereoGroupType::Or,
                StereoGroupKind::And => CanonStereoGroupType::And,
            };
            for atom in group.atoms() {
                atoms[atom.index()].which_stereo_group = group_index + 1;
                atoms[atom.index()].type_of_stereo_group = kind;
            }
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: initCanonAtoms
    Ok(atoms)
}

fn init_fragment_canon_atoms<'a>(
    view: &CanonRankReadView<'a>,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    include_chirality: bool,
    atom_symbols: Option<&'a [String]>,
    bond_symbols: Option<&'a [String]>,
) -> Result<Vec<CanonAtom<'a>>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: detail::initFragmentCanonAtoms
    // RDKit✔️✔️: void initFragmentCanonAtoms(const ROMol &mol,
    // RDKit✔️✔️:                             std::vector<Canon::canon_atom> &atoms,
    // RDKit✔️✔️:                             bool includeChirality,
    // RDKit✔️✔️:                             const std::vector<std::string> *atomSymbols,
    // RDKit✔️✔️:                             const std::vector<std::string> *bondSymbols,
    // RDKit✔️✔️:                             const boost::dynamic_bitset<> &atomsInPlay,
    // RDKit✔️✔️:                             const boost::dynamic_bitset<> &bondsInPlay,
    // RDKit✔️✔️:                             bool needsInit) {
    // RDKit✔️✔️:   needsInit = true;
    // Symbol references borrow caller storage without per-atom or per-bond
    // string allocation. Each included bond assigns both holder references
    // in its existing O(B) initialization pass.
    let mut atoms = view
        .atoms
        .iter()
        .map(empty_canon_atom_from_source_atom)
        .collect::<Vec<_>>();
    // RDKit✔️✔️:   // start by initializing the atoms
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     auto i = atom->getIdx();
    // RDKit✔️✔️:     auto &atomsi = atoms[i];
    // RDKit✔️✔️:     atomsi.atom = atom;
    // RDKit✔️✔️:     atomsi.index = i;
    // RDKit✔️✔️:     atomsi.degree = 0;
    // RDKit✔️✔️:     if (atomsInPlay[i]) {
    // RDKit✔️✔️:       if (atomSymbols) {
    // RDKit✔️✔️:         atomsi.p_symbol = &(*atomSymbols)[i];
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         atomsi.p_symbol = nullptr;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (needsInit) {
    // RDKit✔️✔️:         atomsi.nbrIds = std::make_unique<int[]>(atom->getDegree());
    // RDKit✔️✔️:         advancedInitCanonAtom(mol, atomsi, i);
    // RDKit✔️✔️:         atomsi.bonds.reserve(4);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for atom_idx in 0..view.num_atoms() {
        atoms[atom_idx].is_in_play = atoms_in_play[atom_idx];
        atoms[atom_idx].degree = 0;
        atoms[atom_idx].all_nbr_ids.clear();
        atoms[atom_idx].nbr_ids.clear();
        atoms[atom_idx].bonds.clear();
        if !atoms_in_play[atom_idx] {
            continue;
        }
        atoms[atom_idx].p_symbol = atom_symbols.map(|symbols| symbols[atom_idx].as_str());
        let atom = &view.atoms[atom_idx];
        // RDKit✔️✔️: void advancedInitCanonAtom(const ROMol &mol, Canon::canon_atom &atom,
        // RDKit✔️✔️:                            const int &) {
        // RDKit✔️✔️:   atom.totalNumHs = atom.atom->getTotalNumHs();
        // RDKit✔️✔️:   atom.isRingStereoAtom =
        // RDKit✔️✔️:       (atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
        // RDKit✔️✔️:        atom.atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
        // RDKit✔️✔️:       atom.atom->hasProp(common_properties::_ringStereoAtoms);
        // RDKit✔️✔️:   atom.hasRingNbr = hasRingNbr(mol, atom.atom);
        // RDKit✔️✔️: }
        // The same implicit getter applies to selected prepared fragment rows;
        // NoImplicit and signed storage are source reads, never max(0) defaults.
        // Explicit uint8 + nonnegative int8 is bounded by 382; O(1), no clone.
        atoms[atom_idx].total_num_hs = usize::from(atom.explicit_hydrogens())
            + crate::hcount::implicit_hydrogen_count(atom, &view.valence).map_err(|_| {
                CanonicalRankError::PreparedValenceInvalid {
                    atom_index: atom_idx,
                }
            })? as usize;
        atoms[atom_idx].is_ring_atom = view.rings.num_atom_rings(atom.id()) > 0;
        atoms[atom_idx].is_ring_stereo_atom = matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) && atom.prop("_ringStereoAtoms").is_some();
        atoms[atom_idx].has_ring_nbr = has_ring_nbr_for_kekulize(view, atom_idx);
        atoms[atom_idx].bonds.reserve(4);
    }

    // RDKit✔️✔️:   // now deal with the bonds in the fragment.
    // RDKit✔️✔️:   if (needsInit) {
    // RDKit✔️✔️:     for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:       if (!bondsInPlay[bond->getIdx()] ||
    // RDKit✔️✔️:           !atomsInPlay[bond->getBeginAtomIdx()] ||
    // RDKit✔️✔️:           !atomsInPlay[bond->getEndAtomIdx()]) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Canon::canon_atom &begAt = atoms[bond->getBeginAtomIdx()];
    // RDKit✔️✔️:       Canon::canon_atom &endAt = atoms[bond->getEndAtomIdx()];
    // RDKit✔️✔️:       begAt.nbrIds[begAt.degree++] = bond->getEndAtomIdx();
    // RDKit✔️✔️:       endAt.nbrIds[endAt.degree++] = bond->getBeginAtomIdx();
    // RDKit✔️✔️:       begAt.bonds.push_back(
    // RDKit✔️✔️:           makeBondHolder(bond, bond->getEndAtomIdx(), includeChirality, atoms));
    // RDKit✔️✔️:       endAt.bonds.push_back(makeBondHolder(bond, bond->getBeginAtomIdx(),
    // RDKit✔️✔️:                                            includeChirality, atoms));
    // RDKit✔️✔️:       if (bondSymbols) {
    // RDKit✔️✔️:         begAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
    // RDKit✔️✔️:         endAt.bonds.back().p_symbol = &(*bondSymbols)[bond->getIdx()];
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for (bond_idx, bond) in view.bonds.iter().enumerate() {
        let begin = bond.begin().index();
        let end = bond.end().index();
        if !bonds_in_play[bond_idx] || !atoms_in_play[begin] || !atoms_in_play[end] {
            continue;
        }
        atoms[begin].degree += 1;
        atoms[begin].all_nbr_ids.push(end);
        atoms[begin].nbr_ids.push(end);
        atoms[end].degree += 1;
        atoms[end].all_nbr_ids.push(begin);
        atoms[end].nbr_ids.push(begin);
        let mut begin_holder = make_canon_bond_holder(view, bond_idx, end, include_chirality)?;
        let mut end_holder = make_canon_bond_holder(view, bond_idx, begin, include_chirality)?;
        if let Some(symbols) = bond_symbols {
            begin_holder.p_symbol = Some(symbols[bond_idx].as_str());
            end_holder.p_symbol = Some(symbols[bond_idx].as_str());
        }
        atoms[begin].bonds.push(begin_holder);
        atoms[end].bonds.push(end_holder);
    }

    // RDKit✔️✔️:   for (size_t i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     if (!atomsInPlay[i]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto &atomsi = atoms[i];
    // RDKit✔️✔️:     if (needsInit) {
    // RDKit✔️✔️:       atomsi.totalNumHs += (mol.getAtomWithIdx(i)->getDegree() - atomsi.degree);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::sort(atomsi.bonds.begin(), atomsi.bonds.end(), bondholder::greater);
    // RDKit✔️✔️:   }
    for atom_idx in 0..view.num_atoms() {
        if !atoms_in_play[atom_idx] {
            continue;
        }
        atoms[atom_idx].total_num_hs +=
            view.atom_degree(AtomId::new(atom_idx)) - atoms[atom_idx].degree;
        let initial_ranks = canon_atom_rank_snapshot(&atoms);
        sort_canon_bonds_descending(&mut atoms[atom_idx].bonds, &initial_ranks);
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/new_canon.cpp :: detail::initFragmentCanonAtoms
    Ok(atoms)
}

fn has_ring_nbr_for_kekulize(view: &CanonRankReadView<'_>, atom_idx: usize) -> bool {
    // BEGIN RDKIT CPP FUNCTION hasRingNbr
    // RDKit✔️✔️: bool hasRingNbr(const ROMol &mol, const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad pointer");
    // RDKit✔️✔️:   for (const auto nbr : mol.atomNeighbors(at)) {
    // RDKit✔️✔️:     if ((nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:          nbr->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) &&
    // RDKit✔️✔️:         nbr->hasProp(common_properties::_ringStereoAtoms)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hasRingNbr
    view.adjacency
        .neighbors_of(atom_idx)
        .iter()
        .any(|neighbor_ref| {
            let neighbor = &view.atoms[neighbor_ref.atom_index];
            matches!(
                neighbor.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) && neighbor.prop("_ringStereoAtoms").is_some()
        })
}

fn make_canon_bond_holder<'a>(
    view: &CanonRankReadView<'_>,
    bond_idx: usize,
    other_idx: usize,
    include_chirality: bool,
) -> Result<CanonBondHolder<'a>, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION makeBondHolder
    // RDKit✔️✔️: bondholder makeBondHolder(const Bond *bond, unsigned int otherIdx,
    // RDKit✔️✔️:                           bool includeChirality,
    // RDKit✔️✔️:                           const std::vector<Canon::canon_atom> &atoms) {
    // RDKit✔️✔️:   PRECONDITION(bond, "bad pointer");
    // RDKit✔️✔️:   Bond::BondStereo stereo = Bond::STEREONONE;
    // RDKit✔️✔️:   if (includeChirality) {
    // RDKit✔️✔️:     stereo = bond->getStereo();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Bond::BondType bt =
    // RDKit✔️✔️:       bond->getIsAromatic() ? Bond::AROMATIC : bond->getBondType();
    // RDKit✔️✔️:   bondholder res(bt, stereo, otherIdx, 0, bond->getIdx());
    let bond = &view.bonds[bond_idx];
    let bond_type = if bond.is_aromatic() {
        BondOrder::Aromatic
    } else {
        bond.order()
    };
    let mut controlling_atoms = [None; 4];
    // RDKit✔️✔️:   if (includeChirality) {
    // RDKit✔️✔️:     res.stype = bond->getStereo();
    // RDKit✔️✔️:     if (res.stype == Bond::BondStereo::STEREOCIS ||
    // RDKit✔️✔️:         res.stype == Bond::BondStereo::STEREOTRANS) {
    // RDKit✔️✔️:       res.controllingAtoms[0] = &atoms[bond->getStereoAtoms()[0]];
    // RDKit✔️✔️:       res.controllingAtoms[2] = &atoms[bond->getStereoAtoms()[1]];
    if include_chirality && matches!(bond.stereo(), BondStereo::Cis | BondStereo::Trans) {
        let Some(stereo_atoms) = bond.stereo_atoms() else {
            return Err(CanonicalRankError::ProtocolDebt {
                branch: "makeBondHolder cis/trans stereo atoms",
                reason: "cis/trans bond stereo requires the two RDKit stereo atom references",
            });
        };
        controlling_atoms[0] = Some(stereo_atoms[0].index());
        controlling_atoms[2] = Some(stereo_atoms[1].index());
        // RDKit✔️✔️:       if (bond->getBeginAtom()->getDegree() > 2) {
        // RDKit✔️✔️:         for (const auto nbr :
        // RDKit✔️✔️:              bond->getOwningMol().atomNeighbors(bond->getBeginAtom())) {
        // RDKit✔️✔️:           if (nbr->getIdx() != bond->getEndAtomIdx() &&
        // RDKit✔️✔️:               nbr->getIdx() !=
        // RDKit✔️✔️:                   static_cast<unsigned int>(bond->getStereoAtoms()[0])) {
        // RDKit✔️✔️:             res.controllingAtoms[1] = &atoms[nbr->getIdx()];
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        if view.atom_degree(bond.begin()) > 2 {
            for neighbor in view.atom_neighbors(bond.begin()) {
                if neighbor != bond.end().index() && neighbor != stereo_atoms[0].index() {
                    controlling_atoms[1] = Some(neighbor);
                }
            }
        }
        // RDKit✔️✔️:       if (bond->getEndAtom()->getDegree() > 2) {
        // RDKit✔️✔️:         for (const auto nbr :
        // RDKit✔️✔️:              bond->getOwningMol().atomNeighbors(bond->getEndAtom())) {
        // RDKit✔️✔️:           if (nbr->getIdx() != bond->getBeginAtomIdx() &&
        // RDKit✔️✔️:               nbr->getIdx() !=
        // RDKit✔️✔️:                   static_cast<unsigned int>(bond->getStereoAtoms()[1])) {
        // RDKit✔️✔️:             res.controllingAtoms[3] = &atoms[nbr->getIdx()];
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        if view.atom_degree(bond.end()) > 2 {
            for neighbor in view.atom_neighbors(bond.end()) {
                if neighbor != bond.begin().index() && neighbor != stereo_atoms[1].index() {
                    controlling_atoms[3] = Some(neighbor);
                }
            }
        }
    }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (res.stype == Bond::BondStereo::STEREOATROPCCW ||
    // RDKit✔️✔️:         res.stype == Bond::BondStereo::STEREOATROPCW) {
    if include_chirality && matches!(bond.stereo(), BondStereo::AtropCcw | BondStereo::AtropCw) {
        // RDKit✔️✔️:       Atropisomers::AtropAtomAndBondVec atropAtomAndBondVecs[2];
        // RDKit✔️✔️:       CHECK_INVARIANT(Atropisomers::getAtropisomerAtomsAndBonds(
        // RDKit✔️✔️:                           bond, atropAtomAndBondVecs, bond->getOwningMol()),
        // RDKit✔️✔️:                       "Could not find atropisomer controlling atoms")
        let Some(atrop) = atropisomer_atoms_and_bonds_for_canonical_rank(view, bond.id()) else {
            return Err(CanonicalRankError::ProtocolDebt {
                branch: "Atropisomers::getAtropisomerAtomsAndBonds invariant",
                reason: "could not find atropisomer controlling atoms",
            });
        };
        // RDKit✔️✔️:       res.controllingAtoms[0] =
        // RDKit✔️✔️:           &atoms[atropAtomAndBondVecs[0]
        // RDKit✔️✔️:                      .second[0]
        // RDKit✔️✔️:                      ->getOtherAtom(atropAtomAndBondVecs[0].first)
        // RDKit✔️✔️:                      ->getIdx()];
        controlling_atoms[0] = Some(
            bond_other_atom_index(view, atrop[0].1[0], atrop[0].0)
                .expect("atropisomer neighbor bond must include focus atom"),
        );
        // RDKit✔️✔️:       res.controllingAtoms[2] =
        // RDKit✔️✔️:           &atoms[atropAtomAndBondVecs[1]
        // RDKit✔️✔️:                      .second[0]
        // RDKit✔️✔️:                      ->getOtherAtom(atropAtomAndBondVecs[1].first)
        // RDKit✔️✔️:                      ->getIdx()];
        controlling_atoms[2] = Some(
            bond_other_atom_index(view, atrop[1].1[0], atrop[1].0)
                .expect("atropisomer neighbor bond must include focus atom"),
        );
        // RDKit✔️✔️:       if (atropAtomAndBondVecs[0].second.size() > 1) {
        // RDKit✔️✔️:         res.controllingAtoms[1] =
        // RDKit✔️✔️:             &atoms[atropAtomAndBondVecs[0]
        // RDKit✔️✔️:                        .second[1]
        // RDKit✔️✔️:                        ->getOtherAtom(atropAtomAndBondVecs[0].first)
        // RDKit✔️✔️:                        ->getIdx()];
        // RDKit✔️✔️:       }
        if atrop[0].1.len() > 1 {
            controlling_atoms[1] = Some(
                bond_other_atom_index(view, atrop[0].1[1], atrop[0].0)
                    .expect("atropisomer neighbor bond must include focus atom"),
            );
        }
        // RDKit✔️✔️:       if (atropAtomAndBondVecs[1].second.size() > 1) {
        // RDKit✔️✔️:         res.controllingAtoms[3] =
        // RDKit✔️✔️:             &atoms[atropAtomAndBondVecs[1]
        // RDKit✔️✔️:                        .second[1]
        // RDKit✔️✔️:                        ->getOtherAtom(atropAtomAndBondVecs[1].first)
        // RDKit✔️✔️:                        ->getIdx()];
        // RDKit✔️✔️:     }
        if atrop[1].1.len() > 1 {
            controlling_atoms[3] = Some(
                bond_other_atom_index(view, atrop[1].1[1], atrop[1].0)
                    .expect("atropisomer neighbor bond must include focus atom"),
            );
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION makeBondHolder
    let stype = if include_chirality {
        bond.stereo()
    } else {
        BondStereo::None
    };
    Ok(CanonBondHolder {
        bond_type,
        bond_stereo: rdkit_bond_stereo_rank(stype),
        stype,
        controlling_atoms,
        nbr_sym_class: 0,
        nbr_idx: other_idx,
        p_symbol: None,
        bond_idx,
    })
}

fn new_hanoi_scratch_for_canonical_rank(len: usize, top_level: bool) -> Vec<usize> {
    // RDKit❗❌:   std::vector<int> hanoiTemp(nAts);
    // RDKit❗❌:     localHanoiTemp.resize(nAtoms);
    // Exact constructor scope and O(n) zero-fill. Existing usize indices use
    // 8-byte slots on64-bit targets versus C++ int4; preserve that baseline
    // representation and qualify its cache/space cost, not a new width rewrite.
    let scratch = vec![0; len];
    #[cfg(test)]
    chem31_scratch_trace::record(chem31_scratch_trace::Event::Construct(
        top_level,
        len,
        scratch.as_ptr() as usize,
    ));
    #[cfg(not(test))]
    let _ = top_level;
    scratch
}

fn rank_with_atom_compare_functor_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    break_ties: bool,
    include_ring_stereo: bool,
    flags: CanonRankFlags,
    order: &mut [usize],
) -> Result<(), CanonicalRankError> {
    // BEGIN COMPLETE RDKit .6 CHEM31 rankWithFunctor
    // RDKit❗❌: template <typename T>
    // RDKit❗❌: void rankWithFunctor(T &ftor, bool breakTies, std::vector<int> &order,
    // RDKit❗❌:                      bool useSpecial, bool useChirality, bool includeRingStereo,
    // RDKit❗❌:                      const boost::dynamic_bitset<> *atomsInPlay,
    // RDKit❗❌:                      const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit❗❌:   PRECONDITION(!order.empty(), "order should not be empty");
    // RDKit❗❌:   const ROMol &mol = *ftor.dp_mol;
    // RDKit❗❌:   canon_atom *atoms = ftor.dp_atoms;
    // RDKit❗❌:   const unsigned int nAts = mol.getNumAtoms();
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<int> count(nAts);
    // RDKit❗❌:   std::vector<int> next(nAts);
    // RDKit❗❌:   std::vector<int> changed(nAts, 1);
    // RDKit❗❌:   std::vector<char> touched(nAts, 0);
    // RDKit❗❌:   std::vector<int> hanoiTemp(nAts);
    // RDKit❗❌:   int activeset;
    // RDKit❗❌:   CreateSinglePartition(nAts, order, count, atoms);
    // RDKit❗❌: // ActivatePartitions(nAts,order,count,activeset,next,changed);
    // RDKit❗❌: // RefinePartitions(mol,atoms,ftor,false,order,count,activeset,next,changed,touched);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:   std::cerr << "1--------" << std::endl;
    // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:               << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌: #endif
    // RDKit❗❌:   ftor.df_useNbrs = true;
    // RDKit❗❌:   ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:   std::cerr << "1a--------" << std::endl;
    // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:               << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌: #endif
    // RDKit❗❌:   RefinePartitions(mol, atoms, ftor, true, order, count, activeset, next,
    // RDKit❗❌:                    changed, touched, &hanoiTemp);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:   std::cerr << "2--------" << std::endl;
    // RDKit❗❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:     std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:               << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:   }
    // RDKit❗❌: #endif
    // RDKit❗❌:   bool ties = false;
    // RDKit❗❌:   for (unsigned i = 0; i < nAts; ++i) {
    // RDKit❗❌:     if (!count[i]) {
    // RDKit❗❌:       ties = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (useChirality && ties && includeRingStereo) {
    // RDKit❗❌:     SpecialChiralityAtomCompareFunctor scftor(atoms, mol, atomsInPlay,
    // RDKit❗❌:                                               bondsInPlay);
    // RDKit❗❌:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit❗❌:     RefinePartitions(mol, atoms, scftor, true, order, count, activeset, next,
    // RDKit❗❌:                      changed, touched, &hanoiTemp);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:     std::cerr << "2a--------" << std::endl;
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:       std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:                 << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // RDKit❗❌:   }
    // RDKit❗❌:   ties = false;
    // RDKit❗❌:   unsigned symRingAtoms = 0;
    // RDKit❗❌:   unsigned ringAtoms = 0;
    // RDKit❗❌:   bool branchingRingAtom = false;
    // RDKit❗❌:   RingInfo *ringInfo = mol.getRingInfo();
    // RDKit❗❌:   for (unsigned i = 0; i < nAts; ++i) {
    // RDKit❗❌:     if (ringInfo->isInitialized() && ringInfo->numAtomRings(order[i])) {
    // RDKit❗❌:       if (count[order[i]] > 2) {
    // RDKit❗❌:         symRingAtoms += count[order[i]];
    // RDKit❗❌:       }
    // RDKit❗❌:       ringAtoms++;
    // RDKit❗❌:       if (ringInfo->isInitialized() && ringInfo->numAtomRings(order[i]) > 1 &&
    // RDKit❗❌:           count[order[i]] > 1) {
    // RDKit❗❌:         branchingRingAtom = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!count[i]) {
    // RDKit❗❌:       ties = true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   //      std::cout << " " << ringAtoms << " "  << symRingAtoms << std::endl;
    // RDKit❗❌:   if (useSpecial && ties && ringAtoms > 0 &&
    // RDKit❗❌:       static_cast<float>(symRingAtoms) / ringAtoms > 0.5 && branchingRingAtom) {
    // RDKit❗❌:     SpecialSymmetryAtomCompareFunctor sftor(atoms, mol, atomsInPlay,
    // RDKit❗❌:                                             bondsInPlay);
    // RDKit❗❌:     compareRingAtomsConcerningNumNeighbors(atoms, nAts, mol);
    // RDKit❗❌:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit❗❌:     RefinePartitions(mol, atoms, sftor, true, order, count, activeset, next,
    // RDKit❗❌:                      changed, touched, &hanoiTemp);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:     std::cerr << "2b--------" << std::endl;
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:       std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:                 << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // RDKit❗❌:   }
    // RDKit❗❌:   if (breakTies) {
    // RDKit❗❌:     BreakTies(mol, atoms, ftor, true, order, count, activeset, next, changed,
    // RDKit❗❌:               touched, &hanoiTemp);
    // RDKit❗❌: #ifdef VERBOSE_CANON
    // RDKit❗❌:     std::cerr << "3--------" << std::endl;
    // RDKit❗❌:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗❌:       std::cerr << order[i] + 1 << " " << " index: " << atoms[order[i]].index
    // RDKit❗❌:                 << " count: " << count[order[i]] << std::endl;
    // RDKit❗❌:     }
    // RDKit❗❌: #endif
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END COMPLETE RDKit .6 CHEM31 rankWithFunctor
    // Existing comparer/helper allocations and useSpecial/typed-state baseline
    // remain qualified; this delta shares one scratch across all reached modes.

    let n_atoms = view.num_atoms();
    let mut count = vec![0usize; n_atoms];
    let mut next = vec![-2isize; n_atoms];
    let mut changed = vec![true; n_atoms];
    let mut touched = vec![false; n_atoms];
    let mut hanoi_temp = new_hanoi_scratch_for_canonical_rank(n_atoms, true);
    let mut active_set = -1isize;
    create_single_partition_for_kekulize(n_atoms, order, &mut count, atoms);
    activate_partitions_for_kekulize(
        n_atoms,
        order,
        &count,
        &mut active_set,
        &mut next,
        &mut changed,
    );
    refine_partitions_for_kekulize(
        view,
        atoms,
        true,
        CanonCompareMode::Atom,
        flags,
        order,
        &mut count,
        &mut active_set,
        &mut next,
        &mut changed,
        &mut touched,
        Some(&mut hanoi_temp),
    )?;
    let ties = count.iter().any(|&value| value == 0);
    if flags.use_chirality && ties && include_ring_stereo {
        activate_partitions_for_kekulize(
            n_atoms,
            order,
            &count,
            &mut active_set,
            &mut next,
            &mut changed,
        );
        refine_partitions_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::SpecialChirality,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
            Some(&mut hanoi_temp),
        )?;
    }
    let use_special_symmetry = special_symmetry_rank_refinement_required(view, order, &count);
    if use_special_symmetry {
        compare_ring_atoms_concerning_num_neighbors_for_kekulize(view, atoms);
        activate_partitions_for_kekulize(
            n_atoms,
            order,
            &count,
            &mut active_set,
            &mut next,
            &mut changed,
        );
        refine_partitions_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::SpecialSymmetry,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
            Some(&mut hanoi_temp),
        )?;
    }
    if break_ties {
        break_ties_for_kekulize(
            view,
            atoms,
            true,
            CanonCompareMode::Atom,
            flags,
            order,
            &mut count,
            &mut active_set,
            &mut next,
            &mut changed,
            &mut touched,
            Some(&mut hanoi_temp),
        )?;
    }
    Ok(())
}

fn create_single_partition_for_kekulize(
    n_atoms: usize,
    order: &mut [usize],
    count: &mut [usize],
    atoms: &mut [CanonAtom<'_>],
) {
    // BEGIN RDKIT CPP FUNCTION CreateSinglePartition
    // RDKit✔️✔️: void CreateSinglePartition(unsigned int nAtoms, std::vector<int> &order,
    // RDKit✔️✔️:                            std::vector<int> &count, canon_atom *atoms) {
    // RDKit✔️✔️:   PRECONDITION(!order.empty(), "order should not be empty");
    // RDKit✔️✔️:   PRECONDITION(!count.empty(), "count should not be empty");
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     atoms[i].index = 0;
    // RDKit✔️✔️:     order[i] = i;
    // RDKit✔️✔️:     count[i] = 0;
    // RDKit✔️✔️:   }
    for i in 0..n_atoms {
        atoms[i].index = 0;
        order[i] = i;
        count[i] = 0;
    }
    // RDKit✔️✔️:   count[0] = nAtoms;
    // RDKit✔️✔️: }
    count[0] = n_atoms;
    // END RDKIT CPP FUNCTION CreateSinglePartition
}

fn activate_partitions_for_kekulize(
    n_atoms: usize,
    order: &[usize],
    count: &[usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
) {
    // BEGIN RDKIT CPP FUNCTION ActivatePartitions
    // RDKit✔️✔️: void ActivatePartitions(unsigned int nAtoms, std::vector<int> &order,
    // RDKit✔️✔️:                         std::vector<int> &count, int &activeset,
    // RDKit✔️✔️:                         std::vector<int> &next, std::vector<int> &changed) {
    // RDKit✔️✔️:   unsigned int i, j;
    // RDKit✔️✔️:   activeset = -1;
    *active_set = -1;
    // RDKit✔️✔️:   for (i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     next[i] = -2;
    // RDKit✔️✔️:   }
    next.fill(-2);
    // RDKit✔️✔️:   i = 0;
    // RDKit✔️✔️:   do {
    // RDKit✔️✔️:     j = order[i];
    // RDKit✔️✔️:     if (count[j] > 1) {
    // RDKit✔️✔️:       next[j] = activeset;
    // RDKit✔️✔️:       activeset = j;
    // RDKit✔️✔️:       i += count[j];
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       i++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } while (i < nAtoms);
    let mut i = 0usize;
    while i < n_atoms {
        let j = order[i];
        if count[j] > 1 {
            next[j] = *active_set;
            *active_set = isize::try_from(j).unwrap_or(isize::MAX);
            i += count[j];
        } else {
            i += 1;
        }
    }
    // RDKit✔️✔️:   for (i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:     j = order[i];
    // RDKit✔️✔️:     int flag = 1;
    // RDKit✔️✔️:     changed[j] = flag;
    // RDKit✔️✔️:   }
    for &j in order.iter().take(n_atoms) {
        changed[j] = true;
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION ActivatePartitions
}

#[allow(clippy::too_many_arguments)]
fn refine_partitions_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    mode: bool,
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
    order: &mut [usize],
    count: &mut [usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched_partitions: &mut [bool],
    hanoi_temp: Option<&mut [usize]>,
) -> Result<(), CanonicalRankError> {
    // BEGIN COMPLETE RDKit .6 CHEM31 RefinePartitions
    // RDKit❗❌: template <typename CompareFunc>
    // RDKit❗❌: void RefinePartitions(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit❗❌:                       int mode, std::vector<int> &order,
    // RDKit❗❌:                       std::vector<int> &count, int &activeset,
    // RDKit❗❌:                       std::vector<int> &next, std::vector<int> &changed,
    // RDKit❗❌:                       std::vector<char> &touchedPartitions,
    // RDKit❗❌:                       std::vector<int> *hanoiTemp = nullptr) {
    // RDKit❗❌:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:   std::vector<int> localHanoiTemp;
    // RDKit❗❌:   if (!hanoiTemp) {
    // RDKit❗❌:     localHanoiTemp.resize(nAtoms);
    // RDKit❗❌:     hanoiTemp = &localHanoiTemp;
    // RDKit❗❌:   }
    // RDKit❗❌:   int partition;
    // RDKit❗❌:   int symclass = 0;
    // RDKit❗❌:   int offset;
    // RDKit❗❌:   int index;
    // RDKit❗❌:   int len;
    // RDKit❗❌:   int i;
    // RDKit❗❌:   PRECONDITION(hanoiTemp->size() >= nAtoms, "hanoi scratch is too small");
    // RDKit❗❌:   // std::vector<char> touchedPartitions(mol.getNumAtoms(),0);
    // RDKit❗❌:
    // RDKit❗❌:   // std::cerr<<"&&&&&&&&&&&&&&&& RP"<<std::endl;
    // RDKit❗❌:   while (activeset != -1) {
    // RDKit❗❌:     // std::cerr<<"ITER: "<<activeset<<" next: "<<next[activeset]<<std::endl;
    // RDKit❗❌:     // std::cerr<<" next: ";
    // RDKit❗❌:     // for(unsigned int ii=0;ii<nAtoms;++ii){
    // RDKit❗❌:     //   std::cerr<<ii<<":"<<next[ii]<<" ";
    // RDKit❗❌:     // }
    // RDKit❗❌:     // std::cerr<<std::endl;
    // RDKit❗❌:     // for(unsigned int ii=0;ii<nAtoms;++ii){
    // RDKit❗❌:     //   std::cerr<<order[ii]<<" count: "<<count[order[ii]]<<" index:
    // RDKit❗❌:     //   "<<atoms[order[ii]].index<<std::endl;
    // RDKit❗❌:     // }
    // RDKit❗❌:
    // RDKit❗❌:     partition = activeset;
    // RDKit❗❌:     activeset = next[partition];
    // RDKit❗❌:     next[partition] = -2;
    // RDKit❗❌:
    // RDKit❗❌:     len = count[partition];
    // RDKit❗❌:     offset = atoms[partition].index;
    // RDKit❗❌:     auto start = std::span<int>(&order[offset], len);
    // RDKit❗❌:     // std::cerr<<"\n\n**************************************************************"<<std::endl;
    // RDKit❗❌:     // std::cerr<<"  sort - class:"<<atoms[partition].index<<" len:
    // RDKit❗❌:     // "<<len<<":"; for(unsigned int ii=0;ii<len;++ii){
    // RDKit❗❌:     //   std::cerr<<" "<<order[offset+ii]+1;
    // RDKit❗❌:     // }
    // RDKit❗❌:     // std::cerr<<std::endl;
    // RDKit❗❌:     // for(unsigned int ii=0;ii<nAtoms;++ii){
    // RDKit❗❌:     //   std::cerr<<order[ii]+1<<" count: "<<count[order[ii]]<<" index:
    // RDKit❗❌:     //   "<<atoms[order[ii]].index<<std::endl;
    // RDKit❗❌:     // }
    // RDKit❗❌:     if (RDKit::detail::hanoi(start.data(), len, hanoiTemp->data(), count.data(),
    // RDKit❗❌:                              changed.data(), compar)) {
    // RDKit❗❌:       std::copy_n(hanoiTemp->begin(), len, start.begin());
    // RDKit❗❌:     }
    // RDKit❗❌:     // std::cerr<<"*_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*"<<std::endl;
    // RDKit❗❌:     // std::cerr<<"  result:";
    // RDKit❗❌:     // for(unsigned int ii=0;ii<nAtoms;++ii){
    // RDKit❗❌:     //    std::cerr<<order[ii]+1<<" count: "<<count[order[ii]]<<" index:
    // RDKit❗❌:     //    "<<atoms[order[ii]].index<<std::endl;
    // RDKit❗❌:     //  }
    // RDKit❗❌:     for (int k = 0; k < len; ++k) {
    // RDKit❗❌:       changed[start[k]] = 0;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     index = start[0];
    // RDKit❗❌:     // std::cerr<<"  len:"<<len<<" index:"<<index<<"
    // RDKit❗❌:     // count:"<<count[index]<<std::endl;
    // RDKit❗❌:     for (i = count[index]; i < len; i++) {
    // RDKit❗❌:       index = start[i];
    // RDKit❗❌:       if (count[index]) {
    // RDKit❗❌:         symclass = offset + i;
    // RDKit❗❌:       }
    // RDKit❗❌:       atoms[index].index = symclass;
    // RDKit❗❌:       // std::cerr<<" "<<index+1<<"("<<symclass<<")";
    // RDKit❗❌:       // if(mode && (activeset<0 || count[index]>count[activeset]) ){
    // RDKit❗❌:       //  activeset=index;
    // RDKit❗❌:       //}
    // RDKit❗❌:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
    // RDKit❗❌:         changed[atoms[index].nbrIds[j]] = 1;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // std::cerr<<std::endl;
    // RDKit❗❌:
    // RDKit❗❌:     if (mode) {
    // RDKit❗❌:       index = start[0];
    // RDKit❗❌:       for (i = count[index]; i < len; i++) {
    // RDKit❗❌:         index = start[i];
    // RDKit❗❌:         for (unsigned j = 0; j < atoms[index].degree; ++j) {
    // RDKit❗❌:           unsigned int nbor = atoms[index].nbrIds[j];
    // RDKit❗❌:           touchedPartitions[atoms[nbor].index] = 1;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
    // RDKit❗❌:         if (touchedPartitions[ii]) {
    // RDKit❗❌:           partition = order[ii];
    // RDKit❗❌:           if ((count[partition] > 1) && (next[partition] == -2)) {
    // RDKit❗❌:             next[partition] = activeset;
    // RDKit❗❌:             activeset = partition;
    // RDKit❗❌:           }
    // RDKit❗❌:           touchedPartitions[ii] = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }  // end of RefinePartitions()
    // END COMPLETE RDKit .6 CHEM31 RefinePartitions

    let n_atoms = view.num_atoms();
    let mut local_hanoi_temp;
    let hanoi_temp = match hanoi_temp {
        Some(scratch) => scratch,
        None => {
            local_hanoi_temp = new_hanoi_scratch_for_canonical_rank(n_atoms, false);
            &mut local_hanoi_temp
        }
    };
    // Source validates full-molecule scratch before even inactive activeset.
    if hanoi_temp.len() < n_atoms {
        return Err(CanonicalRankError::HanoiScratchTooSmall);
    }
    #[cfg(test)]
    chem31_scratch_trace::record(chem31_scratch_trace::Event::Refine(
        compare_mode,
        hanoi_temp.as_ptr() as usize,
        hanoi_temp.len(),
    ));
    while *active_set != -1 {
        let partition = usize::try_from(*active_set).expect("active partition is non-negative");
        *active_set = next[partition];
        next[partition] = -2;
        let len = count[partition];
        let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX);
        let result_in_temp = hanoi_order_for_kekulize(
            &mut order[offset..offset + len],
            &mut hanoi_temp[..len],
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        )?;
        #[cfg(test)]
        chem31_scratch_trace::record(chem31_scratch_trace::Event::Sort(
            offset,
            len,
            result_in_temp,
            hanoi_temp.as_ptr() as usize,
        ));
        if result_in_temp {
            order[offset..offset + len].copy_from_slice(&hanoi_temp[..len]);
        }
        for k in 0..len {
            changed[order[offset + k]] = false;
        }
        let mut index = order[offset];
        let mut sym_class = 0usize;
        let mut i = count[index];
        while i < len {
            index = order[offset + i];
            if count[index] > 0 {
                sym_class = offset + i;
            }
            atoms[index].index = i32::try_from(sym_class).unwrap_or(i32::MAX);
            for nbr in atoms[index].nbr_ids.iter().copied() {
                changed[nbr] = true;
            }
            i += 1;
        }
        if mode {
            index = order[offset];
            let mut i = count[index];
            while i < len {
                index = order[offset + i];
                for nbr in atoms[index].nbr_ids.iter().copied() {
                    let partition_idx = usize::try_from(atoms[nbr].index).unwrap_or(usize::MAX);
                    if partition_idx < touched_partitions.len() {
                        touched_partitions[partition_idx] = true;
                    }
                }
                i += 1;
            }
            for ii in 0..n_atoms {
                if touched_partitions[ii] {
                    let partition = order[ii];
                    if count[partition] > 1 && next[partition] == -2 {
                        next[partition] = *active_set;
                        *active_set = isize::try_from(partition).unwrap_or(isize::MAX);
                    }
                    touched_partitions[ii] = false;
                }
            }
        }
    }
    Ok(())
}

fn break_ties_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    mode: bool,
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
    order: &mut [usize],
    count: &mut [usize],
    active_set: &mut isize,
    next: &mut [isize],
    changed: &mut [bool],
    touched_partitions: &mut [bool],
    mut hanoi_temp: Option<&mut [usize]>,
) -> Result<(), CanonicalRankError> {
    // BEGIN COMPLETE RDKit .6 CHEM31 BreakTies
    // RDKit❗❌: template <typename CompareFunc>
    // RDKit❗❌: void BreakTies(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit❗❌:                int mode, std::vector<int> &order, std::vector<int> &count,
    // RDKit❗❌:                int &activeset, std::vector<int> &next,
    // RDKit❗❌:                std::vector<int> &changed,
    // RDKit❗❌:                std::vector<char> &touchedPartitions,
    // RDKit❗❌:                std::vector<int> *hanoiTemp = nullptr) {
    // RDKit❗❌:   unsigned int nAtoms = mol.getNumAtoms();
    // RDKit❗❌:   int partition;
    // RDKit❗❌:   int offset;
    // RDKit❗❌:   int index;
    // RDKit❗❌:   int len;
    // RDKit❗❌:   int oldPart = 0;
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int i = 0; i < nAtoms; i++) {
    // RDKit❗❌:     partition = order[i];
    // RDKit❗❌:     oldPart = atoms[partition].index;
    // RDKit❗❌:     while (count[partition] > 1) {
    // RDKit❗❌:       len = count[partition];
    // RDKit❗❌:       offset = atoms[partition].index + len - 1;
    // RDKit❗❌:       index = order[offset];
    // RDKit❗❌:       atoms[index].index = offset;
    // RDKit❗❌:       count[partition] = len - 1;
    // RDKit❗❌:       count[index] = 1;
    // RDKit❗❌:
    // RDKit❗❌:       // test for ions, water molecules with no
    // RDKit❗❌:       if (atoms[index].degree < 1) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
    // RDKit❗❌:         unsigned int nbor = atoms[index].nbrIds[j];
    // RDKit❗❌:         touchedPartitions[atoms[nbor].index] = 1;
    // RDKit❗❌:         changed[nbor] = 1;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
    // RDKit❗❌:         if (touchedPartitions[ii]) {
    // RDKit❗❌:           int npart = order[ii];
    // RDKit❗❌:           if ((count[npart] > 1) && (next[npart] == -2)) {
    // RDKit❗❌:             next[npart] = activeset;
    // RDKit❗❌:             activeset = npart;
    // RDKit❗❌:           }
    // RDKit❗❌:           touchedPartitions[ii] = 0;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       RefinePartitions(mol, atoms, compar, mode, order, count, activeset, next,
    // RDKit❗❌:                        changed, touchedPartitions, hanoiTemp);
    // RDKit❗❌:     }
    // RDKit❗❌:     // not sure if this works each time
    // RDKit❗❌:     if (atoms[partition].index != oldPart) {
    // RDKit❗❌:       i -= 1;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }  // end of BreakTies()
    // END COMPLETE RDKit .6 CHEM31 BreakTies

    let n_atoms = view.num_atoms();
    let mut i = 0usize;
    while i < n_atoms {
        let partition = order[i];
        let old_part = atoms[partition].index;
        while count[partition] > 1 {
            let len = count[partition];
            let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX) + len - 1;
            let index = order[offset];
            atoms[index].index = i32::try_from(offset).unwrap_or(i32::MAX);
            count[partition] = len - 1;
            count[index] = 1;
            if atoms[index].degree < 1 {
                continue;
            }
            for nbr in atoms[index].nbr_ids.iter().copied() {
                let partition_idx = usize::try_from(atoms[nbr].index).unwrap_or(usize::MAX);
                if partition_idx < touched_partitions.len() {
                    touched_partitions[partition_idx] = true;
                }
                changed[nbr] = true;
            }
            for ii in 0..n_atoms {
                if touched_partitions[ii] {
                    let npart = order[ii];
                    if count[npart] > 1 && next[npart] == -2 {
                        next[npart] = *active_set;
                        *active_set = isize::try_from(npart).unwrap_or(isize::MAX);
                    }
                    touched_partitions[ii] = false;
                }
            }
            refine_partitions_for_kekulize(
                view,
                atoms,
                mode,
                compare_mode,
                flags,
                order,
                count,
                active_set,
                next,
                changed,
                touched_partitions,
                hanoi_temp.as_deref_mut(),
            )?;
        }
        if atoms[partition].index != old_part {
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }
    Ok(())
}

fn hanoi_sort_order_for_kekulize(
    order: &mut [usize],
    offset: usize,
    len: usize,
    count: &mut [usize],
    changed: &[bool],
    atoms: &mut [CanonAtom<'_>],
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> Result<(), CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION hanoisort
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void hanoisort(std::span<int> &base, std::vector<int> &count,
    // RDKit✔️✔️:                std::vector<int> &changed, CompareFunc compar) {
    // RDKit✔️✔️:   std::vector<int> tempVec(base.size());
    let mut temp = vec![0usize; len];
    // RDKit✔️✔️:   if (detail::hanoi(base.data(), base.size(), tempVec.data(), count.data(),
    // RDKit✔️✔️:                     changed.data(), compar)) {
    if hanoi_order_for_kekulize(
        &mut order[offset..offset + len],
        &mut temp,
        count,
        changed,
        atoms,
        compare_mode,
        flags,
    )? {
        // RDKit✔️✔️:     std::copy(tempVec.begin(), tempVec.end(), base.begin());
        order[offset..offset + len].copy_from_slice(&temp);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hanoisort
    Ok(())
}

fn hanoi_order_for_kekulize(
    base: &mut [usize],
    temp: &mut [usize],
    count: &mut [usize],
    changed: &[bool],
    atoms: &mut [CanonAtom<'_>],
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> Result<bool, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION detail::hanoi
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: bool hanoi(int *base, int nel, int *temp, int *count, int *changed,
    // RDKit✔️✔️:            CompareFunc compar) {
    // RDKit✔️✔️:   assert(base);
    // RDKit✔️✔️:   assert(temp);
    // RDKit✔️✔️:   assert(count);
    // RDKit✔️✔️:   assert(changed);
    debug_assert_eq!(base.len(), temp.len());
    // RDKit✔️✔️:   int *b1, *b2;
    // RDKit✔️✔️:   int *t1, *t2;
    // RDKit✔️✔️:   int *s1, *s2;
    // RDKit✔️✔️:   int n1, n2;
    // RDKit✔️✔️:   int result;
    // RDKit✔️✔️:   int *ptr;
    let nel = base.len();

    // RDKit✔️✔️:   if (nel == 1) {
    // RDKit✔️✔️:     count[base[0]] = 1;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   } else if (nel == 2) {
    if nel == 1 {
        count[base[0]] = 1;
        return Ok(false);
    } else if nel == 2 {
        // RDKit✔️✔️:     n1 = base[0];
        // RDKit✔️✔️:     n2 = base[1];
        let n1 = base[0];
        let n2 = base[1];
        // RDKit✔️✔️:     int stat =
        // RDKit✔️✔️:         (/*!changed || */ changed[n1] || changed[n2]) ? compar(n1, n2) : 0;
        let stat = if changed[n1] || changed[n2] {
            compare_canon_atoms_for_kekulize(atoms, n1, n2, compare_mode, flags)?
        } else {
            Ordering::Equal
        };
        // RDKit✔️✔️:     if (stat == 0) {
        // RDKit✔️✔️:       count[n1] = 2;
        // RDKit✔️✔️:       count[n2] = 0;
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     } else if (stat < 0) {
        // RDKit✔️✔️:       count[n1] = 1;
        // RDKit✔️✔️:       count[n2] = 1;
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     } else /* stat > 0 */ {
        // RDKit✔️✔️:       count[n1] = 1;
        // RDKit✔️✔️:       count[n2] = 1;
        // RDKit✔️✔️:       base[0] = n2; /* temp[0] = n2; */
        // RDKit✔️✔️:       base[1] = n1; /* temp[1] = n1; */
        // RDKit✔️✔️:       return false; /* return True;  */
        // RDKit✔️✔️:     }
        match stat {
            Ordering::Equal => {
                count[n1] = 2;
                count[n2] = 0;
            }
            Ordering::Less => {
                count[n1] = 1;
                count[n2] = 1;
            }
            Ordering::Greater => {
                count[n1] = 1;
                count[n2] = 1;
                base[0] = n2;
                base[1] = n1;
            }
        }
        return Ok(false);
    }

    // RDKit✔️✔️:   n1 = nel / 2;
    // RDKit✔️✔️:   n2 = nel - n1;
    let left_len = nel / 2;
    let right_len = nel - left_len;
    // RDKit✔️✔️:   b1 = base;
    // RDKit✔️✔️:   t1 = temp;
    // RDKit✔️✔️:   b2 = base + n1;
    // RDKit✔️✔️:   t2 = temp + n1;
    let (left_in_temp, right_in_temp) = {
        let (base_left, base_right) = base.split_at_mut(left_len);
        let (temp_left, temp_right) = temp.split_at_mut(left_len);

        // RDKit✔️✔️:   if (hanoi(b1, n1, t1, count, changed, compar)) {
        let left_in_temp = hanoi_order_for_kekulize(
            base_left,
            temp_left,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        )?;
        // RDKit✔️✔️:     if (hanoi(b2, n2, t2, count, changed, compar)) {
        let right_in_temp = hanoi_order_for_kekulize(
            base_right,
            temp_right,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        )?;
        (left_in_temp, right_in_temp)
    };
    // RDKit✔️✔️:     s1 = t1;
    // RDKit✔️✔️:     s1 = b1;
    // RDKit✔️✔️:       s2 = t2;
    // RDKit✔️✔️:       s2 = b2;
    // RDKit✔️✔️:     result = false;
    // RDKit✔️✔️:     ptr = base;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     result = true;
    // RDKit✔️✔️:     ptr = temp;
    // RDKit✔️✔️:   }
    let result = !left_in_temp;
    let base_ptr = base.as_ptr();
    let temp_ptr = temp.as_ptr();
    let left_source_ptr = if left_in_temp { temp_ptr } else { base_ptr };
    let right_source_ptr = if right_in_temp {
        // SAFETY: right half starts at left_len within the same backing buffer.
        unsafe { temp_ptr.add(left_len) }
    } else {
        // SAFETY: right half starts at left_len within the same backing buffer.
        unsafe { base_ptr.add(left_len) }
    };
    let output_ptr = if result {
        temp.as_mut_ptr()
    } else {
        base.as_mut_ptr()
    };
    let mut left_pos = 0usize;
    let mut right_pos = 0usize;
    let mut out_pos = 0usize;
    let mut remaining_left = left_len;
    let mut remaining_right = right_len;

    // RDKit✔️✔️:   while (true) {
    loop {
        let left_atom = unsafe { *left_source_ptr.add(left_pos) };
        let right_atom = unsafe { *right_source_ptr.add(right_pos) };
        // RDKit✔️✔️:     assert(*s1 != *s2);
        debug_assert_ne!(left_atom, right_atom);
        // RDKit✔️✔️:     int stat =
        // RDKit✔️✔️:         (/*!changed || */ changed[*s1] || changed[*s2]) ? compar(*s1, *s2) : 0;
        let stat = if changed[left_atom] || changed[right_atom] {
            compare_canon_atoms_for_kekulize(atoms, left_atom, right_atom, compare_mode, flags)?
        } else {
            Ordering::Equal
        };
        // RDKit✔️✔️:     int len1 = count[*s1];
        // RDKit✔️✔️:     int len2 = count[*s2];
        // RDKit✔️✔️:     assert(len1 > 0);
        // RDKit✔️✔️:     assert(len2 > 0);
        let class_left = count[left_atom];
        let class_right = count[right_atom];
        debug_assert!(class_left > 0);
        debug_assert!(class_right > 0);
        // RDKit✔️✔️:     if (stat == 0) {
        if stat == Ordering::Equal {
            // RDKit✔️✔️:       count[*s1] = len1 + len2;
            // RDKit✔️✔️:       count[*s2] = 0;
            count[left_atom] = class_left + class_right;
            count[right_atom] = 0;
            // RDKit✔️✔️:       memmove(ptr, s1, len1 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    left_source_ptr.add(left_pos),
                    output_ptr.add(out_pos),
                    class_left,
                );
            }
            // RDKit✔️✔️:       ptr += len1;
            // RDKit✔️✔️:       n1 -= len1;
            out_pos += class_left;
            remaining_left -= class_left;
            // RDKit✔️✔️:       if (n1 == 0) {
            // RDKit✔️✔️:         if (ptr != s2) {
            // RDKit✔️✔️:           memmove(ptr, s2, n2 * sizeof(int));
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_left == 0 {
                unsafe {
                    std::ptr::copy(
                        right_source_ptr.add(right_pos),
                        output_ptr.add(out_pos),
                        remaining_right,
                    );
                }
                return Ok(result);
            }
            // RDKit✔️✔️:       s1 += len1;
            left_pos += class_left;
            // RDKit✔️✔️:       memmove(ptr, s2, len2 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    right_source_ptr.add(right_pos),
                    output_ptr.add(out_pos),
                    class_right,
                );
            }
            // RDKit✔️✔️:       ptr += len2;
            // RDKit✔️✔️:       n2 -= len2;
            out_pos += class_right;
            remaining_right -= class_right;
            // RDKit✔️✔️:       if (n2 == 0) {
            // RDKit✔️✔️:         memmove(ptr, s1, n1 * sizeof(int));
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_right == 0 {
                unsafe {
                    std::ptr::copy(
                        left_source_ptr.add(left_pos),
                        output_ptr.add(out_pos),
                        remaining_left,
                    );
                }
                return Ok(result);
            }
            // RDKit✔️✔️:       s2 += len2;
            right_pos += class_right;
            // RDKit✔️✔️:     } else if (stat < 0 && len1 > 0) {
        } else if stat == Ordering::Less {
            // RDKit✔️✔️:       memmove(ptr, s1, len1 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    left_source_ptr.add(left_pos),
                    output_ptr.add(out_pos),
                    class_left,
                );
            }
            // RDKit✔️✔️:       ptr += len1;
            // RDKit✔️✔️:       n1 -= len1;
            out_pos += class_left;
            remaining_left -= class_left;
            // RDKit✔️✔️:       if (n1 == 0) {
            // RDKit✔️✔️:         if (ptr != s2) {
            // RDKit✔️✔️:           memmove(ptr, s2, n2 * sizeof(int));
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_left == 0 {
                unsafe {
                    std::ptr::copy(
                        right_source_ptr.add(right_pos),
                        output_ptr.add(out_pos),
                        remaining_right,
                    );
                }
                return Ok(result);
            }
            // RDKit✔️✔️:       s1 += len1;
            left_pos += class_left;
            // RDKit✔️✔️:     } else if (stat > 0 && len2 > 0) /* stat > 0 */ {
        } else {
            // RDKit✔️✔️:       memmove(ptr, s2, len2 * sizeof(int));
            unsafe {
                std::ptr::copy(
                    right_source_ptr.add(right_pos),
                    output_ptr.add(out_pos),
                    class_right,
                );
            }
            // RDKit✔️✔️:       ptr += len2;
            // RDKit✔️✔️:       n2 -= len2;
            out_pos += class_right;
            remaining_right -= class_right;
            // RDKit✔️✔️:       if (n2 == 0) {
            // RDKit✔️✔️:         memmove(ptr, s1, n1 * sizeof(int));
            // RDKit✔️✔️:         return result;
            // RDKit✔️✔️:       }
            if remaining_right == 0 {
                unsafe {
                    std::ptr::copy(
                        left_source_ptr.add(left_pos),
                        output_ptr.add(out_pos),
                        remaining_left,
                    );
                }
                return Ok(result);
            }
            // RDKit✔️✔️:       s2 += len2;
            right_pos += class_right;
            // RDKit✔️✔️:     } else {
            // RDKit✔️✔️:       assert(0);
            // RDKit✔️✔️:     }
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detail::hanoi
}

fn compare_canon_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
    mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> Result<Ordering, CanonicalRankError> {
    if matches!(mode, CanonCompareMode::SpecialSymmetry) {
        return Ok(compare_special_symmetry_atoms_for_kekulize(
            atoms, left, right,
        ));
    }
    if matches!(mode, CanonCompareMode::SpecialChirality) {
        return Ok(compare_special_chirality_atoms_for_kekulize(
            atoms, left, right,
        ));
    }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ok(Ordering::Equal);
    }
    // RDKit✔️✔️:     int v = basecomp(i, j);
    // RDKit✔️✔️:     if (v) {
    // RDKit✔️✔️:       return v;
    // RDKit✔️✔️:     }
    let base_cmp = compare_canon_atom_base_for_kekulize(atoms, left, right, flags)?;
    if base_cmp != Ordering::Equal {
        return Ok(base_cmp);
    }
    // RDKit✔️✔️:     if (df_useNbrs) {
    // RDKit✔️✔️:       if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    // RDKit✔️✔️:       }
    let ranks = canon_atom_rank_snapshot(atoms);
    // RDKit✔️✔️:       for (unsigned int ii = 0;
    // RDKit✔️✔️:            ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size();
    // RDKit✔️✔️:            ++ii) {
    // RDKit✔️✔️:         int cmp =
    // RDKit✔️✔️:             bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:         if (cmp) {
    // RDKit✔️✔️:           return cmp;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }
    // RDKit✔️✔️:       if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    Ok(atoms[left].bonds.len().cmp(&atoms[right].bonds.len()))
}

fn compare_special_chirality_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION SpecialChiralityAtomCompareFunctor::operator()
    // RDKit✔️✔️:   int operator()(int i, int j) const {
    // RDKit✔️✔️:     PRECONDITION(dp_atoms, "no atoms");
    // RDKit✔️✔️:     PRECONDITION(dp_mol, "no molecule");
    // RDKit✔️✔️:     PRECONDITION(i != j, "bad call");
    // RDKit✔️✔️:     if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    // RDKit✔️✔️:     }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    // RDKit✔️✔️:     for (unsigned int ii = 0;
    // RDKit✔️✔️:          ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
    // RDKit✔️✔️:       int cmp =
    // RDKit✔️✔️:           bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    let ranks = canon_atom_rank_snapshot(atoms);
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> swapsi;
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> swapsj;
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[i].bonds, i, swapsi);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborNumSwaps(dp_atoms, dp_atoms[j].bonds, j, swapsj);
    // RDKit✔️✔️:     }
    let swaps_left = if atoms[left].is_in_play {
        update_atom_neighbor_num_swaps_for_kekulize(atoms, left)
    } else {
        Vec::new()
    };
    let swaps_right = if atoms[right].is_in_play {
        update_atom_neighbor_num_swaps_for_kekulize(atoms, right)
    } else {
        Vec::new()
    };
    // RDKit✔️✔️:     for (unsigned int ii = 0; ii < swapsi.size() && ii < swapsj.size(); ++ii) {
    // RDKit✔️✔️:       int cmp = swapsi[ii].second - swapsj[ii].second;
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for idx in 0..swaps_left.len().min(swaps_right.len()) {
        let cmp = swaps_left[idx].1.cmp(&swaps_right[idx].1);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SpecialChiralityAtomCompareFunctor::operator()
    Ordering::Equal
}

fn compare_special_symmetry_atoms_for_kekulize(
    atoms: &mut [CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION SpecialSymmetryAtomCompareFunctor::operator()
    // RDKit✔️✔️:   int operator()(int i, int j) const {
    // RDKit✔️✔️:     PRECONDITION(dp_atoms, "no atoms");
    // RDKit✔️✔️:     PRECONDITION(dp_mol, "no molecule");
    // RDKit✔️✔️:     PRECONDITION(i != j, "bad call");
    // RDKit✔️✔️:     if (dp_atomsInPlay && !((*dp_atomsInPlay)[i] || (*dp_atomsInPlay)[j])) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (dp_atoms[i].neighborNum < dp_atoms[j].neighborNum) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].neighborNum > dp_atoms[j].neighborNum) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    let neighbor_cmp = atoms[left].neighbor_num.cmp(&atoms[right].neighbor_num);
    if neighbor_cmp != Ordering::Equal {
        return neighbor_cmp;
    }
    // RDKit✔️✔️:     if (dp_atoms[i].revistedNeighbors < dp_atoms[j].revistedNeighbors) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].revistedNeighbors > dp_atoms[j].revistedNeighbors) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    let revisited_cmp = atoms[left]
        .revisted_neighbors
        .cmp(&atoms[right].revisted_neighbors);
    if revisited_cmp != Ordering::Equal {
        return revisited_cmp;
    }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[i]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[i].bonds);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!dp_atomsInPlay || (*dp_atomsInPlay)[j]) {
    // RDKit✔️✔️:       updateAtomNeighborIndex(dp_atoms, dp_atoms[j].bonds);
    // RDKit✔️✔️:     }
    if atoms[left].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, left);
    }
    if atoms[right].is_in_play {
        update_atom_neighbor_index_for_kekulize(atoms, right);
    }
    let ranks = canon_atom_rank_snapshot(atoms);
    // RDKit✔️✔️:     for (unsigned int ii = 0;
    // RDKit✔️✔️:          ii < dp_atoms[i].bonds.size() && ii < dp_atoms[j].bonds.size(); ++ii) {
    // RDKit✔️✔️:       int cmp =
    // RDKit✔️✔️:           bondholder::compare(dp_atoms[i].bonds[ii], dp_atoms[j].bonds[ii]);
    // RDKit✔️✔️:       if (cmp) {
    // RDKit✔️✔️:         return cmp;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    for idx in 0..atoms[left].bonds.len().min(atoms[right].bonds.len()) {
        let cmp =
            compare_canon_bond_holder(&atoms[left].bonds[idx], &atoms[right].bonds[idx], &ranks);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    // RDKit✔️✔️:     if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION SpecialSymmetryAtomCompareFunctor::operator()
    atoms[left].bonds.len().cmp(&atoms[right].bonds.len())
}

fn canonical_rank_property_to_int(
    atom: &Atom,
    atom_index: usize,
) -> Result<i32, CanonicalRankError> {
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:348-363
    // Boost❗✔️:   template<class Traits,class OverflowHandler>
    // Boost❗✔️:   struct GetRC_Int2Int
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef GetRC_Sig2Sig_or_Unsig2Unsig<Traits,OverflowHandler> Sig2SigQ     ;
    // Boost❗✔️:     typedef GetRC_Sig2Unsig             <Traits,OverflowHandler> Sig2UnsigQ   ;
    // Boost❗✔️:     typedef GetRC_Unsig2Sig             <Traits,OverflowHandler> Unsig2SigQ   ;
    // Boost❗✔️:     typedef Sig2SigQ                                             Unsig2UnsigQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::sign_mixture sign_mixture ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename
    // Boost❗✔️:       for_sign_mixture<sign_mixture,Sig2SigQ,Sig2UnsigQ,Unsig2SigQ,Unsig2UnsigQ>::type
    // Boost❗✔️:         selector ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename selector::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:348-363
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:406-419
    // Boost❗✔️:   template<class Traits, class OverflowHandler, class Float2IntRounder>
    // Boost❗✔️:   struct GetRC_BuiltIn2BuiltIn
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef GetRC_Int2Int<Traits,OverflowHandler>                    Int2IntQ ;
    // Boost❗✔️:     typedef GetRC_Int2Float<Traits>                                  Int2FloatQ ;
    // Boost❗✔️:     typedef GetRC_Float2Int<Traits,OverflowHandler,Float2IntRounder> Float2IntQ ;
    // Boost❗✔️:     typedef GetRC_Float2Float<Traits,OverflowHandler>                Float2FloatQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::int_float_mixture int_float_mixture ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename for_int_float_mixture<int_float_mixture, Int2IntQ, Int2FloatQ, Float2IntQ, Float2FloatQ>::type selector ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename selector::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:406-419
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:421-435
    // Boost❗✔️:   template<class Traits, class OverflowHandler, class Float2IntRounder>
    // Boost❗✔️:   struct GetRC
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef GetRC_BuiltIn2BuiltIn<Traits,OverflowHandler,Float2IntRounder> BuiltIn2BuiltInQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef dummy_range_checker<Traits> Dummy ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef mpl::identity<Dummy> DummyQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::udt_builtin_mixture udt_builtin_mixture ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename for_udt_builtin_mixture<udt_builtin_mixture,BuiltIn2BuiltInQ,DummyQ,DummyQ,DummyQ>::type selector ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename selector::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:421-435
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:528-554
    // Boost❗✔️:   template<class Traits,class OverflowHandler,class Float2IntRounder,class RawConverter, class UserRangeChecker>
    // Boost❗✔️:   struct get_non_trivial_converter
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef GetRC<Traits,OverflowHandler,Float2IntRounder> InternalRangeCheckerQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef is_same<UserRangeChecker,UseInternalRangeChecker> use_internal_RC ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef mpl::identity<UserRangeChecker> UserRangeCheckerQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename
    // Boost❗✔️:       mpl::eval_if<use_internal_RC,InternalRangeCheckerQ,UserRangeCheckerQ>::type
    // Boost❗✔️:         RangeChecker ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef non_rounding_converter<Traits,RangeChecker,RawConverter>              NonRounding ;
    // Boost❗✔️:     typedef rounding_converter<Traits,RangeChecker,RawConverter,Float2IntRounder> Rounding ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef mpl::identity<NonRounding> NonRoundingQ ;
    // Boost❗✔️:     typedef mpl::identity<Rounding>    RoundingQ    ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::int_float_mixture int_float_mixture ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename
    // Boost❗✔️:       for_int_float_mixture<int_float_mixture, NonRoundingQ, NonRoundingQ, RoundingQ, NonRoundingQ>::type
    // Boost❗✔️:         selector ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename selector::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:528-554
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:556-587
    // Boost❗✔️:   template< class Traits
    // Boost❗✔️:            ,class OverflowHandler
    // Boost❗✔️:            ,class Float2IntRounder
    // Boost❗✔️:            ,class RawConverter
    // Boost❗✔️:            ,class UserRangeChecker
    // Boost❗✔️:           >
    // Boost❗✔️:   struct get_converter_impl
    // Boost❗✔️:   {
    // Boost❗✔️: #if BOOST_WORKAROUND(BOOST_BORLANDC, BOOST_TESTED_AT( 0x0561 ) )
    // Boost❗✔️:     // bcc55 prefers sometimes template parameters to be explicit local types.
    // Boost❗✔️:     // (notice that is is illegal to reuse the names like this)
    // Boost❗✔️:     typedef Traits           Traits ;
    // Boost❗✔️:     typedef OverflowHandler  OverflowHandler ;
    // Boost❗✔️:     typedef Float2IntRounder Float2IntRounder ;
    // Boost❗✔️:     typedef RawConverter     RawConverter ;
    // Boost❗✔️:     typedef UserRangeChecker UserRangeChecker ;
    // Boost❗✔️: #endif
    // Boost❗✔️:
    // Boost❗✔️:     typedef trivial_converter_impl<Traits> Trivial ;
    // Boost❗✔️:     typedef mpl::identity        <Trivial> TrivialQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef get_non_trivial_converter< Traits
    // Boost❗✔️:                                       ,OverflowHandler
    // Boost❗✔️:                                       ,Float2IntRounder
    // Boost❗✔️:                                       ,RawConverter
    // Boost❗✔️:                                       ,UserRangeChecker
    // Boost❗✔️:                                      > NonTrivialQ ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::trivial trivial ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename mpl::eval_if<trivial,TrivialQ,NonTrivialQ>::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:556-587
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:181-188
    // Boost❗✔️: template<class Traits>
    // Boost❗✔️: struct raw_converter
    // Boost❗✔️: {
    // Boost❗✔️:   typedef typename Traits::result_type   result_type   ;
    // Boost❗✔️:   typedef typename Traits::argument_type argument_type ;
    // Boost❗✔️:
    // Boost❗✔️:   static result_type low_level_convert ( argument_type s ) { return static_cast<result_type>(s) ; }
    // Boost❗✔️: } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:181-188
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/conversion_traits.hpp:30-49
    // Boost❗✔️:   template<class T,class S>
    // Boost❗✔️:   struct non_trivial_traits_impl
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef typename get_int_float_mixture   <T,S>::type int_float_mixture ;
    // Boost❗✔️:     typedef typename get_sign_mixture        <T,S>::type sign_mixture ;
    // Boost❗✔️:     typedef typename get_udt_builtin_mixture <T,S>::type udt_builtin_mixture ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename get_is_subranged<T,S>::type subranged ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef mpl::false_ trivial ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef T target_type ;
    // Boost❗✔️:     typedef S source_type ;
    // Boost❗✔️:     typedef T result_type ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename mpl::if_< is_arithmetic<S>, S, S const&>::type argument_type ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename mpl::if_<subranged,S,T>::type supertype ;
    // Boost❗✔️:     typedef typename mpl::if_<subranged,T,S>::type subtype   ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/conversion_traits.hpp:30-49
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/conversion_traits.hpp:79-91
    // Boost❗✔️:   template<class T, class S>
    // Boost❗✔️:   struct get_conversion_traits
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef typename remove_cv<T>::type target_type ;
    // Boost❗✔️:     typedef typename remove_cv<S>::type source_type ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename is_same<target_type,source_type>::type is_trivial ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef trivial_traits_impl    <target_type>             trivial_imp ;
    // Boost❗✔️:     typedef non_trivial_traits_impl<target_type,source_type> non_trivial_imp ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename mpl::if_<is_trivial,trivial_imp,non_trivial_imp>::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/conversion_traits.hpp:79-91
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/bounds.hpp:43-52
    // Boost❗✔️:   template<class N>
    // Boost❗✔️:   struct get_impl
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef mpl::bool_< ::std::numeric_limits<N>::is_integer > is_int ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef Integral<N> impl_int   ;
    // Boost❗✔️:     typedef Float   <N> impl_float ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename mpl::if_<is_int,impl_int,impl_float>::type type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/bounds.hpp:43-52

    // Proposed domain: UInt/Int/String/Double/Bool/IntVector and absent property,
    // with the default classical C++ locale ratified by ROOT. The Uint tag
    // is proposed below; custom global C++ locale facets remain unmodeled.
    // UInt source bounds and getter order were validated by the actual strict release cases;
    // implementation and overflow projection require approval and tests.
    // Behavior: missing initializes0, exact Int returns unchanged; String
    // right-trims only six C spaces before a complete signed decimal parse;
    // every reached invalid cast remains a structured source bad_any_cast.
    // Cost: a borrowed BTreeMap lookup is O(log properties), versus source
    // Dict's linear search; no property clone or eager initialization parse.
    // String conversion is O(bytes) with constant extra storage, avoiding
    // the source String copy and LocaleSwitcher allocation under this fixed
    // locale. These savings justify the improved cost markers below.
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDProps.h:126-129
    // RDKit✔️🔝:   template <typename T>
    // RDKit✔️🔝:   bool getPropIfPresent(const std::string_view key, T &res) const {
    // RDKit✔️🔝:     return d_props.getValIfPresent(key, res);
    // RDKit✔️🔝:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDProps.h:126-129
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/Dict.h:255-264
    // RDKit✔️🔝:   template <typename T>
    // RDKit✔️🔝:   bool getValIfPresent(const std::string_view what, T &res) const {
    // RDKit✔️🔝:     for (const auto &data : _data) {
    // RDKit✔️🔝:       if (data.key == what) {
    // RDKit✔️🔝:         res = from_rdvalue<T>(data.val);
    // RDKit✔️🔝:         return true;
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     return false;
    // RDKit✔️🔝:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/Dict.h:255-264
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:268-292
    // RDKit✔️🔝: // from_rdvalue -> converts string values to appropriate types
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: typename boost::enable_if<boost::is_arithmetic<T>, T>::type from_rdvalue(
    // RDKit✔️🔝:     RDValue_cast_t arg) {
    // RDKit✔️🔝:   T res;
    // RDKit✔️🔝:   if (arg.getTag() == RDTypeTag::StringTag) {
    // RDKit✔️🔝:     Utils::LocaleSwitcher ls;
    // RDKit✔️🔝:     try {
    // RDKit✔️🔝:       res = rdvalue_cast<T>(arg);
    // RDKit✔️🔝:     } catch (const std::bad_any_cast &exc) {
    // RDKit✔️🔝:       try {
    // RDKit✔️🔝: 	std::string val = rdvalue_cast<std::string>(arg);
    // RDKit✔️🔝: 	// trim only the right characters, this mimics how SD values
    // RDKit✔️🔝: 	//  work on read, they will be trimmed by the MolFile parser
    // RDKit✔️🔝: 	boost::trim_right(val);
    // RDKit✔️🔝:         res = boost::lexical_cast<T>(val);
    // RDKit✔️🔝:       } catch (...) {
    // RDKit✔️🔝:         throw exc;
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     res = rdvalue_cast<T>(arg);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:268-292
    // Int tag extraction is O(1). The UInt conversion checks the
    // signed upper bound and preserves typed positive_overflow; actual cases cover all six boundaries.
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450
    // RDKit✔️✔️: template <>
    // RDKit✔️✔️: inline int rdvalue_cast<int>(RDValue_cast_t v) {
    // RDKit✔️✔️:   if (rdvalue_is<int>(v)) {
    // RDKit✔️✔️:     return v.value.i;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (rdvalue_is<unsigned int>(v)) {
    // RDKit✔️✔️:     return boost::numeric_cast<int>(v.value.u);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   throw std::bad_any_cast();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue-taggedunion.h:441-450
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast.hpp:36-46
    // Boost✔️✔️:     template <typename Target, typename Source>
    // Boost✔️✔️:     inline Target lexical_cast(const Source &arg)
    // Boost✔️✔️:     {
    // Boost✔️✔️:         Target result = Target();
    // Boost✔️✔️:
    // Boost✔️✔️:         if (!boost::conversion::detail::try_lexical_convert(arg, result)) {
    // Boost✔️✔️:             boost::conversion::detail::throw_bad_cast<Source, Target>();
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         return result;
    // Boost✔️✔️:     }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast.hpp:36-46
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/try_lexical_convert.hpp:164-202
    // Boost✔️✔️:         template <typename Target, typename Source>
    // Boost✔️✔️:         inline bool try_lexical_convert(const Source& arg, Target& result)
    // Boost✔️✔️:         {
    // Boost✔️✔️:             typedef BOOST_DEDUCED_TYPENAME boost::detail::array_to_pointer_decay<Source>::type src;
    // Boost✔️✔️:
    // Boost✔️✔️:             typedef boost::integral_constant<
    // Boost✔️✔️:                 bool,
    // Boost✔️✔️:                 boost::detail::is_xchar_to_xchar<Target, src >::value ||
    // Boost✔️✔️:                 boost::detail::is_char_array_to_stdstring<Target, src >::value ||
    // Boost✔️✔️:                 boost::detail::is_char_array_to_booststring<Target, src >::value ||
    // Boost✔️✔️:                 (
    // Boost✔️✔️:                      boost::is_same<Target, src >::value &&
    // Boost✔️✔️:                      (boost::detail::is_stdstring<Target >::value || boost::detail::is_booststring<Target >::value)
    // Boost✔️✔️:                 ) ||
    // Boost✔️✔️:                 (
    // Boost✔️✔️:                      boost::is_same<Target, src >::value &&
    // Boost✔️✔️:                      boost::detail::is_character<Target >::value
    // Boost✔️✔️:                 )
    // Boost✔️✔️:             > shall_we_copy_t;
    // Boost✔️✔️:
    // Boost✔️✔️:             typedef boost::detail::is_arithmetic_and_not_xchars<Target, src >
    // Boost✔️✔️:                 shall_we_copy_with_dynamic_check_t;
    // Boost✔️✔️:
    // Boost✔️✔️:             // We do evaluate second `if_` lazily to avoid unnecessary instantiations
    // Boost✔️✔️:             // of `shall_we_copy_with_dynamic_check_t` and improve compilation times.
    // Boost✔️✔️:             typedef BOOST_DEDUCED_TYPENAME boost::conditional<
    // Boost✔️✔️:                 shall_we_copy_t::value,
    // Boost✔️✔️:                 boost::type_identity<boost::detail::copy_converter_impl<Target, src > >,
    // Boost✔️✔️:                 boost::conditional<
    // Boost✔️✔️:                      shall_we_copy_with_dynamic_check_t::value,
    // Boost✔️✔️:                      boost::detail::dynamic_num_converter_impl<Target, src >,
    // Boost✔️✔️:                      boost::detail::lexical_converter_impl<Target, src >
    // Boost✔️✔️:                 >
    // Boost✔️✔️:             >::type caster_type_lazy;
    // Boost✔️✔️:
    // Boost✔️✔️:             typedef BOOST_DEDUCED_TYPENAME caster_type_lazy::type caster_type;
    // Boost✔️✔️:
    // Boost✔️✔️:             return caster_type::try_convert(arg, result);
    // Boost✔️✔️:         }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/try_lexical_convert.hpp:164-202
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical.hpp:458-490
    // Boost✔️✔️:         template<typename Target, typename Source>
    // Boost✔️✔️:         struct lexical_converter_impl
    // Boost✔️✔️:         {
    // Boost✔️✔️:             typedef lexical_cast_stream_traits<Source, Target>  stream_trait;
    // Boost✔️✔️:
    // Boost✔️✔️:             typedef detail::lexical_istream_limited_src<
    // Boost✔️✔️:                 BOOST_DEDUCED_TYPENAME stream_trait::char_type,
    // Boost✔️✔️:                 BOOST_DEDUCED_TYPENAME stream_trait::traits,
    // Boost✔️✔️:                 stream_trait::requires_stringbuf,
    // Boost✔️✔️:                 stream_trait::len_t::value + 1
    // Boost✔️✔️:             > i_interpreter_type;
    // Boost✔️✔️:
    // Boost✔️✔️:             typedef detail::lexical_ostream_limited_src<
    // Boost✔️✔️:                 BOOST_DEDUCED_TYPENAME stream_trait::char_type,
    // Boost✔️✔️:                 BOOST_DEDUCED_TYPENAME stream_trait::traits
    // Boost✔️✔️:             > o_interpreter_type;
    // Boost✔️✔️:
    // Boost✔️✔️:             static inline bool try_convert(const Source& arg, Target& result) {
    // Boost✔️✔️:                 i_interpreter_type i_interpreter;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Disabling ADL, by directly specifying operators.
    // Boost✔️✔️:                 if (!(i_interpreter.operator <<(arg)))
    // Boost✔️✔️:                     return false;
    // Boost✔️✔️:
    // Boost✔️✔️:                 o_interpreter_type out(i_interpreter.cbegin(), i_interpreter.cend());
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Disabling ADL, by directly specifying operators.
    // Boost✔️✔️:                 if(!(out.operator >>(result)))
    // Boost✔️✔️:                     return false;
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         };
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical.hpp:458-490
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:356-361
    // Boost✔️✔️:             template<class Alloc>
    // Boost✔️✔️:             bool operator<<(std::basic_string<CharT,Traits,Alloc> const& str) BOOST_NOEXCEPT {
    // Boost✔️✔️:                 start = str.data();
    // Boost✔️✔️:                 finish = start + str.length();
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:356-361
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:536-561
    // Boost✔️✔️:             template <typename Type>
    // Boost✔️✔️:             bool shr_signed(Type& output) {
    // Boost✔️✔️:                 if (start == finish) return false;
    // Boost✔️✔️:                 CharT const minus = lcast_char_constants<CharT>::minus;
    // Boost✔️✔️:                 CharT const plus = lcast_char_constants<CharT>::plus;
    // Boost✔️✔️:                 typedef BOOST_DEDUCED_TYPENAME make_unsigned<Type>::type utype;
    // Boost✔️✔️:                 utype out_tmp = 0;
    // Boost✔️✔️:                 bool const has_minus = Traits::eq(minus, *start);
    // Boost✔️✔️:
    // Boost✔️✔️:                 /* We won`t use `start' any more, so no need in decrementing it after */
    // Boost✔️✔️:                 if (has_minus || Traits::eq(plus, *start)) {
    // Boost✔️✔️:                     ++start;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 bool succeed = lcast_ret_unsigned<Traits, utype, CharT>(out_tmp, start, finish).convert();
    // Boost✔️✔️:                 if (has_minus) {
    // Boost✔️✔️:                     utype const comp_val = (static_cast<utype>(1) << std::numeric_limits<Type>::digits);
    // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
    // Boost✔️✔️:                     output = static_cast<Type>(0u - out_tmp);
    // Boost✔️✔️:                 } else {
    // Boost✔️✔️:                     utype const comp_val = static_cast<utype>((std::numeric_limits<Type>::max)());
    // Boost✔️✔️:                     succeed = succeed && out_tmp<=comp_val;
    // Boost✔️✔️:                     output = static_cast<Type>(out_tmp);
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 return succeed;
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:536-561
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/lcast_unsigned_converters.hpp:158-289
    // Boost✔️✔️:         template <class Traits, class T, class CharT>
    // Boost✔️✔️:         class lcast_ret_unsigned: boost::noncopyable {
    // Boost✔️✔️:             bool m_multiplier_overflowed;
    // Boost✔️✔️:             T m_multiplier;
    // Boost✔️✔️:             T& m_value;
    // Boost✔️✔️:             const CharT* const m_begin;
    // Boost✔️✔️:             const CharT* m_end;
    // Boost✔️✔️:
    // Boost✔️✔️:         public:
    // Boost✔️✔️:             lcast_ret_unsigned(T& value, const CharT* const begin, const CharT* end) BOOST_NOEXCEPT
    // Boost✔️✔️:                 : m_multiplier_overflowed(false), m_multiplier(1), m_value(value), m_begin(begin), m_end(end)
    // Boost✔️✔️:             {
    // Boost✔️✔️: #ifndef BOOST_NO_LIMITS_COMPILE_TIME_CONSTANTS
    // Boost✔️✔️:                 BOOST_STATIC_ASSERT(!std::numeric_limits<T>::is_signed);
    // Boost✔️✔️:
    // Boost✔️✔️:                 // GCC when used with flag -std=c++0x may not have std::numeric_limits
    // Boost✔️✔️:                 // specializations for __int128 and unsigned __int128 types.
    // Boost✔️✔️:                 // Try compilation with -std=gnu++0x or -std=gnu++11.
    // Boost✔️✔️:                 //
    // Boost✔️✔️:                 // http://gcc.gnu.org/bugzilla/show_bug.cgi?id=40856
    // Boost✔️✔️:                 BOOST_STATIC_ASSERT_MSG(std::numeric_limits<T>::is_specialized,
    // Boost✔️✔️:                     "std::numeric_limits are not specialized for integral type passed to boost::lexical_cast"
    // Boost✔️✔️:                 );
    // Boost✔️✔️: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             inline bool convert() {
    // Boost✔️✔️:                 CharT const czero = lcast_char_constants<CharT>::zero;
    // Boost✔️✔️:                 --m_end;
    // Boost✔️✔️:                 m_value = static_cast<T>(0);
    // Boost✔️✔️:
    // Boost✔️✔️:                 if (m_begin > m_end || *m_end < czero || *m_end >= czero + 10)
    // Boost✔️✔️:                     return false;
    // Boost✔️✔️:                 m_value = static_cast<T>(*m_end - czero);
    // Boost✔️✔️:                 --m_end;
    // Boost✔️✔️:
    // Boost✔️✔️: #ifdef BOOST_LEXICAL_CAST_ASSUME_C_LOCALE
    // Boost✔️✔️:                 return main_convert_loop();
    // Boost✔️✔️: #else
    // Boost✔️✔️:                 std::locale loc;
    // Boost✔️✔️:                 if (loc == std::locale::classic()) {
    // Boost✔️✔️:                     return main_convert_loop();
    // Boost✔️✔️:                 }
    // Boost❌❌:
    // Boost❌❌:                 typedef std::numpunct<CharT> numpunct;
    // Boost❌❌:                 numpunct const& np = BOOST_USE_FACET(numpunct, loc);
    // Boost❌❌:                 std::string const& grouping = np.grouping();
    // Boost❌❌:                 std::string::size_type const grouping_size = grouping.size();
    // Boost❌❌:
    // Boost❌❌:                 /* According to Programming languages - C++
    // Boost❌❌:                  * we MUST check for correct grouping
    // Boost❌❌:                  */
    // Boost❌❌:                 if (!grouping_size || grouping[0] <= 0) {
    // Boost❌❌:                     return main_convert_loop();
    // Boost❌❌:                 }
    // Boost❌❌:
    // Boost❌❌:                 unsigned char current_grouping = 0;
    // Boost❌❌:                 CharT const thousands_sep = np.thousands_sep();
    // Boost❌❌:                 char remained = static_cast<char>(grouping[current_grouping] - 1);
    // Boost❌❌:
    // Boost❌❌:                 for (;m_end >= m_begin; --m_end)
    // Boost❌❌:                 {
    // Boost❌❌:                     if (remained) {
    // Boost❌❌:                         if (!main_convert_iteration()) {
    // Boost❌❌:                             return false;
    // Boost❌❌:                         }
    // Boost❌❌:                         --remained;
    // Boost❌❌:                     } else {
    // Boost❌❌:                         if ( !Traits::eq(*m_end, thousands_sep) ) //|| begin == end ) return false;
    // Boost❌❌:                         {
    // Boost❌❌:                             /*
    // Boost❌❌:                              * According to Programming languages - C++
    // Boost❌❌:                              * Digit grouping is checked. That is, the positions of discarded
    // Boost❌❌:                              * separators is examined for consistency with
    // Boost❌❌:                              * use_facet<numpunct<charT> >(loc ).grouping()
    // Boost❌❌:                              *
    // Boost❌❌:                              * BUT what if there is no separators at all and grouping()
    // Boost❌❌:                              * is not empty? Well, we have no extraced separators, so we
    // Boost❌❌:                              * won`t check them for consistency. This will allow us to
    // Boost❌❌:                              * work with "C" locale from other locales
    // Boost❌❌:                              */
    // Boost❌❌:                             return main_convert_loop();
    // Boost❌❌:                         } else {
    // Boost❌❌:                             if (m_begin == m_end) return false;
    // Boost❌❌:                             if (current_grouping < grouping_size - 1) ++current_grouping;
    // Boost❌❌:                             remained = grouping[current_grouping];
    // Boost❌❌:                         }
    // Boost❌❌:                     }
    // Boost❌❌:                 } /*for*/
    // Boost❌❌:
    // Boost❌❌:                 return true;
    // Boost✔️✔️: #endif
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:         private:
    // Boost✔️✔️:             // Iteration that does not care about grouping/separators and assumes that all
    // Boost✔️✔️:             // input characters are digits
    // Boost✔️✔️:             inline bool main_convert_iteration() BOOST_NOEXCEPT {
    // Boost✔️✔️:                 CharT const czero = lcast_char_constants<CharT>::zero;
    // Boost✔️✔️:                 T const maxv = (std::numeric_limits<T>::max)();
    // Boost✔️✔️:
    // Boost✔️✔️:                 m_multiplier_overflowed = m_multiplier_overflowed || (maxv/10 < m_multiplier);
    // Boost✔️✔️:                 m_multiplier = static_cast<T>(m_multiplier * 10);
    // Boost✔️✔️:
    // Boost✔️✔️:                 T const dig_value = static_cast<T>(*m_end - czero);
    // Boost✔️✔️:                 T const new_sub_value = static_cast<T>(m_multiplier * dig_value);
    // Boost✔️✔️:
    // Boost✔️✔️:                 // We must correctly handle situations like `000000000000000000000000000001`.
    // Boost✔️✔️:                 // So we take care of overflow only if `dig_value` is not '0'.
    // Boost✔️✔️:                 if (*m_end < czero || *m_end >= czero + 10  // checking for correct digit
    // Boost✔️✔️:                     || (dig_value && (                      // checking for overflow of ...
    // Boost✔️✔️:                         m_multiplier_overflowed                             // ... multiplier
    // Boost✔️✔️:                         || static_cast<T>(maxv / dig_value) < m_multiplier  // ... subvalue
    // Boost✔️✔️:                         || static_cast<T>(maxv - new_sub_value) < m_value   // ... whole expression
    // Boost✔️✔️:                     ))
    // Boost✔️✔️:                 ) return false;
    // Boost✔️✔️:
    // Boost✔️✔️:                 m_value = static_cast<T>(m_value + new_sub_value);
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:
    // Boost✔️✔️:             bool main_convert_loop() BOOST_NOEXCEPT {
    // Boost✔️✔️:                 for ( ; m_end >= m_begin; --m_end) {
    // Boost✔️✔️:                     if (!main_convert_iteration()) {
    // Boost✔️✔️:                         return false;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return true;
    // Boost✔️✔️:             }
    // Boost✔️✔️:         };
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/lcast_unsigned_converters.hpp:158-289
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/trim.hpp:233-243
    // Boost✔️✔️:         template<typename SequenceT, typename PredicateT>
    // Boost✔️✔️:         inline void trim_right_if(SequenceT& Input, PredicateT IsSpace)
    // Boost✔️✔️:         {
    // Boost✔️✔️:             Input.erase(
    //                 ::boost::algorithm::detail::trim_end(  // Boost✔️✔️:
    //                     ::boost::begin(Input),  // Boost✔️✔️:
    //                     ::boost::end(Input),  // Boost✔️✔️:
    // Boost✔️✔️:                     IsSpace ),
    // Boost✔️✔️:                 ::boost::end(Input)
    // Boost✔️✔️:                 );
    // Boost✔️✔️:         }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/trim.hpp:233-243
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/trim.hpp:254-260
    // Boost✔️✔️:         template<typename SequenceT>
    // Boost✔️✔️:         inline void trim_right(SequenceT& Input, const std::locale& Loc=std::locale())
    // Boost✔️✔️:         {
    // Boost✔️✔️:             ::boost::algorithm::trim_right_if(
    //                 Input,  // Boost✔️✔️:
    // Boost✔️✔️:                 is_space(Loc) );
    // Boost✔️✔️:         }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/trim.hpp:254-260
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/trim.hpp:44-58
    // Boost✔️✔️:             template< typename ForwardIteratorT, typename PredicateT >
    //             inline ForwardIteratorT trim_end_iter_select(  // Boost✔️✔️:
    //                 ForwardIteratorT InBegin,  // Boost✔️✔️:
    //                 ForwardIteratorT InEnd,  // Boost✔️✔️:
    // Boost✔️✔️:                 PredicateT IsSpace,
    // Boost✔️✔️:                 std::bidirectional_iterator_tag )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 for( ForwardIteratorT It=InEnd; It!=InBegin;  )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     if ( !IsSpace(*(--It)) )
    // Boost✔️✔️:                         return ++It;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return InBegin;
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/trim.hpp:44-58
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/classification.hpp:41-42
    // Boost✔️✔️:                 is_classifiedF(std::ctype_base::mask Type, std::locale const & Loc = std::locale()) :
    // Boost✔️✔️:                     m_Type(Type), m_Locale(Loc) {}
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/classification.hpp:41-42
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/classification.hpp:44-48
    // Boost✔️✔️:                 template<typename CharT>
    // Boost✔️✔️:                 bool operator()( CharT Ch ) const
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     return std::use_facet< std::ctype<CharT> >(m_Locale).is( m_Type, Ch );
    // Boost✔️✔️:                 }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/classification.hpp:44-48
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:47-53
    // RDKit✔️🔝: template <>
    // RDKit✔️🔝: inline std::string rdvalue_cast<std::string>(RDValue_cast_t v) {
    // RDKit✔️🔝:   if (rdvalue_is<std::string>(v)) {
    // RDKit✔️🔝:     return *v.ptrCast<std::string>();
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   throw std::bad_any_cast();
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/RDValue.h:47-53
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/classification.hpp:56-60
    //         inline detail::is_classifiedF  // Boost✔️✔️:
    // Boost✔️✔️:         is_space(const std::locale& Loc=std::locale())
    // Boost✔️✔️:         {
    // Boost✔️✔️:             return detail::is_classifiedF(std::ctype_base::space, Loc);
    // Boost✔️✔️:         }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/classification.hpp:56-60
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/trim.hpp:77-87
    // Boost✔️✔️:             template< typename ForwardIteratorT, typename PredicateT >
    //             inline ForwardIteratorT trim_end(  // Boost✔️✔️:
    //                 ForwardIteratorT InBegin,  // Boost✔️✔️:
    //                 ForwardIteratorT InEnd,  // Boost✔️✔️:
    // Boost✔️✔️:                 PredicateT IsSpace )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️✔️:                     std::iterator_traits<ForwardIteratorT>::iterator_category category;
    // Boost✔️✔️:
    // Boost✔️✔️:                 return ::boost::algorithm::detail::trim_end_iter_select( InBegin, InEnd, IsSpace, category() );
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/algorithm/string/detail/trim.hpp:77-87
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/lcast_char_constants.hpp:30-40
    // Boost✔️✔️:         template < typename Char >
    // Boost✔️✔️:         struct lcast_char_constants {
    // Boost✔️✔️:             // We check in tests assumption that static casted character is
    // Boost✔️✔️:             // equal to correctly written C++ literal: U'0' == static_cast<char32_t>('0')
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, zero  = static_cast<Char>('0'));
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, minus = static_cast<Char>('-'));
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, plus = static_cast<Char>('+'));
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, lowercase_e = static_cast<Char>('e'));
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, capital_e = static_cast<Char>('E'));
    // Boost✔️✔️:             BOOST_STATIC_CONSTANT(Char, c_decimal_separator = static_cast<Char>('.'));
    // Boost✔️✔️:         };
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/lcast_char_constants.hpp:30-40
    // The called String/int specialization borrows the exact input range;
    // these complete accessors/constructor/dispatcher close its source path.
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:169-171
    // Boost✔️✔️:             const CharT* cbegin() const BOOST_NOEXCEPT {
    // Boost✔️✔️:                 return start;
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:169-171
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:173-175
    // Boost✔️✔️:             const CharT* cend() const BOOST_NOEXCEPT {
    // Boost✔️✔️:                 return finish;
    // Boost✔️✔️:             }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:173-175
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:507-511
    // Boost✔️✔️:         public:
    // Boost✔️✔️:             lexical_ostream_limited_src(const CharT* begin, const CharT* end) BOOST_NOEXCEPT
    // Boost✔️✔️:               : start(begin)
    // Boost✔️✔️:               , finish(end)
    // Boost✔️✔️:             {}
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:507-511
    // BEGIN BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:643-643
    // Boost✔️✔️:             bool operator>>(int& output)                        { return shr_signed(output); }
    // END BOOST CPP FUNCTION: target/agent-handoff/Q01/B1/scalar_dependency/boost_1_81_sources/boost/lexical_cast/detail/converter_lexical_streams.hpp:643-643
    // Reviewed UInt source closure: numeric bounds are checked before
    // the low-level cast. O(1), no allocation; errors retain input context.
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/cast.hpp:38-54
    //     template <typename Target, typename Source>  // Boost❗✔️:
    // Boost❗✔️:     inline Target numeric_cast( Source arg )
    // Boost❗✔️:     {
    // Boost❗✔️:         typedef numeric::conversion_traits<Target, Source>   conv_traits;
    // Boost❗✔️:         typedef numeric::numeric_cast_traits<Target, Source> cast_traits;
    // Boost❗✔️:         typedef boost::numeric::converter
    // Boost❗✔️:             <
    // Boost❗✔️:                 Target,
    //                 Source,  // Boost❗✔️:
    // Boost❗✔️:                 conv_traits,
    //                 typename cast_traits::overflow_policy,  // Boost❗✔️:
    //                 typename cast_traits::rounding_policy,  // Boost❗✔️:
    // Boost❗✔️:                 boost::numeric::raw_converter< conv_traits >,
    // Boost❗✔️:                 typename cast_traits::range_checking_policy
    // Boost❗✔️:             > converter;
    // Boost❗✔️:         return converter::convert(arg);
    // Boost❗✔️:     }
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/cast.hpp:38-54
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:340-346
    // Boost❗✔️:   template<class Traits, class OverflowHandler>
    // Boost❗✔️:   struct GetRC_Unsig2Sig
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef GT_HiT<Traits> Pred1 ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef generic_range_checker<Traits,non_applicable,Pred1,OverflowHandler> type ;
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:340-346
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:157-169
    // Boost❗✔️:     template<class Traits>
    // Boost❗✔️:     struct GT_HiT : applicable
    // Boost❗✔️:     {
    // Boost❗✔️:       typedef typename Traits::target_type T ;
    // Boost❗✔️:       typedef typename Traits::source_type S ;
    // Boost❗✔️:       typedef typename Traits::argument_type argument_type ;
    // Boost❗✔️:
    // Boost❗✔️:       static range_check_result apply ( argument_type s )
    // Boost❗✔️:       {
    // Boost❗✔️:         return s > static_cast<S>(bounds<T>::highest())
    // Boost❗✔️:                  ? cPosOverflow : cInRange ;
    // Boost❗✔️:       }
    // Boost❗✔️:     } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:157-169
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:279-295
    // Boost❗✔️:   template<class Traits, class IsNegOverflow, class IsPosOverflow, class OverflowHandler>
    // Boost❗✔️:   struct generic_range_checker
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef OverflowHandler overflow_handler ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::argument_type argument_type ;
    // Boost❗✔️:
    // Boost❗✔️:     static range_check_result out_of_range ( argument_type s )
    // Boost❗✔️:     {
    // Boost❗✔️:       typedef typename combine<IsNegOverflow,IsPosOverflow>::type Predicate ;
    // Boost❗✔️:
    // Boost❗✔️:       return Predicate::apply(s);
    // Boost❗✔️:     }
    // Boost❗✔️:
    // Boost❗✔️:     static void validate_range ( argument_type s )
    // Boost❗✔️:       { OverflowHandler()( out_of_range(s) ) ; }
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:279-295
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:497-517
    // Boost❗✔️:   template<class Traits,class RangeChecker,class RawConverter>
    // Boost❗✔️:   struct non_rounding_converter : public RangeChecker
    // Boost❗✔️:                                  ,public RawConverter
    // Boost❗✔️:   {
    // Boost❗✔️:     typedef RangeChecker RangeCheckerBase ;
    // Boost❗✔️:     typedef RawConverter RawConverterBase ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef Traits traits ;
    // Boost❗✔️:
    // Boost❗✔️:     typedef typename Traits::source_type   source_type   ;
    // Boost❗✔️:     typedef typename Traits::argument_type argument_type ;
    // Boost❗✔️:     typedef typename Traits::result_type   result_type   ;
    // Boost❗✔️:
    // Boost❗✔️:     static source_type nearbyint ( argument_type s ) { return s ; }
    // Boost❗✔️:
    // Boost❗✔️:     static result_type convert ( argument_type s )
    // Boost❗✔️:     {
    // Boost❗✔️:       RangeCheckerBase::validate_range(s);
    // Boost❗✔️:       return RawConverterBase::low_level_convert(s);
    // Boost❗✔️:     }
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/converter.hpp:497-517
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:150-156
    // Boost❗✔️: class positive_overflow : public bad_numeric_cast
    // Boost❗✔️: {
    // Boost❗✔️:   public:
    // Boost❗✔️:
    // Boost❗✔️:     const char * what() const BOOST_NOEXCEPT_OR_NOTHROW BOOST_OVERRIDE
    // Boost❗✔️:       { return "bad numeric conversion: positive overflow"; }
    // Boost❗✔️: };
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:150-156
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:158-174
    // Boost❗✔️: struct def_overflow_handler
    // Boost❗✔️: {
    // Boost❗✔️:   void operator() ( range_check_result r ) // throw(negative_overflow,positive_overflow)
    // Boost❗✔️:   {
    // Boost❗✔️: #ifndef BOOST_NO_EXCEPTIONS
    // Boost❗✔️:     if ( r == cNegOverflow )
    // Boost❗✔️:       throw negative_overflow() ;
    // Boost❗✔️:     else if ( r == cPosOverflow )
    // Boost❗✔️:            throw positive_overflow() ;
    // Boost❗✔️: #else
    // Boost❗✔️:     if ( r == cNegOverflow )
    // Boost❗✔️:       ::boost::throw_exception(negative_overflow()) ;
    // Boost❗✔️:     else if ( r == cPosOverflow )
    // Boost❗✔️:            ::boost::throw_exception(positive_overflow()) ;
    // Boost❗✔️: #endif
    // Boost❗✔️:   }
    // Boost❗✔️: } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/converter_policies.hpp:158-174
    // BEGIN BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/bounds.hpp:19-29
    // Boost❗✔️:   template<class N>
    // Boost❗✔️:   class Integral
    // Boost❗✔️:   {
    // Boost❗✔️:       typedef std::numeric_limits<N> limits ;
    // Boost❗✔️:
    // Boost❗✔️:     public :
    //      // Boost❗✔️:
    // Boost❗✔️:       static N lowest  () { return limits::min BOOST_PREVENT_MACRO_SUBSTITUTION (); }
    // Boost❗✔️:       static N highest () { return limits::max BOOST_PREVENT_MACRO_SUBSTITUTION (); }
    // Boost❗✔️:       static N smallest() { return static_cast<N>(1); }
    // Boost❗✔️:   } ;
    // END BOOST COMPLETE PROPOSED CPP FUNCTION: target/agent-handoff/Q01/B1/uint_complete_preparation_v1/official_source/include/boost/numeric/conversion/detail/bounds.hpp:19-29
    let property = "_CanonicalRankingNumber";
    let Some(value) = atom.prop(property) else {
        return Ok(0);
    };
    crate::property_value_to_int(value).map_err(|error| match error {
        crate::PropertyIntReadError::UnsignedOverflow { value } => {
            CanonicalRankError::UnsignedRankOverflow {
                atom_index,
                property,
                value,
            }
        }
        crate::PropertyIntReadError::InvalidKind { kind } => {
            CanonicalRankError::InvalidPropertyKind {
                atom_index,
                property,
                kind,
            }
        }
        crate::PropertyIntReadError::Lexical { .. }
        | crate::PropertyIntReadError::SignedTextOverflow { .. } => {
            CanonicalRankError::InvalidPropertyKind {
                atom_index,
                property,
                kind: cosmolkit_model::PropertyValueKind::String,
            }
        }
    })
}

fn compare_canon_atom_base_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
    flags: CanonRankFlags,
) -> Result<Ordering, CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION AtomCompareFunctor::basecomp
    // RDKit✔️✔️:   ivi = dp_atoms[i].index;
    // RDKit✔️✔️:   ivj = dp_atoms[j].index;
    // RDKit✔️✔️:   if (df_useNonStereoRanks) {
    // RDKit✔️✔️:     int rankingNumber_i = 0;
    // RDKit✔️✔️:     int rankingNumber_j = 0;
    // RDKit✔️✔️:     dp_atoms[i].atom->getPropIfPresent(
    // RDKit✔️✔️:         common_properties::_CanonicalRankingNumber, rankingNumber_i);
    // RDKit✔️✔️:     dp_atoms[j].atom->getPropIfPresent(
    // RDKit✔️✔️:         common_properties::_CanonicalRankingNumber, rankingNumber_j);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useAtomMaps || df_useAtomMapsOnDummies) {
    // RDKit✔️✔️:     int molAtomMapNumber_i = 0;
    // RDKit✔️✔️:     int molAtomMapNumber_j = 0;
    // RDKit✔️✔️:     if (df_useAtomMaps ||
    // RDKit✔️✔️:         (df_useAtomMapsOnDummies && dp_atoms[i].atom->getAtomicNum() == 0)) {
    // RDKit✔️✔️:       dp_atoms[i].atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit✔️✔️:                                          molAtomMapNumber_i);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (df_useAtomMaps ||
    // RDKit✔️✔️:         (df_useAtomMapsOnDummies && dp_atoms[j].atom->getAtomicNum() == 0)) {
    // RDKit✔️✔️:       dp_atoms[j].atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit✔️✔️:                                          molAtomMapNumber_j);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].degree;
    // RDKit✔️✔️:   ivj = dp_atoms[j].degree;
    // RDKit✔️✔️:   if (dp_atoms[i].p_symbol && dp_atoms[j].p_symbol) {
    // RDKit✔️✔️:     if (*(dp_atoms[i].p_symbol) < *(dp_atoms[j].p_symbol)) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (*(dp_atoms[i].p_symbol) > *(dp_atoms[j].p_symbol)) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].atom->getAtomicNum();
    // RDKit✔️✔️:   ivj = dp_atoms[j].atom->getAtomicNum();
    // RDKit✔️✔️:   if (df_useIsotopes) {
    // RDKit✔️✔️:     ivi = dp_atoms[i].atom->getIsotope();
    // RDKit✔️✔️:     ivj = dp_atoms[j].atom->getIsotope();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ivi = dp_atoms[i].totalNumHs;
    // RDKit✔️✔️:   ivj = dp_atoms[j].totalNumHs;
    // RDKit✔️✔️:   ivi = dp_atoms[i].atom->getFormalCharge();
    // RDKit✔️✔️:   ivj = dp_atoms[j].atom->getFormalCharge();
    // RDKit✔️✔️:   if (df_useChiralPresence) {
    // RDKit✔️✔️:     ivi =
    // RDKit✔️✔️:         dp_atoms[i].atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED;
    // RDKit✔️✔️:     ivj =
    // RDKit✔️✔️:         dp_atoms[j].atom->getChiralTag() != Atom::ChiralType::CHI_UNSPECIFIED;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useChirality) {
    // RDKit✔️✔️:     ivi = dp_atoms[i].whichStereoGroup;
    // RDKit✔️✔️:     ivj = dp_atoms[j].whichStereoGroup;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (df_useChiralityRings) {
    // RDKit✔️✔️:     ivi = getAtomRingNbrCode(i);
    // RDKit✔️✔️:     ivj = getAtomRingNbrCode(j);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION AtomCompareFunctor::basecomp
    let left_atom = &atoms[left];
    let right_atom = &atoms[right];
    let mut cmp = left_atom.index.cmp(&right_atom.index);
    if cmp != Ordering::Equal {
        return Ok(cmp);
    }

    if flags.use_non_stereo_ranks {
        // Same source guards and getter order, with O(log P) borrowed map
        // lookups instead of O(P) Dict scans. String conversion keeps O(B)
        // time and avoids the source String copy and locale object; tagged
        // Int extraction remains O(1). The fixed classical locale preserves
        // every modeled cast result while reducing getter storage/work.
        // RDKit✔️🔝: dp_atoms[i].atom->getPropIfPresent(
        // RDKit✔️🔝:     common_properties::_CanonicalRankingNumber, rankingNumber_i);
        // RDKit✔️🔝: dp_atoms[j].atom->getPropIfPresent(
        // RDKit✔️🔝:     common_properties::_CanonicalRankingNumber, rankingNumber_j);
        // Current class and source flag precede both lookups. The left getter
        // can throw before the right; numeric comparison follows both getters.
        let left_rank = canonical_rank_property_to_int(left_atom.source_atom, left)?;
        let right_rank = canonical_rank_property_to_int(right_atom.source_atom, right)?;
        cmp = left_rank.cmp(&right_rank);
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }

    if flags.use_atom_maps || flags.use_atom_maps_on_dummies {
        let left_map = if flags.use_atom_maps
            || (flags.use_atom_maps_on_dummies && left_atom.atomic_number == 0)
        {
            left_atom.atom_map
        } else {
            0
        };
        let right_map = if flags.use_atom_maps
            || (flags.use_atom_maps_on_dummies && right_atom.atomic_number == 0)
        {
            right_atom.atom_map
        } else {
            0
        };
        cmp = left_map.cmp(&right_map);
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }

    cmp = left_atom.degree.cmp(&right_atom.degree);
    if cmp != Ordering::Equal {
        return Ok(cmp);
    }

    if let (Some(left_symbol), Some(right_symbol)) = (left_atom.p_symbol, right_atom.p_symbol) {
        return Ok(left_symbol.cmp(right_symbol));
    }

    cmp = left_atom.atomic_number.cmp(&right_atom.atomic_number);
    if cmp != Ordering::Equal {
        return Ok(cmp);
    }

    if flags.use_isotopes {
        cmp = left_atom.isotope.cmp(&right_atom.isotope);
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }

    cmp = left_atom.total_num_hs.cmp(&right_atom.total_num_hs);
    if cmp != Ordering::Equal {
        return Ok(cmp);
    }

    // RDKit basecomp stores comparison temporaries as unsigned int.
    // Preserve that C++ conversion behavior here so negative formal charges
    // wrap and order after non-negative charges (e.g. -1 > 0 in this step).
    let left_charge = left_atom.formal_charge as i32 as u32;
    let right_charge = right_atom.formal_charge as i32 as u32;
    cmp = left_charge.cmp(&right_charge);
    if cmp != Ordering::Equal {
        return Ok(cmp);
    }

    if flags.use_chiral_presence {
        cmp = chiral_presence_for_kekulize(left_atom.chiral_tag)
            .cmp(&chiral_presence_for_kekulize(right_atom.chiral_tag));
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }
    if flags.use_chirality {
        // RDKit✔️✔️:     ivi = dp_atoms[i].whichStereoGroup;
        // RDKit✔️✔️:     ivj = dp_atoms[j].whichStereoGroup;
        cmp = compare_stereo_group_state_for_kekulize(atoms, left, right);
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }
    if flags.use_chirality_rings {
        // RDKit✔️✔️:     ivi = getAtomRingNbrCode(i);
        // RDKit✔️✔️:     ivj = getAtomRingNbrCode(j);
        cmp = get_atom_ring_nbr_code_for_kekulize(atoms, left)
            .cmp(&get_atom_ring_nbr_code_for_kekulize(atoms, right));
        if cmp != Ordering::Equal {
            return Ok(cmp);
        }
    }
    Ok(Ordering::Equal)
}

fn atom_pair_has_any_in_play_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
) -> bool {
    atoms[left].is_in_play || atoms[right].is_in_play
}

fn compare_stereo_group_state_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
) -> Ordering {
    let left_group = atoms[left].which_stereo_group;
    let right_group = atoms[right].which_stereo_group;
    match (left_group, right_group) {
        (0, 0) => Ordering::Equal,
        (_, 0) => Ordering::Greater,
        (0, _) => Ordering::Less,
        _ => atoms[left]
            .type_of_stereo_group
            .cmp(&atoms[right].type_of_stereo_group)
            .then_with(|| {
                if left_group == right_group {
                    if atoms[left].type_of_stereo_group == CanonStereoGroupType::Absolute {
                        get_chiral_rank_for_kekulize(atoms, left)
                            .cmp(&get_chiral_rank_for_kekulize(atoms, right))
                    } else {
                        Ordering::Equal
                    }
                } else {
                    stereo_group_symmetry_set_for_kekulize(atoms, left_group)
                        .cmp(&stereo_group_symmetry_set_for_kekulize(atoms, right_group))
                }
            }),
    }
    .then_with(|| {
        if left_group == 0 && right_group == 0 {
            chiral_presence_for_kekulize(atoms[left].chiral_tag)
                .cmp(&chiral_presence_for_kekulize(atoms[right].chiral_tag))
                .then_with(|| {
                    if chiral_presence_for_kekulize(atoms[left].chiral_tag)
                        && chiral_presence_for_kekulize(atoms[right].chiral_tag)
                    {
                        get_chiral_rank_for_kekulize(atoms, left)
                            .cmp(&get_chiral_rank_for_kekulize(atoms, right))
                    } else {
                        Ordering::Equal
                    }
                })
        } else {
            Ordering::Equal
        }
    })
}

fn chiral_presence_for_kekulize(chiral_tag: ChiralTag) -> bool {
    chiral_tag != ChiralTag::Unspecified
}

fn get_chiral_rank_for_kekulize(atoms: &[CanonAtom<'_>], atom_idx: usize) -> u32 {
    // BEGIN RDKIT CPP FUNCTION getChiralRank
    // RDKit✔️✔️: unsigned int getChiralRank(const ROMol *dp_mol, canon_atom *dp_atoms,
    // RDKit✔️✔️:                            unsigned int i) {
    // RDKit✔️✔️:   unsigned int res = 0;
    // RDKit✔️✔️:   std::vector<unsigned int> perm;
    // RDKit✔️✔️:   perm.reserve(dp_atoms[i].atom->getDegree());
    let mut res = 0u32;
    let mut perm = Vec::<i32>::with_capacity(atoms[atom_idx].all_nbr_ids.len());
    // RDKit✔️✔️:   for (const auto nbr : dp_mol->atomNeighbors(dp_atoms[i].atom)) {
    // RDKit✔️✔️:     auto rnk = dp_atoms[nbr->getIdx()].index;
    // RDKit✔️✔️:     if (std::find(perm.begin(), perm.end(), rnk) != perm.end()) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       perm.push_back(rnk);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for neighbor in atoms[atom_idx].all_nbr_ids.iter().copied() {
        let rank = atoms[neighbor].index;
        if perm.contains(&rank) {
            break;
        }
        perm.push(rank);
    }
    // RDKit✔️✔️:   if (perm.size() == dp_atoms[i].atom->getDegree()) {
    if perm.len() == atoms[atom_idx].all_nbr_ids.len() {
        // RDKit✔️✔️:     auto ctag = dp_atoms[i].atom->getChiralTag();
        // RDKit✔️✔️:     if (ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CW ||
        // RDKit✔️✔️:         ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CCW) {
        if matches!(
            atoms[atom_idx].chiral_tag,
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) {
            // RDKit✔️✔️:       auto sortedPerm = perm;
            // RDKit✔️✔️:       std::sort(sortedPerm.begin(), sortedPerm.end());
            // RDKit✔️✔️:       auto nswaps = countSwapsToInterconvert(perm, sortedPerm);
            let mut sorted_perm = perm.clone();
            sorted_perm.sort_unstable();
            let swaps = count_swaps_to_interconvert(&perm, sorted_perm);
            // RDKit✔️✔️:       res = ctag == Atom::ChiralType::CHI_TETRAHEDRAL_CW ? 2 : 1;
            res = if atoms[atom_idx].chiral_tag == ChiralTag::TetrahedralCw {
                2
            } else {
                1
            };
            // RDKit✔️✔️:       if (nswaps % 2) {
            // RDKit✔️✔️:         res = res == 2 ? 1 : 2;
            // RDKit✔️✔️:       }
            if swaps % 2 == 1 {
                res = if res == 2 { 1 } else { 2 };
            }
        }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION getChiralRank
    res
}

fn stereo_group_symmetry_set_for_kekulize(
    atoms: &[CanonAtom<'_>],
    group_idx: usize,
) -> BTreeSet<i32> {
    atoms
        .iter()
        .filter(|atom| atom.which_stereo_group == group_idx)
        .map(|atom| atom.index)
        .collect()
}

fn get_atom_ring_nbr_code_for_kekulize(atoms: &[CanonAtom<'_>], atom_idx: usize) -> i32 {
    // BEGIN RDKIT CPP FUNCTION AtomCompareFunctor::getAtomRingNbrCode
    // RDKit✔️✔️:   unsigned int getAtomRingNbrCode(unsigned int i) const {
    // RDKit✔️✔️:     if (!dp_atoms[i].hasRingNbr) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     }
    if !atoms[atom_idx].has_ring_nbr {
        return 0;
    }
    // RDKit✔️✔️:     auto nbrs = dp_atoms[i].nbrIds.get();
    // RDKit✔️✔️:     unsigned int code = 0;
    // RDKit✔️✔️:     for (unsigned j = 0; j < dp_atoms[i].degree; ++j) {
    // RDKit✔️✔️:       if (dp_atoms[nbrs[j]].isRingStereoAtom) {
    // RDKit✔️✔️:         code += dp_atoms[nbrs[j]].index * 10000 + 1;  // j;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return code;
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION AtomCompareFunctor::getAtomRingNbrCode
    atoms[atom_idx]
        .nbr_ids
        .iter()
        .filter(|&&neighbor| atoms[neighbor].is_ring_stereo_atom)
        .map(|&neighbor| {
            atoms[neighbor]
                .index
                .saturating_mul(10_000)
                .saturating_add(1)
        })
        .sum()
}

fn special_symmetry_rank_refinement_required(
    view: &CanonRankReadView<'_>,
    order: &[usize],
    count: &[usize],
) -> bool {
    let mut ties = false;
    let mut sym_ring_atoms = 0usize;
    let mut ring_atoms = 0usize;
    let mut branching_ring_atom = false;
    for &atom_idx in order.iter().take(view.num_atoms()) {
        let atom = AtomId::new(atom_idx);
        if view.rings.num_atom_rings(atom) > 0 {
            if count[atom_idx] > 2 {
                sym_ring_atoms += count[atom_idx];
            }
            ring_atoms += 1;
            if view.rings.num_atom_rings(atom) > 1 && count[atom_idx] > 1 {
                branching_ring_atom = true;
            }
        }
        if count[atom_idx] == 0 {
            ties = true;
        }
    }
    ties && ring_atoms > 0
        && (sym_ring_atoms as f32) / (ring_atoms as f32) > 0.5
        && branching_ring_atom
}

fn compare_ring_atoms_concerning_num_neighbors_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
) {
    // BEGIN RDKIT CPP FUNCTION compareRingAtomsConcerningNumNeighbors
    // RDKit✔️✔️: void compareRingAtomsConcerningNumNeighbors(Canon::canon_atom *atoms,
    // RDKit✔️✔️:                                             unsigned int nAtoms,
    // RDKit✔️✔️:                                             const ROMol &mol) {
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   RingInfo *ringInfo = mol.getRingInfo();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   std::vector<char> visited(nAtoms);
    // RDKit✔️✔️:   std::vector<char> lastLevelNbrs(nAtoms);
    // RDKit✔️✔️:   std::vector<char> currentLevelNbrs(nAtoms);
    // RDKit✔️✔️:   std::vector<int> revisitedNeighbors(nAtoms);
    let mut visited = vec![false; n_atoms];
    let mut last_level_nbrs = vec![false; n_atoms];
    let mut current_level_nbrs = vec![false; n_atoms];
    let mut revisited_neighbors = vec![0i32; n_atoms];
    // RDKit✔️✔️:   std::vector<int> visitedIndices;
    // RDKit✔️✔️:   std::vector<int> lastLevelIndices;
    // RDKit✔️✔️:   std::vector<int> modifiedRevisited;
    let mut visited_indices = Vec::<usize>::new();
    let mut last_level_indices = Vec::<usize>::new();
    let mut modified_revisited = Vec::<usize>::new();
    // RDKit✔️✔️:   for (unsigned idx = 0; idx < nAtoms; ++idx) {
    // RDKit✔️✔️:     const Canon::canon_atom &a = atoms[idx];
    // RDKit✔️✔️:     if (!ringInfo->isInitialized() ||
    // RDKit✔️✔️:         ringInfo->numAtomRings(a.atom->getIdx()) < 1) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for idx in 0..n_atoms {
        if view.rings.num_atom_rings(AtomId::new(idx)) < 1 {
            continue;
        }
        // RDKit✔️✔️:     std::deque<int> neighbors;
        // RDKit✔️✔️:     neighbors.push_back(idx);
        // RDKit✔️✔️:     unsigned currentRNIdx = 0;
        // RDKit✔️✔️:     atoms[idx].neighborNum.reserve(1000);
        // RDKit✔️✔️:     atoms[idx].revistedNeighbors.assign(1000, 0);
        let mut neighbors = VecDeque::new();
        neighbors.push_back(idx);
        let mut current_rn_idx = 0usize;
        atoms[idx].neighbor_num.clear();
        atoms[idx].neighbor_num.reserve(1000);
        atoms[idx].revisted_neighbors.clear();
        atoms[idx].revisted_neighbors.resize(1000, 0);
        // RDKit✔️✔️:     for (int i : visitedIndices) {
        // RDKit✔️✔️:       visited[i] = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     visitedIndices.clear();
        for i in visited_indices.drain(..) {
            visited[i] = false;
        }
        // RDKit✔️✔️:     for (int i : lastLevelIndices) {
        // RDKit✔️✔️:       lastLevelNbrs[i] = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     lastLevelIndices.clear();
        for i in last_level_indices.drain(..) {
            last_level_nbrs[i] = false;
        }
        // RDKit✔️✔️:     std::vector<int> nextLevelNbrs;
        let mut next_level_nbrs = Vec::<usize>::new();
        // RDKit✔️✔️:     while (!neighbors.empty()) {
        while !neighbors.is_empty() {
            // RDKit✔️✔️:       unsigned int numLevelNbrs = 0;
            // RDKit✔️✔️:       nextLevelNbrs.resize(0);
            let mut num_level_nbrs = 0i32;
            next_level_nbrs.clear();
            // RDKit✔️✔️:       while (!neighbors.empty()) {
            while let Some(nidx) = neighbors.pop_front() {
                // RDKit✔️✔️:         int nidx = neighbors.front();
                // RDKit✔️✔️:         neighbors.pop_front();
                // RDKit✔️✔️:         const Canon::canon_atom &atom = atoms[nidx];
                // RDKit✔️✔️:         if (!ringInfo->isInitialized() ||
                // RDKit✔️✔️:             ringInfo->numAtomRings(atom.atom->getIdx()) < 1) {
                // RDKit✔️✔️:           continue;
                // RDKit✔️✔️:         }
                // symmetrize_sssr() always returns initialized RingInfo here,
                // so only the ring-membership branch remains.
                if view.rings.num_atom_rings(AtomId::new(nidx)) < 1 {
                    continue;
                }
                // RDKit✔️✔️:         lastLevelNbrs[nidx] = 1;
                // RDKit✔️✔️:         lastLevelIndices.push_back(nidx);
                // RDKit✔️✔️:         visited[nidx] = 1;
                // RDKit✔️✔️:         visitedIndices.push_back(nidx);
                last_level_nbrs[nidx] = true;
                last_level_indices.push(nidx);
                visited[nidx] = true;
                visited_indices.push(nidx);
                // RDKit✔️✔️:         for (unsigned int j = 0; j < atom.degree; j++) {
                // RDKit✔️✔️:           int iidx = atom.nbrIds[j];
                // RDKit✔️✔️:           if (!visited[iidx]) {
                // RDKit✔️✔️:             currentLevelNbrs[iidx] = 1;
                // RDKit✔️✔️:             numLevelNbrs++;
                // RDKit✔️✔️:             visited[iidx] = 1;
                // RDKit✔️✔️:             visitedIndices.push_back(iidx);
                // RDKit✔️✔️:             nextLevelNbrs.push_back(iidx);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                for iidx in atoms[nidx].nbr_ids.iter().copied() {
                    if !visited[iidx] {
                        current_level_nbrs[iidx] = true;
                        num_level_nbrs += 1;
                        visited[iidx] = true;
                        visited_indices.push(iidx);
                        next_level_nbrs.push(iidx);
                    }
                }
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         const Canon::canon_atom &natom = atoms[i];
            // RDKit✔️✔️:         for (unsigned int k = 0; k < natom.degree; k++) {
            // RDKit✔️✔️:           int jidx = natom.nbrIds[k];
            // RDKit✔️✔️:           if (currentLevelNbrs[jidx] || lastLevelNbrs[jidx]) {
            // RDKit✔️✔️:             if (revisitedNeighbors[jidx] == 0) {
            // RDKit✔️✔️:               modifiedRevisited.push_back(jidx);
            // RDKit✔️✔️:             }
            // RDKit✔️✔️:             revisitedNeighbors[jidx] += 1;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                for jidx in atoms[i].nbr_ids.iter().copied() {
                    if current_level_nbrs[jidx] || last_level_nbrs[jidx] {
                        if revisited_neighbors[jidx] == 0 {
                            modified_revisited.push(jidx);
                        }
                        revisited_neighbors[jidx] += 1;
                    }
                }
            }
            // RDKit✔️✔️:       for (int i : lastLevelIndices) {
            // RDKit✔️✔️:         lastLevelNbrs[i] = 0;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       lastLevelIndices.clear();
            for i in last_level_indices.drain(..) {
                last_level_nbrs[i] = false;
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         lastLevelNbrs[i] = 1;
            // RDKit✔️✔️:         lastLevelIndices.push_back(i);
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                last_level_nbrs[i] = true;
                last_level_indices.push(i);
            }
            // RDKit✔️✔️:       for (int i : nextLevelNbrs) {
            // RDKit✔️✔️:         currentLevelNbrs[i] = 0;
            // RDKit✔️✔️:       }
            for i in next_level_nbrs.iter().copied() {
                current_level_nbrs[i] = false;
            }
            // RDKit✔️✔️:       std::vector<int> tmp;
            // RDKit✔️✔️:       tmp.reserve(30);
            // RDKit✔️✔️:       for (int i : modifiedRevisited) {
            // RDKit✔️✔️:         tmp.push_back(revisitedNeighbors[i]);
            // RDKit✔️✔️:       }
            let mut tmp = Vec::with_capacity(30);
            for i in modified_revisited.iter().copied() {
                tmp.push(revisited_neighbors[i]);
            }
            // RDKit✔️✔️:       std::sort(tmp.begin(), tmp.end());
            // RDKit✔️✔️:       tmp.push_back(-1);
            tmp.sort_unstable();
            tmp.push(-1);
            // RDKit✔️✔️:       for (int i : tmp) {
            // RDKit✔️✔️:         if (currentRNIdx >= atoms[idx].revistedNeighbors.size()) {
            // RDKit✔️✔️:           atoms[idx].revistedNeighbors.resize(
            // RDKit✔️✔️:               atoms[idx].revistedNeighbors.size() + 1000);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         atoms[idx].revistedNeighbors[currentRNIdx] = i;
            // RDKit✔️✔️:         currentRNIdx++;
            // RDKit✔️✔️:       }
            for value in tmp {
                if current_rn_idx >= atoms[idx].revisted_neighbors.len() {
                    atoms[idx]
                        .revisted_neighbors
                        .resize(atoms[idx].revisted_neighbors.len() + 1000, 0);
                }
                atoms[idx].revisted_neighbors[current_rn_idx] = value;
                current_rn_idx += 1;
            }
            // RDKit✔️✔️:       for (int i : modifiedRevisited) {
            // RDKit✔️✔️:         revisitedNeighbors[i] = 0;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       modifiedRevisited.clear();
            for i in modified_revisited.drain(..) {
                revisited_neighbors[i] = 0;
            }
            // RDKit✔️✔️:       atoms[idx].neighborNum.push_back(numLevelNbrs);
            // RDKit✔️✔️:       atoms[idx].neighborNum.push_back(-1);
            atoms[idx].neighbor_num.push(num_level_nbrs);
            atoms[idx].neighbor_num.push(-1);
            // RDKit✔️✔️:       neighbors.insert(neighbors.end(), nextLevelNbrs.begin(),
            // RDKit✔️✔️:                        nextLevelNbrs.end());
            // RDKit✔️✔️:     }
            neighbors.extend(next_level_nbrs.iter().copied());
        }
        // RDKit✔️✔️:     atoms[idx].revistedNeighbors.resize(currentRNIdx);
        atoms[idx].revisted_neighbors.truncate(current_rn_idx);
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION compareRingAtomsConcerningNumNeighbors
}

fn update_atom_neighbor_index_for_kekulize(atoms: &mut [CanonAtom<'_>], atom_idx: usize) {
    // BEGIN COMPLETE RDKit .6 CHEM32 updateAtomNeighborIndex
    // RDKit❗❌: void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs) {
    // RDKit❗❌:   PRECONDITION(atoms, "bad pointer");
    // RDKit❗❌:   for (auto &nbr : nbrs) {
    // RDKit❗❌:     unsigned nbrIdx = nbr.nbrIdx;
    // RDKit❗❌:     unsigned newSymClass = atoms[nbrIdx].index;
    // RDKit❗❌:     nbr.nbrSymClass = newSymClass;
    // RDKit❗❌:   }
    // RDKit❗❌:   // Neighbor lists are normally very short and already close to sorted after
    // RDKit❗❌:   // partition refinement. Insertion sort avoids std::sort's setup overhead and
    // RDKit❗❌:   // minimizes movement in that common case.
    // RDKit❗❌:   for (size_t i = 1; i < nbrs.size(); ++i) {
    // RDKit❗❌:     if (!bondholder::greater(nbrs[i], nbrs[i - 1])) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:     auto value = std::move(nbrs[i]);
    // RDKit❗❌:     size_t j = i;
    // RDKit❗❌:     do {
    // RDKit❗❌:       nbrs[j] = std::move(nbrs[j - 1]);
    // RDKit❗❌:       --j;
    // RDKit❗❌:     } while (j && bondholder::greater(value, nbrs[j - 1]));
    // RDKit❗❌:     nbrs[j] = std::move(value);
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END COMPLETE RDKit .6 CHEM32 updateAtomNeighborIndex
    // Only this actual repeated-update owner changes sort. Existing update Vec
    // and whole-atom rank snapshot costs remain qualified (❌); the shared
    // initial sort owner for both source initialization paths remains untouched.

    let updates = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| usize::try_from(atoms[bond.nbr_idx].index).unwrap_or(usize::MAX))
        .collect::<Vec<_>>();
    for (bond, nbr_sym_class) in atoms[atom_idx].bonds.iter_mut().zip(updates) {
        bond.nbr_sym_class = nbr_sym_class;
    }
    let ranks = canon_atom_rank_snapshot(atoms);
    insertion_sort_canon_bonds_for_update(&mut atoms[atom_idx].bonds, |left, right| {
        compare_canon_bond_holder(left, right, &ranks) == Ordering::Greater
    });
}

fn insertion_sort_canon_bonds_for_update<'a, F>(nbrs: &mut [CanonBondHolder<'a>], mut greater: F)
where
    F: FnMut(&CanonBondHolder<'a>, &CanonBondHolder<'a>) -> bool,
{
    // BEGIN COMPLETE RDKit .6 CHEM32 updateAtomNeighborIndex_insertion_loop
    // RDKit✔️✔️:   for (size_t i = 1; i < nbrs.size(); ++i) {
    // RDKit✔️✔️:     if (!bondholder::greater(nbrs[i], nbrs[i - 1])) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto value = std::move(nbrs[i]);
    // RDKit✔️✔️:     size_t j = i;
    // RDKit✔️✔️:     do {
    // RDKit✔️✔️:       nbrs[j] = std::move(nbrs[j - 1]);
    // RDKit✔️✔️:       --j;
    // RDKit✔️✔️:     } while (j && bondholder::greater(value, nbrs[j - 1]));
    // RDKit✔️✔️:     nbrs[j] = std::move(value);
    // RDKit✔️✔️:   }
    // END COMPLETE RDKit .6 CHEM32 updateAtomNeighborIndex_insertion_loop
    // Copy holder contains only scalars, fixed arrays, and borrowed str. One
    // shallow local record and fixed-size shifts reproduce std::move with no
    // allocation. Equal keys never shift; sorted n-1, reverse n(n-1)/2 compares,
    // O(1) auxiliary state and O(n^2) worst case, matching the source loop.
    for i in 1..nbrs.len() {
        if !greater(&nbrs[i], &nbrs[i - 1]) {
            continue;
        }
        let value = nbrs[i];
        let mut j = i;
        loop {
            nbrs[j] = nbrs[j - 1];
            j -= 1;
            if j == 0 || !greater(&value, &nbrs[j - 1]) {
                break;
            }
        }
        nbrs[j] = value;
    }
}

fn update_atom_neighbor_num_swaps_for_kekulize(
    atoms: &[CanonAtom<'_>],
    atom_idx: usize,
) -> Vec<(usize, u32)> {
    // BEGIN RDKIT CPP FUNCTION updateAtomNeighborNumSwaps
    // RDKit✔️✔️: void updateAtomNeighborNumSwaps(
    // RDKit✔️✔️:     canon_atom *atoms, std::vector<bondholder> &nbrs, unsigned int atomIdx,
    // RDKit✔️✔️:     std::vector<std::pair<unsigned int, unsigned int>> &result) {
    // RDKit✔️✔️:   bool isRingAtom = queryIsAtomInRing(atoms[atomIdx].atom);
    let is_ring_atom = atoms[atom_idx].is_ring_atom;
    let mut result = Vec::<(usize, u32)>::new();
    // RDKit✔️✔️:   for (auto &nbr : nbrs) {
    for nbr in &atoms[atom_idx].bonds {
        // RDKit✔️✔️:     unsigned nbrIdx = nbr.nbrIdx;
        let nbr_idx = nbr.nbr_idx;
        // RDKit✔️✔️:     std::list<unsigned int> neighborsSeen;
        // RDKit✔️✔️:     bool tooManySimilarNbrs = false;
        let mut neighbors_seen = Vec::<i32>::new();
        let mut too_many_similar_neighbors = false;
        // RDKit✔️✔️:     if (isRingAtom && atoms[nbrIdx].atom->getChiralTag() != 0) {
        if is_ring_atom && atoms[nbr_idx].chiral_tag != ChiralTag::Unspecified {
            // RDKit✔️✔️:       std::vector<int> ref, probe;
            let mut reference = Vec::<usize>::new();
            // RDKit✔️✔️:       for (unsigned i = 0; i < atoms[nbrIdx].degree; ++i) {
            for nbr_nbr_id in atoms[nbr_idx].nbr_ids.iter().copied() {
                // RDKit✔️✔️:         auto nbrNbrId =
                // RDKit✔️✔️:             atoms[nbrIdx].nbrIds[i];
                // RDKit✔️✔️:         ref.push_back(nbrNbrId);
                reference.push(nbr_nbr_id);
                // RDKit✔️✔️:         if ((int)atomIdx != nbrNbrId) {
                if atom_idx != nbr_nbr_id {
                    // RDKit✔️✔️:           if ((std::find(neighborsSeen.begin(), neighborsSeen.end(),
                    // RDKit✔️✔️:                          atoms[nbrNbrId].index) != neighborsSeen.end())) {
                    // RDKit✔️✔️:             tooManySimilarNbrs = true;
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             neighborsSeen.push_back(atoms[nbrNbrId].index);
                    // RDKit✔️✔️:           }
                    let neighbor_rank = atoms[nbr_nbr_id].index;
                    if neighbors_seen.contains(&neighbor_rank) {
                        too_many_similar_neighbors = true;
                    } else {
                        neighbors_seen.push(neighbor_rank);
                    }
                }
            }
            // RDKit✔️✔️:       probe.push_back(atomIdx);
            let mut probe = vec![atom_idx];
            // RDKit✔️✔️:       for (auto &bond : atoms[nbrIdx].bonds) {
            // RDKit✔️✔️:         if (bond.nbrIdx != atomIdx) {
            // RDKit✔️✔️:           probe.push_back(bond.nbrIdx);
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
            for bond in &atoms[nbr_idx].bonds {
                if bond.nbr_idx != atom_idx {
                    probe.push(bond.nbr_idx);
                }
            }
            // RDKit✔️✔️:       if (tooManySimilarNbrs) {
            // RDKit✔️✔️:         result.emplace_back(nbr.nbrSymClass, 0);
            // RDKit✔️✔️:       } else {
            if too_many_similar_neighbors {
                result.push((nbr.nbr_sym_class, 0));
            } else {
                // RDKit✔️✔️:         int nSwaps = static_cast<int>(countSwapsToInterconvert(ref, probe));
                let swaps = count_swaps_to_interconvert(&reference, probe);
                // RDKit✔️✔️:         if (atoms[nbrIdx].atom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW) {
                // RDKit✔️✔️:           if (nSwaps % 2) {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 2);
                // RDKit✔️✔️:           } else {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 1);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         } else if (atoms[nbrIdx].atom->getChiralTag() ==
                // RDKit✔️✔️:                    Atom::CHI_TETRAHEDRAL_CCW) {
                // RDKit✔️✔️:           if (nSwaps % 2) {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 1);
                // RDKit✔️✔️:           } else {
                // RDKit✔️✔️:             result.emplace_back(nbr.nbrSymClass, 2);
                // RDKit✔️✔️:           }
                // RDKit✔️✔️:         }
                match atoms[nbr_idx].chiral_tag {
                    ChiralTag::TetrahedralCw => {
                        result.push((nbr.nbr_sym_class, if swaps % 2 == 1 { 2 } else { 1 }));
                    }
                    ChiralTag::TetrahedralCcw => {
                        result.push((nbr.nbr_sym_class, if swaps % 2 == 1 { 1 } else { 2 }));
                    }
                    _ => {}
                }
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     } else {
        } else {
            // RDKit✔️✔️:       result.emplace_back(nbr.nbrSymClass, 0);
            result.push((nbr.nbr_sym_class, 0));
        }
        // RDKit✔️✔️:     }
    }
    // RDKit✔️✔️:   sort(result.begin(), result.end());
    // RDKit✔️✔️: }
    result.sort_unstable();
    // END RDKIT CPP FUNCTION updateAtomNeighborNumSwaps
    result
}

fn canon_atom_rank_snapshot(atoms: &[CanonAtom<'_>]) -> Vec<i32> {
    atoms.iter().map(|atom| atom.index).collect()
}

fn sort_canon_bonds_descending(bonds: &mut [CanonBondHolder<'_>], atom_ranks: &[i32]) {
    bonds.sort_by(|left, right| compare_canon_bond_holder(right, left, atom_ranks));
}

fn compare_canon_bond_holder(
    left: &CanonBondHolder<'_>,
    right: &CanonBondHolder<'_>,
    atom_ranks: &[i32],
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION bondholder::compare
    // RDKit✔️✔️: static int compare(const bondholder &x, const bondholder &y,
    // RDKit✔️✔️:                    unsigned int div = 1) {
    // RDKit✔️✔️:   if (x.p_symbol && y.p_symbol) {
    // RDKit✔️✔️:     if ((*x.p_symbol) < (*y.p_symbol)) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if ((*x.p_symbol) > (*y.p_symbol)) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (x.bondType < y.bondType) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   } else if (x.bondType > y.bondType) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if let (Some(left_symbol), Some(right_symbol)) = (left.p_symbol, right.p_symbol) {
        let symbol_cmp = left_symbol.cmp(right_symbol);
        if symbol_cmp != Ordering::Equal {
            return symbol_cmp;
        }
    }
    let base_cmp = rdkit_bond_order_rank(left.bond_type)
        .cmp(&rdkit_bond_order_rank(right.bond_type))
        // RDKit✔️✔️:   if (x.bondStereo < y.bondStereo) {
        // RDKit✔️✔️:     return -1;
        // RDKit✔️✔️:   } else if (x.bondStereo > y.bondStereo) {
        // RDKit✔️✔️:     return 1;
        // RDKit✔️✔️:   }
        .then_with(|| left.bond_stereo.cmp(&right.bond_stereo))
        // RDKit✔️✔️:   auto scdiv = x.nbrSymClass / div - y.nbrSymClass / div;
        // RDKit✔️✔️:   if (scdiv) {
        // RDKit✔️✔️:     return scdiv;
        // RDKit✔️✔️:   }
        .then_with(|| left.nbr_sym_class.cmp(&right.nbr_sym_class));
    if base_cmp != Ordering::Equal {
        return base_cmp;
    }
    // RDKit✔️✔️:   if (x.bondStereo && y.bondStereo) {
    // RDKit✔️✔️:     auto cs = x.compareStereo(y);
    // RDKit✔️✔️:     if (cs) {
    // RDKit✔️✔️:       return cs;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if left.bond_stereo != 0 && right.bond_stereo != 0 {
        let stereo_cmp = compare_canon_bond_stereo(left, right, atom_ranks);
        if stereo_cmp != Ordering::Equal {
            return stereo_cmp;
        }
    }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondholder::compare
    Ordering::Equal
}

fn compare_canon_bond_stereo(
    left: &CanonBondHolder<'_>,
    right: &CanonBondHolder<'_>,
    atom_ranks: &[i32],
) -> Ordering {
    // BEGIN RDKIT CPP FUNCTION bondholder::compareStereo
    // RDKit✔️✔️: int bondholder::compareStereo(const bondholder &o) const {
    // RDKit✔️✔️:   auto st1 = stype;
    // RDKit✔️✔️:   auto st2 = o.stype;
    let mut st1 = left.stype;
    let mut st2 = right.stype;
    // RDKit✔️✔️:   if (st1 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:     if (st2 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if st1 == BondStereo::None {
        return if st2 == BondStereo::None {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    // RDKit✔️✔️:   if (st2 == Bond::BondStereo::STEREONONE) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if st2 == BondStereo::None {
        return Ordering::Greater;
    }
    // RDKit✔️✔️:   if (st1 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:     if (st2 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if st1 == BondStereo::Any {
        return if st2 == BondStereo::Any {
            Ordering::Equal
        } else {
            Ordering::Less
        };
    }
    // RDKit✔️✔️:   if (st2 == Bond::BondStereo::STEREOANY) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    if st2 == BondStereo::Any {
        return Ordering::Greater;
    }
    // RDKit✔️✔️:   // we have some kind of specified stereo on both bonds, work is required
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // if both have absolute stereo labels we can compare them directly
    // RDKit✔️✔️:   if ((st1 == Bond::BondStereo::STEREOE || st1 == Bond::BondStereo::STEREOZ) &&
    // RDKit✔️✔️:       (st2 == Bond::BondStereo::STEREOE || st2 == Bond::BondStereo::STEREOZ)) {
    // RDKit✔️✔️:     if (st1 < st2) {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     } else if (st1 > st2) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    if matches!(st1, BondStereo::E | BondStereo::Z) && matches!(st2, BondStereo::E | BondStereo::Z)
    {
        return rdkit_bond_stereo_rank(st1).cmp(&rdkit_bond_stereo_rank(st2));
    }
    // RDKit✔️✔️:   // check to see if we need to flip the controlling atoms due to atom ranks
    // RDKit✔️✔️:   flipIfNeeded(st1, controllingAtoms);
    // RDKit✔️✔️:   flipIfNeeded(st2, o.controllingAtoms);
    st1 = flip_bond_stereo_if_needed(st1, left.controlling_atoms, atom_ranks);
    st2 = flip_bond_stereo_if_needed(st2, right.controlling_atoms, atom_ranks);
    // RDKit✔️✔️:   if (st1 < st2) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   } else if (st1 > st2) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondholder::compareStereo
    rdkit_bond_stereo_rank(st1).cmp(&rdkit_bond_stereo_rank(st2))
}

fn flip_bond_stereo_if_needed(
    mut stereo: BondStereo,
    controlling_atoms: [Option<usize>; 4],
    atom_ranks: &[i32],
) -> BondStereo {
    // BEGIN RDKIT CPP FUNCTION flipIfNeeded
    // RDKit✔️✔️: void flipIfNeeded(Bond::BondStereo &st1,
    // RDKit✔️✔️:                   const canon_atom *const *controllingAtoms) {
    // RDKit✔️✔️:   CHECK_INVARIANT(controllingAtoms[0], "missing controlling atom");
    // RDKit✔️✔️:   CHECK_INVARIANT(controllingAtoms[2], "missing controlling atom");
    let controlling_0 = controlling_atoms[0].expect("missing controlling atom");
    let controlling_2 = controlling_atoms[2].expect("missing controlling atom");
    // RDKit✔️✔️:   bool flip = false;
    let mut flip = false;
    // RDKit✔️✔️:   if (controllingAtoms[1] &&
    // RDKit✔️✔️:       controllingAtoms[1]->index > controllingAtoms[0]->index) {
    // RDKit✔️✔️:     flip = !flip;
    // RDKit✔️✔️:   }
    if controlling_atoms[1].is_some_and(|idx| atom_ranks[idx] > atom_ranks[controlling_0]) {
        flip = !flip;
    }
    // RDKit✔️✔️:   if (controllingAtoms[3] &&
    // RDKit✔️✔️:       controllingAtoms[3]->index > controllingAtoms[2]->index) {
    // RDKit✔️✔️:     flip = !flip;
    // RDKit✔️✔️:   }
    if controlling_atoms[3].is_some_and(|idx| atom_ranks[idx] > atom_ranks[controlling_2]) {
        flip = !flip;
    }
    // RDKit✔️✔️:   if (flip) {
    // RDKit✔️✔️:     if (st1 == Bond::BondStereo::STEREOCIS) {
    // RDKit✔️✔️:       st1 = Bond::BondStereo::STEREOTRANS;
    // RDKit✔️✔️:     } else if (st1 == Bond::BondStereo::STEREOTRANS) {
    // RDKit✔️✔️:       st1 = Bond::BondStereo::STEREOCIS;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if flip {
        stereo = match stereo {
            BondStereo::Cis => BondStereo::Trans,
            BondStereo::Trans => BondStereo::Cis,
            other => other,
        };
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION flipIfNeeded
    stereo
}

fn rdkit_bond_order_rank(order: BondOrder) -> u8 {
    match order {
        BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Hydrogen => 14,
        BondOrder::ThreeCenter => 15,
        BondOrder::DativeOne => 16,
        BondOrder::Dative => 17,
        BondOrder::DativeLeft => 18,
        BondOrder::DativeRight => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    }
}

fn rdkit_bond_stereo_rank(stereo: BondStereo) -> u8 {
    match stereo {
        BondStereo::None => 0,
        BondStereo::Any => 1,
        BondStereo::Z => 2,
        BondStereo::E => 3,
        BondStereo::Cis => 4,
        BondStereo::Trans => 5,
        BondStereo::AtropCw => 6,
        BondStereo::AtropCcw => 7,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{
        AtomQueryPredicate, AtomSpec, BondQueryPredicate, BondSpec, QueryAtom, QueryBond,
        QueryNode, QueryStateRef,
    };
    use cosmolkit_types::Element;

    fn atom(id: usize, spec: AtomSpec) -> Atom {
        Atom::from_spec(AtomId::new(id), spec)
    }

    fn bond(id: usize, begin: usize, end: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
        )
    }

    fn bond_with_spec(id: usize, spec: BondSpec) -> Bond {
        Bond::from_spec(BondId::new(id), spec)
    }

    fn aromatic_bond(id: usize, begin: usize, end: usize) -> Bond {
        bond_with_spec(
            id,
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                .with_aromatic(true),
        )
    }

    fn kekulize_ready_bond(id: usize, begin: usize, end: usize) -> Bond {
        bond_with_spec(
            id,
            BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single)
                .with_aromatic(true),
        )
    }

    fn valence(explicit_valence: Vec<i32>, implicit_hydrogens: Vec<i32>) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence,
            implicit_hydrogens,
        }
    }

    fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }

    mod drawing_ring_input_tests {
        use super::*;

        fn cycle(n: usize) -> TopologyBlock {
            topology(
                (0..n)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..n)
                    .map(|id| aromatic_bond(id, id, (id + 1) % n))
                    .collect(),
            )
        }

        fn rows(n: usize, quality: RingFindType) -> RingInfo {
            let mut rings = RingInfo::new(quality, n, n);
            rings
                .add_ring(&(0..n).collect::<Vec<_>>(), &(0..n).collect::<Vec<_>>())
                .unwrap();
            rings
        }

        fn prerequisites(input: &TopologyBlock, n: usize, supplied: Option<&RingInfo>) {
            assert_eq!(input.atoms.len(), n);
            assert_eq!(input.bonds.len(), n);
            for id in 0..n {
                assert_eq!(
                    input.atoms[id],
                    atom(id, AtomSpec::new(Element::C).with_aromatic(true))
                );
                assert_eq!(input.bonds[id], aromatic_bond(id, id, (id + 1) % n));
                assert_eq!(input.adjacency.neighbors_of(id).len(), 2);
            }
            if let Some(rings) = supplied {
                assert!(rings.is_sssr_or_better());
                assert_eq!(rings.atom_row_count(), n);
                assert_eq!(rings.bond_row_count(), n);
                assert_eq!(
                    rings.atom_rings(),
                    &[(0..n).map(AtomId::new).collect::<Vec<_>>()]
                );
                assert_eq!(
                    rings.bond_rings(),
                    &[(0..n).map(BondId::new).collect::<Vec<_>>()]
                );
                for id in 0..n {
                    assert_eq!(rings.num_atom_rings(AtomId::new(id)), 1);
                    assert_eq!(rings.num_bond_rings(BondId::new(id)), 1);
                }
            }
        }

        fn forwarded(
            base: usize,
            input: &TopologyBlock,
            params: &KekulizeParams,
            query_present: bool,
            supplied: Option<&RingInfo>,
        ) {
            let calls = drawing_ring_if_possible_probe::calls();
            assert_eq!(calls.len(), base + 1);
            assert_eq!(
                calls[base],
                drawing_ring_if_possible_probe::Forward {
                    topology: input as *const _ as usize,
                    params: *params,
                    query_present,
                    rings: supplied.map(|rings| rings as *const _ as usize),
                }
            );
        }

        #[test]
        fn drawing_ring_kekulize_seventeen_actual_calls() {
            // Frozen K1: six unique inputs repeated twice. Expectations come
            // from retained K04/K06 and pinned source, never this test's output.
            let mut calls = 0usize;
            for _repeat in 0..2 {
                for mark in [false, true] {
                    for quality in [None, Some(RingFindType::Sssr), Some(RingFindType::SymmSssr)] {
                        let input = cycle(6);
                        let supplied = quality.map(|quality| rows(6, quality));
                        prerequisites(&input, 6, supplied.as_ref());
                        let snapshot = input.clone();
                        let ring_snapshot = supplied.clone();
                        let params = KekulizeParams {
                            mark_atoms_bonds: mark,
                            canonical: false,
                            max_backtracks: 100,
                        };
                        let params_snapshot = params;
                        let mut expected = snapshot.clone();
                        for atom in &mut expected.atoms {
                            atom.set_aromatic(!mark);
                        }
                        for (bond, order) in expected.bonds.iter_mut().zip([
                            BondOrder::Double,
                            BondOrder::Single,
                            BondOrder::Double,
                            BondOrder::Single,
                            BondOrder::Double,
                            BondOrder::Single,
                        ]) {
                            bond.set_order(order);
                            bond.set_aromatic(!mark);
                        }
                        let base = drawing_ring_if_possible_probe::calls().len();
                        let acquisitions = ring_transport_probe::sssr_calls();
                        let ranks = ring_transport_probe::rank_calls();
                        let result = kekulize_if_possible_with_query_state_and_ring_info(
                            &input,
                            &params,
                            None,
                            supplied.as_ref(),
                        );
                        calls += 1;
                        assert_eq!(input, snapshot);
                        assert_eq!(supplied, ring_snapshot);
                        assert_eq!(params, params_snapshot);
                        forwarded(base, &input, &params, false, supplied.as_ref());
                        assert_eq!(ring_transport_probe::rank_calls(), ranks);
                        assert_eq!(
                            ring_transport_probe::sssr_calls() - acquisitions,
                            if supplied.is_none() { 1 } else { 0 }
                        );
                        let KekulizeAttempt::Applied(assignment) = result.unwrap() else {
                            panic!("K1 must apply");
                        };
                        assert_eq!(assignment.topology, expected);
                        if supplied.is_some() {
                            assert!(assignment.ring_update.is_none());
                        } else {
                            let update = assignment.ring_update.unwrap();
                            assert_eq!(update.find_type(), RingFindType::Sssr);
                            assert_eq!(update.atom_row_count(), 6);
                            assert_eq!(update.bond_row_count(), 6);
                            assert_eq!(update.atom_rings().len(), 1);
                            assert_eq!(update.bond_rings().len(), 1);
                            assert_eq!(update.atom_rings()[0].len(), 6);
                            assert_eq!(update.bond_rings()[0].len(), 6);
                            for id in 0..6 {
                                assert_eq!(update.num_atom_rings(AtomId::new(id)), 1);
                                assert_eq!(update.num_bond_rings(BondId::new(id)), 1);
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 12);
            // Frozen K2: preserve original fallback; no failure-side ring-update claim.
            for quality in [None, Some(RingFindType::Sssr)] {
                let input = cycle(5);
                let supplied = quality.map(|quality| rows(5, quality));
                prerequisites(&input, 5, supplied.as_ref());
                let snapshot = input.clone();
                let ring_snapshot = supplied.clone();
                let params = KekulizeParams {
                    mark_atoms_bonds: false,
                    canonical: false,
                    max_backtracks: 100,
                };
                let base = drawing_ring_if_possible_probe::calls().len();
                let result = kekulize_if_possible_with_query_state_and_ring_info(
                    &input,
                    &params,
                    None,
                    supplied.as_ref(),
                );
                calls += 1;
                assert_eq!(input, snapshot);
                assert_eq!(supplied, ring_snapshot);
                forwarded(base, &input, &params, false, supplied.as_ref());
                assert_eq!(
                    result,
                    Ok(KekulizeAttempt::NotKekulizable {
                        topology: snapshot,
                        problem_atoms: vec![
                            AtomId::new(0),
                            AtomId::new(1),
                            AtomId::new(2),
                            AtomId::new(3),
                            AtomId::new(4)
                        ]
                    })
                );
            }
            assert_eq!(calls, 14);
            // Frozen K3: actual malformed topology, foreign query, consumed dimensions.
            let params = KekulizeParams {
                mark_atoms_bonds: false,
                canonical: false,
                max_backtracks: 100,
            };
            let mut malformed = cycle(6);
            malformed.atoms[0] = atom(1, AtomSpec::new(Element::C).with_aromatic(true));
            let snapshot = malformed.clone();
            let base = drawing_ring_if_possible_probe::calls().len();
            let result = kekulize_if_possible_with_query_state_and_ring_info(
                &malformed, &params, None, None,
            );
            calls += 1;
            assert_eq!(malformed, snapshot);
            forwarded(base, &malformed, &params, false, None);
            assert_eq!(
                result,
                Err(KekulizeError::InvalidTopology(
                    TopologyValidationError::AtomIdMismatch {
                        position: 0,
                        id: AtomId::new(1)
                    }
                ))
            );

            let input = cycle(6);
            let snapshot = input.clone();
            let foreign = cycle(5);
            let foreign_snapshot = foreign.clone();
            let query_atoms = foreign
                .atoms
                .iter()
                .map(|carrier| {
                    QueryAtom::from_carrier_parts(
                        carrier.clone(),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect::<Vec<_>>();
            let query_bonds = foreign
                .bonds
                .iter()
                .map(|carrier| {
                    QueryBond::from_carrier_parts(
                        carrier.clone(),
                        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
                    )
                })
                .collect::<Vec<_>>();
            let atom_snapshot = query_atoms.clone();
            let bond_snapshot = query_bonds.clone();
            let query =
                QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &foreign).unwrap();
            let base = drawing_ring_if_possible_probe::calls().len();
            let result = kekulize_if_possible_with_query_state_and_ring_info(
                &input,
                &params,
                Some(query),
                None,
            );
            calls += 1;
            assert_eq!(input, snapshot);
            assert_eq!(foreign, foreign_snapshot);
            assert_eq!(query_atoms, atom_snapshot);
            assert_eq!(query_bonds, bond_snapshot);
            forwarded(base, &input, &params, true, None);
            assert_eq!(
                result,
                Err(KekulizeError::InvalidQueryState(
                    QueryStateError::AtomCount {
                        actual: 5,
                        expected: 6
                    }
                ))
            );

            let supplied = crate::ring_info_from_selected_rows(6, 5, &[], &[]).unwrap();
            let ring_snapshot = supplied.clone();
            assert!(supplied.is_initialized());
            assert_eq!(supplied.find_type(), RingFindType::OtherOrUnknown);
            assert_eq!(supplied.atom_row_count(), 6);
            assert_eq!(supplied.bond_row_count(), 5);
            let base = drawing_ring_if_possible_probe::calls().len();
            let result = kekulize_if_possible_with_query_state_and_ring_info(
                &input,
                &params,
                None,
                Some(&supplied),
            );
            calls += 1;
            assert_eq!(input, snapshot);
            assert_eq!(supplied, ring_snapshot);
            forwarded(base, &input, &params, false, Some(&supplied));
            assert_eq!(
                result,
                Err(KekulizeError::CanonicalRank(
                    CanonicalRankError::PreparedRingLength {
                        expected_atoms: 6,
                        actual_atoms: 6,
                        expected_bonds: 6,
                        actual_bonds: 5
                    }
                ))
            );
            assert_eq!(calls, 17);
        }
    }

    #[test]
    fn kekulize_source_selection_independent_flags_and_masks_product() {
        // Frozen SOURCE-K B32; Atom.cpp isAromaticAtom inspects incident
        // order/flags even when that bond is outside the selected bond mask.
        let mut calls = 0;
        for order in [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Aromatic,
        ] {
            for bond_flag in [false, true] {
                for atom_flag in [false, true] {
                    for selected in [false, true] {
                        let input = topology(
                            vec![
                                atom(0, AtomSpec::new(Element::C).with_aromatic(atom_flag)),
                                atom(1, AtomSpec::new(Element::C).with_aromatic(atom_flag)),
                            ],
                            vec![bond_with_spec(
                                0,
                                BondSpec::new(AtomId::new(0), AtomId::new(1), order)
                                    .with_aromatic(bond_flag),
                            )],
                        );
                        for (id, row) in input.atoms.iter().enumerate() {
                            assert_eq!(row.id(), AtomId::new(id));
                            assert_eq!(row.element(), Element::C);
                            assert_eq!(row.is_aromatic(), atom_flag);
                        }
                        assert_eq!(input.bonds[0].id(), BondId::new(0));
                        assert_eq!(input.bonds[0].begin(), AtomId::new(0));
                        assert_eq!(input.bonds[0].end(), AtomId::new(1));
                        assert_eq!(input.bonds[0].order(), order);
                        assert_eq!(input.bonds[0].is_aromatic(), bond_flag);
                        let snapshot = input.clone();
                        let result =
                            prepare_kekulize_core(&input, &[true; 2], &[selected], None, None);
                        calls += 1;
                        assert_eq!(input, snapshot, "retained input changed at call {calls}");
                        let prepared = result.unwrap();
                        assert_eq!(
                            prepared.found_aromatic,
                            atom_flag || bond_flag || order == BondOrder::Aromatic,
                            "order={order:?}, bond_flag={bond_flag}, atom_flag={atom_flag}, selected={selected}"
                        );
                        assert_eq!(prepared.bonds_in_play, vec![selected]);
                        assert_eq!(prepared.atoms_in_play, vec![true; 2]);
                    }
                }
            }
        }
        assert_eq!(calls, 32);
        eprintln!("SOURCE-K B32 actual selection calls: {calls}");
    }

    #[test]
    fn selected_kekulize_contract_checks_both_full_index_mask_lengths_before_empty_return() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let before = graph.clone();
        let atoms = [false, false];
        let bonds = [true];
        assert_eq!(
            kekulize_selected_fragment(&graph, &[false], &[], &KekulizeParams::default()),
            Err(KekulizeError::AtomSelectionLength {
                expected: 2,
                actual: 1,
            })
        );
        assert_eq!(
            kekulize_selected_fragment(&graph, &atoms, &[], &KekulizeParams::default()),
            Err(KekulizeError::BondSelectionLength {
                expected: 1,
                actual: 0,
            })
        );
        assert_eq!(graph, before);
        assert_eq!(atoms, [false, false]);
        assert_eq!(bonds, [true]);
    }

    fn line130_cached_source_graph() -> TopologyBlock {
        let elements = [
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
            Element::N,
            Element::C,
            Element::C,
            Element::C,
            Element::C,
            Element::O,
            Element::O,
        ];
        let aromatic_atoms = [0, 1, 2, 3, 4, 9, 10];
        let atoms = elements
            .into_iter()
            .enumerate()
            .map(|(index, element)| {
                atom(
                    index,
                    AtomSpec::new(element).with_aromatic(aromatic_atoms.contains(&index)),
                )
            })
            .collect();
        let bond_rows = [
            (0, 1, BondOrder::Aromatic),
            (1, 2, BondOrder::Aromatic),
            (2, 3, BondOrder::Aromatic),
            (3, 4, BondOrder::Aromatic),
            (4, 5, BondOrder::Single),
            (5, 6, BondOrder::Single),
            (6, 7, BondOrder::Single),
            (7, 8, BondOrder::Single),
            (8, 9, BondOrder::Double),
            (9, 10, BondOrder::Aromatic),
            (9, 11, BondOrder::Single),
            (9, 12, BondOrder::Single),
            (6, 13, BondOrder::Double),
            (1, 14, BondOrder::Single),
            (10, 0, BondOrder::Aromatic),
            (10, 4, BondOrder::Aromatic),
        ];
        let bonds = bond_rows
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end, order))| {
                bond_with_spec(
                    index,
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), order)
                        .with_aromatic(order == BondOrder::Aromatic),
                )
            })
            .collect();
        topology(atoms, bonds)
    }

    fn five_member_n_p_cache_graph(element: Element) -> TopologyBlock {
        let atoms = (0..5)
            .map(|index| {
                let spec = if index == 0 {
                    AtomSpec::new(element)
                        .with_aromatic(true)
                        .with_explicit_hydrogens(1)
                        .with_no_implicit(true)
                } else {
                    AtomSpec::new(Element::C).with_aromatic(true)
                };
                atom(index, spec)
            })
            .collect();
        let bonds = (0..5)
            .map(|index| aromatic_bond(index, index, (index + 1) % 5))
            .collect();
        topology(atoms, bonds)
    }

    #[test]
    fn kekulize_source_cached_postcondition_() {
        let cases = [
            (
                "line130",
                line130_cached_source_graph(),
                false,
                vec![3, 4, 3, 3, 4, 2, 4, 2, 3, 4, 4, 1, 1, 2, 1],
                vec![1, 0, 1, 1, 0, 1, 0, 2, 0, 0, 0, 3, 3, 0, 1],
            ),
            (
                "pyrrole-n",
                five_member_n_p_cache_graph(Element::N),
                true,
                vec![3, 3, 3, 3, 3],
                vec![0, 1, 1, 1, 1],
            ),
            (
                "phosphole-p",
                five_member_n_p_cache_graph(Element::P),
                true,
                vec![3, 3, 3, 3, 3],
                vec![0, 1, 1, 1, 1],
            ),
        ];
        let mut calls = 0;
        let mut mismatches = Vec::new();

        for (name, graph, refresh_center, initial_explicit, initial_implicit) in cases {
            let input_snapshot = graph.clone();
            for mark_atoms_bonds in [false, true] {
                for canonical in [false, true] {
                    calls += 1;
                    let params = KekulizeParams {
                        mark_atoms_bonds,
                        canonical,
                        ..KekulizeParams::default()
                    };
                    match kekulize(&graph, &params) {
                        Err(error) => mismatches.push(format!(
                            "{name}, mark_atoms_bonds={mark_atoms_bonds}, canonical={canonical}: unexpected {error:?}"
                        )),
                        Ok(assignment) => {
                            let expected_explicit = if refresh_center && mark_atoms_bonds {
                                vec![2, 3, 3, 3, 3]
                            } else {
                                initial_explicit.clone()
                            };
                            let expected_implicit = if refresh_center && mark_atoms_bonds {
                                vec![1, 1, 1, 1, 1]
                            } else {
                                initial_implicit.clone()
                            };
                            let expected_refreshed = if refresh_center && mark_atoms_bonds {
                                vec![AtomId::new(0)]
                            } else {
                                Vec::new()
                            };
                            match assignment.final_valence {
                                Some(actual)
                                    if actual.explicit_valence == expected_explicit
                                        && actual.implicit_hydrogens == expected_implicit => {}
                                other => mismatches.push(format!(
                                    "{name}, mark_atoms_bonds={mark_atoms_bonds}, canonical={canonical}: expected final cache {expected_explicit:?}/{expected_implicit:?}, observed {other:?}"
                                )),
                            }
                            if assignment.refreshed_valence_atoms != expected_refreshed {
                                mismatches.push(format!(
                                    "{name}, mark_atoms_bonds={mark_atoms_bonds}, canonical={canonical}: expected refreshed IDs {expected_refreshed:?}, observed {:?}",
                                    assignment.refreshed_valence_atoms
                                ));
                            }
                            for (index, (before, after)) in input_snapshot
                                .atoms
                                .iter()
                                .zip(&assignment.topology.atoms)
                                .enumerate()
                            {
                                let expected_aromatic = before.is_aromatic() && !mark_atoms_bonds;
                                let expected_hydrogens = if index == 0 && refresh_center {
                                    if mark_atoms_bonds { 0 } else { 1 }
                                } else {
                                    before.explicit_hydrogens()
                                };
                                let expected_no_implicit = if index == 0 && refresh_center {
                                    !mark_atoms_bonds
                                } else {
                                    before.no_implicit()
                                };
                                if after.is_aromatic() != expected_aromatic
                                    || after.explicit_hydrogens() != expected_hydrogens
                                    || after.no_implicit() != expected_no_implicit
                                {
                                    mismatches.push(format!(
                                        "{name}, mark_atoms_bonds={mark_atoms_bonds}, canonical={canonical}: atom {index} flags expected aromatic={expected_aromatic}, H={expected_hydrogens}, noImplicit={expected_no_implicit}; observed aromatic={}, H={}, noImplicit={}",
                                        after.is_aromatic(),
                                        after.explicit_hydrogens(),
                                        after.no_implicit()
                                    ));
                                }
                            }
                            if graph != input_snapshot {
                                mismatches.push(format!(
                                    "{name}, mark_atoms_bonds={mark_atoms_bonds}, canonical={canonical}: input graph changed"
                                ));
                            }
                        }
                    }
                }
            }
        }

        assert_eq!(calls, 12, "all frozen owner calls were attempted");
        assert!(
            mismatches.is_empty(),
            "source cache matrix mismatches after all {calls} calls:\n{}",
            mismatches.join("\n")
        );
    }

    #[test]
    fn selected_kekulize_contract_keeps_original_ids_and_caller_state() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let before = graph.clone();
        let none = [false, false];
        let selected_atom = [true, false];
        let selected_bond = [true];
        let unselected_bond = [false];
        for (atoms, bonds) in [
            (&none[..], &selected_bond[..]),
            (&selected_atom[..], &selected_bond[..]),
            (&selected_atom[..], &unselected_bond[..]),
        ] {
            let original_atoms = atoms.to_vec();
            let original_bonds = bonds.to_vec();
            let result = kekulize_selected_fragment(&graph, atoms, bonds, &Default::default())
                .expect("source accepts independent original-index selection masks");
            assert_eq!(result.topology, before);
            assert_eq!(atoms, original_atoms);
            assert_eq!(bonds, original_bonds);
        }
        assert_eq!(graph, before);
    }

    fn selected_kekulize_behavior_ring_pair() -> TopologyBlock {
        topology(
            (0..12)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 0),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 10),
                (10, 11),
                (11, 6),
            ]
            .into_iter()
            .enumerate()
            .map(|(id, (begin, end))| aromatic_bond(id, begin, end))
            .collect(),
        )
    }

    #[test]
    fn selected_kekulize_behavior_full_masks_match_whole_entry_for_both_flags() {
        let graph = selected_kekulize_behavior_ring_pair();
        let before = graph.clone();
        for mark_atoms_bonds in [false, true] {
            for canonical in [false, true] {
                let params = KekulizeParams {
                    mark_atoms_bonds,
                    canonical,
                    ..Default::default()
                };
                let selected =
                    kekulize_selected_fragment(&graph, &[true; 12], &[true; 12], &params).unwrap();
                assert_eq!(selected, kekulize(&graph, &params).unwrap());
            }
        }
        assert_eq!(graph, before);
    }

    #[test]
    fn selected_kekulize_behavior_partial_ring_preserves_original_unselected_rows() {
        let graph = selected_kekulize_behavior_ring_pair();
        let before = graph.clone();
        let mut atoms = [false; 12];
        atoms[0] = true;
        atoms[1] = true;
        let mut bonds = [false; 12];
        bonds[0] = true;
        let result = kekulize_selected_fragment(&graph, &atoms, &bonds, &Default::default())
            .expect("pinned partial ring has no wholly selected candidate ring");
        assert_eq!(result.topology.atoms.len(), 12);
        assert_eq!(result.topology.bonds.len(), 12);
        assert!(!result.topology.atoms[0].is_aromatic());
        assert!(!result.topology.atoms[1].is_aromatic());
        assert!(!result.topology.bonds[0].is_aromatic());
        assert_eq!(result.topology.bonds[0].order(), BondOrder::Aromatic);
        assert_eq!(&result.topology.atoms[2..], &before.atoms[2..]);
        assert_eq!(&result.topology.bonds[1..], &before.bonds[1..]);
        assert_ne!(
            result.topology,
            kekulize(&graph, &Default::default()).unwrap().topology
        );
        assert_eq!(graph, before);
        assert_eq!(atoms.iter().filter(|selected| **selected).count(), 2);
        assert_eq!(bonds.iter().filter(|selected| **selected).count(), 1);
    }

    #[test]
    fn selected_kekulize_behavior_disconnected_ring_keeps_unselected_component() {
        let graph = selected_kekulize_behavior_ring_pair();
        let before = graph.clone();
        let atoms = [
            true, true, true, true, true, true, false, false, false, false, false, false,
        ];
        let bonds = atoms;
        let result = kekulize_selected_fragment(&graph, &atoms, &bonds, &Default::default())
            .expect("one fully selected aromatic component is kekulizable");
        assert!(
            result.topology.atoms[..6]
                .iter()
                .all(|atom| !atom.is_aromatic())
        );
        assert!(
            result.topology.bonds[..6]
                .iter()
                .all(|bond| !bond.is_aromatic())
        );
        assert_eq!(&result.topology.atoms[6..], &before.atoms[6..]);
        assert_eq!(&result.topology.bonds[6..], &before.bonds[6..]);
        assert_eq!(graph, before);
    }

    #[test]
    fn selected_kekulize_behavior_full_nonring_aromatic_atom_reports_source_error() {
        let graph = topology(
            vec![atom(0, AtomSpec::new(Element::C).with_aromatic(true))],
            vec![],
        );
        let before = graph.clone();
        assert_eq!(
            kekulize_selected_fragment(&graph, &[true], &[], &Default::default()),
            Err(KekulizeError::AromaticAtomOutsideRing {
                atom: AtomId::new(0)
            })
        );
        assert_eq!(graph, before);
    }

    #[test]
    fn selected_kekulize_behavior_uninitialized_rings_use_pinned_sssr_rows() {
        let graph = topology(
            (0..10)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 6),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 0),
                (8, 3),
            ]
            .into_iter()
            .enumerate()
            .map(|(id, (begin, end))| aromatic_bond(id, begin, end))
            .collect(),
        );
        let selected = [true; 10];
        let prepared = prepare_kekulize_selection(&graph, &selected, &[true; 11], None).unwrap();
        let atom_rows = prepared
            .candidate_atom_rings
            .iter()
            .map(|ring| ring.iter().map(|atom| atom.index()).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let bond_rows = prepared
            .candidate_bond_rings
            .iter()
            .map(|ring| ring.iter().map(|bond| bond.index()).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(
            atom_rows,
            vec![vec![0, 9, 8, 3, 2, 1], vec![4, 5, 6, 7, 8, 3]]
        );
        assert_eq!(
            bond_rows,
            vec![vec![9, 8, 10, 2, 1, 0], vec![4, 5, 6, 7, 10, 3]]
        );
    }

    mod kekulize_ring_k02_candidate_owner {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene_rows() -> RingInfo {
            ring_info_from_selected_rows(
                6,
                6,
                &[vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                    AtomId::new(4),
                    AtomId::new(5),
                ]],
                &[vec![
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(5),
                ]],
            )
            .unwrap()
        }

        fn ids(rows: &[Vec<AtomId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|atom| atom.index()).collect())
                .collect()
        }

        fn bond_ids(rows: &[Vec<BondId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|bond| bond.index()).collect())
                .collect()
        }

        fn wedged(topology: &TopologyBlock, selected_bonds: &[bool]) -> Vec<bool> {
            let mut wedged_atoms = vec![false; topology.atoms.len()];
            for bond in &topology.bonds {
                if selected_bonds[bond.id().index()]
                    && matches!(
                        bond.direction(),
                        BondDirection::BeginWedge | BondDirection::BeginDash
                    )
                {
                    wedged_atoms[bond.begin().index()] = true;
                }
            }
            wedged_atoms
        }

        // K02 call 1/2 (all-atoms/all-bonds, wedge bond b0 at ring position 0):
        // rotation by startPos 0 keeps the stored order; the wedge marks the
        // row as front-inserted. BEGINWEDGE and BEGINDASH are the same source
        // wedgedAtoms signal (Kekulize.cpp:634-637).
        #[test]
        fn kekulize_ring_k02_all_selection_rotation_identity_both_directions() {
            for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                let graph = topology(
                    (0..6)
                        .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                        .collect(),
                    six_cycle_with(direction, 0),
                );
                let selected = [true; 6];
                let wedged_atoms = wedged(&graph, &selected);
                assert!(wedged_atoms[0]);
                let rings = benzene_rows();
                let before = rings.clone();
                let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                    &graph,
                    &rings,
                    &selected,
                    &selected,
                    &[false; 6],
                    &wedged_atoms,
                )
                .unwrap();
                assert_eq!(ids(&atom_rows), vec![vec![0, 1, 2, 3, 4, 5]]);
                assert_eq!(bond_ids(&bond_rows), vec![vec![0, 1, 2, 3, 4, 5]]);
                // Borrowed assignment is never modified by candidate copying.
                assert_eq!(rings, before);
            }
        }

        fn six_cycle_with(direction: BondDirection, wedge_bond: usize) -> Vec<Bond> {
            (0..6)
                .map(|id| {
                    let spec = BondSpec::new(
                        AtomId::new(id),
                        AtomId::new((id + 1) % 6),
                        BondOrder::Aromatic,
                    )
                    .with_aromatic(true);
                    let spec = if id == wedge_bond {
                        spec.with_direction(direction)
                    } else {
                        spec
                    };
                    bond_with_spec(id, spec)
                })
                .collect()
        }

        // K02 calls 3/4 (partial atom: disconnected two-cycle graph, only the
        // first cycle's atoms selected): only ring-1 rows are copied; ring-2
        // rows are excluded by the atomsInPlay member check, never by graph
        // substitution (containsNonDummy loop, Kekulize.cpp:659-670).
        #[test]
        fn kekulize_ring_k02_partial_atom_excludes_whole_second_ring() {
            for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                let mut atoms = Vec::new();
                for id in 0..12 {
                    atoms.push(atom(id, AtomSpec::new(Element::C).with_aromatic(true)));
                }
                let mut bonds = six_cycle_with(direction, 0);
                for id in 6..12 {
                    bonds.push(bond_with_spec(
                        id,
                        BondSpec::new(
                            AtomId::new(id),
                            AtomId::new(6 + (id - 6 + 1) % 6),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true),
                    ));
                }
                let graph = topology(atoms, bonds);
                let mut atoms_in_play = [false; 12];
                for selected in &mut atoms_in_play[..6] {
                    *selected = true;
                }
                let bonds_in_play = [true; 12];
                let wedged_atoms = wedged(&graph, &bonds_in_play);
                assert!(wedged_atoms[0]);
                assert!(!wedged_atoms[6..].iter().any(|wedge| *wedge));
                let rings = ring_info_from_selected_rows(
                    12,
                    12,
                    &[
                        vec![
                            AtomId::new(0),
                            AtomId::new(1),
                            AtomId::new(2),
                            AtomId::new(3),
                            AtomId::new(4),
                            AtomId::new(5),
                        ],
                        vec![
                            AtomId::new(6),
                            AtomId::new(7),
                            AtomId::new(8),
                            AtomId::new(9),
                            AtomId::new(10),
                            AtomId::new(11),
                        ],
                    ],
                    &[
                        vec![
                            BondId::new(0),
                            BondId::new(1),
                            BondId::new(2),
                            BondId::new(3),
                            BondId::new(4),
                            BondId::new(5),
                        ],
                        vec![
                            BondId::new(6),
                            BondId::new(7),
                            BondId::new(8),
                            BondId::new(9),
                            BondId::new(10),
                            BondId::new(11),
                        ],
                    ],
                )
                .unwrap();
                let before = rings.clone();
                let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                    &graph,
                    &rings,
                    &atoms_in_play,
                    &bonds_in_play,
                    &[false; 12],
                    &wedged_atoms,
                )
                .unwrap();
                assert_eq!(ids(&atom_rows), vec![vec![0, 1, 2, 3, 4, 5]]);
                assert_eq!(bond_ids(&bond_rows), vec![vec![0, 1, 2, 3, 4, 5]]);
                assert_eq!(rings, before);
            }
        }

        // K02 calls 5/6 (partial bond: the ring's b0 deselected so the single
        // stored row fails copyBondRingsWithinFragment): no candidates, and
        // the deselected b0 can no longer contribute a wedge start.
        #[test]
        fn kekulize_ring_k02_partial_bond_excludes_the_ring_row() {
            for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                let graph = topology(
                    (0..6)
                        .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                        .collect(),
                    six_cycle_with(direction, 5),
                );
                let mut bonds_in_play = [true; 6];
                bonds_in_play[0] = false;
                let wedged_atoms = wedged(&graph, &bonds_in_play);
                assert!(wedged_atoms[5]);
                let rings = benzene_rows();
                let before = rings.clone();
                let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                    &graph,
                    &rings,
                    &[true; 6],
                    &bonds_in_play,
                    &[false; 6],
                    &wedged_atoms,
                )
                .unwrap();
                assert!(atom_rows.is_empty());
                assert!(bond_rows.is_empty());
                assert_eq!(rings, before);
            }
        }

        // K02 calls 7/8 (wedge start at ring position 2 via bond b2 whose
        // begin atom is atom 2): the stored row rotates so the wedged atom
        // leads and the row is front-inserted (Kekulize.cpp:675-688).
        #[test]
        fn kekulize_ring_k02_wedge_start_rotates_wedged_atom_to_front() {
            for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                let graph = topology(
                    (0..6)
                        .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                        .collect(),
                    six_cycle_with(direction, 2),
                );
                let selected = [true; 6];
                let wedged_atoms = wedged(&graph, &selected);
                assert!(wedged_atoms[2]);
                let rings = benzene_rows();
                let before = rings.clone();
                let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                    &graph,
                    &rings,
                    &selected,
                    &selected,
                    &[false; 6],
                    &wedged_atoms,
                )
                .unwrap();
                assert_eq!(ids(&atom_rows), vec![vec![2, 3, 4, 5, 0, 1]]);
                assert_eq!(bond_ids(&bond_rows), vec![vec![2, 3, 4, 5, 0, 1]]);
                assert_eq!(rings, before);
            }
        }
    }

    mod kekulize_ring_k03_transition {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene() -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                    .collect(),
            )
        }

        fn full_rows() -> RingInfo {
            ring_info_from_selected_rows(
                6,
                6,
                &[vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                    AtomId::new(4),
                    AtomId::new(5),
                ]],
                &[vec![
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(5),
                ]],
            )
            .unwrap()
        }

        fn empty_rows(find_type: RingFindType) -> RingInfo {
            RingInfo::new(find_type, 6, 6)
        }

        fn reset_state() -> RingInfo {
            let mut rings = RingInfo::new(RingFindType::OtherOrUnknown, 6, 6);
            rings.reset();
            rings
        }

        fn typed_full(find_type: RingFindType) -> RingInfo {
            let mut rings = full_rows();
            rings.initialize(find_type);
            rings
        }

        #[derive(Clone, Copy, PartialEq, Eq, Debug)]
        enum Expect {
            /// Supplied state unchanged: ring_update None.
            PreserveNone,
            /// Fresh SSSR acquired: ring_update Some(initialized SSSR rows).
            AcquireSssr,
            /// Canonical ranking reset an initialized Other state and no
            /// acquisition ran: ring_update Some(uninitialized).
            ResetOnly,
        }

        fn run_case(
            graph: &TopologyBlock,
            valence: &ValenceAssignment,
            name: &str,
            rings: Option<&RingInfo>,
            canonical: bool,
            acquisition: bool,
            expect: Expect,
        ) {
            let supplied_snapshot = rings.cloned();
            let rank_base = ring_transport_probe::rank_calls();
            let sssr_base = ring_transport_probe::sssr_calls();
            let (rows, replaced_by_reset, _atom_ranks) = kekulize_ring_state_transition(
                graph,
                valence,
                &[true; 6],
                &[true; 6],
                canonical,
                acquisition,
                rings,
            )
            .unwrap_or_else(|error| panic!("{name}: unexpected error {error:?}"));
            // Per-call input identity checkpoint: the borrowed supplied
            // assignment is fully unchanged after the transition.
            assert_eq!(rings.cloned(), supplied_snapshot, "{name}: input changed");
            // Counter deltas captured at the actual rank/SSSR sites.
            assert_eq!(
                ring_transport_probe::rank_calls() - rank_base,
                u64::from(canonical),
                "{name}: rank counter delta"
            );
            let expect_sssr_delta = matches!(expect, Expect::AcquireSssr);
            assert_eq!(
                ring_transport_probe::sssr_calls() - sssr_base,
                u64::from(expect_sssr_delta),
                "{name}: SSSR counter delta"
            );
            match expect {
                Expect::PreserveNone => {
                    assert!(!replaced_by_reset, "{name}: unexpected reset disposition");
                    match rows {
                        KekulizeRingRows::Borrowed(borrowed) => {
                            let Some(supplied) = rings else {
                                panic!("{name}: borrowed rows without supplied state")
                            };
                            assert!(std::ptr::eq(borrowed, supplied), "{name}: identity");
                        }
                        KekulizeRingRows::Uninitialized => {
                            assert!(
                                rings.is_none() || !rings.unwrap().is_initialized(),
                                "{name}: uninitialized rows with initialized supply"
                            );
                        }
                        KekulizeRingRows::Acquired(_) => {
                            panic!("{name}: unexpected acquisition")
                        }
                        KekulizeRingRows::Live(_) => {
                            panic!("readonly transition never yields live rows")
                        }
                    }
                }
                Expect::AcquireSssr => {
                    // The one acquired state is owned here: consume it to
                    // inspect the moved buffer, exactly as the engine does.
                    let KekulizeRingRows::Acquired(update) = rows else {
                        panic!("{name}: expected acquired state")
                    };
                    assert!(update.is_initialized(), "{name}: update uninitialized");
                    assert_eq!(update.find_type(), RingFindType::Sssr, "{name}: type");
                    assert_eq!(update.atom_rings().len(), 1, "{name}: row count");
                    assert_eq!(update.atom_rings()[0].len(), 6, "{name}: ring size");
                }
                Expect::ResetOnly => {
                    assert!(replaced_by_reset, "{name}: missing reset disposition");
                    assert!(
                        matches!(rows, KekulizeRingRows::Uninitialized),
                        "{name}: rows"
                    );
                }
            }
        }

        #[test]
        fn kekulize_ring_k03_quality_predicates_on_initialized_empty_states() {
            for find_type in [
                RingFindType::Fast,
                RingFindType::Sssr,
                RingFindType::SymmSssr,
            ] {
                let empty = empty_rows(find_type);
                assert!(empty.is_initialized());
                assert_eq!(empty.atom_rings().len(), 0);
                assert_eq!(empty.bond_rings().len(), 0);
                assert!(empty.is_find_fast_or_better());
                assert_eq!(empty.is_sssr_or_better(), find_type != RingFindType::Fast);
                assert_eq!(empty.is_symm_sssr(), find_type == RingFindType::SymmSssr);
            }
            let reset = reset_state();
            assert!(!reset.is_initialized());
            assert!(!reset.is_find_fast_or_better());
            assert!(!reset.is_sssr_or_better());
            assert!(!reset.is_symm_sssr());
            assert_eq!(reset.atom_row_count(), 0);
            assert_eq!(reset.bond_row_count(), 0);
        }

        // The frozen K01 40-call disposition table: 10 ring inputs x
        // canonical {false,true} x acquisition {false,true}. Each loop
        // iteration is one actual private-owner call; the dispositions come
        // from the K01 literal table, never from observed output.
        #[test]
        fn kekulize_ring_k03_forty_transition_dispositions() {
            let graph = benzene();
            let valence = crate::assign_valence_with_options_from_parts(
                &graph.atoms,
                &graph.bonds,
                &graph.adjacency,
                ValenceModel::RdkitLike,
                false,
            )
            .unwrap();
            let atoms_in_play = [true; 6];
            let bonds_in_play = [true; 6];
            let _ = (&atoms_in_play, &bonds_in_play);
            for canonical in [false, true] {
                for acquisition in [false, true] {
                    // absent (1-10 in the table): acquires only when the
                    // caller guard requires it; canonical ranking of absent
                    // storage changes nothing observable.
                    run_case(
                        &graph,
                        &valence,
                        "absent",
                        None,
                        canonical,
                        acquisition,
                        if acquisition {
                            Expect::AcquireSssr
                        } else {
                            Expect::PreserveNone
                        },
                    );
                    // reset (11-20): reset of uninitialized storage is a
                    // no-op; same family as absent.
                    let reset = reset_state();
                    run_case(
                        &graph,
                        &valence,
                        "reset",
                        Some(&reset),
                        canonical,
                        acquisition,
                        if acquisition {
                            Expect::AcquireSssr
                        } else {
                            Expect::PreserveNone
                        },
                    );
                    // Other-empty / Other-full (21-30): initialized Other
                    // survives noncanonical; canonical ranking resets it,
                    // then acquisition replaces the reset with fresh SSSR.
                    for (name, other) in [
                        ("other-empty", empty_rows(RingFindType::OtherOrUnknown)),
                        ("other-full", typed_full(RingFindType::OtherOrUnknown)),
                    ] {
                        let expect = match (canonical, acquisition) {
                            (false, _) => Expect::PreserveNone,
                            (true, false) => Expect::ResetOnly,
                            (true, true) => Expect::AcquireSssr,
                        };
                        run_case(
                            &graph,
                            &valence,
                            name,
                            Some(&other),
                            canonical,
                            acquisition,
                            expect,
                        );
                    }
                    // Fast/SSSR/Symm, empty and full (31-40 + repeated
                    // fulls): survive ranking and are borrowed unchanged on
                    // every combination; the acquisition guard never fires
                    // for an initialized state.
                    for (name, state) in [
                        ("fast-empty", empty_rows(RingFindType::Fast)),
                        ("fast-full", typed_full(RingFindType::Fast)),
                        ("sssr-empty", empty_rows(RingFindType::Sssr)),
                        ("sssr-full", typed_full(RingFindType::Sssr)),
                        ("summ-empty", empty_rows(RingFindType::SymmSssr)),
                        ("summ-full", typed_full(RingFindType::SymmSssr)),
                    ] {
                        run_case(
                            &graph,
                            &valence,
                            name,
                            Some(&state),
                            canonical,
                            acquisition,
                            Expect::PreserveNone,
                        );
                    }
                }
            }
        }
    }

    mod kekulize_ring_k04_engine {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene() -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                    .collect(),
            )
        }

        fn full_rows(find_type: RingFindType) -> RingInfo {
            let mut rows = ring_info_from_selected_rows(
                6,
                6,
                &[vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                    AtomId::new(4),
                    AtomId::new(5),
                ]],
                &[vec![
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(5),
                ]],
            )
            .unwrap();
            rows.initialize(find_type);
            rows
        }

        fn empty_rows(find_type: RingFindType) -> RingInfo {
            RingInfo::new(find_type, 6, 6)
        }

        fn reset_state() -> RingInfo {
            let mut rings = RingInfo::new(RingFindType::OtherOrUnknown, 6, 6);
            rings.reset();
            rings
        }

        #[derive(Clone, Copy, PartialEq, Eq, Debug)]
        enum Update {
            None,
            SssrRows,
            ResetState,
        }

        #[derive(Clone, Copy, PartialEq, Eq, Debug)]
        enum Outcome {
            /// Topology unchanged, final_valence None (early return).
            EarlyReturn,
            /// Topology unchanged, final_valence Some (aromatic, no fused
            /// call, no marking or marking cleared nothing observable in
            /// bond state; atoms keep aromatic flags).
            UnchangedSomeValence,
            /// Aromatic atom flags cleared, bonds untouched.
            AtomsCleared,
            /// Kekulized: doubles {b0,b2,b4}, singles {b1,b3,b5}; mark
            /// additionally clears aromatic atom and bond flags.
            Kekulized,
            /// AromaticAtomOutsideRing{atom:0} from initialized-empty rows.
            OutsideRing,
        }

        // The complete frozen K01 160-call product: 10 ring inputs x atoms
        // mask {none,all} x bonds mask {none,all} x mark {F,T} x canonical
        // {F,T}. Literal expectations are frozen from the source-derived K01
        // table before this test ran; observing different values is a
        // contradiction for ROOT review, never a value substitution.
        #[test]
        fn kekulize_ring_k04_selected_fragment_160_call_product() {
            let graph = benzene();
            let graph_snapshot = graph.clone();
            let mut calls = 0usize;
            let ring_inputs: [(&str, Option<RingInfo>); 10] = [
                ("absent", None),
                ("reset", Some(reset_state())),
                (
                    "other-empty",
                    Some(empty_rows(RingFindType::OtherOrUnknown)),
                ),
                ("other-full", Some(full_rows(RingFindType::OtherOrUnknown))),
                ("fast-empty", Some(empty_rows(RingFindType::Fast))),
                ("fast-full", Some(full_rows(RingFindType::Fast))),
                ("sssr-empty", Some(empty_rows(RingFindType::Sssr))),
                ("sssr-full", Some(full_rows(RingFindType::Sssr))),
                ("symm-empty", Some(empty_rows(RingFindType::SymmSssr))),
                ("symm-full", Some(full_rows(RingFindType::SymmSssr))),
            ];
            for (ring_name, ring_input) in &ring_inputs {
                for atoms_none in [true, false] {
                    for bonds_none in [true, false] {
                        for mark in [false, true] {
                            for canonical in [false, true] {
                                calls += 1;
                                let label = format!(
                                    "{ring_name}/atoms={}/bonds={}/mark={mark}/canon={canonical}",
                                    if atoms_none { "none" } else { "all" },
                                    if bonds_none { "none" } else { "all" },
                                );
                                let atoms_mask = [!atoms_none; 6];
                                let bonds_mask = [!bonds_none; 6];
                                let params = KekulizeParams {
                                    mark_atoms_bonds: mark,
                                    canonical,
                                    ..KekulizeParams::default()
                                };
                                let supplied_snapshot = ring_input.clone();
                                let rank_base = ring_transport_probe::rank_calls();
                                let sssr_base = ring_transport_probe::sssr_calls();
                                let borrow_base = ring_transport_probe::borrow_sites().len();
                                let result = kekulize_fragment(
                                    &graph,
                                    &atoms_mask,
                                    &bonds_mask,
                                    &params,
                                    None,
                                    ring_input.as_ref(),
                                    None,
                                );
                                // Input immutability on success and error.
                                assert_eq!(&graph, &graph_snapshot, "{label}: graph mutated");
                                assert_eq!(
                                    ring_input, &supplied_snapshot,
                                    "{label}: supplied ring mutated"
                                );

                                // Frozen literal expectation for this combo.
                                let (outcome, update, rank_delta, sssr_delta, borrow_delta) =
                                    expectation(ring_name, atoms_none, bonds_none, mark, canonical);
                                assert_eq!(
                                    ring_transport_probe::rank_calls() - rank_base,
                                    rank_delta,
                                    "{label}: rank delta"
                                );
                                assert_eq!(
                                    ring_transport_probe::sssr_calls() - sssr_base,
                                    sssr_delta,
                                    "{label}: SSSR delta"
                                );
                                assert_eq!(
                                    ring_transport_probe::borrow_sites().len() - borrow_base,
                                    borrow_delta,
                                    "{label}: borrow delta"
                                );
                                match outcome {
                                    Outcome::OutsideRing => {
                                        assert_eq!(
                                            result.err(),
                                            Some(KekulizeError::AromaticAtomOutsideRing {
                                                atom: AtomId::new(0)
                                            }),
                                            "{label}: error category"
                                        );
                                        continue;
                                    }
                                    _ => {}
                                }
                                let assignment =
                                    result.unwrap_or_else(|e| panic!("{label}: {e:?}"));
                                match outcome {
                                    Outcome::EarlyReturn => {
                                        assert_eq!(assignment.topology, graph, "{label}");
                                        assert!(assignment.final_valence.is_none(), "{label}");
                                        assert!(assignment.refreshed_valence_atoms.is_empty());
                                    }
                                    Outcome::UnchangedSomeValence => {
                                        assert_eq!(assignment.topology, graph, "{label}");
                                        assert!(assignment.final_valence.is_some(), "{label}");
                                        assert!(assignment.refreshed_valence_atoms.is_empty());
                                    }
                                    Outcome::AtomsCleared => {
                                        assert!(assignment.final_valence.is_some(), "{label}");
                                        for atom in &assignment.topology.atoms {
                                            assert!(!atom.is_aromatic(), "{label}: atom flag");
                                        }
                                        for bond in &assignment.topology.bonds {
                                            assert_eq!(bond.order(), BondOrder::Aromatic);
                                            assert!(bond.is_aromatic(), "{label}: bond flag");
                                        }
                                    }
                                    Outcome::Kekulized => {
                                        assert!(assignment.final_valence.is_some(), "{label}");
                                        for (index, bond) in
                                            assignment.topology.bonds.iter().enumerate()
                                        {
                                            // Source-derived from the two
                                            // audited owners: the SSSR row
                                            // is atoms [0,5,4,3,2,1] bonds
                                            // [b5..b0]; the matching starts
                                            // at the minimum-rank atom and
                                            // steps to the lower-ranked
                                            // neighbor. Noncanonical iota
                                            // ranks [0..5] start at atom 0
                                            // and step to atom 1, doubling
                                            // b0 -> {b0,b2,b4}; canonical
                                            // ranks [3,1,0,2,4,5] start at
                                            // atom 2 and step to atom 1,
                                            // doubling b1 -> {b1,b3,b5}.
                                            let double_at_even = !canonical;
                                            let expected = if (index % 2 == 0) == double_at_even {
                                                BondOrder::Double
                                            } else {
                                                BondOrder::Single
                                            };
                                            assert_eq!(
                                                bond.order(),
                                                expected,
                                                "{label}: bond {index} order"
                                            );
                                            assert_eq!(
                                                bond.is_aromatic(),
                                                !mark,
                                                "{label}: bond {index} flag"
                                            );
                                        }
                                        for atom in &assignment.topology.atoms {
                                            assert_eq!(atom.is_aromatic(), !mark, "{label}");
                                        }
                                    }
                                    Outcome::OutsideRing => unreachable!(),
                                }
                                match update {
                                    Update::None => {
                                        assert!(assignment.ring_update.is_none(), "{label}")
                                    }
                                    Update::SssrRows => {
                                        let value = assignment
                                            .ring_update
                                            .as_ref()
                                            .unwrap_or_else(|| panic!("{label}: no update"));
                                        assert!(value.is_initialized(), "{label}");
                                        assert_eq!(value.find_type(), RingFindType::Sssr);
                                        assert_eq!(value.atom_rings().len(), 1, "{label}");
                                        assert_eq!(value.atom_rings()[0].len(), 6, "{label}");
                                    }
                                    Update::ResetState => {
                                        let value = assignment
                                            .ring_update
                                            .as_ref()
                                            .unwrap_or_else(|| panic!("{label}: no update"));
                                        assert!(!value.is_initialized(), "{label}");
                                    }
                                }
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 160, "exact census");
        }

        fn expectation(
            ring_name: &str,
            atoms_none: bool,
            bonds_none: bool,
            mark: bool,
            canonical: bool,
        ) -> (Outcome, Update, u64, u64, usize) {
            if atoms_none {
                // Early return before any ring work for every combination.
                return (Outcome::EarlyReturn, Update::None, 0, 0, 0);
            }
            let typed = matches!(
                ring_name,
                "fast-empty"
                    | "fast-full"
                    | "sssr-empty"
                    | "sssr-full"
                    | "symm-empty"
                    | "symm-full"
            );
            let other = ring_name.starts_with("other");
            let rank_delta = u64::from(canonical);
            if bonds_none {
                // No candidate acquisition; the marking guard is the only
                // possible acquisition (mark=true with uninitialized rows).
                let acquires =
                    mark && (ring_name == "absent" || ring_name == "reset" || (other && canonical));
                let update = if acquires {
                    Update::SssrRows
                } else if other && canonical {
                    Update::ResetState
                } else {
                    Update::None
                };
                let sssr_delta = u64::from(acquires);
                let outcome = if mark {
                    Outcome::AtomsCleared
                } else {
                    Outcome::UnchangedSomeValence
                };
                return (outcome, update, rank_delta, sssr_delta, 0);
            }
            // bonds = all: the candidate acquisition guard decides.
            let empty = ring_name.ends_with("-empty");
            if typed {
                if empty {
                    // Initialized-empty rows survive ranking, produce no
                    // candidates and no fused call, but the candidate-read
                    // site still consumes the borrowed rows once; marking
                    // with all/all masks raises the source exception on the
                    // empty membership.
                    let outcome = if mark {
                        Outcome::OutsideRing
                    } else {
                        Outcome::UnchangedSomeValence
                    };
                    return (outcome, Update::None, rank_delta, 0, 1);
                }
                // Initialized full rows are borrowed; one borrow-site read
                // feeds the candidates and the fused dispatch.
                return (Outcome::Kekulized, Update::None, rank_delta, 0, 1);
            }
            if other && !canonical && empty {
                let outcome = if mark {
                    Outcome::OutsideRing
                } else {
                    Outcome::UnchangedSomeValence
                };
                return (outcome, Update::None, 0, 0, 1);
            }
            if other && !canonical {
                // Initialized full Other rows borrowed unchanged.
                return (Outcome::Kekulized, Update::None, 0, 0, 1);
            }
            // absent/reset, any Other under canonical (post-ranking reset),
            // and canonical Other-empty: acquisition at the candidate guard,
            // fresh SSSR rows, fused dispatch, marking reuse.
            (Outcome::Kekulized, Update::SssrRows, rank_delta, 1, 1)
        }
    }

    mod kekulize_ring_k05_query_bonds {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene() -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                    .collect(),
            )
        }

        fn other_full_rows() -> RingInfo {
            let mut rows = ring_info_from_selected_rows(
                6,
                6,
                &[vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                    AtomId::new(4),
                    AtomId::new(5),
                ]],
                &[vec![
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(5),
                ]],
            )
            .unwrap();
            rows.initialize(RingFindType::OtherOrUnknown);
            rows
        }

        fn query_rows(graph: &TopologyBlock, explicit: bool) -> (Vec<QueryAtom>, Vec<QueryBond>) {
            let query_atoms = graph
                .atoms
                .iter()
                .map(|carrier| {
                    QueryAtom::from_carrier_parts(
                        carrier.clone(),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect();
            let query_bonds = graph
                .bonds
                .iter()
                .map(|carrier| {
                    let predicate =
                        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic));
                    if explicit {
                        QueryBond::from_parts(carrier.clone(), predicate)
                    } else {
                        QueryBond::from_carrier_parts(carrier.clone(), predicate)
                    }
                })
                .collect();
            (query_atoms, query_bonds)
        }

        // Eight calls: explicit bond-type predicates remove ALL effective
        // bonds (bondsToUse.none()); carrier-derived predicates do not.
        // Initialized Other FULL rows distinguish reset from preserve and
        // acquisition. Original masks and query values stay unchanged.
        #[test]
        fn kekulize_ring_k05_query_filtered_effective_bonds_none() {
            let graph = benzene();
            let graph_snapshot = graph.clone();
            let rows = other_full_rows();
            let mut calls = 0usize;
            for explicit in [true, false] {
                for mark in [false, true] {
                    for canonical in [false, true] {
                        calls += 1;
                        let label = format!("explicit={explicit}/mark={mark}/canon={canonical}");
                        let (query_atoms, query_bonds) = query_rows(&graph, explicit);
                        let state =
                            QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &graph)
                                .unwrap();
                        let atoms_mask = [true; 6];
                        let bonds_mask = [true; 6];
                        let params = KekulizeParams {
                            mark_atoms_bonds: mark,
                            canonical,
                            ..KekulizeParams::default()
                        };
                        let rows_snapshot = rows.clone();
                        let rank_base = ring_transport_probe::rank_calls();
                        let sssr_base = ring_transport_probe::sssr_calls();
                        let result = kekulize_fragment(
                            &graph,
                            &atoms_mask,
                            &bonds_mask,
                            &params,
                            Some(state),
                            Some(&rows),
                            None,
                        );
                        assert_eq!(&graph, &graph_snapshot, "{label}: graph mutated");
                        assert_eq!(rows, rows_snapshot, "{label}: rows mutated");
                        assert_eq!(query_atoms.len(), 6, "{label}: query rows");
                        assert_eq!(query_bonds.len(), 6, "{label}: query rows");

                        let assignment =
                            result.unwrap_or_else(|error| panic!("{label}: {error:?}"));
                        if explicit {
                            // Effective bonds none: no candidate acquisition;
                            // the marking guard is the only possible find.
                            assert_eq!(
                                ring_transport_probe::rank_calls() - rank_base,
                                u64::from(canonical),
                                "{label}: rank delta"
                            );
                            match (mark, canonical) {
                                (false, false) => {
                                    assert!(assignment.ring_update.is_none(), "{label}");
                                    assert_eq!(assignment.topology, graph, "{label}");
                                }
                                (false, true) => {
                                    let update = assignment.ring_update.unwrap();
                                    assert!(!update.is_initialized(), "{label}: reset");
                                    assert_eq!(assignment.topology, graph, "{label}");
                                }
                                (true, false) => {
                                    // Initialized Other survives and is read
                                    // by nothing: preserved untouched.
                                    assert!(assignment.ring_update.is_none(), "{label}");
                                    for atom in &assignment.topology.atoms {
                                        assert!(!atom.is_aromatic(), "{label}");
                                    }
                                    for bond in &assignment.topology.bonds {
                                        assert_eq!(bond.order(), BondOrder::Aromatic);
                                        assert!(bond.is_aromatic(), "{label}");
                                    }
                                }
                                (true, true) => {
                                    // Ranking reset the Other state, then the
                                    // marking guard acquired fresh SSSR rows.
                                    let update = assignment.ring_update.unwrap();
                                    assert!(update.is_initialized(), "{label}");
                                    assert_eq!(update.find_type(), RingFindType::Sssr);
                                    assert_eq!(update.atom_rings().len(), 1, "{label}");
                                    for atom in &assignment.topology.atoms {
                                        assert!(!atom.is_aromatic(), "{label}");
                                    }
                                }
                            }
                            assert_eq!(
                                ring_transport_probe::sssr_calls() - sssr_base,
                                u64::from(mark && canonical),
                                "{label}: SSSR delta"
                            );
                        } else {
                            // Carrier-derived predicates are ordinary bonds:
                            // full Other-row behavior with a fused dispatch.
                            assert_eq!(
                                ring_transport_probe::rank_calls() - rank_base,
                                u64::from(canonical),
                                "{label}: rank delta"
                            );
                            assert_eq!(
                                ring_transport_probe::sssr_calls() - sssr_base,
                                u64::from(canonical),
                                "{label}: SSSR delta"
                            );
                            for (index, bond) in assignment.topology.bonds.iter().enumerate() {
                                let double_at_even = !canonical;
                                let expected = if (index % 2 == 0) == double_at_even {
                                    BondOrder::Double
                                } else {
                                    BondOrder::Single
                                };
                                assert_eq!(bond.order(), expected, "{label}: b{index}");
                                assert_eq!(bond.is_aromatic(), !mark, "{label}");
                            }
                            match canonical {
                                false => assert!(assignment.ring_update.is_none(), "{label}"),
                                true => {
                                    let update = assignment.ring_update.unwrap();
                                    assert!(update.is_initialized(), "{label}");
                                    assert_eq!(update.find_type(), RingFindType::Sssr);
                                }
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 8, "exact census");
        }
    }

    mod kekulize_ring_k06_failures {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene() -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                    .collect(),
            )
        }

        fn reset_state() -> RingInfo {
            let mut rings = RingInfo::new(RingFindType::OtherOrUnknown, 6, 6);
            rings.reset();
            rings
        }

        fn rows_with_dimensions(atoms: usize, bonds: usize) -> RingInfo {
            ring_info_from_selected_rows(atoms, bonds, &[], &[]).unwrap()
        }

        #[test]
        fn kekulize_ring_k06_typed_failure_table() {
            // Calls 1-2: empty graph, both canonical settings. The
            // atoms-none early return precedes every ring computation.
            let empty = topology(vec![], vec![]);
            for canonical in [false, true] {
                let rank_base = ring_transport_probe::rank_calls();
                let sssr_base = ring_transport_probe::sssr_calls();
                let assignment = kekulize_fragment(
                    &empty,
                    &[],
                    &[],
                    &KekulizeParams {
                        canonical,
                        ..KekulizeParams::default()
                    },
                    None,
                    None,
                    None,
                )
                .unwrap();
                assert_eq!(assignment.topology, empty);
                assert!(assignment.final_valence.is_none());
                assert!(assignment.ring_update.is_none());
                assert_eq!(ring_transport_probe::rank_calls(), rank_base);
                assert_eq!(ring_transport_probe::sssr_calls(), sssr_base);
            }

            // Calls 3-4: nonaromatic acyclic two-carbon graph with an
            // initialized EMPTY Fast assignment supplied. The no-aromatic
            // early return runs AFTER selection/valence work but never
            // validates, fetches or mutates the unused assignment; the
            // initialized-empty quality predicates hold directly.
            let acyclic = topology(
                vec![
                    atom(0, AtomSpec::new(Element::C)),
                    atom(1, AtomSpec::new(Element::C)),
                ],
                vec![bond(0, 0, 1)],
            );
            let supplied = RingInfo::new(RingFindType::Fast, 2, 1);
            assert!(supplied.is_initialized());
            assert!(supplied.is_find_fast_or_better());
            assert_eq!(supplied.atom_rings().len(), 0);
            let supplied_snapshot = supplied.clone();
            for canonical in [false, true] {
                let assignment = kekulize_fragment(
                    &acyclic,
                    &[true, true],
                    &[true],
                    &KekulizeParams {
                        canonical,
                        ..KekulizeParams::default()
                    },
                    None,
                    Some(&supplied),
                    None,
                )
                .unwrap();
                assert_eq!(assignment.topology, acyclic);
                // KekulizeFragment calculates selected I before !foundAromatic.
                assert_eq!(
                    assignment.final_valence,
                    Some(ValenceAssignment {
                        explicit_valence: vec![1, 1],
                        implicit_hydrogens: vec![3, 3],
                    })
                );
                assert!(assignment.ring_update.is_none());
            }
            assert_eq!(supplied, supplied_snapshot);

            // Calls 5-6: aromatic graph with atoms mask none and
            // deliberately mismatched initialized rows (7 atom dimensions).
            // The early return precedes any validation of the unused rows.
            let graph = benzene();
            let mismatched = rows_with_dimensions(7, 6);
            let mismatched_snapshot = mismatched.clone();
            assert_eq!(mismatched.atom_row_count(), 7);
            for canonical in [false, true] {
                let assignment = kekulize_fragment(
                    &graph,
                    &[false; 6],
                    &[true; 6],
                    &KekulizeParams {
                        canonical,
                        ..KekulizeParams::default()
                    },
                    None,
                    Some(&mismatched),
                    None,
                )
                .unwrap();
                assert_eq!(assignment.topology, graph);
                assert!(assignment.ring_update.is_none());
            }
            assert_eq!(mismatched, mismatched_snapshot);

            // Call 7: wrong atom mask length.
            assert_eq!(
                kekulize_fragment(
                    &graph,
                    &[true; 5],
                    &[true; 6],
                    &KekulizeParams::default(),
                    None,
                    None,
                    None
                ),
                Err(KekulizeError::AtomSelectionLength {
                    expected: 6,
                    actual: 5,
                })
            );

            // Call 8: wrong bond mask length.
            assert_eq!(
                kekulize_fragment(
                    &graph,
                    &[true; 6],
                    &[true; 5],
                    &KekulizeParams::default(),
                    None,
                    None,
                    None
                ),
                Err(KekulizeError::BondSelectionLength {
                    expected: 6,
                    actual: 5,
                })
            );

            // Call 9: malformed query state (rows built for a five-atom
            // topology) rejected before any ring work.
            let five_atom_graph = topology(
                (0..5)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..5)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 5))
                    .collect(),
            );
            let query_atoms = five_atom_graph
                .atoms
                .iter()
                .map(|carrier| {
                    QueryAtom::from_carrier_parts(
                        carrier.clone(),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect::<Vec<_>>();
            let query_bonds = five_atom_graph
                .bonds
                .iter()
                .map(|carrier| {
                    QueryBond::from_carrier_parts(
                        carrier.clone(),
                        QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
                    )
                })
                .collect::<Vec<_>>();
            let foreign_state =
                QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &five_atom_graph)
                    .unwrap();
            let result = kekulize_with_query_state_and_ring_info(
                &graph,
                &KekulizeParams::default(),
                Some(foreign_state),
                None,
                None,
            );
            assert!(matches!(
                result,
                Err(KekulizeError::InvalidQueryState(
                    QueryStateError::AtomCount { .. }
                ))
            ));

            // Call 10: initialized wrong ATOM dimensions consumed by the
            // canonical ranking helper.
            let wrong_atoms = rows_with_dimensions(7, 6);
            let wrong_atoms_snapshot = wrong_atoms.clone();
            assert_eq!(
                kekulize_fragment(
                    &graph,
                    &[true; 6],
                    &[true; 6],
                    &KekulizeParams {
                        canonical: true,
                        ..KekulizeParams::default()
                    },
                    None,
                    Some(&wrong_atoms),
                    None
                ),
                Err(KekulizeError::CanonicalRank(
                    CanonicalRankError::PreparedRingLength {
                        expected_atoms: 6,
                        actual_atoms: 7,
                        expected_bonds: 6,
                        actual_bonds: 6,
                    }
                ))
            );
            assert_eq!(wrong_atoms, wrong_atoms_snapshot);

            // Call 11: initialized wrong BOND dimensions consumed by the
            // NONCANONICAL candidate read (ranking skipped).
            let wrong_bonds = rows_with_dimensions(6, 5);
            let wrong_bonds_snapshot = wrong_bonds.clone();
            assert_eq!(
                kekulize_fragment(
                    &graph,
                    &[true; 6],
                    &[true; 6],
                    &KekulizeParams {
                        canonical: false,
                        ..KekulizeParams::default()
                    },
                    None,
                    Some(&wrong_bonds),
                    None
                ),
                Err(KekulizeError::CanonicalRank(
                    CanonicalRankError::PreparedRingLength {
                        expected_atoms: 6,
                        actual_atoms: 6,
                        expected_bonds: 6,
                        actual_bonds: 5,
                    }
                ))
            );
            assert_eq!(wrong_bonds, wrong_bonds_snapshot);

            // Call 12: reset input with an aromatic full selection is a
            // valid reacquisition: fresh SSSR rows, noncanonical doubles at
            // {b0,b2,b4}, all aromatic flags cleared under mark=true.
            let reset = reset_state();
            let reset_snapshot = reset.clone();
            let assignment = kekulize_fragment(
                &graph,
                &[true; 6],
                &[true; 6],
                &KekulizeParams {
                    mark_atoms_bonds: true,
                    canonical: false,
                    ..KekulizeParams::default()
                },
                None,
                Some(&reset),
                None,
            )
            .unwrap();
            assert_eq!(reset, reset_snapshot);
            let update = assignment.ring_update.unwrap();
            assert!(update.is_initialized());
            assert_eq!(update.find_type(), RingFindType::Sssr);
            for (index, bond) in assignment.topology.bonds.iter().enumerate() {
                let expected = if index % 2 == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                };
                assert_eq!(bond.order(), expected);
                assert!(!bond.is_aromatic());
            }
            for atom in &assignment.topology.atoms {
                assert!(!atom.is_aromatic());
            }
        }
    }

    mod kekulize_ring_k07_partial_fragment {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn two_cycles() -> TopologyBlock {
            let atoms = (0..12)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect::<Vec<_>>();
            let bonds = (0..12)
                .map(|id| {
                    let (begin, end) = if id < 6 {
                        (id, (id + 1) % 6)
                    } else {
                        (id, 6 + (id - 6 + 1) % 6)
                    };
                    aromatic_bond(id, begin, end)
                })
                .collect::<Vec<_>>();
            topology(atoms, bonds)
        }

        fn both_cycle_rows(reversed: bool) -> RingInfo {
            let mut first = vec![
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(4),
                AtomId::new(5),
            ];
            let mut first_bonds = vec![
                BondId::new(0),
                BondId::new(1),
                BondId::new(2),
                BondId::new(3),
                BondId::new(4),
                BondId::new(5),
            ];
            let mut second = vec![
                AtomId::new(6),
                AtomId::new(7),
                AtomId::new(8),
                AtomId::new(9),
                AtomId::new(10),
                AtomId::new(11),
            ];
            let mut second_bonds = vec![
                BondId::new(6),
                BondId::new(7),
                BondId::new(8),
                BondId::new(9),
                BondId::new(10),
                BondId::new(11),
            ];
            if reversed {
                first.reverse();
                first_bonds.reverse();
                second.reverse();
                second_bonds.reverse();
            }
            let mut rows = ring_info_from_selected_rows(
                12,
                12,
                &[first, second],
                &[first_bonds, second_bonds],
            )
            .unwrap();
            rows.initialize(RingFindType::OtherOrUnknown);
            rows
        }

        // Eight calls: forward/reversed full-cycle supplied rows x mark x
        // canonical on a disconnected two-cycle topology selecting only the
        // first cycle. Excluded second-cycle atoms/bonds keep full snapshots;
        // every supplied row, membership, type and the initialized flag stay
        // unchanged; canonical=false proves preservation by the ACTUAL
        // borrowed reference identity at the fused consumer.
        #[test]
        fn kekulize_ring_k07_borrowed_row_preservation_on_partial_fragment() {
            let graph = two_cycles();
            let graph_snapshot = graph.clone();
            let mut atoms_mask = [false; 12];
            for selected in &mut atoms_mask[..6] {
                *selected = true;
            }
            let mut bonds_mask = [false; 12];
            for selected in &mut bonds_mask[..6] {
                *selected = true;
            }
            let mut calls = 0usize;
            for reversed in [false, true] {
                let rows = both_cycle_rows(reversed);
                for mark in [false, true] {
                    for canonical in [false, true] {
                        calls += 1;
                        let label = format!("reversed={reversed}/mark={mark}/canon={canonical}");
                        let rows_snapshot = rows.clone();
                        let borrow_base = ring_transport_probe::borrow_sites().len();
                        let assignment = kekulize_fragment(
                            &graph,
                            &atoms_mask,
                            &bonds_mask,
                            &KekulizeParams {
                                mark_atoms_bonds: mark,
                                canonical,
                                ..KekulizeParams::default()
                            },
                            None,
                            Some(&rows),
                            None,
                        )
                        .unwrap();
                        assert_eq!(&graph, &graph_snapshot, "{label}: graph mutated");
                        assert_eq!(rows, rows_snapshot, "{label}: rows mutated");
                        assert_eq!(rows.is_initialized(), true);
                        assert_eq!(rows.find_type(), RingFindType::OtherOrUnknown);
                        assert_eq!(rows.atom_rings().len(), 2, "{label}: row count");

                        // Excluded second cycle: full snapshots unchanged.
                        for atom in &assignment.topology.atoms[6..] {
                            assert!(atom.is_aromatic(), "{label}: excluded atom");
                        }
                        for bond in &assignment.topology.bonds[6..] {
                            assert_eq!(bond.order(), BondOrder::Aromatic, "{label}");
                            assert!(bond.is_aromatic(), "{label}: excluded bond");
                        }
                        // Selected first cycle kekulizes: exactly three
                        // alternating doubles (any valid Kekulé form of the
                        // six-cycle), flags follow mark.
                        let mut doubles = 0;
                        for (index, bond) in assignment.topology.bonds[..6].iter().enumerate() {
                            match bond.order() {
                                BondOrder::Double => doubles += 1,
                                BondOrder::Single => {}
                                other => panic!("{label}: bond {index} order {other:?}"),
                            }
                            assert_eq!(bond.is_aromatic(), !mark, "{label}: flag");
                        }
                        assert_eq!(doubles, 3, "{label}: double count");
                        for (index, window) in assignment.topology.bonds[..6].windows(2).enumerate()
                        {
                            assert_ne!(
                                window[0].order(),
                                window[1].order(),
                                "{label}: alternation at {index}"
                            );
                        }
                        assert_ne!(
                            assignment.topology.bonds[0].order(),
                            assignment.topology.bonds[5].order(),
                            "{label}: wrap alternation"
                        );
                        for atom in &assignment.topology.atoms[..6] {
                            assert_eq!(atom.is_aromatic(), !mark, "{label}");
                        }

                        match canonical {
                            false => {
                                // Preserved untouched: ring_update None and
                                // the fused consumer saw the ACTUAL borrowed
                                // reference, not a clone.
                                assert!(assignment.ring_update.is_none(), "{label}");
                                let sites = ring_transport_probe::borrow_sites();
                                assert_eq!(sites.len(), borrow_base + 1, "{label}");
                                assert_eq!(
                                    sites[borrow_base], &rows as *const RingInfo as usize,
                                    "{label}: borrowed identity"
                                );
                            }
                            true => {
                                // Canonical ranking reset the Other state;
                                // the candidate guard reacquired SSSR rows
                                // (both cycles found, first selected).
                                let update = assignment.ring_update.unwrap();
                                assert!(update.is_initialized(), "{label}");
                                assert_eq!(update.find_type(), RingFindType::Sssr);
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 8, "exact census");
        }
    }

    mod kekulize_ring_move_proof {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn benzene() -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| aromatic_bond(id, id, (id + 1) % 6))
                    .collect(),
            )
        }

        fn other_full_rows() -> RingInfo {
            let mut rows = ring_info_from_selected_rows(
                6,
                6,
                &[vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                    AtomId::new(4),
                    AtomId::new(5),
                ]],
                &[vec![
                    BondId::new(0),
                    BondId::new(1),
                    BondId::new(2),
                    BondId::new(3),
                    BondId::new(4),
                    BondId::new(5),
                ]],
            )
            .unwrap();
            rows.initialize(RingFindType::OtherOrUnknown);
            rows
        }

        #[derive(Clone, Copy, PartialEq, Eq, Debug)]
        enum Expect {
            /// No acquisition: ring_update is None (absent input preserved).
            None,
            /// No acquisition, but the canonical ranking reset the supplied
            /// Other state: ring_update is Some(uninitialized).
            ResetOnly,
            /// Exactly one acquisition moved into ring_update.
            MovedOnce,
        }

        // The exact MOVE-frozen eight-call product: {None, Other-full} x
        // {bonds-none, bonds-all} x {mark false,true} with canonical=true.
        // The observation at the real findSSSR completion captures the
        // acquired buffer's outer row pointers; for every Some(initialized)
        // the returned update carries the SAME buffers (moved, not copied).
        #[test]
        fn kekulize_ring_move_acquired_buffer_identity() {
            let graph = benzene();
            let graph_snapshot = graph.clone();
            let other_full = other_full_rows();
            let mut calls = 0usize;
            for supplied in ["absent", "other-full"] {
                for bonds_none in [true, false] {
                    for mark in [false, true] {
                        calls += 1;
                        let label = format!(
                            "{supplied}/bonds={}/mark={mark}",
                            if bonds_none { "none" } else { "all" }
                        );
                        let rings: Option<&RingInfo> = match supplied {
                            "absent" => None,
                            _ => Some(&other_full),
                        };
                        let supplied_snapshot = rings.cloned();
                        let expect = match (supplied, bonds_none, mark) {
                            ("absent", true, false) => Expect::None,
                            ("other-full", true, false) => Expect::ResetOnly,
                            _ => Expect::MovedOnce,
                        };
                        let bonds_mask = [!bonds_none; 6];
                        let base = ring_transport_probe::acquired_buffers().len();
                        let assignment = kekulize_fragment(
                            &graph,
                            &[true; 6],
                            &bonds_mask,
                            &KekulizeParams {
                                mark_atoms_bonds: mark,
                                canonical: true,
                                ..KekulizeParams::default()
                            },
                            None,
                            rings,
                            None,
                        )
                        .unwrap_or_else(|error| panic!("{label}: {error:?}"));
                        // Immutability on every call.
                        assert_eq!(&graph, &graph_snapshot, "{label}: graph mutated");
                        assert_eq!(
                            rings.cloned(),
                            supplied_snapshot,
                            "{label}: supplied mutated"
                        );
                        let buffers = ring_transport_probe::acquired_buffers();
                        let delta = buffers.len() - base;
                        match expect {
                            Expect::None => {
                                assert_eq!(delta, 0, "{label}: unexpected acquisition");
                                assert!(assignment.ring_update.is_none(), "{label}");
                            }
                            Expect::ResetOnly => {
                                assert_eq!(delta, 0, "{label}: unexpected acquisition");
                                let update = assignment
                                    .ring_update
                                    .as_ref()
                                    .unwrap_or_else(|| panic!("{label}: missing reset"));
                                assert!(!update.is_initialized(), "{label}: reset state");
                            }
                            Expect::MovedOnce => {
                                assert_eq!(delta, 1, "{label}: acquisition count");
                                let update = assignment
                                    .ring_update
                                    .as_ref()
                                    .unwrap_or_else(|| panic!("{label}: missing update"));
                                assert!(update.is_initialized(), "{label}");
                                assert_eq!(update.find_type(), RingFindType::Sssr, "{label}");
                                // The nonempty fixture buffers prove genuine
                                // identity: returned outer pointers equal the
                                // captured finder pointers.
                                assert_eq!(
                                    update.atom_rings().as_ptr() as usize,
                                    buffers[base].0,
                                    "{label}: atom buffer identity"
                                );
                                assert_eq!(
                                    update.atom_rings().len(),
                                    buffers[base].1,
                                    "{label}: atom row count"
                                );
                                assert_eq!(
                                    update.bond_rings().as_ptr() as usize,
                                    buffers[base].2,
                                    "{label}: bond buffer identity"
                                );
                                assert_eq!(
                                    update.bond_rings().len(),
                                    buffers[base].3,
                                    "{label}: bond row count"
                                );
                                assert_eq!(buffers[base].1, 1, "{label}: rows nonempty");
                                assert_eq!(buffers[base].3, 1, "{label}: rows nonempty");
                                let literal_atoms: Vec<usize> =
                                    update.atom_rings()[0].iter().map(|a| a.index()).collect();
                                let literal_bonds: Vec<usize> =
                                    update.bond_rings()[0].iter().map(|b| b.index()).collect();
                                assert_eq!(literal_atoms, vec![0, 5, 4, 3, 2, 1], "{label}");
                                assert_eq!(literal_bonds, vec![5, 4, 3, 2, 1, 0], "{label}");
                                for index in 0..6 {
                                    assert_eq!(
                                        update.atom_members(AtomId::new(index)),
                                        &[0],
                                        "{label}: membership"
                                    );
                                    assert_eq!(
                                        update.bond_members(BondId::new(index)),
                                        &[0],
                                        "{label}: membership"
                                    );
                                }
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 8, "exact census");
        }
    }

    mod kekulize_ring_convert_c1 {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn aromatic_six_cycle(direction: BondDirection, wedge_begin: usize) -> TopologyBlock {
            topology(
                (0..6)
                    .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                    .collect(),
                (0..6)
                    .map(|id| {
                        let spec = BondSpec::new(
                            AtomId::new(id),
                            AtomId::new((id + 1) % 6),
                            BondOrder::Aromatic,
                        )
                        .with_aromatic(true);
                        let spec = if id == wedge_begin {
                            spec.with_direction(direction)
                        } else {
                            spec
                        };
                        bond_with_spec(id, spec)
                    })
                    .collect(),
            )
        }

        fn rows(atom_order: &[usize], bond_order: &[usize]) -> RingInfo {
            ring_info_from_selected_rows(
                6,
                6,
                &[atom_order
                    .iter()
                    .map(|i| AtomId::new(*i))
                    .collect::<Vec<_>>()],
                &[bond_order
                    .iter()
                    .map(|i| BondId::new(*i))
                    .collect::<Vec<_>>()],
            )
            .unwrap()
        }

        fn ids(rows: &[Vec<AtomId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|atom| atom.index()).collect())
                .collect()
        }

        fn bond_ids(rows: &[Vec<BondId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|bond| bond.index()).collect())
                .collect()
        }

        fn wedged_from_graph(
            graph: &TopologyBlock,
            wedge_begin: usize,
            selected: &[bool],
        ) -> Vec<bool> {
            // The ACTUAL selected directed graph bond whose begin atom is
            // the wedge begin (b_begin=(begin,begin+1)) marks wedged_atoms.
            let mut wedged = vec![false; graph.atoms.len()];
            for bond in &graph.bonds {
                if selected[bond.id().index()] && bond.begin().index() == wedge_begin {
                    wedged[bond.begin().index()] = true;
                }
            }
            wedged
        }

        // C1: 32 orientation/stored-order/wedge calls against the CURRENT
        // helper signature. Literal expectations are derived ONLY from the
        // pinned source (rotate atom rows; derive bonds from consecutive
        // graph pairs + closing edge) and are INDEPENDENT of the stored
        // bond order and the wedge direction. This test is RED before the
        // production fix whenever stored rows diverge from graph adjacency.
        #[test]
        fn kekulize_ring_convert_orientation_stored_order_and_wedge() {
            let atom_orders: [(&str, [usize; 6]); 2] =
                [("F", [0, 1, 2, 3, 4, 5]), ("R", [5, 4, 3, 2, 1, 0])];
            let stored_orders: [[usize; 6]; 4] = [
                [0, 1, 2, 3, 4, 5],
                [5, 4, 3, 2, 1, 0],
                [2, 3, 4, 5, 0, 1],
                [0, 2, 4, 1, 3, 5],
            ];
            let mut calls = 0usize;
            for (orientation, atom_order) in &atom_orders {
                for stored in &stored_orders {
                    for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                        for wedge_begin in [0usize, 2usize] {
                            let rings = rows(atom_order, stored);
                            let graph = aromatic_six_cycle(direction, wedge_begin);
                            let graph_snapshot = graph.clone();
                            let rings_snapshot = rings.clone();
                            let atoms_selection = [true; 6];
                            let bonds_selection = [true; 6];
                            let dummy_input = [false; 6];
                            let atoms_snapshot = atoms_selection;
                            let bonds_snapshot = bonds_selection;
                            let dummy_snapshot = dummy_input;
                            let label = format!(
                                "{orientation}/stored{stored:?}/{direction:?}/w{wedge_begin}"
                            );
                            // Prerequisite: the ACTUAL graph bond selected as
                            // the wedge carrier carries the current direction
                            // and its begin atom is the current wedge begin.
                            assert_eq!(
                                graph.bonds[wedge_begin].direction(),
                                direction,
                                "{label}: wedge bond direction"
                            );
                            assert_eq!(
                                graph.bonds[wedge_begin].begin().index(),
                                wedge_begin,
                                "{label}: wedge bond begin"
                            );
                            let wedged = wedged_from_graph(&graph, wedge_begin, &atoms_selection);
                            let wedged_snapshot = wedged.clone();
                            assert!(wedged[wedge_begin], "{label}: graph wedge");
                            let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                                &graph,
                                &rings,
                                &atoms_selection,
                                &bonds_selection,
                                &dummy_input,
                                &wedged,
                            )
                            .unwrap();
                            // Census incremented AFTER the real invocation.
                            calls += 1;
                            // Post-invocation input preservation for THIS
                            // call, before any result assertion.
                            assert_eq!(&graph, &graph_snapshot, "{label}: graph mutated");
                            assert_eq!(rings, rings_snapshot, "{label}: input mutated");
                            assert_eq!(
                                atoms_selection, atoms_snapshot,
                                "{label}: atoms selection mutated"
                            );
                            assert_eq!(
                                bonds_selection, bonds_snapshot,
                                "{label}: bonds selection mutated"
                            );
                            assert_eq!(dummy_input, dummy_snapshot, "{label}: dummy mutated");
                            assert_eq!(wedged, wedged_snapshot, "{label}: wedge array mutated");
                            let (expected_atoms, expected_bonds): (&[usize], &[usize]) =
                                match (*orientation, wedge_begin) {
                                    ("F", 0) => (&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 4, 5]),
                                    ("F", 2) => (&[2, 3, 4, 5, 0, 1], &[2, 3, 4, 5, 0, 1]),
                                    ("R", 0) => (&[0, 5, 4, 3, 2, 1], &[5, 4, 3, 2, 1, 0]),
                                    _ => (&[2, 1, 0, 5, 4, 3], &[1, 0, 5, 4, 3, 2]),
                                };
                            assert_eq!(
                                ids(&atom_rows),
                                vec![expected_atoms.to_vec()],
                                "{label}: atoms"
                            );
                            assert_eq!(
                                bond_ids(&bond_rows),
                                vec![expected_bonds.to_vec()],
                                "{label}: derived bonds"
                            );
                        }
                    }
                }
            }
            assert_eq!(calls, 32, "exact census");
        }
    }

    mod kekulize_ring_convert_c234 {
        use super::*;
        use crate::ring_info_from_selected_rows;

        fn two_cycle_graph(wedge_bond: usize, direction: BondDirection) -> TopologyBlock {
            let atoms = (0..12)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect::<Vec<_>>();
            let bonds = (0..12)
                .map(|id| {
                    let (begin, end) = if id < 6 {
                        (id, (id + 1) % 6)
                    } else {
                        (id, 6 + (id - 6 + 1) % 6)
                    };
                    let spec =
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Aromatic)
                            .with_aromatic(true);
                    let spec = if id == wedge_bond {
                        spec.with_direction(direction)
                    } else {
                        spec
                    };
                    bond_with_spec(id, spec)
                })
                .collect::<Vec<_>>();
            topology(atoms, bonds)
        }

        fn ids(rows: &[Vec<AtomId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|atom| atom.index()).collect())
                .collect()
        }

        fn bond_ids(rows: &[Vec<BondId>]) -> Vec<Vec<usize>> {
            rows.iter()
                .map(|row| row.iter().map(|bond| bond.index()).collect())
                .collect()
        }

        fn wedge_array(graph: &TopologyBlock, wedge_begin: usize) -> Vec<bool> {
            let mut wedged = vec![false; graph.atoms.len()];
            for bond in &graph.bonds {
                if matches!(
                    bond.direction(),
                    BondDirection::BeginWedge | BondDirection::BeginDash
                ) && bond.begin().index() == wedge_begin
                {
                    wedged[bond.begin().index()] = true;
                }
            }
            wedged
        }

        // C2: 32 wrong-stored-set selection calls. Row include/exclude is
        // decided by DERIVED graph bonds only; stored bond order never
        // changes the output.
        #[test]
        fn kekulize_ring_convert_selection_derived_bond_masks() {
            let a_fwd = (0..6).map(AtomId::new).collect::<Vec<_>>();
            let a_bonds = (0..6).map(BondId::new).collect::<Vec<_>>();
            let b_fwd = (6..12).map(AtomId::new).collect::<Vec<_>>();
            let b_bonds = (6..12).map(BondId::new).collect::<Vec<_>>();
            let mut a_rev = a_fwd.clone();
            a_rev.reverse();
            let mut b_rev = b_fwd.clone();
            b_rev.reverse();
            let mut a_bonds_rev = a_bonds.clone();
            a_bonds_rev.reverse();
            let mut b_bonds_rev = b_bonds.clone();
            b_bonds_rev.reverse();
            struct Expected {
                atoms: Vec<Vec<usize>>,
                bonds: Vec<Vec<usize>>,
            }
            let f_a = Expected {
                atoms: vec![vec![0, 1, 2, 3, 4, 5]],
                bonds: vec![vec![0, 1, 2, 3, 4, 5]],
            };
            let f_b = Expected {
                atoms: vec![vec![6, 7, 8, 9, 10, 11]],
                bonds: vec![vec![6, 7, 8, 9, 10, 11]],
            };
            let r_a = Expected {
                atoms: vec![vec![5, 4, 3, 2, 1, 0]],
                bonds: vec![vec![4, 3, 2, 1, 0, 5]],
            };
            let r_b = Expected {
                atoms: vec![vec![11, 10, 9, 8, 7, 6]],
                bonds: vec![vec![10, 9, 8, 7, 6, 11]],
            };
            let mut calls = 0usize;
            for rings_ab_first in [true, false] {
                for reverse_rows in [false, true] {
                    for swapped_stored in [false, true] {
                        // Row order preserved from input.
                        let (row_a, row_b) = if reverse_rows {
                            (a_rev.clone(), b_rev.clone())
                        } else {
                            (a_fwd.clone(), b_fwd.clone())
                        };
                        let (stored_a, stored_b) = if swapped_stored {
                            (b_bonds_rev.clone(), a_bonds_rev.clone())
                        } else {
                            (a_bonds_rev.clone(), b_bonds_rev.clone())
                        };
                        let (first_a, first_b, s_first, s_second) = if rings_ab_first {
                            (
                                row_a.clone(),
                                row_b.clone(),
                                stored_a.clone(),
                                stored_b.clone(),
                            )
                        } else {
                            (
                                row_b.clone(),
                                row_a.clone(),
                                stored_b.clone(),
                                stored_a.clone(),
                            )
                        };
                        let rings = ring_info_from_selected_rows(
                            12,
                            12,
                            &[first_a.clone(), first_b.clone()],
                            &[s_first.clone(), s_second.clone()],
                        )
                        .unwrap();
                        for mask in ["onlyA", "onlyB", "all", "none"] {
                            let graph = two_cycle_graph(usize::MAX, BondDirection::None);
                            let graph_snapshot = graph.clone();
                            let rings_snapshot = rings.clone();
                            let atoms_all = [true; 12];
                            let atoms_snapshot = atoms_all;
                            let dummy = [false; 12];
                            let dummy_snapshot = dummy;
                            let wedged = vec![false; 12];
                            let wedged_snapshot = wedged.clone();
                            let mut bonds_mask = [false; 12];
                            match mask {
                                "onlyA" => bonds_mask[..6].fill(true),
                                "onlyB" => bonds_mask[6..].fill(true),
                                "all" => bonds_mask.fill(true),
                                _ => {}
                            }
                            let mask_snapshot = bonds_mask;
                            let label = format!(
                                "ab={rings_ab_first}/rev={reverse_rows}/swap={swapped_stored}/{mask}"
                            );
                            let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                                &graph,
                                &rings,
                                &atoms_all,
                                &bonds_mask,
                                &dummy,
                                &wedged,
                            )
                            .unwrap();
                            // Post-invocation six-input comparison for THIS
                            // call, before any output assertion; census after.
                            calls += 1;
                            assert_eq!(&graph, &graph_snapshot, "{label}: graph");
                            assert_eq!(rings, rings_snapshot, "{label}: rings");
                            assert_eq!(atoms_all, atoms_snapshot, "{label}: atoms");
                            assert_eq!(bonds_mask, mask_snapshot, "{label}: mask");
                            assert_eq!(dummy, dummy_snapshot, "{label}: dummy");
                            assert_eq!(wedged, wedged_snapshot, "{label}: wedged");
                            let expect_a = if reverse_rows { &r_a } else { &f_a };
                            let expect_b = if reverse_rows { &r_b } else { &f_b };
                            match mask {
                                "onlyA" => {
                                    assert_eq!(ids(&atom_rows), expect_a.atoms, "{label}");
                                    assert_eq!(bond_ids(&bond_rows), expect_a.bonds, "{label}");
                                }
                                "onlyB" => {
                                    assert_eq!(ids(&atom_rows), expect_b.atoms, "{label}");
                                    assert_eq!(bond_ids(&bond_rows), expect_b.bonds, "{label}");
                                }
                                "all" => {
                                    if rings_ab_first {
                                        assert_eq!(
                                            ids(&atom_rows),
                                            [expect_a.atoms.clone(), expect_b.atoms.clone()]
                                                .concat()
                                        );
                                        assert_eq!(
                                            bond_ids(&bond_rows),
                                            [expect_a.bonds.clone(), expect_b.bonds.clone()]
                                                .concat()
                                        );
                                    } else {
                                        assert_eq!(
                                            ids(&atom_rows),
                                            [expect_b.atoms.clone(), expect_a.atoms.clone()]
                                                .concat()
                                        );
                                        assert_eq!(
                                            bond_ids(&bond_rows),
                                            [expect_b.bonds.clone(), expect_a.bonds.clone()]
                                                .concat()
                                        );
                                    }
                                }
                                _ => {
                                    assert!(atom_rows.is_empty(), "{label}");
                                    assert!(bond_rows.is_empty(), "{label}");
                                }
                            }
                        }
                    }
                }
            }
            assert_eq!(calls, 32, "exact census");
        }

        // C3: 8 exact typed errors including REAL engine calls. The broken
        // atom row is detected at derived conversion BEFORE any bond filter.
        #[test]
        fn kekulize_ring_convert_errors_missing_edges() {
            let mut calls = 0usize;
            for row in [[0usize, 1, 3], [0usize, 1, 2]] {
                // Frozen stored bond row is the LITERAL [0,1,2] for BOTH
                // broken atom rows (original contract), never derived from
                // the atom row.
                let stored: Vec<BondId> = [0usize, 1, 2].iter().map(|i| BondId::new(*i)).collect();
                let rings = ring_info_from_selected_rows(
                    12,
                    12,
                    &[row.iter().map(|i| AtomId::new(*i)).collect::<Vec<_>>()],
                    &[stored],
                )
                .unwrap();
                let rings_snapshot = rings.clone();
                let expected = match row[2] {
                    3 => (1usize, 3usize),
                    _ => (2usize, 0usize),
                };
                for exclude_b0 in [false, true] {
                    let label = format!("row{row:?}/excl{exclude_b0}");
                    // Route 1: private helper — fresh graph/rings/atom-mask/
                    // bond-mask snapshots and EVERY argument supplied to THIS
                    // route (dummy/wedged) captured immediately before the
                    // call and compared immediately after THIS Result,
                    // including Err, BEFORE any error check or route 2.
                    let graph = two_cycle_graph(usize::MAX, BondDirection::None);
                    let graph_snapshot = graph.clone();
                    let rings_snapshot = rings.clone();
                    assert!(rings.is_initialized());
                    assert_eq!(rings.find_type(), RingFindType::OtherOrUnknown);
                    assert_eq!(rings.atom_row_count(), 12);
                    assert_eq!(rings.bond_row_count(), 12);
                    let atoms_mask = [true; 12];
                    let atoms_snapshot = atoms_mask;
                    let mut bonds_mask = [true; 12];
                    if exclude_b0 {
                        bonds_mask[0] = false;
                    }
                    let mask_snapshot = bonds_mask;
                    let dummy_input = [false; 12];
                    let dummy_snapshot = dummy_input;
                    let wedged = vec![false; 12];
                    let wedged_snapshot = wedged.clone();
                    let helper = collect_kekulize_candidate_rings(
                        &graph,
                        &rings,
                        &atoms_mask,
                        &bonds_mask,
                        &dummy_input,
                        &wedged,
                    );
                    calls += 1;
                    assert_eq!(&graph, &graph_snapshot, "{label}: helper graph");
                    assert_eq!(rings, rings_snapshot, "{label}: helper rings");
                    assert_eq!(atoms_mask, atoms_snapshot, "{label}: helper atoms");
                    assert_eq!(bonds_mask, mask_snapshot, "{label}: helper mask");
                    assert_eq!(dummy_input, dummy_snapshot, "{label}: helper dummy");
                    assert_eq!(wedged, wedged_snapshot, "{label}: helper wedged");
                    assert_eq!(
                        helper.err(),
                        Some(KekulizeError::RingFinding(
                            crate::RingFindingError::ExpectedBondNotFound {
                                begin: AtomId::new(expected.0),
                                end: AtomId::new(expected.1),
                            }
                        )),
                        "{label}: helper"
                    );
                    // Route 2: REAL engine with supplied initialized Other —
                    // fresh per-route baselines including the params value.
                    let graph = two_cycle_graph(usize::MAX, BondDirection::None);
                    let graph_snapshot = graph.clone();
                    let rings_snapshot = rings.clone();
                    assert!(rings.is_initialized());
                    assert_eq!(rings.find_type(), RingFindType::OtherOrUnknown);
                    assert_eq!(rings.atom_row_count(), 12);
                    assert_eq!(rings.bond_row_count(), 12);
                    let atoms_mask = [true; 12];
                    let atoms_snapshot = atoms_mask;
                    let mask_snapshot = bonds_mask;
                    let params = KekulizeParams {
                        canonical: false,
                        mark_atoms_bonds: false,
                        ..KekulizeParams::default()
                    };
                    let params_snapshot = params.clone();
                    let engine = kekulize_fragment(
                        &graph,
                        &atoms_mask,
                        &bonds_mask,
                        &params,
                        None,
                        Some(&rings),
                        None,
                    );
                    calls += 1;
                    assert_eq!(&graph, &graph_snapshot, "{label}: engine graph");
                    assert_eq!(rings, rings_snapshot, "{label}: engine rings");
                    assert_eq!(atoms_mask, atoms_snapshot, "{label}: engine atoms");
                    assert_eq!(bonds_mask, mask_snapshot, "{label}: engine mask");
                    assert_eq!(params, params_snapshot, "{label}: engine params");
                    assert_eq!(
                        engine.err(),
                        Some(KekulizeError::RingFinding(
                            crate::RingFindingError::ExpectedBondNotFound {
                                begin: AtomId::new(expected.0),
                                end: AtomId::new(expected.1),
                            }
                        )),
                        "{label}: engine"
                    );
                }
            }
            assert_eq!(calls, 8, "exact census");
        }

        // C4: 16 deque front-insertion calls. The wedge sits ONLY on the
        // second graph cycle, so B MUST precede A even for {A,B} input.
        #[test]
        fn kekulize_ring_convert_front_insertion() {
            let mut calls = 0usize;
            for rings_ab_first in [true, false] {
                for reverse_rows in [false, true] {
                    for direction in [BondDirection::BeginWedge, BondDirection::BeginDash] {
                        for wedge_begin in [6usize, 8usize] {
                            let graph = two_cycle_graph(wedge_begin, direction);
                            let graph_snapshot = graph.clone();
                            let a_fwd = (0..6).map(AtomId::new).collect::<Vec<_>>();
                            let b_fwd = (6..12).map(AtomId::new).collect::<Vec<_>>();
                            let a_bonds = (0..6).map(BondId::new).collect::<Vec<_>>();
                            let b_bonds = (6..12).map(BondId::new).collect::<Vec<_>>();
                            let mut a_rev = a_fwd.clone();
                            a_rev.reverse();
                            let mut b_rev = b_fwd.clone();
                            b_rev.reverse();
                            let mut a_bonds_rev = a_bonds.clone();
                            a_bonds_rev.reverse();
                            let mut b_bonds_rev = b_bonds.clone();
                            b_bonds_rev.reverse();
                            // F/R orientation applies to BOTH cycles per
                            // the frozen product axes.
                            let (row_a, row_b, sb_a, sb_b) = if reverse_rows {
                                (
                                    a_rev.clone(),
                                    b_rev.clone(),
                                    a_bonds_rev.clone(),
                                    b_bonds_rev.clone(),
                                )
                            } else {
                                (
                                    a_fwd.clone(),
                                    b_fwd.clone(),
                                    a_bonds.clone(),
                                    b_bonds.clone(),
                                )
                            };
                            let (first_a, first_b, s_first, s_second) = if rings_ab_first {
                                (row_a.clone(), row_b.clone(), sb_a.clone(), sb_b.clone())
                            } else {
                                (row_b.clone(), row_a.clone(), sb_b.clone(), sb_a.clone())
                            };
                            let rings = ring_info_from_selected_rows(
                                12,
                                12,
                                &[first_a, first_b],
                                &[s_first, s_second],
                            )
                            .unwrap();
                            let rings_snapshot = rings.clone();
                            let wedged = wedge_array(&graph, wedge_begin);
                            let wedged_snapshot = wedged.clone();
                            let label = format!(
                                "ab={rings_ab_first}/rev={reverse_rows}/{direction:?}/w{wedge_begin}"
                            );
                            // Actual graph wedge-bond prerequisites precede
                            // the call, on named arrays.
                            let atoms_selection = [true; 12];
                            let atoms_snapshot = atoms_selection;
                            let bonds_selection = [true; 12];
                            let bonds_snapshot = bonds_selection;
                            let dummy_input = [false; 12];
                            let dummy_snapshot = dummy_input;
                            assert_eq!(
                                graph.bonds[wedge_begin].direction(),
                                direction,
                                "{label}: wedge bond direction"
                            );
                            assert_eq!(
                                graph.bonds[wedge_begin].begin().index(),
                                wedge_begin,
                                "{label}: wedge bond begin"
                            );
                            assert!(wedged[wedge_begin], "{label}: graph wedge");
                            let (atom_rows, bond_rows) = collect_kekulize_candidate_rings(
                                &graph,
                                &rings,
                                &atoms_selection,
                                &bonds_selection,
                                &dummy_input,
                                &wedged,
                            )
                            .unwrap();
                            // Fresh six-input comparison immediately after
                            // THIS call, before any output; census after.
                            calls += 1;
                            assert_eq!(&graph, &graph_snapshot, "{label}: graph");
                            assert_eq!(rings, rings_snapshot, "{label}: rings");
                            assert_eq!(atoms_selection, atoms_snapshot, "{label}: atoms");
                            assert_eq!(bonds_selection, bonds_snapshot, "{label}: bonds");
                            assert_eq!(dummy_input, dummy_snapshot, "{label}: dummy");
                            assert_eq!(wedged, wedged_snapshot, "{label}: wedged");
                            // B (wedged) first, A second.
                            let (b_atoms, b_bonds_out): (Vec<usize>, Vec<usize>) =
                                match (reverse_rows, wedge_begin) {
                                    (false, 6) => {
                                        (vec![6, 7, 8, 9, 10, 11], vec![6, 7, 8, 9, 10, 11])
                                    }
                                    (false, _) => {
                                        (vec![8, 9, 10, 11, 6, 7], vec![8, 9, 10, 11, 6, 7])
                                    }
                                    (true, 6) => {
                                        (vec![6, 11, 10, 9, 8, 7], vec![11, 10, 9, 8, 7, 6])
                                    }
                                    (true, _) => {
                                        (vec![8, 7, 6, 11, 10, 9], vec![7, 6, 11, 10, 9, 8])
                                    }
                                };
                            let (a_atoms, a_bonds_out): (Vec<usize>, Vec<usize>) = if reverse_rows {
                                (vec![5, 4, 3, 2, 1, 0], vec![4, 3, 2, 1, 0, 5])
                            } else {
                                (vec![0, 1, 2, 3, 4, 5], vec![0, 1, 2, 3, 4, 5])
                            };
                            assert_eq!(
                                ids(&atom_rows),
                                vec![b_atoms.clone(), a_atoms.clone()],
                                "{label}: atoms"
                            );
                            assert_eq!(
                                bond_ids(&bond_rows),
                                vec![b_bonds_out.clone(), a_bonds_out.clone()],
                                "{label}: bonds"
                            );
                        }
                    }
                }
            }
            assert_eq!(calls, 16, "exact census");
        }
    }

    #[test]
    fn q05_core_sanitize_query_consumers_kekulize_recurses_through_bond_type_queries() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let query_atoms = graph
            .atoms
            .iter()
            .map(|carrier| {
                QueryAtom::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                )
            })
            .collect::<Vec<_>>();
        let explicit_bonds = vec![QueryBond::from_parts(
            graph.bonds[0].clone(),
            QueryNode::and(vec![
                QueryNode::or(vec![
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
                ]),
                QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
            ]),
        )];
        let explicit_state =
            QueryStateRef::try_for_topology(&query_atoms, &explicit_bonds, &graph).unwrap();
        assert!(explicit_state.bond_has_query(BondId::new(0)));
        assert_eq!(
            selected_bond_has_type_query(&graph.bonds[0], Some(explicit_state)),
            Ok(true)
        );

        let carrier_bonds = vec![QueryBond::from_carrier_parts(
            graph.bonds[0].clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        )];
        let carrier_state =
            QueryStateRef::try_for_topology(&query_atoms, &carrier_bonds, &graph).unwrap();
        assert!(!carrier_state.bond_has_query(BondId::new(0)));
        assert_eq!(
            selected_bond_has_type_query(&graph.bonds[0], Some(carrier_state)),
            Ok(false)
        );
    }

    #[test]
    fn canonical_rank_empty_and_disconnected_fragments_are_stable() {
        assert_eq!(
            rank_fragment_atoms(&TopologyBlock::default(), &[], &[]).unwrap(),
            Vec::<usize>::new()
        );

        let disconnected = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::O)),
                atom(2, AtomSpec::new(Element::N)),
            ],
            vec![],
        );
        let first = rank_fragment_atoms(&disconnected, &[true, true, true], &[]).unwrap();
        let second = rank_fragment_atoms(&disconnected, &[true, true, true], &[]).unwrap();
        assert_eq!(first, second);
        assert_eq!(first.len(), 3);
        let mut sorted = first;
        sorted.sort_unstable();
        assert_eq!(sorted, vec![0, 1, 2]);
    }

    #[test]
    fn fragment_rank_symbols_reverse_original_index_order_without_changing_text() {
        // Pinned MolFragmentToSmiles(CCO) with N/O/S and S/O/N produces
        // identical NOS text but output atom orders [0,1,2] and [2,1,0].
        // The direct pinned CanonicalRankAtomsInFragment oracle returns
        // [0,2,1] and [1,2,0]: traversal order is not the rank vector.
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
                atom(2, AtomSpec::new(Element::O)),
            ],
            vec![bond(0, 0, 1), bond(1, 1, 2)],
        );
        let options = CanonicalRankParams::kekulize_fragment_default();
        let first = ["N".to_owned(), "O".to_owned(), "S".to_owned()];
        let reverse = ["S".to_owned(), "O".to_owned(), "N".to_owned()];
        let plain =
            rank_fragment_atoms_with_params(&graph, &[true; 3], &[true; 2], None, None, &options)
                .unwrap();
        assert_eq!(
            plain,
            rank_fragment_atoms(&graph, &[true; 3], &[true; 2]).unwrap()
        );
        assert_eq!(
            rank_fragment_atoms_with_params(
                &graph,
                &[true; 3],
                &[true; 2],
                Some(&first),
                None,
                &options,
            ),
            Ok(vec![0, 2, 1])
        );
        assert_eq!(
            rank_fragment_atoms_with_params(
                &graph,
                &[true; 3],
                &[true; 2],
                Some(&reverse),
                None,
                &options,
            ),
            Ok(vec![1, 2, 0])
        );
    }

    #[test]
    fn fragment_rank_symbols_borrow_both_bond_holders_and_respect_independent_masks() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
                atom(2, AtomSpec::new(Element::O)),
            ],
            vec![bond(0, 0, 1), bond(1, 1, 2)],
        );
        let view = CanonRankReadView::from_topology(&graph).unwrap();
        let atoms = ["N".to_owned(), "O".to_owned(), "S".to_owned()];
        let first = ["-".to_owned(), "=".to_owned()];
        let reverse = ["=".to_owned(), "-".to_owned()];
        for symbols in [&first[..], &reverse[..]] {
            let initialized = init_fragment_canon_atoms(
                &view,
                &[true; 3],
                &[true; 2],
                true,
                Some(&atoms),
                Some(symbols),
            )
            .unwrap();
            assert_eq!(initialized[0].p_symbol, Some("N"));
            assert_eq!(initialized[1].p_symbol, Some("O"));
            assert_eq!(initialized[2].p_symbol, Some("S"));
            for (bond_index, symbol) in symbols.iter().enumerate() {
                let bond = &graph.bonds[bond_index];
                for endpoint in [bond.begin().index(), bond.end().index()] {
                    assert!(
                        initialized[endpoint]
                            .bonds
                            .iter()
                            .any(|holder| holder.bond_idx == bond_index
                                && holder.p_symbol == Some(symbol.as_str()))
                    );
                }
            }
            let masked_bond = init_fragment_canon_atoms(
                &view,
                &[true; 3],
                &[true, false],
                true,
                Some(&atoms),
                Some(symbols),
            )
            .unwrap();
            assert!(masked_bond[2].bonds.is_empty());
            let masked_atom = init_fragment_canon_atoms(
                &view,
                &[true, true, false],
                &[true; 2],
                true,
                Some(&atoms),
                Some(symbols),
            )
            .unwrap();
            assert_eq!(masked_atom[2].p_symbol, None);
            assert!(masked_atom[2].bonds.is_empty());
            assert_eq!(masked_atom[1].bonds.len(), 1);
        }
        // Pinned writer emits C-C=O and C=C-O for these two bond tables.
        let options = CanonicalRankParams::kekulize_fragment_default();
        for symbols in [&first[..], &reverse[..]] {
            assert_eq!(
                rank_fragment_atoms_with_params(
                    &graph,
                    &[true; 3],
                    &[true; 2],
                    None,
                    Some(symbols),
                    &options,
                )
                .unwrap()
                .len(),
                3
            );
        }
    }

    #[test]
    fn fragment_rank_symbols_validate_exact_lengths_before_empty_return() {
        let graph = topology(vec![atom(0, AtomSpec::new(Element::C))], vec![]);
        let options = CanonicalRankParams::kekulize_fragment_default();
        for size in [0, 2] {
            let symbols = vec!["C".to_owned(); size];
            assert_eq!(
                rank_fragment_atoms_with_params(
                    &graph,
                    &[true],
                    &[],
                    Some(&symbols),
                    None,
                    &options,
                ),
                Err(CanonicalRankError::AtomSymbolLength {
                    expected: 1,
                    actual: size,
                })
            );
        }
        let cco = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        for size in [0, 2] {
            let symbols = vec!["-".to_owned(); size];
            assert_eq!(
                rank_fragment_atoms_with_params(
                    &cco,
                    &[true; 2],
                    &[true],
                    None,
                    Some(&symbols),
                    &options,
                ),
                Err(CanonicalRankError::BondSymbolLength {
                    expected: 1,
                    actual: size,
                })
            );
        }
        let empty = TopologyBlock::default();
        let one = ["X".to_owned()];
        assert_eq!(
            rank_fragment_atoms_with_params(&empty, &[true], &[], Some(&one), None, &options),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 0,
                actual: 1,
            })
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&empty, &[], &[], Some(&one), None, &options),
            Err(CanonicalRankError::AtomSymbolLength {
                expected: 0,
                actual: 1,
            })
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&empty, &[], &[], None, Some(&one), &options),
            Err(CanonicalRankError::BondSymbolLength {
                expected: 0,
                actual: 1,
            })
        );
        assert_eq!(
            rank_fragment_atoms_with_params(&empty, &[], &[], None, None, &options),
            Ok(vec![])
        );
    }

    #[test]
    fn fragment_rank_flags_isolate_ties_isotopes_maps_and_chiral_presence() {
        let ethane = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let mut flags = CanonicalRankParams::default();
        flags.break_ties = false;
        let rank = |graph: &TopologyBlock, options: &CanonicalRankParams| {
            rank_fragment_atoms_with_params(
                graph,
                &vec![true; graph.atoms.len()],
                &vec![true; graph.bonds.len()],
                None,
                None,
                options,
            )
            .unwrap()
        };
        assert_eq!(rank(&ethane, &flags), vec![0, 0]);
        flags.break_ties = true;
        assert_eq!(rank(&ethane, &flags), vec![0, 1]);
        flags.break_ties = false;

        // Pinned CanonicalRankAtomsInFragment([13CH4].[12CH4]):
        // includeIsotopes=true -> [1,0], false -> [0,0].
        let isotopes = topology(
            vec![
                atom(0, AtomSpec::new(Element::C).with_isotope(13)),
                atom(1, AtomSpec::new(Element::C).with_isotope(12)),
            ],
            vec![],
        );
        assert_eq!(rank(&isotopes, &flags), vec![1, 0]);
        flags.include_isotopes = false;
        assert_eq!(rank(&isotopes, &flags), vec![0, 0]);
        flags.include_isotopes = true;

        let maps = topology(
            vec![
                atom(0, AtomSpec::new(Element::C).with_atom_map(9)),
                atom(1, AtomSpec::new(Element::C).with_atom_map(3)),
            ],
            vec![],
        );
        assert_eq!(rank(&maps, &flags), vec![1, 0]);
        flags.include_atom_maps = false;
        assert_eq!(rank(&maps, &flags), vec![0, 0]);
        flags.include_atom_maps = true;

        let chiral_presence = topology(
            vec![
                atom(
                    0,
                    AtomSpec::new(Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
                ),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![],
        );
        flags.include_chirality = false;
        flags.include_chiral_presence = false;
        assert_eq!(rank(&chiral_presence, &flags), vec![0, 0]);
        flags.include_chiral_presence = true;
        assert_eq!(rank(&chiral_presence, &flags), vec![1, 0]);
    }

    #[test]
    fn fragment_rank_flags_keep_fragment_chirality_ring_rule_and_masks() {
        let mut fragment = CanonicalRankParams::default();
        fragment.include_chirality = true;
        fragment.include_ring_stereo = false;
        fragment.chirality_rings_use_ring_stereo = false;
        let fragment_flags = CanonRankFlags::from_fragment_options(fragment);
        assert!(fragment_flags.use_chirality_rings);
        let whole_flags = CanonRankFlags::from_fragment_options(CanonicalRankParams {
            chirality_rings_use_ring_stereo: true,
            ..fragment
        });
        assert!(!whole_flags.use_chirality_rings);
        fragment.include_chirality = false;
        assert!(!CanonRankFlags::from_fragment_options(fragment).use_chirality_rings);

        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
                atom(2, AtomSpec::new(Element::C).with_isotope(13)),
            ],
            vec![bond(0, 0, 1), bond(1, 1, 2)],
        );
        let mut options = CanonicalRankParams::default();
        options.break_ties = false;
        let first = rank_fragment_atoms_with_params(
            &graph,
            &[true, true, false],
            &[true, true],
            None,
            None,
            &options,
        )
        .unwrap();
        assert_eq!(first.len(), 3);
        let mut changed = graph.clone();
        changed.atoms[2] = atom(2, AtomSpec::new(Element::C).with_isotope(14));
        let second = rank_fragment_atoms_with_params(
            &changed,
            &[true, true, false],
            &[true, false],
            None,
            None,
            &options,
        )
        .unwrap();
        assert_eq!(first, second);
    }

    #[test]
    fn fragment_rank_prepared_borrows_fast_rings_and_valence_without_mutation() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
                atom(2, AtomSpec::new(Element::O)),
            ],
            vec![bond(0, 0, 1), bond(1, 1, 2)],
        );
        let valence =
            crate::assign_valence_with_options_for_topology(&graph, ValenceModel::RdkitLike, false)
                .unwrap();
        let rings =
            fast_find_rings_from_parts(graph.atoms.len(), &graph.bonds, &graph.adjacency).unwrap();
        let original_graph = graph.clone();
        let original_valence = valence.clone();
        let original_rings = rings.clone();
        let view =
            CanonRankReadView::from_prepared_state(&graph, &valence, Some(&rings), None).unwrap();
        assert!(matches!(view.valence, Cow::Borrowed(_)));
        assert!(matches!(view.rings, Cow::Borrowed(_)));
        let options = CanonicalRankParams::kekulize_fragment_default();
        let prepared = rank_fragment_atoms_with_prepared_state(
            &graph,
            &valence,
            Some(&rings),
            &[true; 3],
            &[true; 2],
            None,
            None,
            &options,
        )
        .unwrap();
        assert_eq!(
            prepared,
            rank_fragment_atoms_with_params(&graph, &[true; 3], &[true; 2], None, None, &options,)
                .unwrap()
        );
        assert_eq!(graph, original_graph);
        assert_eq!(valence, original_valence);
        assert_eq!(rings, original_rings);
    }

    #[test]
    fn fragment_rank_prepared_recomputes_only_missing_or_unknown_rings() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let valence =
            crate::assign_valence_with_options_for_topology(&graph, ValenceModel::RdkitLike, false)
                .unwrap();
        let unknown = RingInfo::new(RingFindType::OtherOrUnknown, 2, 1);
        let original = unknown.clone();
        for supplied in [None, Some(&unknown)] {
            let view =
                CanonRankReadView::from_prepared_state(&graph, &valence, supplied, None).unwrap();
            assert!(matches!(view.valence, Cow::Borrowed(_)));
            assert!(matches!(view.rings, Cow::Owned(_)));
            assert!(view.rings.is_find_fast_or_better());
            assert_eq!(
                rank_fragment_atoms_with_prepared_state(
                    &graph,
                    &valence,
                    supplied,
                    &[true; 2],
                    &[true],
                    None,
                    None,
                    &CanonicalRankParams::kekulize_fragment_default(),
                ),
                rank_fragment_atoms(&graph, &[true; 2], &[true])
            );
        }
        assert_eq!(unknown, original);
    }

    #[test]
    fn fragment_rank_prepared_rejects_invalid_rows_and_accepts_no_implicit() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C).with_no_implicit(true)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let mut valence =
            crate::assign_valence_with_options_for_topology(&graph, ValenceModel::RdkitLike, false)
                .unwrap();
        valence.implicit_hydrogens[0] = -1;
        assert!(CanonRankReadView::from_prepared_state(&graph, &valence, None, None).is_ok());
        let mut bad = valence.clone();
        bad.explicit_valence.pop();
        assert!(matches!(
            CanonRankReadView::from_prepared_state(&graph, &bad, None, None),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 1,
                implicit_len: 2,
            })
        ));
        bad = valence.clone();
        bad.implicit_hydrogens[1] = -1;
        assert!(matches!(
            CanonRankReadView::from_prepared_state(&graph, &bad, None, None),
            Err(CanonicalRankError::PreparedValenceInvalid { atom_index: 1 })
        ));
        let wrong_rings = RingInfo::new(RingFindType::Fast, 1, 0);
        assert!(matches!(
            CanonRankReadView::from_prepared_state(&graph, &valence, Some(&wrong_rings), None),
            Err(CanonicalRankError::PreparedRingLength {
                expected_atoms: 2,
                actual_atoms: 1,
                expected_bonds: 1,
                actual_bonds: 0,
            })
        ));
    }

    #[test]
    fn canonical_rank_masks_control_fragment_initialization() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
                atom(2, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1), bond(1, 1, 2)],
        );
        let view = CanonRankReadView::from_topology(&graph).unwrap();
        let atoms =
            init_fragment_canon_atoms(&view, &[true, true, false], &[true, true], true, None, None)
                .unwrap();
        assert_eq!(
            atoms.iter().map(|atom| atom.is_in_play).collect::<Vec<_>>(),
            vec![true, true, false]
        );
        assert_eq!(
            atoms.iter().map(|atom| atom.degree).collect::<Vec<_>>(),
            vec![1, 1, 0]
        );
        assert_eq!(atoms[0].nbr_ids, vec![1]);
        assert_eq!(atoms[1].nbr_ids, vec![0]);
        assert!(atoms[2].nbr_ids.is_empty());
        assert_eq!(atoms[0].bonds.len(), 1);
        assert_eq!(atoms[1].bonds.len(), 1);
        assert!(atoms[2].bonds.is_empty());

        let bond_masked =
            init_fragment_canon_atoms(&view, &[true, true, true], &[true, false], true, None, None)
                .unwrap();
        assert_eq!(
            bond_masked
                .iter()
                .map(|atom| atom.degree)
                .collect::<Vec<_>>(),
            vec![1, 1, 0]
        );
        assert_eq!(bond_masked[1].total_num_hs, 3);
        assert_eq!(bond_masked[2].total_num_hs, 4);
    }

    #[test]
    fn canonical_rank_preserves_symmetry_classes_until_stable_tie_breaking() {
        let ethane = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let view = CanonRankReadView::from_topology(&ethane).unwrap();
        let mut atoms =
            init_fragment_canon_atoms(&view, &[true, true], &[true], true, None, None).unwrap();
        let mut options = CanonicalRankParams::kekulize_fragment_default();
        options.break_ties = false;
        assert_eq!(
            rank_initialized_atoms(&view, &mut atoms, options).unwrap(),
            vec![0, 0]
        );
        assert_eq!(
            rank_fragment_atoms(&ethane, &[true, true], &[true]).unwrap(),
            vec![0, 1]
        );
    }

    #[test]
    fn canonical_rank_excludes_inactive_isotope_map_and_chirality_dimensions() {
        let graph = topology(
            vec![
                atom(
                    0,
                    AtomSpec::new(Element::C)
                        .with_isotope(13)
                        .with_atom_map(9)
                        .with_chiral_tag(ChiralTag::TetrahedralCw),
                ),
                atom(
                    1,
                    AtomSpec::new(Element::C)
                        .with_isotope(12)
                        .with_atom_map(3)
                        .with_chiral_tag(ChiralTag::TetrahedralCcw),
                ),
            ],
            vec![],
        );
        let view = CanonRankReadView::from_topology(&graph).unwrap();
        let mut inactive =
            init_fragment_canon_atoms(&view, &[false, false], &[], true, None, None).unwrap();
        let mut no_tie_break = CanonicalRankParams::kekulize_fragment_default();
        no_tie_break.break_ties = false;
        assert_eq!(
            rank_initialized_atoms(&view, &mut inactive, no_tie_break).unwrap(),
            vec![0, 0]
        );

        let active = rank_fragment_atoms(&graph, &[true, true], &[]).unwrap();
        assert_ne!(active[0], active[1]);
    }

    #[test]
    fn canonical_rank_rejects_malformed_masks() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        assert!(matches!(
            rank_fragment_atoms(&graph, &[true], &[true]),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        ));
        assert!(matches!(
            rank_fragment_atoms(&graph, &[true, true], &[]),
            Err(CanonicalRankError::BondMaskLength {
                expected: 1,
                actual: 0
            })
        ));
    }

    #[test]
    fn candidate_selection_validates_masks_and_filters_queries_and_ring_boundaries() {
        let aromatic_ring = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        assert!(matches!(
            prepare_kekulize_selection(&aromatic_ring, &[true, true], &[true; 3], None),
            Err(KekulizeError::AtomSelectionLength {
                expected: 3,
                actual: 2
            })
        ));
        assert!(matches!(
            prepare_kekulize_selection(&aromatic_ring, &[true; 3], &[true; 2], None),
            Err(KekulizeError::BondSelectionLength {
                expected: 3,
                actual: 2
            })
        ));

        let empty =
            prepare_kekulize_selection(&aromatic_ring, &[false; 3], &[true; 3], None).unwrap();
        assert!(!empty.found_aromatic);
        assert!(empty.candidate_atom_rings.is_empty());
        assert_eq!(empty.original_total_valences, vec![0; 3]);

        let crossing =
            prepare_kekulize_selection(&aromatic_ring, &[true, true, false], &[true; 3], None)
                .unwrap();
        assert!(crossing.found_aromatic);
        assert!(crossing.candidate_atom_rings.is_empty());
        assert!(crossing.candidate_bond_rings.is_empty());

        let query_ring = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let query_atoms = query_ring
            .atoms
            .iter()
            .map(|carrier| {
                QueryAtom::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                )
            })
            .collect::<Vec<_>>();
        let mut query_bonds = query_ring
            .bonds
            .iter()
            .map(|carrier| {
                QueryBond::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
                )
            })
            .collect::<Vec<_>>();
        query_bonds[0] = QueryBond::from_parts(
            query_ring.bonds[0].clone(),
            QueryNode::or(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            ]),
        );
        let query_state =
            QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &query_ring).unwrap();
        let query =
            prepare_kekulize_selection(&query_ring, &[true; 3], &[true; 3], Some(query_state))
                .unwrap();
        assert_eq!(query.bonds_in_play, vec![false, true, true]);
        assert!(query.candidate_atom_rings.is_empty());

        let complex = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond(0, 0, 1)],
        );
        let complex_atoms = complex
            .atoms
            .iter()
            .map(|carrier| {
                QueryAtom::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                )
            })
            .collect::<Vec<_>>();
        let complex_bonds = vec![QueryBond::from_parts(
            complex.bonds[0].clone(),
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Any),
                QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
            ]),
        )];
        let complex_state =
            QueryStateRef::try_for_topology(&complex_atoms, &complex_bonds, &complex).unwrap();
        let complex =
            prepare_kekulize_selection(&complex, &[true; 2], &[true], Some(complex_state)).unwrap();
        assert_eq!(complex.bonds_in_play, vec![true]);
    }

    #[test]
    fn candidate_aromatic_bonds_are_accounted_before_staged_single_conversion() {
        let aromatic_ring = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let rings = fast_find_rings_from_parts(
            aromatic_ring.atoms.len(),
            &aromatic_ring.bonds,
            &aromatic_ring.adjacency,
        )
        .unwrap();
        let state = mark_double_bond_candidates(
            &aromatic_ring,
            &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            &rings,
            &valence(vec![3; 3], vec![1; 3]),
        )
        .unwrap();
        assert_eq!(state.double_bond_candidates, vec![true; 3]);
        assert!(state.questions.is_empty());
        assert!(state.done.is_empty());
        assert!(
            state
                .topology
                .bonds
                .iter()
                .all(|bond| bond.order() == BondOrder::Single && bond.is_aromatic())
        );
        assert!(
            aromatic_ring
                .bonds
                .iter()
                .all(|bond| bond.order() == BondOrder::Aromatic)
        );
    }

    #[test]
    fn candidate_dummy_question_respects_candidate_and_non_candidate_rings() {
        let mixed_ring = topology(
            vec![
                atom(0, AtomSpec::new(Element::DUMMY)),
                atom(1, AtomSpec::new(Element::C).with_aromatic(true)),
                atom(2, AtomSpec::new(Element::C).with_aromatic(true)),
            ],
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let mixed_rings = fast_find_rings_from_parts(
            mixed_ring.atoms.len(),
            &mixed_ring.bonds,
            &mixed_ring.adjacency,
        )
        .unwrap();
        let mixed = mark_double_bond_candidates(
            &mixed_ring,
            &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            &mixed_rings,
            &valence(vec![0, 3, 3], vec![0, 1, 1]),
        )
        .unwrap();
        assert!(mixed.double_bond_candidates[0]);
        assert_eq!(mixed.questions, vec![AtomId::new(0)]);

        let all_dummy_ring = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::DUMMY)))
                .collect(),
            vec![bond(0, 0, 1), bond(1, 1, 2), bond(2, 2, 0)],
        );
        let all_dummy_rings = fast_find_rings_from_parts(
            all_dummy_ring.atoms.len(),
            &all_dummy_ring.bonds,
            &all_dummy_ring.adjacency,
        )
        .unwrap();
        let all_dummy = mark_double_bond_candidates(
            &all_dummy_ring,
            &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            &all_dummy_rings,
            &valence(vec![0; 3], vec![0; 3]),
        )
        .unwrap();
        assert_eq!(all_dummy.double_bond_candidates, vec![false; 3]);
        assert!(all_dummy.questions.is_empty());
    }

    #[test]
    fn candidate_charge_radical_hydrogen_no_implicit_and_degree_branches_are_exact() {
        let graph = topology(
            vec![
                atom(
                    0,
                    AtomSpec::new(Element::B)
                        .with_aromatic(true)
                        .with_formal_charge(1)
                        .with_explicit_hydrogens(1),
                ),
                atom(
                    1,
                    AtomSpec::new(Element::C)
                        .with_aromatic(true)
                        .with_formal_charge(1)
                        .with_explicit_hydrogens(2),
                ),
                atom(
                    2,
                    AtomSpec::new(Element::C)
                        .with_aromatic(true)
                        .with_explicit_hydrogens(2)
                        .with_radical_electrons(1),
                ),
                atom(
                    3,
                    AtomSpec::new(Element::C)
                        .with_aromatic(true)
                        .with_explicit_hydrogens(2)
                        .with_no_implicit(true),
                ),
                atom(4, AtomSpec::new(Element::C).with_aromatic(true)),
            ],
            vec![],
        );
        let rings = fast_find_rings_from_parts(5, &graph.bonds, &graph.adjacency).unwrap();
        let state = mark_double_bond_candidates(
            &graph,
            &[
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
                AtomId::new(4),
            ],
            &rings,
            &valence(vec![1, 2, 2, 2, 0], vec![0, 0, 0, 0, 4]),
        )
        .unwrap();
        assert_eq!(
            state.double_bond_candidates,
            vec![true, true, true, true, false]
        );
        assert!(is_early_atom_for_kekulize(Element::B.atomic_number()));
        assert!(!is_early_atom_for_kekulize(Element::C.atomic_number()));
    }

    #[test]
    fn candidate_n_p_as_n_oxide_source_cases_take_a_double_bond() {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        for (group, element) in [Element::N, Element::P, Element::AS]
            .into_iter()
            .enumerate()
        {
            let center = group * 4;
            atoms.push(atom(center, AtomSpec::new(element).with_aromatic(true)));
            for offset in 1..=3 {
                atoms.push(atom(center + offset, AtomSpec::new(Element::C)));
                bonds.push(bond_with_spec(
                    bonds.len(),
                    BondSpec::new(
                        AtomId::new(center),
                        AtomId::new(center + offset),
                        if offset == 1 {
                            BondOrder::Double
                        } else {
                            BondOrder::Single
                        },
                    ),
                ));
            }
        }
        let graph = topology(atoms, bonds);
        let rings =
            fast_find_rings_from_parts(graph.atoms.len(), &graph.bonds, &graph.adjacency).unwrap();
        let mut explicit = vec![0; graph.atoms.len()];
        for center in [0, 4, 8] {
            explicit[center] = 5;
        }
        let state = mark_double_bond_candidates(
            &graph,
            &[AtomId::new(0), AtomId::new(4), AtomId::new(8)],
            &rings,
            &valence(explicit, vec![0; graph.atoms.len()]),
        )
        .unwrap();
        assert!(state.double_bond_candidates[0]);
        assert!(state.double_bond_candidates[4]);
        assert!(state.double_bond_candidates[8]);
    }

    #[test]
    fn candidate_zero_contribution_is_ignored_in_total_degree() {
        let graph = topology(
            vec![
                atom(
                    0,
                    AtomSpec::new(Element::C)
                        .with_aromatic(true)
                        .with_explicit_hydrogens(3),
                ),
                atom(1, AtomSpec::new(Element::H)),
            ],
            vec![bond_with_spec(
                0,
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Hydrogen),
            )],
        );
        let rings = fast_find_rings_from_parts(2, &graph.bonds, &graph.adjacency).unwrap();
        let state = mark_double_bond_candidates(
            &graph,
            &[AtomId::new(0)],
            &rings,
            &valence(vec![3, 0], vec![0, 0]),
        )
        .unwrap();
        assert!(state.double_bond_candidates[0]);
    }

    #[test]
    fn candidate_rejects_invalid_duplicate_and_overflow_state_without_editing_input() {
        let graph = topology(
            vec![atom(
                0,
                AtomSpec::new(Element::C)
                    .with_aromatic(true)
                    .with_explicit_hydrogens(1),
            )],
            vec![],
        );
        let rings = fast_find_rings_from_parts(1, &graph.bonds, &graph.adjacency).unwrap();
        assert!(matches!(
            mark_double_bond_candidates(
                &graph,
                &[AtomId::new(1)],
                &rings,
                &valence(vec![0], vec![0])
            ),
            Err(KekulizeError::CandidateAtomOutOfRange {
                atom,
                atom_count: 1
            }) if atom == AtomId::new(1)
        ));
        assert!(matches!(
            mark_double_bond_candidates(
                &graph,
                &[AtomId::new(0), AtomId::new(0)],
                &rings,
                &valence(vec![0], vec![0])
            ),
            Err(KekulizeError::DuplicateCandidateAtom { atom })
                if atom == AtomId::new(0)
        ));
        assert!(matches!(
            mark_double_bond_candidates(
                &graph,
                &[AtomId::new(0)],
                &rings,
                &valence(vec![0], vec![i32::MAX])
            ),
            Err(KekulizeError::IntegerOverflow {
                atom,
                field: "total hydrogen count"
            }) if atom == AtomId::new(0)
        ));
        assert!(graph.bonds.is_empty());
        assert!(graph.atoms[0].is_aromatic());
    }

    #[test]
    fn matching_uses_rank_then_index_for_equal_rank_options() {
        let graph = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::DUMMY)))
                .collect(),
            vec![kekulize_ready_bond(0, 0, 1), kekulize_ready_bond(1, 0, 2)],
        );
        let result = kekulize_matching_worker(
            graph,
            &[AtomId::new(0), AtomId::new(1), AtomId::new(2)],
            &[true; 3],
            &[false; 2],
            &[],
            &[0; 3],
            0,
        )
        .unwrap();
        assert!(result.succeeded);
        assert_eq!(result.topology.bonds[0].order(), BondOrder::Double);
        assert_eq!(result.topology.bonds[1].order(), BondOrder::Single);
        assert_eq!(
            result.done,
            vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]
        );
    }

    #[test]
    fn matching_prioritizes_wedge_end_and_non_wedged_option_and_clears_selected_direction() {
        let wedge = bond_with_spec(
            0,
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_aromatic(true)
                .with_direction(BondDirection::BeginWedge),
        );
        let graph = topology(
            (0..4)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                wedge,
                kekulize_ready_bond(1, 1, 2),
                kekulize_ready_bond(2, 0, 3),
            ],
        );
        let result = kekulize_matching_worker(
            graph,
            &[
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
            ],
            &[true; 4],
            &[false; 3],
            &[],
            &[0; 4],
            0,
        )
        .unwrap();
        assert!(result.succeeded);
        assert_eq!(result.done[0], AtomId::new(1));
        assert_eq!(result.topology.bonds[0].order(), BondOrder::Single);
        assert_eq!(
            result.topology.bonds[0].direction(),
            BondDirection::BeginWedge
        );
        assert_eq!(result.topology.bonds[1].order(), BondOrder::Double);
        assert_eq!(result.topology.bonds[2].order(), BondOrder::Double);

        let selected_wedge = topology(
            vec![
                atom(0, AtomSpec::new(Element::C).with_aromatic(true)),
                atom(1, AtomSpec::new(Element::C).with_aromatic(true)),
            ],
            vec![bond_with_spec(
                0,
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_aromatic(true)
                    .with_direction(BondDirection::BeginDash),
            )],
        );
        let selected = kekulize_matching_worker(
            selected_wedge,
            &[AtomId::new(0), AtomId::new(1)],
            &[true; 2],
            &[false],
            &[],
            &[0; 2],
            0,
        )
        .unwrap();
        assert_eq!(selected.topology.bonds[0].order(), BondOrder::Double);
        assert_eq!(selected.topology.bonds[0].direction(), BondDirection::None);
    }

    #[test]
    fn matching_restores_state_and_obeys_below_at_above_backtrack_limit() {
        let graph = topology(
            (0..4)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                kekulize_ready_bond(0, 0, 1),
                kekulize_ready_bond(1, 0, 2),
                kekulize_ready_bond(2, 1, 3),
            ],
        );
        let run = |max_backtracks| {
            kekulize_matching_worker(
                graph.clone(),
                &[
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                ],
                &[true; 4],
                &[false; 3],
                &[],
                &[0, 1, 2, 3],
                max_backtracks,
            )
            .unwrap()
        };

        let below = run(0);
        assert!(!below.succeeded);
        assert_eq!(below.backtracks, 0);
        assert_eq!(below.problem_atoms, vec![AtomId::new(2), AtomId::new(3)]);
        assert!(
            below
                .topology
                .bonds
                .iter()
                .all(|bond| bond.order() == BondOrder::Single)
        );

        let at = run(1);
        assert!(at.succeeded);
        assert_eq!(at.backtracks, 1);
        assert_eq!(at.topology.bonds[0].order(), BondOrder::Single);
        assert_eq!(at.topology.bonds[1].order(), BondOrder::Double);
        assert_eq!(at.topology.bonds[2].order(), BondOrder::Double);
        assert_eq!(at.double_bond_candidates, vec![false; 4]);
        assert_eq!(at.double_bonds_added, vec![false, true, true]);

        let above = run(2);
        assert!(above.succeeded);
        assert_eq!(above.backtracks, 1);
        assert_eq!(above.done, at.done);
        assert_eq!(
            above
                .topology
                .bonds
                .iter()
                .map(Bond::order)
                .collect::<Vec<_>>(),
            at.topology
                .bonds
                .iter()
                .map(Bond::order)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn matching_rejects_malformed_vectors_and_duplicate_rows() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::C).with_aromatic(true)),
                atom(1, AtomSpec::new(Element::C).with_aromatic(true)),
            ],
            vec![kekulize_ready_bond(0, 0, 1)],
        );
        assert!(matches!(
            kekulize_matching_worker(
                graph.clone(),
                &[AtomId::new(0), AtomId::new(1)],
                &[true],
                &[false],
                &[],
                &[0, 1],
                0
            ),
            Err(KekulizeError::MatchingStateLength {
                field: "double-bond candidates",
                expected: 2,
                actual: 1
            })
        ));
        assert!(matches!(
            kekulize_matching_worker(
                graph,
                &[AtomId::new(0), AtomId::new(1)],
                &[true; 2],
                &[false],
                &[AtomId::new(0), AtomId::new(0)],
                &[0, 1],
                0
            ),
            Err(KekulizeError::DuplicateDoneAtom { atom }) if atom == AtomId::new(0)
        ));
    }

    #[test]
    fn fused_dummy_subsets_follow_source_bit_order_without_counter_limit() {
        let questions = vec![AtomId::new(2), AtomId::new(4), AtomId::new(6)];
        let mut enumerator = QuestionEnumerator::new(questions).unwrap();
        assert_eq!(enumerator.next(), vec![AtomId::new(2)]);
        assert_eq!(enumerator.next(), vec![AtomId::new(4)]);
        assert_eq!(enumerator.next(), vec![AtomId::new(2), AtomId::new(4)]);
        assert_eq!(enumerator.next(), vec![AtomId::new(6)]);
        assert_eq!(enumerator.next(), vec![AtomId::new(2), AtomId::new(6)]);
        assert_eq!(enumerator.next(), vec![AtomId::new(4), AtomId::new(6)]);
        assert_eq!(
            enumerator.next(),
            vec![AtomId::new(2), AtomId::new(4), AtomId::new(6)]
        );
        assert!(enumerator.next().is_empty());

        let mut wide = QuestionEnumerator::new(vec![AtomId::new(0); 32]).unwrap();
        assert_eq!(wide.next(), vec![AtomId::new(0)]);
    }

    #[test]
    fn fused_neighbor_grouping_covers_isolated_fused_and_disconnected_ring_systems() {
        let bond_rings = vec![
            vec![BondId::new(0), BondId::new(1), BondId::new(2)],
            vec![BondId::new(2), BondId::new(3), BondId::new(4)],
            vec![BondId::new(5), BondId::new(6), BondId::new(7)],
            vec![BondId::new(7), BondId::new(8), BondId::new(9)],
            vec![BondId::new(10), BondId::new(11), BondId::new(12)],
        ];
        let neighbors = make_ring_neighbor_map(&bond_rings);
        assert_eq!(neighbors, vec![vec![1], vec![0], vec![3], vec![2], vec![]]);

        let mut done = vec![false; bond_rings.len()];
        assert_eq!(pick_fused_rings(0, &neighbors, &mut done), vec![0, 1]);
        assert_eq!(pick_fused_rings(2, &neighbors, &mut done), vec![2, 3]);
        assert_eq!(pick_fused_rings(4, &neighbors, &mut done), vec![4]);
        assert_eq!(done, vec![true; 5]);
    }

    #[test]
    fn fused_all_dummy_and_crossing_rings_are_excluded_before_dispatch() {
        let all_dummy = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::DUMMY).with_aromatic(true)))
                .collect(),
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let dummy_selection =
            prepare_kekulize_selection(&all_dummy, &[true; 3], &[true; 3], None).unwrap();
        assert!(dummy_selection.found_aromatic);
        assert!(dummy_selection.candidate_atom_rings.is_empty());
        assert!(dummy_selection.candidate_bond_rings.is_empty());
        let dummy_result = kekulize_fused_components(
            all_dummy.clone(),
            &dummy_selection.candidate_atom_rings,
            &dummy_selection.candidate_bond_rings,
            &dummy_selection.rings,
            &dummy_selection.valence,
            &[0, 1, 2],
            0,
        )
        .unwrap();
        assert!(dummy_result.succeeded);
        assert_eq!(dummy_result.topology, all_dummy);

        let crossing = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                aromatic_bond(0, 0, 1),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let crossing_selection =
            prepare_kekulize_selection(&crossing, &[true, true, false], &[true; 3], None).unwrap();
        assert!(crossing_selection.found_aromatic);
        assert!(crossing_selection.candidate_atom_rings.is_empty());
        assert!(crossing_selection.candidate_bond_rings.is_empty());
    }

    #[test]
    fn fused_component_dispatch_handles_shared_and_disconnected_ring_systems() {
        let graph = topology(
            (0..16)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            [
                (0, 1),
                (1, 2),
                (2, 3),
                (3, 4),
                (4, 5),
                (5, 0),
                (5, 6),
                (6, 7),
                (7, 8),
                (8, 9),
                (9, 4),
                (10, 11),
                (11, 12),
                (12, 13),
                (13, 14),
                (14, 15),
                (15, 10),
            ]
            .into_iter()
            .enumerate()
            .map(|(id, (begin, end))| aromatic_bond(id, begin, end))
            .collect(),
        );
        let selection = prepare_kekulize_selection(&graph, &[true; 16], &[true; 17], None).unwrap();
        assert_eq!(selection.candidate_atom_rings.len(), 3);
        let neighbors = make_ring_neighbor_map(&selection.candidate_bond_rings);
        assert_eq!(neighbors.iter().map(Vec::len).sum::<usize>(), 2);

        let result = kekulize_fused_components(
            graph,
            &selection.candidate_atom_rings,
            &selection.candidate_bond_rings,
            &selection.rings,
            &selection.valence,
            &(0..16).collect::<Vec<_>>(),
            100,
        )
        .unwrap();
        assert!(
            result.succeeded,
            "problem atoms: {:?}",
            result.problem_atoms
        );
        assert_eq!(
            result
                .topology
                .bonds
                .iter()
                .filter(|bond| bond.order() == BondOrder::Double)
                .count(),
            8
        );
    }

    #[test]
    fn fused_dummy_retry_resets_aromatic_in_play_bonds_before_each_attempt() {
        let graph = topology(
            vec![
                atom(0, AtomSpec::new(Element::DUMMY).with_aromatic(true)),
                atom(1, AtomSpec::new(Element::DUMMY).with_aromatic(true)),
                atom(2, AtomSpec::new(Element::C).with_aromatic(true)),
                atom(3, AtomSpec::new(Element::DUMMY).with_aromatic(true)),
            ],
            vec![
                kekulize_ready_bond(0, 0, 2),
                bond_with_spec(
                    1,
                    BondSpec::new(AtomId::new(1), AtomId::new(3), BondOrder::Double)
                        .with_aromatic(true),
                ),
            ],
        );
        let result = permute_dummies_and_kekulize(
            graph,
            &[
                AtomId::new(0),
                AtomId::new(1),
                AtomId::new(2),
                AtomId::new(3),
            ],
            &[true, true, true, false],
            &[AtomId::new(0), AtomId::new(1)],
            &[0, 1, 2, 3],
            0,
        )
        .unwrap();
        assert!(result.succeeded);
        assert_eq!(result.topology.bonds[0].order(), BondOrder::Double);
        assert_eq!(result.topology.bonds[1].order(), BondOrder::Single);
    }
}

#[cfg(test)]
mod q01_b1_scalar_rank_dependency_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, PropertyValueKind};
    use cosmolkit_types::Element;

    fn source_atom(index: usize, value: Option<PropertyValue>) -> Atom {
        let mut atom = Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C));
        if let Some(value) = value {
            atom.set_prop("_CanonicalRankingNumber", value).unwrap();
        }
        atom
    }
    fn flags(enabled: bool) -> CanonRankFlags {
        let mut params = CanonicalRankParams::default();
        params.use_non_stereo_ranks = enabled;
        CanonRankFlags::from_fragment_options(params)
    }
    fn cast_error(index: usize, kind: PropertyValueKind) -> CanonicalRankError {
        CanonicalRankError::InvalidPropertyKind {
            atom_index: index,
            property: "_CanonicalRankingNumber",
            kind,
        }
    }

    #[test]
    fn q01_b1_rank_class_and_flag_guards_precede_integer_getters() {
        let source = [
            source_atom(0, Some(PropertyValue::String(("bad".to_owned()).into()))),
            source_atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Ok(Ordering::Greater)
        );
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(cast_error(0, PropertyValueKind::String))
        );
    }

    #[test]
    fn q01_b1_rank_getter_order_uses_actual_comparator_side() {
        let source = [
            source_atom(0, Some(PropertyValue::String(("bad".to_owned()).into()))),
            source_atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        atoms[0].index = 0;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(cast_error(0, PropertyValueKind::String))
        );
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(cast_error(1, PropertyValueKind::Double))
        );
        let source = [
            source_atom(0, Some(PropertyValue::Int(i32::MAX))),
            source_atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        atoms[0].index = 0;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(cast_error(1, PropertyValueKind::Bool))
        );
    }

    #[test]
    fn q01_b1_rank_modes_and_in_play_masks_preserve_source_no_read() {
        let source = [
            source_atom(0, Some(PropertyValue::String(("bad".to_owned()).into()))),
            source_atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        atoms[0].index = 0;
        atoms[1].index = 0;
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Err(cast_error(0, PropertyValueKind::String))
        );
    }

    #[test]
    fn q01_b1_rank_hanoi_guard_and_recursive_error_order_are_source_exact() {
        let source = [
            source_atom(0, Some(PropertyValue::String(("bad".to_owned()).into()))),
            source_atom(1, Some(PropertyValue::Int(1))),
            source_atom(2, Some(PropertyValue::Bool(true))),
            source_atom(3, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for atom in &mut atoms {
            atom.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts[0], 1);
        counts.fill(0);
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 2],
                &mut [0; 2],
                &mut counts,
                &[false; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0, 0, 0]);
        counts.fill(0);
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 2],
                &mut [0; 2],
                &mut counts,
                &[true, false, false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(cast_error(0, PropertyValueKind::String))
        );
        counts.fill(0);
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(cast_error(2, PropertyValueKind::Bool))
        );
        counts.fill(0);
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1, 2, 3],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(cast_error(0, PropertyValueKind::String))
        );
    }
}

#[cfg(test)]
mod uint_source_comparator_proposed_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, PropertyValueKind};
    use cosmolkit_types::Element;
    fn atom(index: usize, value: PropertyValue) -> Atom {
        Atom::from_spec(
            AtomId::new(index),
            AtomSpec::new(Element::C)
                .with_computed_prop("_CanonicalRankingNumber", value)
                .unwrap(),
        )
    }
    fn flags(on: bool) -> CanonRankFlags {
        let mut p = CanonicalRankParams::default();
        p.use_non_stereo_ranks = on;
        CanonRankFlags::from_fragment_options(p)
    }
    fn overflow(index: usize, value: u32) -> CanonicalRankError {
        CanonicalRankError::UnsignedRankOverflow {
            atom_index: index,
            property: "_CanonicalRankingNumber",
            value,
        }
    }
    #[test]
    fn proposed_uint_source_int_limits_and_decimal_projection() {
        // Independent Boost GT_HiT / numeric_limits<int>::max conditions.
        for (value, expected) in [
            (0_u32, Some(0)),
            (1, Some(1)),
            (2147483646, Some(2147483646)),
            (2147483647, Some(2147483647)),
            (2147483648, None),
            (4294967295, None),
        ] {
            let source = atom(7, PropertyValue::UInt(value));
            let before = source.clone();
            assert_eq!(
                canonical_rank_property_to_int(&source, 7),
                expected.ok_or_else(|| overflow(7, value))
            );
            assert_eq!(source, before);
        }
    }
    #[test]
    fn proposed_uint_competing_getter_errors_follow_actual_left_side() {
        for value in [2147483648_u32, 4294967295] {
            let source = [
                atom(0, PropertyValue::UInt(value)),
                atom(1, PropertyValue::String("bad".into())),
            ];
            let mut atoms = source
                .iter()
                .map(empty_canon_atom_from_source_atom)
                .collect::<Vec<_>>();
            for atom in &mut atoms {
                atom.index = 0;
            }
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
                Err(overflow(0, value))
            );
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
                Err(CanonicalRankError::InvalidPropertyKind {
                    atom_index: 1,
                    property: "_CanonicalRankingNumber",
                    kind: PropertyValueKind::String
                })
            );
            let source = [
                atom(0, PropertyValue::Int(0)),
                atom(1, PropertyValue::UInt(value)),
            ];
            let mut atoms = source
                .iter()
                .map(empty_canon_atom_from_source_atom)
                .collect::<Vec<_>>();
            for atom in &mut atoms {
                atom.index = 0;
            }
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
                Err(overflow(1, value))
            );
        }
    }
    #[test]
    fn proposed_uint_class_flag_modes_masks_hanoi_and_recursion() {
        for value in [0_u32, 1, 2147483646, 2147483647, 2147483648, 4294967295] {
            let source = [
                atom(0, PropertyValue::UInt(value)),
                atom(1, PropertyValue::UInt(value)),
            ];
            let mut atoms = source
                .iter()
                .map(empty_canon_atom_from_source_atom)
                .collect::<Vec<_>>();
            atoms[1].index = 1;
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
                Ok(Ordering::Less)
            );
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
                Ok(Ordering::Greater)
            );
            atoms[1].index = 0;
            assert_eq!(
                compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
                Ok(Ordering::Equal)
            );
            for mode in [
                CanonCompareMode::SpecialChirality,
                CanonCompareMode::SpecialSymmetry,
            ] {
                assert_eq!(
                    compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                    Ok(Ordering::Equal)
                );
            }
            atoms[0].is_in_play = false;
            atoms[1].is_in_play = false;
            assert_eq!(
                compare_canon_atoms_for_kekulize(
                    &mut atoms,
                    0,
                    1,
                    CanonCompareMode::Atom,
                    flags(true)
                ),
                Ok(Ordering::Equal)
            );
            let mut counts = [0; 2];
            assert_eq!(
                hanoi_order_for_kekulize(
                    &mut [0],
                    &mut [0],
                    &mut counts,
                    &[true; 2],
                    &mut atoms,
                    CanonCompareMode::Atom,
                    flags(true)
                ),
                Ok(false)
            );
            counts.fill(0);
            assert_eq!(
                hanoi_order_for_kekulize(
                    &mut [0, 1],
                    &mut [0; 2],
                    &mut counts,
                    &[false; 2],
                    &mut atoms,
                    CanonCompareMode::Atom,
                    flags(true)
                ),
                Ok(false)
            );
            atoms[1].is_in_play = true;
            let reached = if value <= 2147483647 {
                Ok(Ordering::Equal)
            } else {
                Err(overflow(0, value))
            };
            assert_eq!(
                compare_canon_atoms_for_kekulize(
                    &mut atoms,
                    0,
                    1,
                    CanonCompareMode::Atom,
                    flags(true)
                ),
                reached
            );
        }
        let source = [
            atom(0, PropertyValue::String("bad".into())),
            atom(1, PropertyValue::Int(1)),
            atom(2, PropertyValue::UInt(4294967295)),
            atom(3, PropertyValue::UInt(2147483648)),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut [0; 4],
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(overflow(2, 4294967295))
        );
    }
}

#[cfg(test)]
mod source_selected_scalar_transport {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Element, TopologyBlock};
    fn graph() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..3)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn selected_stale_explicit_fields_and_unselected_sentinels_survive_nonaromatic_return() {
        let graph = graph();
        for (stored, explicit, implicit) in [
            (1, 1, 3),
            (3, 3, 1),
            (66, 66, 0),
            (-1, 0, 4),
            (255, 0, 4),
            (256, 0, 4),
        ] {
            let cache = ValenceAssignment {
                explicit_valence: vec![stored, -128, 77],
                implicit_hydrogens: vec![-1, 255, 22],
            };
            let before = cache.clone();
            let out = kekulize_fragment(
                &graph,
                &[true, false, false],
                &[],
                &KekulizeParams::default(),
                None,
                None,
                Some(&cache),
            )
            .unwrap();
            assert_eq!(out.topology, graph);
            assert_eq!(
                out.final_valence,
                Some(ValenceAssignment {
                    explicit_valence: vec![explicit, -128, 77],
                    implicit_hydrogens: vec![implicit, 255, 22]
                })
            );
            assert_eq!(cache, before);
            assert!(out.ring_update.is_none());
        }
    }
    #[test]
    fn atoms_none_preserves_actual_cache_and_selected_negative_storage_errors() {
        let graph = graph();
        let cache = ValenceAssignment {
            explicit_valence: vec![128, -2, 66],
            implicit_hydrogens: vec![-1, -1, 17],
        };
        let before = cache.clone();
        let none = kekulize_fragment(
            &graph,
            &[false; 3],
            &[],
            &KekulizeParams::default(),
            None,
            None,
            Some(&cache),
        )
        .unwrap();
        assert_eq!(none.final_valence, Some(before.clone()));
        assert!(
            matches!(kekulize_fragment(&graph, &[true, false, false], &[], &KekulizeParams::default(), None, None, Some(&cache)), Err(KekulizeError::Valence(ValenceError::ExplicitValenceCacheNotInitialized { atom })) if atom == AtomId::new(0))
        );
        assert_eq!(cache, before);
        assert!(
            kekulize_fragment(
                &graph,
                &[false; 3],
                &[],
                &KekulizeParams::default(),
                None,
                None,
                None
            )
            .unwrap()
            .final_valence
            .is_none()
        );
    }
}

#[cfg(test)]
mod source_prepared_fragment_hydrogen_getter {
    use super::*;
    #[test]
    fn no_implicit_and_signed_width_are_observed_before_fragment_ranking() {
        for (no_implicit, field) in [
            (true, 77),
            (true, -1),
            (true, 128),
            (false, 256),
            (false, -256),
        ] {
            let topology = TopologyBlock::try_from_parts(
                vec![
                    Atom::from_spec(
                        AtomId::new(0),
                        cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::C)
                            .with_no_implicit(no_implicit),
                    ),
                    Atom::from_spec(
                        AtomId::new(1),
                        cosmolkit_model::AtomSpec::new(cosmolkit_model::Element::C),
                    ),
                ],
                vec![],
                vec![],
                vec![],
            )
            .unwrap();
            let value = ValenceAssignment {
                explicit_valence: vec![0; 2],
                implicit_hydrogens: vec![field, 0],
            };
            let before = value.clone();
            let view =
                CanonRankReadView::from_prepared_state(&topology, &value, None, Some(&[true; 2]))
                    .unwrap();
            let atoms =
                init_fragment_canon_atoms(&view, &[true; 2], &[], false, None, None).unwrap();
            assert_eq!(
                atoms
                    .iter()
                    .map(|atom| atom.total_num_hs)
                    .collect::<Vec<_>>(),
                vec![0; 2]
            );
            // Source BreakTies distinguishes identical disconnected atoms.
            assert_eq!(
                rank_fragment_atoms_with_prepared_state(
                    &topology,
                    &value,
                    None,
                    &[true; 2],
                    &[],
                    None,
                    None,
                    &CanonicalRankParams::kekulize_fragment_default()
                )
                .unwrap(),
                vec![0, 1]
            );
            let mut symmetric = CanonicalRankParams::kekulize_fragment_default();
            symmetric.break_ties = false;
            assert_eq!(
                rank_fragment_atoms_with_prepared_state(
                    &topology,
                    &value,
                    None,
                    &[true; 2],
                    &[],
                    None,
                    None,
                    &symmetric
                )
                .unwrap(),
                vec![0; 2]
            );
            assert_eq!(value, before);
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, PropertyValueKind};
    use cosmolkit_types::Element;
    fn atom(i: usize, p: Option<PropertyValue>) -> Atom {
        let mut a = Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C));
        if let Some(p) = p {
            a.set_computed_prop("_CanonicalRankingNumber", p).unwrap();
        }
        a
    }
    fn flags(b: bool) -> CanonRankFlags {
        let mut p = CanonicalRankParams::default();
        p.use_non_stereo_ranks = b;
        CanonRankFlags::from_fragment_options(p)
    }
    fn graph(v: Vec<Atom>) -> TopologyBlock {
        TopologyBlock::try_from_parts(v, vec![], vec![], vec![]).unwrap()
    }
    fn overflow(i: usize, v: u32) -> CanonicalRankError {
        CanonicalRankError::UnsignedRankOverflow {
            atom_index: i,
            property: "_CanonicalRankingNumber",
            value: v,
        }
    }
    fn bad(i: usize, k: PropertyValueKind) -> CanonicalRankError {
        CanonicalRankError::InvalidPropertyKind {
            atom_index: i,
            property: "_CanonicalRankingNumber",
            kind: k,
        }
    }

    // FROZEN UINT CONDITION: RANK_SIGNED_0_SIDE0
    #[test]
    fn uint_cell_rank_signed_0_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::Int(0))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_0_SIDE1
    #[test]
    fn uint_cell_rank_signed_0_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(0))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_1_SIDE0
    #[test]
    fn uint_cell_rank_signed_1_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::Int(1))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_1_SIDE1
    #[test]
    fn uint_cell_rank_signed_1_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(1))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483646_SIDE0
    #[test]
    fn uint_cell_rank_signed_2147483646_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::Int(2147483646))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483646_SIDE1
    #[test]
    fn uint_cell_rank_signed_2147483646_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(2147483646))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483647_SIDE0
    #[test]
    fn uint_cell_rank_signed_2147483647_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::Int(2147483647))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483647_SIDE1
    #[test]
    fn uint_cell_rank_signed_2147483647_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483648_SIDE0
    #[test]
    fn uint_cell_rank_signed_2147483648_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::Int(0))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(
            rank_mol_atoms_with_params(&g, &p),
            Err(overflow(0, 2147483648_u32))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_2147483648_SIDE1
    #[test]
    fn uint_cell_rank_signed_2147483648_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(0))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(
            rank_mol_atoms_with_params(&g, &p),
            Err(overflow(1, 2147483648_u32))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_4294967295_SIDE0
    #[test]
    fn uint_cell_rank_signed_4294967295_side0() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::Int(0))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(
            rank_mol_atoms_with_params(&g, &p),
            Err(overflow(0, 4294967295_u32))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: RANK_SIGNED_4294967295_SIDE1
    #[test]
    fn uint_cell_rank_signed_4294967295_side1() {
        let g = graph(vec![
            atom(0, Some(PropertyValue::Int(0))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        assert_eq!(
            rank_mol_atoms_with_params(&g, &p),
            Err(overflow(1, 4294967295_u32))
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_0
    #[test]
    fn uint_cell_uint_guard_class_lt_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_0
    #[test]
    fn uint_cell_uint_guard_class_gt_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_0
    #[test]
    fn uint_cell_uint_guard_flag_false_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_0
    #[test]
    fn uint_cell_uint_guard_left_first_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_0
    #[test]
    fn uint_cell_uint_guard_right_after_left_0() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_0
    #[test]
    fn uint_cell_uint_guard_missing_left_0() {
        let source = [atom(0, None), atom(1, Some(PropertyValue::UInt(0_u32)))];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_0
    #[test]
    fn uint_cell_uint_guard_mask_none_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_0
    #[test]
    fn uint_cell_uint_guard_mask_one_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_0
    #[test]
    fn uint_cell_uint_guard_special_modes_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_0
    #[test]
    fn uint_cell_uint_guard_hanoi_one_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_0
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_0
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_0() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_0
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_0() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(0_u32))),
            atom(3, Some(PropertyValue::UInt(0_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(bad(0, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_0
    #[test]
    fn uint_cell_uint_guard_public_singleton_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(0_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_0
    #[test]
    fn uint_cell_uint_guard_public_flag_false_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_0
    #[test]
    fn uint_cell_uint_guard_public_fragment_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_0
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_0
    #[test]
    fn uint_cell_uint_guard_empty_whole_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(0_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(0_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_0
    #[test]
    fn uint_cell_uint_guard_mask_validation_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_0
    #[test]
    fn uint_cell_uint_guard_prepared_validation_0() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::UInt(0_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_0_String("bad")
    #[test]
    fn uint_cell_getter_order_0_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_0_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_0_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_0_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_1_String("bad")
    #[test]
    fn uint_cell_getter_order_0_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_0_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_0_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_0_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_0_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(0_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_1
    #[test]
    fn uint_cell_uint_guard_class_lt_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_1
    #[test]
    fn uint_cell_uint_guard_class_gt_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_1
    #[test]
    fn uint_cell_uint_guard_flag_false_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_1
    #[test]
    fn uint_cell_uint_guard_left_first_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_1
    #[test]
    fn uint_cell_uint_guard_right_after_left_1() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_1
    #[test]
    fn uint_cell_uint_guard_missing_left_1() {
        let source = [atom(0, None), atom(1, Some(PropertyValue::UInt(1_u32)))];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_1
    #[test]
    fn uint_cell_uint_guard_mask_none_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_1
    #[test]
    fn uint_cell_uint_guard_mask_one_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_1
    #[test]
    fn uint_cell_uint_guard_special_modes_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_1
    #[test]
    fn uint_cell_uint_guard_hanoi_one_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_1
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_1
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_1() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_1
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_1() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(1_u32))),
            atom(3, Some(PropertyValue::UInt(1_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(bad(0, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_1
    #[test]
    fn uint_cell_uint_guard_public_singleton_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(1_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_1
    #[test]
    fn uint_cell_uint_guard_public_flag_false_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_1
    #[test]
    fn uint_cell_uint_guard_public_fragment_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_1
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_1
    #[test]
    fn uint_cell_uint_guard_empty_whole_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(1_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(1_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_1
    #[test]
    fn uint_cell_uint_guard_mask_validation_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_1
    #[test]
    fn uint_cell_uint_guard_prepared_validation_1() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::UInt(1_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_0_String("bad")
    #[test]
    fn uint_cell_getter_order_1_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_1_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_1_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_1_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_1_String("bad")
    #[test]
    fn uint_cell_getter_order_1_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_1_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_1_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_1_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_1_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(1_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_2147483646
    #[test]
    fn uint_cell_uint_guard_class_lt_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_2147483646
    #[test]
    fn uint_cell_uint_guard_class_gt_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_2147483646
    #[test]
    fn uint_cell_uint_guard_flag_false_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_2147483646
    #[test]
    fn uint_cell_uint_guard_left_first_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_2147483646
    #[test]
    fn uint_cell_uint_guard_right_after_left_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_2147483646
    #[test]
    fn uint_cell_uint_guard_missing_left_2147483646() {
        let source = [
            atom(0, None),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_2147483646
    #[test]
    fn uint_cell_uint_guard_mask_none_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_2147483646
    #[test]
    fn uint_cell_uint_guard_mask_one_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_2147483646
    #[test]
    fn uint_cell_uint_guard_special_modes_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_2147483646
    #[test]
    fn uint_cell_uint_guard_hanoi_one_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_2147483646
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_2147483646
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_2147483646
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_2147483646() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(2147483646_u32))),
            atom(3, Some(PropertyValue::UInt(2147483646_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(bad(0, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_2147483646
    #[test]
    fn uint_cell_uint_guard_public_singleton_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(2147483646_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_2147483646
    #[test]
    fn uint_cell_uint_guard_public_flag_false_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_2147483646
    #[test]
    fn uint_cell_uint_guard_public_fragment_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_2147483646
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_2147483646
    #[test]
    fn uint_cell_uint_guard_empty_whole_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(2147483646_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(2147483646_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_2147483646
    #[test]
    fn uint_cell_uint_guard_mask_validation_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_2147483646
    #[test]
    fn uint_cell_uint_guard_prepared_validation_2147483646() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::UInt(2147483646_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_0_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483646_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483646_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483646_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483646_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_1_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483646_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483646_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483646_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483646_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483646_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483646_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_2147483647
    #[test]
    fn uint_cell_uint_guard_class_lt_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_2147483647
    #[test]
    fn uint_cell_uint_guard_class_gt_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_2147483647
    #[test]
    fn uint_cell_uint_guard_flag_false_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_2147483647
    #[test]
    fn uint_cell_uint_guard_left_first_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_2147483647
    #[test]
    fn uint_cell_uint_guard_right_after_left_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_2147483647
    #[test]
    fn uint_cell_uint_guard_missing_left_2147483647() {
        let source = [
            atom(0, None),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_2147483647
    #[test]
    fn uint_cell_uint_guard_mask_none_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_2147483647
    #[test]
    fn uint_cell_uint_guard_mask_one_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_2147483647
    #[test]
    fn uint_cell_uint_guard_special_modes_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_2147483647
    #[test]
    fn uint_cell_uint_guard_hanoi_one_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_2147483647
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_2147483647
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_2147483647
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_2147483647() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(2147483647_u32))),
            atom(3, Some(PropertyValue::UInt(2147483647_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(bad(0, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_2147483647
    #[test]
    fn uint_cell_uint_guard_public_singleton_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(2147483647_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_2147483647
    #[test]
    fn uint_cell_uint_guard_public_flag_false_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_2147483647
    #[test]
    fn uint_cell_uint_guard_public_fragment_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_2147483647
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_2147483647
    #[test]
    fn uint_cell_uint_guard_empty_whole_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(2147483647_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(2147483647_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_2147483647
    #[test]
    fn uint_cell_uint_guard_mask_validation_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_2147483647
    #[test]
    fn uint_cell_uint_guard_prepared_validation_2147483647() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::UInt(2147483647_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_0_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483647_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483647_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483647_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483647_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_1_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483647_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483647_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483647_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483647_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483647_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483647_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_2147483648
    #[test]
    fn uint_cell_uint_guard_class_lt_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_2147483648
    #[test]
    fn uint_cell_uint_guard_class_gt_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_2147483648
    #[test]
    fn uint_cell_uint_guard_flag_false_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_2147483648
    #[test]
    fn uint_cell_uint_guard_left_first_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_2147483648
    #[test]
    fn uint_cell_uint_guard_right_after_left_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(1, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_2147483648
    #[test]
    fn uint_cell_uint_guard_missing_left_2147483648() {
        let source = [
            atom(0, None),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(1, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_2147483648
    #[test]
    fn uint_cell_uint_guard_mask_none_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_2147483648
    #[test]
    fn uint_cell_uint_guard_mask_one_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_2147483648
    #[test]
    fn uint_cell_uint_guard_special_modes_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_2147483648
    #[test]
    fn uint_cell_uint_guard_hanoi_one_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_2147483648
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_2147483648
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_2147483648
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_2147483648() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(2147483648_u32))),
            atom(3, Some(PropertyValue::UInt(2147483648_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(overflow(2, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_2147483648
    #[test]
    fn uint_cell_uint_guard_public_singleton_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(2147483648_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_2147483648
    #[test]
    fn uint_cell_uint_guard_public_flag_false_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_2147483648
    #[test]
    fn uint_cell_uint_guard_public_fragment_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_2147483648
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_2147483648
    #[test]
    fn uint_cell_uint_guard_empty_whole_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(2147483648_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_2147483648
    #[test]
    fn uint_cell_uint_guard_mask_validation_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_2147483648
    #[test]
    fn uint_cell_uint_guard_prepared_validation_2147483648() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::UInt(2147483648_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_0_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483648_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483648_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483648_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483648_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 2147483648_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_1_String("bad")
    #[test]
    fn uint_cell_getter_order_2147483648_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_2147483648_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_2147483648_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_2147483648_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_2147483648_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(2147483648_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_LT_4294967295
    #[test]
    fn uint_cell_uint_guard_class_lt_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 0;
        atoms[1].index = 1;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Less)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_CLASS_GT_4294967295
    #[test]
    fn uint_cell_uint_guard_class_gt_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].index = 1;
        atoms[1].index = 0;
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Ok(Ordering::Greater)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_FLAG_FALSE_4294967295
    #[test]
    fn uint_cell_uint_guard_flag_false_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(false)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_LEFT_FIRST_4294967295
    #[test]
    fn uint_cell_uint_guard_left_first_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_RIGHT_AFTER_LEFT_4294967295
    #[test]
    fn uint_cell_uint_guard_right_after_left_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::Int(2147483647))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(1, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MISSING_LEFT_4294967295
    #[test]
    fn uint_cell_uint_guard_missing_left_4294967295() {
        let source = [
            atom(0, None),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(1, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_NONE_4294967295
    #[test]
    fn uint_cell_uint_guard_mask_none_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = false;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Ok(Ordering::Equal)
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_ONE_4294967295
    #[test]
    fn uint_cell_uint_guard_mask_one_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        atoms[0].is_in_play = false;
        atoms[1].is_in_play = true;
        assert_eq!(
            compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, CanonCompareMode::Atom, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_SPECIAL_MODES_4294967295
    #[test]
    fn uint_cell_uint_guard_special_modes_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        for mode in [
            CanonCompareMode::SpecialChirality,
            CanonCompareMode::SpecialSymmetry,
        ] {
            assert_eq!(
                compare_canon_atoms_for_kekulize(&mut atoms, 0, 1, mode, flags(true)),
                Ok(Ordering::Equal)
            );
        }
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_ONE_4294967295
    #[test]
    fn uint_cell_uint_guard_hanoi_one_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0],
                &mut [0; 1],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [1, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_UNCHANGED_4294967295
    #[test]
    fn uint_cell_uint_guard_hanoi_unchanged_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[false, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Ok(false)
        );
        assert_eq!(counts, [2, 0]);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_CHANGED_4294967295
    #[test]
    fn uint_cell_uint_guard_hanoi_changed_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 2];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [0, 1],
                &mut [0; 2],
                &mut counts,
                &[true, false],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_HANOI_RECURSION_4294967295
    #[test]
    fn uint_cell_uint_guard_hanoi_recursion_4294967295() {
        let source = [
            atom(0, Some(PropertyValue::String("bad".into()))),
            atom(1, Some(PropertyValue::Int(1))),
            atom(2, Some(PropertyValue::UInt(4294967295_u32))),
            atom(3, Some(PropertyValue::UInt(4294967295_u32))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        for a in &mut atoms {
            a.index = 0;
        }
        let mut counts = [0; 4];
        assert_eq!(
            hanoi_order_for_kekulize(
                &mut [2, 3, 0, 1],
                &mut [0; 4],
                &mut counts,
                &[true; 4],
                &mut atoms,
                CanonCompareMode::Atom,
                flags(true)
            ),
            Err(overflow(2, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_SINGLETON_4294967295
    #[test]
    fn uint_cell_uint_guard_public_singleton_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![atom(0, Some(PropertyValue::UInt(4294967295_u32)))]);
        let before = g.clone();
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FLAG_FALSE_4294967295
    #[test]
    fn uint_cell_uint_guard_public_flag_false_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        p.use_non_stereo_ranks = false;
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![0, 0]));
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PUBLIC_FRAGMENT_4294967295
    #[test]
    fn uint_cell_uint_guard_public_fragment_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true, true], &[], None, None, &p),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_FRAGMENT_4294967295
    #[test]
    fn uint_cell_uint_guard_prepared_fragment_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Ok(vec![0, 0])
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_EMPTY_WHOLE_4294967295
    #[test]
    fn uint_cell_uint_guard_empty_whole_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![]);
        assert_eq!(rank_mol_atoms_with_params(&g, &p), Ok(vec![]));
        let populated = graph(vec![atom(0, Some(PropertyValue::UInt(4294967295_u32)))]);
        assert_eq!(
            populated.atoms[0].prop("_CanonicalRankingNumber"),
            Some(&PropertyValue::UInt(4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: UINT_GUARD_MASK_VALIDATION_4294967295
    #[test]
    fn uint_cell_uint_guard_mask_validation_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        assert_eq!(
            rank_fragment_atoms_with_params(&g, &[true], &[], None, None, &p),
            Err(CanonicalRankError::AtomMaskLength {
                expected: 2,
                actual: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: UINT_GUARD_PREPARED_VALIDATION_4294967295
    #[test]
    fn uint_cell_uint_guard_prepared_validation_4294967295() {
        let mut p = CanonicalRankParams::default();
        p.break_ties = false;
        p.use_non_stereo_ranks = true;
        let g = graph(vec![
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::UInt(4294967295_u32))),
        ]);
        let before = g.clone();
        let mut valence = crate::assign_valence(&g, &crate::ValenceParams::default()).unwrap();
        valence.implicit_hydrogens.pop();
        assert_eq!(
            rank_fragment_atoms_with_prepared_state(
                &g,
                &valence,
                None,
                &[true, true],
                &[],
                None,
                None,
                &p
            ),
            Err(CanonicalRankError::PreparedValenceLength {
                atom_count: 2,
                explicit_len: 2,
                implicit_len: 1
            })
        );
        assert_eq!(g, before);
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_0_String("bad")
    #[test]
    fn uint_cell_getter_order_4294967295_0_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_0_Bool(false)
    #[test]
    fn uint_cell_getter_order_4294967295_0_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_0_Double(1.0)
    #[test]
    fn uint_cell_getter_order_4294967295_0_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_0_IntVector([])
    #[test]
    fn uint_cell_getter_order_4294967295_0_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 0, 1, flags(true)),
            Err(overflow(0, 4294967295_u32))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_1_String("bad")
    #[test]
    fn uint_cell_getter_order_4294967295_1_string__bad__() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::String("bad".into()))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::String))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_1_Bool(false)
    #[test]
    fn uint_cell_getter_order_4294967295_1_bool_false_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::Bool(false))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Bool))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_1_Double(1.0)
    #[test]
    fn uint_cell_getter_order_4294967295_1_double_1_0_() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::Double(1.0))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::Double))
        );
    }
    // FROZEN UINT CONDITION: GETTER_ORDER_4294967295_1_IntVector([])
    #[test]
    fn uint_cell_getter_order_4294967295_1_intvector____() {
        let source = [
            atom(0, Some(PropertyValue::UInt(4294967295_u32))),
            atom(1, Some(PropertyValue::IntVector(vec![]))),
        ];
        let mut atoms = source
            .iter()
            .map(empty_canon_atom_from_source_atom)
            .collect::<Vec<_>>();
        // Frozen input uses the same current class on both sides.
        for atom in &mut atoms {
            atom.index = 0;
        }
        assert_eq!(
            compare_canon_atom_base_for_kekulize(&atoms, 1, 0, flags(true)),
            Err(bad(1, PropertyValueKind::IntVector))
        );
    }
}

#[cfg(test)]
mod recovery_chem16 {
    use super::*;
    fn enumerator(n: usize) -> QuestionEnumerator {
        QuestionEnumerator::new((0..n).map(AtomId::new).collect()).unwrap()
    }
    #[test]
    fn zero_question_first_next_finishes_and_repeats_empty() {
        let mut e = enumerator(0);
        assert!(!e.done);
        assert_eq!(e.state.len(), 0);
        assert!(e.next().is_empty());
        assert!(e.done);
        assert!(e.next().is_empty());
    }
    #[test]
    fn complete_small_binary_order_matches_source_subsets() {
        for n in 1..=6 {
            let mut e = enumerator(n);
            for bits in 1usize..(1usize << n) {
                let expected = (0..n)
                    .filter(|i| bits & (1 << i) != 0)
                    .map(AtomId::new)
                    .collect::<Vec<_>>();
                assert_eq!(e.next(), expected);
                assert_eq!(e.done, bits == (1 << n) - 1);
            }
            assert!(e.next().is_empty());
            assert!(e.next().is_empty());
        }
    }
    #[test]
    fn wide_questions_start_lazily_without_integer_width_rejection() {
        for n in [31, 32, 33, 63, 64, 65, 130] {
            let mut e = enumerator(n);
            assert_eq!(e.state.len(), n.div_ceil(64));
            assert_eq!(e.next(), vec![AtomId::new(0)]);
            assert_eq!(e.next(), vec![AtomId::new(1)]);
            assert_eq!(e.next(), vec![AtomId::new(0), AtomId::new(1)]);
            assert!(!e.done);
        }
    }
    #[test]
    fn packed_carry_crosses_word_boundary_and_terminal_tail() {
        let mut e = enumerator(130);
        e.state = vec![u64::MAX, 0, 0];
        assert_eq!(e.next(), (0..64).map(AtomId::new).collect::<Vec<_>>());
        assert_eq!(e.state, vec![0, 1, 0]);
        assert_eq!(e.next(), vec![AtomId::new(64)]);
        e.state = vec![u64::MAX, u64::MAX, 1];
        assert_eq!(e.next(), (0..129).map(AtomId::new).collect::<Vec<_>>());
        assert_eq!(e.state, vec![0, 0, 2]);
        e.state = vec![u64::MAX, u64::MAX, 3];
        assert_eq!(e.next(), (0..130).map(AtomId::new).collect::<Vec<_>>());
        assert!(e.done);
        assert_eq!(e.state, vec![0, 0, 0]);
        assert!(e.next().is_empty());
    }
}
#[cfg(test)]
mod recovery_chem17 {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn graph(specs: Vec<AtomSpec>, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            specs
                .into_iter()
                .enumerate()
                .map(|(i, s)| Atom::from_spec(AtomId::new(i), s))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b, o))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn ring(n: usize) -> TopologyBlock {
        let mut t = graph(
            vec![AtomSpec::new(Element::C).with_aromatic(true); n],
            &(0..n)
                .map(|i| (i, (i + 1) % n, BondOrder::Aromatic))
                .collect::<Vec<_>>(),
        );
        for b in &mut t.bonds {
            b.set_aromatic(true);
        }
        t
    }
    fn cold(n: usize) -> ValenceAssignment {
        ValenceAssignment {
            explicit_valence: vec![-1; n],
            implicit_hydrogens: vec![-1; n],
        }
    }
    fn reset() -> RingInfo {
        let mut r = RingInfo::new(RingFindType::OtherOrUnknown, 0, 0);
        r.reset();
        r
    }
    #[test]
    fn direct_nonaromatic_attempt_retains_cache_without_manufacturing_rings() {
        let mut t = graph(vec![AtomSpec::new(Element::C)], &[]);
        let before = t.clone();
        let mut v = ValenceAssignment {
            explicit_valence: vec![3],
            implicit_hydrogens: vec![-1],
        };
        let mut r = reset();
        source_kekulize_attempt(&mut t, &mut v, &mut r, &KekulizeParams::default()).unwrap();
        assert_eq!(t, before);
        assert_eq!(v.explicit_valence, vec![3]);
        assert_eq!(v.implicit_hydrogens, vec![1]);
        assert!(!r.is_initialized());
    }
    #[test]
    fn atoms_none_precedes_cache_length_and_preserves_all_three_states() {
        let mut t = graph(vec![AtomSpec::new(Element::C)], &[]);
        let before = t.clone();
        let mut v = cold(0);
        let old = v.clone();
        let mut r = reset();
        let mut rows = KekulizeRingRows::Live(&mut r);
        let result = kekulize_fragment_attempt(
            &mut t,
            &mut v,
            &mut rows,
            &[false],
            &[],
            &KekulizeParams::default(),
            None,
        )
        .unwrap();
        assert_eq!(result, (vec![], false));
        assert_eq!(t, before);
        assert_eq!(v, old);
        assert!(!r.is_initialized());
    }
    #[test]
    fn earlier_cache_rows_survive_later_stale_explicit_getter_error() {
        let t = graph(vec![AtomSpec::new(Element::C); 2], &[]);
        let mut v = ValenceAssignment {
            explicit_valence: vec![-1, -2],
            implicit_hydrogens: vec![-1, -1],
        };
        let err = prepare_kekulize_inputs(&t, &[true; 2], &[], None, &mut v).unwrap_err();
        assert!(
            matches!(err,KekulizeError::Valence(ValenceError::ExplicitValenceCacheNotInitialized{atom}) if atom==AtomId::new(1))
        );
        assert_eq!(v.explicit_valence, vec![0, -2]);
        assert_eq!(v.implicit_hydrogens, vec![4, 6]);
    }
    #[test]
    fn failed_matching_keeps_direction_clears_and_single_reset_in_live_graph() {
        let mut t = ring(5);
        for b in &mut t.bonds {
            b.set_direction(BondDirection::EndUpRight);
        }
        let before = t.clone();
        let mut v = cold(5);
        let mut r = reset();
        let err = source_kekulize_attempt(
            &mut t,
            &mut v,
            &mut r,
            &KekulizeParams {
                mark_atoms_bonds: false,
                canonical: false,
                max_backtracks: 0,
            },
        )
        .unwrap_err();
        assert!(matches!(err, KekulizeError::NotKekulizable { .. }));
        assert!(
            t.bonds
                .iter()
                .all(|b| b.order() == BondOrder::Single && b.is_aromatic())
        );
        assert!(t.bonds.iter().any(|b| b.direction() == BondDirection::None));
        assert!(t.atoms.iter().all(Atom::is_aromatic));
        assert!(r.is_initialized());
        assert_eq!(r.find_type(), RingFindType::Sssr);
        assert!(v.explicit_valence.iter().all(|&x| x == 3));
        assert!(v.implicit_hydrogens.iter().all(|&x| x == 1));
        assert_ne!(t, before);
    }
    #[test]
    fn readonly_failure_keeps_original_graph_cache_and_borrowed_ring_rows() {
        let t = ring(5);
        let v = cold(5);
        let r = crate::find_sssr_from_parts(5, &t.bonds, &t.adjacency).unwrap();
        let before = (t.clone(), v.clone(), r.clone());
        assert!(
            kekulize_with_query_state_and_ring_info(
                &t,
                &KekulizeParams {
                    canonical: false,
                    ..KekulizeParams::default()
                },
                None,
                Some(&r),
                Some(&v)
            )
            .is_err()
        );
        assert_eq!((t, v, r), before);
    }
    #[test]
    fn source_recovery_restores_membership_not_old_order_hydrogen_or_direction() {
        let mut t = ring(6);
        t.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1),
        );
        for b in &mut t.bonds {
            b.set_order(BondOrder::Single);
            b.set_direction(BondDirection::EndUpRight);
        }
        let aromatic_atoms = t.atoms.iter().map(Atom::is_aromatic).collect::<Vec<_>>();
        let aromatic_bonds = vec![true; 6];
        t.atoms[0].set_aromatic(false);
        t.atoms[0].set_no_implicit(false);
        t.atoms[0].set_explicit_hydrogens(0);
        t.bonds[0].set_aromatic(false);
        t.bonds[0].set_order(BondOrder::Double);
        t.bonds[0].set_direction(BondDirection::None);
        restore_source_aromatic_membership(&mut t, &aromatic_atoms, &aromatic_bonds);
        assert!(t.atoms.iter().all(Atom::is_aromatic));
        assert!(
            t.bonds
                .iter()
                .all(|b| b.is_aromatic() && b.order() == BondOrder::Aromatic)
        );
        assert_eq!(t.atoms[0].explicit_hydrogens(), 0);
        assert!(!t.atoms[0].no_implicit());
        assert_eq!(t.bonds[0].direction(), BondDirection::None);
    }
    #[test]
    fn late_np_and_stray_aromatic_atom_failure_exposes_source_prefix() {
        // Valid detached graph for source-derived [nH]1cccc1.c before sanitize.
        let mut t = ring(5);
        t.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::N)
                .with_aromatic(true)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1),
        );
        t.atoms.push(Atom::from_spec(
            AtomId::new(5),
            AtomSpec::new(Element::C).with_aromatic(true),
        ));
        t.adjacency = AdjacencyList::try_from_topology(6, &t.bonds).unwrap();
        t.validate().unwrap();
        let aromatic_atoms = vec![true; 6];
        let aromatic_bonds = vec![true; 5];
        let mut v = cold(6);
        let mut r = reset();
        let error = source_kekulize_attempt(
            &mut t,
            &mut v,
            &mut r,
            &KekulizeParams {
                canonical: false,
                ..KekulizeParams::default()
            },
        )
        .unwrap_err();
        assert!(
            matches!(error,KekulizeError::AromaticAtomOutsideRing{atom} if atom==AtomId::new(5))
        );
        assert_eq!(t.atoms[0].explicit_hydrogens(), 0);
        assert!(!t.atoms[0].no_implicit());
        assert!(!t.atoms[0].is_aromatic());
        assert_eq!(v.explicit_valence[0], 2);
        assert_eq!(v.implicit_hydrogens[0], 1);
        assert!(t.bonds.iter().all(|b| !b.is_aromatic()));
        restore_source_aromatic_membership(&mut t, &aromatic_atoms, &aromatic_bonds);
        assert!(t.atoms[0].is_aromatic());
        assert_eq!(t.atoms[0].explicit_hydrogens(), 0);
        assert_eq!(v.explicit_valence[0], 2);
        assert!(r.is_initialized());
    }
    #[test]
    fn borrowed_invalid_dimensions_keep_prepared_valence_error_priority() {
        let t = ring(6);
        let v = cold(6);
        let r = RingInfo::new(RingFindType::OtherOrUnknown, 7, 6);
        let snapshot = r.clone();
        let mut rows = KekulizeRingRows::Borrowed(&r);
        let error = kekulize_ring_state_transition_mut(
            &t, &v, &[true; 6], &[true; 6], true, true, &mut rows,
        )
        .unwrap_err();
        assert!(matches!(
            error,
            KekulizeError::CanonicalRank(CanonicalRankError::PreparedValenceInvalid {
                atom_index: 0
            })
        ));
        assert_eq!(r, snapshot);
    }

    #[test]
    fn rank_error_retains_published_fast_and_success_replaces_with_sssr() {
        let t = ring(6);
        let mut invalid = cold(6);
        let mut r = reset();
        {
            let mut rows = KekulizeRingRows::Live(&mut r);
            assert!(
                kekulize_ring_state_transition_mut(
                    &t, &invalid, &[true; 6], &[true; 6], true, true, &mut rows
                )
                .is_err()
            );
        }
        assert!(r.is_find_fast_or_better());
        assert_eq!(r.find_type(), RingFindType::Fast);
        invalid =
            crate::assign_valence_with_options_for_topology(&t, ValenceModel::RdkitLike, false)
                .unwrap();
        let mut r = reset();
        let mut rows = KekulizeRingRows::Live(&mut r);
        kekulize_ring_state_transition_mut(
            &t, &invalid, &[true; 6], &[true; 6], true, true, &mut rows,
        )
        .unwrap();
        assert_eq!(r.find_type(), RingFindType::Sssr);
        let ptr = r.atom_rings().as_ptr();
        let mut rows = KekulizeRingRows::Live(&mut r);
        kekulize_ring_state_transition_mut(
            &t, &invalid, &[true; 6], &[true; 6], true, true, &mut rows,
        )
        .unwrap();
        assert_eq!(r.atom_rings().as_ptr(), ptr);
    }
    #[test]
    fn retry_only_true_source_sanitize_errors_and_preserve_final_failed_state() {
        assert!(source_sanitize_kekulize_error(
            &KekulizeError::NotKekulizable {
                problem_atoms: vec![]
            }
        ));
        assert!(source_sanitize_kekulize_error(
            &KekulizeError::AromaticAtomOutsideRing {
                atom: AtomId::new(0)
            }
        ));
        assert!(source_sanitize_kekulize_error(
            &KekulizeError::PostconditionValenceMismatch {
                atom: AtomId::new(0),
                before: 3,
                after: 4
            }
        ));
        assert!(!source_sanitize_kekulize_error(
            &KekulizeError::MatchingStateLength {
                field: "cache",
                expected: 1,
                actual: 0
            }
        ));
        assert!(!source_sanitize_kekulize_error(
            &KekulizeError::CanonicalRank(CanonicalRankError::PreparedValenceInvalid {
                atom_index: 0
            })
        ));
        assert!(!source_sanitize_kekulize_error(&KekulizeError::Valence(
            ValenceError::InvalidExplicitValenceInput {
                atom: AtomId::new(0),
                value: -2
            }
        )));
        let mut t = ring(5);
        let mut v = cold(5);
        let mut rings = None;
        assert!(matches!(
            source_kekulize_for_sanitize(&mut t, &mut v, &mut rings, None),
            Err(KekulizeError::NotKekulizable { .. })
        ));
        // Second canonical attempt has failed: first-only recovery did not
        // restore its final single bond state back to aromatic order.
        assert!(
            t.bonds
                .iter()
                .all(|b| b.order() == BondOrder::Single && b.is_aromatic())
        );
        assert!(rings.unwrap().is_initialized());
    }
}

// Private test-only scoped observations of scratch constructor/refinement calls. after rebasing CHEM17.
#[cfg(test)]
mod chem31_scratch_trace {
    use super::CanonCompareMode;
    use std::{cell::RefCell, marker::PhantomData, rc::Rc};
    #[derive(Debug, Clone)]
    pub(super) enum Event {
        Construct(bool, usize, usize), // top-level?, constructor length, data pointer
        Refine(CanonCompareMode, usize, usize), // actual mode, scratch pointer, length
        Sort(usize, usize, bool, usize), // base offset, length, result-in-temp, temp pointer
    }
    std::thread_local! {
        static EVENTS: RefCell<Option<Vec<Event>>> = const { RefCell::new(None) };
    }
    pub(super) struct Capture {
        _thread_bound: PhantomData<Rc<()>>,
    }
    pub(super) fn capture() -> Capture {
        EVENTS.with(|events| {
            let mut events = events.borrow_mut();
            assert!(events.is_none(), "CHEM31 scratch capture cannot nest");
            *events = Some(Vec::new());
        });
        Capture {
            _thread_bound: PhantomData,
        }
    }
    impl Capture {
        pub(super) fn events(&self) -> Vec<Event> {
            EVENTS.with(|events| events.borrow().as_ref().expect("active capture").clone())
        }
    }
    impl Drop for Capture {
        fn drop(&mut self) {
            EVENTS.with(|events| {
                events.borrow_mut().take();
            });
        }
    }
    pub(super) fn record(event: Event) {
        EVENTS.with(|events| {
            if let Some(events) = events.borrow_mut().as_mut() {
                events.push(event);
            }
        });
    }
}

#[cfg(test)]
mod recovery_chem31 {

    use super::*;
    use cosmolkit_model::{AtomSpec, BondSpec};
    use cosmolkit_types::Element;
    fn chem31_graph(n: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = (0..n)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap()
    }
    fn chem31_isolated(n: usize) -> TopologyBlock {
        chem31_graph(n, &[])
    }
    fn chem31_atoms<'a>(view: &CanonRankReadView<'a>) -> Vec<CanonAtom<'a>> {
        init_fragment_canon_atoms(
            view,
            &vec![true; view.num_atoms()],
            &vec![true; view.bonds.len()],
            true,
            None,
            None,
        )
        .unwrap()
    }
    fn chem31_flags() -> CanonRankFlags {
        CanonRankFlags::from_fragment_options(CanonicalRankParams::kekulize_fragment_default())
    }
    #[test]
    fn rdkit_2026_03_6_chem31_inactive_guard_preserves_all_state_and_dirty_supplied_scratch() {
        use super::chem31_scratch_trace::{Event, capture};
        let m = chem31_graph(2, &[(0, 1)]);
        let view = CanonRankReadView::from_topology(&m).unwrap();
        for scratch_len in [1, 2, 5] {
            let mut atoms = chem31_atoms(&view);
            atoms[0].index = 0;
            atoms[1].index = 1;
            let mut order = vec![0, 1];
            let mut count = vec![1, 1];
            let mut active = -1;
            let mut next = vec![-2; 2];
            let mut changed = vec![false; 2];
            let mut touched = vec![false; 2];
            let mut scratch = vec![usize::MAX; scratch_len];
            let before = (
                format!("{atoms:?}"),
                order.clone(),
                count.clone(),
                active,
                next.clone(),
                changed.clone(),
                touched.clone(),
                scratch.clone(),
            );
            let trace = capture();
            let result = super::refine_partitions_for_kekulize(
                &view,
                &mut atoms,
                true,
                CanonCompareMode::Atom,
                chem31_flags(),
                &mut order,
                &mut count,
                &mut active,
                &mut next,
                &mut changed,
                &mut touched,
                Some(&mut scratch),
            );
            if scratch_len < 2 {
                let error = result.unwrap_err();
                assert!(matches!(error, CanonicalRankError::HanoiScratchTooSmall));
                assert_eq!(error.to_string(), "hanoi scratch is too small");
                assert!(trace.events().is_empty());
            } else {
                result.unwrap();
                assert!(
                    matches!(trace.events().as_slice(),[Event::Refine(CanonCompareMode::Atom,_,n)] if *n==scratch_len)
                );
            }
            assert_eq!(
                (
                    format!("{atoms:?}"),
                    order,
                    count,
                    active,
                    next,
                    changed,
                    touched,
                    scratch
                ),
                before
            );
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem31_inactive_fallback_and_empty_constructor_are_distinct_from_heap_count()
    {
        use super::chem31_scratch_trace::{Event, capture};
        for n in [0, 2] {
            let m = chem31_isolated(n);
            let view = CanonRankReadView::from_topology(&m).unwrap();
            let mut atoms = chem31_atoms(&view);
            let mut order = (0..n).collect::<Vec<_>>();
            let mut count = vec![1; n];
            let mut active = -1;
            let mut next = vec![-2; n];
            let mut changed = vec![false; n];
            let mut touched = vec![false; n];
            let trace = capture();
            super::refine_partitions_for_kekulize(
                &view,
                &mut atoms,
                false,
                CanonCompareMode::Atom,
                chem31_flags(),
                &mut order,
                &mut count,
                &mut active,
                &mut next,
                &mut changed,
                &mut touched,
                None,
            )
            .unwrap();
            assert!(
                matches!(trace.events().as_slice(),[Event::Construct(false,len,p),Event::Refine(CanonCompareMode::Atom,r,len2)] if *len==n && *len2==n && p==r)
            );
            // A zero-length constructor invocation does not allocate a heap buffer.
        }
        let trace = capture();
        assert!(
            rank_mol_atoms(&TopologyBlock::default())
                .unwrap()
                .is_empty()
        );
        assert!(trace.events().is_empty());
    }

    #[test]
    fn rdkit_2026_03_6_chem31_nonzero_partition_uses_scratch_prefix_and_false_does_not_copy() {
        use super::chem31_scratch_trace::{Event, capture};
        for (n, symbols, expected, expect_temp, len) in [
            (
                5,
                vec!["0", "z", "x", "y", "4"],
                vec![0, 2, 3, 1, 4],
                true,
                3,
            ),
            (4, vec!["0", "z", "x", "3"], vec![0, 2, 1, 3], false, 2),
        ] {
            let m = chem31_isolated(n);
            let view = CanonRankReadView::from_topology(&m).unwrap();
            let mut atoms = chem31_atoms(&view);
            for i in 0..n {
                atoms[i].index = if i > 0 && i <= len { 1 } else { i as i32 };
                atoms[i].p_symbol = Some(symbols[i]);
            }
            let mut order = (0..n).collect::<Vec<_>>();
            let mut count = vec![1; n];
            count[1] = len;
            for c in &mut count[2..=len] {
                *c = 0;
            }
            let mut active = 1;
            let mut next = vec![-2; n];
            next[1] = -1;
            let mut changed = vec![true; n];
            let mut touched = vec![false; n];
            let mut scratch = vec![usize::MAX; n + 2];
            let ptr = scratch.as_ptr() as usize;
            let trace = capture();
            super::refine_partitions_for_kekulize(
                &view,
                &mut atoms,
                false,
                CanonCompareMode::Atom,
                chem31_flags(),
                &mut order,
                &mut count,
                &mut active,
                &mut next,
                &mut changed,
                &mut touched,
                Some(&mut scratch),
            )
            .unwrap();
            assert_eq!(order, expected);
            assert!(scratch[len..].iter().all(|&x| x == usize::MAX));
            assert!(trace.events().iter().any(
                |e| matches!(e,Event::Sort(1,k,temp,p) if *k==len && *temp==expect_temp && *p==ptr)
            ));
            if expect_temp {
                assert_eq!(&scratch[..len], &order[1..1 + len]);
            } else {
                assert!(scratch.iter().all(|&x| x == usize::MAX));
            }
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem31_breakties_none_degree_zero_skips_and_short_borrowed_keeps_partial_state()
     {
        use super::chem31_scratch_trace::capture;
        let m = chem31_isolated(2);
        let view = CanonRankReadView::from_topology(&m).unwrap();
        let mut atoms = chem31_atoms(&view);
        let mut order = vec![0, 1];
        let mut count = vec![2, 0];
        let mut active = -1;
        let mut next = vec![-2; 2];
        let mut changed = vec![false; 2];
        let mut touched = vec![false; 2];
        {
            let trace = capture();
            super::break_ties_for_kekulize(
                &view,
                &mut atoms,
                true,
                CanonCompareMode::Atom,
                chem31_flags(),
                &mut order,
                &mut count,
                &mut active,
                &mut next,
                &mut changed,
                &mut touched,
                None,
            )
            .unwrap();
            assert!(trace.events().is_empty());
        }
        assert_eq!(count, vec![1, 1]);
        assert_eq!(
            atoms.iter().map(|a| a.index).collect::<Vec<_>>(),
            vec![0, 1]
        );
        let m = chem31_graph(2, &[(0, 1)]);
        let view = CanonRankReadView::from_topology(&m).unwrap();
        let mut atoms = chem31_atoms(&view);
        for a in &mut atoms {
            a.index = 0;
        }
        let mut order = vec![0, 1];
        let mut count = vec![2, 0];
        let mut active = -1;
        let mut next = vec![-2; 2];
        let mut changed = vec![false; 2];
        let mut touched = vec![false; 2];
        let mut scratch = vec![usize::MAX; 1];
        let trace = capture();
        let error = super::break_ties_for_kekulize(
            &view,
            &mut atoms,
            true,
            CanonCompareMode::Atom,
            chem31_flags(),
            &mut order,
            &mut count,
            &mut active,
            &mut next,
            &mut changed,
            &mut touched,
            Some(&mut scratch),
        )
        .unwrap_err();
        assert!(matches!(error, CanonicalRankError::HanoiScratchTooSmall));
        assert_eq!(count, vec![1, 1]);
        assert_eq!(
            atoms.iter().map(|a| a.index).collect::<Vec<_>>(),
            vec![0, 1]
        );
        assert_eq!(changed, vec![true, false]);
        assert_eq!(scratch, vec![usize::MAX]);
        assert!(trace.events().is_empty());
        // BreakTies changes partition state before Refine validates; no rollback claim.
    }

    #[test]
    fn rdkit_2026_03_6_chem31_breakties_none_allocates_once_per_actual_refine() {
        use super::chem31_scratch_trace::{Event, capture};
        let m = chem31_graph(3, &[(0, 1), (1, 2), (2, 0)]);
        let view = CanonRankReadView::from_topology(&m).unwrap();
        let mut atoms = chem31_atoms(&view);
        for a in &mut atoms {
            a.index = 0;
        }
        let mut order = vec![0, 1, 2];
        let mut count = vec![3, 0, 0];
        let mut active = -1;
        let mut next = vec![-2; 3];
        let mut changed = vec![false; 3];
        let mut touched = vec![false; 3];
        let trace = capture();
        super::break_ties_for_kekulize(
            &view,
            &mut atoms,
            true,
            CanonCompareMode::Atom,
            chem31_flags(),
            &mut order,
            &mut count,
            &mut active,
            &mut next,
            &mut changed,
            &mut touched,
            None,
        )
        .unwrap();
        let events = trace.events();
        let fallback = events
            .iter()
            .filter(|e| matches!(e, Event::Construct(false, 3, _)))
            .count();
        let calls = events
            .iter()
            .filter(|e| matches!(e, Event::Refine(_, _, 3)))
            .count();
        assert!(calls >= 2, "{events:?}");
        assert_eq!(fallback, calls);
        assert!(
            !events
                .iter()
                .any(|e| matches!(e, Event::Construct(true, _, _)))
        );
    }

    #[test]
    fn rdkit_2026_03_6_chem31_top_rank_reuses_one_buffer_across_reached_modes_and_fragment_scope() {
        use super::chem31_scratch_trace::{Event, capture};
        // This regular cubane graph reaches the special-symmetry gate; the trace
        // assertion below observes actual routing rather than inferring it from ranks.
        let cube = chem31_graph(
            8,
            &[
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
            ],
        );
        for (m, expect_symmetry, fragment) in [
            (chem31_graph(2, &[(0, 1)]), false, false),
            (cube, true, false),
            (chem31_graph(3, &[(0, 1), (1, 2)]), false, true),
        ] {
            let trace = capture();
            let ranks = if fragment {
                let atoms = vec![true, false, true];
                let bonds = vec![false; m.bonds.len()];
                rank_fragment_atoms(&m, &atoms, &bonds).unwrap()
            } else {
                rank_mol_atoms(&m).unwrap()
            };
            assert_eq!(ranks.len(), m.atoms.len());
            let events = trace.events();
            let constructions = events
                .iter()
                .filter_map(|e| {
                    if let Event::Construct(top, n, p) = e {
                        Some((*top, *n, *p))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>();
            assert_eq!(constructions.len(), 1, "{events:?}");
            assert!(constructions[0].0);
            assert_eq!(constructions[0].1, m.atoms.len());
            let ptr = constructions[0].2;
            assert!(
                events
                    .iter()
                    .any(|e| matches!(e, Event::Refine(CanonCompareMode::Atom, _, _)))
            );
            assert!(
                events
                    .iter()
                    .any(|e| matches!(e, Event::Refine(CanonCompareMode::SpecialChirality, _, _)))
            );
            if expect_symmetry {
                assert!(
                    events.iter().any(|e| matches!(
                        e,
                        Event::Refine(CanonCompareMode::SpecialSymmetry, _, _)
                    )),
                    "{events:?}"
                );
            }
            assert!(
                events
                    .iter()
                    .filter_map(|e| if let Event::Refine(_, p, n) = e {
                        Some((*p, *n))
                    } else {
                        None
                    })
                    .all(|(p, n)| p == ptr && n == m.atoms.len())
            );
        }
    }
}

#[cfg(test)]
mod recovery_chem32 {
    use super::*;
    // The comparison-count closure observes the actual private insertion owner.
    fn chem32_holder(id: usize, key: usize) -> CanonBondHolder<'static> {
        CanonBondHolder {
            bond_type: BondOrder::Single,
            bond_stereo: 0,
            stype: BondStereo::None,
            controlling_atoms: [Some(id), None, Some(id + 1), None],
            nbr_sym_class: key,
            nbr_idx: id + 10,
            p_symbol: Some("shared-borrowed-symbol"),
            bond_idx: id,
        }
    }
    fn chem32_full_record<'a>(
        b: &CanonBondHolder<'a>,
    ) -> (
        BondOrder,
        u8,
        BondStereo,
        [Option<usize>; 4],
        usize,
        usize,
        Option<&'a str>,
        usize,
    ) {
        (
            b.bond_type,
            b.bond_stereo,
            b.stype,
            b.controlling_atoms,
            b.nbr_sym_class,
            b.nbr_idx,
            b.p_symbol,
            b.bond_idx,
        )
    }

    #[test]
    fn rdkit_2026_03_6_chem32_insertion_stability_and_exact_source_comparison_counts() {
        for keys in [
            vec![],
            vec![2],
            vec![3, 2, 1],
            vec![2, 2, 2],
            vec![0, 1, 2, 3],
            vec![4, 2, 3, 1],
            vec![5, 1, 4, 2, 3, 0],
            vec![2, 2, 3, 2],
            (0..128).collect(),
            (0..1024).rev().collect(),
        ] {
            let mut holders = keys
                .iter()
                .enumerate()
                .map(|(id, &key)| chem32_holder(id, key))
                .collect::<Vec<_>>();
            let original = holders.clone();
            let mut expected = holders.clone();
            expected.sort_by(|l, r| r.nbr_sym_class.cmp(&l.nbr_sym_class));
            let mut comparisons = 0;
            insertion_sort_canon_bonds_for_update(&mut holders, |l, r| {
                comparisons += 1;
                compare_canon_bond_holder(l, r, &[]) == Ordering::Greater
            });
            assert_eq!(
                holders.iter().map(chem32_full_record).collect::<Vec<_>>(),
                expected.iter().map(chem32_full_record).collect::<Vec<_>>()
            );
            // Expected ordering uses only source comparator keys; complete record IDs
            // remain payload and prove ties are not broken by neighbor/bond IDs.
            for holder in &holders {
                assert_eq!(
                    chem32_full_record(holder),
                    chem32_full_record(&original[holder.bond_idx])
                );
            }
            if keys.windows(2).all(|w| w[0] >= w[1]) {
                assert_eq!(comparisons, keys.len().saturating_sub(1));
            }
            if keys.windows(2).all(|w| w[0] < w[1]) {
                assert_eq!(comparisons, keys.len() * keys.len().saturating_sub(1) / 2);
            }
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem32_insertion_uses_complete_source_greater_precedence() {
        let ranks = [0, 3, 1];
        let mut cases = Vec::new();
        let mut a = chem32_holder(0, 0);
        let mut b = chem32_holder(1, 0);
        a.p_symbol = Some("A");
        b.p_symbol = Some("B");
        a.bond_type = BondOrder::Triple;
        cases.push((a, b));
        let mut a = chem32_holder(0, 0);
        let b = chem32_holder(1, 0);
        a.p_symbol = None;
        a.bond_type = BondOrder::Triple;
        cases.push((b, a));
        let a = chem32_holder(0, 0);
        let mut b = chem32_holder(1, 0);
        b.bond_type = BondOrder::Double;
        cases.push((a, b));
        for (low, high) in [
            (BondStereo::None, BondStereo::Any),
            (BondStereo::Z, BondStereo::E),
            (BondStereo::E, BondStereo::Cis),
        ] {
            let mut a = chem32_holder(0, 0);
            let mut b = chem32_holder(1, 0);
            a.stype = low;
            a.bond_stereo = rdkit_bond_stereo_rank(low);
            b.stype = high;
            b.bond_stereo = rdkit_bond_stereo_rank(high);
            cases.push((a, b));
        }
        let a = chem32_holder(0, 1);
        let b = chem32_holder(1, 2);
        cases.push((a, b));
        let mut a = chem32_holder(0, 4);
        let mut b = chem32_holder(1, 4);
        a.stype = BondStereo::Cis;
        a.bond_stereo = rdkit_bond_stereo_rank(a.stype);
        b.stype = a.stype;
        b.bond_stereo = a.bond_stereo;
        a.controlling_atoms = [Some(0), None, Some(2), None];
        b.controlling_atoms = [Some(0), Some(1), Some(2), None];
        cases.push((a, b));
        for (low, high) in cases {
            assert_eq!(
                compare_canon_bond_holder(&low, &high, &ranks),
                Ordering::Less
            );
            let mut values = [low, high];
            insertion_sort_canon_bonds_for_update(&mut values, |l, r| {
                compare_canon_bond_holder(l, r, &ranks) == Ordering::Greater
            });
            assert_eq!(values[0].bond_idx, high.bond_idx);
            assert_eq!(values[1].bond_idx, low.bond_idx);
        }
    }

    #[test]
    fn rdkit_2026_03_6_chem32_actual_update_refreshes_all_classes_before_insertion() {
        let source_atoms = (0..5)
            .map(|i| {
                Atom::from_spec(
                    AtomId::new(i),
                    cosmolkit_model::AtomSpec::new(cosmolkit_types::Element::C),
                )
            })
            .collect();
        let source_bonds = [(0, 1), (1, 2), (1, 3), (1, 4)]
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    cosmolkit_model::BondSpec::new(
                        AtomId::new(a),
                        AtomId::new(b),
                        BondOrder::Single,
                    ),
                )
            })
            .collect();
        let graph =
            TopologyBlock::try_from_parts(source_atoms, source_bonds, vec![], vec![]).unwrap();
        let view = CanonRankReadView::from_topology(&graph).unwrap();
        let mut atoms =
            init_fragment_canon_atoms(&view, &[true; 5], &[true; 4], true, None, None).unwrap();
        for (atom, index) in atoms.iter_mut().zip([8, 0, 2, 7, 4]) {
            atom.index = index;
        }
        for (i, bond) in atoms[1].bonds.iter_mut().enumerate() {
            bond.nbr_sym_class = 100 + i;
        }
        let original = atoms[1].bonds.clone();
        update_atom_neighbor_index_for_kekulize(&mut atoms, 1);
        assert_eq!(
            atoms[1]
                .bonds
                .iter()
                .map(|b| b.nbr_sym_class)
                .collect::<Vec<_>>(),
            vec![8, 7, 4, 2]
        );
        for bond in &atoms[1].bonds {
            assert_eq!(bond.nbr_sym_class, atoms[bond.nbr_idx].index as usize);
            let mut expected = *original
                .iter()
                .find(|b| b.bond_idx == bond.bond_idx)
                .unwrap();
            expected.nbr_sym_class = bond.nbr_sym_class;
            assert_eq!(chem32_full_record(bond), chem32_full_record(&expected));
        }
    }
}
