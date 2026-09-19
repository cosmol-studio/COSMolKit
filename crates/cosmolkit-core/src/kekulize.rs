// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    cmp::Ordering,
    collections::{BTreeMap, BTreeSet, VecDeque},
};

use crate::{
    RingFindType, RingFindingError, RingInfo, ValenceAssignment, ValenceError, ValenceModel,
    fast_find_rings_from_parts,
};
use cosmolkit_model::{
    AdjacencyList, Atom, AtomId, Bond, BondId, StereoGroupKind, TemplateAttachmentOrderError,
    TopologyBlock, TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

#[derive(Clone, Debug, Eq, PartialEq, thiserror::Error)]
pub enum CanonicalRankError {
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
    #[error(
        "bond {bond} has inconsistent aromatic state: order {order:?}, aromatic flag {is_aromatic}"
    )]
    AromaticBondStateMismatch {
        bond: BondId,
        order: BondOrder,
        is_aromatic: bool,
    },
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
struct QuestionEnumerator {
    questions: Vec<AtomId>,
    position: u32,
    end: u32,
}

fn atom_is_aromatic_for_kekulize(topology: &TopologyBlock, atom: AtomId) -> bool {
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

fn selected_bond_has_type_query(bond: &Bond) -> Result<bool, KekulizeError> {
    if bond.prop("_MolFileBondQueryComplex").is_some() {
        return Err(KekulizeError::UnsupportedQueryState {
            bond: bond.id(),
            detail: "concrete TopologyBlock cannot represent a recursive or composite bond query",
        });
    }
    Ok(bond.prop("_MolFileBondQuery").is_some())
}

fn checked_total_valence(valence: &ValenceAssignment, atom: AtomId) -> Result<i32, KekulizeError> {
    valence.explicit_valence[atom.index()]
        .checked_add(valence.implicit_hydrogens[atom.index()])
        .ok_or(KekulizeError::IntegerOverflow {
            atom,
            field: "total valence",
        })
}

fn prepare_kekulize_selection(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
) -> Result<PreparedKekulizeSelection, KekulizeError> {
    topology.validate()?;
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment selection
    // RDKit✔️✔️: PRECONDITION(atomsToUse.size() == mol.getNumAtoms(),
    // RDKit✔️✔️:              "atomsToUse is wrong size");
    // RDKit✔️✔️: PRECONDITION(bondsToUse.size() == mol.getNumBonds(),
    // RDKit✔️✔️:              "bondsToUse is wrong size");
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
    // RDKit✔️✔️: // if there are no atoms to use we can directly return
    // RDKit✔️✔️: if (atomsToUse.none()) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    if !atoms_in_play.iter().any(|selected| *selected) {
        return Ok(PreparedKekulizeSelection {
            atoms_in_play: atoms_in_play.to_vec(),
            bonds_in_play: bonds_in_play.to_vec(),
            original_total_valences: vec![0; topology.atoms.len()],
            dummy_atoms: vec![false; topology.atoms.len()],
            rings: RingInfo::new(
                RingFindType::Fast,
                topology.atoms.len(),
                topology.bonds.len(),
            ),
            candidate_atom_rings: Vec::new(),
            candidate_bond_rings: Vec::new(),
            found_aromatic: false,
            valence: ValenceAssignment {
                explicit_valence: vec![0; topology.atoms.len()],
                implicit_hydrogens: vec![0; topology.atoms.len()],
            },
        });
    }

    let mut selected_bonds = bonds_in_play.to_vec();
    let mut found_aromatic = false;
    // RDKit✔️✔️: bool foundAromatic = false;
    // RDKit✔️✔️: for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:   if (bondsToUse[bond->getIdx()]) {
    // RDKit✔️✔️:     if (QueryOps::hasBondTypeQuery(*bond)) {
    // RDKit✔️✔️:       bondsToUse[bond->getIdx()] = 0;
    // RDKit✔️✔️:     } else if (bond->getIsAromatic()) {
    // RDKit✔️✔️:       foundAromatic = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for bond in &topology.bonds {
        if !selected_bonds[bond.id().index()] {
            continue;
        }
        if selected_bond_has_type_query(bond)? {
            selected_bonds[bond.id().index()] = false;
            continue;
        }
        if bond.order() == BondOrder::Aromatic && !bond.is_aromatic()
            || bond.is_aromatic()
                && !matches!(
                    bond.order(),
                    BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                )
        {
            return Err(KekulizeError::AromaticBondStateMismatch {
                bond: bond.id(),
                order: bond.order(),
                is_aromatic: bond.is_aromatic(),
            });
        }
        if bond.is_aromatic() {
            found_aromatic = true;
        }
    }

    let valence = crate::assign_valence_with_options_from_parts(
        &topology.atoms,
        &topology.bonds,
        &topology.adjacency,
        ValenceModel::RdkitLike,
        false,
    )?;
    let mut original_total_valences = vec![0; topology.atoms.len()];
    let mut dummy_atoms = vec![false; topology.atoms.len()];
    // RDKit✔️✔️: auto numAtoms = mol.getNumAtoms();
    // RDKit✔️✔️: INT_VECT valences(numAtoms);
    // RDKit✔️✔️: boost::dynamic_bitset<> dummyAts(numAtoms);
    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (!atomsToUse[atom->getIdx()]) {
    // RDKit✔️✔️:     continue;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   atom->calcImplicitValence(false);
    // RDKit✔️✔️:   valences[atom->getIdx()] = atom->getTotalValence();
    // RDKit✔️✔️:   if (isAromaticAtom(*atom)) {
    // RDKit✔️✔️:     foundAromatic = true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atom->getAtomicNum()) {
    // RDKit✔️✔️:     dummyAts[atom->getIdx()] = 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for (atom_idx, selected) in atoms_in_play.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom_id = AtomId::new(atom_idx);
        original_total_valences[atom_idx] = checked_total_valence(&valence, atom_id)?;
        found_aromatic |= atom_is_aromatic_for_kekulize(topology, atom_id);
        dummy_atoms[atom_idx] = topology.atoms[atom_idx].atomic_number() == 0;
    }

    // RDKit✔️✔️: if (!foundAromatic) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    let rings = if found_aromatic {
        fast_find_rings_from_parts(topology.atoms.len(), &topology.bonds, &topology.adjacency)?
    } else {
        RingInfo::new(
            RingFindType::Fast,
            topology.atoms.len(),
            topology.bonds.len(),
        )
    };

    let mut candidate_rings = VecDeque::new();
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
    if found_aromatic {
        for (atom_ring, bond_ring) in rings.atom_rings().iter().zip(rings.bond_rings()) {
            let atoms_are_selected = atom_ring.iter().all(|atom| atoms_in_play[atom.index()]);
            let contains_non_dummy = atom_ring.iter().any(|atom| !dummy_atoms[atom.index()]);
            let bonds_are_selected = bond_ring.iter().all(|bond| selected_bonds[bond.index()]);
            if atoms_are_selected && contains_non_dummy && bonds_are_selected {
                let wedge_start = atom_ring.iter().position(|atom| wedged_atoms[atom.index()]);
                let start = wedge_start.unwrap_or(0);
                let mut rotated_atoms = atom_ring.clone();
                let mut rotated_bonds = bond_ring.clone();
                rotated_atoms.rotate_left(start);
                rotated_bonds.rotate_left(start);
                if wedge_start.is_some() {
                    candidate_rings.push_front((rotated_atoms, rotated_bonds));
                } else {
                    candidate_rings.push_back((rotated_atoms, rotated_bonds));
                }
            }
        }
    }
    let (candidate_atom_rings, candidate_bond_rings) = candidate_rings.into_iter().unzip();
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment selection
    Ok(PreparedKekulizeSelection {
        atoms_in_play: atoms_in_play.to_vec(),
        bonds_in_play: selected_bonds,
        original_total_valences,
        dummy_atoms,
        rings,
        candidate_atom_rings,
        candidate_bond_rings,
        found_aromatic,
        valence,
    })
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

    let mut working = topology.clone();
    let mut double_bond_candidates = vec![false; topology.atoms.len()];
    let mut questions = Vec::new();
    let mut done = Vec::new();
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: markDbondCands
    // RDKit✔️✔️: bool hasAromaticOrDummyAtom =
    // RDKit✔️✔️:     std::any_of(allAtms.begin(), allAtms.end(), [&mol](int allAtm) {
    // RDKit✔️✔️:       return (!mol.getAtomWithIdx(allAtm)->getAtomicNum() ||
    // RDKit✔️✔️:               isAromaticAtom(*mol.getAtomWithIdx(allAtm)));
    // RDKit✔️✔️:     });
    let has_aromatic_or_dummy_atom = all_atoms.iter().any(|&atom| {
        topology.atoms[atom.index()].atomic_number() == 0
            || atom_is_aromatic_for_kekulize(topology, atom)
    });
    // RDKit✔️✔️: if (!hasAromaticOrDummyAtom) {
    // RDKit✔️✔️:   return;
    // RDKit✔️✔️: }
    if !has_aromatic_or_dummy_atom {
        return Ok(CandidateState {
            topology: working,
            double_bond_candidates,
            questions,
            done,
        });
    }

    // RDKit✔️✔️: boost::dynamic_bitset<> isRingNotCand(mol.getRingInfo()->numRings());
    // RDKit✔️✔️: unsigned int ri = 0;
    // RDKit✔️✔️: for (const auto &aring : mol.getRingInfo()->atomRings()) {
    // RDKit✔️✔️:   isRingNotCand.set(ri);
    // RDKit✔️✔️:   for (auto ai : aring) {
    // RDKit✔️✔️:     const auto at = mol.getAtomWithIdx(ai);
    // RDKit✔️✔️:     if (isAromaticAtom(*at) && mol.getRingInfo()->numAtomRings(ai) == 1) {
    // RDKit✔️✔️:       isRingNotCand.reset(ri);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   ++ri;
    // RDKit✔️✔️: }
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

    // RDKit✔️✔️: boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
    // RDKit✔️✔️: for (int allAtm : allAtms) {
    for &atom_id in all_atoms {
        // RDKit✔️✔️:   inAllAtms.set(allAtm);
        // RDKit✔️✔️:   Atom *at = mol.getAtomWithIdx(allAtm);
        in_all_atoms[atom_id.index()] = true;
        let atom = &topology.atoms[atom_id.index()];
        // RDKit✔️✔️:   if (at->getAtomicNum() && !isAromaticAtom(*at)) {
        // RDKit✔️✔️:     done.push_back(allAtm);
        // RDKit✔️✔️:     continue;
        // RDKit✔️✔️:   }
        if atom.atomic_number() != 0 && !atom_is_aromatic_for_kekulize(topology, atom_id) {
            done.push(atom_id);
            continue;
        }

        // RDKit✔️✔️:   int sbo = 0;
        // RDKit✔️✔️:   unsigned nToIgnore = 0;
        // RDKit✔️✔️:   unsigned int nonArNonDummyNbr = 0;
        let mut single_bond_order = 0i32;
        let mut neighbors_to_ignore = 0usize;
        let mut non_aromatic_non_dummy_neighbors = 0usize;
        // RDKit✔️✔️:   for (const auto bond : mol.atomBonds(at)) {
        for neighbor in topology.adjacency.neighbors_of(atom_id.index()) {
            let bond = &topology.bonds[neighbor.bond.index()];
            let other = &topology.atoms[neighbor.atom_index];
            // RDKit✔️✔️:     auto otherAt = bond->getOtherAtom(at);
            // RDKit✔️✔️:     if (otherAt->getAtomicNum() && !otherAt->getIsAromatic() &&
            // RDKit✔️✔️:         inAllAtms.test(otherAt->getIdx())) {
            // RDKit✔️✔️:       ++nonArNonDummyNbr;
            // RDKit✔️✔️:     }
            if other.atomic_number() != 0
                && !other.is_aromatic()
                && in_all_atoms[neighbor.atom_index]
            {
                non_aromatic_non_dummy_neighbors += 1;
            }
            // RDKit✔️✔️:     if (bond->getIsAromatic() && (bond->getBondType() == Bond::SINGLE ||
            // RDKit✔️✔️:                                   bond->getBondType() == Bond::DOUBLE ||
            // RDKit✔️✔️:                                   bond->getBondType() == Bond::AROMATIC)) {
            // RDKit✔️✔️:       ++sbo;
            // RDKit✔️✔️:       makeSingle.push_back(bond);
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
                // RDKit✔️✔️:     } else {
                // RDKit✔️✔️:       int bondContrib = std::lround(bond->getValenceContrib(at));
                // RDKit✔️✔️:       sbo += bondContrib;
                // RDKit✔️✔️:       if (!bondContrib) {
                // RDKit✔️✔️:         ++nToIgnore;
                // RDKit✔️✔️:       }
                // RDKit✔️✔️:     }
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

        // RDKit✔️✔️:   auto numAtomRings = mol.getRingInfo()->numAtomRings(at->getIdx());
        // RDKit✔️✔️:   const auto &riVect = mol.getRingInfo()->atomMembers(at->getIdx());
        // RDKit✔️✔️:   size_t numNonCandRings = std::count_if(
        // RDKit✔️✔️:       riVect.begin(), riVect.end(),
        // RDKit✔️✔️:       [&isRingNotCand](int ri) { return isRingNotCand.test(ri); });
        let atom_ring_count = rings.num_atom_rings(atom_id);
        let non_candidate_ring_count = rings
            .atom_members(atom_id)
            .iter()
            .filter(|&&ring| is_ring_not_candidate[ring])
            .count();
        // RDKit✔️✔️:   if (!at->getAtomicNum() && nonArNonDummyNbr < numAtomRings &&
        // RDKit✔️✔️:       numNonCandRings < numAtomRings) {
        if atom.atomic_number() == 0
            && non_aromatic_non_dummy_neighbors < atom_ring_count
            && non_candidate_ring_count < atom_ring_count
        {
            // RDKit✔️✔️:     dBndCands[allAtm] = 1;
            // RDKit✔️✔️:     questions.push_back(allAtm);
            double_bond_candidates[atom_id.index()] = true;
            questions.push(atom_id);
        } else {
            // RDKit✔️✔️:     sbo += at->getTotalNumHs();
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
            // RDKit✔️✔️:     auto dv =
            // RDKit✔️✔️:         PeriodicTable::getTable()->getDefaultValence(at->getAtomicNum());
            // RDKit✔️✔️:     auto chrg = at->getFormalCharge();
            // RDKit✔️✔️:     if (isEarlyAtom(at->getAtomicNum())) {
            // RDKit✔️✔️:       chrg = -chrg;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     if (at->getAtomicNum() == 6 && chrg > 0) {
            // RDKit✔️✔️:       chrg = -chrg;
            // RDKit✔️✔️:     }
            // RDKit✔️✔️:     dv += chrg;
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
            // RDKit✔️✔️:     int tbo = at->getTotalValence();
            // RDKit✔️✔️:     int nRadicals = at->getNumRadicalElectrons();
            // RDKit✔️✔️:     int totalDegree = at->getDegree() +
            // RDKit✔️✔️:                       at->getValence(Atom::ValenceType::IMPLICIT) - nToIgnore;
            let total_bond_order = checked_total_valence(valence, atom_id)?;
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
            // RDKit✔️✔️:     const auto &valList =
            // RDKit✔️✔️:         PeriodicTable::getTable()->getValenceList(at->getAtomicNum());
            // RDKit✔️✔️:     unsigned int vi = 1;
            // RDKit✔️✔️:     while (tbo > dv && vi < valList.size() && valList[vi] > 0) {
            // RDKit✔️✔️:       dv = valList[vi] + chrg;
            // RDKit✔️✔️:       ++vi;
            // RDKit✔️✔️:     }
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
            // RDKit✔️✔️:     if (tbo == 5 && sbo == 4 && dv == 3 && totalDegree == 3 &&
            // RDKit✔️✔️:         nRadicals == 0 && chrg == 0 && at->getTotalNumHs() == 0) {
            // RDKit✔️✔️:       switch (at->getAtomicNum()) {
            // RDKit✔️✔️:         case 7:
            // RDKit✔️✔️:         case 15:
            // RDKit✔️✔️:         case 33:
            // RDKit✔️✔️:           dv = 5;
            // RDKit✔️✔️:           break;
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
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
            // RDKit✔️✔️:     if (totalDegree + nRadicals >= dv) {
            // RDKit✔️✔️:       continue;
            // RDKit✔️✔️:     }
            if total_degree + radical_electrons >= default_valence {
                continue;
            }
            // RDKit✔️✔️:     if (dv == (sbo + 1 + nRadicals)) {
            // RDKit✔️✔️:       dBndCands[allAtm] = 1;
            // RDKit✔️✔️:     } else if (!nRadicals && at->getNoImplicit() && dv == (sbo + 2)) {
            // RDKit✔️✔️:       dBndCands[allAtm] = 1;
            // RDKit✔️✔️:     }
            if default_valence == single_bond_order + 1 + radical_electrons
                || (radical_electrons == 0
                    && atom.no_implicit()
                    && default_valence == single_bond_order + 2)
            {
                double_bond_candidates[atom_id.index()] = true;
            }
        }
    }
    // RDKit✔️✔️: }  // loop over all atoms in the fused system
    // RDKit✔️✔️: for (auto &bi : makeSingle) {
    // RDKit✔️✔️:   bi->setBondType(Bond::SINGLE);
    // RDKit✔️✔️: }
    for (bond_idx, should_make_single) in make_single.into_iter().enumerate() {
        if should_make_single {
            working.bonds[bond_idx].set_order(BondOrder::Single);
        }
    }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: markDbondCands
    Ok(CandidateState {
        topology: working,
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

    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeWorker
    // RDKit✔️✔️: bool kekulizeWorker(RWMol &mol, const INT_VECT &allAtms,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> dBndCands,
    // RDKit✔️✔️:                     boost::dynamic_bitset<> dBndAdds, INT_VECT done,
    // RDKit✔️✔️:                     const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit✔️✔️:   INT_DEQUE astack;
    // RDKit✔️✔️:   INT_INT_DEQ_MAP options;
    // RDKit✔️✔️:   int lastOpt = -1;
    // RDKit✔️✔️:   boost::dynamic_bitset<> localBondsAdded(mol.getNumBonds());
    let mut atom_stack = VecDeque::new();
    let mut options = BTreeMap::<AtomId, VecDeque<AtomId>>::new();
    let mut last_option = None;
    let mut local_bonds_added = vec![false; topology.bonds.len()];
    let mut double_bond_candidates = initial_candidates.to_vec();
    let mut double_bonds_added = initial_bonds_added.to_vec();
    let mut done = initial_done.to_vec();

    // RDKit✔️✔️:   boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
    // RDKit✔️✔️:   for (int allAtm : allAtms) {
    // RDKit✔️✔️:     inAllAtms.set(allAtm);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   auto lessByRank = [&atomRanks](int a, int b) {
    // RDKit✔️✔️:     const auto ra = atomRanks.at(static_cast<unsigned int>(a));
    // RDKit✔️✔️:     const auto rb = atomRanks.at(static_cast<unsigned int>(b));
    // RDKit✔️✔️:     return (ra < rb) || (ra == rb && a < b);
    // RDKit✔️✔️:   };
    // `in_all_atoms` was populated during checked input validation above.

    // RDKit✔️✔️:   boost::dynamic_bitset<> wedgeEndAtoms(mol.getNumAtoms());
    // RDKit✔️✔️:   for (const auto bond : mol.bonds()) {
    // RDKit✔️✔️:     if (bond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit✔️✔️:         bond->getBondDir() == Bond::BondDir::BEGINDASH) {
    // RDKit✔️✔️:       const auto endIdx = bond->getEndAtomIdx();
    // RDKit✔️✔️:       if (inAllAtms.test(endIdx)) {
    // RDKit✔️✔️:         wedgeEndAtoms.set(endIdx);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
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

    // RDKit✔️✔️:   INT_VECT sortedAtms(allAtms);
    // RDKit✔️✔️:   std::sort(sortedAtms.begin(), sortedAtms.end(),
    // RDKit✔️✔️:             [&wedgeEndAtoms, &lessByRank](int a, int b) {
    // RDKit✔️✔️:               const bool wa = wedgeEndAtoms.test(a);
    // RDKit✔️✔️:               const bool wb = wedgeEndAtoms.test(b);
    // RDKit✔️✔️:               if (wa != wb) {
    // RDKit✔️✔️:                 return wa;
    // RDKit✔️✔️:               }
    // RDKit✔️✔️:               return lessByRank(a, b);
    // RDKit✔️✔️:             });
    let mut sorted_atoms = all_atoms.to_vec();
    sorted_atoms.sort_by_key(|atom| {
        (
            !wedge_end_atoms[atom.index()],
            atom_ranks[atom.index()],
            atom.index(),
        )
    });

    // RDKit✔️✔️:   int curr = -1;
    // RDKit✔️✔️:   INT_DEQUE btmoves;
    // RDKit✔️✔️:   unsigned int numBT = 0;
    // RDKit✔️✔️:   while ((done.size() < sortedAtms.size()) || !astack.empty()) {
    let mut backtrack_moves = Vec::new();
    let mut backtracks = 0u32;
    while done.len() < sorted_atoms.len() || !atom_stack.is_empty() {
        // RDKit✔️✔️:     if (astack.size() > 0) {
        // RDKit✔️✔️:       curr = astack.front();
        // RDKit✔️✔️:       astack.pop_front();
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       curr = -1;
        // RDKit✔️✔️:       for (int allAtm : sortedAtms) {
        // RDKit✔️✔️:         if (std::find(done.begin(), done.end(), allAtm) == done.end()) {
        // RDKit✔️✔️:           curr = allAtm;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     CHECK_INVARIANT(curr >= 0, "starting point not found");
        // RDKit✔️✔️:     done.push_back(curr);
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

        // RDKit✔️✔️:     INT_DEQUE opts;
        // RDKit✔️✔️:     bool cCand = false;
        // RDKit✔️✔️:     if (dBndCands[curr]) {
        // RDKit✔️✔️:       cCand = true;
        // RDKit✔️✔️:     }
        let current_is_candidate = double_bond_candidates[current.index()];
        // RDKit✔️✔️:     if (options.find(curr) != options.end()) {
        // RDKit✔️✔️:       opts = options[curr];
        // RDKit✔️✔️:       CHECK_INVARIANT(opts.size() > 0, "");
        // RDKit✔️✔️:     } else {
        let mut current_options = if let Some(stored) = options.get(&current) {
            stored.clone()
        } else {
            // RDKit✔️✔️:       INT_DEQUE lstack;
            // RDKit✔️✔️:       std::vector<int> optsV;
            // RDKit✔️✔️:       std::vector<int> wedgedOptsV;
            // RDKit✔️✔️:       std::vector<int> nbrs;
            let mut local_stack = VecDeque::new();
            let mut ordinary_options = VecDeque::new();
            let mut wedged_options = VecDeque::new();
            // RDKit✔️✔️:       for (auto nbrAtom : mol.atomNeighbors(mol.getAtomWithIdx(curr))) {
            // RDKit✔️✔️:         const auto nbrIdx = static_cast<int>(nbrAtom->getIdx());
            // RDKit✔️✔️:         if (!inAllAtms.test(nbrIdx)) {
            // RDKit✔️✔️:           continue;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         if (std::find(done.begin(), done.end(), nbrIdx) != done.end()) {
            // RDKit✔️✔️:           continue;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:         nbrs.push_back(nbrIdx);
            // RDKit✔️✔️:       }
            let mut neighbors = topology
                .adjacency
                .neighbors_of(current.index())
                .iter()
                .map(|neighbor| AtomId::new(neighbor.atom_index))
                .filter(|neighbor| in_all_atoms[neighbor.index()] && !done.contains(neighbor))
                .collect::<Vec<_>>();
            // RDKit✔️✔️:       std::sort(nbrs.begin(), nbrs.end(), lessByRank);
            neighbors.sort_by_key(|atom| (atom_ranks[atom.index()], atom.index()));

            // RDKit✔️✔️:       for (int nbrIdx : nbrs) {
            // RDKit✔️✔️:         auto nbrBond = mol.getBondBetweenAtoms(curr, nbrIdx);
            for neighbor in neighbors {
                let bond_id = bond_between_atoms(&topology, current, neighbor)?;
                let bond = &topology.bonds[bond_id.index()];
                // RDKit✔️✔️:         if (std::find(astack.begin(), astack.end(), nbrIdx) == astack.end()) {
                // RDKit✔️✔️:           lstack.push_back(nbrIdx);
                // RDKit✔️✔️:         }
                if !atom_stack.contains(&neighbor) {
                    local_stack.push_back(neighbor);
                }
                // RDKit✔️✔️:         if (cCand && dBndCands[nbrIdx] &&
                // RDKit✔️✔️:             (nbrBond->getIsAromatic() ||
                // RDKit✔️✔️:              mol.getAtomWithIdx(curr)->getAtomicNum() == 0 ||
                // RDKit✔️✔️:              mol.getAtomWithIdx(nbrIdx)->getAtomicNum() == 0)) {
                if current_is_candidate
                    && double_bond_candidates[neighbor.index()]
                    && (bond.is_aromatic()
                        || topology.atoms[current.index()].atomic_number() == 0
                        || topology.atoms[neighbor.index()].atomic_number() == 0)
                {
                    // RDKit✔️✔️:           if (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
                    // RDKit✔️✔️:               nbrBond->getBondDir() == Bond::BondDir::BEGINDASH) {
                    // RDKit✔️✔️:             wedgedOptsV.push_back(nbrIdx);
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             optsV.push_back(nbrIdx);
                    // RDKit✔️✔️:           }
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
            // RDKit✔️✔️:       for (int v : optsV) {
            // RDKit✔️✔️:         opts.push_back(v);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       for (int v : wedgedOptsV) {
            // RDKit✔️✔️:         opts.push_back(v);
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:       astack.insert(astack.end(), lstack.begin(), lstack.end());
            ordinary_options.append(&mut wedged_options);
            atom_stack.append(&mut local_stack);
            ordinary_options
        };

        // RDKit✔️✔️:     if (cCand) {
        if current_is_candidate {
            // RDKit✔️✔️:       if (!opts.empty()) {
            if let Some(neighbor) = current_options.pop_front() {
                // RDKit✔️✔️:         ncnd = opts.front();
                // RDKit✔️✔️:         opts.pop_front();
                // RDKit✔️✔️:         auto bnd = mol.getBondBetweenAtoms(curr, ncnd);
                // RDKit✔️✔️:         bnd->setBondType(Bond::DOUBLE);
                // RDKit✔️✔️:         if (bnd->getBondDir() != Bond::BondDir::NONE) {
                // RDKit✔️✔️:           bnd->setBondDir(Bond::BondDir::NONE);
                // RDKit✔️✔️:         }
                let bond_id = bond_between_atoms(&topology, current, neighbor)?;
                topology.bonds[bond_id.index()].set_order(BondOrder::Double);
                if topology.bonds[bond_id.index()].direction() != BondDirection::None {
                    topology.bonds[bond_id.index()].set_direction(BondDirection::None);
                }
                // RDKit✔️✔️:         dBndCands[curr] = 0;
                // RDKit✔️✔️:         dBndCands[ncnd] = 0;
                // RDKit✔️✔️:         dBndAdds[bnd->getIdx()] = 1;
                // RDKit✔️✔️:         localBondsAdded[bnd->getIdx()] = 1;
                double_bond_candidates[current.index()] = false;
                double_bond_candidates[neighbor.index()] = false;
                double_bonds_added[bond_id.index()] = true;
                local_bonds_added[bond_id.index()] = true;

                // RDKit✔️✔️:         if (options.find(curr) != options.end()) {
                if options.contains_key(&current) {
                    // RDKit✔️✔️:           if (opts.size() == 0) {
                    // RDKit✔️✔️:             options.erase(curr);
                    // RDKit✔️✔️:             btmoves.pop_back();
                    // RDKit✔️✔️:             if (btmoves.size() > 0) {
                    // RDKit✔️✔️:               lastOpt = btmoves.back();
                    // RDKit✔️✔️:             } else {
                    // RDKit✔️✔️:               lastOpt = -1;
                    // RDKit✔️✔️:             }
                    // RDKit✔️✔️:           } else {
                    // RDKit✔️✔️:             options[curr] = opts;
                    // RDKit✔️✔️:           }
                    if current_options.is_empty() {
                        options.remove(&current);
                        backtrack_moves.pop();
                        last_option = backtrack_moves.last().copied();
                    } else {
                        options.insert(current, current_options);
                    }
                } else {
                    // RDKit✔️✔️:         } else {
                    // RDKit✔️✔️:           if (opts.size() > 0) {
                    // RDKit✔️✔️:             lastOpt = curr;
                    // RDKit✔️✔️:             btmoves.push_back(lastOpt);
                    // RDKit✔️✔️:             options[curr] = opts;
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:         }
                    if !current_options.is_empty() {
                        last_option = Some(current);
                        backtrack_moves.push(current);
                        options.insert(current, current_options);
                    }
                }
            } else if topology.atoms[current.index()].atomic_number() != 0 {
                // RDKit✔️✔️:       } else if (mol.getAtomWithIdx(curr)->getAtomicNum()) {
                // RDKit✔️✔️:         if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
                if let Some(anchor) = last_option.filter(|_| backtracks < max_backtracks) {
                    // RDKit✔️✔️:           backTrack(mol, options, lastOpt, done, astack,
                    // RDKit✔️✔️:                     dBndCands, dBndAdds);
                    // RDKit✔️✔️:           ++numBT;
                    backtrack_kekulize(
                        &mut topology,
                        anchor,
                        &mut done,
                        &mut atom_stack,
                        &mut double_bond_candidates,
                        &mut double_bonds_added,
                    )?;
                    backtracks += 1;
                } else {
                    // RDKit✔️✔️:         } else {
                    // RDKit✔️✔️:           for (unsigned int bidx = 0; bidx < mol.getNumBonds(); ++bidx) {
                    // RDKit✔️✔️:             if (localBondsAdded[bidx]) {
                    // RDKit✔️✔️:               mol.getBondWithIdx(bidx)->setBondType(Bond::SINGLE);
                    // RDKit✔️✔️:             }
                    // RDKit✔️✔️:           }
                    // RDKit✔️✔️:           return false;
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
                    // RDKit✔️✔️:         }
                    // RDKit✔️✔️:       }
                    // RDKit✔️✔️:     }
                    // RDKit✔️✔️:   }
                    // RDKit✔️✔️:   return true;
                    // RDKit✔️✔️: }
                    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeWorker
                    return Ok(MatchingState {
                        topology,
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
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeWorker
    Ok(MatchingState {
        topology,
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
        // RDKit✔️✔️: class QuestionEnumerator {
        // RDKit✔️✔️:  public:
        // RDKit✔️✔️:   QuestionEnumerator(INT_VECT questions)
        // RDKit✔️✔️:       : d_questions(std::move(questions)), d_pos(1) {};
        // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::QuestionEnumerator
        let question_count = questions.len();
        if question_count >= u32::BITS as usize {
            return Err(KekulizeError::QuestionSubsetOverflow {
                questions: question_count,
                bit_width: u32::BITS,
            });
        }
        let end = 1u32 << question_count;
        Ok(Self {
            questions,
            position: 1,
            end,
        })
    }

    fn next(&mut self) -> Vec<AtomId> {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::next
        // RDKit✔️✔️:   INT_VECT next() {
        // RDKit✔️✔️:     INT_VECT res;
        // RDKit✔️✔️:     if (d_pos >= (0x1u << d_questions.size())) {
        // RDKit✔️✔️:       return res;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     for (unsigned int i = 0; i < d_questions.size(); ++i) {
        // RDKit✔️✔️:       if (d_pos & (0x1u << i)) {
        // RDKit✔️✔️:         res.push_back(d_questions[i]);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++d_pos;
        // RDKit✔️✔️:     return res;
        // RDKit✔️✔️:   };
        // RDKit✔️✔️:
        // RDKit✔️✔️:  private:
        // RDKit✔️✔️:   INT_VECT d_questions;
        // RDKit✔️✔️:   unsigned int d_pos;
        // RDKit✔️✔️: };
        // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: QuestionEnumerator::next
        if self.position >= self.end {
            return Vec::new();
        }
        let mut selected = Vec::new();
        for (index, &question) in self.questions.iter().enumerate() {
            if self.position & (1u32 << index) != 0 {
                selected.push(question);
            }
        }
        self.position += 1;
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
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: permuteDummiesAndKekulize
    // RDKit✔️✔️: bool permuteDummiesAndKekulize(RWMol &mol, const INT_VECT &allAtms,
    // RDKit✔️✔️:                                boost::dynamic_bitset<> dBndCands,
    // RDKit✔️✔️:                                INT_VECT &questions,
    // RDKit✔️✔️:                                const UINT_VECT &atomRanks,
    // RDKit✔️✔️:                                unsigned int maxBackTracks) {
    // RDKit✔️✔️:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    // RDKit✔️✔️:   for (int allAtm : allAtms) {
    // RDKit✔️✔️:     atomsInPlay[allAtm] = 1;
    // RDKit✔️✔️:   }
    let mut atoms_in_play = vec![false; topology.atoms.len()];
    for &atom in all_atoms {
        atoms_in_play[atom.index()] = true;
    }
    // RDKit✔️✔️:   bool kekulized = false;
    // RDKit✔️✔️:   QuestionEnumerator qEnum(questions);
    let mut question_enumerator = QuestionEnumerator::new(questions.to_vec())?;
    // RDKit✔️✔️:   while (!kekulized && questions.size()) {
    while !questions.is_empty() {
        // RDKit✔️✔️:     boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
        // RDKit✔️✔️:     INT_VECT done;
        let double_bonds_added = vec![false; topology.bonds.len()];
        // RDKit✔️✔️:     // reset the state: all aromatic bonds are remarked to single:
        // RDKit✔️✔️:     for (const auto bond : mol.bonds()) {
        // RDKit✔️✔️:       if (bond->getIsAromatic() && bond->getBondType() != Bond::SINGLE &&
        // RDKit✔️✔️:           atomsInPlay[bond->getBeginAtomIdx()] &&
        // RDKit✔️✔️:           atomsInPlay[bond->getEndAtomIdx()]) {
        // RDKit✔️✔️:         bond->setBondType(Bond::SINGLE);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for bond in &mut topology.bonds {
            if bond.is_aromatic()
                && bond.order() != BondOrder::Single
                && atoms_in_play[bond.begin().index()]
                && atoms_in_play[bond.end().index()]
            {
                bond.set_order(BondOrder::Single);
            }
        }
        // RDKit✔️✔️:     // pick a new permutation of the questionable atoms:
        // RDKit✔️✔️:     const auto &switchOff = qEnum.next();
        // RDKit✔️✔️:     if (!switchOff.size()) {
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        let switch_off = question_enumerator.next();
        if switch_off.is_empty() {
            break;
        }
        // RDKit✔️✔️:     auto tCands = dBndCands;
        // RDKit✔️✔️:     for (int it : switchOff) {
        // RDKit✔️✔️:       tCands[it] = 0;
        // RDKit✔️✔️:     }
        let mut trial_candidates = initial_candidates.to_vec();
        for atom in switch_off {
            trial_candidates[atom.index()] = false;
        }
        // RDKit✔️✔️:     // try kekulizing again:
        // RDKit✔️✔️:     kekulized =
        // RDKit✔️✔️:         kekulizeWorker(mol, allAtms, tCands, dBndAdds, done, atomRanks,
        // RDKit✔️✔️:                        maxBackTracks);
        let trial = kekulize_matching_worker(
            topology,
            all_atoms,
            &trial_candidates,
            &double_bonds_added,
            &[],
            atom_ranks,
            max_backtracks,
        )?;
        topology = trial.topology;
        if trial.succeeded {
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   return kekulized;
            // RDKit✔️✔️: }
            // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: permuteDummiesAndKekulize
            return Ok(FusedKekulizeState {
                topology,
                succeeded: true,
                problem_atoms: Vec::new(),
            });
        }
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return kekulized;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: permuteDummiesAndKekulize
    Ok(FusedKekulizeState {
        topology,
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
    topology: TopologyBlock,
    atom_rings: &[Vec<AtomId>],
    all_rings: &RingInfo,
    valence: &ValenceAssignment,
    atom_ranks: &[usize],
    max_backtracks: u32,
) -> Result<FusedKekulizeState, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeFused
    // RDKit✔️✔️: void kekulizeFused(RWMol &mol, const VECT_INT_VECT &arings,
    // RDKit✔️✔️:                    const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit✔️✔️:   // get all the atoms in the ring system
    // RDKit✔️✔️:   INT_VECT allAtms;
    // RDKit✔️✔️:   Union(arings, allAtms);
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/types.cpp :: Union(VECT_INT_VECT)
    // RDKit✔️✔️: void Union(const VECT_INT_VECT &rings, INT_VECT &res, const INT_VECT *exclude) {
    // RDKit✔️✔️:   res.resize(0);
    // RDKit✔️✔️:   INT_VECT ring;
    // RDKit✔️✔️:   unsigned int id;
    // RDKit✔️✔️:   auto nrings = static_cast<unsigned int>(rings.size());
    // RDKit✔️✔️:   INT_VECT_CI ri;
    // RDKit✔️✔️:   for (id = 0; id < nrings; id++) {
    // RDKit✔️✔️:     if (exclude) {
    // RDKit✔️✔️:       if (std::find(exclude->begin(), exclude->end(), static_cast<int>(id)) !=
    // RDKit✔️✔️:           exclude->end()) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ring = rings[id];
    // RDKit✔️✔️:     for (ri = ring.begin(); ri != ring.end(); ri++) {
    // RDKit✔️✔️:       if (std::find(res.begin(), res.end(), (*ri)) == res.end()) {
    // RDKit✔️✔️:         res.push_back(*ri);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/RDGeneral/types.cpp :: Union(VECT_INT_VECT)
    let mut all_atoms = Vec::new();
    for ring in atom_rings {
        for &atom in ring {
            if !all_atoms.contains(&atom) {
                all_atoms.push(atom);
            }
        }
    }
    // RDKit✔️✔️:   // get all the atoms that are candidates to receive a double bond
    // RDKit✔️✔️:   // also mark atoms in the fused system that are not aromatic to begin with
    // RDKit✔️✔️:   // as done. Mark all the bonds that are part of the aromatic system
    // RDKit✔️✔️:   // to be single bonds
    // RDKit✔️✔️:   INT_VECT done;
    // RDKit✔️✔️:   INT_VECT questions;
    // RDKit✔️✔️:   auto nats = mol.getNumAtoms();
    // RDKit✔️✔️:   auto nbnds = mol.getNumBonds();
    // RDKit✔️✔️:   boost::dynamic_bitset<> dBndCands(nats);
    // RDKit✔️✔️:   boost::dynamic_bitset<> dBndAdds(nbnds);
    // RDKit✔️✔️:   markDbondCands(mol, allAtms, dBndCands, questions, done);
    let candidate_state = mark_double_bond_candidates(&topology, &all_atoms, all_rings, valence)?;
    let initial_candidates = candidate_state.double_bond_candidates.clone();
    let questions = candidate_state.questions.clone();
    // RDKit✔️✔️:   auto kekulized =
    // RDKit✔️✔️:       kekulizeWorker(mol, allAtms, dBndCands, dBndAdds, done, atomRanks,
    // RDKit✔️✔️:                      maxBackTracks);
    let first_attempt = kekulize_matching_worker(
        candidate_state.topology,
        &all_atoms,
        &initial_candidates,
        &vec![false; topology.bonds.len()],
        &candidate_state.done,
        atom_ranks,
        max_backtracks,
    )?;
    if first_attempt.succeeded {
        return Ok(FusedKekulizeState {
            topology: first_attempt.topology,
            succeeded: true,
            problem_atoms: Vec::new(),
        });
    }
    // RDKit✔️✔️:   if (!kekulized && questions.size()) {
    // RDKit✔️✔️:     // we failed, but there are some dummy atoms we can try permuting.
    // RDKit✔️✔️:     kekulized = permuteDummiesAndKekulize(mol, allAtms, dBndCands, questions,
    // RDKit✔️✔️:                                           atomRanks, maxBackTracks);
    // RDKit✔️✔️:   }
    if !questions.is_empty() {
        let permuted = permute_dummies_and_kekulize(
            first_attempt.topology,
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
    // RDKit✔️✔️:   if (!kekulized) {
    // RDKit✔️✔️:     // we exhausted all option (or crossed the allowed
    // RDKit✔️✔️:     // number of backTracks) and we still need to backtrack
    // RDKit✔️✔️:     // can't kekulize this thing
    // RDKit✔️✔️:     std::vector<unsigned int> problemAtoms;
    // RDKit✔️✔️:     std::ostringstream errout;
    // RDKit✔️✔️:     errout << "Can't kekulize mol.";
    // RDKit✔️✔️:     errout << "  Unkekulized atoms:";
    // RDKit✔️✔️:     for (unsigned int i = 0; i < nats; ++i) {
    // RDKit✔️✔️:       if (dBndCands[i]) {
    // RDKit✔️✔️:         errout << " " << i;
    // RDKit✔️✔️:         problemAtoms.push_back(i);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::string msg = errout.str();
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << msg << std::endl;
    // RDKit✔️✔️:     throw KekulizeException(msg, problemAtoms);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: kekulizeFused
    Ok(FusedKekulizeState {
        topology: first_attempt.topology,
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
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment fused-system dispatch
    // RDKit✔️✔️:     // make a neighbor map for the rings i.e. a ring is a
    // RDKit✔️✔️:     // neighbor to another candidate ring if it shares at least
    // RDKit✔️✔️:     // one bond
    // RDKit✔️✔️:     // useful to figure out fused systems
    // RDKit✔️✔️:     INT_INT_VECT_MAP neighMap;
    // RDKit✔️✔️:     RingUtils::makeRingNeighborMap(brings, neighMap);
    let neighbor_map = make_ring_neighbor_map(bond_rings);
    // RDKit✔️✔️:     int curr = 0;
    // RDKit✔️✔️:     int cnrs = rdcast<int>(arings.size());
    // RDKit✔️✔️:     boost::dynamic_bitset<> fusDone(cnrs);
    let mut done = vec![false; atom_rings.len()];
    // RDKit✔️✔️:     while (curr < cnrs) {
    for current in 0..atom_rings.len() {
        if done[current] {
            continue;
        }
        // RDKit✔️✔️:       INT_VECT fused;
        // RDKit✔️✔️:       RingUtils::pickFusedRings(curr, neighMap, fused, fusDone);
        let fused = pick_fused_rings(current, &neighbor_map, &mut done);
        // RDKit✔️✔️:       VECT_INT_VECT frings(fused.size());
        // RDKit✔️✔️:       std::transform(fused.begin(), fused.end(), frings.begin(),
        // RDKit✔️✔️:                      [&arings](const int ri) { return arings[ri]; });
        let fused_atom_rings = fused
            .into_iter()
            .map(|ring| atom_rings[ring].clone())
            .collect::<Vec<_>>();
        // RDKit✔️✔️:       kekulizeFused(mol, frings, atomRanks, maxBackTracks);
        let state = kekulize_fused_system(
            topology,
            &fused_atom_rings,
            all_rings,
            valence,
            atom_ranks,
            max_backtracks,
        )?;
        topology = state.topology;
        if !state.succeeded {
            return Ok(FusedKekulizeState {
                topology,
                succeeded: false,
                problem_atoms: state.problem_atoms,
            });
        }
        // RDKit✔️✔️:       int rix;
        // RDKit✔️✔️:       for (rix = 0; rix < cnrs; ++rix) {
        // RDKit✔️✔️:         if (!fusDone[rix]) {
        // RDKit✔️✔️:           curr = rix;
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (rix == cnrs) {
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
    }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment fused-system dispatch
    Ok(FusedKekulizeState {
        topology,
        succeeded: true,
        problem_atoms: Vec::new(),
    })
}

fn kekulize_fragment(
    topology: &TopologyBlock,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    params: &KekulizeParams,
) -> Result<TopologyBlock, KekulizeError> {
    let prepared = prepare_kekulize_selection(topology, atoms_in_play, bonds_in_play)?;
    if !prepared.found_aromatic {
        return Ok(topology.clone());
    }

    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment ranking and dispatch
    // RDKit✔️✔️:   UINT_VECT atomRanks(mol.getNumAtoms());
    // RDKit✔️✔️:   if (canonical) {
    // RDKit✔️✔️:     Canon::rankFragmentAtoms(mol, atomRanks, atomsToUse, bondsToUse);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // When canonical=false (e.g. during sanitization), we skip the
    // RDKit✔️✔️:     // expensive ranking step and use atom indices directly.  This is
    // RDKit✔️✔️:     // appropriate because sanitization runs *before* stereo perception:
    // RDKit✔️✔️:     // canonical ranking would be based on incomplete chemistry and the
    // RDKit✔️✔️:     // "deterministic" result would be meaningless.  Callers who need a
    // RDKit✔️✔️:     // canonical Kekulé form should call Kekulize() with canonical=true
    // RDKit✔️✔️:     // after the molecule is fully sanitized and stereo has been assigned.
    // RDKit✔️✔️:     std::iota(atomRanks.begin(), atomRanks.end(), 0u);
    // RDKit✔️✔️:   }
    let atom_ranks = if params.canonical {
        rank_fragment_atoms(topology, &prepared.atoms_in_play, &prepared.bonds_in_play)?
    } else {
        (0..topology.atoms.len()).collect()
    };

    // RDKit✔️✔️:   // if any bonds to kekulize then give it a try:
    // RDKit✔️✔️:   if (bondsToUse.any()) {
    let mut working = topology.clone();
    if prepared.bonds_in_play.iter().any(|selected| *selected)
        && !prepared.candidate_atom_rings.is_empty()
    {
        let fused = kekulize_fused_components(
            working,
            &prepared.candidate_atom_rings,
            &prepared.candidate_bond_rings,
            &prepared.rings,
            &prepared.valence,
            &atom_ranks,
            params.max_backtracks,
        )?;
        working = fused.topology;
        if !fused.succeeded {
            return Err(KekulizeError::NotKekulizable {
                problem_atoms: fused.problem_atoms,
            });
        }
    }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment ranking and dispatch

    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment finalization
    // RDKit✔️✔️:   if (markAtomsBonds) {
    // RDKit✔️✔️:     // if we want the atoms and bonds to be marked non-aromatic do
    // RDKit✔️✔️:     // that here.
    if params.mark_atoms_bonds {
        // RDKit✔️✔️:     if (!mol.getRingInfo()->isInitialized()) {
        // RDKit✔️✔️:       MolOps::findSSSR(mol);
        // RDKit✔️✔️:     }
        // Ring membership was initialized during checked selection.
        // RDKit✔️✔️:     for (auto bond : mol.bonds()) {
        // RDKit✔️✔️:       if (bondsToUse[bond->getIdx()]) {
        // RDKit✔️✔️:         bond->setIsAromatic(false);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        for (bond_idx, selected) in prepared.bonds_in_play.iter().copied().enumerate() {
            if selected {
                working.bonds[bond_idx].set_aromatic(false);
            }
        }
        // RDKit✔️✔️:     for (auto atom : mol.atoms()) {
        // RDKit✔️✔️:       if (atomsToUse[atom->getIdx()] && atom->getIsAromatic()) {
        for (atom_idx, selected) in prepared.atoms_in_play.iter().copied().enumerate() {
            if !selected || !working.atoms[atom_idx].is_aromatic() {
                continue;
            }
            let atom_id = AtomId::new(atom_idx);
            // RDKit✔️✔️:         // if we're doing the full molecule and there are aromatic atoms not in
            // RDKit✔️✔️:         // a ring, throw an exception
            // RDKit✔️✔️:         if (atomsToUse.all() && bondsToUse.all() &&
            // RDKit✔️✔️:             !mol.getRingInfo()->numAtomRings(atom->getIdx())) {
            // RDKit✔️✔️:           std::ostringstream errout;
            // RDKit✔️✔️:           errout << "non-ring atom " << atom->getIdx() << " marked aromatic";
            // RDKit✔️✔️:           auto msg = errout.str();
            // RDKit✔️✔️:           BOOST_LOG(rdErrorLog) << msg << std::endl;
            // RDKit✔️✔️:           throw AtomKekulizeException(msg, atom->getIdx());
            // RDKit✔️✔️:         }
            if prepared.atoms_in_play.iter().all(|selected| *selected)
                && prepared.bonds_in_play.iter().all(|selected| *selected)
                && prepared.rings.num_atom_rings(atom_id) == 0
            {
                return Err(KekulizeError::AromaticAtomOutsideRing { atom: atom_id });
            }
            // RDKit✔️✔️:         atom->setIsAromatic(false);
            working.atoms[atom_idx].set_aromatic(false);
            // RDKit✔️✔️:         // make sure "explicit" Hs on things like pyrroles don't hang around
            // RDKit✔️✔️:         // this was Github Issue 141
            // RDKit✔️✔️:         if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
            // RDKit✔️✔️:             atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
            // RDKit✔️✔️:           atom->setNoImplicit(false);
            // RDKit✔️✔️:           atom->setNumExplicitHs(0);
            // RDKit✔️✔️:           atom->updatePropertyCache(false);
            // RDKit✔️✔️:         }
            if matches!(working.atoms[atom_idx].atomic_number(), 7 | 15)
                && working.atoms[atom_idx].formal_charge() == 0
                && working.atoms[atom_idx].explicit_hydrogens() == 1
            {
                working.atoms[atom_idx].set_no_implicit(false);
                working.atoms[atom_idx].set_explicit_hydrogens(0);
            }
            // RDKit✔️✔️:       }
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }

    // RDKit✔️✔️:   // ok some error checking here force a implicit valence
    // RDKit✔️✔️:   // calculation that should do some error checking by itself. In
    // RDKit✔️✔️:   // addition compare them to what they were before kekulizing
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     if (!atomsToUse[atom->getIdx()]) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    let final_valence = crate::assign_valence_with_options_from_parts(
        &working.atoms,
        &working.bonds,
        &working.adjacency,
        ValenceModel::RdkitLike,
        false,
    )?;
    for (atom_idx, selected) in prepared.atoms_in_play.iter().copied().enumerate() {
        if !selected {
            continue;
        }
        let atom_id = AtomId::new(atom_idx);
        // RDKit✔️✔️:     int val = atom->getTotalValence();
        // RDKit✔️✔️:     if (val != valences[atom->getIdx()]) {
        let after = checked_total_valence(&final_valence, atom_id)?;
        let before = prepared.original_total_valences[atom_idx];
        if after != before {
            // RDKit✔️✔️:       std::ostringstream errout;
            // RDKit✔️✔️:       errout << "Kekulization somehow screwed up valence on " << atom->getIdx()
            // RDKit✔️✔️:              << ": " << val << "!=" << valences[atom->getIdx()] << std::endl;
            // RDKit✔️✔️:       auto msg = errout.str();
            // RDKit✔️✔️:       BOOST_LOG(rdErrorLog) << msg << std::endl;
            // RDKit✔️✔️:       throw AtomKekulizeException(msg, atom->getIdx());
            return Err(KekulizeError::PostconditionValenceMismatch {
                atom: atom_id,
                before,
                after,
            });
            // RDKit✔️✔️:     }
        }
        // RDKit✔️✔️:   }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Kekulize.cpp :: KekulizeFragment finalization
    working.validate()?;
    Ok(working)
}

pub fn kekulize(
    topology: &TopologyBlock,
    params: &KekulizeParams,
) -> Result<KekulizeAssignment, KekulizeError> {
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
    Ok(KekulizeAssignment {
        topology: kekulize_fragment(topology, &atoms_in_play, &bonds_in_play, params)?,
    })
}

pub fn kekulize_if_possible(
    topology: &TopologyBlock,
    params: &KekulizeParams,
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
    // The detached input is immutable, so its complete state is the source
    // snapshot and is returned unchanged for the audited sanitize failures.
    // RDKit✔️✔️:   bool res = true;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     Kekulize(mol, markAtomsBonds, canonical, maxBackTracks);
    match kekulize(topology, params) {
        Ok(assignment) => Ok(KekulizeAttempt::Applied(assignment)),
        // RDKit✔️✔️:   } catch (const MolSanitizeException &) {
        Err(KekulizeError::NotKekulizable { problem_atoms }) => {
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

struct CanonRankReadView<'a> {
    atoms: &'a [Atom],
    bonds: &'a [Bond],
    adjacency: &'a AdjacencyList,
    rings: RingInfo,
    valence: ValenceAssignment,
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
            rings,
            valence,
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
    let view = CanonRankReadView::from_topology(topology)?;
    let mut atoms = init_fragment_canon_atoms(
        &view,
        atoms_in_play,
        bonds_in_play,
        CanonicalRankParams::kekulize_fragment_default().include_chirality,
    )?;
    rank_initialized_atoms(
        &view,
        &mut atoms,
        CanonicalRankParams::kekulize_fragment_default(),
    )
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

#[derive(Debug, Clone)]
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
    canonical_ranking_number: i32,
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

fn empty_canon_atom_from_source_atom(atom: &Atom) -> CanonAtom<'static> {
    CanonAtom {
        index: i32::try_from(atom.id().index()).unwrap_or(i32::MAX),
        is_in_play: true,
        degree: 0,
        atomic_number: atom.atomic_number(),
        isotope: atom.isotope().unwrap_or(0),
        atom_map: atom.atom_map().unwrap_or(0),
        canonical_ranking_number: atom
            .prop("_CanonicalRankingNumber")
            .and_then(|value| value.parse::<i32>().ok())
            .unwrap_or(0),
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

fn init_canon_atoms(
    view: &CanonRankReadView<'_>,
    topology: &TopologyBlock,
    include_chirality: bool,
    include_stereo_groups: bool,
) -> Result<Vec<CanonAtom<'static>>, CanonicalRankError> {
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
    let mut atoms =
        init_fragment_canon_atoms(view, &atoms_in_play, &bonds_in_play, include_chirality)?;
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

fn init_fragment_canon_atoms(
    view: &CanonRankReadView<'_>,
    atoms_in_play: &[bool],
    bonds_in_play: &[bool],
    include_chirality: bool,
) -> Result<Vec<CanonAtom<'static>>, CanonicalRankError> {
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
    // RDKit✔️✔️:       atomsi.p_symbol = nullptr;
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
        let atom = &view.atoms[atom_idx];
        atoms[atom_idx].total_num_hs = usize::from(atom.explicit_hydrogens())
            + usize::try_from(view.valence.implicit_hydrogens[atom_idx].max(0))
                .unwrap_or(usize::MAX);
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
        let begin_holder = make_canon_bond_holder(view, bond_idx, end, include_chirality)?;
        let end_holder = make_canon_bond_holder(view, bond_idx, begin, include_chirality)?;
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

fn rank_with_atom_compare_functor_for_kekulize(
    view: &CanonRankReadView<'_>,
    atoms: &mut [CanonAtom<'_>],
    break_ties: bool,
    include_ring_stereo: bool,
    flags: CanonRankFlags,
    order: &mut [usize],
) -> Result<(), CanonicalRankError> {
    // BEGIN RDKIT CPP FUNCTION detail::rankWithFunctor
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void rankWithFunctor(T &ftor, bool breakTies, std::vector<int> &order,
    // RDKit✔️✔️:                      bool useSpecial, bool useChirality, bool includeRingStereo,
    // RDKit✔️✔️:                      const boost::dynamic_bitset<> *atomsInPlay,
    // RDKit✔️✔️:                      const boost::dynamic_bitset<> *bondsInPlay) {
    // RDKit✔️✔️:   PRECONDITION(!order.empty(), "order should not be empty");
    // RDKit✔️✔️:   const ROMol &mol = *ftor.dp_mol;
    // RDKit✔️✔️:   canon_atom *atoms = ftor.dp_atoms;
    // RDKit✔️✔️:   const unsigned int nAts = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   std::vector<int> count(nAts);
    // RDKit✔️✔️:   std::vector<int> next(nAts);
    // RDKit✔️✔️:   std::vector<int> changed(nAts, 1);
    // RDKit✔️✔️:   std::vector<char> touched(nAts, 0);
    // RDKit✔️✔️:   int activeset;
    let mut count = vec![0usize; n_atoms];
    let mut next = vec![-2isize; n_atoms];
    let mut changed = vec![true; n_atoms];
    let mut touched = vec![false; n_atoms];
    let mut active_set = -1isize;
    // RDKit✔️✔️:   CreateSinglePartition(nAts, order, count, atoms);
    create_single_partition_for_kekulize(n_atoms, order, &mut count, atoms);
    // RDKit✔️✔️:   ftor.df_useNbrs = true;
    // RDKit✔️✔️:   ActivatePartitions(nAts, order, count, activeset, next, changed);
    activate_partitions_for_kekulize(
        n_atoms,
        order,
        &count,
        &mut active_set,
        &mut next,
        &mut changed,
    );
    // RDKit✔️✔️:   RefinePartitions(mol, atoms, ftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                    changed, touched);
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
    );
    // RDKit✔️✔️:   bool ties = false;
    // RDKit✔️✔️:   for (unsigned i = 0; i < nAts; ++i) {
    // RDKit✔️✔️:     if (!count[i]) {
    // RDKit✔️✔️:       ties = true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let ties = count.iter().any(|&value| value == 0);
    // RDKit✔️✔️:   if (useChirality && ties && includeRingStereo) {
    // RDKit✔️✔️:     SpecialChiralityAtomCompareFunctor scftor(atoms, mol, atomsInPlay,
    // RDKit✔️✔️:                                               bondsInPlay);
    // RDKit✔️✔️:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit✔️✔️:     RefinePartitions(mol, atoms, scftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                      changed, touched);
    // RDKit✔️✔️:   }
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
        );
    }
    // RDKit✔️✔️:   ties = false;
    // RDKit✔️✔️:   unsigned symRingAtoms = 0;
    // RDKit✔️✔️:   unsigned ringAtoms = 0;
    // RDKit✔️✔️:   bool branchingRingAtom = false;
    // RDKit✔️✔️:   RingInfo *ringInfo = mol.getRingInfo();
    let use_special_symmetry = special_symmetry_rank_refinement_required(view, order, &count);
    // RDKit✔️✔️:   if (useSpecial && ties && ringAtoms > 0 &&
    // RDKit✔️✔️:       static_cast<float>(symRingAtoms) / ringAtoms > 0.5 && branchingRingAtom) {
    // RDKit✔️✔️:     SpecialSymmetryAtomCompareFunctor sftor(atoms, mol, atomsInPlay,
    // RDKit✔️✔️:                                             bondsInPlay);
    // RDKit✔️✔️:     compareRingAtomsConcerningNumNeighbors(atoms, nAts, mol);
    // RDKit✔️✔️:     ActivatePartitions(nAts, order, count, activeset, next, changed);
    // RDKit✔️✔️:     RefinePartitions(mol, atoms, sftor, true, order, count, activeset, next,
    // RDKit✔️✔️:                      changed, touched);
    // RDKit✔️✔️:   }
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
        );
    }
    // RDKit✔️✔️:   if (breakTies) {
    // RDKit✔️✔️:     BreakTies(mol, atoms, ftor, true, order, count, activeset, next, changed,
    // RDKit✔️✔️:               touched);
    // RDKit✔️✔️:   }
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
        );
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION detail::rankWithFunctor
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
) {
    // BEGIN RDKIT CPP FUNCTION RefinePartitions
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void RefinePartitions(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit✔️✔️:                       int mode, std::vector<int> &order,
    // RDKit✔️✔️:                       std::vector<int> &count, int &activeset,
    // RDKit✔️✔️:                       std::vector<int> &next, std::vector<int> &changed,
    // RDKit✔️✔️:                       std::vector<char> &touchedPartitions) {
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   while (activeset != -1) {
    while *active_set != -1 {
        // RDKit✔️✔️:     partition = activeset;
        // RDKit✔️✔️:     activeset = next[partition];
        // RDKit✔️✔️:     next[partition] = -2;
        let partition = usize::try_from(*active_set).expect("active partition is non-negative");
        *active_set = next[partition];
        next[partition] = -2;
        // RDKit✔️✔️:     len = count[partition];
        // RDKit✔️✔️:     offset = atoms[partition].index;
        let len = count[partition];
        let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX);
        // RDKit✔️✔️:     auto start = std::span<int>(&order[offset], len);
        hanoi_sort_order_for_kekulize(
            order,
            offset,
            len,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        );
        // RDKit✔️✔️:     for (int k = 0; k < len; ++k) {
        // RDKit✔️✔️:       changed[start[k]] = 0;
        // RDKit✔️✔️:     }
        for k in 0..len {
            changed[order[offset + k]] = false;
        }
        // RDKit✔️✔️:     index = start[0];
        // RDKit✔️✔️:     for (i = count[index]; i < len; i++) {
        // RDKit✔️✔️:       index = start[i];
        // RDKit✔️✔️:       if (count[index]) {
        // RDKit✔️✔️:         symclass = offset + i;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       atoms[index].index = symclass;
        // RDKit✔️✔️:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
        // RDKit✔️✔️:         changed[atoms[index].nbrIds[j]] = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
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
        // RDKit✔️✔️:     if (mode) {
        // RDKit✔️✔️:       index = start[0];
        // RDKit✔️✔️:       for (i = count[index]; i < len; i++) {
        // RDKit✔️✔️:         index = start[i];
        // RDKit✔️✔️:         for (unsigned j = 0; j < atoms[index].degree; ++j) {
        // RDKit✔️✔️:           unsigned int nbor = atoms[index].nbrIds[j];
        // RDKit✔️✔️:           touchedPartitions[atoms[nbor].index] = 1;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
        // RDKit✔️✔️:         if (touchedPartitions[ii]) {
        // RDKit✔️✔️:           partition = order[ii];
        // RDKit✔️✔️:           if ((count[partition] > 1) && (next[partition] == -2)) {
        // RDKit✔️✔️:             next[partition] = activeset;
        // RDKit✔️✔️:             activeset = partition;
        // RDKit✔️✔️:           }
        // RDKit✔️✔️:           touchedPartitions[ii] = 0;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
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
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RefinePartitions
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
) {
    // BEGIN RDKIT CPP FUNCTION BreakTies
    // RDKit✔️✔️: template <typename CompareFunc>
    // RDKit✔️✔️: void BreakTies(const ROMol &mol, canon_atom *atoms, CompareFunc compar,
    // RDKit✔️✔️:                int mode, std::vector<int> &order, std::vector<int> &count,
    // RDKit✔️✔️:                int &activeset, std::vector<int> &next,
    // RDKit✔️✔️:                std::vector<int> &changed,
    // RDKit✔️✔️:                std::vector<char> &touchedPartitions) {
    // RDKit✔️✔️:   unsigned int nAtoms = mol.getNumAtoms();
    let n_atoms = view.num_atoms();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < nAtoms; i++) {
    let mut i = 0usize;
    while i < n_atoms {
        // RDKit✔️✔️:     partition = order[i];
        // RDKit✔️✔️:     oldPart = atoms[partition].index;
        let partition = order[i];
        let old_part = atoms[partition].index;
        // RDKit✔️✔️:     while (count[partition] > 1) {
        while count[partition] > 1 {
            // RDKit✔️✔️:       len = count[partition];
            // RDKit✔️✔️:       offset = atoms[partition].index + len - 1;
            // RDKit✔️✔️:       index = order[offset];
            // RDKit✔️✔️:       atoms[index].index = offset;
            // RDKit✔️✔️:       count[partition] = len - 1;
            // RDKit✔️✔️:       count[index] = 1;
            let len = count[partition];
            let offset = usize::try_from(atoms[partition].index).unwrap_or(usize::MAX) + len - 1;
            let index = order[offset];
            atoms[index].index = i32::try_from(offset).unwrap_or(i32::MAX);
            count[partition] = len - 1;
            count[index] = 1;
            // RDKit✔️✔️:       if (atoms[index].degree < 1) {
            // RDKit✔️✔️:         continue;
            // RDKit✔️✔️:       }
            if atoms[index].degree < 1 {
                continue;
            }
            // RDKit✔️✔️:       for (unsigned j = 0; j < atoms[index].degree; ++j) {
            // RDKit✔️✔️:         unsigned int nbor = atoms[index].nbrIds[j];
            // RDKit✔️✔️:         touchedPartitions[atoms[nbor].index] = 1;
            // RDKit✔️✔️:         changed[nbor] = 1;
            // RDKit✔️✔️:       }
            for nbr in atoms[index].nbr_ids.iter().copied() {
                let partition_idx = usize::try_from(atoms[nbr].index).unwrap_or(usize::MAX);
                if partition_idx < touched_partitions.len() {
                    touched_partitions[partition_idx] = true;
                }
                changed[nbr] = true;
            }
            // RDKit✔️✔️:       for (unsigned int ii = 0; ii < nAtoms; ++ii) {
            // RDKit✔️✔️:         if (touchedPartitions[ii]) {
            // RDKit✔️✔️:           int npart = order[ii];
            // RDKit✔️✔️:           if ((count[npart] > 1) && (next[npart] == -2)) {
            // RDKit✔️✔️:             next[npart] = activeset;
            // RDKit✔️✔️:             activeset = npart;
            // RDKit✔️✔️:           }
            // RDKit✔️✔️:           touchedPartitions[ii] = 0;
            // RDKit✔️✔️:         }
            // RDKit✔️✔️:       }
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
            // RDKit✔️✔️:       RefinePartitions(mol, atoms, compar, mode, order, count, activeset, next,
            // RDKit✔️✔️:                        changed, touchedPartitions);
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
            );
        }
        // RDKit✔️✔️:     if (atoms[partition].index != oldPart) {
        // RDKit✔️✔️:       i -= 1;
        // RDKit✔️✔️:     }
        if atoms[partition].index != old_part {
            if i > 0 {
                i -= 1;
            }
        } else {
            i += 1;
        }
    }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION BreakTies
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
) {
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
    ) {
        // RDKit✔️✔️:     std::copy(tempVec.begin(), tempVec.end(), base.begin());
        order[offset..offset + len].copy_from_slice(&temp);
    }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hanoisort
}

fn hanoi_order_for_kekulize(
    base: &mut [usize],
    temp: &mut [usize],
    count: &mut [usize],
    changed: &[bool],
    atoms: &mut [CanonAtom<'_>],
    compare_mode: CanonCompareMode,
    flags: CanonRankFlags,
) -> bool {
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
        return false;
    } else if nel == 2 {
        // RDKit✔️✔️:     n1 = base[0];
        // RDKit✔️✔️:     n2 = base[1];
        let n1 = base[0];
        let n2 = base[1];
        // RDKit✔️✔️:     int stat =
        // RDKit✔️✔️:         (/*!changed || */ changed[n1] || changed[n2]) ? compar(n1, n2) : 0;
        let stat = if changed[n1] || changed[n2] {
            compare_canon_atoms_for_kekulize(atoms, n1, n2, compare_mode, flags)
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
        return false;
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
        );
        // RDKit✔️✔️:     if (hanoi(b2, n2, t2, count, changed, compar)) {
        let right_in_temp = hanoi_order_for_kekulize(
            base_right,
            temp_right,
            count,
            changed,
            atoms,
            compare_mode,
            flags,
        );
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
            compare_canon_atoms_for_kekulize(atoms, left_atom, right_atom, compare_mode, flags)
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
                return result;
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
                return result;
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
                return result;
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
                return result;
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
) -> Ordering {
    if matches!(mode, CanonCompareMode::SpecialSymmetry) {
        return compare_special_symmetry_atoms_for_kekulize(atoms, left, right);
    }
    if matches!(mode, CanonCompareMode::SpecialChirality) {
        return compare_special_chirality_atoms_for_kekulize(atoms, left, right);
    }
    if !atom_pair_has_any_in_play_for_kekulize(atoms, left, right) {
        return Ordering::Equal;
    }
    // RDKit✔️✔️:     int v = basecomp(i, j);
    // RDKit✔️✔️:     if (v) {
    // RDKit✔️✔️:       return v;
    // RDKit✔️✔️:     }
    let base_cmp = compare_canon_atom_base_for_kekulize(atoms, left, right, flags);
    if base_cmp != Ordering::Equal {
        return base_cmp;
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
            return cmp;
        }
    }
    // RDKit✔️✔️:       if (dp_atoms[i].bonds.size() < dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return -1;
    // RDKit✔️✔️:       } else if (dp_atoms[i].bonds.size() > dp_atoms[j].bonds.size()) {
    // RDKit✔️✔️:         return 1;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return 0;
    atoms[left].bonds.len().cmp(&atoms[right].bonds.len())
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

fn compare_canon_atom_base_for_kekulize(
    atoms: &[CanonAtom<'_>],
    left: usize,
    right: usize,
    flags: CanonRankFlags,
) -> Ordering {
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
        return cmp;
    }

    if flags.use_non_stereo_ranks {
        cmp = left_atom
            .canonical_ranking_number
            .cmp(&right_atom.canonical_ranking_number);
        if cmp != Ordering::Equal {
            return cmp;
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
            return cmp;
        }
    }

    cmp = left_atom.degree.cmp(&right_atom.degree);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if let (Some(left_symbol), Some(right_symbol)) = (left_atom.p_symbol, right_atom.p_symbol) {
        return left_symbol.cmp(right_symbol);
    }

    cmp = left_atom.atomic_number.cmp(&right_atom.atomic_number);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if flags.use_isotopes {
        cmp = left_atom.isotope.cmp(&right_atom.isotope);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }

    cmp = left_atom.total_num_hs.cmp(&right_atom.total_num_hs);
    if cmp != Ordering::Equal {
        return cmp;
    }

    // RDKit basecomp stores comparison temporaries as unsigned int.
    // Preserve that C++ conversion behavior here so negative formal charges
    // wrap and order after non-negative charges (e.g. -1 > 0 in this step).
    let left_charge = left_atom.formal_charge as i32 as u32;
    let right_charge = right_atom.formal_charge as i32 as u32;
    cmp = left_charge.cmp(&right_charge);
    if cmp != Ordering::Equal {
        return cmp;
    }

    if flags.use_chiral_presence {
        cmp = chiral_presence_for_kekulize(left_atom.chiral_tag)
            .cmp(&chiral_presence_for_kekulize(right_atom.chiral_tag));
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    if flags.use_chirality {
        // RDKit✔️✔️:     ivi = dp_atoms[i].whichStereoGroup;
        // RDKit✔️✔️:     ivj = dp_atoms[j].whichStereoGroup;
        cmp = compare_stereo_group_state_for_kekulize(atoms, left, right);
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    if flags.use_chirality_rings {
        // RDKit✔️✔️:     ivi = getAtomRingNbrCode(i);
        // RDKit✔️✔️:     ivj = getAtomRingNbrCode(j);
        cmp = get_atom_ring_nbr_code_for_kekulize(atoms, left)
            .cmp(&get_atom_ring_nbr_code_for_kekulize(atoms, right));
        if cmp != Ordering::Equal {
            return cmp;
        }
    }
    Ordering::Equal
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
    // BEGIN RDKIT CPP FUNCTION updateAtomNeighborIndex
    // RDKit✔️✔️: void updateAtomNeighborIndex(canon_atom *atoms, std::vector<bondholder> &nbrs) {
    // RDKit✔️✔️:   PRECONDITION(atoms, "bad pointer");
    // RDKit✔️✔️:   for (auto &nbr : nbrs) {
    // RDKit✔️✔️:     unsigned nbrIdx = nbr.nbrIdx;
    // RDKit✔️✔️:     unsigned newSymClass = atoms[nbrIdx].index;
    // RDKit✔️✔️:     nbr.nbrSymClass = newSymClass;
    // RDKit✔️✔️:   }
    let updates = atoms[atom_idx]
        .bonds
        .iter()
        .map(|bond| usize::try_from(atoms[bond.nbr_idx].index).unwrap_or(usize::MAX))
        .collect::<Vec<_>>();
    for (bond, nbr_sym_class) in atoms[atom_idx].bonds.iter_mut().zip(updates) {
        bond.nbr_sym_class = nbr_sym_class;
    }
    // RDKit✔️✔️:   std::sort(nbrs.begin(), nbrs.end(), bondholder::greater);
    // RDKit✔️✔️: }
    let ranks = canon_atom_rank_snapshot(atoms);
    sort_canon_bonds_descending(&mut atoms[atom_idx].bonds, &ranks);
    // END RDKIT CPP FUNCTION updateAtomNeighborIndex
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
    use cosmolkit_model::{AtomSpec, BondSpec};
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
            init_fragment_canon_atoms(&view, &[true, true, false], &[true, true], true).unwrap();
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
            init_fragment_canon_atoms(&view, &[true, true, true], &[true, false], true).unwrap();
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
        let mut atoms = init_fragment_canon_atoms(&view, &[true, true], &[true], true).unwrap();
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
        let mut inactive = init_fragment_canon_atoms(&view, &[false, false], &[], true).unwrap();
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
            prepare_kekulize_selection(&aromatic_ring, &[true, true], &[true; 3]),
            Err(KekulizeError::AtomSelectionLength {
                expected: 3,
                actual: 2
            })
        ));
        assert!(matches!(
            prepare_kekulize_selection(&aromatic_ring, &[true; 3], &[true; 2]),
            Err(KekulizeError::BondSelectionLength {
                expected: 3,
                actual: 2
            })
        ));

        let empty = prepare_kekulize_selection(&aromatic_ring, &[false; 3], &[true; 3]).unwrap();
        assert!(!empty.found_aromatic);
        assert!(empty.candidate_atom_rings.is_empty());
        assert_eq!(empty.original_total_valences, vec![0; 3]);

        let crossing =
            prepare_kekulize_selection(&aromatic_ring, &[true, true, false], &[true; 3]).unwrap();
        assert!(crossing.found_aromatic);
        assert!(crossing.candidate_atom_rings.is_empty());
        assert!(crossing.candidate_bond_rings.is_empty());

        let query_bond = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Aromatic)
            .with_aromatic(true)
            .with_prop("_MolFileBondQuery", "1")
            .unwrap();
        let query_ring = topology(
            (0..3)
                .map(|id| atom(id, AtomSpec::new(Element::C).with_aromatic(true)))
                .collect(),
            vec![
                bond_with_spec(0, query_bond),
                aromatic_bond(1, 1, 2),
                aromatic_bond(2, 2, 0),
            ],
        );
        let query = prepare_kekulize_selection(&query_ring, &[true; 3], &[true; 3]).unwrap();
        assert_eq!(query.bonds_in_play, vec![false, true, true]);
        assert!(query.candidate_atom_rings.is_empty());

        let complex_query = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("_MolFileBondQueryComplex", "recursive")
            .unwrap();
        let complex = topology(
            vec![
                atom(0, AtomSpec::new(Element::C)),
                atom(1, AtomSpec::new(Element::C)),
            ],
            vec![bond_with_spec(0, complex_query)],
        );
        assert!(matches!(
            prepare_kekulize_selection(&complex, &[true; 2], &[true]),
            Err(KekulizeError::UnsupportedQueryState {
                bond,
                detail: "concrete TopologyBlock cannot represent a recursive or composite bond query"
            }) if bond == BondId::new(0)
        ));
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
    fn fused_dummy_subsets_follow_source_bit_order_and_reject_counter_overflow() {
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

        assert!(matches!(
            QuestionEnumerator::new(vec![AtomId::new(0); u32::BITS as usize]),
            Err(KekulizeError::QuestionSubsetOverflow {
                questions: 32,
                bit_width: 32
            })
        ));
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
            prepare_kekulize_selection(&all_dummy, &[true; 3], &[true; 3]).unwrap();
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
            prepare_kekulize_selection(&crossing, &[true, true, false], &[true; 3]).unwrap();
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
        let selection = prepare_kekulize_selection(&graph, &[true; 16], &[true; 17]).unwrap();
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
