use std::collections::{BTreeMap, BTreeSet, VecDeque};

use crate::{
    AtomId, BondDirection, BondId, BondOrder, Molecule, RingFindingError, RingInfo, ValenceError,
    ValenceModel, assign_valence_with_options, canon_rank::FragmentRankScope,
    canon_rank::rank_fragment_atoms_for_kekulize, molecule::TopologyBlock, symmetrize_sssr,
};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum KekulizeError {
    #[error("can't kekulize mol. unkekulized atoms: {atoms:?}")]
    UnkekulizedAtoms { atoms: Vec<AtomId> },
    #[error("non-ring atom {atom} marked aromatic")]
    NonRingAromaticAtom { atom: AtomId },
    #[error("incomplete RDKit kekulize port for {branch}: {reason}")]
    ProtocolDebt {
        branch: &'static str,
        reason: &'static str,
    },
    #[error("kekulization changed valence on atom {atom}: {before}!={after}")]
    ValenceChanged {
        atom: AtomId,
        before: i32,
        after: i32,
    },
    #[error("kekulize fragment bitset size mismatch: atoms={atoms}, bonds={bonds}")]
    FragmentBitsetSizeMismatch { atoms: usize, bonds: usize },
    #[error("canonical rank {kind} symbol size mismatch: expected={expected}, actual={actual}")]
    CanonicalRankSymbolSizeMismatch {
        kind: &'static str,
        expected: usize,
        actual: usize,
    },
    #[error(transparent)]
    RingFinding(#[from] RingFindingError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct KekulizeAssignment {
    bond_orders: Vec<Option<BondOrder>>,
    bond_aromatic_flags: Vec<Option<bool>>,
    bond_directions: Vec<Option<BondDirection>>,
    atom_aromatic_flags: Vec<Option<bool>>,
    atom_explicit_hydrogens: Vec<Option<u8>>,
    atom_no_implicit: Vec<Option<bool>>,
}

struct KekulizeCandidateState {
    all_atoms: Vec<AtomId>,
    candidate_atoms: BTreeSet<AtomId>,
    aromatic_edges: Vec<(BondId, AtomId, AtomId)>,
    questions: Vec<AtomId>,
    done: Vec<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct KekulizeRing {
    atoms: Vec<AtomId>,
    bonds: Vec<BondId>,
    source_ring: usize,
}

impl KekulizeAssignment {
    fn empty(molecule: &Molecule) -> Self {
        Self {
            bond_orders: vec![None; molecule.num_bonds()],
            bond_aromatic_flags: vec![None; molecule.num_bonds()],
            bond_directions: vec![None; molecule.num_bonds()],
            atom_aromatic_flags: vec![None; molecule.num_atoms()],
            atom_explicit_hydrogens: vec![None; molecule.num_atoms()],
            atom_no_implicit: vec![None; molecule.num_atoms()],
        }
    }

    pub(crate) fn bond_order(&self, bond: BondId) -> Option<BondOrder> {
        self.bond_orders[bond.index()]
    }

    pub(crate) fn bond_aromatic_flag(&self, bond: BondId) -> Option<bool> {
        self.bond_aromatic_flags[bond.index()]
    }

    pub(crate) fn bond_direction(&self, bond: BondId) -> Option<BondDirection> {
        self.bond_directions[bond.index()]
    }

    pub(crate) fn atom_aromatic_flag(&self, atom: AtomId) -> Option<bool> {
        self.atom_aromatic_flags[atom.index()]
    }

    pub(crate) fn atom_explicit_hydrogens(&self, atom: AtomId) -> Option<u8> {
        self.atom_explicit_hydrogens[atom.index()]
    }

    pub(crate) fn atom_no_implicit(&self, atom: AtomId) -> Option<bool> {
        self.atom_no_implicit[atom.index()]
    }

    #[cfg(test)]
    fn is_empty(&self) -> bool {
        self.bond_orders.iter().all(Option::is_none)
            && self.bond_aromatic_flags.iter().all(Option::is_none)
            && self.bond_directions.iter().all(Option::is_none)
            && self.atom_aromatic_flags.iter().all(Option::is_none)
            && self.atom_explicit_hydrogens.iter().all(Option::is_none)
            && self.atom_no_implicit.iter().all(Option::is_none)
    }
}

pub(crate) fn apply_kekulize_assignment(
    topology: &mut TopologyBlock,
    assignment: &KekulizeAssignment,
) -> bool {
    let mut changed = false;
    for bond in &mut topology.bonds {
        if let Some(next) = assignment.bond_order(bond.id()) {
            if bond.order() != next {
                bond.set_order(next);
                changed = true;
            }
        }
        if let Some(next) = assignment.bond_aromatic_flag(bond.id()) {
            if bond.is_aromatic() != next {
                bond.set_aromatic(next);
                changed = true;
            }
        }
        if let Some(next) = assignment.bond_direction(bond.id()) {
            if bond.direction() != next {
                bond.set_direction(next);
                changed = true;
            }
        }
    }
    for atom in &mut topology.atoms {
        if let Some(next) = assignment.atom_aromatic_flag(atom.id()) {
            if atom.is_aromatic() != next {
                atom.set_aromatic(next);
                changed = true;
            }
        }
        if let Some(next) = assignment.atom_explicit_hydrogens(atom.id()) {
            if atom.explicit_hydrogens() != next {
                atom.set_explicit_hydrogens(next);
                changed = true;
            }
        }
        if let Some(next) = assignment.atom_no_implicit(atom.id()) {
            if atom.no_implicit() != next {
                atom.set_no_implicit(next);
                changed = true;
            }
        }
    }
    changed
}

pub(crate) fn kekulize_assignment(
    molecule: &Molecule,
    rings: Option<&RingInfo>,
    clear_aromatic_flags: bool,
    canonical: bool,
    max_backtracks: usize,
) -> Result<KekulizeAssignment, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::Kekulize
    // RDKit❗✔️: void Kekulize(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit❗✔️:               unsigned int maxBackTracks) {
    // RDKit❗✔️:   boost::dynamic_bitset<> atomsToUse(mol.getNumAtoms());
    // RDKit❗✔️:   atomsToUse.set();
    // RDKit❗✔️:   boost::dynamic_bitset<> bondsToUse(mol.getNumBonds());
    // RDKit❗✔️:   bondsToUse.set();
    // RDKit❗✔️:   details::KekulizeFragment(mol, atomsToUse, bondsToUse, markAtomsBonds,
    // RDKit❗✔️:                             canonical, maxBackTracks);
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::Kekulize
    let owned_rings;
    let rings = if let Some(rings) = rings {
        rings
    } else {
        owned_rings = symmetrize_sssr(molecule)?;
        &owned_rings
    };
    let atoms_to_use = vec![true; molecule.num_atoms()];
    let bonds_to_use = vec![true; molecule.num_bonds()];
    kekulize_fragment_assignment(
        molecule,
        rings,
        &atoms_to_use,
        bonds_to_use,
        clear_aromatic_flags,
        canonical,
        max_backtracks,
    )
}

fn kekulize_fragment_assignment(
    molecule: &Molecule,
    rings: &RingInfo,
    atoms_to_use: &[bool],
    mut bonds_to_use: Vec<bool>,
    clear_aromatic_flags: bool,
    canonical: bool,
    max_backtracks: usize,
) -> Result<KekulizeAssignment, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION details::KekulizeFragment
    // RDKit❗✔️: void KekulizeFragment(RWMol &mol, const boost::dynamic_bitset<> &atomsToUse,
    // RDKit❗✔️:                       boost::dynamic_bitset<> bondsToUse, bool markAtomsBonds,
    // RDKit❗✔️:                       bool canonical, unsigned int maxBackTracks) {
    // RDKit❗✔️:   PRECONDITION(atomsToUse.size() == mol.getNumAtoms(),
    // RDKit❗✔️:                "atomsToUse is wrong size");
    // RDKit❗✔️:   PRECONDITION(bondsToUse.size() == mol.getNumBonds(),
    // RDKit❗✔️:                "bondsToUse is wrong size");
    // RDKit❗✔️:   if (atomsToUse.none()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    if atoms_to_use.len() != molecule.num_atoms() || bonds_to_use.len() != molecule.num_bonds() {
        return Err(KekulizeError::FragmentBitsetSizeMismatch {
            atoms: atoms_to_use.len(),
            bonds: bonds_to_use.len(),
        });
    }
    let mut assignment = KekulizeAssignment::empty(molecule);
    if !atoms_to_use.iter().any(|selected| *selected) {
        return Ok(assignment);
    }

    // RDKit❗✔️:   bool foundAromatic = false;
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bondsToUse[bond->getIdx()]) {
    // RDKit❗✔️:       if (QueryOps::hasBondTypeQuery(*bond)) {
    // RDKit❗✔️:         bondsToUse[bond->getIdx()] = 0;
    // RDKit❗✔️:       } else if (bond->getIsAromatic()) {
    // RDKit❗✔️:         foundAromatic = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    let mut found_aromatic = false;
    for bond in molecule.bonds() {
        if bonds_to_use[bond.id().index()] {
            if bond
                .query()
                .is_some_and(crate::valence::has_bond_type_query)
            {
                bonds_to_use[bond.id().index()] = false;
            } else if bond.is_aromatic() {
                found_aromatic = true;
            }
        }
    }

    // RDKit❗✔️:   for (auto atom : mol.atoms()) {
    // RDKit❗✔️:     atom->calcImplicitValence(false);
    // RDKit❗✔️:     valences[atom->getIdx()] = atom->getTotalValence();
    // RDKit❗✔️:     if (isAromaticAtom(*atom)) {
    // RDKit❗✔️:       foundAromatic = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    let pre_kek_valence = assign_valence_with_options(molecule, ValenceModel::RdkitLike, false)?;
    if molecule.atoms().iter().any(crate::Atom::is_aromatic) {
        found_aromatic = true;
    }
    if !found_aromatic {
        return Ok(assignment);
    }

    // RDKit❗✔️:   UINT_VECT atomRanks(mol.getNumAtoms());
    // RDKit❗✔️:   if (canonical) {
    // RDKit❗✔️:     Canon::rankFragmentAtoms(mol, atomRanks, atomsToUse, bondsToUse);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     std::iota(atomRanks.begin(), atomRanks.end(), 0u);
    // RDKit❗✔️:   }
    let atom_ranks = if canonical {
        rank_fragment_atoms_for_kekulize(
            molecule,
            FragmentRankScope::new(atoms_to_use, &bonds_to_use),
        )?
    } else {
        (0..molecule.num_atoms()).collect::<Vec<_>>()
    };

    // RDKit❗✔️:   if (bondsToUse.any()) {
    // RDKit❗✔️:     if (!mol.getRingInfo()->isInitialized()) {
    // RDKit❗✔️:       MolOps::findSSSR(mol, allringsSSSR);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   if (bondsToUse.any()) {
    // RDKit❗✔️:     boost::dynamic_bitset<> wedgedAtoms(numAtoms);
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       if (bondsToUse[bond->getIdx()] &&
    // RDKit❗✔️:           (bond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit❗✔️:            bond->getBondDir() == Bond::BEGINDASH)) {
    // RDKit❗✔️:         wedgedAtoms.set(bond->getBeginAtomIdx());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    let mut wedged_atoms = BTreeSet::new();
    for bond in molecule.bonds() {
        if bonds_to_use[bond.id().index()]
            && matches!(
                bond.direction(),
                BondDirection::BeginWedge | BondDirection::BeginDash
            )
        {
            wedged_atoms.insert(bond.begin());
        }
    }

    // RDKit❗✔️:     const VECT_INT_VECT &allrings =
    // RDKit❗✔️:         allringsSSSR.empty() ? mol.getRingInfo()->atomRings() : allringsSSSR;
    // RDKit❗✔️:     std::deque<INT_VECT> tmpRings;
    // RDKit❗✔️:     auto containsNonDummy = [&atomsToUse, &dummyAts](const INT_VECT &ring) {
    // RDKit❗✔️:       bool ringOk = false;
    // RDKit❗✔️:       for (auto ai : ring) {
    // RDKit❗✔️:         if (!atomsToUse[ai]) {
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (!dummyAts[ai]) {
    // RDKit❗✔️:           ringOk = true;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       return ringOk;
    // RDKit❗✔️:     };
    // RDKit❗✔️:     for (const auto &ring : allrings) {
    // RDKit❗✔️:       if (containsNonDummy(ring)) {
    // RDKit❗✔️:         unsigned int startPos = 0;
    // RDKit❗✔️:         bool hasWedge = false;
    // RDKit❗✔️:         for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit❗✔️:           if (wedgedAtoms[ring[ri]]) {
    // RDKit❗✔️:             startPos = ri;
    // RDKit❗✔️:             hasWedge = true;
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         INT_VECT nring(ring.size());
    // RDKit❗✔️:         for (auto ri = 0u; ri < ring.size(); ++ri) {
    // RDKit❗✔️:           nring[ri] = ring.at((ri + startPos) % ring.size());
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (!hasWedge) {
    // RDKit❗✔️:           tmpRings.push_back(nring);
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           tmpRings.push_front(nring);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     VECT_INT_VECT arings;
    // RDKit❗✔️:     arings.reserve(allrings.size());
    // RDKit❗✔️:     arings.insert(arings.end(), tmpRings.begin(), tmpRings.end());
    // RDKit❗✔️:     VECT_INT_VECT allbrings;
    // RDKit❗✔️:     RingUtils::convertToBonds(arings, allbrings, mol);
    // RDKit❗✔️:     VECT_INT_VECT brings;
    // RDKit❗✔️:     brings.reserve(allbrings.size());
    // RDKit❗✔️:     auto copyBondRingsWithinFragment = [&bondsToUse](const INT_VECT &ring) {
    // RDKit❗✔️:       return std::all_of(ring.begin(), ring.end(), [&bondsToUse](const int bi) {
    // RDKit❗✔️:         return bondsToUse[bi];
    // RDKit❗✔️:       });
    // RDKit❗✔️:     };
    // RDKit❗✔️:     VECT_INT_VECT aringsRemaining;
    // RDKit❗✔️:     aringsRemaining.reserve(arings.size());
    // RDKit❗✔️:     for (unsigned i = 0; i < allbrings.size(); ++i) {
    // RDKit❗✔️:       if (copyBondRingsWithinFragment(allbrings[i])) {
    // RDKit❗✔️:         brings.push_back(allbrings[i]);
    // RDKit❗✔️:         aringsRemaining.push_back(arings[i]);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     arings = std::move(aringsRemaining);
    let kekulize_rings =
        ordered_kekulize_rings(molecule, rings, atoms_to_use, &bonds_to_use, &wedged_atoms);
    let ring_systems = fused_ring_systems(&kekulize_rings);
    for system in ring_systems {
        let ring_bond_set = system
            .iter()
            .flat_map(|&ring_idx| kekulize_rings[ring_idx].bonds.iter().copied())
            .collect::<BTreeSet<_>>();
        let ring_atom_set = system
            .iter()
            .flat_map(|&ring_idx| kekulize_rings[ring_idx].atoms.iter().copied())
            .collect::<BTreeSet<_>>();

        if ring_bond_set
            .iter()
            .any(|bond| molecule.bonds()[bond.index()].is_aromatic())
            || ring_atom_set
                .iter()
                .any(|atom| molecule.atoms()[atom.index()].is_aromatic())
        {
            kekulize_fused_assignment(
                molecule,
                rings,
                &kekulize_rings,
                &system,
                &atom_ranks,
                max_backtracks,
                &mut assignment,
            )?;
        }
    }

    // RDKit❗✔️:   if (markAtomsBonds) {
    // RDKit❗✔️:     for (auto bond : mol.bonds()) {
    // RDKit❗✔️:       if (bondsToUse[bond->getIdx()]) {
    // RDKit❗✔️:         bond->setIsAromatic(false);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (auto atom : mol.atoms()) {
    // RDKit❗✔️:       if (atomsToUse[atom->getIdx()] && atom->getIsAromatic()) {
    if clear_aromatic_flags {
        for bond in molecule.bonds() {
            if bonds_to_use[bond.id().index()] && bond.is_aromatic() {
                assignment.bond_aromatic_flags[bond.id().index()] = Some(false);
            }
        }
        for atom in molecule.atoms() {
            if atoms_to_use[atom.id().index()] && atom.is_aromatic() {
                if rings.num_atom_rings(atom.id()) == 0 {
                    return Err(KekulizeError::NonRingAromaticAtom { atom: atom.id() });
                }
                assignment.atom_aromatic_flags[atom.id().index()] = Some(false);
                // RDKit❗✔️:         if ((atom->getAtomicNum() == 7 || atom->getAtomicNum() == 15) &&
                // RDKit❗✔️:             atom->getFormalCharge() == 0 && atom->getNumExplicitHs() == 1) {
                // RDKit❗✔️:           atom->setNoImplicit(false);
                // RDKit❗✔️:           atom->setNumExplicitHs(0);
                if matches!(atom.atomic_number(), 7 | 15)
                    && atom.formal_charge() == 0
                    && atom.explicit_hydrogens() == 1
                {
                    assignment.atom_no_implicit[atom.id().index()] = Some(false);
                    assignment.atom_explicit_hydrogens[atom.id().index()] = Some(0);
                }
            }
        }
    }
    let mut kekulized = molecule.clone();
    apply_kekulize_assignment(kekulized.topology_block_mut(), &assignment);
    let post_kek_valence = assign_valence_with_options(&kekulized, ValenceModel::RdkitLike, false)?;
    for atom in molecule.atoms() {
        let idx = atom.id().index();
        let before =
            pre_kek_valence.explicit_valence[idx] + pre_kek_valence.implicit_hydrogens[idx];
        let after =
            post_kek_valence.explicit_valence[idx] + post_kek_valence.implicit_hydrogens[idx];
        if before != after {
            return Err(KekulizeError::ValenceChanged {
                atom: atom.id(),
                before,
                after,
            });
        }
    }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION details::KekulizeFragment
    Ok(assignment)
}

fn kekulize_fused_assignment(
    molecule: &Molecule,
    rings: &RingInfo,
    kekulize_rings: &[KekulizeRing],
    fused_ring_indices: &[usize],
    atom_ranks: &[usize],
    max_backtracks: usize,
    assignment: &mut KekulizeAssignment,
) -> Result<(), KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION kekulizeFused
    // RDKit❗✔️: void kekulizeFused(RWMol &mol, const VECT_INT_VECT &arings,
    // RDKit❗✔️:                    const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit❗✔️:   INT_VECT allAtms;
    // RDKit❗✔️:   Union(arings, allAtms);
    // RDKit❗✔️:   INT_VECT done;
    // RDKit❗✔️:   INT_VECT questions;
    // RDKit❗✔️:   boost::dynamic_bitset<> dBndCands(nats);
    // RDKit❗✔️:   boost::dynamic_bitset<> dBndAdds(nbnds);
    // RDKit❗✔️:   markDbondCands(mol, allAtms, dBndCands, questions, done);
    let state = mark_double_bond_candidates(
        molecule,
        rings,
        kekulize_rings,
        fused_ring_indices,
        assignment,
    )?;

    // RDKit❗✔️:   auto kekulized =
    // RDKit❗✔️:       kekulizeWorker(mol, allAtms, dBndCands, dBndAdds, done, atomRanks,
    // RDKit❗✔️:                      maxBackTracks);
    let matched_bonds = if let Some(matched_bonds) = kekulize_worker_matching(
        molecule,
        &state,
        &state.done,
        atom_ranks,
        max_backtracks,
        &BTreeSet::new(),
    ) {
        matched_bonds
    } else if let Some(matched_bonds) = permute_dummies_and_match(
        molecule,
        &state,
        &state.questions,
        atom_ranks,
        max_backtracks,
    ) {
        matched_bonds
    } else {
        return Err(KekulizeError::UnkekulizedAtoms {
            atoms: state.candidate_atoms.iter().copied().collect(),
        });
    };

    for (bond, _, _) in state.aromatic_edges {
        assignment.bond_orders[bond.index()] = Some(if matched_bonds.contains(&bond) {
            BondOrder::Double
        } else {
            BondOrder::Single
        });
        if matched_bonds.contains(&bond)
            && molecule.bonds()[bond.index()].direction() != BondDirection::None
        {
            assignment.bond_directions[bond.index()] = Some(BondDirection::None);
        }
    }
    // RDKit❗✔️:   if (!kekulized) {
    // RDKit❗✔️:     throw KekulizeException(msg, problemAtoms);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION kekulizeFused
    Ok(())
}

fn mark_double_bond_candidates(
    molecule: &Molecule,
    rings: &RingInfo,
    kekulize_rings: &[KekulizeRing],
    fused_ring_indices: &[usize],
    assignment: &mut KekulizeAssignment,
) -> Result<KekulizeCandidateState, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION markDbondCands
    // RDKit❗✔️: void markDbondCands(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗✔️:                     boost::dynamic_bitset<> &dBndCands, INT_VECT &questions,
    // RDKit❗✔️:                     INT_VECT &done) {
    // RDKit❗✔️:   bool hasAromaticOrDummyAtom =
    // RDKit❗✔️:       std::any_of(allAtms.begin(), allAtms.end(), [&mol](int allAtm) {
    // RDKit❗✔️:         return (!mol.getAtomWithIdx(allAtm)->getAtomicNum() ||
    // RDKit❗✔️:                 isAromaticAtom(*mol.getAtomWithIdx(allAtm)));
    // RDKit❗✔️:       });
    let mut all_atoms = Vec::new();
    let mut seen_atoms = BTreeSet::new();
    for &ring_idx in fused_ring_indices {
        for &atom_id in &kekulize_rings[ring_idx].atoms {
            if seen_atoms.insert(atom_id) {
                all_atoms.push(atom_id);
            }
        }
    }
    let mut ring_bonds = Vec::new();
    let mut seen_bonds = BTreeSet::new();
    for &ring_idx in fused_ring_indices {
        for &bond_id in &kekulize_rings[ring_idx].bonds {
            if seen_bonds.insert(bond_id) {
                ring_bonds.push(bond_id);
            }
        }
    }
    let all_atom_set = all_atoms.iter().copied().collect::<BTreeSet<_>>();
    if !all_atoms.iter().any(|atom| {
        molecule.atoms()[atom.index()].atomic_number() == 0
            || molecule.atoms()[atom.index()].is_aromatic()
    }) {
        return Ok(KekulizeCandidateState {
            all_atoms,
            candidate_atoms: BTreeSet::new(),
            aromatic_edges: Vec::new(),
            questions: Vec::new(),
            done: Vec::new(),
        });
    }
    let valence = assign_valence_with_options(molecule, ValenceModel::RdkitLike, false)?;
    let fused_source_ring_set = fused_ring_indices
        .iter()
        .map(|&ring_idx| kekulize_rings[ring_idx].source_ring)
        .collect::<BTreeSet<_>>();
    let mut is_ring_not_cand = BTreeSet::new();
    for &ring_idx in fused_ring_indices {
        let mut ring_is_candidate = false;
        let source_ring = kekulize_rings[ring_idx].source_ring;
        for atom in &rings.atom_rings()[source_ring] {
            if molecule.atoms()[atom.index()].is_aromatic() && rings.num_atom_rings(*atom) == 1 {
                ring_is_candidate = true;
                break;
            }
        }
        if !ring_is_candidate {
            is_ring_not_cand.insert(source_ring);
        }
    }
    let mut candidate_atoms = BTreeSet::<AtomId>::new();
    let mut aromatic_edges = Vec::<(BondId, AtomId, AtomId)>::new();
    let mut questions = Vec::<AtomId>::new();
    let mut done = Vec::<AtomId>::new();
    // RDKit❗✔️:   std::vector<Bond *> makeSingle;
    for bond_id in ring_bonds {
        let bond = &molecule.bonds()[bond_id.index()];
        if bond.is_aromatic()
            && matches!(
                bond.order(),
                BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
            )
        {
            assignment.bond_orders[bond_id.index()] = Some(BondOrder::Single);
            aromatic_edges.push((bond_id, bond.begin(), bond.end()));
        }
    }
    for &atom_id in &all_atoms {
        let atom = &molecule.atoms()[atom_id.index()];
        if atom.atomic_number() != 0 && !atom.is_aromatic() {
            done.push(atom_id);
            continue;
        }
        let mut sbo = 0i32;
        let mut n_to_ignore = 0usize;
        let mut non_ar_non_dummy_nbr = 0usize;
        for bond in molecule.bonds().iter().filter(|bond| {
            (bond.begin() == atom_id && all_atom_set.contains(&bond.end()))
                || (bond.end() == atom_id && all_atom_set.contains(&bond.begin()))
        }) {
            let other = if bond.begin() == atom_id {
                bond.end()
            } else {
                bond.begin()
            };
            let other_atom = &molecule.atoms()[other.index()];
            if other_atom.atomic_number() != 0 && !other_atom.is_aromatic() {
                non_ar_non_dummy_nbr += 1;
            }
            if bond.is_aromatic()
                && matches!(
                    bond.order(),
                    BondOrder::Single | BondOrder::Double | BondOrder::Aromatic
                )
            {
                sbo += 1;
            } else {
                let bond_contrib =
                    crate::valence::bond_valence_contrib(bond, atom_id)?.round() as i32;
                sbo += bond_contrib;
                if bond_contrib == 0 {
                    n_to_ignore += 1;
                }
            }
        }
        let num_atom_rings = rings.num_atom_rings(atom_id);
        let num_non_cand_rings = rings
            .atom_members(atom_id)
            .iter()
            .filter(|ring_idx| {
                fused_source_ring_set.contains(ring_idx) && is_ring_not_cand.contains(ring_idx)
            })
            .count();
        if atom.atomic_number() == 0
            && non_ar_non_dummy_nbr < num_atom_rings
            && num_non_cand_rings < num_atom_rings
        {
            candidate_atoms.insert(atom_id);
            questions.push(atom_id);
            continue;
        }
        if atom.atomic_number() == 0 {
            continue;
        }
        sbo += i32::from(atom.explicit_hydrogens()) + valence.implicit_hydrogens[atom_id.index()];
        let mut dv = rdkit_kekulize_default_valence(atom.atomic_number())?;
        let mut chrg = atom.formal_charge();
        if rdkit_is_early_atom(atom.atomic_number()) {
            chrg = -chrg;
        }
        if atom.atomic_number() == 6 && chrg > 0 {
            chrg = -chrg;
        }
        dv += i32::from(chrg);
        let tbo =
            valence.explicit_valence[atom_id.index()] + valence.implicit_hydrogens[atom_id.index()];
        let n_radicals = i32::from(atom.radical_electrons());
        let total_degree = i32::try_from(
            molecule
                .bonds()
                .iter()
                .filter(|bond| bond.begin() == atom_id || bond.end() == atom_id)
                .count(),
        )
        .unwrap_or(i32::MAX)
            + valence.implicit_hydrogens[atom_id.index()]
            - i32::try_from(n_to_ignore).unwrap_or(i32::MAX);
        let valence_list = crate::rdkit_valence_list(atom.atomic_number())?.ok_or(
            ValenceError::UnsupportedBranch {
                reason: "kekulize valence list atomic number out of range",
            },
        )?;
        let mut vi = 1usize;
        while tbo > dv && vi < valence_list.len() && valence_list[vi] > 0 {
            dv = valence_list[vi] + i32::from(chrg);
            vi += 1;
        }
        if tbo == 5
            && sbo == 4
            && dv == 3
            && total_degree == 3
            && n_radicals == 0
            && chrg == 0
            && atom.explicit_hydrogens() == 0
            && valence.implicit_hydrogens[atom_id.index()] == 0
            && matches!(atom.atomic_number(), 7 | 15 | 33)
        {
            dv = 5;
        }
        if total_degree + n_radicals >= dv {
            continue;
        }
        if dv == (sbo + 1 + n_radicals)
            || (n_radicals == 0 && atom.no_implicit() && dv == (sbo + 2))
        {
            candidate_atoms.insert(atom_id);
        }
    }
    // RDKit❗✔️:   for (auto &bi : makeSingle) {
    // RDKit❗✔️:     bi->setBondType(Bond::SINGLE);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION markDbondCands
    Ok(KekulizeCandidateState {
        all_atoms,
        candidate_atoms,
        aromatic_edges,
        questions,
        done,
    })
}

fn kekulize_worker_matching(
    molecule: &Molecule,
    state: &KekulizeCandidateState,
    initial_done: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: usize,
    switched_off: &BTreeSet<AtomId>,
) -> Option<BTreeSet<BondId>> {
    // BEGIN RDKIT CPP FUNCTION kekulizeWorker
    // RDKit❗✔️: bool kekulizeWorker(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗✔️:                     boost::dynamic_bitset<> dBndCands,
    // RDKit❗✔️:                     boost::dynamic_bitset<> dBndAdds, INT_VECT done,
    // RDKit❗✔️:                     const UINT_VECT &atomRanks, unsigned int maxBackTracks) {
    // RDKit❗✔️:   INT_DEQUE astack;
    // RDKit❗✔️:   INT_INT_DEQ_MAP options;
    // RDKit❗✔️:   int lastOpt = -1;
    // RDKit❗✔️:   boost::dynamic_bitset<> localBondsAdded(mol.getNumBonds());
    // RDKit❗✔️:   boost::dynamic_bitset<> inAllAtms(mol.getNumAtoms());
    // RDKit❗✔️:   for (int allAtm : allAtms) {
    // RDKit❗✔️:     inAllAtms.set(allAtm);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   auto lessByRank = [&atomRanks](int a, int b) {
    // RDKit❗✔️:     const auto ra = atomRanks.at(static_cast<unsigned int>(a));
    // RDKit❗✔️:     const auto rb = atomRanks.at(static_cast<unsigned int>(b));
    // RDKit❗✔️:     return (ra < rb) || (ra == rb && a < b);
    // RDKit❗✔️:   };
    // RDKit❗✔️:   boost::dynamic_bitset<> wedgeEndAtoms(mol.getNumAtoms());
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit❗✔️:         bond->getBondDir() == Bond::BondDir::BEGINDASH) {
    // RDKit❗✔️:       const auto endIdx = bond->getEndAtomIdx();
    // RDKit❗✔️:       if (inAllAtms.test(endIdx)) {
    // RDKit❗✔️:         wedgeEndAtoms.set(endIdx);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   INT_VECT sortedAtms(allAtms);
    // RDKit❗✔️:   std::sort(sortedAtms.begin(), sortedAtms.end(),
    // RDKit❗✔️:             [&wedgeEndAtoms, &lessByRank](int a, int b) {
    // RDKit❗✔️:               const bool wa = wedgeEndAtoms.test(a);
    // RDKit❗✔️:               const bool wb = wedgeEndAtoms.test(b);
    // RDKit❗✔️:               if (wa != wb) {
    // RDKit❗✔️:                 return wa;
    // RDKit❗✔️:               }
    // RDKit❗✔️:               return lessByRank(a, b);
    // RDKit❗✔️:             });
    // RDKit❗✔️:   int curr = -1;
    // RDKit❗✔️:   INT_DEQUE btmoves;
    // RDKit❗✔️:   unsigned int numBT = 0;
    // RDKit❗✔️:   while ((done.size() < sortedAtms.size()) || !astack.empty()) {
    // RDKit❗✔️:     if (astack.size() > 0) {
    // RDKit❗✔️:       curr = astack.front();
    // RDKit❗✔️:       astack.pop_front();
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       curr = -1;
    // RDKit❗✔️:       for (int allAtm : sortedAtms) {
    // RDKit❗✔️:         if (std::find(done.begin(), done.end(), allAtm) == done.end()) {
    // RDKit❗✔️:           curr = allAtm;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     CHECK_INVARIANT(curr >= 0, "starting point not found");
    // RDKit❗✔️:     done.push_back(curr);
    // RDKit❗✔️:     INT_DEQUE opts;
    // RDKit❗✔️:     bool cCand = false;
    // RDKit❗✔️:     if (dBndCands[curr]) {
    // RDKit❗✔️:       cCand = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     int ncnd;
    // RDKit❗✔️:     if (options.find(curr) != options.end()) {
    // RDKit❗✔️:       opts = options[curr];
    // RDKit❗✔️:       CHECK_INVARIANT(opts.size() > 0, "");
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       INT_DEQUE lstack;
    // RDKit❗✔️:       std::vector<int> optsV;
    // RDKit❗✔️:       std::vector<int> wedgedOptsV;
    // RDKit❗✔️:       std::vector<int> nbrs;
    // RDKit❗✔️:       for (auto nbrAtom : mol.atomNeighbors(mol.getAtomWithIdx(curr))) {
    // RDKit❗✔️:         const auto nbrIdx = static_cast<int>(nbrAtom->getIdx());
    // RDKit❗✔️:         if (!inAllAtms.test(nbrIdx)) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (std::find(done.begin(), done.end(), nbrIdx) != done.end()) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         nbrs.push_back(nbrIdx);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       std::sort(nbrs.begin(), nbrs.end(), lessByRank);
    // RDKit❗✔️:       for (int nbrIdx : nbrs) {
    // RDKit❗✔️:         auto nbrBond = mol.getBondBetweenAtoms(curr, nbrIdx);
    // RDKit❗✔️:         if (std::find(astack.begin(), astack.end(), nbrIdx) == astack.end()) {
    // RDKit❗✔️:           lstack.push_back(nbrIdx);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (cCand && dBndCands[nbrIdx] &&
    // RDKit❗✔️:             (nbrBond->getIsAromatic() ||
    // RDKit❗✔️:              mol.getAtomWithIdx(curr)->getAtomicNum() == 0 ||
    // RDKit❗✔️:              mol.getAtomWithIdx(nbrIdx)->getAtomicNum() == 0)) {
    // RDKit❗✔️:           if (nbrBond->getBondDir() == Bond::BondDir::BEGINWEDGE ||
    // RDKit❗✔️:               nbrBond->getBondDir() == Bond::BondDir::BEGINDASH) {
    // RDKit❗✔️:             wedgedOptsV.push_back(nbrIdx);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             optsV.push_back(nbrIdx);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (int v : optsV) {
    // RDKit❗✔️:         opts.push_back(v);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       for (int v : wedgedOptsV) {
    // RDKit❗✔️:         opts.push_back(v);
    // RDKit❗✔️:       }
    // RDKit❗✔️:       astack.insert(astack.end(), lstack.begin(), lstack.end());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (cCand) {
    // RDKit❗✔️:       if (!opts.empty()) {
    // RDKit❗✔️:         ncnd = opts.front();
    // RDKit❗✔️:         opts.pop_front();
    // RDKit❗✔️:         auto bnd = mol.getBondBetweenAtoms(curr, ncnd);
    // RDKit❗✔️:         bnd->setBondType(Bond::DOUBLE);
    // RDKit❗✔️:         if (bnd->getBondDir() != Bond::BondDir::NONE) {
    // RDKit❗✔️:           bnd->setBondDir(Bond::BondDir::NONE);
    // RDKit❗✔️:         }
    // RDKit❗✔️:         dBndCands[curr] = 0;
    // RDKit❗✔️:         dBndCands[ncnd] = 0;
    // RDKit❗✔️:         dBndAdds[bnd->getIdx()] = 1;
    // RDKit❗✔️:         localBondsAdded[bnd->getIdx()] = 1;
    // RDKit❗✔️:         if (options.find(curr) != options.end()) {
    // RDKit❗✔️:           if (opts.size() == 0) {
    // RDKit❗✔️:             options.erase(curr);
    // RDKit❗✔️:             btmoves.pop_back();
    // RDKit❗✔️:             if (btmoves.size() > 0) {
    // RDKit❗✔️:               lastOpt = btmoves.back();
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               lastOpt = -1;
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             options[curr] = opts;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           if (opts.size() > 0) {
    // RDKit❗✔️:             lastOpt = curr;
    // RDKit❗✔️:             btmoves.push_back(lastOpt);
    // RDKit❗✔️:             options[curr] = opts;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else if (mol.getAtomWithIdx(curr)->getAtomicNum()) {
    // RDKit❗✔️:         if ((lastOpt >= 0) && (numBT < maxBackTracks)) {
    // RDKit❗✔️:           backTrack(mol, options, lastOpt, done, astack, dBndCands, dBndAdds);
    // RDKit❗✔️:           ++numBT;
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           for (unsigned int bidx = 0; bidx < mol.getNumBonds(); ++bidx) {
    // RDKit❗✔️:             if (localBondsAdded[bidx]) {
    // RDKit❗✔️:               mol.getBondWithIdx(bidx)->setBondType(Bond::SINGLE);
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:           return false;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return true;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION kekulizeWorker
    let all_atom_set = state.all_atoms.iter().copied().collect::<BTreeSet<_>>();
    let mut adjacency = BTreeMap::<AtomId, Vec<(AtomId, BondId)>>::new();
    for &(bond, begin, end) in &state.aromatic_edges {
        if all_atom_set.contains(&begin) && all_atom_set.contains(&end) {
            adjacency.entry(begin).or_default().push((end, bond));
            adjacency.entry(end).or_default().push((begin, bond));
        }
    }
    for neighbors in adjacency.values_mut() {
        neighbors.sort_by_key(|(atom, _)| (atom_ranks[atom.index()], atom.index()));
    }

    let mut wedge_end_atoms = BTreeSet::new();
    for bond in molecule.bonds() {
        if matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) && all_atom_set.contains(&bond.end())
        {
            wedge_end_atoms.insert(bond.end());
        }
    }

    let mut sorted_atoms = state.all_atoms.clone();
    sorted_atoms.sort_by(|left, right| {
        match (
            wedge_end_atoms.contains(left),
            wedge_end_atoms.contains(right),
        ) {
            (true, false) => std::cmp::Ordering::Less,
            (false, true) => std::cmp::Ordering::Greater,
            _ => (atom_ranks[left.index()], left.index())
                .cmp(&(atom_ranks[right.index()], right.index())),
        }
    });

    let mut d_bnd_cands = vec![false; molecule.num_atoms()];
    for atom in &state.candidate_atoms {
        if !switched_off.contains(atom) {
            d_bnd_cands[atom.index()] = true;
        }
    }
    let mut d_bnd_adds = vec![false; molecule.num_bonds()];
    let mut local_bonds_added = vec![false; molecule.num_bonds()];
    let mut done = initial_done.to_vec();
    let mut astack = VecDeque::<AtomId>::new();
    let mut options = BTreeMap::<AtomId, VecDeque<AtomId>>::new();
    let mut last_opt = None::<AtomId>;
    let mut btmoves = Vec::<AtomId>::new();
    let mut matched = BTreeSet::<BondId>::new();
    let mut num_bt = 0usize;

    while done.len() < sorted_atoms.len() || !astack.is_empty() {
        let curr = if let Some(curr) = astack.pop_front() {
            curr
        } else {
            sorted_atoms
                .iter()
                .copied()
                .find(|atom| !done.contains(atom))
                .expect("starting point not found")
        };
        done.push(curr);
        let c_cand = d_bnd_cands[curr.index()];
        let mut opts = if let Some(saved) = options.get(&curr) {
            saved.clone()
        } else {
            let mut lstack = VecDeque::new();
            let mut opts_v = Vec::new();
            let mut wedged_opts_v = Vec::new();
            let mut nbrs = adjacency
                .get(&curr)
                .map(|neighbors| {
                    neighbors
                        .iter()
                        .filter_map(|(nbr, _)| (!done.contains(nbr)).then_some(*nbr))
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            nbrs.sort_by_key(|atom| (atom_ranks[atom.index()], atom.index()));
            for nbr_idx in nbrs {
                let nbr_bond = adjacency
                    .get(&curr)
                    .and_then(|neighbors| {
                        neighbors
                            .iter()
                            .find_map(|(nbr, bond)| (*nbr == nbr_idx).then_some(*bond))
                    })
                    .expect("bond between neighbor atoms not found");
                if !astack.contains(&nbr_idx) {
                    lstack.push_back(nbr_idx);
                }
                let bond = &molecule.bonds()[nbr_bond.index()];
                if c_cand
                    && d_bnd_cands[nbr_idx.index()]
                    && (bond.is_aromatic()
                        || molecule.atoms()[curr.index()].atomic_number() == 0
                        || molecule.atoms()[nbr_idx.index()].atomic_number() == 0)
                {
                    if matches!(
                        bond.direction(),
                        BondDirection::BeginWedge | BondDirection::BeginDash
                    ) {
                        wedged_opts_v.push(nbr_idx);
                    } else {
                        opts_v.push(nbr_idx);
                    }
                }
            }
            let mut computed = VecDeque::new();
            computed.extend(opts_v);
            computed.extend(wedged_opts_v);
            astack.extend(lstack);
            computed
        };
        if c_cand {
            if let Some(ncnd) = opts.pop_front() {
                let bond_id = adjacency
                    .get(&curr)
                    .and_then(|neighbors| {
                        neighbors
                            .iter()
                            .find_map(|(nbr, bond)| (*nbr == ncnd).then_some(*bond))
                    })
                    .expect("bond between current atom and option atom not found");
                d_bnd_cands[curr.index()] = false;
                d_bnd_cands[ncnd.index()] = false;
                d_bnd_adds[bond_id.index()] = true;
                local_bonds_added[bond_id.index()] = true;
                matched.insert(bond_id);
                if options.contains_key(&curr) {
                    if opts.is_empty() {
                        options.remove(&curr);
                        let _ = btmoves.pop();
                        last_opt = btmoves.last().copied();
                    } else {
                        options.insert(curr, opts);
                    }
                } else if !opts.is_empty() {
                    last_opt = Some(curr);
                    btmoves.push(curr);
                    options.insert(curr, opts);
                }
            } else if molecule.atoms()[curr.index()].atomic_number() != 0 {
                if let Some(last_opt_atom) = last_opt {
                    if num_bt < max_backtracks {
                        back_track(
                            molecule,
                            &mut options,
                            last_opt_atom,
                            &mut done,
                            &mut astack,
                            &mut d_bnd_cands,
                            &mut d_bnd_adds,
                            &mut matched,
                        );
                        num_bt += 1;
                    } else {
                        for (bond_idx, was_added) in local_bonds_added.iter().enumerate() {
                            if *was_added {
                                matched.remove(&BondId::new(bond_idx));
                            }
                        }
                        return None;
                    }
                } else {
                    for (bond_idx, was_added) in local_bonds_added.iter().enumerate() {
                        if *was_added {
                            matched.remove(&BondId::new(bond_idx));
                        }
                    }
                    return None;
                }
            }
        }
    }
    Some(matched)
}

fn permute_dummies_and_match(
    molecule: &Molecule,
    state: &KekulizeCandidateState,
    questions: &[AtomId],
    atom_ranks: &[usize],
    max_backtracks: usize,
) -> Option<BTreeSet<BondId>> {
    // BEGIN RDKIT CPP FUNCTION QuestionEnumerator
    // RDKit✔️✔️: class QuestionEnumerator {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   QuestionEnumerator(INT_VECT questions)
    // RDKit✔️✔️:       : d_questions(std::move(questions)), d_pos(1) {};
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
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION QuestionEnumerator
    // BEGIN RDKIT CPP FUNCTION permuteDummiesAndKekulize
    // RDKit❗✔️: bool permuteDummiesAndKekulize(RWMol &mol, const INT_VECT &allAtms,
    // RDKit❗✔️:                                boost::dynamic_bitset<> dBndCands,
    // RDKit❗✔️:                                INT_VECT &questions,
    // RDKit❗✔️:                                const UINT_VECT &atomRanks,
    // RDKit❗✔️:                                unsigned int maxBackTracks) {
    // RDKit❗✔️:   boost::dynamic_bitset<> atomsInPlay(mol.getNumAtoms());
    // RDKit❗✔️:   for (int allAtm : allAtms) {
    // RDKit❗✔️:     atomsInPlay[allAtm] = 1;
    // RDKit❗✔️:   }
    // RDKit✔️✔️:   bool kekulized = false;
    // RDKit✔️✔️:   QuestionEnumerator qEnum(questions);
    // RDKit✔️✔️:   while (!kekulized && questions.size()) {
    // RDKit❗✔️:     boost::dynamic_bitset<> dBndAdds(mol.getNumBonds());
    // RDKit✔️✔️:     INT_VECT done;
    // RDKit❗✔️:     // reset the state: all aromatic bonds are remarked to single:
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       if (bond->getIsAromatic() && bond->getBondType() != Bond::SINGLE &&
    // RDKit❗✔️:           atomsInPlay[bond->getBeginAtomIdx()] &&
    // RDKit❗✔️:           atomsInPlay[bond->getEndAtomIdx()]) {
    // RDKit❗✔️:         bond->setBondType(Bond::SINGLE);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit✔️✔️:     const auto &switchOff = qEnum.next();
    // RDKit✔️✔️:     if (!switchOff.size()) {
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     auto tCands = dBndCands;
    // RDKit✔️✔️:     for (int it : switchOff) {
    // RDKit✔️✔️:       tCands[it] = 0;
    // RDKit✔️✔️:     }
    // RDKit❗✔️:     kekulized =
    // RDKit❗✔️:         kekulizeWorker(mol, allAtms, tCands, dBndAdds, done, atomRanks,
    // RDKit❗✔️:                        maxBackTracks);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return kekulized;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION permuteDummiesAndKekulize
    if questions.is_empty() {
        return None;
    }
    let question_count = questions.len();
    for mask in 1usize..(1usize << question_count) {
        let switched_off = questions
            .iter()
            .enumerate()
            .filter_map(|(idx, atom)| ((mask & (1usize << idx)) != 0).then_some(*atom))
            .collect::<BTreeSet<_>>();
        if let Some(matched) = kekulize_worker_matching(
            molecule,
            state,
            &[],
            atom_ranks,
            max_backtracks,
            &switched_off,
        ) {
            return Some(matched);
        }
    }
    None
}

fn back_track(
    molecule: &Molecule,
    options: &mut BTreeMap<AtomId, VecDeque<AtomId>>,
    last_opt: AtomId,
    done: &mut Vec<AtomId>,
    aqueue: &mut VecDeque<AtomId>,
    d_bnd_cands: &mut [bool],
    d_bnd_adds: &mut [bool],
    matched: &mut BTreeSet<BondId>,
) {
    // BEGIN RDKIT CPP FUNCTION backTrack
    // RDKit✔️✔️: void backTrack(RWMol &mol, INT_INT_DEQ_MAP &, int lastOpt, INT_VECT &done,
    // RDKit✔️✔️:                INT_DEQUE &aqueue, boost::dynamic_bitset<> &dBndCands,
    // RDKit✔️✔️:                boost::dynamic_bitset<> &dBndAdds) {
    // RDKit✔️✔️:   auto ei = std::find(done.begin(), done.end(), lastOpt);
    // RDKit✔️✔️:   INT_VECT tdone;
    // RDKit✔️✔️:   tdone.insert(tdone.end(), done.begin(), ei);
    let split = done
        .iter()
        .position(|atom| *atom == last_opt)
        .expect("lastOpt must exist in done");
    let tdone = done[..split].to_vec();
    // RDKit✔️✔️:   INT_VECT_CRI eri = std::find(done.rbegin(), done.rend(), lastOpt);
    // RDKit✔️✔️:   ++eri;
    // RDKit✔️✔️:   for (INT_VECT_CRI ri = done.rbegin(); ri != eri; ++ri) {
    // RDKit✔️✔️:     aqueue.push_front(*ri);
    // RDKit✔️✔️:   }
    for atom in done[split..].iter().rev() {
        aqueue.push_front(*atom);
    }
    // RDKit✔️✔️:   for (unsigned int bi = 0; bi < nbnds; ++bi) {
    for bond in molecule.bonds() {
        let bi = bond.id().index();
        if d_bnd_adds[bi] {
            let aid1 = bond.begin();
            let aid2 = bond.end();
            if !tdone.contains(&aid1) && !tdone.contains(&aid2) {
                d_bnd_adds[bi] = false;
                matched.remove(&bond.id());
                d_bnd_cands[aid1.index()] = true;
                d_bnd_cands[aid2.index()] = true;
            }
        }
    }
    *done = tdone;
    options.retain(|atom, _| done.contains(atom));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION backTrack
}

fn rdkit_kekulize_default_valence(atomic_number: u8) -> Result<i32, ValenceError> {
    crate::valence::rdkit_default_valence(atomic_number)
}

// Ported RDKit entrypoint kept private until the public operation surface has a
// deliberate `KekulizeIfPossible` API. Tests exercise it directly.
#[allow(dead_code)]
fn kekulize_if_possible_assignment(
    molecule: &Molecule,
    clear_aromatic_flags: bool,
    canonical: bool,
    max_backtracks: usize,
) -> Result<Option<KekulizeAssignment>, KekulizeError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::KekulizeIfPossible
    // RDKit❗✔️: bool KekulizeIfPossible(RWMol &mol, bool markAtomsBonds, bool canonical,
    // RDKit❗✔️:                         unsigned int maxBackTracks) {
    // RDKit❗✔️:   boost::dynamic_bitset<> aromaticBonds(mol.getNumBonds());
    // RDKit❗✔️:   for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:     if (bond->getIsAromatic()) {
    // RDKit❗✔️:       aromaticBonds.set(bond->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   boost::dynamic_bitset<> aromaticAtoms(mol.getNumAtoms());
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (isAromaticAtom(*atom)) {
    // RDKit❗✔️:       aromaticAtoms.set(atom->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool res = true;
    // RDKit❗✔️:   try {
    // RDKit❗✔️:     Kekulize(mol, markAtomsBonds, canonical, maxBackTracks);
    // RDKit❗✔️:   } catch (const MolSanitizeException &) {
    // RDKit❗✔️:     res = false;
    // RDKit❗✔️:     for (unsigned int i = 0; i < mol.getNumBonds(); ++i) {
    // RDKit❗✔️:       if (aromaticBonds[i]) {
    // RDKit❗✔️:         auto bond = mol.getBondWithIdx(i);
    // RDKit❗✔️:         bond->setIsAromatic(true);
    // RDKit❗✔️:         bond->setBondType(Bond::BondType::AROMATIC);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit❗✔️:       if (aromaticAtoms[i]) {
    // RDKit❗✔️:         mol.getAtomWithIdx(i)->setIsAromatic(true);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION MolOps::KekulizeIfPossible
    match kekulize_assignment(
        molecule,
        None,
        clear_aromatic_flags,
        canonical,
        max_backtracks,
    ) {
        Ok(assignment) => Ok(Some(assignment)),
        Err(KekulizeError::ProtocolDebt { branch, reason }) => {
            Err(KekulizeError::ProtocolDebt { branch, reason })
        }
        Err(
            KekulizeError::UnkekulizedAtoms { .. }
            | KekulizeError::NonRingAromaticAtom { .. }
            | KekulizeError::ValenceChanged { .. }
            | KekulizeError::FragmentBitsetSizeMismatch { .. }
            | KekulizeError::CanonicalRankSymbolSizeMismatch { .. }
            | KekulizeError::Valence(_)
            | KekulizeError::RingFinding(_),
        ) => Ok(None),
    }
}

fn rdkit_is_early_atom(atomic_number: u8) -> bool {
    // BEGIN RDKIT CPP FUNCTION isEarlyAtom
    // RDKit✔️✔️: bool isEarlyAtom(int atomicNum) {
    // RDKit✔️✔️:   static const bool table[119] = {
    // RDKit✔️✔️:       false,  // #0 *
    // RDKit✔️✔️:       false,  // #1 H
    // RDKit✔️✔️:       false,  // #2 He
    // RDKit✔️✔️:       true,   // #3 Li
    // RDKit✔️✔️:       true,   // #4 Be
    // RDKit✔️✔️:       true,   // #5 B
    // RDKit✔️✔️:       false,  // #6 C
    // RDKit✔️✔️:       false,  // #7 N
    // RDKit✔️✔️:       false,  // #8 O
    // RDKit✔️✔️:       false,  // #9 F
    // RDKit✔️✔️:       false,  // #10 Ne
    // RDKit✔️✔️:       true,   // #11 Na
    // RDKit✔️✔️:       true,   // #12 Mg
    // RDKit✔️✔️:       true,   // #13 Al
    // RDKit✔️✔️:       false,  // #14 Si
    // RDKit✔️✔️:       false,  // #15 P
    // RDKit✔️✔️:       false,  // #16 S
    // RDKit✔️✔️:       false,  // #17 Cl
    // RDKit✔️✔️:       false,  // #18 Ar
    // RDKit✔️✔️:       true,   // #19 K
    // RDKit✔️✔️:       true,   // #20 Ca
    // RDKit✔️✔️:       true,   // #21 Sc
    // RDKit✔️✔️:       true,   // #22 Ti
    // RDKit✔️✔️:       false,  // #23 V
    // RDKit✔️✔️:       false,  // #24 Cr
    // RDKit✔️✔️:       false,  // #25 Mn
    // RDKit✔️✔️:       false,  // #26 Fe
    // RDKit✔️✔️:       false,  // #27 Co
    // RDKit✔️✔️:       false,  // #28 Ni
    // RDKit✔️✔️:       false,  // #29 Cu
    // RDKit✔️✔️:       true,   // #30 Zn
    // RDKit✔️✔️:       true,   // #31 Ga
    // RDKit✔️✔️:       true,   // #32 Ge  see github #2606
    // RDKit✔️✔️:       false,  // #33 As
    // RDKit✔️✔️:       false,  // #34 Se
    // RDKit✔️✔️:       false,  // #35 Br
    // RDKit✔️✔️:       false,  // #36 Kr
    // RDKit✔️✔️:       true,   // #37 Rb
    // RDKit✔️✔️:       true,   // #38 Sr
    // RDKit✔️✔️:       true,   // #39 Y
    // RDKit✔️✔️:       true,   // #40 Zr
    // RDKit✔️✔️:       true,   // #41 Nb
    // RDKit✔️✔️:       false,  // #42 Mo
    // RDKit✔️✔️:       false,  // #43 Tc
    // RDKit✔️✔️:       false,  // #44 Ru
    // RDKit✔️✔️:       false,  // #45 Rh
    // RDKit✔️✔️:       false,  // #46 Pd
    // RDKit✔️✔️:       false,  // #47 Ag
    // RDKit✔️✔️:       true,   // #48 Cd
    // RDKit✔️✔️:       true,   // #49 In
    // RDKit✔️✔️:       true,   // #50 Sn  see github #2606
    // RDKit✔️✔️:       true,   // #51 Sb  see github #2775
    // RDKit✔️✔️:       false,  // #52 Te
    // RDKit✔️✔️:       false,  // #53 I
    // RDKit✔️✔️:       false,  // #54 Xe
    // RDKit✔️✔️:       true,   // #55 Cs
    // RDKit✔️✔️:       true,   // #56 Ba
    // RDKit✔️✔️:       true,   // #57 La
    // RDKit✔️✔️:       true,   // #58 Ce
    // RDKit✔️✔️:       true,   // #59 Pr
    // RDKit✔️✔️:       true,   // #60 Nd
    // RDKit✔️✔️:       true,   // #61 Pm
    // RDKit✔️✔️:       false,  // #62 Sm
    // RDKit✔️✔️:       false,  // #63 Eu
    // RDKit✔️✔️:       false,  // #64 Gd
    // RDKit✔️✔️:       false,  // #65 Tb
    // RDKit✔️✔️:       false,  // #66 Dy
    // RDKit✔️✔️:       false,  // #67 Ho
    // RDKit✔️✔️:       false,  // #68 Er
    // RDKit✔️✔️:       false,  // #69 Tm
    // RDKit✔️✔️:       false,  // #70 Yb
    // RDKit✔️✔️:       false,  // #71 Lu
    // RDKit✔️✔️:       true,   // #72 Hf
    // RDKit✔️✔️:       true,   // #73 Ta
    // RDKit✔️✔️:       false,  // #74 W
    // RDKit✔️✔️:       false,  // #75 Re
    // RDKit✔️✔️:       false,  // #76 Os
    // RDKit✔️✔️:       false,  // #77 Ir
    // RDKit✔️✔️:       false,  // #78 Pt
    // RDKit✔️✔️:       false,  // #79 Au
    // RDKit✔️✔️:       true,   // #80 Hg
    // RDKit✔️✔️:       true,   // #81 Tl
    // RDKit✔️✔️:       true,   // #82 Pb  see github #2606
    // RDKit✔️✔️:       true,   // #83 Bi  see github #2775
    // RDKit✔️✔️:       false,  // #84 Po
    // RDKit✔️✔️:       false,  // #85 At
    // RDKit✔️✔️:       false,  // #86 Rn
    // RDKit✔️✔️:       true,   // #87 Fr
    // RDKit✔️✔️:       true,   // #88 Ra
    // RDKit✔️✔️:       true,   // #89 Ac
    // RDKit✔️✔️:       true,   // #90 Th
    // RDKit✔️✔️:       true,   // #91 Pa
    // RDKit✔️✔️:       true,   // #92 U
    // RDKit✔️✔️:       true,   // #93 Np
    // RDKit✔️✔️:       false,  // #94 Pu
    // RDKit✔️✔️:       false,  // #95 Am
    // RDKit✔️✔️:       false,  // #96 Cm
    // RDKit✔️✔️:       false,  // #97 Bk
    // RDKit✔️✔️:       false,  // #98 Cf
    // RDKit✔️✔️:       false,  // #99 Es
    // RDKit✔️✔️:       false,  // #100 Fm
    // RDKit✔️✔️:       false,  // #101 Md
    // RDKit✔️✔️:       false,  // #102 No
    // RDKit✔️✔️:       false,  // #103 Lr
    // RDKit✔️✔️:       true,   // #104 Rf
    // RDKit✔️✔️:       true,   // #105 Db
    // RDKit✔️✔️:       true,   // #106 Sg
    // RDKit✔️✔️:       true,   // #107 Bh
    // RDKit✔️✔️:       true,   // #108 Hs
    // RDKit✔️✔️:       true,   // #109 Mt
    // RDKit✔️✔️:       true,   // #110 Ds
    // RDKit✔️✔️:       true,   // #111 Rg
    // RDKit✔️✔️:       true,   // #112 Cn
    // RDKit✔️✔️:       true,   // #113 Nh
    // RDKit✔️✔️:       true,   // #114 Fl
    // RDKit✔️✔️:       true,   // #115 Mc
    // RDKit✔️✔️:       true,   // #116 Lv
    // RDKit✔️✔️:       true,   // #117 Ts
    // RDKit✔️✔️:       true,   // #118 Og
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   return ((unsigned int)atomicNum < 119) && table[atomicNum];
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isEarlyAtom
    matches!(
        atomic_number,
        3 | 4
            | 5
            | 11
            | 12
            | 13
            | 19
            | 20
            | 21
            | 22
            | 30
            | 31
            | 32
            | 37
            | 38
            | 39
            | 40
            | 41
            | 48
            | 49
            | 50
            | 51
            | 55
            | 56
            | 57
            | 58
            | 59
            | 60
            | 61
            | 72
            | 73
            | 80
            | 81
            | 82
            | 83
            | 87
            | 88
            | 89
            | 90
            | 91
            | 92
            | 93
            | 104
            | 105
            | 106
            | 107
            | 108
            | 109
            | 110
            | 111
            | 112
            | 113
            | 114
            | 115
            | 116
            | 117
            | 118
    )
}

fn ordered_kekulize_rings(
    molecule: &Molecule,
    rings: &RingInfo,
    atoms_to_use: &[bool],
    bonds_to_use: &[bool],
    wedged_atoms: &BTreeSet<AtomId>,
) -> Vec<KekulizeRing> {
    let mut front = Vec::new();
    let mut back = Vec::new();
    for ring_idx in 0..rings.num_rings() {
        let atoms = &rings.atom_rings()[ring_idx];
        let bonds = &rings.bond_rings()[ring_idx];
        let contains_non_dummy = atoms
            .iter()
            .all(|atom| atom.index() < molecule.num_atoms() && atoms_to_use[atom.index()])
            && atoms
                .iter()
                .any(|atom| molecule.atoms()[atom.index()].atomic_number() != 0);
        if !contains_non_dummy {
            continue;
        }
        if !bonds.iter().all(|bond| bonds_to_use[bond.index()]) {
            continue;
        }
        let start_pos = atoms
            .iter()
            .position(|atom| wedged_atoms.contains(atom))
            .unwrap_or(0);
        let rotated_atoms = (0..atoms.len())
            .map(|offset| atoms[(offset + start_pos) % atoms.len()])
            .collect::<Vec<_>>();
        let entry = KekulizeRing {
            atoms: rotated_atoms,
            bonds: bonds.clone(),
            source_ring: ring_idx,
        };
        if start_pos == 0 && !atoms.iter().any(|atom| wedged_atoms.contains(atom)) {
            back.push(entry);
        } else {
            front.push(entry);
        }
    }
    front.extend(back);
    front
}

fn fused_ring_systems(rings: &[KekulizeRing]) -> Vec<Vec<usize>> {
    let mut systems = Vec::<Vec<usize>>::new();
    let mut seen = vec![false; rings.len()];
    for start in 0..rings.len() {
        if seen[start] {
            continue;
        }
        let mut stack = vec![start];
        let mut system = Vec::new();
        seen[start] = true;
        while let Some(ring) = stack.pop() {
            system.push(ring);
            for next in 0..rings.len() {
                if !seen[next] && rings_share_bond(rings, ring, next) {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        systems.push(system);
    }
    systems
}

fn rings_share_bond(rings: &[KekulizeRing], left: usize, right: usize) -> bool {
    rings[left]
        .bonds
        .iter()
        .any(|bond| rings[right].bonds.contains(bond))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        AtomSpec, BondOrder, BondQueryPredicate, BondSpec, Element, MoleculeBuilder, QueryNode,
    };

    #[test]
    fn kekulize_assignment_returns_empty_for_empty_molecule() {
        let molecule = Molecule::new();

        let assignment = kekulize_assignment(&molecule, None, true, false, 100).unwrap();

        assert!(assignment.is_empty());
    }

    #[test]
    fn kekulize_assignment_kekulizes_benzene_like_aromatic_cycle() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let assignment = kekulize_assignment(&molecule, None, true, false, 100).unwrap();
        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        let single_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Single))
            .count();

        assert_eq!(double_bonds, 3);
        assert_eq!(single_bonds, 3);
        assert!(
            molecule
                .bonds()
                .iter()
                .all(|bond| assignment.bond_aromatic_flag(bond.id()) == Some(false))
        );
        assert!(
            molecule
                .atoms()
                .iter()
                .all(|atom| assignment.atom_aromatic_flag(atom.id()) == Some(false))
        );
    }

    #[test]
    fn kekulize_fragment_returns_empty_when_atoms_to_use_is_empty_like_rdkit() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();
        let rings = symmetrize_sssr(&molecule).unwrap();

        let assignment = kekulize_fragment_assignment(
            &molecule,
            &rings,
            &[false; 6],
            vec![true; molecule.num_bonds()],
            true,
            false,
            100,
        )
        .unwrap();

        assert!(assignment.is_empty());
    }

    #[test]
    fn kekulize_fragment_rejects_bitset_size_mismatch() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();
        let rings = symmetrize_sssr(&molecule).unwrap();

        let error = kekulize_fragment_assignment(
            &molecule,
            &rings,
            &[true; 5],
            vec![true; molecule.num_bonds()],
            true,
            false,
            100,
        )
        .unwrap_err();

        assert!(matches!(
            error,
            KekulizeError::FragmentBitsetSizeMismatch { atoms: 5, bonds: 6 }
        ));
    }

    #[test]
    fn kekulize_assignment_runs_canonical_fragment_ranking_for_plain_aromatic_cycle() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let assignment = kekulize_assignment(&molecule, None, true, true, 100).unwrap();

        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        let single_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Single))
            .count();
        assert_eq!(double_bonds, 3);
        assert_eq!(single_bonds, 3);
    }

    #[test]
    fn kekulize_if_possible_returns_assignment_for_kekulizable_molecule() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let assignment = kekulize_if_possible_assignment(&molecule, true, false, 100)
            .unwrap()
            .unwrap();

        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        assert_eq!(double_bonds, 3);
    }

    #[test]
    fn kekulize_if_possible_runs_canonical_fragment_ranking_for_plain_aromatic_cycle() {
        let molecule = Molecule::from_smiles_with_sanitize("c1ccccc1", false).unwrap();

        let assignment = kekulize_if_possible_assignment(&molecule, true, true, 100)
            .unwrap()
            .unwrap();

        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        assert_eq!(double_bonds, 3);
    }

    #[test]
    fn canonical_kekulize_handles_stereo_ranking_branch() {
        let mut builder = MoleculeBuilder::new();
        let atoms = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true)))
            .collect::<Vec<_>>();
        for idx in 0..6 {
            builder
                .add_bond(
                    BondSpec::new(atoms[idx], atoms[(idx + 1) % 6], BondOrder::Aromatic)
                        .with_aromatic(true),
                )
                .unwrap();
        }
        let mut molecule = builder.build().unwrap();
        molecule.topology_block_mut().atoms[0].set_chiral_tag(crate::ChiralTag::TetrahedralCw);

        let assignment = kekulize_assignment(&molecule, None, true, true, 100).unwrap();

        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        assert_eq!(double_bonds, 3);
    }

    #[test]
    fn kekulize_assignment_handles_dummy_atom_permutation_branch() {
        let molecule = Molecule::from_smiles_with_sanitize("*1cccc1", false).unwrap();

        let assignment = kekulize_assignment(&molecule, None, true, false, 100).unwrap();

        let double_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Double))
            .count();
        let single_bonds = molecule
            .bonds()
            .iter()
            .filter(|bond| assignment.bond_order(bond.id()) == Some(BondOrder::Single))
            .count();

        assert_eq!(double_bonds, 2, "{:?}", assignment.bond_orders);
        assert_eq!(single_bonds, 1, "{:?}", assignment.bond_orders);
    }

    #[test]
    fn kekulize_assignment_excludes_bond_type_query_from_bonds_to_use() {
        let mut builder = MoleculeBuilder::new();
        let a0 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let a1 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        let a2 = builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true));
        builder
            .add_bond(
                BondSpec::new(a0, a1, BondOrder::Aromatic)
                    .with_aromatic(true)
                    .with_query(QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                        BondOrder::Single,
                        BondOrder::Aromatic,
                    ]))),
            )
            .unwrap();
        builder
            .add_bond(BondSpec::new(a1, a2, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        builder
            .add_bond(BondSpec::new(a2, a0, BondOrder::Aromatic).with_aromatic(true))
            .unwrap();
        let molecule = builder.build().unwrap();

        let assignment = kekulize_assignment(&molecule, None, true, false, 100).unwrap();

        assert_eq!(assignment.bond_order(BondId::new(0)), None);
        assert_eq!(assignment.bond_aromatic_flag(BondId::new(0)), None);
        assert_eq!(assignment.bond_aromatic_flag(BondId::new(1)), Some(false));
        assert_eq!(assignment.bond_aromatic_flag(BondId::new(2)), Some(false));
    }

    #[test]
    fn kekulize_assignment_keeps_non_bond_type_query_in_bonds_to_use() {
        let mut builder = MoleculeBuilder::new();
        let atoms = (0..6)
            .map(|_| builder.add_atom(AtomSpec::new(Element::C).with_aromatic(true)))
            .collect::<Vec<_>>();
        builder
            .add_bond(
                BondSpec::new(atoms[0], atoms[1], BondOrder::Aromatic)
                    .with_aromatic(true)
                    .with_query(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
            )
            .unwrap();
        for idx in 1..6 {
            builder
                .add_bond(
                    BondSpec::new(atoms[idx], atoms[(idx + 1) % 6], BondOrder::Aromatic)
                        .with_aromatic(true),
                )
                .unwrap();
        }
        let molecule = builder.build().unwrap();

        let assignment = kekulize_assignment(&molecule, None, true, false, 100).unwrap();

        assert_eq!(assignment.bond_aromatic_flag(BondId::new(0)), Some(false));
    }

    #[test]
    fn kekulize_source_markers_do_not_hide_elided_rdkit_lines() {
        let source = include_str!("kekulize.rs");

        for line in source.lines().filter(|line| line.contains("RDKit")) {
            assert!(
                !line.contains("..."),
                "RDKit source marker must not contain elided source: {line}"
            );
        }
    }
}
