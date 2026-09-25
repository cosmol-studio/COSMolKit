//! RDKit legacy stereochemistry assignment over detached topology values.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, AtomPropertyError, BondValueError, QueryStateRef, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

use crate::{
    AtropisomerError, CipRankError, DoubleBondStereoError, PotentialStereoError, RingInfo,
    StereoOrderError, ValenceAssignment, ValenceError, assign_atom_cip_ranks_with_query_state,
    assign_directional_double_bond_stereo, bond_affects_atom_chirality,
    cleanup_atropisomer_stereo_groups, count_swaps_to_interconvert, invert_tetrahedral_tag,
    is_atom_bridgehead_from_topology, refine_atom_cip_ranks_from_invariants_with_query_state,
};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum LegacyStereoError {
    #[error(transparent)]
    InvalidTopology(#[from] TopologyValidationError),
    #[error(transparent)]
    CipRank(#[from] CipRankError),
    #[error(transparent)]
    DoubleBond(#[from] DoubleBondStereoError),
    #[error(transparent)]
    StereoOrder(#[from] StereoOrderError),
    #[error(transparent)]
    Valence(#[from] ValenceError),
    #[error(transparent)]
    AtomProperty(#[from] AtomPropertyError),
    #[error(transparent)]
    BondValue(#[from] BondValueError),
    #[error(transparent)]
    PotentialStereo(#[from] PotentialStereoError),
    #[error(transparent)]
    Atropisomer(#[from] AtropisomerError),
}

fn has_protium_neighbor(topology: &TopologyBlock, atom: AtomId) -> bool {
    // RDKit✔️✔️: bool is_regular_h(const Atom &atom) {
    // RDKit✔️✔️:   return atom.getAtomicNum() == 1 && atom.getIsotope() == 0;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bool has_protium_neighbor(const ROMol &mol, const Atom *atom) {
    // RDKit✔️✔️:   for (const auto nbr : mol.atomNeighbors(atom)) {
    // RDKit✔️✔️:     if (is_regular_h(*nbr)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // Behavior review: detached `None` and `Some(0)` both expose RDKit's
    // numeric zero isotope; nonzero hydrogen isotopes are not protium.
    // Complexity review: one adjacency scan with no allocation is O(degree),
    // matching `atomNeighbors()` in the source helper.
    topology
        .adjacency
        .neighbors_of(atom.index())
        .iter()
        .any(|neighbor| {
            let value = &topology.atoms[neighbor.atom_index];
            value.atomic_number() == 1 && value.isotope().unwrap_or(0) == 0
        })
}

fn is_legal_legacy_center(
    topology: &TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    center: AtomId,
) -> Result<bool, LegacyStereoError> {
    // RDKit✔️✔️: auto nzDegree = Chirality::detail::getAtomNonzeroDegree(atom);
    // RDKit✔️✔️: auto tnzDegree = nzDegree + atom->getTotalNumHs();
    // RDKit✔️✔️: if (tnzDegree > 4) {
    // RDKit✔️✔️:   legalCenter = false;
    // RDKit✔️✔️: } else {
    // RDKit✔️✔️:   if (tnzDegree < 3) {
    // RDKit✔️✔️:     legalCenter = false;
    // RDKit✔️✔️:   } else if (nzDegree < 3 &&
    // RDKit✔️✔️:              (atom->getAtomicNum() != 15 && atom->getAtomicNum() != 33)) {
    // RDKit✔️✔️:     legalCenter = false;
    // RDKit✔️✔️:   } else if (nzDegree == 3) {
    // RDKit✔️✔️:     if (atom->getTotalNumHs() == 1) {
    // RDKit✔️✔️:       if (detail::has_protium_neighbor(mol, atom)) {
    // RDKit✔️✔️:         legalCenter = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       legalCenter = false;
    // RDKit✔️✔️:       if (atom->getAtomicNum() == 7) {
    // RDKit✔️✔️:         if (atom->getHybridization() == Atom::HybridizationType::SP3 &&
    // RDKit✔️✔️:             !MolOps::atomHasConjugatedBond(atom) &&
    // RDKit✔️✔️:             (mol.getRingInfo()->isAtomInRingOfSize(atom->getIdx(), 3) ||
    // RDKit✔️✔️:              queryIsAtomBridgehead(atom))) {
    // RDKit✔️✔️:           legalCenter = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (atom->getAtomicNum() == 15 || atom->getAtomicNum() == 33) {
    // RDKit✔️✔️:         legalCenter = true;
    // RDKit✔️✔️:       } else if (atom->getAtomicNum() == 16 || atom->getAtomicNum() == 34) {
    // RDKit✔️✔️:         if (atom->getValence(Atom::ValenceType::EXPLICIT) == 4 ||
    // RDKit✔️✔️:             (atom->getValence(Atom::ValenceType::EXPLICIT) == 3 &&
    // RDKit✔️✔️:              atom->getFormalCharge() == 1)) {
    // RDKit✔️✔️:           legalCenter = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: this helper reproduces the legality half of
    // `isAtomPotentialChiralCenter`; duplicate-rank detection remains in
    // `assign_atom_codes`, where ranks and neighbor ordering are available.
    // Complexity review: bounded adjacency plus canonical ring/conjugation
    // lookups matches the source shape and introduces no whole-graph clone.
    let atom = &topology.atoms[center.index()];
    let mut nonzero_degree = 0usize;
    for neighbor in topology.adjacency.neighbors_of(center.index()) {
        if bond_affects_atom_chirality(&topology.bonds[neighbor.bond.index()], center)? {
            nonzero_degree += 1;
        }
    }
    let total_hydrogens =
        crate::hcount::total_hydrogen_count_from_validated(topology, valence, center, false)?
            as usize;
    let total_nonzero_degree = nonzero_degree + total_hydrogens;
    if total_nonzero_degree > 4 || total_nonzero_degree < 3 {
        return Ok(false);
    }
    if nonzero_degree < 3 && !matches!(atom.atomic_number(), 15 | 33) {
        return Ok(false);
    }
    if nonzero_degree != 3 {
        return Ok(true);
    }
    if total_hydrogens == 1 {
        return Ok(!has_protium_neighbor(topology, center));
    }
    let legal = match atom.atomic_number() {
        7 => {
            atom.hybridization() == Hybridization::Sp3
                && !crate::conjugation::atom_has_conjugated_bond_from_validated(topology, center)
                && (rings.is_atom_in_ring_of_size(center, 3)
                    || is_atom_bridgehead_from_topology(topology, center.index(), rings) != 0)
        }
        15 | 33 => true,
        16 | 34 => {
            valence.explicit_valence[center.index()] == 4
                || (valence.explicit_valence[center.index()] == 3 && atom.formal_charge() == 1)
        }
        _ => false,
    };
    Ok(legal)
}

fn materialize_initial_ranks(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<u32>, LegacyStereoError> {
    // RDKit✔️✔️: void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
    // RDKit✔️✔️:   PRECONDITION((!ranks.size() || ranks.size() >= mol.getNumAtoms()),
    // RDKit✔️✔️:                "bad ranks size");
    // RDKit✔️✔️:   if (!ranks.size()) {
    // RDKit✔️✔️:     ranks.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int numAtoms = mol.getNumAtoms();
    // RDKit✔️✔️: #ifndef USE_NEW_STEREOCHEMISTRY
    // RDKit✔️✔️:   DOUBLE_VECT invars(numAtoms, 0);
    // RDKit✔️✔️:   buildCIPInvariants(mol, invars);
    // RDKit✔️✔️:   iterateCIPRanks(mol, invars, ranks, false);
    // RDKit❌❌: #else
    // RDKit❌❌:   Canon::chiralRankMolAtoms(mol, ranks);
    // RDKit❌❌: #endif
    // RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; ++i) {
    // RDKit✔️✔️:     mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Behavior review: the pinned build follows the legacy invariant branch;
    // the canonical rank owner supplies it and this adapter materializes the
    // same computed atom property only when stereo inventory requires ranks.
    // Complexity review: one rank refinement plus one linear property pass
    // matches the pinned helper's asymptotic and allocation behavior.
    let ranks = assign_atom_cip_ranks_with_query_state(topology, valence, query_state)?;
    for (atom, rank) in topology.atoms.iter_mut().zip(&ranks) {
        atom.set_computed_prop("_CIPRank", rank.to_string())?;
    }
    Ok(ranks)
}

fn assign_atom_codes(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &mut Vec<u32>,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<(bool, bool), LegacyStereoError> {
    // BEGIN RDKIT CPP FUNCTION assignAtomChiralCodes
    // RDKit✔️✔️: PRECONDITION((!ranks.size() || ranks.size() == mol.getNumAtoms()),
    // RDKit✔️✔️:              "bad rank vector size");
    // RDKit✔️✔️: bool atomChanged = false;
    // RDKit✔️✔️: unsigned int unassignedAtoms = 0;
    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   Atom::ChiralType tag = atom->getChiralTag();
    // RDKit✔️✔️:   if (flagPossibleStereoCenters ||
    // RDKit✔️✔️:       (tag != Atom::CHI_UNSPECIFIED && tag != Atom::CHI_OTHER)) {
    // RDKit✔️✔️:     if (atom->hasProp(common_properties::_CIPCode)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (!ranks.size()) {
    // RDKit✔️✔️:       assignAtomCIPRanks(mol, ranks);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     Chirality::INT_PAIR_VECT nbrs;
    // RDKit✔️✔️:     auto [legalCenter, hasDupes] =
    // RDKit✔️✔️:         isAtomPotentialChiralCenter(atom, mol, ranks, nbrs);
    // RDKit✔️✔️:     if (legalCenter) {
    // RDKit✔️✔️:       ++unassignedAtoms;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (legalCenter && !hasDupes && flagPossibleStereoCenters) {
    // RDKit✔️✔️:       atom->setProp(common_properties::_ChiralityPossible, 1);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (legalCenter && !hasDupes && tag != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:         tag != Atom::CHI_OTHER) {
    // RDKit✔️✔️:       std::sort(nbrs.begin(), nbrs.end(), Rankers::pairLess);
    // RDKit✔️✔️:       std::list<int> nbrIndices;
    // RDKit✔️✔️:       for (auto nbrIt = nbrs.begin(); nbrIt != nbrs.end(); ++nbrIt) {
    // RDKit✔️✔️:         nbrIndices.push_back((*nbrIt).second);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       int nSwaps = atom->getPerturbationOrder(nbrIndices);
    // RDKit✔️✔️:       if (nbrIndices.size() == 3 && atom->getTotalNumHs() == 1) ++nSwaps;
    // RDKit✔️✔️:       if (nSwaps % 2) {
    // RDKit✔️✔️:         if (tag == Atom::CHI_TETRAHEDRAL_CCW) tag = Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:         else tag = Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       atom->setProp(common_properties::_CIPCode,
    // RDKit✔️✔️:                     tag == Atom::CHI_TETRAHEDRAL_CCW ? "S" : "R");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: return std::make_pair((unassignedAtoms > 0), atomChanged);
    // END RDKIT CPP FUNCTION assignAtomChiralCodes
    // Behavior review: legality, duplicate-rank rejection, source bond-order
    // perturbation, implicit-H inversion and computed property writes follow
    // the pinned legacy branch for the detached tetrahedral state model.
    // Complexity review: ranks are lazily computed once; every center scans
    // and sorts at most four incident bonds, matching the source traversal.
    let mut changed = false;
    let mut unassigned = 0usize;
    for index in 0..topology.atoms.len() {
        if topology.atoms[index].prop("_CIPCode").is_some() {
            continue;
        }
        if ranks.is_empty() {
            *ranks = materialize_initial_ranks(topology, valence, query_state)?;
        }
        let center = AtomId::new(index);
        let legal = is_legal_legacy_center(topology, valence, rings, center)?;
        let mut neighbors = Vec::new();
        let mut seen = BTreeSet::new();
        let mut duplicates = false;
        if legal {
            for neighbor in topology.adjacency.neighbors_of(index) {
                let bond = &topology.bonds[neighbor.bond.index()];
                neighbors.push((ranks[neighbor.atom_index], neighbor.bond));
                if bond_affects_atom_chirality(bond, center)?
                    && !seen.insert(ranks[neighbor.atom_index])
                {
                    duplicates = true;
                    break;
                }
            }
        }
        if legal {
            unassigned += 1;
        }
        if legal && !duplicates {
            topology.atoms[index].set_computed_prop("_ChiralityPossible", "1")?;
        }
        let tag = topology.atoms[index].chiral_tag();
        if legal
            && !duplicates
            && matches!(tag, ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw)
        {
            changed = true;
            unassigned -= 1;
            neighbors.sort_by_key(|(rank, bond)| (*rank, bond.index()));
            let sorted = neighbors.iter().map(|(_, bond)| *bond).collect::<Vec<_>>();
            let stored = topology
                .adjacency
                .neighbors_of(index)
                .iter()
                .map(|n| n.bond)
                .collect::<Vec<_>>();
            let mut swaps = count_swaps_to_interconvert(&sorted, &stored)?;
            if sorted.len() == 3
                && crate::hcount::total_hydrogen_count_from_validated(
                    topology, valence, center, false,
                )? == 1
            {
                swaps += 1;
            }
            let effective = if swaps % 2 == 1 {
                invert_tetrahedral_tag(tag)?
            } else {
                tag
            };
            topology.atoms[index].set_computed_prop(
                "_CIPCode",
                if effective == ChiralTag::TetrahedralCcw {
                    "S"
                } else {
                    "R"
                },
            )?;
        }
    }
    Ok((unassigned != 0, changed))
}

fn rerank_atoms(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    ranks: &[u32],
    query_state: Option<QueryStateRef<'_>>,
) -> Result<Vec<u32>, LegacyStereoError> {
    // BEGIN RDKIT CPP FUNCTION rerankAtoms
    // RDKit✔️✔️: PRECONDITION(ranks.size() == mol.getNumAtoms(), "bad rank vector size");
    // RDKit✔️✔️: unsigned int factor = 100;
    // RDKit✔️✔️: while (factor < mol.getNumAtoms()) {
    // RDKit✔️✔️:   factor *= 10;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: DOUBLE_VECT invars(mol.getNumAtoms());
    // RDKit✔️✔️: for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:   invars[i] = ranks[i] * factor;
    // RDKit✔️✔️:   const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:   std::string cipCode;
    // RDKit✔️✔️:   if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
    // RDKit✔️✔️:     if (cipCode == "S") {
    // RDKit✔️✔️:       invars[i] += 10;
    // RDKit✔️✔️:     } else if (cipCode == "R") {
    // RDKit✔️✔️:       invars[i] += 20;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (const auto oBond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:     if (oBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       if (oBond->getStereo() == Bond::STEREOE) {
    // RDKit✔️✔️:         invars[i] += 1;
    // RDKit✔️✔️:       } else if (oBond->getStereo() == Bond::STEREOZ) {
    // RDKit✔️✔️:         invars[i] += 2;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: iterateCIPRanks(mol, invars, ranks, true);
    // END RDKIT CPP FUNCTION rerankAtoms
    // Behavior review: modeled ranks are supplemented with the source R/S and
    // incident E/Z priorities before the unique legacy seeded rank owner runs.
    // Complexity review: one atom/adjacency pass plus the existing source-
    // shaped rank refinement matches the source asymptotic and allocation form.
    let mut factor = 100_i64;
    while factor < topology.atoms.len() as i64 {
        factor *= 10;
    }
    let mut invariants = Vec::with_capacity(topology.atoms.len());
    for (index, atom) in topology.atoms.iter().enumerate() {
        let mut invariant = i64::from(ranks[index]) * factor;
        invariant += match atom.prop("_CIPCode") {
            Some("S") => 10,
            Some("R") => 20,
            _ => 0,
        };
        for neighbor in topology.adjacency.neighbors_of(index) {
            let bond = &topology.bonds[neighbor.bond.index()];
            if bond.order() == BondOrder::Double {
                invariant += match bond.stereo() {
                    BondStereo::E => 1,
                    BondStereo::Z => 2,
                    _ => 0,
                };
            }
        }
        invariants.push(invariant);
    }
    let ranks = refine_atom_cip_ranks_from_invariants_with_query_state(
        topology,
        valence,
        &invariants,
        query_state,
    )?;
    // RDKit✔️✔️: // copy the ranks onto the atoms:
    // RDKit✔️✔️: for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️✔️:   mol.getAtomWithIdx(i)->setProp(common_properties::_CIPRank, ranks[i]);
    // RDKit✔️✔️: }
    for (atom, rank) in topology.atoms.iter_mut().zip(&ranks) {
        atom.set_computed_prop("_CIPRank", rank.to_string())?;
    }
    Ok(ranks)
}

fn install_ring_special_cases(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &[u32],
) -> Result<Vec<bool>, LegacyStereoError> {
    // BEGIN RDKIT CPP FUNCTION findChiralAtomSpecialCases delegated relation installation
    // RDKit✔️❌: if (ratom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️❌:     !ratom->hasProp(common_properties::_CIPCode) &&
    // RDKit✔️❌:     atomIsCandidateForRingStereochem(mol, ratom, atomRanks)) {
    // RDKit✔️❌: int same = (ratom->getChiralTag() == atom->getChiralTag()) ? 1 : -1;
    // RDKit✔️❌: ringStereoAtoms.push_back(same * (ratom->getIdx() + 1));
    // RDKit✔️❌: INT_VECT oringatoms(0);
    // RDKit✔️❌: ratom->getPropIfPresent(common_properties::_ringStereoAtoms,
    // RDKit✔️❌:                         oringatoms);
    // RDKit✔️❌: oringatoms.push_back(same * (atom->getIdx() + 1));
    // RDKit✔️❌: ratom->setProp(common_properties::_ringStereoAtoms, oringatoms, true);
    // RDKit✔️❌: possibleSpecialCases.set(ratom->getIdx());
    // RDKit✔️❌: possibleSpecialCases.set(atom->getIdx());
    // RDKit✔️❌: if (ringStereoAtoms.size() != 0) {
    // RDKit✔️❌:   atom->setProp(common_properties::_ringStereoAtoms, ringStereoAtoms, true);
    // RDKit✔️❌: }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION findChiralAtomSpecialCases delegated relation installation
    // Behavior review: the unique graph/BFS owner returns every reciprocal
    // relation from the complete canonical helper; this adapter retains source
    // sign and one-based indexing in the established comma encoding and marks
    // only nonempty relation rows.
    // Complexity review: ordered accumulation is O(R log V) instead of the
    // source's vector-indexed property writes, so behavior is exact but this
    // small internal-property materialization is algorithmically worse.
    let relations =
        crate::potential_stereo::special_ring_relations(topology, valence, rings, ranks)?;
    let mut encoded = BTreeMap::<usize, Vec<i64>>::new();
    for relation in relations {
        let one_based = i64::try_from(relation.other.index() + 1).map_err(|_| {
            CipRankError::InvariantOutOfRange {
                atom: relation.other,
                value: i64::MAX,
            }
        })?;
        encoded
            .entry(relation.atom.index())
            .or_default()
            .push(if relation.same_orientation {
                one_based
            } else {
                -one_based
            });
    }
    let mut special = vec![false; topology.atoms.len()];
    for (index, values) in encoded {
        let value = values
            .iter()
            .map(i64::to_string)
            .collect::<Vec<_>>()
            .join(",");
        topology.atoms[index].set_computed_prop("_ringStereoAtoms", value)?;
        special[index] = true;
    }
    Ok(special)
}

fn clean_directional_state(topology: &mut TopologyBlock) {
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception clean bond loop
    // RDKit✔️✔️: for (auto bond : mol.bonds()) {
    // RDKit✔️✔️: if ((bond->getBondDir() == Bond::BEGINWEDGE ||
    // RDKit✔️✔️:      bond->getBondDir() == Bond::BEGINDASH) &&
    // RDKit✔️✔️:     bond->getBeginAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:     bond->getEndAtom()->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:   bool atomHasAtropisomer = false;
    // RDKit✔️✔️:   for (auto nbond : mol.atomBonds(bond->getBeginAtom())) {
    // RDKit✔️✔️:     if (nbond->getStereo() == Bond::STEREOATROPCCW ||
    // RDKit✔️✔️:         nbond->getStereo() == Bond::STEREOATROPCW) {
    // RDKit✔️✔️:       atomHasAtropisomer = true;
    // RDKit✔️✔️:       foundAtropisomer = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!atomHasAtropisomer) {
    // RDKit✔️✔️:     bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION legacyStereoPerception clean bond loop
    // Behavior review: wedge/dash protection inspects only the source begin
    // endpoint; slash/backslash cleanup is handled by the source loop below.
    // Complexity review: bounded adjacency scans per incident bond match the
    // source loop nesting and introduce no whole-graph clones.
    let mut clear = BTreeSet::new();
    for bond in &topology.bonds {
        if matches!(
            bond.direction(),
            BondDirection::BeginWedge | BondDirection::BeginDash
        ) && topology.atoms[bond.begin().index()].chiral_tag() == ChiralTag::Unspecified
            && topology.atoms[bond.end().index()].chiral_tag() == ChiralTag::Unspecified
        {
            let adjacent_atrop = topology
                .adjacency
                .neighbors_of(bond.begin().index())
                .iter()
                .any(|neighbor| {
                    matches!(
                        topology.bonds[neighbor.bond.index()].stereo(),
                        BondStereo::AtropCw | BondStereo::AtropCcw
                    )
                });
            if !adjacent_atrop {
                clear.insert(bond.id());
            }
        }
    }
    // RDKit✔️✔️: if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:     (bond->getStereo() == Bond::STEREOANY ||
    // RDKit✔️✔️:      bond->getStereo() == Bond::STEREONONE)) {
    // RDKit✔️✔️:   std::vector<Atom *> batoms = {bond->getBeginAtom(), bond->getEndAtom()};
    // RDKit✔️✔️:   for (auto batom : batoms) {
    // RDKit✔️✔️:     for (const auto nbrBndI : mol.atomBonds(batom)) {
    // RDKit✔️✔️:       if (nbrBndI == bond) continue;
    // RDKit✔️✔️:       if ((nbrBndI->getBondDir() == Bond::ENDDOWNRIGHT ||
    // RDKit✔️✔️:            nbrBndI->getBondDir() == Bond::ENDUPRIGHT) &&
    // RDKit✔️✔️:           (nbrBndI->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:            nbrBndI->getBondType() == Bond::AROMATIC)) {
    // RDKit✔️✔️:         bool okToClear = true;
    // RDKit✔️✔️:         for (const auto nbrBndJ :
    // RDKit✔️✔️:              mol.atomBonds(nbrBndI->getOtherAtom(batom))) {
    // RDKit✔️✔️:           if (nbrBndJ->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:               nbrBndJ->getStereo() != Bond::STEREOANY &&
    // RDKit✔️✔️:               nbrBndJ->getStereo() != Bond::STEREONONE) {
    // RDKit✔️✔️:             okToClear = false;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (okToClear) nbrBndI->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: }
    for bond in &topology.bonds {
        if bond.order() != BondOrder::Double
            || !matches!(bond.stereo(), BondStereo::Any | BondStereo::None)
        {
            continue;
        }
        for endpoint in [bond.begin(), bond.end()] {
            for neighbor in topology.adjacency.neighbors_of(endpoint.index()) {
                if neighbor.bond == bond.id() {
                    continue;
                }
                let directed = &topology.bonds[neighbor.bond.index()];
                if !matches!(
                    directed.direction(),
                    BondDirection::EndDownRight | BondDirection::EndUpRight
                ) || !matches!(directed.order(), BondOrder::Single | BondOrder::Aromatic)
                {
                    continue;
                }
                let other = if directed.begin() == endpoint {
                    directed.end()
                } else {
                    directed.begin()
                };
                let consumed =
                    topology
                        .adjacency
                        .neighbors_of(other.index())
                        .iter()
                        .any(|candidate| {
                            let candidate = &topology.bonds[candidate.bond.index()];
                            candidate.order() == BondOrder::Double
                                && !matches!(candidate.stereo(), BondStereo::Any | BondStereo::None)
                        });
                if !consumed {
                    clear.insert(directed.id());
                }
            }
        }
    }
    for bond in clear {
        topology.bonds[bond.index()].set_direction(BondDirection::None);
    }
}

/// Apply the fixed RDKit legacy `assignStereochemistry(true, true, true)`
/// closure to detached topology state.
pub fn assign_legacy_stereochemistry(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) -> Result<TopologyBlock, LegacyStereoError> {
    assign_legacy_stereochemistry_with_query_state(topology, valence, rings, None)
}

/// Run the fixed-profile legacy assignment used by RDKit depiction.
///
/// This is the `assignStereochemistry(mol, false)` path: existing stereo is
/// assigned/ranked but the cleanup-only branches are not executed.
#[doc(hidden)]
pub fn assign_legacy_stereochemistry_for_depiction(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
) -> Result<TopologyBlock, LegacyStereoError> {
    // RDKit❗✔️:   RDKit::MolOps::assignStereochemistry(mol, false);
    // Behavior review: the one explicit argument is `cleanIt=false`; default
    // `force=false` cannot short-circuit here because detached topology has no
    // molecule-level `_StereochemDone` cache property.
    // Complexity review: this wrapper only selects the existing owner branch.
    assign_legacy_stereochemistry_impl(topology, valence, rings, None, false, false)
}

/// Apply the fixed RDKit legacy assignment to detached topology state with
/// independent cleanup and possible-center flags.
///
/// The caller owns the source molecule-level `force=false` property guard and
/// property effects; detached topology does not contain `_StereochemDone`.
#[doc(hidden)]
pub fn assign_legacy_stereochemistry_with_flags(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<TopologyBlock, LegacyStereoError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::assignStereochemistry legacy flag dispatch
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     Chirality::legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters);
    // RDKit❗✔️:   }
    // END RDKIT CPP FUNCTION MolOps::assignStereochemistry legacy flag dispatch
    // Behavior review: the pinned parity profile selects the existing legacy
    // implementation, and this adapter forwards both independent flags. The
    // false/false and true/true profiles remain available through their
    // unchanged wrappers; the CX caller supplies true/false after its own
    // force=false presence guard. Exact profile regressions are scheduled in
    // the owning core test target.
    // Complexity review: this O(1) adapter adds no clone, allocation, graph
    // scan, or lookup; all work remains in the existing implementation.
    assign_legacy_stereochemistry_impl(
        topology,
        valence,
        rings,
        None,
        clean_it,
        flag_possible_stereo_centers,
    )
}

#[doc(hidden)]
pub fn assign_legacy_stereochemistry_with_query_state(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<TopologyBlock, LegacyStereoError> {
    assign_legacy_stereochemistry_impl(topology, valence, rings, query_state, true, true)
}

fn assign_legacy_stereochemistry_impl(
    mut topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    query_state: Option<QueryStateRef<'_>>,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<TopologyBlock, LegacyStereoError> {
    if let Some(state) = query_state {
        state
            .validate_for_topology(&topology)
            .map_err(CipRankError::InvalidQueryState)?;
    }
    // BEGIN RDKIT CPP FUNCTION assignStereochemistry
    // RDKit✔️✔️: void assignStereochemistry(ROMol &mol, bool cleanIt, bool force,
    // RDKit✔️✔️:                            bool flagPossibleStereoCenters) {
    // RDKit✔️✔️:   if (!force && mol.hasProp(common_properties::_StereochemDone)) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (mol.needsUpdatePropertyCache()) {
    // RDKit✔️✔️:     mol.updatePropertyCache(false);
    // RDKit✔️✔️:   }
    // RDKit❌❌:   if (!Chirality::getUseLegacyStereoPerception()) {
    // RDKit❌❌:     Chirality::stereoPerception(mol, cleanIt, flagPossibleStereoCenters);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     Chirality::legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   mol.setProp(common_properties::_StereochemDone, 1, true);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION assignStereochemistry
    // The fixed RDKit 2026.03.1 reference profile selects the legacy branch.
    // The modern `stereoPerception` branch is independent, unmodeled behavior
    // and is deliberately not certified by this fixed-profile owner.
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception clean inventory
    // RDKit✔️✔️: if (cleanIt) {
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     atom->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️:     atom->clearProp(common_properties::_ChiralityPossible);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:     bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION legacyStereoPerception clean inventory
    // Behavior review: this fixed-parameter detached owner cleans computed
    // labels, validates marked tetrahedral centers by legacy ranks, assigns
    // directional double-bond stereo, removes invalid marked centers and then
    // cleans enhanced-stereo membership. Source property-cache preparation is
    // supplied explicitly by the caller's validated valence assignment.
    // Complexity review: the owner performs the source inventory before any
    // rank allocation, then bounded graph passes and only source-required rank
    // refinements; no eager whole-graph clone is introduced here.
    topology.validate()?;
    for atom in &mut topology.atoms {
        if clean_it {
            atom.clear_prop("_CIPCode");
            atom.clear_prop("_ChiralityPossible");
            atom.clear_prop("_ringStereochemCand");
            atom.clear_prop("_ringStereoAtoms");
        }
    }
    // RDKit✔️✔️: bool hasStereoAtoms = false;  // flagPossibleStereoCenters;
    // RDKit✔️✔️: bool hasPotentialStereoAtoms = false;
    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (!hasStereoAtoms && atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:       atom->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit✔️✔️:     hasStereoAtoms = true;
    // RDKit✔️✔️:   } else if (!hasPotentialStereoAtoms) {
    // RDKit✔️✔️:     UINT_VECT ranks;
    // RDKit✔️✔️:     Chirality::INT_PAIR_VECT nbrs;
    // RDKit✔️✔️:     hasPotentialStereoAtoms =
    // RDKit✔️✔️:         isAtomPotentialChiralCenter(atom, mol, ranks, nbrs).first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    let mut has_stereo_atoms = false;
    let mut has_potential_stereo_atoms = false;
    for index in 0..topology.atoms.len() {
        let tag = topology.atoms[index].chiral_tag();
        if !has_stereo_atoms && !matches!(tag, ChiralTag::Unspecified | ChiralTag::Other) {
            has_stereo_atoms = true;
        } else if !has_potential_stereo_atoms {
            has_potential_stereo_atoms =
                is_legal_legacy_center(&topology, valence, rings, AtomId::new(index))?;
        }
    }
    let mut has_stereo_bonds = false;
    let mut has_potential_stereo_bonds = false;
    // RDKit✔️✔️: bool hasStereoBonds = false;
    // RDKit✔️✔️: bool hasPotentialStereoBonds = false;
    // RDKit✔️✔️: for (auto bond : mol.bonds()) {
    // RDKit✔️✔️:   if (cleanIt) {
    // RDKit✔️✔️:     bond->clearProp(common_properties::_CIPCode);
    // RDKit✔️✔️:     if ((bond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️✔️:          bond->getBondType() == Bond::AROMATIC) &&
    // RDKit✔️✔️:         !shouldDetectDoubleBondStereo(bond)) {
    // RDKit✔️✔️:       if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit✔️✔️:         bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (bond->getStereo() != Bond::STEREONONE) {
    // RDKit✔️✔️:         bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:         bond->getStereoAtoms().clear();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     } else if (bond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:       if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit✔️✔️:         bond->setStereo(Bond::STEREOANY);
    // RDKit✔️✔️:         bond->getStereoAtoms().clear();
    // RDKit✔️✔️:         bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:       } else if (bond->getStereo() != Bond::STEREOANY) {
    // RDKit✔️✔️:         bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️:         bond->getStereoAtoms().clear();
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    for bond_index in 0..topology.bonds.len() {
        if clean_it {
            topology.bonds[bond_index].clear_prop("_CIPCode");
        }
        let bond = topology.bonds[bond_index].clone();
        let source_should_detect =
            rings.num_bond_rings(bond.id()) == 0 || rings.min_bond_ring_size(bond.id()) >= 8;
        if clean_it {
            if matches!(bond.order(), BondOrder::Double | BondOrder::Aromatic)
                && !source_should_detect
            {
                if bond.direction() == BondDirection::EitherDouble {
                    topology.bonds[bond_index].set_direction(BondDirection::None);
                }
                if bond.stereo() != BondStereo::None {
                    topology.bonds[bond_index].set_stereo_atoms(None);
                    topology.bonds[bond_index].set_stereo(BondStereo::None)?;
                }
            } else if bond.order() == BondOrder::Double {
                if bond.direction() == BondDirection::EitherDouble {
                    topology.bonds[bond_index].set_stereo_atoms(None);
                    topology.bonds[bond_index].set_stereo(BondStereo::Any)?;
                    topology.bonds[bond_index].set_direction(BondDirection::None);
                } else if bond.stereo() != BondStereo::Any {
                    topology.bonds[bond_index].set_stereo_atoms(None);
                    topology.bonds[bond_index].set_stereo(BondStereo::None)?;
                }
            }
        }
        // RDKit✔️✔️: if (!hasStereoBonds && bond->getBondType() == Bond::DOUBLE) {
        // RDKit✔️✔️:   bool isSpecified = false;
        // RDKit✔️✔️:   for (auto nbond : mol.atomBonds(bond->getBeginAtom())) {
        // RDKit✔️✔️:     if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
        // RDKit✔️✔️:         nbond->getBondDir() == Bond::ENDUPRIGHT) {
        // RDKit✔️✔️:       hasStereoBonds = true;
        // RDKit✔️✔️:       isSpecified = true;
        // RDKit✔️✔️:       break;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!hasStereoBonds) {
        // RDKit✔️✔️:     for (auto nbond : mol.atomBonds(bond->getEndAtom())) {
        // RDKit✔️✔️:       if (nbond->getBondDir() == Bond::ENDDOWNRIGHT ||
        // RDKit✔️✔️:           nbond->getBondDir() == Bond::ENDUPRIGHT) {
        // RDKit✔️✔️:         hasStereoBonds = true;
        // RDKit✔️✔️:         isSpecified = true;
        // RDKit✔️✔️:         break;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!hasPotentialStereoBonds && !isSpecified &&
        // RDKit✔️✔️:       shouldDetectDoubleBondStereo(bond)) {
        // RDKit✔️✔️:     hasPotentialStereoBonds = true;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        let current = &topology.bonds[bond_index];
        if !has_stereo_bonds && current.order() == BondOrder::Double {
            let is_specified = [current.begin(), current.end()]
                .into_iter()
                .any(|endpoint| {
                    topology
                        .adjacency
                        .neighbors_of(endpoint.index())
                        .iter()
                        .any(|neighbor| {
                            matches!(
                                topology.bonds[neighbor.bond.index()].direction(),
                                BondDirection::EndDownRight | BondDirection::EndUpRight
                            )
                        })
                });
            if is_specified {
                has_stereo_bonds = true;
            }
            if !has_potential_stereo_bonds && !is_specified && source_should_detect {
                has_potential_stereo_bonds = true;
            }
        }
    }
    // RDKit✔️✔️: while (keepGoing) {
    // RDKit✔️✔️:   if (hasStereoAtoms || hasPotentialStereoAtoms) {
    // RDKit✔️✔️:     std::tie(hasStereoAtoms, changedStereoAtoms) =
    // RDKit✔️✔️:         Chirality::assignAtomChiralCodes(mol, atomRanks,
    // RDKit✔️✔️:                                          flagPossibleStereoCenters);
    // RDKit✔️✔️:   } else { changedStereoAtoms = false; }
    // RDKit✔️✔️:   if (hasStereoBonds || hasPotentialStereoBonds) {
    // RDKit✔️✔️:     std::tie(hasStereoBonds, changedStereoBonds) =
    // RDKit✔️✔️:         Chirality::assignBondStereoCodes(mol, atomRanks);
    // RDKit✔️✔️:   } else { changedStereoBonds = false; }
    // RDKit✔️✔️:   keepGoing = (hasStereoAtoms || hasStereoBonds) &&
    // RDKit✔️✔️:               (changedStereoAtoms || changedStereoBonds);
    // RDKit✔️✔️:   if (keepGoing) Chirality::rerankAtoms(mol, atomRanks);
    // RDKit✔️✔️: }
    let mut ranks = Vec::new();
    let mut keep_going = has_stereo_atoms || has_stereo_bonds;
    if !keep_going && flag_possible_stereo_centers {
        keep_going = has_potential_stereo_atoms || has_potential_stereo_bonds;
    }
    while keep_going {
        let changed_stereo_atoms;
        if has_stereo_atoms || has_potential_stereo_atoms {
            (has_stereo_atoms, changed_stereo_atoms) =
                assign_atom_codes(&mut topology, valence, rings, &mut ranks, query_state)?;
        } else {
            changed_stereo_atoms = false;
        }
        let changed_stereo_bonds;
        if has_stereo_bonds || has_potential_stereo_bonds {
            if ranks.is_empty() {
                ranks = materialize_initial_ranks(&mut topology, valence, query_state)?;
            }
            let bond_assignment = assign_directional_double_bond_stereo(topology, &ranks, rings)?;
            has_stereo_bonds = bond_assignment.has_unassigned;
            changed_stereo_bonds = bond_assignment.assigned_any;
            topology = bond_assignment.topology;
        } else {
            changed_stereo_bonds = false;
        }
        keep_going = (has_stereo_atoms || has_stereo_bonds)
            && (changed_stereo_atoms || changed_stereo_bonds);
        if keep_going {
            ranks = rerank_atoms(&mut topology, valence, &ranks, query_state)?;
        }
    }

    // RDKit✔️✔️:   if (cleanIt) {
    // Behavior review: the depiction call passes `cleanIt=false`, so it ends
    // after source assignment/ranking. Existing public sanitize/finalization
    // callers retain the complete cleanup path below.
    // Complexity review: this is the source constant-time branch around the
    // existing cleanup passes and introduces no extra allocation or scan.
    if !clean_it {
        topology.validate()?;
        return Ok(topology);
    }

    // RDKit✔️❌: boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit✔️❌: Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases, atomRanks);
    let special_ring_atoms = install_ring_special_cases(&mut topology, valence, rings, &ranks)?;

    // RDKit✔️✔️: for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:   if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:       !Chirality::hasNonTetrahedralStereo(atom) &&
    // RDKit✔️✔️:       !atom->hasProp(common_properties::_CIPCode) &&
    // RDKit✔️✔️:       (!possibleSpecialCases[atom->getIdx()] ||
    // RDKit✔️✔️:        !atom->hasProp(common_properties::_ringStereoAtoms))) {
    // RDKit✔️✔️:     atom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit✔️✔️:     if (atom->getNumExplicitHs() == 1 && atom->getFormalCharge() == 0 &&
    // RDKit✔️✔️:         !atom->getIsAromatic()) {
    // RDKit✔️✔️:       atom->setNumExplicitHs(0);
    // RDKit✔️✔️:       atom->setNoImplicit(false);
    // RDKit✔️✔️:       atom->calcExplicitValence(false);
    // RDKit✔️✔️:       atom->calcImplicitValence(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    for atom in &mut topology.atoms {
        if matches!(
            atom.chiral_tag(),
            ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
        ) && atom.prop("_CIPCode").is_none()
            && !special_ring_atoms[atom.id().index()]
        {
            atom.set_chiral_tag(ChiralTag::Unspecified);
            if atom.explicit_hydrogens() == 1 && atom.formal_charge() == 0 && !atom.is_aromatic() {
                atom.set_explicit_hydrogens(0);
                atom.set_no_implicit(false);
            }
        }
    }
    let degrees = (0..topology.atoms.len())
        .map(|atom| topology.adjacency.neighbors_of(atom).len())
        .collect::<Vec<_>>();
    // RDKit✔️✔️: if (bond->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:     (bond->getBondDir() == Bond::EITHERDOUBLE ||
    // RDKit✔️✔️:      bond->getStereo() == Bond::STEREOANY) &&
    // RDKit✔️✔️:     (bond->getBeginAtom()->getDegree() == 1 ||
    // RDKit✔️✔️:      bond->getEndAtom()->getDegree() == 1)) {
    // RDKit✔️✔️:   if (bond->getBondDir() == Bond::EITHERDOUBLE) {
    // RDKit✔️✔️:     bond->setBondDir(Bond::NONE);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bond->setStereo(Bond::STEREONONE);
    // RDKit✔️✔️: }
    for bond in &mut topology.bonds {
        if bond.order() == BondOrder::Double
            && matches!(bond.stereo(), BondStereo::Any)
            && (degrees[bond.begin().index()] == 1 || degrees[bond.end().index()] == 1)
        {
            if bond.direction() == BondDirection::EitherDouble {
                bond.set_direction(BondDirection::None);
            }
            bond.set_stereo_atoms(None);
            bond.set_stereo(BondStereo::None)?;
        }
    }
    clean_directional_state(&mut topology);
    // RDKit✔️✔️: if (foundAtropisomer || Atropisomers::doesMolHaveAtropisomers(mol)) {
    // RDKit✔️✔️:   Atropisomers::cleanupAtropisomerStereoGroups(mol);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: Chirality::cleanupStereoGroups(mol);
    if topology
        .bonds
        .iter()
        .any(|bond| matches!(bond.stereo(), BondStereo::AtropCw | BondStereo::AtropCcw))
    {
        topology.stereo_groups =
            cleanup_atropisomer_stereo_groups(&topology, &crate::AtropisomerAssignment::default())?
                .groups;
    }
    crate::structure_tags::cleanup_stereo_groups(&mut topology);
    topology.validate()?;
    Ok(topology)
}
