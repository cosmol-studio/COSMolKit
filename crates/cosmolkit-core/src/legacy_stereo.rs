//! RDKit legacy stereochemistry assignment over detached topology values.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_model::{
    AtomId, AtomPropertyError, BondValueError, QueryStateRef, TopologyBlock,
    TopologyValidationError,
};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

use crate::{
    AtropisomerError, CipRankError, DoubleBondStereoError, PotentialStereoError, RingFindingError,
    RingInfo, StereoOrderError, ValenceAssignment, ValenceError,
    assign_atom_cip_ranks_with_query_state, assign_directional_double_bond_stereo,
    bond_affects_atom_chirality, cleanup_atropisomer_stereo_groups, count_swaps_to_interconvert,
    invert_tetrahedral_tag, is_atom_bridgehead_from_topology,
    refine_atom_cip_ranks_from_invariants_with_query_state,
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
    RingFinding(#[from] RingFindingError),
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
            // Source getValence observes the stored signed-width E field.
            // Reuse its existing owner; no fresh assignment or allocation.
            let explicit = crate::valence::cached_explicit_valence(atom, Some(valence))?;
            explicit == 4 || (explicit == 3 && atom.formal_charge() == 1)
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
    // BEGIN RDKIT CPP FUNCTION materialize_initial_ranks pinned Chirality.cpp:1325-1345
    // RDKit✔️✔️: void assignAtomCIPRanks(const ROMol &mol, UINT_VECT &ranks) {
    // RDKit✔️✔️:   PRECONDITION((!ranks.size() || ranks.size() >= mol.getNumAtoms()),
    // RDKit✔️✔️:                "bad ranks size");
    // RDKit✔️✔️:   if (!ranks.size()) {
    // RDKit✔️✔️:     ranks.resize(mol.getNumAtoms());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   unsigned int numAtoms = mol.getNumAtoms();
    // RDKit✔️✔️: #ifndef USE_NEW_STEREOCHEMISTRY
    // RDKit✔️✔️:   // get the initial invariants:
    // RDKit✔️✔️:   DOUBLE_VECT invars(numAtoms, 0);
    // RDKit✔️✔️:   buildCIPInvariants(mol, invars);
    // RDKit✔️✔️:   iterateCIPRanks(mol, invars, ranks, false);
    // RDKit❌❌: #else
    // RDKit❌❌:   Canon::chiralRankMolAtoms(mol, ranks);
    // RDKit❌❌: #endif
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // copy the ranks onto the atoms:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < numAtoms; ++i) {
    // RDKit❗✔️:     mol[i]->setProp(common_properties::_CIPRank, ranks[i], 1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION materialize_initial_ranks
    // Behavior: the pinned legacy rank engine returns segmented u32 ranks;
    // proposed exact UInt preserves unsigned width and computed membership.
    // Modern compile-time ranking is unmodeled. Complexity: the same single
    // refinement plus O(V) writes; numeric conversion removes temporary text.
    let ranks = assign_atom_cip_ranks_with_query_state(topology, valence, query_state)?;
    for (atom, rank) in topology.atoms.iter_mut().zip(&ranks) {
        let value = cosmolkit_model::PropertyValue::UInt(*rank);
        atom.set_computed_prop("_CIPRank", value)?;
    }
    Ok(ranks)
}

fn assign_atom_codes(
    topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    ranks: &mut Vec<u32>,
    flag_possible_stereo_centers: bool,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<(bool, bool), LegacyStereoError> {
    // BEGIN RDKIT CPP FUNCTION assign_atom_codes pinned Chirality.cpp:1741-1821
    // RDKit✔️✔️: std::pair<bool, bool> assignAtomChiralCodes(ROMol &mol, UINT_VECT &ranks,
    // RDKit✔️✔️:                                             bool flagPossibleStereoCenters) {
    // RDKit✔️✔️:   PRECONDITION((!ranks.size() || ranks.size() == mol.getNumAtoms()),
    // RDKit✔️✔️:                "bad rank vector size");
    // RDKit✔️✔️:   bool atomChanged = false;
    // RDKit✔️✔️:   unsigned int unassignedAtoms = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // ------------------
    // RDKit✔️✔️:   // now loop over each atom and, if it's marked as chiral,
    // RDKit✔️✔️:   //  figure out the appropriate CIP label:
    // RDKit✔️✔️:   for (auto atom : mol.atoms()) {
    // RDKit✔️✔️:     Atom::ChiralType tag = atom->getChiralTag();
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // only worry about this atom if it has a marked chirality
    // RDKit✔️✔️:     // we understand:
    // RDKit✔️✔️:     if (flagPossibleStereoCenters ||
    // RDKit✔️✔️:         (tag != Atom::CHI_UNSPECIFIED && tag != Atom::CHI_OTHER)) {
    // RDKit✔️✔️:       if (atom->hasProp(common_properties::_CIPCode)) {
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (!ranks.size()) {
    // RDKit✔️✔️:         //  if we need to, get the "CIP" ranking of each atom:
    // RDKit✔️✔️:         assignAtomCIPRanks(mol, ranks);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       Chirality::INT_PAIR_VECT nbrs;
    // RDKit✔️✔️:       // note that hasDupes is only evaluated if legalCenter==true
    // RDKit✔️✔️:       auto [legalCenter, hasDupes] =
    // RDKit✔️✔️:           isAtomPotentialChiralCenter(atom, mol, ranks, nbrs);
    // RDKit✔️✔️:       if (legalCenter) {
    // RDKit✔️✔️:         ++unassignedAtoms;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (legalCenter && !hasDupes && flagPossibleStereoCenters) {
    // RDKit✔️✔️:         atom->setProp(common_properties::_ChiralityPossible, 1);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       if (legalCenter && !hasDupes && tag != Atom::CHI_UNSPECIFIED &&
    // RDKit✔️✔️:           tag != Atom::CHI_OTHER) {
    // RDKit✔️✔️:         // stereochem is possible and we have no duplicate neighbors, assign
    // RDKit✔️✔️:         // a CIP code:
    // RDKit✔️✔️:         atomChanged = true;
    // RDKit✔️✔️:         --unassignedAtoms;
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // sort the list of neighbors by their CIP ranks:
    // RDKit✔️✔️:         std::sort(nbrs.begin(), nbrs.end(), Rankers::pairLess);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // collect the list of neighbor indices:
    // RDKit✔️✔️:         std::list<int> nbrIndices;
    // RDKit✔️✔️:         for (Chirality::INT_PAIR_VECT_CI nbrIt = nbrs.begin();
    // RDKit✔️✔️:              nbrIt != nbrs.end(); ++nbrIt) {
    // RDKit✔️✔️:           nbrIndices.push_back((*nbrIt).second);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // ask the atom how many swaps we have to make:
    // RDKit✔️✔️:         int nSwaps = atom->getPerturbationOrder(nbrIndices);
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if the atom has 3 neighbors and a hydrogen, add a swap:
    // RDKit✔️✔️:         if (nbrIndices.size() == 3 && atom->getTotalNumHs() == 1) {
    // RDKit✔️✔️:           ++nSwaps;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:
    // RDKit✔️✔️:         // if that number is odd, we'll change our chirality:
    // RDKit✔️✔️:         if (nSwaps % 2) {
    // RDKit✔️✔️:           if (tag == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:             tag = Atom::CHI_TETRAHEDRAL_CW;
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             tag = Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         // now assign the CIP code:
    // RDKit✔️✔️:         std::string cipCode;
    // RDKit✔️✔️:         if (tag == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:           cipCode = "S";
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           cipCode = "R";
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         atom->setProp(common_properties::_CIPCode, cipCode);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return std::make_pair((unassignedAtoms > 0), atomChanged);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION assign_atom_codes
    // Behavior: the independent flag/tag guard runs before CIP presence,
    // lazy ranking and legality. Possible is ordinary Int1 when requested;
    // CIP is ordinary String. Existing setters retain computed membership.
    // Complexity: one lazy rank vector, bounded-degree neighbor collection/
    // sorting and linear atom scan. The existing tree property storage and
    // BTreeSet duplicate check differ from source constant/vector lookups;
    // no additional whole-graph clone or production observer is introduced.
    let mut changed = false;
    let mut unassigned = 0usize;
    for index in 0..topology.atoms.len() {
        let tag = topology.atoms[index].chiral_tag();
        if !flag_possible_stereo_centers && matches!(tag, ChiralTag::Unspecified | ChiralTag::Other)
        {
            continue;
        }
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
        if legal && !duplicates && flag_possible_stereo_centers {
            topology.atoms[index].set_prop("_ChiralityPossible", 1_i32)?;
        }
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
            topology.atoms[index].set_prop(
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
    // BEGIN RDKIT CPP FUNCTION rerank_atoms pinned Chirality.cpp:2067-2117
    // RDKit✔️✔️: void rerankAtoms(const ROMol &mol, UINT_VECT &ranks) {
    // RDKit✔️✔️:   PRECONDITION(ranks.size() == mol.getNumAtoms(), "bad rank vector size");
    // RDKit✔️✔️:   unsigned int factor = 100;
    // RDKit✔️✔️:   while (factor < mol.getNumAtoms()) {
    // RDKit✔️✔️:     factor *= 10;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit❌❌: #ifdef VERBOSE_CANON
    // RDKit❌❌:   BOOST_LOG(rdDebugLog) << "rerank PRE: " << std::endl;
    // RDKit❌❌:   for (int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❌❌:     BOOST_LOG(rdDebugLog) << "  " << i << ": " << ranks[i] << std::endl;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️:
    // RDKit✔️✔️:   DOUBLE_VECT invars(mol.getNumAtoms());
    // RDKit✔️✔️:   // and now supplement them:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     invars[i] = ranks[i] * factor;
    // RDKit✔️✔️:     const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:     // Priority order: R > S > nothing
    // RDKit✔️✔️:     std::string cipCode;
    // RDKit✔️✔️:     if (atom->getPropIfPresent(common_properties::_CIPCode, cipCode)) {
    // RDKit✔️✔️:       if (cipCode == "S") {
    // RDKit✔️✔️:         invars[i] += 10;
    // RDKit✔️✔️:       } else if (cipCode == "R") {
    // RDKit✔️✔️:         invars[i] += 20;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (const auto oBond : mol.atomBonds(atom)) {
    // RDKit✔️✔️:       if (oBond->getBondType() == Bond::DOUBLE) {
    // RDKit✔️✔️:         if (oBond->getStereo() == Bond::STEREOE) {
    // RDKit✔️✔️:           invars[i] += 1;
    // RDKit✔️✔️:         } else if (oBond->getStereo() == Bond::STEREOZ) {
    // RDKit✔️✔️:           invars[i] += 2;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   iterateCIPRanks(mol, invars, ranks, true);
    // RDKit✔️✔️:   // copy the ranks onto the atoms:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❗✔️:     mol.getAtomWithIdx(i)->setProp(common_properties::_CIPRank, ranks[i]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit❌❌: #ifdef VERBOSE_CANON
    // RDKit❌❌:   BOOST_LOG(rdDebugLog) << "   post: " << std::endl;
    // RDKit❌❌:   for (int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit❌❌:     BOOST_LOG(rdDebugLog) << "  " << i << ": " << ranks[i] << std::endl;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION rerank_atoms
    // Behavior: the existing seeded rank engine consumes exact source CIP
    // and E/Z supplements; ordinary numeric overwrite preserves pre-existing
    // computed membership. VERBOSE_CANON logging is not modeled. Complexity:
    // O(V+E) invariant preparation plus the same rank refinement and O(V)
    // writes; no additional ranking or temporary numeric String allocation.
    let mut factor = 100_i64;
    while factor < topology.atoms.len() as i64 {
        factor *= 10;
    }
    let mut invariants = Vec::with_capacity(topology.atoms.len());
    for (index, atom) in topology.atoms.iter().enumerate() {
        let mut invariant = i64::from(ranks[index]) * factor;
        invariant += match atom
            .prop("_CIPCode")
            .and_then(|value| value.as_string().ok())
            .map(|value| value.as_bytes())
        {
            Some(b"S") => 10,
            Some(b"R") => 20,
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
    for (atom, rank) in topology.atoms.iter_mut().zip(&ranks) {
        let value = cosmolkit_model::PropertyValue::UInt(*rank);
        atom.set_prop("_CIPRank", value)?;
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
    // Apply the one engine's detached changes in source property insertion
    // order; ordinary writes retain computed membership in the model store.
    let result = crate::potential_stereo::special_ring_cases(topology, valence, rings, ranks)?;
    for update in result.updates {
        let atom = &mut topology.atoms[update.atom.index()];
        if update.computed {
            atom.set_computed_prop(update.key, update.value)?;
        } else {
            atom.set_prop(update.key, update.value)?;
        }
    }
    Ok(result.flags)
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
    // `force=false` presence guard belongs to the caller that borrows molecule
    // properties. This topology-only delegator runs only when that caller
    // dispatches; it cannot inspect molecule-level `_StereochemDone` itself.
    // Complexity review: this wrapper only selects the existing owner branch.
    assign_legacy_stereochemistry_impl(topology, valence, rings, None, false, false)
        .map(|assignment| assignment.topology)
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
    .map(|assignment| assignment.topology)
}

#[doc(hidden)]
pub fn assign_legacy_stereochemistry_with_query_state(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    query_state: Option<QueryStateRef<'_>>,
) -> Result<TopologyBlock, LegacyStereoError> {
    assign_legacy_stereochemistry_impl(topology, valence, rings, query_state, true, true)
        .map(|assignment| assignment.topology)
}

/// Detached state effects of the fixed legacy stereochemistry profile.
#[derive(Debug, Clone, PartialEq)]
pub struct LegacyStereoAssignment {
    pub topology: TopologyBlock,
    /// Source-required ring preparation, absent when the input cache is reused.
    pub ring_update: Option<RingInfo>,
}

/// Run legacy stereochemistry and retain its detached ring-state effects.
#[doc(hidden)]
pub fn assign_legacy_stereochemistry_with_assignments(
    topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<LegacyStereoAssignment, LegacyStereoError> {
    // RDKit❗✔️:     Chirality::legacyStereoPerception(mol, cleanIt, flagPossibleStereoCenters);
    // One owner returns both the topology and its source ring preparation;
    // this adapter neither copies a cache nor executes a second algorithm.
    assign_legacy_stereochemistry_impl(
        topology,
        valence,
        rings,
        None,
        clean_it,
        flag_possible_stereo_centers,
    )
}

fn assign_legacy_stereochemistry_impl(
    mut topology: TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    query_state: Option<QueryStateRef<'_>>,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
) -> Result<LegacyStereoAssignment, LegacyStereoError> {
    let mut ring_update = None;
    assign_legacy_stereochemistry_source(
        &mut topology,
        valence,
        rings,
        query_state,
        clean_it,
        flag_possible_stereo_centers,
        &mut ring_update,
    )?;
    Ok(LegacyStereoAssignment {
        topology,
        ring_update,
    })
}

/// Borrow the actual detached graph for reached native source operations.
/// Property/stereo mutations preceding an error remain observable to the caller;
/// no empty replacement graph or copied working graph stands in for that state.
/// The output ring preparation is retained even on a later property/stereo
/// failure; a source caller moves that effect to its actual cache before
/// propagating the error, while old owning APIs expose effects on success.
#[doc(hidden)]
pub fn assign_legacy_stereochemistry_source(
    mut topology: &mut TopologyBlock,
    valence: &ValenceAssignment,
    rings: &RingInfo,
    query_state: Option<QueryStateRef<'_>>,
    clean_it: bool,
    flag_possible_stereo_centers: bool,
    ring_update: &mut Option<RingInfo>,
) -> Result<(), LegacyStereoError> {
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
    // BEGIN RDKIT CPP FUNCTION legacyStereoPerception ring preparation
    // RDKit✔️✔️:   // later we're going to need ring information, get it now if we don't
    // RDKit✔️✔️:   // have it already:
    // RDKit✔️✔️:   // NOTE, if called from the SMART code, the ring info will be DUMMY, and
    // RDKit✔️✔️:   // contains no information
    // RDKit✔️✔️:   if (!mol.getRingInfo()->isFindFastOrBetter()) {
    // RDKit✔️✔️:     MolOps::fastFindRings(mol);
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION legacyStereoPerception ring preparation
    // A borrowed Fast-or-better cache is reused without allocation. The
    // OtherOrUnknown/absent path invokes the existing full-topology core
    // owner and keeps its result for every subsequent ring consumer; this
    // adds only the source-required O(V+E) search on that branch.
    let source_rings = rings;
    if !rings.is_find_fast_or_better() {
        *ring_update = Some(crate::fast_find_rings(&topology)?);
    }
    let rings = ring_update.as_ref().unwrap_or(source_rings);
    for atom in &mut topology.atoms {
        if clean_it {
            atom.clear_prop("_CIPCode")?;
            atom.clear_prop("_ChiralityPossible")?;
            atom.clear_prop("_ringStereochemCand")?;
            atom.clear_prop("_ringStereoAtoms")?;
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
    let mut has_unassigned_double_bond = false;
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
            topology.bonds[bond_index].clear_prop("_CIPCode")?;
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
        has_unassigned_double_bond |=
            current.order() == BondOrder::Double && current.stereo() == BondStereo::None;
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
            #[cfg(test)]
            let observed_before = state_owner_tests::before(&topology);
            (has_stereo_atoms, changed_stereo_atoms) = assign_atom_codes(
                &mut topology,
                valence,
                rings,
                &mut ranks,
                flag_possible_stereo_centers,
                query_state,
            )?;
            #[cfg(test)]
            state_owner_tests::after(
                observed_before,
                &topology,
                &ranks,
                (has_stereo_atoms, changed_stereo_atoms),
            );
        } else {
            changed_stereo_atoms = false;
        }
        let changed_stereo_bonds;
        // BEGIN RDKIT CPP FUNCTION assignBondStereoCodes lazy-rank boundary
        // RDKit✔️✔️:   for (auto dblBond : mol.bonds()) {
        // RDKit✔️✔️:     if (dblBond->getBondType() == Bond::BondType::DOUBLE) {
        // RDKit✔️✔️:       if (dblBond->getStereo() != Bond::BondStereo::STEREONONE) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (!ranks.size()) {
        // RDKit✔️✔️:         assignAtomCIPRanks(mol, ranks);
        // RDKit✔️✔️:       }
        // END RDKIT CPP FUNCTION assignBondStereoCodes lazy-rank boundary
        // The existing inventory records whether any double bond reaches
        // this source branch. With no ranks and all double bonds already
        // assigned, the source bond pass changes no state and returns false,
        // false. Preserve the absent rank properties used by depiction.
        // Reuse the inventory boolean; no additional graph traversal or rank
        // allocation occurs on this branch. The existing owner handles every
        // branch requiring ranks, with unchanged validation and effects.
        if has_stereo_bonds || has_potential_stereo_bonds {
            if ranks.is_empty() && !has_unassigned_double_bond {
                has_stereo_bonds = false;
                changed_stereo_bonds = false;
            } else {
                if ranks.is_empty() {
                    ranks = materialize_initial_ranks(&mut topology, valence, query_state)?;
                }
                (has_stereo_bonds, changed_stereo_bonds) =
                    crate::double_stereo::assign_directional_double_bond_stereo_source(
                        topology, &ranks, rings,
                    )?;
            }
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
        return Ok(());
    }

    // RDKit✔️❌: boost::dynamic_bitset<> possibleSpecialCases(mol.getNumAtoms());
    // RDKit✔️❌: Chirality::findChiralAtomSpecialCases(mol, possibleSpecialCases, atomRanks);
    // BEGIN RDKIT CPP FUNCTION findChiralAtomSpecialCases ring preparation
    // RDKit✔️❌:   if (!mol.getRingInfo()->isSymmSssr()) {
    // RDKit✔️❌:     VECT_INT_VECT sssrs;
    // RDKit✔️❌:     MolOps::symmetrizeSSSR(mol, sssrs);
    // RDKit✔️❌:   }
    // END RDKIT CPP FUNCTION findChiralAtomSpecialCases ring preparation
    // This source transition precedes every chiral-atom guard. Retain its
    // ring update for the caller as well as all cleanup consumers below.
    // The existing detached symmetrization owner reconstructs its SSSR
    // context instead of reusing the source molecule's cached extra rings.
    if !rings.is_symm_sssr() {
        *ring_update = Some(crate::symmetrized_sssr(
            &topology,
            &crate::RingSearchParams::default(),
        )?);
    }
    let rings = ring_update.as_ref().unwrap_or(source_rings);
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
    Ok(())
}

#[cfg(test)]
mod legacy_ring_prepass_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::Element;

    fn two_six_membered_rings() -> TopologyBlock {
        let atoms = (0..12)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = (0..12)
            .map(|index| {
                let base = if index < 6 { 0 } else { 6 };
                let offset = index - base;
                let next = base + (offset + 1) % 6;
                let order = if index == 0 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                };
                let mut spec = BondSpec::new(AtomId::new(index), AtomId::new(next), order);
                if index == 0 {
                    spec = spec.with_stereo(BondStereo::Any);
                }
                Bond::from_spec(BondId::new(index), spec)
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("two source-aligned six-membered cycles")
    }

    #[test]
    fn legacy_depiction_already_assigned_bonds_do_not_materialize_ranks() {
        let atoms = (0..4)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        let mut bonds = (0..3)
            .map(|i| {
                let order = if i == 1 {
                    BondOrder::Double
                } else {
                    BondOrder::Single
                };
                let mut bond = Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(i), AtomId::new(i + 1), order),
                );
                if i != 1 {
                    bond.set_direction(BondDirection::EndUpRight);
                }
                bond
            })
            .collect::<Vec<_>>();
        bonds[1].set_stereo_atoms(Some([AtomId::new(0), AtomId::new(3)]));
        bonds[1].set_stereo(BondStereo::E).unwrap();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
        let rings =
            crate::symmetrized_sssr(&topology, &crate::RingSearchParams::default()).unwrap();
        let valence = crate::assign_valence_with_options_for_topology(
            &topology,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let result =
            assign_legacy_stereochemistry_for_depiction(topology.clone(), &valence, &rings)
                .unwrap();
        // Pinned assignBondStereoCodes continues past already assigned double
        // bonds before its lazy assignAtomCIPRanks call. No new rank property
        // may influence subsequent depiction ordering on this source branch.
        assert!(result.atoms.iter().all(|a| a.prop("_CIPRank").is_none()));
        assert_eq!(result, topology);
    }

    fn assign_with_rings(topology: &TopologyBlock, rings: &RingInfo) -> TopologyBlock {
        let valence = crate::assign_valence_with_options_for_topology(
            topology,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .expect("non-strict valence for fixed ring topology");
        assign_legacy_stereochemistry_with_flags(topology.clone(), &valence, rings, true, false)
            .expect("pinned legacy clean=true possible=false profile")
    }

    #[test]
    fn legacy_ring_prepass_recomputes_initialized_empty_other_or_unknown() {
        let topology = two_six_membered_rings();
        let before = topology.clone();
        let rings = crate::ring_info_from_selected_rows(12, 12, &[], &[]).unwrap();
        let rings_before = rings.clone();
        assert!(rings.is_initialized());
        assert!(!rings.is_find_fast_or_better());

        let result = assign_with_rings(&topology, &rings);
        // Pinned RDKit 2026.03.1 legacy AssignStereochemistry on
        // C1=CCCCC1.C1CCCCC1 changes bond 0 STEREOANY to STEREONONE.
        assert_eq!(result.bonds[0].stereo(), BondStereo::None);
        assert_eq!(topology, before);
        assert_eq!(rings, rings_before);
    }

    #[test]
    fn legacy_ring_prepass_recomputes_nonempty_retained_other_or_unknown() {
        let topology = two_six_membered_rings();
        let source_rings = crate::fast_find_rings(&topology).unwrap();
        assert_eq!(source_rings.num_rings(), 2);
        let retained_index = source_rings
            .atom_rings()
            .iter()
            .position(|row| row.iter().all(|atom| atom.index() >= 6))
            .expect("second source ring exists");
        let retained = crate::ring_info_from_selected_rows(
            12,
            12,
            &[source_rings.atom_rings()[retained_index].clone()],
            &[source_rings.bond_rings()[retained_index].clone()],
        )
        .unwrap();
        let retained_before = retained.clone();
        assert_eq!(retained.num_bond_rings(BondId::new(0)), 0);

        let result = assign_with_rings(&topology, &retained);
        assert_eq!(result.bonds[0].stereo(), BondStereo::None);
        assert_eq!(
            retained, retained_before,
            "source retained rows are immutable"
        );
    }

    #[test]
    fn legacy_ring_prepass_reuses_fast_or_better_source_membership() {
        let topology = two_six_membered_rings();
        let before = topology.clone();
        let rings = crate::fast_find_rings(&topology).unwrap();
        let rings_before = rings.clone();
        assert!(rings.is_find_fast_or_better());
        assert_eq!(rings.num_rings(), 2);

        let result = assign_with_rings(&topology, &rings);
        assert_eq!(result.bonds[0].stereo(), BondStereo::None);
        assert_eq!(topology, before);
        assert_eq!(rings, rings_before);
    }
}

#[cfg(test)]
mod state_owner_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondId, BondSpec, PropertyValue};
    use cosmolkit_types::Element;
    use std::cell::RefCell;

    #[derive(Debug)]
    struct Observation {
        before: TopologyBlock,
        after: TopologyBlock,
        ranks: Vec<u32>,
        flags: (bool, bool),
    }
    thread_local! {
        static OBSERVATIONS: RefCell<Option<Vec<Observation>>> = const { RefCell::new(None) };
    }
    pub(super) fn before(topology: &TopologyBlock) -> Option<TopologyBlock> {
        OBSERVATIONS.with(|slot| slot.borrow().as_ref().map(|_| topology.clone()))
    }
    pub(super) fn after(
        before: Option<TopologyBlock>,
        topology: &TopologyBlock,
        ranks: &[u32],
        flags: (bool, bool),
    ) {
        if let Some(before) = before {
            OBSERVATIONS.with(|slot| {
                slot.borrow_mut().as_mut().unwrap().push(Observation {
                    before,
                    after: topology.clone(),
                    ranks: ranks.to_vec(),
                    flags,
                })
            });
        }
    }

    struct Fixture {
        source_pin: String,
        states: Vec<Vec<Vec<String>>>,
        cells: BTreeMap<String, Cell>,
    }
    struct Cell {
        returns: Vec<Returned>,
        native_inner_calls: usize,
        states: BTreeMap<String, usize>,
    }
    struct Returned {
        ordinal: usize,
        possible: bool,
        has_unassigned: bool,
        changed: bool,
        ranks: Vec<u32>,
    }
    fn parse<T: std::str::FromStr>(text: &str) -> T
    where
        T::Err: std::fmt::Debug,
    {
        text.parse().unwrap()
    }
    fn unquote(text: &str) -> String {
        assert!(
            !text.contains('\\'),
            "these fixed source keys/values contain no escapes"
        );
        text.strip_prefix('"')
            .unwrap()
            .strip_suffix('"')
            .unwrap()
            .to_owned()
    }
    fn topology(rows: &[Vec<String>]) -> TopologyBlock {
        let mut atoms = Vec::new();
        let mut bonds = Vec::new();
        for row in rows {
            if row[0] == "ATOM" {
                let index = parse(&row[1]);
                assert_eq!(index, atoms.len());
                let mut spec = AtomSpec::new(Element::from_atomic_number(parse(&row[2])).unwrap())
                    .with_isotope(parse(&row[3]))
                    .with_formal_charge(parse(&row[4]))
                    .with_explicit_hydrogens(parse(&row[5]))
                    .with_no_implicit(row[6] == "1")
                    .with_radical_electrons(parse(&row[7]))
                    .with_aromatic(row[8] == "1")
                    .with_hybridization(Hybridization::from_rdkit_code(parse(&row[9])).unwrap())
                    .with_chiral_tag(ChiralTag::from_rdkit_code(parse(&row[10])).unwrap());
                if row[11] != "0" {
                    spec = spec.with_atom_map(parse(&row[11]));
                }
                atoms.push(Atom::from_spec(AtomId::new(index), spec));
            } else if row[0] == "BOND" {
                let index = parse(&row[1]);
                assert_eq!(index, bonds.len());
                let mut spec = BondSpec::new(
                    AtomId::new(parse(&row[2])),
                    AtomId::new(parse(&row[3])),
                    BondOrder::from_rdkit_code(parse(&row[4])).unwrap(),
                )
                .with_aromatic(row[5] == "1")
                .with_conjugated(row[6] == "1")
                .with_direction(BondDirection::from_rdkit_code(parse(&row[7])).unwrap())
                .with_stereo(BondStereo::from_rdkit_code(parse(&row[8])).unwrap());
                let stereo = row[9]
                    .split(',')
                    .filter(|s| !s.is_empty())
                    .map(|s| AtomId::new(parse(s)))
                    .collect::<Vec<_>>();
                if !stereo.is_empty() {
                    assert_eq!(stereo.len(), 2);
                    spec = spec.with_stereo_atoms(stereo[0], stereo[1]);
                }
                bonds.push(Bond::from_spec(BondId::new(index), spec));
            }
        }
        let mut result = TopologyBlock::try_from_parts(atoms, bonds, vec![], vec![]).unwrap();
        for row in rows.iter().filter(|r| r[0] == "PROP" && r[1] != "MOL") {
            let key = unquote(&row[2]);
            if key == "__computedProps" {
                continue;
            }
            let value = match row[3].as_str() {
                "1" => PropertyValue::Int(parse(&row[5])),
                "6" => PropertyValue::UInt(parse::<u32>(&row[5])),
                "3" => PropertyValue::String((unquote(&row[5])).into()),
                kind => panic!("unmodeled fixed atom/bond kind {kind}"),
            };
            let index: usize = parse(&row[1][1..]);
            if row[1].starts_with('A') {
                if row[4] == "1" {
                    result.atoms[index].set_computed_prop(key, value).unwrap();
                } else {
                    result.atoms[index].set_prop(key, value).unwrap();
                }
            } else {
                assert!(row[1].starts_with('B'));
                if row[4] == "1" {
                    result.bonds[index].set_computed_prop(key, value).unwrap();
                } else {
                    result.bonds[index].set_prop(key, value).unwrap();
                }
            }
        }
        // MOL UInt/Int/computed fields are preserved in the source fixture but
        // deliberately not admitted as detached core values: their destination
        // still requires the separate human model decision.
        for row in rows
            .iter()
            .filter(|r| r[0] == "PROP" && r[1].starts_with('A') && r[2] == "\"__computedProps\"")
        {
            let atom = &result.atoms[parse::<usize>(&row[1][1..])];
            let source_keys = row[5]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(unquote)
                .collect::<Vec<_>>();
            for key in ["_CIPRank", "_CIPCode", "_ChiralityPossible"] {
                assert_eq!(
                    atom.is_prop_computed(key).unwrap(),
                    source_keys.iter().any(|k| k == key)
                );
            }
        }
        result
    }
    fn discrepancy<T: PartialEq + std::fmt::Debug>(
        label: &str,
        field: &str,
        actual: &T,
        expected: &T,
        failures: &mut Vec<String>,
    ) {
        if actual != expected {
            failures.push(format!("{label} {field}: {actual:?} != {expected:?}"));
        }
    }

    #[test]
    fn d2_state_owner_native_24_calls_and_two_typed_errors() {
        let fixture: Fixture =
            include!("../../../testdata/depict_2d/expected/rdkit/stereo_property_owners_core.rs");
        assert_eq!(
            fixture.source_pin,
            "351f8f378f8ad6bbd517980c38896e66bf907af8c"
        );
        let mut attempted = 0;
        let mut unexpected_errors = 0;
        let mut preserved = 0;
        let mut inner_calls = 0;
        let mut failures = Vec::new();
        for (label, cell) in fixture.cells.iter().filter(|(k, _)| k.contains("_OWNER_")) {
            let flags = label.rsplit('_').next().unwrap().as_bytes();
            let input = topology(&fixture.states[cell.states["OUTER_PRE"]]);
            let expected = topology(&fixture.states[cell.states["OUTER_POST"]]);
            let valence = crate::assign_valence_with_options_for_topology(
                &input,
                crate::ValenceModel::RdkitLike,
                false,
            )
            .unwrap();
            let rings =
                crate::symmetrized_sssr(&input, &crate::RingSearchParams::default()).unwrap();
            for repeat in 0..2 {
                let before = (input.clone(), valence.clone(), rings.clone());
                OBSERVATIONS.with(|slot| *slot.borrow_mut() = Some(Vec::new()));
                let result = assign_legacy_stereochemistry_with_flags(
                    input.clone(),
                    &valence,
                    &rings,
                    flags[0] == b'1',
                    flags[1] == b'1',
                );
                attempted += 1;
                let observed = OBSERVATIONS.with(|slot| slot.borrow_mut().take().unwrap());
                inner_calls += observed.len();
                let current = (input.clone(), valence.clone(), rings.clone());
                if current == before {
                    preserved += 1;
                }
                discrepancy(
                    label,
                    "whole-input preservation",
                    &current,
                    &before,
                    &mut failures,
                );
                match result {
                    Ok(output) => discrepancy(
                        label,
                        "whole represented output",
                        &output,
                        &expected,
                        &mut failures,
                    ),
                    Err(error) => {
                        unexpected_errors += 1;
                        failures.push(format!("{label}/{repeat}: {error:?}"));
                    }
                }
                discrepancy(
                    label,
                    "actual inner call count",
                    &observed.len(),
                    &cell.native_inner_calls,
                    &mut failures,
                );
                for (index, observation) in observed.iter().enumerate() {
                    if let Some(source) = cell.returns.get(index) {
                        discrepancy(
                            label,
                            "source ordinal",
                            &index,
                            &source.ordinal,
                            &mut failures,
                        );
                        discrepancy(
                            label,
                            "forwarded flag",
                            &(flags[1] == b'1'),
                            &source.possible,
                            &mut failures,
                        );
                        discrepancy(
                            label,
                            "actual returned flags",
                            &observation.flags,
                            &(source.has_unassigned, source.changed),
                            &mut failures,
                        );
                        discrepancy(
                            label,
                            "actual returned ranks",
                            &observation.ranks,
                            &source.ranks,
                            &mut failures,
                        );
                        discrepancy(
                            label,
                            "whole inner PRE",
                            &observation.before,
                            &topology(&fixture.states[cell.states[&format!("INNER_PRE_{index}")]]),
                            &mut failures,
                        );
                        discrepancy(
                            label,
                            "whole inner POST",
                            &observation.after,
                            &topology(&fixture.states[cell.states[&format!("INNER_POST_{index}")]]),
                            &mut failures,
                        );
                    } else {
                        failures.push(format!("{label} unexpected inner call {index}"));
                    }
                }
            }
        }
        let cell = &fixture.cells["2508_OWNER_11"];
        let valid = topology(&fixture.states[cell.states["OUTER_PRE"]]);
        let valid_valence = crate::assign_valence_with_options_for_topology(
            &valid,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        let rings = crate::symmetrized_sssr(&valid, &crate::RingSearchParams::default()).unwrap();
        let mut errors = 0;
        let mut error_preserved = 0;
        for topology_error in [true, false] {
            let mut input = valid.clone();
            let mut valence = valid_valence.clone();
            if topology_error {
                input.atoms.pop();
            } else {
                valence.explicit_valence.clear();
            }
            let before = (input.clone(), valence.clone(), rings.clone());
            let result = assign_legacy_stereochemistry_with_flags(
                input.clone(),
                &valence,
                &rings,
                false,
                true,
            );
            errors += 1;
            let typed = if topology_error {
                matches!(result, Err(LegacyStereoError::InvalidTopology(_)))
            } else {
                matches!(
                    result,
                    Err(LegacyStereoError::CipRank(CipRankError::ValenceRowCount {
                        field: "explicit_valence",
                        ..
                    }))
                )
            };
            if !typed {
                failures.push(format!("typed control {topology_error}: {result:?}"));
            }
            let after = (input, valence, rings.clone());
            if after == before {
                error_preserved += 1;
            }
            discrepancy(
                "typed error",
                "whole-input preservation",
                &after,
                &before,
                &mut failures,
            );
        }
        eprintln!(
            "STATE_OWNER attempted={attempted}/24 unexpected_errors={unexpected_errors} preservation={preserved}/24 inner_calls={inner_calls}/30 typed_controls={errors}/2 error_preservation={error_preserved}/2 discrepancies={}",
            failures.len()
        );
        assert_eq!(
            (attempted, errors, preserved, error_preserved),
            (24, 2, 24, 2)
        );
        assert!(failures.is_empty(), "{}", failures.join("\n"));
    }
}

#[cfg(test)]
mod source590_borrow_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondId, BondSpec, Element, PropertyValue};
    fn graph() -> TopologyBlock {
        TopologyBlock::try_from_parts(
            (0..2)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn source590_borrowed_legacy_is_same_owner_as_owning_entry_without_state_copy() {
        let g = graph();
        let v = crate::assign_valence_with_options_for_topology(
            &g,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        for clean in [false, true] {
            for flag in [false, true] {
                let r = crate::symmetrized_sssr(&g, &crate::RingSearchParams::default()).unwrap();
                let expected =
                    assign_legacy_stereochemistry_with_assignments(g.clone(), &v, &r, clean, flag)
                        .unwrap();
                let mut actual = g.clone();
                let atoms = actual.atoms.as_ptr();
                let bonds = actual.bonds.as_ptr();
                let mut update = None;
                assign_legacy_stereochemistry_source(
                    &mut actual,
                    &v,
                    &r,
                    None,
                    clean,
                    flag,
                    &mut update,
                )
                .unwrap();
                assert_eq!(actual, expected.topology);
                assert_eq!(update, expected.ring_update);
                assert_eq!(actual.atoms.as_ptr(), atoms);
                assert_eq!(actual.bonds.as_ptr(), bonds);
            }
        }
    }
    #[test]
    fn source590_borrowed_legacy_error_preserves_prior_property_clear_and_ring_effect() {
        let mut g = graph();
        let v = crate::assign_valence_with_options_for_topology(
            &g,
            crate::ValenceModel::RdkitLike,
            false,
        )
        .unwrap();
        g.atoms[0].set_prop("_CIPCode", "old").unwrap();
        g.atoms[1]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let r = RingInfo::new(crate::RingFindType::OtherOrUnknown, 2, 1);
        let mut update = None;
        assert!(
            assign_legacy_stereochemistry_source(&mut g, &v, &r, None, true, false, &mut update)
                .is_err()
        );
        assert_eq!(g.atoms[0].prop("_CIPCode"), None);
        assert_eq!(
            g.atoms[1].prop("__computedProps"),
            Some(&PropertyValue::Bool(false))
        );
        assert_eq!(
            update.as_ref().unwrap().find_type(),
            crate::RingFindType::Fast
        );
        assert_eq!(r.find_type(), crate::RingFindType::OtherOrUnknown);
    }
}
