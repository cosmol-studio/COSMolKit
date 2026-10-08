use crate::materialize::{
    ProductBuilder, ReactantProductMapping, int_prop, invariant, inversion_flag, uint_prop,
};
use crate::{ReactionCoordinateSelection, ReactionInput, ReactionProductError, ReactionRole};
use cosmolkit_core::{
    count_swaps_to_interconvert, find_double_bond_stereo_atoms_with_rank_reader,
    invert_atom_chirality,
};
use cosmolkit_model::{AtomId, BondId, QueryGraph, StereoGroup};
use cosmolkit_search::{SearchTarget, atom_total_degree_with_context, build_query_match_context};
use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag};

struct StereoBondEndCap {
    anchor: usize,
    non_anchor: Option<usize>,
}
impl StereoBondEndCap {
    fn new(
        input: &ReactionInput<'_>,
        atom: usize,
        other: usize,
        anchor: usize,
    ) -> Result<Self, ReactionProductError> {
        // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: StereoBondEndCap::StereoBondEndCap
        // RDKit❗✔️:   StereoBondEndCap(const ROMol &mol, const Atom *atom,
        // RDKit❗✔️:                    const Atom *otherDblBndAtom, const unsigned stereoAtomIdx)
        // RDKit❗✔️:       : m_anchor(stereoAtomIdx) {
        // RDKit❗✔️:     PRECONDITION(atom, "no atom");
        // RDKit❗✔️:     PRECONDITION(otherDblBndAtom, "no atom");
        // RDKit❗✔️:     PRECONDITION(atom->getTotalDegree() <= 3,
        // RDKit❗✔️:                  "Stereo Bond extremes must have less than four neighbors");
        // RDKit❗✔️:
        // RDKit❗✔️:     const auto nbrIdxItr = mol.getAtomNeighbors(atom);
        // RDKit❗✔️:     const unsigned otherIdx = otherDblBndAtom->getIdx();
        // RDKit❗✔️:
        // RDKit❗✔️:     auto isNonAnchor = [otherIdx, stereoAtomIdx](const unsigned &nbrIdx) {
        // RDKit❗✔️:       return nbrIdx != otherIdx && nbrIdx != stereoAtomIdx;
        // RDKit❗✔️:     };
        // RDKit❗✔️:
        // RDKit❗✔️:     auto nonAnchorItr =
        // RDKit❗✔️:         std::find_if(nbrIdxItr.first, nbrIdxItr.second, isNonAnchor);
        // RDKit❗✔️:     if (nonAnchorItr != nbrIdxItr.second) {
        // RDKit❗✔️:       mp_nonAnchor = mol.getAtomWithIdx(*nonAnchorItr);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // END RDKIT COMPLETE CPP FUNCTION
        // RDKit❗✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
        // RDKit❗✔️:   PRECONDITION(at, "no atom");
        // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
        // RDKit❗✔️:                "atom not associated with this molecule");
        // RDKit❗✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
        // RDKit❗✔️: }
        // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
        // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
        // RDKit❗✔️:
        // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
        // RDKit❗✔️:   const auto res = d_graph[vd];
        // RDKit❗✔️:
        // RDKit❗✔️:   POSTCONDITION(res, "");
        // RDKit❗✔️:   return res;
        // RDKit❗✔️: }
        // Constant-cost borrowed getter context followed by the source's single
        // O(degree) find_if; no graph copies, candidate buffers or sorting.
        let anchor = crate::materialize::source_u32("stereo anchor", anchor)? as usize;
        let source_atom = input
            .topology
            .atoms
            .get(atom)
            .ok_or_else(|| invariant("StereoBondEndCap", "no atom", Some(atom), None, None))?;
        let other_atom = input.topology.atoms.get(other).ok_or_else(|| {
            invariant(
                "StereoBondEndCap",
                "no other double-bond atom",
                Some(other),
                None,
                None,
            )
        })?;
        let target = SearchTarget::new(
            input.topology,
            input.coordinates,
            &input.topology.stereo_groups,
            input.rings,
            input.valence,
        );
        let context = build_query_match_context(&target);
        let source_index = source_atom.id().index();
        let total_degree =
            atom_total_degree_with_context(source_atom, &context).map_err(|source| {
                ReactionProductError::StereoGetter {
                    atom: source_atom.id(),
                    source,
                }
            })?;
        if total_degree > 3 {
            return Err(invariant(
                "StereoBondEndCap",
                "Stereo Bond extremes must have less than four neighbors",
                Some(source_index),
                None,
                None,
            ));
        }
        let neighbors = input
            .topology
            .adjacency
            .try_neighbors_of(source_index)
            .ok_or_else(|| {
                invariant(
                    "StereoBondEndCap",
                    "atom adjacency row missing",
                    Some(source_index),
                    None,
                    None,
                )
            })?;
        let other_index = other_atom.id().index();
        let non_anchor = if let Some(neighbor) = neighbors
            .iter()
            .find(|n| n.atom_index != other_index && n.atom_index != anchor)
        {
            let selected = input
                .topology
                .atoms
                .get(neighbor.atom_index)
                .ok_or_else(|| {
                    invariant(
                        "StereoBondEndCap",
                        "non-anchor atom row out of range",
                        Some(neighbor.atom_index),
                        None,
                        None,
                    )
                })?;
            // Native stores the pointer and getNonAnchorIdx reads its field.
            // The private caller retains this immutable input for both calls.
            Some(selected.id().index())
        } else {
            None
        };
        Ok(Self { anchor, non_anchor })
    }

    fn has_non_anchor(&self) -> bool {
        // RDKit❗✔️:   bool hasNonAnchor() const { return mp_nonAnchor != nullptr; }
        self.non_anchor.is_some()
    }

    fn anchor_idx(&self) -> usize {
        // RDKit❗✔️:   unsigned getAnchorIdx() const { return m_anchor; }
        self.anchor
    }

    fn non_anchor_idx(&self) -> Result<usize, ReactionProductError> {
        // RDKit❗✔️:   unsigned getNonAnchorIdx() const { return mp_nonAnchor->getIdx(); }
        // Native dereference requires a present non-anchor. Null is undefined,
        // not a zero-index fallback; retain an explicit structural error.
        self.non_anchor.ok_or_else(|| {
            invariant(
                "StereoBondEndCap::getNonAnchorIdx",
                "null non-anchor pointer",
                None,
                None,
                None,
            )
        })
    }

    fn product_candidates<'a>(
        &self,
        mapping: &'a ReactantProductMapping,
    ) -> Result<(&'a [usize], bool), ReactionProductError> {
        // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: getProductAnchorCandidates
        // RDKit❗🔝:   std::pair<UINT_VECT, bool> getProductAnchorCandidates(
        // RDKit❗🔝:       ReactantProductAtomMapping *mapping) {
        // RDKit❗🔝:     auto &react2Prod = mapping->reactProdAtomMap;
        // RDKit❗🔝:
        // RDKit❗🔝:     bool swapStereo = false;
        // RDKit❗🔝:     auto newAnchorMatches = react2Prod.find(getAnchorIdx());
        // RDKit❗🔝:     if (newAnchorMatches != react2Prod.end()) {
        // RDKit❗🔝:       // The corresponding StereoAtom exists in the product
        // RDKit❗🔝:       return {newAnchorMatches->second, swapStereo};
        // RDKit❗🔝:
        // RDKit❗🔝:     } else if (hasNonAnchor()) {
        // RDKit❗🔝:       // The non-StereoAtom neighbor exists in the product
        // RDKit❗🔝:       newAnchorMatches = react2Prod.find(getNonAnchorIdx());
        // RDKit❗🔝:       if (newAnchorMatches != react2Prod.end()) {
        // RDKit❗🔝:         swapStereo = true;
        // RDKit❗🔝:         return {newAnchorMatches->second, swapStereo};
        // RDKit❗🔝:       }
        // RDKit❗🔝:     }
        // RDKit❗🔝:     // None of the neighbors survived the reaction
        // RDKit❗🔝:     return {{}, swapStereo};
        // RDKit❗🔝:   }
        // END RDKIT COMPLETE CPP FUNCTION
        // Cost improvement: the source copies its selected UINT_VECT. This
        // private flow consumes a borrowed slice while mapping is immutable,
        // preserving encounter order and present-empty key distinctions without
        // that vector allocation/copy or extra temporary ownership.
        if let Some(matches) = mapping.reactant_to_product.get(&self.anchor_idx()) {
            return Ok((matches, false));
        }
        if self.has_non_anchor() {
            if let Some(matches) = mapping.reactant_to_product.get(&self.non_anchor_idx()?) {
                return Ok((matches, true));
            }
        }
        Ok((&[], false))
    }
}

fn anchor_index(
    product: &ProductBuilder,
    atom: usize,
    candidates: &[usize],
) -> Result<usize, ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: reactProdMapAnchorIdx
    // RDKit❗✔️: unsigned reactProdMapAnchorIdx(Atom *atom, const RDKit::UINT_VECT &pMatches) {
    // RDKit❗✔️:   PRECONDITION(atom, "no atom");
    // RDKit❗✔️:   if (pMatches.size() == 1) {
    // RDKit❗✔️:     return pMatches[0];
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   const auto &pMol = atom->getOwningMol();
    // RDKit❗✔️:   const unsigned atomIdx = atom->getIdx();
    // RDKit❗✔️:
    // RDKit❗✔️:   auto areAtomsBonded = [&pMol, atomIdx](const unsigned &pAnchor) {
    // RDKit❗✔️:     return pMol.getBondBetweenAtoms(atomIdx, pAnchor) != nullptr;
    // RDKit❗✔️:   };
    // RDKit❗✔️:
    // RDKit❗✔️:   auto match = std::find_if(pMatches.begin(), pMatches.end(), areAtomsBonded);
    // RDKit❗✔️:
    // RDKit❗✔️:   CHECK_INVARIANT(match != pMatches.end(), "match not found");
    // RDKit❗✔️:   return *match;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Source direct singleton access precedes its O(candidates * degree) scan.
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    let source_atom = product.topology.atoms.get(atom).ok_or_else(|| {
        invariant(
            "reactProdMapAnchorIdx",
            "atom row out of range",
            None,
            Some(atom),
            None,
        )
    })?;
    let source_index = source_atom.id().index();
    for &candidate in candidates {
        // Use the actual source Atom::getIdx and canonical checked source edge
        // lookup. Every lookup occurs only when std::find_if would reach it;
        // a first match prevents inspection of later invalid candidate indices.
        if cosmolkit_model::source_bond_between_atoms(
            product.topology.atoms.len(),
            AtomId::new(source_index),
            AtomId::new(candidate),
            || product.neighbors[source_index].iter(),
        )?
        .is_some()
        {
            return Ok(candidate);
        }
    }
    Err(invariant(
        "reactProdMapAnchorIdx",
        "match not found",
        None,
        Some(source_index),
        None,
    ))
}

fn forward_bond_stereo(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &mut ReactantProductMapping,
    product_bond: BondId,
    reactant_bond: BondId,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: forwardReactantBondStereo
    // RDKit❗❌: void forwardReactantBondStereo(ReactantProductAtomMapping *mapping, Bond *pBond,
    // RDKit❗❌:                                const ROMol &reactant, const Bond *rBond) {
    // RDKit❗❌:   PRECONDITION(mapping, "no mapping");
    // RDKit❗❌:   PRECONDITION(pBond, "no bond");
    // RDKit❗❌:   PRECONDITION(rBond, "no bond");
    // RDKit❗❌:   PRECONDITION(rBond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit❗❌:                "bond in reactant must have defined stereo");
    // RDKit❗❌:
    // RDKit❗❌:   auto &prod2React = mapping->prodReactAtomMap;
    // RDKit❗❌:
    // RDKit❗❌:   const Atom *rStart = rBond->getBeginAtom();
    // RDKit❗❌:   const Atom *rEnd = rBond->getEndAtom();
    // RDKit❗❌:   const auto rStereoAtoms = Chirality::findStereoAtoms(rBond);
    // RDKit❗❌:   if (rStereoAtoms.size() != 2) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog)
    // RDKit❗❌:         << "WARNING: neither stereo atoms nor CIP codes found for double bond. "
    // RDKit❗❌:            "Stereochemistry info will not be propagated to product."
    // RDKit❗❌:         << std::endl;
    // RDKit❗❌:     pBond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   StereoBondEndCap start(reactant, rStart, rEnd, rStereoAtoms[0]);
    // RDKit❗❌:   StereoBondEndCap end(reactant, rEnd, rStart, rStereoAtoms[1]);
    // RDKit❗❌:
    // RDKit❗❌:   // The bond might be matched backwards in the reaction
    // RDKit❗❌:   if (prod2React[pBond->getBeginAtomIdx()] == rEnd->getIdx()) {
    // RDKit❗❌:     std::swap(start, end);
    // RDKit❗❌:   } else if (prod2React[pBond->getBeginAtomIdx()] != rStart->getIdx()) {
    // RDKit❗❌:     throw std::logic_error("Reactant and Product bond ends do not match");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   /**
    // RDKit❗❌:    *  The reactants stereo can be transmitted in three similar ways:
    // RDKit❗❌:    *
    // RDKit❗❌:    * 1. Survival of both stereoatoms: direct forwarding happens, i.e.,
    // RDKit❗❌:    *
    // RDKit❗❌:    *    C/C=C/[Br] in reaction [C:1]=[C:2]>>[Si:1]=[C:2]:
    // RDKit❗❌:    *
    // RDKit❗❌:    *    C/C=C/[Br] >> C/Si=C/[Br], C/C=Si/[Br] (2 product sets)
    // RDKit❗❌:    *
    // RDKit❗❌:    *    Both stereoatoms exist unaltered in both product sets, so we can forward
    // RDKit❗❌:    *    the same bond stereochemistry (trans) and set the stereoatoms in the
    // RDKit❗❌:    *    product to the mapped indexes of the stereoatoms in the reactant.
    // RDKit❗❌:    *
    // RDKit❗❌:    * 2. Survival of both anti-stereoatoms: as this pair is symmetric to the
    // RDKit❗❌:    *    stereoatoms, direct forwarding also happens in this case, i.e.,
    // RDKit❗❌:    *
    // RDKit❗❌:    *    Cl/C(C)=C(/Br)F in reaction
    // RDKit❗❌:    *        [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3].[Cl:4]:
    // RDKit❗❌:    *      Cl/C(C)=C(/Br)F >> C/C=C/F + Br + Cl
    // RDKit❗❌:    *
    // RDKit❗❌:    *    Both stereoatoms in the reactant are split from the molecule,
    // RDKit❗❌:    *    but the anti-stereoatoms remain in it. Since these have symmetrical
    // RDKit❗❌:    *    orientation to the stereoatoms, we can use these (their mapped
    // RDKit❗❌:    *    equivalents) as stereoatoms in the product and use the same
    // RDKit❗❌:    *    stereochemistry label (trans).
    // RDKit❗❌:    *
    // RDKit❗❌:    * 3. Survival of a mixed pair stereoatom-anti-stereoatom: such a pair
    // RDKit❗❌:    *    defines the opposite stereochemistry to the one labeled on the
    // RDKit❗❌:    *    reactant, but it is also valid, as long ase we use the properly mapped
    // RDKit❗❌:    *    indexes:
    // RDKit❗❌:    *
    // RDKit❗❌:    *    Cl/C(C)=C(/Br)F in reaction [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3]:
    // RDKit❗❌:    *
    // RDKit❗❌:    *        Cl/C(C)=C(/Br)F >> C/C=C/F + Br
    // RDKit❗❌:    *
    // RDKit❗❌:    *    In this case, one of the stereoatoms is conserved, and the other one is
    // RDKit❗❌:    *    switched to the other neighbor at the same end of the bond as the
    // RDKit❗❌:    *    non-conserved stereoatom. Since the reference changed, the
    // RDKit❗❌:    *    stereochemistry label needs to be flipped too: in this case, the
    // RDKit❗❌:    *    reactant was trans, and the product will be cis.
    // RDKit❗❌:    *
    // RDKit❗❌:    *    Reaction [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Cl:4] would have the same
    // RDKit❗❌:    *    effect, with the only difference that the non-conserved stereoatom would
    // RDKit❗❌:    *    be the one at the opposite end of the reactant.
    // RDKit❗❌:    */
    // RDKit❗❌:   auto pStartAnchorCandidates = start.getProductAnchorCandidates(mapping);
    // RDKit❗❌:   auto pEndAnchorCandidates = end.getProductAnchorCandidates(mapping);
    // RDKit❗❌:
    // RDKit❗❌:   // The reaction has invalidated the reactant's stereochemistry
    // RDKit❗❌:   if (pStartAnchorCandidates.first.empty() ||
    // RDKit❗❌:       pEndAnchorCandidates.first.empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   unsigned pStartAnchorIdx = reactProdMapAnchorIdx(
    // RDKit❗❌:       pBond->getBeginAtom(), pStartAnchorCandidates.first);
    // RDKit❗❌:   unsigned pEndAnchorIdx =
    // RDKit❗❌:       reactProdMapAnchorIdx(pBond->getEndAtom(), pEndAnchorCandidates.first);
    // RDKit❗❌:
    // RDKit❗❌:   const ROMol &m = pBond->getOwningMol();
    // RDKit❗❌:   if (m.getBondBetweenAtoms(pBond->getBeginAtomIdx(), pStartAnchorIdx) ==
    // RDKit❗❌:           nullptr ||
    // RDKit❗❌:       m.getBondBetweenAtoms(pBond->getEndAtomIdx(), pEndAnchorIdx) == nullptr) {
    // RDKit❗❌:     BOOST_LOG(rdWarningLog) << "stereo atoms in input cannot be mapped to "
    // RDKit❗❌:                                "output (atoms are no longer bonded)\n";
    // RDKit❗❌:   } else {
    // RDKit❗❌:     pBond->setStereoAtoms(pStartAnchorIdx, pEndAnchorIdx);
    // RDKit❗❌:     bool flipStereo =
    // RDKit❗❌:         (pStartAnchorCandidates.second + pEndAnchorCandidates.second) % 2;
    // RDKit❗❌:
    // RDKit❗❌:     if (rBond->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗❌:         rBond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit❗❌:       if (flipStereo) {
    // RDKit❗❌:         pBond->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         pBond->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       if (flipStereo) {
    // RDKit❗❌:         pBond->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗❌:       } else {
    // RDKit❗❌:         pBond->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗❌: void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
    // RDKit❗❌:   PRECONDITION(
    // RDKit❗❌:       getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
    // RDKit❗❌:       "bgnIdx not connected to begin atom of bond");
    // RDKit❗❌:   PRECONDITION(
    // RDKit❗❌:       getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
    // RDKit❗❌:       "endIdx not connected to end atom of bond");
    // RDKit❗❌:
    // RDKit❗❌:   auto &atoms = getStereoAtoms();
    // RDKit❗❌:   atoms.clear();
    // RDKit❗❌:   atoms.push_back(bgnIdx);
    // RDKit❗❌:   atoms.push_back(endIdx);
    // RDKit❗❌: }
    // Source body scans only reached endpoints/candidates. The sole CORE
    // findStereoAtoms currently additionally validates O(V+E) detached state;
    // this transitive cost/invalid-state difference remains explicit, not hidden
    // by a duplicated local stereo or rank algorithm.
    product
        .topology
        .bonds
        .get(product_bond.index())
        .ok_or_else(|| {
            invariant(
                "forwardReactantBondStereo",
                "no product bond",
                None,
                None,
                Some(product_bond),
            )
        })?;
    let r = input
        .topology
        .bonds
        .get(reactant_bond.index())
        .ok_or_else(|| {
            invariant(
                "forwardReactantBondStereo",
                "no reactant bond",
                None,
                None,
                Some(reactant_bond),
            )
        })?;
    if r.stereo().rdkit_code() <= BondStereo::Any.rdkit_code() {
        return Err(invariant(
            "forwardReactantBondStereo",
            "bond in reactant must have defined stereo",
            None,
            None,
            Some(reactant_bond),
        ));
    }
    let r_start = input.topology.atoms.get(r.begin().index()).ok_or_else(|| {
        invariant(
            "forwardReactantBondStereo",
            "reactant begin atom missing",
            Some(r.begin().index()),
            None,
            Some(reactant_bond),
        )
    })?;
    let r_end = input.topology.atoms.get(r.end().index()).ok_or_else(|| {
        invariant(
            "forwardReactantBondStereo",
            "reactant end atom missing",
            Some(r.end().index()),
            None,
            Some(reactant_bond),
        )
    })?;
    let found = find_double_bond_stereo_atoms_with_rank_reader(
        input.topology,
        reactant_bond,
        |atom| uint_prop(&input.topology.atoms[atom.index()], "_CIPRank"),
        |warning| eprintln!("{warning}"),
    )?;
    let Some([a, b]) = found.atoms else {
        eprintln!(
            "WARNING: neither stereo atoms nor CIP codes found for double bond. Stereochemistry info will not be propagated to product."
        );
        product.topology.bonds[product_bond.index()].set_stereo(BondStereo::None)?;
        return Ok(());
    };
    let mut start = StereoBondEndCap::new(input, r.begin().index(), r.end().index(), a.index())?;
    let mut end = StereoBondEndCap::new(input, r.end().index(), r.begin().index(), b.index())?;
    let p = &product.topology.bonds[product_bond.index()];
    let begin = p.begin().index();
    let finish = p.end().index();
    // std::map::operator[] value-initializes a missing unsigned entry before
    // either comparison. Keep that observable mutation even on later failure.
    let source_begin = *mapping.product_to_reactant.entry(begin).or_insert(0);
    if source_begin == r_end.id().index() {
        std::mem::swap(&mut start, &mut end);
    } else if source_begin != r_start.id().index() {
        return Err(invariant(
            "forwardReactantBondStereo",
            "Reactant and Product bond ends do not match",
            Some(source_begin),
            Some(begin),
            Some(product_bond),
        ));
    }
    let (start_matches, start_swap) = start.product_candidates(mapping)?;
    let (end_matches, end_swap) = end.product_candidates(mapping)?;
    if start_matches.is_empty() || end_matches.is_empty() {
        return Ok(());
    }
    product.topology.atoms.get(begin).ok_or_else(|| {
        invariant(
            "forwardReactantBondStereo",
            "product begin atom missing",
            None,
            Some(begin),
            Some(product_bond),
        )
    })?;
    let start_anchor = anchor_index(product, begin, start_matches)?;
    product.topology.atoms.get(finish).ok_or_else(|| {
        invariant(
            "forwardReactantBondStereo",
            "product end atom missing",
            None,
            Some(finish),
            Some(product_bond),
        )
    })?;
    let end_anchor = anchor_index(product, finish, end_matches)?;
    // Preserve the source OR short circuit: disconnected start suppresses the
    // end lookup, including its source range check; neither warning path writes.
    if cosmolkit_model::source_bond_between_atoms(
        product.topology.atoms.len(),
        AtomId::new(begin),
        AtomId::new(start_anchor),
        || product.neighbors[begin].iter(),
    )?
    .is_none()
        || cosmolkit_model::source_bond_between_atoms(
            product.topology.atoms.len(),
            AtomId::new(finish),
            AtomId::new(end_anchor),
            || product.neighbors[finish].iter(),
        )?
        .is_none()
    {
        eprintln!("stereo atoms in input cannot be mapped to output (atoms are no longer bonded)");
        return Ok(());
    }
    // Source setStereoAtoms repeats the connectedness checks before any write;
    // this immutable graph cannot change between them and the two above.
    let p = &mut product.topology.bonds[product_bond.index()];
    p.set_stereo_atoms(Some([AtomId::new(start_anchor), AtomId::new(end_anchor)]));
    let flip = start_swap != end_swap;
    let cis = matches!(r.stereo(), BondStereo::Cis | BondStereo::Z) != flip;
    p.set_stereo(if cis {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    })?;
    Ok(())
}

fn other_atom_index(
    bond: &cosmolkit_model::Bond,
    atom: AtomId,
) -> Result<AtomId, ReactionProductError> {
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // Shared reaction source helper: constant endpoint comparisons, no graph
    // work or allocation; nonincident values retain the source postcondition.
    if bond.begin() == atom {
        Ok(bond.end())
    } else if bond.end() == atom {
        Ok(bond.begin())
    } else {
        Err(invariant(
            "Bond::getOtherAtomIdx",
            "bad index",
            None,
            Some(atom.index()),
            Some(bond.id()),
        ))
    }
}

fn translate_directions(
    product: &mut ProductBuilder,
    bond: BondId,
    start: BondId,
    end: BondId,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: translateProductStereoBondDirections
    // RDKit❗✔️: void translateProductStereoBondDirections(Bond *pBond, const Bond *start,
    // RDKit❗✔️:                                           const Bond *end) {
    // RDKit❗✔️:   PRECONDITION(pBond, "no bond");
    // RDKit❗✔️:   PRECONDITION(start && end && Chirality::hasStereoBondDir(start) &&
    // RDKit❗✔️:                    Chirality::hasStereoBondDir(end),
    // RDKit❗✔️:                "Both neighboring bonds must have bond directions");
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned pStartAnchorIdx = start->getOtherAtomIdx(pBond->getBeginAtomIdx());
    // RDKit❗✔️:   unsigned pEndAnchorIdx = end->getOtherAtomIdx(pBond->getEndAtomIdx());
    // RDKit❗✔️:
    // RDKit❗✔️:   pBond->setStereoAtoms(pStartAnchorIdx, pEndAnchorIdx);
    // RDKit❗✔️:
    // RDKit❗✔️:   bool sameDir = start->getBondDir() == end->getBondDir();
    // RDKit❗✔️:
    // RDKit❗✔️:   if (start->getBeginAtom() == pBond->getBeginAtom()) {
    // RDKit❗✔️:     sameDir = !sameDir;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (end->getBeginAtom() != pBond->getEndAtom()) {
    // RDKit❗✔️:     sameDir = !sameDir;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (sameDir) {
    // RDKit❗✔️:     pBond->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     pBond->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗✔️: unsigned int Bond::getOtherAtomIdx(const unsigned int thisIdx) const {
    // RDKit❗✔️:   if (d_beginAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_endAtomIdx;
    // RDKit❗✔️:   } else if (d_endAtomIdx == thisIdx) {
    // RDKit❗✔️:     return d_beginAtomIdx;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // This "precondition" would check exactly the same that is checked
    // RDKit❗✔️:   // above, but no need to be redundant, so just throw.
    // RDKit❗✔️:   POSTCONDITION(false, "bad index");
    // RDKit❗✔️: }
    // RDKit❗✔️: void Bond::setStereoAtoms(unsigned int bgnIdx, unsigned int endIdx) {
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       getOwningMol().getBondBetweenAtoms(getBeginAtomIdx(), bgnIdx) != nullptr,
    // RDKit❗✔️:       "bgnIdx not connected to begin atom of bond");
    // RDKit❗✔️:   PRECONDITION(
    // RDKit❗✔️:       getOwningMol().getBondBetweenAtoms(getEndAtomIdx(), endIdx) != nullptr,
    // RDKit❗✔️:       "endIdx not connected to end atom of bond");
    // RDKit❗✔️:
    // RDKit❗✔️:   auto &atoms = getStereoAtoms();
    // RDKit❗✔️:   atoms.clear();
    // RDKit❗✔️:   atoms.push_back(bgnIdx);
    // RDKit❗✔️:   atoms.push_back(endIdx);
    // RDKit❗✔️: }
    // Same source constant scalar work and two O(degree) connectedness checks;
    // no allocated collections, cloned bond/graph state, or inferred directions.
    let p = product.topology.bonds.get(bond.index()).ok_or_else(|| {
        invariant(
            "translateProductStereoBondDirections",
            "no product bond",
            None,
            None,
            Some(bond),
        )
    })?;
    let begin = p.begin();
    let finish = p.end();
    let start = product.topology.bonds.get(start.index()).ok_or_else(|| {
        invariant(
            "translateProductStereoBondDirections",
            "Both neighboring bonds must have bond directions",
            None,
            None,
            Some(bond),
        )
    })?;
    let end = product.topology.bonds.get(end.index()).ok_or_else(|| {
        invariant(
            "translateProductStereoBondDirections",
            "Both neighboring bonds must have bond directions",
            None,
            None,
            Some(bond),
        )
    })?;
    if !cosmolkit_core::has_stereo_bond_direction(start.direction())
        || !cosmolkit_core::has_stereo_bond_direction(end.direction())
    {
        return Err(invariant(
            "translateProductStereoBondDirections",
            "Both neighboring bonds must have bond directions",
            None,
            None,
            Some(bond),
        ));
    }
    // This source helper is checked, rather than treating every non-begin input
    // as the end. Keep start's POSTCONDITION before evaluating end's helper.
    let start_anchor = other_atom_index(start, begin)?;
    let end_anchor = other_atom_index(end, finish)?;
    if cosmolkit_model::source_bond_between_atoms(
        product.topology.atoms.len(),
        begin,
        start_anchor,
        || product.neighbors[begin.index()].iter(),
    )?
    .is_none()
    {
        return Err(invariant(
            "Bond::setStereoAtoms",
            "bgnIdx not connected to begin atom of bond",
            None,
            Some(start_anchor.index()),
            Some(bond),
        ));
    }
    if cosmolkit_model::source_bond_between_atoms(
        product.topology.atoms.len(),
        finish,
        end_anchor,
        || product.neighbors[finish.index()].iter(),
    )?
    .is_none()
    {
        return Err(invariant(
            "Bond::setStereoAtoms",
            "endIdx not connected to end atom of bond",
            None,
            Some(end_anchor.index()),
            Some(bond),
        ));
    }
    let start_begin = start.begin();
    let end_begin = end.begin();
    let mut same = start.direction() == end.direction();
    // Source pointer equality compares rows in this sole owning product, not
    // arbitrary AtomId fields stored in those rows.
    let p = &mut product.topology.bonds[bond.index()];
    p.set_stereo_atoms(Some([start_anchor, end_anchor]));
    if start_begin == begin {
        same = !same;
    }
    if end_begin != finish {
        same = !same;
    }
    p.set_stereo(if same {
        BondStereo::Trans
    } else {
        BondStereo::Cis
    })?;
    Ok(())
}

fn directed_neighbor(
    product: &ProductBuilder,
    atom: usize,
) -> Result<Option<BondId>, ReactionProductError> {
    // RDKit❗✔️: bool hasStereoBondDir(const Bond *bond) {
    // RDKit❗✔️:   PRECONDITION(bond, "no bond");
    // RDKit❗✔️:   return bond->getBondDir() == Bond::BondDir::ENDDOWNRIGHT ||
    // RDKit❗✔️:          bond->getBondDir() == Bond::BondDir::ENDUPRIGHT;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: const Bond *getNeighboringDirectedBond(const ROMol &mol, const Atom *atom) {
    // RDKit❗✔️:   PRECONDITION(atom, "no atom");
    // RDKit❗✔️:   for (const auto &bondIdx :
    // RDKit❗✔️:        boost::make_iterator_range(mol.getAtomBonds(atom))) {
    // RDKit❗✔️:     const Bond *bond = mol[bondIdx];
    // RDKit❗✔️:
    // RDKit❗✔️:     if (bond->getBondType() != Bond::BondType::DOUBLE &&
    // RDKit❗✔️:         hasStereoBondDir(bond)) {
    // RDKit❗✔️:       return bond;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return nullptr;
    // RDKit❗✔️: }
    // RDKit❗✔️:
    let source_atom = product.topology.atoms.get(atom).ok_or_else(|| {
        invariant(
            "getNeighboringDirectedBond",
            "no atom",
            None,
            Some(atom),
            None,
        )
    })?;
    let source_index = source_atom.id().index();
    let neighbors = product.neighbors.get(source_index).ok_or_else(|| {
        invariant(
            "getNeighboringDirectedBond",
            "atom adjacency row missing",
            None,
            Some(source_index),
            None,
        )
    })?;
    cosmolkit_core::neighboring_directed_bond_from_incident(neighbors.iter().map(|neighbor| {
        product
            .topology
            .bonds
            .get(neighbor.bond.index())
            .ok_or_else(|| {
                invariant(
                    "getNeighboringDirectedBond",
                    "incident bond row missing",
                    None,
                    Some(source_index),
                    Some(neighbor.bond),
                )
            })
    }))
    .map(|bond| bond.map(cosmolkit_model::Bond::id))
}

pub(crate) fn update_stereo_bonds(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &mut ReactantProductMapping,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateStereoBonds
    // RDKit❗❌: void updateStereoBonds(RWMOL_SPTR product, const ROMol &reactant,
    // RDKit❗❌:                        ReactantProductAtomMapping *mapping) {
    // RDKit❗❌:   for (Bond *pBond : product->bonds()) {
    // RDKit❗❌:     // We are only interested in double bonds
    // RDKit❗❌:     if (pBond->getBondType() != Bond::BondType::DOUBLE) {
    // RDKit❗❌:       continue;
    // RDKit❗❌:     } else if (pBond->hasProp(_UnknownStereoRxnBond)) {
    // RDKit❗❌:       pBond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗❌:       pBond->clearProp(_UnknownStereoRxnBond);
    // RDKit❗❌:       continue;
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // Check if the reaction defined the stereo for the bond: SMARTS can only
    // RDKit❗❌:     // use bond directions for this, and both sides of the double bond must have
    // RDKit❗❌:     // them, else they will be ignored, as there is no reference to decide the
    // RDKit❗❌:     // stereo.
    // RDKit❗❌:     const auto *pBondStartDirBond =
    // RDKit❗❌:         Chirality::getNeighboringDirectedBond(*product, pBond->getBeginAtom());
    // RDKit❗❌:     const auto *pBondEndDirBond =
    // RDKit❗❌:         Chirality::getNeighboringDirectedBond(*product, pBond->getEndAtom());
    // RDKit❗❌:     if (pBondStartDirBond != nullptr && pBondEndDirBond != nullptr) {
    // RDKit❗❌:       translateProductStereoBondDirections(pBond, pBondStartDirBond,
    // RDKit❗❌:                                            pBondEndDirBond);
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // If the reaction did not specify the stereo, then we need to rely on the
    // RDKit❗❌:       // atom mapping and use the reactant's stereo.
    // RDKit❗❌:
    // RDKit❗❌:       // The atoms and the bond might have been added in the reaction
    // RDKit❗❌:       const auto begIdxItr =
    // RDKit❗❌:           mapping->prodReactAtomMap.find(pBond->getBeginAtomIdx());
    // RDKit❗❌:       if (begIdxItr == mapping->prodReactAtomMap.end()) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:       const auto endIdxItr =
    // RDKit❗❌:           mapping->prodReactAtomMap.find(pBond->getEndAtomIdx());
    // RDKit❗❌:       if (endIdxItr == mapping->prodReactAtomMap.end()) {
    // RDKit❗❌:         continue;
    // RDKit❗❌:       }
    // RDKit❗❌:
    // RDKit❗❌:       const Bond *rBond =
    // RDKit❗❌:           reactant.getBondBetweenAtoms(begIdxItr->second, endIdxItr->second);
    // RDKit❗❌:
    // RDKit❗❌:       if (rBond && rBond->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗❌:         // The bond might not have been present in the reactant, or its order
    // RDKit❗❌:         // might have changed
    // RDKit❗❌:         if (rBond->getStereo() > Bond::BondStereo::STEREOANY) {
    // RDKit❗❌:           // If the bond had stereo, forward it
    // RDKit❗❌:           forwardReactantBondStereo(mapping, pBond, reactant, rBond);
    // RDKit❗❌:         } else if (rBond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit❗❌:           pBond->setStereo(Bond::BondStereo::STEREOANY);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // No stereo: Bond::BondStereo::STEREONONE
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗❌:   void clearProp(const std::string_view key) const {
    // RDKit❗❌:     STR_VECT compLst;
    // RDKit❗❌:     if (getPropIfPresent(RDKit::detail::computedPropName, compLst)) {
    // RDKit❗❌:       auto svi = std::find(compLst.begin(), compLst.end(), key);
    // RDKit❗❌:       if (svi != compLst.end()) {
    // RDKit❗❌:         compLst.erase(svi);
    // RDKit❗❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     d_props.clearVal(key);
    // RDKit❗❌:   }
    // Source physical bond traversal and conditional local reads. Forwarding
    // retains the sole CORE findStereoAtoms whole-topology validation overhead;
    // detached property-map operations also differ from native dictionaries.
    for row in 0..product.topology.bonds.len() {
        let p = &product.topology.bonds[row];
        if p.order() != BondOrder::Double {
            continue;
        }
        if p.prop("_UnknownStereoRxnBond").is_some() {
            let p = &mut product.topology.bonds[row];
            p.set_stereo(BondStereo::None)?;
            p.clear_prop("_UnknownStereoRxnBond")?;
            continue;
        }
        let begin = p.begin().index();
        let end = p.end().index();
        let id = p.id();
        // Both source helper calls are evaluated, in order, even when the
        // first returns no directed bond. Each helper itself stops at first hit.
        let start = directed_neighbor(product, begin)?;
        let finish = directed_neighbor(product, end)?;
        if let (Some(start), Some(finish)) = (start, finish) {
            translate_directions(product, id, start, finish)?;
        } else {
            let Some(&r_begin) = mapping.product_to_reactant.get(&begin) else {
                continue;
            };
            let Some(&r_end) = mapping.product_to_reactant.get(&end) else {
                continue;
            };
            let neighbors = input.topology.adjacency.try_neighbors_of(r_begin);
            // Canonical source endpoint range checks run before the borrowed
            // row is consumed. A missing row is rejected explicitly below,
            // never accepted as a zero-degree source atom.
            let r_bond = cosmolkit_model::source_bond_between_atoms(
                input.topology.atoms.len(),
                AtomId::new(r_begin),
                AtomId::new(r_end),
                || neighbors.into_iter().flatten(),
            )?;
            neighbors.ok_or_else(|| {
                invariant(
                    "updateStereoBonds",
                    "reactant adjacency row missing",
                    Some(r_begin),
                    None,
                    None,
                )
            })?;
            let Some(r_bond) = r_bond else {
                continue;
            };
            let r = input.topology.bonds.get(r_bond.index()).ok_or_else(|| {
                invariant(
                    "updateStereoBonds",
                    "reactant bond row missing",
                    None,
                    None,
                    Some(r_bond),
                )
            })?;
            if r.order() != BondOrder::Double {
                continue;
            }
            if r.stereo().rdkit_code() > BondStereo::Any.rdkit_code() {
                forward_bond_stereo(product, input, mapping, id, r_bond)?;
            } else if r.stereo() == BondStereo::Any {
                product.topology.bonds[row].set_stereo(BondStereo::Any)?;
            }
        }
    }
    Ok(())
}

pub(crate) fn correct_matched_chirality(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &mut ReactantProductMapping,
    reactant_atom: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: checkAndCorrectChiralityOfMatchingAtomsInProduct
    // RDKit❗🔝: void checkAndCorrectChiralityOfMatchingAtomsInProduct(
    // RDKit❗🔝:     const ROMol &reactant, unsigned reactantAtomIdx, const Atom &reactantAtom,
    // RDKit❗🔝:     RWMOL_SPTR product, ReactantProductAtomMapping *mapping) {
    // RDKit❗🔝:   for (unsigned i = 0; i < mapping->reactProdAtomMap[reactantAtomIdx].size();
    // RDKit❗🔝:        i++) {
    // RDKit❗🔝:     unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
    // RDKit❗🔝:     Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗🔝:
    // RDKit❗🔝:     int inversionFlag = 0;
    // RDKit❗🔝:     productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗🔝:                                   inversionFlag);
    // RDKit❗🔝:     // if stereochemistry wasn't present in the reactant or if we're
    // RDKit❗🔝:     // either creating or destroying stereo we don't mess with this
    // RDKit❗🔝:     if (reactantAtom.getChiralTag() == Atom::CHI_UNSPECIFIED ||
    // RDKit❗🔝:         reactantAtom.getChiralTag() == Atom::CHI_OTHER || inversionFlag > 2) {
    // RDKit❗🔝:       continue;
    // RDKit❗🔝:     }
    // RDKit❗🔝:
    // RDKit❗🔝:     // we can only do something sensible here if the degree in the reactants
    // RDKit❗🔝:     // and products differs by at most one
    // RDKit❗🔝:     if (reactantAtom.getDegree() < 3 || productAtom->getDegree() < 3 ||
    // RDKit❗🔝:         std::abs(static_cast<int>(reactantAtom.getDegree()) -
    // RDKit❗🔝:                  static_cast<int>(productAtom->getDegree())) > 1) {
    // RDKit❗🔝:       continue;
    // RDKit❗🔝:     }
    // RDKit❗🔝:     unsigned int nUnknown = 0;
    // RDKit❗🔝:     // get the order of the bonds around the atom in the reactant:
    // RDKit❗🔝:     INT_LIST rOrder;
    // RDKit❗🔝:     for (const auto &nbri :
    // RDKit❗🔝:          boost::make_iterator_range(reactant.getAtomBonds(&reactantAtom))) {
    // RDKit❗🔝:       rOrder.push_back(reactant[nbri]->getIdx());
    // RDKit❗🔝:     }
    // RDKit❗🔝:     INT_LIST pOrder;
    // RDKit❗🔝:     for (const auto &nbri :
    // RDKit❗🔝:          boost::make_iterator_range(product->getAtomNeighbors(productAtom))) {
    // RDKit❗🔝:       if (mapping->prodReactAtomMap.find(nbri) ==
    // RDKit❗🔝:               mapping->prodReactAtomMap.end() ||
    // RDKit❗🔝:           !reactant.getBondBetweenAtoms(reactantAtom.getIdx(),
    // RDKit❗🔝:                                         mapping->prodReactAtomMap[nbri])) {
    // RDKit❗🔝:         ++nUnknown;
    // RDKit❗🔝:         // if there's more than one bond in the product that doesn't
    // RDKit❗🔝:         // correspond to anything in the reactant, we're also doomed
    // RDKit❗🔝:         if (nUnknown > 1) {
    // RDKit❗🔝:           break;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         // otherwise, add a -1 to the bond order that we'll fill in later
    // RDKit❗🔝:         pOrder.push_back(-1);
    // RDKit❗🔝:       } else {
    // RDKit❗🔝:         const Bond *rBond = reactant.getBondBetweenAtoms(
    // RDKit❗🔝:             reactantAtom.getIdx(), mapping->prodReactAtomMap[nbri]);
    // RDKit❗🔝:         CHECK_INVARIANT(rBond, "expected reactant bond not found");
    // RDKit❗🔝:         pOrder.push_back(rBond->getIdx());
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (nUnknown == 1) {
    // RDKit❗🔝:       if (reactantAtom.getDegree() == productAtom->getDegree()) {
    // RDKit❗🔝:         // there's a reactant bond that hasn't yet been accounted for:
    // RDKit❗🔝:         int unmatchedBond = -1;
    // RDKit❗🔝:
    // RDKit❗🔝:         for (const auto rBond : reactant.atomBonds(&reactantAtom)) {
    // RDKit❗🔝:           if (std::find(pOrder.begin(), pOrder.end(), rBond->getIdx()) ==
    // RDKit❗🔝:               pOrder.end()) {
    // RDKit❗🔝:             unmatchedBond = rBond->getIdx();
    // RDKit❗🔝:             break;
    // RDKit❗🔝:           }
    // RDKit❗🔝:         }
    // RDKit❗🔝:         // what must be true at this point:
    // RDKit❗🔝:         //  1) there's a -1 in pOrder that we'll substitute for
    // RDKit❗🔝:         //  2) unmatchedBond contains the index of the substitution
    // RDKit❗🔝:         auto bPos = std::find(pOrder.begin(), pOrder.end(), -1);
    // RDKit❗🔝:         if (unmatchedBond >= 0 && bPos != pOrder.end()) {
    // RDKit❗🔝:           *bPos = unmatchedBond;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         nUnknown = 0;
    // RDKit❗🔝:         CHECK_INVARIANT(
    // RDKit❗🔝:             std::find(pOrder.begin(), pOrder.end(), -1) == pOrder.end(),
    // RDKit❗🔝:             "extra unmapped atom");
    // RDKit❗🔝:       } else if (productAtom->getDegree() > reactantAtom.getDegree()) {
    // RDKit❗🔝:         // the product has an extra bond. we can just remove the -1 from the
    // RDKit❗🔝:         // list:
    // RDKit❗🔝:         auto bPos = std::find(pOrder.begin(), pOrder.end(), -1);
    // RDKit❗🔝:         pOrder.erase(bPos);
    // RDKit❗🔝:         nUnknown = 0;
    // RDKit❗🔝:         CHECK_INVARIANT(
    // RDKit❗🔝:             std::find(pOrder.begin(), pOrder.end(), -1) == pOrder.end(),
    // RDKit❗🔝:             "extra unmapped atom");
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (!nUnknown) {
    // RDKit❗🔝:       if (reactantAtom.getDegree() > productAtom->getDegree()) {
    // RDKit❗🔝:         // we lost a bond from the reactant.
    // RDKit❗🔝:         // we can just remove the unmatched reactant bond from the list
    // RDKit❗🔝:         INT_LIST::iterator rOrderIter = rOrder.begin();
    // RDKit❗🔝:         while (rOrderIter != rOrder.end() && rOrder.size() > pOrder.size()) {
    // RDKit❗🔝:           // we may invalidate the iterator so keep track of what comes next:
    // RDKit❗🔝:           auto thisOne = rOrderIter++;
    // RDKit❗🔝:           if (std::find(pOrder.begin(), pOrder.end(), *thisOne) ==
    // RDKit❗🔝:               pOrder.end()) {
    // RDKit❗🔝:             // not in the products:
    // RDKit❗🔝:             rOrder.erase(thisOne);
    // RDKit❗🔝:           }
    // RDKit❗🔝:         }
    // RDKit❗🔝:       }
    // RDKit❗🔝:       productAtom->setChiralTag(reactantAtom.getChiralTag());
    // RDKit❗🔝:       int nSwaps = countSwapsToInterconvert(rOrder, pOrder);
    // RDKit❗🔝:       bool invert = false;
    // RDKit❗🔝:       if (nSwaps % 2) {
    // RDKit❗🔝:         invert = true;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       int inversionFlag;
    // RDKit❗🔝:       if (productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗🔝:                                         inversionFlag) &&
    // RDKit❗🔝:           inversionFlag == 1) {
    // RDKit❗🔝:         invert = !invert;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       if (invert) {
    // RDKit❗🔝:         productAtom->invertChirality();
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Cost improvement: preserve source physical list order and first-match
    // semantics in contiguous Vec lists instead of allocating one std::list
    // node per bond. Same degree-local O(d^2) matching; at most one unknown
    // substitution/removal, no topology copies or independent stereo logic.
    crate::materialize::source_u32("reactant atom", reactant_atom)?;
    let matches = mapping
        .reactant_to_product
        .entry(reactant_atom)
        .or_default();
    let mut source_i = 0u32;
    while (source_i as usize) < matches.len() {
        let row = matches[source_i as usize];
        source_i = source_i.wrapping_add(1);
        let p = product.topology.atoms.get(row).ok_or_else(|| {
            invariant(
                "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                "product atom row missing",
                Some(reactant_atom),
                Some(row),
                None,
            )
        })?;
        let r = input.topology.atoms.get(reactant_atom).ok_or_else(|| {
            invariant(
                "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                "reactant atom row missing",
                Some(reactant_atom),
                Some(row),
                None,
            )
        })?;
        let flag = inversion_flag(p)?.unwrap_or(0);
        if matches!(r.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) || flag > 2 {
            continue;
        }
        let reactant_index = r.id().index();
        let product_index = p.id().index();
        let reactant_neighbors = input
            .topology
            .adjacency
            .try_neighbors_of(reactant_index)
            .ok_or_else(|| {
                invariant(
                    "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                    "reactant adjacency row missing",
                    Some(reactant_index),
                    Some(row),
                    None,
                )
            })?;
        if reactant_neighbors.len() < 3 {
            continue;
        }
        let pn = product.neighbors.get(product_index).ok_or_else(|| {
            invariant(
                "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                "product adjacency row missing",
                Some(reactant_index),
                Some(product_index),
                None,
            )
        })?;
        if pn.len() < 3 || reactant_neighbors.len().abs_diff(pn.len()) > 1 {
            continue;
        }
        let mut unknown = 0u32;
        let mut r_order = Vec::with_capacity(reactant_neighbors.len());
        for neighbor in reactant_neighbors {
            let bond = input
                .topology
                .bonds
                .get(neighbor.bond.index())
                .ok_or_else(|| {
                    invariant(
                        "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                        "reactant bond row missing",
                        Some(reactant_index),
                        Some(row),
                        Some(neighbor.bond),
                    )
                })?;
            r_order
                .push(crate::materialize::source_u32("reactant bond", bond.id().index())? as i32);
        }
        let mut p_order = Vec::with_capacity(pn.len());
        for neighbor in pn {
            let bond = if let Some(&other) = mapping.product_to_reactant.get(&neighbor.atom_index) {
                cosmolkit_model::source_bond_between_atoms(
                    input.topology.atoms.len(),
                    AtomId::new(reactant_index),
                    AtomId::new(other),
                    || reactant_neighbors.iter(),
                )?
            } else {
                None
            };
            if let Some(bond) = bond {
                // Native repeats the immutable source edge lookup in its else
                // branch. Reuse that exact found row without a second scan.
                let bond = input.topology.bonds.get(bond.index()).ok_or_else(|| {
                    invariant(
                        "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                        "expected reactant bond not found",
                        Some(reactant_index),
                        Some(row),
                        Some(bond),
                    )
                })?;
                p_order.push(
                    crate::materialize::source_u32("reactant bond", bond.id().index())? as i32,
                );
            } else {
                unknown = unknown.wrapping_add(1);
                if unknown > 1 {
                    break;
                }
                p_order.push(-1);
            }
        }
        if unknown == 1 {
            if reactant_neighbors.len() == pn.len() {
                let unmatched = r_order.iter().copied().find(|b| !p_order.contains(b));
                if let (Some(bond), Some(position)) =
                    (unmatched, p_order.iter().position(|b| *b == -1))
                {
                    if bond >= 0 {
                        p_order[position] = bond;
                    }
                }
                unknown = 0;
                if p_order.contains(&-1) {
                    return Err(invariant(
                        "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                        "extra unmapped atom",
                        Some(reactant_atom),
                        Some(row),
                        None,
                    ));
                }
            } else if pn.len() > reactant_neighbors.len() {
                let position = p_order.iter().position(|b| *b == -1).ok_or_else(|| {
                    invariant(
                        "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                        "extra unmapped atom",
                        Some(reactant_atom),
                        Some(row),
                        None,
                    )
                })?;
                p_order.remove(position);
                unknown = 0;
                if p_order.contains(&-1) {
                    return Err(invariant(
                        "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                        "extra unmapped atom",
                        Some(reactant_atom),
                        Some(row),
                        None,
                    ));
                }
            }
        }
        if unknown == 0 {
            if reactant_neighbors.len() > pn.len() {
                let mut position = 0;
                while position < r_order.len() && r_order.len() > p_order.len() {
                    if !p_order.contains(&r_order[position]) {
                        r_order.remove(position);
                    } else {
                        position += 1;
                    }
                }
            }
            product.topology.atoms[row].set_chiral_tag(r.chiral_tag());
            let mut invert = count_swaps_to_interconvert(&r_order, &p_order)? % 2 != 0;
            if inversion_flag(&product.topology.atoms[row])? == Some(1) {
                invert = !invert;
            }
            if invert {
                invert_atom_chirality(&mut product.topology.atoms[row])?;
            }
        }
    }
    Ok(())
}

pub(crate) fn correct_unmatched_chirality(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &mut ReactantProductMapping,
    atoms: &[usize],
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: checkAndCorrectChiralityOfProduct
    // RDKit❗🔝: void checkAndCorrectChiralityOfProduct(
    // RDKit❗🔝:     const std::vector<const Atom *> &chiralAtomsToCheck, RWMOL_SPTR product,
    // RDKit❗🔝:     ReactantProductAtomMapping *mapping) {
    // RDKit❗🔝:   for (auto reactantAtom : chiralAtomsToCheck) {
    // RDKit❗🔝:     CHECK_INVARIANT(reactantAtom->getChiralTag() != Atom::CHI_UNSPECIFIED,
    // RDKit❗🔝:                     "missing atom chirality.");
    // RDKit❗🔝:     const auto reactAtomDegree =
    // RDKit❗🔝:         reactantAtom->getOwningMol().getAtomDegree(reactantAtom);
    // RDKit❗🔝:     for (unsigned i = 0;
    // RDKit❗🔝:          i < mapping->reactProdAtomMap[reactantAtom->getIdx()].size(); i++) {
    // RDKit❗🔝:       unsigned productAtomIdx =
    // RDKit❗🔝:           mapping->reactProdAtomMap[reactantAtom->getIdx()][i];
    // RDKit❗🔝:       Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗🔝:       CHECK_INVARIANT(
    // RDKit❗🔝:           reactantAtom->getChiralTag() == productAtom->getChiralTag(),
    // RDKit❗🔝:           "invalid product chirality.");
    // RDKit❗🔝:
    // RDKit❗🔝:       if (reactAtomDegree != product->getAtomDegree(productAtom)) {
    // RDKit❗🔝:         // If the number of bonds to the atom has changed in the course of the
    // RDKit❗🔝:         // reaction we're lost, so remove chirality.
    // RDKit❗🔝:         //  A word of explanation here: the atoms in the chiralAtomsToCheck
    // RDKit❗🔝:         //  set are not explicitly mapped atoms of the reaction, so we really
    // RDKit❗🔝:         //  have no idea what to do with this case. At the moment I'm not even
    // RDKit❗🔝:         //  really sure how this could happen, but better safe than sorry.
    // RDKit❗🔝:         productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗🔝:       } else if (reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗🔝:                  reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗🔝:         // this will contain the indices of product bonds in the
    // RDKit❗🔝:         // reactant order:
    // RDKit❗🔝:         INT_LIST newOrder;
    // RDKit❗🔝:         ROMol::OEDGE_ITER beg, end;
    // RDKit❗🔝:         boost::tie(beg, end) =
    // RDKit❗🔝:             reactantAtom->getOwningMol().getAtomBonds(reactantAtom);
    // RDKit❗🔝:         while (beg != end) {
    // RDKit❗🔝:           const Bond *reactantBond = reactantAtom->getOwningMol()[*beg];
    // RDKit❗🔝:           unsigned int oAtomIdx =
    // RDKit❗🔝:               reactantBond->getOtherAtomIdx(reactantAtom->getIdx());
    // RDKit❗🔝:           CHECK_INVARIANT(mapping->reactProdAtomMap.find(oAtomIdx) !=
    // RDKit❗🔝:                               mapping->reactProdAtomMap.end(),
    // RDKit❗🔝:                           "other atom from bond not mapped.");
    // RDKit❗🔝:           const Bond *productBond;
    // RDKit❗🔝:           unsigned neighborBondIdx = mapping->reactProdAtomMap[oAtomIdx][i];
    // RDKit❗🔝:           productBond = product->getBondBetweenAtoms(productAtom->getIdx(),
    // RDKit❗🔝:                                                      neighborBondIdx);
    // RDKit❗🔝:           CHECK_INVARIANT(productBond, "no matching bond found in product");
    // RDKit❗🔝:           newOrder.push_back(productBond->getIdx());
    // RDKit❗🔝:           ++beg;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         int nSwaps = productAtom->getPerturbationOrder(newOrder);
    // RDKit❗🔝:         if (nSwaps % 2) {
    // RDKit❗🔝:           productAtom->invertChirality();
    // RDKit❗🔝:         }
    // RDKit❗🔝:       } else {
    // RDKit❗🔝:         // not tetrahedral chirality, don't do anything.
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }  // end of loop over chiralAtomsToCheck
    // RDKit❗🔝: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Cost improvement: source INT_LIST newOrder/reference/probe use linked
    // nodes. Contiguous degree-local vectors retain order with fewer allocations
    // and pointer traversals. Canonical CORE perturbation/inversion owns parity;
    // no graph copies, inferred chirality or extra degree restrictions.
    for &reactant_atom in atoms {
        let r = input.topology.atoms.get(reactant_atom).ok_or_else(|| {
            invariant(
                "checkAndCorrectChiralityOfProduct",
                "reactant atom row missing",
                Some(reactant_atom),
                None,
                None,
            )
        })?;
        let reactant_index = r.id().index();
        if r.chiral_tag() == ChiralTag::Unspecified {
            return Err(invariant(
                "checkAndCorrectChiralityOfProduct",
                "missing atom chirality.",
                Some(reactant_index),
                None,
                None,
            ));
        }
        let rn = input
            .topology
            .adjacency
            .try_neighbors_of(reactant_index)
            .ok_or_else(|| {
                invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "reactant adjacency row missing",
                    Some(reactant_index),
                    None,
                    None,
                )
            })?;
        mapping
            .reactant_to_product
            .entry(reactant_index)
            .or_default();
        let rows = &mapping.reactant_to_product[&reactant_index];
        let mut duplicate = 0u32;
        while (duplicate as usize) < rows.len() {
            let copy_index = duplicate as usize;
            let row = rows[copy_index];
            duplicate = duplicate.wrapping_add(1);
            let p = product.topology.atoms.get(row).ok_or_else(|| {
                invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "product atom row missing",
                    Some(reactant_index),
                    Some(row),
                    None,
                )
            })?;
            if r.chiral_tag() != p.chiral_tag() {
                return Err(invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "invalid product chirality.",
                    Some(reactant_index),
                    Some(row),
                    None,
                ));
            }
            let product_index = p.id().index();
            let pn = product.neighbors.get(product_index).ok_or_else(|| {
                invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "product adjacency row missing",
                    Some(reactant_index),
                    Some(product_index),
                    None,
                )
            })?;
            if rn.len() != pn.len() {
                product.topology.atoms[row].set_chiral_tag(ChiralTag::Unspecified);
            } else if matches!(
                r.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                let mut new_order = Vec::with_capacity(rn.len());
                for neighbor in rn {
                    let reactant_bond = input
                        .topology
                        .bonds
                        .get(neighbor.bond.index())
                        .ok_or_else(|| {
                            invariant(
                                "checkAndCorrectChiralityOfProduct",
                                "reactant bond row missing",
                                Some(reactant_index),
                                Some(row),
                                Some(neighbor.bond),
                            )
                        })?;
                    let other = other_atom_index(reactant_bond, r.id())?;
                    // The current center rows stay borrowed. Native [] for an
                    // already-found other key cannot insert or modify them.
                    let other_rows = if other.index() == reactant_index {
                        rows.as_slice()
                    } else {
                        mapping
                            .reactant_to_product
                            .get(&other.index())
                            .ok_or_else(|| {
                                invariant(
                                    "checkAndCorrectChiralityOfProduct",
                                    "other atom from bond not mapped.",
                                    Some(other.index()),
                                    Some(row),
                                    Some(neighbor.bond),
                                )
                            })?
                            .as_slice()
                    };
                    let other_product = *other_rows.get(copy_index).ok_or_else(|| {
                        invariant(
                            "checkAndCorrectChiralityOfProduct",
                            "other atom product copy missing",
                            Some(other.index()),
                            Some(row),
                            Some(neighbor.bond),
                        )
                    })?;
                    let product_bond = cosmolkit_model::source_bond_between_atoms(
                        product.topology.atoms.len(),
                        AtomId::new(product_index),
                        AtomId::new(other_product),
                        || pn.iter(),
                    )?
                    .ok_or_else(|| {
                        invariant(
                            "checkAndCorrectChiralityOfProduct",
                            "no matching bond found in product",
                            Some(other.index()),
                            Some(row),
                            None,
                        )
                    })?;
                    let product_bond = product
                        .topology
                        .bonds
                        .get(product_bond.index())
                        .ok_or_else(|| {
                            invariant(
                                "checkAndCorrectChiralityOfProduct",
                                "product bond row missing",
                                Some(other.index()),
                                Some(row),
                                Some(product_bond),
                            )
                        })?;
                    new_order.push(crate::materialize::source_u32(
                        "product bond",
                        product_bond.id().index(),
                    )? as i32);
                }
                // Read actual source Bond::getIdx fields in product incidence
                // order, then let the sole CORE getPerturbationOrder owner do
                // source signed-width conversion and probe/reference matching.
                let current_order = pn
                    .iter()
                    .map(|neighbor| {
                        product
                            .topology
                            .bonds
                            .get(neighbor.bond.index())
                            .map(|b| b.id().index())
                            .ok_or_else(|| {
                                invariant(
                                    "Atom::getPerturbationOrder",
                                    "product bond row missing",
                                    None,
                                    Some(product_index),
                                    Some(neighbor.bond),
                                )
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                if cosmolkit_core::atom_perturbation_order(&new_order, current_order)? % 2 != 0 {
                    invert_atom_chirality(&mut product.topology.atoms[row])?;
                }
            }
        }
    }
    Ok(())
}

fn new_group(source: &StereoGroup, atoms: Vec<AtomId>) -> StereoGroup {
    // RDKit❗🔝: StereoGroup::StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
    // RDKit❗🔝:                          std::vector<Bond *> &&bonds, unsigned readId)
    // RDKit❗🔝:     : d_grouptype(grouptype),
    // RDKit❗🔝:       d_atoms(atoms),
    // RDKit❗🔝:       d_bonds(bonds),
    // RDKit❗🔝:       d_readId{readId} {}
    // Moving detached member vectors avoids the source rvalue constructor
    // copying named lvalue parameters into d_atoms/d_bonds.
    // The modeled optional source ID remains distinct; constructor write ID 0
    // and empty bond members reproduce the runner's explicit constructors.
    let group = StereoGroup::new(source.kind(), atoms, Vec::new());
    match source.id() {
        Some(id) => group.with_id(id),
        None => group,
    }
}

pub(crate) fn copy_enhanced_groups(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: copyEnhancedStereoGroups
    // RDKit❗🔝: void copyEnhancedStereoGroups(const ROMol &reactant, RWMOL_SPTR product,
    // RDKit❗🔝:                               const ReactantProductAtomMapping &mapping) {
    // RDKit❗🔝:   std::vector<StereoGroup> new_stereo_groups;
    // RDKit❗🔝:   for (const auto &sg : reactant.getStereoGroups()) {
    // RDKit❗🔝:     std::vector<Atom *> atoms;
    // RDKit❗🔝:     std::vector<Bond *> bonds;
    // RDKit❗🔝:     for (const auto &reactantAtom : sg.getAtoms()) {
    // RDKit❗🔝:       auto productAtoms = mapping.reactProdAtomMap.find(reactantAtom->getIdx());
    // RDKit❗🔝:       if (productAtoms == mapping.reactProdAtomMap.end()) {
    // RDKit❗🔝:         continue;
    // RDKit❗🔝:       }
    // RDKit❗🔝:
    // RDKit❗🔝:       for (auto &productAtomIdx : productAtoms->second) {
    // RDKit❗🔝:         auto productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗🔝:         // If chirality destroyed by the reaction, skip the atom
    // RDKit❗🔝:         if (productAtom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    // RDKit❗🔝:           continue;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         // If chirality defined explicitly by the reaction, skip the atom
    // RDKit❗🔝:         int flagVal = 0;
    // RDKit❗🔝:         productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗🔝:                                       flagVal);
    // RDKit❗🔝:         if (flagVal == 4) {
    // RDKit❗🔝:           continue;
    // RDKit❗🔝:         }
    // RDKit❗🔝:         atoms.push_back(productAtom);
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     if (!atoms.empty()) {
    // RDKit❗🔝:       new_stereo_groups.emplace_back(sg.getGroupType(), std::move(atoms),
    // RDKit❗🔝:                                      std::move(bonds), sg.getReadId());
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:
    // RDKit❗🔝:   // Although we have added storage, and canonicalization of Atropisomers,
    // RDKit❗🔝:   // searching is not yet supported.  When it is, we will need to copy
    // RDKit❗🔝:   // bond-part of the SG groups to the products as appropriate.
    // RDKit❗🔝:
    // RDKit❗🔝:   if (!new_stereo_groups.empty()) {
    // RDKit❗🔝:     auto &existing_sg = product->getStereoGroups();
    // RDKit❗🔝:     new_stereo_groups.insert(new_stereo_groups.end(), existing_sg.begin(),
    // RDKit❗🔝:                              existing_sg.end());
    // RDKit❗🔝:     product->setStereoGroups(std::move(new_stereo_groups));
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Cost improvement: new_group moves member vectors and the final append
    // moves existing groups rather than deep-copying their vectors as source.
    // Preserve physical group/member/copy order, duplicate pointers and staging
    // before the sole MODEL setStereoGroups merge, with no sorted membership.
    let mut groups = Vec::new();
    for group in &input.topology.stereo_groups {
        let mut atoms = Vec::new();
        for atom in group.atoms() {
            let source_atom = input.topology.atoms.get(atom.index()).ok_or_else(|| {
                invariant(
                    "copyEnhancedStereoGroups",
                    "reactant stereo-group atom row missing",
                    Some(atom.index()),
                    None,
                    None,
                )
            })?;
            let Some(rows) = mapping.reactant_to_product.get(&source_atom.id().index()) else {
                continue;
            };
            for &row in rows {
                let p = product.topology.atoms.get(row).ok_or_else(|| {
                    invariant(
                        "copyEnhancedStereoGroups",
                        "product stereo-group atom row missing",
                        Some(source_atom.id().index()),
                        Some(row),
                        None,
                    )
                })?;
                if p.chiral_tag() == ChiralTag::Unspecified {
                    continue;
                }
                if inversion_flag(p)? == Some(4) {
                    continue;
                }
                atoms.push(p.id());
            }
        }
        if !atoms.is_empty() {
            groups.push(new_group(group, atoms));
        }
    }
    if !groups.is_empty() {
        // Moving existing rows avoids the source's SG deep copies and retains
        // exact source prepend order, IDs and bond membership.
        groups.append(&mut product.topology.stereo_groups);
        product.topology.stereo_groups = cosmolkit_model::merge_absolute_stereo_groups(groups);
    }
    Ok(())
}

pub(crate) fn copy_template_groups(
    product: &mut ProductBuilder,
    template: &QueryGraph,
    template_index: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: copyTemplateStereoGroupsToMol
    // RDKit❗❌: void copyTemplateStereoGroupsToMol(const ROMol &templateMol,
    // RDKit❗❌:                                    RWMOL_SPTR product) {
    // RDKit❗❌:   const auto &stereoGroups = templateMol.getStereoGroups();
    // RDKit❗❌:   if (stereoGroups.empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> atomsInTemplateStereoGroups(product->getNumAtoms());
    // RDKit❗❌:   std::vector<StereoGroup> newStereoGroups;
    // RDKit❗❌:   for (const auto &sg : stereoGroups) {
    // RDKit❗❌:     bool keepIt = true;
    // RDKit❗❌:     std::vector<Atom *> atoms;
    // RDKit❗❌:     for (const auto &atom : sg.getAtoms()) {
    // RDKit❗❌:       if (auto mapNum = atom->getAtomMapNum()) {
    // RDKit❗❌:         for (auto productAtom : product->atoms()) {
    // RDKit❗❌:           int oldMapNum = 0;
    // RDKit❗❌:           if (productAtom->getPropIfPresent(common_properties::reactionMapNum,
    // RDKit❗❌:                                             oldMapNum) &&
    // RDKit❗❌:               oldMapNum == mapNum) {
    // RDKit❗❌:             atoms.push_back(productAtom);
    // RDKit❗❌:             atomsInTemplateStereoGroups.set(productAtom->getIdx());
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       } else {
    // RDKit❗❌:         keepIt = false;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (keepIt && !atoms.empty()) {
    // RDKit❗❌:       std::vector<Bond *> bonds;
    // RDKit❗❌:       newStereoGroups.emplace_back(sg.getGroupType(), std::move(atoms),
    // RDKit❗❌:                                    std::move(bonds), sg.getReadId());
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   if (!newStereoGroups.empty()) {
    // RDKit❗❌:     // remove any stereo groups that are already present in the product (these
    // RDKit❗❌:     // were copied over from the reactant in copyEnhancedStereoGroups()) and
    // RDKit❗❌:     // that overlap with the added ones
    // RDKit❗❌:     for (const auto &productSG : product->getStereoGroups()) {
    // RDKit❗❌:       unsigned int nOverlappingAtoms = 0;
    // RDKit❗❌:       for (const auto atom : productSG.getAtoms()) {
    // RDKit❗❌:         if (atomsInTemplateStereoGroups[atom->getIdx()]) {
    // RDKit❗❌:           ++nOverlappingAtoms;
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       if (!nOverlappingAtoms) {
    // RDKit❗❌:         // no overlapping atoms, we can just keep the stereogroup.
    // RDKit❗❌:         newStereoGroups.push_back(productSG);
    // RDKit❗❌:       } else if (nOverlappingAtoms < productSG.getAtoms().size()) {
    // RDKit❗❌:         // some of the atoms in the stereo group are not already there
    // RDKit❗❌:         // in the product, we need to split the stereo group
    // RDKit❗❌:         std::vector<Atom *> newAtoms;
    // RDKit❗❌:         for (const auto atom : productSG.getAtoms()) {
    // RDKit❗❌:           if (!atomsInTemplateStereoGroups[atom->getIdx()]) {
    // RDKit❗❌:             newAtoms.push_back(atom);
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:         std::vector<Bond *> newBonds;
    // RDKit❗❌:         newStereoGroups.emplace_back(productSG.getGroupType(),
    // RDKit❗❌:                                      std::move(newAtoms), std::move(newBonds),
    // RDKit❗❌:                                      productSG.getReadId());
    // RDKit❗❌:       }
    // RDKit❗❌:       // else: all atoms in the stereo group are already there, we can skip it
    // RDKit❗❌:     }
    // RDKit❗❌:     product->setStereoGroups(std::move(newStereoGroups));
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    if template.stereo_groups().is_empty() {
        return Ok(());
    }
    // Source physical scans and copied existing-group staging are preserved.
    // Known cost gap: Vec<bool> stores byte flags whereas Boost dynamic_bitset
    // packs bits; keep ❌ rather than hide that extra temporary mask footprint.
    let mut marked = vec![false; product.topology.atoms.len()];
    let mut groups = Vec::new();
    for group in template.stereo_groups() {
        let mut keep = true;
        let mut atoms = Vec::new();
        for atom in group.atoms() {
            let t = template.atom(atom.index()).ok_or_else(|| {
                invariant(
                    "copyTemplateStereoGroupsToMol",
                    "template stereo atom missing",
                    None,
                    None,
                    None,
                )
            })?;
            let map =
                crate::validation::atom_map(t, ReactionRole::Product, template_index)?.unwrap_or(0);
            if map != 0 {
                for p in &product.topology.atoms {
                    if int_prop(p, "old_mapno")? == Some(map) {
                        atoms.push(p.id());
                        *marked.get_mut(p.id().index()).ok_or_else(|| {
                            invariant(
                                "copyTemplateStereoGroupsToMol",
                                "stereo-group atom bit index out of range",
                                None,
                                Some(p.id().index()),
                                None,
                            )
                        })? = true;
                    }
                }
            } else {
                keep = false;
                break;
            }
        }
        // Marks from a dropped group deliberately persist, as in source.
        if keep && !atoms.is_empty() {
            groups.push(new_group(group, atoms));
        }
    }
    if !groups.is_empty() {
        // Source iterates existing groups by const reference and replaces the
        // complete group set only at the final setStereoGroups statement. Keep
        // this staging even when a reached detached row/bit read fails.
        let is_marked = |atom: AtomId| {
            let source_atom = product.topology.atoms.get(atom.index()).ok_or_else(|| {
                invariant(
                    "copyTemplateStereoGroupsToMol",
                    "existing stereo-group atom row missing",
                    None,
                    Some(atom.index()),
                    None,
                )
            })?;
            marked
                .get(source_atom.id().index())
                .copied()
                .ok_or_else(|| {
                    invariant(
                        "copyTemplateStereoGroupsToMol",
                        "stereo-group atom bit index out of range",
                        None,
                        Some(source_atom.id().index()),
                        None,
                    )
                })
        };
        for group in &product.topology.stereo_groups {
            let mut overlaps = 0u32;
            for &atom in group.atoms() {
                if is_marked(atom)? {
                    overlaps = overlaps.wrapping_add(1);
                }
            }
            if overlaps == 0 {
                groups.push(group.clone());
            } else if (overlaps as usize) < group.atoms().len() {
                let mut atoms = Vec::new();
                for &atom in group.atoms() {
                    if !is_marked(atom)? {
                        atoms.push(atom);
                    }
                }
                groups.push(new_group(group, atoms));
            }
        }
        product.topology.stereo_groups = cosmolkit_model::merge_absolute_stereo_groups(groups);
    }
    Ok(())
}

pub(crate) fn propagate_coordinates(
    points: &mut Vec<[f64; 3]>,
    is_3d: &mut bool,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
    selection: ReactionCoordinateSelection,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: generateProductConformers
    // RDKit❗🔝: void generateProductConformers(Conformer *productConf, const ROMol &reactant,
    // RDKit❗🔝:                                ReactantProductAtomMapping *mapping) {
    // RDKit❗🔝:   if (!reactant.getNumConformers()) {
    // RDKit❗🔝:     return;
    // RDKit❗🔝:   }
    // RDKit❗🔝:   const Conformer &reactConf = reactant.getConformer();
    // RDKit❗🔝:   if (reactConf.is3D()) {
    // RDKit❗🔝:     productConf->set3D(true);
    // RDKit❗🔝:   }
    // RDKit❗🔝:   for (std::map<unsigned int, std::vector<unsigned int>>::const_iterator pr =
    // RDKit❗🔝:            mapping->reactProdAtomMap.begin();
    // RDKit❗🔝:        pr != mapping->reactProdAtomMap.end(); ++pr) {
    // RDKit❗🔝:     std::vector<unsigned> prodIdxs = pr->second;
    // RDKit❗🔝:     if (prodIdxs.size() > 1) {
    // RDKit❗🔝:       BOOST_LOG(rdWarningLog) << "reactant atom match more than one product "
    // RDKit❗🔝:                                  "atom, coordinates need to be revised\n";
    // RDKit❗🔝:     }
    // RDKit❗🔝:     // is this reliable when multiple product atom mapping occurs????
    // RDKit❗🔝:     for (unsigned int prodIdx : prodIdxs) {
    // RDKit❗🔝:       productConf->setAtomPos(prodIdx, reactConf.getAtomPos(pr->first));
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗🔝: const RDGeom::Point3D &Conformer::getAtomPos(unsigned int atomId) const {
    // RDKit❗🔝:   if (dp_mol) {
    // RDKit❗🔝:     PRECONDITION(dp_mol->getNumAtoms() == d_positions.size(), "");
    // RDKit❗🔝:   }
    // RDKit❗🔝:   URANGE_CHECK(atomId, d_positions.size());
    // RDKit❗🔝:   return d_positions.at(atomId);
    // RDKit❗🔝: }
    // The sole caller performs Native productConf->resize before this function.
    // Source absent conformers return before selection/flags/mapping reads.
    if input.coordinates.conformers_2d.is_empty() && input.coordinates.conformers_3d.is_empty() {
        return Ok(());
    }
    let Some(source) =
        cosmolkit_smiles::select_cx_coordinates(input.coordinates, selection.into())?
    else {
        return Ok(());
    };
    if let cosmolkit_smiles::CoordinateSource::ThreeD(conf) = source {
        if conf.is_3d() {
            *is_3d = true;
        }
    }
    // Cost improvement: Native copies every mapped product-index vector. The
    // immutable private mapping is borrowed for the identical ordered writes,
    // avoiding those allocations/copies while preserving duplicates and order.
    for (&reactant_atom, rows) in &mapping.reactant_to_product {
        if rows.len() > 1 {
            eprintln!(
                "reactant atom match more than one product atom, coordinates need to be revised"
            );
        }
        for &row in rows {
            // Argument conversion precedes the Native getter body; its owning
            // conformer shape precondition is checked only on reached reads.
            let source_index =
                crate::materialize::source_u32("reactant coordinate atom", reactant_atom)? as usize;
            let point = match source {
                cosmolkit_smiles::CoordinateSource::ThreeD(conf) => {
                    conf.validate_for_atom_count(input.topology.atoms.len())?;
                    conf.coordinates().get(source_index).copied()
                }
                cosmolkit_smiles::CoordinateSource::TwoD(conf) => {
                    conf.validate_for_atom_count(input.topology.atoms.len())?;
                    conf.coordinates()
                        .get(source_index)
                        .map(|p| [p[0], p[1], 0.0])
                }
            }
            .ok_or_else(|| {
                invariant(
                    "generateProductConformers",
                    "reactant conformer atom out of range",
                    Some(source_index),
                    None,
                    None,
                )
            })?;
            // Reuse the sole MODEL source Conformer::setAtomPos implementation:
            // max-unsigned throws, gaps zero-fill, stored IEEE bits are untouched.
            cosmolkit_model::source_set_atom_position(points, row, point)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod complete_product_anchor_candidates_source_tests {
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};
    fn mapping(rows: &[(usize, Vec<usize>)]) -> ReactantProductMapping {
        ReactantProductMapping {
            mapped: vec![],
            skipped: vec![],
            reactant_to_product: rows.iter().cloned().collect(),
            product_to_reactant: BTreeMap::new(),
            product_atom_bond: BTreeMap::new(),
            template_bonds: BTreeSet::new(),
        }
    }
    #[test]
    fn primary_anchor_takes_precedence_and_preserves_duplicate_encounter_order() {
        for non_anchor in [None, Some(0)] {
            let cap = StereoBondEndCap {
                anchor: 7,
                non_anchor,
            };
            let mapping = mapping(&[(7, vec![3, 0, 3]), (0, vec![9])]);
            let (rows, swap) = cap.product_candidates(&mapping).unwrap();
            assert_eq!(rows, [3, 0, 3]);
            assert!(!swap);
        }
    }
    #[test]
    fn present_empty_primary_key_still_wins_over_present_nonempty_secondary_key() {
        let cap = StereoBondEndCap {
            anchor: 7,
            non_anchor: Some(0),
        };
        let mapping = mapping(&[(7, vec![]), (0, vec![9])]);
        let (rows, swap) = cap.product_candidates(&mapping).unwrap();
        assert!(rows.is_empty());
        assert!(!swap);
    }
    #[test]
    fn secondary_index_zero_is_present_and_sets_swap_even_for_empty_matches() {
        for matches in [vec![], vec![4, 1, 4]] {
            let cap = StereoBondEndCap {
                anchor: 7,
                non_anchor: Some(0),
            };
            let mapping = mapping(&[(0, matches.clone())]);
            let (rows, swap) = cap.product_candidates(&mapping).unwrap();
            assert_eq!(rows, matches);
            assert!(swap);
            assert!(cap.has_non_anchor());
            assert_eq!(cap.non_anchor_idx().unwrap(), 0);
        }
    }
    #[test]
    fn absent_secondary_pointer_does_not_use_zero_key_or_dereference_it() {
        let cap = StereoBondEndCap {
            anchor: 7,
            non_anchor: None,
        };
        let mapping = mapping(&[(0, vec![9])]);
        let (rows, swap) = cap.product_candidates(&mapping).unwrap();
        assert!(rows.is_empty());
        assert!(!swap);
        assert!(!cap.has_non_anchor());
    }
    #[test]
    fn absent_both_keys_and_identical_anchor_indices_preserve_source_false_flag() {
        for rows in [vec![], vec![(0, vec![2, 1])]] {
            let cap = StereoBondEndCap {
                anchor: 0,
                non_anchor: Some(0),
            };
            let mapping = mapping(&rows);
            let (selected, swap) = cap.product_candidates(&mapping).unwrap();
            assert_eq!(
                selected,
                rows.first().map_or(&[][..], |(_, v)| v.as_slice())
            );
            assert!(!swap);
            assert_eq!(cap.anchor_idx(), 0);
        }
    }
    #[test]
    fn selection_is_read_only_and_borrows_the_exact_selected_map_vector() {
        for secondary in [false, true] {
            let cap = StereoBondEndCap {
                anchor: 7,
                non_anchor: Some(0),
            };
            let selected_key = if secondary { 0 } else { 7 };
            let mapping = mapping(&[(selected_key, vec![4, 2, 4])]);
            let before = mapping.reactant_to_product.clone();
            let (rows, swap) = cap.product_candidates(&mapping).unwrap();
            assert_eq!(
                rows.as_ptr(),
                mapping.reactant_to_product[&selected_key].as_ptr()
            );
            assert_eq!(swap, secondary);
            assert_eq!(mapping.reactant_to_product, before);
        }
    }
}

#[cfg(test)]
mod complete_anchor_index_source_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, TopologyBlock};
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;
    fn product(count: usize) -> ProductBuilder {
        ProductBuilder {
            topology: TopologyBlock {
                atoms: (0..count)
                    .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                    .collect(),
                ..TopologyBlock::default()
            },
            neighbors: vec![vec![]; count],
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; count],
            bond_origins: vec![],
        }
    }
    fn connect(product: &mut ProductBuilder, begin: usize, end: usize) {
        cosmolkit_model::add_source_bond_order(
            &mut product.topology.atoms,
            &mut product.topology.bonds,
            cosmolkit_model::SourceBondNeighbors {
                original: None,
                appended: &mut product.neighbors,
            },
            None,
            AtomId::new(begin),
            AtomId::new(end),
            BondOrder::Single,
        )
        .unwrap();
        product.bond_origins.push(None);
    }
    #[test]
    fn singleton_returns_candidate_without_owner_atom_or_connectivity_lookup() {
        assert_eq!(anchor_index(&product(0), 99, &[77]).unwrap(), 77);
    }
    #[test]
    fn empty_candidates_never_read_adjacency_but_still_report_source_match_invariant() {
        let mut product = product(1);
        product.neighbors.clear();
        assert!(matches!(
            anchor_index(&product, 0, &[]),
            Err(ReactionProductError::Invariant {
                stage: "reactProdMapAnchorIdx",
                detail: "match not found",
                ..
            })
        ));
    }
    #[test]
    fn first_bonded_candidate_in_given_order_wins_over_neighbor_order_and_duplicates() {
        let mut product = product(4);
        connect(&mut product, 0, 1);
        connect(&mut product, 0, 2);
        assert_eq!(anchor_index(&product, 0, &[3, 2, 1, 2]).unwrap(), 2);
    }
    #[test]
    fn reached_out_of_range_candidate_errors_before_later_match_but_unreached_bad_candidate_is_ignored()
     {
        let mut product = product(3);
        connect(&mut product, 0, 1);
        assert!(
            matches!(anchor_index(&product,0,&[9,1]),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(9))
        );
        assert_eq!(anchor_index(&product, 0, &[1, 9]).unwrap(), 1);
    }
    #[test]
    fn actual_atom_index_is_used_and_its_source_range_check_precedes_candidate_range_check() {
        let mut product = product(3);
        connect(&mut product, 2, 1);
        product.topology.atoms[0] = product.topology.atoms[0].clone().with_id(AtomId::new(2));
        assert_eq!(anchor_index(&product, 0, &[1, 0]).unwrap(), 1);
        product.topology.atoms[0] = product.topology.atoms[0].clone().with_id(AtomId::new(99));
        assert!(
            matches!(anchor_index(&product,0,&[88,1]),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
    }
    #[test]
    fn no_connected_candidate_reports_match_failure_without_zero_or_first_candidate_fallback() {
        assert!(matches!(
            anchor_index(&product(3), 0, &[1, 2]),
            Err(ReactionProductError::Invariant {
                detail: "match not found",
                ..
            })
        ));
    }
}

#[cfg(test)]
mod complete_stereo_end_cap_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondSpec, CoordinateBlock, MoleculeProperties,
        TopologyBlock,
    };
    use cosmolkit_types::Element;
    fn topology(count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let bonds: Vec<_> = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock {
            atoms: (0..count)
                .map(|i| {
                    Atom::from_spec(
                        AtomId::new(i),
                        AtomSpec::new(Element::C).with_no_implicit(true),
                    )
                })
                .collect(),
            adjacency: AdjacencyList::from_topology(count, &bonds),
            bonds,
            ..TopologyBlock::default()
        }
    }
    fn input<'a>(
        topology: &'a TopologyBlock,
        coordinates: &'a CoordinateBlock,
        properties: &'a MoleculeProperties,
    ) -> ReactionInput<'a> {
        ReactionInput {
            topology,
            coordinates,
            properties,
            rings: None,
            valence: None,
        }
    }
    #[test]
    fn first_non_anchor_neighbor_wins_in_physical_order_and_zero_is_present() {
        let t = topology(4, &[(1, 3), (1, 0), (1, 2)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let cap = StereoBondEndCap::new(&input(&t, &c, &p), 1, 3, 2).unwrap();
        assert_eq!(cap.anchor_idx(), 2);
        assert!(cap.has_non_anchor());
        assert_eq!(cap.non_anchor_idx().unwrap(), 0);
        let cap = StereoBondEndCap::new(&input(&t, &c, &p), 1, 3, 99).unwrap();
        assert_eq!(cap.non_anchor_idx().unwrap(), 0);
        assert_eq!(cap.anchor_idx(), 99);
    }
    #[test]
    fn only_other_and_anchor_or_no_neighbors_leaves_non_anchor_absent() {
        for edges in [vec![], vec![(0, 1)], vec![(0, 1), (0, 2)]] {
            let t = topology(3, &edges);
            let c = CoordinateBlock::default();
            let p = MoleculeProperties::default();
            let cap = StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 2).unwrap();
            assert!(!cap.has_non_anchor());
        }
    }
    #[test]
    fn explicit_and_implicit_hydrogens_participate_in_source_total_degree_precondition() {
        let mut t = topology(3, &[(0, 1), (0, 2)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        t.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(1),
        );
        assert!(StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 2).is_ok());
        t.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(2),
        );
        assert!(matches!(
            StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 2),
            Err(ReactionProductError::Invariant {
                detail: "Stereo Bond extremes must have less than four neighbors",
                ..
            })
        ));
        t.atoms[0] = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        let v = cosmolkit_core::ValenceAssignment {
            explicit_valence: vec![2; 3],
            implicit_hydrogens: vec![2, 0, 0],
        };
        let mut i = input(&t, &c, &p);
        i.valence = Some(&v);
        assert!(matches!(
            StereoBondEndCap::new(&i, 0, 1, 2),
            Err(ReactionProductError::Invariant {
                detail: "Stereo Bond extremes must have less than four neighbors",
                ..
            })
        ));
    }
    #[test]
    fn missing_implicit_valence_is_propagated_and_no_implicit_bypasses_it() {
        let mut t = topology(2, &[(0, 1)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        assert!(StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 8).is_ok());
        t.atoms[0] = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        assert!(
            matches!(StereoBondEndCap::new(&input(&t,&c,&p),0,1,8),Err(ReactionProductError::StereoGetter {atom,..}) if atom==AtomId::new(0))
        );
    }
    #[test]
    fn actual_atom_other_and_selected_atom_ids_are_used_instead_of_projected_rows() {
        let mut t = topology(4, &[(2, 1), (2, 3)]);
        t.atoms[0] = t.atoms[0].clone().with_id(AtomId::new(2));
        t.atoms[1] = t.atoms[1].clone().with_id(AtomId::new(3));
        t.atoms[3] = t.atoms[3].clone().with_id(AtomId::new(27));
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let cap = StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 1).unwrap();
        assert!(!cap.has_non_anchor());
        let cap = StereoBondEndCap::new(&input(&t, &c, &p), 0, 2, 1).unwrap();
        assert_eq!(cap.non_anchor_idx().unwrap(), 27);
    }
    #[test]
    fn missing_rows_and_missing_adjacency_are_structural_errors_without_empty_fallback() {
        let mut t = topology(2, &[]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        for (atom, other, detail) in [(8, 9, "no atom"), (0, 9, "no other double-bond atom")] {
            assert!(
                matches!(StereoBondEndCap::new(&input(&t,&c,&p),atom,other,1),Err(ReactionProductError::Invariant {detail:d,..}) if d==detail)
            );
        }
        t.adjacency = AdjacencyList::default();
        assert!(matches!(
            StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 1),
            Err(ReactionProductError::Invariant {
                detail: "atom adjacency row missing",
                ..
            })
        ));
    }
    #[test]
    fn only_selected_neighbor_pointer_is_resolved_and_source_degree_check_precedes_it() {
        let mut t = topology(4, &[(0, 1), (0, 3)]);
        t.atoms.truncate(2);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        assert!(
            !StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 3)
                .unwrap()
                .has_non_anchor()
        );
        assert!(matches!(
            StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 9),
            Err(ReactionProductError::Invariant {
                detail: "non-anchor atom row out of range",
                ..
            })
        ));
        t.atoms[0] = Atom::from_spec(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(2),
        );
        assert!(matches!(
            StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, 9),
            Err(ReactionProductError::Invariant {
                detail: "Stereo Bond extremes must have less than four neighbors",
                ..
            })
        ));
    }
    #[test]
    fn unsigned_anchor_boundary_preserves_maximum_and_rejects_unrepresentable_projection() {
        let t = topology(2, &[]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        assert_eq!(
            StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, u32::MAX as usize)
                .unwrap()
                .anchor_idx(),
            u32::MAX as usize
        );
        if usize::BITS > 32 {
            assert!(matches!(
                StereoBondEndCap::new(&input(&t, &c, &p), 0, 1, u32::MAX as usize + 1),
                Err(ReactionProductError::RowOverflow {
                    kind: "stereo anchor",
                    ..
                })
            ));
        }
    }
}

#[cfg(test)]
mod complete_forward_bond_stereo_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondSpec, CoordinateBlock, MoleculeProperties,
        PropertyValue, TopologyBlock,
    };
    use cosmolkit_types::Element;
    use std::collections::{BTreeMap, BTreeSet};
    pub(super) fn reactant(stereo: BondStereo) -> TopologyBlock {
        let edges = [(1, 2), (1, 0), (2, 3), (1, 4), (2, 5)];
        let bonds: Vec<_> = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                let spec = BondSpec::new(
                    AtomId::new(a),
                    AtomId::new(b),
                    if i == 0 {
                        BondOrder::Double
                    } else {
                        BondOrder::Single
                    },
                );
                Bond::from_spec(
                    BondId::new(i),
                    if i == 0 {
                        spec.with_stereo_atoms(AtomId::new(0), AtomId::new(3))
                            .with_stereo(stereo)
                    } else {
                        spec
                    },
                )
            })
            .collect();
        TopologyBlock {
            atoms: (0..6)
                .map(|i| {
                    Atom::from_spec(
                        AtomId::new(i),
                        AtomSpec::new(Element::C).with_no_implicit(true),
                    )
                })
                .collect(),
            adjacency: AdjacencyList::from_topology(6, &bonds),
            bonds,
            ..TopologyBlock::default()
        }
    }
    pub(super) fn product(reverse: bool) -> ProductBuilder {
        let mut t = reactant(BondStereo::Any);
        if reverse {
            t.bonds[0] = Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(2), AtomId::new(1), BondOrder::Double)
                    .with_stereo_atoms(AtomId::new(5), AtomId::new(4))
                    .with_stereo(BondStereo::Any),
            );
        }
        let neighbors = (0..6)
            .map(|i| t.adjacency.neighbors_of(i).to_vec())
            .collect();
        ProductBuilder {
            topology: t,
            neighbors,
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; 6],
            bond_origins: vec![None; 5],
        }
    }
    pub(super) fn mapping() -> ReactantProductMapping {
        ReactantProductMapping {
            mapped: vec![true; 6],
            skipped: vec![false; 6],
            reactant_to_product: (0..6).map(|i| (i, vec![i])).collect(),
            product_to_reactant: (0..6).map(|i| (i, i)).collect(),
            product_atom_bond: BTreeMap::new(),
            template_bonds: BTreeSet::new(),
        }
    }
    pub(super) fn input<'a>(
        t: &'a TopologyBlock,
        c: &'a CoordinateBlock,
        p: &'a MoleculeProperties,
    ) -> ReactionInput<'a> {
        ReactionInput {
            topology: t,
            coordinates: c,
            properties: p,
            rings: None,
            valence: None,
        }
    }
    #[test]
    fn all_defined_stereo_codes_forward_with_primary_secondary_parity_and_reversed_bond_ends() {
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        for stereo in [
            BondStereo::Z,
            BondStereo::E,
            BondStereo::Cis,
            BondStereo::Trans,
            BondStereo::AtropCw,
            BondStereo::AtropCcw,
        ] {
            for mask in 0..4 {
                for reverse in [false, true] {
                    let r = reactant(stereo);
                    let mut p = product(reverse);
                    let mut m = mapping();
                    if mask & 1 != 0 {
                        m.reactant_to_product.remove(&0);
                    }
                    if mask & 2 != 0 {
                        m.reactant_to_product.remove(&3);
                    }
                    forward_bond_stereo(
                        &mut p,
                        &input(&r, &c, &props),
                        &mut m,
                        BondId::new(0),
                        BondId::new(0),
                    )
                    .unwrap();
                    let left = if mask & 1 != 0 { 4 } else { 0 };
                    let right = if mask & 2 != 0 { 5 } else { 3 };
                    assert_eq!(
                        p.topology.bonds[0].stereo_atoms(),
                        Some(if reverse {
                            [AtomId::new(right), AtomId::new(left)]
                        } else {
                            [AtomId::new(left), AtomId::new(right)]
                        })
                    );
                    let flip = (mask == 1) || (mask == 2);
                    let cis = matches!(stereo, BondStereo::Z | BondStereo::Cis) != flip;
                    assert_eq!(
                        p.topology.bonds[0].stereo(),
                        if cis {
                            BondStereo::Cis
                        } else {
                            BondStereo::Trans
                        }
                    );
                }
            }
        }
    }
    #[test]
    fn absent_begin_map_key_inserts_zero_before_endpoint_mismatch_error() {
        let r = reactant(BondStereo::E);
        let mut p = product(false);
        let before = p.topology.bonds[0].clone();
        let mut m = mapping();
        m.product_to_reactant.remove(&1);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        assert!(matches!(
            forward_bond_stereo(
                &mut p,
                &input(&r, &c, &props),
                &mut m,
                BondId::new(0),
                BondId::new(0)
            ),
            Err(ReactionProductError::Invariant {
                detail: "Reactant and Product bond ends do not match",
                ..
            })
        ));
        assert_eq!(m.product_to_reactant[&1], 0);
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn default_zero_begin_mapping_is_a_real_source_endpoint_and_can_forward() {
        let mut r = reactant(BondStereo::E);
        let swap = |v: usize| {
            if v == 0 {
                1
            } else if v == 1 {
                0
            } else {
                v
            }
        };
        r.bonds = r
            .bonds
            .iter()
            .map(|b| {
                let spec = BondSpec::new(
                    AtomId::new(swap(b.begin().index())),
                    AtomId::new(swap(b.end().index())),
                    b.order(),
                );
                let spec = if b.id() == BondId::new(0) {
                    spec.with_stereo_atoms(AtomId::new(1), AtomId::new(3))
                        .with_stereo(BondStereo::E)
                } else {
                    spec
                };
                Bond::from_spec(b.id(), spec)
            })
            .collect();
        r.adjacency = AdjacencyList::from_topology(6, &r.bonds);
        let mut p = product(false);
        let mut m = mapping();
        m.product_to_reactant.remove(&1);
        m.reactant_to_product.insert(1, vec![0]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        forward_bond_stereo(
            &mut p,
            &input(&r, &c, &props),
            &mut m,
            BondId::new(0),
            BondId::new(0),
        )
        .unwrap();
        assert_eq!(m.product_to_reactant[&1], 0);
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::Trans);
    }
    #[test]
    fn present_empty_primary_or_absent_both_candidates_retains_entire_product_bond() {
        let r = reactant(BondStereo::Z);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        for empty_primary in [false, true] {
            let mut p = product(false);
            let before = p.topology.bonds[0].clone();
            let mut m = mapping();
            if empty_primary {
                m.reactant_to_product.insert(0, vec![]);
            } else {
                m.reactant_to_product.remove(&0);
                m.reactant_to_product.remove(&4);
            }
            forward_bond_stereo(
                &mut p,
                &input(&r, &c, &props),
                &mut m,
                BondId::new(0),
                BondId::new(0),
            )
            .unwrap();
            assert_eq!(p.topology.bonds[0], before);
        }
    }
    #[test]
    fn disconnected_start_warning_short_circuits_invalid_end_index_without_writes() {
        let r = reactant(BondStereo::E);
        let mut p = product(false);
        let before = p.topology.bonds[0].clone();
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![5]);
        m.reactant_to_product.insert(3, vec![99]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        forward_bond_stereo(
            &mut p,
            &input(&r, &c, &props),
            &mut m,
            BondId::new(0),
            BondId::new(0),
        )
        .unwrap();
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn connected_start_reaches_invalid_end_range_check_before_any_product_write() {
        let r = reactant(BondStereo::E);
        let mut p = product(false);
        let before = p.topology.bonds[0].clone();
        let mut m = mapping();
        m.reactant_to_product.insert(3, vec![99]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        assert!(
            matches!(forward_bond_stereo(&mut p,&input(&r,&c,&props),&mut m,BondId::new(0),BondId::new(0)),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn multiple_candidates_choose_first_connected_in_given_order_without_relabeling_swap() {
        let r = reactant(BondStereo::Z);
        let mut p = product(false);
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![5, 4, 0]);
        m.reactant_to_product.insert(3, vec![4, 5, 3]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        forward_bond_stereo(
            &mut p,
            &input(&r, &c, &props),
            &mut m,
            BondId::new(0),
            BondId::new(0),
        )
        .unwrap();
        assert_eq!(
            p.topology.bonds[0].stereo_atoms(),
            Some([AtomId::new(4), AtomId::new(5)])
        );
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::Cis);
    }
    #[test]
    fn undefined_stereo_precondition_precedes_lazy_valence_or_mapping_work() {
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        for stereo in [BondStereo::None, BondStereo::Any] {
            let mut r = reactant(stereo);
            r.atoms[1] = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C));
            let mut p = product(false);
            let mut m = mapping();
            m.product_to_reactant.remove(&1);
            assert!(matches!(
                forward_bond_stereo(
                    &mut p,
                    &input(&r, &c, &props),
                    &mut m,
                    BondId::new(0),
                    BondId::new(0)
                ),
                Err(ReactionProductError::Invariant {
                    detail: "bond in reactant must have defined stereo",
                    ..
                })
            ));
            assert!(!m.product_to_reactant.contains_key(&1));
        }
    }
    #[test]
    fn no_source_reference_atoms_clears_only_stereo_before_constructor_and_map_access() {
        let mut r = reactant(BondStereo::AtropCw);
        r.bonds[0].set_stereo_atoms(None);
        r.atoms[1] = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C));
        let mut p = product(false);
        let refs = p.topology.bonds[0].stereo_atoms();
        let mut m = mapping();
        m.product_to_reactant.clear();
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        forward_bond_stereo(
            &mut p,
            &input(&r, &c, &props),
            &mut m,
            BondId::new(0),
            BondId::new(0),
        )
        .unwrap();
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::None);
        assert_eq!(p.topology.bonds[0].stereo_atoms(), refs);
        assert!(m.product_to_reactant.is_empty());
    }
    #[test]
    fn stored_references_skip_bad_cip_rank_but_reached_constructor_fact_error_precedes_map_insertion()
     {
        let mut r = reactant(BondStereo::E);
        r.atoms[0]
            .set_prop("_CIPRank", PropertyValue::String("not-a-rank".into()))
            .unwrap();
        let mut p = product(false);
        let mut m = mapping();
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        forward_bond_stereo(
            &mut p,
            &input(&r, &c, &props),
            &mut m,
            BondId::new(0),
            BondId::new(0),
        )
        .unwrap();
        r.atoms[1] = Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C));
        m.product_to_reactant.remove(&1);
        assert!(matches!(
            forward_bond_stereo(
                &mut p,
                &input(&r, &c, &props),
                &mut m,
                BondId::new(0),
                BondId::new(0)
            ),
            Err(ReactionProductError::StereoGetter { .. })
        ));
        assert!(!m.product_to_reactant.contains_key(&1));
    }
}

#[cfg(test)]
mod complete_translate_directions_source_tests {
    use super::*;
    use cosmolkit_model::{Atom, AtomSpec, Bond, BondSpec, NeighborRef, TopologyBlock};
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;
    fn product(
        start_dir: BondDirection,
        end_dir: BondDirection,
        start_reverse: bool,
        end_reverse: bool,
    ) -> ProductBuilder {
        let (a, b) = if start_reverse { (0, 1) } else { (1, 0) };
        let (c, d) = if end_reverse { (3, 2) } else { (2, 3) };
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double)
                    .with_stereo_atoms(AtomId::new(0), AtomId::new(3))
                    .with_stereo(BondStereo::Any),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single)
                    .with_direction(start_dir),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(AtomId::new(c), AtomId::new(d), BondOrder::Single)
                    .with_direction(end_dir),
            ),
        ];
        let mut neighbors = vec![vec![]; 4];
        for bond in &bonds {
            neighbors[bond.begin().index()].push(NeighborRef {
                atom_index: bond.end().index(),
                bond: bond.id(),
            });
            neighbors[bond.end().index()].push(NeighborRef {
                atom_index: bond.begin().index(),
                bond: bond.id(),
            });
        }
        ProductBuilder {
            topology: TopologyBlock {
                atoms: (0..4)
                    .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                    .collect(),
                bonds,
                ..TopologyBlock::default()
            },
            neighbors,
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; 4],
            bond_origins: vec![None; 3],
        }
    }
    #[test]
    fn sixteen_direction_and_endpoint_orientations_follow_source_two_toggles() {
        for a in [BondDirection::EndDownRight, BondDirection::EndUpRight] {
            for b in [BondDirection::EndDownRight, BondDirection::EndUpRight] {
                for ar in [false, true] {
                    for br in [false, true] {
                        let mut p = product(a, b, ar, br);
                        translate_directions(
                            &mut p,
                            BondId::new(0),
                            BondId::new(1),
                            BondId::new(2),
                        )
                        .unwrap();
                        let trans = (a == b) ^ (!ar) ^ br;
                        assert_eq!(
                            p.topology.bonds[0].stereo(),
                            if trans {
                                BondStereo::Trans
                            } else {
                                BondStereo::Cis
                            }
                        );
                        assert_eq!(
                            p.topology.bonds[0].stereo_atoms(),
                            Some([AtomId::new(0), AtomId::new(3)])
                        );
                    }
                }
            }
        }
    }
    #[test]
    fn all_non_stereo_directions_fail_before_anchor_graph_access_or_writes() {
        for direction in [
            BondDirection::None,
            BondDirection::BeginWedge,
            BondDirection::BeginDash,
            BondDirection::Unknown,
            BondDirection::EitherDouble,
        ] {
            for bad_start in [false, true] {
                let good = BondDirection::EndUpRight;
                let mut p = product(
                    if bad_start { direction } else { good },
                    if bad_start { good } else { direction },
                    false,
                    false,
                );
                p.neighbors.clear();
                let before = p.topology.bonds[0].clone();
                assert!(matches!(
                    translate_directions(&mut p, BondId::new(0), BondId::new(1), BondId::new(2)),
                    Err(ReactionProductError::Invariant {
                        detail: "Both neighboring bonds must have bond directions",
                        ..
                    })
                ));
                assert_eq!(p.topology.bonds[0], before);
            }
        }
    }
    #[test]
    fn nonincident_start_or_end_reports_source_bad_index_before_stereo_write() {
        for row in [1, 2] {
            let mut p = product(
                BondDirection::EndUpRight,
                BondDirection::EndDownRight,
                false,
                false,
            );
            p.topology.bonds[row] = Bond::from_spec(
                BondId::new(row),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single)
                    .with_direction(BondDirection::EndUpRight),
            );
            let before = p.topology.bonds[0].clone();
            assert!(matches!(
                translate_directions(&mut p, BondId::new(0), BondId::new(1), BondId::new(2)),
                Err(ReactionProductError::Invariant {
                    stage: "Bond::getOtherAtomIdx",
                    detail: "bad index",
                    ..
                })
            ));
            assert_eq!(p.topology.bonds[0], before);
        }
    }
    #[test]
    fn start_connectedness_failure_precedes_end_and_retains_previous_reference_atoms() {
        let mut p = product(
            BondDirection::EndUpRight,
            BondDirection::EndDownRight,
            false,
            false,
        );
        p.neighbors[1].retain(|n| n.atom_index != 0);
        p.neighbors[2].clear();
        let before = p.topology.bonds[0].clone();
        assert!(matches!(
            translate_directions(&mut p, BondId::new(0), BondId::new(1), BondId::new(2)),
            Err(ReactionProductError::Invariant {
                detail: "bgnIdx not connected to begin atom of bond",
                ..
            })
        ));
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn connected_start_reaches_end_connectedness_and_range_errors_before_writes() {
        let mut p = product(
            BondDirection::EndUpRight,
            BondDirection::EndDownRight,
            false,
            false,
        );
        p.neighbors[2].retain(|n| n.atom_index != 3);
        let before = p.topology.bonds[0].clone();
        assert!(matches!(
            translate_directions(&mut p, BondId::new(0), BondId::new(1), BondId::new(2)),
            Err(ReactionProductError::Invariant {
                detail: "endIdx not connected to end atom of bond",
                ..
            })
        ));
        assert_eq!(p.topology.bonds[0], before);
        p.topology.bonds[2] = Bond::from_spec(
            BondId::new(2),
            BondSpec::new(AtomId::new(2), AtomId::new(99), BondOrder::Single)
                .with_direction(BondDirection::EndUpRight),
        );
        assert!(
            matches!(translate_directions(&mut p,BondId::new(0),BondId::new(1),BondId::new(2)),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn source_does_not_require_central_or_neighbor_bond_order_to_be_double_or_single() {
        let mut p = product(
            BondDirection::EndUpRight,
            BondDirection::EndUpRight,
            false,
            false,
        );
        p.topology.bonds[0].set_order(BondOrder::Single);
        p.topology.bonds[1].set_order(BondOrder::Triple);
        translate_directions(&mut p, BondId::new(0), BondId::new(1), BondId::new(2)).unwrap();
        assert_eq!(p.topology.bonds[0].order(), BondOrder::Single);
        assert_eq!(p.topology.bonds[1].order(), BondOrder::Triple);
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::Cis);
    }
    #[test]
    fn row_projection_errors_are_structural_and_product_pointer_check_is_first() {
        let mut p = product(
            BondDirection::EndUpRight,
            BondDirection::EndUpRight,
            false,
            false,
        );
        assert!(matches!(
            translate_directions(&mut p, BondId::new(99), BondId::new(98), BondId::new(97)),
            Err(ReactionProductError::Invariant {
                detail: "no product bond",
                ..
            })
        ));
        assert!(matches!(
            translate_directions(&mut p, BondId::new(0), BondId::new(98), BondId::new(97)),
            Err(ReactionProductError::Invariant {
                detail: "Both neighboring bonds must have bond directions",
                ..
            })
        ));
    }
}

#[cfg(test)]
mod complete_update_stereo_bonds_source_tests {
    use super::complete_forward_bond_stereo_source_tests::{input, mapping, product, reactant};
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, CoordinateBlock, MoleculeProperties, NeighborRef, PropertyValue,
    };
    fn run(
        p: &mut ProductBuilder,
        r: &cosmolkit_model::TopologyBlock,
        m: &mut ReactantProductMapping,
    ) -> Result<(), ReactionProductError> {
        update_stereo_bonds(
            p,
            &input(
                r,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            m,
        )
    }
    #[test]
    fn non_double_bond_skips_unknown_property_adjacency_and_input_mapping_access() {
        let mut p = product(false);
        p.topology.bonds[0].set_order(BondOrder::Single);
        p.topology.bonds[0]
            .set_prop("_UnknownStereoRxnBond", PropertyValue::Int(1))
            .unwrap();
        p.neighbors.clear();
        let before = p.topology.clone();
        let mut m = mapping();
        m.product_to_reactant.clear();
        run(&mut p, &reactant(BondStereo::E), &mut m).unwrap();
        assert_eq!(p.topology, before);
    }
    #[test]
    fn unknown_property_presence_clears_stereo_and_property_without_reading_its_value_or_graph() {
        for value in [
            PropertyValue::Bool(false),
            PropertyValue::Int(0),
            PropertyValue::String("any-tag".into()),
        ] {
            let mut p = product(false);
            p.topology.bonds[0]
                .set_prop("_UnknownStereoRxnBond", value)
                .unwrap();
            let refs = p.topology.bonds[0].stereo_atoms();
            p.neighbors.clear();
            run(&mut p, &reactant(BondStereo::E), &mut mapping()).unwrap();
            assert_eq!(p.topology.bonds[0].stereo(), BondStereo::None);
            assert_eq!(p.topology.bonds[0].stereo_atoms(), refs);
            assert!(p.topology.bonds[0].prop("_UnknownStereoRxnBond").is_none());
        }
    }
    #[test]
    fn reached_clear_error_retains_stereo_none_prefix_and_original_property_then_stops() {
        let mut p = product(false);
        p.topology.bonds[0]
            .set_prop("_UnknownStereoRxnBond", PropertyValue::Int(1))
            .unwrap();
        p.topology.bonds[0]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        p.neighbors.clear();
        assert!(matches!(
            run(&mut p, &reactant(BondStereo::E), &mut mapping()),
            Err(ReactionProductError::BondValue(_))
        ));
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::None);
        assert_eq!(
            p.topology.bonds[0].prop("_UnknownStereoRxnBond"),
            Some(&PropertyValue::Int(1))
        );
    }
    #[test]
    fn absent_unknown_property_does_not_read_bad_computed_membership_metadata() {
        let mut p = product(false);
        p.topology.bonds[0]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        let before = p.topology.bonds[0].clone();
        let mut m = mapping();
        m.product_to_reactant.clear();
        run(&mut p, &reactant(BondStereo::E), &mut m).unwrap();
        assert_eq!(p.topology.bonds[0], before);
    }
    #[test]
    fn two_product_directions_take_priority_over_all_reactant_mapping_and_fact_work() {
        let mut p = product(false);
        p.topology.bonds[1].set_direction(BondDirection::EndUpRight);
        p.topology.bonds[2].set_direction(BondDirection::EndUpRight);
        let mut m = mapping();
        m.product_to_reactant.clear();
        let r = cosmolkit_model::TopologyBlock::default();
        run(&mut p, &r, &mut m).unwrap();
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::Cis);
    }
    #[test]
    fn end_directed_neighbor_is_evaluated_even_when_start_has_no_match() {
        let mut p = product(false);
        p.neighbors[2].push(NeighborRef {
            atom_index: 9,
            bond: BondId::new(99),
        });
        let mut m = mapping();
        m.product_to_reactant.clear();
        assert!(matches!(
            run(&mut p, &reactant(BondStereo::E), &mut m),
            Err(ReactionProductError::Invariant {
                detail: "incident bond row missing",
                ..
            })
        ));
    }
    #[test]
    fn source_neighbor_helper_stops_at_first_directed_non_double_and_skips_later_bad_rows() {
        let mut p = product(false);
        p.topology.bonds[1].set_direction(BondDirection::EndUpRight);
        p.topology.bonds[2].set_direction(BondDirection::EndDownRight);
        p.neighbors[1].push(NeighborRef {
            atom_index: 9,
            bond: BondId::new(99),
        });
        p.neighbors[2].push(NeighborRef {
            atom_index: 9,
            bond: BondId::new(98),
        });
        run(&mut p, &reactant(BondStereo::E), &mut mapping()).unwrap();
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::Trans);
    }
    #[test]
    fn missing_endpoint_mapping_skips_source_graph_access_without_inserting_keys() {
        for absent_begin in [false, true] {
            let mut p = product(false);
            let before = p.topology.clone();
            let mut m = mapping();
            if absent_begin {
                m.product_to_reactant.remove(&1);
                m.product_to_reactant.insert(2, 99);
            } else {
                m.product_to_reactant.remove(&2);
                m.product_to_reactant.insert(1, 99);
            }
            let before_map = m.product_to_reactant.clone();
            run(&mut p, &cosmolkit_model::TopologyBlock::default(), &mut m).unwrap();
            assert_eq!(p.topology, before);
            assert_eq!(m.product_to_reactant, before_map);
        }
    }
    #[test]
    fn reached_source_range_and_missing_csr_are_errors_instead_of_missing_bond_fallback() {
        let mut p = product(false);
        let mut m = mapping();
        m.product_to_reactant.insert(1, 99);
        assert!(
            matches!(run(&mut p,&reactant(BondStereo::E),&mut m),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        let mut r = reactant(BondStereo::E);
        r.adjacency = AdjacencyList::default();
        assert!(matches!(
            run(&mut p, &r, &mut mapping()),
            Err(ReactionProductError::Invariant {
                detail: "reactant adjacency row missing",
                ..
            })
        ));
    }
    #[test]
    fn absent_or_non_double_or_no_stereo_source_bond_retains_existing_product_stereo() {
        for case in 0..3 {
            let mut p = product(false);
            p.topology.bonds[0].set_stereo(BondStereo::Cis).unwrap();
            let before = p.topology.bonds[0].clone();
            let mut r = reactant(BondStereo::None);
            if case == 0 {
                r.adjacency = AdjacencyList::from_topology(6, &[]);
            } else if case == 1 {
                r.bonds[0].set_order(BondOrder::Single);
            }
            run(&mut p, &r, &mut mapping()).unwrap();
            assert_eq!(p.topology.bonds[0], before);
        }
    }
    #[test]
    fn source_any_changes_only_label_while_defined_stereo_dispatches_full_forwarder() {
        for stereo in [
            BondStereo::Any,
            BondStereo::E,
            BondStereo::Z,
            BondStereo::Cis,
            BondStereo::Trans,
            BondStereo::AtropCw,
            BondStereo::AtropCcw,
        ] {
            let mut p = product(false);
            p.topology.bonds[0].set_stereo(BondStereo::Cis).unwrap();
            run(&mut p, &reactant(stereo), &mut mapping()).unwrap();
            assert_eq!(
                p.topology.bonds[0].stereo(),
                match stereo {
                    BondStereo::Any => BondStereo::Any,
                    BondStereo::Z | BondStereo::Cis => BondStereo::Cis,
                    _ => BondStereo::Trans,
                }
            );
        }
    }
    #[test]
    fn earlier_bond_writes_remain_when_a_later_double_bond_fails_its_clear() {
        let mut p = product(false);
        p.topology.bonds[0]
            .set_prop("_UnknownStereoRxnBond", PropertyValue::Int(1))
            .unwrap();
        p.topology.bonds[1].set_order(BondOrder::Double);
        p.topology.bonds[1]
            .set_prop("_UnknownStereoRxnBond", PropertyValue::Int(1))
            .unwrap();
        p.topology.bonds[1]
            .set_prop("__computedProps", PropertyValue::Bool(false))
            .unwrap();
        assert!(run(&mut p, &reactant(BondStereo::E), &mut mapping()).is_err());
        assert!(p.topology.bonds[0].prop("_UnknownStereoRxnBond").is_none());
        assert_eq!(p.topology.bonds[0].stereo(), BondStereo::None);
        assert!(p.topology.bonds[1].prop("_UnknownStereoRxnBond").is_some());
        assert_eq!(p.topology.bonds[1].stereo(), BondStereo::None);
    }
}

#[cfg(test)]
mod complete_correct_matched_chirality_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondSpec, CoordinateBlock, MoleculeProperties,
        PropertyValue, TopologyBlock,
    };
    use cosmolkit_types::Element;
    use std::collections::{BTreeMap, BTreeSet};
    pub(super) fn topology(neighbors: &[usize], tag: ChiralTag) -> TopologyBlock {
        let bonds: Vec<_> = neighbors
            .iter()
            .enumerate()
            .map(|(i, &n)| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(0), AtomId::new(n), BondOrder::Single),
                )
            })
            .collect();
        let mut atoms: Vec<_> = (0..10)
            .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
            .collect();
        atoms[0].set_chiral_tag(tag);
        TopologyBlock {
            atoms,
            adjacency: AdjacencyList::from_topology(10, &bonds),
            bonds,
            ..TopologyBlock::default()
        }
    }
    pub(super) fn product(neighbors: &[usize], tag: ChiralTag) -> ProductBuilder {
        let t = topology(neighbors, tag);
        let neighbors = (0..10)
            .map(|i| t.adjacency.neighbors_of(i).to_vec())
            .collect();
        ProductBuilder {
            bond_origins: vec![None; t.bonds.len()],
            topology: t,
            neighbors,
            bookmarks: BTreeMap::new(),
            atom_origins: vec![None; 10],
        }
    }
    pub(super) fn mapping() -> ReactantProductMapping {
        ReactantProductMapping {
            mapped: vec![true; 10],
            skipped: vec![false; 10],
            reactant_to_product: BTreeMap::from([(0, vec![0])]),
            product_to_reactant: (0..10).map(|i| (i, i)).collect(),
            product_atom_bond: BTreeMap::new(),
            template_bonds: BTreeSet::new(),
        }
    }
    fn run(
        r: &TopologyBlock,
        p: &mut ProductBuilder,
        m: &mut ReactantProductMapping,
    ) -> Result<(), ReactionProductError> {
        correct_matched_chirality(
            p,
            &ReactionInput {
                topology: r,
                coordinates: &CoordinateBlock::default(),
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            m,
            0,
        )
    }
    #[test]
    fn permutation_parity_and_source_inversion_flags_combine_by_xor() {
        for (order, odd) in [
            (vec![1, 2, 3], false),
            (vec![2, 1, 3], true),
            (vec![3, 1, 2], false),
        ] {
            for flag in [None, Some(-1), Some(0), Some(1), Some(2)] {
                for tag in [ChiralTag::TetrahedralCw, ChiralTag::TetrahedralCcw] {
                    let r = topology(&[1, 2, 3], tag);
                    let mut p = product(&order, ChiralTag::Other);
                    if let Some(flag) = flag {
                        p.topology.atoms[0]
                            .set_prop("molInversionFlag", PropertyValue::Int(flag))
                            .unwrap();
                    }
                    run(&r, &mut p, &mut mapping()).unwrap();
                    let invert = odd ^ (flag == Some(1));
                    assert_eq!(
                        p.topology.atoms[0].chiral_tag(),
                        if invert {
                            if tag == ChiralTag::TetrahedralCw {
                                ChiralTag::TetrahedralCcw
                            } else {
                                ChiralTag::TetrahedralCw
                            }
                        } else {
                            tag
                        }
                    );
                }
            }
        }
    }
    #[test]
    fn source_tag_and_flag_guards_skip_missing_adjacency_but_flag_cast_is_first() {
        for (tag, flag) in [
            (ChiralTag::Unspecified, 0),
            (ChiralTag::Other, 0),
            (ChiralTag::TetrahedralCw, 3),
        ] {
            let mut r = topology(&[1, 2, 3], tag);
            r.adjacency = AdjacencyList::default();
            let mut p = product(&[1, 2, 3], ChiralTag::Allene);
            p.neighbors.clear();
            p.topology.atoms[0]
                .set_prop("molInversionFlag", PropertyValue::Int(flag))
                .unwrap();
            run(&r, &mut p, &mut mapping()).unwrap();
            assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Allene);
            p.topology.atoms[0]
                .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
                .unwrap();
            assert!(matches!(
                run(&r, &mut p, &mut mapping()),
                Err(ReactionProductError::PropertyInt { .. })
            ));
        }
    }
    #[test]
    fn degree_guards_preserve_source_short_circuit_and_skip_chiral_assignment() {
        let r = topology(&[1, 2], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 3], ChiralTag::Other);
        p.neighbors.clear();
        run(&r, &mut p, &mut mapping()).unwrap();
        for (rn, pn) in [
            (vec![1, 2, 3], vec![1, 2]),
            (vec![1, 2, 3], vec![1, 2, 3, 4, 5]),
        ] {
            let r = topology(&rn, ChiralTag::TetrahedralCw);
            let mut p = product(&pn, ChiralTag::Other);
            run(&r, &mut p, &mut mapping()).unwrap();
            assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Other);
        }
    }
    #[test]
    fn one_unknown_same_degree_substitutes_first_unmatched_source_bond_in_place() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[3, 2, 1], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.remove(&3);
        run(&r, &mut p, &mut m).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn one_extra_unknown_product_bond_is_removed_before_permutation() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 4, 3, 2], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.remove(&4);
        run(&r, &mut p, &mut m).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn lost_reactant_bond_is_removed_in_source_order_only_until_lengths_agree() {
        let r = topology(&[1, 2, 3, 4], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 4, 2], ChiralTag::Other);
        run(&r, &mut p, &mut mapping()).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn two_unknowns_break_before_later_bad_mapping_and_one_unknown_with_lost_bond_skips() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 3, 4], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.remove(&1);
        m.product_to_reactant.remove(&2);
        m.product_to_reactant.insert(3, 99);
        run(&r, &mut p, &mut m).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Other);
        let r = topology(&[1, 2, 3, 4], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 5], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.remove(&5);
        run(&r, &mut p, &mut m).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Other);
    }
    #[test]
    fn missing_mapped_key_is_value_initialized_empty_without_any_product_read() {
        let r = topology(&[], ChiralTag::TetrahedralCw);
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms.clear();
        p.neighbors.clear();
        let mut m = mapping();
        m.reactant_to_product.clear();
        run(&r, &mut p, &mut m).unwrap();
        assert_eq!(m.reactant_to_product.get(&0), Some(&vec![]));
    }
    #[test]
    fn reached_original_edge_range_error_does_not_become_an_unknown_bond() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 3], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.insert(1, 99);
        assert!(
            matches!(run(&r,&mut p,&mut m),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Other);
    }
    #[test]
    fn chiral_tag_assignment_precedes_permutation_failure_without_swallowing_error() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 3], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.insert(2, 1);
        m.product_to_reactant.insert(3, 1);
        assert!(run(&r, &mut p, &mut m).is_err());
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
    }
    #[test]
    fn signed_unmatched_bond_guard_rejects_high_bit_source_index_before_substitution() {
        let mut r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        r.bonds[2] = Bond::from_spec(
            BondId::new(0x8000_0000),
            BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
        );
        let mut p = product(&[1, 2, 3], ChiralTag::Other);
        let mut m = mapping();
        m.product_to_reactant.remove(&3);
        assert!(matches!(
            run(&r, &mut p, &mut m),
            Err(ReactionProductError::Invariant {
                detail: "extra unmapped atom",
                ..
            })
        ));
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::Other);
    }
    #[test]
    fn one_to_many_rows_preserve_earlier_chiral_write_when_later_flag_read_fails() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[2, 1, 3], ChiralTag::Other);
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![0, 6]);
        p.topology.atoms[6].set_chiral_tag(ChiralTag::Allene);
        p.topology.atoms[6]
            .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(run(&r, &mut p, &mut m).is_err());
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(p.topology.atoms[6].chiral_tag(), ChiralTag::Allene);
    }
}

#[cfg(test)]
mod complete_correct_unmatched_chirality_source_tests {
    use super::complete_correct_matched_chirality_source_tests::{mapping, product, topology};
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Bond, BondSpec, CoordinateBlock, MoleculeProperties, NeighborRef,
        PropertyValue,
    };
    fn mapping_with_neighbors() -> ReactantProductMapping {
        let mut m = mapping();
        for i in 1..10 {
            m.reactant_to_product.insert(i, vec![i]);
        }
        m
    }
    fn run(
        r: &cosmolkit_model::TopologyBlock,
        p: &mut ProductBuilder,
        m: &mut ReactantProductMapping,
        rows: &[usize],
    ) -> Result<(), ReactionProductError> {
        correct_unmatched_chirality(
            p,
            &ReactionInput {
                topology: r,
                coordinates: &CoordinateBlock::default(),
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            m,
            rows,
        )
    }
    #[test]
    fn physical_product_bond_order_uses_source_perturbation_parity_for_both_cw_and_ccw() {
        for (order, odd) in [
            (vec![1, 2, 3], false),
            (vec![2, 1, 3], true),
            (vec![3, 1, 2], false),
        ] {
            for tag in [ChiralTag::TetrahedralCw, ChiralTag::TetrahedralCcw] {
                let r = topology(&[1, 2, 3], tag);
                let mut p = product(&order, tag);
                run(&r, &mut p, &mut mapping_with_neighbors(), &[0]).unwrap();
                assert_eq!(
                    p.topology.atoms[0].chiral_tag(),
                    if odd {
                        if tag == ChiralTag::TetrahedralCw {
                            ChiralTag::TetrahedralCcw
                        } else {
                            ChiralTag::TetrahedralCw
                        }
                    } else {
                        tag
                    }
                );
            }
        }
    }
    #[test]
    fn degree_change_clears_only_chiral_tag_for_every_nonunspecified_source_tag() {
        for tag in [
            ChiralTag::TetrahedralCw,
            ChiralTag::TetrahedralCcw,
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            let r = topology(&[1, 2, 3], tag);
            let mut p = product(&[1, 2], tag);
            p.topology.atoms[0].set_chiral_permutation(Some(17));
            p.topology.atoms[0]
                .set_prop("molInversionFlag", PropertyValue::String("ignored".into()))
                .unwrap();
            let mut expected = p.topology.atoms[0].clone();
            expected.set_chiral_tag(ChiralTag::Unspecified);
            run(&r, &mut p, &mut mapping_with_neighbors(), &[0]).unwrap();
            assert_eq!(p.topology.atoms[0], expected);
        }
    }
    #[test]
    fn equal_degree_non_cw_ccw_tags_skip_source_bonds_other_mapping_and_inversion_metadata() {
        for tag in [
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            let mut r = topology(&[1, 2, 3], tag);
            r.bonds.clear();
            let mut p = product(&[2, 1, 3], tag);
            p.topology.atoms[0]
                .set_prop("molInversionFlag", PropertyValue::String("ignored".into()))
                .unwrap();
            let before = p.topology.clone();
            run(&r, &mut p, &mut mapping(), &[0]).unwrap();
            assert_eq!(p.topology, before);
        }
    }
    #[test]
    fn missing_reactant_chirality_precedes_degree_or_mapping_and_mismatched_product_precedes_degree()
     {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.adjacency = AdjacencyList::default();
        let mut p = product(&[], ChiralTag::Other);
        let mut m = mapping();
        m.reactant_to_product.clear();
        assert!(matches!(
            run(&r, &mut p, &mut m, &[0]),
            Err(ReactionProductError::Invariant {
                detail: "missing atom chirality.",
                ..
            })
        ));
        assert!(m.reactant_to_product.is_empty());
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        p.neighbors.clear();
        assert!(matches!(
            run(&r, &mut p, &mut mapping(), &[0]),
            Err(ReactionProductError::Invariant {
                detail: "invalid product chirality.",
                ..
            })
        ));
    }
    #[test]
    fn absent_center_mapping_inserts_empty_after_reactant_degree_read_and_skips_product() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms.clear();
        let mut m = mapping();
        m.reactant_to_product.clear();
        run(&r, &mut p, &mut m, &[0]).unwrap();
        assert_eq!(m.reactant_to_product.get(&0), Some(&vec![]));
        let mut r = r;
        r.adjacency = AdjacencyList::default();
        let mut m = mapping();
        m.reactant_to_product.clear();
        assert!(run(&r, &mut p, &mut m, &[0]).is_err());
        assert!(m.reactant_to_product.is_empty());
    }
    #[test]
    fn other_neighbor_index_comes_from_actual_reactant_bond_not_csr_atom_projection() {
        let mut r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        r.adjacency = AdjacencyList::from_topology(
            10,
            &[
                Bond::from_spec(
                    BondId::new(0),
                    BondSpec::new(AtomId::new(0), AtomId::new(7), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(1),
                    BondSpec::new(AtomId::new(0), AtomId::new(8), BondOrder::Single),
                ),
                Bond::from_spec(
                    BondId::new(2),
                    BondSpec::new(AtomId::new(0), AtomId::new(9), BondOrder::Single),
                ),
            ],
        );
        let mut p = product(&[2, 1, 3], ChiralTag::TetrahedralCw);
        run(&r, &mut p, &mut mapping_with_neighbors(), &[0]).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn nonincident_source_bond_fails_checked_other_atom_helper_before_mapping() {
        let mut r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        r.bonds[0] = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
        );
        let mut p = product(&[1, 2, 3], ChiralTag::TetrahedralCw);
        assert!(matches!(
            run(&r, &mut p, &mut mapping(), &[0]),
            Err(ReactionProductError::Invariant {
                stage: "Bond::getOtherAtomIdx",
                detail: "bad index",
                ..
            })
        ));
    }
    #[test]
    fn missing_neighbor_key_and_short_copy_vector_remain_distinct_structural_failures() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        for present in [false, true] {
            let mut p = product(&[1, 2, 3], ChiralTag::TetrahedralCw);
            let mut m = mapping_with_neighbors();
            if present {
                m.reactant_to_product.insert(1, vec![]);
            } else {
                m.reactant_to_product.remove(&1);
            }
            assert!(
                matches!(run(&r,&mut p,&mut m,&[0]),Err(ReactionProductError::Invariant {detail,..}) if detail==if present {"other atom product copy missing"}else{"other atom from bond not mapped."})
            );
            assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        }
    }
    #[test]
    fn missing_product_edge_or_out_of_range_neighbor_is_not_silently_ignored() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut m = mapping_with_neighbors();
        m.reactant_to_product.insert(1, vec![4]);
        assert!(matches!(
            run(&r, &mut p, &mut m, &[0]),
            Err(ReactionProductError::Invariant {
                detail: "no matching bond found in product",
                ..
            })
        ));
        m.reactant_to_product.insert(1, vec![99]);
        assert!(
            matches!(run(&r,&mut p,&mut m,&[0]),Err(ReactionProductError::TopologyEdit(cosmolkit_model::TopologyEditError::AtomOutOfRange {atom,..})) if atom==AtomId::new(99))
        );
    }
    #[test]
    fn one_to_many_uses_aligned_copy_index_and_preserves_earlier_inversion_on_later_error() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[2, 1, 3], ChiralTag::TetrahedralCw);
        p.topology.atoms[6].set_chiral_tag(ChiralTag::TetrahedralCw);
        p.neighbors[6] = vec![
            NeighborRef {
                atom_index: 7,
                bond: BondId::new(3),
            },
            NeighborRef {
                atom_index: 8,
                bond: BondId::new(4),
            },
            NeighborRef {
                atom_index: 9,
                bond: BondId::new(5),
            },
        ];
        for (i, n) in [7, 8, 9].into_iter().enumerate() {
            p.topology.bonds.push(Bond::from_spec(
                BondId::new(i + 3),
                BondSpec::new(AtomId::new(6), AtomId::new(n), BondOrder::Single),
            ));
        }
        let mut m = mapping_with_neighbors();
        m.reactant_to_product.insert(0, vec![0, 6]);
        m.reactant_to_product.insert(2, vec![2, 8]);
        m.reactant_to_product.insert(3, vec![3, 9]);
        assert!(matches!(
            run(&r, &mut p, &mut m, &[0]),
            Err(ReactionProductError::Invariant {
                detail: "other atom product copy missing",
                ..
            })
        ));
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
        assert_eq!(p.topology.atoms[6].chiral_tag(), ChiralTag::TetrahedralCw);
    }
    #[test]
    fn repeated_source_atom_is_processed_again_and_observes_first_inversion() {
        let r = topology(&[1, 2, 3], ChiralTag::TetrahedralCw);
        let mut p = product(&[2, 1, 3], ChiralTag::TetrahedralCw);
        assert!(matches!(
            run(&r, &mut p, &mut mapping_with_neighbors(), &[0, 0]),
            Err(ReactionProductError::Invariant {
                detail: "invalid product chirality.",
                ..
            })
        ));
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn zero_degree_source_cw_atom_has_empty_perturbation_without_added_degree_guard() {
        let r = topology(&[], ChiralTag::TetrahedralCw);
        let mut p = product(&[], ChiralTag::TetrahedralCw);
        run(&r, &mut p, &mut mapping(), &[0]).unwrap();
        assert_eq!(p.topology.atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
    }
}

#[cfg(test)]
mod complete_copy_enhanced_groups_source_tests {
    use super::complete_correct_matched_chirality_source_tests::{mapping, product, topology};
    use super::*;
    use cosmolkit_model::{CoordinateBlock, MoleculeProperties, PropertyValue, StereoGroupKind};
    fn group(
        kind: StereoGroupKind,
        atoms: &[usize],
        bonds: &[usize],
        id: Option<u32>,
        write: u32,
    ) -> StereoGroup {
        let g = StereoGroup::new(
            kind,
            atoms.iter().map(|&i| AtomId::new(i)).collect(),
            bonds.iter().map(|&i| BondId::new(i)).collect(),
        )
        .with_write_id(write);
        if let Some(id) = id { g.with_id(id) } else { g }
    }
    fn run(
        r: &cosmolkit_model::TopologyBlock,
        p: &mut ProductBuilder,
        m: &ReactantProductMapping,
    ) -> Result<(), ReactionProductError> {
        copy_enhanced_groups(
            p,
            &ReactionInput {
                topology: r,
                coordinates: &CoordinateBlock::default(),
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            m,
        )
    }
    #[test]
    fn mapped_group_members_preserve_source_atom_copy_order_duplicates_and_literal_flag_four_filter()
     {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![group(StereoGroupKind::Or, &[0, 1, 0], &[7], Some(13), 77)];
        let mut p = product(&[], ChiralTag::Other);
        for a in &mut p.topology.atoms {
            a.set_chiral_tag(ChiralTag::TetrahedralCw);
        }
        p.topology.atoms[2].set_chiral_tag(ChiralTag::Unspecified);
        p.topology.atoms[2]
            .set_prop(
                "molInversionFlag",
                PropertyValue::String("unreached".into()),
            )
            .unwrap();
        p.topology.atoms[3]
            .set_prop("molInversionFlag", PropertyValue::Int(4))
            .unwrap();
        p.topology.atoms[4]
            .set_prop("molInversionFlag", PropertyValue::Int(3))
            .unwrap();
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![4, 2, 3, 4]);
        m.reactant_to_product.insert(1, vec![0]);
        run(&r, &mut p, &m).unwrap();
        assert_eq!(p.topology.stereo_groups.len(), 1);
        let g = &p.topology.stereo_groups[0];
        assert_eq!(
            g.atoms(),
            [
                AtomId::new(4),
                AtomId::new(4),
                AtomId::new(0),
                AtomId::new(4),
                AtomId::new(4)
            ]
        );
        assert!(g.bonds().is_empty());
        assert_eq!(g.id(), Some(13));
        assert_eq!(g.write_id(), 0);
    }
    #[test]
    fn every_flag_other_than_four_keeps_nonunspecified_atom_without_tetrahedral_restriction() {
        for flag in [
            None,
            Some(-1),
            Some(0),
            Some(1),
            Some(2),
            Some(3),
            Some(4),
            Some(5),
        ] {
            let mut r = topology(&[], ChiralTag::Unspecified);
            r.stereo_groups = vec![group(StereoGroupKind::And, &[0], &[], None, 0)];
            let mut p = product(&[], ChiralTag::Allene);
            if let Some(flag) = flag {
                p.topology.atoms[0]
                    .set_prop("molInversionFlag", PropertyValue::Int(flag))
                    .unwrap();
            }
            run(&r, &mut p, &mapping()).unwrap();
            assert_eq!(p.topology.stereo_groups.len(), usize::from(flag != Some(4)));
        }
    }
    #[test]
    fn absent_and_empty_mapping_or_bond_only_source_groups_leave_existing_groups_untouched() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![
            group(StereoGroupKind::And, &[1], &[], None, 0),
            group(StereoGroupKind::Or, &[2], &[], None, 0),
            group(StereoGroupKind::Absolute, &[], &[1], None, 0),
        ];
        let mut p = product(&[], ChiralTag::Other);
        p.topology.stereo_groups = vec![
            group(StereoGroupKind::Absolute, &[7], &[2], None, 3),
            group(StereoGroupKind::Absolute, &[8], &[1], None, 4),
        ];
        let before = p.topology.stereo_groups.clone();
        let mut m = mapping();
        m.reactant_to_product.insert(1, vec![]);
        run(&r, &mut p, &m).unwrap();
        assert_eq!(p.topology.stereo_groups, before);
    }
    #[test]
    fn new_nonabsolute_groups_prepend_and_existing_groups_keep_full_read_write_and_bond_state() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![
            group(StereoGroupKind::Or, &[0], &[2], Some(7), 99),
            group(StereoGroupKind::And, &[0], &[], Some(8), 88),
        ];
        let mut p = product(&[], ChiralTag::SquarePlanar);
        p.topology.stereo_groups =
            vec![group(StereoGroupKind::And, &[5, 5], &[4, 2], Some(91), 73)];
        let existing = p.topology.stereo_groups[0].clone();
        run(&r, &mut p, &mapping()).unwrap();
        assert_eq!(
            p.topology
                .stereo_groups
                .iter()
                .map(StereoGroup::kind)
                .collect::<Vec<_>>(),
            [
                StereoGroupKind::Or,
                StereoGroupKind::And,
                StereoGroupKind::And
            ]
        );
        assert_eq!(p.topology.stereo_groups[2], existing);
        assert_eq!(p.topology.stereo_groups[0].id(), Some(7));
        assert_eq!(p.topology.stereo_groups[0].write_id(), 0);
        assert!(p.topology.stereo_groups[0].bonds().is_empty());
    }
    #[test]
    fn absolute_groups_use_sole_source_reverse_concatenation_without_sorting_or_deduplication() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![group(StereoGroupKind::Absolute, &[0, 1], &[], Some(7), 9)];
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms[1].set_chiral_tag(ChiralTag::Other);
        p.topology.stereo_groups = vec![
            group(StereoGroupKind::Or, &[8], &[2], Some(12), 33),
            group(StereoGroupKind::Absolute, &[7, 7, 6], &[2, 1], Some(9), 22),
        ];
        let mut m = mapping();
        m.reactant_to_product.insert(1, vec![1]);
        run(&r, &mut p, &m).unwrap();
        assert_eq!(p.topology.stereo_groups.len(), 2);
        let g = &p.topology.stereo_groups[1];
        assert_eq!(
            g.atoms(),
            [
                AtomId::new(7),
                AtomId::new(7),
                AtomId::new(6),
                AtomId::new(0),
                AtomId::new(1)
            ]
        );
        assert_eq!(g.bonds(), [BondId::new(2), BondId::new(1)]);
        assert_eq!(g.id(), None);
        assert_eq!(g.write_id(), 0);
    }
    #[test]
    fn flag_conversion_failure_after_staged_groups_preserves_entire_existing_group_set() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![
            group(StereoGroupKind::Or, &[0], &[], None, 0),
            group(StereoGroupKind::And, &[1], &[], None, 0),
        ];
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms[1].set_chiral_tag(ChiralTag::Other);
        p.topology.atoms[1]
            .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
            .unwrap();
        p.topology.stereo_groups = vec![group(StereoGroupKind::And, &[9], &[3], Some(14), 19)];
        let before = p.topology.stereo_groups.clone();
        let mut m = mapping();
        m.reactant_to_product.insert(1, vec![1]);
        assert!(matches!(
            run(&r, &mut p, &m),
            Err(ReactionProductError::PropertyInt { .. })
        ));
        assert_eq!(p.topology.stereo_groups, before);
    }
    #[test]
    fn read_id_none_explicit_zero_and_nonzero_preserve_existing_model_projection_and_new_write_id_zero()
     {
        for id in [None, Some(0), Some(42)] {
            let mut r = topology(&[], ChiralTag::Unspecified);
            r.stereo_groups = vec![group(StereoGroupKind::And, &[0], &[], id, 88)];
            let mut p = product(&[], ChiralTag::Other);
            run(&r, &mut p, &mapping()).unwrap();
            assert_eq!(p.topology.stereo_groups[0].id(), id);
            assert_eq!(p.topology.stereo_groups[0].write_id(), 0);
        }
    }
    #[test]
    fn mapping_key_and_encoded_product_pointer_use_actual_atom_index_fields() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.atoms[1] = r.atoms[1].clone().with_id(AtomId::new(5));
        r.stereo_groups = vec![group(StereoGroupKind::And, &[1], &[], None, 0)];
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms[6].set_chiral_tag(ChiralTag::Other);
        p.topology.atoms[6] = p.topology.atoms[6].clone().with_id(AtomId::new(9));
        let mut m = mapping();
        m.reactant_to_product.insert(5, vec![6]);
        run(&r, &mut p, &m).unwrap();
        assert_eq!(p.topology.stereo_groups[0].atoms(), [AtomId::new(9)]);
    }
    #[test]
    fn reached_missing_source_or_product_group_atom_is_structural_error_before_assignment() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![group(StereoGroupKind::And, &[99], &[], None, 0)];
        let mut p = product(&[], ChiralTag::Other);
        assert!(matches!(
            run(&r, &mut p, &mapping()),
            Err(ReactionProductError::Invariant {
                detail: "reactant stereo-group atom row missing",
                ..
            })
        ));
        r.stereo_groups = vec![group(StereoGroupKind::And, &[0], &[], None, 0)];
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![99]);
        assert!(matches!(
            run(&r, &mut p, &m),
            Err(ReactionProductError::Invariant {
                detail: "product stereo-group atom row missing",
                ..
            })
        ));
        assert!(p.topology.stereo_groups.is_empty());
    }
    #[test]
    fn source_signed_flag_conversion_is_not_lossy_unsigned_cast_or_presence_only_filter() {
        let mut r = topology(&[], ChiralTag::Unspecified);
        r.stereo_groups = vec![group(StereoGroupKind::And, &[0], &[], None, 0)];
        let mut p = product(&[], ChiralTag::Other);
        p.topology.atoms[0]
            .set_prop("molInversionFlag", PropertyValue::UInt(u32::MAX))
            .unwrap();
        assert!(matches!(
            run(&r, &mut p, &mapping()),
            Err(ReactionProductError::PropertyInt { .. })
        ));
        assert!(p.topology.stereo_groups.is_empty());
    }
}

#[cfg(test)]
mod complete_propagate_coordinates_source_tests {
    use super::complete_correct_matched_chirality_source_tests::mapping;
    use super::*;
    use cosmolkit_model::{
        Atom, AtomSpec, Conformer2D, Conformer3D, CoordinateBlock, CoordinateDimension,
        CoordinateValidationError, MoleculeProperties, TopologyBlock,
    };
    use cosmolkit_types::Element;
    fn coordinates(points: Vec<[f64; 3]>, is_3d: bool) -> CoordinateBlock {
        CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(9, points, is_3d)],
            ..CoordinateBlock::default()
        }
    }
    fn run(
        c: &CoordinateBlock,
        atom_count: usize,
        points: &mut Vec<[f64; 3]>,
        is_3d: &mut bool,
        m: &ReactantProductMapping,
        selection: ReactionCoordinateSelection,
    ) -> Result<(), ReactionProductError> {
        let t = TopologyBlock {
            atoms: (0..atom_count)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            ..TopologyBlock::default()
        };
        propagate_coordinates(
            points,
            is_3d,
            &ReactionInput {
                topology: &t,
                coordinates: c,
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            m,
            selection,
        )
    }
    #[test]
    fn no_source_conformers_return_before_explicit_selection_mapping_flags_or_resize() {
        let mut points = vec![[7.0; 3]];
        let mut flag = false;
        let mut m = mapping();
        m.reactant_to_product.insert(99, vec![u32::MAX as usize]);
        run(
            &CoordinateBlock::default(),
            0,
            &mut points,
            &mut flag,
            &m,
            ReactionCoordinateSelection::ThreeD { id: 123 },
        )
        .unwrap();
        assert_eq!(points, [[7.0; 3]]);
        assert!(!flag);
    }
    #[test]
    fn ascending_source_keys_and_product_copy_order_overwrite_and_zero_extend_exactly() {
        let c = coordinates(vec![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]], true);
        let mut m = mapping();
        m.reactant_to_product =
            std::collections::BTreeMap::from([(1, vec![0, 3]), (0, vec![3, 1, 3])]);
        let mut points = vec![[9.0; 3]];
        let mut flag = false;
        run(
            &c,
            2,
            &mut points,
            &mut flag,
            &m,
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert_eq!(
            points,
            [[4.0, 5.0, 6.0], [1.0, 2.0, 3.0], [0.0; 3], [4.0, 5.0, 6.0]]
        );
        assert!(flag);
    }
    #[test]
    fn source_three_d_flag_is_written_before_first_reached_position_error() {
        let c = coordinates(vec![[1.0; 3]], true);
        let mut m = mapping();
        m.reactant_to_product = std::collections::BTreeMap::from([(99, vec![0])]);
        let mut points = vec![[9.0; 3]];
        let mut flag = false;
        assert!(matches!(
            run(
                &c,
                1,
                &mut points,
                &mut flag,
                &m,
                ReactionCoordinateSelection::Auto
            ),
            Err(ReactionProductError::Invariant {
                detail: "reactant conformer atom out of range",
                ..
            })
        ));
        assert!(flag);
        assert_eq!(points, [[9.0; 3]]);
    }
    #[test]
    fn non_three_d_sources_never_clear_existing_three_d_flag_and_two_d_adds_zero_z() {
        let c = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[-0.0, 2.0]])],
            ..CoordinateBlock::default()
        };
        let mut points = vec![];
        let mut flag = true;
        run(
            &c,
            1,
            &mut points,
            &mut flag,
            &mapping(),
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert!(flag);
        assert_eq!(points[0][0].to_bits(), (-0.0f64).to_bits());
        assert_eq!(points[0][2].to_bits(), 0.0f64.to_bits());
        let c = coordinates(vec![[1.0, 2.0, 7.0]], false);
        let mut flag = false;
        run(
            &c,
            1,
            &mut points,
            &mut flag,
            &mapping(),
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert!(!flag);
        assert_eq!(points[0], [1.0, 2.0, 7.0]);
    }
    #[test]
    fn empty_product_copy_vector_never_reads_source_index_or_owning_shape() {
        let c = coordinates(vec![[1.0; 3]], true);
        let mut m = mapping();
        m.reactant_to_product = std::collections::BTreeMap::from([(99, vec![])]);
        let mut points = vec![[9.0; 3]];
        let mut flag = false;
        run(
            &c,
            2,
            &mut points,
            &mut flag,
            &m,
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert!(flag);
        assert_eq!(points, [[9.0; 3]]);
    }
    #[test]
    fn owning_conformer_shape_precondition_precedes_source_index_and_destination_overflow() {
        let c = coordinates(vec![[1.0; 3]], true);
        let mut m = mapping();
        m.reactant_to_product = std::collections::BTreeMap::from([(99, vec![u32::MAX as usize])]);
        let mut points = vec![];
        let mut flag = false;
        assert!(matches!(
            run(
                &c,
                2,
                &mut points,
                &mut flag,
                &m,
                ReactionCoordinateSelection::Auto
            ),
            Err(ReactionProductError::Coordinate(
                CoordinateValidationError::RowCount { .. }
            ))
        ));
        assert!(flag);
        assert!(points.is_empty());
    }
    #[test]
    fn source_getter_error_precedes_destination_setter_overflow() {
        let c = coordinates(vec![[1.0; 3]], false);
        let mut m = mapping();
        m.reactant_to_product = std::collections::BTreeMap::from([(99, vec![u32::MAX as usize])]);
        assert!(matches!(
            run(
                &c,
                1,
                &mut vec![],
                &mut false,
                &m,
                ReactionCoordinateSelection::Auto
            ),
            Err(ReactionProductError::Invariant {
                detail: "reactant conformer atom out of range",
                ..
            })
        ));
    }
    #[test]
    fn destination_unsigned_max_error_preserves_earlier_copy_and_flag_prefix_without_allocation() {
        let c = coordinates(vec![[1.0, 2.0, 3.0]], true);
        let mut m = mapping();
        m.reactant_to_product.insert(0, vec![0, u32::MAX as usize]);
        let mut points = vec![[9.0; 3]];
        let mut flag = false;
        assert!(
            matches!(run(&c,1,&mut points,&mut flag,&m,ReactionCoordinateSelection::Auto),Err(ReactionProductError::Coordinate(CoordinateValidationError::AtomIndexOverflow {atom})) if atom==u32::MAX as usize)
        );
        assert_eq!(points, [[1.0, 2.0, 3.0]]);
        assert!(flag);
    }
    #[test]
    fn source_position_copy_preserves_nan_payload_infinity_and_signed_zero_bits() {
        let point = [
            f64::from_bits(0x7ff8_0000_0000_1234),
            f64::NEG_INFINITY,
            -0.0,
        ];
        let c = coordinates(vec![point], true);
        let mut points = vec![];
        run(
            &c,
            1,
            &mut points,
            &mut false,
            &mapping(),
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert_eq!(points[0].map(f64::to_bits), point.map(f64::to_bits));
    }
    #[test]
    fn auto_uses_actual_source_front_not_minimum_conformer_identifier() {
        let c = CoordinateBlock {
            conformers_3d: vec![
                Conformer3D::new(9, vec![[9.0; 3]], true),
                Conformer3D::new(1, vec![[1.0; 3]], true),
            ],
            ..CoordinateBlock::default()
        };
        let mut points = vec![];
        run(
            &c,
            1,
            &mut points,
            &mut false,
            &mapping(),
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert_eq!(points, [[9.0; 3]]);
    }
    #[test]
    fn mixed_source_order_and_explicit_dimension_selection_reuse_sole_conformer_selector() {
        let c = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(7, vec![[7.0, 2.0]])],
            conformers_3d: vec![Conformer3D::new(9, vec![[9.0; 3]], true)],
            source_conformer_order: Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]),
            ..CoordinateBlock::default()
        };
        let mut points = vec![];
        run(
            &c,
            1,
            &mut points,
            &mut false,
            &mapping(),
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert_eq!(points, [[9.0; 3]]);
        run(
            &c,
            1,
            &mut points,
            &mut false,
            &mapping(),
            ReactionCoordinateSelection::TwoD { id: 7 },
        )
        .unwrap();
        assert_eq!(points, [[7.0, 2.0, 0.0]]);
    }
    #[test]
    fn explicit_missing_selection_errors_before_flag_or_coordinate_writes_on_nonempty_input() {
        let c = coordinates(vec![[1.0; 3]], true);
        let mut points = vec![[9.0; 3]];
        let mut flag = false;
        assert!(matches!(
            run(
                &c,
                1,
                &mut points,
                &mut flag,
                &mapping(),
                ReactionCoordinateSelection::ThreeD { id: 77 }
            ),
            Err(ReactionProductError::CoordinateSelection(_))
        ));
        assert!(!flag);
        assert_eq!(points, [[9.0; 3]]);
    }
}

#[cfg(test)]
mod complete_copy_template_groups_source_tests {
    use super::complete_correct_matched_chirality_source_tests::product;
    use super::*;
    use cosmolkit_model::{AtomSpec, Element, PropertyValue, QueryAtom, StereoGroupKind};
    fn group(
        kind: StereoGroupKind,
        atoms: &[usize],
        bonds: &[usize],
        id: u32,
        write: u32,
    ) -> StereoGroup {
        StereoGroup::new(
            kind,
            atoms.iter().map(|&i| AtomId::new(i)).collect(),
            bonds.iter().map(|&i| BondId::new(i)).collect(),
        )
        .with_id(id)
        .with_write_id(write)
    }
    fn template(maps: &[Option<i32>], groups: Vec<StereoGroup>) -> QueryGraph {
        let atoms = maps
            .iter()
            .enumerate()
            .map(|(i, map)| {
                let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
                if let Some(map) = map {
                    a.set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                        .unwrap();
                }
                a
            })
            .collect();
        let mut q = QueryGraph::from_parts(atoms, vec![], [], vec![], vec![], vec![]).unwrap();
        for group in groups {
            q.add_stereo_group(group);
        }
        q
    }
    fn mapped_product(maps: &[Option<i32>]) -> ProductBuilder {
        let mut p = product(&[], ChiralTag::Other);
        for (atom, map) in p.topology.atoms.iter_mut().zip(maps) {
            if let Some(map) = map {
                atom.set_prop("old_mapno", PropertyValue::Int(*map))
                    .unwrap();
            }
        }
        p
    }
    #[test]
    fn empty_template_groups_skip_all_product_rows_properties_and_existing_members() {
        let q = template(&[], vec![]);
        let mut p = mapped_product(&[]);
        p.topology.stereo_groups = vec![group(StereoGroupKind::And, &[99], &[98], 7, 9)];
        let before = p.topology.clone();
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(p.topology, before);
    }
    #[test]
    fn template_member_order_expands_all_matching_product_atoms_without_filtering_or_deduplication()
    {
        let q = template(
            &[Some(7), Some(8)],
            vec![group(StereoGroupKind::Or, &[1, 0, 1], &[4], 13, 77)],
        );
        let mut p = mapped_product(&[Some(7), Some(8), Some(7)]);
        p.topology.atoms[0].set_chiral_tag(ChiralTag::Unspecified);
        p.topology.atoms[1]
            .set_prop("molInversionFlag", PropertyValue::String("ignored".into()))
            .unwrap();
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(
            p.topology.stereo_groups[0].atoms(),
            [
                AtomId::new(1),
                AtomId::new(0),
                AtomId::new(2),
                AtomId::new(1)
            ]
        );
        assert!(p.topology.stereo_groups[0].bonds().is_empty());
        assert_eq!(p.topology.stereo_groups[0].id(), Some(13));
        assert_eq!(p.topology.stereo_groups[0].write_id(), 0);
    }
    #[test]
    fn zero_map_abandons_group_but_earlier_marks_still_remove_overlapping_existing_group() {
        let q = template(
            &[Some(7), None, Some(8)],
            vec![
                group(StereoGroupKind::Or, &[0, 1], &[], 1, 0),
                group(StereoGroupKind::And, &[2], &[], 2, 0),
            ],
        );
        let mut p = mapped_product(&[Some(7), Some(8)]);
        p.topology.stereo_groups = vec![group(StereoGroupKind::And, &[0], &[2], 99, 44)];
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(p.topology.stereo_groups.len(), 1);
        assert_eq!(p.topology.stereo_groups[0].atoms(), [AtomId::new(1)]);
        assert_eq!(p.topology.stereo_groups[0].id(), Some(2));
    }
    #[test]
    fn zero_or_absent_first_map_skips_product_property_casts_and_keeps_existing_groups_when_none_new()
     {
        for map in [None, Some(0)] {
            let q = template(&[map], vec![group(StereoGroupKind::And, &[0], &[], 1, 0)]);
            let mut p = mapped_product(&[]);
            p.topology.atoms[0]
                .set_prop("old_mapno", PropertyValue::String("unreached".into()))
                .unwrap();
            p.topology.stereo_groups = vec![
                group(StereoGroupKind::Absolute, &[8], &[], 1, 3),
                group(StereoGroupKind::Absolute, &[9], &[], 2, 4),
            ];
            let before = p.topology.clone();
            copy_template_groups(&mut p, &q, 0).unwrap();
            assert_eq!(p.topology, before);
        }
    }
    #[test]
    fn no_product_matches_or_empty_template_members_do_not_trigger_existing_group_merge() {
        let q = template(
            &[Some(7)],
            vec![
                group(StereoGroupKind::And, &[0], &[], 1, 0),
                group(StereoGroupKind::Or, &[], &[2], 2, 0),
            ],
        );
        let mut p = mapped_product(&[Some(8)]);
        p.topology.stereo_groups = vec![
            group(StereoGroupKind::Absolute, &[8], &[], 1, 3),
            group(StereoGroupKind::Absolute, &[9], &[], 2, 4),
        ];
        let before = p.topology.stereo_groups.clone();
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(p.topology.stereo_groups, before);
    }
    #[test]
    fn existing_no_overlap_is_copied_partial_overlap_is_split_and_full_overlap_is_removed() {
        let q = template(
            &[Some(7)],
            vec![group(StereoGroupKind::Or, &[0], &[], 13, 77)],
        );
        let mut p = mapped_product(&[Some(7)]);
        let no_overlap = group(StereoGroupKind::And, &[3, 3], &[2, 1], 21, 88);
        p.topology.stereo_groups = vec![
            no_overlap.clone(),
            group(StereoGroupKind::And, &[0, 0, 2, 2], &[4, 3], 22, 99),
            group(StereoGroupKind::And, &[0, 0], &[5], 23, 66),
            group(StereoGroupKind::Or, &[], &[6], 24, 55),
        ];
        let bond_only = p.topology.stereo_groups[3].clone();
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(p.topology.stereo_groups.len(), 4);
        assert_eq!(p.topology.stereo_groups[1], no_overlap);
        assert_eq!(
            p.topology.stereo_groups[2].atoms(),
            [AtomId::new(2), AtomId::new(2)]
        );
        assert!(p.topology.stereo_groups[2].bonds().is_empty());
        assert_eq!(p.topology.stereo_groups[2].id(), Some(22));
        assert_eq!(p.topology.stereo_groups[2].write_id(), 0);
        assert_eq!(p.topology.stereo_groups[3], bond_only);
    }
    #[test]
    fn negative_nonzero_map_numbers_are_compared_as_signed_source_ints() {
        let q = template(
            &[Some(-7)],
            vec![group(StereoGroupKind::And, &[0], &[], 1, 0)],
        );
        let mut p = mapped_product(&[Some(-7), Some(7), Some(-7)]);
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(
            p.topology.stereo_groups[0].atoms(),
            [AtomId::new(0), AtomId::new(2)]
        );
    }
    #[test]
    fn reached_bad_product_property_after_staging_keeps_original_groups_without_partial_commit() {
        let q = template(
            &[Some(7), Some(8)],
            vec![
                group(StereoGroupKind::And, &[0], &[], 1, 0),
                group(StereoGroupKind::Or, &[1], &[], 2, 0),
            ],
        );
        let mut p = mapped_product(&[Some(7)]);
        p.topology.atoms[1]
            .set_prop("old_mapno", PropertyValue::String("bad".into()))
            .unwrap();
        p.topology.stereo_groups = vec![group(StereoGroupKind::And, &[4], &[2], 9, 33)];
        let before = p.topology.stereo_groups.clone();
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::PropertyInt { .. })
        ));
        assert_eq!(p.topology.stereo_groups, before);
    }
    #[test]
    fn template_member_and_matched_product_bit_range_failures_are_structural_before_assignment() {
        let q = template(
            &[Some(7)],
            vec![group(StereoGroupKind::And, &[99], &[], 1, 0)],
        );
        let mut p = mapped_product(&[Some(7)]);
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::Invariant {
                detail: "template stereo atom missing",
                ..
            })
        ));
        let q = template(
            &[Some(7)],
            vec![group(StereoGroupKind::And, &[0], &[], 1, 0)],
        );
        p.topology.atoms[0] = p.topology.atoms[0].clone().with_id(AtomId::new(99));
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::Invariant {
                detail: "stereo-group atom bit index out of range",
                ..
            })
        ));
        assert!(p.topology.stereo_groups.is_empty());
    }
    #[test]
    fn existing_group_read_error_after_new_group_staging_preserves_original_group_set() {
        let q = template(
            &[Some(7)],
            vec![group(StereoGroupKind::And, &[0], &[], 1, 0)],
        );
        let mut p = mapped_product(&[Some(7)]);
        p.topology.stereo_groups = vec![
            group(StereoGroupKind::And, &[3], &[], 1, 2),
            group(StereoGroupKind::Or, &[99], &[2], 9, 33),
        ];
        let before = p.topology.stereo_groups.clone();
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::Invariant {
                detail: "existing stereo-group atom row missing",
                ..
            })
        ));
        assert_eq!(p.topology.stereo_groups, before);
        p.topology.stereo_groups = vec![group(StereoGroupKind::And, &[3], &[], 1, 2)];
        p.topology.atoms[3] = p.topology.atoms[3].clone().with_id(AtomId::new(99));
        let before = p.topology.stereo_groups.clone();
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::Invariant {
                detail: "stereo-group atom bit index out of range",
                ..
            })
        ));
        assert_eq!(p.topology.stereo_groups, before);
    }
    #[test]
    fn final_absolute_merge_retains_source_reverse_concatenation_and_duplicate_members() {
        let q = template(
            &[Some(7), Some(8)],
            vec![
                group(StereoGroupKind::Absolute, &[0], &[], 1, 7),
                group(StereoGroupKind::Absolute, &[1], &[], 2, 8),
            ],
        );
        let mut p = mapped_product(&[Some(7), Some(8)]);
        p.topology.stereo_groups =
            vec![group(StereoGroupKind::Absolute, &[8, 8, 7], &[2, 1], 3, 9)];
        copy_template_groups(&mut p, &q, 0).unwrap();
        assert_eq!(p.topology.stereo_groups.len(), 1);
        let g = &p.topology.stereo_groups[0];
        assert_eq!(
            g.atoms(),
            [
                AtomId::new(8),
                AtomId::new(8),
                AtomId::new(7),
                AtomId::new(1),
                AtomId::new(0)
            ]
        );
        assert_eq!(g.bonds(), [BondId::new(2), BondId::new(1)]);
        assert_eq!(g.id(), None);
        assert_eq!(g.write_id(), 0);
    }
    #[test]
    fn out_of_signed_range_template_map_is_not_wrapped_or_treated_as_absent() {
        let mut q = template(&[None], vec![group(StereoGroupKind::And, &[0], &[], 1, 0)]);
        q.atom_mut(0)
            .unwrap()
            .set_prop("molAtomMapNumber", PropertyValue::UInt(u32::MAX))
            .unwrap();
        let mut p = mapped_product(&[]);
        assert!(matches!(
            copy_template_groups(&mut p, &q, 0),
            Err(ReactionProductError::TemplateProperty(_))
        ));
        assert!(p.topology.stereo_groups.is_empty());
    }
}
