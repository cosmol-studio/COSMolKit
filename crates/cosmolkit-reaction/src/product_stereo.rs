use crate::materialize::{
    ProductBuilder, ReactantProductMapping, int_prop, invariant, inversion_flag, uint_prop,
};
use crate::{ReactionCoordinateSelection, ReactionInput, ReactionProductError, ReactionRole};
use cosmolkit_core::{
    StereoAtomSearchWarning, count_swaps_to_interconvert,
    find_double_bond_stereo_atoms_with_rank_reader, invert_atom_chirality,
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
        // Borrowed canonical context is constant-cost; only the reached getter
        // checks the current valence fact. Neighbor search is O(degree), as source.
        let target = SearchTarget::new(
            input.topology,
            input.coordinates,
            &input.topology.stereo_groups,
            input.rings,
            input.valence,
        );
        let context = build_query_match_context(&target);
        let total_degree = atom_total_degree_with_context(&input.topology.atoms[atom], &context)
            .map_err(|source| ReactionProductError::StereoGetter {
                atom: AtomId::new(atom),
                source,
            })?;
        if total_degree > 3 {
            return Err(invariant(
                "StereoBondEndCap",
                "Stereo Bond extremes must have less than four neighbors",
                Some(atom),
                None,
                None,
            ));
        }
        let non_anchor = input
            .topology
            .adjacency
            .neighbors_of(atom)
            .iter()
            .find(|n| n.atom_index != other && n.atom_index != anchor)
            .map(|n| n.atom_index);
        Ok(Self { anchor, non_anchor })
    }

    fn product_candidates<'a>(&self, mapping: &'a ReactantProductMapping) -> (&'a [usize], bool) {
        // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: getProductAnchorCandidates
        // RDKit❗✔️:   std::pair<UINT_VECT, bool> getProductAnchorCandidates(
        // RDKit❗✔️:       ReactantProductAtomMapping *mapping) {
        // RDKit❗✔️:     auto &react2Prod = mapping->reactProdAtomMap;
        // RDKit❗✔️:
        // RDKit❗✔️:     bool swapStereo = false;
        // RDKit❗✔️:     auto newAnchorMatches = react2Prod.find(getAnchorIdx());
        // RDKit❗✔️:     if (newAnchorMatches != react2Prod.end()) {
        // RDKit❗✔️:       // The corresponding StereoAtom exists in the product
        // RDKit❗✔️:       return {newAnchorMatches->second, swapStereo};
        // RDKit❗✔️:
        // RDKit❗✔️:     } else if (hasNonAnchor()) {
        // RDKit❗✔️:       // The non-StereoAtom neighbor exists in the product
        // RDKit❗✔️:       newAnchorMatches = react2Prod.find(getNonAnchorIdx());
        // RDKit❗✔️:       if (newAnchorMatches != react2Prod.end()) {
        // RDKit❗✔️:         swapStereo = true;
        // RDKit❗✔️:         return {newAnchorMatches->second, swapStereo};
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     // None of the neighbors survived the reaction
        // RDKit❗✔️:     return {{}, swapStereo};
        // RDKit❗✔️:   }
        // END RDKIT COMPLETE CPP FUNCTION
        // Source copies UINT_VECT; a borrowed slice preserves its order and empty
        // present-key distinction while avoiding that allocation.
        if let Some(matches) = mapping.reactant_to_product.get(&self.anchor) {
            return (matches, false);
        }
        if let Some(other) = self.non_anchor {
            if let Some(matches) = mapping.reactant_to_product.get(&other) {
                return (matches, true);
            }
        }
        (&[], false)
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
    candidates
        .iter()
        .copied()
        .find(|row| product.bond_between(atom, *row).is_some())
        .ok_or_else(|| {
            invariant(
                "reactProdMapAnchorIdx",
                "match not found",
                None,
                Some(atom),
                None,
            )
        })
}

fn forward_bond_stereo(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
    product_bond: BondId,
    reactant_bond: BondId,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: forwardReactantBondStereo
    // RDKit❗✔️: void forwardReactantBondStereo(ReactantProductAtomMapping *mapping, Bond *pBond,
    // RDKit❗✔️:                                const ROMol &reactant, const Bond *rBond) {
    // RDKit❗✔️:   PRECONDITION(mapping, "no mapping");
    // RDKit❗✔️:   PRECONDITION(pBond, "no bond");
    // RDKit❗✔️:   PRECONDITION(rBond, "no bond");
    // RDKit❗✔️:   PRECONDITION(rBond->getStereo() > Bond::BondStereo::STEREOANY,
    // RDKit❗✔️:                "bond in reactant must have defined stereo");
    // RDKit❗✔️:
    // RDKit❗✔️:   auto &prod2React = mapping->prodReactAtomMap;
    // RDKit❗✔️:
    // RDKit❗✔️:   const Atom *rStart = rBond->getBeginAtom();
    // RDKit❗✔️:   const Atom *rEnd = rBond->getEndAtom();
    // RDKit❗✔️:   const auto rStereoAtoms = Chirality::findStereoAtoms(rBond);
    // RDKit❗✔️:   if (rStereoAtoms.size() != 2) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:         << "WARNING: neither stereo atoms nor CIP codes found for double bond. "
    // RDKit❗✔️:            "Stereochemistry info will not be propagated to product."
    // RDKit❗✔️:         << std::endl;
    // RDKit❗✔️:     pBond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   StereoBondEndCap start(reactant, rStart, rEnd, rStereoAtoms[0]);
    // RDKit❗✔️:   StereoBondEndCap end(reactant, rEnd, rStart, rStereoAtoms[1]);
    // RDKit❗✔️:
    // RDKit❗✔️:   // The bond might be matched backwards in the reaction
    // RDKit❗✔️:   if (prod2React[pBond->getBeginAtomIdx()] == rEnd->getIdx()) {
    // RDKit❗✔️:     std::swap(start, end);
    // RDKit❗✔️:   } else if (prod2React[pBond->getBeginAtomIdx()] != rStart->getIdx()) {
    // RDKit❗✔️:     throw std::logic_error("Reactant and Product bond ends do not match");
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   /**
    // RDKit❗✔️:    *  The reactants stereo can be transmitted in three similar ways:
    // RDKit❗✔️:    *
    // RDKit❗✔️:    * 1. Survival of both stereoatoms: direct forwarding happens, i.e.,
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    C/C=C/[Br] in reaction [C:1]=[C:2]>>[Si:1]=[C:2]:
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    C/C=C/[Br] >> C/Si=C/[Br], C/C=Si/[Br] (2 product sets)
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    Both stereoatoms exist unaltered in both product sets, so we can forward
    // RDKit❗✔️:    *    the same bond stereochemistry (trans) and set the stereoatoms in the
    // RDKit❗✔️:    *    product to the mapped indexes of the stereoatoms in the reactant.
    // RDKit❗✔️:    *
    // RDKit❗✔️:    * 2. Survival of both anti-stereoatoms: as this pair is symmetric to the
    // RDKit❗✔️:    *    stereoatoms, direct forwarding also happens in this case, i.e.,
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    Cl/C(C)=C(/Br)F in reaction
    // RDKit❗✔️:    *        [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3].[Cl:4]:
    // RDKit❗✔️:    *      Cl/C(C)=C(/Br)F >> C/C=C/F + Br + Cl
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    Both stereoatoms in the reactant are split from the molecule,
    // RDKit❗✔️:    *    but the anti-stereoatoms remain in it. Since these have symmetrical
    // RDKit❗✔️:    *    orientation to the stereoatoms, we can use these (their mapped
    // RDKit❗✔️:    *    equivalents) as stereoatoms in the product and use the same
    // RDKit❗✔️:    *    stereochemistry label (trans).
    // RDKit❗✔️:    *
    // RDKit❗✔️:    * 3. Survival of a mixed pair stereoatom-anti-stereoatom: such a pair
    // RDKit❗✔️:    *    defines the opposite stereochemistry to the one labeled on the
    // RDKit❗✔️:    *    reactant, but it is also valid, as long ase we use the properly mapped
    // RDKit❗✔️:    *    indexes:
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    Cl/C(C)=C(/Br)F in reaction [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Br:3]:
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *        Cl/C(C)=C(/Br)F >> C/C=C/F + Br
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    In this case, one of the stereoatoms is conserved, and the other one is
    // RDKit❗✔️:    *    switched to the other neighbor at the same end of the bond as the
    // RDKit❗✔️:    *    non-conserved stereoatom. Since the reference changed, the
    // RDKit❗✔️:    *    stereochemistry label needs to be flipped too: in this case, the
    // RDKit❗✔️:    *    reactant was trans, and the product will be cis.
    // RDKit❗✔️:    *
    // RDKit❗✔️:    *    Reaction [Cl:4][C:1]=[C:2][Br:3]>>[C:1]=[C:2].[Cl:4] would have the same
    // RDKit❗✔️:    *    effect, with the only difference that the non-conserved stereoatom would
    // RDKit❗✔️:    *    be the one at the opposite end of the reactant.
    // RDKit❗✔️:    */
    // RDKit❗✔️:   auto pStartAnchorCandidates = start.getProductAnchorCandidates(mapping);
    // RDKit❗✔️:   auto pEndAnchorCandidates = end.getProductAnchorCandidates(mapping);
    // RDKit❗✔️:
    // RDKit❗✔️:   // The reaction has invalidated the reactant's stereochemistry
    // RDKit❗✔️:   if (pStartAnchorCandidates.first.empty() ||
    // RDKit❗✔️:       pEndAnchorCandidates.first.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   unsigned pStartAnchorIdx = reactProdMapAnchorIdx(
    // RDKit❗✔️:       pBond->getBeginAtom(), pStartAnchorCandidates.first);
    // RDKit❗✔️:   unsigned pEndAnchorIdx =
    // RDKit❗✔️:       reactProdMapAnchorIdx(pBond->getEndAtom(), pEndAnchorCandidates.first);
    // RDKit❗✔️:
    // RDKit❗✔️:   const ROMol &m = pBond->getOwningMol();
    // RDKit❗✔️:   if (m.getBondBetweenAtoms(pBond->getBeginAtomIdx(), pStartAnchorIdx) ==
    // RDKit❗✔️:           nullptr ||
    // RDKit❗✔️:       m.getBondBetweenAtoms(pBond->getEndAtomIdx(), pEndAnchorIdx) == nullptr) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "stereo atoms in input cannot be mapped to "
    // RDKit❗✔️:                                "output (atoms are no longer bonded)\n";
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     pBond->setStereoAtoms(pStartAnchorIdx, pEndAnchorIdx);
    // RDKit❗✔️:     bool flipStereo =
    // RDKit❗✔️:         (pStartAnchorCandidates.second + pEndAnchorCandidates.second) % 2;
    // RDKit❗✔️:
    // RDKit❗✔️:     if (rBond->getStereo() == Bond::BondStereo::STEREOCIS ||
    // RDKit❗✔️:         rBond->getStereo() == Bond::BondStereo::STEREOZ) {
    // RDKit❗✔️:       if (flipStereo) {
    // RDKit❗✔️:         pBond->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         pBond->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       if (flipStereo) {
    // RDKit❗✔️:         pBond->setStereo(Bond::BondStereo::STEREOCIS);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         pBond->setStereo(Bond::BondStereo::STEREOTRANS);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // No rank assignment: the sole CORE reader first uses stored reference
    // atoms, then reaches unsigned property conversion only where source does.
    let found =
        find_double_bond_stereo_atoms_with_rank_reader(input.topology, reactant_bond, |atom| {
            uint_prop(&input.topology.atoms[atom.index()], "_CIPRank")
        })?;
    for warning in found.warnings {
        match warning {
            StereoAtomSearchWarning::DuplicateCipRank { center: _ } => {
                eprintln!("Warning: duplicate CIP ranks found in findHighestCIPNeighbor()")
            }
            StereoAtomSearchWarning::UnableToAssign { bond } => {
                eprintln!("Unable to assign stereo atoms for bond {}", bond.index())
            }
        }
    }
    let Some([a, b]) = found.atoms else {
        eprintln!(
            "WARNING: neither stereo atoms nor CIP codes found for double bond. Stereochemistry info will not be propagated to product."
        );
        product.topology.bonds[product_bond.index()].set_stereo(BondStereo::None)?;
        return Ok(());
    };
    let r = &input.topology.bonds[reactant_bond.index()];
    let mut start = StereoBondEndCap::new(input, r.begin().index(), r.end().index(), a.index())?;
    let mut end = StereoBondEndCap::new(input, r.end().index(), r.begin().index(), b.index())?;
    let p = &product.topology.bonds[product_bond.index()];
    let begin = p.begin().index();
    let finish = p.end().index();
    let source_begin = mapping
        .product_to_reactant
        .get(&begin)
        .copied()
        .ok_or_else(|| {
            invariant(
                "forwardReactantBondStereo",
                "Reactant and Product bond ends do not match",
                None,
                Some(begin),
                Some(product_bond),
            )
        })?;
    if source_begin == r.end().index() {
        std::mem::swap(&mut start, &mut end);
    } else if source_begin != r.begin().index() {
        return Err(invariant(
            "forwardReactantBondStereo",
            "Reactant and Product bond ends do not match",
            Some(source_begin),
            Some(begin),
            Some(product_bond),
        ));
    }
    let (start_matches, start_swap) = start.product_candidates(mapping);
    let (end_matches, end_swap) = end.product_candidates(mapping);
    if start_matches.is_empty() || end_matches.is_empty() {
        return Ok(());
    }
    let start_anchor = anchor_index(product, begin, start_matches)?;
    let end_anchor = anchor_index(product, finish, end_matches)?;
    if product.bond_between(begin, start_anchor).is_none()
        || product.bond_between(finish, end_anchor).is_none()
    {
        eprintln!("stereo atoms in input cannot be mapped to output (atoms are no longer bonded)");
        return Ok(());
    }
    let flip = start_swap != end_swap;
    let cis = matches!(r.stereo(), BondStereo::Cis | BondStereo::Z) != flip;
    let p = &mut product.topology.bonds[product_bond.index()];
    p.set_stereo_atoms(Some([AtomId::new(start_anchor), AtomId::new(end_anchor)]));
    p.set_stereo(if cis {
        BondStereo::Cis
    } else {
        BondStereo::Trans
    })?;
    Ok(())
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
    let p = &product.topology.bonds[bond.index()];
    let begin = p.begin();
    let finish = p.end();
    let start = &product.topology.bonds[start.index()];
    let end = &product.topology.bonds[end.index()];
    let start_anchor = if start.begin() == begin {
        start.end()
    } else {
        start.begin()
    };
    let end_anchor = if end.begin() == finish {
        end.end()
    } else {
        end.begin()
    };
    let mut same = start.direction() == end.direction();
    if start.begin() == begin {
        same = !same;
    }
    if end.begin() != finish {
        same = !same;
    }
    let p = &mut product.topology.bonds[bond.index()];
    p.set_stereo_atoms(Some([start_anchor, end_anchor]));
    p.set_stereo(if same {
        BondStereo::Trans
    } else {
        BondStereo::Cis
    })?;
    Ok(())
}

fn directed_neighbor(product: &ProductBuilder, atom: usize) -> Option<BondId> {
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
    product.neighbors[atom].iter().find_map(|neighbor| {
        let bond = &product.topology.bonds[neighbor.bond.index()];
        (bond.order() != BondOrder::Double
            && matches!(
                bond.direction(),
                BondDirection::EndDownRight | BondDirection::EndUpRight
            ))
        .then_some(neighbor.bond)
    })
}

pub(crate) fn update_stereo_bonds(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateStereoBonds
    // RDKit❗✔️: void updateStereoBonds(RWMOL_SPTR product, const ROMol &reactant,
    // RDKit❗✔️:                        ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   for (Bond *pBond : product->bonds()) {
    // RDKit❗✔️:     // We are only interested in double bonds
    // RDKit❗✔️:     if (pBond->getBondType() != Bond::BondType::DOUBLE) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     } else if (pBond->hasProp(_UnknownStereoRxnBond)) {
    // RDKit❗✔️:       pBond->setStereo(Bond::BondStereo::STEREONONE);
    // RDKit❗✔️:       pBond->clearProp(_UnknownStereoRxnBond);
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // Check if the reaction defined the stereo for the bond: SMARTS can only
    // RDKit❗✔️:     // use bond directions for this, and both sides of the double bond must have
    // RDKit❗✔️:     // them, else they will be ignored, as there is no reference to decide the
    // RDKit❗✔️:     // stereo.
    // RDKit❗✔️:     const auto *pBondStartDirBond =
    // RDKit❗✔️:         Chirality::getNeighboringDirectedBond(*product, pBond->getBeginAtom());
    // RDKit❗✔️:     const auto *pBondEndDirBond =
    // RDKit❗✔️:         Chirality::getNeighboringDirectedBond(*product, pBond->getEndAtom());
    // RDKit❗✔️:     if (pBondStartDirBond != nullptr && pBondEndDirBond != nullptr) {
    // RDKit❗✔️:       translateProductStereoBondDirections(pBond, pBondStartDirBond,
    // RDKit❗✔️:                                            pBondEndDirBond);
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // If the reaction did not specify the stereo, then we need to rely on the
    // RDKit❗✔️:       // atom mapping and use the reactant's stereo.
    // RDKit❗✔️:
    // RDKit❗✔️:       // The atoms and the bond might have been added in the reaction
    // RDKit❗✔️:       const auto begIdxItr =
    // RDKit❗✔️:           mapping->prodReactAtomMap.find(pBond->getBeginAtomIdx());
    // RDKit❗✔️:       if (begIdxItr == mapping->prodReactAtomMap.end()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       const auto endIdxItr =
    // RDKit❗✔️:           mapping->prodReactAtomMap.find(pBond->getEndAtomIdx());
    // RDKit❗✔️:       if (endIdxItr == mapping->prodReactAtomMap.end()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       const Bond *rBond =
    // RDKit❗✔️:           reactant.getBondBetweenAtoms(begIdxItr->second, endIdxItr->second);
    // RDKit❗✔️:
    // RDKit❗✔️:       if (rBond && rBond->getBondType() == Bond::BondType::DOUBLE) {
    // RDKit❗✔️:         // The bond might not have been present in the reactant, or its order
    // RDKit❗✔️:         // might have changed
    // RDKit❗✔️:         if (rBond->getStereo() > Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:           // If the bond had stereo, forward it
    // RDKit❗✔️:           forwardReactantBondStereo(mapping, pBond, reactant, rBond);
    // RDKit❗✔️:         } else if (rBond->getStereo() == Bond::BondStereo::STEREOANY) {
    // RDKit❗✔️:           pBond->setStereo(Bond::BondStereo::STEREOANY);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // No stereo: Bond::BondStereo::STEREONONE
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Source bond-order traversal and local adjacency lookups; no graph copies.
    for row in 0..product.topology.bonds.len() {
        let p = &product.topology.bonds[row];
        if p.order() != BondOrder::Double {
            continue;
        }
        if p.prop("_UnknownStereoRxnBond").is_some() {
            let p = &mut product.topology.bonds[row];
            p.set_stereo(BondStereo::None)?;
            p.clear_prop("_UnknownStereoRxnBond");
            continue;
        }
        let begin = p.begin().index();
        let end = p.end().index();
        let id = p.id();
        if let (Some(start), Some(finish)) = (
            directed_neighbor(product, begin),
            directed_neighbor(product, end),
        ) {
            translate_directions(product, id, start, finish)?;
        } else {
            let Some(&r_begin) = mapping.product_to_reactant.get(&begin) else {
                continue;
            };
            let Some(&r_end) = mapping.product_to_reactant.get(&end) else {
                continue;
            };
            let Some(r_bond) = input
                .topology
                .adjacency
                .neighbors_of(r_begin)
                .iter()
                .find(|n| n.atom_index == r_end)
                .map(|n| n.bond)
            else {
                continue;
            };
            let r = &input.topology.bonds[r_bond.index()];
            if r.order() != BondOrder::Double {
                continue;
            }
            match r.stereo() {
                BondStereo::None => {}
                BondStereo::Any => product.topology.bonds[row].set_stereo(BondStereo::Any)?,
                _ => forward_bond_stereo(product, input, mapping, id, r_bond)?,
            }
        }
    }
    Ok(())
}

pub(crate) fn correct_matched_chirality(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
    reactant_atom: usize,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: checkAndCorrectChiralityOfMatchingAtomsInProduct
    // RDKit❗✔️: void checkAndCorrectChiralityOfMatchingAtomsInProduct(
    // RDKit❗✔️:     const ROMol &reactant, unsigned reactantAtomIdx, const Atom &reactantAtom,
    // RDKit❗✔️:     RWMOL_SPTR product, ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   for (unsigned i = 0; i < mapping->reactProdAtomMap[reactantAtomIdx].size();
    // RDKit❗✔️:        i++) {
    // RDKit❗✔️:     unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
    // RDKit❗✔️:     Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗✔️:
    // RDKit❗✔️:     int inversionFlag = 0;
    // RDKit❗✔️:     productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗✔️:                                   inversionFlag);
    // RDKit❗✔️:     // if stereochemistry wasn't present in the reactant or if we're
    // RDKit❗✔️:     // either creating or destroying stereo we don't mess with this
    // RDKit❗✔️:     if (reactantAtom.getChiralTag() == Atom::CHI_UNSPECIFIED ||
    // RDKit❗✔️:         reactantAtom.getChiralTag() == Atom::CHI_OTHER || inversionFlag > 2) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // we can only do something sensible here if the degree in the reactants
    // RDKit❗✔️:     // and products differs by at most one
    // RDKit❗✔️:     if (reactantAtom.getDegree() < 3 || productAtom->getDegree() < 3 ||
    // RDKit❗✔️:         std::abs(static_cast<int>(reactantAtom.getDegree()) -
    // RDKit❗✔️:                  static_cast<int>(productAtom->getDegree())) > 1) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     unsigned int nUnknown = 0;
    // RDKit❗✔️:     // get the order of the bonds around the atom in the reactant:
    // RDKit❗✔️:     INT_LIST rOrder;
    // RDKit❗✔️:     for (const auto &nbri :
    // RDKit❗✔️:          boost::make_iterator_range(reactant.getAtomBonds(&reactantAtom))) {
    // RDKit❗✔️:       rOrder.push_back(reactant[nbri]->getIdx());
    // RDKit❗✔️:     }
    // RDKit❗✔️:     INT_LIST pOrder;
    // RDKit❗✔️:     for (const auto &nbri :
    // RDKit❗✔️:          boost::make_iterator_range(product->getAtomNeighbors(productAtom))) {
    // RDKit❗✔️:       if (mapping->prodReactAtomMap.find(nbri) ==
    // RDKit❗✔️:               mapping->prodReactAtomMap.end() ||
    // RDKit❗✔️:           !reactant.getBondBetweenAtoms(reactantAtom.getIdx(),
    // RDKit❗✔️:                                         mapping->prodReactAtomMap[nbri])) {
    // RDKit❗✔️:         ++nUnknown;
    // RDKit❗✔️:         // if there's more than one bond in the product that doesn't
    // RDKit❗✔️:         // correspond to anything in the reactant, we're also doomed
    // RDKit❗✔️:         if (nUnknown > 1) {
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         // otherwise, add a -1 to the bond order that we'll fill in later
    // RDKit❗✔️:         pOrder.push_back(-1);
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         const Bond *rBond = reactant.getBondBetweenAtoms(
    // RDKit❗✔️:             reactantAtom.getIdx(), mapping->prodReactAtomMap[nbri]);
    // RDKit❗✔️:         CHECK_INVARIANT(rBond, "expected reactant bond not found");
    // RDKit❗✔️:         pOrder.push_back(rBond->getIdx());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (nUnknown == 1) {
    // RDKit❗✔️:       if (reactantAtom.getDegree() == productAtom->getDegree()) {
    // RDKit❗✔️:         // there's a reactant bond that hasn't yet been accounted for:
    // RDKit❗✔️:         int unmatchedBond = -1;
    // RDKit❗✔️:
    // RDKit❗✔️:         for (const auto rBond : reactant.atomBonds(&reactantAtom)) {
    // RDKit❗✔️:           if (std::find(pOrder.begin(), pOrder.end(), rBond->getIdx()) ==
    // RDKit❗✔️:               pOrder.end()) {
    // RDKit❗✔️:             unmatchedBond = rBond->getIdx();
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         // what must be true at this point:
    // RDKit❗✔️:         //  1) there's a -1 in pOrder that we'll substitute for
    // RDKit❗✔️:         //  2) unmatchedBond contains the index of the substitution
    // RDKit❗✔️:         auto bPos = std::find(pOrder.begin(), pOrder.end(), -1);
    // RDKit❗✔️:         if (unmatchedBond >= 0 && bPos != pOrder.end()) {
    // RDKit❗✔️:           *bPos = unmatchedBond;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         nUnknown = 0;
    // RDKit❗✔️:         CHECK_INVARIANT(
    // RDKit❗✔️:             std::find(pOrder.begin(), pOrder.end(), -1) == pOrder.end(),
    // RDKit❗✔️:             "extra unmapped atom");
    // RDKit❗✔️:       } else if (productAtom->getDegree() > reactantAtom.getDegree()) {
    // RDKit❗✔️:         // the product has an extra bond. we can just remove the -1 from the
    // RDKit❗✔️:         // list:
    // RDKit❗✔️:         auto bPos = std::find(pOrder.begin(), pOrder.end(), -1);
    // RDKit❗✔️:         pOrder.erase(bPos);
    // RDKit❗✔️:         nUnknown = 0;
    // RDKit❗✔️:         CHECK_INVARIANT(
    // RDKit❗✔️:             std::find(pOrder.begin(), pOrder.end(), -1) == pOrder.end(),
    // RDKit❗✔️:             "extra unmapped atom");
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!nUnknown) {
    // RDKit❗✔️:       if (reactantAtom.getDegree() > productAtom->getDegree()) {
    // RDKit❗✔️:         // we lost a bond from the reactant.
    // RDKit❗✔️:         // we can just remove the unmatched reactant bond from the list
    // RDKit❗✔️:         INT_LIST::iterator rOrderIter = rOrder.begin();
    // RDKit❗✔️:         while (rOrderIter != rOrder.end() && rOrder.size() > pOrder.size()) {
    // RDKit❗✔️:           // we may invalidate the iterator so keep track of what comes next:
    // RDKit❗✔️:           auto thisOne = rOrderIter++;
    // RDKit❗✔️:           if (std::find(pOrder.begin(), pOrder.end(), *thisOne) ==
    // RDKit❗✔️:               pOrder.end()) {
    // RDKit❗✔️:             // not in the products:
    // RDKit❗✔️:             rOrder.erase(thisOne);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       productAtom->setChiralTag(reactantAtom.getChiralTag());
    // RDKit❗✔️:       int nSwaps = countSwapsToInterconvert(rOrder, pOrder);
    // RDKit❗✔️:       bool invert = false;
    // RDKit❗✔️:       if (nSwaps % 2) {
    // RDKit❗✔️:         invert = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       int inversionFlag;
    // RDKit❗✔️:       if (productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗✔️:                                         inversionFlag) &&
    // RDKit❗✔️:           inversionFlag == 1) {
    // RDKit❗✔️:         invert = !invert;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (invert) {
    // RDKit❗✔️:         productAtom->invertChirality();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Same small degree-local ordered lists and source quadratic permutation
    // matching, with one mutable product atom at a time, no topology clone.
    let r = &input.topology.atoms[reactant_atom];
    let reactant_neighbors = input.topology.adjacency.neighbors_of(reactant_atom);
    let matches = mapping
        .reactant_to_product
        .get(&reactant_atom)
        .ok_or_else(|| {
            invariant(
                "checkAndCorrectChiralityOfMatchingAtomsInProduct",
                "mapped atom missing",
                Some(reactant_atom),
                None,
                None,
            )
        })?;
    for &row in matches {
        let flag = inversion_flag(&product.topology.atoms[row])?.unwrap_or(0);
        if matches!(r.chiral_tag(), ChiralTag::Unspecified | ChiralTag::Other) || flag > 2 {
            continue;
        }
        let pn = &product.neighbors[row];
        if reactant_neighbors.len() < 3
            || pn.len() < 3
            || reactant_neighbors.len().abs_diff(pn.len()) > 1
        {
            continue;
        }
        let mut unknown = 0;
        let mut r_order: Vec<i32> = reactant_neighbors
            .iter()
            .map(|n| n.bond.index() as i32)
            .collect();
        let mut p_order = Vec::new();
        for n in pn {
            let bond = mapping
                .product_to_reactant
                .get(&n.atom_index)
                .and_then(|other| reactant_neighbors.iter().find(|r| r.atom_index == *other))
                .map(|n| n.bond);
            if let Some(bond) = bond {
                p_order.push(bond.index() as i32);
            } else {
                unknown += 1;
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
                    p_order[position] = bond;
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
                invert_atom_chirality(&mut product.topology.atoms[row]);
            }
        }
    }
    Ok(())
}

pub(crate) fn correct_unmatched_chirality(
    product: &mut ProductBuilder,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
    atoms: &[usize],
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: checkAndCorrectChiralityOfProduct
    // RDKit❗✔️: void checkAndCorrectChiralityOfProduct(
    // RDKit❗✔️:     const std::vector<const Atom *> &chiralAtomsToCheck, RWMOL_SPTR product,
    // RDKit❗✔️:     ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   for (auto reactantAtom : chiralAtomsToCheck) {
    // RDKit❗✔️:     CHECK_INVARIANT(reactantAtom->getChiralTag() != Atom::CHI_UNSPECIFIED,
    // RDKit❗✔️:                     "missing atom chirality.");
    // RDKit❗✔️:     const auto reactAtomDegree =
    // RDKit❗✔️:         reactantAtom->getOwningMol().getAtomDegree(reactantAtom);
    // RDKit❗✔️:     for (unsigned i = 0;
    // RDKit❗✔️:          i < mapping->reactProdAtomMap[reactantAtom->getIdx()].size(); i++) {
    // RDKit❗✔️:       unsigned productAtomIdx =
    // RDKit❗✔️:           mapping->reactProdAtomMap[reactantAtom->getIdx()][i];
    // RDKit❗✔️:       Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗✔️:       CHECK_INVARIANT(
    // RDKit❗✔️:           reactantAtom->getChiralTag() == productAtom->getChiralTag(),
    // RDKit❗✔️:           "invalid product chirality.");
    // RDKit❗✔️:
    // RDKit❗✔️:       if (reactAtomDegree != product->getAtomDegree(productAtom)) {
    // RDKit❗✔️:         // If the number of bonds to the atom has changed in the course of the
    // RDKit❗✔️:         // reaction we're lost, so remove chirality.
    // RDKit❗✔️:         //  A word of explanation here: the atoms in the chiralAtomsToCheck
    // RDKit❗✔️:         //  set are not explicitly mapped atoms of the reaction, so we really
    // RDKit❗✔️:         //  have no idea what to do with this case. At the moment I'm not even
    // RDKit❗✔️:         //  really sure how this could happen, but better safe than sorry.
    // RDKit❗✔️:         productAtom->setChiralTag(Atom::CHI_UNSPECIFIED);
    // RDKit❗✔️:       } else if (reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit❗✔️:                  reactantAtom->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit❗✔️:         // this will contain the indices of product bonds in the
    // RDKit❗✔️:         // reactant order:
    // RDKit❗✔️:         INT_LIST newOrder;
    // RDKit❗✔️:         ROMol::OEDGE_ITER beg, end;
    // RDKit❗✔️:         boost::tie(beg, end) =
    // RDKit❗✔️:             reactantAtom->getOwningMol().getAtomBonds(reactantAtom);
    // RDKit❗✔️:         while (beg != end) {
    // RDKit❗✔️:           const Bond *reactantBond = reactantAtom->getOwningMol()[*beg];
    // RDKit❗✔️:           unsigned int oAtomIdx =
    // RDKit❗✔️:               reactantBond->getOtherAtomIdx(reactantAtom->getIdx());
    // RDKit❗✔️:           CHECK_INVARIANT(mapping->reactProdAtomMap.find(oAtomIdx) !=
    // RDKit❗✔️:                               mapping->reactProdAtomMap.end(),
    // RDKit❗✔️:                           "other atom from bond not mapped.");
    // RDKit❗✔️:           const Bond *productBond;
    // RDKit❗✔️:           unsigned neighborBondIdx = mapping->reactProdAtomMap[oAtomIdx][i];
    // RDKit❗✔️:           productBond = product->getBondBetweenAtoms(productAtom->getIdx(),
    // RDKit❗✔️:                                                      neighborBondIdx);
    // RDKit❗✔️:           CHECK_INVARIANT(productBond, "no matching bond found in product");
    // RDKit❗✔️:           newOrder.push_back(productBond->getIdx());
    // RDKit❗✔️:           ++beg;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         int nSwaps = productAtom->getPerturbationOrder(newOrder);
    // RDKit❗✔️:         if (nSwaps % 2) {
    // RDKit❗✔️:           productAtom->invertChirality();
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // not tetrahedral chirality, don't do anything.
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }  // end of loop over chiralAtomsToCheck
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    for &reactant_atom in atoms {
        let r = &input.topology.atoms[reactant_atom];
        if r.chiral_tag() == ChiralTag::Unspecified {
            return Err(invariant(
                "checkAndCorrectChiralityOfProduct",
                "missing atom chirality.",
                Some(reactant_atom),
                None,
                None,
            ));
        }
        let rn = input.topology.adjacency.neighbors_of(reactant_atom);
        let rows = mapping
            .reactant_to_product
            .get(&reactant_atom)
            .ok_or_else(|| {
                invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "mapped atom missing",
                    Some(reactant_atom),
                    None,
                    None,
                )
            })?;
        for (duplicate, &row) in rows.iter().enumerate() {
            if product.topology.atoms[row].chiral_tag() != r.chiral_tag() {
                return Err(invariant(
                    "checkAndCorrectChiralityOfProduct",
                    "invalid product chirality.",
                    Some(reactant_atom),
                    Some(row),
                    None,
                ));
            }
            if rn.len() != product.neighbors[row].len() {
                product.topology.atoms[row].set_chiral_tag(ChiralTag::Unspecified);
            } else if matches!(
                r.chiral_tag(),
                ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
            ) {
                let mut new_order = Vec::with_capacity(rn.len());
                for n in rn {
                    let other = mapping
                        .reactant_to_product
                        .get(&n.atom_index)
                        .and_then(|rows| rows.get(duplicate))
                        .copied()
                        .ok_or_else(|| {
                            invariant(
                                "checkAndCorrectChiralityOfProduct",
                                "other atom from bond not mapped.",
                                Some(n.atom_index),
                                Some(row),
                                Some(n.bond),
                            )
                        })?;
                    let bond = product.bond_between(row, other).ok_or_else(|| {
                        invariant(
                            "checkAndCorrectChiralityOfProduct",
                            "no matching bond found in product",
                            Some(n.atom_index),
                            Some(row),
                            None,
                        )
                    })?;
                    new_order.push(bond);
                }
                let current_order: Vec<BondId> =
                    product.neighbors[row].iter().map(|n| n.bond).collect();
                // Atom::getPerturbationOrder calls countSwaps(probe, ref), in
                // this argument order; use the unique CORE permutation owner.
                if count_swaps_to_interconvert(&new_order, &current_order)? % 2 != 0 {
                    invert_atom_chirality(&mut product.topology.atoms[row]);
                }
            }
        }
    }
    Ok(())
}

fn new_group(source: &StereoGroup, atoms: Vec<AtomId>) -> StereoGroup {
    // RDKit❗✔️: StereoGroup(StereoGroupType grouptype, std::vector<Atom *> &&atoms,
    // RDKit❗✔️:             std::vector<Bond *> &&bonds, unsigned readId = 0);
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
    // RDKit❗✔️: void copyEnhancedStereoGroups(const ROMol &reactant, RWMOL_SPTR product,
    // RDKit❗✔️:                               const ReactantProductAtomMapping &mapping) {
    // RDKit❗✔️:   std::vector<StereoGroup> new_stereo_groups;
    // RDKit❗✔️:   for (const auto &sg : reactant.getStereoGroups()) {
    // RDKit❗✔️:     std::vector<Atom *> atoms;
    // RDKit❗✔️:     std::vector<Bond *> bonds;
    // RDKit❗✔️:     for (const auto &reactantAtom : sg.getAtoms()) {
    // RDKit❗✔️:       auto productAtoms = mapping.reactProdAtomMap.find(reactantAtom->getIdx());
    // RDKit❗✔️:       if (productAtoms == mapping.reactProdAtomMap.end()) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       for (auto &productAtomIdx : productAtoms->second) {
    // RDKit❗✔️:         auto productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗✔️:         // If chirality destroyed by the reaction, skip the atom
    // RDKit❗✔️:         if (productAtom->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         // If chirality defined explicitly by the reaction, skip the atom
    // RDKit❗✔️:         int flagVal = 0;
    // RDKit❗✔️:         productAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗✔️:                                       flagVal);
    // RDKit❗✔️:         if (flagVal == 4) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         atoms.push_back(productAtom);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (!atoms.empty()) {
    // RDKit❗✔️:       new_stereo_groups.emplace_back(sg.getGroupType(), std::move(atoms),
    // RDKit❗✔️:                                      std::move(bonds), sg.getReadId());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // Although we have added storage, and canonicalization of Atropisomers,
    // RDKit❗✔️:   // searching is not yet supported.  When it is, we will need to copy
    // RDKit❗✔️:   // bond-part of the SG groups to the products as appropriate.
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!new_stereo_groups.empty()) {
    // RDKit❗✔️:     auto &existing_sg = product->getStereoGroups();
    // RDKit❗✔️:     new_stereo_groups.insert(new_stereo_groups.end(), existing_sg.begin(),
    // RDKit❗✔️:                              existing_sg.end());
    // RDKit❗✔️:     product->setStereoGroups(std::move(new_stereo_groups));
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut groups = Vec::new();
    for group in &input.topology.stereo_groups {
        let mut atoms = Vec::new();
        for atom in group.atoms() {
            let Some(rows) = mapping.reactant_to_product.get(&atom.index()) else {
                continue;
            };
            for &row in rows {
                let p = &product.topology.atoms[row];
                if p.chiral_tag() == ChiralTag::Unspecified {
                    continue;
                }
                if inversion_flag(p)? == Some(4) {
                    continue;
                }
                atoms.push(AtomId::new(row));
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
    // RDKit❗✔️: void copyTemplateStereoGroupsToMol(const ROMol &templateMol,
    // RDKit❗✔️:                                    RWMOL_SPTR product) {
    // RDKit❗✔️:   const auto &stereoGroups = templateMol.getStereoGroups();
    // RDKit❗✔️:   if (stereoGroups.empty()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> atomsInTemplateStereoGroups(product->getNumAtoms());
    // RDKit❗✔️:   std::vector<StereoGroup> newStereoGroups;
    // RDKit❗✔️:   for (const auto &sg : stereoGroups) {
    // RDKit❗✔️:     bool keepIt = true;
    // RDKit❗✔️:     std::vector<Atom *> atoms;
    // RDKit❗✔️:     for (const auto &atom : sg.getAtoms()) {
    // RDKit❗✔️:       if (auto mapNum = atom->getAtomMapNum()) {
    // RDKit❗✔️:         for (auto productAtom : product->atoms()) {
    // RDKit❗✔️:           int oldMapNum = 0;
    // RDKit❗✔️:           if (productAtom->getPropIfPresent(common_properties::reactionMapNum,
    // RDKit❗✔️:                                             oldMapNum) &&
    // RDKit❗✔️:               oldMapNum == mapNum) {
    // RDKit❗✔️:             atoms.push_back(productAtom);
    // RDKit❗✔️:             atomsInTemplateStereoGroups.set(productAtom->getIdx());
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         keepIt = false;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keepIt && !atoms.empty()) {
    // RDKit❗✔️:       std::vector<Bond *> bonds;
    // RDKit❗✔️:       newStereoGroups.emplace_back(sg.getGroupType(), std::move(atoms),
    // RDKit❗✔️:                                    std::move(bonds), sg.getReadId());
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (!newStereoGroups.empty()) {
    // RDKit❗✔️:     // remove any stereo groups that are already present in the product (these
    // RDKit❗✔️:     // were copied over from the reactant in copyEnhancedStereoGroups()) and
    // RDKit❗✔️:     // that overlap with the added ones
    // RDKit❗✔️:     for (const auto &productSG : product->getStereoGroups()) {
    // RDKit❗✔️:       unsigned int nOverlappingAtoms = 0;
    // RDKit❗✔️:       for (const auto atom : productSG.getAtoms()) {
    // RDKit❗✔️:         if (atomsInTemplateStereoGroups[atom->getIdx()]) {
    // RDKit❗✔️:           ++nOverlappingAtoms;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (!nOverlappingAtoms) {
    // RDKit❗✔️:         // no overlapping atoms, we can just keep the stereogroup.
    // RDKit❗✔️:         newStereoGroups.push_back(productSG);
    // RDKit❗✔️:       } else if (nOverlappingAtoms < productSG.getAtoms().size()) {
    // RDKit❗✔️:         // some of the atoms in the stereo group are not already there
    // RDKit❗✔️:         // in the product, we need to split the stereo group
    // RDKit❗✔️:         std::vector<Atom *> newAtoms;
    // RDKit❗✔️:         for (const auto atom : productSG.getAtoms()) {
    // RDKit❗✔️:           if (!atomsInTemplateStereoGroups[atom->getIdx()]) {
    // RDKit❗✔️:             newAtoms.push_back(atom);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:         std::vector<Bond *> newBonds;
    // RDKit❗✔️:         newStereoGroups.emplace_back(productSG.getGroupType(),
    // RDKit❗✔️:                                      std::move(newAtoms), std::move(newBonds),
    // RDKit❗✔️:                                      productSG.getReadId());
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // else: all atoms in the stereo group are already there, we can skip it
    // RDKit❗✔️:     }
    // RDKit❗✔️:     product->setStereoGroups(std::move(newStereoGroups));
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    if template.stereo_groups().is_empty() {
        return Ok(());
    }
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
                        marked[p.id().index()] = true;
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
        for group in std::mem::take(&mut product.topology.stereo_groups) {
            let overlaps = group.atoms().iter().filter(|a| marked[a.index()]).count();
            if overlaps == 0 {
                groups.push(group);
            } else if overlaps < group.atoms().len() {
                let atoms = group
                    .atoms()
                    .iter()
                    .copied()
                    .filter(|a| !marked[a.index()])
                    .collect();
                groups.push(new_group(&group, atoms));
            }
        }
        product.topology.stereo_groups = cosmolkit_model::merge_absolute_stereo_groups(groups);
    }
    Ok(())
}

pub(crate) fn propagate_coordinates(
    points: &mut Vec<[f64; 3]>,
    is_3d: &mut bool,
    atom_count: usize,
    input: &ReactionInput<'_>,
    mapping: &ReactantProductMapping,
    selection: ReactionCoordinateSelection,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: generateProductConformers
    // RDKit❗✔️: void generateProductConformers(Conformer *productConf, const ROMol &reactant,
    // RDKit❗✔️:                                ReactantProductAtomMapping *mapping) {
    // RDKit❗✔️:   if (!reactant.getNumConformers()) {
    // RDKit❗✔️:     return;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const Conformer &reactConf = reactant.getConformer();
    // RDKit❗✔️:   if (reactConf.is3D()) {
    // RDKit❗✔️:     productConf->set3D(true);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (std::map<unsigned int, std::vector<unsigned int>>::const_iterator pr =
    // RDKit❗✔️:            mapping->reactProdAtomMap.begin();
    // RDKit❗✔️:        pr != mapping->reactProdAtomMap.end(); ++pr) {
    // RDKit❗✔️:     std::vector<unsigned> prodIdxs = pr->second;
    // RDKit❗✔️:     if (prodIdxs.size() > 1) {
    // RDKit❗✔️:       BOOST_LOG(rdWarningLog) << "reactant atom match more than one product "
    // RDKit❗✔️:                                  "atom, coordinates need to be revised\n";
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // is this reliable when multiple product atom mapping occurs????
    // RDKit❗✔️:     for (unsigned int prodIdx : prodIdxs) {
    // RDKit❗✔️:       productConf->setAtomPos(prodIdx, reactConf.getAtomPos(pr->first));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Caller performs the source resize before each reagent; new rows remain
    // exactly zero. D4 resolution reuses the existing unique CX selector.
    points.resize(atom_count, [0.0; 3]);
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
    for (&reactant_atom, rows) in &mapping.reactant_to_product {
        if rows.len() > 1 {
            eprintln!(
                "reactant atom match more than one product atom, coordinates need to be revised"
            );
        }
        let point = match source {
            cosmolkit_smiles::CoordinateSource::ThreeD(conf) => {
                conf.coordinates().get(reactant_atom).copied()
            }
            cosmolkit_smiles::CoordinateSource::TwoD(conf) => conf
                .coordinates()
                .get(reactant_atom)
                .map(|p| [p[0], p[1], 0.0]),
        }
        .ok_or_else(|| {
            invariant(
                "generateProductConformers",
                "reactant conformer atom out of range",
                Some(reactant_atom),
                None,
                None,
            )
        })?;
        for &row in rows {
            let dest = points.get_mut(row).ok_or_else(|| {
                invariant(
                    "generateProductConformers",
                    "product conformer atom out of range",
                    Some(reactant_atom),
                    Some(row),
                    None,
                )
            })?;
            *dest = point;
        }
    }
    Ok(())
}
