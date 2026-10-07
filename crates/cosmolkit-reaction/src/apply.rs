use crate::materialize::{invariant, inversion_flag, update_from_template};
use crate::{
    Reaction, ReactionApplyError, ReactionApplyParams, ReactionInput, ReactionRole,
    ReactionValidationParams,
};
use cosmolkit_model::{
    AtomId, BondSpec, QueryGraph, TopologyBatchEdit, TopologyBlock, TopologyMapping,
};
use cosmolkit_types::{BondOrder, ChiralTag};
use std::{
    borrow::Cow,
    collections::{BTreeMap, VecDeque},
};

/// Detached changes only: the sole runtime remaps coordinates/properties and
/// marks stale cache facts before atomic commit. No-match carries no block copy.
#[doc(hidden)]
pub struct ReactionApplyChanges {
    pub change: Option<(TopologyBlock, TopologyMapping)>,
    /// Exact source operation bool, independent of graph equality.
    pub changed: bool,
    /// RWMol::commitBatchEdit clears computed molecule properties on removal.
    pub clears_computed_properties: bool,
}

fn map_number(
    template: &QueryGraph,
    atom: usize,
    role: ReactionRole,
) -> Result<u32, ReactionApplyError> {
    let value = template.atom(atom).ok_or_else(|| {
        invariant(
            "run_Reactant",
            "template atom out of range",
            None,
            Some(atom),
            None,
        )
    })?;
    Ok(crate::validation::atom_map(value, role, 0)
        .map_err(crate::ReactionProductError::from)?
        .unwrap_or(0) as u32)
}
fn template_bond(template: &QueryGraph, begin: usize, end: usize) -> Option<usize> {
    template.adjacency()[begin]
        .iter()
        .find(|(other, _)| *other == end)
        .map(|(_, bond)| *bond)
}
fn identify_removed(
    template: &QueryGraph,
    matches: &[usize],
    preserved: &BTreeMap<u32, usize>,
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: identifyAtomsInReactantTemplateNotProductTemplate
    // RDKit❗✔️: void identifyAtomsInReactantTemplateNotProductTemplate(
    // RDKit❗✔️:     const ROMol &reactant, boost::dynamic_bitset<> &atoms,
    // RDKit❗✔️:     std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗✔️:     const MatchVectType &reactantMatch) {
    // RDKit❗✔️:   for (const auto atom : reactant.atoms()) {
    // RDKit❗✔️:     if (atom->getAtomMapNum()) {
    // RDKit❗✔️:       if (reactantProductMap.find(atom->getAtomMapNum()) ==
    // RDKit❗✔️:           reactantProductMap.end()) {
    // RDKit❗✔️:         // atom map not present in product
    // RDKit❗✔️:         atoms.set(reactantMatch[atom->getIdx()].second);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       // unmapped atoms in the reactants are lost in the products:
    // RDKit❗✔️:       atoms.set(reactantMatch[atom->getIdx()].second);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    for row in 0..template.num_atoms() {
        let map = map_number(template, row, ReactionRole::Reactant)?;
        if map == 0 || !preserved.contains_key(&map) {
            remove[matches[row]] = true;
        }
    }
    Ok(())
}

fn traverse_removed(
    input: ReactionInput<'_>,
    template: &QueryGraph,
    matches: &[usize],
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: traverseToFindAtomsToRemove
    // RDKit❗✔️: void traverseToFindAtomsToRemove(const ROMol &reactant, const ROMol &templ,
    // RDKit❗✔️:                                  boost::dynamic_bitset<> &atoms,
    // RDKit❗✔️:                                  const MatchVectType &reactantMatch) {
    // RDKit❗✔️:   // toRemove marks both atoms that need to be removed and those we can traverse
    // RDKit❗✔️:   // to
    // RDKit❗✔️:   boost::dynamic_bitset<> toRemove = ~atoms;
    // RDKit❗✔️:   for (const auto &tpl : reactantMatch) {
    // RDKit❗✔️:     toRemove.reset(tpl.second);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (const auto &tpl : reactantMatch) {
    // RDKit❗✔️:     std::deque<const Atom *> toConsider;
    // RDKit❗✔️:     if (templ.getAtomWithIdx(tpl.first)->getAtomMapNum() &&
    // RDKit❗✔️:         !atoms[tpl.second]) {
    // RDKit❗✔️:       toConsider.push_back(reactant.getAtomWithIdx(tpl.second));
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     while (!toConsider.empty()) {
    // RDKit❗✔️:       auto atom = toConsider.back();
    // RDKit❗✔️:       toConsider.pop_back();
    // RDKit❗✔️:       toRemove.reset(atom->getIdx());
    // RDKit❗✔️:       for (const auto nbr : reactant.atomNeighbors(atom)) {
    // RDKit❗✔️:         if (toRemove[nbr->getIdx()]) {
    // RDKit❗✔️:           toConsider.push_front(nbr);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   atoms |= toRemove;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut to_remove: Vec<bool> = remove.iter().map(|marked| !marked).collect();
    for &row in matches {
        to_remove[row] = false;
    }
    for (template_atom, &matched) in matches.iter().enumerate() {
        let mut pending = VecDeque::new();
        if map_number(template, template_atom, ReactionRole::Reactant)? != 0 && !remove[matched] {
            pending.push_back(matched);
        }
        while let Some(atom) = pending.pop_back() {
            to_remove[atom] = false;
            for neighbor in input.topology.adjacency.neighbors_of(atom) {
                if to_remove[neighbor.atom_index] {
                    pending.push_front(neighbor.atom_index);
                }
            }
        }
    }
    for (marked, unretained) in remove.iter_mut().zip(to_remove) {
        *marked |= unretained;
    }
    Ok(())
}

fn update_atoms(
    edit: &mut TopologyBatchEdit,
    input: ReactionInput<'_>,
    reactant_template: &QueryGraph,
    product_template: &QueryGraph,
    product_maps: &BTreeMap<u32, usize>,
    preserved: &BTreeMap<u32, usize>,
    matches: &[usize],
) -> Result<bool, ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateAtomsModifiedByReaction
    // RDKit❗✔️: bool updateAtomsModifiedByReaction(
    // RDKit❗✔️:     RWMol &reactant, const ROMOL_SPTR reactantTemplate,
    // RDKit❗✔️:     const ROMOL_SPTR productTemplate,
    // RDKit❗✔️:     const std::map<unsigned int, unsigned int> &productAtomMap,
    // RDKit❗✔️:     const std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗✔️:     const MatchVectType &match) {
    // RDKit❗✔️:   bool molModified = false;
    // RDKit❗✔️:   for (const auto &pr : reactantProductMap) {
    // RDKit❗✔️:     const auto rAtom = reactantTemplate->getAtomWithIdx(pr.second);
    // RDKit❗✔️:     const auto pAtom =
    // RDKit❗✔️:         productTemplate->getAtomWithIdx(productAtomMap.at(pr.first));
    // RDKit❗✔️:     const auto atom = reactant.getAtomWithIdx(match[pr.second].second);
    // RDKit❗✔️:     if (rAtom->getAtomicNum() != pAtom->getAtomicNum() &&
    // RDKit❗✔️:         (pAtom->getAtomicNum() || !pAtom->hasQuery())) {
    // RDKit❗✔️:       atom->setAtomicNum(pAtom->getAtomicNum());
    // RDKit❗✔️:       molModified = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (ReactionRunnerUtils::updatePropsFromImplicitProps(pAtom, atom)) {
    // RDKit❗✔️:       molModified = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // check if we need to modify stereo
    // RDKit❗✔️:     int molInversionFlag;
    // RDKit❗✔️:     if (pAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗✔️:                                 molInversionFlag)) {
    // RDKit❗✔️:       auto atomTag = atom->getChiralTag();
    // RDKit❗✔️:       switch (molInversionFlag) {
    // RDKit❗✔️:         case 0:  // no chiral impact, do nothing
    // RDKit❗✔️:         case 2:  // retention, do nothing
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 1:
    // RDKit❗✔️:           // inversion
    // RDKit❗✔️:           if (atomTag != Atom::ChiralType::CHI_OTHER &&
    // RDKit❗✔️:               atomTag != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit❗✔️:             atom->invertChirality();
    // RDKit❗✔️:             molModified = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 3:
    // RDKit❗✔️:           // destroy
    // RDKit❗✔️:           atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit❗✔️:           molModified = true;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         case 4:
    // RDKit❗✔️:           // create
    // RDKit❗✔️:           atom->setChiralTag(pAtom->getChiralTag());
    // RDKit❗✔️:           molModified = true;
    // RDKit❗✔️:           // check swaps
    // RDKit❗✔️:           {
    // RDKit❗✔️:             std::vector<int> porder;
    // RDKit❗✔️:             for (const auto nbrAtom : productTemplate->atomNeighbors(pAtom)) {
    // RDKit❗✔️:               if (nbrAtom->getAtomMapNum()) {
    // RDKit❗✔️:                 porder.push_back(nbrAtom->getAtomMapNum());
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:             // get the ordered vect of atom map numbers for the neighbors
    // RDKit❗✔️:             // of atom
    // RDKit❗✔️:             std::vector<int> aorder;
    // RDKit❗✔️:             for (auto aidx :
    // RDKit❗✔️:                  boost::make_iterator_range(reactant.getAtomNeighbors(atom))) {
    // RDKit❗✔️:               auto miter = std::find_if(
    // RDKit❗✔️:                   match.begin(), match.end(), [aidx](const auto &pr) {
    // RDKit❗✔️:                     return static_cast<unsigned int>(pr.second) == aidx;
    // RDKit❗✔️:                   });
    // RDKit❗✔️:               if (miter != match.end()) {
    // RDKit❗✔️:                 auto rNbr = reactantTemplate->getAtomWithIdx(miter->first);
    // RDKit❗✔️:                 if (rNbr->getAtomMapNum()) {
    // RDKit❗✔️:                   aorder.push_back(rNbr->getAtomMapNum());
    // RDKit❗✔️:                 }
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:             if (porder.size() == aorder.size()) {
    // RDKit❗✔️:               auto nswaps = countSwapsToInterconvert(aorder, porder);
    // RDKit❗✔️:               if (nswaps % 2) {
    // RDKit❗✔️:                 atom->invertChirality();
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           }
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         default:
    // RDKit❗✔️:           BOOST_LOG(rdWarningLog)
    // RDKit❗✔️:               << "unrecognized chiral inversion/retention flag "
    // RDKit❗✔️:                  "on product atom ignored\n";
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return molModified;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut changed = false;
    for (&map, &reactant_row) in preserved {
        let r = &reactant_template.atoms()[reactant_row];
        let p = &product_template.atoms()[product_maps[&map]];
        let atom_row = matches[reactant_row];
        let atom = edit.atom_mut(AtomId::new(atom_row))?;
        if r.atomic_number() != p.atomic_number()
            && (p.atomic_number() != 0 || p.predicate_is_carrier_derived())
        {
            let element = p.element().ok_or_else(|| {
                invariant(
                    "updateAtomsModifiedByReaction",
                    "product atomic number is not a modeled element",
                    Some(atom_row),
                    Some(p.id().index()),
                    None,
                )
            })?;
            atom.set_element(element);
            changed = true;
        }
        if update_from_template(p, atom)? {
            changed = true;
        }
        if let Some(flag) = inversion_flag(p)? {
            match flag {
                0 | 2 => {}
                1 => {
                    if !matches!(atom.chiral_tag(), ChiralTag::Other | ChiralTag::Unspecified) {
                        cosmolkit_core::invert_atom_chirality(atom);
                        changed = true;
                    }
                }
                3 => {
                    atom.set_chiral_tag(ChiralTag::Unspecified);
                    changed = true;
                }
                4 => {
                    atom.set_chiral_tag(p.chiral_tag());
                    changed = true;
                    let mut p_order = Vec::new();
                    for &(neighbor, _) in &product_template.adjacency()[p.id().index()] {
                        let map = map_number(product_template, neighbor, ReactionRole::Product)?;
                        if map != 0 {
                            p_order.push(map as i32);
                        }
                    }
                    let mut a_order = Vec::new();
                    for neighbor in input.topology.adjacency.neighbors_of(atom_row) {
                        if let Some(query_row) =
                            matches.iter().position(|row| *row == neighbor.atom_index)
                        {
                            let map =
                                map_number(reactant_template, query_row, ReactionRole::Reactant)?;
                            if map != 0 {
                                a_order.push(map as i32);
                            }
                        }
                    }
                    if p_order.len() == a_order.len()
                        && cosmolkit_core::count_swaps_to_interconvert(&a_order, &p_order)
                            .map_err(crate::ReactionProductError::from)?
                            % 2
                            != 0
                    {
                        cosmolkit_core::invert_atom_chirality(atom);
                    }
                }
                _ => eprintln!(
                    "unrecognized chiral inversion/retention flag on product atom ignored"
                ),
            }
        }
    }
    Ok(changed)
}

fn update_bonds(
    edit: &mut TopologyBatchEdit,
    reactant_template: &QueryGraph,
    product_template: &QueryGraph,
    product_maps: &BTreeMap<u32, usize>,
    preserved: &BTreeMap<u32, usize>,
    matches: &[usize],
) -> Result<(bool, bool), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateBondsModifiedByReaction
    // RDKit❗✔️: bool updateBondsModifiedByReaction(
    // RDKit❗✔️:     RWMol &reactant, const ROMOL_SPTR reactantTemplate,
    // RDKit❗✔️:     const ROMOL_SPTR productTemplate,
    // RDKit❗✔️:     const std::map<unsigned int, unsigned int> &productAtomMap,
    // RDKit❗✔️:     const std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗✔️:     const MatchVectType &match) {
    // RDKit❗✔️:   bool molModified = false;
    // RDKit❗✔️:   for (const auto &pr : reactantProductMap) {
    // RDKit❗✔️:     const auto rAtom = reactantTemplate->getAtomWithIdx(pr.second);
    // RDKit❗✔️:     const auto pAtom =
    // RDKit❗✔️:         productTemplate->getAtomWithIdx(productAtomMap.at(pr.first));
    // RDKit❗✔️:     const auto atom = reactant.getAtomWithIdx(match[pr.second].second);
    // RDKit❗✔️:     for (const auto nbr : productTemplate->atomNeighbors(pAtom)) {
    // RDKit❗✔️:       if (nbr->getAtomMapNum() &&
    // RDKit❗✔️:           reactantProductMap.find(nbr->getAtomMapNum()) !=
    // RDKit❗✔️:               reactantProductMap.end()) {
    // RDKit❗✔️:         const auto pBond = productTemplate->getBondBetweenAtoms(pAtom->getIdx(),
    // RDKit❗✔️:                                                                 nbr->getIdx());
    // RDKit❗✔️:         ASSERT_INVARIANT(pBond,
    // RDKit❗✔️:                          "missing bond between known neighbors in product");
    // RDKit❗✔️:         const auto rBond = reactantTemplate->getBondBetweenAtoms(
    // RDKit❗✔️:             rAtom->getIdx(), reactantProductMap.at(nbr->getAtomMapNum()));
    // RDKit❗✔️:         if (rBond) {
    // RDKit❗✔️:           if (pBond->getBondType() != Bond::BondType::UNSPECIFIED &&
    // RDKit❗✔️:               pBond->getBondType() != rBond->getBondType()) {
    // RDKit❗✔️:             const auto bond = reactant.getBondBetweenAtoms(
    // RDKit❗✔️:                 match[rBond->getBeginAtomIdx()].second,
    // RDKit❗✔️:                 match[rBond->getEndAtomIdx()].second);
    // RDKit❗✔️:             ASSERT_INVARIANT(
    // RDKit❗✔️:                 bond, "missing bond between known neighbors in reactant");
    // RDKit❗✔️:             bond->setBondType(pBond->getBondType());
    // RDKit❗✔️:             molModified = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           // there was no corresponding bond in the reactant template, was there
    // RDKit❗✔️:           // one in the reactant?
    // RDKit❗✔️:           const auto bond = reactant.getBondBetweenAtoms(
    // RDKit❗✔️:               match[rAtom->getIdx()].second,
    // RDKit❗✔️:               match[reactantProductMap.at(nbr->getAtomMapNum())].second);
    // RDKit❗✔️:           if (!bond) {
    // RDKit❗✔️:             auto begIdx = match[reactantProductMap.at(
    // RDKit❗✔️:                                     pBond->getBeginAtom()->getAtomMapNum())]
    // RDKit❗✔️:                               .second;
    // RDKit❗✔️:             auto endIdx = match[reactantProductMap.at(
    // RDKit❗✔️:                                     pBond->getEndAtom()->getAtomMapNum())]
    // RDKit❗✔️:                               .second;
    // RDKit❗✔️:
    // RDKit❗✔️:             ReactionRunnerUtils::addBondToProduct(*pBond, reactant, begIdx,
    // RDKit❗✔️:                                                   endIdx);
    // RDKit❗✔️:             molModified = true;
    // RDKit❗✔️:           } else if (bond->getBondType() != pBond->getBondType()) {
    // RDKit❗✔️:             bond->setBondType(pBond->getBondType());
    // RDKit❗✔️:             molModified = true;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     // now look for bonds which were in the reactant template but are not in the
    // RDKit❗✔️:     // product template
    // RDKit❗✔️:     for (const auto nbr : reactantTemplate->atomNeighbors(rAtom)) {
    // RDKit❗✔️:       if (nbr->getAtomMapNum() &&
    // RDKit❗✔️:           productAtomMap.find(nbr->getAtomMapNum()) != productAtomMap.end() &&
    // RDKit❗✔️:           !productTemplate->getBondBetweenAtoms(
    // RDKit❗✔️:               pAtom->getIdx(), productAtomMap.at(nbr->getAtomMapNum()))) {
    // RDKit❗✔️:         // remove the bond in the reactant
    // RDKit❗✔️:         reactant.removeBond(atom->getIdx(), match[nbr->getIdx()].second);
    // RDKit❗✔️:         molModified = true;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return molModified;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut changed = false;
    let mut removed = false;
    for (&map, &reactant_row) in preserved {
        let product_row = product_maps[&map];
        let atom_row = matches[reactant_row];
        for &(neighbor, bond_row) in &product_template.adjacency()[product_row] {
            let neighbor_map = map_number(product_template, neighbor, ReactionRole::Product)?;
            if neighbor_map == 0 {
                continue;
            }
            let Some(&neighbor_reactant_row) = preserved.get(&neighbor_map) else {
                continue;
            };
            let p = &product_template.bonds()[bond_row];
            if let Some(r_bond_row) =
                template_bond(reactant_template, reactant_row, neighbor_reactant_row)
            {
                let r = &reactant_template.bonds()[r_bond_row];
                if p.bond().order() != BondOrder::Unspecified
                    && p.bond().order() != r.bond().order()
                {
                    let bond = edit
                        .bond_between_atoms(
                            AtomId::new(matches[r.begin().index()]),
                            AtomId::new(matches[r.end().index()]),
                        )?
                        .ok_or_else(|| {
                            invariant(
                                "updateBondsModifiedByReaction",
                                "missing bond between known neighbors in reactant",
                                Some(atom_row),
                                None,
                                None,
                            )
                        })?;
                    edit.bond_mut(bond)?.set_order(p.bond().order());
                    changed = true;
                }
            } else {
                let neighbor_row = matches[neighbor_reactant_row];
                if let Some(bond) =
                    edit.bond_between_atoms(AtomId::new(atom_row), AtomId::new(neighbor_row))?
                {
                    let b = edit.bond_mut(bond)?;
                    if b.order() != p.bond().order() {
                        b.set_order(p.bond().order());
                        changed = true;
                    }
                } else {
                    let begin_map =
                        map_number(product_template, p.begin().index(), ReactionRole::Product)?;
                    let end_map =
                        map_number(product_template, p.end().index(), ReactionRole::Product)?;
                    let begin = matches[*preserved.get(&begin_map).ok_or_else(|| {
                        invariant(
                            "updateBondsModifiedByReaction",
                            "begin product map not preserved",
                            None,
                            Some(p.begin().index()),
                            None,
                        )
                    })?];
                    let end = matches[*preserved.get(&end_map).ok_or_else(|| {
                        invariant(
                            "updateBondsModifiedByReaction",
                            "end product map not preserved",
                            None,
                            Some(p.end().index()),
                            None,
                        )
                    })?];
                    // RDKit❗✔️: Bond *addBondToProduct(const Bond &origB, RWMol &product,
                    // RDKit❗✔️:                        unsigned int begAtomIdx, unsigned int endAtomIdx) {
                    // RDKit❗✔️:   if (!origB.hasQuery()) {
                    // RDKit❗✔️:     auto idx = product.addBond(begAtomIdx, endAtomIdx, origB.getBondType());
                    // RDKit❗✔️:     return product.getBondWithIdx(idx - 1);
                    // RDKit❗✔️:   } else {
                    // RDKit❗✔️:     QueryBond *qbond = new QueryBond(origB.getBondType());
                    // RDKit❗✔️:     qbond->setBeginAtomIdx(begAtomIdx);
                    // RDKit❗✔️:     qbond->setEndAtomIdx(endAtomIdx);
                    // RDKit❗✔️:     qbond->setQuery(origB.getQuery()->copy());
                    // RDKit❗✔️:     bool takeOwnership = true;
                    // RDKit❗✔️:     product.addBond(qbond, takeOwnership);
                    // RDKit❗✔️:     return qbond;
                    // RDKit❗✔️:   }
                    // RDKit❗✔️: }
                    let spec =
                        BondSpec::new(AtomId::new(begin), AtomId::new(end), p.bond().order());
                    let spec = if p.predicate_is_carrier_derived() {
                        spec
                    } else {
                        spec.with_query(p.predicate().clone())
                    };
                    edit.add_bond(spec)?;
                    changed = true;
                }
            }
        }
        for &(neighbor, _) in &reactant_template.adjacency()[reactant_row] {
            let neighbor_map = map_number(reactant_template, neighbor, ReactionRole::Reactant)?;
            if neighbor_map != 0
                && let Some(&neighbor_product_row) = product_maps.get(&neighbor_map)
                && template_bond(product_template, product_row, neighbor_product_row).is_none()
            {
                // RWMol keeps pending removed bonds visible until batch commit.
                if let Some(bond) =
                    edit.bond_between_atoms(AtomId::new(atom_row), AtomId::new(matches[neighbor]))?
                {
                    edit.remove_bond(bond)?;
                    removed = true;
                }
                changed = true;
            }
        }
    }
    Ok((changed, removed))
}

#[doc(hidden)]
pub fn apply_reaction(
    reaction: &Reaction,
    input: ReactionInput<'_>,
    params: &ReactionApplyParams,
) -> Result<ReactionApplyChanges, ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactant
    // RDKit❗✔️: bool run_Reactant(const ChemicalReaction &rxn, RWMol &reactant,
    // RDKit❗✔️:                   bool removeUnmatchedAtoms) {
    // RDKit❗✔️:   PRECONDITION(rxn.getNumReactantTemplates() == 1,
    // RDKit❗✔️:                "only one reactant supported");
    // RDKit❗✔️:   PRECONDITION(rxn.getNumProductTemplates() == 1, "only one product supported");
    // RDKit❗✔️:   if (!rxn.isInitialized()) {
    // RDKit❗✔️:     throw ChemicalReactionException(
    // RDKit❗✔️:         "initMatchers() must be called before runReactants()");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const unsigned int reactantIdx = 0;
    // RDKit❗✔️:   const auto reactantTemplate = rxn.getReactants()[reactantIdx];
    // RDKit❗✔️:   const auto productTemplate = rxn.getProducts()[0];
    // RDKit❗✔️:
    // RDKit❗✔️:   std::map<unsigned int, unsigned int>
    // RDKit❗✔️:       productAtomMap;  // atom mapnum -> product atom index
    // RDKit❗✔️:   for (const auto atom : productTemplate->atoms()) {
    // RDKit❗✔️:     if (atom->getAtomMapNum()) {
    // RDKit❗✔️:       productAtomMap[atom->getAtomMapNum()] = atom->getIdx();
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::map<unsigned int, unsigned int>
    // RDKit❗✔️:       reactantProductMap;  // atom mapnum -> reactant atom index, for atoms
    // RDKit❗✔️:                            // which are also mapped in the product
    // RDKit❗✔️:   for (const auto atom : reactantTemplate->atoms()) {
    // RDKit❗✔️:     if (atom->getAtomMapNum()) {
    // RDKit❗✔️:       if (productAtomMap.find(atom->getAtomMapNum()) != productAtomMap.end()) {
    // RDKit❗✔️:         reactantProductMap[atom->getAtomMapNum()] = atom->getIdx();
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // we don't support reactions with unmapped or new atoms in the products
    // RDKit❗✔️:   for (const auto atom : productTemplate->atoms()) {
    // RDKit❗✔️:     if (!atom->getAtomMapNum() ||
    // RDKit❗✔️:         reactantProductMap.find(atom->getAtomMapNum()) ==
    // RDKit❗✔️:             reactantProductMap.end()) {
    // RDKit❗✔️:       throw ChemicalReactionException(
    // RDKit❗✔️:           "single component reactions which add atoms in the product "
    // RDKit❗✔️:           "are not supported");
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   auto reactantMatch = ReactionRunnerUtils::getReactantMatchesToTemplate(
    // RDKit❗✔️:       reactant, *reactantTemplate, 1, rxn.getSubstructParams());
    // RDKit❗✔️:   if (reactantMatch.empty()) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   const auto &match = reactantMatch[0];
    // RDKit❗✔️:
    // RDKit❗✔️:   // we now have a match for the reactant, so we can work on it
    // RDKit❗✔️:   // start by marking atoms which are in the reactant template, but not in the
    // RDKit❗✔️:   // product template for removal
    // RDKit❗✔️:   boost::dynamic_bitset<> atomsToRemove(reactant.getNumAtoms());
    // RDKit❗✔️:   // finds atoms in the reactantTemplate which aren't in the productTemplate
    // RDKit❗✔️:   ReactionRunnerUtils::identifyAtomsInReactantTemplateNotProductTemplate(
    // RDKit❗✔️:       *reactantTemplate, atomsToRemove, reactantProductMap, match);
    // RDKit❗✔️:   if (removeUnmatchedAtoms) {
    // RDKit❗✔️:     // identify atoms which did not match something in the reactant template but
    // RDKit❗✔️:     // which should be removed from the molecule
    // RDKit❗✔️:     ReactionRunnerUtils::traverseToFindAtomsToRemove(
    // RDKit❗✔️:         reactant, *reactantTemplate, atomsToRemove, match);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   bool molModified = false;
    // RDKit❗✔️:   reactant.beginBatchEdit();
    // RDKit❗✔️:
    // RDKit❗✔️:   if (updateAtomsModifiedByReaction(reactant, reactantTemplate, productTemplate,
    // RDKit❗✔️:                                     productAtomMap, reactantProductMap,
    // RDKit❗✔️:                                     match)) {
    // RDKit❗✔️:     molModified = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (updateBondsModifiedByReaction(reactant, reactantTemplate, productTemplate,
    // RDKit❗✔️:                                     productAtomMap, reactantProductMap,
    // RDKit❗✔️:                                     match)) {
    // RDKit❗✔️:     molModified = true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // remove atoms which aren't transferred to the products (marked above)
    // RDKit❗✔️:   if (atomsToRemove.count()) {
    // RDKit❗✔️:     molModified = true;
    // RDKit❗✔️:     for (unsigned int i = 0; i < atomsToRemove.size(); ++i) {
    // RDKit❗✔️:       if (atomsToRemove[i]) {
    // RDKit❗✔️:         reactant.removeAtom(i);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   reactant.commitBatchEdit();
    // RDKit❗✔️:   return molModified;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    if reaction.num_reactant_templates() != 1 || reaction.num_product_templates() != 1 {
        return Err(ReactionApplyError::ApplicabilityArity {
            reactants: reaction.num_reactant_templates(),
            products: reaction.num_product_templates(),
        });
    }
    let reaction = if reaction.is_initialized() {
        Cow::Borrowed(reaction)
    } else {
        Cow::Owned(crate::initialize_reaction(
            reaction,
            &ReactionValidationParams::default(),
        )?)
    };
    let r = &reaction.reactant_templates()[0];
    let p = &reaction.product_templates()[0];
    let mut product_maps = BTreeMap::new();
    for row in 0..p.num_atoms() {
        let map = map_number(p, row, ReactionRole::Product)?;
        if map != 0 {
            product_maps.insert(map, row);
        }
    }
    let mut preserved = BTreeMap::new();
    for row in 0..r.num_atoms() {
        let map = map_number(r, row, ReactionRole::Reactant)?;
        if map != 0 && product_maps.contains_key(&map) {
            preserved.insert(map, row);
        }
    }
    for row in 0..p.num_atoms() {
        let map = map_number(p, row, ReactionRole::Product)?;
        if map == 0 || !preserved.contains_key(&map) {
            return Err(ReactionApplyError::AddsProductAtom {
                atom: p.atoms()[row].id(),
            });
        }
    }
    let matches = crate::matching::matches_to_template(input, r, 1, reaction.match_params(), 0, 0)?;
    let Some(matches) = matches.first() else {
        return Ok(ReactionApplyChanges {
            change: None,
            changed: false,
            clears_computed_properties: false,
        });
    };
    let mut remove = vec![false; input.topology.atoms.len()];
    identify_removed(r, matches, &preserved, &mut remove)?;
    if params.remove_unmatched_atoms {
        traverse_removed(input, r, matches, &mut remove)?;
    }
    // One owned topology copy after a successful match, never all molecule
    // blocks. The existing canonical MODEL batch editor performs compaction.
    let mut edit = input.topology.clone().into_batch_edit()?;
    let mut changed = update_atoms(&mut edit, input, r, p, &product_maps, &preserved, matches)?;
    let (bonds_changed, bonds_removed) =
        update_bonds(&mut edit, r, p, &product_maps, &preserved, matches)?;
    changed |= bonds_changed;
    let atoms_removed = remove.iter().any(|marked| *marked);
    if atoms_removed {
        changed = true;
        for (row, marked) in remove.into_iter().enumerate() {
            if marked {
                edit.remove_atom(AtomId::new(row))?;
            }
        }
    }
    if !changed {
        edit.abort();
        return Ok(ReactionApplyChanges {
            change: None,
            changed: false,
            clears_computed_properties: false,
        });
    }
    let change = edit.finish()?;
    // The source does not updatePropertyCache here. The runtime invalidates
    // stale derived facts and remaps every coordinate set on validated commit.
    Ok(ReactionApplyChanges {
        change: Some(change),
        changed,
        clears_computed_properties: atoms_removed || bonds_removed,
    })
}
