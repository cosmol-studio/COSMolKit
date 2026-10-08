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

fn identify_removed(
    template: &QueryGraph,
    matches: &[usize],
    preserved: &BTreeMap<u32, usize>,
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    identify_removed_source(
        template,
        |index| {
            let target = *matches.get(index).ok_or_else(|| {
                invariant(
                    "identifyAtomsInReactantTemplateNotProductTemplate",
                    "template match index out of range",
                    None,
                    Some(index),
                    None,
                )
            })?;
            Ok(crate::materialize::source_u32("matched reactant row", target)? as i32)
        },
        preserved,
        remove,
    )
}

fn identify_removed_source(
    template: &QueryGraph,
    mut match_second: impl FnMut(usize) -> Result<i32, ReactionApplyError>,
    preserved: &BTreeMap<u32, usize>,
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: identifyAtomsInReactantTemplateNotProductTemplate
    // RDKit❗❌: void identifyAtomsInReactantTemplateNotProductTemplate(
    // RDKit❗❌:     const ROMol &reactant, boost::dynamic_bitset<> &atoms,
    // RDKit❗❌:     std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗❌:     const MatchVectType &reactantMatch) {
    // RDKit❗❌:   for (const auto atom : reactant.atoms()) {
    // RDKit❗❌:     if (atom->getAtomMapNum()) {
    // RDKit❗❌:       if (reactantProductMap.find(atom->getAtomMapNum()) ==
    // RDKit❗❌:           reactantProductMap.end()) {
    // RDKit❗❌:         // atom map not present in product
    // RDKit❗❌:         atoms.set(reactantMatch[atom->getIdx()].second);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       // unmapped atoms in the reactants are lost in the products:
    // RDKit❗❌:       atoms.set(reactantMatch[atom->getIdx()].second);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗✔️:   int getAtomMapNum() const {
    // RDKit❗✔️:     int mapno = 0;
    // RDKit❗✔️:     getPropIfPresent(common_properties::molAtomMapNumber, mapno);
    // RDKit❗✔️:     return mapno;
    // RDKit❗✔️:   }
    // Iterate physical template pointers, but index the source match vector
    // using each Atom::getIdx, only for atoms actually marked for removal.
    // Native calls getAtomMapNum twice for nonzero maps; preserve both reads.
    // No pair first-index read, canonicalization or graph validation is added.
    for (row, atom) in template.atoms().iter().enumerate() {
        let map = map_number(template, row, ReactionRole::Reactant)?;
        if map != 0 {
            let key = map_number(template, row, ReactionRole::Reactant)?;
            if preserved.contains_key(&key) {
                continue;
            }
        }
        let index = crate::materialize::source_u32("template atom", atom.id().index())? as usize;
        let second = match_second(index)?;
        let bit = second as usize;
        *remove.get_mut(bit).ok_or_else(|| {
            invariant(
                "identifyAtomsInReactantTemplateNotProductTemplate",
                "removed atom bit index out of range",
                Some(bit),
                None,
                None,
            )
        })? = true;
    }
    Ok(())
}

fn traverse_removed(
    input: ReactionInput<'_>,
    template: &QueryGraph,
    matches: &[usize],
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    let matched = matches.iter().enumerate().map(|(query, &target)| {
        Ok((
            crate::materialize::source_u32("matched query row", query)? as i32,
            crate::materialize::source_u32("matched reactant row", target)? as i32,
        ))
    });
    traverse_removed_source(input, template, matched, remove)
}

fn traverse_removed_source(
    input: ReactionInput<'_>,
    template: &QueryGraph,
    matched: impl IntoIterator<Item = Result<(i32, i32), ReactionApplyError>> + Clone,
    remove: &mut [bool],
) -> Result<(), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: traverseToFindAtomsToRemove
    // RDKit❗❌: void traverseToFindAtomsToRemove(const ROMol &reactant, const ROMol &templ,
    // RDKit❗❌:                                  boost::dynamic_bitset<> &atoms,
    // RDKit❗❌:                                  const MatchVectType &reactantMatch) {
    // RDKit❗❌:   // toRemove marks both atoms that need to be removed and those we can traverse
    // RDKit❗❌:   // to
    // RDKit❗❌:   boost::dynamic_bitset<> toRemove = ~atoms;
    // RDKit❗❌:   for (const auto &tpl : reactantMatch) {
    // RDKit❗❌:     toRemove.reset(tpl.second);
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto &tpl : reactantMatch) {
    // RDKit❗❌:     std::deque<const Atom *> toConsider;
    // RDKit❗❌:     if (templ.getAtomWithIdx(tpl.first)->getAtomMapNum() &&
    // RDKit❗❌:         !atoms[tpl.second]) {
    // RDKit❗❌:       toConsider.push_back(reactant.getAtomWithIdx(tpl.second));
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     while (!toConsider.empty()) {
    // RDKit❗❌:       auto atom = toConsider.back();
    // RDKit❗❌:       toConsider.pop_back();
    // RDKit❗❌:       toRemove.reset(atom->getIdx());
    // RDKit❗❌:       for (const auto nbr : reactant.atomNeighbors(atom)) {
    // RDKit❗❌:         if (toRemove[nbr->getIdx()]) {
    // RDKit❗❌:           toConsider.push_front(nbr);
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   atoms |= toRemove;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   const auto res = d_graph[vd];
    // RDKit❗✔️:
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit❗✔️: };
    // Preserve Native FIFO deque operations and mark-on-pop. Multiple pending
    // references to the same atom are deliberately not deduplicated.
    // Byte flags retain the known packed-mask footprint gap; pointers project
    // directly to borrowed Atom references, without topology/atom clones.
    let mut to_remove: Vec<bool> = remove.iter().map(|marked| !marked).collect();
    for pair in matched.clone() {
        let (_, second) = pair?;
        let bit = second as usize;
        *to_remove.get_mut(bit).ok_or_else(|| {
            invariant(
                "traverseToFindAtomsToRemove",
                "matched reset bit index out of range",
                Some(bit),
                None,
                None,
            )
        })? = false;
    }
    for pair in matched {
        let (first, second) = pair?;
        let query = (first as u32) as usize;
        let mut pending = VecDeque::new();
        if map_number(template, query, ReactionRole::Reactant)? != 0
            && !*remove.get(second as usize).ok_or_else(|| {
                invariant(
                    "traverseToFindAtomsToRemove",
                    "seed bit index out of range",
                    Some(second as usize),
                    None,
                    None,
                )
            })?
        {
            let row = (second as u32) as usize;
            let atom = input.topology.atoms.get(row).ok_or_else(|| {
                invariant(
                    "traverseToFindAtomsToRemove",
                    "seed reactant atom row missing",
                    Some(row),
                    None,
                    None,
                )
            })?;
            pending.push_back(atom);
        }
        while let Some(atom) = pending.pop_back() {
            let index =
                crate::materialize::source_u32("reactant atom", atom.id().index())? as usize;
            *to_remove.get_mut(index).ok_or_else(|| {
                invariant(
                    "traverseToFindAtomsToRemove",
                    "visited reset bit index out of range",
                    Some(index),
                    None,
                    None,
                )
            })? = false;
            let neighbors = input
                .topology
                .adjacency
                .try_neighbors_of(index)
                .ok_or_else(|| {
                    invariant(
                        "traverseToFindAtomsToRemove",
                        "reactant adjacency row missing",
                        Some(index),
                        None,
                        None,
                    )
                })?;
            for neighbor in neighbors {
                let atom = input
                    .topology
                    .atoms
                    .get(neighbor.atom_index)
                    .ok_or_else(|| {
                        invariant(
                            "traverseToFindAtomsToRemove",
                            "neighbor reactant atom row missing",
                            Some(neighbor.atom_index),
                            None,
                            Some(neighbor.bond),
                        )
                    })?;
                let index =
                    crate::materialize::source_u32("reactant atom", atom.id().index())? as usize;
                if *to_remove.get(index).ok_or_else(|| {
                    invariant(
                        "traverseToFindAtomsToRemove",
                        "neighbor bit index out of range",
                        Some(index),
                        None,
                        Some(neighbor.bond),
                    )
                })? {
                    pending.push_front(atom);
                }
            }
        }
    }
    // Both masks retain identical length by construction; only this final
    // source OR changes caller state, including when no matches were present.
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
    update_atoms_source(
        edit,
        input,
        reactant_template,
        product_template,
        product_maps,
        preserved,
        matches.len(),
        |index| {
            let target = *matches.get(index).ok_or_else(|| {
                invariant(
                    "updateAtomsModifiedByReaction",
                    "template match index out of range",
                    None,
                    Some(index),
                    None,
                )
            })?;
            Ok((
                crate::materialize::source_u32("matched query row", index)? as i32,
                crate::materialize::source_u32("matched reactant row", target)? as i32,
            ))
        },
    )
}

fn update_atoms_source(
    edit: &mut TopologyBatchEdit,
    input: ReactionInput<'_>,
    reactant_template: &QueryGraph,
    product_template: &QueryGraph,
    product_maps: &BTreeMap<u32, usize>,
    preserved: &BTreeMap<u32, usize>,
    match_count: usize,
    mut match_at: impl FnMut(usize) -> Result<(i32, i32), ReactionApplyError>,
) -> Result<bool, ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateAtomsModifiedByReaction
    // RDKit❗❌: bool updateAtomsModifiedByReaction(
    // RDKit❗❌:     RWMol &reactant, const ROMOL_SPTR reactantTemplate,
    // RDKit❗❌:     const ROMOL_SPTR productTemplate,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &productAtomMap,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗❌:     const MatchVectType &match) {
    // RDKit❗❌:   bool molModified = false;
    // RDKit❗❌:   for (const auto &pr : reactantProductMap) {
    // RDKit❗❌:     const auto rAtom = reactantTemplate->getAtomWithIdx(pr.second);
    // RDKit❗❌:     const auto pAtom =
    // RDKit❗❌:         productTemplate->getAtomWithIdx(productAtomMap.at(pr.first));
    // RDKit❗❌:     const auto atom = reactant.getAtomWithIdx(match[pr.second].second);
    // RDKit❗❌:     if (rAtom->getAtomicNum() != pAtom->getAtomicNum() &&
    // RDKit❗❌:         (pAtom->getAtomicNum() || !pAtom->hasQuery())) {
    // RDKit❗❌:       atom->setAtomicNum(pAtom->getAtomicNum());
    // RDKit❗❌:       molModified = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     if (ReactionRunnerUtils::updatePropsFromImplicitProps(pAtom, atom)) {
    // RDKit❗❌:       molModified = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     // check if we need to modify stereo
    // RDKit❗❌:     int molInversionFlag;
    // RDKit❗❌:     if (pAtom->getPropIfPresent(common_properties::molInversionFlag,
    // RDKit❗❌:                                 molInversionFlag)) {
    // RDKit❗❌:       auto atomTag = atom->getChiralTag();
    // RDKit❗❌:       switch (molInversionFlag) {
    // RDKit❗❌:         case 0:  // no chiral impact, do nothing
    // RDKit❗❌:         case 2:  // retention, do nothing
    // RDKit❗❌:           break;
    // RDKit❗❌:         case 1:
    // RDKit❗❌:           // inversion
    // RDKit❗❌:           if (atomTag != Atom::ChiralType::CHI_OTHER &&
    // RDKit❗❌:               atomTag != Atom::ChiralType::CHI_UNSPECIFIED) {
    // RDKit❗❌:             atom->invertChirality();
    // RDKit❗❌:             molModified = true;
    // RDKit❗❌:           }
    // RDKit❗❌:           break;
    // RDKit❗❌:         case 3:
    // RDKit❗❌:           // destroy
    // RDKit❗❌:           atom->setChiralTag(Atom::ChiralType::CHI_UNSPECIFIED);
    // RDKit❗❌:           molModified = true;
    // RDKit❗❌:           break;
    // RDKit❗❌:         case 4:
    // RDKit❗❌:           // create
    // RDKit❗❌:           atom->setChiralTag(pAtom->getChiralTag());
    // RDKit❗❌:           molModified = true;
    // RDKit❗❌:           // check swaps
    // RDKit❗❌:           {
    // RDKit❗❌:             std::vector<int> porder;
    // RDKit❗❌:             for (const auto nbrAtom : productTemplate->atomNeighbors(pAtom)) {
    // RDKit❗❌:               if (nbrAtom->getAtomMapNum()) {
    // RDKit❗❌:                 porder.push_back(nbrAtom->getAtomMapNum());
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             // get the ordered vect of atom map numbers for the neighbors
    // RDKit❗❌:             // of atom
    // RDKit❗❌:             std::vector<int> aorder;
    // RDKit❗❌:             for (auto aidx :
    // RDKit❗❌:                  boost::make_iterator_range(reactant.getAtomNeighbors(atom))) {
    // RDKit❗❌:               auto miter = std::find_if(
    // RDKit❗❌:                   match.begin(), match.end(), [aidx](const auto &pr) {
    // RDKit❗❌:                     return static_cast<unsigned int>(pr.second) == aidx;
    // RDKit❗❌:                   });
    // RDKit❗❌:               if (miter != match.end()) {
    // RDKit❗❌:                 auto rNbr = reactantTemplate->getAtomWithIdx(miter->first);
    // RDKit❗❌:                 if (rNbr->getAtomMapNum()) {
    // RDKit❗❌:                   aorder.push_back(rNbr->getAtomMapNum());
    // RDKit❗❌:                 }
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:             if (porder.size() == aorder.size()) {
    // RDKit❗❌:               auto nswaps = countSwapsToInterconvert(aorder, porder);
    // RDKit❗❌:               if (nswaps % 2) {
    // RDKit❗❌:                 atom->invertChirality();
    // RDKit❗❌:               }
    // RDKit❗❌:             }
    // RDKit❗❌:           }
    // RDKit❗❌:           break;
    // RDKit❗❌:         default:
    // RDKit❗❌:           BOOST_LOG(rdWarningLog)
    // RDKit❗❌:               << "unrecognized chiral inversion/retention flag "
    // RDKit❗❌:                  "on product atom ignored\n";
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return molModified;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
    // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
    // RDKit❗✔️:   const auto res = d_graph[vd];
    // RDKit❗✔️:
    // RDKit❗✔️:   POSTCONDITION(res, "");
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // RDKit❗✔️: ROMol::ADJ_ITER_PAIR ROMol::getAtomNeighbors(Atom const *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return boost::adjacent_vertices(at->getIdx(), d_graph);
    // RDKit❗✔️: };
    // Source map order and per-atom mutation order are preserved. Reuse sole
    // template-property, countSwaps and invertChirality owners; never infer
    // molModified from equality or the bool returned by invertChirality.
    // Known detached editor/property helper costs remain ❌; the source
    // neighbor * physical-match find_if loops are intentionally retained.
    let mut changed = false;
    for (&map, &reactant_row) in preserved {
        let reactant_row =
            crate::materialize::source_u32("reactant template atom", reactant_row)? as usize;
        let r = reactant_template.atoms().get(reactant_row).ok_or_else(|| {
            invariant(
                "updateAtomsModifiedByReaction",
                "reactant template atom row missing",
                Some(reactant_row),
                None,
                None,
            )
        })?;
        let product_row = *product_maps.get(&map).ok_or_else(|| {
            invariant(
                "updateAtomsModifiedByReaction",
                "product atom map key missing",
                Some(reactant_row),
                None,
                None,
            )
        })?;
        let product_row =
            crate::materialize::source_u32("product template atom", product_row)? as usize;
        let p = product_template.atoms().get(product_row).ok_or_else(|| {
            invariant(
                "updateAtomsModifiedByReaction",
                "product template atom row missing",
                Some(reactant_row),
                Some(product_row),
                None,
            )
        })?;
        let (_, second) = match_at(reactant_row)?;
        let atom_row = (second as u32) as usize;
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
                        cosmolkit_core::invert_atom_chirality(atom)
                            .map_err(crate::ReactionProductError::from)?;
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
                    let product_id =
                        crate::materialize::source_u32("product template atom", p.id().index())?
                            as usize;
                    let product_neighbors = product_template
                        .adjacency()
                        .get(product_id)
                        .ok_or_else(|| {
                            invariant(
                                "updateAtomsModifiedByReaction",
                                "product template adjacency row missing",
                                Some(atom_row),
                                Some(product_id),
                                None,
                            )
                        })?;
                    for &(neighbor, _) in product_neighbors {
                        if map_number(product_template, neighbor, ReactionRole::Product)? != 0 {
                            p_order.push(map_number(
                                product_template,
                                neighbor,
                                ReactionRole::Product,
                            )? as i32);
                        }
                    }
                    let mut a_order = Vec::new();
                    let actual_id =
                        crate::materialize::source_u32("reactant atom", atom.id().index())?
                            as usize;
                    let reactant_neighbors = input
                        .topology
                        .adjacency
                        .try_neighbors_of(actual_id)
                        .ok_or_else(|| {
                        invariant(
                            "updateAtomsModifiedByReaction",
                            "reactant adjacency row missing",
                            Some(actual_id),
                            None,
                            None,
                        )
                    })?;
                    for neighbor in reactant_neighbors {
                        // Source find_if uses physical first-match order and
                        // compares unsigned(second) with the graph vertex;
                        // only the winning pair's first indexes the template.
                        let mut found = None;
                        for index in 0..match_count {
                            let pair = match_at(index)?;
                            if (pair.1 as u32) as usize == neighbor.atom_index {
                                found = Some(pair.0);
                                break;
                            }
                        }
                        if let Some(first) = found {
                            let query = (first as u32) as usize;
                            if map_number(reactant_template, query, ReactionRole::Reactant)? != 0 {
                                a_order.push(map_number(
                                    reactant_template,
                                    query,
                                    ReactionRole::Reactant,
                                )? as i32);
                            }
                        }
                    }
                    if p_order.len() == a_order.len()
                        && cosmolkit_core::count_swaps_to_interconvert(&a_order, &p_order)
                            .map_err(crate::ReactionProductError::from)?
                            % 2
                            != 0
                    {
                        cosmolkit_core::invert_atom_chirality(atom)
                            .map_err(crate::ReactionProductError::from)?;
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

fn apply_template_atom(
    q: &QueryGraph,
    row: usize,
) -> Result<&cosmolkit_model::QueryAtom, ReactionApplyError> {
    q.atoms().get(row).ok_or_else(|| {
        invariant(
            "updateBondsModifiedByReaction",
            "template atom row missing",
            None,
            Some(row),
            None,
        )
        .into()
    })
}
fn apply_template_neighbors(
    q: &QueryGraph,
    row: usize,
) -> Result<&[(usize, usize)], ReactionApplyError> {
    q.adjacency().get(row).map(Vec::as_slice).ok_or_else(|| {
        invariant(
            "updateBondsModifiedByReaction",
            "template adjacency row missing",
            None,
            Some(row),
            None,
        )
        .into()
    })
}
fn apply_template_edge(
    q: &QueryGraph,
    begin: usize,
    end: usize,
) -> Result<Option<&cosmolkit_model::QueryBond>, ReactionApplyError> {
    // Borrow only the queried adjacency row, after the canonical source range
    // checks. Owned NeighborRef transport introduces no edge-search algorithm.
    let mut missing = false;
    let id = cosmolkit_model::source_bond_between_atoms(
        q.atoms().len(),
        AtomId::new(begin),
        AtomId::new(end),
        || {
            let row = q.adjacency().get(begin);
            missing = row.is_none();
            row.into_iter()
                .flatten()
                .map(|&(atom_index, bond)| cosmolkit_model::NeighborRef {
                    atom_index,
                    bond: cosmolkit_model::BondId::new(bond),
                })
        },
    )?;
    if missing {
        return Err(invariant(
            "updateBondsModifiedByReaction",
            "template adjacency row missing",
            None,
            Some(begin),
            None,
        )
        .into());
    }
    id.map(|id| {
        q.bonds().get(id.index()).ok_or_else(|| {
            invariant(
                "updateBondsModifiedByReaction",
                "template edge bond row missing",
                None,
                Some(id.index()),
                None,
            )
            .into()
        })
    })
    .transpose()
}

fn update_bonds(
    edit: &mut TopologyBatchEdit,
    reactant_template: &QueryGraph,
    product_template: &QueryGraph,
    product_maps: &BTreeMap<u32, usize>,
    preserved: &BTreeMap<u32, usize>,
    matches: &[usize],
) -> Result<(bool, bool), ReactionApplyError> {
    update_bonds_source(
        edit,
        reactant_template,
        product_template,
        product_maps,
        preserved,
        |index| {
            let target = *matches.get(index).ok_or_else(|| {
                invariant(
                    "updateBondsModifiedByReaction",
                    "template match index out of range",
                    None,
                    Some(index),
                    None,
                )
            })?;
            Ok(crate::materialize::source_u32("matched reactant row", target)? as i32)
        },
    )
}

fn update_bonds_source(
    edit: &mut TopologyBatchEdit,
    reactant_template: &QueryGraph,
    product_template: &QueryGraph,
    product_maps: &BTreeMap<u32, usize>,
    preserved: &BTreeMap<u32, usize>,
    mut match_second: impl FnMut(usize) -> Result<i32, ReactionApplyError>,
) -> Result<(bool, bool), ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: updateBondsModifiedByReaction
    // RDKit❗❌: bool updateBondsModifiedByReaction(
    // RDKit❗❌:     RWMol &reactant, const ROMOL_SPTR reactantTemplate,
    // RDKit❗❌:     const ROMOL_SPTR productTemplate,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &productAtomMap,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &reactantProductMap,
    // RDKit❗❌:     const MatchVectType &match) {
    // RDKit❗❌:   bool molModified = false;
    // RDKit❗❌:   for (const auto &pr : reactantProductMap) {
    // RDKit❗❌:     const auto rAtom = reactantTemplate->getAtomWithIdx(pr.second);
    // RDKit❗❌:     const auto pAtom =
    // RDKit❗❌:         productTemplate->getAtomWithIdx(productAtomMap.at(pr.first));
    // RDKit❗❌:     const auto atom = reactant.getAtomWithIdx(match[pr.second].second);
    // RDKit❗❌:     for (const auto nbr : productTemplate->atomNeighbors(pAtom)) {
    // RDKit❗❌:       if (nbr->getAtomMapNum() &&
    // RDKit❗❌:           reactantProductMap.find(nbr->getAtomMapNum()) !=
    // RDKit❗❌:               reactantProductMap.end()) {
    // RDKit❗❌:         const auto pBond = productTemplate->getBondBetweenAtoms(pAtom->getIdx(),
    // RDKit❗❌:                                                                 nbr->getIdx());
    // RDKit❗❌:         ASSERT_INVARIANT(pBond,
    // RDKit❗❌:                          "missing bond between known neighbors in product");
    // RDKit❗❌:         const auto rBond = reactantTemplate->getBondBetweenAtoms(
    // RDKit❗❌:             rAtom->getIdx(), reactantProductMap.at(nbr->getAtomMapNum()));
    // RDKit❗❌:         if (rBond) {
    // RDKit❗❌:           if (pBond->getBondType() != Bond::BondType::UNSPECIFIED &&
    // RDKit❗❌:               pBond->getBondType() != rBond->getBondType()) {
    // RDKit❗❌:             const auto bond = reactant.getBondBetweenAtoms(
    // RDKit❗❌:                 match[rBond->getBeginAtomIdx()].second,
    // RDKit❗❌:                 match[rBond->getEndAtomIdx()].second);
    // RDKit❗❌:             ASSERT_INVARIANT(
    // RDKit❗❌:                 bond, "missing bond between known neighbors in reactant");
    // RDKit❗❌:             bond->setBondType(pBond->getBondType());
    // RDKit❗❌:             molModified = true;
    // RDKit❗❌:           }
    // RDKit❗❌:         } else {
    // RDKit❗❌:           // there was no corresponding bond in the reactant template, was there
    // RDKit❗❌:           // one in the reactant?
    // RDKit❗❌:           const auto bond = reactant.getBondBetweenAtoms(
    // RDKit❗❌:               match[rAtom->getIdx()].second,
    // RDKit❗❌:               match[reactantProductMap.at(nbr->getAtomMapNum())].second);
    // RDKit❗❌:           if (!bond) {
    // RDKit❗❌:             auto begIdx = match[reactantProductMap.at(
    // RDKit❗❌:                                     pBond->getBeginAtom()->getAtomMapNum())]
    // RDKit❗❌:                               .second;
    // RDKit❗❌:             auto endIdx = match[reactantProductMap.at(
    // RDKit❗❌:                                     pBond->getEndAtom()->getAtomMapNum())]
    // RDKit❗❌:                               .second;
    // RDKit❗❌:
    // RDKit❗❌:             ReactionRunnerUtils::addBondToProduct(*pBond, reactant, begIdx,
    // RDKit❗❌:                                                   endIdx);
    // RDKit❗❌:             molModified = true;
    // RDKit❗❌:           } else if (bond->getBondType() != pBond->getBondType()) {
    // RDKit❗❌:             bond->setBondType(pBond->getBondType());
    // RDKit❗❌:             molModified = true;
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     // now look for bonds which were in the reactant template but are not in the
    // RDKit❗❌:     // product template
    // RDKit❗❌:     for (const auto nbr : reactantTemplate->atomNeighbors(rAtom)) {
    // RDKit❗❌:       if (nbr->getAtomMapNum() &&
    // RDKit❗❌:           productAtomMap.find(nbr->getAtomMapNum()) != productAtomMap.end() &&
    // RDKit❗❌:           !productTemplate->getBondBetweenAtoms(
    // RDKit❗❌:               pAtom->getIdx(), productAtomMap.at(nbr->getAtomMapNum()))) {
    // RDKit❗❌:         // remove the bond in the reactant
    // RDKit❗❌:         reactant.removeBond(atom->getIdx(), match[nbr->getIdx()].second);
    // RDKit❗❌:         molModified = true;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return molModified;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // RDKit❗❌: Bond *addBondToProduct(const Bond &origB, RWMol &product,
    // RDKit❗❌:                        unsigned int begAtomIdx, unsigned int endAtomIdx) {
    // RDKit❗❌:   if (!origB.hasQuery()) {
    // RDKit❗❌:     auto idx = product.addBond(begAtomIdx, endAtomIdx, origB.getBondType());
    // RDKit❗❌:     return product.getBondWithIdx(idx - 1);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     QueryBond *qbond = new QueryBond(origB.getBondType());
    // RDKit❗❌:     qbond->setBeginAtomIdx(begAtomIdx);
    // RDKit❗❌:     qbond->setEndAtomIdx(endAtomIdx);
    // RDKit❗❌:     qbond->setQuery(origB.getQuery()->copy());
    // RDKit❗❌:     bool takeOwnership = true;
    // RDKit❗❌:     product.addBond(qbond, takeOwnership);
    // RDKit❗❌:     return qbond;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // RDKit❗✔️: QueryBond::QueryBond(BondType bT) : Bond(bT) {
    // RDKit❗✔️:   if (bT != Bond::UNSPECIFIED) {
    // RDKit❗✔️:     dp_query = makeBondOrderEqualsQuery(bT);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     dp_query = makeBondNullQuery();
    // RDKit❗✔️:   }
    // RDKit❗✔️: };
    // RDKit❗✔️: Bond::Bond(BondType bT) : RDProps() {
    // RDKit❗✔️:   initBond();
    // RDKit❗✔️:   d_bondType = bT;
    // RDKit❗✔️: };
    // RDKit❗✔️: void Bond::initBond() {
    // RDKit❗✔️:   d_bondType = UNSPECIFIED;
    // RDKit❗✔️:   d_dirTag = NONE;
    // RDKit❗✔️:   d_stereo = STEREONONE;
    // RDKit❗✔️:   dp_mol = nullptr;
    // RDKit❗✔️:   d_beginAtomIdx = 0;
    // RDKit❗✔️:   d_endAtomIdx = 0;
    // RDKit❗✔️:   df_isAromatic = 0;
    // RDKit❗✔️:   d_index = 0;
    // RDKit❗✔️:   df_isConjugated = 0;
    // RDKit❗✔️:   dp_stereoAtoms = nullptr;
    // RDKit❗✔️: };
    // RDKit❗✔️:   void setQuery(QUERYBOND_QUERY *what) override {
    // RDKit❗✔️:     // free up any existing query (Issue255):
    // RDKit❗✔️:     delete dp_query;
    // RDKit❗✔️:     dp_query = what;
    // RDKit❗✔️:   }
    // Keep native ordered unsigned maps, physical neighbor order and reached
    // accesses. Source bond lookup/add/remove algorithms have one MODEL owner.
    // Cost ❌: checked detached rows/property transport and byte removal masks
    // retain known gaps; the query CSR adapter streams without buffering.
    let stage = "updateBondsModifiedByReaction";
    let failure = |detail, row| invariant(stage, detail, None, Some(row), None);
    let product_row = |map: u32| -> Result<usize, ReactionApplyError> {
        Ok(crate::materialize::source_u32(
            "product template atom",
            *product_maps
                .get(&map)
                .ok_or_else(|| failure("product atom map key missing", map as usize))?,
        )? as usize)
    };
    let kept_row = |map: u32| -> Result<usize, ReactionApplyError> {
        Ok(crate::materialize::source_u32(
            "reactant template atom",
            *preserved
                .get(&map)
                .ok_or_else(|| failure("reactant product map key missing", map as usize))?,
        )? as usize)
    };
    let id = |kind, row| crate::materialize::source_u32(kind, row).map(|v| v as usize);
    let mut changed = false;
    let mut removed = false;
    for (&map, &reactant_row) in preserved {
        let reactant_row = id("reactant template atom", reactant_row)?;
        let r = apply_template_atom(reactant_template, reactant_row)?;
        let p = apply_template_atom(product_template, product_row(map)?)?;
        let actual_row = (match_second(reactant_row)? as u32) as usize;
        let actual_id = edit.atom_mut(AtomId::new(actual_row))?.id().index();
        let p_id = id("product template atom", p.id().index())?;
        for &(neighbor, _) in apply_template_neighbors(product_template, p_id)? {
            let nbr = apply_template_atom(product_template, neighbor)?;
            if map_number(product_template, neighbor, ReactionRole::Product)? == 0 {
                continue;
            }
            if !preserved.contains_key(&map_number(
                product_template,
                neighbor,
                ReactionRole::Product,
            )?) {
                continue;
            }
            let nbr_id = id("product template atom", nbr.id().index())?;
            let pb = apply_template_edge(product_template, p_id, nbr_id)?
                .ok_or_else(|| failure("missing bond between known neighbors in product", p_id))?;
            let r_id = id("reactant template atom", r.id().index())?;
            let rb = apply_template_edge(
                reactant_template,
                r_id,
                kept_row(map_number(
                    product_template,
                    neighbor,
                    ReactionRole::Product,
                )?)?,
            )?;
            if let Some(rb) = rb {
                if pb.bond().order() != BondOrder::Unspecified
                    && pb.bond().order() != rb.bond().order()
                {
                    let begin = (match_second(id("reactant template atom", rb.begin().index())?)?
                        as u32) as usize;
                    let end = (match_second(id("reactant template atom", rb.end().index())?)?
                        as u32) as usize;
                    let bond = edit
                        .bond_between_atoms(AtomId::new(begin), AtomId::new(end))?
                        .ok_or_else(|| {
                            failure(
                                "missing bond between known neighbors in reactant",
                                actual_row,
                            )
                        })?;
                    edit.bond_mut(bond)?.set_order(pb.bond().order());
                    changed = true;
                }
            } else {
                let begin = (match_second(r_id)? as u32) as usize;
                let end = (match_second(kept_row(map_number(
                    product_template,
                    neighbor,
                    ReactionRole::Product,
                )?)?)? as u32) as usize;
                if let Some(bond) = edit.bond_between_atoms(AtomId::new(begin), AtomId::new(end))? {
                    let bond = edit.bond_mut(bond)?;
                    if bond.order() != pb.bond().order() {
                        bond.set_order(pb.bond().order());
                        changed = true;
                    }
                } else {
                    let beg_atom = id("product template atom", pb.begin().index())?;
                    let begin = (match_second(kept_row(map_number(
                        product_template,
                        beg_atom,
                        ReactionRole::Product,
                    )?)?)? as u32) as usize;
                    let end_atom = id("product template atom", pb.end().index())?;
                    let end = (match_second(kept_row(map_number(
                        product_template,
                        end_atom,
                        ReactionRole::Product,
                    )?)?)? as u32) as usize;
                    if pb.predicate_is_carrier_derived() {
                        edit.add_bond_by_order(
                            AtomId::new(begin),
                            AtomId::new(end),
                            pb.bond().order(),
                        )?;
                    } else {
                        let bond = cosmolkit_model::Bond::from_spec(
                            cosmolkit_model::BondId::new(0),
                            BondSpec::new(AtomId::new(begin), AtomId::new(end), pb.bond().order())
                                .with_query(pb.predicate().clone()),
                        );
                        edit.add_bond_value_source(Cow::Owned(bond))?;
                    }
                    changed = true;
                }
            }
        }
        for &(neighbor, _) in apply_template_neighbors(
            reactant_template,
            id("reactant template atom", r.id().index())?,
        )? {
            let nbr = apply_template_atom(reactant_template, neighbor)?;
            if map_number(reactant_template, neighbor, ReactionRole::Reactant)? != 0
                && product_maps.contains_key(&map_number(
                    reactant_template,
                    neighbor,
                    ReactionRole::Reactant,
                )?)
                && apply_template_edge(
                    product_template,
                    p_id,
                    product_row(map_number(
                        reactant_template,
                        neighbor,
                        ReactionRole::Reactant,
                    )?)?,
                )?
                .is_none()
            {
                let end = (match_second(id("reactant template atom", nbr.id().index())?)? as u32)
                    as usize;
                let begin = AtomId::new(id("reactant atom", actual_id)?);
                let end = AtomId::new(end);
                let existed = edit.bond_between_atoms(begin, end)?.is_some();
                edit.remove_bond_between_atoms(begin, end)?;
                changed = true;
                removed |= existed;
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
    // ROOT-approved D1 immutable initialization projection. Source arity gates
    // precede preparation; the source function keeps its own Native init gate.
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
    apply_reaction_source(&reaction, input, params)
}

fn apply_reaction_source(
    reaction: &Reaction,
    input: ReactionInput<'_>,
    params: &ReactionApplyParams,
) -> Result<ReactionApplyChanges, ReactionApplyError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactant
    // RDKit❗❌: bool run_Reactant(const ChemicalReaction &rxn, RWMol &reactant,
    // RDKit❗❌:                   bool removeUnmatchedAtoms) {
    // RDKit❗❌:   PRECONDITION(rxn.getNumReactantTemplates() == 1,
    // RDKit❗❌:                "only one reactant supported");
    // RDKit❗❌:   PRECONDITION(rxn.getNumProductTemplates() == 1, "only one product supported");
    // RDKit❗❌:   if (!rxn.isInitialized()) {
    // RDKit❗❌:     throw ChemicalReactionException(
    // RDKit❗❌:         "initMatchers() must be called before runReactants()");
    // RDKit❗❌:   }
    // RDKit❗❌:   const unsigned int reactantIdx = 0;
    // RDKit❗❌:   const auto reactantTemplate = rxn.getReactants()[reactantIdx];
    // RDKit❗❌:   const auto productTemplate = rxn.getProducts()[0];
    // RDKit❗❌:
    // RDKit❗❌:   std::map<unsigned int, unsigned int>
    // RDKit❗❌:       productAtomMap;  // atom mapnum -> product atom index
    // RDKit❗❌:   for (const auto atom : productTemplate->atoms()) {
    // RDKit❗❌:     if (atom->getAtomMapNum()) {
    // RDKit❗❌:       productAtomMap[atom->getAtomMapNum()] = atom->getIdx();
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   std::map<unsigned int, unsigned int>
    // RDKit❗❌:       reactantProductMap;  // atom mapnum -> reactant atom index, for atoms
    // RDKit❗❌:                            // which are also mapped in the product
    // RDKit❗❌:   for (const auto atom : reactantTemplate->atoms()) {
    // RDKit❗❌:     if (atom->getAtomMapNum()) {
    // RDKit❗❌:       if (productAtomMap.find(atom->getAtomMapNum()) != productAtomMap.end()) {
    // RDKit❗❌:         reactantProductMap[atom->getAtomMapNum()] = atom->getIdx();
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // we don't support reactions with unmapped or new atoms in the products
    // RDKit❗❌:   for (const auto atom : productTemplate->atoms()) {
    // RDKit❗❌:     if (!atom->getAtomMapNum() ||
    // RDKit❗❌:         reactantProductMap.find(atom->getAtomMapNum()) ==
    // RDKit❗❌:             reactantProductMap.end()) {
    // RDKit❗❌:       throw ChemicalReactionException(
    // RDKit❗❌:           "single component reactions which add atoms in the product "
    // RDKit❗❌:           "are not supported");
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   auto reactantMatch = ReactionRunnerUtils::getReactantMatchesToTemplate(
    // RDKit❗❌:       reactant, *reactantTemplate, 1, rxn.getSubstructParams());
    // RDKit❗❌:   if (reactantMatch.empty()) {
    // RDKit❗❌:     return false;
    // RDKit❗❌:   }
    // RDKit❗❌:   const auto &match = reactantMatch[0];
    // RDKit❗❌:
    // RDKit❗❌:   // we now have a match for the reactant, so we can work on it
    // RDKit❗❌:   // start by marking atoms which are in the reactant template, but not in the
    // RDKit❗❌:   // product template for removal
    // RDKit❗❌:   boost::dynamic_bitset<> atomsToRemove(reactant.getNumAtoms());
    // RDKit❗❌:   // finds atoms in the reactantTemplate which aren't in the productTemplate
    // RDKit❗❌:   ReactionRunnerUtils::identifyAtomsInReactantTemplateNotProductTemplate(
    // RDKit❗❌:       *reactantTemplate, atomsToRemove, reactantProductMap, match);
    // RDKit❗❌:   if (removeUnmatchedAtoms) {
    // RDKit❗❌:     // identify atoms which did not match something in the reactant template but
    // RDKit❗❌:     // which should be removed from the molecule
    // RDKit❗❌:     ReactionRunnerUtils::traverseToFindAtomsToRemove(
    // RDKit❗❌:         reactant, *reactantTemplate, atomsToRemove, match);
    // RDKit❗❌:   }
    // RDKit❗❌:   bool molModified = false;
    // RDKit❗❌:   reactant.beginBatchEdit();
    // RDKit❗❌:
    // RDKit❗❌:   if (updateAtomsModifiedByReaction(reactant, reactantTemplate, productTemplate,
    // RDKit❗❌:                                     productAtomMap, reactantProductMap,
    // RDKit❗❌:                                     match)) {
    // RDKit❗❌:     molModified = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (updateBondsModifiedByReaction(reactant, reactantTemplate, productTemplate,
    // RDKit❗❌:                                     productAtomMap, reactantProductMap,
    // RDKit❗❌:                                     match)) {
    // RDKit❗❌:     molModified = true;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // remove atoms which aren't transferred to the products (marked above)
    // RDKit❗❌:   if (atomsToRemove.count()) {
    // RDKit❗❌:     molModified = true;
    // RDKit❗❌:     for (unsigned int i = 0; i < atomsToRemove.size(); ++i) {
    // RDKit❗❌:       if (atomsToRemove[i]) {
    // RDKit❗❌:         reactant.removeAtom(i);
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   reactant.commitBatchEdit();
    // RDKit❗❌:   return molModified;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    if reaction.num_reactant_templates() != 1 || reaction.num_product_templates() != 1 {
        return Err(ReactionApplyError::ApplicabilityArity {
            reactants: reaction.num_reactant_templates(),
            products: reaction.num_product_templates(),
        });
    }
    if !reaction.is_initialized() {
        return Err(crate::ReactionRunError::NeedsInitialization.into());
    }
    let r = &reaction.reactant_templates()[0];
    let p = &reaction.product_templates()[0];
    let mut product_maps = BTreeMap::new();
    for (row, atom) in p.atoms().iter().enumerate() {
        if map_number(p, row, ReactionRole::Product)? != 0 {
            let key = map_number(p, row, ReactionRole::Product)?;
            product_maps.insert(
                key,
                crate::materialize::source_u32("product template atom", atom.id().index())?
                    as usize,
            );
        }
    }
    let mut preserved = BTreeMap::new();
    for (row, atom) in r.atoms().iter().enumerate() {
        if map_number(r, row, ReactionRole::Reactant)? != 0
            && product_maps.contains_key(&map_number(r, row, ReactionRole::Reactant)?)
        {
            let key = map_number(r, row, ReactionRole::Reactant)?;
            preserved.insert(
                key,
                crate::materialize::source_u32("reactant template atom", atom.id().index())?
                    as usize,
            );
        }
    }
    for (row, atom) in p.atoms().iter().enumerate() {
        if map_number(p, row, ReactionRole::Product)? == 0
            || !preserved.contains_key(&map_number(p, row, ReactionRole::Product)?)
        {
            return Err(ReactionApplyError::AddsProductAtom { atom: atom.id() });
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
    // Source commit runs after every match, including an identity reaction.
    // Template-property writes contribute to the bool through the sole helper.
    // Copy only side blocks reached by actual removals,
    // execute their canonical source validation, then let the runtime apply
    // the returned row mapping/property-clear effect at its single boundary.
    let mut coordinates = if atoms_removed {
        input.coordinates.clone()
    } else {
        Default::default()
    };
    let mut properties = if atoms_removed || bonds_removed {
        input.properties.clone()
    } else {
        Default::default()
    };
    let change = edit.finish_source::<ReactionApplyError>(
        &mut coordinates,
        &mut properties,
        &mut |value| {
            cosmolkit_core::property_value_to_uint(value)
                .map_err(|e| crate::ReactionProductError::from(e).into())
        },
        &mut |value| {
            cosmolkit_core::property_value_to_string(value).map_err(ReactionApplyError::from)
        },
    )?;
    // The source does not updatePropertyCache here. The runtime invalidates
    // stale derived facts and remaps every coordinate set on validated commit.
    Ok(ReactionApplyChanges {
        change: changed.then_some(change),
        changed,
        clears_computed_properties: atoms_removed || bonds_removed,
    })
}

#[cfg(test)]
mod complete_identify_removed_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;
    fn graph(maps: &[Option<i32>]) -> QueryGraph {
        QueryGraph::from_parts(
            maps.iter()
                .enumerate()
                .map(|(i, map)| {
                    let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
                    if let Some(map) = map {
                        a.set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                            .unwrap();
                    }
                    a
                })
                .collect(),
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn run(
        q: &QueryGraph,
        pairs: &[(i32, i32)],
        preserved: &BTreeMap<u32, usize>,
        remove: &mut [bool],
    ) -> Result<(), ReactionApplyError> {
        identify_removed_source(
            q,
            |index| {
                Ok(pairs
                    .get(index)
                    .ok_or_else(|| {
                        invariant(
                            "identifyAtomsInReactantTemplateNotProductTemplate",
                            "template match index out of range",
                            None,
                            Some(index),
                            None,
                        )
                    })?
                    .1)
            },
            preserved,
            remove,
        )
    }
    #[test]
    fn empty_template_skips_all_match_access_and_preserves_existing_mask() {
        let mut remove = [true, false, true];
        identify_removed_source(
            &graph(&[]),
            |_| panic!("unreached"),
            &BTreeMap::new(),
            &mut remove,
        )
        .unwrap();
        assert_eq!(remove, [true, false, true]);
    }
    #[test]
    fn physical_template_order_indexes_match_by_actual_atom_id_and_ignores_pair_first_fields() {
        let mut q = graph(&[None, Some(0), Some(9), Some(10)]);
        for (a, id) in q.atoms_mut().iter_mut().zip([2, 0, 3, 1]) {
            *a = a.clone().with_id(AtomId::new(id));
        }
        let pairs = [(-900, 4), (999, 2), (3, 1), (0, 0)];
        let mut seen = vec![];
        let mut remove = [true, false, false, false, false];
        identify_removed_source(
            &q,
            |index| {
                seen.push(index);
                Ok(pairs[index].1)
            },
            &BTreeMap::from([(9, 77)]),
            &mut remove,
        )
        .unwrap();
        assert_eq!(seen, [2, 0, 1]);
        assert_eq!(remove, [true, true, true, false, true]);
    }
    #[test]
    fn preserved_positive_and_negative_maps_use_unsigned_map_keys_without_reached_match_reads() {
        let mut q = graph(&[Some(1), Some(-1)]);
        for a in q.atoms_mut() {
            *a = a.clone().with_id(AtomId::new(99));
        }
        let mut remove = [false];
        identify_removed_source(
            &q,
            |_| panic!("preserved atoms never read match"),
            &BTreeMap::from([(1, 0), (u32::MAX, 999)]),
            &mut remove,
        )
        .unwrap();
        assert_eq!(remove, [false]);
    }
    #[test]
    fn duplicate_target_bits_are_reached_in_each_physical_template_iteration() {
        let q = graph(&[None, None]);
        let mut seen = vec![];
        let mut remove = [false; 3];
        identify_removed_source(
            &q,
            |index| {
                seen.push(index);
                Ok(2)
            },
            &BTreeMap::new(),
            &mut remove,
        )
        .unwrap();
        assert_eq!(seen, [0, 1]);
        assert_eq!(remove, [false, false, true]);
    }
    #[test]
    fn property_conversion_error_keeps_prior_removed_bits_and_precedes_later_match_read() {
        let mut q = graph(&[None, Some(11)]);
        q.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::String("bad".into()))
            .unwrap();
        let mut seen = vec![];
        let mut remove = [false; 2];
        assert!(matches!(
            identify_removed_source(
                &q,
                |index| {
                    seen.push(index);
                    Ok(index as i32)
                },
                &BTreeMap::new(),
                &mut remove
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::TemplateProperty(_)
            ))
        ));
        assert_eq!(seen, [0]);
        assert_eq!(remove, [true, false]);
    }
    #[test]
    fn reached_missing_match_slot_is_structural_after_prior_bit_writes() {
        let mut q = graph(&[None, None]);
        q.atoms_mut()[1] = q.atoms()[1].clone().with_id(AtomId::new(99));
        let mut remove = [false; 2];
        assert!(matches!(
            run(&q, &[(123, 0)], &BTreeMap::new(), &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "template match index out of range",
                    product_atom: Some(99),
                    ..
                }
            ))
        ));
        assert_eq!(remove, [true, false]);
    }
    #[test]
    fn negative_and_positive_out_of_range_bit_indices_fail_at_reached_set_with_prefix_retained() {
        let q = graph(&[None, None]);
        for second in [-1, 2] {
            let mut remove = [false; 2];
            let bit = second as usize;
            assert!(
                matches!(run(&q,&[(0,0),(1,second)],&BTreeMap::new(),&mut remove),Err(ReactionApplyError::Product(crate::ReactionProductError::Invariant{stage:"identifyAtomsInReactantTemplateNotProductTemplate",detail:"removed atom bit index out of range",reactant_atom:Some(actual),..})) if actual==bit)
            );
            assert_eq!(remove, [true, false]);
        }
    }
    #[test]
    fn canonical_index_adapter_converts_only_reached_removed_atoms_and_rejects_source_uint_overflow()
     {
        let q = graph(&[Some(11), None]);
        let mut remove = [false; 2];
        identify_removed(
            &q,
            &[usize::MAX, 1],
            &BTreeMap::from([(11, 0)]),
            &mut remove,
        )
        .unwrap();
        assert_eq!(remove, [false, true]);
        let before = remove;
        assert!(matches!(
            identify_removed(&q, &[usize::MAX, 1], &BTreeMap::new(), &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::RowOverflow {
                    kind: "matched reactant row",
                    ..
                }
            ))
        ));
        assert_eq!(remove, before);
    }
}

#[cfg(test)]
mod complete_traverse_removed_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondId, CoordinateBlock, MoleculeProperties,
        PropertyValue, QueryAtom,
    };
    use cosmolkit_types::Element;
    fn graph(maps: &[Option<i32>]) -> QueryGraph {
        QueryGraph::from_parts(
            maps.iter()
                .enumerate()
                .map(|(i, map)| {
                    let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
                    if let Some(map) = map {
                        a.set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                            .unwrap();
                    }
                    a
                })
                .collect(),
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn topology(count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            atoms: (0..count)
                .map(|i| Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C)))
                .collect(),
            adjacency: AdjacencyList::from_topology(count, &bonds),
            bonds,
            ..Default::default()
        }
    }
    fn input<'a>(
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
    fn run(
        t: &TopologyBlock,
        q: &QueryGraph,
        pairs: &[(i32, i32)],
        remove: &mut [bool],
    ) -> Result<(), ReactionApplyError> {
        traverse_removed_source(
            input(
                t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            q,
            pairs.iter().copied().map(Ok),
            remove,
        )
    }
    #[test]
    fn no_matches_marks_entire_mask_even_when_source_graph_shape_is_unreached() {
        let mut remove = [false, true, false, true];
        run(&TopologyBlock::default(), &graph(&[]), &[], &mut remove).unwrap();
        assert_eq!(remove, [true; 4]);
    }
    #[test]
    fn retained_mapped_component_cycle_and_duplicate_paths_preserve_removed_barrier_and_drop_disconnected_atoms()
     {
        let t = topology(6, &[(0, 1), (0, 2), (1, 3), (2, 3), (3, 0), (0, 5)]);
        let mut remove = [false, false, false, false, false, true];
        run(&t, &graph(&[Some(11)]), &[(0, 0)], &mut remove).unwrap();
        assert_eq!(remove, [false, false, false, false, true, true]);
    }
    #[test]
    fn signed_sparse_reordered_query_pairs_control_which_matched_component_can_seed() {
        let t = topology(4, &[(0, 1), (2, 3)]);
        let q = graph(&[None, None, Some(-7)]);
        let mut remove = [false; 4];
        run(&t, &q, &[(2, 2), (0, 0)], &mut remove).unwrap();
        assert_eq!(remove, [false, true, false, false]);
    }
    #[test]
    fn actual_source_atom_ids_drive_reset_and_neighbor_identity_without_eager_id_canonicalization()
    {
        let q = graph(&[Some(1)]);
        let mut t = topology(2, &[]);
        t.atoms[1] = t.atoms[1].clone().with_id(AtomId::new(0));
        let mut remove = [false; 2];
        run(&t, &q, &[(0, 1)], &mut remove).unwrap();
        assert_eq!(remove, [false, false]);
        let mut t = topology(3, &[(0, 1)]);
        t.atoms[1] = t.atoms[1].clone().with_id(AtomId::new(2));
        let mut remove = [false; 3];
        run(&t, &q, &[(0, 0)], &mut remove).unwrap();
        assert_eq!(remove, [false, true, false]);
    }
    #[test]
    fn all_match_resets_precede_any_template_lookup_and_fail_without_committing_original_mask() {
        let mut remove = [false];
        assert!(matches!(
            run(
                &TopologyBlock::default(),
                &graph(&[]),
                &[(-1, 0), (0, 99)],
                &mut remove
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "matched reset bit index out of range",
                    reactant_atom: Some(99),
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false]);
    }
    #[test]
    fn absent_map_and_already_removed_seeds_skip_missing_input_pointers_and_csr_rows() {
        for (map, original) in [(None, false), (Some(1), true)] {
            let mut remove = [original];
            run(
                &TopologyBlock::default(),
                &graph(&[map]),
                &[(0, 0)],
                &mut remove,
            )
            .unwrap();
            assert_eq!(remove, [original]);
        }
    }
    #[test]
    fn reached_missing_csr_row_is_error_instead_of_empty_neighbor_fallback() {
        let mut t = topology(1, &[]);
        t.adjacency = AdjacencyList::default();
        let mut remove = [false];
        assert!(matches!(
            run(&t, &graph(&[Some(1)]), &[(0, 0)], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "reactant adjacency row missing",
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false]);
    }
    #[test]
    fn reached_missing_seed_pointer_preserves_mask_before_final_or() {
        let mut remove = [false];
        assert!(matches!(
            run(
                &TopologyBlock::default(),
                &graph(&[Some(1)]),
                &[(0, 0)],
                &mut remove
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "seed reactant atom row missing",
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false]);
    }
    #[test]
    fn reached_missing_neighbor_pointer_preserves_mask_before_final_or() {
        let mut t = topology(2, &[(0, 1)]);
        t.atoms.truncate(1);
        let mut remove = [false; 2];
        assert!(matches!(
            run(&t, &graph(&[Some(1)]), &[(0, 0)], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "neighbor reactant atom row missing",
                    reactant_atom: Some(1),
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false; 2]);
    }
    #[test]
    fn fifo_physical_neighbor_order_selects_first_reached_branch_error_without_mask_commit() {
        let mut t = topology(5, &[(0, 1), (0, 2), (1, 3), (2, 4)]);
        t.atoms[3] = t.atoms[3].clone().with_id(AtomId::new(99));
        t.atoms[4] = t.atoms[4].clone().with_id(AtomId::new(98));
        let mut remove = [false; 5];
        assert!(matches!(
            run(&t, &graph(&[Some(1)]), &[(0, 0)], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "neighbor bit index out of range",
                    reactant_atom: Some(99),
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false; 5]);
    }
    #[test]
    fn later_template_property_error_discards_all_temporary_traversal_changes() {
        let t = topology(4, &[(0, 1)]);
        let mut q = graph(&[Some(1), Some(2)]);
        q.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::String("bad".into()))
            .unwrap();
        let mut remove = [false; 4];
        assert!(matches!(
            run(&t, &q, &[(0, 0), (1, 2)], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::TemplateProperty(_)
            ))
        ));
        assert_eq!(remove, [false; 4]);
    }
    #[test]
    fn canonical_adapter_preserves_source_signed_bit_width_and_reports_larger_usize_overflow() {
        let t = topology(0, &[]);
        let q = graph(&[Some(1)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let mut remove = [false];
        assert!(matches!(
            traverse_removed(input(&t, &c, &p), &q, &[u32::MAX as usize], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "matched reset bit index out of range",
                    reactant_atom: Some(usize::MAX),
                    ..
                }
            ))
        ));
        assert!(matches!(
            traverse_removed(input(&t, &c, &p), &q, &[usize::MAX], &mut remove),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::RowOverflow {
                    kind: "matched reactant row",
                    ..
                }
            ))
        ));
        assert_eq!(remove, [false]);
    }
}

#[cfg(test)]
mod complete_update_atoms_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomQueryPredicate, AtomSpec, Bond, BondId, CoordinateBlock,
        MoleculeProperties, PropertyValue, QueryAtom, QueryBond, QueryNode,
    };
    use cosmolkit_types::Element;
    fn topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            atoms: elements
                .iter()
                .enumerate()
                .map(|(i, &e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
                .collect(),
            adjacency: AdjacencyList::from_topology(elements.len(), &bonds),
            bonds,
            ..Default::default()
        }
    }
    fn graph(elements: &[Element], maps: &[Option<i32>], edges: &[(usize, usize)]) -> QueryGraph {
        let atoms = elements
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(e));
                if let Some(Some(map)) = maps.get(i) {
                    a.set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                        .unwrap();
                }
                a
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b))| {
                QueryBond::new(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![]).unwrap()
    }
    fn input<'a>(
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
    fn run(
        edit: &mut TopologyBatchEdit,
        t: &TopologyBlock,
        r: &QueryGraph,
        p: &QueryGraph,
        pm: &BTreeMap<u32, usize>,
        kept: &BTreeMap<u32, usize>,
        pairs: &[(i32, i32)],
    ) -> Result<bool, ReactionApplyError> {
        update_atoms_source(
            edit,
            input(
                t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            r,
            p,
            pm,
            kept,
            pairs.len(),
            |index| {
                Ok(*pairs.get(index).ok_or_else(|| {
                    invariant(
                        "updateAtomsModifiedByReaction",
                        "template match index out of range",
                        None,
                        Some(index),
                        None,
                    )
                })?)
            },
        )
    }
    fn one(element: Element) -> QueryGraph {
        graph(&[element], &[Some(1)], &[])
    }
    fn raw(edit: TopologyBatchEdit) -> TopologyBlock {
        edit.into_source_batch_parts().0
    }
    #[test]
    fn empty_preserved_map_never_reads_templates_matches_or_working_atoms() {
        let t = topology(&[], &[]);
        let mut edit = t.clone().into_batch_edit().unwrap();
        let q = one(Element::C);
        assert!(
            !update_atoms_source(
                &mut edit,
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                ),
                &q,
                &q,
                &BTreeMap::new(),
                &BTreeMap::new(),
                usize::MAX,
                |_| panic!("unreached")
            )
            .unwrap()
        );
        assert!(raw(edit).atoms.is_empty());
    }
    #[test]
    fn checked_lookup_order_is_reactant_then_product_map_then_product_atom_then_match_then_actual_atom()
     {
        let r = one(Element::C);
        let p = one(Element::N);
        let t = topology(&[Element::C], &[]);
        for case in 0..5 {
            let mut edit = t.clone().into_batch_edit().unwrap();
            let kept = BTreeMap::from([(1, if case == 0 { 99 } else { 0 })]);
            let pm = if case == 1 {
                BTreeMap::new()
            } else {
                BTreeMap::from([(1, if case == 2 { 99 } else { 0 })])
            };
            let mut reads = 0;
            let result = update_atoms_source(
                &mut edit,
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default(),
                ),
                &r,
                &p,
                &pm,
                &kept,
                1,
                |_| {
                    reads += 1;
                    if case == 3 {
                        Err(invariant(
                            "updateAtomsModifiedByReaction",
                            "template match index out of range",
                            None,
                            Some(0),
                            None,
                        )
                        .into())
                    } else {
                        Ok((999, if case == 4 { -1 } else { 0 }))
                    }
                },
            );
            assert!(result.is_err());
            assert_eq!(reads, if case < 3 { 0 } else { 1 });
            assert_eq!(raw(edit).atoms[0].element(), Element::C);
        }
    }
    #[test]
    fn element_change_decision_compares_template_atoms_and_reports_assignment_even_if_actual_element_same()
     {
        for (actual, product, expected, changed) in [
            (Element::N, Element::N, Element::N, true),
            (Element::O, Element::C, Element::O, false),
        ] {
            let t = topology(&[actual], &[]);
            let mut edit = t.clone().into_batch_edit().unwrap();
            let result = run(
                &mut edit,
                &t,
                &one(Element::C),
                &one(product),
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[(123, 0)],
            )
            .unwrap();
            assert_eq!(result, changed);
            assert_eq!(raw(edit).atoms[0].element(), expected);
        }
    }
    #[test]
    fn dummy_query_preserves_actual_element_but_ordinary_dummy_carrier_sets_zero() {
        let r = one(Element::C);
        for query in [true, false] {
            let mut p = one(Element::DUMMY);
            if !query {
                p.atoms_mut()[0] = QueryAtom::from_carrier_parts(
                    Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::DUMMY)),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0)),
                );
            }
            let t = topology(&[Element::C], &[]);
            let mut edit = t.clone().into_batch_edit().unwrap();
            assert_eq!(
                run(
                    &mut edit,
                    &t,
                    &r,
                    &p,
                    &BTreeMap::from([(1, 0)]),
                    &BTreeMap::from([(1, 0)]),
                    &[(0, 0)]
                )
                .unwrap(),
                !query
            );
            assert_eq!(
                raw(edit).atoms[0].element(),
                if query { Element::C } else { Element::DUMMY }
            );
        }
    }
    #[test]
    fn source_template_properties_and_errors_preserve_scalar_write_order_before_chirality() {
        for failure in [0, 1] {
            let mut p = one(Element::N);
            for (key, value) in [
                ("_QueryFormalCharge", PropertyValue::Int(-2)),
                ("_QueryHCount", PropertyValue::UInt(3)),
                (
                    "_QueryMass",
                    if failure == 0 {
                        PropertyValue::String("bad".into())
                    } else {
                        PropertyValue::UInt(12)
                    },
                ),
                ("_QueryIsotope", PropertyValue::UInt(13)),
            ] {
                p.atoms_mut()[0].set_prop(key, value).unwrap();
            }
            if failure == 1 {
                p.atoms_mut()[0]
                    .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
                    .unwrap();
            }
            let mut t = topology(&[Element::C], &[]);
            t.atoms[0].set_isotope(Some(1));
            let mut edit = t.clone().into_batch_edit().unwrap();
            assert!(
                run(
                    &mut edit,
                    &t,
                    &one(Element::C),
                    &p,
                    &BTreeMap::from([(1, 0)]),
                    &BTreeMap::from([(1, 0)]),
                    &[(0, 0)]
                )
                .is_err()
            );
            let result = raw(edit);
            let a = &result.atoms[0];
            assert_eq!(
                (
                    a.element(),
                    a.formal_charge(),
                    a.explicit_hydrogens(),
                    a.no_implicit(),
                    a.isotope()
                ),
                (
                    Element::N,
                    -2,
                    3,
                    true,
                    Some(if failure == 0 { 1 } else { 13 })
                )
            );
        }
    }
    #[test]
    fn inversion_flags_cover_all_chiral_tags_with_source_modified_bool_independent_of_equality() {
        let tags = [
            ChiralTag::Unspecified,
            ChiralTag::Other,
            ChiralTag::TetrahedralCw,
            ChiralTag::TetrahedralCcw,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ];
        for tag in tags {
            for flag in [-1, 0, 1, 2, 3, 4, 5] {
                let mut t = topology(&[Element::C], &[]);
                t.atoms[0].set_chiral_tag(tag);
                let mut p = one(Element::C);
                p.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
                p.atoms_mut()[0]
                    .set_prop("molInversionFlag", PropertyValue::Int(flag))
                    .unwrap();
                let mut edit = t.clone().into_batch_edit().unwrap();
                let changed = run(
                    &mut edit,
                    &t,
                    &one(Element::C),
                    &p,
                    &BTreeMap::from([(1, 0)]),
                    &BTreeMap::from([(1, 0)]),
                    &[(0, 0)],
                )
                .unwrap();
                let expected = match flag {
                    1 => match tag {
                        ChiralTag::TetrahedralCw => ChiralTag::TetrahedralCcw,
                        ChiralTag::TetrahedralCcw => ChiralTag::TetrahedralCw,
                        _ => tag,
                    },
                    3 => ChiralTag::Unspecified,
                    4 => ChiralTag::TetrahedralCw,
                    _ => tag,
                };
                assert_eq!(raw(edit).atoms[0].chiral_tag(), expected);
                assert_eq!(
                    changed,
                    flag == 3
                        || flag == 4
                        || (flag == 1 && !matches!(tag, ChiralTag::Other | ChiralTag::Unspecified))
                );
            }
        }
    }
    #[test]
    fn non_tetrahedral_permutation_inversion_uses_sole_core_tables() {
        for tag in [
            ChiralTag::Tetrahedral,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            let mut t = topology(&[Element::C], &[]);
            t.atoms[0] = Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_chiral_tag(tag)
                    .with_chiral_permutation(1),
            );
            let mut p = one(Element::C);
            p.atoms_mut()[0]
                .set_prop("molInversionFlag", PropertyValue::Int(1))
                .unwrap();
            let mut edit = t.clone().into_batch_edit().unwrap();
            assert!(
                run(
                    &mut edit,
                    &t,
                    &one(Element::C),
                    &p,
                    &BTreeMap::from([(1, 0)]),
                    &BTreeMap::from([(1, 0)]),
                    &[(0, 0)]
                )
                .unwrap()
            );
            assert_eq!(raw(edit).atoms[0].chiral_permutation(), Some(2));
        }
    }
    fn star_p(maps: &[Option<i32>]) -> QueryGraph {
        let mut p = graph(&[Element::C; 3], maps, &[(0, 1), (0, 2)]);
        p.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        p.atoms_mut()[0]
            .set_prop("molInversionFlag", PropertyValue::Int(4))
            .unwrap();
        p
    }
    #[test]
    fn creation_uses_signed_pair_first_for_neighbor_template_lookup_and_swaps_odd_permutation() {
        let t = topology(&[Element::C; 3], &[(0, 1), (0, 2)]);
        let r = graph(
            &[Element::C; 3],
            &[Some(1), Some(-7), Some(8)],
            &[(0, 1), (0, 2)],
        );
        let p = star_p(&[Some(1), Some(-7), Some(8)]);
        let mut edit = t.clone().into_batch_edit().unwrap();
        assert!(
            run(
                &mut edit,
                &t,
                &r,
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[(0, 0), (2, 1), (1, 2)]
            )
            .unwrap()
        );
        assert_eq!(raw(edit).atoms[0].chiral_tag(), ChiralTag::TetrahedralCcw);
    }
    #[test]
    fn creation_find_if_keeps_first_duplicate_and_never_dereferences_later_invalid_query_index() {
        let t = topology(&[Element::C; 3], &[(0, 1), (0, 2)]);
        let r = graph(
            &[Element::C; 3],
            &[Some(1), Some(7), Some(8)],
            &[(0, 1), (0, 2)],
        );
        let p = star_p(&[Some(1), Some(7), Some(8)]);
        let mut edit = t.clone().into_batch_edit().unwrap();
        run(
            &mut edit,
            &t,
            &r,
            &p,
            &BTreeMap::from([(1, 0)]),
            &BTreeMap::from([(1, 0)]),
            &[(0, 0), (1, 1), (-1, 1), (2, 2)],
        )
        .unwrap();
        assert_eq!(raw(edit).atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
    }
    #[test]
    fn creation_skips_unmapped_neighbors_and_unequal_lengths_but_propagates_equal_length_permutation_error_after_set_tag()
     {
        let t = topology(&[Element::C; 3], &[(0, 1), (0, 2)]);
        let r = graph(
            &[Element::C; 3],
            &[Some(1), Some(7), Some(8)],
            &[(0, 1), (0, 2)],
        );
        for last in [None, Some(9)] {
            let p = star_p(&[Some(1), Some(7), last]);
            let mut edit = t.clone().into_batch_edit().unwrap();
            let result = run(
                &mut edit,
                &t,
                &r,
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[(0, 0), (1, 1), (2, 2)],
            );
            if last.is_none() {
                assert!(result.unwrap());
            } else {
                assert!(matches!(
                    result,
                    Err(ReactionApplyError::Product(
                        crate::ReactionProductError::StereoOrder(_)
                    ))
                ));
            }
            assert_eq!(raw(edit).atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        }
    }
    #[test]
    fn creation_checks_product_neighbors_then_actual_working_atom_csr_after_assigning_chiral_tag() {
        let t = topology(&[Element::C], &[]);
        let r = one(Element::C);
        for bad_product in [true, false] {
            let mut p = one(Element::C);
            p.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
            p.atoms_mut()[0]
                .set_prop("molInversionFlag", PropertyValue::Int(4))
                .unwrap();
            let mut edit = t.clone().into_batch_edit().unwrap();
            if bad_product {
                p.atoms_mut()[0] = p.atoms()[0].clone().with_id(AtomId::new(99));
            } else {
                let a = edit.atom_mut(AtomId::new(0)).unwrap();
                *a = a.clone().with_id(AtomId::new(99));
            }
            let result = run(
                &mut edit,
                &t,
                &r,
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[(0, 0)],
            );
            assert!(
                matches!(result,Err(ReactionApplyError::Product(crate::ReactionProductError::Invariant{detail,..})) if detail==if bad_product {"product template adjacency row missing"}else{"reactant adjacency row missing"})
            );
            assert_eq!(raw(edit).atoms[0].chiral_tag(), ChiralTag::TetrahedralCw);
        }
    }
    #[test]
    fn unsigned_map_order_preserves_earlier_atom_write_before_later_missing_map_key() {
        let t = topology(&[Element::C; 2], &[]);
        let r = graph(&[Element::C; 2], &[], &[]);
        let p = graph(&[Element::C, Element::N], &[], &[]);
        let mut edit = t.clone().into_batch_edit().unwrap();
        assert!(matches!(
            run(
                &mut edit,
                &t,
                &r,
                &p,
                &BTreeMap::from([(1, 1)]),
                &BTreeMap::from([(u32::MAX, 0), (1, 1)]),
                &[(0, 0), (1, 1)]
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "product atom map key missing",
                    ..
                }
            ))
        ));
        let result = raw(edit);
        assert_eq!(result.atoms[0].element(), Element::C);
        assert_eq!(result.atoms[1].element(), Element::N);
    }
    #[test]
    fn canonical_adapter_leaves_unused_match_rows_unconverted_and_inversion_error_retains_prior_atomic_change()
     {
        let t = topology(&[Element::C], &[]);
        let r = one(Element::C);
        let p = one(Element::C);
        let mut edit = t.clone().into_batch_edit().unwrap();
        assert!(
            !update_atoms(
                &mut edit,
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                ),
                &r,
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[0, usize::MAX]
            )
            .unwrap()
        );
        let mut p = one(Element::N);
        p.atoms_mut()[0]
            .set_prop("molInversionFlag", PropertyValue::Int(1))
            .unwrap();
        let a = edit.atom_mut(AtomId::new(0)).unwrap();
        a.set_chiral_tag(ChiralTag::Tetrahedral);
        a.set_prop("_chiralPermutation", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            run(
                &mut edit,
                &t,
                &r,
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[(0, 0)]
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::StereoOrder(_)
            ))
        ));
        let result = raw(edit);
        assert_eq!(result.atoms[0].element(), Element::N);
        assert_eq!(result.atoms[0].chiral_tag(), ChiralTag::Tetrahedral);
    }
}

#[cfg(test)]
mod complete_update_bonds_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondId, BondQueryPredicate, PropertyValue, QueryAtom,
        QueryBond, QueryNode, SourceAtomValenceFacts,
    };
    use cosmolkit_types::{BondDirection, BondStereo, Element};
    fn graph(maps: &[i32], edges: &[(usize, usize, BondOrder)], query: bool) -> QueryGraph {
        let atoms = maps
            .iter()
            .enumerate()
            .map(|(i, &map)| {
                let mut a = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
                a.set_prop("molAtomMapNumber", PropertyValue::Int(map))
                    .unwrap();
                a
            })
            .collect();
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b, o))| {
                let bond = Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o)
                        .with_aromatic(true)
                        .with_conjugated(true)
                        .with_direction(BondDirection::BeginWedge)
                        .with_stereo(BondStereo::Any)
                        .with_prop("template-only", PropertyValue::Int(1))
                        .unwrap(),
                );
                let predicate =
                    QueryNode::Not(Box::new(QueryNode::predicate(BondQueryPredicate::Any)));
                if query {
                    QueryBond::from_parts(bond, predicate)
                } else {
                    QueryBond::from_carrier_parts(bond, predicate)
                }
            })
            .collect();
        QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![]).unwrap()
    }
    fn topology(n: usize, edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            atoms: (0..n)
                .map(|i| {
                    let mut a = Atom::from_spec(AtomId::new(i), AtomSpec::new(Element::C));
                    a.set_source_valence_facts(SourceAtomValenceFacts {
                        explicit_valence: 3,
                        implicit_valence: 1,
                    });
                    a
                })
                .collect(),
            adjacency: AdjacencyList::from_topology(n, &bonds),
            bonds,
            ..Default::default()
        }
    }
    fn maps() -> BTreeMap<u32, usize> {
        BTreeMap::from([(1, 0), (2, 1)])
    }
    fn run(
        edit: &mut TopologyBatchEdit,
        r: &QueryGraph,
        p: &QueryGraph,
        pm: &BTreeMap<u32, usize>,
        kept: &BTreeMap<u32, usize>,
        matches: &[i32],
    ) -> Result<(bool, bool), ReactionApplyError> {
        update_bonds_source(edit, r, p, pm, kept, |i| {
            Ok(*matches.get(i).ok_or_else(|| {
                invariant(
                    "updateBondsModifiedByReaction",
                    "template match index out of range",
                    None,
                    Some(i),
                    None,
                )
            })?)
        })
    }
    #[test]
    fn empty_map_reads_no_templates_matches_or_editor() {
        let mut e = topology(0, &[]).into_batch_edit().unwrap();
        let q = graph(&[], &[], false);
        assert_eq!(
            update_bonds_source(
                &mut e,
                &q,
                &q,
                &BTreeMap::new(),
                &BTreeMap::new(),
                |_| panic!("unreached")
            )
            .unwrap(),
            (false, false)
        );
    }
    #[test]
    fn lookup_precedence_and_signed_target_failures_do_not_panic_or_write() {
        let q = graph(&[1], &[], false);
        for case in 0..5 {
            let mut e = topology(1, &[]).into_batch_edit().unwrap();
            let kept = BTreeMap::from([(1, if case == 0 { 99 } else { 0 })]);
            let pm = if case == 1 {
                BTreeMap::new()
            } else {
                BTreeMap::from([(1, if case == 2 { 99 } else { 0 })])
            };
            let mut reads = 0;
            assert!(
                update_bonds_source(&mut e, &q, &q, &pm, &kept, |_| {
                    reads += 1;
                    if case == 3 {
                        Err(invariant("test", "match failure", None, None, None).into())
                    } else {
                        Ok(-1)
                    }
                })
                .is_err()
            );
            assert_eq!(reads, usize::from(case >= 3));
            assert!(e.into_source_batch_parts().0.bonds.is_empty());
        }
    }
    #[test]
    fn existing_template_bond_compares_templates_and_modified_bool_even_if_actual_order_equal() {
        let r = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
        let p = graph(&[1, 2], &[(0, 1, BondOrder::Double)], false);
        let mut e = topology(2, &[(0, 1, BondOrder::Double)])
            .into_batch_edit()
            .unwrap();
        assert_eq!(
            run(&mut e, &r, &p, &maps(), &maps(), &[0, 1]).unwrap(),
            (true, false)
        );
        assert_eq!(
            e.into_source_batch_parts().0.bonds[0].order(),
            BondOrder::Double
        );
    }
    #[test]
    fn same_template_order_and_unspecified_product_never_read_actual_bond_endpoints() {
        for o in [BondOrder::Single, BondOrder::Unspecified] {
            let r = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
            let p = graph(&[1, 2], &[(0, 1, o)], false);
            let mut e = topology(1, &[]).into_batch_edit().unwrap();
            assert_eq!(
                run(
                    &mut e,
                    &r,
                    &p,
                    &maps(),
                    &BTreeMap::from([(1, 0), (2, 1)]),
                    &[0, 0]
                )
                .unwrap(),
                (false, false)
            );
            assert!(e.into_source_batch_parts().0.bonds.is_empty());
        }
    }
    #[test]
    fn differing_template_bond_requires_actual_edge_and_retains_prior_updates_on_later_error() {
        let r = graph(
            &[1, 2, 3],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
            false,
        );
        let p = graph(
            &[1, 2, 3],
            &[(0, 1, BondOrder::Double), (0, 2, BondOrder::Double)],
            false,
        );
        let mut e = topology(3, &[(0, 1, BondOrder::Single)])
            .into_batch_edit()
            .unwrap();
        assert!(matches!(
            run(
                &mut e,
                &r,
                &p,
                &BTreeMap::from([(1, 0), (2, 1), (3, 2)]),
                &BTreeMap::from([(1, 0), (2, 1), (3, 2)]),
                &[0, 1, 2]
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "missing bond between known neighbors in reactant",
                    ..
                }
            ))
        ));
        assert_eq!(
            e.into_source_batch_parts().0.bonds[0].order(),
            BondOrder::Double
        );
    }
    #[test]
    fn no_template_edge_updates_existing_actual_edge_including_unspecified_order_without_copying_query()
     {
        for o in [BondOrder::Unspecified, BondOrder::Double] {
            let r = graph(&[1, 2], &[], false);
            let p = graph(&[1, 2], &[(0, 1, o)], true);
            let mut e = topology(2, &[(0, 1, BondOrder::Single)])
                .into_batch_edit()
                .unwrap();
            assert_eq!(
                run(&mut e, &r, &p, &maps(), &maps(), &[0, 1]).unwrap(),
                (true, false)
            );
            let t = e.into_source_batch_parts().0;
            assert_eq!(t.bonds[0].order(), o);
            assert!(t.bonds[0].query().is_none());
            assert!(!t.atoms[0].is_aromatic());
            assert_eq!(t.atoms[0].source_valence_facts().explicit_valence, 3);
        }
    }
    #[test]
    fn new_ordinary_bond_uses_product_endpoint_order_and_source_cache_aromatic_and_pending_masks() {
        let r = graph(&[1, 2], &[], false);
        let p = graph(&[1, 2], &[(1, 0, BondOrder::Aromatic)], false);
        let mut e = topology(3, &[]).into_batch_edit().unwrap();
        e.remove_atom(AtomId::new(1)).unwrap();
        assert_eq!(
            run(&mut e, &r, &p, &maps(), &maps(), &[1, 2]).unwrap(),
            (true, false)
        );
        let (t, _, mask) = e.into_source_batch_parts();
        assert_eq!(
            (t.bonds[0].begin(), t.bonds[0].end()),
            (AtomId::new(2), AtomId::new(1))
        );
        assert!(t.bonds[0].is_aromatic());
        assert!(t.bonds[0].query().is_none());
        assert!(t.atoms[1].is_aromatic() && t.atoms[2].is_aromatic());
        assert_eq!(
            t.atoms[1].source_valence_facts(),
            SourceAtomValenceFacts::UNINITIALIZED
        );
        assert_eq!(mask.unwrap(), vec![true]);
        assert_eq!(t.bonds[0].direction(), BondDirection::None);
        assert_eq!(t.bonds[0].stereo(), BondStereo::None);
        assert!(t.bonds[0].prop("template-only").is_none());
    }
    #[test]
    fn new_query_bond_deep_copies_only_predicate_preserves_cache_and_does_not_inherit_atom_removal()
    {
        let r = graph(&[1, 2], &[], false);
        let p = graph(&[1, 2], &[(1, 0, BondOrder::Aromatic)], true);
        let mut e = topology(2, &[]).into_batch_edit().unwrap();
        e.remove_atom(AtomId::new(0)).unwrap();
        assert_eq!(
            run(&mut e, &r, &p, &maps(), &maps(), &[0, 1]).unwrap(),
            (true, false)
        );
        let (t, _, mask) = e.into_source_batch_parts();
        let b = &t.bonds[0];
        assert_eq!((b.begin(), b.end()), (AtomId::new(1), AtomId::new(0)));
        assert_eq!(b.query(), Some(p.bonds()[0].predicate()));
        assert!(!std::ptr::eq(b.query().unwrap(), p.bonds()[0].predicate()));
        assert!(!b.is_aromatic() && !b.is_conjugated());
        assert_eq!(b.direction(), BondDirection::None);
        assert_eq!(b.stereo(), BondStereo::None);
        assert!(b.prop("template-only").is_none());
        assert_eq!(t.atoms[0].source_valence_facts().explicit_valence, 3);
        assert!(!t.atoms[0].is_aromatic());
        assert_eq!(mask.unwrap(), vec![false]);
    }
    #[test]
    fn removal_keeps_pending_edge_visible_duplicate_visits_and_modified_bool_for_absent_edge() {
        let r = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
        let p = graph(&[1, 2], &[], false);
        for exists in [false, true] {
            let t = topology(
                2,
                if exists {
                    &[(0, 1, BondOrder::Single)]
                } else {
                    &[]
                },
            );
            let mut e = t.into_batch_edit().unwrap();
            assert_eq!(
                run(&mut e, &r, &p, &maps(), &maps(), &[0, 1]).unwrap(),
                (true, exists)
            );
            assert_eq!(
                e.bond_between_atoms(AtomId::new(0), AtomId::new(1))
                    .unwrap()
                    .is_some(),
                exists
            );
            let (t, _, mask) = e.into_source_batch_parts();
            assert_eq!(t.bonds.len(), usize::from(exists));
            assert_eq!(mask.unwrap(), if exists { vec![true] } else { vec![] });
        }
    }
    #[test]
    fn unmapped_or_unpreserved_neighbors_skip_later_match_and_edge_reads_but_property_errors_propagate()
     {
        for map in [0, 3] {
            let r = graph(&[1], &[], false);
            let p = graph(&[1, map], &[(0, 1, BondOrder::Double)], false);
            let mut e = topology(1, &[]).into_batch_edit().unwrap();
            assert_eq!(
                run(
                    &mut e,
                    &r,
                    &p,
                    &BTreeMap::from([(1, 0)]),
                    &BTreeMap::from([(1, 0)]),
                    &[0]
                )
                .unwrap(),
                (false, false)
            );
        }
        let r = graph(&[1, 2], &[], false);
        let mut p = graph(&[1, 2], &[(0, 1, BondOrder::Double)], false);
        p.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::String("bad".into()))
            .unwrap();
        let mut e = topology(2, &[]).into_batch_edit().unwrap();
        assert!(run(&mut e, &r, &p, &maps(), &maps(), &[0, 1]).is_err());
        assert!(e.into_source_batch_parts().0.bonds.is_empty());
    }
    #[test]
    fn actual_atom_id_controls_removal_endpoint_and_source_template_ids_control_neighbors() {
        let r = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
        let p = graph(&[1, 2], &[], false);
        let mut e = topology(3, &[(1, 2, BondOrder::Single)])
            .into_batch_edit()
            .unwrap();
        let a = e.atom_mut(AtomId::new(0)).unwrap();
        *a = a.clone().with_id(AtomId::new(2));
        assert_eq!(
            run(&mut e, &r, &p, &maps(), &BTreeMap::from([(1, 0)]), &[0, 1]).unwrap(),
            (true, true)
        );
        assert_eq!(e.into_source_batch_parts().2.unwrap(), vec![true]);
        let mut p = graph(&[1], &[], false);
        p.atoms_mut()[0] = p.atoms()[0].clone().with_id(AtomId::new(99));
        let mut e = topology(1, &[]).into_batch_edit().unwrap();
        assert!(matches!(
            run(
                &mut e,
                &graph(&[1], &[], false),
                &p,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[0]
            ),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::Invariant {
                    detail: "template adjacency row missing",
                    ..
                }
            ))
        ));
    }
    #[test]
    fn canonical_match_adapter_converts_only_reached_rows_and_unsigned_map_order_keeps_write_prefix()
     {
        let q = graph(&[1], &[], false);
        let mut e = topology(1, &[]).into_batch_edit().unwrap();
        assert_eq!(
            update_bonds(
                &mut e,
                &q,
                &q,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[0, usize::MAX]
            )
            .unwrap(),
            (false, false)
        );
        let r = graph(&[1, 2], &[], false);
        let p = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
        let mut e = topology(2, &[]).into_batch_edit().unwrap();
        assert!(
            run(
                &mut e,
                &r,
                &p,
                &maps(),
                &BTreeMap::from([(u32::MAX, 99), (1, 0), (2, 1)]),
                &[0, 1]
            )
            .is_err()
        );
        assert_eq!(e.into_source_batch_parts().0.bonds.len(), 1);
    }
    #[test]
    fn actual_id_conversion_is_reached_only_for_removal_and_query_edge_adapter_keeps_range_check_order()
     {
        let q = graph(&[1], &[], false);
        let mut e = topology(1, &[]).into_batch_edit().unwrap();
        let a = e.atom_mut(AtomId::new(0)).unwrap();
        *a = a.clone().with_id(AtomId::new(usize::MAX));
        assert_eq!(
            run(
                &mut e,
                &q,
                &q,
                &BTreeMap::from([(1, 0)]),
                &BTreeMap::from([(1, 0)]),
                &[0]
            )
            .unwrap(),
            (false, false)
        );
        for (begin, end, expected) in [(99, 100, 99), (0, 100, 100)] {
            assert!(
                matches!(apply_template_edge(&q,begin,end),Err(ReactionApplyError::Edit(cosmolkit_model::TopologyEditError::AtomOutOfRange{atom,..})) if atom.index()==expected)
            );
        }
        let q = graph(&[1, 2], &[(0, 1, BondOrder::Single)], false);
        assert_eq!(
            apply_template_edge(&q, 0, 1).unwrap().unwrap().id(),
            BondId::new(0)
        );
        assert_eq!(
            apply_template_edge(&q, 1, 0).unwrap().unwrap().id(),
            BondId::new(0)
        );
    }
}

#[cfg(test)]
mod complete_apply_reaction_source_tests {
    use super::*;
    use cosmolkit_model::{
        AdjacencyList, Atom, AtomSpec, Bond, BondId, Conformer3D, CoordinateBlock,
        MoleculeProperties, PropertyValue, SourceAtomValenceFacts,
    };
    use cosmolkit_types::Element;
    fn topology(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> TopologyBlock {
        let bonds = edges
            .iter()
            .enumerate()
            .map(|(i, &(a, b, o))| {
                Bond::from_spec(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            atoms: elements
                .iter()
                .enumerate()
                .map(|(i, &e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
                .collect(),
            adjacency: AdjacencyList::from_topology(elements.len(), &bonds),
            bonds,
            ..Default::default()
        }
    }
    fn reaction(text: &str) -> Reaction {
        let mut r = crate::parse_smirks(text).unwrap();
        r.needs_init = false;
        r
    }
    fn input<'a>(
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
    fn run(
        r: &Reaction,
        t: &TopologyBlock,
        remove: bool,
    ) -> Result<ReactionApplyChanges, ReactionApplyError> {
        apply_reaction_source(
            r,
            input(
                t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            &ReactionApplyParams {
                remove_unmatched_atoms: remove,
            },
        )
    }
    #[test]
    fn arity_precedes_initialization_and_initialization_precedes_template_properties() {
        let t = topology(&[], &[]);
        assert!(matches!(
            run(&Reaction::new(), &t, true),
            Err(ReactionApplyError::ApplicabilityArity {
                reactants: 0,
                products: 0
            })
        ));
        let mut r = reaction("[C:1]>>[C:1]");
        r.needs_init = true;
        r.products[0].atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::String("bad".into()))
            .unwrap();
        assert!(matches!(
            run(&r, &t, true),
            Err(ReactionApplyError::Matching(
                crate::ReactionRunError::NeedsInitialization
            ))
        ));
    }
    #[test]
    fn approved_public_initialization_uses_private_copy_and_source_keeps_native_gate() {
        let mut r = reaction("[C:1]>>[N:1]");
        r.needs_init = true;
        let t = topology(&[Element::C], &[]);
        assert!(run(&r, &t, true).is_err());
        let out = apply_reaction(
            &r,
            input(
                &t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            &Default::default(),
        )
        .unwrap();
        assert!(out.changed);
        assert_eq!(out.change.unwrap().0.atoms[0].element(), Element::N);
        assert!(!r.is_initialized());
        assert_eq!(t.atoms[0].element(), Element::C);
    }
    #[test]
    fn unmapped_and_new_product_atoms_reject_before_input_matching() {
        for text in ["[C:1]>>C", "[C:1]>>[C:2]"] {
            let r = reaction(text);
            let t = topology(&[Element::N], &[]);
            assert!(
                matches!(run(&r,&t,true),Err(ReactionApplyError::AddsProductAtom{atom}) if atom.index()==0)
            );
        }
    }
    #[test]
    fn product_map_reads_precede_reactant_map_reads_and_duplicate_product_keys_keep_last_atom() {
        let mut r = reaction("[C:1]>>[C:1]");
        // Fresh query atoms have no parser-created authoritative map slot, so
        // these raw wrong-tag values actually reach native getProp<int>.
        r.reactants[0].atoms_mut()[0] =
            cosmolkit_model::QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        r.products[0].atoms_mut()[0] =
            cosmolkit_model::QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C));
        r.reactants[0].atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::String("reactant".into()))
            .unwrap();
        r.products[0].atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::String("product".into()))
            .unwrap();
        assert!(matches!(
            run(&r, &topology(&[Element::C], &[]), true),
            Err(ReactionApplyError::Product(
                crate::ReactionProductError::TemplateProperty(
                    crate::ReactionValidationError::Property {
                        role: ReactionRole::Product,
                        ..
                    }
                )
            ))
        ));
        let r = reaction("[C:1]>>([C:1].[N:1])");
        let out = run(&r, &topology(&[Element::C], &[]), true).unwrap();
        assert!(out.changed);
        assert_eq!(out.change.unwrap().0.atoms[0].element(), Element::N);
    }
    #[test]
    fn no_match_returns_no_detached_copy_and_never_reads_unreached_coordinate_shape() {
        let r = reaction("[C:1]>>[N:1]");
        let t = topology(&[Element::O], &[]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(7, vec![], true)],
            ..Default::default()
        };
        let out = apply_reaction_source(
            &r,
            input(&t, &c, &MoleculeProperties::default()),
            &Default::default(),
        )
        .unwrap();
        assert!(!out.changed && !out.clears_computed_properties && out.change.is_none());
        assert_eq!(c.conformers_3d[0].coordinates().len(), 0);
    }
    #[test]
    fn template_properties_determine_source_modified_bool_and_identity_keeps_no_change() {
        let mut r = reaction("[C:1]>>[C:1]");
        let p = &mut r.products[0].atoms_mut()[0];
        for (k, v) in [
            ("_QueryFormalCharge", PropertyValue::Int(-2)),
            ("_QueryHCount", PropertyValue::UInt(3)),
            ("_QueryIsotope", PropertyValue::UInt(13)),
        ] {
            p.set_prop(k, v).unwrap();
        }
        for already_set in [false, true] {
            let mut t = topology(&[Element::C], &[]);
            t.atoms[0].set_source_valence_facts(SourceAtomValenceFacts {
                explicit_valence: 4,
                implicit_valence: 0,
            });
            t.atoms[0]
                .set_computed_prop("kept", PropertyValue::Int(7))
                .unwrap();
            if already_set {
                t.atoms[0].set_formal_charge(-2);
                t.atoms[0].set_explicit_hydrogens(3);
                t.atoms[0].set_no_implicit(true);
                t.atoms[0].set_isotope(Some(13));
            }
            let before = t.clone();
            let out = run(&r, &t, true).unwrap();
            assert_eq!(out.changed, !already_set);
            assert!(!out.clears_computed_properties);
            if already_set {
                assert!(out.change.is_none());
            } else {
                let (t2, m) = out.change.unwrap();
                assert_eq!(
                    (
                        t2.atoms[0].formal_charge(),
                        t2.atoms[0].explicit_hydrogens(),
                        t2.atoms[0].isotope()
                    ),
                    (-2, 3, Some(13))
                );
                assert_eq!(
                    t2.atoms[0].source_valence_facts(),
                    t.atoms[0].source_valence_facts()
                );
                assert!(t2.atoms[0].prop("kept").is_some());
                assert_eq!(m, TopologyMapping::identity(1, 0));
            }
            assert_eq!(t, before);
        }
    }
    #[test]
    fn actual_atomic_change_returns_source_bool_and_identity_mapping_without_cache_recompute() {
        let r = reaction("[C:1]>>[N:1]");
        let mut t = topology(&[Element::C], &[]);
        t.atoms[0].set_source_valence_facts(SourceAtomValenceFacts {
            explicit_valence: 4,
            implicit_valence: 0,
        });
        let out = run(&r, &t, false).unwrap();
        assert!(out.changed && !out.clears_computed_properties);
        let (t2, m) = out.change.unwrap();
        assert_eq!(t2.atoms[0].element(), Element::N);
        assert_eq!(t2.atoms[0].source_valence_facts().explicit_valence, 4);
        assert_eq!(m, TopologyMapping::identity(1, 0));
    }
    #[test]
    fn source_removal_commits_dense_mapping_clears_computed_props_and_keeps_original_immutable() {
        let r = reaction("[C:1][O:2]>>[C:1]");
        let mut t = topology(
            &[Element::C, Element::O, Element::N],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        for a in &mut t.atoms {
            a.set_computed_prop("computed", PropertyValue::Int(4))
                .unwrap();
        }
        let before = t.clone();
        let mut props = MoleculeProperties::default();
        props
            .set_computed_prop("molecule", PropertyValue::Int(3))
            .unwrap();
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(
                8,
                vec![[0.0; 3], [1.0; 3], [2.0; 3]],
                true,
            )],
            ..Default::default()
        };
        let out = apply_reaction_source(
            &r,
            input(&t, &c, &props),
            &ReactionApplyParams {
                remove_unmatched_atoms: false,
            },
        )
        .unwrap();
        assert!(out.changed && out.clears_computed_properties);
        let (t2, m) = out.change.unwrap();
        assert_eq!(
            t2.atoms.iter().map(Atom::element).collect::<Vec<_>>(),
            [Element::C, Element::N]
        );
        assert_eq!(
            m.atoms.old_to_new,
            [Some(AtomId::new(0)), None, Some(AtomId::new(1))]
        );
        assert_eq!(m.bonds.old_to_new, [None, Some(BondId::new(0))]);
        assert!(t2.atoms[0].prop("computed").is_none());
        assert_eq!(t, before);
        assert!(props.prop("molecule").is_some());
        assert_eq!(c.conformers_3d[0].coordinates().len(), 3);
    }
    #[test]
    fn remove_unmatched_option_keeps_reachable_sidechain_and_controls_disconnected_fragment() {
        let r = reaction("[C:1][O:2]>>[C:1]");
        let t = topology(
            &[Element::C, Element::O, Element::N, Element::F],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        for remove in [false, true] {
            let out = run(&r, &t, remove).unwrap();
            let (t2, m) = out.change.unwrap();
            assert_eq!(t2.atoms.len(), if remove { 2 } else { 3 });
            assert_eq!(m.atoms.old_to_new[3].is_none(), remove);
            assert_eq!(t2.atoms[1].element(), Element::N);
        }
    }
    #[test]
    fn new_bond_then_atom_removal_preserves_appended_bond_mapping_and_endpoint_renumbering() {
        let r = reaction("([C:1].[O:3].[C:2])>>[C:1][C:2]");
        let t = topology(&[Element::C, Element::O, Element::C], &[]);
        let out = run(&r, &t, false).unwrap();
        assert!(out.changed && out.clears_computed_properties);
        let (t, m) = out.change.unwrap();
        assert_eq!(t.atoms.len(), 2);
        assert_eq!(
            (t.bonds[0].begin(), t.bonds[0].end()),
            (AtomId::new(0), AtomId::new(1))
        );
        assert_eq!(m.bonds.new_to_old, [None]);
        assert!(m.bonds.old_to_new.is_empty());
    }
    #[test]
    fn bond_type_update_keeps_cache_and_computed_state_without_a_removal_effect() {
        let r = reaction("[C:1][C:2]>>[C:1]=[C:2]");
        let mut t = topology(&[Element::C; 2], &[(0, 1, BondOrder::Single)]);
        t.bonds[0]
            .set_computed_prop("kept", PropertyValue::Int(1))
            .unwrap();
        t.atoms[0].set_source_valence_facts(SourceAtomValenceFacts {
            explicit_valence: 4,
            implicit_valence: 0,
        });
        let out = run(&r, &t, false).unwrap();
        assert!(out.changed && !out.clears_computed_properties);
        let (t, m) = out.change.unwrap();
        assert_eq!(t.bonds[0].order(), BondOrder::Double);
        assert!(t.bonds[0].prop("kept").is_some());
        assert_eq!(t.atoms[0].source_valence_facts().explicit_valence, 4);
        assert_eq!(m, TopologyMapping::identity(2, 1));
    }
    #[test]
    fn native_one_match_cap_precedes_protection_filter_without_refilling() {
        let mut r = reaction("[C:1]>>[N:1]");
        r.match_params.max_matches = 0;
        r.match_params.uniquify = true;
        let mut t = topology(&[Element::C; 2], &[]);
        t.atoms[0]
            .set_prop("_protected", PropertyValue::Bool(false))
            .unwrap();
        let out = run(&r, &t, false).unwrap();
        assert!(!out.changed && out.change.is_none());
        assert_eq!(r.match_params.max_matches, 0);
        assert!(r.match_params.uniquify);
    }
    #[test]
    fn helper_failure_after_private_atom_write_keeps_original_input_unchanged() {
        let mut r = reaction("[C:1]>>[N:1]");
        r.products[0].atoms_mut()[0]
            .set_prop("_QueryMass", PropertyValue::String("bad".into()))
            .unwrap();
        let t = topology(&[Element::C], &[]);
        let before = t.clone();
        assert!(run(&r, &t, false).is_err());
        assert_eq!(t, before);
    }
    #[test]
    fn source_commit_endpoint_parse_error_is_structural_and_preserves_caller_blocks() {
        let r = reaction("[C:1][O:2]>>[C:1]");
        let mut t = topology(
            &[Element::C, Element::O, Element::N],
            &[(0, 1, BondOrder::Single), (0, 2, BondOrder::Single)],
        );
        t.bonds[1]
            .set_prop(
                "_MolFileBondEndPts",
                PropertyValue::String("(2 1 nope)".into()),
            )
            .unwrap();
        let before = t.clone();
        assert!(matches!(
            run(&r, &t, false),
            Err(ReactionApplyError::Edit(
                cosmolkit_model::TopologyEditError::BondEndPointsParse { .. }
            ))
        ));
        assert_eq!(t, before);
    }
}
