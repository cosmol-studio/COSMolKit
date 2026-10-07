use crate::{Reaction, ReactionParseError, ReactionRole};
use cosmolkit_model::{AtomId, QueryGraph};
use cosmolkit_types::ChiralTag;
use std::collections::BTreeMap;

fn neighbor_order(
    graph: &QueryGraph,
    atom: AtomId,
    other_degree: usize,
    role: ReactionRole,
    template: usize,
) -> Result<(usize, Vec<i32>), ReactionParseError> {
    // RDKit❗✔️: std::pair<unsigned int, std::vector<int>> getNbrOrder(const Atom *atom1,
    // RDKit❗✔️:                                                       const Atom *atom2) {
    // RDKit❗✔️:   std::vector<int> order;
    // RDKit❗✔️:   order.reserve(atom1->getDegree());
    // RDKit❗✔️:   unsigned nUnmapped = 0;
    // RDKit❗✔️:   for (const auto nbrAtom : atom1->getOwningMol().atomNeighbors(atom1)) {
    // RDKit❗✔️:     if (nbrAtom->getAtomMapNum() > 0) {
    // RDKit❗✔️:       order.push_back(nbrAtom->getAtomMapNum());
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       order.push_back(-1);
    // RDKit❗✔️:       ++nUnmapped;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (atom1->getDegree() < atom2->getDegree()) {
    // RDKit❗✔️:     order.push_back(-1);
    // RDKit❗✔️:     ++nUnmapped;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return {nUnmapped, order};
    // RDKit❗✔️: }
    let neighbors = &graph.adjacency()[atom.index()];
    let mut order = Vec::with_capacity(neighbors.len());
    let mut unmapped = 0;
    for &(neighbor, _) in neighbors {
        match crate::validation::atom_map(&graph.atoms()[neighbor], role, template)
            .map_err(|source| ReactionParseError::TemplateProperty { source })?
            .filter(|&map| map > 0)
        {
            Some(map) => order.push(map),
            None => {
                order.push(-1);
                unmapped += 1;
            }
        }
    }
    if neighbors.len() < other_degree {
        order.push(-1);
        unmapped += 1;
    }
    Ok((unmapped, order))
}

fn check_order_overlap(order: &mut [i32], unmapped: usize, reference: &[i32]) -> bool {
    // RDKit❗✔️: bool checkOrderOverlap(std::vector<int> &order, unsigned int nUnmapped,
    // RDKit❗✔️:                        const std::vector<int> &refOrder) {
    // RDKit❗✔️:   bool allFound = true;
    // RDKit❗✔️:   for (auto elem : refOrder) {
    // RDKit❗✔️:     if (elem >= 0) {
    // RDKit❗✔️:       if (std::find(order.begin(), order.end(), elem) == order.end()) {
    // RDKit❗✔️:         // this one was not there, is there an unmapped slot for
    // RDKit❗✔️:         // it (i.e. a -1 value in the order)?
    // RDKit❗✔️:         if (nUnmapped) {
    // RDKit❗✔️:           auto negOne = std::find(order.begin(), order.end(), -1);
    // RDKit❗✔️:           if (negOne != order.end()) {
    // RDKit❗✔️:             *negOne = elem;
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             allFound = false;
    // RDKit❗✔️:             break;
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else {
    // RDKit❗✔️:           allFound = false;
    // RDKit❗✔️:           break;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return allFound;
    // RDKit❗✔️: }
    for &elem in reference {
        if elem >= 0 && !order.contains(&elem) {
            if unmapped != 0 {
                if let Some(slot) = order.iter_mut().find(|slot| **slot == -1) {
                    *slot = elem;
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
    }
    true
}

fn swaps(
    reactant: &QueryGraph,
    react_atom: AtomId,
    product: &QueryGraph,
    prod_atom: AtomId,
    react_template: usize,
    prod_template: usize,
) -> Result<Option<usize>, ReactionParseError> {
    // RDKit❗✔️: int countSwapsBetweenReactantAndProduct(const Atom *reactAtom,
    // RDKit❗✔️:                                         const Atom *prodAtom) {
    // RDKit❗✔️:   PRECONDITION(reactAtom, "bad atom");
    // RDKit❗✔️:   PRECONDITION(prodAtom, "bad atom");
    // RDKit❗✔️:   if (reactAtom->getDegree() >= 3 && prodAtom->getDegree() >= 3 &&
    // RDKit❗✔️:       std::abs(static_cast<int>(prodAtom->getDegree()) -
    // RDKit❗✔️:                static_cast<int>(reactAtom->getDegree())) <= 1) {
    // RDKit❗✔️:     std::vector<int> reactOrder;
    // RDKit❗✔️:     unsigned int nReactUnmapped;
    // RDKit❗✔️:     std::tie(nReactUnmapped, reactOrder) = getNbrOrder(reactAtom, prodAtom);
    // RDKit❗✔️:     if (nReactUnmapped <= 1) {
    // RDKit❗✔️:       std::vector<int> prodOrder;
    // RDKit❗✔️:       unsigned int nProdUnmapped;
    // RDKit❗✔️:       std::tie(nProdUnmapped, prodOrder) = getNbrOrder(prodAtom, reactAtom);
    // RDKit❗✔️:       if (nProdUnmapped <= 1) {
    // RDKit❗✔️:         // check that each element of the product mappings is
    // RDKit❗✔️:         // in the reactant mappings
    // RDKit❗✔️:         if (checkOrderOverlap(reactOrder, nReactUnmapped, prodOrder)) {
    // RDKit❗✔️:           // found a match for all the product atoms, what about all
    // RDKit❗✔️:           // the reactant atoms?
    // RDKit❗✔️:           if (checkOrderOverlap(prodOrder, nProdUnmapped, reactOrder)) {
    // RDKit❗✔️:             return countSwapsToInterconvert(reactOrder, prodOrder);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return -1;
    // RDKit❗✔️: }
    let r_degree = reactant.adjacency()[react_atom.index()].len();
    let p_degree = product.adjacency()[prod_atom.index()].len();
    if r_degree < 3 || p_degree < 3 || r_degree.abs_diff(p_degree) > 1 {
        return Ok(None);
    }
    let (r_unmapped, mut r_order) = neighbor_order(
        reactant,
        react_atom,
        p_degree,
        ReactionRole::Reactant,
        react_template,
    )?;
    if r_unmapped > 1 {
        return Ok(None);
    }
    let (p_unmapped, mut p_order) = neighbor_order(
        product,
        prod_atom,
        r_degree,
        ReactionRole::Product,
        prod_template,
    )?;
    if p_unmapped > 1 {
        return Ok(None);
    }
    if check_order_overlap(&mut r_order, r_unmapped, &p_order)
        && check_order_overlap(&mut p_order, p_unmapped, &r_order)
    {
        return cosmolkit_core::count_swaps_to_interconvert(&r_order, &p_order)
            .map(Some)
            .map_err(|source| ReactionParseError::StereoOrder {
                product: prod_template,
                atom: prod_atom,
                source,
            });
    }
    Ok(None)
}

pub(crate) fn update_products_stereochem(
    reaction: &mut Reaction,
) -> Result<(), ReactionParseError> {
    // RDKit❗✔️: void updateProductsStereochem(ChemicalReaction *rxn) {
    // RDKit❗✔️:   std::map<int, Atom *> reactantMapping;
    // RDKit❗✔️:   getMappingNumAtomIdxMapReactants(*rxn, reactantMapping);
    // RDKit❗✔️:   for (MOL_SPTR_VECT::const_iterator prodIt = rxn->beginProductTemplates();
    // RDKit❗✔️:        prodIt != rxn->endProductTemplates(); ++prodIt) {
    // RDKit❗✔️:     for (auto prodAtom : (*prodIt)->atoms()) {
    // RDKit❗✔️:       if (prodAtom->hasProp(common_properties::molInversionFlag)) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (!prodAtom->hasProp(common_properties::molAtomMapNumber)) {
    // RDKit❗✔️:         // if we have stereochemistry specified, it's automatically
    // RDKit❗✔️:         // creating stereochem:
    // RDKit❗✔️:         prodAtom->setProp(common_properties::molInversionFlag, 4);
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       int mapNum;
    // RDKit❗✔️:       prodAtom->getProp(common_properties::molAtomMapNumber, mapNum);
    // RDKit❗✔️:       if (reactantMapping.find(mapNum) != reactantMapping.end()) {
    // RDKit❗✔️:         const auto reactAtom = reactantMapping[mapNum];
    // RDKit❗✔️:         if (prodAtom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:             prodAtom->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit❗✔️:           if (reactAtom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:               reactAtom->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit❗✔️:             // both have stereochem specified, we're either preserving
    // RDKit❗✔️:             // or inverting
    // RDKit❗✔️:             if (reactAtom->getChiralTag() == prodAtom->getChiralTag()) {
    // RDKit❗✔️:               prodAtom->setProp(common_properties::molInversionFlag, 2);
    // RDKit❗✔️:             } else {
    // RDKit❗✔️:               // FIX: this is technically fragile: it should be checking
    // RDKit❗✔️:               // if the atoms both have tetrahedral chirality. However,
    // RDKit❗✔️:               // at the moment that's the only chirality available, so
    // RDKit❗✔️:               // there's no need to go monkeying around.
    // RDKit❗✔️:               prodAtom->setProp(common_properties::molInversionFlag, 1);
    // RDKit❗✔️:             }
    // RDKit❗✔️:
    // RDKit❗✔️:             // FIX this should move out into a separate function
    // RDKit❗✔️:             // last thing to check here: if the ordering of the bonds
    // RDKit❗✔️:             // around the atom changed from reactants->products then we
    // RDKit❗✔️:             // may need to adjust the inversion flag
    // RDKit❗✔️:             int nSwaps =
    // RDKit❗✔️:                 countSwapsBetweenReactantAndProduct(reactAtom, prodAtom);
    // RDKit❗✔️:             if (nSwaps >= 0 && nSwaps % 2) {
    // RDKit❗✔️:               auto mival =
    // RDKit❗✔️:                   prodAtom->getProp<int>(common_properties::molInversionFlag);
    // RDKit❗✔️:               if (mival == 1) {
    // RDKit❗✔️:                 prodAtom->setProp(common_properties::molInversionFlag, 2);
    // RDKit❗✔️:               } else if (mival == 2) {
    // RDKit❗✔️:                 prodAtom->setProp(common_properties::molInversionFlag, 1);
    // RDKit❗✔️:               } else {
    // RDKit❗✔️:                 CHECK_INVARIANT(false, "inconsistent molInversionFlag");
    // RDKit❗✔️:               }
    // RDKit❗✔️:             }
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             // stereochem in the product, but not in the reactant
    // RDKit❗✔️:             prodAtom->setProp(common_properties::molInversionFlag, 4);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         } else if (reactantMapping[mapNum]->getChiralTag() !=
    // RDKit❗✔️:                        Atom::CHI_UNSPECIFIED &&
    // RDKit❗✔️:                    reactantMapping[mapNum]->getChiralTag() != Atom::CHI_OTHER) {
    // RDKit❗✔️:           // stereochem in the reactant, but not the product:
    // RDKit❗✔️:           prodAtom->setProp(common_properties::molInversionFlag, 3);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       } else {
    // RDKit❗✔️:         // introduction of new stereocenter by the reaction
    // RDKit❗✔️:         prodAtom->setProp(common_properties::molInversionFlag, 4);
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:
    // RDKit❗✔️: void getMappingNumAtomIdxMapReactants(
    // RDKit❗✔️:     const ChemicalReaction &rxn, std::map<int, Atom *> &reactantAtomMapping) {
    // RDKit❗✔️:   for (auto reactIt = rxn.beginReactantTemplates();
    // RDKit❗✔️:        reactIt != rxn.endReactantTemplates(); ++reactIt) {
    // RDKit❗✔️:     for (const auto atom : (*reactIt)->atoms()) {
    // RDKit❗✔️:       int reactMapNum;
    // RDKit❗✔️:       if (atom->getPropIfPresent(common_properties::molAtomMapNumber,
    // RDKit❗✔️:                                  reactMapNum)) {
    // RDKit❗✔️:         reactantAtomMapping[reactMapNum] = atom;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // RDKit❗✔️:
    let mut mapping = BTreeMap::new();
    for (template, graph) in reaction.reactants.iter().enumerate() {
        for atom in graph.atoms() {
            if let Some(map) = crate::validation::atom_map(atom, ReactionRole::Reactant, template)
                .map_err(|source| ReactionParseError::TemplateProperty { source })?
            {
                mapping.insert(map, (template, atom.id()));
            }
        }
    }
    for (product, graph) in reaction.products.iter_mut().enumerate() {
        for row in 0..graph.num_atoms() {
            let atom = &graph.atoms()[row];
            if atom.mol_inversion_flag().is_some() {
                continue;
            }
            let Some(map) = crate::validation::atom_map(atom, ReactionRole::Product, product)
                .map_err(|source| ReactionParseError::TemplateProperty { source })?
            else {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(4));
                continue;
            };
            let Some(&(react_template, react_atom)) = mapping.get(&map) else {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(4));
                continue;
            };
            let reactant = &reaction.reactants[react_template];
            let r_tag = reactant.atoms()[react_atom.index()].chiral_tag();
            let p_tag = atom.chiral_tag();
            let defined = |tag| !matches!(tag, ChiralTag::Unspecified | ChiralTag::Other);
            let flag = if defined(p_tag) {
                if defined(r_tag) {
                    let mut flag = if r_tag == p_tag { 2 } else { 1 };
                    if swaps(
                        reactant,
                        react_atom,
                        graph,
                        AtomId::new(row),
                        react_template,
                        product,
                    )?
                    .is_some_and(|n| n % 2 != 0)
                    {
                        flag = if flag == 1 { 2 } else { 1 };
                    }
                    Some(flag)
                } else {
                    Some(4)
                }
            } else if defined(r_tag) {
                Some(3)
            } else {
                None
            };
            if let Some(flag) = flag {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(flag));
            }
        }
    }
    Ok(())
}
