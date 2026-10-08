use crate::{Reaction, ReactionParseError, ReactionRole};
use cosmolkit_model::{AtomId, QueryGraph};
use cosmolkit_types::ChiralTag;
use std::collections::BTreeMap;

fn neighbor_order(
    graph: &QueryGraph,
    atom: AtomId,
    other: (&QueryGraph, AtomId, ReactionRole, usize),
    role: ReactionRole,
    template: usize,
) -> Result<(u32, Vec<i32>), ReactionParseError> {
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
    let degree = template_atom_degree_source(graph, atom, role, template)?;
    let neighbors = &graph.adjacency()[atom.index()];
    let mut order = Vec::with_capacity(degree as usize);
    let mut unmapped = 0u32;
    for &(neighbor, _) in neighbors {
        let map = crate::validation::atom_map(&graph.atoms()[neighbor], role, template)
            .map_err(|source| ReactionParseError::TemplateProperty { source })?
            .unwrap_or(0); // Atom::getAtomMapNum's source-defined absent value.
        if map > 0 {
            order.push(map);
        } else {
            order.push(-1);
            unmapped = unmapped.wrapping_add(1);
        }
    }
    // Native reads atom2's degree only after all neighbor properties.
    if degree < template_atom_degree_source(other.0, other.1, other.2, other.3)? {
        order.push(-1);
        unmapped = unmapped.wrapping_add(1);
    }
    Ok((unmapped, order))
}

fn template_atom_degree_source(
    graph: &QueryGraph,
    atom: AtomId,
    role: ReactionRole,
    template: usize,
) -> Result<u32, ReactionParseError> {
    // RDKit❗✔️: unsigned int Atom::getDegree() const {
    // RDKit❗✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
    // RDKit❗✔️: }
    // RDKit❗✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
    // RDKit❗✔️:   PRECONDITION(at, "no atom");
    // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
    // RDKit❗✔️:                "atom not associated with this molecule");
    // RDKit❗✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
    // RDKit❗✔️: };
    graph
        .try_atom_degree(atom)
        .ok_or(ReactionParseError::TemplateAtomBounds {
            role,
            template,
            atom,
            atom_count: graph.num_atoms(),
        })
}

fn check_order_overlap(order: &mut [i32], unmapped: u32, reference: &[i32]) -> bool {
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
) -> Result<i32, ReactionParseError> {
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
    let r_degree =
        template_atom_degree_source(reactant, react_atom, ReactionRole::Reactant, react_template)?;
    if r_degree < 3 {
        return Ok(-1);
    }
    let p_degree =
        template_atom_degree_source(product, prod_atom, ReactionRole::Product, prod_template)?;
    if p_degree < 3 || source_degree_difference(r_degree, p_degree)? > 1 {
        return Ok(-1);
    }
    let (r_unmapped, mut r_order) = neighbor_order(
        reactant,
        react_atom,
        (product, prod_atom, ReactionRole::Product, prod_template),
        ReactionRole::Reactant,
        react_template,
    )?;
    if r_unmapped > 1 {
        return Ok(-1);
    }
    let (p_unmapped, mut p_order) = neighbor_order(
        product,
        prod_atom,
        (reactant, react_atom, ReactionRole::Reactant, react_template),
        ReactionRole::Product,
        prod_template,
    )?;
    if p_unmapped > 1 {
        return Ok(-1);
    }
    if check_order_overlap(&mut r_order, r_unmapped, &p_order)
        && check_order_overlap(&mut p_order, p_unmapped, &r_order)
    {
        return cosmolkit_core::count_swaps_to_interconvert(&r_order, &p_order)
            .map(|count| count as u32 as i32)
            .map_err(|source| ReactionParseError::StereoOrder {
                product: prod_template,
                atom: prod_atom,
                source,
            });
    }
    Ok(-1)
}

fn source_degree_difference(reactant: u32, product: u32) -> Result<i32, ReactionParseError> {
    // RDKit❗✔️:       std::abs(static_cast<int>(prodAtom->getDegree()) -
    // RDKit❗✔️:                static_cast<int>(reactAtom->getDegree())) <= 1) {
    // The pinned ABI uses signed32 int conversion. Signed subtraction and
    // abs(INT_MIN) are undefined in C++; retain a structural error for those
    // states rather than silently replacing native arithmetic with abs_diff.
    (product as i32)
        .checked_sub(reactant as i32)
        .and_then(i32::checked_abs)
        .ok_or(ReactionParseError::StereoDegreeArithmetic { reactant, product })
}

fn map_reactant_atoms_source(
    reaction: &Reaction,
    mapping: &mut BTreeMap<i32, (usize, usize)>,
) -> Result<(), ReactionParseError> {
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
    // A template and row locate the actual borrowed atom object; do not use
    // a copied atom ID as an indirection that could select another object.
    for (template, graph) in reaction.reactants.iter().enumerate() {
        for (row, atom) in graph.atoms().iter().enumerate() {
            if let Some(map) = crate::validation::atom_map(atom, ReactionRole::Reactant, template)
                .map_err(|source| ReactionParseError::TemplateProperty { source })?
            {
                mapping.insert(map, (template, row));
            }
        }
    }
    Ok(())
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
    let mut mapping = BTreeMap::new();
    map_reactant_atoms_source(reaction, &mut mapping)?;
    for (product, graph) in reaction.products.iter_mut().enumerate() {
        for row in 0..graph.num_atoms() {
            let atom = &graph.atoms()[row];
            if crate::management::atom_has_property_source(atom, b"molInversionFlag") {
                continue;
            }
            let Some(map) = crate::validation::atom_map(atom, ReactionRole::Product, product)
                .map_err(|source| ReactionParseError::TemplateProperty { source })?
            else {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(4));
                continue;
            };
            let Some(&(react_template, react_row)) = mapping.get(&map) else {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(4));
                continue;
            };
            let reactant = &reaction.reactants[react_template];
            let react_atom = &reactant.atoms()[react_row];
            let r_tag = react_atom.chiral_tag();
            let p_tag = atom.chiral_tag();
            let defined = |tag| !matches!(tag, ChiralTag::Unspecified | ChiralTag::Other);
            if defined(p_tag) {
                if defined(r_tag) {
                    // Source mutation precedes the possibly failing neighbor/property
                    // reads. Preserve that write even when swaps returns an error.
                    graph.atoms_mut()[row].set_mol_inversion_flag(Some(if r_tag == p_tag {
                        2
                    } else {
                        1
                    }));
                    let n_swaps = swaps(
                        reactant,
                        react_atom.id(),
                        graph,
                        AtomId::new(row),
                        react_template,
                        product,
                    )?;
                    if n_swaps >= 0 && n_swaps % 2 != 0 {
                        let value = graph.atoms()[row].mol_inversion_flag();
                        let flag = match value {
                            Some(1) => 2,
                            Some(2) => 1,
                            _ => {
                                return Err(ReactionParseError::StereoInversionFlag {
                                    product,
                                    atom: AtomId::new(row),
                                    flag: value,
                                });
                            }
                        };
                        graph.atoms_mut()[row].set_mol_inversion_flag(Some(flag));
                    }
                } else {
                    graph.atoms_mut()[row].set_mol_inversion_flag(Some(4));
                }
            } else if defined(r_tag) {
                graph.atoms_mut()[row].set_mol_inversion_flag(Some(3));
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod complete_reactant_mapping_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;

    fn graph(maps: &[Option<i32>]) -> QueryGraph {
        let atoms = maps
            .iter()
            .enumerate()
            .map(|(row, map)| {
                let mut atom = QueryAtom::new(AtomId::new(row), AtomSpec::new(Element::C));
                if let Some(map) = map {
                    atom.set_prop("molAtomMapNumber", PropertyValue::Int(*map))
                        .unwrap();
                }
                atom
            })
            .collect();
        QueryGraph::from_parts(atoms, vec![], [], vec![], vec![], vec![]).unwrap()
    }

    #[test]
    fn retains_prefix_and_last_encountered_atom_for_all_present_maps() {
        let mut reaction = Reaction::new();
        reaction.reactants = vec![
            graph(&[Some(7), None, Some(0), Some(-2), Some(7)]),
            graph(&[Some(7), Some(3)]),
        ];
        reaction.products = vec![graph(&[Some(99)])];
        reaction.agents = vec![graph(&[Some(100)])];
        let mut mapping = BTreeMap::from([(80, (9, 9)), (7, (8, 8))]);
        map_reactant_atoms_source(&reaction, &mut mapping).unwrap();
        assert_eq!(
            mapping,
            BTreeMap::from([
                (-2, (0, 3)),
                (0, (0, 2)),
                (3, (1, 1)),
                (7, (1, 0)),
                (80, (9, 9))
            ])
        );
        assert!(reaction.needs_init);
    }

    #[test]
    fn property_failure_preserves_preceding_updates_without_clearing_target() {
        let mut reaction = Reaction::new();
        reaction.reactants = vec![graph(&[Some(5), Some(6), Some(8)])];
        reaction.reactants[0].atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        let mut mapping = BTreeMap::from([(80, (9, 9))]);
        assert!(matches!(
            map_reactant_atoms_source(&reaction, &mut mapping),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
        assert_eq!(mapping, BTreeMap::from([(5, (0, 0)), (80, (9, 9))]));
    }

    #[test]
    fn empty_reactants_leave_target_untouched_and_typed_overflow_propagates() {
        let mut reaction = Reaction::new();
        let mut mapping = BTreeMap::from([(80, (9, 9))]);
        map_reactant_atoms_source(&reaction, &mut mapping).unwrap();
        assert_eq!(mapping.len(), 1);
        reaction.reactants = vec![graph(&[None])];
        reaction.reactants[0].atoms_mut()[0].set_atom_map(Some(u32::MAX));
        assert!(matches!(
            map_reactant_atoms_source(&reaction, &mut mapping),
            Err(ReactionParseError::TemplateProperty {
                source: crate::ReactionValidationError::MapOverflow { .. }
            })
        ));
        assert_eq!(mapping, BTreeMap::from([(80, (9, 9))]));
    }
}

#[cfg(test)]
mod complete_neighbor_order_source_tests {
    use super::*;
    use cosmolkit_model::{AtomSpec, BondId, BondSpec, PropertyValue, QueryAtom, QueryBond};
    use cosmolkit_types::{BondOrder, Element};

    pub(super) fn star(maps: &[Option<i32>], order: &[usize]) -> QueryGraph {
        let atoms = std::iter::once(None)
            .chain(maps.iter().copied())
            .enumerate()
            .map(|(i, map)| {
                let mut atom = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
                if let Some(map) = map {
                    atom.set_prop("molAtomMapNumber", PropertyValue::Int(map))
                        .unwrap();
                }
                atom
            })
            .collect();
        let bonds = order
            .iter()
            .enumerate()
            .map(|(i, &nbr)| {
                QueryBond::new(
                    BondId::new(i),
                    BondSpec::new(AtomId::new(0), AtomId::new(nbr), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(atoms, bonds, [], vec![], vec![], vec![]).unwrap()
    }

    #[test]
    fn follows_incident_edge_order_and_retains_duplicate_positive_maps() {
        let graph = star(&[Some(7), Some(4), Some(7)], &[3, 1, 2]);
        let result = neighbor_order(
            &graph,
            AtomId::new(0),
            (&graph, AtomId::new(0), ReactionRole::Product, 8),
            ReactionRole::Reactant,
            3,
        )
        .unwrap();
        assert_eq!(result, (0, vec![7, 7, 4]));
    }

    #[test]
    fn absent_zero_negative_and_degree_difference_add_source_minus_one_slots() {
        let graph = star(&[None, Some(0), Some(-2), Some(9)], &[1, 2, 3, 4]);
        let larger = star(&[Some(1); 7], &[1, 2, 3, 4, 5, 6, 7]);
        assert_eq!(
            neighbor_order(
                &graph,
                AtomId::new(0),
                (&larger, AtomId::new(0), ReactionRole::Product, 8),
                ReactionRole::Reactant,
                3
            )
            .unwrap(),
            (4, vec![-1, -1, -1, 9, -1])
        );
        assert_eq!(
            neighbor_order(
                &larger,
                AtomId::new(0),
                (&graph, AtomId::new(0), ReactionRole::Product, 8),
                ReactionRole::Reactant,
                3
            )
            .unwrap(),
            (0, vec![1; 7])
        );
    }

    #[test]
    fn neighbor_property_error_precedes_other_atom_degree_read() {
        let mut graph = star(&[Some(1)], &[1]);
        graph.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        assert!(matches!(
            neighbor_order(
                &graph,
                AtomId::new(0),
                (&graph, AtomId::new(99), ReactionRole::Product, 8),
                ReactionRole::Reactant,
                3
            ),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
        graph.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::Int(1))
            .unwrap();
        assert!(matches!(
            neighbor_order(
                &graph,
                AtomId::new(0),
                (&graph, AtomId::new(99), ReactionRole::Product, 8),
                ReactionRole::Reactant,
                3
            ),
            Err(ReactionParseError::TemplateAtomBounds {
                role: ReactionRole::Product,
                template: 8,
                ..
            })
        ));
        assert!(matches!(
            neighbor_order(
                &graph,
                AtomId::new(99),
                (&graph, AtomId::new(0), ReactionRole::Product, 8),
                ReactionRole::Reactant,
                3
            ),
            Err(ReactionParseError::TemplateAtomBounds {
                role: ReactionRole::Reactant,
                template: 3,
                ..
            })
        ));
    }
}

#[cfg(test)]
mod complete_order_overlap_source_tests {
    use super::*;

    #[test]
    fn source_nonzero_unmapped_flag_is_not_a_consumable_slot_budget() {
        let mut order = [-1, 7, -1, -1];
        assert!(check_order_overlap(&mut order, 1, &[0, 9, 9, -3, 4, 7]));
        assert_eq!(order, [0, 7, 9, 4]);
        let mut order = [-1, -1];
        assert!(check_order_overlap(&mut order, u32::MAX, &[2, 3]));
        assert_eq!(order, [2, 3]);
    }

    #[test]
    fn failed_overlap_keeps_prior_replacements_and_stops_at_first_missing() {
        let mut order = [-1, 7];
        assert!(!check_order_overlap(&mut order, 1, &[4, 5, 6]));
        assert_eq!(order, [4, 7]);
        let mut order = [-1, 7];
        assert!(!check_order_overlap(&mut order, 0, &[4, 7]));
        assert_eq!(order, [-1, 7]);
    }

    #[test]
    fn negative_reference_values_are_ignored_and_existing_duplicates_are_membership() {
        let mut order = [2, 2, -1];
        assert!(check_order_overlap(&mut order, 0, &[-1, -9, 2, 2]));
        assert_eq!(order, [2, 2, -1]);
        assert!(check_order_overlap(&mut [], 0, &[-2, -1]));
        assert!(!check_order_overlap(&mut [], 1, &[0]));
        assert!(check_order_overlap(&mut order, 0, &[]));
    }
}

#[cfg(test)]
mod complete_template_swaps_source_tests {
    use super::complete_neighbor_order_source_tests::star;
    use super::*;
    use cosmolkit_model::PropertyValue;

    pub(super) fn count(
        reactant: &QueryGraph,
        product: &QueryGraph,
    ) -> Result<i32, ReactionParseError> {
        swaps(reactant, AtomId::new(0), product, AtomId::new(0), 2, 8)
    }

    #[test]
    fn delegates_first_match_swap_count_without_sorting_or_parity_shortcut() {
        let r = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
        for (order, expected) in [
            ([1, 2, 3], 0),
            ([2, 1, 3], 1),
            ([2, 3, 1], 2),
            ([3, 2, 1], 1),
        ] {
            let p = star(&[Some(1), Some(2), Some(3)], &order);
            assert_eq!(count(&r, &p).unwrap(), expected);
        }
    }

    #[test]
    fn degree_change_uses_one_placeholder_and_two_way_overlap() {
        let three = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
        let four = star(&[Some(1), Some(2), Some(3), Some(4)], &[1, 2, 3, 4]);
        assert_eq!(count(&three, &four).unwrap(), 0);
        assert_eq!(count(&four, &three).unwrap(), 0);
        let missing = star(&[Some(4), Some(5), Some(6)], &[1, 2, 3]);
        assert_eq!(count(&three, &missing).unwrap(), -1);
        let five = star(&[Some(1); 5], &[1, 2, 3, 4, 5]);
        assert_eq!(count(&three, &five).unwrap(), -1);
    }

    #[test]
    fn source_short_circuits_degree_then_reactant_unmapped_before_product_properties() {
        let small = star(&[Some(1), Some(2)], &[1, 2]);
        assert_eq!(
            swaps(&small, AtomId::new(0), &small, AtomId::new(99), 2, 8).unwrap(),
            -1
        );
        let r = star(&[None, None, Some(3)], &[1, 2, 3]);
        let mut p = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
        p.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        assert_eq!(count(&r, &p).unwrap(), -1);
        let mapped = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
        assert!(matches!(
            count(&mapped, &p),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
    }

    #[test]
    fn duplicate_membership_does_not_hide_source_swap_invariant_failure() {
        let r = star(&[Some(1), Some(1), Some(2)], &[1, 2, 3]);
        let p = star(&[Some(1), Some(2), Some(2)], &[1, 2, 3]);
        assert!(matches!(
            count(&r, &p),
            Err(ReactionParseError::StereoOrder { product: 8, .. })
        ));
    }

    #[test]
    fn signed_source_degree_arithmetic_keeps_pinned_cast_and_undefined_state_errors() {
        assert_eq!(source_degree_difference(3, 4).unwrap(), 1);
        assert_eq!(source_degree_difference(4, 3).unwrap(), 1);
        assert_eq!(source_degree_difference(3, u32::MAX).unwrap(), 4);
        assert_eq!(source_degree_difference(1 << 31, 1 << 31).unwrap(), 0);
        assert!(matches!(
            source_degree_difference(3, 1 << 31),
            Err(ReactionParseError::StereoDegreeArithmetic { .. })
        ));
        assert!(matches!(
            source_degree_difference(0, 1 << 31),
            Err(ReactionParseError::StereoDegreeArithmetic { .. })
        ));
    }
}

#[cfg(test)]
mod complete_products_stereo_source_tests {
    use super::complete_neighbor_order_source_tests::star;
    use super::*;
    use cosmolkit_model::PropertyValue;

    fn tagged(map: Option<i32>, tag: ChiralTag) -> QueryGraph {
        let mut graph = star(&[], &[]);
        if let Some(map) = map {
            graph.atoms_mut()[0]
                .set_prop("molAtomMapNumber", PropertyValue::Int(map))
                .unwrap();
        }
        graph.atoms_mut()[0].set_chiral_tag(tag);
        graph
    }

    #[test]
    fn every_modeled_chiral_tag_follows_literal_source_unspecified_other_exclusions() {
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
        for r in tags {
            for p in tags {
                let mut reaction = Reaction::new();
                reaction.reactants = vec![tagged(Some(10), r)];
                reaction.products = vec![tagged(Some(10), p)];
                reaction.needs_init = false;
                update_products_stereochem(&mut reaction).unwrap();
                let defined = |t| !matches!(t, ChiralTag::Unspecified | ChiralTag::Other);
                let expected = if defined(p) {
                    Some(if defined(r) {
                        if r == p { 2 } else { 1 }
                    } else {
                        4
                    })
                } else if defined(r) {
                    Some(3)
                } else {
                    None
                };
                assert_eq!(
                    reaction.products[0].atoms()[0].mol_inversion_flag(),
                    expected,
                    "{r:?}->{p:?}"
                );
                assert!(!reaction.needs_init);
                assert_eq!(reaction.products[0].atoms()[0].chiral_tag(), p);
            }
        }
    }

    #[test]
    fn absent_or_new_mapping_creates_flag_even_without_product_stereo() {
        let mut reaction = Reaction::new();
        reaction.products = vec![
            tagged(None, ChiralTag::Unspecified),
            tagged(Some(99), ChiralTag::Other),
        ];
        update_products_stereochem(&mut reaction).unwrap();
        assert!(
            reaction
                .products
                .iter()
                .all(|g| g.atoms()[0].mol_inversion_flag() == Some(4))
        );
    }

    #[test]
    fn existing_typed_or_wrong_tag_raw_flag_skips_map_reads_and_normalization() {
        let mut reaction = Reaction::new();
        let mut raw = tagged(None, ChiralTag::TetrahedralCw);
        raw.atoms_mut()[0]
            .set_prop("molInversionFlag", PropertyValue::Bool(false))
            .unwrap();
        raw.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        let mut typed = raw.clone();
        typed.atoms_mut()[0].set_mol_inversion_flag(Some(-9));
        reaction.products = vec![raw, typed];
        update_products_stereochem(&mut reaction).unwrap();
        assert_eq!(reaction.products[0].atoms()[0].mol_inversion_flag(), None);
        assert_eq!(
            reaction.products[1].atoms()[0].mol_inversion_flag(),
            Some(-9)
        );
        assert!(reaction.products.iter().all(|g| g.atoms()[0].prop("molInversionFlag") == Some(&PropertyValue::Bool(false))));
    }

    #[test]
    fn swap_parity_toggles_equal_and_unequal_tags_only_for_nonnegative_odd_counts() {
        for (p_tag, order, expected) in [
            (ChiralTag::TetrahedralCw, [1, 2, 3], 2),
            (ChiralTag::TetrahedralCw, [2, 1, 3], 1),
            (ChiralTag::TetrahedralCw, [2, 3, 1], 2),
            (ChiralTag::TetrahedralCcw, [2, 1, 3], 2),
        ] {
            let mut r = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
            let mut p = star(&[Some(1), Some(2), Some(3)], &order);
            r.atoms_mut()[0].set_atom_map(Some(10));
            p.atoms_mut()[0].set_atom_map(Some(10));
            r.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
            p.atoms_mut()[0].set_chiral_tag(p_tag);
            let mut reaction = Reaction::new();
            reaction.reactants = vec![r];
            reaction.products = vec![p];
            update_products_stereochem(&mut reaction).unwrap();
            assert_eq!(
                reaction.products[0].atoms()[0].mol_inversion_flag(),
                Some(expected)
            );
        }
    }

    #[test]
    fn initial_inversion_write_survives_later_swap_property_failure() {
        let mut r = star(&[Some(1), Some(2), Some(3)], &[1, 2, 3]);
        let mut p = r.clone();
        for graph in [&mut r, &mut p] {
            graph.atoms_mut()[0].set_atom_map(Some(10));
            graph.atoms_mut()[0].set_chiral_tag(ChiralTag::TetrahedralCw);
        }
        p.atoms_mut()[1]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        let mut reaction = Reaction::new();
        reaction.reactants = vec![r];
        reaction.products = vec![p];
        assert!(matches!(
            update_products_stereochem(&mut reaction),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
        assert_eq!(
            reaction.products[0].atoms()[0].mol_inversion_flag(),
            Some(2)
        );
        assert!(
            reaction.products[0].atoms()[1..]
                .iter()
                .all(|a| a.mol_inversion_flag().is_none())
        );
    }

    #[test]
    fn reactant_map_failure_precedes_any_product_mutation() {
        let mut reaction = Reaction::new();
        let mut r = tagged(None, ChiralTag::Unspecified);
        r.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        reaction.reactants = vec![r];
        reaction.products = vec![tagged(None, ChiralTag::Unspecified)];
        assert!(matches!(
            update_products_stereochem(&mut reaction),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
        assert_eq!(reaction.products[0].atoms()[0].mol_inversion_flag(), None);
    }

    #[test]
    fn last_duplicate_reactant_mapping_including_zero_and_negative_wins() {
        for map in [0, -2, 10] {
            let mut reaction = Reaction::new();
            reaction.reactants = vec![
                tagged(Some(map), ChiralTag::TetrahedralCw),
                tagged(Some(map), ChiralTag::TetrahedralCcw),
            ];
            reaction.products = vec![tagged(Some(map), ChiralTag::TetrahedralCw)];
            update_products_stereochem(&mut reaction).unwrap();
            assert_eq!(
                reaction.products[0].atoms()[0].mol_inversion_flag(),
                Some(1)
            );
        }
    }

    #[test]
    fn product_encounter_order_retains_prefix_writes_on_first_map_error() {
        let mut reaction = Reaction::new();
        let good = tagged(None, ChiralTag::Unspecified);
        let mut bad = good.clone();
        bad.atoms_mut()[0]
            .set_prop("molAtomMapNumber", PropertyValue::Bool(true))
            .unwrap();
        reaction.products = vec![good.clone(), bad, good];
        assert!(matches!(
            update_products_stereochem(&mut reaction),
            Err(ReactionParseError::TemplateProperty { .. })
        ));
        assert_eq!(
            reaction.products[0].atoms()[0].mol_inversion_flag(),
            Some(4)
        );
        assert_eq!(reaction.products[1].atoms()[0].mol_inversion_flag(), None);
        assert_eq!(reaction.products[2].atoms()[0].mol_inversion_flag(), None);
    }
}
