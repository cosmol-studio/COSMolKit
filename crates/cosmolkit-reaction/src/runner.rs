use crate::materialize::{
    ProductBuilder, add_neighbors, atom_mappings, convert_template, invariant,
    transfer_atom_properties, transfer_bond_properties,
};
use crate::product_stereo::{
    copy_enhanced_groups, copy_template_groups, correct_matched_chirality,
    correct_unmatched_chirality, propagate_coordinates, update_stereo_bonds,
};
use crate::{
    Reaction, ReactionCoordinateSelection, ReactionInput, ReactionProduct, ReactionProductError,
    ReactionRunError, ReactionRunParams, ReactionSingleRunParams, ReactionValidationParams,
};
use cosmolkit_core::ValenceModel;
use cosmolkit_model::{
    AdjacencyList, AtomId, Conformer3D, CoordinateBlock, MoleculeProperties, QueryGraph,
    TopologyBlock,
};
use cosmolkit_types::{BondDirection, BondOrder};
use std::borrow::Cow;

fn initialized(reaction: &Reaction) -> Result<Cow<'_, Reaction>, ReactionRunError> {
    // ROOT-approved D1: run initializes a private temporary only when needed.
    // Immutable caller state and warning/error reporting remain intact.
    if reaction.is_initialized() {
        Ok(Cow::Borrowed(reaction))
    } else {
        Ok(Cow::Owned(crate::initialize_reaction(
            reaction,
            &ReactionValidationParams::default(),
        )?))
    }
}

fn add_reactant(
    reaction: &Reaction,
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    matches: &[usize],
    template: &QueryGraph,
    template_index: usize,
    mut conformer: Option<(&mut Vec<[f64; 3]>, &mut bool)>,
    selection: ReactionCoordinateSelection,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addReactantAtomsAndBonds
    // RDKit❗✔️: void addReactantAtomsAndBonds(const ChemicalReaction &rxn, RWMOL_SPTR product,
    // RDKit❗✔️:                               const ROMOL_SPTR reactantSptr,
    // RDKit❗✔️:                               const MatchVectType &match,
    // RDKit❗✔️:                               const ROMOL_SPTR reactantTemplate,
    // RDKit❗✔️:                               Conformer *productConf,
    // RDKit❗✔️: 			      unsigned int reactantId) {
    // RDKit❗✔️:   // start by looping over all matches and marking the reactant atoms that
    // RDKit❗✔️:   // have already been "added" by virtue of being in the product. We'll also
    // RDKit❗✔️:   // mark "skipped" atoms: those that are in the match, but not in this
    // RDKit❗✔️:   // particular product (or, perhaps, not in any product)
    // RDKit❗✔️:   // At the same time we'll set up a map between the indices of those
    // RDKit❗✔️:   // atoms and their index in the product.
    // RDKit❗✔️:   ReactantProductAtomMapping *mapping = getAtomMappingsReactantProduct(
    // RDKit❗✔️:       match, *reactantTemplate, product, reactantSptr->getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   boost::dynamic_bitset<> visitedAtoms(reactantSptr->getNumAtoms());
    // RDKit❗✔️:
    // RDKit❗✔️:   const ROMol *reactant = reactantSptr.get();
    // RDKit❗✔️:
    // RDKit❗✔️:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗✔️:   // Loop over the bonds in the product and look for those that have
    // RDKit❗✔️:   // the NullBond property set. These are bonds for which no information
    // RDKit❗✔️:   // (other than their existence) was provided in the template
    // RDKit❗✔️:   setReactantBondPropertiesToProduct(product, *reactant, mapping);
    // RDKit❗✔️:
    // RDKit❗✔️:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗✔️:   // Loop over the atoms in the match that were added to the product
    // RDKit❗✔️:   // From the corresponding atom in the reactant, do a graph traversal
    // RDKit❗✔️:   // to find other connected atoms that should be added:
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<const Atom *> chiralAtomsToCheck;
    // RDKit❗✔️:   for (const auto &matchIdx : match) {
    // RDKit❗✔️:     int reactantAtomIdx = matchIdx.second;
    // RDKit❗✔️:     if (mapping->mappedAtoms[reactantAtomIdx]) {
    // RDKit❗✔️:       CHECK_INVARIANT(mapping->reactProdAtomMap.find(reactantAtomIdx) !=
    // RDKit❗✔️:                           mapping->reactProdAtomMap.end(),
    // RDKit❗✔️:                       "mapped reactant atom not present in product.");
    // RDKit❗✔️:
    // RDKit❗✔️:       const Atom *reactantAtom = reactant->getAtomWithIdx(reactantAtomIdx);
    // RDKit❗✔️:       for (unsigned i = 0;
    // RDKit❗✔️:            i < mapping->reactProdAtomMap[reactantAtomIdx].size(); i++) {
    // RDKit❗✔️:         // here's a pointer to the atom in the product:
    // RDKit❗✔️:         unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
    // RDKit❗✔️:         Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗✔️:         setReactantAtomPropertiesToProduct(productAtom, *reactantAtom,
    // RDKit❗✔️:                                            rxn.getImplicitPropertiesFlag(), reactantId);
    // RDKit❗✔️:         if (reactantAtom->hasQuery()) {
    // RDKit❗✔️:           // finally: if the reactant atom is a query we should copy over the
    // RDKit❗✔️:           // query information. We need to replace the atom to do this
    // RDKit❗✔️:           QueryAtom newAtom(*productAtom);
    // RDKit❗✔️:           newAtom.setQuery(reactantAtom->getQuery()->copy());
    // RDKit❗✔️:           // replaceAtom copies
    // RDKit❗✔️:           product->replaceAtom(productAtomIdx, &newAtom);
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:       // now traverse:
    // RDKit❗✔️:       addReactantNeighborsToProduct(*reactant, *reactantAtom, product,
    // RDKit❗✔️:                                     visitedAtoms, chiralAtomsToCheck, mapping, reactantId);
    // RDKit❗✔️:
    // RDKit❗✔️:       // now that we've added all the reactant's neighbors, check to see if
    // RDKit❗✔️:       // it is chiral in the reactant but is not in the reaction. If so
    // RDKit❗✔️:       // we need to worry about its chirality
    // RDKit❗✔️:       checkAndCorrectChiralityOfMatchingAtomsInProduct(
    // RDKit❗✔️:           *reactant, reactantAtomIdx, *reactantAtom, product, mapping);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }  // end of loop over matched atoms
    // RDKit❗✔️:
    // RDKit❗✔️:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗✔️:   // now we need to loop over atoms from the reactants that were chiral but
    // RDKit❗✔️:   // not directly involved in the reaction in order to make sure their
    // RDKit❗✔️:   // chirality hasn't been disturbed
    // RDKit❗✔️:   checkAndCorrectChiralityOfProduct(chiralAtomsToCheck, product, mapping);
    // RDKit❗✔️:
    // RDKit❗✔️:   updateStereoBonds(product, *reactant, mapping);
    // RDKit❗✔️:
    // RDKit❗✔️:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗✔️:   // Copy enhanced StereoGroup data from reactant to product if it is
    // RDKit❗✔️:   // still valid. Uses ChiralTag checks above.
    // RDKit❗✔️:   copyEnhancedStereoGroups(*reactant, product, *mapping);
    // RDKit❗✔️:
    // RDKit❗✔️:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗✔️:   // finally we may need to set the coordinates in the product conformer:
    // RDKit❗✔️:   if (productConf) {
    // RDKit❗✔️:     productConf->resize(product->getNumAtoms());
    // RDKit❗✔️:     generateProductConformers(productConf, *reactant, mapping);
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   delete (mapping);
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let mut mapping = atom_mappings(
        matches,
        template,
        template_index,
        product,
        input.topology.atoms.len(),
    )?;
    let mut visited = vec![false; input.topology.atoms.len()];
    transfer_bond_properties(product, input, template_index, &mapping)?;
    let mut chiral_to_check = Vec::new();
    for &reactant_atom in matches {
        if mapping.mapped[reactant_atom] {
            let rows = mapping
                .reactant_to_product
                .get(&reactant_atom)
                .ok_or_else(|| {
                    invariant(
                        "addReactantAtomsAndBonds",
                        "mapped reactant atom not present in product.",
                        Some(reactant_atom),
                        None,
                        None,
                    )
                })?;
            let r = &input.topology.atoms[reactant_atom];
            for &row in rows {
                transfer_atom_properties(
                    &mut product.topology.atoms[row],
                    r,
                    reaction.implicit_properties(),
                    template_index,
                )?;
                product.atom_origins[row] = Some(crate::ReactionRowOrigin {
                    input: template_index,
                    row: AtomId::new(reactant_atom),
                });
            }
            // Independent G7 query-bearing reagent inputs are excluded; the
            // borrowed input carries canonical concrete Atom/Bond rows only.
            add_neighbors(
                input,
                reactant_atom,
                product,
                &mut visited,
                &mut chiral_to_check,
                &mut mapping,
                template_index,
            )?;
            correct_matched_chirality(product, &input, &mapping, reactant_atom)?;
        }
    }
    correct_unmatched_chirality(product, &input, &mapping, &chiral_to_check)?;
    update_stereo_bonds(product, &input, &mapping)?;
    copy_enhanced_groups(product, &input, &mapping)?;
    if let Some((points, is_3d)) = conformer.as_mut() {
        propagate_coordinates(
            points,
            is_3d,
            product.topology.atoms.len(),
            &input,
            &mapping,
            selection,
        )?;
    }
    Ok(())
}

fn one_product_set(
    reaction: &Reaction,
    inputs: &[ReactionInput<'_>],
    matches: &[Vec<usize>],
    selections: &[ReactionCoordinateSelection],
    set: usize,
) -> Result<Vec<ReactionProduct>, ReactionRunError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: generateOneProductSet
    // RDKit❗✔️: generateOneProductSet(const ChemicalReaction &rxn,
    // RDKit❗✔️:                       const MOL_SPTR_VECT &reactants,
    // RDKit❗✔️:                       const std::vector<MatchVectType> &reactantsMatch) {
    // RDKit❗✔️:   PRECONDITION(reactants.size() == reactantsMatch.size(),
    // RDKit❗✔️:                "vector size mismatch");
    // RDKit❗✔️:
    // RDKit❗✔️:   // if any of the reactants have a conformer, we'll go ahead and
    // RDKit❗✔️:   // generate conformers for the products:
    // RDKit❗✔️:   bool doConfs = false;
    // RDKit❗✔️:   // if any of the reactants have a single bond with directionality specified,
    // RDKit❗✔️:   // we will make sure that the output molecules have directionality
    // RDKit❗✔️:   // specified.
    // RDKit❗✔️:   bool doBondDirs = false;
    // RDKit❗✔️:   for (const auto &reactant : reactants) {
    // RDKit❗✔️:     if (reactant->getNumConformers()) {
    // RDKit❗✔️:       doConfs = true;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     for (const auto bnd : reactant->bonds()) {
    // RDKit❗✔️:       if (bnd->getBondType() == Bond::SINGLE &&
    // RDKit❗✔️:           bnd->getBondDir() > Bond::NONE) {
    // RDKit❗✔️:         doBondDirs = true;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (doConfs && doBondDirs) {
    // RDKit❗✔️:       break;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   MOL_SPTR_VECT res;
    // RDKit❗✔️:   res.resize(rxn.getNumProductTemplates());
    // RDKit❗✔️:   unsigned int prodId = 0;
    // RDKit❗✔️:   for (auto pTemplIt = rxn.beginProductTemplates();
    // RDKit❗✔️:        pTemplIt != rxn.endProductTemplates(); ++pTemplIt) {
    // RDKit❗✔️:     // copy product template and its properties to a new product RWMol
    // RDKit❗✔️:     RWMOL_SPTR product = convertTemplateToMol(*pTemplIt);
    // RDKit❗✔️:     Conformer *conf = nullptr;
    // RDKit❗✔️:     if (doConfs) {
    // RDKit❗✔️:       conf = new Conformer();
    // RDKit❗✔️:       conf->set3D(false);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     unsigned int reactantId = 0;
    // RDKit❗✔️:     for (auto iter = rxn.beginReactantTemplates();
    // RDKit❗✔️:          iter != rxn.endReactantTemplates(); ++iter, reactantId++) {
    // RDKit❗✔️:       addReactantAtomsAndBonds(rxn, product, reactants.at(reactantId),
    // RDKit❗✔️:                                reactantsMatch.at(reactantId), *iter, conf, reactantId);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     if (doConfs) {
    // RDKit❗✔️:       product->addConformer(conf, true);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // if there was bond direction information in any reactant, it has been
    // RDKit❗✔️:     // lost, add it back.
    // RDKit❗✔️:     if (doBondDirs) {
    // RDKit❗✔️:       MolOps::setDoubleBondNeighborDirections(*product);
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // if the product template has stereo groups, copy them over now
    // RDKit❗✔️:     if (!(*pTemplIt)->getStereoGroups().empty()) {
    // RDKit❗✔️:       copyTemplateStereoGroupsToMol(**pTemplIt, product);
    // RDKit❗✔️:     }
    // RDKit❗✔️:     product->updatePropertyCache(false);
    // RDKit❗✔️:     res[prodId] = product;
    // RDKit❗✔️:     ++prodId;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Ordered output allocation matches source; construction owns only fresh
    // product rows, one adjacency finalization and source-defined cache work.
    if inputs.len() != matches.len() || inputs.len() != selections.len() {
        return Err(ReactionRunError::ReactantArity {
            expected: inputs.len(),
            actual: matches.len(),
        });
    }
    let mut do_confs = false;
    let mut do_bond_dirs = false;
    for input in inputs {
        if !input.coordinates.conformers_2d.is_empty()
            || !input.coordinates.conformers_3d.is_empty()
        {
            do_confs = true;
        }
        if input
            .topology
            .bonds
            .iter()
            .any(|b| b.order() == BondOrder::Single && b.direction() != BondDirection::None)
        {
            do_bond_dirs = true;
        }
        if do_confs && do_bond_dirs {
            break;
        }
    }
    // Explicit absent selections must report D4's error even when no source
    // conformer exists. Resolution stays in reached product construction.
    let explicit = selections
        .iter()
        .any(|s| *s != ReactionCoordinateSelection::Auto);
    let mut result = Vec::with_capacity(reaction.num_product_templates());
    for (template_index, template) in reaction.product_templates().iter().enumerate() {
        let build = (|| -> Result<ReactionProduct, ReactionProductError> {
            let mut product = convert_template(template, template_index)?;
            let mut points = Vec::new();
            let mut is_3d = false;
            for (reactant_index, reactant_template) in
                reaction.reactant_templates().iter().enumerate()
            {
                let conformer = (do_confs || explicit).then_some((&mut points, &mut is_3d));
                add_reactant(
                    reaction,
                    &mut product,
                    inputs[reactant_index],
                    &matches[reactant_index],
                    reactant_template,
                    reactant_index,
                    conformer,
                    selections[reactant_index],
                )?;
            }
            let mut coordinates = CoordinateBlock::default();
            if do_confs {
                coordinates
                    .conformers_3d
                    .push(Conformer3D::new(0, points, is_3d));
            }
            // Canonical adjacency is finalized once, after ordered insertion.
            product.topology.adjacency = AdjacencyList::try_from_topology(
                product.topology.atoms.len(),
                &product.topology.bonds,
            )?;
            let mut properties = MoleculeProperties::default();
            let rings = if do_bond_dirs {
                // RDKit❗✔️: if (!mol.getRingInfo()->isSymmSssr()) {
                // RDKit❗✔️:   RDKit::MolOps::symmetrizeSSSR(mol);
                // RDKit❗✔️: }
                // Chirality.cpp source performs this preparation only after
                // generateOneProductSet's doBondDirs dispatch. No uncalled
                // sanitization, cache repair or rank assignment is introduced.
                let rings = cosmolkit_core::symmetrized_sssr(
                    &product.topology,
                    &cosmolkit_core::RingSearchParams::default(),
                )?;
                let update = cosmolkit_core::set_double_bond_neighbor_directions(
                    product.topology,
                    &rings,
                    None,
                )?;
                product.topology = update.topology;
                if update.needs_detect_bond_stereo {
                    // RDKit❗✔️:     mol.setProp("_needsDetectBondStereo", 1);
                    properties.set_prop("_needsDetectBondStereo", 1_i32)?;
                }
                Some(rings)
            } else {
                None
            };
            copy_template_groups(&mut product, template, template_index)?;
            let valence = cosmolkit_core::assign_valence_with_options_for_topology(
                &product.topology,
                ValenceModel::RdkitLike,
                false,
            )?;
            Ok(ReactionProduct {
                topology: product.topology,
                coordinates,
                properties,
                atom_origins: product.atom_origins,
                bond_origins: product.bond_origins,
                valence,
                rings,
            })
        })();
        result.push(build.map_err(|source| ReactionRunError::Product {
            set,
            template: template_index,
            source,
        })?);
    }
    Ok(result)
}

#[doc(hidden)]
pub fn run_reactants(
    reaction: &Reaction,
    inputs: &[ReactionInput<'_>],
    params: &ReactionRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactants
    // RDKit❗✔️: std::vector<MOL_SPTR_VECT> run_Reactants(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                          const MOL_SPTR_VECT &reactants,
    // RDKit❗✔️:                                          unsigned int maxProducts) {
    // RDKit❗✔️:   if (!rxn.isInitialized()) {
    // RDKit❗✔️:     throw ChemicalReactionException(
    // RDKit❗✔️:         "initMatchers() must be called before runReactants()");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (reactants.size() != rxn.getNumReactantTemplates()) {
    // RDKit❗✔️:     throw ChemicalReactionException(
    // RDKit❗✔️:         "Number of reactants provided does not match number of reactant "
    // RDKit❗✔️:         "templates.");
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (auto msptr : reactants) {
    // RDKit❗✔️:     CHECK_INVARIANT(msptr, "bad molecule in reactants");
    // RDKit❗✔️:     msptr->clearAllAtomBookmarks();  // we use this as scratch space
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   std::vector<MOL_SPTR_VECT> productMols;
    // RDKit❗✔️:   productMols.clear();
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we have no products, return now:
    // RDKit❗✔️:   if (!rxn.getNumProductTemplates()) {
    // RDKit❗✔️:     return productMols;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // find the matches for each reactant:
    // RDKit❗✔️:   VectVectMatchVectType matchesByReactant;
    // RDKit❗✔️:   if (!ReactionRunnerUtils::getReactantMatches(
    // RDKit❗✔️:           reactants, rxn, matchesByReactant, maxProducts)) {
    // RDKit❗✔️:     // some reactants didn't find a match, return an empty product list:
    // RDKit❗✔️:     return productMols;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   // -------------------------------------------------------
    // RDKit❗✔️:   // we now have matches for each reactant, so we can start creating products:
    // RDKit❗✔️:   // start by doing the combinatorics on the matches:
    // RDKit❗✔️:   VectVectMatchVectType reactantMatchesPerProduct;
    // RDKit❗✔️:   ReactionRunnerUtils::generateReactantCombinations(
    // RDKit❗✔️:       matchesByReactant, reactantMatchesPerProduct, maxProducts);
    // RDKit❗✔️:   productMols.resize(reactantMatchesPerProduct.size());
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int productId = 0; productId != productMols.size();
    // RDKit❗✔️:        ++productId) {
    // RDKit❗✔️:     MOL_SPTR_VECT lProds = ReactionRunnerUtils::generateOneProductSet(
    // RDKit❗✔️:         rxn, reactants, reactantMatchesPerProduct[productId]);
    // RDKit❗✔️:     productMols[productId] = lProds;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return productMols;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let reaction = initialized(reaction)?;
    if inputs.len() != reaction.num_reactant_templates() {
        return Err(ReactionRunError::ReactantArity {
            expected: reaction.num_reactant_templates(),
            actual: inputs.len(),
        });
    }
    if !params.coordinate_selections.is_empty()
        && params.coordinate_selections.len() != inputs.len()
    {
        return Err(ReactionRunError::CoordinateSelectionArity {
            expected: inputs.len(),
            actual: params.coordinate_selections.len(),
        });
    }
    if reaction.num_product_templates() == 0 {
        return Ok(Vec::new());
    }
    let Some(matches) =
        crate::matching::reactant_matches(inputs, &reaction, params.max_products, None)?
    else {
        return Ok(Vec::new());
    };
    let combinations = crate::matching::reactant_combinations(&matches, params.max_products)?;
    let default_selections;
    let selections = if params.coordinate_selections.is_empty() {
        default_selections = vec![ReactionCoordinateSelection::Auto; inputs.len()];
        &default_selections
    } else {
        &params.coordinate_selections
    };
    let mut products = Vec::with_capacity(combinations.len());
    for (set, combination) in combinations.iter().enumerate() {
        products.push(one_product_set(
            &reaction,
            inputs,
            combination,
            selections,
            set,
        )?);
    }
    Ok(products)
}

#[doc(hidden)]
pub fn run_reactant(
    reaction: &Reaction,
    input: ReactionInput<'_>,
    reactant_template: usize,
    params: &ReactionSingleRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactant
    // RDKit❗✔️: std::vector<MOL_SPTR_VECT> run_Reactant(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                         const ROMOL_SPTR &reactant,
    // RDKit❗✔️:                                         unsigned int reactantIdx) {
    // RDKit❗✔️:   if (!rxn.isInitialized()) {
    // RDKit❗✔️:     throw ChemicalReactionException(
    // RDKit❗✔️:         "initMatchers() must be called before runReactants()");
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   PRECONDITION(reactant, "bad molecule in reactants");
    // RDKit❗✔️:   reactant->clearAllAtomBookmarks();  // we use this as scratch space
    // RDKit❗✔️:   std::vector<MOL_SPTR_VECT> productMols;
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we have no products, return now:
    // RDKit❗✔️:   if (!rxn.getNumProductTemplates()) {
    // RDKit❗✔️:     return productMols;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   PRECONDITION(static_cast<size_t>(reactantIdx) < rxn.getReactants().size(),
    // RDKit❗✔️:                "reactantIdx out of bounds");
    // RDKit❗✔️:   // find the matches for each reactant:
    // RDKit❗✔️:   VectVectMatchVectType matchesByReactant;
    // RDKit❗✔️:
    // RDKit❗✔️:   // assemble the reactants (use an empty mol for missing reactants)
    // RDKit❗✔️:   MOL_SPTR_VECT reactants(rxn.getNumReactantTemplates());
    // RDKit❗✔️:   for (size_t i = 0; i < rxn.getNumReactantTemplates(); ++i) {
    // RDKit❗✔️:     if (i == reactantIdx) {
    // RDKit❗✔️:       reactants[i] = reactant;
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       reactants[i] = ROMOL_SPTR(new ROMol);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!ReactionRunnerUtils::getReactantMatches(
    // RDKit❗✔️:           reactants, rxn, matchesByReactant, 1000, reactantIdx)) {
    // RDKit❗✔️:     return productMols;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   VectMatchVectType &matches = matchesByReactant[reactantIdx];
    // RDKit❗✔️:   // each match on a reactant is a separate product
    // RDKit❗✔️:   VectVectMatchVectType matchesAtReactants(matches.size());
    // RDKit❗✔️:   for (size_t i = 0; i < matches.size(); ++i) {
    // RDKit❗✔️:     matchesAtReactants[i].resize(rxn.getReactants().size());
    // RDKit❗✔️:     matchesAtReactants[i][reactantIdx] = matches[i];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   productMols.resize(matches.size());
    // RDKit❗✔️:
    // RDKit❗✔️:   for (unsigned int productId = 0; productId != productMols.size();
    // RDKit❗✔️:        ++productId) {
    // RDKit❗✔️:     MOL_SPTR_VECT lProds = ReactionRunnerUtils::generateOneProductSet(
    // RDKit❗✔️:         rxn, reactants, matchesAtReactants[productId]);
    // RDKit❗✔️:     productMols[productId] = lProds;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return productMols;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let reaction = initialized(reaction)?;
    // SOURCE returns for zero products before validating the template index.
    if reaction.num_product_templates() == 0 {
        return Ok(Vec::new());
    }
    if reactant_template >= reaction.num_reactant_templates() {
        return Err(ReactionRunError::ReactantTemplateIndex {
            index: reactant_template,
            count: reaction.num_reactant_templates(),
        });
    }
    // Source unselected slots are empty detached graphs, with no Molecule or
    // synthetic runtime/capabilities and no fabricated cache facts.
    let empty_topology = TopologyBlock::default();
    let empty_coordinates = CoordinateBlock::default();
    let empty_properties = MoleculeProperties::default();
    let empty = ReactionInput {
        topology: &empty_topology,
        coordinates: &empty_coordinates,
        properties: &empty_properties,
        rings: None,
        valence: None,
    };
    let mut inputs = vec![empty; reaction.num_reactant_templates()];
    inputs[reactant_template] = input;
    let Some(mut matches) =
        crate::matching::reactant_matches(&inputs, &reaction, 1000, Some(reactant_template))?
    else {
        return Ok(Vec::new());
    };
    let selected = std::mem::take(&mut matches[reactant_template]);
    let mut selections = vec![ReactionCoordinateSelection::Auto; inputs.len()];
    selections[reactant_template] = params.coordinate_selection;
    let mut result = Vec::with_capacity(selected.len());
    for (set, selected_match) in selected.into_iter().enumerate() {
        let mut combination = vec![Vec::new(); inputs.len()];
        combination[reactant_template] = selected_match;
        let mut products = one_product_set(&reaction, &inputs, &combination, &selections, set)?;
        for (template, product) in products.iter_mut().enumerate() {
            // SOURCE react_idx remains its reaction-template index. Private
            // detached origins instead address the sole actual runtime input.
            for origin in product.atom_origins.iter_mut().flatten() {
                if origin.input != reactant_template {
                    return Err(ReactionRunError::Product {
                        set,
                        template,
                        source: invariant(
                            "run_Reactant",
                            "row came from an unselected empty reagent",
                            Some(origin.row.index()),
                            None,
                            None,
                        ),
                    });
                }
                origin.input = 0;
            }
            for origin in product.bond_origins.iter_mut().flatten() {
                if origin.input != reactant_template {
                    return Err(ReactionRunError::Product {
                        set,
                        template,
                        source: invariant(
                            "run_Reactant",
                            "bond came from an unselected empty reagent",
                            None,
                            None,
                            Some(origin.row),
                        ),
                    });
                }
                origin.input = 0;
            }
        }
        result.push(products);
    }
    Ok(result)
}
