use crate::materialize::{
    ProductBuilder, add_neighbors, convert_template, invariant, transfer_atom_properties,
    transfer_bond_properties,
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
pub(crate) fn initialize_for_run(reaction: &mut Reaction) -> Result<(), ReactionRunError> {
    // RDKit❗✔️:   if (!self->isInitialized()) {
    // RDKit❗✔️:     NOGIL gil;
    // RDKit❗✔️:     self->initReactantMatchers();
    // RDKit❗✔️:   }
    // Same-instance source initialization preserves validation mutation/error
    // prefixes. The source's false validation result returns normally; the
    // reached runner retains its own initialization precondition and order.
    // Cost: reuse the unique initializer, with no Reaction/template copy.
    if !reaction.is_initialized() {
        crate::validation::init_reactant_matchers_source(
            reaction,
            &ReactionValidationParams::default(),
            &mut crate::ReactionValidationReport::default(),
        )
        .map_err(|source| crate::ReactionInitializationError::Validation { source })?;
    }
    Ok(())
}

fn add_reactant(
    reaction: &Reaction,
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    matches: &[usize],
    template: &QueryGraph,
    template_index: usize,
    conformer: Option<(&mut Vec<[f64; 3]>, &mut bool)>,
    selection: ReactionCoordinateSelection,
) -> Result<(), ReactionProductError> {
    // Canonical matcher output is indexed by query row. Preserve that source
    // pair sequence without sorting, deduplication, or another chemistry path.
    let matched = matches.iter().enumerate().map(|(query, &target)| {
        Ok((
            crate::materialize::source_u32("matched query row", query)? as i32,
            crate::materialize::source_u32("matched reactant row", target)? as i32,
        ))
    });
    add_reactant_source(
        reaction,
        product,
        input,
        matched,
        template,
        template_index,
        conformer,
        selection,
    )
}

fn add_reactant_source(
    reaction: &Reaction,
    product: &mut ProductBuilder,
    input: ReactionInput<'_>,
    matched: impl IntoIterator<Item = Result<(i32, i32), ReactionProductError>> + Clone,
    template: &QueryGraph,
    template_index: usize,
    mut conformer: Option<(&mut Vec<[f64; 3]>, &mut bool)>,
    selection: ReactionCoordinateSelection,
) -> Result<(), ReactionProductError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: addReactantAtomsAndBonds
    // RDKit❗❌: void addReactantAtomsAndBonds(const ChemicalReaction &rxn, RWMOL_SPTR product,
    // RDKit❗❌:                               const ROMOL_SPTR reactantSptr,
    // RDKit❗❌:                               const MatchVectType &match,
    // RDKit❗❌:                               const ROMOL_SPTR reactantTemplate,
    // RDKit❗❌:                               Conformer *productConf,
    // RDKit❗❌: 			      unsigned int reactantId) {
    // RDKit❗❌:   // start by looping over all matches and marking the reactant atoms that
    // RDKit❗❌:   // have already been "added" by virtue of being in the product. We'll also
    // RDKit❗❌:   // mark "skipped" atoms: those that are in the match, but not in this
    // RDKit❗❌:   // particular product (or, perhaps, not in any product)
    // RDKit❗❌:   // At the same time we'll set up a map between the indices of those
    // RDKit❗❌:   // atoms and their index in the product.
    // RDKit❗❌:   ReactantProductAtomMapping *mapping = getAtomMappingsReactantProduct(
    // RDKit❗❌:       match, *reactantTemplate, product, reactantSptr->getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:   boost::dynamic_bitset<> visitedAtoms(reactantSptr->getNumAtoms());
    // RDKit❗❌:
    // RDKit❗❌:   const ROMol *reactant = reactantSptr.get();
    // RDKit❗❌:
    // RDKit❗❌:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗❌:   // Loop over the bonds in the product and look for those that have
    // RDKit❗❌:   // the NullBond property set. These are bonds for which no information
    // RDKit❗❌:   // (other than their existence) was provided in the template
    // RDKit❗❌:   setReactantBondPropertiesToProduct(product, *reactant, mapping);
    // RDKit❗❌:
    // RDKit❗❌:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗❌:   // Loop over the atoms in the match that were added to the product
    // RDKit❗❌:   // From the corresponding atom in the reactant, do a graph traversal
    // RDKit❗❌:   // to find other connected atoms that should be added:
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<const Atom *> chiralAtomsToCheck;
    // RDKit❗❌:   for (const auto &matchIdx : match) {
    // RDKit❗❌:     int reactantAtomIdx = matchIdx.second;
    // RDKit❗❌:     if (mapping->mappedAtoms[reactantAtomIdx]) {
    // RDKit❗❌:       CHECK_INVARIANT(mapping->reactProdAtomMap.find(reactantAtomIdx) !=
    // RDKit❗❌:                           mapping->reactProdAtomMap.end(),
    // RDKit❗❌:                       "mapped reactant atom not present in product.");
    // RDKit❗❌:
    // RDKit❗❌:       const Atom *reactantAtom = reactant->getAtomWithIdx(reactantAtomIdx);
    // RDKit❗❌:       for (unsigned i = 0;
    // RDKit❗❌:            i < mapping->reactProdAtomMap[reactantAtomIdx].size(); i++) {
    // RDKit❗❌:         // here's a pointer to the atom in the product:
    // RDKit❗❌:         unsigned productAtomIdx = mapping->reactProdAtomMap[reactantAtomIdx][i];
    // RDKit❗❌:         Atom *productAtom = product->getAtomWithIdx(productAtomIdx);
    // RDKit❗❌:         setReactantAtomPropertiesToProduct(productAtom, *reactantAtom,
    // RDKit❗❌:                                            rxn.getImplicitPropertiesFlag(), reactantId);
    // RDKit❌❌:         if (reactantAtom->hasQuery()) {
    // RDKit❌❌:           // finally: if the reactant atom is a query we should copy over the
    // RDKit❌❌:           // query information. We need to replace the atom to do this
    // RDKit❌❌:           QueryAtom newAtom(*productAtom);
    // RDKit❌❌:           newAtom.setQuery(reactantAtom->getQuery()->copy());
    // RDKit❌❌:           // replaceAtom copies
    // RDKit❌❌:           product->replaceAtom(productAtomIdx, &newAtom);
    // RDKit❌❌:         }
    // RDKit❗❌:       }
    // RDKit❗❌:       // now traverse:
    // RDKit❗❌:       addReactantNeighborsToProduct(*reactant, *reactantAtom, product,
    // RDKit❗❌:                                     visitedAtoms, chiralAtomsToCheck, mapping, reactantId);
    // RDKit❗❌:
    // RDKit❗❌:       // now that we've added all the reactant's neighbors, check to see if
    // RDKit❗❌:       // it is chiral in the reactant but is not in the reaction. If so
    // RDKit❗❌:       // we need to worry about its chirality
    // RDKit❗❌:       checkAndCorrectChiralityOfMatchingAtomsInProduct(
    // RDKit❗❌:           *reactant, reactantAtomIdx, *reactantAtom, product, mapping);
    // RDKit❗❌:     }
    // RDKit❗❌:   }  // end of loop over matched atoms
    // RDKit❗❌:
    // RDKit❗❌:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗❌:   // now we need to loop over atoms from the reactants that were chiral but
    // RDKit❗❌:   // not directly involved in the reaction in order to make sure their
    // RDKit❗❌:   // chirality hasn't been disturbed
    // RDKit❗❌:   checkAndCorrectChiralityOfProduct(chiralAtomsToCheck, product, mapping);
    // RDKit❗❌:
    // RDKit❗❌:   updateStereoBonds(product, *reactant, mapping);
    // RDKit❗❌:
    // RDKit❗❌:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗❌:   // Copy enhanced StereoGroup data from reactant to product if it is
    // RDKit❗❌:   // still valid. Uses ChiralTag checks above.
    // RDKit❗❌:   copyEnhancedStereoGroups(*reactant, product, *mapping);
    // RDKit❗❌:
    // RDKit❗❌:   // ---------- ---------- ---------- ---------- ---------- ----------
    // RDKit❗❌:   // finally we may need to set the coordinates in the product conformer:
    // RDKit❗❌:   if (productConf) {
    // RDKit❗❌:     productConf->resize(product->getNumAtoms());
    // RDKit❗❌:     generateProductConformers(productConf, *reactant, mapping);
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   delete (mapping);
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Native signed pair order is preserved without allocating a second match
    // vector. Clone only the lazy iterator for the second source traversal.
    // Cost gap: the visited and mapping flags use byte Vec<bool> instead of
    // packed dynamic_bitset; detached provenance and checked row reads add work.
    let mut mapping = crate::materialize::atom_mappings_source(
        matched.clone(),
        template,
        template_index,
        product,
        crate::materialize::source_u32("reactant atom count", input.topology.atoms.len())?,
    )?;
    let mut visited = vec![false; input.topology.atoms.len()];
    transfer_bond_properties(product, input, template_index, &mapping)?;
    let mut chiral_to_check = Vec::new();
    for matched in matched {
        let (_, reactant_signed) = matched?;
        // Native bitset indexes use size_t, while getAtomWithIdx and mapping
        // keys convert the same source int to unsigned int.
        let bit = reactant_signed as usize;
        let reactant_atom = (reactant_signed as u32) as usize;
        if *mapping.mapped.get(bit).ok_or_else(|| {
            invariant(
                "addReactantAtomsAndBonds",
                "mapped bit index out of range",
                Some(bit),
                None,
                None,
            )
        })? {
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
            let r = input.topology.atoms.get(reactant_atom).ok_or_else(|| {
                invariant(
                    "addReactantAtomsAndBonds",
                    "matched reactant atom row missing",
                    Some(reactant_atom),
                    None,
                    None,
                )
            })?;
            let mut copy = 0u32;
            while (copy as usize) < rows.len() {
                let row = rows[copy as usize];
                let product_atom = product.topology.atoms.get_mut(row).ok_or_else(|| {
                    invariant(
                        "addReactantAtomsAndBonds",
                        "mapped product atom row missing",
                        Some(reactant_atom),
                        Some(row),
                        None,
                    )
                })?;
                transfer_atom_properties(
                    product_atom,
                    r,
                    reaction.implicit_properties(),
                    template_index,
                )?;
                *product.atom_origins.get_mut(row).ok_or_else(|| {
                    invariant(
                        "addReactantAtomsAndBonds",
                        "mapped product provenance row missing",
                        Some(reactant_atom),
                        Some(row),
                        None,
                    )
                })? = Some(crate::ReactionRowOrigin {
                    input: template_index,
                    row: AtomId::new(reactant_atom),
                });
                copy = copy.wrapping_add(1);
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
            correct_matched_chirality(product, &input, &mut mapping, reactant_atom)?;
        }
    }
    correct_unmatched_chirality(product, &input, &mut mapping, &chiral_to_check)?;
    update_stereo_bonds(product, &input, &mut mapping)?;
    copy_enhanced_groups(product, &input, &mapping)?;
    if let Some((points, is_3d)) = conformer.as_mut() {
        points.resize(product.topology.atoms.len(), [0.0; 3]);
        propagate_coordinates(points, is_3d, &input, &mapping, selection)?;
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
    // RDKit❗❌: generateOneProductSet(const ChemicalReaction &rxn,
    // RDKit❗❌:                       const MOL_SPTR_VECT &reactants,
    // RDKit❗❌:                       const std::vector<MatchVectType> &reactantsMatch) {
    // RDKit❗❌:   PRECONDITION(reactants.size() == reactantsMatch.size(),
    // RDKit❗❌:                "vector size mismatch");
    // RDKit❗❌:
    // RDKit❗❌:   // if any of the reactants have a conformer, we'll go ahead and
    // RDKit❗❌:   // generate conformers for the products:
    // RDKit❗❌:   bool doConfs = false;
    // RDKit❗❌:   // if any of the reactants have a single bond with directionality specified,
    // RDKit❗❌:   // we will make sure that the output molecules have directionality
    // RDKit❗❌:   // specified.
    // RDKit❗❌:   bool doBondDirs = false;
    // RDKit❗❌:   for (const auto &reactant : reactants) {
    // RDKit❗❌:     if (reactant->getNumConformers()) {
    // RDKit❗❌:       doConfs = true;
    // RDKit❗❌:     }
    // RDKit❗❌:     for (const auto bnd : reactant->bonds()) {
    // RDKit❗❌:       if (bnd->getBondType() == Bond::SINGLE &&
    // RDKit❗❌:           bnd->getBondDir() > Bond::NONE) {
    // RDKit❗❌:         doBondDirs = true;
    // RDKit❗❌:         break;
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:     if (doConfs && doBondDirs) {
    // RDKit❗❌:       break;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   MOL_SPTR_VECT res;
    // RDKit❗❌:   res.resize(rxn.getNumProductTemplates());
    // RDKit❗❌:   unsigned int prodId = 0;
    // RDKit❗❌:   for (auto pTemplIt = rxn.beginProductTemplates();
    // RDKit❗❌:        pTemplIt != rxn.endProductTemplates(); ++pTemplIt) {
    // RDKit❗❌:     // copy product template and its properties to a new product RWMol
    // RDKit❗❌:     RWMOL_SPTR product = convertTemplateToMol(*pTemplIt);
    // RDKit❗❌:     Conformer *conf = nullptr;
    // RDKit❗❌:     if (doConfs) {
    // RDKit❗❌:       conf = new Conformer();
    // RDKit❗❌:       conf->set3D(false);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     unsigned int reactantId = 0;
    // RDKit❗❌:     for (auto iter = rxn.beginReactantTemplates();
    // RDKit❗❌:          iter != rxn.endReactantTemplates(); ++iter, reactantId++) {
    // RDKit❗❌:       addReactantAtomsAndBonds(rxn, product, reactants.at(reactantId),
    // RDKit❗❌:                                reactantsMatch.at(reactantId), *iter, conf, reactantId);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     if (doConfs) {
    // RDKit❗❌:       product->addConformer(conf, true);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // if there was bond direction information in any reactant, it has been
    // RDKit❗❌:     // lost, add it back.
    // RDKit❗❌:     if (doBondDirs) {
    // RDKit❗❌:       MolOps::setDoubleBondNeighborDirections(*product);
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     // if the product template has stereo groups, copy them over now
    // RDKit❗❌:     if (!(*pTemplIt)->getStereoGroups().empty()) {
    // RDKit❗❌:       copyTemplateStereoGroupsToMol(**pTemplIt, product);
    // RDKit❗❌:     }
    // RDKit❗❌:     product->updatePropertyCache(false);
    // RDKit❗❌:     res[prodId] = product;
    // RDKit❗❌:     ++prodId;
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Source ordered product construction is retained. Detached adjacency
    // finalization and returned cache vectors allocate O(V+E)/O(V) extra state
    // relative to source graph/cache mutation; inherited helper costs remain ❌.
    if inputs.len() != matches.len() {
        return Err(ReactionRunError::ReactantArity {
            expected: inputs.len(),
            actual: matches.len(),
        });
    }
    if inputs.len() != selections.len() {
        return Err(ReactionRunError::CoordinateSelectionArity {
            expected: inputs.len(),
            actual: selections.len(),
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
    let mut result = Vec::with_capacity(reaction.num_product_templates());
    for (template_index, template) in reaction.product_templates().iter().enumerate() {
        let build = (|| -> Result<ReactionProduct, ReactionProductError> {
            let mut product = convert_template(template, template_index)?;
            let mut points = Vec::new();
            let mut is_3d = false;
            for (reactant_index, reactant_template) in
                reaction.reactant_templates().iter().enumerate()
            {
                // Source vector::at is reached after product conversion, not
                // a whole-reaction arity/graph validation at function entry.
                let input = *inputs.get(reactant_index).ok_or_else(|| {
                    invariant(
                        "generateOneProductSet",
                        "reactant input index out of range",
                        Some(reactant_index),
                        None,
                        None,
                    )
                })?;
                let matched = matches.get(reactant_index).ok_or_else(|| {
                    invariant(
                        "generateOneProductSet",
                        "reactant match index out of range",
                        Some(reactant_index),
                        None,
                        None,
                    )
                })?;
                let selection = selections[reactant_index];
                let conformer = do_confs.then_some((&mut points, &mut is_3d));
                add_reactant(
                    reaction,
                    &mut product,
                    input,
                    matched,
                    reactant_template,
                    reactant_index,
                    conformer,
                    selection,
                )?;
            }
            let mut coordinates = CoordinateBlock::default();
            if do_confs {
                // RDKit❗✔️: unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
                // RDKit❗✔️:   PRECONDITION(conf, "bad conformer");
                // RDKit❗✔️:   PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
                // RDKit❗✔️:                "Number of atom mismatch");
                // RDKit❗✔️:   if (assignId) {
                // RDKit❗✔️:     int maxId = -1;
                // RDKit❗✔️:     for (auto cptr : d_confs) {
                // RDKit❗✔️:       maxId = std::max((int)(cptr->getId()), maxId);
                // RDKit❗✔️:     }
                // RDKit❗✔️:     maxId++;
                // RDKit❗✔️:     conf->setId((unsigned int)maxId);
                // RDKit❗✔️:   }
                // RDKit❗✔️:   conf->setOwningMol(this);
                // RDKit❗✔️:   CONFORMER_SPTR nConf(conf);
                // RDKit❗✔️:   d_confs.push_back(nConf);
                // RDKit❗✔️:   return conf->getId();
                // RDKit❗✔️: }
                // convertTemplateToMol creates a fresh molecule with no source
                // conformers. The exact assignId loop therefore assigns zero.
                let conformer = Conformer3D::new(0, points, is_3d);
                conformer.validate_for_atom_count(product.topology.atoms.len())?;
                coordinates.conformers_3d.push(conformer);
            }
            // Canonical adjacency is finalized once, after ordered insertion.
            product.topology.adjacency = AdjacencyList::try_from_topology(
                product.topology.atoms.len(),
                &product.topology.bonds,
            )?;
            let mut properties = MoleculeProperties::default();
            let rings = if do_bond_dirs {
                // RDKit❗❌: if (!mol.getRingInfo()->isSymmSssr()) {
                // RDKit❗❌:   RDKit::MolOps::symmetrizeSSSR(mol);
                // RDKit❗❌: }
                // Chirality.cpp source performs this preparation only after
                // generateOneProductSet's doBondDirs dispatch. No uncalled
                // sanitization, cache repair or rank assignment is introduced.
                let (rings, prepared_properties) =
                    cosmolkit_core::symmetrized_sssr_with_properties(
                        &product.topology,
                        properties,
                        &cosmolkit_core::RingSearchParams::default(),
                    )?;
                properties = prepared_properties;
                let update = cosmolkit_core::set_double_bond_neighbor_directions(
                    product.topology,
                    &rings,
                    None,
                )?;
                product.topology = update.topology;
                if update.needs_detect_bond_stereo {
                    // RDKit❗❌:     mol.setProp("_needsDetectBondStereo", 1);
                    properties.set_prop("_needsDetectBondStereo", 1_i32)?;
                }
                Some(rings)
            } else {
                None
            };
            if !template.stereo_groups().is_empty() {
                copy_template_groups(&mut product, template, template_index)?;
            }
            // Use the sole source updatePropertyCache traversal directly on
            // borrowed rows; no unrelated whole-topology validation is added.
            let valence = cosmolkit_core::assign_valence_with_options_from_parts(
                &product.topology.atoms,
                &product.topology.bonds,
                &product.topology.adjacency,
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
    reaction: &mut Reaction,
    inputs: &[ReactionInput<'_>],
    params: &ReactionRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    initialize_for_run(reaction)?;
    run_reactants_source(reaction, inputs, params)
}

fn run_reactants_source(
    reaction: &Reaction,
    inputs: &[ReactionInput<'_>],
    params: &ReactionRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactants
    // RDKit❗❌: std::vector<MOL_SPTR_VECT> run_Reactants(const ChemicalReaction &rxn,
    // RDKit❗❌:                                          const MOL_SPTR_VECT &reactants,
    // RDKit❗❌:                                          unsigned int maxProducts) {
    // RDKit❗❌:   if (!rxn.isInitialized()) {
    // RDKit❗❌:     throw ChemicalReactionException(
    // RDKit❗❌:         "initMatchers() must be called before runReactants()");
    // RDKit❗❌:   }
    // RDKit❗❌:   if (reactants.size() != rxn.getNumReactantTemplates()) {
    // RDKit❗❌:     throw ChemicalReactionException(
    // RDKit❗❌:         "Number of reactants provided does not match number of reactant "
    // RDKit❗❌:         "templates.");
    // RDKit❗❌:   }
    // RDKit❗❌:   for (auto msptr : reactants) {
    // RDKit❗❌:     CHECK_INVARIANT(msptr, "bad molecule in reactants");
    // RDKit❗❌:     msptr->clearAllAtomBookmarks();  // we use this as scratch space
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::vector<MOL_SPTR_VECT> productMols;
    // RDKit❗❌:   productMols.clear();
    // RDKit❗❌:
    // RDKit❗❌:   // if we have no products, return now:
    // RDKit❗❌:   if (!rxn.getNumProductTemplates()) {
    // RDKit❗❌:     return productMols;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   // find the matches for each reactant:
    // RDKit❗❌:   VectVectMatchVectType matchesByReactant;
    // RDKit❗❌:   if (!ReactionRunnerUtils::getReactantMatches(
    // RDKit❗❌:           reactants, rxn, matchesByReactant, maxProducts)) {
    // RDKit❗❌:     // some reactants didn't find a match, return an empty product list:
    // RDKit❗❌:     return productMols;
    // RDKit❗❌:   }
    // RDKit❗❌:   // -------------------------------------------------------
    // RDKit❗❌:   // we now have matches for each reactant, so we can start creating products:
    // RDKit❗❌:   // start by doing the combinatorics on the matches:
    // RDKit❗❌:   VectVectMatchVectType reactantMatchesPerProduct;
    // RDKit❗❌:   ReactionRunnerUtils::generateReactantCombinations(
    // RDKit❗❌:       matchesByReactant, reactantMatchesPerProduct, maxProducts);
    // RDKit❗❌:   productMols.resize(reactantMatchesPerProduct.size());
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int productId = 0; productId != productMols.size();
    // RDKit❗❌:        ++productId) {
    // RDKit❗❌:     MOL_SPTR_VECT lProds = ReactionRunnerUtils::generateOneProductSet(
    // RDKit❗❌:         rxn, reactants, reactantMatchesPerProduct[productId]);
    // RDKit❗❌:     productMols[productId] = lProds;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return productMols;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    if !reaction.is_initialized() {
        return Err(ReactionRunError::NeedsInitialization);
    }
    // Borrowed inputs are non-null and have no atom-bookmark scratch carrier.
    // All source reaction scratch bookmarks live only in each fresh private
    // ProductBuilder; this function never mutates caller input/model state.
    // Keep the explicit source owner/bookmark projection gap in final review.
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
        crate::matching::reactant_matches(inputs, reaction, params.max_products, None)?
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
    // Source productId is uint32, including its native wrap behavior.
    let mut product_id = 0u32;
    while (product_id as usize) != combinations.len() {
        let set = product_id as usize;
        let combination = &combinations[set];
        let mut product_set = one_product_set(reaction, inputs, combination, selections, set)?;
        if params.copy_atom_properties {
            for (template, product) in product_set.iter_mut().enumerate() {
                copy_user_atom_properties(inputs, product).map_err(|source| {
                    ReactionRunError::Product {
                        set,
                        template,
                        source,
                    }
                })?;
            }
        }
        products.push(product_set);
        product_id = product_id.wrapping_add(1);
    }
    Ok(products)
}

/// CK opt-in metadata propagation, not an alternative reaction algorithm.
/// Existing typed origins include the reagent index and all duplicated rows.
/// Only missing user keys are filled, so explicit template values always win.
fn copy_user_atom_properties(
    inputs: &[ReactionInput<'_>],
    product: &mut ReactionProduct,
) -> Result<(), ReactionProductError> {
    if product.atom_origins.len() != product.topology.atoms.len() {
        return Err(invariant(
            "copy_atom_properties",
            "origin count differs from atom count",
            None,
            None,
            None,
        ));
    }
    for (destination, origin) in product.topology.atoms.iter_mut().zip(&product.atom_origins) {
        let Some(origin) = origin else { continue };
        let source = inputs
            .get(origin.input)
            .and_then(|input| input.topology.atoms.get(origin.row.index()))
            .ok_or_else(|| {
                invariant(
                    "copy_atom_properties",
                    "source origin is out of range",
                    Some(origin.row.index()),
                    Some(destination.id().index()),
                    None,
                )
            })?;
        for (key, value) in source.property_records(false, false)? {
            if cosmolkit_model::is_user_atom_property(key.as_bytes())
                && destination.prop(key).is_none()
            {
                destination.set_prop(key.clone(), value)?;
            }
        }
    }
    Ok(())
}

#[doc(hidden)]
pub fn run_reactant(
    reaction: &mut Reaction,
    input: ReactionInput<'_>,
    reactant_template: usize,
    params: &ReactionSingleRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    initialize_for_run(reaction)?;
    let mut result = run_reactant_source(reaction, input, reactant_template, params)?;
    // Detached output addresses the sole actual input. Keep source react_idx
    // atom properties unchanged; only private runtime origins are projected.
    for (set, products) in result.iter_mut().enumerate() {
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
    }
    Ok(result)
}

fn run_reactant_source(
    reaction: &Reaction,
    input: ReactionInput<'_>,
    reactant_template: usize,
    params: &ReactionSingleRunParams,
) -> Result<Vec<Vec<ReactionProduct>>, ReactionRunError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionRunner.cpp :: run_Reactant
    // RDKit❗❌: std::vector<MOL_SPTR_VECT> run_Reactant(const ChemicalReaction &rxn,
    // RDKit❗❌:                                         const ROMOL_SPTR &reactant,
    // RDKit❗❌:                                         unsigned int reactantIdx) {
    // RDKit❗❌:   if (!rxn.isInitialized()) {
    // RDKit❗❌:     throw ChemicalReactionException(
    // RDKit❗❌:         "initMatchers() must be called before runReactants()");
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   PRECONDITION(reactant, "bad molecule in reactants");
    // RDKit❗❌:   reactant->clearAllAtomBookmarks();  // we use this as scratch space
    // RDKit❗❌:   std::vector<MOL_SPTR_VECT> productMols;
    // RDKit❗❌:
    // RDKit❗❌:   // if we have no products, return now:
    // RDKit❗❌:   if (!rxn.getNumProductTemplates()) {
    // RDKit❗❌:     return productMols;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   PRECONDITION(static_cast<size_t>(reactantIdx) < rxn.getReactants().size(),
    // RDKit❗❌:                "reactantIdx out of bounds");
    // RDKit❗❌:   // find the matches for each reactant:
    // RDKit❗❌:   VectVectMatchVectType matchesByReactant;
    // RDKit❗❌:
    // RDKit❗❌:   // assemble the reactants (use an empty mol for missing reactants)
    // RDKit❗❌:   MOL_SPTR_VECT reactants(rxn.getNumReactantTemplates());
    // RDKit❗❌:   for (size_t i = 0; i < rxn.getNumReactantTemplates(); ++i) {
    // RDKit❗❌:     if (i == reactantIdx) {
    // RDKit❗❌:       reactants[i] = reactant;
    // RDKit❗❌:     } else {
    // RDKit❗❌:       reactants[i] = ROMOL_SPTR(new ROMol);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   if (!ReactionRunnerUtils::getReactantMatches(
    // RDKit❗❌:           reactants, rxn, matchesByReactant, 1000, reactantIdx)) {
    // RDKit❗❌:     return productMols;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   VectMatchVectType &matches = matchesByReactant[reactantIdx];
    // RDKit❗❌:   // each match on a reactant is a separate product
    // RDKit❗❌:   VectVectMatchVectType matchesAtReactants(matches.size());
    // RDKit❗❌:   for (size_t i = 0; i < matches.size(); ++i) {
    // RDKit❗❌:     matchesAtReactants[i].resize(rxn.getReactants().size());
    // RDKit❗❌:     matchesAtReactants[i][reactantIdx] = matches[i];
    // RDKit❗❌:   }
    // RDKit❗❌:   productMols.resize(matches.size());
    // RDKit❗❌:
    // RDKit❗❌:   for (unsigned int productId = 0; productId != productMols.size();
    // RDKit❗❌:        ++productId) {
    // RDKit❗❌:     MOL_SPTR_VECT lProds = ReactionRunnerUtils::generateOneProductSet(
    // RDKit❗❌:         rxn, reactants, matchesAtReactants[productId]);
    // RDKit❗❌:     productMols[productId] = lProds;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return productMols;
    // RDKit❗❌: }
    // END RDKIT COMPLETE CPP FUNCTION
    if !reaction.is_initialized() {
        return Err(ReactionRunError::NeedsInitialization);
    }
    // Source input pointer is non-null by the borrowed detached type. Native
    // input atom-bookmark scratch is absent; each product owns private scratch.
    // This owner/side-effect projection remains explicit final comparison.
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
        crate::matching::reactant_matches(&inputs, reaction, 1000, Some(reactant_template))?
    else {
        return Ok(Vec::new());
    };
    let selected = std::mem::take(&mut matches[reactant_template]);
    let mut selections = vec![ReactionCoordinateSelection::Auto; inputs.len()];
    selections[reactant_template] = params.coordinate_selection;
    // Stage all per-match vectors before building products, as source does.
    let mut matches_at_reactants = Vec::with_capacity(selected.len());
    for selected_match in &selected {
        let mut combination = vec![Vec::new(); inputs.len()];
        combination[reactant_template] = selected_match.clone();
        matches_at_reactants.push(combination);
    }
    let mut result = vec![Vec::new(); selected.len()];
    let mut product_id = 0u32;
    while (product_id as usize) != result.len() {
        let set = product_id as usize;
        result[set] = one_product_set(
            reaction,
            &inputs,
            &matches_at_reactants[set],
            &selections,
            set,
        )?;
        product_id = product_id.wrapping_add(1);
    }
    Ok(result)
}

#[cfg(test)]
mod complete_add_reactant_source_tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomSpec, Bond, BondId, BondSpec, PropertyValue, QueryAtom, StereoGroup,
        StereoGroupKind,
    };
    use cosmolkit_types::{ChiralTag, Element};
    use std::collections::BTreeMap;
    pub(super) fn graph(maps: &[Option<u32>]) -> QueryGraph {
        QueryGraph::from_parts(
            maps.iter()
                .enumerate()
                .map(|(i, map)| {
                    let mut spec = AtomSpec::new(Element::C);
                    if let Some(map) = map {
                        spec = spec.with_atom_map(*map);
                    }
                    QueryAtom::new(AtomId::new(i), spec)
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
    pub(super) fn topology(count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
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
    fn product(maps: &[u32]) -> ProductBuilder {
        convert_template(
            &graph(&maps.iter().map(|&i| Some(i)).collect::<Vec<_>>()),
            0,
        )
        .unwrap()
    }
    fn run(
        p: &mut ProductBuilder,
        t: &TopologyBlock,
        q: &QueryGraph,
        pairs: &[(i32, i32)],
        c: &CoordinateBlock,
        conf: Option<(&mut Vec<[f64; 3]>, &mut bool)>,
    ) -> Result<(), ReactionProductError> {
        add_reactant_source(
            &Reaction::new(),
            p,
            ReactionInput {
                topology: t,
                coordinates: c,
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            pairs.iter().copied().map(Ok),
            q,
            7,
            conf,
            ReactionCoordinateSelection::Auto,
        )
    }
    #[test]
    fn reordered_sparse_signed_pairs_keep_query_keys_and_source_target_origins() {
        let q = graph(&[Some(11), Some(12), Some(13)]);
        let t = topology(4, &[]);
        let mut p = product(&[11, 13]);
        run(
            &mut p,
            &t,
            &q,
            &[(2, 3), (0, 1)],
            &CoordinateBlock::default(),
            None,
        )
        .unwrap();
        assert_eq!(
            p.atom_origins
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [1, 3]
        );
        for (a, id) in p.topology.atoms.iter().zip([1u32, 3]) {
            assert_eq!(a.prop("react_atom_idx"), Some(&PropertyValue::UInt(id)));
            assert_eq!(a.prop("react_idx"), Some(&PropertyValue::UInt(7)));
        }
    }
    #[test]
    fn traversal_adds_neighbors_once_before_later_matched_atom_and_keeps_source_bond_origins() {
        let q = graph(&[Some(11), Some(12)]);
        let t = topology(3, &[(0, 2), (2, 1)]);
        let mut p = product(&[11, 12]);
        run(
            &mut p,
            &t,
            &q,
            &[(0, 0), (1, 1)],
            &CoordinateBlock::default(),
            None,
        )
        .unwrap();
        assert_eq!(p.topology.atoms.len(), 3);
        assert_eq!(p.topology.bonds.len(), 2);
        assert_eq!(
            p.atom_origins
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(
            p.bond_origins
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [0, 1]
        );
    }
    #[test]
    fn skipped_matched_atoms_and_their_attached_components_are_not_added() {
        let q = graph(&[Some(11), Some(99)]);
        let t = topology(4, &[(0, 2), (2, 1), (1, 3)]);
        let mut p = product(&[11]);
        run(
            &mut p,
            &t,
            &q,
            &[(0, 0), (1, 1)],
            &CoordinateBlock::default(),
            None,
        )
        .unwrap();
        assert_eq!(
            p.atom_origins
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [0, 2]
        );
        assert_eq!(p.topology.bonds.len(), 1);
    }
    #[test]
    fn replicated_product_atoms_receive_properties_and_one_coordinate_per_copy_after_growth() {
        let q = graph(&[Some(11)]);
        let mut t = topology(2, &[(0, 1)]);
        t.atoms[0].set_isotope(Some(13));
        let mut p = product(&[11, 11]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[1., 2., 3.], [4., 5., 6.]], true)],
            ..Default::default()
        };
        let mut points = vec![[9.; 3]];
        let mut is3d = false;
        run(
            &mut p,
            &t,
            &q,
            &[(0, 0)],
            &c,
            Some((&mut points, &mut is3d)),
        )
        .unwrap();
        assert_eq!(p.topology.atoms.len(), 4);
        assert_eq!(
            points,
            [[1., 2., 3.], [1., 2., 3.], [4., 5., 6.], [4., 5., 6.]]
        );
        assert!(is3d);
        assert_eq!(p.topology.atoms[0].isotope(), None);
        assert_eq!(p.topology.atoms[1].isotope(), None);
        // ChemicalReaction's default implicit flag is false. Ordinary template
        // C keeps its isotope until the source caller explicitly enables it.
        let mut reaction = Reaction::new();
        reaction.implicit_properties = true;
        let mut implicit = product(&[11, 11]);
        add_reactant_source(
            &reaction,
            &mut implicit,
            ReactionInput {
                topology: &t,
                coordinates: &c,
                properties: &MoleculeProperties::default(),
                rings: None,
                valence: None,
            },
            [(0, 0)].into_iter().map(Ok),
            &q,
            7,
            None,
            ReactionCoordinateSelection::Auto,
        )
        .unwrap();
        assert_eq!(implicit.topology.atoms[0].isotope(), Some(13));
        assert_eq!(implicit.topology.atoms[1].isotope(), Some(13));
    }
    #[test]
    fn no_source_conformers_still_resizes_existing_product_buffer_only_at_final_phase() {
        let q = graph(&[Some(11)]);
        let t = topology(2, &[(0, 1)]);
        let mut p = product(&[11]);
        let mut points = vec![[9.; 3], [8.; 3], [7.; 3]];
        let mut is3d = true;
        run(
            &mut p,
            &t,
            &q,
            &[(0, 0)],
            &CoordinateBlock::default(),
            Some((&mut points, &mut is3d)),
        )
        .unwrap();
        assert_eq!(points, [[9.; 3], [8.; 3]]);
        assert!(is3d);
    }
    #[test]
    fn absent_output_conformer_does_not_read_malformed_source_coordinate_rows() {
        let q = graph(&[Some(11)]);
        let t = topology(1, &[]);
        let mut p = product(&[11]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![], true)],
            ..Default::default()
        };
        run(&mut p, &t, &q, &[(0, 0)], &c, None).unwrap();
        assert_eq!(p.atom_origins[0].unwrap().row.index(), 0);
    }
    #[test]
    fn mapping_failure_precedes_atom_transfer_and_final_coordinate_resize() {
        let q = graph(&[Some(11)]);
        let t = topology(1, &[]);
        let mut p = product(&[11]);
        let before = p.topology.clone();
        let mut points = vec![[9.; 3], [8.; 3]];
        let mut is3d = false;
        assert!(
            run(
                &mut p,
                &t,
                &q,
                &[(-1, 0)],
                &CoordinateBlock::default(),
                Some((&mut points, &mut is3d))
            )
            .is_err()
        );
        assert_eq!(p.topology, before);
        assert_eq!(points.len(), 2);
        assert!(p.atom_origins[0].is_none());
    }
    #[test]
    fn earlier_match_properties_and_provenance_survive_later_transfer_error_without_coordinate_resize()
     {
        let q = graph(&[Some(11), Some(12)]);
        let mut t = topology(2, &[]);
        t.atoms[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        let mut p = product(&[11, 12]);
        p.topology.atoms[1]
            .set_prop("molInversionFlag", PropertyValue::String("bad".into()))
            .unwrap();
        let mut points = vec![[9.; 3]];
        let mut is3d = false;
        assert!(
            run(
                &mut p,
                &t,
                &q,
                &[(0, 0), (1, 1)],
                &CoordinateBlock::default(),
                Some((&mut points, &mut is3d))
            )
            .is_err()
        );
        assert_eq!(p.atom_origins[0].unwrap().row.index(), 0);
        assert!(p.atom_origins[1].is_none());
        assert_eq!(
            p.topology.atoms[1].prop("react_atom_idx"),
            Some(&PropertyValue::UInt(1))
        );
        assert_eq!(points, [[9.; 3]]);
    }
    #[test]
    fn signed_bit_index_with_empty_bookmark_fails_when_loop_reaches_bit_not_in_mapping() {
        let q = graph(&[Some(11)]);
        let t = topology(0, &[]);
        let mut p = product(&[]);
        p.bookmarks = BTreeMap::from([(11, vec![])]);
        assert!(matches!(
            run(
                &mut p,
                &t,
                &q,
                &[(0, -1)],
                &CoordinateBlock::default(),
                None
            ),
            Err(ReactionProductError::Invariant {
                stage: "addReactantAtomsAndBonds",
                detail: "mapped bit index out of range",
                ..
            })
        ));
    }
    #[test]
    fn duplicate_match_pairs_are_not_deduplicated_or_rejected_and_share_reached_mapping_rows() {
        let q = graph(&[Some(11)]);
        let t = topology(2, &[(0, 1)]);
        let mut p = product(&[11]);
        run(
            &mut p,
            &t,
            &q,
            &[(0, 0), (0, 0)],
            &CoordinateBlock::default(),
            None,
        )
        .unwrap();
        assert_eq!(p.topology.atoms.len(), 3);
        assert_eq!(p.topology.bonds.len(), 2);
        // Native mapping push_back keeps both occurrences; BFS replicates the
        // neighbor for both mapping entries, even when both entries name atom0.
        assert_eq!(
            p.atom_origins
                .iter()
                .map(|x| x.unwrap().row.index())
                .collect::<Vec<_>>(),
            [0, 1, 1]
        );
    }
    #[test]
    fn enhanced_groups_copy_after_chirality_and_before_coordinate_failure() {
        let q = graph(&[Some(11)]);
        let mut t = topology(1, &[]);
        t.atoms[0].set_chiral_tag(ChiralTag::Other);
        t.stereo_groups = vec![
            StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(0)], vec![])
                .expect("valid distinct stereo members")
                .with_id(3),
        ];
        let mut p = product(&[11]);
        p.topology.atoms[0].set_chiral_tag(ChiralTag::Other);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![], true)],
            ..Default::default()
        };
        let mut points = vec![];
        let mut is3d = false;
        assert!(
            run(
                &mut p,
                &t,
                &q,
                &[(0, 0)],
                &c,
                Some((&mut points, &mut is3d))
            )
            .is_err()
        );
        assert_eq!(p.topology.stereo_groups.len(), 1);
        assert_eq!(p.topology.stereo_groups[0].atoms(), [AtomId::new(0)]);
        assert_eq!(points.len(), 1);
        assert!(is3d);
    }
}

#[cfg(test)]
mod complete_one_product_set_source_tests {
    use super::complete_add_reactant_source_tests::{graph, topology};
    use super::*;
    use cosmolkit_model::{Conformer2D, PropertyValue, StereoGroup, StereoGroupKind};
    fn reaction(reactants: Vec<QueryGraph>, products: Vec<QueryGraph>) -> Reaction {
        let mut r = Reaction::new();
        r.reactants = reactants;
        r.products = products;
        r
    }
    fn input<'a>(
        t: &'a TopologyBlock,
        c: &'a CoordinateBlock,
        props: &'a MoleculeProperties,
    ) -> ReactionInput<'a> {
        ReactionInput {
            topology: t,
            coordinates: c,
            properties: props,
            rings: None,
            valence: None,
        }
    }
    #[test]
    fn source_input_match_vector_precondition_precedes_zero_products() {
        let t = topology(0, &[]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        assert!(matches!(
            one_product_set(&Reaction::new(), &[input(&t, &c, &props)], &[], &[], 5),
            Err(ReactionRunError::ReactantArity {
                expected: 1,
                actual: 0
            })
        ));
    }
    #[test]
    fn selector_projection_has_its_own_arity_error_after_source_vector_precondition() {
        let t = topology(0, &[]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        assert!(matches!(
            one_product_set(
                &Reaction::new(),
                &[input(&t, &c, &props)],
                &[vec![]],
                &[],
                5
            ),
            Err(ReactionRunError::CoordinateSelectionArity {
                expected: 1,
                actual: 0
            })
        ));
    }
    #[test]
    fn no_products_never_reaches_missing_reaction_template_inputs() {
        let r = reaction(vec![graph(&[Some(11)])], vec![]);
        assert!(one_product_set(&r, &[], &[], &[], 1).unwrap().is_empty());
    }
    #[test]
    fn missing_source_vector_slot_reports_reached_product_context_without_panicking() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        assert!(matches!(
            one_product_set(&r, &[], &[], &[], 23),
            Err(ReactionRunError::Product {
                set: 23,
                template: 0,
                source: ReactionProductError::Invariant {
                    stage: "generateOneProductSet",
                    detail: "reactant input index out of range",
                    ..
                }
            })
        ));
    }
    #[test]
    fn template_conversion_error_precedes_reached_missing_input_slot() {
        let mut q = graph(&[Some(11)]);
        q.atoms_mut()[0]
            .set_prop("_QueryMass", PropertyValue::String("bad".into()))
            .unwrap();
        let r = reaction(vec![graph(&[Some(11)])], vec![q]);
        assert!(matches!(
            one_product_set(&r, &[], &[], &[], 3),
            Err(ReactionRunError::Product {
                set: 3,
                template: 0,
                source: ReactionProductError::PropertyUInt {
                    key: "_QueryMass",
                    ..
                }
            })
        ));
    }
    #[test]
    fn absent_source_coordinates_do_not_create_conformer_for_explicit_selection() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let t = topology(1, &[]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        let result = one_product_set(
            &r,
            &[input(&t, &c, &props)],
            &[vec![0]],
            &[ReactionCoordinateSelection::ThreeD { id: 99 }],
            0,
        )
        .unwrap();
        assert!(result[0].coordinates.conformers_3d.is_empty());
        assert!(result[0].coordinates.conformers_2d.is_empty());
        assert!(result[0].rings.is_none());
    }
    #[test]
    fn source_dimension_is_projected_into_one_fresh_id_zero_conformer() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let t = topology(1, &[]);
        let props = MoleculeProperties::default();
        for three in [false, true] {
            let c = if three {
                CoordinateBlock {
                    conformers_3d: vec![Conformer3D::new(17, vec![[1., 2., 3.]], true)],
                    ..Default::default()
                }
            } else {
                CoordinateBlock {
                    conformers_2d: vec![Conformer2D::new(19, vec![[1., 2.]])],
                    ..Default::default()
                }
            };
            let result = one_product_set(
                &r,
                &[input(&t, &c, &props)],
                &[vec![0]],
                &[ReactionCoordinateSelection::Auto],
                0,
            )
            .unwrap();
            assert_eq!(result[0].coordinates.conformers_3d.len(), 1);
            let conf = &result[0].coordinates.conformers_3d[0];
            assert_eq!(conf.id(), 0);
            assert_eq!(conf.is_3d(), three);
            assert_eq!(conf.coordinates(), [[1., 2., if three { 3. } else { 0. }]]);
        }
    }
    #[test]
    fn later_reactant_growth_preserves_earlier_coordinates_and_new_template_rows_stay_zero() {
        let r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)])],
            vec![graph(&[Some(11), Some(12), None])],
        );
        let t0 = topology(2, &[(0, 1)]);
        let t1 = topology(2, &[(0, 1)]);
        let c0 = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[1.; 3], [2.; 3]], true)],
            ..Default::default()
        };
        let c1 = CoordinateBlock {
            conformers_2d: vec![Conformer2D::new(0, vec![[3.; 2], [4.; 2]])],
            ..Default::default()
        };
        let props = MoleculeProperties::default();
        let result = one_product_set(
            &r,
            &[input(&t0, &c0, &props), input(&t1, &c1, &props)],
            &[vec![0], vec![0]],
            &[ReactionCoordinateSelection::Auto; 2],
            0,
        )
        .unwrap();
        assert_eq!(
            result[0].coordinates.conformers_3d[0].coordinates(),
            [[1.; 3], [3., 3., 0.], [0.; 3], [2.; 3], [4., 4., 0.]]
        );
        assert!(result[0].coordinates.conformers_3d[0].is_3d());
        assert!(result[0].atom_origins[2].is_none());
        assert_eq!(result[0].atom_origins[3].unwrap().input, 0);
        assert_eq!(result[0].atom_origins[4].unwrap().input, 1);
    }
    #[test]
    fn add_conformer_atom_count_precondition_runs_even_without_reaction_template_iterations() {
        let r = reaction(vec![], vec![graph(&[None])]);
        let t = topology(1, &[]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[1.; 3]], true)],
            ..Default::default()
        };
        let props = MoleculeProperties::default();
        assert!(matches!(
            one_product_set(
                &r,
                &[input(&t, &c, &props)],
                &[vec![]],
                &[ReactionCoordinateSelection::Auto],
                4
            ),
            Err(ReactionRunError::Product {
                set: 4,
                template: 0,
                source: ReactionProductError::Coordinate(_)
            })
        ));
    }
    #[test]
    fn only_directed_single_input_bonds_trigger_source_direction_restoration() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        for order in [BondOrder::Single, BondOrder::Double] {
            for direction in [
                BondDirection::None,
                BondDirection::BeginWedge,
                BondDirection::BeginDash,
                BondDirection::EndDownRight,
                BondDirection::EndUpRight,
                BondDirection::EitherDouble,
                BondDirection::Unknown,
            ] {
                let mut t = topology(2, &[(0, 1)]);
                t.bonds[0].set_order(order);
                t.bonds[0].set_direction(direction);
                let result = one_product_set(
                    &r,
                    &[input(&t, &c, &props)],
                    &[vec![0]],
                    &[ReactionCoordinateSelection::Auto],
                    0,
                )
                .unwrap();
                assert_eq!(
                    result[0].rings.is_some(),
                    order == BondOrder::Single && direction != BondDirection::None
                );
            }
        }
    }
    #[test]
    fn template_stereo_groups_are_copied_after_all_reactants_with_map_order_and_read_id() {
        let mut q = graph(&[Some(11), Some(12)]);
        q.add_stereo_group(
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(1), AtomId::new(0)],
                vec![],
            )
            .expect("valid distinct stereo members")
            .with_id(7),
        );
        let r = reaction(vec![graph(&[Some(11)]), graph(&[Some(12)])], vec![q]);
        let t = topology(1, &[]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        let result = one_product_set(
            &r,
            &[input(&t, &c, &props); 2],
            &[vec![0], vec![0]],
            &[ReactionCoordinateSelection::Auto; 2],
            0,
        )
        .unwrap();
        assert_eq!(
            result[0].topology.stereo_groups[0].atoms(),
            [AtomId::new(1), AtomId::new(0)]
        );
        assert_eq!(result[0].topology.stereo_groups[0].id(), Some(7));
    }
    #[test]
    fn final_property_cache_uses_source_nonstrict_valence_and_canonical_adjacency() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let t = topology(6, &[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]);
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        let result = one_product_set(
            &r,
            &[input(&t, &c, &props)],
            &[vec![0]],
            &[ReactionCoordinateSelection::Auto],
            0,
        )
        .unwrap();
        assert_eq!(result[0].valence.explicit_valence, [5, 1, 1, 1, 1, 1]);
        assert_eq!(
            result[0].topology.adjacency,
            AdjacencyList::from_topology(6, &result[0].topology.bonds)
        );
    }
    #[test]
    fn product_template_order_and_input_origins_are_independent_and_inputs_unchanged() {
        let r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)])],
            vec![graph(&[Some(12)]), graph(&[Some(11)])],
        );
        let t = topology(1, &[]);
        let before = t.clone();
        let c = CoordinateBlock::default();
        let mut props = MoleculeProperties::default();
        props.set_prop("not_copied", 17i32).unwrap();
        let result = one_product_set(
            &r,
            &[input(&t, &c, &props); 2],
            &[vec![0], vec![0]],
            &[ReactionCoordinateSelection::Auto; 2],
            0,
        )
        .unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(result[0].atom_origins[0].unwrap().input, 1);
        assert_eq!(result[1].atom_origins[0].unwrap().input, 0);
        assert!(
            result
                .iter()
                .all(|p| p.properties.prop("not_copied").is_none())
        );
        assert_eq!(t, before);
    }
    #[test]
    fn later_product_failure_retains_exact_set_and_template_context_and_input_state() {
        let mut bad = graph(&[Some(11)]);
        bad.atoms_mut()[0]
            .set_prop("_QueryMass", PropertyValue::String("bad".into()))
            .unwrap();
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)]), bad]);
        let t = topology(1, &[]);
        let before = t.clone();
        let c = CoordinateBlock::default();
        let props = MoleculeProperties::default();
        assert!(matches!(
            one_product_set(
                &r,
                &[input(&t, &c, &props)],
                &[vec![0]],
                &[ReactionCoordinateSelection::Auto],
                9
            ),
            Err(ReactionRunError::Product {
                set: 9,
                template: 1,
                source: ReactionProductError::PropertyUInt { .. }
            })
        ));
        assert_eq!(t, before);
    }
}

#[cfg(test)]
mod complete_run_reactants_source_tests {
    use super::complete_add_reactant_source_tests::{graph, topology};
    use super::*;
    use cosmolkit_model::PropertyValue;
    use cosmolkit_types::Element;
    #[test]
    fn copied_user_atom_properties_use_origins_and_preserve_template_priority() {
        let mut product = graph(&[Some(11), Some(11), Some(12), None]);
        product
            .atom_mut(0)
            .unwrap()
            .set_prop("tracking_id", PropertyValue::Int(99))
            .unwrap();
        let r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12), Some(13)])],
            vec![product],
        );
        let mut first = topology(1, &[]);
        let mut second = topology(2, &[]);
        first.atoms[0].set_prop("tracking_id", 42).unwrap();
        first.atoms[0].set_prop("note", "first").unwrap();
        first.atoms[0].set_prop("_CIPCode", "R").unwrap();
        first.atoms[0].set_computed_prop("derived", 7).unwrap();
        second.atoms[0].set_prop("tracking_id", 84).unwrap();
        second.atoms[1].set_prop("deleted", true).unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let inputs = [
            input(&first, &coordinates, &properties),
            input(&second, &coordinates, &properties),
        ];
        let native = run_reactants_source(&r, &inputs, &ReactionRunParams::default()).unwrap();
        let copied = run_reactants_source(
            &r,
            &inputs,
            &ReactionRunParams {
                copy_atom_properties: true,
                ..Default::default()
            },
        )
        .unwrap();
        let atoms = &copied[0][0].topology.atoms;
        assert_eq!(atoms.len(), 4);
        assert_eq!(atoms[0].prop("tracking_id"), Some(&PropertyValue::Int(99)));
        assert_eq!(atoms[1].prop("tracking_id"), Some(&PropertyValue::Int(42)));
        assert_eq!(atoms[2].prop("tracking_id"), Some(&PropertyValue::Int(84)));
        assert!(atoms[3].prop("tracking_id").is_none());
        for atom in atoms {
            assert!(atom.prop("deleted").is_none());
            assert!(atom.prop("_CIPCode").is_none());
            assert!(atom.prop("derived").is_none());
        }
        assert!(native[0][0].topology.atoms[1].prop("tracking_id").is_none());
        assert_eq!(
            native[0][0].topology.atoms[0].prop("tracking_id"),
            Some(&PropertyValue::Int(99))
        );
        assert_eq!(native[0][0].atom_origins, copied[0][0].atom_origins);
        assert_eq!(native[0][0].valence, copied[0][0].valence);
        assert_eq!(
            first.atoms[0].prop("tracking_id"),
            Some(&PropertyValue::Int(42))
        );
    }
    fn reaction(reactants: Vec<QueryGraph>, products: Vec<QueryGraph>) -> Reaction {
        let mut r = Reaction::new();
        r.reactants = reactants;
        r.products = products;
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
    #[test]
    fn source_initialization_precondition_precedes_arity_and_zero_products() {
        let r = Reaction::new();
        let t = topology(0, &[]);
        assert!(matches!(
            run_reactants_source(
                &r,
                &[input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                )],
                &ReactionRunParams::default()
            ),
            Err(ReactionRunError::NeedsInitialization)
        ));
    }
    #[test]
    fn source_reactant_arity_precedes_zero_products_and_selector_projection_validation() {
        let r = reaction(vec![graph(&[Some(11)])], vec![]);
        let params = ReactionRunParams {
            coordinate_selections: vec![ReactionCoordinateSelection::Auto; 2],
            ..Default::default()
        };
        assert!(matches!(
            run_reactants_source(&r, &[], &params),
            Err(ReactionRunError::ReactantArity {
                expected: 1,
                actual: 0
            })
        ));
    }
    #[test]
    fn no_product_templates_return_before_search_and_coordinate_reads() {
        let r = reaction(vec![graph(&[Some(11)])], vec![]);
        let mut t = topology(1, &[]);
        t.atoms[0] = t.atoms[0].clone().with_id(AtomId::new(99));
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![], true)],
            ..Default::default()
        };
        let p = MoleculeProperties::default();
        assert!(
            run_reactants_source(&r, &[input(&t, &c, &p)], &ReactionRunParams::default())
                .unwrap()
                .is_empty()
        );
    }
    #[test]
    fn selector_projection_arity_reports_its_exact_count_after_source_arity() {
        let r = reaction(vec![graph(&[Some(11)])], vec![]);
        let t = topology(1, &[]);
        let params = ReactionRunParams {
            coordinate_selections: vec![ReactionCoordinateSelection::Auto; 2],
            ..Default::default()
        };
        assert!(matches!(
            run_reactants_source(
                &r,
                &[input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                )],
                &params
            ),
            Err(ReactionRunError::CoordinateSelectionArity {
                expected: 1,
                actual: 2
            })
        ));
    }
    #[test]
    fn first_unmatched_reactant_prevents_later_search_and_coordinate_selection() {
        let r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)])],
            vec![graph(&[Some(11), Some(12)])],
        );
        let mut t0 = topology(1, &[]);
        t0.atoms[0].set_element(Element::O);
        let mut t1 = topology(1, &[]);
        t1.atoms[0] = t1.atoms[0].clone().with_id(AtomId::new(99));
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let params = ReactionRunParams {
            coordinate_selections: vec![ReactionCoordinateSelection::ThreeD { id: 99 }; 2],
            ..Default::default()
        };
        assert!(
            run_reactants_source(&r, &[input(&t0, &c, &p), input(&t1, &c, &p)], &params)
                .unwrap()
                .is_empty()
        );
    }
    #[test]
    fn physical_match_order_and_cartesian_product_order_survive_all_finite_and_unlimited_limits() {
        let mut r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)])],
            vec![graph(&[Some(11), Some(12)])],
        );
        r.match_params.uniquify = true;
        r.match_params.max_matches = 1;
        let t0 = topology(2, &[]);
        let t1 = topology(3, &[]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let inputs = [input(&t0, &c, &p), input(&t1, &c, &p)];
        for max in [0, 1, 2, 5, 6, 99] {
            let result = run_reactants_source(
                &r,
                &inputs,
                &ReactionRunParams {
                    max_products: max,
                    ..Default::default()
                },
            )
            .unwrap();
            let expected = if max == 0 { 6 } else { (max as usize).min(6) };
            assert_eq!(result.len(), expected);
            let origins = result
                .iter()
                .map(|set| {
                    assert_eq!(set.len(), 1);
                    let rows = &set[0].atom_origins;
                    assert_eq!(rows[0].unwrap().input, 0);
                    assert_eq!(rows[1].unwrap().input, 1);
                    (rows[0].unwrap().row.index(), rows[1].unwrap().row.index())
                })
                .collect::<Vec<_>>();
            assert_eq!(
                origins,
                [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)][..expected]
            );
        }
        assert!(r.match_params.uniquify);
        assert_eq!(r.match_params.max_matches, 1);
    }
    #[test]
    fn protected_presence_filters_after_source_matching_limit_without_refilling() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        for protected in [
            PropertyValue::Bool(false),
            PropertyValue::Int(0),
            PropertyValue::String("raw".into()),
        ] {
            let mut t = topology(3, &[]);
            t.atoms[0].set_prop("_protected", protected).unwrap();
            for max in [0, 2] {
                let result = run_reactants_source(
                    &r,
                    &[input(&t, &c, &p)],
                    &ReactionRunParams {
                        max_products: max,
                        ..Default::default()
                    },
                )
                .unwrap();
                let rows = result
                    .iter()
                    .map(|s| s[0].atom_origins[0].unwrap().row.index())
                    .collect::<Vec<_>>();
                assert_eq!(rows, if max == 0 { vec![1, 2] } else { vec![1] });
            }
            assert!(t.atoms[0].prop("_protected").is_some());
        }
    }
    #[test]
    fn reached_coordinate_selection_error_keeps_product_set_and_template_context() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let t = topology(1, &[]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![[1.; 3]], true)],
            ..Default::default()
        };
        let p = MoleculeProperties::default();
        let params = ReactionRunParams {
            coordinate_selections: vec![ReactionCoordinateSelection::ThreeD { id: 99 }],
            ..Default::default()
        };
        assert!(matches!(
            run_reactants_source(&r, &[input(&t, &c, &p)], &params),
            Err(ReactionRunError::Product {
                set: 0,
                template: 0,
                source: ReactionProductError::CoordinateSelection(_)
            })
        ));
    }
    #[test]
    fn approved_d1_public_projection_initializes_only_private_temporary_and_preserves_inputs() {
        let mut r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        r.needs_init = true;
        let before_r = r.reactants.clone();
        let before_p = r.products.clone();
        let t = topology(2, &[]);
        let before = t.clone();
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let result = run_reactants(
            &mut r,
            &[input(&t, &c, &p)],
            &ReactionRunParams {
                max_products: 0,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(result.len(), 2);
        assert!(r.is_initialized());
        assert_eq!(r.reactants, before_r);
        assert_eq!(r.products, before_p);
        assert_eq!(t, before);
    }
    #[test]
    fn reached_empty_reactant_combination_error_is_propagated_without_fabricating_product() {
        let r = reaction(vec![], vec![graph(&[None])]);
        assert!(matches!(
            run_reactants_source(&r, &[], &ReactionRunParams::default()),
            Err(ReactionRunError::CombinationLevel { level: 0, count: 0 })
        ));
    }
}

#[cfg(test)]
mod complete_run_reactant_source_tests {
    use super::complete_add_reactant_source_tests::{graph, topology};
    use super::*;
    use cosmolkit_model::PropertyValue;
    use cosmolkit_types::Element;
    fn reaction(reactants: Vec<QueryGraph>, products: Vec<QueryGraph>) -> Reaction {
        let mut r = Reaction::new();
        r.reactants = reactants;
        r.products = products;
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
    #[test]
    fn initialization_precondition_precedes_zero_products_and_selected_index() {
        let t = topology(0, &[]);
        assert!(matches!(
            run_reactant_source(
                &Reaction::new(),
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                ),
                usize::MAX,
                &ReactionSingleRunParams::default()
            ),
            Err(ReactionRunError::NeedsInitialization)
        ));
    }
    #[test]
    fn zero_products_return_before_index_search_or_coordinate_reads() {
        let r = reaction(vec![], vec![]);
        let mut t = topology(1, &[]);
        t.atoms[0] = t.atoms[0].clone().with_id(AtomId::new(99));
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(0, vec![], true)],
            ..Default::default()
        };
        assert!(
            run_reactant_source(
                &r,
                input(&t, &c, &MoleculeProperties::default()),
                usize::MAX,
                &ReactionSingleRunParams::default()
            )
            .unwrap()
            .is_empty()
        );
    }
    #[test]
    fn source_index_precondition_precedes_search_of_invalid_input() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let mut t = topology(1, &[]);
        t.atoms[0] = t.atoms[0].clone().with_id(AtomId::new(99));
        assert!(matches!(
            run_reactant_source(
                &r,
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                ),
                1,
                &ReactionSingleRunParams::default()
            ),
            Err(ReactionRunError::ReactantTemplateIndex { index: 1, count: 1 })
        ));
    }
    #[test]
    fn only_selected_template_is_matched_and_unselected_slots_materialize_as_empty_reagents() {
        let mut unused = graph(&[Some(11)]);
        unused.atoms_mut()[0]
            .set_prop(
                "molAtomMapNumber",
                PropertyValue::String("unreached".into()),
            )
            .unwrap();
        let r = reaction(
            vec![unused, graph(&[Some(12)])],
            vec![graph(&[Some(11), Some(12)])],
        );
        let t = topology(1, &[]);
        let result = run_reactant_source(
            &r,
            input(
                &t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            1,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(result.len(), 1);
        assert!(result[0][0].atom_origins[0].is_none());
        assert_eq!(result[0][0].atom_origins[1].unwrap().input, 1);
        assert_eq!(
            result[0][0].topology.atoms[1].prop("react_idx"),
            Some(&PropertyValue::UInt(1))
        );
    }
    #[test]
    fn source_fixed_thousand_match_limit_overrides_stored_match_options_without_reordering() {
        let mut r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        r.match_params.max_matches = 1;
        r.match_params.uniquify = true;
        let t = topology(1001, &[]);
        let result = run_reactant_source(
            &r,
            input(
                &t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            0,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(result.len(), 1000);
        for (row, set) in result.iter().enumerate() {
            assert_eq!(set.len(), 1);
            assert_eq!(set[0].atom_origins[0].unwrap().row.index(), row);
        }
        assert_eq!(r.match_params.max_matches, 1);
        assert!(r.match_params.uniquify);
    }
    #[test]
    fn detached_origin_projection_changes_atom_and_bond_input_slots_but_preserves_source_properties()
     {
        let mut r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)])],
            vec![graph(&[Some(12)])],
        );
        let t = topology(2, &[(0, 1)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let source = run_reactant_source(
            &r,
            input(&t, &c, &p),
            1,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        let projected = run_reactant(
            &mut r,
            input(&t, &c, &p),
            1,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(source.len(), 2);
        for (raw, public) in source.iter().zip(&projected) {
            assert_eq!(raw[0].topology, public[0].topology);
            assert!(raw[0].atom_origins.iter().flatten().all(|o| o.input == 1));
            assert!(raw[0].bond_origins.iter().flatten().all(|o| o.input == 1));
            assert!(
                public[0]
                    .atom_origins
                    .iter()
                    .flatten()
                    .all(|o| o.input == 0)
            );
            assert!(
                public[0]
                    .bond_origins
                    .iter()
                    .flatten()
                    .all(|o| o.input == 0)
            );
            assert!(
                public[0]
                    .topology
                    .atoms
                    .iter()
                    .all(|a| a.prop("react_idx") == Some(&PropertyValue::UInt(1)))
            );
        }
    }
    #[test]
    fn empty_slots_before_and_after_selected_reagent_preserve_source_coordinates_and_zero_new_rows()
    {
        let r = reaction(
            vec![graph(&[Some(11)]), graph(&[Some(12)]), graph(&[Some(13)])],
            vec![graph(&[Some(11), Some(12), Some(13)])],
        );
        let t = topology(2, &[(0, 1)]);
        let c = CoordinateBlock {
            conformers_3d: vec![Conformer3D::new(9, vec![[1.; 3], [2.; 3]], true)],
            ..Default::default()
        };
        let p = MoleculeProperties::default();
        let result = run_reactant_source(
            &r,
            input(&t, &c, &p),
            1,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(
            result[0][0].coordinates.conformers_3d[0].coordinates(),
            [[0.; 3], [1.; 3], [0.; 3], [2.; 3]]
        );
        assert!(result[0][0].coordinates.conformers_3d[0].is_3d());
        assert_eq!(result[0][0].coordinates.conformers_3d[0].id(), 0);
    }
    #[test]
    fn selected_no_match_returns_before_template_conversion_and_coordinate_selection() {
        let mut bad = graph(&[Some(11)]);
        bad.atoms_mut()[0]
            .set_prop("_QueryMass", PropertyValue::String("unreached".into()))
            .unwrap();
        let r = reaction(vec![graph(&[Some(11)])], vec![bad]);
        let mut t = topology(1, &[]);
        t.atoms[0].set_element(Element::O);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        assert!(
            run_reactant_source(
                &r,
                input(&t, &c, &p),
                0,
                &ReactionSingleRunParams {
                    coordinate_selection: ReactionCoordinateSelection::ThreeD { id: 99 }
                }
            )
            .unwrap()
            .is_empty()
        );
    }
    #[test]
    fn symmetric_matching_products_are_kept_in_source_order_without_deduplication() {
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        let t = topology(2, &[(0, 1)]);
        let c = CoordinateBlock::default();
        let p = MoleculeProperties::default();
        let result = run_reactant_source(
            &r,
            input(&t, &c, &p),
            0,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(result.len(), 2);
        assert_eq!(
            result
                .iter()
                .map(|s| s[0].atom_origins[0].unwrap().row.index())
                .collect::<Vec<_>>(),
            [0, 1]
        );
        assert!(
            result
                .iter()
                .all(|s| s[0].topology.atoms.len() == 2 && s[0].topology.bonds.len() == 1)
        );
    }
    #[test]
    fn later_product_failure_preserves_source_product_context_and_caller_input() {
        let mut bad = graph(&[Some(11)]);
        bad.atoms_mut()[0]
            .set_prop("_QueryMass", PropertyValue::String("bad".into()))
            .unwrap();
        let r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)]), bad]);
        let t = topology(2, &[]);
        let before = t.clone();
        assert!(matches!(
            run_reactant_source(
                &r,
                input(
                    &t,
                    &CoordinateBlock::default(),
                    &MoleculeProperties::default()
                ),
                0,
                &ReactionSingleRunParams::default()
            ),
            Err(ReactionRunError::Product {
                set: 0,
                template: 1,
                source: ReactionProductError::PropertyUInt { .. }
            })
        ));
        assert_eq!(t, before);
    }
    #[test]
    fn approved_d1_initialization_does_not_change_original_reaction_or_source_topology() {
        let mut r = reaction(vec![graph(&[Some(11)])], vec![graph(&[Some(11)])]);
        r.needs_init = true;
        let reactants = r.reactants.clone();
        let products = r.products.clone();
        let t = topology(1, &[]);
        let before = t.clone();
        let result = run_reactant(
            &mut r,
            input(
                &t,
                &CoordinateBlock::default(),
                &MoleculeProperties::default(),
            ),
            0,
            &ReactionSingleRunParams::default(),
        )
        .unwrap();
        assert_eq!(result.len(), 1);
        assert!(r.is_initialized());
        assert_eq!(r.reactants, reactants);
        assert_eq!(r.products, products);
        assert_eq!(t, before);
    }
}
