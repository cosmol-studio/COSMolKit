use crate::{Reaction, ReactionTemplateRemovalParams};
use cosmolkit_model::QueryGraph;

#[derive(Debug, Clone)]
pub struct ReactionTemplateRemoval {
    pub reaction: Reaction,
    pub removed_templates: Vec<QueryGraph>,
}

fn count_atoms_with_property(graph: &QueryGraph, property: &[u8]) -> u32 {
    // RDKit❗🔝: unsigned getNumAtomsWithDistinctProperty(const ROMol &mol,
    // RDKit❗🔝:                                          const std::string_view &prop) {
    // RDKit❗🔝:   unsigned numPropAtoms = 0;
    // RDKit❗🔝:   for (const auto atom : mol.atoms()) {
    // RDKit❗🔝:     if (atom->hasProp(prop)) {
    // RDKit❗🔝:       ++numPropAtoms;
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return numPropAtoms;
    // RDKit❗🔝: }
    let mut count = 0u32;
    for atom in graph.atoms() {
        if atom_has_property_source(atom, property) {
            count = count.wrapping_add(1);
        }
    }
    count
}

pub(crate) fn atom_has_property_source(
    atom: &impl crate::materialize::ReactionAtomPropertyRead,
    property: &[u8],
) -> bool {
    // RDKit❗🔝: bool hasProp(const std::string_view key) const { return d_props.hasVal(key); }
    // RDKit❗🔝:   bool hasVal(const std::string_view what) const {
    // RDKit❗🔝:     for (const auto &data : _data) {
    // RDKit❗🔝:       if (data.key == what) {
    // RDKit❗🔝:         return true;
    // RDKit❗🔝:       }
    // RDKit❗🔝:     }
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // RDKit❗🔝: inline constexpr std::string_view molAtomMapNumber = "molAtomMapNumber";
    // RDKit❗🔝: inline constexpr std::string_view _chiralPermutation = "_chiralPermutation";
    // RDKit❗🔝: inline constexpr std::string_view molParity = "molParity";
    // RDKit❗🔝: inline constexpr std::string_view molInversionFlag = "molInversionFlag";
    // Behavior: count each atom with a named property, not unique values. No
    // numeric/string conversion occurs, including for zero or wrong-tag raw
    // values. Byte keys preserve native string_view's non-UTF8/NUL semantics.
    // Explicit Option fields project known source dictionary presence; they
    // count once together with a raw key. Other semantic bool/list fields do
    // not carry independent dictionary presence and are never guessed here.
    // That independent MODEL transport capability remains a source-review gap;
    // this reaction caller requests only the represented atom-map property.
    // Complexity: one atom pass without allocation, copied graph, value set or
    // deduplication. BTreeMap O(log P) lookup replaces native Dict O(P) scan and
    // preserves exact key-presence membership; four fixed Option projections
    // add constant work and no owned-key allocation.
    let projected_presence = match property {
        b"molAtomMapNumber" => atom.typed_map().is_some(),
        b"_chiralPermutation" => atom.typed_permutation().is_some(),
        b"molParity" => atom.typed_parity().is_some(),
        b"molInversionFlag" => atom.typed_inversion().is_some(),
        _ => false,
    };
    atom.raw_property_exists(property) || projected_presence
}

fn is_agent(graph: &QueryGraph, threshold: f64) -> bool {
    // RDKit❗✔️: bool isReactionTemplateMoleculeAgent(const ROMol &mol, double agentThreshold) {
    // RDKit❗✔️:   unsigned numMappedAtoms = MolOps::getNumAtomsWithDistinctProperty(
    // RDKit❗✔️:       mol, common_properties::molAtomMapNumber);
    // RDKit❗✔️:   unsigned numAtoms = mol.getNumHeavyAtoms();
    // RDKit❗✔️:   return !(numAtoms > 0u && static_cast<double>(numMappedAtoms) /
    // RDKit❗✔️:                                     static_cast<double>(numAtoms) >=
    // RDKit❗✔️:                                 agentThreshold);
    // RDKit❗✔️: }
    // RDKit❗✔️: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit❗✔️:   unsigned int res = 0;
    // RDKit❗✔️:   for (const auto atom : atoms()) {
    // RDKit❗✔️:     if (atom->getAtomicNum() > 1) {
    // RDKit❗✔️:       ++res;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: };
    // Source counts map-property presence on every atom, but divides by
    // heavy atoms only. The heavy-count zero guard precedes the division;
    // preserve the exact negated >= comparison, including NaN/infinities.
    // Reuse canonical property-presence owner; two allocation-free atom passes.
    let mapped = count_atoms_with_property(graph, b"molAtomMapNumber");
    let heavy = graph
        .atoms()
        .iter()
        .filter(|atom| atom.atomic_number() > 1)
        .fold(0u32, |count, _| count.wrapping_add(1));
    !(heavy > 0 && (mapped as f64) / (heavy as f64) >= threshold)
}

fn remove_unmapped_reactants_source(
    reaction: &mut Reaction,
    threshold: f64,
    move_to_agents: bool,
    mut target: Option<&mut Vec<QueryGraph>>,
) {
    // RDKit❗❌: void ChemicalReaction::removeUnmappedReactantTemplates(
    // RDKit❗❌:     double thresholdUnmappedAtoms, bool moveToAgentTemplates,
    // RDKit❗❌:     MOL_SPTR_VECT *targetVector) {
    // RDKit❗❌:   MOL_SPTR_VECT res_reactantTemplates;
    // RDKit❗❌:   for (auto iter = beginReactantTemplates(); iter != endReactantTemplates();
    // RDKit❗❌:        ++iter) {
    // RDKit❗❌:     if (isReactionTemplateMoleculeAgent(*iter->get(), thresholdUnmappedAtoms)) {
    // RDKit❗❌:       if (moveToAgentTemplates) {
    // RDKit❗❌:         m_agentTemplates.push_back(*iter);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (targetVector) {
    // RDKit❗❌:         targetVector->push_back(*iter);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res_reactantTemplates.push_back(*iter);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   m_reactantTemplates.clear();
    // RDKit❗❌:   m_reactantTemplates.insert(m_reactantTemplates.begin(),
    // RDKit❗❌:                              res_reactantTemplates.begin(),
    // RDKit❗❌:                              res_reactantTemplates.end());
    // RDKit❗❌:   res_reactantTemplates.clear();
    // RDKit❗❌: }
    // Native's temporary vector contains references to retained templates.
    // The owned detached representation stores their indexes, then stably
    // installs those same graph values in the original vector allocation.
    // This preserves encounter order, original capacity and kept graph data,
    // without deep-copying every kept template or taking the role buffer.
    // All per-template agent/target appends happen before changing this role,
    // and agents are appended before the optional target, as in source.
    // Cost gap: duplicate outputs clone owned QueryGraph values, whereas native
    // copies shared pointers. O(total duplicated graph data) instead of O(M)
    // pointer copies; identity/alias transport remains a final comparison.
    let mut retained = Vec::new();
    for (index, template) in reaction.reactants.iter().enumerate() {
        if is_agent(template, threshold) {
            if move_to_agents {
                reaction.agents.push(template.clone());
            }
            if let Some(output) = target.as_deref_mut() {
                output.push(template.clone());
            }
        } else {
            retained.push(index);
        }
    }
    let mut positions = retained.iter().copied();
    let mut next = positions.next();
    let mut index = 0;
    reaction.reactants.retain(|_| {
        let keep = next == Some(index);
        if keep {
            next = positions.next();
        }
        index += 1;
        keep
    });
    retained.clear();
}

#[doc(hidden)]
pub fn without_unmapped_reactants(
    reaction: &Reaction,
    params: &ReactionTemplateRemovalParams,
) -> ReactionTemplateRemoval {
    // Immutable value projection owns the full source-shaped reaction copy.
    let mut result = reaction.clone();
    let mut removed_templates = Vec::new();
    remove_unmapped_reactants_source(
        &mut result,
        params.threshold_unmapped_atoms,
        params.move_to_agent_templates,
        Some(&mut removed_templates),
    );
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates,
    }
}

fn remove_unmapped_products_source(
    reaction: &mut Reaction,
    threshold: f64,
    move_to_agents: bool,
    mut target: Option<&mut Vec<QueryGraph>>,
) {
    // RDKit❗❌: void ChemicalReaction::removeUnmappedProductTemplates(
    // RDKit❗❌:     double thresholdUnmappedAtoms, bool moveToAgentTemplates,
    // RDKit❗❌:     MOL_SPTR_VECT *targetVector) {
    // RDKit❗❌:   MOL_SPTR_VECT res_productTemplates;
    // RDKit❗❌:   for (auto iter = beginProductTemplates(); iter != endProductTemplates();
    // RDKit❗❌:        ++iter) {
    // RDKit❗❌:     if (isReactionTemplateMoleculeAgent(*iter->get(), thresholdUnmappedAtoms)) {
    // RDKit❗❌:       if (moveToAgentTemplates) {
    // RDKit❗❌:         m_agentTemplates.push_back(*iter);
    // RDKit❗❌:       }
    // RDKit❗❌:       if (targetVector) {
    // RDKit❗❌:         targetVector->push_back(*iter);
    // RDKit❗❌:       }
    // RDKit❗❌:     } else {
    // RDKit❗❌:       res_productTemplates.push_back(*iter);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   m_productTemplates.clear();
    // RDKit❗❌:   m_productTemplates.insert(m_productTemplates.begin(),
    // RDKit❗❌:                             res_productTemplates.begin(),
    // RDKit❗❌:                             res_productTemplates.end());
    // RDKit❗❌:   res_productTemplates.clear();
    // RDKit❗❌: }
    // Native's temporary vector contains references to retained templates.
    // The owned detached representation stores their indexes, then stably
    // installs those same graph values in the original vector allocation.
    // This preserves encounter order, original capacity and kept graph data,
    // without deep-copying every kept template or taking the role buffer.
    // All per-template agent/target appends happen before changing this role,
    // and agents are appended before the optional target, as in source.
    // Cost gap: duplicate outputs clone owned QueryGraph values, whereas native
    // copies shared pointers. O(total duplicated graph data) instead of O(M)
    // pointer copies; identity/alias transport remains a final comparison.
    let mut retained = Vec::new();
    for (index, template) in reaction.products.iter().enumerate() {
        if is_agent(template, threshold) {
            if move_to_agents {
                reaction.agents.push(template.clone());
            }
            if let Some(output) = target.as_deref_mut() {
                output.push(template.clone());
            }
        } else {
            retained.push(index);
        }
    }
    let mut positions = retained.iter().copied();
    let mut next = positions.next();
    let mut index = 0;
    reaction.products.retain(|_| {
        let keep = next == Some(index);
        if keep {
            next = positions.next();
        }
        index += 1;
        keep
    });
    retained.clear();
}

#[doc(hidden)]
pub fn without_unmapped_products(
    reaction: &Reaction,
    params: &ReactionTemplateRemovalParams,
) -> ReactionTemplateRemoval {
    // Immutable value projection owns the full source-shaped reaction copy.
    let mut result = reaction.clone();
    let mut removed_templates = Vec::new();
    remove_unmapped_products_source(
        &mut result,
        params.threshold_unmapped_atoms,
        params.move_to_agent_templates,
        Some(&mut removed_templates),
    );
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates,
    }
}

fn remove_agents_source(reaction: &mut Reaction, target: Option<&mut Vec<QueryGraph>>) {
    // RDKit❗❌: void ChemicalReaction::removeAgentTemplates(MOL_SPTR_VECT *targetVector) {
    // RDKit❗❌:   if (targetVector) {
    // RDKit❗❌:     for (auto iter = beginAgentTemplates(); iter != endAgentTemplates();
    // RDKit❗❌:          ++iter) {
    // RDKit❗❌:       targetVector->push_back(*iter);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   m_agentTemplates.clear();
    // RDKit❗❌: }
    // Append to the optional existing target in encounter order before the
    // actual source vector clear; retain its allocation and all other state.
    // Owned detached copies cost graph-data size where native copies shared
    // pointers. That explicit cost/alias transport gap is not a chemistry rule.
    if let Some(output) = target {
        for template in &reaction.agents {
            output.push(template.clone());
        }
    }
    reaction.agents.clear();
}

#[doc(hidden)]
pub fn without_agents(reaction: &Reaction) -> ReactionTemplateRemoval {
    let mut result = reaction.clone();
    let mut removed_templates = Vec::new();
    remove_agents_source(&mut result, Some(&mut removed_templates));
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates,
    }
}

#[cfg(test)]
mod source574_property_presence_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;

    fn graph(specs: Vec<AtomSpec>) -> QueryGraph {
        QueryGraph::from_parts(
            specs
                .into_iter()
                .enumerate()
                .map(|(i, spec)| QueryAtom::new(AtomId::new(i), spec))
                .collect(),
            vec![],
            vec![],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }

    #[test]
    fn source574_counts_atoms_even_when_values_repeat_or_have_different_tags() {
        let mut atoms = vec![];
        for (i, value) in [
            PropertyValue::Int(7),
            PropertyValue::Int(7),
            PropertyValue::Bool(false),
        ]
        .into_iter()
        .enumerate()
        {
            let mut atom = QueryAtom::new(AtomId::new(i), AtomSpec::new(Element::C));
            atom.set_prop("ordinary", value).unwrap();
            atoms.push(atom);
        }
        let graph = QueryGraph::from_parts(atoms, vec![], vec![], vec![], vec![], vec![]).unwrap();
        assert_eq!(count_atoms_with_property(&graph, b"ordinary"), 3);
        assert_eq!(count_atoms_with_property(&graph, b"absent"), 0);
    }

    #[test]
    fn source574_zero_map_and_duplicate_raw_typed_presence_count_once_per_atom() {
        let mut graph = graph(vec![
            AtomSpec::new(Element::C).with_atom_map(0),
            AtomSpec::new(Element::C).with_atom_map(4),
            AtomSpec::new(Element::C).with_atom_map(4),
        ]);
        graph
            .atom_mut(0)
            .unwrap()
            .set_prop("molAtomMapNumber", PropertyValue::Bool(false))
            .unwrap();
        assert_eq!(count_atoms_with_property(&graph, b"molAtomMapNumber"), 3);
        assert!(!is_agent(&graph, 1.0));
        assert!(is_agent(&graph, f64::NAN));
    }

    #[test]
    fn source574_byte_property_names_and_explicit_option_presence() {
        let mut atom = QueryAtom::new(
            AtomId::new(0),
            AtomSpec::new(Element::C)
                .with_chiral_permutation(0)
                .with_mol_parity(0)
                .with_mol_inversion_flag(0),
        );
        atom.set_prop(vec![0, 255], PropertyValue::UInt(0)).unwrap();
        let graph =
            QueryGraph::from_parts(vec![atom], vec![], vec![], vec![], vec![], vec![]).unwrap();
        for key in [
            b"_chiralPermutation".as_slice(),
            b"molParity",
            b"molInversionFlag",
            &[0, 255],
        ] {
            assert_eq!(count_atoms_with_property(&graph, key), 1);
        }
        assert_eq!(count_atoms_with_property(&graph, b"molAtomMapNumber"), 0);
    }
}

#[cfg(test)]
mod complete_is_agent_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    fn graph(elements: &[Element]) -> QueryGraph {
        QueryGraph::from_parts(
            elements
                .iter()
                .enumerate()
                .map(|(i, e)| {
                    QueryAtom::new(AtomId::new(i), AtomSpec::new(*e).with_atom_map(i as u32))
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
    #[test]
    fn heavy_denominator_mapped_hydrogens_zero_guard_and_nonfinite_thresholds_follow_source() {
        for elements in [&[][..], &[Element::H, Element::H][..]] {
            let molecule = graph(elements);
            for threshold in [f64::NEG_INFINITY, 0.0, 1.0, f64::INFINITY, f64::NAN] {
                assert!(is_agent(&molecule, threshold));
            }
        }
        let mut molecule = graph(&[Element::C, Element::H, Element::H]);
        for threshold in [f64::NEG_INFINITY, 0.0, 1.0, 3.0] {
            assert!(!is_agent(&molecule, threshold));
        }
        for threshold in [4.0, f64::INFINITY, f64::NAN] {
            assert!(is_agent(&molecule, threshold));
        }
        // Atom-index validation is not a source precondition of this count.
        molecule.atoms_mut()[0] =
            QueryAtom::new(AtomId::new(42), AtomSpec::new(Element::C).with_atom_map(0));
        assert!(molecule.validate().is_err());
        assert!(!is_agent(&molecule, 3.0));
    }
}

#[cfg(test)]
mod complete_remove_unmapped_reactants_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;
    fn template(label: &str, mapped: bool) -> QueryGraph {
        let spec = if mapped {
            AtomSpec::new(Element::C).with_atom_map(0)
        } else {
            AtomSpec::new(Element::C)
        };
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), spec)],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn labels(rows: &[QueryGraph]) -> Vec<PropertyValue> {
        rows.iter()
            .map(|g| g.prop("label").unwrap().clone())
            .collect()
    }
    fn expected(rows: &[&str]) -> Vec<PropertyValue> {
        rows.iter().map(|s| PropertyValue::from(*s)).collect()
    }
    fn reaction() -> Reaction {
        let templates = vec![
            template("k0", true),
            template("r0", false),
            template("k1", true),
            template("r1", false),
        ];
        let mut reaction = Reaction::from_templates(
            templates,
            vec![template("other", true)],
            vec![template("existing", true)],
        )
        .unwrap();
        reaction.implicit_properties = true;
        reaction.match_params.max_matches = 0;
        reaction
    }
    #[test]
    fn source_options_append_to_existing_outputs_and_preserve_kept_storage_flags_and_order() {
        for needs_init in [false, true] {
            for move_to_agents in [false, true] {
                for output in [false, true] {
                    let mut reaction = reaction();
                    reaction.needs_init = needs_init;
                    let storage = reaction.reactants.as_ptr();
                    let capacity = reaction.reactants.capacity();
                    let kept0 = reaction.reactants[0].atoms().as_ptr();
                    let kept1 = reaction.reactants[2].atoms().as_ptr();
                    let other = reaction.products.as_ptr();
                    let mut target = vec![template("prefix", true)];
                    remove_unmapped_reactants_source(
                        &mut reaction,
                        0.5,
                        move_to_agents,
                        if output { Some(&mut target) } else { None },
                    );
                    assert_eq!(labels(&reaction.reactants), expected(&["k0", "k1"]));
                    assert_eq!(reaction.reactants.as_ptr(), storage);
                    assert_eq!(reaction.reactants.capacity(), capacity);
                    assert_eq!(reaction.reactants[0].atoms().as_ptr(), kept0);
                    assert_eq!(reaction.reactants[1].atoms().as_ptr(), kept1);
                    assert_eq!(reaction.products.as_ptr(), other);
                    assert_eq!(reaction.needs_init, needs_init);
                    assert!(reaction.implicit_properties);
                    assert_eq!(reaction.match_params.max_matches, 0);
                    assert_eq!(
                        labels(&reaction.agents),
                        expected(if move_to_agents {
                            &["existing", "r0", "r1"]
                        } else {
                            &["existing"]
                        })
                    );
                    assert_eq!(
                        labels(&target),
                        expected(if output {
                            &["prefix", "r0", "r1"]
                        } else {
                            &["prefix"]
                        })
                    );
                }
            }
        }
    }
    #[test]
    fn immutable_projection_keeps_original_and_delegates_complete_source_removal() {
        let mut source = reaction();
        source.needs_init = false;
        let result = without_unmapped_reactants(
            &source,
            &ReactionTemplateRemovalParams {
                threshold_unmapped_atoms: 0.5,
                move_to_agent_templates: true,
            },
        );
        assert_eq!(
            labels(&source.reactants),
            expected(&["k0", "r0", "k1", "r1"])
        );
        assert_eq!(labels(&source.agents), expected(&["existing"]));
        assert_eq!(labels(&result.reaction.reactants), expected(&["k0", "k1"]));
        assert_eq!(
            labels(&result.reaction.agents),
            expected(&["existing", "r0", "r1"])
        );
        assert_eq!(labels(&result.removed_templates), expected(&["r0", "r1"]));
        assert!(!result.reaction.needs_init);
        assert!(!source.needs_init);
    }
    #[test]
    fn source_accepts_nonfinite_threshold_and_does_not_add_graph_validation_or_init_reset() {
        let mut reaction = reaction();
        reaction.needs_init = false;
        reaction.reactants[1].atoms_mut()[0] =
            QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        assert!(reaction.reactants[1].validate().is_err());
        let storage = reaction.reactants.as_ptr();
        let capacity = reaction.reactants.capacity();
        remove_unmapped_reactants_source(&mut reaction, f64::NAN, false, None);
        assert!(reaction.reactants.is_empty());
        assert_eq!(reaction.reactants.as_ptr(), storage);
        assert_eq!(reaction.reactants.capacity(), capacity);
        assert_eq!(labels(&reaction.agents), expected(&["existing"]));
        assert!(!reaction.needs_init);
    }
}

#[cfg(test)]
mod complete_remove_unmapped_products_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;
    fn template(label: &str, mapped: bool) -> QueryGraph {
        let spec = if mapped {
            AtomSpec::new(Element::C).with_atom_map(0)
        } else {
            AtomSpec::new(Element::C)
        };
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), spec)],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn labels(rows: &[QueryGraph]) -> Vec<PropertyValue> {
        rows.iter()
            .map(|g| g.prop("label").unwrap().clone())
            .collect()
    }
    fn expected(rows: &[&str]) -> Vec<PropertyValue> {
        rows.iter().map(|s| PropertyValue::from(*s)).collect()
    }
    fn reaction() -> Reaction {
        let templates = vec![
            template("k0", true),
            template("r0", false),
            template("k1", true),
            template("r1", false),
        ];
        let mut reaction = Reaction::from_templates(
            vec![template("other", true)],
            templates,
            vec![template("existing", true)],
        )
        .unwrap();
        reaction.implicit_properties = true;
        reaction.match_params.max_matches = 0;
        reaction
    }
    #[test]
    fn source_options_append_to_existing_outputs_and_preserve_kept_storage_flags_and_order() {
        for needs_init in [false, true] {
            for move_to_agents in [false, true] {
                for output in [false, true] {
                    let mut reaction = reaction();
                    reaction.needs_init = needs_init;
                    let storage = reaction.products.as_ptr();
                    let capacity = reaction.products.capacity();
                    let kept0 = reaction.products[0].atoms().as_ptr();
                    let kept1 = reaction.products[2].atoms().as_ptr();
                    let other = reaction.reactants.as_ptr();
                    let mut target = vec![template("prefix", true)];
                    remove_unmapped_products_source(
                        &mut reaction,
                        0.5,
                        move_to_agents,
                        if output { Some(&mut target) } else { None },
                    );
                    assert_eq!(labels(&reaction.products), expected(&["k0", "k1"]));
                    assert_eq!(reaction.products.as_ptr(), storage);
                    assert_eq!(reaction.products.capacity(), capacity);
                    assert_eq!(reaction.products[0].atoms().as_ptr(), kept0);
                    assert_eq!(reaction.products[1].atoms().as_ptr(), kept1);
                    assert_eq!(reaction.reactants.as_ptr(), other);
                    assert_eq!(reaction.needs_init, needs_init);
                    assert!(reaction.implicit_properties);
                    assert_eq!(reaction.match_params.max_matches, 0);
                    assert_eq!(
                        labels(&reaction.agents),
                        expected(if move_to_agents {
                            &["existing", "r0", "r1"]
                        } else {
                            &["existing"]
                        })
                    );
                    assert_eq!(
                        labels(&target),
                        expected(if output {
                            &["prefix", "r0", "r1"]
                        } else {
                            &["prefix"]
                        })
                    );
                }
            }
        }
    }
    #[test]
    fn immutable_projection_keeps_original_and_delegates_complete_source_removal() {
        let mut source = reaction();
        source.needs_init = false;
        let result = without_unmapped_products(
            &source,
            &ReactionTemplateRemovalParams {
                threshold_unmapped_atoms: 0.5,
                move_to_agent_templates: true,
            },
        );
        assert_eq!(
            labels(&source.products),
            expected(&["k0", "r0", "k1", "r1"])
        );
        assert_eq!(labels(&source.agents), expected(&["existing"]));
        assert_eq!(labels(&result.reaction.products), expected(&["k0", "k1"]));
        assert_eq!(
            labels(&result.reaction.agents),
            expected(&["existing", "r0", "r1"])
        );
        assert_eq!(labels(&result.removed_templates), expected(&["r0", "r1"]));
        assert!(!result.reaction.needs_init);
        assert!(!source.needs_init);
    }
    #[test]
    fn source_accepts_nonfinite_threshold_and_does_not_add_graph_validation_or_init_reset() {
        let mut reaction = reaction();
        reaction.needs_init = false;
        reaction.products[1].atoms_mut()[0] =
            QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        assert!(reaction.products[1].validate().is_err());
        let storage = reaction.products.as_ptr();
        let capacity = reaction.products.capacity();
        remove_unmapped_products_source(&mut reaction, f64::NAN, false, None);
        assert!(reaction.products.is_empty());
        assert_eq!(reaction.products.as_ptr(), storage);
        assert_eq!(reaction.products.capacity(), capacity);
        assert_eq!(labels(&reaction.agents), expected(&["existing"]));
        assert!(!reaction.needs_init);
    }
}

#[cfg(test)]
mod complete_remove_agents_source {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, PropertyText, PropertyValue, QueryAtom};
    use cosmolkit_types::Element;
    fn template(label: &str) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [(PropertyText::from("label"), PropertyValue::from(label))],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn labels(rows: &[QueryGraph]) -> Vec<PropertyValue> {
        rows.iter()
            .map(|g| g.prop("label").unwrap().clone())
            .collect()
    }
    fn expected(rows: &[&str]) -> Vec<PropertyValue> {
        rows.iter().map(|s| PropertyValue::from(*s)).collect()
    }
    #[test]
    fn source_optional_target_appends_before_clear_preserving_vector_capacity_and_other_state() {
        for needs_init in [false, true] {
            for output in [false, true] {
                let mut reaction = Reaction::from_templates(
                    vec![template("r")],
                    vec![template("p")],
                    vec![template("a0"), template("a1")],
                )
                .unwrap();
                reaction.needs_init = needs_init;
                reaction.implicit_properties = true;
                reaction.match_params.max_matches = 0;
                let storage = reaction.agents.as_ptr();
                let capacity = reaction.agents.capacity();
                let reactants = reaction.reactants.as_ptr();
                let products = reaction.products.as_ptr();
                let mut target = vec![template("prefix")];
                remove_agents_source(&mut reaction, if output { Some(&mut target) } else { None });
                assert!(reaction.agents.is_empty());
                assert_eq!(reaction.agents.as_ptr(), storage);
                assert_eq!(reaction.agents.capacity(), capacity);
                assert_eq!(reaction.reactants.as_ptr(), reactants);
                assert_eq!(reaction.products.as_ptr(), products);
                assert_eq!(reaction.needs_init, needs_init);
                assert!(reaction.implicit_properties);
                assert_eq!(reaction.match_params.max_matches, 0);
                assert_eq!(
                    labels(&target),
                    expected(if output {
                        &["prefix", "a0", "a1"]
                    } else {
                        &["prefix"]
                    })
                );
            }
        }
    }
    #[test]
    fn immutable_projection_preserves_original_templates_and_initialization() {
        let mut source = Reaction::from_templates(
            vec![template("r")],
            vec![template("p")],
            vec![template("a0"), template("a1")],
        )
        .unwrap();
        source.needs_init = false;
        let result = without_agents(&source);
        assert_eq!(labels(&source.agents), expected(&["a0", "a1"]));
        assert!(result.reaction.agents.is_empty());
        assert_eq!(labels(&result.removed_templates), expected(&["a0", "a1"]));
        assert_eq!(result.reaction.reactants, source.reactants);
        assert_eq!(result.reaction.products, source.products);
        assert!(!source.needs_init);
        assert!(!result.reaction.needs_init);
    }
    #[test]
    fn source_empty_role_keeps_target_prefix_and_unvalidated_template_state_is_copied() {
        let mut reaction = Reaction::new();
        let mut target = vec![template("prefix")];
        remove_agents_source(&mut reaction, Some(&mut target));
        assert_eq!(labels(&target), expected(&["prefix"]));
        let mut graph = template("invalid");
        graph.atoms_mut()[0] = QueryAtom::new(AtomId::new(99), AtomSpec::new(Element::C));
        reaction.agents.push(graph);
        remove_agents_source(&mut reaction, Some(&mut target));
        assert!(reaction.agents.is_empty());
        assert_eq!(target[1].atoms()[0].id().index(), 99);
        assert!(target[1].validate().is_err());
    }
}
