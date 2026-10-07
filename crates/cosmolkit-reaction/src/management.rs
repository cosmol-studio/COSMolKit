use crate::{Reaction, ReactionTemplateRemovalParams};
use cosmolkit_model::QueryGraph;

#[derive(Debug, Clone)]
pub struct ReactionTemplateRemoval {
    pub reaction: Reaction,
    pub removed_templates: Vec<QueryGraph>,
}

fn is_agent(graph: &QueryGraph, threshold: f64) -> bool {
    // RDKit❗✔️: bool isReactionTemplateMoleculeAgent(const ROMol &mol, double agentThreshold) {
    // RDKit❗✔️:   unsigned numMappedAtoms = MolOps::getNumAtomsWithDistinctProperty(
    // RDKit❗✔️:       mol, common_properties::molAtomMapNumber);
    // RDKit❗✔️:   unsigned numAtoms = mol.getNumHeavyAtoms();
    // RDKit❗✔️:   return !(numAtoms > 0u && static_cast<double>(numMappedAtoms) /
    // RDKit❗✔️:                                     static_cast<double>(numAtoms) >=
    // RDKit❗✔️:                                 agentThreshold);
    // RDKit❗✔️: unsigned getNumAtomsWithDistinctProperty(const ROMol &mol,
    // RDKit❗✔️:                                          const std::string_view &prop) {
    // RDKit❗✔️:   unsigned numPropAtoms = 0;
    // RDKit❗✔️:   for (const auto atom : mol.atoms()) {
    // RDKit❗✔️:     if (atom->hasProp(prop)) {
    // RDKit❗✔️:       ++numPropAtoms;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return numPropAtoms;
    // RDKit❗✔️: unsigned int ROMol::getNumHeavyAtoms() const {
    // RDKit❗✔️:   unsigned int res = 0;
    // RDKit❗✔️:   for (const auto atom : atoms()) {
    // RDKit❗✔️:     if (atom->getAtomicNum() > 1) {
    // RDKit❗✔️:       ++res;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // Presence includes a canonical map of zero and a raw property of any type.
    // Count each atom once even if both carriers store the name.
    let mapped = graph
        .atoms()
        .iter()
        .filter(|atom| atom.atom_map().is_some() || atom.prop("molAtomMapNumber").is_some())
        .count();
    let heavy = graph
        .atoms()
        .iter()
        .filter(|atom| atom.atomic_number() > 1)
        .count();
    !(heavy > 0 && (mapped as f64) / (heavy as f64) >= threshold)
}

#[doc(hidden)]
pub fn without_unmapped_reactants(
    reaction: &Reaction,
    params: &ReactionTemplateRemovalParams,
) -> ReactionTemplateRemoval {
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
    // D1 immutable value projection; retain the source needsInit state.
    let mut result = reaction.clone();
    let mut removed = Vec::new();
    let templates = std::mem::take(&mut result.reactants);
    for template in templates {
        if is_agent(&template, params.threshold_unmapped_atoms) {
            if params.move_to_agent_templates {
                result.agents.push(template.clone());
            }
            removed.push(template);
        } else {
            result.reactants.push(template);
        }
    }
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates: removed,
    }
}

#[doc(hidden)]
pub fn without_unmapped_products(
    reaction: &Reaction,
    params: &ReactionTemplateRemovalParams,
) -> ReactionTemplateRemoval {
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
    // D1 immutable value projection; retain the source needsInit state.
    let mut result = reaction.clone();
    let mut removed = Vec::new();
    let templates = std::mem::take(&mut result.products);
    for template in templates {
        if is_agent(&template, params.threshold_unmapped_atoms) {
            if params.move_to_agent_templates {
                result.agents.push(template.clone());
            }
            removed.push(template);
        } else {
            result.products.push(template);
        }
    }
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates: removed,
    }
}

#[doc(hidden)]
pub fn without_agents(reaction: &Reaction) -> ReactionTemplateRemoval {
    // RDKit❗❌: void ChemicalReaction::removeAgentTemplates(MOL_SPTR_VECT *targetVector) {
    // RDKit❗❌:   if (targetVector) {
    // RDKit❗❌:     for (auto iter = beginAgentTemplates(); iter != endAgentTemplates();
    // RDKit❗❌:          ++iter) {
    // RDKit❗❌:       targetVector->push_back(*iter);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   m_agentTemplates.clear();
    let mut result = reaction.clone();
    let removed_templates = std::mem::take(&mut result.agents);
    ReactionTemplateRemoval {
        reaction: result,
        removed_templates,
    }
}
