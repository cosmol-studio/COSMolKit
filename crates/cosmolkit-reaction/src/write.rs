use crate::{Reaction, ReactionRole, ReactionWriteError, ReactionWriteParams};
use cosmolkit_model::{PropertyText, QueryGraph};
use cosmolkit_search::{SmartsWriteOutput, SmartsWriteParams};

pub(crate) struct SmirksBase {
    pub text: PropertyText,
    /// Writer-produced atom/bond orders retained in original R>A>P order.
    pub outputs: Vec<SmartsWriteOutput>,
}

fn template_to_string(
    graph: &QueryGraph,
    role: ReactionRole,
    index: usize,
    params: &ReactionWriteParams,
) -> Result<SmartsWriteOutput, ReactionWriteError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionWriter.cpp :: molToString
    // RDKit❗✔️: std::string molToString(RDKit::ROMol &mol, bool toSmiles,
    // RDKit❗✔️:                         const RDKit::SmilesWriteParams &params) {
    // RDKit❗✔️:   std::string res = "";
    // RDKit❗✔️:   if (toSmiles) {
    // RDKit❗✔️:     res = MolToSmiles(mol, params);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     res = MolToSmarts(mol, params);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   std::vector<int> mapping;
    // RDKit❗✔️:   if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    // RDKit❗✔️:     res = "(" + res + ")";
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    let smarts_params = SmartsWriteParams {
        include_atom_maps: true,
        do_isomeric_smiles: params.do_isomeric_smiles,
        include_dative_bonds: params.include_dative_bonds,
        rooted_at_atom: params.rooted_at_atom,
    };
    let mut output = cosmolkit_search::query_graph_to_smarts_output(graph, &smarts_params)
        .map_err(|source| ReactionWriteError::Template {
            role,
            template: index,
            source,
        })?;
    let components = cosmolkit_core::query_connected_components(graph).map_err(|source| {
        ReactionWriteError::Connectivity {
            role,
            template: index,
            source,
        }
    })?;
    if components.components.len() > 1 {
        let mut grouped = PropertyText::with_capacity(output.text.len() + 2);
        grouped.push_byte(b'(');
        grouped.extend_bytes(output.text.as_bytes());
        grouped.push_byte(b')');
        output.text = grouped;
    }
    Ok(output)
}

fn role_to_string(
    graphs: &[QueryGraph],
    role: ReactionRole,
    params: &ReactionWriteParams,
) -> Result<(PropertyText, Vec<SmartsWriteOutput>), ReactionWriteError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionWriter.cpp :: chemicalReactionTemplatesToString
    // RDKit❗✔️: std::string chemicalReactionTemplatesToString(
    // RDKit❗✔️:     const RDKit::ChemicalReaction &rxn, RDKit::ReactionMoleculeType type,
    // RDKit❗✔️:     bool toSmiles, const RDKit::SmilesWriteParams &params) {
    // RDKit❗✔️:   std::string res = "";
    // RDKit❗✔️:   std::vector<std::string> vfragsmi;
    // RDKit❗✔️:   auto begin = getStartIterator(rxn, type);
    // RDKit❗✔️:   auto end = getEndIterator(rxn, type);
    // RDKit❗✔️:   for (; begin != end; ++begin) {
    // RDKit❗✔️:     vfragsmi.push_back(molToString(**begin, toSmiles, params));
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (params.canonical) {
    // RDKit❗✔️:     std::sort(vfragsmi.begin(), vfragsmi.end());
    // RDKit❗✔️:   }
    // RDKit❗✔️:   for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit❗✔️:     res += vfragsmi[i];
    // RDKit❗✔️:     if (i < vfragsmi.size() - 1) {
    // RDKit❗✔️:       res += ".";
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Source same ordered writes, copied per-template strings and optional
    // lexical role sorting. Output evidence is retained in insertion order.
    let mut outputs = Vec::with_capacity(graphs.len());
    let mut text = Vec::with_capacity(graphs.len());
    for (index, graph) in graphs.iter().enumerate() {
        let output = template_to_string(graph, role, index, params)?;
        text.push(output.text.clone());
        outputs.push(output);
    }
    if params.canonical {
        text.sort_unstable();
    }
    // Counted byte concatenation preserves source sorting and every payload
    // byte; no Unicode conversion or Display adapter participates in chemistry.
    let capacity = text.iter().map(PropertyText::len).sum::<usize>() + text.len().saturating_sub(1);
    let mut joined = PropertyText::with_capacity(capacity);
    for (index, fragment) in text.iter().enumerate() {
        if index != 0 {
            joined.push_byte(b'.');
        }
        joined.extend_bytes(fragment.as_bytes());
    }
    Ok((joined, outputs))
}

pub(crate) fn smirks_base(
    reaction: &Reaction,
    params: &ReactionWriteParams,
) -> Result<SmirksBase, ReactionWriteError> {
    // BEGIN RDKIT COMPLETE CPP FUNCTION: Code/GraphMol/ChemReactions/ReactionWriter.cpp :: chemicalReactionToRxnToString
    // RDKit❗✔️: std::string chemicalReactionToRxnToString(
    // RDKit❗✔️:     const RDKit::ChemicalReaction &rxn, bool toSmiles,
    // RDKit❗✔️:     const RDKit::SmilesWriteParams &params, bool includeCX, std::uint32_t flags = RDKit::SmilesWrite::CXSmilesFields::CX_ALL) {
    // RDKit❗✔️:   std::string res = "";
    // RDKit❗✔️:   res +=
    // RDKit❗✔️:       chemicalReactionTemplatesToString(rxn, RDKit::Reactant, toSmiles, params);
    // RDKit❗✔️:   res += ">";
    // RDKit❗✔️:   res += chemicalReactionTemplatesToString(rxn, RDKit::Agent, toSmiles, params);
    // RDKit❗✔️:   res += ">";
    // RDKit❗✔️:   res +=
    // RDKit❗✔️:       chemicalReactionTemplatesToString(rxn, RDKit::Product, toSmiles, params);
    // RDKit❗✔️:
    // RDKit❗✔️:
    // RDKit❗✔️:   if (includeCX) {
    // RDKit❗✔️:     std::vector<RDKit::ROMol *> mols;
    // RDKit❗✔️:
    // RDKit❗✔️:     // Collect reactants, agents, and products into mols vector
    // RDKit❗✔️:     for (auto type : {RDKit::Reactant, RDKit::Agent, RDKit::Product}) {
    // RDKit❗✔️:       for (auto begin = getStartIterator(rxn, type);
    // RDKit❗✔️:            begin != getEndIterator(rxn, type); ++begin) {
    // RDKit❗✔️:         mols.push_back((*begin).get());
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     auto ext = RDKit::SmilesWrite::getCXExtensions(mols, flags);
    // RDKit❗✔️:     if (!ext.empty()) {
    // RDKit❗✔️:       res += " ";
    // RDKit❗✔️:       res += ext;
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT COMPLETE CPP FUNCTION
    // Source SMARTS half: read templates without lifecycle initialization or
    // a query-to-molecule conversion. The CX vector half consumes these exact
    // original-order outputs after role-string sorting, at its owner boundary.
    let (reactants, mut outputs) = role_to_string(
        reaction.reactant_templates(),
        ReactionRole::Reactant,
        params,
    )?;
    let (agents, mut agent_outputs) =
        role_to_string(reaction.agent_templates(), ReactionRole::Agent, params)?;
    let (products, mut product_outputs) =
        role_to_string(reaction.product_templates(), ReactionRole::Product, params)?;
    outputs.append(&mut agent_outputs);
    outputs.append(&mut product_outputs);
    let mut text = PropertyText::with_capacity(reactants.len() + agents.len() + products.len() + 2);
    text.extend_bytes(reactants.as_bytes());
    text.push_byte(b'>');
    text.extend_bytes(agents.as_bytes());
    text.push_byte(b'>');
    text.extend_bytes(products.as_bytes());
    Ok(SmirksBase { text, outputs })
}

#[cfg(test)]
mod byte_caller_tests {
    use super::smirks_base;
    use crate::{Reaction, ReactionWriteParams, parse_smirks};

    #[test]
    fn source_role_sorting_preserves_original_writer_evidence_order() {
        let reaction = parse_smirks("N.C>O>C").unwrap();
        let original = reaction.clone();
        let plain = smirks_base(&reaction, &ReactionWriteParams::default()).unwrap();
        assert_eq!(plain.text.as_bytes(), b"N.C>O>C");
        let canonical = smirks_base(
            &reaction,
            &ReactionWriteParams {
                canonical: true,
                ..ReactionWriteParams::default()
            },
        )
        .unwrap();
        assert_eq!(canonical.text.as_bytes(), b"C.N>O>C");
        assert_eq!(canonical.outputs.len(), 4);
        assert_eq!(canonical.outputs[0].text.as_bytes(), b"N");
        assert_eq!(canonical.outputs[1].text.as_bytes(), b"C");
        assert_eq!(canonical.outputs[2].text.as_bytes(), b"O");
        assert_eq!(canonical.outputs[3].text.as_bytes(), b"C");
        assert_eq!(reaction.reactants, original.reactants);
        assert_eq!(reaction.agents, original.agents);
        assert_eq!(reaction.products, original.products);
        assert_eq!(reaction.needs_init, original.needs_init);
        assert_eq!(reaction.implicit_properties, original.implicit_properties);
        assert_eq!(reaction.properties, original.properties);
    }

    #[test]
    fn source_disconnected_template_grouping_and_empty_roles_keep_delimiters() {
        let reaction = parse_smirks("(C.N)>>C").unwrap();
        let output = smirks_base(&reaction, &ReactionWriteParams::default()).unwrap();
        assert_eq!(output.text.as_bytes(), b"(C.N)>>C");
        assert_eq!(output.outputs.len(), 2);
        assert_eq!(output.outputs[0].atom_order.len(), 2);
        let empty = smirks_base(&Reaction::new(), &ReactionWriteParams::default()).unwrap();
        assert_eq!(empty.text.as_bytes(), b">>");
        assert!(empty.outputs.is_empty());
    }
}
