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
    // RDKit❗❌: std::string molToString(RDKit::ROMol &mol, bool toSmiles,
    // RDKit❗❌:                         const RDKit::SmilesWriteParams &params) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   if (toSmiles) {
    // RDKit❌❌:     res = MolToSmiles(mol, params);
    // RDKit❗❌:   } else {
    // RDKit❗❌:     res = MolToSmarts(mol, params);
    // RDKit❗❌:   }
    // RDKit❗❌:   std::vector<int> mapping;
    // RDKit❗❌:   if (RDKit::MolOps::getMolFrags(mol, mapping) > 1) {
    // RDKit❗❌:     res = "(" + res + ")";
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // This existing SMIRKS-only type boundary fixes source toSmiles=false.
    // The independent MolToSmiles capability is unmodeled (❌❌ above), not
    // accepted and silently redirected to SMARTS. No new format/API is added.
    // Reuse the now-complete SMARTS output owner before source getMolFrags.
    // Count actual graph components, never delimiter bytes in emitted text.
    // Cost ❌: CORE checked query component projection validates model state
    // and materializes component member vectors beyond Native's label vector;
    // Vec-backed short text versus Native SSO also costs extra allocations.
    let smarts_params = SmartsWriteParams {
        include_atom_maps: true,
        isomeric_smiles: params.isomeric_smiles,
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
        // Behavior: source concatenation copies counted bytes, including NUL
        // and non-UTF-8 payloads, between the two literal delimiters.
        // One O(bytes) counted copy. Native chained concatenation has mixed
        // SSO/long-buffer costs; do not claim identical allocation behavior.
        let mut text = PropertyText::with_capacity(output.text.len() + 2);
        text.push_byte(b'(');
        text.extend_bytes(output.text.as_bytes());
        text.push_byte(b')');
        output.text = text;
    }
    Ok(output)
}

fn get_start_iterator(reaction: &Reaction, role: ReactionRole) -> std::slice::Iter<'_, QueryGraph> {
    // RDKit❗✔️: MOL_SPTR_VECT::const_iterator getStartIterator(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                                ReactionMoleculeType t) {
    // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator begin;
    // RDKit❗✔️:   if (t == Reactant) {
    // RDKit❗✔️:     begin = rxn.beginReactantTemplates();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (t == Product) {
    // RDKit❗✔️:     begin = rxn.beginProductTemplates();
    // RDKit❗✔️:     ;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (t == Agent) {
    // RDKit❗✔️:     begin = rxn.beginAgentTemplates();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return begin;
    // RDKit❗✔️: }
    // Each valid source role selects exactly its own borrowed vector begin.
    // Rust's closed role enum excludes invalid C++ enum values for which the
    // source returns an uninitialized iterator. No default role is inferred.
    // Cost: O(1), no query copies, allocations, or lifecycle changes.
    match role {
        ReactionRole::Reactant => reaction.reactant_templates().iter(),
        ReactionRole::Product => reaction.product_templates().iter(),
        ReactionRole::Agent => reaction.agent_templates().iter(),
    }
}

fn get_end_iterator(reaction: &Reaction, role: ReactionRole) -> usize {
    // RDKit❗✔️: MOL_SPTR_VECT::const_iterator getEndIterator(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                              ReactionMoleculeType t) {
    // RDKit❗✔️:   MOL_SPTR_VECT::const_iterator end;
    // RDKit❗✔️:   if (t == Reactant) {
    // RDKit❗✔️:     end = rxn.endReactantTemplates();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (t == Product) {
    // RDKit❗✔️:     end = rxn.endProductTemplates();
    // RDKit❗✔️:     ;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (t == Agent) {
    // RDKit❗✔️:     end = rxn.endAgentTemplates();
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return end;
    // RDKit❗✔️: }
    // The end sentinel is the selected vector's one-past-last ordinal.
    // get_start_iterator starts at ordinal zero; no dereference of the end
    // sentinel occurs. Closed valid role dispatch has no source-absent default.
    // Cost: O(1), no iteration, allocation, query clone, or initialization.
    match role {
        ReactionRole::Reactant => reaction.reactant_templates().len(),
        ReactionRole::Product => reaction.product_templates().len(),
        ReactionRole::Agent => reaction.agent_templates().len(),
    }
}

fn role_to_string(
    reaction: &Reaction,
    role: ReactionRole,
    params: &ReactionWriteParams,
) -> Result<(PropertyText, Vec<SmartsWriteOutput>), ReactionWriteError> {
    // RDKit❗❌: std::string chemicalReactionTemplatesToString(
    // RDKit❗❌:     const RDKit::ChemicalReaction &rxn, RDKit::ReactionMoleculeType type,
    // RDKit❗❌:     bool toSmiles, const RDKit::SmilesWriteParams &params) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   std::vector<std::string> vfragsmi;
    // RDKit❗❌:   auto begin = getStartIterator(rxn, type);
    // RDKit❗❌:   auto end = getEndIterator(rxn, type);
    // RDKit❗❌:   for (; begin != end; ++begin) {
    // RDKit❗❌:     vfragsmi.push_back(molToString(**begin, toSmiles, params));
    // RDKit❗❌:   }
    // RDKit❗❌:   if (params.canonical) {
    // RDKit❗❌:     std::sort(vfragsmi.begin(), vfragsmi.end());
    // RDKit❗❌:   }
    // RDKit❗❌:   for (unsigned i = 0; i < vfragsmi.size(); ++i) {
    // RDKit❗❌:     res += vfragsmi[i];
    // RDKit❗❌:     if (i < vfragsmi.size() - 1) {
    // RDKit❗❌:       res += ".";
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // The existing SMIRKS boundary fixes source toSmiles=false. Both native
    // role iterators select the same vector; the Rust end is its one-past-last
    // ordinal. Emit every template before sorting and preserve error order.
    let begin = get_start_iterator(reaction, role);
    let end = get_end_iterator(reaction, role);
    let mut outputs = Vec::with_capacity(end);
    let mut text = Vec::with_capacity(end);
    for (index, graph) in begin.enumerate().take(end) {
        let output = template_to_string(graph, role, index, params)?;
        text.push(output.text.clone());
        outputs.push(output);
    }
    // PropertyText compares unsigned counted bytes like std::string traits;
    // embedded NUL and non-UTF-8 payloads neither terminate nor decode text.
    // Only text is sorted: Native stores output-order props on the original
    // templates, represented here by original-order detached writer evidence.
    if params.canonical {
        text.sort_unstable();
    }
    let mut joined = PropertyText::new();
    for (index, fragment) in text.iter().enumerate() {
        joined.extend_bytes(fragment.as_bytes());
        if index + 1 < text.len() {
            joined.push_byte(b'.');
        }
    }
    // Cost ❌: retaining output evidence adds a second template-text copy and
    // vector beyond Native vfragsmi; Vec short text also lacks Native SSO.
    // Sorting stays O(n log n) byte comparisons, joining O(total bytes).
    Ok((joined, outputs))
}

pub(crate) fn smirks_base(
    reaction: &Reaction,
    params: &ReactionWriteParams,
    include_cx: bool,
    flags: cosmolkit_smiles::CxSmilesFields,
) -> Result<SmirksBase, ReactionWriteError> {
    // RDKit❗❌: std::string chemicalReactionToRxnToString(
    // RDKit❗❌:     const RDKit::ChemicalReaction &rxn, bool toSmiles,
    // RDKit❗❌:     const RDKit::SmilesWriteParams &params, bool includeCX, std::uint32_t flags = RDKit::SmilesWrite::CXSmilesFields::CX_ALL) {
    // RDKit❗❌:   std::string res = "";
    // RDKit❗❌:   res +=
    // RDKit❗❌:       chemicalReactionTemplatesToString(rxn, RDKit::Reactant, toSmiles, params);
    // RDKit❗❌:   res += ">";
    // RDKit❗❌:   res += chemicalReactionTemplatesToString(rxn, RDKit::Agent, toSmiles, params);
    // RDKit❗❌:   res += ">";
    // RDKit❗❌:   res +=
    // RDKit❗❌:       chemicalReactionTemplatesToString(rxn, RDKit::Product, toSmiles, params);
    // RDKit❗❌:
    // RDKit❗❌:
    // RDKit❗❌:   if (includeCX) {
    // RDKit❗❌:     std::vector<RDKit::ROMol *> mols;
    // RDKit❗❌:
    // RDKit❗❌:     // Collect reactants, agents, and products into mols vector
    // RDKit❗❌:     for (auto type : {RDKit::Reactant, RDKit::Agent, RDKit::Product}) {
    // RDKit❗❌:       for (auto begin = getStartIterator(rxn, type);
    // RDKit❗❌:            begin != getEndIterator(rxn, type); ++begin) {
    // RDKit❗❌:         mols.push_back((*begin).get());
    // RDKit❗❌:       }
    // RDKit❗❌:     }
    // RDKit❗❌:
    // RDKit❗❌:     auto ext = RDKit::SmilesWrite::getCXExtensions(mols, flags);
    // RDKit❗❌:     if (!ext.empty()) {
    // RDKit❗❌:       res += " ";
    // RDKit❗❌:       res += ext;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // This established SMIRKS-only boundary specializes toSmiles=false;
    // independent reaction SMILES stays unmodeled, never silently redirected.
    // Source role strings are completed in R>A>P before ANY CX preflight.
    let (reactants, mut outputs) = role_to_string(reaction, ReactionRole::Reactant, params)?;
    let (agents, mut agent_outputs) = role_to_string(reaction, ReactionRole::Agent, params)?;
    let (products, mut product_outputs) = role_to_string(reaction, ReactionRole::Product, params)?;
    outputs.append(&mut agent_outputs);
    outputs.append(&mut product_outputs);
    let mut text = PropertyText::with_capacity(reactants.len() + agents.len() + products.len() + 2);
    text.extend_bytes(reactants.as_bytes());
    text.push_byte(b'>');
    text.extend_bytes(agents.as_bytes());
    text.push_byte(b'>');
    text.extend_bytes(products.as_bytes());
    if include_cx {
        let mut templates = Vec::with_capacity(outputs.len());
        for role in [
            ReactionRole::Reactant,
            ReactionRole::Agent,
            ReactionRole::Product,
        ] {
            for graph in get_start_iterator(reaction, role).take(get_end_iterator(reaction, role)) {
                templates.push(graph);
            }
        }
        // Original template order is independent of canonical role-text sort.
        // Consume exact typed writer getter values; the vector CX owner alone
        // inserts molecules, offsets orders and invokes the complete writer.
        let templates: Vec<_> = templates.into_iter().zip(&outputs).collect();
        let composition = cosmolkit_search::compose_query_cx_templates(&templates, flags)
            .map_err(ReactionWriteError::Cx)?;
        if !composition.extension.is_empty() {
            text.push_byte(b' ');
            text.extend_bytes(composition.extension.as_bytes());
        }
    }
    // Cost ❌: detached per-template output evidence and zipped reference Vec
    // add buffers beyond Native's raw pointers; no query/graph/AST clone here.
    // Role sort and joining retain source asymptotic costs and counted bytes.
    Ok(SmirksBase { text, outputs })
}

#[cfg(test)]
mod template_to_string_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondQueryPredicate, BondSpec, PropertyValue, QueryAtom,
        QueryBond, QueryNode,
    };
    use cosmolkit_types::{BondOrder, ChiralTag, Element};
    fn graph(elements: &[Element], edges: &[(usize, usize, BondOrder)]) -> QueryGraph {
        QueryGraph::from_parts(
            elements
                .iter()
                .enumerate()
                .map(|(i, &e)| QueryAtom::new(AtomId::new(i), AtomSpec::new(e)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b, o))| {
                    QueryBond::new(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), o),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn write(q: &QueryGraph) -> SmartsWriteOutput {
        template_to_string(q, ReactionRole::Reactant, 0, &Default::default()).unwrap()
    }
    #[test]
    fn empty_graph_stays_empty_and_does_not_claim_source_output_order_properties() {
        let out = write(&graph(&[], &[]));
        assert!(out.text.is_empty());
        assert!(!out.source_orders_written);
        assert!(out.atom_order.is_empty());
    }
    #[test]
    fn disconnected_components_wrap_once_and_keep_exact_source_orders() {
        let q = graph(&[Element::C, Element::O], &[]);
        let before = q.clone();
        let out = write(&q);
        assert_eq!(out.text, PropertyText::from("([#6].[#8])"));
        assert_eq!(out.atom_order, vec![AtomId::new(0), AtomId::new(1)]);
        assert!(out.source_orders_written);
        assert_eq!(q, before);
    }
    #[test]
    fn connected_template_uses_no_outer_template_parentheses() {
        let out = write(&graph(
            &[Element::C, Element::O],
            &[(0, 1, BondOrder::Single)],
        ));
        assert_eq!(out.text, PropertyText::from("[#6]-[#8]"));
        assert_eq!(out.bond_order, vec![BondId::new(0)]);
    }
    #[test]
    fn dot_bytes_in_atom_symbol_never_infer_disconnected_topology() {
        let mut q = graph(&[Element::C], &[]);
        q.atom_mut(0)
            .unwrap()
            .set_prop("smilesSymbol", PropertyValue::String("a.b".into()))
            .unwrap();
        assert_eq!(write(&q).text, PropertyText::from("[a.b]"));
    }
    #[test]
    fn binary_non_utf8_and_nul_symbol_is_wrapped_as_counted_bytes() {
        let mut q = graph(&[Element::C, Element::O], &[]);
        let mut symbol = PropertyText::new();
        symbol.extend_bytes(&[255, 0, b'.']);
        q.atom_mut(0)
            .unwrap()
            .set_prop("smilesSymbol", PropertyValue::String(symbol))
            .unwrap();
        let out = write(&q);
        let mut expected = PropertyText::new();
        expected.extend_bytes(b"([");
        expected.extend_bytes(&[255, 0, b'.']);
        expected.extend_bytes(b"].[#8])");
        assert_eq!(out.text, expected);
        assert_eq!(out.atom_order.len(), 2);
    }
    #[test]
    fn rooted_atom_parameter_changes_source_order_before_template_wrapping() {
        let q = graph(&[Element::C, Element::O], &[]);
        let mut p = ReactionWriteParams::default();
        p.rooted_at_atom = Some(1);
        let out = template_to_string(&q, ReactionRole::Product, 7, &p).unwrap();
        assert_eq!(out.text, PropertyText::from("([#8].[#6])"));
        assert_eq!(out.atom_order, vec![AtomId::new(1), AtomId::new(0)]);
    }
    #[test]
    fn writer_flags_reach_sole_smarts_writer_without_canonicalizing_template_components() {
        let mut q = graph(&[Element::O, Element::C], &[]);
        q.atom_mut(1)
            .unwrap()
            .set_chiral_tag(ChiralTag::TetrahedralCw);
        q.set_prop("_StereochemDone", PropertyValue::Bool(false))
            .unwrap();
        let mut p = ReactionWriteParams::default();
        p.isomeric_smiles = false;
        p.canonical = true;
        let out = template_to_string(&q, ReactionRole::Agent, 2, &p).unwrap();
        assert_eq!(out.text, PropertyText::from("([#8].[#6])"));
    }
    #[test]
    fn source_writer_error_retains_role_and_template_index() {
        let mut q = graph(&[Element::C, Element::O], &[(0, 1, BondOrder::Single)]);
        q.bonds_mut()[0].set_predicate(QueryNode::predicate(BondQueryPredicate::HasStereo));
        assert!(matches!(
            template_to_string(&q, ReactionRole::Agent, 19, &Default::default()),
            Err(ReactionWriteError::Template {
                role: ReactionRole::Agent,
                template: 19,
                source: cosmolkit_search::SmartsWriteError::UnwritableBondQuery {
                    description: "BondStereo"
                }
            })
        ));
    }
    #[test]
    fn source_dative_option_and_atom_to_left_order_are_forwarded() {
        let q = graph(&[Element::C, Element::O], &[(0, 1, BondOrder::Dative)]);
        let mut p = ReactionWriteParams::default();
        p.rooted_at_atom = Some(1);
        assert_eq!(
            template_to_string(&q, ReactionRole::Reactant, 0, &p)
                .unwrap()
                .text,
            PropertyText::from("[#8]<-[#6]")
        );
        p.include_dative_bonds = false;
        assert_eq!(
            template_to_string(&q, ReactionRole::Reactant, 0, &p)
                .unwrap()
                .text,
            PropertyText::from("[#8]-[#6]")
        );
    }
}

#[cfg(test)]
mod get_start_iterator_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    fn graph(e: Element) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(e))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn reaction() -> Reaction {
        Reaction::from_templates(
            vec![graph(Element::C), graph(Element::N)],
            vec![graph(Element::O)],
            vec![graph(Element::F), graph(Element::N), graph(Element::O)],
        )
        .unwrap()
    }
    #[test]
    fn reactant_begin_borrows_selected_vector() {
        let r = reaction();
        let mut i = get_start_iterator(&r, ReactionRole::Reactant);
        assert_eq!(i.len(), 2);
        assert!(std::ptr::eq(i.next().unwrap(), &r.reactant_templates()[0]));
    }
    #[test]
    fn product_begin_borrows_selected_vector() {
        let r = reaction();
        let mut i = get_start_iterator(&r, ReactionRole::Product);
        assert_eq!(i.len(), 1);
        assert!(std::ptr::eq(i.next().unwrap(), &r.product_templates()[0]));
    }
    #[test]
    fn agent_begin_borrows_selected_vector() {
        let r = reaction();
        let mut i = get_start_iterator(&r, ReactionRole::Agent);
        assert_eq!(i.len(), 3);
        assert!(std::ptr::eq(i.next().unwrap(), &r.agent_templates()[0]));
    }
    #[test]
    fn empty_selected_role_does_not_fall_back_to_other_roles() {
        let r = Reaction::from_templates(vec![], vec![graph(Element::C)], vec![]).unwrap();
        assert!(
            get_start_iterator(&r, ReactionRole::Reactant)
                .next()
                .is_none()
        );
        assert!(get_start_iterator(&r, ReactionRole::Agent).next().is_none());
    }
    #[test]
    fn iteration_keeps_physical_order_and_does_not_change_initialization() {
        let r = reaction();
        let before = r.clone();
        let actual: Vec<_> = get_start_iterator(&r, ReactionRole::Reactant).collect();
        assert!(std::ptr::eq(actual[0], &r.reactant_templates()[0]));
        assert!(std::ptr::eq(actual[1], &r.reactant_templates()[1]));
        assert_eq!(r.reactant_templates(), before.reactant_templates());
        assert_eq!(r.product_templates(), before.product_templates());
        assert_eq!(r.agent_templates(), before.agent_templates());
        assert_eq!(r.is_initialized(), before.is_initialized());
        assert_eq!(r.implicit_properties(), before.implicit_properties());
    }
}

#[cfg(test)]
mod get_end_iterator_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, QueryAtom};
    use cosmolkit_types::Element;
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn each_end_matches_its_selected_template_vector() {
        let r = Reaction::from_templates(
            vec![graph()],
            vec![graph(), graph()],
            vec![graph(), graph(), graph()],
        )
        .unwrap();
        for (role, count) in [
            (ReactionRole::Reactant, 1),
            (ReactionRole::Product, 2),
            (ReactionRole::Agent, 3),
        ] {
            assert_eq!(get_end_iterator(&r, role), count);
            assert_eq!(get_start_iterator(&r, role).count(), count);
        }
    }
    #[test]
    fn empty_reaction_has_all_begin_end_pairs_equal() {
        let r = Reaction::new();
        for role in [
            ReactionRole::Reactant,
            ReactionRole::Product,
            ReactionRole::Agent,
        ] {
            assert_eq!(get_end_iterator(&r, role), 0);
            assert!(get_start_iterator(&r, role).next().is_none());
        }
    }
    #[test]
    fn empty_role_end_is_not_another_roles_length() {
        let r = Reaction::from_templates(vec![], vec![graph(), graph()], vec![]).unwrap();
        assert_eq!(get_end_iterator(&r, ReactionRole::Reactant), 0);
        assert_eq!(get_end_iterator(&r, ReactionRole::Agent), 0);
    }
    #[test]
    fn end_lookup_keeps_template_storage_and_lifecycle_unchanged() {
        let r = Reaction::from_templates(vec![graph()], vec![], vec![]).unwrap();
        let before = r.reactant_templates().as_ptr();
        let initialized = r.is_initialized();
        assert_eq!(get_end_iterator(&r, ReactionRole::Reactant), 1);
        assert_eq!(r.reactant_templates().as_ptr(), before);
        assert_eq!(r.is_initialized(), initialized);
    }
}

#[cfg(test)]
mod role_to_string_source_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondQueryPredicate, BondSpec, PropertyValue, QueryAtom,
        QueryBond, QueryNode,
    };
    use cosmolkit_types::{BondOrder, Element};
    fn graph(e: Element) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(e))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn empty() -> QueryGraph {
        QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap()
    }
    fn r(v: Vec<QueryGraph>) -> Reaction {
        Reaction::from_templates(v, vec![], vec![]).unwrap()
    }
    fn write(r: &Reaction, p: &ReactionWriteParams) -> (PropertyText, Vec<SmartsWriteOutput>) {
        role_to_string(r, ReactionRole::Reactant, p).unwrap()
    }
    fn bad() -> QueryGraph {
        let mut b = QueryBond::new(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        b.set_predicate(QueryNode::predicate(BondQueryPredicate::HasStereo));
        QueryGraph::from_parts(
            vec![
                QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
                QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::N)),
            ],
            vec![b],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn every_role_uses_its_canonical_begin_end_pair() {
        let r = Reaction::from_templates(
            vec![graph(Element::C)],
            vec![graph(Element::O)],
            vec![graph(Element::F)],
        )
        .unwrap();
        for (role, text) in [
            (ReactionRole::Reactant, "[#6]"),
            (ReactionRole::Product, "[#8]"),
            (ReactionRole::Agent, "[#9]"),
        ] {
            assert_eq!(
                role_to_string(&r, role, &Default::default()).unwrap().0,
                PropertyText::from(text)
            );
        }
    }
    #[test]
    fn canonical_sort_keeps_duplicate_text_and_original_output_evidence() {
        let r = r(vec![
            graph(Element::O),
            graph(Element::C),
            graph(Element::C),
        ]);
        let mut p = ReactionWriteParams::default();
        p.canonical = true;
        let (text, outputs) = write(&r, &p);
        assert_eq!(text, PropertyText::from("[#6].[#6].[#8]"));
        assert_eq!(
            outputs.iter().map(|o| o.text.clone()).collect::<Vec<_>>(),
            vec![
                PropertyText::from("[#8]"),
                PropertyText::from("[#6]"),
                PropertyText::from("[#6]")
            ]
        );
    }
    #[test]
    fn noncanonical_role_preserves_insertion_order() {
        let r = r(vec![graph(Element::O), graph(Element::C)]);
        let mut p = ReactionWriteParams::default();
        p.canonical = false;
        assert_eq!(write(&r, &p).0, PropertyText::from("[#8].[#6]"));
    }
    #[test]
    fn empty_templates_keep_their_source_separator_positions() {
        let r = r(vec![graph(Element::C), empty(), graph(Element::O)]);
        let mut p = ReactionWriteParams::default();
        p.canonical = false;
        let (text, outputs) = write(&r, &p);
        assert_eq!(text, PropertyText::from("[#6]..[#8]"));
        assert_eq!(outputs.len(), 3);
        assert!(!outputs[1].source_orders_written);
    }
    #[test]
    fn canonical_text_order_compares_unsigned_bytes_after_nul() {
        let symbol = |bytes: &[u8]| {
            let mut q = graph(Element::C);
            let mut text = PropertyText::new();
            text.extend_bytes(bytes);
            q.atom_mut(0)
                .unwrap()
                .set_prop("smilesSymbol", PropertyValue::String(text))
                .unwrap();
            q
        };
        let r = r(vec![symbol(&[0, 255]), symbol(&[0, 127])]);
        let mut p = ReactionWriteParams::default();
        p.canonical = true;
        let mut expected = PropertyText::new();
        expected.extend_bytes(&[b'[', 0, 127, b']', b'.', b'[', 0, 255, b']']);
        assert_eq!(write(&r, &p).0, expected);
    }
    #[test]
    fn template_errors_precede_sorting_and_keep_original_role_index() {
        let r = Reaction::from_templates(vec![], vec![], vec![graph(Element::O), bad(), bad()])
            .unwrap();
        let mut p = ReactionWriteParams::default();
        p.canonical = true;
        assert!(matches!(
            role_to_string(&r, ReactionRole::Agent, &p),
            Err(ReactionWriteError::Template {
                role: ReactionRole::Agent,
                template: 1,
                ..
            })
        ));
    }
    #[test]
    fn writer_options_are_forwarded_before_role_joining() {
        let q = QueryGraph::from_parts(
            vec![
                QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
                QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::O)),
            ],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        let r = r(vec![q]);
        let mut p = ReactionWriteParams::default();
        p.rooted_at_atom = Some(1);
        let (text, outputs) = write(&r, &p);
        assert_eq!(text, PropertyText::from("([#8].[#6])"));
        assert_eq!(outputs[0].atom_order, vec![AtomId::new(1), AtomId::new(0)]);
    }
    #[test]
    fn empty_selected_role_does_not_visit_other_role_errors() {
        let r = Reaction::from_templates(vec![], vec![bad()], vec![]).unwrap();
        let initialized = r.is_initialized();
        let (text, outputs) = write(&r, &Default::default());
        assert!(text.is_empty());
        assert!(outputs.is_empty());
        assert_eq!(r.is_initialized(), initialized);
    }
}
#[cfg(test)]
mod reaction_writer_complete_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, Element, PropertyValue};
    use cosmolkit_smiles::CxSmilesFields as F;
    fn graph(e: Element, label: &str) -> QueryGraph {
        let mut q = QueryGraph::from_parts(
            vec![cosmolkit_model::QueryAtom::new(
                AtomId::new(0),
                AtomSpec::new(e),
            )],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        if !label.is_empty() {
            q.atom_mut(0).unwrap().set_prop("atomLabel", label).unwrap();
        }
        q
    }
    fn run(r: &Reaction, cx: bool, flags: F) -> Result<SmirksBase, ReactionWriteError> {
        smirks_base(r, &Default::default(), cx, flags)
    }
    #[test]
    fn empty_roles_keep_two_separators_with_or_without_cx() {
        for cx in [false, true] {
            let r = run(&Reaction::new(), cx, F::ALL).unwrap();
            assert_eq!(r.text, PropertyText::from(">>"));
            assert!(r.outputs.is_empty());
        }
    }
    #[test]
    fn roles_emit_r_then_a_then_p_without_initialization() {
        let r = Reaction::from_templates(
            vec![graph(Element::C, "")],
            vec![graph(Element::O, "")],
            vec![graph(Element::F, "")],
        )
        .unwrap();
        assert!(!r.is_initialized());
        let o = run(&r, false, F::ALL).unwrap();
        assert_eq!(o.text, PropertyText::from("[#6]>[#9]>[#8]"));
        assert_eq!(
            o.outputs.iter().map(|x| x.text.clone()).collect::<Vec<_>>(),
            ["[#6]", "[#9]", "[#8]"].map(PropertyText::from)
        );
        assert!(!r.is_initialized());
    }
    #[test]
    fn canonical_role_sort_does_not_sort_cx_template_order() {
        let r = Reaction::from_templates(
            vec![graph(Element::O, "oxygen"), graph(Element::C, "carbon")],
            vec![],
            vec![],
        )
        .unwrap();
        let mut p = ReactionWriteParams::default();
        p.canonical = true;
        let o = smirks_base(&r, &p, true, F::ATOM_LABELS).unwrap();
        assert_eq!(o.text, PropertyText::from("[#6].[#8]>> |$oxygen;carbon$|"));
        assert_eq!(o.outputs[0].text, PropertyText::from("[#8]"));
        assert_eq!(o.outputs[1].text, PropertyText::from("[#6]"));
    }
    #[test]
    fn cx_collection_keeps_reactant_agent_product_order() {
        let r = Reaction::from_templates(
            vec![graph(Element::C, "R")],
            vec![graph(Element::O, "P")],
            vec![graph(Element::F, "A")],
        )
        .unwrap();
        assert_eq!(
            run(&r, true, F::ATOM_LABELS).unwrap().text,
            PropertyText::from("[#6]>[#9]>[#8] |$R;A;P$|")
        );
    }
    #[test]
    fn flags_none_does_not_append_space_or_skip_cx_presence_checks() {
        let r = Reaction::from_templates(vec![graph(Element::C, "R")], vec![], vec![]).unwrap();
        assert_eq!(
            run(&r, true, F::NONE).unwrap().text,
            PropertyText::from("[#6]>>")
        );
        let e = QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap();
        let r = Reaction::from_templates(vec![e], vec![], vec![]).unwrap();
        assert_eq!(
            run(&r, false, F::NONE).unwrap().text,
            PropertyText::from(">>")
        );
        assert!(matches!(
            run(&r, true, F::NONE),
            Err(ReactionWriteError::Cx(
                cosmolkit_search::SmartsWriteError::CxMissingOutputOrder { template: 0 }
            ))
        ));
    }
    #[test]
    fn false_include_cx_does_not_read_molecule_cx_properties() {
        let mut q = graph(Element::C, "");
        q.set_prop("_molLinkNodes", PropertyValue::Double(f64::NAN))
            .unwrap();
        let r = Reaction::from_templates(vec![q], vec![], vec![]).unwrap();
        assert_eq!(
            run(&r, false, F::ALL).unwrap().text,
            PropertyText::from("[#6]>>")
        );
    }
    #[test]
    fn source_counted_atom_labels_survive_role_and_cx_joining() {
        let mut q = graph(Element::C, "");
        q.atom_mut(0)
            .unwrap()
            .set_prop(
                "atomLabel",
                PropertyValue::String(PropertyText::from(vec![b'x', 0, 255])),
            )
            .unwrap();
        let r = Reaction::from_templates(vec![q], vec![], vec![]).unwrap();
        assert_eq!(
            run(&r, true, F::ATOM_LABELS).unwrap().text.as_bytes(),
            b"[#6]>> |$x\0\xff$|"
        );
    }
    #[test]
    fn coordinate_flags_forward_to_single_native_front_frame_writer() {
        let mut q = graph(Element::C, "");
        q.add_conformer_3d(cosmolkit_model::Conformer3D::new(
            9,
            vec![[1., 2., 3.]],
            true,
        ))
        .unwrap();
        q.add_conformer_3d(cosmolkit_model::Conformer3D::new(
            1,
            vec![[4., 5., 6.]],
            true,
        ))
        .unwrap();
        let r = Reaction::from_templates(vec![q], vec![], vec![]).unwrap();
        let o = run(&r, true, F::COORDS).unwrap();
        assert_eq!(o.text, PropertyText::from("[#6]>> |(1,2,3)|"));
        assert_eq!(r.reactant_templates()[0].conformers_3d().len(), 2);
    }
}

pub(crate) fn write_smirks_with_params(
    reaction: &Reaction,
    params: &ReactionWriteParams,
) -> Result<PropertyText, ReactionWriteError> {
    // RDKit❗✔️: std::string ChemicalReactionToRxnSmarts(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                         const SmilesWriteParams &params) {
    // RDKit❗✔️:   return chemicalReactionToRxnToString(rxn, false, params, false);
    // RDKit❗✔️: }
    // Source SMARTS wrapper fixes toSmiles=false and includeCX=false.
    // Caller options are borrowed unchanged; one sole complete writer call.
    // O(1) wrapper overhead, move the source byte result without another copy.
    smirks_base(
        reaction,
        params,
        false,
        cosmolkit_smiles::CxSmilesFields::ALL,
    )
    .map(|output| output.text)
}
#[cfg(test)]
mod reaction_smarts_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, Element, QueryAtom};
    fn g(e: Element) -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(e))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn empty_source_wrapper_returns_both_separators() {
        assert_eq!(
            write_smirks_with_params(&Reaction::new(), &Default::default()).unwrap(),
            PropertyText::from(">>")
        );
    }
    #[test]
    fn wrapper_forwards_canonical_role_sort() {
        let r =
            Reaction::from_templates(vec![g(Element::O), g(Element::C)], vec![], vec![]).unwrap();
        let mut p = ReactionWriteParams::default();
        p.canonical = true;
        assert_eq!(
            write_smirks_with_params(&r, &p).unwrap(),
            PropertyText::from("[#6].[#8]>>")
        );
    }
    #[test]
    fn source_literal_false_skips_empty_template_cx_order_requirement() {
        let q = QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap();
        let r = Reaction::from_templates(vec![q], vec![], vec![]).unwrap();
        let mut p = ReactionWriteParams::default();
        p.include_cx = true;
        assert_eq!(
            write_smirks_with_params(&r, &p).unwrap(),
            PropertyText::from(">>")
        );
    }
    #[test]
    fn source_wrapper_propagates_role_and_index_of_first_writer_error() {
        let r = Reaction::from_templates(vec![], vec![g(Element::C)], vec![]).unwrap();
        let mut p = ReactionWriteParams::default();
        p.rooted_at_atom = Some(9);
        assert!(matches!(
            write_smirks_with_params(&r, &p),
            Err(ReactionWriteError::Template {
                role: ReactionRole::Product,
                template: 0,
                source: cosmolkit_search::SmartsWriteError::RootedAtomOutOfRange { atom: 9 }
            })
        ));
        assert!(!r.is_initialized());
    }
}

pub(crate) fn write_cx_smirks_with_params(
    reaction: &Reaction,
    params: &ReactionWriteParams,
    flags: cosmolkit_smiles::CxSmilesFields,
) -> Result<PropertyText, ReactionWriteError> {
    // RDKit❗✔️: std::string ChemicalReactionToCXRxnSmarts(const ChemicalReaction &rxn,
    // RDKit❗✔️:                                         const SmilesWriteParams &params,
    // RDKit❗✔️:                                         std::uint32_t flags) {
    // RDKit❗✔️:   return chemicalReactionToRxnToString(rxn, false, params, true, flags);
    // RDKit❗✔️: }
    // This distinct source overload fixes false/true, preserving explicit
    // caller flags. No duplicated reaction, molecule writer or CX algorithm.
    // O(1) wrapper overhead; move the canonical counted-byte result.
    smirks_base(reaction, params, true, flags).map(|output| output.text)
}
#[cfg(test)]
mod reaction_cx_wrapper_source_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, Element, QueryAtom};
    use cosmolkit_smiles::CxSmilesFields as F;
    fn r() -> Reaction {
        let mut q = QueryGraph::from_parts(
            vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap();
        q.atom_mut(0).unwrap().set_prop("atomLabel", "R").unwrap();
        Reaction::from_templates(vec![q], vec![], vec![]).unwrap()
    }
    #[test]
    fn literal_true_emits_cx_with_default_non_cx_parameter_flag() {
        let p = ReactionWriteParams::default();
        assert!(!p.include_cx);
        assert_eq!(
            write_cx_smirks_with_params(&r(), &p, F::ATOM_LABELS).unwrap(),
            PropertyText::from("[#6]>> |$R$|")
        );
    }
    #[test]
    fn explicit_none_flags_retain_base_without_trailing_space() {
        assert_eq!(
            write_cx_smirks_with_params(&r(), &Default::default(), F::NONE).unwrap(),
            PropertyText::from("[#6]>>")
        );
    }
    #[test]
    fn empty_reaction_reaches_cx_without_inventing_a_template() {
        assert_eq!(
            write_cx_smirks_with_params(&Reaction::new(), &Default::default(), F::ALL).unwrap(),
            PropertyText::from(">>")
        );
    }
    #[test]
    fn empty_template_failure_is_forwarded_as_structural_cx_error() {
        let e = QueryGraph::from_parts(vec![], vec![], [], vec![], vec![], vec![]).unwrap();
        let rx = Reaction::from_templates(vec![], vec![e], vec![]).unwrap();
        assert!(matches!(
            write_cx_smirks_with_params(&rx, &Default::default(), F::NONE),
            Err(ReactionWriteError::Cx(
                cosmolkit_search::SmartsWriteError::CxMissingOutputOrder { template: 0 }
            ))
        ));
    }
    #[test]
    fn template_writer_error_precedes_requested_cx_and_leaves_lifecycle_unchanged() {
        let rx = r();
        let mut p = ReactionWriteParams::default();
        p.rooted_at_atom = Some(7);
        assert!(matches!(
            write_cx_smirks_with_params(&rx, &p, F::ALL),
            Err(ReactionWriteError::Template {
                role: ReactionRole::Reactant,
                template: 0,
                source: cosmolkit_search::SmartsWriteError::RootedAtomOutOfRange { atom: 7 }
            })
        ));
        assert!(!rx.is_initialized());
        assert_eq!(
            rx.reactant_templates()[0]
                .atom(0)
                .unwrap()
                .prop("atomLabel"),
            Some(&cosmolkit_model::PropertyValue::String(PropertyText::from(
                "R"
            )))
        );
    }
}
