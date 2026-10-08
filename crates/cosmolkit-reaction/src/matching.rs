use crate::{Reaction, ReactionInput, ReactionRunError};
use cosmolkit_model::QueryGraph;
use cosmolkit_search::{
    SearchTarget, SubstructMatchParams, build_query_match_context,
    try_get_substruct_atom_matches_with_params_and_context,
};

type AtomMatches = Vec<Vec<usize>>;

pub(crate) fn matches_to_template(
    input: ReactionInput<'_>,
    template: &QueryGraph,
    max_matches: u32,
    params: &SubstructMatchParams,
    reactant_index: usize,
    template_index: usize,
) -> Result<AtomMatches, ReactionRunError> {
    // BEGIN RDKIT CPP FUNCTION: ReactionRunner.cpp 156-198
    // RDKit❗✔️: VectMatchVectType getReactantMatchesToTemplate(
    // RDKit❗✔️:     const ROMol &reactant, const ROMol &templ, unsigned int maxMatches,
    // RDKit❗✔️:     const SubstructMatchParameters &ssparams) {
    // RDKit❗✔️:   // NOTE that we are *not* uniquifying the results.
    // RDKit❗✔️:   //   This is because we need multiple matches in reactions. For example,
    // RDKit❗✔️:   //   The ring-closure coded as:
    // RDKit❗✔️:   //     [C:1]=[C:2] + [C:3]=[C:4][C:5]=[C:6] ->
    // RDKit❗✔️:   //     [C:1]1[C:2][C:3][C:4]=[C:5][C:6]1
    // RDKit❗✔️:   //   should give 4 products here:
    // RDKit❗✔️:   //     [Cl]C=C + [Br]C=CC=C ->
    // RDKit❗✔️:   //       [Cl]C1C([Br])C=CCC1
    // RDKit❗✔️:   //       [Cl]C1CC(Br)C=CC1
    // RDKit❗✔️:   //       C1C([Br])C=CCC1[Cl]
    // RDKit❗✔️:   //       C1CC([Br])C=CC1[Cl]
    // RDKit❗✔️:   //   Yes, in this case there are only 2 unique products, but that's
    // RDKit❗✔️:   //   a factor of the reactants' symmetry.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //   There's no particularly straightforward way of solving this problem
    // RDKit❗✔️:   //   of recognizing cases where we should give all matches and cases where we
    // RDKit❗✔️:   //   shouldn't; it's safer to just produce everything and let the caller deal
    // RDKit❗✔️:   //   with uniquifying their results.
    // RDKit❗✔️:   VectMatchVectType res;
    // RDKit❗✔️:
    // RDKit❗✔️:   SubstructMatchParameters ssps = ssparams;
    // RDKit❗✔️:   ssps.uniquify = false;
    // RDKit❗✔️:   ssps.maxMatches = maxMatches;
    // RDKit❗✔️:   auto matchesHere = SubstructMatch(reactant, templ, ssps);
    // RDKit❗✔️:   res.reserve(matchesHere.size());
    // RDKit❗✔️:   for (const auto &match : matchesHere) {
    // RDKit❗✔️:     bool keep = true;
    // RDKit❗✔️:     for (const auto &pr : match) {
    // RDKit❗✔️:       if (reactant.getAtomWithIdx(pr.second)->hasProp(
    // RDKit❗✔️:               common_properties::_protected)) {
    // RDKit❗✔️:         keep = false;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:     if (keep) {
    // RDKit❗✔️:       res.push_back(std::move(match));
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION
    // Same SOURCE SEARCH owner, borrowed facts, atom-only projection: no query
    // or cache copies, no unused bond mappings. Source match order is retained.
    let target = SearchTarget::new(
        input.topology,
        input.coordinates,
        &input.topology.stereo_groups,
        input.rings,
        input.valence,
    );
    let context = build_query_match_context(&target);
    let mut source_params = params.clone();
    source_params.uniquify = false;
    source_params.max_matches = max_matches as usize;
    let matches = try_get_substruct_atom_matches_with_params_and_context(
        &target,
        template,
        &source_params,
        &context,
    )
    .map_err(|source| ReactionRunError::Matching {
        reactant: reactant_index,
        template: template_index,
        source,
    })?;
    let mut retained = Vec::with_capacity(matches.len());
    for matched in matches {
        let mut keep = true;
        for &target_atom in &matched {
            // RDKit❗✔️: const Atom *ROMol::getAtomWithIdx(unsigned int idx) const {
            // RDKit❗✔️:   URANGE_CHECK(idx, getNumAtoms());
            // RDKit❗✔️:
            // RDKit❗✔️:   auto vd = boost::vertex(idx, d_graph);
            // RDKit❗✔️:   const auto res = d_graph[vd];
            // RDKit❗✔️:
            // RDKit❗✔️:   POSTCONDITION(res, "");
            // RDKit❗✔️:   return res;
            // RDKit❗✔️: }
            let atom =
                input
                    .topology
                    .atoms
                    .get(target_atom)
                    .ok_or(ReactionRunError::MatchAtomIndex {
                        reactant: reactant_index,
                        template: template_index,
                        atom: target_atom,
                        atom_count: input.topology.atoms.len(),
                    })?;
            // RDKit❗🔝: bool hasProp(const std::string_view key) const { return d_props.hasVal(key); }
            // RDKit❗🔝:   bool hasVal(const std::string_view what) const {
            // RDKit❗🔝:     for (const auto &data : _data) {
            // RDKit❗🔝:       if (data.key == what) {
            // RDKit❗🔝:         return true;
            // RDKit❗🔝:       }
            // RDKit❗🔝:     }
            // RDKit❗🔝:     return false;
            // RDKit❗🔝:   }
            // Presence only: raw wrong-tag, zero and false values protect too.
            // Canonical BTreeMap lookup avoids the native linear Dict scan.
            if atom.prop("_protected").is_some() {
                keep = false;
                break;
            }
        }
        if keep {
            retained.push(matched);
        }
    }
    Ok(retained)
}

fn get_reactant_matches_source(
    inputs: &[ReactionInput<'_>],
    reaction: &Reaction,
    by_reactant: &mut Vec<AtomMatches>,
    max_matches: u32,
    selected: u32,
) -> Result<bool, ReactionRunError> {
    // RDKit❗✔️: bool getReactantMatches(const MOL_SPTR_VECT &reactants,
    // RDKit❗✔️:                         const ChemicalReaction &rxn,
    // RDKit❗✔️:                         VectVectMatchVectType &matchesByReactant,
    // RDKit❗✔️:                         unsigned int maxMatches,
    // RDKit❗✔️:                         unsigned int matchSingleReactant = MatchAll) {
    // RDKit❗✔️:   PRECONDITION(reactants.size() == rxn.getNumReactantTemplates(),
    // RDKit❗✔️:                "reactant size mismatch");
    // RDKit❗✔️:
    // RDKit❗✔️:   matchesByReactant.clear();
    // RDKit❗✔️:   matchesByReactant.resize(reactants.size());
    // RDKit❗✔️:
    // RDKit❗✔️:   bool res = true;
    // RDKit❗✔️:   unsigned int i = 0;
    // RDKit❗✔️:   for (auto iter = rxn.beginReactantTemplates();
    // RDKit❗✔️:        iter != rxn.endReactantTemplates(); ++iter, i++) {
    // RDKit❗✔️:     if (matchSingleReactant == MatchAll || matchSingleReactant == i) {
    // RDKit❗✔️:       auto matches =
    // RDKit❗✔️:           getReactantMatchesToTemplate(*reactants[i].get(), *iter->get(),
    // RDKit❗✔️:                                        maxMatches, rxn.getSubstructParams());
    // RDKit❗✔️:       if (matches.empty()) {
    // RDKit❗✔️:         // no point continuing if we don't match one of the reactants:
    // RDKit❗✔️:         res = false;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       matchesByReactant[i] = std::move(matches);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // Native MatchAll is UINT_MAX. Initialization/index checks belong to the
    // reached runner, not this source helper; unselected slots stay empty.
    if inputs.len() != reaction.num_reactant_templates() {
        return Err(ReactionRunError::ReactantArity {
            expected: reaction.num_reactant_templates(),
            actual: inputs.len(),
        });
    }
    by_reactant.clear();
    by_reactant.resize_with(inputs.len(), Vec::new);
    let mut source_index = 0u32;
    for template in &reaction.reactants {
        let index = source_index as usize;
        if selected == u32::MAX || selected == source_index {
            let matches = matches_to_template(
                inputs[index],
                template,
                max_matches,
                &reaction.match_params,
                index,
                index,
            )?;
            if matches.is_empty() {
                return Ok(false);
            }
            by_reactant[index] = matches;
        }
        source_index = source_index.wrapping_add(1);
    }
    Ok(true)
}

pub(crate) fn reactant_matches(
    inputs: &[ReactionInput<'_>],
    reaction: &Reaction,
    max_matches: u32,
    single: Option<usize>,
) -> Result<Option<Vec<AtomMatches>>, ReactionRunError> {
    let mut by_reactant = Vec::new();
    let selected = single.map_or(u32::MAX, |index| index as u32);
    if !get_reactant_matches_source(inputs, reaction, &mut by_reactant, max_matches, selected)? {
        return Ok(None);
    }
    Ok(Some(by_reactant))
}

fn recurse_combinations(
    by_reactant: &[AtomMatches],
    per_product: &mut Vec<Vec<Vec<usize>>>,
    level: u32,
    combination: Vec<Vec<usize>>,
    max_products: u32,
) -> Result<bool, ReactionRunError> {
    // RDKit❗🔝: bool recurseOverReactantCombinations(
    // RDKit❗🔝:     const VectVectMatchVectType &matchesByReactant,
    // RDKit❗🔝:     VectVectMatchVectType &matchesPerProduct, unsigned int level,
    // RDKit❗🔝:     VectMatchVectType combination, unsigned int maxProducts) {
    // RDKit❗🔝:   unsigned int nReactants = matchesByReactant.size();
    // RDKit❗🔝:   URANGE_CHECK(level, nReactants);
    // RDKit❗🔝:   PRECONDITION(combination.size() == nReactants, "bad combination size");
    // RDKit❗🔝:
    // RDKit❗🔝:   if (maxProducts && matchesPerProduct.size() >= maxProducts) {
    // RDKit❗🔝:     return false;
    // RDKit❗🔝:   }
    // RDKit❗🔝:
    // RDKit❗🔝:   bool keepGoing = true;
    // RDKit❗🔝:   for (auto reactIt = matchesByReactant[level].begin();
    // RDKit❗🔝:        reactIt != matchesByReactant[level].end(); ++reactIt) {
    // RDKit❗🔝:     VectMatchVectType prod = combination;
    // RDKit❗🔝:     prod[level] = *reactIt;
    // RDKit❗🔝:     if (level == nReactants - 1) {
    // RDKit❗🔝:       // this is the bottom of the recursion:
    // RDKit❗🔝:       if (maxProducts && matchesPerProduct.size() >= maxProducts) {
    // RDKit❗🔝:         keepGoing = false;
    // RDKit❗🔝:         break;
    // RDKit❗🔝:       }
    // RDKit❗🔝:       matchesPerProduct.push_back(prod);
    // RDKit❗🔝:
    // RDKit❗🔝:     } else {
    // RDKit❗🔝:       keepGoing = recurseOverReactantCombinations(
    // RDKit❗🔝:           matchesByReactant, matchesPerProduct, level + 1, prod, maxProducts);
    // RDKit❗🔝:     }
    // RDKit❗🔝:   }
    // RDKit❗🔝:   return keepGoing;
    // RDKit❗🔝: }
    // Preserve every branch/duplicate and last recursive return; do not add
    // an early break after recursion. Native prod lvalues copy into the next
    // by-value parameter/output; moving owned vectors here removes those
    // copies without changing any input/output alias or encounter order.
    let count = by_reactant.len() as u32;
    if level >= count {
        return Err(ReactionRunError::CombinationLevel { level, count });
    }
    if combination.len() != count as usize {
        return Err(ReactionRunError::CombinationSize {
            expected: count,
            actual: combination.len(),
        });
    }
    if max_products != 0 && per_product.len() >= max_products as usize {
        return Ok(false);
    }
    let index = level as usize;
    let mut keep_going = true;
    for matched in &by_reactant[index] {
        let mut product = combination.clone();
        product[index] = matched.clone();
        if level == count - 1 {
            if max_products != 0 && per_product.len() >= max_products as usize {
                keep_going = false;
                break;
            }
            per_product.push(product);
        } else {
            keep_going = recurse_combinations(
                by_reactant,
                per_product,
                level.wrapping_add(1),
                product,
                max_products,
            )?;
        }
    }
    Ok(keep_going)
}

fn generate_reactant_combinations_source(
    by_reactant: &[AtomMatches],
    products: &mut Vec<Vec<Vec<usize>>>,
    max_products: u32,
) -> Result<(), ReactionRunError> {
    // RDKit❗✔️: void generateReactantCombinations(
    // RDKit❗✔️:     const VectVectMatchVectType &matchesByReactant,
    // RDKit❗✔️:     VectVectMatchVectType &matchesPerProduct, unsigned int maxProducts) {
    // RDKit❗✔️:   matchesPerProduct.clear();
    // RDKit❗✔️:   VectMatchVectType tmp;
    // RDKit❗✔️:   tmp.clear();
    // RDKit❗✔️:   tmp.resize(matchesByReactant.size());
    // RDKit❗✔️:   if (!recurseOverReactantCombinations(matchesByReactant, matchesPerProduct, 0,
    // RDKit❗✔️:                                        tmp, maxProducts)) {
    // RDKit❗✔️:     BOOST_LOG(rdWarningLog) << "Maximum product count hit " << maxProducts
    // RDKit❗✔️:                             << ", stopping reaction early...\n";
    // RDKit❗✔️:   }
    // RDKit❗✔️: }
    // Clear the real caller output before the reached recursion preconditions,
    // including zero levels. Do not substitute a successful empty product set.
    products.clear();
    let combination = vec![Vec::new(); by_reactant.len()];
    if !recurse_combinations(by_reactant, products, 0, combination, max_products)? {
        eprintln!("Maximum product count hit {max_products}, stopping reaction early...");
    }
    Ok(())
}

pub(crate) fn reactant_combinations(
    by_reactant: &[AtomMatches],
    max_products: u32,
) -> Result<Vec<Vec<Vec<usize>>>, ReactionRunError> {
    let mut products = Vec::new();
    generate_reactant_combinations_source(by_reactant, &mut products, max_products)?;
    Ok(products)
}

#[cfg(test)]
mod complete_matches_to_template_source_tests {
    use super::*;
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, MoleculeProperties,
        PropertyValue, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    fn target(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            elements
                .iter()
                .enumerate()
                .map(|(i, &e)| Atom::from_spec(AtomId::new(i), AtomSpec::new(e)))
                .collect(),
            edges
                .iter()
                .enumerate()
                .map(|(i, &(a, b))| {
                    Bond::from_spec(
                        BondId::new(i),
                        BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn template(text: &str) -> QueryGraph {
        crate::parse_smirks(&format!("{text}>>C"))
            .unwrap()
            .reactants
            .remove(0)
    }
    fn run(
        topology: &TopologyBlock,
        template: &QueryGraph,
        max: u32,
        params: &SubstructMatchParams,
    ) -> Result<AtomMatches, ReactionRunError> {
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        matches_to_template(
            ReactionInput {
                topology,
                coordinates: &coordinates,
                properties: &properties,
                rings: None,
                valence: None,
            },
            template,
            max,
            params,
            7,
            9,
        )
    }

    #[test]
    fn symmetry_is_not_uniquified_and_native_max_overrides_caller_options() {
        let topology = target(&[Element::C, Element::C], &[(0, 1)]);
        let query = template("CC");
        let params = SubstructMatchParams {
            uniquify: true,
            max_matches: 1,
            ..Default::default()
        };
        assert_eq!(
            run(&topology, &query, 0, &params).unwrap(),
            [vec![0, 1], vec![1, 0]]
        );
        assert_eq!(run(&topology, &query, 1, &params).unwrap(), [vec![0, 1]]);
        assert!(params.uniquify);
        assert_eq!(params.max_matches, 1);
    }

    #[test]
    fn protected_is_exact_property_presence_for_false_zero_and_wrong_tag_values() {
        for value in [
            PropertyValue::Bool(false),
            PropertyValue::Int(0),
            PropertyValue::String("raw".into()),
        ] {
            let mut topology = target(&[Element::C, Element::C, Element::C], &[]);
            topology.atoms[1].set_prop("_protected", value).unwrap();
            assert_eq!(
                run(&topology, &template("C"), 0, &Default::default()).unwrap(),
                [vec![0], vec![2]]
            );
            assert!(topology.atoms[1].prop("_protected").is_some());
        }
    }

    #[test]
    fn native_match_cap_precedes_protection_filter_without_refilling() {
        let mut topology = target(&[Element::C, Element::C], &[]);
        topology.atoms[0]
            .set_prop("_protected", PropertyValue::Bool(false))
            .unwrap();
        assert!(
            run(&topology, &template("C"), 1, &Default::default())
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            run(&topology, &template("C"), 2, &Default::default()).unwrap(),
            [vec![1]]
        );
    }

    #[test]
    fn any_protected_atom_rejects_the_whole_mapping_without_reordering_survivors() {
        let mut topology = target(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
        topology.atoms[0]
            .set_prop("_protected", PropertyValue::Int(0))
            .unwrap();
        let matches = run(&topology, &template("CC"), 0, &Default::default()).unwrap();
        assert_eq!(matches, [vec![1, 2], vec![2, 1]]);
        assert!(
            run(&topology, &template("N"), 0, &Default::default())
                .unwrap()
                .is_empty()
        );
    }

    #[test]
    fn remaining_match_options_and_final_callback_are_forwarded_to_the_single_search_owner() {
        let topology = target(&[Element::C, Element::C, Element::C], &[]);
        let calls = Arc::new(AtomicUsize::new(0));
        let observed = calls.clone();
        let params = SubstructMatchParams {
            extra_final_check: Some(Arc::new(move |_, m| {
                observed.fetch_add(1, Ordering::SeqCst);
                m == [2]
            })),
            ..Default::default()
        };
        assert_eq!(
            run(&topology, &template("C"), 0, &params).unwrap(),
            [vec![2]]
        );
        assert_eq!(calls.load(Ordering::SeqCst), 3);
        assert!(params.extra_final_check.is_some());
    }
    mod complete_get_reactant_matches_source_tests {
        use super::*;

        fn input<'a>(
            topology: &'a TopologyBlock,
            coordinates: &'a CoordinateBlock,
            properties: &'a MoleculeProperties,
        ) -> ReactionInput<'a> {
            ReactionInput {
                topology,
                coordinates,
                properties,
                rings: None,
                valence: None,
            }
        }
        fn reaction(templates: &[&str]) -> Reaction {
            let mut reaction = Reaction::new();
            reaction.reactants = templates.iter().map(|text| template(text)).collect();
            reaction
        }

        #[test]
        fn arity_precondition_precedes_output_clear() {
            let reaction = reaction(&["C"]);
            let mut output = vec![vec![vec![42]]];
            assert!(matches!(
                get_reactant_matches_source(&[], &reaction, &mut output, 0, u32::MAX),
                Err(ReactionRunError::ReactantArity {
                    expected: 1,
                    actual: 0
                })
            ));
            assert_eq!(output, [vec![vec![42]]]);
        }

        #[test]
        fn output_resets_and_matches_without_extra_initialization_gate() {
            let reaction = reaction(&["C", "N"]);
            assert!(reaction.needs_init);
            let c = target(&[Element::C, Element::C], &[]);
            let n = target(&[Element::N], &[]);
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let inputs = [
                input(&c, &coordinates, &properties),
                input(&n, &coordinates, &properties),
            ];
            let mut output = vec![vec![vec![42]]; 4];
            assert!(
                get_reactant_matches_source(&inputs, &reaction, &mut output, 0, u32::MAX).unwrap()
            );
            assert_eq!(output, [vec![vec![0], vec![1]], vec![vec![0]]]);
            assert!(reaction.needs_init);
        }

        #[test]
        fn single_selection_skips_all_other_templates_and_out_of_range_is_no_selection() {
            let reaction = reaction(&["N", "C"]);
            let c = target(&[Element::C], &[]);
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let inputs = [input(&c, &coordinates, &properties); 2];
            let mut output = vec![vec![vec![42]]];
            assert!(get_reactant_matches_source(&inputs, &reaction, &mut output, 0, 1).unwrap());
            assert_eq!(output, [vec![], vec![vec![0]]]);
            assert!(get_reactant_matches_source(&inputs, &reaction, &mut output, 0, 9).unwrap());
            assert_eq!(output.len(), 2);
            assert!(output.iter().all(Vec::is_empty));
            assert!(
                !get_reactant_matches_source(&inputs, &reaction, &mut output, 0, u32::MAX).unwrap()
            );
        }

        #[test]
        fn first_empty_selected_match_keeps_prior_output_and_does_not_evaluate_later_template() {
            let mut reaction = reaction(&["C", "N", "C"]);
            let calls = Arc::new(AtomicUsize::new(0));
            let observed = calls.clone();
            reaction.match_params.extra_final_check = Some(Arc::new(move |_, _| {
                observed.fetch_add(1, Ordering::SeqCst);
                true
            }));
            let c = target(&[Element::C, Element::C], &[]);
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let inputs = [input(&c, &coordinates, &properties); 3];
            let mut output = vec![vec![vec![42]]];
            assert!(
                !get_reactant_matches_source(&inputs, &reaction, &mut output, 0, u32::MAX).unwrap()
            );
            assert_eq!(output, [vec![vec![0], vec![1]], vec![], vec![]]);
            assert_eq!(calls.load(Ordering::SeqCst), 2);
        }

        #[test]
        fn reached_search_error_keeps_prior_slots_and_exact_reactant_context() {
            let reaction = reaction(&["C", "[C;H1]"]);
            let c = target(&[Element::C], &[]);
            let coordinates = CoordinateBlock::default();
            let properties = MoleculeProperties::default();
            let inputs = [input(&c, &coordinates, &properties); 2];
            let mut output = vec![vec![vec![42]]];
            assert!(matches!(
                get_reactant_matches_source(&inputs, &reaction, &mut output, 0, u32::MAX),
                Err(ReactionRunError::Matching {
                    reactant: 1,
                    template: 1,
                    ..
                })
            ));
            assert_eq!(output, [vec![vec![0]], vec![]]);
        }
    }
}

#[cfg(test)]
mod complete_recurse_combinations_source_tests {
    use super::*;

    #[test]
    fn source_preconditions_precede_existing_output_cap_and_keep_output_unchanged() {
        let mut output = vec![vec![vec![42]]];
        assert!(matches!(
            recurse_combinations(&[], &mut output, 0, vec![], 1),
            Err(ReactionRunError::CombinationLevel { level: 0, count: 0 })
        ));
        assert!(matches!(
            recurse_combinations(&[vec![]], &mut output, 1, vec![], 1),
            Err(ReactionRunError::CombinationLevel { level: 1, count: 1 })
        ));
        assert!(matches!(
            recurse_combinations(&[vec![]], &mut output, 0, vec![], 1),
            Err(ReactionRunError::CombinationSize {
                expected: 1,
                actual: 0
            })
        ));
        assert_eq!(output, [vec![vec![42]]]);
    }

    #[test]
    fn encounter_order_duplicates_and_empty_mappings_survive_cartesian_recursion() {
        let levels = vec![vec![vec![1], vec![1]], vec![vec![], vec![3, 4]]];
        let mut output = vec![];
        assert!(recurse_combinations(&levels, &mut output, 0, vec![vec![]; 2], 0).unwrap());
        assert_eq!(
            output,
            [
                vec![vec![1], vec![]],
                vec![vec![1], vec![3, 4]],
                vec![vec![1], vec![]],
                vec![vec![1], vec![3, 4]]
            ]
        );
        assert_eq!(levels[0], [vec![1], vec![1]]);
    }

    #[test]
    fn source_keep_going_reports_attempted_extra_candidate_not_merely_reaching_cap() {
        let levels = vec![vec![vec![1], vec![2]], vec![vec![3], vec![4]]];
        for (cap, keep, expected_len) in [(1, false, 1), (3, false, 3), (4, true, 4), (0, true, 4)]
        {
            let mut output = vec![];
            assert_eq!(
                recurse_combinations(&levels, &mut output, 0, vec![vec![]; 2], cap).unwrap(),
                keep
            );
            assert_eq!(output.len(), expected_len);
            assert_eq!(output[0], [vec![1], vec![3]]);
        }
        let mut output = vec![vec![vec![42]]];
        assert!(!recurse_combinations(&[vec![vec![1]]], &mut output, 0, vec![vec![]], 1).unwrap());
        assert_eq!(output, [vec![vec![42]]]);
        assert!(recurse_combinations(&[vec![vec![1]]], &mut output, 0, vec![vec![]], 0).unwrap());
        assert_eq!(output, [vec![vec![42]], vec![vec![1]]]);
    }

    #[test]
    fn empty_level_produces_no_product_and_returns_true_without_clearing_prior_output() {
        for levels in [vec![vec![], vec![vec![1]]], vec![vec![vec![1]], vec![]]] {
            let mut output = vec![vec![vec![42]]];
            assert!(recurse_combinations(&levels, &mut output, 0, vec![vec![]; 2], 0).unwrap());
            assert_eq!(output, [vec![vec![42]]]);
        }
    }

    #[test]
    fn nonzero_start_level_preserves_supplied_combination_prefix() {
        let levels = vec![vec![], vec![vec![3], vec![4]]];
        let mut output = vec![];
        assert!(
            recurse_combinations(&levels, &mut output, 1, vec![vec![42], vec![99]], 0).unwrap()
        );
        assert_eq!(output, [vec![vec![42], vec![3]], vec![vec![42], vec![4]]]);
    }
}

#[cfg(test)]
mod complete_generate_reactant_combinations_source_tests {
    use super::*;

    #[test]
    fn native_clear_precedes_zero_level_recursion_error_and_preserves_output_allocation() {
        let mut output = Vec::with_capacity(32);
        output.push(vec![vec![42]]);
        let storage = output.as_ptr();
        assert!(matches!(
            generate_reactant_combinations_source(&[], &mut output, 1),
            Err(ReactionRunError::CombinationLevel { level: 0, count: 0 })
        ));
        assert!(output.is_empty());
        assert_eq!(output.capacity(), 32);
        assert_eq!(output.as_ptr(), storage);
    }

    #[test]
    fn native_output_reset_retains_order_and_duplicates_in_real_caller_storage() {
        let levels = vec![vec![vec![1], vec![1]], vec![vec![2], vec![3]]];
        let mut output = Vec::with_capacity(32);
        output.push(vec![vec![42]]);
        let storage = output.as_ptr();
        generate_reactant_combinations_source(&levels, &mut output, 0).unwrap();
        assert_eq!(
            output,
            [
                vec![vec![1], vec![2]],
                vec![vec![1], vec![3]],
                vec![vec![1], vec![2]],
                vec![vec![1], vec![3]]
            ]
        );
        assert_eq!(output.as_ptr(), storage);
        assert_eq!(levels[0], [vec![1], vec![1]]);
    }

    #[test]
    fn cap_truncation_returns_normally_with_exact_prefix_and_no_deduplication() {
        let levels = vec![vec![vec![1], vec![2]], vec![vec![3], vec![4]]];
        let all = reactant_combinations(&levels, 0).unwrap();
        for cap in [1, 2, 3, 4, 5] {
            let actual = reactant_combinations(&levels, cap).unwrap();
            assert_eq!(actual, all[..(cap as usize).min(all.len())]);
        }
    }

    #[test]
    fn empty_match_level_and_empty_mapping_value_have_distinct_native_results() {
        let mut output = vec![vec![vec![42]]];
        generate_reactant_combinations_source(&[vec![]], &mut output, 0).unwrap();
        assert!(output.is_empty());
        generate_reactant_combinations_source(&[vec![vec![]]], &mut output, 0).unwrap();
        assert_eq!(output, [vec![Vec::<usize>::new()]]);
        assert!(matches!(
            reactant_combinations(&[], 0),
            Err(ReactionRunError::CombinationLevel { level: 0, count: 0 })
        ));
    }
}
