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
            if input.topology.atoms[target_atom]
                .prop("_protected")
                .is_some()
            {
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

pub(crate) fn reactant_matches(
    inputs: &[ReactionInput<'_>],
    reaction: &Reaction,
    max_matches: u32,
    single: Option<usize>,
) -> Result<Option<Vec<AtomMatches>>, ReactionRunError> {
    // BEGIN RDKIT CPP FUNCTION: ReactionRunner.cpp 200-228
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
    // RDKit❗✔️: }  // end of getReactantMatches()
    // END RDKIT CPP FUNCTION
    // One ordered pass; stop after first empty selected reactant as source.
    if inputs.len() != reaction.reactants.len() {
        return Err(ReactionRunError::ReactantArity {
            expected: reaction.reactants.len(),
            actual: inputs.len(),
        });
    }
    if reaction.needs_init {
        return Err(ReactionRunError::NeedsInitialization);
    }
    if let Some(index) = single
        && index >= reaction.reactants.len()
    {
        return Err(ReactionRunError::ReactantTemplateIndex {
            index,
            count: reaction.reactants.len(),
        });
    }
    let mut by_reactant = vec![Vec::new(); inputs.len()];
    for (index, template) in reaction.reactants.iter().enumerate() {
        if single.is_none_or(|selected| selected == index) {
            let matches = matches_to_template(
                inputs[index],
                template,
                max_matches,
                &reaction.match_params,
                index,
                index,
            )?;
            if matches.is_empty() {
                return Ok(None);
            }
            by_reactant[index] = matches;
        }
    }
    Ok(Some(by_reactant))
}

fn recurse_combinations(
    by_reactant: &[AtomMatches],
    per_product: &mut Vec<Vec<Vec<usize>>>,
    level: usize,
    combination: Vec<Vec<usize>>,
    max_products: u32,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: ReactionRunner.cpp 233-264
    // RDKit❗✔️: bool recurseOverReactantCombinations(
    // RDKit❗✔️:     const VectVectMatchVectType &matchesByReactant,
    // RDKit❗✔️:     VectVectMatchVectType &matchesPerProduct, unsigned int level,
    // RDKit❗✔️:     VectMatchVectType combination, unsigned int maxProducts) {
    // RDKit❗✔️:   unsigned int nReactants = matchesByReactant.size();
    // RDKit❗✔️:   URANGE_CHECK(level, nReactants);
    // RDKit❗✔️:   PRECONDITION(combination.size() == nReactants, "bad combination size");
    // RDKit❗✔️:
    // RDKit❗✔️:   if (maxProducts && matchesPerProduct.size() >= maxProducts) {
    // RDKit❗✔️:     return false;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   bool keepGoing = true;
    // RDKit❗✔️:   for (auto reactIt = matchesByReactant[level].begin();
    // RDKit❗✔️:        reactIt != matchesByReactant[level].end(); ++reactIt) {
    // RDKit❗✔️:     VectMatchVectType prod = combination;
    // RDKit❗✔️:     prod[level] = *reactIt;
    // RDKit❗✔️:     if (level == nReactants - 1) {
    // RDKit❗✔️:       // this is the bottom of the recursion:
    // RDKit❗✔️:       if (maxProducts && matchesPerProduct.size() >= maxProducts) {
    // RDKit❗✔️:         keepGoing = false;
    // RDKit❗✔️:         break;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       matchesPerProduct.push_back(prod);
    // RDKit❗✔️:
    // RDKit❗✔️:     } else {
    // RDKit❗✔️:       keepGoing = recurseOverReactantCombinations(
    // RDKit❗✔️:           matchesByReactant, matchesPerProduct, level + 1, prod, maxProducts);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return keepGoing;
    // RDKit❗✔️: }  // end of recurseOverReactantCombinations
    // END RDKIT CPP FUNCTION
    // Source copies the partial combination at each branch and returns the
    // last recursive keepGoing value; do not add an early-break heuristic.
    if max_products != 0 && per_product.len() >= max_products as usize {
        return false;
    }
    let mut keep_going = true;
    for matched in &by_reactant[level] {
        let mut product = combination.clone();
        product[level] = matched.clone();
        if level == by_reactant.len() - 1 {
            if max_products != 0 && per_product.len() >= max_products as usize {
                keep_going = false;
                break;
            }
            per_product.push(product);
        } else {
            keep_going =
                recurse_combinations(by_reactant, per_product, level + 1, product, max_products);
        }
    }
    keep_going
}

pub(crate) fn reactant_combinations(
    by_reactant: &[AtomMatches],
    max_products: u32,
) -> Result<Vec<Vec<Vec<usize>>>, ReactionRunError> {
    // BEGIN RDKIT CPP FUNCTION: ReactionRunner.cpp 288-300
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
    // RDKit❗✔️: }  // end of generateReactantCombinations()
    // END RDKIT CPP FUNCTION
    if by_reactant.is_empty() {
        return Err(ReactionRunError::EmptyCombinationLevels);
    }
    let mut products = Vec::new();
    let combination = vec![Vec::new(); by_reactant.len()];
    if !recurse_combinations(by_reactant, &mut products, 0, combination, max_products) {
        eprintln!("Maximum product count hit {max_products}, stopping reaction early...");
    }
    Ok(products)
}
