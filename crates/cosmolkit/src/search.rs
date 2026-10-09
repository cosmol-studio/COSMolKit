//! Canonical public search over detached query values and read-only molecule inputs.

use crate::{Molecule, QueryGraph, SmartsParseError, SmartsParseParams};
use cosmolkit_search::{MatchResult, SubstructMatchError, SubstructMatchParams};

/// Find a maximum common query substructure using the source owner's defaults.
///
/// This experimental search borrows at least two molecules without modifying
/// their topology, coordinates, properties or derived caches. The result is a
/// query, not a concrete molecule. Full FMCS parity/performance is not claimed.
pub fn maximum_common_substructure(
    inputs: &[&Molecule],
) -> Result<crate::McsResult, crate::McsError> {
    maximum_common_substructure_with_params(inputs, &crate::McsParameters::default())
}

/// Find an MCS with explicit options. `timeout` is in seconds; an interrupted
/// search returns its best partial result with `completed == false`.
///
/// Options requiring ring or valence state use only valid existing assignments;
/// absent state produces the owner's typed error, never an implicit write-back.
pub fn maximum_common_substructure_with_params(
    inputs: &[&Molecule],
    params: &crate::McsParameters,
) -> Result<crate::McsResult, crate::McsError> {
    let targets: Vec<_> = inputs
        .iter()
        .map(|mol| mol.detached_search_target())
        .collect();
    cosmolkit_search::find_mcs(&targets, params)
}

/// Parse SMARTS with the source owner's unchanged defaults.
pub fn parse_smarts(text: &str) -> Result<QueryGraph, SmartsParseError> {
    parse_smarts_with_params(text, &SmartsParseParams::default())
}

/// Parse SMARTS without creating a live molecule or changing parser options.
pub fn parse_smarts_with_params(
    text: &str,
    params: &SmartsParseParams,
) -> Result<QueryGraph, SmartsParseError> {
    cosmolkit_search::parse_smarts(text, params)
}

/// Registered callable behind the bound QueryGraph factory; uses the sole parser.
pub fn from_smarts(text: &str) -> Result<QueryGraph, SmartsParseError> {
    parse_smarts(text)
}

/// Explicit-parameter factory without adding parser ownership to the model value.
pub fn from_smarts_with_params(
    text: &str,
    params: &SmartsParseParams,
) -> Result<QueryGraph, SmartsParseError> {
    parse_smarts_with_params(text, params)
}

/// Compile an owned detached query for reuse against independent targets.
pub fn compile_query(query: &QueryGraph) -> Result<crate::CompiledQuery, crate::QueryCompileError> {
    cosmolkit_search::compile_query(query)
}

/// Serialize the canonical detached query through its existing owner.
pub fn write_smarts(
    query: &QueryGraph,
    params: &crate::SmartsWriteParams,
) -> Result<crate::PropertyText, crate::SmartsWriteError> {
    cosmolkit_search::query_graph_to_smarts(query, params)
}

/// Serialize query annotations through the same SMARTS owner.
pub fn write_cx_smarts(
    query: &QueryGraph,
    params: &crate::SmartsWriteParams,
) -> Result<crate::PropertyText, crate::SmartsWriteError> {
    cosmolkit_search::query_graph_to_cx_smarts(query, params)
}

impl Molecule {
    /// Serialize concrete rows without changing live molecule state.
    pub fn to_smarts(&self) -> Result<crate::PropertyText, crate::SmartsWriteError> {
        self.to_smarts_with_params(&crate::SmartsWriteParams::default())
    }

    pub fn to_smarts_with_params(
        &self,
        params: &crate::SmartsWriteParams,
    ) -> Result<crate::PropertyText, crate::SmartsWriteError> {
        cosmolkit_search::topology_to_smarts(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            params,
            false,
        )
    }

    /// Serialize concrete rows without changing live molecule state.
    pub fn to_cx_smarts(&self) -> Result<crate::PropertyText, crate::SmartsWriteError> {
        self.to_cx_smarts_with_params(&crate::SmartsWriteParams::default())
    }

    pub fn to_cx_smarts_with_params(
        &self,
        params: &crate::SmartsWriteParams,
    ) -> Result<crate::PropertyText, crate::SmartsWriteError> {
        cosmolkit_search::topology_to_smarts(
            self.topology(),
            self.coordinate_block_runtime(),
            self.properties(),
            params,
            true,
        )
    }

    fn detached_search_target(&self) -> cosmolkit_search::SearchTarget<'_> {
        #[cfg(any(
            feature = "cap-valence",
            feature = "cap-fingerprints",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors",
            feature = "cap-forcefields"
        ))]
        let cache = self.derived_cache_runtime();
        #[cfg(any(
            feature = "cap-valence",
            feature = "cap-fingerprints",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors",
            feature = "cap-forcefields"
        ))]
        let valence = cache.valence_assignment();
        #[cfg(not(any(
            feature = "cap-valence",
            feature = "cap-fingerprints",
            feature = "cap-hydrogens",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-descriptors",
            feature = "cap-forcefields"
        )))]
        let valence = None;
        #[cfg(any(
            feature = "cap-rings",
            feature = "cap-descriptors",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-hydrogens",
            feature = "cap-kekulize",
            feature = "cap-aromaticity",
            feature = "cap-fingerprints"
        ))]
        let rings = self.derived_cache_runtime().valid_ring_info();
        #[cfg(not(any(
            feature = "cap-rings",
            feature = "cap-descriptors",
            feature = "cap-smiles",
            feature = "cap-sanitize",
            feature = "cap-hydrogens",
            feature = "cap-kekulize",
            feature = "cap-aromaticity",
            feature = "cap-fingerprints"
        )))]
        let rings = None;
        cosmolkit_search::SearchTarget::new(
            self.topology(),
            self.coordinate_block_runtime(),
            &self.topology().stereo_groups,
            rings,
            valence,
        )
    }

    /// Return the first source-ordered match, propagating matcher failures.
    pub fn substruct_match(
        &self,
        query: &QueryGraph,
    ) -> Result<Option<MatchResult>, SubstructMatchError> {
        let params = SubstructMatchParams {
            max_matches: 1,
            ..Default::default()
        };
        Ok(self
            .substruct_matches_with_params(query, &params)?
            .into_iter()
            .next())
    }

    /// Return all matches with the owner's exact default options.
    pub fn substruct_matches(
        &self,
        query: &QueryGraph,
    ) -> Result<Vec<MatchResult>, SubstructMatchError> {
        self.substruct_matches_with_params(query, &SubstructMatchParams::default())
    }

    /// Preserve all options, atom/bond mappings, ordering and typed owner errors.
    pub fn substruct_matches_with_params(
        &self,
        query: &QueryGraph,
        params: &SubstructMatchParams,
    ) -> Result<Vec<MatchResult>, SubstructMatchError> {
        let target = self.detached_search_target();
        if let Some(context) = prepared_live_search_context(&target)? {
            cosmolkit_search::try_get_substruct_matches_with_params_and_context(
                &target, query, params, &context,
            )
        } else {
            cosmolkit_search::try_get_substruct_matches_with_params(&target, query, params)
        }
    }

    /// Test for a match without converting an unsupported branch into false.
    pub fn has_substruct_match(&self, query: &QueryGraph) -> Result<bool, SubstructMatchError> {
        Ok(self.substruct_match(query)?.is_some())
    }

    /// Return the first match with the complete configured matcher policy.
    pub fn substruct_match_with_params(
        &self,
        query: &QueryGraph,
        params: &SubstructMatchParams,
    ) -> Result<Option<MatchResult>, SubstructMatchError> {
        // RDKit✔️✔️: SubstructMatchParameters ps = params;
        // RDKit✔️✔️: ps.maxMatches = 1;
        // RDKit✔️✔️: std::vector<MatchVectType> matches;
        // RDKit✔️✔️: pySubstructHelper(mol, query, params, matches);
        // RDKit✔️✔️: MatchVectType match;
        // RDKit✔️✔️: if (matches.size()) {
        // RDKit✔️✔️:   match = matches[0];
        // RDKit✔️✔️: }
        // RDKit✔️✔️: return convertMatches(match);
        // Wrap/substructmethods.h::helpGetSubstructMatch passes original params,
        // not its unused ps copy. Preserve max_matches and callback visits by
        // using the configured owner call, then its first ordered result.
        // No additional molecule or unused parameter copy is needed.
        Ok(self
            .substruct_matches_with_params(query, params)?
            .into_iter()
            .next())
    }

    /// Test for a match with the same policy as the configured all-match query.
    pub fn has_substruct_match_with_params(
        &self,
        query: &QueryGraph,
        params: &SubstructMatchParams,
    ) -> Result<bool, SubstructMatchError> {
        // RDKit✔️✔️: SubstructMatchParameters ps = params;
        // RDKit✔️✔️: ps.maxMatches = 1;
        // RDKit✔️✔️: std::vector<MatchVectType> matches;
        // RDKit✔️✔️: pySubstructHelper(mol, query, params, matches);
        // RDKit✔️✔️: return matches.size() != 0;
        // Same source helper policy; do not truncate callback visits.
        Ok(self.substruct_match_with_params(query, params)?.is_some())
    }

    /// Match a reusable compiled query without recompiling its execution plan.
    pub fn substruct_matches_compiled(
        &self,
        query: &crate::CompiledQuery,
    ) -> Result<Vec<MatchResult>, crate::MatchError> {
        let target = self.detached_search_target();
        if let Some(context) = prepared_live_search_context(&target)? {
            query.matches_prepared_target(&target, &context)
        } else {
            query.matches_target(&target)
        }
    }
}

// The facade owns topology correspondence and reads validity-checked cache
// assignments. Passing its exact prepared context preserves source ring
// quality, memberships and initialized-empty state through direct, recursive
// and compiled matching; no domain algorithm or cache mutation occurs here.
fn prepared_live_search_context<'a>(
    target: &'a cosmolkit_search::SearchTarget<'_>,
) -> Result<Option<cosmolkit_search::QueryMatchContext<'a>>, SubstructMatchError> {
    use cosmolkit_search::SearchTargetAccess;
    match target.ring_info() {
        Some(rings) => cosmolkit_search::build_ring_query_match_context(
            target.topology_block(),
            rings,
            target.valence(),
        )
        .map(Some)
        .map_err(SubstructMatchError::from),
        None => Ok(None),
    }
}

#[cfg(all(test, feature = "cap-smiles", feature = "cap-hydrogens"))]
mod original_smarts_public_regressions {
    use super::*;
    use cosmolkit_search::SearchTargetAccess;
    use std::sync::Arc;

    #[test]
    fn mcs_public_default_query_and_inputs_are_preserved() {
        let a = Molecule::from_smiles("CCO").unwrap();
        let b = Molecule::from_smiles("CCN").unwrap();
        let topology = a.topology_arc_runtime();
        let coordinates = a.coordinates_arc_runtime();
        let cache = a.derived_cache_arc_runtime();
        let result = maximum_common_substructure(&[&a, &b]).unwrap();
        assert_eq!(
            (result.atom_count, result.bond_count, result.completed),
            (2, 1, true)
        );
        assert_eq!(result.smarts, "[#6]-[#6]".into());
        let query = result.query.as_ref().unwrap();
        assert!(a.has_substruct_match(query).unwrap());
        assert!(b.has_substruct_match(query).unwrap());
        assert!(Arc::ptr_eq(&topology, &a.topology_arc_runtime()));
        assert!(Arc::ptr_eq(&coordinates, &a.coordinates_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &a.derived_cache_arc_runtime()));
        assert_eq!(a.to_smiles().unwrap(), "CCO".into());
        assert_eq!(b.to_smiles().unwrap(), "CCN".into());
    }

    #[test]
    fn mcs_public_options_empty_result_and_errors_are_forwarded() {
        let a = Molecule::from_smiles("CCO").unwrap();
        let b = Molecule::from_smiles("CCN").unwrap();
        let params = crate::McsParameters {
            atom_comparator: crate::McsAtomComparator::AtomCompareAny,
            ..Default::default()
        };
        let result = maximum_common_substructure_with_params(&[&a, &b], &params).unwrap();
        assert_eq!((result.atom_count, result.bond_count), (3, 2));
        let c = Molecule::from_smiles("Cl").unwrap();
        let d = Molecule::from_smiles("Br").unwrap();
        let empty = maximum_common_substructure(&[&c, &d]).unwrap();
        assert_eq!(
            (empty.atom_count, empty.bond_count, empty.completed),
            (0, 0, true)
        );
        assert!(empty.query.is_none());
        assert!(empty.smarts.is_empty());
        assert!(matches!(
            maximum_common_substructure(&[&a]),
            Err(crate::McsError::State(
                cosmolkit_search::McsError::TooFewInputs { count: 1 }
            ))
        ));
        let invalid = crate::McsParameters {
            threshold: 1.1,
            ..Default::default()
        };
        assert!(matches!(
            maximum_common_substructure_with_params(&[&a, &b], &invalid),
            Err(crate::McsError::State(
                cosmolkit_search::McsError::ThresholdAboveOne
            ))
        ));
    }

    #[test]
    fn configured_single_and_boolean_matching_preserve_callback_policy_and_inputs() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        let molecule = Molecule::from_smiles("CCC").unwrap();
        let query = crate::parse_smarts("C").unwrap();
        let topology = molecule.topology_arc_runtime();
        let visits = Arc::new(AtomicUsize::new(0));
        let observed = Arc::clone(&visits);
        let params = SubstructMatchParams {
            max_matches: 2,
            extra_final_check: Some(Arc::new(move |_, _| {
                observed.fetch_add(1, Ordering::Relaxed);
                true
            })),
            ..Default::default()
        };
        let first = molecule
            .substruct_match_with_params(&query, &params)
            .unwrap()
            .unwrap();
        assert_eq!(first.atom_mapping, vec![0]);
        // The pinned configured RDKit wrapper passes the original params,
        // rather than its unused maxMatches=1 copy, to the matcher.
        assert_eq!(visits.swap(0, Ordering::Relaxed), 2);
        assert!(
            molecule
                .has_substruct_match_with_params(&query, &params)
                .unwrap()
        );
        assert_eq!(visits.load(Ordering::Relaxed), 2);
        assert_eq!(params.max_matches, 2);
        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert_eq!(molecule.to_smiles().unwrap(), "CCC".into());
    }

    #[test]
    fn original_smarts_prepared_rings_are_borrowed_unchanged() {
        let molecule = Molecule::from_smiles("C12C3C4C1C5C2C3C45").unwrap();
        let topology = molecule.topology_arc_runtime();
        let cache = molecule.derived_cache_arc_runtime();
        let target = molecule.detached_search_target();
        let rings = target.ring_info().unwrap();
        assert_eq!(rings.num_rings(), 6);
        let context = prepared_live_search_context(&target).unwrap().unwrap();
        let query = crate::parse_smarts("[R3]").unwrap();
        let expected = (0..8).map(|i| vec![i]).collect::<Vec<_>>();
        let matches = cosmolkit_search::try_get_substruct_matches_with_params_and_context(
            &target,
            &query,
            &SubstructMatchParams::default(),
            &context,
        )
        .unwrap();
        assert_eq!(
            matches
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(
            molecule
                .substruct_matches(&query)
                .unwrap()
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        let compiled = crate::compile_query(&query).unwrap();
        assert_eq!(
            molecule
                .substruct_matches_compiled(&compiled)
                .unwrap()
                .iter()
                .map(|m| m.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert!(std::ptr::eq(
            rings,
            molecule.detached_search_target().ring_info().unwrap()
        ));
    }

    #[test]
    fn original_smarts_add_hs_preserves_rings_without_valence_cache() {
        let source = Molecule::from_smiles("C12C3C4C1C5C2C3C45").unwrap();
        let molecule = source.with_hydrogens().unwrap();
        let cache = molecule.derived_cache_arc_runtime();
        let topology = molecule.topology_arc_runtime();
        let target = molecule.detached_search_target();
        let actual = target
            .valence()
            .expect("AddHs moves refreshed source scalar rows");
        let expected_valence =
            cosmolkit_core::assign_valence(molecule.topology(), &Default::default()).unwrap();
        assert_eq!(actual, &expected_valence);
        let rings = target.ring_info().unwrap();
        assert_eq!(rings.num_rings(), 6);
        assert_eq!(
            rings.atom_rings(),
            source
                .detached_search_target()
                .ring_info()
                .unwrap()
                .atom_rings()
        );
        assert!(prepared_live_search_context(&target).unwrap().is_some());
        let query = crate::parse_smarts("[C;H1]-[C;R3]").unwrap();
        let expected = vec![
            vec![0, 1],
            vec![0, 3],
            vec![0, 5],
            vec![1, 2],
            vec![1, 6],
            vec![2, 3],
            vec![2, 7],
            vec![3, 4],
            vec![4, 5],
            vec![4, 7],
            vec![5, 6],
            vec![6, 7],
        ];
        assert_eq!(
            molecule
                .substruct_matches(&query)
                .unwrap()
                .iter()
                .map(|r| r.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        let compiled = crate::compile_query(&query).unwrap();
        assert_eq!(
            molecule
                .substruct_matches_compiled(&compiled)
                .unwrap()
                .iter()
                .map(|r| r.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert!(Arc::ptr_eq(&cache, &molecule.derived_cache_arc_runtime()));
        assert!(Arc::ptr_eq(&topology, &molecule.topology_arc_runtime()));
        assert!(std::ptr::eq(
            rings,
            molecule.detached_search_target().ring_info().unwrap()
        ));
        assert_eq!(
            molecule.detached_search_target().valence(),
            Some(&expected_valence)
        );
        // Independent valence preparation must also retain sparse ring rows.
        let assigned = molecule.with_assigned_valence().unwrap();
        assert!(assigned.detached_search_target().valence().is_some());
        assert_eq!(
            assigned
                .substruct_matches(&query)
                .unwrap()
                .iter()
                .map(|r| r.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
        assert_eq!(
            assigned
                .substruct_matches_compiled(&compiled)
                .unwrap()
                .iter()
                .map(|r| r.atom_mapping.clone())
                .collect::<Vec<_>>(),
            expected
        );
    }

    #[test]
    fn original_smarts_all_142_ordered_public_observations() {
        let fixture: serde_json::Value = serde_json::from_str(include_str!(
            "../../../testdata/regression/smarts_user_original104/observations142.json"
        ))
        .unwrap();
        assert_eq!(fixture["groups"].as_array().unwrap().len(), 104);
        let observations = fixture["observations"].as_array().unwrap();
        assert_eq!(observations.len(), 142);
        for o in observations {
            let label = o["id"].as_str().unwrap();
            let mut molecule = Molecule::from_smiles(o["smiles"].as_str().unwrap()).unwrap();
            if o["add_hydrogens"].as_bool().unwrap() {
                molecule = molecule.with_hydrogens().unwrap();
            }
            let query = crate::parse_smarts(o["query"].as_str().unwrap()).unwrap();
            let params = SubstructMatchParams {
                max_matches: o["max_matches"].as_u64().unwrap() as usize,
                uniquify: o["uniquify"].as_bool().unwrap(),
                use_chirality: o["use_chirality"].as_bool().unwrap(),
                ..Default::default()
            };
            let expected: Vec<Vec<usize>> =
                serde_json::from_value(o["rdkit_maps"].clone()).unwrap();
            let actual = molecule
                .substruct_matches_with_params(&query, &params)
                .unwrap();
            assert_eq!(
                actual
                    .iter()
                    .map(|m| m.atom_mapping.clone())
                    .collect::<Vec<_>>(),
                expected,
                "{label}"
            );
            if !params.use_chirality {
                let compiled = crate::compile_query(&query).unwrap();
                let actual = molecule.substruct_matches_compiled(&compiled).unwrap();
                assert_eq!(
                    actual
                        .iter()
                        .map(|m| m.atom_mapping.clone())
                        .collect::<Vec<_>>(),
                    expected,
                    "compiled {label}"
                );
            }
        }
    }
}
