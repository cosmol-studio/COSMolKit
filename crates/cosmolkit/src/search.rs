//! Canonical public search over detached query values and read-only molecule inputs.

use crate::{Molecule, QueryGraph, SmartsParseError, SmartsParseParams};
use cosmolkit_search::{MatchResult, SubstructMatchError, SubstructMatchParams};

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
) -> Result<String, crate::SmartsWriteError> {
    cosmolkit_search::query_graph_to_smarts(query, params)
}

/// Serialize query annotations through the same SMARTS owner.
pub fn write_cx_smarts(
    query: &QueryGraph,
    params: &crate::SmartsWriteParams,
) -> Result<String, crate::SmartsWriteError> {
    cosmolkit_search::query_graph_to_cx_smarts(query, params)
}

impl Molecule {
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
        assert!(target.valence().is_none());
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
        assert!(molecule.detached_search_target().valence().is_none());
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
