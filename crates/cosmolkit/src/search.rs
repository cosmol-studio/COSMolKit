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
        cosmolkit_search::try_get_substruct_matches_with_params(
            &self.detached_search_target(),
            query,
            params,
        )
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
        query.matches_target(&self.detached_search_target())
    }
}
