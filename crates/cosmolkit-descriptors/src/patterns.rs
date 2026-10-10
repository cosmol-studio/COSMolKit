//! Retained fixed SMARTS patterns for descriptor count functions.
//!
//! Source: RDKit `Lipinski.cpp` anonymous-namespace `ss_matcher` plus the
//! `pattern_flyweight` (lines 29-74). Each fixed pattern string is parsed
//! exactly once per process and the compiled query is retained and shared by
//! every later call; descriptor count units never re-parse per call.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use cosmolkit_search::{
    MatchResult, QueryGraph, QueryMatchContext, SearchTarget, SmartsParseParams,
    SubstructMatchParams, build_ring_query_match_context, build_topology_query_match_context,
    build_valence_query_match_context, parse_smarts,
    try_get_substruct_matches_with_params_and_context,
};

use crate::{DescriptorError, DescriptorInput, DescriptorResult, DescriptorSearchCause};

/// Successfully compiled fixed patterns, keyed by the pattern text.
///
/// This is the `no_tracking` flyweight: entries are inserted once and never
/// evicted. Failed compiles are not retained; the source aborts on a failed
/// pattern compile (`POSTCONDITION`), so it defines no retention behavior for
/// failures, and COSMolKit surfaces the typed compile error on each attempt.
type PatternStore = HashMap<&'static str, Arc<QueryGraph>>;

fn pattern_store() -> &'static Mutex<PatternStore> {
    static PATTERNS: OnceLock<Mutex<PatternStore>> = OnceLock::new();
    PATTERNS.get_or_init(|| Mutex::new(HashMap::new()))
}

#[cfg(test)]
thread_local! {
    /// Per-thread count of pattern compilations, mirroring the search
    /// crate's cold-context counter. Only observable in tests.
    pub(crate) static PATTERN_COMPILES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };

    /// Per-thread count of prepared query-context BUILDS performed by this
    /// crate. Only observable in tests; used to prove that a whole
    /// multi-pattern evaluation (e.g. StrictLinkages) builds exactly ONE
    /// context and reuses it for every pattern through the `_with_context`
    /// path instead of preparing per pattern.
    pub(crate) static CONTEXT_BUILDS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

/// Build ONE prepared shared query context from the input's final rows.
///
/// Every descriptor-side context construction goes through this single
/// entry so the test-only [`CONTEXT_BUILDS`] counter observes all of them;
/// multi-pattern evaluations call this once and pass the result to the
/// `_with_context` forms for every pattern.
pub(crate) fn prepared_context<'a>(
    input: &DescriptorInput<'a>,
    function: &'static str,
) -> DescriptorResult<QueryMatchContext<'a>> {
    #[cfg(test)]
    CONTEXT_BUILDS.with(|count| count.set(count.get() + 1));
    // AddHs preserves authoritative ring memberships without appending empty
    // rows for new leaf hydrogens. Reuse the source-semantic sparse adapter;
    // Some(valence) keeps the supplied final assignment borrowed, never rebuilt.
    // Exact-size detached validation remains unchanged in the search owner.
    build_ring_query_match_context(input.topology(), input.ring_info(), Some(input.valence()))
        .map_err(|source| DescriptorError::Search {
            function,
            source: DescriptorSearchCause::Context(source),
        })
}

/// Acquire the retained, once-compiled query for a fixed pattern string.
///
/// The registry keeps one `Arc<QueryGraph>` per distinct pattern for the
/// process lifetime (flyweight `no_tracking`), so the first caller pays the
/// SMARTS parse and every later caller pays only a map lookup plus a
/// reference-count clone.
pub(crate) fn retained_pattern(
    function: &'static str,
    pattern: &'static str,
) -> DescriptorResult<Arc<QueryGraph>> {
    // RDKit✔️✔️: class ss_matcher {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   ss_matcher(const std::string &pattern) : m_pattern(pattern) {
    // RDKit✔️✔️:     m_needCopies = (pattern.find_first_of("$") != std::string::npos);
    // RDKit✔️✔️:     RDKit::RWMol *p = RDKit::SmartsToMol(pattern);
    // RDKit✔️✔️:     m_matcher = p;
    // RDKit✔️✔️:     POSTCONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:   };
    // RDKit✔️✔️:   const RDKit::ROMol *getMatcher() const { return m_matcher; };
    // RDKit✔️✔️:
    // RDKit✔️✔️:  private:
    // RDKit✔️✔️:   ss_matcher() : m_pattern("") {};
    // RDKit✔️✔️:   std::string m_pattern;
    // RDKit✔️✔️:   bool m_needCopies{false};
    // RDKit✔️✔️:   const RDKit::ROMol *m_matcher{nullptr};
    // RDKit✔️✔️: };
    // RDKit✔️✔️: }  // namespace
    // RDKit✔️✔️:
    // RDKit✔️✔️: typedef boost::flyweight<boost::flyweights::key_value<std::string, ss_matcher>,
    // RDKit✔️✔️:                          boost::flyweights::no_tracking>
    // RDKit✔️✔️:     pattern_flyweight;
    // RDKit✔️✔️: #define SMARTSCOUNTFUNC(nm, pattern, vers)         \
    // RDKit✔️✔️:   const std::string nm##Version = vers;            \
    // RDKit✔️✔️:   unsigned int calc##nm(const RDKit::ROMol &mol) { \
    // RDKit✔️✔️:     pattern_flyweight m(pattern);                  \
    // RDKit✔️✔️:     return m.get().countMatches(mol);              \
    // RDKit✔️✔️:   }                                                \
    // RDKit✔️✔️:   extern int no_such_variable
    //
    // Behavioral notes: the source parses once per distinct pattern
    // (flyweight) and aborts the process on a failed parse
    // (`POSTCONDITION`); COSMolKit parses once per process on success and
    // returns the typed `DescriptorSearchCause::Compile` error instead of
    // aborting. Complexity review: a mutex-guarded hash lookup plus an
    // `Arc` clone per call is the same cost class as the flyweight handle;
    // no repeated scanning, no pattern-size allocations, no hot-path loss.
    let mut store = pattern_store()
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner());
    if let Some(query) = store.get(pattern) {
        return Ok(Arc::clone(query));
    }
    #[cfg(test)]
    PATTERN_COMPILES.with(|count| count.set(count.get() + 1));
    let parsed = parse_smarts(pattern, &SmartsParseParams::default()).map_err(|source| {
        DescriptorError::Search {
            function,
            source: DescriptorSearchCause::Compile(source),
        }
    })?;
    let query = Arc::new(parsed);
    store.insert(pattern, Arc::clone(&query));
    Ok(query)
}

/// Ordered pattern matches using one already-prepared shared query context.
///
/// Multi-pattern descriptor units call this form with a single context so
/// every pattern (including recursive ones) reuses the same target-derived
/// preparation, matching the source's reuse of ROMol-internal query state.
/// The returned matches keep the matcher's deterministic order and default
/// `SubstructMatchParameters` uniquify semantics; each `MatchResult`
/// `atom_mapping` is indexed by query atom and holds target atom indices
/// (the projection of RDKit `MatchVectType`'s `mIt.second`).
pub(crate) fn pattern_matches_with_context(
    input: &DescriptorInput<'_>,
    function: &'static str,
    pattern: &'static str,
    context: &QueryMatchContext<'_>,
) -> DescriptorResult<Vec<MatchResult>> {
    // Prepared-input wrapper: constructs the SearchTarget from the
    // prepared DescriptorInput rows and DELEGATES the actual fixed-pattern
    // matching to the ONE common owner `pattern_matches_on_target`, where
    // the complete pinned countMatches body, its behavior review and the
    // recursive-query-copy improvement note live. This wrapper only
    // constructs the target view and delegates; it owns no matching loop.
    let target = SearchTarget::new(
        input.topology(),
        input.coordinates(),
        &input.topology().stereo_groups,
        Some(input.ring_info()),
        Some(input.valence()),
    );
    pattern_matches_on_target(&target, function, pattern, context)
}

/// The ONE common fixed-pattern match owner over explicit target rows.
///
/// Both the prepared-input path ([`pattern_matches_with_context`]) and the
/// topology-only path route here: acquire the retained query and run the
/// default-parameter match over the supplied `SearchTarget`. The caller owns
/// context construction (prepared chemistry rows vs the narrow
/// topology-only context) and ring/valence row supply.
pub(crate) fn pattern_matches_on_target(
    target: &SearchTarget<'_>,
    function: &'static str,
    pattern: &'static str,
    context: &QueryMatchContext<'_>,
) -> DescriptorResult<Vec<MatchResult>> {
    // BEGIN RDKIT CPP FUNCTION: Lipinski.cpp ss_matcher::countMatches — the
    // ONE actual fixed-pattern matching owner reached by BOTH the prepared
    // path and the topology-only path. Callers construct the SearchTarget
    // and the context; this owner acquires the retained query and runs the
    // default-parameter match.
    // RDKit✔️✔️:   unsigned int countMatches(const RDKit::ROMol &mol) const {
    // RDKit✔️✔️:     PRECONDITION(m_matcher, "no matcher");
    // RDKit✔️✔️:     std::vector<RDKit::MatchVectType> matches;
    // RDKit✔️✔️:     // This is an ugly one. Recursive queries aren't thread safe.
    // RDKit✔️✔️:     // Unfortunately we have to take a performance hit here in order
    // RDKit✔️✔️:     // to guarantee thread safety
    // RDKit✔️✔️:     if (m_needCopies) {
    // RDKit✔️✔️:       const RDKit::ROMol nm(*(m_matcher), true);
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       const RDKit::ROMol &nm = *m_matcher;
    // RDKit✔️✔️:       RDKit::SubstructMatch(mol, nm, matches);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return matches.size();
    // RDKit✔️✔️:   }
    //
    // Behavior review (complete pinned body above): the three-argument
    // SubstructMatch overload runs with default SubstructMatchParameters —
    // reproduced by SubstructMatchParams::default() (uniquify=true,
    // maxMatches=1000, recursionPossible=true, useChirality=false,
    // useQueryQueryMatches=false); the result is the ordered match vector
    // whose length the counting callers use. The caller supplies the
    // target-derived context (prepared chemistry rows OR the narrow
    // topology-only context); this owner performs no context construction
    // itself and no chemistry computation.
    //
    // Cost review: one retained-query acquisition (map hit after the first
    // process call, Arc clone) plus one matcher pass over the query
    // pattern against the caller's target. The source's m_needCopies
    // deep-copies the whole query mol on EVERY call for recursive-query
    // thread safety (Lipinski.cpp:43-47); the Rust matcher shares one
    // immutable Arc<QueryGraph> and builds its per-call recursive-query
    // cache instead — this is ONLY a narrow recursive-query-copy
    // improvement: no query copy is ever made
    // while the observable match set is identical. No other whole-operation
    // cost-equivalence claim is made here; matcher-internal costs belong to
    // the search crate's own reviews.
    let query = retained_pattern(function, pattern)?;
    try_get_substruct_matches_with_params_and_context(
        target,
        &query,
        &SubstructMatchParams::default(),
        context,
    )
    .map_err(|source| DescriptorError::Search {
        function,
        source: DescriptorSearchCause::Match(source),
    })
}

/// Topology-only fixed-pattern match count for element-only patterns.
///
/// Builds the ONE narrow topology query context (validation + borrowed
/// adjacency, ring_info None, valence None) and delegates to the ONE
/// common match owner with default `SubstructMatchParameters`. The
/// caller-supplied topology is structurally validated; NO ring/valence
/// rows are computed, borrowed or fabricated.
pub(crate) fn count_pattern_matches_topology_only(
    topology: &cosmolkit_model::TopologyBlock,
    function: &'static str,
    pattern: &'static str,
) -> DescriptorResult<u32> {
    let context =
        build_topology_query_match_context(topology).map_err(|source| DescriptorError::Search {
            function,
            source: DescriptorSearchCause::Context(source),
        })?;
    let matches = {
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let target = SearchTarget::new(topology, &coordinates, &topology.stereo_groups, None, None);
        pattern_matches_on_target(&target, function, pattern, &context)?
    };
    u32::try_from(matches.len()).map_err(|_| DescriptorError::CountOverflow {
        function,
        field: "matches",
    })
}

/// Narrow borrowed-valence fixed-pattern match count (HBD-PUBLIC).
///
/// Builds the ONE narrow borrowed-valence query context (topology +
/// valence-row validation, borrowed adjacency/valence, ring_info None)
/// and delegates to the ONE common match owner with default
/// `SubstructMatchParameters`. The supplied prepared valence rows are
/// BORROWED — never recomputed, extended or cloned — and NO ring rows are
/// computed, borrowed or fabricated (the HBD pattern reads H/valence
/// predicates only). The empty coordinate block is an UNREAD interface
/// value of `SearchTarget`, not fabricated chemistry or cache state.
pub(crate) fn count_pattern_matches_with_valence(
    topology: &cosmolkit_model::TopologyBlock,
    valence: &cosmolkit_core::ValenceAssignment,
    function: &'static str,
    pattern: &'static str,
) -> DescriptorResult<u32> {
    let context = build_valence_query_match_context(topology, valence).map_err(|source| {
        DescriptorError::Search {
            function,
            source: DescriptorSearchCause::Context(source),
        }
    })?;
    let matches = {
        let coordinates = cosmolkit_model::CoordinateBlock::default();
        let target = SearchTarget::new(
            topology,
            &coordinates,
            &topology.stereo_groups,
            None,
            Some(valence),
        );
        pattern_matches_on_target(&target, function, pattern, &context)?
    };
    u32::try_from(matches.len()).map_err(|_| DescriptorError::CountOverflow {
        function,
        field: "matches",
    })
}

/// Count pattern matches using one already-prepared shared query context.
///
/// Thin `.len()` form of [`pattern_matches_with_context`]; multi-pattern
/// units share one context across every pattern through these forms.
pub(crate) fn count_pattern_matches_with_context(
    input: &DescriptorInput<'_>,
    function: &'static str,
    pattern: &'static str,
    context: &QueryMatchContext<'_>,
) -> DescriptorResult<u32> {
    let matches = pattern_matches_with_context(input, function, pattern, context)?;
    u32::try_from(matches.len()).map_err(|_| DescriptorError::CountOverflow {
        function,
        field: "matches",
    })
}

/// Count pattern matches, preparing the shared query context once.
///
/// Single-pattern units use this convenience form; it builds one prepared
/// context from the input's final rows and delegates to the shared-context
/// entry. Units counting several patterns against one input build the
/// context once themselves and call the `_with_context` form for every
/// pattern.
pub(crate) fn count_pattern_matches(
    input: &DescriptorInput<'_>,
    function: &'static str,
    pattern: &'static str,
) -> DescriptorResult<u32> {
    let context = prepared_context(input, function)?;
    count_pattern_matches_with_context(input, function, pattern, &context)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_core::{RingFindType, RingInfo, ValenceModel};
    use cosmolkit_model::{CoordinateBlock, MoleculeProperties};
    use cosmolkit_search::QueryMatchContextError;

    #[test]
    fn descriptor_sparse_ring_context_retains_invalid_state_errors() {
        let topology = cosmolkit_smiles::parse_smiles("CCO", &Default::default())
            .unwrap()
            .topology;
        let valence =
            cosmolkit_core::assign_valence_for_topology(&topology, ValenceModel::RdkitLike)
                .unwrap();
        let coordinates = CoordinateBlock::default();
        let properties = MoleculeProperties::default();
        let count = |rings: &RingInfo, valence: &cosmolkit_core::ValenceAssignment| {
            crate::num_hba_prepared(&DescriptorInput::new(
                &topology,
                &coordinates,
                &properties,
                valence,
                rings,
            ))
        };
        let mut sparse = RingInfo::new(RingFindType::Sssr, 0, 0);
        assert_eq!(count(&sparse, &valence).unwrap(), 1);
        sparse.reset();
        assert!(matches!(
            count(&sparse, &valence),
            Err(DescriptorError::Search {
                source: DescriptorSearchCause::Context(QueryMatchContextError::UninitializedRings),
                ..
            })
        ));
        let too_long = RingInfo::new(RingFindType::Sssr, 4, 2);
        assert!(matches!(
            count(&too_long, &valence),
            Err(DescriptorError::Search {
                source: DescriptorSearchCause::Context(
                    QueryMatchContextError::RingMembershipRows {
                        field: "atoms",
                        expected: 3,
                        actual: 4,
                    }
                ),
                ..
            })
        ));
        let sparse = RingInfo::new(RingFindType::Sssr, 0, 0);
        let mut incomplete = valence;
        incomplete.implicit_hydrogens.pop();
        assert!(matches!(
            count(&sparse, &incomplete),
            Err(DescriptorError::Search {
                source: DescriptorSearchCause::Context(QueryMatchContextError::ValenceRows {
                    field: "implicit_hydrogens",
                    expected: 3,
                    actual: 2,
                }),
                ..
            })
        ));
    }
}
