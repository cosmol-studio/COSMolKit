//! Subgraph isomorphism matching (VF2) for molecule pattern matching.
//!
//! ## RDKit provenance (protocol: dev/source_reproduction_protocol.md)
//!
//! This module reproduces RDKit's substructure matching from:
//! - `third_party/rdkit/Code/GraphMol/Substruct/vf2.hpp` (~682 lines C++)
//! - `third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp` (~735 lines C++)
//!
//! The VF2 algorithm implementation is adapted from vflib-2.0 by P. Foggia,
//! extensively modified by Greg Landrum, ported to Rust with depth-based
//! term_1/term_2 tracking (BackTrack decrements counters instead of
//! recomputing from scratch).
//!
//! ## Marker convention
//!
//! Each copied C++ block below uses the two-axis status marker:
//! - RDKit✔️✔️: fully reproduced behavior and performance
//! - RDKit✔️❌: functionally correct, but with a known performance gap
//! - RDKit❗✔️: unfinished behavior that must not be presented as parity
//! - RDKit❌❌: not yet ported

use crate::query_behavior::{
    QueryMatchContext, and_query_match, atom_predicate_matches_with_context,
    bond_predicate_matches_with_context, build_query_match_context, or_query_match,
    query_atom_query_match, query_bond_query_match, xor_query_match,
};
use crate::{AtomQueryPredicate, BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryNode};
use crate::{SearchTarget, SearchTargetAccess};
use cosmolkit_core::PeriodicTableError;
use cosmolkit_model::{
    Atom, Bond, Conformer3D, NeighborRef, PropertyValue, StereoGroupKind, TopologyBlock,
};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag};
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of a single substructure match.
pub type SubstructMatchResult = crate::MatchResult;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SubstructMatchError {
    #[error(
        "RDKit substructure matching branch {branch} is unsupported until {rdkit_function} is source-ported"
    )]
    Unsupported {
        branch: &'static str,
        rdkit_function: &'static str,
    },
    #[error(transparent)]
    PeriodicTable(#[from] PeriodicTableError),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstructMatchOverload {
    DetachedTarget,
    MolBundle,
    ResonanceMolSupplier,
    SubstructLibrary,
}

pub fn check_substruct_match_overload_support(
    overload: SubstructMatchOverload,
) -> Result<(), SubstructMatchError> {
    match overload {
        SubstructMatchOverload::DetachedTarget => Ok(()),
        SubstructMatchOverload::MolBundle => Err(SubstructMatchError::Unsupported {
            branch: "MolBundle substructure-match overloads",
            rdkit_function: "SubstructMatch(MolBundle, ROMol/MolBundle, params)",
        }),
        SubstructMatchOverload::ResonanceMolSupplier => Err(SubstructMatchError::Unsupported {
            branch: "resonance substructure-match overload",
            rdkit_function: "SubstructMatch(ResonanceMolSupplier, ROMol, params)",
        }),
        SubstructMatchOverload::SubstructLibrary => Err(SubstructMatchError::Unsupported {
            branch: "SubstructLibrary search overloads",
            rdkit_function: "SubstructLibrary::getMatches/hasMatch/countMatches",
        }),
    }
}

#[derive(Debug, thiserror::Error)]
pub enum SubstructMatchParamsJsonError {
    #[error("invalid substructure match parameter JSON: {0}")]
    InvalidJson(#[from] serde_json::Error),
    #[error("invalid JSON value for substructure match parameter '{field}'")]
    InvalidField { field: &'static str },
}

type SubstructMatchResultList = Result<Vec<SubstructMatchResult>, SubstructMatchError>;

/// Query input accepted by the matcher boundary.
///
/// `QueryGraph` is the canonical representation. Concrete target values are
/// deliberately not accepted as query input: lowering belongs at the facade.
pub trait QueryInput {
    fn query_graph(&self) -> Result<Cow<'_, QueryGraph>, SubstructMatchError>;
}

impl QueryInput for QueryGraph {
    fn query_graph(&self) -> Result<Cow<'_, QueryGraph>, SubstructMatchError> {
        Ok(Cow::Borrowed(self))
    }
}

/// Parameters controlling substructure matching behaviour.
pub type ExtraAtomCheck =
    Arc<dyn for<'a> Fn(&QueryGraph, &QueryAtom, &SearchTarget<'a>, &Atom) -> bool + Send + Sync>;
pub type ExtraBondCheck = Arc<dyn Fn(&Bond, &Bond) -> bool + Send + Sync>;
pub type ExtraFinalCheck = Arc<dyn for<'a> Fn(&SearchTarget<'a>, &[usize]) -> bool + Send + Sync>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AtomCoordsMatchFunctor {
    pub ref_conf_id: i32,
    pub query_conf_id: i32,
    pub tol2: f64,
}

impl AtomCoordsMatchFunctor {
    #[must_use]
    pub fn new(ref_conf_id: i32, query_conf_id: i32, tolerance: f64) -> Self {
        Self {
            ref_conf_id,
            query_conf_id,
            tol2: tolerance * tolerance,
        }
    }

    #[must_use]
    pub fn matches(
        &self,
        query_mol: &QueryGraph,
        query_atom: &QueryAtom,
        target_mol: &SearchTarget<'_>,
        target_atom: &Atom,
    ) -> bool {
        // RDKit✔️✔️: bool AtomCoordsMatchFunctor::operator()(const Atom &queryAtom,
        // RDKit✔️✔️:                                         const Atom &targetAtom) const {
        // RDKit✔️✔️:   if (!queryAtom.getOwningMol().getNumConformers() ||
        // RDKit✔️✔️:       !targetAtom.getOwningMol().getNumConformers()) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   const auto &queryPos = queryAtom.getOwningMol()
        // RDKit✔️✔️:                              .getConformer(d_queryConfId)
        // RDKit✔️✔️:                              .getAtomPos(queryAtom.getIdx());
        // RDKit✔️✔️:   const auto &targetPos = targetAtom.getOwningMol()
        // RDKit✔️✔️:                               .getConformer(d_refConfId)
        // RDKit✔️✔️:                               .getAtomPos(targetAtom.getIdx());
        // RDKit✔️✔️:   return (queryPos - targetPos).lengthSq() <= d_tol2;
        // RDKit✔️✔️: };
        // Complexity review: both versions select two conformers, index two
        // coordinate rows, and compare three squared deltas in O(1) after the
        // conformer-id lookup. No coordinate data is cloned or allocated.
        fn conformer<T>(molecule: &T, id: i32) -> Option<&Conformer3D>
        where
            T: ConformerSource,
        {
            let conformers = molecule.conformers_3d();
            if id < 0 {
                conformers.first()
            } else {
                conformers
                    .iter()
                    .find(|conformer| conformer.id() == id as usize)
            }
        }

        let Some(query_conformer) = conformer(query_mol, self.query_conf_id) else {
            return false;
        };
        let Some(target_conformer) = conformer(target_mol, self.ref_conf_id) else {
            return false;
        };
        let Some(query_position) = query_conformer.coordinates().get(query_atom.id().index())
        else {
            return false;
        };
        let Some(target_position) = target_conformer.coordinates().get(target_atom.id().index())
        else {
            return false;
        };
        query_position
            .iter()
            .zip(target_position)
            .map(|(query, target)| {
                let delta = query - target;
                delta * delta
            })
            .sum::<f64>()
            <= self.tol2
    }
}

trait ConformerSource {
    fn conformers_3d(&self) -> &[Conformer3D];
}

impl ConformerSource for SearchTarget<'_> {
    fn conformers_3d(&self) -> &[Conformer3D] {
        self.coordinate_block().conformers_3d.as_slice()
    }
}

impl ConformerSource for QueryGraph {
    fn conformers_3d(&self) -> &[Conformer3D] {
        self.conformers_3d()
    }
}

impl Default for AtomCoordsMatchFunctor {
    fn default() -> Self {
        Self::new(-1, -1, 1e-4)
    }
}

#[derive(Clone)]
pub struct SubstructMatchParams {
    /// Maximum number of matches to return (default: 1000).
    pub max_matches: usize,
    /// Whether to uniquify results (default: true).
    pub uniquify: bool,
    /// Whether atom/bond stereochemistry participates in matching.
    pub use_chirality: bool,
    /// Whether enhanced stereo groups participate in final matching.
    pub use_enhanced_stereo: bool,
    /// Whether specified query stereo may match unspecified molecule stereo.
    pub specified_stereo_query_matches_unspecified: bool,
    /// Whether two query atoms are compared as query trees.
    pub use_query_query_matches: bool,
    /// Whether recursive query nodes may be evaluated.
    pub recursion_possible: bool,
    /// Maximum matches used while evaluating recursive query nodes.
    pub max_recursive_matches: usize,
    /// Requested matcher thread count; matching is currently single-threaded.
    pub num_threads: i32,
    /// Whether aromatic bonds may match conjugated single or double bonds.
    pub aromatic_matches_conjugated: bool,
    /// Whether aromatic bonds may match any single or double bond.
    pub aromatic_matches_single_or_double: bool,
    /// Atom property names that must have equal string values on both atoms.
    pub atom_properties: Vec<String>,
    /// Bond property names that must have equal string values on both bonds.
    pub bond_properties: Vec<String>,
    /// Optional caller-provided atom compatibility check.
    pub extra_atom_check: Option<ExtraAtomCheck>,
    /// Whether `extra_atom_check` replaces the default atom comparison.
    pub extra_atom_check_overrides_default_check: bool,
    /// Optional caller-provided bond compatibility check.
    pub extra_bond_check: Option<ExtraBondCheck>,
    /// Whether `extra_bond_check` replaces the default bond comparison.
    pub extra_bond_check_overrides_default_check: bool,
    /// Whether generic-group labels participate in final matching.
    pub use_generic_matchers: bool,
    /// Optional caller-provided final check over target atom indices.
    pub extra_final_check: Option<ExtraFinalCheck>,
}

impl fmt::Debug for SubstructMatchParams {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SubstructMatchParams")
            .field("max_matches", &self.max_matches)
            .field("uniquify", &self.uniquify)
            .field("use_chirality", &self.use_chirality)
            .field("use_enhanced_stereo", &self.use_enhanced_stereo)
            .field(
                "specified_stereo_query_matches_unspecified",
                &self.specified_stereo_query_matches_unspecified,
            )
            .field("use_query_query_matches", &self.use_query_query_matches)
            .field("recursion_possible", &self.recursion_possible)
            .field("max_recursive_matches", &self.max_recursive_matches)
            .field("num_threads", &self.num_threads)
            .field(
                "aromatic_matches_conjugated",
                &self.aromatic_matches_conjugated,
            )
            .field(
                "aromatic_matches_single_or_double",
                &self.aromatic_matches_single_or_double,
            )
            .field("atom_properties", &self.atom_properties)
            .field("bond_properties", &self.bond_properties)
            .field("extra_atom_check", &self.extra_atom_check.is_some())
            .field(
                "extra_atom_check_overrides_default_check",
                &self.extra_atom_check_overrides_default_check,
            )
            .field("extra_bond_check", &self.extra_bond_check.is_some())
            .field(
                "extra_bond_check_overrides_default_check",
                &self.extra_bond_check_overrides_default_check,
            )
            .field("use_generic_matchers", &self.use_generic_matchers)
            .field("extra_final_check", &self.extra_final_check.is_some())
            .finish()
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum RecursiveQueryCacheKey {
    Serial(u32),
    OwnedQuery(usize),
}

type RecursiveQueryMatchCache = BTreeMap<RecursiveQueryCacheKey, Vec<bool>>;

struct RecursiveLocker {
    cache: RecursiveQueryMatchCache,
}

impl RecursiveLocker {
    fn new(query: &QueryGraph, recursion_possible: bool) -> Self {
        // RDKit✔️🔝: RecursiveLocker(const ROMol &query, const bool recursionPossible) {
        // RDKit✔️🔝:   if (recursionPossible) {
        // RDKit✔️🔝:     locked.reserve(query.getNumAtoms());
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // Rust keeps recursive match state in this call-local cache instead of
        // mutating and locking query nodes. This preserves the source lifetime
        // semantics while avoiding the O(query atoms) pointer-vector reserve
        // and every mutex operation. The query and flag remain inputs here so
        // this constructor is the canonical source boundary.
        let _ = (query, recursion_possible);
        Self {
            cache: RecursiveQueryMatchCache::new(),
        }
    }
}

impl Drop for RecursiveLocker {
    fn drop(&mut self) {
        // RDKit✔️✔️: ~RecursiveLocker() {
        // RDKit✔️✔️:   for (auto v : locked) {
        // RDKit✔️✔️:     v->clear();
        // RDKit✔️✔️: #ifdef RDK_BUILD_THREADSAFE_SSS
        // RDKit✔️✔️:     v->d_mutex.unlock();
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Complexity review: dropping the call-local cache clears each stored
        // atom-membership vector once, matching RDKit's linear clear pass.
        // No unlock is required because immutable query nodes are never shared
        // mutably; ownership enforces the same cleanup on every return path.
        self.cache.clear();
    }
}

fn recursive_query_cache_key(
    query: &crate::query_behavior::RecursiveStructureQuery,
) -> RecursiveQueryCacheKey {
    if query.serial_number() != 0 {
        RecursiveQueryCacheKey::Serial(query.serial_number())
    } else {
        RecursiveQueryCacheKey::OwnedQuery(query as *const _ as usize)
    }
}

impl Default for SubstructMatchParams {
    fn default() -> Self {
        // RDKit✔️✔️: bool useChirality = false;  //!< Use chirality in determining whether or not
        // RDKit✔️✔️:                             //!< atoms/bonds match
        // RDKit✔️✔️: bool uniquify = true;            //!< uniquify (by atom index) match results
        // RDKit✔️✔️: unsigned int maxMatches = 1000;  //!< maximum number of matches to return
        // RDKit✔️✔️: bool specifiedStereoQueryMatchesUnspecified =
        // RDKit✔️✔️:     false;  //!< If set, query atoms and bonds with specified stereochemistry
        // RDKit✔️✔️:             //!< will match atoms and bonds with unspecified stereochemistry
        // RDKit✔️✔️: bool useEnhancedStereo = false;
        // RDKit✔️✔️: bool aromaticMatchesConjugated = false;
        // RDKit✔️✔️: bool useQueryQueryMatches = false;
        // RDKit✔️✔️: bool useGenericMatchers = false;
        // RDKit✔️✔️: bool recursionPossible = true;
        // RDKit✔️✔️: int numThreads = 1;
        // RDKit✔️✔️: std::vector<std::string> atomProperties;
        // RDKit✔️✔️: std::vector<std::string> bondProperties;
        // RDKit✔️✔️: std::function<bool(const ROMol &, std::span<const unsigned int>)>
        // RDKit✔️✔️:     extraFinalCheck;
        // RDKit✔️✔️: unsigned int maxRecursiveMatches = 1000;
        // RDKit✔️✔️: bool aromaticMatchesSingleOrDouble = false;
        // RDKit✔️✔️: std::function<bool(const Atom &, const Atom &)> extraAtomCheck;
        // RDKit✔️✔️: bool extraAtomCheckOverridesDefaultCheck = false;
        // RDKit✔️✔️: std::function<bool(const Bond &, const Bond &)> extraBondCheck;
        // RDKit✔️✔️: bool extraBondCheckOverridesDefaultCheck = false;
        // Complexity review: initialization is O(1) and allocates only empty
        // Vec headers and absent callback slots, matching the C++ defaults.
        Self {
            max_matches: 1000,
            uniquify: true,
            use_chirality: false,
            use_enhanced_stereo: false,
            specified_stereo_query_matches_unspecified: false,
            use_query_query_matches: false,
            recursion_possible: true,
            max_recursive_matches: 1000,
            num_threads: 1,
            aromatic_matches_conjugated: false,
            aromatic_matches_single_or_double: false,
            atom_properties: Vec::new(),
            bond_properties: Vec::new(),
            extra_atom_check: None,
            extra_atom_check_overrides_default_check: false,
            extra_bond_check: None,
            extra_bond_check_overrides_default_check: false,
            use_generic_matchers: false,
            extra_final_check: None,
        }
    }
}

fn json_param_bool(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<bool>, SubstructMatchParamsJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    match value {
        serde_json::Value::Bool(value) => Ok(Some(*value)),
        serde_json::Value::Number(value) if value.as_u64() == Some(1) => Ok(Some(true)),
        serde_json::Value::Number(value) if value.as_u64() == Some(0) => Ok(Some(false)),
        serde_json::Value::String(value) if value == "true" || value == "1" => Ok(Some(true)),
        serde_json::Value::String(value) if value == "false" || value == "0" => Ok(Some(false)),
        _ => Err(SubstructMatchParamsJsonError::InvalidField { field }),
    }
}

fn json_param_usize(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<usize>, SubstructMatchParamsJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    let parsed = match value {
        serde_json::Value::Number(value) => value
            .as_u64()
            .and_then(|value| u32::try_from(value).ok())
            .map(|value| value as usize),
        serde_json::Value::String(value) => value.parse::<u32>().ok().map(|value| value as usize),
        _ => None,
    };
    parsed
        .map(Some)
        .ok_or(SubstructMatchParamsJsonError::InvalidField { field })
}

fn json_param_i32(
    object: &serde_json::Map<String, serde_json::Value>,
    field: &'static str,
) -> Result<Option<i32>, SubstructMatchParamsJsonError> {
    let Some(value) = object.get(field) else {
        return Ok(None);
    };
    let parsed = match value {
        serde_json::Value::Number(value) => {
            value.as_i64().and_then(|value| i32::try_from(value).ok())
        }
        serde_json::Value::String(value) => value.parse().ok(),
        _ => None,
    };
    parsed
        .map(Some)
        .ok_or(SubstructMatchParamsJsonError::InvalidField { field })
}

pub fn update_substruct_match_params_from_json(
    params: &mut SubstructMatchParams,
    json: &str,
) -> Result<(), SubstructMatchParamsJsonError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: updateSubstructMatchParamsFromJSON
    // RDKit✔️✔️: void updateSubstructMatchParamsFromJSON(SubstructMatchParameters &params,
    // RDKit✔️✔️:                                         const std::string &json) {
    // RDKit✔️✔️:   if (json.empty()) {
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::istringstream ss;
    // RDKit✔️✔️:   ss.str(json);
    // RDKit✔️✔️:   boost::property_tree::ptree pt;
    // RDKit✔️✔️:   boost::property_tree::read_json(ss, pt);
    // RDKit✔️✔️:   PT_OPT_GET(useChirality);
    // RDKit✔️✔️:   PT_OPT_GET(useEnhancedStereo);
    // RDKit✔️✔️:   PT_OPT_GET(aromaticMatchesConjugated);
    // RDKit✔️✔️:   PT_OPT_GET(useQueryQueryMatches);
    // RDKit✔️✔️:   PT_OPT_GET(recursionPossible);
    // RDKit✔️✔️:   PT_OPT_GET(uniquify);
    // RDKit✔️✔️:   PT_OPT_GET(maxMatches);
    // RDKit✔️✔️:   PT_OPT_GET(maxRecursiveMatches);
    // RDKit✔️✔️:   PT_OPT_GET(numThreads);
    // RDKit✔️✔️:   PT_OPT_GET(specifiedStereoQueryMatchesUnspecified);
    // RDKit✔️✔️:   PT_OPT_GET(aromaticMatchesSingleOrDouble);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Local complexity review: both parsers are O(input length), followed by
    // eleven expected O(1) object lookups. No molecule/query data is touched
    // or cloned. Each successful field is assigned immediately in the source
    // PT_OPT_GET order, so a later conversion error preserves earlier updates.
    if json.is_empty() {
        return Ok(());
    }
    let value: serde_json::Value = serde_json::from_str(json)?;
    let object = value
        .as_object()
        .ok_or(SubstructMatchParamsJsonError::InvalidField { field: "root" })?;
    if let Some(value) = json_param_bool(object, "useChirality")? {
        params.use_chirality = value;
    }
    if let Some(value) = json_param_bool(object, "useEnhancedStereo")? {
        params.use_enhanced_stereo = value;
    }
    if let Some(value) = json_param_bool(object, "aromaticMatchesConjugated")? {
        params.aromatic_matches_conjugated = value;
    }
    if let Some(value) = json_param_bool(object, "useQueryQueryMatches")? {
        params.use_query_query_matches = value;
    }
    if let Some(value) = json_param_bool(object, "recursionPossible")? {
        params.recursion_possible = value;
    }
    if let Some(value) = json_param_bool(object, "uniquify")? {
        params.uniquify = value;
    }
    if let Some(value) = json_param_usize(object, "maxMatches")? {
        params.max_matches = value;
    }
    if let Some(value) = json_param_usize(object, "maxRecursiveMatches")? {
        params.max_recursive_matches = value;
    }
    if let Some(value) = json_param_i32(object, "numThreads")? {
        params.num_threads = value;
    }
    if let Some(value) = json_param_bool(object, "specifiedStereoQueryMatchesUnspecified")? {
        params.specified_stereo_query_matches_unspecified = value;
    }
    if let Some(value) = json_param_bool(object, "aromaticMatchesSingleOrDouble")? {
        params.aromatic_matches_single_or_double = value;
    }
    Ok(())
}

pub fn substruct_match_params_to_json(params: &SubstructMatchParams) -> String {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: substructMatchParamsToJSON
    // RDKit✔️✔️: std::string substructMatchParamsToJSON(const SubstructMatchParameters &params) {
    // RDKit✔️✔️:   boost::property_tree::ptree pt;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PT_OPT_PUT(useChirality);
    // RDKit✔️✔️:   PT_OPT_PUT(useEnhancedStereo);
    // RDKit✔️✔️:   PT_OPT_PUT(aromaticMatchesConjugated);
    // RDKit✔️✔️:   PT_OPT_PUT(useQueryQueryMatches);
    // RDKit✔️✔️:   PT_OPT_PUT(recursionPossible);
    // RDKit✔️✔️:   PT_OPT_PUT(uniquify);
    // RDKit✔️✔️:   PT_OPT_PUT(maxMatches);
    // RDKit✔️✔️:   PT_OPT_PUT(maxRecursiveMatches);
    // RDKit✔️✔️:   PT_OPT_PUT(numThreads);
    // RDKit✔️✔️:   PT_OPT_PUT(specifiedStereoQueryMatchesUnspecified);
    // RDKit✔️✔️:   PT_OPT_PUT(aromaticMatchesSingleOrDouble);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::stringstream ss;
    // RDKit✔️✔️:   boost::property_tree::json_parser::write_json(ss, pt);
    // RDKit✔️✔️:   return ss.str();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Local complexity review: both implementations serialize the same fixed
    // eleven scalar fields in O(output length), without molecule/query work.
    format!(
        concat!(
            "{{\n",
            "    \"useChirality\": \"{}\",\n",
            "    \"useEnhancedStereo\": \"{}\",\n",
            "    \"aromaticMatchesConjugated\": \"{}\",\n",
            "    \"useQueryQueryMatches\": \"{}\",\n",
            "    \"recursionPossible\": \"{}\",\n",
            "    \"uniquify\": \"{}\",\n",
            "    \"maxMatches\": \"{}\",\n",
            "    \"maxRecursiveMatches\": \"{}\",\n",
            "    \"numThreads\": \"{}\",\n",
            "    \"specifiedStereoQueryMatchesUnspecified\": \"{}\",\n",
            "    \"aromaticMatchesSingleOrDouble\": \"{}\"\n",
            "}}\n",
        ),
        params.use_chirality,
        params.use_enhanced_stereo,
        params.aromatic_matches_conjugated,
        params.use_query_query_matches,
        params.recursion_possible,
        params.uniquify,
        params.max_matches,
        params.max_recursive_matches,
        params.num_threads,
        params.specified_stereo_query_matches_unspecified,
        params.aromatic_matches_single_or_double,
    )
}

// ---------------------------------------------------------------------------
// Internal minimum-degree graph representation
// ---------------------------------------------------------------------------

/// Minimal adjacency info needed for VF2.
///
/// `nbrs[i]` is a slice into `edges`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct Vf2Graph {
    n_atoms: usize,
    n_bonds: usize,
    /// Canonical source/target endpoints indexed by bond id.
    edge_endpoints: Vec<(usize, usize)>,
    /// For each atom index, the neighbor indices and bond ids.
    adjacency: Vec<Vec<(usize, usize)>>, // (neighbor_atom_index, bond_index)
}

/// A borrowed ordered adjacency row from one of the existing graph owners.
/// Query and compiled rows store `(neighbor, bond_index)` tuples; topology
/// rows store the same observation in `NeighborRef` values.
#[derive(Debug, Clone, Copy)]
enum Vf2NeighborRow<'a> {
    Pairs(&'a [(usize, usize)]),
    NeighborRefs(&'a [NeighborRef]),
}

impl<'a> Vf2NeighborRow<'a> {
    fn len(self) -> usize {
        match self {
            Self::Pairs(row) => row.len(),
            Self::NeighborRefs(row) => row.len(),
        }
    }

    fn get(self, index: usize) -> Option<(usize, usize)> {
        match self {
            Self::Pairs(row) => row.get(index).copied(),
            Self::NeighborRefs(row) => row
                .get(index)
                .map(|neighbor| (neighbor.atom_index, neighbor.bond.index())),
        }
    }

    fn iter(self) -> Vf2NeighborIter<'a> {
        // RDKit❗✔️: typename Graph::out_edge_iterator bNbrs, eNbrs;
        // RDKit❗✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit❗✔️: while (bNbrs != eNbrs) {
        // RDKit❗✔️:       ++bNbrs;
        // RDKit❗✔️: }
        // The Rust observation borrows each owner's existing row and advances
        // in stored order. Tuple rows copy their two indices; topology rows
        // project NeighborRef's canonical atom and bond indices on demand.
        // Complexity: iterator construction/advance does not allocate or copy
        // a row; indexed lookup remains constant time on either slice.
        match self {
            Self::Pairs(row) => Vf2NeighborIter::Pairs(row.iter()),
            Self::NeighborRefs(row) => Vf2NeighborIter::NeighborRefs(row.iter()),
        }
    }
}

enum Vf2NeighborIter<'a> {
    Pairs(std::slice::Iter<'a, (usize, usize)>),
    NeighborRefs(std::slice::Iter<'a, NeighborRef>),
}

impl Iterator for Vf2NeighborIter<'_> {
    type Item = (usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        match self {
            Self::Pairs(iter) => iter.next().copied(),
            Self::NeighborRefs(iter) => iter
                .next()
                .map(|neighbor| (neighbor.atom_index, neighbor.bond.index())),
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        match self {
            Self::Pairs(iter) => iter.size_hint(),
            Self::NeighborRefs(iter) => iter.size_hint(),
        }
    }
}

impl ExactSizeIterator for Vf2NeighborIter<'_> {}

/// Query-side graph plan retained by [`CompiledQuery`].
///
/// This is an intentionally opaque alias at the public boundary: callers can
/// inspect the resulting atom order, but cannot depend on VF2 storage details.
pub(crate) type CompiledQueryGraph = Vf2Graph;

/// Build a VF2-compatible adjacency view from a molecule.
///
/// The C++ code iterates `out_edges` via Boost graph. We pre-build adjacency
/// once and use raw index lookups.
trait Vf2GraphSource {
    fn vf2_num_atoms(&self) -> usize;
    fn vf2_num_bonds(&self) -> usize;
    fn vf2_bond_endpoints(&self, index: usize) -> (usize, usize);
}

impl Vf2GraphSource for QueryGraph {
    fn vf2_num_atoms(&self) -> usize {
        self.num_atoms()
    }
    fn vf2_num_bonds(&self) -> usize {
        self.num_bonds()
    }
    fn vf2_bond_endpoints(&self, index: usize) -> (usize, usize) {
        self.bonds()[index].endpoints()
    }
}

#[cfg(test)]
std::thread_local! {
    static VF2_GRAPH_BUILD_ENTRIES: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}

fn build_vf2_graph<G: Vf2GraphSource>(mol: &G) -> Vf2Graph {
    #[cfg(test)]
    VF2_GRAPH_BUILD_ENTRIES.with(|entries| entries.set(entries.get() + 1));

    // RDKit source (implicit in vf2.hpp usage of out_edges):
    //   The VF2 state stores Graph *g1, *g2 and calls:
    //     boost::out_edges(node, *g)
    //     boost::out_degree(node, *g)
    //     boost::adjacent_vertices(node, *g)
    //   These are all O(1) in Boost adjacency_list.
    //
    // Rust-only compiled-query representation:
    //   This one-time O(V+E) materialization is limited to explicit compiled-
    //   query construction. Ordinary and prepared matching borrow the existing
    //   QueryGraph/TopologyBlock rows through Vf2GraphRef; compiled matching
    //   borrows the retained graph. Indexed row lookup remains O(degree).
    let n_atoms = mol.vf2_num_atoms();
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_atoms];
    let mut edge_endpoints = Vec::with_capacity(mol.vf2_num_bonds());
    for bond_idx in 0..mol.vf2_num_bonds() {
        let (b, e) = mol.vf2_bond_endpoints(bond_idx);
        edge_endpoints.push((b, e));
        adjacency[b].push((e, bond_idx));
        adjacency[e].push((b, bond_idx));
    }
    Vf2Graph {
        n_atoms,
        n_bonds: mol.vf2_num_bonds(),
        edge_endpoints,
        adjacency,
    }
}

fn get_other_idx(g: Vf2GraphRef<'_>, edge: usize, vertex: NodeId) -> NodeId {
    // RDKit✔️✔️: template <class Graph, class VertexDescr, class EdgeDescr>
    // RDKit✔️✔️: VertexDescr getOtherIdx(const Graph &g, const EdgeDescr &edge,
    // RDKit✔️✔️:                         const VertexDescr &vertex) {
    // RDKit✔️✔️:   VertexDescr tmp = boost::source(edge, g);
    // RDKit✔️✔️:   if (tmp == vertex) {
    // RDKit✔️✔️:     tmp = boost::target(edge, g);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return tmp;
    // RDKit✔️✔️: }
    // Complexity review: the endpoint table provides the same O(1) source and
    // target lookup as the Boost edge descriptor, with no per-call allocation.
    let (source, target) = g.bond_endpoints(edge);
    if source == vertex { target } else { source }
}

impl Vf2Graph {
    pub(crate) fn num_atoms(&self) -> usize {
        self.n_atoms
    }

    pub(crate) fn num_bonds(&self) -> usize {
        self.n_bonds
    }

    fn neighbor_row(&self, node: usize) -> Vf2NeighborRow<'_> {
        Vf2NeighborRow::Pairs(&self.adjacency[node])
    }
}

/// Allocation-free observations over the three existing graph owners.
/// Unlike RDKit's concrete Boost graph pointer, this Rust-only enum borrows
/// the owner-specific storage directly and carries no reconstructed graph.
#[derive(Debug, Clone, Copy)]
enum Vf2GraphRef<'a> {
    Query(&'a QueryGraph),
    Target(&'a TopologyBlock),
    Compiled(&'a Vf2Graph),
}

impl<'a> Vf2GraphRef<'a> {
    fn query(graph: &'a QueryGraph) -> Self {
        Self::Query(graph)
    }

    fn target(graph: &'a TopologyBlock) -> Self {
        Self::Target(graph)
    }

    fn compiled(graph: &'a Vf2Graph) -> Self {
        Self::Compiled(graph)
    }

    fn num_atoms(self) -> usize {
        // RDKit❗✔️:   Graph *g1, *g2;
        // RDKit❗✔️:         n1(num_vertices(*ag1)),
        // RDKit❗✔️:         n2(num_vertices(*ag2)) {
        // The variants retain those graph owners by reference; each atom count
        // is a direct field/slice-length observation with no graph creation.
        match self {
            Self::Query(graph) => graph.num_atoms(),
            Self::Target(graph) => graph.atoms.len(),
            Self::Compiled(graph) => graph.num_atoms(),
        }
    }

    fn num_bonds(self) -> usize {
        match self {
            Self::Query(graph) => graph.num_bonds(),
            Self::Target(graph) => graph.bonds.len(),
            Self::Compiled(graph) => graph.num_bonds(),
        }
    }

    fn out_degree(self, node: usize) -> usize {
        // RDKit❗✔️: if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // Each existing row's length is its undirected out-degree, observed
        // directly in O(1) without a temporary row or graph.
        self.neighbor_row(node).len()
    }

    fn neighbor_row(self, node: usize) -> Vf2NeighborRow<'a> {
        // RDKit❗✔️: typename Graph::out_edge_iterator bNbrs, eNbrs;
        // RDKit❗✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit❗✔️: RDK_ADJ_ITER n1iter_beg, n1iter_end;
        // RDKit❗✔️:             boost::adjacent_vertices(pair.n1, *g1);
        // Query and compiled rows already store ordered pairs. A validated
        // target topology borrows ordered NeighborRef rows whose BondId index
        // is the canonical bond-table edge descriptor.
        match self {
            Self::Query(graph) => Vf2NeighborRow::Pairs(&graph.adjacency()[node]),
            Self::Target(graph) => Vf2NeighborRow::NeighborRefs(graph.adjacency.neighbors_of(node)),
            Self::Compiled(graph) => graph.neighbor_row(node),
        }
    }

    fn bond_endpoints(self, edge: usize) -> (usize, usize) {
        // RDKit❗✔️:   VertexDescr tmp = boost::source(edge, g);
        // RDKit❗✔️:   if (tmp == vertex) {
        // RDKit❗✔️:     tmp = boost::target(edge, g);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return tmp;
        // The Rust view observes the descriptor's canonical endpoint pair in
        // O(1); the caller applies the same opposite-endpoint selection.
        match self {
            Self::Query(graph) => graph.bonds()[edge].endpoints(),
            Self::Target(graph) => {
                let bond = &graph.bonds[edge];
                (bond.begin().index(), bond.end().index())
            }
            Self::Compiled(graph) => graph.edge_endpoints[edge],
        }
    }
}

// ---------------------------------------------------------------------------
// Atom and bond matching functors
// ---------------------------------------------------------------------------

fn property_compat(
    properties1: &BTreeMap<String, PropertyValue>,
    properties2: &BTreeMap<String, PropertyValue>,
    properties: &[String],
) -> bool {
    // RDKit✔️🔝: bool propertyCompat(const RDProps *r1, const RDProps *r2,
    // RDKit✔️🔝:                     const std::vector<std::string> &properties) {
    // RDKit✔️🔝:   PRECONDITION(r1, "bad RDProps");
    // RDKit✔️🔝:   PRECONDITION(r2, "bad RDProps");
    // RDKit✔️🔝:
    // RDKit✔️🔝:   for (const auto &prop : properties) {
    // RDKit✔️🔝:     std::string prop1;
    // RDKit✔️🔝:     bool hasprop1 = r1->getPropIfPresent<std::string>(prop, prop1);
    // RDKit✔️🔝:     std::string prop2;
    // RDKit✔️🔝:     bool hasprop2 = r2->getPropIfPresent<std::string>(prop, prop2);
    // RDKit✔️🔝:     if (hasprop1 && hasprop2) {
    // RDKit✔️🔝:       if (prop1 != prop2) {
    // RDKit✔️🔝:         return false;
    // RDKit✔️🔝:       }
    // RDKit✔️🔝:     } else if (hasprop1 || hasprop2) {
    // RDKit✔️🔝:       // only one has the property
    // RDKit✔️🔝:       return false;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    //
    // Typed references make both source pointer preconditions
    // unrepresentable. COSMolKit's canonical atom/bond property maps store
    // only strings, exactly the type requested by the source function, so a
    // pair of Option<&String> values preserves the source's present/missing
    // cases without temporary string copies. Local complexity review: both
    // implementations scan the requested property list once and short-circuit
    // at the first mismatch without cloning or allocating. RDKit's Dict scans
    // its vector of entries for each lookup (O(P*N)); BTreeMap lookup is
    // O(log N), making this O(P*log N) while preserving lookup semantics.
    for property in properties {
        if properties1.get(property) != properties2.get(property) {
            return false;
        }
    }
    true
}

// RDKit source (SubstructMatch.cpp):
//   class AtomLabelFunctor {
//    public:
//     AtomLabelFunctor(const ROMol &query, const ROMol &mol,
//                      const SubstructMatchParameters &ps)
//         : d_query(query), d_mol(mol), d_params(ps) {};
//     bool operator()(unsigned int i, unsigned int j) const {
//       bool res = false;
//       if (d_params.useChirality) {
//         const Atom *qAt = d_query.getAtomWithIdx(i);
//         if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
//             qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
//           const Atom *mAt = d_mol.getAtomWithIdx(j);
//           if (!d_params.specifiedStereoQueryMatchesUnspecified &&
//               mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
//               mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) {
//             return false;
//           }
//         }
//       }
//       res = atomCompat(d_query[i], d_mol[j], d_params);
//       return res;
//     }
//    private:
//     const ROMol &d_query;
//     const ROMol &d_mol;
//     const SubstructMatchParameters &d_params;
//   };
//
// RDKit❗✔️: AtomLabelFunctor is ported as plain functions. The
//   useChirality specified/unspecified precheck is wired below; the final
//   tetrahedral parity check remains in MolMatchFinalCheckFunctor.

fn has_chiral_label(chiral_tag: ChiralTag) -> bool {
    // RDKit✔️✔️: bool hasChiralLabel(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   return at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:          at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️: }
    // Rust's reference type enforces the non-null precondition. Complexity
    // review: both implementations read one enum and perform at most two O(1)
    // comparisons without allocation.
    matches!(
        chiral_tag,
        ChiralTag::TetrahedralCw | ChiralTag::TetrahedralCcw
    )
}

type MatchVect = Vec<(i32, i32)>;

fn insert_if_needed(matches: &mut BTreeSet<MatchVect>, candidate: MatchVect) -> bool {
    // RDKit✔️✔️: bool insertIfNeeded(std::set<MatchVectType> &matches, const MatchVectType &m) {
    // RDKit✔️✔️:   bool shouldInsert = true;
    // RDKit✔️✔️:   std::unordered_set<int> matchAsSet;
    // RDKit✔️✔️:   std::transform(m.begin(), m.end(),
    // RDKit✔️✔️:                  std::inserter(matchAsSet, matchAsSet.begin()),
    // RDKit✔️✔️:                  [](const std::pair<int, int> &p) { return p.second; });
    // RDKit✔️✔️:   for (auto it = matches.begin(); it != matches.end(); ++it) {
    // RDKit✔️✔️:     std::unordered_set<int> existingMatchAsSet;
    // RDKit✔️✔️:     std::transform(
    // RDKit✔️✔️:         it->begin(), it->end(),
    // RDKit✔️✔️:         std::inserter(existingMatchAsSet, existingMatchAsSet.begin()),
    // RDKit✔️✔️:         [](const std::pair<int, int> &p) { return p.second; });
    // RDKit✔️✔️:     if (matchAsSet == existingMatchAsSet) {
    // RDKit✔️✔️:       if (m < *it) {
    // RDKit✔️✔️:         matches.erase(it);
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         shouldInsert = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (shouldInsert) {
    // RDKit✔️✔️:     matches.insert(m);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return shouldInsert;
    // RDKit✔️✔️: }
    // Complexity review: both scan O(number of matches), build one O(match
    // length) hash set per comparison, and use a logarithmic ordered-set erase
    // and insert. Rust retains no temporary sets after the call.
    let candidate_atoms: HashSet<i32> = candidate.iter().map(|pair| pair.1).collect();
    let existing = matches.iter().find(|existing| {
        existing.iter().map(|pair| pair.1).collect::<HashSet<_>>() == candidate_atoms
    });
    let mut should_insert = true;
    if let Some(existing) = existing.cloned() {
        if candidate < existing {
            matches.remove(&existing);
        } else {
            should_insert = false;
        }
    }
    if should_insert {
        matches.insert(candidate);
    }
    should_insert
}

fn try_to_insert(
    matches: &mut BTreeSet<MatchVect>,
    candidate: MatchVect,
    params: &SubstructMatchParams,
) -> bool {
    // RDKit✔️✔️: bool tryToInsert(std::set<MatchVectType> &matches, const MatchVectType &match,
    // RDKit✔️✔️:                  const SubstructMatchParameters &params) {
    // RDKit✔️✔️:   if (matches.size() == params.maxMatches) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!params.uniquify) {
    // RDKit✔️✔️:     matches.insert(match);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     insertIfNeeded(matches, match);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Complexity review: the limit check is O(1), ordinary insertion is
    // O(log M), and the uniquify branch delegates to the source-equivalent
    // O(M * match length) canonical helper without additional copying.
    if matches.len() == params.max_matches {
        return false;
    }
    if !params.uniquify {
        matches.insert(candidate);
    } else {
        insert_if_needed(matches, candidate);
    }
    true
}

fn atom_label_matches(
    query: &QueryGraph,
    mol: &SearchTarget<'_>,
    query_index: usize,
    mol_index: usize,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> Result<bool, PeriodicTableError> {
    // RDKit✔️✔️: bool operator()(unsigned int i, unsigned int j) const {
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:     if (d_params.useChirality) {
    // RDKit✔️✔️:       const Atom *qAt = d_query.getAtomWithIdx(i);
    // RDKit✔️✔️:       if (qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:           qAt->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:         const Atom *mAt = d_mol.getAtomWithIdx(j);
    // RDKit✔️✔️:         if (!d_params.specifiedStereoQueryMatchesUnspecified &&
    // RDKit✔️✔️:             mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW &&
    // RDKit✔️✔️:             mAt->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW) {
    // RDKit✔️✔️:           return false;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   res = atomCompat(d_query[i], d_mol[j], d_params);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: the precheck is O(1), then this delegates exactly once
    // to canonical atom_compat; it introduces no allocation or repeated query
    // evaluation beyond the source functor.
    let query_atom = &query.atoms()[query_index];
    let mol_atom = &mol.atoms()[mol_index];
    if params.use_chirality
        && has_chiral_label(query_atom.chiral_tag())
        && !params.specified_stereo_query_matches_unspecified
        && !has_chiral_label(mol_atom.chiral_tag())
    {
        return Ok(false);
    }
    atom_compat(
        query_atom,
        query,
        mol_atom,
        mol,
        params,
        recursive_cache,
        query_ctx,
    )
}

fn atom_matches(query_atom: &QueryAtom, mol_atom: &Atom, mol: &SearchTarget<'_>) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Atom.cpp :: Atom::Match
    // RDKit✔️✔️: bool Atom::Match(Atom const *what) const {
    // RDKit✔️✔️:   PRECONDITION(what, "bad query atom");
    // RDKit✔️✔️:   bool res = getAtomicNum() == what->getAtomicNum();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // special dummy--dummy match case:
    // RDKit✔️✔️:   //   [*] matches [*],[1*],[2*],etc.
    // RDKit✔️✔️:   //   [1*] only matches [*] and [1*]
    // RDKit✔️✔️:   if (res) {
    // RDKit✔️✔️:     if (!this->getAtomicNum()) {
    // RDKit✔️✔️:       // this is the new behavior, based on the isotopes:
    // RDKit✔️✔️:       int tgt = this->getIsotope();
    // RDKit✔️✔️:       int test = what->getIsotope();
    // RDKit✔️✔️:       if (tgt && test && tgt != test) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       // standard atom-atom match: The general rule here is that if this atom
    // RDKit✔️✔️:       // has a property that
    // RDKit✔️✔️:       // deviates from the default, then the other atom should match that value.
    // RDKit✔️✔️:       if ((this->getFormalCharge() &&
    // RDKit✔️✔️:            this->getFormalCharge() != what->getFormalCharge()) ||
    // RDKit✔️✔️:           (this->getIsotope() && this->getIsotope() != what->getIsotope()) ||
    // RDKit✔️✔️:           (this->getNumRadicalElectrons() &&
    // RDKit✔️✔️:            this->getNumRadicalElectrons() != what->getNumRadicalElectrons())) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Behavior: `None` and an explicit zero isotope both represent RDKit's
    // `getIsotope() == 0`. The target atomic-number accessor also preserves the
    // effective source identity used by detached matching overrides.
    // Complexity: the source performs constant-time scalar reads. This path
    // adds one O(1) target-identity access (an optional indexed override read),
    // with no allocation, cloning, molecule scan, keyed lookup, or collection.
    let query_atomic_number = query_atom.atomic_number();
    let target_atomic_number = mol.query_atomic_number(mol_atom);
    if query_atomic_number != target_atomic_number {
        return false;
    }
    let query_isotope = query_atom.isotope().unwrap_or(0);
    let target_isotope = mol_atom.isotope().unwrap_or(0);
    if query_atomic_number == 0 {
        return query_isotope == 0 || target_isotope == 0 || query_isotope == target_isotope;
    }
    (query_atom.formal_charge() == 0 || query_atom.formal_charge() == mol_atom.formal_charge())
        && (query_isotope == 0 || query_isotope == target_isotope)
        && (query_atom.radical_electrons() == 0
            || query_atom.radical_electrons() == mol_atom.radical_electrons())
}

fn recursive_smarts_root_matches(
    atom: &Atom,
    recursive_query: &crate::query_behavior::RecursiveStructureQuery,
    mol: &SearchTarget<'_>,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.h :: RecursiveStructureQuery
    // RDKit✔️✔️: class RDKIT_GRAPHMOL_EXPORT RecursiveStructureQuery
    // RDKit✔️✔️:     : public Queries::SetQuery<int, Atom const *, true> {
    // RDKit✔️✔️:   RecursiveStructureQuery(ROMol const *query, unsigned int serialNumber = 0)
    // RDKit✔️✔️:       : Queries::SetQuery<int, Atom const *, true>(),
    // RDKit✔️✔️:         d_serialNumber(serialNumber) {
    // RDKit✔️✔️:     setQueryMol(query);
    // RDKit✔️✔️:     setDataFunc(getAtIdx);
    // RDKit✔️✔️:     setDescription("RecursiveStructure");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   static inline int getAtIdx(Atom const *at) {
    // RDKit✔️✔️:     PRECONDITION(at, "bad atom argument");
    // RDKit✔️✔️:     return at->getIdx();
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp :: detail::RecursiveMatcher
    // RDKit✔️✔️:   if (!query.hasProp(common_properties::_queryRootAtom)) {
    // RDKit✔️✔️:     matches.push_back(pairs.begin()->second);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     int rootIdx;
    // RDKit✔️✔️:     query.getProp(common_properties::_queryRootAtom, rootIdx);
    // RDKit✔️✔️:     bool found = false;
    // RDKit✔️✔️:     for (const auto &pairIter : pairs) {
    // RDKit✔️✔️:       if (pairIter.first == static_cast<unsigned int>(rootIdx)) {
    // RDKit✔️✔️:         matches.push_back(pairIter.second);
    // RDKit✔️✔️:         found = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // END RDKIT CPP FUNCTION
    //
    // COSMolKit currently parses the recursive SMARTS used by Lipinski NumHBA
    // without `_queryRootAtom`; matching therefore uses RDKit's first mapped
    // query atom as the recursive root and tests membership in the cached
    // RecursiveStructureQuery atom-index set.
    if let Some(cache) = recursive_cache {
        return cache
            .get(&recursive_query_cache_key(recursive_query))
            .and_then(|match_starts| match_starts.get(atom.id().index()))
            .copied()
            .unwrap_or(false);
    }

    let Some(query) = recursive_query.query_graph() else {
        return false;
    };
    substruct_match_impl(
        mol,
        query,
        &SubstructMatchParams {
            max_matches: 1000,
            uniquify: false,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
            ..Default::default()
        },
    )
    .unwrap_or_default()
    .into_iter()
    .any(|matched| matched.atom_mapping.first().copied() == Some(atom.id().index()))
}

fn atom_query_predicate_matches_for_substruct(
    atom: &Atom,
    pred: &AtomQueryPredicate,
    mol: &SearchTarget<'_>,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> Result<bool, PeriodicTableError> {
    match pred {
        // RDKit✔️✔️: Chiral SMARTS labels are not ordinary atom-compatibility
        // constraints when `useChirality` is false. AtomLabelFunctor and
        // MolMatchFinalCheckFunctor handle stereochemistry explicitly.
        AtomQueryPredicate::ChiralTagMatch(_) | AtomQueryPredicate::ChiralPermutationMatch(_)
            if !params.use_chirality =>
        {
            Ok(true)
        }
        AtomQueryPredicate::RecursiveSmarts(recursive_query) => Ok(recursive_smarts_root_matches(
            atom,
            recursive_query,
            mol,
            recursive_cache,
        )),
        _ => atom_predicate_matches_with_context(atom, pred, mol, query_ctx),
    }
}

/// RDKit❗✔️: Evaluation of an atom query node for the SMARTS subset currently
/// modeled by COSMolKit.
///
/// Recursive SMARTS are evaluated through the recursive match cache used by
/// SubstructMatch; unsupported predicate leaves still evaluate false.
fn evaluate_atom_query(
    query: &crate::QueryNode<AtomQueryPredicate>,
    atom: &Atom,
    mol: &SearchTarget<'_>,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> Result<bool, PeriodicTableError> {
    match query {
        crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
            atom,
            pred,
            mol,
            params,
            recursive_cache,
            query_ctx,
        ),
        crate::QueryNode::And(children) => {
            for child in children {
                if !evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)? {
                    return Ok(false);
                }
            }
            Ok(true)
        }
        crate::QueryNode::Or(children) => {
            for child in children {
                if evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)? {
                    return Ok(true);
                }
            }
            Ok(false)
        }
        crate::QueryNode::Xor(children) => {
            let mut matched = false;
            for child in children {
                if evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)? {
                    if matched {
                        return Ok(false);
                    }
                    matched = true;
                }
            }
            Ok(matched)
        }
        crate::QueryNode::Not(child) => Ok(!evaluate_atom_query(
            child,
            atom,
            mol,
            params,
            recursive_cache,
            query_ctx,
        )?),
    }
}

// RDKit source (SubstructMatch.cpp):
//   class BondLabelFunctor {
//    public:
//     BondLabelFunctor(const ROMol &query, const ROMol &mol,
//                      const SubstructMatchParameters &ps)
//         : d_query(query), d_mol(mol), d_params(ps) {};
//     bool operator()(MolGraph::edge_descriptor i,
//                     MolGraph::edge_descriptor j) const {
//       if (d_params.useChirality) {
//         const Bond *qBnd = d_query[i];
//         if (qBnd->getBondType() == Bond::DOUBLE &&
//             qBnd->getStereo() > Bond::STEREOANY) {
//           const Bond *mBnd = d_mol[j];
//           if (mBnd->getBondType() == Bond::DOUBLE &&
//               !d_params.specifiedStereoQueryMatchesUnspecified &&
//               mBnd->getStereo() <= Bond::STEREOANY) {
//             return false;
//           }
//         }
//       }
//       bool res = bondCompat(d_query[i], d_mol[j], d_params);
//       return res;
//     }
//    private:
//     const ROMol &d_query;
//     const ROMol &d_mol;
//     const SubstructMatchParameters &d_params;
//   };

fn rdkit_bond_stereo_is_above_any(stereo: BondStereo) -> bool {
    !matches!(stereo, BondStereo::None | BondStereo::Any)
}

fn bond_label_matches(
    query: &QueryGraph,
    mol: &SearchTarget<'_>,
    query_index: usize,
    mol_index: usize,
    params: &SubstructMatchParams,
    query_ctx: &QueryMatchContext,
) -> bool {
    // RDKit✔️✔️: bool operator()(MolGraph::edge_descriptor i,
    // RDKit✔️✔️:                 MolGraph::edge_descriptor j) const {
    // RDKit✔️✔️:   if (d_params.useChirality) {
    // RDKit✔️✔️:     const Bond *qBnd = d_query[i];
    // RDKit✔️✔️:     if (qBnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:         qBnd->getStereo() > Bond::STEREOANY) {
    // RDKit✔️✔️:       const Bond *mBnd = d_mol[j];
    // RDKit✔️✔️:       if (mBnd->getBondType() == Bond::DOUBLE &&
    // RDKit✔️✔️:           !d_params.specifiedStereoQueryMatchesUnspecified &&
    // RDKit✔️✔️:           mBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool res = bondCompat(d_query[i], d_mol[j], d_params);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Complexity review: the stereo precheck is O(1), then this delegates
    // exactly once to canonical bond_compat. It adds no allocation, scan,
    // cloning, or repeated query evaluation beyond the source functor.
    let query_bond = &query.bonds()[query_index];
    let mol_bond = &mol.bonds()[mol_index];
    if params.use_chirality
        && query_bond.bond().order() == BondOrder::Double
        && rdkit_bond_stereo_is_above_any(query_bond.bond().stereo())
        && mol_bond.order() == BondOrder::Double
        && !params.specified_stereo_query_matches_unspecified
        && !rdkit_bond_stereo_is_above_any(mol_bond.stereo())
    {
        return false;
    }
    bond_compat(query_bond, query, mol_bond, mol, params, query_ctx)
}

/// RDKit❗✔️: Evaluation of a bond query node for the currently modeled SMARTS
/// bond predicate subset.
fn evaluate_bond_query(
    query: &crate::QueryNode<BondQueryPredicate>,
    bond: &Bond,
    mol: &SearchTarget<'_>,
    query_ctx: &QueryMatchContext,
) -> bool {
    match query {
        crate::QueryNode::Predicate(pred) => {
            bond_predicate_matches_with_context(bond, pred, mol, query_ctx)
        }
        crate::QueryNode::And(children) => and_query_match(children, false, |child| {
            evaluate_bond_query(child, bond, mol, query_ctx)
        }),
        crate::QueryNode::Or(children) => or_query_match(children, false, |child| {
            evaluate_bond_query(child, bond, mol, query_ctx)
        }),
        crate::QueryNode::Xor(children) => xor_query_match(children, false, |child| {
            evaluate_bond_query(child, bond, mol, query_ctx)
        }),
        crate::QueryNode::Not(child) => !evaluate_bond_query(child, bond, mol, query_ctx),
    }
}

// ---------------------------------------------------------------------------
// VF2 State Machine
// ---------------------------------------------------------------------------
//
// ## RDKit source reproduction: vf2.hpp
//
// The following section reproduces the VF2SubState class from vf2.hpp.
// The C++ code is shown as verbatim comments with RDKit markers.
//
// ### Key design differences from RDKit:
//
// 1. `core_1`/`core_2`: Same role — mapping from query atom idx → mol atom idx
//    and vice versa. Both use the `NULL_NODE` sentinel, matching vf2.hpp.
//
// 2. `term_1`/`term_2`: Stores the core_len *depth* at which each atom was
//    added to the terminal set, exactly as in vf2.hpp. BackTrack decrements
//    counters keyed by depth, not recomputes from scratch.
//
// 3. RDKit supports COW state copies through `share_count`, but its canonical
//    `vf2`/`vf2_all` entry constructs one state and recursively mutates and
//    backtracks that same state. Rust matching follows that same single-state
//    traversal. The explicit `clone_state`/`Clone` path deep-copies its vectors
//    for an independent caller; canonical recursion does not call it.
//
// 4. `Vf2GraphRef` borrows QueryGraph and TopologyBlock adjacency directly.
//    Only an explicit compiled-query plan owns a materialized `Vf2Graph`.

// RDKit source (vf2.hpp):
//   typedef std::uint32_t node_id;
//   const node_id NULL_NODE = 0xFFFFFFFF;

type NodeId = usize;
const NULL_NODE: NodeId = usize::MAX;

// RDKit source (vf2.hpp):
//   template <class Graph>
//   struct Pair {
//     node_id n1, n2;
//     bool hasiter{false};
//     RDK_ADJ_ITER nbrbeg, nbrend;
//     Pair() : n1(NULL_NODE), n2(NULL_NODE) {}
//   };

#[derive(Debug, Clone)]
struct Vf2Pair {
    n1: NodeId,
    n2: NodeId,
    hasiter: bool,
    /// VF2+ source atom in the mol graph (g2) whose adjacency drives the
    /// neighbor iterator.
    nbr_node: NodeId,
    /// VF2+ neighbor iterator over mol graph (g2) neighbors.
    nbr_cursor: usize,
    nbr_end: usize,
}

impl Vf2Pair {
    fn new() -> Self {
        Self {
            n1: NULL_NODE,
            n2: NULL_NODE,
            hasiter: false,
            nbr_node: NULL_NODE,
            nbr_cursor: 0,
            nbr_end: 0,
        }
    }
}

#[derive(Debug, Clone, Copy)]
struct NodeInfo {
    id: u32,
    in_deg: u32,
    out_deg: u32,
}

fn node_info_cmp1(a: &NodeInfo, b: &NodeInfo) -> std::cmp::Ordering {
    // RDKit✔️✔️: static bool nodeInfoComp1(const NodeInfo &a, const NodeInfo &b) {
    // RDKit✔️✔️:   if (a.out < b.out) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.out > b.out) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.in < b.in) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.in > b.in) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: return false;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform at most two integer
    // comparisons in O(1) time without allocation or temporary collections.
    a.out_deg
        .cmp(&b.out_deg)
        .then_with(|| a.in_deg.cmp(&b.in_deg))
}

fn node_info_cmp2(a: &NodeInfo, b: &NodeInfo) -> std::cmp::Ordering {
    // RDKit✔️✔️: static int nodeInfoComp2(const NodeInfo &a, const NodeInfo &b) {
    // RDKit✔️✔️:   if (!a.in && b.in) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.in && !b.in) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.out < b.out) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.out > b.out) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.in < b.in) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a.in > b.in) {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: }
    // Complexity review: both implementations perform a bounded sequence of
    // integer comparisons in O(1) time without allocation or cloning.
    if a.in_deg == 0 && b.in_deg != 0 {
        return std::cmp::Ordering::Greater;
    }
    if a.in_deg != 0 && b.in_deg == 0 {
        return std::cmp::Ordering::Less;
    }
    a.out_deg
        .cmp(&b.out_deg)
        .then_with(|| a.in_deg.cmp(&b.in_deg))
}

// RDKit source (vf2.hpp), SortNodesByFrequency:
//   Sorts the nodes of a graphs, returning a heap-allocated vector
//   with the node ids in the proper orders.
//   The sorting criterion takes into account:
//     1 - The number of nodes with the same in/out degree.
//     2 - The valence of the nodes.
//   The nodes at the beginning of the vector are the most singular,
//   from which the matching should start.

fn sort_nodes_by_frequency(g: Vf2GraphRef<'_>) -> Vec<NodeId> {
    // RDKit✔️✔️: template <class Graph>
    // RDKit✔️✔️: node_id *SortNodesByFrequency(const Graph *g) {
    // RDKit✔️✔️:   std::vector<NodeInfo> vect;
    // RDKit✔️✔️:   vect.reserve(boost::num_vertices(*g));
    // RDKit✔️✔️:   typename Graph::vertex_iterator bNode, eNode;
    // RDKit✔️✔️:   boost::tie(bNode, eNode) = boost::vertices(*g);
    // RDKit✔️✔️:   while (bNode != eNode) {
    // RDKit✔️✔️:     NodeInfo t;
    // RDKit✔️✔️:     t.id = vect.size();
    // RDKit✔️✔️:     t.in = boost::out_degree(*bNode, *g);  // <- assuming undirected graph
    // RDKit✔️✔️:     t.out = boost::out_degree(*bNode, *g);
    // RDKit✔️✔️:     vect.push_back(t);
    // RDKit✔️✔️:     ++bNode;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::sort(vect.begin(), vect.end(), nodeInfoComp1);
    let mut vect: Vec<NodeInfo> = (0..g.num_atoms())
        .map(|i| {
            // RDKit's NodeInfo uses node_id (uint32_t) for all three fields.
            // The detached graph uses usize indices, so convert at this
            // source-width metadata boundary and widen IDs when returning.
            let deg = g.out_degree(i) as u32;
            NodeInfo {
                id: i as u32,
                in_deg: deg,
                out_deg: deg,
            }
        })
        .collect();
    vect.sort_unstable_by(node_info_cmp1);

    // RDKit✔️✔️:   unsigned int run = 1;
    // RDKit✔️✔️:   for (unsigned int i = 0; i < vect.size(); i += run) {
    // RDKit✔️✔️:     for (run = 1; i + run < vect.size() && vect[i + run].in == vect[i].in &&
    // RDKit✔️✔️:                   vect[i + run].out == vect[i].out;
    // RDKit✔️✔️:          ++run) {
    // RDKit✔️✔️:       ;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     for (unsigned int j = 0; j < run; ++j) {
    // RDKit✔️✔️:       vect[i + j].in += vect[i + j].out;
    // RDKit✔️✔️:       vect[i + j].out = run;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut i = 0;
    while i < vect.len() {
        let mut run = 1;
        while i + run < vect.len()
            && vect[i + run].in_deg == vect[i].in_deg
            && vect[i + run].out_deg == vect[i].out_deg
        {
            run += 1;
        }
        for j in 0..run {
            // `NodeInfo::in` is uint32_t upstream, so unsigned overflow wraps
            // even in debug builds where ordinary Rust addition would panic.
            vect[i + j].in_deg = vect[i + j].in_deg.wrapping_add(vect[i + j].out_deg); // valence sum
            vect[i + j].out_deg = run as u32; // frequency
        }
        i += run;
    }

    // RDKit✔️✔️:   std::sort(vect.begin(), vect.end(), nodeInfoComp2);
    // The source comparator has no node-ID tiebreak; equal records remain
    // unordered just as with std::sort, so do not add a stable ID key here.
    vect.sort_unstable_by(node_info_cmp2);

    // RDKit✔️✔️:   node_id *nodes = new node_id[vect.size()];
    // RDKit✔️✔️:   for (unsigned int i = 0; i < vect.size(); ++i) {
    // RDKit✔️✔️:     nodes[i] = vect[i].id;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return nodes;
    // RDKit✔️✔️: }
    // Complexity review: both versions allocate O(V) node metadata and an
    // O(V) result, perform two O(V log V) unstable sorts, and scan runs in
    // O(V). Degree lookup and all loop bodies remain O(1) per visited node.
    vect.iter().map(|ni| ni.id as usize).collect()
}

// RDKit source (vf2.hpp), VF2SubState class:
//   template <class Graph, class VertexCompatible, class EdgeCompatible,
//             class MatchChecking>
//   class VF2SubState {
//    private:
//     Graph *g1, *g2;
//     VertexCompatible &vc;
//     EdgeCompatible &ec;
//     MatchChecking &mc;
//     unsigned int n1, n2;
//     unsigned int core_len;
//     unsigned int t1_len;
//     unsigned int t2_len;  // Core nodes are also counted by these...
//     node_id *core_1;
//     node_id *core_2;
//     node_id *term_1;
//     node_id *term_2;
//     node_id *order;
//     long *share_count;
//     int *vs_compared;

/// RDKit❗✔️: VF2 subgraph isomorphism state.
///
/// g1 = query graph, g2 = molecule graph.
/// core_1[i] = mapping from query atom i -> mol atom j (or NULL_NODE).
/// core_2[j] = mapping from mol atom j -> query atom i (or NULL_NODE).
/// term_1[i] = depth (core_len) when atom i entered terminal set (0 = not terminal).
/// term_2[j] = same for mol atoms.
struct Vf2SubState<'a> {
    g1: Vf2GraphRef<'a>,
    g2: Vf2GraphRef<'a>,
    n1: usize,
    n2: usize,
    core_len: usize,
    t1_len: usize,
    t2_len: usize,
    core_1: Vec<NodeId>,
    core_2: Vec<NodeId>,
    term_1: Vec<usize>,
    term_2: Vec<usize>,
    order: Option<Vec<NodeId>>,
}

impl<'a> Vf2SubState<'a> {
    fn new(g1: Vf2GraphRef<'a>, g2: Vf2GraphRef<'a>, sort_nodes: bool) -> Self {
        // RDKit✔️✔️: VF2SubState(Graph *ag1, Graph *ag2, VertexCompatible &avc,
        // RDKit✔️✔️:             EdgeCompatible &aec, MatchChecking &amc, bool sortNodes = false)
        // RDKit✔️✔️:     : g1(ag1),
        // RDKit✔️✔️:       g2(ag2),
        // RDKit✔️✔️:       vc(avc),
        // RDKit✔️✔️:       ec(aec),
        // RDKit✔️✔️:       mc(amc),
        // RDKit✔️✔️:       n1(num_vertices(*ag1)),
        // RDKit✔️✔️:       n2(num_vertices(*ag2)) {
        // RDKit✔️✔️:   if (sortNodes) {
        // RDKit✔️✔️:     order = SortNodesByFrequency(ag1);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     order = nullptr;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   core_len = 0;
        // RDKit✔️✔️:   t1_len = 0;
        // RDKit✔️✔️:   t2_len = 0;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   core_1 = new node_id[n1];
        // RDKit✔️✔️:   core_2 = new node_id[n2];
        // RDKit✔️✔️:   term_1 = new node_id[n1];
        // RDKit✔️✔️:   term_2 = new node_id[n2];
        // RDKit✔️✔️:   share_count = new long;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   for (unsigned int i = 0; i < n1; i++) {
        // RDKit✔️✔️:     core_1[i] = NULL_NODE;
        // RDKit✔️✔️:     term_1[i] = 0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   for (unsigned int i = 0; i < n2; i++) {
        // RDKit✔️✔️:     core_2[i] = NULL_NODE;
        // RDKit✔️✔️:     term_2[i] = 0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   vs_compared = nullptr;
        // RDKit✔️✔️:   // vs_compared = new int[n1*n2];
        // RDKit✔️✔️:   // memset((void *)vs_compared,0,n1*n2*sizeof(int));
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // es_compared = new std::map<unsigned int,bool>();
        // RDKit✔️✔️:   *share_count = 1;
        // RDKit✔️✔️: }
        // The compatibility functors remain explicit arguments to Rust match
        // methods, so the state stores only the source fields those methods use.
        // Complexity review: both implementations initialize four O(V) arrays
        // and optionally run the same O(V log V) ordering routine. Vec uses the
        // same contiguous storage and does not add asymptotic or hot-path work.
        let n1 = g1.num_atoms();
        let n2 = g2.num_atoms();
        let order = if sort_nodes {
            Some(sort_nodes_by_frequency(g1))
        } else {
            None
        };

        // RDKit✔️✔️: core_len = 0; t1_len = 0; t2_len = 0;
        // RDKit✔️✔️: core_1[i] = NULL_NODE; term_1[i] = 0;
        // RDKit✔️✔️: core_2[j] = NULL_NODE; term_2[j] = 0;
        Self {
            g1,
            g2,
            n1,
            n2,
            core_len: 0,
            t1_len: 0,
            t2_len: 0,
            core_1: vec![NULL_NODE; n1],
            core_2: vec![NULL_NODE; n2],
            term_1: vec![0usize; n1],
            term_2: vec![0usize; n2],
            order,
        }
    }

    fn with_order(g1: Vf2GraphRef<'a>, g2: Vf2GraphRef<'a>, order: &[usize]) -> Self {
        let mut state = Self::new(g1, g2, false);
        state.order = Some(order.to_vec());
        state
    }

    fn clone_state(&self) -> Self {
        // RDKit✔️❌: VF2SubState(const VF2SubState &state)
        // RDKit✔️❌:     : g1(state.g1),
        // RDKit✔️❌:       g2(state.g2),
        // RDKit✔️❌:       vc(state.vc),
        // RDKit✔️❌:       ec(state.ec),
        // RDKit✔️❌:       mc(state.mc),
        // RDKit✔️❌:       n1(state.n1),
        // RDKit✔️❌:       n2(state.n2),
        // RDKit✔️❌:       order(state.order),
        // RDKit✔️❌:       vs_compared(state.vs_compared)
        // RDKit✔️❌:   // es_compared(state.es_compared)
        // RDKit✔️❌: {
        // RDKit✔️❌:   core_len = state.core_len;
        // RDKit✔️❌:   t1_len = state.t1_len;
        // RDKit✔️❌:   t2_len = state.t2_len;
        // RDKit✔️❌:
        // RDKit✔️❌:   core_1 = state.core_1;
        // RDKit✔️❌:   core_2 = state.core_2;
        // RDKit✔️❌:   term_1 = state.term_1;
        // RDKit✔️❌:   term_2 = state.term_2;
        // RDKit✔️❌:   share_count = state.share_count;
        // RDKit✔️❌:
        // RDKit✔️❌:   ++(*share_count);
        // RDKit✔️❌: }
        // Compatibility callbacks are passed to Rust match calls rather than
        // stored in the state. Deep-copying Vec state preserves the copied
        // values and makes subsequent mutation independent. Complexity review:
        // this is O(V) with five allocations, while RDKit shares the arrays and
        // increments one reference count in O(1).
        Self {
            g1: self.g1,
            g2: self.g2,
            n1: self.n1,
            n2: self.n2,
            core_len: self.core_len,
            t1_len: self.t1_len,
            t2_len: self.t2_len,
            core_1: self.core_1.clone(),
            core_2: self.core_2.clone(),
            term_1: self.term_1.clone(),
            term_2: self.term_2.clone(),
            order: self.order.clone(),
        }
    }

    fn clone(&self) -> Self {
        // RDKit✔️❌: VF2SubState *Clone() { return new VF2SubState(*this); }
        // Complexity review: this forwards to the single O(V) Rust state-copy
        // implementation, while RDKit's shared-array copy is O(1). No second
        // clone path is introduced.
        self.clone_state()
    }

    fn debug_order(&self) -> Option<&[NodeId]> {
        self.order.as_deref()
    }

    fn is_goal(&self) -> bool {
        // RDKit✔️✔️: bool IsGoal() { return core_len == n1; }
        // Complexity review: one integer equality in O(1), without allocation.
        self.core_len == self.n1
    }

    fn match_checks(
        &self,
        c1: &[NodeId],
        c2: &[NodeId],
        check: &mut impl FnMut(&[NodeId], &[NodeId]) -> bool,
    ) -> bool {
        // RDKit✔️✔️: bool MatchChecks(const node_id c1[], const node_id c2[]) {
        // RDKit✔️✔️:   return mc(c1, c2);
        // RDKit✔️✔️: }
        // Complexity review: both forms make one callback invocation and pass
        // existing mapping storage by reference without allocation or cloning.
        check(c1, c2)
    }

    fn is_dead(&self) -> bool {
        // RDKit✔️✔️: bool IsDead() { return n1 > n2 || t1_len > t2_len; }
        // Complexity review: at most two integer comparisons in O(1), without
        // allocation or temporary collections.
        self.n1 > self.n2 || self.t1_len > self.t2_len
    }

    fn core_len(&self) -> usize {
        // RDKit✔️✔️: unsigned int CoreLen() { return core_len; }
        // Complexity review: one field read in O(1), without allocation.
        self.core_len
    }

    // RDKit source (vf2.hpp):
    //   bool NextPair(Pair<Graph> &pair) {
    //     if (pair.n1 == NULL_NODE) { pair.n1 = 0; }
    //     if (pair.n2 == NULL_NODE) { pair.n2 = 0; }
    //     else { pair.n2++; }
    //     ...
    //     if (t1_len > core_len && t2_len > core_len) {
    //       while (pair.n1 < n1 &&
    //              (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
    //         pair.n1++; pair.n2 = 0;
    //       }
    //       ...
    //     } else if (pair.n1 == 0 && order != nullptr) {
    //       // Optimisation: ...
    //       unsigned int i = 0;
    //       while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) { i++; }
    //       ...
    //     } else {
    //       while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
    //         pair.n1++; pair.n2 = 0;
    //       }
    //     }
    //     // VF2 Plus iterator ...
    //     if (pair.hasiter) { ... }
    //     else if (t1_len > core_len && t2_len > core_len) {
    //       while (pair.n2 < n2 &&
    //              (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
    //         pair.n2++;
    //       }
    //     } else {
    //       while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) { pair.n2++; }
    //     }
    //     return pair.n1 < n1 && pair.n2 < n2;
    //   }

    /// RDKit✔️❌: NextPair — find the next candidate pair (n1 from query,
    ///   n2 from mol) to try matching.
    ///
    /// Uses terminal-set-based iteration from vf2.hpp, including the VF2+
    /// neighbor iterator that restricts mol-side candidates to neighbors of
    /// the already-mapped terminal predecessor.
    fn next_pair(&self, pair: &mut Vf2Pair) -> bool {
        // RDKit✔️✔️: bool NextPair(Pair<Graph> &pair) {
        // RDKit✔️✔️:   if (pair.n1 == NULL_NODE) {
        // RDKit✔️✔️:     pair.n1 = 0;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (pair.n2 == NULL_NODE) {
        // RDKit✔️✔️:     pair.n2 = 0;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     pair.n2++;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️: #if 0
        // RDKit✔️✔️:   std::cerr<<" **** np: "<< prev_n1<<","<<prev_n2<<std::endl;
        // RDKit✔️✔️:   std::cerr<<"in_1 ";
        // RDKit✔️✔️:   for(unsigned int i=0;i<n1;++i){
        // RDKit✔️✔️:     std::cerr<<"("<<in_1[i]<<","<<out_1[i]<<"), ";
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::cerr<<std::endl;
        // RDKit✔️✔️:   std::cerr<<"in_2 ";
        // RDKit✔️✔️:   for(unsigned int i=0;i<n2;++i){
        // RDKit✔️✔️:     std::cerr<<"("<<in_2[i]<<","<<out_2[i]<<"), ";
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::cerr<<std::endl;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:   if (t1_len > core_len && t2_len > core_len) {
        // RDKit✔️✔️:     while (pair.n1 < n1 &&
        // RDKit✔️✔️:            (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
        // RDKit✔️✔️:       pair.n1++;
        // RDKit✔️✔️:       pair.n2 = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     /* Initialize VF2 Plus neighbor iterator.
        // RDKit✔️✔️:      * The next query node (pair.n1) has been selected from the terminal
        // RDKit✔️✔️:      * set and is therefore adjacent to an already mapped atom (in
        // RDKit✔️✔️:      * core_1). Rather than select pair.n2 from all atoms (0...n2) we can
        // RDKit✔️✔️:      * select it from the neighbors of this mapped atom (0...deg(nbor))
        // RDKit✔️✔️:      * since it must also be adajcent to this mapped atom!
        // RDKit✔️✔️:      */
        // RDKit✔️✔️:     if (!pair.hasiter) {
        // RDKit✔️✔️:       RDK_ADJ_ITER n1iter_beg, n1iter_end;
        // RDKit✔️✔️:       boost::tie(n1iter_beg, n1iter_end) =
        // RDKit✔️✔️:           boost::adjacent_vertices(pair.n1, *g1);
        // RDKit✔️✔️:
        // RDKit✔️✔️:       while (n1iter_beg != n1iter_end && core_1[*n1iter_beg] == NULL_NODE) {
        // RDKit✔️✔️:         ++n1iter_beg;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:
        // RDKit✔️✔️:       assert(n1iter_beg != n1iter_end);
        // RDKit✔️✔️:
        // RDKit✔️✔️:       boost::tie(pair.nbrbeg, pair.nbrend) =
        // RDKit✔️✔️:           boost::adjacent_vertices(core_1[*n1iter_beg], *g2);
        // RDKit✔️✔️:       pair.hasiter = true;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (pair.n1 == 0 && order != nullptr) {
        // RDKit✔️✔️:     // Optimisation: if the order vector is laid out in a DFS/BFS then this
        // RDKit✔️✔️:     // loop can be replaced with:
        // RDKit✔️✔️:     //   pair.n1=order[core_len];
        // RDKit✔️✔️:     // :)
        // RDKit✔️✔️:     unsigned int i = 0;
        // RDKit✔️✔️:     while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) {
        // RDKit✔️✔️:       i++;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (i == n1) {
        // RDKit✔️✔️:       pair.n1 = n1;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
        // RDKit✔️✔️:       pair.n1++;
        // RDKit✔️✔️:       pair.n2 = 0;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   /* VF2 Plus iterator available? */
        // RDKit✔️✔️:   if (pair.hasiter) {
        // RDKit✔️✔️:     while (pair.nbrbeg < pair.nbrend && core_2[*pair.nbrbeg] != NULL_NODE) {
        // RDKit✔️✔️:       ++pair.nbrbeg;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     if (pair.nbrbeg < pair.nbrend) {
        // RDKit✔️✔️:       pair.n2 = *pair.nbrbeg;
        // RDKit✔️✔️:       ++pair.nbrbeg;
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       pair.n2 = n2;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else if (t1_len > core_len && t2_len > core_len) {
        // RDKit✔️✔️:     while (pair.n2 < n2 &&
        // RDKit✔️✔️:            (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
        // RDKit✔️✔️:       pair.n2++;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) {
        // RDKit✔️✔️:       pair.n2++;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return pair.n1 < n1 && pair.n2 < n2;
        // RDKit✔️✔️: }
        // Complexity review: both versions scan at most O(V) unmapped nodes
        // outside the terminal branch and O(degree) adjacency entries in the
        // VF2+ branch, with no allocation per candidate pair.
        // RDKit✔️✔️: if (pair.n1 == NULL_NODE) pair.n1 = 0;
        // RDKit✔️✔️: if (pair.n2 == NULL_NODE) pair.n2 = 0;
        // RDKit✔️✔️: else pair.n2++;
        if pair.n1 == NULL_NODE {
            pair.n1 = 0;
        }
        if pair.n2 == NULL_NODE {
            pair.n2 = 0;
        } else {
            pair.n2 += 1;
        }

        // --- Select query node (n1) ---
        // RDKit✔️✔️: if (t1_len > core_len && t2_len > core_len) {
        if self.t1_len > self.core_len && self.t2_len > self.core_len {
            // RDKit✔️✔️: while (pair.n1 < n1 &&
            // RDKit✔️✔️:   (core_1[pair.n1] != NULL_NODE || term_1[pair.n1] == 0)) {
            // RDKit✔️✔️:   pair.n1++; pair.n2 = 0;
            // RDKit✔️✔️: }
            while pair.n1 < self.n1
                && (self.core_1[pair.n1] != NULL_NODE || self.term_1[pair.n1] == 0)
            {
                pair.n1 += 1;
                pair.n2 = 0;
            }
            // RDKit✔️✔️: /* Initialize VF2 Plus neighbor iterator.
            // RDKit✔️✔️:  * The next query node (pair.n1) has been selected from the terminal
            // RDKit✔️✔️:  * set and is therefore adjacent to an already mapped atom (in
            // RDKit✔️✔️:  * core_1). Rather than select pair.n2 from all atoms (0...n2) we can
            // RDKit✔️✔️:  * select it from the neighbors of this mapped atom (0...deg(nbor))
            // RDKit✔️✔️:  * since it must also be adajcent to this mapped atom!
            // RDKit✔️✔️:  */
            // RDKit✔️✔️: if (!pair.hasiter) {
            // RDKit✔️✔️:   boost::tie(n1iter_beg, n1iter_end) =
            // RDKit✔️✔️:       boost::adjacent_vertices(pair.n1, *g1);
            // RDKit✔️✔️:   while (n1iter_beg != n1iter_end && core_1[*n1iter_beg] == NULL_NODE) {
            // RDKit✔️✔️:     ++n1iter_beg;
            // RDKit✔️✔️:   }
            // RDKit✔️✔️:   assert(n1iter_beg != n1iter_end);
            // RDKit✔️✔️:   boost::tie(pair.nbrbeg, pair.nbrend) =
            // RDKit✔️✔️:       boost::adjacent_vertices(core_1[*n1iter_beg], *g2);
            // RDKit✔️✔️:   pair.hasiter = true;
            // RDKit✔️✔️: }
            if !pair.hasiter {
                let mut mapped_terminal_neighbor = NULL_NODE;
                for (query_neighbor, _) in self.g1.neighbor_row(pair.n1).iter() {
                    if self.core_1[query_neighbor] != NULL_NODE {
                        mapped_terminal_neighbor = self.core_1[query_neighbor];
                        break;
                    }
                }
                debug_assert_ne!(mapped_terminal_neighbor, NULL_NODE);
                if mapped_terminal_neighbor != NULL_NODE {
                    pair.nbr_node = mapped_terminal_neighbor;
                    pair.nbr_cursor = 0;
                    pair.nbr_end = self.g2.neighbor_row(mapped_terminal_neighbor).len();
                    pair.hasiter = true;
                }
            }
        } else if pair.n1 == 0 {
            // RDKit✔️✔️: } else if (pair.n1 == 0 && order != nullptr) {
            if let Some(order) = &self.order {
                // RDKit✔️✔️:   unsigned int i = 0;
                // RDKit✔️✔️:   while (i < n1 && core_1[pair.n1 = order[i]] != NULL_NODE) { i++; }
                // RDKit✔️✔️:   if (i == n1) pair.n1 = n1;
                let mut i = 0;
                while i < self.n1 {
                    let candidate = order[i];
                    if self.core_1[candidate] == NULL_NODE {
                        pair.n1 = candidate;
                        break;
                    }
                    i += 1;
                }
                if i == self.n1 {
                    pair.n1 = self.n1;
                }
            } else {
                // RDKit✔️✔️: } else {
                // RDKit✔️✔️:   while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
                // RDKit✔️✔️:     pair.n1++; pair.n2 = 0;
                // RDKit✔️✔️:   }
                while pair.n1 < self.n1 && self.core_1[pair.n1] != NULL_NODE {
                    pair.n1 += 1;
                    pair.n2 = 0;
                }
            }
        } else {
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   while (pair.n1 < n1 && core_1[pair.n1] != NULL_NODE) {
            // RDKit✔️✔️:     pair.n1++; pair.n2 = 0;
            // RDKit✔️✔️:   }
            while pair.n1 < self.n1 && self.core_1[pair.n1] != NULL_NODE {
                pair.n1 += 1;
                pair.n2 = 0;
            }
        }

        // --- Select mol node (n2) ---
        // RDKit✔️✔️: if (pair.hasiter) { ... }
        if pair.hasiter {
            let neighbors = self.g2.neighbor_row(pair.nbr_node);
            // RDKit✔️✔️: while (pair.nbrbeg < pair.nbrend && core_2[*pair.nbrbeg] != NULL_NODE) {
            // RDKit✔️✔️:   ++pair.nbrbeg;
            // RDKit✔️✔️: }
            while pair.nbr_cursor < pair.nbr_end
                && self.core_2[neighbors
                    .get(pair.nbr_cursor)
                    .expect("VF2+ cursor is within its borrowed neighbor row")
                    .0]
                    != NULL_NODE
            {
                pair.nbr_cursor += 1;
            }
            // RDKit✔️✔️: if (pair.nbrbeg < pair.nbrend) {
            // RDKit✔️✔️:   pair.n2 = *pair.nbrbeg;
            // RDKit✔️✔️:   ++pair.nbrbeg;
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   pair.n2 = n2;
            // RDKit✔️✔️: }
            if pair.nbr_cursor < pair.nbr_end {
                pair.n2 = neighbors
                    .get(pair.nbr_cursor)
                    .expect("VF2+ cursor is within its borrowed neighbor row")
                    .0;
                pair.nbr_cursor += 1;
            } else {
                pair.n2 = self.n2;
            }
        } else if self.t1_len > self.core_len && self.t2_len > self.core_len {
            // RDKit✔️✔️: } else if (t1_len > core_len && t2_len > core_len) {
            // RDKit✔️✔️:   while (pair.n2 < n2 &&
            // RDKit✔️✔️:     (core_2[pair.n2] != NULL_NODE || term_2[pair.n2] == 0)) {
            // RDKit✔️✔️:     pair.n2++;
            // RDKit✔️✔️:   }
            while pair.n2 < self.n2
                && (self.core_2[pair.n2] != NULL_NODE || self.term_2[pair.n2] == 0)
            {
                pair.n2 += 1;
            }
        } else {
            // RDKit✔️✔️: } else {
            // RDKit✔️✔️:   while (pair.n2 < n2 && core_2[pair.n2] != NULL_NODE) { pair.n2++; }
            // RDKit✔️✔️: }
            while pair.n2 < self.n2 && self.core_2[pair.n2] != NULL_NODE {
                pair.n2 += 1;
            }
        }

        // RDKit✔️✔️: return pair.n1 < n1 && pair.n2 < n2;
        pair.n1 < self.n1 && pair.n2 < self.n2
    }

    // RDKit source (vf2.hpp), IsFeasiblePair:
    //   bool IsFeasiblePair(node_id node1, node_id node2) {
    //     assert(node1 < n1); assert(node2 < n2);
    //     assert(core_1[node1] == NULL_NODE); assert(core_2[node2] == NULL_NODE);
    //
    //     // O(1) check for adjacency list
    //     if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
    //       return false;
    //     }
    //     if (!vc(node1, node2)) { return false; }
    //
    //     unsigned int other1, other2;
    //     // Check the out edges of node1
    //     typename Graph::out_edge_iterator bNbrs, eNbrs;
    //     boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
    //     while (bNbrs != eNbrs) {
    //       other1 = getOtherIdx(*g1, *bNbrs, node1);
    //       if (core_1[other1] != NULL_NODE) {
    //         other2 = core_1[other1];
    //         typename Graph::edge_descriptor oEdge;
    //         bool found;
    //         boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
    //         if (!found || !ec(*bNbrs, oEdge)) { return false; }
    //       }
    //       ++bNbrs;
    //     }
    //     return true;
    //   }

    /// RDKit✔️❌: IsFeasiblePair — check if (node1, node2) can be added.
    ///
    /// Performs degree check, vertex compatibility, and edge compatibility
    /// for already-matched neighbors. RDK_VF2_PRUNING (terminal count
    /// pre-check) is not enabled — the C++ code also has it behind an
    /// ifdef that is not defined at the top of vf2.hpp.
    fn is_feasible_pair(
        &self,
        node1: NodeId,
        node2: NodeId,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
    ) -> bool {
        // RDKit✔️✔️: bool IsFeasiblePair(node_id node1, node_id node2) {
        // RDKit✔️✔️:   assert(node1 < n1);
        // RDKit✔️✔️:   assert(node2 < n2);
        // RDKit✔️✔️:   assert(core_1[node1] == NULL_NODE);
        // RDKit✔️✔️:   assert(core_2[node2] == NULL_NODE);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // std::cerr<<"  ifp:"<<node1<<"-"<<node2<<"
        // RDKit✔️✔️:   // "<<vs_compared->size()<<std::endl;
        // RDKit✔️✔️:   // int &isCompat=vs_compared[node1*n2+node2];
        // RDKit✔️✔️:   // if(isCompat==0){
        // RDKit✔️✔️:   //   isCompat=vc(node1,node2)?1:-1;
        // RDKit✔️✔️:   // }
        // RDKit✔️✔️:   // if( isCompat<0 ){
        // RDKit✔️✔️:   //   //std::cerr<<"  short1"<<std::endl;
        // RDKit✔️✔️:   //   return false;
        // RDKit✔️✔️:   // }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // O(1) check for adjacency list
        // RDKit✔️✔️:   if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   if (!vc(node1, node2)) {
        // RDKit✔️✔️:     return false;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   unsigned int other1, other2;
        // RDKit✔️✔️: #ifdef RDK_VF2_PRUNING
        // RDKit✔️✔️:   unsigned int term1 = 0, term2 = 0;
        // RDKit✔️✔️:   unsigned int new1 = 0, new2 = 0;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // Check the out edges of node1
        // RDKit✔️✔️:   typename Graph::out_edge_iterator bNbrs, eNbrs;
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     other1 = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:     if (core_1[other1] != NULL_NODE) {
        // RDKit✔️✔️:       other2 = core_1[other1];
        // RDKit✔️✔️:       typename Graph::edge_descriptor oEdge;
        // RDKit✔️✔️:       bool found;
        // RDKit✔️✔️:       boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
        // RDKit✔️✔️:       if (!found || !ec(*bNbrs, oEdge)) {
        // RDKit✔️✔️:         // std::cerr<<"  short2"<<std::endl;
        // RDKit✔️✔️:         return false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #ifdef RDK_VF2_PRUNING
        // RDKit✔️✔️:     else {
        // RDKit✔️✔️:       if (term_1[other1]) ++term1;
        // RDKit✔️✔️:       if (!term_1[other1]) ++new1;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️: #ifdef RDK_VF2_PRUNING
        // RDKit✔️✔️:   // Check the out edges of node2
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     other2 = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:     if (core_2[other2] != NULL_NODE) {
        // RDKit✔️✔️:       // do nothing
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       if (term_2[other2]) ++term2;
        // RDKit✔️✔️:       if (!term_2[other2]) ++new2;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   // std::cerr<<(termin1 <= termin2 && termout1 <= termout2 &&
        // RDKit✔️✔️:   // (termin1+termout1+new1)<=(termin2+termout2+new2))<<std::endl;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // n.b. term1+new1 == boost::out_degree(node1) and
        // RDKit✔️✔️:   //      term2+new2 == boost::out_degree(node2)
        // RDKit✔️✔️:   return term1 <= term2 && (term1 + new1) <= (term2 + new2);
        // RDKit✔️✔️: #else
        // RDKit✔️✔️:   return true;
        // RDKit✔️✔️: #endif
        // RDKit✔️✔️: }
        // Complexity review: both active builds do O(1) degree and vertex
        // checks, scan O(degree(node1)) query edges, and perform target edge
        // lookup in O(degree(node2)); neither allocates per candidate.
        // RDKit✔️✔️: assert(node1 < n1); assert(node2 < n2);
        // RDKit✔️✔️: assert(core_1[node1] == NULL_NODE);
        // RDKit✔️✔️: assert(core_2[node2] == NULL_NODE);
        debug_assert!(node1 < self.n1);
        debug_assert!(node2 < self.n2);
        debug_assert_eq!(self.core_1[node1], NULL_NODE);
        debug_assert_eq!(self.core_2[node2], NULL_NODE);
        if self.core_1[node1] != NULL_NODE || self.core_2[node2] != NULL_NODE {
            return false;
        }

        // RDKit✔️✔️: if (boost::out_degree(node1, *g1) > boost::out_degree(node2, *g2)) {
        // RDKit✔️✔️:   return false;
        // RDKit✔️✔️: }
        if self.g1.out_degree(node1) > self.g2.out_degree(node2) {
            return false;
        }

        // RDKit✔️✔️: if (!vc(node1, node2)) { return false; }
        if !atom_fn(node1, node2) {
            return false;
        }

        // RDKit✔️✔️: // Check the out edges of node1
        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   other1 = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (core_1[other1] != NULL_NODE) {
        // RDKit✔️✔️:     other2 = core_1[other1];
        // RDKit✔️✔️:     if (!found || !ec(*bNbrs, oEdge)) { return false; }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for (_, edge_idx1) in self.g1.neighbor_row(node1).iter() {
            let other1 = get_other_idx(self.g1, edge_idx1, node1);
            if other1 == node1 {
                continue;
            }
            if self.core_1[other1] != NULL_NODE {
                let other2 = self.core_1[other1];
                // Check that (node2, other2) has a matching bond.
                let bond_found = self.find_bond(node2, other2);
                match bond_found {
                    Some(edge_idx2) => {
                        if !bond_fn(edge_idx1, edge_idx2) {
                            return false;
                        }
                    }
                    None => return false,
                }
            }
        }

        true
    }

    /// Find a bond between atom `a` and `b` in the molecule graph (g2).
    fn find_bond(&self, a: NodeId, b: NodeId) -> Option<usize> {
        // RDKit✔️✔️:         boost::tie(oEdge, found) = boost::edge(node2, other2, *g2);
        // The existing target edge descriptor is found by the same ordered
        // incident-row scan; the simple validated topology has one edge for
        // the endpoint pair and this returns its canonical bond index.
        for (nbr, bond_idx) in self.g2.neighbor_row(a).iter() {
            if nbr == b {
                return Some(bond_idx);
            }
        }
        None
    }

    fn add_pair(&mut self, node1: NodeId, node2: NodeId) {
        // RDKit✔️✔️: void AddPair(node_id node1, node_id node2) {
        // RDKit✔️✔️:   assert(node1 < n1);
        // RDKit✔️✔️:   assert(node2 < n2);
        // RDKit✔️✔️:   assert(core_len < n1);
        // RDKit✔️✔️:   assert(core_len < n2);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   ++core_len;
        // RDKit✔️✔️:   if (!term_1[node1]) {
        // RDKit✔️✔️:     term_1[node1] = core_len;
        // RDKit✔️✔️:     ++t1_len;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (!term_2[node2]) {
        // RDKit✔️✔️:     term_2[node2] = core_len;
        // RDKit✔️✔️:     ++t2_len;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   core_1[node1] = node2;
        // RDKit✔️✔️:   core_2[node2] = node1;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   typename Graph::out_edge_iterator bNbrs, eNbrs;
        // RDKit✔️✔️:   // FIX: this is explicitly ignoring directionality
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:     if (!term_1[other]) {
        // RDKit✔️✔️:       term_1[other] = core_len;
        // RDKit✔️✔️:       ++t1_len;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   // FIX: this is explicitly ignoring directionality
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:     if (!term_2[other]) {
        // RDKit✔️✔️:       term_2[other] = core_len;
        // RDKit✔️✔️:       ++t2_len;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Complexity review: both versions update O(1) mapping fields and scan
        // each selected node's adjacency once in O(degree1 + degree2), without
        // allocation or whole-graph rescanning.
        debug_assert!(node1 < self.n1);
        debug_assert!(node2 < self.n2);
        debug_assert!(self.core_len < self.n1);
        debug_assert!(self.core_len < self.n2);
        // RDKit✔️✔️: ++core_len;
        self.core_len += 1;
        let depth = self.core_len;

        // RDKit✔️✔️: if (!term_1[node1]) { term_1[node1] = core_len; ++t1_len; }
        if self.term_1[node1] == 0 {
            self.term_1[node1] = depth;
            self.t1_len += 1;
        }

        // RDKit✔️✔️: if (!term_2[node2]) { term_2[node2] = core_len; ++t2_len; }
        if self.term_2[node2] == 0 {
            self.term_2[node2] = depth;
            self.t2_len += 1;
        }

        // RDKit✔️✔️: core_1[node1] = node2; core_2[node2] = node1;
        self.core_1[node1] = node2;
        self.core_2[node2] = node1;

        // RDKit✔️✔️: // FIX: explicitly ignoring directionality
        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (!term_1[other]) { term_1[other] = core_len; ++t1_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for (_, edge) in self.g1.neighbor_row(node1).iter() {
            let other = get_other_idx(self.g1, edge, node1);
            if other == node1 {
                continue;
            }
            if self.term_1[other] == 0 {
                self.term_1[other] = depth;
                self.t1_len += 1;
            }
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:   if (!term_2[other]) { term_2[other] = core_len; ++t2_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for (_, edge) in self.g2.neighbor_row(node2).iter() {
            let other = get_other_idx(self.g2, edge, node2);
            if other == node2 {
                continue;
            }
            if self.term_2[other] == 0 {
                self.term_2[other] = depth;
                self.t2_len += 1;
            }
        }
    }

    fn back_track(&mut self, node1: NodeId, node2: NodeId) {
        // RDKit✔️✔️: void BackTrack(node_id node1, node_id node2) {
        // RDKit✔️✔️:   if (term_1[node1] == core_len) {
        // RDKit✔️✔️:     term_1[node1] = 0;
        // RDKit✔️✔️:     --t1_len;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   typename Graph::out_edge_iterator bNbrs, eNbrs;
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:     if (term_1[other] == core_len) {
        // RDKit✔️✔️:       term_1[other] = 0;
        // RDKit✔️✔️:       --t1_len;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   if (term_2[node2] == core_len) {
        // RDKit✔️✔️:     term_2[node2] = 0;
        // RDKit✔️✔️:     --t2_len;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️:   while (bNbrs != eNbrs) {
        // RDKit✔️✔️:     unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:     if (term_2[other] == core_len) {
        // RDKit✔️✔️:       term_2[other] = 0;
        // RDKit✔️✔️:       --t2_len;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     ++bNbrs;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   core_1[node1] = NULL_NODE;
        // RDKit✔️✔️:   core_2[node2] = NULL_NODE;
        // RDKit✔️✔️:   --core_len;
        // RDKit✔️✔️: }
        // Complexity review: both versions scan each removed node's adjacency
        // once in O(degree1 + degree2), mutate depth-tagged entries in place,
        // and allocate no temporary collections.
        let depth = self.core_len;

        // RDKit✔️✔️: if (term_1[node1] == core_len) { term_1[node1] = 0; --t1_len; }
        if self.term_1[node1] == depth {
            self.term_1[node1] = 0;
            self.t1_len -= 1;
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node1, *g1);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g1, *bNbrs, node1);
        // RDKit✔️✔️:   if (term_1[other] == core_len) { term_1[other] = 0; --t1_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for (_, edge) in self.g1.neighbor_row(node1).iter() {
            let other = get_other_idx(self.g1, edge, node1);
            if other == node1 {
                continue;
            }
            if self.term_1[other] == depth {
                self.term_1[other] = 0;
                self.t1_len -= 1;
            }
        }

        // RDKit✔️✔️: if (term_2[node2] == core_len) { term_2[node2] = 0; --t2_len; }
        if self.term_2[node2] == depth {
            self.term_2[node2] = 0;
            self.t2_len -= 1;
        }

        // RDKit✔️✔️: boost::tie(bNbrs, eNbrs) = boost::out_edges(node2, *g2);
        // RDKit✔️✔️: while (bNbrs != eNbrs) {
        // RDKit✔️✔️:   unsigned int other = getOtherIdx(*g2, *bNbrs, node2);
        // RDKit✔️✔️:   if (term_2[other] == core_len) { term_2[other] = 0; --t2_len; }
        // RDKit✔️✔️:   ++bNbrs;
        // RDKit✔️✔️: }
        for (_, edge) in self.g2.neighbor_row(node2).iter() {
            let other = get_other_idx(self.g2, edge, node2);
            if other == node2 {
                continue;
            }
            if self.term_2[other] == depth {
                self.term_2[other] = 0;
                self.t2_len -= 1;
            }
        }

        // RDKit✔️✔️: core_1[node1] = NULL_NODE;
        // RDKit✔️✔️: core_2[node2] = NULL_NODE;
        // RDKit✔️✔️: --core_len;
        self.core_1[node1] = NULL_NODE;
        self.core_2[node2] = NULL_NODE;
        self.core_len -= 1;
    }

    fn get_core_set_into(&self, c1: &mut [NodeId], c2: &mut [NodeId]) -> usize {
        // RDKit❗✔️: void GetCoreSet(node_id c1[], node_id c2[]) {
        // RDKit❗✔️:   unsigned int i, j;
        // RDKit❗✔️:   for (i = 0, j = 0; i < n1; ++i) {
        // RDKit❗✔️:     if (core_1[i] != NULL_NODE) {
        // RDKit❗✔️:       c1[j] = i;
        // RDKit❗✔️:       c2[j] = core_1[i];
        // RDKit❗✔️:       ++j;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️: }
        // Behavior: query indices are scanned ascending; only the populated
        // prefix is written, and the returned count identifies that prefix.
        // Complexity: one O(n1) scan, core_len paired indexed writes, and no
        // allocation, matching RDKit's caller-array contract.
        debug_assert!(c1.len() >= self.core_len);
        debug_assert!(c2.len() >= self.core_len);
        let mut written = 0;
        for (query_index, &target_index) in self.core_1.iter().enumerate() {
            if target_index != NULL_NODE {
                c1[written] = query_index;
                c2[written] = target_index;
                written += 1;
            }
        }
        written
    }

    fn match_one(
        &mut self,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
        mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
        c1: &mut [NodeId],
        c2: &mut [NodeId],
    ) -> bool {
        // RDKit❗✔️: bool Match(node_id c1[], node_id c2[]) {
        // RDKit❗✔️:   if (IsGoal()) {
        // RDKit❗✔️:     GetCoreSet(c1, c2);
        // RDKit❗✔️:     if (MatchChecks(c1, c2)) {
        // RDKit❗✔️:       return true;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (IsDead()) {
        // RDKit❗✔️:     return false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   Pair<Graph> pair;
        // RDKit❗✔️:   while (NextPair(pair)) {
        // RDKit❗✔️:     if (IsFeasiblePair(pair.n1, pair.n2)) {
        // RDKit❗✔️:       AddPair(pair.n1, pair.n2);
        // RDKit❗✔️:       if (Match(c1, c2)) {  // recurse
        // RDKit❗✔️:         return true;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       BackTrack(pair.n1, pair.n2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // Behavior: each goal writes the caller-owned prefix, its final check
        // sees only that prefix, and success returns on the same DFS branch.
        // Complexity: query-sized buffers are allocated once by the outer
        // invocation and reused here, with no mapping allocation at each goal.
        if self.is_goal() {
            let written = self.get_core_set_into(c1, c2);
            debug_assert_eq!(written, self.core_len);
            let accepted = match match_check.as_mut() {
                Some(check) => self.match_checks(&c1[..written], &c2[..written], check),
                None => true,
            };
            if accepted {
                return true;
            }
        }
        if self.is_dead() {
            return false;
        }
        let mut pair = Vf2Pair::new();
        while self.next_pair(&mut pair) {
            if self.is_feasible_pair(pair.n1, pair.n2, atom_fn, bond_fn) {
                self.add_pair(pair.n1, pair.n2);
                if self.match_one(atom_fn, bond_fn, match_check.as_deref_mut(), c1, c2) {
                    return true;
                }
                self.back_track(pair.n1, pair.n2);
            }
        }
        false
    }

    fn match_all(
        &mut self,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
        mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
        c1: &mut [NodeId],
        c2: &mut [NodeId],
        results: &mut Vec<Vec<(NodeId, NodeId)>>,
        max_matches: usize,
    ) -> bool {
        // RDKit❗✔️: template <class DoubleBackInsertionSequence>
        // RDKit❗✔️: bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
        // RDKit❗✔️:               unsigned int lim = 0) {
        // RDKit❗✔️:   if (IsGoal()) {
        // RDKit❗✔️:     GetCoreSet(c1, c2);
        // RDKit❗✔️:     if (MatchChecks(c1, c2)) {
        // RDKit❗✔️:       typename DoubleBackInsertionSequence::value_type newSeq;
        // RDKit❗✔️:       newSeq.reserve(core_len);
        // RDKit❗✔️:       for (unsigned int i = 0; i < core_len; ++i) {
        // RDKit❗✔️:         newSeq.emplace_back(c1[i], c2[i]);
        // RDKit❗✔️:       }
        // RDKit❗✔️:       res.push_back(newSeq);
        // RDKit❗✔️:       return lim && res.size() >= lim;
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (IsDead()) {
        // RDKit❗✔️:     return false;
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   Pair<Graph> pair;
        // RDKit❗✔️:   while (NextPair(pair)) {
        // RDKit❗✔️:     if (IsFeasiblePair(pair.n1, pair.n2)) {
        // RDKit❗✔️:       AddPair(pair.n1, pair.n2);
        // RDKit❗✔️:       if (MatchAll(c1, c2, res, lim)) {  // recurse
        // RDKit❗✔️:         return true;
        // RDKit❗✔️:       }
        // RDKit❗✔️:       BackTrack(pair.n1, pair.n2);
        // RDKit❗✔️:     }
        // RDKit❗✔️:   }
        // RDKit❗✔️:   return false;
        // RDKit❗✔️: }
        // Behavior: goals are checked before dead-state handling; each accepted
        // mapping is appended before the limit check and source-order recursion
        // backtracks only after an unaccepted/continuing child returns false.
        // Complexity: the caller-owned arrays are reused at every goal; one
        // capacity-sized pair vector is created only for each accepted mapping.
        if self.is_goal() {
            let written = self.get_core_set_into(c1, c2);
            debug_assert_eq!(written, self.core_len);
            let accepted = match match_check.as_mut() {
                Some(check) => self.match_checks(&c1[..written], &c2[..written], check),
                None => true,
            };
            if accepted {
                let mut new_sequence = Vec::with_capacity(written);
                new_sequence.extend(
                    c1[..written]
                        .iter()
                        .copied()
                        .zip(c2[..written].iter().copied()),
                );
                results.push(new_sequence);
                return max_matches > 0 && results.len() >= max_matches;
            }
        }
        if self.is_dead() {
            return false;
        }
        let mut pair = Vf2Pair::new();
        while self.next_pair(&mut pair) {
            if self.is_feasible_pair(pair.n1, pair.n2, atom_fn, bond_fn) {
                self.add_pair(pair.n1, pair.n2);
                if self.match_all(
                    atom_fn,
                    bond_fn,
                    match_check.as_deref_mut(),
                    c1,
                    c2,
                    results,
                    max_matches,
                ) {
                    return true;
                }
                self.back_track(pair.n1, pair.n2);
            }
        }
        false
    }
}

// ---------------------------------------------------------------------------
// VF2 recursive matching
// ---------------------------------------------------------------------------
//
// RDKit source (vf2.hpp):
//   bool Match(node_id c1[], node_id c2[]) {
//     if (IsGoal()) { GetCoreSet(c1, c2); if (MatchChecks(c1, c2)) return true; }
//     if (IsDead()) return false;
//     Pair<Graph> pair;
//     while (NextPair(pair)) {
//       if (IsFeasiblePair(pair.n1, pair.n2)) {
//         AddPair(pair.n1, pair.n2);
//         if (Match(c1, c2)) return true;  // recurse
//         BackTrack(pair.n1, pair.n2);
//       }
//     }
//     return false;
//   }

/// RDKit❗✔️: Match — find first match via VF2 recursion.
///
/// Matches RDKit's `Match(c1, c2)` entry point. `match_check` allows
/// final verification (like MolMatchFinalCheckFunctor). If None, all
/// completed matches are accepted.
fn vf2_match(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    c1: &mut [NodeId],
    c2: &mut [NodeId],
) -> bool {
    // RDKit❗✔️: template <class SubState>
    // RDKit❗✔️: bool match(int *pn, node_id c1[], node_id c2[], SubState &s) {
    // RDKit❗✔️:   if (s.Match(c1, c2)) {
    // RDKit❗✔️:     // not needed, pn = num query atoms (n1)...
    // RDKit❗✔️:     *pn = s.CoreLen();
    // RDKit❗✔️:     return true;
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return false;
    // RDKit❗✔️: }
    // Rust forwards the invocation-owned arrays through the same recursive
    // member; vf2_entry_one projects the accepted prefix after success.
    state.match_one(atom_fn, bond_fn, match_check, c1, c2)
}

// RDKit source (vf2.hpp), MatchAll:
//   template <class DoubleBackInsertionSequence>
//   bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
//                 unsigned int lim = 0) {
//     if (IsGoal()) {
//       GetCoreSet(c1, c2);
//       if (MatchChecks(c1, c2)) {
//         typename DoubleBackInsertionSequence::value_type newSeq;
//         newSeq.reserve(core_len);
//         for (unsigned int i = 0; i < core_len; ++i) {
//           newSeq.emplace_back(c1[i], c2[i]);
//         }
//         res.push_back(newSeq);
//         return lim && res.size() >= lim;
//       }
//     }
//     if (IsDead()) return false;
//     Pair<Graph> pair;
//     while (NextPair(pair)) {
//       if (IsFeasiblePair(pair.n1, pair.n2)) {
//         AddPair(pair.n1, pair.n2);
//         if (MatchAll(c1, c2, res, lim)) return true;  // recurse
//         BackTrack(pair.n1, pair.n2);
//       }
//     }
//     return false;
//   }

/// RDKit❗✔️: MatchAll — find all matches up to `max_matches`.
///
/// Collects each accepted mapping into `results` as one ordered paired sequence.
/// Returns true when the limit has been reached, signaling the caller
/// to stop.
fn vf2_match_all(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    c1: &mut [NodeId],
    c2: &mut [NodeId],
    results: &mut Vec<Vec<(NodeId, NodeId)>>,
    max_matches: usize,
) -> bool {
    // RDKit❗✔️: template <class SubState, class DoubleBackInsertionSequence>
    // RDKit❗✔️: bool match(node_id c1[], node_id c2[], SubState &s,
    // RDKit❗✔️:            DoubleBackInsertionSequence &res, unsigned int max_results) {
    // RDKit❗✔️:   s.MatchAll(c1, c2, res, max_results);
    // RDKit❗✔️:   return !res.empty();
    // RDKit❗✔️: }
    // Complexity review: this wrapper adds one emptiness check after invoking
    // the same member recursion and forwards caller-owned scratch without
    // copying or re-enumerating results.
    state.match_all(atom_fn, bond_fn, match_check, c1, c2, results, max_matches);
    !results.is_empty()
}

fn vf2_entry_one(
    g1: Vf2GraphRef<'_>,
    g2: Vf2GraphRef<'_>,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    result: &mut Vec<(NodeId, NodeId)>,
) -> bool {
    // RDKit✔️✔️: template <
    // RDKit✔️✔️:     class Graph, class VertexLabeling  // binary predicate
    // RDKit✔️✔️:     ,
    // RDKit✔️✔️:     class EdgeLabeling  // binary predicate
    // RDKit✔️✔️:     ,
    // RDKit✔️✔️:     class MatchChecking  // binary predicate
    // RDKit✔️✔️:     ,
    // RDKit✔️✔️:     class
    // RDKit✔️✔️:     BackInsertionSequence  // contains
    // RDKit✔️✔️:                            // std::pair<vertex_descriptor,vertex_descriptor>
    // RDKit✔️✔️:     >
    // RDKit✔️✔️: bool vf2(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
    // RDKit✔️✔️:          EdgeLabeling &edge_labeling, MatchChecking &match_checking,
    // RDKit✔️✔️:          BackInsertionSequence &F) {
    // RDKit✔️✔️:   detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling, MatchChecking>
    // RDKit✔️✔️:       s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false);
    // RDKit✔️✔️:   detail::node_id *ni1 = new detail::node_id[num_vertices(g1)];
    // RDKit✔️✔️:   detail::node_id *ni2 = new detail::node_id[num_vertices(g2)];
    // RDKit✔️✔️:   int n = 0;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   F.clear();
    // RDKit✔️✔️:   if (match(&n, ni1, ni2, s0)) {
    // RDKit✔️✔️:     auto sz = num_vertices(g1);
    // RDKit✔️✔️:     F.reserve(sz);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < sz; ++i) {
    // RDKit✔️✔️:       F.emplace_back(ni1[i], ni2[i]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   delete[] ni1;
    // RDKit✔️✔️:   delete[] ni2;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return !F.empty();
    // RDKit✔️✔️: };
    // Behavior: the invocation keeps RDKit's unsorted state and first-match
    // DFS, then projects exactly the accepted query-order prefix once.
    // Complexity: these two O(V) arrays are the sole mapping scratch buffers
    // for the recursion; the result is materialized only after acceptance.
    let mut state = Vf2SubState::new(g1, g2, false);
    let mut c1 = vec![NULL_NODE; g1.num_atoms()];
    let mut c2 = vec![NULL_NODE; g1.num_atoms()];
    result.clear();
    if vf2_match(&mut state, atom_fn, bond_fn, match_check, &mut c1, &mut c2) {
        let matched = state.core_len;
        result.reserve(g1.num_atoms());
        result.extend(
            c1[..matched]
                .iter()
                .copied()
                .zip(c2[..matched].iter().copied()),
        );
    }
    !result.is_empty()
}

fn vf2_entry_all(
    g1: Vf2GraphRef<'_>,
    g2: Vf2GraphRef<'_>,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    results: &mut Vec<Vec<(NodeId, NodeId)>>,
    max_results: usize,
) -> bool {
    vf2_entry_all_ordered(
        g1,
        g2,
        atom_fn,
        bond_fn,
        match_check,
        results,
        max_results,
        None,
    )
}

fn vf2_entry_all_ordered(
    g1: Vf2GraphRef<'_>,
    g2: Vf2GraphRef<'_>,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    results: &mut Vec<Vec<(NodeId, NodeId)>>,
    max_results: usize,
    order: Option<&[usize]>,
) -> bool {
    // RDKit❗✔️: template <class Graph, class VertexLabeling  // binary predicate
    // RDKit❗✔️:           ,
    // RDKit❗✔️:           class EdgeLabeling  // binary predicate
    // RDKit❗✔️:           ,
    // RDKit❗✔️:           class MatchChecking  // binary predicate
    // RDKit❗✔️:           ,
    // RDKit❗✔️:           class DoubleBackInsertionSequence  // contains a back insertion
    // RDKit❗✔️:                                              // sequence
    // RDKit❗✔️:           >
    // RDKit❗✔️: bool vf2_all(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
    // RDKit❗✔️:              EdgeLabeling &edge_labeling, MatchChecking &match_checking,
    // RDKit❗✔️:              DoubleBackInsertionSequence &F, unsigned int max_results = 1000) {
    // RDKit❗✔️:   detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling, MatchChecking>
    // RDKit❗✔️:       s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false);
    // RDKit❗✔️:   std::unique_ptr<detail::node_id[]> ni1(new detail::node_id[num_vertices(g1)]);
    // RDKit❗✔️:   std::unique_ptr<detail::node_id[]> ni2(new detail::node_id[num_vertices(g2)]);
    // RDKit❗✔️:
    // RDKit❗✔️:   F.clear();
    // RDKit❗✔️:   F.resize(0);
    // RDKit❗✔️:
    // RDKit❗✔️:   match(ni1.get(), ni2.get(), s0, F, max_results);
    // RDKit❗✔️:
    // RDKit❗✔️:   return !F.empty();
    // RDKit❗✔️: };
    // Behavior: the same unsorted state and recursive order are retained; the
    // limit stops only after an accepted mapping has entered the result list.
    // Complexity: the two query-sized scratch arrays are reused at all goals;
    // each accepted mapping gets one reserved vector of paired indices.
    let mut state = match order {
        Some(order) => Vf2SubState::with_order(g1, g2, order),
        None => Vf2SubState::new(g1, g2, false),
    };
    let mut c1 = vec![NULL_NODE; g1.num_atoms()];
    let mut c2 = vec![NULL_NODE; g1.num_atoms()];
    results.clear();
    vf2_match_all(
        &mut state,
        atom_fn,
        bond_fn,
        match_check,
        &mut c1,
        &mut c2,
        results,
        max_results,
    )
}

// ---------------------------------------------------------------------------
// Final match check (simplified MolMatchFinalCheckFunctor)
// ---------------------------------------------------------------------------
//
// RDKit source (SubstructMatch.cpp):
//   bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
//                                              const std::uint32_t m_c[]) {
//     if (d_params.extraFinalCheck || d_params.useGenericMatchers) { ... }
//     HashedStorageType match;
//     if (d_params.uniquify) {
//       match.resize(d_mol.getNumAtoms());
//       std::fill(match.begin(), match.end(), 0);
//       for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
//         match[m_c[i]] = 1;
//       }
//       if (matchesSeen.find(match) != matchesSeen.end()) { return false; }
//     }
//     if (!d_params.useChirality) {
//       if (d_params.uniquify) { matchesSeen.insert(match); }
//       return true;
//     }
//     // ... chirality checks ...
//   }

/// RDKit✔️✔️: Final match atom-set mask used for uniquification.
fn match_mask(atom_mapping: &[usize], mol_num_atoms: usize) -> Vec<bool> {
    let mut mask = vec![false; mol_num_atoms];
    for &ma in atom_mapping {
        if ma < mol_num_atoms {
            mask[ma] = true;
        }
    }
    mask
}

fn count_swaps_to_interconvert_i32(reference: &[i32], probe: &[i32]) -> Option<u32> {
    if reference.len() != probe.len() {
        return None;
    }
    let mut probe = probe.to_vec();
    let mut swaps = 0_u32;
    for (index, expected) in reference.iter().copied().enumerate() {
        if probe[index] == expected {
            continue;
        }
        let found = probe[index..]
            .iter()
            .position(|value| *value == expected)
            .map(|offset| index + offset)?;
        probe.swap(index, found);
        swaps = swaps.checked_add(1)?;
    }
    Some(swaps)
}

fn rdkit_atom_perturbation_order_from_bond_indices(
    mol: &QueryGraph,
    atom_idx: usize,
    probe: &[i32],
) -> Result<u32, SubstructMatchError> {
    // BEGIN RDKIT CPP FUNCTION Atom::getPerturbationOrder
    // RDKit✔️✔️: int Atom::getPerturbationOrder(const INT_LIST &probe) const {
    // RDKit✔️✔️:   INT_LIST ref;
    // RDKit✔️✔️:   for (const auto bond : getOwningMol().atomBonds(this)) {
    // RDKit✔️✔️:     ref.push_back(bond->getIdx());
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return static_cast<int>(countSwapsToInterconvert(probe, ref));
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    let reference: Vec<i32> = mol
        .adjacency()
        .get(atom_idx)
        .into_iter()
        .flatten()
        .map(|(_, bond)| i32::try_from(*bond))
        .collect::<Result<_, _>>()
        .map_err(|_| SubstructMatchError::Unsupported {
            branch: "MolMatchFinalCheckFunctor/Atom::getPerturbationOrder/bond-index-overflow",
            rdkit_function: "Atom::getPerturbationOrder",
        })?;
    count_swaps_to_interconvert_i32(probe, &reference).ok_or(SubstructMatchError::Unsupported {
        branch: "MolMatchFinalCheckFunctor/Atom::getPerturbationOrder/unmodeled-bond-ordering",
        rdkit_function: "Atom::getPerturbationOrder/countSwapsToInterconvert",
    })
}

fn rdkit_translate_ez_label_to_cis_trans(stereo: BondStereo) -> BondStereo {
    match stereo {
        BondStereo::E => BondStereo::Trans,
        BondStereo::Z => BondStereo::Cis,
        other => other,
    }
}

fn enhanced_stereo_is_ok(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    q_to_mol: &[NodeId],
    mol_stereo_groups: &[Option<usize>],
    matches: &[Option<bool>],
) -> bool {
    // RDKit✔️✔️: bool enhancedStereoIsOK(
    // RDKit✔️✔️:     const ROMol &mol, const ROMol &query,
    // RDKit✔️✔️:     std::unordered_map<unsigned int, unsigned int> &q_to_mol,
    // RDKit✔️✔️:     const std::unordered_map<unsigned int, StereoGroup const *>
    // RDKit✔️✔️:         &molStereoGroups,
    // RDKit✔️✔️:     const std::unordered_map<unsigned int, bool> &matches) {
    // RDKit✔️✔️:   std::unordered_map<unsigned int, StereoGroup const *> molAtomsToQueryGroups;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // If the query has stereo groups:
    // RDKit✔️✔️:   // * OR only matches AND or OR (not absolute)
    // RDKit✔️✔️:   // * AND only matches OR
    // RDKit✔️✔️:   for (const auto &sg : query.getStereoGroups()) {
    // RDKit✔️✔️:     if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // StereoGroup const* matched_mol_group = nullptr;
    // RDKit✔️✔️:     const bool is_and = sg.getGroupType() == StereoGroupType::STEREO_AND;
    // RDKit✔️✔️:     for (const auto a : sg.getAtoms()) {
    // RDKit✔️✔️:       const auto mol_group = molStereoGroups.find(q_to_mol[a->getIdx()]);
    // RDKit✔️✔️:       if (mol_group == molStereoGroups.end()) {
    // RDKit✔️✔️:         // group matching absolute. not ok.
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       } else if (is_and && mol_group->second->getGroupType() !=
    // RDKit✔️✔️:                                StereoGroupType::STEREO_AND) {
    // RDKit✔️✔️:         // AND matching OR. not ok.
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       molAtomsToQueryGroups[q_to_mol[a->getIdx()]] = &sg;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   // If the mol has stereo groups:
    // RDKit✔️✔️:   // * All atoms must either be the same or opposite, you can't mix
    // RDKit✔️✔️:   // * Only one stereogroup must cover all matched atoms in the mol stereo group
    // RDKit✔️✔️:   for (const auto &sg : mol.getStereoGroups()) {
    // RDKit✔️✔️:     if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     bool doesMatch = false;
    // RDKit✔️✔️:     bool seen = false;
    // RDKit✔️✔️:     StereoGroup const *QGroup = nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     for (const auto &a : sg.getAtoms()) {
    // RDKit✔️✔️:       auto thisDoesMatch = matches.find(a->getIdx());
    // RDKit✔️✔️:       if (thisDoesMatch == matches.end()) {
    // RDKit✔️✔️:         // not matched
    // RDKit✔️✔️:         continue;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       auto pos = molAtomsToQueryGroups.find(a->getIdx());
    // RDKit✔️✔️:       auto thisQGroup =
    // RDKit✔️✔️:           pos == molAtomsToQueryGroups.end() ? nullptr : pos->second;
    // RDKit✔️✔️:       if (!seen) {
    // RDKit✔️✔️:         doesMatch = thisDoesMatch->second;
    // RDKit✔️✔️:         QGroup = thisQGroup;
    // RDKit✔️✔️:         seen = true;
    // RDKit✔️✔️:       } else if (doesMatch != thisDoesMatch->second) {
    // RDKit✔️✔️:         // diastereomer. not ok.
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       } else if (thisQGroup != QGroup) {
    // RDKit✔️✔️:         // mix of groups in query. not ok.
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Complexity review: both implementations allocate O(mol atoms) lookup
    // state and scan each query/target group member once. Vec indexing replaces
    // unordered-map lookup with O(1) direct indexing and no worse allocation.
    let mut mol_atoms_to_query_groups = vec![None; mol.num_atoms()];
    for (query_group_idx, group) in query.stereo_groups().iter().enumerate() {
        if group.kind() == StereoGroupKind::Absolute {
            continue;
        }
        let is_and = group.kind() == StereoGroupKind::And;
        for atom in group.atoms() {
            let mol_atom = q_to_mol[atom.index()];
            let Some(mol_group_idx) = mol_stereo_groups[mol_atom] else {
                return false;
            };
            if is_and && mol.stereo_groups()[mol_group_idx].kind() != StereoGroupKind::And {
                return false;
            }
            mol_atoms_to_query_groups[mol_atom] = Some(query_group_idx);
        }
    }

    for group in mol.stereo_groups() {
        if group.kind() == StereoGroupKind::Absolute {
            continue;
        }
        let mut first: Option<(bool, Option<usize>)> = None;
        for atom in group.atoms() {
            let mol_atom = atom.index();
            let Some(does_match) = matches[mol_atom] else {
                continue;
            };
            let query_group = mol_atoms_to_query_groups[mol_atom];
            match first {
                None => first = Some((does_match, query_group)),
                Some((first_match, _)) if first_match != does_match => return false,
                Some((_, first_group)) if first_group != query_group => return false,
                Some(_) => {}
            }
        }
    }
    true
}

struct MolMatchFinalCheckSetup {
    mol_stereo_groups: Vec<Option<usize>>,
}

impl MolMatchFinalCheckSetup {
    fn new(_query: &QueryGraph, mol: &SearchTarget<'_>, params: &SubstructMatchParams) -> Self {
        // RDKit✔️✔️: MolMatchFinalCheckFunctor::MolMatchFinalCheckFunctor(
        // RDKit✔️✔️:     const ROMol &query, const ROMol &mol, const SubstructMatchParameters &ps)
        // RDKit✔️✔️:     : d_query(query), d_mol(mol), d_params(ps) {
        // RDKit✔️✔️:   if (d_params.useEnhancedStereo) {
        // RDKit✔️✔️:     for (const auto &sg : d_mol.getStereoGroups()) {
        // RDKit✔️✔️:       if (sg.getGroupType() == StereoGroupType::STEREO_ABSOLUTE) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       for (const auto a : sg.getAtoms()) {
        // RDKit✔️✔️:         d_molStereoGroups[a->getIdx()] = &sg;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: }
        // Complexity review: both build the group lookup once in O(mol atoms
        // plus non-absolute group members), then reuse it across goal checks.
        let mut mol_stereo_groups = vec![None; mol.num_atoms()];
        if params.use_enhanced_stereo {
            for (group_idx, group) in mol.stereo_groups().iter().enumerate() {
                if group.kind() == StereoGroupKind::Absolute {
                    continue;
                }
                for atom in group.atoms() {
                    mol_stereo_groups[atom.index()] = Some(group_idx);
                }
            }
        }
        Self { mol_stereo_groups }
    }
}

fn find_bond_between<'a>(mol: &'a SearchTarget<'_>, begin: usize, end: usize) -> Option<&'a Bond> {
    mol.bonds().iter().find(|bond| {
        let b = bond.begin().index();
        let e = bond.end().index();
        (b == begin && e == end) || (b == end && e == begin)
    })
}

fn find_query_bond_between(query: &QueryGraph, begin: usize, end: usize) -> Option<&QueryBond> {
    query.bonds().iter().find(|bond| {
        let (b, e) = bond.endpoints();
        (b == begin && e == end) || (b == end && e == begin)
    })
}

fn rdkit_match_final_check(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    c1: &[NodeId],
    c2: &[NodeId],
    setup: &MolMatchFinalCheckSetup,
    matches_seen: &mut HashSet<Vec<bool>>,
) -> Result<bool, SubstructMatchError> {
    // BEGIN RDKIT CPP FUNCTION MolMatchFinalCheckFunctor::operator()
    // RDKit✔️✔️: bool MolMatchFinalCheckFunctor::operator()(const std::uint32_t q_c[],
    // RDKit✔️✔️:                                            const std::uint32_t m_c[]) {
    // RDKit✔️✔️:   if (d_params.extraFinalCheck || d_params.useGenericMatchers) {
    // RDKit✔️✔️:     const std::span<const std::uint32_t> aids(m_c, d_query.getNumAtoms());
    // RDKit✔️✔️:     if (d_params.useGenericMatchers &&
    // RDKit✔️✔️:         !GenericGroups::genericAtomMatcher(d_mol, d_query, aids)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (d_params.extraFinalCheck && !d_params.extraFinalCheck(d_mol, aids)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // Complexity review: this adds the source-equivalent O(Q) dispatcher and
    // its selected generic-group matcher only when the option is enabled.
    if params.use_generic_matchers {
        if !crate::generic_groups::generic_atom_matcher(mol, query, c2) {
            return Ok(false);
        }
    }
    if let Some(extra_final_check) = &params.extra_final_check
        && !extra_final_check(mol, c2)
    {
        return Ok(false);
    }
    // RDKit✔️✔️:   HashedStorageType match;
    // RDKit✔️✔️:   if (d_params.uniquify) {
    // RDKit✔️✔️:     match.resize(d_mol.getNumAtoms());
    // RDKit✔️✔️:     std::fill(match.begin(), match.end(), 0);
    // RDKit✔️✔️:     for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    // RDKit✔️✔️:       match[m_c[i]] = 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (matchesSeen.find(match) != matchesSeen.end()) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut q_to_mol = vec![NULL_NODE; query.num_atoms()];
    for (&qa, &ma) in c1.iter().zip(c2.iter()) {
        if qa < q_to_mol.len() {
            q_to_mol[qa] = ma;
        }
    }
    let match_key = if params.uniquify {
        let mask = match_mask(&q_to_mol, mol.num_atoms());
        // RDKit's unordered_set gives expected constant-time membership after
        // hashing the target-atom mask. HashSet mirrors that lookup; only
        // membership is observed, while raw_matches retains VF2 result order.
        if matches_seen.contains(&mask) {
            return Ok(false);
        }
        Some(mask)
    } else {
        None
    };

    // RDKit✔️✔️:   if (!d_params.useChirality) {
    // RDKit✔️✔️:     if (d_params.uniquify) {
    // RDKit✔️✔️:       matchesSeen.insert(match);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    if !params.use_chirality {
        if let Some(mask) = match_key {
            matches_seen.insert(mask);
        }
        return Ok(true);
    }

    // RDKit✔️✔️:   std::unordered_map<unsigned int, bool> matches;
    let mol_stereo_groups = &setup.mol_stereo_groups;
    let mut stereo_matches = vec![None; mol.num_atoms()];

    // RDKit✔️✔️:   // check chiral atoms:
    // RDKit✔️✔️:   for (unsigned int i = 0; i < d_query.getNumAtoms(); ++i) {
    // RDKit✔️✔️:     const Atom *qAt = d_query.getAtomWithIdx(q_c[i]);
    // RDKit✔️✔️:     if (qAt->getDegree() < 3 || !detail::hasChiralLabel(qAt)) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for qi in 0..query.num_atoms() {
        let q_at = &query.atoms()[qi];
        if query.adjacency().get(qi).map_or(0, Vec::len) < 3 || !has_chiral_label(q_at.chiral_tag())
        {
            continue;
        }
        let mi = q_to_mol[qi];
        let m_at = &mol.atoms()[mi];
        // RDKit✔️✔️:     if (!detail::hasChiralLabel(mAt)) {
        // RDKit✔️✔️:       if (d_params.specifiedStereoQueryMatchesUnspecified) {
        // RDKit✔️✔️:         continue;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if !has_chiral_label(m_at.chiral_tag()) {
            if params.specified_stereo_query_matches_unspecified {
                continue;
            }
            return Ok(false);
        }
        // RDKit✔️✔️:     if (qAt->getDegree() > mAt->getDegree()) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if query.adjacency().get(qi).map_or(0, Vec::len)
            > mol.topology_block().adjacency.neighbors_of(mi).len()
        {
            return Ok(false);
        }

        // RDKit✔️✔️:     INT_LIST qOrder;
        // RDKit✔️✔️:     INT_LIST mOrder;
        // RDKit✔️✔️:     for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
        // RDKit✔️✔️:       const Bond *qB = d_query.getBondBetweenAtoms(q_c[i], q_c[j]);
        // RDKit✔️✔️:       const Bond *mB = d_mol.getBondBetweenAtoms(m_c[i], m_c[j]);
        // RDKit✔️✔️:       if (qB && mB) {
        // RDKit✔️✔️:         mOrder.push_back(mB->getIdx());
        // RDKit✔️✔️:         qOrder.push_back(qB->getIdx());
        // RDKit✔️✔️:         if (mOrder.size() == qAt->getDegree()) {
        // RDKit✔️✔️:           break;
        // RDKit✔️✔️:         }
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut q_order: Vec<i32> = Vec::new();
        let mut m_order: Vec<i32> = Vec::new();
        for qj in 0..query.num_atoms() {
            let Some(q_bond) = find_query_bond_between(query, qi, qj) else {
                continue;
            };
            let mj = q_to_mol[qj];
            let Some(m_bond) = find_bond_between(mol, mi, mj) else {
                continue;
            };
            q_order.push(i32::try_from(q_bond.id().index()).map_err(|_| {
                SubstructMatchError::Unsupported {
                    branch: "MolMatchFinalCheckFunctor/qOrder/bond-index-overflow",
                    rdkit_function: "Atom::getPerturbationOrder",
                }
            })?);
            m_order.push(i32::try_from(m_bond.id().index()).map_err(|_| {
                SubstructMatchError::Unsupported {
                    branch: "MolMatchFinalCheckFunctor/mOrder/bond-index-overflow",
                    rdkit_function: "countSwapsToInterconvert",
                }
            })?);
            if m_order.len() == query.adjacency().get(qi).map_or(0, Vec::len) {
                break;
            }
        }
        if q_order.len() != query.adjacency().get(qi).map_or(0, Vec::len)
            || q_order.len() != m_order.len()
        {
            return Err(SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/chiral-atom-missing-matched-neighbors",
                rdkit_function: "MolMatchFinalCheckFunctor::operator()",
            });
        }
        // RDKit✔️✔️:     int qPermCount = qAt->getPerturbationOrder(qOrder);
        let q_perm_count = rdkit_atom_perturbation_order_from_bond_indices(query, qi, &q_order)?;

        // RDKit✔️✔️:     unsigned unmatchedNeighbors = mAt->getDegree() - mOrder.size();
        // RDKit✔️✔️:     mOrder.insert(mOrder.end(), unmatchedNeighbors, -1);
        let unmatched_neighbors = mol
            .topology_block()
            .adjacency
            .neighbors_of(mi)
            .len()
            .saturating_sub(m_order.len());
        m_order.extend(std::iter::repeat_n(-1, unmatched_neighbors));

        // RDKit✔️✔️:     INT_LIST moOrder;
        // RDKit✔️✔️:     for (const auto &bond : d_mol.atomBonds(mAt)) {
        // RDKit✔️✔️:       const int dbidx = bond->getIdx();
        // RDKit✔️✔️:       if (std::find(mOrder.begin(), mOrder.end(), dbidx) != mOrder.end()) {
        // RDKit✔️✔️:         moOrder.push_back(dbidx);
        // RDKit✔️✔️:       } else {
        // RDKit✔️✔️:         moOrder.push_back(-1);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mo_order: Vec<i32> = mol
            .topology_block()
            .adjacency
            .neighbors_of(mi)
            .iter()
            .map(|neighbor| {
                i32::try_from(neighbor.bond.index()).map(|bond_idx| {
                    if m_order.contains(&bond_idx) {
                        bond_idx
                    } else {
                        -1
                    }
                })
            })
            .collect::<Result<_, _>>()
            .map_err(|_| SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/moOrder/bond-index-overflow",
                rdkit_function: "countSwapsToInterconvert",
            })?;
        // RDKit✔️✔️:     const int mPermCount =
        // RDKit✔️✔️:         static_cast<int>(countSwapsToInterconvert(moOrder, mOrder));
        let m_perm_count = count_swaps_to_interconvert_i32(&mo_order, &m_order).ok_or(
            SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/mPermCount/unmodeled-bond-ordering",
                rdkit_function: "countSwapsToInterconvert",
            },
        )?;

        // RDKit✔️✔️:     const bool requireMatch = qPermCount % 2 == mPermCount % 2;
        // RDKit✔️✔️:     const bool labelsMatch = qAt->getChiralTag() == mAt->getChiralTag();
        // RDKit✔️✔️:     const bool matchOK = requireMatch == labelsMatch;
        // RDKit✔️✔️:     // if this is not part of a stereogroup and doesn't match, return false
        // RDKit✔️✔️:     const auto msg = d_molStereoGroups.find(m_c[i]);
        // RDKit✔️✔️:     if (msg == d_molStereoGroups.end()) {
        // RDKit✔️✔️:       if (!matchOK) {
        // RDKit✔️✔️:         return false;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       matches[m_c[i]] = matchOK;
        // RDKit✔️✔️:     }
        let require_match = q_perm_count % 2 == m_perm_count % 2;
        let labels_match = q_at.chiral_tag() == m_at.chiral_tag();
        let match_ok = require_match == labels_match;
        if mol_stereo_groups[mi].is_some() {
            stereo_matches[mi] = Some(match_ok);
        } else if !match_ok {
            return Ok(false);
        }
    }

    // RDKit✔️✔️:   std::unordered_map<unsigned int, unsigned int> q_to_mol;
    // RDKit✔️✔️:   for (unsigned int j = 0; j < d_query.getNumAtoms(); ++j) {
    // RDKit✔️✔️:     q_to_mol[q_c[j]] = m_c[j];
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (d_params.useEnhancedStereo) {
    // RDKit✔️✔️:     if (!detail::enhancedStereoIsOK(d_mol, d_query, q_to_mol, d_molStereoGroups,
    // RDKit✔️✔️:                                     matches)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if params.use_enhanced_stereo
        && !enhanced_stereo_is_ok(mol, query, &q_to_mol, mol_stereo_groups, &stereo_matches)
    {
        return Ok(false);
    }

    // RDKit✔️✔️:   // now check double bonds
    // RDKit✔️✔️:   for (const auto &qBnd : d_query.bonds()) {
    // RDKit✔️✔️:     if (qBnd->getBondType() != Bond::DOUBLE ||
    // RDKit✔️✔️:         qBnd->getStereo() <= Bond::STEREOANY) {
    // RDKit✔️✔️:       continue;
    // RDKit✔️✔️:     }
    for q_bnd in query.bonds() {
        if q_bnd.bond().order() != BondOrder::Double
            || !rdkit_bond_stereo_is_above_any(q_bnd.bond().stereo())
        {
            continue;
        }
        // RDKit✔️✔️:     if (qBnd->getStereoAtoms().size() != 2) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(q_stereo_atoms) = q_bnd.bond().stereo_atoms() else {
            continue;
        };
        // RDKit✔️✔️:     const Bond *mBnd = d_mol.getBondBetweenAtoms(
        // RDKit✔️✔️:         q_to_mol[qBnd->getBeginAtomIdx()], q_to_mol[qBnd->getEndAtomIdx()]);
        let q_begin_mol = q_to_mol[q_bnd.begin().index()];
        let q_end_mol = q_to_mol[q_bnd.end().index()];
        let Some(m_bnd) = find_bond_between(mol, q_begin_mol, q_end_mol) else {
            return Err(SubstructMatchError::Unsupported {
                branch: "MolMatchFinalCheckFunctor/double-bond-matching-bond-missing",
                rdkit_function: "MolMatchFinalCheckFunctor::operator()",
            });
        };
        // RDKit✔️✔️:     if (mBnd->getBondType() != Bond::DOUBLE) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        if m_bnd.order() != BondOrder::Double {
            continue;
        }
        // RDKit✔️✔️:     if (!d_params.specifiedStereoQueryMatchesUnspecified &&
        // RDKit✔️✔️:         mBnd->getStereo() <= Bond::STEREOANY) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if !params.specified_stereo_query_matches_unspecified
            && !rdkit_bond_stereo_is_above_any(m_bnd.stereo())
        {
            return Ok(false);
        }
        // RDKit✔️✔️:     if (mBnd->getStereoAtoms().size() != 2) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(m_stereo_atoms) = m_bnd.stereo_atoms() else {
            continue;
        };

        // RDKit✔️✔️:     unsigned int end1Matches = 0;
        // RDKit✔️✔️:     unsigned int end2Matches = 0;
        // RDKit✔️✔️:     if (q_to_mol[qBnd->getBeginAtomIdx()] == mBnd->getBeginAtomIdx()) {
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        // RDKit✔️✔️:         end1Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        // RDKit✔️✔️:         end2Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     } else {
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[0]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[1])) {
        // RDKit✔️✔️:         end1Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       if (q_to_mol[qBnd->getStereoAtoms()[1]] ==
        // RDKit✔️✔️:           static_cast<unsigned>(mBnd->getStereoAtoms()[0])) {
        // RDKit✔️✔️:         end2Matches = 1;
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        let mut end1_matches = 0_u32;
        let mut end2_matches = 0_u32;
        if q_begin_mol == m_bnd.begin().index() {
            if q_to_mol[q_stereo_atoms[0].index()] == m_stereo_atoms[0].index() {
                end1_matches = 1;
            }
            if q_to_mol[q_stereo_atoms[1].index()] == m_stereo_atoms[1].index() {
                end2_matches = 1;
            }
        } else {
            if q_to_mol[q_stereo_atoms[0].index()] == m_stereo_atoms[1].index() {
                end1_matches = 1;
            }
            if q_to_mol[q_stereo_atoms[1].index()] == m_stereo_atoms[0].index() {
                end2_matches = 1;
            }
        }

        // RDKit✔️✔️:     const unsigned totalMatches = end1Matches + end2Matches;
        // RDKit✔️✔️:     const auto mStereo =
        // RDKit✔️✔️:         Chirality::translateEZLabelToCisTrans(mBnd->getStereo());
        // RDKit✔️✔️:     const auto qStereo =
        // RDKit✔️✔️:         Chirality::translateEZLabelToCisTrans(qBnd->getStereo());
        // RDKit✔️✔️:     if (mStereo == qStereo && totalMatches == 1) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (mStereo != qStereo && totalMatches != 1) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        let total_matches = end1_matches + end2_matches;
        let m_stereo = rdkit_translate_ez_label_to_cis_trans(m_bnd.stereo());
        let q_stereo = rdkit_translate_ez_label_to_cis_trans(q_bnd.bond().stereo());
        if m_stereo == q_stereo && total_matches == 1 {
            return Ok(false);
        }
        if m_stereo != q_stereo && total_matches != 1 {
            return Ok(false);
        }
    }

    // RDKit✔️✔️:   if (d_params.uniquify) {
    // RDKit✔️✔️:     matchesSeen.insert(match);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return true;
    if let Some(mask) = match_key {
        matches_seen.insert(mask);
    }
    Ok(true)
}

// ---------------------------------------------------------------------------
// Bond mapping builder
// ---------------------------------------------------------------------------

/// Build the bond mapping for a match result.
///
/// For each query bond (by index), find the corresponding molecular bond
/// that connects the matched query endpoints.
#[allow(dead_code)]
fn build_bond_mapping(
    query_atom_to_mol: &[Option<usize>],
    query: &Vf2Graph,
    mol: &Vf2Graph,
) -> Vec<usize> {
    let mut bond_mapping = Vec::with_capacity(query.n_bonds);
    for bond_idx in 0..query.n_bonds {
        // Find the query atoms connected by this bond.
        let mut q_begin = NULL_NODE;
        let mut q_end = NULL_NODE;
        for qa in 0..query.n_atoms {
            for &(nbr, eidx) in &query.adjacency[qa] {
                if eidx == bond_idx {
                    q_begin = qa;
                    q_end = nbr;
                    break;
                }
            }
            if q_begin != NULL_NODE {
                break;
            }
        }

        if q_begin != NULL_NODE {
            let m_begin = query_atom_to_mol[q_begin];
            let m_end = query_atom_to_mol[q_end];
            if let (Some(mb), Some(me)) = (m_begin, m_end) {
                // Find bond between mb and me in mol.
                let mut mol_bond_idx = NULL_NODE;
                for &(nbr, eidx) in &mol.adjacency[mb] {
                    if nbr == me {
                        mol_bond_idx = eidx;
                        break;
                    }
                }
                bond_mapping.push(mol_bond_idx);
            } else {
                bond_mapping.push(NULL_NODE);
            }
        } else {
            bond_mapping.push(NULL_NODE);
        }
    }
    bond_mapping
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

fn preflight_atom_query(
    query: &crate::QueryNode<AtomQueryPredicate>,
) -> Result<(), SubstructMatchError> {
    match query {
        crate::QueryNode::Predicate(AtomQueryPredicate::UnsupportedFeature(branch)) => {
            Err(SubstructMatchError::Unsupported {
                branch,
                rdkit_function: "QueryAtom::Match",
            })
        }
        crate::QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => {
            let inner_query = query
                .query_graph()
                .ok_or(SubstructMatchError::Unsupported {
                    branch: "recursive SMARTS without a compiled query molecule",
                    rdkit_function: "RecursiveStructureQuery::getQueryMol",
                })?;
            preflight_query_molecule(inner_query)
        }
        crate::QueryNode::Predicate(_) => Ok(()),
        crate::QueryNode::And(children)
        | crate::QueryNode::Or(children)
        | crate::QueryNode::Xor(children) => {
            for child in children {
                preflight_atom_query(child)?;
            }
            Ok(())
        }
        crate::QueryNode::Not(child) => preflight_atom_query(child),
    }
}

fn preflight_bond_query(
    query: &crate::QueryNode<BondQueryPredicate>,
) -> Result<(), SubstructMatchError> {
    match query {
        crate::QueryNode::Predicate(BondQueryPredicate::UnsupportedFeature(branch)) => {
            Err(SubstructMatchError::Unsupported {
                branch,
                rdkit_function: "QueryBond::Match",
            })
        }
        crate::QueryNode::Predicate(_) => Ok(()),
        crate::QueryNode::And(children)
        | crate::QueryNode::Or(children)
        | crate::QueryNode::Xor(children) => {
            for child in children {
                preflight_bond_query(child)?;
            }
            Ok(())
        }
        crate::QueryNode::Not(child) => preflight_bond_query(child),
    }
}

fn preflight_query_molecule(query: &QueryGraph) -> Result<(), SubstructMatchError> {
    // This fail-closed preflight has no RDKit counterpart: RDKit query leaves
    // are executable, while COSMolKit can preserve explicitly unsupported
    // leaves imported from other formats. Inspecting every leaf before VF2
    // prevents AND/OR short-circuiting from turning unsupported chemistry into
    // a plausible match or mismatch.
    //
    // Local complexity review: this is O(A + B + Q), where Q includes all
    // owned recursive query trees. It allocates no collections and performs no
    // molecule or query clones. Each supported leaf is visited once before the
    // existing matcher traversal; failure returns at the first unsupported
    // leaf.
    for atom in query.atoms() {
        preflight_atom_query(atom.predicate())?;
    }
    for bond in query.bonds() {
        preflight_bond_query(bond.predicate())?;
    }
    Ok(())
}

fn recursive_matcher(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
    query_context: Option<&QueryMatchContext<'_>>,
) -> Result<Vec<bool>, SubstructMatchError> {
    // RDKit✔️❌: unsigned int RecursiveMatcher(const ROMol &mol, const ROMol &query,
    // RDKit✔️❌:                               std::vector<int> &matches,
    // RDKit✔️❌:                               SUBQUERY_MAP &subqueryMap,
    // RDKit✔️❌:                               const SubstructMatchParameters &params,
    // RDKit✔️❌:                               std::vector<RecursiveStructureQuery *> &locked) {
    // RDKit✔️❌:   SubstructMatchParameters lparams = params;
    // RDKit✔️❌:   lparams.maxMatches = std::max(params.maxRecursiveMatches, params.maxMatches);
    // RDKit✔️❌:   lparams.uniquify = false;
    // RDKit✔️❌:   for (auto qAtom : query.atoms()) {
    // RDKit✔️❌:     if (qAtom->hasQuery()) {
    // RDKit✔️❌:       MatchSubqueries(mol, qAtom->getQuery(), lparams, subqueryMap, locked);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   detail::AtomLabelFunctor atomLabeler(query, mol, lparams);
    // RDKit✔️❌:   detail::BondLabelFunctor bondLabeler(query, mol, lparams);
    // RDKit✔️❌:   MolMatchFinalCheckFunctor matchChecker(query, mol, lparams);
    // RDKit✔️❌:
    // RDKit✔️❌:   matches.clear();
    // RDKit✔️❌:   matches.resize(0);
    // RDKit✔️❌:   std::vector<detail::ssPairType> pms;
    // RDKit✔️❌:   bool found =
    // RDKit✔️❌:       boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
    // RDKit✔️❌:                      bondLabeler, matchChecker, pms, lparams.maxMatches);
    // RDKit✔️❌:   unsigned int res = 0;
    // RDKit✔️❌:   if (found) {
    // RDKit✔️❌:     matches.reserve(pms.size());
    // RDKit✔️❌:     for (const auto &pairs : pms) {
    // RDKit✔️❌:       if (!query.hasProp(common_properties::_queryRootAtom)) {
    // RDKit✔️❌:         matches.push_back(pairs.begin()->second);
    // RDKit✔️❌:       } else {
    // RDKit✔️❌:         int rootIdx;
    // RDKit✔️❌:         query.getProp(common_properties::_queryRootAtom, rootIdx);
    // RDKit✔️❌:         bool found = false;
    // RDKit✔️❌:         for (const auto &pairIter : pairs) {
    // RDKit✔️❌:           if (pairIter.first == static_cast<unsigned int>(rootIdx)) {
    // RDKit✔️❌:             matches.push_back(pairIter.second);
    // RDKit✔️❌:             found = true;
    // RDKit✔️❌:             break;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (!found) {
    // RDKit✔️❌:           BOOST_LOG(rdErrorLog)
    // RDKit✔️❌:               << "no match found for queryRootAtom" << std::endl;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (matches.size() == lparams.maxMatches) {
    // RDKit✔️❌:         break;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     res = matches.size();
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: }
    // Complexity review: nested preparation and VF2 follow the source. The
    // membership result is one O(target atoms) bool Vec in place of RDKit's
    // ordered set; the root-only projection avoids constructing unused bond
    // rows while the canonical VF2 goal-mapping allocations remain.
    let mut local_params = params.clone();
    local_params.max_matches = params.max_recursive_matches.max(params.max_matches);
    local_params.uniquify = false;
    for atom in query.atoms() {
        match_subqueries(
            mol,
            atom.predicate(),
            &local_params,
            recursive_cache,
            query_context,
        )?;
    }

    // Recursive queries see the same owning target chemistry; preserve both
    // supplied prepared state and the atom-only projection through every depth.
    let matches = match query_context {
        Some(context) => {
            substruct_match_impl_with_recursive_cache_and_context::<AtomOnlyMatchResultProjection>(
                mol,
                query,
                &local_params,
                Some(recursive_cache),
                context,
                None,
                None,
            )
        }
        None => substruct_match_impl_with_recursive_cache::<AtomOnlyMatchResultProjection>(
            mol,
            query,
            &local_params,
            Some(recursive_cache),
        ),
    }?;
    let root_index = query
        .prop("_queryRootAtom")
        .and_then(|value| value.parse::<i32>().ok())
        .map_or(0, |root_index| root_index as u32 as usize);
    let mut match_starts = vec![false; mol.num_atoms()];
    for matched in matches.into_iter().take(local_params.max_matches) {
        if let Some(&root_atom_idx) = matched.get(root_index)
            && root_atom_idx != NULL_NODE
            && root_atom_idx < match_starts.len()
        {
            match_starts[root_atom_idx] = true;
        }
    }
    Ok(match_starts)
}

fn match_subqueries(
    mol: &SearchTarget<'_>,
    query: &crate::QueryNode<AtomQueryPredicate>,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
    query_context: Option<&QueryMatchContext<'_>>,
) -> Result<(), SubstructMatchError> {
    // RDKit✔️❌: void MatchSubqueries(const ROMol &mol, QueryAtom::QUERYATOM_QUERY *query,
    // RDKit✔️❌:                      const SubstructMatchParameters &params,
    // RDKit✔️❌:                      SUBQUERY_MAP &subqueryMap,
    // RDKit✔️❌:                      std::vector<RecursiveStructureQuery *> &locked) {
    // RDKit✔️❌:   PRECONDITION(query, "bad query");
    // RDKit✔️❌:   if (query->getDescription() == "RecursiveStructure") {
    // RDKit✔️❌:     auto *rsq = (RecursiveStructureQuery *)query;
    // RDKit✔️❌: #ifdef RDK_BUILD_THREADSAFE_SSS
    // RDKit✔️❌:     rsq->d_mutex.lock();
    // RDKit✔️❌: #endif
    // RDKit✔️❌:     locked.push_back(rsq);
    // RDKit✔️❌:     rsq->clear();
    // RDKit✔️❌:     bool matchDone = false;
    // RDKit✔️❌:     if (rsq->getSerialNumber() &&
    // RDKit✔️❌:         subqueryMap.find(rsq->getSerialNumber()) != subqueryMap.end()) {
    // RDKit✔️❌:       matchDone = true;
    // RDKit✔️❌:       auto orsq =
    // RDKit✔️❌:           (const RecursiveStructureQuery *)subqueryMap[rsq->getSerialNumber()];
    // RDKit✔️❌:       for (auto setIter = orsq->beginSet(); setIter != orsq->endSet();
    // RDKit✔️❌:            ++setIter) {
    // RDKit✔️❌:         rsq->insert(*setIter);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:
    // RDKit✔️❌:     if (!matchDone) {
    // RDKit✔️❌:       ROMol const *queryMol = rsq->getQueryMol();
    // RDKit✔️❌:       if (queryMol) {
    // RDKit✔️❌:         std::vector<int> matchStarts;
    // RDKit✔️❌:         unsigned int res = RecursiveMatcher(mol, *queryMol, matchStarts,
    // RDKit✔️❌:                                             subqueryMap, params, locked);
    // RDKit✔️❌:         if (res) {
    // RDKit✔️❌:           for (int &matchStart : matchStarts) {
    // RDKit✔️❌:             rsq->insert(matchStart);
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (rsq->getSerialNumber()) {
    // RDKit✔️❌:         subqueryMap[rsq->getSerialNumber()] = query;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   for (auto childIt = query->beginChildren(); childIt != query->endChildren();
    // RDKit✔️❌:        ++childIt) {
    // RDKit✔️❌:     MatchSubqueries(mol, childIt->get(), params, subqueryMap, locked);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Complexity review: every query node is visited once unless a serial-key
    // cache hit skips recursive VF2, matching RDKit. BTreeMap lookup is O(log
    // R) instead of unordered-map average O(1), but recursive VF2 dominates;
    // query trees and match vectors are never cloned here.
    match query {
        crate::QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive_query)) => {
            let cache_key = recursive_query_cache_key(recursive_query);
            if !recursive_cache.contains_key(&cache_key) {
                let match_starts = match recursive_query.query_graph() {
                    Some(inner_query) => {
                        recursive_matcher(mol, inner_query, params, recursive_cache, query_context)?
                    }
                    None => vec![false; mol.num_atoms()],
                };
                recursive_cache.insert(cache_key, match_starts);
            }
        }
        crate::QueryNode::Predicate(_) => {}
        crate::QueryNode::And(children)
        | crate::QueryNode::Or(children)
        | crate::QueryNode::Xor(children) => {
            for child in children {
                match_subqueries(mol, child, params, recursive_cache, query_context)?;
            }
        }
        crate::QueryNode::Not(child) => {
            match_subqueries(mol, child, params, recursive_cache, query_context)?;
        }
    }
    Ok(())
}

fn populate_recursive_query_match_cache(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
    query_context: Option<&QueryMatchContext<'_>>,
) -> Result<(), SubstructMatchError> {
    for atom in query.atoms() {
        match_subqueries(
            mol,
            atom.predicate(),
            params,
            recursive_cache,
            query_context,
        )?;
    }
    Ok(())
}

fn substruct_match_impl_with_recursive_cache<P: MatchResultProjection>(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
) -> Result<Vec<P::Output>, SubstructMatchError> {
    let query_ctx = build_query_match_context(mol);
    substruct_match_impl_with_recursive_cache_and_context::<P>(
        mol,
        query,
        params,
        recursive_cache,
        &query_ctx,
        None,
        None,
    )
}

trait MatchResultProjection {
    type Output;
    fn project(
        query: &QueryGraph,
        target_graph: Vf2GraphRef<'_>,
        atom_mapping: Vec<usize>,
    ) -> Self::Output;
}
struct FullMatchResultProjection;
impl MatchResultProjection for FullMatchResultProjection {
    type Output = SubstructMatchResult;
    fn project(
        query: &QueryGraph,
        target_graph: Vf2GraphRef<'_>,
        atom_mapping: Vec<usize>,
    ) -> Self::Output {
        let bond_mapping = materialize_bond_mapping(query, target_graph, &atom_mapping);
        SubstructMatchResult {
            atom_mapping,
            bond_mapping,
        }
    }
}
struct AtomOnlyMatchResultProjection;
impl MatchResultProjection for AtomOnlyMatchResultProjection {
    type Output = Vec<usize>;
    fn project(
        _query: &QueryGraph,
        _target_graph: Vf2GraphRef<'_>,
        atom_mapping: Vec<usize>,
    ) -> Self::Output {
        atom_mapping
    }
}
#[cfg(test)]
thread_local! {
    static BOND_MAPPING_MATERIALIZATION_COUNT: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
}
fn materialize_bond_mapping(
    query: &QueryGraph,
    target_graph: Vf2GraphRef<'_>,
    atom_mapping: &[usize],
) -> Vec<usize> {
    // RDKit atom matches are projected above; the existing public full result also requires target bond indices.
    #[cfg(test)]
    BOND_MAPPING_MATERIALIZATION_COUNT.with(|count| count.set(count.get() + 1));
    let mut bond_mapping = Vec::with_capacity(query.num_bonds());
    for qbond in query.bonds() {
        let q_begin = qbond.begin().index();
        let q_end = qbond.end().index();
        let m_begin = atom_mapping[q_begin];
        let m_end = atom_mapping[q_end];
        if m_begin != NULL_NODE && m_end != NULL_NODE {
            let found = target_graph
                .neighbor_row(m_begin)
                .iter()
                .find(|(neighbor, _)| *neighbor == m_end)
                .map(|(_, edge)| edge);
            bond_mapping.push(found.unwrap_or(NULL_NODE));
        } else {
            bond_mapping.push(NULL_NODE);
        }
    }

    bond_mapping
}

fn substruct_match_impl_with_recursive_cache_and_context<P: MatchResultProjection>(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
    query_order: Option<&[usize]>,
    compiled_graph: Option<&CompiledQueryGraph>,
) -> Result<Vec<P::Output>, SubstructMatchError> {
    let m_num_atoms = mol.num_atoms();
    let q_num_atoms = query.num_atoms();

    // RDKit source (SubstructMatch.cpp):
    //   if (!mNumAtoms || !qNumAtoms || qNumAtoms > mNumAtoms) {
    //     return matches;
    //   }
    if m_num_atoms == 0 || q_num_atoms == 0 || q_num_atoms > m_num_atoms {
        return Ok(Vec::new());
    }

    // RDKit's vf2_all receives the existing topology objects. The ordinary
    // path borrows QueryGraph and SearchTarget adjacency; only a compiled
    // query uses its already-retained owned graph. No graph is rebuilt here.
    let q_graph = match compiled_graph {
        Some(compiled_graph) => Vf2GraphRef::compiled(compiled_graph),
        None => Vf2GraphRef::query(query),
    };
    let m_graph = Vf2GraphRef::target(mol.topology_block());

    // Build atom matching closure.
    // RDKit source:
    //   detail::AtomLabelFunctor atomLabeler(query, mol, params);
    //   detail::BondLabelFunctor bondLabeler(query, mol, params);
    //   MolMatchFinalCheckFunctor matchChecker(query, mol, params);
    let atom_match_error = std::cell::Cell::new(None);
    let atom_fn = |qi: usize, mj: usize| -> bool {
        match atom_label_matches(query, mol, qi, mj, params, recursive_cache, query_ctx) {
            Ok(matched) => matched,
            Err(error) => {
                if atom_match_error.get().is_none() {
                    atom_match_error.set(Some(error));
                }
                false
            }
        }
    };

    let bond_fn = |qei: usize, mei: usize| -> bool {
        bond_label_matches(query, mol, qei, mei, params, query_ctx)
    };

    // RDKit source:
    //   bool found = boost::vf2_all(query.getTopology(), mol.getTopology(),
    //                               atomLabeler, bondLabeler, matchChecker,
    //                               pms, params.maxMatches);
    let mut raw_matches: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
    let mut matches_seen: HashSet<Vec<bool>> = HashSet::new();
    let final_check_setup = MolMatchFinalCheckSetup::new(query, mol, params);
    let mut final_check_error: Option<SubstructMatchError> = None;
    let mut check_fn = |c1: &[NodeId], c2: &[NodeId]| -> bool {
        match rdkit_match_final_check(
            mol,
            query,
            params,
            c1,
            c2,
            &final_check_setup,
            &mut matches_seen,
        ) {
            Ok(accepted) => accepted,
            Err(err) => {
                final_check_error = Some(err);
                false
            }
        }
    };

    vf2_entry_all_ordered(
        q_graph,
        m_graph,
        &atom_fn,
        &bond_fn,
        Some(&mut check_fn),
        &mut raw_matches,
        params.max_matches,
        query_order,
    );
    if let Some(error) = atom_match_error.get() {
        return Err(error.into());
    }
    if let Some(err) = final_check_error {
        return Err(err);
    }

    Ok(project_match_results::<P>(query, m_graph, &raw_matches))
}
fn project_match_results<P: MatchResultProjection>(
    query: &QueryGraph,
    target_graph: Vf2GraphRef<'_>,
    raw_matches: &[Vec<(NodeId, NodeId)>],
) -> Vec<P::Output> {
    // RDKit source: third_party/rdkit/Code/GraphMol/Substruct/SubstructMatch.cpp
    // RDKit✔️❌:   if (found) {
    // RDKit✔️❌:     const unsigned int nQueryAtoms = query.getNumAtoms();
    // RDKit✔️❌:     matches.reserve(pms.size());
    // RDKit✔️❌:     MatchVectType matchVect(nQueryAtoms);
    // RDKit✔️❌:     for (const auto &pairs : pms) {
    // RDKit✔️❌:       for (const auto &pair : pairs) {
    // RDKit✔️❌:         matchVect[pair.first] = pair;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       matches.push_back(matchVect);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // Behavior review: source-produced pairs populate their query-index slot,
    // each accepted row remains in order, and the prior bounds guard and
    // NULL_NODE output are preserved. P01 and the frozen S06/S08 routes pass.
    // Complexity review: P01 verifies one reserved output outer vector, one
    // direct final atom map per result, no reallocations, and no conversion
    // scratch/copy. RDKit's MatchVectType is atom-only; this result also needs
    // a bond map and incident-neighbor search per query bond. That required
    // extra allocation/search remains source-relative overhead; this marker
    // does not claim global matcher or projection allocation parity.
    let q_num_atoms = query.num_atoms();
    let mut results = Vec::with_capacity(raw_matches.len());

    for pairs in raw_matches {
        let mut atom_mapping = vec![NULL_NODE; q_num_atoms];
        for &(qa, ma) in pairs {
            if qa < q_num_atoms {
                atom_mapping[qa] = ma;
            }
        }

        results.push(P::project(query, target_graph, atom_mapping));
    }

    results
}
#[cfg(test)]
fn project_substruct_matches(
    query: &QueryGraph,
    target_graph: Vf2GraphRef<'_>,
    raw_matches: &[Vec<(NodeId, NodeId)>],
) -> Vec<SubstructMatchResult> {
    project_match_results::<FullMatchResultProjection>(query, target_graph, raw_matches)
}

fn substruct_matches_with_compiled_query<P: MatchResultProjection>(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    compiled_graph: &CompiledQueryGraph,
) -> Result<Vec<P::Output>, SubstructMatchError> {
    substruct_matches_with_compiled_query_and_context::<P>(mol, query, params, compiled_graph, None)
}

fn substruct_matches_with_compiled_query_and_context<P: MatchResultProjection>(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    compiled_graph: &CompiledQueryGraph,
    query_context: Option<&QueryMatchContext>,
) -> Result<Vec<P::Output>, SubstructMatchError> {
    preflight_query_molecule(query)?;
    if mol.num_atoms() == 0 || query.num_atoms() == 0 || query.num_atoms() > mol.num_atoms() {
        return Ok(Vec::new());
    }
    let mut recursive_locker = RecursiveLocker::new(query, params.recursion_possible);
    if params.recursion_possible {
        populate_recursive_query_match_cache(
            mol,
            query,
            params,
            &mut recursive_locker.cache,
            query_context,
        )?;
    }
    let owned_query_context;
    let query_ctx = match query_context {
        Some(query_context) => query_context,
        None => {
            owned_query_context = build_query_match_context(mol);
            &owned_query_context
        }
    };
    // RDKit's `vf2_all` creates its initial state with `sortNodes=false`.
    // Keep source graph order for deterministic enumeration even though the
    // compiled plan retains its separate atom-order metadata.
    substruct_match_impl_with_recursive_cache_and_context::<P>(
        mol,
        query,
        params,
        Some(&recursive_locker.cache),
        &query_ctx,
        None,
        Some(compiled_graph),
    )
}

#[allow(dead_code)] // Made cross-crate in the authorized Step 10 export.
fn substruct_atom_matches_with_compiled_query(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    compiled_graph: &CompiledQueryGraph,
) -> Result<Vec<Vec<usize>>, SubstructMatchError> {
    substruct_matches_with_compiled_query::<AtomOnlyMatchResultProjection>(
        mol,
        query,
        params,
        compiled_graph,
    )
}

pub(crate) fn try_get_substruct_atom_matches_with_compiled_query_and_context(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    compiled_graph: &CompiledQueryGraph,
    query_context: &QueryMatchContext,
) -> Result<Vec<Vec<usize>>, SubstructMatchError> {
    substruct_matches_with_compiled_query_and_context::<AtomOnlyMatchResultProjection>(
        mol,
        query,
        params,
        compiled_graph,
        Some(query_context),
    )
}

fn atom_compat(
    query_atom: &QueryAtom,
    query_mol: &QueryGraph,
    mol_atom: &Atom,
    mol: &SearchTarget<'_>,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> Result<bool, PeriodicTableError> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: atomCompat
    // RDKit✔️✔️: bool atomCompat(const Atom *a1, const Atom *a2,
    // RDKit✔️✔️:                 const SubstructMatchParameters &ps) {
    // RDKit✔️✔️:   PRECONDITION(a1, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(a2, "bad atom");
    // RDKit✔️✔️:   // std::cerr << "\t\tatomCompat: "<< a1 << " " << a1->getIdx() << "-" << a2 <<
    // RDKit✔️✔️:   // " " << a2->getIdx() << std::endl;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ps.extraAtomCheckOverridesDefaultCheck && ps.extraAtomCheck) {
    // RDKit✔️✔️:     return ps.extraAtomCheck(*a1, *a2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool res;
    // RDKit✔️✔️:   if (ps.useQueryQueryMatches && a1->hasQuery() && a2->hasQuery()) {
    // RDKit✔️✔️:     res = static_cast<const QueryAtom *>(a1)->QueryMatch(
    // RDKit✔️✔️:         static_cast<const QueryAtom *>(a2));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = a1->Match(a2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!res) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!ps.atomProperties.empty()) {
    // RDKit✔️✔️:     if (!propertyCompat(a1, a2, ps.atomProperties)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (ps.extraAtomCheck && !ps.extraAtomCheck(*a1, *a2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Typed references make the source pointer preconditions
    // unrepresentable. Local complexity review: both implementations perform
    // the same constant-time origin/option dispatch, one default atom or query
    // match, the requested property scan, and at most one callback invocation
    // after the default match. Arc callback dispatch is the Rust equivalent of
    // std::function dispatch and allocates nothing per match. Query-tree
    // traversal and recursive-cache lookup retain their existing complexity;
    // property_compat has the separately documented BTreeMap improvement. No
    // atom, molecule, query, property map, or callback is cloned in this hot
    // path.
    if params.extra_atom_check_overrides_default_check
        && let Some(extra_atom_check) = &params.extra_atom_check
    {
        return Ok(extra_atom_check(query_mol, query_atom, mol, mol_atom));
    }

    let target_query = (params.use_query_query_matches
        && !query_atom.predicate_is_carrier_derived()
        && mol.atom_has_query(mol_atom.id()))
    .then(|| mol.atom_query_predicate(mol_atom.id()))
    .flatten();
    let matches = if let Some(target_query) = target_query {
        query_atom_query_match(
            query_atom.predicate(),
            Some(target_query),
            mol_atom,
            mol,
            query_ctx,
        )?
    } else if query_atom.predicate_is_carrier_derived() {
        atom_matches(query_atom, mol_atom, mol)
    } else {
        // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryAtom.cpp :: QueryAtom::Match
        // RDKit❗✔️: bool QueryAtom::Match(Atom const *what) const {
        // RDKit❗✔️:   PRECONDITION(what, "bad query atom");
        // RDKit❗✔️:   PRECONDITION(dp_query, "no query set");
        // RDKit❗✔️:   return dp_query->Match(what);
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION
        // The existing predicate evaluator applies that explicit query tree
        // to this ordinary target atom. Carrier-derived rows use Atom::Match
        // above.
        evaluate_atom_query(
            query_atom.predicate(),
            mol_atom,
            mol,
            params,
            recursive_cache,
            query_ctx,
        )?
    };
    if !matches {
        return Ok(false);
    }
    if !params.atom_properties.is_empty()
        && !property_compat(
            query_atom.props(),
            mol_atom.props(),
            &params.atom_properties,
        )
    {
        return Ok(false);
    }
    if let Some(extra_atom_check) = &params.extra_atom_check
        && !extra_atom_check(query_mol, query_atom, mol, mol_atom)
    {
        return Ok(false);
    }
    Ok(matches)
}

#[allow(deprecated)]
fn chiral_atom_compat(
    query_atom: &QueryAtom,
    _query_mol: &QueryGraph,
    mol_atom: &Atom,
    mol: &SearchTarget<'_>,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: chiralAtomCompat
    // RDKit✔️✔️: bool chiralAtomCompat(const Atom *&a1, const Atom *&a2) {
    // RDKit✔️✔️:   /// DEPRECATED
    // RDKit✔️✔️:   PRECONDITION(a1, "bad atom");
    // RDKit✔️✔️:   PRECONDITION(a2, "bad atom");
    // RDKit✔️✔️:   bool res = a1->Match(a2);
    // RDKit✔️✔️:   if (res) {
    // RDKit✔️✔️:     std::string s1, s2;
    // RDKit✔️✔️:     bool hascode1 = a1->getPropIfPresent(common_properties::_CIPCode, s1);
    // RDKit✔️✔️:     bool hascode2 = a2->getPropIfPresent(common_properties::_CIPCode, s2);
    // RDKit✔️✔️:     if (hascode1 || hascode2) {
    // RDKit✔️✔️:       res = hascode1 && hascode2 && s1 == s2;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::cerr << "\t\tchiralAtomCompat: " << a1 << " " << a1->getIdx() << "-"
    // RDKit✔️✔️:             << a2 << " " << a2->getIdx() << std::endl;
    // RDKit✔️✔️:   std::cerr << "\t\t    " << res << std::endl;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Rust references make both pointer preconditions unrepresentable. Local
    // complexity review: the shared atom matcher has the same source-defined
    // atom-query cost, followed by two property lookups and one string
    // comparison only after a successful atom match. No molecule, atom,
    // property map, or string is cloned. BTreeMap lookup retains the canonical
    // atom property representation and has the same logarithmic lookup class
    // as RDKit's property dictionary for the modeled state.
    let mut matches = atom_matches(query_atom, mol_atom, mol);
    if matches {
        let query_cip = query_atom.prop("_CIPCode");
        let mol_cip = mol_atom.prop("_CIPCode");
        if query_cip.is_some() || mol_cip.is_some() {
            matches = query_cip.is_some() && mol_cip.is_some() && query_cip == mol_cip;
        }
    }
    eprintln!(
        "\t\tchiralAtomCompat: {:p} {}-{:p} {}",
        query_atom,
        query_atom.id().index(),
        mol_atom,
        mol_atom.id().index()
    );
    eprintln!("\t\t    {}", u8::from(matches));
    matches
}

fn bond_compat(
    query_bond: &QueryBond,
    query_mol: &QueryGraph,
    mol_bond: &Bond,
    mol: &SearchTarget<'_>,
    params: &SubstructMatchParams,
    query_ctx: &QueryMatchContext,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: bondCompat
    // RDKit✔️✔️: bool bondCompat(const Bond *b1, const Bond *b2,
    // RDKit✔️✔️:                 const SubstructMatchParameters &ps) {
    // RDKit✔️✔️:   PRECONDITION(b1, "bad bond");
    // RDKit✔️✔️:   PRECONDITION(b2, "bad bond");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ps.extraBondCheckOverridesDefaultCheck && ps.extraBondCheck) {
    // RDKit✔️✔️:     return ps.extraBondCheck(*b1, *b2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto isConjugatedSingleOrDoubleBond([](const Bond *bond) {
    // RDKit✔️✔️:     return bond->getIsConjugated() && (bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                                        bond->getBondType() == Bond::DOUBLE);
    // RDKit✔️✔️:   });
    // RDKit✔️✔️:   auto isSingleOrDoubleBond([](const Bond *bond) {
    // RDKit✔️✔️:     return (bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:             bond->getBondType() == Bond::DOUBLE);
    // RDKit✔️✔️:   });
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (ps.useQueryQueryMatches && b1->hasQuery() && b2->hasQuery()) {
    // RDKit✔️✔️:     res = static_cast<const QueryBond *>(b1)->QueryMatch(
    // RDKit✔️✔️:         static_cast<const QueryBond *>(b2));
    // RDKit✔️✔️:   } else if (ps.aromaticMatchesConjugated && !b1->hasQuery() &&
    // RDKit✔️✔️:              !b2->hasQuery() &&
    // RDKit✔️✔️:              ((b1->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                b2->getBondType() == Bond::AROMATIC) ||
    // RDKit✔️✔️:               (b1->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                isConjugatedSingleOrDoubleBond(b2)) ||
    // RDKit✔️✔️:               (b2->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                isConjugatedSingleOrDoubleBond(b1)))) {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:   } else if (ps.aromaticMatchesSingleOrDouble && !b1->hasQuery() &&
    // RDKit✔️✔️:              !b2->hasQuery() &&
    // RDKit✔️✔️:              ((b1->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                b2->getBondType() == Bond::AROMATIC) ||
    // RDKit✔️✔️:               (b1->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                isSingleOrDoubleBond(b2)) ||
    // RDKit✔️✔️:               (b2->getBondType() == Bond::AROMATIC &&
    // RDKit✔️✔️:                isSingleOrDoubleBond(b1)))) {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     res = b1->Match(b2);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!res) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (b1->getBondType() == Bond::DATIVE && b2->getBondType() == Bond::DATIVE) {
    // RDKit✔️✔️:     // for dative bonds we need to make sure that the direction also matches:
    // RDKit✔️✔️:     if (!b1->getBeginAtom()->Match(b2->getBeginAtom()) ||
    // RDKit✔️✔️:         !b1->getEndAtom()->Match(b2->getEndAtom())) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (!ps.bondProperties.empty()) {
    // RDKit✔️✔️:     if (!propertyCompat(b1, b2, ps.bondProperties)) {
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (ps.extraBondCheck && !ps.extraBondCheck(*b1, *b2)) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Rust references make both pointer preconditions unrepresentable. Local
    // complexity review: flag/order checks and dative endpoint lookups are
    // constant time; query-tree matching reuses the canonical evaluator with
    // its source-equivalent tree complexity. The property scan is linear in
    // the requested names with logarithmic canonical BTreeMap lookup, and at
    // most one Arc callback dispatch occurs. Nothing is cloned or allocated.
    if params.extra_bond_check_overrides_default_check
        && let Some(extra_bond_check) = &params.extra_bond_check
    {
        return extra_bond_check(query_bond.bond(), mol_bond);
    }

    let is_conjugated_single_or_double = |bond: &Bond| {
        bond.is_conjugated() && matches!(bond.order(), BondOrder::Single | BondOrder::Double)
    };
    let is_single_or_double =
        |bond: &Bond| matches!(bond.order(), BondOrder::Single | BondOrder::Double);
    let aromatic_pair_matches = |other_matches: &dyn Fn(&Bond) -> bool| {
        (query_bond.bond().order() == BondOrder::Aromatic
            && mol_bond.order() == BondOrder::Aromatic)
            || (query_bond.bond().order() == BondOrder::Aromatic && other_matches(mol_bond))
            || (mol_bond.order() == BondOrder::Aromatic && other_matches(query_bond.bond()))
    };

    let query_has_query = !query_bond.predicate_is_carrier_derived();
    let target_has_query = mol.bond_has_query(mol_bond.id());
    let matches = if params.use_query_query_matches && query_has_query && target_has_query {
        let target_query = mol
            .bond_query_predicate(mol_bond.id())
            .expect("validated query target exposes each explicit bond predicate");
        query_bond_query_match(
            query_bond.predicate(),
            Some(target_query),
            mol_bond,
            mol,
            query_ctx,
        )
    } else if params.aromatic_matches_conjugated
        && !query_has_query
        && !target_has_query
        && aromatic_pair_matches(&is_conjugated_single_or_double)
    {
        true
    } else if params.aromatic_matches_single_or_double
        && !query_has_query
        && !target_has_query
        && aromatic_pair_matches(&is_single_or_double)
    {
        true
    } else {
        evaluate_bond_query(query_bond.predicate(), mol_bond, mol, query_ctx)
    };
    if !matches {
        return false;
    }

    if query_bond.bond().order() == BondOrder::Dative && mol_bond.order() == BondOrder::Dative {
        let query_begin = &query_mol.atoms()[query_bond.begin().index()];
        let query_end = &query_mol.atoms()[query_bond.end().index()];
        let mol_begin = &mol.atoms()[mol_bond.begin().index()];
        let mol_end = &mol.atoms()[mol_bond.end().index()];
        if !atom_matches(query_begin, mol_begin, mol) || !atom_matches(query_end, mol_end, mol) {
            return false;
        }
    }
    if !params.bond_properties.is_empty()
        && !property_compat(
            query_bond.bond().props(),
            mol_bond.props(),
            &params.bond_properties,
        )
    {
        return false;
    }
    if let Some(extra_bond_check) = &params.extra_bond_check
        && !extra_bond_check(query_bond.bond(), mol_bond)
    {
        return false;
    }
    matches
}

fn remove_duplicates(matches: &mut Vec<SubstructMatchResult>, atom_count: usize) {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: removeDuplicates
    // RDKit✔️✔️: void removeDuplicates(std::vector<MatchVectType> &matches,
    // RDKit✔️✔️:                       unsigned int nAtoms) {
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //  This works by tracking the indices of the atoms in each match vector.
    // RDKit✔️✔️:   //  This can lead to unexpected behavior when looking at rings and queries
    // RDKit✔️✔️:   //  that don't specify bond orders.  For example querying this molecule:
    // RDKit✔️✔️:   //    C1CCC=1
    // RDKit✔️✔️:   //  with the pattern constructed from SMARTS C~C~C~C will return a
    // RDKit✔️✔️:   //  single match, despite the fact that there are 4 different paths
    // RDKit✔️✔️:   //  when valence is considered.  The defense of this behavior is
    // RDKit✔️✔️:   //  that the 4 paths are equivalent in the semantics of the query.
    // RDKit✔️✔️:   //  Also, OELib returns the same results
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   std::unordered_set<std::string> seen;
    // RDKit✔️✔️:   std::vector<MatchVectType> res;
    // RDKit✔️✔️:   res.reserve(matches.size());
    // RDKit✔️✔️:   seen.reserve(matches.size());
    // RDKit✔️✔️:   for (const auto &match : matches) {
    // RDKit✔️✔️:     std::string val(nAtoms, '0');
    // RDKit✔️✔️:     for (const auto &ci : match) {
    // RDKit✔️✔️:       val[ci.second] = '1';
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     const bool inserted = seen.insert(std::move(val)).second;
    // RDKit✔️✔️:     if (inserted) {
    // RDKit✔️✔️:       res.push_back(match);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res.shrink_to_fit();
    // RDKit✔️✔️:   matches = std::move(res);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Local complexity review: both versions allocate one atom-count-sized
    // signature per examined match and use expected O(1) hash insertion, for
    // O(matches * atom_count) time and space bounded by unique signatures.
    // Vec<bool> packs the same binary information as the source string. Moving
    // accepted Rust match values avoids the source copy and preserves order.
    let mut seen = HashSet::with_capacity(matches.len());
    let mut unique = Vec::with_capacity(matches.len());
    for matched in matches.drain(..) {
        let mut signature = vec![false; atom_count];
        for &atom_index in &matched.atom_mapping {
            signature[atom_index] = true;
        }
        if seen.insert(signature) {
            unique.push(matched);
        }
    }
    unique.shrink_to_fit();
    *matches = unique;
}

fn query_contains_atomic_number(
    query: &crate::QueryNode<AtomQueryPredicate>,
    atomic_number: u8,
) -> bool {
    match query {
        crate::QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(value)) => {
            *value == atomic_number
        }
        crate::QueryNode::And(children)
        | crate::QueryNode::Or(children)
        | crate::QueryNode::Xor(children) => children
            .iter()
            .any(|child| query_contains_atomic_number(child, atomic_number)),
        crate::QueryNode::Not(child) => query_contains_atomic_number(child, atomic_number),
        crate::QueryNode::Predicate(_) => false,
    }
}

pub(crate) fn is_atom_terminal_r_group_or_query_hydrogen(
    molecule: &SearchTarget<'_>,
    atom_index: usize,
) -> bool {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: isAtomTerminalRGroupOrQueryHydrogen
    // RDKit✔️✔️: bool isAtomTerminalRGroupOrQueryHydrogen(const Atom *atom) {
    // RDKit✔️✔️:   return (atom->getDegree() == 1 && isAtomDummy(atom)) ||
    // RDKit✔️✔️:          (atom->hasQuery() &&
    // RDKit✔️✔️:           describeQuery(atom).find("AtomAtomicNum 1 = val") !=
    // RDKit✔️✔️:               std::string::npos);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.h :: isAtomDummy
    // RDKit✔️✔️: inline bool isAtomDummy(const Atom *a) {
    // RDKit✔️✔️:   return (!a->hasQuery() && a->getAtomicNum() == 0) ||
    // RDKit✔️✔️:          (a->hasQuery() && !a->getQuery()->getNegation() &&
    // RDKit✔️✔️:           a->getQuery()->getDescription() == "AtomNull");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Local complexity review: degree is an indexed adjacency-slice length;
    // dummy classification is O(1), and the typed query traversal is O(n)
    // time/O(h) stack, matching describeQuery's traversal without allocating
    // its intermediate string. No molecule state or query node is cloned.
    let atom = &molecule.atoms()[atom_index];
    let is_dummy = atom.atomic_number() == 0;
    (molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom_index)
        .len()
        == 1
        && is_dummy)
        || false
}

fn core_substitution_score(
    molecule: &SearchTarget<'_>,
    query: &SearchTarget<'_>,
    matched: &SubstructMatchResult,
) -> f64 {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: detail::ScoreMatchesByDegreeOfCoreSubstitution
    // RDKit✔️✔️: class ScoreMatchesByDegreeOfCoreSubstitution {
    // RDKit✔️✔️:  public:
    // RDKit✔️✔️:   typedef std::pair<unsigned int, double> IdxScorePair;
    // RDKit✔️✔️:   ScoreMatchesByDegreeOfCoreSubstitution(
    // RDKit✔️✔️:       const RDKit::ROMol &mol, const RDKit::ROMol &query,
    // RDKit✔️✔️:       const std::vector<RDKit::MatchVectType> &matches)
    // RDKit✔️✔️:       : d_mol(mol),
    // RDKit✔️✔️:         d_query(query),
    // RDKit✔️✔️:         d_matches(matches),
    // RDKit✔️✔️:         d_sumIndices(0.0),
    // RDKit✔️✔️:         d_minIdx(-1),
    // RDKit✔️✔️:         d_isSorted(false) {
    // RDKit✔️✔️:     PRECONDITION(!matches.empty(), "matches must not be empty");
    // RDKit✔️✔️:     auto na = d_mol.getNumAtoms();
    // RDKit✔️✔️:     d_sumIndices = static_cast<double>(na * (na + 1) / 2);
    // RDKit✔️✔️:     unsigned int i = 0;
    // RDKit✔️✔️:     d_matchIdxVsScore.reserve(d_matches.size());
    // RDKit✔️✔️:     for (const auto &match : d_matches) {
    // RDKit✔️✔️:       d_matchIdxVsScore.emplace_back(i++, computeScore(match));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const RDKit::MatchVectType &getMostSubstitutedCoreMatch() {
    // RDKit✔️✔️:     if (d_minIdx == -1) {
    // RDKit✔️✔️:       d_minIdx = std::min_element(d_matchIdxVsScore.begin(),
    // RDKit✔️✔️:                                   d_matchIdxVsScore.end(), compare)
    // RDKit✔️✔️:                      ->first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return d_matches.at(d_minIdx);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution() {
    // RDKit✔️✔️:     if (!d_isSorted) {
    // RDKit✔️✔️:       std::sort(d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), compare);
    // RDKit✔️✔️:       d_isSorted = true;
    // RDKit✔️✔️:       d_minIdx = d_matchIdxVsScore.front().first;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     std::vector<MatchVectType> res(d_matches.size());
    // RDKit✔️✔️:     std::transform(
    // RDKit✔️✔️:         d_matchIdxVsScore.begin(), d_matchIdxVsScore.end(), res.begin(),
    // RDKit✔️✔️:         [this](const IdxScorePair &pair) { return d_matches.at(pair.first); });
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:  private:
    // RDKit✔️✔️:   static bool compare(const IdxScorePair &aPair, const IdxScorePair &bPair) {
    // RDKit✔️✔️:     return (aPair.second < bPair.second);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   bool doesRGroupMatchHydrogen(const std::pair<int, int> &pair) const {
    // RDKit✔️✔️:     const auto queryAtom = d_query.getAtomWithIdx(pair.first);
    // RDKit✔️✔️:     const auto molAtom = d_mol.getAtomWithIdx(pair.second);
    // RDKit✔️✔️:     return (molAtom->getAtomicNum() == 1 &&
    // RDKit✔️✔️:             isAtomTerminalRGroupOrQueryHydrogen(queryAtom));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   double computeScore(const RDKit::MatchVectType &match) const {
    // RDKit✔️✔️:     double penalty = 0.0;
    // RDKit✔️✔️:     double i = 0.0;
    // RDKit✔️✔️:     for (const auto &pair : match) {
    // RDKit✔️✔️:       i += static_cast<double>(pair.second);
    // RDKit✔️✔️:       if (doesRGroupMatchHydrogen(pair)) {
    // RDKit✔️✔️:         penalty += 1.0;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     penalty += i / d_sumIndices;
    // RDKit✔️✔️:     return penalty;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   const RDKit::ROMol &d_mol;
    // RDKit✔️✔️:   const RDKit::ROMol &d_query;
    // RDKit✔️✔️:   const std::vector<RDKit::MatchVectType> &d_matches;
    // RDKit✔️✔️:   std::vector<IdxScorePair> d_matchIdxVsScore;
    // RDKit✔️✔️:   double d_sumIndices;
    // RDKit✔️✔️:   int d_minIdx;
    // RDKit✔️✔️:   bool d_isSorted;
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION
    //
    // The Rust wrappers compute and retain the same per-match scores without
    // materializing a stateful scorer object.
    let atom_count = molecule.num_atoms();
    let sum_indices = (atom_count * (atom_count + 1) / 2) as f64;
    let mut penalty = 0.0;
    let mut index_sum = 0.0;
    for (query_index, &molecule_index) in matched.atom_mapping.iter().enumerate() {
        index_sum += molecule_index as f64;
        if molecule.atoms()[molecule_index].atomic_number() == 1
            && is_atom_terminal_r_group_or_query_hydrogen(query, query_index)
        {
            penalty += 1.0;
        }
    }
    penalty + index_sum / sum_indices
}

fn get_most_substituted_core_match<'a>(
    molecule: &SearchTarget<'_>,
    query: &SearchTarget<'_>,
    matches: &'a [SubstructMatchResult],
) -> &'a SubstructMatchResult {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: getMostSubstitutedCoreMatch
    // RDKit✔️✔️: const MatchVectType &getMostSubstitutedCoreMatch(
    // RDKit✔️✔️:     const ROMol &mol, const ROMol &core,
    // RDKit✔️✔️:     const std::vector<MatchVectType> &matches) {
    // RDKit✔️✔️:   detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
    // RDKit✔️✔️:                                                              matches);
    // RDKit✔️✔️:   return matchScorer.getMostSubstitutedCoreMatch();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // The canonical scorer above reproduces the complete source helper class.
    // Local complexity review: one linear score pass and min selection gives
    // O(matches * query_atoms) time and O(1) auxiliary space, equivalent to
    // constructing and scanning RDKit's score vector, with fewer allocations.
    assert!(!matches.is_empty(), "matches must not be empty");
    matches
        .iter()
        .min_by(|left, right| {
            core_substitution_score(molecule, query, left)
                .total_cmp(&core_substitution_score(molecule, query, right))
        })
        .expect("non-empty matches")
}

#[cfg(test)]
mod q86_atom_dispatch_tests {
    use super::*;
    use cosmolkit_model::{AtomId, AtomSpec, CoordinateBlock, QueryStateRef, TopologyBlock};
    use cosmolkit_types::Element;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn atom(element: Element, charge: i8, isotope: Option<u16>) -> Atom {
        let mut spec = AtomSpec::new(element).with_formal_charge(charge);
        if let Some(isotope) = isotope {
            spec = spec.with_isotope(isotope);
        }
        Atom::from_spec(AtomId::new(0), spec)
    }

    fn topology(atom: Atom) -> TopologyBlock {
        TopologyBlock::try_from_parts(vec![atom], Vec::new(), Vec::new(), Vec::new()).unwrap()
    }

    fn graph(atom: QueryAtom) -> QueryGraph {
        QueryGraph::from_parts(
            vec![atom],
            Vec::new(),
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn compat(
        query: &QueryGraph,
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        rows: &[QueryAtom],
        params: &SubstructMatchParams,
    ) -> bool {
        let state = QueryStateRef::try_for_topology(rows, &[], topology).unwrap();
        let target = SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None)
            .try_with_query_state(state)
            .unwrap();
        let context = build_query_match_context(&target);
        atom_compat(
            &query.atoms()[0],
            query,
            &topology.atoms[0],
            &target,
            params,
            None,
            &context,
        )
        .unwrap()
    }

    #[test]
    fn q86_atom_dispatch_requires_option_and_both_explicit_origins() {
        let carrier = atom(Element::C, 0, None);
        let topology = topology(carrier.clone());
        let coordinates = CoordinateBlock::default();
        let oxygen = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));
        let query = graph(QueryAtom::from_parts(carrier.clone(), oxygen.clone()));
        let explicit_rows = vec![QueryAtom::from_parts(carrier.clone(), oxygen.clone())];
        let carrier_rows = vec![QueryAtom::from_carrier_parts(
            carrier.clone(),
            oxygen.clone(),
        )];

        let mut params = SubstructMatchParams::default();
        assert!(!compat(
            &query,
            &topology,
            &coordinates,
            &explicit_rows,
            &params,
        ));
        params.use_query_query_matches = true;
        assert!(compat(
            &query,
            &topology,
            &coordinates,
            &explicit_rows,
            &params,
        ));
        assert!(!compat(
            &query,
            &topology,
            &coordinates,
            &carrier_rows,
            &params,
        ));

        let carrier_query = graph(QueryAtom::from_carrier_parts(
            carrier.clone(),
            oxygen.clone(),
        ));
        assert!(compat(
            &carrier_query,
            &topology,
            &coordinates,
            &explicit_rows,
            &params,
        ));
    }

    #[test]
    fn q86_atom_dispatch_uses_current_target_carrier_for_fallbacks() {
        let current_oxygen = atom(Element::O, 1, Some(14));
        let oxygen_topology = topology(current_oxygen.clone());
        let coordinates = CoordinateBlock::default();
        let old_carrier = atom(Element::C, 0, Some(13));
        let target_rows = vec![QueryAtom::from_parts(
            old_carrier,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        )];
        let oxygen_query = graph(QueryAtom::from_parts(
            atom(Element::C, 0, None),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        ));

        let mut params = SubstructMatchParams::default();
        assert!(compat(
            &oxygen_query,
            &oxygen_topology,
            &coordinates,
            &target_rows,
            &params,
        ));
        params.use_query_query_matches = true;
        assert!(!compat(
            &oxygen_query,
            &oxygen_topology,
            &coordinates,
            &target_rows,
            &params,
        ));

        let constrained_carrier = atom(Element::C, 1, Some(13));
        let carrier_query = graph(QueryAtom::from_carrier_parts(
            constrained_carrier.clone(),
            QueryNode::predicate(AtomQueryPredicate::Any),
        ));
        let target_rows = vec![QueryAtom::from_parts(
            constrained_carrier.clone(),
            QueryNode::predicate(AtomQueryPredicate::Any),
        )];
        let matching_topology = topology(constrained_carrier);
        assert!(compat(
            &carrier_query,
            &matching_topology,
            &coordinates,
            &target_rows,
            &params,
        ));
        let mismatching_topology = topology(atom(Element::C, 0, Some(14)));
        let mismatching_rows = vec![QueryAtom::from_parts(
            atom(Element::C, 0, None),
            QueryNode::predicate(AtomQueryPredicate::Any),
        )];
        assert!(!compat(
            &carrier_query,
            &mismatching_topology,
            &coordinates,
            &mismatching_rows,
            &params,
        ));
    }

    #[test]
    fn q86_atom_dispatch_preserves_override_property_and_post_callback_order() {
        let mut query_atom = QueryAtom::from_parts(
            atom(Element::C, 0, None),
            QueryNode::predicate(AtomQueryPredicate::Any),
        );
        query_atom.set_prop("gate", "query").unwrap();
        let query = graph(query_atom);
        let mut target_atom = atom(Element::C, 0, None);
        target_atom.set_prop("gate", "target").unwrap();
        let topology = topology(target_atom.clone());
        let coordinates = CoordinateBlock::default();
        let rows = vec![QueryAtom::from_parts(
            target_atom,
            QueryNode::predicate(AtomQueryPredicate::Any),
        )];

        let calls = Arc::new(AtomicUsize::new(0));
        let callback_calls = Arc::clone(&calls);
        let mut params = SubstructMatchParams::default();
        params.use_query_query_matches = true;
        params.atom_properties = vec!["gate".to_owned()];
        params.extra_atom_check = Some(Arc::new(move |_, _, _, _| {
            callback_calls.fetch_add(1, Ordering::SeqCst);
            true
        }));
        assert!(!compat(&query, &topology, &coordinates, &rows, &params));
        assert_eq!(calls.load(Ordering::SeqCst), 0);

        params.extra_atom_check_overrides_default_check = true;
        assert!(compat(&query, &topology, &coordinates, &rows, &params));
        assert_eq!(calls.load(Ordering::SeqCst), 1);

        params.extra_atom_check_overrides_default_check = false;
        params.atom_properties.clear();
        let rejecting_calls = Arc::new(AtomicUsize::new(0));
        let callback_calls = Arc::clone(&rejecting_calls);
        params.extra_atom_check = Some(Arc::new(move |_, _, _, _| {
            callback_calls.fetch_add(1, Ordering::SeqCst);
            false
        }));
        assert!(!compat(&query, &topology, &coordinates, &rows, &params));
        assert_eq!(rejecting_calls.load(Ordering::SeqCst), 1);
    }
}

#[cfg(test)]
mod q86_bond_dispatch_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondSpec, CoordinateBlock, QueryStateRef, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn atom(index: usize) -> Atom {
        Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C))
    }

    fn bond(order: BondOrder, conjugated: bool) -> Bond {
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), order).with_conjugated(conjugated),
        )
    }

    fn topology(bond: Bond) -> TopologyBlock {
        TopologyBlock::try_from_parts(vec![atom(0), atom(1)], vec![bond], Vec::new(), Vec::new())
            .unwrap()
    }

    fn graph(bond: QueryBond) -> QueryGraph {
        let atoms = vec![
            QueryAtom::from_carrier_parts(
                atom(0),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            ),
            QueryAtom::from_carrier_parts(
                atom(1),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            ),
        ];
        QueryGraph::from_parts(
            atoms,
            vec![bond],
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    fn compat(
        query: &QueryGraph,
        topology: &TopologyBlock,
        coordinates: &CoordinateBlock,
        target_bond: QueryBond,
        params: &SubstructMatchParams,
    ) -> bool {
        let atom_rows: Vec<_> = topology
            .atoms
            .iter()
            .cloned()
            .map(|atom| {
                let atomic_number = atom.atomic_number();
                QueryAtom::from_carrier_parts(
                    atom,
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
                )
            })
            .collect();
        let bond_rows = [target_bond];
        let state = QueryStateRef::try_for_topology(&atom_rows, &bond_rows, topology).unwrap();
        let target = SearchTarget::new(topology, coordinates, &topology.stereo_groups, None, None)
            .try_with_query_state(state)
            .unwrap();
        let context = build_query_match_context(&target);
        bond_compat(
            &query.bonds()[0],
            query,
            &topology.bonds[0],
            &target,
            params,
            &context,
        )
    }

    #[test]
    fn q86_bond_dispatch_requires_option_and_both_explicit_origins() {
        let current = bond(BondOrder::Single, false);
        let topology = topology(current.clone());
        let coordinates = CoordinateBlock::default();
        let query = graph(QueryBond::from_parts(
            current.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        ));
        let explicit_double = QueryBond::from_parts(
            current.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        );
        let carrier_double = QueryBond::from_carrier_parts(
            current.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        );

        let mut params = SubstructMatchParams::default();
        assert!(compat(
            &query,
            &topology,
            &coordinates,
            explicit_double.clone(),
            &params,
        ));
        params.use_query_query_matches = true;
        assert!(!compat(
            &query,
            &topology,
            &coordinates,
            explicit_double,
            &params,
        ));
        assert!(compat(
            &query,
            &topology,
            &coordinates,
            carrier_double,
            &params,
        ));

        let carrier_query = graph(QueryBond::from_carrier_parts(
            current.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
        ));
        let explicit_target =
            QueryBond::from_parts(current, QueryNode::predicate(BondQueryPredicate::Any));
        assert!(!compat(
            &carrier_query,
            &topology,
            &coordinates,
            explicit_target,
            &params,
        ));
    }

    #[test]
    fn q86_bond_dispatch_uses_bond_null_and_source_second_and_relation() {
        let current = bond(BondOrder::Single, false);
        let topology = topology(current.clone());
        let coordinates = CoordinateBlock::default();
        let query = graph(QueryBond::from_parts(
            current.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        ));
        let mut params = SubstructMatchParams::default();
        params.use_query_query_matches = true;

        for (target_predicate, expected) in [
            (QueryNode::predicate(BondQueryPredicate::Any), true),
            (
                QueryNode::and(vec![QueryNode::predicate(BondQueryPredicate::Order(
                    BondOrder::Double,
                ))]),
                true,
            ),
            (
                QueryNode::and(vec![QueryNode::predicate(BondQueryPredicate::Order(
                    BondOrder::Single,
                ))]),
                false,
            ),
        ] {
            let target_bond = QueryBond::from_parts(current.clone(), target_predicate);
            assert_eq!(
                compat(&query, &topology, &coordinates, target_bond, &params),
                expected
            );
        }
    }

    #[test]
    fn q86_bond_dispatch_preserves_aromatic_and_current_property_precedence() {
        let aromatic = bond(BondOrder::Aromatic, false);
        let query = graph(QueryBond::from_carrier_parts(
            aromatic.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
        ));
        let conjugated_single = bond(BondOrder::Single, true);
        let conjugated_topology = topology(conjugated_single.clone());
        let coordinates = CoordinateBlock::default();
        let carrier_target = QueryBond::from_carrier_parts(
            conjugated_single.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        let mut params = SubstructMatchParams::default();
        params.aromatic_matches_conjugated = true;
        assert!(compat(
            &query,
            &conjugated_topology,
            &coordinates,
            carrier_target,
            &params,
        ));
        let explicit_target = QueryBond::from_parts(
            conjugated_single,
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert!(!compat(
            &query,
            &conjugated_topology,
            &coordinates,
            explicit_target,
            &params,
        ));

        let mut current_single = bond(BondOrder::Single, false);
        current_single.set_prop("gate", "current").unwrap();
        let current_topology = topology(current_single);
        let mut query_double = QueryBond::from_parts(
            bond(BondOrder::Double, false),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        );
        query_double.bond_mut().set_prop("gate", "query").unwrap();
        let query = graph(query_double);
        let mut old_double = bond(BondOrder::Double, false);
        old_double.set_prop("gate", "query").unwrap();
        let target_row = QueryBond::from_parts(
            old_double,
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
        );
        params.use_query_query_matches = true;
        params.bond_properties.clear();
        assert!(compat(
            &query,
            &current_topology,
            &coordinates,
            target_row.clone(),
            &params,
        ));
        params.bond_properties = vec!["gate".to_owned()];
        assert!(!compat(
            &query,
            &current_topology,
            &coordinates,
            target_row,
            &params,
        ));
    }
}

fn sort_matches_by_degree_of_core_substitution(
    molecule: &SearchTarget<'_>,
    query: &SearchTarget<'_>,
    matches: &[SubstructMatchResult],
) -> Vec<SubstructMatchResult> {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/Substruct/SubstructUtils.cpp :: sortMatchesByDegreeOfCoreSubstitution
    // RDKit✔️✔️: std::vector<MatchVectType> sortMatchesByDegreeOfCoreSubstitution(
    // RDKit✔️✔️:     const ROMol &mol, const ROMol &core,
    // RDKit✔️✔️:     const std::vector<MatchVectType> &matches) {
    // RDKit✔️✔️:   detail::ScoreMatchesByDegreeOfCoreSubstitution matchScorer(mol, core,
    // RDKit✔️✔️:                                                              matches);
    // RDKit✔️✔️:   return matchScorer.sortMatchesByDegreeOfCoreSubstitution();
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION
    //
    // Local complexity review: scores are computed once and the indexed rows
    // are sorted in O(matches log matches), matching the source helper. The
    // returned mappings are cloned once, as in RDKit's result transform.
    assert!(!matches.is_empty(), "matches must not be empty");
    let mut scored = matches
        .iter()
        .enumerate()
        .map(|(index, matched)| (index, core_substitution_score(molecule, query, matched)))
        .collect::<Vec<_>>();
    scored.sort_by(|left, right| left.1.total_cmp(&right.1));
    scored
        .into_iter()
        .map(|(index, _)| matches[index].clone())
        .collect()
}

fn substruct_match_impl(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
) -> SubstructMatchResultList {
    // RDKit✔️❌: std::vector<MatchVectType> SubstructMatch(
    // RDKit✔️❌:     const ROMol &mol, const ROMol &query,
    // RDKit✔️❌:     const SubstructMatchParameters &params) {
    // RDKit✔️❌:   std::vector<MatchVectType> matches;
    // RDKit✔️❌:   const auto &mNumAtoms = mol.getNumAtoms();
    // RDKit✔️❌:   const auto &qNumAtoms = query.getNumAtoms();
    // RDKit✔️❌:   if (!mNumAtoms || !qNumAtoms || qNumAtoms > mNumAtoms) {
    // RDKit✔️❌:     return matches;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   detail::RecursiveLocker locker(query, params.recursionPossible);
    // RDKit✔️❌:
    // RDKit✔️❌:   if (params.recursionPossible) {
    // RDKit✔️❌:     detail::SUBQUERY_MAP subqueryMap;
    // RDKit✔️❌:     ROMol::ConstAtomIterator atIt;
    // RDKit✔️❌:     for (const auto atom : query.atoms()) {
    // RDKit✔️❌:       if (atom->hasQuery()) {
    // RDKit✔️❌:         detail::MatchSubqueries(mol, atom->getQuery(), params, subqueryMap,
    // RDKit✔️❌:                                 locker.locked);
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   detail::AtomLabelFunctor atomLabeler(query, mol, params);
    // RDKit✔️❌:   detail::BondLabelFunctor bondLabeler(query, mol, params);
    // RDKit✔️❌:   MolMatchFinalCheckFunctor matchChecker(query, mol, params);
    // RDKit✔️❌:
    // RDKit✔️❌:   std::vector<detail::ssPairType> pms;
    // RDKit✔️❌:   bool found =
    // RDKit✔️❌:       boost::vf2_all(query.getTopology(), mol.getTopology(), atomLabeler,
    // RDKit✔️❌:                      bondLabeler, matchChecker, pms, params.maxMatches);
    // RDKit✔️❌:   if (found) {
    // RDKit✔️❌:     const unsigned int nQueryAtoms = query.getNumAtoms();
    // RDKit✔️❌:     matches.reserve(pms.size());
    // RDKit✔️❌:     MatchVectType matchVect(nQueryAtoms);
    // RDKit✔️❌:     for (const auto &pairs : pms) {
    // RDKit✔️❌:       for (const auto &pair : pairs) {
    // RDKit✔️❌:         matchVect[pair.first] = pair;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       matches.push_back(matchVect);
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return matches;
    // RDKit✔️❌: }
    // Complexity review: the fail-closed preflight adds one O(A+B+Q) query
    // scan. The VF2 core now allocates its two query-sized mapping buffers
    // once per invocation, reuses them at every goal, allocates no mapping
    // output for rejected goals, and appends one paired vector per accepted
    // all-match result. Ordinary/context calls borrow existing query/target
    // adjacency; explicit compiled calls borrow their retained graph. These
    // remove the prior per-goal mapping vectors and per-call graph rebuilds.
    // The second marker remains ❌ for the distinct accepted-result projection
    // cost: each row currently builds a temporary Vec<Option<usize>> and then
    // collects a separate output atom map (along with the result bond map),
    // while RDKit reuses its MatchVectType before copying each output row.
    preflight_query_molecule(query)?;
    if mol.num_atoms() == 0 || query.num_atoms() == 0 || query.num_atoms() > mol.num_atoms() {
        return Ok(Vec::new());
    }
    let mut recursive_locker = RecursiveLocker::new(query, params.recursion_possible);
    if params.recursion_possible {
        populate_recursive_query_match_cache(
            mol,
            query,
            params,
            &mut recursive_locker.cache,
            None,
        )?;
    }
    substruct_match_impl_with_recursive_cache::<FullMatchResultProjection>(
        mol,
        query,
        params,
        Some(&recursive_locker.cache),
    )
}

/// Check if a molecule contains a substructure match for the given query.
///
/// This is the public API for `has_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn has_substruct_match<Q: QueryInput + ?Sized>(mol: &SearchTarget<'_>, query: &Q) -> bool {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    let Ok(query) = query.query_graph() else {
        return false;
    };
    substruct_match_impl(mol, &query, &params)
        .map(|matches| !matches.is_empty())
        .unwrap_or(false)
}

/// Get the first substructure match, if any.
///
/// This is the public API for `get_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_match<Q: QueryInput + ?Sized>(
    mol: &SearchTarget<'_>,
    query: &Q,
) -> Option<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    let query = query.query_graph().ok()?;
    substruct_match_impl(mol, &query, &params)
        .ok()
        .and_then(|matches| matches.into_iter().next())
}

/// Get all substructure matches with default parameters.
///
/// This is the public API for `get_substruct_matches`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches<Q: QueryInput + ?Sized>(
    mol: &SearchTarget<'_>,
    query: &Q,
) -> Vec<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    let Ok(query) = query.query_graph() else {
        return Vec::new();
    };
    substruct_match_impl(mol, &query, &params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters.
///
/// This is the public API for `get_substruct_matches_with_params`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches_with_params<Q: QueryInput + ?Sized>(
    mol: &SearchTarget<'_>,
    query: &Q,
    params: &SubstructMatchParams,
) -> Vec<SubstructMatchResult> {
    let Ok(query) = query.query_graph() else {
        return Vec::new();
    };
    substruct_match_impl(mol, &query, params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters and structured
/// unsupported-feature errors for source-porting callers.
pub fn try_get_substruct_matches_with_params<Q: QueryInput + ?Sized>(
    mol: &SearchTarget<'_>,
    query: &Q,
    params: &SubstructMatchParams,
) -> SubstructMatchResultList {
    let query = query.query_graph()?;
    substruct_match_impl(mol, &query, params)
}

#[doc(hidden)]
pub fn try_get_substruct_matches_with_params_and_context(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    query_context: &QueryMatchContext,
) -> SubstructMatchResultList {
    // This narrow entry retains the canonical preflight, recursive-query
    // preparation, VF2 implementation, final checks, and result ordering. It
    // only lets callers that run several immutable queries against one target
    // reuse the target-derived match context, as RDKit reuses ROMol state.
    preflight_query_molecule(query)?;
    if mol.num_atoms() == 0 || query.num_atoms() == 0 || query.num_atoms() > mol.num_atoms() {
        return Ok(Vec::new());
    }
    let mut recursive_locker = RecursiveLocker::new(query, params.recursion_possible);
    if params.recursion_possible {
        populate_recursive_query_match_cache(
            mol,
            query,
            params,
            &mut recursive_locker.cache,
            Some(query_context),
        )?;
    }
    substruct_match_impl_with_recursive_cache_and_context::<FullMatchResultProjection>(
        mol,
        query,
        params,
        Some(&recursive_locker.cache),
        query_context,
        None,
        None,
    )
}

pub(crate) fn compile_query_order_from_graph(query: &CompiledQueryGraph) -> Vec<usize> {
    sort_nodes_by_frequency(Vf2GraphRef::compiled(query))
}

pub(crate) fn compile_query_graph(query: &QueryGraph) -> CompiledQueryGraph {
    build_vf2_graph(query)
}

pub(crate) fn get_substruct_matches_with_compiled_query(
    mol: &SearchTarget<'_>,
    query: &QueryGraph,
    params: &SubstructMatchParams,
    compiled_graph: &CompiledQueryGraph,
) -> SubstructMatchResultList {
    substruct_matches_with_compiled_query::<FullMatchResultProjection>(
        mol,
        query,
        params,
        compiled_graph,
    )
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod q66_node_info_comparator_tests {
    use super::{NodeInfo, node_info_cmp1, node_info_cmp2};
    use std::cmp::Ordering;

    const fn info(id: u32, in_deg: u32, out_deg: u32) -> NodeInfo {
        NodeInfo {
            id,
            in_deg,
            out_deg,
        }
    }

    #[test]
    fn q66_node_info_cmp1_orders_out_then_in_and_preserves_ties() {
        assert_eq!(
            node_info_cmp1(&info(0, u32::MAX, 1), &info(1, 0, 2)),
            Ordering::Less
        );
        assert_eq!(
            node_info_cmp1(&info(0, 0, 2), &info(1, u32::MAX, 1)),
            Ordering::Greater
        );
        assert_eq!(
            node_info_cmp1(&info(0, 1, u32::MAX), &info(1, 2, u32::MAX)),
            Ordering::Less
        );
        assert_eq!(
            node_info_cmp1(&info(0, 2, u32::MAX), &info(1, 1, u32::MAX)),
            Ordering::Greater
        );
        assert_eq!(
            node_info_cmp1(
                &info(u32::MAX, u32::MAX, u32::MAX),
                &info(0, u32::MAX, u32::MAX)
            ),
            Ordering::Equal
        );
    }

    #[test]
    fn q66_node_info_cmp2_orders_zero_frequency_then_out_then_in() {
        assert_eq!(
            node_info_cmp2(&info(0, 0, 0), &info(1, 1, u32::MAX)),
            Ordering::Greater
        );
        assert_eq!(
            node_info_cmp2(&info(0, 1, u32::MAX), &info(1, 0, 0)),
            Ordering::Less
        );
        assert_eq!(
            node_info_cmp2(&info(0, 0, 1), &info(1, 0, 2)),
            Ordering::Less
        );
        assert_eq!(
            node_info_cmp2(&info(0, 0, u32::MAX), &info(1, 0, 0)),
            Ordering::Greater
        );
        assert_eq!(
            node_info_cmp2(&info(0, 1, 1), &info(1, 2, 1)),
            Ordering::Less
        );
        assert_eq!(
            node_info_cmp2(&info(0, 2, 1), &info(1, 1, 1)),
            Ordering::Greater
        );
        assert_eq!(
            node_info_cmp2(
                &info(u32::MAX, u32::MAX, u32::MAX),
                &info(0, u32::MAX, u32::MAX)
            ),
            Ordering::Equal
        );
    }
}

#[cfg(test)]
mod q33_plain_atom_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondSpec, CoordinateBlock, QueryAtomIdentity, QueryBond,
        TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn atom(
        index: usize,
        element: Element,
        formal_charge: i8,
        isotope: Option<u16>,
        radical_electrons: u8,
        aromatic: bool,
        explicit_hydrogens: u8,
    ) -> Atom {
        let mut spec = AtomSpec::new(element)
            .with_formal_charge(formal_charge)
            .with_radical_electrons(radical_electrons)
            .with_aromatic(aromatic)
            .with_explicit_hydrogens(explicit_hydrogens);
        if let Some(isotope) = isotope {
            spec = spec.with_isotope(isotope);
        }
        Atom::from_spec(AtomId::new(index), spec)
    }

    fn carrier_query(atom: Atom) -> QueryAtom {
        QueryAtom::from_carrier_parts(
            atom.clone(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
        )
    }

    fn topology(atoms: Vec<Atom>, bonds: Vec<Bond>) -> TopologyBlock {
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed query-matching topology is valid")
    }

    fn plain_match(
        query_atom: &QueryAtom,
        target_atom: Atom,
        target_atomic_number: Option<u8>,
    ) -> bool {
        let topology = topology(vec![target_atom], Vec::new());
        let coordinates = CoordinateBlock::default();
        if let Some(atomic_number) = target_atomic_number {
            let overrides = [Some(atomic_number)];
            let target =
                SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None)
                    .with_atomic_number_overrides(&overrides);
            atom_matches(query_atom, &topology.atoms[0], &target)
        } else {
            let target =
                SearchTarget::new(&topology, &coordinates, &topology.stereo_groups, None, None);
            atom_matches(query_atom, &topology.atoms[0], &target)
        }
    }

    fn query_graph(atoms: Vec<QueryAtom>, bonds: Vec<QueryBond>) -> QueryGraph {
        QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed detached query graph is valid")
    }

    #[test]
    fn q33_plain_atom_uses_current_query_and_effective_target_identity() {
        let query = carrier_query(atom(0, Element::C, 0, None, 0, false, 0));
        assert!(plain_match(
            &query,
            atom(0, Element::C, 0, None, 0, false, 0),
            None
        ));
        assert!(!plain_match(
            &query,
            atom(0, Element::N, 0, None, 0, false, 0),
            None
        ));
        assert!(!plain_match(
            &query,
            atom(0, Element::C, 0, None, 0, false, 0),
            Some(7)
        ));

        // The ordinary carrier's current identity controls Atom::Match even
        // when its preserved query-tree snapshot still says atomic number 6.
        let reidentified = query
            .clone()
            .with_identity(QueryAtomIdentity::AtomicNumber(7));
        assert!(reidentified.predicate_is_carrier_derived());
        assert_eq!(
            reidentified.predicate(),
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6))
        );
        assert!(plain_match(
            &reidentified,
            atom(0, Element::C, 0, None, 0, false, 0),
            Some(7)
        ));
    }

    #[test]
    fn q33_plain_atom_dummy_isotopes_use_nonzero_source_defaults() {
        let mut zero_query = carrier_query(atom(0, Element::DUMMY, 0, None, 0, false, 0));
        zero_query.set_isotope(Some(0));
        assert_eq!(zero_query.isotope(), None);
        let mut zero_target = atom(0, Element::DUMMY, 0, None, 0, false, 0);
        zero_target.set_isotope(Some(0));
        assert_eq!(zero_target.isotope(), None);

        let wildcard_query = carrier_query(atom(0, Element::DUMMY, 0, None, 0, false, 0));
        for target_isotope in [None, Some(13)] {
            assert!(plain_match(
                &wildcard_query,
                atom(0, Element::DUMMY, 0, target_isotope, 0, false, 0),
                None
            ));
        }

        let labeled_query = carrier_query(atom(0, Element::DUMMY, 0, Some(12), 0, false, 0));
        for (target_isotope, expected) in [(None, true), (Some(12), true), (Some(13), false)] {
            assert_eq!(
                plain_match(
                    &labeled_query,
                    atom(0, Element::DUMMY, 0, target_isotope, 0, false, 0),
                    None
                ),
                expected,
                "target isotope {target_isotope:?}"
            );
        }
    }

    #[test]
    fn q33_plain_atom_nondefault_charge_isotope_and_radical_constraints_are_asymmetric() {
        for (query_charge, target_charge, expected) in [(1, 1, true), (1, 0, false), (0, 1, true)] {
            let query = carrier_query(atom(0, Element::C, query_charge, None, 0, false, 0));
            assert_eq!(
                plain_match(
                    &query,
                    atom(0, Element::C, target_charge, None, 0, false, 0),
                    None
                ),
                expected,
                "charge query={query_charge}, target={target_charge}"
            );
        }

        for (query_isotope, target_isotope, expected) in [
            (Some(13), Some(13), true),
            (Some(13), Some(14), false),
            (None, Some(13), true),
        ] {
            let query = carrier_query(atom(0, Element::C, 0, query_isotope, 0, false, 0));
            assert_eq!(
                plain_match(
                    &query,
                    atom(0, Element::C, 0, target_isotope, 0, false, 0),
                    None
                ),
                expected,
                "isotope query={query_isotope:?}, target={target_isotope:?}"
            );
        }

        for (query_radicals, target_radicals, expected) in
            [(2, 2, true), (2, 1, false), (0, 1, true)]
        {
            let query = carrier_query(atom(0, Element::C, 0, None, query_radicals, false, 0));
            assert_eq!(
                plain_match(
                    &query,
                    atom(0, Element::C, 0, None, target_radicals, false, 0),
                    None
                ),
                expected,
                "radicals query={query_radicals}, target={target_radicals}"
            );
        }

        let query = carrier_query(atom(0, Element::C, 0, None, 0, true, 4));
        assert!(plain_match(
            &query,
            atom(0, Element::C, 0, None, 0, false, 0),
            None
        ));
    }

    #[test]
    fn q33_plain_atom_current_identity_reaches_chiral_and_dative_callers() {
        let query_atom = carrier_query(atom(0, Element::C, 0, None, 0, false, 0));
        let query = query_graph(vec![query_atom.clone()], Vec::new());
        let single_topology = topology(vec![atom(0, Element::C, 0, None, 0, false, 0)], Vec::new());
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &single_topology,
            &coordinates,
            &single_topology.stereo_groups,
            None,
            None,
        );
        assert!(chiral_atom_compat(
            &query_atom,
            &query,
            &single_topology.atoms[0],
            &target
        ));
        let overrides = [Some(7)];
        let target_with_override = SearchTarget::new(
            &single_topology,
            &coordinates,
            &single_topology.stereo_groups,
            None,
            None,
        )
        .with_atomic_number_overrides(&overrides);
        assert!(!chiral_atom_compat(
            &query_atom,
            &query,
            &single_topology.atoms[0],
            &target_with_override
        ));

        let query_atoms = vec![
            carrier_query(atom(0, Element::C, 0, None, 0, false, 0)),
            carrier_query(atom(1, Element::N, 0, None, 0, false, 0)),
        ];
        let query_dative = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Dative),
        );
        let query_bond = QueryBond::from_carrier_parts(
            query_dative,
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Dative)),
        );
        let query_graph = query_graph(query_atoms, vec![query_bond]);
        let target_bond = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Dative),
        );
        let target_topology = topology(
            vec![
                atom(0, Element::C, 0, None, 0, false, 0),
                atom(1, Element::N, 0, None, 0, false, 0),
            ],
            vec![target_bond],
        );
        let target_coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            &target_topology,
            &target_coordinates,
            &target_topology.stereo_groups,
            None,
            None,
        );
        let target_context = build_query_match_context(&target);
        assert!(bond_compat(
            query_graph.bond(0).expect("query dative bond"),
            &query_graph,
            &target_topology.bonds[0],
            &target,
            &SubstructMatchParams::default(),
            &target_context,
        ));

        let endpoint_overrides = [Some(7), None];
        let target_with_override = SearchTarget::new(
            &target_topology,
            &target_coordinates,
            &target_topology.stereo_groups,
            None,
            None,
        )
        .with_atomic_number_overrides(&endpoint_overrides);
        let target_context = build_query_match_context(&target_with_override);
        assert!(!bond_compat(
            query_graph.bond(0).expect("query dative bond"),
            &query_graph,
            &target_topology.bonds[0],
            &target_with_override,
            &SubstructMatchParams::default(),
            &target_context,
        ));
    }
}

#[cfg(test)]
mod uff_one_fix_result_projection_tests {
    use super::*;
    use cosmolkit_model::{
        AtomId, AtomSpec, BondId, BondSpec, CoordinateBlock, QueryAtom, QueryBond,
        RecursiveStructureQuery, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn atom(index: usize, element: Element) -> Atom {
        Atom::from_spec(AtomId::new(index), AtomSpec::new(element))
    }

    fn target_topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| atom(index, element))
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed Search projection target is valid")
    }

    fn query_graph(elements: &[Element], edges: &[(usize, usize)]) -> QueryGraph {
        let atoms = elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| {
                let carrier = atom(index, element);
                QueryAtom::from_carrier_parts(
                    carrier.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(carrier.atomic_number())),
                )
            })
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end))| {
                let carrier = Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                );
                QueryBond::from_carrier_parts(
                    carrier,
                    QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                )
            })
            .collect();
        QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed Search projection query is valid")
    }

    fn assert_projection_matches_full(
        target_topology: &TopologyBlock,
        query: &QueryGraph,
        params: &SubstructMatchParams,
    ) -> Vec<Vec<usize>> {
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(
            target_topology,
            &coordinates,
            &target_topology.stereo_groups,
            None,
            None,
        );
        let compiled_graph = compile_query_graph(query);
        let query_context = build_query_match_context(&target);

        BOND_MAPPING_MATERIALIZATION_COUNT.with(|count| count.set(0));
        // Exercise the borrowed-context entry wired through the hidden
        // cross-crate adapter used by the forcefields caller.
        let atom_rows = try_get_substruct_atom_matches_with_compiled_query_and_context(
            &target,
            query,
            params,
            &compiled_graph,
            &query_context,
        )
        .expect("atom-only projection keeps canonical matcher errors");
        BOND_MAPPING_MATERIALIZATION_COUNT
            .with(|count| assert_eq!(count.get(), 0, "atom-only projection built a bond map"));

        let full_rows =
            get_substruct_matches_with_compiled_query(&target, query, params, &compiled_graph)
                .expect("full projection keeps canonical matcher errors");
        BOND_MAPPING_MATERIALIZATION_COUNT.with(|count| {
            assert_eq!(
                count.get(),
                full_rows.len(),
                "each full hit materializes exactly one bond map"
            )
        });
        let full_atom_rows = full_rows
            .iter()
            .map(|matched| matched.atom_mapping.clone())
            .collect::<Vec<_>>();
        assert_eq!(atom_rows, full_atom_rows);
        atom_rows
    }

    #[test]
    fn uff_one_fix_atom_projection_preserves_empty_and_no_hit_results() {
        let target = target_topology(&[Element::C, Element::C], &[(0, 1)]);
        let empty_query = query_graph(&[], &[]);
        assert!(assert_projection_matches_full(
            &target,
            &empty_query,
            &SubstructMatchParams::default(),
        )
        .is_empty());

        let nitrogen_query = query_graph(&[Element::N], &[]);
        assert!(
            assert_projection_matches_full(
                &target,
                &nitrogen_query,
                &SubstructMatchParams::default(),
            )
            .is_empty()
        );
    }

    #[test]
    fn uff_one_fix_atom_projection_preserves_multiple_hits_and_max_matches() {
        let target = target_topology(&[Element::C, Element::C, Element::C], &[]);
        let query = query_graph(&[Element::C], &[]);
        for uniquify in [true, false] {
            let params = SubstructMatchParams {
                uniquify,
                ..SubstructMatchParams::default()
            };
            assert_eq!(
                assert_projection_matches_full(&target, &query, &params),
                [vec![0], vec![1], vec![2]],
                "uniquify={uniquify}"
            );
        }

        let params = SubstructMatchParams {
            max_matches: 1,
            ..SubstructMatchParams::default()
        };
        assert_eq!(
            assert_projection_matches_full(&target, &query, &params),
            [vec![0]],
        );
    }

    #[test]
    fn uff_one_fix_atom_projection_preserves_symmetric_hit_order_and_uniqueness() {
        let target = target_topology(&[Element::C, Element::C, Element::C], &[(0, 1), (1, 2)]);
        let query = query_graph(&[Element::C, Element::C], &[(0, 1)]);

        assert_eq!(
            assert_projection_matches_full(&target, &query, &SubstructMatchParams::default(),),
            [vec![0, 1], vec![1, 2]],
        );
        let params = SubstructMatchParams {
            uniquify: false,
            ..SubstructMatchParams::default()
        };
        assert_eq!(
            assert_projection_matches_full(&target, &query, &params),
            [vec![0, 1], vec![1, 0], vec![1, 2], vec![2, 1]],
        );
    }

    #[test]
    fn uff_one_fix_atom_projection_preserves_recursive_query_results() {
        let target = target_topology(&[Element::C, Element::O, Element::C], &[(0, 1), (1, 2)]);
        let mut inner_query = crate::parse_smarts("C-O", &crate::SmartsParseParams::default())
            .expect("fixed recursive inner SMARTS parses");
        inner_query.set_prop("_queryRootAtom", "0");
        let recursive = RecursiveStructureQuery::from_query_graph(inner_query, 101);
        let recursive_atom = QueryAtom::from_parts(
            atom(0, Element::C),
            QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(recursive)),
        );
        let query = QueryGraph::from_parts(
            vec![recursive_atom],
            Vec::new(),
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed recursive SMARTS graph is valid");

        assert_eq!(
            assert_projection_matches_full(&target, &query, &SubstructMatchParams::default(),),
            [vec![0], vec![2]],
        );
    }
}

#[cfg(test)]
mod search_shared_perf_allocator {
    use std::alloc::{GlobalAlloc, Layout, System};
    use std::cell::Cell;

    std::thread_local! {
        static TRACK_ALLOCATIONS: Cell<bool> = const { Cell::new(false) };
        static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
        static REALLOCATIONS: Cell<usize> = const { Cell::new(0) };
    }

    struct TestCountingAllocator;

    // SAFETY: every allocation operation is forwarded to `System` with the
    // exact layout and pointer supplied by the allocator caller. The counters
    // observe only thread-local Cell state and do not inspect allocated bytes.
    unsafe impl GlobalAlloc for TestCountingAllocator {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            record_allocation();
            // SAFETY: this forwards the original layout unchanged to System.
            unsafe { System.alloc(layout) }
        }

        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            record_allocation();
            // SAFETY: this forwards the original layout unchanged to System.
            unsafe { System.alloc_zeroed(layout) }
        }

        unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
            // SAFETY: this forwards the original pointer/layout pair to System.
            unsafe { System.dealloc(pointer, layout) }
        }

        unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
            record_reallocation();
            // SAFETY: this forwards the original pointer/layout and requested
            // size unchanged to System.
            unsafe { System.realloc(pointer, layout, new_size) }
        }
    }

    #[global_allocator]
    static TEST_ALLOCATOR: TestCountingAllocator = TestCountingAllocator;

    #[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
    pub(super) struct AllocationCounts {
        pub(super) allocations: usize,
        pub(super) reallocations: usize,
    }

    fn record_allocation() {
        if TRACK_ALLOCATIONS.try_with(Cell::get).unwrap_or(false) {
            let _ = ALLOCATIONS.try_with(|count| count.set(count.get().saturating_add(1)));
        }
    }

    fn record_reallocation() {
        if TRACK_ALLOCATIONS.try_with(Cell::get).unwrap_or(false) {
            let _ = REALLOCATIONS.try_with(|count| count.set(count.get().saturating_add(1)));
        }
    }

    struct TrackingGuard;

    impl Drop for TrackingGuard {
        fn drop(&mut self) {
            let _ = TRACK_ALLOCATIONS.try_with(|tracking| tracking.set(false));
        }
    }

    pub(super) fn measure_allocations(operation: impl FnOnce()) -> AllocationCounts {
        // Warm all thread-local keys while tracking is disabled, so TLS setup
        // cannot contaminate the measured region.
        TRACK_ALLOCATIONS.with(|tracking| tracking.set(false));
        ALLOCATIONS.with(|count| count.set(0));
        REALLOCATIONS.with(|count| count.set(0));

        let guard = TrackingGuard;
        TRACK_ALLOCATIONS.with(|tracking| tracking.set(true));
        operation();
        drop(guard);

        AllocationCounts {
            allocations: ALLOCATIONS.with(Cell::get),
            reallocations: REALLOCATIONS.with(Cell::get),
        }
    }
}

#[cfg(test)]
mod search_projection_p01_tests {
    use super::search_shared_perf_allocator::measure_allocations;
    use super::{
        NULL_NODE, NodeId, VF2_GRAPH_BUILD_ENTRIES, Vf2GraphRef, project_substruct_matches,
    };
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, QueryAtom, QueryBond, QueryGraph,
        TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    struct FixedCase {
        query_atoms: usize,
        query_edges: &'static [(usize, usize)],
        target_atoms: usize,
        target_edges: &'static [(usize, usize)],
        raw_rows: &'static [&'static [(usize, usize)]],
        expected_atoms: &'static [&'static [usize]],
        expected_bonds: &'static [&'static [usize]],
    }

    const CASES: [FixedCase; 4] = [
        FixedCase {
            query_atoms: 1,
            query_edges: &[],
            target_atoms: 3,
            target_edges: &[(0, 1), (1, 2)],
            raw_rows: &[&[(0, 0)], &[(0, 2)]],
            expected_atoms: &[&[0], &[2]],
            expected_bonds: &[&[], &[]],
        },
        FixedCase {
            query_atoms: 2,
            query_edges: &[(0, 1)],
            target_atoms: 3,
            target_edges: &[(0, 1), (1, 2)],
            raw_rows: &[&[(0, 0), (1, 1)], &[(0, 2), (1, 1)]],
            expected_atoms: &[&[0, 1], &[2, 1]],
            expected_bonds: &[&[0], &[1]],
        },
        FixedCase {
            query_atoms: 3,
            query_edges: &[(0, 1), (1, 2)],
            target_atoms: 4,
            target_edges: &[(0, 1), (1, 2), (2, 3)],
            raw_rows: &[&[(0, 0), (1, 1), (2, 2)], &[(0, 3), (1, 2), (2, 1)]],
            expected_atoms: &[&[0, 1, 2], &[3, 2, 1]],
            expected_bonds: &[&[0, 1], &[2, 1]],
        },
        FixedCase {
            query_atoms: 2,
            query_edges: &[],
            target_atoms: 3,
            target_edges: &[(0, 1), (1, 2)],
            raw_rows: &[&[(0, 0), (1, 2)], &[(0, 2), (1, 0)]],
            expected_atoms: &[&[0, 2], &[2, 0]],
            expected_bonds: &[&[], &[]],
        },
    ];

    // Frozen from the contract: A/D [0,2,3], B/C [0,3,5].
    const EXPECTED_ALLOCATIONS: [[usize; 3]; 4] = [[0, 2, 3], [0, 3, 5], [0, 3, 5], [0, 2, 3]];

    fn query(case: &FixedCase) -> QueryGraph {
        let atoms = (0..case.query_atoms)
            .map(|index| QueryAtom::new(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = case
            .query_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                QueryBond::new(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("P01 literal query graph is valid")
    }

    fn target(case: &FixedCase) -> TopologyBlock {
        let atoms = (0..case.target_atoms)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = case
            .target_edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("P01 literal target topology is valid")
    }

    #[test]
    fn search_projection_p01_literal_rows_allocations_and_immutability() {
        VF2_GRAPH_BUILD_ENTRIES.with(|entries| entries.set(0));
        let mut actual_calls = 0;

        for (case_index, case) in CASES.iter().enumerate() {
            let query = query(case);
            let target = target(case);

            for goal_count in 0..=2 {
                for reverse_pairs in [false, true] {
                    let raw_matches: Vec<Vec<(NodeId, NodeId)>> = case
                        .raw_rows
                        .iter()
                        .take(goal_count)
                        .map(|row| {
                            if reverse_pairs {
                                row.iter().copied().rev().collect()
                            } else {
                                row.to_vec()
                            }
                        })
                        .collect();
                    let query_before = query.clone();
                    let target_before = target.clone();
                    let raw_matches_before = raw_matches.clone();
                    let target_graph = Vf2GraphRef::target(&target);
                    let builds_before = VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get);

                    let mut actual = Vec::new();
                    let allocation_counts = measure_allocations(|| {
                        actual = project_substruct_matches(&query, target_graph, &raw_matches);
                    });

                    let builds_after = VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get);
                    assert_eq!(
                        builds_after, builds_before,
                        "case {case_index}, goals {goal_count}, reversed {reverse_pairs} built a graph"
                    );
                    assert_eq!(
                        allocation_counts.allocations, EXPECTED_ALLOCATIONS[case_index][goal_count],
                        "case {case_index}, goals {goal_count}, reversed {reverse_pairs} allocation count"
                    );
                    assert_eq!(
                        allocation_counts.reallocations, 0,
                        "case {case_index}, goals {goal_count}, reversed {reverse_pairs} reallocated"
                    );
                    assert_eq!(actual.len(), goal_count);
                    for row_index in 0..goal_count {
                        assert_eq!(
                            actual[row_index].atom_mapping.as_slice(),
                            case.expected_atoms[row_index],
                            "case {case_index}, row {row_index}, reversed {reverse_pairs} atom map"
                        );
                        assert_eq!(
                            actual[row_index].bond_mapping.as_slice(),
                            case.expected_bonds[row_index],
                            "case {case_index}, row {row_index}, reversed {reverse_pairs} bond map"
                        );
                        assert!(
                            actual[row_index]
                                .atom_mapping
                                .iter()
                                .all(|&node| node != NULL_NODE),
                            "all frozen accepted query atoms are mapped"
                        );
                    }
                    assert_eq!(query, query_before, "P01 query was mutated");
                    assert_eq!(target, target_before, "P01 target was mutated");
                    assert_eq!(raw_matches, raw_matches_before, "P01 raw rows were mutated");
                    actual_calls += 1;
                }
            }
        }

        assert_eq!(actual_calls, 24);
    }
}

#[cfg(test)]
mod search_shared_perf_s01_tests {
    use super::search_shared_perf_allocator::measure_allocations;
    use super::{NULL_NODE, NodeId, Vf2Graph, Vf2GraphRef, Vf2SubState};

    #[derive(Clone, Copy)]
    enum MappingShape {
        None,
        FirstOnly,
        All,
    }

    fn empty_graph(n_atoms: usize) -> Vf2Graph {
        Vf2Graph {
            n_atoms,
            n_bonds: 0,
            edge_endpoints: Vec::new(),
            adjacency: vec![Vec::new(); n_atoms],
        }
    }

    fn state_for<'a>(
        query: &'a Vf2Graph,
        target: &'a Vf2Graph,
        shape: MappingShape,
    ) -> Vf2SubState<'a> {
        const FIRST: &[(NodeId, NodeId)] = &[(0, 2)];
        const ALL: &[(NodeId, NodeId)] = &[(0, 2), (1, 0), (2, 1)];

        let mapping = match (query.n_atoms, shape) {
            (0, _) | (1, MappingShape::None) | (3, MappingShape::None) => &[][..],
            (1, MappingShape::FirstOnly | MappingShape::All) | (3, MappingShape::FirstOnly) => {
                FIRST
            }
            (3, MappingShape::All) => ALL,
            _ => unreachable!("S01 fixes query sizes to 0, 1, and 3"),
        };

        let mut state = Vf2SubState::new(
            Vf2GraphRef::compiled(query),
            Vf2GraphRef::compiled(target),
            false,
        );
        for (depth, &(query_index, target_index)) in mapping.iter().enumerate() {
            state.core_1[query_index] = target_index;
            state.core_2[target_index] = query_index;
            state.term_1[query_index] = depth + 1;
            state.term_2[target_index] = depth + 1;
        }
        state.core_len = mapping.len();
        state.t1_len = mapping.len();
        state.t2_len = mapping.len();
        state
    }

    #[derive(Clone, Debug, PartialEq, Eq)]
    struct StateSnapshot {
        core_len: usize,
        t1_len: usize,
        t2_len: usize,
        core_1: Vec<NodeId>,
        core_2: Vec<NodeId>,
        term_1: Vec<usize>,
        term_2: Vec<usize>,
    }

    fn snapshot(state: &Vf2SubState<'_>) -> StateSnapshot {
        StateSnapshot {
            core_len: state.core_len,
            t1_len: state.t1_len,
            t2_len: state.t2_len,
            core_1: state.core_1.clone(),
            core_2: state.core_2.clone(),
            term_1: state.term_1.clone(),
            term_2: state.term_2.clone(),
        }
    }

    #[test]
    fn search_shared_perf_s01_extraction_writes_literal_prefix_without_allocating() {
        let query_graphs = [empty_graph(0), empty_graph(1), empty_graph(3)];
        let target_graph = empty_graph(3);
        let shapes = [
            MappingShape::None,
            MappingShape::FirstOnly,
            MappingShape::All,
        ];
        let mut states = Vec::with_capacity(9);
        for query in &query_graphs {
            for shape in shapes {
                states.push(state_for(query, &target_graph, shape));
            }
        }

        let sentinel = NULL_NODE - 1;
        let mut scratch: Vec<(Vec<NodeId>, Vec<NodeId>)> = states
            .iter()
            .map(|state| (vec![sentinel; state.n1], vec![sentinel; state.n1]))
            .collect();
        for (state, (c1, c2)) in states.iter().zip(&scratch) {
            assert_eq!(c1.capacity(), state.n1);
            assert_eq!(c2.capacity(), state.n1);
        }

        // Literal source-order rows; do not infer expected values from core_1.
        const EXPECTED: [(usize, &[NodeId], &[NodeId]); 9] = [
            (0, &[], &[]),
            (0, &[], &[]),
            (0, &[], &[]),
            (0, &[], &[]),
            (1, &[0], &[2]),
            (1, &[0], &[2]),
            (0, &[], &[]),
            (1, &[0], &[2]),
            (3, &[0, 1, 2], &[2, 0, 1]),
        ];
        let before: Vec<_> = states.iter().map(snapshot).collect();
        let identities: Vec<_> = scratch
            .iter()
            .map(|(c1, c2)| (c1.as_ptr(), c1.capacity(), c2.as_ptr(), c2.capacity()))
            .collect();
        let mut written = [usize::MAX; 9];
        let mut calls = 0;

        let allocations = measure_allocations(|| {
            for index in 0..9 {
                let (c1, c2) = &mut scratch[index];
                written[index] = states[index].get_core_set_into(c1, c2);
                calls += 1;
            }
        });

        assert_eq!(calls, 9);
        assert_eq!(allocations.allocations, 0);
        assert_eq!(allocations.reallocations, 0);
        for index in 0..9 {
            let (expected_len, expected_c1, expected_c2) = EXPECTED[index];
            assert_eq!(written[index], expected_len, "case {index}");
            assert_eq!(
                &scratch[index].0[..expected_len],
                expected_c1,
                "case {index}"
            );
            assert_eq!(
                &scratch[index].1[..expected_len],
                expected_c2,
                "case {index}"
            );
            assert!(
                scratch[index].0[expected_len..]
                    .iter()
                    .all(|&value| value == sentinel)
            );
            assert!(
                scratch[index].1[expected_len..]
                    .iter()
                    .all(|&value| value == sentinel)
            );

            let (c1_ptr, c1_capacity, c2_ptr, c2_capacity) = identities[index];
            assert_eq!(scratch[index].0.as_ptr(), c1_ptr, "case {index}");
            assert_eq!(scratch[index].0.capacity(), c1_capacity, "case {index}");
            assert_eq!(scratch[index].1.as_ptr(), c2_ptr, "case {index}");
            assert_eq!(scratch[index].1.capacity(), c2_capacity, "case {index}");
            assert_eq!(snapshot(&states[index]), before[index], "case {index}");
        }
    }
}

#[cfg(test)]
mod search_shared_perf_s02_tests {
    use super::{NodeId, Vf2Graph, Vf2GraphRef, vf2_entry_one};

    type Mapping = &'static [(NodeId, NodeId)];
    type Trace = &'static [Mapping];

    const ONE_FIRST: Trace = &[&[(0, 0)]];
    const ONE_ALL: Trace = &[&[(0, 0)], &[(0, 1)], &[(0, 2)]];
    const TWO_FIRST: Trace = &[&[(0, 0), (1, 1)]];
    const TWO_ALL: Trace = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
        &[(0, 2), (1, 1)],
    ];
    const TWO_THROUGH_FIRST_TARGET_TWO: Trace = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
    ];

    fn empty_graph(n_atoms: usize) -> Vf2Graph {
        Vf2Graph {
            n_atoms,
            n_bonds: 0,
            edge_endpoints: Vec::new(),
            adjacency: vec![Vec::new(); n_atoms],
        }
    }

    fn expected_case(query_atoms: usize, policy: usize) -> (Trace, Option<Mapping>) {
        match (query_atoms, policy) {
            (1, 0) => (ONE_FIRST, Some(&[(0, 0)])),
            (1, 1) => (ONE_ALL, None),
            (1, 2) => (ONE_ALL, Some(&[(0, 2)])),
            (1, 3) => (ONE_FIRST, Some(&[(0, 0)])),
            (2, 0) => (TWO_FIRST, Some(&[(0, 0), (1, 1)])),
            (2, 1) => (TWO_ALL, None),
            (2, 2) => (TWO_THROUGH_FIRST_TARGET_TWO, Some(&[(0, 2), (1, 0)])),
            (2, 3) => (TWO_FIRST, Some(&[(0, 0), (1, 1)])),
            _ => unreachable!("S02 freezes query sizes one/two and four policies"),
        }
    }

    #[test]
    fn search_shared_perf_s02_first_match_preserves_literal_callback_order() {
        let query_graphs = [empty_graph(1), empty_graph(2)];
        let target_graph = empty_graph(3);
        let atom_fn = |_: usize, _: usize| true;
        let bond_fn = |_: usize, _: usize| true;
        let mut actual_calls = 0;

        for (query_index, query) in query_graphs.iter().enumerate() {
            let query_atoms = query_index + 1;
            for policy in 0..4 {
                let (expected_trace, expected_result) = expected_case(query_atoms, policy);
                let mut trace: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
                let mut goal_ordinal = 0;
                let mut result = Vec::new();
                let found = {
                    let mut callback = |c1: &[NodeId], c2: &[NodeId]| {
                        trace.push(c1.iter().copied().zip(c2.iter().copied()).collect());
                        let ordinal = goal_ordinal;
                        goal_ordinal += 1;
                        match policy {
                            0 => true,
                            1 => false,
                            2 => c2.first() == Some(&2),
                            3 => ordinal % 2 == 0,
                            _ => unreachable!("policy was frozen to four cases"),
                        }
                    };
                    vf2_entry_one(
                        Vf2GraphRef::compiled(query),
                        Vf2GraphRef::compiled(&target_graph),
                        &atom_fn,
                        &bond_fn,
                        Some(&mut callback),
                        &mut result,
                    )
                };
                actual_calls += 1;

                assert_eq!(
                    trace.len(),
                    expected_trace.len(),
                    "query/policy {query_atoms}/{policy}"
                );
                for (observed, expected) in trace.iter().zip(expected_trace) {
                    assert_eq!(
                        observed.as_slice(),
                        *expected,
                        "query/policy {query_atoms}/{policy}"
                    );
                }
                assert_eq!(
                    found,
                    expected_result.is_some(),
                    "query/policy {query_atoms}/{policy}"
                );
                assert_eq!(
                    result.as_slice(),
                    expected_result.unwrap_or(&[]),
                    "query/policy {query_atoms}/{policy}"
                );
            }
        }

        assert_eq!(actual_calls, 8);
    }
}

#[cfg(test)]
mod search_shared_perf_s03_tests {
    use super::search_shared_perf_allocator::measure_allocations;
    use super::{NULL_NODE, NodeId, Vf2Graph, Vf2GraphRef, Vf2SubState, vf2_match_all};

    type Mapping = &'static [(NodeId, NodeId)];
    type Sequences = &'static [Mapping];

    const ONE_GOALS: Sequences = &[&[(0, 0)], &[(0, 1)], &[(0, 2)]];
    const TWO_GOALS: Sequences = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
        &[(0, 2), (1, 1)],
    ];

    const NO_SEQUENCES: Sequences = &[];
    const ONE_TARGET_TWO: Sequences = &[&[(0, 2)]];
    const ONE_EVEN_ORDINALS: Sequences = &[&[(0, 0)], &[(0, 2)]];
    const TWO_TARGET_TWO: Sequences = &[&[(0, 2), (1, 0)], &[(0, 2), (1, 1)]];
    const TWO_EVEN_ORDINALS: Sequences = &[&[(0, 0), (1, 1)], &[(0, 1), (1, 0)], &[(0, 2), (1, 0)]];

    // Columns follow the literal limit list [0, 1, 2, 9]; rows are the four
    // callback policies in their frozen order.
    const ONE_TRACE_COUNTS: [[usize; 4]; 4] =
        [[3, 1, 2, 3], [3, 3, 3, 3], [3, 3, 3, 3], [3, 1, 3, 3]];
    const TWO_TRACE_COUNTS: [[usize; 4]; 4] =
        [[6, 1, 2, 6], [6, 6, 6, 6], [6, 5, 6, 6], [6, 1, 3, 6]];

    fn empty_graph(n_atoms: usize) -> Vf2Graph {
        Vf2Graph {
            n_atoms,
            n_bonds: 0,
            edge_endpoints: Vec::new(),
            adjacency: vec![Vec::new(); n_atoms],
        }
    }

    fn expected_sequences(query_index: usize, policy: usize) -> Sequences {
        match (query_index, policy) {
            (0, 0) => ONE_GOALS,
            (0, 1) => NO_SEQUENCES,
            (0, 2) => ONE_TARGET_TWO,
            (0, 3) => ONE_EVEN_ORDINALS,
            (1, 0) => TWO_GOALS,
            (1, 1) => NO_SEQUENCES,
            (1, 2) => TWO_TARGET_TWO,
            (1, 3) => TWO_EVEN_ORDINALS,
            _ => unreachable!("S03 freezes two query shapes and four policies"),
        }
    }

    #[test]
    fn search_shared_perf_s03_limits_order_backtracking_and_reject_allocations() {
        let query_graphs = [empty_graph(1), empty_graph(2)];
        let target_graph = empty_graph(3);
        let atom_fn = |_: usize, _: usize| true;
        let bond_fn = |_: usize, _: usize| true;
        let limits = [0, 1, 2, 9];
        let mut actual_calls = 0;
        let mut measured_reject_calls = 0;

        for (query_index, query) in query_graphs.iter().enumerate() {
            let query_atoms = query_index + 1;
            let all_goals = if query_index == 0 {
                ONE_GOALS
            } else {
                TWO_GOALS
            };
            let trace_counts = if query_index == 0 {
                ONE_TRACE_COUNTS
            } else {
                TWO_TRACE_COUNTS
            };

            for policy in 0..4 {
                let accepted_sequences = expected_sequences(query_index, policy);
                for (limit_index, limit) in limits.into_iter().enumerate() {
                    let mut state = Vf2SubState::new(
                        Vf2GraphRef::compiled(query),
                        Vf2GraphRef::compiled(&target_graph),
                        false,
                    );
                    let mut c1 = vec![NULL_NODE; query_atoms];
                    let mut c2 = vec![NULL_NODE; query_atoms];
                    let mut results: Vec<Vec<(NodeId, NodeId)>> = Vec::with_capacity(6);
                    let mut trace = [[(NULL_NODE, NULL_NODE); 2]; 6];
                    let mut trace_len = 0;
                    let mut goal_ordinal = 0;
                    let mut callback = |mapped_query: &[NodeId], mapped_target: &[NodeId]| {
                        for pair_index in 0..mapped_query.len() {
                            trace[trace_len][pair_index] =
                                (mapped_query[pair_index], mapped_target[pair_index]);
                        }
                        trace_len += 1;
                        let ordinal = goal_ordinal;
                        goal_ordinal += 1;
                        match policy {
                            0 => true,
                            1 => false,
                            2 => mapped_target.first() == Some(&2),
                            3 => ordinal % 2 == 0,
                            _ => unreachable!("policy was frozen to four cases"),
                        }
                    };
                    let measure_reject_case = query_index == 1 && policy == 1 && limit == 0;
                    let mut found = false;
                    let allocation_counts = if measure_reject_case {
                        measured_reject_calls += 1;
                        Some(measure_allocations(|| {
                            found = vf2_match_all(
                                &mut state,
                                &atom_fn,
                                &bond_fn,
                                Some(&mut callback),
                                &mut c1,
                                &mut c2,
                                &mut results,
                                limit,
                            );
                        }))
                    } else {
                        found = vf2_match_all(
                            &mut state,
                            &atom_fn,
                            &bond_fn,
                            Some(&mut callback),
                            &mut c1,
                            &mut c2,
                            &mut results,
                            limit,
                        );
                        None
                    };
                    actual_calls += 1;

                    let expected_trace_len = trace_counts[policy][limit_index];
                    let expected_result_count = if limit == 0 {
                        accepted_sequences.len()
                    } else {
                        accepted_sequences.len().min(limit)
                    };
                    assert_eq!(
                        trace_len, expected_trace_len,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    for (observed, expected) in trace[..trace_len]
                        .iter()
                        .zip(all_goals.iter().take(expected_trace_len))
                    {
                        assert_eq!(
                            &observed[..query_atoms],
                            *expected,
                            "query/policy/limit {query_atoms}/{policy}/{limit}"
                        );
                    }
                    assert_eq!(
                        found,
                        expected_result_count != 0,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    assert_eq!(
                        results.len(),
                        expected_result_count,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    for (observed, expected) in results
                        .iter()
                        .zip(accepted_sequences.iter().take(expected_result_count))
                    {
                        assert_eq!(
                            observed.as_slice(),
                            *expected,
                            "query/policy/limit {query_atoms}/{policy}/{limit}"
                        );
                    }

                    let reached_positive_limit = limit > 0 && expected_result_count >= limit;
                    if !reached_positive_limit {
                        assert_eq!(state.core_len, 0);
                        assert_eq!(state.t1_len, 0);
                        assert_eq!(state.t2_len, 0);
                        assert!(state.core_1.iter().all(|&value| value == NULL_NODE));
                        assert!(state.core_2.iter().all(|&value| value == NULL_NODE));
                        assert!(state.term_1.iter().all(|&value| value == 0));
                        assert!(state.term_2.iter().all(|&value| value == 0));
                    }

                    if measure_reject_case {
                        let counts = allocation_counts.expect("reject-all case is measured");
                        assert_eq!(counts.allocations, 0);
                        assert_eq!(counts.reallocations, 0);
                    } else {
                        assert!(allocation_counts.is_none());
                    }
                }
            }
        }

        assert_eq!(actual_calls, 32);
        assert_eq!(measured_reject_calls, 1);
    }
}

#[cfg(test)]
mod search_shared_perf_s04_tests {
    use super::search_shared_perf_allocator::measure_allocations;
    use super::{Vf2Graph, Vf2GraphRef, Vf2NeighborRow, build_vf2_graph};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, QueryAtom, QueryBond, QueryGraph,
        TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    const MAX_ATOMS: usize = 4;
    const MAX_BONDS: usize = 3;
    const MAX_DEGREE: usize = 3;
    const UNUSED_PAIR: (usize, usize) = (usize::MAX, usize::MAX);

    struct FixedShape {
        atom_count: usize,
        edges: &'static [(usize, usize)],
        neighbors: &'static [&'static [(usize, usize)]],
    }

    const EMPTY: FixedShape = FixedShape {
        atom_count: 0,
        edges: &[],
        neighbors: &[],
    };
    const CHAIN: FixedShape = FixedShape {
        atom_count: 3,
        edges: &[(0, 1), (1, 2)],
        neighbors: &[&[(1, 0)], &[(0, 0), (2, 1)], &[(1, 1)]],
    };
    const BRANCHED: FixedShape = FixedShape {
        atom_count: 4,
        edges: &[(0, 1), (0, 2), (0, 3)],
        neighbors: &[&[(1, 0), (2, 1), (3, 2)], &[(0, 0)], &[(0, 1)], &[(0, 2)]],
    };
    const TRIANGLE: FixedShape = FixedShape {
        atom_count: 3,
        edges: &[(0, 1), (1, 2), (2, 0)],
        neighbors: &[&[(1, 0), (2, 2)], &[(0, 0), (2, 1)], &[(1, 1), (0, 2)]],
    };
    const SHAPES: [FixedShape; 4] = [EMPTY, CHAIN, BRANCHED, TRIANGLE];

    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    struct ViewObservation {
        atom_count: usize,
        bond_count: usize,
        degrees: [usize; MAX_ATOMS],
        row_lengths: [usize; MAX_ATOMS],
        iterator_lengths_before: [usize; MAX_ATOMS],
        iterator_lengths_after: [usize; MAX_ATOMS],
        iterated: [[(usize, usize); MAX_DEGREE]; MAX_ATOMS],
        indexed: [[(usize, usize); MAX_DEGREE]; MAX_ATOMS],
        row_pointers: [*const (); MAX_ATOMS],
        endpoints: [(usize, usize); MAX_BONDS],
    }

    impl Default for ViewObservation {
        fn default() -> Self {
            Self {
                atom_count: 0,
                bond_count: 0,
                degrees: [0; MAX_ATOMS],
                row_lengths: [0; MAX_ATOMS],
                iterator_lengths_before: [0; MAX_ATOMS],
                iterator_lengths_after: [0; MAX_ATOMS],
                iterated: [[UNUSED_PAIR; MAX_DEGREE]; MAX_ATOMS],
                indexed: [[UNUSED_PAIR; MAX_DEGREE]; MAX_ATOMS],
                row_pointers: [std::ptr::null(); MAX_ATOMS],
                endpoints: [UNUSED_PAIR; MAX_BONDS],
            }
        }
    }

    fn atoms(count: usize) -> Vec<Atom> {
        (0..count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect()
    }

    fn bonds(shape: &FixedShape) -> Vec<Bond> {
        shape
            .edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect()
    }

    fn query(shape: &FixedShape) -> QueryGraph {
        let query_atoms = (0..shape.atom_count)
            .map(|index| QueryAtom::new(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let query_bonds = shape
            .edges
            .iter()
            .enumerate()
            .map(|(index, &(begin, end))| {
                QueryBond::new(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        QueryGraph::from_parts(
            query_atoms,
            query_bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed S04 query graph is valid")
    }

    fn target(shape: &FixedShape) -> TopologyBlock {
        TopologyBlock::try_from_parts(
            atoms(shape.atom_count),
            bonds(shape),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed S04 target topology is valid")
    }

    fn make_view<'a>(
        shape_index: usize,
        representation: usize,
        queries: &'a [QueryGraph; 4],
        targets: &'a [TopologyBlock; 4],
        compiled: &'a [Vf2Graph; 4],
    ) -> Vf2GraphRef<'a> {
        match representation {
            0 => Vf2GraphRef::query(&queries[shape_index]),
            1 => Vf2GraphRef::target(&targets[shape_index]),
            2 => Vf2GraphRef::compiled(&compiled[shape_index]),
            _ => unreachable!("S04 freezes Query, Target and Compiled views"),
        }
    }

    fn borrowed_row_pointer(row: Vf2NeighborRow<'_>) -> *const () {
        match row {
            Vf2NeighborRow::Pairs(row) => row.as_ptr().cast(),
            Vf2NeighborRow::NeighborRefs(row) => row.as_ptr().cast(),
        }
    }

    fn owner_row_pointer(view: Vf2GraphRef<'_>, node: usize) -> *const () {
        match view {
            Vf2GraphRef::Query(graph) => graph.adjacency()[node].as_ptr().cast(),
            Vf2GraphRef::Target(graph) => graph.adjacency.neighbors_of(node).as_ptr().cast(),
            Vf2GraphRef::Compiled(graph) => graph.adjacency[node].as_ptr().cast(),
        }
    }

    fn observe(view: Vf2GraphRef<'_>, observation: &mut ViewObservation) {
        observation.atom_count = view.num_atoms();
        observation.bond_count = view.num_bonds();
        for node in 0..observation.atom_count {
            let row = view.neighbor_row(node);
            observation.degrees[node] = view.out_degree(node);
            observation.row_lengths[node] = row.len();
            observation.row_pointers[node] = borrowed_row_pointer(row);

            let mut iterator = row.iter();
            observation.iterator_lengths_before[node] = iterator.len();
            let mut index = 0;
            while let Some(pair) = iterator.next() {
                observation.iterated[node][index] = pair;
                index += 1;
            }
            observation.iterator_lengths_after[node] = iterator.len();

            for index in 0..row.len() {
                observation.indexed[node][index] = row.get(index).unwrap_or(UNUSED_PAIR);
            }
        }
        for edge in 0..observation.bond_count {
            observation.endpoints[edge] = view.bond_endpoints(edge);
        }
    }

    #[test]
    fn search_shared_perf_s04_borrowed_views_match_literal_owners_without_allocation() {
        let queries: [QueryGraph; 4] = std::array::from_fn(|index| query(&SHAPES[index]));
        let targets: [TopologyBlock; 4] = std::array::from_fn(|index| target(&SHAPES[index]));
        let compiled: [Vf2Graph; 4] = std::array::from_fn(|index| build_vf2_graph(&queries[index]));
        let queries_before = queries.clone();
        let targets_before = targets.clone();
        let compiled_before = compiled.clone();

        let mut expected_row_pointers = [[std::ptr::null(); MAX_ATOMS]; 12];
        for shape_index in 0..SHAPES.len() {
            for representation in 0..3 {
                let case_index = shape_index * 3 + representation;
                let view = make_view(shape_index, representation, &queries, &targets, &compiled);
                for node in 0..SHAPES[shape_index].atom_count {
                    expected_row_pointers[case_index][node] = owner_row_pointer(view, node);
                }
            }
        }

        let mut observations = [ViewObservation::default(); 12];
        let mut actual_views = 0;
        let allocations = measure_allocations(|| {
            for shape_index in 0..SHAPES.len() {
                for representation in 0..3 {
                    let view =
                        make_view(shape_index, representation, &queries, &targets, &compiled);
                    observe(view, &mut observations[actual_views]);
                    actual_views += 1;
                }
            }
        });

        assert_eq!(actual_views, 12);
        assert_eq!(allocations.allocations, 0);
        assert_eq!(allocations.reallocations, 0);
        assert_eq!(queries, queries_before);
        assert_eq!(targets, targets_before);
        assert_eq!(compiled, compiled_before);

        for (case_index, observation) in observations.iter().enumerate() {
            let shape_index = case_index / 3;
            let shape = &SHAPES[shape_index];
            assert_eq!(
                observation.atom_count, shape.atom_count,
                "case {case_index}"
            );
            assert_eq!(
                observation.bond_count,
                shape.edges.len(),
                "case {case_index}"
            );
            for node in 0..shape.atom_count {
                let expected_row = shape.neighbors[node];
                let degree = expected_row.len();
                assert_eq!(
                    observation.degrees[node], degree,
                    "case/node {case_index}/{node}"
                );
                assert_eq!(
                    observation.row_lengths[node], degree,
                    "case/node {case_index}/{node}"
                );
                assert_eq!(
                    observation.iterator_lengths_before[node], degree,
                    "case/node {case_index}/{node}"
                );
                assert_eq!(
                    observation.iterator_lengths_after[node], 0,
                    "case/node {case_index}/{node}"
                );
                assert_eq!(
                    &observation.iterated[node][..degree],
                    expected_row,
                    "ordered iterator case/node {case_index}/{node}"
                );
                assert_eq!(
                    &observation.indexed[node][..degree],
                    expected_row,
                    "indexed row case/node {case_index}/{node}"
                );
                assert!(
                    observation.iterated[node][degree..]
                        .iter()
                        .all(|pair| *pair == UNUSED_PAIR)
                );
                assert!(
                    observation.indexed[node][degree..]
                        .iter()
                        .all(|pair| *pair == UNUSED_PAIR)
                );
                assert_eq!(
                    observation.row_pointers[node], expected_row_pointers[case_index][node],
                    "borrowed identity case/node {case_index}/{node}"
                );

                for &(neighbor, edge) in expected_row {
                    let (begin, end) = observation.endpoints[edge];
                    assert!(
                        (begin == node && end == neighbor) || (begin == neighbor && end == node),
                        "endpoint orientations case/node/edge {case_index}/{node}/{edge}"
                    );
                }
            }
            assert_eq!(
                &observation.endpoints[..shape.edges.len()],
                shape.edges,
                "literal endpoint order case {case_index}"
            );
            assert!(
                observation.endpoints[shape.edges.len()..]
                    .iter()
                    .all(|pair| *pair == UNUSED_PAIR)
            );
        }
    }
}

#[cfg(test)]
mod search_shared_perf_s05_tests {
    use super::{NULL_NODE, NodeId, Vf2GraphRef, Vf2SubState, vf2_entry_one, vf2_match_all};
    use cosmolkit_model::{Atom, AtomId, AtomSpec, QueryAtom, QueryGraph, TopologyBlock};
    use cosmolkit_types::Element;

    type Mapping = &'static [(NodeId, NodeId)];
    type Sequences = &'static [Mapping];
    type Trace = &'static [Mapping];

    const ONE_GOALS: Sequences = &[&[(0, 0)], &[(0, 1)], &[(0, 2)]];
    const TWO_GOALS: Sequences = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
        &[(0, 2), (1, 1)],
    ];
    const NO_SEQUENCES: Sequences = &[];
    const ONE_TARGET_TWO: Sequences = &[&[(0, 2)]];
    const ONE_EVEN_ORDINALS: Sequences = &[&[(0, 0)], &[(0, 2)]];
    const TWO_TARGET_TWO: Sequences = &[&[(0, 2), (1, 0)], &[(0, 2), (1, 1)]];
    const TWO_EVEN_ORDINALS: Sequences = &[&[(0, 0), (1, 1)], &[(0, 1), (1, 0)], &[(0, 2), (1, 0)]];

    const ONE_FIRST: Trace = &[&[(0, 0)]];
    const ONE_ALL: Trace = &[&[(0, 0)], &[(0, 1)], &[(0, 2)]];
    const TWO_FIRST: Trace = &[&[(0, 0), (1, 1)]];
    const TWO_ALL: Trace = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
        &[(0, 2), (1, 1)],
    ];
    const TWO_THROUGH_FIRST_TARGET_TWO: Trace = &[
        &[(0, 0), (1, 1)],
        &[(0, 0), (1, 2)],
        &[(0, 1), (1, 0)],
        &[(0, 1), (1, 2)],
        &[(0, 2), (1, 0)],
    ];

    const LIMITS: [usize; 4] = [0, 1, 2, 9];
    const TRACE_COUNTS: [[[usize; 4]; 4]; 2] = [
        [[3, 1, 2, 3], [3, 3, 3, 3], [3, 3, 3, 3], [3, 1, 3, 3]],
        [[6, 1, 2, 6], [6, 6, 6, 6], [6, 5, 6, 6], [6, 1, 3, 6]],
    ];

    fn empty_query(atom_count: usize) -> QueryGraph {
        let atoms = (0..atom_count)
            .map(|index| QueryAtom::new(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        QueryGraph::from_parts(
            atoms,
            Vec::new(),
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed S05 query is valid")
    }

    fn empty_target(atom_count: usize) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        TopologyBlock::try_from_parts(atoms, Vec::new(), Vec::new(), Vec::new())
            .expect("fixed S05 target topology is valid")
    }

    fn first_expected(query_index: usize, policy: usize) -> (Trace, Option<Mapping>) {
        match (query_index, policy) {
            (0, 0) => (ONE_FIRST, Some(&[(0, 0)])),
            (0, 1) => (ONE_ALL, None),
            (0, 2) => (ONE_ALL, Some(&[(0, 2)])),
            (0, 3) => (ONE_FIRST, Some(&[(0, 0)])),
            (1, 0) => (TWO_FIRST, Some(&[(0, 0), (1, 1)])),
            (1, 1) => (TWO_ALL, None),
            (1, 2) => (TWO_THROUGH_FIRST_TARGET_TWO, Some(&[(0, 2), (1, 0)])),
            (1, 3) => (TWO_FIRST, Some(&[(0, 0), (1, 1)])),
            _ => unreachable!("S05 freezes two query shapes and four policies"),
        }
    }

    fn all_expected(query_index: usize, policy: usize) -> Sequences {
        match (query_index, policy) {
            (0, 0) => ONE_GOALS,
            (0, 1) => NO_SEQUENCES,
            (0, 2) => ONE_TARGET_TWO,
            (0, 3) => ONE_EVEN_ORDINALS,
            (1, 0) => TWO_GOALS,
            (1, 1) => NO_SEQUENCES,
            (1, 2) => TWO_TARGET_TWO,
            (1, 3) => TWO_EVEN_ORDINALS,
            _ => unreachable!("S05 freezes two query shapes and four policies"),
        }
    }

    // These are the same 8 first-match and 32 all-match literal cases as
    // S02/S03, rerun through QueryGraph and TopologyBlock borrowed storage.
    #[test]
    fn search_shared_perf_s05_first_match_uses_borrowed_query_target_views() {
        let query_graphs = [empty_query(1), empty_query(2)];
        let target_graph = empty_target(3);
        let query_graphs_before = query_graphs.clone();
        let target_graph_before = target_graph.clone();
        let atom_fn = |_: usize, _: usize| true;
        let bond_fn = |_: usize, _: usize| true;
        let mut actual_calls = 0;

        for (query_index, query) in query_graphs.iter().enumerate() {
            let query_view = Vf2GraphRef::query(query);
            let target_view = Vf2GraphRef::target(&target_graph);
            for policy in 0..4 {
                let (expected_trace, expected_result) = first_expected(query_index, policy);
                let mut trace: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
                let mut goal_ordinal = 0;
                let mut result = Vec::new();
                let found = {
                    let mut callback = |c1: &[NodeId], c2: &[NodeId]| {
                        trace.push(c1.iter().copied().zip(c2.iter().copied()).collect());
                        let ordinal = goal_ordinal;
                        goal_ordinal += 1;
                        match policy {
                            0 => true,
                            1 => false,
                            2 => c2.first() == Some(&2),
                            3 => ordinal % 2 == 0,
                            _ => unreachable!("policy was frozen to four cases"),
                        }
                    };
                    vf2_entry_one(
                        query_view,
                        target_view,
                        &atom_fn,
                        &bond_fn,
                        Some(&mut callback),
                        &mut result,
                    )
                };
                actual_calls += 1;

                assert_eq!(
                    trace.len(),
                    expected_trace.len(),
                    "query/policy {}/{policy}",
                    query_index + 1
                );
                for (observed, expected) in trace.iter().zip(expected_trace) {
                    assert_eq!(
                        observed.as_slice(),
                        *expected,
                        "query/policy {}/{policy}",
                        query_index + 1
                    );
                }
                assert_eq!(
                    found,
                    expected_result.is_some(),
                    "query/policy {}/{policy}",
                    query_index + 1
                );
                assert_eq!(
                    result.as_slice(),
                    expected_result.unwrap_or(&[]),
                    "query/policy {}/{policy}",
                    query_index + 1
                );
            }
        }

        assert_eq!(actual_calls, 8);
        assert_eq!(query_graphs, query_graphs_before);
        assert_eq!(target_graph, target_graph_before);
    }

    #[test]
    fn search_shared_perf_s05_all_match_uses_borrowed_query_target_views() {
        let query_graphs = [empty_query(1), empty_query(2)];
        let target_graph = empty_target(3);
        let query_graphs_before = query_graphs.clone();
        let target_graph_before = target_graph.clone();
        let atom_fn = |_: usize, _: usize| true;
        let bond_fn = |_: usize, _: usize| true;
        let mut actual_calls = 0;

        for (query_index, query) in query_graphs.iter().enumerate() {
            let query_atoms = query_index + 1;
            let query_view = Vf2GraphRef::query(query);
            let target_view = Vf2GraphRef::target(&target_graph);
            let all_goals = if query_index == 0 {
                ONE_GOALS
            } else {
                TWO_GOALS
            };

            for policy in 0..4 {
                let accepted_sequences = all_expected(query_index, policy);
                for (limit_index, limit) in LIMITS.into_iter().enumerate() {
                    let mut state = Vf2SubState::new(query_view, target_view, false);
                    let mut c1 = vec![NULL_NODE; query_atoms];
                    let mut c2 = vec![NULL_NODE; query_atoms];
                    let mut results: Vec<Vec<(NodeId, NodeId)>> = Vec::new();
                    let mut trace = [[(NULL_NODE, NULL_NODE); 2]; 6];
                    let mut trace_len = 0;
                    let mut goal_ordinal = 0;
                    let mut callback = |mapped_query: &[NodeId], mapped_target: &[NodeId]| {
                        for pair_index in 0..mapped_query.len() {
                            trace[trace_len][pair_index] =
                                (mapped_query[pair_index], mapped_target[pair_index]);
                        }
                        trace_len += 1;
                        let ordinal = goal_ordinal;
                        goal_ordinal += 1;
                        match policy {
                            0 => true,
                            1 => false,
                            2 => mapped_target.first() == Some(&2),
                            3 => ordinal % 2 == 0,
                            _ => unreachable!("policy was frozen to four cases"),
                        }
                    };
                    let found = vf2_match_all(
                        &mut state,
                        &atom_fn,
                        &bond_fn,
                        Some(&mut callback),
                        &mut c1,
                        &mut c2,
                        &mut results,
                        limit,
                    );
                    actual_calls += 1;

                    let expected_trace_len = TRACE_COUNTS[query_index][policy][limit_index];
                    let expected_result_count = if limit == 0 {
                        accepted_sequences.len()
                    } else {
                        accepted_sequences.len().min(limit)
                    };
                    assert_eq!(
                        trace_len, expected_trace_len,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    for (observed, expected) in trace[..trace_len]
                        .iter()
                        .zip(all_goals.iter().take(expected_trace_len))
                    {
                        assert_eq!(
                            &observed[..query_atoms],
                            *expected,
                            "query/policy/limit {query_atoms}/{policy}/{limit}"
                        );
                    }
                    assert_eq!(
                        found,
                        expected_result_count != 0,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    assert_eq!(
                        results.len(),
                        expected_result_count,
                        "query/policy/limit {query_atoms}/{policy}/{limit}"
                    );
                    for (observed, expected) in results
                        .iter()
                        .zip(accepted_sequences.iter().take(expected_result_count))
                    {
                        assert_eq!(
                            observed.as_slice(),
                            *expected,
                            "query/policy/limit {query_atoms}/{policy}/{limit}"
                        );
                    }

                    let reached_positive_limit = limit > 0 && expected_result_count >= limit;
                    if !reached_positive_limit {
                        assert_eq!(state.core_len, 0);
                        assert_eq!(state.t1_len, 0);
                        assert_eq!(state.t2_len, 0);
                        assert!(state.core_1.iter().all(|&value| value == NULL_NODE));
                        assert!(state.core_2.iter().all(|&value| value == NULL_NODE));
                        assert!(state.term_1.iter().all(|&value| value == 0));
                        assert!(state.term_2.iter().all(|&value| value == 0));
                    }
                }
            }
        }

        assert_eq!(actual_calls, 32);
        assert_eq!(query_graphs, query_graphs_before);
        assert_eq!(target_graph, target_graph_before);
    }
}

#[cfg(test)]
mod search_shared_perf_s06_tests {
    use super::{
        SubstructMatchParams, VF2_GRAPH_BUILD_ENTRIES, compile_query_graph,
        get_substruct_matches_with_compiled_query, substruct_match_impl,
        try_get_substruct_matches_with_params_and_context,
    };
    use crate::{
        SearchTarget, SearchTargetAccess, SmartsParseParams, build_prepared_query_match_context,
        parse_smarts,
    };
    use cosmolkit_core::{ValenceModel, assign_valence_with_options_for_topology, fast_find_rings};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, QueryGraph, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    type ExpectedRows = &'static [(&'static [usize], &'static [usize])];

    const Q76_EXPECTED: ExpectedRows = &[
        (&[1, 2, 3], &[1, 2]),
        (&[4, 5, 6], &[4, 5]),
        (&[5, 4, 3], &[4, 3]),
    ];
    const NO_MATCHES: ExpectedRows = &[];
    const DISCONNECTED_EXPECTED: ExpectedRows = &[
        (&[0, 1], &[]),
        (&[0, 2], &[]),
        (&[1, 0], &[]),
        (&[1, 2], &[]),
        (&[2, 0], &[]),
        (&[2, 1], &[]),
    ];
    const RECURSIVE_EXPECTED: ExpectedRows = &[(&[1], &[])];
    const LIMITS: [usize; 4] = [0, 1, 2, 9];

    struct Fixture {
        label: &'static str,
        query: QueryGraph,
        topology: TopologyBlock,
        expected: ExpectedRows,
        uniquify: bool,
    }

    fn topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed S06 target topology is valid")
    }

    fn query(smarts: &str) -> QueryGraph {
        parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("fixed S06 query {smarts:?} parses: {error}"))
    }

    fn graph_build_count() -> usize {
        VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get)
    }

    #[test]
    fn search_shared_perf_s06_canonical_routes_borrow_prepared_graphs() {
        VF2_GRAPH_BUILD_ENTRIES.with(|entries| entries.set(0));

        let q76_target = topology(
            &[
                Element::C,
                Element::C,
                Element::C,
                Element::O,
                Element::C,
                Element::C,
                Element::O,
            ],
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
        );
        let disconnected_target = topology(&[Element::C, Element::C, Element::C], &[]);
        let cco_target = topology(&[Element::C, Element::C, Element::O], &[(0, 1), (1, 2)]);
        let fixtures = vec![
            Fixture {
                label: "Q76 C-C-O",
                query: query("C-C-O"),
                topology: q76_target.clone(),
                expected: Q76_EXPECTED,
                uniquify: true,
            },
            Fixture {
                label: "Q76 N no-match",
                query: query("N"),
                topology: q76_target,
                expected: NO_MATCHES,
                uniquify: true,
            },
            Fixture {
                label: "disconnected C.C",
                query: query("C.C"),
                topology: disconnected_target,
                expected: DISCONNECTED_EXPECTED,
                uniquify: false,
            },
            Fixture {
                label: "recursive rooted carbon",
                query: query("[C;$([C]-[O])]"),
                topology: cco_target,
                expected: RECURSIVE_EXPECTED,
                uniquify: true,
            },
        ];
        let compiled_graphs = fixtures
            .iter()
            .map(|fixture| compile_query_graph(&fixture.query))
            .collect::<Vec<_>>();
        assert_eq!(graph_build_count(), fixtures.len());

        let mut actual_calls = 0;
        for (fixture_index, fixture) in fixtures.iter().enumerate() {
            let query_before = fixture.query.clone();
            let topology_before = fixture.topology.clone();
            let coordinates = CoordinateBlock::default();
            let coordinates_before = coordinates.clone();
            let rings = fast_find_rings(&fixture.topology).expect("fixed S06 rings prepare");
            let rings_before = rings.clone();
            let valence = assign_valence_with_options_for_topology(
                &fixture.topology,
                ValenceModel::RdkitLike,
                false,
            )
            .expect("fixed S06 valence prepares");
            let valence_before = valence.clone();
            let query_context =
                build_prepared_query_match_context(&fixture.topology, &rings, &valence)
                    .expect("fixed S06 prepared query context validates");
            let target = SearchTarget::new(
                &fixture.topology,
                &coordinates,
                &fixture.topology.stereo_groups,
                Some(&rings),
                Some(&valence),
            );

            for limit in LIMITS {
                let params = SubstructMatchParams {
                    max_matches: limit,
                    uniquify: fixture.uniquify,
                    recursion_possible: true,
                    ..SubstructMatchParams::default()
                };
                for route in 0..3 {
                    let builds_before = graph_build_count();
                    let result = match route {
                        0 => substruct_match_impl(&target, &fixture.query, &params),
                        1 => try_get_substruct_matches_with_params_and_context(
                            &target,
                            &fixture.query,
                            &params,
                            &query_context,
                        ),
                        2 => get_substruct_matches_with_compiled_query(
                            &target,
                            &fixture.query,
                            &params,
                            &compiled_graphs[fixture_index],
                        ),
                        _ => unreachable!("S06 freezes ordinary, prepared and compiled routes"),
                    };
                    let builds_after = graph_build_count();
                    assert_eq!(
                        builds_after, builds_before,
                        "{} route {route} limit {limit} rebuilt a graph",
                        fixture.label
                    );
                    actual_calls += 1;

                    let actual = result.unwrap_or_else(|error| {
                        panic!("{} route {route} limit {limit}: {error}", fixture.label)
                    });
                    let expected_count = if limit == 0 {
                        fixture.expected.len()
                    } else {
                        fixture.expected.len().min(limit)
                    };
                    assert_eq!(
                        actual.len(),
                        expected_count,
                        "{} route {route} limit {limit}",
                        fixture.label
                    );
                    for (result, (expected_atoms, expected_bonds)) in actual
                        .iter()
                        .zip(fixture.expected.iter().take(expected_count))
                    {
                        assert_eq!(
                            result.atom_mapping, *expected_atoms,
                            "{} route {route} limit {limit} atom mapping",
                            fixture.label
                        );
                        assert_eq!(
                            result.bond_mapping, *expected_bonds,
                            "{} route {route} limit {limit} bond mapping",
                            fixture.label
                        );
                    }
                    assert_eq!(
                        fixture.query, query_before,
                        "{} query mutated",
                        fixture.label
                    );
                    assert_eq!(
                        fixture.topology, topology_before,
                        "{} topology mutated",
                        fixture.label
                    );
                    assert_eq!(
                        coordinates, coordinates_before,
                        "{} coordinates mutated",
                        fixture.label
                    );
                    assert_eq!(rings, rings_before, "{} rings mutated", fixture.label);
                    assert_eq!(valence, valence_before, "{} valence mutated", fixture.label);
                }
            }
        }

        assert_eq!(actual_calls, 48);
        assert_eq!(graph_build_count(), fixtures.len());
    }
}

#[cfg(test)]
mod search_shared_perf_s07_tests {
    use super::{
        SubstructMatchParams, VF2_GRAPH_BUILD_ENTRIES,
        try_get_substruct_matches_with_params_and_context,
    };
    use crate::{
        SearchTarget, SmartsParseParams, build_prepared_query_match_context, parse_smarts,
    };
    use cosmolkit_core::{ValenceModel, assign_valence_with_options_for_topology, fast_find_rings};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, CoordinateBlock, QueryGraph, TopologyBlock,
    };
    use cosmolkit_types::{BondOrder, Element};

    type ExpectedRows = &'static [(&'static [usize], &'static [usize])];

    const CARBON_EXPECTED: ExpectedRows = &[(&[0], &[]), (&[1], &[])];
    const OXYGEN_EXPECTED: ExpectedRows = &[(&[2], &[])];
    const NO_MATCHES: ExpectedRows = &[];
    const DISCONNECTED_EXPECTED: ExpectedRows = &[(&[0, 1], &[])];
    const CARBON_OXYGEN_EXPECTED: ExpectedRows = &[(&[1, 2], &[1])];
    const RECURSIVE_EXPECTED: ExpectedRows = &[(&[1], &[])];

    struct FeatureQuery {
        label: &'static str,
        graph: QueryGraph,
        expected: ExpectedRows,
    }

    fn topology(elements: &[Element], edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = elements
            .iter()
            .copied()
            .enumerate()
            .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
            .collect();
        let bonds = edges
            .iter()
            .copied()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new())
            .expect("fixed S07 CCO target topology is valid")
    }

    fn query(smarts: &str) -> QueryGraph {
        parse_smarts(smarts, &SmartsParseParams::default())
            .unwrap_or_else(|error| panic!("fixed S07 query {smarts:?} parses: {error}"))
    }

    fn graph_build_count() -> usize {
        VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get)
    }

    #[test]
    fn search_shared_perf_s07_prepared_context_reuses_six_fixed_queries() {
        VF2_GRAPH_BUILD_ENTRIES.with(|entries| entries.set(0));
        let queries = [
            FeatureQuery {
                label: "C",
                graph: query("C"),
                expected: CARBON_EXPECTED,
            },
            FeatureQuery {
                label: "O",
                graph: query("O"),
                expected: OXYGEN_EXPECTED,
            },
            FeatureQuery {
                label: "N",
                graph: query("N"),
                expected: NO_MATCHES,
            },
            FeatureQuery {
                label: "C.C",
                graph: query("C.C"),
                expected: DISCONNECTED_EXPECTED,
            },
            FeatureQuery {
                label: "C-O",
                graph: query("C-O"),
                expected: CARBON_OXYGEN_EXPECTED,
            },
            FeatureQuery {
                label: "recursive C-O root",
                graph: query("[C;$([C]-[O])]"),
                expected: RECURSIVE_EXPECTED,
            },
        ];
        let query_snapshots = queries
            .iter()
            .map(|feature| feature.graph.clone())
            .collect::<Vec<_>>();
        let topology = topology(&[Element::C, Element::C, Element::O], &[(0, 1), (1, 2)]);
        let topology_before = topology.clone();
        let coordinates = CoordinateBlock::default();
        let coordinates_before = coordinates.clone();
        let rings = fast_find_rings(&topology).expect("fixed S07 rings prepare");
        let rings_before = rings.clone();
        let valence =
            assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
                .expect("fixed S07 valence prepares");
        let valence_before = valence.clone();
        let context = build_prepared_query_match_context(&topology, &rings, &valence)
            .expect("fixed S07 prepared context validates");
        let target = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );
        let params = SubstructMatchParams::default();

        let mut actual_calls = 0;
        for sweep in 0..2 {
            for feature in &queries {
                let builds_before = graph_build_count();
                let actual = try_get_substruct_matches_with_params_and_context(
                    &target,
                    &feature.graph,
                    &params,
                    &context,
                )
                .unwrap_or_else(|error| panic!("S07 {} sweep {sweep}: {error}", feature.label));
                let builds_after = graph_build_count();
                assert_eq!(
                    builds_after, builds_before,
                    "{} sweep {sweep} rebuilt a matching graph",
                    feature.label
                );
                actual_calls += 1;

                assert_eq!(
                    actual.len(),
                    feature.expected.len(),
                    "{} sweep {sweep} result count",
                    feature.label
                );
                for (result, (expected_atoms, expected_bonds)) in
                    actual.iter().zip(feature.expected)
                {
                    assert_eq!(
                        result.atom_mapping, *expected_atoms,
                        "{} sweep {sweep} atom mapping",
                        feature.label
                    );
                    assert_eq!(
                        result.bond_mapping, *expected_bonds,
                        "{} sweep {sweep} bond mapping",
                        feature.label
                    );
                }
                assert_eq!(
                    queries
                        .iter()
                        .map(|candidate| candidate.graph.clone())
                        .collect::<Vec<_>>(),
                    query_snapshots,
                    "query inputs mutated during {} sweep {sweep}",
                    feature.label
                );
                assert_eq!(topology, topology_before, "S07 topology mutated");
                assert_eq!(coordinates, coordinates_before, "S07 coordinates mutated");
                assert_eq!(rings, rings_before, "S07 rings mutated");
                assert_eq!(valence, valence_before, "S07 valence mutated");
            }
        }

        assert_eq!(actual_calls, 12);
        assert_eq!(graph_build_count(), 0);
    }
}

#[cfg(test)]
mod search_shared_perf_s08_tests {
    use super::{
        SubstructMatchError, SubstructMatchParams, VF2_GRAPH_BUILD_ENTRIES, compile_query_graph,
        get_substruct_matches_with_compiled_query, substruct_match_impl,
        try_get_substruct_matches_with_params_and_context,
    };
    use crate::{AtomQueryPredicate, SearchTarget, build_prepared_query_match_context};
    use cosmolkit_core::{ValenceModel, assign_valence_with_options_for_topology, fast_find_rings};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, CoordinateBlock, QueryAtom, QueryGraph, QueryNode, TopologyBlock,
    };
    use cosmolkit_types::Element;
    use std::collections::BTreeMap;
    use std::sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    };

    #[test]
    fn search_shared_perf_s08_unsupported_leaf_precedes_all_private_routes() {
        let query = QueryGraph::from_parts(
            vec![QueryAtom::from_parts(
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                QueryNode::and(vec![
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
                    QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature(
                        "S08 fixed unsupported atom leaf",
                    )),
                ]),
            )],
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed S08 query with ordered unsupported leaf is valid");
        let query_before = query.clone();
        let topology = TopologyBlock::try_from_parts(
            vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("fixed S08 target topology is valid");
        let topology_before = topology.clone();
        let coordinates = CoordinateBlock::default();
        let coordinates_before = coordinates.clone();
        let rings = fast_find_rings(&topology).expect("fixed S08 rings prepare");
        let rings_before = rings.clone();
        let valence =
            assign_valence_with_options_for_topology(&topology, ValenceModel::RdkitLike, false)
                .expect("fixed S08 non-strict valence prepares");
        let valence_before = valence.clone();
        let context = build_prepared_query_match_context(&topology, &rings, &valence)
            .expect("fixed S08 prepared query context validates");
        let target = SearchTarget::new(
            &topology,
            &coordinates,
            &topology.stereo_groups,
            Some(&rings),
            Some(&valence),
        );

        let atom_callback_calls = Arc::new(AtomicUsize::new(0));
        let atom_callback_calls_in_check = Arc::clone(&atom_callback_calls);
        let final_callback_calls = Arc::new(AtomicUsize::new(0));
        let final_callback_calls_in_check = Arc::clone(&final_callback_calls);
        let params = SubstructMatchParams {
            extra_atom_check: Some(Arc::new(move |_, _, _, _| {
                atom_callback_calls_in_check.fetch_add(1, Ordering::SeqCst);
                true
            })),
            extra_atom_check_overrides_default_check: true,
            extra_final_check: Some(Arc::new(move |_, _| {
                final_callback_calls_in_check.fetch_add(1, Ordering::SeqCst);
                true
            })),
            ..SubstructMatchParams::default()
        };

        VF2_GRAPH_BUILD_ENTRIES.with(|entries| entries.set(0));
        let compiled_graph = compile_query_graph(&query);
        assert_eq!(
            VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get),
            1,
            "the compiled route's graph is built exactly once before matching"
        );

        let expected_error = SubstructMatchError::Unsupported {
            branch: "S08 fixed unsupported atom leaf",
            rdkit_function: "QueryAtom::Match",
        };
        let mut actual_calls = 0;
        for route in 0..3 {
            let builds_before = VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get);
            let result = match route {
                0 => substruct_match_impl(&target, &query, &params),
                1 => try_get_substruct_matches_with_params_and_context(
                    &target, &query, &params, &context,
                ),
                2 => get_substruct_matches_with_compiled_query(
                    &target,
                    &query,
                    &params,
                    &compiled_graph,
                ),
                _ => unreachable!("S08 freezes ordinary, prepared and compiled routes"),
            };
            let builds_after = VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get);
            assert_eq!(
                builds_after, builds_before,
                "route {route} built a matching graph after unsupported preflight"
            );
            assert_eq!(result, Err(expected_error.clone()), "route {route} error");
            assert_eq!(
                atom_callback_calls.load(Ordering::SeqCst),
                0,
                "route {route} reached an atom goal callback"
            );
            assert_eq!(
                final_callback_calls.load(Ordering::SeqCst),
                0,
                "route {route} reached a completed-goal callback"
            );
            assert_eq!(query, query_before, "route {route} mutated the query");
            assert_eq!(
                topology, topology_before,
                "route {route} mutated the topology"
            );
            assert_eq!(
                coordinates, coordinates_before,
                "route {route} mutated the coordinates"
            );
            assert_eq!(
                rings, rings_before,
                "route {route} mutated the ring assignment"
            );
            assert_eq!(
                valence, valence_before,
                "route {route} mutated the valence assignment"
            );
            actual_calls += 1;
        }

        assert_eq!(actual_calls, 3);
        assert_eq!(VF2_GRAPH_BUILD_ENTRIES.with(std::cell::Cell::get), 1);
    }
}
