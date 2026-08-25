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

use crate::search::query::{
    QueryMatchContext, and_query_match, atom_predicate_matches_with_context, atom_queries_match,
    bond_predicate_matches_with_context, bond_queries_match, build_query_match_context,
    or_query_match, xor_query_match,
};
use crate::{
    Atom, AtomQueryPredicate, Bond, BondOrder, BondQueryPredicate, BondStereo, ChiralTag, Molecule,
    StereoGroupKind,
};
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fmt;
use std::sync::Arc;

// ---------------------------------------------------------------------------
// Result types
// ---------------------------------------------------------------------------

/// Result of a single substructure match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SubstructMatchResult {
    /// Mapping from query atom index to molecule atom index.
    pub atom_mapping: Vec<usize>,
    /// Mapping from query bond index to molecule bond index.
    pub bond_mapping: Vec<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SubstructMatchError {
    #[error(
        "RDKit substructure matching branch {branch} is unsupported until {rdkit_function} is source-ported"
    )]
    Unsupported {
        branch: &'static str,
        rdkit_function: &'static str,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SubstructMatchOverload {
    Molecule,
    MolBundle,
    ResonanceMolSupplier,
    SubstructLibrary,
}

pub fn check_substruct_match_overload_support(
    overload: SubstructMatchOverload,
) -> Result<(), SubstructMatchError> {
    match overload {
        SubstructMatchOverload::Molecule => Ok(()),
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

/// Parameters controlling substructure matching behaviour.
pub type ExtraAtomCheck = Arc<dyn Fn(&Molecule, &Atom, &Molecule, &Atom) -> bool + Send + Sync>;
pub type ExtraBondCheck = Arc<dyn Fn(&Bond, &Bond) -> bool + Send + Sync>;
pub type ExtraFinalCheck = Arc<dyn Fn(&Molecule, &[usize]) -> bool + Send + Sync>;

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
        query_mol: &Molecule,
        query_atom: &Atom,
        target_mol: &Molecule,
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
        fn conformer(molecule: &Molecule, id: i32) -> Option<&crate::Conformer3D> {
            if id < 0 {
                molecule.conformers_3d().first()
            } else {
                molecule
                    .conformers_3d()
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
    fn new(query: &Molecule, recursion_possible: bool) -> Self {
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
    query: &crate::search::query::RecursiveStructureQuery,
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
        serde_json::Value::Number(value) => {
            value.as_u64().and_then(|value| usize::try_from(value).ok())
        }
        serde_json::Value::String(value) => value.parse().ok(),
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
    // or cloned. Staging prevents partial mutation on malformed input.
    if json.is_empty() {
        return Ok(());
    }
    let value: serde_json::Value = serde_json::from_str(json)?;
    let object = value
        .as_object()
        .ok_or(SubstructMatchParamsJsonError::InvalidField { field: "root" })?;
    let use_chirality = json_param_bool(object, "useChirality")?;
    let use_enhanced_stereo = json_param_bool(object, "useEnhancedStereo")?;
    let aromatic_matches_conjugated = json_param_bool(object, "aromaticMatchesConjugated")?;
    let use_query_query_matches = json_param_bool(object, "useQueryQueryMatches")?;
    let recursion_possible = json_param_bool(object, "recursionPossible")?;
    let uniquify = json_param_bool(object, "uniquify")?;
    let max_matches = json_param_usize(object, "maxMatches")?;
    let max_recursive_matches = json_param_usize(object, "maxRecursiveMatches")?;
    let num_threads = json_param_i32(object, "numThreads")?;
    let specified_stereo = json_param_bool(object, "specifiedStereoQueryMatchesUnspecified")?;
    let aromatic_matches_single_or_double =
        json_param_bool(object, "aromaticMatchesSingleOrDouble")?;

    if let Some(value) = use_chirality {
        params.use_chirality = value;
    }
    if let Some(value) = use_enhanced_stereo {
        params.use_enhanced_stereo = value;
    }
    if let Some(value) = aromatic_matches_conjugated {
        params.aromatic_matches_conjugated = value;
    }
    if let Some(value) = use_query_query_matches {
        params.use_query_query_matches = value;
    }
    if let Some(value) = recursion_possible {
        params.recursion_possible = value;
    }
    if let Some(value) = uniquify {
        params.uniquify = value;
    }
    if let Some(value) = max_matches {
        params.max_matches = value;
    }
    if let Some(value) = max_recursive_matches {
        params.max_recursive_matches = value;
    }
    if let Some(value) = num_threads {
        params.num_threads = value;
    }
    if let Some(value) = specified_stereo {
        params.specified_stereo_query_matches_unspecified = value;
    }
    if let Some(value) = aromatic_matches_single_or_double {
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
    let fields = serde_json::json!({
        "useChirality": params.use_chirality.to_string(),
        "useEnhancedStereo": params.use_enhanced_stereo.to_string(),
        "aromaticMatchesConjugated": params.aromatic_matches_conjugated.to_string(),
        "useQueryQueryMatches": params.use_query_query_matches.to_string(),
        "recursionPossible": params.recursion_possible.to_string(),
        "uniquify": params.uniquify.to_string(),
        "maxMatches": params.max_matches.to_string(),
        "maxRecursiveMatches": params.max_recursive_matches.to_string(),
        "numThreads": params.num_threads.to_string(),
        "specifiedStereoQueryMatchesUnspecified": params.specified_stereo_query_matches_unspecified.to_string(),
        "aromaticMatchesSingleOrDouble": params.aromatic_matches_single_or_double.to_string(),
    });
    serde_json::to_string_pretty(&fields).expect("fixed scalar JSON serialization cannot fail")
        + "\n"
}

// ---------------------------------------------------------------------------
// Internal minimum-degree graph representation
// ---------------------------------------------------------------------------

/// Minimal adjacency info needed for VF2.
///
/// `nbrs[i]` is a slice into `edges`.
#[derive(Debug, Clone)]
struct Vf2Graph {
    n_atoms: usize,
    n_bonds: usize,
    /// Canonical source/target endpoints indexed by bond id.
    edge_endpoints: Vec<(usize, usize)>,
    /// For each atom index, the neighbor indices and bond ids.
    adjacency: Vec<Vec<(usize, usize)>>, // (neighbor_atom_index, bond_index)
}

/// Build a VF2-compatible adjacency view from a molecule.
///
/// The C++ code iterates `out_edges` via Boost graph. We pre-build adjacency
/// once and use raw index lookups.
fn build_vf2_graph(mol: &Molecule) -> Vf2Graph {
    // RDKit source (implicit in vf2.hpp usage of out_edges):
    //   The VF2 state stores Graph *g1, *g2 and calls:
    //     boost::out_edges(node, *g)
    //     boost::out_degree(node, *g)
    //     boost::adjacent_vertices(node, *g)
    //   These are all O(1) in Boost adjacency_list.
    //
    // RDKit✔️❌: We build a flat adjacency Vec<(usize, usize)> per atom.
    //   This adds a one-time O(V+E) allocation vs the Boost inline storage,
    //   but lookups are O(degree) which matches the original hot-path cost.
    let n_atoms = mol.num_atoms();
    let mut adjacency: Vec<Vec<(usize, usize)>> = vec![Vec::new(); n_atoms];
    let mut edge_endpoints = Vec::with_capacity(mol.num_bonds());
    for (bond_idx, bond) in mol.bonds().iter().enumerate() {
        let b = bond.begin().index();
        let e = bond.end().index();
        edge_endpoints.push((b, e));
        adjacency[b].push((e, bond_idx));
        adjacency[e].push((b, bond_idx));
    }
    Vf2Graph {
        n_atoms,
        n_bonds: mol.num_bonds(),
        edge_endpoints,
        adjacency,
    }
}

fn get_other_idx(g: &Vf2Graph, edge: usize, vertex: NodeId) -> NodeId {
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
    let (source, target) = g.edge_endpoints[edge];
    if source == vertex { target } else { source }
}

impl Vf2Graph {
    fn out_degree(&self, node: usize) -> usize {
        self.adjacency[node].len()
    }

    fn out_edges(&self, node: usize) -> &[(usize, usize)] {
        &self.adjacency[node]
    }
}

// ---------------------------------------------------------------------------
// Atom and bond matching functors
// ---------------------------------------------------------------------------

fn property_compat(
    properties1: &BTreeMap<String, String>,
    properties2: &BTreeMap<String, String>,
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

fn has_chiral_label(atom: &Atom) -> bool {
    // RDKit✔️✔️: bool hasChiralLabel(const Atom *at) {
    // RDKit✔️✔️:   PRECONDITION(at, "bad atom");
    // RDKit✔️✔️:   return at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CW ||
    // RDKit✔️✔️:          at->getChiralTag() == Atom::CHI_TETRAHEDRAL_CCW;
    // RDKit✔️✔️: }
    // Rust's reference type enforces the non-null precondition. Complexity
    // review: both implementations read one enum and perform at most two O(1)
    // comparisons without allocation.
    matches!(
        atom.chiral_tag(),
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
    query: &Molecule,
    mol: &Molecule,
    query_index: usize,
    mol_index: usize,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
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
        && has_chiral_label(query_atom)
        && !params.specified_stereo_query_matches_unspecified
        && !has_chiral_label(mol_atom)
    {
        return false;
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

fn atom_matches(query_atom: &Atom, query_mol: &Molecule, mol_atom: &Atom, mol: &Molecule) -> bool {
    if let Some(query_node) = query_atom.query() {
        let query_ctx = build_query_match_context(mol);
        return evaluate_atom_query(
            query_node,
            mol_atom,
            mol,
            &SubstructMatchParams::default(),
            None,
            &query_ctx,
        );
    }

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
    // Local complexity review: the plain-atom path is constant time and uses
    // only scalar field reads, exactly as the source. No allocation, cloning,
    // molecule scan, keyed lookup, or temporary collection is introduced.
    if query_atom.atomic_number() != mol_atom.atomic_number() {
        return false;
    }
    let _ = query_mol;
    if query_atom.atomic_number() == 0 {
        return match (query_atom.isotope(), mol_atom.isotope()) {
            (Some(query_isotope), Some(mol_isotope)) => query_isotope == mol_isotope,
            _ => true,
        };
    }
    (query_atom.formal_charge() == 0 || query_atom.formal_charge() == mol_atom.formal_charge())
        && (query_atom.isotope().is_none() || query_atom.isotope() == mol_atom.isotope())
        && (query_atom.radical_electrons() == 0
            || query_atom.radical_electrons() == mol_atom.radical_electrons())
}

fn recursive_smarts_root_matches(
    atom: &Atom,
    recursive_query: &crate::search::query::RecursiveStructureQuery,
    mol: &Molecule,
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

    let Some(query) = recursive_query.query_mol() else {
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
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
    match pred {
        // RDKit✔️✔️: Chiral SMARTS labels are not ordinary atom-compatibility
        // constraints when `useChirality` is false. AtomLabelFunctor and
        // MolMatchFinalCheckFunctor handle stereochemistry explicitly.
        AtomQueryPredicate::ChiralTagMatch(_) | AtomQueryPredicate::ChiralPermutationMatch(_)
            if !params.use_chirality =>
        {
            true
        }
        AtomQueryPredicate::RecursiveSmarts(recursive_query) => {
            recursive_smarts_root_matches(atom, recursive_query, mol, recursive_cache)
        }
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
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
    match query {
        crate::QueryNode::Predicate(pred) => atom_query_predicate_matches_for_substruct(
            atom,
            pred,
            mol,
            params,
            recursive_cache,
            query_ctx,
        ),
        crate::QueryNode::And(children) => and_query_match(children, false, |child| {
            evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)
        }),
        crate::QueryNode::Or(children) => or_query_match(children, false, |child| {
            evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)
        }),
        crate::QueryNode::Xor(children) => xor_query_match(children, false, |child| {
            evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)
        }),
        crate::QueryNode::Not(child) => {
            !evaluate_atom_query(child, atom, mol, params, recursive_cache, query_ctx)
        }
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
    query: &Molecule,
    mol: &Molecule,
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
        && query_bond.order() == BondOrder::Double
        && rdkit_bond_stereo_is_above_any(query_bond.stereo())
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
    mol: &Molecule,
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
//    and vice versa. Uses `Option<usize>` instead of NULL_NODE sentinel.
//
// 2. `term_1`/`term_2`: Stores the core_len *depth* at which each atom was
//    added to the terminal set, exactly as in vf2.hpp. BackTrack decrements
//    counters keyed by depth, not recomputes from scratch.
//
// 3. No shared_ptr copy semantics: VF2SubState in RDKit uses COW with
//    `share_count`. Rust's Clone+Vf2State avoids raw pointer sharing.
//    This means each VF2 recursive branch owns its state, which is
//    semantically correct but allocates O(depth * n) instead of
//    O(n) shared storage. For typical molecule sizes (<1000 atoms) this
//    is negligible; for very large searches the COW approach could be
//    reinstated with Arc<Vec<NodeId>>.
//
// 4. No boost graph: hand-rolled Vf2Graph adjacency.

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
    id: usize,
    in_deg: usize,
    out_deg: usize,
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

fn sort_nodes_by_frequency(g: &Vf2Graph) -> Vec<NodeId> {
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
    let mut vect: Vec<NodeInfo> = (0..g.n_atoms)
        .map(|i| {
            let deg = g.out_degree(i);
            NodeInfo {
                id: i,
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
            vect[i + j].in_deg += vect[i + j].out_deg; // valence sum
            vect[i + j].out_deg = run; // frequency
        }
        i += run;
    }

    // RDKit✔️✔️:   std::sort(vect.begin(), vect.end(), nodeInfoComp2);
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
    vect.iter().map(|ni| ni.id).collect()
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
/// core_1[i] = mapping from query atom i -> mol atom j (or None).
/// core_2[j] = mapping from mol atom j -> query atom i (or None).
/// term_1[i] = depth (core_len) when atom i entered terminal set (0 = not terminal).
/// term_2[j] = same for mol atoms.
struct Vf2SubState<'a> {
    g1: &'a Vf2Graph,
    g2: &'a Vf2Graph,
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
    fn new(g1: &'a Vf2Graph, g2: &'a Vf2Graph, sort_nodes: bool) -> Self {
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
        let n1 = g1.n_atoms;
        let n2 = g2.n_atoms;
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
                for &(query_neighbor, _) in self.g1.out_edges(pair.n1) {
                    if self.core_1[query_neighbor] != NULL_NODE {
                        mapped_terminal_neighbor = self.core_1[query_neighbor];
                        break;
                    }
                }
                debug_assert_ne!(mapped_terminal_neighbor, NULL_NODE);
                if mapped_terminal_neighbor != NULL_NODE {
                    pair.nbr_node = mapped_terminal_neighbor;
                    pair.nbr_cursor = 0;
                    pair.nbr_end = self.g2.out_edges(mapped_terminal_neighbor).len();
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
            // RDKit✔️✔️: while (pair.nbrbeg < pair.nbrend && core_2[*pair.nbrbeg] != NULL_NODE) {
            // RDKit✔️✔️:   ++pair.nbrbeg;
            // RDKit✔️✔️: }
            let neighbors = self.g2.out_edges(pair.nbr_node);
            while pair.nbr_cursor < pair.nbr_end
                && self.core_2[neighbors[pair.nbr_cursor].0] != NULL_NODE
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
                pair.n2 = neighbors[pair.nbr_cursor].0;
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
        for &(_, edge_idx1) in self.g1.out_edges(node1) {
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
        for &(nbr, bond_idx) in self.g2.out_edges(a) {
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
        for &(_, edge) in self.g1.out_edges(node1) {
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
        for &(_, edge) in self.g2.out_edges(node2) {
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
        for &(_, edge) in self.g1.out_edges(node1) {
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
        for &(_, edge) in self.g2.out_edges(node2) {
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

    fn get_core_set(&self) -> (Vec<NodeId>, Vec<NodeId>) {
        // RDKit✔️❌: void GetCoreSet(node_id c1[], node_id c2[]) {
        // RDKit✔️❌:   unsigned int i, j;
        // RDKit✔️❌:   for (i = 0, j = 0; i < n1; ++i) {
        // RDKit✔️❌:     if (core_1[i] != NULL_NODE) {
        // RDKit✔️❌:       c1[j] = i;
        // RDKit✔️❌:       c2[j] = core_1[i];
        // RDKit✔️❌:       ++j;
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌: }
        // Complexity review: both scan n1 entries in O(V) and write core_len
        // outputs. Rust allocates two result Vecs here, whereas RDKit writes
        // into caller-provided arrays, so repeated goal checks pay two extra
        // allocations despite identical mapping order and asymptotic cost.
        let mut c1 = Vec::with_capacity(self.core_len);
        let mut c2 = Vec::with_capacity(self.core_len);
        for i in 0..self.n1 {
            if self.core_1[i] != NULL_NODE {
                c1.push(i);
                c2.push(self.core_1[i]);
            }
        }
        (c1, c2)
    }

    fn match_one(
        &mut self,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
        mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    ) -> Option<(Vec<NodeId>, Vec<NodeId>)> {
        // RDKit✔️❌: bool Match(node_id c1[], node_id c2[]) {
        // RDKit✔️❌:   if (IsGoal()) {
        // RDKit✔️❌:     GetCoreSet(c1, c2);
        // RDKit✔️❌:     if (MatchChecks(c1, c2)) {
        // RDKit✔️❌:       return true;
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   if (IsDead()) {
        // RDKit✔️❌:     return false;
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   Pair<Graph> pair;
        // RDKit✔️❌:   while (NextPair(pair)) {
        // RDKit✔️❌:     if (IsFeasiblePair(pair.n1, pair.n2)) {
        // RDKit✔️❌:       AddPair(pair.n1, pair.n2);
        // RDKit✔️❌:       if (Match(c1, c2)) {  // recurse
        // RDKit✔️❌:         return true;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       BackTrack(pair.n1, pair.n2);
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:   return false;
        // RDKit✔️❌: }
        // Complexity review: candidate generation, feasibility checks, and
        // depth-first recursion match RDKit's search tree. The known gap is
        // inherited from get_core_set(), which allocates two Vecs at each goal.
        if self.is_goal() {
            let (c1, c2) = self.get_core_set();
            let accepted = match match_check.as_mut() {
                Some(check) => self.match_checks(&c1, &c2, check),
                None => true,
            };
            if accepted {
                return Some((c1, c2));
            }
        }
        if self.is_dead() {
            return None;
        }
        let mut pair = Vf2Pair::new();
        while self.next_pair(&mut pair) {
            if self.is_feasible_pair(pair.n1, pair.n2, atom_fn, bond_fn) {
                self.add_pair(pair.n1, pair.n2);
                if let Some(result) = self.match_one(atom_fn, bond_fn, match_check.as_deref_mut()) {
                    return Some(result);
                }
                self.back_track(pair.n1, pair.n2);
            }
        }
        None
    }

    fn match_all(
        &mut self,
        atom_fn: &impl Fn(usize, usize) -> bool,
        bond_fn: &impl Fn(usize, usize) -> bool,
        mut match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
        results: &mut Vec<(Vec<NodeId>, Vec<NodeId>)>,
        max_matches: usize,
    ) -> bool {
        // RDKit✔️❌: template <class DoubleBackInsertionSequence>
        // RDKit✔️❌: bool MatchAll(node_id c1[], node_id c2[], DoubleBackInsertionSequence &res,
        // RDKit✔️❌:               unsigned int lim = 0) {
        // RDKit✔️❌:   if (IsGoal()) {
        // RDKit✔️❌:     GetCoreSet(c1, c2);
        // RDKit✔️❌:     if (MatchChecks(c1, c2)) {
        // RDKit✔️❌:       typename DoubleBackInsertionSequence::value_type newSeq;
        // RDKit✔️❌:       newSeq.reserve(core_len);
        // RDKit✔️❌:       for (unsigned int i = 0; i < core_len; ++i) {
        // RDKit✔️❌:         newSeq.emplace_back(c1[i], c2[i]);
        // RDKit✔️❌:       }
        // RDKit✔️❌:       res.push_back(newSeq);
        // RDKit✔️❌:       return lim && res.size() >= lim;
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   if (IsDead()) {
        // RDKit✔️❌:     return false;
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   Pair<Graph> pair;
        // RDKit✔️❌:   while (NextPair(pair)) {
        // RDKit✔️❌:     if (IsFeasiblePair(pair.n1, pair.n2)) {
        // RDKit✔️❌:       AddPair(pair.n1, pair.n2);
        // RDKit✔️❌:       if (MatchAll(c1, c2, res, lim)) {  // recurse
        // RDKit✔️❌:         return true;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       BackTrack(pair.n1, pair.n2);
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:   return false;
        // RDKit✔️❌: }
        // Complexity review: the DFS search tree, candidate order, early limit,
        // and per-result O(core_len) storage match RDKit. The known extra cost
        // is get_core_set() allocating two Vecs before each final check.
        if self.is_goal() {
            let (c1, c2) = self.get_core_set();
            let accepted = match match_check.as_mut() {
                Some(check) => self.match_checks(&c1, &c2, check),
                None => true,
            };
            if accepted {
                results.push((c1, c2));
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

/// RDKit✔️❌: Match — find first match via VF2 recursion.
///
/// Matches RDKit's `Match(c1, c2)` entry point. `match_check` allows
/// final verification (like MolMatchFinalCheckFunctor). If None, all
/// completed matches are accepted.
fn vf2_match(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
) -> Option<(Vec<NodeId>, Vec<NodeId>)> {
    // RDKit✔️❌: template <class SubState>
    // RDKit✔️❌: bool match(int *pn, node_id c1[], node_id c2[], SubState &s) {
    // RDKit✔️❌:   if (s.Match(c1, c2)) {
    // RDKit✔️❌:     // not needed, pn = num query atoms (n1)...
    // RDKit✔️❌:     *pn = s.CoreLen();
    // RDKit✔️❌:     return true;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return false;
    // RDKit✔️❌: }
    // Rust returns the mapping and its length is available directly. Complexity
    // and allocation behavior are exactly those of the single member core.
    state.match_one(atom_fn, bond_fn, match_check)
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

/// RDKit✔️❌: MatchAll — find all matches up to `max_matches`.
///
/// Collects matches into `results` as (c1, c2) pairs.
/// Returns true when the limit has been reached, signaling the caller
/// to stop.
fn vf2_match_all(
    state: &mut Vf2SubState,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    results: &mut Vec<(Vec<NodeId>, Vec<NodeId>)>,
    max_matches: usize,
) -> bool {
    // RDKit✔️❌: template <class SubState, class DoubleBackInsertionSequence>
    // RDKit✔️❌: bool match(node_id c1[], node_id c2[], SubState &s,
    // RDKit✔️❌:            DoubleBackInsertionSequence &res, unsigned int max_results) {
    // RDKit✔️❌:   s.MatchAll(c1, c2, res, max_results);
    // RDKit✔️❌:   return !res.empty();
    // RDKit✔️❌: }
    // Complexity review: this wrapper adds one emptiness check after invoking
    // the single member recursion core; it does not copy or re-enumerate results.
    state.match_all(atom_fn, bond_fn, match_check, results, max_matches);
    !results.is_empty()
}

fn vf2_entry_one(
    g1: &Vf2Graph,
    g2: &Vf2Graph,
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
    // Complexity review: both allocate two O(V) mapping buffers, construct one
    // unsorted state, run the same first-match DFS, and fill one O(V) result.
    let mut state = Vf2SubState::new(g1, g2, false);
    result.clear();
    if let Some((c1, c2)) = vf2_match(&mut state, atom_fn, bond_fn, match_check) {
        result.reserve(c1.len());
        result.extend(c1.into_iter().zip(c2));
    }
    !result.is_empty()
}

fn vf2_entry_all(
    g1: &Vf2Graph,
    g2: &Vf2Graph,
    atom_fn: &impl Fn(usize, usize) -> bool,
    bond_fn: &impl Fn(usize, usize) -> bool,
    match_check: Option<&mut impl FnMut(&[NodeId], &[NodeId]) -> bool>,
    results: &mut Vec<(Vec<NodeId>, Vec<NodeId>)>,
    max_results: usize,
) -> bool {
    // RDKit✔️❌: template <class Graph, class VertexLabeling  // binary predicate
    // RDKit✔️❌:           ,
    // RDKit✔️❌:           class EdgeLabeling  // binary predicate
    // RDKit✔️❌:           ,
    // RDKit✔️❌:           class MatchChecking  // binary predicate
    // RDKit✔️❌:           ,
    // RDKit✔️❌:           class DoubleBackInsertionSequence  // contains a back insertion
    // RDKit✔️❌:                                              // sequence
    // RDKit✔️❌:           >
    // RDKit✔️❌: bool vf2_all(const Graph &g1, const Graph &g2, VertexLabeling &vertex_labeling,
    // RDKit✔️❌:              EdgeLabeling &edge_labeling, MatchChecking &match_checking,
    // RDKit✔️❌:              DoubleBackInsertionSequence &F, unsigned int max_results = 1000) {
    // RDKit✔️❌:   detail::VF2SubState<const Graph, VertexLabeling, EdgeLabeling, MatchChecking>
    // RDKit✔️❌:       s0(&g1, &g2, vertex_labeling, edge_labeling, match_checking, false);
    // RDKit✔️❌:   std::unique_ptr<detail::node_id[]> ni1(new detail::node_id[num_vertices(g1)]);
    // RDKit✔️❌:   std::unique_ptr<detail::node_id[]> ni2(new detail::node_id[num_vertices(g2)]);
    // RDKit✔️❌:
    // RDKit✔️❌:   F.clear();
    // RDKit✔️❌:   F.resize(0);
    // RDKit✔️❌:
    // RDKit✔️❌:   match(ni1.get(), ni2.get(), s0, F, max_results);
    // RDKit✔️❌:
    // RDKit✔️❌:   return !F.empty();
    // RDKit✔️❌: };
    // Complexity review: search order and result storage match RDKit. The known
    // gap is the member core allocating mapping Vecs before each final check,
    // while RDKit reuses ni1/ni2 across goal states.
    let mut state = Vf2SubState::new(g1, g2, false);
    results.clear();
    vf2_match_all(
        &mut state,
        atom_fn,
        bond_fn,
        match_check,
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
    crate::source_port_helpers::count_swaps_to_interconvert(reference, probe)
        .ok()
        .and_then(|swaps| u32::try_from(swaps).ok())
}

fn rdkit_atom_perturbation_order_from_bond_indices(
    mol: &Molecule,
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
        .topology_block()
        .adjacency
        .neighbors_of(atom_idx)
        .iter()
        .map(|neighbor| i32::try_from(neighbor.bond.index()))
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
    mol: &Molecule,
    query: &Molecule,
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
    fn new(_query: &Molecule, mol: &Molecule, params: &SubstructMatchParams) -> Self {
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

fn find_bond_between(mol: &Molecule, begin: usize, end: usize) -> Option<&Bond> {
    mol.bonds().iter().find(|bond| {
        let b = bond.begin().index();
        let e = bond.end().index();
        (b == begin && e == end) || (b == end && e == begin)
    })
}

fn rdkit_match_final_check(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    c1: &[NodeId],
    c2: &[NodeId],
    setup: &MolMatchFinalCheckSetup,
    matches_seen: &mut Vec<Vec<bool>>,
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
    if params.use_generic_matchers && !super::generic_groups::generic_atom_matcher(mol, query, c2) {
        return Ok(false);
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
        if matches_seen.iter().any(|existing| *existing == mask) {
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
            matches_seen.push(mask);
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
        if query.topology_block().adjacency.neighbors_of(qi).len() < 3 || !has_chiral_label(q_at) {
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
        if !has_chiral_label(m_at) {
            if params.specified_stereo_query_matches_unspecified {
                continue;
            }
            return Ok(false);
        }
        // RDKit✔️✔️:     if (qAt->getDegree() > mAt->getDegree()) {
        // RDKit✔️✔️:       return false;
        // RDKit✔️✔️:     }
        if query.topology_block().adjacency.neighbors_of(qi).len()
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
            let Some(q_bond) = find_bond_between(query, qi, qj) else {
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
            if m_order.len() == query.topology_block().adjacency.neighbors_of(qi).len() {
                break;
            }
        }
        if q_order.len() != query.topology_block().adjacency.neighbors_of(qi).len()
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
        if q_bnd.order() != BondOrder::Double || !rdkit_bond_stereo_is_above_any(q_bnd.stereo()) {
            continue;
        }
        // RDKit✔️✔️:     if (qBnd->getStereoAtoms().size() != 2) {
        // RDKit✔️✔️:       continue;
        // RDKit✔️✔️:     }
        let Some(q_stereo_atoms) = q_bnd.stereo_atoms() else {
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
        let q_stereo = rdkit_translate_ez_label_to_cis_trans(q_bnd.stereo());
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
        matches_seen.push(mask);
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
            let inner_query = query.query_mol().ok_or(SubstructMatchError::Unsupported {
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

fn preflight_query_molecule(query: &Molecule) -> Result<(), SubstructMatchError> {
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
        if let Some(query) = atom.query() {
            preflight_atom_query(query)?;
        }
    }
    for bond in query.bonds() {
        if let Some(query) = bond.query() {
            preflight_bond_query(query)?;
        }
    }
    Ok(())
}

fn recursive_matcher(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
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
    // ordered set, while the canonical VF2 mapping-allocation gap remains.
    let mut local_params = params.clone();
    local_params.max_matches = params.max_recursive_matches.max(params.max_matches);
    local_params.uniquify = false;
    for atom in query.atoms() {
        if let Some(query_node) = atom.query() {
            match_subqueries(mol, query_node, &local_params, recursive_cache)?;
        }
    }

    let matches = substruct_match_impl_with_recursive_cache(
        mol,
        query,
        &local_params,
        Some(recursive_cache),
    )?;
    let root_index = query
        .prop("_queryRootAtom")
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(0);
    let mut match_starts = vec![false; mol.num_atoms()];
    for matched in matches.into_iter().take(local_params.max_matches) {
        if let Some(&root_atom_idx) = matched.atom_mapping.get(root_index)
            && root_atom_idx != NULL_NODE
            && root_atom_idx < match_starts.len()
        {
            match_starts[root_atom_idx] = true;
        }
    }
    Ok(match_starts)
}

fn match_subqueries(
    mol: &Molecule,
    query: &crate::QueryNode<AtomQueryPredicate>,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
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
                let match_starts = match recursive_query.query_mol() {
                    Some(inner_query) => {
                        recursive_matcher(mol, inner_query, params, recursive_cache)?
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
                match_subqueries(mol, child, params, recursive_cache)?;
            }
        }
        crate::QueryNode::Not(child) => {
            match_subqueries(mol, child, params, recursive_cache)?;
        }
    }
    Ok(())
}

fn populate_recursive_query_match_cache(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: &mut RecursiveQueryMatchCache,
) -> Result<(), SubstructMatchError> {
    for atom in query.atoms() {
        if let Some(query_node) = atom.query() {
            match_subqueries(mol, query_node, params, recursive_cache)?;
        }
    }
    Ok(())
}

fn substruct_match_impl_with_recursive_cache(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
) -> SubstructMatchResultList {
    let m_num_atoms = mol.num_atoms();
    let q_num_atoms = query.num_atoms();
    let query_ctx = build_query_match_context(mol);

    // RDKit source (SubstructMatch.cpp):
    //   if (!mNumAtoms || !qNumAtoms || qNumAtoms > mNumAtoms) {
    //     return matches;
    //   }
    if m_num_atoms == 0 || q_num_atoms == 0 || q_num_atoms > m_num_atoms {
        return Ok(Vec::new());
    }

    // Build VF2 graphs.
    let q_graph = build_vf2_graph(query);
    let m_graph = build_vf2_graph(mol);

    // Build atom matching closure.
    // RDKit source:
    //   detail::AtomLabelFunctor atomLabeler(query, mol, params);
    //   detail::BondLabelFunctor bondLabeler(query, mol, params);
    //   MolMatchFinalCheckFunctor matchChecker(query, mol, params);
    let atom_fn = |qi: usize, mj: usize| -> bool {
        atom_label_matches(query, mol, qi, mj, params, recursive_cache, &query_ctx)
    };

    let bond_fn = |qei: usize, mei: usize| -> bool {
        bond_label_matches(query, mol, qei, mei, params, &query_ctx)
    };

    // RDKit source:
    //   bool found = boost::vf2_all(query.getTopology(), mol.getTopology(),
    //                               atomLabeler, bondLabeler, matchChecker,
    //                               pms, params.maxMatches);
    let mut raw_matches: Vec<(Vec<NodeId>, Vec<NodeId>)> = Vec::new();
    let mut matches_seen: Vec<Vec<bool>> = Vec::new();
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

    vf2_entry_all(
        &q_graph,
        &m_graph,
        &atom_fn,
        &bond_fn,
        Some(&mut check_fn),
        &mut raw_matches,
        params.max_matches,
    );
    if let Some(err) = final_check_error {
        return Err(err);
    }

    // RDKit source (SubstructMatch.cpp):
    //   if (found) {
    //     const unsigned int nQueryAtoms = query.getNumAtoms();
    //     matches.reserve(pms.size());
    //     MatchVectType matchVect(nQueryAtoms);
    //     for (const auto &pairs : pms) {
    //       for (const auto &pair : pairs) {
    //         matchVect[pair.first] = pair;
    //       }
    //       matches.push_back(matchVect);
    //     }
    //   }
    let mut results: Vec<SubstructMatchResult> = Vec::new();

    for (c1, c2) in &raw_matches {
        // Build atom_mapping: query_atom_index -> mol_atom_index.
        // RDKit uses MatchVectType (vector<pair<int,int>>) where
        // pair.second is the mol atom index and pair.first is query atom index.
        let mut atom_to_mol: Vec<Option<usize>> = vec![None; q_num_atoms];
        for (&qa, &ma) in c1.iter().zip(c2.iter()) {
            if qa < q_num_atoms {
                atom_to_mol[qa] = Some(ma);
            }
        }

        // Build bond mapping by looking up bonds between matched atoms.
        let mut bond_mapping = Vec::with_capacity(query.num_bonds());
        for qbond in query.bonds() {
            let q_begin = qbond.begin().index();
            let q_end = qbond.end().index();
            let m_begin = atom_to_mol[q_begin];
            let m_end = atom_to_mol[q_end];
            match (m_begin, m_end) {
                (Some(mb), Some(me)) => {
                    // Find bond between mb and me in mol.
                    let found = m_graph.adjacency[mb]
                        .iter()
                        .find(|&&(nbr, _)| nbr == me)
                        .map(|&(_, eidx)| eidx);
                    bond_mapping.push(found.unwrap_or(NULL_NODE));
                }
                _ => {
                    bond_mapping.push(NULL_NODE);
                }
            }
        }

        results.push(SubstructMatchResult {
            atom_mapping: atom_to_mol
                .into_iter()
                .map(|x| x.unwrap_or(NULL_NODE))
                .collect(),
            bond_mapping,
        });
    }

    Ok(results)
}

fn atom_compat(
    query_atom: &Atom,
    query_mol: &Molecule,
    mol_atom: &Atom,
    mol: &Molecule,
    params: &SubstructMatchParams,
    recursive_cache: Option<&RecursiveQueryMatchCache>,
    query_ctx: &QueryMatchContext,
) -> bool {
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
    // the same constant-time option/flag dispatch, one default query or atom
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
        return extra_atom_check(query_mol, query_atom, mol, mol_atom);
    }

    let matches = if params.use_query_query_matches
        && let (Some(query), Some(mol_query)) = (query_atom.query(), mol_atom.query())
    {
        atom_queries_match(query, mol_query)
    } else if let Some(query_node) = query_atom.query() {
        evaluate_atom_query(
            query_node,
            mol_atom,
            mol,
            params,
            recursive_cache,
            query_ctx,
        )
    } else {
        atom_matches(query_atom, query_mol, mol_atom, mol)
    };
    if !matches {
        return false;
    }
    if !params.atom_properties.is_empty()
        && !property_compat(
            query_atom.props(),
            mol_atom.props(),
            &params.atom_properties,
        )
    {
        return false;
    }
    if let Some(extra_atom_check) = &params.extra_atom_check
        && !extra_atom_check(query_mol, query_atom, mol, mol_atom)
    {
        return false;
    }
    matches
}

#[allow(deprecated)]
fn chiral_atom_compat(
    query_atom: &Atom,
    query_mol: &Molecule,
    mol_atom: &Atom,
    mol: &Molecule,
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
    let mut matches = atom_matches(query_atom, query_mol, mol_atom, mol);
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
    query_bond: &Bond,
    query_mol: &Molecule,
    mol_bond: &Bond,
    mol: &Molecule,
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
        return extra_bond_check(query_bond, mol_bond);
    }

    let is_conjugated_single_or_double = |bond: &Bond| {
        bond.is_conjugated() && matches!(bond.order(), BondOrder::Single | BondOrder::Double)
    };
    let is_single_or_double =
        |bond: &Bond| matches!(bond.order(), BondOrder::Single | BondOrder::Double);
    let aromatic_pair_matches = |other_matches: &dyn Fn(&Bond) -> bool| {
        (query_bond.order() == BondOrder::Aromatic && mol_bond.order() == BondOrder::Aromatic)
            || (query_bond.order() == BondOrder::Aromatic && other_matches(mol_bond))
            || (mol_bond.order() == BondOrder::Aromatic && other_matches(query_bond))
    };

    let matches = if params.use_query_query_matches
        && let (Some(query), Some(mol_query)) = (query_bond.query(), mol_bond.query())
    {
        bond_queries_match(query, mol_query)
    } else if params.aromatic_matches_conjugated
        && query_bond.query().is_none()
        && mol_bond.query().is_none()
        && aromatic_pair_matches(&is_conjugated_single_or_double)
    {
        true
    } else if params.aromatic_matches_single_or_double
        && query_bond.query().is_none()
        && mol_bond.query().is_none()
        && aromatic_pair_matches(&is_single_or_double)
    {
        true
    } else if let Some(query) = query_bond.query() {
        evaluate_bond_query(query, mol_bond, mol, query_ctx)
    } else {
        query_bond.order() == BondOrder::Unspecified
            || mol_bond.order() == BondOrder::Unspecified
            || query_bond.order() == mol_bond.order()
    };
    if !matches {
        return false;
    }

    if query_bond.order() == BondOrder::Dative && mol_bond.order() == BondOrder::Dative {
        let query_begin = &query_mol.atoms()[query_bond.begin().index()];
        let query_end = &query_mol.atoms()[query_bond.end().index()];
        let mol_begin = &mol.atoms()[mol_bond.begin().index()];
        let mol_end = &mol.atoms()[mol_bond.end().index()];
        if !atom_matches(query_begin, query_mol, mol_begin, mol)
            || !atom_matches(query_end, query_mol, mol_end, mol)
        {
            return false;
        }
    }
    if !params.bond_properties.is_empty()
        && !property_compat(
            query_bond.props(),
            mol_bond.props(),
            &params.bond_properties,
        )
    {
        return false;
    }
    if let Some(extra_bond_check) = &params.extra_bond_check
        && !extra_bond_check(query_bond, mol_bond)
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
    molecule: &Molecule,
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
    let is_dummy = match atom.query() {
        None => atom.atomic_number() == 0,
        Some(crate::QueryNode::Predicate(AtomQueryPredicate::Any)) => true,
        Some(_) => false,
    };
    (molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom_index)
        .len()
        == 1
        && is_dummy)
        || atom
            .query()
            .is_some_and(|query| query_contains_atomic_number(query, 1))
}

fn core_substitution_score(
    molecule: &Molecule,
    query: &Molecule,
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
    molecule: &Molecule,
    query: &Molecule,
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

fn sort_matches_by_degree_of_core_substitution(
    molecule: &Molecule,
    query: &Molecule,
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
    mol: &Molecule,
    query: &Molecule,
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
    // Complexity review: preflight adds one linear query-tree scan for
    // fail-closed unsupported leaves. Recursive preparation, VF2 search, and
    // result materialization otherwise retain RDKit's asymptotic behavior.
    // The second marker remains ❌ because Rust's VF2 result path allocates
    // mapping Vecs at goal checks, as documented on the canonical VF2 core.
    preflight_query_molecule(query)?;
    if mol.num_atoms() == 0 || query.num_atoms() == 0 || query.num_atoms() > mol.num_atoms() {
        return Ok(Vec::new());
    }
    let mut recursive_locker = RecursiveLocker::new(query, params.recursion_possible);
    if params.recursion_possible {
        populate_recursive_query_match_cache(mol, query, params, &mut recursive_locker.cache)?;
    }
    substruct_match_impl_with_recursive_cache(mol, query, params, Some(&recursive_locker.cache))
}

/// Check if a molecule contains a substructure match for the given query.
///
/// This is the public API for `has_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn has_substruct_match(mol: &Molecule, query: &Molecule) -> bool {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    substruct_match_impl(mol, query, &params)
        .map(|matches| !matches.is_empty())
        .unwrap_or(false)
}

/// Get the first substructure match, if any.
///
/// This is the public API for `get_substruct_match`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_match(mol: &Molecule, query: &Molecule) -> Option<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    let mut params = params;
    params.max_matches = 1;
    substruct_match_impl(mol, query, &params)
        .ok()
        .and_then(|matches| matches.into_iter().next())
}

/// Get all substructure matches with default parameters.
///
/// This is the public API for `get_substruct_matches`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches(mol: &Molecule, query: &Molecule) -> Vec<SubstructMatchResult> {
    let params = SubstructMatchParams::default();
    substruct_match_impl(mol, query, &params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters.
///
/// This is the public API for `get_substruct_matches_with_params`.
/// RDKit✔️❌: VF2-based substructure matching ported from vf2.hpp + SubstructMatch.cpp.
pub fn get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> Vec<SubstructMatchResult> {
    substruct_match_impl(mol, query, params).unwrap_or_default()
}

/// Get all substructure matches with custom parameters and structured
/// unsupported-feature errors for source-porting callers.
pub fn try_get_substruct_matches_with_params(
    mol: &Molecule,
    query: &Molecule,
    params: &SubstructMatchParams,
) -> SubstructMatchResultList {
    substruct_match_impl(mol, query, params)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::MoleculeBuilder;
    use crate::search::smarts_parse::compile_query_fixture;

    #[test]
    fn smarts_ring_connectivity_zero_rejects_ring_atoms() {
        let chain = Molecule::from_smiles("CCC").expect("chain");
        let ring = Molecule::from_smiles("C1CCCCC1").expect("ring");
        let no_ring_bonds = compile_query_fixture("[Cx0]").expect("x0 query");
        let has_ring_bond = compile_query_fixture("[Cx]").expect("x query");

        assert!(!get_substruct_matches(&chain, &no_ring_bonds).is_empty());
        assert!(get_substruct_matches(&ring, &no_ring_bonds).is_empty());
        assert!(get_substruct_matches(&chain, &has_ring_bond).is_empty());
        assert!(!get_substruct_matches(&ring, &has_ring_bond).is_empty());
    }

    #[test]
    fn shared_count_swaps_substruct_preserves_none_failure_mapping() {
        assert_eq!(count_swaps_to_interconvert_i32(&[1, 2], &[1]), None);
        assert_eq!(count_swaps_to_interconvert_i32(&[1, 2], &[1, 3]), None);
    }

    fn make_mol_c() -> Molecule {
        // Methane: C
        let mut builder = MoleculeBuilder::new();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder.build().expect("build methane")
    }

    fn make_mol_cc() -> Molecule {
        // Ethane: CC
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        builder.build().expect("build ethane")
    }

    fn make_mol_cco() -> Molecule {
        // Ethanol: CCO
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add C-C bond");
        builder
            .add_bond(crate::BondSpec::new(c1, o, BondOrder::Single))
            .expect("add C-O bond");
        builder.build().expect("build ethanol")
    }

    fn make_mol_coc() -> Molecule {
        // Dimethyl ether: COC
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, o, BondOrder::Single))
            .expect("add C-O bond 1");
        builder
            .add_bond(crate::BondSpec::new(o, c1, BondOrder::Single))
            .expect("add O-C bond 2");
        builder.build().expect("build dimethyl ether")
    }

    #[test]
    fn test_has_substruct_match_self() {
        let c = make_mol_c();
        assert!(
            has_substruct_match(&c, &c),
            "a molecule should match itself"
        );
    }

    #[test]
    fn test_has_substruct_match_cc_in_cco() {
        let cc = make_mol_cc();
        let cco = make_mol_cco();
        assert!(
            has_substruct_match(&cco, &cc),
            "CCO should contain CC as substructure"
        );
    }

    #[test]
    fn test_has_substruct_match_no_match() {
        let c = make_mol_c();
        let cco = make_mol_cco();
        assert!(
            !has_substruct_match(&c, &cco),
            "a single carbon should not contain CCO"
        );
    }

    #[test]
    fn test_get_substruct_match_self() {
        let cco = make_mol_cco();
        let result = get_substruct_match(&cco, &cco);
        assert!(result.is_some(), "self-match should return Some");
        let result = result.unwrap();
        assert_eq!(result.atom_mapping.len(), 3);
        // Identity mapping: 0->0, 1->1, 2->2
        for (qa, ma) in result.atom_mapping.iter().enumerate() {
            assert_eq!(*ma, qa, "self-match should have identity mapping");
        }
    }

    #[test]
    fn test_get_substruct_match_cc_in_cco() {
        let cc = make_mol_cc();
        let cco = make_mol_cco();
        let result = get_substruct_match(&cco, &cc);
        assert!(result.is_some(), "CC should match in CCO");
    }

    #[test]
    fn test_get_substruct_match_no_match() {
        let c = make_mol_c();
        let cco = make_mol_cco();
        let result = get_substruct_match(&c, &cco);
        assert!(
            result.is_none(),
            "C should not match CCO (query larger than mol)"
        );
    }

    #[test]
    fn test_get_substruct_matches_cco_in_cco() {
        let cco = make_mol_cco();
        let matches = get_substruct_matches(&cco, &cco);
        assert!(!matches.is_empty(), "should find at least self-match");
    }

    #[test]
    fn test_substruct_coc_matches_cco() {
        // COC (dimethyl ether) should not match CCO (ethanol) — different topology.
        let coc = make_mol_coc();
        let cco = make_mol_cco();
        assert!(
            !has_substruct_match(&cco, &coc),
            "CCO should not match COC topology"
        );
        // But CO should match CCO (CO is a substructure of CCO).
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c, o, BondOrder::Single))
            .expect("add CO bond");
        let co = builder.build().expect("build CO");
        assert!(has_substruct_match(&cco, &co), "CCO should match CO");
    }

    #[test]
    fn test_has_substruct_match_empty_mol() {
        let empty = Molecule::new();
        let c = make_mol_c();
        assert!(
            !has_substruct_match(&empty, &c),
            "empty molecule should not match anything"
        );
        assert!(
            !has_substruct_match(&c, &empty),
            "molecule should not match empty query"
        );
    }

    #[test]
    fn test_substruct_match_params_max_matches() {
        let c = make_mol_c();
        let params = SubstructMatchParams {
            max_matches: 1,
            uniquify: true,
            use_chirality: false,
            specified_stereo_query_matches_unspecified: false,
            ..Default::default()
        };
        let matches = get_substruct_matches_with_params(&c, &c, &params);
        assert_eq!(matches.len(), 1, "max_matches=1 should return one match");
    }

    #[test]
    fn feature_smarts_substruct_matches_required_query_semantics() {
        let cases = [
            (
                "Donor",
                "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
                "CCO",
                vec![2],
                vec![2],
            ),
            (
                "Acceptor",
                "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
                "CC(=O)C",
                vec![2],
                vec![2],
            ),
            (
                "Aromatic",
                "[a]",
                "c1ccccc1",
                vec![0],
                vec![0, 1, 2, 3, 4, 5],
            ),
            ("Halogen", "[F,Cl,Br,I]", "CCl", vec![1], vec![1]),
            (
                "Basic",
                "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
                "[NH4+]",
                vec![0],
                vec![0],
            ),
            (
                "Acidic",
                "[$([C,S](=[O,S,P])-[O;H1,-1])]",
                "CC(=O)O",
                vec![1],
                vec![1],
            ),
        ];

        for (name, smarts, smiles, expected_first, expected_atoms) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, false)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = compile_query_fixture(smarts)
                .unwrap_or_else(|_| panic!("{name} SMARTS should build query molecule"));
            let matches = get_substruct_matches(&mol, &query);
            assert!(
                !matches.is_empty(),
                "{name} should produce at least one match"
            );
            assert_eq!(
                matches[0].atom_mapping, expected_first,
                "{name} first match atom mapping"
            );
            let mut atom_indices: Vec<usize> = matches
                .iter()
                .flat_map(|matched| matched.atom_mapping.iter().copied())
                .filter(|idx| *idx != NULL_NODE)
                .collect();
            atom_indices.sort_unstable();
            atom_indices.dedup();
            assert_eq!(atom_indices, expected_atoms, "{name} feature SMARTS");
        }
    }

    #[test]
    fn lipinski_hba_recursive_smarts_matches_rdkit_root_semantics() {
        const HBA: &str = "[$([O,S;H1;v2]-[!$(*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=!@[O,N,P,S])]),$([nH0X2,o,s;+0])]";
        let cases = [
            (
                "alcohol oxygen recursive branch",
                "[O,S;H1;v2]-[!$(*=[O,N,P,S])]",
                "CCO",
                vec![vec![2, 1]],
            ),
            (
                "carboxylic acid oxygen rejected by negated recursive neighbor",
                "[O,S;H1;v2]-[!$(*=[O,N,P,S])]",
                "CC(=O)O",
                Vec::<Vec<usize>>::new(),
            ),
            (
                "amine nitrogen recursive branch",
                "[N;v3;!$(N-*=!@[O,N,P,S])]",
                "CCN",
                vec![vec![2]],
            ),
            (
                "amide nitrogen rejected by mixed bond recursive query",
                "[N;v3;!$(N-*=!@[O,N,P,S])]",
                "CC(=O)N",
                Vec::<Vec<usize>>::new(),
            ),
            (
                "mixed single-double non-ring bond query",
                "N-*=!@[O,N,P,S]",
                "CC(=O)N",
                vec![vec![3, 1, 2]],
            ),
            ("full HBA ethanol", HBA, "CCO", vec![vec![2]]),
            ("full HBA carboxylic acid", HBA, "CC(=O)O", vec![vec![2]]),
            ("full HBA amide", HBA, "CC(=O)N", vec![vec![2]]),
            ("full HBA pyridine", HBA, "c1ccncc1", vec![vec![3]]),
            ("full HBA furan", HBA, "c1ccoc1", vec![vec![3]]),
        ];

        for (name, smarts, smiles, expected) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, true)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = compile_query_fixture(smarts)
                .unwrap_or_else(|_| panic!("{name} SMARTS should build"));
            let matches = get_substruct_matches(&mol, &query);
            let atom_mappings = matches
                .iter()
                .map(|matched| matched.atom_mapping.clone())
                .collect::<Vec<_>>();
            assert_eq!(atom_mappings, expected, "{name}");
        }
    }

    #[test]
    fn smarts_recursive_compiled_query() {
        let query = compile_query_fixture("[$(C=O)_101,$(C=O)_101]")
            .expect("recursive SMARTS should compile once during parsing");
        let molecule = Molecule::from_smiles("CC(=O)C").expect("acetone fixture");
        let mut cache = RecursiveQueryMatchCache::new();
        populate_recursive_query_match_cache(
            &molecule,
            &query,
            &SubstructMatchParams::default(),
            &mut cache,
        )
        .expect("compiled recursive queries should populate the match cache");

        assert_eq!(
            cache.len(),
            1,
            "equal serial numbers share one compiled result"
        );
        assert_eq!(
            get_substruct_matches(&molecule, &query)[0].atom_mapping,
            vec![1]
        );
    }

    #[test]
    fn smarts_match_recursive() {
        let molecule = Molecule::from_smiles("CC(=O)C").expect("acetone fixture");
        let rooted_query = compile_query_fixture("C=O")
            .expect("inner query")
            .with_prop("_queryRootAtom", "1");
        let mut cache = RecursiveQueryMatchCache::new();
        let starts = recursive_matcher(
            &molecule,
            &rooted_query,
            &SubstructMatchParams::default(),
            &mut cache,
        )
        .expect("rooted recursive matcher");
        assert_eq!(
            starts
                .iter()
                .enumerate()
                .filter_map(|(index, matched)| matched.then_some(index))
                .collect::<Vec<_>>(),
            vec![2]
        );

        let propane = Molecule::from_smiles("CCC").expect("propane fixture");
        let carbon = compile_query_fixture("C").expect("carbon query");
        let params = SubstructMatchParams {
            max_matches: 1,
            max_recursive_matches: 3,
            ..SubstructMatchParams::default()
        };
        let starts = recursive_matcher(
            &propane,
            &carbon,
            &params,
            &mut RecursiveQueryMatchCache::new(),
        )
        .expect("recursive match limit");
        assert_eq!(starts, vec![true, true, true]);
    }

    #[test]
    fn smarts_match_subqueries_execute() {
        let molecule = Molecule::from_smiles("CC(=O)C").expect("acetone fixture");
        let query =
            compile_query_fixture("[$(C=O)_101,$(C=O)_101]").expect("serial recursive query");
        let query_node = query.atoms()[0].query().expect("atom query tree");
        let mut cache = RecursiveQueryMatchCache::new();
        match_subqueries(
            &molecule,
            query_node,
            &SubstructMatchParams::default(),
            &mut cache,
        )
        .expect("execute recursive query tree");

        assert_eq!(cache.len(), 1, "equal serials reuse one result");
        assert_eq!(
            cache.get(&RecursiveQueryCacheKey::Serial(101)),
            Some(&vec![false, true, false, false])
        );
    }

    #[test]
    fn smarts_match_recursive_lock() {
        let molecule = Molecule::from_smiles("CC(=O)C").expect("acetone fixture");
        let query = compile_query_fixture("[$(C=O)]").expect("recursive SMARTS query");
        let enabled = SubstructMatchParams::default();
        assert_eq!(
            try_get_substruct_matches_with_params(&molecule, &query, &enabled)
                .expect("enabled recursive match")[0]
                .atom_mapping,
            vec![1]
        );

        let disabled = SubstructMatchParams {
            recursion_possible: false,
            ..SubstructMatchParams::default()
        };
        assert!(
            try_get_substruct_matches_with_params(&molecule, &query, &disabled)
                .expect("disabled recursion is a non-match")
                .is_empty()
        );

        assert_eq!(
            try_get_substruct_matches_with_params(&molecule, &query, &enabled)
                .expect("recursive state is rebuilt after scoped cleanup")[0]
                .atom_mapping,
            vec![1]
        );
    }

    #[test]
    fn smarts_match_entry() {
        let empty = Molecule::new();
        let carbon = Molecule::from_smiles("C").expect("carbon fixture");
        assert!(
            try_get_substruct_matches_with_params(
                &empty,
                &carbon,
                &SubstructMatchParams::default(),
            )
            .expect("empty target")
            .is_empty()
        );
        assert!(
            try_get_substruct_matches_with_params(
                &carbon,
                &empty,
                &SubstructMatchParams::default(),
            )
            .expect("empty query")
            .is_empty()
        );

        let ethane = Molecule::from_smiles("CC").expect("ethane fixture");
        assert!(
            try_get_substruct_matches_with_params(
                &carbon,
                &ethane,
                &SubstructMatchParams::default(),
            )
            .expect("oversized query")
            .is_empty()
        );

        let propane = Molecule::from_smiles("CCC").expect("propane fixture");
        let params = SubstructMatchParams {
            max_matches: 1,
            uniquify: false,
            ..SubstructMatchParams::default()
        };
        let matches = try_get_substruct_matches_with_params(&propane, &ethane, &params)
            .expect("bounded entry match");
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].atom_mapping.len(), ethane.num_atoms());
        assert_eq!(matches[0].atom_mapping, vec![0, 1]);
    }

    #[test]
    fn smarts_unsupported_query_errors() {
        let target = Molecule::from_smiles("CC").expect("target fixture");
        let params = SubstructMatchParams::default();

        let mut atom_builder = MoleculeBuilder::new();
        atom_builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_query(
            crate::QueryNode::and(vec![
                crate::QueryNode::predicate(AtomQueryPredicate::Any),
                crate::QueryNode::not(crate::QueryNode::predicate(
                    AtomQueryPredicate::UnsupportedFeature("unsupported atom leaf"),
                )),
            ]),
        ));
        let atom_query = atom_builder.build().expect("atom query fixture");
        assert_eq!(
            try_get_substruct_matches_with_params(&target, &atom_query, &params),
            Err(SubstructMatchError::Unsupported {
                branch: "unsupported atom leaf",
                rdkit_function: "QueryAtom::Match",
            })
        );

        let mut bond_builder = MoleculeBuilder::new();
        let begin = bond_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let end = bond_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        bond_builder
            .add_bond(
                crate::BondSpec::new(begin, end, BondOrder::Single).with_query(
                    crate::QueryNode::or(vec![
                        crate::QueryNode::predicate(BondQueryPredicate::Any),
                        crate::QueryNode::predicate(BondQueryPredicate::UnsupportedFeature(
                            "unsupported bond leaf",
                        )),
                    ]),
                ),
            )
            .expect("bond query edge");
        let bond_query = bond_builder.build().expect("bond query fixture");
        assert_eq!(
            try_get_substruct_matches_with_params(&target, &bond_query, &params),
            Err(SubstructMatchError::Unsupported {
                branch: "unsupported bond leaf",
                rdkit_function: "QueryBond::Match",
            })
        );

        let mut inner_builder = MoleculeBuilder::new();
        inner_builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_query(
            crate::QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature(
                "unsupported recursive leaf",
            )),
        ));
        let recursive_query = crate::search::query::RecursiveStructureQuery::from_molecule(
            inner_builder.build().expect("inner query fixture"),
            0,
        );
        let mut outer_builder = MoleculeBuilder::new();
        outer_builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY).with_query(
            crate::QueryNode::predicate(AtomQueryPredicate::RecursiveSmarts(recursive_query)),
        ));
        let outer_query = outer_builder.build().expect("outer query fixture");
        assert_eq!(
            try_get_substruct_matches_with_params(&target, &outer_query, &params),
            Err(SubstructMatchError::Unsupported {
                branch: "unsupported recursive leaf",
                rdkit_function: "QueryAtom::Match",
            })
        );
    }

    #[test]
    fn smarts_substruct_property_compat() {
        fn properties(entries: &[(&str, &str)]) -> BTreeMap<String, String> {
            entries
                .iter()
                .map(|(key, value)| ((*key).to_owned(), (*value).to_owned()))
                .collect()
        }

        let requested = vec!["test_prop".to_owned()];
        let empty = properties(&[]);
        let one = properties(&[("test_prop", "1")]);
        let same = properties(&[("test_prop", "1"), ("ignored", "left")]);
        let different = properties(&[("test_prop", "2")]);
        let unrequested_difference = properties(&[("ignored", "right")]);

        assert!(property_compat(&empty, &empty, &requested));
        assert!(property_compat(&one, &same, &requested));
        assert!(!property_compat(&one, &different, &requested));
        assert!(!property_compat(&one, &empty, &requested));
        assert!(!property_compat(&empty, &one, &requested));
        assert!(property_compat(&empty, &unrequested_difference, &requested));
        assert!(property_compat(&one, &different, &[]));
        assert!(!property_compat(
            &properties(&[("first", "same"), ("second", "left")]),
            &properties(&[("first", "same"), ("second", "right")]),
            &["first".to_owned(), "second".to_owned()],
        ));
    }

    #[test]
    fn smarts_substruct_atom_compat() {
        fn carbon_chain(atom_count: usize, first_property: Option<&str>) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            let mut atoms = Vec::with_capacity(atom_count);
            for atom_index in 0..atom_count {
                let mut atom = crate::AtomSpec::new(crate::Element::C);
                if atom_index == 0
                    && let Some(value) = first_property
                {
                    atom = atom.with_prop("test_prop", value);
                }
                atoms.push(builder.add_atom(atom));
            }
            for pair in atoms.windows(2) {
                builder
                    .add_bond(crate::BondSpec::new(pair[0], pair[1], BondOrder::Single))
                    .expect("chain bond");
            }
            builder.build().expect("carbon chain")
        }

        let mut property_params = SubstructMatchParams::default();
        property_params.atom_properties = vec!["test_prop".to_owned()];
        let cases = [
            (None, None, 7),
            (Some("1"), Some("1"), 1),
            (Some("1"), None, 6),
            (None, Some("1"), 0),
            (Some("1"), Some("2"), 0),
        ];
        for (target_property, query_property, expected) in cases {
            let target = carbon_chain(9, target_property);
            let query = carbon_chain(3, query_property);
            assert_eq!(
                get_substruct_matches_with_params(&target, &query, &property_params).len(),
                expected,
                "target={target_property:?}, query={query_property:?}"
            );
        }

        let query_query_molecule = |predicate| {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(
                crate::AtomSpec::new(crate::Element::C)
                    .with_query(crate::QueryNode::predicate(predicate)),
            );
            builder.build().expect("single query atom")
        };
        let query = query_query_molecule(AtomQueryPredicate::AtomicNumber(6));
        let target = query_query_molecule(AtomQueryPredicate::AtomicNumber(8));
        assert!(has_substruct_match(&target, &query));
        let mut query_query_params = SubstructMatchParams::default();
        query_query_params.use_query_query_matches = true;
        assert!(get_substruct_matches_with_params(&target, &query, &query_query_params).is_empty());

        let carbon = carbon_chain(1, None);
        let mut oxygen_builder = MoleculeBuilder::new();
        oxygen_builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let oxygen = oxygen_builder.build().expect("oxygen query");

        let mut callback_params = SubstructMatchParams::default();
        let expected_query_index = 0;
        callback_params.extra_atom_check = Some(Arc::new(move |_, query_atom, _, mol_atom| {
            query_atom.id().index() == expected_query_index && mol_atom.atomic_number() == 6
        }));
        callback_params.extra_atom_check_overrides_default_check = true;
        assert_eq!(
            get_substruct_matches_with_params(&carbon, &oxygen, &callback_params).len(),
            1
        );

        callback_params.extra_atom_check_overrides_default_check = false;
        assert!(get_substruct_matches_with_params(&carbon, &oxygen, &callback_params).is_empty());

        callback_params.extra_atom_check = Some(Arc::new(|_, _, _, _| false));
        assert!(get_substruct_matches_with_params(&carbon, &carbon, &callback_params).is_empty());
    }

    #[test]
    fn smarts_match_atom_coords() {
        fn one_atom_with_conformers(conformers: &[(usize, [f64; 3])]) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            for &(id, position) in conformers {
                builder
                    .add_conformer(crate::Conformer3D::new(id, vec![position], true))
                    .expect("add conformer");
            }
            builder.build().expect("coordinate fixture")
        }

        let query = one_atom_with_conformers(&[(0, [0.0, 0.0, 0.1]), (7, [5.0, 0.0, 0.0])]);
        let target = one_atom_with_conformers(&[(0, [0.0, 0.0, 0.0]), (9, [5.1, 0.0, 0.0])]);
        let missing = one_atom_with_conformers(&[]);

        let default_matcher = AtomCoordsMatchFunctor::default();
        assert!(!default_matcher.matches(&query, &query.atoms()[0], &target, &target.atoms()[0],));
        assert!(
            !default_matcher.matches(&query, &query.atoms()[0], &missing, &missing.atoms()[0],)
        );

        let matcher = AtomCoordsMatchFunctor::new(9, 7, 0.15);
        assert!(matcher.matches(&query, &query.atoms()[0], &target, &target.atoms()[0],));
        let mut params = SubstructMatchParams::default();
        params.extra_atom_check = Some(Arc::new(move |query_mol, query_atom, mol, mol_atom| {
            matcher.matches(query_mol, query_atom, mol, mol_atom)
        }));
        assert_eq!(
            try_get_substruct_matches_with_params(&target, &query, &params)
                .expect("coordinate-constrained match")
                .len(),
            1
        );
    }

    #[test]
    fn smarts_substruct_chiral_atom_compat() {
        fn atom(element: crate::Element, cip: Option<&str>) -> Molecule {
            let mut spec = crate::AtomSpec::new(element);
            if let Some(cip) = cip {
                spec = spec.with_prop("_CIPCode", cip);
            }
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(spec);
            builder.build().expect("single atom fixture")
        }

        let carbon = atom(crate::Element::C, None);
        let oxygen = atom(crate::Element::O, None);
        assert!(!chiral_atom_compat(
            &carbon.atoms()[0],
            &carbon,
            &oxygen.atoms()[0],
            &oxygen,
        ));

        assert!(chiral_atom_compat(
            &carbon.atoms()[0],
            &carbon,
            &carbon.atoms()[0],
            &carbon,
        ));

        let carbon_r = atom(crate::Element::C, Some("R"));
        let another_carbon_r = atom(crate::Element::C, Some("R"));
        let carbon_s = atom(crate::Element::C, Some("S"));
        assert!(chiral_atom_compat(
            &carbon_r.atoms()[0],
            &carbon_r,
            &another_carbon_r.atoms()[0],
            &another_carbon_r,
        ));
        assert!(!chiral_atom_compat(
            &carbon_r.atoms()[0],
            &carbon_r,
            &carbon_s.atoms()[0],
            &carbon_s,
        ));
        assert!(!chiral_atom_compat(
            &carbon_r.atoms()[0],
            &carbon_r,
            &carbon.atoms()[0],
            &carbon,
        ));
        assert!(!chiral_atom_compat(
            &carbon.atoms()[0],
            &carbon,
            &carbon_r.atoms()[0],
            &carbon_r,
        ));
    }

    #[test]
    fn smarts_substruct_bond_compat() {
        fn two_atom_molecule(
            begin: crate::Element,
            end: crate::Element,
            bond: crate::BondSpec,
        ) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(crate::AtomSpec::new(begin));
            builder.add_atom(crate::AtomSpec::new(end));
            builder.add_bond(bond).expect("two-atom bond");
            builder.build().expect("two-atom molecule")
        }

        fn compatible(query: &Molecule, target: &Molecule, params: &SubstructMatchParams) -> bool {
            bond_compat(
                &query.bonds()[0],
                query,
                &target.bonds()[0],
                target,
                params,
                &build_query_match_context(target),
            )
        }

        let single = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            ),
        );
        let double = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Double,
            ),
        );
        let unspecified = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Unspecified,
            ),
        );
        assert!(!compatible(
            &single,
            &double,
            &SubstructMatchParams::default()
        ));
        assert!(compatible(
            &unspecified,
            &double,
            &SubstructMatchParams::default()
        ));

        let aromatic = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Aromatic,
            )
            .with_aromatic(true),
        );
        let conjugated_single = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            )
            .with_conjugated(true),
        );
        let mut params = SubstructMatchParams::default();
        params.aromatic_matches_conjugated = true;
        assert!(compatible(&aromatic, &conjugated_single, &params));
        assert!(!compatible(&aromatic, &single, &params));
        params.aromatic_matches_conjugated = false;
        params.aromatic_matches_single_or_double = true;
        assert!(compatible(&aromatic, &single, &params));
        assert!(compatible(&aromatic, &double, &params));

        let query_single = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            )
            .with_query(crate::QueryNode::predicate(BondQueryPredicate::Order(
                BondOrder::Single,
            ))),
        );
        let query_double = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Double,
            )
            .with_query(crate::QueryNode::predicate(BondQueryPredicate::Order(
                BondOrder::Double,
            ))),
        );
        let mut query_query_params = SubstructMatchParams::default();
        query_query_params.use_query_query_matches = true;
        assert!(!compatible(
            &query_single,
            &query_double,
            &query_query_params
        ));

        let property_single = two_atom_molecule(
            crate::Element::C,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Single,
            )
            .with_prop("test_prop", "left"),
        );
        let mut property_params = SubstructMatchParams::default();
        property_params.bond_properties = vec!["test_prop".to_owned()];
        assert!(!compatible(&property_single, &single, &property_params));

        let mut callback_params = SubstructMatchParams::default();
        callback_params.extra_bond_check = Some(Arc::new(|query, target| {
            query.order() == BondOrder::Single && target.order() == BondOrder::Double
        }));
        callback_params.extra_bond_check_overrides_default_check = true;
        assert!(compatible(&single, &double, &callback_params));
        callback_params.extra_bond_check_overrides_default_check = false;
        assert!(!compatible(&single, &double, &callback_params));
        callback_params.extra_bond_check = Some(Arc::new(|_, _| false));
        assert!(!compatible(&single, &single, &callback_params));

        let dative_cn = two_atom_molecule(
            crate::Element::C,
            crate::Element::N,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Dative,
            ),
        );
        let dative_nc = two_atom_molecule(
            crate::Element::N,
            crate::Element::C,
            crate::BondSpec::new(
                crate::AtomId::new(0),
                crate::AtomId::new(1),
                BondOrder::Dative,
            ),
        );
        assert!(compatible(
            &dative_cn,
            &dative_cn,
            &SubstructMatchParams::default()
        ));
        assert!(!compatible(
            &dative_cn,
            &dative_nc,
            &SubstructMatchParams::default()
        ));
    }

    #[test]
    fn smarts_substruct_remove_duplicates() {
        let result = |atom_mapping: &[usize], bond_mapping: &[usize]| SubstructMatchResult {
            atom_mapping: atom_mapping.to_vec(),
            bond_mapping: bond_mapping.to_vec(),
        };
        let first = result(&[0, 1, 2, 3], &[10, 11, 12]);
        let same_atom_set_different_path = result(&[3, 2, 1, 0], &[20, 21, 22]);
        let distinct = result(&[0, 1, 2, 4], &[30, 31, 32]);
        let repeated_distinct = result(&[4, 2, 1, 0], &[40, 41, 42]);
        let mut matches = vec![
            first.clone(),
            same_atom_set_different_path,
            distinct.clone(),
            repeated_distinct,
        ];

        remove_duplicates(&mut matches, 5);

        assert_eq!(matches, vec![first, distinct]);
        assert_eq!(matches.capacity(), matches.len());
    }

    fn core_substitution_fixtures() -> (Molecule, Molecule, Vec<SubstructMatchResult>) {
        let mut molecule_builder = MoleculeBuilder::new();
        let carbon_zero = molecule_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = molecule_builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let carbon_two = molecule_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        molecule_builder
            .add_bond(crate::BondSpec::new(
                carbon_zero,
                hydrogen,
                BondOrder::Single,
            ))
            .expect("C-H bond");
        molecule_builder
            .add_bond(crate::BondSpec::new(
                carbon_zero,
                carbon_two,
                BondOrder::Single,
            ))
            .expect("C-C bond");
        let molecule = molecule_builder.build().expect("target fixture");

        let mut query_builder = MoleculeBuilder::new();
        let dummy = query_builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let carbon = query_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        query_builder
            .add_bond(crate::BondSpec::new(dummy, carbon, BondOrder::Single))
            .expect("query bond");
        let query = query_builder.build().expect("query fixture");

        let hydrogen_match = SubstructMatchResult {
            atom_mapping: vec![1, 0],
            bond_mapping: vec![0],
        };
        let substituted_match = SubstructMatchResult {
            atom_mapping: vec![2, 0],
            bond_mapping: vec![1],
        };
        (molecule, query, vec![hydrogen_match, substituted_match])
    }

    #[test]
    fn smarts_substruct_get_most_substituted_core_match() {
        let (molecule, query, matches) = core_substitution_fixtures();
        assert_eq!(
            get_most_substituted_core_match(&molecule, &query, &matches),
            &matches[1]
        );
        assert!(core_substitution_score(&molecule, &query, &matches[1]) < 1.0);
        assert!(core_substitution_score(&molecule, &query, &matches[0]) >= 1.0);
    }

    #[test]
    #[should_panic(expected = "matches must not be empty")]
    fn smarts_substruct_get_most_substituted_core_match_rejects_empty() {
        let (molecule, query, _) = core_substitution_fixtures();
        let _ = get_most_substituted_core_match(&molecule, &query, &[]);
    }

    #[test]
    fn smarts_substruct_sort_matches_by_degree_of_core_substitution() {
        let (molecule, query, matches) = core_substitution_fixtures();
        let sorted = sort_matches_by_degree_of_core_substitution(&molecule, &query, &matches);
        assert_eq!(sorted, vec![matches[1].clone(), matches[0].clone()]);
        assert_eq!(matches[0].atom_mapping, vec![1, 0]);
    }

    #[test]
    fn smarts_substruct_is_atom_terminal_r_group_or_query_hydrogen() {
        let (_, terminal_dummy_query, _) = core_substitution_fixtures();
        assert!(is_atom_terminal_r_group_or_query_hydrogen(
            &terminal_dummy_query,
            0
        ));
        assert!(!is_atom_terminal_r_group_or_query_hydrogen(
            &terminal_dummy_query,
            1
        ));

        let mut hydrogen_query_builder = MoleculeBuilder::new();
        hydrogen_query_builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY).with_query(
            crate::QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
        ));
        let hydrogen_query = hydrogen_query_builder
            .build()
            .expect("hydrogen query fixture");
        assert!(is_atom_terminal_r_group_or_query_hydrogen(
            &hydrogen_query,
            0
        ));

        let mut nonterminal_dummy_builder = MoleculeBuilder::new();
        let dummy = nonterminal_dummy_builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let carbon_one =
            nonterminal_dummy_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_two =
            nonterminal_dummy_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        nonterminal_dummy_builder
            .add_bond(crate::BondSpec::new(dummy, carbon_one, BondOrder::Single))
            .expect("first dummy bond");
        nonterminal_dummy_builder
            .add_bond(crate::BondSpec::new(dummy, carbon_two, BondOrder::Single))
            .expect("second dummy bond");
        let nonterminal_dummy = nonterminal_dummy_builder
            .build()
            .expect("nonterminal dummy fixture");
        assert!(!is_atom_terminal_r_group_or_query_hydrogen(
            &nonterminal_dummy,
            0
        ));
    }

    #[test]
    fn smarts_substruct_update_substruct_match_params_from_j_s_o_n() {
        let mut params = SubstructMatchParams::default();
        params.max_matches = 77;
        update_substruct_match_params_from_json(&mut params, "").expect("empty JSON no-op");
        assert_eq!(params.max_matches, 77);

        update_substruct_match_params_from_json(
            &mut params,
            r#"{
                "useChirality": true,
                "useEnhancedStereo": "true",
                "aromaticMatchesConjugated": true,
                "useQueryQueryMatches": "1",
                "recursionPossible": false,
                "uniquify": "false",
                "maxMatches": 12,
                "maxRecursiveMatches": "34",
                "numThreads": -2,
                "specifiedStereoQueryMatchesUnspecified": true,
                "aromaticMatchesSingleOrDouble": "true",
                "unknownOption": "ignored"
            }"#,
        )
        .expect("source JSON fields");
        assert!(params.use_chirality);
        assert!(params.use_enhanced_stereo);
        assert!(params.aromatic_matches_conjugated);
        assert!(params.use_query_query_matches);
        assert!(!params.recursion_possible);
        assert!(!params.uniquify);
        assert_eq!(params.max_matches, 12);
        assert_eq!(params.max_recursive_matches, 34);
        assert_eq!(params.num_threads, -2);
        assert!(params.specified_stereo_query_matches_unspecified);
        assert!(params.aromatic_matches_single_or_double);

        let before = params.clone();
        assert!(
            update_substruct_match_params_from_json(
                &mut params,
                r#"{"useChirality": false, "maxMatches": "invalid"}"#,
            )
            .is_err()
        );
        assert_eq!(params.use_chirality, before.use_chirality);
        assert_eq!(params.max_matches, before.max_matches);
    }

    #[test]
    fn smarts_substruct_substruct_match_params_to_j_s_o_n() {
        let mut params = SubstructMatchParams::default();
        params.use_chirality = true;
        params.use_enhanced_stereo = true;
        params.aromatic_matches_conjugated = true;
        params.use_query_query_matches = true;
        params.recursion_possible = false;
        params.uniquify = false;
        params.max_matches = 12;
        params.max_recursive_matches = 34;
        params.num_threads = -2;
        params.specified_stereo_query_matches_unspecified = true;
        params.aromatic_matches_single_or_double = true;
        params.atom_properties.push("notSerialized".to_owned());
        params.bond_properties.push("notSerialized".to_owned());

        let json = substruct_match_params_to_json(&params);
        let value: serde_json::Value = serde_json::from_str(&json).expect("writer JSON");
        let object = value.as_object().expect("JSON object");
        assert_eq!(object.len(), 11);
        assert_eq!(object["useChirality"], "true");
        assert_eq!(object["maxMatches"], "12");
        assert_eq!(object["numThreads"], "-2");
        assert!(!object.contains_key("atomProperties"));
        assert!(!object.contains_key("bondProperties"));

        let mut roundtrip = SubstructMatchParams::default();
        update_substruct_match_params_from_json(&mut roundtrip, &json).expect("writer roundtrip");
        assert_eq!(roundtrip.use_chirality, params.use_chirality);
        assert_eq!(roundtrip.use_enhanced_stereo, params.use_enhanced_stereo);
        assert_eq!(
            roundtrip.aromatic_matches_conjugated,
            params.aromatic_matches_conjugated
        );
        assert_eq!(
            roundtrip.use_query_query_matches,
            params.use_query_query_matches
        );
        assert_eq!(roundtrip.recursion_possible, params.recursion_possible);
        assert_eq!(roundtrip.uniquify, params.uniquify);
        assert_eq!(roundtrip.max_matches, params.max_matches);
        assert_eq!(
            roundtrip.max_recursive_matches,
            params.max_recursive_matches
        );
        assert_eq!(roundtrip.num_threads, params.num_threads);
        assert_eq!(
            roundtrip.specified_stereo_query_matches_unspecified,
            params.specified_stereo_query_matches_unspecified
        );
        assert_eq!(
            roundtrip.aromatic_matches_single_or_double,
            params.aromatic_matches_single_or_double
        );
    }

    #[test]
    fn smarts_substruct_match_parameters() {
        let params = SubstructMatchParams::default();
        assert!(!params.use_chirality);
        assert!(!params.use_enhanced_stereo);
        assert!(!params.use_generic_matchers);
        assert!(params.recursion_possible);
        assert!(params.uniquify);
        assert_eq!(params.max_matches, 1000);
        assert_eq!(params.max_recursive_matches, 1000);
        assert_eq!(params.num_threads, 1);
        assert!(params.atom_properties.is_empty());
        assert!(params.bond_properties.is_empty());
        assert!(params.extra_atom_check.is_none());
        assert!(params.extra_bond_check.is_none());
        assert!(params.extra_final_check.is_none());

        assert_eq!(
            check_substruct_match_overload_support(SubstructMatchOverload::Molecule),
            Ok(())
        );
        for (overload, expected_branch) in [
            (
                SubstructMatchOverload::MolBundle,
                "MolBundle substructure-match overloads",
            ),
            (
                SubstructMatchOverload::ResonanceMolSupplier,
                "resonance substructure-match overload",
            ),
            (
                SubstructMatchOverload::SubstructLibrary,
                "SubstructLibrary search overloads",
            ),
        ] {
            assert!(matches!(
                check_substruct_match_overload_support(overload),
                Err(SubstructMatchError::Unsupported { branch, .. }) if branch == expected_branch
            ));
        }
    }

    #[test]
    fn maccs_patterns_substruct_matches_required_topology_semantics() {
        let cases = [
            (
                "four-membered ring",
                "*1~*~*~*~1",
                "C1CCC1",
                true,
                vec![0, 1, 2, 3],
            ),
            (
                "four-membered ring rejects chain",
                "*1~*~*~*~1",
                "CCCC",
                false,
                Vec::new(),
            ),
            ("ring bond", "*@*(@*)@*", "C12CC1C2", true, vec![0, 2, 1, 3]),
            (
                "non-ring oxygen bridge",
                "*!@[#8]!@*",
                "COC",
                true,
                vec![0, 1, 2],
            ),
            (
                "branch degree",
                "*~*(~*)(~*)~*",
                "CC(C)(C)C",
                true,
                vec![0, 1, 2, 3, 4],
            ),
            (
                "recursive ring closure",
                "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
                "CC1CC1",
                true,
                vec![0],
            ),
        ];

        for (name, smarts, smiles, expected_match, expected_first) in cases {
            let mol = Molecule::from_smiles_with_sanitize(smiles, true)
                .unwrap_or_else(|_| panic!("{name} molecule should parse"));
            let query = compile_query_fixture(smarts)
                .unwrap_or_else(|_| panic!("{name} MACCS SMARTS should build query"));
            let matches = get_substruct_matches(&mol, &query);
            assert_eq!(
                !matches.is_empty(),
                expected_match,
                "{name} match truth value"
            );
            if expected_match {
                assert_eq!(
                    matches[0].atom_mapping, expected_first,
                    "{name} first match"
                );
            }
        }
    }

    struct MaccsPatternGolden {
        bit: u16,
        smarts: &'static str,
        smiles: &'static str,
        first_match: &'static [usize],
    }

    fn rdkit_maccs_pattern_positive_goldens() -> &'static [MaccsPatternGolden] {
        // RDKit source: MACCS.cpp::Patterns initializes these SMARTS strings
        // with `RDKit::SmartsToMol(...)`. `first_match` values were generated
        // from pinned RDKit 2026.03.1 using `Mol.GetSubstructMatch()`.
        &[
            MaccsPatternGolden {
                bit: 8,
                smarts: "[!#6!#1]1~*~*~*~1",
                smiles: "O1CCC1",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 11,
                smarts: "*1~*~*~*~1",
                smiles: "C1CCC1",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 13,
                smarts: "[#8]~[#7](~[#6])~[#6]",
                smiles: "ON(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 14,
                smarts: "[#16]-[#16]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 15,
                smarts: "[#8]~[#6](~[#8])~[#8]",
                smiles: "O=C(O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 16,
                smarts: "[!#6!#1]1~*~*~1",
                smiles: "O1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 17,
                smarts: "[#6]#[#6]",
                smiles: "C#C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 19,
                smarts: "*1~*~*~*~*~*~*~1",
                smiles: "C1CCCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5, 6],
            },
            MaccsPatternGolden {
                bit: 20,
                smarts: "[#14]",
                smiles: "[SiH4]",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 21,
                smarts: "[#6]=[#6](~[!#6!#1])~[!#6!#1]",
                smiles: "C=C(O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 22,
                smarts: "*1~*~*~1",
                smiles: "C1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 23,
                smarts: "[#7]~[#6](~[#8])~[#8]",
                smiles: "NC(=O)O",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 24,
                smarts: "[#7]-[#8]",
                smiles: "ON(C)C",
                first_match: &[1, 0],
            },
            MaccsPatternGolden {
                bit: 25,
                smarts: "[#7]~[#6](~[#7])~[#7]",
                smiles: "NC(N)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 26,
                smarts: "[#6]=@[#6](@*)@*",
                smiles: "C1=C2CCCC2C1",
                first_match: &[0, 1, 2, 5],
            },
            MaccsPatternGolden {
                bit: 28,
                smarts: "[!#6!#1]~[CH2]~[!#6!#1]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 30,
                smarts: "[#6]~[!#6!#1](~[#6])(~[#6])~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 31,
                smarts: "[!#6!#1]~[F,Cl,Br,I]",
                smiles: "N[Pt](Cl)(Cl)N",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 32,
                smarts: "[#6]~[#16]~[#7]",
                smiles: "CSN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 33,
                smarts: "[#7]~[#16]",
                smiles: "CSN",
                first_match: &[2, 1],
            },
            MaccsPatternGolden {
                bit: 34,
                smarts: "[CH2]=*",
                smiles: "C=C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 36,
                smarts: "[#16R]",
                smiles: "S1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 37,
                smarts: "[#7]~[#6](~[#8])~[#7]",
                smiles: "NC(=O)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 38,
                smarts: "[#7]~[#6](~[#6])~[#7]",
                smiles: "NC(C)N",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 39,
                smarts: "[#8]~[#16](~[#8])~[#8]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 40,
                smarts: "[#16]-[#8]",
                smiles: "CSO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 41,
                smarts: "[#6]#[#7]",
                smiles: "C#N",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 43,
                smarts: "[!#6!#1!H0]~*~[!#6!#1!H0]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 44,
                smarts: "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]",
                smiles: "[SeH2]",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 45,
                smarts: "[#6]=[#6]~[#7]",
                smiles: "C=CN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 47,
                smarts: "[#16]~*~[#7]",
                smiles: "SCN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 48,
                smarts: "[#8]~[!#6!#1](~[#8])~[#8]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 49,
                smarts: "[!+0]",
                smiles: "O=N(=O)O",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 50,
                smarts: "[#6]=[#6](~[#6])~[#6]",
                smiles: "CC(C)=C",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 51,
                smarts: "[#6]~[#16]~[#8]",
                smiles: "CSO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 52,
                smarts: "[#7]~[#7]",
                smiles: "NNO",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 53,
                smarts: "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]",
                smiles: "NCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 54,
                smarts: "[!#6!#1!H0]~*~*~[!#6!#1!H0]",
                smiles: "NCCN",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 55,
                smarts: "[#8]~[#16]~[#8]",
                smiles: "CS(=O)(=O)C",
                first_match: &[2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 56,
                smarts: "[#8]~[#7](~[#8])~[#6]",
                smiles: "ON(O)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 57,
                smarts: "[#8R]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 58,
                smarts: "[!#6!#1]~[#16]~[!#6!#1]",
                smiles: "CS(=O)(=O)C",
                first_match: &[2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 59,
                smarts: "[#16]!:*:*",
                smiles: "Sc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 60,
                smarts: "[#16]=[#8]",
                smiles: "CS(=O)C",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 61,
                smarts: "*~[#16](~*)~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 62,
                smarts: "*@*!@*@*",
                smiles: "C1CC1C1CC1",
                first_match: &[0, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 63,
                smarts: "[#7]=[#8]",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 64,
                smarts: "*@*!@[#16]",
                smiles: "Sc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 65,
                smarts: "c:n",
                smiles: "c1ncccc1",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 66,
                smarts: "[#6]~[#6](~[#6])(~[#6])~*",
                smiles: "CC(C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 67,
                smarts: "[!#6!#1]~[#16]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 68,
                smarts: "[!#6!#1!H0]~[!#6!#1!H0]",
                smiles: "NNO",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 69,
                smarts: "[!#6!#1]~[!#6!#1!H0]",
                smiles: "CSN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 70,
                smarts: "[!#6!#1]~[#7]~[!#6!#1]",
                smiles: "NNO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 71,
                smarts: "[#7]~[#8]",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 72,
                smarts: "[#8]~*~*~[#8]",
                smiles: "OCCO",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 73,
                smarts: "[#16]=*",
                smiles: "CS(=O)C",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 74,
                smarts: "[CH3]~*~[CH3]",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 75,
                smarts: "*!@[#7]@*",
                smiles: "CN1CC1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 76,
                smarts: "[#6]=[#6](~*)~*",
                smiles: "CC(C)=C",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 77,
                smarts: "[#7]~*~[#7]",
                smiles: "NC(=O)N",
                first_match: &[0, 1, 3],
            },
            MaccsPatternGolden {
                bit: 78,
                smarts: "[#6]=[#7]",
                smiles: "C=N",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 79,
                smarts: "[#7]~*~*~[#7]",
                smiles: "NCCN",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 80,
                smarts: "[#7]~*~*~*~[#7]",
                smiles: "NCCCN",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 81,
                smarts: "[#16]~*(~*)~*",
                smiles: "Sc1ccccc1",
                first_match: &[0, 1, 2, 6],
            },
            MaccsPatternGolden {
                bit: 82,
                smarts: "*~[CH2]~[!#6!#1!H0]",
                smiles: "CCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 83,
                smarts: "[!#6!#1]1~*~*~*~*~1",
                smiles: "O1CCCC1",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 84,
                smarts: "[NH2]",
                smiles: "CCN",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 85,
                smarts: "[#6]~[#7](~[#6])~[#6]",
                smiles: "CN(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 86,
                smarts: "[C;H2,H3][!#6!#1][C;H2,H3]",
                smiles: "COC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 87,
                smarts: "[F,Cl,Br,I]!@*@*",
                smiles: "Clc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 89,
                smarts: "[#8]~*~*~*~[#8]",
                smiles: "OCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 90,
                smarts: "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
                smiles: "N1CCC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 91,
                smarts: "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
                smiles: "NCCCCCN",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 92,
                smarts: "[#8]~[#6](~[#7])~[#6]",
                smiles: "CC(=O)N",
                first_match: &[2, 1, 3, 0],
            },
            MaccsPatternGolden {
                bit: 93,
                smarts: "[!#6!#1]~[CH3]",
                smiles: "COC",
                first_match: &[1, 0],
            },
            MaccsPatternGolden {
                bit: 94,
                smarts: "[!#6!#1]~[#7]",
                smiles: "CSN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 95,
                smarts: "[#7]~*~*~[#8]",
                smiles: "NCCO",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 96,
                smarts: "*1~*~*~*~*~1",
                smiles: "C1CCCC1",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 97,
                smarts: "[#7]~*~*~*~[#8]",
                smiles: "NCCCO",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 98,
                smarts: "[!#6!#1]1~*~*~*~*~*~1",
                smiles: "O1CCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 99,
                smarts: "[#6]=[#6]",
                smiles: "C=C",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 100,
                smarts: "*~[CH2]~[#7]",
                smiles: "CCN",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 101,
                smarts: "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
                smiles: "C1CCCCCCC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 102,
                smarts: "[!#6!#1]~[#8]",
                smiles: "CSO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 104,
                smarts: "[!#6!#1!H0]~*~[CH2]~*",
                smiles: "CCCO",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 105,
                smarts: "*@*(@*)@*",
                smiles: "C12CC1C2",
                first_match: &[0, 2, 1, 3],
            },
            MaccsPatternGolden {
                bit: 106,
                smarts: "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]",
                smiles: "COS(=O)(=O)O",
                first_match: &[1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 107,
                smarts: "[F,Cl,Br,I]~*(~*)~*",
                smiles: "CC(C)(C)Cl",
                first_match: &[4, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 108,
                smarts: "[CH3]~*~*~*~[CH2]~*",
                smiles: "CCCCCC",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 109,
                smarts: "*~[CH2]~[#8]",
                smiles: "CCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 110,
                smarts: "[#7]~[#6]~[#8]",
                smiles: "CC(=O)N",
                first_match: &[3, 1, 2],
            },
            MaccsPatternGolden {
                bit: 111,
                smarts: "[#7]~*~[CH2]~*",
                smiles: "CCCN",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 112,
                smarts: "*~*(~*)(~*)~*",
                smiles: "CC(C)(C)C",
                first_match: &[0, 1, 2, 3, 4],
            },
            MaccsPatternGolden {
                bit: 113,
                smarts: "[#8]!:*:*",
                smiles: "Oc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 114,
                smarts: "[CH3]~[CH2]~*",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 115,
                smarts: "[CH3]~*~[CH2]~*",
                smiles: "CCCC",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 116,
                smarts: "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]",
                smiles: "CCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 117,
                smarts: "[#7]~*~[#8]",
                smiles: "CC(=O)N",
                first_match: &[3, 1, 2],
            },
            MaccsPatternGolden {
                bit: 118,
                smarts: "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]",
                smiles: "CCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 119,
                smarts: "[#7]=*",
                smiles: "N=O",
                first_match: &[0, 1],
            },
            MaccsPatternGolden {
                bit: 120,
                smarts: "[!#6R]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 121,
                smarts: "[#7R]",
                smiles: "N1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 122,
                smarts: "*~[#7](~*)~*",
                smiles: "ON(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 123,
                smarts: "[#8]~[#6]~[#8]",
                smiles: "OCO",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 124,
                smarts: "[!#6!#1]~[!#6!#1]",
                smiles: "CSSC",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 126,
                smarts: "*!@[#8]!@*",
                smiles: "COC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 127,
                smarts: "*@*!@[#8]",
                smiles: "Oc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 128,
                smarts: "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
                smiles: "CCCCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 129,
                smarts: "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
                smiles: "CCCCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 131,
                smarts: "[!#6!#1!H0]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 132,
                smarts: "[#8]~*~[CH2]~*",
                smiles: "CCCO",
                first_match: &[3, 2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 133,
                smarts: "*@*!@[#7]",
                smiles: "Nc1ccccc1",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 135,
                smarts: "[#7]!:*:*",
                smiles: "Nc1ccccc1",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 136,
                smarts: "[#8]=*",
                smiles: "CS(=O)C",
                first_match: &[2, 1],
            },
            MaccsPatternGolden {
                bit: 137,
                smarts: "[!C!cR]",
                smiles: "O1CC1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 138,
                smarts: "[!#6!#1]~[CH2]~*",
                smiles: "CCO",
                first_match: &[2, 1, 0],
            },
            MaccsPatternGolden {
                bit: 139,
                smarts: "[O!H0]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 140,
                smarts: "[#8]",
                smiles: "CCO",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 141,
                smarts: "[CH3]",
                smiles: "CC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 142,
                smarts: "[#7]",
                smiles: "C#N",
                first_match: &[1],
            },
            MaccsPatternGolden {
                bit: 144,
                smarts: "*!:*:*!:*",
                smiles: "Cc1ccccc1C",
                first_match: &[0, 1, 6, 7],
            },
            MaccsPatternGolden {
                bit: 145,
                smarts: "*1~*~*~*~*~*~1",
                smiles: "C1CCCCC1",
                first_match: &[0, 1, 2, 3, 4, 5],
            },
            MaccsPatternGolden {
                bit: 147,
                smarts: "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]",
                smiles: "CCCC",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 148,
                smarts: "*~[!#6!#1](~*)~*",
                smiles: "C[S](C)(C)C",
                first_match: &[0, 1, 2, 3],
            },
            MaccsPatternGolden {
                bit: 149,
                smarts: "[C;H3,H4]",
                smiles: "C",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 150,
                smarts: "*!@*@*!@*",
                smiles: "Cc1ccccc1C",
                first_match: &[0, 1, 6, 7],
            },
            MaccsPatternGolden {
                bit: 151,
                smarts: "[#7!H0]",
                smiles: "CCN",
                first_match: &[2],
            },
            MaccsPatternGolden {
                bit: 152,
                smarts: "[#8]~[#6](~[#6])~[#6]",
                smiles: "CC(C)(C)O",
                first_match: &[4, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 154,
                smarts: "[#6]=[#8]",
                smiles: "CC(=O)O",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 155,
                smarts: "*!@[CH2]!@*",
                smiles: "CCC",
                first_match: &[0, 1, 2],
            },
            MaccsPatternGolden {
                bit: 156,
                smarts: "[#7]~*(~*)~*",
                smiles: "CC(C)N",
                first_match: &[3, 1, 0, 2],
            },
            MaccsPatternGolden {
                bit: 157,
                smarts: "[#6]-[#8]",
                smiles: "CCO",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 158,
                smarts: "[#6]-[#7]",
                smiles: "CCN",
                first_match: &[1, 2],
            },
            MaccsPatternGolden {
                bit: 162,
                smarts: "a",
                smiles: "c1ccccc1",
                first_match: &[0],
            },
            MaccsPatternGolden {
                bit: 165,
                smarts: "[R]",
                smiles: "C1CC1",
                first_match: &[0],
            },
        ]
    }

    #[test]
    fn maccs_patterns_match_rdkit_positive_truth_and_first_atom_maps() {
        let goldens = rdkit_maccs_pattern_positive_goldens();
        assert_eq!(goldens.len(), 136);

        for golden in goldens {
            let mol =
                Molecule::from_smiles_with_sanitize(golden.smiles, true).unwrap_or_else(|error| {
                    panic!("MACCS bit {} target SMILES failed: {error}", golden.bit)
                });
            let query = compile_query_fixture(golden.smarts)
                .unwrap_or_else(|error| panic!("MACCS bit {} SMARTS failed: {error}", golden.bit));
            let matches = get_substruct_matches(&mol, &query);
            assert!(
                !matches.is_empty(),
                "MACCS bit {} should match RDKit positive target {} with SMARTS {}",
                golden.bit,
                golden.smiles,
                golden.smarts
            );
            assert_eq!(
                matches[0].atom_mapping, golden.first_match,
                "MACCS bit {} first RDKit atom map for {}",
                golden.bit, golden.smiles
            );
        }
    }

    #[test]
    fn maccs_bit_030_requires_four_explicit_neighbors_like_rdkit() {
        let mol = Molecule::from_smiles("ON(C)C").expect("fixture should parse");
        let query = compile_query_fixture("[#6]~[!#6!#1](~[#6])(~[#6])~*")
            .expect("MACCS bit 30 SMARTS should build");

        assert_eq!(query.num_atoms(), 5);
        assert_eq!(query.num_bonds(), 4);
        assert!(
            !has_substruct_match(&mol, &query),
            "RDKit does not match MACCS bit 30 against ON(C)C through the max-matches=1 path"
        );
        assert!(
            get_substruct_matches(&mol, &query).is_empty(),
            "RDKit does not match MACCS bit 30 against ON(C)C"
        );
    }

    #[test]
    fn test_atom_matches_basic() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        let q_atom = &mol.atoms()[0];
        let m_atom = &mol.atoms()[1];
        assert!(
            atom_matches(q_atom, &mol, m_atom, &mol),
            "two carbons should match"
        );
    }

    #[test]
    fn test_atom_matches_different_elements() {
        let mut builder = MoleculeBuilder::new();
        let c = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let o = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(c, o, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        // Query atom = C, target atom = O — should not match (C != O, and
        // atomic number check should reject since 6 != 8).
        // But our atom_matches uses query atomic number: query=6, mol=8.
        // If query atomic number != 0, it must match. So C does not match O.
        let c_atom = &mol.atoms()[0]; // C
        let o_atom = &mol.atoms()[1]; // O
        assert!(
            !atom_matches(c_atom, &mol, o_atom, &mol),
            "C should not match O via basic atomic number check"
        );
    }

    #[test]
    fn atom_matches_reproduces_rdkit_plain_atom_defaults() {
        fn one_atom(spec: crate::AtomSpec) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(spec);
            builder.build().expect("one-atom fixture")
        }

        let dummy = one_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let isotope_one = one_atom(crate::AtomSpec::new(crate::Element::DUMMY).with_isotope(1));
        let isotope_two = one_atom(crate::AtomSpec::new(crate::Element::DUMMY).with_isotope(2));
        assert!(atom_matches(
            &dummy.atoms()[0],
            &dummy,
            &isotope_one.atoms()[0],
            &isotope_one,
        ));
        assert!(atom_matches(
            &isotope_one.atoms()[0],
            &isotope_one,
            &dummy.atoms()[0],
            &dummy,
        ));
        assert!(!atom_matches(
            &isotope_one.atoms()[0],
            &isotope_one,
            &isotope_two.atoms()[0],
            &isotope_two,
        ));

        let neutral_carbon = one_atom(crate::AtomSpec::new(crate::Element::C));
        let charged_carbon =
            one_atom(crate::AtomSpec::new(crate::Element::C).with_formal_charge(1));
        assert!(atom_matches(
            &neutral_carbon.atoms()[0],
            &neutral_carbon,
            &charged_carbon.atoms()[0],
            &charged_carbon,
        ));
        assert!(!atom_matches(
            &charged_carbon.atoms()[0],
            &charged_carbon,
            &neutral_carbon.atoms()[0],
            &neutral_carbon,
        ));

        let radical_carbon =
            one_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(1));
        assert!(!atom_matches(
            &radical_carbon.atoms()[0],
            &radical_carbon,
            &neutral_carbon.atoms()[0],
            &neutral_carbon,
        ));
    }

    #[test]
    fn test_bond_matches_single() {
        let mut builder = MoleculeBuilder::new();
        let c0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let c1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(c0, c1, BondOrder::Single))
            .expect("add bond");
        let mol = builder.build().expect("build");
        assert!(bond_compat(
            &mol.bonds()[0],
            &mol,
            &mol.bonds()[0],
            &mol,
            &SubstructMatchParams::default(),
            &build_query_match_context(&mol),
        ));
    }

    #[test]
    fn test_vf2_graph_building() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        assert_eq!(g.n_atoms, 2);
        assert_eq!(g.n_bonds, 1);
        assert_eq!(g.out_degree(0), 1);
        assert_eq!(g.out_degree(1), 1);
        assert_eq!(g.out_edges(0)[0].0, 1);
        assert_eq!(g.out_edges(1)[0].0, 0);
    }

    #[test]
    fn smarts_vf2_other_index() {
        let cc = make_mol_cc();
        let graph = build_vf2_graph(&cc);
        assert_eq!(get_other_idx(&graph, 0, 0), 1);
        assert_eq!(get_other_idx(&graph, 0, 1), 0);

        // RDKit returns source(edge) when `vertex` is not the source; it does
        // not validate that `vertex` is an endpoint.
        assert_eq!(get_other_idx(&graph, 0, 99), 0);
    }

    #[test]
    fn test_sort_nodes_by_frequency_small() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let order = sort_nodes_by_frequency(&g);
        // Both nodes have degree 1, so order depends on sort stability.
        assert_eq!(order.len(), 2);
        // Both should be present.
        assert!(order.contains(&0));
        assert!(order.contains(&1));
    }

    #[test]
    fn smarts_vf2_node_order() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let leaf1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let leaf2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let leaf3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for leaf in [leaf1, leaf2, leaf3] {
            builder
                .add_bond(crate::BondSpec::new(center, leaf, BondOrder::Single))
                .expect("add star bond");
        }
        let molecule = builder.build().expect("build star and isolated atom");

        let order = sort_nodes_by_frequency(&build_vf2_graph(&molecule));
        assert_eq!(order[0], center.index());
        assert_eq!(order[4], isolated.index());
        let mut middle = order[1..4].to_vec();
        middle.sort_unstable();
        assert_eq!(middle, vec![leaf1.index(), leaf2.index(), leaf3.index()]);
    }

    #[test]
    fn smarts_vf2_node_compare_degree() {
        let lower_out = NodeInfo {
            id: 0,
            in_deg: 9,
            out_deg: 1,
        };
        let higher_out = NodeInfo {
            id: 1,
            in_deg: 0,
            out_deg: 2,
        };
        assert_eq!(
            node_info_cmp1(&lower_out, &higher_out),
            std::cmp::Ordering::Less
        );

        let lower_in = NodeInfo {
            id: 2,
            in_deg: 1,
            out_deg: 3,
        };
        let higher_in = NodeInfo {
            id: 3,
            in_deg: 2,
            out_deg: 3,
        };
        assert_eq!(
            node_info_cmp1(&lower_in, &higher_in),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            node_info_cmp1(&lower_in, &NodeInfo { id: 4, ..lower_in }),
            std::cmp::Ordering::Equal
        );
    }

    #[test]
    fn smarts_vf2_node_compare_frequency() {
        let isolated = NodeInfo {
            id: 0,
            in_deg: 0,
            out_deg: 0,
        };
        let connected = NodeInfo {
            id: 1,
            in_deg: 2,
            out_deg: 9,
        };
        assert_eq!(
            node_info_cmp2(&isolated, &connected),
            std::cmp::Ordering::Greater
        );
        assert_eq!(
            node_info_cmp2(&connected, &isolated),
            std::cmp::Ordering::Less
        );

        let rarer = NodeInfo {
            id: 2,
            in_deg: 8,
            out_deg: 1,
        };
        let common = NodeInfo {
            id: 3,
            in_deg: 1,
            out_deg: 2,
        };
        assert_eq!(node_info_cmp2(&rarer, &common), std::cmp::Ordering::Less);

        let lower_valence = NodeInfo {
            id: 4,
            in_deg: 2,
            out_deg: 3,
        };
        let higher_valence = NodeInfo {
            id: 5,
            in_deg: 4,
            out_deg: 3,
        };
        assert_eq!(
            node_info_cmp2(&lower_valence, &higher_valence),
            std::cmp::Ordering::Less
        );
        assert_eq!(
            node_info_cmp2(
                &lower_valence,
                &NodeInfo {
                    id: 6,
                    ..lower_valence
                }
            ),
            std::cmp::Ordering::Equal
        );
    }

    #[test]
    fn test_vf2_state_initial() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let state = Vf2SubState::new(&g, &g, true);
        assert!(!state.is_goal());
        assert!(!state.is_dead());
        assert_eq!(state.core_len, 0);
    }

    #[test]
    fn smarts_vf2_state_new() {
        let query = build_vf2_graph(&make_mol_cc());
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let target = build_vf2_graph(&target_molecule);

        let unsorted = Vf2SubState::new(&query, &target, false);
        assert_eq!((unsorted.n1, unsorted.n2), (2, 3));
        assert_eq!(
            (unsorted.core_len, unsorted.t1_len, unsorted.t2_len),
            (0, 0, 0)
        );
        assert_eq!(unsorted.core_1, vec![NULL_NODE; 2]);
        assert_eq!(unsorted.core_2, vec![NULL_NODE; 3]);
        assert_eq!(unsorted.term_1, vec![0; 2]);
        assert_eq!(unsorted.term_2, vec![0; 3]);
        assert!(unsorted.debug_order().is_none());

        let sorted = Vf2SubState::new(&query, &target, true);
        assert_eq!(
            sorted.debug_order(),
            Some(sort_nodes_by_frequency(&query).as_slice())
        );
    }

    #[test]
    fn smarts_vf2_state_clone() {
        let molecule = make_mol_cc();
        let graph = build_vf2_graph(&molecule);
        let mut original = Vf2SubState::new(&graph, &graph, true);
        original.add_pair(0, 0);

        let mut cloned = original.clone_state();
        assert_eq!(cloned.core_len, original.core_len);
        assert_eq!(cloned.core_1, original.core_1);
        assert_eq!(cloned.core_2, original.core_2);
        assert_eq!(cloned.term_1, original.term_1);
        assert_eq!(cloned.term_2, original.term_2);
        assert_eq!(cloned.order, original.order);

        cloned.add_pair(1, 1);
        assert_eq!(cloned.core_len, 2);
        assert_eq!(original.core_len, 1);
        assert_eq!(original.core_1[1], NULL_NODE);
    }

    #[test]
    fn smarts_vf2_goal() {
        let molecule = make_mol_cc();
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, false);
        assert!(!state.is_goal());
        state.add_pair(0, 0);
        assert!(!state.is_goal());
        state.add_pair(1, 1);
        assert!(state.is_goal());
    }

    #[test]
    fn smarts_vf2_match_checks() {
        let molecule = make_mol_cc();
        let graph = build_vf2_graph(&molecule);
        let state = Vf2SubState::new(&graph, &graph, false);
        let mut seen = None;
        let mut check = |c1: &[NodeId], c2: &[NodeId]| {
            seen = Some((c1.to_vec(), c2.to_vec()));
            c1 == [0, 1] && c2 == [1, 0]
        };
        assert!(state.match_checks(&[0, 1], &[1, 0], &mut check));
        assert_eq!(seen, Some((vec![0, 1], vec![1, 0])));
    }

    #[test]
    fn smarts_vf2_dead() {
        let query_molecule = Molecule::from_smiles("CCC").expect("parse query");
        let target_molecule = make_mol_cc();
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        assert!(Vf2SubState::new(&query, &target, false).is_dead());

        let mut terminal_dead = Vf2SubState::new(&target, &query, false);
        terminal_dead.t1_len = 2;
        terminal_dead.t2_len = 1;
        assert!(terminal_dead.is_dead());
        terminal_dead.t2_len = 2;
        assert!(!terminal_dead.is_dead());
    }

    #[test]
    fn smarts_vf2_core_len() {
        let molecule = make_mol_cc();
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, false);
        assert_eq!(state.core_len(), 0);
        state.add_pair(0, 0);
        assert_eq!(state.core_len(), 1);
        state.add_pair(1, 1);
        assert_eq!(state.core_len(), 2);
    }

    #[test]
    fn test_vf2_next_pair_initial() {
        let cc = make_mol_cc();
        let g = build_vf2_graph(&cc);
        let state = Vf2SubState::new(&g, &g, true);
        let mut pair = Vf2Pair::new();
        let has_next = state.next_pair(&mut pair);
        assert!(has_next, "should find a next pair");
    }

    #[test]
    fn smarts_vf2_next_pair() {
        let mut builder = MoleculeBuilder::new();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let leaf1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let leaf2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for leaf in [leaf1, leaf2] {
            builder
                .add_bond(crate::BondSpec::new(center, leaf, BondOrder::Single))
                .expect("add star bond");
        }
        let molecule = builder.build().expect("build star");
        let graph = build_vf2_graph(&molecule);

        let state = Vf2SubState::new(&graph, &graph, true);
        let mut pair = Vf2Pair::new();
        let mut initial_pairs = Vec::new();
        while state.next_pair(&mut pair) {
            initial_pairs.push((pair.n1, pair.n2));
        }
        assert_eq!(
            initial_pairs,
            vec![
                (center.index(), 0),
                (center.index(), 1),
                (center.index(), 2)
            ]
        );

        let mut terminal_state = Vf2SubState::new(&graph, &graph, true);
        terminal_state.add_pair(center.index(), center.index());
        let mut terminal_pair = Vf2Pair::new();
        let mut target_neighbors = Vec::new();
        while terminal_state.next_pair(&mut terminal_pair) {
            assert!(terminal_pair.n1 == leaf1.index() || terminal_pair.n1 == leaf2.index());
            target_neighbors.push(terminal_pair.n2);
        }
        target_neighbors.sort_unstable();
        assert_eq!(target_neighbors, vec![leaf1.index(), leaf2.index()]);
    }

    #[test]
    fn smarts_vf2_feasible_pair() {
        let query_molecule = make_mol_cc();
        let query = build_vf2_graph(&query_molecule);

        let single_atom = Molecule::from_smiles("C").expect("parse single atom");
        let single = build_vf2_graph(&single_atom);
        let state = Vf2SubState::new(&query, &single, false);
        assert!(!state.is_feasible_pair(0, 0, &|_, _| true, &|_, _| true));

        let target_molecule = Molecule::from_smiles("CC").expect("parse target");
        let target = build_vf2_graph(&target_molecule);
        let state = Vf2SubState::new(&query, &target, false);
        assert!(!state.is_feasible_pair(0, 0, &|_, _| false, &|_, _| true));
        assert!(state.is_feasible_pair(0, 0, &|_, _| true, &|_, _| true));

        let mut mapped = Vf2SubState::new(&query, &target, false);
        mapped.add_pair(0, 0);
        assert!(!mapped.is_feasible_pair(1, 1, &|_, _| true, &|_, _| false));
        assert!(mapped.is_feasible_pair(1, 1, &|_, _| true, &|_, _| true));

        let disconnected_molecule = Molecule::from_smiles("CC.CC").expect("parse fragments");
        let disconnected = build_vf2_graph(&disconnected_molecule);
        let mut missing_edge = Vf2SubState::new(&query, &disconnected, false);
        missing_edge.add_pair(0, 0);
        assert!(!missing_edge.is_feasible_pair(1, 2, &|_, _| true, &|_, _| true));
    }

    #[test]
    fn smarts_vf2_add_pair() {
        let molecule = Molecule::from_smiles("CCC").expect("parse chain");
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, false);

        state.add_pair(1, 1);
        assert_eq!(state.core_len, 1);
        assert_eq!(state.core_1, vec![NULL_NODE, 1, NULL_NODE]);
        assert_eq!(state.core_2, vec![NULL_NODE, 1, NULL_NODE]);
        assert_eq!(state.term_1, vec![1, 1, 1]);
        assert_eq!(state.term_2, vec![1, 1, 1]);
        assert_eq!((state.t1_len, state.t2_len), (3, 3));

        state.add_pair(0, 0);
        assert_eq!(state.core_len, 2);
        assert_eq!(state.core_1, vec![0, 1, NULL_NODE]);
        assert_eq!(state.core_2, vec![0, 1, NULL_NODE]);
        assert_eq!(state.term_1, vec![1, 1, 1]);
        assert_eq!(state.term_2, vec![1, 1, 1]);
        assert_eq!((state.t1_len, state.t2_len), (3, 3));
    }

    #[test]
    fn smarts_vf2_core_set() {
        let molecule = Molecule::from_smiles("CCC").expect("parse chain");
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, false);
        state.add_pair(2, 0);
        state.add_pair(0, 2);

        let (c1, c2) = state.get_core_set();
        assert_eq!(c1, vec![0, 2]);
        assert_eq!(c2, vec![2, 0]);
    }

    #[test]
    fn smarts_vf2_clone() {
        let molecule = Molecule::from_smiles("CCC").expect("parse chain");
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, true);
        state.add_pair(1, 1);
        let mut cloned = state.clone();

        cloned.back_track(1, 1);
        assert_eq!(cloned.core_len, 0);
        assert_eq!(state.core_len, 1);
        assert_eq!(state.core_1[1], 1);
        assert_eq!(state.core_2[1], 1);
    }

    #[test]
    fn smarts_vf2_backtrack() {
        let molecule = Molecule::from_smiles("CCC").expect("parse chain");
        let graph = build_vf2_graph(&molecule);
        let mut state = Vf2SubState::new(&graph, &graph, false);
        state.add_pair(1, 1);
        state.add_pair(0, 0);
        assert_eq!(state.core_len, 2);

        state.back_track(0, 0);
        assert_eq!(state.core_len, 1);
        assert_eq!(state.core_1, vec![NULL_NODE, 1, NULL_NODE]);
        assert_eq!(state.core_2, vec![NULL_NODE, 1, NULL_NODE]);
        assert_eq!(state.term_1, vec![1, 1, 1]);
        assert_eq!(state.term_2, vec![1, 1, 1]);
        assert_eq!((state.t1_len, state.t2_len), (3, 3));

        state.back_track(1, 1);
        assert_eq!(state.core_len, 0);
        assert_eq!(state.term_1, vec![0, 0, 0]);
        assert_eq!(state.term_2, vec![0, 0, 0]);
        assert_eq!((state.t1_len, state.t2_len), (0, 0));
    }

    #[test]
    fn smarts_vf2_match_one() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        let mut state = Vf2SubState::new(&query, &target, true);
        let mut checks = 0;
        let mut accept_second = |_: &[NodeId], _: &[NodeId]| {
            checks += 1;
            checks == 2
        };
        let (c1, c2) = state
            .match_one(&|_, _| true, &|_, _| true, Some(&mut accept_second))
            .expect("second complete mapping accepted");
        assert_eq!(c1, vec![0, 1]);
        assert_eq!(c2.len(), 2);
        assert_eq!(checks, 2);

        let mut dead = Vf2SubState::new(&target, &query, true);
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(
            dead.match_one(&|_, _| true, &|_, _| true, Some(&mut accept_all))
                .is_none()
        );
    }

    #[test]
    fn smarts_vf2_match_all() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);

        let mut state = Vf2SubState::new(&query, &target, true);
        let mut results = Vec::new();
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(!state.match_all(
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut results,
            0,
        ));
        assert_eq!(results.len(), 4);

        let mut limited_state = Vf2SubState::new(&query, &target, true);
        let mut limited = Vec::new();
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(limited_state.match_all(
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut limited,
            2,
        ));
        assert_eq!(limited.len(), 2);
    }

    #[test]
    fn smarts_vf2_free_match_one() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        let mut state = Vf2SubState::new(&query, &target, false);
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        let (c1, c2) = vf2_match(
            &mut state,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
        )
        .expect("free match finds mapping");
        assert_eq!(c1.len(), state.core_len());
        assert_eq!(c1, vec![0, 1]);
        assert_eq!(c2.len(), 2);

        let mut dead = Vf2SubState::new(&target, &query, false);
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(vf2_match(&mut dead, &|_, _| true, &|_, _| true, Some(&mut accept_all),).is_none());
    }

    #[test]
    fn smarts_vf2_free_match_all() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        let mut state = Vf2SubState::new(&query, &target, false);
        let mut results = Vec::new();
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(vf2_match_all(
            &mut state,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut results,
            2,
        ));
        assert_eq!(results.len(), 2);

        let mut dead = Vf2SubState::new(&target, &query, false);
        let mut no_results = Vec::new();
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(!vf2_match_all(
            &mut dead,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut no_results,
            2,
        ));
        assert!(no_results.is_empty());
    }

    #[test]
    fn smarts_vf2_entry_one() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        let mut result = vec![(99, 99)];
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(vf2_entry_one(
            &query,
            &target,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut result,
        ));
        assert_eq!(result.len(), 2);
        assert!(!result.contains(&(99, 99)));

        let mut no_result = vec![(99, 99)];
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(!vf2_entry_one(
            &target,
            &query,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut no_result,
        ));
        assert!(no_result.is_empty());
    }

    #[test]
    fn smarts_vf2_entry_all() {
        let query_molecule = make_mol_cc();
        let target_molecule = Molecule::from_smiles("CCC").expect("parse target");
        let query = build_vf2_graph(&query_molecule);
        let target = build_vf2_graph(&target_molecule);
        let mut results = vec![(vec![99], vec![99])];
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(vf2_entry_all(
            &query,
            &target,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut results,
            3,
        ));
        assert_eq!(results.len(), 3);
        assert!(!results.iter().any(|mapping| mapping.0 == [99]));

        let mut no_results = vec![(vec![99], vec![99])];
        let mut accept_all = |_: &[NodeId], _: &[NodeId]| true;
        assert!(!vf2_entry_all(
            &target,
            &query,
            &|_, _| true,
            &|_, _| true,
            Some(&mut accept_all),
            &mut no_results,
            3,
        ));
        assert!(no_results.is_empty());
    }

    #[test]
    fn smarts_match_chiral_label() {
        let atom = |tag| {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_chiral_tag(tag));
            let molecule = builder.build().expect("build one atom");
            molecule.atoms()[0].clone()
        };
        assert!(has_chiral_label(&atom(ChiralTag::TetrahedralCw)));
        assert!(has_chiral_label(&atom(ChiralTag::TetrahedralCcw)));
        for tag in [
            ChiralTag::Unspecified,
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ] {
            assert!(
                !has_chiral_label(&atom(tag)),
                "unexpected label for {tag:?}"
            );
        }
    }

    #[test]
    fn smarts_match_enhanced_stereo() {
        fn grouped(groups: &[(StereoGroupKind, &[usize])]) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            let atoms: Vec<_> = (0..3)
                .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::C)))
                .collect();
            for (kind, members) in groups {
                builder
                    .add_stereo_group(crate::StereoGroup::new(
                        *kind,
                        members.iter().map(|&idx| atoms[idx]).collect(),
                        Vec::new(),
                    ))
                    .expect("add stereo group");
            }
            builder.build().expect("build grouped molecule")
        }

        let query_or = grouped(&[(StereoGroupKind::Or, &[0, 1])]);
        let mol_or = grouped(&[(StereoGroupKind::Or, &[0, 1])]);
        let mol_and = grouped(&[(StereoGroupKind::And, &[0, 1])]);
        let identity = [0, 1, 2];
        let or_membership = [Some(0), Some(0), None];
        let same = [Some(true), Some(true), None];
        assert!(enhanced_stereo_is_ok(
            &mol_or,
            &query_or,
            &identity,
            &or_membership,
            &same,
        ));
        assert!(enhanced_stereo_is_ok(
            &mol_and,
            &query_or,
            &identity,
            &or_membership,
            &same,
        ));

        let query_and = grouped(&[(StereoGroupKind::And, &[0, 1])]);
        assert!(!enhanced_stereo_is_ok(
            &mol_or,
            &query_and,
            &identity,
            &or_membership,
            &same,
        ));
        assert!(enhanced_stereo_is_ok(
            &mol_and,
            &query_and,
            &identity,
            &or_membership,
            &same,
        ));

        let absolute = grouped(&[]);
        assert!(!enhanced_stereo_is_ok(
            &absolute,
            &query_or,
            &identity,
            &[None, None, None],
            &[None, None, None],
        ));
        assert!(!enhanced_stereo_is_ok(
            &mol_or,
            &query_or,
            &identity,
            &or_membership,
            &[Some(true), Some(false), None],
        ));

        let split_query = grouped(&[(StereoGroupKind::Or, &[0]), (StereoGroupKind::Or, &[1])]);
        assert!(!enhanced_stereo_is_ok(
            &mol_or,
            &split_query,
            &identity,
            &or_membership,
            &same,
        ));
    }

    #[test]
    fn smarts_match_insert_unique() {
        let mut matches = BTreeSet::new();
        let first = vec![(0, 2), (1, 1)];
        assert!(insert_if_needed(&mut matches, first.clone()));
        assert!(!insert_if_needed(&mut matches, first.clone()));

        let lexicographically_smaller = vec![(0, 1), (1, 2)];
        assert!(insert_if_needed(
            &mut matches,
            lexicographically_smaller.clone()
        ));
        assert_eq!(matches, BTreeSet::from([lexicographically_smaller]));

        let distinct_atom_set = vec![(0, 3), (1, 2)];
        assert!(insert_if_needed(&mut matches, distinct_atom_set.clone()));
        assert!(matches.contains(&distinct_atom_set));
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn smarts_match_try_insert() {
        let first = vec![(0, 2), (1, 1)];
        let duplicate_atom_set = vec![(0, 1), (1, 2)];

        let mut params = SubstructMatchParams {
            max_matches: 1,
            uniquify: false,
            ..SubstructMatchParams::default()
        };
        let mut matches = BTreeSet::new();
        assert!(try_to_insert(&mut matches, first.clone(), &params));
        assert!(!try_to_insert(
            &mut matches,
            duplicate_atom_set.clone(),
            &params
        ));
        assert_eq!(matches, BTreeSet::from([first.clone()]));

        params.max_matches = 10;
        params.uniquify = true;
        assert!(try_to_insert(
            &mut matches,
            duplicate_atom_set.clone(),
            &params
        ));
        assert_eq!(matches, BTreeSet::from([duplicate_atom_set]));

        params.uniquify = false;
        assert!(try_to_insert(&mut matches, first.clone(), &params));
        assert!(try_to_insert(&mut matches, first, &params));
        assert_eq!(matches.len(), 2);
    }

    #[test]
    fn smarts_match_final_check_setup() {
        let mut builder = MoleculeBuilder::new();
        let atoms: Vec<_> = (0..4)
            .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::C)))
            .collect();
        builder
            .add_stereo_group(crate::StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![atoms[0]],
                Vec::new(),
            ))
            .expect("add absolute group");
        builder
            .add_stereo_group(crate::StereoGroup::new(
                StereoGroupKind::Or,
                vec![atoms[1], atoms[2]],
                Vec::new(),
            ))
            .expect("add OR group");
        builder
            .add_stereo_group(crate::StereoGroup::new(
                StereoGroupKind::And,
                vec![atoms[3]],
                Vec::new(),
            ))
            .expect("add AND group");
        let molecule = builder.build().expect("build grouped molecule");

        let disabled =
            MolMatchFinalCheckSetup::new(&molecule, &molecule, &SubstructMatchParams::default());
        assert_eq!(disabled.mol_stereo_groups, vec![None; 4]);

        let enabled = MolMatchFinalCheckSetup::new(
            &molecule,
            &molecule,
            &SubstructMatchParams {
                use_enhanced_stereo: true,
                ..SubstructMatchParams::default()
            },
        );
        assert_eq!(
            enabled.mol_stereo_groups,
            vec![None, Some(1), Some(1), Some(2)]
        );
    }

    #[test]
    fn smarts_match_final_check() {
        let query = make_mol_cc();
        let molecule = Molecule::from_smiles("CCC").expect("parse target");
        let mapping = ([0, 1], [0, 1]);

        let mut params = SubstructMatchParams::default();
        let setup = MolMatchFinalCheckSetup::new(&query, &molecule, &params);
        let mut seen = Vec::new();
        assert!(
            rdkit_match_final_check(
                &molecule, &query, &params, &mapping.0, &mapping.1, &setup, &mut seen,
            )
            .expect("final check")
        );
        assert!(
            !rdkit_match_final_check(
                &molecule, &query, &params, &mapping.0, &mapping.1, &setup, &mut seen,
            )
            .expect("duplicate final check")
        );

        params.uniquify = false;
        params.extra_final_check = Some(Arc::new(|_, atom_ids| atom_ids == [0, 1]));
        let setup = MolMatchFinalCheckSetup::new(&query, &molecule, &params);
        assert!(
            rdkit_match_final_check(
                &molecule,
                &query,
                &params,
                &mapping.0,
                &mapping.1,
                &setup,
                &mut Vec::new(),
            )
            .expect("accepted callback")
        );
        params.extra_final_check = Some(Arc::new(|_, _| false));
        assert!(
            !rdkit_match_final_check(
                &molecule,
                &query,
                &params,
                &mapping.0,
                &mapping.1,
                &setup,
                &mut Vec::new(),
            )
            .expect("rejected callback")
        );

        params.extra_final_check = None;
        params.use_generic_matchers = true;
        assert!(
            rdkit_match_final_check(
                &molecule,
                &query,
                &params,
                &mapping.0,
                &mapping.1,
                &setup,
                &mut Vec::new(),
            )
            .expect("generic matcher without labels accepts the match")
        );
    }

    #[test]
    fn smarts_match_atom_label() {
        fn one_atom(element: crate::Element, tag: ChiralTag) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            builder.add_atom(crate::AtomSpec::new(element).with_chiral_tag(tag));
            builder.build().expect("build one atom")
        }

        let query = one_atom(crate::Element::C, ChiralTag::TetrahedralCw);
        let unspecified = one_atom(crate::Element::C, ChiralTag::Unspecified);
        let specified = one_atom(crate::Element::C, ChiralTag::TetrahedralCcw);
        let oxygen = one_atom(crate::Element::O, ChiralTag::TetrahedralCcw);
        let context = build_query_match_context(&unspecified);
        let mut params = SubstructMatchParams {
            use_chirality: true,
            ..SubstructMatchParams::default()
        };
        assert!(!atom_label_matches(
            &query,
            &unspecified,
            0,
            0,
            &params,
            None,
            &context,
        ));
        params.specified_stereo_query_matches_unspecified = true;
        assert!(atom_label_matches(
            &query,
            &unspecified,
            0,
            0,
            &params,
            None,
            &context,
        ));
        let context = build_query_match_context(&specified);
        assert!(atom_label_matches(
            &query, &specified, 0, 0, &params, None, &context,
        ));
        let context = build_query_match_context(&oxygen);
        assert!(!atom_label_matches(
            &query, &oxygen, 0, 0, &params, None, &context,
        ));
    }

    #[test]
    fn smarts_match_bond_label() {
        fn double_bond(stereo: BondStereo) -> Molecule {
            let mut builder = MoleculeBuilder::new();
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Double).with_stereo(stereo))
                .expect("add double bond");
            builder.build().expect("build double bond")
        }

        let query = double_bond(BondStereo::E);
        let unspecified = double_bond(BondStereo::None);
        let specified = double_bond(BondStereo::Z);
        let mut params = SubstructMatchParams {
            use_chirality: true,
            ..SubstructMatchParams::default()
        };
        assert!(!bond_label_matches(
            &query,
            &unspecified,
            0,
            0,
            &params,
            &build_query_match_context(&unspecified),
        ));
        params.specified_stereo_query_matches_unspecified = true;
        assert!(bond_label_matches(
            &query,
            &unspecified,
            0,
            0,
            &params,
            &build_query_match_context(&unspecified),
        ));
        assert!(bond_label_matches(
            &query,
            &specified,
            0,
            0,
            &params,
            &build_query_match_context(&specified),
        ));

        let single = Molecule::from_smiles("CC").expect("parse single bond");
        assert!(!bond_label_matches(
            &query,
            &single,
            0,
            0,
            &params,
            &build_query_match_context(&single),
        ));
    }
}
