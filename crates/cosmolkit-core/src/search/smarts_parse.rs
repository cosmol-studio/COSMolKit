//! SMARTS parser — recursive-descent parser producing query-predicate trees.
//!
//! ## RDKit provenance (protocol: dev/source_reproduction_protocol.md)
//!
//! The SMARTS parser corresponds to RDKit's `GraphMol/SmilesParse/SmilesParse.cpp`
//! (MolFromSmarts entry point, labelRecursivePatterns helper) and the bison/flex
//! grammars `smarts.yy` / `smarts.ll`.  The flex/bison lexer and parser are
//! replaced here by a hand-written recursive-descent parser that produces the
//! same semantic query trees.
//!
//! C++ source lines are copied verbatim as commented blocks with two-axis
//! RDKit status markers per `dev/source_reproduction_protocol.md`:
//!   // RDKit✔️✔️: <C++ line>   — fully ported, behaviour identical
//!   // RDKit❗✔️: <C++ line>   — adapted for Rust / COSMolKit differences
//!
//! ## Design
//!
//! This module is the sole canonical SMARTS parser/compiler owner.
//! The shared typed query graph remains owned by `search::query`; callers must
//! migrate into this module and that model instead of adding another parser,
//! compatibility graph, or consumer-local decoder. The parser tokenizes the
//! SMARTS string, recursively produces a private graph, and immediately
//! materializes the sole query-bearing `Molecule` representation.
//!
//! The parser's private graph is materialized immediately as the canonical
//! query-bearing `Molecule` consumed by the substructure matcher.

use std::collections::BTreeMap;

use crate::search::query::{
    AtomRangeDataFunction, CompositeQueryType, make_atom_null_query,
    make_atom_possible_range_query, make_bond_is_in_ring_query, make_bond_null_query,
    make_bond_order_equals_query, query_bond_expand_query,
};
use crate::{
    AtomQueryPredicate, AtomSpec, BondDirection, BondOrder, BondQueryPredicate, BondSpec,
    ChiralTag, Element, Molecule, MoleculeBuilder, QueryNode, SmartsParseError,
};

// ---------------------------------------------------------------------------
// ParsedSmartsGraph - private parser output
// ---------------------------------------------------------------------------

/// The result of parsing a SMARTS string.
///
/// RDKit✔️✔️: RDKit returns an RWMol with QueryAtom / QueryBond objects.
/// COSMolKit returns a separate struct of query-predicate trees paired via
/// indexing (`atom_queries[i]` is the i-th atom, `bond_queries[i]` is the bond
/// between atoms i and i+1 in SMARTS order).
#[derive(Debug, Clone)]
struct ParsedSmartsGraph {
    /// Query trees for each atom in the pattern.
    pub atom_queries: Vec<QueryNode<AtomQueryPredicate>>,
    /// Atom-map properties aligned with `atom_queries`.
    pub atom_maps: Vec<Option<u32>>,
    /// Query trees for each bond in the pattern (length = atom_queries.len() - 1).
    pub bond_queries: Vec<QueryNode<BondQueryPredicate>>,
    /// Directional state aligned with `bond_queries`.
    pub bond_directions: Vec<BondDirection>,
    /// Query bond endpoints in SMARTS atom-index space.
    pub bond_edges: Vec<(usize, usize)>,
    /// Ring-closure specifications: (closure_number, atom_index_in_pattern)
    pub ring_closures: Vec<(u32, usize)>,
    /// Ring-closure bond query specifications: (closure_number, atom_index, bond_query).
    pub ring_closure_bonds: Vec<(u32, usize, QueryNode<BondQueryPredicate>)>,
}

impl ParsedSmartsGraph {
    #[must_use]
    pub fn new(
        atom_queries: Vec<QueryNode<AtomQueryPredicate>>,
        atom_maps: Vec<Option<u32>>,
        bond_queries: Vec<QueryNode<BondQueryPredicate>>,
        bond_directions: Vec<BondDirection>,
        bond_edges: Vec<(usize, usize)>,
        ring_closures: Vec<(u32, usize)>,
        ring_closure_bonds: Vec<(u32, usize, QueryNode<BondQueryPredicate>)>,
    ) -> Self {
        Self {
            atom_queries,
            atom_maps,
            bond_queries,
            bond_directions,
            bond_edges,
            ring_closures,
            ring_closure_bonds,
        }
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atom_queries.len()
    }

    #[must_use]
    pub fn atom_query(&self, idx: usize) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.atom_queries.get(idx)
    }

    #[must_use]
    pub fn bond_query(&self, idx: usize) -> Option<&QueryNode<BondQueryPredicate>> {
        self.bond_queries.get(idx)
    }
}

fn to_mol(inp: &str) -> Result<Molecule, String> {
    // BEGIN RDKIT CPP FUNCTION toMol
    // RDKit✔️✔️: std::unique_ptr<RWMol> toMol(const std::string &inp,
    // RDKit✔️✔️:                              int func(const std::string &,
    // RDKit✔️✔️:                                       std::vector<RDKit::RWMol *> &),
    // RDKit✔️✔️:                              const std::string &origInp) {
    // RDKit✔️✔️:   // empty strings produce empty molecules:
    // RDKit✔️✔️:   if (inp.empty()) {
    // RDKit✔️✔️:     return std::make_unique<RWMol>();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::unique_ptr<RWMol> res;
    // RDKit✔️✔️:   std::vector<RDKit::RWMol *> molVect;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     func(inp, molVect);
    // RDKit✔️✔️:     if (!molVect.empty()) {
    // RDKit✔️✔️:       res.reset(molVect[0]);
    // RDKit✔️✔️:       SmilesParseOps::CloseMolRings(res.get(), false);
    // RDKit✔️✔️:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
    // RDKit✔️✔️:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
    // RDKit✔️✔️:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
    // RDKit✔️✔️:       // No sense leaving this bookmark intact:
    // RDKit✔️✔️:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
    // RDKit✔️✔️:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
    // RDKit✔️✔️:                              // occur last in this if.
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } catch (SmilesParseException &e) {
    // RDKit✔️✔️:     std::string nm = "SMILES";
    // RDKit✔️✔️:     if (func == smarts_parse) {
    // RDKit✔️✔️:       nm = "SMARTS";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
    // RDKit✔️✔️:                           << " for input: '" << origInp << "'" << std::endl;
    // RDKit✔️✔️:
    // RDKit✔️✔️:     // reset res so that we return a nullptr. We don't want to reset(),
    // RDKit✔️✔️:     // because that would delete the mol and leak any unmatched
    // RDKit✔️✔️:     // ring closure bonds. These will be cleaned up in the loop below.
    // RDKit✔️✔️:     res.release();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto *molPtr : molVect) {
    // RDKit✔️✔️:     if (molPtr) {
    // RDKit✔️✔️:       // Clean-up the bond bookmarks when not calling CloseMolRings
    // RDKit✔️✔️:       SmilesParseOps::CleanupAfterParseError(molPtr);
    // RDKit✔️✔️:       delete molPtr;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION toMol
    // Local complexity review: the canonical parser is linear in the SMARTS
    // input, and materialization makes one pass over atoms and one over bonds.
    // Rust ownership replaces RDKit's raw-pointer cleanup without an extra
    // graph scan, temporary molecule list, matcher-time decode, or clone.
    if inp.is_empty() {
        return Ok(Molecule::new());
    }

    let parsed = parse_smarts(inp).map_err(|error| error.to_string())?;
    let mut builder = MoleculeBuilder::new();
    let atom_ids: Vec<_> = parsed
        .atom_queries
        .iter()
        .zip(parsed.atom_maps.iter().copied())
        .map(|(query, atom_map)| {
            let (query, chiral_tag, chiral_permutation) =
                materialize_smarts_atom_state(query.clone()).map_err(|error| error.to_string())?;
            let isotope = atom_isotope(&query);
            let (atomic_number, aromatic) = source_query_atom_identity(&query);
            let mut spec = AtomSpec::new(
                Element::from_atomic_number(atomic_number)
                    .expect("source SMARTS atomic number is a valid element"),
            )
            .with_aromatic(aromatic)
            .with_query(query)
            .with_chiral_tag(chiral_tag);
            if let Some(chiral_permutation) = chiral_permutation {
                spec = spec.with_chiral_permutation(chiral_permutation);
            }
            if let Some(isotope) = isotope {
                spec = spec.with_isotope(isotope);
            }
            if let Some(atom_map) = atom_map {
                spec = spec.with_atom_map(atom_map);
            }
            Ok(builder.add_atom(spec))
        })
        .collect::<Result<_, String>>()?;
    for (bond_index, (((begin_atom_index, end_atom_index), query), direction)) in parsed
        .bond_edges
        .iter()
        .copied()
        .zip(parsed.bond_queries.iter())
        .zip(parsed.bond_directions.iter().copied())
        .enumerate()
    {
        // RDKit's SMARTS grammar assigns this parser-order index to every
        // materialized bond. CX coordinate-bond extensions consume it before
        // CleanupAfterParsing removes the temporary property.
        builder
            .add_bond(
                BondSpec::new(
                    atom_ids[begin_atom_index],
                    atom_ids[end_atom_index],
                    representative_bond_order(query),
                )
                .with_direction(direction)
                .with_prop(
                    crate::notation::smiles::CXSMILES_BOND_IDX_PROP,
                    bond_index.to_string(),
                )
                .with_query(query.clone()),
            )
            .map_err(|error| error.to_string())?;
    }
    builder.build().map_err(|error| error.to_string())
}

fn source_query_atom_identity(query: &QueryNode<AtomQueryPredicate>) -> (u8, bool) {
    // BEGIN RDKIT CPP BLOCK SMARTS QueryAtom identity retention
    // RDKit✔️✔️: atom_expr: atom_expr AND_TOKEN atom_expr {
    // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
    // RDKit✔️✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
    // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
    // RDKit✔️✔️:   delete $3;
    // RDKit✔️✔️:   $$ = $1;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | atom_expr OR_TOKEN atom_expr {
    // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
    // RDKit✔️✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
    // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
    // RDKit✔️✔️:   $1->setAtomicNum(0);
    // RDKit✔️✔️:   delete $3;
    // RDKit✔️✔️:   $$ = $1;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: point_query: NOT_TOKEN point_query {
    // RDKit✔️✔️:   $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
    // RDKit✔️✔️:   $2->setAtomicNum(0);
    // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($2);
    // RDKit✔️✔️:   $$ = $2;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | HASH_TOKEN number { $$ = new QueryAtom($2); }
    // END RDKIT CPP BLOCK SMARTS QueryAtom identity retention
    // `QueryNode` stores the final query tree instead of RDKit's mutable
    // QueryAtom construction object. AND retains the left QueryAtom identity;
    // OR/XOR and NOT clear only its atomic number, while its aromatic flag is
    // retained. This single ordered traversal reconstructs those observable
    // source fields in O(query height), without allocating or cloning.
    match query {
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_number)) => {
            (*atomic_number, false)
        }
        QueryNode::Predicate(AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        }) => (*atomic_number, *aromatic),
        QueryNode::Predicate(AtomQueryPredicate::IsAromatic(aromatic)) => (0, *aromatic),
        QueryNode::And(children) => children
            .first()
            .map(source_query_atom_identity)
            .unwrap_or((0, false)),
        QueryNode::Or(children) | QueryNode::Xor(children) => (
            0,
            children
                .first()
                .map(source_query_atom_identity)
                .unwrap_or((0, false))
                .1,
        ),
        QueryNode::Not(child) => (0, source_query_atom_identity(child).1),
        QueryNode::Predicate(_) => (0, false),
    }
}

fn materialize_smarts_atom_state(
    query: QueryNode<AtomQueryPredicate>,
) -> Result<(QueryNode<AtomQueryPredicate>, ChiralTag, Option<u32>), SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION atom_expr_and_point_query / atom_expr reductions
    // RDKit✔️✔️: atom_expr->expandQuery(point_query->getQuery()->copy(), Queries::COMPOSITE_AND, true);
    // RDKit✔️✔️: if (atom_expr->getChiralTag() == Atom::CHI_UNSPECIFIED) {
    // RDKit✔️✔️:   atom_expr->setChiralTag(point_query->getChiralTag());
    // RDKit✔️✔️:   int perm;
    // RDKit✔️✔️:   if (point_query->getPropIfPresent(common_properties::_chiralPermutation, perm)) {
    // RDKit✔️✔️:     atom_expr->setProp(common_properties::_chiralPermutation, perm);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // RDKit✔️✔️: $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
    // RDKit✔️✔️: if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
    // END RDKIT CPP FUNCTION atom_expr_and_point_query / atom_expr reductions
    // RDKit's SMARTS grammar stores chirality on QueryAtom, independently of
    // its query tree. The recursive parser initially represents every grammar
    // reduction as a typed node; strip the two temporary chirality nodes here
    // while rebuilding composites through QueryAtom::expandQuery's null
    // algebra, then apply CheckChiralitySpecifications exactly once.
    fn strip(
        query: QueryNode<AtomQueryPredicate>,
        chiral_tag: &mut ChiralTag,
        chiral_permutation: &mut Option<u32>,
        accept_permutation: &mut bool,
    ) -> QueryNode<AtomQueryPredicate> {
        fn rebuild(
            children: Vec<QueryNode<AtomQueryPredicate>>,
            how: CompositeQueryType,
            chiral_tag: &mut ChiralTag,
            chiral_permutation: &mut Option<u32>,
            accept_permutation: &mut bool,
        ) -> QueryNode<AtomQueryPredicate> {
            let mut children = children.into_iter();
            let mut rebuilt = children.next().map_or_else(make_atom_null_query, |child| {
                strip(child, chiral_tag, chiral_permutation, accept_permutation)
            });
            for child in children {
                let child = strip(child, chiral_tag, chiral_permutation, accept_permutation);
                super::query::query_atom_expand_query(&mut rebuilt, child, how, true);
            }
            rebuilt
        }

        match query {
            QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(tag)) => {
                *accept_permutation = *chiral_tag == ChiralTag::Unspecified;
                if *accept_permutation {
                    *chiral_tag = tag;
                }
                make_atom_null_query()
            }
            QueryNode::Predicate(AtomQueryPredicate::ChiralPermutationMatch(permutation)) => {
                if *accept_permutation {
                    *chiral_permutation = Some(permutation);
                }
                *accept_permutation = false;
                make_atom_null_query()
            }
            QueryNode::And(children) => rebuild(
                children,
                CompositeQueryType::And,
                chiral_tag,
                chiral_permutation,
                accept_permutation,
            ),
            QueryNode::Or(children) => rebuild(
                children,
                CompositeQueryType::Or,
                chiral_tag,
                chiral_permutation,
                accept_permutation,
            ),
            QueryNode::Xor(children) => rebuild(
                children,
                CompositeQueryType::Xor,
                chiral_tag,
                chiral_permutation,
                accept_permutation,
            ),
            QueryNode::Not(child) => QueryNode::not(strip(
                *child,
                chiral_tag,
                chiral_permutation,
                accept_permutation,
            )),
            query => {
                *accept_permutation = false;
                query
            }
        }
    }

    let mut chiral_tag = ChiralTag::Unspecified;
    let mut chiral_permutation = None;
    let mut accept_permutation = false;
    let query = strip(
        query,
        &mut chiral_tag,
        &mut chiral_permutation,
        &mut accept_permutation,
    );
    if let Some(permutation) = chiral_permutation {
        if !crate::notation::smiles::check_chiral_permutation(chiral_tag, permutation as i32) {
            return Err(SmartsParseError::Parse(format!(
                "invalid chiral permutation {permutation} for {}",
                chiral_tag.rdkit_name()
            )));
        }
        if chiral_tag == ChiralTag::Tetrahedral {
            if permutation <= 1 {
                chiral_tag = ChiralTag::TetrahedralCcw;
                chiral_permutation = None;
            } else if permutation == 2 {
                chiral_tag = ChiralTag::TetrahedralCw;
                chiral_permutation = None;
            }
        }
    }
    Ok((query, chiral_tag, chiral_permutation))
}

fn atom_isotope(query: &QueryNode<AtomQueryPredicate>) -> Option<u16> {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::Isotope(isotope)) => Some(*isotope),
        QueryNode::And(children) => children.iter().find_map(atom_isotope),
        QueryNode::Not(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Predicate(_) => None,
    }
}

fn representative_bond_order(query: &QueryNode<BondQueryPredicate>) -> BondOrder {
    match query {
        QueryNode::Predicate(BondQueryPredicate::Order(order)) => *order,
        QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)) => BondOrder::Aromatic,
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children
            .iter()
            .find_map(|child| match child {
                QueryNode::Predicate(BondQueryPredicate::Order(order)) => Some(*order),
                QueryNode::Predicate(BondQueryPredicate::IsAromatic(true)) => {
                    Some(BondOrder::Aromatic)
                }
                _ => None,
            })
            .unwrap_or(BondOrder::Single),
        _ => BondOrder::Single,
    }
}

#[cfg(test)]
pub(crate) fn compile_query_fixture(smarts: &str) -> Result<Molecule, String> {
    mol_from_smarts(smarts, &SmartsParseParams::default()).map_err(|error| error.to_string())
}

// ---------------------------------------------------------------------------
// SmartsParseParams
// ---------------------------------------------------------------------------

/// RDKit source: SmilesParse.h lines 56-67
/// RDKit✔️✔️: struct RDKIT_SMILESPARSE_EXPORT SmartsParserParams {
/// RDKit✔️✔️:   bool allowCXSMILES = true;
/// RDKit✔️✔️:   bool strictCXSMILES = true;
/// RDKit✔️✔️:   bool parseName = true;
/// RDKit✔️✔️:   bool mergeHs = false;
/// RDKit✔️✔️:   bool skipCleanup = false;
/// RDKit✔️✔️:   bool debugParse = false;
/// RDKit✔️✔️:   std::map<std::string, std::string> replacements;
/// RDKit✔️✔️: };
#[derive(Debug, Clone)]
pub struct SmartsParseParams {
    pub allow_cxsmiles: bool,
    pub strict_cxsmiles: bool,
    pub parse_name: bool,
    pub merge_hs: bool,
    pub skip_cleanup: bool,
    pub debug_parse: bool,
    pub replacements: BTreeMap<String, String>,
}

impl Default for SmartsParseParams {
    fn default() -> Self {
        Self {
            allow_cxsmiles: true,
            strict_cxsmiles: true,
            parse_name: true,
            merge_hs: false,
            skip_cleanup: false,
            debug_parse: false,
            replacements: BTreeMap::new(),
        }
    }
}

// ---------------------------------------------------------------------------
// Top-level parse entry point
// ---------------------------------------------------------------------------

// RDKit source: SmilesParse.cpp lines 548-576
// RDKit✔️✔️: std::unique_ptr<RWMol> MolFromSmarts(
// RDKit✔️✔️:     const std::string &smarts,
// RDKit✔️✔️:     const SmartsParserParams &params) {
// RDKit❌❌:   if (yysmarts_debug != params.debugParse) {
// RDKit❌❌:     yysmarts_debug = params.debugParse;
// RDKit❌❌:   }
// RDKit✔️✔️:   std::string lsmarts, name, cxPart;
// RDKit✔️✔️:   preprocessSmiles(smarts, params, lsmarts, name, cxPart);
// RDKit✔️✔️:   auto res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
// RDKit✔️✔️:   handleCXPartAndName(res.get(), params, cxPart, name);
// RDKit✔️✔️:   return res;
// RDKit✔️✔️: }
fn parse_smarts(smarts: &str) -> Result<ParsedSmartsGraph, SmartsParseError> {
    parse_smarts_with_params(smarts, &SmartsParseParams::default())
}

/// Parse a SMARTS string with custom parameters.
///
/// This private syntax-tree entry applies the same preprocessing and recursive
/// labeling as the public molecule compiler, but intentionally leaves molecule
/// postprocessing to `mol_from_smarts`.
fn parse_smarts_with_params(
    smarts: &str,
    params: &SmartsParseParams,
) -> Result<ParsedSmartsGraph, SmartsParseError> {
    // RDKit✔️✔️: preprocessSmiles — trim whitespace, handle replacements
    let preprocessed = preprocess_smarts(smarts, params);
    let input = label_recursive_patterns(&preprocessed.smarts);

    // RDKit✔️✔️: meta_start:
    // RDKit✔️✔️: START_MOL mol {
    // RDKit✔️✔️: // the molList has already been updated, no need to do anything
    // RDKit✔️✔️: }
    // RDKit✔️✔️: bad_atom_def:
    // RDKit✔️✔️: ATOM_OPEN_TOKEN bad_atom_def
    // RDKit✔️✔️: | ATOM_CLOSE_TOKEN bad_atom_def
    // RDKit✔️✔️: | COLON_TOKEN bad_atom_def
    // RDKit✔️✔️: | atom_expr {
    // RDKit✔️✔️:   delete $1;
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: ;
    // RDKit✔️✔️: | START_ATOM atomd EOS_TOKEN {
    // RDKit✔️✔️:   lastAtom = $2;
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_ATOM bad_atom_def {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_ATOM {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_BOND bond_expr EOS_TOKEN {
    // RDKit✔️✔️:   lastBond = $2;
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_BOND bond_expr {
    // RDKit✔️✔️:   delete $2;
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_BOND {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | meta_start BAD_CHARACTER {
    // RDKit✔️✔️:   yyerrok;
    // RDKit✔️✔️:   yyErrorCleanup(molList);
    // RDKit✔️✔️:   yyerror(input, molList, current_token_position, "syntax error");
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | meta_start error EOS_TOKEN{
    // RDKit✔️✔️:   yyerrok;
    // RDKit✔️✔️:   yyErrorCleanup(molList);
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | meta_start EOS_TOKEN {
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | error EOS_TOKEN {
    // RDKit✔️✔️:   yyerrok;
    // RDKit✔️✔️:   yyErrorCleanup(molList);
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | meta_start EOS_TOKEN {
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    smarts_parse_entry(&input)
}

pub fn mol_from_smarts(
    smarts: &str,
    params: &SmartsParseParams,
) -> Result<Molecule, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION MolFromSmarts
    // RDKit✔️❌: std::unique_ptr<RWMol> MolFromSmarts(const std::string &smarts,
    // RDKit✔️❌:                                      const SmartsParserParams &params) {
    // RDKit✔️✔️:   // Calling MolFromSmarts in a multithreaded context is generally safe *unless*
    // RDKit✔️✔️:   // the value of debugParse is different for different threads. The if
    // RDKit✔️✔️:   // statement below avoids a TSAN warning in the case where multiple threads
    // RDKit✔️✔️:   // all use the same value for debugParse.
    // RDKit❌❌:   if (yysmarts_debug != params.debugParse) {
    // RDKit❌❌:     yysmarts_debug = params.debugParse;
    // RDKit❌❌:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   std::string lsmarts, name, cxPart;
    // RDKit✔️✔️:   preprocessSmiles(smarts, params, lsmarts, name, cxPart);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   auto res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
    // RDKit✔️❌:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit✔️❌:   if (res) {
    // RDKit✔️❌:     if (params.mergeHs) {
    // RDKit✔️❌:       MolOps::mergeQueryHs(*res);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     MolOps::setBondStereoFromDirections(*res);
    // RDKit✔️❌:     if (!params.skipCleanup) {
    // RDKit✔️❌:       SmilesParseOps::CleanupAfterParsing(res.get());
    // RDKit✔️❌:     }
    // RDKit✔️✔️:     if (!name.empty()) {
    // RDKit✔️✔️:       res->setProp(common_properties::_Name, name);
    // RDKit✔️✔️:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return res;
    // RDKit✔️❌: };
    // END RDKIT CPP FUNCTION MolFromSmarts
    // Local complexity review: preprocessing, recursive labeling, parsing,
    // stereo assignment, and cleanup remain linear in their inputs or graph.
    // Shared CX and cleanup adapters rebuild through MoleculeBuilder, adding
    // material O(V + E) cloning versus RDKit's in-place edits. Query-H merging
    // uses the same canonical typed query graph and one builder compaction.
    if params.debug_parse {
        return Err(SmartsParseError::UnsupportedFeature(
            "Bison debug_parse diagnostic output",
        ));
    }
    let preprocessed = preprocess_smarts(smarts, params);
    let labeled = label_recursive_patterns(&preprocessed.smarts);
    let mut molecule = to_mol(&labeled).map_err(SmartsParseError::Parse)?;
    let mut name = preprocessed.name;
    handle_c_x_part_and_name(&mut molecule, params, &preprocessed.cx_part, &mut name)?;
    if params.merge_hs {
        merge_query_hs_in_place(&mut molecule, false, false)?;
    }
    crate::notation::smiles::apply_bond_stereo_from_directions_to_molecule(&mut molecule)
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    if !params.skip_cleanup {
        crate::notation::smiles::cleanup_after_parsing_molecule(&mut molecule)
            .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    }
    if !name.is_empty() {
        molecule = molecule.with_name(name);
    }
    Ok(molecule)
}

fn smarts_parse_helper(input: &str) -> Result<ParsedSmartsGraph, SmartsParseError> {
    // RDKit✔️✔️: int smarts_parse_helper(const std::string &inp, ...,
    // RDKit✔️✔️:   return generic_parse_helper<yysmarts_lex_init,
    // RDKit✔️✔️:     setup_smarts_string, yysmarts_lex_destroy>(yysmarts_parse,
    // RDKit✔️✔️:     inp, molVect, atom, bond, start_tok, "SMARTS");
    // Local complexity review: preprocessing has already produced one input
    // buffer; molecule tokenization and recursive descent each make one linear
    // pass. This wrapper adds no copy, rescan, alternate parser, or consumer
    // local decoding and keeps SMARTS on the sole canonical query path.
    let tokens = tokenize(input)?;
    let mut parser = SmartsParser::new(&tokens, input);
    parser.parse_smarts_molecule()
}

fn parse_atom_entry(input: &str) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
    // RDKit✔️✔️: | START_ATOM atomd EOS_TOKEN {
    // RDKit✔️✔️:   lastAtom = $2;
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_ATOM bad_atom_def {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_ATOM {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    let tokens = generic_parse_helper(input, ScannerStart::Atom)?;
    let mut parser = SmartsParser::new(&tokens, input);
    let atom = parser.parse_atomd()?;
    parser.require_end("atom SMARTS")?;
    Ok(atom.query)
}

fn parse_bond_entry(input: &str) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
    // RDKit✔️✔️: | START_BOND bond_expr EOS_TOKEN {
    // RDKit✔️✔️:   lastBond = $2;
    // RDKit✔️✔️:   YYACCEPT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_BOND bond_expr {
    // RDKit✔️✔️:   delete $2;
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | START_BOND {
    // RDKit✔️✔️:   YYABORT;
    // RDKit✔️✔️: }
    let tokens = generic_parse_helper(input, ScannerStart::Bond)?;
    let mut parser = SmartsParser::new(&tokens, input);
    let bond = parser.parse_bond_expr()?;
    parser.require_end("bond SMARTS")?;
    Ok(bond)
}

fn smarts_bond_parse(input: &str) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
    // RDKit✔️✔️: int smarts_bond_parse(const std::string &inp, Bond *&bond) {
    // RDKit✔️✔️:   auto start_tok = static_cast<int>(START_BOND);
    // RDKit✔️✔️:   std::vector<RWMol *> molVect;
    // RDKit✔️✔️:   Atom *atom = nullptr;
    // RDKit✔️✔️:   return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
    // RDKit✔️✔️: }
    // Local complexity review: dispatch is O(1) and delegates to the sole
    // scanner/parser path; no molecule allocation or alternate bond decoder is
    // introduced.
    parse_bond_entry(input)
}

fn smarts_atom_parse(input: &str) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
    // RDKit✔️✔️: int smarts_atom_parse(const std::string &inp, Atom *&atom) {
    // RDKit✔️✔️:   auto start_tok = static_cast<int>(START_ATOM);
    // RDKit✔️✔️:   std::vector<RWMol *> molVect;
    // RDKit✔️✔️:   Bond *bond = nullptr;
    // RDKit✔️✔️:   return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
    // RDKit✔️✔️: }
    // Local complexity review: dispatch is O(1) and delegates to the sole
    // scanner/parser path; no molecule allocation or alternate atom decoder is
    // introduced.
    parse_atom_entry(input)
}

fn to_atom(inp: &str) -> Result<Option<QueryNode<AtomQueryPredicate>>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION toAtom
    // RDKit✔️✔️: std::unique_ptr<Atom> toAtom(const std::string &inp,
    // RDKit✔️✔️:                              int func(const std::string &, Atom *&)) {
    // RDKit✔️✔️:   // empty strings produce nullptrs:
    // RDKit✔️✔️:   if (inp.empty()) {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Atom *res = nullptr;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     func(inp, res);
    // RDKit✔️✔️:   } catch (SmilesParseException &e) {
    // RDKit✔️✔️:     std::string nm = "SMILES";
    // RDKit✔️✔️:     if (func != smiles_atom_parse) {
    // RDKit✔️✔️:       nm = "SMARTS";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
    // RDKit✔️✔️:                           << " for input: '" << inp << "'" << std::endl;
    // RDKit✔️✔️:     res = nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return std::unique_ptr<Atom>(res);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION toAtom
    // Local complexity review: the empty branch is O(1); non-empty input is
    // parsed once by the canonical atom entry point. Option and Result replace
    // nullable ownership and exception logging without allocation, cloning,
    // rescanning, or a second atom decoder.
    if inp.is_empty() {
        return Ok(None);
    }
    smarts_atom_parse(inp).map(Some)
}

fn atom_from_smarts(
    smarts: &str,
) -> Result<Option<QueryNode<AtomQueryPredicate>>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION AtomFromSmarts
    // RDKit✔️✔️: std::unique_ptr<Atom> AtomFromSmarts(const std::string &smiles) {
    // RDKit✔️✔️:   yysmarts_debug = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return toAtom(smiles, smarts_atom_parse);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION AtomFromSmarts
    // Local complexity review: this constant-time entry wrapper delegates once
    // to `to_atom`, which performs the sole linear canonical atom parse. Rust
    // has no mutable bison debug global to reset, so the source assignment has
    // no allocation, scan, branch, or synchronization analogue.
    to_atom(smarts)
}

fn to_bond(inp: &str) -> Result<Option<QueryNode<BondQueryPredicate>>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION toBond
    // RDKit✔️✔️: std::unique_ptr<Bond> toBond(const std::string &inp,
    // RDKit✔️✔️:                              int func(const std::string &, Bond *&)) {
    // RDKit✔️✔️:   // empty strings produce nullptrs:
    // RDKit✔️✔️:   if (inp.empty()) {
    // RDKit✔️✔️:     return nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   Bond *res = nullptr;
    // RDKit✔️✔️:   try {
    // RDKit✔️✔️:     func(inp, res);
    // RDKit✔️✔️:   } catch (SmilesParseException &e) {
    // RDKit✔️✔️:     std::string nm = "SMILES";
    // RDKit✔️✔️:     if (func != smiles_bond_parse) {
    // RDKit✔️✔️:       nm = "SMARTS";
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
    // RDKit✔️✔️:                           << " for input: '" << inp << "'" << std::endl;
    // RDKit✔️✔️:     res = nullptr;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return std::unique_ptr<Bond>(res);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION toBond
    // Local complexity review: the empty branch is O(1); non-empty input is
    // parsed once by the canonical bond entry point. Option and Result replace
    // nullable ownership and exception logging without allocation, cloning,
    // rescanning, or a second bond decoder.
    if inp.is_empty() {
        return Ok(None);
    }
    smarts_bond_parse(inp).map(Some)
}

fn bond_from_smarts(
    smarts: &str,
) -> Result<Option<QueryNode<BondQueryPredicate>>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION BondFromSmarts
    // RDKit✔️✔️: std::unique_ptr<Bond> BondFromSmarts(const std::string &smiles) {
    // RDKit✔️✔️:   yysmarts_debug = false;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return toBond(smiles, smarts_bond_parse);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION BondFromSmarts
    // Local complexity review: this constant-time wrapper delegates once to
    // the sole linear bond parser. The absent bison debug global requires no
    // Rust state write and introduces no scan, allocation, or synchronization.
    to_bond(smarts)
}

fn smarts_parse_entry(input: &str) -> Result<ParsedSmartsGraph, SmartsParseError> {
    // RDKit✔️✔️: int smarts_parse(const std::string &inp, std::vector<RDKit::RWMol *> &molVect) {
    // RDKit✔️✔️:   auto start_tok = static_cast<int>(START_MOL);
    // RDKit✔️✔️:   Atom *atom = nullptr;
    // RDKit✔️✔️:   Bond *bond = nullptr;
    // RDKit✔️✔️:   return smarts_parse_helper(inp, molVect, atom, bond, start_tok);
    // RDKit✔️✔️: }
    // Local complexity review: dispatch is O(1) and delegates non-empty input
    // to the sole scanner/parser path. Accepting empty molecule SMARTS avoids
    // token allocation, matching the grammar's START_MOL/EOS branch without a
    // second parser, rescan, or molecule conversion.
    if input.is_empty() {
        return Ok(ParsedSmartsGraph::new(
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        ));
    }
    smarts_parse_helper(input)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PreprocessedSmarts {
    smarts: String,
    name: String,
    cx_part: String,
}

fn preprocess_smarts(smarts: &str, params: &SmartsParseParams) -> PreprocessedSmarts {
    // BEGIN RDKIT CPP FUNCTION preprocessSmiles<SmartsParserParams>
    // RDKit✔️✔️: // despite the name: works for both SMILES and SMARTS
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: void preprocessSmiles(const std::string &smiles, const T &params,
    // RDKit✔️✔️:                       std::string &lsmiles, std::string &name,
    // RDKit✔️✔️:                       std::string &cxPart) {
    // RDKit✔️✔️:   cxPart = "";
    // RDKit✔️✔️:   name = "";
    // RDKit✔️✔️:   if (params.parseName && !params.allowCXSMILES) {
    // RDKit✔️✔️:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️✔️:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️✔️:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️✔️:       name = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (params.allowCXSMILES) {
    // RDKit✔️✔️:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️✔️:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️✔️:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️✔️:       cxPart = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (lsmiles.empty()) {
    // RDKit✔️✔️:     lsmiles = smiles;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!params.replacements.empty()) {
    // RDKit✔️✔️:     std::string smi = lsmiles;
    // RDKit✔️✔️:     for (auto loopAgain = true; loopAgain;) {
    // RDKit✔️✔️:       loopAgain = false;
    // RDKit✔️✔️:       for (const auto &pr : params.replacements) {
    // RDKit✔️✔️:         if (smi.find(pr.first) != std::string::npos) {
    // RDKit✔️✔️:           loopAgain = true;
    // RDKit✔️✔️:           boost::replace_all(smi, pr.first, pr.second);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     lsmiles = smi;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION preprocessSmiles<SmartsParserParams>
    // Local complexity review: splitting performs one byte scan and at most
    // one owned copy per output field. Replacement uses the same ordered-map,
    // repeat-until-fixed-point scans and whole-string replacements as RDKit;
    // no extra parser, tokenization pass, or consumer-local preprocessing is
    // introduced.
    let mut processed = PreprocessedSmarts {
        smarts: String::new(),
        name: String::new(),
        cx_part: String::new(),
    };
    if params.parse_name && !params.allow_cxsmiles {
        if let Some(split_index) = smarts.bytes().position(|byte| matches!(byte, b' ' | b'\t'))
            && split_index != 0
        {
            processed.smarts = smarts[..split_index].to_string();
            processed.name = smarts[split_index..].trim().to_string();
        }
    } else if params.allow_cxsmiles
        && let Some(split_index) = smarts.bytes().position(|byte| matches!(byte, b' ' | b'\t'))
        && split_index != 0
    {
        processed.smarts = smarts[..split_index].to_string();
        processed.cx_part = smarts[split_index..].trim().to_string();
    }

    if processed.smarts.is_empty() {
        processed.smarts = smarts.to_string();
    }

    if !params.replacements.is_empty() {
        loop {
            let mut loop_again = false;
            for (key, value) in &params.replacements {
                if processed.smarts.contains(key) {
                    loop_again = true;
                    processed.smarts = processed.smarts.replace(key, value);
                }
            }
            if !loop_again {
                break;
            }
        }
    }
    processed
}

fn handle_c_x_part_and_name(
    molecule: &mut Molecule,
    params: &SmartsParseParams,
    cx_part: &str,
    name: &mut String,
) -> Result<(), SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION handleCXPartAndName<SmartsParserParams>
    // RDKit✔️❌: template <typename T>
    // RDKit✔️❌: void handleCXPartAndName(RWMol *res, const T &params, const std::string &cxPart,
    // RDKit✔️❌:                          std::string &name) {
    // RDKit✔️❌:   if (!res || cxPart.empty()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   std::string::const_iterator pos = cxPart.cbegin();
    // RDKit✔️❌:   bool cxfailed = false;
    // RDKit✔️❌:   if (params.allowCXSMILES) {
    // RDKit✔️❌:     if (*pos == '|') {
    // RDKit✔️❌:       try {
    // RDKit✔️❌:         SmilesParseOps::parseCXExtensions(*res, cxPart, pos);
    // RDKit✔️❌:       } catch (...) {
    // RDKit✔️❌:         cxfailed = true;
    // RDKit✔️❌:         if (params.strictCXSMILES) {
    // RDKit✔️❌:           throw;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:       res->setProp("_CXSMILES_Data", std::string(cxPart.cbegin(), pos));
    // RDKit✔️❌:     } else if (params.strictCXSMILES && !params.parseName &&
    // RDKit✔️❌:                pos != cxPart.cend()) {
    // RDKit✔️❌:       throw RDKit::SmilesParseException(
    // RDKit✔️❌:           "CXSMILES extension does not start with | and parseName=false");
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (!cxfailed && params.parseName && pos != cxPart.end()) {
    // RDKit✔️❌:     std::string nmpart(pos, cxPart.cend());
    // RDKit✔️❌:     name = boost::trim_copy(nmpart);
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION handleCXPartAndName<SmartsParserParams>
    // Local complexity review: the shared CX parser performs the same single
    // extension scan as RDKit, and name handling remains linear in the suffix.
    // The adapter currently rebuilds the molecule through MoleculeBuilder,
    // cloning O(atoms + bonds + attached state) where RDKit mutates RWMol in
    // place; this known material allocation cost accounts for the second-axis
    // performance gap. Keeping this adapter preserves one CX parser and all
    // query predicates while the canonical builder gains a narrowed mutation
    // path.
    crate::notation::smiles::apply_cx_part_and_name_to_molecule(
        molecule,
        params.allow_cxsmiles,
        params.strict_cxsmiles,
        params.parse_name,
        cx_part,
        name,
    )
    .map_err(|error| SmartsParseError::CxSmiles(error.to_string()))
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QueryHydrogenType {
    NotAHydrogen,
    UnmergableQueryHydrogen,
    QueryHydrogen,
}

fn query_has_hydrogen(query: &QueryNode<AtomQueryPredicate>, in_or: bool) -> (bool, bool) {
    // BEGIN RDKIT CPP FUNCTION queryHasHs
    // RDKit✔️✔️: template <class Q>
    // RDKit✔️✔️: std::pair<bool, bool> queryHasHs(Q queryAtom, bool inor = false) {
    // RDKit✔️✔️:   for (auto childit = queryAtom->beginChildren();
    // RDKit✔️✔️:        childit != queryAtom->endChildren(); ++childit) {
    // RDKit✔️✔️:     QueryAtom::QUERYATOM_QUERY::CHILD_TYPE query = *childit;
    // RDKit✔️✔️:     if (query->getDescription() == "AtomOr") {
    // RDKit✔️✔️:       return queryHasHs(query, true);
    // RDKit✔️✔️:     } else if (query->getDescription() == "AtomAtomicNum") {
    // RDKit✔️✔️:       if (static_cast<ATOM_EQUALS_QUERY *>(query.get())->getVal() == 1 &&
    // RDKit✔️✔️:           !query->getNegation()) {
    // RDKit✔️✔️:         return std::make_pair(true, inor);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (query->getDescription() == "AtomType") {
    // RDKit✔️✔️:       auto val = static_cast<ATOM_EQUALS_QUERY *>(query.get())->getVal();
    // RDKit✔️✔️:       // 1001 == aromtic hydrogen (not a thing, really)
    // RDKit✔️✔️:       // 1 == aliphatic hydrogen
    // RDKit✔️✔️:       if ((val == 1001 || val == 1) && !query->getNegation()) {
    // RDKit✔️✔️:         return std::make_pair(true, inor);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return std::make_pair(false, inor);
    // RDKit✔️✔️:   ;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION queryHasHs
    // Local complexity review: both versions inspect immediate children in
    // order and recurse only through an AtomOr child. They use O(depth) stack
    // space for nested ORs, allocate no collections, and short-circuit at the
    // same first OR or positive hydrogen predicate.
    let children = match query {
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children,
        QueryNode::Not(_) | QueryNode::Predicate(_) => return (false, in_or),
    };
    for child in children {
        match child {
            QueryNode::Or(_) => return query_has_hydrogen(child, true),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1))
            | QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 1, ..
            }) => return (true, in_or),
            QueryNode::Not(_) | QueryNode::And(_) | QueryNode::Xor(_) | QueryNode::Predicate(_) => {
            }
        }
    }
    (false, in_or)
}

fn is_query_hydrogen(atom: &crate::Atom, degree: usize) -> QueryHydrogenType {
    // BEGIN RDKIT CPP FUNCTION isQueryH
    // RDKit✔️✔️: HydrogenType isQueryH(const Atom *atom) {
    // RDKit✔️✔️:   PRECONDITION(atom, "bogus atom");
    // RDKit✔️✔️:   if (atom->getAtomicNum() == 1) {
    // RDKit✔️✔️:     // the simple case: the atom is flagged as being an H and
    // RDKit✔️✔️:     // has no query
    // RDKit✔️✔️:     if (!atom->hasQuery() ||
    // RDKit✔️✔️:         (!atom->getQuery()->getNegation() &&
    // RDKit✔️✔️:          atom->getQuery()->getDescription() == "AtomAtomicNum")) {
    // RDKit✔️✔️:       return HydrogenType::QueryHydrogen;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (!(atom->getDegree() <= 1)) {
    // RDKit✔️✔️:     // bonded and unbonded H atoms will continue rest will be returned
    // RDKit✔️✔️:     return HydrogenType::NotAHydrogen;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atom->hasQuery() && atom->getQuery()->getNegation()) {
    // RDKit✔️✔️:     // we will not merge negated queries
    // RDKit✔️✔️:     return HydrogenType::NotAHydrogen;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (atom->hasQuery()) {
    // RDKit✔️✔️:     std::pair<bool, bool> res = std::make_pair(false, false);
    // RDKit✔️✔️:     if (atom->getQuery()->getDescription() == "AtomOr") {
    // RDKit✔️✔️:       res = queryHasHs(atom->getQuery(), true);
    // RDKit✔️✔️:     } else if (atom->getQuery()->getDescription() == "AtomAnd") {
    // RDKit✔️✔️:       res = queryHasHs(atom->getQuery(), false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (res.first) {     // hasH
    // RDKit✔️✔️:       if (res.second) {  // inOr
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog)
    // RDKit✔️✔️:             << "WARNING: merging explicit H queries involved "
    // RDKit✔️✔️:                "in ORs is not supported. This query will not "
    // RDKit✔️✔️:                "be merged"
    // RDKit✔️✔️:             << std::endl;
    // RDKit✔️✔️:         return HydrogenType::UnMergableQueryHydrogen;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         return HydrogenType::QueryHydrogen;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return HydrogenType::NotAHydrogen;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isQueryH
    // Local complexity review: both implementations perform constant-time
    // atom/degree checks and at most one queryHasHs traversal. Rust's enum
    // matching adds no allocation, cloning, repeated graph scan, or lookup.
    let root_atomic_number_query = matches!(
        atom.query(),
        Some(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_)))
    );
    let canonical_atomic_number = match atom.query() {
        Some(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(value))) => *value,
        _ => atom.atomic_number(),
    };
    if canonical_atomic_number == 1 && (atom.query().is_none() || root_atomic_number_query) {
        return QueryHydrogenType::QueryHydrogen;
    }
    if degree > 1 || matches!(atom.query(), Some(QueryNode::Not(_))) {
        return QueryHydrogenType::NotAHydrogen;
    }
    if let Some(query @ (QueryNode::Or(_) | QueryNode::And(_))) = atom.query() {
        let (has_hydrogen, in_or) = query_has_hydrogen(query, matches!(query, QueryNode::Or(_)));
        if has_hydrogen {
            return if in_or {
                QueryHydrogenType::UnmergableQueryHydrogen
            } else {
                QueryHydrogenType::QueryHydrogen
            };
        }
    }
    QueryHydrogenType::NotAHydrogen
}

fn needs_hs(molecule: &Molecule) -> Result<bool, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION needsHs
    // RDKit✔️✔️: bool needsHs(const ROMol &mol) {
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     bool includeNeighbors = false;
    // RDKit✔️✔️:     if (atom->getTotalNumHs(includeNeighbors)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION needsHs
    // Local complexity review: both versions make one O(V) atom pass and
    // short-circuit on the first nonzero total-H count. COSMolKit reuses a
    // cached O(1)-per-atom valence assignment when present; otherwise its
    // canonical valence calculation is O(V + E), with no repeated graph scan,
    // clone, or temporary molecule.
    let computed;
    let valence = if let Some(cached) = crate::cached_valence_assignment(molecule) {
        cached
    } else {
        computed =
            crate::assign_valence_with_options(molecule, crate::ValenceModel::RdkitLike, false)
                .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
        &computed
    };
    Ok(molecule.atoms().iter().any(|atom| {
        let implicit = valence
            .implicit_hydrogens
            .get(atom.id().index())
            .copied()
            .unwrap_or(0)
            .max(0) as usize;
        usize::from(atom.explicit_hydrogens()) + implicit != 0
    }))
}

fn merge_recursive_query_hydrogens(
    query: &mut QueryNode<AtomQueryPredicate>,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<(), SmartsParseError> {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
            if let Some(query_molecule) = recursive.query_mol_mut() {
                merge_query_hs_in_place(query_molecule, merge_unmapped_only, merge_isotopes)?;
            }
        }
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            for child in children {
                merge_recursive_query_hydrogens(child, merge_unmapped_only, merge_isotopes)?;
            }
        }
        QueryNode::Not(child) => {
            merge_recursive_query_hydrogens(child, merge_unmapped_only, merge_isotopes)?;
        }
        QueryNode::Predicate(_) => {}
    }
    Ok(())
}

fn merge_query_hs_in_place(
    molecule: &mut Molecule,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<(), SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION mergeQueryHs(RWMol)
    // RDKit✔️❌: void mergeQueryHs(RWMol &mol, bool mergeUnmappedOnly, bool mergeIsotopes) {
    // RDKit✔️❌:   std::vector<unsigned int> atomsToRemove;
    // RDKit✔️❌:
    // RDKit✔️❌:   boost::dynamic_bitset<> hatoms(mol.getNumAtoms());
    // RDKit✔️❌:   for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
    // RDKit✔️❌:     hatoms[i] = isQueryH(mol.getAtomWithIdx(i)) == HydrogenType::QueryHydrogen;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   unsigned int currIdx = 0, stopIdx = mol.getNumAtoms();
    // RDKit✔️❌:   while (currIdx < stopIdx) {
    // RDKit✔️❌:     Atom *atom = mol.getAtomWithIdx(currIdx);
    // RDKit✔️❌:     if (!hatoms[currIdx]) {
    // RDKit✔️❌:       unsigned int numHsToRemove = 0;
    // RDKit✔️❌:       ROMol::ADJ_ITER begin, end;
    // RDKit✔️❌:       boost::tie(begin, end) = mol.getAtomNeighbors(atom);
    // RDKit✔️❌:
    // RDKit✔️❌:       while (begin != end) {
    // RDKit✔️❌:         if (hatoms[*begin]) {
    // RDKit✔️❌:           Atom &bgn = *mol.getAtomWithIdx(*begin);
    // RDKit✔️❌:           bool checkUnmapped =
    // RDKit✔️❌:               !mergeUnmappedOnly ||
    // RDKit✔️❌:               !bgn.hasProp(common_properties::molAtomMapNumber);
    // RDKit✔️❌:           bool checkIsotope = mergeIsotopes || bgn.getIsotope() == 0;
    // RDKit✔️❌:           if (checkUnmapped && checkIsotope) {
    // RDKit✔️❌:             atomsToRemove.push_back(rdcast<unsigned int>(*begin));
    // RDKit✔️❌:             ++numHsToRemove;
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:         ++begin;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       if (numHsToRemove) {
    // RDKit✔️❌:         //
    // RDKit✔️❌:         //  We have H neighbors:
    // RDKit✔️❌:         //   Add the appropriate queries to compensate for their removal.
    // RDKit✔️❌:         //
    // RDKit✔️❌:         //  Examples:
    // RDKit✔️❌:         //    C[H] -> [C;!H0]
    // RDKit✔️❌:         //    C([H])[H] -> [C;!H0;!H1]
    // RDKit✔️❌:         //
    // RDKit✔️❌:         //  It would be more efficient to do this using range queries like:
    // RDKit✔️❌:         //    C([H])[H] -> [C;H{2-}]
    // RDKit✔️❌:         //  but that would produce non-standard SMARTS without the user
    // RDKit✔️❌:         //  having started with a non-standard SMARTS.
    // RDKit✔️❌:         //
    // RDKit✔️❌:         if (!atom->hasQuery()) {
    // RDKit✔️❌:           // it wasn't a query atom, we need to replace it so that we can add
    // RDKit✔️❌:           // a query:
    // RDKit✔️❌:           ATOM_EQUALS_QUERY *tmp = makeAtomNumQuery(atom->getAtomicNum());
    // RDKit✔️❌:           auto *newAt = new QueryAtom;
    // RDKit✔️❌:           newAt->setQuery(tmp);
    // RDKit✔️❌:           newAt->updateProps(*atom);
    // RDKit✔️❌:           mol.replaceAtom(atom->getIdx(), newAt);
    // RDKit✔️❌:           delete newAt;
    // RDKit✔️❌:           atom = mol.getAtomWithIdx(currIdx);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         for (unsigned int i = 0; i < numHsToRemove; ++i) {
    // RDKit✔️❌:           ATOM_EQUALS_QUERY *tmp = makeAtomHCountQuery(i);
    // RDKit✔️❌:           tmp->setNegation(true);
    // RDKit✔️❌:           atom->expandQuery(tmp);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }  // end of numHsToRemove test
    // RDKit✔️❌:
    // RDKit✔️❌:       // recurse if needed (was github isusue 544)
    // RDKit✔️❌:       if (atom->hasQuery()) {
    // RDKit✔️❌:         if (atom->getQuery()->getDescription() == "RecursiveStructure") {
    // RDKit✔️❌:           auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
    // RDKit✔️❌:           CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit✔️❌:           RWMol *rqm = new RWMol(*rsq->getQueryMol());
    // RDKit✔️❌:           mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
    // RDKit✔️❌:           rsq->setQueryMol(rqm);
    // RDKit✔️❌:         }
    // RDKit✔️❌:
    // RDKit✔️❌:         // FIX: shouldn't be repeating this code here
    // RDKit✔️❌:         std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
    // RDKit✔️❌:             atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
    // RDKit✔️❌:         while (childStack.size()) {
    // RDKit✔️❌:           QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
    // RDKit✔️❌:           childStack.pop_front();
    // RDKit✔️❌:           if (qry->getDescription() == "RecursiveStructure") {
    // RDKit✔️❌:             auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
    // RDKit✔️❌:             CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit✔️❌:             RWMol *rqm = new RWMol(*rsq->getQueryMol());
    // RDKit✔️❌:             mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
    // RDKit✔️❌:             rsq->setQueryMol(rqm);
    // RDKit✔️❌:           } else if (qry->beginChildren() != qry->endChildren()) {
    // RDKit✔️❌:             childStack.insert(childStack.end(), qry->beginChildren(),
    // RDKit✔️❌:                               qry->endChildren());
    // RDKit✔️❌:           }
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }  // end of recursion loop
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++currIdx;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol.beginBatchEdit();
    // RDKit✔️❌:   for (auto aidx : atomsToRemove) {
    // RDKit✔️❌:     mol.removeAtom(aidx);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   mol.commitBatchEdit();
    // RDKit✔️❌: };
    // END RDKIT CPP FUNCTION mergeQueryHs(RWMol)
    // Local complexity review: classification, neighbor traversal, query
    // expansion, recursive traversal, and final compaction remain O(V + E + Q)
    // for each recursively owned query molecule. Rust uses Vec bit/count state
    // and one MoleculeBuilder rebuild, adding a material O(V + E) clone versus
    // RDKit's batch-edit mutation; this accounts for the second-axis gap.
    let hydrogen_types = molecule
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            is_query_hydrogen(
                atom,
                molecule
                    .topology_block()
                    .adjacency
                    .neighbors_of(index)
                    .len(),
            )
        })
        .collect::<Vec<_>>();
    let mut removals = Vec::new();
    let mut hydrogen_counts = vec![0_u8; molecule.num_atoms()];
    for atom_index in 0..molecule.num_atoms() {
        if hydrogen_types[atom_index] == QueryHydrogenType::QueryHydrogen {
            continue;
        }
        for neighbor in molecule.topology_block().adjacency.neighbors_of(atom_index) {
            if hydrogen_types[neighbor.atom_index] != QueryHydrogenType::QueryHydrogen {
                continue;
            }
            let hydrogen = &molecule.atoms()[neighbor.atom_index];
            let map_ok = !merge_unmapped_only || hydrogen.atom_map().is_none();
            let isotope_ok =
                merge_isotopes || hydrogen.isotope().is_none_or(|isotope| isotope == 0);
            if map_ok && isotope_ok {
                removals.push(hydrogen.id());
                hydrogen_counts[atom_index] = hydrogen_counts[atom_index].saturating_add(1);
            }
        }
    }
    removals.sort_unstable_by_key(|atom| atom.index());
    removals.dedup();
    let mut builder = molecule.to_builder();
    for (atom_index, count) in hydrogen_counts.into_iter().enumerate() {
        let atom = builder
            .atom_mut(crate::AtomId::new(atom_index))
            .expect("existing query atom");
        if count != 0 {
            let mut children = vec![atom.query().cloned().unwrap_or_else(|| {
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()))
            })];
            for hydrogen_count in 0..count {
                children.push(QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::HydrogenCount(hydrogen_count),
                ))));
            }
            atom.set_query(Some(QueryNode::And(children)));
        }
        if let Some(query) = atom.query_mut() {
            merge_recursive_query_hydrogens(query, merge_unmapped_only, merge_isotopes)?;
        }
    }
    builder.remove_atoms_for_construction(&removals);
    *molecule = builder
        .build()
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    Ok(())
}

fn merge_query_hs(
    molecule: &Molecule,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<Molecule, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION mergeQueryHs(ROMol)
    // RDKit✔️✔️: ROMol *mergeQueryHs(const ROMol &mol, bool mergeUnmappedOnly,
    // RDKit✔️✔️:                     bool mergeIsotopes) {
    // RDKit✔️✔️:   auto *res = new RWMol(mol);
    // RDKit✔️✔️:   mergeQueryHs(*res, mergeUnmappedOnly, mergeIsotopes);
    // RDKit✔️✔️:   return static_cast<ROMol *>(res);
    // RDKit✔️✔️: };
    // END RDKIT CPP FUNCTION mergeQueryHs(ROMol)
    // Local complexity review: both implementations clone the input molecule
    // once and delegate to the sole in-place implementation. COSMolKit's
    // Arc-backed clone is O(1) until the delegated builder materializes its
    // graph; no duplicate traversal, query logic, or extra molecule copy is
    // introduced here.
    let mut result = molecule.clone();
    merge_query_hs_in_place(&mut result, merge_unmapped_only, merge_isotopes)?;
    Ok(result)
}

fn query_node_h_status(query: &QueryNode<AtomQueryPredicate>) -> (bool, bool) {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => recursive
            .query_mol()
            .map(has_query_hs)
            .unwrap_or((false, false)),
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            let mut query_hs = false;
            for child in children {
                let result = query_node_h_status(child);
                if result.1 {
                    return result;
                }
                query_hs |= result.0;
            }
            (query_hs, false)
        }
        QueryNode::Not(child) => query_node_h_status(child),
        QueryNode::Predicate(_) => (false, false),
    }
}

fn has_query_hs(molecule: &Molecule) -> (bool, bool) {
    // BEGIN RDKIT CPP FUNCTION hasQueryHs
    // RDKit✔️✔️: std::pair<bool, bool> hasQueryHs(const ROMol &mol) {
    // RDKit✔️✔️:   bool queryHs = false;
    // RDKit✔️✔️:   // We don't care about announcing ORs or other items during isQueryH
    // RDKit✔️✔️:   RDLog::LogStateSetter blocker;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   for (const auto atom : mol.atoms()) {
    // RDKit✔️✔️:     switch (isQueryH(atom)) {
    // RDKit✔️✔️:       case HydrogenType::UnMergableQueryHydrogen:
    // RDKit✔️✔️:         return std::make_pair(true, true);
    // RDKit✔️✔️:       case HydrogenType::QueryHydrogen:
    // RDKit✔️✔️:         queryHs = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       default:  // HydrogenType::NotAHydrogen:
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     if (atom->hasQuery()) {
    // RDKit✔️✔️:       if (atom->getQuery()->getDescription() == "RecursiveStructure") {
    // RDKit✔️✔️:         auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
    // RDKit✔️✔️:         CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit✔️✔️:         auto res = hasQueryHs(*rsq->getQueryMol());
    // RDKit✔️✔️:         if (res.second) {  // unmergableH implies queryH
    // RDKit✔️✔️:           return res;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         queryHs |= res.first;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:
    // RDKit✔️✔️:       // FIX: shouldn't be repeating this code here -- yet again!
    // RDKit✔️✔️:       std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
    // RDKit✔️✔️:           atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
    // RDKit✔️✔️:       while (!childStack.empty()) {
    // RDKit✔️✔️:         QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
    // RDKit✔️✔️:         childStack.pop_front();
    // RDKit✔️✔️:         if (qry->getDescription() == "RecursiveStructure") {
    // RDKit✔️✔️:           auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
    // RDKit✔️✔️:           CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit✔️✔️:           auto res = hasQueryHs(*rsq->getQueryMol());
    // RDKit✔️✔️:           if (res.second) {
    // RDKit✔️✔️:             return res;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           queryHs |= res.first;
    // RDKit✔️✔️:         } else {
    // RDKit✔️✔️:           childStack.insert(childStack.end(), qry->beginChildren(),
    // RDKit✔️✔️:                             qry->endChildren());
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }  // end of recursion loop
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return std::make_pair(queryHs, false);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION hasQueryHs
    // Local complexity review: both versions traverse atoms and each owned
    // recursive query tree once, short-circuiting on the first unmergeable H.
    // Runtime is O(V + Q) over the recursive query closure and stack depth is
    // O(query depth); Rust allocates no child queue or reparsed SMARTS cache.
    let mut query_hs = false;
    for (index, atom) in molecule.atoms().iter().enumerate() {
        let degree = molecule
            .topology_block()
            .adjacency
            .neighbors_of(index)
            .len();
        match is_query_hydrogen(atom, degree) {
            QueryHydrogenType::UnmergableQueryHydrogen => return (true, true),
            QueryHydrogenType::QueryHydrogen => query_hs = true,
            QueryHydrogenType::NotAHydrogen => {}
        }
        if let Some(query) = atom.query() {
            let result = query_node_h_status(query);
            if result.1 {
                return result;
            }
            query_hs |= result.0;
        }
    }
    (query_hs, false)
}

fn label_recursive_patterns(sma: &str) -> String {
    // RDKit✔️✔️: std::string labelRecursivePatterns(const std::string &sma) {
    // RDKit✔️✔️: #ifndef NO_AUTOMATIC_SMARTS_RELABELLING
    // RDKit✔️✔️:   std::list<SmaState> state;
    // RDKit✔️✔️:   std::list<unsigned int> startRecurse;
    // RDKit✔️✔️:   std::map<std::string, std::string> patterns;
    // RDKit✔️✔️:   std::string res;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   state.push_back(BASE);
    // RDKit✔️✔️:
    // RDKit✔️✔️:   unsigned int pos = 0;
    // RDKit✔️✔️:   while (pos < sma.size()) {
    // RDKit✔️✔️:     res += sma[pos];
    // RDKit✔️✔️:     if (sma[pos] == '$' && pos + 1 < sma.size() && sma[pos + 1] == '(') {
    // RDKit✔️✔️:       state.push_back(RECURSE);
    // RDKit✔️✔️:       startRecurse.push_back(pos);
    // RDKit✔️✔️:       ++pos;
    // RDKit✔️✔️:       res += sma[pos];
    // RDKit✔️✔️:     } else if (sma[pos] == '(') {
    // RDKit✔️✔️:       state.push_back(BRANCH);
    // RDKit✔️✔️:     } else if (sma[pos] == ')') {
    // RDKit✔️✔️:       if (state.empty() || state.back() == BASE) {
    // RDKit✔️✔️:         // seriously bogus input. Just return the input
    // RDKit✔️✔️:         // and let the SMARTS parser itself report the error
    // RDKit✔️✔️:         return sma;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       SmaState currState = state.back();
    // RDKit✔️✔️:       state.pop_back();
    // RDKit✔️✔️:       if (currState == RECURSE) {
    // RDKit✔️✔️:         unsigned int dollarPos = startRecurse.back();
    // RDKit✔️✔️:         startRecurse.pop_back();
    // RDKit✔️✔️:         if (pos + 1 >= sma.size() || sma[pos + 1] != '_') {
    // RDKit✔️✔️:           std::string recurs = sma.substr(dollarPos, pos - dollarPos + 1);
    // RDKit✔️✔️:           std::string label;
    // RDKit✔️✔️:           if (patterns.find(recurs) != patterns.end()) {
    // RDKit✔️✔️:             // seen this one before, add the label
    // RDKit✔️✔️:             label = patterns[recurs];
    // RDKit✔️✔️:           } else {
    // RDKit✔️✔️:             label = std::to_string(patterns.size() + 100);
    // RDKit✔️✔️:             patterns[recurs] = label;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:           res += "_" + label;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else if (currState == BRANCH) {
    // RDKit✔️✔️:         // no need to do anything here.
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++pos;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // std::cerr<< " >"<<sma<<"->"<<res<<std::endl;
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: #else
    // RDKit✔️✔️:   return sma;
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️: }
    // Local complexity review: input and output are scanned once in O(n),
    // stack operations are amortized O(1), and recursive-pattern lookup is
    // O(log p), matching the source. Vec stacks avoid std::list node
    // allocations; BTreeMap preserves std::map lookup complexity. Byte indexing
    // reproduces std::string behavior and introduces no parser or reparse path.
    #[derive(Clone, Copy, PartialEq)]
    enum SmaState {
        Base,
        Branch,
        Recurse,
    }
    use SmaState::*;

    let mut state: Vec<SmaState> = vec![Base];
    let mut start_recurse: Vec<usize> = Vec::new();
    let mut patterns: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    let mut res = Vec::with_capacity(sma.len());
    let bytes = sma.as_bytes();

    let mut pos: usize = 0;
    while pos < bytes.len() {
        res.push(bytes[pos]);
        if bytes[pos] == b'$' && pos + 1 < bytes.len() && bytes[pos + 1] == b'(' {
            state.push(Recurse);
            start_recurse.push(pos);
            pos += 1;
            res.push(bytes[pos]);
        } else if bytes[pos] == b'(' {
            state.push(Branch);
        } else if bytes[pos] == b')' {
            if state.is_empty() || state.last() == Some(&Base) {
                return sma.to_string();
            }
            let curr_state = state.pop().expect("non-base SMARTS state");
            if curr_state == Recurse {
                let dollar_pos = start_recurse.pop().expect("recursive SMARTS start");
                if pos + 1 >= bytes.len() || bytes[pos + 1] != b'_' {
                    let recurs = &bytes[dollar_pos..=pos];
                    let label = if let Some(lbl) = patterns.get(recurs) {
                        lbl.clone()
                    } else {
                        let lbl = format!("{}", patterns.len() + 100);
                        patterns.insert(recurs.to_vec(), lbl.clone());
                        lbl
                    };
                    res.push(b'_');
                    res.extend_from_slice(label.as_bytes());
                }
            }
        }
        pos += 1;
    }
    String::from_utf8(res).expect("SMARTS relabeling preserves UTF-8")
}

// ---------------------------------------------------------------------------
// Tokenizer
// ---------------------------------------------------------------------------

/// RDKit❗✔️: Comparison of token types between SMILES/smarts.ll and our
/// tokenizer. The flex/generated lexer produces tokens like ORGANIC_ATOM_TOKEN,
/// AROMATIC_ATOM_TOKEN, ATOM_TOKEN, BOND_TOKEN, etc. We collapse those into
/// a simpler enum since our parser uses recursive descent rather than bison.
#[derive(Debug, Clone, PartialEq)]
enum Token {
    /// Organic element symbol: B, C, N, O, S, P, F, Cl, Br, I, *
    OrganicElement(String),
    /// Aromatic element: c, n, o, s, p
    AromaticElement(String),
    /// Bracket atom content: the raw text between [ and ] (excluding brackets)
    BracketContent(String),
    /// Bond specifier: -, =, #, :, ~, /, \\
    BondSpec(BondLexeme),
    /// Open parenthesis (branch)
    OpenParen,
    /// Close parenthesis
    CloseParen,
    /// Ring closure digit 0-9
    RingClosureDigit(u32),
    /// Ring closure %NN
    RingClosurePercent(u32),
    /// Logical AND operator &
    And,
    /// Low-precedence logical AND operator ;
    Semi,
    /// Logical OR operator (comma)
    Or,
    /// Logical NOT operator !
    Not,
    /// Low-order bit separator (.)
    Dot,
    /// End of token stream
    EndOfStream,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BondLexeme {
    Symbol(char),
    DativeRight,
    DativeLeft,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScannerStart {
    Molecule,
    Atom,
    Bond,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ScannerState {
    Initial,
    Atom,
    Branch,
    Recursion,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum ScannerToken {
    Start(ScannerStart),
    OrganicElement(String),
    AromaticElement(String),
    AtomElement(String),
    AtomPrimitive(char),
    BondSpec(char),
    DativeRight,
    DativeLeft,
    ChiralClass(String),
    At,
    Hybridization(u8),
    GroupOpen,
    GroupClose,
    BeginRecurse,
    EndRecurse,
    AtomOpen,
    AtomClose,
    RangeOpen,
    RangeClose,
    Colon,
    Underscore,
    Hash,
    Minus,
    Plus,
    Separator,
    Percent,
    Digit(u8),
    Not,
    Semi,
    And,
    Or,
    EndOfStream,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ScannedToken {
    token: ScannerToken,
    start: usize,
    end: usize,
}

struct SmartsScanner {
    chars: Vec<char>,
    start: ScannerStart,
    states: Vec<ScannerState>,
    pos: usize,
}

impl SmartsScanner {
    fn new(input: &str, start: ScannerStart) -> Self {
        Self {
            chars: input.chars().collect(),
            start,
            states: vec![ScannerState::Initial],
            pos: 0,
        }
    }

    fn state(&self) -> ScannerState {
        *self
            .states
            .last()
            .expect("scanner state stack is non-empty")
    }

    fn emit(&mut self, token: ScannerToken, width: usize) -> ScannedToken {
        let start = self.pos;
        self.pos += width;
        ScannedToken {
            token,
            start,
            end: self.pos,
        }
    }

    fn scan(mut self) -> Result<Vec<ScannedToken>, SmartsParseError> {
        let mut tokens = vec![ScannedToken {
            token: ScannerToken::Start(self.start),
            start: 0,
            end: 0,
        }];

        while self.pos < self.chars.len() {
            let state = self.state();
            let ch = self.chars[self.pos];

            // RDKit✔️✔️: \n\t\treturn EOS_TOKEN;
            if ch == '\n' {
                tokens.push(self.emit(ScannerToken::EndOfStream, 1));
                return Ok(tokens);
            }

            // RDKit✔️✔️: @[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }
            if ch == '@' {
                let mut cursor = self.pos + 1;
                while self.chars.get(cursor) == Some(&' ') {
                    cursor += 1;
                }
                if cursor + 1 < self.chars.len() {
                    let class: String = self.chars[cursor..cursor + 2].iter().collect();
                    if matches!(class.as_str(), "TH" | "AL" | "SP" | "TB" | "OH") {
                        tokens.push(
                            self.emit(ScannerToken::ChiralClass(class), cursor + 2 - self.pos),
                        );
                        continue;
                    }
                }
                // RDKit✔️✔️: @\t\t{ return AT_TOKEN; }
                tokens.push(self.emit(ScannerToken::At, 1));
                continue;
            }

            // RDKit✔️✔️: <IN_ATOM_STATE>\$\( { yy_push_state(IN_RECURSION_STATE,yyscanner); return BEGIN_RECURSE; }
            if state == ScannerState::Atom
                && ch == '$'
                && self.chars.get(self.pos + 1) == Some(&'(')
            {
                self.states.push(ScannerState::Recursion);
                tokens.push(self.emit(ScannerToken::BeginRecurse, 2));
                continue;
            }

            // RDKit✔️✔️: \( { yy_push_state(IN_BRANCH_STATE,yyscanner); return GROUP_OPEN_TOKEN; }
            // RDKit✔️✔️: <IN_BRANCH_STATE>\) { yy_pop_state(yyscanner); return GROUP_CLOSE_TOKEN; }
            // RDKit✔️✔️: <IN_RECURSION_STATE>\) { yy_pop_state(yyscanner); return END_RECURSE; }
            if ch == '(' {
                self.states.push(ScannerState::Branch);
                tokens.push(self.emit(ScannerToken::GroupOpen, 1));
                continue;
            }
            if ch == ')' {
                let token = match state {
                    ScannerState::Branch => {
                        self.states.pop();
                        ScannerToken::GroupClose
                    }
                    ScannerState::Recursion => {
                        self.states.pop();
                        ScannerToken::EndRecurse
                    }
                    _ => ScannerToken::GroupClose,
                };
                tokens.push(self.emit(token, 1));
                continue;
            }

            // RDKit✔️✔️: \[ { yy_push_state(IN_ATOM_STATE,yyscanner); return ATOM_OPEN_TOKEN; }
            // RDKit✔️✔️: <IN_ATOM_STATE>\] { yy_pop_state(yyscanner); return ATOM_CLOSE_TOKEN; }
            // RDKit✔️✔️: \] { return ATOM_CLOSE_TOKEN; }
            if ch == '[' {
                self.states.push(ScannerState::Atom);
                tokens.push(self.emit(ScannerToken::AtomOpen, 1));
                continue;
            }
            if ch == ']' {
                if state == ScannerState::Atom {
                    self.states.pop();
                }
                tokens.push(self.emit(ScannerToken::AtomClose, 1));
                continue;
            }

            if state == ScannerState::Atom {
                if let Some(token) = self.scan_atom_token()? {
                    tokens.push(token);
                    continue;
                }
            }

            if let Some(token) = self.scan_common_token()? {
                tokens.push(token);
                continue;
            }

            // RDKit✔️✔️: .\t\treturn BAD_CHARACTER;
            return Err(SmartsParseError::UnexpectedCharacter {
                position: self.pos,
                character: ch,
                context: "unexpected character in SMARTS string".to_string(),
            });
        }

        if self.states.contains(&ScannerState::Atom) {
            let start = tokens
                .iter()
                .rev()
                .find(|token| token.token == ScannerToken::AtomOpen)
                .map_or(0, |token| token.start);
            return Err(SmartsParseError::UnclosedBracket(start));
        }

        // RDKit✔️✔️: <<EOF>>\t\t{ return EOS_TOKEN; }
        tokens.push(ScannedToken {
            token: ScannerToken::EndOfStream,
            start: self.pos,
            end: self.pos,
        });
        Ok(tokens)
    }

    fn scan_atom_token(&mut self) -> Result<Option<ScannedToken>, SmartsParseError> {
        let ch = self.chars[self.pos];
        let rest: String = self.chars[self.pos..].iter().collect();

        // RDKit✔️✔️: <IN_ATOM_STATE>He |
        // RDKit✔️✔️: <IN_ATOM_STATE>Li |
        // RDKit✔️✔️: <IN_ATOM_STATE>Be |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ne |
        // RDKit✔️✔️: <IN_ATOM_STATE>Na |
        // RDKit✔️✔️: <IN_ATOM_STATE>Mg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Al |
        // RDKit✔️✔️: <IN_ATOM_STATE>Si |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ar |
        // RDKit✔️✔️: <IN_ATOM_STATE>K |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ca |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sc |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ti |
        // RDKit✔️✔️: <IN_ATOM_STATE>V |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Mn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Co |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fe |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ni |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cu |
        // RDKit✔️✔️: <IN_ATOM_STATE>Zn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ga |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ge |
        // RDKit✔️✔️: <IN_ATOM_STATE>As |
        // RDKit✔️✔️: <IN_ATOM_STATE>Se |
        // RDKit✔️✔️: <IN_ATOM_STATE>Kr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Y |
        // RDKit✔️✔️: <IN_ATOM_STATE>Zr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Nb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Mo |
        // RDKit✔️✔️: <IN_ATOM_STATE>Tc |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ru |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rh |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pd |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ag |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cd |
        // RDKit✔️✔️: <IN_ATOM_STATE>In |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Te |
        // RDKit✔️✔️: <IN_ATOM_STATE>Xe |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cs |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ba |
        // RDKit✔️✔️: <IN_ATOM_STATE>La |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ce |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Nd |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Eu |
        // RDKit✔️✔️: <IN_ATOM_STATE>Gd |
        // RDKit✔️✔️: <IN_ATOM_STATE>Tb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Dy |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ho |
        // RDKit✔️✔️: <IN_ATOM_STATE>Er |
        // RDKit✔️✔️: <IN_ATOM_STATE>Tm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Yb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Lu |
        // RDKit✔️✔️: <IN_ATOM_STATE>Hf |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ta |
        // RDKit✔️✔️: <IN_ATOM_STATE>W |
        // RDKit✔️✔️: <IN_ATOM_STATE>Re |
        // RDKit✔️✔️: <IN_ATOM_STATE>Os |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ir |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pt |
        // RDKit✔️✔️: <IN_ATOM_STATE>Au |
        // RDKit✔️✔️: <IN_ATOM_STATE>Hg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Tl |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bi |
        // RDKit✔️✔️: <IN_ATOM_STATE>Po |
        // RDKit✔️✔️: <IN_ATOM_STATE>At |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ra |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ac |
        // RDKit✔️✔️: <IN_ATOM_STATE>Th |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pa |
        // RDKit✔️✔️: <IN_ATOM_STATE>U |
        // RDKit✔️✔️: <IN_ATOM_STATE>Np |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pu |
        // RDKit✔️✔️: <IN_ATOM_STATE>Am |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bk |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cf |
        // RDKit✔️✔️: <IN_ATOM_STATE>Es |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Md |
        // RDKit✔️✔️: <IN_ATOM_STATE>No |
        // RDKit✔️✔️: <IN_ATOM_STATE>Lr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rf |
        // RDKit✔️✔️: <IN_ATOM_STATE>Db |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bh |
        // RDKit✔️✔️: <IN_ATOM_STATE>Hs |
        // RDKit✔️✔️: <IN_ATOM_STATE>Mt |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ds |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Uut |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fl |
        // RDKit✔️✔️: <IN_ATOM_STATE>Uup |
        // RDKit✔️✔️: <IN_ATOM_STATE>Lv	{   yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
        // RDKit✔️✔️: 				return ATOM_TOKEN;
        // RDKit✔️✔️: 			}
        // RDKit✔️✔️: <IN_ATOM_STATE>D {
        // Flex uses longest-match selection, so three-letter temporary element
        // names and then two-letter names are checked before one-letter names.
        for symbol in ELEMENT_SYMBOLS {
            if rest.starts_with(symbol) {
                return Ok(Some(self.emit(
                    ScannerToken::AtomElement((*symbol).to_string()),
                    symbol.chars().count(),
                )));
            }
        }

        // RDKit✔️✔️: <IN_ATOM_STATE>si { yylval->ival = 14; return AROMATIC_ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>as { yylval->ival = 33; return AROMATIC_ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>se { yylval->ival = 34; return AROMATIC_ATOM_TOKEN; }
        // RDKit✔️✔️: <IN_ATOM_STATE>te { yylval->ival = 52; return AROMATIC_ATOM_TOKEN; }
        for symbol in ["si", "as", "se", "te"] {
            if rest.starts_with(symbol) {
                return Ok(Some(
                    self.emit(ScannerToken::AromaticElement(symbol.to_string()), 2),
                ));
            }
        }

        // RDKit✔️✔️: <IN_ATOM_STATE>D {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomExplicitDegreeQuery(1));
        // RDKit✔️✔️: \treturn COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>d {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomNonHydrogenDegreeQuery(1));
        // RDKit✔️✔️: \treturn COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>X {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomTotalDegreeQuery(1));
        // RDKit✔️✔️: \treturn COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>x {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHasRingBondQuery());
        // RDKit✔️✔️: \treturn RINGBOND_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>v {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomTotalValenceQuery(1));
        // RDKit✔️✔️: \treturn COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>z {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHasHeteroatomNbrsQuery());
        // RDKit✔️✔️: \treturn HETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>Z {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHasAliphaticHeteroatomNbrsQuery());
        // RDKit✔️✔️: \treturn ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>h {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️:         yylval->atom->setQuery(makeAtomHasImplicitHQuery());
        // RDKit✔️✔️: \treturn IMPLICIT_H_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>R {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(new AtomRingQuery(-1));
        // RDKit✔️✔️: \treturn COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>r {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit✔️✔️: \treturn MIN_RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: <IN_ATOM_STATE>k {
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit✔️✔️: \treturn RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        if matches!(
            ch,
            'D' | 'd' | 'X' | 'x' | 'v' | 'z' | 'Z' | 'h' | 'R' | 'r' | 'k'
        ) {
            return Ok(Some(self.emit(ScannerToken::AtomPrimitive(ch), 1)));
        }

        // RDKit✔️✔️: \^0\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::S));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \^1\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \^2\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP2));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \^3\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \^4\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \^5\t\t{
        // RDKit✔️✔️: \tyylval->atom = new QueryAtom();
        // RDKit✔️✔️: \tyylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D2));
        // RDKit✔️✔️: \treturn HYB_TOKEN;
        // RDKit✔️✔️: }
        if ch == '^' {
            if let Some(value) = self
                .chars
                .get(self.pos + 1)
                .and_then(|digit| digit.to_digit(10))
                .filter(|value| *value <= 5)
            {
                return Ok(Some(self.emit(ScannerToken::Hybridization(value as u8), 2)));
            }
        }
        Ok(None)
    }

    fn scan_common_token(&mut self) -> Result<Option<ScannedToken>, SmartsParseError> {
        let ch = self.chars[self.pos];
        let rest: String = self.chars[self.pos..].iter().collect();

        // RDKit✔️✔️: B ... I { yylval->ival = ...; return ORGANIC_ATOM_TOKEN; }
        for symbol in ["Cl", "Br", "B", "C", "N", "O", "F", "P", "S", "I"] {
            if rest.starts_with(symbol) {
                return Ok(Some(self.emit(
                    ScannerToken::OrganicElement(symbol.to_string()),
                    symbol.chars().count(),
                )));
            }
        }
        // RDKit✔️✔️: b ... s { yylval->ival = ...; return AROMATIC_ATOM_TOKEN; }
        if matches!(ch, 'b' | 'c' | 'n' | 'o' | 'p' | 's') {
            return Ok(Some(
                self.emit(ScannerToken::AromaticElement(ch.to_string()), 1),
            ));
        }
        // RDKit✔️✔️: \* { ... return SIMPLE_ATOM_QUERY_TOKEN; }
        // RDKit✔️✔️: a { ... return SIMPLE_ATOM_QUERY_TOKEN; }
        // RDKit✔️✔️: A { ... return SIMPLE_ATOM_QUERY_TOKEN; }
        if ch == '*' || ch == 'A' || ch == 'a' {
            let token = if ch == 'a' {
                ScannerToken::AromaticElement("a".to_string())
            } else {
                ScannerToken::OrganicElement(ch.to_string())
            };
            return Ok(Some(self.emit(token, 1)));
        }
        // RDKit✔️✔️: H { return H_TOKEN; }
        if ch == 'H' {
            return Ok(Some(self.emit(ScannerToken::AtomPrimitive('H'), 1)));
        }

        // RDKit✔️✔️: \-\> { yylval->bond = new QueryBond(Bond::DATIVER); return BOND_TOKEN; }
        // RDKit✔️✔️: \<\- { yylval->bond = new QueryBond(Bond::DATIVEL); return BOND_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \*			{
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomNullQuery());
        // RDKit✔️✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: a			{
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomAromaticQuery());
        // RDKit✔️✔️: 	yylval->atom->setIsAromatic(true);
        // RDKit✔️✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: A			{
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomAliphaticQuery());
        // RDKit✔️✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️:
        // RDKit✔️✔️: \: 			{ return COLON_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \_ 			{ return UNDERSCORE_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \#			{ return HASH_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \=	{ yylval->bond = new QueryBond(Bond::DOUBLE);
        // RDKit✔️✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::DOUBLE));
        // RDKit✔️✔️: 	return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \~	{ yylval->bond = new QueryBond();
        // RDKit✔️✔️: 	yylval->bond->setQuery(makeBondNullQuery());
        // RDKit✔️✔️: 	return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \$	{ yylval->bond = new QueryBond(Bond::QUADRUPLE);
        // RDKit✔️✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::QUADRUPLE));
        // RDKit✔️✔️:     return BOND_TOKEN; }
        // RDKit✔️✔️:
        // RDKit✔️✔️: [\\]{1,2}    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit✔️✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
        // RDKit✔️✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit✔️✔️: 	return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: [\/]    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit✔️✔️: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
        // RDKit✔️✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit✔️✔️: 	return BOND_TOKEN;  }
        // RDKit✔️✔️:
        // RDKit✔️✔️: \-\> {
        // RDKit✔️✔️:     yylval->bond = new QueryBond(Bond::DATIVER);
        // RDKit✔️✔️:     return BOND_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: \<\- {
        // RDKit✔️✔️:     yylval->bond = new QueryBond(Bond::DATIVEL);
        // RDKit✔️✔️:     return BOND_TOKEN;
        // RDKit✔️✔️: }
        if rest.starts_with("->") {
            return Ok(Some(self.emit(ScannerToken::DativeRight, 2)));
        }
        if rest.starts_with("<-") {
            return Ok(Some(self.emit(ScannerToken::DativeLeft, 2)));
        }
        // RDKit✔️✔️: \= ... \~ ... \$ ... [\\]{1,2} ... [\/] { ... return BOND_TOKEN; }
        if matches!(ch, '=' | '~' | '$' | '/' | '\\') {
            let width = if ch == '\\' && self.chars.get(self.pos + 1) == Some(&'\\') {
                2
            } else {
                1
            };
            return Ok(Some(self.emit(ScannerToken::BondSpec(ch), width)));
        }

        // RDKit✔️✔️: \:\t\t\t{ return COLON_TOKEN; }
        // RDKit✔️✔️: \_\t\t\t{ return UNDERSCORE_TOKEN; }
        // RDKit✔️✔️: \#\t\t\t{ return HASH_TOKEN; }
        // RDKit✔️✔️: \-\t\t\t{ return MINUS_TOKEN; }
        // RDKit✔️✔️: \+\t\t\t{ return PLUS_TOKEN; }
        // RDKit✔️✔️: \{       \t{ return RANGE_OPEN_TOKEN; }
        // RDKit✔️✔️: \}       \t{ return RANGE_CLOSE_TOKEN; }
        // RDKit✔️✔️: \.       \t{ return SEPARATOR_TOKEN; }
        // RDKit✔️✔️: \%              { return PERCENT_TOKEN; }
        // RDKit✔️✔️: [0]\t\t{ yylval->ival = 0;  return ZERO_TOKEN; }
        // RDKit✔️✔️: [1-9]\t\t{ yylval->ival = yytext[0]-'0';  return NONZERO_DIGIT_TOKEN; }
        // RDKit✔️✔️: \!\t\t\t{ return NOT_TOKEN; }
        // RDKit✔️✔️: \;\t\t\t{ return SEMI_TOKEN; }
        // RDKit✔️✔️: \&\t\t\t{ return AND_TOKEN; }
        // RDKit✔️✔️: \,\t\t\t{ return OR_TOKEN; }
        let token = match ch {
            ':' => ScannerToken::Colon,
            '_' => ScannerToken::Underscore,
            '#' => ScannerToken::Hash,
            '-' => ScannerToken::Minus,
            '+' => ScannerToken::Plus,
            '{' => ScannerToken::RangeOpen,
            '}' => ScannerToken::RangeClose,
            '.' => ScannerToken::Separator,
            '%' => ScannerToken::Percent,
            '0'..='9' => ScannerToken::Digit(ch.to_digit(10).expect("ASCII digit") as u8),
            '!' => ScannerToken::Not,
            ';' => ScannerToken::Semi,
            '&' => ScannerToken::And,
            ',' => ScannerToken::Or,
            _ => return Ok(None),
        };
        Ok(Some(self.emit(token, 1)))
    }
}

// This is ordered longest-first to reproduce flex longest-match behavior.
const ELEMENT_SYMBOLS: &[&str] = &[
    "Uut", "Uup", "He", "Li", "Be", "Ne", "Na", "Mg", "Al", "Si", "Ar", "Ca", "Sc", "Ti", "Cr",
    "Mn", "Co", "Fe", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Kr", "Rb", "Sr", "Zr", "Nb", "Mo",
    "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "Re", "Os",
    "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "Np",
    "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Fl", "Lv", "K", "V", "Y", "W", "U",
];

/// Tokenize a molecule SMARTS with the sole stateful scanner.
///
/// Local complexity review: scanning and compaction are each linear in input
/// length. State operations are O(1); token storage is O(n). The implementation
/// creates the same O(n) input/token buffers as flex and does not rescan an
/// already compacted token range. Fixed element-rule lookup has constant size.
fn tokenize(input: &str) -> Result<Vec<(Token, usize)>, SmartsParseError> {
    generic_parse_helper(input, ScannerStart::Molecule)
}

fn generic_parse_helper(
    input: &str,
    start: ScannerStart,
) -> Result<Vec<(Token, usize)>, SmartsParseError> {
    // RDKit✔️✔️: template<int(*lex_init)(void**), size_t(*string_setup)(...),
    // RDKit✔️✔️: int generic_parse_helper(T parser, const std::string &inp, ...)
    // RDKit✔️✔️: TEST_ASSERT(!lex_init(&scanner));
    // RDKit✔️✔️: res = parser(inp.c_str() + ltrim, ... start_tok, ...);
    // RDKit✔️✔️: lex_destroy(scanner);
    // Local complexity review: scanner setup and token compaction each make
    // one linear pass and retain O(n) token storage. Rust owns the scanner
    // directly, so there is no FFI lifetime or destroy call, no second parse
    // path, and no rescan after compaction.
    let scanned = SmartsScanner::new(input, start).scan()?;
    compact_scanned_tokens(input, &scanned)
}

fn compact_scanned_tokens(
    input: &str,
    scanned: &[ScannedToken],
) -> Result<Vec<(Token, usize)>, SmartsParseError> {
    let chars: Vec<char> = input.chars().collect();
    let mut tokens = Vec::new();
    let mut i = 1usize;
    while i < scanned.len() {
        let current = &scanned[i];
        match &current.token {
            ScannerToken::Start(_) => unreachable!("start token is first"),
            ScannerToken::OrganicElement(symbol) | ScannerToken::AtomElement(symbol) => {
                tokens.push((Token::OrganicElement(symbol.clone()), current.start));
            }
            ScannerToken::AromaticElement(symbol) => {
                tokens.push((Token::AromaticElement(symbol.clone()), current.start));
            }
            ScannerToken::AtomOpen => {
                let content_start = current.end;
                let mut depth = 1usize;
                let mut cursor = i + 1;
                while cursor < scanned.len() && depth > 0 {
                    match scanned[cursor].token {
                        ScannerToken::AtomOpen => depth += 1,
                        ScannerToken::AtomClose => depth -= 1,
                        _ => {}
                    }
                    cursor += 1;
                }
                if depth != 0 {
                    return Err(SmartsParseError::UnclosedBracket(current.start));
                }
                let close = &scanned[cursor - 1];
                let content: String = chars[content_start..close.start].iter().collect();
                tokens.push((Token::BracketContent(content), current.start));
                i = cursor;
                continue;
            }
            ScannerToken::BondSpec(ch) => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(*ch)), current.start))
            }
            ScannerToken::DativeRight => {
                tokens.push((Token::BondSpec(BondLexeme::DativeRight), current.start));
            }
            ScannerToken::DativeLeft => {
                tokens.push((Token::BondSpec(BondLexeme::DativeLeft), current.start));
            }
            ScannerToken::At => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol('@')), current.start))
            }
            ScannerToken::Colon => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(':')), current.start))
            }
            ScannerToken::Hash => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol('#')), current.start))
            }
            ScannerToken::Minus => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol('-')), current.start))
            }
            ScannerToken::GroupOpen => tokens.push((Token::OpenParen, current.start)),
            ScannerToken::GroupClose => tokens.push((Token::CloseParen, current.start)),
            ScannerToken::Separator => tokens.push((Token::Dot, current.start)),
            ScannerToken::Digit(value) => {
                tokens.push((Token::RingClosureDigit(u32::from(*value)), current.start));
            }
            ScannerToken::Percent => {
                // RDKit✔️✔️: ring_number: digit
                // RDKit✔️✔️: | PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit
                // RDKit✔️✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit ... GROUP_CLOSE_TOKEN
                let (number, consumed) = compact_ring_number(scanned, i)?;
                tokens.push((Token::RingClosurePercent(number), current.start));
                i = consumed;
                continue;
            }
            ScannerToken::Not => tokens.push((Token::Not, current.start)),
            ScannerToken::Semi => tokens.push((Token::Semi, current.start)),
            ScannerToken::And => tokens.push((Token::And, current.start)),
            ScannerToken::Or => tokens.push((Token::Or, current.start)),
            ScannerToken::EndOfStream => tokens.push((Token::EndOfStream, current.start)),
            token => {
                return Err(SmartsParseError::UnexpectedCharacter {
                    position: current.start,
                    character: chars.get(current.start).copied().unwrap_or('?'),
                    context: format!("unexpected {token:?} token in molecule SMARTS"),
                });
            }
        }
        i += 1;
    }
    Ok(tokens)
}

fn invalid_percent(position: usize) -> SmartsParseError {
    SmartsParseError::UnexpectedCharacter {
        position,
        character: '%',
        context: "expected two digits after %".to_string(),
    }
}

fn compact_ring_number(
    scanned: &[ScannedToken],
    percent_index: usize,
) -> Result<(u32, usize), SmartsParseError> {
    let Some(next) = scanned.get(percent_index + 1) else {
        return Err(invalid_percent(scanned[percent_index].start));
    };
    let (digits_start, close_required) = match next.token {
        ScannerToken::Digit(value) if value != 0 => {
            let _ = value;
            (percent_index + 1, false)
        }
        ScannerToken::GroupOpen => (percent_index + 2, true),
        _ => return Err(invalid_percent(scanned[percent_index].start)),
    };
    let mut cursor = digits_start;
    let mut value = 0u32;
    let mut count = 0usize;
    while let Some(ScannedToken {
        token: ScannerToken::Digit(digit),
        ..
    }) = scanned.get(cursor)
    {
        value = value * 10 + u32::from(*digit);
        count += 1;
        cursor += 1;
        if !close_required && count == 2 {
            break;
        }
        if close_required && count == 5 {
            break;
        }
    }
    if (!close_required && count != 2)
        || (close_required
            && (count == 0
                || !matches!(
                    scanned.get(cursor).map(|token| &token.token),
                    Some(ScannerToken::GroupClose)
                )))
    {
        return Err(invalid_percent(scanned[percent_index].start));
    }
    let consumed = if close_required { cursor + 1 } else { cursor };
    Ok((value, consumed))
}

fn invalid_atom_operator(position: usize, operator: char) -> SmartsParseError {
    SmartsParseError::InvalidAtomPrimitive {
        position,
        detail: format!("operator '{operator}' has no left operand"),
    }
}

// ---------------------------------------------------------------------------
// Recursive-descent SMARTS Parser
// ---------------------------------------------------------------------------

/// RDKit❗✔️: Our recursive-descent SMARTS parser. In RDKit, the parser is
/// generated by bison from smarts.yy. We implement the same grammar logic
/// by hand.
struct SmartsParser<'a> {
    tokens: &'a [(Token, usize)],
    input: &'a str,
    pos: usize,
    /// Track ring closures: closure number to atom, query, direction, and input position.
    ring_closure_targets:
        BTreeMap<u32, (usize, QueryNode<BondQueryPredicate>, BondDirection, usize)>,
}

struct ParsedSmartsAtom {
    query: QueryNode<AtomQueryPredicate>,
    atom_map: Option<u32>,
}

fn split_atom_map_suffix(content: &str) -> Result<(&str, Option<u32>), SmartsParseError> {
    let bytes = content.as_bytes();
    let mut digit_start = bytes.len();
    while digit_start > 0 && bytes[digit_start - 1].is_ascii_digit() {
        digit_start -= 1;
    }
    if digit_start == bytes.len() || digit_start == 0 || bytes[digit_start - 1] != b':' {
        return Ok((content, None));
    }
    let colon = digit_start - 1;
    let atom_map = content[digit_start..].parse::<u32>().map_err(|_| {
        SmartsParseError::InvalidAtomPrimitive {
            position: digit_start,
            detail: "atom map number is out of range".to_string(),
        }
    })?;
    Ok((&content[..colon], Some(atom_map)))
}

impl<'a> SmartsParser<'a> {
    fn new(tokens: &'a [(Token, usize)], input: &'a str) -> Self {
        Self {
            tokens,
            input,
            pos: 0,
            ring_closure_targets: BTreeMap::new(),
        }
    }

    fn peek(&self) -> &(Token, usize) {
        &self.tokens[self.pos]
    }

    fn advance(&mut self) {
        self.pos += 1;
    }

    fn pos_info(&self) -> usize {
        self.tokens[self.pos].1
    }

    fn require_end(&self, context: &str) -> Result<(), SmartsParseError> {
        match self.peek() {
            (Token::EndOfStream, _) => Ok(()),
            (token, position) => Err(SmartsParseError::UnexpectedCharacter {
                position: *position,
                character: self
                    .input
                    .chars()
                    .nth(*position)
                    .unwrap_or_else(|| format!("{token:?}").chars().next().unwrap_or('?')),
                context: format!("unexpected trailing token in {context}"),
            }),
        }
    }

    /// Parse the full SMARTS pattern into a private graph.
    ///
    /// RDKit source: smarts.yy — the top-level grammar rule produces a molecule.
    fn parse_smarts_molecule(&mut self) -> Result<ParsedSmartsGraph, SmartsParseError> {
        // RDKit✔️✔️: mol: atomd {
        // RDKit✔️✔️:   int sz     = molList->size();
        // RDKit✔️✔️:   molList->resize( sz + 1);
        // RDKit✔️✔️:   (*molList)[ sz ] = new RWMol();
        // RDKit✔️✔️:   $1->setProp(RDKit::common_properties::_SmilesStart,1);
        // RDKit✔️✔️:   (*molList)[ sz ]->addAtom($1,true,true);
        // RDKit✔️✔️:   //delete $1;
        // RDKit✔️✔️:   $$ = sz;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | mol atomd       {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit✔️✔️:   int atomIdx1=a1->getIdx();
        // RDKit✔️✔️:   int atomIdx2=mp->addAtom($2,true,true);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
        // RDKit✔️✔️:   newB->setOwningMol(mp);
        // RDKit✔️✔️:   newB->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:   newB->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:   newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond(newB,true);
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol bond_expr atomd  {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2 = mp->addAtom($3,true,true);
        // RDKit✔️✔️:   if( $2->getBondType() == Bond::DATIVER ){
        // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   }else if ( $2->getBondType() == Bond::DATIVEL ){
        // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond($2);
        // RDKit✔️✔️:   delete $2;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol SEPARATOR_TOKEN atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
        // RDKit✔️✔️:   mp->addAtom($3,true,true);
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol ring_number {
        // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(atom, nullptr);
        // RDKit✔️✔️:   newB->setOwningMol(mp);
        // RDKit✔️✔️:   newB->setBeginAtomIdx(atom->getIdx());
        // RDKit✔️✔️:   mp->setBondBookmark(newB,$2);
        // RDKit✔️✔️:   if(!(mp->getAllBondsWithBookmark($2).size()%2)){
        // RDKit✔️✔️:     newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   mp->setAtomBookmark(atom,$2);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   INT_VECT tmp;
        // RDKit✔️✔️:   if(atom->hasProp(RDKit::common_properties::_RingClosures)){
        // RDKit✔️✔️:     atom->getProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   tmp.push_back(-($2+1));
        // RDKit✔️✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol bond_expr ring_number {
        // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
        // RDKit✔️✔️:
        // RDKit✔️✔️:   mp->setBondBookmark($2,$3);
        // RDKit✔️✔️:   $2->setOwningMol(mp);
        // RDKit✔️✔️:   $2->setBeginAtomIdx(atom->getIdx());
        // RDKit✔️✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->setAtomBookmark(atom,$3);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   SmilesParseOps::CheckRingClosureBranchStatus(atom,mp);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   INT_VECT tmp;
        // RDKit✔️✔️:   if(atom->hasProp(RDKit::common_properties::_RingClosures)){
        // RDKit✔️✔️:     atom->getProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   tmp.push_back(-($3+1));
        // RDKit✔️✔️:   atom->setProp(RDKit::common_properties::_RingClosures,tmp);
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol branch_open_token atomd {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   Atom *a1 = mp->getActiveAtom();
        // RDKit✔️✔️:   int atomIdx1=a1->getIdx();
        // RDKit✔️✔️:   int atomIdx2=mp->addAtom($3,true,true);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   QueryBond *newB = SmilesParseOps::getUnspecifiedQueryBond(a1,mp->getAtomWithIdx(atomIdx2));
        // RDKit✔️✔️:   newB->setOwningMol(mp);
        // RDKit✔️✔️:   newB->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:   newB->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:   newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond(newB);
        // RDKit✔️✔️:   delete newB;
        // RDKit✔️✔️:
        // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol branch_open_token bond_expr atomd  {
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   int atomIdx1 = mp->getActiveAtom()->getIdx();
        // RDKit✔️✔️:   int atomIdx2 = mp->addAtom($4,true,true);
        // RDKit✔️✔️:   if( $3->getBondType() == Bond::DATIVER ){
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $3->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   }else if ( $3->getBondType() == Bond::DATIVEL ){
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx2);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setBondType(Bond::DATIVE);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     $3->setBeginAtomIdx(atomIdx1);
        // RDKit✔️✔️:     $3->setEndAtomIdx(atomIdx2);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   $3->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️:   mp->addBond($3,true);
        // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
        // RDKit✔️✔️:
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️:
        // RDKit✔️✔️: | mol GROUP_CLOSE_TOKEN {
        // RDKit✔️✔️:   if(branchPoints.empty()){
        // RDKit✔️✔️:      yyerror(input,molList,branchPoints,scanner,start_token, current_token_position, "extra close parentheses");
        // RDKit✔️✔️:      yyErrorCleanup(molList);
        // RDKit✔️✔️:      YYABORT;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RWMol *mp = (*molList)[$$];
        // RDKit✔️✔️:   mp->setActiveAtom(branchPoints.back().first);
        // RDKit✔️✔️:   branchPoints.pop_back();
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: ;
        // RDKit✔️✔️:
        let mut atom_queries: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();
        let mut atom_maps: Vec<Option<u32>> = Vec::new();
        let mut bond_queries: Vec<QueryNode<BondQueryPredicate>> = Vec::new();
        let mut bond_directions: Vec<BondDirection> = Vec::new();
        let mut bond_edges: Vec<(usize, usize)> = Vec::new();
        let mut ring_closures: Vec<(u32, usize)> = Vec::new();
        let mut ring_closure_bonds: Vec<(u32, usize, QueryNode<BondQueryPredicate>)> = Vec::new();

        // RDKit✔️✔️: Parse the first atom
        let first = self.parse_atomd()?;
        atom_queries.push(first.query);
        atom_maps.push(first.atom_map);

        // RDKit✔️✔️: Parse the rest of the pattern
        let _ = self.parse_smarts_chain(
            &mut atom_queries,
            &mut atom_maps,
            &mut bond_queries,
            &mut bond_directions,
            &mut bond_edges,
            &mut ring_closures,
            &mut ring_closure_bonds,
            0,
        )?;

        if let Some(number) = self.ring_closure_targets.keys().next().copied() {
            return Err(SmartsParseError::UnbalancedRingClosure(number));
        }
        self.require_end("molecule SMARTS")?;

        Ok(ParsedSmartsGraph::new(
            atom_queries,
            atom_maps,
            bond_queries,
            bond_directions,
            bond_edges,
            ring_closures,
            ring_closure_bonds,
        ))
    }

    /// RDKit source: smarts.yy — mol → atom atom_list
    /// RDKit✔️✔️: Parse the chain of atoms, bonds, branches, and ring closures.
    fn parse_smarts_chain(
        &mut self,
        atom_queries: &mut Vec<QueryNode<AtomQueryPredicate>>,
        atom_maps: &mut Vec<Option<u32>>,
        bond_queries: &mut Vec<QueryNode<BondQueryPredicate>>,
        bond_directions: &mut Vec<BondDirection>,
        bond_edges: &mut Vec<(usize, usize)>,
        ring_closures: &mut Vec<(u32, usize)>,
        ring_closure_bonds: &mut Vec<(u32, usize, QueryNode<BondQueryPredicate>)>,
        mut active_atom_idx: usize,
    ) -> Result<usize, SmartsParseError> {
        loop {
            match self.peek() {
                (Token::EndOfStream, _) => break,
                (Token::CloseParen, _) => break,

                // Bond spec followed by atom
                (Token::BondSpec(_), _) | (Token::Not, _) | (Token::And, _) | (Token::Semi, _) => {
                    let direction = self.current_bond_direction();
                    let bond = self.parse_bond_expr()?;
                    match self.peek() {
                        (Token::RingClosureDigit(n), pos) | (Token::RingClosurePercent(n), pos) => {
                            let num = *n;
                            let bond_pos = *pos;
                            self.advance();
                            // RDKit source: smarts.yy lines 321-337
                            // RDKit✔️✔️: | mol bond_expr ring_number {
                            // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
                            // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
                            // RDKit✔️✔️:   mp->setBondBookmark($2,$3);
                            // RDKit✔️✔️:   $2->setOwningMol(mp);
                            // RDKit✔️✔️:   $2->setBeginAtomIdx(atom->getIdx());
                            // RDKit✔️✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
                            // RDKit✔️✔️:   mp->setAtomBookmark(atom,$3);
                            self.record_ring_closure(
                                num,
                                active_atom_idx,
                                bond,
                                direction,
                                bond_pos,
                                bond_queries,
                                bond_directions,
                                bond_edges,
                                ring_closures,
                                ring_closure_bonds,
                            );
                        }
                        _ => {
                            let (bond, reverse_endpoints) = normalize_dative_bond(bond);
                            let atom = self.parse_atomd()?;
                            atom_queries.push(atom.query);
                            atom_maps.push(atom.atom_map);
                            let end_atom_idx = atom_queries.len() - 1;
                            bond_queries.push(bond);
                            bond_directions.push(direction);
                            if reverse_endpoints {
                                bond_edges.push((end_atom_idx, active_atom_idx));
                            } else {
                                bond_edges.push((active_atom_idx, end_atom_idx));
                            }
                            active_atom_idx = end_atom_idx;
                        }
                    }
                }

                // No bond spec — implicit single/aromatic bond (SMARTS semantics)
                _ => {
                    // Check if next is a ring closure or branch first
                    match self.peek() {
                        (Token::RingClosureDigit(n), pos) | (Token::RingClosurePercent(n), pos) => {
                            let num = *n;
                            let bond_pos = *pos;
                            self.advance();
                            // Record ring closure on RDKit's current active atom.
                            self.record_ring_closure(
                                num,
                                active_atom_idx,
                                unspecified_smarts_bond_query(),
                                BondDirection::None,
                                bond_pos,
                                bond_queries,
                                bond_directions,
                                bond_edges,
                                ring_closures,
                                ring_closure_bonds,
                            );
                        }
                        (Token::OpenParen, _) => {
                            let _branch_position = self.parse_branch_open_token()?;
                            // RDKit✔️✔️: branchPoints.push_back({atomIdx1, $2});
                            // RDKit✔️✔️: GROUP_CLOSE_TOKEN restores the active atom
                            // RDKit✔️✔️: with mp->setActiveAtom(branchPoints.back().first).
                            let _branch_active = self.parse_smarts_chain(
                                atom_queries,
                                atom_maps,
                                bond_queries,
                                bond_directions,
                                bond_edges,
                                ring_closures,
                                ring_closure_bonds,
                                active_atom_idx,
                            )?;
                            match self.peek() {
                                (Token::CloseParen, _) => {
                                    self.advance();
                                }
                                (tok, pos) => {
                                    return Err(SmartsParseError::UnexpectedCharacter {
                                        position: *pos,
                                        character: format!("{:?}", tok)
                                            .chars()
                                            .next()
                                            .unwrap_or('?'),
                                        context: "expected close parenthesis".to_string(),
                                    });
                                }
                            }
                        }
                        (Token::Dot, _) => {
                            // RDKit✔️✔️: Dot separates disconnected fragments
                            self.advance();
                            let atom = self.parse_atomd()?;
                            atom_queries.push(atom.query);
                            atom_maps.push(atom.atom_map);
                            active_atom_idx = atom_queries.len() - 1;
                        }
                        // Atom follows implicitly with default bond
                        _ => {
                            bond_queries.push(unspecified_smarts_bond_query());
                            bond_directions.push(BondDirection::None);
                            let atom = self.parse_atomd()?;
                            atom_queries.push(atom.query);
                            atom_maps.push(atom.atom_map);
                            let end_atom_idx = atom_queries.len() - 1;
                            bond_edges.push((active_atom_idx, end_atom_idx));
                            active_atom_idx = end_atom_idx;
                        }
                    }
                }
            }
        }

        Ok(active_atom_idx)
    }

    fn record_ring_closure(
        &mut self,
        num: u32,
        atom_idx: usize,
        bond: QueryNode<BondQueryPredicate>,
        direction: BondDirection,
        bond_pos: usize,
        bond_queries: &mut Vec<QueryNode<BondQueryPredicate>>,
        bond_directions: &mut Vec<BondDirection>,
        bond_edges: &mut Vec<(usize, usize)>,
        ring_closures: &mut Vec<(u32, usize)>,
        ring_closure_bonds: &mut Vec<(u32, usize, QueryNode<BondQueryPredicate>)>,
    ) {
        ring_closures.push((num, atom_idx));
        ring_closure_bonds.push((num, atom_idx, bond.clone()));
        if let Some((open_atom_idx, open_bond, open_direction, _open_pos)) =
            self.ring_closure_targets.remove(&num)
        {
            let closing_is_unspecified = bond == unspecified_smarts_bond_query();
            let resolved_bond = if closing_is_unspecified {
                open_bond
            } else {
                bond
            };
            let resolved_direction = if direction != BondDirection::None {
                direction
            } else {
                open_direction
            };
            bond_queries.push(resolved_bond);
            bond_directions.push(resolved_direction);
            bond_edges.push((open_atom_idx, atom_idx));
        } else {
            self.ring_closure_targets
                .insert(num, (atom_idx, bond, direction, bond_pos));
        }
    }

    fn current_bond_direction(&self) -> BondDirection {
        // RDKit✔️✔️: [\\]{1,2} { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit✔️✔️:   yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
        // RDKit✔️✔️:   yylval->bond->setQuery(makeSingleOrAromaticBondQuery()); }
        // RDKit✔️✔️: [\/] { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit✔️✔️:   yylval->bond->setBondDir(Bond::ENDUPRIGHT);
        // RDKit✔️✔️:   yylval->bond->setQuery(makeSingleOrAromaticBondQuery()); }
        // RDKit's composite grammar retains the left QueryBond and expands its
        // query, so direction comes from the first primitive in the expression.
        self.tokens[self.pos..]
            .iter()
            .find_map(|(token, _)| match token {
                Token::BondSpec(BondLexeme::Symbol('/')) => Some(BondDirection::EndUpRight),
                Token::BondSpec(BondLexeme::Symbol('\\')) => Some(BondDirection::EndDownRight),
                Token::BondSpec(_) => Some(BondDirection::None),
                Token::Not => None,
                _ => Some(BondDirection::None),
            })
            .unwrap_or(BondDirection::None)
    }

    /// Parse the `atomd` production through the sole typed-query parser.
    ///
    /// Local complexity review: token dispatch is O(1). Bracket content is
    /// parsed once in O(n) time with O(n) query storage; no branch reparses or
    /// clones the SMARTS input.
    fn parse_atomd(&mut self) -> Result<ParsedSmartsAtom, SmartsParseError> {
        // RDKit✔️✔️: atomd:\tsimple_atom
        // RDKit✔️✔️: | hydrogen_atom
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN atom_expr ATOM_CLOSE_TOKEN
        // RDKit✔️✔️: {
        // RDKit✔️✔️:   $$ = $2;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN atom_expr COLON_TOKEN number ATOM_CLOSE_TOKEN
        // RDKit✔️✔️: {
        // RDKit✔️✔️:   $$ = $2;
        // RDKit✔️✔️:   $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: ;
        let (token, _pos) = self.peek().clone();
        match token {
            Token::OrganicElement(name) | Token::AromaticElement(name) => {
                let query = parse_simple_atom(&name).ok_or_else(|| {
                    SmartsParseError::InvalidAtomPrimitive {
                        position: self.pos_info(),
                        detail: format!("invalid simple atom '{name}'"),
                    }
                })?;
                self.advance();
                Ok(ParsedSmartsAtom {
                    query,
                    atom_map: None,
                })
            }
            Token::BracketContent(content) => {
                self.advance();
                self.parse_bracket_atom_content(&content)
            }
            Token::EndOfStream => Err(SmartsParseError::UnexpectedEnd(
                "expected atom but reached end".to_string(),
            )),
            _ => {
                let pos = self.pos_info();
                Err(SmartsParseError::UnexpectedCharacter {
                    position: pos,
                    character: '?',
                    context: "expected atom expression".to_string(),
                })
            }
        }
    }

    /// Parse the content inside a bracket atom: e.g. "C", "C@@H", "N+", "O-", "#6", "6X4"
    ///
    /// RDKit source: smarts.yy — the ATOM_TOKEN production and its associated actions
    /// RDKit✔️✔️: Bracket atom content is parsed as a sequence of primitives AND-ed together.
    fn parse_bracket_atom_content(
        &mut self,
        content: &str,
    ) -> Result<ParsedSmartsAtom, SmartsParseError> {
        // RDKit✔️✔️: atom_expr: atom_expr AND_TOKEN atom_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | atom_expr OR_TOKEN atom_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
        // RDKit✔️✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit✔️✔️:   $1->setAtomicNum(0);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | atom_expr SEMI_TOKEN atom_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | atom_expr point_query {
        // RDKit✔️✔️:   atom_expr_and_point_query($1, $2);
        // RDKit✔️✔️:   delete $2;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | atom_expr AND_TOKEN point_query {
        // RDKit✔️✔️:   atom_expr_and_point_query($1, $3);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | point_query
        // RDKit✔️✔️: ;
        //
        // COSMolKit stores only typed query predicates here, so RDKit's
        // mutable QueryAtom chemical-property cleanup and chiral-tag copying
        // have no separate state to migrate: chirality remains a predicate in
        // the composed tree. Local complexity review: the expression is read
        // once in O(n) time and stored in O(n) query nodes. Finalization moves
        // vectors without rescanning or cloning their children.
        let (content, atom_map) = split_atom_map_suffix(content)?;
        let chars: Vec<char> = content.chars().collect();
        let len = chars.len();
        if len == 0 {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: 0,
                detail: "empty atom expression".to_string(),
            });
        }
        if let Some(query) = self.try_parse_hydrogen_atom(&chars, len)? {
            return Ok(ParsedSmartsAtom { query, atom_map });
        }
        let mut i = 0;
        let mut needs_operand = true;
        let mut clauses: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();
        let mut current_or_terms: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();
        let mut current_term: Vec<QueryNode<AtomQueryPredicate>> = Vec::new();

        fn finalize_term(
            current_term: &mut Vec<QueryNode<AtomQueryPredicate>>,
            current_or_terms: &mut Vec<QueryNode<AtomQueryPredicate>>,
        ) {
            if current_term.is_empty() {
                return;
            }
            let mut terms = std::mem::take(current_term).into_iter();
            let mut term = terms.next().expect("nonempty atom-query term");
            for query in terms {
                super::query::query_atom_expand_query(
                    &mut term,
                    query,
                    CompositeQueryType::And,
                    true,
                );
            }
            current_or_terms.push(term);
        }

        fn finalize_clause(
            current_term: &mut Vec<QueryNode<AtomQueryPredicate>>,
            current_or_terms: &mut Vec<QueryNode<AtomQueryPredicate>>,
            clauses: &mut Vec<QueryNode<AtomQueryPredicate>>,
        ) {
            finalize_term(current_term, current_or_terms);
            if current_or_terms.is_empty() {
                return;
            }
            let mut terms = std::mem::take(current_or_terms).into_iter();
            let mut clause = terms.next().expect("nonempty atom-query clause");
            for query in terms {
                super::query::query_atom_expand_query(
                    &mut clause,
                    query,
                    CompositeQueryType::Or,
                    true,
                );
            }
            clauses.push(clause);
        }

        while i < len {
            let ch = chars[i];

            // Handle logical OR (comma)
            if ch == ',' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                finalize_term(&mut current_term, &mut current_or_terms);
                needs_operand = true;
                i += 1;
                continue;
            }

            // Handle logical AND (ampersand)
            if ch == '&' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                needs_operand = true;
                i += 1;
                continue;
            }

            // Handle semicolon (AND)
            if ch == ';' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                finalize_clause(&mut current_term, &mut current_or_terms, &mut clauses);
                needs_operand = true;
                i += 1;
                continue;
            }

            let (pred, consumed) = self.parse_point_query(&chars, i, len)?;

            current_term.push(pred);
            needs_operand = false;

            i = consumed;
        }

        if needs_operand {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: len.saturating_sub(1),
                detail: "atom expression ends without an operand".to_string(),
            });
        }

        finalize_clause(&mut current_term, &mut current_or_terms, &mut clauses);

        // Combine predicates
        // RDKit source: smarts.yy precedence gives implicit/`&` high-precedence
        // AND inside each comma term, comma OR inside a clause, and `;`
        // low-precedence AND across clauses. Every reduction uses the canonical
        // QueryAtom::expandQuery port so AtomNull algebra is preserved.
        let mut clauses = clauses.into_iter();
        let mut query = clauses.next().expect("at least one bracket clause");
        for clause in clauses {
            super::query::query_atom_expand_query(
                &mut query,
                clause,
                CompositeQueryType::And,
                true,
            );
        }
        Ok(ParsedSmartsAtom { query, atom_map })
    }

    fn parse_point_query(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        // RDKit✔️✔️: point_query: NOT_TOKEN point_query {
        // RDKit✔️✔️:   $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
        // RDKit✔️✔️:   $2->setAtomicNum(0);
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($2);
        // RDKit✔️✔️:   $$ = $2;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | recursive_query
        // RDKit✔️✔️: | atom_query
        // RDKit✔️✔️: ;
        //
        // Typed predicates contain no separate mutable atomic number or
        // chemical-property cache, so negation wraps the query node directly.
        // Local complexity review: each leading NOT is consumed once and the
        // selected recursive/atom query is parsed once, for O(n) time and
        // O(number of effective NOT nodes) storage.
        let mut pos = start;
        let mut negate = false;
        while chars.get(pos) == Some(&'!') {
            negate = !negate;
            pos += 1;
        }
        if pos >= len {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: "NOT has no point query".to_string(),
            });
        }
        let (query, consumed) = self.parse_atom_primitive(chars, pos, len)?;
        Ok((if negate { QueryNode::not(query) } else { query }, consumed))
    }

    fn try_parse_hydrogen_atom(
        &self,
        chars: &[char],
        len: usize,
    ) -> Result<Option<QueryNode<AtomQueryPredicate>>, SmartsParseError> {
        // RDKit✔️✔️: hydrogen_atom:\tATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN
        // RDKit✔️✔️: {
        // RDKit✔️✔️:   $$ = new QueryAtom(1);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN
        // RDKit✔️✔️: {
        // RDKit✔️✔️:   $$ = new QueryAtom(1);
        // RDKit✔️✔️:   $$->setProp(RDKit::common_properties::molAtomMapNumber,$4);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN number H_TOKEN ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setIsotope($2);
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN number H_TOKEN COLON_TOKEN number ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setIsotope($2);
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   newQ->setProp(RDKit::common_properties::molAtomMapNumber,$5);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️:
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN H_TOKEN charge_spec ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setFormalCharge($3);
        // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomFormalChargeQuery($3),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setFormalCharge($3);
        // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomFormalChargeQuery($3),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   newQ->setProp(RDKit::common_properties::molAtomMapNumber,$5);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN number H_TOKEN charge_spec ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setIsotope($2);
        // RDKit✔️✔️:   newQ->setFormalCharge($4);
        // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomFormalChargeQuery($4),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_OPEN_TOKEN number H_TOKEN charge_spec COLON_TOKEN number ATOM_CLOSE_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom(1);
        // RDKit✔️✔️:   newQ->setIsotope($2);
        // RDKit✔️✔️:   newQ->setFormalCharge($4);
        // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomIsotopeQuery($2),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   newQ->expandQuery(makeAtomFormalChargeQuery($4),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   newQ->setProp(RDKit::common_properties::molAtomMapNumber,$6);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: ;
        //
        // Local complexity review: one left-to-right pass over bracket content,
        // O(n) time and O(1) parser state. The resulting query has at most four
        // nodes, matching RDKit's bounded expansion and allocation behavior.
        if len == 0 {
            return Ok(None);
        }
        let mut pos = 0usize;
        let mut isotope = None;
        if chars[pos].is_ascii_digit() {
            let (num, consumed) = self.parse_number(chars, pos, len)?;
            isotope = Some(num as u16);
            pos = consumed;
        }
        if pos >= len || chars[pos] != 'H' {
            return Ok(None);
        }
        pos += 1;
        if pos < len && chars[pos].is_ascii_digit() {
            return Ok(None);
        }

        let mut formal_charge = None;
        if pos < len && matches!(chars[pos], '+' | '-') {
            let (pred, consumed) = self.parse_atom_primitive(chars, pos, len)?;
            match pred {
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(charge)) => {
                    formal_charge = Some(charge);
                    pos = consumed;
                }
                _ => return Ok(None),
            }
        }

        if pos != len {
            return Ok(None);
        }

        let mut clauses = vec![QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1))];
        if let Some(isotope) = isotope {
            clauses.push(super::query::make_atom_isotope_query(isotope));
        }
        if let Some(formal_charge) = formal_charge {
            clauses.push(super::query::make_atom_formal_charge_query(formal_charge));
        }
        Ok(Some(if clauses.len() == 1 {
            clauses.pop().expect("single hydrogen atom clause")
        } else {
            QueryNode::And(clauses)
        }))
    }

    /// Parse a single atom primitive from the bracket content starting at position `i`.
    ///
    /// RDKit source: smarts.yy — primitives within ATOM_TOKEN
    /// RDKit✔️✔️: Handles all SMARTS atom primitives.
    fn parse_atom_primitive(
        &self,
        chars: &[char],
        i: usize,
        len: usize,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        // RDKit✔️✔️: atom_query:\tsimple_atom
        // RDKit✔️✔️: | number simple_atom {
        // RDKit✔️✔️:   $2->setIsotope($1);
        // RDKit✔️✔️:   $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   $$=$2;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ATOM_TOKEN
        // RDKit✔️✔️: | number ATOM_TOKEN {
        // RDKit✔️✔️:   $2->setIsotope($1);
        // RDKit✔️✔️:   $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   $$=$2;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | HASH_TOKEN number { $$ = new QueryAtom($2); }
        // RDKit✔️✔️: | number HASH_TOKEN number {
        // RDKit✔️✔️:   $$ = new QueryAtom($3);
        // RDKit✔️✔️:   $$->setIsotope($1);
        // RDKit✔️✔️:   $$->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | HETERONEIGHBOR_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | MIN_RINGSIZE_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | RINGSIZE_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | RINGBOND_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | IMPLICIT_H_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | HETERONEIGHBOR_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomNumHeteroatomNbrsQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | MIN_RINGSIZE_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomMinRingSizeQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGSIZE_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomInRingOfSizeQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGBOND_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomRingBondCountQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | IMPLICIT_H_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomImplicitHCountQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        //
        // Local complexity review: every non-recursive branch consumes one
        // primitive and an optional decimal in O(token length), constructs a
        // bounded number of typed nodes, and never rescans earlier input.
        if i >= len {
            return Err(SmartsParseError::UnexpectedEnd(
                "expected atom primitive".to_string(),
            ));
        }

        let ch = chars[i];

        if let Some(range_query) = self.parse_possible_range_query(chars, i, len)? {
            return Ok(range_query);
        }

        // Atomic number: #N
        // RDKit✔️✔️: smarts.yy — HASH_TOKEN NUMBER
        if ch == '#' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(num as u8)),
                consumed,
            ));
        }

        // Recursive SMARTS: $(...)
        if ch == '$' {
            return self.parse_recursive_query(chars, i, len);
        }

        if let Some(charge) = self.parse_charge_spec(chars, i, len)? {
            return Ok(charge);
        }

        // RDKit✔️✔️: | AT_TOKEN AT_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomNullQuery());
        // RDKit✔️✔️:   newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | AT_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomNullQuery());
        // RDKit✔️✔️:   newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | CHI_CLASS_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomNullQuery());
        // RDKit✔️✔️:   newQ->setChiralTag($1);
        // RDKit✔️✔️:   newQ->setProp(common_properties::_chiralPermutation,0);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | CHI_CLASS_TOKEN number {
        // RDKit✔️✔️:   if($2==0){
        // RDKit✔️✔️:     yyerror(input,molList,branchPoints,scanner,start_token, current_token_position,
        // RDKit✔️✔️:             "chiral permutation cannot be zero");
        // RDKit✔️✔️:     yyErrorCleanup(molList);
        // RDKit✔️✔️:     YYABORT;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomNullQuery());
        // RDKit✔️✔️:   newQ->setChiralTag($1);
        // RDKit✔️✔️:   newQ->setProp(common_properties::_chiralPermutation,$2);
        // RDKit✔️✔️:   $$=newQ;
        // RDKit✔️✔️: }
        if ch == '@' {
            let start = i + 1;
            if start < len && chars[start] == '@' {
                // @@
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(
                        crate::ChiralTag::TetrahedralCw,
                    )),
                    start + 1,
                ));
            }
            let class = if start + 1 < len {
                Some(chars[start..=start + 1].iter().collect::<String>())
            } else {
                None
            };
            if let Some((tag, class_end)) = class.as_deref().and_then(|class| {
                let tag = match class {
                    "TH" => crate::ChiralTag::Tetrahedral,
                    "AL" => crate::ChiralTag::Allene,
                    "SP" => crate::ChiralTag::SquarePlanar,
                    "TB" => crate::ChiralTag::TrigonalBipyramidal,
                    "OH" => crate::ChiralTag::Octahedral,
                    _ => return None,
                };
                Some((tag, start + 2))
            }) {
                let (permutation, consumed) = self.parse_optional_number(chars, class_end, len);
                if consumed != class_end && permutation == 0 {
                    return Err(SmartsParseError::InvalidAtomPrimitive {
                        position: class_end,
                        detail: "chiral permutation cannot be zero".to_string(),
                    });
                }
                return Ok((
                    QueryNode::And(vec![
                        QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(tag)),
                        QueryNode::Predicate(AtomQueryPredicate::ChiralPermutationMatch(
                            permutation,
                        )),
                    ]),
                    consumed,
                ));
            }
            // @
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(
                    crate::ChiralTag::TetrahedralCcw,
                )),
                start,
            ));
        }

        // Element symbol
        // RDKit✔️✔️: <IN_ATOM_STATE>Hg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Tl |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pb |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bi |
        // RDKit✔️✔️: <IN_ATOM_STATE>Po |
        // RDKit✔️✔️: <IN_ATOM_STATE>At |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ra |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ac |
        // RDKit✔️✔️: <IN_ATOM_STATE>Th |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pa |
        // RDKit✔️✔️: <IN_ATOM_STATE>U |
        // RDKit✔️✔️: <IN_ATOM_STATE>Np |
        // RDKit✔️✔️: <IN_ATOM_STATE>Pu |
        // RDKit✔️✔️: <IN_ATOM_STATE>Am |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bk |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cf |
        // RDKit✔️✔️: <IN_ATOM_STATE>Es |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fm |
        // RDKit✔️✔️: <IN_ATOM_STATE>Md |
        // RDKit✔️✔️: <IN_ATOM_STATE>No |
        // RDKit✔️✔️: <IN_ATOM_STATE>Lr |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rf |
        // RDKit✔️✔️: <IN_ATOM_STATE>Db |
        // RDKit✔️✔️: <IN_ATOM_STATE>Sg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Bh |
        // RDKit✔️✔️: <IN_ATOM_STATE>Hs |
        // RDKit✔️✔️: <IN_ATOM_STATE>Mt |
        // RDKit✔️✔️: <IN_ATOM_STATE>Ds |
        // RDKit✔️✔️: <IN_ATOM_STATE>Rg |
        // RDKit✔️✔️: <IN_ATOM_STATE>Cn |
        // RDKit✔️✔️: <IN_ATOM_STATE>Uut |
        // RDKit✔️✔️: <IN_ATOM_STATE>Fl |
        // RDKit✔️✔️: <IN_ATOM_STATE>Uup |
        // RDKit✔️✔️: <IN_ATOM_STATE>Lv  { yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
        // RDKit✔️✔️:                       return ATOM_TOKEN;
        // RDKit✔️✔️:                    }
        //
        // RDKit flex selects the longest matching token in IN_ATOM_STATE, so
        // two-letter elements such as Hg must be consumed before the H_TOKEN
        // hydrogen-count rule below.
        if ch.is_ascii_uppercase() {
            let start = i;
            let end = i + 1;
            if end < len && chars[end].is_ascii_lowercase() {
                let two_char: String = chars[start..=end].iter().collect();
                if let Some(query) = parse_atom_token(&two_char) {
                    return Ok((query, end + 1));
                }
            }
            if ch != 'H' {
                let one_char: String = chars[start..end].iter().collect();
                if let Some(query) = parse_atom_token(&one_char) {
                    return Ok((query, end));
                }
            }
        }

        // Hydrogen-count SMARTS queries: `h` or `h<N>`, `H` or `H<N>`
        // RDKit✔️✔️: smarts.ll / smarts.yy split `h` into AtomHasImplicitH /
        // RDKit✔️✔️: AtomImplicitHCount and `H` into AtomHCount.
        if ch == 'h' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((super::query::make_atom_has_implicit_h_query(), consumed));
            }
            return Ok((
                super::query::make_atom_implicit_h_count_query(num as u8),
                consumed,
            ));
        }
        if ch == 'H' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((super::query::make_atom_h_count_query(1), consumed));
            }
            return Ok((super::query::make_atom_h_count_query(num as u8), consumed));
        }

        // Ring membership: R or R<N>
        // RDKit✔️✔️: smarts.yy — R_TOKEN (optional NUMBER)
        if ch == 'R' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::InRing), consumed));
            }
            if num == 0 {
                // RDKit✔️✔️: <IN_ATOM_STATE>R {
                // RDKit✔️✔️:   yylval->atom = new QueryAtom();
                // RDKit✔️✔️:   yylval->atom->setQuery(new AtomRingQuery(-1));
                // RDKit✔️✔️:   return COMPLEX_ATOM_QUERY_TOKEN;
                // RDKit✔️✔️: }
                //
                // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
                // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
                // RDKit✔️✔️:   $$ = $1;
                // RDKit✔️✔️: }
                //
                // AtomRingQuery(0) is equivalent to "not in any ring".
                return Ok((
                    QueryNode::not(QueryNode::Predicate(AtomQueryPredicate::InRing)),
                    consumed,
                ));
            }
            return Ok((
                // RDKit✔️✔️: <IN_ATOM_STATE>R {
                // RDKit✔️✔️:   yylval->atom = new QueryAtom();
                // RDKit✔️✔️:   yylval->atom->setQuery(new AtomRingQuery(-1));
                // RDKit✔️✔️:   return COMPLEX_ATOM_QUERY_TOKEN;
                // RDKit✔️✔️: }
                // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
                // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
                // RDKit✔️✔️:   $$ = $1;
                // RDKit✔️✔️: }
                QueryNode::Predicate(AtomQueryPredicate::NumAtomRings(num as i32)),
                consumed,
            ));
        }
        if ch == 'r' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::InRing), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::SmallestRingSize(num as u8)),
                consumed,
            ));
        }
        if ch == 'k' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((QueryNode::Predicate(AtomQueryPredicate::InRing), consumed));
            }
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::InRingOfSize(num as u8)),
                consumed,
            ));
        }

        // Connectivity/degree: X or X<N>
        // RDKit✔️✔️: smarts.yy — X_TOKEN (optional NUMBER)
        if ch == 'X' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                super::query::make_atom_total_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num as u8
                }),
                consumed,
            ));
        }

        // Non-hydrogen degree: d or d<N>
        if ch == 'd' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                super::query::make_atom_non_hydrogen_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num
                }),
                consumed,
            ));
        }

        // Heteroatom-neighbor queries: z/Z and their explicit counts.
        if ch == 'z' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                if consumed == i + 1 {
                    super::query::make_atom_has_heteroatom_nbrs_query()
                } else {
                    super::query::make_atom_num_heteroatom_nbrs_query(num as u8)
                },
                consumed,
            ));
        }
        if ch == 'Z' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                if consumed == i + 1 {
                    super::query::make_atom_has_aliphatic_heteroatom_nbrs_query()
                } else {
                    super::query::make_atom_num_aliphatic_heteroatom_nbrs_query(num as u8)
                },
                consumed,
            ));
        }

        // Ring connectivity: x or x<N>
        // RDKit✔️✔️: <IN_ATOM_STATE>x {
        // RDKit✔️✔️:   yylval->atom = new QueryAtom();
        // RDKit✔️✔️:   yylval->atom->setQuery(makeAtomHasRingBondQuery());
        // RDKit✔️✔️:   return RINGBOND_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGBOND_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomRingBondCountQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        if ch == 'x' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            if consumed == i + 1 {
                return Ok((super::query::make_atom_has_ring_bond_query(), consumed));
            }
            return Ok((
                super::query::make_atom_ring_bond_count_query(num as u8),
                consumed,
            ));
        }

        // Degree: D or D<N>
        // RDKit✔️✔️: smarts.yy — D_TOKEN (optional NUMBER)
        if ch == 'D' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                super::query::make_atom_explicit_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num as u8
                }),
                consumed,
            ));
        }

        // Hybridization: ^1, ^2, ^3, ...
        // RDKit✔️✔️: caret hybridization primitives map to explicit hybridization predicates.
        if ch == '^' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            let hybridization = match num {
                0 => crate::Hybridization::S,
                1 => crate::Hybridization::Sp,
                2 => crate::Hybridization::Sp2,
                3 => crate::Hybridization::Sp3,
                4 => crate::Hybridization::Sp3d,
                5 => crate::Hybridization::Sp3d2,
                _ => {
                    return Err(SmartsParseError::InvalidAtomPrimitive {
                        position: i,
                        detail: "hybridization must be in ^0..^5".to_string(),
                    });
                }
            };
            return Ok((
                super::query::make_atom_hybridization_query(hybridization),
                consumed,
            ));
        }

        // Valence: v or v<N>
        // RDKit✔️✔️: smarts.yy — v_TOKEN (optional NUMBER)
        if ch == 'v' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len);
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::TotalValence(if consumed == i + 1 {
                    1
                } else {
                    num as u8
                })),
                consumed,
            ));
        }

        // Aromatic query: a / A
        // RDKit✔️✔️: smarts.yy — a_TOKEN / A_TOKEN
        if ch == 'a' {
            return Ok((
                parse_simple_atom("a").expect("aromatic simple query token"),
                i + 1,
            ));
        }
        if ch == 'A' {
            return Ok((
                parse_simple_atom("A").expect("aliphatic simple query token"),
                i + 1,
            ));
        }

        // Unsaturated: u
        // RDKit✔️✔️: smarts.yy — u_TOKEN
        if ch == 'u' {
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::IsUnsaturated),
                i + 1,
            ));
        }

        // RDKit's `number simple_atom`, `number ATOM_TOKEN`, and
        // `number HASH_TOKEN number` reductions retain the atom query as the
        // left child and append the isotope query. Other numeric prefixes
        // remain standalone isotope primitives and are combined by atom_expr
        // in input order.
        if ch.is_ascii_digit() {
            let (num, consumed) = self.parse_number(chars, i, len)?;
            let isotope =
                u16::try_from(num).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: "isotope is out of range".to_string(),
                })?;

            // RDKit✔️✔️: | number H_TOKEN {
            // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit✔️✔️:   newQ->setQuery(makeAtomIsotopeQuery($1));
            // RDKit✔️✔️:   newQ->setIsotope($1);
            // RDKit✔️✔️:   newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
            // RDKit✔️✔️:   newQ->setNumExplicitHs(1);
            // RDKit✔️✔️:   newQ->setNoImplicit(true);
            // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_H_MASK;
            // RDKit✔️✔️:   $$=newQ;
            // RDKit✔️✔️: }
            // RDKit✔️✔️: | number H_TOKEN number {
            // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit✔️✔️:   newQ->setQuery(makeAtomIsotopeQuery($1));
            // RDKit✔️✔️:   newQ->setIsotope($1);
            // RDKit✔️✔️:   newQ->expandQuery(makeAtomHCountQuery($3),Queries::COMPOSITE_AND,true);
            // RDKit✔️✔️:   newQ->setNumExplicitHs($3);
            // RDKit✔️✔️:   newQ->setNoImplicit(true);
            // RDKit✔️✔️:   newQ->getFlags() |= SMARTS_H_MASK;
            // RDKit✔️✔️:   $$=newQ;
            // RDKit✔️✔️: }
            // Local complexity review: both implementations parse the isotope
            // and optional H count once, then allocate one two-leaf query.
            if chars.get(consumed) == Some(&'H') {
                let (hydrogen_count, end) = self.parse_optional_number(chars, consumed + 1, len);
                let hydrogen_count = if end == consumed + 1 {
                    1
                } else {
                    u8::try_from(hydrogen_count).map_err(|_| {
                        SmartsParseError::InvalidAtomPrimitive {
                            position: consumed + 1,
                            detail: "hydrogen count is out of range".to_string(),
                        }
                    })?
                };
                let mut query = super::query::make_atom_isotope_query(isotope);
                super::query::query_atom_expand_query(
                    &mut query,
                    super::query::make_atom_h_count_query(hydrogen_count),
                    CompositeQueryType::And,
                    true,
                );
                return Ok((query, end));
            }

            let mut atom_and_end = None;
            if chars.get(consumed) == Some(&'#') {
                let (atomic_number, end) = self.parse_number(chars, consumed + 1, len)?;
                let atomic_number = u8::try_from(atomic_number).map_err(|_| {
                    SmartsParseError::InvalidAtomPrimitive {
                        position: consumed + 1,
                        detail: "atomic number is out of range".to_string(),
                    }
                })?;
                atom_and_end = Some((
                    QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
                    end,
                ));
            } else if let Some(next) = chars.get(consumed) {
                let mut end = consumed + 1;
                if chars.get(end).is_some_and(char::is_ascii_lowercase) {
                    let two_char = chars[consumed..=end].iter().collect::<String>();
                    if parse_atom_token(&two_char).is_some() {
                        end += 1;
                    }
                }
                let symbol = chars[consumed..end].iter().collect::<String>();
                if *next == '*' || parse_atom_token(&symbol).is_some() {
                    atom_and_end = parse_atom_token(&symbol).map(|atom| (atom, end));
                }
            }
            if let Some((mut atom, end)) = atom_and_end {
                super::query::query_atom_expand_query(
                    &mut atom,
                    super::query::make_atom_isotope_query(isotope),
                    CompositeQueryType::And,
                    true,
                );
                return Ok((atom, end));
            }
            return Ok((super::query::make_atom_isotope_query(isotope), consumed));
        }

        // Lowercase aromatic element inside bracket (e.g. [c], [se])
        if ch.is_ascii_lowercase() && ch != 'a' && ch != 'u' && ch != 'v' && ch != 'r' && ch != 'h'
        {
            let mut consumed = i + 1;
            if i + 1 < len {
                let two_char: String = chars[i..=i + 1].iter().collect();
                match two_char.as_str() {
                    // RDKit✔️✔️: <IN_ATOM_STATE>si	{  yylval->ival = 14;  return AROMATIC_ATOM_TOKEN;  }
                    // RDKit✔️✔️: <IN_ATOM_STATE>as	{  yylval->ival = 33;  return AROMATIC_ATOM_TOKEN;  }
                    // RDKit✔️✔️: <IN_ATOM_STATE>se	{  yylval->ival = 34;  return AROMATIC_ATOM_TOKEN;  }
                    // RDKit✔️✔️: <IN_ATOM_STATE>te	{  yylval->ival = 52;  return AROMATIC_ATOM_TOKEN;  }
                    "si" | "as" | "se" | "te" => {
                        consumed = i + 2;
                    }
                    _ => {}
                }
            }
            let name: String = chars[i..consumed].iter().collect();
            let query =
                parse_simple_atom(&name).ok_or_else(|| SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: format!("invalid aromatic simple atom '{name}'"),
                })?;
            return Ok((query, consumed));
        }

        // Wildcard inside a bracket uses the same SIMPLE_ATOM_QUERY_TOKEN path.
        if ch == '*' {
            return Ok((
                parse_simple_atom("*").expect("wildcard simple query token"),
                i + 1,
            ));
        }

        Err(SmartsParseError::InvalidAtomPrimitive {
            position: i,
            detail: format!("unexpected character '{}'", ch),
        })
    }

    fn parse_charge_spec(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<Option<(QueryNode<AtomQueryPredicate>, usize)>, SmartsParseError> {
        // RDKit✔️✔️: charge_spec: PLUS_TOKEN PLUS_TOKEN { $$=2; }
        // RDKit✔️✔️: | PLUS_TOKEN number { $$=$2; }
        // RDKit✔️✔️: | PLUS_TOKEN { $$=1; }
        // RDKit✔️✔️: | MINUS_TOKEN MINUS_TOKEN { $$=-2; }
        // RDKit✔️✔️: | MINUS_TOKEN number { $$=-$2; }
        // RDKit✔️✔️: | MINUS_TOKEN { $$=-1; }
        // Local complexity review: one sign, one optional repeated sign, and
        // one decimal scan are consumed once. The helper is O(number length)
        // with O(1) auxiliary state and creates one typed charge leaf, matching
        // the bounded bison reduction without rescans or temporary collections.
        let Some(sign) = chars.get(start).copied() else {
            return Ok(None);
        };
        if !matches!(sign, '+' | '-') {
            return Ok(None);
        }
        let next = start + 1;
        if chars.get(next) == Some(&sign) {
            let charge = if sign == '+' { 2 } else { -2 };
            return Ok(Some((
                super::query::make_atom_formal_charge_query(charge),
                next + 1,
            )));
        }
        let (magnitude, consumed) = self.parse_optional_number(chars, next, len);
        let magnitude = if consumed == next { 1 } else { magnitude };
        let signed = if sign == '+' {
            i64::from(magnitude)
        } else {
            -i64::from(magnitude)
        };
        let charge = i8::try_from(signed).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
            position: start,
            detail: format!(
                "formal charge {signed} is outside COSMolKit's modeled i8 charge range"
            ),
        })?;
        Ok(Some((
            super::query::make_atom_formal_charge_query(charge),
            consumed,
        )))
    }

    fn parse_possible_range_query(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<Option<(QueryNode<AtomQueryPredicate>, usize)>, SmartsParseError> {
        // RDKit✔️✔️: possible_range_query : COMPLEX_ATOM_QUERY_TOKEN
        // RDKit✔️✔️: | HETERONEIGHBOR_ATOM_QUERY_TOKEN {
        // RDKit✔️✔️:   $1->setQuery(makeAtomNumHeteroatomNbrsQuery(0));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN {
        // RDKit✔️✔️:   $1->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery(0));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | MIN_RINGSIZE_ATOM_QUERY_TOKEN {
        // RDKit✔️✔️:   $1->setQuery(makeAtomMinRingSizeQuery(5)); // this is going to be ignored anyway
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGBOND_ATOM_QUERY_TOKEN {
        // RDKit✔️✔️:   $1->setQuery(makeAtomRingBondCountQuery(0));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | IMPLICIT_H_ATOM_QUERY_TOKEN {
        // RDKit✔️✔️:   $1->setQuery(makeAtomImplicitHCountQuery(0));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | PLUS_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomFormalChargeQuery(0));
        // RDKit✔️✔️:   $$ = newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | MINUS_TOKEN {
        // RDKit✔️✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit✔️✔️:   newQ->setQuery(makeAtomNegativeFormalChargeQuery(0));
        // RDKit✔️✔️:   $$ = newQ;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: ;
        // RDKit✔️✔️: | possible_range_query RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
        // RDKit✔️✔️:   ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>($4,oq->getDataFunc(),
        // RDKit✔️✔️:     std::string("greater_")+oq->getDescription());
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
        // RDKit✔️✔️:   ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>($3,oq->getDataFunc(),
        // RDKit✔️✔️:     std::string("less_")+oq->getDescription());
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | possible_range_query RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   ATOM_EQUALS_QUERY *oq = static_cast<ATOM_EQUALS_QUERY *>($1->getQuery());
        // RDKit✔️✔️:   ATOM_RANGE_QUERY *nq = makeAtomRangeQuery($3,$5,false,false,
        // RDKit✔️✔️:     oq->getDataFunc(),
        // RDKit✔️✔️:     std::string("range_")+oq->getDescription());
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: /* "k" queries have to be handled differently */
        // RDKit✔️✔️: | RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN MINUS_TOKEN number RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   int lv = -1;
        // RDKit✔️✔️:   int uv = $4;
        // RDKit✔️✔️:   ATOM_GREATEREQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_GREATEREQUAL_QUERY>(uv,[lv,uv](Atom const *at) {
        // RDKit✔️✔️:             return queryAtomIsInRingOfSize(at, lv, uv);
        // RDKit✔️✔️:           },std::string("greater_AtomRingSize"));
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN number MINUS_TOKEN RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   int lv = $3;
        // RDKit✔️✔️:   int uv = -1;
        // RDKit✔️✔️:   ATOM_LESSEQUAL_QUERY *nq = makeAtomSimpleQuery<ATOM_LESSEQUAL_QUERY>(lv,[lv,uv](Atom const *at) {
        // RDKit✔️✔️:             return queryAtomIsInRingOfSize(at, lv, uv);
        // RDKit✔️✔️:           },std::string("less_AtomRingSize"));
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGSIZE_ATOM_QUERY_TOKEN RANGE_OPEN_TOKEN number MINUS_TOKEN number RANGE_CLOSE_TOKEN {
        // RDKit✔️✔️:   int lv = $3;
        // RDKit✔️✔️:   int uv = $5;
        // RDKit✔️✔️:   ATOM_RANGE_QUERY *nq = makeAtomRangeQuery(lv,uv,false,false,[lv,uv](Atom const *at) {
        // RDKit✔️✔️:             return queryAtomIsInRingOfSize(at, lv, uv);
        // RDKit✔️✔️:           },std::string("range_AtomRingSize"));
        // RDKit✔️✔️:   $1->setQuery(nq);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        //
        // Local complexity review: the parser selects one data-function tag,
        // consumes each range character once, and creates one bounded typed
        // leaf. This is O(range-token length) time and O(1) auxiliary space,
        // matching bison's bounded reductions without rescanning, cloning,
        // keyed lookup, or a second query representation.
        let data_function = match chars.get(start) {
            Some('D') => AtomRangeDataFunction::ExplicitDegree,
            Some('d') => AtomRangeDataFunction::NonHydrogenDegree,
            Some('X') => AtomRangeDataFunction::TotalDegree,
            Some('v') => AtomRangeDataFunction::TotalValence,
            Some('R') => AtomRangeDataFunction::NumAtomRings,
            Some('z') => AtomRangeDataFunction::NumHeteroatomNeighbors,
            Some('Z') => AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors,
            Some('r') => AtomRangeDataFunction::MinRingSize,
            Some('x') => AtomRangeDataFunction::RingBondCount,
            Some('h') => AtomRangeDataFunction::ImplicitHydrogenCount,
            Some('+') => AtomRangeDataFunction::FormalCharge,
            Some('-') => AtomRangeDataFunction::NegativeFormalCharge,
            Some('k') => AtomRangeDataFunction::AtomRingSize {
                lower: 0,
                upper: 0,
                lower_open: false,
                upper_open: false,
            },
            _ => return Ok(None),
        };
        if chars.get(start + 1) != Some(&'{') {
            return Ok(None);
        }

        let (lower, upper, consumed) = self.parse_possible_range_bounds(chars, start + 2, len)?;
        let data_function = if matches!(data_function, AtomRangeDataFunction::AtomRingSize { .. }) {
            AtomRangeDataFunction::AtomRingSize {
                lower: lower.unwrap_or(-1),
                upper: upper.unwrap_or(-1),
                lower_open: false,
                upper_open: false,
            }
        } else {
            data_function
        };
        let query =
            make_atom_possible_range_query(lower, upper, data_function).ok_or_else(|| {
                SmartsParseError::InvalidAtomPrimitive {
                    position: start,
                    detail: "empty atom range".to_string(),
                }
            })?;
        Ok(Some((query, consumed)))
    }

    fn parse_possible_range_bounds(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<(Option<i32>, Option<i32>, usize), SmartsParseError> {
        let mut pos = start;
        let lower = if chars.get(pos) == Some(&'-') {
            None
        } else {
            let (value, consumed) = self.parse_number(chars, pos, len)?;
            pos = consumed;
            Some(
                i32::try_from(value).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
                    position: start,
                    detail: "atom range bound is out of range".to_string(),
                })?,
            )
        };
        if chars.get(pos) != Some(&'-') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '-' in atom range".to_string(),
            });
        }
        pos += 1;
        let upper = if chars.get(pos).is_some_and(char::is_ascii_digit) {
            let (value, consumed) = self.parse_number(chars, pos, len)?;
            pos = consumed;
            Some(
                i32::try_from(value).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
                    position: start,
                    detail: "atom range bound is out of range".to_string(),
                })?,
            )
        } else {
            None
        };
        if lower.is_none() && upper.is_none() {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "empty atom range".to_string(),
            });
        }
        if chars.get(pos) != Some(&'}') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '}' to close atom range".to_string(),
            });
        }
        Ok((lower, upper, pos + 1))
    }

    fn parse_recursive_query(
        &self,
        chars: &[char],
        start: usize,
        len: usize,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        // RDKit✔️✔️: recursive_query: BEGIN_RECURSE mol END_RECURSE {
        // RDKit✔️✔️:   // this is a recursive SMARTS expression
        // RDKit✔️✔️:   QueryAtom *qA = new QueryAtom();
        // RDKit✔️✔️:   //  FIX: there's maybe a leak here
        // RDKit✔️✔️:   RWMol *molP = (*molList)[$2];
        // RDKit✔️✔️:   // close any rings in the molecule:
        // RDKit✔️✔️:   SmilesParseOps::CloseMolRings(molP,0);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   //molP->debugMol(std::cout);
        // RDKit✔️✔️:   qA->setQuery(new RecursiveStructureQuery(molP));
        // RDKit✔️✔️:   //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
        // RDKit✔️✔️:   int sz = molList->size();
        // RDKit✔️✔️:   if ( sz==$2+1) {
        // RDKit✔️✔️:     molList->resize( sz-1 );
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   $$ = qA;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | BEGIN_RECURSE mol END_RECURSE UNDERSCORE_TOKEN  nonzero_number{
        // RDKit✔️✔️:   // UNDOCUMENTED EXTENSION:
        // RDKit✔️✔️:   // this is a recursive SMARTS expression with a serial number
        // RDKit✔️✔️:   // please don't write your own SMARTS that include this extension:
        // RDKit✔️✔️:   // the RDKit smarts parsing code will automatically insert serial
        // RDKit✔️✔️:   // numbers for recursive smarts patterns.
        // RDKit✔️✔️:   QueryAtom *qA = new QueryAtom();
        // RDKit✔️✔️:   //  FIX: there's maybe a leak here
        // RDKit✔️✔️:   RWMol *molP = (*molList)[$2];
        // RDKit✔️✔️:   // close any rings in the molecule:
        // RDKit✔️✔️:   SmilesParseOps::CloseMolRings(molP,0);
        // RDKit✔️✔️:
        // RDKit✔️✔️:   //molP->debugMol(std::cout);
        // RDKit✔️✔️:   qA->setQuery(new RecursiveStructureQuery(molP,$5));
        // RDKit✔️✔️:   //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
        // RDKit✔️✔️:   int sz = molList->size();
        // RDKit✔️✔️:   if ( sz==$2+1) {
        // RDKit✔️✔️:     molList->resize( sz-1 );
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   $$ = qA;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: ;
        //
        // Local complexity review: delimiter discovery and compilation each
        // scan the recursive SMARTS once, O(n), and graph construction is
        // O(V+E). The compiled query molecule is stored on the predicate, so
        // matching does not reparse the source string.
        if chars.get(start + 1) != Some(&'(') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: "expected '(' after '$'".to_string(),
            });
        }
        let mut depth = 1usize;
        let mut end = start + 2;
        while end < len && depth > 0 {
            match chars[end] {
                '(' => depth += 1,
                ')' => depth -= 1,
                _ => {}
            }
            end += 1;
        }
        if depth != 0 {
            return Err(SmartsParseError::UnclosedParenthesis(start + 1));
        }
        if end == start + 3 {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: "recursive SMARTS requires a molecule".to_string(),
            });
        }

        let recursive_smarts: String = chars[start..end].iter().collect();
        let mut consumed = end;
        let mut serial_number = 0;
        if chars.get(consumed) == Some(&'_') {
            consumed += 1;
            let serial_start = consumed;
            if !chars
                .get(consumed)
                .is_some_and(|digit| matches!(*digit, '1'..='9'))
            {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: consumed,
                    detail: "recursive SMARTS serial number must be nonzero".to_string(),
                });
            }
            while consumed < len && chars[consumed].is_ascii_digit() {
                consumed += 1;
            }
            serial_number = chars[serial_start..consumed]
                .iter()
                .collect::<String>()
                .parse()
                .map_err(|_| SmartsParseError::InvalidAtomPrimitive {
                    position: serial_start,
                    detail: "recursive SMARTS serial number is out of range".to_string(),
                })?;
        }
        let inner = recursive_smarts
            .strip_prefix("$(")
            .and_then(|value| value.strip_suffix(')'))
            .expect("balanced recursive SMARTS has delimiters");
        let query_mol = to_mol(inner).map_err(|detail| SmartsParseError::InvalidAtomPrimitive {
            position: start,
            detail,
        })?;
        Ok((
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(
                super::query::RecursiveStructureQuery::from_smarts(
                    recursive_smarts,
                    query_mol,
                    serial_number,
                ),
            )),
            consumed,
        ))
    }

    /// Parse a number from position i. Returns (value, consumed_index).
    fn parse_number(
        &self,
        chars: &[char],
        i: usize,
        len: usize,
    ) -> Result<(u32, usize), SmartsParseError> {
        // RDKit✔️✔️: number: ZERO_TOKEN | nonzero_number
        // RDKit✔️✔️: nonzero_number: NONZERO_DIGIT_TOKEN
        // RDKit✔️✔️: | nonzero_number digit { ... number too large ... }
        // RDKit✔️✔️: digit: NONZERO_DIGIT_TOKEN | ZERO_TOKEN
        // Local complexity review: one left-to-right digit fold is O(n) time
        // and O(1) state, matching the source reduction without reparsing.
        if i >= len || !chars[i].is_ascii_digit() {
            return Err(SmartsParseError::UnexpectedEnd(
                "expected number".to_string(),
            ));
        }
        let mut val = 0u32;
        let mut pos = i;
        while pos < len && chars[pos].is_ascii_digit() {
            let digit = chars[pos].to_digit(10).expect("ASCII digit");
            val = val
                .checked_mul(10)
                .and_then(|value| value.checked_add(digit))
                .ok_or_else(|| SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: "number too large".to_string(),
                })?;
            if val > i32::MAX as u32 {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: "number too large".to_string(),
                });
            }
            pos += 1;
        }
        Ok((val, pos))
    }

    /// Parse an optional number from position i. Returns (value if present else 0, consumed_index).
    fn parse_optional_number(&self, chars: &[char], i: usize, len: usize) -> (u32, usize) {
        // RDKit✔️✔️: nonzero_number: NONZERO_DIGIT_TOKEN ... digit
        // RDKit✔️✔️: digit: NONZERO_DIGIT_TOKEN | ZERO_TOKEN
        // Local complexity review: the optional fold is one linear pass with
        // constant state and no allocation; callers decide whether absence is
        // the grammar's implicit default.
        if i >= len || !chars[i].is_ascii_digit() {
            return (0, i);
        }
        let mut val = 0u32;
        let mut pos = i;
        while pos < len && chars[pos].is_ascii_digit() {
            val = val * 10 + chars[pos].to_digit(10).unwrap();
            pos += 1;
        }
        (val, pos)
    }

    /// Parse `bond_expr` with bison's declared operator precedence.
    fn parse_bond_expr(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        // RDKit✔️✔️: bond_expr:bond_expr AND_TOKEN bond_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | bond_expr OR_TOKEN bond_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | bond_expr SEMI_TOKEN bond_expr {
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   delete $3;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | bond_query
        // RDKit✔️✔️: ;
        // Local complexity review: each token is consumed once through four
        // fixed precedence levels. Every source expandQuery maps to one O(1)
        // typed query combination with the same child order and null-query
        // algebra. There is no rescan, clone, temporary token collection, or
        // second bond-expression parser; time and query storage are O(n).
        self.parse_bond_semi_expr()
    }

    fn parse_branch_open_token(&mut self) -> Result<usize, SmartsParseError> {
        // RDKit✔️✔️: branch_open_token: GROUP_OPEN_TOKEN { $$ = current_token_position; };
        // Local complexity review: this helper performs one token check and
        // advance in O(1), returning the original source position without a
        // scan, allocation, or duplicate branch parser.
        let (Token::OpenParen, position) = self.peek() else {
            return Err(SmartsParseError::UnexpectedCharacter {
                position: self.pos_info(),
                character: '?',
                context: "expected branch opening token".to_string(),
            });
        };
        let position = *position;
        self.advance();
        Ok(position)
    }

    fn parse_bond_semi_expr(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        let mut query = self.parse_bond_or_expr()?;
        while matches!(self.peek(), (Token::Semi, _)) {
            self.advance();
            let rhs = self.parse_bond_or_expr()?;
            query_bond_expand_query(&mut query, rhs, CompositeQueryType::And, true);
        }
        Ok(query)
    }

    fn parse_bond_or_expr(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        let mut query = self.parse_bond_and_expr()?;
        while matches!(self.peek(), (Token::Or, _)) {
            self.advance();
            let rhs = self.parse_bond_and_expr()?;
            query_bond_expand_query(&mut query, rhs, CompositeQueryType::Or, true);
        }
        Ok(query)
    }

    fn parse_bond_and_expr(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        let mut query = self.parse_bond_query()?;
        while matches!(self.peek(), (Token::And, _)) {
            self.advance();
            let rhs = self.parse_bond_query()?;
            query_bond_expand_query(&mut query, rhs, CompositeQueryType::And, true);
        }
        Ok(query)
    }

    fn parse_bond_query(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        // RDKit✔️✔️: bond_query: bondd
        // RDKit✔️✔️: | bond_query bondd {
        // RDKit✔️✔️:   $1->expandQuery($2->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   delete $2;
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: ;
        // Local complexity review: each adjacent primitive is consumed once
        // and contributes one O(1) ordered expandQuery operation. The loop is
        // O(n) time and query storage, with no token rescan, subtree clone,
        // lookup, temporary collection, or alternate bond-query decoder.
        let mut query = self.parse_bondd()?;
        while matches!(self.peek(), (Token::BondSpec(_), _) | (Token::Not, _)) {
            let rhs = self.parse_bondd()?;
            query_bond_expand_query(&mut query, rhs, CompositeQueryType::And, true);
        }
        Ok(query)
    }

    fn parse_bondd(&mut self) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
        // RDKit✔️❌: bondd: BOND_TOKEN
        // RDKit✔️❌: | MINUS_TOKEN {
        // RDKit✔️❌:   QueryBond *newB= new QueryBond();
        // RDKit✔️❌:   newB->setBondType(Bond::SINGLE);
        // RDKit✔️❌:   newB->setQuery(makeBondOrderEqualsQuery(Bond::SINGLE));
        // RDKit✔️❌:   $$ = newB;
        // RDKit✔️❌: }
        // RDKit✔️❌: | HASH_TOKEN {
        // RDKit✔️❌:   QueryBond *newB= new QueryBond();
        // RDKit✔️❌:   newB->setBondType(Bond::TRIPLE);
        // RDKit✔️❌:   newB->setQuery(makeBondOrderEqualsQuery(Bond::TRIPLE));
        // RDKit✔️❌:   $$ = newB;
        // RDKit✔️❌: }
        // RDKit✔️❌: | COLON_TOKEN {
        // RDKit✔️❌:   QueryBond *newB= new QueryBond();
        // RDKit✔️❌:   newB->setBondType(Bond::AROMATIC);
        // RDKit✔️❌:   newB->setQuery(makeBondOrderEqualsQuery(Bond::AROMATIC));
        // RDKit✔️❌:   $$ = newB;
        // RDKit✔️❌: }
        // RDKit✔️❌: | AT_TOKEN {
        // RDKit✔️❌:   QueryBond *newB= new QueryBond();
        // RDKit✔️❌:   newB->setQuery(makeBondIsInRingQuery());
        // RDKit✔️❌:   $$ = newB;
        // RDKit✔️❌: }
        // RDKit✔️❌: | NOT_TOKEN bondd {
        // RDKit✔️❌:   $2->getQuery()->setNegation(!($2->getQuery()->getNegation()));
        // RDKit✔️❌:   $$ = $2;
        // RDKit✔️❌: }
        // RDKit✔️❌: ;
        // Local complexity review: primitive dispatch and each NOT level are
        // O(1), so a run of NOT tokens is O(n) time and O(n) parser stack,
        // matching the right-recursive grammar without scans or clones. Typed
        // primitives avoid RDKit's QueryBond allocation, but an effective Rust
        // negation uses one Box where RDKit flips an inline bool; repeated NOT
        // may allocate/free that Box, so allocation behavior is materially
        // worse even though the final query semantics are exact.
        match self.peek() {
            (Token::Not, _) => {
                self.advance();
                let mut query = self.parse_bondd()?;
                let is_negated = matches!(&query, QueryNode::Not(_));
                query.set_negation(!is_negated);
                Ok(query)
            }
            (Token::BondSpec(lexeme), _) => {
                let query = bond_spec_to_query(*lexeme);
                self.advance();
                Ok(query)
            }
            (_, position) => Err(SmartsParseError::UnexpectedCharacter {
                position: *position,
                character: self.input.chars().nth(*position).unwrap_or('?'),
                context: "expected bond query primitive".to_string(),
            }),
        }
    }
}

// ---------------------------------------------------------------------------
// SMARTS primitive helpers
// ---------------------------------------------------------------------------

/// Reduce RDKit's `simple_atom` production into the sole typed atom-query leaf.
fn parse_simple_atom(name: &str) -> Option<QueryNode<AtomQueryPredicate>> {
    // RDKit✔️✔️: simple_atom: 	ORGANIC_ATOM_TOKEN {
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // This construction (and some others) may seem odd, but the
    // RDKit✔️✔️:   // SMARTS definition requires that an atom which is aliphatic on
    // RDKit✔️✔️:   // input (i.e. something in the "organic subset" that is given with
    // RDKit✔️✔️:   // a capital letter) only match aliphatic atoms.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   // The following rule applies a similar logic to aromatic atoms.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   $$ = new QueryAtom($1);
    // RDKit✔️✔️:   $$->setQuery(makeAtomTypeQuery($1,false));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | AROMATIC_ATOM_TOKEN {
    // RDKit✔️✔️:   $$ = new QueryAtom($1);
    // RDKit✔️✔️:   $$->setIsAromatic(true);
    // RDKit✔️✔️:   $$->setQuery(makeAtomTypeQuery($1,true));
    // RDKit✔️✔️: }
    // RDKit✔️✔️: | SIMPLE_ATOM_QUERY_TOKEN
    // RDKit✔️✔️: ;
    // Local complexity review: symbol dispatch is over a fixed-size set and
    // constructs one inline typed leaf in O(1) time/space. This removes the
    // source QueryAtom/query-object allocations without adding traversal,
    // cloning, temporary collections, or a second simple-atom decoder.
    let atom_type = |atomic_number, aromatic| {
        QueryNode::Predicate(AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        })
    };
    Some(match name {
        "B" => atom_type(5, false),
        "C" => atom_type(6, false),
        "N" => atom_type(7, false),
        "O" => atom_type(8, false),
        "F" => atom_type(9, false),
        "P" => atom_type(15, false),
        "S" => atom_type(16, false),
        "Cl" => atom_type(17, false),
        "Br" => atom_type(35, false),
        "I" => atom_type(53, false),
        "b" => atom_type(5, true),
        "c" => atom_type(6, true),
        "n" => atom_type(7, true),
        "o" => atom_type(8, true),
        "p" => atom_type(15, true),
        "s" => atom_type(16, true),
        "si" => atom_type(14, true),
        "as" => atom_type(33, true),
        "se" => atom_type(34, true),
        "te" => atom_type(52, true),
        "*" => make_atom_null_query(),
        "a" => super::query::make_atom_aromatic_query(),
        "A" => super::query::make_atom_aliphatic_query(),
        _ => return None,
    })
}

/// Decode either RDKit's `simple_atom` production or an `ATOM_TOKEN` emitted
/// by the bracket-atom lexer.
fn parse_atom_token(name: &str) -> Option<QueryNode<AtomQueryPredicate>> {
    parse_simple_atom(name).or_else(|| {
        element_symbol_to_atomic_number(name).map(|atomic_number| {
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_number))
        })
    })
}

/// Convert a bond specifier character to a bond query predicate.
///
/// RDKit source: smarts.yy / smarts.ll bond handling.
/// RDKit✔️✔️: `-`, `=`, `#`, `:`, `~`, `@`, and dative tokens become their
/// RDKit✔️✔️: typed bond predicates. SMARTS `/` and `\\` additionally store
/// RDKit✔️✔️: directional state on RDKit's QueryBond, while their graph query is
/// RDKit✔️✔️: the same `SingleOrAromaticBond` predicate as an unspecified bond.
fn bond_spec_to_query(lexeme: BondLexeme) -> QueryNode<BondQueryPredicate> {
    match lexeme {
        BondLexeme::DativeRight => make_bond_order_equals_query(BondOrder::DativeRight),
        BondLexeme::DativeLeft => make_bond_order_equals_query(BondOrder::DativeLeft),
        BondLexeme::Symbol(ch) => match ch {
            '-' => make_bond_order_equals_query(BondOrder::Single),
            '=' => make_bond_order_equals_query(BondOrder::Double),
            '#' => make_bond_order_equals_query(BondOrder::Triple),
            ':' => make_bond_order_equals_query(BondOrder::Aromatic),
            '$' => make_bond_order_equals_query(BondOrder::Quadruple),
            '@' => make_bond_is_in_ring_query(),
            '~' => make_bond_null_query(),
            '/' | '\\' => unspecified_smarts_bond_query(),
            _ => QueryNode::Predicate(BondQueryPredicate::Any),
        },
    }
}

fn unspecified_smarts_bond_query() -> QueryNode<BondQueryPredicate> {
    crate::smiles::get_unspecified_query_bond(false, None).1
}

fn normalize_dative_bond(
    query: QueryNode<BondQueryPredicate>,
) -> (QueryNode<BondQueryPredicate>, bool) {
    // RDKit✔️✔️: if( $2->getBondType() == Bond::DATIVER ){
    // RDKit✔️✔️:   $2->setBeginAtomIdx(atomIdx1);
    // RDKit✔️✔️:   $2->setEndAtomIdx(atomIdx2);
    // RDKit✔️✔️:   $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️: }else if ( $2->getBondType() == Bond::DATIVEL ){
    // RDKit✔️✔️:   $2->setBeginAtomIdx(atomIdx2);
    // RDKit✔️✔️:   $2->setEndAtomIdx(atomIdx1);
    // RDKit✔️✔️:   $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️: }
    match query {
        QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::DativeRight)) => {
            (make_bond_order_equals_query(BondOrder::Dative), false)
        }
        QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::DativeLeft)) => {
            (make_bond_order_equals_query(BondOrder::Dative), true)
        }
        query => (query, false),
    }
}

/// Look up the atomic number for an element symbol.
///
/// RDKit✔️✔️: Standard periodic table mapping.
fn element_symbol_to_atomic_number(symbol: &str) -> Option<u8> {
    match symbol {
        "H" => Some(1),
        "He" => Some(2),
        "Li" => Some(3),
        "Be" => Some(4),
        "B" => Some(5),
        "C" => Some(6),
        "N" => Some(7),
        "O" => Some(8),
        "F" => Some(9),
        "Ne" => Some(10),
        "Na" => Some(11),
        "Mg" => Some(12),
        "Al" => Some(13),
        "Si" => Some(14),
        "P" => Some(15),
        "S" => Some(16),
        "Cl" => Some(17),
        "Ar" => Some(18),
        "K" => Some(19),
        "Ca" => Some(20),
        "Sc" => Some(21),
        "Ti" => Some(22),
        "V" => Some(23),
        "Cr" => Some(24),
        "Mn" => Some(25),
        "Fe" => Some(26),
        "Co" => Some(27),
        "Ni" => Some(28),
        "Cu" => Some(29),
        "Zn" => Some(30),
        "Ga" => Some(31),
        "Ge" => Some(32),
        "As" => Some(33),
        "Se" => Some(34),
        "Br" => Some(35),
        "Kr" => Some(36),
        "Rb" => Some(37),
        "Sr" => Some(38),
        "Y" => Some(39),
        "Zr" => Some(40),
        "Nb" => Some(41),
        "Mo" => Some(42),
        "Tc" => Some(43),
        "Ru" => Some(44),
        "Rh" => Some(45),
        "Pd" => Some(46),
        "Ag" => Some(47),
        "Cd" => Some(48),
        "In" => Some(49),
        "Sn" => Some(50),
        "Sb" => Some(51),
        "Te" => Some(52),
        "I" => Some(53),
        "Xe" => Some(54),
        "Cs" => Some(55),
        "Ba" => Some(56),
        "La" => Some(57),
        "Ce" => Some(58),
        "Pr" => Some(59),
        "Nd" => Some(60),
        "Pm" => Some(61),
        "Sm" => Some(62),
        "Eu" => Some(63),
        "Gd" => Some(64),
        "Tb" => Some(65),
        "Dy" => Some(66),
        "Ho" => Some(67),
        "Er" => Some(68),
        "Tm" => Some(69),
        "Yb" => Some(70),
        "Lu" => Some(71),
        "Hf" => Some(72),
        "Ta" => Some(73),
        "W" => Some(74),
        "Re" => Some(75),
        "Os" => Some(76),
        "Ir" => Some(77),
        "Pt" => Some(78),
        "Au" => Some(79),
        "Hg" => Some(80),
        "Tl" => Some(81),
        "Pb" => Some(82),
        "Bi" => Some(83),
        "Po" => Some(84),
        "At" => Some(85),
        "Rn" => Some(86),
        "Fr" => Some(87),
        "Ra" => Some(88),
        "Ac" => Some(89),
        "Th" => Some(90),
        "Pa" => Some(91),
        "U" => Some(92),
        "Np" => Some(93),
        "Pu" => Some(94),
        "Am" => Some(95),
        "Cm" => Some(96),
        "Bk" => Some(97),
        "Cf" => Some(98),
        "Es" => Some(99),
        "Fm" => Some(100),
        "Md" => Some(101),
        "No" => Some(102),
        "Lr" => Some(103),
        "Rf" => Some(104),
        "Db" => Some(105),
        "Sg" => Some(106),
        "Bh" => Some(107),
        "Hs" => Some(108),
        "Mt" => Some(109),
        "Ds" => Some(110),
        "Rg" => Some(111),
        "Cn" => Some(112),
        "Nh" => Some(113),
        "Fl" => Some(114),
        "Mc" => Some(115),
        "Lv" => Some(116),
        "Ts" => Some(117),
        "Og" => Some(118),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tokenize_simple_smarts() {
        // CC — two organic carbons
        let tokens = tokenize("CC").unwrap();
        assert_eq!(tokens.len(), 3); // C, C, EOS
        assert_eq!(tokens[0].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[1].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[2].0, Token::EndOfStream);
    }

    #[test]
    fn test_tokenize_bracket_atom() {
        let tokens = tokenize("[N+]").unwrap();
        assert_eq!(tokens.len(), 2); // BracketContent, EOS
        match &tokens[0].0 {
            Token::BracketContent(content) => {
                assert_eq!(content, "N+");
            }
            _ => panic!("expected bracket content"),
        }
    }

    #[test]
    fn test_tokenize_ring_closure() {
        let tokens = tokenize("C1CC1").unwrap();
        assert_eq!(tokens.len(), 6); // C, 1, C, C, 1, EOS
        assert_eq!(tokens[0].0, Token::OrganicElement("C".to_string()));
        assert_eq!(tokens[1].0, Token::RingClosureDigit(1));
        assert_eq!(tokens[4].0, Token::RingClosureDigit(1));
    }

    #[test]
    fn test_tokenize_percent_ring_closure() {
        let tokens = tokenize("C%10CC%10").unwrap();
        assert_eq!(tokens.len(), 6);
        assert_eq!(tokens[1].0, Token::RingClosurePercent(10));
        assert_eq!(tokens[4].0, Token::RingClosurePercent(10));
    }

    #[test]
    fn test_tokenize_bond_specs() {
        let tokens = tokenize("C=O").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec(BondLexeme::Symbol('=')));
        let tokens = tokenize("C#N").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec(BondLexeme::Symbol('#')));
        let tokens = tokenize("C~C").unwrap();
        assert_eq!(tokens[1].0, Token::BondSpec(BondLexeme::Symbol('~')));
    }

    fn scanner_tokens(input: &str, start: ScannerStart) -> Vec<ScannerToken> {
        SmartsScanner::new(input, start)
            .scan()
            .expect("SMARTS scanner input")
            .into_iter()
            .map(|token| token.token)
            .collect()
    }

    #[test]
    fn smarts_lexer_initial() {
        assert_eq!(
            scanner_tokens("C.N", ScannerStart::Molecule),
            vec![
                ScannerToken::Start(ScannerStart::Molecule),
                ScannerToken::OrganicElement("C".to_string()),
                ScannerToken::Separator,
                ScannerToken::OrganicElement("N".to_string()),
                ScannerToken::EndOfStream,
            ]
        );
        assert_eq!(
            scanner_tokens("-", ScannerStart::Bond)[0],
            ScannerToken::Start(ScannerStart::Bond)
        );
    }

    #[test]
    fn smarts_lexer_atom_state() {
        let tokens = scanner_tokens("[FeDdxXvzZhRrk]", ScannerStart::Atom);
        assert!(tokens.contains(&ScannerToken::AtomElement("Fe".to_string())));
        for primitive in ['D', 'd', 'x', 'X', 'v', 'z', 'Z', 'h', 'R', 'r', 'k'] {
            assert!(tokens.contains(&ScannerToken::AtomPrimitive(primitive)));
        }
        assert!(matches!(
            tokens.first(),
            Some(ScannerToken::Start(ScannerStart::Atom))
        ));
    }

    #[test]
    fn smarts_lexer_branch_state() {
        assert_eq!(
            scanner_tokens("C(=O)N", ScannerStart::Molecule),
            vec![
                ScannerToken::Start(ScannerStart::Molecule),
                ScannerToken::OrganicElement("C".to_string()),
                ScannerToken::GroupOpen,
                ScannerToken::BondSpec('='),
                ScannerToken::OrganicElement("O".to_string()),
                ScannerToken::GroupClose,
                ScannerToken::OrganicElement("N".to_string()),
                ScannerToken::EndOfStream,
            ]
        );
    }

    #[test]
    fn smarts_lexer_recursion_state() {
        let tokens = scanner_tokens("[$(C(=O)N)_100]", ScannerStart::Atom);
        assert!(tokens.contains(&ScannerToken::BeginRecurse));
        assert!(tokens.contains(&ScannerToken::EndRecurse));
        assert_eq!(
            tokens
                .iter()
                .filter(|t| **t == ScannerToken::GroupOpen)
                .count(),
            1
        );
        assert_eq!(
            tokens
                .iter()
                .filter(|t| **t == ScannerToken::GroupClose)
                .count(),
            1
        );
    }

    #[test]
    fn smarts_lexer_elements() {
        let tokens = scanner_tokens("[Uut,C,si,as,se,te]", ScannerStart::Atom);
        assert!(tokens.contains(&ScannerToken::AtomElement("Uut".to_string())));
        assert!(tokens.contains(&ScannerToken::OrganicElement("C".to_string())));
        for aromatic in ["si", "as", "se", "te"] {
            assert!(tokens.contains(&ScannerToken::AromaticElement(aromatic.to_string())));
        }
        assert_eq!(
            scanner_tokens("Brcn", ScannerStart::Molecule)[1..5],
            [
                ScannerToken::OrganicElement("Br".to_string()),
                ScannerToken::AromaticElement("c".to_string()),
                ScannerToken::AromaticElement("n".to_string()),
                ScannerToken::EndOfStream,
            ]
        );
    }

    #[test]
    fn smarts_lexer_atom_primitives() {
        let tokens = scanner_tokens("[*aA;!&,H:+-#_]", ScannerStart::Atom);
        assert!(tokens.contains(&ScannerToken::OrganicElement("*".to_string())));
        assert!(tokens.contains(&ScannerToken::AromaticElement("a".to_string())));
        assert!(tokens.contains(&ScannerToken::OrganicElement("A".to_string())));
        assert!(tokens.contains(&ScannerToken::AtomPrimitive('H')));
        for expected in [
            ScannerToken::Semi,
            ScannerToken::Not,
            ScannerToken::And,
            ScannerToken::Or,
            ScannerToken::Colon,
            ScannerToken::Plus,
            ScannerToken::Minus,
            ScannerToken::Hash,
            ScannerToken::Underscore,
        ] {
            assert!(tokens.contains(&expected));
        }
    }

    #[test]
    fn smarts_lexer_bonds() {
        let tokens = scanner_tokens("C-=~$\\\\/-><-N", ScannerStart::Molecule);
        for symbol in ['=', '~', '$', '\\', '/'] {
            assert!(tokens.contains(&ScannerToken::BondSpec(symbol)));
        }
        assert!(tokens.contains(&ScannerToken::Minus));
        assert!(tokens.contains(&ScannerToken::DativeRight));
        assert!(tokens.contains(&ScannerToken::DativeLeft));
    }

    #[test]
    fn smarts_lexer_stereo() {
        let tokens = scanner_tokens("[@@TH@ AL^0^1^2^3^4^5]", ScannerStart::Atom);
        assert!(tokens.contains(&ScannerToken::At));
        assert!(tokens.contains(&ScannerToken::ChiralClass("TH".to_string())));
        assert!(tokens.contains(&ScannerToken::ChiralClass("AL".to_string())));
        for value in 0..=5 {
            assert!(tokens.contains(&ScannerToken::Hybridization(value)));
        }
    }

    #[test]
    fn smarts_lexer_ranges_rings() {
        let tokens = scanner_tokens("[r{2-5}]%12", ScannerStart::Molecule);
        assert!(tokens.contains(&ScannerToken::RangeOpen));
        assert!(tokens.contains(&ScannerToken::RangeClose));
        assert!(tokens.contains(&ScannerToken::Percent));
        assert_eq!(
            tokens
                .iter()
                .filter(|t| matches!(t, ScannerToken::Digit(_)))
                .count(),
            4
        );
    }

    #[test]
    fn smarts_lexer_errors() {
        assert!(matches!(
            SmartsScanner::new("C C", ScannerStart::Molecule).scan(),
            Err(SmartsParseError::UnexpectedCharacter {
                position: 1,
                character: ' ',
                ..
            })
        ));
        assert!(matches!(
            SmartsScanner::new("[C", ScannerStart::Atom).scan(),
            Err(SmartsParseError::UnclosedBracket(0))
        ));
        assert_eq!(
            scanner_tokens("C", ScannerStart::Molecule).last(),
            Some(&ScannerToken::EndOfStream)
        );
        assert_eq!(
            scanner_tokens("C\nN", ScannerStart::Molecule),
            vec![
                ScannerToken::Start(ScannerStart::Molecule),
                ScannerToken::OrganicElement("C".to_string()),
                ScannerToken::EndOfStream,
            ]
        );
    }

    #[test]
    fn smarts_grammar_meta_start() {
        assert_eq!(
            parse_smarts("").expect("empty molecule start").num_atoms(),
            0
        );
        assert!(parse_atom_entry("[C]").is_ok());
        assert_eq!(
            parse_bond_entry("=").expect("bond start"),
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double))
        );
        assert!(parse_atom_entry("").is_err());
        assert!(parse_bond_entry("").is_err());
    }

    #[test]
    fn smarts_grammar_bad_atom_def() {
        for invalid in ["[", "[]C", "[C]C", "[:]"] {
            assert!(parse_atom_entry(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn smarts_grammar_mol() {
        let mol = parse_smarts("C(=O).N").expect("branch and fragment productions");
        assert_eq!(mol.num_atoms(), 3);
        assert_eq!(mol.bond_edges, vec![(0, 1)]);

        let quadruple = parse_smarts("C$C").expect("quadruple bond production");
        assert_eq!(
            quadruple.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Quadruple))
        );

        let right = parse_smarts("C->N").expect("right dative production");
        assert_eq!(right.bond_edges, vec![(0, 1)]);
        assert_eq!(
            right.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Dative))
        );
        let left = parse_smarts("C<-N").expect("left dative production");
        assert_eq!(left.bond_edges, vec![(1, 0)]);
        assert!(parse_smarts("C)").is_err());
        assert!(matches!(
            parse_smarts("C1CC"),
            Err(SmartsParseError::UnbalancedRingClosure(1))
        ));
    }

    #[test]
    fn smarts_parser_baseline_regressions() {
        for (smarts, atomic_number, aromatic) in [
            ("C", 6, false),
            ("N", 7, false),
            ("O", 8, false),
            ("F", 9, false),
            ("Cl", 17, false),
            ("Br", 35, false),
            ("c", 6, true),
            ("n", 7, true),
        ] {
            let parsed = parse_smarts(smarts).expect("simple atom regression");
            assert_eq!(
                parsed.atom_queries,
                vec![QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number,
                    aromatic,
                })],
                "{smarts}"
            );
        }

        assert_eq!(
            parse_smarts("*").expect("wildcard regression").atom_queries,
            vec![QueryNode::Predicate(AtomQueryPredicate::Any)]
        );
        for (smarts, order) in [
            ("C-C", BondOrder::Single),
            ("C=O", BondOrder::Double),
            ("C#N", BondOrder::Triple),
            ("c:c", BondOrder::Aromatic),
        ] {
            assert_eq!(
                parse_smarts(smarts)
                    .expect("explicit bond regression")
                    .bond_queries,
                vec![QueryNode::Predicate(BondQueryPredicate::Order(order))],
                "{smarts}"
            );
        }
        assert_eq!(
            parse_smarts("C~N")
                .expect("any-bond regression")
                .bond_queries,
            vec![QueryNode::Predicate(BondQueryPredicate::Any)]
        );
    }

    #[test]
    fn smarts_grammar_atomd() {
        assert_eq!(
            parse_atom_entry("C").expect("simple_atom branch"),
            QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            })
        );
        assert_eq!(
            parse_atom_entry("[H]").expect("hydrogen_atom branch"),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1))
        );

        let expression = parse_atom_entry("[#6]").expect("atom_expr branch");
        assert_eq!(
            expression,
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
        );

        let mapped = parse_smarts("[#6:17]").expect("mapped atom_expr branch");
        assert_eq!(mapped.atom_maps, vec![Some(17)]);
        assert_eq!(
            mapped.atom_queries,
            vec![QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))]
        );
    }

    #[test]
    fn smarts_grammar_hydrogen_atom() {
        fn contains(query: &QueryNode<AtomQueryPredicate>, expected: &AtomQueryPredicate) -> bool {
            match query {
                QueryNode::Predicate(predicate) => predicate == expected,
                QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                    children.iter().any(|child| contains(child, expected))
                }
                QueryNode::Not(child) => contains(child, expected),
            }
        }

        let cases = [
            ("[H]", None, None, None),
            ("[H:1]", None, None, Some(1)),
            ("[2H]", Some(2), None, None),
            ("[2H:2]", Some(2), None, Some(2)),
            ("[H+]", None, Some(1), None),
            ("[H+:3]", None, Some(1), Some(3)),
            ("[2H-]", Some(2), Some(-1), None),
            ("[2H-:4]", Some(2), Some(-1), Some(4)),
        ];
        for (smarts, isotope, charge, atom_map) in cases {
            let parsed = parse_smarts(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}"));
            let query = &parsed.atom_queries[0];
            assert!(
                contains(query, &AtomQueryPredicate::AtomicNumber(1)),
                "{smarts} must be a hydrogen atom"
            );
            if let Some(isotope) = isotope {
                assert!(contains(query, &AtomQueryPredicate::Isotope(isotope)));
            }
            if let Some(charge) = charge {
                assert!(contains(query, &AtomQueryPredicate::FormalCharge(charge)));
            }
            assert_eq!(parsed.atom_maps[0], atom_map, "{smarts}");
        }

        assert_eq!(
            parse_atom_entry("[H2]").expect("hydrogen-count expression"),
            QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(2))
        );
    }

    #[test]
    fn smarts_atom_map_property() {
        let parsed = parse_smarts("[C:7]-[O:8]").expect("mapped SMARTS parse tree");
        assert_eq!(parsed.atom_maps, vec![Some(7), Some(8)]);
        assert_eq!(
            parsed.atom_queries,
            vec![
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 8,
                    aromatic: false,
                }),
            ]
        );

        let query = to_mol("[C:7]-[O:8]").expect("mapped query molecule");
        assert_eq!(query.atoms()[0].atom_map(), Some(7));
        assert_eq!(query.atoms()[1].atom_map(), Some(8));

        let target = Molecule::from_smiles("CO").expect("unmapped target molecule");
        assert!(
            !crate::search::substruct::get_substruct_matches(&target, &query).is_empty(),
            "atom-map properties must not constrain substructure matching"
        );
    }

    #[test]
    fn smarts_grammar_atom_expr() {
        fn atom_type(atomic_number: u8) -> QueryNode<AtomQueryPredicate> {
            QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number,
                aromatic: false,
            })
        }

        let hydrogen_count = QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(1));
        assert_eq!(
            parse_atom_entry("[C,N;H1]").expect("semicolon is lower precedence than comma"),
            QueryNode::And(vec![
                QueryNode::Or(vec![atom_type(6), atom_type(7)]),
                hydrogen_count.clone(),
            ])
        );
        assert_eq!(
            parse_atom_entry("[C;N,H1]").expect("comma is higher precedence than semicolon"),
            QueryNode::And(vec![
                atom_type(6),
                QueryNode::Or(vec![atom_type(7), hydrogen_count.clone()]),
            ])
        );
        assert_eq!(
            parse_atom_entry("[C&N,H1]").expect("ampersand is higher precedence than comma"),
            QueryNode::Or(vec![
                QueryNode::And(vec![atom_type(6), atom_type(7)]),
                hydrogen_count,
            ])
        );
        assert_eq!(
            parse_atom_entry("[!C]").expect("point-query negation"),
            QueryNode::not(atom_type(6))
        );

        for invalid in ["[]", "[!]", "[C,]", "[C&]", "[C;]", "[,C]", "[&C]"] {
            assert!(parse_atom_entry(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn smarts_grammar_point_query() {
        let carbon = QueryNode::Predicate(AtomQueryPredicate::AtomType {
            atomic_number: 6,
            aromatic: false,
        });
        assert_eq!(parse_atom_entry("[C]").expect("atom_query branch"), carbon);
        assert_eq!(
            parse_atom_entry("[!C]").expect("NOT point_query branch"),
            QueryNode::not(carbon.clone())
        );
        assert_eq!(
            parse_atom_entry("[!!C]").expect("nested NOT point_query branch"),
            carbon
        );

        let recursive = parse_atom_entry("[$(CO)]").expect("recursive_query branch");
        match recursive {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => {
                assert_eq!(query.source_smarts(), Some("$(CO)"));
                assert_eq!(query.query_mol().map(Molecule::num_atoms), Some(2));
            }
            other => panic!("expected recursive point query, got {other:?}"),
        }
        assert!(parse_atom_entry("[!]").is_err());
    }

    #[test]
    fn smarts_grammar_recursive_query() {
        fn recursive(smarts: &str) -> super::super::query::RecursiveStructureQuery {
            match parse_atom_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")) {
                QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(query)) => query,
                other => panic!("expected recursive query for {smarts}, got {other:?}"),
            }
        }

        let plain = recursive("[$(C1CC1)]");
        assert_eq!(plain.source_smarts(), Some("$(C1CC1)"));
        assert_eq!(plain.serial_number(), 0);
        assert_eq!(plain.query_mol().map(Molecule::num_atoms), Some(3));

        let numbered = recursive("[$([N])_12]");
        assert_eq!(numbered.source_smarts(), Some("$([N])"));
        assert_eq!(numbered.serial_number(), 12);
        assert_eq!(numbered.query_mol().map(Molecule::num_atoms), Some(1));

        for invalid in ["[$()]", "[$(C)_]", "[$(C)_0]", "[$(C1CC)]"] {
            assert!(parse_atom_entry(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn smarts_grammar_atom_query() {
        fn contains(query: &QueryNode<AtomQueryPredicate>, expected: &AtomQueryPredicate) -> bool {
            match query {
                QueryNode::Predicate(predicate) => predicate == expected,
                QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                    children.iter().any(|child| contains(child, expected))
                }
                QueryNode::Not(child) => contains(child, expected),
            }
        }

        let isotope_carbon = parse_atom_entry("[13C]").expect("number simple_atom");
        assert!(contains(&isotope_carbon, &AtomQueryPredicate::Isotope(13)));
        assert!(contains(
            &isotope_carbon,
            &AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }
        ));
        let isotope_iron = parse_atom_entry("[57Fe]").expect("number ATOM_TOKEN");
        assert!(contains(&isotope_iron, &AtomQueryPredicate::Isotope(57)));
        assert!(contains(
            &isotope_iron,
            &AtomQueryPredicate::AtomicNumber(26)
        ));
        let isotope_atomic_number = parse_atom_entry("[13#6]").expect("number HASH number");
        assert!(contains(
            &isotope_atomic_number,
            &AtomQueryPredicate::Isotope(13)
        ));
        assert!(contains(
            &isotope_atomic_number,
            &AtomQueryPredicate::AtomicNumber(6)
        ));
        assert_eq!(
            parse_atom_entry("[2H,13C]").expect("number H_TOKEN in atom_expr"),
            QueryNode::Or(vec![
                QueryNode::And(vec![
                    QueryNode::Predicate(AtomQueryPredicate::Isotope(2)),
                    QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(1)),
                ]),
                QueryNode::And(vec![
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: false,
                    }),
                    QueryNode::Predicate(AtomQueryPredicate::Isotope(13)),
                ]),
            ])
        );

        let expected = [
            ("[D]", AtomQueryPredicate::ExplicitDegree(1)),
            ("[D0]", AtomQueryPredicate::ExplicitDegree(0)),
            ("[X]", AtomQueryPredicate::TotalDegree(1)),
            ("[X0]", AtomQueryPredicate::TotalDegree(0)),
            ("[v]", AtomQueryPredicate::TotalValence(1)),
            ("[v0]", AtomQueryPredicate::TotalValence(0)),
            ("[d]", AtomQueryPredicate::NonHydrogenDegree(1)),
            ("[d2]", AtomQueryPredicate::NonHydrogenDegree(2)),
            ("[z]", AtomQueryPredicate::HasHeteroatomNeighbors),
            ("[z2]", AtomQueryPredicate::NumHeteroatomNeighbors(2)),
            ("[Z]", AtomQueryPredicate::HasAliphaticHeteroatomNeighbors),
            (
                "[Z2]",
                AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(2),
            ),
            ("[r0]", AtomQueryPredicate::SmallestRingSize(0)),
            ("[k5]", AtomQueryPredicate::InRingOfSize(5)),
            ("[x0]", AtomQueryPredicate::RingBondCount(0)),
            ("[h0]", AtomQueryPredicate::ImplicitHydrogenCount(0)),
        ];
        for (smarts, predicate) in expected {
            assert_eq!(
                parse_atom_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")),
                QueryNode::Predicate(predicate),
                "{smarts}"
            );
        }
        assert_eq!(
            parse_atom_entry("[k]").expect("bare ring-size token"),
            QueryNode::Predicate(AtomQueryPredicate::InRing)
        );

        let square_planar = parse_atom_entry("[@SP3]").expect("chiral class and permutation");
        assert!(contains(
            &square_planar,
            &AtomQueryPredicate::ChiralTagMatch(crate::ChiralTag::SquarePlanar)
        ));
        assert!(contains(
            &square_planar,
            &AtomQueryPredicate::ChiralPermutationMatch(3)
        ));
        assert!(parse_atom_entry("[@SP0]").is_err());
        assert!(parse_atom_entry("[^6]").is_err());
    }

    #[test]
    fn smarts_grammar_possible_range_query() {
        for smarts in [
            "[D{-2}]", "[d{1-}]", "[X{2-4}]", "[v{2-4}]", "[R{1-}]", "[z{-2}]", "[Z{1-3}]",
            "[r{3-6}]", "[x{-2}]", "[h{1-}]", "[+{-2}]", "[-{1-2}]", "[k{3-6}]",
        ] {
            assert!(
                matches!(
                    parse_atom_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")),
                    QueryNode::Predicate(AtomQueryPredicate::Range(_))
                ),
                "{smarts} must compile to the canonical range leaf"
            );
        }

        let mut degree_builder = Molecule::builder();
        let center = degree_builder.add_atom(AtomSpec::new(Element::C));
        let terminals = [
            degree_builder.add_atom(AtomSpec::new(Element::C)),
            degree_builder.add_atom(AtomSpec::new(Element::C)),
            degree_builder.add_atom(AtomSpec::new(Element::C)),
        ];
        for terminal in terminals {
            degree_builder
                .add_bond(BondSpec::new(center, terminal, BondOrder::Single))
                .expect("degree fixture bond");
        }
        let degree_mol = degree_builder.build().expect("degree fixture");
        let center_atom = &degree_mol.atoms()[center.index()];
        let terminal_atom = &degree_mol.atoms()[terminals[0].index()];
        for (smarts, center_matches, terminal_matches) in [
            ("[D{-2}]", false, true),
            ("[D{2-}]", true, false),
            ("[D{2-3}]", true, false),
        ] {
            let query = parse_atom_entry(smarts).expect("degree range query");
            assert_eq!(
                crate::search::query::atom_matches_query(center_atom, &query, &degree_mol),
                center_matches,
                "{smarts} center"
            );
            assert_eq!(
                crate::search::query::atom_matches_query(terminal_atom, &query, &degree_mol),
                terminal_matches,
                "{smarts} terminal"
            );
        }

        let mut ring_builder = Molecule::builder();
        let triangle = [
            ring_builder.add_atom(AtomSpec::new(Element::C)),
            ring_builder.add_atom(AtomSpec::new(Element::C)),
            ring_builder.add_atom(AtomSpec::new(Element::C)),
        ];
        let square = [
            ring_builder.add_atom(AtomSpec::new(Element::C)),
            ring_builder.add_atom(AtomSpec::new(Element::C)),
            ring_builder.add_atom(AtomSpec::new(Element::C)),
            ring_builder.add_atom(AtomSpec::new(Element::C)),
        ];
        for ring in [&triangle[..], &square[..]] {
            for index in 0..ring.len() {
                ring_builder
                    .add_bond(BondSpec::new(
                        ring[index],
                        ring[(index + 1) % ring.len()],
                        BondOrder::Single,
                    ))
                    .expect("ring fixture bond");
            }
        }
        let ring_mol = ring_builder.build().expect("ring fixture");
        for (smarts, triangle_matches, square_matches) in [
            ("[k{-3}]", true, false),
            ("[k{4-}]", false, true),
            ("[k{3-4}]", true, true),
        ] {
            let query = parse_atom_entry(smarts).expect("ring-size range query");
            assert_eq!(
                crate::search::query::atom_matches_query(
                    &ring_mol.atoms()[triangle[0].index()],
                    &query,
                    &ring_mol,
                ),
                triangle_matches,
                "{smarts} triangle"
            );
            assert_eq!(
                crate::search::query::atom_matches_query(
                    &ring_mol.atoms()[square[0].index()],
                    &query,
                    &ring_mol,
                ),
                square_matches,
                "{smarts} square"
            );
        }

        for invalid in ["[D{}]", "[D{-}]", "[D{2}]", "[D{2--3}]", "[k{3-4]"] {
            assert!(parse_atom_entry(invalid).is_err(), "{invalid}");
        }
    }

    #[test]
    fn smarts_grammar_simple_atom() {
        for (smarts, expected) in [
            (
                "C",
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
            ),
            (
                "[c]",
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: true,
                }),
            ),
            (
                "[si]",
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 14,
                    aromatic: true,
                }),
            ),
            (
                "[Cl]",
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 17,
                    aromatic: false,
                }),
            ),
        ] {
            assert_eq!(
                parse_atom_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")),
                expected,
                "{smarts}"
            );
        }
        assert_eq!(
            parse_atom_entry("*").expect("wildcard simple atom"),
            QueryNode::Predicate(AtomQueryPredicate::Any)
        );
        assert_eq!(
            parse_atom_entry("[a]").expect("aromatic wildcard simple atom"),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true))
        );
        assert_eq!(
            parse_atom_entry("[A]").expect("aliphatic wildcard simple atom"),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(false))
        );
    }

    #[test]
    fn smarts_grammar_bond_expr() {
        let single = QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single));
        let double = QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double));
        let aromatic = QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic));
        let in_ring = QueryNode::Predicate(BondQueryPredicate::IsInRing(true));

        assert_eq!(
            parse_bond_entry("-&@,=;!:").expect("bond-expression precedence"),
            QueryNode::And(vec![
                QueryNode::Or(vec![
                    QueryNode::And(vec![single.clone(), in_ring.clone()]),
                    double.clone(),
                ]),
                QueryNode::not(aromatic.clone()),
            ])
        );
        assert_eq!(
            parse_bond_entry("-@").expect("adjacent bondd primitives"),
            QueryNode::And(vec![single.clone(), in_ring])
        );
        assert_eq!(
            parse_bond_entry("!!:").expect("right-associative bond negation"),
            aromatic
        );
        assert_eq!(
            parse_bond_entry("~,-").expect("RDKit null-query OR algebra"),
            QueryNode::Predicate(BondQueryPredicate::Any)
        );
        assert_eq!(
            parse_smarts("C-,=O")
                .expect("molecule bond_expr")
                .bond_queries,
            vec![QueryNode::Or(vec![single, double])]
        );

        for invalid in ["", "&-", ",-", ";-", "-&", "-,", "-;"] {
            assert!(parse_bond_entry(invalid).is_err(), "{invalid:?}");
        }
    }

    #[test]
    fn smarts_grammar_bond_query() {
        let single = QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single));
        let in_ring = QueryNode::Predicate(BondQueryPredicate::IsInRing(true));
        let aromatic = QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic));

        assert_eq!(parse_bond_entry("-").expect("one bondd"), single.clone());
        assert_eq!(
            parse_bond_entry("-@!:").expect("left-associated adjacent bondd primitives"),
            QueryNode::And(vec![
                QueryNode::And(vec![single, in_ring]),
                QueryNode::not(aromatic),
            ])
        );
        assert_eq!(
            parse_bond_entry("~-").expect("RDKit null-query implicit AND algebra"),
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single))
        );
    }

    #[test]
    fn smarts_grammar_bondd() {
        for (smarts, expected) in [
            (
                "-",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Single)),
            ),
            (
                "=",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double)),
            ),
            (
                "#",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Triple)),
            ),
            (
                ":",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic)),
            ),
            (
                "$",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Quadruple)),
            ),
            (
                "@",
                QueryNode::Predicate(BondQueryPredicate::IsInRing(true)),
            ),
            ("~", QueryNode::Predicate(BondQueryPredicate::Any)),
            (
                "->",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::DativeRight)),
            ),
            (
                "<-",
                QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::DativeLeft)),
            ),
        ] {
            assert_eq!(
                parse_bond_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")),
                expected,
                "{smarts}"
            );
        }

        for directional in ["/", "\\"] {
            assert_eq!(
                parse_bond_entry(directional).expect("directional BOND_TOKEN"),
                unspecified_smarts_bond_query(),
                "{directional}"
            );
        }
        assert_eq!(
            parse_bond_entry("!@").expect("negated bondd"),
            QueryNode::not(QueryNode::Predicate(BondQueryPredicate::IsInRing(true)))
        );
        assert!(parse_bond_entry("!").is_err());
    }

    #[test]
    fn smarts_grammar_charge_spec() {
        for (smarts, charge) in [
            ("[+]", 1),
            ("[++]", 2),
            ("[+2]", 2),
            ("[+0]", 0),
            ("[-]", -1),
            ("[--]", -2),
            ("[-2]", -2),
            ("[-0]", 0),
            ("[+127]", 127),
            ("[-128]", -128),
        ] {
            assert_eq!(
                parse_atom_entry(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}")),
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(charge)),
                "{smarts}"
            );
        }
        assert_eq!(
            parse_atom_entry("[+++]").expect("third plus belongs to the next atom_query"),
            QueryNode::And(vec![
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(2)),
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(1)),
            ])
        );
        assert!(parse_atom_entry("[+128]").is_err());
        assert!(parse_atom_entry("[-129]").is_err());
    }

    #[test]
    fn smarts_grammar_ring_number() {
        assert_eq!(
            tokenize("C%12CC%12").expect("two-digit percent ring"),
            vec![
                (Token::OrganicElement("C".to_string()), 0),
                (Token::RingClosurePercent(12), 1),
                (Token::OrganicElement("C".to_string()), 4),
                (Token::OrganicElement("C".to_string()), 5),
                (Token::RingClosurePercent(12), 6),
                (Token::EndOfStream, 9),
            ]
        );
        assert_eq!(
            parse_smarts("C%(123)CC%(123)")
                .expect("parenthesized percent ring")
                .ring_closures,
            vec![(123, 0), (123, 2)]
        );
        assert!(parse_smarts("C%0C").is_err());
        assert!(parse_smarts("C%123C").is_err());
        assert!(parse_smarts("C%(123456)C").is_err());
    }

    #[test]
    fn smarts_grammar_number() {
        assert_eq!(
            parse_atom_entry("[#0]").expect("zero number"),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(0))
        );
        assert_eq!(
            parse_atom_entry("[#17]").expect("multi-digit number"),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(17))
        );
        assert!(parse_atom_entry("[#2147483648]").is_err());
    }

    #[test]
    fn smarts_grammar_nonzero_number() {
        assert!(parse_atom_entry("[13C]").is_ok());
        assert!(parse_atom_entry("[0C]").is_ok());
        assert!(parse_atom_entry("[C1]").is_ok());
    }

    #[test]
    fn smarts_grammar_digit() {
        assert!(parse_atom_entry("[C0]").is_ok());
        assert!(parse_atom_entry("[C9]").is_ok());
        assert!(parse_atom_entry("[C10]").is_ok());
    }

    #[test]
    fn smarts_grammar_branch_open_token() {
        let tokens = tokenize("C(C)").expect("branch tokens");
        let mut parser = SmartsParser::new(&tokens, "C(C)");
        parser.advance();
        assert_eq!(parser.parse_branch_open_token().expect("branch open"), 1);
        assert!(matches!(parser.peek(), (Token::OrganicElement(symbol), 2) if symbol == "C"));

        let parsed = parse_smarts("C(C)(N)O").expect("sibling branches");
        assert_eq!(parsed.num_atoms(), 4);
        assert_eq!(parsed.bond_edges, vec![(0, 1), (0, 2), (0, 3)]);
    }

    #[test]
    fn smarts_parse_generic_parse_helper() {
        assert!(matches!(
            generic_parse_helper("CC", ScannerStart::Molecule),
            Ok(tokens) if tokens.iter().any(|(token, _)| matches!(token, Token::EndOfStream))
        ));
        assert!(generic_parse_helper("=", ScannerStart::Bond).is_ok());
        assert!(generic_parse_helper("[C", ScannerStart::Atom).is_err());
    }

    #[test]
    fn smarts_parse_smarts_parse_helper() {
        let parsed = smarts_parse_helper("C=O").expect("SMARTS helper dispatch");
        assert_eq!(parsed.num_atoms(), 2);
        assert_eq!(parsed.bond_queries.len(), 1);
        assert!(smarts_parse_helper("C(").is_err());
    }

    #[test]
    fn smarts_parse_smarts_bond_parse() {
        assert_eq!(
            smarts_bond_parse("=").expect("bond entry dispatch"),
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double))
        );
        assert!(smarts_bond_parse("").is_err());
        assert!(smarts_bond_parse("C").is_err());
    }

    #[test]
    fn smarts_parse_smarts_atom_parse() {
        assert_eq!(
            smarts_atom_parse("C").expect("atom entry dispatch"),
            QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            })
        );
        assert_eq!(
            smarts_atom_parse("[#7]").expect("bracket atom entry dispatch"),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(7))
        );
        assert!(smarts_atom_parse("").is_err());
        assert!(smarts_atom_parse("=").is_err());
        assert!(smarts_atom_parse("CC").is_err());
    }

    #[test]
    fn smarts_parse_to_atom() {
        assert_eq!(to_atom("").expect("empty atom SMARTS"), None);
        assert_eq!(
            to_atom("C").expect("organic atom SMARTS"),
            Some(QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }))
        );
        assert_eq!(
            to_atom("[#7]").expect("bracket atom SMARTS"),
            Some(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(7)))
        );
        assert!(to_atom("CC").is_err());
        assert!(to_atom("=").is_err());
    }

    #[test]
    fn smarts_parse__atom_from_smarts() {
        assert_eq!(atom_from_smarts("").expect("empty atom SMARTS"), None);
        assert_eq!(
            atom_from_smarts("[13C;H2;+1]").expect("query atom"),
            to_atom("[13C;H2;+1]").expect("canonical atom parser")
        );
        assert!(atom_from_smarts("CC").is_err());
    }

    #[test]
    fn smarts_parse_to_bond() {
        assert_eq!(to_bond("").expect("empty bond SMARTS"), None);
        assert_eq!(
            to_bond("=").expect("double bond SMARTS"),
            Some(QueryNode::Predicate(BondQueryPredicate::Order(
                BondOrder::Double
            )))
        );
        assert_eq!(
            to_bond("~").expect("any bond SMARTS"),
            Some(QueryNode::Predicate(BondQueryPredicate::Any))
        );
        assert!(to_bond("C").is_err());
    }

    #[test]
    fn smarts_parse__bond_from_smarts() {
        assert_eq!(bond_from_smarts("").expect("empty bond SMARTS"), None);
        assert_eq!(
            bond_from_smarts("-,=").expect("composite bond query"),
            to_bond("-,=").expect("canonical bond parser")
        );
        assert!(bond_from_smarts("C").is_err());
    }

    #[test]
    fn smarts_is_query_h() {
        fn classify(smarts: &str, atom_index: usize) -> QueryHydrogenType {
            let molecule = to_mol(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}"));
            let degree = molecule
                .topology_block()
                .adjacency
                .neighbors_of(atom_index)
                .len();
            is_query_hydrogen(&molecule.atoms()[atom_index], degree)
        }

        assert_eq!(classify("[#1]", 0), QueryHydrogenType::QueryHydrogen);
        assert_eq!(classify("[H]", 0), QueryHydrogenType::QueryHydrogen);
        assert_eq!(classify("[!#1]", 0), QueryHydrogenType::NotAHydrogen);
        assert_eq!(
            classify("[#1,#6,#7]", 0),
            QueryHydrogenType::UnmergableQueryHydrogen
        );
        assert_eq!(classify("[#6;#1]", 0), QueryHydrogenType::QueryHydrogen);
        assert_eq!(
            classify("[#6]-[#1;H0](-[#6])-[#6]", 1),
            QueryHydrogenType::NotAHydrogen
        );
    }

    #[test]
    fn smarts_needs_hs() {
        assert!(needs_hs(&Molecule::from_smiles("CC").expect("ethane")).expect("ethane Hs"));
        assert!(
            !needs_hs(&Molecule::from_smiles("[O][O]").expect("oxygen radical"))
                .expect("oxygen radical Hs")
        );
        assert!(!needs_hs(&Molecule::from_smiles("FF").expect("fluorine")).expect("fluorine Hs"));

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(AtomSpec::new(Element::C));
        for _ in 0..4 {
            let hydrogen = builder.add_atom(AtomSpec::new(Element::H));
            builder
                .add_bond(BondSpec::new(carbon, hydrogen, BondOrder::Single))
                .expect("explicit C-H bond");
        }
        let explicit_methane = builder.build().expect("explicit methane");
        assert!(!needs_hs(&explicit_methane).expect("explicit methane Hs"));
    }

    #[test]
    fn smarts_merge_query_hs_in_place() {
        fn has_negated_h_count(query: &QueryNode<AtomQueryPredicate>, count: u8) -> bool {
            match query {
                QueryNode::Not(child) => matches!(
                    child.as_ref(),
                    QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(value)) if *value == count
                ),
                QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                    children
                        .iter()
                        .any(|child| has_negated_h_count(child, count))
                }
                QueryNode::Predicate(_) => false,
            }
        }

        let mut one_h = to_mol("C[H]").expect("one explicit query H");
        merge_query_hs_in_place(&mut one_h, false, false).expect("merge one H");
        assert_eq!(one_h.num_atoms(), 1);
        assert!(has_negated_h_count(
            one_h.atoms()[0].query().expect("expanded carbon query"),
            0
        ));

        let mut two_h = to_mol("C([H])[H]").expect("two explicit query Hs");
        merge_query_hs_in_place(&mut two_h, false, false).expect("merge two Hs");
        assert_eq!(two_h.num_atoms(), 1);
        let carbon_query = two_h.atoms()[0].query().expect("expanded carbon query");
        assert!(has_negated_h_count(carbon_query, 0));
        assert!(has_negated_h_count(carbon_query, 1));

        let mut hydrogen_molecule = to_mol("[H][H]").expect("hydrogen molecule");
        merge_query_hs_in_place(&mut hydrogen_molecule, false, false)
            .expect("preserve unconnected-to-heavy Hs");
        assert_eq!(hydrogen_molecule.num_atoms(), 2);

        let mut negated = to_mol("[!#1]-[#1]").expect("negated heavy query");
        merge_query_hs_in_place(&mut negated, false, false).expect("merge into negated query");
        assert_eq!(negated.num_atoms(), 1);

        let mut unmergeable = to_mol("[#6]-[#1,#6]").expect("OR query H");
        merge_query_hs_in_place(&mut unmergeable, false, false).expect("preserve OR query H");
        assert_eq!(unmergeable.num_atoms(), 2);

        let mut mapped = to_mol("[#6]-[#1:7]").expect("mapped query H");
        merge_query_hs_in_place(&mut mapped, true, false).expect("preserve mapped H");
        assert_eq!(mapped.num_atoms(), 2);
        merge_query_hs_in_place(&mut mapped, false, false).expect("merge mapped H");
        assert_eq!(mapped.num_atoms(), 1);

        let mut isotopic = to_mol("[#6]-[2#1]").expect("isotopic query H");
        merge_query_hs_in_place(&mut isotopic, false, false).expect("preserve isotopic H");
        assert_eq!(isotopic.num_atoms(), 2);
        merge_query_hs_in_place(&mut isotopic, false, true).expect("merge isotopic H");
        assert_eq!(isotopic.num_atoms(), 1);
    }

    #[test]
    fn smarts_merge_query_hs_copy() {
        let source = to_mol("C([H])[H]").expect("source query");
        let source_snapshot = source.clone();
        let merged = merge_query_hs(&source, false, false).expect("copied merge");
        assert_eq!(source, source_snapshot);
        assert_eq!(source.num_atoms(), 3);
        assert_eq!(merged.num_atoms(), 1);
        assert!(merged.atoms()[0].query().is_some());
    }

    #[test]
    fn smarts_has_query_hs() {
        assert_eq!(
            has_query_hs(&to_mol("CCCC").expect("no query Hs")),
            (false, false)
        );
        assert_eq!(
            has_query_hs(&to_mol("[#1]").expect("query H")),
            (true, false)
        );
        assert_eq!(
            has_query_hs(&to_mol("[#1,N]").expect("unmergeable query H")),
            (true, true)
        );
        assert_eq!(
            has_query_hs(&to_mol("[$(C-[H])]").expect("recursive query H")),
            (true, false)
        );
        assert_eq!(
            has_query_hs(&to_mol("[$([C,#1])]").expect("recursive OR query H")),
            (true, true)
        );
        for smarts in ["[#1,#6,#7]", "[1;#7,#1,#6]", "[1&#7,#1,#6]"] {
            assert_eq!(
                has_query_hs(&to_mol(smarts).unwrap_or_else(|error| panic!("{smarts}: {error}"))),
                (true, true),
                "{smarts}"
            );
        }
    }

    #[test]
    fn smarts_parse__mol_from_smarts() {
        let molecule = mol_from_smarts(
            "[C:7]=O |$carbon;oxygen$| carbonyl query",
            &SmartsParseParams::default(),
        )
        .expect("SMARTS molecule entry");
        assert_eq!(molecule.num_atoms(), 2);
        assert_eq!(molecule.num_bonds(), 1);
        assert_eq!(molecule.atoms()[0].atom_map(), Some(7));
        assert_eq!(molecule.atoms()[0].prop("atomLabel"), Some("carbon"));
        assert_eq!(molecule.atoms()[1].prop("atomLabel"), Some("oxygen"));
        assert_eq!(molecule.properties().name(), Some("carbonyl query"));
        assert!(molecule.atoms().iter().all(|atom| atom.query().is_some()));
        assert!(molecule.bonds().iter().all(|bond| bond.query().is_some()));

        let merge_params = SmartsParseParams {
            merge_hs: true,
            ..SmartsParseParams::default()
        };
        let merged = mol_from_smarts("[H]C", &merge_params).expect("merge query hydrogen");
        assert_eq!(merged.num_atoms(), 1);
        assert_eq!(merged.num_bonds(), 0);

        let recursive =
            mol_from_smarts("[$(C-[H])]", &merge_params).expect("merge recursive query hydrogen");
        assert_eq!(has_query_hs(&recursive), (false, false));

        let mapped = mol_from_smarts("[#6]-[#1:7]", &merge_params)
            .expect("MolFromSmarts merges mapped H with default options");
        assert_eq!(mapped.num_atoms(), 1);

        let isotopic = mol_from_smarts("[#6]-[2#1]", &merge_params)
            .expect("MolFromSmarts preserves isotopic H with default options");
        assert_eq!(isotopic.num_atoms(), 2);
    }

    #[test]
    fn smarts_parse_smarts_parse() {
        assert_eq!(
            smarts_parse_entry("")
                .expect("empty molecule entry")
                .num_atoms(),
            0
        );

        let connected = smarts_parse_entry("C=O").expect("connected molecule entry");
        assert_eq!(connected.num_atoms(), 2);
        assert_eq!(connected.bond_edges, vec![(0, 1)]);

        let fragments = smarts_parse_entry("C.N").expect("fragment molecule entry");
        assert_eq!(fragments.num_atoms(), 2);
        assert!(fragments.bond_edges.is_empty());

        assert!(smarts_parse_entry("C(").is_err());
    }

    #[test]
    fn smarts_parse_to_mol() {
        let empty = to_mol("").expect("empty SMARTS molecule");
        assert_eq!(empty.num_atoms(), 0);
        assert_eq!(empty.num_bonds(), 0);

        let connected = to_mol("[C:7]=O").expect("connected query molecule");
        assert_eq!(connected.num_atoms(), 2);
        assert_eq!(connected.num_bonds(), 1);
        assert_eq!(connected.atoms()[0].atom_map(), Some(7));
        assert!(connected.atoms().iter().all(|atom| atom.query().is_some()));
        assert!(connected.bonds().iter().all(|bond| bond.query().is_some()));

        let ring = to_mol("C1CC1").expect("ring query molecule");
        assert_eq!(ring.num_atoms(), 3);
        assert_eq!(ring.num_bonds(), 3);

        let fragments = to_mol("C.N").expect("disconnected query molecule");
        assert_eq!(fragments.num_atoms(), 2);
        assert_eq!(fragments.num_bonds(), 0);

        assert!(to_mol("C(").is_err());
        assert!(to_mol("C1CC").is_err());
    }

    #[test]
    fn test_tokenize_unclosed_bracket() {
        let result = tokenize("[NH");
        assert!(result.is_err());
        assert!(matches!(
            result.unwrap_err(),
            SmartsParseError::UnclosedBracket(_)
        ));
    }

    #[test]
    fn test_parse_simple_atom() {
        assert_eq!(
            parse_simple_atom("C"),
            Some(QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            }))
        );
        assert_eq!(
            parse_simple_atom("*"),
            Some(QueryNode::Predicate(AtomQueryPredicate::Any))
        );
        assert_eq!(parse_simple_atom("Xe"), None);
    }

    #[test]
    fn test_bracket_atom_element() {
        // [C] — carbon in brackets
        let mol = parse_smarts("[C]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
        assert_eq!(mol.bond_queries.len(), 0);
    }

    #[test]
    fn test_parse_simple_smarts() {
        // CC
        let mol = parse_smarts("CC").unwrap();
        assert_eq!(mol.atom_queries.len(), 2);
        assert_eq!(mol.bond_queries.len(), 1);
        assert_eq!(
            mol.atom_queries[0],
            QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            })
        );
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))
        );
    }

    #[test]
    fn test_parse_bonded_smarts() {
        // C=O
        let mol = parse_smarts("C=O").unwrap();
        assert_eq!(mol.atom_queries.len(), 2);
        assert_eq!(mol.bond_queries.len(), 1);
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Double))
        );
    }

    #[test]
    fn test_bracket_with_charge() {
        let mol = parse_smarts("[N+]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
        assert_eq!(
            mol.atom_queries[0],
            QueryNode::And(vec![
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 7,
                    aromatic: false,
                }),
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(1)),
            ])
        );
    }

    #[test]
    fn test_bracket_with_negative_charge_defaults_to_minus_one() {
        let mol = parse_smarts("[O-]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
        assert_eq!(
            mol.atom_queries[0],
            QueryNode::And(vec![
                QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 8,
                    aromatic: false,
                }),
                QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-1)),
            ])
        );
    }

    #[test]
    fn test_bracket_with_chirality() {
        let mol = parse_smarts("[C@@H]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
    }

    #[test]
    fn test_parse_ring_closure() {
        let mol = parse_smarts("C1CC1").unwrap();
        assert_eq!(mol.atom_queries.len(), 3);
        assert_eq!(mol.ring_closures.len(), 2);
    }

    #[test]
    fn test_parse_explicit_bond_ring_closure_like_rdkit() {
        let mol = parse_smarts("*1~*~*~*~1").unwrap();
        assert_eq!(mol.atom_queries.len(), 4);
        assert_eq!(mol.bond_queries.len(), 4);
        assert_eq!(mol.bond_edges, vec![(0, 1), (1, 2), (2, 3), (0, 3)]);
        assert_eq!(mol.ring_closures, vec![(1, 0), (1, 3)]);
        assert_eq!(mol.ring_closure_bonds.len(), 2);
        assert_eq!(
            mol.ring_closure_bonds[1],
            (1, 3, QueryNode::Predicate(BondQueryPredicate::Any))
        );
    }

    #[test]
    fn test_parse_branch() {
        let mol = parse_smarts("C(C)C").unwrap();
        assert_eq!(mol.atom_queries.len(), 3);
    }

    #[test]
    fn test_parse_branch_continues_from_branch_point_like_rdkit() {
        let mol = parse_smarts("[#6]~[!#6!#1](~[#6])(~[#6])~*").unwrap();
        assert_eq!(mol.atom_queries.len(), 5);
        assert_eq!(mol.bond_queries.len(), 4);
        assert_eq!(mol.bond_edges, vec![(0, 1), (1, 2), (1, 3), (1, 4)]);
    }

    #[test]
    fn test_bracket_atomic_number_primitive() {
        let mol = parse_smarts("[#6]").unwrap();
        assert_eq!(mol.atom_queries.len(), 1);
    }

    #[test]
    fn test_label_recursive_patterns_noop() {
        // No $() in input
        assert_eq!(label_recursive_patterns("CCO"), "CCO");
    }

    #[test]
    fn test_label_recursive_patterns_simple() {
        let result = label_recursive_patterns("[$([N])]");
        // Should append label after closing paren
        assert!(result.contains("_100") || result == "[$([N])]");
    }

    #[test]
    fn smarts_parse_label_recursive_patterns() {
        assert_eq!(label_recursive_patterns("CCO"), "CCO");
        assert_eq!(label_recursive_patterns("C(C)N"), "C(C)N");
        assert_eq!(
            label_recursive_patterns("[$(C),$(C),$(N)]"),
            "[$(C)_100,$(C)_100,$(N)_101]"
        );
        assert_eq!(
            label_recursive_patterns("[$(C)_777,$(N)]"),
            "[$(C)_777,$(N)_100]"
        );
        assert_eq!(label_recursive_patterns("[$($(C))]"), "[$($(C)_100)_101]");
        assert_eq!(label_recursive_patterns("C)($(N))"), "C)($(N))");
        assert_eq!(label_recursive_patterns("[é,$(N)]"), "[é,$(N)_100]");
    }

    #[test]
    fn test_element_symbol_lookup() {
        assert_eq!(element_symbol_to_atomic_number("C"), Some(6));
        assert_eq!(element_symbol_to_atomic_number("O"), Some(8));
        assert_eq!(element_symbol_to_atomic_number("Cl"), Some(17));
        assert_eq!(element_symbol_to_atomic_number("Br"), Some(35));
        assert_eq!(element_symbol_to_atomic_number("Xx"), None);
    }

    #[test]
    fn test_parse_aromatic_smarts() {
        // c1ccccc1 — benzene SMARTS
        let mol = parse_smarts("c1ccccc1").unwrap();
        assert_eq!(mol.atom_queries.len(), 6);
        // Every atom is aromatic carbon
        for aq in &mol.atom_queries {
            assert_eq!(
                aq,
                &QueryNode::Predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: true,
                })
            );
        }
    }

    #[test]
    fn test_smarts_molecule_num_atoms() {
        let mol = parse_smarts("CCO").unwrap();
        assert_eq!(mol.num_atoms(), 3);
        assert!(mol.atom_query(0).is_some());
        assert!(mol.bond_query(0).is_some());
        assert!(mol.bond_query(1).is_some());
    }

    #[test]
    fn ring_connectivity_distinguishes_bare_x_from_explicit_x0() {
        let has_ring_bond = parse_smarts("[Cx]").expect("bare x SMARTS");
        assert!(atom_query_contains(
            &has_ring_bond.atom_queries[0],
            &AtomQueryPredicate::HasRingBond
        ));

        let no_ring_bonds = parse_smarts("[Cx0]").expect("x0 SMARTS");
        assert!(atom_query_contains(
            &no_ring_bonds.atom_queries[0],
            &AtomQueryPredicate::RingBondCount(0)
        ));
    }

    #[test]
    fn default_feature_smarts_parse_into_expected_query_shapes() {
        let cases = [
            (
                "Donor",
                "[$([N;!H0;v3,v4&+1]),$([O,S;H1;+0]),n&H1&+0]",
                "or",
            ),
            (
                "Acceptor",
                "[$([O,S;H1;v2;!$(*-*=[O,N,P,S])]),$([O,S;H0;v2]),$([O,S;-]),$([N;v3;!$(N-*=[O,N,P,S])]),n&H0&+0,$([o,s;+0;!$([o,s]:n);!$([o,s]:c:n)])]",
                "or",
            ),
            ("Aromatic", "[a]", "aromatic"),
            ("Halogen", "[F,Cl,Br,I]", "or"),
            (
                "Basic",
                "[#7;+,$([N;H2&+0][$([C,a]);!$([C,a](=O))]),$([N;H1&+0]([$([C,a]);!$([C,a](=O))])[$([C,a]);!$([C,a](=O))]),$([N;H0&+0]([C;!$(C(=O))])([C;!$(C(=O))])[C;!$(C(=O))])]",
                "and",
            ),
            ("Acidic", "[$([C,S](=[O,S,P])-[O;H1,-1])]", "recursive"),
        ];

        for (name, pattern, expected_shape) in cases {
            let parsed =
                parse_smarts(pattern).unwrap_or_else(|_| panic!("{name} SMARTS should parse"));
            assert_eq!(parsed.atom_queries.len(), 1, "{name} atom query count");
            let atom_query = &parsed.atom_queries[0];
            match expected_shape {
                "or" => assert!(matches!(atom_query, QueryNode::Or(_)), "{name} root"),
                "and" => assert!(matches!(atom_query, QueryNode::And(_)), "{name} root"),
                "aromatic" => assert!(
                    matches!(
                        atom_query,
                        QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true))
                    ),
                    "{name} root"
                ),
                "recursive" => assert!(
                    matches!(
                        atom_query,
                        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(_))
                    ),
                    "{name} root"
                ),
                _ => panic!("unexpected expected shape for {name}"),
            }
        }
    }

    fn atom_query_contains(
        query: &QueryNode<AtomQueryPredicate>,
        predicate: &AtomQueryPredicate,
    ) -> bool {
        match query {
            QueryNode::Predicate(candidate) => candidate == predicate,
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children
                    .iter()
                    .any(|child| atom_query_contains(child, predicate))
            }
            QueryNode::Not(child) => atom_query_contains(child, predicate),
        }
    }

    fn atom_query_contains_recursive_smarts(query: &QueryNode<AtomQueryPredicate>) -> bool {
        match query {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(_)) => true,
            QueryNode::Predicate(_) => false,
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children.iter().any(atom_query_contains_recursive_smarts)
            }
            QueryNode::Not(child) => atom_query_contains_recursive_smarts(child),
        }
    }

    fn bond_query_contains(
        query: &QueryNode<BondQueryPredicate>,
        predicate: &BondQueryPredicate,
    ) -> bool {
        match query {
            QueryNode::Predicate(candidate) => candidate == predicate,
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children
                    .iter()
                    .any(|child| bond_query_contains(child, predicate))
            }
            QueryNode::Not(child) => bond_query_contains(child, predicate),
        }
    }

    #[test]
    fn maccs_patterns_parse_targeted_source_smarts_categories() {
        // RDKit source: MACCS.cpp::Patterns constructs these SMARTS strings
        // with RDKit::SmartsToMol before matching them in GenerateFP().
        let recursive = parse_smarts(
            "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
        )
        .expect("MACCS bit 90 recursive SMARTS should parse");
        assert_eq!(recursive.num_atoms(), 1);
        assert!(atom_query_contains_recursive_smarts(
            &recursive.atom_queries[0]
        ));

        let ring_atom = parse_smarts("[R]").expect("MACCS bit 165 ring atom should parse");
        assert!(atom_query_contains(
            &ring_atom.atom_queries[0],
            &AtomQueryPredicate::InRing
        ));

        let ring_bond = parse_smarts("*@*(@*)@*").expect("MACCS bit 105 ring bonds should parse");
        assert!(
            ring_bond
                .bond_queries
                .iter()
                .any(|query| bond_query_contains(query, &BondQueryPredicate::IsInRing(true))),
            "MACCS bit 105 should preserve @ ring-bond queries"
        );

        let non_ring_bond =
            parse_smarts("*!@[#8]!@*").expect("MACCS bit 126 non-ring bonds should parse");
        assert!(
            non_ring_bond.bond_queries.iter().any(|query| matches!(
                query,
                QueryNode::Not(child)
                    if **child == QueryNode::Predicate(BondQueryPredicate::IsInRing(true))
            )),
            "MACCS bit 126 should preserve !@ non-ring-bond queries"
        );

        let wildcard = parse_smarts("*~[CH2]~[#7]").expect("MACCS bit 100 wildcard should parse");
        assert_eq!(
            wildcard.atom_queries[0],
            QueryNode::Predicate(AtomQueryPredicate::Any)
        );

        let negation = parse_smarts("[!#6!#1!H0]").expect("MACCS bit 131 negation should parse");
        assert!(
            matches!(&negation.atom_queries[0], QueryNode::And(children) if children.iter().any(|child| matches!(child, QueryNode::Not(_)))),
            "MACCS bit 131 should preserve atom-query negation"
        );

        let alternatives =
            parse_smarts("[F,Cl,Br,I]").expect("MACCS bit 31 OR alternatives should parse");

        fn or_leaf_count(query: &QueryNode<AtomQueryPredicate>) -> Option<usize> {
            match query {
                QueryNode::Or(children) => children.iter().try_fold(0, |count, child| {
                    or_leaf_count(child).map(|child_count| count + child_count)
                }),
                QueryNode::Predicate(_) => Some(1),
                QueryNode::And(_) | QueryNode::Xor(_) | QueryNode::Not(_) => None,
            }
        }
        assert!(
            matches!(&alternatives.atom_queries[0], QueryNode::Or(_))
                && or_leaf_count(&alternatives.atom_queries[0]) == Some(4),
            "MACCS bit 31 should parse four halogen alternatives"
        );

        let hydrogen_count =
            parse_smarts("[C;H3,H4]").expect("MACCS bit 149 hydrogen counts should parse");
        assert!(atom_query_contains(
            &hydrogen_count.atom_queries[0],
            &AtomQueryPredicate::HydrogenCount(3)
        ));
        assert!(atom_query_contains(
            &hydrogen_count.atom_queries[0],
            &AtomQueryPredicate::HydrogenCount(4)
        ));

        let branch_ring_closure = parse_smarts("[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]")
            .expect("MACCS bit 116 branch/ring-closure SMARTS should parse");
        assert_eq!(branch_ring_closure.num_atoms(), 1);
        assert!(atom_query_contains_recursive_smarts(
            &branch_ring_closure.atom_queries[0]
        ));

        let explicit_ring_closure =
            parse_smarts("*1~*~*~*~1").expect("MACCS bit 11 ring closure should parse");
        assert_eq!(
            explicit_ring_closure.bond_edges,
            vec![(0, 1), (1, 2), (2, 3), (0, 3)]
        );
        assert_eq!(explicit_ring_closure.ring_closures, vec![(1, 0), (1, 3)]);
        assert_eq!(explicit_ring_closure.ring_closure_bonds.len(), 2);
    }

    #[test]
    fn maccs_patterns_parse_source_smarts() {
        let cases = [
            (8, "[!#6!#1]1~*~*~*~1"),
            (11, "*1~*~*~*~1"),
            (13, "[#8]~[#7](~[#6])~[#6]"),
            (14, "[#16]-[#16]"),
            (15, "[#8]~[#6](~[#8])~[#8]"),
            (16, "[!#6!#1]1~*~*~1"),
            (17, "[#6]#[#6]"),
            (19, "*1~*~*~*~*~*~*~1"),
            (20, "[#14]"),
            (21, "[#6]=[#6](~[!#6!#1])~[!#6!#1]"),
            (22, "*1~*~*~1"),
            (23, "[#7]~[#6](~[#8])~[#8]"),
            (24, "[#7]-[#8]"),
            (25, "[#7]~[#6](~[#7])~[#7]"),
            (26, "[#6]=@[#6](@*)@*"),
            (28, "[!#6!#1]~[CH2]~[!#6!#1]"),
            (30, "[#6]~[!#6!#1](~[#6])(~[#6])~*"),
            (31, "[!#6!#1]~[F,Cl,Br,I]"),
            (32, "[#6]~[#16]~[#7]"),
            (33, "[#7]~[#16]"),
            (34, "[CH2]=*"),
            (36, "[#16R]"),
            (37, "[#7]~[#6](~[#8])~[#7]"),
            (38, "[#7]~[#6](~[#6])~[#7]"),
            (39, "[#8]~[#16](~[#8])~[#8]"),
            (40, "[#16]-[#8]"),
            (41, "[#6]#[#7]"),
            (43, "[!#6!#1!H0]~*~[!#6!#1!H0]"),
            (44, "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]"),
            (45, "[#6]=[#6]~[#7]"),
            (47, "[#16]~*~[#7]"),
            (48, "[#8]~[!#6!#1](~[#8])~[#8]"),
            (49, "[!+0]"),
            (50, "[#6]=[#6](~[#6])~[#6]"),
            (51, "[#6]~[#16]~[#8]"),
            (52, "[#7]~[#7]"),
            (53, "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]"),
            (54, "[!#6!#1!H0]~*~*~[!#6!#1!H0]"),
            (55, "[#8]~[#16]~[#8]"),
            (56, "[#8]~[#7](~[#8])~[#6]"),
            (57, "[#8R]"),
            (58, "[!#6!#1]~[#16]~[!#6!#1]"),
            (59, "[#16]!:*:*"),
            (60, "[#16]=[#8]"),
            (61, "*~[#16](~*)~*"),
            (62, "*@*!@*@*"),
            (63, "[#7]=[#8]"),
            (64, "*@*!@[#16]"),
            (65, "c:n"),
            (66, "[#6]~[#6](~[#6])(~[#6])~*"),
            (67, "[!#6!#1]~[#16]"),
            (68, "[!#6!#1!H0]~[!#6!#1!H0]"),
            (69, "[!#6!#1]~[!#6!#1!H0]"),
            (70, "[!#6!#1]~[#7]~[!#6!#1]"),
            (71, "[#7]~[#8]"),
            (72, "[#8]~*~*~[#8]"),
            (73, "[#16]=*"),
            (74, "[CH3]~*~[CH3]"),
            (75, "*!@[#7]@*"),
            (76, "[#6]=[#6](~*)~*"),
            (77, "[#7]~*~[#7]"),
            (78, "[#6]=[#7]"),
            (79, "[#7]~*~*~[#7]"),
            (80, "[#7]~*~*~*~[#7]"),
            (81, "[#16]~*(~*)~*"),
            (82, "*~[CH2]~[!#6!#1!H0]"),
            (83, "[!#6!#1]1~*~*~*~*~1"),
            (84, "[NH2]"),
            (85, "[#6]~[#7](~[#6])~[#6]"),
            (86, "[C;H2,H3][!#6!#1][C;H2,H3]"),
            (87, "[F,Cl,Br,I]!@*@*"),
            (89, "[#8]~*~*~*~[#8]"),
            (
                90,
                "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]",
            ),
            (
                91,
                "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (92, "[#8]~[#6](~[#7])~[#6]"),
            (93, "[!#6!#1]~[CH3]"),
            (94, "[!#6!#1]~[#7]"),
            (95, "[#7]~*~*~[#8]"),
            (96, "*1~*~*~*~*~1"),
            (97, "[#7]~*~*~*~[#8]"),
            (98, "[!#6!#1]1~*~*~*~*~*~1"),
            (99, "[#6]=[#6]"),
            (100, "*~[CH2]~[#7]"),
            (
                101,
                "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]",
            ),
            (102, "[!#6!#1]~[#8]"),
            (104, "[!#6!#1!H0]~*~[CH2]~*"),
            (105, "*@*(@*)@*"),
            (106, "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]"),
            (107, "[F,Cl,Br,I]~*(~*)~*"),
            (108, "[CH3]~*~*~*~[CH2]~*"),
            (109, "*~[CH2]~[#8]"),
            (110, "[#7]~[#6]~[#8]"),
            (111, "[#7]~*~[CH2]~*"),
            (112, "*~*(~*)(~*)~*"),
            (113, "[#8]!:*:*"),
            (114, "[CH3]~[CH2]~*"),
            (115, "[CH3]~*~[CH2]~*"),
            (116, "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]"),
            (117, "[#7]~*~[#8]"),
            (118, "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]"),
            (119, "[#7]=*"),
            (120, "[!#6R]"),
            (121, "[#7R]"),
            (122, "*~[#7](~*)~*"),
            (123, "[#8]~[#6]~[#8]"),
            (124, "[!#6!#1]~[!#6!#1]"),
            (126, "*!@[#8]!@*"),
            (127, "*@*!@[#8]"),
            (
                128,
                "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]",
            ),
            (
                129,
                "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]",
            ),
            (131, "[!#6!#1!H0]"),
            (132, "[#8]~*~[CH2]~*"),
            (133, "*@*!@[#7]"),
            (135, "[#7]!:*:*"),
            (136, "[#8]=*"),
            (137, "[!C!cR]"),
            (138, "[!#6!#1]~[CH2]~*"),
            (139, "[O!H0]"),
            (140, "[#8]"),
            (141, "[CH3]"),
            (142, "[#7]"),
            (144, "*!:*:*!:*"),
            (145, "*1~*~*~*~*~*~1"),
            (147, "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]"),
            (148, "*~[!#6!#1](~*)~*"),
            (149, "[C;H3,H4]"),
            (150, "*!@*@*!@*"),
            (151, "[#7!H0]"),
            (152, "[#8]~[#6](~[#6])~[#6]"),
            (154, "[#6]=[#8]"),
            (155, "*!@[CH2]!@*"),
            (156, "[#7]~*(~*)~*"),
            (157, "[#6]-[#8]"),
            (158, "[#6]-[#7]"),
            (162, "a"),
            (165, "[R]"),
        ];

        assert_eq!(cases.len(), 136);
        for (bit, smarts) in cases {
            parse_smarts(smarts)
                .unwrap_or_else(|error| panic!("MACCS bit {bit} SMARTS failed: {error}"));
        }
    }

    #[test]
    fn test_smarts_parse_params_default() {
        let params = SmartsParseParams::default();
        assert!(params.allow_cxsmiles);
        assert!(!params.merge_hs);
    }

    #[test]
    fn smarts_debug_parse_is_explicitly_unsupported() {
        let params = SmartsParseParams {
            debug_parse: true,
            ..SmartsParseParams::default()
        };
        assert!(matches!(
            mol_from_smarts("C", &params),
            Err(SmartsParseError::UnsupportedFeature(
                "Bison debug_parse diagnostic output"
            ))
        ));
    }

    #[test]
    fn smarts_parse_preprocess_smiles() {
        let default_params = SmartsParseParams::default();
        assert_eq!(
            preprocess_smarts("CC |$;;;_R1$| molecule name", &default_params),
            PreprocessedSmarts {
                smarts: "CC".to_string(),
                name: String::new(),
                cx_part: "|$;;;_R1$| molecule name".to_string(),
            }
        );

        let mut name_params = SmartsParseParams {
            allow_cxsmiles: false,
            ..SmartsParseParams::default()
        };
        assert_eq!(
            preprocess_smarts("C\t  query name  ", &name_params),
            PreprocessedSmarts {
                smarts: "C".to_string(),
                name: "query name".to_string(),
                cx_part: String::new(),
            }
        );

        name_params.parse_name = false;
        assert_eq!(
            preprocess_smarts("C query name", &name_params).smarts,
            "C query name"
        );
        assert_eq!(
            preprocess_smarts(" C", &default_params).smarts,
            " C",
            "a separator at byte zero leaves the input intact"
        );

        let mut replacement_params = SmartsParseParams {
            allow_cxsmiles: false,
            parse_name: false,
            ..SmartsParseParams::default()
        };
        replacement_params
            .replacements
            .insert("{A}".to_string(), "N".to_string());
        replacement_params
            .replacements
            .insert("{B}".to_string(), "{A}".to_string());
        assert_eq!(
            preprocess_smarts("C{B}", &replacement_params).smarts,
            "CN",
            "replacements repeat until no key remains"
        );
    }

    #[test]
    fn smarts_parse_handle_c_x_part_and_name() {
        let params = SmartsParseParams::default();
        let mut molecule = to_mol("C=O").expect("query molecule");
        let atom_queries = molecule
            .atoms()
            .iter()
            .map(|atom| atom.query().cloned())
            .collect::<Vec<_>>();
        let bond_queries = molecule
            .bonds()
            .iter()
            .map(|bond| bond.query().cloned())
            .collect::<Vec<_>>();

        let mut unchanged_name = "existing".to_string();
        handle_c_x_part_and_name(&mut molecule, &params, "", &mut unchanged_name)
            .expect("empty CX part");
        assert_eq!(unchanged_name, "existing");

        let mut plain_name = String::new();
        handle_c_x_part_and_name(&mut molecule, &params, "  query name  ", &mut plain_name)
            .expect("plain trailing name");
        assert_eq!(plain_name, "query name");

        let strict_no_name = SmartsParseParams {
            parse_name: false,
            ..SmartsParseParams::default()
        };
        let mut ignored_name = String::new();
        assert_eq!(
            handle_c_x_part_and_name(
                &mut molecule,
                &strict_no_name,
                "not-a-cx-extension",
                &mut ignored_name,
            ),
            Err(SmartsParseError::CxSmiles(
                "CXSMILES extension does not start with | and parseName=false".to_string()
            ))
        );

        let mut cx_name = String::new();
        handle_c_x_part_and_name(
            &mut molecule,
            &params,
            "|$carbonyl-carbon;oxygen$,atomProp:1.note.carbonyl| named query",
            &mut cx_name,
        )
        .expect("shared CX parser");
        assert_eq!(
            molecule.atoms()[0].prop("atomLabel"),
            Some("carbonyl-carbon")
        );
        assert_eq!(molecule.atoms()[1].prop("atomLabel"), Some("oxygen"));
        assert_eq!(molecule.atoms()[1].prop("note"), Some("carbonyl"));
        assert_eq!(
            molecule.properties().prop("_CXSMILES_Data"),
            Some("|$carbonyl-carbon;oxygen$,atomProp:1.note.carbonyl|")
        );
        assert_eq!(cx_name, "named query");
        assert_eq!(
            molecule
                .atoms()
                .iter()
                .map(|atom| atom.query().cloned())
                .collect::<Vec<_>>(),
            atom_queries
        );
        assert_eq!(
            molecule
                .bonds()
                .iter()
                .map(|bond| bond.query().cloned())
                .collect::<Vec<_>>(),
            bond_queries
        );

        let non_strict = SmartsParseParams {
            strict_cxsmiles: false,
            ..SmartsParseParams::default()
        };
        let mut failed_name = String::new();
        handle_c_x_part_and_name(
            &mut molecule,
            &non_strict,
            "|unterminated",
            &mut failed_name,
        )
        .expect("non-strict CX failure");
        assert_eq!(molecule.properties().prop("_CXSMILES_Data"), Some(""));
        assert!(failed_name.is_empty());
    }

    #[test]
    fn test_varied_bonds() {
        let mol = parse_smarts("C#N").unwrap();
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Triple))
        );

        let mol = parse_smarts("C:N").unwrap();
        assert_eq!(
            mol.bond_queries[0],
            QueryNode::Predicate(BondQueryPredicate::Order(BondOrder::Aromatic))
        );
    }

    #[test]
    fn test_empty_smarts_molecule() {
        let result = parse_smarts("").expect("meta_start EOS accepts an empty molecule");
        assert_eq!(result.num_atoms(), 0);
    }
}
