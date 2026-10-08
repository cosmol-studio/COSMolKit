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
//! This module is the sole canonical SMARTS parser/compiler owner.  Parsing
//! produces the first-class [`QueryGraph`] model consumed by matching and
//! serialization; concrete molecules remain query-free at the public API
//! boundary.

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::rc::Rc;

use cosmolkit_cx::{
    CxCoordinateBondKind, CxParseProgress, CxProgressPhase, CxRecord, CxSGroupHierarchy,
};

use crate::query_behavior::SmartsParseError;
use crate::query_behavior::{
    AtomRangeDataFunction, CompositeQueryType, make_atom_null_query,
    make_atom_possible_range_query, make_atom_possible_ring_range_query,
    make_bond_is_in_ring_query, make_bond_null_query, make_bond_order_equals_query,
    query_bond_expand_query,
};
use crate::{QueryAtom, QueryBond, QueryGraph};
use cosmolkit_model::{
    AtomId, AtomMapping, AtomQueryPredicate, Bond, BondId, BondMapping, BondQueryPredicate,
    BondSpec, PropertyText, QueryAtomIdentity, QueryNode, StereoGroup, SubstanceGroup,
    SubstanceGroupId, TopologyMapping, query_substance_groups, replace_query_substance_groups,
};
use cosmolkit_types::{BondDirection, BondOrder, ChiralTag, Element, Hybridization};

#[cfg(test)]
use cosmolkit_model::AtomSpec;

// ---------------------------------------------------------------------------
// QueryGraphBuilder - private parser construction state
// ---------------------------------------------------------------------------

/// Parser-owned construction state for a [`QueryGraph`].
///
/// This parser keeps source row indexes private and lowers them once into a
/// canonical `QueryGraph`.
/// COSMolKit keeps parser indexes and ring bookkeeping private, then lowers
/// them once into the independent `QueryGraph` value. This is not a second
/// public graph model and never becomes a query-bearing `Molecule`.
#[derive(Debug, Clone, Default)]
struct QueryGraphBuilder {
    /// One complete source query carrier per atom in the pattern.
    atoms: Vec<QueryAtom>,
    /// Query trees for source-ordered bonds, including reconciled ring closures.
    bond_queries: Vec<QueryNode<BondQueryPredicate>>,
    /// Source QueryBond types aligned with `bond_queries`.
    bond_orders: Vec<BondOrder>,
    /// Directional state aligned with `bond_queries`.
    bond_directions: Vec<BondDirection>,
    /// Source grammar encounter indexes, independent of final graph rows.
    bond_source_indexes: Vec<u32>,
    /// Source unspecified-order properties aligned with emitted bond rows.
    bond_unspecified_orders: Vec<bool>,
    /// Query bond endpoints in SMARTS atom-index space.
    bond_edges: Vec<(usize, usize)>,
    /// Undirected endpoint index for source `getBondBetweenAtoms` checks.
    bond_pairs: BTreeSet<(usize, usize)>,
}

fn ordered_atom_pair(first: usize, second: usize) -> (usize, usize) {
    // Constant-time canonicalization keeps the endpoint index undirected.
    if first <= second {
        (first, second)
    } else {
        (second, first)
    }
}

impl QueryGraphBuilder {
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    #[must_use]
    pub fn atom_query(&self, idx: usize) -> Option<&QueryNode<AtomQueryPredicate>> {
        self.atoms.get(idx).map(QueryAtom::predicate)
    }

    #[must_use]
    pub fn bond_query(&self, idx: usize) -> Option<&QueryNode<BondQueryPredicate>> {
        self.bond_queries.get(idx)
    }

    fn push_atom(&mut self, atom: ParsedSmartsAtom) -> usize {
        let index = self.atoms.len();
        let mut carrier = atom.carrier.with_id(AtomId::new(index));
        carrier.set_atom_map(atom.atom_map);
        self.atoms.push(carrier);
        index
    }

    fn push_bond(
        &mut self,
        begin: usize,
        end: usize,
        bond: ParsedSmartsBond,
        direction: BondDirection,
        source_index: u32,
    ) {
        // Local complexity: endpoint lookup is O(log E) in the BTreeSet and
        // keeps one O(E) duplicate-edge index alongside the source rows.
        self.bond_pairs.insert(ordered_atom_pair(begin, end));
        self.bond_queries.push(bond.query);
        self.bond_orders.push(bond.carrier_order);
        self.bond_directions.push(direction);
        self.bond_source_indexes.push(source_index);
        self.bond_unspecified_orders.push(bond.unspecified_order);
        self.bond_edges.push((begin, end));
    }

    fn has_bond_between(&self, first: usize, second: usize) -> bool {
        // One canonicalization and one O(log E) lookup; no row rescan.
        self.bond_pairs.contains(&ordered_atom_pair(first, second))
    }

    fn implicit_bond_order(&self, begin: usize, end: usize) -> BondOrder {
        // BEGIN RDKIT CPP FUNCTION GetUnspecifiedBondType
        // RDKit✔️✔️:   if (atom1->getIsAromatic() && atom2->getIsAromatic()) {
        // RDKit✔️✔️:     res = Bond::AROMATIC;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     res = Bond::SINGLE;
        // RDKit✔️✔️:   }
        // END RDKIT CPP FUNCTION GetUnspecifiedBondType
        // SMARTS implicit bonds are created through getUnspecifiedQueryBond;
        // its source order is determined by the two source atom flags. Atom
        // carrier specs are indexed in the same parser row space, so this is
        // one lookup per endpoint and does not inspect the bond predicate.
        let both_aromatic = [begin, end]
            .into_iter()
            .all(|index| self.atoms.get(index).is_some_and(QueryAtom::is_aromatic));
        if both_aromatic {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        }
    }

    /// Lower parser indexes and query state into the canonical query value.
    ///
    /// Parser row storage is deliberately consumed here: query carriers,
    /// bond directions, and endpoints are installed on their final
    /// `QueryAtom`/`QueryBond` values. Ring-closure records remain in
    /// `SmartsParser` and are consumed when their ordinary bond is emitted;
    /// the builder does not keep a cloned mirror of them.
    fn finish(self) -> Result<QueryGraph, SmartsParseError> {
        // BEGIN COMPLETE PINNED SF250 toMol
        // RDKit✔️❌: std::unique_ptr<RWMol> toMol(const std::string &inp,
        // RDKit✔️❌:                              int func(const std::string &,
        // RDKit✔️❌:                                       std::vector<RDKit::RWMol *> &),
        // RDKit✔️❌:                              const std::string &origInp) {
        // RDKit✔️❌:   // empty strings produce empty molecules:
        // RDKit✔️❌:   if (inp.empty()) {
        // RDKit✔️❌:     return std::make_unique<RWMol>();
        // RDKit✔️❌:   }
        // RDKit✔️❌:   std::unique_ptr<RWMol> res;
        // RDKit✔️❌:   std::vector<RDKit::RWMol *> molVect;
        // RDKit✔️❌:   try {
        // RDKit✔️❌:     func(inp, molVect);
        // RDKit✔️❌:     if (!molVect.empty()) {
        // RDKit✔️❌:       res.reset(molVect[0]);
        // RDKit✔️❌:       SmilesParseOps::CloseMolRings(res.get(), false);
        // RDKit✔️❌:       SmilesParseOps::CheckChiralitySpecifications(res.get(), true);
        // RDKit✔️❌:       SmilesParseOps::SetUnspecifiedBondTypes(res.get());
        // RDKit✔️❌:       SmilesParseOps::AdjustAtomChiralityFlags(res.get());
        // RDKit✔️❌:       // No sense leaving this bookmark intact:
        // RDKit✔️❌:       if (res->hasAtomBookmark(ci_RIGHTMOST_ATOM)) {
        // RDKit✔️❌:         res->clearAtomBookmark(ci_RIGHTMOST_ATOM);
        // RDKit✔️❌:       }
        // RDKit✔️❌:       molVect[0] = nullptr;  // NOTE: to avoid leaks on failures, this should
        // RDKit✔️❌:                              // occur last in this if.
        // RDKit✔️❌:     }
        // RDKit✔️❌:   } catch (SmilesParseException &e) {
        // RDKit✔️❌:     std::string nm = "SMILES";
        // RDKit✔️❌:     if (func == smarts_parse) {
        // RDKit✔️❌:       nm = "SMARTS";
        // RDKit✔️❌:     }
        // RDKit✔️❌:     BOOST_LOG(rdErrorLog) << nm << " Parse Error: " << e.what()
        // RDKit✔️❌:                           << " for input: '" << origInp << "'" << std::endl;
        // RDKit✔️❌:
        // RDKit✔️❌:     // reset res so that we return a nullptr. We don't want to reset(),
        // RDKit✔️❌:     // because that would delete the mol and leak any unmatched
        // RDKit✔️❌:     // ring closure bonds. These will be cleaned up in the loop below.
        // RDKit✔️❌:     res.release();
        // RDKit✔️❌:   }
        // RDKit✔️❌:   for (auto *molPtr : molVect) {
        // RDKit✔️❌:     if (molPtr) {
        // RDKit✔️❌:       // Clean-up the bond bookmarks when not calling CloseMolRings
        // RDKit✔️❌:       SmilesParseOps::CleanupAfterParseError(molPtr);
        // RDKit✔️❌:       delete molPtr;
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌:
        // RDKit✔️❌:   return res;
        // RDKit✔️❌: }
        // END COMPLETE PINNED SF250 toMol
        // BEGIN COMPLETE PINNED reached SetUnspecifiedBondTypes
        // RDKit✔️❌: void SetUnspecifiedBondTypes(RWMol *mol) {
        // RDKit✔️❌:   PRECONDITION(mol, "no molecule");
        // RDKit✔️❌:   for (auto bond : mol->bonds()) {
        // RDKit✔️❌:     if (bond->hasProp(RDKit::common_properties::_unspecifiedOrder)) {
        // RDKit✔️❌:       bond->setBondType(GetUnspecifiedBondType(mol, bond->getBeginAtom(),
        // RDKit✔️❌:                                                bond->getEndAtom()));
        // RDKit✔️❌:       if (bond->getBondType() == Bond::AROMATIC) {
        // RDKit✔️❌:         bond->setIsAromatic(true);
        // RDKit✔️❌:       } else {
        // RDKit✔️❌:         bond->setIsAromatic(false);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // RDKit✔️❌: }
        // END COMPLETE PINNED reached SetUnspecifiedBondTypes
        // RDKit✔️✔️: Bond::BondType GetUnspecifiedBondType(const RWMol *mol, const Atom *atom1,
        // RDKit✔️✔️:                                       const Atom *atom2) {
        // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
        // RDKit✔️✔️:   PRECONDITION(atom1, "no atom1");
        // RDKit✔️✔️:   PRECONDITION(atom2, "no atom2");
        // RDKit✔️✔️:   Bond::BondType res;
        // RDKit✔️✔️:   if (atom1->getIsAromatic() && atom2->getIsAromatic()) {
        // RDKit✔️✔️:     res = Bond::AROMATIC;
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     res = Bond::SINGLE;
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return res;
        // RDKit✔️✔️: }
        // The sole scanner/parser has already consumed EOS and closed every
        // pending ring in source bookmark order before this graph is finished.
        // Empty input produced the canonical empty builder, without scanning.
        // Recursive reductions use finish_recursive and skip these root passes.
        // Source strict chirality check precedes unspecified bond flags and
        // shared chirality adjustment. Setter dispatch here is Bond*, so its
        // type/aromatic flags change without replacing the query predicate.
        // No live rightmost bookmark escapes the parser's owned local state;
        // Result errors drop all owned root/recursive rows rather than leaking
        // or returning a partly finalized graph. No exception is suppressed.
        // Cost: graph assembly plus one V chirality pass, one E bond pass and
        // canonical chirality ordering. Detached validation and assignment
        // buffers add O(V+E) storage beyond the source in-place finalization.
        let mut graph = self.finish_recursive()?;
        for atom in graph.atoms_mut() {
            materialize_smarts_atom_state(atom)?;
        }
        for bond_index in 0..graph.num_bonds() {
            let bond = &graph.bonds()[bond_index];
            if bond
                .bond()
                .prop(crate::query_graph_behavior::UNSPECIFIED_ORDER_PROP)
                .is_some()
            {
                let both_aromatic = graph.atoms()[bond.begin().index()].is_aromatic()
                    && graph.atoms()[bond.end().index()].is_aromatic();
                let order = if both_aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                };
                let carrier = graph.bonds_mut()[bond_index].bond_mut();
                carrier.set_order(order);
                carrier.set_aromatic(both_aromatic);
            }
        }
        crate::query_graph_behavior::finalize_query_parser_chirality(&mut graph)?;
        Ok(graph)
    }

    fn finish_recursive(self) -> Result<QueryGraph, SmartsParseError> {
        // RDKit❗✔️: recursive_query: BEGIN_RECURSE mol END_RECURSE {
        // RDKit❗✔️:   // this is a recursive SMARTS expression
        // RDKit❗✔️:   QueryAtom *qA = new QueryAtom();
        // RDKit❗✔️:   //  FIX: there's maybe a leak here
        // RDKit❗✔️:   RWMol *molP = (*molList)[$2];
        // RDKit❗✔️:   // close any rings in the molecule:
        // RDKit❗✔️:   SmilesParseOps::CloseMolRings(molP,0);
        // RDKit❗✔️:
        // RDKit❗✔️:   //molP->debugMol(std::cout);
        // RDKit❗✔️:   qA->setQuery(new RecursiveStructureQuery(molP));
        // RDKit❗✔️:   //std::cout << "qA: " << qA << " " << qA->getQuery() << std::endl;
        let Self {
            mut atoms,
            bond_queries,
            bond_orders,
            bond_directions,
            bond_source_indexes,
            bond_unspecified_orders,
            bond_edges,
            bond_pairs: _,
        } = self;

        if bond_queries.len() != bond_directions.len()
            || bond_queries.len() != bond_orders.len()
            || bond_queries.len() != bond_edges.len()
            || bond_queries.len() != bond_source_indexes.len()
            || bond_queries.len() != bond_unspecified_orders.len()
            || atoms
                .iter()
                .enumerate()
                .any(|(index, atom)| atom.id() != AtomId::new(index))
        {
            return Err(SmartsParseError::Parse(
                "SMARTS parser state has misaligned graph arrays".to_owned(),
            ));
        }

        // BEGIN RDKIT CPP BLOCK SMARTS graph row creation
        // RDKit❗❗: int atomIdx2 = mp->addAtom($3,true,true);
        // RDKit❗❗: $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit❗❗: mp->addBond($2);
        // END RDKIT CPP BLOCK SMARTS graph row creation
        // The exact grammar anchors above record parser-order insertion and
        // the source bond-row property. Each complete QueryAtom is carried
        // directly into the one validated QueryGraph; concrete Atom/AtomSpec
        // constraints are not part of this parser identity boundary.
        // Complexity review: row indexing is O(V+E); model construction and
        // validation add allocation/work beyond RDKit's incremental RWMol
        // insertion, and their relative cost is unresolved.
        let mut bonds = Vec::with_capacity(bond_queries.len());
        for (bond_index, ((endpoints, query), direction)) in bond_edges
            .into_iter()
            .zip(bond_queries)
            .zip(bond_directions)
            .enumerate()
        {
            let (begin, end) = endpoints;
            if begin >= atoms.len() || end >= atoms.len() {
                return Err(SmartsParseError::Parse(format!(
                    "SMARTS bond {bond_index} references an atom outside the graph"
                )));
            }
            let mut bond = Bond::from_spec(
                cosmolkit_model::BondId::new(bond_index),
                BondSpec::new(
                    cosmolkit_model::AtomId::new(begin),
                    cosmolkit_model::AtomId::new(end),
                    bond_orders[bond_index],
                )
                .with_direction(direction)
                .with_prop(
                    crate::query_graph_behavior::CXSMILES_BOND_IDX_PROP,
                    cosmolkit_model::PropertyValue::UInt(bond_source_indexes[bond_index]),
                )
                .expect("the internal CXSMARTS bond-index property key is non-empty"),
            );
            // RDKit✔️✔️: newB->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
            if bond_unspecified_orders[bond_index] {
                bond.set_prop(crate::query_graph_behavior::UNSPECIFIED_ORDER_PROP, 1_i32)
                    .expect("the internal unspecified-order property key is non-empty");
            }
            bonds.push(QueryBond::from_parts(bond, query));
        }

        let mut graph = QueryGraph::from_parts(
            atoms,
            bonds,
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
        Ok(graph)
    }
}

#[cfg(test)]
fn query_graph_for_test(inp: &str) -> Result<QueryGraph, String> {
    // Test-only parser helper: keep the result in the canonical QueryGraph
    // representation rather than recreating the removed query-bearing
    // Molecule projection used by the historical tests.
    if inp.is_empty() {
        return QueryGraph::from_parts(
            Vec::new(),
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .map_err(|error| error.to_string());
    }
    let graph = parse_smarts_graph(inp)
        .map_err(|error| error.to_string())?
        .finish()
        .map_err(|error| error.to_string())?;
    Ok(graph)
}

fn materialize_smarts_atom_state(atom: &mut QueryAtom) -> Result<(), SmartsParseError> {
    // RDKit✔️✔️: void CheckChiralitySpecifications(RDKit::RWMol *mol, bool strict) {
    // RDKit✔️✔️:   PRECONDITION(mol, "no molecule");
    // RDKit✔️✔️:   for (const auto atom : mol->atoms()) {
    // RDKit✔️✔️:     int permutation;
    // RDKit✔️✔️:     if (atom->getChiralTag() > RDKit::Atom::ChiralType::CHI_OTHER &&
    // RDKit✔️✔️:         permutationLimits.find(atom->getChiralTag()) !=
    // RDKit✔️✔️:             permutationLimits.end() &&
    // RDKit✔️✔️:         atom->getPropIfPresent(common_properties::_chiralPermutation,
    // RDKit✔️✔️:                                permutation)) {
    // RDKit✔️✔️:       if (!checkChiralPermutation(atom->getChiralTag(), permutation)) {
    // RDKit✔️✔️:         std::string error =
    // RDKit✔️✔️:             (boost::format("Invalid chiral specification on atom %d") %
    // RDKit✔️✔️:              atom->getIdx())
    // RDKit✔️✔️:                 .str();
    // RDKit✔️✔️:         BOOST_LOG(rdWarningLog) << error << std::endl;
    // RDKit✔️✔️:         if (strict) {
    // RDKit✔️✔️:           throw SmilesParseException(error);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // directly convert @TH1 -> @ and @TH2 -> @@
    // RDKit✔️✔️:       if (atom->getChiralTag() == RDKit::Atom::ChiralType::CHI_TETRAHEDRAL) {
    // RDKit✔️✔️:         if (permutation == 0 || permutation == 1) {
    // RDKit✔️✔️:           atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CCW);
    // RDKit✔️✔️:           atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:         } else if (permutation == 2) {
    // RDKit✔️✔️:           atom->setChiralTag(RDKit::Atom::ChiralType::CHI_TETRAHEDRAL_CW);
    // RDKit✔️✔️:           atom->clearProp(common_properties::_chiralPermutation);
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // The grammar supplies its source int32 permutation in the existing
    // typed field. Standard tags short-circuit before the scalar read. Only
    // a present mapped-class permutation is checked; TH0/TH1/TH2 are converted
    // and their permutation cleared. Constant class dispatch and scalar access.
    let tag = atom.chiral_tag();
    if !matches!(
        tag,
        ChiralTag::Tetrahedral
            | ChiralTag::Allene
            | ChiralTag::SquarePlanar
            | ChiralTag::TrigonalBipyramidal
            | ChiralTag::Octahedral
    ) {
        return Ok(());
    }
    if let Some(permutation) = atom.chiral_permutation() {
        if !crate::query_graph_behavior::check_chiral_permutation(tag, permutation as i32) {
            return Err(SmartsParseError::Parse(format!(
                "Invalid chiral specification on atom {}",
                atom.id().index()
            )));
        }
        if tag == ChiralTag::Tetrahedral {
            atom.set_chiral_tag(if permutation <= 1 {
                ChiralTag::TetrahedralCcw
            } else {
                ChiralTag::TetrahedralCw
            });
            atom.set_chiral_permutation(None);
        }
    }
    Ok(())
}

#[doc(hidden)]
pub fn compile_query_fixture(smarts: &str) -> Result<QueryGraph, String> {
    parse_smarts(smarts, &SmartsParseParams::default()).map_err(|error| error.to_string())
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
    pub replacements: BTreeMap<PropertyText, PropertyText>,
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
fn parse_smarts_graph(smarts: impl AsRef<[u8]>) -> Result<QueryGraphBuilder, SmartsParseError> {
    let smarts = smarts.as_ref();
    parse_smarts_with_params(smarts, &SmartsParseParams::default())
}

/// Parse a SMARTS string with custom parameters.
///
/// This private syntax-tree entry applies the same preprocessing and recursive
/// labeling as the public molecule compiler, but intentionally leaves molecule
/// postprocessing to `mol_from_smarts`.
fn parse_smarts_with_params(
    smarts: impl AsRef<[u8]>,
    params: &SmartsParseParams,
) -> Result<QueryGraphBuilder, SmartsParseError> {
    let smarts = smarts.as_ref();
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

/// Apply representation-independent CX records directly to the canonical
/// SMARTS query graph. This is deliberately separate from the SMILES lowerer:
/// query predicates must never be materialized through a concrete molecule.
fn apply_cx_to_query(
    graph: &mut QueryGraph,
    records: &[CxRecord],
    stereo_tracker: &mut crate::cx_lowering::CxStereoGroupTracker,
    cx_sequence_id: &mut u32,
) -> Result<(), SmartsParseError> {
    for record in records {
        match record {
            CxRecord::Coordinates(coordinates) => {
                append_query_conformer(graph, coordinates)?;
            }
            CxRecord::AtomLabels(values) => {
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value
                        && let Some(atom) = graph.atom_mut(index)
                    {
                        atom.set_prop("atomLabel", value)?;
                    }
                }
            }
            CxRecord::AtomValues(values) => {
                for (index, value) in values.iter().enumerate() {
                    if let Some(value) = value
                        && let Some(atom) = graph.atom_mut(index)
                    {
                        atom.set_prop("molFileValue", value)?;
                    }
                }
            }
            CxRecord::AtomProperties(properties) => {
                for property in properties {
                    if let Some(atom) = graph.atom_mut(property.atom) {
                        atom.set_prop(property.name.clone(), property.value.clone())
                            .map_err(SmartsParseError::from)?;
                    }
                }
            }
            CxRecord::CoordinateBonds(annotation) => {
                for reference in &annotation.bonds {
                    apply_cx_coordinate_bond_to_query(graph, *reference, annotation.kind)?;
                }
            }
            CxRecord::ZeroBonds(indices) => {
                for index in indices {
                    apply_cx_zero_bond_to_query(graph, *index)?;
                }
            }
            CxRecord::Unsaturation(indices) => {
                for item_index in 0..indices.len() {
                    crate::cx_lowering::apply_cx_query_constraint_item(graph, record, item_index)
                        .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::RingBonds(constraints) => {
                for item_index in 0..constraints.len() {
                    crate::cx_lowering::apply_cx_query_constraint_item(graph, record, item_index)
                        .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::Substitution(constraints) => {
                for item_index in 0..constraints.len() {
                    crate::cx_lowering::apply_cx_query_constraint_item(graph, record, item_index)
                        .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::EnhancedStereo(stereo) => {
                crate::cx_lowering::merge_cx_enhanced_stereo(graph, stereo_tracker, stereo)
                    .map_err(SmartsParseError::from)?;
            }
            CxRecord::WedgedBonds(wedges) => {
                for wedge in wedges {
                    crate::cx_lowering::apply_cx_wedge_bond_to_query(graph, wedge)
                        .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::DoubleBondStereo(stereo) => {
                for &index in &stereo.bonds {
                    crate::cx_lowering::apply_cx_double_bond_stereo_to_query(
                        graph,
                        index,
                        stereo.stereo,
                    )
                    .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::Radicals(radicals) => {
                for radical in radicals {
                    apply_cx_radical_to_query(graph, *radical)?;
                }
            }
            CxRecord::LinkNodes(nodes) => {
                crate::cx_lowering::apply_cx_link_nodes_to_query(graph, nodes)
                    .map_err(SmartsParseError::from)?;
            }
            CxRecord::DataSGroup(data) => {
                crate::cx_lowering::apply_cx_data_sgroup_to_query(graph, data, *cx_sequence_id)
                    .map_err(SmartsParseError::from)?;
                *cx_sequence_id = cx_sequence_id.wrapping_add(1);
            }
            CxRecord::SGroupHierarchy(hierarchies) => {
                crate::cx_lowering::apply_cx_sgroup_hierarchy_to_query(graph, hierarchies)
                    .map_err(SmartsParseError::from)?;
            }
            CxRecord::PolymerSGroup(polymer) => {
                crate::cx_lowering::apply_cx_polymer_sgroup_to_query(
                    graph,
                    polymer,
                    *cx_sequence_id,
                )
                .map_err(SmartsParseError::from)?;
                *cx_sequence_id = cx_sequence_id.wrapping_add(1);
            }
            CxRecord::VariableAttachments(attachments) => {
                for attachment in attachments {
                    crate::cx_lowering::apply_cx_variable_attachment_to_query(graph, attachment)
                        .map_err(SmartsParseError::from)?;
                }
            }
            CxRecord::Unknown(_) => {}
        }
    }
    Ok(())
}

struct PendingCxSGroupHierarchy {
    record_index: usize,
    next_item_index: usize,
    hierarchy_index: usize,
    child_index: usize,
    parent_resolved: bool,
    next_parent_index: usize,
    resolved_parent: Option<(SubstanceGroupId, u32)>,
    groups: Vec<SubstanceGroup>,
    dirty: bool,
}

fn apply_cx_sgroup_hierarchy_progress_item(
    pending: &mut PendingCxSGroupHierarchy,
    hierarchies: &[CxSGroupHierarchy],
    item_index: usize,
) -> Result<(), crate::cx_lowering::CxQueryLoweringError> {
    if pending.next_item_index != item_index {
        return Err(crate::cx_lowering::CxQueryLoweringError::InvalidGraph(
            "CX SGroup hierarchy item checkpoints are out of order".to_owned(),
        ));
    }
    let hierarchy = hierarchies.get(pending.hierarchy_index).ok_or_else(|| {
        crate::cx_lowering::CxQueryLoweringError::InvalidGraph(
            "CX SGroup hierarchy item checkpoint references a missing parent".to_owned(),
        )
    })?;
    if !pending.parent_resolved || pending.child_index >= hierarchy.children.len() {
        return Err(crate::cx_lowering::CxQueryLoweringError::InvalidGraph(
            "CX SGroup hierarchy child has no matching parent checkpoint".to_owned(),
        ));
    }
    let child_id = hierarchy.children[pending.child_index];
    let changed = crate::cx_lowering::apply_cx_sgroup_hierarchy_child(
        &mut pending.groups,
        pending.resolved_parent,
        child_id,
    )?;
    pending.dirty |= changed;
    pending.next_item_index += 1;
    pending.child_index += 1;
    Ok(())
}

fn commit_pending_cx_sgroup_hierarchy(
    graph: &mut QueryGraph,
    pending: PendingCxSGroupHierarchy,
) -> Result<(), SmartsParseError> {
    if pending.dirty {
        replace_query_substance_groups(graph, pending.groups).map_err(SmartsParseError::from)?;
    }
    Ok(())
}

#[cfg(test)]
fn apply_cx_progress_to_query(
    graph: &mut QueryGraph,
    progress: &CxParseProgress,
) -> Result<(), SmartsParseError> {
    let mut source_cursor = progress.consumed();
    apply_cx_progress_to_query_with_cursor(graph, progress, &mut source_cursor)
}

fn apply_cx_progress_to_query_with_cursor(
    graph: &mut QueryGraph,
    progress: &CxParseProgress,
    source_cursor: &mut usize,
) -> Result<(), SmartsParseError> {
    let mut pending_coordinates: Option<(usize, cosmolkit_cx::CxCoordinates)> = None;
    let mut pending_stereo: Option<(usize, usize)> = None;
    let mut pending_query_constraints: Option<(usize, usize)> = None;
    let mut pending_link_nodes: Option<(usize, usize)> = None;
    let mut pending_data_sgroup: Option<(usize, usize)> = None;
    let mut pending_polymer_sgroup: Option<(usize, usize)> = None;
    let mut pending_variable_attachments: Option<(usize, usize)> = None;
    let mut pending_wedges: Option<(usize, usize)> = None;
    let mut pending_double_bond_stereo: Option<(usize, usize)> = None;
    let mut pending_sgroup_hierarchy: Option<PendingCxSGroupHierarchy> = None;
    let mut cx_sequence_id = 0_u32;
    let mut stereo_tracker = crate::cx_lowering::CxStereoGroupTracker::new(graph);
    for checkpoint in progress.checkpoints() {
        // Source helpers commit the current item before consuming its next
        // delimiter. Preserve that iterator position if graph lowering fails.
        *source_cursor = checkpoint.cursor;
        let record = progress
            .records()
            .get(checkpoint.record_index)
            .ok_or_else(|| {
                SmartsParseError::CxSmiles(
                    "CX progress checkpoint references a missing record".to_owned(),
                )
            })?;
        if let CxRecord::DataSGroup(data) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_data_sgroup.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX data SGroup progress has an invalid begin checkpoint".to_owned(),
                        ));
                    }
                    pending_data_sgroup = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX data SGroup item checkpoint has no field index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_data_sgroup else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX data SGroup item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index
                        || *committed != item_index
                        || item_index > 5
                        || (item_index == 5 && data.coordinates.is_none())
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX data SGroup field checkpoints are out of order".to_owned(),
                        ));
                    }
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_data_sgroup.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX data SGroup completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    let expected_fields = 5 + usize::from(data.coordinates.is_some());
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != expected_fields
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX data SGroup completion does not match parsed fields".to_owned(),
                        ));
                    }
                    crate::cx_lowering::apply_cx_data_sgroup_to_query(graph, data, cx_sequence_id)
                        .map_err(SmartsParseError::from)?;
                    cx_sequence_id = cx_sequence_id.wrapping_add(1);
                }
            }
            continue;
        }
        if let CxRecord::PolymerSGroup(polymer) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_polymer_sgroup.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX polymer SGroup progress has an invalid begin checkpoint".to_owned(),
                        ));
                    }
                    pending_polymer_sgroup = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX polymer SGroup item checkpoint has no field index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_polymer_sgroup else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX polymer SGroup item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index || *committed != item_index {
                        return Err(SmartsParseError::CxSmiles(
                            "CX polymer SGroup field checkpoints are out of order".to_owned(),
                        ));
                    }
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_polymer_sgroup.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX polymer SGroup completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    let expected_items = polymer.atoms.len()
                        + polymer.head_crossings.len()
                        + polymer.tail_crossings.len()
                        + usize::from(!polymer.label.is_empty())
                        + usize::from(!polymer.connect.is_empty());
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != expected_items
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX polymer SGroup completion does not match parsed fields".to_owned(),
                        ));
                    }
                    crate::cx_lowering::apply_cx_polymer_sgroup_to_query(
                        graph,
                        polymer,
                        cx_sequence_id,
                    )
                    .map_err(SmartsParseError::from)?;
                    cx_sequence_id = cx_sequence_id.wrapping_add(1);
                }
            }
            continue;
        }
        if let CxRecord::WedgedBonds(wedges) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_wedges.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX wedge progress has an invalid begin checkpoint".to_owned(),
                        ));
                    }
                    pending_wedges = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX wedge item checkpoint has no pair index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_wedges else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX wedge item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index || *committed != item_index {
                        return Err(SmartsParseError::CxSmiles(
                            "CX wedge item checkpoints are out of order".to_owned(),
                        ));
                    }
                    let wedge = wedges.get(item_index).ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX wedge item references a missing pair".to_owned(),
                        )
                    })?;
                    crate::cx_lowering::apply_cx_wedge_bond_to_query(graph, wedge)
                        .map_err(SmartsParseError::from)?;
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_wedges.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX wedge completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != wedges.len()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX wedge completion does not match parsed pairs".to_owned(),
                        ));
                    }
                }
            }
            continue;
        }
        if let CxRecord::DoubleBondStereo(stereo) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_double_bond_stereo.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX double-bond stereo progress has an invalid begin checkpoint"
                                .to_owned(),
                        ));
                    }
                    pending_double_bond_stereo = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX double-bond stereo item checkpoint has no bond index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_double_bond_stereo else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX double-bond stereo item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index || *committed != item_index {
                        return Err(SmartsParseError::CxSmiles(
                            "CX double-bond stereo item checkpoints are out of order".to_owned(),
                        ));
                    }
                    let bond_index = *stereo.bonds.get(item_index).ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX double-bond stereo item references a missing bond".to_owned(),
                        )
                    })?;
                    crate::cx_lowering::apply_cx_double_bond_stereo_to_query(
                        graph,
                        bond_index,
                        stereo.stereo,
                    )
                    .map_err(SmartsParseError::from)?;
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_double_bond_stereo.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX double-bond stereo completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != stereo.bonds.len()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX double-bond stereo completion does not match parsed bonds"
                                .to_owned(),
                        ));
                    }
                }
            }
            continue;
        }
        if let CxRecord::VariableAttachments(attachments) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_variable_attachments.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX variable-attachment progress has an invalid begin checkpoint"
                                .to_owned(),
                        ));
                    }
                    pending_variable_attachments = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX variable-attachment item checkpoint has no field index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_variable_attachments else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX variable-attachment item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index || *committed != item_index {
                        return Err(SmartsParseError::CxSmiles(
                            "CX variable-attachment item checkpoints are out of order".to_owned(),
                        ));
                    }
                    let attachment = attachments.get(item_index / 2).ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX variable-attachment item references a missing row".to_owned(),
                        )
                    })?;
                    if item_index % 2 == 0 {
                        // Source degree validation occurs after at1idx and before its colon.
                        crate::cx_lowering::validate_cx_variable_attachment_atom_to_query(
                            graph, attachment,
                        )
                        .map_err(SmartsParseError::from)?;
                    } else {
                        // The source installs bond properties only after the full endpoint list.
                        crate::cx_lowering::apply_cx_variable_attachment_effect_to_query(
                            graph, attachment,
                        )
                        .map_err(SmartsParseError::from)?;
                    }
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_variable_attachments.take()
                    else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX variable-attachment completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != attachments.len() * 2
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX variable-attachment completion does not match parsed rows"
                                .to_owned(),
                        ));
                    }
                }
            }
            continue;
        }
        if let CxRecord::SGroupHierarchy(hierarchies) = record {
            // BEGIN COMPLETE PINNED SF197 hierarchy replay
            // RDKit✔️❌: bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
            // RDKit✔️❌:   // these look like: |SgH:1:0|
            // RDKit✔️❌:   // from CXSMILES docs:
            // RDKit✔️❌:   //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
            // RDKit✔️❌:   if (first >= last || *first != 'S' || first + 3 >= last ||
            // RDKit✔️❌:       *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
            // RDKit✔️❌:     return false;
            // RDKit✔️❌:   }
            // RDKit✔️❌:   first += 4;
            // RDKit✔️❌:   auto &sgs = getSubstanceGroups(mol);
            // RDKit✔️❌:   while (1) {
            // RDKit✔️❌:     unsigned int parentId;
            // RDKit✔️❌:     if (!read_int(first, last, parentId)) {
            // RDKit✔️❌:       return false;
            // RDKit✔️❌:     }
            // RDKit✔️❌:
            // RDKit✔️❌:     bool validParent = true;
            // RDKit✔️❌:     auto psg = find_matching_sgroup(sgs, parentId);
            // RDKit✔️❌:     if (psg == sgs.end()) {
            // RDKit✔️❌:       validParent = false;
            // RDKit✔️❌:     } else {
            // RDKit✔️❌:       psg->getPropIfPresent("index", parentId);
            // RDKit✔️❌:     }
            // RDKit✔️❌:     if (first <= last && *first == ':') {
            // RDKit✔️❌:       ++first;
            // RDKit✔️❌:       std::vector<unsigned int> children;
            // RDKit✔️❌:       if (!read_int_list(first, last, children, '.')) {
            // RDKit✔️❌:         return false;
            // RDKit✔️❌:       }
            // RDKit✔️❌:       if (validParent) {
            // RDKit✔️❌:         for (auto childId : children) {
            // RDKit✔️❌:           if (childId >= sgs.size()) {
            // RDKit✔️❌:             throw SmilesParseException(
            // RDKit✔️❌:                 "child id references non-existent SGroup");
            // RDKit✔️❌:           }
            // RDKit✔️❌:           auto csg = find_matching_sgroup(sgs, childId);
            // RDKit✔️❌:           if (csg != sgs.end()) {
            // RDKit✔️❌:             unsigned int cid;
            // RDKit✔️❌:             csg->getProp("index", cid);
            // RDKit✔️❌:             csg->setProp("PARENT", parentId);
            // RDKit✔️❌:           }
            // RDKit✔️❌:         }
            // RDKit✔️❌:       }
            // RDKit✔️❌:       if (first <= last && *first == ',') {
            // RDKit✔️❌:         ++first;
            // RDKit✔️❌:       } else {
            // RDKit✔️❌:         break;
            // RDKit✔️❌:       }
            // RDKit✔️❌:     } else {
            // RDKit✔️❌:       return false;
            // RDKit✔️❌:     }
            // RDKit✔️❌:   }
            // RDKit✔️❌:
            // RDKit✔️❌:   return true;
            // RDKit✔️❌: }
            // END COMPLETE PINNED SF197 hierarchy replay
            // Behavior: replay parent conversion before colon/list grammar, then
            // source-ordered children only after the whole list succeeds. Retain
            // prior writes before a later parent/child property failure.
            // Complexity: pending group copy/checkpoints add source-relative cost;
            // the canonical finder/converters are reused without reparsing input.
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if let Some(hierarchy_index) = checkpoint.item_index {
                        let Some(pending) = pending_sgroup_hierarchy.as_mut() else {
                            return Err(SmartsParseError::CxSmiles(
                                "CX SGroup hierarchy parent has no record begin".to_owned(),
                            ));
                        };
                        let hierarchy = hierarchies.get(hierarchy_index).ok_or_else(|| {
                            SmartsParseError::CxSmiles(
                                "CX SGroup hierarchy parent references a missing row".to_owned(),
                            )
                        })?;
                        if pending.record_index != checkpoint.record_index
                            || pending.next_parent_index != hierarchy_index
                            || (pending.parent_resolved
                                && pending.child_index
                                    != hierarchies[pending.hierarchy_index].children.len())
                        {
                            let pending = pending_sgroup_hierarchy
                                .take()
                                .expect("pending hierarchy was just checked");
                            commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                            return Err(SmartsParseError::CxSmiles(
                                "CX SGroup hierarchy parent checkpoints are out of order"
                                    .to_owned(),
                            ));
                        }
                        // Pinned parent getPropIfPresent precedes colon/list
                        // parsing, including empty children and later failure.
                        match crate::cx_lowering::resolve_cx_sgroup_hierarchy_parent(
                            &pending.groups,
                            hierarchy.parent,
                        ) {
                            Ok(resolved_parent) => {
                                pending.hierarchy_index = hierarchy_index;
                                pending.child_index = 0;
                                pending.parent_resolved = true;
                                pending.next_parent_index += 1;
                                pending.resolved_parent = resolved_parent;
                            }
                            Err(error) => {
                                let pending = pending_sgroup_hierarchy
                                    .take()
                                    .expect("pending hierarchy was just checked");
                                commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                                return Err(SmartsParseError::CxSmiles(error.to_string()));
                            }
                        }
                        continue;
                    }
                    if pending_sgroup_hierarchy.is_some() {
                        if let Some(pending) = pending_sgroup_hierarchy.take() {
                            commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                        }
                        return Err(SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy progress has an invalid begin checkpoint"
                                .to_owned(),
                        ));
                    }
                    pending_sgroup_hierarchy = Some(PendingCxSGroupHierarchy {
                        record_index: checkpoint.record_index,
                        next_item_index: 0,
                        hierarchy_index: 0,
                        child_index: 0,
                        parent_resolved: false,
                        next_parent_index: 0,
                        resolved_parent: None,
                        groups: query_substance_groups(graph).to_vec(),
                        dirty: false,
                    });
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy item checkpoint has no child index".to_owned(),
                        )
                    })?;
                    let Some(pending) = pending_sgroup_hierarchy.as_mut() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if pending.record_index != checkpoint.record_index {
                        let pending = pending_sgroup_hierarchy
                            .take()
                            .expect("pending hierarchy was just checked");
                        commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                        return Err(SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy item references another record".to_owned(),
                        ));
                    }
                    let apply_result =
                        apply_cx_sgroup_hierarchy_progress_item(pending, hierarchies, item_index);
                    if let Err(error) = apply_result {
                        let pending = pending_sgroup_hierarchy
                            .take()
                            .expect("pending hierarchy was just checked");
                        commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                        return Err(SmartsParseError::CxSmiles(error.to_string()));
                    }
                }
                CxProgressPhase::Complete => {
                    let Some(pending) = pending_sgroup_hierarchy.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    let expected_items = hierarchies
                        .iter()
                        .map(|hierarchy| hierarchy.children.len())
                        .sum::<usize>();
                    if pending.record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || pending.next_item_index != expected_items
                        || pending.next_parent_index != hierarchies.len()
                    {
                        commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                        return Err(SmartsParseError::CxSmiles(
                            "CX SGroup hierarchy completion does not match parsed children"
                                .to_owned(),
                        ));
                    }
                    commit_pending_cx_sgroup_hierarchy(graph, pending)?;
                }
            }
            continue;
        }
        if let CxRecord::EnhancedStereo(stereo) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_stereo.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX enhanced-stereo progress has an invalid begin checkpoint"
                                .to_owned(),
                        ));
                    }
                    pending_stereo = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX enhanced-stereo item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_stereo else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX enhanced-stereo item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index
                        || *committed != item_index
                        || stereo.atoms.get(item_index).is_none()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX enhanced-stereo item checkpoints are out of order".to_owned(),
                        ));
                    }
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_stereo.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX enhanced-stereo completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != stereo.atoms.len()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX enhanced-stereo completion does not match committed indices"
                                .to_owned(),
                        ));
                    }
                    crate::cx_lowering::merge_cx_enhanced_stereo(
                        graph,
                        &mut stereo_tracker,
                        stereo,
                    )
                    .map_err(SmartsParseError::from)?;
                }
            }
            continue;
        }
        if let CxRecord::LinkNodes(nodes) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_link_nodes.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX link-node progress has an invalid begin checkpoint".to_owned(),
                        ));
                    }
                    pending_link_nodes = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX link-node item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_link_nodes else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX link-node item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index
                        || *committed != item_index
                        || nodes.get(item_index).is_none()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX link-node item checkpoints are out of order".to_owned(),
                        ));
                    }
                    // Native omitted-outer degree failure precedes comma and
                    // later syntax failure; do not write the whole property yet.
                    crate::cx_lowering::cx_link_node_outer_atoms(graph, &nodes[item_index])
                        .map_err(SmartsParseError::from)?;
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_link_nodes.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX link-node completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != nodes.len()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX link-node completion does not match committed items".to_owned(),
                        ));
                    }
                    crate::cx_lowering::apply_cx_link_nodes_to_query(graph, nodes)
                        .map_err(SmartsParseError::from)?;
                }
            }
            continue;
        }
        if matches!(
            record,
            CxRecord::Unsaturation(_) | CxRecord::RingBonds(_) | CxRecord::Substitution(_)
        ) {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if checkpoint.item_index.is_some() || pending_query_constraints.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "CX query-constraint progress has an invalid begin checkpoint"
                                .to_owned(),
                        ));
                    }
                    pending_query_constraints = Some((checkpoint.record_index, 0));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX query-constraint item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    let Some((record_index, committed)) = &mut pending_query_constraints else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX query-constraint item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index || *committed != item_index {
                        return Err(SmartsParseError::CxSmiles(
                            "CX query-constraint item checkpoints are out of order".to_owned(),
                        ));
                    }
                    crate::cx_lowering::apply_cx_query_constraint_item(graph, record, item_index)
                        .map_err(SmartsParseError::from)?;
                    *committed += 1;
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, committed)) = pending_query_constraints.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX query-constraint completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    let item_count = match record {
                        CxRecord::Unsaturation(indices) => indices.len(),
                        CxRecord::RingBonds(constraints) => constraints.len(),
                        CxRecord::Substitution(constraints) => constraints.len(),
                        _ => unreachable!("query-constraint progress kind changed"),
                    };
                    if record_index != checkpoint.record_index
                        || checkpoint.item_index.is_some()
                        || committed != item_count
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX query-constraint completion does not match committed items"
                                .to_owned(),
                        ));
                    }
                }
            }
            continue;
        }
        if let CxRecord::Coordinates(coordinates) = record {
            match checkpoint.phase {
                CxProgressPhase::Begin => {
                    if pending_coordinates.is_some() {
                        return Err(SmartsParseError::CxSmiles(
                            "nested CX coordinate progress record".to_owned(),
                        ));
                    }
                    pending_coordinates = Some((
                        checkpoint.record_index,
                        cosmolkit_cx::CxCoordinates {
                            conformer: coordinates.conformer,
                            values: Vec::new(),
                            is_3d: true,
                        },
                    ));
                }
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX coordinate item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    let Some((record_index, pending)) = &mut pending_coordinates else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX coordinate item has no begin checkpoint".to_owned(),
                        ));
                    };
                    if *record_index != checkpoint.record_index
                        || pending.values.len() != item_index
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX coordinate item checkpoints are out of order".to_owned(),
                        ));
                    }
                    let value = coordinates.values.get(item_index).copied().ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX coordinate item checkpoint references a missing row".to_owned(),
                        )
                    })?;
                    pending.values.push(value);
                }
                CxProgressPhase::Complete => {
                    let Some((record_index, mut pending)) = pending_coordinates.take() else {
                        return Err(SmartsParseError::CxSmiles(
                            "CX coordinate completion has no begin checkpoint".to_owned(),
                        ));
                    };
                    if record_index != checkpoint.record_index
                        || pending.values.len() != coordinates.values.len()
                    {
                        return Err(SmartsParseError::CxSmiles(
                            "CX coordinate completion does not match committed rows".to_owned(),
                        ));
                    }
                    pending.is_3d = coordinates.is_3d;
                    append_query_conformer(graph, &pending)?;
                }
            }
            continue;
        }
        if matches!(record, CxRecord::AtomLabels(_) | CxRecord::AtomValues(_)) {
            match checkpoint.phase {
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX atom slot checkpoint has no atom index".to_owned(),
                        )
                    })?;
                    apply_cx_atom_slot_to_query(graph, record, item_index)?;
                }
                CxProgressPhase::Complete if checkpoint.item_index.is_none() => {}
                CxProgressPhase::Begin | CxProgressPhase::Complete => {
                    return Err(SmartsParseError::CxSmiles(
                        "CX atom slot progress has an invalid checkpoint".to_owned(),
                    ));
                }
            }
            continue;
        }
        if matches!(record, CxRecord::AtomProperties(_)) {
            match checkpoint.phase {
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX atomProp item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    apply_cx_atom_property_to_query(graph, record, item_index)?;
                }
                CxProgressPhase::Complete if checkpoint.item_index.is_none() => {}
                CxProgressPhase::Begin | CxProgressPhase::Complete => {
                    return Err(SmartsParseError::CxSmiles(
                        "CX atomProp progress has an invalid checkpoint".to_owned(),
                    ));
                }
            }
            continue;
        }
        if matches!(
            record,
            CxRecord::CoordinateBonds(_) | CxRecord::ZeroBonds(_)
        ) {
            match checkpoint.phase {
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX bond item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    apply_cx_bond_item_to_query(graph, record, item_index)?;
                }
                CxProgressPhase::Complete if checkpoint.item_index.is_none() => {}
                CxProgressPhase::Begin | CxProgressPhase::Complete => {
                    return Err(SmartsParseError::CxSmiles(
                        "CX bond progress has an invalid checkpoint".to_owned(),
                    ));
                }
            }
            continue;
        }
        if matches!(record, CxRecord::Radicals(_)) {
            match checkpoint.phase {
                CxProgressPhase::Item => {
                    let item_index = checkpoint.item_index.ok_or_else(|| {
                        SmartsParseError::CxSmiles(
                            "CX radical item checkpoint has no item index".to_owned(),
                        )
                    })?;
                    apply_cx_radical_item_to_query(graph, record, item_index)?;
                }
                CxProgressPhase::Complete if checkpoint.item_index.is_none() => {}
                CxProgressPhase::Begin | CxProgressPhase::Complete => {
                    return Err(SmartsParseError::CxSmiles(
                        "CX radical progress has an invalid checkpoint".to_owned(),
                    ));
                }
            }
            continue;
        }
        match checkpoint.phase {
            CxProgressPhase::Complete if checkpoint.item_index.is_none() => {
                if pending_coordinates.is_some() {
                    return Err(SmartsParseError::CxSmiles(
                        "CX record completed before coordinate progress".to_owned(),
                    ));
                }
                apply_cx_to_query(
                    graph,
                    std::slice::from_ref(record),
                    &mut stereo_tracker,
                    &mut cx_sequence_id,
                )?;
            }
            CxProgressPhase::Begin | CxProgressPhase::Item | CxProgressPhase::Complete => {
                unreachable!("the dispatcher currently emits record-complete checkpoints only")
            }
        }
    }
    if let Some((record_index, mut pending)) = pending_coordinates {
        let Some(CxRecord::Coordinates(coordinates)) = progress.records().get(record_index) else {
            return Err(SmartsParseError::CxSmiles(
                "CX coordinate progress references a missing record".to_owned(),
            ));
        };
        if pending.values.len() != coordinates.values.len() {
            return Err(SmartsParseError::CxSmiles(
                "CX coordinate progress lost a committed row".to_owned(),
            ));
        }
        pending.is_3d = coordinates.is_3d;
        append_query_conformer(graph, &pending)?;
    }
    if let Some(pending) = pending_sgroup_hierarchy {
        commit_pending_cx_sgroup_hierarchy(graph, pending)?;
    }
    // BEGIN COMPLETE PINNED SF207 actual progress caller
    // RDKit❗❌: void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
    // RDKit❗❌:                        std::string::const_iterator &first,
    // RDKit❗❌:                        unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit❗❌:   // BOOST_LOG(rdWarningLog) << "parseCXNExtensions: " << extText << std::endl;
    // RDKit❗❌:   if (extText.empty()) {
    // RDKit❗❌:     return;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (extText[0] != '|') {
    // RDKit❗❌:     throw RDKit::SmilesParseException(
    // RDKit❗❌:         "CXSMILES extension does not start with |");
    // RDKit❗❌:   }
    // RDKit❗❌:   first = extText.begin();
    // RDKit❗❌:   bool ok =
    // RDKit❗❌:       parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
    // RDKit❗❌:   if (!ok) {
    // RDKit❗❌:     throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
    // RDKit❗❌:   }
    // RDKit❗❌:   processCXSmilesLabels(mol);
    // RDKit❗❌:   mol.clearProp("_cxsmilesLabelsProcessed");
    // RDKit❗❌:   mol.clearProp(cxsgTracker);
    // RDKit❗❌: }
    // END COMPLETE PINNED SF207 actual progress caller
    // Empty source extText returns before label/property reads. Complete
    // nonempty syntax has advanced its source iterator past the final pipe
    // before label processing and both ordered clearProp operations. Failed
    // syntax keeps earlier graph effects and bypasses the success epilogue.
    if progress.is_complete() && progress.consumed() != 0 {
        // parseCXExtensions advances past the closing pipe before its
        // successful-only CX label processing and tracker cleanup.
        *source_cursor = progress.consumed();
        crate::cx_lowering::finish_cx_smiles_labels(graph).map_err(SmartsParseError::from)?;
    }
    if !progress.is_complete() {
        // Source scanning captured a terminal warning, but only successful
        // replay reaches it. Any earlier graph/property/query failure returned
        // above before this side effect, matching native sequential parsing.
        progress.emit_source_warning();
    }
    if !progress.is_complete() {
        // Source scanning captured a terminal warning, but only successful
        // replay reaches it. Any earlier graph/property/query failure returned
        // above before this side effect, matching native sequential parsing.
        progress.emit_source_warning();
    }
    Ok(())
}

fn apply_cx_atom_slot_to_query(
    graph: &mut QueryGraph,
    record: &CxRecord,
    atom_index: usize,
) -> Result<(), SmartsParseError> {
    let (property, value) = match record {
        CxRecord::AtomLabels(values) => ("atomLabel", values.get(atom_index)),
        CxRecord::AtomValues(values) => ("molFileValue", values.get(atom_index)),
        _ => {
            return Err(SmartsParseError::CxSmiles(
                "CX atom slot checkpoint references a non-slot record".to_owned(),
            ));
        }
    };
    let Some(value) = value.and_then(Option::as_ref) else {
        return Err(SmartsParseError::CxSmiles(
            "CX atom slot checkpoint references an empty slot".to_owned(),
        ));
    };
    // RDKit: the pinned helper writes the nonempty label/value immediately
    // after `read_text_to`; `VALID_ATIDX` skips indices outside the graph.
    if let Some(atom) = graph.atom_mut(atom_index) {
        atom.set_prop(property, value)
            .map_err(SmartsParseError::from)?;
    }
    Ok(())
}

fn apply_cx_atom_property_to_query(
    graph: &mut QueryGraph,
    record: &CxRecord,
    item_index: usize,
) -> Result<(), SmartsParseError> {
    let CxRecord::AtomProperties(properties) = record else {
        return Err(SmartsParseError::CxSmiles(
            "CX atomProp checkpoint references a non-property record".to_owned(),
        ));
    };
    let property = properties.get(item_index).ok_or_else(|| {
        SmartsParseError::CxSmiles(
            "CX atomProp checkpoint references a missing property".to_owned(),
        )
    })?;
    // RDKit: valid atom indices commit nonempty properties before the helper
    // advances a following colon; out-of-range indices have no effect.
    if let Some(atom) = graph.atom_mut(property.atom) {
        atom.set_prop(property.name.clone(), property.value.clone())
            .map_err(SmartsParseError::from)?;
    }
    Ok(())
}

fn apply_cx_bond_item_to_query(
    graph: &mut QueryGraph,
    record: &CxRecord,
    item_index: usize,
) -> Result<(), SmartsParseError> {
    match record {
        CxRecord::CoordinateBonds(annotation) => {
            let reference = annotation.bonds.get(item_index).ok_or_else(|| {
                SmartsParseError::CxSmiles(
                    "CX coordinate bond checkpoint references a missing item".to_owned(),
                )
            })?;
            apply_cx_coordinate_bond_to_query(graph, *reference, annotation.kind)
        }
        CxRecord::ZeroBonds(indices) => {
            let index = indices.get(item_index).copied().ok_or_else(|| {
                SmartsParseError::CxSmiles(
                    "CX zero-bond checkpoint references a missing item".to_owned(),
                )
            })?;
            apply_cx_zero_bond_to_query(graph, index)
        }
        _ => Err(SmartsParseError::CxSmiles(
            "CX bond item checkpoint references a non-bond record".to_owned(),
        )),
    }
}

fn apply_cx_coordinate_bond_to_query(
    graph: &mut QueryGraph,
    reference: cosmolkit_cx::CxBondReference,
    kind: CxCoordinateBondKind,
) -> Result<(), SmartsParseError> {
    // BEGIN COMPLETE PINNED SF189
    // RDKit✔️❌: bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                             Bond::BondType typ, unsigned int startAtomIdx,
    // RDKit✔️❌:                             unsigned int startBondIdx) {
    // RDKit✔️❌:   if (first >= last || (*first != 'C' && *first != 'H')) {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   if (first >= last || *first != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   while (first <= last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int aidx;
    // RDKit✔️❌:     unsigned int bidx;
    // RDKit✔️❌:     if (read_int_pair(first, last, aidx, bidx)) {
    // RDKit✔️❌:       if (VALID_ATIDX(aidx) && VALID_BNDIDX(bidx)) {
    // RDKit✔️❌:         auto bnd = get_bond_with_smiles_idx(mol, bidx - startBondIdx);
    // RDKit✔️❌:         if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
    // RDKit✔️❌:                      bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
    // RDKit✔️❌:           BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
    // RDKit✔️❌:                                   << " involving atom " << aidx << std::endl;
    // RDKit✔️❌:           return false;
    // RDKit✔️❌:         }
    // RDKit✔️❌:         bnd->setBondType(typ);
    // RDKit✔️❌:         if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
    // RDKit✔️❌:           unsigned int tmp = bnd->getBeginAtomIdx();
    // RDKit✔️❌:           bnd->setBeginAtomIdx(aidx - startAtomIdx);
    // RDKit✔️❌:           bnd->setEndAtomIdx(tmp);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     } else {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF189
    // The same query destination body owns complete-record and checkpoint
    // writes. This adapter only translates its structured error; the caller
    // already retained the pair's source cursor before invoking it.
    crate::cx_lowering::apply_cx_coordinate_bond(graph, reference, kind)
        .map_err(SmartsParseError::from)
}

fn apply_cx_zero_bond_to_query(
    graph: &mut QueryGraph,
    index: usize,
) -> Result<(), SmartsParseError> {
    // BEGIN COMPLETE PINNED SF190
    // RDKit✔️❌: bool parse_zero_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                       unsigned int, unsigned int startBondIdx) {
    // RDKit✔️❌:   // these look like: C1CCCCC~CCCC1 |Z:5|
    // RDKit✔️❌:   if (first >= last || *first != 'Z') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   if (first >= last || *first != ':') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:
    // RDKit✔️❌:   while (first < last && *first >= '0' && *first <= '9') {
    // RDKit✔️❌:     unsigned int bondIdx;
    // RDKit✔️❌:     if (!read_int(first, last, bondIdx)) {
    // RDKit✔️❌:       return false;
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (VALID_BNDIDX(bondIdx)) {
    // RDKit✔️❌:       auto bond = get_bond_with_smiles_idx(mol, bondIdx - startBondIdx);
    // RDKit✔️❌:
    // RDKit✔️❌:       if (!bond) {
    // RDKit✔️❌:         BOOST_LOG(rdWarningLog)
    // RDKit✔️❌:             << "bond " << bondIdx
    // RDKit✔️❌:             << " not found, cannot mark as zero order bond." << std::endl;
    // RDKit✔️❌:         return false;
    // RDKit✔️❌:       }
    // RDKit✔️❌:       bond->setBondType(Bond::ZERO);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     if (first < last && *first == ',') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF190
    // Ordered checkpoint replay already retained this index's source cursor.
    // Use the same destination body as complete records; translate only error.
    crate::cx_lowering::apply_cx_zero_bond(graph, index).map_err(SmartsParseError::from)
}

fn apply_cx_radical_item_to_query(
    graph: &mut QueryGraph,
    record: &CxRecord,
    item_index: usize,
) -> Result<(), SmartsParseError> {
    let CxRecord::Radicals(radicals) = record else {
        return Err(SmartsParseError::CxSmiles(
            "CX radical checkpoint references a non-radical record".to_owned(),
        ));
    };
    let radical = radicals.get(item_index).ok_or_else(|| {
        SmartsParseError::CxSmiles("CX radical checkpoint references a missing item".to_owned())
    })?;
    apply_cx_radical_to_query(graph, *radical)
}

fn apply_cx_radical_to_query(
    graph: &mut QueryGraph,
    radical: cosmolkit_cx::CxRadical,
) -> Result<(), SmartsParseError> {
    // Canonical SEARCH effect is shared with complete records. This private
    // parser adapter translates the domain error only; no second field effect.
    crate::cx_lowering::apply_cx_radical_to_query(graph, radical).map_err(SmartsParseError::from)
}

fn append_query_conformer(
    graph: &mut QueryGraph,
    coordinates: &cosmolkit_cx::CxCoordinates,
) -> Result<(), SmartsParseError> {
    let mut values = vec![[0.0; 3]; graph.num_atoms()];
    // BEGIN COMPLETE PINNED SF188
    // RDKit✔️❌: bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
    // RDKit✔️❌:                   unsigned int startAtomIdx, unsigned int confIdx) {
    // RDKit✔️❌:   if (first >= last || *first != '(') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   auto *conf = new Conformer(mol.getNumAtoms());
    // RDKit✔️❌:   mol.addConformer(conf);
    // RDKit✔️❌:   conf->setId(confIdx);
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   unsigned int atIdx = 0;
    // RDKit✔️❌:   bool is3D = false;
    // RDKit✔️❌:   while (first <= last && *first != ')') {
    // RDKit✔️❌:     RDGeom::Point3D pt;
    // RDKit✔️❌:     std::string tkn = read_text_to(first, last, ";)");
    // RDKit✔️❌:     if (VALID_ATIDX(atIdx)) {
    // RDKit✔️❌:       if (!tkn.empty()) {
    // RDKit✔️❌:         std::vector<std::string> tokens;
    // RDKit✔️❌:         boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
    // RDKit✔️❌:         if (tokens.size() >= 1 && tokens[0].size()) {
    // RDKit✔️❌:           pt.x = boost::lexical_cast<double>(tokens[0]);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (tokens.size() >= 2 && tokens[1].size()) {
    // RDKit✔️❌:           pt.y = boost::lexical_cast<double>(tokens[1]);
    // RDKit✔️❌:         }
    // RDKit✔️❌:         if (tokens.size() >= 3 && tokens[2].size()) {
    // RDKit✔️❌:           pt.z = boost::lexical_cast<double>(tokens[2]);
    // RDKit✔️❌:           is3D = true;
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:
    // RDKit✔️❌:       conf->setAtomPos(atIdx - startAtomIdx, pt);
    // RDKit✔️❌:     }
    // RDKit✔️❌:     ++atIdx;
    // RDKit✔️❌:     if (first <= last && *first != ')') {
    // RDKit✔️❌:       ++first;
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   // make sure that the conformer really is 3D!
    // RDKit✔️❌:   if (is3D && hasNonZeroZCoords(*conf)) {
    // RDKit✔️❌:     conf->set3D(true);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     conf->set3D(false);
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (first >= last || *first != ')') {
    // RDKit✔️❌:     return false;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   ++first;
    // RDKit✔️❌:   return true;
    // RDKit✔️❌: }
    // END COMPLETE PINNED SF188
    let atom_count = u32::try_from(graph.num_atoms()).map_err(|_| {
        SmartsParseError::CxSmiles("CX atom count exceeds source unsigned32 domain".to_owned())
    })?;
    for (slot, source) in coordinates.values.iter().enumerate() {
        let index = slot as u32;
        if index < atom_count {
            values[index as usize] = source.unwrap_or([0.0; 3]);
        }
    }
    graph
        .add_conformer_3d(cosmolkit_model::Conformer3D::new(
            coordinates.conformer,
            values,
            coordinates.is_3d,
        ))
        .map_err(SmartsParseError::from)
}

pub fn parse_smarts(
    smarts: impl AsRef<[u8]>,
    params: &SmartsParseParams,
) -> Result<QueryGraph, SmartsParseError> {
    let smarts = smarts.as_ref();
    // BEGIN RDKIT CPP FUNCTION MolFromSmarts
    // RDKit❗❌: std::unique_ptr<RWMol> MolFromSmarts(const std::string &smarts,
    // RDKit❗❌:                                      const SmartsParserParams &params) {
    // RDKit❗❌:   // Calling MolFromSmarts in a multithreaded context is generally safe *unless*
    // RDKit❗❌:   // the value of debugParse is different for different threads. The if
    // RDKit❗❌:   // statement below avoids a TSAN warning in the case where multiple threads
    // RDKit❗❌:   // all use the same value for debugParse.
    // RDKit❌❌:   if (yysmarts_debug != params.debugParse) {
    // RDKit❌❌:     yysmarts_debug = params.debugParse;
    // RDKit❗❌:   }
    // RDKit❗❌:
    // RDKit❗❌:   std::string lsmarts, name, cxPart;
    // RDKit❗❌:   preprocessSmiles(smarts, params, lsmarts, name, cxPart);
    // RDKit❗❌:
    // RDKit❗❌:   auto res = toMol(labelRecursivePatterns(lsmarts), smarts_parse, lsmarts);
    // RDKit❗❌:   handleCXPartAndName(res.get(), params, cxPart, name);
    // RDKit❗❌:   if (res) {
    // RDKit❗❌:     if (params.mergeHs) {
    // RDKit❗❌:       MolOps::mergeQueryHs(*res);
    // RDKit❗❌:     }
    // RDKit❗❌:     MolOps::setBondStereoFromDirections(*res);
    // RDKit❗❌:     if (!params.skipCleanup) {
    // RDKit❗❌:       SmilesParseOps::CleanupAfterParsing(res.get());
    // RDKit❗❌:     }
    // RDKit❗❌:     if (!name.empty()) {
    // RDKit❗❌:       res->setProp(common_properties::_Name, name);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return res;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION MolFromSmarts
    // Local complexity review: current lexer borrows each remaining tail;
    // token buffering remains an extra O(n) allocation. CX source-index
    // lookup scans O(E) properties per reference, matching RDKit. Query-H
    // merging compacts the detached graph; no MoleculeBuilder is used here.
    if params.debug_parse {
        return Err(SmartsParseError::UnsupportedFeature(
            "Bison debug_parse diagnostic output",
        ));
    }
    let preprocessed = preprocess_smarts(smarts, params);
    let labeled = label_recursive_patterns(&preprocessed.smarts);
    // Lower parser state to the canonical query value before any source
    // post-processing. All production post-processing below stays in the
    // query model; test-only compatibility adapters are outside this path.
    // Preserve typed reached property/carrier errors through the root entry.
    let mut parsed_graph = smarts_parse_entry(&labeled)?.finish()?;
    // BEGIN RDKIT CPP FUNCTION handleCXPartAndName
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
    // END RDKIT CPP FUNCTION handleCXPartAndName
    // BEGIN RDKIT CPP FUNCTION parseCXExtensions
    // RDKit✔️❌: void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
    // RDKit✔️❌:                        std::string::const_iterator &first,
    // RDKit✔️❌:                        unsigned int startAtomIdx, unsigned int startBondIdx) {
    // RDKit✔️❌:   // BOOST_LOG(rdWarningLog) << "parseCXNExtensions: " << extText << std::endl;
    // RDKit✔️❌:   if (extText.empty()) {
    // RDKit✔️❌:     return;
    // RDKit✔️❌:   }
    // RDKit✔️❌:   if (extText[0] != '|') {
    // RDKit✔️❌:     throw RDKit::SmilesParseException(
    // RDKit✔️❌:         "CXSMILES extension does not start with |");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   first = extText.begin();
    // RDKit✔️❌:   bool ok =
    // RDKit✔️❌:       parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
    // RDKit✔️❌:   if (!ok) {
    // RDKit✔️❌:     throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
    // RDKit✔️❌:   }
    // RDKit✔️❌:   processCXSmilesLabels(mol);
    // RDKit✔️❌:   mol.clearProp("_cxsmilesLabelsProcessed");
    // RDKit✔️❌:   mol.clearProp(cxsgTracker);
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION parseCXExtensions
    // Behavior: empty CX bypasses every CX/name read. With CX disabled, the
    // entire nonempty suffix becomes the name only when parseName is enabled.
    // Native catch(...) retains earlier effects and its actual iterator under
    // strict=false, writes the consumed raw prefix, and skips suffix naming.
    // Strict=true propagates reached typed errors before prefix/name writes.
    // A syntax false return uses the native wrapper's fixed exception text.
    // Successful parsing processes labels and clears both markers before the
    // prefix property is written; only then does source trimming assign name.
    // Local cost: prefix/name each copy their counted source bytes once. The
    // reached progress parser buffers records and checkpoints before replay,
    // unlike native streaming mutation, so the known buffering cost is ❌.
    // No clone of the whole query graph, error decoding or heuristic recovery.
    let mut name = preprocessed.name;
    if !preprocessed.cx_part.is_empty() {
        let mut cx_failed = false;
        let mut consumed = 0;
        if params.allow_cxsmiles {
            if preprocessed.cx_part.starts_with(b"|") {
                let atom_count = u32::try_from(parsed_graph.num_atoms()).map_err(|_| {
                    SmartsParseError::CxSmiles(
                        "CX atom count exceeds source unsigned32 domain".to_owned(),
                    )
                })?;
                let progress = cosmolkit_cx::parse_cx_extensions_progress_with_atom_window(
                    &preprocessed.cx_part,
                    0,
                    atom_count,
                );
                consumed = progress.consumed();
                let mut lowering_cursor = consumed;
                if let Err(error) = apply_cx_progress_to_query_with_cursor(
                    &mut parsed_graph,
                    &progress,
                    &mut lowering_cursor,
                ) {
                    cx_failed = true;
                    if params.strict_cxsmiles {
                        return Err(error);
                    }
                    consumed = lowering_cursor;
                } else if !progress.is_complete() {
                    cx_failed = true;
                    if params.strict_cxsmiles {
                        // parse_it returned false: the native wrapper throws
                        // this fixed text, rather than its internal diagnostic.
                        return Err(SmartsParseError::CxSmiles(
                            "failure parsing CXSMILES extensions".to_owned(),
                        ));
                    }
                }
                let prefix = preprocessed.cx_part.get(..consumed).ok_or_else(|| {
                    SmartsParseError::CxSmiles(
                        "CX source cursor exceeds the counted extension".to_owned(),
                    )
                })?;
                parsed_graph.set_prop("_CXSMILES_Data", PropertyText::from_bytes(prefix))?;
            } else if params.strict_cxsmiles && !params.parse_name {
                return Err(SmartsParseError::CxSmiles(
                    "CXSMILES extension does not start with | and parseName=false".to_owned(),
                ));
            }
        }
        if !cx_failed && params.parse_name && consumed < preprocessed.cx_part.len() {
            // Native assignment also replaces the local name with an empty
            // trimmed suffix; only the later _Name assignment tests emptiness.
            name = trim_source_whitespace(&preprocessed.cx_part[consumed..]).to_owned();
        }
    }
    if params.merge_hs {
        merge_query_hs_in_place(&mut parsed_graph, false, false)?;
    }
    crate::query_graph_behavior::set_bond_stereo_from_directions(&mut parsed_graph)?;
    if !params.skip_cleanup {
        crate::query_graph_behavior::cleanup_query_graph_parser_state(&mut parsed_graph)?;
    }
    if !name.is_empty() {
        parsed_graph = parsed_graph.with_name(PropertyText::from(name));
    }
    Ok(parsed_graph)
}

fn smarts_parse_helper(input: impl AsRef<[u8]>) -> Result<QueryGraphBuilder, SmartsParseError> {
    let input = input.as_ref();
    // RDKit✔️✔️: int smarts_parse_helper(const std::string &inp, ...,
    // RDKit✔️✔️:   return generic_parse_helper<yysmarts_lex_init,
    // RDKit✔️✔️:     setup_smarts_string, yysmarts_lex_destroy>(yysmarts_parse,
    // RDKit✔️✔️:     inp, molVect, atom, bond, start_tok, "SMARTS");
    // Local complexity review: preprocessing has already produced one input
    // buffer. Recursive descent follows source reductions, but the lexer
    // still buffers scanned and compacted tokens and recursive subparses;
    // raw suffix lookup removes the former O(n^2) suffix copying.
    let tokens = tokenize(input)?;
    let mut parser = SmartsParser::new(&tokens, input);
    parser.parse_smarts_molecule()
}

fn parse_atom_entry(
    input: impl AsRef<[u8]>,
) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
    let input = input.as_ref();
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
    Ok(atom.carrier.predicate().clone())
}

fn parse_bond_entry(
    input: impl AsRef<[u8]>,
) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
    let input = input.as_ref();
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
    Ok(bond.query)
}

fn smarts_bond_parse(
    input: impl AsRef<[u8]>,
) -> Result<QueryNode<BondQueryPredicate>, SmartsParseError> {
    let input = input.as_ref();
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

fn smarts_atom_parse(
    input: impl AsRef<[u8]>,
) -> Result<QueryNode<AtomQueryPredicate>, SmartsParseError> {
    let input = input.as_ref();
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

fn to_atom(
    inp: impl AsRef<[u8]>,
) -> Result<Option<QueryNode<AtomQueryPredicate>>, SmartsParseError> {
    let inp = inp.as_ref();
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
    smarts: impl AsRef<[u8]>,
) -> Result<Option<QueryNode<AtomQueryPredicate>>, SmartsParseError> {
    let smarts = smarts.as_ref();
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

fn to_bond(
    inp: impl AsRef<[u8]>,
) -> Result<Option<QueryNode<BondQueryPredicate>>, SmartsParseError> {
    let inp = inp.as_ref();
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
    smarts: impl AsRef<[u8]>,
) -> Result<Option<QueryNode<BondQueryPredicate>>, SmartsParseError> {
    let smarts = smarts.as_ref();
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

fn smarts_parse_entry(input: impl AsRef<[u8]>) -> Result<QueryGraphBuilder, SmartsParseError> {
    let input = input.as_ref();
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
    let input = input.as_ref();
    if input.is_empty() {
        return Ok(QueryGraphBuilder::default());
    }
    smarts_parse_helper(input)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct PreprocessedSmarts {
    smarts: Vec<u8>,
    name: Vec<u8>,
    cx_part: Vec<u8>,
}

#[doc(hidden)]
pub fn trim_source_whitespace(value: &[u8]) -> &[u8] {
    // Boost✔️✔️: inline SequenceT trim_copy( const SequenceT& Input, const std::locale& Loc=std::locale() )
    // Boost✔️✔️:         {
    // Boost✔️✔️:             return
    // Boost✔️✔️:                 ::boost::algorithm::trim_copy_if(
    // Boost✔️✔️:                     Input,
    // Boost✔️✔️:                     is_space(Loc) );
    // Boost✔️✔️:         }
    // Boost✔️✔️: inline SequenceT trim_copy_if(const SequenceT& Input, PredicateT IsSpace)
    // Boost✔️✔️:         {
    // Boost✔️✔️:             BOOST_STRING_TYPENAME
    // Boost✔️✔️:                 range_const_iterator<SequenceT>::type TrimEnd=
    // Boost✔️✔️:                     ::boost::algorithm::detail::trim_end(
    // Boost✔️✔️:                         ::boost::begin(Input),
    // Boost✔️✔️:                         ::boost::end(Input),
    // Boost✔️✔️:                         IsSpace);
    // Boost✔️✔️:
    // Boost✔️✔️:             return SequenceT(
    // Boost✔️✔️:                 detail::trim_begin(
    // Boost✔️✔️:                     ::boost::begin(Input),
    // Boost✔️✔️:                     TrimEnd,
    // Boost✔️✔️:                     IsSpace),
    // Boost✔️✔️:                 TrimEnd
    // Boost✔️✔️:                 );
    // Boost✔️✔️:         }
    // Boost✔️✔️: inline ForwardIteratorT trim_begin(
    // Boost✔️✔️:                 ForwardIteratorT InBegin,
    // Boost✔️✔️:                 ForwardIteratorT InEnd,
    // Boost✔️✔️:                 PredicateT IsSpace )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 ForwardIteratorT It=InBegin;
    // Boost✔️✔️:                 for(; It!=InEnd; ++It )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     if (!IsSpace(*It))
    // Boost✔️✔️:                         return It;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return It;
    // Boost✔️✔️:             }
    // Boost✔️✔️: inline ForwardIteratorT trim_end(
    // Boost✔️✔️:                 ForwardIteratorT InBegin,
    // Boost✔️✔️:                 ForwardIteratorT InEnd,
    // Boost✔️✔️:                 PredicateT IsSpace )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️✔️:                     std::iterator_traits<ForwardIteratorT>::iterator_category category;
    // Boost✔️✔️:
    // Boost✔️✔️:                 return ::boost::algorithm::detail::trim_end_iter_select( InBegin, InEnd, IsSpace, category() );
    // Boost✔️✔️:             }
    // Boost✔️✔️: inline ForwardIteratorT trim_end_iter_select(
    // Boost✔️✔️:                 ForwardIteratorT InBegin,
    // Boost✔️✔️:                 ForwardIteratorT InEnd,
    // Boost✔️✔️:                 PredicateT IsSpace,
    // Boost✔️✔️:                 std::bidirectional_iterator_tag )
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 for( ForwardIteratorT It=InEnd; It!=InBegin;  )
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     if ( !IsSpace(*(--It)) )
    // Boost✔️✔️:                         return ++It;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 return InBegin;
    // Boost✔️✔️:             }
    // Boost✔️✔️: struct is_classifiedF :
    // Boost✔️✔️:                 public predicate_facade<is_classifiedF>
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 // Boost.ResultOf support
    // Boost✔️✔️:                 typedef bool result_type;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Constructor from a locale
    // Boost✔️✔️:                 is_classifiedF(std::ctype_base::mask Type, std::locale const & Loc = std::locale()) :
    // Boost✔️✔️:                     m_Type(Type), m_Locale(Loc) {}
    // Boost✔️✔️:                 // Operation
    // Boost✔️✔️:                 template<typename CharT>
    // Boost✔️✔️:                 bool operator()( CharT Ch ) const
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     return std::use_facet< std::ctype<CharT> >(m_Locale).is( m_Type, Ch );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:                 #if defined(BOOST_BORLANDC) && (BOOST_BORLANDC >= 0x560) && (BOOST_BORLANDC <= 0x582) && !defined(_USE_OLD_RW_STL)
    // Boost✔️✔️:                     template<>
    // Boost✔️✔️:                     bool operator()( char const Ch ) const
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         return std::use_facet< std::ctype<char> >(m_Locale).is( m_Type, Ch );
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 #endif
    // Boost✔️✔️:
    // Boost✔️✔️:             private:
    // Boost✔️✔️:                 std::ctype_base::mask m_Type;
    // Boost✔️✔️:                 std::locale m_Locale;
    // Boost✔️✔️:             }
    // Exact pinned C-locale ctype<char> space bytes, without decoding. Boost
    // obtains TrimEnd first, then searches TrimBegin only before that bound.
    // Each edge scan is linear in trimmed bytes with constant borrowed storage;
    // the caller owns the single counted copy, like trim_copy's constructor.
    let end = match value
        .iter()
        .rposition(|b| !matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 12 | 11))
    {
        Some(last) => last + 1,
        None => 0,
    };
    let start = value[..end]
        .iter()
        .position(|b| !matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 12 | 11))
        .unwrap_or(end);
    &value[start..end]
}

fn preprocess_smarts(smarts: impl AsRef<[u8]>, params: &SmartsParseParams) -> PreprocessedSmarts {
    let smarts = smarts.as_ref();
    // BEGIN RDKIT CPP FUNCTION preprocessSmiles<SmartsParserParams>
    // RDKit✔️❌: // despite the name: works for both SMILES and SMARTS
    // RDKit✔️❌: template <typename T>
    // RDKit✔️❌: void preprocessSmiles(const std::string &smiles, const T &params,
    // RDKit✔️❌:                       std::string &lsmiles, std::string &name,
    // RDKit✔️❌:                       std::string &cxPart) {
    // RDKit✔️❌:   cxPart = "";
    // RDKit✔️❌:   name = "";
    // RDKit✔️❌:   if (params.parseName && !params.allowCXSMILES) {
    // RDKit✔️❌:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️❌:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️❌:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️❌:       name = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️❌:     }
    // RDKit✔️❌:   } else if (params.allowCXSMILES) {
    // RDKit✔️❌:     size_t sidx = smiles.find_first_of(" \t");
    // RDKit✔️❌:     if (sidx != std::string::npos && sidx != 0) {
    // RDKit✔️❌:       lsmiles = smiles.substr(0, sidx);
    // RDKit✔️❌:       cxPart = boost::trim_copy(smiles.substr(sidx, smiles.size() - sidx));
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (lsmiles.empty()) {
    // RDKit✔️❌:     lsmiles = smiles;
    // RDKit✔️❌:   }
    // RDKit✔️❌:
    // RDKit✔️❌:   if (!params.replacements.empty()) {
    // RDKit✔️❌:     std::string smi = lsmiles;
    // RDKit✔️❌:     for (auto loopAgain = true; loopAgain;) {
    // RDKit✔️❌:       loopAgain = false;
    // RDKit✔️❌:       for (const auto &pr : params.replacements) {
    // RDKit✔️❌:         if (smi.find(pr.first) != std::string::npos) {
    // RDKit✔️❌:           loopAgain = true;
    // RDKit✔️❌:           boost::replace_all(smi, pr.first, pr.second);
    // RDKit✔️❌:         }
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:     lsmiles = smi;
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // END RDKIT CPP FUNCTION preprocessSmiles<SmartsParserParams>
    // BEGIN COMPLETE REACHED replace_all
    // Boost✔️❌: inline void replace_all(
    // Boost✔️❌:             SequenceT& Input,
    // Boost✔️❌:             const Range1T& Search,
    // Boost✔️❌:             const Range2T& Format )
    // Boost✔️❌:         {
    // Boost✔️❌:             ::boost::algorithm::find_format_all(
    // Boost✔️❌:                 Input,
    // Boost✔️❌:                 ::boost::algorithm::first_finder(Search),
    // Boost✔️❌:                 ::boost::algorithm::const_formatter(Format) );
    // Boost✔️❌:         }
    // END COMPLETE REACHED replace_all
    // BEGIN COMPLETE REACHED first_finderF
    // Boost✔️✔️: struct first_finderF
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 typedef SearchIteratorT search_iterator_type;
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Construction
    // Boost✔️✔️:                 template< typename SearchT >
    // Boost✔️✔️:                 first_finderF( const SearchT& Search, PredicateT Comp ) :
    // Boost✔️✔️:                     m_Search(::boost::begin(Search), ::boost::end(Search)), m_Comp(Comp) {}
    // Boost✔️✔️:                 first_finderF(
    // Boost✔️✔️:                         search_iterator_type SearchBegin,
    // Boost✔️✔️:                         search_iterator_type SearchEnd,
    // Boost✔️✔️:                         PredicateT Comp ) :
    // Boost✔️✔️:                     m_Search(SearchBegin, SearchEnd), m_Comp(Comp) {}
    // Boost✔️✔️:
    // Boost✔️✔️:                 // Operation
    // Boost✔️✔️:                 template< typename ForwardIteratorT >
    // Boost✔️✔️:                 iterator_range<ForwardIteratorT>
    // Boost✔️✔️:                 operator()(
    // Boost✔️✔️:                     ForwardIteratorT Begin,
    // Boost✔️✔️:                     ForwardIteratorT End ) const
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     typedef iterator_range<ForwardIteratorT> result_type;
    // Boost✔️✔️:                     typedef ForwardIteratorT input_iterator_type;
    // Boost✔️✔️:
    // Boost✔️✔️:                     // Outer loop
    // Boost✔️✔️:                     for(input_iterator_type OuterIt=Begin;
    // Boost✔️✔️:                         OuterIt!=End;
    // Boost✔️✔️:                         ++OuterIt)
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         // Sanity check
    // Boost✔️✔️:                         if( boost::empty(m_Search) )
    // Boost✔️✔️:                             return result_type( End, End );
    // Boost✔️✔️:
    // Boost✔️✔️:                         input_iterator_type InnerIt=OuterIt;
    // Boost✔️✔️:                         search_iterator_type SubstrIt=m_Search.begin();
    // Boost✔️✔️:                         for(;
    // Boost✔️✔️:                             InnerIt!=End && SubstrIt!=m_Search.end();
    // Boost✔️✔️:                             ++InnerIt,++SubstrIt)
    // Boost✔️✔️:                         {
    // Boost✔️✔️:                             if( !( m_Comp(*InnerIt,*SubstrIt) ) )
    // Boost✔️✔️:                                 break;
    // Boost✔️✔️:                         }
    // Boost✔️✔️:
    // Boost✔️✔️:                         // Substring matching succeeded
    // Boost✔️✔️:                         if ( SubstrIt==m_Search.end() )
    // Boost✔️✔️:                             return result_type( OuterIt, InnerIt );
    // Boost✔️✔️:                     }
    // Boost✔️✔️:
    // Boost✔️✔️:                     return result_type( End, End );
    // Boost✔️✔️:                 }
    // Boost✔️✔️:
    // Boost✔️✔️:             private:
    // Boost✔️✔️:                 iterator_range<search_iterator_type> m_Search;
    // Boost✔️✔️:                 PredicateT m_Comp;
    // Boost✔️✔️:             }
    // END COMPLETE REACHED first_finderF
    // BEGIN COMPLETE REACHED nonoverlapping in-place find_format_all_impl2
    // Boost✔️❌: inline void find_format_all_impl2(
    // Boost✔️❌:                 InputT& Input,
    // Boost✔️❌:                 FinderT Finder,
    // Boost✔️❌:                 FormatterT Formatter,
    // Boost✔️❌:                 FindResultT FindResult,
    // Boost✔️❌:                 FormatResultT FormatResult)
    // Boost✔️❌:             {
    // Boost✔️❌:                 typedef BOOST_STRING_TYPENAME
    // Boost✔️❌:                     range_iterator<InputT>::type input_iterator_type;
    // Boost✔️❌:                 typedef find_format_store<
    // Boost✔️❌:                         input_iterator_type,
    // Boost✔️❌:                         FormatterT,
    // Boost✔️❌:                         FormatResultT > store_type;
    // Boost✔️❌:
    // Boost✔️❌:                 // Create store for the find result
    // Boost✔️❌:                 store_type M( FindResult, FormatResult, Formatter );
    // Boost✔️❌:
    // Boost✔️❌:                 // Instantiate replacement storage
    // Boost✔️❌:                 std::deque<
    // Boost✔️❌:                     BOOST_STRING_TYPENAME range_value<InputT>::type> Storage;
    // Boost✔️❌:
    // Boost✔️❌:                 // Initialize replacement iterators
    // Boost✔️❌:                 input_iterator_type InsertIt=::boost::begin(Input);
    // Boost✔️❌:                 input_iterator_type SearchIt=::boost::begin(Input);
    // Boost✔️❌:
    // Boost✔️❌:                 while( M )
    // Boost✔️❌:                 {
    // Boost✔️❌:                     // process the segment
    // Boost✔️❌:                     InsertIt=process_segment(
    // Boost✔️❌:                         Storage,
    // Boost✔️❌:                         Input,
    // Boost✔️❌:                         InsertIt,
    // Boost✔️❌:                         SearchIt,
    // Boost✔️❌:                         M.begin() );
    // Boost✔️❌:
    // Boost✔️❌:                     // Adjust search iterator
    // Boost✔️❌:                     SearchIt=M.end();
    // Boost✔️❌:
    // Boost✔️❌:                     // Copy formatted replace to the storage
    // Boost✔️❌:                     ::boost::algorithm::detail::copy_to_storage( Storage, M.format_result() );
    // Boost✔️❌:
    // Boost✔️❌:                     // Find range for a next match
    // Boost✔️❌:                     M=Finder( SearchIt, ::boost::end(Input) );
    // Boost✔️❌:                 }
    // Boost✔️❌:
    // Boost✔️❌:                 // process the last segment
    // Boost✔️❌:                 InsertIt=::boost::algorithm::detail::process_segment(
    // Boost✔️❌:                     Storage,
    // Boost✔️❌:                     Input,
    // Boost✔️❌:                     InsertIt,
    // Boost✔️❌:                     SearchIt,
    // Boost✔️❌:                     ::boost::end(Input) );
    // Boost✔️❌:
    // Boost✔️❌:                 if ( Storage.empty() )
    // Boost✔️❌:                 {
    // Boost✔️❌:                     // Truncate input
    // Boost✔️❌:                     ::boost::algorithm::detail::erase( Input, InsertIt, ::boost::end(Input) );
    // Boost✔️❌:                 }
    // Boost✔️❌:                 else
    // Boost✔️❌:                 {
    // Boost✔️❌:                     // Copy remaining data to the end of input
    // Boost✔️❌:                     ::boost::algorithm::detail::insert( Input, ::boost::end(Input), Storage.begin(), Storage.end() );
    // Boost✔️❌:                 }
    // Boost✔️❌:             }
    // END COMPLETE REACHED nonoverlapping in-place find_format_all_impl2
    // Behavior: source split recognizes space/tab only and ignores a split
    // at byte zero. Name or CX suffix is trimmed with the C-locale byte class;
    // an empty core after the split uses the complete original counted bytes.
    // PropertyText map order is unsigned-byte lexicographic, matching native
    // std::string compare. Every rule checks for a source match, then replaces
    // nonoverlapping matches of that pass's original buffer. A later map rule
    // sees earlier output; matching any rule repeats the ordered map. Empty
    // search sets loopAgain exactly like std::string::find, but Boost's finder
    // has no empty match, so native repetition is retained without a guard.
    // NUL/high bytes remain input data. No UTF8 validation, replacement cycle
    // heuristic, iteration limit, success rejection or consumer preprocessing.
    // Cost: same ordered fixed-point scans and substring comparisons, but the
    // Rust matched rule allocates a whole output Vec while Boost mutates the
    // input and buffers insertion overflow in a deque. That extra output
    // allocation/copy is a known memory/material constant cost, marked ❌.
    let mut processed = PreprocessedSmarts {
        smarts: Vec::new(),
        name: Vec::new(),
        cx_part: Vec::new(),
    };
    if params.parse_name && !params.allow_cxsmiles {
        if let Some(split_index) = smarts
            .iter()
            .copied()
            .position(|byte| matches!(byte, b' ' | b'\t'))
            && split_index != 0
        {
            processed.smarts = smarts[..split_index].to_vec();
            processed.name = trim_source_whitespace(&smarts[split_index..]).to_vec();
        }
    } else if params.allow_cxsmiles
        && let Some(split_index) = smarts
            .iter()
            .copied()
            .position(|byte| matches!(byte, b' ' | b'\t'))
        && split_index != 0
    {
        processed.smarts = smarts[..split_index].to_vec();
        processed.cx_part = trim_source_whitespace(&smarts[split_index..]).to_vec();
    }

    if processed.smarts.is_empty() {
        processed.smarts = smarts.to_vec();
    }

    if !params.replacements.is_empty() {
        loop {
            let mut loop_again = false;
            for (key, value) in &params.replacements {
                let key = key.as_bytes();
                let value = value.as_bytes();
                // std::string::find(empty) succeeds; boost replace_all has no
                // match for an empty search range. Source still loops again.
                if key.is_empty() {
                    loop_again = true;
                    continue;
                }
                if processed
                    .smarts
                    .windows(key.len())
                    .any(|window| window == key)
                {
                    loop_again = true;
                    let mut replaced = Vec::new();
                    let mut cursor = 0;
                    while cursor < processed.smarts.len() {
                        if processed.smarts[cursor..].starts_with(key) {
                            replaced.extend_from_slice(value);
                            cursor += key.len();
                        } else {
                            replaced.push(processed.smarts[cursor]);
                            cursor += 1;
                        }
                    }
                    processed.smarts = replaced;
                }
            }
            if !loop_again {
                break;
            }
        }
    }
    processed
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

fn is_query_hydrogen(atom: &QueryAtom, degree: usize) -> QueryHydrogenType {
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
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_))
    );
    if atom.atomic_number() == 1 && root_atomic_number_query {
        return QueryHydrogenType::QueryHydrogen;
    }
    if degree > 1 || matches!(atom.predicate(), QueryNode::Not(_)) {
        return QueryHydrogenType::NotAHydrogen;
    }
    if let query @ (QueryNode::Or(_) | QueryNode::And(_)) = atom.predicate() {
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

fn merge_recursive_query_hydrogens(
    query: &mut QueryNode<AtomQueryPredicate>,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<(), SmartsParseError> {
    // RDKit❗❌:       // recurse if needed (was github isusue 544)
    // RDKit❗❌:       if (atom->hasQuery()) {
    // RDKit❗❌:         if (atom->getQuery()->getDescription() == "RecursiveStructure") {
    // RDKit❗❌:           auto *rsq = dynamic_cast<RecursiveStructureQuery *>(atom->getQuery());
    // RDKit❗❌:           CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit❗❌:           RWMol *rqm = new RWMol(*rsq->getQueryMol());
    // RDKit❗❌:           mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
    // RDKit❗❌:           rsq->setQueryMol(rqm);
    // RDKit❗❌:         }
    // RDKit❗❌:
    // RDKit❗❌:         // FIX: shouldn't be repeating this code here
    // RDKit❗❌:         std::list<QueryAtom::QUERYATOM_QUERY::CHILD_TYPE> childStack(
    // RDKit❗❌:             atom->getQuery()->beginChildren(), atom->getQuery()->endChildren());
    // RDKit❗❌:         while (childStack.size()) {
    // RDKit❗❌:           QueryAtom::QUERYATOM_QUERY::CHILD_TYPE qry = childStack.front();
    // RDKit❗❌:           childStack.pop_front();
    // RDKit❗❌:           if (qry->getDescription() == "RecursiveStructure") {
    // RDKit❗❌:             auto *rsq = dynamic_cast<RecursiveStructureQuery *>(qry.get());
    // RDKit❗❌:             CHECK_INVARIANT(rsq, "could not convert recursive structure query");
    // RDKit❗❌:             RWMol *rqm = new RWMol(*rsq->getQueryMol());
    // RDKit❗❌:             mergeQueryHs(*rqm, mergeUnmappedOnly, mergeIsotopes);
    // RDKit❗❌:             rsq->setQueryMol(rqm);
    // RDKit❗❌:           } else if (qry->beginChildren() != qry->endChildren()) {
    // RDKit❗❌:             childStack.insert(childStack.end(), qry->beginChildren(),
    // RDKit❗❌:                               qry->endChildren());
    // RDKit❗❌:           }
    // RDKit❗❌:         }
    // RDKit❗❌:       }  // end of recursion loop
    // Native root dispatch precedes its FIFO childStack. A canonical Not
    // represents a source negation bit, not an extra native child edge.
    // Clone the actual owned recursive graph just as new RWMol(*getQueryMol)
    // does; a missing graph is an invariant error, never silent omission.
    // Queue visits each source node once with O(width) live references; no
    // query-tree cloning or guessed textual reconstruction is performed.
    let mut pending = VecDeque::from([query]);
    while let Some(mut node) = pending.pop_front() {
        while let QueryNode::Not(child) = node {
            node = child.as_mut();
        }
        match node {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                let mut nested = recursive
                    .query_graph()
                    .cloned()
                    .ok_or(SmartsParseError::MissingRecursiveQueryGraph)?;
                merge_query_hs_in_place(&mut nested, merge_unmapped_only, merge_isotopes)?;
                recursive.set_query_graph(nested);
            }
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                pending.extend(children.iter_mut());
            }
            QueryNode::Predicate(_) => {}
            QueryNode::Not(_) => unreachable!("source negation layers consumed above"),
        }
    }
    Ok(())
}

pub(crate) fn remap_query_substance_groups_after_removal(
    groups: &[cosmolkit_model::SubstanceGroup],
    mapping: &TopologyMapping,
) -> Result<Vec<cosmolkit_model::SubstanceGroup>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION removedParentInHierarchy
    // RDKit❗❌: bool removedParentInHierarchy(
    // RDKit❗❌:     unsigned int idx, const std::vector<SubstanceGroup> &sgs,
    // RDKit❗❌:     const boost::dynamic_bitset<> &toRemove,
    // RDKit❗❌:     const std::map<unsigned int, unsigned int> &indexLookup) {
    // RDKit❗❌:   PRECONDITION(idx < sgs.size(), "cannot find SubstanceGroup");
    // RDKit❗❌:   if (toRemove[idx]) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   unsigned int parent;
    // RDKit❗❌:   if (sgs[idx].getPropIfPresent("PARENT", parent)) {
    // RDKit❗❌:     auto piter = indexLookup.find(parent);
    // RDKit❗❌:     if (piter != indexLookup.end()) {
    // RDKit❗❌:       return removedParentInHierarchy(piter->second, sgs, toRemove,
    // RDKit❗❌:                                       indexLookup);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION removedParentInHierarchy
    // Source hierarchy closure is equivalent here because QueryGraph validation
    // requires dense group IDs and in-range parents. The queue avoids recursive
    // stack growth while preserving the source's ancestor-removal result.
    let atom_map = mapping.atoms().old_to_new();
    let bond_map = mapping.bonds().old_to_new();
    let mut removed_groups = vec![false; groups.len()];

    // BEGIN RDKIT CPP FUNCTION SubstanceGroup::includesAtom and includesBond
    // RDKit❗❌: bool SubstanceGroup::includesAtom(unsigned int atomIdx) const {
    // RDKit❗❌:   if (std::find(d_atoms.begin(), d_atoms.end(), atomIdx) != d_atoms.end()) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   if (std::find(d_patoms.begin(), d_patoms.end(), atomIdx) != d_patoms.end()) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto &ap : d_saps) {
    // RDKit❗❌:     if (ap.aIdx == atomIdx || ap.lvIdx == rdcast<int>(atomIdx)) {
    // RDKit❗❌:       return true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   return false;
    // RDKit❗❌: }
    // RDKit❗❌: bool SubstanceGroup::includesBond(unsigned int bondIdx) const {
    // RDKit❗❌:   if (std::find(d_bonds.begin(), d_bonds.end(), bondIdx) != d_bonds.end()) {
    // RDKit❗❌:     return true;
    // RDKit❗❌:   }
    // RDKit❗❌:   for (const auto &cs : d_cstates) {
    // RDKit❗❌:     if (cs.bondIdx == bondIdx) {
    // RDKit❗❌:       return true;
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❗:   return false;
    // RDKit❗❗: }
    // END RDKIT CPP FUNCTION SubstanceGroup::includesAtom and includesBond
    for (index, mapped) in atom_map.iter().enumerate() {
        if mapped.is_some() {
            continue;
        }
        let atom = AtomId::new(index);
        for group in groups {
            if group.includes_atom(atom) {
                removed_groups[group.id().index()] = true;
            }
        }
    }
    for (index, mapped) in bond_map.iter().enumerate() {
        if mapped.is_some() {
            continue;
        }
        let bond = BondId::new(index);
        for group in groups {
            // The pinned includesBond helper checks only d_bonds and cstates;
            // model includes_bond also sees XBHEAD/XBCORR and is wider here.
            if group.bonds().contains(&bond)
                || group.cstates().iter().any(|cstate| cstate.bond == bond)
            {
                removed_groups[group.id().index()] = true;
            }
        }
    }

    // BEGIN RDKIT CPP FUNCTION removeSubstanceGroupsReferencing parent propagation
    // RDKit❗❌: // now go through and keep everything that shouldn't be removed
    // RDKit❗❌: // and who doesn't have a PARENT that should be removed in their hierarchy
    // RDKit❗❌: if (piter != indexLookup.end() && !toRemove[piter->second]) {
    // RDKit❗❌:   if (!removedParentInHierarchy(piter->second, sgs, toRemove,
    // RDKit❗❌:                               indexLookup)) {
    // RDKit❗❌:     keepIt = true;
    // RDKit❗❌:   }
    // RDKit❗❌: }
    // END RDKIT CPP FUNCTION removeSubstanceGroupsReferencing parent propagation
    let mut children = vec![Vec::new(); groups.len()];
    for group in groups {
        if let Some(parent) = group.parent() {
            children[parent.index()].push(group.id().index());
        }
    }
    let mut pending = removed_groups
        .iter()
        .enumerate()
        .filter_map(|(index, is_removed)| is_removed.then_some(index))
        .collect::<VecDeque<_>>();
    while let Some(parent) = pending.pop_front() {
        for &child in &children[parent] {
            if !removed_groups[child] {
                removed_groups[child] = true;
                pending.push_back(child);
            }
        }
    }

    let mut group_map = vec![None; groups.len()];
    let mut next_id = 0;
    for group in groups {
        let index = group.id().index();
        if !removed_groups[index] {
            group_map[index] = Some(SubstanceGroupId::new(next_id));
            next_id += 1;
        }
    }

    let mut remapped = Vec::with_capacity(next_id);
    for group in groups {
        let index = group.id().index();
        let Some(new_id) = group_map[index] else {
            continue;
        };
        let Some(group) = group.remapped(new_id, atom_map, bond_map, &group_map) else {
            return Err(SmartsParseError::Parse(format!(
                "query SGroup {} retains a reference to a removed graph row",
                group.id().index()
            )));
        };
        remapped.push(group);
    }
    Ok(remapped)
}

pub(crate) fn remap_query_stereo_groups_after_removal(
    groups: &[StereoGroup],
    mapping: &TopologyMapping,
) -> Result<Vec<StereoGroup>, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION removeAtomFromGroups/removeBondFromGroups
    // RDKit❗❌: auto atomPos = findAtom(group);
    // RDKit❗❌: if (atomPos != group.d_atoms.end()) {
    // RDKit❗❌:   group.d_atoms.erase(atomPos);
    // RDKit❗❌: }
    // RDKit❗❌: auto bondPos = findBond(group);
    // RDKit❗❌: if (bondPos != group.d_bonds.end()) {
    // RDKit❗❌:   group.d_bonds.erase(bondPos);
    // RDKit❗❌: }
    // RDKit❗❌: groups.erase(std::remove_if(groups.begin(), groups.end(),
    // RDKit❗❌:                             [](const auto &gp) {
    // RDKit❗❌:                               return gp.getAtoms().empty() &&
    // RDKit❗❌:                                      gp.getBonds().empty();
    // RDKit❗❌:                             }),
    // RDKit❗❌:              groups.end());
    // END RDKIT CPP FUNCTION removeAtomFromGroups/removeBondFromGroups
    let mut groups = groups.to_vec();
    for (index, mapped) in mapping.bonds().old_to_new().iter().enumerate().rev() {
        if mapped.is_none() {
            let removed = BondId::new(index);
            for group in &mut groups {
                group.remove_bond(removed);
            }
            groups.retain(|group| !group.is_empty());
        }
    }
    for (index, mapped) in mapping.atoms().old_to_new().iter().enumerate().rev() {
        if mapped.is_none() {
            let removed = AtomId::new(index);
            for group in &mut groups {
                group.remove_atom(removed);
            }
            groups.retain(|group| !group.is_empty());
        }
    }

    groups
        .into_iter()
        .map(|group| {
            let atoms = group
                .atoms()
                .iter()
                .map(|atom| {
                    mapping.atoms().old_to_new()[atom.index()].ok_or_else(|| {
                        SmartsParseError::Parse(format!(
                            "stereo group retains removed atom {} after source removal",
                            atom.index()
                        ))
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let bonds = group
                .bonds()
                .iter()
                .map(|bond| {
                    mapping.bonds().old_to_new()[bond.index()].ok_or_else(|| {
                        SmartsParseError::Parse(format!(
                            "stereo group retains removed bond {} after source removal",
                            bond.index()
                        ))
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mut remapped = StereoGroup::new(group.kind(), atoms, bonds);
            // RDKit❗✔️: extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
            // A removal retains the original group object and both IDs.
            cosmolkit_model::set_stereo_group_write_id(
                &mut remapped,
                cosmolkit_model::stereo_group_write_id(&group),
            );
            Ok(if let Some(id) = group.id() {
                remapped.with_id(id)
            } else {
                remapped
            })
        })
        .collect()
}

fn merge_query_hs_in_place(
    molecule: &mut QueryGraph,
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
    // plus O(R log R) removal sorting and one canonical QueryGraph rebuild,
    // adding a material O(V + E) clone versus
    // RDKit's batch-edit mutation; this accounts for the second-axis gap.
    let hydrogen_types = molecule
        .atoms()
        .iter()
        .enumerate()
        .map(|(index, atom)| {
            is_query_hydrogen(atom, molecule.adjacency().get(index).map_or(0, Vec::len))
        })
        .collect::<Vec<_>>();
    let mut removals = Vec::new();
    // RDKit uses an unsigned int for numHsToRemove; keep counts wider than a
    // byte so every accepted neighboring H contributes its source !H<i> term.
    let mut hydrogen_counts = vec![0_u32; molecule.num_atoms()];
    for atom_index in 0..molecule.num_atoms() {
        if hydrogen_types[atom_index] == QueryHydrogenType::QueryHydrogen {
            continue;
        }
        for &(neighbor_index, _) in &molecule.adjacency()[atom_index] {
            if hydrogen_types[neighbor_index] != QueryHydrogenType::QueryHydrogen {
                continue;
            }
            let hydrogen = &molecule.atoms()[neighbor_index];
            let map_ok = !merge_unmapped_only || hydrogen.atom_map().is_none();
            let isotope_ok =
                merge_isotopes || hydrogen.isotope().is_none_or(|isotope| isotope == 0);
            if map_ok && isotope_ok {
                removals.push(hydrogen.id());
                hydrogen_counts[atom_index] = hydrogen_counts[atom_index].wrapping_add(1);
            }
        }
    }
    removals.sort_unstable_by_key(|atom| atom.index());
    removals.dedup();
    let mut removed_atoms = vec![false; molecule.num_atoms()];
    for atom in &removals {
        removed_atoms[atom.index()] = true;
    }
    // First pass: build the complete old-to-new map for every surviving atom
    // before any predicate is touched or any carrier is remapped. A carrier
    // may reference a later surviving atom, so remapping against a partially
    // built map would falsely report that target as removed.
    let mut atom_old_to_new = vec![None; molecule.num_atoms()];
    let mut next_index = 0_usize;
    for (atom_index, mapping) in atom_old_to_new.iter_mut().enumerate() {
        if !removed_atoms[atom_index] {
            *mapping = Some(AtomId::new(next_index));
            next_index += 1;
        }
    }
    let mut atom_new_to_old = vec![None; next_index];
    for (old_index, new_id) in atom_old_to_new.iter().enumerate() {
        if let Some(new_id) = new_id {
            atom_new_to_old[new_id.index()] = Some(AtomId::new(old_index));
        }
    }

    let mut bond_old_to_new = vec![None; molecule.num_bonds()];
    let mut bond_new_to_old = Vec::new();
    for bond in molecule.bonds() {
        if atom_old_to_new[bond.begin().index()].is_none()
            || atom_old_to_new[bond.end().index()].is_none()
        {
            continue;
        }
        let new_id = BondId::new(bond_new_to_old.len());
        bond_old_to_new[bond.id().index()] = Some(new_id);
        bond_new_to_old.push(Some(bond.id()));
    }
    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    mapping
        .validate_for_counts(
            molecule.num_atoms(),
            next_index,
            molecule.num_bonds(),
            mapping.bonds().new_to_old().len(),
        )
        .map_err(|error| {
            SmartsParseError::Parse(format!("query-H removal produced invalid mapping: {error}"))
        })?;

    // Second pass: every predicate edit and remap is applied to a clone. The
    // caller's graph is assigned only after the whole rebuild succeeds, so a
    // failed remap leaves `molecule` unchanged.
    let mut atoms = Vec::with_capacity(next_index);
    for atom in molecule.atoms() {
        let Some(new_id) = mapping.atoms().old_to_new()[atom.index()] else {
            continue;
        };
        let mut atom_value = atom.clone().with_id(new_id);
        let mut predicate =
            std::mem::replace(atom_value.predicate_mut(), QueryNode::and(Vec::new()));
        let count = hydrogen_counts[atom.index()];
        for hydrogen_count in 0..count {
            crate::query_behavior::query_atom_expand_query(
                &mut predicate,
                QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::HydrogenCount(hydrogen_count as i32),
                ))),
                CompositeQueryType::And,
                true,
            );
        }
        if hydrogen_types[atom.index()] != QueryHydrogenType::QueryHydrogen {
            merge_recursive_query_hydrogens(&mut predicate, merge_unmapped_only, merge_isotopes)?;
        }
        // This compaction renumbers the query atom table, so surviving
        // carriers remap their template attachment targets through the one
        // shared primitive; a carrier whose referenced atom is a removed
        // hydrogen fails under the established lost-target policy instead of
        // keeping a stale row.
        atom_value
            .remap_template_attachment_order(mapping.atoms().old_to_new())
            .map_err(|source| SmartsParseError::TemplateAttachmentRemap {
                carrier: atom.index(),
                source,
            })?;
        *atom_value.predicate_mut() = predicate;
        atoms.push(atom_value);
    }
    let mut bonds = Vec::new();
    for bond in molecule.bonds() {
        let Some(new_id) = mapping.bonds().old_to_new()[bond.id().index()] else {
            continue;
        };
        let begin = mapping.atoms().old_to_new()[bond.begin().index()].ok_or_else(|| {
            SmartsParseError::Parse(format!(
                "retained query bond {} has removed begin atom",
                bond.id().index()
            ))
        })?;
        let end = mapping.atoms().old_to_new()[bond.end().index()].ok_or_else(|| {
            SmartsParseError::Parse(format!(
                "retained query bond {} has removed end atom",
                bond.id().index()
            ))
        })?;
        let stereo_atoms = bond.bond().stereo_atoms().and_then(|[first, second]| {
            Some([
                mapping.atoms().old_to_new()[first.index()]?,
                mapping.atoms().old_to_new()[second.index()]?,
            ])
        });
        bonds.push(QueryBond::from_parts(
            bond.bond()
                .clone()
                .remapped(new_id, begin, end, stereo_atoms),
            bond.predicate().clone(),
        ));
    }

    let substance_groups =
        remap_query_substance_groups_after_removal(query_substance_groups(molecule), &mapping)?;
    let stereo_groups =
        remap_query_stereo_groups_after_removal(molecule.stereo_groups(), &mapping)?;
    let props = molecule
        .ordered_props()
        .map(|(key, value)| (key.clone(), value.clone()));
    // RDKit source RWMol.cpp::batchRemoveAtoms preserves every conformer and
    // removes the coordinates whose old atoms were deleted:
    // RDKit❗❌: // do the same with the coordinates in the conformations
    // RDKit❗❌: for (auto conf : d_confs) {
    // RDKit❗❌:   RDGeom::POINT3D_VECT &positions = conf->getPositions();
    // RDKit❗❌:   RDGeom::POINT3D_VECT newPositions;
    // RDKit❗❌:   newPositions.reserve(getNumAtoms());
    // RDKit❗❌:   for (RDGeom::POINT3D_VECT::size_type i = 0; i < positions.size(); ++i) {
    // RDKit❗❌:     if (oldIndices[i] != nullptr) {
    // RDKit❗❌:       newPositions.push_back(positions[i]);
    // RDKit❗❌:     }
    // RDKit❗❌:   }
    // RDKit❗❌:   CHECK_INVARIANT(newPositions.size() == getNumAtoms(), "Lost coordinates!");
    // RDKit❗❌:   positions.swap(newPositions);
    // RDKit❗❌: }
    let coordinates = molecule.coordinate_block(None);
    let conformers_2d = coordinates
        .conformers_2d
        .iter()
        .map(|conformer| {
            let coordinates = mapping
                .atoms()
                .old_to_new()
                .iter()
                .enumerate()
                .filter_map(|(index, new_id)| new_id.map(|_| conformer.coordinates()[index]))
                .collect();
            let mut remapped = cosmolkit_model::Conformer2D::new(conformer.id(), coordinates);
            for (key, value) in conformer.props() {
                remapped = remapped.with_prop(key.clone(), value.clone());
            }
            remapped
        })
        .collect();
    let conformers_3d = coordinates
        .conformers_3d
        .iter()
        .map(|conformer| {
            let coordinates = mapping
                .atoms()
                .old_to_new()
                .iter()
                .enumerate()
                .filter_map(|(index, new_id)| new_id.map(|_| conformer.coordinates()[index]))
                .collect();
            let mut remapped =
                cosmolkit_model::Conformer3D::new(conformer.id(), coordinates, conformer.is_3d());
            for (key, value) in conformer.props() {
                remapped = remapped.with_prop(key.clone(), value.clone());
            }
            remapped
        })
        .collect();
    let mut rebuilt = QueryGraph::from_parts(
        atoms,
        bonds,
        props,
        conformers_2d,
        conformers_3d,
        stereo_groups,
    )
    .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    rebuilt
        .set_source_conformer_order(coordinates.source_conformer_order)
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    replace_query_substance_groups(&mut rebuilt, substance_groups)
        .map_err(|error| SmartsParseError::Parse(error.to_string()))?;
    *molecule = rebuilt;
    Ok(())
}

fn merge_query_hs(
    molecule: &QueryGraph,
    merge_unmapped_only: bool,
    merge_isotopes: bool,
) -> Result<QueryGraph, SmartsParseError> {
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
            .query_graph()
            .map(has_query_hs_graph)
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

fn has_query_hs(molecule: &QueryGraph) -> (bool, bool) {
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
    has_query_hs_graph(molecule)
}

fn has_query_hs_graph(graph: &QueryGraph) -> (bool, bool) {
    let mut query_hs = false;
    for (index, atom) in graph.atoms().iter().enumerate() {
        let degree = graph.adjacency().get(index).map_or(0, Vec::len);
        match is_query_hydrogen(atom, degree) {
            QueryHydrogenType::UnmergableQueryHydrogen => return (true, true),
            QueryHydrogenType::QueryHydrogen => query_hs = true,
            QueryHydrogenType::NotAHydrogen => {}
        }
        let result = query_node_h_status(atom.predicate());
        if result.1 {
            return result;
        }
        query_hs |= result.0;
    }
    (query_hs, false)
}

fn label_recursive_patterns(sma: impl AsRef<[u8]>) -> PropertyText {
    let sma = sma.as_ref();
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
    // Behavior: the pinned default automatic-relabeling branch uses exact
    // counted source bytes. Recursive keys include the original `$(` and `)`,
    // including nested source text rather than previously inserted labels.
    // Existing underscore suppresses lookup/insertion; repeated keys reuse
    // their first decimal label. A stray close discards all generated labels
    // and returns every original byte. Both source unsigned32 cursor increments
    // and size_t label addition wrap explicitly, without a size heuristic.
    // The independent NO_AUTOMATIC_SMARTS_RELABELLING build branch is outside
    // this pinned default configuration; no parser policy switch is invented.
    // Cost: one cursor pass and amortized output/stack append. Let L be the
    // sum of recursive substring lengths and p the number of distinct keys;
    // copied keys and ordered-map byte comparisons retain source O(L log p)
    // work, including nested long keys. Vec stacks avoid list nodes; no
    // UTF-8 validation, reparse, fallback decoder or additional input copy.
    #[derive(Clone, Copy, PartialEq)]
    enum SmaState {
        Base,
        Branch,
        Recurse,
    }
    use SmaState::*;

    let sma = sma.as_ref();
    let mut state: Vec<SmaState> = vec![Base];
    let mut start_recurse: Vec<u32> = Vec::new();
    let mut patterns: BTreeMap<Vec<u8>, String> = BTreeMap::new();
    let mut res = Vec::with_capacity(sma.len());
    let bytes = sma;

    let mut pos: u32 = 0;
    while (pos as usize) < bytes.len() {
        res.push(bytes[pos as usize]);
        if bytes[pos as usize] == b'$'
            && (pos.wrapping_add(1) as usize) < bytes.len()
            && bytes[pos.wrapping_add(1) as usize] == b'('
        {
            state.push(Recurse);
            start_recurse.push(pos);
            pos = pos.wrapping_add(1);
            res.push(bytes[pos as usize]);
        } else if bytes[pos as usize] == b'(' {
            state.push(Branch);
        } else if bytes[pos as usize] == b')' {
            if state.is_empty() || state.last() == Some(&Base) {
                return PropertyText::from_bytes(sma);
            }
            let curr_state = state.pop().expect("non-base SMARTS state");
            if curr_state == Recurse {
                let dollar_pos = start_recurse.pop().expect("recursive SMARTS start");
                if (pos.wrapping_add(1) as usize) >= bytes.len()
                    || bytes[pos.wrapping_add(1) as usize] != b'_'
                {
                    let length = pos.wrapping_sub(dollar_pos).wrapping_add(1) as usize;
                    let recurs = &bytes[dollar_pos as usize..dollar_pos as usize + length];
                    let label = if let Some(lbl) = patterns.get(recurs) {
                        lbl.clone()
                    } else {
                        let lbl = format!("{}", patterns.len().wrapping_add(100));
                        patterns.insert(recurs.to_vec(), lbl.clone());
                        lbl
                    };
                    res.push(b'_');
                    res.extend_from_slice(label.as_bytes());
                }
            }
        }
        pos = pos.wrapping_add(1);
    }
    PropertyText::from(res)
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
    /// Organic element symbol: B, C, N, O, S, P, F, Cl, Br, I
    OrganicElement(String),
    /// Aromatic element: c, n, o, s, p
    AromaticElement(String),
    /// Generic aromatic/aliphatic/wildcard atom query from SIMPLE_ATOM_QUERY_TOKEN.
    SimpleAtomQuery(u8),
    /// Flex's one-byte BAD_CHARACTER token, preserved until parser dispatch.
    BadCharacter(u8),
    /// Bracket atom text and the lexical atom tokens selected inside it.
    BracketContent(BracketContent),
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
    Symbol(u8),
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
    SimpleAtomQuery(u8),
    AtomElement(String),
    AtomPrimitive(u8),
    BondSpec(u8),
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
    BadCharacter(u8),
    Not,
    Semi,
    And,
    Or,
    EndOfStream,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ScannedToken {
    token: ScannerToken,
    span: SmartsTokenSpan,
}

#[derive(Debug, Clone, PartialEq)]
struct BracketContent {
    text: Vec<u8>,
    span: SmartsTokenSpan,
    lexical_tokens: Vec<ScannedToken>,
}

/// A scanner/parser span carries byte offsets in the original helper input
/// and the trimmed parser buffer. Bracket tokens use byte offsets relative
/// to the bracket-content origin; no text decoding or scalar indexing occurs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SmartsTokenSpan {
    input_char_start: usize,
    input_char_end: usize,
    parser_byte_start: usize,
    parser_byte_end: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SmartsScannerInputWindow {
    byte_start: usize,
    byte_end: usize,
}

fn setup_smarts_input(input: impl AsRef<[u8]>) -> SmartsScannerInputWindow {
    let input = input.as_ref();
    // BEGIN RDKIT CPP FUNCTION setup_smarts_string
    // RDKit✔️❌: size_t setup_smarts_string(const std::string &text,yyscan_t yyscanner){
    // RDKit✔️❌:   yyconst char * yybytes = text.c_str();
    // RDKit✔️❌:   yy_size_t _yybytes_len=text.size(), n, start, end;
    // RDKit✔️❌:   for(start = 0 ; start < _yybytes_len; ++start) {
    // RDKit✔️❌:     if (yybytes[start] > 32) { break; }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   for(end = _yybytes_len ; end > start; --end) {
    // RDKit✔️❌:     if (yybytes[end] > 32) { break; }
    // RDKit✔️❌:   }
    // RDKit✔️❌:   _yybytes_len = end-start+1;
    // RDKit✔️❌:   memcpy(buf, yybytes+start, _yybytes_len);
    // RDKit✔️❌:   buf[_yybytes_len] = buf[_yybytes_len+1] = YY_END_OF_BUFFER_CHAR;
    // RDKit✔️❌:   return start;
    // RDKit✔️❌: }
    // Local complexity review: setup keeps the same two O(n) byte passes and
    // uses a borrowed window; the canonical raw scanner in
    // `SmartsScanner::new` owns one O(n) raw byte buffer, matching Flex storage.
    let bytes = input;
    let mut start = 0;
    while start < bytes.len() {
        if i8::from_ne_bytes([bytes[start]]) > 32 {
            break;
        }
        start += 1;
    }

    let mut end = bytes.len();
    while end > start {
        // setup_smarts_string starts at text.size(), so this first read is
        // the NUL terminator supplied by std::string::c_str().
        let byte = bytes.get(end).copied().unwrap_or(0);
        if i8::from_ne_bytes([byte]) > 32 {
            break;
        }
        end -= 1;
    }

    // In the all-trimmed case, source copies only that terminal NUL. The
    // equivalent Rust scanner window is empty; otherwise end identifies the
    // last copied input byte and is inclusive in RDKit's length expression.
    let byte_end = end + 1;
    SmartsScannerInputWindow {
        byte_start: start,
        byte_end,
    }
}

struct SmartsScanner {
    bytes: Vec<u8>,
    parser_byte_base: usize,
    start: ScannerStart,
    states: Vec<ScannerState>,
    pos: usize,
    scan_end: usize,
}

impl SmartsScanner {
    fn new(input: &[u8], start: ScannerStart, window: SmartsScannerInputWindow) -> Self {
        // RDKit❗✔️:     ltrim = string_setup(inp, scanner);
        // RDKit❗✔️:     res = parser(inp.c_str() + ltrim, &molVect, atom, bond,
        // RDKit❗✔️:                          numAtomsParsed, numBondsParsed, branchPoints, scanner,
        // RDKit❗✔️:                          start_tok, current_token_position);
        // Every scanner position is a byte offset into the same counted input.
        // One O(n) byte copy matches Flex buffer ownership; no Unicode map.
        Self {
            bytes: input.to_vec(),
            parser_byte_base: window.byte_start,
            start,
            states: vec![ScannerState::Initial],
            pos: window.byte_start,
            scan_end: window.byte_end,
        }
    }

    fn state(&self) -> ScannerState {
        *self
            .states
            .last()
            .expect("scanner state stack is non-empty")
    }

    fn source_span(&self, input_char_start: usize, input_char_end: usize) -> SmartsTokenSpan {
        // RDKit✔️✔️: #define YY_USER_ACTION current_token_position += yyleng;
        // Existing span names are retained privately; both origins index bytes.
        debug_assert!(input_char_start <= input_char_end && input_char_end <= self.bytes.len());
        SmartsTokenSpan {
            input_char_start,
            input_char_end,
            parser_byte_start: input_char_start - self.parser_byte_base,
            parser_byte_end: input_char_end - self.parser_byte_base,
        }
    }

    fn emit(&mut self, token: ScannerToken, width: usize) -> ScannedToken {
        // RDKit❗✔️: #define YY_USER_ACTION current_token_position += yyleng;
        let start = self.pos;
        self.pos += width;
        ScannedToken {
            token,
            span: self.source_span(start, self.pos),
        }
    }

    fn bad_character_error_position(&self) -> usize {
        // RDKit❗✔️: #define YY_USER_ACTION current_token_position += yyleng;
        // RDKit❗✔️: .		return BAD_CHARACTER;
        // Flex's `.` consumes one source byte even when it begins a multibyte
        // UTF-8 scalar, so the parser counter is post-consumption by one byte.
        self.pos - self.parser_byte_base + 1
    }

    fn emit_bad_character(&mut self, character: u8) -> ScannedToken {
        // RDKit❗✔️: #define YY_USER_ACTION current_token_position += yyleng;
        // RDKit❗✔️: .		return BAD_CHARACTER;
        // Flex consumes exactly one byte. The input span and parser counter
        // advance by that byte, including NUL and invalid UTF-8. Scanning stops
        // at this token for the existing parser's error priority.
        let start = self.pos;
        let parser_byte_end = self.bad_character_error_position();
        self.pos += 1;
        let span = self.source_span(start, self.pos);
        ScannedToken {
            token: ScannerToken::BadCharacter(character),
            span: SmartsTokenSpan {
                parser_byte_end,
                ..span
            },
        }
    }

    fn scan(mut self) -> Result<Vec<ScannedToken>, SmartsParseError> {
        // BEGIN COMPLETE PINNED SMARTS FLEX RULES
        // RDKit❗✔️: %%
        // RDKit❗✔️:
        // RDKit❗✔️: %{
        // RDKit❗✔️:   if (start_token)
        // RDKit❗✔️:     {
        // RDKit❗✔️:       int t = start_token;
        // RDKit❗✔️:       start_token = 0;
        // RDKit❗✔️:       return t;
        // RDKit❗✔️:     }
        // RDKit❗✔️: %}
        // RDKit❗✔️:
        // RDKit❗✔️: @[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️: @[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: @		{ return AT_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>He |
        // RDKit❗✔️: <IN_ATOM_STATE>Li |
        // RDKit❗✔️: <IN_ATOM_STATE>Be |
        // RDKit❗✔️: <IN_ATOM_STATE>Ne |
        // RDKit❗✔️: <IN_ATOM_STATE>Na |
        // RDKit❗✔️: <IN_ATOM_STATE>Mg |
        // RDKit❗✔️: <IN_ATOM_STATE>Al |
        // RDKit❗✔️: <IN_ATOM_STATE>Si |
        // RDKit❗✔️: <IN_ATOM_STATE>Ar |
        // RDKit❗✔️: <IN_ATOM_STATE>K |
        // RDKit❗✔️: <IN_ATOM_STATE>Ca |
        // RDKit❗✔️: <IN_ATOM_STATE>Sc |
        // RDKit❗✔️: <IN_ATOM_STATE>Ti |
        // RDKit❗✔️: <IN_ATOM_STATE>V |
        // RDKit❗✔️: <IN_ATOM_STATE>Cr |
        // RDKit❗✔️: <IN_ATOM_STATE>Mn |
        // RDKit❗✔️: <IN_ATOM_STATE>Co |
        // RDKit❗✔️: <IN_ATOM_STATE>Fe |
        // RDKit❗✔️: <IN_ATOM_STATE>Ni |
        // RDKit❗✔️: <IN_ATOM_STATE>Cu |
        // RDKit❗✔️: <IN_ATOM_STATE>Zn |
        // RDKit❗✔️: <IN_ATOM_STATE>Ga |
        // RDKit❗✔️: <IN_ATOM_STATE>Ge |
        // RDKit❗✔️: <IN_ATOM_STATE>As |
        // RDKit❗✔️: <IN_ATOM_STATE>Se |
        // RDKit❗✔️: <IN_ATOM_STATE>Kr |
        // RDKit❗✔️: <IN_ATOM_STATE>Rb |
        // RDKit❗✔️: <IN_ATOM_STATE>Sr |
        // RDKit❗✔️: <IN_ATOM_STATE>Y |
        // RDKit❗✔️: <IN_ATOM_STATE>Zr |
        // RDKit❗✔️: <IN_ATOM_STATE>Nb |
        // RDKit❗✔️: <IN_ATOM_STATE>Mo |
        // RDKit❗✔️: <IN_ATOM_STATE>Tc |
        // RDKit❗✔️: <IN_ATOM_STATE>Ru |
        // RDKit❗✔️: <IN_ATOM_STATE>Rh |
        // RDKit❗✔️: <IN_ATOM_STATE>Pd |
        // RDKit❗✔️: <IN_ATOM_STATE>Ag |
        // RDKit❗✔️: <IN_ATOM_STATE>Cd |
        // RDKit❗✔️: <IN_ATOM_STATE>In |
        // RDKit❗✔️: <IN_ATOM_STATE>Sn |
        // RDKit❗✔️: <IN_ATOM_STATE>Sb |
        // RDKit❗✔️: <IN_ATOM_STATE>Te |
        // RDKit❗✔️: <IN_ATOM_STATE>Xe |
        // RDKit❗✔️: <IN_ATOM_STATE>Cs |
        // RDKit❗✔️: <IN_ATOM_STATE>Ba |
        // RDKit❗✔️: <IN_ATOM_STATE>La |
        // RDKit❗✔️: <IN_ATOM_STATE>Ce |
        // RDKit❗✔️: <IN_ATOM_STATE>Pr |
        // RDKit❗✔️: <IN_ATOM_STATE>Nd |
        // RDKit❗✔️: <IN_ATOM_STATE>Pm |
        // RDKit❗✔️: <IN_ATOM_STATE>Sm |
        // RDKit❗✔️: <IN_ATOM_STATE>Eu |
        // RDKit❗✔️: <IN_ATOM_STATE>Gd |
        // RDKit❗✔️: <IN_ATOM_STATE>Tb |
        // RDKit❗✔️: <IN_ATOM_STATE>Dy |
        // RDKit❗✔️: <IN_ATOM_STATE>Ho |
        // RDKit❗✔️: <IN_ATOM_STATE>Er |
        // RDKit❗✔️: <IN_ATOM_STATE>Tm |
        // RDKit❗✔️: <IN_ATOM_STATE>Yb |
        // RDKit❗✔️: <IN_ATOM_STATE>Lu |
        // RDKit❗✔️: <IN_ATOM_STATE>Hf |
        // RDKit❗✔️: <IN_ATOM_STATE>Ta |
        // RDKit❗✔️: <IN_ATOM_STATE>W |
        // RDKit❗✔️: <IN_ATOM_STATE>Re |
        // RDKit❗✔️: <IN_ATOM_STATE>Os |
        // RDKit❗✔️: <IN_ATOM_STATE>Ir |
        // RDKit❗✔️: <IN_ATOM_STATE>Pt |
        // RDKit❗✔️: <IN_ATOM_STATE>Au |
        // RDKit❗✔️: <IN_ATOM_STATE>Hg |
        // RDKit❗✔️: <IN_ATOM_STATE>Tl |
        // RDKit❗✔️: <IN_ATOM_STATE>Pb |
        // RDKit❗✔️: <IN_ATOM_STATE>Bi |
        // RDKit❗✔️: <IN_ATOM_STATE>Po |
        // RDKit❗✔️: <IN_ATOM_STATE>At |
        // RDKit❗✔️: <IN_ATOM_STATE>Rn |
        // RDKit❗✔️: <IN_ATOM_STATE>Fr |
        // RDKit❗✔️: <IN_ATOM_STATE>Ra |
        // RDKit❗✔️: <IN_ATOM_STATE>Ac |
        // RDKit❗✔️: <IN_ATOM_STATE>Th |
        // RDKit❗✔️: <IN_ATOM_STATE>Pa |
        // RDKit❗✔️: <IN_ATOM_STATE>U |
        // RDKit❗✔️: <IN_ATOM_STATE>Np |
        // RDKit❗✔️: <IN_ATOM_STATE>Pu |
        // RDKit❗✔️: <IN_ATOM_STATE>Am |
        // RDKit❗✔️: <IN_ATOM_STATE>Cm |
        // RDKit❗✔️: <IN_ATOM_STATE>Bk |
        // RDKit❗✔️: <IN_ATOM_STATE>Cf |
        // RDKit❗✔️: <IN_ATOM_STATE>Es |
        // RDKit❗✔️: <IN_ATOM_STATE>Fm |
        // RDKit❗✔️: <IN_ATOM_STATE>Md |
        // RDKit❗✔️: <IN_ATOM_STATE>No |
        // RDKit❗✔️: <IN_ATOM_STATE>Lr |
        // RDKit❗✔️: <IN_ATOM_STATE>Rf |
        // RDKit❗✔️: <IN_ATOM_STATE>Db |
        // RDKit❗✔️: <IN_ATOM_STATE>Sg |
        // RDKit❗✔️: <IN_ATOM_STATE>Bh |
        // RDKit❗✔️: <IN_ATOM_STATE>Hs |
        // RDKit❗✔️: <IN_ATOM_STATE>Mt |
        // RDKit❗✔️: <IN_ATOM_STATE>Ds |
        // RDKit❗✔️: <IN_ATOM_STATE>Rg |
        // RDKit❗✔️: <IN_ATOM_STATE>Cn |
        // RDKit❗✔️: <IN_ATOM_STATE>Uut |
        // RDKit❗✔️: <IN_ATOM_STATE>Fl |
        // RDKit❗✔️: <IN_ATOM_STATE>Uup |
        // RDKit❗✔️: <IN_ATOM_STATE>Lv	{   yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
        // RDKit❗✔️: 				return ATOM_TOKEN;
        // RDKit❗✔️: 			}
        // RDKit❗✔️: <IN_ATOM_STATE>D {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomExplicitDegreeQuery(1));
        // RDKit❗✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: <IN_ATOM_STATE>d {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomNonHydrogenDegreeQuery(1));
        // RDKit❗✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>X {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomTotalDegreeQuery(1));
        // RDKit❗✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>x {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHasRingBondQuery());
        // RDKit❗✔️: 	return RINGBOND_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>v {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomTotalValenceQuery(1));
        // RDKit❗✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>z {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHasHeteroatomNbrsQuery());
        // RDKit❗✔️: 	return HETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>Z {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHasAliphaticHeteroatomNbrsQuery());
        // RDKit❗✔️: 	return ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>h {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️:         yylval->atom->setQuery(makeAtomHasImplicitHQuery());
        // RDKit❗✔️: 	return IMPLICIT_H_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>R {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(new AtomRingQuery(-1));
        // RDKit❗✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>r {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit❗✔️: 	return MIN_RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>k {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit❗✔️: 	return RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: H			{  return H_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: B			{  yylval->ival = 5;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: C			{  yylval->ival = 6;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: N			{  yylval->ival = 7;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: O			{  yylval->ival = 8;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: F			{  yylval->ival = 9;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: P			{  yylval->ival = 15;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: S			{  yylval->ival = 16;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: Cl			{  yylval->ival = 17;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: Br			{  yylval->ival = 35;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: I			{  yylval->ival = 53;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: b			{  yylval->ival = 5;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: c			{  yylval->ival = 6;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: n			{  yylval->ival = 7;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: o			{  yylval->ival = 8;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: p			{  yylval->ival = 15;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: s			{  yylval->ival = 16;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>si	{  yylval->ival = 14;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>as	{  yylval->ival = 33;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>se	{  yylval->ival = 34;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>te	{  yylval->ival = 52;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: \*			{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomNullQuery());
        // RDKit❗✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: a			{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomAromaticQuery());
        // RDKit❗✔️: 	yylval->atom->setIsAromatic(true);
        // RDKit❗✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: A			{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomAliphaticQuery());
        // RDKit❗✔️: 	return SIMPLE_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: \: 			{ return COLON_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \_ 			{ return UNDERSCORE_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \#			{ return HASH_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \=	{ yylval->bond = new QueryBond(Bond::DOUBLE);
        // RDKit❗✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::DOUBLE));
        // RDKit❗✔️: 	return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: \~	{ yylval->bond = new QueryBond();
        // RDKit❗✔️: 	yylval->bond->setQuery(makeBondNullQuery());
        // RDKit❗✔️: 	return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: \$	{ yylval->bond = new QueryBond(Bond::QUADRUPLE);
        // RDKit❗✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::QUADRUPLE));
        // RDKit❗✔️:     return BOND_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: [\\]{1,2}    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
        // RDKit❗✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗✔️: 	return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: [\/]    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
        // RDKit❗✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗✔️: 	return BOND_TOKEN;  }
        // RDKit❗✔️:
        // RDKit❗✔️: \-\> {
        // RDKit❗✔️:     yylval->bond = new QueryBond(Bond::DATIVER);
        // RDKit❗✔️:     return BOND_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: \<\- {
        // RDKit❗✔️:     yylval->bond = new QueryBond(Bond::DATIVEL);
        // RDKit❗✔️:     return BOND_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: \-			{ return MINUS_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \+			{ return PLUS_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: <IN_ATOM_STATE>\$\(              { yy_push_state(IN_RECURSION_STATE,yyscanner); return BEGIN_RECURSE; }
        // RDKit❗✔️:
        // RDKit❗✔️: \(       	{ yy_push_state(IN_BRANCH_STATE,yyscanner); return GROUP_OPEN_TOKEN; }
        // RDKit❗✔️: <IN_BRANCH_STATE>\)       	{ yy_pop_state(yyscanner); return GROUP_CLOSE_TOKEN; }
        // RDKit❗✔️: <IN_RECURSION_STATE>\)       	{ yy_pop_state(yyscanner); return END_RECURSE; }
        // RDKit❗✔️:
        // RDKit❗✔️: \{       	{ return RANGE_OPEN_TOKEN; }
        // RDKit❗✔️: \}       	{ return RANGE_CLOSE_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️:
        // RDKit❗✔️: \[			{ yy_push_state(IN_ATOM_STATE,yyscanner); return ATOM_OPEN_TOKEN; }
        // RDKit❗✔️: <IN_ATOM_STATE>\]	{ yy_pop_state(yyscanner); return ATOM_CLOSE_TOKEN; }
        // RDKit❗✔️: \]			{ /* FIX: ???
        // RDKit❗✔️:                            This rule is here because otherwise recursive SMARTS queries like:
        // RDKit❗✔️: 	                   [$(C(=O)[O,N])] lex improperly (no ATOM_CLOSE token is returned).
        // RDKit❗✔️:  			   I am not 100% sure that the approach we're using here will work
        // RDKit❗✔️:                            all the time, but I'm hoping that any problems caused here in
        // RDKit❗✔️:                            the lexer will get caught in the parser.
        // RDKit❗✔️: 			  */
        // RDKit❗✔️:                           return ATOM_CLOSE_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \.       	{ return SEPARATOR_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \%              { return PERCENT_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: [0]		{ yylval->ival = 0;  return ZERO_TOKEN; }
        // RDKit❗✔️: [1-9]		{ yylval->ival = yytext[0]-'0';  return NONZERO_DIGIT_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \!			{ return NOT_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \;			{ return SEMI_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \&			{ return AND_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \,			{ return OR_TOKEN; }
        // RDKit❗✔️:
        // RDKit❗✔️: \^0		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::S));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: \^1		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: \^2		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP2));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️:
        // RDKit❗✔️: \^3		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: \^4		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: \^5		{
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D2));
        // RDKit❗✔️: 	return HYB_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: \n		return EOS_TOKEN;
        // RDKit❗✔️:
        // RDKit❗✔️: <<EOF>>		{ return EOS_TOKEN; }
        // RDKit❗✔️: .		return BAD_CHARACTER;
        // RDKit❗✔️:
        // RDKit❗✔️: %%
        // END COMPLETE PINNED SMARTS FLEX RULES

        let mut tokens = vec![ScannedToken {
            token: ScannerToken::Start(self.start),
            span: self.source_span(self.pos, self.pos),
        }];

        while self.pos < self.scan_end {
            let state = self.state();
            let ch = self.bytes[self.pos];

            // RDKit❗✔️: \n		return EOS_TOKEN;
            if ch == b'\n' {
                tokens.push(self.emit(ScannerToken::EndOfStream, 1));
                return Ok(tokens);
            }

            // RDKit✔️✔️: @[' ']*TH { yylval->chiraltype = Atom::ChiralType::CHI_TETRAHEDRAL; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*AL { yylval->chiraltype = Atom::ChiralType::CHI_ALLENE; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*SP { yylval->chiraltype = Atom::ChiralType::CHI_SQUAREPLANAR; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*TB { yylval->chiraltype = Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL; return CHI_CLASS_TOKEN; }
            // RDKit✔️✔️: @[' ']*OH { yylval->chiraltype = Atom::ChiralType::CHI_OCTAHEDRAL; return CHI_CLASS_TOKEN; }
            if ch == b'@' {
                let mut cursor = self.pos + 1;
                while self.bytes.get(cursor) == Some(&b' ') {
                    cursor += 1;
                }
                if cursor + 1 < self.scan_end {
                    let class = match &self.bytes[cursor..cursor + 2] {
                        b"TH" => "TH",
                        b"AL" => "AL",
                        b"SP" => "SP",
                        b"TB" => "TB",
                        b"OH" => "OH",
                        _ => "",
                    };
                    if matches!(class, "TH" | "AL" | "SP" | "TB" | "OH") {
                        tokens.push(self.emit(
                            ScannerToken::ChiralClass(class.to_owned()),
                            cursor + 2 - self.pos,
                        ));
                        continue;
                    }
                }
                // RDKit✔️✔️: @		{ return AT_TOKEN; }
                tokens.push(self.emit(ScannerToken::At, 1));
                continue;
            }

            // RDKit✔️✔️: <IN_ATOM_STATE>\$\(              { yy_push_state(IN_RECURSION_STATE,yyscanner); return BEGIN_RECURSE; }
            if state == ScannerState::Atom
                && ch == b'$'
                && self.bytes.get(self.pos + 1) == Some(&b'(')
            {
                self.states.push(ScannerState::Recursion);
                tokens.push(self.emit(ScannerToken::BeginRecurse, 2));
                continue;
            }

            // RDKit✔️✔️: \(       	{ yy_push_state(IN_BRANCH_STATE,yyscanner); return GROUP_OPEN_TOKEN; }
            // RDKit✔️✔️: <IN_BRANCH_STATE>\)       	{ yy_pop_state(yyscanner); return GROUP_CLOSE_TOKEN; }
            // RDKit✔️✔️: <IN_RECURSION_STATE>\)       	{ yy_pop_state(yyscanner); return END_RECURSE; }
            if ch == b'(' {
                self.states.push(ScannerState::Branch);
                tokens.push(self.emit(ScannerToken::GroupOpen, 1));
                continue;
            }
            if ch == b')' {
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

            // RDKit✔️✔️: \[			{ yy_push_state(IN_ATOM_STATE,yyscanner); return ATOM_OPEN_TOKEN; }
            // RDKit✔️✔️: <IN_ATOM_STATE>\]	{ yy_pop_state(yyscanner); return ATOM_CLOSE_TOKEN; }
            // RDKit✔️✔️: \]			{ /* FIX: ???
            // RDKit✔️✔️:                            This rule is here because otherwise recursive SMARTS queries like:
            // RDKit✔️✔️: 	                   [$(C(=O)[O,N])] lex improperly (no ATOM_CLOSE token is returned).
            // RDKit✔️✔️:  			   I am not 100% sure that the approach we're using here will work
            // RDKit✔️✔️:                            all the time, but I'm hoping that any problems caused here in
            // RDKit✔️✔️:                            the lexer will get caught in the parser.
            // RDKit✔️✔️: 			  */
            // RDKit✔️✔️:                           return ATOM_CLOSE_TOKEN; }
            if ch == b'[' {
                self.states.push(ScannerState::Atom);
                tokens.push(self.emit(ScannerToken::AtomOpen, 1));
                continue;
            }
            if ch == b']' {
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

            tokens.push(self.emit_bad_character(ch));
            return Ok(tokens);
        }

        // RDKit❗✔️: <<EOF>>		{ return EOS_TOKEN; }
        // The legacy CK variant points at the opening bracket; retain that
        // projection while expressing its position in parser-buffer bytes.
        if self.states.contains(&ScannerState::Atom) {
            let start = tokens
                .iter()
                .rev()
                .find(|token| token.token == ScannerToken::AtomOpen)
                .map_or(0, |token| token.span.parser_byte_start);
            return Err(SmartsParseError::UnclosedBracket(start));
        }

        // RDKit❗✔️: <<EOF>>		{ return EOS_TOKEN; }
        tokens.push(ScannedToken {
            token: ScannerToken::EndOfStream,
            span: self.source_span(self.pos, self.pos),
        });
        Ok(tokens)
    }

    fn scan_atom_token(&mut self) -> Result<Option<ScannedToken>, SmartsParseError> {
        let ch = self.bytes[self.pos];
        // Borrow the source byte suffix; literal matching has no suffix copy,
        // text decoding, scalar index or token-attempt allocation.
        let rest = &self.bytes[self.pos..self.scan_end];

        // RDKit✔️❌: <IN_ATOM_STATE>He |
        // RDKit✔️❌: <IN_ATOM_STATE>Li |
        // RDKit✔️❌: <IN_ATOM_STATE>Be |
        // RDKit✔️❌: <IN_ATOM_STATE>Ne |
        // RDKit✔️❌: <IN_ATOM_STATE>Na |
        // RDKit✔️❌: <IN_ATOM_STATE>Mg |
        // RDKit✔️❌: <IN_ATOM_STATE>Al |
        // RDKit✔️❌: <IN_ATOM_STATE>Si |
        // RDKit✔️❌: <IN_ATOM_STATE>Ar |
        // RDKit✔️❌: <IN_ATOM_STATE>K |
        // RDKit✔️❌: <IN_ATOM_STATE>Ca |
        // RDKit✔️❌: <IN_ATOM_STATE>Sc |
        // RDKit✔️❌: <IN_ATOM_STATE>Ti |
        // RDKit✔️❌: <IN_ATOM_STATE>V |
        // RDKit✔️❌: <IN_ATOM_STATE>Cr |
        // RDKit✔️❌: <IN_ATOM_STATE>Mn |
        // RDKit✔️❌: <IN_ATOM_STATE>Co |
        // RDKit✔️❌: <IN_ATOM_STATE>Fe |
        // RDKit✔️❌: <IN_ATOM_STATE>Ni |
        // RDKit✔️❌: <IN_ATOM_STATE>Cu |
        // RDKit✔️❌: <IN_ATOM_STATE>Zn |
        // RDKit✔️❌: <IN_ATOM_STATE>Ga |
        // RDKit✔️❌: <IN_ATOM_STATE>Ge |
        // RDKit✔️❌: <IN_ATOM_STATE>As |
        // RDKit✔️❌: <IN_ATOM_STATE>Se |
        // RDKit✔️❌: <IN_ATOM_STATE>Kr |
        // RDKit✔️❌: <IN_ATOM_STATE>Rb |
        // RDKit✔️❌: <IN_ATOM_STATE>Sr |
        // RDKit✔️❌: <IN_ATOM_STATE>Y |
        // RDKit✔️❌: <IN_ATOM_STATE>Zr |
        // RDKit✔️❌: <IN_ATOM_STATE>Nb |
        // RDKit✔️❌: <IN_ATOM_STATE>Mo |
        // RDKit✔️❌: <IN_ATOM_STATE>Tc |
        // RDKit✔️❌: <IN_ATOM_STATE>Ru |
        // RDKit✔️❌: <IN_ATOM_STATE>Rh |
        // RDKit✔️❌: <IN_ATOM_STATE>Pd |
        // RDKit✔️❌: <IN_ATOM_STATE>Ag |
        // RDKit✔️❌: <IN_ATOM_STATE>Cd |
        // RDKit✔️❌: <IN_ATOM_STATE>In |
        // RDKit✔️❌: <IN_ATOM_STATE>Sn |
        // RDKit✔️❌: <IN_ATOM_STATE>Sb |
        // RDKit✔️❌: <IN_ATOM_STATE>Te |
        // RDKit✔️❌: <IN_ATOM_STATE>Xe |
        // RDKit✔️❌: <IN_ATOM_STATE>Cs |
        // RDKit✔️❌: <IN_ATOM_STATE>Ba |
        // RDKit✔️❌: <IN_ATOM_STATE>La |
        // RDKit✔️❌: <IN_ATOM_STATE>Ce |
        // RDKit✔️❌: <IN_ATOM_STATE>Pr |
        // RDKit✔️❌: <IN_ATOM_STATE>Nd |
        // RDKit✔️❌: <IN_ATOM_STATE>Pm |
        // RDKit✔️❌: <IN_ATOM_STATE>Sm |
        // RDKit✔️❌: <IN_ATOM_STATE>Eu |
        // RDKit✔️❌: <IN_ATOM_STATE>Gd |
        // RDKit✔️❌: <IN_ATOM_STATE>Tb |
        // RDKit✔️❌: <IN_ATOM_STATE>Dy |
        // RDKit✔️❌: <IN_ATOM_STATE>Ho |
        // RDKit✔️❌: <IN_ATOM_STATE>Er |
        // RDKit✔️❌: <IN_ATOM_STATE>Tm |
        // RDKit✔️❌: <IN_ATOM_STATE>Yb |
        // RDKit✔️❌: <IN_ATOM_STATE>Lu |
        // RDKit✔️❌: <IN_ATOM_STATE>Hf |
        // RDKit✔️❌: <IN_ATOM_STATE>Ta |
        // RDKit✔️❌: <IN_ATOM_STATE>W |
        // RDKit✔️❌: <IN_ATOM_STATE>Re |
        // RDKit✔️❌: <IN_ATOM_STATE>Os |
        // RDKit✔️❌: <IN_ATOM_STATE>Ir |
        // RDKit✔️❌: <IN_ATOM_STATE>Pt |
        // RDKit✔️❌: <IN_ATOM_STATE>Au |
        // RDKit✔️❌: <IN_ATOM_STATE>Hg |
        // RDKit✔️❌: <IN_ATOM_STATE>Tl |
        // RDKit✔️❌: <IN_ATOM_STATE>Pb |
        // RDKit✔️❌: <IN_ATOM_STATE>Bi |
        // RDKit✔️❌: <IN_ATOM_STATE>Po |
        // RDKit✔️❌: <IN_ATOM_STATE>At |
        // RDKit✔️❌: <IN_ATOM_STATE>Rn |
        // RDKit✔️❌: <IN_ATOM_STATE>Fr |
        // RDKit✔️❌: <IN_ATOM_STATE>Ra |
        // RDKit✔️❌: <IN_ATOM_STATE>Ac |
        // RDKit✔️❌: <IN_ATOM_STATE>Th |
        // RDKit✔️❌: <IN_ATOM_STATE>Pa |
        // RDKit✔️❌: <IN_ATOM_STATE>U |
        // RDKit✔️❌: <IN_ATOM_STATE>Np |
        // RDKit✔️❌: <IN_ATOM_STATE>Pu |
        // RDKit✔️❌: <IN_ATOM_STATE>Am |
        // RDKit✔️❌: <IN_ATOM_STATE>Cm |
        // RDKit✔️❌: <IN_ATOM_STATE>Bk |
        // RDKit✔️❌: <IN_ATOM_STATE>Cf |
        // RDKit✔️❌: <IN_ATOM_STATE>Es |
        // RDKit✔️❌: <IN_ATOM_STATE>Fm |
        // RDKit✔️❌: <IN_ATOM_STATE>Md |
        // RDKit✔️❌: <IN_ATOM_STATE>No |
        // RDKit✔️❌: <IN_ATOM_STATE>Lr |
        // RDKit✔️❌: <IN_ATOM_STATE>Rf |
        // RDKit✔️❌: <IN_ATOM_STATE>Db |
        // RDKit✔️❌: <IN_ATOM_STATE>Sg |
        // RDKit✔️❌: <IN_ATOM_STATE>Bh |
        // RDKit✔️❌: <IN_ATOM_STATE>Hs |
        // RDKit✔️❌: <IN_ATOM_STATE>Mt |
        // RDKit✔️❌: <IN_ATOM_STATE>Ds |
        // RDKit✔️❌: <IN_ATOM_STATE>Rg |
        // RDKit✔️❌: <IN_ATOM_STATE>Cn |
        // RDKit✔️❌: <IN_ATOM_STATE>Uut |
        // RDKit✔️❌: <IN_ATOM_STATE>Fl |
        // RDKit✔️❌: <IN_ATOM_STATE>Uup |
        // RDKit✔️❌: <IN_ATOM_STATE>Lv	{   yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
        // RDKit✔️❌: 				return ATOM_TOKEN;
        // RDKit✔️❌: 			}
        // RDKit✔️❌: <IN_ATOM_STATE>D {
        // Flex uses longest-match selection, so three-letter temporary element
        // names and then two-letter names are checked before one-letter names.
        for symbol in ELEMENT_SYMBOLS {
            if rest.starts_with(symbol.as_bytes()) {
                return Ok(Some(self.emit(
                    ScannerToken::AtomElement((*symbol).to_string()),
                    symbol.len(),
                )));
            }
        }

        // RDKit✔️❌: <IN_ATOM_STATE>si	{  yylval->ival = 14;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: <IN_ATOM_STATE>as	{  yylval->ival = 33;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: <IN_ATOM_STATE>se	{  yylval->ival = 34;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: <IN_ATOM_STATE>te	{  yylval->ival = 52;  return AROMATIC_ATOM_TOKEN;  }
        for symbol in ["si", "as", "se", "te"] {
            if rest.starts_with(symbol.as_bytes()) {
                return Ok(Some(
                    self.emit(ScannerToken::AromaticElement(symbol.to_string()), 2),
                ));
            }
        }

        // RDKit✔️❌: <IN_ATOM_STATE>D {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomExplicitDegreeQuery(1));
        // RDKit✔️❌: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>d {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomNonHydrogenDegreeQuery(1));
        // RDKit✔️❌: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>X {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomTotalDegreeQuery(1));
        // RDKit✔️❌: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>x {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHasRingBondQuery());
        // RDKit✔️❌: 	return RINGBOND_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>v {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomTotalValenceQuery(1));
        // RDKit✔️❌: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>z {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHasHeteroatomNbrsQuery());
        // RDKit✔️❌: 	return HETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>Z {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHasAliphaticHeteroatomNbrsQuery());
        // RDKit✔️❌: 	return ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>h {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌:         yylval->atom->setQuery(makeAtomHasImplicitHQuery());
        // RDKit✔️❌: 	return IMPLICIT_H_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>R {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(new AtomRingQuery(-1));
        // RDKit✔️❌: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>r {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit✔️❌: 	return MIN_RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: <IN_ATOM_STATE>k {
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomInRingQuery());
        // RDKit✔️❌: 	return RINGSIZE_ATOM_QUERY_TOKEN;
        // RDKit✔️❌: }
        if matches!(
            ch,
            b'D' | b'd' | b'X' | b'x' | b'v' | b'z' | b'Z' | b'h' | b'R' | b'r' | b'k'
        ) {
            return Ok(Some(self.emit(ScannerToken::AtomPrimitive(ch), 1)));
        }

        // RDKit✔️❌: \^0		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::S));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: \^1		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: \^2		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP2));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: \^3		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: \^4		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        // RDKit✔️❌: \^5		{
        // RDKit✔️❌: 	yylval->atom = new QueryAtom();
        // RDKit✔️❌: 	yylval->atom->setQuery(makeAtomHybridizationQuery(Atom::SP3D2));
        // RDKit✔️❌: 	return HYB_TOKEN;
        // RDKit✔️❌: }
        if ch == b'^' {
            if let Some(value) = self
                .bytes
                .get(self.pos + 1)
                .filter(|digit| digit.is_ascii_digit())
                .map(|digit| u32::from(*digit - b'0'))
                .filter(|value| *value <= 5)
            {
                return Ok(Some(self.emit(ScannerToken::Hybridization(value as u8), 2)));
            }
        }
        Ok(None)
    }

    fn scan_common_token(&mut self) -> Result<Option<ScannedToken>, SmartsParseError> {
        let ch = self.bytes[self.pos];
        // Borrow the source byte suffix; literal matching has no suffix copy,
        // text decoding, scalar index or token-attempt allocation.
        let rest = &self.bytes[self.pos..self.scan_end];

        // RDKit✔️❌: B			{  yylval->ival = 5;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: C			{  yylval->ival = 6;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: N			{  yylval->ival = 7;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: O			{  yylval->ival = 8;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: F			{  yylval->ival = 9;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: P			{  yylval->ival = 15;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: S			{  yylval->ival = 16;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: Cl			{  yylval->ival = 17;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: Br			{  yylval->ival = 35;  return ORGANIC_ATOM_TOKEN;  }
        // RDKit✔️❌: I			{  yylval->ival = 53;  return ORGANIC_ATOM_TOKEN;  }
        for symbol in ["Cl", "Br", "B", "C", "N", "O", "F", "P", "S", "I"] {
            if rest.starts_with(symbol.as_bytes()) {
                return Ok(Some(self.emit(
                    ScannerToken::OrganicElement(symbol.to_string()),
                    symbol.len(),
                )));
            }
        }
        // RDKit✔️❌: b			{  yylval->ival = 5;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: c			{  yylval->ival = 6;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: n			{  yylval->ival = 7;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: o			{  yylval->ival = 8;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: p			{  yylval->ival = 15;  return AROMATIC_ATOM_TOKEN;  }
        // RDKit✔️❌: s			{  yylval->ival = 16;  return AROMATIC_ATOM_TOKEN;  }
        if matches!(ch, b'b' | b'c' | b'n' | b'o' | b'p' | b's') {
            return Ok(Some(
                self.emit(
                    ScannerToken::AromaticElement(
                        match ch {
                            b'b' => "b",
                            b'c' => "c",
                            b'n' => "n",
                            b'o' => "o",
                            b'p' => "p",
                            b's' => "s",
                            _ => unreachable!("matched aromatic lexer literal"),
                        }
                        .to_owned(),
                    ),
                    1,
                ),
            ));
        }
        // Preserve the simple-query token kind; parse_simple_atom constructs its value.
        if ch == b'*' || ch == b'A' || ch == b'a' {
            return Ok(Some(self.emit(ScannerToken::SimpleAtomQuery(ch), 1)));
        }
        // RDKit✔️❌: H			{  return H_TOKEN;  }
        if ch == b'H' {
            return Ok(Some(self.emit(ScannerToken::AtomPrimitive(b'H'), 1)));
        }

        // RDKit✔️❌: \: 			{ return COLON_TOKEN; }
        // RDKit✔️❌: \_ 			{ return UNDERSCORE_TOKEN; }
        // RDKit✔️❌: \#			{ return HASH_TOKEN; }
        // Bond query and carrier construction is anchored in parsed_bond_spec;
        // direction transport is handled by current_bond_direction.
        if rest.starts_with(b"->") {
            return Ok(Some(self.emit(ScannerToken::DativeRight, 2)));
        }
        if rest.starts_with(b"<-") {
            return Ok(Some(self.emit(ScannerToken::DativeLeft, 2)));
        }
        if matches!(ch, b'=' | b'~' | b'$' | b'/' | b'\\') {
            let width = if ch == b'\\' && self.bytes.get(self.pos + 1) == Some(&b'\\') {
                2
            } else {
                1
            };
            return Ok(Some(self.emit(ScannerToken::BondSpec(ch), width)));
        }

        // RDKit✔️❌: \: 			{ return COLON_TOKEN; }
        // RDKit✔️❌: \_ 			{ return UNDERSCORE_TOKEN; }
        // RDKit✔️❌: \#			{ return HASH_TOKEN; }
        // RDKit✔️❌: \-			{ return MINUS_TOKEN; }
        // RDKit✔️❌: \+			{ return PLUS_TOKEN; }
        // RDKit✔️❌: \{       	{ return RANGE_OPEN_TOKEN; }
        // RDKit✔️❌: \}       	{ return RANGE_CLOSE_TOKEN; }
        // RDKit✔️❌: \.       	{ return SEPARATOR_TOKEN; }
        // RDKit✔️❌: \%              { return PERCENT_TOKEN; }
        // RDKit✔️❌: [0]		{ yylval->ival = 0;  return ZERO_TOKEN; }
        // RDKit✔️❌: [1-9]		{ yylval->ival = yytext[0]-'0';  return NONZERO_DIGIT_TOKEN; }
        // RDKit✔️❌: \!			{ return NOT_TOKEN; }
        // RDKit✔️❌: \;			{ return SEMI_TOKEN; }
        // RDKit✔️❌: \&			{ return AND_TOKEN; }
        // RDKit✔️❌: \,			{ return OR_TOKEN; }
        let token = match ch {
            b':' => ScannerToken::Colon,
            b'_' => ScannerToken::Underscore,
            b'#' => ScannerToken::Hash,
            b'-' => ScannerToken::Minus,
            b'+' => ScannerToken::Plus,
            b'{' => ScannerToken::RangeOpen,
            b'}' => ScannerToken::RangeClose,
            b'.' => ScannerToken::Separator,
            b'%' => ScannerToken::Percent,
            b'0'..=b'9' => ScannerToken::Digit(ch - b'0'),
            b'!' => ScannerToken::Not,
            b';' => ScannerToken::Semi,
            b'&' => ScannerToken::And,
            b',' => ScannerToken::Or,
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
/// length. State operations are O(1); token storage is O(n), and the byte map
/// adds O(n) usize storage beyond the source scanner. Fixed element-rule
/// lookup has constant size, and compacted token ranges are not rescanned.
fn tokenize(input: impl AsRef<[u8]>) -> Result<Vec<(Token, SmartsTokenSpan)>, SmartsParseError> {
    let input = input.as_ref();
    generic_parse_helper(input, ScannerStart::Molecule)
}

fn generic_parse_helper(
    input: impl AsRef<[u8]>,
    start: ScannerStart,
) -> Result<Vec<(Token, SmartsTokenSpan)>, SmartsParseError> {
    let input = input.as_ref();
    // RDKit❗❌: int generic_parse_helper(T parser,
    // RDKit❗❌:                          const std::string &inp,
    // RDKit❗❌:                          std::vector<RDKit::RWMol *> &molVect,
    // RDKit❗❌:                          Atom *&atom,
    // RDKit❗❌:                          Bond *&bond,
    // RDKit❗❌:                          int start_tok,
    // RDKit❗❌:                          const std::string& input_type) {
    // RDKit❗❌:   TEST_ASSERT(!lex_init(&scanner));
    // RDKit❗❌:     ltrim = string_setup(inp, scanner);
    // RDKit❗❌:     unsigned int current_token_position = 0;
    // RDKit❗❌:     res = parser(inp.c_str() + ltrim, &molVect, atom, bond,
    // RDKit❗❌:                          numAtomsParsed, numBondsParsed, branchPoints, scanner,
    // RDKit❗❌:                          start_tok, current_token_position);
    // RDKit❗❌:   lex_destroy(scanner);
    // Local complexity review: setup, scanner position mapping and token
    // compaction are linear. Rust retains O(n) characters, byte boundaries,
    // and tokens; the boundary map is extra O(n) usize storage versus the
    // source's streaming scanner, while avoiding repeated prefix rescans.
    let input = input.as_ref();
    let window = setup_smarts_input(input);
    let scanned = SmartsScanner::new(input, start, window).scan()?;
    compact_scanned_tokens(input, &scanned)
}

fn compact_scanned_tokens(
    input: &[u8],
    scanned: &[ScannedToken],
) -> Result<Vec<(Token, SmartsTokenSpan)>, SmartsParseError> {
    // RDKit❗❌:     res = parser(inp.c_str() + ltrim, &molVect, atom, bond,
    // RDKit❗❌:                          numAtomsParsed, numBondsParsed, branchPoints, scanner,
    // RDKit❗❌:                          start_tok, current_token_position);
    // Local complexity review: token compaction is linear, but Rust
    // also collects an O(n) character vector for safe bracket slices;
    // RDKit consumes its scanner buffer directly. The byte boundary
    // map prevents repeated prefix scans but adds another O(n) vector.
    let chars = input;
    let mut tokens = Vec::new();
    let mut i = 1usize;
    while i < scanned.len() {
        let current = &scanned[i];
        match &current.token {
            ScannerToken::Start(_) => unreachable!("start token is first"),
            ScannerToken::OrganicElement(symbol) | ScannerToken::AtomElement(symbol) => {
                tokens.push((Token::OrganicElement(symbol.clone()), current.span));
            }
            ScannerToken::AromaticElement(symbol) => {
                tokens.push((Token::AromaticElement(symbol.clone()), current.span));
            }
            ScannerToken::SimpleAtomQuery(ch) => {
                tokens.push((Token::SimpleAtomQuery(*ch), current.span));
            }
            ScannerToken::BadCharacter(character) => {
                tokens.push((Token::BadCharacter(*character), current.span));
            }
            ScannerToken::AtomOpen => {
                let content_char_start = current.span.input_char_end;
                let content_parser_byte_start = current.span.parser_byte_end;
                let mut depth = 1usize;
                let mut cursor = i + 1;
                let mut bad_character_index = None;
                while cursor < scanned.len() && depth > 0 {
                    match scanned[cursor].token {
                        ScannerToken::AtomOpen => depth += 1,
                        ScannerToken::AtomClose => depth -= 1,
                        ScannerToken::BadCharacter(_) => {
                            bad_character_index = Some(cursor);
                            break;
                        }
                        _ => {}
                    }
                    cursor += 1;
                }
                if depth != 0 {
                    if let Some(bad_index) = bad_character_index {
                        let bad_character = &scanned[bad_index];
                        let content_char_end = bad_character.span.input_char_start;
                        if content_char_start < content_char_end {
                            let content = &chars[content_char_start..content_char_end];
                            let lexical_tokens = scanned[i + 1..bad_index]
                                .iter()
                                .filter(|token| {
                                    matches!(
                                        &token.token,
                                        ScannerToken::SimpleAtomQuery(_)
                                            | ScannerToken::AromaticElement(_)
                                    )
                                })
                                .map(|token| ScannedToken {
                                    token: token.token.clone(),
                                    span: SmartsTokenSpan {
                                        input_char_start: token.span.input_char_start
                                            - content_char_start,
                                        input_char_end: token.span.input_char_end
                                            - content_char_start,
                                        parser_byte_start: token.span.parser_byte_start
                                            - content_parser_byte_start,
                                        parser_byte_end: token.span.parser_byte_end
                                            - content_parser_byte_start,
                                    },
                                })
                                .collect();
                            let content_span = SmartsTokenSpan {
                                input_char_start: content_char_start,
                                input_char_end: content_char_end,
                                parser_byte_start: content_parser_byte_start,
                                parser_byte_end: bad_character.span.parser_byte_start,
                            };
                            let bracket_span = SmartsTokenSpan {
                                input_char_start: current.span.input_char_start,
                                input_char_end: content_char_end,
                                parser_byte_start: current.span.parser_byte_start,
                                parser_byte_end: bad_character.span.parser_byte_start,
                            };
                            tokens.push((
                                Token::BracketContent(BracketContent {
                                    text: content.to_vec(),
                                    span: content_span,
                                    lexical_tokens,
                                }),
                                bracket_span,
                            ));
                        }
                        let ScannerToken::BadCharacter(character) = &bad_character.token else {
                            unreachable!("the recorded terminal token is BAD_CHARACTER")
                        };
                        tokens.push((Token::BadCharacter(*character), bad_character.span));
                        return Ok(tokens);
                    }
                    return Err(SmartsParseError::UnclosedBracket(
                        current.span.parser_byte_start,
                    ));
                }
                let close = &scanned[cursor - 1];
                let content_char_end = close.span.input_char_start;
                let content = &chars[content_char_start..content_char_end];
                let lexical_tokens = scanned[i + 1..cursor - 1]
                    .iter()
                    .filter(|token| {
                        matches!(
                            &token.token,
                            ScannerToken::SimpleAtomQuery(_) | ScannerToken::AromaticElement(_)
                        )
                    })
                    .map(|token| ScannedToken {
                        token: token.token.clone(),
                        span: SmartsTokenSpan {
                            input_char_start: token.span.input_char_start - content_char_start,
                            input_char_end: token.span.input_char_end - content_char_start,
                            parser_byte_start: token.span.parser_byte_start
                                - content_parser_byte_start,
                            parser_byte_end: token.span.parser_byte_end - content_parser_byte_start,
                        },
                    })
                    .collect();
                let content_span = SmartsTokenSpan {
                    input_char_start: content_char_start,
                    input_char_end: content_char_end,
                    parser_byte_start: content_parser_byte_start,
                    parser_byte_end: close.span.parser_byte_start,
                };
                let bracket_span = SmartsTokenSpan {
                    input_char_start: current.span.input_char_start,
                    input_char_end: close.span.input_char_end,
                    parser_byte_start: current.span.parser_byte_start,
                    parser_byte_end: close.span.parser_byte_end,
                };
                tokens.push((
                    Token::BracketContent(BracketContent {
                        text: content.to_vec(),
                        span: content_span,
                        lexical_tokens,
                    }),
                    bracket_span,
                ));
                i = cursor;
                continue;
            }
            ScannerToken::BondSpec(ch) => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(*ch)), current.span))
            }
            ScannerToken::DativeRight => {
                tokens.push((Token::BondSpec(BondLexeme::DativeRight), current.span));
            }
            ScannerToken::DativeLeft => {
                tokens.push((Token::BondSpec(BondLexeme::DativeLeft), current.span));
            }
            ScannerToken::At => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(b'@')), current.span))
            }
            ScannerToken::Colon => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(b':')), current.span))
            }
            ScannerToken::Hash => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(b'#')), current.span))
            }
            ScannerToken::Minus => {
                tokens.push((Token::BondSpec(BondLexeme::Symbol(b'-')), current.span))
            }
            ScannerToken::GroupOpen => tokens.push((Token::OpenParen, current.span)),
            ScannerToken::GroupClose => tokens.push((Token::CloseParen, current.span)),
            ScannerToken::Separator => tokens.push((Token::Dot, current.span)),
            ScannerToken::Digit(value) => {
                tokens.push((Token::RingClosureDigit(u32::from(*value)), current.span));
            }
            ScannerToken::Percent => {
                // RDKit❗✔️: ring_number:  digit
                // RDKit❗✔️: | PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
                // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
                // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
                // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
                // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
                // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
                let (number, consumed) = compact_ring_number(&chars, scanned, i)?;
                let last_span = scanned[consumed - 1].span;
                let mut number_span = current.span;
                number_span.input_char_end = last_span.input_char_end;
                number_span.parser_byte_end = last_span.parser_byte_end;
                tokens.push((Token::RingClosurePercent(number), number_span));
                i = consumed;
                continue;
            }
            ScannerToken::Not => tokens.push((Token::Not, current.span)),
            ScannerToken::Semi => tokens.push((Token::Semi, current.span)),
            ScannerToken::And => tokens.push((Token::And, current.span)),
            ScannerToken::Or => tokens.push((Token::Or, current.span)),
            ScannerToken::EndOfStream => tokens.push((Token::EndOfStream, current.span)),
            token => {
                return Err(SmartsParseError::UnexpectedCharacter {
                    position: current.span.parser_byte_end,
                    character: chars
                        .get(current.span.input_char_start)
                        .copied()
                        .expect("non-EOS scanner token owns a source byte"),
                    context: format!("unexpected {token:?} token in molecule SMARTS"),
                });
            }
        }
        i += 1;
    }
    Ok(tokens)
}

fn invalid_percent(
    chars: &[u8],
    scanned: &[ScannedToken],
    unexpected_index: usize,
) -> SmartsParseError {
    // RDKit❗✔️: #define YY_USER_ACTION current_token_position += yyleng;
    // The parser position is the consumed byte endpoint of the first token
    // that cannot continue a `ring_number` production. An EOS token has no
    // character to report; BAD_CHARACTER uses the scanner's shared mapping.
    let Some(unexpected) = scanned.get(unexpected_index) else {
        return SmartsParseError::UnexpectedEnd("expected ring closure number".to_string());
    };
    match unexpected.token {
        ScannerToken::EndOfStream => {
            SmartsParseError::UnexpectedEnd("expected ring closure number".to_string())
        }
        ScannerToken::BadCharacter(character) => {
            SmartsParser::bad_character_error(unexpected.span, character)
        }
        _ => SmartsParseError::UnexpectedCharacter {
            position: unexpected.span.parser_byte_end,
            character: *chars
                .get(unexpected.span.input_char_start)
                .expect("non-EOS scanner token has a source character"),
            context: "invalid ring closure number".to_string(),
        },
    }
}

fn compact_ring_number(
    chars: &[u8],
    scanned: &[ScannedToken],
    percent_index: usize,
) -> Result<(u32, usize), SmartsParseError> {
    // RDKit❗✔️: ring_number:  digit
    // RDKit❗✔️: | PERCENT_TOKEN NONZERO_DIGIT_TOKEN digit { $$ = $2*10+$3; }
    // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit GROUP_CLOSE_TOKEN { $$ = $3; }
    // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit GROUP_CLOSE_TOKEN { $$ = $3*10+$4; }
    // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*100+$4*10+$5; }
    // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*1000+$4*100+$5*10+$6; }
    // RDKit❗✔️: | PERCENT_TOKEN GROUP_OPEN_TOKEN digit digit digit digit digit GROUP_CLOSE_TOKEN { $$ = $3*10000+$4*1000+$5*100+$6*10+$7; }
    // Local complexity review: inspect at most five grouped digits or two
    // shorthand digits; valid parsing stays allocation-free and O(1) bounded.
    let Some(next) = scanned.get(percent_index + 1) else {
        return Err(invalid_percent(chars, scanned, percent_index + 1));
    };
    let (digits_start, close_required) = match next.token {
        ScannerToken::Digit(value) if value != 0 => {
            let _ = value;
            (percent_index + 1, false)
        }
        ScannerToken::GroupOpen => (percent_index + 2, true),
        _ => {
            return Err(invalid_percent(chars, scanned, percent_index + 1));
        }
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
        return Err(invalid_percent(chars, scanned, cursor));
    }
    let consumed = if close_required { cursor + 1 } else { cursor };
    Ok((value, consumed))
}

fn invalid_atom_operator(position: usize, operator: u8) -> SmartsParseError {
    SmartsParseError::InvalidAtomPrimitive {
        position,
        detail: format!("operator '{operator}' has no left operand"),
    }
}

// ---------------------------------------------------------------------------
// Recursive-descent SMARTS Parser
// ---------------------------------------------------------------------------

/// Recursive-descent SMARTS parser for the currently modeled grammar.
struct SmartsParser<'a> {
    tokens: &'a [(Token, SmartsTokenSpan)],
    input: &'a [u8],
    pos: usize,
    /// Preserve every occurrence in the source's sorted-label/token order.
    ring_closure_targets: BTreeMap<u32, Vec<RingClosureOccurrence>>,
    /// One grammar invocation shares this counter through recursive atomd.
    num_bonds_parsed: Rc<Cell<u32>>,
}

struct ParsedSmartsAtom {
    carrier: QueryAtom,
    atom_map: Option<u32>,
}

struct SimpleAtom {
    query: QueryNode<AtomQueryPredicate>,
    atomic_number: u8,
    aromatic: bool,
}

#[derive(Debug, Clone)]
struct ParsedSmartsBond {
    query: QueryNode<BondQueryPredicate>,
    carrier_order: BondOrder,
    /// Mirrors RDKit's `_unspecifiedOrder` marker independently of type/query.
    unspecified_order: bool,
}

#[derive(Debug, Clone)]
struct RingClosureOccurrence {
    atom_idx: usize,
    bond: ParsedSmartsBond,
    direction: BondDirection,
    source_index: Option<u32>,
}

impl ParsedSmartsBond {
    fn expand_query(&mut self, other: Self, how: CompositeQueryType) {
        // BEGIN RDKIT CPP FUNCTION bond_query reduction
        // RDKit❗✔️: bond_query: bondd
        // RDKit❗✔️: | bond_query bondd {
        // RDKit❗✔️:   $1->expandQuery($2->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit❗✔️:   delete $2;
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION bond_query reduction
        // Each grammar reduction keeps the left QueryBond and only expands
        // its predicate. The typed query helper mirrors that action; moving
        // the left carrier order unchanged preserves the source field.
        query_bond_expand_query(&mut self.query, other.query, how, true);
    }
}

impl SimpleAtom {
    fn into_parsed_atom(self, atom_map: Option<u32>) -> ParsedSmartsAtom {
        ParsedAtomExpr::with_identity(self.query, self.atomic_number, self.aromatic)
            .into_parsed_atom(atom_map)
    }
}

#[derive(Debug, Clone)]
struct ParsedAtomExpr {
    carrier: QueryAtom,
    hydrogen_mask: bool,
    charge_mask: bool,
}

impl ParsedAtomExpr {
    fn query_only(query: QueryNode<AtomQueryPredicate>) -> Self {
        Self {
            carrier: QueryAtom::from_identity_parts(
                AtomId::new(0),
                QueryAtomIdentity::Element(Element::DUMMY),
                query,
            ),
            hydrogen_mask: false,
            charge_mask: false,
        }
    }

    fn with_identity(
        query: QueryNode<AtomQueryPredicate>,
        atomic_number: u8,
        aromatic: bool,
    ) -> Self {
        // BEGIN RDKIT CPP FUNCTION QueryAtom::QueryAtom(int num)
        // RDKit❗✔️: explicit QueryAtom(int num) : Atom(num), dp_query(makeAtomNumQuery(num)) {}
        // END RDKIT CPP FUNCTION QueryAtom::QueryAtom(int num)
        let mut carrier = QueryAtom::from_identity_parts(
            AtomId::new(0),
            QueryAtomIdentity::from_atomic_number(atomic_number),
            query,
        );
        carrier.set_aromatic(aromatic);
        Self {
            carrier,
            hydrogen_mask: false,
            charge_mask: false,
        }
    }

    fn clear_chemical_properties(&mut self) {
        // BEGIN RDKIT CPP FUNCTION SmilesParseOps::ClearAtomChemicalProps
        // RDKit✔️✔️:   atom->setIsotope(0);
        // RDKit✔️✔️:   atom->setFormalCharge(0);
        // RDKit✔️✔️:   atom->setNumExplicitHs(0);
        // END RDKIT CPP FUNCTION SmilesParseOps::ClearAtomChemicalProps
        // The source clear deliberately leaves noImplicit and parser masks
        // untouched. Apply only the three named carrier writes at each
        // reduction, without scanning the final predicate tree.
        self.carrier.set_isotope(None);
        self.carrier.set_formal_charge(0);
        self.carrier.set_explicit_hydrogens(0);
    }

    fn reset_atomic_number(mut self) -> Self {
        // BEGIN RDKIT CPP FUNCTION atom_expr OR carrier action
        // RDKit✔️✔️: $1->setAtomicNum(0);
        // END RDKIT CPP FUNCTION atom_expr OR carrier action
        self.carrier = self
            .carrier
            .with_identity(QueryAtomIdentity::Element(Element::DUMMY));
        self
    }

    fn reduce_atom_expr(mut self, other: Self, how: CompositeQueryType) -> Self {
        // BEGIN RDKIT CPP FUNCTION atom_expr reduction
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
        // RDKit✔️✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit✔️✔️:   $1->setAtomicNum(0);
        // RDKit✔️✔️:   $$ = $1;
        // END RDKIT CPP FUNCTION atom_expr reduction
        // All atom_expr reductions keep the left QueryAtom. AND and SEMI
        // clear its three chemical properties; OR additionally resets only
        // atomic number and retains aromaticity and noImplicit.
        let ParsedAtomExpr {
            carrier: mut other_carrier,
            ..
        } = other;
        let other_query =
            std::mem::replace(other_carrier.predicate_mut(), QueryNode::and(Vec::new()));
        crate::query_behavior::query_atom_expand_query(
            self.carrier.predicate_mut(),
            other_query,
            how,
            true,
        );
        // RDKit❗✔️: atom_expr: atom_expr AND_TOKEN atom_expr {
        // RDKit❗✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit❗✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit❗✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit❗✔️:   delete $3;
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        // RDKit❗✔️: | atom_expr OR_TOKEN atom_expr {
        // RDKit❗✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_OR,true);
        // RDKit❗✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit❗✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit❗✔️:   $1->setAtomicNum(0);
        // RDKit❗✔️:   delete $3;
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        // RDKit❗✔️: | atom_expr SEMI_TOKEN atom_expr {
        // RDKit❗✔️:   $1->expandQuery($3->getQuery()->copy(),Queries::COMPOSITE_AND,true);
        // RDKit❗✔️:   if ($1->getChiralTag()==Atom::CHI_UNSPECIFIED) { $1->setChiralTag($3->getChiralTag()); }
        // RDKit❗✔️:   SmilesParseOps::ClearAtomChemicalProps($1);
        // RDKit❗✔️:   delete $3;
        if self.carrier.chiral_tag() == ChiralTag::Unspecified {
            self.carrier.set_chiral_tag(other_carrier.chiral_tag());
        }
        // atom_expr reductions deliberately keep only the left permutation.
        self.clear_chemical_properties();
        if how == CompositeQueryType::Or {
            self = self.reset_atomic_number();
        }
        // The source carries the left QueryAtom's flags through reduction.
        self
    }

    fn and_point_query(mut self, point_query: Self) -> Self {
        // BEGIN RDKIT CPP FUNCTION atom_expr_and_point_query
        // RDKit✔️✔️:     atom_expr->expandQuery(point_query->getQuery()->copy(), Queries::COMPOSITE_AND, true);
        // RDKit✔️✔️:     if (point_query->getFlags() & SMARTS_H_MASK) {
        // RDKit✔️✔️:       if (!(atom_expr->getFlags() & SMARTS_H_MASK)) {
        // RDKit✔️✔️:         atom_expr->setNumExplicitHs(point_query->getNumExplicitHs());
        // RDKit✔️✔️:         atom_expr->setNoImplicit(true);
        // RDKit✔️✔️:         atom_expr->getFlags() |= SMARTS_H_MASK;
        // RDKit✔️✔️:       } else if (atom_expr->getNumExplicitHs() != point_query->getNumExplicitHs()) {
        // RDKit✔️✔️:         atom_expr->setNumExplicitHs(0);
        // RDKit✔️✔️:         atom_expr->setNoImplicit(true);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (point_query->getFlags() & SMARTS_CHARGE_MASK) {
        // RDKit✔️✔️:       if (!(atom_expr->getFlags() & SMARTS_CHARGE_MASK)) {
        // RDKit✔️✔️:         atom_expr->setFormalCharge(point_query->getFormalCharge());
        // RDKit✔️✔️:         atom_expr->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit✔️✔️:       } else if (atom_expr->getFormalCharge() != point_query->getFormalCharge()) {
        // RDKit✔️✔️:         atom_expr->setFormalCharge(0);
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // END RDKIT CPP FUNCTION atom_expr_and_point_query
        // This grammar action transfers only right-side H/charge carrier
        // fields indicated by its source masks; conflicts retain both query
        // leaves but clear the carrier field. Predicate construction remains
        // delegated to the ordered QueryAtom::expandQuery port.
        let ParsedAtomExpr {
            carrier: mut point_carrier,
            hydrogen_mask: point_hydrogen_mask,
            charge_mask: point_charge_mask,
        } = point_query;
        let point_hydrogens = point_carrier.explicit_hydrogens();
        let point_charge = point_carrier.formal_charge();
        let point_query =
            std::mem::replace(point_carrier.predicate_mut(), QueryNode::and(Vec::new()));
        crate::query_behavior::query_atom_expand_query(
            self.carrier.predicate_mut(),
            point_query,
            CompositeQueryType::And,
            true,
        );
        // RDKit❗✔️:     atom_expr->expandQuery(point_query->getQuery()->copy(), Queries::COMPOSITE_AND, true);
        // RDKit❗✔️:     if (atom_expr->getChiralTag() == Atom::CHI_UNSPECIFIED) {
        // RDKit❗✔️:       atom_expr->setChiralTag(point_query->getChiralTag());
        // RDKit❗✔️:       int perm;
        // RDKit❗✔️:       if (point_query->getPropIfPresent(common_properties::_chiralPermutation, perm)) {
        // RDKit❗✔️:         atom_expr->setProp(common_properties::_chiralPermutation, perm);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     if (point_query->getFlags() & SMARTS_H_MASK) {
        if self.carrier.chiral_tag() == ChiralTag::Unspecified {
            self.carrier.set_chiral_tag(point_carrier.chiral_tag());
            if let Some(permutation) = point_carrier.chiral_permutation() {
                self.carrier.set_chiral_permutation(Some(permutation));
            }
        }
        if point_hydrogen_mask {
            if !self.hydrogen_mask {
                self.carrier.set_explicit_hydrogens(point_hydrogens);
                self.carrier.set_no_implicit(true);
                self.hydrogen_mask = true;
            } else if self.carrier.explicit_hydrogens() != point_hydrogens {
                self.carrier.set_explicit_hydrogens(0);
                self.carrier.set_no_implicit(true);
            }
        }
        if point_charge_mask {
            if !self.charge_mask {
                self.carrier.set_formal_charge(point_charge);
                self.charge_mask = true;
            } else if self.carrier.formal_charge() != point_charge {
                self.carrier.set_formal_charge(0);
            }
        }
        self
    }

    fn into_parsed_atom(self, atom_map: Option<u32>) -> ParsedSmartsAtom {
        ParsedSmartsAtom {
            carrier: self.carrier,
            atom_map,
        }
    }
}

fn split_atom_map_suffix(content: &[u8]) -> Result<(&[u8], Option<u32>), SmartsParseError> {
    let bytes = content;
    let mut digit_start = bytes.len();
    while digit_start > 0 && bytes[digit_start - 1].is_ascii_digit() {
        digit_start -= 1;
    }
    if digit_start == bytes.len() || digit_start == 0 || bytes[digit_start - 1] != b':' {
        return Ok((content, None));
    }
    let colon = digit_start - 1;
    let atom_map = content[digit_start..]
        .iter()
        .try_fold(0_u32, |value, digit| {
            value
                .checked_mul(10)
                .and_then(|value| value.checked_add(u32::from(*digit - b'0')))
                .ok_or(())
        })
        .map_err(|_| SmartsParseError::InvalidAtomPrimitive {
            position: digit_start,
            detail: "atom map number is out of range".to_string(),
        })?;
    Ok((&content[..colon], Some(atom_map)))
}

impl<'a> SmartsParser<'a> {
    fn new(tokens: &'a [(Token, SmartsTokenSpan)], input: &'a [u8]) -> Self {
        Self {
            tokens,
            input: input.as_ref(),
            pos: 0,
            ring_closure_targets: BTreeMap::new(),
            num_bonds_parsed: Rc::new(Cell::new(0)),
        }
    }

    fn next_bond_source_index(&self) -> Result<u32, SmartsParseError> {
        // RDKit❗✔️: $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // One shared counter follows grammar reductions, including recursive
        // atomd before its parent bond. Allocation is O(1), with no graph scan.
        // Source unsigned32 postincrement wraps at UINT_MAX; this counter is
        // source encounter metadata, never a persistent BondId.
        let index = self.num_bonds_parsed.get();
        let next = index.wrapping_add(1);
        self.num_bonds_parsed.set(next);
        Ok(index)
    }

    fn peek(&self) -> &(Token, SmartsTokenSpan) {
        &self.tokens[self.pos]
    }

    fn advance(&mut self) {
        self.pos += 1;
    }

    fn input_character_position(&self) -> usize {
        // These custom Rust atom-primitive diagnostics retain their existing
        // helper-input character-index projection; lexer/parser diagnostics
        // use the separately named `source_error_position` byte coordinate.
        self.tokens[self.pos].1.input_char_start
    }

    fn source_error_position(&self) -> usize {
        // RDKit❗✔️:   yyerror(input, molList, current_token_position, "syntax error");
        // Parser diagnostics use the scanner's post-consumption byte counter,
        // while `input_character_position` remains the helper-input character
        // projection used by the existing custom atom-primitive errors.
        self.tokens[self.pos].1.parser_byte_end
    }

    fn bad_character_error(span: SmartsTokenSpan, character: u8) -> SmartsParseError {
        // RDKit❗✔️: | meta_start BAD_CHARACTER {
        // RDKit❗✔️:   yyerrok;
        // RDKit❗✔️:   yyErrorCleanup(molList);
        // RDKit❗✔️:   yyerror(input, molList, current_token_position, "syntax error");
        // RDKit❗✔️:   YYABORT;
        // RDKit❗✔️: }
        SmartsParseError::UnexpectedCharacter {
            position: span.parser_byte_end,
            character,
            context: "unexpected character in SMARTS string".to_string(),
        }
    }

    fn unexpected_token_error(&self, span: SmartsTokenSpan, context: &str) -> SmartsParseError {
        // RDKit❗✔️: void
        // RDKit❗✔️: yysmarts_error( const char *input,
        // RDKit❗✔️:                 std::vector<RDKit::RWMol *> *ms,
        // RDKit❗✔️:                 RDKit::Atom* &,
        // RDKit❗✔️:                 RDKit::Bond* &,
        // RDKit❗✔️:                 unsigned int &,
        // RDKit❗✔️:                 unsigned int &,
        // RDKit❗✔️:                 std::vector<std::pair<unsigned int, unsigned int>>&,
        // RDKit❗✔️:                 void *,
        // RDKit❗✔️:                 int,
        // RDKit❗✔️:                 unsigned int bad_token_position,
        // RDKit❗✔️:                 const char *msg  )
        // RDKit❗✔️: {
        // RDKit❗✔️:   yyErrorCleanup(ms);
        // RDKit❗✔️:   SmilesParseOps::detail::printSyntaxErrorMessage(input,
        // RDKit❗✔️:                                                   msg,
        // RDKit❗✔️:                                                   bad_token_position,
        // RDKit❗✔️:                                                   "SMARTS");
        // RDKit❗✔️: }
        // A typed diagnostic preserves EOS absence separately from the exact
        // token byte. It never substitutes Debug/Display text or a default byte.
        if matches!(self.peek().0, Token::EndOfStream) {
            SmartsParseError::UnexpectedEnd(context.to_owned())
        } else {
            SmartsParseError::UnexpectedCharacter {
                position: span.parser_byte_end,
                character: self.input[span.input_char_start],
                context: context.to_owned(),
            }
        }
    }

    fn require_end(&self, context: &str) -> Result<(), SmartsParseError> {
        // RDKit❗✔️:   yyerror(input, molList, current_token_position, "syntax error");
        match self.peek() {
            (Token::EndOfStream, _) => Ok(()),
            (Token::BadCharacter(character), span) => {
                Err(Self::bad_character_error(*span, *character))
            }
            (_token, span) => Err(SmartsParseError::UnexpectedCharacter {
                position: span.parser_byte_end,
                character: self.input[span.input_char_start],
                context: format!("unexpected trailing token in {context}"),
            }),
        }
    }

    /// Parse the full SMARTS pattern into a private graph.
    ///
    /// The top-level parser builds a molecule from source-ordered productions.
    fn parse_smarts_molecule(&mut self) -> Result<QueryGraphBuilder, SmartsParseError> {
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
        let mut graph = QueryGraphBuilder::default();

        // RDKit✔️✔️: mol: atomd {
        let first = self.parse_atomd()?;
        graph.push_atom(first);
        // RDKit✔️✔️: $1->setProp(RDKit::common_properties::_SmilesStart,1);
        graph.atoms[0]
            .set_prop(crate::query_graph_behavior::SMILES_START_PROP, 1_i32)
            .expect("the internal source-start property key is non-empty");

        // RDKit✔️✔️: | mol atomd       {
        let _ = self.parse_smarts_chain(&mut graph, 0)?;

        self.require_end("molecule SMARTS")?;
        self.close_ring_closures(&mut graph)?;

        Ok(graph)
    }

    /// Parse source molecule reductions beginning with the first atom.
    /// Parse atom, bond, branch, and ring-closure tokens into the graph.
    fn parse_smarts_chain(
        &mut self,
        graph: &mut QueryGraphBuilder,
        mut active_atom_idx: usize,
    ) -> Result<usize, SmartsParseError> {
        // Local complexity review: each token is consumed once by this loop
        // or one nested branch call. Branch-start validation is O(1); total
        // time is O(n), graph storage O(n), and recursion space O(branch depth).
        loop {
            match self.peek() {
                (Token::EndOfStream, _) => break,
                (Token::CloseParen, _) => break,

                // Bond spec followed by atom
                (Token::BondSpec(_), _) | (Token::Not, _) | (Token::And, _) | (Token::Semi, _) => {
                    let direction = self.current_bond_direction();
                    let bond = self.parse_bond_expr()?;
                    match self.peek() {
                        (Token::RingClosureDigit(n), _) | (Token::RingClosurePercent(n), _) => {
                            let num = *n;
                            self.advance();
                            // RDKit✔️✔️: | mol bond_expr ring_number {
                            // RDKit✔️✔️:   RWMol * mp = (*molList)[$$];
                            // RDKit✔️✔️:   Atom *atom=mp->getActiveAtom();
                            // RDKit✔️✔️:   mp->setBondBookmark($2,$3);
                            // RDKit✔️✔️:   $2->setOwningMol(mp);
                            // RDKit✔️✔️:   $2->setBeginAtomIdx(atom->getIdx());
                            // RDKit✔️✔️:   $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
                            // RDKit✔️✔️:   mp->setAtomBookmark(atom,$3);
                            self.record_ring_closure(
                                graph,
                                num,
                                active_atom_idx,
                                bond,
                                direction,
                                true,
                            )?;
                        }
                        _ => {
                            let (bond, reverse_endpoints) = normalize_dative_bond(bond);
                            let atom = self.parse_atomd()?;
                            let end_atom_idx = graph.push_atom(atom);
                            let source_index = self.next_bond_source_index()?;
                            if reverse_endpoints {
                                graph.push_bond(
                                    end_atom_idx,
                                    active_atom_idx,
                                    bond,
                                    direction,
                                    source_index,
                                );
                            } else {
                                graph.push_bond(
                                    active_atom_idx,
                                    end_atom_idx,
                                    bond,
                                    direction,
                                    source_index,
                                );
                            }
                            active_atom_idx = end_atom_idx;
                        }
                    }
                }

                // No bond spec — implicit single/aromatic bond (SMARTS semantics)
                _ => {
                    // Check if next is a ring closure or branch first
                    match self.peek() {
                        (Token::RingClosureDigit(n), _) | (Token::RingClosurePercent(n), _) => {
                            let num = *n;
                            self.advance();
                            // Record ring closure on the parser's active atom.
                            self.record_ring_closure(
                                graph,
                                num,
                                active_atom_idx,
                                ParsedSmartsBond {
                                    query: unspecified_smarts_bond_query(),
                                    carrier_order: BondOrder::Unspecified,
                                    unspecified_order: true,
                                },
                                BondDirection::None,
                                false,
                            )?;
                        }
                        (Token::OpenParen, _) => {
                            let _branch_position = self.parse_branch_open_token()?;
                            // RDKit source: smarts.yy branch productions begin
                            // with atomd or bond_expr atomd; a branch cannot
                            // start with a separator, ring closure, or paren.
                            // RDKit✔️✔️: | mol branch_open_token atomd {
                            // RDKit✔️✔️: | mol branch_open_token bond_expr atomd {
                            if matches!(
                                self.peek(),
                                (
                                    Token::OpenParen
                                        | Token::CloseParen
                                        | Token::RingClosureDigit(_)
                                        | Token::RingClosurePercent(_)
                                        | Token::Dot
                                        | Token::EndOfStream,
                                    _
                                )
                            ) {
                                return Err(self.unexpected_token_error(
                                    self.peek().1,
                                    "expected atom expression",
                                ));
                            }
                            // RDKit✔️✔️:   branchPoints.push_back({atomIdx1, $2});
                            // RDKit✔️✔️: | mol GROUP_CLOSE_TOKEN {
                            // RDKit✔️✔️:   mp->setActiveAtom(branchPoints.back().first);
                            // RDKit✔️✔️:   branchPoints.pop_back();
                            let _branch_active = self.parse_smarts_chain(graph, active_atom_idx)?;
                            match self.peek() {
                                (Token::CloseParen, _) => {
                                    self.advance();
                                }
                                (_tok, pos) => {
                                    // RDKit❗✔️:   yyerror(input, molList, current_token_position, "syntax error");
                                    return Err(self.unexpected_token_error(
                                        *pos,
                                        "expected close parenthesis",
                                    ));
                                }
                            }
                        }
                        (Token::Dot, _) => {
                            // RDKit✔️✔️: | mol SEPARATOR_TOKEN atomd {
                            self.advance();
                            let atom = self.parse_atomd()?;
                            active_atom_idx = graph.push_atom(atom);
                            // RDKit✔️✔️: $3->setProp(RDKit::common_properties::_SmilesStart,1,true);
                            graph.atoms[active_atom_idx]
                                .set_computed_prop(
                                    crate::query_graph_behavior::SMILES_START_PROP,
                                    1_i32,
                                )
                                .expect("the internal source-start property key is non-empty");
                        }
                        // Atom follows implicitly with default bond
                        _ => {
                            let atom = self.parse_atomd()?;
                            let end_atom_idx = graph.push_atom(atom);
                            graph.push_bond(
                                active_atom_idx,
                                end_atom_idx,
                                ParsedSmartsBond {
                                    query: unspecified_smarts_bond_query(),
                                    carrier_order: graph
                                        .implicit_bond_order(active_atom_idx, end_atom_idx),
                                    unspecified_order: true,
                                },
                                BondDirection::None,
                                self.next_bond_source_index()?,
                            );
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
        graph: &mut QueryGraphBuilder,
        num: u32,
        atom_idx: usize,
        bond: ParsedSmartsBond,
        direction: BondDirection,
        explicit: bool,
    ) -> Result<(), SmartsParseError> {
        // RDKit✔️✔️: mp->setBondBookmark(newB,$2);
        // RDKit✔️✔️: mp->setAtomBookmark(atom,$2);
        // RDKit✔️✔️: mp->setBondBookmark($2,$3);
        // RDKit✔️✔️: mp->setAtomBookmark(atom,$3);
        // Each token occurrence is appended once to its source-label bucket;
        // BTreeMap preserves CloseMolRings label order and Vec preserves the
        // source atom occurrence order without rescanning parser input.
        // Local complexity: O(log L) label lookup and amortized O(1) append.
        // RDKit✔️✔️: if(!(mp->getAllBondsWithBookmark($2).size()%2)){
        // RDKit✔️✔️:   newB->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: $2->setProp("_cxsmilesBondIdx",numBondsParsed++);
        // The implicit action checks bookmark count after appending; explicit
        // bond_expr allocates at every occurrence, including opening tokens.
        // RDKit❗✔️: tmp.push_back(-($2+1));
        // Preserve the source token-position placeholders on the atom carrier;
        // explicit wrapping spells out the source signed carrier conversion.
        let atom = &mut graph.atoms[atom_idx];
        let mut closures = match atom.prop("_RingClosures") {
            Some(cosmolkit_model::PropertyValue::IntVector(ids)) => ids.clone(),
            None => Vec::new(),
            Some(_) => {
                return Err(SmartsParseError::Parse(
                    "ring closure property has the wrong type".into(),
                ));
            }
        };
        closures.push((num as i32).wrapping_add(1).wrapping_neg());
        atom.set_prop("_RingClosures", closures)
            .expect("the internal ring-closure property key is non-empty");
        let count_before = self.ring_closure_targets.get(&num).map_or(0, Vec::len);
        let source_index = if explicit || count_before % 2 == 1 {
            Some(self.next_bond_source_index()?)
        } else {
            None
        };
        self.ring_closure_targets
            .entry(num)
            .or_default()
            .push(RingClosureOccurrence {
                atom_idx,
                bond,
                direction,
                source_index,
            });
        Ok(())
    }

    fn close_ring_closures(
        &mut self,
        graph: &mut QueryGraphBuilder,
    ) -> Result<(), SmartsParseError> {
        // BEGIN RDKIT CPP FUNCTION SmilesParseOps::CloseMolRings
        // RDKit❗❌: void CloseMolRings(RWMol *mol, bool toleratePartials) {
        // RDKit❗❌:   auto bookmarkIt = mol->getAtomBookmarks()->begin();
        // RDKit❗❌:   while (bookmarkIt != mol->getAtomBookmarks()->end()) {
        // RDKit❗❌:     auto &bookmark = *bookmarkIt;
        // RDKit❗❌:     auto atomIt = bookmark.second.begin();
        // RDKit❗❌:     auto atomsEnd = bookmark.second.end();
        // RDKit❗❌:     while (atomIt != atomsEnd) {
        // RDKit❗❌:       Atom *atom1 = *atomIt;
        // RDKit❗❌:       ++atomIt;
        // RDKit❗❌:       if (!toleratePartials && atomIt == atomsEnd) {
        // RDKit❗❌:         ReportParseError("unclosed ring");
        // RDKit❗❌:       } else if (atomIt != atomsEnd && *atomIt == atom1) {
        // RDKit❗❌:         auto fmt =
        // RDKit❗❌:             boost::format{
        // RDKit❗❌:                 "duplicated ring closure %1% bonds atom %2% to itself"} %
        // RDKit❗❌:             bookmark.first % atom1->getIdx();
        // RDKit❗❌:         std::string msg = fmt.str();
        // RDKit❗❌:         ReportParseError(msg.c_str(), true);
        // RDKit❗❌:       } else if (mol->getBondBetweenAtoms(atom1->getIdx(),
        // RDKit❗❌:                                             (*atomIt)->getIdx()) != nullptr) {
        // RDKit❗❌:         auto fmt =
        // RDKit❗❌:             boost::format{
        // RDKit❗❌:                 "ring closure %1% duplicates bond between atom %2% and atom "
        // RDKit❗❌:                 "%3%"} %
        // RDKit❗❌:             bookmark.first % atom1->getIdx() % (*atomIt)->getIdx();
        // RDKit❗❌:         std::string msg = fmt.str();
        // RDKit❗❌:         ReportParseError(msg.c_str(), true);
        // RDKit❗❌:       } else if (atomIt != atomsEnd) {
        // RDKit❗❌:         Atom *atom2 = *atomIt;
        // RDKit❗❌:         ++atomIt;
        // RDKit❗❌:         int bondIdx = -1;
        // RDKit❗❌:         // We're guaranteed two partial bonds, one for each time
        // RDKit❗❌:         // the ring index was used.  We give the first specification
        // RDKit❗❌:         // priority.
        // RDKit❗❌:         CHECK_INVARIANT(mol->hasBondBookmark(bookmark.first),
        // RDKit❗❌:                         "Missing bond bookmark");
        // RDKit❗❌:         RWMol::BOND_PTR_LIST bonds =
        // RDKit❗❌:             mol->getAllBondsWithBookmark(bookmark.first);
        // RDKit❗❌:         auto bondIt = bonds.begin();
        // RDKit❗❌:         CHECK_INVARIANT(bonds.size() >= 2, "Missing bond");
        // RDKit❗❌:         Bond *bond1 = *bondIt;
        // RDKit❗❌:         ++bondIt;
        // RDKit❗❌:         Bond *bond2 = *bondIt;
        // RDKit❗❌:         CHECK_INVARIANT(bond1->getBeginAtomIdx() == atom1->getIdx(),
        // RDKit❗❌:                         "bad begin atom");
        // RDKit❗❌:         CHECK_INVARIANT(bond2->getBeginAtomIdx() == atom2->getIdx(),
        // RDKit❗❌:                         "bad begin atom");
        // RDKit❗❌:       Bond *matchedBond;
        // RDKit❗❌:       if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
        // RDKit❗❌:         matchedBond = bond1;
        // RDKit❗❌:         if (matchedBond->getBondType() == Bond::DATIVEL) {
        // RDKit❗❌:           matchedBond->setBeginAtomIdx(atom2->getIdx());
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom1->getIdx());
        // RDKit❗❌:           matchedBond->setBondType(Bond::DATIVE);
        // RDKit❗❌:         } else if (matchedBond->getBondType() == Bond::DATIVER) {
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit❗❌:           matchedBond->setBondType(Bond::DATIVE);
        // RDKit❗❌:         } else {
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit❗❌:         }
        // RDKit❗❌:         swapBondDirIfNeeded(bond1, bond2);
        // RDKit❗❌:         delete bond2;
        // RDKit❗❌:       } else {
        // RDKit❗❌:         matchedBond = bond2;
        // RDKit❗❌:         if (matchedBond->getBondType() == Bond::DATIVEL) {
        // RDKit❗❌:           matchedBond->setBeginAtomIdx(atom1->getIdx());
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom2->getIdx());
        // RDKit❗❌:           matchedBond->setBondType(Bond::DATIVE);
        // RDKit❗❌:         } else if (matchedBond->getBondType() == Bond::DATIVER) {
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom1->getIdx());
        // RDKit❗❌:           matchedBond->setBondType(Bond::DATIVE);
        // RDKit❗❌:         } else {
        // RDKit❗❌:           matchedBond->setEndAtomIdx(atom1->getIdx());
        // RDKit❗❌:         }
        // RDKit❗❌:         swapBondDirIfNeeded(bond2, bond1);
        // RDKit❗❌:         delete bond1;
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // END RDKIT CPP FUNCTION SmilesParseOps::CloseMolRings
        // BEGIN RDKIT CPP HELPER swapBondDirIfNeeded
        // RDKit❗❌: void swapBondDirIfNeeded(Bond *bond1, const Bond *bond2) {
        // RDKit❗❌:   PRECONDITION(bond1, "bad bond1");
        // RDKit❗❌:   PRECONDITION(bond2, "bad bond2");
        // RDKit❗❌:   if (bond1->getBondDir() == Bond::NONE && bond2->getBondDir() != Bond::NONE) {
        // RDKit❗❌:     bond1->setBondDir(bond2->getBondDir());
        // RDKit❗❌:     if (bond1->getBeginAtom() != bond2->getBeginAtom()) {
        // RDKit❗❌:       switch (bond1->getBondDir()) {
        // RDKit❗❌:         case Bond::ENDDOWNRIGHT:
        // RDKit❗❌:           bond1->setBondDir(Bond::ENDUPRIGHT);
        // RDKit❗❌:           break;
        // RDKit❗❌:         case Bond::ENDUPRIGHT:
        // RDKit❗❌:           bond1->setBondDir(Bond::ENDDOWNRIGHT);
        // RDKit❗❌:           break;
        // RDKit❗❌:         default:
        // RDKit❗❌:           break;
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // END RDKIT CPP HELPER swapBondDirIfNeeded
        // BEGIN RDKIT CPP HELPER getUnspecifiedQueryBond
        // RDKit❗❌: RDKit::QueryBond *getUnspecifiedQueryBond(const RDKit::Atom *a1,
        // RDKit❗❌:                                           const RDKit::Atom *a2) {
        // RDKit❗❌:   PRECONDITION(a1, "bad atom pointer");
        // RDKit❗❌:   QueryBond *newB;
        // RDKit❗❌:   if (!a1->getIsAromatic() || (a2 && !a2->getIsAromatic())) {
        // RDKit❗❌:     newB = new QueryBond(Bond::SINGLE);
        // RDKit❗❌:     newB->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗❌:   } else {
        // RDKit❗❌:     newB = new QueryBond(Bond::AROMATIC);
        // RDKit❗❌:     newB->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗❌:   }
        // RDKit❗❌:   newB->setProp(RDKit::common_properties::_unspecifiedOrder, 1);
        // RDKit❗❌:   return newB;
        // RDKit❗❌: }
        // END RDKIT CPP HELPER getUnspecifiedQueryBond
        // BEGIN RDKIT CPP HELPER GetUnspecifiedBondType
        // RDKit❗❌: Bond::BondType GetUnspecifiedBondType(const RWMol *mol, const Atom *atom1,
        // RDKit❗❌:                                       const Atom *atom2) {
        // RDKit❗❌:   PRECONDITION(mol, "no molecule");
        // RDKit❗❌:   PRECONDITION(atom1, "no atom1");
        // RDKit❗❌:   PRECONDITION(atom2, "no atom2");
        // RDKit❗❌:   Bond::BondType res;
        // RDKit❗❌:   if (atom1->getIsAromatic() && atom2->getIsAromatic()) {
        // RDKit❗❌:     res = Bond::AROMATIC;
        // RDKit❗❌:   } else {
        // RDKit❗❌:     res = Bond::SINGLE;
        // RDKit❗❌:   }
        // RDKit❗❌:   return res;
        // RDKit❗❌: }
        // END RDKIT CPP HELPER GetUnspecifiedBondType
        // BEGIN RDKIT CPP HELPER SetUnspecifiedBondTypes
        // RDKit❗❌: void SetUnspecifiedBondTypes(RWMol *mol) {
        // RDKit❗❌:   PRECONDITION(mol, "no molecule");
        // RDKit❗❌:   for (auto bond : mol->bonds()) {
        // RDKit❗❌:     if (bond->hasProp(RDKit::common_properties::_unspecifiedOrder)) {
        // RDKit❗❌:       bond->setBondType(GetUnspecifiedBondType(mol, bond->getBeginAtom(),
        // RDKit❗❌:                                                bond->getEndAtom()));
        // RDKit❗❌:       if (bond->getBondType() == Bond::AROMATIC) {
        // RDKit❗❌:         bond->setIsAromatic(true);
        // RDKit❗❌:       } else {
        // RDKit❗❌:         bond->setIsAromatic(false);
        // RDKit❗❌:       }
        // RDKit❗❌:     }
        // RDKit❗❌:   }
        // RDKit❗❌: }
        // END RDKIT CPP HELPER SetUnspecifiedBondTypes
        // Local complexity review: occurrences remain grouped in a sorted
        // BTreeMap and are paired in source order in one pass. The endpoint
        // BTreeSet makes duplicate lookup O(log E) with an additional O(E)
        // index; source uses molecule adjacency lookup instead. Reconciliation
        // otherwise appends source-sorted closure rows without rescanning tokens.
        let ring_closures = std::mem::take(&mut self.ring_closure_targets);
        for (number, occurrences) in ring_closures {
            let mut occurrences = occurrences.into_iter();
            while let Some(open) = occurrences.next() {
                let Some(close) = occurrences.next() else {
                    return Err(SmartsParseError::Parse("unclosed ring".to_owned()));
                };
                let open_atom_idx = open.atom_idx;
                let close_atom_idx = close.atom_idx;
                if open_atom_idx == close_atom_idx {
                    return Err(SmartsParseError::SelfRingClosure {
                        ring: number,
                        atom: open_atom_idx,
                    });
                }
                if graph.has_bond_between(open_atom_idx, close_atom_idx) {
                    return Err(SmartsParseError::DuplicateRingBond {
                        ring: number,
                        begin_atom: open_atom_idx,
                        end_atom: close_atom_idx,
                    });
                }

                // RDKit✔️✔️: if (bond2->hasProp("_cxsmilesBondIdx")) {
                // RDKit✔️✔️:   bond1->setProp("_cxsmilesBondIdx",
                // RDKit✔️✔️:                bond2->getProp<unsigned int>("_cxsmilesBondIdx"));
                // RDKit✔️✔️: }
                // Both implicit closing occurrences and explicit tokens carry
                // the second index; final ring row order does not renumber it.
                let source_index = close.source_index.or(open.source_index).ok_or_else(|| {
                    SmartsParseError::Parse(
                        "closed SMARTS ring has no source bond index".to_owned(),
                    )
                })?;
                let open_is_unspecified = ring_closure_is_unspecified(&open.bond);
                let selected_is_open = !open_is_unspecified;
                let (mut selected_bond, selected_direction, other_direction) = if selected_is_open {
                    (open.bond, open.direction, close.direction)
                } else {
                    (close.bond, close.direction, open.direction)
                };

                let direction = if selected_direction != BondDirection::None {
                    selected_direction
                } else {
                    match other_direction {
                        BondDirection::EndUpRight => BondDirection::EndDownRight,
                        BondDirection::EndDownRight => BondDirection::EndUpRight,
                        other => other,
                    }
                };

                let (begin, end) = match (selected_bond.carrier_order, selected_is_open) {
                    (BondOrder::DativeRight, true) => (open_atom_idx, close_atom_idx),
                    (BondOrder::DativeLeft, true) => (close_atom_idx, open_atom_idx),
                    (BondOrder::DativeRight, false) => (close_atom_idx, open_atom_idx),
                    (BondOrder::DativeLeft, false) => (open_atom_idx, close_atom_idx),
                    (_, true) => (open_atom_idx, close_atom_idx),
                    (_, false) => (close_atom_idx, open_atom_idx),
                };
                if matches!(
                    selected_bond.carrier_order,
                    BondOrder::DativeRight | BondOrder::DativeLeft
                ) {
                    selected_bond.carrier_order = BondOrder::Dative;
                } else if selected_bond.unspecified_order {
                    // RDKit✔️✔️: bond->setBondType(GetUnspecifiedBondType(mol, atom1, atom2));
                    selected_bond.carrier_order =
                        graph.implicit_bond_order(open_atom_idx, close_atom_idx);
                }
                let bond_index = i32::try_from(graph.bond_queries.len()).map_err(|_| {
                    SmartsParseError::Parse(
                        "ring bond index exceeds the source signed carrier".into(),
                    )
                })?;
                graph.push_bond(begin, end, selected_bond, direction, source_index);
                // RDKit❗✔️: *closurePos = bondIdx - 1;
                for index in [open_atom_idx, close_atom_idx] {
                    let atom = &mut graph.atoms[index];
                    let mut ids = atom
                        .prop("_RingClosures")
                        .and_then(|value| value.as_int_vector().ok())
                        .ok_or_else(|| {
                            SmartsParseError::Parse("missing ring atom-order placeholder".into())
                        })?
                        .to_vec();
                    let position = ids
                        .iter()
                        .position(|&id| id == (number as i32).wrapping_add(1).wrapping_neg())
                        .ok_or_else(|| {
                            SmartsParseError::Parse(
                                "missing unresolved ring atom-order placeholder".into(),
                            )
                        })?;
                    ids[position] = bond_index;
                    atom.set_prop("_RingClosures", ids)
                        .expect("the internal ring-closure property key is non-empty");
                }
            }
        }
        Ok(())
    }

    fn current_bond_direction(&self) -> BondDirection {
        // BEGIN RDKIT CPP FUNCTION SMARTS directional bond token actions
        // RDKit❗❌: [\\]{1,2}    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit❗❌: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
        // RDKit❗❌: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗❌: 	return BOND_TOKEN;  }
        // RDKit❗❌: [\/]    { yylval->bond = new QueryBond(Bond::SINGLE);
        // RDKit❗❌: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
        // RDKit❗❌: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
        // RDKit❗❌: 	return BOND_TOKEN;  }
        // END RDKIT CPP FUNCTION SMARTS directional bond token actions
        // The source selects direction on the left QueryBond. This adapter
        // scans the remaining token slice for that first direction token.
        self.tokens[self.pos..]
            .iter()
            .find_map(|(token, _)| match token {
                Token::BondSpec(BondLexeme::Symbol(b'/')) => Some(BondDirection::EndUpRight),
                Token::BondSpec(BondLexeme::Symbol(b'\\')) => Some(BondDirection::EndDownRight),
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
        // RDKit✔️✔️: atomd:	simple_atom
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
        let (token, pos) = self.peek().clone();
        match token {
            Token::OrganicElement(name) | Token::AromaticElement(name) => {
                let atom = parse_simple_atom(&name).ok_or_else(|| {
                    SmartsParseError::InvalidAtomPrimitive {
                        position: self.input_character_position(),
                        detail: format!("invalid simple atom {name:?}"),
                    }
                })?;
                self.advance();
                Ok(atom.into_parsed_atom(None))
            }
            Token::SimpleAtomQuery(name) => {
                let name = [name];
                let atom = parse_simple_atom(&name).ok_or_else(|| {
                    SmartsParseError::InvalidAtomPrimitive {
                        position: self.input_character_position(),
                        detail: format!("invalid simple atom query {name:?}"),
                    }
                })?;
                self.advance();
                Ok(atom.into_parsed_atom(None))
            }
            Token::BracketContent(content) => {
                self.advance();
                self.parse_bracket_atom_content(&content)
            }
            Token::EndOfStream => Err(SmartsParseError::UnexpectedEnd(
                "expected atom but reached end".to_string(),
            )),
            Token::BadCharacter(character) => Err(Self::bad_character_error(pos, character)),
            _ => {
                // RDKit❗✔️:   yyerror(input, molList, current_token_position, "syntax error");
                Err(self.unexpected_token_error(self.peek().1, "expected atom expression"))
            }
        }
    }

    /// Parse bracket atom content in source grammar order, for example
    /// `C@@H`, `N+`, `O-`, `#6`, or `6X4`.
    fn parse_bracket_atom_content(
        &mut self,
        content: &BracketContent,
    ) -> Result<ParsedSmartsAtom, SmartsParseError> {
        // Atom carrier state and predicate trees are reduced together by
        // `ParsedAtomExpr`; the carrier is not reconstructed from the result.
        // One QueryAtom carrier and its predicate move together through each
        // source precedence reduction; no final identity is inferred from the
        // predicate. Local complexity review: content is scanned once in
        // O(n), while query nodes and common carrier state move in source order.
        let (content_text, atom_map) = split_atom_map_suffix(&content.text)?;
        let chars = content_text;
        let len = chars.len();
        if len == 0 {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: 0,
                detail: "empty atom expression".to_string(),
            });
        }
        if let Some(atom) = self.try_parse_hydrogen_atom(&chars, len)? {
            return Ok(atom.into_parsed_atom(atom_map));
        }
        let mut i = 0;
        let mut needs_operand = true;
        let mut clauses: Vec<ParsedAtomExpr> = Vec::new();
        let mut current_or_terms: Vec<ParsedAtomExpr> = Vec::new();
        let mut current_term: Vec<ParsedAtomExpr> = Vec::new();
        let mut lexical_index = 0usize;

        fn finalize_term(
            current_term: &mut Vec<ParsedAtomExpr>,
            current_or_terms: &mut Vec<ParsedAtomExpr>,
        ) {
            if current_term.is_empty() {
                return;
            }
            let mut terms = std::mem::take(current_term).into_iter();
            let mut term = terms.next().expect("nonempty atom-query term");
            for point_query in terms {
                term = term.and_point_query(point_query);
            }
            current_or_terms.push(term);
        }

        fn finalize_clause(
            current_term: &mut Vec<ParsedAtomExpr>,
            current_or_terms: &mut Vec<ParsedAtomExpr>,
            clauses: &mut Vec<ParsedAtomExpr>,
        ) {
            finalize_term(current_term, current_or_terms);
            if current_or_terms.is_empty() {
                return;
            }
            let mut terms = std::mem::take(current_or_terms).into_iter();
            let mut clause = terms.next().expect("nonempty atom-query clause");
            for expression in terms {
                clause = clause.reduce_atom_expr(expression, CompositeQueryType::Or);
            }
            clauses.push(clause);
        }

        while i < len {
            let ch = chars[i];

            // Handle logical OR (comma)
            if ch == b',' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                finalize_term(&mut current_term, &mut current_or_terms);
                needs_operand = true;
                i += 1;
                continue;
            }

            // Handle logical AND (ampersand)
            if ch == b'&' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                needs_operand = true;
                i += 1;
                continue;
            }

            // Handle semicolon (AND)
            if ch == b';' {
                if needs_operand {
                    return Err(invalid_atom_operator(i, ch));
                }
                finalize_clause(&mut current_term, &mut current_or_terms, &mut clauses);
                needs_operand = true;
                i += 1;
                continue;
            }

            let (point_query, consumed) = self.parse_point_query(
                &chars,
                i,
                len,
                &content.lexical_tokens,
                &mut lexical_index,
            )?;

            current_term.push(point_query);
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

        // Combine clauses
        // RDKit source: smarts.yy precedence gives implicit/`&` high-precedence
        // AND inside each comma term, comma OR inside a clause, and `;`
        // low-precedence AND across clauses. Every reduction uses the canonical
        // QueryAtom::expandQuery port so AtomNull algebra is preserved.
        let mut clauses = clauses.into_iter();
        let mut expression = clauses.next().expect("at least one bracket clause");
        for clause in clauses {
            expression = expression.reduce_atom_expr(clause, CompositeQueryType::And);
        }
        Ok(expression.into_parsed_atom(atom_map))
    }

    fn parse_point_query(
        &self,
        chars: &[u8],
        start: usize,
        len: usize,
        lexical_tokens: &[ScannedToken],
        lexical_index: &mut usize,
    ) -> Result<(ParsedAtomExpr, usize), SmartsParseError> {
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
        // Local complexity review: each leading NOT is consumed once, the
        // selected recursive/atom query is parsed once, and carrier effects
        // are applied per source reduction. This is O(n) time with O(1)
        // auxiliary state; query negation is reduced by its source parity.
        let mut pos = start;
        while lexical_tokens
            .get(*lexical_index)
            .is_some_and(|token| token.span.input_char_start < pos)
        {
            *lexical_index += 1;
        }
        let mut negation_count = 0usize;
        let mut negate = false;
        while chars.get(pos) == Some(&b'!') {
            negation_count += 1;
            negate = !negate;
            pos += 1;
        }
        if pos >= len {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: "NOT has no point query".to_string(),
            });
        }
        while lexical_tokens
            .get(*lexical_index)
            .is_some_and(|token| token.span.input_char_start < pos)
        {
            *lexical_index += 1;
        }
        let lexical_token = lexical_tokens
            .get(*lexical_index)
            .filter(|token| token.span.input_char_start == pos);
        let (query, consumed) = self.parse_atom_primitive(chars, pos, len, lexical_token)?;
        let mut atom = self.apply_source_atom_carrier(
            chars,
            pos,
            consumed,
            lexical_token,
            ParsedAtomExpr::query_only(query),
        )?;
        for action_index in 0..negation_count {
            if action_index == 0 && negate {
                let negated = atom.carrier.predicate().is_negated();
                atom.carrier.predicate_mut().set_negation(!negated);
            }
            atom = atom.reset_atomic_number();
            atom.clear_chemical_properties();
        }
        while lexical_tokens
            .get(*lexical_index)
            .is_some_and(|token| token.span.input_char_start < consumed)
        {
            *lexical_index += 1;
        }
        Ok((atom, consumed))
    }

    fn apply_source_atom_carrier(
        &self,
        chars: &[u8],
        start: usize,
        end: usize,
        lexical_token: Option<&ScannedToken>,
        mut atom: ParsedAtomExpr,
    ) -> Result<ParsedAtomExpr, SmartsParseError> {
        // BEGIN RDKIT CPP FUNCTION atom_query carrier actions
        // RDKit❗❌: | number simple_atom {
        // RDKit❗❌:   $2->setIsotope($1);
        // RDKit❗❌:   $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit❗❌:   $$=$2;
        // RDKit❗❌: }
        // RDKit❗❌: | number ATOM_TOKEN {
        // RDKit❗❌:   $2->setIsotope($1);
        // RDKit❗❌:   $2->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit❗❌:   $$=$2;
        // RDKit❗❌: }
        // RDKit❗❌: | HASH_TOKEN number { $$ = new QueryAtom($2); }
        // RDKit❗❌: | number HASH_TOKEN number {
        // RDKit❗❌:   $$ = new QueryAtom($3);
        // RDKit❗❌:   $$->setIsotope($1);
        // RDKit❗❌:   $$->expandQuery(makeAtomIsotopeQuery($1),Queries::COMPOSITE_AND,true);
        // RDKit❗❌: }
        // RDKit❗❌: | number H_TOKEN {
        // RDKit❗❌:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗❌:   newQ->setQuery(makeAtomIsotopeQuery($1));
        // RDKit❗❌:   newQ->setIsotope($1);
        // RDKit❗❌:   newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
        // RDKit❗❌:   newQ->setNumExplicitHs(1);
        // RDKit❗❌:   newQ->setNoImplicit(true);
        // RDKit❗❌:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗❌:   $$=newQ;
        // RDKit❗❌: }
        // RDKit❗❌: | number H_TOKEN number {
        // RDKit❗❌:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗❌:   newQ->setQuery(makeAtomIsotopeQuery($1));
        // RDKit❗❌:   newQ->setIsotope($1);
        // RDKit❗❌:   newQ->expandQuery(makeAtomHCountQuery($3),Queries::COMPOSITE_AND,true);
        // RDKit❗❌:   newQ->setNumExplicitHs($3);
        // RDKit❗❌:   newQ->setNoImplicit(true);
        // RDKit❗❌:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗❌:   $$=newQ;
        // RDKit❗❌: }
        // RDKit❗❌: | H_TOKEN number {
        // RDKit❗❌:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗❌:   newQ->setQuery(makeAtomHCountQuery($2));
        // RDKit❗❌:   newQ->setNumExplicitHs($2);
        // RDKit❗❌:   newQ->setNoImplicit(true);
        // RDKit❗❌:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗❌:   $$=newQ;
        // RDKit❗❌: }
        // RDKit❗❌: | H_TOKEN {
        // RDKit❗❌:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗❌:   newQ->setQuery(makeAtomHCountQuery(1));
        // RDKit❗❌:   newQ->setNumExplicitHs(1);
        // RDKit❗❌:   newQ->setNoImplicit(true);
        // RDKit❗❌:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗❌:   $$=newQ;
        // RDKit❗❌: }
        // RDKit❗❌: | charge_spec {
        // RDKit❗❌:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗❌:   newQ->setQuery(makeAtomFormalChargeQuery($1));
        // RDKit❗❌:   newQ->setFormalCharge($1);
        // RDKit❗❌:   newQ->getFlags() |= SMARTS_CHARGE_MASK;
        // RDKit❗❌:   $$=newQ;
        // RDKit❗❌: }
        // END RDKIT CPP FUNCTION atom_query carrier actions
        // The query leaf is built in `parse_atom_primitive`; this bounded
        // character pass applies the carrier writes from the corresponding
        // source reductions. It does not inspect the resulting predicate tree.
        if start >= end {
            return Ok(atom);
        }
        if let Some(lexical_token) = lexical_token {
            match &lexical_token.token {
                ScannerToken::SimpleAtomQuery(b'a') => {
                    atom.carrier.set_aromatic(true);
                    return Ok(atom);
                }
                ScannerToken::AromaticElement(name) => {
                    let simple = parse_simple_atom(name).ok_or_else(|| {
                        SmartsParseError::InvalidAtomPrimitive {
                            position: start,
                            detail: format!("invalid aromatic element token '{name}'"),
                        }
                    })?;
                    atom = set_atom_carrier_identity(
                        atom,
                        u32::from(simple.atomic_number),
                        simple.aromatic,
                        lexical_token.span.input_char_start,
                    )?;
                    return Ok(atom);
                }
                _ => {}
            }
        }
        let ch = chars[start];
        if ch == b'@' {
            // RDKit❗✔️: | AT_TOKEN AT_TOKEN {
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomNullQuery());
            // RDKit❗✔️:   newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CW);
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // RDKit❗✔️: | AT_TOKEN {
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomNullQuery());
            // RDKit❗✔️:   newQ->setChiralTag(Atom::CHI_TETRAHEDRAL_CCW);
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // RDKit❗✔️: | CHI_CLASS_TOKEN {
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomNullQuery());
            // RDKit❗✔️:   newQ->setChiralTag($1);
            // RDKit❗✔️:   newQ->setProp(common_properties::_chiralPermutation,0);
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // RDKit❗✔️: | CHI_CLASS_TOKEN number {
            // RDKit❗✔️:   if($2==0){
            // RDKit❗✔️:     yyerror(input,molList,branchPoints,scanner,start_token, current_token_position,
            // RDKit❗✔️:             "chiral permutation cannot be zero");
            // RDKit❗✔️:     yyErrorCleanup(molList);
            // RDKit❗✔️:     YYABORT;
            // RDKit❗✔️:   }
            // RDKit❗✔️:
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomNullQuery());
            // RDKit❗✔️:   newQ->setChiralTag($1);
            // RDKit❗✔️:   newQ->setProp(common_properties::_chiralPermutation,$2);
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // RDKit❗✔️: | HYB_TOKEN
            // RDKit❗✔️: | number {
            // Only one bounded grammar primitive is decoded here, before
            // negation and atom_expr reductions. Carrier state stays separate
            // from predicate algebra, as in the source QueryAtom actions.
            let query = std::mem::replace(atom.carrier.predicate_mut(), make_atom_null_query());
            match query {
                QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(tag)) => {
                    atom.carrier.set_chiral_tag(tag);
                }
                QueryNode::And(children) => {
                    let mut children = children.into_iter();
                    match (children.next(), children.next(), children.next()) {
                        (
                            Some(QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(tag))),
                            Some(QueryNode::Predicate(AtomQueryPredicate::ChiralPermutationMatch(
                                permutation,
                            ))),
                            None,
                        ) => {
                            atom.carrier.set_chiral_tag(tag);
                            atom.carrier.set_chiral_permutation(Some(permutation));
                        }
                        _ => {
                            return Err(SmartsParseError::Parse(
                                "invalid chiral point-query carrier".into(),
                            ));
                        }
                    }
                }
                _ => {
                    return Err(SmartsParseError::Parse(
                        "invalid chiral point-query primitive".into(),
                    ));
                }
            }
            return Ok(atom);
        }
        if ch == b'#' {
            let (number, _) = self.parse_number(chars, start + 1, end)?;
            atom = set_atom_carrier_identity(atom, number, false, start)?;
            return Ok(atom);
        }
        // RDKit source: third_party/rdkit/Code/GraphMol/SmilesParse/smarts.ll
        // RDKit❗❌: <IN_ATOM_STATE>He |
        // RDKit❗❌: <IN_ATOM_STATE>Ho |
        // RDKit❗❌: <IN_ATOM_STATE>Hf |
        // RDKit❗❌: <IN_ATOM_STATE>Hg |
        // RDKit❗❌: <IN_ATOM_STATE>Hs |
        // RDKit❗❌: H			{  return H_TOKEN;  }
        // These full element tokens win by length before the single H token.
        // Check only the bounded two-character candidate, then retain the H
        // count carrier action below when that candidate is not an element.
        if ch == b'H' && start + 2 <= end {
            let name = &chars[start..start + 2];
            if let Some(simple) = parse_atom_token(&name) {
                atom = set_atom_carrier_identity(
                    atom,
                    u32::from(simple.atomic_number),
                    simple.aromatic,
                    start,
                )?;
                return Ok(atom);
            }
        }
        if matches!(ch, b'+' | b'-') {
            if let Some((charge, _)) = self.parse_charge_spec(chars, start, end)? {
                // RDKit's `setFormalCharge(int)` writes into its signed
                // `int8_t` carrier while the query leaf retains the `int`.
                // The pinned signed-char projection wraps at this boundary.
                atom.carrier.set_formal_charge(charge as i8);
                atom.charge_mask = true;
            }
            return Ok(atom);
        }
        if ch == b'H' {
            let (count, consumed) = self.parse_optional_number(chars, start + 1, end)?;
            let count = if consumed == start + 1 {
                1
            } else {
                count as u8
            };
            // RDKit❗✔️: void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
            // RDKit❗✔️: std::uint8_t d_numExplicitHs;
            atom.carrier.set_explicit_hydrogens(count);
            atom.carrier.set_no_implicit(true);
            atom.hydrogen_mask = true;
            return Ok(atom);
        }
        if ch.is_ascii_digit() {
            let (number, consumed) = self.parse_number(chars, start, end)?;
            // RDKit's unsigned-int isotope assignment narrows into its
            // uint16_t carrier; this carrier projection is independent of the
            // full signed-int predicate target constructed by the parser.
            let isotope = number as u16;
            if chars.get(consumed) == Some(&b'H') {
                let (count, count_end) = self.parse_optional_number(chars, consumed + 1, end)?;
                let count = if count_end == consumed + 1 {
                    1
                } else {
                    count as u8
                };
                // RDKit❗✔️: void setNumExplicitHs(unsigned int what) { d_numExplicitHs = what; }
                // RDKit❗✔️: std::uint8_t d_numExplicitHs;
                atom.carrier.set_isotope(Some(isotope));
                atom.carrier.set_explicit_hydrogens(count);
                atom.carrier.set_no_implicit(true);
                atom.hydrogen_mask = true;
                return Ok(atom);
            }
            if chars.get(consumed) == Some(&b'#') {
                let (atomic_number, _) = self.parse_number(chars, consumed + 1, end)?;
                atom = set_atom_carrier_identity(atom, atomic_number, false, consumed + 1)?;
                atom.carrier.set_isotope(Some(isotope));
                return Ok(atom);
            }
            let symbol_start = consumed;
            if symbol_start < end {
                let mut symbol_end = symbol_start + 1;
                if chars.get(symbol_end).is_some_and(u8::is_ascii_lowercase) {
                    let name = &chars[symbol_start..=symbol_end];
                    if parse_atom_token(&name).is_some() {
                        symbol_end += 1;
                    }
                }
                let name = &chars[symbol_start..symbol_end];
                if let Some(simple) = parse_atom_token(&name) {
                    atom = set_atom_carrier_identity(
                        atom,
                        u32::from(simple.atomic_number),
                        simple.aromatic,
                        symbol_start,
                    )?;
                    atom.carrier.set_isotope(Some(isotope));
                }
            }
            return Ok(atom);
        }

        if ch.is_ascii_uppercase() && ch != b'H' {
            let mut symbol_end = start + 1;
            if start + 3 <= end {
                let name = &chars[start..start + 3];
                if parse_atom_token(&name).is_some() {
                    symbol_end = start + 3;
                }
            }
            if chars.get(symbol_end).is_some_and(u8::is_ascii_lowercase) {
                if symbol_end == start + 1 {
                    let name = &chars[start..=symbol_end];
                    if parse_atom_token(&name).is_some() {
                        symbol_end += 1;
                    }
                }
            }
            let name = &chars[start..symbol_end];
            if let Some(simple) = parse_atom_token(&name) {
                atom = set_atom_carrier_identity(
                    atom,
                    u32::from(simple.atomic_number),
                    simple.aromatic,
                    start,
                )?;
            }
            return Ok(atom);
        }
        Ok(atom)
    }

    fn try_parse_hydrogen_atom(
        &self,
        chars: &[u8],
        len: usize,
    ) -> Result<Option<ParsedAtomExpr>, SmartsParseError> {
        // RDKit✔️✔️: hydrogen_atom:	ATOM_OPEN_TOKEN H_TOKEN ATOM_CLOSE_TOKEN
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
            isotope = Some(
                i32::try_from(num)
                    .expect("SMARTS number is bounded to the source nonnegative int32 range"),
            );
            pos = consumed;
        }
        if pos >= len || chars[pos] != b'H' {
            return Ok(None);
        }
        pos += 1;
        if pos < len && chars[pos].is_ascii_digit() {
            return Ok(None);
        }

        let mut formal_charge = None;
        if pos < len && matches!(chars[pos], b'+' | b'-') {
            let Some((charge, consumed)) = self.parse_charge_spec(chars, pos, len)? else {
                return Ok(None);
            };
            formal_charge = Some(charge);
            pos = consumed;
        }

        if pos != len {
            return Ok(None);
        }

        let mut atom = ParsedAtomExpr::with_identity(
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(1)),
            1,
            false,
        );
        if let Some(isotope) = isotope {
            atom.carrier.set_isotope(Some(isotope as u16));
            crate::query_behavior::query_atom_expand_query(
                atom.carrier.predicate_mut(),
                crate::query_behavior::make_atom_isotope_query(isotope),
                CompositeQueryType::And,
                true,
            );
        }
        if let Some(formal_charge) = formal_charge {
            // Keep the source int query value separate from Atom's int8_t
            // carrier projection.
            atom.carrier.set_formal_charge(formal_charge as i8);
            atom.charge_mask = true;
            crate::query_behavior::query_atom_expand_query(
                atom.carrier.predicate_mut(),
                crate::query_behavior::make_atom_formal_charge_query(formal_charge),
                CompositeQueryType::And,
                true,
            );
        }
        Ok(Some(atom))
    }

    /// Parse one currently modeled atom primitive from bracket content at `i`.
    ///
    /// Grammar anchors below identify the source productions represented here.
    fn parse_atom_primitive(
        &self,
        chars: &[u8],
        i: usize,
        len: usize,
        lexical_token: Option<&ScannedToken>,
    ) -> Result<(QueryNode<AtomQueryPredicate>, usize), SmartsParseError> {
        // RDKit✔️✔️: atom_query:	simple_atom
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
        // RDKit❗✔️: | HASH_TOKEN number { $$ = new QueryAtom($2); }
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

        if let Some(lexical_token) = lexical_token {
            match &lexical_token.token {
                ScannerToken::SimpleAtomQuery(b'a') => {
                    let atom = parse_simple_atom("a").expect("generic aromatic query token");
                    return Ok((atom.query, lexical_token.span.input_char_end));
                }
                ScannerToken::AromaticElement(name) => {
                    let atom = parse_simple_atom(name).ok_or_else(|| {
                        SmartsParseError::InvalidAtomPrimitive {
                            position: i,
                            detail: format!("invalid aromatic simple atom '{name}'"),
                        }
                    })?;
                    return Ok((atom.query, lexical_token.span.input_char_end));
                }
                _ => {}
            }
        }

        // Atomic number: #N
        // RDKit❗✔️: | HASH_TOKEN number { $$ = new QueryAtom($2); }
        // RDKit❗✔️: ATOM_EQUALS_QUERY *makeAtomNumQuery(int what) {
        // RDKit❗✔️:   return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNum,
        // RDKit❗✔️:                                                 "AtomAtomicNum");
        // RDKit❗✔️: }
        if ch == b'#' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            let atomic_number =
                u8::try_from(num).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
                    position: i + 1,
                    detail: "atomic number is outside the parser carrier's u8 range".to_string(),
                })?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
                consumed,
            ));
        }

        // Recursive SMARTS: $(...)
        if ch == b'$' {
            return self.parse_recursive_query(chars, i, len);
        }

        if let Some((charge, consumed)) = self.parse_charge_spec(chars, i, len)? {
            return Ok((
                crate::query_behavior::make_atom_formal_charge_query(charge),
                consumed,
            ));
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
        if ch == b'@' {
            let start = i + 1;
            if start < len && chars[start] == b'@' {
                // @@
                return Ok((
                    QueryNode::Predicate(AtomQueryPredicate::ChiralTagMatch(
                        ChiralTag::TetrahedralCw,
                    )),
                    start + 1,
                ));
            }
            let class = if start + 1 < len {
                Some(&chars[start..=start + 1])
            } else {
                None
            };
            if let Some((tag, class_end)) = class.and_then(|class| {
                let tag = match class {
                    b"TH" => ChiralTag::Tetrahedral,
                    b"AL" => ChiralTag::Allene,
                    b"SP" => ChiralTag::SquarePlanar,
                    b"TB" => ChiralTag::TrigonalBipyramidal,
                    b"OH" => ChiralTag::Octahedral,
                    _ => return None,
                };
                Some((tag, start + 2))
            }) {
                let (permutation, consumed) = self.parse_optional_number(chars, class_end, len)?;
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
                    ChiralTag::TetrahedralCcw,
                )),
                start,
            ));
        }

        // Element symbol
        // The lexer source below contains the full element-symbol union. This
        // character path applies longest-token matching to one- and two-letter
        // symbols and the source's three-letter temporary symbols.
        // RDKit❗❌: <IN_ATOM_STATE>He |
        // RDKit❗❌: <IN_ATOM_STATE>Li |
        // RDKit❗❌: <IN_ATOM_STATE>Be |
        // RDKit❗❌: <IN_ATOM_STATE>Ne |
        // RDKit❗❌: <IN_ATOM_STATE>Na |
        // RDKit❗❌: <IN_ATOM_STATE>Mg |
        // RDKit❗❌: <IN_ATOM_STATE>Al |
        // RDKit❗❌: <IN_ATOM_STATE>Si |
        // RDKit❗❌: <IN_ATOM_STATE>Ar |
        // RDKit❗❌: <IN_ATOM_STATE>K |
        // RDKit❗❌: <IN_ATOM_STATE>Ca |
        // RDKit❗❌: <IN_ATOM_STATE>Sc |
        // RDKit❗❌: <IN_ATOM_STATE>Ti |
        // RDKit❗❌: <IN_ATOM_STATE>V |
        // RDKit❗❌: <IN_ATOM_STATE>Cr |
        // RDKit❗❌: <IN_ATOM_STATE>Mn |
        // RDKit❗❌: <IN_ATOM_STATE>Co |
        // RDKit❗❌: <IN_ATOM_STATE>Fe |
        // RDKit❗❌: <IN_ATOM_STATE>Ni |
        // RDKit❗❌: <IN_ATOM_STATE>Cu |
        // RDKit❗❌: <IN_ATOM_STATE>Zn |
        // RDKit❗❌: <IN_ATOM_STATE>Ga |
        // RDKit❗❌: <IN_ATOM_STATE>Ge |
        // RDKit❗❌: <IN_ATOM_STATE>As |
        // RDKit❗❌: <IN_ATOM_STATE>Se |
        // RDKit❗❌: <IN_ATOM_STATE>Kr |
        // RDKit❗❌: <IN_ATOM_STATE>Rb |
        // RDKit❗❌: <IN_ATOM_STATE>Sr |
        // RDKit❗❌: <IN_ATOM_STATE>Y |
        // RDKit❗❌: <IN_ATOM_STATE>Zr |
        // RDKit❗❌: <IN_ATOM_STATE>Nb |
        // RDKit❗❌: <IN_ATOM_STATE>Mo |
        // RDKit❗❌: <IN_ATOM_STATE>Tc |
        // RDKit❗❌: <IN_ATOM_STATE>Ru |
        // RDKit❗❌: <IN_ATOM_STATE>Rh |
        // RDKit❗❌: <IN_ATOM_STATE>Pd |
        // RDKit❗❌: <IN_ATOM_STATE>Ag |
        // RDKit❗❌: <IN_ATOM_STATE>Cd |
        // RDKit❗❌: <IN_ATOM_STATE>In |
        // RDKit❗❌: <IN_ATOM_STATE>Sn |
        // RDKit❗❌: <IN_ATOM_STATE>Sb |
        // RDKit❗❌: <IN_ATOM_STATE>Te |
        // RDKit❗❌: <IN_ATOM_STATE>Xe |
        // RDKit❗❌: <IN_ATOM_STATE>Cs |
        // RDKit❗❌: <IN_ATOM_STATE>Ba |
        // RDKit❗❌: <IN_ATOM_STATE>La |
        // RDKit❗❌: <IN_ATOM_STATE>Ce |
        // RDKit❗❌: <IN_ATOM_STATE>Pr |
        // RDKit❗❌: <IN_ATOM_STATE>Nd |
        // RDKit❗❌: <IN_ATOM_STATE>Pm |
        // RDKit❗❌: <IN_ATOM_STATE>Sm |
        // RDKit❗❌: <IN_ATOM_STATE>Eu |
        // RDKit❗❌: <IN_ATOM_STATE>Gd |
        // RDKit❗❌: <IN_ATOM_STATE>Tb |
        // RDKit❗❌: <IN_ATOM_STATE>Dy |
        // RDKit❗❌: <IN_ATOM_STATE>Ho |
        // RDKit❗❌: <IN_ATOM_STATE>Er |
        // RDKit❗❌: <IN_ATOM_STATE>Tm |
        // RDKit❗❌: <IN_ATOM_STATE>Yb |
        // RDKit❗❌: <IN_ATOM_STATE>Lu |
        // RDKit❗❌: <IN_ATOM_STATE>Hf |
        // RDKit❗❌: <IN_ATOM_STATE>Ta |
        // RDKit❗❌: <IN_ATOM_STATE>W |
        // RDKit❗❌: <IN_ATOM_STATE>Re |
        // RDKit❗❌: <IN_ATOM_STATE>Os |
        // RDKit❗❌: <IN_ATOM_STATE>Ir |
        // RDKit❗❌: <IN_ATOM_STATE>Pt |
        // RDKit❗❌: <IN_ATOM_STATE>Au |
        // RDKit❗❌: <IN_ATOM_STATE>Hg |
        // RDKit❗❌: <IN_ATOM_STATE>Tl |
        // RDKit❗❌: <IN_ATOM_STATE>Pb |
        // RDKit❗❌: <IN_ATOM_STATE>Bi |
        // RDKit❗❌: <IN_ATOM_STATE>Po |
        // RDKit❗❌: <IN_ATOM_STATE>At |
        // RDKit❗❌: <IN_ATOM_STATE>Rn |
        // RDKit❗❌: <IN_ATOM_STATE>Fr |
        // RDKit❗❌: <IN_ATOM_STATE>Ra |
        // RDKit❗❌: <IN_ATOM_STATE>Ac |
        // RDKit❗❌: <IN_ATOM_STATE>Th |
        // RDKit❗❌: <IN_ATOM_STATE>Pa |
        // RDKit❗❌: <IN_ATOM_STATE>U |
        // RDKit❗❌: <IN_ATOM_STATE>Np |
        // RDKit❗❌: <IN_ATOM_STATE>Pu |
        // RDKit❗❌: <IN_ATOM_STATE>Am |
        // RDKit❗❌: <IN_ATOM_STATE>Cm |
        // RDKit❗❌: <IN_ATOM_STATE>Bk |
        // RDKit❗❌: <IN_ATOM_STATE>Cf |
        // RDKit❗❌: <IN_ATOM_STATE>Es |
        // RDKit❗❌: <IN_ATOM_STATE>Fm |
        // RDKit❗❌: <IN_ATOM_STATE>Md |
        // RDKit❗❌: <IN_ATOM_STATE>No |
        // RDKit❗❌: <IN_ATOM_STATE>Lr |
        // RDKit❗❌: <IN_ATOM_STATE>Rf |
        // RDKit❗❌: <IN_ATOM_STATE>Db |
        // RDKit❗❌: <IN_ATOM_STATE>Sg |
        // RDKit❗❌: <IN_ATOM_STATE>Bh |
        // RDKit❗❌: <IN_ATOM_STATE>Hs |
        // RDKit❗❌: <IN_ATOM_STATE>Mt |
        // RDKit❗❌: <IN_ATOM_STATE>Ds |
        // RDKit❗❌: <IN_ATOM_STATE>Rg |
        // RDKit❗❌: <IN_ATOM_STATE>Cn |
        // RDKit❗❌: <IN_ATOM_STATE>Uut |
        // RDKit❗❌: <IN_ATOM_STATE>Fl |
        // RDKit❗❌: <IN_ATOM_STATE>Uup |
        // RDKit❗❌: <IN_ATOM_STATE>Lv	{   yylval->atom = new QueryAtom( PeriodicTable::getTable()->getAtomicNumber( yytext ) );
        // RDKit❗❌: 				return ATOM_TOKEN;
        // RDKit❗❌: 			}
        // Local complexity review: this path creates temporary String values
        // before fixed-symbol lookup; the source lexer supplies its matched token
        // and atomic number directly, so this adapter has materially more allocation.
        // RDKit flex selects the longest matching token in IN_ATOM_STATE, so
        // two-letter elements such as Hg must be consumed before the H_TOKEN
        // hydrogen-count rule below.
        if ch.is_ascii_uppercase() {
            let start = i;
            let three_end = start + 3;
            if three_end <= len {
                let three_char = &chars[start..three_end];
                if let Some(atom) = parse_atom_token(&three_char) {
                    return Ok((atom.query, three_end));
                }
            }
            let end = i + 1;
            if end < len && chars[end].is_ascii_lowercase() {
                let two_char = &chars[start..=end];
                if let Some(atom) = parse_atom_token(&two_char) {
                    return Ok((atom.query, end + 1));
                }
            }
            if ch != b'H' {
                let one_char = &chars[start..end];
                if let Some(atom) = parse_atom_token(&one_char) {
                    return Ok((atom.query, end));
                }
            }
        }

        // Hydrogen-count SMARTS queries: `h` or `h<N>`, `H` or `H<N>`
        // RDKit❗✔️: <IN_ATOM_STATE>h {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️:         yylval->atom->setQuery(makeAtomHasImplicitHQuery());
        // RDKit❗✔️: 	return IMPLICIT_H_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: | IMPLICIT_H_ATOM_QUERY_TOKEN number {
        // RDKit❗✔️:   $1->setQuery(makeAtomImplicitHCountQuery($2));
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        if ch == b'h' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((
                    crate::query_behavior::make_atom_has_implicit_h_query(),
                    consumed,
                ));
            }
            return Ok((
                crate::query_behavior::make_atom_implicit_h_count_query(num as i32),
                consumed,
            ));
        }
        // RDKit❗✔️: | H_TOKEN number {
        // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗✔️:   newQ->setQuery(makeAtomHCountQuery($2));
        // RDKit❗✔️:   newQ->setNumExplicitHs($2);
        // RDKit❗✔️:   newQ->setNoImplicit(true);
        // RDKit❗✔️:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗✔️:   $$=newQ;
        // RDKit❗✔️: }
        // RDKit❗✔️: | H_TOKEN {
        // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
        // RDKit❗✔️:   newQ->setQuery(makeAtomHCountQuery(1));
        // RDKit❗✔️:   newQ->setNumExplicitHs(1);
        // RDKit❗✔️:   newQ->setNoImplicit(true);
        // RDKit❗✔️:   newQ->getFlags() |= SMARTS_H_MASK;
        // RDKit❗✔️:   $$=newQ;
        // RDKit❗✔️: }
        if ch == b'H' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((crate::query_behavior::make_atom_h_count_query(1), consumed));
            }
            return Ok((
                crate::query_behavior::make_atom_h_count_query(num as i32),
                consumed,
            ));
        }

        // Ring membership: R or R<N>
        if ch == b'R' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((crate::query_behavior::make_atom_ring_query(-1), consumed));
            }
            if num == 0 {
                // RDKit✔️✔️: <IN_ATOM_STATE>R {
                // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
                // RDKit✔️✔️: 	yylval->atom->setQuery(new AtomRingQuery(-1));
                // RDKit✔️✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
                // RDKit✔️✔️: }
                //
                // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
                // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
                // RDKit✔️✔️:   $$ = $1;
                // RDKit✔️✔️: }
                //
                // Keep the source AtomRingQuery(0) predicate identity.
                return Ok((crate::query_behavior::make_atom_ring_query(0), consumed));
            }
            return Ok((
                // RDKit✔️✔️: <IN_ATOM_STATE>R {
                // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
                // RDKit✔️✔️: 	yylval->atom->setQuery(new AtomRingQuery(-1));
                // RDKit✔️✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
                // RDKit✔️✔️: }
                // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
                // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
                // RDKit✔️✔️:   $$ = $1;
                // RDKit✔️✔️: }
                crate::query_behavior::make_atom_ring_query(num as i32),
                consumed,
            ));
        }
        if ch == b'r' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((crate::query_behavior::make_atom_in_ring_query(), consumed));
            }
            return Ok((
                crate::query_behavior::make_atom_min_ring_size_query(num as i32),
                consumed,
            ));
        }
        if ch == b'k' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((crate::query_behavior::make_atom_in_ring_query(), consumed));
            }
            return Ok((
                crate::query_behavior::make_atom_in_ring_of_size_query(num as i32),
                consumed,
            ));
        }

        // Connectivity/degree: X or X<N>
        // RDKit✔️✔️: <IN_ATOM_STATE>X {
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomTotalDegreeQuery(1));
        // RDKit✔️✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | COMPLEX_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   static_cast<ATOM_EQUALS_QUERY *>($1->getQuery())->setVal($2);
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        if ch == b'X' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                crate::query_behavior::make_atom_total_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num as i32
                }),
                consumed,
            ));
        }

        // Non-hydrogen degree: d or d<N>
        if ch == b'd' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                crate::query_behavior::make_atom_non_hydrogen_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num
                }),
                consumed,
            ));
        }

        // RDKit❗✔️: <IN_ATOM_STATE>z {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHasHeteroatomNbrsQuery());
        // RDKit❗✔️: 	return HETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: <IN_ATOM_STATE>Z {
        // RDKit❗✔️: 	yylval->atom = new QueryAtom();
        // RDKit❗✔️: 	yylval->atom->setQuery(makeAtomHasAliphaticHeteroatomNbrsQuery());
        // RDKit❗✔️: 	return ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN;
        // RDKit❗✔️: }
        // RDKit❗✔️: | HETERONEIGHBOR_ATOM_QUERY_TOKEN number {
        // RDKit❗✔️:   $1->setQuery(makeAtomNumHeteroatomNbrsQuery($2));
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        // RDKit❗✔️: | ALIPHATICHETERONEIGHBOR_ATOM_QUERY_TOKEN number {
        // RDKit❗✔️:   $1->setQuery(makeAtomNumAliphaticHeteroatomNbrsQuery($2));
        // RDKit❗✔️:   $$ = $1;
        // RDKit❗✔️: }
        // Heteroatom-neighbor queries retain distinct source factories for z
        // and Z; parse_optional_number keeps the shared signed-int guard.
        if ch == b'z' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                if consumed == i + 1 {
                    crate::query_behavior::make_atom_has_heteroatom_nbrs_query()
                } else {
                    crate::query_behavior::make_atom_num_heteroatom_nbrs_query(num as i32)
                },
                consumed,
            ));
        }
        if ch == b'Z' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                if consumed == i + 1 {
                    crate::query_behavior::make_atom_has_aliphatic_heteroatom_nbrs_query()
                } else {
                    crate::query_behavior::make_atom_num_aliphatic_heteroatom_nbrs_query(num as i32)
                },
                consumed,
            ));
        }

        // Ring connectivity: x or x<N>
        // RDKit✔️✔️: <IN_ATOM_STATE>x {
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomHasRingBondQuery());
        // RDKit✔️✔️: 	return RINGBOND_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        // RDKit✔️✔️: | RINGBOND_ATOM_QUERY_TOKEN number {
        // RDKit✔️✔️:   $1->setQuery(makeAtomRingBondCountQuery($2));
        // RDKit✔️✔️:   $$ = $1;
        // RDKit✔️✔️: }
        if ch == b'x' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            if consumed == i + 1 {
                return Ok((
                    crate::query_behavior::make_atom_has_ring_bond_query(),
                    consumed,
                ));
            }
            return Ok((
                crate::query_behavior::make_atom_ring_bond_count_query(num as i32),
                consumed,
            ));
        }

        // Degree: D or D<N>
        // RDKit✔️✔️: <IN_ATOM_STATE>D {
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomExplicitDegreeQuery(1));
        // RDKit✔️✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        if ch == b'D' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                crate::query_behavior::make_atom_explicit_degree_query(if consumed == i + 1 {
                    1
                } else {
                    num as i32
                }),
                consumed,
            ));
        }

        // Hybridization: ^1, ^2, ^3, ...
        if ch == b'^' {
            let (num, consumed) = self.parse_number(chars, i + 1, len)?;
            let hybridization = match num {
                0 => Hybridization::S,
                1 => Hybridization::Sp,
                2 => Hybridization::Sp2,
                3 => Hybridization::Sp3,
                4 => Hybridization::Sp3d,
                5 => Hybridization::Sp3d2,
                _ => {
                    return Err(SmartsParseError::InvalidAtomPrimitive {
                        position: i,
                        detail: "hybridization must be in ^0..^5".to_string(),
                    });
                }
            };
            return Ok((
                crate::query_behavior::make_atom_hybridization_query(hybridization),
                consumed,
            ));
        }

        // Valence: v or v<N>
        // RDKit✔️✔️: <IN_ATOM_STATE>v {
        // RDKit✔️✔️: 	yylval->atom = new QueryAtom();
        // RDKit✔️✔️: 	yylval->atom->setQuery(makeAtomTotalValenceQuery(1));
        // RDKit✔️✔️: 	return COMPLEX_ATOM_QUERY_TOKEN;
        // RDKit✔️✔️: }
        if ch == b'v' {
            let (num, consumed) = self.parse_optional_number(chars, i + 1, len)?;
            return Ok((
                QueryNode::Predicate(AtomQueryPredicate::TotalValence(if consumed == i + 1 {
                    1
                } else {
                    num as i32
                })),
                consumed,
            ));
        }

        // The scanner's SIMPLE_ATOM_QUERY_TOKEN branch above identifies `a`;
        // this branch handles the distinct aliphatic `A` query token.
        if ch == b'A' {
            return Ok((
                parse_simple_atom("A")
                    .expect("aliphatic simple query token")
                    .query,
                i + 1,
            ));
        }

        // Unsaturated: u
        if ch == b'u' {
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
            let isotope = i32::try_from(num)
                .expect("SMARTS number is bounded to the source nonnegative int32 range");

            // RDKit❗✔️: | number H_TOKEN {
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomIsotopeQuery($1));
            // RDKit❗✔️:   newQ->setIsotope($1);
            // RDKit❗✔️:   newQ->expandQuery(makeAtomHCountQuery(1),Queries::COMPOSITE_AND,true);
            // RDKit❗✔️:   newQ->setNumExplicitHs(1);
            // RDKit❗✔️:   newQ->setNoImplicit(true);
            // RDKit❗✔️:   newQ->getFlags() |= SMARTS_H_MASK;
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // RDKit❗✔️: | number H_TOKEN number {
            // RDKit❗✔️:   QueryAtom *newQ = new QueryAtom();
            // RDKit❗✔️:   newQ->setQuery(makeAtomIsotopeQuery($1));
            // RDKit❗✔️:   newQ->setIsotope($1);
            // RDKit❗✔️:   newQ->expandQuery(makeAtomHCountQuery($3),Queries::COMPOSITE_AND,true);
            // RDKit❗✔️:   newQ->setNumExplicitHs($3);
            // RDKit❗✔️:   newQ->setNoImplicit(true);
            // RDKit❗✔️:   newQ->getFlags() |= SMARTS_H_MASK;
            // RDKit❗✔️:   $$=newQ;
            // RDKit❗✔️: }
            // Local complexity review: the source number is parsed once per
            // query leaf; the bounded numeric suffix is read again only for
            // the separate carrier write, without allocation or rescanning
            // any unrelated query input.
            if chars.get(consumed) == Some(&b'H') {
                let (hydrogen_count, end) = self.parse_optional_number(chars, consumed + 1, len)?;
                let hydrogen_count = if end == consumed + 1 {
                    1
                } else {
                    hydrogen_count
                };
                let mut query = crate::query_behavior::make_atom_isotope_query(isotope);
                crate::query_behavior::query_atom_expand_query(
                    &mut query,
                    crate::query_behavior::make_atom_h_count_query(hydrogen_count as i32),
                    CompositeQueryType::And,
                    true,
                );
                return Ok((query, end));
            }

            let mut atom_and_end = None;
            if chars.get(consumed) == Some(&b'#') {
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
                if chars.get(end).is_some_and(u8::is_ascii_lowercase) {
                    let two_char = &chars[consumed..=end];
                    if parse_atom_token(&two_char).is_some() {
                        end += 1;
                    }
                }
                let symbol = &chars[consumed..end];
                if *next == b'*' || parse_atom_token(&symbol).is_some() {
                    atom_and_end = parse_atom_token(&symbol).map(|atom| (atom.query, end));
                }
            }
            if let Some((mut atom, end)) = atom_and_end {
                crate::query_behavior::query_atom_expand_query(
                    &mut atom,
                    crate::query_behavior::make_atom_isotope_query(isotope),
                    CompositeQueryType::And,
                    true,
                );
                return Ok((atom, end));
            }
            return Ok((
                crate::query_behavior::make_atom_isotope_query(isotope),
                consumed,
            ));
        }

        // Wildcard inside a bracket uses the same SIMPLE_ATOM_QUERY_TOKEN path.
        if ch == b'*' {
            return Ok((
                parse_simple_atom("*")
                    .expect("wildcard simple query token")
                    .query,
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
        chars: &[u8],
        start: usize,
        len: usize,
    ) -> Result<Option<(i32, usize)>, SmartsParseError> {
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
        if !matches!(sign, b'+' | b'-') {
            return Ok(None);
        }
        let next = start + 1;
        if chars.get(next) == Some(&sign) {
            let charge = if sign == b'+' { 2 } else { -2 };
            return Ok(Some((charge, next + 1)));
        }
        let (magnitude, consumed) = self.parse_optional_number(chars, next, len)?;
        let magnitude = if consumed == next { 1 } else { magnitude };
        let magnitude = i32::try_from(magnitude)
            .expect("parse_optional_number enforces the source nonnegative int32 range");
        let charge = if sign == b'+' { magnitude } else { -magnitude };
        Ok(Some((charge, consumed)))
    }

    fn parse_possible_range_query(
        &self,
        chars: &[u8],
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
            Some(b'D') => AtomRangeDataFunction::ExplicitDegree,
            Some(b'd') => AtomRangeDataFunction::NonHydrogenDegree,
            Some(b'X') => AtomRangeDataFunction::TotalDegree,
            Some(b'v') => AtomRangeDataFunction::TotalValence,
            Some(b'R') => AtomRangeDataFunction::NumAtomRings,
            Some(b'z') => AtomRangeDataFunction::NumHeteroatomNeighbors,
            Some(b'Z') => AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors,
            Some(b'r') => AtomRangeDataFunction::MinRingSize,
            Some(b'x') => AtomRangeDataFunction::RingBondCount,
            Some(b'h') => AtomRangeDataFunction::ImplicitHydrogenCount,
            Some(b'+') => AtomRangeDataFunction::FormalCharge,
            Some(b'-') => AtomRangeDataFunction::NegativeFormalCharge,
            Some(b'k') => AtomRangeDataFunction::AtomRingSize {
                lower: 0,
                upper: 0,
                lower_open: false,
                upper_open: false,
            },
            _ => return Ok(None),
        };
        if chars.get(start + 1) != Some(&b'{') {
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
        let query = (if matches!(
            &data_function,
            AtomRangeDataFunction::NumAtomRings
                | AtomRangeDataFunction::MinRingSize
                | AtomRangeDataFunction::RingBondCount
                | AtomRangeDataFunction::AtomRingSize { .. }
        ) {
            make_atom_possible_ring_range_query(lower, upper, data_function)
        } else {
            make_atom_possible_range_query(lower, upper, data_function)
        })
        .ok_or_else(|| SmartsParseError::InvalidAtomPrimitive {
            position: start,
            detail: "empty atom range".to_string(),
        })?;
        Ok(Some((query, consumed)))
    }

    fn parse_possible_range_bounds(
        &self,
        chars: &[u8],
        start: usize,
        len: usize,
    ) -> Result<(Option<i32>, Option<i32>, usize), SmartsParseError> {
        let mut pos = start;
        let lower = if chars.get(pos) == Some(&b'-') {
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
        if chars.get(pos) != Some(&b'-') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '-' in atom range".to_string(),
            });
        }
        pos += 1;
        let upper = if chars.get(pos).is_some_and(u8::is_ascii_digit) {
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
        if chars.get(pos) != Some(&b'}') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start.saturating_sub(2),
                detail: "expected '}' to close atom range".to_string(),
            });
        }
        Ok((lower, upper, pos + 1))
    }

    fn parse_recursive_query(
        &self,
        chars: &[u8],
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
        if chars.get(start + 1) != Some(&b'(') {
            return Err(SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: "expected '(' after '$'".to_string(),
            });
        }
        let mut depth = 1usize;
        let mut end = start + 2;
        while end < len && depth > 0 {
            match chars[end] {
                b'(' => depth += 1,
                b')' => depth -= 1,
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

        let recursive_smarts = &chars[start..end];
        let mut consumed = end;
        let mut serial_number = 0;
        if chars.get(consumed) == Some(&b'_') {
            consumed += 1;
            let serial_start = consumed;
            if !chars
                .get(consumed)
                .is_some_and(|digit| matches!(*digit, b'1'..=b'9'))
            {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: consumed,
                    detail: "recursive SMARTS serial number must be nonzero".to_string(),
                });
            }
            // The explicit suffix reduces through the same bounded
            // `nonzero_number` rule as every other source integer token.
            (serial_number, consumed) = self.parse_number(chars, serial_start, len)?;
        }
        let inner = recursive_smarts
            .strip_prefix(b"$(")
            .and_then(|value| value.strip_suffix(b")"))
            .expect("balanced recursive SMARTS has delimiters");
        // RDKit✔️✔️: SmilesParseOps::CloseMolRings(molP,0);
        // RDKit✔️✔️: qA->setQuery(new RecursiveStructureQuery(molP));
        // Recursive mol reductions are part of the same grammar invocation.
        // Do not run public preprocessing/CX/cleanup or reset the counter.
        let query_mol = (|| {
            let tokens = tokenize(inner)?;
            let mut parser = SmartsParser::new(&tokens, inner);
            parser.num_bonds_parsed = Rc::clone(&self.num_bonds_parsed);
            parser.parse_smarts_molecule()?.finish_recursive()
        })()
        .map_err(
            |error: SmartsParseError| SmartsParseError::InvalidAtomPrimitive {
                position: start,
                detail: error.to_string(),
            },
        )?;
        Ok((
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(
                crate::query_behavior::RecursiveStructureQuery::from_query_graph(
                    query_mol,
                    serial_number,
                )
                .with_source_smarts(recursive_smarts),
            )),
            consumed,
        ))
    }

    /// Parse a number from position i. Returns (value, consumed_index).
    fn parse_number(
        &self,
        chars: &[u8],
        i: usize,
        len: usize,
    ) -> Result<(u32, usize), SmartsParseError> {
        // RDKit❗✔️: number:  ZERO_TOKEN
        // RDKit❗✔️: | nonzero_number
        // RDKit❗✔️: nonzero_number:  NONZERO_DIGIT_TOKEN
        // RDKit❗✔️: | nonzero_number digit {
        // RDKit❗✔️:     if($1 >= std::numeric_limits<std::int32_t>::max()/10 ||
        // RDKit❗✔️:      $1*10 >= std::numeric_limits<std::int32_t>::max()-$2 ){
        // RDKit❗✔️:      yysmarts_error(input,molList,lastAtom,lastBond,numAtomsParsed,numBondsParsed,branchPoints,scanner,start_token, current_token_position, "number too large");
        // RDKit❗✔️:      YYABORT;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   $$ = $1*10 + $2; }
        // RDKit❗✔️: digit: NONZERO_DIGIT_TOKEN
        // RDKit❗✔️: | ZERO_TOKEN
        // Local complexity review: one left-to-right digit fold is O(n) time
        // and O(1) state; each source reduction is modeled by one bounded
        // guard and update, without reparsing or allocation.
        if i >= len || !chars[i].is_ascii_digit() {
            return Err(SmartsParseError::UnexpectedEnd(
                "expected number".to_string(),
            ));
        }
        if chars[i] == b'0' {
            return Ok((0, i + 1));
        }
        let mut val = u32::from(chars[i] - b'0');
        let mut pos = i + 1;
        while pos < len && chars[pos].is_ascii_digit() {
            let digit = u32::from(chars[pos] - b'0');
            if val >= i32::MAX as u32 / 10 || val * 10 >= i32::MAX as u32 - digit {
                return Err(SmartsParseError::InvalidAtomPrimitive {
                    position: i,
                    detail: "number too large".to_string(),
                });
            }
            val = val * 10 + digit;
            pos += 1;
        }
        Ok((val, pos))
    }

    /// Parse an optional number from position i. Returns (value if present else 0, consumed_index).
    fn parse_optional_number(
        &self,
        chars: &[u8],
        i: usize,
        len: usize,
    ) -> Result<(u32, usize), SmartsParseError> {
        // RDKit❗✔️: [0]		{ yylval->ival = 0;  return ZERO_TOKEN; }
        // RDKit❗✔️: [1-9]		{ yylval->ival = yytext[0]-'0';  return NONZERO_DIGIT_TOKEN; }
        // An absent token leaves the parser index unchanged; callers apply
        // the token-specific default. Present numbers share `parse_number` so
        // zero token width and signed-source bounds stay identical.
        // Local complexity review: absence is O(1); a present value delegates
        // to one O(n), O(1)-state fold with no intermediate allocation.
        if i >= len || !chars[i].is_ascii_digit() {
            return Ok((0, i));
        }
        self.parse_number(chars, i, len)
    }

    /// Parse `bond_expr` with bison's declared operator precedence.
    fn parse_bond_expr(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
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
        // RDKit❗✔️: branch_open_token: GROUP_OPEN_TOKEN { $$ = current_token_position; };
        // Local complexity review: this helper performs one token check and
        // advance in O(1), returning the original source position without a
        // scan, allocation, or duplicate branch parser.
        let (Token::OpenParen, position) = self.peek() else {
            return Err(self.unexpected_token_error(self.peek().1, "expected branch opening token"));
        };
        let position = position.parser_byte_end;
        self.advance();
        Ok(position)
    }

    fn parse_bond_semi_expr(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
        let mut bond = self.parse_bond_or_expr()?;
        while matches!(self.peek(), (Token::Semi, _)) {
            self.advance();
            let rhs = self.parse_bond_or_expr()?;
            bond.expand_query(rhs, CompositeQueryType::And);
        }
        Ok(bond)
    }

    fn parse_bond_or_expr(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
        let mut bond = self.parse_bond_and_expr()?;
        while matches!(self.peek(), (Token::Or, _)) {
            self.advance();
            let rhs = self.parse_bond_and_expr()?;
            bond.expand_query(rhs, CompositeQueryType::Or);
        }
        Ok(bond)
    }

    fn parse_bond_and_expr(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
        let mut bond = self.parse_bond_query()?;
        while matches!(self.peek(), (Token::And, _)) {
            self.advance();
            let rhs = self.parse_bond_query()?;
            bond.expand_query(rhs, CompositeQueryType::And);
        }
        Ok(bond)
    }

    fn parse_bond_query(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
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
        let mut bond = self.parse_bondd()?;
        while matches!(self.peek(), (Token::BondSpec(_), _) | (Token::Not, _)) {
            let rhs = self.parse_bondd()?;
            bond.expand_query(rhs, CompositeQueryType::And);
        }
        Ok(bond)
    }

    fn parse_bondd(&mut self) -> Result<ParsedSmartsBond, SmartsParseError> {
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
            (Token::BadCharacter(character), span) => {
                Err(Self::bad_character_error(*span, *character))
            }
            (Token::Not, _) => {
                self.advance();
                let mut bond = self.parse_bondd()?;
                let is_negated = matches!(&bond.query, QueryNode::Not(_));
                bond.query.set_negation(!is_negated);
                Ok(bond)
            }
            (Token::BondSpec(lexeme), _) => {
                let bond = parsed_bond_spec(*lexeme);
                self.advance();
                Ok(bond)
            }
            (_, span) => Err(self.unexpected_token_error(*span, "expected bond query primitive")),
        }
    }
}

// ---------------------------------------------------------------------------
// SMARTS primitive helpers
// ---------------------------------------------------------------------------

/// Reduce RDKit's `simple_atom` production into its query leaf and carrier.
fn parse_simple_atom(name: impl AsRef<[u8]>) -> Option<SimpleAtom> {
    let name = name.as_ref();
    // These lexer actions map to the same query leaves and carrier flags in
    // SimpleAtom below. The inline representation avoids heap allocations
    // for RDKit QueryAtom and query objects without changing those values.
    // RDKit✔️🔝: \*			{
    // RDKit✔️🔝: 	yylval->atom = new QueryAtom();
    // RDKit✔️🔝: 	yylval->atom->setQuery(makeAtomNullQuery());
    // RDKit✔️🔝: 	return SIMPLE_ATOM_QUERY_TOKEN;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: a			{
    // RDKit✔️🔝: 	yylval->atom = new QueryAtom();
    // RDKit✔️🔝: 	yylval->atom->setQuery(makeAtomAromaticQuery());
    // RDKit✔️🔝: 	yylval->atom->setIsAromatic(true);
    // RDKit✔️🔝: 	return SIMPLE_ATOM_QUERY_TOKEN;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: A			{
    // RDKit✔️🔝: 	yylval->atom = new QueryAtom();
    // RDKit✔️🔝: 	yylval->atom->setQuery(makeAtomAliphaticQuery());
    // RDKit✔️🔝: 	return SIMPLE_ATOM_QUERY_TOKEN;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: simple_atom: 	ORGANIC_ATOM_TOKEN {
    // RDKit✔️🔝:   //
    // RDKit✔️🔝:   // This construction (and some others) may seem odd, but the
    // RDKit✔️🔝:   // SMARTS definition requires that an atom which is aliphatic on
    // RDKit✔️🔝:   // input (i.e. something in the "organic subset" that is given with
    // RDKit✔️🔝:   // a capital letter) only match aliphatic atoms.
    // RDKit✔️🔝:   //
    // RDKit✔️🔝:   // The following rule applies a similar logic to aromatic atoms.
    // RDKit✔️🔝:   //
    // RDKit✔️🔝:   $$ = new QueryAtom($1);
    // RDKit✔️🔝:   $$->setQuery(makeAtomTypeQuery($1,false));
    // RDKit✔️🔝: }
    // RDKit✔️🔝: | AROMATIC_ATOM_TOKEN {
    // RDKit✔️🔝:   $$ = new QueryAtom($1);
    // RDKit✔️🔝:   $$->setIsAromatic(true);
    // RDKit✔️🔝:   $$->setQuery(makeAtomTypeQuery($1,true));
    // RDKit✔️🔝: }
    // RDKit✔️🔝: | SIMPLE_ATOM_QUERY_TOKEN
    // RDKit✔️🔝: ;
    // Local complexity review: symbol dispatch is over a fixed-size set and
    // constructs one inline typed leaf in O(1) time/space. This removes the
    // source QueryAtom/query-object allocations without adding traversal,
    // cloning, temporary collections, or a second simple-atom decoder.
    let atom_type = |atomic_number, aromatic| SimpleAtom {
        query: QueryNode::Predicate(AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        }),
        atomic_number,
        aromatic,
    };
    Some(match name {
        b"B" => atom_type(5, false),
        b"C" => atom_type(6, false),
        b"N" => atom_type(7, false),
        b"O" => atom_type(8, false),
        b"F" => atom_type(9, false),
        b"P" => atom_type(15, false),
        b"S" => atom_type(16, false),
        b"Cl" => atom_type(17, false),
        b"Br" => atom_type(35, false),
        b"I" => atom_type(53, false),
        b"b" => atom_type(5, true),
        b"c" => atom_type(6, true),
        b"n" => atom_type(7, true),
        b"o" => atom_type(8, true),
        b"p" => atom_type(15, true),
        b"s" => atom_type(16, true),
        b"si" => atom_type(14, true),
        b"as" => atom_type(33, true),
        b"se" => atom_type(34, true),
        b"te" => atom_type(52, true),
        b"*" => SimpleAtom {
            query: make_atom_null_query(),
            atomic_number: 0,
            aromatic: false,
        },
        b"a" => SimpleAtom {
            query: crate::query_behavior::make_atom_aromatic_query(),
            atomic_number: 0,
            aromatic: true,
        },
        b"A" => SimpleAtom {
            query: crate::query_behavior::make_atom_aliphatic_query(),
            atomic_number: 0,
            aromatic: false,
        },
        _ => return None,
    })
}

/// Decode either RDKit's `simple_atom` production or an `ATOM_TOKEN` emitted
/// by the bracket-atom lexer.
fn parse_atom_token(name: impl AsRef<[u8]>) -> Option<SimpleAtom> {
    let name = name.as_ref();
    parse_simple_atom(name).or_else(|| {
        element_symbol_to_atomic_number(name).map(|atomic_number| SimpleAtom {
            query: QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(atomic_number)),
            atomic_number,
            aromatic: false,
        })
    })
}

fn set_atom_carrier_identity(
    mut atom: ParsedAtomExpr,
    atomic_number: u32,
    aromatic: bool,
    position: usize,
) -> Result<ParsedAtomExpr, SmartsParseError> {
    // BEGIN RDKIT CPP FUNCTION Atom::setAtomicNum
    // RDKit❗✔️: void setAtomicNum(int newNum) { d_atomicNum = newNum; }
    // END RDKIT CPP FUNCTION Atom::setAtomicNum
    let atomic_number =
        u8::try_from(atomic_number).map_err(|_| SmartsParseError::InvalidAtomPrimitive {
            position,
            detail: "atomic number is outside the parser carrier's u8 range".to_string(),
        })?;
    atom.carrier = atom
        .carrier
        .with_identity(QueryAtomIdentity::from_atomic_number(atomic_number));
    atom.carrier.set_aromatic(aromatic);
    Ok(atom)
}

/// Convert a bond specifier character to a bond query predicate. Carrier type
/// and direction are tracked separately by parsed-bond construction and parser
/// actions.
fn bond_spec_to_query(lexeme: BondLexeme) -> QueryNode<BondQueryPredicate> {
    match lexeme {
        BondLexeme::DativeRight => make_bond_order_equals_query(BondOrder::DativeRight),
        BondLexeme::DativeLeft => make_bond_order_equals_query(BondOrder::DativeLeft),
        BondLexeme::Symbol(ch) => match ch {
            b'-' => make_bond_order_equals_query(BondOrder::Single),
            b'=' => make_bond_order_equals_query(BondOrder::Double),
            b'#' => make_bond_order_equals_query(BondOrder::Triple),
            b':' => make_bond_order_equals_query(BondOrder::Aromatic),
            b'$' => make_bond_order_equals_query(BondOrder::Quadruple),
            b'@' => make_bond_is_in_ring_query(),
            b'~' => make_bond_null_query(),
            b'/' | b'\\' => unspecified_smarts_bond_query(),
            _ => QueryNode::Predicate(BondQueryPredicate::Any),
        },
    }
}

fn parsed_bond_spec(lexeme: BondLexeme) -> ParsedSmartsBond {
    // BEGIN RDKIT CPP FUNCTION SMARTS BOND_TOKEN carrier actions
    // RDKit❗✔️: \=	{ yylval->bond = new QueryBond(Bond::DOUBLE);
    // RDKit❗✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::DOUBLE));
    // RDKit❗✔️: 	return BOND_TOKEN;  }
    // RDKit❗✔️: \~	{ yylval->bond = new QueryBond();
    // RDKit❗✔️: 	yylval->bond->setQuery(makeBondNullQuery());
    // RDKit❗✔️: 	return BOND_TOKEN;  }
    // RDKit❗✔️: \$	{ yylval->bond = new QueryBond(Bond::QUADRUPLE);
    // RDKit❗✔️: 	yylval->bond->setQuery(makeBondOrderEqualsQuery(Bond::QUADRUPLE));
    // RDKit❗✔️:     return BOND_TOKEN; }
    // RDKit❗✔️: [\\]{1,2}    { yylval->bond = new QueryBond(Bond::SINGLE);
    // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDDOWNRIGHT);
    // RDKit❗✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
    // RDKit❗✔️: 	return BOND_TOKEN;  }
    // RDKit❗✔️: [\/]    { yylval->bond = new QueryBond(Bond::SINGLE);
    // RDKit❗✔️: 	yylval->bond->setBondDir(Bond::ENDUPRIGHT);
    // RDKit❗✔️: 	yylval->bond->setQuery(makeSingleOrAromaticBondQuery());
    // RDKit❗✔️: 	return BOND_TOKEN;  }
    // RDKit❗✔️: \-\> {
    // RDKit❗✔️:     yylval->bond = new QueryBond(Bond::DATIVER);
    // RDKit❗✔️:     return BOND_TOKEN;
    // RDKit❗✔️: }
    // RDKit❗✔️: \<\- {
    // RDKit❗✔️:     yylval->bond = new QueryBond(Bond::DATIVEL);
    // RDKit❗✔️:     return BOND_TOKEN;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION SMARTS BOND_TOKEN carrier actions
    let carrier_order = match lexeme {
        BondLexeme::DativeRight => BondOrder::DativeRight,
        BondLexeme::DativeLeft => BondOrder::DativeLeft,
        BondLexeme::Symbol(ch) => match ch {
            b'-' | b'/' | b'\\' => BondOrder::Single,
            b'=' => BondOrder::Double,
            b'#' => BondOrder::Triple,
            b':' => BondOrder::Aromatic,
            b'$' => BondOrder::Quadruple,
            b'~' | b'@' => BondOrder::Unspecified,
            _ => BondOrder::Unspecified,
        },
    };
    ParsedSmartsBond {
        query: bond_spec_to_query(lexeme),
        carrier_order,
        unspecified_order: false,
    }
}

fn ring_closure_is_unspecified(bond: &ParsedSmartsBond) -> bool {
    // BEGIN RDKIT CPP FUNCTION CloseMolRings unspecified-order selection
    // RDKit❗❌: if (!bond1->hasProp(common_properties::_unspecifiedOrder)) {
    // END RDKIT CPP FUNCTION CloseMolRings unspecified-order selection
    // This bit mirrors the source property directly; query equality and the
    // ordinary carrier order are independent state and cannot stand in for it.
    // Local complexity: one field read, O(1).
    bond.unspecified_order
}

fn unspecified_smarts_bond_query() -> QueryNode<BondQueryPredicate> {
    crate::query_behavior::make_single_or_aromatic_bond_query()
}

fn normalize_dative_bond(bond: ParsedSmartsBond) -> (ParsedSmartsBond, bool) {
    // BEGIN RDKIT CPP FUNCTION QueryBond constructor and QueryBond::setBondType
    // RDKit✔️✔️: QueryBond::QueryBond(BondType bT) : Bond(bT) {
    // RDKit✔️✔️:   if (bT != Bond::UNSPECIFIED) {
    // RDKit✔️✔️:     dp_query = makeBondOrderEqualsQuery(bT);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     dp_query = makeBondNullQuery();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // RDKit✔️✔️: void QueryBond::setBondType(BondType bT) {
    // RDKit✔️✔️:   // NOTE: calling this blows out any existing query
    // RDKit✔️✔️:   d_bondType = bT;
    // RDKit✔️✔️:   delete dp_query;
    // RDKit✔️✔️:   dp_query = nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   dp_query = makeBondOrderEqualsQuery(bT);
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION QueryBond constructor and QueryBond::setBondType
    // RDKit✔️✔️:   if( $2->getBondType() == Bond::DATIVER ){
    // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx1);
    // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx2);
    // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️:   }else if ( $2->getBondType() == Bond::DATIVEL ){
    // RDKit✔️✔️:     $2->setBeginAtomIdx(atomIdx2);
    // RDKit✔️✔️:     $2->setEndAtomIdx(atomIdx1);
    // RDKit✔️✔️:     $2->setBondType(Bond::DATIVE);
    // RDKit✔️✔️: }
    // The ordinary and branch reductions hold QueryBond*, so normalization
    // replaces even a composite or negated query. CloseMolRings holds Bond*;
    // its nonvirtual setter changes only the carrier and is handled separately.
    // Complexity: constant-size dispatch and one equality leaf, like upstream;
    // dropping the previous tree visits its nodes, like delete dp_query.
    match bond.carrier_order {
        BondOrder::DativeRight => (
            ParsedSmartsBond {
                query: make_bond_order_equals_query(BondOrder::Dative),
                carrier_order: BondOrder::Dative,
                ..bond
            },
            false,
        ),
        BondOrder::DativeLeft => (
            ParsedSmartsBond {
                query: make_bond_order_equals_query(BondOrder::Dative),
                carrier_order: BondOrder::Dative,
                ..bond
            },
            true,
        ),
        _ => (bond, false),
    }
}

/// Look up the atomic number for an element symbol.
///
/// RDKit✔️✔️: Standard periodic table mapping.
fn element_symbol_to_atomic_number(symbol: impl AsRef<[u8]>) -> Option<u8> {
    let symbol = symbol.as_ref();
    // RDKit source: third_party/rdkit/Code/GraphMol/atomic_data.cpp
    // RDKit❗✔️: // we leave Uut and Uup in here for backwards
    // RDKit❗✔️: // compatibility. Nh and Mc (the entries appearing first
    // RDKit❗✔️: // for a particular atomic number) will be the values returned
    // RDKit❗✔️: // when looking an atomic symbol up using atomic number.
    // RDKit❗✔️: 113 Nh	7	1.36	0	2.0	284	2	284	284.17873	-1
    // RDKit❗✔️: 113 Uut	7	1.36	0	2.0	284	2	284	284.17873	-1
    // RDKit❗✔️: 115 Mc	7	1.62	0	2.0	288	2	288	288.19274	-1
    // RDKit❗✔️: 115 Uup	7	1.62	0	2.0	288	2	288	288.19274	-1
    match symbol {
        b"H" => Some(1),
        b"He" => Some(2),
        b"Li" => Some(3),
        b"Be" => Some(4),
        b"B" => Some(5),
        b"C" => Some(6),
        b"N" => Some(7),
        b"O" => Some(8),
        b"F" => Some(9),
        b"Ne" => Some(10),
        b"Na" => Some(11),
        b"Mg" => Some(12),
        b"Al" => Some(13),
        b"Si" => Some(14),
        b"P" => Some(15),
        b"S" => Some(16),
        b"Cl" => Some(17),
        b"Ar" => Some(18),
        b"K" => Some(19),
        b"Ca" => Some(20),
        b"Sc" => Some(21),
        b"Ti" => Some(22),
        b"V" => Some(23),
        b"Cr" => Some(24),
        b"Mn" => Some(25),
        b"Fe" => Some(26),
        b"Co" => Some(27),
        b"Ni" => Some(28),
        b"Cu" => Some(29),
        b"Zn" => Some(30),
        b"Ga" => Some(31),
        b"Ge" => Some(32),
        b"As" => Some(33),
        b"Se" => Some(34),
        b"Br" => Some(35),
        b"Kr" => Some(36),
        b"Rb" => Some(37),
        b"Sr" => Some(38),
        b"Y" => Some(39),
        b"Zr" => Some(40),
        b"Nb" => Some(41),
        b"Mo" => Some(42),
        b"Tc" => Some(43),
        b"Ru" => Some(44),
        b"Rh" => Some(45),
        b"Pd" => Some(46),
        b"Ag" => Some(47),
        b"Cd" => Some(48),
        b"In" => Some(49),
        b"Sn" => Some(50),
        b"Sb" => Some(51),
        b"Te" => Some(52),
        b"I" => Some(53),
        b"Xe" => Some(54),
        b"Cs" => Some(55),
        b"Ba" => Some(56),
        b"La" => Some(57),
        b"Ce" => Some(58),
        b"Pr" => Some(59),
        b"Nd" => Some(60),
        b"Pm" => Some(61),
        b"Sm" => Some(62),
        b"Eu" => Some(63),
        b"Gd" => Some(64),
        b"Tb" => Some(65),
        b"Dy" => Some(66),
        b"Ho" => Some(67),
        b"Er" => Some(68),
        b"Tm" => Some(69),
        b"Yb" => Some(70),
        b"Lu" => Some(71),
        b"Hf" => Some(72),
        b"Ta" => Some(73),
        b"W" => Some(74),
        b"Re" => Some(75),
        b"Os" => Some(76),
        b"Ir" => Some(77),
        b"Pt" => Some(78),
        b"Au" => Some(79),
        b"Hg" => Some(80),
        b"Tl" => Some(81),
        b"Pb" => Some(82),
        b"Bi" => Some(83),
        b"Po" => Some(84),
        b"At" => Some(85),
        b"Rn" => Some(86),
        b"Fr" => Some(87),
        b"Ra" => Some(88),
        b"Ac" => Some(89),
        b"Th" => Some(90),
        b"Pa" => Some(91),
        b"U" => Some(92),
        b"Np" => Some(93),
        b"Pu" => Some(94),
        b"Am" => Some(95),
        b"Cm" => Some(96),
        b"Bk" => Some(97),
        b"Cf" => Some(98),
        b"Es" => Some(99),
        b"Fm" => Some(100),
        b"Md" => Some(101),
        b"No" => Some(102),
        b"Lr" => Some(103),
        b"Rf" => Some(104),
        b"Db" => Some(105),
        b"Sg" => Some(106),
        b"Bh" => Some(107),
        b"Hs" => Some(108),
        b"Mt" => Some(109),
        b"Ds" => Some(110),
        b"Rg" => Some(111),
        b"Cn" => Some(112),
        b"Nh" => Some(113),
        b"Uut" => Some(113),
        b"Fl" => Some(114),
        b"Mc" => Some(115),
        b"Uup" => Some(115),
        b"Lv" => Some(116),
        b"Ts" => Some(117),
        b"Og" => Some(118),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod query_hydrogen_merge_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::*;
    use cosmolkit_model::{
        AtomId, BondId, Conformer2D, Conformer3D, SGroupAttachPoint, SGroupBondRole, SGroupBracket,
        SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData, SGroupDisplay, StereoGroup,
        StereoGroupKind, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind, TemplateAttachment,
        TemplateAttachmentOrder, TemplateAttachmentOrderError, query_substance_groups,
        replace_query_substance_groups,
    };

    // H-C-C chain: atom 0 is the carrier, atom 1 is a removable query
    // hydrogen, atom 2 is a later surviving atom.
    fn chain_h_c_c(carrier: Option<TemplateAttachmentOrder>) -> QueryGraph {
        let atom_0 = match carrier {
            Some(order) => QueryAtom::new(
                AtomId::new(0),
                AtomSpec::new(Element::C).with_template_attachment_order(order),
            ),
            None => QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
        };
        let atom_1 = QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::H));
        let atom_2 = QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::C));
        let bonds = vec![
            QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryBond::new(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
        ];
        QueryGraph::from_parts(
            vec![atom_0, atom_1, atom_2],
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("valid H-C-C query graph")
    }

    // H-C(<)-C(-<C) with two later surviving targets for ordered labels.
    fn branched_h_c_c_c(carrier: TemplateAttachmentOrder) -> QueryGraph {
        let atom_0 = QueryAtom::new(
            AtomId::new(0),
            AtomSpec::new(Element::C).with_template_attachment_order(carrier),
        );
        let atom_1 = QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::H));
        let atom_2 = QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::C));
        let atom_3 = QueryAtom::new(AtomId::new(3), AtomSpec::new(Element::C));
        let bonds = vec![
            QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryBond::new(
                BondId::new(1),
                BondSpec::new(AtomId::new(0), AtomId::new(2), BondOrder::Single),
            ),
            QueryBond::new(
                BondId::new(2),
                BondSpec::new(AtomId::new(0), AtomId::new(3), BondOrder::Single),
            ),
        ];
        QueryGraph::from_parts(
            vec![atom_0, atom_1, atom_2, atom_3],
            bonds,
            Vec::<(
                cosmolkit_model::PropertyText,
                cosmolkit_model::PropertyValue,
            )>::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("valid branched query graph")
    }

    fn query_hydrogen_graph_with_conformers() -> QueryGraph {
        let source =
            parse_smarts("[C]([H])[C]", &SmartsParseParams::default()).expect("source query graph");
        let mut props = source.props().clone();
        props.insert("_CXSMILES_Data".into(), "|source-prefix|".into());
        QueryGraph::from_parts(
            source.atoms().to_vec(),
            source.bonds().to_vec(),
            props,
            vec![
                Conformer2D::new(9, vec![[0.0, 0.0], [1.0, 0.5], [2.0, 1.0]])
                    .with_prop("frame", "first-2d"),
                Conformer2D::new(11, vec![[0.0, 1.0], [1.0, 1.5], [2.0, 2.0]])
                    .with_prop("frame", "second-2d"),
            ],
            vec![
                Conformer3D::new(
                    17,
                    vec![[0.0, 0.0, 0.0], [1.0, 0.5, 0.0], [2.0, 1.0, 0.0]],
                    true,
                )
                .with_prop("frame", "first-3d"),
                Conformer3D::new(
                    19,
                    vec![[0.0, 1.0, 0.0], [1.0, 1.5, 0.0], [2.0, 2.0, 0.0]],
                    false,
                )
                .with_prop("frame", "second-3d"),
            ],
            vec![
                StereoGroup::new(
                    StereoGroupKind::Or,
                    vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)],
                    vec![BondId::new(0), BondId::new(1)],
                )
                .with_id(17),
            ],
        )
        .expect("query graph with complete conformer and stereo state")
    }

    #[test]
    fn carrier_referencing_a_later_surviving_atom_is_remapped() {
        let order =
            TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(2), "Al")])
                .expect("attachment order");
        let graph = chain_h_c_c(Some(order));
        let merged = merge_query_hs(&graph, false, false).expect("merge query hydrogens");
        assert_eq!(merged.num_atoms(), 2);
        let carrier = merged.atoms()[0]
            .template_attachment_order()
            .expect("carrier keeps its attachment order");
        assert_eq!(carrier.entries().len(), 1);
        assert_eq!(carrier.entries()[0].target(), AtomId::new(1));
        assert_eq!(carrier.entries()[0].label(), "Al");
    }

    #[test]
    fn carrier_referencing_a_removed_target_fails_loudly() {
        let order =
            TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(1), "H")])
                .expect("attachment order");
        let graph = chain_h_c_c(Some(order));
        match merge_query_hs(&graph, false, false) {
            Err(SmartsParseError::TemplateAttachmentRemap {
                carrier,
                source: TemplateAttachmentOrderError::TargetRemoved { target, .. },
            }) => {
                assert_eq!(carrier, 0);
                assert_eq!(target, AtomId::new(1));
            }
            other => panic!("expected TargetRemoved for carrier 0, got {other:?}"),
        }
    }

    #[test]
    fn merge_without_attachments_keeps_the_source_predicate_expansion() {
        let graph = chain_h_c_c(None);
        let merged = merge_query_hs(&graph, false, false).expect("merge query hydrogens");
        assert_eq!(merged.num_atoms(), 2);
        assert_eq!(
            merged.atoms()[0].predicate(),
            &QueryNode::And(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::HydrogenCount(0)
                ))),
            ])
        );
    }

    #[test]
    fn ordered_attachment_labels_are_preserved_through_renumbering() {
        let order = TemplateAttachmentOrder::new(vec![
            TemplateAttachment::new(AtomId::new(2), "Al"),
            TemplateAttachment::new(AtomId::new(3), "Br"),
        ])
        .expect("attachment order");
        let graph = branched_h_c_c_c(order);
        let merged = merge_query_hs(&graph, false, false).expect("merge query hydrogens");
        assert_eq!(merged.num_atoms(), 3);
        let carrier = merged.atoms()[0]
            .template_attachment_order()
            .expect("carrier keeps its attachment order");
        assert_eq!(carrier.entries().len(), 2);
        assert_eq!(carrier.entries()[0].target(), AtomId::new(1));
        assert_eq!(carrier.entries()[0].label(), "Al");
        assert_eq!(carrier.entries()[1].target(), AtomId::new(2));
        assert_eq!(carrier.entries()[1].label(), "Br");
    }

    #[test]
    fn failed_merge_leaves_the_in_place_graph_unchanged() {
        let order =
            TemplateAttachmentOrder::new(vec![TemplateAttachment::new(AtomId::new(1), "H")])
                .expect("attachment order");
        let mut graph = chain_h_c_c(Some(order));
        let before = graph.clone();
        assert!(merge_query_hs_in_place(&mut graph, false, false).is_err());
        assert_eq!(graph, before);
    }

    #[test]
    fn parse_smarts_merges_plain_query_hydrogens() {
        let params = SmartsParseParams {
            merge_hs: true,
            ..SmartsParseParams::default()
        };
        let merged = parse_smarts("[C][H]", &params).expect("parse and merge");
        assert_eq!(merged.num_atoms(), 1);
        // The parser lowers a bracketed aliphatic atom to `AtomType`, so the
        // merged predicate retains that node as the first `And` child.
        assert_eq!(
            merged.atoms()[0].predicate(),
            &QueryNode::And(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 6,
                    aromatic: false,
                }),
                QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::HydrogenCount(0)
                ))),
            ])
        );
    }

    #[test]
    fn merge_query_hydrogens_match_map_presence_and_isotope_flags() {
        let source = parse_smarts(
            "[C]([H])([H:0])([H:5])([2H])([2H:0])",
            &SmartsParseParams::default(),
        )
        .expect("parse mapped and isotopic query H atoms without merging");
        assert_eq!(source.num_atoms(), 6);
        assert_eq!(source.atoms()[1].atom_map(), None);
        assert_eq!(source.atoms()[2].atom_map(), Some(0));
        assert_eq!(source.atoms()[3].atom_map(), Some(5));
        assert_eq!(source.atoms()[4].isotope(), Some(2));
        assert_eq!(source.atoms()[4].atom_map(), None);
        assert_eq!(source.atoms()[5].isotope(), Some(2));
        assert_eq!(source.atoms()[5].atom_map(), Some(0));
        let source_before = source.clone();

        let cases = [
            (false, false, vec![0, 4, 5]),
            (false, true, vec![0]),
            (true, false, vec![0, 2, 3, 4, 5]),
            (true, true, vec![0, 2, 3, 5]),
        ];
        for (merge_unmapped_only, merge_isotopes, expected_source_atoms) in cases {
            let merged = merge_query_hs(&source, merge_unmapped_only, merge_isotopes)
                .expect("copying wrapper merges source-selected H atoms");
            assert_eq!(
                source, source_before,
                "copying wrapper input must be unchanged"
            );
            assert_eq!(merged.num_atoms(), expected_source_atoms.len());
            assert_eq!(merged.num_bonds(), expected_source_atoms.len() - 1);

            let removed_hydrogen_count = source.num_atoms() - expected_source_atoms.len();
            let mut expected_predicate = QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            });
            // MolOps::mergeQueryHs calls expandQuery for EACH removed H.
            // QueryAtom::expandQuery builds a binary COMPOSITE_AND each time;
            // it does not flatten the already existing composite tree.
            for hydrogen_count in 0..removed_hydrogen_count {
                expected_predicate = QueryNode::And(vec![
                    expected_predicate,
                    QueryNode::Not(Box::new(QueryNode::Predicate(
                        AtomQueryPredicate::HydrogenCount(hydrogen_count as i32),
                    ))),
                ]);
            }
            assert_eq!(
                merged.atoms()[0].predicate(),
                &expected_predicate,
                "flags: merge_unmapped_only={merge_unmapped_only}, merge_isotopes={merge_isotopes}"
            );

            for (new_index, &old_index) in expected_source_atoms.iter().enumerate() {
                let actual = &merged.atoms()[new_index];
                let original = &source.atoms()[old_index];
                assert_eq!(actual.id(), AtomId::new(new_index));
                assert_eq!(actual.identity(), original.identity());
                assert_eq!(actual.atomic_number(), original.atomic_number());
                assert_eq!(actual.isotope(), original.isotope());
                assert_eq!(actual.atom_map(), original.atom_map());
                if old_index != 0 {
                    assert_eq!(actual.predicate(), original.predicate());
                    let mapped_bond = &merged.bonds()[new_index - 1];
                    let source_bond = source
                        .bonds()
                        .iter()
                        .find(|bond| {
                            bond.begin() == AtomId::new(0) && bond.end() == AtomId::new(old_index)
                        })
                        .expect("source bond to retained hydrogen");
                    assert_eq!(mapped_bond.id(), BondId::new(new_index - 1));
                    assert_eq!(mapped_bond.begin(), AtomId::new(0));
                    assert_eq!(mapped_bond.end(), AtomId::new(new_index));
                    assert_eq!(mapped_bond.predicate(), source_bond.predicate());
                }
            }

            let mut in_place = source.clone();
            merge_query_hs_in_place(&mut in_place, merge_unmapped_only, merge_isotopes)
                .expect("in-place helper applies identical source flags");
            assert_eq!(in_place, merged);
        }
    }

    #[test]
    fn query_sgroups_hydrogen_remap_preserves_typed_references_and_all_graph_state() {
        let mut graph = query_hydrogen_graph_with_conformers();
        let group = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
            .with_rdkit_sequence_id(23)
            .with_external_id(41)
            .with_atoms(vec![AtomId::new(0), AtomId::new(2), AtomId::new(0)])
            .with_bonds(vec![BondId::new(1), BondId::new(1)])
            .with_bond_role(BondId::new(1), SGroupBondRole::Contained)
            .with_head_crossing_bonds(vec![BondId::new(1), BondId::new(1)])
            .with_crossing_bond_correspondence(vec![BondId::new(1)])
            .with_parent_atoms(vec![AtomId::new(2)])
            .with_label("typed polymer data")
            .with_connection(SGroupConnection::HeadToTail)
            .with_subtype("SUP")
            .with_bracket_style(SGroupBracketStyle::Bracket)
            .with_display(SGroupDisplay {
                brackets: vec![SGroupBracket::new([
                    [0.0, 0.0, 0.0],
                    [1.0, 0.0, 0.0],
                    [0.0, 1.0, 0.0],
                ])],
                field_position: Some([0.25, 0.75]),
                display_tag: Some("typed-display".into()),
            })
            .with_expansion_state("expanded")
            .with_class("polymer-class")
            .with_component_number(4)
            .with_data(SGroupData {
                field_name: Some("FIELD".into()),
                field_type: Some("S".into()),
                values: vec!["one".into(), "two".into()],
                ..SGroupData::default()
            })
            .with_attach_points(vec![SGroupAttachPoint {
                atom: AtomId::new(0),
                leaving_atom: Some(AtomId::new(2)),
                label: Some("R".into()),
                order: Some(2),
            }])
            .with_cstates(vec![SGroupCState::new(BondId::new(1), [0.25, 0.5, 0.75])])
            .with_prop("custom", "preserved")
            .expect("original valid fixture constructor must succeed")
            .with_data_field("first raw row")
            .with_data_field("second raw row");
        replace_query_substance_groups(&mut graph, vec![group])
            .expect("typed query SGroup is valid for the source graph");
        let source = graph.clone();

        let merged = merge_query_hs(&graph, false, false).expect("merge query hydrogens");

        let expected_group =
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Superatom)
                .with_rdkit_sequence_id(23)
                .with_external_id(41)
                .with_atoms(vec![AtomId::new(0), AtomId::new(1), AtomId::new(0)])
                .with_bonds(vec![BondId::new(0), BondId::new(0)])
                .with_bond_role(BondId::new(0), SGroupBondRole::Contained)
                .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
                .with_crossing_bond_correspondence(vec![BondId::new(0)])
                .with_parent_atoms(vec![AtomId::new(1)])
                .with_label("typed polymer data")
                .with_connection(SGroupConnection::HeadToTail)
                .with_subtype("SUP")
                .with_bracket_style(SGroupBracketStyle::Bracket)
                .with_display(SGroupDisplay {
                    brackets: vec![SGroupBracket::new([
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                    ])],
                    field_position: Some([0.25, 0.75]),
                    display_tag: Some("typed-display".into()),
                })
                .with_expansion_state("expanded")
                .with_class("polymer-class")
                .with_component_number(4)
                .with_data(SGroupData {
                    field_name: Some("FIELD".into()),
                    field_type: Some("S".into()),
                    values: vec!["one".into(), "two".into()],
                    ..SGroupData::default()
                })
                .with_attach_points(vec![SGroupAttachPoint {
                    atom: AtomId::new(0),
                    leaving_atom: Some(AtomId::new(1)),
                    label: Some("R".into()),
                    order: Some(2),
                }])
                .with_cstates(vec![SGroupCState::new(BondId::new(0), [0.25, 0.5, 0.75])])
                .with_prop("custom", "preserved")
                .expect("original valid fixture constructor must succeed")
                .with_data_field("first raw row")
                .with_data_field("second raw row");

        assert_eq!(merged.num_atoms(), 2);
        assert_eq!(merged.num_bonds(), 1);
        assert_eq!(query_substance_groups(&merged), &[expected_group]);
        assert_eq!(
            merged.stereo_groups(),
            &[StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(1)],
                vec![BondId::new(0)],
            )
            .with_id(17),]
        );
        assert_eq!(merged.props(), source.props());
        for (new, old) in [(0, 0), (1, 2)] {
            assert_eq!(
                merged.atoms()[new].predicate_is_carrier_derived(),
                source.atoms()[old].predicate_is_carrier_derived()
            );
        }
        assert_eq!(
            merged.bonds()[0].predicate_is_carrier_derived(),
            source.bonds()[1].predicate_is_carrier_derived()
        );
        assert_eq!(merged.bonds()[0].predicate(), source.bonds()[1].predicate());

        let coordinates = merged.coordinate_block(None);
        assert_eq!(coordinates.conformers_2d.len(), 2);
        assert_eq!(coordinates.conformers_2d[0].id(), 9);
        assert_eq!(
            coordinates.conformers_2d[0]
                .props()
                .get(b"frame".as_slice())
                .map(crate::smarts_parse::fixture_text),
            Some("first-2d")
        );
        assert_eq!(
            coordinates.conformers_2d[0].coordinates(),
            &[[0.0, 0.0], [2.0, 1.0]]
        );
        assert_eq!(coordinates.conformers_2d[1].id(), 11);
        assert_eq!(
            coordinates.conformers_2d[1]
                .props()
                .get(b"frame".as_slice())
                .map(crate::smarts_parse::fixture_text),
            Some("second-2d")
        );
        assert_eq!(
            coordinates.conformers_2d[1].coordinates(),
            &[[0.0, 1.0], [2.0, 2.0]]
        );
        assert_eq!(coordinates.conformers_3d.len(), 2);
        assert_eq!(coordinates.conformers_3d[0].id(), 17);
        assert!(coordinates.conformers_3d[0].is_3d());
        assert_eq!(
            coordinates.conformers_3d[0]
                .props()
                .get(b"frame".as_slice())
                .map(crate::smarts_parse::fixture_text),
            Some("first-3d")
        );
        assert_eq!(
            coordinates.conformers_3d[0].coordinates(),
            &[[0.0, 0.0, 0.0], [2.0, 1.0, 0.0]]
        );
        assert_eq!(coordinates.conformers_3d[1].id(), 19);
        assert!(!coordinates.conformers_3d[1].is_3d());
        assert_eq!(
            coordinates.conformers_3d[1]
                .props()
                .get(b"frame".as_slice())
                .map(crate::smarts_parse::fixture_text),
            Some("second-3d")
        );
        assert_eq!(
            coordinates.conformers_3d[1].coordinates(),
            &[[0.0, 1.0, 0.0], [2.0, 2.0, 0.0]]
        );
    }

    #[test]
    fn query_sgroups_hydrogen_deletion_cascades_parents_and_compacts_ids() {
        let mut graph = chain_h_c_c(None);
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(1)]),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_parent(SubstanceGroupId::new(2))
                .with_atoms(vec![AtomId::new(0)]),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(2)]),
            SubstanceGroup::new(SubstanceGroupId::new(3), SubstanceGroupKind::Data)
                .with_bonds(vec![BondId::new(0)]),
            SubstanceGroup::new(SubstanceGroupId::new(4), SubstanceGroupKind::Data)
                .with_parent(SubstanceGroupId::new(0))
                .with_atoms(vec![AtomId::new(0)]),
        ];
        replace_query_substance_groups(&mut graph, groups)
            .expect("source graph SGroup references are valid");

        let merged = merge_query_hs(&graph, false, false).expect("merge query hydrogens");

        let groups = query_substance_groups(&merged);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].id(), SubstanceGroupId::new(0));
        assert_eq!(groups[0].parent(), Some(SubstanceGroupId::new(1)));
        assert_eq!(groups[0].atoms(), &[AtomId::new(0)]);
        assert_eq!(groups[1].id(), SubstanceGroupId::new(1));
        assert_eq!(groups[1].parent(), None);
        assert_eq!(groups[1].atoms(), &[AtomId::new(1)]);
    }
}

#[cfg(test)]
mod q03_setup_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::*;

    #[test]
    fn q03_setup_window_matches_signed_plain_char_trim_and_terminal_read() {
        assert_eq!(
            setup_smarts_input(""),
            SmartsScannerInputWindow {
                byte_start: 0,
                byte_end: 1,
            }
        );
        assert_eq!(
            setup_smarts_input(" \t\r\n"),
            SmartsScannerInputWindow {
                byte_start: 4,
                byte_end: 5,
            }
        );
        assert_eq!(
            setup_smarts_input("\t C\t "),
            SmartsScannerInputWindow {
                byte_start: 2,
                byte_end: 3,
            }
        );
        // The unchanged U+2603 fixture supplies three UTF-8 bytes. Source
        // signed-char trimming skips them at the edges, then copies the
        // counted std::string terminator as one input NUL before Flex sentinels.
        assert_eq!(
            setup_smarts_input("☃"),
            SmartsScannerInputWindow {
                byte_start: 3,
                byte_end: 4,
            }
        );
        assert_eq!(
            setup_smarts_input("C☃C"),
            SmartsScannerInputWindow {
                byte_start: 0,
                byte_end: 5,
            }
        );
    }

    #[test]
    fn q03_scanner_keeps_trimmed_window_in_original_character_coordinates() {
        let input = "\t C\t ";
        let scanned = SmartsScanner::new(
            input.as_bytes(),
            ScannerStart::Molecule,
            setup_smarts_input(input),
        )
        .scan()
        .expect("edge controls are trimmed before scanning");
        assert_eq!(
            scanned[1].token,
            ScannerToken::OrganicElement("C".to_owned())
        );
        assert_eq!(
            (
                scanned[1].span.input_char_start,
                scanned[1].span.input_char_end
            ),
            (2, 3)
        );
        assert_eq!(
            scanned.last().expect("terminal EOS").token,
            ScannerToken::EndOfStream
        );
        assert_eq!(
            (
                scanned.last().unwrap().span.input_char_start,
                scanned.last().unwrap().span.input_char_end
            ),
            (3, 3)
        );
    }

    fn q03_scan_and_compact(input: &str) -> Vec<(Token, SmartsTokenSpan)> {
        let scanned = SmartsScanner::new(
            input.as_bytes(),
            ScannerStart::Molecule,
            setup_smarts_input(input),
        )
        .scan()
        .expect("valid scanner input");
        compact_scanned_tokens(input.as_bytes(), &scanned).expect("valid token stream")
    }

    #[test]
    fn q03_ascii_multichar_and_bracket_spans_keep_their_coordinate_origins() {
        let input = "ClBr";
        let scanned = SmartsScanner::new(
            input.as_bytes(),
            ScannerStart::Molecule,
            setup_smarts_input(input),
        )
        .scan()
        .expect("two source organic atoms");
        assert_eq!(
            scanned[1].token,
            ScannerToken::OrganicElement("Cl".to_owned())
        );
        assert_eq!(
            scanned[1].span,
            SmartsTokenSpan {
                input_char_start: 0,
                input_char_end: 2,
                parser_byte_start: 0,
                parser_byte_end: 2,
            }
        );
        assert_eq!(
            scanned[2].token,
            ScannerToken::OrganicElement("Br".to_owned())
        );
        assert_eq!(
            scanned[2].span,
            SmartsTokenSpan {
                input_char_start: 2,
                input_char_end: 4,
                parser_byte_start: 2,
                parser_byte_end: 4,
            }
        );

        let bracketed = q03_scan_and_compact("[a]");
        let (Token::BracketContent(content), bracket_span) = &bracketed[0] else {
            panic!("bracket content token")
        };
        assert_eq!(content.text, b"a".to_vec());
        assert_eq!(
            *bracket_span,
            SmartsTokenSpan {
                input_char_start: 0,
                input_char_end: 3,
                parser_byte_start: 0,
                parser_byte_end: 3,
            }
        );
        assert_eq!(
            content.span,
            SmartsTokenSpan {
                input_char_start: 1,
                input_char_end: 2,
                parser_byte_start: 1,
                parser_byte_end: 2,
            }
        );
        assert_eq!(content.lexical_tokens.len(), 1);
        assert_eq!(
            content.lexical_tokens[0].token,
            ScannerToken::SimpleAtomQuery(b'a')
        );
        assert_eq!(
            content.lexical_tokens[0].span,
            SmartsTokenSpan {
                input_char_start: 0,
                input_char_end: 1,
                parser_byte_start: 0,
                parser_byte_end: 1,
            }
        );

        let recursive = q03_scan_and_compact("[$([a])]");
        let (Token::BracketContent(content), _) = &recursive[0] else {
            panic!("recursive bracket content token")
        };
        assert_eq!(content.text, b"$([a])".to_vec());
        assert_eq!(content.lexical_tokens.len(), 1);
        assert_eq!(
            content.lexical_tokens[0].token,
            ScannerToken::SimpleAtomQuery(b'a')
        );
        assert_eq!(
            content.lexical_tokens[0].span,
            SmartsTokenSpan {
                input_char_start: 3,
                input_char_end: 4,
                parser_byte_start: 3,
                parser_byte_end: 4,
            }
        );
    }

    #[test]
    fn q03_trimmed_tokens_keep_original_character_and_parser_byte_offsets() {
        let input = "\t ClBr\t ";
        let scanned = SmartsScanner::new(
            input.as_bytes(),
            ScannerStart::Molecule,
            setup_smarts_input(input),
        )
        .scan()
        .expect("edge controls are trimmed before scanning");
        assert_eq!(
            scanned[1].span,
            SmartsTokenSpan {
                input_char_start: 2,
                input_char_end: 4,
                parser_byte_start: 0,
                parser_byte_end: 2,
            }
        );
        assert_eq!(
            scanned[2].span,
            SmartsTokenSpan {
                input_char_start: 4,
                input_char_end: 6,
                parser_byte_start: 2,
                parser_byte_end: 4,
            }
        );
        assert_eq!(
            scanned.last().expect("EOS token").span,
            SmartsTokenSpan {
                input_char_start: 6,
                input_char_end: 6,
                parser_byte_start: 4,
                parser_byte_end: 4,
            }
        );
    }

    #[test]
    fn q03_dispatch_tokens_keep_source_consumption_positions() {
        let newline =
            SmartsScanner::new(b"C\nC", ScannerStart::Molecule, setup_smarts_input("C\nC"))
                .scan()
                .expect("newline returns EOS without scanning the suffix");
        assert_eq!(newline.len(), 3);
        assert_eq!(newline[2].token, ScannerToken::EndOfStream);
        assert_eq!(
            newline[2].span,
            SmartsTokenSpan {
                input_char_start: 1,
                input_char_end: 2,
                parser_byte_start: 1,
                parser_byte_end: 2,
            }
        );

        let ascii = SmartsScanner::new(b"C?C", ScannerStart::Molecule, setup_smarts_input("C?C"))
            .scan()
            .expect("BAD_CHARACTER is transported to parser dispatch");
        assert_eq!(ascii.len(), 3);
        assert_eq!(ascii[2].token, ScannerToken::BadCharacter(b'?'));
        assert_eq!(
            ascii[2].span,
            SmartsTokenSpan {
                input_char_start: 1,
                input_char_end: 2,
                parser_byte_start: 1,
                parser_byte_end: 2,
            }
        );

        let multibyte = SmartsScanner::new(
            "C☃C".as_bytes(),
            ScannerStart::Molecule,
            setup_smarts_input("C☃C"),
        )
        .scan()
        .expect("first BAD_CHARACTER stops Flex-style scanning");
        assert_eq!(multibyte.len(), 3);
        assert_eq!(multibyte[2].token, ScannerToken::BadCharacter(0xE2));
        assert_eq!(
            multibyte[2].span,
            SmartsTokenSpan {
                input_char_start: 1,
                input_char_end: 2,
                parser_byte_start: 1,
                parser_byte_end: 2,
            }
        );
    }

    #[test]
    fn q03_bad_character_dispatch_uses_existing_atom_and_bond_starts() {
        assert_eq!(
            parse_atom_entry("C?").expect_err("atom start preserves BAD_CHARACTER"),
            SmartsParseError::UnexpectedCharacter {
                position: 2,
                character: b'?',
                context: "unexpected character in SMARTS string".to_owned(),
            }
        );
        assert_eq!(
            parse_bond_entry("-?").expect_err("bond start preserves BAD_CHARACTER"),
            SmartsParseError::UnexpectedCharacter {
                position: 2,
                character: b'?',
                context: "unexpected character in SMARTS string".to_owned(),
            }
        );
    }
}

#[cfg(test)]
mod cx_progress_coordinate_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{apply_cx_progress_to_query, parse_smarts_graph};

    #[test]
    fn cx_progress_coordinates_retains_committed_rows_and_pads_graph_atoms() {
        let mut graph = parse_smarts_graph("CC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|(1,2,3;bad,0,0)|");
        assert!(!progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("committed CX effects");

        assert_eq!(graph.conformers_3d().len(), 1);
        let conformer = &graph.conformers_3d()[0];
        assert_eq!(conformer.coordinates(), &[[1.0, 2.0, 3.0], [0.0, 0.0, 0.0]]);
        assert!(conformer.is_3d());
    }

    #[test]
    fn cx_progress_coordinates_projects_complete_rows_to_graph_atom_count() {
        let mut graph = parse_smarts_graph("CC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|(1,2,3;4,5,6;7,8,9)|");
        assert!(progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("complete CX effects");

        assert_eq!(graph.conformers_3d().len(), 1);
        assert_eq!(
            graph.conformers_3d()[0].coordinates(),
            &[[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        );
    }
}

#[cfg(test)]
mod cx_progress_label_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_model::{
        AtomId, BondId, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
        query_substance_groups, replace_query_substance_groups,
    };

    #[test]
    fn cx_progress_labels_and_values_apply_written_slots_and_decoded_text() {
        let mut graph = parse_smarts_graph("CCC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let labels = cosmolkit_cx::parse_cx_extensions_progress("|$a&#321;b;;last$|");
        assert!(labels.is_complete());
        apply_cx_progress_to_query(&mut graph, &labels).expect("label effects");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "aAb".to_owned().into()
            ))
        );
        assert_eq!(graph.atom(1).and_then(|atom| atom.prop("atomLabel")), None);
        assert_eq!(
            graph.atom(2).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "last".to_owned().into()
            ))
        );

        let values = cosmolkit_cx::parse_cx_extensions_progress("|$_AV:v0;;v2$|");
        assert!(values.is_complete());
        apply_cx_progress_to_query(&mut graph, &values).expect("value effects");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("molFileValue")),
            Some(&cosmolkit_model::PropertyValue::String(
                "v0".to_owned().into()
            ))
        );
        assert_eq!(
            graph.atom(1).and_then(|atom| atom.prop("molFileValue")),
            None
        );
        assert_eq!(
            graph.atom(2).and_then(|atom| atom.prop("molFileValue")),
            Some(&cosmolkit_model::PropertyValue::String(
                "v2".to_owned().into()
            ))
        );
    }

    #[test]
    fn cx_progress_labels_apply_items_before_a_missing_closing_delimiter() {
        let mut graph = parse_smarts_graph("CC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|$first;second");
        assert!(!progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("committed label effects");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "first".to_owned().into()
            ))
        );
        assert_eq!(
            graph.atom(1).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "second".to_owned().into()
            ))
        );
    }

    #[test]
    fn query_sgroups_cx_progress_labels_replace_atom_and_preserve_existing_groups() {
        let mut graph = parse_smarts_graph("CC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let substance_group = SubstanceGroup::new(
            SubstanceGroupId::new(0),
            SubstanceGroupKind::StructuralRepeatUnit,
        )
        .with_rdkit_sequence_id(31)
        .with_atoms(vec![AtomId::new(0), AtomId::new(1), AtomId::new(0)])
        .with_bonds(vec![BondId::new(0), BondId::new(0)])
        .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
        .with_crossing_bond_correspondence(vec![BondId::new(0)])
        .with_label("existing polymer group");
        replace_query_substance_groups(&mut graph, vec![substance_group.clone()])
            .expect("typed group references are valid");
        let stereo_group = cosmolkit_model::StereoGroup::new(
            cosmolkit_model::StereoGroupKind::And,
            vec![AtomId::new(0), AtomId::new(1)],
            vec![BondId::new(0)],
        )
        .with_id(73);
        graph.add_stereo_group(stereo_group.clone());

        let progress = cosmolkit_cx::parse_cx_extensions_progress("|$Q_e;$|");
        assert!(progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("source progress label pass");

        let atom = graph.atom(0).expect("replaced query atom");
        assert_eq!(
            atom.predicate(),
            &crate::query_behavior::make_q_atom_query()
        );
        assert_eq!(
            atom.identity(),
            cosmolkit_model::QueryAtomIdentity::Element(cosmolkit_types::Element::DUMMY)
        );
        assert!(atom.no_implicit());
        assert_eq!(
            atom.prop("atomLabel"),
            Some(&cosmolkit_model::PropertyValue::String(
                "Q_e".to_owned().into()
            ))
        );
        assert_eq!(query_substance_groups(&graph), &[substance_group]);
        assert_eq!(graph.stereo_groups(), &[stereo_group]);
        assert_eq!(graph.prop("_cxsmilesLabelsProcessed"), None);
    }

    #[test]
    fn query_sgroups_cx_progress_keeps_special_label_unprocessed_until_block_completion() {
        let mut graph = parse_smarts_graph("C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let original_predicate = graph.atom(0).expect("source atom").predicate().clone();
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|$Q_e;");
        assert!(!progress.is_complete());

        apply_cx_progress_to_query(&mut graph, &progress).expect("committed label property");

        let atom = graph.atom(0).expect("partial query atom");
        assert_eq!(
            atom.prop("atomLabel"),
            Some(&cosmolkit_model::PropertyValue::String(
                "Q_e".to_owned().into()
            ))
        );
        assert_eq!(atom.predicate(), &original_predicate);
        assert_eq!(atom.element(), Some(cosmolkit_types::Element::C));
        assert_eq!(graph.prop("_cxsmilesLabelsProcessed"), None);
    }
}

#[cfg(test)]
mod cx_progress_properties_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{apply_cx_progress_to_query, parse_smarts_graph};

    #[test]
    fn cx_progress_properties_applies_items_and_skips_out_of_graph_indices() {
        let mut graph = parse_smarts_graph("CC")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(
            "|atomProp:0.label.first:1.kind.second:9.outside.skip|",
        );
        assert!(progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("atomProp effects");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("label")),
            Some(&cosmolkit_model::PropertyValue::String(
                "first".to_owned().into()
            ))
        );
        assert_eq!(
            graph.atom(1).and_then(|atom| atom.prop("kind")),
            Some(&cosmolkit_model::PropertyValue::String(
                "second".to_owned().into()
            ))
        );
    }

    #[test]
    fn cx_progress_properties_keeps_prior_item_after_later_parse_failure() {
        let mut graph = parse_smarts_graph("C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress =
            cosmolkit_cx::parse_cx_extensions_progress("|atomProp:0.kept.value:0.later.&#oops|");
        assert!(!progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("committed property effect");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("kept")),
            Some(&cosmolkit_model::PropertyValue::String(
                "value".to_owned().into()
            ))
        );
        assert_eq!(graph.atom(0).and_then(|atom| atom.prop("later")), None);
    }
}

#[cfg(test)]
mod cx_progress_bond_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{BondOrder, apply_cx_progress_to_query, parse_smarts_graph};

    #[test]
    fn cx_progress_bonds_applies_kinds_orients_pairs_and_skips_out_of_graph_indices() {
        let mut graph = parse_smarts_graph("C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");

        let coordinate = cosmolkit_cx::parse_cx_extensions_progress("|C:1.0,9.0,0.9|");
        assert!(coordinate.is_complete());
        apply_cx_progress_to_query(&mut graph, &coordinate).expect("coordinate-bond effects");
        let bond = graph.bonds_mut().first().expect("one query bond");
        assert_eq!(bond.bond().order(), BondOrder::Dative);
        assert_eq!(bond.endpoints(), (1, 0));

        let hydrogen = cosmolkit_cx::parse_cx_extensions_progress("|H:0.0|");
        apply_cx_progress_to_query(&mut graph, &hydrogen).expect("hydrogen-bond effect");
        let bond = graph.bonds_mut().first().expect("one query bond");
        assert_eq!(bond.bond().order(), BondOrder::Hydrogen);
        assert_eq!(bond.endpoints(), (0, 1));

        let zero = cosmolkit_cx::parse_cx_extensions_progress("|Z:0,9|");
        apply_cx_progress_to_query(&mut graph, &zero).expect("zero-bond effect");
        let bond = graph.bonds_mut().first().expect("one query bond");
        assert_eq!(bond.bond().order(), BondOrder::Zero);
    }

    #[test]
    fn cx_progress_bonds_keeps_prior_effects_before_syntax_or_target_failure() {
        let mut graph = parse_smarts_graph("C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let malformed = cosmolkit_cx::parse_cx_extensions_progress("|C:1.0,0x.1|");
        assert!(!malformed.is_complete());
        apply_cx_progress_to_query(&mut graph, &malformed).expect("prior pair effect");
        let bond = graph.bonds_mut().first().expect("one query bond");
        assert_eq!(bond.bond().order(), BondOrder::Dative);
        assert_eq!(bond.endpoints(), (1, 0));

        let mut graph = parse_smarts_graph("C-C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let mismatch = cosmolkit_cx::parse_cx_extensions_progress("|C:1.0,2.0|");
        assert!(mismatch.is_complete());
        assert!(apply_cx_progress_to_query(&mut graph, &mismatch).is_err());
        let bond = graph.bonds_mut().first().expect("first query bond");
        assert_eq!(bond.bond().order(), BondOrder::Dative);
        assert_eq!(bond.endpoints(), (1, 0));

        let mut graph = parse_smarts_graph("C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let zero_overflow = cosmolkit_cx::parse_cx_extensions_progress("|Z:0,4294967296|");
        assert!(!zero_overflow.is_complete());
        apply_cx_progress_to_query(&mut graph, &zero_overflow).expect("prior zero-bond effect");
        assert_eq!(graph.bonds_mut()[0].bond().order(), BondOrder::Zero);
    }
}

#[cfg(test)]
mod cx_progress_radical_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{apply_cx_progress_to_query, parse_smarts_graph};

    #[test]
    fn cx_progress_radicals_applies_electron_classes_and_skips_invalid_atoms() {
        let mut graph = parse_smarts_graph("C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|^1:0,9^5:1|");
        assert!(progress.is_complete());
        apply_cx_progress_to_query(&mut graph, &progress).expect("radical effects");
        assert_eq!(graph.atom(0).expect("first atom").radical_electrons(), 1);
        assert_eq!(graph.atom(1).expect("second atom").radical_electrons(), 3);
    }

    #[test]
    fn cx_progress_radicals_keeps_prior_effect_after_later_parse_failure() {
        let text = "|^1:0,4294967296|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert!(progress.error().is_some());
        assert_eq!(progress.consumed(), text.rfind('|').expect("closing pipe"));

        let mut graph = parse_smarts_graph("C-C")
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph");
        apply_cx_progress_to_query(&mut graph, &progress).expect("committed radical effect");
        assert_eq!(graph.atom(0).expect("first atom").radical_electrons(), 1);
        assert_eq!(graph.atom(1).expect("second atom").radical_electrons(), 0);
    }
}

#[cfg(test)]
mod cx_progress_stereo_merge_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_model::{AtomId, StereoGroup, StereoGroupKind, replace_query_stereo_groups};

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    #[test]
    fn cx_progress_stereo_merge_matches_both_consumers_and_source_order() {
        let input = "|a:0,1,0,o1:1,2,o1:2,2,&3:0,99,&3:1,a:99,o4:,o4:2,&858993459:2|";
        let parsed = cosmolkit_cx::parse_cx_extensions(input).expect("complete CX records");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph("CCO");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct record lowering");
        let mut progressed = initial.clone();
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress record lowering");

        let expected_groups = vec![
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(0),
                    AtomId::new(2),
                ],
                Vec::new(),
            )
            .with_id(858_993_459),
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(2),
                    AtomId::new(2),
                ],
                Vec::new(),
            )
            .with_id(1),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            )
            .with_id(3),
            StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(2)], Vec::new()).with_id(4),
        ];
        assert_eq!(direct.stereo_groups(), expected_groups);
        assert_eq!(progressed.stereo_groups(), expected_groups);

        let mut expected_graph = initial;
        replace_query_stereo_groups(&mut expected_graph, expected_groups)
            .expect("validated expected groups");
        expected_graph
            .set_prop(
                "__computedProps",
                cosmolkit_model::PropertyValue::StringVector(Vec::new()),
            )
            .unwrap();
        assert_eq!(direct, expected_graph);
        assert_eq!(progressed, expected_graph);
    }

    #[test]
    fn cx_progress_second_stereo_session_reports_native_tracker_invariant() {
        // CXSmilesOps.cpp clears _sgTracker after a successful parse, but
        // parse_enhanced_stereo asserts groups.size() == sgTracker.size().
        // A second session on a graph with existing groups is NOT an append API.
        let input = "|o7:0|";
        let parsed = cosmolkit_cx::parse_cx_extensions(input).unwrap();
        let initial = query_graph("CC");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).unwrap();
        let before = direct.clone();
        assert_eq!(
            crate::apply_cx_to_query_graph(&mut direct, &parsed),
            Err(crate::cx_lowering::CxQueryLoweringError::StereoGroupMerge(
                cosmolkit_core::parser_stereo_order::CxStereoGroupMergeError::TrackerCount {
                    groups: 1,
                    hashes: 0,
                }
            ))
        );
        assert_eq!(direct, before);
        let mut progressed = initial;
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        apply_cx_progress_to_query(&mut progressed, &progress).unwrap();
        let before = progressed.clone();
        let error = apply_cx_progress_to_query(&mut progressed, &progress).unwrap_err();
        assert!(
            error
                .to_string()
                .contains("source stereo group/tracker count invariant failed: 1 != 0")
        );
        assert_eq!(progressed, before);
        assert_eq!(direct, progressed);
        assert_eq!(direct.stereo_groups()[0].atoms(), &[AtomId::new(0)]);
    }

    #[test]
    fn cx_progress_stereo_keeps_completed_group_before_later_index_failure() {
        let input = "|o7:0,&3:1,4294967296|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        assert!(!progress.is_complete());
        let number_start = input.find("4294967296").expect("later index overflow");
        assert_eq!(progress.consumed(), number_start + "4294967296".len());
        assert!(progress.error().is_some());

        let initial = query_graph("CC");
        let mut graph = initial.clone();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("completed source group before failed group");

        let expected_group =
            StereoGroup::new(StereoGroupKind::Or, vec![AtomId::new(0)], Vec::new()).with_id(7);
        assert_eq!(graph.stereo_groups(), &[expected_group.clone()]);
        let mut expected = initial;
        replace_query_stereo_groups(&mut expected, vec![expected_group])
            .expect("valid completed group");
        assert_eq!(graph, expected);
    }

    #[test]
    fn cx_progress_stereo_commits_completed_helper_before_missing_pipe() {
        let input = "|&3:0,1";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), input.len());
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(input.len())
        );

        let initial = query_graph("CC");
        let mut graph = initial.clone();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("completed enhanced-stereo helper before outer delimiter failure");

        let expected_group = StereoGroup::new(
            StereoGroupKind::And,
            vec![AtomId::new(0), AtomId::new(1)],
            Vec::new(),
        )
        .with_id(3);
        let mut expected = initial;
        replace_query_stereo_groups(&mut expected, vec![expected_group])
            .expect("valid completed group");
        assert_eq!(graph, expected);
    }
}

#[cfg(test)]
mod cx_progress_constraints_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_model::{AtomQueryPredicate, QueryNode};

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    fn contains_predicate(
        node: &QueryNode<AtomQueryPredicate>,
        expected: &AtomQueryPredicate,
    ) -> bool {
        match node {
            QueryNode::Predicate(predicate) => predicate == expected,
            QueryNode::And(children) | QueryNode::Or(children) => children
                .iter()
                .any(|child| contains_predicate(child, expected)),
            _ => false,
        }
    }

    fn apply_complete(graph: &mut QueryGraph, input: &str) {
        let parsed = cosmolkit_cx::parse_cx_extensions(input).expect("complete CX record");
        crate::apply_cx_to_query_graph(graph, &parsed).expect("direct CX lowering");
    }

    #[test]
    fn cx_progress_constraints_apply_query_scan_and_skip_source_invalid_atoms() {
        let input = "|u:0,99,rb:1:4,1:*,99:3,s:0:2,1:*|";
        let parsed = cosmolkit_cx::parse_cx_extensions(input).expect("complete CX records");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph("CC");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress lowering");
        assert_eq!(progressed, direct);

        assert!(contains_predicate(
            progressed.atom(0).expect("first atom").predicate(),
            &AtomQueryPredicate::IsUnsaturated
        ));
        assert!(contains_predicate(
            progressed.atom(0).expect("first atom").predicate(),
            &AtomQueryPredicate::NonHydrogenDegree(2)
        ));
        assert!(contains_predicate(
            progressed.atom(1).expect("second atom").predicate(),
            &AtomQueryPredicate::RingBondCountLessEqual(4)
        ));
        assert!(contains_predicate(
            progressed.atom(1).expect("second atom").predicate(),
            &AtomQueryPredicate::RingBondCount(
                crate::query_behavior::QUERY_SCAN_MAGIC_VALUE as i32
            )
        ));
        assert!(contains_predicate(
            progressed.atom(1).expect("second atom").predicate(),
            &AtomQueryPredicate::NonHydrogenDegree(crate::query_behavior::QUERY_SCAN_MAGIC_VALUE)
        ));
    }

    #[test]
    fn cx_progress_constraints_keep_prior_items_before_later_failure() {
        for (input, valid_prefix) in [
            ("|u:0,4294967296|", "|u:0|"),
            ("|rb:0:3,1:1|", "|rb:0:3|"),
            ("|s:0:2,1:x|", "|s:0:2|"),
        ] {
            let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
            assert!(!progress.is_complete(), "{input:?}");
            assert!(progress.error().is_some(), "{input:?}");

            let initial = query_graph("CC");
            let mut expected = initial.clone();
            apply_complete(&mut expected, valid_prefix);
            assert_eq!(
                expected.prop("__computedProps"),
                Some(&cosmolkit_model::PropertyValue::StringVector(Vec::new()))
            );
            expected.clear_prop("__computedProps").unwrap();
            let mut progressed = initial;
            apply_cx_progress_to_query(&mut progressed, &progress)
                .expect("apply source-committed query items");
            assert_eq!(progressed, expected, "{input:?}");
        }
    }

    #[test]
    fn cx_progress_constraints_recognizes_malformed_empty_substitution_without_effect() {
        let input = "|s:0:x|";
        let failure = input.find('x').expect("malformed substitution value");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(input);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        assert_eq!(progress.error().map(|error| error.offset), Some(failure));
        assert!(matches!(
            progress.records(),
            [cosmolkit_cx::CxRecord::Substitution(constraints)] if constraints.is_empty()
        ));

        let initial = query_graph("CC");
        let mut graph = initial.clone();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("recognized incomplete substitution has no completed item");
        assert_eq!(graph, initial);
    }
}

#[cfg(test)]
mod cx_progress_linknodes_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    #[test]
    fn cx_progress_linknodes_direct_and_progress_lowering_match_source_property() {
        let text = "|LN:0:1.3.1.2,1:2.4.2.0|";
        let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete link-node record");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph("C(C)C");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress lowering");

        assert_eq!(progressed, direct);
        assert_eq!(
            direct.prop("_molLinkNodes"),
            Some(&cosmolkit_model::PropertyValue::String(
                "1 3 2 1 2 1 3|2 4 2 2 3 2 1".into()
            ))
        );
    }

    #[test]
    fn cx_progress_linknodes_does_not_commit_a_partial_helper() {
        let text = "|$left$LN:0:1.2.0.0,4294967296:1.2.0.0|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert!(progress.error().is_some());
        assert!(
            progress
                .checkpoints()
                .iter()
                .any(|checkpoint| checkpoint.phase == cosmolkit_cx::CxProgressPhase::Item)
        );
        assert!(
            progress
                .checkpoints()
                .iter()
                .all(
                    |checkpoint| checkpoint.phase != cosmolkit_cx::CxProgressPhase::Complete
                        || checkpoint.record_index == 0
                )
        );

        let mut graph = query_graph("C");
        apply_cx_progress_to_query(&mut graph, &progress).expect("source partial effects");
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "left".to_owned().into()
            ))
        );
        assert_eq!(graph.prop("_molLinkNodes"), None);
    }

    #[test]
    fn cx_progress_linknodes_commits_complete_helper_before_missing_pipe() {
        let text = "|LN:0:1.2.0.0";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert!(
            progress
                .checkpoints()
                .iter()
                .any(|checkpoint| { checkpoint.phase == cosmolkit_cx::CxProgressPhase::Complete })
        );

        let mut graph = query_graph("C");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("link-node helper completed before outer delimiter failure");
        assert_eq!(
            graph.prop("_molLinkNodes"),
            Some(&cosmolkit_model::PropertyValue::String(
                "1 2 2 1 1 1 1".into()
            ))
        );
    }

    #[test]
    fn cx_progress_linknodes_degree_error_preserves_existing_property() {
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|LN:0:1.2|");
        assert!(progress.is_complete());

        let mut graph = query_graph("CC")
            .with_prop("_molLinkNodes", "prior")
            .expect("original valid fixture property must succeed");
        assert!(apply_cx_progress_to_query(&mut graph, &progress).is_err());
        assert_eq!(
            graph.prop("_molLinkNodes"),
            Some(&cosmolkit_model::PropertyValue::String("prior".into()))
        );
    }
}

#[cfg(test)]
mod cx_progress_sgroups_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_cx::{CxRecord, CxSGroupHierarchy};
    use cosmolkit_model::{
        AtomId, QueryAtomIdentity, SubstanceGroup, SubstanceGroupId, SubstanceGroupKind,
        query_substance_groups, replace_query_substance_groups,
    };
    use cosmolkit_types::Element;

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    fn query_graph_with_existing_groups() -> QueryGraph {
        let mut graph = query_graph("CC");
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0)]),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(1)]),
        ];
        replace_query_substance_groups(&mut graph, groups).expect("valid initial groups");
        graph
    }

    fn query_graph_with_hierarchy_groups() -> QueryGraph {
        let mut graph = query_graph("CC");
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(5)
                .with_prop("_cxsmilesindex", "5")
                .expect("original valid fixture constructor must succeed")
                .with_prop("index", "71")
                .expect("original valid fixture property must succeed"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(0)
                .with_prop("_cxsmilesindex", "0")
                .expect("original valid fixture constructor must succeed")
                .with_prop("index", "92")
                .expect("original valid fixture property must succeed"),
            SubstanceGroup::new(SubstanceGroupId::new(2), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(1)
                .with_prop("_cxsmilesindex", "1")
                .expect("original valid fixture constructor must succeed")
                .with_prop("index", "93")
                .expect("original valid fixture property must succeed"),
        ];
        replace_query_substance_groups(&mut graph, groups).expect("valid hierarchy groups");
        graph
    }

    #[test]
    fn cx_progress_sgroups_direct_and_progress_lowering_keep_source_and_storage_ids() {
        let text = "|SgD:9:IGNORED:::::SgD:1,0:FIELD:value,with,comma:=:unit:tag:(1,2)|";
        let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete data SGroups");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph_with_existing_groups();
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress lowering");
        assert_eq!(progressed, direct);

        let groups = query_substance_groups(&progressed);
        assert_eq!(groups.len(), 3);
        let data_group = &groups[2];
        assert_eq!(data_group.id(), SubstanceGroupId::new(2));
        assert_eq!(data_group.rdkit_sequence_id(), Some(1));
        assert_eq!(data_group.atoms(), &[AtomId::new(1), AtomId::new(0)]);
        assert_eq!(
            data_group.props().get(b"_cxsmilesindex".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(1))
        );
        assert_eq!(
            data_group.props().get(b"index".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(3))
        );
        assert_eq!(
            data_group.props().get(b"DATAFIELDS".as_slice()),
            Some(&cosmolkit_model::PropertyValue::StringVector(vec![
                "value,with,comma".into()
            ]))
        );
        assert_eq!(
            data_group.props().get(b"COORDS".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("(1,2".into()))
        );
        assert_eq!(
            data_group.data_fields(),
            &[cosmolkit_model::PropertyText::from("value,with,comma")]
        );
        let typed_data = data_group.data().expect("typed DAT payload");
        assert_eq!(
            typed_data
                .field_name
                .as_ref()
                .map(crate::smarts_parse::fixture_text),
            Some("FIELD")
        );
        assert_eq!(
            typed_data
                .query_op
                .as_ref()
                .map(crate::smarts_parse::fixture_text),
            Some("=")
        );
        assert_eq!(
            typed_data
                .field_info
                .as_ref()
                .map(crate::smarts_parse::fixture_text),
            Some("unit")
        );
        assert_eq!(
            typed_data.values,
            [cosmolkit_model::PropertyText::from("value,with,comma")]
        );
        assert_eq!(
            typed_data
                .field_display
                .as_ref()
                .map(crate::smarts_parse::fixture_text),
            Some("    0.0000    0.0000    DR    ALL  0       0")
        );
    }

    #[test]
    fn cx_progress_sgroups_failure_before_complete_has_no_search_effect() {
        let text = "|SgD:0:FIELD:bad&#oops:QUERY:INFO:TAG:|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert!(progress.error().is_some());

        let mut graph = query_graph("C");
        graph
            .atom_mut(0)
            .expect("first atom")
            .set_prop("atomLabel", "Q_e")
            .expect("atom label property");
        let initial = graph.clone();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("incomplete local DAT state has no destination effect");

        assert_eq!(graph, initial);
        assert!(query_substance_groups(&graph).is_empty());
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "Q_e".to_owned().into()
            ))
        );
    }

    #[test]
    fn cx_progress_sgroups_attach_before_later_failure_after_source_label_processing() {
        let text = "|SgD:0:FIELD:one::::SgD:0:NEXT:bad&#oops:QUERY:INFO:TAG:|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert!(progress.error().is_some());

        let mut graph = query_graph("C");
        graph
            .atom_mut(0)
            .expect("first atom")
            .set_prop("atomLabel", "Q_e")
            .expect("atom label property");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("first completed data group commits before the later parse error");

        assert_eq!(
            graph.atom(0).expect("first atom").identity(),
            QueryAtomIdentity::Element(Element::DUMMY)
        );
        assert_eq!(
            graph.prop("_cxsmilesLabelsProcessed"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
        let groups = query_substance_groups(&graph);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].id(), SubstanceGroupId::new(0));
        assert_eq!(groups[0].rdkit_sequence_id(), Some(0));
        assert_eq!(
            groups[0].props().get(b"index".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(1))
        );
        assert_eq!(
            groups[0].data_fields(),
            &[cosmolkit_model::PropertyText::from("one")]
        );
    }

    #[test]
    fn cx_progress_sgroups_commit_completed_helper_before_missing_outer_pipe() {
        let text = "|SgD:0:FIELD:DATA:QUERY:INFO:TAG:(raw|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), text.len());
        assert!(
            progress
                .checkpoints()
                .iter()
                .any(|checkpoint| { checkpoint.phase == cosmolkit_cx::CxProgressPhase::Complete })
        );

        let mut graph = query_graph("C");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("complete source helper commits before outer delimiter failure");
        let groups = query_substance_groups(&graph);
        assert_eq!(groups.len(), 1);
        assert_eq!(
            groups[0].props().get(b"COORDS".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("(raw|".into()))
        );
    }

    #[test]
    fn cx_progress_hierarchy_direct_and_progress_paths_preserve_source_identity() {
        let text = "|SgH:5:0.1.0|";
        let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete hierarchy");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph_with_hierarchy_groups();
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct hierarchy lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress)
            .expect("progress hierarchy lowering");

        assert_eq!(progressed, direct);
        let groups = query_substance_groups(&progressed);
        assert_eq!(groups[0].parent(), None);
        assert_eq!(groups[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[1].props().get(b"PARENT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(71))
        );
        assert_eq!(groups[2].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[2].props().get(b"PARENT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(71))
        );
        assert_eq!(groups[1].rdkit_sequence_id(), Some(0));
        assert_eq!(groups[2].rdkit_sequence_id(), Some(1));
    }

    #[test]
    fn cx_progress_hierarchy_keeps_prior_mutation_after_later_syntax_error() {
        let text = "|SgH:5:0,77x:1|";
        let error_cursor = text.find('x').expect("bad later parent delimiter");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), error_cursor);
        assert_eq!(
            progress.error().map(|error| error.offset),
            Some(error_cursor)
        );
        assert!(matches!(
            progress.records(),
            [CxRecord::SGroupHierarchy(hierarchies)] if hierarchies == &[
                CxSGroupHierarchy { parent: 5, children: vec![0] },
                CxSGroupHierarchy { parent: 77, children: Vec::new() },
            ]
        ));
        assert_eq!(progress.checkpoints().len(), 4);
        assert_eq!(
            progress.checkpoints()[1].phase,
            cosmolkit_cx::CxProgressPhase::Begin
        );
        assert_eq!(progress.checkpoints()[2].item_index, Some(0));
        assert_eq!(
            progress.checkpoints()[3].phase,
            cosmolkit_cx::CxProgressPhase::Begin
        );

        let mut graph = query_graph_with_hierarchy_groups();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("prior source relationship commits before later syntax failure");

        let groups = query_substance_groups(&graph);
        assert_eq!(groups[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[1].props().get(b"PARENT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(71))
        );
        assert_eq!(groups[2].parent(), None);
    }

    #[test]
    fn cx_progress_hierarchy_keeps_prior_mutation_after_lowering_error() {
        let progress = cosmolkit_cx::parse_cx_extensions_progress("|SgH:5:0.9|");
        assert!(progress.is_complete());

        let mut graph = query_graph_with_hierarchy_groups();
        let error = apply_cx_progress_to_query(&mut graph, &progress)
            .expect_err("valid parent with out-of-range later child fails");

        assert!(
            error
                .to_string()
                .contains("child id references non-existent SGroup")
        );
        let groups = query_substance_groups(&graph);
        assert_eq!(groups[1].parent(), Some(SubstanceGroupId::new(0)));
        assert_eq!(
            groups[1].props().get(b"PARENT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(71))
        );
        assert_eq!(groups[2].parent(), None);
    }
}

#[cfg(test)]
mod cx_progress_polymer_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_model::{
        AtomId, BondId, QueryAtomIdentity, SGroupConnection, SubstanceGroup, SubstanceGroupId,
        SubstanceGroupKind, query_substance_groups, replace_query_substance_groups,
    };
    use cosmolkit_types::Element;

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    fn set_query_label(graph: &mut QueryGraph) {
        graph
            .atom_mut(0)
            .expect("first query atom")
            .set_prop("atomLabel", "Q_e")
            .expect("query atom label property");
    }

    #[test]
    fn cx_progress_polymer_direct_and_progress_preserve_typed_order_and_ids() {
        let text = "|Sg:alt:1,0,1:repeat:hh&#44;f:0,2:1:|";
        let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete polymer SGroup");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let mut initial = query_graph("CCCC");
        let existing = SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
            .with_rdkit_sequence_id(77)
            .with_atoms(vec![AtomId::new(3)]);
        replace_query_substance_groups(&mut initial, vec![existing])
            .expect("valid preexisting query SGroup");

        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct polymer lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress polymer lowering");
        assert_eq!(progressed, direct);

        let groups = query_substance_groups(&progressed);
        assert_eq!(groups.len(), 2);
        let group = &groups[1];
        assert_eq!(group.id(), SubstanceGroupId::new(1));
        assert_eq!(group.rdkit_sequence_id(), Some(0));
        assert_eq!(group.kind(), &SubstanceGroupKind::Copolymer);
        assert_eq!(
            group.subtype().map(crate::smarts_parse::fixture_text),
            Some("ALT")
        );
        assert_eq!(
            group.label().map(crate::smarts_parse::fixture_text),
            Some("repeat")
        );
        assert_eq!(group.connection(), Some(&SGroupConnection::HeadToHead));
        assert_eq!(
            group.atoms(),
            &[AtomId::new(1), AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(
            group.bonds(),
            &[BondId::new(0), BondId::new(2), BondId::new(1)]
        );
        assert_eq!(
            group.head_crossing_bonds(),
            &[BondId::new(0), BondId::new(2)]
        );
        assert_eq!(
            group.crossing_bond_correspondence(),
            &[BondId::new(0), BondId::new(1)]
        );
        assert_eq!(
            group.props().get(b"_cxsmilesindex".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(0))
        );
        assert_eq!(
            group.props().get(b"index".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(2))
        );
        assert_eq!(
            group.props().get(b"CONNECT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("HH".into()))
        );
        assert_eq!(
            group.props().get(b"SUBTYPE".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("ALT".into()))
        );
    }

    #[test]
    fn cx_progress_polymer_maps_all_source_types_and_keeps_skipped_sequence_ids() {
        let cases = [
            ("n", SubstanceGroupKind::StructuralRepeatUnit, None),
            ("mon", SubstanceGroupKind::Monomer, None),
            ("mer", SubstanceGroupKind::Mer, None),
            ("co", SubstanceGroupKind::Copolymer, None),
            ("xl", SubstanceGroupKind::Crosslink, None),
            ("mod", SubstanceGroupKind::Modification, None),
            ("mix", SubstanceGroupKind::MixtureComponent, None),
            ("f", SubstanceGroupKind::Formulation, None),
            ("any", SubstanceGroupKind::AnyPolymer, None),
            (
                "gen",
                SubstanceGroupKind::Generic("GEN".to_owned().into()),
                None,
            ),
            (
                "c",
                SubstanceGroupKind::Generic("COM".to_owned().into()),
                None,
            ),
            ("grf", SubstanceGroupKind::Graft, None),
            ("alt", SubstanceGroupKind::Copolymer, Some("ALT")),
            ("ran", SubstanceGroupKind::Copolymer, Some("RAN")),
            ("blk", SubstanceGroupKind::Copolymer, Some("BLO")),
        ];
        let mut source_records = vec!["Sg:n:99".to_owned()];
        source_records.extend(
            cases
                .iter()
                .map(|(type_code, _, _)| format!("Sg:{type_code}:0")),
        );
        let text = format!("|{}|", source_records.join(","));
        let progress = cosmolkit_cx::parse_cx_extensions_progress(&text);
        assert!(progress.is_complete());

        let mut graph = query_graph("CC");
        apply_cx_progress_to_query(&mut graph, &progress).expect("source polymer records lower");
        let groups = query_substance_groups(&graph);
        assert_eq!(groups.len(), cases.len());
        for (index, (group, (type_code, kind, subtype))) in groups.iter().zip(cases).enumerate() {
            assert_eq!(group.id(), SubstanceGroupId::new(index));
            assert_eq!(group.rdkit_sequence_id(), Some(index as u32 + 1));
            assert_eq!(group.kind(), &kind, "{type_code}");
            assert_eq!(
                group.subtype().map(crate::smarts_parse::fixture_text),
                subtype,
                "{type_code}"
            );
            assert_eq!(
                group.props().get(b"_cxsmilesindex".as_slice()),
                Some(&cosmolkit_model::PropertyValue::UInt((index + 1) as u32))
            );
            assert_eq!(
                group.props().get(b"index".as_slice()),
                Some(&cosmolkit_model::PropertyValue::UInt((index + 1) as u32))
            );
            assert_eq!(group.atoms(), &[AtomId::new(0)], "{type_code}");
        }
    }

    #[test]
    fn cx_progress_polymer_partial_syntax_has_no_search_effect() {
        let text = "|Sg:n:0:repeat:eu:1:4294967296|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(!progress.is_complete());
        assert!(progress.error().is_some());

        let mut graph = query_graph("CCC");
        set_query_label(&mut graph);
        let initial = graph.clone();
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("incomplete local polymer state has no graph effect");
        assert_eq!(graph, initial);
        assert!(query_substance_groups(&graph).is_empty());
        assert_eq!(
            graph.atom(0).and_then(|atom| atom.prop("atomLabel")),
            Some(&cosmolkit_model::PropertyValue::String(
                "Q_e".to_owned().into()
            ))
        );
    }

    #[test]
    fn cx_progress_polymer_skip_and_error_follow_source_commit_order() {
        let skipped_then_bad = "|Sg:n:0::eu:3,Sg:unknown:1|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(skipped_then_bad);
        assert!(!progress.is_complete());
        let mut skipped_graph = query_graph("CCC");
        set_query_label(&mut skipped_graph);
        let skipped_initial = skipped_graph.clone();
        apply_cx_progress_to_query(&mut skipped_graph, &progress)
            .expect("invalid crossing atom skips group before source label processing");
        assert_eq!(skipped_graph, skipped_initial);
        assert!(query_substance_groups(&skipped_graph).is_empty());

        let invalid_bond = cosmolkit_cx::parse_cx_extensions_progress("|Sg:n:0:label:eu:2|");
        assert!(invalid_bond.is_complete());
        let mut graph = query_graph("CCC");
        set_query_label(&mut graph);
        let error = apply_cx_progress_to_query(&mut graph, &invalid_bond)
            .expect_err("source-valid crossing atom reaches checked bond insertion");
        assert!(
            error
                .to_string()
                .contains("outside the topology with 2 bonds")
        );
        assert!(query_substance_groups(&graph).is_empty());
        assert_eq!(
            graph.atom(0).unwrap().identity(),
            QueryAtomIdentity::Element(Element::DUMMY)
        );
        assert_eq!(
            graph.prop("_cxsmilesLabelsProcessed"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );
    }

    #[test]
    fn cx_progress_polymer_commits_completed_group_before_later_or_outer_failure() {
        let later_error = "|Sg:n:0,Sg:unknown:1|";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(later_error);
        assert!(!progress.is_complete());
        let mut graph = query_graph("CC");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("first complete helper commits before later parser error");
        let groups = query_substance_groups(&graph);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].rdkit_sequence_id(), Some(0));
        assert_eq!(
            groups[0].props().get(b"index".as_slice()),
            Some(&cosmolkit_model::PropertyValue::UInt(1))
        );

        let missing_outer_pipe = "|Sg:n:0:repeat:hh";
        let progress = cosmolkit_cx::parse_cx_extensions_progress(missing_outer_pipe);
        assert!(!progress.is_complete());
        assert!(
            progress
                .checkpoints()
                .iter()
                .any(|checkpoint| { checkpoint.phase == cosmolkit_cx::CxProgressPhase::Complete })
        );
        let mut graph = query_graph("CC");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("completed polymer helper commits before missing outer pipe");
        let groups = query_substance_groups(&graph);
        assert_eq!(groups.len(), 1);
        assert_eq!(
            groups[0].label().map(crate::smarts_parse::fixture_text),
            Some("repeat")
        );
        assert_eq!(
            groups[0].props().get(b"CONNECT".as_slice()),
            Some(&cosmolkit_model::PropertyValue::String("HH".into()))
        );
        assert_eq!(groups[0].head_crossing_bonds(), &[BondId::new(0)]);
    }
}

#[cfg(test)]
mod cx_progress_attachments_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    fn bond_prop<'a>(graph: &'a QueryGraph, bond_index: usize, key: &str) -> Option<&'a str> {
        graph
            .bond(bond_index)
            .and_then(|bond| bond.bond().prop(key))
            .map(|value| {
                crate::smarts_parse::fixture_text(
                    value.as_string().expect("string property fixture"),
                )
            })
    }

    #[test]
    fn cx_progress_attachments_direct_and_progress_preserve_source_values() {
        let text = "|m:0:1.1.99,3:0.|";
        let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete attachments");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph("CCCC");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct attachment lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress)
            .expect("progress attachment lowering");
        assert_eq!(progressed, direct);

        assert_eq!(
            bond_prop(&progressed, 0, "_MolFileBondEndPts"),
            Some("(2 2 2)")
        );
        assert_eq!(bond_prop(&progressed, 0, "_MolFileBondAttach"), Some("ANY"));
        assert_eq!(
            bond_prop(&progressed, 2, "_MolFileBondEndPts"),
            Some("(1 1)")
        );
        assert_eq!(bond_prop(&progressed, 2, "_MolFileBondAttach"), Some("ANY"));
        assert_eq!(bond_prop(&progressed, 1, "_MolFileBondEndPts"), None);
    }

    #[test]
    fn cx_progress_attachments_skip_invalid_atoms_and_keep_empty_endpoint_count() {
        let empty_endpoints = cosmolkit_cx::parse_cx_extensions_progress("|m:0:|");
        assert!(empty_endpoints.is_complete());
        let mut graph = query_graph("CC");
        apply_cx_progress_to_query(&mut graph, &empty_endpoints)
            .expect("empty nested list still commits an empty source endpoint list");
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondEndPts"), Some("(0)"));
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondAttach"), Some("ANY"));

        let invalid_primary = cosmolkit_cx::parse_cx_extensions_progress("|m:99:0|");
        assert!(invalid_primary.is_complete());
        let mut graph = query_graph("CC");
        apply_cx_progress_to_query(&mut graph, &invalid_primary)
            .expect("source-invalid primary atom is skipped");
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondEndPts"), None);
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondAttach"), None);
    }

    #[test]
    fn cx_progress_attachments_keep_prior_rows_and_reject_later_source_degree() {
        let later_syntax_error = "|m:0:2,3x:1|";
        let error_offset = later_syntax_error.find('x').unwrap();
        let progress = cosmolkit_cx::parse_cx_extensions_progress(later_syntax_error);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), error_offset);
        let mut graph = query_graph("CCCC");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("the first complete attachment row survives later syntax failure");
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondEndPts"), Some("(1 3)"));
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondAttach"), Some("ANY"));
        assert_eq!(bond_prop(&graph, 2, "_MolFileBondEndPts"), None);

        let partial_endpoint = cosmolkit_cx::parse_cx_extensions_progress("|m:0:2.4294967296|");
        assert!(!partial_endpoint.is_complete());
        let mut graph = query_graph("CCC");
        apply_cx_progress_to_query(&mut graph, &partial_endpoint)
            .expect("an incomplete nested list does not commit its row");
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondEndPts"), None);
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondAttach"), None);

        let later_degree_error = cosmolkit_cx::parse_cx_extensions_progress("|m:0:2,1:0|");
        assert!(later_degree_error.is_complete());
        let mut graph = query_graph("CCC");
        let error = apply_cx_progress_to_query(&mut graph, &later_degree_error)
            .expect_err("source rejects a position-variation atom with degree two");
        assert!(
            error
                .to_string()
                .contains("position variation bond to atom with more than one bond")
        );
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondEndPts"), Some("(1 3)"));
        assert_eq!(bond_prop(&graph, 0, "_MolFileBondAttach"), Some("ANY"));
        assert_eq!(bond_prop(&graph, 1, "_MolFileBondEndPts"), None);
    }
}

#[cfg(test)]
mod cx_progress_directions_tests {
    fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }
    fn fixture_value(value: &cosmolkit_model::PropertyValue) -> &str {
        fixture_text(value.as_string().expect("original string fixture kind"))
    }

    use super::{QueryGraph, apply_cx_progress_to_query, parse_smarts_graph};
    use cosmolkit_model::{AtomId, PropertyValue};
    use cosmolkit_types::{BondDirection, BondStereo, ChiralTag};

    fn query_graph(smarts: &str) -> QueryGraph {
        parse_smarts_graph(smarts)
            .expect("SMARTS syntax")
            .finish()
            .expect("query graph")
    }

    fn bond_prop<'a>(graph: &'a QueryGraph, bond_index: usize, key: &str) -> Option<&'a str> {
        graph
            .bond(bond_index)
            .and_then(|bond| bond.bond().prop(key))
            .map(|value| {
                crate::smarts_parse::fixture_text(
                    value.as_string().expect("string property fixture"),
                )
            })
    }

    #[test]
    fn cx_progress_directions_direct_and_progress_preserve_effects() {
        let wedge_text = "|wU:1.0,wD:2.1|";
        let parsed = cosmolkit_cx::parse_cx_extensions(wedge_text).expect("complete wedges");
        let progress = cosmolkit_cx::parse_cx_extensions_progress(wedge_text);
        assert!(progress.is_complete());
        assert_eq!(parsed.records(), progress.records());

        let initial = query_graph("C-C-C");
        let mut direct = initial.clone();
        crate::apply_cx_to_query_graph(&mut direct, &parsed).expect("direct wedge lowering");
        let mut progressed = initial;
        apply_cx_progress_to_query(&mut progressed, &progress).expect("progress wedge lowering");
        assert_eq!(progressed, direct);
        assert_eq!(progressed.bond(0).unwrap().endpoints(), (1, 0));
        assert_eq!(progressed.bond(1).unwrap().endpoints(), (2, 1));
        assert_eq!(
            progressed.bond(0).unwrap().bond().direction(),
            BondDirection::BeginWedge
        );
        assert_eq!(
            progressed.bond(1).unwrap().bond().direction(),
            BondDirection::BeginDash
        );
        // CXSmilesOps::parse_wedged_bonds stores unsigned bond_cfg values.
        assert_eq!(
            progressed.bonds()[0].bond().prop("_MolFileBondCfg"),
            Some(&PropertyValue::UInt(1))
        );
        assert_eq!(
            progressed.bonds()[1].bond().prop("_MolFileBondCfg"),
            Some(&PropertyValue::UInt(3))
        );
        assert_eq!(
            progressed.prop("_needsDetectAtomStereo"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );

        let mut unknown = query_graph("C-C");
        unknown.atoms_mut()[1].set_chiral_tag(ChiralTag::TetrahedralCw);
        let unknown_progress = cosmolkit_cx::parse_cx_extensions_progress("|w:1.0|");
        apply_cx_progress_to_query(&mut unknown, &unknown_progress)
            .expect("unknown wedge lowering");
        assert_eq!(
            unknown.bond(0).unwrap().bond().direction(),
            BondDirection::Unknown
        );
        assert_eq!(
            unknown.atom(1).unwrap().chiral_tag(),
            ChiralTag::Unspecified
        );
        assert_eq!(
            unknown.prop("_needsDetectBondStereo"),
            Some(&cosmolkit_model::PropertyValue::Int(1))
        );

        for (text, expected) in [
            ("|ctu:1|", BondStereo::Any),
            ("|c:1|", BondStereo::Cis),
            ("|t:1|", BondStereo::Trans),
        ] {
            let parsed = cosmolkit_cx::parse_cx_extensions(text).expect("complete stereo record");
            let progress = cosmolkit_cx::parse_cx_extensions_progress(text);
            let mut initial = query_graph("FC=CF");
            initial.bonds_mut()[1]
                .bond_mut()
                .set_endpoints(AtomId::new(2), AtomId::new(1));
            let mut direct = initial.clone();
            crate::apply_cx_to_query_graph(&mut direct, &parsed)
                .expect("direct double-bond stereo lowering");
            let mut progressed = initial;
            apply_cx_progress_to_query(&mut progressed, &progress)
                .expect("progress double-bond stereo lowering");
            assert_eq!(progressed, direct);
            let bond = progressed.bond(1).unwrap().bond();
            assert_eq!(bond.stereo(), expected);
            assert_eq!(bond.stereo_atoms(), Some([AtomId::new(3), AtomId::new(0)]));
            assert_eq!(
                progressed.prop("_needsDetectBondStereo"),
                Some(&cosmolkit_model::PropertyValue::Int(1))
            );
        }
    }

    #[test]
    fn cx_progress_directions_keeps_prior_effects_on_parse_and_target_failures() {
        let malformed = "|wU:0.0,1x.1|";
        let failure = malformed.find('x').unwrap();
        let progress = cosmolkit_cx::parse_cx_extensions_progress(malformed);
        assert!(!progress.is_complete());
        assert_eq!(progress.consumed(), failure);
        let mut graph = query_graph("C-C");
        apply_cx_progress_to_query(&mut graph, &progress)
            .expect("the complete first wedge pair survives later syntax failure");
        assert_eq!(
            graph.bond(0).unwrap().bond().prop("_MolFileBondCfg"),
            Some(&PropertyValue::UInt(1))
        );
        assert_eq!(
            graph.bond(0).unwrap().bond().direction(),
            BondDirection::BeginWedge
        );

        let duplicate = cosmolkit_cx::parse_cx_extensions_progress("|wU:0.0,1.0|");
        assert!(duplicate.is_complete());
        let mut graph = query_graph("C-C");
        assert!(apply_cx_progress_to_query(&mut graph, &duplicate).is_err());
        assert_eq!(
            graph.bond(0).unwrap().bond().prop("_MolFileBondCfg"),
            Some(&PropertyValue::UInt(1))
        );

        let mut graph = query_graph("C(C)(C)C");
        assert_eq!(graph.bond(0).unwrap().endpoints(), (0, 1));
        assert_eq!(graph.bond(1).unwrap().endpoints(), (0, 2));
        assert_eq!(graph.bond(2).unwrap().endpoints(), (0, 3));
        let mismatch = cosmolkit_cx::parse_cx_extensions_progress("|wU:0.0,3.1|");
        assert!(apply_cx_progress_to_query(&mut graph, &mismatch).is_err());
        assert_eq!(
            graph.bond(0).unwrap().bond().prop("_MolFileBondCfg"),
            Some(&PropertyValue::UInt(1))
        );
        assert_eq!(graph.bond(1).unwrap().bond().prop("_MolFileBondCfg"), None);

        let stereo_overflow = cosmolkit_cx::parse_cx_extensions_progress("|ctu:1,4294967296|");
        assert!(!stereo_overflow.is_complete());
        let mut graph = query_graph("FC=CF");
        apply_cx_progress_to_query(&mut graph, &stereo_overflow)
            .expect("the completed double-bond item survives later integer overflow");
        assert_eq!(graph.bond(1).unwrap().bond().stereo(), BondStereo::Any);
        assert_eq!(
            graph.bond(1).unwrap().bond().stereo_atoms(),
            Some([AtomId::new(0), AtomId::new(3)])
        );

        let invalid = cosmolkit_cx::parse_cx_extensions_progress("|wD:99.0ctu:99|");
        let mut graph = query_graph("C-C");
        apply_cx_progress_to_query(&mut graph, &invalid)
            .expect("source-invalid atom and bond indices are skipped");
        assert_eq!(graph.bond(0).unwrap().bond().prop("_MolFileBondCfg"), None);
        assert_eq!(graph.bond(0).unwrap().bond().stereo(), BondStereo::None);

        let degree_limited = cosmolkit_cx::parse_cx_extensions_progress("|c:0|");
        let mut graph = query_graph("C=C");
        apply_cx_progress_to_query(&mut graph, &degree_limited)
            .expect("source does not set stereo without two endpoint neighbors");
        assert_eq!(graph.bond(0).unwrap().bond().stereo(), BondStereo::None);
        assert_eq!(graph.bond(0).unwrap().bond().stereo_atoms(), None);
    }
}

#[cfg(test)]
mod uint_counter_proposed_tests {
    use super::*;
    #[test]
    fn proposed_uint_source_counter_preserves_unsigned_postincrement() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        for (value, next) in [
            (0_u32, 1_u32),
            (1, 2),
            (2147483646, 2147483647),
            (2147483647, 2147483648),
            (2147483648, 2147483649),
            (4294967295, 0),
        ] {
            parser.num_bonds_parsed.set(value);
            assert_eq!(parser.next_bond_source_index(), Ok(value));
            assert_eq!(parser.num_bonds_parsed.get(), next);
        }
    }
}

#[cfg(test)]
mod uint_complete_source_condition_cells {
    use super::*;
    // FROZEN UINT CONDITION: COUNTER_0
    #[test]
    fn uint_cell_counter_0_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(0_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(0_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 1_u32);
    }
    // FROZEN UINT CONDITION: COUNTER_1
    #[test]
    fn uint_cell_counter_1_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(1_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(1_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 2_u32);
    }
    // FROZEN UINT CONDITION: COUNTER_2147483646
    #[test]
    fn uint_cell_counter_2147483646_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(2147483646_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(2147483646_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 2147483647_u32);
    }
    // FROZEN UINT CONDITION: COUNTER_2147483647
    #[test]
    fn uint_cell_counter_2147483647_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(2147483647_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(2147483647_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 2147483648_u32);
    }
    // FROZEN UINT CONDITION: COUNTER_2147483648
    #[test]
    fn uint_cell_counter_2147483648_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(2147483648_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(2147483648_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 2147483649_u32);
    }
    // FROZEN UINT CONDITION: COUNTER_4294967295
    #[test]
    fn uint_cell_counter_4294967295_smarts_parse() {
        let tokens = [];
        let parser = SmartsParser::new(&tokens, b"");
        parser.num_bonds_parsed.set(4294967295_u32);
        assert_eq!(parser.next_bond_source_index(), Ok(4294967295_u32));
        assert_eq!(parser.num_bonds_parsed.get(), 0_u32);
    }
}

#[cfg(test)]
mod original_smarts_boundary_regressions {
    use super::*;

    #[test]
    fn original_smarts_grammar_permutation_transfer() {
        // Actual pinned native yysmarts_parse/CheckChiralitySpecifications
        // observations, including right/left invalid-permutation asymmetry.
        for (text, permutation, final_tag) in [
            ("[O;@TH2]", None, ChiralTag::Tetrahedral),
            ("[C,@TH1]", None, ChiralTag::Tetrahedral),
            ("[C&@TH1]", Some(1), ChiralTag::TetrahedralCcw),
            ("[C@TH1]", Some(1), ChiralTag::TetrahedralCcw),
            ("[@TH1;C]", Some(1), ChiralTag::TetrahedralCcw),
            ("[C;@TH9]", None, ChiralTag::Tetrahedral),
            ("[C@TH1@TH2]", Some(1), ChiralTag::TetrahedralCcw),
        ] {
            let tokens = generic_parse_helper(text, ScannerStart::Atom).unwrap();
            let mut parser = SmartsParser::new(&tokens, text.as_bytes());
            let mut atom = parser.parse_atomd().unwrap().carrier;
            assert_eq!(atom.chiral_tag(), ChiralTag::Tetrahedral, "grammar {text}");
            assert_eq!(atom.chiral_permutation(), permutation, "grammar {text}");
            materialize_smarts_atom_state(&mut atom).unwrap();
            assert_eq!(atom.chiral_tag(), final_tag, "checked {text}");
            assert_eq!(atom.chiral_permutation(), None, "checked {text}");
        }
        let text = "[@TH9;C]";
        let tokens = generic_parse_helper(text, ScannerStart::Atom).unwrap();
        let mut parser = SmartsParser::new(&tokens, text.as_bytes());
        let mut atom = parser.parse_atomd().unwrap().carrier;
        assert_eq!(atom.chiral_tag(), ChiralTag::Tetrahedral);
        assert_eq!(atom.chiral_permutation(), Some(9));
        assert!(materialize_smarts_atom_state(&mut atom).is_err());
    }

    fn recursive_graph(query: &QueryNode<AtomQueryPredicate>) -> Option<&QueryGraph> {
        match query {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                recursive.query_graph()
            }
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children.iter().find_map(recursive_graph)
            }
            QueryNode::Not(child) => recursive_graph(child),
            _ => None,
        }
    }

    #[test]
    fn original_smarts_recursive_mol_keeps_raw_chirality() {
        for (text, permutation) in [
            ("[c;$([@TH1])]", 1),
            ("[O;$([@TH2])]", 2),
            ("[C;$([@TH2])]", 2),
            ("[$([@TH9])]", 9),
        ] {
            let graph = parse_smarts(text, &SmartsParseParams::default()).unwrap();
            let recursive = recursive_graph(graph.atoms()[0].predicate()).unwrap();
            assert_eq!(
                recursive.atoms()[0].chiral_tag(),
                ChiralTag::Tetrahedral,
                "{text}"
            );
            assert_eq!(
                recursive.atoms()[0].chiral_permutation(),
                Some(permutation),
                "{text}"
            );
            assert!(
                recursive.atoms()[0].prop("_SmilesStart").is_some(),
                "{text}"
            );
        }
        assert!(parse_smarts("[@TH9]", &SmartsParseParams::default()).is_err());
    }
}

#[cfg(test)]
mod smarts_source_state_tests {
    use super::*;

    #[test]
    fn querybond_normalization_replaces_composite_or_negated_predicate() {
        let previous = QueryNode::not(QueryNode::and(vec![
            make_bond_order_equals_query(BondOrder::Single),
            make_bond_is_in_ring_query(),
        ]));
        for order in [BondOrder::DativeRight, BondOrder::DativeLeft] {
            let (normalized, reverse) = normalize_dative_bond(ParsedSmartsBond {
                query: previous.clone(),
                carrier_order: order,
                unspecified_order: false,
            });
            assert_eq!(normalized.carrier_order, BondOrder::Dative);
            assert_eq!(
                normalized.query,
                make_bond_order_equals_query(BondOrder::Dative)
            );
            assert_eq!(reverse, order == BondOrder::DativeLeft);
        }
        let original = ParsedSmartsBond {
            query: previous.clone(),
            carrier_order: BondOrder::Single,
            unspecified_order: false,
        };
        let (unchanged, reverse) = normalize_dative_bond(original);
        assert_eq!(unchanged.query, previous);
        assert_eq!(unchanged.carrier_order, BondOrder::Single);
        assert!(!reverse);
    }

    #[test]
    fn ring_closure_base_setter_preserves_directional_query() {
        // Source grammar input constructed to exercise CloseMolRings' Bond*
        // dispatch; it is not a stored external counterexample or fixture.
        for (text, order, endpoints) in [
            ("C->1CCC1", BondOrder::DativeRight, (0, 3)),
            ("C1CCC<-1", BondOrder::DativeLeft, (0, 3)),
        ] {
            let query = parse_smarts_graph(text).unwrap().finish().unwrap();
            let bond = query.bond(3).unwrap();
            assert_eq!(bond.bond().order(), BondOrder::Dative);
            assert_eq!(bond.endpoints(), endpoints);
            assert_eq!(bond.predicate(), &make_bond_order_equals_query(order));
        }
    }
}

#[cfg(test)]
mod source_complete_smarts_composition_tests {
    use super::*;
    use cosmolkit_model::{BondStereo, PropertyValue, RecursiveStructureQuery};

    fn first_recursive(query: &QueryNode<AtomQueryPredicate>) -> &QueryGraph {
        match query {
            QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                recursive.query_graph().unwrap()
            }
            QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
                children
                    .iter()
                    .find_map(|child| match child {
                        QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(recursive)) => {
                            recursive.query_graph()
                        }
                        _ => None,
                    })
                    .unwrap()
            }
            QueryNode::Not(child) => first_recursive(child),
            _ => panic!("fixture has no source recursive query"),
        }
    }

    #[test]
    fn each_source_expand_query_call_retains_its_own_nested_and_node() {
        let graph = parse_smarts(
            "[C]([H])[H]",
            &SmartsParseParams {
                merge_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(graph.num_atoms(), 1);
        assert_eq!(
            graph.atoms()[0].predicate(),
            &QueryNode::And(vec![
                QueryNode::And(vec![
                    QueryNode::Predicate(AtomQueryPredicate::AtomType {
                        atomic_number: 6,
                        aromatic: false
                    }),
                    QueryNode::Not(Box::new(QueryNode::Predicate(
                        AtomQueryPredicate::HydrogenCount(0)
                    ))),
                ]),
                QueryNode::Not(Box::new(QueryNode::Predicate(
                    AtomQueryPredicate::HydrogenCount(1)
                ))),
            ])
        );
    }

    #[test]
    fn source_query_hydrogen_guard_skips_recursive_merging_for_retained_hydrogen_rows() {
        let input = "[#1;$([C][H])]";
        let raw = parse_smarts(input, &Default::default()).unwrap();
        assert_eq!(raw.num_atoms(), 1);
        assert_eq!(first_recursive(raw.atoms()[0].predicate()).num_atoms(), 2);
        let merged = parse_smarts(
            input,
            &SmartsParseParams {
                merge_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(merged.num_atoms(), 1);
        assert_eq!(
            first_recursive(merged.atoms()[0].predicate()).num_atoms(),
            2
        );
        // The native copying wrapper also quick-copies recursive query molecules.
        assert_eq!(merged, raw.clone());
    }

    #[test]
    fn source_h_token_in_an_expression_is_a_count_query_and_reaches_recursion() {
        // smarts.yy H_TOKEN action installs AtomHCount(1); HASH_TOKEN number
        // installs AtomAtomicNum. Retain the original failing input as its own
        // exact source grammar case rather than discarding its coverage.
        let input = "[H;$([C][H])]";
        let raw = parse_smarts(input, &Default::default()).unwrap();
        assert_eq!(first_recursive(raw.atoms()[0].predicate()).num_atoms(), 2);
        let merged = parse_smarts(
            input,
            &SmartsParseParams {
                merge_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        assert_eq!(
            first_recursive(merged.atoms()[0].predicate()).num_atoms(),
            1
        );
        assert!(
            matches!(raw.atoms()[0].predicate(), QueryNode::And(children)
            if children.iter().any(|node| matches!(node, QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(1)))))
        );
    }

    #[test]
    fn source_non_hydrogen_recursive_rows_reach_the_same_hydrogen_owner() {
        let input = "[C;$([C][H])]";
        let raw = parse_smarts(input, &Default::default()).unwrap();
        assert_eq!(first_recursive(raw.atoms()[0].predicate()).num_atoms(), 2);
        let merged = parse_smarts(
            input,
            &SmartsParseParams {
                merge_hs: true,
                ..Default::default()
            },
        )
        .unwrap();
        let nested = first_recursive(merged.atoms()[0].predicate());
        assert_eq!(nested.num_atoms(), 1);
        assert_eq!(nested.num_bonds(), 0);
        assert!(matches!(nested.atoms()[0].predicate(), QueryNode::And(_)));
        assert_eq!(first_recursive(raw.atoms()[0].predicate()).num_atoms(), 2);
    }

    #[test]
    fn source_direction_clear_errors_propagate_before_any_bond_changes() {
        let mut graph = parse_smarts("F/C=C/Cl", &Default::default()).unwrap();
        graph.set_prop("_needsDetectBondStereo", 1_i32).unwrap();
        graph.set_prop("__computedProps", 7_i32).unwrap();
        let before = graph.clone();
        assert!(matches!(
            crate::query_graph_behavior::set_bond_stereo_from_directions(&mut graph),
            Err(SmartsParseError::MoleculeProperty(_))
        ));
        assert_eq!(graph, before);
        assert_eq!(
            graph.prop("_needsDetectBondStereo"),
            Some(&PropertyValue::Int(1))
        );
    }

    #[test]
    fn missing_recursive_source_graph_is_a_structural_error() {
        let mut query = QueryNode::Predicate(AtomQueryPredicate::RecursiveSmarts(
            RecursiveStructureQuery::new(),
        ));
        assert_eq!(
            merge_recursive_query_hydrogens(&mut query, false, false),
            Err(SmartsParseError::MissingRecursiveQueryGraph)
        );
    }

    #[test]
    fn source_direction_postprocessing_sets_carriers_before_cis_trans_and_names_last() {
        for (input, stereo) in [
            ("F/C=C/Cl sample", BondStereo::Trans),
            ("F/C=C\\Cl sample", BondStereo::Cis),
        ] {
            let graph = parse_smarts(input, &Default::default()).unwrap();
            let bond = graph.bonds()[1].bond();
            assert_eq!(bond.stereo(), stereo);
            assert_eq!(bond.stereo_atoms(), Some([AtomId::new(0), AtomId::new(3)]));
            assert_eq!(graph.name().unwrap(), Some(&"sample".into()));
            assert_eq!(graph.prop("_needsDetectBondStereo"), None);
        }
    }
}

#[cfg(test)]
pub(crate) fn fixture_text(value: &cosmolkit_model::PropertyText) -> &str {
    std::str::from_utf8(value.as_bytes()).expect("original text fixture has valid UTF8 bytes")
}
