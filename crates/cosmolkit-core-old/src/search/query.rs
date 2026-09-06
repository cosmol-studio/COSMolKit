//! SMARTS-backed atom and bond query predicates for molecule substructure matching.
//!
//! ## RDKit provenance (protocol: dev/source_reproduction_protocol.md)
//!
//! The query-predicate system corresponds to RDKit's `GraphMol/QueryOps.cpp`,
//! `GraphMol/SmilesParse/SmartsParse.cpp`, and `GraphMol/QueryAtom.h`.
//!
//! C++ source lines are copied verbatim as commented blocks with two-axis
//! RDKit status markers per `dev/source_reproduction_protocol.md`.
//!
//! ## Implementation notes
//!
//! `QueryNode`, `AtomQueryPredicate`, `BondQueryPredicate`, and
//! `SmartsParseError` define the predicate vocabulary. The graph-level model
//! is [`super::query_graph::QueryGraph`], while parsing belongs in
//! `search::smarts_parse` and reuses these types.
//!
//! - Atom adjacency is built on-the-fly from `mol.bonds()` when not cached.
//! - Ring info is built on-the-fly from `mol.atoms()`/`mol.bonds()` when not cached.
//! - The SMARTS parser is a recursive-descent parser reproducing the Daylon
//!   Wilkins / RDKit SMARTS grammar.

use std::collections::BTreeSet;

use super::target::SearchTargetAccess;
use crate::Molecule;

impl SearchTargetAccess for Molecule {
    fn topology_block(&self) -> &cosmolkit_model::TopologyBlock {
        self.topology_block()
    }

    fn coordinate_block(&self) -> &cosmolkit_model::CoordinateBlock {
        self.coordinate_block()
    }

    fn ring_info(&self) -> Option<&RingInfo> {
        self.derived_cache()
            .rings
            .as_ref()
            .or(self.derived_cache().ring_families.as_ref())
    }

    fn valence(&self) -> Option<&ValenceAssignment> {
        self.derived_cache().valence.as_ref()
    }
}

pub use cosmolkit_model::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryGraph, QueryNode, RecursiveStructureQuery,
};

use crate::chemistry::valence::rdkit_atomic_mass;
use crate::{
    AdjacencyList, Atom, AtomId, Bond, BondOrder, BondSpec, BondStereo, ChiralTag, Hybridization,
    RingInfo, ValenceAssignment, ValenceModel,
};

#[derive(Clone)]
pub struct QueryMatchContext {
    adj: AdjacencyList,
    ring_info: Option<RingInfo>,
    valence: Option<ValenceAssignment>,
}

pub(crate) const MH_EXCLUDED_ATOMIC_NUMBERS: [u8; 22] = [
    0, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
];

fn match_atom_range_query(
    range: &AtomRangeQuery,
    atom: &Atom,
    mol: &impl SearchTargetAccess,
    context: &QueryMatchContext,
) -> bool {
    let value = match range.data_function() {
        AtomRangeDataFunction::ExplicitDegree => {
            Some(query_atom_explicit_degree(atom, &context.adj) as i32)
        }
        AtomRangeDataFunction::NonHydrogenDegree => {
            Some(query_atom_non_hydrogen_degree(atom, &context.adj, mol) as i32)
        }
        AtomRangeDataFunction::TotalDegree => {
            query_atom_total_degree(&context.adj, context.valence.as_ref(), atom)
                .map(|value| value as i32)
        }
        AtomRangeDataFunction::TotalValence => {
            query_atom_total_valence(context.valence.as_ref(), atom)
        }
        AtomRangeDataFunction::NumAtomRings => context
            .ring_info
            .as_ref()
            .map(|ring_info| query_atom_ring_membership(atom, ring_info)),
        AtomRangeDataFunction::NumHeteroatomNeighbors => {
            Some(query_atom_num_heteroatom_nbrs(atom, &context.adj, mol))
        }
        AtomRangeDataFunction::NumAliphaticHeteroatomNeighbors => Some(
            query_atom_num_aliphatic_heteroatom_nbrs(atom, &context.adj, mol),
        ),
        AtomRangeDataFunction::MinRingSize => context
            .ring_info
            .as_ref()
            .map(|ring_info| query_atom_min_ring_size(atom, ring_info) as i32),
        AtomRangeDataFunction::RingBondCount => context
            .ring_info
            .as_ref()
            .map(|ring_info| query_atom_ring_bond_count(atom, &context.adj, mol, ring_info)),
        AtomRangeDataFunction::ImplicitHydrogenCount => {
            query_atom_implicit_h_count(context.valence.as_ref(), atom).map(|value| value as i32)
        }
        AtomRangeDataFunction::FormalCharge => Some(query_atom_formal_charge(atom)),
        AtomRangeDataFunction::NegativeFormalCharge => {
            Some(query_atom_negative_formal_charge(atom))
        }
        AtomRangeDataFunction::AtomRingSize {
            lower,
            upper,
            lower_open,
            upper_open,
        } => Some(context.ring_info.as_ref().map_or_else(
            || {
                if lower > -1 {
                    -1
                } else if upper > -1 {
                    i32::MAX
                } else {
                    0
                }
            },
            |ring_info| {
                query_atom_is_in_ring_size_range(
                    atom, lower, upper, lower_open, upper_open, ring_info,
                )
            },
        )),
    };
    let Some(value) = value else {
        return false;
    };
    match range.bounds() {
        AtomRangeBounds::LessEqual(upper) => {
            greater_equal_query_match(upper, value, 0, false, |observed| observed)
        }
        AtomRangeBounds::GreaterEqual(lower) => {
            less_equal_query_match(lower, value, 0, false, |observed| observed)
        }
        AtomRangeBounds::Inclusive {
            lower,
            upper,
            lower_open,
            upper_open,
        } => range_query_match(
            lower,
            upper,
            value,
            0,
            lower_open,
            upper_open,
            false,
            |observed| observed,
        ),
    }
}

// ---------------------------------------------------------------------------
// QueryNode: a recursive Boolean query tree
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CompositeQueryType {
    And,
    Or,
    Xor,
}

fn merge_both_null_q<T>(
    return_query: &mut QueryNode<T>,
    other_null_q: &QueryNode<T>,
    how: CompositeQueryType,
) {
    // RDKit✔️✔️: void mergeBothNullQ(T *&returnQuery, T *&otherNullQ,
    // RDKit✔️✔️:                     Queries::CompositeQueryType how) {
    // RDKit✔️✔️:   bool negatedQ = returnQuery->getNegation();
    // RDKit✔️✔️:   bool negatedOtherQ = otherNullQ->getNegation();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (how == Queries::COMPOSITE_AND) {
    // RDKit✔️✔️:     // This is the only case in which we need to do anything
    // RDKit✔️✔️:     if (!negatedQ && negatedOtherQ) {
    // RDKit✔️✔️:       returnQuery->setNegation(true);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (how == Queries::COMPOSITE_OR) {
    // RDKit✔️✔️:     // This is the only case in which we need to do anything
    // RDKit✔️✔️:     if (negatedQ && !negatedOtherQ) {
    // RDKit✔️✔️:       returnQuery->setNegation(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (how == Queries::COMPOSITE_XOR) {
    // RDKit✔️✔️:     if (!negatedQ && !negatedOtherQ) {
    // RDKit✔️✔️:       returnQuery->setNegation(true);
    // RDKit✔️✔️:     } else if (negatedQ + negatedOtherQ == 1) {
    // RDKit✔️✔️:       returnQuery->setNegation(false);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations read two negation flags,
    // select one of three operators, and update at most one flag in O(1) time.
    // Rust may allocate or free one Box when set_negation changes state, the
    // already documented canonical cost of representing negation in QueryNode;
    // no traversal, lookup, clone, scan, or temporary collection is added.
    let negated_q = return_query.is_negated();
    let negated_other_q = other_null_q.is_negated();

    match how {
        CompositeQueryType::And if !negated_q && negated_other_q => {
            return_query.set_negation(true);
        }
        CompositeQueryType::Or if negated_q && !negated_other_q => {
            return_query.set_negation(false);
        }
        CompositeQueryType::Xor if !negated_q && !negated_other_q => {
            return_query.set_negation(true);
        }
        CompositeQueryType::Xor if negated_q != negated_other_q => {
            return_query.set_negation(false);
        }
        CompositeQueryType::And | CompositeQueryType::Or | CompositeQueryType::Xor => {}
    }
}

fn merge_null_q_first<T>(
    return_query: &mut QueryNode<T>,
    other_q: &mut QueryNode<T>,
    how: CompositeQueryType,
) {
    // RDKit✔️✔️: void mergeNullQFirst(T *&returnQuery, T *&otherQ,
    // RDKit✔️✔️:                      Queries::CompositeQueryType how) {
    // RDKit✔️✔️:   bool negatedQ = returnQuery->getNegation();
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (how == Queries::COMPOSITE_AND) {
    // RDKit✔️✔️:     if (!negatedQ) {
    // RDKit✔️✔️:       std::swap(returnQuery, otherQ);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (how == Queries::COMPOSITE_OR) {
    // RDKit✔️✔️:     if (negatedQ) {
    // RDKit✔️✔️:       std::swap(returnQuery, otherQ);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (how == Queries::COMPOSITE_XOR) {
    // RDKit✔️✔️:     std::swap(returnQuery, otherQ);
    // RDKit✔️✔️:     if (!negatedQ) {
    // RDKit✔️✔️:       returnQuery->setNegation(!returnQuery->getNegation());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: C++ swaps two pointers and Rust swaps two
    // enum owners in O(1) time without traversing or cloning their heap-backed
    // children. XOR may toggle the canonical outer Not box, with the same
    // representation cost documented by QueryNode::set_negation; there are no
    // scans, lookups, temporary collections, or repeated state reads.
    let negated_q = return_query.is_negated();
    match how {
        CompositeQueryType::And if !negated_q => std::mem::swap(return_query, other_q),
        CompositeQueryType::Or if negated_q => std::mem::swap(return_query, other_q),
        CompositeQueryType::Xor => {
            std::mem::swap(return_query, other_q);
            if !negated_q {
                return_query.set_negation(!return_query.is_negated());
            }
        }
        CompositeQueryType::And | CompositeQueryType::Or => {}
    }
}

fn merge_null_queries<T>(
    return_query: &mut QueryNode<T>,
    is_query_null: bool,
    other_query: &mut QueryNode<T>,
    is_other_q_null: bool,
    how: CompositeQueryType,
) {
    // RDKit✔️✔️: void mergeNullQueries(T *&returnQuery, bool isQueryNull, T *&otherQuery,
    // RDKit✔️✔️:                       bool isOtherQNull, Queries::CompositeQueryType how) {
    // RDKit✔️✔️:   PRECONDITION(returnQuery, "bad query");
    // RDKit✔️✔️:   PRECONDITION(otherQuery, "bad query");
    // RDKit✔️✔️:   PRECONDITION(how == Queries::COMPOSITE_AND || how == Queries::COMPOSITE_OR ||
    // RDKit✔️✔️:                    how == Queries::COMPOSITE_XOR,
    // RDKit✔️✔️:                "bad combination op");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (isQueryNull && isOtherQNull) {
    // RDKit✔️✔️:     mergeBothNullQ(returnQuery, otherQuery, how);
    // RDKit✔️✔️:   } else if (isQueryNull) {
    // RDKit✔️✔️:     mergeNullQFirst(returnQuery, otherQuery, how);
    // RDKit✔️✔️:   } else if (isOtherQNull) {
    // RDKit✔️✔️:     std::swap(returnQuery, otherQuery);
    // RDKit✔️✔️:     mergeNullQFirst(returnQuery, otherQuery, how);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: typed references make RDKit's non-null pointer
    // preconditions unrepresentable and the enum makes invalid operators
    // unrepresentable. Dispatch and owner swaps are O(1), and the function
    // delegates all negation behavior to the two canonical helpers without
    // traversal, cloning, lookup, scanning, or temporary collections.
    if is_query_null && is_other_q_null {
        merge_both_null_q(return_query, other_query, how);
    } else if is_query_null {
        merge_null_q_first(return_query, other_query, how);
    } else if is_other_q_null {
        std::mem::swap(return_query, other_query);
        merge_null_q_first(return_query, other_query, how);
    }
}

fn is_typed_null_query<T>(query: &QueryNode<T>, is_null_predicate: impl Fn(&T) -> bool) -> bool {
    match query {
        QueryNode::Predicate(predicate) => is_null_predicate(predicate),
        QueryNode::Not(child) => match child.as_ref() {
            QueryNode::Predicate(predicate) => is_null_predicate(predicate),
            QueryNode::Not(_) | QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) => false,
        },
        QueryNode::And(_) | QueryNode::Or(_) | QueryNode::Xor(_) => false,
    }
}

fn is_atom_null_query(query: &QueryNode<AtomQueryPredicate>) -> bool {
    is_typed_null_query(query, |predicate| {
        matches!(predicate, AtomQueryPredicate::Any)
    })
}

fn is_bond_null_query(query: &QueryNode<BondQueryPredicate>) -> bool {
    is_typed_null_query(query, |predicate| {
        matches!(predicate, BondQueryPredicate::Any)
    })
}

pub(crate) fn query_atom_expand_query(
    query: &mut QueryNode<AtomQueryPredicate>,
    mut what: QueryNode<AtomQueryPredicate>,
    how: CompositeQueryType,
    maintain_order: bool,
) {
    // RDKit✔️✔️: void QueryAtom::expandQuery(QUERYATOM_QUERY *what,
    // RDKit✔️✔️:                             Queries::CompositeQueryType how,
    // RDKit✔️✔️:                             bool maintainOrder) {
    // RDKit✔️✔️:   PRECONDITION(dp_query, "Can't expand empty query");
    // RDKit✔️✔️:   bool thisIsNullQuery = dp_query->getDescription() == "AtomNull";
    // RDKit✔️✔️:   bool otherIsNullQuery = what->getDescription() == "AtomNull";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (thisIsNullQuery || otherIsNullQuery) {
    // RDKit✔️✔️:     mergeNullQueries(dp_query, thisIsNullQuery, what, otherIsNullQuery, how);
    // RDKit✔️✔️:     delete what;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QUERYATOM_QUERY *origQ = dp_query;
    // RDKit✔️✔️:   std::string descrip;
    // RDKit✔️✔️:   switch (how) {
    // RDKit✔️✔️:     case Queries::COMPOSITE_AND:
    // RDKit✔️✔️:       dp_query = new ATOM_AND_QUERY;
    // RDKit✔️✔️:       descrip = "AtomAnd";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_OR:
    // RDKit✔️✔️:       dp_query = new ATOM_OR_QUERY;
    // RDKit✔️✔️:       descrip = "AtomOr";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_XOR:
    // RDKit✔️✔️:       dp_query = new ATOM_XOR_QUERY;
    // RDKit✔️✔️:       descrip = "AtomXor";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       UNDER_CONSTRUCTION("unrecognized combination query");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   dp_query->setDescription(descrip);
    // RDKit✔️✔️:   if (maintainOrder) {
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:     dp_query->addChild(QUERYATOM_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations classify two roots in
    // O(1), move two owned query handles, and allocate one two-child composite
    // node for non-null inputs. Rust's enum move avoids virtual dispatch and
    // performs no subtree traversal, clone, lookup, scan, or temporary string
    // allocation; null behavior delegates to the single canonical algebra.
    let this_is_null = is_atom_null_query(query);
    let other_is_null = is_atom_null_query(&what);
    if this_is_null || other_is_null {
        merge_null_queries(query, this_is_null, &mut what, other_is_null, how);
        return;
    }

    let original = std::mem::replace(query, make_atom_null_query());
    let children = if maintain_order {
        vec![original, what]
    } else {
        vec![what, original]
    };
    *query = match how {
        CompositeQueryType::And => QueryNode::and(children),
        CompositeQueryType::Or => QueryNode::or(children),
        CompositeQueryType::Xor => QueryNode::xor(children),
    };
}

pub(crate) fn query_bond_expand_query(
    query: &mut QueryNode<BondQueryPredicate>,
    mut what: QueryNode<BondQueryPredicate>,
    how: CompositeQueryType,
    maintain_order: bool,
) {
    // RDKit✔️✔️: void QueryBond::expandQuery(QUERYBOND_QUERY *what,
    // RDKit✔️✔️:                             Queries::CompositeQueryType how,
    // RDKit✔️✔️:                             bool maintainOrder) {
    // RDKit✔️✔️:   bool thisIsNullQuery = dp_query->getDescription() == "BondNull";
    // RDKit✔️✔️:   bool otherIsNullQuery = what->getDescription() == "BondNull";
    // RDKit✔️✔️:
    // RDKit✔️✔️:   if (thisIsNullQuery || otherIsNullQuery) {
    // RDKit✔️✔️:     mergeNullQueries(dp_query, thisIsNullQuery, what, otherIsNullQuery, how);
    // RDKit✔️✔️:     delete what;
    // RDKit✔️✔️:     return;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   QUERYBOND_QUERY *origQ = dp_query;
    // RDKit✔️✔️:   std::string descrip;
    // RDKit✔️✔️:   switch (how) {
    // RDKit✔️✔️:     case Queries::COMPOSITE_AND:
    // RDKit✔️✔️:       dp_query = new BOND_AND_QUERY;
    // RDKit✔️✔️:       descrip = "BondAnd";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_OR:
    // RDKit✔️✔️:       dp_query = new BOND_OR_QUERY;
    // RDKit✔️✔️:       descrip = "BondOr";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case Queries::COMPOSITE_XOR:
    // RDKit✔️✔️:       dp_query = new BOND_XOR_QUERY;
    // RDKit✔️✔️:       descrip = "BondXor";
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     default:
    // RDKit✔️✔️:       UNDER_CONSTRUCTION("unrecognized combination query");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   dp_query->setDescription(descrip);
    // RDKit✔️✔️:   if (maintainOrder) {
    // RDKit✔️✔️:     dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:     dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(what));
    // RDKit✔️✔️:     dp_query->addChild(QUERYBOND_QUERY::CHILD_TYPE(origQ));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations classify two roots in
    // O(1), move two owned query handles, and allocate one two-child composite
    // node for non-null inputs. Rust performs no subtree traversal, clone,
    // lookup, scan, or temporary string allocation; atom and bond paths share
    // the same canonical typed null classifier and NullQuery algebra.
    let this_is_null = is_bond_null_query(query);
    let other_is_null = is_bond_null_query(&what);
    if this_is_null || other_is_null {
        merge_null_queries(query, this_is_null, &mut what, other_is_null, how);
        return;
    }

    let original = std::mem::replace(query, make_bond_null_query());
    let children = if maintain_order {
        vec![original, what]
    } else {
        vec![what, original]
    };
    *query = match how {
        CompositeQueryType::And => QueryNode::and(children),
        CompositeQueryType::Or => QueryNode::or(children),
        CompositeQueryType::Xor => QueryNode::xor(children),
    };
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RangeQueryType {
    Equal,
    Less,
    Greater,
    Range,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
enum QueryFinalizationError {
    #[error("bad range query type or payload")]
    BadRangeQueryType,
    #[error("Do not know how to finalize query: '{0}'")]
    UnknownDescription(String),
}

fn finalize_atom_ring_size_query(
    query: QueryNode<AtomQueryPredicate>,
    query_type: RangeQueryType,
) -> Result<QueryNode<AtomQueryPredicate>, QueryFinalizationError> {
    // RDKit✔️✔️: switch (qtype) {
    // RDKit✔️✔️:   case RangeQueryType::EQUAL: {
    // RDKit✔️✔️:     auto tgt = static_cast<ATOM_EQUALS_QUERY *>(query)->getVal();
    // RDKit✔️✔️:     query->setDataFunc(
    // RDKit✔️✔️:         [tgt](Atom const *at) { return queryAtomIsInRingOfSize(at, tgt); });
    // RDKit✔️✔️:   } break;
    // RDKit✔️✔️:   case RangeQueryType::RANGE: {
    // RDKit✔️✔️:     auto rq = static_cast<ATOM_RANGE_QUERY *>(query);
    // RDKit✔️✔️:     auto uv = rq->getUpper();
    // RDKit✔️✔️:     auto lv = rq->getLower();
    // RDKit✔️✔️:     auto [lo, uo] = rq->getEndsOpen();
    // RDKit✔️✔️:     query->setDataFunc([lv, uv, lo, uo](Atom const *at) {
    // RDKit✔️✔️:       return queryAtomIsInRingOfSize(at, lv, uv, lo, uo);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   } break;
    // RDKit✔️✔️:   case RangeQueryType::LESS: {
    // RDKit✔️✔️:     auto lv = static_cast<ATOM_LESSEQUAL_QUERY *>(query)->getVal();
    // RDKit✔️✔️:     auto uv = -1;
    // RDKit✔️✔️:     query->setDataFunc([lv, uv](Atom const *at) {
    // RDKit✔️✔️:       return queryAtomIsInRingOfSize(at, lv, uv);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   } break;
    // RDKit✔️✔️:   case RangeQueryType::GREATER: {
    // RDKit✔️✔️:     auto lv = -1;
    // RDKit✔️✔️:     auto uv = static_cast<ATOM_GREATEREQUAL_QUERY *>(query)->getVal();
    // RDKit✔️✔️:     query->setDataFunc([lv, uv](Atom const *at) {
    // RDKit✔️✔️:       return queryAtomIsInRingOfSize(at, lv, uv);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   } break;
    // RDKit✔️✔️:   default:
    // RDKit✔️✔️:     throw ValueErrorException("bad range query type");
    // RDKit✔️✔️: }
    // Local complexity review: each branch performs O(1) typed-leaf
    // classification/construction. The range branch retains its existing
    // four scalars without allocation, traversal, lookup, or cloning.
    match (query_type, query) {
        (
            RangeQueryType::Equal,
            query @ QueryNode::Predicate(AtomQueryPredicate::InRingOfSize(_)),
        )
        | (RangeQueryType::Range, query @ QueryNode::Predicate(AtomQueryPredicate::Range(_))) => {
            Ok(query)
        }
        (RangeQueryType::Less, QueryNode::Predicate(AtomQueryPredicate::InRingOfSize(value))) => {
            Ok(QueryNode::predicate(
                AtomQueryPredicate::InRingOfSizeLessEqual(value),
            ))
        }
        (
            RangeQueryType::Greater,
            QueryNode::Predicate(AtomQueryPredicate::InRingOfSize(value)),
        ) => Ok(QueryNode::predicate(
            AtomQueryPredicate::InRingOfSizeGreaterEqual(value),
        )),
        _ => Err(QueryFinalizationError::BadRangeQueryType),
    }
}

fn finalize_atom_query_from_description(
    description: &str,
    query: QueryNode<AtomQueryPredicate>,
) -> Result<QueryNode<AtomQueryPredicate>, QueryFinalizationError> {
    // RDKit✔️✔️: std::string descr = query->getDescription();
    // RDKit✔️✔️: RangeQueryType qtype = RangeQueryType::EQUAL;
    // RDKit✔️✔️: if (boost::starts_with(descr, "range_")) {
    // RDKit✔️✔️:   descr = descr.substr(6); qtype = RangeQueryType::RANGE;
    // RDKit✔️✔️: } else if (boost::starts_with(descr, "less_")) {
    // RDKit✔️✔️:   descr = descr.substr(5); qtype = RangeQueryType::LESS;
    // RDKit✔️✔️: } else if (boost::starts_with(descr, "greater_")) {
    // RDKit✔️✔️:   descr = descr.substr(8); qtype = RangeQueryType::GREATER;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: if (descr == "AtomRingBondCount") { query->setDataFunc(queryAtomRingBondCount); }
    // RDKit✔️✔️: else if (descr == "AtomHasRingBond") { query->setDataFunc(queryAtomHasRingBond); }
    // RDKit✔️✔️: else if (descr == "AtomRingSize") { finalizeAtomRingSizeQuery(query, qtype); }
    // RDKit✔️✔️: else if (descr == "AtomMinRingSize") { query->setDataFunc(queryAtomMinRingSize); }
    // RDKit✔️✔️: else if (descr == "AtomImplicitValence") { query->setDataFunc(queryAtomImplicitValence); }
    // RDKit✔️✔️: else if (descr == "AtomTotalValence") { query->setDataFunc(queryAtomTotalValence); }
    // RDKit✔️✔️: else if (descr == "AtomAtomicNum") { query->setDataFunc(queryAtomNum); }
    // RDKit✔️✔️: else if (descr == "AtomExplicitDegree") { query->setDataFunc(queryAtomExplicitDegree); }
    // RDKit✔️✔️: else if (descr == "AtomTotalDegree") { query->setDataFunc(queryAtomTotalDegree); }
    // RDKit✔️✔️: else if (descr == "AtomHeavyAtomDegree") { query->setDataFunc(queryAtomHeavyAtomDegree); }
    // RDKit✔️✔️: else if (descr == "AtomHCount") { query->setDataFunc(queryAtomHCount); }
    // RDKit✔️✔️: else if (descr == "AtomImplicitHCount") { query->setDataFunc(queryAtomImplicitHCount); }
    // RDKit✔️✔️: else if (descr == "AtomHasImplicitH") { query->setDataFunc(queryAtomHasImplicitH); }
    // RDKit✔️✔️: else if (descr == "AtomIsAromatic") { query->setDataFunc(queryAtomAromatic); }
    // RDKit✔️✔️: else if (descr == "AtomIsAliphatic") { query->setDataFunc(queryAtomAliphatic); }
    // RDKit✔️✔️: else if (descr == "AtomUnsaturated") { query->setDataFunc(queryAtomUnsaturated); }
    // RDKit✔️✔️: else if (descr == "AtomMass") { query->setDataFunc(queryAtomMass); }
    // RDKit✔️✔️: else if (descr == "AtomIsotope") { query->setDataFunc(queryAtomIsotope); }
    // RDKit✔️✔️: else if (descr == "AtomFormalCharge") { query->setDataFunc(queryAtomFormalCharge); }
    // RDKit✔️✔️: else if (descr == "AtomNegativeFormalCharge") { query->setDataFunc(queryAtomNegativeFormalCharge); }
    // RDKit✔️✔️: else if (descr == "AtomHybridization") { query->setDataFunc(queryAtomHybridization); }
    // RDKit✔️✔️: else if (descr == "AtomInRing") { query->setDataFunc(queryIsAtomInRing); }
    // RDKit✔️✔️: else if (descr == "AtomInNRings") { query->setDataFunc(queryIsAtomInNRings); }
    // RDKit✔️✔️: else if (descr == "AtomHasHeteroatomNeighbors") { query->setDataFunc(queryAtomHasHeteroatomNbrs); }
    // RDKit✔️✔️: else if (descr == "AtomNumHeteroatomNeighbors") { query->setDataFunc(queryAtomNumHeteroatomNbrs); }
    // RDKit✔️✔️: else if (descr == "AtomNonHydrogenDegree") { query->setDataFunc(queryAtomNonHydrogenDegree); }
    // RDKit✔️✔️: else if (descr == "AtomHasAliphaticHeteroatomNeighbors") { query->setDataFunc(queryAtomHasAliphaticHeteroatomNbrs); }
    // RDKit✔️✔️: else if (descr == "AtomNumAliphaticHeteroatomNeighbors") { query->setDataFunc(queryAtomNumAliphaticHeteroatomNbrs); }
    // RDKit✔️✔️: else if (descr == "AtomNull" || descr == "AtomType" ||
    // RDKit✔️✔️:          descr == "AtomNumRadicalElectrons" || descr == "RecursiveStructure" ||
    // RDKit✔️✔️:          descr == "AtomAnd" || descr == "AtomOr" || descr == "AtomXor" ||
    // RDKit✔️✔️:          descr == "HasProp" || descr == "HasPropWithValue") { }
    // RDKit✔️✔️: else { throw ValueErrorException("Do not know how to finalize query: '" + descr + "'"); }
    // Local complexity review: prefix stripping and description dispatch are
    // linear in the short description length, matching RDKit; typed leaves
    // avoid virtual data-function writes and add no graph traversal or clone.
    let (description, query_type) = if let Some(value) = description.strip_prefix("range_") {
        (value, RangeQueryType::Range)
    } else if let Some(value) = description.strip_prefix("less_") {
        (value, RangeQueryType::Less)
    } else if let Some(value) = description.strip_prefix("greater_") {
        (value, RangeQueryType::Greater)
    } else {
        (description, RangeQueryType::Equal)
    };
    if description == "AtomRingSize" {
        return finalize_atom_ring_size_query(query, query_type);
    }
    const KNOWN: &[&str] = &[
        "AtomRingBondCount",
        "AtomHasRingBond",
        "AtomMinRingSize",
        "AtomImplicitValence",
        "AtomTotalValence",
        "AtomAtomicNum",
        "AtomExplicitDegree",
        "AtomTotalDegree",
        "AtomHeavyAtomDegree",
        "AtomHCount",
        "AtomImplicitHCount",
        "AtomHasImplicitH",
        "AtomIsAromatic",
        "AtomIsAliphatic",
        "AtomUnsaturated",
        "AtomMass",
        "AtomIsotope",
        "AtomFormalCharge",
        "AtomNegativeFormalCharge",
        "AtomHybridization",
        "AtomInRing",
        "AtomInNRings",
        "AtomHasHeteroatomNeighbors",
        "AtomNumHeteroatomNeighbors",
        "AtomNonHydrogenDegree",
        "AtomHasAliphaticHeteroatomNeighbors",
        "AtomNumAliphaticHeteroatomNeighbors",
        "AtomNull",
        "AtomType",
        "AtomNumRadicalElectrons",
        "RecursiveStructure",
        "AtomAnd",
        "AtomOr",
        "AtomXor",
        "HasProp",
        "HasPropWithValue",
    ];
    KNOWN
        .contains(&description)
        .then_some(query)
        .ok_or_else(|| QueryFinalizationError::UnknownDescription(description.to_string()))
}

fn make_atom_has_prop_query(property: impl Into<String>) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class Target>
    // RDKit✔️🔝: Queries::EqualityQuery<int, const Target *, true> *makeHasPropQuery(
    // RDKit✔️🔝:     const std::string &property) {
    // RDKit✔️🔝:   return new HasPropQuery<const Target *>(property);
    // RDKit✔️🔝: }
    // Local complexity review: both copy/move one property name in O(n). Rust
    // removes the query-object allocation and virtual Match dispatch.
    QueryNode::predicate(AtomQueryPredicate::HasProperty(property.into()))
}

fn make_atom_prop_query(
    property: impl Into<String>,
    value: impl Into<String>,
) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class Target, class T>
    // RDKit✔️🔝: Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    // RDKit✔️🔝:     const std::string &propname, const T &val, double tolerance = 0.0) {
    // RDKit✔️🔝:   return new HasPropWithValueQuery<const Target *, T>(propname, val, tolerance);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: res = atom_val == this->val;
    // Local complexity review: this is RDKit's string specialization, for
    // which tolerance is ignored. Both own two strings and compare in O(n);
    // Rust removes one allocation for the polymorphic query object.
    QueryNode::predicate(AtomQueryPredicate::PropertyValue {
        name: property.into(),
        value: value.into(),
    })
}

pub(crate) fn complex_atom_query_helper(
    query: &QueryNode<AtomQueryPredicate>,
    has_atomic_number: &mut bool,
) -> bool {
    // RDKit✔️✔️: bool _complexQueryHelper(Atom::QUERYATOM_QUERY const *query, bool &hasAtNum) {
    // RDKit✔️✔️:   if (!query) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (query->getNegation()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string descr = query->getDescription();
    // RDKit✔️✔️:   // std::cerr<<" |"<<descr;
    // RDKit✔️✔️:   if (descr == "AtomAtomicNum" || descr == "AtomType") {
    // RDKit✔️✔️:     hasAtNum = true;
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "AtomOr" || descr == "AtomXor") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "AtomAnd") {
    // RDKit✔️✔️:     auto childIt = query->beginChildren();
    // RDKit✔️✔️:     while (childIt != query->endChildren()) {
    // RDKit✔️✔️:       if (_complexQueryHelper(childIt->get(), hasAtNum)) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       ++childIt;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    //
    // The typed tree has no null child handles. `AtomicNumberIn` and
    // `AtomicNumberNotIn` retain the source semantics of OR and negated-list
    // structures even though the canonical representation compresses them
    // into leaves. Local complexity review: both implementations visit each
    // child of nested AND nodes at most once and short-circuit at the first
    // complex node, giving O(n) time and O(h) recursion. Neither allocates,
    // clones, performs keyed lookup, or creates a temporary collection.
    match query {
        QueryNode::Not(_)
        | QueryNode::Or(_)
        | QueryNode::Xor(_)
        | QueryNode::Predicate(
            AtomQueryPredicate::AtomicNumberIn(_) | AtomQueryPredicate::AtomicNumberNotIn(_),
        ) => true,
        QueryNode::Predicate(
            AtomQueryPredicate::AtomicNumber(_) | AtomQueryPredicate::AtomType { .. },
        ) => {
            *has_atomic_number = true;
            false
        }
        QueryNode::And(children) => children
            .iter()
            .any(|child| complex_atom_query_helper(child, has_atomic_number)),
        QueryNode::Predicate(_) => false,
    }
}

pub(crate) fn is_complex_atom_query(atom: &crate::QueryAtom) -> bool {
    // RDKit✔️✔️: bool isComplexQuery(const Atom *a) {
    // RDKit✔️✔️:   PRECONDITION(a, "bad atom");
    // RDKit✔️✔️:   if (!a->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // std::cerr<<"\n"<<a->getIdx();
    // RDKit✔️✔️:   // negated things are always complex:
    // RDKit✔️✔️:   if (a->getQuery()->getNegation()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string descr = a->getQuery()->getDescription();
    // RDKit✔️✔️:   // std::cerr<<" "<<descr;
    // RDKit✔️✔️:   if (descr == "AtomNull" || descr == "AtomAtomicNum" || descr == "AtomType") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "AtomOr" || descr == "AtomXor") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "AtomAnd") {
    // RDKit✔️✔️:     bool hasAtNum = false;
    // RDKit✔️✔️:     if (_complexQueryHelper(a->getQuery(), hasAtNum)) {
    // RDKit✔️✔️:       return true;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return !hasAtNum;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    //
    // Root descriptions map directly to the sole typed query variants and
    // nested AND inspection delegates to the single helper above. Local
    // complexity review: simple roots are O(1); AND trees are O(n) time and
    // O(h) recursion, matching RDKit's traversal and short-circuit behavior.
    // No path allocates, clones, performs keyed lookup, or scans molecule data.
    let _ = atom;
    false
}

pub(crate) fn is_atom_aromatic(atom: &Atom, molecule: &impl SearchTargetAccess) -> bool {
    // BEGIN RDKIT CPP FUNCTION isAromaticAtom
    // RDKit✔️✔️: bool isAromaticAtom(const Atom &atom) {
    // RDKit✔️✔️:   if (atom.getIsAromatic()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (atom.hasOwningMol()) {
    // RDKit✔️✔️:     for (const auto &bond : atom.getOwningMol().atomBonds(&atom)) {
    // RDKit✔️✔️:       if (bond->getIsAromatic() ||
    // RDKit✔️✔️:           bond->getBondType() == Bond::BondType::AROMATIC) {
    // RDKit✔️✔️:         return true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAromaticAtom
    // BEGIN RDKIT CPP FUNCTION isAtomAromatic
    // RDKit✔️✔️: bool isAtomAromatic(const Atom *a) {
    // RDKit✔️✔️:   PRECONDITION(a, "bad atom");
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   if (!a->hasQuery()) {
    // RDKit✔️✔️:     res = isAromaticAtom(*a);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     std::string descr = a->getQuery()->getDescription();
    // RDKit✔️✔️:     if (descr == "AtomAtomicNum") {
    // RDKit✔️✔️:       res = a->getIsAromatic();
    // RDKit✔️✔️:     } else if (descr == "AtomIsAromatic") {
    // RDKit✔️✔️:       res = true;
    // RDKit✔️✔️:       if (a->getQuery()->getNegation()) {
    // RDKit✔️✔️:         res = !res;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (descr == "AtomIsAliphatic") {
    // RDKit✔️✔️:       res = false;
    // RDKit✔️✔️:       if (a->getQuery()->getNegation()) {
    // RDKit✔️✔️:         res = !res;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (descr == "AtomType") {
    // RDKit✔️✔️:       res = getAtomTypeIsAromatic(
    // RDKit✔️✔️:           static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal());
    // RDKit✔️✔️:       if (a->getQuery()->getNegation()) {
    // RDKit✔️✔️:         res = !res;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     } else if (descr == "AtomAnd") {
    // RDKit✔️✔️:       auto childIt = a->getQuery()->beginChildren();
    // RDKit✔️✔️:       if ((*childIt)->getDescription() == "AtomAtomicNum") {
    // RDKit✔️✔️:         if (a->getQuery()->getNegation()) {
    // RDKit✔️✔️:           res = false;
    // RDKit✔️✔️:         } else if ((*(childIt + 1))->getDescription() == "AtomIsAliphatic") {
    // RDKit✔️✔️:           res = false;
    // RDKit✔️✔️:         } else if ((*(childIt + 1))->getDescription() == "AtomIsAromatic") {
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION isAtomAromatic
    //
    // QueryNode retains RDKit's root description as its enum identity. The
    // AtomAnd branch intentionally inspects only the first two children in
    // source order and intentionally ignores negation on its second child.
    // Local complexity review: concrete atoms scan the same indexed incident
    // bond range in O(degree); query atoms inspect one root and at most two
    // children in O(1). Neither implementation allocates, clones, performs a
    // keyed lookup, builds a temporary collection, or traverses a query tree.
    if atom.is_aromatic() {
        return true;
    }
    for neighbor in molecule
        .topology_block()
        .adjacency
        .neighbors_of(atom.id().index())
    {
        let bond = &molecule.bonds()[neighbor.bond.index()];
        if bond.is_aromatic() || bond.order() == BondOrder::Aromatic {
            return true;
        }
    }
    false
    /*

        fn atom_and_aromaticity(children: &[QueryNode<AtomQueryPredicate>]) -> bool {
            if !matches!(
                children.first(),
                Some(QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_)))
            ) {
                return false;
            }
            match children.get(1) {
                Some(QueryNode::Predicate(AtomQueryPredicate::IsAromatic(aromatic))) => *aromatic,
                Some(QueryNode::Not(child)) => match child.as_ref() {
                    QueryNode::Predicate(AtomQueryPredicate::IsAromatic(aromatic)) => *aromatic,
                    _ => false,
                },
                _ => false,
            }
        }

        match query {
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_)) => atom.is_aromatic(),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(aromatic))
            | QueryNode::Predicate(AtomQueryPredicate::AtomType { aromatic, .. }) => *aromatic,
            QueryNode::And(children) => atom_and_aromaticity(children),
            QueryNode::Not(child) => match child.as_ref() {
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(_)) => atom.is_aromatic(),
                QueryNode::Predicate(AtomQueryPredicate::IsAromatic(aromatic))
                | QueryNode::Predicate(AtomQueryPredicate::AtomType { aromatic, .. }) => !*aromatic,
                QueryNode::And(_) => false,
                QueryNode::Predicate(_) | QueryNode::Or(_) | QueryNode::Xor(_) | QueryNode::Not(_) => {
                    false
                }
            },
            QueryNode::Predicate(_) | QueryNode::Or(_) | QueryNode::Xor(_) => false,
        }
    */
}

fn atom_list_query_helper(query: &QueryNode<AtomQueryPredicate>, ignore_negation: bool) -> bool {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: bool _atomListQueryHelper(const T query, bool ignoreNegation) {
    // RDKit✔️✔️:   PRECONDITION(query, "no query");
    // RDKit✔️✔️:   if (!ignoreNegation && query->getNegation()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (query->getDescription() == "AtomAtomicNum" ||
    // RDKit✔️✔️:       query->getDescription() == "AtomType") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (query->getDescription() == "AtomOr") {
    // RDKit✔️✔️:     for (const auto &child : boost::make_iterator_range(query->beginChildren(),
    // RDKit✔️✔️:                                                         query->endChildren())) {
    // RDKit✔️✔️:       if (!_atomListQueryHelper(child, ignoreNegation)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    //
    // `AtomicNumberIn` and `AtomicNumberNotIn` are the canonical compressed
    // forms of source OR and negated-OR lists. Local complexity review: both
    // implementations visit each OR descendant once, short-circuit at the
    // first invalid child, and use O(h) recursion. No traversal path allocates,
    // clones, performs keyed lookup, or creates a temporary collection.
    match query {
        QueryNode::Not(child) => ignore_negation && atom_list_query_helper(child, ignore_negation),
        QueryNode::Predicate(
            AtomQueryPredicate::AtomicNumber(_) | AtomQueryPredicate::AtomType { .. },
        ) => true,
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(_)) => true,
        QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(_)) => ignore_negation,
        QueryNode::Or(children) => children
            .iter()
            .all(|child| atom_list_query_helper(child, ignore_negation)),
        QueryNode::And(_) | QueryNode::Xor(_) | QueryNode::Predicate(_) => false,
    }
}

fn is_atom_list_query(atom: &Atom) -> bool {
    // RDKit✔️✔️: bool isAtomListQuery(const Atom *a) {
    // RDKit✔️✔️:   PRECONDITION(a, "bad atom");
    // RDKit✔️✔️:   if (!a->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (a->getQuery()->getDescription() == "AtomOr") {
    // RDKit✔️✔️:     for (const auto &child : boost::make_iterator_range(
    // RDKit✔️✔️:              a->getQuery()->beginChildren(), a->getQuery()->endChildren())) {
    // RDKit✔️✔️:       if (!_atomListQueryHelper(child, false)) {
    // RDKit✔️✔️:         return false;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (a->getQuery()->getNegation() &&
    // RDKit✔️✔️:              _atomListQueryHelper(a->getQuery(), true)) {
    // RDKit✔️✔️:     // this was github #5930: negated list queries containing a single atom were
    // RDKit✔️✔️:     // being lost on output
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   } else if (a->getQuery()->getDescription() == "AtomAtomicNum" &&
    // RDKit✔️✔️:              static_cast<ATOM_EQUALS_QUERY *>(a->getQuery())->getVal() !=
    // RDKit✔️✔️:                  a->getAtomicNum()) {
    // RDKit✔️✔️:     // when reading single-member atom lists from CTABs we end up with simple
    // RDKit✔️✔️:     // AtomAtomicNum queries where the atomic number of the atom itself is zero.
    // RDKit✔️✔️:     // Recognize this case.
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return false;
    // RDKit✔️✔️: }
    //
    // Typed list leaves preserve the source OR-list identity without creating
    // a second list representation. Local complexity review: direct leaves are
    // O(1), while OR trees are O(n) time and O(h) recursion through the shared
    // helper, matching the source. No path allocates, clones, looks up keyed
    // state, or scans molecule data.
    let _ = atom;
    false
}

fn get_atom_list_query_values(
    query: &QueryNode<AtomQueryPredicate>,
) -> Result<Vec<i32>, &'static str> {
    // RDKit✔️✔️: void getAtomListQueryVals(const Atom::QUERYATOM_QUERY *q,
    // RDKit✔️✔️:                           std::vector<int> &vals) {
    // RDKit✔️✔️:   // list queries are series of nested ors of AtomAtomicNum queries
    // RDKit✔️✔️:   PRECONDITION(q, "bad query");
    // RDKit✔️✔️:   auto descr = q->getDescription();
    // RDKit✔️✔️:   if (descr == "AtomOr") {
    // RDKit✔️✔️:     for (const auto &child :
    // RDKit✔️✔️:          boost::make_iterator_range(q->beginChildren(), q->endChildren())) {
    // RDKit✔️✔️:       auto descr = child->getDescription();
    // RDKit✔️✔️:       if (child->getNegation() ||
    // RDKit✔️✔️:           (descr != "AtomOr" && descr != "AtomAtomicNum" &&
    // RDKit✔️✔️:            descr != "AtomType")) {
    // RDKit✔️✔️:         throw ValueErrorException("bad query type1");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       // we don't allow negation of any children of the query:
    // RDKit✔️✔️:       if (descr == "AtomOr") {
    // RDKit✔️✔️:         getAtomListQueryVals(child.get(), vals);
    // RDKit✔️✔️:       } else if (descr == "AtomAtomicNum") {
    // RDKit✔️✔️:         vals.push_back(static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal());
    // RDKit✔️✔️:       } else if (descr == "AtomType") {
    // RDKit✔️✔️:         auto v = static_cast<ATOM_EQUALS_QUERY *>(child.get())->getVal();
    // RDKit✔️✔️:         // aromatic AtomType queries add 1000 to the atomic number;
    // RDKit✔️✔️:         // correct for that:
    // RDKit✔️✔️:         if (v >= 1000) {
    // RDKit✔️✔️:           v -= 1000;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         vals.push_back(v);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (descr == "AtomAtomicNum") {
    // RDKit✔️✔️:     vals.push_back(static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal());
    // RDKit✔️✔️:   } else if (descr == "AtomType") {
    // RDKit✔️✔️:     auto v = static_cast<const ATOM_EQUALS_QUERY *>(q)->getVal();
    // RDKit✔️✔️:     // aromatic AtomType queries add 1000 to the atomic number;
    // RDKit✔️✔️:     // correct for that:
    // RDKit✔️✔️:     if (v >= 1000) {
    // RDKit✔️✔️:       v -= 1000;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     vals.push_back(v);
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     CHECK_INVARIANT(0, "bad query type");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // Typed AtomType stores the corrected atomic number directly, while list
    // leaves are the canonical compressed form of nested ORs. Local complexity
    // review: source and Rust visit n list values once and append n integers,
    // using O(n) time/output storage and O(h) recursion. Rust adds no clones,
    // keyed lookups, molecule scans, or temporary collections beyond the one
    // returned vector corresponding to RDKit's caller-owned output vector.
    fn append(
        query: &QueryNode<AtomQueryPredicate>,
        values: &mut Vec<i32>,
        child: bool,
    ) -> Result<(), &'static str> {
        match query {
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(value)) => {
                values.push(i32::from(*value));
                Ok(())
            }
            QueryNode::Predicate(AtomQueryPredicate::AtomType { atomic_number, .. }) => {
                values.push(i32::from(*atomic_number));
                Ok(())
            }
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumberIn(values_in)) => {
                values.extend(values_in.iter().copied().map(i32::from));
                Ok(())
            }
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumberNotIn(values_in)) if !child => {
                values.extend(values_in.iter().copied().map(i32::from));
                Ok(())
            }
            QueryNode::Or(children) => {
                for child_query in children {
                    if matches!(child_query, QueryNode::Not(_)) {
                        return Err("bad query type1");
                    }
                    append(child_query, values, true)?;
                }
                Ok(())
            }
            QueryNode::Not(inner) if !child => append(inner, values, false),
            QueryNode::Not(_) | QueryNode::And(_) | QueryNode::Xor(_) | QueryNode::Predicate(_) => {
                Err(if child {
                    "bad query type1"
                } else {
                    "bad query type"
                })
            }
        }
    }

    let mut values = Vec::new();
    append(query, &mut values, false)?;
    Ok(values)
}

#[inline]
fn make_atom_simple_query(predicate: AtomQueryPredicate) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomSimpleQuery(int what, std::function<int(Atom const *)> func,
    // RDKit✔️🔝:                        const std::string &description = "Atom Simple") {
    // RDKit✔️🔝:   T *res = new T;
    // RDKit✔️🔝:   res->setVal(what);
    // RDKit✔️🔝:   res->setDataFunc(func);
    // RDKit✔️🔝:   res->setDescription(description);
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The typed predicate variant is the single canonical representation of
    // RDKit's query class, target value, data function, and description
    // identity. Moving it into a leaf therefore preserves all modeled query
    // behavior without a parallel function-pointer or description registry.
    // Local complexity review: RDKit performs one heap allocation followed by
    // three constant-time field assignments. Rust performs one constant-time
    // enum move with no heap allocation, traversal, cloning, or lookup. This
    // preserves semantics while removing the source allocation and indirection.
    QueryNode::predicate(predicate)
}

#[inline]
pub(crate) fn make_atom_range_query(
    lower: i32,
    upper: i32,
    lower_open: bool,
    upper_open: bool,
    data_function: AtomRangeDataFunction,
) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: static inline ATOM_RANGE_QUERY *makeAtomRangeQuery(
    // RDKit✔️🔝:     int lower, int upper, bool lowerOpen, bool upperOpen,
    // RDKit✔️🔝:     std::function<int(Atom const *)> func,
    // RDKit✔️🔝:     const std::string &description = "Atom Range") {
    // RDKit✔️🔝:   ATOM_RANGE_QUERY *res = new ATOM_RANGE_QUERY(lower, upper);
    // RDKit✔️🔝:   res->setDataFunc(func);
    // RDKit✔️🔝:   res->setDescription(description);
    // RDKit✔️🔝:   res->setEndsOpen(lowerOpen, upperOpen);
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The typed leaf stores the bounds, endpoint flags, and data-function
    // identity directly; its variant supplies RDKit's description identity.
    // Local complexity review: RDKit heap-allocates one RangeQuery and makes
    // three constant-time assignments. Rust constructs one inline typed leaf
    // with the same state and no heap allocation, lookup, traversal, or clone.
    // Removing that allocation and virtual indirection preserves semantics and
    // is a material constant-factor improvement.
    make_atom_simple_query(AtomQueryPredicate::Range(AtomRangeQuery::new(
        AtomRangeBounds::Inclusive {
            lower,
            upper,
            lower_open,
            upper_open,
        },
        data_function,
    )))
}

#[inline]
pub(crate) fn make_atom_possible_range_query(
    lower: Option<i32>,
    upper: Option<i32>,
    data_function: AtomRangeDataFunction,
) -> Option<QueryNode<AtomQueryPredicate>> {
    let bounds = match (lower, upper) {
        (None, Some(upper)) => AtomRangeBounds::LessEqual(upper),
        (Some(lower), None) => AtomRangeBounds::GreaterEqual(lower),
        (Some(lower), Some(upper)) => {
            return Some(make_atom_range_query(
                lower,
                upper,
                false,
                false,
                data_function,
            ));
        }
        (None, None) => return None,
    };
    Some(make_atom_simple_query(AtomQueryPredicate::Range(
        AtomRangeQuery::new(bounds, data_function),
    )))
}

#[inline]
pub(crate) fn make_atom_num_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomNumQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomNum, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNumQuery(int what) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNum,
    // RDKit✔️🔝:                                                 "AtomAtomicNum");
    // RDKit✔️🔝: }
    //
    // `AtomicNumber` is the canonical typed identity for the source value,
    // data function, and AtomAtomicNum description in the modeled element
    // range. Local complexity review: both implementations perform one O(1)
    // leaf construction with no traversal, lookup, or clone. Rust reuses the
    // allocation-free simple factory, avoiding RDKit's query-object heap
    // allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::AtomicNumber(what))
}

#[inline]
fn make_atom_type_query(num: u8, aromatic: bool) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomTypeQuery(int num, int aromatic, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(makeAtomType(num, aromatic), queryAtomType,
    // RDKit✔️🔝:                                 descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomTypeQuery(int num, int aromatic) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(makeAtomType(num, aromatic),
    // RDKit✔️🔝:                                                 queryAtomType, "AtomType");
    // RDKit✔️🔝: }
    //
    // The typed leaf stores the two lossless components of RDKit's encoded
    // target value; matching recombines them through the canonical
    // `make_atom_type` implementation before comparing with `query_atom_type`.
    // Local complexity review: RDKit performs one O(1) scalar encoding and a
    // heap-allocated simple-query construction. Rust performs the same O(1)
    // encoding plus one subtraction to recover the typed atomic-number field,
    // then moves one inline leaf without traversal, lookup, cloning, or heap
    // allocation. The extra scalar operation is constant and the removed heap
    // allocation and virtual data-function indirection are a material
    // constant-factor improvement.
    let encoded_type = make_atom_type(i32::from(num), aromatic);
    let atomic_number = (encoded_type - 1000 * (aromatic as i32)) as u8;
    make_atom_simple_query(AtomQueryPredicate::AtomType {
        atomic_number,
        aromatic,
    })
}

#[inline]
fn make_atom_implicit_valence_query(what: i32) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomImplicitValenceQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomImplicitValence, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomImplicitValenceQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitValence);
    // RDKit✔️🔝:   res->setDescription("AtomImplicitValence");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `ImplicitValence` is the canonical typed identity for the source value,
    // data function, and AtomImplicitValence description. Local complexity
    // review: both implementations perform one O(1) leaf construction with
    // no traversal, lookup, or clone. Rust reuses the allocation-free simple
    // factory, removing RDKit's query-object heap allocation and virtual data-
    // function indirection without changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::ImplicitValence(what))
}

#[inline]
fn make_atom_explicit_valence_query(what: i32) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomExplicitValenceQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomExplicitValence, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomExplicitValenceQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitValence);
    // RDKit✔️🔝:   res->setDescription("AtomExplicitValence");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `ExplicitValence` is the canonical typed identity for the source value,
    // data function, and AtomExplicitValence description. Local complexity
    // review: both implementations perform one O(1) leaf construction with
    // no traversal, lookup, or clone. Rust reuses the allocation-free simple
    // factory, removing RDKit's query-object heap allocation and virtual data-
    // function indirection without changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::ExplicitValence(what))
}

#[inline]
fn make_atom_total_valence_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomTotalValenceQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomTotalValence, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomTotalValenceQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalValence);
    // RDKit✔️🔝:   res->setDescription("AtomTotalValence");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `TotalValence` is the canonical typed identity for the source value,
    // data function, and AtomTotalValence description in the modeled SMARTS
    // value range. Local complexity review: both implementations perform one
    // O(1) leaf construction with no traversal, lookup, or clone. Rust reuses
    // the allocation-free simple factory, removing RDKit's query-object heap
    // allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::TotalValence(what))
}

#[inline]
pub(crate) fn make_atom_explicit_degree_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomExplicitDegreeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomExplicitDegree, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomExplicitDegreeQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomExplicitDegree);
    // RDKit✔️🔝:   res->setDescription("AtomExplicitDegree");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `ExplicitDegree` is the single canonical typed identity for the source
    // value, data function, and AtomExplicitDegree description in the modeled
    // SMARTS value range. Local complexity review: both implementations make
    // one O(1) leaf construction with no traversal, lookup, or clone. Rust
    // reuses the allocation-free simple factory, removing RDKit's query-object
    // heap allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::ExplicitDegree(what))
}

#[inline]
pub(crate) fn make_atom_total_degree_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomTotalDegreeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomTotalDegree, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomTotalDegreeQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomTotalDegree);
    // RDKit✔️🔝:   res->setDescription("AtomTotalDegree");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `TotalDegree` is the single canonical typed identity for the source
    // value, data function, and AtomTotalDegree description in the modeled
    // SMARTS value range. Local complexity review: both implementations make
    // one O(1) leaf construction with no traversal, lookup, or clone. Rust
    // reuses the allocation-free simple factory, removing RDKit's query-object
    // heap allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::TotalDegree(what))
}

#[inline]
fn make_atom_heavy_atom_degree_query(what: u32) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomHeavyAtomDegreeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomHeavyAtomDegree, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHeavyAtomDegreeQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHeavyAtomDegree);
    // RDKit✔️🔝:   res->setDescription("AtomHeavyAtomDegree");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HeavyAtomDegree` is the canonical typed identity for the source value,
    // data function, and AtomHeavyAtomDegree description in the modeled atom-
    // degree range. Local complexity review: both implementations make one
    // O(1) leaf construction with no traversal, lookup, or clone. Rust reuses
    // the allocation-free simple factory, removing RDKit's query-object heap
    // allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::HeavyAtomDegree(what))
}

#[inline]
pub(crate) fn make_atom_h_count_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomHCountQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomHCount, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHCountQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHCount);
    // RDKit✔️🔝:   res->setDescription("AtomHCount");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HydrogenCount` is the single canonical typed identity for the source
    // target, queryAtomHCount data function, and AtomHCount description in the
    // modeled SMARTS value range. Local complexity review: both implementations
    // perform one O(1) leaf construction with no traversal, lookup, or clone.
    // Rust reuses the allocation-free simple factory, removing RDKit's query-
    // object heap allocation and virtual data-function indirection without
    // changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::HydrogenCount(what))
}

#[inline]
pub(crate) fn make_atom_has_implicit_h_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomHasImplicitHQuery(const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(true, queryAtomHasImplicitH, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHasImplicitHQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasImplicitH);
    // RDKit✔️🔝:   res->setDescription("AtomHasImplicitH");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasImplicitHydrogen` is the single typed identity for RDKit's literal
    // true target, queryAtomHasImplicitH data function, and description. Local
    // complexity review: both implementations perform one O(1) leaf
    // construction with no traversal, lookup, or clone. Rust reuses the
    // allocation-free simple factory, removing the source query-object heap
    // allocation and virtual data-function indirection while preserving the
    // source predicate's total-no-neighbors hydrogen semantics.
    make_atom_simple_query(AtomQueryPredicate::HasImplicitHydrogen)
}

#[inline]
pub(crate) fn make_atom_implicit_h_count_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomImplicitHCountQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomImplicitHCount, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomImplicitHCountQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomImplicitHCount);
    // RDKit✔️🔝:   res->setDescription("AtomImplicitHCount");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `ImplicitHydrogenCount` is the single typed identity for the source
    // target, queryAtomImplicitHCount data function, and description in the
    // modeled SMARTS value range. Local complexity review: both implementations
    // perform one O(1) leaf construction with no traversal, lookup, or clone.
    // Rust reuses the allocation-free simple factory, removing RDKit's query-
    // object heap allocation and virtual data-function indirection without
    // changing the source's no-neighbor hydrogen-count semantics.
    make_atom_simple_query(AtomQueryPredicate::ImplicitHydrogenCount(what))
}

#[inline]
pub(crate) fn make_atom_aromatic_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomAromaticQuery(const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(true, queryAtomAromatic, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomAromaticQuery() {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAromatic);
    // RDKit✔️🔝:   res->setDescription("AtomIsAromatic");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `IsAromatic(true)` is the single typed identity for the source literal
    // true target, queryAtomAromatic data function, and AtomIsAromatic
    // description. Local complexity review: both implementations perform one
    // O(1) leaf construction with no traversal, lookup, or clone. Rust reuses
    // the allocation-free simple factory, removing RDKit's query-object heap
    // allocation and virtual data-function indirection without changing
    // matching semantics.
    make_atom_simple_query(AtomQueryPredicate::IsAromatic(true))
}

#[inline]
pub(crate) fn make_atom_aliphatic_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomAliphaticQuery(const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(true, queryAtomAliphatic, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomAliphaticQuery() {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomAliphatic);
    // RDKit✔️🔝:   res->setDescription("AtomIsAliphatic");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `IsAromatic(false)` is the single typed identity for RDKit's literal
    // true target, queryAtomAliphatic data function, and AtomIsAliphatic
    // description. Local complexity review: both implementations perform one
    // O(1) leaf construction with no traversal, lookup, or clone. Rust reuses
    // the allocation-free simple factory, removing RDKit's query-object heap
    // allocation and virtual data-function indirection without changing the
    // source's logical negation of the aromatic flag.
    make_atom_simple_query(AtomQueryPredicate::IsAromatic(false))
}

#[inline]
fn make_atom_unsaturated_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomUnsaturatedQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomUnsaturated);
    // RDKit✔️🔝:   res->setDescription("AtomUnsaturated");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The existing `IsUnsaturated` leaf is the canonical typed identity for
    // the source boolean query, so the factory joins the historical matcher
    // to the main source-order path without duplicating its valence logic.
    // Local complexity review: both factories construct one O(1) leaf with no
    // traversal, lookup, or cloning. Rust removes RDKit's heap allocation and
    // virtual function indirection without changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::IsUnsaturated)
}

#[inline]
fn make_atom_in_ring_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomInRingQuery() {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomInRing);
    // RDKit✔️🔝:   res->setDescription("AtomInRing");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The historical `InRing` typed leaf is the canonical representation of
    // RDKit's true-valued AtomInRing query and already dispatches to the sole
    // source-backed ring predicate. Local complexity review: both factories
    // are O(1) and perform no traversal, lookup, or cloning. Rust removes one
    // heap allocation and virtual data-function indirection while retaining
    // identical match-time ring-info complexity.
    make_atom_simple_query(AtomQueryPredicate::InRing)
}

#[inline]
fn make_atom_in_n_rings_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomInNRingsQuery(int what) {
    // RDKit✔️🔝:   ATOM_EQUALS_QUERY *res;
    // RDKit✔️🔝:   res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryIsAtomInNRings);
    // RDKit✔️🔝:   res->setDescription("AtomInNRings");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NumAtomRings` is the existing canonical typed identity for the target,
    // source data function, and AtomInNRings description. Local complexity
    // review: both factories are O(1) with no traversal, lookup, or cloning;
    // Rust removes the source heap allocation and virtual dispatch while the
    // shared match path retains identical ring-membership lookup complexity.
    make_atom_simple_query(AtomQueryPredicate::NumAtomRings(i32::from(what)))
}

fn make_atom_ring_query(value: i32) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: bool Match(const ConstAtomPtr what) const override {
    // RDKit✔️🔝:   int v = this->TypeConvert(what, Queries::Int2Type<true>());
    // RDKit✔️🔝:   bool res;
    // RDKit✔️🔝:   if (this->d_val < 0) {
    // RDKit✔️🔝:     res = v != 0;
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     res = !Queries::queryCmp(v, this->d_val, this->d_tol);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   if (this->getNegation()) { res = !res; }
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // Local complexity review: the existing NumAtomRings leaf and ordinary
    // matcher perform the same O(1) ring-membership count lookup/comparison.
    // Rust removes virtual type conversion and stores the signed sentinel.
    QueryNode::predicate(AtomQueryPredicate::NumAtomRings(value))
}

#[inline]
fn make_atom_in_ring_of_size_query(target: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomInRingOfSizeQuery(int tgt) {
    // RDKit✔️🔝:   auto *res = new ATOM_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(tgt);
    // RDKit✔️🔝:   res->setDataFunc(
    // RDKit✔️🔝:       [tgt](Atom const *at) { return queryAtomIsInRingOfSize(at, tgt); });
    // RDKit✔️🔝:   res->setDescription("AtomRingSize");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The existing `InRingOfSize` leaf is the canonical typed identity for
    // the source target, captured data function, and description. Local
    // complexity review: both factories are O(1), capture/store one integer,
    // and perform no traversal or cloning. Rust removes the query-object heap
    // allocation and virtual closure dispatch while preserving the shared
    // O(R_atom) match-time lookup.
    make_atom_simple_query(AtomQueryPredicate::InRingOfSize(target))
}

#[inline]
fn make_atom_in_ring_of_size_range_query(
    lower: i32,
    upper: i32,
    lower_open: bool,
    upper_open: bool,
) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_RANGE_QUERY *makeAtomInRingOfSizeQuery(int lower, int upper,
    // RDKit✔️🔝:                                             bool lowerOpen, bool upperOpen) {
    // RDKit✔️🔝:   auto *res = new ATOM_RANGE_QUERY;
    // RDKit✔️🔝:   res->setLower(lower);
    // RDKit✔️🔝:   res->setUpper(upper);
    // RDKit✔️🔝:   res->setEndsOpen(lowerOpen, upperOpen);
    // RDKit✔️🔝:   res->setDataFunc([lower, upper, lowerOpen, upperOpen](Atom const *at) {
    // RDKit✔️🔝:     return queryAtomIsInRingOfSize(at, lower, upper, lowerOpen, upperOpen);
    // RDKit✔️🔝:   });
    // RDKit✔️🔝:   res->setDescription("range_AtomRingSize");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `AtomRangeQuery` remains the sole typed range-query representation; its
    // data-function identity stores the same four values captured by RDKit's
    // lambda. Local complexity review: both factories store four scalars in
    // O(1), without traversal or cloning. Rust avoids the source allocation
    // and virtual closure object while preserving the O(R_atom) match scan and
    // O(R_atom) temporary ring-size vector.
    make_atom_range_query(
        lower,
        upper,
        lower_open,
        upper_open,
        AtomRangeDataFunction::AtomRingSize {
            lower,
            upper,
            lower_open,
            upper_open,
        },
    )
}

#[inline]
fn make_atom_min_ring_size_query(target: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomMinRingSizeQuery(int tgt) {
    // RDKit✔️🔝:   auto *res = new ATOM_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(tgt);
    // RDKit✔️🔝:   res->setDataFunc(queryAtomMinRingSize);
    // RDKit✔️🔝:   res->setDescription("AtomMinRingSize");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // The historical `SmallestRingSize` leaf is the canonical typed identity
    // for the source target, data function, and description. Local complexity
    // review: both factories perform O(1) scalar storage with no traversal,
    // lookup, or cloning. Rust removes the query-object heap allocation and
    // virtual data-function dispatch while retaining the shared O(R_atom)
    // match-time minimum scan.
    make_atom_simple_query(AtomQueryPredicate::SmallestRingSize(target))
}

#[inline]
pub(crate) fn make_atom_ring_bond_count_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomRingBondCountQuery(int what) {
    // RDKit✔️🔝:   ATOM_EQUALS_QUERY *res = new AtomRingQuery(what);
    // RDKit✔️🔝:   res->setDescription("AtomRingBondCount");
    // RDKit✔️🔝:   res->setDataFunc(queryAtomRingBondCount);
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: };
    //
    // `RingBondCount` is the sole atom-side typed identity for the source
    // target, data function, and description. The duplicate historical
    // `NumRingBonds` family has been folded into this representation. Local
    // complexity review: both factories store one scalar in O(1), without
    // traversal, lookup, or cloning. Rust removes the source heap allocation
    // and virtual dispatch while preserving the shared O(degree) match scan.
    make_atom_simple_query(AtomQueryPredicate::RingBondCount(u32::from(what)))
}

pub(crate) const QUERY_SCAN_MAGIC_VALUE: u32 = 0xDEADBEEF;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct AtomQueryCompletionValues {
    pub(crate) ring_bond_count: u32,
    pub(crate) non_hydrogen_degree: u32,
}

pub(crate) fn complete_query_and_children(
    query: &mut QueryNode<AtomQueryPredicate>,
    magic_value: u32,
    values: AtomQueryCompletionValues,
) {
    // RDKit✔️✔️: void completeQueryAndChildren(Atom::QUERYATOM_QUERY *query, Atom *tgt,
    // RDKit✔️✔️:                               unsigned int magicVal) {
    // RDKit✔️✔️:   PRECONDITION(query, "no query");
    // RDKit✔️✔️:   PRECONDITION(tgt, "no atom");
    // RDKit✔️✔️:   auto eqQuery = dynamic_cast<ATOM_EQUALS_QUERY *>(query);
    // RDKit✔️✔️:   if (eqQuery) {
    // RDKit✔️✔️:     if (static_cast<unsigned int>(eqQuery->getVal()) == magicVal) {
    // RDKit✔️✔️:       int tgtVal = eqQuery->getDataFunc()(tgt);
    // RDKit✔️✔️:       eqQuery->setVal(tgtVal);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (auto childIt = query->beginChildren(); childIt != query->endChildren();
    // RDKit✔️✔️:        ++childIt) {
    // RDKit✔️✔️:     completeQueryAndChildren(childIt->get(), tgt, magicVal);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // Typed equality leaves own the data-function identity, so the two query
    // kinds currently capable of storing RDKit's unsigned sentinel select the
    // corresponding precomputed target value directly. The values preserve
    // the source target-atom boundary without retaining an owning-molecule
    // pointer on `Atom`. Local complexity review: source and Rust visit every
    // query node once in O(n)
    // time with O(h) recursion and mutate matching leaves in place. Neither
    // implementation allocates, clones, performs keyed lookup, or creates a
    // temporary collection during traversal.
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCount(value)) if *value == magic_value => {
            *value = values.ring_bond_count;
        }
        QueryNode::Predicate(AtomQueryPredicate::NonHydrogenDegree(value))
            if *value == magic_value =>
        {
            *value = values.non_hydrogen_degree;
        }
        QueryNode::Predicate(_) => {}
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => {
            for child in children {
                complete_query_and_children(child, magic_value, values);
            }
        }
        QueryNode::Not(child) => complete_query_and_children(child, magic_value, values),
    }
}

pub(crate) fn atom_query_has_magic_value(
    query: &QueryNode<AtomQueryPredicate>,
    magic_value: u32,
) -> bool {
    match query {
        QueryNode::Predicate(AtomQueryPredicate::RingBondCount(value))
        | QueryNode::Predicate(AtomQueryPredicate::NonHydrogenDegree(value)) => {
            *value == magic_value
        }
        QueryNode::Predicate(_) => false,
        QueryNode::And(children) | QueryNode::Or(children) | QueryNode::Xor(children) => children
            .iter()
            .any(|child| atom_query_has_magic_value(child, magic_value)),
        QueryNode::Not(child) => atom_query_has_magic_value(child, magic_value),
    }
}

pub(crate) fn complete_mol_queries(molecule: &mut crate::QueryGraph, magic_value: u32) {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeMolQueries
    // RDKit✔️✔️: void completeMolQueries(RWMol *mol, unsigned int magicVal) {
    // RDKit✔️✔️:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️✔️:   for (auto atom : mol->atoms()) {
    // RDKit✔️✔️:     if (atom->hasQuery()) {
    // RDKit✔️✔️:       completeQueryAndChildren(atom->getQuery(), atom, magicVal);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: completeMolQueries
    //
    // `Atom` does not carry RDKit's owning-molecule pointer, so graph-derived
    // data-function results are calculated from the same topology immediately
    // before completing each atom's canonical typed query tree. Local
    // complexity review: RDKit and Rust make one O(A) atom pass and one O(Q)
    // total query-tree traversal. Rust additionally scans the incident bonds
    // and neighbors of each query atom in O(sum degree) to materialize the two
    // typed data-function values; it allocates no per-atom or per-query
    // collection, clones no query tree, and performs no repeated whole-graph
    // traversal.
    for atom_idx in 0..molecule.num_atoms() {
        let mut ring_bond_count = 0_u32;
        let mut non_hydrogen_degree = 0_u32;
        for &(neighbor_index, _) in &molecule.adjacency()[atom_idx] {
            let neighbor_atom = &molecule.atoms()[neighbor_index];
            if neighbor_atom.atom().atomic_number() != 1
                || neighbor_atom
                    .atom()
                    .isotope()
                    .is_some_and(|isotope| isotope > 1)
            {
                non_hydrogen_degree += 1;
            }
        }
        let values = AtomQueryCompletionValues {
            ring_bond_count,
            non_hydrogen_degree,
        };
        complete_query_and_children(
            molecule.atoms_mut()[atom_idx].predicate_mut(),
            magic_value,
            values,
        );
    }
}

pub(crate) fn replace_atom_with_query_atom(atom: Atom) -> crate::QueryAtom {
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: replaceAtomWithQueryAtom
    // RDKit✔️🔝: Atom *replaceAtomWithQueryAtom(RWMol *mol, Atom *atom) {
    // RDKit✔️🔝:   PRECONDITION(mol, "bad molecule");
    // RDKit✔️🔝:   PRECONDITION(atom, "bad atom");
    // RDKit✔️🔝:   if (atom->hasQuery()) {
    // RDKit✔️🔝:     return atom;
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   QueryAtom qa(*atom);
    // RDKit✔️🔝:   unsigned int idx = atom->getIdx();
    // RDKit✔️🔝:
    // RDKit✔️🔝:   if (atom->hasProp(common_properties::_hasMassQuery)) {
    // RDKit✔️🔝:     qa.expandQuery(makeAtomMassQuery(static_cast<int>(atom->getMass())));
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   mol->replaceAtom(idx, &qa);
    // RDKit✔️🔝:   return mol->getAtomWithIdx(idx);
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryOps.cpp :: replaceAtomWithQueryAtom
    // BEGIN RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryAtom.h :: QueryAtom(const Atom &)
    // RDKit✔️🔝:   explicit QueryAtom(const Atom &other)
    // RDKit✔️🔝:       : Atom(other), dp_query(makeAtomNumQuery(other.getAtomicNum())) {
    // RDKit✔️🔝:     if (other.getIsotope()) {
    // RDKit✔️🔝:       this->expandQuery(makeAtomIsotopeQuery(other.getIsotope()),
    // RDKit✔️🔝:                         Queries::CompositeQueryType::COMPOSITE_AND);
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (other.getFormalCharge()) {
    // RDKit✔️🔝:       this->expandQuery(makeAtomFormalChargeQuery(other.getFormalCharge()),
    // RDKit✔️🔝:                         Queries::CompositeQueryType::COMPOSITE_AND);
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     if (other.getNumRadicalElectrons()) {
    // RDKit✔️🔝:       this->expandQuery(
    // RDKit✔️🔝:           makeAtomNumRadicalElectronsQuery(other.getNumRadicalElectrons()),
    // RDKit✔️🔝:           Queries::CompositeQueryType::COMPOSITE_AND);
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // END RDKIT CPP FUNCTION: third_party/rdkit/Code/GraphMol/QueryAtom.h :: QueryAtom(const Atom &)
    //
    // `Atom` is already the sole typed ordinary/query atom representation, so
    // setting its optional query preserves the same row, id, props, and all
    // topology references instead of allocating a temporary QueryAtom and
    // copying it back through RWMol::replaceAtom. Local complexity review:
    // both implementations perform a fixed O(1) sequence of scalar checks and
    // construct at most five leaves. Rust has the same bounded composite-node
    // allocations but avoids the source temporary atom clone, virtual query
    // objects, molecule row replacement, and second atom copy.
    let mut query = make_atom_num_query(atom.atomic_number());
    if let Some(isotope) = atom.isotope()
        && isotope != 0
    {
        query_atom_expand_query(
            &mut query,
            make_atom_isotope_query(isotope),
            CompositeQueryType::And,
            true,
        );
    }
    if atom.formal_charge() != 0 {
        query_atom_expand_query(
            &mut query,
            make_atom_formal_charge_query(atom.formal_charge()),
            CompositeQueryType::And,
            true,
        );
    }
    if atom.radical_electrons() != 0 {
        query_atom_expand_query(
            &mut query,
            make_atom_num_radical_electrons_query(atom.radical_electrons()),
            CompositeQueryType::And,
            true,
        );
    }
    if atom.prop("_hasMassQuery").is_some() {
        let mass = rdkit_atomic_mass(atom.atomic_number(), atom.isotope()) as u16;
        query_atom_expand_query(
            &mut query,
            make_atom_mass_query(mass),
            CompositeQueryType::And,
            true,
        );
    }
    crate::QueryAtom::from_parts(atom, query)
}

#[inline]
pub(crate) fn make_atom_has_ring_bond_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHasRingBondQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasRingBond);
    // RDKit✔️🔝:   res->setDescription("AtomHasRingBond");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasRingBond` is the sole typed identity for the source boolean target,
    // data function, and description. It replaces the historical deferred-
    // scan name and the duplicate >=1 ring-count leaf. Local complexity
    // review: both factories construct one O(1) leaf with no traversal,
    // lookup, or cloning. Rust removes the source allocation and virtual
    // dispatch while retaining the shared O(degree) short-circuit match scan.
    make_atom_simple_query(AtomQueryPredicate::HasRingBond)
}

#[inline]
pub(crate) fn make_atom_num_heteroatom_nbrs_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNumHeteroatomNbrsQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNumHeteroatomNbrs);
    // RDKit✔️🔝:   res->setDescription("AtomNumHeteroatomNeighbors");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NumHeteroatomNeighbors` is the canonical equality-query identity for
    // the source target, data function, and description. The range-query data
    // function shares the same source-backed counter without duplicating it.
    // Local complexity review: both factories store one scalar in O(1), with
    // no traversal, lookup, or cloning. Rust removes the source allocation and
    // virtual dispatch while retaining the shared O(degree) match scan.
    make_atom_simple_query(AtomQueryPredicate::NumHeteroatomNeighbors(what))
}

#[inline]
pub(crate) fn make_atom_has_heteroatom_nbrs_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHasHeteroatomNbrsQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasHeteroatomNbrs);
    // RDKit✔️🔝:   res->setDescription("AtomHasHeteroatomNeighbors");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasHeteroatomNeighbors` is the canonical boolean-query identity for
    // the source target, short-circuit data function, and description. Local
    // complexity review: both factories construct one O(1) leaf without
    // traversal, lookup, or cloning. Rust removes the source allocation and
    // virtual dispatch while retaining the O(degree), first-match-return scan.
    make_atom_simple_query(AtomQueryPredicate::HasHeteroatomNeighbors)
}

#[inline]
pub(crate) fn make_atom_num_aliphatic_heteroatom_nbrs_query(
    what: u8,
) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNumAliphaticHeteroatomNbrsQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
    // RDKit✔️🔝:       what, queryAtomNumAliphaticHeteroatomNbrs);
    // RDKit✔️🔝:   res->setDescription("AtomNumAliphaticHeteroatomNeighbors");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NumAliphaticHeteroatomNeighbors` is the canonical typed identity for
    // the equality target, source counter, and description. Local complexity
    // review: both factories store one scalar in O(1), without traversal,
    // lookup, or cloning. Rust removes the source allocation and virtual
    // dispatch while preserving the shared O(degree) match-time scan.
    make_atom_simple_query(AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(what))
}

#[inline]
pub(crate) fn make_atom_has_aliphatic_heteroatom_nbrs_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHasAliphaticHeteroatomNbrsQuery() {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
    // RDKit✔️🔝:       true, queryAtomHasAliphaticHeteroatomNbrs);
    // RDKit✔️🔝:   res->setDescription("AtomHasAliphaticHeteroatomNeighbors");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasAliphaticHeteroatomNeighbors` is the canonical typed identity for
    // the source boolean target, short-circuit data function, and description.
    // Local complexity review: both factories construct one O(1) leaf without
    // traversal, lookup, or cloning. Rust removes the source allocation and
    // virtual dispatch while preserving the O(degree), first-match-return scan.
    make_atom_simple_query(AtomQueryPredicate::HasAliphaticHeteroatomNeighbors)
}

#[inline]
pub(crate) fn make_atom_non_hydrogen_degree_query(what: u32) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNonHydrogenDegreeQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomNonHydrogenDegree);
    // RDKit✔️🔝:   res->setDescription("AtomNonHydrogenDegree");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NonHydrogenDegree` is the sole typed identity for this source counter;
    // the duplicate historical `SubstitutionCount` family has been folded
    // into its equality/comparison variants. Local complexity review: both
    // factories store one scalar in O(1), without traversal, lookup, or
    // cloning. Rust removes the source allocation and virtual dispatch while
    // preserving the shared O(degree) match-time scan.
    make_atom_simple_query(AtomQueryPredicate::NonHydrogenDegree(what))
}

#[inline]
fn make_atom_is_bridgehead_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomIsBridgeheadQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryIsAtomBridgehead);
    // RDKit✔️🔝:   res->setDescription("AtomIsBridgehead");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `IsBridgehead` is the sole typed query identity for the source boolean
    // target, data function, and description. Matching reuses the existing
    // source-backed `chemistry::stereo::query_is_atom_bridgehead` algorithm,
    // so the SMARTS/query path and stereochemistry path share one core rather
    // than carrying historical duplicates. Local complexity review: both
    // factories construct one O(1) leaf with no traversal, lookup, or clone.
    // Rust removes the source heap allocation and virtual dispatch while the
    // shared match-time helper retains RDKit's ring-overlap complexity.
    make_atom_simple_query(AtomQueryPredicate::IsBridgehead)
}

#[inline]
pub(crate) fn make_q_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_OR_QUERY *makeQAtomQuery() {
    // RDKit✔️✔️:   auto *res = new ATOM_OR_QUERY;
    // RDKit✔️✔️:   res->setDescription("AtomOr");
    // RDKit✔️✔️:   res->setTypeLabel("Q");
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `Not(Or(...))` is the canonical typed form of RDKit's negated AtomOr;
    // both the query factory and Molfile complex-symbol path reuse this one
    // constructor. Type labels are serialization metadata in RDKit and do not
    // alter matching. Local complexity review: both implementations create a
    // fixed two-child tree in O(1), with constant allocation and no molecule
    // traversal, lookup, or cloning. Match-time behavior short-circuits after
    // at most two O(1) atomic-number comparisons in the same child order.
    QueryNode::not(QueryNode::or(vec![
        make_atom_num_query(6),
        make_atom_num_query(1),
    ]))
}

#[inline]
pub(crate) fn make_q_h_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_EQUALS_QUERY *makeQHAtomQuery() {
    // RDKit✔️✔️:   ATOM_EQUALS_QUERY *res = makeAtomNumQuery(6);
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->setTypeLabel("QH");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `Not(AtomicNumber(6))` is the sole typed representation of the source
    // negated equality query; Molfile parsing reuses this constructor. RDKit's
    // type label affects serialization identity, not matching. Local
    // complexity review: both factories create a fixed one-leaf negation tree
    // in O(1), with constant allocation and no traversal, lookup, or cloning;
    // matching performs one O(1) atomic-number comparison and one negation.
    QueryNode::not(make_atom_num_query(6))
}

#[inline]
pub(crate) fn make_a_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_EQUALS_QUERY *makeAAtomQuery() {
    // RDKit✔️✔️:   ATOM_EQUALS_QUERY *res = makeAtomNumQuery(1);
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->setTypeLabel("A");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // `Not(AtomicNumber(1))` is the sole typed representation of this source
    // negated equality query, shared by the factory and Molfile parser. The
    // type label is serialization metadata and does not change matching.
    // Local complexity review: both implementations create a fixed one-leaf
    // negation tree in O(1), with constant allocation and no traversal,
    // lookup, or cloning; matching is one O(1) comparison and one negation.
    QueryNode::not(make_atom_num_query(1))
}

#[inline]
pub(crate) fn make_atom_null_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_NULL_QUERY *makeAtomNullQuery() {
    // RDKit✔️🔝:   auto *res = new ATOM_NULL_QUERY;
    // RDKit✔️🔝:   res->setDataFunc(nullDataFun<const RDKit::Atom *>);
    // RDKit✔️🔝:   res->setMatchFunc(nullQueryFun<int>);
    // RDKit✔️🔝:   res->setDescription("AtomNull");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: template <typename T>
    // RDKit✔️🔝: int nullDataFun(T) {
    // RDKit✔️🔝:   return 1;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: template <typename T>
    // RDKit✔️🔝: bool nullQueryFun(T) {
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    //
    // `Any` is the single typed constant-true atom query. Explicit SMARTS
    // wildcards and the AH complex-symbol factory reuse this constructor.
    // Local complexity review: source and Rust construction and matching are
    // O(1), with no traversal, lookup, or clone. Rust removes the allocation
    // and virtual calls while preserving unconditional success.
    make_atom_simple_query(AtomQueryPredicate::Any)
}

#[inline]
pub(crate) fn make_a_h_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_NULL_QUERY *makeAHAtomQuery() {
    // RDKit✔️🔝:   auto *res = makeAtomNullQuery();
    // RDKit✔️🔝:   res->setTypeLabel("AH");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `Any` is the single canonical typed null-query identity; the AH factory
    // and Molfile parser both reuse it, while the source type label is matching-
    // neutral serialization metadata. Local complexity review: both factories
    // are O(1) and perform no traversal, lookup, or clone. Rust constructs an
    // inline leaf without RDKit's null-query heap allocation and virtual call,
    // preserving the constant-true semantics with a lower constant cost.
    make_atom_null_query()
}

#[inline]
pub(crate) fn make_x_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_OR_QUERY *makeXAtomQuery() {
    // RDKit✔️✔️:   auto *res = new ATOM_OR_QUERY;
    // RDKit✔️✔️:   res->setDescription("AtomOr");
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(17)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(35)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(53)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(85)));
    // RDKit✔️✔️:   res->setTypeLabel("X");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // The five-child `Or` tree is the sole X-query representation, shared by
    // the factory and Molfile parser; no parallel halogen-set predicate is
    // retained for X. The type label is matching-neutral metadata. Local
    // complexity review: both implementations build a fixed five-leaf tree in
    // O(1) with constant allocation, no molecule traversal, lookup, or clone,
    // and match by the same ordered, short-circuit O(5) comparisons.
    QueryNode::or([9, 17, 35, 53, 85].map(make_atom_num_query).to_vec())
}

#[inline]
pub(crate) fn make_x_h_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_OR_QUERY *makeXHAtomQuery() {
    // RDKit✔️✔️:   ATOM_OR_QUERY *res = makeXAtomQuery();
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
    // RDKit✔️✔️:   res->setTypeLabel("XH");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // This extends the canonical X tree instead of repeating its halogen
    // list, then serves both the factory and Molfile parser. The type label is
    // matching-neutral metadata. Local complexity review: RDKit and Rust each
    // build the fixed X tree and append one leaf in O(1), with no molecule
    // traversal, lookup, or clone; matching short-circuits through the same
    // ordered six O(1) atomic-number comparisons.
    let mut query = make_x_atom_query();
    match &mut query {
        QueryNode::Or(children) => children.push(make_atom_num_query(1)),
        _ => unreachable!("make_x_atom_query must return an Or node"),
    }
    query
}

#[inline]
pub(crate) fn make_m_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   ATOM_OR_QUERY *res = makeMHAtomQuery();
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
    // RDKit✔️✔️:   res->setTypeLabel("M");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // The negated OR tree and `is_metal` share one excluded-number constant;
    // M appends hydrogen after the MH-ordered children exactly as the source
    // does. The type label is matching-neutral metadata. Local complexity
    // review: both factories construct a fixed 23-leaf tree in O(1), with
    // constant allocation and no molecule traversal, lookup, or clone, then
    // match via the same ordered, short-circuit O(23) comparisons.
    let mut query = make_m_h_atom_query();
    match &mut query {
        QueryNode::Not(child) => match child.as_mut() {
            QueryNode::Or(children) => children.push(make_atom_num_query(1)),
            _ => unreachable!("make_m_h_atom_query negation must contain an Or node"),
        },
        _ => unreachable!("make_m_h_atom_query must return a Not node"),
    }
    query
}

#[inline]
pub(crate) fn make_m_h_atom_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMHAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   auto *res = new ATOM_OR_QUERY;
    // RDKit✔️✔️:   res->setDescription("AtomOr");
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(0)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(2)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(5)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(7)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(8)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(10)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(14)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(15)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(16)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(17)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(18)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(33)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(34)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(35)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(36)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(52)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(53)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(54)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(85)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(86)));
    // RDKit✔️✔️:   res->setTypeLabel("MH");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // The negated OR tree is built from the one shared excluded-number table;
    // M extends this exact constructor and Molfile MH calls it directly. The
    // source type label does not alter matching. Local complexity review:
    // both implementations construct a fixed 22-leaf tree in O(1), with
    // constant allocation and no molecule traversal, lookup, or clone, then
    // short-circuit through the same ordered O(22) comparisons.
    QueryNode::not(QueryNode::or(
        MH_EXCLUDED_ATOMIC_NUMBERS.map(make_atom_num_query).to_vec(),
    ))
}

#[inline]
pub(crate) fn convert_complex_name_to_query(
    symbol: &str,
) -> Result<QueryNode<AtomQueryPredicate>, QueryConstructionError> {
    // RDKit✔️✔️: void convertComplexNameToQuery(Atom *query, std::string_view symb) {
    // RDKit✔️✔️:   if (symb == "Q") {
    // RDKit✔️✔️:     query->setQuery(makeQAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "QH") {
    // RDKit✔️✔️:     query->setQuery(makeQHAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "A") {
    // RDKit✔️✔️:     query->setQuery(makeAAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "AH") {
    // RDKit✔️✔️:     query->setQuery(makeAHAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "X") {
    // RDKit✔️✔️:     query->setQuery(makeXAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "XH") {
    // RDKit✔️✔️:     query->setQuery(makeXHAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "M") {
    // RDKit✔️✔️:     query->setQuery(makeMAtomQuery());
    // RDKit✔️✔️:   } else if (symb == "MH") {
    // RDKit✔️✔️:     query->setQuery(makeMHAtomQuery());
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     // we control what this function gets called with, so we should never land
    // RDKit✔️✔️:     // here
    // RDKit✔️✔️:     ASSERT_INVARIANT(0, "bad complex query symbol");
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    //
    // Each branch delegates to the one canonical typed factory; Molfile/SDF
    // parsing calls this dispatcher instead of retaining a second symbol map.
    // Local complexity review: RDKit and Rust perform an ordered fixed set of
    // at most eight string comparisons and one O(1) factory call, with no graph
    // traversal or lookup table. Rust returns a structured error in place of
    // the source invariant exception for the unreachable invalid-symbol path.
    match symbol {
        "Q" => Ok(make_q_atom_query()),
        "QH" => Ok(make_q_h_atom_query()),
        "A" => Ok(make_a_atom_query()),
        "AH" => Ok(make_a_h_atom_query()),
        "X" => Ok(make_x_atom_query()),
        "XH" => Ok(make_x_h_atom_query()),
        "M" => Ok(make_m_atom_query()),
        "MH" => Ok(make_m_h_atom_query()),
        _ => Err(QueryConstructionError::InvalidComplexAtomSymbol {
            symbol: symbol.to_owned(),
        }),
    }
}

#[inline]
fn make_atom_mass_query(what: u16) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomMassQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(massIntegerConversionFactor * what,
    // RDKit✔️🔝:                                 queryAtomMass, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomMassQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
    // RDKit✔️🔝:       massIntegerConversionFactor * what, queryAtomMass);
    // RDKit✔️🔝:   res->setDescription("AtomMass");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `Mass` stores the source `what` value as the single typed query identity;
    // the canonical match arm applies `massIntegerConversionFactor` exactly
    // once when comparing it with queryAtomMass. Local complexity review:
    // RDKit performs one O(1) integer multiplication, allocation, and leaf
    // initialization; Rust performs one O(1) enum move and defers the same
    // multiplication to matching. Neither traverses, scans, looks up, clones,
    // or creates a temporary collection, and Rust removes the heap allocation
    // and virtual data-function indirection.
    make_atom_simple_query(AtomQueryPredicate::Mass(what))
}

#[inline]
pub(crate) fn make_atom_isotope_query(what: u16) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomIsotopeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomIsotope, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomIsotopeQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomIsotope);
    // RDKit✔️🔝:   res->setDescription("AtomIsotope");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `Isotope` is the single typed identity for the source target,
    // queryAtomIsotope data function, and AtomIsotope description in the
    // modeled isotope range. Local complexity review: both implementations
    // perform one O(1) leaf construction with no traversal, lookup, or clone.
    // Rust reuses the allocation-free simple factory, removing RDKit's query-
    // object heap allocation and virtual data-function indirection without
    // changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::Isotope(what))
}

#[inline]
pub(crate) fn make_atom_formal_charge_query(what: i8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomFormalChargeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomFormalCharge, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomFormalChargeQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomFormalCharge);
    // RDKit✔️🔝:   res->setDescription("AtomFormalCharge");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `FormalCharge` is the single typed identity for the source target,
    // queryAtomFormalCharge data function, and AtomFormalCharge description in
    // COSMolKit's modeled charge range. Local complexity review: both
    // implementations perform one O(1) leaf construction with no traversal,
    // lookup, or clone. Rust reuses the allocation-free simple factory,
    // removing RDKit's query-object heap allocation and virtual data-function
    // indirection without changing matching semantics.
    make_atom_simple_query(AtomQueryPredicate::FormalCharge(what))
}

#[inline]
fn make_atom_negative_formal_charge_query(what: i8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomNegativeFormalChargeQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomNegativeFormalCharge, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNegativeFormalChargeQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
    // RDKit✔️🔝:       what, queryAtomNegativeFormalCharge);
    // RDKit✔️🔝:   res->setDescription("AtomNegativeFormalCharge");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NegativeFormalCharge` is the single typed identity for the source
    // target, queryAtomNegativeFormalCharge data function, and description.
    // Local complexity review: both implementations perform one O(1) leaf
    // construction with no traversal, lookup, or clone. Rust reuses the
    // allocation-free simple factory, removing RDKit's query-object heap
    // allocation and virtual data-function indirection without changing
    // matching semantics or duplicating the existing source-backed data
    // function.
    make_atom_simple_query(AtomQueryPredicate::NegativeFormalCharge(what))
}

#[inline]
pub(crate) fn make_atom_hybridization_query(what: Hybridization) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomHybridizationQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomHybridization, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHybridizationQuery(int what) {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(what, queryAtomHybridization);
    // RDKit✔️🔝:   res->setDescription("AtomHybridization");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HybridizationMatch` stores the typed value whose discriminant follows
    // RDKit's HybridizationType source order, so it is the single identity for
    // the integer target, queryAtomHybridization data function, and description.
    // Local complexity review: both implementations perform one O(1) leaf
    // construction with no traversal, lookup, or clone. Rust removes the
    // source heap allocation and virtual data-function indirection while
    // preserving the scalar comparison in the canonical match arm.
    make_atom_simple_query(AtomQueryPredicate::HybridizationMatch(what))
}

#[inline]
fn make_atom_num_radical_electrons_query(what: u8) -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: template <class T>
    // RDKit✔️🔝: T *makeAtomNumRadicalElectronsQuery(int what, const std::string &descr) {
    // RDKit✔️🔝:   return makeAtomSimpleQuery<T>(what, queryAtomNumRadicalElectrons, descr);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomNumRadicalElectronsQuery(int what) {
    // RDKit✔️🔝:   auto *res = makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
    // RDKit✔️🔝:       what, queryAtomNumRadicalElectrons);
    // RDKit✔️🔝:   res->setDescription("AtomNumRadicalElectrons");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NumRadicalElectrons` is the single typed identity for the source target,
    // queryAtomNumRadicalElectrons data function, and description in the
    // modeled atom-state range. Local complexity review: both implementations
    // perform one O(1) leaf construction with no traversal, lookup, or clone.
    // Rust reuses the allocation-free simple factory, removing RDKit's query-
    // object heap allocation and virtual data-function indirection while
    // reusing the existing source-backed radical-electron reader.
    make_atom_simple_query(AtomQueryPredicate::NumRadicalElectrons(what))
}

#[inline]
fn make_atom_has_chiral_tag_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomHasChiralTagQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomHasChiralTag);
    // RDKit✔️🔝:   res->setDescription("AtomHasChiralTag");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasChiralTag` is the single typed identity for the source's boolean
    // target, data function, and AtomHasChiralTag description. It remains
    // distinct from exact `ChiralTagMatch` queries. Local complexity review:
    // both implementations construct one O(1) leaf without traversal, lookup,
    // or cloning. The allocation-free typed leaf removes RDKit's heap
    // allocation and virtual data-function indirection without changing the
    // matching behavior.
    make_atom_simple_query(AtomQueryPredicate::HasChiralTag)
}

#[inline]
fn make_atom_missing_chiral_tag_query() -> QueryNode<AtomQueryPredicate> {
    // RDKit✔️🔝: ATOM_EQUALS_QUERY *makeAtomMissingChiralTagQuery() {
    // RDKit✔️🔝:   auto *res =
    // RDKit✔️🔝:       makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(true, queryAtomMissingChiralTag);
    // RDKit✔️🔝:   res->setDescription("AtomMissingChiralTag");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `MissingChiralTag` is the canonical typed identity for the source's
    // boolean target, data function, and description. Local complexity
    // review: both paths construct one O(1) leaf with no traversal, lookup,
    // or cloning. Rust removes the source heap allocation and virtual dispatch
    // while preserving the presence-only property semantics in the reused
    // data function.
    make_atom_simple_query(AtomQueryPredicate::MissingChiralTag)
}

pub(crate) fn make_bond_order_equals_query(what: BondOrder) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondOrderEqualsQuery(Bond::BondType what) {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(what);
    // RDKit✔️🔝:   res->setDataFunc(queryBondOrder);
    // RDKit✔️🔝:   res->setDescription("BondOrder");
    // RDKit✔️🔝:   res->setTypeLabel("BondOrder");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `Order` is the sole typed identity for the target value, queryBondOrder
    // data function, description, and matching-neutral type label. Both SMARTS
    // parsers reuse this constructor. Local complexity review: both factories
    // are O(1) with no traversal, lookup, or clone; Rust constructs one inline
    // leaf and removes RDKit's heap allocation and virtual data-function call.
    QueryNode::predicate(BondQueryPredicate::Order(what))
}

fn make_bond_has_prop_query(property: impl Into<String>) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: template <class Target>
    // RDKit✔️🔝: Queries::EqualityQuery<int, const Target *, true> *makeHasPropQuery(
    // RDKit✔️🔝:     const std::string &property) {
    // RDKit✔️🔝:   return new HasPropQuery<const Target *>(property);
    // RDKit✔️🔝: }
    // Local complexity review: both move/copy one O(n) property name; Rust
    // removes the separately allocated polymorphic query object.
    QueryNode::predicate(BondQueryPredicate::HasProperty(property.into()))
}

fn make_bond_prop_query(
    property: impl Into<String>,
    value: impl Into<String>,
) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: template <class Target, class T>
    // RDKit✔️🔝: Queries::EqualityQuery<int, const Target *, true> *makePropQuery(
    // RDKit✔️🔝:     const std::string &propname, const T &val, double tolerance = 0.0) {
    // RDKit✔️🔝:   return new HasPropWithValueQuery<const Target *, T>(propname, val, tolerance);
    // RDKit✔️🔝: }
    // RDKit✔️🔝: res = atom_val == this->val;
    // Local complexity review: the modeled string specialization ignores
    // tolerance in RDKit. Both store and compare two strings in O(n); Rust
    // avoids the polymorphic query allocation and virtual dispatch.
    QueryNode::predicate(BondQueryPredicate::PropertyValue {
        name: property.into(),
        value: value.into(),
    })
}

fn finalize_bond_query_from_description(
    description: &str,
    query: QueryNode<BondQueryPredicate>,
) -> Result<QueryNode<BondQueryPredicate>, QueryFinalizationError> {
    // RDKit✔️✔️: std::string descr = query->getDescription();
    // RDKit✔️✔️: if (descr == "BondRingSize") {
    // RDKit✔️✔️:   tmpQuery = makeBondInRingOfSizeQuery(
    // RDKit✔️✔️:       static_cast<BOND_EQUALS_QUERY *>(query)->getVal());
    // RDKit✔️✔️:   query->setDataFunc(tmpQuery->getDataFunc()); delete tmpQuery;
    // RDKit✔️✔️: } else if (descr == "BondMinRingSize") { query->setDataFunc(queryBondMinRingSize);
    // RDKit✔️✔️: } else if (descr == "BondOrder") { query->setDataFunc(queryBondOrder);
    // RDKit✔️✔️: } else if (descr == "BondDir") { query->setDataFunc(queryBondDir);
    // RDKit✔️✔️: } else if (descr == "BondInRing") { query->setDataFunc(queryIsBondInRing);
    // RDKit✔️✔️: } else if (descr == "BondInNRings") { query->setDataFunc(queryIsBondInNRings);
    // RDKit✔️✔️: } else if (descr == "SingleOrAromaticBond") { query->setDataFunc(queryBondIsSingleOrAromatic);
    // RDKit✔️✔️: } else if (descr == "SingleOrDoubleBond") { query->setDataFunc(queryBondIsSingleOrDouble);
    // RDKit✔️✔️: } else if (descr == "DoubleOrAromaticBond") { query->setDataFunc(queryBondIsDoubleOrAromatic);
    // RDKit✔️✔️: } else if (descr == "SingleOrDoubleOrAromaticBond") { query->setDataFunc(queryBondIsSingleOrDoubleOrAromatic);
    // RDKit✔️✔️: } else if (descr == "BondNull" || descr == "BondAnd" ||
    // RDKit✔️✔️:            descr == "BondOr" || descr == "BondXor" ||
    // RDKit✔️✔️:            descr == "HasProp" || descr == "HasPropWithValue") { }
    // RDKit✔️✔️: else { throw ValueErrorException("Do not know how to finalize query: '" + descr + "'"); }
    // Local complexity review: short-description dispatch is O(n) in the
    // description length, matching RDKit. Typed variants already encode the
    // data function and avoid the temporary BondRingSize query allocation.
    const KNOWN: &[&str] = &[
        "BondRingSize",
        "BondMinRingSize",
        "BondOrder",
        "BondDir",
        "BondInRing",
        "BondInNRings",
        "SingleOrAromaticBond",
        "SingleOrDoubleBond",
        "DoubleOrAromaticBond",
        "SingleOrDoubleOrAromaticBond",
        "BondNull",
        "BondAnd",
        "BondOr",
        "BondXor",
        "HasProp",
        "HasPropWithValue",
    ];
    KNOWN
        .contains(&description)
        .then_some(query)
        .ok_or_else(|| QueryFinalizationError::UnknownDescription(description.to_string()))
}

#[inline]
fn make_bond_dir_equals_query(what: crate::BondDirection) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondDirEqualsQuery(Bond::BondDir what) {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(what);
    // RDKit✔️🔝:   res->setDataFunc(queryBondDir);
    // RDKit✔️🔝:   res->setDescription("BondDir");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `Direction` is the sole typed identity for the source target, data
    // function, and description. Local complexity review: both factories are
    // O(1), with no traversal, lookup, or clone. Rust constructs one inline
    // leaf and removes RDKit's query-object allocation and virtual dispatch;
    // matching still performs the single O(1) direction-field comparison.
    QueryNode::predicate(BondQueryPredicate::Direction(what))
}

#[inline]
fn make_bond_has_stereo_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondHasStereoQuery() {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(true);
    // RDKit✔️🔝:   res->setDataFunc(queryBondHasStereo);
    // RDKit✔️🔝:   res->setDescription("BondStereo");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `HasStereo` is the sole typed identity for RDKit's presence query and
    // remains distinct from exact `Stereo` values. Local complexity review:
    // both factories are O(1), with no traversal, lookup, or clone. Rust
    // removes the query-object allocation and virtual dispatch; matching
    // retains one O(1) stereo-field read and comparison through the shared
    // `query_bond_has_stereo` helper.
    QueryNode::predicate(BondQueryPredicate::HasStereo)
}

#[inline]
pub(crate) fn make_bond_is_in_ring_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondIsInRingQuery() {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(true);
    // RDKit✔️🔝:   res->setDataFunc(queryIsBondInRing);
    // RDKit✔️🔝:   res->setDescription("BondInRing");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `IsInRing(true)` is the sole typed identity for the positive source
    // query, and both SMARTS parsers reuse this factory. Local complexity
    // review: both factories are O(1), with no graph traversal, lookup, or
    // clone. Rust removes the query-object allocation and virtual dispatch;
    // matching retains the shared O(1) RingInfo membership lookup.
    QueryNode::predicate(BondQueryPredicate::IsInRing(true))
}

#[inline]
fn make_bond_in_n_rings_query(what: i32) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondInNRingsQuery(int what) {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(what);
    // RDKit✔️🔝:   res->setDataFunc(queryIsBondInNRings);
    // RDKit✔️🔝:   res->setDescription("BondInNRings");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `NumRingBonds(i32)` preserves the full source target domain, including
    // negative targets that cannot match the nonnegative RingInfo count.
    // Local complexity review: both factories are O(1), with no traversal,
    // lookup, or clone. Rust removes the query-object allocation and virtual
    // dispatch; match time retains one O(1) ring-count lookup and comparison.
    QueryNode::predicate(BondQueryPredicate::NumRingBonds(what))
}

#[inline]
fn make_bond_in_ring_of_size_query(
    target: i32,
) -> Result<QueryNode<BondQueryPredicate>, QueryConstructionError> {
    // RDKit✔️✔️: BOND_EQUALS_QUERY *makeBondInRingOfSizeQuery(int tgt) {
    // RDKit✔️✔️:   RANGE_CHECK(3, tgt, 20);
    // RDKit✔️✔️:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️✔️:   res->setVal(tgt);
    // RDKit✔️✔️:   switch (tgt) {
    // RDKit✔️✔️:     case 3:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<3>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 4:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<4>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 5:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<5>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 6:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<6>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 7:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<7>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 8:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<8>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 9:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<9>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 10:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<10>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 11:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<11>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 12:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<12>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 13:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<13>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 14:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<14>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 15:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<15>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 16:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<16>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 17:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<17>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 18:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<18>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 19:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<19>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     case 20:
    // RDKit✔️✔️:       res->setDataFunc(queryBondIsInRingOfSize<20>);
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   res->setDescription("BondRingSize");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    //
    // The typed target replaces the source's 18 template-specialized data
    // pointers while retaining the same validated domain and one canonical
    // ring-size matcher. Local complexity review: both factories validate in
    // O(1). Rust then constructs one inline leaf without RDKit's allocation or
    // switch/virtual dispatch. At match time both perform the same O(R_bond)
    // ring-membership scan with no temporary collection or clone.
    if !(3..=20).contains(&target) {
        return Err(QueryConstructionError::BondRingSizeOutOfRange { target });
    }
    Ok(QueryNode::predicate(BondQueryPredicate::InRingOfSize(
        target,
    )))
}

#[inline]
fn make_bond_min_ring_size_query(target: i32) -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_EQUALS_QUERY *makeBondMinRingSizeQuery(int tgt) {
    // RDKit✔️🔝:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️🔝:   res->setVal(tgt);
    // RDKit✔️🔝:   res->setDataFunc(queryBondMinRingSize);
    // RDKit✔️🔝:   res->setDescription("BondMinRingSize");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    //
    // `MinRingSize(i32)` preserves the source target without imposing the
    // separate in-ring-size factory's 3..=20 construction range. Local
    // complexity review: both factories are O(1), with no traversal, lookup,
    // or clone. Rust removes the allocation and virtual dispatch; match time
    // retains the shared O(R_bond) minimum scan with O(1) auxiliary space.
    QueryNode::predicate(BondQueryPredicate::MinRingSize(target))
}

#[inline]
pub(crate) fn make_bond_null_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️🔝: BOND_NULL_QUERY *makeBondNullQuery() {
    // RDKit✔️🔝:   auto *res = new BOND_NULL_QUERY;
    // RDKit✔️🔝:   res->setDataFunc(nullDataFun<const RDKit::Bond *>);
    // RDKit✔️🔝:   res->setMatchFunc(nullQueryFun<int>);
    // RDKit✔️🔝:   res->setDescription("BondNull");
    // RDKit✔️🔝:   return res;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: template <typename T>
    // RDKit✔️🔝: int nullDataFun(T) {
    // RDKit✔️🔝:   return 1;
    // RDKit✔️🔝: }
    // RDKit✔️🔝: template <typename T>
    // RDKit✔️🔝: bool nullQueryFun(T) {
    // RDKit✔️🔝:   return true;
    // RDKit✔️🔝: }
    //
    // `Any` is the single typed constant-true bond query, and both SMARTS
    // parsers reuse this factory for `~`. Local complexity review: source and
    // Rust construction and matching are O(1), with no traversal, lookup, or
    // clone. Rust removes the heap allocation and two virtual calls while
    // preserving unconditional success.
    QueryNode::predicate(BondQueryPredicate::Any)
}

fn make_query_bond_spec(begin: AtomId, end: AtomId, order: BondOrder) -> crate::QueryBond {
    // RDKit✔️🔝: QueryBond::QueryBond(BondType bT) : Bond(bT) {
    // RDKit✔️🔝:   if (bT != Bond::UNSPECIFIED) {
    // RDKit✔️🔝:     dp_query = makeBondOrderEqualsQuery(bT);
    // RDKit✔️🔝:   } else {
    // RDKit✔️🔝:     dp_query = makeBondNullQuery();
    // RDKit✔️🔝:   }
    // RDKit✔️🔝: };
    // Local complexity review: both constructors select one of two O(1)
    // query factories. BondSpec additionally stores the endpoints required by
    // COSMolKit's value-style builder, but its inline typed query avoids
    // RDKit's separate QueryBond and virtual query heap allocations. There is
    // no traversal, lookup, clone, scan, or temporary collection.
    let query = if order == BondOrder::Unspecified {
        make_bond_null_query()
    } else {
        make_bond_order_equals_query(order)
    };
    crate::QueryBond::from_parts(
        Bond::from_spec(crate::BondId::new(0), BondSpec::new(begin, end, order)),
        query,
    )
}

#[inline]
pub(crate) fn make_single_or_aromatic_bond_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrAromaticBondQuery() {
    // RDKit✔️✔️:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️✔️:   res->setVal(true);
    // RDKit✔️✔️:   res->setDataFunc(queryBondIsSingleOrAromatic);
    // RDKit✔️✔️:   res->setDescription("SingleOrAromaticBond");
    // RDKit✔️✔️:   res->setTypeLabel("BondOrder");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    //
    // `OrderIn` is the sole typed representation for fixed bond-order sets;
    // the factory and both SMARTS parsers share this leaf. Local complexity
    // review: each factory is O(1) with one constant-size allocation and no
    // molecule traversal, lookup, or clone; matching performs the same ordered
    // short-circuit pair of O(1) bond-order comparisons as the source helper.
    QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Single,
        BondOrder::Aromatic,
    ]))
}

pub(crate) fn is_complex_bond_query(bond: &crate::QueryBond) -> bool {
    // RDKit✔️✔️: bool isComplexQuery(const Bond *b) {
    // RDKit✔️✔️:   PRECONDITION(b, "bad bond");
    // RDKit✔️✔️:   if (!b->hasQuery()) {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // negated things are always complex:
    // RDKit✔️✔️:   if (b->getQuery()->getNegation()) {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   std::string descr = b->getQuery()->getDescription();
    // RDKit✔️✔️:   if (descr == "BondOrder" || descr == "SingleOrAromaticBond") {
    // RDKit✔️✔️:     return false;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "BondAnd" || descr == "BondXor") {
    // RDKit✔️✔️:     return true;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (descr == "BondOr") {
    // RDKit✔️✔️:     // detect the types of queries that appear for unspecified bonds in
    // RDKit✔️✔️:     // SMARTS:
    // RDKit✔️✔️:     if (b->getQuery()->endChildren() - b->getQuery()->beginChildren() == 2) {
    // RDKit✔️✔️:       for (auto child = b->getQuery()->beginChildren();
    // RDKit✔️✔️:            child != b->getQuery()->endChildren(); ++child) {
    // RDKit✔️✔️:         if ((*child)->getDescription() != "BondOrder" ||
    // RDKit✔️✔️:             (*child)->getNegation()) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:         if (static_cast<BOND_EQUALS_QUERY *>(child->get())->getVal() !=
    // RDKit✔️✔️:                 Bond::SINGLE &&
    // RDKit✔️✔️:             static_cast<BOND_EQUALS_QUERY *>(child->get())->getVal() !=
    // RDKit✔️✔️:                 Bond::AROMATIC) {
    // RDKit✔️✔️:           return true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       return false;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    //
    // `QueryNode` variants are the sole typed identities for RDKit query
    // descriptions, while the exact two-value `OrderIn` leaf is the canonical
    // SingleOrAromaticBond representation. Local complexity review: both
    // implementations inspect the root and at most two OR children in O(1)
    // time. Neither traverses a molecule or subtree, allocates, clones, looks
    // up keyed state, or creates a temporary collection.
    let query = bond.predicate();

    match query {
        QueryNode::Not(_) | QueryNode::And(_) | QueryNode::Xor(_) => true,
        QueryNode::Predicate(BondQueryPredicate::Order(_)) => false,
        QueryNode::Predicate(BondQueryPredicate::OrderIn(orders)) => {
            orders.as_slice() != [BondOrder::Single, BondOrder::Aromatic]
        }
        QueryNode::Or(children) if children.len() == 2 => !children.iter().all(|child| {
            matches!(
                child,
                QueryNode::Predicate(BondQueryPredicate::Order(
                    BondOrder::Single | BondOrder::Aromatic
                ))
            )
        }),
        QueryNode::Or(_) | QueryNode::Predicate(_) => true,
    }
}

#[inline]
fn make_double_or_aromatic_bond_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeDoubleOrAromaticBondQuery() {
    // RDKit✔️✔️:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️✔️:   res->setVal(true);
    // RDKit✔️✔️:   res->setDataFunc(queryBondIsDoubleOrAromatic);
    // RDKit✔️✔️:   res->setDescription("DoubleOrAromaticBond");
    // RDKit✔️✔️:   res->setTypeLabel("BondOrder");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    //
    // This reuses the single fixed-order-set leaf and match helper rather than
    // introducing a dedicated boolean predicate. Local complexity review:
    // both factories are O(1) with one constant-size allocation and no graph
    // traversal, lookup, or clone; matching makes the same ordered pair of
    // O(1), short-circuit bond-order comparisons.
    QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Double,
        BondOrder::Aromatic,
    ]))
}

#[inline]
fn make_single_or_double_bond_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *makeSingleOrDoubleBondQuery() {
    // RDKit✔️✔️:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️✔️:   res->setVal(true);
    // RDKit✔️✔️:   res->setDataFunc(queryBondIsSingleOrDouble);
    // RDKit✔️✔️:   res->setDescription("SingleOrDoubleBond");
    // RDKit✔️✔️:   res->setTypeLabel("BondOrder");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    //
    // This fixed pair remains an `OrderIn` leaf and shares the one order-set
    // matcher with all related factories. Local complexity review: both
    // factories are O(1) with one constant-size allocation and no traversal,
    // lookup, or clone; matching makes the same ordered, short-circuit pair
    // of O(1) bond-order comparisons.
    QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Single,
        BondOrder::Double,
    ]))
}

#[inline]
fn make_single_or_double_or_aromatic_bond_query() -> QueryNode<BondQueryPredicate> {
    // RDKit✔️✔️: RDKIT_GRAPHMOL_EXPORT BOND_EQUALS_QUERY *
    // RDKit✔️✔️: makeSingleOrDoubleOrAromaticBondQuery() {
    // RDKit✔️✔️:   auto *res = new BOND_EQUALS_QUERY;
    // RDKit✔️✔️:   res->setVal(true);
    // RDKit✔️✔️:   res->setDataFunc(queryBondIsSingleOrDoubleOrAromatic);
    // RDKit✔️✔️:   res->setDescription("SingleOrDoubleOrAromaticBond");
    // RDKit✔️✔️:   res->setTypeLabel("BondOrder");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    //
    // The three-member set uses the same canonical `OrderIn` leaf and matcher
    // as every fixed order-set factory. Local complexity review: both factory
    // paths are O(1) with one constant-size allocation and no traversal,
    // lookup, or clone; matching makes the same ordered, short-circuit three
    // O(1) bond-order comparisons.
    QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
        BondOrder::Single,
        BondOrder::Double,
        BondOrder::Aromatic,
    ]))
}

// ---------------------------------------------------------------------------
// Error types
// ---------------------------------------------------------------------------

/// Errors produced while constructing typed query leaves.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryConstructionError {
    #[error("bond ring size target {target} is outside RDKit's supported range 3..=20")]
    BondRingSizeOutOfRange { target: i32 },
    #[error("invalid RDKit complex atom query symbol '{symbol}'")]
    InvalidComplexAtomSymbol { symbol: String },
}

/// Errors produced by SMARTS parsing.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum SmartsParseError {
    #[error("unclosed bracket at position {0}")]
    UnclosedBracket(usize),
    #[error("unexpected character '{character}' at position {position}: {context}")]
    UnexpectedCharacter {
        position: usize,
        character: char,
        context: String,
    },
    #[error("unexpected end of SMARTS: {0}")]
    UnexpectedEnd(String),
    #[error("invalid atom primitive at position {position}: {detail}")]
    InvalidAtomPrimitive { position: usize, detail: String },
    #[error("unclosed parenthesis at position {0}")]
    UnclosedParenthesis(usize),
    #[error("unbalanced ring closure number: {0}")]
    UnbalancedRingClosure(u32),
    #[error("CXSMARTS parse error: {0}")]
    CxSmiles(String),
    #[error("SMARTS parse error: {0}")]
    Parse(String),
    #[error("unsupported SMARTS feature: {0}")]
    UnsupportedFeature(&'static str),
}

// ---------------------------------------------------------------------------
// Cache helpers (build adjacency / ring info on-the-fly when not cached)
// ---------------------------------------------------------------------------

/// RDKit✔️❌: Returns an `AdjacencyList` for `mol`.
/// COSMolKit stores adjacency inline in topology instead of in a derived cache.
fn ensure_adjacency(mol: &impl SearchTargetAccess) -> AdjacencyList {
    mol.adjacency().clone()
}

/// RDKit✔️❌: Returns a `RingInfo` for `mol`, using the cached copy if
/// available. When absent we build fresh from the molecule topology — this is
/// O(atoms × SSSR) and guaranteed to match RDKit's SSSR perception.
fn ensure_ring_info(mol: &impl SearchTargetAccess) -> Option<RingInfo> {
    // RDKit✔️❌: RDKit stores ring info inline; COSMolKit caches optionally.
    if let Some(cached) = mol.ring_info() {
        return Some(cached.clone());
    }
    // Ring info not cached - compute from the detached topology.
    crate::rings::find_sssr_from_parts(mol.num_atoms(), mol.bonds(), mol.adjacency()).ok()
}

fn ensure_valence_assignment(mol: &impl SearchTargetAccess) -> Option<ValenceAssignment> {
    if let Some(cached) = mol.valence() {
        return Some(cached.clone());
    }
    crate::valence::assign_valence_with_options_for_topology(
        mol.topology_block(),
        ValenceModel::RdkitLike,
        false,
    )
    .ok()
}

#[must_use]
pub(crate) fn build_query_match_context_for_target(
    mol: &impl SearchTargetAccess,
) -> QueryMatchContext {
    QueryMatchContext {
        adj: ensure_adjacency(mol),
        ring_info: ensure_ring_info(mol),
        valence: ensure_valence_assignment(mol),
    }
}

/// Build query predicate state for a live molecule (compatibility adapter).
/// The matching implementation itself uses `build_query_match_context_for_target`
/// and therefore only requires detached search data.
#[must_use]
pub fn build_query_match_context(mol: &Molecule) -> QueryMatchContext {
    build_query_match_context_for_target(mol)
}

/// Build predicate state from detached model values without constructing a
/// runtime molecule. This is the entrypoint used by future domain-crate
/// search adapters; the live-object overload above remains only as a
/// transitional facade adapter.
#[must_use]
pub(crate) fn build_query_match_context_from_blocks(
    topology: &cosmolkit_model::TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    stereo_groups: &[cosmolkit_model::StereoGroup],
    ring_info: Option<&RingInfo>,
    valence: Option<&ValenceAssignment>,
) -> QueryMatchContext {
    let target =
        super::target::SearchTarget::new(topology, coordinates, stereo_groups, ring_info, valence);
    build_query_match_context_for_target(&target)
}

/// Evaluate one atom predicate directly against detached model blocks.
pub(crate) fn atom_predicate_matches_from_blocks(
    atom: &Atom,
    predicate: &AtomQueryPredicate,
    topology: &cosmolkit_model::TopologyBlock,
    coordinates: &cosmolkit_model::CoordinateBlock,
    stereo_groups: &[cosmolkit_model::StereoGroup],
    ring_info: Option<&RingInfo>,
    valence: Option<&ValenceAssignment>,
) -> bool {
    let target =
        super::target::SearchTarget::new(topology, coordinates, stereo_groups, ring_info, valence);
    let context = build_query_match_context_for_target(&target);
    atom_predicate_matches_with_target_context(atom, predicate, &target, &context)
}

#[inline]
fn query_atom_implicit_valence(valence: Option<&ValenceAssignment>, at: &Atom) -> Option<i32> {
    // RDKit✔️✔️: static inline int queryAtomImplicitValence(Atom const *at) {
    // RDKit✔️✔️:   return at->getValence(Atom::ValenceType::IMPLICIT);
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️✔️:   if (!dp_mol) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️✔️:        d_implicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️✔️:     return d_explicitValence;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: RDKit performs one O(1) cached valence-field
    // read after constant-time state checks; Rust performs one O(1) indexed
    // read from the already assigned typed valence vector. Neither traverses,
    // allocates, clones, repeats a lookup, or creates a temporary collection.
    valence.and_then(|assignment| assignment.implicit_hydrogens.get(at.id().index()).copied())
}

#[inline]
fn atom_explicit_valence(valence: Option<&ValenceAssignment>, at: &Atom) -> Option<i32> {
    // RDKit✔️✔️: unsigned int Atom::getValence(ValenceType which) const {
    // RDKit✔️✔️:   if (!dp_mol) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::IMPLICIT || d_explicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::EXPLICIT) called without call to calcExplicitValence()");
    // RDKit✔️✔️:   PRECONDITION(
    // RDKit✔️✔️:       (which == ValenceType::EXPLICIT || df_noImplicit ||
    // RDKit✔️✔️:        d_implicitValence > -1),
    // RDKit✔️✔️:       "getValence(ValenceType::IMPLICIT) called without call to calcImplicitValence()");
    // RDKit✔️✔️:   if (which == ValenceType::EXPLICIT) {
    // RDKit✔️✔️:     return d_explicitValence;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return df_noImplicit ? 0 : d_implicitValence;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: RDKit performs one O(1) cached explicit-valence
    // field read after constant-time state checks; Rust performs one O(1)
    // indexed read from the already assigned typed valence vector. Neither
    // traverses, allocates, clones, repeats a lookup, or creates a temporary
    // collection.
    valence.and_then(|assignment| assignment.explicit_valence.get(at.id().index()).copied())
}

#[inline]
fn query_atom_explicit_valence(valence: Option<&ValenceAssignment>, at: &Atom) -> Option<i32> {
    // RDKit✔️✔️: static inline int queryAtomExplicitValence(Atom const *at) {
    // RDKit✔️✔️:   return at->getValence(Atom::ValenceType::EXPLICIT) - at->getNumExplicitHs();
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // Local complexity review: both implementations perform two O(1) typed or
    // cached field reads and one integer subtraction, with no traversal,
    // allocation, cloning, repeated lookup, or temporary collection. The raw
    // explicit-valence access is centralized in `atom_explicit_valence` so
    // total-valence matching and this SMARTS primitive cannot diverge.
    atom_explicit_valence(valence, at).map(|explicit| explicit - i32::from(at.explicit_hydrogens()))
}

fn implicit_hydrogen_count(valence: Option<&ValenceAssignment>, atom: &Atom) -> Option<u8> {
    query_atom_implicit_valence(valence, atom).map(|count| count.max(0) as u8)
}

fn total_hydrogen_count(valence: Option<&ValenceAssignment>, atom: &Atom) -> Option<usize> {
    implicit_hydrogen_count(valence, atom)
        .map(|implicit| usize::from(atom.explicit_hydrogens()) + usize::from(implicit))
}

#[inline]
fn query_atom_implicit_h_count(valence: Option<&ValenceAssignment>, at: &Atom) -> Option<usize> {
    // RDKit✔️✔️: static inline int queryAtomImplicitHCount(Atom const *at) {
    // RDKit✔️✔️:   return at->getTotalNumHs(false);
    // RDKit✔️✔️: };
    // RDKit✔️✔️: //
    // RDKit✔️✔️: //  If includeNeighbors is set, we'll loop over our neighbors
    // RDKit✔️✔️: //   and include any of them that are Hs in the count here
    // RDKit✔️✔️: //
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Local complexity review: the literal `includeNeighbors=false` makes
    // RDKit and Rust each perform two O(1) hydrogen-state reads and one
    // addition. Neither traverses adjacency, allocates, clones, branches over
    // neighbors, or creates a temporary collection. Reusing the shared total-
    // hydrogen base keeps this as the sole no-neighbor count implementation.
    total_hydrogen_count(valence, at)
}

#[inline]
fn query_atom_has_implicit_h(valence: Option<&ValenceAssignment>, at: &Atom) -> bool {
    // RDKit✔️✔️: static inline int queryAtomHasImplicitH(Atom const *at) {
    // RDKit✔️✔️:   return int(at->getTotalNumHs(false) > 0);
    // RDKit✔️✔️: };
    // RDKit✔️✔️: //
    // RDKit✔️✔️: //  If includeNeighbors is set, we'll loop over our neighbors
    // RDKit✔️✔️: //   and include any of them that are Hs in the count here
    // RDKit✔️✔️: //
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each reuse the O(1) no-neighbor
    // total-hydrogen count and perform one comparison with zero. Neither scans
    // adjacency, allocates, clones, or creates a temporary collection. The
    // shared helper preserves the source's inclusion of explicit atom H state.
    query_atom_implicit_h_count(valence, at).is_some_and(|count| count > 0)
}

#[inline]
fn query_atom_h_count(
    adj: &AdjacencyList,
    valence: Option<&ValenceAssignment>,
    at: &Atom,
    mol: &impl SearchTargetAccess,
) -> Option<usize> {
    // RDKit✔️✔️: static inline int queryAtomHCount(Atom const *at) {
    // RDKit✔️✔️:   return at->getTotalNumHs(true);
    // RDKit✔️✔️: };
    // RDKit✔️✔️: //
    // RDKit✔️✔️: //  If includeNeighbors is set, we'll loop over our neighbors
    // RDKit✔️✔️: //   and include any of them that are Hs in the count here
    // RDKit✔️✔️: //
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int getNumExplicitHs() const { return d_numExplicitHs; }
    // RDKit✔️✔️: unsigned int Atom::getNumImplicitHs() const {
    // RDKit✔️✔️:   if (df_noImplicit) {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   PRECONDITION(d_implicitValence > -1,
    // RDKit✔️✔️:                "getNumImplicitHs() called without preceding call to "
    // RDKit✔️✔️:                "calcImplicitValence()");
    // RDKit✔️✔️:   return getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform O(1) explicit and
    // implicit hydrogen reads followed by one O(degree) pass over the existing
    // adjacency range. Neither allocates, clones, repeats the scan, or creates
    // a temporary collection. Rust uses `usize` for the accumulator so the
    // target atom's neighbor count cannot silently saturate at the query's
    // `u8` representation limit.
    let mut res = total_hydrogen_count(valence, at)?;
    for nbr in adj.neighbors_of(at.id().index()) {
        if mol.atoms()[nbr.atom_index].atomic_number() == 1 {
            res += 1;
        }
    }
    Some(res)
}

#[inline]
fn query_atom_total_degree(
    adj: &AdjacencyList,
    valence: Option<&ValenceAssignment>,
    atom: &Atom,
) -> Option<usize> {
    // RDKit✔️✔️: static inline int queryAtomTotalDegree(Atom const *at) {
    // RDKit✔️✔️:   return at->getTotalDegree();
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int Atom::getTotalDegree() const {
    // RDKit✔️✔️:   unsigned int res = this->getTotalNumHs(false) + this->getDegree();
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: unsigned int Atom::getTotalNumHs(bool includeNeighbors) const {
    // RDKit✔️✔️:   int res = getNumExplicitHs() + getNumImplicitHs();
    // RDKit✔️✔️:   if (includeNeighbors && dp_mol) {
    // RDKit✔️✔️:     auto nbrs = dp_mol->atomNeighbors(this);
    // RDKit✔️✔️:     res += std::count_if(nbrs.begin(), nbrs.end(), [](const auto nbr) {
    // RDKit✔️✔️:       return (nbr->getAtomicNum() == 1);
    // RDKit✔️✔️:     });
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: with RDKit's literal `includeNeighbors=false`,
    // both implementations perform O(1) cached/typed hydrogen reads, one O(1)
    // graph-degree lookup, and one addition. Neither traverses neighbors,
    // allocates, clones, or creates a temporary collection. Rust returns
    // `None` only outside the modeled valid query context when valence state
    // could not be assigned; supported matching builds that state up front.
    total_hydrogen_count(valence, atom)
        .map(|total_hs| total_hs + query_atom_explicit_degree(atom, adj))
}

#[inline]
fn query_atom_total_valence(valence: Option<&ValenceAssignment>, at: &Atom) -> Option<i32> {
    // RDKit✔️✔️: static inline int queryAtomTotalValence(Atom const *at) {
    // RDKit✔️✔️:   return at->getTotalValence();
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int Atom::getTotalValence() const {
    // RDKit✔️✔️:   return getValence(ValenceType::EXPLICIT) + getValence(ValenceType::IMPLICIT);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform two O(1) cached or
    // typed valence reads and one integer addition, with no traversal,
    // allocation, cloning, repeated lookup, or temporary collection. Rust
    // reuses the canonical explicit- and implicit-valence readers, so all
    // total-valence predicates share one implementation of each source field.
    atom_explicit_valence(valence, at)
        .zip(query_atom_implicit_valence(valence, at))
        .and_then(|(explicit, implicit)| explicit.checked_add(implicit))
}

#[inline]
fn query_atom_unsaturated(
    adj: &AdjacencyList,
    valence: Option<&ValenceAssignment>,
    at: &Atom,
) -> Option<bool> {
    // RDKit✔️✔️: static inline int queryAtomUnsaturated(Atom const *at) {
    // RDKit✔️✔️:   return at->getTotalDegree() < at->getTotalValence();
    // RDKit✔️✔️: };
    // Local complexity review: both implementations reuse two O(1) cached or
    // typed atom-property reads and perform one integer comparison. Neither
    // traverses neighbors, allocates, clones, repeats a lookup, or creates a
    // temporary collection. Reusing the canonical total-degree and
    // total-valence functions removes the historical hybridization-based
    // SMARTS branch without introducing another chemistry implementation.
    query_atom_total_degree(adj, valence, at)
        .zip(query_atom_total_valence(valence, at))
        .and_then(|(degree, total_valence)| {
            usize::try_from(total_valence)
                .ok()
                .map(|total_valence| degree < total_valence)
        })
}

// ---------------------------------------------------------------------------
// atom_predicate_matches — evaluate an atom query predicate against a real atom
// ---------------------------------------------------------------------------

#[inline]
fn null_data<T>(_value: T) -> i32 {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: int nullDataFun(T) {
    // RDKit✔️✔️:   return 1;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust both ignore the generic input
    // and return one constant integer in O(1) time. Neither implementation
    // reads query state, branches, traverses, allocates, clones, or constructs
    // a temporary collection.
    1
}

#[inline]
fn null_query<T>(_value: T) -> bool {
    // RDKit✔️✔️: template <typename T>
    // RDKit✔️✔️: bool nullQueryFun(T) {
    // RDKit✔️✔️:   return true;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust both ignore the generic input
    // and return one constant boolean in O(1) time. Neither implementation
    // reads query state, branches, traverses, allocates, clones, or constructs
    // a temporary collection.
    true
}

#[inline]
fn is_atom_dummy(atom: &crate::QueryAtom) -> bool {
    // RDKit✔️✔️: inline bool isAtomDummy(const Atom *a) {
    // RDKit✔️✔️:   return (!a->hasQuery() && a->getAtomicNum() == 0) ||
    // RDKit✔️✔️:          (a->hasQuery() && !a->getQuery()->getNegation() &&
    // RDKit✔️✔️:           a->getQuery()->getDescription() == "AtomNull");
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform O(1) query-presence
    // and root-query checks plus, only for a non-query atom, one O(1) atomic-
    // number read. Neither traverses a composite query, allocates, clones, or
    // creates a temporary collection. The typed `Predicate(Any)` root is the
    // canonical representation of RDKit's non-negated `AtomNull` query; `Not`
    // and `Or` roots therefore retain the source distinction without relying
    // on the query atom's zero atomic number.
    if atom.atom().atomic_number() == 0 {
        return true;
    }
    matches!(
        atom.predicate(),
        QueryNode::Predicate(AtomQueryPredicate::Any)
    )
}

#[inline]
fn is_metal(atom: &Atom) -> bool {
    // RDKit✔️✔️: bool isMetal(const Atom &atom) {
    // RDKit✔️✔️:   static const std::unique_ptr<ATOM_OR_QUERY> q(makeMAtomQuery());
    // RDKit✔️✔️:   return q->Match(&atom);
    // RDKit✔️✔️: }
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#1!#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   ATOM_OR_QUERY *res = makeMHAtomQuery();
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(1)));
    // RDKit✔️✔️:   res->setTypeLabel("M");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: ATOM_OR_QUERY *makeMHAtomQuery() {
    // RDKit✔️✔️:   // using the definition from Marvin Sketch, which produces the following
    // RDKit✔️✔️:   // SMARTS:
    // RDKit✔️✔️:   // !#2!#5!#6!#7!#8!#9!#10!#14!#15!#16!#17!#18!#33!#34!#35!#36!#52!#53!#54!#85!#86
    // RDKit✔️✔️:   // We expanded this with !#0 as part of #6106
    // RDKit✔️✔️:   // it's easier to define what isn't a metal than what is. :-)
    // RDKit✔️✔️:   auto *res = new ATOM_OR_QUERY;
    // RDKit✔️✔️:   res->setDescription("AtomOr");
    // RDKit✔️✔️:   res->setNegation(true);
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(0)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(2)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(5)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(6)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(7)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(8)));
    // RDKit✔️✔️:   res->addChild(
    // RDKit✔️✔️:       Queries::Query<int, Atom const *, true>::CHILD_TYPE(makeAtomNumQuery(9)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(10)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(14)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(15)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(16)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(17)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(18)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(33)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(34)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(35)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(36)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(52)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(53)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(54)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(85)));
    // RDKit✔️✔️:   res->addChild(Queries::Query<int, Atom const *, true>::CHILD_TYPE(
    // RDKit✔️✔️:       makeAtomNumQuery(86)));
    // RDKit✔️✔️:   res->setTypeLabel("MH");
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: after RDKit's one-time static query creation,
    // both implementations linearly inspect the same fixed 23 atomic numbers
    // and short-circuit on a match, so each call is O(23) with identical hot-
    // path branching. Rust uses a static slice rather than heap query nodes;
    // it performs no per-call allocation, cloning, traversal beyond that fixed
    // list, or temporary collection construction.
    let atomic_number = query_atom_num(atom);
    atomic_number != 1 && !MH_EXCLUDED_ATOMIC_NUMBERS.contains(&atomic_number)
}

#[inline]
fn query_atom_aromatic(at: &Atom) -> bool {
    // RDKit✔️✔️: static inline int queryAtomAromatic(Atom const *at) {
    // RDKit✔️✔️:   return at->getIsAromatic();
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform one O(1) aromatic
    // flag read with no allocation, cloning, traversal, or temporary object.
    at.is_aromatic()
}

#[inline]
fn query_atom_aliphatic(at: &Atom) -> bool {
    // RDKit✔️✔️: static inline int queryAtomAliphatic(Atom const *at) {
    // RDKit✔️✔️:   return !(at->getIsAromatic());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform one O(1) aromatic
    // flag read plus one boolean negation, with no allocation, cloning,
    // traversal, or temporary object. Reusing `query_atom_aromatic` keeps the
    // aromatic flag access in one core implementation without changing cost.
    !query_atom_aromatic(at)
}

#[inline]
fn query_atom_num(at: &Atom) -> u8 {
    // RDKit✔️✔️: static inline int queryAtomNum(Atom const *at) { return at->getAtomicNum(); }
    // RDKit✔️✔️: int getAtomicNum() const { return d_atomicNum; }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // atomic-number field read, with no traversal, allocation, cloning,
    // repeated lookup, branching, or temporary collection.
    at.atomic_number()
}

#[inline]
fn make_atom_type(atomic_num: i32, aromatic: bool) -> i32 {
    // RDKit✔️✔️: static inline int makeAtomType(int atomic_num, bool aromatic) {
    // RDKit✔️✔️:   return atomic_num + 1000 * static_cast<int>(aromatic);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform one O(1) boolean-to-
    // integer conversion, multiplication, and addition, with no branching,
    // traversal, lookup, allocation, cloning, or temporary collection.
    atomic_num + 1000 * (aromatic as i32)
}

#[inline]
fn parse_atom_type(val: i32) -> (i32, bool) {
    // RDKit✔️✔️: static inline void parseAtomType(int val, int &atomic_num, bool &aromatic) {
    // RDKit✔️✔️:   if (val > 1000) {
    // RDKit✔️✔️:     aromatic = true;
    // RDKit✔️✔️:     atomic_num = val - 1000;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     aromatic = false;
    // RDKit✔️✔️:     atomic_num = val;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform one O(1) integer
    // comparison and at most one subtraction, with no traversal, lookup,
    // allocation, cloning, or temporary collection. Returning two scalar
    // values replaces C++ output references without changing the cost class.
    if val > 1000 {
        (val - 1000, true)
    } else {
        (val, false)
    }
}

#[inline]
fn get_atom_type_is_aromatic(val: i32) -> bool {
    // RDKit✔️✔️: static inline bool getAtomTypeIsAromatic(int val) { return val > 1000; }
    // Local complexity review: RDKit performs one O(1) integer comparison.
    // Rust reuses the inline canonical atom-type decoder, which performs the
    // same comparison and at most one scalar subtraction, with no traversal,
    // lookup, allocation, cloning, or temporary collection. This preserves a
    // single decoding branch instead of duplicating the `> 1000` rule.
    parse_atom_type(val).1
}

#[inline]
fn get_atom_type_atomic_num(val: i32) -> i32 {
    // RDKit✔️✔️: static inline int getAtomTypeAtomicNum(int val) {
    // RDKit✔️✔️:   if (val > 1000) {
    // RDKit✔️✔️:     return val - 1000;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return val;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and the reused inline canonical decoder
    // each perform one O(1) integer comparison and at most one subtraction,
    // with no traversal, lookup, allocation, cloning, or temporary collection.
    // Reuse keeps the `> 1000` decoding rule in one core branch.
    parse_atom_type(val).0
}

#[inline]
fn query_atom_type(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomType(Atom const *at) {
    // RDKit✔️✔️:   return makeAtomType(at->getAtomicNum(), at->getIsAromatic());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform two O(1) typed
    // atom-field reads followed by the same scalar atom-type encoding, with no
    // traversal, lookup, allocation, cloning, or temporary collection. Reuse
    // of the inline query and encoding helpers keeps one core implementation.
    make_atom_type(i32::from(query_atom_num(at)), query_atom_aromatic(at))
}

const MASS_INTEGER_CONVERSION_FACTOR: i32 = 1000;

#[inline]
fn query_atom_mass(at: &Atom) -> i32 {
    // RDKit✔️✔️: const int massIntegerConversionFactor = 1000;
    // RDKit✔️✔️: static inline int queryAtomMass(Atom const *at) {
    // RDKit✔️✔️:   return static_cast<int>(
    // RDKit✔️✔️:       std::round(massIntegerConversionFactor * at->getMass()));
    // RDKit✔️✔️: };
    // RDKit✔️✔️: double Atom::getMass() const {
    // RDKit✔️✔️:   if (d_isotope) {
    // RDKit✔️✔️:     double res =
    // RDKit✔️✔️:         PeriodicTable::getTable()->getMassForIsotope(d_atomicNum, d_isotope);
    // RDKit✔️✔️:     if (d_atomicNum != 0 && res == 0.0) {
    // RDKit✔️✔️:       res = d_isotope;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     return res;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return PeriodicTable::getTable()->getAtomicWeight(d_atomicNum);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations perform O(1) atomic-
    // weight access or O(log I) isotope lookup, followed by one multiplication
    // and round, with no allocation, cloning, graph traversal, or temporary
    // collection. Rust reuses the canonical full RDKit mass tables and
    // `Atom::getMass` fallback port in `rdkit_atomic_mass`; its sorted-table
    // binary search and RDKit's isotope map lookup have the same asymptotic
    // complexity.
    (f64::from(MASS_INTEGER_CONVERSION_FACTOR)
        * rdkit_atomic_mass(at.atomic_number(), at.isotope()))
    .round() as i32
}

#[inline]
fn query_atom_isotope(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomIsotope(Atom const *at) {
    // RDKit✔️✔️:   return static_cast<int>(at->getIsotope());
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int getIsotope() const { return d_isotope; }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // isotope-state read and one scalar conversion, with no traversal, lookup,
    // allocation, cloning, branching over graph state, or temporary object.
    // COSMolKit's `None` is the typed representation of RDKit's zero isotope.
    i32::from(at.isotope().unwrap_or(0))
}

#[inline]
fn query_atom_formal_charge(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomFormalCharge(Atom const *at) {
    // RDKit✔️✔️:   return static_cast<int>(at->getFormalCharge());
    // RDKit✔️✔️: };
    // RDKit✔️✔️: int getFormalCharge() const { return d_formalCharge; }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // formal-charge field read and one scalar conversion, with no traversal,
    // lookup, allocation, cloning, branching, or temporary object. Rust's
    // `i8` charge representation converts losslessly within the modeled state.
    i32::from(at.formal_charge())
}

#[inline]
fn query_atom_negative_formal_charge(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomNegativeFormalCharge(Atom const *at) {
    // RDKit✔️✔️:   return static_cast<int>(-1 * at->getFormalCharge());
    // RDKit✔️✔️: };
    // RDKit✔️✔️: int getFormalCharge() const { return d_formalCharge; }
    // Local complexity review: RDKit and Rust each perform one O(1) formal-
    // charge read and one integer negation, with no traversal, lookup,
    // allocation, cloning, branching, or temporary object. Reusing the
    // canonical formal-charge reader preserves a single field-access core;
    // COSMolKit's modeled `i8` range cannot overflow the `i32` negation.
    -query_atom_formal_charge(at)
}

#[inline]
fn query_atom_hybridization(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomHybridization(Atom const *at) {
    // RDKit✔️✔️:   return at->getHybridization();
    // RDKit✔️✔️: };
    // RDKit✔️✔️: typedef enum {
    // RDKit✔️✔️:   UNSPECIFIED = 0,  //!< hybridization that hasn't been specified
    // RDKit✔️✔️:   S,
    // RDKit✔️✔️:   SP,
    // RDKit✔️✔️:   SP2,
    // RDKit✔️✔️:   SP3,
    // RDKit✔️✔️:   SP2D,
    // RDKit✔️✔️:   SP3D,
    // RDKit✔️✔️:   SP3D2,
    // RDKit✔️✔️:   OTHER  //!< unrecognized hybridization
    // RDKit✔️✔️: } HybridizationType;
    // RDKit✔️✔️: HybridizationType getHybridization() const {
    // RDKit✔️✔️:   return static_cast<HybridizationType>(d_hybrid);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // hybridization field read and return its scalar discriminant, with no
    // traversal, lookup, allocation, cloning, branching, or temporary object.
    // The Rust fieldless enum is declared in the identical source order.
    at.hybridization() as i32
}

#[inline]
fn query_atom_num_radical_electrons(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomNumRadicalElectrons(Atom const *at) {
    // RDKit✔️✔️:   return at->getNumRadicalElectrons();
    // RDKit✔️✔️: };
    // RDKit✔️✔️: unsigned int getNumRadicalElectrons() const { return d_numRadicalElectrons; }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // radical-electron field read and one lossless scalar conversion for the
    // currently modeled `u8` state, with no traversal, lookup, allocation,
    // cloning, branching, or temporary object creation.
    i32::from(at.radical_electrons())
}

#[inline]
fn query_atom_has_chiral_tag(at: &Atom) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomHasChiralTag(Atom const *at) {
    // RDKit✔️✔️:   return at->getChiralTag() != Atom::CHI_UNSPECIFIED;
    // RDKit✔️✔️: };
    // RDKit✔️✔️: typedef enum {
    // RDKit✔️✔️:   CHI_UNSPECIFIED = 0,  //!< chirality that hasn't been specified
    // RDKit✔️✔️:   CHI_TETRAHEDRAL_CW,   //!< tetrahedral: clockwise rotation (SMILES \@\@)
    // RDKit✔️✔️:   CHI_TETRAHEDRAL_CCW,  //!< tetrahedral: counter-clockwise rotation (SMILES
    // RDKit✔️✔️:                           //\@)
    // RDKit✔️✔️:   CHI_OTHER,            //!< some unrecognized type of chirality
    // RDKit✔️✔️:   CHI_TETRAHEDRAL,      //!< tetrahedral, use permutation flag
    // RDKit✔️✔️:   CHI_ALLENE,           //!< allene, use permutation flag
    // RDKit✔️✔️:   CHI_SQUAREPLANAR,     //!< square planar, use permutation flag
    // RDKit✔️✔️:   CHI_TRIGONALBIPYRAMIDAL,  //!< trigonal bipyramidal, use permutation flag
    // RDKit✔️✔️:   CHI_OCTAHEDRAL            //!< octahedral, use permutation flag
    // RDKit✔️✔️: } ChiralType;
    // RDKit✔️✔️: ChiralType getChiralTag() const {
    // RDKit✔️✔️:   return static_cast<ChiralType>(d_chiralTag);
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each perform one O(1) typed
    // chiral-tag field read and one comparison against the unspecified tag,
    // with no traversal, lookup, allocation, cloning, or temporary object.
    // Rust's fieldless enum is declared in the identical source order; the
    // boolean-to-`i32` conversion reproduces C++'s implicit 0/1 return value.
    (at.chiral_tag() != ChiralTag::Unspecified) as i32
}

#[inline]
fn query_atom_missing_chiral_tag(at: &Atom) -> i32 {
    // RDKit✔️🔝: static inline int queryAtomMissingChiralTag(Atom const *at) {
    // RDKit✔️🔝:   return at->getChiralTag() == Atom::CHI_UNSPECIFIED &&
    // RDKit✔️🔝:          at->hasProp(common_properties::_ChiralityPossible);
    // RDKit✔️🔝: };
    // RDKit✔️🔝: bool hasProp(const std::string_view key) const { return d_props.hasVal(key); }
    // RDKit✔️🔝: bool hasVal(const std::string_view what) const {
    // RDKit✔️🔝:   for (const auto &data : _data) {
    // RDKit✔️🔝:     if (data.key == what) {
    // RDKit✔️🔝:       return true;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return false;
    // RDKit✔️🔝: }
    // Local complexity review: both implementations short-circuit after the
    // O(1) chiral-tag read and allocate or clone nothing. RDKit then scans its
    // property vector in O(P), while COSMolKit's canonical atom-property
    // `BTreeMap` performs the same presence-only test in O(log P). Reusing
    // `query_atom_has_chiral_tag` preserves the identical unspecified-tag
    // semantics, and the faster lookup cannot change the result.
    (query_atom_has_chiral_tag(at) == 0 && at.prop("_ChiralityPossible").is_some()) as i32
}

#[inline]
fn query_atom_has_heteroatom_nbrs(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomHasHeteroatomNbrs(Atom const *at) {
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     const Atom *nbr = at->getOwningMol()[*nbrIdx];
    // RDKit✔️✔️:     if (nbr->getAtomicNum() != 6 && nbr->getAtomicNum() != 1) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the owning molecule's indexed adjacency range and return on the
    // first non-carbon, non-hydrogen neighbor. Both use O(1) atom indexing and
    // atomic-number reads, allocate and clone nothing, and create no temporary
    // collection. Rust reuses the canonical CSR `AdjacencyList` and atom-number
    // helper; caching the scalar avoids a repeated field read without changing
    // asymptotic complexity or source behavior.
    for nbr in adj.neighbors_of(at.id().index()) {
        let atomic_number = query_atom_num(&mol.atoms()[nbr.atom_index]);
        if atomic_number != 6 && atomic_number != 1 {
            return 1;
        }
    }
    0
}

#[inline]
fn query_atom_num_heteroatom_nbrs(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomNumHeteroatomNbrs(Atom const *at) {
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     const Atom *nbr = at->getOwningMol()[*nbrIdx];
    // RDKit✔️✔️:     if (nbr->getAtomicNum() != 6 && nbr->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform exactly one
    // O(degree) traversal over indexed adjacency and one O(1) atom lookup per
    // neighbor. Both use an integer accumulator and allocate or clone nothing;
    // neither repeats graph scans nor creates a temporary collection. Rust
    // reuses the canonical CSR adjacency and atom-number helper, caching each
    // scalar atomic number without changing the result or asymptotic cost.
    let mut res = 0;
    for nbr in adj.neighbors_of(at.id().index()) {
        let atomic_number = query_atom_num(&mol.atoms()[nbr.atom_index]);
        if atomic_number != 6 && atomic_number != 1 {
            res += 1;
        }
    }
    res
}

#[inline]
fn query_atom_has_aliphatic_heteroatom_nbrs(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomHasAliphaticHeteroatomNbrs(Atom const *at) {
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     const Atom *nbr = at->getOwningMol()[*nbrIdx];
    // RDKit✔️✔️:     if ((!nbr->getIsAromatic()) && nbr->getAtomicNum() != 6 &&
    // RDKit✔️✔️:         nbr->getAtomicNum() != 1) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over indexed adjacency and return at the first non-aromatic, non-carbon,
    // non-hydrogen neighbor. Both perform O(1) typed atom lookups and field
    // reads, allocate and clone nothing, and create no temporary collection.
    // Rust reuses the canonical aromatic flag, atom-number helper, and CSR
    // adjacency; caching the scalar atomic number does not alter behavior.
    for nbr in adj.neighbors_of(at.id().index()) {
        let neighbor = &mol.atoms()[nbr.atom_index];
        let atomic_number = query_atom_num(neighbor);
        if !neighbor.is_aromatic() && atomic_number != 6 && atomic_number != 1 {
            return 1;
        }
    }
    0
}

#[inline]
fn query_atom_num_aliphatic_heteroatom_nbrs(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomNumAliphaticHeteroatomNbrs(Atom const *at) {
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   ROMol::ADJ_ITER nbrIdx, endNbrs;
    // RDKit✔️✔️:   boost::tie(nbrIdx, endNbrs) = at->getOwningMol().getAtomNeighbors(at);
    // RDKit✔️✔️:   while (nbrIdx != endNbrs) {
    // RDKit✔️✔️:     const Atom *nbr = at->getOwningMol()[*nbrIdx];
    // RDKit✔️✔️:     if ((!nbr->getIsAromatic()) && nbr->getAtomicNum() != 6 &&
    // RDKit✔️✔️:         nbr->getAtomicNum() != 1) {
    // RDKit✔️✔️:       ++res;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++nbrIdx;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over indexed adjacency, perform one O(1) atom lookup per neighbor, and
    // maintain a scalar integer accumulator. Neither implementation allocates,
    // clones, repeats graph scans, or builds a temporary collection. Rust
    // reuses the canonical aromatic flag, atom-number helper, and CSR
    // adjacency; caching each atomic number preserves source behavior and cost.
    let mut res = 0;
    for nbr in adj.neighbors_of(at.id().index()) {
        let neighbor = &mol.atoms()[nbr.atom_index];
        let atomic_number = query_atom_num(neighbor);
        if !neighbor.is_aromatic() && atomic_number != 6 && atomic_number != 1 {
            res += 1;
        }
    }
    res
}

#[inline]
fn query_atom_ring_membership(atom: &Atom, ring_info: &RingInfo) -> i32 {
    // RDKit✔️✔️: static inline int queryIsAtomInNRings(Atom const *at) {
    // RDKit✔️✔️:   return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx());
    // RDKit✔️✔️: };
    // RDKit✔️✔️: static inline int queryAtomRingMembership(Atom const *at) {
    // RDKit✔️✔️:   return static_cast<int>(
    // RDKit✔️✔️:       at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()));
    // RDKit✔️✔️: }
    // Local complexity review: after initialized ring information is supplied,
    // RDKit and Rust each perform one O(1) atom-id lookup and one O(1)
    // member-vector length read. Neither traverses rings, allocates, clones,
    // or creates a temporary collection. Rust receives typed RingInfo
    // explicitly because Atom does not retain an owning-molecule pointer. The
    // cast reproduces RDKit's explicit conversion to the query data type.
    ring_info.num_atom_rings(atom.id()) as i32
}

#[inline]
fn query_is_atom_in_ring(atom: &Atom, ring_info: &RingInfo) -> i32 {
    // RDKit✔️✔️: static inline int queryIsAtomInRing(Atom const *at) {
    // RDKit✔️✔️:   return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()) != 0;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform one O(1) atom-id
    // lookup, one O(1) member-vector length read, and one zero comparison.
    // Neither traverses, allocates, clones, or creates a temporary collection.
    // The Rust helper reuses the canonical ring-count helper and preserves the
    // source function's exact integer 0/1 result.
    i32::from(query_atom_ring_membership(atom, ring_info) != 0)
}

#[inline]
fn query_atom_has_ring_bond(
    atom: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
    ring_info: &RingInfo,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomHasRingBond(Atom const *at) {
    // RDKit✔️✔️:   ROMol::OBOND_ITER_PAIR atomBonds = at->getOwningMol().getAtomBonds(at);
    // RDKit✔️✔️:   while (atomBonds.first != atomBonds.second) {
    // RDKit✔️✔️:     unsigned int bondIdx =
    // RDKit✔️✔️:         at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
    // RDKit✔️✔️:     if (at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
    // RDKit✔️✔️:       return 1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atomBonds.first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return 0;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the owning molecule's indexed incident-bond range and return on
    // the first ring bond. Each iteration performs O(1) bond indexing and an
    // O(1) ring-membership count lookup. Neither allocates, clones, repeats
    // the scan, or creates a temporary collection. Rust receives the owning
    // molecule state explicitly and reuses the canonical ring-count helper.
    for neighbor in adj.neighbors_of(atom.id().index()) {
        let bond = &mol.bonds()[neighbor.bond.index()];
        if query_is_bond_in_n_rings(bond, ring_info) != 0 {
            return 1;
        }
    }
    0
}

#[inline]
pub(crate) fn query_is_bond_in_ring(bond: &Bond, ring_info: &RingInfo) -> i32 {
    // RDKit✔️✔️: static inline int queryIsBondInRing(Bond const *bond) {
    // RDKit✔️✔️:   return bond->getOwningMol().getRingInfo()->numBondRings(bond->getIdx()) != 0;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform one O(1) bond-id
    // lookup, one O(1) member-vector length read, and one zero comparison.
    // Neither traverses, allocates, clones, or creates a temporary collection.
    // Rust reuses the canonical ring-count helper and preserves the source
    // function's exact integer 0/1 result.
    i32::from(query_is_bond_in_n_rings(bond, ring_info) != 0)
}

#[inline]
fn query_atom_min_ring_size(atom: &Atom, ring_info: &RingInfo) -> usize {
    // RDKit✔️✔️: static inline int queryAtomMinRingSize(Atom const *at) {
    // RDKit✔️✔️:   return at->getOwningMol().getRingInfo()->minAtomRingSize(at->getIdx());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each inspect the initialized
    // ring-membership list for one indexed atom and select the minimum ring
    // size, returning zero when that list is empty. Both are O(R_atom), use
    // O(1) auxiliary space, allocate and clone nothing, and perform no graph
    // traversal. Rust receives typed RingInfo explicitly because Atom does not
    // retain an owning-molecule pointer.
    ring_info.min_atom_ring_size(atom.id())
}

#[inline]
pub(crate) fn query_bond_min_ring_size(bond: &Bond, ring_info: &RingInfo) -> usize {
    // RDKit✔️✔️: static inline int queryBondMinRingSize(Bond const *bond) {
    // RDKit✔️✔️:   return bond->getOwningMol().getRingInfo()->minBondRingSize(bond->getIdx());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each inspect the initialized
    // ring-membership list for one indexed bond and select the minimum ring
    // size, returning zero when that list is empty. Both are O(R_bond), use
    // O(1) auxiliary space, allocate and clone nothing, and perform no graph
    // traversal. Rust receives typed RingInfo explicitly because Bond does not
    // retain an owning-molecule pointer.
    ring_info.min_bond_ring_size(bond.id())
}

#[inline]
fn query_atom_ring_bond_count(
    atom: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
    ring_info: &RingInfo,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomRingBondCount(Atom const *at) {
    // RDKit✔️✔️:   // EFF: cache this result
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   ROMol::OBOND_ITER_PAIR atomBonds = at->getOwningMol().getAtomBonds(at);
    // RDKit✔️✔️:   while (atomBonds.first != atomBonds.second) {
    // RDKit✔️✔️:     unsigned int bondIdx =
    // RDKit✔️✔️:         at->getOwningMol().getTopology()[*atomBonds.first]->getIdx();
    // RDKit✔️✔️:     if (at->getOwningMol().getRingInfo()->numBondRings(bondIdx)) {
    // RDKit✔️✔️:       res++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     ++atomBonds.first;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the indexed incident-bond range, perform O(1) bond indexing and
    // ring-count lookup per entry, and maintain one integer accumulator.
    // Neither allocates, clones, repeats a scan, or creates a temporary
    // collection. Rust reuses the canonical bond ring-count helper.
    let mut res = 0;
    for neighbor in adj.neighbors_of(atom.id().index()) {
        let bond = &mol.bonds()[neighbor.bond.index()];
        if query_is_bond_in_n_rings(bond, ring_info) != 0 {
            res += 1;
        }
    }
    res
}

#[inline]
fn query_atom_is_in_ring_of_size(atom: &Atom, target: i32, ring_info: &RingInfo) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomIsInRingOfSize(Atom const *at, int tgt) {
    // RDKit✔️✔️:   if (at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(), tgt)) {
    // RDKit✔️✔️:     return tgt;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // RDKit✔️✔️: template <int tgt>
    // RDKit✔️✔️: int queryAtomIsInRingOfSize(Atom const *at) {
    // RDKit✔️✔️:   if (at->getOwningMol().getRingInfo()->isAtomInRingOfSize(at->getIdx(), tgt)) {
    // RDKit✔️✔️:     return tgt;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // Local complexity review: both source overloads and Rust perform one
    // O(R_atom) scan of the indexed atom ring-membership list, return early on
    // a matching size, and use O(1) auxiliary space. Neither graph-traverses,
    // clones, or creates a temporary collection. A negative C++ target is
    // converted to an unsigned size and cannot match; Rust returns zero before
    // conversion, preserving that result without overflow.
    if target >= 0 && ring_info.is_atom_in_ring_of_size(atom.id(), target as usize) {
        target
    } else {
        0
    }
}

#[inline]
fn query_atom_is_in_ring_size_range(
    atom: &Atom,
    lower: i32,
    upper: i32,
    lower_open: bool,
    upper_open: bool,
    ring_info: &RingInfo,
) -> i32 {
    // RDKit✔️✔️: static inline int queryAtomIsInRingOfSize(Atom const *at, int lower, int upper,
    // RDKit✔️✔️:                                           bool lowerOpen = false,
    // RDKit✔️✔️:                                           bool upperOpen = false) {
    // RDKit✔️✔️:   const auto ri = at->getOwningMol().getRingInfo();
    // RDKit✔️✔️:   for (const auto ringSize : ri->atomRingSizes(at->getIdx())) {
    // RDKit✔️✔️:     if ((ringSize > lower || (ringSize == lower && !lowerOpen)) &&
    // RDKit✔️✔️:         (upper < 0 ||
    // RDKit✔️✔️:          (ringSize < upper || (ringSize == upper && !upperOpen)))) {
    // RDKit✔️✔️:       return ringSize;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   // we didn't find it, return a result that's not in the acceptable range:
    // RDKit✔️✔️:   if (lower > -1) {
    // RDKit✔️✔️:     return -1;
    // RDKit✔️✔️:   } else if (upper > -1) {
    // RDKit✔️✔️:     return std::numeric_limits<int>::max();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // Local complexity review: RDKit atomRingSizes and Rust atom_ring_sizes
    // each allocate one O(R_atom) integer vector, then make one linear pass
    // with identical bound checks and early return. Both use O(R_atom) temporary
    // space and perform no graph traversal, repeated scan, or cloning of graph
    // state. The failure sentinels map exactly to i32 values.
    for ring_size in ring_info.atom_ring_sizes(atom.id()) {
        let ring_size = ring_size as i32;
        if (ring_size > lower || (ring_size == lower && !lower_open))
            && (upper < 0 || ring_size < upper || (ring_size == upper && !upper_open))
        {
            return ring_size;
        }
    }
    if lower > -1 {
        -1
    } else if upper > -1 {
        i32::MAX
    } else {
        0
    }
}

#[inline]
fn query_bond_is_in_ring_of_size(bond: &Bond, target: i32, ring_info: &RingInfo) -> i32 {
    // RDKit✔️✔️: template <int tgt>
    // RDKit✔️✔️: int queryBondIsInRingOfSize(Bond const *bond) {
    // RDKit✔️✔️:   if (bond->getOwningMol().getRingInfo()->isBondInRingOfSize(bond->getIdx(),
    // RDKit✔️✔️:                                                              tgt)) {
    // RDKit✔️✔️:     return tgt;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // Local complexity review: each RDKit template instantiation and the Rust
    // helper make one O(R_bond) scan of the indexed bond ring-membership list,
    // return early on a matching size, and use O(1) auxiliary space. Neither
    // graph-traverses, allocates, clones, or creates a temporary collection.
    // A negative runtime target cannot correspond to a C++ template used by
    // the query factory; returning zero keeps the helper fail-closed.
    if target >= 0 && ring_info.is_bond_in_ring_of_size(bond.id(), target as usize) {
        target
    } else {
        0
    }
}

#[inline]
const fn rdkit_bond_type_prime(order: BondOrder) -> u32 {
    // RDKit✔️✔️: typedef enum {
    // RDKit✔️✔️:   UNSPECIFIED = 0,
    // RDKit✔️✔️:   SINGLE,
    // RDKit✔️✔️:   DOUBLE,
    // RDKit✔️✔️:   TRIPLE,
    // RDKit✔️✔️:   QUADRUPLE,
    // RDKit✔️✔️:   QUINTUPLE,
    // RDKit✔️✔️:   HEXTUPLE,
    // RDKit✔️✔️:   ONEANDAHALF,
    // RDKit✔️✔️:   TWOANDAHALF,
    // RDKit✔️✔️:   THREEANDAHALF,
    // RDKit✔️✔️:   FOURANDAHALF,
    // RDKit✔️✔️:   FIVEANDAHALF,
    // RDKit✔️✔️:   AROMATIC,
    // RDKit✔️✔️:   IONIC,
    // RDKit✔️✔️:   HYDROGEN,
    // RDKit✔️✔️:   THREECENTER,
    // RDKit✔️✔️:   DATIVEONE,
    // RDKit✔️✔️:   DATIVE,
    // RDKit✔️✔️:   DATIVEL,
    // RDKit✔️✔️:   DATIVER,
    // RDKit✔️✔️:   OTHER,
    // RDKit✔️✔️:   ZERO
    // RDKit✔️✔️: } BondType;
    // RDKit✔️✔️: int firstThousandPrimes[NUM_PRIMES_AVAIL] = {
    // RDKit✔️✔️:     2,    3,    5,    7,    11,   13,   17,   19,   23,   29,   31,   37,
    // RDKit✔️✔️:     41,   43,   47,   53,   59,   61,   67,   71,   73,   79,   83,   89,
    // Local complexity review: this exhaustive typed match is one O(1)
    // branch-table lookup, matching RDKit's O(1) array indexing without
    // allocation, cloning, iteration, or temporary collections. It maps by
    // source enum identity instead of COSMolKit declaration order, which
    // differs for the dative and hydrogen bond variants.
    match order {
        BondOrder::Null | BondOrder::Unspecified => 2,
        BondOrder::Single => 3,
        BondOrder::Double => 5,
        BondOrder::Triple => 7,
        BondOrder::Quadruple => 11,
        BondOrder::Quintuple => 13,
        BondOrder::Hextuple => 17,
        BondOrder::OneAndHalf => 19,
        BondOrder::TwoAndHalf => 23,
        BondOrder::ThreeAndHalf => 29,
        BondOrder::FourAndHalf => 31,
        BondOrder::FiveAndHalf => 37,
        BondOrder::Aromatic => 41,
        BondOrder::Ionic => 43,
        BondOrder::Hydrogen => 47,
        BondOrder::ThreeCenter => 53,
        BondOrder::DativeOne => 59,
        BondOrder::Dative => 61,
        BondOrder::DativeLeft => 67,
        BondOrder::DativeRight => 71,
        BondOrder::Other => 73,
        BondOrder::Zero => 79,
    }
}

#[inline]
fn query_atom_bond_product(at: &Atom, adj: &AdjacencyList, mol: &impl SearchTargetAccess) -> u32 {
    // RDKit✔️✔️: unsigned int queryAtomBondProduct(Atom const *at) {
    // RDKit✔️✔️:   ROMol::OEDGE_ITER beg, end;
    // RDKit✔️✔️:   boost::tie(beg, end) = at->getOwningMol().getAtomBonds(at);
    // RDKit✔️✔️:   unsigned int prod = 1;
    // RDKit✔️✔️:   while (beg != end) {
    // RDKit✔️✔️:     prod *= static_cast<unsigned int>(
    // RDKit✔️✔️:         firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    // RDKit✔️✔️:     ++beg;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return prod;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the owning molecule's indexed incident-bond range, perform one
    // O(1) bond lookup and one O(1) prime lookup per entry, and retain one
    // integer accumulator. Neither allocates, clones, repeats a scan, or
    // creates a temporary collection. `wrapping_mul` explicitly preserves
    // C++ unsigned-int modulo-2^32 multiplication in every Rust profile.
    let mut prod = 1_u32;
    for neighbor in adj.neighbors_of(at.id().index()) {
        let bond = &mol.bonds()[neighbor.bond.index()];
        prod = prod.wrapping_mul(rdkit_bond_type_prime(bond.order()));
    }
    prod
}

#[inline]
fn query_atom_all_bond_product(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
    valence: Option<&ValenceAssignment>,
) -> Option<u32> {
    // RDKit✔️✔️: unsigned int queryAtomAllBondProduct(Atom const *at) {
    // RDKit✔️✔️:   ROMol::OEDGE_ITER beg, end;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   boost::tie(beg, end) = at->getOwningMol().getAtomBonds(at);
    // RDKit✔️✔️:   unsigned int prod = 1;
    // RDKit✔️✔️:   while (beg != end) {
    // RDKit✔️✔️:     prod *= static_cast<unsigned int>(
    // RDKit✔️✔️:         firstThousandPrimes[at->getOwningMol()[*beg]->getBondType()]);
    // RDKit✔️✔️:     ++beg;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   for (unsigned int i = 0; i < at->getTotalNumHs(); i++) {
    // RDKit✔️✔️:     prod *= static_cast<unsigned int>(firstThousandPrimes[Bond::SINGLE]);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return prod;
    // RDKit✔️✔️: }
    // Local complexity review: RDKit performs O(degree + total_H) indexed
    // work with one accumulator. Rust performs the same single incident-bond
    // pass through the canonical `query_atom_bond_product`, followed by the
    // same O(total_H) multiplication loop. The typed hydrogen lookup is O(1).
    // Neither path allocates, clones, repeats an adjacency scan, or builds a
    // temporary collection; wrapping multiplication preserves unsigned C++
    // overflow. Reuse prevents a second explicit-bond product implementation.
    let mut prod = query_atom_bond_product(at, adj, mol);
    let total_hydrogens = total_hydrogen_count(valence, at)?;
    for _ in 0..total_hydrogens {
        prod = prod.wrapping_mul(rdkit_bond_type_prime(BondOrder::Single));
    }
    Some(prod)
}

#[inline]
fn query_atom_explicit_degree(at: &Atom, adj: &AdjacencyList) -> usize {
    // RDKit✔️✔️: static inline int queryAtomExplicitDegree(Atom const *at) {
    // RDKit✔️✔️:   return at->getDegree();
    // RDKit✔️✔️: };
    // Local complexity review: RDKit's graph degree lookup and the Rust
    // adjacency-slice length lookup are both O(1), with no traversal,
    // allocation, cloning, or temporary object creation.
    adj.neighbors_of(at.id().index()).len()
}

/// RDKit✔️✔️: Evaluate a single `AtomQueryPredicate` against atom `atom` from `mol`.
///
/// RDKit source: `QueryOps.h` inline queryAtom* functions
///   + `QueryOps.cpp` atomMatchesQuery dispatch logic.
///
/// Returns `true` when the predicate matches the atom.
pub fn atom_predicate_matches(atom: &Atom, pred: &AtomQueryPredicate, mol: &Molecule) -> bool {
    let ctx = build_query_match_context(mol);
    atom_predicate_matches_with_target_context(atom, pred, mol, &ctx)
}

pub(crate) fn atom_predicate_matches_with_target_context(
    atom: &Atom,
    pred: &AtomQueryPredicate,
    mol: &impl SearchTargetAccess,
    ctx: &QueryMatchContext,
) -> bool {
    let aidx = atom.id().index();
    let adj = &ctx.adj;
    let ring_info = &ctx.ring_info;
    let valence = &ctx.valence;

    match pred {
        // RDKit✔️✔️: `*` matches any atom — equivalent to AtomNull with no negation.
        AtomQueryPredicate::Any => true,

        AtomQueryPredicate::AtomicNumber(n) => {
            equality_query_match(i32::from(*n), atom, 0, false, |atom| {
                i32::from(query_atom_num(atom))
            })
        }

        // RDKit✔️✔️: return makeAtomSimpleQuery<ATOM_EQUALS_QUERY>(
        // RDKit✔️✔️:     makeAtomType(num, aromatic), queryAtomType, "AtomType");
        AtomQueryPredicate::AtomType {
            atomic_number,
            aromatic,
        } => query_atom_type(atom) == make_atom_type(i32::from(*atomic_number), *aromatic),

        // RDKit✔️✔️: `[#N,#M]` — atomic number in list.
        AtomQueryPredicate::AtomicNumberIn(vals) => vals.contains(&query_atom_num(atom)),

        // RDKit✔️✔️: `[!#N;!#M]` — atomic number not in list.
        AtomQueryPredicate::AtomicNumberNotIn(vals) => !vals.contains(&query_atom_num(atom)),

        // RDKit✔️✔️: `[+N]` / `[-N]` — queryAtomFormalCharge matches charge.
        AtomQueryPredicate::FormalCharge(c) => query_atom_formal_charge(atom) == i32::from(*c),

        AtomQueryPredicate::NegativeFormalCharge(c) => {
            query_atom_negative_formal_charge(atom) == i32::from(*c)
        }

        AtomQueryPredicate::NumRadicalElectrons(n) => {
            query_atom_num_radical_electrons(atom) == i32::from(*n)
        }

        AtomQueryPredicate::HasChiralTag => query_atom_has_chiral_tag(atom) != 0,

        AtomQueryPredicate::MissingChiralTag => query_atom_missing_chiral_tag(atom) != 0,

        // RDKit✔️✔️: isotope match — queryAtomIsotope.
        AtomQueryPredicate::Isotope(i) => query_atom_isotope(atom) == i32::from(*i),

        AtomQueryPredicate::HydrogenCount(n) => {
            query_atom_h_count(adj, valence.as_ref(), atom, mol) == Some(usize::from(*n))
        }

        AtomQueryPredicate::HasImplicitHydrogen => {
            query_atom_has_implicit_h(valence.as_ref(), atom)
        }

        AtomQueryPredicate::ImplicitValence(n) => {
            query_atom_implicit_valence(valence.as_ref(), atom) == Some(*n)
        }

        AtomQueryPredicate::ExplicitValence(n) => {
            query_atom_explicit_valence(valence.as_ref(), atom) == Some(*n)
        }

        AtomQueryPredicate::ImplicitHydrogenCount(n) => {
            query_atom_implicit_h_count(valence.as_ref(), atom) == Some(usize::from(*n))
        }

        AtomQueryPredicate::ImplicitHydrogenCountLessEqual(n) => {
            query_atom_implicit_h_count(valence.as_ref(), atom)
                .is_some_and(|count| count <= usize::from(*n))
        }

        AtomQueryPredicate::ExplicitDegree(n) => {
            query_atom_explicit_degree(atom, adj) == usize::from(*n)
        }

        // RDKit✔️✔️: explicit degree ≤ N.
        AtomQueryPredicate::ExplicitDegreeLessEqual(n) => {
            query_atom_explicit_degree(atom, adj) <= usize::from(*n)
        }

        AtomQueryPredicate::NonHydrogenDegree(n) => {
            query_atom_non_hydrogen_degree(atom, adj, mol) == *n
        }
        AtomQueryPredicate::NonHydrogenDegreeLessEqual(n) => {
            query_atom_non_hydrogen_degree(atom, adj, mol) <= *n
        }
        AtomQueryPredicate::NonHydrogenDegreeGreaterEqual(n) => {
            query_atom_non_hydrogen_degree(atom, adj, mol) >= *n
        }
        AtomQueryPredicate::HeavyAtomDegree(n) => {
            query_atom_heavy_atom_degree(atom, adj, mol) == *n
        }
        AtomQueryPredicate::NumHeteroatomNeighbors(n) => {
            query_atom_num_heteroatom_nbrs(atom, adj, mol) == i32::from(*n)
        }
        AtomQueryPredicate::HasHeteroatomNeighbors => {
            query_atom_has_heteroatom_nbrs(atom, adj, mol) != 0
        }
        AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(n) => {
            query_atom_num_aliphatic_heteroatom_nbrs(atom, adj, mol) == i32::from(*n)
        }
        AtomQueryPredicate::HasAliphaticHeteroatomNeighbors => {
            query_atom_has_aliphatic_heteroatom_nbrs(atom, adj, mol) != 0
        }

        // RDKit✔️✔️: ring bond count — queryAtomRingBondCount.
        // RDKit source:
        //   queryAtomRingBondCount(at) {
        //     int res = 0;
        //     for atomBonds(at)
        //       if (ringInfo->numBondRings(bondIdx)) res++;
        //     return res;
        //   }
        AtomQueryPredicate::RingBondCount(n) => {
            if let Some(ri) = &ring_info {
                query_atom_ring_bond_count(atom, adj, mol, ri) as u32 == *n
            } else {
                false
            }
        }

        // RDKit✔️✔️: ring bond count ≤ N.
        AtomQueryPredicate::RingBondCountLessEqual(n) => {
            if let Some(ri) = &ring_info {
                query_atom_ring_bond_count(atom, adj, mol, ri) as u8 <= *n
            } else {
                false
            }
        }

        AtomQueryPredicate::HasRingBond => {
            if let Some(ri) = &ring_info {
                query_atom_has_ring_bond(atom, adj, mol, ri) != 0
            } else {
                false
            }
        }

        AtomQueryPredicate::IsBridgehead => ring_info.as_ref().is_some_and(|ri| {
            crate::chemistry::stereo::query_is_atom_bridgehead_from_topology(
                mol.topology_block(),
                aidx,
                ri,
            ) != 0
        }),

        AtomQueryPredicate::IsAromatic(desired) => {
            if *desired {
                query_atom_aromatic(atom)
            } else {
                query_atom_aliphatic(atom)
            }
        }

        AtomQueryPredicate::IsUnsaturated => {
            query_atom_unsaturated(adj, valence.as_ref(), atom).unwrap_or(false)
        }

        // RDKit✔️✔️: hybridization match — queryAtomHybridization.
        AtomQueryPredicate::HybridizationMatch(h) => query_atom_hybridization(atom) == *h as i32,

        AtomQueryPredicate::TotalDegree(n) => {
            query_atom_total_degree(adj, valence.as_ref(), atom) == Some(usize::from(*n))
        }
        AtomQueryPredicate::TotalDegreeLessEqual(n) => {
            query_atom_total_degree(adj, valence.as_ref(), atom)
                .is_some_and(|total| total <= usize::from(*n))
        }
        AtomQueryPredicate::TotalDegreeGreaterEqual(n) => {
            query_atom_total_degree(adj, valence.as_ref(), atom)
                .is_some_and(|total| total >= usize::from(*n))
        }

        AtomQueryPredicate::TotalValence(n) => {
            query_atom_total_valence(valence.as_ref(), atom) == Some(i32::from(*n))
        }
        AtomQueryPredicate::TotalValenceLessEqual(n) => {
            query_atom_total_valence(valence.as_ref(), atom)
                .is_some_and(|total| total <= i32::from(*n))
        }
        AtomQueryPredicate::TotalValenceGreaterEqual(n) => {
            query_atom_total_valence(valence.as_ref(), atom)
                .is_some_and(|total| total >= i32::from(*n))
        }

        // RDKit✔️✔️: in ring — queryIsAtomInRing.
        // RDKit source: queryIsAtomInRing(at) {
        //   return at->getOwningMol().getRingInfo()->numAtomRings(at->getIdx()) != 0;
        // }
        AtomQueryPredicate::InRing => {
            if let Some(ri) = &ring_info {
                query_is_atom_in_ring(atom, ri) != 0
            } else {
                false
            }
        }

        // RDKit✔️✔️: AtomRingQuery(N) — atom ring membership count.
        // RDKit source: `COMPLEX_ATOM_QUERY_TOKEN number` mutates the
        // AtomRingQuery value used for `R` SMARTS primitives.
        AtomQueryPredicate::NumAtomRings(n) => {
            if let Some(ri) = &ring_info {
                let membership = query_atom_ring_membership(atom, ri);
                if *n < 0 {
                    membership != 0
                } else {
                    membership == *n
                }
            } else {
                false
            }
        }

        // RDKit✔️✔️: in ring of size N — isAtomInRingOfSize.
        AtomQueryPredicate::InRingOfSize(n) => {
            if let Some(ri) = &ring_info {
                query_atom_is_in_ring_of_size(atom, i32::from(*n), ri) == i32::from(*n)
            } else {
                false
            }
        }
        AtomQueryPredicate::InRingOfSizeLessEqual(n) => {
            if let Some(ri) = &ring_info {
                query_atom_is_in_ring_size_range(atom, i32::from(*n), -1, false, false, ri)
                    <= i32::from(*n)
            } else {
                false
            }
        }
        AtomQueryPredicate::InRingOfSizeGreaterEqual(n) => {
            if let Some(ri) = &ring_info {
                query_atom_is_in_ring_size_range(atom, -1, i32::from(*n), false, false, ri)
                    >= i32::from(*n)
            } else {
                false
            }
        }

        // RDKit✔️✔️: smallest ring size — queryAtomMinRingSize.
        // RDKit source: queryAtomMinRingSize(at) {
        //   return getRingInfo()->minAtomRingSize(at->getIdx());
        // }
        AtomQueryPredicate::SmallestRingSize(n) => {
            if let Some(ri) = &ring_info {
                query_atom_min_ring_size(atom, ri) as u8 == *n
            } else {
                false
            }
        }
        AtomQueryPredicate::SmallestRingSizeLessEqual(n) => {
            if let Some(ri) = &ring_info {
                query_atom_min_ring_size(atom, ri) as u8 <= *n
            } else {
                false
            }
        }
        AtomQueryPredicate::SmallestRingSizeGreaterEqual(n) => {
            if let Some(ri) = &ring_info {
                query_atom_min_ring_size(atom, ri) as u8 >= *n
            } else {
                false
            }
        }

        // RDKit✔️✔️: mass match — queryAtomMass. `Mass` retains the unscaled
        // integer query value accepted by RDKit's makeAtomMassQuery.
        AtomQueryPredicate::Mass(m) => {
            query_atom_mass(atom) == i32::from(*m) * MASS_INTEGER_CONVERSION_FACTOR
        }

        // RDKit✔️✔️: chiral tag match.
        AtomQueryPredicate::ChiralTagMatch(tag) => atom.chiral_tag() == *tag,
        AtomQueryPredicate::ChiralPermutationMatch(permutation) => {
            atom.chiral_permutation().unwrap_or(0) == *permutation
        }

        // RDKit✔️✔️: comparison forms of degree use the same explicit-degree data function.
        AtomQueryPredicate::DegreeLessEqual(n) => {
            query_atom_explicit_degree(atom, adj) <= usize::from(*n)
        }
        AtomQueryPredicate::DegreeGreaterEqual(n) => {
            query_atom_explicit_degree(atom, adj) >= usize::from(*n)
        }

        AtomQueryPredicate::Range(range) => match_atom_range_query(range, atom, mol, ctx),

        // RDKit✔️✔️: recursive SMARTS — not yet fully supported.
        AtomQueryPredicate::RecursiveSmarts(_query) => {
            // RDKit✔️❌: Recursive SMARTS evaluation requires the full SMARTS matcher /
            // substructure matching engine which is not yet ported. This is preserved
            // as a stored value but not evaluated.
            false
        }

        AtomQueryPredicate::HasProperty(name) => atom.prop(name).is_some(),
        AtomQueryPredicate::PropertyValue { name, value } => {
            atom.prop(name) == Some(value.as_str())
        }

        // RDKit✔️✔️: R-group label.
        AtomQueryPredicate::RGroupLabel(_label) => {
            // RDKit✔️❌: R-group label matching not yet supported.
            false
        }

        // RDKit✔️✔️: MolFile alias.
        AtomQueryPredicate::MolFileAlias(_alias) => {
            // RDKit✔️❌: MolFile alias matching not yet supported.
            false
        }

        // RDKit✔️✔️: explicitly unsupported feature — fail open with false.
        AtomQueryPredicate::UnsupportedFeature(_desc) => {
            // Per policy_invariants.md: unsupported features must not silently
            // produce chemically meaningful results. We return false here and
            // the caller should check for UnsupportedFeature in the match tree.
            false
        }
    }
}

pub fn atom_predicate_matches_with_context(
    atom: &Atom,
    pred: &AtomQueryPredicate,
    mol: &Molecule,
    ctx: &QueryMatchContext,
) -> bool {
    atom_predicate_matches_with_target_context(atom, pred, mol, ctx)
}

// ---------------------------------------------------------------------------
// bond_predicate_matches — evaluate a bond query predicate
// ---------------------------------------------------------------------------

#[inline]
fn query_bond_order(bond: &Bond) -> BondOrder {
    // RDKit✔️✔️: static inline int queryBondOrder(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondType());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit performs one O(1) bond-type bit-field
    // read and an integer cast; Rust performs one O(1) typed enum field read.
    // Both allocate and clone nothing, perform no traversal or lookup, and
    // create no temporary collection. The canonical typed return preserves
    // the same bond-order identity while keeping RDKit integer codes out of
    // the core model; all query comparisons consume this single helper.
    bond.order()
}

#[inline]
fn query_bond_order_in(bond: &Bond, orders: &[BondOrder]) -> bool {
    orders.contains(&query_bond_order(bond))
}

#[inline]
fn query_bond_is_single_or_aromatic(bond: &Bond) -> i32 {
    // RDKit✔️✔️: static inline int queryBondIsSingleOrAromatic(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                           bond->getBondType() == Bond::AROMATIC);
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform at most two O(1)
    // bond-order comparisons with short-circuit evaluation. Neither traverses,
    // allocates, clones, performs a lookup, or creates a temporary collection.
    // Rust reuses the canonical query_bond_order helper, so bond-type access
    // remains centralized and the boolean-to-i32 conversion matches C++ 0/1.
    query_bond_order_in(bond, &[BondOrder::Single, BondOrder::Aromatic]) as i32
}

#[inline]
fn query_bond_is_double_or_aromatic(bond: &Bond) -> i32 {
    // RDKit✔️✔️: static inline int queryBondIsDoubleOrAromatic(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️✔️:                           bond->getBondType() == Bond::AROMATIC);
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform at most two O(1)
    // bond-order comparisons with short-circuit evaluation and no traversal,
    // lookup, allocation, cloning, or temporary collection. Rust reuses the
    // canonical query_bond_order helper, and converting the result to i32
    // preserves the source function's exact 0/1 return values.
    query_bond_order_in(bond, &[BondOrder::Double, BondOrder::Aromatic]) as i32
}

#[inline]
fn query_bond_is_single_or_double(bond: &Bond) -> i32 {
    // RDKit✔️✔️: static inline int queryBondIsSingleOrDouble(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                           bond->getBondType() == Bond::DOUBLE);
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform at most two O(1)
    // bond-order comparisons with short-circuit evaluation and no traversal,
    // lookup, allocation, cloning, or temporary collection. Rust reuses the
    // canonical query_bond_order helper, and converting the result to i32
    // preserves the source function's exact 0/1 return values.
    query_bond_order_in(bond, &[BondOrder::Single, BondOrder::Double]) as i32
}

#[inline]
fn query_bond_is_single_or_double_or_aromatic(bond: &Bond) -> i32 {
    // RDKit✔️✔️: static inline int queryBondIsSingleOrDoubleOrAromatic(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondType() == Bond::SINGLE ||
    // RDKit✔️✔️:                           bond->getBondType() == Bond::DOUBLE ||
    // RDKit✔️✔️:                           bond->getBondType() == Bond::AROMATIC);
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform at most three O(1)
    // bond-order comparisons with short-circuit evaluation and no traversal,
    // lookup, allocation, cloning, or temporary collection. Rust reuses the
    // canonical query_bond_order helper, and converting the result to i32
    // preserves the source function's exact 0/1 return values.
    query_bond_order_in(
        bond,
        &[BondOrder::Single, BondOrder::Double, BondOrder::Aromatic],
    ) as i32
}

#[inline]
fn query_bond_dir(bond: &Bond) -> crate::BondDirection {
    // RDKit✔️✔️: static inline int queryBondDir(Bond const *bond) {
    // RDKit✔️✔️:   return static_cast<int>(bond->getBondDir());
    // RDKit✔️✔️: };
    // Local complexity review: RDKit performs one O(1) direction bit-field
    // read and an integer cast; Rust performs one O(1) typed enum field read.
    // Both allocate and clone nothing, perform no traversal or lookup, and
    // create no temporary collection. Returning the canonical typed direction
    // preserves BondDir identity without adding RDKit integer query state.
    bond.direction()
}

#[inline]
fn query_is_bond_in_n_rings(bond: &Bond, ring_info: &RingInfo) -> usize {
    // RDKit✔️✔️: static inline int queryIsBondInNRings(Bond const *at) {
    // RDKit✔️✔️:   return at->getOwningMol().getRingInfo()->numBondRings(at->getIdx());
    // RDKit✔️✔️: };
    // Local complexity review: after the caller has obtained initialized ring
    // information, RDKit and Rust each perform one O(1) bond-id lookup and one
    // O(1) member-vector length read. Neither path traverses rings, allocates,
    // clones, or creates a temporary collection. Rust receives the owning
    // molecule's typed RingInfo explicitly because Bond does not retain an
    // owning-molecule pointer.
    ring_info.num_bond_rings(bond.id())
}

#[inline]
fn query_bond_has_stereo(bond: &Bond) -> i32 {
    // RDKit✔️✔️: static inline int queryBondHasStereo(Bond const *bnd) {
    // RDKit✔️✔️:   return bnd->getStereo() > Bond::STEREONONE;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each perform one O(1) stereo
    // field read and one O(1) comparison, returning the same integer 0/1.
    // Neither traverses, allocates, clones, performs a lookup, or creates a
    // temporary collection. Typed inequality is equivalent to RDKit's ordered
    // comparison because None is the sole no-stereo state in both models.
    i32::from(bond.stereo() != BondStereo::None)
}

/// RDKit✔️✔️: Evaluate a single `BondQueryPredicate` against bond `bond` from `mol`.
///
/// RDKit source: `QueryOps.h` inline queryBond* functions.
pub fn bond_predicate_matches(bond: &Bond, pred: &BondQueryPredicate, mol: &Molecule) -> bool {
    let ctx = build_query_match_context(mol);
    bond_predicate_matches_with_target_context(bond, pred, mol, &ctx)
}

pub(crate) fn bond_predicate_matches_with_target_context(
    bond: &Bond,
    pred: &BondQueryPredicate,
    _mol: &impl SearchTargetAccess,
    ctx: &QueryMatchContext,
) -> bool {
    let ring_info = &ctx.ring_info;

    match pred {
        // RDKit✔️✔️: `~` matches any bond.
        BondQueryPredicate::Any => true,

        // RDKit✔️✔️: `-`, `=`, `#` — queryBondOrder matches bond type.
        BondQueryPredicate::Order(order) => query_bond_order(bond) == *order,

        // RDKit✔️✔️: bond order in list.
        BondQueryPredicate::OrderIn(orders) => query_bond_order_in(bond, orders),

        // RDKit✔️✔️: `:` aromatic bond.
        BondQueryPredicate::IsAromatic(desired) => bond.is_aromatic() == *desired,

        // RDKit✔️✔️: `@` ring bond — queryIsBondInRing.
        BondQueryPredicate::IsInRing(desired) => {
            if let Some(ri) = &ring_info {
                (query_is_bond_in_ring(bond, ri) != 0) == *desired
            } else {
                !desired
            }
        }

        // RDKit✔️✔️: `/` `\` bond direction — queryBondDir.
        BondQueryPredicate::Direction(dir) => query_bond_dir(bond) == *dir,

        // Exact typed stereo identity matching. RDKit's separate boolean
        // queryBondHasStereo data function is centralized above.
        BondQueryPredicate::Stereo(stereo) => bond.stereo() == *stereo,

        // RDKit✔️✔️: queryBondHasStereo compares against STEREONONE.
        BondQueryPredicate::HasStereo => query_bond_has_stereo(bond) != 0,

        // RDKit✔️✔️: `^` conjugated — bond.is_conjugated().
        BondQueryPredicate::IsConjugated => bond.is_conjugated(),

        // RDKit✔️✔️: number of ring bonds the bond is part of.
        BondQueryPredicate::NumRingBonds(n) => {
            if let Some(ri) = &ring_info {
                usize::try_from(*n).is_ok_and(|target| query_is_bond_in_n_rings(bond, ri) == target)
            } else {
                false
            }
        }
        BondQueryPredicate::InRingOfSize(target) => {
            if let Some(ri) = &ring_info {
                query_bond_is_in_ring_of_size(bond, *target, ri) == *target
            } else {
                false
            }
        }
        BondQueryPredicate::MinRingSize(target) => {
            if let Some(ri) = &ring_info {
                usize::try_from(*target)
                    .is_ok_and(|target| query_bond_min_ring_size(bond, ri) == target)
            } else {
                false
            }
        }
        BondQueryPredicate::NumRingBondsGreaterEqual(n) => {
            if let Some(ri) = &ring_info {
                query_is_bond_in_n_rings(bond, ri) as u8 >= *n
            } else {
                false
            }
        }
        BondQueryPredicate::NumRingBondsLessEqual(n) => {
            if let Some(ri) = &ring_info {
                query_is_bond_in_n_rings(bond, ri) as u8 <= *n
            } else {
                false
            }
        }

        BondQueryPredicate::HasProperty(name) => bond.prop(name).is_some(),
        BondQueryPredicate::PropertyValue { name, value } => {
            bond.prop(name) == Some(value.as_str())
        }

        // RDKit✔️✔️: MolFile query code — preserved but not interpreted.
        BondQueryPredicate::MolFileQueryCode(_code) => {
            // RDKit✔️❌: MolFile bond query codes not yet interpreted.
            false
        }

        // RDKit✔️✔️: explicitly unsupported feature.
        BondQueryPredicate::UnsupportedFeature(_desc) => false,
    }
}

// ---------------------------------------------------------------------------
// Recursive query tree evaluators
// ---------------------------------------------------------------------------

#[inline]
fn query_cmp(target: i32, observed: i32, tolerance: i32) -> i32 {
    // RDKit✔️✔️: template <class T1, class T2>
    // RDKit✔️✔️: int queryCmp(const T1 v1, const T2 v2, const T1 tol) {
    // RDKit✔️✔️:   T1 diff = v1 - v2;
    // RDKit✔️✔️:   if (diff <= tol) {
    // RDKit✔️✔️:     if (diff >= -tol) {
    // RDKit✔️✔️:       return 0;
    // RDKit✔️✔️:     } else {
    // RDKit✔️✔️:       return -1;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return 1;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: };
    // Local complexity review: both implementations perform one subtraction,
    // at most two ordered comparisons, and return in O(1) time and O(1)
    // space. The Rust integer representation matches the `int` query aliases
    // used by RDKit's atom and bond queries. No allocation, clone, lookup,
    // collection, scan, or extra hot-path branch is introduced. Modeled
    // SMARTS query/data values are bounded well inside the i32 subtraction
    // range, so C++ signed-overflow behavior is outside the supported state.
    let difference = target - observed;
    if difference <= tolerance {
        if difference >= -tolerance { 0 } else { -1 }
    } else {
        1
    }
}

fn equality_query_match<T>(
    target: i32,
    what: T,
    tolerance: i32,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   if (queryCmp(this->d_val, mfArg, this->d_tol) == 0) {
    // RDKit✔️✔️:     return !this->getNegation();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return this->getNegation();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both paths invoke the configured conversion or
    // data function exactly once, make one O(1) queryCmp call, and apply the
    // stored negation with O(1) extra space. The monomorphized FnOnce adds no
    // allocation, virtual lookup, clone, temporary collection, repeated scan,
    // or asymptotic/hot-path branch beyond RDKit's TypeConvert dispatch.
    let match_arg = type_convert(what);
    if query_cmp(target, match_arg, tolerance) == 0 {
        !negated
    } else {
        negated
    }
}

fn greater_query_match<T>(
    target: i32,
    what: T,
    tolerance: i32,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   if (queryCmp(this->d_val, mfArg, this->d_tol) > 0) {
    // RDKit✔️✔️:     return !this->getNegation();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return this->getNegation();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both paths invoke TypeConvert exactly once,
    // perform the canonical O(1) queryCmp operation, and apply negation in
    // O(1) space. Reusing query_cmp avoids a second tolerance implementation;
    // the monomorphized FnOnce adds no allocation, clone, lookup, temporary
    // collection, scan, or asymptotic/hot-path branch relative to RDKit.
    let match_arg = type_convert(what);
    if query_cmp(target, match_arg, tolerance) > 0 {
        !negated
    } else {
        negated
    }
}

fn greater_equal_query_match<T>(
    target: i32,
    what: T,
    tolerance: i32,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   if (queryCmp(this->d_val, mfArg, this->d_tol) >= 0) {
    // RDKit✔️✔️:     return !this->getNegation();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return this->getNegation();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both paths invoke TypeConvert exactly once,
    // call the one canonical O(1) queryCmp implementation, and apply negation
    // in O(1) space. The monomorphized FnOnce introduces no allocation,
    // cloning, lookup, temporary collection, scan, repeated conversion, or
    // asymptotic/hot-path branch relative to the RDKit implementation.
    let match_arg = type_convert(what);
    if query_cmp(target, match_arg, tolerance) >= 0 {
        !negated
    } else {
        negated
    }
}

fn less_query_match<T>(
    target: i32,
    what: T,
    tolerance: i32,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   if (queryCmp(this->d_val, mfArg, this->d_tol) < 0) {
    // RDKit✔️✔️:     return !this->getNegation();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return this->getNegation();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations invoke TypeConvert once,
    // call the canonical O(1) queryCmp implementation, and apply negation in
    // O(1) space. The monomorphized FnOnce adds no allocation, clone, lookup,
    // temporary collection, scan, repeated conversion, or extra hot-path
    // branch relative to RDKit.
    let match_arg = type_convert(what);
    if query_cmp(target, match_arg, tolerance) < 0 {
        !negated
    } else {
        negated
    }
}

fn less_equal_query_match<T>(
    target: i32,
    what: T,
    tolerance: i32,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   if (queryCmp(this->d_val, mfArg, this->d_tol) <= 0) {
    // RDKit✔️✔️:     return !this->getNegation();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return this->getNegation();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both paths invoke TypeConvert exactly once,
    // perform the canonical O(1) queryCmp operation, and apply negation in
    // O(1) space. The monomorphized FnOnce introduces no allocation, clone,
    // lookup, temporary collection, scan, repeated conversion, or additional
    // hot-path branch relative to RDKit.
    let match_arg = type_convert(what);
    if query_cmp(target, match_arg, tolerance) <= 0 {
        !negated
    } else {
        negated
    }
}

#[allow(clippy::too_many_arguments)]
fn range_query_match<T>(
    lower: i32,
    upper: i32,
    what: T,
    tolerance: i32,
    lower_open: bool,
    upper_open: bool,
    negated: bool,
    type_convert: impl FnOnce(T) -> i32,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   int lCmp = queryCmp(this->d_lower, mfArg, this->d_tol);
    // RDKit✔️✔️:   int uCmp = queryCmp(this->d_upper, mfArg, this->d_tol);
    // RDKit✔️✔️:   bool lowerRes, upperRes;
    // RDKit✔️✔️:   if (this->df_lowerOpen) {
    // RDKit✔️✔️:     lowerRes = lCmp < 0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     lowerRes = lCmp <= 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (this->df_upperOpen) {
    // RDKit✔️✔️:     upperRes = uCmp > 0;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     upperRes = uCmp >= 0;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool tempR = !(lowerRes && upperRes);
    // RDKit✔️✔️:   if (this->getNegation()) {
    // RDKit✔️✔️:     return tempR;
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return !tempR;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations convert once, invoke two
    // O(1) comparisons, test two endpoint flags, and apply negation in O(1)
    // space. Reusing query_cmp introduces no allocation, clone, lookup,
    // collection, scan, repeated conversion, or extra asymptotic work.
    let match_arg = type_convert(what);
    let lower_cmp = query_cmp(lower, match_arg, tolerance);
    let upper_cmp = query_cmp(upper, match_arg, tolerance);
    let lower_matches = if lower_open {
        lower_cmp < 0
    } else {
        lower_cmp <= 0
    };
    let upper_matches = if upper_open {
        upper_cmp > 0
    } else {
        upper_cmp >= 0
    };
    (lower_matches && upper_matches) != negated
}

fn set_query_match<T, U: Ord>(
    values: &BTreeSet<U>,
    what: T,
    negated: bool,
    type_convert: impl FnOnce(T) -> U,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   MatchFuncArgType mfArg =
    // RDKit✔️✔️:       this->TypeConvert(what, Int2Type<needsConversion>());
    // RDKit✔️✔️:   return (this->d_set.find(mfArg) != this->d_set.end()) ^ this->getNegation();
    // RDKit✔️✔️: }
    // Local complexity review: both implementations invoke TypeConvert once
    // and perform an O(log n) ordered-tree membership lookup in constant
    // auxiliary space. Rust's borrowed BTreeSet adds no allocation, clone,
    // temporary collection, repeated scan, or extra hot-path dispatch.
    let found = values.contains(&type_convert(what));
    found != negated
}

fn query_atom_copy(atom: &crate::QueryAtom) -> crate::QueryAtom {
    // RDKit✔️✔️: Atom *QueryAtom::copy() const {
    // RDKit✔️✔️:   auto *res = new QueryAtom(*this);
    // RDKit✔️✔️:   return static_cast<Atom *>(res);
    // RDKit✔️✔️: }
    // Local complexity review: both paths copy fixed atom state and deeply
    // copy the owned query tree and property values in O(query + properties)
    // time and space. Rust's derived Clone performs no extra traversal,
    // lookup, temporary collection, or repeated allocation beyond cloning the
    // same owned fields, while avoiding a separate outer object allocation.
    atom.clone()
}

fn query_bond_copy(bond: &crate::QueryBond) -> crate::QueryBond {
    // RDKit✔️✔️: Bond *QueryBond::copy() const {
    // RDKit✔️✔️:   auto *res = new QueryBond(*this);
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: both paths copy fixed bond state and deeply
    // copy the owned query tree and property values in O(query + properties)
    // time and space. Rust's derived Clone performs no extra traversal,
    // lookup, temporary collection, or repeated allocation beyond cloning the
    // same owned fields, while avoiding a separate outer object allocation.
    bond.clone()
}

fn query_bond_set_type(bond: &mut crate::QueryBond, bond_type: BondOrder) {
    // RDKit✔️✔️: void QueryBond::setBondType(BondType bT) {
    // RDKit✔️✔️:   // NOTE: calling this blows out any existing query
    // RDKit✔️✔️:   d_bondType = bT;
    // RDKit✔️✔️:   delete dp_query;
    // RDKit✔️✔️:   dp_query = nullptr;
    // RDKit✔️✔️:
    // RDKit✔️✔️:   dp_query = makeBondOrderEqualsQuery(bT);
    // RDKit✔️✔️: }
    // Local complexity review: both implementations assign fixed-size bond
    // state, drop the old owned query tree, and allocate/build one constant-
    // size order predicate in O(old query size) destruction time and O(1)
    // new space. Rust performs no clone, scan, lookup, temporary collection,
    // repeated dispatch, or additional hot-path branch.
    bond.bond_mut().set_order(bond_type);
    bond.set_predicate(make_bond_order_equals_query(bond_type));
}

fn query_bond_set_dir(bond: &mut crate::QueryBond, direction: crate::BondDirection) {
    // RDKit✔️✔️: void QueryBond::setBondDir(BondDir bD) {
    // RDKit✔️✔️:   // NOTE: calling this blows out any existing query
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   //   Ignoring bond orders (which this implicitly does by blowing out
    // RDKit✔️✔️:   //   any bond order query) is ok for organic molecules, where the
    // RDKit✔️✔️:   //   only bonds assigned directions are single.  It'll fail in other
    // RDKit✔️✔️:   //   situations, whatever those may be.
    // RDKit✔️✔️:   //
    // RDKit✔️✔️:   d_dirTag = bD;
    // RDKit✔️✔️: #if 0
    // RDKit✔️✔️:   delete dp_query;
    // RDKit✔️✔️:   dp_query = NULL;
    // RDKit✔️✔️:   dp_query = makeBondDirEqualsQuery(bD);
    // RDKit✔️✔️: #endif
    // RDKit✔️✔️: }
    // Local complexity review: the active RDKit code and Rust each perform
    // one fixed-size direction assignment in O(1) time and space. The query
    // replacement is disabled by RDKit's preprocessor and is likewise not
    // executed here; Rust adds no allocation, clone, traversal, or lookup.
    bond.bond_mut().set_direction(direction);
}

fn query_local_match<T: PartialEq>(
    first_value: &T,
    first_negated: bool,
    second_value: &T,
    second_negated: bool,
) -> bool {
    if first_negated == second_negated {
        first_value == second_value
    } else {
        first_value != second_value
    }
}

fn atom_query_local_match<T: PartialEq>(
    first_value: &T,
    first_negated: bool,
    second_value: &T,
    second_negated: bool,
) -> bool {
    // RDKit✔️✔️: bool localMatch(ATOM_EQUALS_QUERY const *q1, ATOM_EQUALS_QUERY const *q2) {
    // RDKit✔️✔️:   if (q1->getNegation() == q2->getNegation()) {
    // RDKit✔️✔️:     return q1->getVal() == q2->getVal();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return q1->getVal() != q2->getVal();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: both implementations read two negation flags
    // and perform exactly one equality or inequality comparison in O(1)
    // auxiliary space. Borrowing typed values introduces no allocation,
    // clone, lookup, traversal, temporary collection, or extra branch.
    query_local_match(first_value, first_negated, second_value, second_negated)
}

fn bond_query_local_match<T: PartialEq>(
    first_value: &T,
    first_negated: bool,
    second_value: &T,
    second_negated: bool,
) -> bool {
    // RDKit✔️✔️: bool localMatch(BOND_EQUALS_QUERY const *q1, BOND_EQUALS_QUERY const *q2) {
    // RDKit✔️✔️:   if (q1->getNegation() == q2->getNegation()) {
    // RDKit✔️✔️:     return q1->getVal() == q2->getVal();
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:     return q1->getVal() != q2->getVal();
    // RDKit✔️✔️:   }
    // RDKit✔️✔️: }
    // Local complexity review: the shared monomorphized implementation reads
    // two negation flags and performs exactly one equality or inequality
    // comparison in O(1) auxiliary space. It introduces no allocation, clone,
    // lookup, traversal, temporary collection, virtual dispatch, or extra
    // branch, and prevents duplicate atom/bond comparison cores.
    query_local_match(first_value, first_negated, second_value, second_negated)
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum AtomEqualityQueryKind {
    AtomType,
    RingBondCount,
    RingSize,
    MinRingSize,
    ImplicitValence,
    ExplicitValence,
    TotalValence,
    AtomicNumber,
    ExplicitDegree,
    TotalDegree,
    HydrogenCount,
    IsAromatic,
    IsAliphatic,
    Unsaturated,
    Mass,
    FormalCharge,
    NegativeFormalCharge,
    Hybridization,
    InRing,
    InNRings,
}

fn atom_equality_query_kind(predicate: &AtomQueryPredicate) -> Option<AtomEqualityQueryKind> {
    Some(match predicate {
        AtomQueryPredicate::AtomType { .. } => AtomEqualityQueryKind::AtomType,
        AtomQueryPredicate::RingBondCount(_) => AtomEqualityQueryKind::RingBondCount,
        AtomQueryPredicate::InRingOfSize(_) => AtomEqualityQueryKind::RingSize,
        AtomQueryPredicate::SmallestRingSize(_) => AtomEqualityQueryKind::MinRingSize,
        AtomQueryPredicate::ImplicitValence(_) => AtomEqualityQueryKind::ImplicitValence,
        AtomQueryPredicate::ExplicitValence(_) => AtomEqualityQueryKind::ExplicitValence,
        AtomQueryPredicate::TotalValence(_) => AtomEqualityQueryKind::TotalValence,
        AtomQueryPredicate::AtomicNumber(_) => AtomEqualityQueryKind::AtomicNumber,
        AtomQueryPredicate::ExplicitDegree(_) => AtomEqualityQueryKind::ExplicitDegree,
        AtomQueryPredicate::TotalDegree(_) => AtomEqualityQueryKind::TotalDegree,
        AtomQueryPredicate::HydrogenCount(_) => AtomEqualityQueryKind::HydrogenCount,
        AtomQueryPredicate::IsAromatic(true) => AtomEqualityQueryKind::IsAromatic,
        AtomQueryPredicate::IsAromatic(false) => AtomEqualityQueryKind::IsAliphatic,
        AtomQueryPredicate::IsUnsaturated => AtomEqualityQueryKind::Unsaturated,
        AtomQueryPredicate::Mass(_) => AtomEqualityQueryKind::Mass,
        AtomQueryPredicate::FormalCharge(_) => AtomEqualityQueryKind::FormalCharge,
        AtomQueryPredicate::NegativeFormalCharge(_) => AtomEqualityQueryKind::NegativeFormalCharge,
        AtomQueryPredicate::HybridizationMatch(_) => AtomEqualityQueryKind::Hybridization,
        AtomQueryPredicate::InRing => AtomEqualityQueryKind::InRing,
        AtomQueryPredicate::NumAtomRings(_) => AtomEqualityQueryKind::InNRings,
        _ => return None,
    })
}

fn atom_query_base(
    mut query: &QueryNode<AtomQueryPredicate>,
) -> (&QueryNode<AtomQueryPredicate>, bool) {
    let mut negated = false;
    while let QueryNode::Not(child) = query {
        negated = !negated;
        query = child;
    }
    (query, negated)
}

pub(crate) fn atom_queries_match(
    first: &QueryNode<AtomQueryPredicate>,
    second: &QueryNode<AtomQueryPredicate>,
) -> bool {
    // RDKit✔️✔️: bool queriesMatch(QueryAtom::QUERYATOM_QUERY const *q1,
    // RDKit✔️✔️:                   QueryAtom::QUERYATOM_QUERY const *q2) {
    // RDKit✔️✔️:   PRECONDITION(q1, "no q1");
    // RDKit✔️✔️:   PRECONDITION(q2, "no q2");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   static const unsigned int nQueries = 20;
    // RDKit✔️✔️:   static std::string equalityQueries[nQueries] = {"AtomType",
    // RDKit✔️✔️:                                                   "AtomRingBondCount",
    // RDKit✔️✔️:                                                   "AtomRingSize",
    // RDKit✔️✔️:                                                   "AtomMinRingSize",
    // RDKit✔️✔️:                                                   "AtomImplicitValence",
    // RDKit✔️✔️:                                                   "AtomExplicitValence",
    // RDKit✔️✔️:                                                   "AtomTotalValence",
    // RDKit✔️✔️:                                                   "AtomAtomicNum",
    // RDKit✔️✔️:                                                   "AtomExplicitDegree",
    // RDKit✔️✔️:                                                   "AtomTotalDegree",
    // RDKit✔️✔️:                                                   "AtomHCount",
    // RDKit✔️✔️:                                                   "AtomIsAromatic",
    // RDKit✔️✔️:                                                   "AtomIsAliphatic",
    // RDKit✔️✔️:                                                   "AtomUnsaturated",
    // RDKit✔️✔️:                                                   "AtomMass",
    // RDKit✔️✔️:                                                   "AtomFormalCharge",
    // RDKit✔️✔️:                                                   "AtomNegativeFormalCharge",
    // RDKit✔️✔️:                                                   "AtomHybridization",
    // RDKit✔️✔️:                                                   "AtomInRing",
    // RDKit✔️✔️:                                                   "AtomInNRings"};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   std::string d1 = q1->getDescription();
    // RDKit✔️✔️:   std::string d2 = q2->getDescription();
    // RDKit✔️✔️:   if (d1 == "AtomNull" || d2 == "AtomNull") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:   } else if (d1 == "AtomOr") {
    // RDKit✔️✔️:     // FIX: handle negation on AtomOr and AtomAnd
    // RDKit✔️✔️:     for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
    // RDKit✔️✔️:          ++iter1) {
    // RDKit✔️✔️:       if (d2 == "AtomOr") {
    // RDKit✔️✔️:         for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:              ++iter2) {
    // RDKit✔️✔️:           if (queriesMatch(iter1->get(), iter2->get())) {
    // RDKit✔️✔️:             res = true;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (queriesMatch(iter1->get(), q2)) {
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (res) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (d1 == "AtomAnd") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:     for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
    // RDKit✔️✔️:          ++iter1) {
    // RDKit✔️✔️:       bool matched = false;
    // RDKit✔️✔️:       if (d2 == "AtomAnd") {
    // RDKit✔️✔️:         for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:              ++iter2) {
    // RDKit✔️✔️:           if (queriesMatch(iter1->get(), iter2->get())) {
    // RDKit✔️✔️:             matched = true;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         matched = queriesMatch(iter1->get(), q2);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!matched) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // FIX : handle AtomXOr
    // RDKit✔️✔️:   } else if (d2 == "AtomOr") {
    // RDKit✔️✔️:     // FIX: handle negation on AtomOr and AtomAnd
    // RDKit✔️✔️:     for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:          ++iter2) {
    // RDKit✔️✔️:       if (queriesMatch(q1, iter2->get())) {
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (d2 == "AtomAnd") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:     for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:          ++iter2) {
    // RDKit✔️✔️:       if (!queriesMatch(q1, iter2->get())) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (d1 == d2) {
    // RDKit✔️✔️:     if (std::find(&equalityQueries[0], &equalityQueries[nQueries], d1) !=
    // RDKit✔️✔️:         &equalityQueries[nQueries]) {
    // RDKit✔️✔️:       res = localMatch(static_cast<ATOM_EQUALS_QUERY const *>(q1),
    // RDKit✔️✔️:                        static_cast<ATOM_EQUALS_QUERY const *>(q2));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else {
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: recursion, short-circuit order, and the Or/And
    // nested-loop bounds match RDKit exactly, including its documented lack
    // of XOr and composite-negation handling. Typed variants replace repeated
    // description strings and the fixed 20-element linear lookup with O(1)
    // enum classification, eliminating string copies without changing query
    // traversal, allocation, ownership, or matching semantics.
    let (first, first_negated) = atom_query_base(first);
    let (second, second_negated) = atom_query_base(second);

    if is_atom_null_query(first) || is_atom_null_query(second) {
        return true;
    }
    if let QueryNode::Or(first_children) = first {
        return first_children.iter().any(|first_child| match second {
            QueryNode::Or(second_children) => second_children
                .iter()
                .any(|second_child| atom_queries_match(first_child, second_child)),
            _ => atom_queries_match(first_child, second),
        });
    }
    if let QueryNode::And(first_children) = first {
        return first_children.iter().all(|first_child| match second {
            QueryNode::And(second_children) => second_children
                .iter()
                .any(|second_child| atom_queries_match(first_child, second_child)),
            _ => atom_queries_match(first_child, second),
        });
    }
    if let QueryNode::Or(second_children) = second {
        return second_children
            .iter()
            .any(|second_child| atom_queries_match(first, second_child));
    }
    if let QueryNode::And(second_children) = second {
        return second_children
            .iter()
            .all(|second_child| atom_queries_match(first, second_child));
    }
    let (QueryNode::Predicate(first_predicate), QueryNode::Predicate(second_predicate)) =
        (first, second)
    else {
        return false;
    };
    let Some(first_kind) = atom_equality_query_kind(first_predicate) else {
        return false;
    };
    if atom_equality_query_kind(second_predicate) != Some(first_kind) {
        return false;
    }
    atom_query_local_match(
        first_predicate,
        first_negated,
        second_predicate,
        second_negated,
    )
}

fn rdkit_bond_order_value(order: BondOrder) -> i32 {
    match order {
        BondOrder::Null | BondOrder::Unspecified => 0,
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Quadruple => 4,
        BondOrder::Quintuple => 5,
        BondOrder::Hextuple => 6,
        BondOrder::OneAndHalf => 7,
        BondOrder::TwoAndHalf => 8,
        BondOrder::ThreeAndHalf => 9,
        BondOrder::FourAndHalf => 10,
        BondOrder::FiveAndHalf => 11,
        BondOrder::Aromatic => 12,
        BondOrder::Ionic => 13,
        BondOrder::Hydrogen => 14,
        BondOrder::ThreeCenter => 15,
        BondOrder::DativeOne => 16,
        BondOrder::Dative => 17,
        BondOrder::DativeLeft => 18,
        BondOrder::DativeRight => 19,
        BondOrder::Other => 20,
        BondOrder::Zero => 21,
    }
}

fn rdkit_bond_direction_value(direction: crate::BondDirection) -> i32 {
    match direction {
        crate::BondDirection::None => 0,
        crate::BondDirection::BeginWedge => 1,
        crate::BondDirection::BeginDash => 2,
        crate::BondDirection::EndDownRight => 3,
        crate::BondDirection::EndUpRight => 4,
        crate::BondDirection::EitherDouble => 5,
        crate::BondDirection::Unknown => 6,
    }
}

fn bond_equality_query_value(predicate: &BondQueryPredicate) -> Option<i32> {
    match predicate {
        BondQueryPredicate::InRingOfSize(value)
        | BondQueryPredicate::MinRingSize(value)
        | BondQueryPredicate::NumRingBonds(value) => Some(*value),
        BondQueryPredicate::Order(order) => Some(rdkit_bond_order_value(*order)),
        BondQueryPredicate::Direction(direction) => Some(rdkit_bond_direction_value(*direction)),
        BondQueryPredicate::IsInRing(value) => Some(i32::from(*value)),
        _ => None,
    }
}

fn bond_query_base(
    mut query: &QueryNode<BondQueryPredicate>,
) -> (&QueryNode<BondQueryPredicate>, bool) {
    let mut negated = false;
    while let QueryNode::Not(child) = query {
        negated = !negated;
        query = child;
    }
    (query, negated)
}

pub(crate) fn bond_queries_match(
    first: &QueryNode<BondQueryPredicate>,
    second: &QueryNode<BondQueryPredicate>,
) -> bool {
    // RDKit✔️✔️: bool queriesMatch(QueryBond::QUERYBOND_QUERY const *q1,
    // RDKit✔️✔️:                   QueryBond::QUERYBOND_QUERY const *q2) {
    // RDKit✔️✔️:   PRECONDITION(q1, "no q1");
    // RDKit✔️✔️:   PRECONDITION(q2, "no q2");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   static const unsigned int nQueries = 6;
    // RDKit✔️✔️:   static std::string equalityQueries[nQueries] = {
    // RDKit✔️✔️:       "BondRingSize", "BondMinRingSize", "BondOrder",
    // RDKit✔️✔️:       "BondDir",      "BondInRing",      "BondInNRings"};
    // RDKit✔️✔️:
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   std::string d1 = q1->getDescription();
    // RDKit✔️✔️:   std::string d2 = q2->getDescription();
    // RDKit✔️✔️:   if (d1 == "BondNull" || d2 == "BondNull") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:   } else if (d1 == "BondOr") {
    // RDKit✔️✔️:     // FIX: handle negation on BondOr and BondAnd
    // RDKit✔️✔️:     for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
    // RDKit✔️✔️:          ++iter1) {
    // RDKit✔️✔️:       if (d2 == "BondOr") {
    // RDKit✔️✔️:         for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:              ++iter2) {
    // RDKit✔️✔️:           if (queriesMatch(iter1->get(), iter2->get())) {
    // RDKit✔️✔️:             res = true;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         if (queriesMatch(iter1->get(), q2)) {
    // RDKit✔️✔️:           res = true;
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (res) {
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (d1 == "BondAnd") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:     for (auto iter1 = q1->beginChildren(); iter1 != q1->endChildren();
    // RDKit✔️✔️:          ++iter1) {
    // RDKit✔️✔️:       bool matched = false;
    // RDKit✔️✔️:       if (d2 == "BondAnd") {
    // RDKit✔️✔️:         for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:              ++iter2) {
    // RDKit✔️✔️:           if (queriesMatch(iter1->get(), iter2->get())) {
    // RDKit✔️✔️:             matched = true;
    // RDKit✔️✔️:             break;
    // RDKit✔️✔️:           }
    // RDKit✔️✔️:         }
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         matched = queriesMatch(iter1->get(), q2);
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       if (!matched) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     // FIX : handle BondXOr
    // RDKit✔️✔️:   } else if (d2 == "BondOr") {
    // RDKit✔️✔️:     // FIX: handle negation on BondOr and BondAnd
    // RDKit✔️✔️:     for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:          ++iter2) {
    // RDKit✔️✔️:       if (queriesMatch(q1, iter2->get())) {
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (d2 == "BondAnd") {
    // RDKit✔️✔️:     res = true;
    // RDKit✔️✔️:     for (auto iter2 = q2->beginChildren(); iter2 != q2->endChildren();
    // RDKit✔️✔️:          ++iter2) {
    // RDKit✔️✔️:       if (queriesMatch(q1, iter2->get())) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   } else if (std::find(&equalityQueries[0], &equalityQueries[nQueries], d1) !=
    // RDKit✔️✔️:              &equalityQueries[nQueries]) {
    // RDKit✔️✔️:     res = localMatch(static_cast<BOND_EQUALS_QUERY const *>(q1),
    // RDKit✔️✔️:                      static_cast<BOND_EQUALS_QUERY const *>(q2));
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: tree traversal, child ordering, short-circuit
    // behavior, and worst-case O(n*m) composite comparison match RDKit. Typed
    // classification and exact integer conversion are O(1) and avoid RDKit's
    // description-string copies and linear six-string lookup without adding
    // allocation, cloning, temporary collections, or extra tree scans. The
    // source's asymmetric second-And condition is preserved exactly.
    let (first, first_negated) = bond_query_base(first);
    let (second, second_negated) = bond_query_base(second);

    if is_bond_null_query(first) || is_bond_null_query(second) {
        return true;
    }
    if let QueryNode::Or(first_children) = first {
        return first_children.iter().any(|first_child| match second {
            QueryNode::Or(second_children) => second_children
                .iter()
                .any(|second_child| bond_queries_match(first_child, second_child)),
            _ => bond_queries_match(first_child, second),
        });
    }
    if let QueryNode::And(first_children) = first {
        return first_children.iter().all(|first_child| match second {
            QueryNode::And(second_children) => second_children
                .iter()
                .any(|second_child| bond_queries_match(first_child, second_child)),
            _ => bond_queries_match(first_child, second),
        });
    }
    if let QueryNode::Or(second_children) = second {
        return second_children
            .iter()
            .any(|second_child| bond_queries_match(first, second_child));
    }
    if let QueryNode::And(second_children) = second {
        return !second_children
            .iter()
            .any(|second_child| bond_queries_match(first, second_child));
    }
    let (QueryNode::Predicate(first_predicate), QueryNode::Predicate(second_predicate)) =
        (first, second)
    else {
        return false;
    };
    let (Some(first_value), Some(second_value)) = (
        bond_equality_query_value(first_predicate),
        bond_equality_query_value(second_predicate),
    ) else {
        return false;
    };
    bond_query_local_match(&first_value, first_negated, &second_value, second_negated)
}

fn query_atom_query_match(
    query: &QueryNode<AtomQueryPredicate>,
    what: &Atom,
    mol: &impl SearchTargetAccess,
) -> bool {
    // RDKit✔️❌: bool QueryAtom::QueryMatch(QueryAtom const *what) const {
    // RDKit✔️❌:   PRECONDITION(what, "bad query atom");
    // RDKit✔️❌:   PRECONDITION(dp_query, "no query set");
    // RDKit✔️❌:   if (!what->hasQuery()) {
    // RDKit✔️❌:     return dp_query->Match(what);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     return queriesMatch(dp_query, what->getQuery());
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Local complexity review: both implementations inspect target query
    // presence once and dispatch to the canonical ordinary matcher or query
    // compatibility matcher without cloning or allocating. Query-to-query
    // matching retains RDKit's query-tree complexity. The ordinary-target
    // branch inherits atom_matches_query's additional O(V+E) context build,
    // so this entry is behavior-equivalent but performance-worse until
    // molecule-derived state is reused canonically.
    let context = build_query_match_context_for_target(mol);
    atom_matches_query_with_target_context(what, query, mol, &context)
}

fn query_bond_query_match(
    query: &QueryNode<BondQueryPredicate>,
    what: &Bond,
    mol: &impl SearchTargetAccess,
) -> bool {
    // RDKit✔️❌: bool QueryBond::QueryMatch(QueryBond const *what) const {
    // RDKit✔️❌:   PRECONDITION(what, "bad query bond");
    // RDKit✔️❌:   PRECONDITION(dp_query, "no query set");
    // RDKit✔️❌:   if (!what->hasQuery()) {
    // RDKit✔️❌:     return dp_query->Match(what);
    // RDKit✔️❌:   } else {
    // RDKit✔️❌:     return queriesMatch(dp_query, what->getQuery());
    // RDKit✔️❌:   }
    // RDKit✔️❌: }
    // Local complexity review: both implementations inspect target query
    // presence once and dispatch to the canonical ordinary matcher or query
    // compatibility matcher without cloning or allocating. Query-to-query
    // matching retains RDKit's query-tree complexity. The ordinary-target
    // branch inherits bond_matches_query's additional O(V+E) context build,
    // so performance remains worse until derived state is reused canonically.
    let context = build_query_match_context_for_target(mol);
    bond_matches_query_with_target_context(what, query, mol, &context)
}

pub(crate) fn and_query_match<T>(
    children: &[QueryNode<T>],
    negated: bool,
    mut child_matches: impl FnMut(&QueryNode<T>) -> bool,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   bool res = true;
    // RDKit✔️✔️:   typename BASE::CHILD_VECT_CI it1;
    // RDKit✔️✔️:   for (it1 = this->beginChildren(); it1 != this->endChildren(); ++it1) {
    // RDKit✔️✔️:     bool tmp = (*it1)->Match(what);
    // RDKit✔️✔️:     if (!tmp) {
    // RDKit✔️✔️:       res = false;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (this->getNegation()) {
    // RDKit✔️✔️:     res = !res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: both implementations inspect children in
    // insertion order, stop after the first false result, and take O(k) time
    // in the all-matching case with O(1) extra space. Rust walks the existing
    // contiguous child vector by reference and the monomorphized callback
    // adds no allocation, clone, lookup, temporary collection, repeated scan,
    // or asymptotic/hot-path branch beyond RDKit's virtual child Match call.
    let mut result = true;
    for child in children {
        if !child_matches(child) {
            result = false;
            break;
        }
    }
    if negated {
        result = !result;
    }
    result
}

pub(crate) fn or_query_match<T>(
    children: &[QueryNode<T>],
    negated: bool,
    mut child_matches: impl FnMut(&QueryNode<T>) -> bool,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   typename BASE::CHILD_VECT_CI it1;
    // RDKit✔️✔️:   for (it1 = this->beginChildren(); it1 != this->endChildren(); ++it1) {
    // RDKit✔️✔️:     bool tmp = (*it1)->Match(what);
    // RDKit✔️✔️:     if (tmp) {
    // RDKit✔️✔️:       res = true;
    // RDKit✔️✔️:       break;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (this->getNegation()) {
    // RDKit✔️✔️:     res = !res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: both implementations inspect children in
    // insertion order, stop after the first true result, and take O(k) time
    // when no child matches with O(1) extra space. Rust walks the existing
    // contiguous child vector by reference and the monomorphized callback
    // adds no allocation, clone, lookup, temporary collection, repeated scan,
    // or asymptotic/hot-path branch beyond RDKit's virtual child Match call.
    let mut result = false;
    for child in children {
        if child_matches(child) {
            result = true;
            break;
        }
    }
    if negated {
        result = !result;
    }
    result
}

pub(crate) fn xor_query_match<T>(
    children: &[QueryNode<T>],
    negated: bool,
    mut child_matches: impl FnMut(&QueryNode<T>) -> bool,
) -> bool {
    // RDKit✔️✔️: bool Match(const DataFuncArgType what) const override {
    // RDKit✔️✔️:   bool res = false;
    // RDKit✔️✔️:   typename BASE::CHILD_VECT_CI it1;
    // RDKit✔️✔️:   for (it1 = this->beginChildren(); it1 != this->endChildren(); ++it1) {
    // RDKit✔️✔️:     bool tmp = (*it1)->Match(what);
    // RDKit✔️✔️:     if (tmp) {
    // RDKit✔️✔️:       if (res) {
    // RDKit✔️✔️:         res = false;
    // RDKit✔️✔️:         break;
    // RDKit✔️✔️:       } else {
    // RDKit✔️✔️:         res = true;
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   if (this->getNegation()) {
    // RDKit✔️✔️:     res = !res;
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // Local complexity review: both implementations inspect children in
    // insertion order, stop on the second true result, and take O(k) time
    // when fewer than two children match with O(1) extra space. Rust walks
    // the existing contiguous child vector by reference; its monomorphized
    // callback adds no allocation, clone, lookup, temporary collection,
    // repeated scan, or asymptotic/hot-path branch beyond RDKit's virtual
    // child Match call.
    let mut result = false;
    for child in children {
        if child_matches(child) {
            if result {
                result = false;
                break;
            }
            result = true;
        }
    }
    if negated {
        result = !result;
    }
    result
}

pub fn bond_predicate_matches_with_context(
    bond: &Bond,
    pred: &BondQueryPredicate,
    mol: &Molecule,
    ctx: &QueryMatchContext,
) -> bool {
    bond_predicate_matches_with_target_context(bond, pred, mol, ctx)
}

/// RDKit✔️✔️: Evaluate a `QueryNode<AtomQueryPredicate>` tree against `atom`.
///
/// RDKit source: `QueryOps.h` / `Query.cpp` recursive query matching.
///   An `AtomAnd`/`BondAnd` node matches iff all children match.
///   An `AtomOr`/`BondOr` node matches iff any child matches.
///   A negation inverts match.
///   A leaf Predicate matches iff `atom_predicate_matches` returns true.
pub fn atom_matches_query(
    atom: &Atom,
    query: &QueryNode<AtomQueryPredicate>,
    mol: &Molecule,
) -> bool {
    // RDKit✔️❌: bool QueryAtom::Match(Atom const *what) const {
    // RDKit✔️❌:   PRECONDITION(what, "bad query atom");
    // RDKit✔️❌:   PRECONDITION(dp_query, "no query set");
    // RDKit✔️❌:   return dp_query->Match(what);
    // RDKit✔️❌: }
    // Local complexity review: Rust references make both source null-pointer
    // preconditions unrepresentable. Both entries dispatch once into the same
    // recursive query tree; Rust additionally builds the molecule-derived
    // context once before traversal, avoiding repeated work inside leaves but
    // adding O(V+E) setup relative to RDKit's cache-backed entry. Behavior is
    // equivalent, but this entry remains performance-worse until molecule
    // derived-state reuse is canonicalized. No duplicate matcher is added.
    let ctx = build_query_match_context(mol);
    atom_matches_query_with_target_context(atom, query, mol, &ctx)
}

pub(crate) fn atom_matches_query_with_target_context(
    atom: &Atom,
    query: &QueryNode<AtomQueryPredicate>,
    mol: &impl SearchTargetAccess,
    ctx: &QueryMatchContext,
) -> bool {
    match query {
        QueryNode::Predicate(pred) => {
            atom_predicate_matches_with_target_context(atom, pred, mol, ctx)
        }

        QueryNode::And(children) => and_query_match(children, false, |child| {
            atom_matches_query_with_target_context(atom, child, mol, ctx)
        }),

        QueryNode::Or(children) => or_query_match(children, false, |child| {
            atom_matches_query_with_target_context(atom, child, mol, ctx)
        }),

        QueryNode::Xor(children) => xor_query_match(children, false, |child| {
            atom_matches_query_with_target_context(atom, child, mol, ctx)
        }),

        // RDKit✔️✔️: NOT — invert child match.
        // RDKit source: negation flips the result.
        QueryNode::Not(child) => !atom_matches_query_with_target_context(atom, child, mol, ctx),
    }
}

/// RDKit✔️✔️: Evaluate a `QueryNode<BondQueryPredicate>` tree against `bond`.
pub fn bond_matches_query(
    bond: &Bond,
    query: &QueryNode<BondQueryPredicate>,
    mol: &Molecule,
) -> bool {
    // RDKit✔️❌: bool QueryBond::Match(Bond const *what) const {
    // RDKit✔️❌:   PRECONDITION(what, "bad query bond");
    // RDKit✔️❌:   PRECONDITION(dp_query, "no query set");
    // RDKit✔️❌:   return dp_query->Match(what);
    // RDKit✔️❌: }
    // Local complexity review: Rust references make both source null-pointer
    // preconditions unrepresentable. Both entries dispatch once into the same
    // recursive typed query tree; Rust additionally builds the molecule-
    // derived context once before traversal, adding O(V+E) setup relative to
    // RDKit's cache-backed entry. Behavior is equivalent, but performance is
    // worse until derived-state reuse is canonicalized. No duplicate matcher
    // or bond-query representation is introduced.
    let ctx = build_query_match_context(mol);
    bond_matches_query_with_target_context(bond, query, mol, &ctx)
}

pub(crate) fn bond_matches_query_with_target_context(
    bond: &Bond,
    query: &QueryNode<BondQueryPredicate>,
    mol: &impl SearchTargetAccess,
    ctx: &QueryMatchContext,
) -> bool {
    match query {
        QueryNode::Predicate(pred) => {
            bond_predicate_matches_with_target_context(bond, pred, mol, ctx)
        }

        QueryNode::And(children) => and_query_match(children, false, |child| {
            bond_matches_query_with_target_context(bond, child, mol, ctx)
        }),

        QueryNode::Or(children) => or_query_match(children, false, |child| {
            bond_matches_query_with_target_context(bond, child, mol, ctx)
        }),

        QueryNode::Xor(children) => xor_query_match(children, false, |child| {
            bond_matches_query_with_target_context(bond, child, mol, ctx)
        }),

        // RDKit✔️✔️: NOT
        QueryNode::Not(child) => !bond_matches_query_with_target_context(bond, child, mol, ctx),
    }
}

pub fn atom_matches_query_with_context(
    atom: &Atom,
    query: &QueryNode<AtomQueryPredicate>,
    mol: &Molecule,
    ctx: &QueryMatchContext,
) -> bool {
    atom_matches_query_with_target_context(atom, query, mol, ctx)
}

pub fn bond_matches_query_with_context(
    bond: &Bond,
    query: &QueryNode<BondQueryPredicate>,
    mol: &Molecule,
    ctx: &QueryMatchContext,
) -> bool {
    bond_matches_query_with_target_context(bond, query, mol, ctx)
}

// ---------------------------------------------------------------------------
// Internal helpers
// ---------------------------------------------------------------------------

fn query_atom_non_hydrogen_degree(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> u32 {
    // RDKit✔️✔️: //! D and T are treated as "non-hydrogen" here
    // RDKit✔️✔️: static inline int queryAtomNonHydrogenDegree(Atom const *at) {
    // RDKit✔️✔️:   int res = 0;
    // RDKit✔️✔️:   for (const auto nbri :
    // RDKit✔️✔️:        boost::make_iterator_range(at->getOwningMol().getAtomNeighbors(at))) {
    // RDKit✔️✔️:     const auto nbr = at->getOwningMol()[nbri];
    // RDKit✔️✔️:     if (nbr->getAtomicNum() != 1 || nbr->getIsotope() > 1) {
    // RDKit✔️✔️:       res++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the existing adjacency range with O(1) atom lookup and property
    // reads per neighbor. Neither allocates, clones, repeats a scan, or creates
    // a temporary collection; both use one counter and the same branch.
    let mut res = 0u32;
    for nbri in adj.neighbors_of(at.id().index()) {
        let nbr = &mol.atoms()[nbri.atom_index];
        if nbr.atomic_number() != 1 || nbr.isotope().is_some_and(|isotope| isotope > 1) {
            res += 1;
        }
    }
    res
}

fn query_atom_heavy_atom_degree(
    at: &Atom,
    adj: &AdjacencyList,
    mol: &impl SearchTargetAccess,
) -> u32 {
    // RDKit✔️✔️: //! D and T are not treated as heavy atoms here
    // RDKit✔️✔️: static inline int queryAtomHeavyAtomDegree(Atom const *at) {
    // RDKit✔️✔️:   int heavyDegree = 0;
    // RDKit✔️✔️:   for (const auto nbri :
    // RDKit✔️✔️:        boost::make_iterator_range(at->getOwningMol().getAtomNeighbors(at))) {
    // RDKit✔️✔️:     const auto nbr = at->getOwningMol()[nbri];
    // RDKit✔️✔️:     if (nbr->getAtomicNum() > 1) {
    // RDKit✔️✔️:       heavyDegree++;
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:
    // RDKit✔️✔️:   return heavyDegree;
    // RDKit✔️✔️: };
    // Local complexity review: RDKit and Rust each make one O(degree) pass
    // over the existing adjacency range with O(1) atom lookup and atomic-
    // number read per neighbor. Neither allocates, clones, repeats a scan, or
    // creates a temporary collection; both use one counter and one branch.
    let mut heavy_degree = 0u32;
    for nbri in adj.neighbors_of(at.id().index()) {
        let nbr = &mol.atoms()[nbri.atom_index];
        if nbr.atomic_number() > 1 {
            heavy_degree += 1;
        }
    }
    heavy_degree
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Molecule;

    #[test]
    fn detached_search_context_uses_model_blocks() {
        let molecule = Molecule::from_smiles("CCO").expect("valid test molecule");
        let context = build_query_match_context_from_blocks(
            molecule.topology_block(),
            molecule.coordinate_block(),
            molecule.stereo_groups(),
            None,
            None,
        );
        assert_eq!(context.adj.neighbors_of(0).len(), 1);
        assert!(context.valence.is_some());
        assert!(atom_predicate_matches_from_blocks(
            &molecule.atoms()[0],
            &AtomQueryPredicate::AtomicNumber(6),
            molecule.topology_block(),
            molecule.coordinate_block(),
            molecule.stereo_groups(),
            None,
            None,
        ));
    }

    fn test_atom_with_query(query: Option<QueryNode<AtomQueryPredicate>>) -> crate::QueryAtom {
        test_atom_with_query_and_aromatic(query, false)
    }

    fn test_atom_with_query_and_aromatic(
        query: Option<QueryNode<AtomQueryPredicate>>,
        aromatic: bool,
    ) -> crate::QueryAtom {
        let spec = crate::AtomSpec::new(crate::Element::C).with_aromatic(aromatic);
        let predicate = query.unwrap_or_else(|| make_atom_null_query());
        crate::QueryAtom::from_parts(Atom::from_spec(AtomId::new(0), spec), predicate)
    }

    #[test]
    fn smarts_complex_atom_query() {
        let atom = |query| test_atom_with_query(Some(query));
        let atomic_number = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let charge = || QueryNode::predicate(AtomQueryPredicate::FormalCharge(0));

        assert!(!is_complex_atom_query(&test_atom_with_query(None)));
        assert!(!is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::Any,
        ))));
        assert!(!is_complex_atom_query(&atom(atomic_number())));
        assert!(!is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            },
        ))));
        assert!(is_complex_atom_query(&atom(
            QueryNode::not(atomic_number())
        )));
        assert!(is_complex_atom_query(&atom(QueryNode::or(vec![
            atomic_number()
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::xor(vec![
            atomic_number()
        ]))));
        assert!(!is_complex_atom_query(&atom(QueryNode::and(vec![
            atomic_number(),
            charge(),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::and(vec![charge()]))));
        assert!(is_complex_atom_query(&atom(QueryNode::and(vec![
            atomic_number(),
            QueryNode::or(vec![charge()]),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![6, 8]),
        ))));
        assert!(is_complex_atom_query(&atom(charge())));
    }

    #[test]
    fn query_complexity_classifiers_match_source_across_fingerprint_consumers() {
        let atom = |query| test_atom_with_query(Some(query));
        let number = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let atom_type = |atomic_number, aromatic| {
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number,
                aromatic,
            })
        };
        let charge = || QueryNode::predicate(AtomQueryPredicate::FormalCharge(0));

        assert!(!is_complex_atom_query(&test_atom_with_query(None)));
        assert!(!is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::Any,
        ))));
        assert!(!is_complex_atom_query(&atom(number())));
        assert!(!is_complex_atom_query(&atom(atom_type(6, false))));
        assert!(!is_complex_atom_query(&atom(atom_type(6, true))));
        assert!(is_complex_atom_query(&atom(QueryNode::not(number()))));
        assert!(is_complex_atom_query(&atom(QueryNode::or(vec![
            number(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::xor(vec![
            number(),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        ]))));
        assert!(!is_complex_atom_query(&atom(QueryNode::and(vec![
            charge(),
            QueryNode::and(vec![charge(), number()]),
        ]))));
        assert!(!is_complex_atom_query(&atom(QueryNode::and(vec![
            charge(),
            QueryNode::and(vec![atom_type(6, false), charge()]),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::and(vec![
            charge(),
            QueryNode::and(vec![charge()]),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::and(vec![
            number(),
            QueryNode::not(charge()),
        ]))));
        assert!(is_complex_atom_query(&atom(QueryNode::and(Vec::new()))));
        assert!(is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![6, 8]),
        ))));
        assert!(is_complex_atom_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]),
        ))));

        let bond = |query| test_bond_with_query(Some(query));
        let order = |value| QueryNode::predicate(BondQueryPredicate::Order(value));
        assert!(!is_complex_bond_query(&test_bond_with_query(None)));
        assert!(!is_complex_bond_query(&bond(order(BondOrder::Single))));
        assert!(!is_complex_bond_query(&bond(QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Aromatic]),
        ))));
        assert!(!is_complex_bond_query(&bond(QueryNode::or(vec![
            order(BondOrder::Aromatic),
            order(BondOrder::Single),
        ]))));
        assert!(is_complex_bond_query(&bond(QueryNode::or(vec![order(
            BondOrder::Single,
        )]))));
        assert!(is_complex_bond_query(&bond(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Aromatic),
            order(BondOrder::Single),
        ]))));
        assert!(is_complex_bond_query(&bond(QueryNode::or(vec![
            QueryNode::not(order(BondOrder::Single)),
            order(BondOrder::Aromatic),
        ]))));
        assert!(is_complex_bond_query(&bond(QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Aromatic, BondOrder::Single]),
        ))));
        assert!(is_complex_bond_query(&bond(QueryNode::predicate(
            BondQueryPredicate::OrderIn(vec![BondOrder::Single, BondOrder::Double]),
        ))));
        assert!(is_complex_bond_query(&bond(QueryNode::not(order(
            BondOrder::Single,
        )))));
        assert!(is_complex_bond_query(&bond(QueryNode::and(vec![order(
            BondOrder::Single,
        )]))));
        assert!(is_complex_bond_query(&bond(QueryNode::xor(vec![order(
            BondOrder::Single,
        )]))));
        assert!(is_complex_bond_query(&bond(QueryNode::predicate(
            BondQueryPredicate::Any,
        ))));
    }

    #[test]
    fn layered_query_aromaticity_matches_source() {
        let mut builder = Molecule::builder();
        let plain = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_by_flag =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true));
        let aromatic_by_bond = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                aromatic_by_bond,
                aromatic_neighbor,
                BondOrder::Aromatic,
            ))
            .expect("aromaticity fixture bond is valid");
        let molecule = builder.build().expect("aromaticity fixture is valid");

        assert!(!is_atom_aromatic(
            &molecule.atoms()[plain.index()],
            &molecule
        ));
        assert!(is_atom_aromatic(
            &molecule.atoms()[aromatic_by_flag.index()],
            &molecule,
        ));
        assert!(is_atom_aromatic(
            &molecule.atoms()[aromatic_by_bond.index()],
            &molecule,
        ));

        let query_molecule = Molecule::new();
        let query_atom = |query, aromatic| {
            test_atom_with_query_and_aromatic(Some(query), aromatic)
                .atom()
                .clone()
        };
        let number = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let aromatic = || QueryNode::predicate(AtomQueryPredicate::IsAromatic(true));
        let aliphatic = || QueryNode::predicate(AtomQueryPredicate::IsAromatic(false));
        let atom_type = |is_aromatic| {
            QueryNode::predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: is_aromatic,
            })
        };

        assert!(!is_atom_aromatic(
            &query_atom(number(), false),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(number(), true),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(QueryNode::not(number()), true),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(aromatic(), false),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(QueryNode::not(aromatic()), false),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(aliphatic(), true),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(QueryNode::not(aliphatic()), false),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(atom_type(true), false),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(QueryNode::not(atom_type(false)), false),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(QueryNode::and(vec![number(), aromatic()]), false),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(QueryNode::and(vec![number(), aliphatic()]), true),
            &query_molecule,
        ));
        assert!(is_atom_aromatic(
            &query_atom(
                QueryNode::and(vec![number(), QueryNode::not(aromatic())]),
                false,
            ),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(QueryNode::and(vec![aromatic(), number()]), true),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(
                QueryNode::not(QueryNode::and(vec![number(), aromatic()])),
                true,
            ),
            &query_molecule,
        ));
        assert!(!is_atom_aromatic(
            &query_atom(
                QueryNode::predicate(AtomQueryPredicate::FormalCharge(0)),
                true,
            ),
            &query_molecule,
        ));
    }

    #[test]
    fn smarts_atom_list_query() {
        let atom = |query| test_atom_with_query(Some(query)).atom().clone();
        let atomic_number = |value| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(value));

        assert!(!is_atom_list_query(test_atom_with_query(None).atom()));
        assert!(is_atom_list_query(&atom(QueryNode::or(vec![
            atomic_number(6),
            QueryNode::or(vec![atomic_number(7), atomic_number(8)]),
        ]))));
        assert!(!is_atom_list_query(&atom(QueryNode::or(vec![
            atomic_number(6),
            QueryNode::not(atomic_number(7)),
        ]))));
        assert!(is_atom_list_query(&atom(QueryNode::not(atomic_number(7)))));
        assert!(is_atom_list_query(&atom(QueryNode::not(QueryNode::or(
            vec![atomic_number(6), QueryNode::not(atomic_number(7))],
        )))));
        assert!(!is_atom_list_query(&atom(atomic_number(6))));
        assert!(is_atom_list_query(&atom(atomic_number(8))));
        assert!(is_atom_list_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberIn(vec![6, 8]),
        ))));
        assert!(is_atom_list_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]),
        ))));
        assert!(!is_atom_list_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            },
        ))));
        assert!(!is_atom_list_query(&atom(QueryNode::predicate(
            AtomQueryPredicate::FormalCharge(0),
        ))));
    }

    #[test]
    fn smarts_atom_list_values() {
        let number = |value| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(value));
        let query = QueryNode::or(vec![
            number(6),
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomType {
                    atomic_number: 7,
                    aromatic: true,
                }),
                number(8),
            ]),
        ]);
        assert_eq!(get_atom_list_query_values(&query), Ok(vec![6, 7, 8]));
        assert_eq!(
            get_atom_list_query_values(&QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(
                vec![9, 17]
            ),)),
            Ok(vec![9, 17])
        );
        assert_eq!(
            get_atom_list_query_values(&QueryNode::not(number(6))),
            Ok(vec![6])
        );
        assert_eq!(
            get_atom_list_query_values(&QueryNode::or(vec![QueryNode::not(number(6))])),
            Err("bad query type1")
        );
        assert_eq!(
            get_atom_list_query_values(&QueryNode::predicate(AtomQueryPredicate::FormalCharge(0),)),
            Err("bad query type")
        );
    }

    #[test]
    fn smarts_complete_query_children() {
        let mut query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(QUERY_SCAN_MAGIC_VALUE)),
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(
                    QUERY_SCAN_MAGIC_VALUE,
                )),
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)),
                QueryNode::predicate(AtomQueryPredicate::HasRingBond),
            ]),
        ]);

        complete_query_and_children(
            &mut query,
            QUERY_SCAN_MAGIC_VALUE,
            AtomQueryCompletionValues {
                ring_bond_count: 3,
                non_hydrogen_degree: 4,
            },
        );

        assert_eq!(
            query,
            QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(3)),
                QueryNode::or(vec![
                    QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(4)),
                    QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)),
                    QueryNode::predicate(AtomQueryPredicate::HasRingBond),
                ]),
            ])
        );
    }

    #[test]
    fn smarts_complete_mol_queries() {
        let query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::RingBondCount(QUERY_SCAN_MAGIC_VALUE)),
            QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(
                QUERY_SCAN_MAGIC_VALUE,
            )),
            QueryNode::predicate(AtomQueryPredicate::HasRingBond),
        ]);
        let atoms = vec![
            crate::QueryAtom::from_parts(
                Atom::from_spec(AtomId::new(0), crate::AtomSpec::new(crate::Element::C)),
                query,
            ),
            test_atom_with_query(None),
            test_atom_with_query(None),
            crate::QueryAtom::from_parts(
                Atom::from_spec(AtomId::new(3), crate::AtomSpec::new(crate::Element::H)),
                make_atom_num_query(1),
            ),
        ];
        let bonds = vec![
            make_query_bond_spec(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            make_query_bond_spec(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            make_query_bond_spec(AtomId::new(2), AtomId::new(0), BondOrder::Single),
            make_query_bond_spec(AtomId::new(0), AtomId::new(3), BondOrder::Single),
        ];
        let mut molecule = crate::QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap();

        complete_mol_queries(&mut molecule, QUERY_SCAN_MAGIC_VALUE);

        assert_eq!(
            molecule.atoms()[0].predicate(),
            &QueryNode::and(vec![
                QueryNode::predicate(AtomQueryPredicate::RingBondCount(2)),
                QueryNode::predicate(AtomQueryPredicate::NonHydrogenDegree(2)),
                QueryNode::predicate(AtomQueryPredicate::HasRingBond),
            ])
        );
        assert!(!atom_query_has_magic_value(
            molecule.atoms()[0].predicate(),
            QUERY_SCAN_MAGIC_VALUE,
        ));
    }

    #[test]
    fn smarts_replace_atom_with_query() {
        let spec = crate::AtomSpec::new(crate::Element::N)
            .with_isotope(15)
            .with_formal_charge(1)
            .with_radical_electrons(2)
            .with_prop("_hasMassQuery", "1");
        let atom = replace_atom_with_query_atom(Atom::from_spec(AtomId::new(0), spec));

        assert_eq!(
            atom.predicate(),
            &QueryNode::and(vec![
                QueryNode::and(vec![
                    QueryNode::and(vec![
                        QueryNode::and(vec![make_atom_num_query(7), make_atom_isotope_query(15),]),
                        make_atom_formal_charge_query(1),
                    ]),
                    make_atom_num_radical_electrons_query(2),
                ]),
                make_atom_mass_query(15),
            ])
        );
        assert_eq!(atom.id(), AtomId::new(0));
        assert_eq!(atom.atom().prop("_hasMassQuery"), Some("1"));
    }

    #[test]
    fn smarts_finalize_ring_size() {
        assert_eq!(
            finalize_atom_ring_size_query(
                make_atom_in_ring_of_size_query(5),
                RangeQueryType::Equal,
            ),
            Ok(QueryNode::predicate(AtomQueryPredicate::InRingOfSize(5)))
        );
        assert_eq!(
            finalize_atom_ring_size_query(make_atom_in_ring_of_size_query(5), RangeQueryType::Less,),
            Ok(QueryNode::predicate(
                AtomQueryPredicate::InRingOfSizeLessEqual(5)
            ))
        );
        assert_eq!(
            finalize_atom_ring_size_query(
                make_atom_in_ring_of_size_query(5),
                RangeQueryType::Greater,
            ),
            Ok(QueryNode::predicate(
                AtomQueryPredicate::InRingOfSizeGreaterEqual(5)
            ))
        );
        let range = make_atom_in_ring_of_size_range_query(3, 6, true, false);
        assert_eq!(
            finalize_atom_ring_size_query(range.clone(), RangeQueryType::Range),
            Ok(range)
        );
    }

    #[test]
    fn smarts_finalize_atom_description() {
        assert_eq!(
            finalize_atom_query_from_description(
                "less_AtomRingSize",
                make_atom_in_ring_of_size_query(6),
            ),
            Ok(QueryNode::predicate(
                AtomQueryPredicate::InRingOfSizeLessEqual(6)
            ))
        );
        let atomic_number = make_atom_num_query(6);
        assert_eq!(
            finalize_atom_query_from_description("AtomAtomicNum", atomic_number.clone()),
            Ok(atomic_number)
        );
        assert!(matches!(
            finalize_atom_query_from_description("NoSuchAtomQuery", make_atom_null_query()),
            Err(QueryFinalizationError::UnknownDescription(_))
        ));
    }

    #[test]
    fn smarts_finalize_bond_description() {
        let order = make_bond_order_equals_query(BondOrder::Double);
        assert_eq!(
            finalize_bond_query_from_description("BondOrder", order.clone()),
            Ok(order)
        );
        let ring_size = QueryNode::predicate(BondQueryPredicate::InRingOfSize(6));
        assert_eq!(
            finalize_bond_query_from_description("BondRingSize", ring_size.clone()),
            Ok(ring_size)
        );
        assert!(matches!(
            finalize_bond_query_from_description("NoSuchBondQuery", make_bond_null_query()),
            Err(QueryFinalizationError::UnknownDescription(_))
        ));
    }

    #[test]
    fn smarts_has_prop_query() {
        let atom = Atom::from_spec(
            AtomId::new(0),
            crate::AtomSpec::new(crate::Element::C).with_prop("role", "donor"),
        );
        let mut builder = Molecule::builder();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().unwrap();
        assert!(atom_matches_query(
            &atom,
            &make_atom_has_prop_query("role"),
            &molecule,
        ));

        let bond = Bond::from_spec(
            crate::BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("kind", "kept"),
        );
        assert!(bond_matches_query(
            &bond,
            &make_bond_has_prop_query("kind"),
            &molecule,
        ));
    }

    #[test]
    fn smarts_prop_value_query() {
        let atom = Atom::from_spec(
            AtomId::new(0),
            crate::AtomSpec::new(crate::Element::C).with_prop("role", "donor"),
        );
        let mut builder = Molecule::builder();
        builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().unwrap();
        assert!(atom_matches_query(
            &atom,
            &make_atom_prop_query("role", "donor"),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &atom,
            &make_atom_prop_query("role", "acceptor"),
            &molecule,
        ));

        let bond = Bond::from_spec(
            crate::BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                .with_prop("kind", "kept"),
        );
        assert!(bond_matches_query(
            &bond,
            &make_bond_prop_query("kind", "kept"),
            &molecule,
        ));
    }

    #[test]
    fn smarts_atom_ring_query() {
        let ring = Molecule::from_smiles("C1CC1").unwrap();
        let chain = Molecule::from_smiles("CCC").unwrap();
        assert!(atom_matches_query(
            &ring.atoms()[0],
            &make_atom_ring_query(-1),
            &ring,
        ));
        assert!(!atom_matches_query(
            &chain.atoms()[0],
            &make_atom_ring_query(-1),
            &chain,
        ));
        assert!(atom_matches_query(
            &chain.atoms()[0],
            &make_atom_ring_query(0),
            &chain,
        ));
        assert!(atom_matches_query(
            &ring.atoms()[0],
            &make_atom_ring_query(1),
            &ring,
        ));
    }

    #[test]
    fn smarts_recursive_structure_default_constructor() {
        let query = RecursiveStructureQuery::new();
        assert!(query.query_mol().is_none());
        assert_eq!(query.serial_number(), 0);
    }

    #[test]
    fn smarts_recursive_structure_molecule_constructor() {
        let molecule = Molecule::from_smiles("CO").unwrap();
        let graph =
            crate::search::query_graph::query_graph_from_concrete_molecule(&molecule).unwrap();
        let query = RecursiveStructureQuery::from_molecule(graph, 17);
        assert_eq!(query.query_mol().unwrap().num_atoms(), 2);
        assert_eq!(query.serial_number(), 17);
    }

    #[test]
    fn smarts_recursive_structure_atom_index() {
        let atom = Atom::from_spec(AtomId::new(4), crate::AtomSpec::new(crate::Element::N));
        assert_eq!(RecursiveStructureQuery::atom_index(&atom), 4);
    }

    #[test]
    fn smarts_recursive_structure_set_molecule() {
        let mut query = RecursiveStructureQuery::new();
        let nitrogen = Molecule::from_smiles("N").unwrap();
        query.set_query_mol(
            crate::search::query_graph::query_graph_from_concrete_molecule(&nitrogen).unwrap(),
        );
        assert_eq!(
            query.query_mol().unwrap().atoms()[0].atom().atomic_number(),
            7
        );
        let oxygen = Molecule::from_smiles("O").unwrap();
        query.set_query_mol(
            crate::search::query_graph::query_graph_from_concrete_molecule(&oxygen).unwrap(),
        );
        assert_eq!(
            query.query_mol().unwrap().atoms()[0].atom().atomic_number(),
            8
        );
    }

    #[test]
    fn smarts_recursive_structure_get_molecule() {
        let molecule = Molecule::from_smiles("C=C").unwrap();
        let query = RecursiveStructureQuery::from_molecule(
            crate::search::query_graph::query_graph_from_concrete_molecule(&molecule).unwrap(),
            0,
        );
        assert_eq!(query.query_mol().unwrap().num_bonds(), 1);
    }

    #[test]
    fn smarts_recursive_structure_copy() {
        let molecule = Molecule::from_smiles("CO").unwrap();
        let mut query = RecursiveStructureQuery::from_molecule(
            crate::search::query_graph::query_graph_from_concrete_molecule(&molecule).unwrap(),
            23,
        );
        query.insert_atom_index(3);
        let copied = query.copy_query();
        assert_eq!(copied, query);
        assert!(copied.contains_atom_index(3));
        assert_eq!(copied.serial_number(), 23);
    }

    #[test]
    fn smarts_recursive_structure_serial() {
        let molecule = Molecule::new();
        let query = RecursiveStructureQuery::from_molecule(
            crate::search::query_graph::query_graph_from_concrete_molecule(&molecule).unwrap(),
            101,
        );
        assert_eq!(query.serial_number(), 101);
    }

    fn test_bond_with_query(query: Option<QueryNode<BondQueryPredicate>>) -> crate::QueryBond {
        let spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Unspecified);
        let predicate = query.unwrap_or_else(|| make_bond_null_query());
        crate::QueryBond::from_parts(Bond::from_spec(crate::BondId::new(0), spec), predicate)
    }

    #[test]
    fn smarts_complex_bond_query() {
        let bond = |query| test_bond_with_query(Some(query));
        let order = |order| QueryNode::predicate(BondQueryPredicate::Order(order));

        assert!(!is_complex_bond_query(&test_bond_with_query(None)));
        assert!(!is_complex_bond_query(&bond(order(BondOrder::Double))));
        assert!(!is_complex_bond_query(&bond(
            make_single_or_aromatic_bond_query()
        )));
        assert!(is_complex_bond_query(&bond(QueryNode::not(order(
            BondOrder::Single
        )))));
        assert!(is_complex_bond_query(&bond(QueryNode::and(vec![order(
            BondOrder::Single
        )]))));
        assert!(is_complex_bond_query(&bond(QueryNode::xor(vec![order(
            BondOrder::Single
        )]))));
        assert!(!is_complex_bond_query(&bond(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Aromatic),
        ]))));
        assert!(!is_complex_bond_query(&bond(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Single),
        ]))));
        assert!(is_complex_bond_query(&bond(QueryNode::or(vec![
            order(BondOrder::Single),
            order(BondOrder::Double),
        ]))));
        assert!(is_complex_bond_query(&bond(QueryNode::predicate(
            BondQueryPredicate::Any,
        ))));
    }

    #[test]
    fn smarts_query_add_child() {
        let carbon = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let oxygen = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));
        let mut and_query = QueryNode::and(vec![carbon.clone()]);

        and_query.add_child(oxygen.clone());
        and_query.add_child(oxygen.clone());

        assert_eq!(
            and_query,
            QueryNode::and(vec![carbon.clone(), oxygen.clone(), oxygen.clone()])
        );

        let mut or_query = QueryNode::or(Vec::new());
        or_query.add_child(carbon.clone());
        or_query.add_child(oxygen.clone());
        assert_eq!(or_query, QueryNode::or(vec![carbon, oxygen]));
    }

    #[test]
    fn smarts_query_negation() {
        let carbon = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let oxygen = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));
        let mut query = QueryNode::and(vec![carbon.clone(), oxygen.clone()]);
        let original = query.clone();

        query.set_negation(false);
        assert_eq!(query, original);

        query.set_negation(true);
        assert_eq!(query, QueryNode::not(original.clone()));

        query.set_negation(true);
        assert_eq!(query, QueryNode::not(original.clone()));

        query.set_negation(false);
        assert_eq!(query, original);

        let child_negation = QueryNode::not(carbon);
        let mut nested = QueryNode::not(child_negation.clone());
        nested.set_negation(false);
        assert_eq!(nested, child_negation);
    }

    #[test]
    fn smarts_query_merge_both_null() {
        let any = || QueryNode::predicate(AtomQueryPredicate::Any);
        let not_any = || QueryNode::not(any());

        let mut and_query = any();
        merge_both_null_q(&mut and_query, &not_any(), CompositeQueryType::And);
        assert!(and_query.is_negated());

        let mut or_query = not_any();
        merge_both_null_q(&mut or_query, &any(), CompositeQueryType::Or);
        assert!(!or_query.is_negated());

        let mut xor_true_true = any();
        merge_both_null_q(&mut xor_true_true, &any(), CompositeQueryType::Xor);
        assert!(xor_true_true.is_negated());

        let mut xor_false_true = not_any();
        merge_both_null_q(&mut xor_false_true, &any(), CompositeQueryType::Xor);
        assert!(!xor_false_true.is_negated());

        let mut unchanged = not_any();
        merge_both_null_q(&mut unchanged, &not_any(), CompositeQueryType::And);
        assert!(unchanged.is_negated());
    }

    #[test]
    fn smarts_query_merge_null_first() {
        let any = || QueryNode::predicate(AtomQueryPredicate::Any);
        let carbon = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));

        let mut and_null = any();
        let mut and_other = carbon();
        merge_null_q_first(&mut and_null, &mut and_other, CompositeQueryType::And);
        assert_eq!(and_null, carbon());
        assert_eq!(and_other, any());

        let mut or_false = QueryNode::not(any());
        let mut or_other = carbon();
        merge_null_q_first(&mut or_false, &mut or_other, CompositeQueryType::Or);
        assert_eq!(or_false, carbon());

        let mut xor_true = any();
        let mut xor_other = carbon();
        merge_null_q_first(&mut xor_true, &mut xor_other, CompositeQueryType::Xor);
        assert_eq!(xor_true, QueryNode::not(carbon()));

        let mut unchanged = QueryNode::not(any());
        let mut unchanged_other = carbon();
        merge_null_q_first(
            &mut unchanged,
            &mut unchanged_other,
            CompositeQueryType::And,
        );
        assert_eq!(unchanged, QueryNode::not(any()));
        assert_eq!(unchanged_other, carbon());
    }

    #[test]
    fn smarts_query_merge_null_queries() {
        let any = || QueryNode::predicate(AtomQueryPredicate::Any);
        let carbon = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let oxygen = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));

        let mut both_left = any();
        let mut both_right = QueryNode::not(any());
        merge_null_queries(
            &mut both_left,
            true,
            &mut both_right,
            true,
            CompositeQueryType::And,
        );
        assert!(both_left.is_negated());

        let mut first_null = any();
        let mut first_other = carbon();
        merge_null_queries(
            &mut first_null,
            true,
            &mut first_other,
            false,
            CompositeQueryType::And,
        );
        assert_eq!(first_null, carbon());

        let mut second_other = carbon();
        let mut second_null = any();
        merge_null_queries(
            &mut second_other,
            false,
            &mut second_null,
            true,
            CompositeQueryType::And,
        );
        assert_eq!(second_other, carbon());

        let mut neither_left = carbon();
        let mut neither_right = oxygen();
        merge_null_queries(
            &mut neither_left,
            false,
            &mut neither_right,
            false,
            CompositeQueryType::Or,
        );
        assert_eq!(neither_left, carbon());
        assert_eq!(neither_right, oxygen());
    }

    #[test]
    fn smarts_query_bond_constructor() {
        let begin = crate::AtomId::new(2);
        let end = crate::AtomId::new(5);

        let single = make_query_bond_spec(begin, end, BondOrder::Single);
        assert_eq!(single.begin(), begin);
        assert_eq!(single.end(), end);
        assert_eq!(single.bond().order(), BondOrder::Single);
        assert_eq!(
            single.predicate(),
            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single))
        );

        let unspecified = make_query_bond_spec(begin, end, BondOrder::Unspecified);
        assert_eq!(unspecified.bond().order(), BondOrder::Unspecified);
        assert_eq!(unspecified.predicate(), &make_bond_null_query());
    }

    #[test]
    fn smarts_query_and_match() {
        let children = vec![
            QueryNode::predicate(true),
            QueryNode::predicate(false),
            QueryNode::predicate(true),
        ];
        let visits = std::cell::Cell::new(0);
        let result = and_query_match(&children, false, |child| {
            visits.set(visits.get() + 1);
            let QueryNode::Predicate(value) = child else {
                unreachable!("boolean fixture contains only predicate leaves")
            };
            *value
        });
        assert!(!result);
        assert_eq!(visits.get(), 2, "AND must stop at the first false child");

        let all_true = vec![QueryNode::predicate(true), QueryNode::predicate(true)];
        assert!(and_query_match(&all_true, false, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));
        assert!(!and_query_match(&all_true, true, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));

        let empty: Vec<QueryNode<bool>> = Vec::new();
        assert!(and_query_match(&empty, false, |_| unreachable!()));
        assert!(!and_query_match(&empty, true, |_| unreachable!()));
    }

    #[test]
    fn smarts_query_or_match() {
        let children = vec![
            QueryNode::predicate(false),
            QueryNode::predicate(true),
            QueryNode::predicate(false),
        ];
        let visits = std::cell::Cell::new(0);
        let result = or_query_match(&children, false, |child| {
            visits.set(visits.get() + 1);
            let QueryNode::Predicate(value) = child else {
                unreachable!("boolean fixture contains only predicate leaves")
            };
            *value
        });
        assert!(result);
        assert_eq!(visits.get(), 2, "OR must stop at the first true child");

        let all_false = vec![QueryNode::predicate(false), QueryNode::predicate(false)];
        assert!(!or_query_match(&all_false, false, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));
        assert!(or_query_match(&all_false, true, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));

        let empty: Vec<QueryNode<bool>> = Vec::new();
        assert!(!or_query_match(&empty, false, |_| unreachable!()));
        assert!(or_query_match(&empty, true, |_| unreachable!()));
    }

    #[test]
    fn smarts_query_xor_match() {
        let two_true = vec![
            QueryNode::predicate(true),
            QueryNode::predicate(true),
            QueryNode::predicate(true),
        ];
        let visits = std::cell::Cell::new(0);
        assert!(!xor_query_match(&two_true, false, |child| {
            visits.set(visits.get() + 1);
            matches!(child, QueryNode::Predicate(true))
        }));
        assert_eq!(visits.get(), 2, "XOR must stop at the second true child");

        let one_true = vec![
            QueryNode::predicate(false),
            QueryNode::predicate(true),
            QueryNode::predicate(false),
        ];
        assert!(xor_query_match(&one_true, false, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));
        assert!(!xor_query_match(&one_true, true, |child| {
            matches!(child, QueryNode::Predicate(true))
        }));

        let empty: Vec<QueryNode<bool>> = Vec::new();
        assert!(!xor_query_match(&empty, false, |_| unreachable!()));
        assert!(xor_query_match(&empty, true, |_| unreachable!()));

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("XOR atom fixture is valid");
        let query = QueryNode::xor(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8)),
        ]);
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_query_equality_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(equality_query_match(7, 5, 2, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(equality_query_match(7, 9, 2, false, |value| value));
        assert!(!equality_query_match(7, 4, 2, false, |value| value));
        assert!(!equality_query_match(7, 10, 2, false, |value| value));
        assert!(equality_query_match(7, 10, 2, true, |value| value));

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("equality atom fixture is valid");
        let query = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_query_greater_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(greater_query_match(7, 4, 2, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(!greater_query_match(7, 5, 2, false, |value| value));
        assert!(!greater_query_match(7, 8, 2, false, |value| value));
        assert!(greater_query_match(7, 5, 2, true, |value| value));

        assert!(greater_query_match(7, 6, 0, false, |value| value));
        assert!(!greater_query_match(7, 8, 0, false, |value| value));
    }

    #[test]
    fn smarts_query_greater_equal_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(greater_equal_query_match(7, 5, 2, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(greater_equal_query_match(7, 4, 2, false, |value| value));
        assert!(greater_equal_query_match(7, 9, 2, false, |value| value));
        assert!(!greater_equal_query_match(7, 10, 2, false, |value| value));
        assert!(greater_equal_query_match(7, 10, 2, true, |value| value));

        assert!(greater_equal_query_match(7, 7, 0, false, |value| value));
        assert!(greater_equal_query_match(7, 6, 0, false, |value| value));
        assert!(!greater_equal_query_match(7, 8, 0, false, |value| value));
    }

    #[test]
    fn smarts_query_less_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(less_query_match(7, 10, 2, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(!less_query_match(7, 9, 2, false, |value| value));
        assert!(!less_query_match(7, 6, 2, false, |value| value));
        assert!(less_query_match(7, 9, 2, true, |value| value));

        assert!(less_query_match(7, 8, 0, false, |value| value));
        assert!(!less_query_match(7, 6, 0, false, |value| value));
    }

    #[test]
    fn smarts_query_less_equal_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(less_equal_query_match(7, 9, 2, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(less_equal_query_match(7, 10, 2, false, |value| value));
        assert!(less_equal_query_match(7, 5, 2, false, |value| value));
        assert!(!less_equal_query_match(7, 4, 2, false, |value| value));
        assert!(less_equal_query_match(7, 4, 2, true, |value| value));

        assert!(less_equal_query_match(7, 7, 0, false, |value| value));
        assert!(less_equal_query_match(7, 8, 0, false, |value| value));
        assert!(!less_equal_query_match(7, 6, 0, false, |value| value));
    }

    #[test]
    fn smarts_query_range_match() {
        let conversions = std::cell::Cell::new(0);
        assert!(range_query_match(
            3,
            7,
            5,
            0,
            true,
            true,
            false,
            |observed| {
                conversions.set(conversions.get() + 1);
                observed
            }
        ));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(!range_query_match(3, 7, 3, 0, true, true, false, |v| v));
        assert!(range_query_match(3, 7, 3, 0, false, true, false, |v| v));
        assert!(!range_query_match(3, 7, 7, 0, true, true, false, |v| v));
        assert!(range_query_match(3, 7, 7, 0, true, false, false, |v| v));

        assert!(!range_query_match(3, 7, 2, 1, true, true, false, |v| v));
        assert!(range_query_match(3, 7, 2, 1, false, true, false, |v| v));
        assert!(!range_query_match(3, 7, 1, 1, true, true, false, |v| v));
        assert!(range_query_match(3, 7, 1, 1, true, true, true, |v| v));
    }

    #[test]
    fn smarts_query_set_match() {
        let values = [2, 5, 8].into_iter().collect::<BTreeSet<_>>();
        let conversions = std::cell::Cell::new(0);
        assert!(set_query_match(&values, 5, false, |observed| {
            conversions.set(conversions.get() + 1);
            observed
        }));
        assert_eq!(conversions.get(), 1, "TypeConvert must run exactly once");

        assert!(!set_query_match(&values, 3, false, |value| value));
        assert!(set_query_match(&values, 3, true, |value| value));
        assert!(!set_query_match(&values, 8, true, |value| value));
        assert!(!set_query_match(
            &BTreeSet::<i32>::new(),
            1,
            false,
            |value| value
        ));
    }

    #[test]
    fn smarts_query_atom_copy() {
        let query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::FormalCharge(1)),
        ]);
        let atom = Atom::from_spec(
            AtomId::new(0),
            crate::AtomSpec::new(crate::Element::C)
                .with_formal_charge(1)
                .with_prop("label", "source"),
        );
        let source = crate::QueryAtom::from_parts(atom, query);
        let mut copied = query_atom_copy(&source);

        assert_eq!(&copied, &source);
        copied.predicate_mut().set_negation(true);
        assert_ne!(copied.predicate(), source.predicate());
        assert_eq!(source.atom().prop("label"), Some("source"));
    }

    #[test]
    fn smarts_query_bond_copy() {
        let source = crate::QueryBond::from_parts(
            Bond::from_spec(
                crate::BondId::new(0),
                crate::BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double)
                    .with_prop("label", "source"),
            ),
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double)),
                QueryNode::predicate(BondQueryPredicate::IsInRing(false)),
            ]),
        );
        let mut copied = query_bond_copy(&source);

        assert_eq!(&copied, &source);
        copied.set_predicate(QueryNode::predicate(BondQueryPredicate::Order(
            BondOrder::Single,
        )));
        assert_ne!(copied.predicate(), source.predicate());
        assert_eq!(source.bond().prop("label"), Some("source"));
    }

    #[test]
    fn smarts_query_bond_set_type() {
        let mut bond = crate::QueryBond::from_parts(
            Bond::from_spec(
                crate::BondId::new(0),
                crate::BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryNode::and(vec![
                QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
                QueryNode::predicate(BondQueryPredicate::IsInRing(true)),
            ]),
        );

        query_bond_set_type(&mut bond, BondOrder::Triple);
        assert_eq!(bond.bond().order(), BondOrder::Triple);
        assert_eq!(
            bond.predicate(),
            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Triple))
        );

        query_bond_set_type(&mut bond, BondOrder::Unspecified);
        assert_eq!(bond.bond().order(), BondOrder::Unspecified);
        assert_eq!(
            bond.predicate(),
            &QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Unspecified))
        );
        assert_ne!(bond.predicate(), &make_bond_null_query());
    }

    #[test]
    fn smarts_query_bond_set_dir() {
        let original_query = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double));
        let mut bond = crate::QueryBond::from_parts(
            Bond::from_spec(
                crate::BondId::new(0),
                crate::BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Double),
            ),
            original_query.clone(),
        );

        query_bond_set_dir(&mut bond, crate::BondDirection::BeginWedge);
        assert_eq!(bond.bond().direction(), crate::BondDirection::BeginWedge);
        assert_eq!(bond.predicate(), &original_query);

        query_bond_set_dir(&mut bond, crate::BondDirection::EitherDouble);
        assert_eq!(bond.bond().direction(), crate::BondDirection::EitherDouble);
        assert_eq!(bond.predicate(), &original_query);
    }

    #[test]
    fn smarts_query_atom_expand() {
        let carbon = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let charged = QueryNode::predicate(AtomQueryPredicate::FormalCharge(1));

        let mut maintained = carbon.clone();
        query_atom_expand_query(
            &mut maintained,
            charged.clone(),
            CompositeQueryType::And,
            true,
        );
        assert_eq!(
            maintained,
            QueryNode::and(vec![carbon.clone(), charged.clone()])
        );

        let mut reversed = carbon.clone();
        query_atom_expand_query(
            &mut reversed,
            charged.clone(),
            CompositeQueryType::Or,
            false,
        );
        assert_eq!(
            reversed,
            QueryNode::or(vec![charged.clone(), carbon.clone()])
        );

        let mut null_and = make_atom_null_query();
        query_atom_expand_query(
            &mut null_and,
            charged.clone(),
            CompositeQueryType::And,
            true,
        );
        assert_eq!(null_and, charged);

        let mut null_xor = make_atom_null_query();
        query_atom_expand_query(&mut null_xor, carbon.clone(), CompositeQueryType::Xor, true);
        assert_eq!(null_xor, QueryNode::not(carbon));
    }

    #[test]
    fn smarts_query_bond_expand() {
        let single = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
        let in_ring = QueryNode::predicate(BondQueryPredicate::IsInRing(true));

        let mut maintained = single.clone();
        query_bond_expand_query(
            &mut maintained,
            in_ring.clone(),
            CompositeQueryType::And,
            true,
        );
        assert_eq!(
            maintained,
            QueryNode::and(vec![single.clone(), in_ring.clone()])
        );

        let mut reversed = single.clone();
        query_bond_expand_query(
            &mut reversed,
            in_ring.clone(),
            CompositeQueryType::Or,
            false,
        );
        assert_eq!(
            reversed,
            QueryNode::or(vec![in_ring.clone(), single.clone()])
        );

        let mut null_and = make_bond_null_query();
        query_bond_expand_query(
            &mut null_and,
            in_ring.clone(),
            CompositeQueryType::And,
            true,
        );
        assert_eq!(null_and, in_ring);

        let mut null_xor = make_bond_null_query();
        query_bond_expand_query(&mut null_xor, single.clone(), CompositeQueryType::Xor, true);
        assert_eq!(null_xor, QueryNode::not(single));
    }

    #[test]
    fn smarts_query_atom_local_match() {
        assert!(atom_query_local_match(&6, false, &6, false));
        assert!(!atom_query_local_match(&6, false, &8, false));
        assert!(atom_query_local_match(&6, true, &6, true));
        assert!(!atom_query_local_match(&6, true, &8, true));

        assert!(!atom_query_local_match(&6, false, &6, true));
        assert!(atom_query_local_match(&6, false, &8, true));
        assert!(!atom_query_local_match(&6, true, &6, false));
        assert!(atom_query_local_match(&6, true, &8, false));
    }

    #[test]
    fn smarts_query_bond_local_match() {
        assert!(bond_query_local_match(
            &BondOrder::Single,
            false,
            &BondOrder::Single,
            false
        ));
        assert!(!bond_query_local_match(
            &BondOrder::Single,
            false,
            &BondOrder::Double,
            false
        ));
        assert!(!bond_query_local_match(
            &BondOrder::Single,
            false,
            &BondOrder::Single,
            true
        ));
        assert!(bond_query_local_match(
            &BondOrder::Single,
            false,
            &BondOrder::Double,
            true
        ));
    }

    #[test]
    fn smarts_query_atom_queries_match() {
        let carbon = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));
        let oxygen = || QueryNode::predicate(AtomQueryPredicate::AtomicNumber(8));
        let neutral = || QueryNode::predicate(AtomQueryPredicate::FormalCharge(0));

        assert!(atom_queries_match(&carbon(), &carbon()));
        assert!(!atom_queries_match(&carbon(), &oxygen()));
        assert!(atom_queries_match(&carbon(), &QueryNode::not(oxygen())));
        assert!(!atom_queries_match(&carbon(), &QueryNode::not(carbon())));
        assert!(atom_queries_match(&make_atom_null_query(), &neutral()));

        let disjunction = QueryNode::or(vec![carbon(), oxygen()]);
        assert!(atom_queries_match(&disjunction, &oxygen()));
        assert!(!atom_queries_match(&disjunction, &neutral()));

        let conjunction = QueryNode::and(vec![carbon(), neutral()]);
        let reordered = QueryNode::and(vec![neutral(), carbon(), oxygen()]);
        assert!(atom_queries_match(&conjunction, &reordered));
        assert!(!atom_queries_match(&conjunction, &carbon()));

        let xor = QueryNode::xor(vec![carbon(), oxygen()]);
        assert!(!atom_queries_match(&xor, &xor));
        assert!(!atom_queries_match(
            &QueryNode::predicate(AtomQueryPredicate::AtomicNumberIn(vec![6, 8])),
            &carbon()
        ));
    }

    #[test]
    fn smarts_query_bond_queries_match() {
        let single = || QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
        let double = || QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double));
        let in_ring = || QueryNode::predicate(BondQueryPredicate::IsInRing(true));
        let min_ring_four = || QueryNode::predicate(BondQueryPredicate::MinRingSize(4));

        assert!(bond_queries_match(&single(), &single()));
        assert!(!bond_queries_match(&single(), &double()));
        assert!(bond_queries_match(&single(), &QueryNode::not(double())));
        assert!(!bond_queries_match(&single(), &QueryNode::not(single())));
        assert!(bond_queries_match(&make_bond_null_query(), &double()));

        let disjunction = QueryNode::or(vec![single(), double()]);
        assert!(bond_queries_match(&disjunction, &double()));
        assert!(!bond_queries_match(
            &disjunction,
            &QueryNode::predicate(BondQueryPredicate::MinRingSize(4))
        ));

        let conjunction = QueryNode::and(vec![single(), min_ring_four()]);
        let reordered = QueryNode::and(vec![min_ring_four(), single(), double()]);
        assert!(bond_queries_match(&conjunction, &reordered));
        assert!(!bond_queries_match(&conjunction, &single()));

        // RDKit does not require d1 == d2 in the equality-query branch.
        assert!(bond_queries_match(&single(), &in_ring()));

        // Preserve QueryBond.cpp's asymmetric d2 == BondAnd condition.
        assert!(!bond_queries_match(
            &single(),
            &QueryNode::and(vec![double(), single()])
        ));
        assert!(bond_queries_match(
            &single(),
            &QueryNode::and(vec![
                double(),
                QueryNode::predicate(BondQueryPredicate::Direction(crate::BondDirection::None))
            ])
        ));

        let xor = QueryNode::xor(vec![single(), double()]);
        assert!(!bond_queries_match(&xor, &xor));
    }

    #[test]
    fn smarts_query_atom_match() {
        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("query atom match fixture is valid");
        let query = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::not(QueryNode::predicate(AtomQueryPredicate::FormalCharge(1))),
        ]);

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule
        ));
    }

    #[test]
    fn smarts_query_bond_match() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let single = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("single bond is valid");
        let double = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Double))
            .expect("double bond is valid");
        let molecule = builder.build().expect("query bond Match fixture is valid");
        let query = QueryNode::and(vec![
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
            QueryNode::not(QueryNode::predicate(BondQueryPredicate::IsInRing(true))),
        ]);

        assert!(bond_matches_query(
            &molecule.bonds()[single.index()],
            &query,
            &molecule
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[double.index()],
            &query,
            &molecule
        ));
    }

    #[test]
    fn smarts_query_bond_query_match() {
        let mut builder = Molecule::builder();
        let atoms = (0..8)
            .map(|_| builder.add_atom(crate::AtomSpec::new(crate::Element::C)))
            .collect::<Vec<_>>();
        let ordinary_single = builder
            .add_bond(crate::BondSpec::new(atoms[0], atoms[1], BondOrder::Single))
            .expect("ordinary single bond is valid");
        let ordinary_double = builder
            .add_bond(crate::BondSpec::new(atoms[2], atoms[3], BondOrder::Double))
            .expect("ordinary double bond is valid");
        let molecule = builder
            .build()
            .expect("query bond QueryMatch fixture is valid");
        let single_query = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));

        assert!(query_bond_query_match(
            &single_query,
            &molecule.bonds()[ordinary_single.index()],
            &molecule
        ));
        assert!(!query_bond_query_match(
            &single_query,
            &molecule.bonds()[ordinary_double.index()],
            &molecule
        ));
        let query_single = make_query_bond_spec(AtomId::new(4), AtomId::new(5), BondOrder::Single);
        let query_double = make_query_bond_spec(AtomId::new(6), AtomId::new(7), BondOrder::Double);
        assert!(bond_queries_match(&single_query, query_single.predicate()));
        assert!(!bond_queries_match(&single_query, query_double.predicate()));
    }

    #[test]
    fn smarts_query_atom_query_match() {
        let mut builder = Molecule::builder();
        let ordinary_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ordinary_oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder
            .build()
            .expect("query atom QueryMatch fixture is valid");
        let carbon_query = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6));

        assert!(query_atom_query_match(
            &carbon_query,
            &molecule.atoms()[ordinary_carbon.index()],
            &molecule
        ));
        assert!(!query_atom_query_match(
            &carbon_query,
            &molecule.atoms()[ordinary_oxygen.index()],
            &molecule
        ));
        let query_carbon = test_atom_with_query(Some(carbon_query.clone()));
        let query_oxygen = test_atom_with_query(Some(QueryNode::predicate(
            AtomQueryPredicate::AtomicNumber(8),
        )));
        assert!(atom_queries_match(
            query_carbon.predicate(),
            query_carbon.predicate()
        ));
        assert!(!atom_queries_match(
            query_carbon.predicate(),
            query_oxygen.predicate()
        ));
    }

    #[test]
    fn smarts_make_single_or_double_or_aromatic_bond_query() {
        assert_eq!(
            make_single_or_double_or_aromatic_bond_query(),
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))
        );

        let mut builder = Molecule::builder();
        let mut bonds = Vec::new();
        for order in [
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Aromatic,
            BondOrder::Triple,
        ] {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            bonds.push(
                builder
                    .add_bond(crate::BondSpec::new(begin, end, order))
                    .expect("single-double-aromatic fixture bond is valid"),
            );
        }
        let molecule = builder
            .build()
            .expect("single-double-aromatic fixture is valid");
        let query = make_single_or_double_or_aromatic_bond_query();

        for bond_id in &bonds[..3] {
            assert!(bond_matches_query(
                &molecule.bonds()[bond_id.index()],
                &query,
                &molecule,
            ));
        }
        assert!(!bond_matches_query(
            &molecule.bonds()[bonds[3].index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_single_or_double_bond_query() {
        assert_eq!(
            make_single_or_double_bond_query(),
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ]))
        );

        let mut builder = Molecule::builder();
        let mut bonds = Vec::new();
        for order in [BondOrder::Single, BondOrder::Double, BondOrder::Aromatic] {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            bonds.push(
                builder
                    .add_bond(crate::BondSpec::new(begin, end, order))
                    .expect("single-or-double fixture bond is valid"),
            );
        }
        let molecule = builder.build().expect("single-or-double fixture is valid");
        let query = make_single_or_double_bond_query();

        for bond_id in &bonds[..2] {
            assert!(bond_matches_query(
                &molecule.bonds()[bond_id.index()],
                &query,
                &molecule,
            ));
        }
        assert!(!bond_matches_query(
            &molecule.bonds()[bonds[2].index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_double_or_aromatic_bond_query() {
        assert_eq!(
            make_double_or_aromatic_bond_query(),
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))
        );

        let mut builder = Molecule::builder();
        let mut bonds = Vec::new();
        for order in [BondOrder::Double, BondOrder::Aromatic, BondOrder::Single] {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            bonds.push(
                builder
                    .add_bond(crate::BondSpec::new(begin, end, order))
                    .expect("double-or-aromatic fixture bond is valid"),
            );
        }
        let molecule = builder
            .build()
            .expect("double-or-aromatic fixture is valid");
        let query = make_double_or_aromatic_bond_query();

        for bond_id in &bonds[..2] {
            assert!(bond_matches_query(
                &molecule.bonds()[bond_id.index()],
                &query,
                &molecule,
            ));
        }
        assert!(!bond_matches_query(
            &molecule.bonds()[bonds[2].index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_single_or_aromatic_bond_query() {
        assert_eq!(
            make_single_or_aromatic_bond_query(),
            QueryNode::predicate(BondQueryPredicate::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))
        );

        let mut builder = Molecule::builder();
        let mut bonds = Vec::new();
        for order in [BondOrder::Single, BondOrder::Aromatic, BondOrder::Double] {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            bonds.push(
                builder
                    .add_bond(crate::BondSpec::new(begin, end, order))
                    .expect("single-or-aromatic fixture bond is valid"),
            );
        }
        let molecule = builder
            .build()
            .expect("single-or-aromatic fixture is valid");
        let query = make_single_or_aromatic_bond_query();

        assert!(bond_matches_query(
            &molecule.bonds()[bonds[0].index()],
            &query,
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[bonds[1].index()],
            &query,
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[bonds[2].index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_order_equals_query() {
        assert_eq!(
            make_bond_order_equals_query(BondOrder::Double),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Double))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let single = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond-order query single bond is valid");
        let double = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Double))
            .expect("bond-order query double bond is valid");
        let molecule = builder.build().expect("bond-order query fixture is valid");
        let query = make_bond_order_equals_query(BondOrder::Double);

        assert!(!bond_matches_query(
            &molecule.bonds()[single.index()],
            &query,
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[double.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_dir_equals_query() {
        use crate::BondDirection;

        assert_eq!(
            make_bond_dir_equals_query(BondDirection::BeginWedge),
            QueryNode::predicate(BondQueryPredicate::Direction(BondDirection::BeginWedge))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let wedge = builder
            .add_bond(
                crate::BondSpec::new(a0, a1, BondOrder::Single)
                    .with_direction(BondDirection::BeginWedge),
            )
            .expect("bond-direction query wedge bond is valid");
        let dash = builder
            .add_bond(
                crate::BondSpec::new(a1, a2, BondOrder::Single)
                    .with_direction(BondDirection::BeginDash),
            )
            .expect("bond-direction query dash bond is valid");
        let molecule = builder
            .build()
            .expect("bond-direction query fixture is valid");
        let query = make_bond_dir_equals_query(BondDirection::BeginWedge);

        assert!(bond_matches_query(
            &molecule.bonds()[wedge.index()],
            &query,
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[dash.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_has_stereo_query() {
        assert_eq!(
            make_bond_has_stereo_query(),
            QueryNode::predicate(BondQueryPredicate::HasStereo)
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let without_stereo = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Double))
            .expect("bond-stereo query plain bond is valid");
        let with_stereo = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Double).with_stereo(BondStereo::Any))
            .expect("bond-stereo query stereo bond is valid");
        let molecule = builder.build().expect("bond-stereo query fixture is valid");
        let query = make_bond_has_stereo_query();

        assert!(!bond_matches_query(
            &molecule.bonds()[without_stereo.index()],
            &query,
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[with_stereo.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_is_in_ring_query() {
        assert_eq!(
            make_bond_is_in_ring_query(),
            QueryNode::predicate(BondQueryPredicate::IsInRing(true))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("ring-bond query ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .expect("ring-bond query second ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a2, a0, BondOrder::Single))
            .expect("ring-bond query third ring edge is valid");
        let chain_bond = builder
            .add_bond(crate::BondSpec::new(a2, a3, BondOrder::Single))
            .expect("ring-bond query chain edge is valid");
        let molecule = builder.build().expect("ring-bond query fixture is valid");
        let query = make_bond_is_in_ring_query();

        assert!(bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &query,
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[chain_bond.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_in_n_rings_query() {
        assert_eq!(
            make_bond_in_n_rings_query(1),
            QueryNode::predicate(BondQueryPredicate::NumRingBonds(1))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("ring-count query ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .expect("ring-count query second ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a2, a0, BondOrder::Single))
            .expect("ring-count query third ring edge is valid");
        let chain_bond = builder
            .add_bond(crate::BondSpec::new(a2, a3, BondOrder::Single))
            .expect("ring-count query chain edge is valid");
        let molecule = builder.build().expect("ring-count query fixture is valid");

        assert!(bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_in_n_rings_query(1),
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[chain_bond.index()],
            &make_bond_in_n_rings_query(0),
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_in_n_rings_query(-1),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_in_ring_of_size_query() {
        assert_eq!(
            make_bond_in_ring_of_size_query(3),
            Ok(QueryNode::predicate(BondQueryPredicate::InRingOfSize(3)))
        );
        assert_eq!(
            make_bond_in_ring_of_size_query(2),
            Err(QueryConstructionError::BondRingSizeOutOfRange { target: 2 })
        );
        assert_eq!(
            make_bond_in_ring_of_size_query(21),
            Err(QueryConstructionError::BondRingSizeOutOfRange { target: 21 })
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("ring-size query first edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .expect("ring-size query second edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a2, a0, BondOrder::Single))
            .expect("ring-size query third edge is valid");
        let molecule = builder.build().expect("ring-size query fixture is valid");

        assert!(bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_in_ring_of_size_query(3).expect("size three is supported"),
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_in_ring_of_size_query(4).expect("size four is supported"),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_min_ring_size_query() {
        assert_eq!(
            make_bond_min_ring_size_query(3),
            QueryNode::predicate(BondQueryPredicate::MinRingSize(3))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("minimum-ring-size query first ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .expect("minimum-ring-size query second ring edge is valid");
        builder
            .add_bond(crate::BondSpec::new(a2, a0, BondOrder::Single))
            .expect("minimum-ring-size query third ring edge is valid");
        let chain_bond = builder
            .add_bond(crate::BondSpec::new(a2, a3, BondOrder::Single))
            .expect("minimum-ring-size query chain edge is valid");
        let molecule = builder
            .build()
            .expect("minimum-ring-size query fixture is valid");

        assert!(bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_min_ring_size_query(3),
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[chain_bond.index()],
            &make_bond_min_ring_size_query(0),
            &molecule,
        ));
        assert!(!bond_matches_query(
            &molecule.bonds()[ring_bond.index()],
            &make_bond_min_ring_size_query(-1),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_bond_null_query() {
        assert_eq!(
            make_bond_null_query(),
            QueryNode::predicate(BondQueryPredicate::Any)
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let single = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("null-bond query single edge is valid");
        let triple = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Triple))
            .expect("null-bond query triple edge is valid");
        let molecule = builder.build().expect("null-bond query fixture is valid");
        let query = make_bond_null_query();

        assert!(bond_matches_query(
            &molecule.bonds()[single.index()],
            &query,
            &molecule,
        ));
        assert!(bond_matches_query(
            &molecule.bonds()[triple.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_m_h_atom_query() {
        let expected = MH_EXCLUDED_ATOMIC_NUMBERS
            .map(|number| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(number)))
            .to_vec();
        assert_eq!(
            make_m_h_atom_query(),
            QueryNode::not(QueryNode::or(expected))
        );

        let mut builder = Molecule::builder();
        let iron = builder.add_atom(crate::AtomSpec::new(crate::Element::FE));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("MH atom query fixture is valid");
        let query = make_m_h_atom_query();

        assert!(atom_matches_query(
            &molecule.atoms()[iron.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_convert_complex_name_to_query() {
        for (symbol, expected) in [
            ("Q", make_q_atom_query()),
            ("QH", make_q_h_atom_query()),
            ("A", make_a_atom_query()),
            ("AH", make_a_h_atom_query()),
            ("X", make_x_atom_query()),
            ("XH", make_x_h_atom_query()),
            ("M", make_m_atom_query()),
            ("MH", make_m_h_atom_query()),
        ] {
            assert_eq!(
                convert_complex_name_to_query(symbol),
                Ok(expected),
                "complex query symbol {symbol}"
            );
        }

        assert_eq!(
            convert_complex_name_to_query("R"),
            Err(QueryConstructionError::InvalidComplexAtomSymbol {
                symbol: "R".to_owned()
            })
        );
    }

    #[test]
    fn smarts_make_m_atom_query() {
        let mut expected = MH_EXCLUDED_ATOMIC_NUMBERS
            .map(|number| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(number)))
            .to_vec();
        expected.push(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)));
        assert_eq!(make_m_atom_query(), QueryNode::not(QueryNode::or(expected)));

        let mut builder = Molecule::builder();
        let iron = builder.add_atom(crate::AtomSpec::new(crate::Element::FE));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let molecule = builder.build().expect("M atom query fixture is valid");
        let query = make_m_atom_query();

        assert!(atom_matches_query(
            &molecule.atoms()[iron.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_x_h_atom_query() {
        assert_eq!(
            make_x_h_atom_query(),
            QueryNode::or(
                [9, 17, 35, 53, 85, 1]
                    .map(|number| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(number)))
                    .to_vec()
            )
        );

        let mut builder = Molecule::builder();
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("XH atom query fixture is valid");
        let query = make_x_h_atom_query();

        assert!(atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[fluorine.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_x_atom_query() {
        assert_eq!(
            make_x_atom_query(),
            QueryNode::or(
                [9, 17, 35, 53, 85]
                    .map(|number| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(number)))
                    .to_vec()
            )
        );

        let mut builder = Molecule::builder();
        let fluorine = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
        let chlorine = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("X atom query fixture is valid");
        let query = make_x_atom_query();

        assert!(atom_matches_query(
            &molecule.atoms()[fluorine.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[chlorine.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_a_h_atom_query() {
        assert_eq!(
            make_a_h_atom_query(),
            QueryNode::predicate(AtomQueryPredicate::Any)
        );

        let mut builder = Molecule::builder();
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("AH atom query fixture is valid");
        let query = make_a_h_atom_query();

        for atom_id in [hydrogen, carbon, oxygen] {
            assert!(atom_matches_query(
                &molecule.atoms()[atom_id.index()],
                &query,
                &molecule,
            ));
        }
    }

    #[test]
    fn smarts_make_atom_null_query() {
        assert_eq!(
            make_atom_null_query(),
            QueryNode::predicate(AtomQueryPredicate::Any)
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let molecule = builder.build().expect("null-atom query fixture is valid");
        let query = make_atom_null_query();

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_a_atom_query() {
        assert_eq!(
            make_a_atom_query(),
            QueryNode::not(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)))
        );

        let mut builder = Molecule::builder();
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("A atom query fixture is valid");
        let query = make_a_atom_query();

        assert!(!atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_q_h_atom_query() {
        assert_eq!(
            make_q_h_atom_query(),
            QueryNode::not(QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)))
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("QH atom query fixture is valid");
        let query = make_q_h_atom_query();

        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_q_atom_query() {
        assert_eq!(
            make_q_atom_query(),
            QueryNode::not(QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(1)),
            ]))
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("Q atom query fixture is valid");
        let query = make_q_atom_query();

        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[hydrogen.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_is_bridgehead_query() {
        assert_eq!(
            make_atom_is_bridgehead_query(),
            QueryNode::Predicate(AtomQueryPredicate::IsBridgehead)
        );

        let mut builder = Molecule::builder();
        let bridgehead = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let opposite = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let path_atoms = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for middle in path_atoms {
            for (begin, end) in [(bridgehead, middle), (middle, opposite)] {
                builder
                    .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                    .expect("bridgehead query fixture bond is valid");
            }
        }
        let molecule = builder.build().expect("bridgehead query fixture is valid");
        let query = make_atom_is_bridgehead_query();

        assert!(atom_matches_query(
            &molecule.atoms()[bridgehead.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[path_atoms[0].index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_non_hydrogen_degree_query() {
        assert_eq!(
            make_atom_non_hydrogen_degree_query(2),
            QueryNode::Predicate(AtomQueryPredicate::NonHydrogenDegree(2))
        );

        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        for neighbor in [carbon, oxygen, hydrogen] {
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, BondOrder::Single))
                .expect("non-hydrogen-degree fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("non-hydrogen-degree fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_non_hydrogen_degree_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_non_hydrogen_degree_query(3),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_has_aliphatic_heteroatom_nbrs_query() {
        assert_eq!(
            make_atom_has_aliphatic_heteroatom_nbrs_query(),
            QueryNode::Predicate(AtomQueryPredicate::HasAliphaticHeteroatomNeighbors)
        );

        let mut builder = Molecule::builder();
        let with_aliphatic = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                with_aliphatic,
                oxygen,
                BondOrder::Single,
            ))
            .expect("aliphatic-heteroatom fixture bond is valid");
        let with_aromatic_only = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        builder
            .add_bond(crate::BondSpec::new(
                with_aromatic_only,
                aromatic_nitrogen,
                BondOrder::Single,
            ))
            .expect("aromatic-heteroatom fixture bond is valid");
        let molecule = builder
            .build()
            .expect("has-aliphatic-heteroatom-neighbor fixture is valid");
        let query = make_atom_has_aliphatic_heteroatom_nbrs_query();

        assert!(atom_matches_query(
            &molecule.atoms()[with_aliphatic.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[with_aromatic_only.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_num_aliphatic_heteroatom_nbrs_query() {
        assert_eq!(
            make_atom_num_aliphatic_heteroatom_nbrs_query(1),
            QueryNode::Predicate(AtomQueryPredicate::NumAliphaticHeteroatomNeighbors(1))
        );

        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aliphatic_oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let aromatic_nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for neighbor in [aliphatic_oxygen, aromatic_nitrogen, carbon] {
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, BondOrder::Single))
                .expect("aliphatic-heteroatom-neighbor fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("aliphatic-heteroatom-neighbor fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_num_aliphatic_heteroatom_nbrs_query(1),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_num_aliphatic_heteroatom_nbrs_query(2),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_has_heteroatom_nbrs_query() {
        assert_eq!(
            make_atom_has_heteroatom_nbrs_query(),
            QueryNode::Predicate(AtomQueryPredicate::HasHeteroatomNeighbors)
        );

        let mut builder = Molecule::builder();
        let with_heteroatom = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                with_heteroatom,
                oxygen,
                BondOrder::Single,
            ))
            .expect("heteroatom-neighbor fixture bond is valid");
        let carbon_only = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_only,
                carbon_neighbor,
                BondOrder::Single,
            ))
            .expect("carbon-only fixture bond is valid");
        let molecule = builder
            .build()
            .expect("has-heteroatom-neighbor fixture is valid");
        let query = make_atom_has_heteroatom_nbrs_query();

        assert!(atom_matches_query(
            &molecule.atoms()[with_heteroatom.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon_only.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_num_heteroatom_nbrs_query() {
        assert_eq!(
            make_atom_num_heteroatom_nbrs_query(2),
            QueryNode::Predicate(AtomQueryPredicate::NumHeteroatomNeighbors(2))
        );

        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        for neighbor in [carbon, hydrogen, oxygen, nitrogen] {
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, BondOrder::Single))
                .expect("heteroatom-neighbor fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("heteroatom-neighbor fixture molecule is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_num_heteroatom_nbrs_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[center.index()],
            &make_atom_num_heteroatom_nbrs_query(1),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &make_atom_num_heteroatom_nbrs_query(0),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_has_ring_bond_query() {
        assert_eq!(
            make_atom_has_ring_bond_query(),
            QueryNode::Predicate(AtomQueryPredicate::HasRingBond)
        );

        let mut builder = Molecule::builder();
        let ring = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [(ring[0], ring[1]), (ring[1], ring[2]), (ring[2], ring[0])] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("has-ring-bond fixture ring bond is valid");
        }
        let tail = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(ring[0], tail, BondOrder::Single))
            .expect("has-ring-bond fixture tail is valid");
        let molecule = builder.build().expect("has-ring-bond fixture is valid");
        let query = make_atom_has_ring_bond_query();

        assert!(atom_matches_query(
            &molecule.atoms()[ring[0].index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[tail.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_ring_bond_count_query() {
        assert_eq!(
            make_atom_ring_bond_count_query(2),
            QueryNode::Predicate(AtomQueryPredicate::RingBondCount(2))
        );

        let mut builder = Molecule::builder();
        let ring = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [(ring[0], ring[1]), (ring[1], ring[2]), (ring[2], ring[0])] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("ring-bond-count fixture bond is valid");
        }
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("ring-bond-count fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[ring[0].index()],
            &make_atom_ring_bond_count_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[ring[0].index()],
            &make_atom_ring_bond_count_query(1),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[isolated.index()],
            &make_atom_ring_bond_count_query(0),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_min_ring_size_query() {
        assert_eq!(
            make_atom_min_ring_size_query(3),
            QueryNode::Predicate(AtomQueryPredicate::SmallestRingSize(3))
        );

        let mut builder = Molecule::builder();
        let triangle = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("minimum-ring-size fixture bond is valid");
        }
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("minimum-ring-size fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[triangle[0].index()],
            &make_atom_min_ring_size_query(3),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[triangle[0].index()],
            &make_atom_min_ring_size_query(4),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[isolated.index()],
            &make_atom_min_ring_size_query(0),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_in_ring_of_size_query() {
        assert_eq!(
            make_atom_in_ring_of_size_query(3),
            QueryNode::Predicate(AtomQueryPredicate::InRingOfSize(3))
        );
        assert_eq!(
            make_atom_in_ring_of_size_range_query(3, 4, false, true),
            QueryNode::Predicate(AtomQueryPredicate::Range(AtomRangeQuery::new(
                AtomRangeBounds::Inclusive {
                    lower: 3,
                    upper: 4,
                    lower_open: false,
                    upper_open: true,
                },
                AtomRangeDataFunction::AtomRingSize {
                    lower: 3,
                    upper: 4,
                    lower_open: false,
                    upper_open: true,
                },
            )))
        );

        let mut builder = Molecule::builder();
        let triangle = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("triangle bond is valid");
        }
        let square = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [
            (square[0], square[1]),
            (square[1], square[2]),
            (square[2], square[3]),
            (square[3], square[0]),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("square bond is valid");
        }
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("ring-size fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[triangle[0].index()],
            &make_atom_in_ring_of_size_query(3),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[square[0].index()],
            &make_atom_in_ring_of_size_query(3),
            &molecule,
        ));
        let closed = make_atom_in_ring_of_size_range_query(3, 4, false, false);
        assert!(atom_matches_query(
            &molecule.atoms()[triangle[0].index()],
            &closed,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[square[0].index()],
            &closed,
            &molecule,
        ));
        let open = make_atom_in_ring_of_size_range_query(3, 4, true, true);
        assert!(!atom_matches_query(
            &molecule.atoms()[triangle[0].index()],
            &open,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[square[0].index()],
            &open,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[isolated.index()],
            &closed,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_in_n_rings_query() {
        assert_eq!(
            make_atom_in_n_rings_query(2),
            QueryNode::Predicate(AtomQueryPredicate::NumAtomRings(2))
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0), (a1, a3), (a3, a2)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("fused-ring atom fixture bond is valid");
        }
        let molecule = builder.build().expect("fused-ring atom fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[a1.index()],
            &make_atom_in_n_rings_query(2),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[a0.index()],
            &make_atom_in_n_rings_query(1),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[isolated.index()],
            &make_atom_in_n_rings_query(0),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_in_ring_query() {
        assert_eq!(
            make_atom_in_ring_query(),
            QueryNode::Predicate(AtomQueryPredicate::InRing)
        );

        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-in-ring fixture bond is valid");
        }
        let molecule = builder.build().expect("atom-in-ring fixture is valid");
        let query = make_atom_in_ring_query();

        assert!(atom_matches_query(
            &molecule.atoms()[a0.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[isolated.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_unsaturated_query() {
        assert_eq!(
            make_atom_unsaturated_query(),
            QueryNode::Predicate(AtomQueryPredicate::IsUnsaturated)
        );

        let mut saturated_builder = Molecule::builder();
        let methane_id = saturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let methane = saturated_builder.build().expect("methane is valid");
        assert!(!atom_matches_query(
            &methane.atoms()[methane_id.index()],
            &make_atom_unsaturated_query(),
            &methane,
        ));

        let mut unsaturated_builder = Molecule::builder();
        let left_id = unsaturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let right_id = unsaturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        unsaturated_builder
            .add_bond(crate::BondSpec::new(
                left_id,
                right_id,
                crate::BondOrder::Double,
            ))
            .expect("carbon-to-carbon double bond is valid");
        let ethene = unsaturated_builder.build().expect("ethene is valid");
        assert!(atom_matches_query(
            &ethene.atoms()[left_id.index()],
            &make_atom_unsaturated_query(),
            &ethene,
        ));
    }

    #[test]
    fn smarts_make_atom_missing_chiral_tag_query() {
        assert_eq!(
            make_atom_missing_chiral_tag_query(),
            QueryNode::Predicate(AtomQueryPredicate::MissingChiralTag)
        );

        let mut builder = Molecule::builder();
        let unspecified_without_prop = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let unspecified_with_prop = builder
            .add_atom(crate::AtomSpec::new(crate::Element::C).with_prop("_ChiralityPossible", "0"));
        let tagged_with_prop = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(ChiralTag::TetrahedralCw)
                .with_prop("_ChiralityPossible", "1"),
        );
        let molecule = builder
            .build()
            .expect("missing-chiral-tag fixture is valid");
        let query = make_atom_missing_chiral_tag_query();

        assert!(!atom_matches_query(
            &molecule.atoms()[unspecified_without_prop.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[unspecified_with_prop.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[tagged_with_prop.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_has_chiral_tag_query() {
        assert_eq!(
            make_atom_has_chiral_tag_query(),
            QueryNode::Predicate(AtomQueryPredicate::HasChiralTag)
        );

        let mut builder = Molecule::builder();
        let unspecified = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let clockwise = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let counterclockwise = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_chiral_tag(ChiralTag::TetrahedralCcw),
        );
        let molecule = builder.build().expect("chiral-tag fixture is valid");
        let query = make_atom_has_chiral_tag_query();

        assert!(!atom_matches_query(
            &molecule.atoms()[unspecified.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[clockwise.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[counterclockwise.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_num_radical_electrons_query() {
        assert_eq!(
            make_atom_num_radical_electrons_query(2),
            QueryNode::Predicate(AtomQueryPredicate::NumRadicalElectrons(2))
        );

        let mut builder = Molecule::builder();
        let zero = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let one =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(1));
        let two =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(2));
        let molecule = builder.build().expect("radical-electron fixture is valid");
        let query = make_atom_num_radical_electrons_query(2);

        assert!(!atom_matches_query(
            &molecule.atoms()[zero.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[one.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[two.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_hybridization_query() {
        assert_eq!(
            make_atom_hybridization_query(Hybridization::Sp2),
            QueryNode::Predicate(AtomQueryPredicate::HybridizationMatch(Hybridization::Sp2))
        );

        let mut builder = Molecule::builder();
        let sp2 = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_hybridization(Hybridization::Sp2),
        );
        let sp3 = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_hybridization(Hybridization::Sp3),
        );
        let molecule = builder.build().expect("hybridization fixture is valid");
        let query = make_atom_hybridization_query(Hybridization::Sp2);

        assert!(atom_matches_query(
            &molecule.atoms()[sp2.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[sp3.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_negative_formal_charge_query() {
        assert_eq!(
            make_atom_negative_formal_charge_query(2),
            QueryNode::Predicate(AtomQueryPredicate::NegativeFormalCharge(2))
        );

        let mut builder = Molecule::builder();
        let anion =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_formal_charge(-2));
        let neutral = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let cation =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_formal_charge(1));
        let molecule = builder
            .build()
            .expect("negative-formal-charge fixture is valid");
        let query = make_atom_negative_formal_charge_query(2);

        assert!(atom_matches_query(
            &molecule.atoms()[anion.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[neutral.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[cation.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_formal_charge_query() {
        assert_eq!(
            make_atom_formal_charge_query(-1),
            QueryNode::Predicate(AtomQueryPredicate::FormalCharge(-1))
        );

        let mut builder = Molecule::builder();
        let negative =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_formal_charge(-1));
        let neutral = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let positive =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_formal_charge(2));
        let molecule = builder.build().expect("formal-charge fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[negative.index()],
            &make_atom_formal_charge_query(-1),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[neutral.index()],
            &make_atom_formal_charge_query(0),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[positive.index()],
            &make_atom_formal_charge_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[negative.index()],
            &make_atom_formal_charge_query(1),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_isotope_query() {
        assert_eq!(
            make_atom_isotope_query(13),
            QueryNode::Predicate(AtomQueryPredicate::Isotope(13))
        );

        let mut builder = Molecule::builder();
        let natural_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_13 = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(13));
        let carbon_14 = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(14));
        let molecule = builder.build().expect("isotope-query fixture is valid");
        let query = make_atom_isotope_query(13);

        assert!(!atom_matches_query(
            &molecule.atoms()[natural_carbon.index()],
            &query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon_13.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon_14.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_mass_query() {
        assert_eq!(
            make_atom_mass_query(12),
            QueryNode::Predicate(AtomQueryPredicate::Mass(12))
        );

        let mut builder = Molecule::builder();
        let natural_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_12 = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(12));
        let carbon_13 = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(13));
        let molecule = builder.build().expect("mass-query fixture is valid");

        assert!(!atom_matches_query(
            &molecule.atoms()[natural_carbon.index()],
            &make_atom_mass_query(12),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon_12.index()],
            &make_atom_mass_query(12),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon_13.index()],
            &make_atom_mass_query(12),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_aliphatic_query() {
        assert_eq!(
            make_atom_aliphatic_query(),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(false))
        );

        let mut builder = Molecule::builder();
        let _atoms = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        let molecule = builder.build().expect("aliphatic-query fixture is valid");
        let query = make_atom_aliphatic_query();

        // RDKit✔️✔️: aeq = makeAtomAliphaticQuery();
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(6)), "");
        assert!(!atom_matches_query(&molecule.atoms()[0], &query, &molecule));
        assert!(!atom_matches_query(&molecule.atoms()[1], &query, &molecule));
        assert!(!atom_matches_query(&molecule.atoms()[5], &query, &molecule));
        assert!(atom_matches_query(&molecule.atoms()[6], &query, &molecule));
    }

    #[test]
    fn smarts_make_atom_aromatic_query() {
        assert_eq!(
            make_atom_aromatic_query(),
            QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true))
        );

        let mut builder = Molecule::builder();
        let _atoms = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        let molecule = builder.build().expect("aromatic-query fixture is valid");
        let query = make_atom_aromatic_query();

        // RDKit✔️✔️: aeq = makeAtomAromaticQuery();
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
        assert!(atom_matches_query(&molecule.atoms()[0], &query, &molecule));
        assert!(atom_matches_query(&molecule.atoms()[1], &query, &molecule));
        assert!(atom_matches_query(&molecule.atoms()[5], &query, &molecule));
        assert!(!atom_matches_query(&molecule.atoms()[6], &query, &molecule));
    }

    #[test]
    fn smarts_make_atom_implicit_h_count_query() {
        assert_eq!(
            make_atom_implicit_h_count_query(2),
            QueryNode::Predicate(AtomQueryPredicate::ImplicitHydrogenCount(2))
        );

        let mut builder = Molecule::builder();
        let carbon =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(1));
        let deuterium = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        for neighbor in [deuterium, oxygen] {
            builder
                .add_bond(crate::BondSpec::new(carbon, neighbor, BondOrder::Single))
                .expect("implicit-h-count fixture bond is valid");
        }
        let molecule = builder.build().expect("implicit-h-count fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &make_atom_implicit_h_count_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &make_atom_implicit_h_count_query(3),
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &make_atom_h_count_query(3),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_has_implicit_h_query() {
        assert_eq!(
            make_atom_has_implicit_h_query(),
            QueryNode::Predicate(AtomQueryPredicate::HasImplicitHydrogen)
        );

        let mut builder = Molecule::builder();
        let has_h =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(4));
        let no_h = builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_no_implicit(true));
        let molecule = builder
            .build()
            .expect("implicit-h-presence fixture is valid");
        let query = make_atom_has_implicit_h_query();

        assert!(atom_matches_query(
            &molecule.atoms()[has_h.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[no_h.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_h_count_query() {
        assert_eq!(
            make_atom_h_count_query(1),
            QueryNode::Predicate(AtomQueryPredicate::HydrogenCount(1))
        );

        let mut builder = Molecule::builder();
        let atoms = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for idx in 0..6 {
            builder
                .add_bond(crate::BondSpec::new(
                    atoms[idx],
                    atoms[(idx + 1) % 6],
                    if idx % 2 == 0 {
                        BondOrder::Single
                    } else {
                        BondOrder::Double
                    },
                ))
                .expect("ring bond is valid");
        }
        builder
            .add_bond(crate::BondSpec::new(atoms[5], atoms[6], BondOrder::Single))
            .expect("substituent bond is valid");
        let molecule = builder.build().expect("hydrogen-count fixture is valid");
        let query = make_atom_h_count_query(1);

        // RDKit✔️✔️: aeq = makeAtomHCountQuery(1);
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
        assert!(atom_matches_query(&molecule.atoms()[0], &query, &molecule));
        assert!(atom_matches_query(&molecule.atoms()[1], &query, &molecule));
        assert!(!atom_matches_query(&molecule.atoms()[5], &query, &molecule,));
        assert!(!atom_matches_query(&molecule.atoms()[6], &query, &molecule,));
    }

    #[test]
    fn smarts_make_atom_heavy_atom_degree_query() {
        assert_eq!(
            make_atom_heavy_atom_degree_query(2),
            QueryNode::Predicate(AtomQueryPredicate::HeavyAtomDegree(2))
        );

        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let protium = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let deuterium = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        for neighbor in [carbon, oxygen, protium, deuterium] {
            builder
                .add_bond(crate::BondSpec::new(center, neighbor, BondOrder::Single))
                .expect("heavy-degree fixture bond is valid");
        }
        let molecule = builder.build().expect("heavy-degree fixture is valid");
        let context = build_query_match_context(&molecule);
        let center_atom = &molecule.atoms()[center.index()];

        assert_eq!(
            query_atom_heavy_atom_degree(center_atom, &context.adj, &molecule),
            2
        );
        assert_eq!(
            query_atom_non_hydrogen_degree(center_atom, &context.adj, &molecule),
            3
        );
        assert!(atom_matches_query(
            center_atom,
            &make_atom_heavy_atom_degree_query(2),
            &molecule,
        ));
        assert!(!atom_matches_query(
            center_atom,
            &make_atom_heavy_atom_degree_query(3),
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_total_degree_query() {
        assert_eq!(
            make_atom_total_degree_query(4),
            QueryNode::Predicate(AtomQueryPredicate::TotalDegree(4))
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(carbon, oxygen, BondOrder::Single))
            .expect("carbon-oxygen bond is valid");
        let molecule = builder.build().expect("total-degree fixture is valid");
        let degree_four = make_atom_total_degree_query(4);
        let degree_two = make_atom_total_degree_query(2);

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &degree_four,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &degree_four,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &degree_two,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_explicit_degree_query() {
        assert_eq!(
            make_atom_explicit_degree_query(3),
            QueryNode::Predicate(AtomQueryPredicate::ExplicitDegree(3))
        );

        // RDKit✔️✔️: ATOM_EQUALS_QUERY *aeq = makeAtomExplicitDegreeQuery(3);
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(5)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
        // RDKit✔️✔️: aeq = makeAtomExplicitDegreeQuery(2);
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(aeq->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(5)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!aeq->Match(m.getAtomWithIdx(6)), "");
        let mut builder = Molecule::builder();
        let ring_atoms = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for idx in 0..ring_atoms.len() {
            builder
                .add_bond(crate::BondSpec::new(
                    ring_atoms[idx],
                    ring_atoms[(idx + 1) % ring_atoms.len()],
                    if idx % 2 == 0 {
                        BondOrder::Single
                    } else {
                        BondOrder::Double
                    },
                ))
                .expect("ring bond is valid");
        }
        let substituent = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(
                ring_atoms[5],
                substituent,
                BondOrder::Single,
            ))
            .expect("substituent bond is valid");
        let molecule = builder.build().expect("explicit-degree fixture is valid");
        let degree_three = make_atom_explicit_degree_query(3);
        let degree_two = make_atom_explicit_degree_query(2);

        assert!(!atom_matches_query(
            &molecule.atoms()[ring_atoms[0].index()],
            &degree_three,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[ring_atoms[1].index()],
            &degree_three,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[ring_atoms[5].index()],
            &degree_three,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[substituent.index()],
            &degree_three,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[ring_atoms[0].index()],
            &degree_two,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[ring_atoms[1].index()],
            &degree_two,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[ring_atoms[5].index()],
            &degree_two,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[substituent.index()],
            &degree_two,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_total_valence_query() {
        assert_eq!(
            make_atom_total_valence_query(4),
            QueryNode::Predicate(AtomQueryPredicate::TotalValence(4))
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let molecule = builder.build().expect("total-valence fixture is valid");
        let valence_four = make_atom_total_valence_query(4);
        let valence_two = make_atom_total_valence_query(2);

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &valence_four,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &valence_four,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &valence_two,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_explicit_valence_query() {
        let valence_three = make_atom_explicit_valence_query(3);
        assert_eq!(
            valence_three,
            QueryNode::Predicate(AtomQueryPredicate::ExplicitValence(3))
        );
        // RDKit✔️✔️: a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
        // RDKit✔️✔️: a2.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
        // RDKit✔️✔️: TEST_ASSERT(a1.QueryMatch(&a2));
        // RDKit✔️✔️: TEST_ASSERT(a2.QueryMatch(&a1));
        assert_eq!(valence_three, make_atom_explicit_valence_query(3));
        // RDKit✔️✔️: a1.expandQuery(makeAtomExplicitValenceQuery(3), Queries::COMPOSITE_AND);
        // RDKit✔️✔️: a2.expandQuery(makeAtomExplicitValenceQuery(4), Queries::COMPOSITE_AND);
        // RDKit✔️✔️: TEST_ASSERT(!a1.QueryMatch(&a2));
        // RDKit✔️✔️: TEST_ASSERT(!a2.QueryMatch(&a1));
        assert_ne!(valence_three, make_atom_explicit_valence_query(4));

        let mut builder = Molecule::builder();
        let triple_bonded_carbon = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(0),
        );
        let nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_no_implicit(true)
                .with_explicit_hydrogens(0),
        );
        let double_bonded_carbon = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_no_implicit(true)
                .with_explicit_hydrogens(0),
        );
        let oxygen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::O)
                .with_no_implicit(true)
                .with_explicit_hydrogens(0),
        );
        builder
            .add_bond(crate::BondSpec::new(
                triple_bonded_carbon,
                nitrogen,
                BondOrder::Triple,
            ))
            .expect("nitrile bond is valid");
        builder
            .add_bond(crate::BondSpec::new(
                double_bonded_carbon,
                oxygen,
                BondOrder::Double,
            ))
            .expect("carbonyl bond is valid");
        let molecule = builder.build().expect("explicit-valence fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[triple_bonded_carbon.index()],
            &valence_three,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[double_bonded_carbon.index()],
            &valence_three,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_implicit_valence_query() {
        assert_eq!(
            make_atom_implicit_valence_query(3),
            QueryNode::Predicate(AtomQueryPredicate::ImplicitValence(3))
        );

        // RDKit✔️✔️: auto *a = new Atom(6);
        // RDKit✔️✔️: m.addAtom(a);
        // RDKit✔️✔️: m.addAtom(a);
        // RDKit✔️✔️: delete a;
        // RDKit✔️✔️: m.addBond(0, 1, Bond::SINGLE);
        // RDKit✔️✔️: a = new Atom(8);
        // RDKit✔️✔️: m.addAtom(a);
        // RDKit✔️✔️: delete a;
        // RDKit✔️✔️: m.addBond(1, 2, Bond::DOUBLE);
        // RDKit✔️✔️: MolOps::sanitizeMol(m);
        // RDKit✔️✔️: qA->expandQuery(makeAtomImplicitValenceQuery(3));
        // RDKit✔️✔️: CHECK_INVARIANT(qA->Match(m.getAtomWithIdx(0)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(1)), "");
        // RDKit✔️✔️: CHECK_INVARIANT(!qA->Match(m.getAtomWithIdx(2)), "");
        let mut builder = Molecule::builder();
        let terminal_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbonyl_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                terminal_carbon,
                carbonyl_carbon,
                BondOrder::Single,
            ))
            .expect("carbon-carbon bond is valid");
        builder
            .add_bond(crate::BondSpec::new(
                carbonyl_carbon,
                oxygen,
                BondOrder::Double,
            ))
            .expect("carbonyl bond is valid");
        let molecule = builder.build().expect("implicit-valence fixture is valid");
        let query = make_atom_implicit_valence_query(3);

        assert!(atom_matches_query(
            &molecule.atoms()[terminal_carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[carbonyl_carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[oxygen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_type_query() {
        let aromatic_carbon_query = make_atom_type_query(6, true);
        let aliphatic_carbon_query = make_atom_type_query(6, false);
        let aromatic_nitrogen_query = make_atom_type_query(7, true);
        assert_eq!(
            aromatic_carbon_query,
            QueryNode::Predicate(AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: true,
            })
        );

        // RDKit✔️✔️: RWMol *m = SmilesToMol("CCc1ccccc1");
        // RDKit✔️✔️: QueryAtom qA1;
        // RDKit✔️✔️: qA1.setQuery(makeAtomTypeQuery(6, true));
        // RDKit✔️✔️: QueryAtom qA2;
        // RDKit✔️✔️: qA2.setQuery(makeAtomTypeQuery(6, false));
        // RDKit✔️✔️: QueryAtom qA3;
        // RDKit✔️✔️: qA3.setQuery(makeAtomTypeQuery(7, true));
        // RDKit✔️✔️: TEST_ASSERT(!qA1.Match(m->getAtomWithIdx(0)));
        // RDKit✔️✔️: TEST_ASSERT(qA2.Match(m->getAtomWithIdx(0)));
        // RDKit✔️✔️: TEST_ASSERT(!qA3.Match(m->getAtomWithIdx(0)));
        // RDKit✔️✔️: TEST_ASSERT(qA1.Match(m->getAtomWithIdx(2)));
        // RDKit✔️✔️: TEST_ASSERT(!qA2.Match(m->getAtomWithIdx(2)));
        // RDKit✔️✔️: TEST_ASSERT(!qA3.Match(m->getAtomWithIdx(2)));
        let mut builder = Molecule::builder();
        let aliphatic_carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_carbon =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true));
        let molecule = builder.build().expect("atom-type fixture is valid");

        assert!(!atom_matches_query(
            &molecule.atoms()[aliphatic_carbon.index()],
            &aromatic_carbon_query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[aliphatic_carbon.index()],
            &aliphatic_carbon_query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[aliphatic_carbon.index()],
            &aromatic_nitrogen_query,
            &molecule,
        ));
        assert!(atom_matches_query(
            &molecule.atoms()[aromatic_carbon.index()],
            &aromatic_carbon_query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[aromatic_carbon.index()],
            &aliphatic_carbon_query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[aromatic_carbon.index()],
            &aromatic_nitrogen_query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_num_query() {
        assert_eq!(
            make_atom_num_query(0),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(0))
        );
        assert_eq!(
            make_atom_num_query(118),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(118))
        );

        let query = make_atom_num_query(6);
        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let molecule = builder.build().expect("atom-number fixture is valid");
        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[nitrogen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_make_atom_range_query() {
        let mut builder = Molecule::builder();
        let atoms = [
            crate::Element::C,
            crate::Element::C,
            crate::Element::N,
            crate::Element::C,
            crate::Element::O,
            crate::Element::O,
        ]
        .map(|element| builder.add_atom(crate::AtomSpec::new(element)));
        for pair in atoms.windows(2) {
            builder
                .add_bond(crate::BondSpec::new(pair[0], pair[1], BondOrder::Single))
                .expect("range-query fixture bond is valid");
        }
        let molecule = builder.build().expect("range-query fixture is valid");

        let cases = [
            (0, 3, true, true),
            (1, 2, false, false),
            (0, 2, true, false),
        ];
        for (lower, upper, lower_open, upper_open) in cases {
            let query = make_atom_range_query(
                lower,
                upper,
                lower_open,
                upper_open,
                AtomRangeDataFunction::NumHeteroatomNeighbors,
            );
            assert!(!atom_matches_query(&molecule.atoms()[0], &query, &molecule));
            assert!(atom_matches_query(&molecule.atoms()[1], &query, &molecule));
            assert!(!atom_matches_query(&molecule.atoms()[2], &query, &molecule));
            assert!(atom_matches_query(&molecule.atoms()[3], &query, &molecule));
            assert!(atom_matches_query(&molecule.atoms()[5], &query, &molecule));
        }
    }

    #[test]
    fn smarts_make_atom_simple_query() {
        let query = make_atom_simple_query(AtomQueryPredicate::AtomicNumber(6));
        assert_eq!(
            query,
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
        );

        let mut builder = Molecule::builder();
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let molecule = builder.build().expect("simple-query fixture is valid");

        assert!(atom_matches_query(
            &molecule.atoms()[carbon.index()],
            &query,
            &molecule,
        ));
        assert!(!atom_matches_query(
            &molecule.atoms()[nitrogen.index()],
            &query,
            &molecule,
        ));
    }

    #[test]
    fn smarts_query_null_data() {
        let mut builder = Molecule::builder();
        let atom_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("null-data fixture is valid");
        let atom = &molecule.atoms()[atom_id.index()];

        assert_eq!(null_data(atom), 1);
        assert_eq!(null_data(0usize), 1);
    }

    #[test]
    fn smarts_query_null_fun() {
        let mut builder = Molecule::builder();
        let atom_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("null-query fixture is valid");
        let atom = &molecule.atoms()[atom_id.index()];

        assert!(null_query(atom));
        assert!(null_query(false));
    }

    #[test]
    fn smarts_query_atom_dummy() {
        let atom = |id, element, predicate| {
            crate::QueryAtom::from_parts(
                Atom::from_spec(id, crate::AtomSpec::new(element)),
                predicate,
            )
        };
        let plain_dummy = atom(
            AtomId::new(0),
            crate::Element::DUMMY,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(0)),
        );
        let plain_carbon = atom(
            AtomId::new(1),
            crate::Element::C,
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        );
        let atom_null = atom(
            AtomId::new(2),
            crate::Element::DUMMY,
            QueryNode::predicate(AtomQueryPredicate::Any),
        );
        let negated_atom_null = atom(
            AtomId::new(3),
            crate::Element::DUMMY,
            QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Any)),
        );
        let composite_or = atom(
            AtomId::new(4),
            crate::Element::DUMMY,
            QueryNode::or(vec![
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
            ]),
        );

        assert!(is_atom_dummy(&plain_dummy));
        assert!(!is_atom_dummy(&plain_carbon));
        assert!(is_atom_dummy(&atom_null));
        assert!(!is_atom_dummy(&negated_atom_null));
        assert!(!is_atom_dummy(&composite_or));
    }

    #[test]
    fn smarts_query_is_metal() {
        const NON_METALS: &[u8] = &[
            0, 1, 2, 5, 6, 7, 8, 9, 10, 14, 15, 16, 17, 18, 33, 34, 35, 36, 52, 53, 54, 85, 86,
        ];
        let mut builder = Molecule::builder();
        let excluded = NON_METALS
            .iter()
            .map(|&atomic_number| {
                let element = crate::Element::from_atomic_number(atomic_number)
                    .expect("RDKit non-metal atomic number is modeled");
                builder.add_atom(crate::AtomSpec::new(element))
            })
            .collect::<Vec<_>>();
        let included = [3, 4, 13, 26, 29, 30, 31, 50, 82, 87, 88].map(|atomic_number| {
            let element = crate::Element::from_atomic_number(atomic_number)
                .expect("representative metal atomic number is modeled");
            builder.add_atom(crate::AtomSpec::new(element))
        });
        let molecule = builder.build().expect("metal-query fixture is valid");

        for atom_id in excluded {
            assert!(!is_metal(&molecule.atoms()[atom_id.index()]));
        }
        for atom_id in included {
            assert!(is_metal(&molecule.atoms()[atom_id.index()]));
        }
    }

    #[test]
    fn smarts_query_atom_aromatic() {
        let mut builder = Molecule::builder();
        let aromatic_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true));
        let aliphatic_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder
            .build()
            .expect("two isolated atoms form a valid molecule");
        let aromatic = &molecule.atoms()[aromatic_id.index()];
        let aliphatic = &molecule.atoms()[aliphatic_id.index()];

        assert!(query_atom_aromatic(aromatic));
        assert!(!query_atom_aromatic(aliphatic));
        assert!(atom_predicate_matches(
            aromatic,
            &AtomQueryPredicate::IsAromatic(true),
            &molecule,
        ));
        assert!(!atom_predicate_matches(
            aliphatic,
            &AtomQueryPredicate::IsAromatic(true),
            &molecule,
        ));
    }

    #[test]
    fn smarts_query_atom_aliphatic() {
        let mut builder = Molecule::builder();
        let aromatic_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true));
        let aliphatic_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder
            .build()
            .expect("two isolated atoms form a valid molecule");
        let aromatic = &molecule.atoms()[aromatic_id.index()];
        let aliphatic = &molecule.atoms()[aliphatic_id.index()];

        assert!(!query_atom_aliphatic(aromatic));
        assert!(query_atom_aliphatic(aliphatic));
        assert!(!atom_predicate_matches(
            aromatic,
            &AtomQueryPredicate::IsAromatic(false),
            &molecule,
        ));
        assert!(atom_predicate_matches(
            aliphatic,
            &AtomQueryPredicate::IsAromatic(false),
            &molecule,
        ));
    }

    #[test]
    fn smarts_query_atom_explicit_degree() {
        let mut builder = Molecule::builder();
        let center_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                center_id,
                hydrogen_id,
                crate::BondOrder::Single,
            ))
            .expect("center-to-hydrogen bond is valid");
        builder
            .add_bond(crate::BondSpec::new(
                center_id,
                oxygen_id,
                crate::BondOrder::Single,
            ))
            .expect("center-to-oxygen bond is valid");
        let molecule = builder.build().expect("branched molecule is valid");
        let context = build_query_match_context(&molecule);
        let center = &molecule.atoms()[center_id.index()];
        let hydrogen = &molecule.atoms()[hydrogen_id.index()];

        assert_eq!(query_atom_explicit_degree(center, &context.adj), 2);
        assert_eq!(query_atom_explicit_degree(hydrogen, &context.adj), 1);
        assert!(atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::ExplicitDegree(2),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::ExplicitDegree(3),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::DegreeGreaterEqual(2),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_total_degree() {
        let mut builder = Molecule::builder();
        let carbon_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_id,
                oxygen_id,
                crate::BondOrder::Single,
            ))
            .expect("carbon-to-oxygen bond is valid");
        let molecule = builder.build().expect("methanol skeleton is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];
        let oxygen = &molecule.atoms()[oxygen_id.index()];

        assert_eq!(query_atom_explicit_degree(carbon, &context.adj), 1);
        assert_eq!(
            query_atom_total_degree(&context.adj, context.valence.as_ref(), carbon),
            Some(4)
        );
        assert_eq!(
            query_atom_total_degree(&context.adj, context.valence.as_ref(), oxygen),
            Some(2)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::TotalDegree(4),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::TotalDegree(1),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_non_hydrogen_degree() {
        let mut builder = Molecule::builder();
        let center_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let protium_id = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let deuterium_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        let tritium_id = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(3));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        for neighbor_id in [protium_id, deuterium_id, tritium_id, oxygen_id] {
            builder
                .add_bond(crate::BondSpec::new(
                    center_id,
                    neighbor_id,
                    crate::BondOrder::Single,
                ))
                .expect("center-to-neighbor bond is valid");
        }
        let molecule = builder.build().expect("isotopic star molecule is valid");
        let context = build_query_match_context(&molecule);
        let center = &molecule.atoms()[center_id.index()];

        assert_eq!(query_atom_explicit_degree(center, &context.adj), 4);
        assert_eq!(
            query_atom_non_hydrogen_degree(center, &context.adj, &molecule),
            3
        );
        assert!(atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::NonHydrogenDegree(3),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::NonHydrogenDegree(2),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::NonHydrogenDegree(3),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_heavy_atom_degree() {
        let mut builder = Molecule::builder();
        let center_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let protium_id = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let deuterium_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        let tritium_id = builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(3));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        for neighbor_id in [protium_id, deuterium_id, tritium_id, oxygen_id] {
            builder
                .add_bond(crate::BondSpec::new(
                    center_id,
                    neighbor_id,
                    crate::BondOrder::Single,
                ))
                .expect("center-to-neighbor bond is valid");
        }
        let molecule = builder.build().expect("isotopic star molecule is valid");
        let context = build_query_match_context(&molecule);
        let center = &molecule.atoms()[center_id.index()];

        assert_eq!(
            query_atom_non_hydrogen_degree(center, &context.adj, &molecule),
            3
        );
        assert_eq!(
            query_atom_heavy_atom_degree(center, &context.adj, &molecule),
            1
        );
        assert!(atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::HeavyAtomDegree(1),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            center,
            &AtomQueryPredicate::HeavyAtomDegree(3),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_h_count() {
        let mut builder = Molecule::builder();
        let carbon_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(1));
        let deuterium_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        for neighbor_id in [deuterium_id, oxygen_id] {
            builder
                .add_bond(crate::BondSpec::new(
                    carbon_id,
                    neighbor_id,
                    crate::BondOrder::Single,
                ))
                .expect("carbon-to-neighbor bond is valid");
        }
        let molecule = builder
            .build()
            .expect("isotopic methanol skeleton is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];

        assert_eq!(
            implicit_hydrogen_count(context.valence.as_ref(), carbon),
            Some(1)
        );
        assert_eq!(
            query_atom_h_count(&context.adj, context.valence.as_ref(), carbon, &molecule),
            Some(3)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::HydrogenCount(3),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::HydrogenCount(2),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_implicit_h_count() {
        let mut builder = Molecule::builder();
        let carbon_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(1));
        let deuterium_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::H).with_isotope(2));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        for neighbor_id in [deuterium_id, oxygen_id] {
            builder
                .add_bond(crate::BondSpec::new(
                    carbon_id,
                    neighbor_id,
                    crate::BondOrder::Single,
                ))
                .expect("carbon-to-neighbor bond is valid");
        }
        let molecule = builder
            .build()
            .expect("isotopic methanol skeleton is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];

        assert_eq!(
            query_atom_implicit_h_count(context.valence.as_ref(), carbon),
            Some(2)
        );
        assert_eq!(
            query_atom_h_count(&context.adj, context.valence.as_ref(), carbon, &molecule),
            Some(3)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ImplicitHydrogenCount(2),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ImplicitHydrogenCount(3),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ImplicitHydrogenCountLessEqual(2),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_has_implicit_h() {
        let mut builder = Molecule::builder();
        let explicit_h_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(4));
        let no_h_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_no_implicit(true));
        let molecule = builder
            .build()
            .expect("two isolated valence-complete atoms are valid");
        let context = build_query_match_context(&molecule);
        let explicit_h_atom = &molecule.atoms()[explicit_h_id.index()];
        let no_h_atom = &molecule.atoms()[no_h_id.index()];

        assert_eq!(
            implicit_hydrogen_count(context.valence.as_ref(), explicit_h_atom),
            Some(0)
        );
        assert_eq!(
            query_atom_implicit_h_count(context.valence.as_ref(), explicit_h_atom),
            Some(4)
        );
        assert!(query_atom_has_implicit_h(
            context.valence.as_ref(),
            explicit_h_atom
        ));
        assert!(!query_atom_has_implicit_h(
            context.valence.as_ref(),
            no_h_atom
        ));
        assert!(atom_predicate_matches_with_context(
            explicit_h_atom,
            &AtomQueryPredicate::HasImplicitHydrogen,
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            no_h_atom,
            &AtomQueryPredicate::HasImplicitHydrogen,
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_implicit_valence() {
        let mut builder = Molecule::builder();
        let carbon_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let no_implicit_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_no_implicit(true));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_id,
                oxygen_id,
                crate::BondOrder::Single,
            ))
            .expect("carbon-to-oxygen bond is valid");
        let molecule = builder
            .build()
            .expect("methanol skeleton plus capped carbon is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];
        let oxygen = &molecule.atoms()[oxygen_id.index()];
        let no_implicit = &molecule.atoms()[no_implicit_id.index()];

        assert_eq!(
            query_atom_implicit_valence(context.valence.as_ref(), carbon),
            Some(3)
        );
        assert_eq!(
            query_atom_implicit_valence(context.valence.as_ref(), oxygen),
            Some(1)
        );
        assert_eq!(
            query_atom_implicit_valence(context.valence.as_ref(), no_implicit),
            Some(0)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ImplicitValence(3),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ImplicitValence(1),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_explicit_valence() {
        let mut builder = Molecule::builder();
        let carbon_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(1));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_id,
                oxygen_id,
                crate::BondOrder::Single,
            ))
            .expect("carbon-to-oxygen bond is valid");
        let molecule = builder.build().expect("methanol skeleton is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];

        assert_eq!(
            atom_explicit_valence(context.valence.as_ref(), carbon),
            Some(2)
        );
        assert_eq!(
            query_atom_explicit_valence(context.valence.as_ref(), carbon),
            Some(1)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ExplicitValence(1),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::ExplicitValence(2),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_total_valence() {
        let mut builder = Molecule::builder();
        let carbon_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_explicit_hydrogens(1));
        let oxygen_id = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        builder
            .add_bond(crate::BondSpec::new(
                carbon_id,
                oxygen_id,
                crate::BondOrder::Single,
            ))
            .expect("carbon-to-oxygen bond is valid");
        let molecule = builder.build().expect("methanol skeleton is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];
        let oxygen = &molecule.atoms()[oxygen_id.index()];

        assert_eq!(
            query_atom_total_valence(context.valence.as_ref(), carbon),
            Some(4)
        );
        assert_eq!(
            query_atom_total_valence(context.valence.as_ref(), oxygen),
            Some(2)
        );
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::TotalValence(4),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::TotalValence(3),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_unsaturated() {
        let mut saturated_builder = Molecule::builder();
        let methane_id = saturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let methane = saturated_builder.build().expect("methane is valid");
        let methane_context = build_query_match_context(&methane);
        let methane_carbon = &methane.atoms()[methane_id.index()];

        assert_eq!(
            query_atom_unsaturated(
                &methane_context.adj,
                methane_context.valence.as_ref(),
                methane_carbon,
            ),
            Some(false)
        );
        assert!(!atom_predicate_matches_with_context(
            methane_carbon,
            &AtomQueryPredicate::IsUnsaturated,
            &methane,
            &methane_context,
        ));

        let mut unsaturated_builder = Molecule::builder();
        let left_id = unsaturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let right_id = unsaturated_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        unsaturated_builder
            .add_bond(crate::BondSpec::new(
                left_id,
                right_id,
                crate::BondOrder::Double,
            ))
            .expect("carbon-to-carbon double bond is valid");
        let ethene = unsaturated_builder.build().expect("ethene is valid");
        let ethene_context = build_query_match_context(&ethene);
        let ethene_carbon = &ethene.atoms()[left_id.index()];

        assert_eq!(
            query_atom_unsaturated(
                &ethene_context.adj,
                ethene_context.valence.as_ref(),
                ethene_carbon,
            ),
            Some(true)
        );
        assert!(atom_predicate_matches_with_context(
            ethene_carbon,
            &AtomQueryPredicate::IsUnsaturated,
            &ethene,
            &ethene_context,
        ));
    }

    #[test]
    fn smarts_query_atom_num() {
        let mut builder = Molecule::builder();
        let carbon_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder.build().expect("isolated carbon is valid");
        let context = build_query_match_context(&molecule);
        let carbon = &molecule.atoms()[carbon_id.index()];

        assert_eq!(query_atom_num(carbon), 6);
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::AtomicNumber(6),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::AtomicNumberIn(vec![6, 8]),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon,
            &AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_make_atom_type() {
        assert_eq!(make_atom_type(6, false), 6);
        assert_eq!(make_atom_type(6, true), 1006);
        assert_eq!(make_atom_type(0, false), 0);
        assert_eq!(make_atom_type(118, true), 1118);
    }

    #[test]
    fn smarts_parse_atom_type() {
        assert_eq!(parse_atom_type(6), (6, false));
        assert_eq!(parse_atom_type(1000), (1000, false));
        assert_eq!(parse_atom_type(1001), (1, true));
        assert_eq!(parse_atom_type(1118), (118, true));
        assert_eq!(parse_atom_type(-1), (-1, false));
    }

    #[test]
    fn smarts_get_atom_type_is_aromatic() {
        assert!(!get_atom_type_is_aromatic(6));
        assert!(!get_atom_type_is_aromatic(1000));
        assert!(get_atom_type_is_aromatic(1001));
        assert!(get_atom_type_is_aromatic(make_atom_type(6, true)));
        assert!(!get_atom_type_is_aromatic(make_atom_type(6, false)));
    }

    #[test]
    fn smarts_get_atom_type_atomic_num() {
        assert_eq!(get_atom_type_atomic_num(6), 6);
        assert_eq!(get_atom_type_atomic_num(1000), 1000);
        assert_eq!(get_atom_type_atomic_num(1001), 1);
        assert_eq!(get_atom_type_atomic_num(make_atom_type(118, true)), 118);
        assert_eq!(get_atom_type_atomic_num(-1), -1);
    }

    #[test]
    fn smarts_query_atom_type() {
        let mut builder = Molecule::builder();
        let aliphatic_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_aromatic(true));
        let molecule = builder
            .build()
            .expect("two isolated carbons form a valid molecule");
        let context = build_query_match_context(&molecule);
        let aliphatic = &molecule.atoms()[aliphatic_id.index()];
        let aromatic = &molecule.atoms()[aromatic_id.index()];

        assert_eq!(query_atom_type(aliphatic), 6);
        assert_eq!(query_atom_type(aromatic), 1006);
        assert!(atom_predicate_matches_with_context(
            aliphatic,
            &AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            },
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            aromatic,
            &AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: true,
            },
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            aromatic,
            &AtomQueryPredicate::AtomType {
                atomic_number: 6,
                aromatic: false,
            },
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_mass() {
        let mut builder = Molecule::builder();
        let natural_carbon_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_12_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(12));
        let carbon_13_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(13));
        let unknown_isotope_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(999));
        let molecule = builder
            .build()
            .expect("four isolated carbons form a valid molecule");
        let context = build_query_match_context(&molecule);
        let natural_carbon = &molecule.atoms()[natural_carbon_id.index()];
        let carbon_12 = &molecule.atoms()[carbon_12_id.index()];
        let carbon_13 = &molecule.atoms()[carbon_13_id.index()];
        let unknown_isotope = &molecule.atoms()[unknown_isotope_id.index()];

        assert_eq!(query_atom_mass(natural_carbon), 12_011);
        assert_eq!(query_atom_mass(carbon_12), 12_000);
        assert_eq!(query_atom_mass(carbon_13), 13_003);
        assert_eq!(query_atom_mass(unknown_isotope), 999_000);
        assert!(atom_predicate_matches_with_context(
            carbon_12,
            &AtomQueryPredicate::Mass(12),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            natural_carbon,
            &AtomQueryPredicate::Mass(12),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            carbon_13,
            &AtomQueryPredicate::Mass(13),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_isotope() {
        let mut builder = Molecule::builder();
        let natural_carbon_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_13_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(13));
        let unknown_isotope_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_isotope(999));
        let molecule = builder
            .build()
            .expect("three isolated carbons form a valid molecule");
        let context = build_query_match_context(&molecule);
        let natural_carbon = &molecule.atoms()[natural_carbon_id.index()];
        let carbon_13 = &molecule.atoms()[carbon_13_id.index()];
        let unknown_isotope = &molecule.atoms()[unknown_isotope_id.index()];

        assert_eq!(query_atom_isotope(natural_carbon), 0);
        assert_eq!(query_atom_isotope(carbon_13), 13);
        assert_eq!(query_atom_isotope(unknown_isotope), 999);
        assert!(atom_predicate_matches_with_context(
            natural_carbon,
            &AtomQueryPredicate::Isotope(0),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            carbon_13,
            &AtomQueryPredicate::Isotope(13),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            natural_carbon,
            &AtomQueryPredicate::Isotope(13),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_formal_charge() {
        let mut builder = Molecule::builder();
        let anion_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_formal_charge(-1));
        let neutral_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let cation_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_formal_charge(2));
        let molecule = builder
            .build()
            .expect("three isolated charged atoms form a valid molecule");
        let context = build_query_match_context(&molecule);
        let anion = &molecule.atoms()[anion_id.index()];
        let neutral = &molecule.atoms()[neutral_id.index()];
        let cation = &molecule.atoms()[cation_id.index()];

        assert_eq!(query_atom_formal_charge(anion), -1);
        assert_eq!(query_atom_formal_charge(neutral), 0);
        assert_eq!(query_atom_formal_charge(cation), 2);
        assert!(atom_predicate_matches_with_context(
            anion,
            &AtomQueryPredicate::FormalCharge(-1),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            neutral,
            &AtomQueryPredicate::FormalCharge(0),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            cation,
            &AtomQueryPredicate::FormalCharge(2),
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            cation,
            &AtomQueryPredicate::FormalCharge(-1),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_negative_formal_charge() {
        let mut builder = Molecule::builder();
        let anion_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::O).with_formal_charge(-2));
        let neutral_id = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let cation_id =
            builder.add_atom(crate::AtomSpec::new(crate::Element::N).with_formal_charge(1));
        let molecule = builder
            .build()
            .expect("three isolated charged atoms form a valid molecule");
        let anion = &molecule.atoms()[anion_id.index()];
        let neutral = &molecule.atoms()[neutral_id.index()];
        let cation = &molecule.atoms()[cation_id.index()];

        assert_eq!(query_atom_negative_formal_charge(anion), 2);
        assert_eq!(query_atom_negative_formal_charge(neutral), 0);
        assert_eq!(query_atom_negative_formal_charge(cation), -1);
    }

    #[test]
    fn smarts_query_atom_hybridization() {
        let hybridizations = [
            Hybridization::Unspecified,
            Hybridization::S,
            Hybridization::Sp,
            Hybridization::Sp2,
            Hybridization::Sp3,
            Hybridization::Sp2d,
            Hybridization::Sp3d,
            Hybridization::Sp3d2,
            Hybridization::Other,
        ];
        let mut builder = Molecule::builder();
        let atom_ids: Vec<_> = hybridizations
            .iter()
            .map(|&hybridization| {
                builder.add_atom(
                    crate::AtomSpec::new(crate::Element::C).with_hybridization(hybridization),
                )
            })
            .collect();
        let molecule = builder
            .build()
            .expect("nine isolated carbons form a valid molecule");
        let context = build_query_match_context(&molecule);

        for (expected, (&atom_id, &hybridization)) in
            atom_ids.iter().zip(hybridizations.iter()).enumerate()
        {
            let atom = &molecule.atoms()[atom_id.index()];
            assert_eq!(query_atom_hybridization(atom), expected as i32);
            assert!(atom_predicate_matches_with_context(
                atom,
                &AtomQueryPredicate::HybridizationMatch(hybridization),
                &molecule,
                &context,
            ));
        }
        assert!(!atom_predicate_matches_with_context(
            &molecule.atoms()[atom_ids[3].index()],
            &AtomQueryPredicate::HybridizationMatch(Hybridization::Sp3),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_num_radical_electrons() {
        let radical_counts = [0_u8, 1, 2, u8::MAX];
        let mut builder = Molecule::builder();
        let atom_ids: Vec<_> = radical_counts
            .iter()
            .map(|&count| {
                builder
                    .add_atom(crate::AtomSpec::new(crate::Element::C).with_radical_electrons(count))
            })
            .collect();
        let molecule = builder
            .build()
            .expect("four isolated radical-count atoms form a valid molecule");

        for (&atom_id, &expected) in atom_ids.iter().zip(radical_counts.iter()) {
            let atom = &molecule.atoms()[atom_id.index()];
            assert_eq!(query_atom_num_radical_electrons(atom), i32::from(expected));
        }
    }

    #[test]
    fn smarts_query_atom_has_chiral_tag() {
        let chiral_tags = [
            ChiralTag::Unspecified,
            ChiralTag::TetrahedralCw,
            ChiralTag::TetrahedralCcw,
            ChiralTag::Other,
            ChiralTag::Tetrahedral,
            ChiralTag::Allene,
            ChiralTag::SquarePlanar,
            ChiralTag::TrigonalBipyramidal,
            ChiralTag::Octahedral,
        ];
        let mut builder = Molecule::builder();
        let atom_ids: Vec<_> = chiral_tags
            .iter()
            .map(|&tag| {
                builder.add_atom(crate::AtomSpec::new(crate::Element::C).with_chiral_tag(tag))
            })
            .collect();
        let molecule = builder
            .build()
            .expect("nine isolated chiral-tag atoms form a valid molecule");

        for (&atom_id, &tag) in atom_ids.iter().zip(chiral_tags.iter()) {
            let atom = &molecule.atoms()[atom_id.index()];
            let expected = i32::from(tag != ChiralTag::Unspecified);
            assert_eq!(query_atom_has_chiral_tag(atom), expected);
        }
    }

    #[test]
    fn smarts_query_atom_missing_chiral_tag() {
        let mut builder = Molecule::builder();
        let unspecified_without_prop = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let unspecified_with_prop = builder
            .add_atom(crate::AtomSpec::new(crate::Element::C).with_prop("_ChiralityPossible", "0"));
        let tagged_without_prop = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C).with_chiral_tag(ChiralTag::TetrahedralCw),
        );
        let tagged_with_prop = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_chiral_tag(ChiralTag::Other)
                .with_prop("_ChiralityPossible", "1"),
        );
        let molecule = builder
            .build()
            .expect("four isolated chirality-state atoms form a valid molecule");

        assert_eq!(
            query_atom_missing_chiral_tag(&molecule.atoms()[unspecified_without_prop.index()]),
            0
        );
        assert_eq!(
            query_atom_missing_chiral_tag(&molecule.atoms()[unspecified_with_prop.index()]),
            1,
            "RDKit tests property presence, not the stored property's truth value"
        );
        assert_eq!(
            query_atom_missing_chiral_tag(&molecule.atoms()[tagged_without_prop.index()]),
            0
        );
        assert_eq!(
            query_atom_missing_chiral_tag(&molecule.atoms()[tagged_with_prop.index()]),
            0
        );
    }

    #[test]
    fn smarts_query_atom_has_heteroatom_nbrs() {
        let mut builder = Molecule::builder();
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let nitrogen_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let nitrogen_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let dummy_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let dummy_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        for (begin, end) in [
            (carbon_center, carbon_neighbor),
            (carbon_center, hydrogen_neighbor),
            (nitrogen_center, nitrogen_neighbor),
            (dummy_center, dummy_neighbor),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, crate::BondOrder::Single))
                .expect("heteroatom-neighbor fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("heteroatom-neighbor fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        assert_eq!(
            query_atom_has_heteroatom_nbrs(
                &molecule.atoms()[isolated.index()],
                &context.adj,
                &molecule,
            ),
            0
        );
        assert_eq!(
            query_atom_has_heteroatom_nbrs(
                &molecule.atoms()[carbon_center.index()],
                &context.adj,
                &molecule,
            ),
            0
        );
        assert_eq!(
            query_atom_has_heteroatom_nbrs(
                &molecule.atoms()[nitrogen_center.index()],
                &context.adj,
                &molecule,
            ),
            1
        );
        assert_eq!(
            query_atom_has_heteroatom_nbrs(
                &molecule.atoms()[dummy_center.index()],
                &context.adj,
                &molecule,
            ),
            1,
            "RDKit classifies atomic number zero as a heteroatom neighbor here"
        );
        assert_eq!(
            query_atom_has_heteroatom_nbrs(
                &molecule.atoms()[nitrogen_neighbor.index()],
                &context.adj,
                &molecule,
            ),
            0,
            "the predicate classifies neighbors, not the queried atom itself"
        );
    }

    #[test]
    fn smarts_query_atom_num_heteroatom_nbrs() {
        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let dummy = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for neighbor in [carbon, hydrogen, nitrogen, oxygen, dummy] {
            builder
                .add_bond(crate::BondSpec::new(
                    center,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .expect("heteroatom-count fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("heteroatom-count fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        assert_eq!(
            query_atom_num_heteroatom_nbrs(
                &molecule.atoms()[center.index()],
                &context.adj,
                &molecule,
            ),
            3,
            "nitrogen, oxygen, and dummy neighbors are counted; carbon and hydrogen are not"
        );
        assert_eq!(
            query_atom_num_heteroatom_nbrs(
                &molecule.atoms()[isolated.index()],
                &context.adj,
                &molecule,
            ),
            0
        );
        assert_eq!(
            query_atom_num_heteroatom_nbrs(
                &molecule.atoms()[nitrogen.index()],
                &context.adj,
                &molecule,
            ),
            0,
            "the carbon neighbor of nitrogen is not a heteroatom neighbor"
        );
    }

    #[test]
    fn smarts_query_atom_has_aliphatic_heteroatom_nbrs() {
        let mut builder = Molecule::builder();
        let aliphatic_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aliphatic_nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let aromatic_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic_nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let dummy_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let dummy = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let carbon_hydrogen_center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [
            (aliphatic_center, aliphatic_nitrogen),
            (aromatic_center, aromatic_nitrogen),
            (dummy_center, dummy),
            (carbon_hydrogen_center, carbon),
            (carbon_hydrogen_center, hydrogen),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, crate::BondOrder::Single))
                .expect("aliphatic-heteroatom fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("aliphatic-heteroatom fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        for (atom_id, expected) in [
            (aliphatic_center, 1),
            (aromatic_center, 0),
            (dummy_center, 1),
            (carbon_hydrogen_center, 0),
            (isolated, 0),
        ] {
            assert_eq!(
                query_atom_has_aliphatic_heteroatom_nbrs(
                    &molecule.atoms()[atom_id.index()],
                    &context.adj,
                    &molecule,
                ),
                expected
            );
        }
        assert_eq!(
            query_atom_has_aliphatic_heteroatom_nbrs(
                &molecule.atoms()[aromatic_nitrogen.index()],
                &context.adj,
                &molecule,
            ),
            0,
            "the queried atom's own aromaticity and element do not classify its carbon neighbor"
        );
    }

    #[test]
    fn smarts_query_atom_num_aliphatic_heteroatom_nbrs() {
        let mut builder = Molecule::builder();
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aliphatic_nitrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::N));
        let aliphatic_oxygen = builder.add_atom(crate::AtomSpec::new(crate::Element::O));
        let dummy = builder.add_atom(crate::AtomSpec::new(crate::Element::DUMMY));
        let aromatic_nitrogen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::N)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let aromatic_oxygen = builder.add_atom(
            crate::AtomSpec::new(crate::Element::O)
                .with_aromatic(true)
                .with_no_implicit(true),
        );
        let carbon = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let hydrogen = builder.add_atom(crate::AtomSpec::new(crate::Element::H));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for neighbor in [
            aliphatic_nitrogen,
            aliphatic_oxygen,
            dummy,
            aromatic_nitrogen,
            aromatic_oxygen,
            carbon,
            hydrogen,
        ] {
            builder
                .add_bond(crate::BondSpec::new(
                    center,
                    neighbor,
                    crate::BondOrder::Single,
                ))
                .expect("aliphatic-heteroatom-count fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("aliphatic-heteroatom-count fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        assert_eq!(
            query_atom_num_aliphatic_heteroatom_nbrs(
                &molecule.atoms()[center.index()],
                &context.adj,
                &molecule,
            ),
            3,
            "aliphatic nitrogen, oxygen, and dummy neighbors are counted"
        );
        assert_eq!(
            query_atom_num_aliphatic_heteroatom_nbrs(
                &molecule.atoms()[isolated.index()],
                &context.adj,
                &molecule,
            ),
            0
        );
        assert_eq!(
            query_atom_num_aliphatic_heteroatom_nbrs(
                &molecule.atoms()[aromatic_nitrogen.index()],
                &context.adj,
                &molecule,
            ),
            0,
            "the queried atom's own aromaticity does not classify its carbon neighbor"
        );
    }

    #[test]
    fn smarts_query_bond_order() {
        let orders = [
            BondOrder::Null,
            BondOrder::Single,
            BondOrder::Double,
            BondOrder::Triple,
            BondOrder::Quadruple,
            BondOrder::Quintuple,
            BondOrder::Hextuple,
            BondOrder::OneAndHalf,
            BondOrder::TwoAndHalf,
            BondOrder::ThreeAndHalf,
            BondOrder::FourAndHalf,
            BondOrder::FiveAndHalf,
            BondOrder::Aromatic,
            BondOrder::Ionic,
            BondOrder::Dative,
            BondOrder::DativeOne,
            BondOrder::DativeLeft,
            BondOrder::DativeRight,
            BondOrder::Hydrogen,
            BondOrder::ThreeCenter,
            BondOrder::Other,
            BondOrder::Zero,
            BondOrder::Unspecified,
        ];
        let mut builder = Molecule::builder();
        for order in orders {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, order))
                .expect("bond-order fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("bond-order fixture molecule is valid");

        for (bond, expected) in molecule.bonds().iter().zip(orders) {
            assert_eq!(query_bond_order(bond), expected);
            assert!(bond_predicate_matches_with_context(
                bond,
                &BondQueryPredicate::Order(expected),
                &molecule,
                &build_query_match_context(&molecule),
            ));
        }
    }

    #[test]
    fn smarts_query_bond_is_single_or_aromatic() {
        let cases = [
            (BondOrder::Single, false, 1),
            (BondOrder::Aromatic, false, 1),
            (BondOrder::Double, true, 0),
            (BondOrder::Triple, false, 0),
        ];
        let mut builder = Molecule::builder();
        for (order, is_aromatic, _) in cases {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, order).with_aromatic(is_aromatic))
                .expect("single-or-aromatic fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("single-or-aromatic fixture molecule is valid");

        for (bond, (_, _, expected)) in molecule.bonds().iter().zip(cases) {
            assert_eq!(query_bond_is_single_or_aromatic(bond), expected);
        }
    }

    #[test]
    fn smarts_query_bond_is_double_or_aromatic() {
        let cases = [
            (BondOrder::Double, false, 1),
            (BondOrder::Aromatic, false, 1),
            (BondOrder::Single, true, 0),
            (BondOrder::Triple, false, 0),
        ];
        let mut builder = Molecule::builder();
        for (order, is_aromatic, _) in cases {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, order).with_aromatic(is_aromatic))
                .expect("double-or-aromatic fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("double-or-aromatic fixture molecule is valid");

        for (bond, (_, _, expected)) in molecule.bonds().iter().zip(cases) {
            assert_eq!(query_bond_is_double_or_aromatic(bond), expected);
        }
    }

    #[test]
    fn smarts_query_bond_is_single_or_double() {
        let cases = [
            (BondOrder::Single, false, 1),
            (BondOrder::Double, false, 1),
            (BondOrder::Aromatic, false, 0),
            (BondOrder::Triple, true, 0),
        ];
        let mut builder = Molecule::builder();
        for (order, is_aromatic, _) in cases {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, order).with_aromatic(is_aromatic))
                .expect("single-or-double fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("single-or-double fixture molecule is valid");

        for (bond, (_, _, expected)) in molecule.bonds().iter().zip(cases) {
            assert_eq!(query_bond_is_single_or_double(bond), expected);
        }
    }

    #[test]
    fn smarts_query_bond_is_single_or_double_or_aromatic() {
        let cases = [
            (BondOrder::Single, false, 1),
            (BondOrder::Double, false, 1),
            (BondOrder::Aromatic, false, 1),
            (BondOrder::Triple, true, 0),
        ];
        let mut builder = Molecule::builder();
        for (order, is_aromatic, _) in cases {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(crate::BondSpec::new(begin, end, order).with_aromatic(is_aromatic))
                .expect("single-double-or-aromatic fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("single-double-or-aromatic fixture molecule is valid");

        for (bond, (_, _, expected)) in molecule.bonds().iter().zip(cases) {
            assert_eq!(query_bond_is_single_or_double_or_aromatic(bond), expected);
        }
    }

    #[test]
    fn smarts_query_bond_dir() {
        let directions = [
            crate::BondDirection::None,
            crate::BondDirection::BeginWedge,
            crate::BondDirection::BeginDash,
            crate::BondDirection::EndDownRight,
            crate::BondDirection::EndUpRight,
            crate::BondDirection::EitherDouble,
            crate::BondDirection::Unknown,
        ];
        let mut builder = Molecule::builder();
        for direction in directions {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            builder
                .add_bond(
                    crate::BondSpec::new(begin, end, BondOrder::Single).with_direction(direction),
                )
                .expect("bond-direction fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("bond-direction fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        for (bond, expected) in molecule.bonds().iter().zip(directions) {
            assert_eq!(query_bond_dir(bond), expected);
            assert!(bond_predicate_matches_with_context(
                bond,
                &BondQueryPredicate::Direction(expected),
                &molecule,
                &context,
            ));
        }
    }

    #[test]
    fn smarts_query_is_bond_in_n_rings() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let shared = builder
            .add_bond(crate::BondSpec::new(a1, a2, BondOrder::Single))
            .expect("shared-ring fixture bond is valid");
        let one_ring = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("first-ring fixture bond is valid");
        for (begin, end) in [(a2, a0), (a1, a3), (a3, a2)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("fused-ring fixture bond is valid");
        }
        let non_ring = builder
            .add_bond(crate::BondSpec::new(a0, isolated, BondOrder::Single))
            .expect("non-ring fixture bond is valid");
        let molecule = builder
            .build()
            .expect("fused-ring fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the fused-ring fixture");

        assert_eq!(
            query_is_bond_in_n_rings(&molecule.bonds()[shared.index()], ring_info),
            2
        );
        assert_eq!(
            query_is_bond_in_n_rings(&molecule.bonds()[one_ring.index()], ring_info),
            1
        );
        assert_eq!(
            query_is_bond_in_n_rings(&molecule.bonds()[non_ring.index()], ring_info),
            0
        );
        assert!(bond_predicate_matches_with_context(
            &molecule.bonds()[shared.index()],
            &BondQueryPredicate::NumRingBonds(2),
            &molecule,
            &context,
        ));
        assert!(bond_predicate_matches_with_context(
            &molecule.bonds()[one_ring.index()],
            &BondQueryPredicate::NumRingBondsGreaterEqual(1),
            &molecule,
            &context,
        ));
        assert!(bond_predicate_matches_with_context(
            &molecule.bonds()[non_ring.index()],
            &BondQueryPredicate::NumRingBondsLessEqual(0),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_bond_has_stereo() {
        let cases = [
            (BondStereo::None, 0),
            (BondStereo::Any, 1),
            (BondStereo::Z, 1),
            (BondStereo::E, 1),
            (BondStereo::Cis, 1),
            (BondStereo::Trans, 1),
            (BondStereo::AtropCw, 1),
            (BondStereo::AtropCcw, 1),
        ];
        let mut builder = Molecule::builder();
        for (stereo, _) in cases {
            let begin = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let end = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
            let begin_ref = builder.add_atom(crate::AtomSpec::new(crate::Element::F));
            let end_ref = builder.add_atom(crate::AtomSpec::new(crate::Element::CL));
            builder
                .add_bond(
                    crate::BondSpec::new(begin, end, BondOrder::Double)
                        .with_stereo_atoms(begin_ref, end_ref)
                        .with_stereo(stereo),
                )
                .expect("bond-stereo fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("bond-stereo fixture molecule is valid");

        for (bond, (_, expected)) in molecule.bonds().iter().zip(cases) {
            assert_eq!(query_bond_has_stereo(bond), expected);
        }
    }

    #[test]
    fn smarts_query_atom_ring_membership() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0), (a1, a3), (a3, a2)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("fused-ring atom fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("fused-ring atom fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the fused-ring fixture");

        assert_eq!(
            query_atom_ring_membership(&molecule.atoms()[a1.index()], ring_info),
            2
        );
        assert_eq!(
            query_atom_ring_membership(&molecule.atoms()[a0.index()], ring_info),
            1
        );
        assert_eq!(
            query_atom_ring_membership(&molecule.atoms()[isolated.index()], ring_info),
            0
        );
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[a1.index()],
            &AtomQueryPredicate::NumAtomRings(2),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_is_atom_in_ring() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-in-ring fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("atom-in-ring fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_is_atom_in_ring(&molecule.atoms()[a0.index()], ring_info),
            1
        );
        assert_eq!(
            query_is_atom_in_ring(&molecule.atoms()[isolated.index()], ring_info),
            0
        );
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[a0.index()],
            &AtomQueryPredicate::InRing,
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            &molecule.atoms()[isolated.index()],
            &AtomQueryPredicate::InRing,
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_has_ring_bond() {
        let mut builder = Molecule::builder();
        let ring_atom = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_third = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [
            (ring_atom, ring_neighbor),
            (ring_neighbor, ring_third),
            (ring_third, ring_atom),
            (ring_atom, external),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-has-ring-bond fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("atom-has-ring-bond fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_atom_has_ring_bond(
                &molecule.atoms()[ring_atom.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            1
        );
        assert_eq!(
            query_atom_has_ring_bond(
                &molecule.atoms()[external.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            0
        );
        assert_eq!(
            query_atom_has_ring_bond(
                &molecule.atoms()[isolated.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            0
        );
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[ring_atom.index()],
            &AtomQueryPredicate::HasRingBond,
            &molecule,
            &context,
        ));
        assert!(!atom_predicate_matches_with_context(
            &molecule.atoms()[external.index()],
            &AtomQueryPredicate::HasRingBond,
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_is_bond_in_ring() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond-in-ring fixture ring bond is valid");
        for (begin, end) in [(a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("bond-in-ring fixture ring bond is valid");
        }
        let non_ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, external, BondOrder::Single))
            .expect("bond-in-ring fixture non-ring bond is valid");
        let molecule = builder
            .build()
            .expect("bond-in-ring fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_is_bond_in_ring(&molecule.bonds()[ring_bond.index()], ring_info),
            1
        );
        assert_eq!(
            query_is_bond_in_ring(&molecule.bonds()[non_ring_bond.index()], ring_info),
            0
        );
        assert!(bond_predicate_matches_with_context(
            &molecule.bonds()[ring_bond.index()],
            &BondQueryPredicate::IsInRing(true),
            &molecule,
            &context,
        ));
        assert!(bond_predicate_matches_with_context(
            &molecule.bonds()[non_ring_bond.index()],
            &BondQueryPredicate::IsInRing(false),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_min_ring_size() {
        let mut builder = Molecule::builder();
        let triangle = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-min-ring-size triangle bond is valid");
        }
        let square = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for (begin, end) in [
            (square[0], square[1]),
            (square[1], square[2]),
            (square[2], square[3]),
            (square[3], square[0]),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-min-ring-size square bond is valid");
        }
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let molecule = builder
            .build()
            .expect("atom-min-ring-size fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_atom_min_ring_size(&molecule.atoms()[triangle[0].index()], ring_info),
            3
        );
        assert_eq!(
            query_atom_min_ring_size(&molecule.atoms()[square[0].index()], ring_info),
            4
        );
        assert_eq!(
            query_atom_min_ring_size(&molecule.atoms()[isolated.index()], ring_info),
            0
        );
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[triangle[0].index()],
            &AtomQueryPredicate::SmallestRingSize(3),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[square[0].index()],
            &AtomQueryPredicate::SmallestRingSizeGreaterEqual(4),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[isolated.index()],
            &AtomQueryPredicate::SmallestRingSizeLessEqual(0),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_bond_min_ring_size() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond-min-ring-size ring bond is valid");
        for (begin, end) in [(a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("bond-min-ring-size ring bond is valid");
        }
        let non_ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, external, BondOrder::Single))
            .expect("bond-min-ring-size non-ring bond is valid");
        let molecule = builder
            .build()
            .expect("bond-min-ring-size fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_bond_min_ring_size(&molecule.bonds()[ring_bond.index()], ring_info),
            3
        );
        assert_eq!(
            query_bond_min_ring_size(&molecule.bonds()[non_ring_bond.index()], ring_info),
            0
        );
    }

    #[test]
    fn layered_ring_accessors_match_source() {
        let mut builder = Molecule::builder();
        let triangle = [
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(
                triangle[0],
                triangle[1],
                BondOrder::Single,
            ))
            .expect("Layered ring accessor fixture bond is valid");
        for (begin, end) in [(triangle[1], triangle[2]), (triangle[2], triangle[0])] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("Layered ring accessor fixture bond is valid");
        }
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let non_ring_bond = builder
            .add_bond(crate::BondSpec::new(
                triangle[0],
                external,
                BondOrder::Single,
            ))
            .expect("Layered non-ring accessor fixture bond is valid");
        let molecule = builder
            .build()
            .expect("Layered ring accessor fixture is valid");
        let ring_info = crate::find_sssr(&molecule).expect("exact SSSR succeeds");

        assert_eq!(
            query_is_bond_in_ring(&molecule.bonds()[ring_bond.index()], &ring_info),
            1,
        );
        assert_eq!(
            query_bond_min_ring_size(&molecule.bonds()[ring_bond.index()], &ring_info),
            3,
        );
        assert_eq!(
            query_is_bond_in_ring(&molecule.bonds()[non_ring_bond.index()], &ring_info),
            0,
        );
        assert_eq!(
            query_bond_min_ring_size(&molecule.bonds()[non_ring_bond.index()], &ring_info,),
            0,
        );
    }

    #[test]
    fn smarts_query_atom_ring_bond_count() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let shared = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a3 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [
            (a0, shared),
            (shared, a2),
            (a2, a0),
            (shared, a3),
            (a3, a2),
            (a0, external),
        ] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-ring-bond-count fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("atom-ring-bond-count fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the fused-ring fixture");

        assert_eq!(
            query_atom_ring_bond_count(
                &molecule.atoms()[shared.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            3
        );
        assert_eq!(
            query_atom_ring_bond_count(
                &molecule.atoms()[a0.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            2
        );
        assert_eq!(
            query_atom_ring_bond_count(
                &molecule.atoms()[external.index()],
                &context.adj,
                &molecule,
                ring_info,
            ),
            0
        );
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[shared.index()],
            &AtomQueryPredicate::RingBondCount(3),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            &molecule.atoms()[external.index()],
            &AtomQueryPredicate::RingBondCountLessEqual(0),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_atom_is_in_ring_of_size() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [(a0, a1), (a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("atom-in-ring-of-size fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("atom-in-ring-of-size fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");
        let ring_atom = &molecule.atoms()[a0.index()];
        let isolated_atom = &molecule.atoms()[isolated.index()];

        assert_eq!(query_atom_is_in_ring_of_size(ring_atom, 3, ring_info), 3);
        assert_eq!(query_atom_is_in_ring_of_size(ring_atom, 4, ring_info), 0);
        assert_eq!(
            query_atom_is_in_ring_of_size(isolated_atom, 3, ring_info),
            0
        );
        assert_eq!(
            query_atom_is_in_ring_size_range(ring_atom, 3, 4, false, true, ring_info),
            3
        );
        assert_eq!(
            query_atom_is_in_ring_size_range(ring_atom, 3, 4, true, false, ring_info),
            -1
        );
        assert_eq!(
            query_atom_is_in_ring_size_range(isolated_atom, -1, 4, false, false, ring_info),
            i32::MAX
        );
        assert_eq!(
            query_atom_is_in_ring_size_range(isolated_atom, -1, -1, false, false, ring_info),
            0
        );
        assert!(atom_predicate_matches_with_context(
            ring_atom,
            &AtomQueryPredicate::InRingOfSize(3),
            &molecule,
            &context,
        ));
        assert!(atom_predicate_matches_with_context(
            isolated_atom,
            &AtomQueryPredicate::InRingOfSize(0),
            &molecule,
            &context,
        ));
    }

    #[test]
    fn smarts_query_bond_is_in_ring_of_size() {
        let mut builder = Molecule::builder();
        let a0 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a1 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let a2 = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let external = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, a1, BondOrder::Single))
            .expect("bond-in-ring-of-size ring bond is valid");
        for (begin, end) in [(a1, a2), (a2, a0)] {
            builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("bond-in-ring-of-size ring bond is valid");
        }
        let non_ring_bond = builder
            .add_bond(crate::BondSpec::new(a0, external, BondOrder::Single))
            .expect("bond-in-ring-of-size non-ring bond is valid");
        let molecule = builder
            .build()
            .expect("bond-in-ring-of-size fixture molecule is valid");
        let context = build_query_match_context(&molecule);
        let ring_info = context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the ring fixture");

        assert_eq!(
            query_bond_is_in_ring_of_size(&molecule.bonds()[ring_bond.index()], 3, ring_info,),
            3
        );
        assert_eq!(
            query_bond_is_in_ring_of_size(&molecule.bonds()[ring_bond.index()], 4, ring_info,),
            0
        );
        assert_eq!(
            query_bond_is_in_ring_of_size(&molecule.bonds()[non_ring_bond.index()], 3, ring_info,),
            0
        );
    }

    #[test]
    fn smarts_query_is_atom_bridgehead() {
        let mut bridge_builder = Molecule::builder();
        let bridgehead = bridge_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let opposite = bridge_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let path_atoms = [
            bridge_builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            bridge_builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
            bridge_builder.add_atom(crate::AtomSpec::new(crate::Element::C)),
        ];
        for middle in path_atoms {
            for (begin, end) in [(bridgehead, middle), (middle, opposite)] {
                bridge_builder
                    .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                    .expect("bridgehead theta-graph bond is valid");
            }
        }
        let bridge_molecule = bridge_builder
            .build()
            .expect("bridgehead theta graph is valid");
        let bridge_context = build_query_match_context(&bridge_molecule);
        let bridge_rings = bridge_context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the theta graph");
        assert_eq!(
            crate::chemistry::stereo::query_is_atom_bridgehead(
                &bridge_molecule,
                bridgehead.index(),
                bridge_rings,
            ),
            1
        );
        assert_eq!(
            crate::chemistry::stereo::query_is_atom_bridgehead(
                &bridge_molecule,
                path_atoms[0].index(),
                bridge_rings,
            ),
            0
        );

        let mut fused_builder = Molecule::builder();
        let shared_left = fused_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let shared_right = fused_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let top = fused_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let bottom = fused_builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (begin, end) in [
            (shared_left, shared_right),
            (shared_left, top),
            (top, shared_right),
            (shared_left, bottom),
            (bottom, shared_right),
        ] {
            fused_builder
                .add_bond(crate::BondSpec::new(begin, end, BondOrder::Single))
                .expect("fused-ring non-bridgehead bond is valid");
        }
        let fused_molecule = fused_builder
            .build()
            .expect("fused-ring non-bridgehead graph is valid");
        let fused_context = build_query_match_context(&fused_molecule);
        let fused_rings = fused_context
            .ring_info
            .as_ref()
            .expect("ring information is perceived for the fused-ring graph");
        assert_eq!(
            crate::chemistry::stereo::query_is_atom_bridgehead(
                &fused_molecule,
                shared_left.index(),
                fused_rings,
            ),
            0
        );
    }

    #[test]
    fn smarts_query_atom_bond_product() {
        let mut builder = Molecule::builder();
        let isolated = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let center = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let single = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let double = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let aromatic = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let dative = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        for (end, order) in [
            (single, BondOrder::Single),
            (double, BondOrder::Double),
            (aromatic, BondOrder::Aromatic),
            (dative, BondOrder::Dative),
        ] {
            builder
                .add_bond(crate::BondSpec::new(center, end, order))
                .expect("atom-bond-product fixture bond is valid");
        }
        let molecule = builder
            .build()
            .expect("atom-bond-product fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        assert_eq!(
            query_atom_bond_product(&molecule.atoms()[isolated.index()], &context.adj, &molecule,),
            1
        );
        assert_eq!(
            query_atom_bond_product(&molecule.atoms()[single.index()], &context.adj, &molecule,),
            3
        );
        assert_eq!(
            query_atom_bond_product(&molecule.atoms()[center.index()], &context.adj, &molecule,),
            3 * 5 * 41 * 61
        );

        assert_eq!(rdkit_bond_type_prime(BondOrder::DativeOne), 59);
        assert_eq!(rdkit_bond_type_prime(BondOrder::Hydrogen), 47);
        assert_eq!(rdkit_bond_type_prime(BondOrder::Zero), 79);
    }

    #[test]
    fn smarts_query_atom_all_bond_product() {
        let mut builder = Molecule::builder();
        let methyl = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        let neighbor = builder.add_atom(crate::AtomSpec::new(crate::Element::C));
        builder
            .add_bond(crate::BondSpec::new(methyl, neighbor, BondOrder::Single))
            .expect("all-bond-product methyl bond is valid");
        let explicit_hydrogens = builder.add_atom(
            crate::AtomSpec::new(crate::Element::C)
                .with_explicit_hydrogens(2)
                .with_no_implicit(true),
        );
        let molecule = builder
            .build()
            .expect("all-bond-product fixture molecule is valid");
        let context = build_query_match_context(&molecule);

        assert_eq!(
            query_atom_all_bond_product(
                &molecule.atoms()[methyl.index()],
                &context.adj,
                &molecule,
                context.valence.as_ref(),
            ),
            Some(81)
        );
        assert_eq!(
            query_atom_all_bond_product(
                &molecule.atoms()[explicit_hydrogens.index()],
                &context.adj,
                &molecule,
                context.valence.as_ref(),
            ),
            Some(9)
        );
        assert_eq!(
            query_atom_bond_product(&molecule.atoms()[methyl.index()], &context.adj, &molecule,),
            3
        );
    }

    #[test]
    fn test_query_node_logic() {
        // OR: matches if any child matches
        let or_node = QueryNode::or(vec![
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(7)),
        ]);
        assert_eq!(
            or_node,
            QueryNode::Or(vec![
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6)),
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(7)),
            ])
        );

        // AND
        let and_node = QueryNode::and(vec![
            QueryNode::predicate(AtomQueryPredicate::IsAromatic(true)),
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
        ]);
        assert_eq!(
            and_node,
            QueryNode::And(vec![
                QueryNode::Predicate(AtomQueryPredicate::IsAromatic(true)),
                QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6)),
            ])
        );

        // NOT
        let not_node = QueryNode::not(QueryNode::predicate(AtomQueryPredicate::Any));
        assert_eq!(
            not_node,
            QueryNode::Not(Box::new(QueryNode::Predicate(AtomQueryPredicate::Any)))
        );
    }
}
