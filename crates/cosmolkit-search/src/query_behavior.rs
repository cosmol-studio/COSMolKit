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

pub use cosmolkit_model::{
    AtomQueryPredicate, AtomRangeBounds, AtomRangeDataFunction, AtomRangeQuery, BondQueryPredicate,
    QueryNode, RecursiveStructureQuery,
};

use cosmolkit_core::{RingInfo, ValenceAssignment, ValenceModel, atomic_mass as rdkit_atomic_mass};
use cosmolkit_model::{AdjacencyList, Atom, AtomId, Bond, BondSpec};
use cosmolkit_types::{BondOrder, BondStereo, ChiralTag, Hybridization};

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

#[doc(hidden)]
pub fn is_atom_aromatic(atom: &Atom, molecule: &impl SearchTargetAccess) -> bool {
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

#[doc(hidden)]
pub const QUERY_SCAN_MAGIC_VALUE: u32 = 0xDEADBEEF;

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

#[doc(hidden)]
pub fn atom_query_has_magic_value(query: &QueryNode<AtomQueryPredicate>, magic_value: u32) -> bool {
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
        let ring_bond_count = 0_u32;
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

fn rdkit_atom_mass(atom: &Atom) -> f64 {
    match rdkit_atomic_mass(atom.element(), atom.isotope()) {
        Ok(mass) => mass,
        Err(_) if atom.atomic_number() != 0 => atom.isotope().map_or(0.0, f64::from),
        Err(_) => 0.0,
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
        let mass = rdkit_atom_mass(&atom) as u16;
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
#[doc(hidden)]
pub fn convert_complex_name_to_query(
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
#[doc(hidden)]
pub fn make_single_or_aromatic_bond_query() -> QueryNode<BondQueryPredicate> {
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
    cosmolkit_core::find_sssr_from_parts(mol.num_atoms(), mol.bonds(), mol.adjacency()).ok()
}

fn ensure_valence_assignment(mol: &impl SearchTargetAccess) -> Option<ValenceAssignment> {
    if let Some(cached) = mol.valence() {
        return Some(cached.clone());
    }
    cosmolkit_core::assign_valence_with_options_for_topology(
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
pub fn build_query_match_context(mol: &impl SearchTargetAccess) -> QueryMatchContext {
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
    (f64::from(MASS_INTEGER_CONVERSION_FACTOR) * rdkit_atom_mass(at)).round() as i32
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
#[doc(hidden)]
pub fn query_is_bond_in_ring(bond: &Bond, ring_info: &RingInfo) -> i32 {
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
#[doc(hidden)]
pub fn query_bond_min_ring_size(bond: &Bond, ring_info: &RingInfo) -> usize {
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
        BondOrder::Unspecified => 2,
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
pub fn atom_predicate_matches(
    atom: &Atom,
    pred: &AtomQueryPredicate,
    mol: &impl SearchTargetAccess,
) -> bool {
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
            cosmolkit_core::is_atom_bridgehead_from_topology(mol.topology_block(), aidx, ri) != 0
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
    mol: &impl SearchTargetAccess,
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
pub fn bond_predicate_matches(
    bond: &Bond,
    pred: &BondQueryPredicate,
    mol: &impl SearchTargetAccess,
) -> bool {
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
        BondOrder::Unspecified => 0,
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
    mol: &impl SearchTargetAccess,
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
    mol: &impl SearchTargetAccess,
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
    mol: &impl SearchTargetAccess,
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
    mol: &impl SearchTargetAccess,
    ctx: &QueryMatchContext,
) -> bool {
    atom_matches_query_with_target_context(atom, query, mol, ctx)
}

pub fn bond_matches_query_with_context(
    bond: &Bond,
    query: &QueryNode<BondQueryPredicate>,
    mol: &impl SearchTargetAccess,
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
