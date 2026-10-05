//! Source root dispatch, including ignored second-child negation and child order.
use cosmolkit_model::{AtomId, QueryAtom, QueryAtomIdentity, QueryGraph};
use cosmolkit_search::{AtomQueryPredicate as A, QueryNode as Q, is_query_atom_aromatic};

fn aromaticity(query: Q<A>) -> bool {
    let atom =
        QueryAtom::from_identity_parts(AtomId::new(0), QueryAtomIdentity::AtomicNumber(0), query);
    let graph = QueryGraph::from_parts(
        vec![atom],
        vec![],
        Default::default(),
        vec![],
        vec![],
        vec![],
    )
    .unwrap();
    assert_eq!(graph.atoms()[0].atomic_number(), 0);
    is_query_atom_aromatic(&graph.atoms()[0], &graph)
}
#[test]
fn negated_aromatic_and_aliphatic_roots_follow_source() {
    assert!(aromaticity(Q::predicate(A::IsAromatic(true))));
    assert!(!aromaticity(Q::not(Q::predicate(A::IsAromatic(true)))));
    assert!(!aromaticity(Q::predicate(A::IsAromatic(false))));
    assert!(aromaticity(Q::not(Q::predicate(A::IsAromatic(false)))));
    assert!(aromaticity(Q::predicate(A::AtomType {
        atomic_number: 6,
        aromatic: true
    })));
    assert!(!aromaticity(Q::not(Q::predicate(A::AtomType {
        atomic_number: 6,
        aromatic: true
    }))));
}
#[test]
fn and_checks_source_order_and_ignores_second_child_negation() {
    let number = || Q::predicate(A::AtomicNumber(6));
    let aromatic = || Q::predicate(A::IsAromatic(true));
    assert!(aromaticity(Q::and(vec![number(), aromatic()])));
    assert!(aromaticity(Q::and(vec![number(), Q::not(aromatic())])));
    assert!(!aromaticity(Q::and(vec![aromatic(), number()])));
    assert!(!aromaticity(Q::not(Q::and(vec![number(), aromatic()]))));
    assert!(!aromaticity(Q::and(vec![
        number(),
        Q::predicate(A::FormalCharge(0)),
        aromatic()
    ])));
}
#[test]
fn unrelated_and_boolean_roots_remain_nonaromatic() {
    for query in [
        Q::predicate(A::Any),
        Q::predicate(A::FormalCharge(0)),
        Q::or(vec![Q::predicate(A::IsAromatic(true))]),
        Q::xor(vec![Q::predicate(A::IsAromatic(true))]),
        Q::not(Q::not(Q::predicate(A::IsAromatic(true)))),
    ] {
        assert!(!aromaticity(query));
    }
}
