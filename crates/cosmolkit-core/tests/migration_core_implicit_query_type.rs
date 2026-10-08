use cosmolkit_core::calculate_implicit_valence_for_topology;
use cosmolkit_model::{
    Atom, AtomId, AtomSpec, Bond, BondId, BondQueryPredicate as P, BondSpec, QueryNode as Q,
    TopologyBlock,
};
use cosmolkit_types::{BondOrder, Element};

fn order(value: BondOrder) -> Q<P> {
    Q::predicate(P::Order(value))
}
fn ring() -> Q<P> {
    Q::predicate(P::IsInRing(true))
}
fn topology(query: Option<Q<P>>) -> TopologyBlock {
    let mut spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single);
    if let Some(query) = query {
        spec = spec.with_query(query);
    }
    TopologyBlock::try_from_parts(
        vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::F)),
        ],
        vec![Bond::from_spec(BondId::new(0), spec)],
        vec![],
        vec![],
    )
    .unwrap()
}

#[test]
fn implicit_query_guard_matches_source_root_negation_and_ordered_seen_state() {
    // Exact native QueryOps.cpp::hasComplexBondTypeQueryHelper: root plain
    // BondOrder alone is false; legacy multiorder and negated order are true.
    // Only direct BondOrder child descriptions update seen for later siblings.
    let single = || order(BondOrder::Single);
    let double = || order(BondOrder::Double);
    let cases = vec![
        ("ordinary", None, 3),
        ("singleOrder", Some(single()), 3),
        (
            "singleOrAromaticFactory",
            Some(Q::predicate(P::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Aromatic,
            ]))),
            0,
        ),
        (
            "doubleOrAromaticFactory",
            Some(Q::predicate(P::OrderIn(vec![
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))),
            0,
        ),
        (
            "singleOrDoubleFactory",
            Some(Q::predicate(P::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
            ]))),
            0,
        ),
        (
            "singleDoubleAromaticFactory",
            Some(Q::predicate(P::OrderIn(vec![
                BondOrder::Single,
                BondOrder::Double,
                BondOrder::Aromatic,
            ]))),
            0,
        ),
        ("negatedOrder", Some(Q::not(single())), 0),
        ("andTwoOrders", Some(Q::and(vec![single(), double()])), 0),
        ("orTwoOrders", Some(Q::or(vec![single(), double()])), 0),
        ("xorTwoOrders", Some(Q::xor(vec![single(), double()])), 0),
        ("orderAndRing", Some(Q::and(vec![single(), ring()])), 3),
        ("orderOrRing", Some(Q::or(vec![single(), ring()])), 3),
        (
            "negatedCompositeOrderRing",
            Some(Q::not(Q::and(vec![single(), ring()]))),
            3,
        ),
        (
            "nestedOrderBeforeDirect",
            Some(Q::and(vec![Q::and(vec![single()]), double()])),
            3,
        ),
        (
            "directOrderBeforeNested",
            Some(Q::and(vec![single(), Q::and(vec![double()])])),
            0,
        ),
        ("nonOrderNegation", Some(Q::not(ring())), 3),
    ];
    for (label, query, expected) in cases {
        let graph = topology(query);
        let before = graph.clone();
        for (strict, check_it) in [(false, false), (false, true), (true, false), (true, true)] {
            assert_eq!(
                calculate_implicit_valence_for_topology(
                    &graph,
                    AtomId::new(0),
                    1,
                    strict,
                    check_it
                ),
                Ok(expected),
                "{label},strict={strict},check={check_it}"
            );
        }
        assert_eq!(graph, before);
    }
}
