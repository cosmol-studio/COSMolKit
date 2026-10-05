//! Query complexity branches from pinned RDKit QueryOps.cpp.
use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
use cosmolkit_search::{
    AtomQueryPredicate as A, BondQueryPredicate as B, QueryAtom, QueryBond, QueryNode as Q,
    is_complex_atom_query, is_complex_bond_query,
};
use cosmolkit_types::{BondOrder as O, Element};

fn atom(query: Q<A>) -> QueryAtom {
    QueryAtom::from_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        query,
    )
}
fn bond(query: Q<B>) -> QueryBond {
    QueryBond::from_parts(
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), O::Single),
        ),
        query,
    )
}

#[test]
fn atom_simple_roots_are_distinct_from_standalone_constraints() {
    for query in [
        Q::predicate(A::Any),
        Q::predicate(A::AtomicNumber(6)),
        Q::predicate(A::AtomType {
            atomic_number: 6,
            aromatic: true,
        }),
    ] {
        assert!(!is_complex_atom_query(&atom(query)));
    }
    for query in [
        Q::predicate(A::FormalCharge(0)),
        Q::predicate(A::IsAromatic(true)),
        Q::predicate(A::ExplicitDegree(2)),
        Q::predicate(A::AtomicNumberIn(vec![6, 7])),
        Q::predicate(A::AtomicNumberNotIn(vec![6, 7])),
        Q::not(Q::predicate(A::Any)),
        Q::or(vec![Q::predicate(A::AtomicNumber(6))]),
        Q::xor(vec![Q::predicate(A::AtomicNumber(6))]),
    ] {
        assert!(is_complex_atom_query(&atom(query)));
    }
}

#[test]
fn atom_and_requires_atomic_identity_and_rejects_complex_descendants() {
    for query in [
        Q::and(vec![]),
        Q::and(vec![Q::predicate(A::Any), Q::predicate(A::FormalCharge(0))]),
        Q::and(vec![
            Q::predicate(A::AtomicNumber(6)),
            Q::not(Q::predicate(A::FormalCharge(1))),
        ]),
        Q::and(vec![
            Q::predicate(A::AtomicNumber(6)),
            Q::and(vec![Q::or(vec![])]),
        ]),
        Q::and(vec![
            Q::predicate(A::AtomicNumber(6)),
            Q::predicate(A::AtomicNumberIn(vec![6, 7])),
        ]),
        Q::and(vec![
            Q::predicate(A::AtomicNumber(6)),
            Q::predicate(A::AtomicNumberNotIn(vec![6, 7])),
        ]),
    ] {
        assert!(is_complex_atom_query(&atom(query)));
    }
    for identity in [
        A::AtomicNumber(6),
        A::AtomType {
            atomic_number: 6,
            aromatic: false,
        },
    ] {
        let query = Q::and(vec![
            Q::predicate(A::FormalCharge(0)),
            Q::and(vec![
                Q::predicate(identity),
                Q::predicate(A::ExplicitDegree(2)),
            ]),
        ]);
        assert!(!is_complex_atom_query(&atom(query)));
    }
}

#[test]
fn bond_only_source_simple_order_branches_are_noncomplex() {
    for query in [
        Q::predicate(B::Order(O::Double)),
        Q::predicate(B::OrderIn(vec![O::Single, O::Aromatic])),
        Q::or(vec![
            Q::predicate(B::Order(O::Single)),
            Q::predicate(B::Order(O::Aromatic)),
        ]),
        Q::or(vec![
            Q::predicate(B::Order(O::Aromatic)),
            Q::predicate(B::Order(O::Single)),
        ]),
        // Source checks each child independently; it does not require distinct orders.
        Q::or(vec![
            Q::predicate(B::Order(O::Single)),
            Q::predicate(B::Order(O::Single)),
        ]),
    ] {
        assert!(!is_complex_bond_query(&bond(query)));
    }
    for query in [
        Q::predicate(B::Any),
        Q::predicate(B::IsInRing(true)),
        Q::predicate(B::OrderIn(vec![O::Double, O::Aromatic])),
        Q::not(Q::predicate(B::Order(O::Single))),
        Q::and(vec![Q::predicate(B::Order(O::Single))]),
        Q::xor(vec![Q::predicate(B::Order(O::Single))]),
        Q::or(vec![]),
        Q::or(vec![Q::predicate(B::Order(O::Single))]),
        Q::or(vec![
            Q::predicate(B::Order(O::Single)),
            Q::predicate(B::Order(O::Double)),
        ]),
        Q::or(vec![
            Q::predicate(B::Order(O::Single)),
            Q::not(Q::predicate(B::Order(O::Aromatic))),
        ]),
        Q::or(vec![Q::predicate(B::Order(O::Single)); 3]),
    ] {
        assert!(is_complex_bond_query(&bond(query)));
    }
}

#[test]
fn ordinary_carriers_keep_source_has_query_false() {
    let ordinary_atom = QueryAtom::from_carrier_parts(
        Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
        Q::not(Q::predicate(A::AtomicNumber(6))),
    );
    let ordinary_bond = QueryBond::from_carrier_parts(
        Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), O::Single),
        ),
        Q::predicate(B::Any),
    );
    assert!(!is_complex_atom_query(&ordinary_atom));
    assert!(!is_complex_bond_query(&ordinary_bond));
}
