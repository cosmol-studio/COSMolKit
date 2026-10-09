//! Source path filtering retains real query atomic numbers and borrowed identity.
use cosmolkit_core::{
    GraphPath, SubgraphSearchParams, query_atom_paths_in_range, query_bond_paths_in_range,
    query_subgraphs_in_range,
};
use cosmolkit_model::BondOrder;
use cosmolkit_model::{
    AtomId, AtomQueryPredicate, BondId, BondSpec, QueryAtom, QueryAtomIdentity, QueryBond,
    QueryGraph, QueryNode,
};

fn graph() -> QueryGraph {
    let atoms = [0, 6, 1, 6]
        .into_iter()
        .enumerate()
        .map(|(i, z)| {
            QueryAtom::from_identity_parts(
                AtomId::new(i),
                QueryAtomIdentity::from_atomic_number(z),
                QueryNode::predicate(AtomQueryPredicate::Any),
            )
        })
        .collect();
    let bonds = [(0, 1), (1, 2), (1, 3)]
        .into_iter()
        .enumerate()
        .map(|(i, (a, b))| {
            QueryBond::new(
                BondId::new(i),
                BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
            )
        })
        .collect();
    QueryGraph::from_parts(atoms, bonds, std::collections::BTreeMap::<cosmolkit_model::PropertyText, cosmolkit_model::PropertyValue>::new(), vec![], vec![], vec![]).unwrap()
}
#[test]
fn zero_query_identity_is_retained_and_only_actual_hydrogen_is_filtered() {
    let query = graph();
    let before = query.clone();
    let params = SubgraphSearchParams::default();
    let paths = query_bond_paths_in_range(&query, 1, 2, &params).unwrap();
    assert_eq!(
        paths[&1],
        [
            GraphPath::Bonds(vec![BondId::new(0)]),
            GraphPath::Bonds(vec![BondId::new(2)]),
        ]
    );
    assert_eq!(
        paths[&2],
        [GraphPath::Bonds(vec![BondId::new(0), BondId::new(2)])]
    );
    let subgraphs = query_subgraphs_in_range(&query, 1, 2, &params).unwrap();
    assert_eq!(subgraphs[&1], [vec![BondId::new(0)], vec![BondId::new(2)]]);
    assert_eq!(subgraphs[&2], [vec![BondId::new(0), BondId::new(2)]]);
    assert_eq!(
        query_subgraphs_in_range(
            &query,
            1,
            1,
            &SubgraphSearchParams {
                use_hydrogens: true,
                ..params
            }
        )
        .unwrap()[&1]
            .len(),
        3
    );
    assert_eq!(query.atoms()[0].atomic_number(), 0);
    assert_eq!(query, before);
}
#[test]
fn query_roots_and_source_ranges_use_the_shared_enumerator() {
    let query = graph();
    let rooted = SubgraphSearchParams {
        rooted_at_atom: Some(AtomId::new(3)),
        ..SubgraphSearchParams::default()
    };
    assert_eq!(
        query_subgraphs_in_range(&query, 1, 1, &rooted).unwrap()[&1],
        [vec![BondId::new(2)]]
    );
    assert!(query_bond_paths_in_range(&query, 3, 2, &rooted).is_err());
}

#[test]
fn layered_atom_paths_keep_atom_ids_lengths_and_query_hydrogen_filtering() {
    let query = graph();
    let before = query.clone();
    let params = SubgraphSearchParams::default();
    let paths = query_atom_paths_in_range(&query, 2, 3, &params).unwrap();
    assert_eq!(
        paths[&2],
        [
            GraphPath::Atoms(vec![AtomId::new(0), AtomId::new(1)]),
            GraphPath::Atoms(vec![AtomId::new(1), AtomId::new(3)]),
        ]
    );
    assert_eq!(
        paths[&3],
        [GraphPath::Atoms(vec![
            AtomId::new(0),
            AtomId::new(1),
            AtomId::new(3)
        ])]
    );
    assert!(query_atom_paths_in_range(&query, 3, 2, &params).is_err());
    let with_hs = SubgraphSearchParams {
        use_hydrogens: true,
        ..params
    };
    assert_eq!(
        query_atom_paths_in_range(&query, 2, 2, &with_hs).unwrap()[&2].len(),
        3
    );
    assert_eq!(query, before);
}
#[test]
fn search01_query_mask_keeps_real_hydrogen_and_zero_identity_gates() {
    let query = graph();
    let before = query.clone();
    let mask = [true, false, false, false];
    for use_hydrogens in [false, true] {
        let p = SubgraphSearchParams {
            ignore_atoms: Some(&mask),
            use_hydrogens,
            rooted_at_atom: Some(AtomId::new(1)),
        };
        let expected = if use_hydrogens {
            vec![vec![BondId::new(1)], vec![BondId::new(2)]]
        } else {
            vec![vec![BondId::new(2)]]
        };
        assert_eq!(
            query_subgraphs_in_range(&query, 1, 1, &p).unwrap()[&1],
            expected
        );
        let expected: Vec<_> = expected.into_iter().map(GraphPath::Bonds).collect();
        assert_eq!(
            query_bond_paths_in_range(&query, 1, 1, &p).unwrap()[&1],
            expected
        );
    }
    assert_eq!(query, before);
}
