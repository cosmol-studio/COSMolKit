use std::collections::BTreeMap;

use cosmolkit_model::{Atom, AtomId, AtomSpec, Bond, BondId, BondSpec, QueryGraph, TopologyBlock};
use cosmolkit_search::{SmartsParseParams, compile_query, parse_smarts};
use cosmolkit_types::{BondOrder, Element};

fn topology(elements: &[Element], edge_endpoints: &[(usize, usize)]) -> TopologyBlock {
    let atoms = elements
        .iter()
        .copied()
        .enumerate()
        .map(|(index, element)| Atom::from_spec(AtomId::new(index), AtomSpec::new(element)))
        .collect();
    let bonds = edge_endpoints
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
        .expect("fixed Q-compile target topology is valid")
}

#[test]
fn q65_compiled_graph_handles_empty_and_disconnected_queries() {
    let empty = QueryGraph::from_parts(
        Vec::new(),
        Vec::new(),
        BTreeMap::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    )
    .expect("empty Q65 query graph is valid");
    let compiled_empty = compile_query(&empty).expect("empty query graph compiles");
    assert_eq!(compiled_empty.num_atoms(), 0);
    assert_eq!(compiled_empty.num_bonds(), 0);
    assert!(compiled_empty.atom_order().is_empty());

    let disconnected = parse_smarts("C.O", &SmartsParseParams::default())
        .expect("fixed disconnected Q65 query parses");
    let compiled_disconnected = compile_query(&disconnected).expect("disconnected query compiles");
    assert_eq!(compiled_disconnected.num_atoms(), 2);
    assert_eq!(compiled_disconnected.num_bonds(), 0);
    let mut atom_order = compiled_disconnected.atom_order().to_vec();
    atom_order.sort_unstable();
    assert_eq!(atom_order, [0, 1]);

    let target = topology(&[Element::C, Element::O], &[]);
    assert!(
        compiled_disconnected
            .matches(&target)
            .expect("disconnected query matches detached target")
            .iter()
            .any(|matched| matched.atom_mapping == [0, 1])
    );
}

#[test]
fn q65_compiled_graph_keeps_undirected_endpoints_and_bond_indices() {
    let query = parse_smarts("C-O-N", &SmartsParseParams::default())
        .expect("fixed Q65 linear query parses");
    let compiled = compile_query(&query).expect("linear query compiles");

    // Both target bonds are stored with endpoint order opposite to the query.
    // RDKit's MolGraph is undirected and keeps each edge's bond identity.
    let target = topology(&[Element::C, Element::O, Element::N], &[(1, 0), (2, 1)]);
    let matches = compiled
        .matches(&target)
        .expect("reversed endpoint target is valid");
    assert!(
        matches
            .iter()
            .any(|matched| { matched.atom_mapping == [0, 1, 2] && matched.bond_mapping == [0, 1] })
    );
}

#[test]
fn q66_compiled_order_prioritizes_rare_degrees_and_keeps_equal_degree_ties() {
    let query =
        parse_smarts("C-O-N.C", &SmartsParseParams::default()).expect("fixed Q66 query parses");
    let compiled = compile_query(&query).expect("query compiles");
    let order = compiled.atom_order();

    assert_eq!(order.len(), 4);
    assert_eq!(
        order.first(),
        Some(&1),
        "degree-two node is most frequent-rare"
    );
    assert_eq!(
        order.last(),
        Some(&3),
        "isolated node follows nonzero degree"
    );
    let mut degree_one_tie = order[1..3].to_vec();
    degree_one_tie.sort_unstable();
    assert_eq!(degree_one_tie, [0, 2]);
}

#[test]
fn q67_frequency_order_covers_disconnected_groups_and_source_ties() {
    let query_smarts = "C.C-C.C-C-C.C-C-C.C(-C)(-C)-C.C(-C)(-C)-C";
    let query = parse_smarts(query_smarts, &SmartsParseParams::default())
        .expect("fixed Q67 disconnected query parses");
    let compiled = compile_query(&query).expect("disconnected query compiles");

    assert_eq!(compiled.num_atoms(), 17);
    assert_eq!(compiled.num_bonds(), 11);

    let mut degrees = vec![0; compiled.num_atoms()];
    for bond in compiled.query().bonds() {
        let (begin, end) = bond.endpoints();
        degrees[begin] += 1;
        degrees[end] += 1;
    }

    let order = compiled.atom_order();
    let mut sorted_ids = order.to_vec();
    sorted_ids.sort_unstable();
    assert_eq!(sorted_ids, (0..compiled.num_atoms()).collect::<Vec<_>>());

    let ordered_degrees: Vec<_> = order.iter().map(|&atom| degrees[atom]).collect();
    assert_eq!(&ordered_degrees[..4], &[2, 2, 3, 3]);
    assert!(ordered_degrees[4..16].iter().all(|&degree| degree == 1));
    assert_eq!(ordered_degrees[16], 0);

    // RDKit's comparator does not break ties by node ID; compare the tied
    // source identities as sets while checking their exact rank groups above.
    let mut degree_two_ids = order[..2].to_vec();
    degree_two_ids.sort_unstable();
    assert_eq!(degree_two_ids, [4, 7]);
    let mut degree_three_ids = order[2..4].to_vec();
    degree_three_ids.sort_unstable();
    assert_eq!(degree_three_ids, [9, 13]);
    assert_eq!(order[16], 0);
}

#[test]
fn q68_compiled_plan_owns_query_ids_and_remains_immutable() {
    let compiled = {
        let query =
            parse_smarts("C-O", &SmartsParseParams::default()).expect("fixed Q68 query parses");
        let compiled = compile_query(&query).expect("query compiles");
        assert_eq!(compiled.query(), &query);
        compiled
    };

    assert_eq!(compiled.query().num_atoms(), 2);
    assert_eq!(compiled.query().num_bonds(), 1);
    assert_eq!(compiled.query().atoms()[0].id(), AtomId::new(0));
    assert_eq!(compiled.query().atoms()[1].id(), AtomId::new(1));
    assert_eq!(compiled.query().bonds()[0].id(), BondId::new(0));

    let before_match = compiled.clone();
    let target = topology(&[Element::C, Element::O], &[(0, 1)]);
    assert!(
        !compiled
            .matches(&target)
            .expect("compiled query matches detached topology")
            .is_empty()
    );
    assert_eq!(compiled, before_match);
}

#[test]
fn q69_compiled_recursive_query_resets_between_targets() {
    let query = parse_smarts("[$(C-C)]", &SmartsParseParams::default())
        .expect("fixed Q69 recursive query parses");
    let compiled = compile_query(&query).expect("recursive query compiles");
    let before_matches = compiled.clone();
    let carbon_pair = topology(&[Element::C, Element::C], &[(0, 1)]);
    let carbon_nitrogen_pair = topology(&[Element::C, Element::N], &[(0, 1)]);

    assert!(
        !compiled
            .matches(&carbon_pair)
            .expect("recursive query matches carbon pair")
            .is_empty()
    );
    assert!(
        compiled
            .matches(&carbon_nitrogen_pair)
            .expect("recursive query checks the second target")
            .is_empty()
    );
    assert!(
        !compiled
            .matches(&carbon_pair)
            .expect("recursive cache is recomputed for the original target")
            .is_empty()
    );
    assert_eq!(compiled, before_matches);
}
