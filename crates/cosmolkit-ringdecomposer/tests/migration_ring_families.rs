use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};

use cosmolkit_ringdecomposer::{Graph, RingDecomposerError, RingDecomposition};

#[derive(Debug)]
struct ExpectedCycle {
    declared_size: usize,
    edges: Vec<(usize, usize)>,
}

#[derive(Debug)]
struct Fixture {
    node_count: usize,
    edge_count: usize,
    relevant_cycle_count: usize,
    edges: Vec<(usize, usize)>,
    cycles: Vec<ExpectedCycle>,
}

fn corpus_path(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../testdata/ringdecomposer/corpus/v1.1.3_rdkit")
        .join(name)
}

fn parse_usize(field: Option<&str>, context: &str) -> usize {
    field
        .unwrap_or_else(|| panic!("missing {context}"))
        .parse()
        .unwrap_or_else(|error| panic!("invalid {context}: {error}"))
}

fn zero_based(field: Option<&str>, context: &str) -> usize {
    parse_usize(field, context)
        .checked_sub(1)
        .unwrap_or_else(|| panic!("{context} must be one-based"))
}

fn read_fixture(name: &str) -> Fixture {
    let text = fs::read_to_string(corpus_path(name))
        .unwrap_or_else(|error| panic!("cannot read {name}: {error}"));
    let mut header = None;
    let mut edges = Vec::new();
    let mut cycles = BTreeMap::<usize, ExpectedCycle>::new();

    for (line_number, line) in text.lines().enumerate() {
        let mut fields = line.split_whitespace();
        let Some(record) = fields.next() else {
            continue;
        };
        match record {
            "c" => {}
            "p" => {
                assert_eq!(fields.next(), Some("edge"), "{name}:{}", line_number + 1);
                let node_count = parse_usize(fields.next(), "node count");
                let edge_count = parse_usize(fields.next(), "edge count");
                let cycle_count = parse_usize(fields.next(), "relevant-cycle count");
                assert!(
                    fields.next().is_none(),
                    "{name}:{} trailing header fields",
                    line_number + 1
                );
                assert!(
                    header
                        .replace((node_count, edge_count, cycle_count))
                        .is_none(),
                    "{name}: duplicate header"
                );
            }
            "e" => {
                let from = zero_based(fields.next(), "edge start");
                let to = zero_based(fields.next(), "edge end");
                assert!(
                    fields.next().is_none(),
                    "{name}:{} trailing edge fields",
                    line_number + 1
                );
                edges.push((from, to));
            }
            "r" => {
                let cycle_id = parse_usize(fields.next(), "cycle id");
                let declared_size = parse_usize(fields.next(), "cycle size");
                let from = zero_based(fields.next(), "cycle edge start");
                let to = zero_based(fields.next(), "cycle edge end");
                assert!(
                    fields.next().is_none(),
                    "{name}:{} trailing cycle fields",
                    line_number + 1
                );
                let cycle = cycles.entry(cycle_id).or_insert_with(|| ExpectedCycle {
                    declared_size,
                    edges: Vec::new(),
                });
                assert_eq!(
                    cycle.declared_size, declared_size,
                    "{name}: inconsistent cycle size for {cycle_id}"
                );
                cycle.edges.push((from, to));
            }
            other => panic!("{name}:{} unknown record {other}", line_number + 1),
        }
    }

    let (node_count, edge_count, relevant_cycle_count) =
        header.unwrap_or_else(|| panic!("{name}: missing header"));
    assert_eq!(edges.len(), edge_count, "{name}: edge-row count");
    assert_eq!(
        cycles.len(),
        relevant_cycle_count,
        "{name}: cycle-row groups"
    );
    Fixture {
        node_count,
        edge_count,
        relevant_cycle_count,
        edges,
        cycles: cycles.into_values().collect(),
    }
}

fn build_graph(fixture: &Fixture) -> Graph {
    let mut graph = Graph::new(fixture.node_count);
    for &(from, to) in &fixture.edges {
        graph.add_undirected_edge(from, to).unwrap();
    }
    graph
}

fn validate_expected_cycle(
    name: &str,
    cycle_index: usize,
    cycle: &ExpectedCycle,
    graph: &Graph,
) -> (BTreeSet<usize>, BTreeSet<usize>) {
    assert_eq!(
        cycle.edges.len(),
        cycle.declared_size,
        "{name}: cycle {} row count",
        cycle_index + 1
    );
    let mut edge_ids = BTreeSet::new();
    let mut degree = BTreeMap::<usize, usize>::new();
    let mut adjacency = BTreeMap::<usize, Vec<usize>>::new();
    for &(from, to) in &cycle.edges {
        assert!(from < graph.node_count() && to < graph.node_count());
        let edge = graph.edge_id(from, to).unwrap().index();
        assert!(
            edge_ids.insert(edge),
            "{name}: duplicate edge in cycle {}",
            cycle_index + 1
        );
        *degree.entry(from).or_default() += 1;
        *degree.entry(to).or_default() += 1;
        adjacency.entry(from).or_default().push(to);
        adjacency.entry(to).or_default().push(from);
    }
    assert!(
        degree.values().all(|&value| value == 2),
        "{name}: cycle {} is not degree two",
        cycle_index + 1
    );
    assert_eq!(
        degree.len(),
        cycle.declared_size,
        "{name}: cycle {} vertex count",
        cycle_index + 1
    );
    let start = *degree.keys().next().expect("nonempty source cycle");
    let mut seen = BTreeSet::from([start]);
    let mut stack = vec![start];
    while let Some(node) = stack.pop() {
        for &neighbor in &adjacency[&node] {
            if seen.insert(neighbor) {
                stack.push(neighbor);
            }
        }
    }
    assert_eq!(
        seen.len(),
        degree.len(),
        "{name}: cycle {} is disconnected",
        cycle_index + 1
    );
    (edge_ids, degree.into_keys().collect())
}

fn validate_fixture(name: &str) {
    let fixture = read_fixture(name);
    let graph = build_graph(&fixture);
    assert_eq!(graph.node_count(), fixture.node_count);
    assert_eq!(graph.edge_count(), fixture.edge_count);

    let first = RingDecomposition::calculate(graph.clone()).unwrap();
    let second = RingDecomposition::calculate(graph.clone()).unwrap();
    assert_eq!(first, second, "{name}: decomposition is not deterministic");
    assert_eq!(
        first.relevant_cycle_count(),
        fixture.relevant_cycle_count as f64,
        "{name}: relevant-cycle count"
    );

    let urf_edges = first
        .urfs()
        .iter()
        .map(|urf| {
            urf.edges()
                .iter()
                .map(|edge| edge.index())
                .collect::<BTreeSet<_>>()
        })
        .collect::<Vec<_>>();
    let urf_nodes = first
        .urfs()
        .iter()
        .map(|urf| urf.nodes().iter().copied().collect::<BTreeSet<_>>())
        .collect::<Vec<_>>();
    for urf in first.urfs() {
        assert!(
            urf.edges()
                .windows(2)
                .all(|pair| pair[0].index() < pair[1].index()),
            "{name}: URF edge ids are not in stable input order"
        );
        assert_eq!(
            urf.edges().len(),
            urf.edges().iter().collect::<BTreeSet<_>>().len(),
            "{name}: duplicate URF edge id"
        );
        assert_eq!(
            urf.nodes().len(),
            urf.nodes().iter().collect::<BTreeSet<_>>().len(),
            "{name}: duplicate URF node id"
        );
    }

    let mut expected_cycles = Vec::with_capacity(fixture.cycles.len());
    let mut all_expected_edges = BTreeSet::new();
    let mut all_expected_nodes = BTreeSet::new();
    for (cycle_index, cycle) in fixture.cycles.iter().enumerate() {
        let (cycle_edges, cycle_nodes) = validate_expected_cycle(name, cycle_index, cycle, &graph);
        assert!(
            (0..first.urf_count()).any(|index| {
                cycle_edges.is_subset(&urf_edges[index]) && cycle_nodes.is_subset(&urf_nodes[index])
            }),
            "{name}: source cycle {} is not contained by any returned URF",
            cycle_index + 1
        );
        all_expected_edges.extend(cycle_edges.iter().copied());
        all_expected_nodes.extend(cycle_nodes.iter().copied());
        expected_cycles.push((cycle_edges, cycle_nodes));
    }
    for index in 0..first.urf_count() {
        assert!(
            expected_cycles.iter().any(|(cycle_edges, cycle_nodes)| {
                cycle_edges.is_subset(&urf_edges[index]) && cycle_nodes.is_subset(&urf_nodes[index])
            }),
            "{name}: returned URF {index} contains no complete source cycle"
        );
    }
    let all_returned_edges = urf_edges.iter().flatten().copied().collect::<BTreeSet<_>>();
    let all_returned_nodes = urf_nodes.iter().flatten().copied().collect::<BTreeSet<_>>();
    assert_eq!(
        all_expected_edges, all_returned_edges,
        "{name}: global URF edge membership differs from all source cycles"
    );
    assert_eq!(
        all_expected_nodes, all_returned_nodes,
        "{name}: global URF node membership differs from all source cycles"
    );
}

fn graph(node_count: usize, edges: &[(usize, usize)]) -> Graph {
    let mut graph = Graph::new(node_count);
    for &(from, to) in edges {
        graph.add_undirected_edge(from, to).unwrap();
    }
    graph
}

fn edge_indices(decomposition: &RingDecomposition) -> Vec<Vec<usize>> {
    decomposition
        .urfs()
        .iter()
        .map(|urf| urf.edges().iter().map(|edge| edge.index()).collect())
        .collect()
}

#[test]
fn empty_singleton_and_forest_follow_rdl_result_boundaries() {
    assert_eq!(
        RingDecomposition::calculate(Graph::new(0)),
        Err(RingDecomposerError::EmptyGraph)
    );
    for graph in [Graph::new(1), graph(5, &[(0, 1), (1, 2), (3, 4)])] {
        let result = RingDecomposition::calculate(graph).unwrap();
        assert_eq!(result.urf_count(), 0);
        assert_eq!(result.relevant_cycle_count(), 0.0);
        assert!(result.urfs().is_empty());
    }
}

#[test]
fn odd_and_even_single_cycles_have_complete_stable_membership() {
    for (graph, nodes, edges) in [
        (
            graph(3, &[(2, 0), (0, 1), (1, 2)]),
            vec![0, 2, 1],
            vec![0, 1, 2],
        ),
        (
            graph(4, &[(3, 0), (0, 1), (1, 2), (2, 3)]),
            vec![0, 3, 1, 2],
            vec![0, 1, 2, 3],
        ),
    ] {
        let result = RingDecomposition::calculate(graph).unwrap();
        assert_eq!(result.urf_count(), 1);
        assert_eq!(result.relevant_cycle_count(), 1.0);
        assert_eq!(result.urfs()[0].nodes(), nodes);
        assert_eq!(
            result.urfs()[0]
                .edges()
                .iter()
                .map(|edge| edge.index())
                .collect::<Vec<_>>(),
            edges
        );
    }
}

#[test]
fn disconnected_and_articulation_separated_cycles_keep_component_order() {
    let disconnected =
        RingDecomposition::calculate(graph(7, &[(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)]))
            .unwrap();
    assert_eq!(disconnected.urf_count(), 2);
    assert_eq!(disconnected.relevant_cycle_count(), 2.0);
    assert_eq!(
        edge_indices(&disconnected),
        vec![vec![0, 1, 2], vec![3, 4, 5]]
    );

    let articulation =
        RingDecomposition::calculate(graph(5, &[(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)]))
            .unwrap();
    assert_eq!(articulation.urf_count(), 2);
    assert_eq!(articulation.relevant_cycle_count(), 2.0);
    assert_eq!(
        edge_indices(&articulation),
        vec![vec![3, 4, 5], vec![0, 1, 2]]
    );
}

#[test]
fn fused_spiro_and_bridged_graphs_close_relation_and_path_count_branches() {
    let fused =
        RingDecomposition::calculate(graph(4, &[(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)])).unwrap();
    assert_eq!((fused.urf_count(), fused.relevant_cycle_count()), (2, 2.0));
    assert_eq!(edge_indices(&fused), vec![vec![0, 1, 2], vec![1, 3, 4]]);

    let spiro =
        RingDecomposition::calculate(graph(5, &[(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]))
            .unwrap();
    assert_eq!((spiro.urf_count(), spiro.relevant_cycle_count()), (2, 2.0));

    let bridged =
        RingDecomposition::calculate(graph(5, &[(0, 1), (1, 4), (0, 2), (2, 4), (0, 3), (3, 4)]))
            .unwrap();
    assert_eq!(
        (bridged.urf_count(), bridged.relevant_cycle_count()),
        (3, 3.0)
    );
    assert_eq!(
        edge_indices(&bridged),
        vec![vec![0, 1, 2, 3], vec![0, 1, 4, 5], vec![2, 3, 4, 5]]
    );
}

macro_rules! corpus_cases {
    ($($test:ident => $file:literal),+ $(,)?) => {$(
        #[test]
        fn $test() {
            validate_fixture($file);
        }
    )+};
}

corpus_cases! {
    source_molecule_00 => "molecule_00_0.dimacs",
    source_molecule_01 => "molecule_01_0.dimacs",
    source_molecule_02 => "molecule_02_0.dimacs",
    source_molecule_03 => "molecule_03_0.dimacs",
    source_molecule_04 => "molecule_04_0.dimacs",
    source_molecule_05 => "molecule_05_0.dimacs",
    source_molecule_06 => "molecule_06_0.dimacs",
    source_molecule_07 => "molecule_07_0.dimacs",
    source_molecule_08 => "molecule_08_0.dimacs",
    source_molecule_09 => "molecule_09_0.dimacs",
    source_molecule_10 => "molecule_10_0.dimacs",
    source_molecule_11 => "molecule_11_0.dimacs",
    source_molecule_12 => "molecule_12_0.dimacs",
    source_molecule_13 => "molecule_13_0.dimacs",
    source_molecule_14 => "molecule_14_0.dimacs",
    source_molecule_15 => "molecule_15_0.dimacs",
    source_molecule_16 => "molecule_16_0.dimacs",
    source_molecule_17 => "molecule_17_0.dimacs",
    source_molecule_18 => "molecule_18_0.dimacs",
    source_molecule_19 => "molecule_19_0.dimacs",
    source_molecule_20 => "molecule_20_0.dimacs",
    source_molecule_21 => "molecule_21_0.dimacs",
    source_molecule_22 => "molecule_22_0.dimacs",
    source_molecule_23 => "molecule_23_0.dimacs",
    source_molecule_24 => "molecule_24_0.dimacs",
    source_molecule_25 => "molecule_25_0.dimacs",
}
