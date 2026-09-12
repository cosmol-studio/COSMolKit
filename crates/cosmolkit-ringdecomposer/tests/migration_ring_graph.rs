use cosmolkit_ringdecomposer::{EdgeId, Graph, RingDecomposerError};

fn edge_rows(graph: &Graph) -> Vec<(usize, usize, usize)> {
    graph
        .edges()
        .iter()
        .map(|edge| (edge.id().index(), edge.from(), edge.to()))
        .collect()
}

fn neighbor_row(graph: &Graph, node: usize) -> Vec<(usize, usize)> {
    graph
        .neighbors(node)
        .expect("node is in range")
        .iter()
        .map(|(neighbor, edge)| (*neighbor, edge.index()))
        .collect()
}

#[test]
fn zero_and_single_node_graphs_have_source_shaped_connectivity() {
    let empty = Graph::new(0);
    assert_eq!(empty.node_count(), 0);
    assert_eq!(empty.edge_count(), 0);
    assert!(empty.edges().is_empty());
    assert_eq!(empty.neighbors(0), None);
    assert!(!empty.is_adjacent(0, 0));
    assert!(empty.is_connected());

    let singleton = Graph::new(1);
    assert_eq!(singleton.node_count(), 1);
    assert_eq!(singleton.neighbors(0), Some(&[][..]));
    assert_eq!(singleton.neighbors(1), None);
    assert!(!singleton.is_adjacent(0, 0));
    assert!(singleton.is_connected());
}

#[test]
fn successful_edges_use_dense_input_ids_and_canonical_storage() {
    let mut graph = Graph::new(5);
    assert_eq!(graph.add_undirected_edge(3, 1), Ok(EdgeId::new(0)));
    assert_eq!(graph.add_undirected_edge(0, 3), Ok(EdgeId::new(1)));
    assert_eq!(graph.add_undirected_edge(4, 3), Ok(EdgeId::new(2)));
    assert_eq!(graph.add_undirected_edge(2, 1), Ok(EdgeId::new(3)));

    assert_eq!(
        edge_rows(&graph),
        vec![(0, 1, 3), (1, 0, 3), (2, 3, 4), (3, 1, 2)]
    );
    assert_eq!(neighbor_row(&graph, 0), vec![(3, 1)]);
    assert_eq!(neighbor_row(&graph, 1), vec![(3, 0), (2, 3)]);
    assert_eq!(neighbor_row(&graph, 2), vec![(1, 3)]);
    assert_eq!(neighbor_row(&graph, 3), vec![(1, 0), (0, 1), (4, 2)]);
    assert_eq!(neighbor_row(&graph, 4), vec![(3, 2)]);
}

#[test]
fn every_failed_addition_preserves_graph_and_next_edge_id() {
    let mut graph = Graph::new(3);
    assert_eq!(graph.add_undirected_edge(0, 1), Ok(EdgeId::new(0)));

    let failures = [
        (
            (7, 8),
            RingDecomposerError::NodeOutOfRange {
                node: 7,
                node_count: 3,
            },
        ),
        (
            (2, 9),
            RingDecomposerError::NodeOutOfRange {
                node: 9,
                node_count: 3,
            },
        ),
        ((2, 2), RingDecomposerError::SelfLoop { node: 2 }),
        (
            (0, 1),
            RingDecomposerError::DuplicateEdge { from: 0, to: 1 },
        ),
        (
            (1, 0),
            RingDecomposerError::DuplicateEdge { from: 1, to: 0 },
        ),
    ];

    for ((from, to), expected) in failures {
        let edges_before = edge_rows(&graph);
        let rows_before = (0..graph.node_count())
            .map(|node| neighbor_row(&graph, node))
            .collect::<Vec<_>>();
        assert_eq!(graph.add_undirected_edge(from, to), Err(expected));
        assert_eq!(edge_rows(&graph), edges_before);
        assert_eq!(
            (0..graph.node_count())
                .map(|node| neighbor_row(&graph, node))
                .collect::<Vec<_>>(),
            rows_before
        );
    }

    assert_eq!(graph.add_undirected_edge(1, 2), Ok(EdgeId::new(1)));
}

#[test]
fn lookup_is_orientation_independent_and_errors_keep_attempted_order() {
    let mut graph = Graph::new(3);
    graph.add_undirected_edge(2, 0).unwrap();

    assert_eq!(graph.edge_id(0, 2), Ok(EdgeId::new(0)));
    assert_eq!(graph.edge_id(2, 0), Ok(EdgeId::new(0)));
    assert_eq!(
        graph.edge_id(1, 2),
        Err(RingDecomposerError::EdgeNotFound { from: 1, to: 2 })
    );
    assert_eq!(
        graph.edge_id(9, 0),
        Err(RingDecomposerError::EdgeNotFound { from: 9, to: 0 })
    );
    assert_eq!(graph.neighbors(3), None);
    assert!(!graph.is_adjacent(9, 0));
    assert!(!graph.is_adjacent(0, 9));
}

#[test]
fn chain_branch_ring_and_disconnected_shapes_preserve_complete_rows() {
    let mut chain = Graph::new(4);
    chain.add_undirected_edge(0, 1).unwrap();
    chain.add_undirected_edge(1, 2).unwrap();
    chain.add_undirected_edge(2, 3).unwrap();
    assert!(chain.is_connected());
    assert_eq!(neighbor_row(&chain, 0), vec![(1, 0)]);
    assert_eq!(neighbor_row(&chain, 1), vec![(0, 0), (2, 1)]);
    assert_eq!(neighbor_row(&chain, 2), vec![(1, 1), (3, 2)]);
    assert_eq!(neighbor_row(&chain, 3), vec![(2, 2)]);

    let mut branch = Graph::new(4);
    branch.add_undirected_edge(0, 2).unwrap();
    branch.add_undirected_edge(0, 1).unwrap();
    branch.add_undirected_edge(0, 3).unwrap();
    assert!(branch.is_connected());
    assert_eq!(neighbor_row(&branch, 0), vec![(2, 0), (1, 1), (3, 2)]);
    assert_eq!(neighbor_row(&branch, 1), vec![(0, 1)]);
    assert_eq!(neighbor_row(&branch, 2), vec![(0, 0)]);
    assert_eq!(neighbor_row(&branch, 3), vec![(0, 2)]);

    let mut ring = Graph::new(3);
    ring.add_undirected_edge(2, 0).unwrap();
    ring.add_undirected_edge(1, 2).unwrap();
    ring.add_undirected_edge(0, 1).unwrap();
    assert!(ring.is_connected());
    assert_eq!(neighbor_row(&ring, 0), vec![(2, 0), (1, 2)]);
    assert_eq!(neighbor_row(&ring, 1), vec![(2, 1), (0, 2)]);
    assert_eq!(neighbor_row(&ring, 2), vec![(0, 0), (1, 1)]);

    let mut disconnected = Graph::new(4);
    disconnected.add_undirected_edge(0, 1).unwrap();
    disconnected.add_undirected_edge(2, 3).unwrap();
    assert!(!disconnected.is_connected());
    assert_eq!(neighbor_row(&disconnected, 0), vec![(1, 0)]);
    assert_eq!(neighbor_row(&disconnected, 1), vec![(0, 0)]);
    assert_eq!(neighbor_row(&disconnected, 2), vec![(3, 1)]);
    assert_eq!(neighbor_row(&disconnected, 3), vec![(2, 1)]);
}

#[test]
fn two_cycles_sharing_an_articulation_keep_global_input_order() {
    let mut graph = Graph::new(5);
    for (from, to) in [(2, 0), (0, 1), (1, 2), (4, 2), (2, 3), (3, 4)] {
        graph.add_undirected_edge(from, to).unwrap();
    }

    assert!(graph.is_connected());
    assert_eq!(
        edge_rows(&graph),
        vec![
            (0, 0, 2),
            (1, 0, 1),
            (2, 1, 2),
            (3, 2, 4),
            (4, 2, 3),
            (5, 3, 4),
        ]
    );
    assert_eq!(neighbor_row(&graph, 0), vec![(2, 0), (1, 1)]);
    assert_eq!(neighbor_row(&graph, 1), vec![(0, 1), (2, 2)]);
    assert_eq!(
        neighbor_row(&graph, 2),
        vec![(0, 0), (1, 2), (4, 3), (3, 4)]
    );
    assert_eq!(neighbor_row(&graph, 3), vec![(2, 4), (4, 5)]);
    assert_eq!(neighbor_row(&graph, 4), vec![(2, 3), (3, 5)]);
}
