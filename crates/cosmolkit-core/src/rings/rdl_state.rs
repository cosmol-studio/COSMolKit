//! Ownership conversion at the CK/independent RDL boundary. The large source
//! buffers move once; later cycle queries borrow the CK snapshot directly.

use cosmolkit_model::ring_decomposition::{
    RingComponent, RingCycleFamily, RingDecomposition, RingEdge, RingFamily, RingGraph,
};
use cosmolkit_ringdecomposer::persistence::CycleFamilyView;
use cosmolkit_ringdecomposer::{CycleComponentAccess, EdgeId};

fn graph(graph: cosmolkit_ringdecomposer::Graph) -> RingGraph {
    let graph = graph.into_state();
    RingGraph {
        node_count: graph.node_count,
        adjacency: graph
            .adjacency
            .into_iter()
            .map(|row| {
                row.into_iter()
                    .map(|(node, edge)| (node, edge.index()))
                    .collect()
            })
            .collect(),
        edges: graph
            .edges
            .into_iter()
            .map(|e| RingEdge {
                id: e.id().index(),
                from: e.from(),
                to: e.to(),
            })
            .collect(),
        edge_ids_by_endpoints: graph
            .edge_ids_by_endpoints
            .into_iter()
            .map(|(ends, id)| (ends, id.index()))
            .collect(),
    }
}

/// Consume a source result before introducing Arc sharing. Numeric indices
/// are lowered to CK-owned integers; matrices, DAG rows, prototypes and family
/// memberships transfer their existing buffers, without source recomputation.
pub fn into_model(decomposition: cosmolkit_ringdecomposer::RingDecomposition) -> RingDecomposition {
    let d = decomposition.into_state();
    RingDecomposition::new(
        graph(d.graph),
        d.urfs
            .into_iter()
            .map(|u| {
                let (nodes, edges) = u.into_parts();
                RingFamily {
                    nodes,
                    edges: edges.into_iter().map(EdgeId::index).collect(),
                }
            })
            .collect(),
        d.relevant_cycle_count,
        d.components
            .into_iter()
            .map(|c| RingComponent {
                graph: graph(c.graph),
                original_nodes: c.original_nodes,
                original_edges: c.original_edges.into_iter().map(EdgeId::index).collect(),
                predecessor: c.predecessor,
                distance: c.distance,
                reachable: c.reachable,
                directed_paths: c.directed_paths,
                cycle_families: c
                    .cycle_families
                    .into_iter()
                    .map(|f| RingCycleFamily {
                        weight: f.weight,
                        root: f.root,
                        p: f.p,
                        q: f.q,
                        x: f.x,
                        prototype: f.prototype,
                        relevant: f.relevant,
                    })
                    .collect(),
                prototype_counts: c.prototype_counts,
                relations: c.relations,
                urfs: c.urfs,
            })
            .collect(),
    )
    .expect("the source engine only constructs finite nonnegative integral cycle counts")
}

#[derive(Clone, Copy)]
struct Component<'a>(&'a RingComponent);

impl<'a> CycleComponentAccess<'a> for Component<'a> {
    fn node_count(self) -> usize {
        self.0.graph.node_count
    }
    fn edge_count(self) -> usize {
        self.0.graph.edges.len()
    }
    fn edge_id(self, from: usize, to: usize) -> EdgeId {
        let ends = (from.min(to), from.max(to));
        EdgeId::new(self.0.graph.edge_ids_by_endpoints[&ends])
    }
    fn original_edge(self, edge: usize) -> EdgeId {
        EdgeId::new(self.0.original_edges[edge])
    }
    fn urfs(self) -> &'a [Vec<usize>] {
        &self.0.urfs
    }
    fn family(self, family: usize) -> CycleFamilyView<'a> {
        let f = &self.0.cycle_families[family];
        CycleFamilyView {
            weight: f.weight,
            root: f.root,
            p: f.p,
            q: f.q,
            x: f.x,
            prototype: &f.prototype,
            relevant: f.relevant,
        }
    }
    fn path_neighbors(self, root: usize, node: usize) -> &'a [usize] {
        &self.0.directed_paths[root][node]
    }
}

pub fn smallest_ring_basis(d: &RingDecomposition) -> Vec<Vec<EdgeId>> {
    cosmolkit_ringdecomposer::smallest_ring_basis(d.components.iter().map(Component))
}

pub fn relevant_cycles(d: &RingDecomposition) -> impl Iterator<Item = Vec<EdgeId>> + '_ {
    cosmolkit_ringdecomposer::relevant_cycles(d.components.iter().map(Component))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn owned_lowering_moves_buffers_and_borrowed_queries_preserve_source_order() {
        let mut g = cosmolkit_ringdecomposer::Graph::new(5);
        for (a, b) in [(0, 1), (1, 2), (2, 0), (1, 3), (3, 4), (4, 2)] {
            g.add_undirected_edge(a, b).unwrap();
        }
        let source = cosmolkit_ringdecomposer::RingDecomposition::calculate(g).unwrap();
        let expected_basis = source.smallest_ring_basis();
        let expected_cycles = source.relevant_cycles().collect::<Vec<_>>();
        let c = source.persisted_components().next().unwrap();
        let pointers = (
            c.original_nodes.as_ptr(),
            c.predecessor.as_ptr(),
            c.distance.as_ptr(),
            c.reachable.as_ptr(),
            c.relations.as_ptr(),
            c.urfs.as_ptr(),
            c.directed_paths().next().unwrap().as_ptr(),
            c.cycle_families().next().unwrap().prototype.as_ptr(),
        );
        let detached = into_model(source);
        let c = &detached.components[0];
        assert_eq!(
            pointers,
            (
                c.original_nodes.as_ptr(),
                c.predecessor.as_ptr(),
                c.distance.as_ptr(),
                c.reachable.as_ptr(),
                c.relations.as_ptr(),
                c.urfs.as_ptr(),
                c.directed_paths[0].as_ptr(),
                c.cycle_families[0].prototype.as_ptr()
            )
        );
        let shared = std::sync::Arc::new(detached);
        let other = std::sync::Arc::clone(&shared);
        for _ in 0..2 {
            assert_eq!(smallest_ring_basis(&shared), expected_basis);
            assert_eq!(
                relevant_cycles(&shared).collect::<Vec<_>>(),
                expected_cycles
            );
        }
        assert!(std::sync::Arc::ptr_eq(&shared, &other));
        assert_eq!(pointers.1, shared.components[0].predecessor.as_ptr());
    }
}
