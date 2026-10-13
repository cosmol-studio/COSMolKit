//! Detached persistence boundary. Restoring state does not run ring finding.
use super::*;

pub struct CycleFamilyView<'a> {
    pub weight: usize,
    pub root: usize,
    pub p: usize,
    pub q: usize,
    pub x: Option<usize>,
    pub prototype: &'a [bool],
    pub relevant: bool,
}

pub struct ComponentView<'a> {
    pub graph: &'a Graph,
    pub original_nodes: &'a [usize],
    pub original_edges: &'a [EdgeId],
    pub predecessor: &'a [Vec<Option<usize>>],
    pub distance: &'a [Vec<Option<usize>>],
    pub reachable: &'a [Vec<bool>],
    pub prototype_counts: &'a [usize],
    pub relations: &'a [Vec<Vec<bool>>],
    pub urfs: &'a [Vec<usize>],
    data: &'a DecomposedComponent,
}

impl ComponentView<'_> {
    pub fn directed_paths(&self) -> impl Iterator<Item = &[Vec<usize>]> {
        self.data
            .shortest_paths
            .directed_paths
            .iter()
            .map(|p| p.adjacency.as_slice())
    }

    pub fn cycle_families(&self) -> impl Iterator<Item = CycleFamilyView<'_>> {
        self.data
            .cycle_families
            .families
            .iter()
            .map(|f| CycleFamilyView {
                weight: f.weight,
                root: f.r,
                p: f.p,
                q: f.q,
                x: f.x,
                prototype: &f.prototype,
                relevant: f.relevant,
            })
    }
}

pub struct CycleFamilyState {
    pub weight: usize,
    pub root: usize,
    pub p: usize,
    pub q: usize,
    pub x: Option<usize>,
    pub prototype: Vec<bool>,
    pub relevant: bool,
}

pub struct ComponentState {
    pub graph: Graph,
    pub original_nodes: Vec<usize>,
    pub original_edges: Vec<EdgeId>,
    pub predecessor: Vec<Vec<Option<usize>>>,
    pub distance: Vec<Vec<Option<usize>>>,
    pub reachable: Vec<Vec<bool>>,
    pub directed_paths: Vec<Vec<Vec<usize>>>,
    pub cycle_families: Vec<CycleFamilyState>,
    pub prototype_counts: Vec<usize>,
    pub relations: Vec<Vec<Vec<bool>>>,
    pub urfs: Vec<Vec<usize>>,
}

impl RingDecomposition {
    pub fn persisted_components(&self) -> impl Iterator<Item = ComponentView<'_>> {
        self.components.iter().map(|data| ComponentView {
            graph: &data.component.graph,
            original_nodes: &data.component.original_nodes,
            original_edges: &data.component.original_edges,
            predecessor: &data.shortest_paths.predecessor,
            distance: &data.shortest_paths.distance,
            reachable: &data.shortest_paths.reachable,
            prototype_counts: &data.urf_info.nof_protos,
            relations: &data.urf_info.relations,
            urfs: &data.urf_info.urfs,
            data,
        })
    }

    pub fn from_persisted_state(
        graph: Graph,
        urfs: Vec<(Vec<usize>, Vec<EdgeId>)>,
        relevant_cycle_count: f64,
        components: Vec<ComponentState>,
    ) -> Result<Self, &'static str> {
        if graph.node_count() == 0
            || !relevant_cycle_count.is_finite()
            || relevant_cycle_count < 0.0
            || relevant_cycle_count.fract() != 0.0
        {
            return Err("invalid RDL graph or relevant-cycle count");
        }
        if urfs.iter().any(|(nodes, edges)| {
            nodes.iter().any(|&n| n >= graph.node_count())
                || edges.iter().any(|e| e.index() >= graph.edge_count())
        }) {
            return Err("RDL family index out of range");
        }
        if components.iter().map(|c| c.urfs.len()).sum::<usize>() != urfs.len() {
            return Err("RDL component/global family count mismatch");
        }
        let components = components
            .into_iter()
            .map(|c| c.restore(&graph))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            graph,
            urfs: urfs
                .into_iter()
                .map(|(nodes, edges)| UniqueRingFamily { nodes, edges })
                .collect(),
            relevant_cycle_count,
            components,
        })
    }
}

fn square<T>(matrix: &[Vec<T>], n: usize) -> bool {
    matrix.len() == n && matrix.iter().all(|row| row.len() == n)
}

impl ComponentState {
    fn restore(self, parent: &Graph) -> Result<DecomposedComponent, &'static str> {
        let n = self.graph.node_count();
        let m = self.graph.edge_count();
        if n == 0 || m < n || self.original_nodes.len() != n || self.original_edges.len() != m {
            return Err("invalid RDL cyclic component dimensions");
        }
        // Track only this component's IDs. Allocating parent-sized masks for
        // every BCC would make many disconnected small rings quadratic.
        let mut nodes = std::collections::BTreeSet::new();
        for &node in &self.original_nodes {
            if node >= parent.node_count() {
                return Err("RDL component node out of range");
            }
            if !nodes.insert(node) {
                return Err("duplicate RDL component node");
            }
        }
        let mut edges = std::collections::BTreeSet::new();
        for (edge, original) in self.graph.edges().iter().zip(&self.original_edges) {
            if original.index() >= parent.edge_count() {
                return Err("RDL component edge out of range");
            }
            if !edges.insert(original.index()) {
                return Err("duplicate RDL component edge");
            }
            let original = parent.edges()[original.index()];
            if canonical_edge_key(
                self.original_nodes[edge.from()],
                self.original_nodes[edge.to()],
            ) != canonical_edge_key(original.from(), original.to())
            {
                return Err("RDL component edge mapping mismatch");
            }
        }
        if !square(&self.predecessor, n)
            || !square(&self.distance, n)
            || !square(&self.reachable, n)
            || self.directed_paths.len() != n
        {
            return Err("invalid RDL shortest-path dimensions");
        }
        for root in 0..n {
            if self.predecessor[root][root] != Some(root)
                || self.distance[root][root] != Some(0)
                || self.predecessor[root].iter().flatten().any(|&v| v >= n)
                || self.distance[root].iter().flatten().any(|&v| v >= n)
                || self.directed_paths[root].len() != n
            {
                return Err("invalid RDL shortest-path row");
            }
            for (from, row) in self.directed_paths[root].iter().enumerate() {
                for &to in row {
                    if to >= n
                        || !self.graph.is_adjacent(from, to)
                        || !matches!((self.distance[root][from], self.distance[root][to]),
                            (Some(a), Some(b)) if a == b + 1)
                    {
                        return Err("invalid RDL shortest-path DAG edge");
                    }
                }
            }
        }
        for family in &self.cycle_families {
            if [family.root, family.p, family.q]
                .into_iter()
                .any(|v| v >= n)
                || family.x.is_some_and(|v| v >= n)
                || family.prototype.len() != m
                || family.prototype.iter().filter(|&&b| b).count() != family.weight
            {
                return Err("invalid RDL cycle-family prototype");
            }
            let closing = match family.x {
                Some(x) => {
                    self.graph.is_adjacent(family.p, x) && self.graph.is_adjacent(family.q, x)
                }
                None => self.graph.is_adjacent(family.p, family.q),
            };
            if !closing {
                return Err("missing RDL cycle-family closing edge");
            }
            // Validate only paths used by this family. The source deliberately
            // leaves other DAG rows empty according to its vertex ordering.
            let dag = &self.directed_paths[family.root];
            let mut pending = vec![family.p, family.q];
            let mut seen = vec![false; n];
            while let Some(node) = pending.pop() {
                if node == family.root || std::mem::replace(&mut seen[node], true) {
                    continue;
                }
                if dag[node].is_empty() {
                    return Err("incomplete RDL cycle-family path");
                }
                pending.extend_from_slice(&dag[node]);
            }
        }
        if self
            .prototype_counts
            .iter()
            .try_fold(0usize, |sum, &v| sum.checked_add(v))
            != Some(self.cycle_families.len())
            || self.relations.len() != self.prototype_counts.len()
            || self
                .relations
                .iter()
                .zip(&self.prototype_counts)
                .any(|(matrix, &count)| !square(matrix, count))
        {
            return Err("invalid RDL family-relation dimensions");
        }
        let mut covered = vec![false; self.cycle_families.len()];
        for urf in &self.urfs {
            if urf.is_empty() {
                return Err("empty RDL unique ring family");
            }
            for &family in urf {
                let f = self
                    .cycle_families
                    .get(family)
                    .ok_or("RDL family index out of range")?;
                if !f.relevant || std::mem::replace(&mut covered[family], true) {
                    return Err("invalid RDL relevant family membership");
                }
            }
        }
        if covered
            .iter()
            .zip(&self.cycle_families)
            .any(|(&covered, f)| covered != f.relevant)
        {
            return Err("incomplete RDL relevant family membership");
        }
        Ok(DecomposedComponent {
            component: BiconnectedComponent {
                graph: self.graph,
                original_nodes: self.original_nodes,
                original_edges: self.original_edges,
            },
            shortest_paths: ShortestPathInfo {
                predecessor: self.predecessor,
                distance: self.distance,
                reachable: self.reachable,
                directed_paths: self
                    .directed_paths
                    .into_iter()
                    .map(|adjacency| DirectedPathGraph {
                        node_count: n,
                        adjacency,
                    })
                    .collect(),
            },
            cycle_families: CycleFamilies {
                families: self
                    .cycle_families
                    .into_iter()
                    .map(|f| CycleFamily {
                        weight: f.weight,
                        r: f.root,
                        p: f.p,
                        q: f.q,
                        x: f.x,
                        prototype: f.prototype,
                        relevant: f.relevant,
                    })
                    .collect(),
            },
            urf_info: UrfInfo {
                nof_protos: self.prototype_counts,
                relations: self.relations,
                urfs: self.urfs,
            },
        })
    }
}
