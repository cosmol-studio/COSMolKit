//! Source RDL cycle-basis and streamed relevant-cycle queries.
use super::persistence::CycleFamilyView;
use super::{DecomposedComponent, EdgeId, RingDecomposition};

/// Read-only source data needed by cycle queries. External owners can retain
/// the complete state in their own vocabulary without rebuilding or cloning it.
/// Implementations preserve component, family, path and edge index ordering.
pub trait CycleComponentAccess<'a>: Copy + 'a {
    fn node_count(self) -> usize;
    fn edge_count(self) -> usize;
    fn edge_id(self, from: usize, to: usize) -> EdgeId;
    fn original_edge(self, edge: usize) -> EdgeId;
    fn urfs(self) -> &'a [Vec<usize>];
    fn family(self, family: usize) -> CycleFamilyView<'a>;
    fn path_neighbors(self, root: usize, node: usize) -> &'a [usize];
}

impl<'a> CycleComponentAccess<'a> for &'a DecomposedComponent {
    fn node_count(self) -> usize {
        self.component.graph.node_count()
    }
    fn edge_count(self) -> usize {
        self.component.graph.edge_count()
    }
    fn edge_id(self, from: usize, to: usize) -> EdgeId {
        self.component
            .graph
            .edge_id(from, to)
            .expect("source component edge exists")
    }
    fn original_edge(self, edge: usize) -> EdgeId {
        self.component.original_edges[edge]
    }
    fn urfs(self) -> &'a [Vec<usize>] {
        &self.urf_info.urfs
    }
    fn family(self, family: usize) -> CycleFamilyView<'a> {
        let f = &self.cycle_families.families[family];
        CycleFamilyView {
            weight: f.weight,
            root: f.r,
            p: f.p,
            q: f.q,
            x: f.x,
            prototype: &f.prototype,
            relevant: f.relevant,
        }
    }
    fn path_neighbors(self, root: usize, node: usize) -> &'a [usize] {
        &self.shortest_paths.directed_paths[root].adjacency[node]
    }
}

fn bit(words: &[u64], index: usize) -> bool {
    // RDL✔️✔️: return (bitset[pos / BlockSize] >> (pos % BlockSize)) & 1;
    words[index / 64] & (1_u64 << (index % 64)) != 0
}

fn set(words: &mut [u64], index: usize) {
    // RDL✔️✔️: bitset[pos / BlockSize] |= 1 << (pos % BlockSize);
    words[index / 64] |= 1_u64 << (index % 64);
}

fn swap_columns(rows: &mut [Vec<u64>], left: usize, right: usize) {
    // RDL✔️✔️: if (!val1 != !val2) {
    // RDL✔️✔️:   RDL_bitset_flip(rows[i], col1);
    // RDL✔️✔️:   RDL_bitset_flip(rows[i], col2);
    // RDL✔️✔️: }
    for row in rows {
        if bit(row, left) != bit(row, right) {
            row[left / 64] ^= 1_u64 << (left % 64);
            row[right / 64] ^= 1_u64 << (right % 64);
        }
    }
}

impl RingDecomposition {
    /// Minimum cycle basis in source BCC/URF order, with original edge IDs.
    pub fn smallest_ring_basis(&self) -> Vec<Vec<EdgeId>> {
        smallest_ring_basis(self.components.iter())
    }

    /// Stream relevant cycles without materializing exponentially many paths.
    pub fn relevant_cycles(&self) -> impl Iterator<Item = Vec<EdgeId>> + '_ {
        relevant_cycles(self.components.iter())
    }
}

/// The same source basis algorithm for borrowed, independently owned state.
pub fn smallest_ring_basis<'a, C: CycleComponentAccess<'a>>(
    components: impl Iterator<Item = C>,
) -> Vec<Vec<EdgeId>> {
    // RDL✔️✔️: for (bcc_index = 0; bcc_index < data->bccGraphs->nof_bcc; ++bcc_index) {
    // RDL✔️✔️:   curr_graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL✔️✔️:   size=curr_graph->E - curr_graph->V + 1;
    // RDL✔️✔️:   for (internal_index = 0;
    // RDL✔️✔️:       internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index) {
    // RDL✔️✔️:     compressedSize = RDL_bitset_compressed(&(prototypesForGaussian[internal_index]),
    // RDL✔️✔️:         data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype,
    // RDL✔️✔️:         curr_graph->E);
    // RDL✔️✔️:   }
    // RDL✔️✔️:   for (internal_index = 0;
    // RDL✔️✔️:       internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index, ++urf) {
    // RDL✔️✔️:     for (k = 0; k < currentBasisCyclesSize; ++k) {
    // RDL✔️✔️:       if (RDL_bitset_test(currentCycle, k)) {
    // RDL✔️✔️:         RDL_bitset_xor_inplace(currentCycle, basisCycles[k], compressedSize);
    // RDL✔️✔️:       }
    // RDL✔️✔️:     }
    // RDL✔️✔️:     if (RDL_bitset_empty(currentCycle, empty_cycle, compressedSize)) {
    // RDL✔️✔️:       free(currentCycle);
    // RDL✔️✔️:       continue;
    // RDL✔️✔️:     }
    // RDL✔️✔️:     oldBasisCyclesSize = currentBasisCyclesSize;
    // RDL✔️✔️:     basisCycles[currentBasisCyclesSize] = currentCycle;
    // RDL✔️✔️:     ++currentBasisCyclesSize;
    // RDL✔️✔️:     if (!RDL_bitset_test(currentCycle, oldBasisCyclesSize)) {
    // RDL✔️✔️:       for (k = oldBasisCyclesSize + 1; k < curr_graph->E; ++k) {
    // RDL✔️✔️:         if (RDL_bitset_test(currentCycle, k)) {
    // RDL✔️✔️:           RDL_swap_columns(basisCycles, currentBasisCyclesSize, oldBasisCyclesSize, k);
    // RDL✔️✔️:           RDL_swap_columns(prototypesForGaussian, data->nofURFsPerBCC[bcc_index], oldBasisCyclesSize, k);
    // RDL✔️✔️:           break;
    // RDL✔️✔️:         }
    // RDL✔️✔️:       }
    // RDL✔️✔️:     }
    // RDL✔️✔️:     prototype = data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype;
    // RDL✔️✔️:     if(++added == size) {
    // RDL✔️✔️:       break;
    // RDL✔️✔️:     }
    // Packed u64 Gaussian rows have the same order and asymptotic costs
    // as compressed C byte rows; only independent output rings allocate.
    let mut result = Vec::new();
    for data in components {
        let target = data.edge_count() - data.node_count() + 1;
        let mut prototypes: Vec<Vec<u64>> = data
            .urfs()
            .iter()
            .map(|urf| {
                let mut words = vec![0; data.edge_count().div_ceil(64)];
                for (index, &present) in data.family(urf[0]).prototype.iter().enumerate() {
                    if present {
                        set(&mut words, index);
                    }
                }
                words
            })
            .collect();
        let mut basis: Vec<Vec<u64>> = Vec::with_capacity(target);
        for index in 0..prototypes.len() {
            let mut current = prototypes[index].clone();
            for (pivot, previous) in basis.iter().enumerate() {
                if bit(&current, pivot) {
                    for (word, &other) in current.iter_mut().zip(previous) {
                        *word ^= other;
                    }
                }
            }
            if current.iter().all(|&word| word == 0) {
                continue;
            }
            let pivot = basis.len();
            basis.push(current);
            if !bit(&basis[pivot], pivot) {
                let column = (pivot + 1..data.edge_count())
                    .find(|&column| bit(&basis[pivot], column))
                    .expect("a nonzero reduced source cycle has a pivot");
                swap_columns(&mut basis, pivot, column);
                swap_columns(&mut prototypes, pivot, column);
            }
            let family = data.family(data.urfs()[index][0]);
            result.push(
                family
                    .prototype
                    .iter()
                    .enumerate()
                    .filter_map(|(edge, &present)| present.then_some(data.original_edge(edge)))
                    .collect(),
            );
            if basis.len() == target {
                break;
            }
        }
    }
    result
}

/// The source streaming iterator, borrowing owner state without reconstruction.
pub fn relevant_cycles<'a, C: CycleComponentAccess<'a>>(
    components: impl Iterator<Item = C> + 'a,
) -> impl Iterator<Item = Vec<EdgeId>> + 'a {
    // RDL✔️✔️: it = RDL_initCycleIterator(
    // RDL✔️✔️:     RDL_ALL_IT,
    // RDL✔️✔️:     0, 0,
    // RDL✔️✔️:     0, 0,
    // RDL✔️✔️:     0, data->bccGraphs->nof_bcc-1,
    // RDL✔️✔️:     'b', data);
    // BCC, then URF, then RCF is exactly the source iterator nesting.
    components.flat_map(|data| {
        data.urfs().iter().flat_map(move |urf| {
            urf.iter()
                .flat_map(move |&family| FamilyCycles::new(data, family))
        })
    })
}

struct PathFrame {
    node: usize,
    next: usize,
    split_copy: Option<Vec<u64>>,
}

struct Paths<C> {
    data: C,
    root: usize,
    stack: Vec<PathFrame>,
    words: Vec<u64>,
}

impl<'a, C: CycleComponentAccess<'a>> Paths<C> {
    fn new(data: C, root: usize, target: usize) -> Self {
        // RDL✔️🔝: new_element->r = r;
        // RDL✔️🔝: new_element->p = p;
        // RDL✔️🔝: new_element->ai_index = 0;
        // RDL✔️🔝: new_element->split_copy = NULL;
        // RDL✔️🔝: RDL_stack_push(it->stack, new_element);
        // RDL✔️🔝: if (!RDL_pathIteratorNext(it)) {
        // Borrow existing source state; contiguous frames avoid one malloc
        // per source linked-stack node. Split masks retain identical lifetime.
        let mut paths = Self {
            data,
            root,
            stack: vec![PathFrame {
                node: target,
                next: 0,
                split_copy: None,
            }],
            words: vec![0; data.edge_count().div_ceil(64)],
        };
        assert!(paths.advance(), "a source cycle family has a root path");
        paths
    }

    fn advance(&mut self) -> bool {
        // RDL✔️🔝: if (RDL_stack_empty(it->stack)) {
        // RDL✔️🔝:   it->end = 1;
        // RDL✔️🔝: }
        // RDL✔️🔝: while (!RDL_stack_empty(it->stack)) {
        // RDL✔️🔝:   top = RDL_stack_top(it->stack);
        // RDL✔️🔝:   if(top->r == top->p) {
        // RDL✔️🔝:     RDL_stack_pop(it->stack);
        // RDL✔️🔝:     RDL_lpStackElement_delete(top);
        // RDL✔️🔝:     break;
        // RDL✔️🔝:   }
        // RDL✔️🔝:   if(it->spi->dPaths[top->r]->degree[top->p] > 1)/*split*/
        // RDL✔️🔝:   {
        // RDL✔️🔝:     vertex = it->spi->dPaths[top->r]->adjList[top->p][top->ai_index][0];
        // RDL✔️🔝:     if (top->ai_index == 0) {
        // RDL✔️🔝:       top->split_copy = malloc(it->compressed_number * sizeof(*top->split_copy));
        // RDL✔️🔝:       memcpy(top->split_copy, it->compressed_path,
        // RDL✔️🔝:           it->compressed_number * sizeof(*top->split_copy));
        // RDL✔️🔝:     }
        // RDL✔️🔝:     else /*not continuing the last path*/
        // RDL✔️🔝:     {
        // RDL✔️🔝:       memcpy(it->compressed_path, top->split_copy,
        // RDL✔️🔝:           it->compressed_number * sizeof(*top->split_copy));
        // RDL✔️🔝:     }
        // RDL✔️🔝:     RDL_bitset_set(it->compressed_path, RDL_edgeId(it->gra,top->p,vertex));
        // RDL✔️🔝:     ++top->ai_index;
        // RDL✔️🔝:     if (top->ai_index >= it->spi->dPaths[top->r]->degree[top->p]) {
        // RDL✔️🔝:       RDL_stack_pop(it->stack);
        // RDL✔️🔝:       RDL_lpStackElement_delete(top);
        // RDL✔️🔝:     }
        // RDL✔️🔝:     RDL_stack_push(it->stack, new_element);
        // Packed masks, source adjacency order, constant live path count.
        while let Some(top) = self.stack.last_mut() {
            if top.node == self.root {
                self.stack.pop();
                return true;
            }
            let neighbors = self.data.path_neighbors(self.root, top.node);
            let vertex = neighbors[top.next];
            let edge = self.data.edge_id(top.node, vertex);
            if neighbors.len() > 1 {
                if top.next == 0 {
                    top.split_copy = Some(self.words.clone());
                } else {
                    self.words.copy_from_slice(
                        top.split_copy
                            .as_ref()
                            .expect("source split saved its prefix"),
                    );
                }
                set(&mut self.words, edge.index());
                top.next += 1;
                if top.next == neighbors.len() {
                    self.stack.pop();
                }
                self.stack.push(PathFrame {
                    node: vertex,
                    next: 0,
                    split_copy: None,
                });
            } else {
                set(&mut self.words, edge.index());
                top.node = vertex;
                top.next = 0;
            }
        }
        false
    }
}

struct FamilyCycles<C> {
    data: C,
    family: usize,
    first: Paths<C>,
    second: Paths<C>,
    combined: Vec<u64>,
    ready: bool,
}

impl<'a, C: CycleComponentAccess<'a>> FamilyCycles<C> {
    fn new(data: C, family: usize) -> Self {
        // RDL✔️✔️: it->it2 = RDL_listPaths(cfam->r, cfam->q, it->mode, gra, spi);
        // RDL✔️✔️: it->it1 = RDL_listPaths(cfam->r, cfam->p, it->mode, gra, spi);
        let source = data.family(family);
        Self {
            data,
            family,
            first: Paths::new(data, source.root, source.p),
            second: Paths::new(data, source.root, source.q),
            combined: vec![0; data.edge_count().div_ceil(64)],
            ready: true,
        }
    }
}

impl<'a, C: CycleComponentAccess<'a>> Iterator for FamilyCycles<C> {
    type Item = Vec<EdgeId>;
    fn next(&mut self) -> Option<Self::Item> {
        // RDL✔️✔️: memcpy(compressed_cycle, it1->compressed_path,
        // RDL✔️✔️:     it1->compressed_number * sizeof(*compressed_cycle));
        // RDL✔️✔️: RDL_bitset_or_inplace(compressed_cycle, it2->compressed_path,
        // RDL✔️✔️:     it2->compressed_number);
        // RDL✔️✔️: if(cfam->x < UINT_MAX) {
        // RDL✔️✔️:   /*even cycle*/
        // RDL✔️✔️:   RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->p,cfam->x));
        // RDL✔️✔️:   RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->q,cfam->x));
        // RDL✔️✔️: }
        // RDL✔️✔️: else {
        // RDL✔️✔️:   /*odd cycle*/
        // RDL✔️✔️:   RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->p,cfam->q));
        // RDL✔️✔️: }
        // RDL✔️✔️: RDL_pathIteratorNext(it->it2);
        // RDL✔️✔️: if (it->it2->end) {
        // RDL✔️✔️:   if (!it->it1->end) {
        // RDL✔️✔️:     RDL_pathIteratorNext(it->it1);
        // RDL✔️✔️:     RDL_deletePathIterator(it->it2);
        // RDL✔️✔️:     it->it2 = RDL_listPaths(cfam->r, cfam->q, it->mode, gra, spi);
        // RDL✔️✔️:   }
        // RDL✔️✔️: }
        // Edge conversion scans component edge order, not global edge order.
        // Only two streamed paths and one packed combined mask are live.
        if !self.ready {
            return None;
        }
        let family = self.data.family(self.family);
        for ((word, &first), &second) in self
            .combined
            .iter_mut()
            .zip(&self.first.words)
            .zip(&self.second.words)
        {
            *word = first | second;
        }
        for (from, to) in if let Some(x) = family.x {
            [(family.p, x), (family.q, x)]
        } else {
            [(family.p, family.q), (family.p, family.q)]
        } {
            set(&mut self.combined, self.data.edge_id(from, to).index());
        }
        let result = (0..self.data.edge_count())
            .filter(|&edge| bit(&self.combined, edge))
            .map(|edge| self.data.original_edge(edge))
            .collect();
        if !self.second.advance() {
            if self.first.advance() {
                self.second = Paths::new(self.data, family.root, family.q);
            } else {
                self.ready = false;
            }
        }
        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Graph;
    #[test]
    fn source_basis_and_relevant_iterator_keep_theta_cycle_order() {
        let mut graph = Graph::new(5);
        for (from, to) in [(0, 2), (2, 1), (0, 3), (3, 1), (0, 4), (4, 1)] {
            graph.add_undirected_edge(from, to).unwrap();
        }
        let decomposition = RingDecomposition::calculate(graph).unwrap();
        let ids = |cycles: Vec<Vec<EdgeId>>| {
            cycles
                .into_iter()
                .map(|cycle| cycle.into_iter().map(EdgeId::index).collect::<Vec<_>>())
                .collect::<Vec<_>>()
        };
        // RDLtarjan.c:270–285 maps component edges in original index order;
        // RDL_getSSSR/RDL_cycleIteratorGetCycle scan that order, not path order.
        assert_eq!(
            ids(decomposition.smallest_ring_basis()),
            [vec![0, 1, 2, 3], vec![0, 1, 4, 5]]
        );
        assert_eq!(
            ids(decomposition.relevant_cycles().collect()),
            [vec![0, 1, 2, 3], vec![0, 1, 4, 5], vec![2, 3, 4, 5]]
        );
        assert_eq!(decomposition.relevant_cycles().count(), 3);
    }
}
