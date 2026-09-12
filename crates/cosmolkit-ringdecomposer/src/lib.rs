use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum RingDecomposerError {
    #[error("node index {node} out of range for graph with {node_count} nodes")]
    NodeOutOfRange { node: usize, node_count: usize },
    #[error("self-loop at node {node} is not allowed")]
    SelfLoop { node: usize },
    #[error("duplicate edge between {from} and {to}")]
    DuplicateEdge { from: usize, to: usize },
    #[error("graph has no nodes")]
    EmptyGraph,
    #[error("edge not found between node {from} and node {to}")]
    EdgeNotFound { from: usize, to: usize },
    #[error("relevant-cycle count exceeds the finite RDL counting range")]
    RelevantCycleCountOverflow,
    #[error("invalid shortest-path DAG while counting paths from {root} to {target}")]
    InvalidShortestPathDag { root: usize, target: usize },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct EdgeId(usize);

impl EdgeId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Edge {
    id: EdgeId,
    from: usize,
    to: usize,
}

impl Edge {
    #[must_use]
    pub const fn id(self) -> EdgeId {
        self.id
    }

    #[must_use]
    pub const fn from(self) -> usize {
        self.from
    }

    #[must_use]
    pub const fn to(self) -> usize {
        self.to
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Graph {
    node_count: usize,
    adjacency: Vec<Vec<(usize, EdgeId)>>,
    edges: Vec<Edge>,
    edge_ids_by_endpoints: BTreeMap<(usize, usize), EdgeId>,
}

impl Graph {
    #[must_use]
    pub fn new(node_count: usize) -> Self {
        // BEGIN RDL C FUNCTION RDL_initNewGraph
        // RDL✔️✔️: RDL_graph *RDL_initNewGraph(unsigned V)
        // RDL✔️✔️: {
        // RDL✔️✔️:   return RDL_initNewGraph_g(V, 1);
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_initNewGraph
        //
        // BEGIN RDL C FUNCTION RDL_initNewGraph_g
        // RDL❗✔️: RDL_graph *RDL_initNewGraph_g(unsigned V, char owns_edges)
        // RDL❗✔️: {
        // RDL❗✔️:   RDL_graph *graph = malloc(sizeof(*graph));
        // RDL❗✔️:   unsigned *degree;
        // RDL❗✔️:   unsigned i;
        // RDL❗✔️:   degree = malloc(V * sizeof(*degree));
        // RDL❗✔️:   for(i=0; i<V; ++i)
        // RDL❗✔️:   {
        // RDL❗✔️:     degree[i] = 0;
        // RDL❗✔️:   }
        // RDL❗✔️:   RDL_initGraph(graph,V,degree,owns_edges);
        // RDL❗✔️:
        // RDL❗✔️:   return graph;
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_initNewGraph_g
        Self {
            node_count,
            adjacency: vec![Vec::new(); node_count],
            edges: Vec::new(),
            edge_ids_by_endpoints: BTreeMap::new(),
        }
    }

    #[must_use]
    pub const fn node_count(&self) -> usize {
        self.node_count
    }

    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    #[must_use]
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    #[must_use]
    pub fn neighbors(&self, node: usize) -> Option<&[(usize, EdgeId)]> {
        self.adjacency.get(node).map(Vec::as_slice)
    }

    pub fn add_undirected_edge(
        &mut self,
        from: usize,
        to: usize,
    ) -> Result<EdgeId, RingDecomposerError> {
        // BEGIN RDL C FUNCTION RDL_addUEdge
        // RDL✔️✔️: unsigned RDL_addUEdge(RDL_graph *gra, RDL_node from, RDL_node to)
        // RDL✔️✔️: {
        // RDL✔️✔️:   return RDL_addUEdge_g(gra, from, to);
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_addUEdge
        //
        // BEGIN RDL C FUNCTION RDL_addUEdge_g
        // RDL✔️✔️: unsigned RDL_addUEdge_g(RDL_graph *gra, unsigned from, unsigned to)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned edge_id, i;
        // RDL✔️✔️:
        // RDL✔️✔️:   if(from >= gra->V || to >= gra->V) {
        // RDL✔️✔️:     RDL_outputFunc(RDL_ERROR, "Tried to add an edge with atoms not in range.\n");
        // RDL✔️✔️:     RDL_outputFunc(RDL_ERROR,  "edge (%u,%u) can not be added to graph with %u atoms.\n",
        // RDL✔️✔️:         from, to, gra->V);
        // RDL✔️✔️:     return RDL_INVALID_RESULT;
        // RDL✔️✔️:   }
        if from >= self.node_count {
            return Err(RingDecomposerError::NodeOutOfRange {
                node: from,
                node_count: self.node_count,
            });
        }
        if to >= self.node_count {
            return Err(RingDecomposerError::NodeOutOfRange {
                node: to,
                node_count: self.node_count,
            });
        }
        // RDL✔️✔️:
        // RDL✔️✔️:   if (from == to) {
        // RDL✔️✔️:     RDL_outputFunc(RDL_WARNING, "Adding a loop is not allowed, node %u\n", from);
        // RDL✔️✔️:     return RDL_INVALID_RESULT;
        // RDL✔️✔️:   }
        if from == to {
            return Err(RingDecomposerError::SelfLoop { node: from });
        }
        let (left, right) = canonical_edge_key(from, to);
        // RDL✔️✔️:
        // RDL✔️✔️:   for(i=0; i<gra->degree[from]; ++i) {
        // RDL✔️✔️:     if(gra->adjList[from][i][0] == to) {
        // RDL✔️✔️:       /*edge already exists*/
        // RDL✔️✔️:       return RDL_DUPLICATE_EDGE;
        // RDL✔️✔️:     }
        // RDL✔️✔️:   }
        if self.edge_ids_by_endpoints.contains_key(&(left, right)) {
            return Err(RingDecomposerError::DuplicateEdge { from, to });
        }
        // RDL✔️✔️:   RDL_addEdge(gra, from, to);
        // RDL✔️✔️:   RDL_addEdge(gra, to, from);
        // RDL✔️✔️:   --gra->E; /*was incremented twice*/
        // RDL✔️✔️:
        // RDL✔️✔️:   edge_id = RDL_addToEdgeArray(gra, from, to);
        let edge_id = EdgeId::new(self.edges.len());
        self.edges.push(Edge {
            id: edge_id,
            from: left,
            to: right,
        });
        self.edge_ids_by_endpoints.insert((left, right), edge_id);
        // RDL✔️✔️:
        // RDL✔️✔️:   gra->adjList[from][gra->degree[from]-1][1] = edge_id;
        // RDL✔️✔️:   gra->adjList[to][gra->degree[to]-1][1] = edge_id;
        self.adjacency[from].push((to, edge_id));
        self.adjacency[to].push((from, edge_id));
        // RDL✔️✔️:
        // RDL✔️✔️:   return edge_id;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_addUEdge_g
        Ok(edge_id)
    }

    pub fn edge_id(&self, from: usize, to: usize) -> Result<EdgeId, RingDecomposerError> {
        // BEGIN RDL C FUNCTION RDL_edgeId
        // RDL✔️✔️: unsigned RDL_edgeId(const RDL_graph *gra, unsigned from, unsigned to)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned j, edge;
        // RDL✔️✔️:
        // RDL✔️✔️:   if(from > to) {
        // RDL✔️✔️:     /*swap order to make from < to*/
        // RDL✔️✔️:     edge = to;
        // RDL✔️✔️:     to = from;
        // RDL✔️✔️:     from = edge;
        // RDL✔️✔️:   }
        let (left, right) = canonical_edge_key(from, to);
        // RDL✔️✔️:
        // RDL✔️✔️:   edge = RDL_INVALID_RESULT;
        // RDL✔️✔️:
        // RDL✔️✔️:   for(j=0; j<gra->degree[from]; ++j) {
        // RDL✔️✔️:     if(gra->adjList[from][j][0] == to) {
        // RDL✔️✔️:       edge = gra->adjList[from][j][1];
        // RDL✔️✔️:       break;
        // RDL✔️✔️:     }
        // RDL✔️✔️:   }
        // RDL✔️✔️:
        // RDL✔️✔️:   return edge;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_edgeId
        self.edge_ids_by_endpoints
            .get(&(left, right))
            .copied()
            .ok_or(RingDecomposerError::EdgeNotFound { from, to })
    }

    #[must_use]
    pub fn is_adjacent(&self, from: usize, to: usize) -> bool {
        // BEGIN RDL C FUNCTION RDL_isAdj
        // RDL✔️✔️: int RDL_isAdj(const RDL_graph *graph, unsigned i, unsigned j)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned idx;
        // RDL✔️✔️:   for(idx=0; idx<graph->degree[i]; ++idx)
        // RDL✔️✔️:   {
        // RDL✔️✔️:     if(graph->adjList[i][idx][0] == j)
        // RDL✔️✔️:     {
        // RDL✔️✔️:       return 1;
        // RDL✔️✔️:     }
        // RDL✔️✔️:   }
        // RDL✔️✔️:   return 0;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_isAdj
        self.adjacency
            .get(from)
            .is_some_and(|neighbors| neighbors.iter().any(|(node, _)| *node == to))
    }

    #[must_use]
    pub fn is_connected(&self) -> bool {
        // BEGIN RDL C FUNCTION RDL_checkGraphConnected
        // RDL✔️✔️: char RDL_checkGraphConnected(const RDL_graph *gra)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned i;
        // RDL✔️✔️:   char *visited;
        // RDL✔️✔️:   char result;
        // RDL✔️✔️:   result = 1;
        // RDL✔️✔️:   visited = malloc(gra->V * sizeof(*visited));
        // RDL✔️✔️:   visited[0] = 1; /*start vertex*/
        if self.node_count == 0 {
            return true;
        }
        let mut visited = vec![false; self.node_count];
        let mut stack = vec![0usize];
        visited[0] = true;
        while let Some(node) = stack.pop() {
            for &(neighbor, _) in &self.adjacency[node] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
        // RDL✔️✔️:   for(i=1; i<gra->V; ++i)
        // RDL✔️✔️:   {
        // RDL✔️✔️:     visited[i] = 0; /*unvisited*/
        // RDL✔️✔️:   }
        // RDL✔️✔️:   RDL_DFSvisit(gra, 0, visited);
        // RDL✔️✔️:   for(i=0; i<gra->V; ++i)
        // RDL✔️✔️:   {
        // RDL✔️✔️:     if(visited[i] == 0) /*one was unvisited*/
        // RDL✔️✔️:     {
        // RDL✔️✔️:       result = 0;
        // RDL✔️✔️:       break;
        // RDL✔️✔️:     }
        // RDL✔️✔️:   }
        // RDL✔️✔️:
        // RDL✔️✔️:   free(visited);
        // RDL✔️✔️:   return result;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_checkGraphConnected
        visited.into_iter().all(|is_visited| is_visited)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct DirectedPathGraph {
    node_count: usize,
    adjacency: Vec<Vec<usize>>,
}

impl DirectedPathGraph {
    #[must_use]
    fn new(node_count: usize) -> Self {
        Self {
            node_count,
            adjacency: vec![Vec::new(); node_count],
        }
    }

    fn add_directed_edge(&mut self, from: usize, to: usize) {
        // BEGIN RDL C FUNCTION RDL_addEdge
        // RDL✔️✔️: void RDL_addEdge(RDL_graph *gra, unsigned from, unsigned to)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned i;
        // RDL✔️✔️:
        // RDL✔️✔️:   /* loops */
        // RDL✔️✔️:   if (from == to) {
        // RDL✔️✔️:     return;
        // RDL✔️✔️:   }
        if from == to || from >= self.node_count || to >= self.node_count {
            return;
        }
        // RDL✔️✔️:
        // RDL✔️✔️:   for(i=0; i<gra->degree[from]; ++i) {
        // RDL✔️✔️:     if(gra->adjList[from][i][0] == to) {
        // RDL✔️✔️:       /*edge already exists*/
        // RDL✔️✔️:       return;
        // RDL✔️✔️:     }
        // RDL✔️✔️:   }
        if self.adjacency[from].contains(&to) {
            return;
        }
        // RDL✔️✔️:   ++gra->E;
        // RDL✔️✔️:   ++gra->degree[from];
        // RDL✔️✔️:   ...
        // RDL✔️✔️:   gra->adjList[from][ gra->degree[from]-1 ][0] = to;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_addEdge
        self.adjacency[from].push(to);
    }

    #[must_use]
    const fn node_count(&self) -> usize {
        self.node_count
    }

    #[must_use]
    fn neighbors(&self, node: usize) -> Option<&[usize]> {
        self.adjacency.get(node).map(Vec::as_slice)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BiconnectedComponents {
    components: Vec<BiconnectedComponent>,
    edge_to_component: Vec<Option<(usize, usize)>>,
    node_to_components: Vec<Vec<(usize, usize)>>,
}

impl BiconnectedComponents {
    #[must_use]
    fn calculate(graph: &Graph) -> Self {
        tarjan_biconnected_components(graph)
    }

    #[must_use]
    fn components(&self) -> &[BiconnectedComponent] {
        &self.components
    }

    #[must_use]
    fn component_count(&self) -> usize {
        self.components.len()
    }

    #[must_use]
    fn edge_component(&self, edge: EdgeId) -> Option<(usize, usize)> {
        self.edge_to_component.get(edge.index()).copied().flatten()
    }

    #[must_use]
    fn node_components(&self, node: usize) -> Option<&[(usize, usize)]> {
        self.node_to_components.get(node).map(Vec::as_slice)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct BiconnectedComponent {
    graph: Graph,
    original_nodes: Vec<usize>,
    original_edges: Vec<EdgeId>,
}

impl BiconnectedComponent {
    #[must_use]
    const fn graph(&self) -> &Graph {
        &self.graph
    }

    #[must_use]
    fn original_nodes(&self) -> &[usize] {
        &self.original_nodes
    }

    #[must_use]
    fn original_edges(&self) -> &[EdgeId] {
        &self.original_edges
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ShortestPathInfo {
    predecessor: Vec<Vec<Option<usize>>>,
    distance: Vec<Vec<Option<usize>>>,
    reachable: Vec<Vec<bool>>,
    directed_paths: Vec<DirectedPathGraph>,
}

impl ShortestPathInfo {
    #[must_use]
    fn calculate(graph: &Graph) -> Self {
        // BEGIN RDL C FUNCTION RDL_AllPairsShortestPaths
        // RDL✔️✔️: RDL_sPathInfo *RDL_AllPairsShortestPaths(RDL_graph *gra)
        // RDL✔️✔️: {
        // RDL✔️✔️:   RDL_sPathInfo *info = malloc(sizeof(*info));
        // RDL✔️✔️:
        // RDL✔️✔️:   RDL_initializeSPathInfo(info, gra);
        // RDL✔️✔️:   RDL_findpaths(info, gra);
        // RDL✔️✔️:
        // RDL✔️✔️:   return info;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_AllPairsShortestPaths
        let mut info = Self::initialized_for_graph(graph);
        info.find_paths(graph);
        info
    }

    fn initialized_for_graph(graph: &Graph) -> Self {
        // BEGIN RDL C FUNCTION RDL_initializeSPathInfo
        // RDL✔️✔️: static void RDL_initializeSPathInfo(RDL_sPathInfo *info, RDL_graph *gra)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned i;
        // RDL✔️✔️:   info->pred = RDL_alloc2DUIntArray(gra->V, gra->V);
        // RDL✔️✔️:   info->dist = RDL_alloc2DUIntArray(gra->V, gra->V);
        // RDL✔️✔️:   info->reachable = RDL_alloc2DCharArray(gra->V, gra->V);
        // RDL✔️✔️:   info->dPaths = malloc(gra->V * sizeof(*info->dPaths));
        // RDL✔️✔️:   for(i=0; i<gra->V; ++i)
        // RDL✔️✔️:   {
        // RDL✔️✔️:     info->dPaths[i] = RDL_initNewGraph_g(gra->V, 0);
        // RDL✔️✔️:   }
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_initializeSPathInfo
        let node_count = graph.node_count();
        Self {
            predecessor: vec![vec![None; node_count]; node_count],
            distance: vec![vec![None; node_count]; node_count],
            reachable: vec![vec![false; node_count]; node_count],
            directed_paths: vec![DirectedPathGraph::new(node_count); node_count],
        }
    }

    fn find_paths(&mut self, graph: &Graph) {
        // BEGIN RDL C FUNCTION RDL_findpaths
        // RDL✔️✔️: static void RDL_findpaths(RDL_sPathInfo *spi, RDL_graph *gra)
        // RDL✔️✔️: {
        // RDL✔️✔️:   const unsigned INFINITY = UINT_MAX;
        // RDL✔️✔️:   unsigned i,j,w,adj,run;
        // RDL✔️✔️:   unsigned q_head, q_nextfree, q_size;
        // RDL✔️✔️:   char *color = malloc(gra->V * sizeof(*color));
        // RDL✔️✔️:   unsigned *queue = malloc(gra->V * sizeof(*queue));
        // RDL✔️✔️:
        // RDL✔️✔️:   for(run=1; run<3; ++run) /*2 run throughs to get Vr correct*/
        // RDL✔️✔️:   {
        let node_count = graph.node_count();
        for run in 1..3 {
            // RDL✔️✔️:     for(i=0; i<gra->V; ++i)
            for source in 0..node_count {
                // RDL✔️✔️:       q_head=0;
                // RDL✔️✔️:       q_nextfree=0;
                // RDL✔️✔️:       q_size=0;
                let mut queue = Vec::with_capacity(node_count);
                let mut queue_head = 0usize;
                let mut color = vec!['w'; node_count];
                // RDL✔️✔️:       for(j=0; j<gra->V; ++j)
                // RDL✔️✔️:       {
                // RDL✔️✔️:         if(run==1)
                // RDL✔️✔️:         {
                if run == 1 {
                    for target in 0..node_count {
                        // RDL✔️✔️:           spi->dist[i][j] = INFINITY;
                        // RDL✔️✔️:           spi->pred[i][j] = INFINITY;
                        // RDL✔️✔️:           spi->reachable[i][j] = 0;
                        self.distance[source][target] = None;
                        self.predecessor[source][target] = None;
                        self.reachable[source][target] = false;
                    }
                }
                // RDL✔️✔️:         }
                // RDL✔️✔️:         color[j] = 'w'; /*white*/
                // RDL✔️✔️:       }
                // RDL✔️✔️:       spi->dist[i][i] = 0;
                // RDL✔️✔️:       spi->pred[i][i] = i;
                self.distance[source][source] = Some(0);
                self.predecessor[source][source] = Some(source);
                // RDL✔️✔️:       color[i] = 'b'; /*black*/
                color[source] = 'b';
                // RDL✔️✔️:       queue[q_nextfree]=i; /*enqueue*/
                // RDL✔️✔️:       ++q_nextfree; /*enqueue*/
                // RDL✔️✔️:       ++q_size; /*enqueue*/
                queue.push(source);
                // RDL✔️✔️:
                // RDL✔️✔️:       while(q_size > 0)
                while queue_head < queue.len() {
                    // RDL✔️✔️:         w = queue[q_head]; /*deqeue*/
                    // RDL✔️✔️:         ++q_head; /*dequeue*/
                    // RDL✔️✔️:         --q_size; /*dequeue*/
                    let current = queue[queue_head];
                    queue_head += 1;
                    // RDL✔️✔️:         for(adj=0; adj<gra->degree[w]; ++adj) /*for each node adj to w*/
                    for &(neighbor, _) in &graph.adjacency[current] {
                        // RDL✔️✔️:           j=gra->adjList[w][adj][0];
                        // RDL✔️✔️:           /*if j precedes i in order (only first run)*/
                        // RDL✔️✔️:           if(run==2 || ((gra->degree[j] == gra->degree[i]) && j<i) ||
                        // RDL✔️✔️:           gra->degree[j] < gra->degree[i])
                        if run != 2
                            && !((graph.adjacency[neighbor].len() == graph.adjacency[source].len()
                                && neighbor < source)
                                || graph.adjacency[neighbor].len() < graph.adjacency[source].len())
                        {
                            continue;
                        }
                        // RDL✔️✔️:           if(color[j] == 'w') /*unvisited*/
                        if color[neighbor] != 'w' {
                            continue;
                        }
                        // RDL✔️✔️:           {
                        let candidate_distance =
                            self.distance[source][current].expect("queued vertex is reachable") + 1;
                        // RDL✔️✔️:             if(run==2)
                        if run == 2 {
                            // RDL✔️✔️:             {/*if in the 2nd run a dist gets shorter*/
                            // RDL✔️✔️:               if(spi->dist[i][w]+1 < spi->dist[i][j])
                            // RDL✔️✔️:               {/*not element of Vr*/
                            // RDL✔️✔️:                 spi->reachable[i][j] = 0;
                            // RDL✔️✔️:               }
                            if self.distance[source][neighbor]
                                .is_none_or(|distance| candidate_distance < distance)
                            {
                                self.reachable[source][neighbor] = false;
                            }
                        }
                        // RDL✔️✔️:             }
                        // RDL✔️✔️:             if(run==1 || (spi->dist[i][w]+1 < spi->dist[i][j]))
                        // RDL✔️✔️:             {/*if 2nd run and dist stays the same, pred shouldn't
                        // RDL✔️✔️:             change to keep the shortest path along ordering*/
                        if run == 1
                            || self.distance[source][neighbor]
                                .is_none_or(|distance| candidate_distance < distance)
                        {
                            // RDL✔️✔️:               /*predecessor of j on a shortest path from i to j*/
                            // RDL✔️✔️:               spi->pred[i][j] = w;
                            self.predecessor[source][neighbor] = Some(current);
                        }
                        // RDL✔️✔️:             }
                        // RDL✔️✔️:             if(spi->dist[i][j] > spi->dist[i][w]+1)
                        // RDL✔️✔️:             {
                        // RDL✔️✔️:               spi->dist[i][j] = spi->dist[i][w] + 1;
                        // RDL✔️✔️:             }
                        if self.distance[source][neighbor]
                            .is_none_or(|distance| candidate_distance < distance)
                        {
                            self.distance[source][neighbor] = Some(candidate_distance);
                        }
                        // RDL✔️✔️:             color[j] = 'b';
                        color[neighbor] = 'b';
                        // RDL✔️✔️:             queue[q_nextfree] = j; /*enqueue*/
                        // RDL✔️✔️:             ++q_nextfree; /*enqueue*/
                        // RDL✔️✔️:             ++q_size; /*enqueue*/
                        queue.push(neighbor);
                        // RDL✔️✔️:             if(run==1)
                        // RDL✔️✔️:             {/*reachable should not change to 1 in the 2nd run*/
                        // RDL✔️✔️:               spi->reachable[i][j] = 1;
                        // RDL✔️✔️:             }
                        if run == 1 {
                            self.reachable[source][neighbor] = true;
                        }
                    }
                    // RDL✔️✔️:           }
                    // RDL✔️✔️:         }
                    // RDL✔️✔️:       }
                }
            }
            // RDL✔️✔️:     }
            // RDL✔️✔️:   }
        }
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_findpaths
    }

    #[must_use]
    fn predecessor(&self, from: usize, to: usize) -> Option<usize> {
        self.predecessor
            .get(from)
            .and_then(|row| row.get(to))
            .copied()
            .flatten()
    }

    #[must_use]
    fn distance(&self, from: usize, to: usize) -> Option<usize> {
        self.distance
            .get(from)
            .and_then(|row| row.get(to))
            .copied()
            .flatten()
    }

    #[must_use]
    fn reachable_preceding(&self, from: usize, to: usize) -> bool {
        self.reachable
            .get(from)
            .and_then(|row| row.get(to))
            .copied()
            .unwrap_or(false)
    }

    #[must_use]
    fn directed_path_graphs(&self) -> &[DirectedPathGraph] {
        &self.directed_paths
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CycleFamily {
    weight: usize,
    r: usize,
    p: usize,
    q: usize,
    x: Option<usize>,
    prototype: Vec<bool>,
    relevant: bool,
}

impl CycleFamily {
    #[must_use]
    const fn weight(&self) -> usize {
        self.weight
    }

    #[must_use]
    const fn root(&self) -> usize {
        self.r
    }

    #[must_use]
    const fn p(&self) -> usize {
        self.p
    }

    #[must_use]
    const fn q(&self) -> usize {
        self.q
    }

    #[must_use]
    const fn x(&self) -> Option<usize> {
        self.x
    }

    #[must_use]
    fn prototype(&self) -> &[bool] {
        &self.prototype
    }

    #[must_use]
    const fn is_relevant(&self) -> bool {
        self.relevant
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct CycleFamilies {
    families: Vec<CycleFamily>,
}

impl CycleFamilies {
    #[must_use]
    fn calculate(graph: &mut Graph, shortest_paths: &mut ShortestPathInfo) -> Self {
        find_cycle_families(graph, shortest_paths)
    }

    #[must_use]
    fn families(&self) -> &[CycleFamily] {
        &self.families
    }

    #[must_use]
    fn len(&self) -> usize {
        self.families.len()
    }

    #[must_use]
    fn is_empty(&self) -> bool {
        self.families.is_empty()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct UrfInfo {
    nof_protos: Vec<usize>,
    relations: Vec<Vec<Vec<bool>>>,
    urfs: Vec<Vec<usize>>,
}

impl UrfInfo {
    fn check_urf_relation(
        cycle_families: &mut CycleFamilies,
        graph: &Graph,
        shortest_paths: &ShortestPathInfo,
    ) -> Self {
        // BEGIN RDL C FUNCTION RDL_checkURFRelation
        // RDL❗✔️: RDL_URFinfo *RDL_checkURFRelation(RDL_cfURF *RCFs, RDL_graph *graph, RDL_sPathInfo* spi)
        // RDL❗✔️: {
        // RDL❗✔️:   RDL_URFinfo *uInfo = RDL_initUrfInfo(RCFs, graph);
        let mut info = Self::init_urf_info(cycle_families);
        // RDL❗✔️:   RDL_findRelations(RCFs, graph, uInfo, spi);
        info.find_relations(cycle_families, graph, shortest_paths);
        // RDL❗✔️:
        // RDL❗✔️:   uInfo->nofURFs = RDL_countURFs(uInfo);
        // RDL❗✔️:   RDL_fillURFs(uInfo, RCFs);
        info.fill_urfs();
        // RDL❗✔️:   return uInfo;
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_checkURFRelation
        info
    }

    fn init_urf_info(cycle_families: &CycleFamilies) -> Self {
        // BEGIN RDL C FUNCTION RDL_initUrfInfo
        // RDL❗✔️: RDL_URFinfo *RDL_initUrfInfo(RDL_cfURF *RCFs, RDL_graph *graph)
        // RDL❗✔️: {
        // RDL❗✔️:   RDL_URFinfo *urfInfo;
        // RDL❗✔️:   char ***URFrel;
        // RDL❗✔️:   unsigned numOfWeights = 1;
        // RDL❗✔️:   unsigned *numOfProtos;
        // RDL❗✔️:   unsigned i, j, k, currWeight, currIdx;
        if cycle_families.is_empty() {
            return Self {
                nof_protos: Vec::new(),
                relations: Vec::new(),
                urfs: Vec::new(),
            };
        }
        // RDL❗✔️:
        // RDL❗✔️:   /* count number of different weights occuring*/
        let mut nof_protos = Vec::<usize>::new();
        let mut current_weight = cycle_families.families[0].weight();
        nof_protos.push(0);
        for family in cycle_families.families() {
            if family.weight() != current_weight {
                current_weight = family.weight();
                nof_protos.push(0);
            }
            *nof_protos
                .last_mut()
                .expect("weight bucket exists after initialization") += 1;
        }
        // RDL❗✔️:
        // RDL❗✔️:   /*allocate everything*/
        // RDL❗✔️:   /*initialize with zeros*/
        let relations = nof_protos
            .iter()
            .map(|&count| vec![vec![false; count]; count])
            .collect::<Vec<_>>();
        // RDL❗✔️:
        // RDL❗✔️:   urfInfo->nofWeights = numOfWeights;
        // RDL❗✔️:   urfInfo->nofProtos = numOfProtos;
        // RDL❗✔️:   urfInfo->URFrel = URFrel;
        // RDL❗✔️:
        // RDL❗✔️:   return urfInfo;
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_initUrfInfo
        Self {
            nof_protos,
            relations,
            urfs: Vec::new(),
        }
    }

    fn idx_weight(&self, weight: usize, index: usize) -> usize {
        // BEGIN RDL C FUNCTION RDL_idxWeight
        // RDL✔️✔️: unsigned RDL_idxWeight(RDL_URFinfo *uInfo, unsigned weight, unsigned j)
        // RDL✔️✔️: {
        // RDL✔️✔️:   unsigned i, sum = 0;
        // RDL✔️✔️:   for(i=0; i<weight; ++i)
        // RDL✔️✔️:   {
        // RDL✔️✔️:     sum += uInfo->nofProtos[i];
        // RDL✔️✔️:   }
        // RDL✔️✔️:   return sum + j;
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_idxWeight
        self.nof_protos.iter().take(weight).sum::<usize>() + index
    }

    fn find_relations(
        &mut self,
        cycle_families: &mut CycleFamilies,
        graph: &Graph,
        shortest_paths: &ShortestPathInfo,
    ) {
        // BEGIN RDL C FUNCTION RDL_findRelations
        // RDL✔️✔️: void RDL_findRelations(RDL_cfURF *RCFs, RDL_graph *graph,
        // RDL✔️✔️:     RDL_URFinfo *uInfo, RDL_sPathInfo* spi)
        // RDL✔️✔️: {
        // RDL✔️✔️:   RDL_checkDependencies(RCFs, graph, uInfo);
        self.check_dependencies(cycle_families, graph);
        // RDL✔️✔️:   RDL_checkEdges(RCFs, graph, uInfo, spi);
        self.check_edges(cycle_families, graph, shortest_paths);
        // RDL✔️✔️:   RDL_findTransitiveClosure(uInfo);
        self.find_transitive_closure();
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_findRelations
    }

    fn check_dependencies(&mut self, cycle_families: &mut CycleFamilies, graph: &Graph) {
        // BEGIN RDL C FUNCTION RDL_checkDependencies
        // RDL❗✔️: void RDL_checkDependencies(RDL_cfURF *RCFs, RDL_graph *graph, RDL_URFinfo *uInfo)
        // RDL❗✔️: {
        // RDL❗✔️:   if(RCFs->nofFams < 3) {
        if cycle_families.len() < 3 {
            // RDL❗✔️:     /*if only 0, 1 or 2 families exist, they are all relevant and independent*/
            for weight in 0..self.nof_protos.len() {
                for index in 0..self.nof_protos[weight] {
                    self.relations[weight][index][index] = true;
                }
            }
            for family in &mut cycle_families.families {
                family.relevant = true;
            }
            // RDL❗✔️:     return;
            return;
        }
        // RDL❗✔️:   }
        let mut basis_cycles = Vec::<Vec<bool>>::with_capacity(
            graph
                .edge_count()
                .saturating_sub(graph.node_count())
                .saturating_add(1),
        );
        let mut relevant_cycles = Vec::<Vec<bool>>::new();
        let mut relevant_cycles_map = Vec::<usize>::new();
        let mut prototypes_for_gaussian = cycle_families
            .families()
            .iter()
            .map(|family| family.prototype().to_vec())
            .collect::<Vec<_>>();
        // RDL❗✔️:
        // RDL❗✔️:   /* iterate over all weights */
        for weight in 0..self.nof_protos.len() {
            // RDL❗✔️:     currentBasisCyclesSmallerEnd = currentBasisCyclesSize;
            // RDL❗✔️:     currentRelevantCyclesSmallerEnd = currentRelevantCyclesSize;
            let current_basis_cycles_smaller_end = basis_cycles.len();
            let current_relevant_cycles_smaller_end = relevant_cycles.len();
            for index_in_weight in 0..self.nof_protos[weight] {
                // RDL❗✔️:       index = RDL_idxWeights(uInfo, i, j);
                let family_index = self.idx_weight(weight, index_in_weight);
                // RDL❗✔️:       /* copy current cycle */
                let mut cycle_with_smaller_added = prototypes_for_gaussian[family_index].clone();
                // RDL❗✔️:       /* add smaller cycles of the basis (set B<) */
                for pivot in 0..current_basis_cycles_smaller_end {
                    if cycle_with_smaller_added[pivot] {
                        xor_bool_slice(&mut cycle_with_smaller_added, &basis_cycles[pivot]);
                    }
                }
                // RDL❗✔️:
                // RDL❗✔️:       if (RDL_bitset_empty(compressedCycleWithSmallerAdded, empty_cycle, compressedSize)) {
                if !cycle_with_smaller_added.iter().any(|&bit| bit) {
                    // RDL❗✔️:         continue;
                    continue;
                }
                // RDL❗✔️:
                // RDL❗✔️:       /* this cycle does not depend on strictly smaller cycles => relevant */
                relevant_cycles.push(cycle_with_smaller_added.clone());
                relevant_cycles_map.push(index_in_weight);
                cycle_families.families[family_index].relevant = true;
                self.relations[weight][index_in_weight][index_in_weight] = true;
                // RDL❗✔️:
                // RDL❗✔️:       /* make a copy for adding the equal sized cycles */
                let mut cycle_with_equal_added = cycle_with_smaller_added.clone();
                // RDL❗✔️:       /* add equal sized cycles of the basis (set B=) */
                for pivot in current_basis_cycles_smaller_end..basis_cycles.len() {
                    if cycle_with_equal_added[pivot] {
                        xor_bool_slice(&mut cycle_with_equal_added, &basis_cycles[pivot]);
                    }
                }
                // RDL❗✔️:
                // RDL❗✔️:       if (RDL_bitset_empty(compressedCycleWithEqualAdded, empty_cycle, compressedSize)) {
                if !cycle_with_equal_added.iter().any(|&bit| bit) {
                    // RDL❗✔️:         for (k = currentRelevantCyclesSmallerEnd; k < currentRelevantCyclesSize-1; ++k) {
                    for relevant_index in
                        current_relevant_cycles_smaller_end..relevant_cycles.len() - 1
                    {
                        let mut comparison = cycle_with_smaller_added.clone();
                        xor_bool_slice(&mut comparison, &relevant_cycles[relevant_index]);
                        if !comparison.iter().any(|&bit| bit) {
                            let related_index = relevant_cycles_map[relevant_index];
                            self.relations[weight][index_in_weight][related_index] = true;
                            self.relations[weight][related_index][index_in_weight] = true;
                        }
                    }
                } else {
                    // RDL❗✔️:         oldBasisCyclesSize = currentBasisCyclesSize;
                    let old_basis_cycles_size = basis_cycles.len();
                    // RDL❗✔️:         compressedBasisCycles[currentBasisCyclesSize] = compressedCycleWithEqualAdded;
                    basis_cycles.push(cycle_with_equal_added);
                    // RDL❗✔️:
                    // RDL❗✔️:         if (!RDL_bitset_test(compressedCycleWithEqualAdded, oldBasisCyclesSize)) {
                    if !basis_cycles[old_basis_cycles_size][old_basis_cycles_size] {
                        for column in old_basis_cycles_size + 1..graph.edge_count() {
                            if basis_cycles[old_basis_cycles_size][column] {
                                swap_bool_columns(&mut basis_cycles, old_basis_cycles_size, column);
                                swap_bool_columns(
                                    &mut relevant_cycles,
                                    old_basis_cycles_size,
                                    column,
                                );
                                swap_bool_columns(
                                    &mut prototypes_for_gaussian,
                                    old_basis_cycles_size,
                                    column,
                                );
                                break;
                            }
                        }
                    }
                }
            }
        }
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_checkDependencies
    }

    fn check_edges(
        &mut self,
        cycle_families: &CycleFamilies,
        graph: &Graph,
        shortest_paths: &ShortestPathInfo,
    ) {
        // BEGIN RDL C FUNCTION RDL_checkEdges
        // RDL❗✔️: void RDL_checkEdges(RDL_cfURF *RCFs, RDL_graph *graph, RDL_URFinfo *uInfo, RDL_sPathInfo* spi)
        // RDL❗✔️: {
        // RDL❗✔️:   for(i=0; i<uInfo->nofWeights; ++i) /*go through all matrices*/
        for weight in 0..self.nof_protos.len() {
            let mut edge_lists = vec![None; self.nof_protos[weight]];
            for left in 0..self.nof_protos[weight] {
                for right in left + 1..self.nof_protos[weight] {
                    // RDL❗✔️:         if(uInfo->URFrel[i][j][k] == 1)
                    if !self.relations[weight][left][right] {
                        continue;
                    }
                    if edge_lists[left].is_none() {
                        edge_lists[left] = Some(make_edge_list(
                            &cycle_families.families[self.idx_weight(weight, left)],
                            graph,
                            shortest_paths,
                        ));
                    }
                    if edge_lists[right].is_none() {
                        edge_lists[right] = Some(make_edge_list(
                            &cycle_families.families[self.idx_weight(weight, right)],
                            graph,
                            shortest_paths,
                        ));
                    }
                    let shared = sorted_edge_lists_share_edge(
                        edge_lists[left].as_deref().expect("left edge list exists"),
                        edge_lists[right]
                            .as_deref()
                            .expect("right edge list exists"),
                    );
                    // RDL❗✔️:           /* if shared edge, relation is confirmed */
                    self.relations[weight][left][right] = shared;
                    self.relations[weight][right][left] = shared;
                }
            }
        }
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_checkEdges
    }

    fn find_transitive_closure(&mut self) {
        // BEGIN RDL C FUNCTION RDL_findTransitiveClosure
        // RDL✔️✔️: void RDL_findTransitiveClosure(RDL_URFinfo *uInfo)
        // RDL✔️✔️: {
        // RDL✔️✔️:   for(i = 0; i < uInfo->nofWeights; ++i) {
        for weight in 0..self.nof_protos.len() {
            let mut visited = vec![false; self.nof_protos[weight]];
            let mut connected_components = Vec::<Vec<usize>>::new();
            // RDL✔️✔️:     for(j = 0; j < uInfo->nofProtos[i]; ++j) {
            for start in 0..self.nof_protos[weight] {
                if visited[start] {
                    continue;
                }
                let mut component = Vec::new();
                let mut stack = vec![start];
                visited[start] = true;
                while let Some(current) = stack.pop() {
                    component.push(current);
                    for (next, &related) in self.relations[weight][current].iter().enumerate() {
                        if related && !visited[next] {
                            visited[next] = true;
                            stack.push(next);
                        }
                    }
                }
                connected_components.push(component);
            }
            // RDL✔️✔️:
            // RDL✔️✔️:     /* every pair inside a CC is URF related */
            for component in connected_components {
                for left in 0..component.len() {
                    for right in left + 1..component.len() {
                        let u = component[left];
                        let v = component[right];
                        self.relations[weight][u][v] = true;
                        self.relations[weight][v][u] = true;
                    }
                }
            }
        }
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_findTransitiveClosure
    }

    fn fill_urfs(&mut self) {
        // BEGIN RDL C FUNCTION RDL_countURFs
        // RDL✔️✔️: unsigned RDL_countURFs(RDL_URFinfo *uInfo)
        // RDL✔️✔️: {
        // RDL✔️✔️:   /*count number of 1s indicating URF-relation*/
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_countURFs
        //
        // BEGIN RDL C FUNCTION RDL_fillURFs
        // RDL✔️✔️: void RDL_fillURFs(RDL_URFinfo *uInfo, RDL_cfURF *CFs)
        // RDL✔️✔️: {
        // RDL✔️✔️:   for(i=0; i<uInfo->nofWeights; ++i)
        let mut urfs = Vec::new();
        for weight in 0..self.nof_protos.len() {
            let mut already_in_urf = vec![false; self.nof_protos[weight]];
            for start in 0..self.nof_protos[weight] {
                // RDL✔️✔️:       if((uInfo->URFrel[i][j][j] == 1) && (alreadyInURF[i][j] == 0))
                if !self.relations[weight][start][start] || already_in_urf[start] {
                    continue;
                }
                let mut urf = Vec::new();
                for current in start..self.nof_protos[weight] {
                    if self.relations[weight][start][current] {
                        urf.push(self.idx_weight(weight, current));
                        already_in_urf[current] = true;
                    }
                }
                urfs.push(urf);
            }
        }
        // RDL✔️✔️: }
        // END RDL C FUNCTION RDL_fillURFs
        self.urfs = urfs;
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RingDecomposition {
    graph: Graph,
    urfs: Vec<UniqueRingFamily>,
    relevant_cycle_count: f64,
}

impl RingDecomposition {
    pub fn calculate(graph: Graph) -> Result<Self, RingDecomposerError> {
        // BEGIN RDL C FUNCTION RDL_calculate
        // RDL❗✔️: RDL_data *RDL_calculate(RDL_graph *gra)
        // RDL❗✔️: {
        // RDL❗✔️:   RDL_data *data;
        // RDL❗✔️:   unsigned i, j, k, urf_index, rcf_index;
        // RDL❗✔️:   unsigned nof_relevant_fams, nof_relevant_fams_sum;
        // RDL❗✔️:
        // RDL❗✔️:   if (!gra) {
        // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "The graph is NULL.\n");
        // RDL❗✔️:     return NULL;
        // RDL❗✔️:   }
        // RDL❗✔️:
        // RDL❗✔️:   /* we can't calculate anything if there isn't at least one node */
        // RDL❗✔️:   if(!gra->V) {
        // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "The graph has no nodes.\n");
        // RDL❗✔️:     return NULL;
        // RDL❗✔️:   }
        // RDL❗✔️:   data = malloc(sizeof(*data));
        // RDL❗✔️:
        // RDL❗✔️:   /* FIRST STEP: TARJAN */
        // RDL❗✔️:   data->bccGraphs = RDL_tarjanBCC(gra);
        // RDL❗✔️:   data->nofURFs = 0;
        // RDL❗✔️:   data->nofRCFs = 0;
        // RDL❗✔️:
        // RDL❗✔️:   /* allocate result structures */
        // RDL❗✔️:   data->spiPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->spiPerBCC));
        // RDL❗✔️:   data->CFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->CFsPerBCC));
        // RDL❗✔️:   data->urfInfoPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->urfInfoPerBCC));
        // RDL❗✔️:   data->nofURFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->nofURFsPerBCC));
        // RDL❗✔️:   data->nofRCFsPerBCC = malloc(data->bccGraphs->nof_bcc * sizeof(*data->nofRCFsPerBCC));
        // RDL❗✔️:
        // RDL❗✔️:   /* the followeing steps are performed for each BCC */
        // RDL❗✔️:   for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
        // RDL❗✔️:     /* solve APSP problem */
        // RDL❗✔️:     data->spiPerBCC[i] = RDL_AllPairsShortestPaths(data->bccGraphs->bcc_graphs[i]);
        // RDL❗✔️:     /* calculate RCFs with Vismara's algortihm */
        // RDL❗✔️:     data->CFsPerBCC[i] = RDL_findCycleFams(data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
        // RDL❗✔️:     /* this can fail if the BCC is to large */
        // RDL❗✔️:     if (!data->CFsPerBCC[i]) {
        // RDL❗✔️:       /* delete families and APSP information UPTO current index */
        // RDL❗✔️:       for (j = 0; j < i; ++j) {
        // RDL❗✔️:         RDL_deleteAPSP(data->spiPerBCC[j], data->bccGraphs->bcc_graphs[j]->V);
        // RDL❗✔️:         RDL_deleteCycleFams(data->CFsPerBCC[j]);
        // RDL❗✔️:         if (data->nofURFsPerBCC[j] > 0) {
        // RDL❗✔️:           RDL_deleteURFInfo(data->urfInfoPerBCC[j]);
        // RDL❗✔️:         }
        // RDL❗✔️:       }
        // RDL❗✔️:       /* delete APSP for current index (there is no RCF info) */
        // RDL❗✔️:       RDL_deleteAPSP(data->spiPerBCC[i], data->bccGraphs->bcc_graphs[i]->V);
        // RDL❗✔️:
        // RDL❗✔️:       /* free alloced per BCC structures */
        // RDL❗✔️:       free(data->spiPerBCC);
        // RDL❗✔️:       free(data->CFsPerBCC);
        // RDL❗✔️:       free(data->nofURFsPerBCC);
        // RDL❗✔️:       free(data->nofRCFsPerBCC);
        // RDL❗✔️:       free(data->urfInfoPerBCC);
        // RDL❗✔️:       /* and of course the BCC graph */
        // RDL❗✔️:       RDL_deleteBCCGraph(data->bccGraphs);
        // RDL❗✔️:       /* free resulting struct */
        // RDL❗✔️:       free(data);
        // RDL❗✔️:
        // RDL❗✔️:       return NULL;
        // RDL❗✔️:     }
        // RDL❗✔️:     if(data->CFsPerBCC[i]->nofFams > 0) {
        // RDL❗✔️:       /* if there is at least one RCF, check URF relation */
        // RDL❗✔️:       data->urfInfoPerBCC[i] = RDL_checkURFRelation(data->CFsPerBCC[i],
        // RDL❗✔️:           data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
        // RDL❗✔️:       data->nofURFsPerBCC[i] = data->urfInfoPerBCC[i]->nofURFs;
        // RDL❗✔️:
        // RDL❗✔️:       nof_relevant_fams = 0;
        // RDL❗✔️:       /* count RCFs */
        // RDL❗✔️:       for (j = 0; j < data->CFsPerBCC[i]->nofFams; ++j) {
        // RDL❗✔️:         if (data->CFsPerBCC[i]->fams[j]->mark) {
        // RDL❗✔️:           ++nof_relevant_fams;
        // RDL❗✔️:         }
        // RDL❗✔️:       }
        // RDL❗✔️:
        // RDL❗✔️:       nof_relevant_fams_sum = 0;
        // RDL❗✔️:       for (j = 0; j < data->nofURFsPerBCC[i]; ++j) {
        // RDL❗✔️:         nof_relevant_fams_sum += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
        // RDL❗✔️:       }
        // RDL❗✔️:
        // RDL❗✔️:       if (nof_relevant_fams != nof_relevant_fams_sum) {
        // RDL❗✔️:         RDL_outputFunc(RDL_ERROR, "different number of relevant families!\n");
        // RDL❗✔️:         /* internal check, should never happen */
        // RDL❗✔️:         assert(0);
        // RDL❗✔️:       }
        // RDL❗✔️:       data->nofRCFsPerBCC[i] = nof_relevant_fams;
        // RDL❗✔️:     }
        // RDL❗✔️:     else {
        // RDL❗✔️:       data->nofURFsPerBCC[i] = 0;
        // RDL❗✔️:       data->nofRCFsPerBCC[i] = 0;
        // RDL❗✔️:     }
        // RDL❗✔️:     data->nofURFs += data->nofURFsPerBCC[i];
        // RDL❗✔️:     data->nofRCFs += data->nofRCFsPerBCC[i];
        // RDL❗✔️:   }
        // RDL❗✔️:
        // RDL❗✔️:   /* create a mapping from URFs to BCCs */
        // RDL❗✔️:   data->urf_to_bcc = malloc(data->nofURFs * sizeof(*data->urf_to_bcc));
        // RDL❗✔️:   /* create a mapping from RCFs to URFs */
        // RDL❗✔️:   data->rcf_to_urf = malloc(data->nofRCFs * sizeof(*data->rcf_to_urf));
        // RDL❗✔️:   urf_index = 0;
        // RDL❗✔️:   rcf_index = 0;
        // RDL❗✔️:   for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
        // RDL❗✔️:     for (j = 0; j < data->nofURFsPerBCC[i]; ++j, ++urf_index) {
        // RDL❗✔️:       data->urf_to_bcc[urf_index][0] = i;
        // RDL❗✔️:       data->urf_to_bcc[urf_index][1] = j;
        // RDL❗✔️:       for (k = 0; k < data->urfInfoPerBCC[i]->nofCFsPerURF[j]; ++k, ++rcf_index) {
        // RDL❗✔️:         data->rcf_to_urf[rcf_index][0] = urf_index;
        // RDL❗✔️:         data->rcf_to_urf[rcf_index][1] = k;
        // RDL❗✔️:       }
        // RDL❗✔️:     }
        // RDL❗✔️:   }
        // RDL❗✔️:
        // RDL❗✔️:   data->graph = gra;
        // RDL❗✔️:
        // RDL❗✔️:   return data;
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_calculate

        if graph.node_count == 0 {
            return Err(RingDecomposerError::EmptyGraph);
        }
        let cyclomatic_number =
            graph.edge_count() + connected_component_count(&graph) - graph.node_count();
        if cyclomatic_number == 0 {
            return Ok(Self {
                graph,
                urfs: Vec::new(),
                relevant_cycle_count: 0.0,
            });
        }

        let bcc = BiconnectedComponents::calculate(&graph);
        if bcc.component_count() == 0 {
            return Ok(Self {
                graph,
                urfs: Vec::new(),
                relevant_cycle_count: 0.0,
            });
        }

        let mut urfs = Vec::new();
        let mut relevant_cycle_count = 0.0;
        for component in bcc.components() {
            let mut component_graph = component.graph().clone();
            let mut shortest_paths = ShortestPathInfo::calculate(&component_graph);
            let mut cycle_families =
                CycleFamilies::calculate(&mut component_graph, &mut shortest_paths);
            if cycle_families.is_empty() {
                continue;
            }
            let urf_info =
                UrfInfo::check_urf_relation(&mut cycle_families, &component_graph, &shortest_paths);
            let marked_family_count = cycle_families
                .families()
                .iter()
                .filter(|family| family.is_relevant())
                .count();
            let grouped_family_count = urf_info.urfs.iter().map(Vec::len).sum::<usize>();
            debug_assert_eq!(marked_family_count, grouped_family_count);

            let component_relevant_cycle_count = count_relevant_cycles(
                &urf_info,
                &cycle_families,
                &component_graph,
                &shortest_paths,
            )?;
            relevant_cycle_count += component_relevant_cycle_count;
            if !relevant_cycle_count.is_finite() {
                return Err(RingDecomposerError::RelevantCycleCountOverflow);
            }
            for urf in &urf_info.urfs {
                urfs.push(unique_ring_family_from_component_urf(
                    urf,
                    &cycle_families,
                    component,
                    &shortest_paths,
                ));
            }
        }
        Ok(Self {
            graph,
            urfs,
            relevant_cycle_count,
        })
    }

    #[must_use]
    pub const fn graph(&self) -> &Graph {
        &self.graph
    }

    #[must_use]
    pub fn urfs(&self) -> &[UniqueRingFamily] {
        &self.urfs
    }

    #[must_use]
    pub fn urf_count(&self) -> usize {
        // BEGIN RDL C FUNCTION RDL_getNofURF
        // RDL❗✔️: unsigned RDL_getNofURF(const RDL_data *data)
        // RDL❗✔️: {
        // RDL❗✔️:   if (!data) {
        // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
        // RDL❗✔️:     return RDL_INVALID_RESULT;
        // RDL❗✔️:   }
        // RDL❗✔️:   return data->nofURFs;
        // RDL❗✔️: }
        // END RDL C FUNCTION RDL_getNofURF
        self.urfs.len()
    }

    #[must_use]
    pub const fn relevant_cycle_count(&self) -> f64 {
        self.relevant_cycle_count
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UniqueRingFamily {
    nodes: Vec<usize>,
    edges: Vec<EdgeId>,
}

impl UniqueRingFamily {
    #[must_use]
    pub fn nodes(&self) -> &[usize] {
        &self.nodes
    }

    #[must_use]
    pub fn edges(&self) -> &[EdgeId] {
        &self.edges
    }
}

fn canonical_edge_key(from: usize, to: usize) -> (usize, usize) {
    if from <= to { (from, to) } else { (to, from) }
}

fn connected_component_count(graph: &Graph) -> usize {
    if graph.node_count == 0 {
        return 0;
    }
    let mut seen = vec![false; graph.node_count];
    let mut count = 0usize;
    for start in 0..graph.node_count {
        if seen[start] {
            continue;
        }
        count += 1;
        seen[start] = true;
        let mut stack = vec![start];
        while let Some(node) = stack.pop() {
            for &(neighbor, _) in &graph.adjacency[node] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    count
}

fn find_cycle_families(graph: &mut Graph, shortest_paths: &mut ShortestPathInfo) -> CycleFamilies {
    // BEGIN RDL C FUNCTION RDL_findCycleFams
    // RDL✔️✔️: RDL_cfURF *RDL_findCycleFams(RDL_graph *gra, RDL_sPathInfo *spi)
    // RDL✔️✔️: {
    // RDL✔️✔️:   RDL_cfURF *rc = malloc(sizeof(*rc));
    // RDL✔️✔️:   /*number of CFs is at most 2m^2+vn (Vismara Lemma 3)*/
    // RDL✔️✔️:
    // RDL✔️✔️:   rc->alloced = 64;
    // RDL✔️✔️:   rc->fams = malloc(rc->alloced * sizeof(*rc->fams));
    // RDL✔️✔️:   rc->nofFams = 0;
    // RDL✔️✔️:   RDL_vismara(rc, gra, spi);
    let mut families = Vec::new();
    vismara_cycle_families(&mut families, graph, shortest_paths);
    // RDL✔️✔️:
    // RDL✔️✔️:   if (!rc->fams) {
    // RDL✔️✔️:     RDL_outputFunc(RDL_ERROR, "Graph is too large, can't allocate memory!\n");
    // RDL✔️✔️:     free(rc);
    // RDL✔️✔️:     return NULL;
    // RDL✔️✔️:   }
    // RDL✔️✔️:
    // RDL✔️✔️:   rc->alloced = rc->nofFams;
    // RDL✔️✔️:   rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));
    // RDL✔️✔️:
    // RDL✔️✔️:   /* sort by size */
    // RDL✔️✔️:   qsort(rc->fams, rc->nofFams, sizeof(*rc->fams), RDL_cycleFamsComp);
    families.sort_by_key(CycleFamily::weight);
    // RDL✔️✔️:   return rc;
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_findCycleFams
    CycleFamilies { families }
}

fn vismara_cycle_families(
    families: &mut Vec<CycleFamily>,
    graph: &mut Graph,
    shortest_paths: &mut ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_vismara
    // RDL✔️✔️: static void RDL_vismara(RDL_cfURF *rc, RDL_graph *gra, RDL_sPathInfo *spi)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned i,j;
    // RDL✔️✔️:   unsigned rv,yv,zv,pv,qv; /*variables as in Vismara's algorithm, extended by a 'v'*/
    // RDL✔️✔️:   unsigned *evenCand; /*'S' in Vismara's algorithm*/
    // RDL✔️✔️:   unsigned nofCandidates = 0;
    // RDL✔️✔️:   evenCand = malloc(gra->V * sizeof(*evenCand));
    // RDL✔️✔️:
    // RDL✔️✔️:   for(rv = 0; rv < gra->V; ++rv) {
    for root in 0..graph.node_count() {
        // RDL✔️✔️:     for(yv = 0; yv < gra->V; ++yv) {
        for y in 0..graph.node_count() {
            // RDL✔️✔️:       /*all yv reachable from rv respecting the ordering*/
            // RDL✔️✔️:       if(spi->reachable[rv][yv] == 1) {
            if !shortest_paths.reachable_preceding(root, y) {
                continue;
            }
            // RDL✔️✔️:         nofCandidates = 0;
            let mut even_candidates = Vec::new();
            // RDL✔️✔️:         for(i = 0; i < gra->degree[yv]; ++i) {
            for &(z, _) in &graph.adjacency[y] {
                // RDL✔️✔️:           zv = gra->adjList[yv][i][0];
                // RDL✔️✔️:           /*all zv reachable from rv respecting the ordering and adjacent to yv*/
                // RDL✔️✔️:           if(spi->reachable[rv][zv] == 1) {
                if !shortest_paths.reachable_preceding(root, z) {
                    continue;
                }
                let root_to_z = shortest_paths
                    .distance(root, z)
                    .expect("reachable vertex has distance");
                let root_to_y = shortest_paths
                    .distance(root, y)
                    .expect("reachable vertex has distance");
                // RDL✔️✔️:             if(spi->dist[rv][zv] + 1 == spi->dist[rv][yv]) {
                // RDL✔️✔️:               evenCand[nofCandidates] = zv;
                // RDL✔️✔️:               ++nofCandidates;
                // RDL✔️✔️:             }
                if root_to_z + 1 == root_to_y {
                    even_candidates.push(z);
                }
                // RDL✔️✔️:             else if(spi->dist[rv][zv] != spi->dist[rv][yv] + 1
                // RDL✔️✔️:                 && (gra->degree[zv] < gra->degree[yv] ||
                // RDL✔️✔️:                     (gra->degree[zv] == gra->degree[yv] && zv<yv))
                // RDL✔️✔️:                 && RDL_pathsShareOnlyStart(rv, yv, zv, gra, spi) == 1) {
                else if root_to_z != root_to_y + 1
                    && (graph.adjacency[z].len() < graph.adjacency[y].len()
                        || (graph.adjacency[z].len() == graph.adjacency[y].len() && z < y))
                    && paths_share_only_start(root, y, z, graph, shortest_paths)
                {
                    // RDL✔️✔️:               /*add odd cycle rv-yv rv-zv zv-yv*/
                    // RDL✔️✔️:               RDL_addOdd(rv, yv, zv, gra, spi, rc);
                    add_odd_cycle_family(families, root, y, z, graph, shortest_paths);
                }
                // RDL✔️✔️:             /*to fill dPaths in sPathInfo with the edges to r*/
                // RDL✔️✔️:             if(spi->dist[rv][zv] == 1) {
                // RDL✔️✔️:               RDL_addEdge(spi->dPaths[rv], zv, rv);
                // RDL✔️✔️:             }
                if root_to_z == 1 {
                    shortest_paths.directed_paths[root].add_directed_edge(z, root);
                }
                // RDL✔️✔️:           }
            }
            // RDL✔️✔️:         }
            // RDL✔️✔️:         /*any pair in evenCand*/
            // RDL✔️✔️:         for(i = 0; i < nofCandidates; ++i) {
            for i in 0..even_candidates.len() {
                // RDL✔️✔️:           pv = evenCand[i];
                let p = even_candidates[i];
                // RDL✔️✔️:           for(j = i+1; j < nofCandidates; ++j) {
                for &q in &even_candidates[i + 1..] {
                    // RDL✔️✔️:             qv = evenCand[j];
                    // RDL✔️✔️:             if((RDL_pathsShareOnlyStart(rv, pv, qv, gra, spi) == 1)) {
                    if paths_share_only_start(root, p, q, graph, shortest_paths) {
                        // RDL✔️✔️:               /*add even cycle rv-pv rv-qv pv-yv-qv*/
                        // RDL✔️✔️:               RDL_addEven(rv, pv, yv, qv, gra, spi, rc);
                        add_even_cycle_family(families, root, p, y, q, graph, shortest_paths);
                    }
                    // RDL✔️✔️:             }
                }
                // RDL✔️✔️:           }
            }
            // RDL✔️✔️:         }
            // RDL✔️✔️:         /*to fill dPaths in sPathInfo/fill U_r (see Vismara)*/
            // RDL✔️✔️:         for(i = 0; i < nofCandidates; ++i) {
            for p in even_candidates {
                // RDL✔️✔️:           pv = evenCand[i];
                // RDL✔️✔️:           RDL_addEdge(spi->dPaths[rv], yv, pv);
                shortest_paths.directed_paths[root].add_directed_edge(y, p);
                // RDL✔️✔️:         }
            }
            // RDL✔️✔️:       }
        }
        // RDL✔️✔️:     }
    }
    // RDL✔️✔️:   }
    // RDL✔️✔️:
    // RDL✔️✔️:   free(evenCand);
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_vismara
}

fn paths_share_only_start(
    root: usize,
    y: usize,
    z: usize,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> bool {
    // BEGIN RDL C FUNCTION RDL_pathsShareOnlyStart
    // RDL✔️✔️: static int RDL_pathsShareOnlyStart(unsigned r, unsigned y, unsigned z, RDL_graph *gra, RDL_sPathInfo *spi)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned result = 0, i, pnt, count=0;
    // RDL✔️✔️:   unsigned *vertInRY, *vertInRZ; /*edges in P(r,y) and P(r,z)*/
    let mut vertices_in_root_y = vec![false; graph.node_count()];
    let mut vertices_in_root_z = vec![false; graph.node_count()];
    // RDL✔️✔️:   vertInRY[y] = 1;
    // RDL✔️✔️:   vertInRZ[z] = 1;
    vertices_in_root_y[y] = true;
    vertices_in_root_z[z] = true;
    // RDL✔️✔️:   pnt = y;
    // RDL✔️✔️:   do
    // RDL✔️✔️:   {
    let mut point = y;
    loop {
        // RDL✔️✔️:     pnt = spi->pred[r][pnt];
        point = shortest_paths
            .predecessor(root, point)
            .expect("reachable path has predecessor");
        // RDL✔️✔️:     vertInRY[pnt] = 1;
        vertices_in_root_y[point] = true;
        // RDL✔️✔️:   }while(pnt != r);
        if point == root {
            break;
        }
    }
    // RDL✔️✔️:   pnt = z;
    // RDL✔️✔️:   do
    // RDL✔️✔️:   {
    point = z;
    loop {
        // RDL✔️✔️:     pnt = spi->pred[r][pnt];
        point = shortest_paths
            .predecessor(root, point)
            .expect("reachable path has predecessor");
        // RDL✔️✔️:     vertInRZ[pnt] = 1;
        vertices_in_root_z[point] = true;
        // RDL✔️✔️:   }while(pnt != r);
        if point == root {
            break;
        }
    }
    // RDL✔️✔️:   for(i=0; i<gra->V && count<2; ++i)
    let mut count = 0usize;
    for index in 0..graph.node_count() {
        // RDL✔️✔️:     if(vertInRY[i] == 1 && vertInRZ[i] == 1)
        // RDL✔️✔️:     {
        // RDL✔️✔️:       ++count;
        // RDL✔️✔️:     }
        if vertices_in_root_y[index] && vertices_in_root_z[index] {
            count += 1;
            if count >= 2 {
                break;
            }
        }
    }
    // RDL✔️✔️:   if(count == 1 && (vertInRY[r] == 1) && (vertInRZ[r] == 1))
    // RDL✔️✔️:   {
    // RDL✔️✔️:     result = 1;
    // RDL✔️✔️:   }
    // RDL✔️✔️:
    // RDL✔️✔️:   return result;
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_pathsShareOnlyStart
    count == 1 && vertices_in_root_y[root] && vertices_in_root_z[root]
}

fn find_prototype(
    root: usize,
    y: usize,
    z: usize,
    x: Option<usize>,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> Vec<bool> {
    // BEGIN RDL C FUNCTION RDL_findPrototype
    // RDL✔️✔️: static char *RDL_findPrototype(unsigned r, unsigned y, unsigned z, unsigned x,
    // RDL✔️✔️:     RDL_graph *gra, RDL_sPathInfo *spi)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned i, vert1, vert2;
    // RDL✔️✔️:   char *proto;
    // RDL✔️✔️:
    // RDL✔️✔️:   proto = malloc(gra->E * sizeof(*proto));
    // RDL✔️✔️:   for(i=0; i<gra->E; ++i)
    // RDL✔️✔️:   {
    // RDL✔️✔️:     proto[i] = 0;
    // RDL✔️✔️:   }
    let mut prototype = vec![false; graph.edge_count()];
    // RDL✔️✔️:   /*path from r to y*/
    mark_path_edges(&mut prototype, root, y, graph, shortest_paths);
    // RDL✔️✔️:   /*path from r to z*/
    mark_path_edges(&mut prototype, root, z, graph, shortest_paths);
    // RDL✔️✔️:   if(x == UINT_MAX)/*odd cycle*/
    if let Some(x) = x {
        // RDL✔️✔️:   else /*even cycle*/
        // RDL✔️✔️:   {
        // RDL✔️✔️:     proto[RDL_edgeId(gra,y,x)] = 1;
        // RDL✔️✔️:     proto[RDL_edgeId(gra,z,x)] = 1;
        // RDL✔️✔️:   }
        prototype[graph.edge_id(y, x).expect("cycle edge exists").index()] = true;
        prototype[graph.edge_id(z, x).expect("cycle edge exists").index()] = true;
    } else {
        // RDL✔️✔️:   {
        // RDL✔️✔️:     proto[RDL_edgeId(gra,y,z)] = 1;
        // RDL✔️✔️:   }
        prototype[graph.edge_id(y, z).expect("cycle edge exists").index()] = true;
    }
    // RDL✔️✔️:   return proto;
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_findPrototype
    prototype
}

fn mark_path_edges(
    prototype: &mut [bool],
    root: usize,
    target: usize,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) {
    let mut vertex_1 = target;
    loop {
        let vertex_2 = vertex_1;
        vertex_1 = shortest_paths
            .predecessor(root, vertex_1)
            .expect("reachable path has predecessor");
        prototype[graph
            .edge_id(vertex_1, vertex_2)
            .expect("path edge exists")
            .index()] = true;
        if vertex_1 == root {
            break;
        }
    }
}

fn add_odd_cycle_family(
    families: &mut Vec<CycleFamily>,
    root: usize,
    y: usize,
    z: usize,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_addOdd
    // RDL✔️✔️: static void RDL_addOdd(unsigned r, unsigned y, unsigned z,
    // RDL✔️✔️:     RDL_graph *gra, RDL_sPathInfo *spi, RDL_cfURF *rc)
    // RDL✔️✔️: {
    // RDL✔️✔️:   RDL_cfam *new;
    // RDL✔️✔️:   if (rc->alloced == rc->nofFams) {
    // RDL✔️✔️:     rc->alloced *= 2;
    // RDL✔️✔️:     rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));
    // RDL✔️✔️:   }
    // RDL✔️✔️:   /* skip if we're out of memory... */
    // RDL✔️✔️:   if (rc->fams) {
    // RDL✔️✔️:     new = malloc(sizeof(*new));
    // RDL✔️✔️:     new->r = r;
    // RDL✔️✔️:     new->p = y;
    // RDL✔️✔️:     new->q = z;
    // RDL✔️✔️:     new->x = UINT_MAX; /*odd cycle*/
    // RDL✔️✔️:     new->mark = 0;
    // RDL✔️✔️:     new->prototype = RDL_findPrototype(r, y, z, UINT_MAX, gra, spi);
    // RDL✔️✔️:     new->weight = spi->dist[r][y] + spi->dist[r][z] + 1;
    // RDL✔️✔️:     rc->fams[rc->nofFams++] = new;
    // RDL✔️✔️:   }
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_addOdd
    families.push(CycleFamily {
        weight: shortest_paths.distance(root, y).expect("reachable y")
            + shortest_paths.distance(root, z).expect("reachable z")
            + 1,
        r: root,
        p: y,
        q: z,
        x: None,
        prototype: find_prototype(root, y, z, None, graph, shortest_paths),
        relevant: false,
    });
}

fn add_even_cycle_family(
    families: &mut Vec<CycleFamily>,
    root: usize,
    y: usize,
    x: usize,
    z: usize,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_addEven
    // RDL✔️✔️: static void RDL_addEven(unsigned r, unsigned y, unsigned x,
    // RDL✔️✔️:     unsigned z, RDL_graph *gra, RDL_sPathInfo *spi, RDL_cfURF *rc)
    // RDL✔️✔️: {
    // RDL✔️✔️:   RDL_cfam *new;
    // RDL✔️✔️:   if (rc->alloced == rc->nofFams) {
    // RDL✔️✔️:     rc->alloced *= 2;
    // RDL✔️✔️:     rc->fams = realloc(rc->fams, rc->alloced * sizeof(*rc->fams));
    // RDL✔️✔️:   }
    // RDL✔️✔️:   /* skip if we're out of memory... */
    // RDL✔️✔️:   if (rc->fams) {
    // RDL✔️✔️:     new = malloc(sizeof(**rc->fams));
    // RDL✔️✔️:     new->r = r;
    // RDL✔️✔️:     new->p = y;
    // RDL✔️✔️:     new->q = z;
    // RDL✔️✔️:     new->x = x; /*even cycle*/
    // RDL✔️✔️:     new->mark = 0;
    // RDL✔️✔️:     new->prototype = RDL_findPrototype(r, y, z, x, gra, spi);
    // RDL✔️✔️:     new->weight = spi->dist[r][y] + spi->dist[r][z] + 2;
    // RDL✔️✔️:     rc->fams[rc->nofFams++] = new;
    // RDL✔️✔️:   }
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_addEven
    families.push(CycleFamily {
        weight: shortest_paths.distance(root, y).expect("reachable y")
            + shortest_paths.distance(root, z).expect("reachable z")
            + 2,
        r: root,
        p: y,
        q: z,
        x: Some(x),
        prototype: find_prototype(root, y, z, Some(x), graph, shortest_paths),
        relevant: false,
    });
}

fn xor_bool_slice(dst: &mut [bool], src: &[bool]) {
    // BEGIN RDL C FUNCTION RDL_bitset_xor_inplace
    // RDL✔️✔️: void RDL_bitset_xor_inplace(unsigned char* dst, unsigned const char* src, unsigned size)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned i;
    // RDL✔️✔️:
    // RDL✔️✔️:   for (i = 0; i < size; ++i) {
    // RDL✔️✔️:     dst[i] ^= src[i];
    // RDL✔️✔️:   }
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_bitset_xor_inplace
    for (dst_bit, src_bit) in dst.iter_mut().zip(src) {
        *dst_bit ^= *src_bit;
    }
}

fn swap_bool_columns(rows: &mut [Vec<bool>], left: usize, right: usize) {
    // BEGIN RDL C FUNCTION RDL_swap_columns
    // RDL✔️✔️: void RDL_swap_columns(unsigned char** rows, unsigned nof_rows, unsigned col1, unsigned col2)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned char val1, val2;
    // RDL✔️✔️:   unsigned i;
    // RDL✔️✔️:
    // RDL✔️✔️:   for (i = 0; i < nof_rows; ++i) {
    // RDL✔️✔️:     val1 = RDL_bitset_test(rows[i], col1);
    // RDL✔️✔️:     val2 = RDL_bitset_test(rows[i], col2);
    // RDL✔️✔️:
    // RDL✔️✔️:     if (!val1 != !val2) {
    // RDL✔️✔️:       RDL_bitset_flip(rows[i], col1);
    // RDL✔️✔️:       RDL_bitset_flip(rows[i], col2);
    // RDL✔️✔️:     }
    // RDL✔️✔️:   }
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_swap_columns
    for row in rows {
        if row[left] != row[right] {
            row.swap(left, right);
        }
    }
}

fn find_edges(
    edges: &mut [bool],
    root: usize,
    target: usize,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
    visited: &mut [bool],
) {
    // BEGIN RDL C FUNCTION RDL_giveEdges
    // RDL✔️✔️: void RDL_giveEdges(unsigned a, unsigned b, char *array,
    // RDL✔️✔️:     const RDL_graph *gra, const RDL_sPathInfo *spi, char* visited)
    // RDL✔️✔️: {
    // RDL✔️✔️:   unsigned i, vertex, edge;
    // RDL✔️✔️:
    // RDL✔️✔️:   if(a==b) {
    // RDL✔️✔️:     return;
    // RDL✔️✔️:   }
    if root == target {
        return;
    }
    // RDL✔️✔️:
    // RDL✔️✔️:   visited[b] = 1;
    visited[target] = true;
    // RDL✔️✔️:
    // RDL✔️✔️:   /*for each vertex adjacent to b in U_a*/
    if let Some(neighbors) = shortest_paths.directed_paths[root].neighbors(target) {
        for &vertex in neighbors {
            // RDL✔️✔️:     edge = RDL_edgeId(gra,b,vertex);
            // RDL✔️✔️:     array[edge] = 1;
            edges[graph
                .edge_id(target, vertex)
                .expect("directed path edge exists in graph")
                .index()] = true;
            // RDL✔️✔️:     if (!visited[vertex]) {
            if !visited[vertex] {
                // RDL✔️✔️:       RDL_giveEdges(a, vertex, array, gra, spi, visited);
                find_edges(edges, root, vertex, graph, shortest_paths, visited);
            }
            // RDL✔️✔️:     }
        }
    }
    // RDL✔️✔️: }
    // END RDL C FUNCTION RDL_giveEdges
}

fn make_edge_list(
    family: &CycleFamily,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> Vec<usize> {
    // BEGIN RDL C FUNCTION make_edge_list
    // RDL❗✔️: static void make_edge_list(
    // RDL❗✔️:     unsigned** edge_list,
    // RDL❗✔️:     unsigned* edge_list_size,
    // RDL❗✔️:     char* edges, RDL_graph *graph,
    // RDL❗✔️:     RDL_cfam* rcf, RDL_sPathInfo* spi)
    // RDL❗✔️: {
    // RDL❗✔️:   memset(edges, 0, graph->E * sizeof(*edges));
    let mut edge_membership = vec![false; graph.edge_count()];
    // RDL❗✔️:   RDL_findEdges(edges, rcf, graph, spi);
    find_family_edges(&mut edge_membership, family, graph, shortest_paths);
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < graph->E; ++i) {
    // RDL❗✔️:     if (edges[i]) {
    // RDL❗✔️:       (*edge_list)[*edge_list_size] = i;
    // RDL❗✔️:       ++(*edge_list_size);
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    // RDL❗✔️: }
    // END RDL C FUNCTION make_edge_list
    edge_membership
        .iter()
        .enumerate()
        .filter_map(|(edge, &contains)| contains.then_some(edge))
        .collect()
}

fn find_family_edges(
    edges: &mut [bool],
    family: &CycleFamily,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_getEdges_internal
    // RDL❗✔️: static void RDL_getEdges_internal(const RDL_data *data,
    // RDL❗✔️:     unsigned bcc_index, unsigned urf_internal_index,
    // RDL❗✔️:     unsigned rcf_internal_index, char* edges)
    // RDL❗✔️: {
    // RDL❗✔️:   const RDL_cfam **URF;
    // RDL❗✔️:   char *visited;
    // RDL❗✔️:   const RDL_graph* graph;
    // RDL❗✔️:
    // RDL❗✔️:   graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL❗✔️:
    // RDL❗✔️:   URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];
    // RDL❗✔️:
    // RDL❗✔️:   /* perform BFS on both paths and collect edges */
    // RDL❗✔️:   visited = malloc(graph->V *  sizeof(*visited));
    // RDL❗✔️:   memset(visited, 0, graph->V *  sizeof(*visited));
    // RDL❗✔️:
    // RDL❗✔️:   RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
    // RDL❗✔️:       edges, graph, data->spiPerBCC[bcc_index], visited);
    // RDL❗✔️:
    // RDL❗✔️:   memset(visited, 0, graph->V *  sizeof(*visited));
    // RDL❗✔️:   RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
    // RDL❗✔️:       edges, graph, data->spiPerBCC[bcc_index], visited);
    // RDL❗✔️:
    // RDL❗✔️:   /*
    // RDL❗✔️:    * in contrast to the nodes, we have to add additional edges, because we're
    // RDL❗✔️:    * looking at cycles...
    // RDL❗✔️:    */
    // RDL❗✔️:   if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    // RDL❗✔️:     edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->x)] = 1;
    // RDL❗✔️:     edges[RDL_edgeId(graph,URF[rcf_internal_index]->p,URF[rcf_internal_index]->x)] = 1;
    // RDL❗✔️:   }
    // RDL❗✔️:   else /*odd cycle*/ {
    // RDL❗✔️:     edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->p)] = 1;
    // RDL❗✔️:   }
    // RDL❗✔️:   free(visited);
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getEdges_internal
    //
    // BEGIN RDL C FUNCTION RDL_findEdges
    // RDL❗✔️: void RDL_findEdges(char *edges, RDL_cfam *RCF, RDL_graph *gra, RDL_sPathInfo *spi)
    // RDL❗✔️: {
    // RDL❗✔️:   char* visited;
    let mut visited = vec![false; graph.node_count()];
    // RDL❗✔️:   RDL_giveEdges(RCF->r, RCF->p, edges, gra, spi, visited);
    find_edges(
        edges,
        family.root(),
        family.p(),
        graph,
        shortest_paths,
        &mut visited,
    );
    // RDL❗✔️:   memset(visited, 0, gra->E * sizeof(*visited));
    visited.fill(false);
    // RDL❗✔️:   RDL_giveEdges(RCF->r, RCF->q, edges, gra, spi, visited);
    find_edges(
        edges,
        family.root(),
        family.q(),
        graph,
        shortest_paths,
        &mut visited,
    );
    // RDL❗✔️:   if(RCF->x < UINT_MAX) /*even family*/
    if let Some(x) = family.x() {
        // RDL❗✔️:     edges[RDL_edgeId(gra, RCF->x, RCF->p)] = 1;
        // RDL❗✔️:     edges[RDL_edgeId(gra, RCF->x, RCF->q)] = 1;
        edges[graph
            .edge_id(x, family.p())
            .expect("cycle edge exists")
            .index()] = true;
        edges[graph
            .edge_id(x, family.q())
            .expect("cycle edge exists")
            .index()] = true;
    } else {
        // RDL❗✔️:     edges[RDL_edgeId(gra, RCF->p, RCF->q)] = 1;
        edges[graph
            .edge_id(family.p(), family.q())
            .expect("cycle edge exists")
            .index()] = true;
    }
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_findEdges
}

fn sorted_edge_lists_share_edge(left: &[usize], right: &[usize]) -> bool {
    // BEGIN RDL C FUNCTION RDL_shareEdges
    // RDL❗✔️: char RDL_shareEdges(RDL_cfURF *RCFs, unsigned idx1, unsigned idx2, RDL_graph *graph,
    // RDL❗✔️:                     RDL_sPathInfo *spi)
    // RDL❗✔️: {
    // RDL❗✔️:   for(i=0; i<graph->E; ++i)
    // RDL❗✔️:   {
    // RDL❗✔️:     if(edges1[i] == 1 && edges2[i] == 1)
    // RDL❗✔️:     {
    // RDL❗✔️:       result = 1;
    // RDL❗✔️:       break;
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_shareEdges
    let mut left_index = 0usize;
    let mut right_index = 0usize;
    while left_index < left.len() && right_index < right.len() {
        match left[left_index].cmp(&right[right_index]) {
            std::cmp::Ordering::Less => left_index += 1,
            std::cmp::Ordering::Greater => right_index += 1,
            std::cmp::Ordering::Equal => return true,
        }
    }
    false
}

fn find_indegree(
    root: usize,
    target: usize,
    indegree: &mut [Option<usize>],
    shortest_paths: &ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_findIndegree
    // RDL❗✔️: static void RDL_findIndegree(unsigned a, unsigned b,
    // RDL❗✔️:     unsigned *indegree, const RDL_sPathInfo *spi)
    // RDL❗✔️: {
    // RDL❗✔️: /* similar to the function "List_Paths" from Vismara*/
    // RDL❗✔️:   unsigned i, vertex;
    // RDL❗✔️:
    // RDL❗✔️:   indegree[b] = 0;
    indegree[target] = Some(0);
    // RDL❗✔️:
    // RDL❗✔️:   if(a==b)
    // RDL❗✔️:   {
    // RDL❗✔️:     return;
    // RDL❗✔️:   }
    if root == target {
        return;
    }
    // RDL❗✔️:   /*for each vertex adjacent to b in U_a*/
    // RDL❗✔️:   for(i=0; i<spi->dPaths[a]->degree[b]; ++i)
    // RDL❗✔️:   {
    for &vertex in shortest_paths.directed_paths[root]
        .neighbors(target)
        .unwrap_or_default()
    {
        // RDL❗✔️:     vertex = spi->dPaths[a]->adjList[b][i][0];
        // RDL❗✔️:     if (indegree[vertex] == UINT_MAX) {
        // RDL❗✔️:       RDL_findIndegree(a, vertex, indegree, spi);
        // RDL❗✔️:     }
        if indegree[vertex].is_none() {
            find_indegree(root, vertex, indegree, shortest_paths);
        }
        // RDL❗✔️:     indegree[vertex] += 1;
        let degree = indegree[vertex].expect("recursive visit initializes indegree");
        indegree[vertex] = Some(degree + 1);
        // RDL❗✔️:   }
    }
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_findIndegree
}

fn count_paths(
    root: usize,
    target: usize,
    node_count: usize,
    shortest_paths: &ShortestPathInfo,
) -> Result<f64, RingDecomposerError> {
    // BEGIN RDL C FUNCTION RDL_countPaths
    // RDL❗✔️: double RDL_countPaths(unsigned a, unsigned b,
    // RDL❗✔️:     unsigned V, const RDL_sPathInfo *spi)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned i, curr, vertex, empty_end=0;
    // RDL❗✔️:   unsigned *indegree, *empty_set;
    // RDL❗✔️:   double *nof_paths, result;
    // RDL❗✔️:
    // RDL❗✔️:   nof_paths = malloc(V * sizeof(*nof_paths));
    // RDL❗✔️:   indegree = malloc(V * sizeof(*indegree));
    // RDL❗✔️:   empty_set = malloc(V * sizeof(*empty_set));
    // RDL❗✔️:   for (i = 0; i < V; ++i) {
    // RDL❗✔️:     nof_paths[i] = 0.0;
    // RDL❗✔️:     indegree[i] = UINT_MAX;
    // RDL❗✔️:   }
    let mut number_of_paths = vec![0.0; node_count];
    let mut indegree = vec![None; node_count];
    let mut empty_set = Vec::new();
    // RDL❗✔️:
    // RDL❗✔️:   RDL_findIndegree(a, b, indegree, spi);
    find_indegree(root, target, &mut indegree, shortest_paths);
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < V; ++i) {
    // RDL❗✔️:     if (indegree[i] == 0) {
    // RDL❗✔️:       empty_set[empty_end] = i;
    // RDL❗✔️:       ++empty_end;
    // RDL❗✔️:       nof_paths[i] = 1.0;
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    for (node, degree) in indegree.iter().enumerate() {
        if *degree == Some(0) {
            empty_set.push(node);
            number_of_paths[node] = 1.0;
        }
    }
    // RDL❗✔️:
    // RDL❗✔️:   if (empty_end != 1 || empty_set[0] != b) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "invalid topological sort!");
    // RDL❗✔️:     /* should never happen */
    // RDL❗✔️:     assert(0);
    // RDL❗✔️:     free(nof_paths);
    // RDL❗✔️:     free(indegree);
    // RDL❗✔️:     free(empty_set);
    // RDL❗✔️:     return DBL_MAX;
    // RDL❗✔️:   }
    if empty_set.len() != 1 || empty_set[0] != target {
        return Err(RingDecomposerError::InvalidShortestPathDag { root, target });
    }
    // RDL❗✔️:
    // RDL❗✔️:   /* look at the nodes in an topological sort */
    // RDL❗✔️:   while (empty_end) {
    while let Some(current) = empty_set.pop() {
        // RDL❗✔️:     curr = empty_set[empty_end-1];
        // RDL❗✔️:     --empty_end;
        // RDL❗✔️:     for(i=0; i<spi->dPaths[a]->degree[curr]; ++i) {
        for &vertex in shortest_paths.directed_paths[root]
            .neighbors(current)
            .unwrap_or_default()
        {
            // RDL❗✔️:       vertex = spi->dPaths[a]->adjList[curr][i][0];
            // RDL❗✔️:       --indegree[vertex];
            let degree = indegree[vertex]
                .ok_or(RingDecomposerError::InvalidShortestPathDag { root, target })?;
            let next_degree = degree
                .checked_sub(1)
                .ok_or(RingDecomposerError::InvalidShortestPathDag { root, target })?;
            indegree[vertex] = Some(next_degree);
            // RDL❗✔️:       nof_paths[vertex] += nof_paths[curr];
            number_of_paths[vertex] += number_of_paths[current];
            // RDL❗✔️:       if (!indegree[vertex]) {
            // RDL❗✔️:         empty_set[empty_end] = vertex;
            // RDL❗✔️:         ++empty_end;
            // RDL❗✔️:       }
            if next_degree == 0 {
                empty_set.push(vertex);
            }
            // RDL❗✔️:     }
        }
        // RDL❗✔️:   }
    }
    // RDL❗✔️:
    // RDL❗✔️:   result = nof_paths[a];
    // RDL❗✔️:
    // RDL❗✔️:   free(nof_paths);
    // RDL❗✔️:   free(indegree);
    // RDL❗✔️:   free(empty_set);
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_countPaths
    Ok(number_of_paths[root])
}

fn count_relevant_cycles_for_family(
    family: &CycleFamily,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> Result<f64, RingDecomposerError> {
    // BEGIN RDL C FUNCTION RDL_getNofRCForRCF_internal
    // RDL❗✔️: static double RDL_getNofRCForRCF_internal(const RDL_data *data,
    // RDL❗✔️:     unsigned bcc_index, unsigned urf_internal_index,
    // RDL❗✔️:     unsigned rcf_internal_index)
    // RDL❗✔️: {
    // RDL❗✔️:   double nofPaths1, nofPaths2, result;
    // RDL❗✔️:   const double prod_limit = sqrt(DBL_MAX) - 1.0;
    // RDL❗✔️:   const RDL_cfam** URF;
    // RDL❗✔️:   const RDL_graph* graph;
    // RDL❗✔️:
    // RDL❗✔️:   result = 0.0;
    // RDL❗✔️:
    // RDL❗✔️:   graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL❗✔️:
    // RDL❗✔️:   URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];
    // RDL❗✔️:
    // RDL❗✔️:   /* count the number of paths one both sides */
    // RDL❗✔️:   nofPaths1 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
    // RDL❗✔️:       graph->V, data->spiPerBCC[bcc_index]);
    // RDL❗✔️:   nofPaths2 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
    // RDL❗✔️:       graph->V, data->spiPerBCC[bcc_index]);
    let number_of_paths_1 = count_paths(
        family.root(),
        family.q(),
        graph.node_count(),
        shortest_paths,
    )?;
    let number_of_paths_2 = count_paths(
        family.root(),
        family.p(),
        graph.node_count(),
        shortest_paths,
    )?;
    // RDL❗✔️:
    // RDL❗✔️:   result = nofPaths1 * nofPaths2;
    // RDL❗✔️:
    // RDL❗✔️:   /* check if either of the number paths is larger than what we can multiply */
    // RDL❗✔️:   if (nofPaths1 >= prod_limit || nofPaths2 >= prod_limit) {
    // RDL❗✔️:     RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
    // RDL❗✔️:     return RDL_INVALID_RC_COUNT;
    // RDL❗✔️:   }
    let product_limit = f64::MAX.sqrt() - 1.0;
    if number_of_paths_1 >= product_limit || number_of_paths_2 >= product_limit {
        return Err(RingDecomposerError::RelevantCycleCountOverflow);
    }
    let result = number_of_paths_1 * number_of_paths_2;
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNofRCForRCF_internal
    Ok(result)
}

fn count_relevant_cycles_for_urf(
    urf: &[usize],
    cycle_families: &CycleFamilies,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> Result<f64, RingDecomposerError> {
    // BEGIN RDL C FUNCTION RDL_getNofRCForURF
    // RDL❗✔️: double RDL_getNofRCForURF(const RDL_data *data, unsigned index)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned i,nofFams, bcc_index, internal_index;
    // RDL❗✔️:   double result=0, prod;
    // RDL❗✔️:
    // RDL❗✔️:   const double sum_limit = DBL_MAX/2.0 - 1.0;
    // RDL❗✔️:
    // RDL❗✔️:   if (!data) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    // RDL❗✔️:     return RDL_INVALID_RC_COUNT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   if (index >= data->nofURFs) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    // RDL❗✔️:     return RDL_INVALID_RC_COUNT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   /* map indices */
    // RDL❗✔️:   bcc_index = data->urf_to_bcc[index][0];
    // RDL❗✔️:   internal_index = data->urf_to_bcc[index][1];
    // RDL❗✔️:
    // RDL❗✔️:   nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
    // RDL❗✔️:
    // RDL❗✔️:   /* iterate all RCFs in this URF */
    // RDL❗✔️:   for(i=0; i<nofFams; ++i) {
    let mut result = 0.0;
    let sum_limit = f64::MAX / 2.0 - 1.0;
    for &family_index in urf {
        // RDL❗✔️:     prod = RDL_getNofRCForRCF_internal(data, bcc_index, internal_index, i);
        let product = count_relevant_cycles_for_family(
            &cycle_families.families()[family_index],
            graph,
            shortest_paths,
        )?;
        // RDL❗✔️:
        // RDL❗✔️:     /* check if any summand is larger than limit */
        // RDL❗✔️:     if (prod >= sum_limit || result >= sum_limit) {
        // RDL❗✔️:       RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
        // RDL❗✔️:       return RDL_INVALID_RC_COUNT;
        // RDL❗✔️:     }
        if product >= sum_limit || result >= sum_limit {
            return Err(RingDecomposerError::RelevantCycleCountOverflow);
        }
        // RDL❗✔️:
        // RDL❗✔️:     result += prod;
        result += product;
        // RDL❗✔️:   }
    }
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNofRCForURF
    Ok(result)
}

fn count_relevant_cycles(
    urf_info: &UrfInfo,
    cycle_families: &CycleFamilies,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) -> Result<f64, RingDecomposerError> {
    // BEGIN RDL C FUNCTION RDL_getNofRC
    // RDL❗✔️: double RDL_getNofRC(const RDL_data *data)
    // RDL❗✔️: {
    // RDL❗✔️:   double result = 0, intermediate;
    // RDL❗✔️:   unsigned i;
    // RDL❗✔️:
    // RDL❗✔️:   if (!data) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    // RDL❗✔️:     return RDL_INVALID_RC_COUNT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   /* add up results for URFs */
    // RDL❗✔️:   for (i = 0; i < data->nofURFs; ++i) {
    let mut result = 0.0;
    for urf in &urf_info.urfs {
        // RDL❗✔️:     intermediate = RDL_getNofRCForURF(data, i);
        let intermediate =
            count_relevant_cycles_for_urf(urf, cycle_families, graph, shortest_paths)?;
        // RDL❗✔️:     if (intermediate == RDL_INVALID_RC_COUNT) {
        // RDL❗✔️:       return RDL_INVALID_RC_COUNT;
        // RDL❗✔️:     }
        // RDL❗✔️:
        // RDL❗✔️:     result += intermediate;
        result += intermediate;
        if !result.is_finite() {
            return Err(RingDecomposerError::RelevantCycleCountOverflow);
        }
        // RDL❗✔️:   }
    }
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNofRC
    Ok(result)
}

fn find_vertices(
    vertices: &mut [bool],
    root: usize,
    target: usize,
    shortest_paths: &ShortestPathInfo,
    visited: &mut [bool],
) {
    // BEGIN RDL C FUNCTION RDL_giveVertices
    // RDL❗✔️: void RDL_giveVertices(unsigned a, unsigned b, char *array,
    // RDL❗✔️:     const RDL_sPathInfo *spi, char *visited)
    // RDL❗✔️: {
    // RDL❗✔️:  /* similar to the function "List_Paths" from Vismara*/
    // RDL❗✔️:   unsigned i, vertex;
    // RDL❗✔️:
    // RDL❗✔️:   visited[b] = 1;
    visited[target] = true;
    // RDL❗✔️:
    // RDL❗✔️:   if(a==b)
    // RDL❗✔️:   {
    // RDL❗✔️:     array[a] = 1;
    // RDL❗✔️:     return;
    // RDL❗✔️:   }
    if root == target {
        vertices[root] = true;
        return;
    }
    // RDL❗✔️:   array[b] = 1;
    vertices[target] = true;
    // RDL❗✔️:   /*for each vertex adjacent to b in U_a*/
    // RDL❗✔️:   for(i=0; i<spi->dPaths[a]->degree[b]; ++i)
    // RDL❗✔️:   {
    if let Some(neighbors) = shortest_paths.directed_paths[root].neighbors(target) {
        for &vertex in neighbors {
            // RDL❗✔️:     vertex = spi->dPaths[a]->adjList[b][i][0];
            // RDL❗✔️:     if (!visited[vertex]) {
            if !visited[vertex] {
                // RDL❗✔️:       RDL_giveVertices(a, vertex, array, spi, visited);
                find_vertices(vertices, root, vertex, shortest_paths, visited);
                // RDL❗✔️:     }
            }
            // RDL❗✔️:   }
        }
    }
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_giveVertices
}

fn find_family_nodes(
    nodes: &mut [bool],
    family: &CycleFamily,
    graph: &Graph,
    shortest_paths: &ShortestPathInfo,
) {
    // BEGIN RDL C FUNCTION RDL_getNodes_internal
    // RDL❗✔️: static void RDL_getNodes_internal(const RDL_data *data, unsigned bcc_index,
    // RDL❗✔️:     unsigned urf_internal_index, unsigned rcf_internal_index, char* atoms)
    // RDL❗✔️: {
    // RDL❗✔️:   const RDL_cfam **URF;
    // RDL❗✔️:   char *visited;
    // RDL❗✔️:   const RDL_graph* graph;
    // RDL❗✔️:
    // RDL❗✔️:   graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL❗✔️:
    // RDL❗✔️:   URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];
    // RDL❗✔️:
    // RDL❗✔️:   visited = malloc(graph->V *  sizeof(*visited));
    let mut visited = vec![false; graph.node_count()];
    // RDL❗✔️:
    // RDL❗✔️:   /* BFS on both shortest path graphs */
    // RDL❗✔️:   memset(visited, 0, graph->V *  sizeof(*visited));
    // RDL❗✔️:   RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
    // RDL❗✔️:       atoms, data->spiPerBCC[bcc_index], visited);
    find_vertices(
        nodes,
        family.root(),
        family.q(),
        shortest_paths,
        &mut visited,
    );
    // RDL❗✔️:
    // RDL❗✔️:   memset(visited, 0, graph->V *  sizeof(*visited));
    // RDL❗✔️:   RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
    // RDL❗✔️:       atoms, data->spiPerBCC[bcc_index], visited);
    visited.fill(false);
    find_vertices(
        nodes,
        family.root(),
        family.p(),
        shortest_paths,
        &mut visited,
    );
    // RDL❗✔️:
    // RDL❗✔️:   if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    // RDL❗✔️:     atoms[URF[rcf_internal_index]->x] = 1;
    // RDL❗✔️:   }
    if let Some(x) = family.x() {
        nodes[x] = true;
    }
    // RDL❗✔️:   free(visited);
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNodes_internal
}

fn unique_ring_family_from_component_urf(
    urf: &[usize],
    cycle_families: &CycleFamilies,
    component: &BiconnectedComponent,
    shortest_paths: &ShortestPathInfo,
) -> UniqueRingFamily {
    // BEGIN RDL C FUNCTION RDL_getNodesURF
    // RDL❗✔️: static RDL_node *RDL_getNodesURF(const RDL_data *data, unsigned index)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned i,nofFams,nextfree=0,alloced, bcc_index, internal_index;
    // RDL❗✔️:   char *atoms;
    // RDL❗✔️:   RDL_node *result;
    // RDL❗✔️:   const RDL_graph* graph;
    // RDL❗✔️:
    // RDL❗✔️:   bcc_index = data->urf_to_bcc[index][0];
    // RDL❗✔️:   internal_index = data->urf_to_bcc[index][1];
    // RDL❗✔️:   graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL❗✔️:   atoms = malloc(graph->V * sizeof(*atoms));
    // RDL❗✔️:   memset(atoms, 0, graph->V * sizeof(*atoms));
    // RDL❗✔️:
    // RDL❗✔️:   nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
    // RDL❗✔️:   alloced = RDL_RESERVED_START;
    // RDL❗✔️:   result = malloc(alloced * sizeof(*result));
    // RDL❗✔️:
    // RDL❗✔️:   /* iterate over all RCFs in this URFs and collect data in one bitset (atoms) */
    // RDL❗✔️:   for(i=0; i<nofFams; ++i) {
    // RDL❗✔️:     RDL_getNodes_internal(data, bcc_index, internal_index, i, atoms);
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   /* translate into a dynamic table */
    // RDL❗✔️:   for(i=0; i<graph->V; ++i) {
    // RDL❗✔️:     if(atoms[i] == 1) {
    // RDL❗✔️:       if(nextfree == alloced) {/*double the size*/
    // RDL❗✔️:         alloced *= 2;
    // RDL❗✔️:         result = realloc(result, alloced*sizeof(*result));
    // RDL❗✔️:       }
    // RDL❗✔️:       result[nextfree++] = data->bccGraphs->node_from_bcc_mapping[bcc_index][i];
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   result = realloc(result, (nextfree+1)*sizeof(*result));
    // RDL❗✔️:
    // RDL❗✔️:   result[nextfree] = UINT_MAX;
    // RDL❗✔️:   free(atoms);
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNodesURF
    //
    // BEGIN RDL C FUNCTION RDL_getNodesForURF
    // RDL❗✔️: unsigned RDL_getNodesForURF(const RDL_data *data, unsigned index, RDL_node **ptr)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned i;
    // RDL❗✔️:
    // RDL❗✔️:   if (!data) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    // RDL❗✔️:     (*ptr) = malloc(sizeof(**ptr));
    // RDL❗✔️:     return RDL_INVALID_RESULT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   if (index >= data->nofURFs) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    // RDL❗✔️:     (*ptr) = malloc(sizeof(**ptr));
    // RDL❗✔️:     return RDL_INVALID_RESULT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   (*ptr) = RDL_getNodesURF(data, index);
    // RDL❗✔️:   for(i=0; (*ptr)[i]<UINT_MAX; ++i); /*counts the number of atoms*/
    // RDL❗✔️:   return i;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getNodesForURF
    //
    // BEGIN RDL C FUNCTION RDL_getEdgesURF
    // RDL❗✔️: static unsigned *RDL_getEdgesURF(const RDL_data *data, unsigned index)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned i,nofFams,nextfree=0,alloced, bcc_index, internal_index;
    // RDL❗✔️:   char *edges;
    // RDL❗✔️:   unsigned *result;
    // RDL❗✔️:   const RDL_graph* graph;
    // RDL❗✔️:
    // RDL❗✔️:   /* map indices */
    // RDL❗✔️:   bcc_index = data->urf_to_bcc[index][0];
    // RDL❗✔️:   internal_index = data->urf_to_bcc[index][1];
    // RDL❗✔️:   graph = data->bccGraphs->bcc_graphs[bcc_index];
    // RDL❗✔️:   edges = malloc(graph->E * sizeof(*edges));
    // RDL❗✔️:   memset(edges, 0, graph->E * sizeof(*edges));
    // RDL❗✔️:
    // RDL❗✔️:   nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
    // RDL❗✔️:   alloced = RDL_RESERVED_START;
    // RDL❗✔️:   result = malloc(alloced * sizeof(*result));
    // RDL❗✔️:
    // RDL❗✔️:   /* iterate over all RCFs in this URF and collect results in bitset (edges) */
    // RDL❗✔️:   for(i=0; i<nofFams; ++i) {
    // RDL❗✔️:     RDL_getEdges_internal(data, bcc_index, internal_index, i, edges);
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   /* translate into dynamic table */
    // RDL❗✔️:   for(i=0; i<graph->E; ++i) {
    // RDL❗✔️:     if(edges[i] == 1) {
    // RDL❗✔️:       if(nextfree == alloced)
    // RDL❗✔️:       {/*double the size*/
    // RDL❗✔️:         alloced *= 2;
    // RDL❗✔️:         result = realloc(result, alloced * sizeof(*result));
    // RDL❗✔️:       }
    // RDL❗✔️:       result[nextfree++] = data->bccGraphs->edge_from_bcc_mapping[bcc_index][i];
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    // RDL❗✔️:   result = realloc(result, (nextfree+1)*sizeof(*result));
    // RDL❗✔️:
    // RDL❗✔️:   result[nextfree] = UINT_MAX;
    // RDL❗✔️:
    // RDL❗✔️:   free(edges);
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getEdgesURF
    //
    // BEGIN RDL C FUNCTION RDL_getEdgesForURF
    // RDL❗✔️: unsigned RDL_getEdgesForURF(const RDL_data *data, unsigned index, RDL_edge **ptr)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned nextfree, alloced;
    // RDL❗✔️:   RDL_edge *result;
    // RDL❗✔️:   unsigned *edgeIndices;
    // RDL❗✔️:   unsigned i;
    // RDL❗✔️:
    // RDL❗✔️:   if (!data) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    // RDL❗✔️:     (*ptr) = malloc(sizeof(**ptr));
    // RDL❗✔️:     return RDL_INVALID_RESULT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   if (index >= data->nofURFs) {
    // RDL❗✔️:     RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    // RDL❗✔️:     (*ptr) = malloc(sizeof(**ptr));
    // RDL❗✔️:     return RDL_INVALID_RESULT;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   nextfree = 0;
    // RDL❗✔️:   alloced = RDL_RESERVED_START;
    // RDL❗✔️:   result = malloc(alloced * sizeof(*result));
    // RDL❗✔️:   edgeIndices = RDL_getEdgesURF(data, index);
    // RDL❗✔️:   /* translate result from edge indices to edges */
    // RDL❗✔️:   for(i=0; edgeIndices[i]<UINT_MAX; ++i)
    // RDL❗✔️:   {
    // RDL❗✔️:     if(nextfree == alloced)/*more space needed in result*/
    // RDL❗✔️:     {
    // RDL❗✔️:       alloced *= 2; /* double the space */
    // RDL❗✔️:       result = realloc(result, alloced * sizeof(*result));
    // RDL❗✔️:     }
    // RDL❗✔️:     result[nextfree][0] = data->graph->edges[edgeIndices[i]][0];
    // RDL❗✔️:     result[nextfree][1] = data->graph->edges[edgeIndices[i]][1];
    // RDL❗✔️:     ++nextfree;
    // RDL❗✔️:   }
    // RDL❗✔️:   result = realloc(result, nextfree * sizeof(*result));
    // RDL❗✔️:   free(edgeIndices);
    // RDL❗✔️:   (*ptr) = result;
    // RDL❗✔️:   return nextfree;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_getEdgesForURF

    let mut edge_membership = vec![false; component.graph().edge_count()];
    let mut node_membership = vec![false; component.original_nodes().len()];
    for &family_index in urf {
        let family = &cycle_families.families()[family_index];
        find_family_edges(
            &mut edge_membership,
            family,
            component.graph(),
            shortest_paths,
        );
        find_family_nodes(
            &mut node_membership,
            family,
            component.graph(),
            shortest_paths,
        );
    }
    let mut edges = Vec::new();
    for (local_edge, &contains) in edge_membership.iter().enumerate() {
        if !contains {
            continue;
        }
        edges.push(component.original_edges()[local_edge]);
    }
    let nodes = node_membership
        .iter()
        .enumerate()
        .filter_map(|(local_node, &contains)| {
            contains.then_some(component.original_nodes()[local_node])
        })
        .collect();
    UniqueRingFamily { nodes, edges }
}

fn tarjan_biconnected_components(graph: &Graph) -> BiconnectedComponents {
    // BEGIN RDL C FUNCTION RDL_tarjanBCC
    // RDL❗✔️: RDL_BCCGraph* RDL_tarjanBCC(const RDL_graph* graph)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned *bcc, u, node, i, j, k, ii[2],
    // RDL❗✔️:     nof_bcc, *bcc_size, nof_nontrivial_bcc=0, curr_bcc, found;
    // RDL❗✔️:   unsigned *non_trivial_mapping;
    // RDL❗✔️:   RDL_BCCGraph* result;
    // RDL❗✔️:
    // RDL❗✔️:   nof_bcc = RDL_tarjan(graph, &bcc);
    let edge_components = tarjan_edge_components(graph);
    // RDL❗✔️:
    // RDL❗✔️:   result = malloc(sizeof(*result));
    // RDL❗✔️:   bcc_size = malloc(nof_bcc * sizeof(*bcc_size));
    let max_component = edge_components.iter().copied().max().unwrap_or(0);
    let mut component_sizes = vec![0usize; max_component + 1];
    // RDL❗✔️:   non_trivial_mapping = malloc(nof_bcc * sizeof(*non_trivial_mapping));
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < nof_bcc; ++i) {
    // RDL❗✔️:     bcc_size[i] = 0;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < graph->E; ++i) {
    for &component in &edge_components {
        if component != 0 {
            // RDL❗✔️:     curr_bcc = bcc[i];
            // RDL❗✔️:     ++bcc_size[curr_bcc - 1];
            component_sizes[component] += 1;
        }
    }
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < nof_bcc; ++i) {
    // RDL❗✔️:     if(bcc_size[i] > 1) {
    // RDL❗✔️:       non_trivial_mapping[i] = nof_nontrivial_bcc;
    // RDL❗✔️:       ++nof_nontrivial_bcc;
    // RDL❗✔️:     }
    // RDL❗✔️:     else {
    // RDL❗✔️:       non_trivial_mapping[i] = RDL_NO_RINGSYSTEM;
    // RDL❗✔️:     }
    // RDL❗✔️:   }
    let mut non_trivial_mapping = vec![None; max_component + 1];
    let mut next_component = 0usize;
    for component in 1..=max_component {
        if component_sizes[component] > 1 {
            non_trivial_mapping[component] = Some(next_component);
            next_component += 1;
        }
    }

    let mut original_edges_by_component = vec![Vec::new(); next_component];
    let mut original_nodes_by_component = vec![Vec::<usize>::new(); next_component];
    let mut local_node_by_original_and_component = BTreeMap::<(usize, usize), usize>::new();
    let mut edge_to_component = vec![None; graph.edge_count()];
    let mut node_to_components = vec![Vec::<(usize, usize)>::new(); graph.node_count()];

    // RDL❗✔️:   for (i = 0; i < graph->E; ++i) {
    for (edge_index, &raw_component) in edge_components.iter().enumerate() {
        // RDL❗✔️:     curr_bcc = bcc[i];
        let Some(component_index) = non_trivial_mapping.get(raw_component).copied().flatten()
        else {
            // RDL❗✔️:     /* skip trivial BCCs */
            // RDL❗✔️:     if (bcc_size[curr_bcc-1] <= 1) {
            // RDL❗✔️:       continue;
            // RDL❗✔️:     }
            continue;
        };
        // RDL❗✔️:
        // RDL❗✔️:     /* -1 because of 1-based BCCs */
        // RDL❗✔️:     u = non_trivial_mapping[curr_bcc - 1];
        let local_edge_index = original_edges_by_component[component_index].len();
        edge_to_component[edge_index] = Some((component_index, local_edge_index));
        original_edges_by_component[component_index].push(EdgeId::new(edge_index));
        let edge = graph.edges[edge_index];
        // RDL❗✔️:
        // RDL❗✔️:     for (k = 0; k <= 1; ++k) {
        for node in [edge.from(), edge.to()] {
            // RDL❗✔️:       node = graph->edges[i][k];
            // RDL❗✔️:
            // RDL❗✔️:       found = 0;
            if local_node_by_original_and_component.contains_key(&(component_index, node)) {
                continue;
            }
            // RDL❗✔️:       if (!found) {
            // RDL❗✔️:         ++result->nof_bcc_per_node[node];
            let local_node_index = original_nodes_by_component[component_index].len();
            local_node_by_original_and_component.insert((component_index, node), local_node_index);
            node_to_components[node].push((component_index, local_node_index));
            // RDL❗✔️:         result->node_from_bcc_mapping[u][result->nof_nodes_per_bcc[u]-1] = node;
            original_nodes_by_component[component_index].push(node);
        }
        // RDL❗✔️:     }
    }
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   for (u = 0; u < nof_nontrivial_bcc; ++u) {
    // RDL❗✔️:     result->bcc_graphs[u] = RDL_initNewGraph(result->nof_nodes_per_bcc[u]);
    // RDL❗✔️:   }
    let mut component_graphs = original_nodes_by_component
        .iter()
        .map(|nodes| Graph::new(nodes.len()))
        .collect::<Vec<_>>();
    // RDL❗✔️:
    // RDL❗✔️:   for (i = 0; i < graph->E; ++i) {
    for (edge_index, component_mapping) in edge_to_component.iter().copied().enumerate() {
        // RDL❗✔️:     u = result->edge_to_bcc_mapping[i][0];
        let Some((component_index, _)) = component_mapping else {
            // RDL❗✔️:
            // RDL❗✔️:     /* skip if there is no BCC */
            // RDL❗✔️:     if (u == RDL_NO_RINGSYSTEM) {
            // RDL❗✔️:       continue;
            // RDL❗✔️:     }
            continue;
        };
        let edge = graph.edges[edge_index];
        let left = local_node_by_original_and_component[&(component_index, edge.from())];
        let right = local_node_by_original_and_component[&(component_index, edge.to())];
        // RDL❗✔️:     RDL_addUEdge(result->bcc_graphs[u], ii[0], ii[1]);
        let _ = component_graphs[component_index].add_undirected_edge(left, right);
    }
    // RDL❗✔️:
    // RDL❗✔️:   return result;
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_tarjanBCC

    let components = component_graphs
        .into_iter()
        .enumerate()
        .map(|(index, graph)| BiconnectedComponent {
            graph,
            original_nodes: std::mem::take(&mut original_nodes_by_component[index]),
            original_edges: std::mem::take(&mut original_edges_by_component[index]),
        })
        .collect();

    BiconnectedComponents {
        components,
        edge_to_component,
        node_to_components,
    }
}

fn tarjan_edge_components(graph: &Graph) -> Vec<usize> {
    // BEGIN RDL C FUNCTION RDL_tarjan
    // RDL❗✔️: static unsigned RDL_tarjan(const RDL_graph* graph,
    // RDL❗✔️:                     unsigned** bcc)
    // RDL❗✔️: {
    // RDL❗✔️:   unsigned *d, *low, u, time, curr_bcc, i;
    // RDL❗✔️:   RDL_stack* edge_stack;
    // RDL❗✔️:
    // RDL❗✔️:   d = malloc(graph->V * sizeof(*d));
    // RDL❗✔️:   low = malloc(graph->V * sizeof(*low));
    // RDL❗✔️:   *bcc = malloc(graph->E * sizeof(**bcc));
    let mut discovery = vec![0usize; graph.node_count()];
    let mut low = vec![0usize; graph.node_count()];
    let mut edge_components = vec![0usize; graph.edge_count()];
    let mut edge_stack = Vec::<EdgeId>::new();
    let mut time = 0usize;
    let mut current_component = 0usize;
    // RDL❗✔️:   for (i = 0; i < graph->E; ++i) {
    // RDL❗✔️:     (*bcc)[i] = 0;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   edge_stack = RDL_stack_new();
    // RDL❗✔️:
    // RDL❗✔️:   time = 0;
    // RDL❗✔️:   curr_bcc = 1;
    // RDL❗✔️:
    // RDL❗✔️:   for (u = 0; u < graph->V; ++u) {
    // RDL❗✔️:     d[u] = 0;
    // RDL❗✔️:     low[u] = 0;
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   for (u = 0; u < graph->V; ++u) {
    for node in 0..graph.node_count() {
        // RDL❗✔️:     if (d[u] == 0) {
        // RDL❗✔️:       RDL_tarjanVisit(graph, u, d, low, &time, &curr_bcc, *bcc, edge_stack);
        // RDL❗✔️:     }
        if discovery[node] == 0 {
            tarjan_visit(
                graph,
                node,
                None,
                &mut discovery,
                &mut low,
                &mut time,
                &mut current_component,
                &mut edge_stack,
                &mut edge_components,
            );
        }
    }
    // RDL❗✔️:   }
    // RDL❗✔️:
    // RDL❗✔️:   return (curr_bcc-1);
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_tarjan
    edge_components
}

#[allow(clippy::too_many_arguments)]
fn tarjan_visit(
    graph: &Graph,
    node: usize,
    parent: Option<usize>,
    discovery: &mut [usize],
    low: &mut [usize],
    time: &mut usize,
    current_component: &mut usize,
    edge_stack: &mut Vec<EdgeId>,
    edge_components: &mut [usize],
) {
    // BEGIN RDL C FUNCTION RDL_tarjanVisit
    // RDL❗✔️: static void RDL_tarjanVisit(const RDL_graph* graph, unsigned start, unsigned* d,
    // RDL❗✔️:     unsigned* low, unsigned* time, unsigned* curr_bcc, unsigned* bcc,
    // RDL❗✔️:     RDL_stack* edge_stack)
    // RDL❗✔️: {
    // RDL❗✔️:   ...
    // RDL❗✔️:   ++(*time);
    // RDL❗✔️:   d[start] = *time;
    // RDL❗✔️:   low[start] = d[start];
    *time += 1;
    discovery[node] = *time;
    low[node] = discovery[node];
    // RDL❗✔️:
    // RDL❗✔️:   while (!RDL_stack_empty(dfs_stack)) {
    for &(neighbor, edge_id) in &graph.adjacency[node] {
        // RDL❗✔️:       v = graph->adjList[u][j][0];
        // RDL❗✔️:       uv_edge = RDL_edgeId(graph, u, v);
        // RDL❗✔️:
        // RDL❗✔️:       if (d[v] == 0) {
        if discovery[neighbor] == 0 {
            // RDL❗✔️:         edge_elements[next_free_edge_element] = uv_edge;
            // RDL❗✔️:         RDL_stack_push(edge_stack, &(edge_elements[next_free_edge_element]));
            edge_stack.push(edge_id);
            tarjan_visit(
                graph,
                neighbor,
                Some(node),
                discovery,
                low,
                time,
                current_component,
                edge_stack,
                edge_components,
            );
            // RDL❗✔️:         low[u] = low[u] < low[v] ? low[u] : low[v];
            low[node] = low[node].min(low[neighbor]);
            // RDL❗✔️:         if (low[v] >= d[u]) {
            if low[neighbor] >= discovery[node] {
                *current_component += 1;
                // RDL❗✔️:           do {
                loop {
                    // RDL❗✔️:             element = RDL_stack_top(edge_stack);
                    // RDL❗✔️:             curr_edge = *element;
                    // RDL❗✔️:             RDL_stack_pop(edge_stack);
                    let Some(component_edge) = edge_stack.pop() else {
                        break;
                    };
                    // RDL❗✔️:
                    // RDL❗✔️:             bcc[curr_edge] = *curr_bcc;
                    edge_components[component_edge.index()] = *current_component;
                    // RDL❗✔️:           } while (curr_edge != uv_edge);
                    if component_edge == edge_id {
                        break;
                    }
                }
                // RDL❗✔️:           ++(*curr_bcc);
            }
            // RDL❗✔️:         }
            // RDL❗✔️:       }
            // RDL❗✔️:       else if (d[v] < d[u] && v != parent) {
        } else if Some(neighbor) != parent && discovery[neighbor] < discovery[node] {
            // RDL❗✔️:         edge_elements[next_free_edge_element] = uv_edge;
            // RDL❗✔️:         RDL_stack_push(edge_stack, &(edge_elements[next_free_edge_element]));
            edge_stack.push(edge_id);
            // RDL❗✔️:         low[u] = low[u] < d[v] ? low[u] : d[v];
            low[node] = low[node].min(discovery[neighbor]);
        }
        // RDL❗✔️:     }
    }
    // RDL❗✔️: }
    // END RDL C FUNCTION RDL_tarjanVisit
}

#[cfg(test)]
mod tests {
    use super::{
        BiconnectedComponents, CycleFamilies, CycleFamily, EdgeId, Graph, RingDecomposerError,
        RingDecomposition, ShortestPathInfo, UrfInfo, count_paths,
        count_relevant_cycles_for_family,
    };

    #[test]
    fn graph_add_undirected_edge_matches_rdl_validation() {
        let mut graph = Graph::new(3);

        assert_eq!(graph.add_undirected_edge(0, 1).unwrap().index(), 0);
        assert_eq!(graph.add_undirected_edge(2, 1).unwrap().index(), 1);
        assert!(graph.is_adjacent(0, 1));
        assert!(graph.is_adjacent(1, 0));
        assert!(graph.is_adjacent(1, 2));
        assert_eq!(graph.edge_id(1, 0).unwrap().index(), 0);
        assert_eq!(graph.edge_id(1, 2).unwrap().index(), 1);

        assert!(matches!(
            graph.add_undirected_edge(0, 3),
            Err(RingDecomposerError::NodeOutOfRange {
                node: 3,
                node_count: 3
            })
        ));
        assert!(matches!(
            graph.add_undirected_edge(1, 1),
            Err(RingDecomposerError::SelfLoop { node: 1 })
        ));
        assert!(matches!(
            graph.add_undirected_edge(1, 0),
            Err(RingDecomposerError::DuplicateEdge { from: 1, to: 0 })
        ));
    }

    #[test]
    fn graph_connectivity_matches_rdl_depth_first_logic() {
        let mut graph = Graph::new(4);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();

        assert!(!graph.is_connected());

        graph.add_undirected_edge(2, 3).unwrap();

        assert!(graph.is_connected());
    }

    #[test]
    fn decomposition_handles_acyclic_graph_without_guessing_urf_state() {
        let mut graph = Graph::new(3);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();

        let decomposition = RingDecomposition::calculate(graph).unwrap();

        assert_eq!(decomposition.urf_count(), 0);
        assert_eq!(decomposition.relevant_cycle_count(), 0.0);
        assert!(decomposition.urfs().is_empty());
    }

    #[test]
    fn tarjan_bcc_drops_trivial_edges_like_rdl() {
        let mut graph = Graph::new(4);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 3).unwrap();

        let bcc = BiconnectedComponents::calculate(&graph);

        assert_eq!(bcc.component_count(), 0);
        assert_eq!(bcc.edge_component(EdgeId::new(0)), None);
        assert_eq!(bcc.edge_component(EdgeId::new(1)), None);
        assert_eq!(bcc.edge_component(EdgeId::new(2)), None);
    }

    #[test]
    fn tarjan_bcc_keeps_simple_cycle_as_one_nontrivial_component() {
        let mut graph = Graph::new(3);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 0).unwrap();

        let bcc = BiconnectedComponents::calculate(&graph);

        assert_eq!(bcc.component_count(), 1);
        let component = &bcc.components()[0];
        assert_eq!(component.original_nodes(), &[0, 1, 2]);
        assert_eq!(
            component
                .original_edges()
                .iter()
                .map(|edge| edge.index())
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
        assert_eq!(component.graph().node_count(), 3);
        assert_eq!(component.graph().edge_count(), 3);
        assert_eq!(bcc.edge_component(EdgeId::new(0)), Some((0, 0)));
        assert_eq!(bcc.edge_component(EdgeId::new(1)), Some((0, 1)));
        assert_eq!(bcc.edge_component(EdgeId::new(2)), Some((0, 2)));
        assert_eq!(bcc.node_components(0).unwrap(), &[(0, 0)]);
        assert_eq!(bcc.node_components(1).unwrap(), &[(0, 1)]);
        assert_eq!(bcc.node_components(2).unwrap(), &[(0, 2)]);
    }

    #[test]
    fn tarjan_bcc_splits_two_cycles_sharing_articulation() {
        let mut graph = Graph::new(5);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 0).unwrap();
        graph.add_undirected_edge(2, 3).unwrap();
        graph.add_undirected_edge(3, 4).unwrap();
        graph.add_undirected_edge(4, 2).unwrap();

        let bcc = BiconnectedComponents::calculate(&graph);

        assert_eq!(bcc.component_count(), 2);
        assert_eq!(
            bcc.components()
                .iter()
                .map(|component| component.original_edges().len())
                .collect::<Vec<_>>(),
            vec![3, 3]
        );
        assert_eq!(bcc.node_components(2).unwrap().len(), 2);
        assert_eq!(bcc.components()[0].graph().edge_count(), 3);
        assert_eq!(bcc.components()[1].graph().edge_count(), 3);
    }

    #[test]
    fn apsp_records_shortest_distances_and_predecessors() {
        let mut graph = Graph::new(4);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 3).unwrap();

        let paths = ShortestPathInfo::calculate(&graph);

        assert_eq!(paths.distance(0, 0), Some(0));
        assert_eq!(paths.distance(0, 3), Some(3));
        assert_eq!(paths.predecessor(0, 0), Some(0));
        assert_eq!(paths.predecessor(0, 3), Some(2));
        assert_eq!(paths.predecessor(3, 0), Some(1));
        assert_eq!(paths.distance(0, 99), None);
    }

    #[test]
    fn apsp_reachable_preceding_matches_degree_ordering_rule() {
        let mut graph = Graph::new(3);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 0).unwrap();

        let paths = ShortestPathInfo::calculate(&graph);

        assert!(!paths.reachable_preceding(0, 1));
        assert!(!paths.reachable_preceding(0, 2));
        assert!(paths.reachable_preceding(1, 0));
        assert!(!paths.reachable_preceding(1, 2));
        assert!(paths.reachable_preceding(2, 0));
        assert!(paths.reachable_preceding(2, 1));
        assert_eq!(paths.directed_path_graphs().len(), 3);
    }

    #[test]
    fn directed_path_graph_rejects_invalid_arcs_and_preserves_row_order() {
        let mut graph = super::DirectedPathGraph::new(3);
        graph.add_directed_edge(0, 2);
        graph.add_directed_edge(0, 2);
        graph.add_directed_edge(0, 0);
        graph.add_directed_edge(3, 0);
        graph.add_directed_edge(0, 3);
        graph.add_directed_edge(0, 1);
        graph.add_directed_edge(2, 1);

        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.neighbors(0), Some(&[2, 1][..]));
        assert_eq!(graph.neighbors(1), Some(&[][..]));
        assert_eq!(graph.neighbors(2), Some(&[1][..]));
        assert_eq!(graph.neighbors(3), None);
    }

    #[test]
    fn apsp_preserves_equal_path_predecessor_and_unreachable_rows() {
        let mut graph = Graph::new(5);
        graph.add_undirected_edge(0, 2).unwrap();
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(2, 3).unwrap();
        graph.add_undirected_edge(1, 3).unwrap();

        let paths = ShortestPathInfo::calculate(&graph);

        assert_eq!(paths.distance(0, 0), Some(0));
        assert_eq!(paths.predecessor(0, 0), Some(0));
        assert_eq!(paths.distance(0, 3), Some(2));
        assert_eq!(paths.predecessor(0, 3), Some(2));
        assert_eq!(paths.distance(3, 0), Some(2));
        assert_eq!(paths.predecessor(3, 0), Some(2));
        for connected in 0..4 {
            assert_eq!(paths.distance(connected, 4), None);
            assert_eq!(paths.predecessor(connected, 4), None);
            assert!(!paths.reachable_preceding(connected, 4));
            assert_eq!(paths.distance(4, connected), None);
            assert_eq!(paths.predecessor(4, connected), None);
            assert!(!paths.reachable_preceding(4, connected));
        }
        assert_eq!(paths.distance(4, 4), Some(0));
        assert_eq!(paths.predecessor(4, 4), Some(4));
        assert_eq!(paths.distance(99, 0), None);
        assert_eq!(paths.predecessor(0, 99), None);
        assert!(!paths.reachable_preceding(99, 0));
    }

    #[test]
    fn tarjan_bcc_maps_every_row_for_disconnected_cycles() {
        let mut graph = Graph::new(7);
        for (from, to) in [(0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3)] {
            graph.add_undirected_edge(from, to).unwrap();
        }

        let bcc = BiconnectedComponents::calculate(&graph);

        assert_eq!(bcc.component_count(), 2);
        assert_eq!(bcc.components()[0].original_nodes(), &[0, 1, 2]);
        assert_eq!(
            bcc.components()[0].original_edges(),
            &[EdgeId::new(0), EdgeId::new(1), EdgeId::new(2)]
        );
        assert_eq!(bcc.components()[1].original_nodes(), &[3, 4, 5]);
        assert_eq!(
            bcc.components()[1].original_edges(),
            &[EdgeId::new(3), EdgeId::new(4), EdgeId::new(5)]
        );
        for edge in 0..3 {
            assert_eq!(bcc.edge_component(EdgeId::new(edge)), Some((0, edge)));
        }
        for edge in 3..6 {
            assert_eq!(bcc.edge_component(EdgeId::new(edge)), Some((1, edge - 3)));
        }
        assert_eq!(bcc.edge_component(EdgeId::new(6)), None);
        assert_eq!(bcc.edge_component(EdgeId::new(99)), None);
        assert_eq!(bcc.node_components(0), Some(&[(0, 0)][..]));
        assert_eq!(bcc.node_components(1), Some(&[(0, 1)][..]));
        assert_eq!(bcc.node_components(2), Some(&[(0, 2)][..]));
        assert_eq!(bcc.node_components(3), Some(&[(1, 0)][..]));
        assert_eq!(bcc.node_components(4), Some(&[(1, 1)][..]));
        assert_eq!(bcc.node_components(5), Some(&[(1, 2)][..]));
        assert_eq!(bcc.node_components(6), Some(&[][..]));
        assert_eq!(bcc.node_components(7), None);
    }

    #[test]
    fn tarjan_bcc_maps_articulation_into_both_local_index_spaces() {
        let mut graph = Graph::new(5);
        for (from, to) in [(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 2)] {
            graph.add_undirected_edge(from, to).unwrap();
        }

        let bcc = BiconnectedComponents::calculate(&graph);

        assert_eq!(bcc.component_count(), 2);
        assert_eq!(bcc.components()[0].original_nodes(), &[2, 3, 4]);
        assert_eq!(
            bcc.components()[0].original_edges(),
            &[EdgeId::new(3), EdgeId::new(4), EdgeId::new(5)]
        );
        assert_eq!(bcc.components()[1].original_nodes(), &[0, 1, 2]);
        assert_eq!(
            bcc.components()[1].original_edges(),
            &[EdgeId::new(0), EdgeId::new(1), EdgeId::new(2)]
        );
        assert_eq!(bcc.edge_component(EdgeId::new(0)), Some((1, 0)));
        assert_eq!(bcc.edge_component(EdgeId::new(1)), Some((1, 1)));
        assert_eq!(bcc.edge_component(EdgeId::new(2)), Some((1, 2)));
        assert_eq!(bcc.edge_component(EdgeId::new(3)), Some((0, 0)));
        assert_eq!(bcc.edge_component(EdgeId::new(4)), Some((0, 1)));
        assert_eq!(bcc.edge_component(EdgeId::new(5)), Some((0, 2)));
        assert_eq!(bcc.node_components(0), Some(&[(1, 0)][..]));
        assert_eq!(bcc.node_components(1), Some(&[(1, 1)][..]));
        assert_eq!(bcc.node_components(2), Some(&[(1, 2), (0, 0)][..]));
        assert_eq!(bcc.node_components(3), Some(&[(0, 1)][..]));
        assert_eq!(bcc.node_components(4), Some(&[(0, 2)][..]));
        for component in bcc.components() {
            assert_eq!(component.graph().node_count(), 3);
            assert_eq!(component.graph().edge_count(), 3);
        }
    }

    #[test]
    fn cycle_families_find_triangle_family() {
        let mut graph = Graph::new(3);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 0).unwrap();
        let mut paths = ShortestPathInfo::calculate(&graph);

        let cycle_families = CycleFamilies::calculate(&mut graph, &mut paths);

        assert_eq!(cycle_families.len(), 1);
        let family = &cycle_families.families()[0];
        assert_eq!(family.weight(), 3);
        assert_eq!(family.x(), None);
        assert_eq!(
            family
                .prototype()
                .iter()
                .filter(|&&contains| contains)
                .count(),
            3
        );
    }

    #[test]
    fn decomposition_returns_single_urf_for_triangle() {
        let mut graph = Graph::new(3);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 0).unwrap();

        let decomposition = RingDecomposition::calculate(graph).unwrap();

        assert_eq!(decomposition.urf_count(), 1);
        assert_eq!(decomposition.relevant_cycle_count(), 1.0);
        assert_eq!(decomposition.urfs()[0].nodes(), &[0, 1, 2]);
        assert_eq!(
            decomposition.urfs()[0]
                .edges()
                .iter()
                .map(|edge| edge.index())
                .collect::<Vec<_>>(),
            vec![0, 1, 2]
        );
    }

    #[test]
    fn decomposition_returns_single_urf_for_square() {
        let mut graph = Graph::new(4);
        graph.add_undirected_edge(0, 1).unwrap();
        graph.add_undirected_edge(1, 2).unwrap();
        graph.add_undirected_edge(2, 3).unwrap();
        graph.add_undirected_edge(3, 0).unwrap();

        let decomposition = RingDecomposition::calculate(graph).unwrap();

        assert_eq!(decomposition.urf_count(), 1);
        assert_eq!(decomposition.relevant_cycle_count(), 1.0);
        assert_eq!(decomposition.urfs()[0].nodes(), &[0, 1, 2, 3]);
        assert_eq!(
            decomposition.urfs()[0]
                .edges()
                .iter()
                .map(|edge| edge.index())
                .collect::<Vec<_>>(),
            vec![0, 1, 2, 3]
        );
    }

    fn test_family(weight: usize, prototype: &[bool]) -> CycleFamily {
        CycleFamily {
            weight,
            r: 0,
            p: 0,
            q: 0,
            x: None,
            prototype: prototype.to_vec(),
            relevant: false,
        }
    }

    #[test]
    fn urf_fast_path_and_empty_initialization_keep_exact_private_state() {
        let empty = CycleFamilies {
            families: Vec::new(),
        };
        assert_eq!(
            UrfInfo::init_urf_info(&empty),
            UrfInfo {
                nof_protos: Vec::new(),
                relations: Vec::new(),
                urfs: Vec::new(),
            }
        );

        for family_count in 1..=2 {
            let mut families = CycleFamilies {
                families: (0..family_count)
                    .map(|index| {
                        let mut prototype = vec![false; 2];
                        prototype[index] = true;
                        test_family(3, &prototype)
                    })
                    .collect(),
            };
            let mut info = UrfInfo::init_urf_info(&families);
            info.check_dependencies(&mut families, &Graph::new(2));
            info.fill_urfs();
            assert!(families.families().iter().all(CycleFamily::is_relevant));
            assert_eq!(
                info.urfs,
                (0..family_count)
                    .map(|index| vec![index])
                    .collect::<Vec<_>>()
            );
        }
    }

    #[test]
    fn transitive_closure_relates_complete_components_in_source_order() {
        let mut info = UrfInfo {
            nof_protos: vec![4],
            relations: vec![vec![vec![false; 4]; 4]],
            urfs: Vec::new(),
        };
        for index in 0..4 {
            info.relations[0][index][index] = true;
        }
        info.relations[0][0][1] = true;
        info.relations[0][1][0] = true;
        info.relations[0][1][2] = true;
        info.relations[0][2][1] = true;

        info.find_transitive_closure();
        assert!(info.relations[0][0][2]);
        assert!(info.relations[0][2][0]);
        assert!(!info.relations[0][0][3]);
        info.fill_urfs();
        assert_eq!(info.urfs, vec![vec![0, 1, 2], vec![3]]);
    }

    #[test]
    fn gaussian_dependency_distinguishes_smaller_equal_and_new_pivots() {
        let mut smaller_dependent = CycleFamilies {
            families: vec![
                test_family(3, &[true, true, false, false]),
                test_family(3, &[true, false, true, false]),
                test_family(4, &[false, true, true, false]),
            ],
        };
        let mut smaller_info = UrfInfo::init_urf_info(&smaller_dependent);
        smaller_info.check_dependencies(&mut smaller_dependent, &Graph::new(4));
        assert_eq!(
            smaller_dependent
                .families()
                .iter()
                .map(CycleFamily::is_relevant)
                .collect::<Vec<_>>(),
            vec![true, true, false]
        );

        let mut equal_dependent = CycleFamilies {
            families: vec![
                test_family(3, &[true, true, false, false]),
                test_family(3, &[true, false, true, false]),
                test_family(3, &[false, true, true, false]),
            ],
        };
        let mut equal_info = UrfInfo::init_urf_info(&equal_dependent);
        equal_info.check_dependencies(&mut equal_dependent, &Graph::new(4));
        assert!(
            equal_dependent
                .families()
                .iter()
                .all(CycleFamily::is_relevant)
        );

        let mut swapped_pivots = CycleFamilies {
            families: vec![
                test_family(3, &[false, true, false, false]),
                test_family(3, &[false, false, true, false]),
                test_family(3, &[false, false, false, true]),
            ],
        };
        let mut swapped_info = UrfInfo::init_urf_info(&swapped_pivots);
        swapped_info.check_dependencies(&mut swapped_pivots, &Graph::new(4));
        assert!(
            swapped_pivots
                .families()
                .iter()
                .all(CycleFamily::is_relevant)
        );
    }

    #[test]
    fn complete_shortest_path_dag_counts_two_paths_and_rejects_a_cycle() {
        let graph = Graph::new(4);
        let mut paths = ShortestPathInfo::initialized_for_graph(&graph);
        paths.directed_paths[0].add_directed_edge(3, 1);
        paths.directed_paths[0].add_directed_edge(3, 2);
        paths.directed_paths[0].add_directed_edge(1, 0);
        paths.directed_paths[0].add_directed_edge(2, 0);
        assert_eq!(count_paths(0, 3, 4, &paths), Ok(2.0));

        let mut invalid = ShortestPathInfo::initialized_for_graph(&graph);
        invalid.directed_paths[0].add_directed_edge(3, 1);
        invalid.directed_paths[0].add_directed_edge(1, 3);
        assert_eq!(
            count_paths(0, 3, 4, &invalid),
            Err(RingDecomposerError::InvalidShortestPathDag { root: 0, target: 3 })
        );
    }

    #[test]
    fn relevant_cycle_product_overflow_is_a_structured_error() {
        let layer_count = 513usize;
        let node_count = layer_count * 2 + 2;
        let target = node_count - 1;
        let graph = Graph::new(node_count);
        let mut paths = ShortestPathInfo::initialized_for_graph(&graph);
        let mut previous = vec![0usize];
        for layer in 0..layer_count {
            let current = vec![layer * 2 + 1, layer * 2 + 2];
            for &node in &current {
                for &predecessor in &previous {
                    paths.directed_paths[0].add_directed_edge(node, predecessor);
                }
            }
            previous = current;
        }
        for predecessor in previous {
            paths.directed_paths[0].add_directed_edge(target, predecessor);
        }
        let family = CycleFamily {
            weight: 3,
            r: 0,
            p: target,
            q: target,
            x: None,
            prototype: Vec::new(),
            relevant: true,
        };

        assert_eq!(
            count_relevant_cycles_for_family(&family, &graph, &paths),
            Err(RingDecomposerError::RelevantCycleCountOverflow)
        );
    }

    #[test]
    fn private_relation_pipeline_covers_smaller_and_equal_weight_cycles() {
        let mut fused = Graph::new(4);
        for (from, to) in [(0, 1), (1, 2), (2, 0), (1, 3), (3, 2)] {
            fused.add_undirected_edge(from, to).unwrap();
        }
        let mut fused_paths = ShortestPathInfo::calculate(&fused);
        let mut fused_families = CycleFamilies::calculate(&mut fused, &mut fused_paths);
        assert_eq!(
            fused_families
                .families()
                .iter()
                .map(CycleFamily::weight)
                .collect::<Vec<_>>(),
            vec![3, 3]
        );
        let fused_info = UrfInfo::check_urf_relation(&mut fused_families, &fused, &fused_paths);
        assert_eq!(
            fused_families
                .families()
                .iter()
                .map(CycleFamily::is_relevant)
                .collect::<Vec<_>>(),
            vec![true, true]
        );
        assert_eq!(fused_info.urfs, vec![vec![0], vec![1]]);

        let mut theta = Graph::new(5);
        for (from, to) in [(0, 1), (1, 4), (0, 2), (2, 4), (0, 3), (3, 4)] {
            theta.add_undirected_edge(from, to).unwrap();
        }
        let mut theta_paths = ShortestPathInfo::calculate(&theta);
        let mut theta_families = CycleFamilies::calculate(&mut theta, &mut theta_paths);
        assert_eq!(theta_families.len(), 3);
        assert!(
            theta_families
                .families()
                .iter()
                .all(|family| family.weight() == 4)
        );
        let theta_info = UrfInfo::check_urf_relation(&mut theta_families, &theta, &theta_paths);
        assert!(
            theta_families
                .families()
                .iter()
                .all(CycleFamily::is_relevant)
        );
        assert_eq!(theta_info.urfs, vec![vec![0], vec![1], vec![2]]);
    }
}
