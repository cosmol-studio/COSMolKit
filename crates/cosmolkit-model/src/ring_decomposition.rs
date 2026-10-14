//! Complete detached RDL state in CK-owned vocabulary, without an engine
//! dependency. CORE supplies source-produced state; MODEL does not find rings.

use std::collections::BTreeMap;

#[doc(hidden)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RingEdge {
    pub id: usize,
    pub from: usize,
    pub to: usize,
}

impl RingEdge {
    pub const fn from(&self) -> usize {
        self.from
    }
    pub const fn to(&self) -> usize {
        self.to
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingGraph {
    pub node_count: usize,
    pub adjacency: Vec<Vec<(usize, usize)>>,
    pub edges: Vec<RingEdge>,
    pub edge_ids_by_endpoints: BTreeMap<(usize, usize), usize>,
}

impl RingGraph {
    pub const fn node_count(&self) -> usize {
        self.node_count
    }
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }
    pub fn edges(&self) -> &[RingEdge] {
        &self.edges
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingCycleFamily {
    pub weight: usize,
    pub root: usize,
    pub p: usize,
    pub q: usize,
    pub x: Option<usize>,
    pub prototype: Vec<bool>,
    pub relevant: bool,
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingComponent {
    pub graph: RingGraph,
    pub original_nodes: Vec<usize>,
    pub original_edges: Vec<usize>,
    pub predecessor: Vec<Vec<Option<usize>>>,
    pub distance: Vec<Vec<Option<usize>>>,
    pub reachable: Vec<Vec<bool>>,
    pub directed_paths: Vec<Vec<Vec<usize>>>,
    pub cycle_families: Vec<RingCycleFamily>,
    pub prototype_counts: Vec<usize>,
    pub relations: Vec<Vec<Vec<bool>>>,
    pub urfs: Vec<Vec<usize>>,
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RingFamily {
    pub nodes: Vec<usize>,
    pub edges: Vec<usize>,
}

impl RingFamily {
    pub fn nodes(&self) -> &[usize] {
        &self.nodes
    }
    pub fn edges(&self) -> &[usize] {
        &self.edges
    }
}

#[doc(hidden)]
#[derive(Debug, Clone, PartialEq)]
pub struct RingDecomposition {
    pub graph: RingGraph,
    pub urfs: Vec<RingFamily>,
    relevant_cycle_count: f64,
    pub components: Vec<RingComponent>,
}

// As in the source engine, construction excludes NaN and leaves structural
// f64 equality unchanged (including +0 == -0). Never use Arc identity for Eq.
impl Eq for RingDecomposition {}

impl RingDecomposition {
    pub fn new(
        graph: RingGraph,
        urfs: Vec<RingFamily>,
        relevant_cycle_count: f64,
        components: Vec<RingComponent>,
    ) -> Result<Self, &'static str> {
        if !relevant_cycle_count.is_finite()
            || relevant_cycle_count < 0.0
            || relevant_cycle_count.fract() != 0.0
        {
            return Err("invalid RDL relevant-cycle count");
        }
        Ok(Self {
            graph,
            urfs,
            relevant_cycle_count,
            components,
        })
    }
    pub const fn graph(&self) -> &RingGraph {
        &self.graph
    }
    pub fn urfs(&self) -> &[RingFamily] {
        &self.urfs
    }
    pub fn urf_count(&self) -> usize {
        self.urfs.len()
    }
    pub const fn relevant_cycle_count(&self) -> f64 {
        self.relevant_cycle_count
    }
}
