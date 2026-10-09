//! Source-backed graph paths and detached topology subsets.
//!
//! This module owns graph algorithms over validated detached model values. It
//! never accepts a live molecule or runtime capability.

use std::collections::{BTreeMap, VecDeque};

use cosmolkit_model::{
    AtomId, AtomMapping, AtomQueryPredicate, BondId, BondMapping, BondQueryPredicate, BondStereo,
    MappingValidationError, QueryAtom, QueryBond, QueryGraph, QueryGraphError, QueryNode,
    StereoGroup, SubstanceGroup, SubstanceGroupId, TemplateAttachmentOrderError, TopologyBlock,
    TopologyMapping, TopologyValidationError,
};

use crate::{PeriodicTableError, atomic_mass};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PathRepresentation {
    Bonds,
    Atoms,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GraphPath {
    Bonds(Vec<BondId>),
    Atoms(Vec<AtomId>),
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PathSearchParams<'a> {
    /// Borrowed source ignoreAtoms bitmap; None differs from Some(empty).
    pub ignore_atoms: Option<&'a [bool]>,
    pub representation: PathRepresentation,
    pub use_hydrogens: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub only_shortest_paths: bool,
}

impl Default for PathSearchParams<'_> {
    fn default() -> Self {
        Self {
            ignore_atoms: None,
            representation: PathRepresentation::Bonds,
            use_hydrogens: false,
            rooted_at_atom: None,
            only_shortest_paths: false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SubgraphSearchParams<'a> {
    /// Borrowed source ignoreAtoms bitmap; None differs from Some(empty).
    pub ignore_atoms: Option<&'a [bool]>,
    pub use_hydrogens: bool,
    pub rooted_at_atom: Option<AtomId>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UniqueSubgraphParams<'a> {
    /// Borrowed source ignoreAtoms bitmap; None differs from Some(empty).
    pub ignore_atoms: Option<&'a [bool]>,
    pub use_hydrogens: bool,
    pub use_bond_orders: bool,
    pub rooted_at_atom: Option<AtomId>,
    pub extra_atom_invariants: Option<Vec<u32>>,
}

impl Default for UniqueSubgraphParams<'_> {
    fn default() -> Self {
        Self {
            use_hydrogens: false,
            ignore_atoms: None,
            use_bond_orders: true,
            rooted_at_atom: None,
            extra_atom_invariants: None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AtomEnvironmentParams {
    pub use_hydrogens: bool,
    pub enforce_radius: bool,
}

impl Default for AtomEnvironmentParams {
    fn default() -> Self {
        Self {
            use_hydrogens: false,
            enforce_radius: true,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConnectedComponents {
    pub atom_to_component: Vec<usize>,
    pub components: Vec<Vec<AtomId>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomEnvironment {
    pub bonds: Vec<BondId>,
    pub atom_distances: Vec<Option<usize>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct SubtopologyParams {
    pub copy_as_query: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum DetachedPathSubgraph {
    Concrete(TopologyBlock),
    Query {
        graph: QueryGraph,
        substance_groups: Vec<SubstanceGroup>,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SubtopologyResult {
    pub subgraph: DetachedPathSubgraph,
    pub mapping: TopologyMapping,
}

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum PathError {
    #[error("{0}")]
    StereoGroup(#[from] cosmolkit_model::StereoGroupError),

    #[error("molecule property operation failed: {0}")]
    MoleculeProperty(#[from] cosmolkit_model::MoleculePropertyError),
    #[error("bond property operation failed: {0}")]
    BondProperty(#[from] cosmolkit_model::BondValueError),
    #[error("atom property operation failed: {0}")]
    AtomProperty(#[from] cosmolkit_model::AtomPropertyError),
    #[error(transparent)]
    InvalidTopology(TopologyValidationError),
    #[error("{role} atom {atom} is out of range for {atom_count} atoms")]
    AtomOutOfRange {
        role: &'static str,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("shortest-path endpoints must be distinct (both were atom {atom})")]
    EqualShortestPathEndpoints { atom: AtomId },
    #[error("path length range {lower}..={upper} is invalid")]
    InvalidLengthRange { lower: usize, upper: usize },
    #[error("path length {length} cannot be represented by the source unsigned range")]
    LengthOverflow { length: usize },
    #[error("atom path row {position} refers to atom {atom}, outside {atom_count} atoms")]
    AtomPathOutOfRange {
        position: usize,
        atom: AtomId,
        atom_count: usize,
    },
    #[error("bad ignoreAtoms size: {actual} rows, expected {expected}")]
    IgnoredAtomMaskLength { actual: usize, expected: usize },
    #[error("extra atom invariants have {actual} rows, expected {expected}")]
    ExtraInvariantLength { actual: usize, expected: usize },
    #[error(transparent)]
    PeriodicTable(PeriodicTableError),
    #[error("detached subset topology is invalid: {0}")]
    InvalidSubsetTopology(TopologyValidationError),
    #[error(transparent)]
    QueryGraph(QueryGraphError),
    #[error(transparent)]
    Mapping(MappingValidationError),
    #[error("atom {carrier} template attachment remap failed: {source}")]
    TemplateAttachmentRemap {
        carrier: AtomId,
        source: TemplateAttachmentOrderError,
    },
    #[error(transparent)]
    Matrix(crate::MatrixError),
}

/// Borrowed connectivity only. Query atomic number 0 remains 0. All loops
/// below are shared; access performs O(1) indexed reads and allocates nothing.
#[derive(Clone, Copy)]
enum PathGraphAccess<'a> {
    Concrete(&'a TopologyBlock),
    Query(&'a QueryGraph),
}
impl PathGraphAccess<'_> {
    fn atom_count(self) -> usize {
        match self {
            Self::Concrete(t) => t.atoms.len(),
            Self::Query(q) => q.num_atoms(),
        }
    }
    fn bond_count(self) -> usize {
        match self {
            Self::Concrete(t) => t.bonds.len(),
            Self::Query(q) => q.num_bonds(),
        }
    }
    fn atomic_number(self, atom: usize) -> u8 {
        match self {
            Self::Concrete(t) => t.atoms[atom].atomic_number(),
            Self::Query(q) => q.atoms()[atom].atomic_number(),
        }
    }
    fn bond(&self, index: usize) -> &cosmolkit_model::Bond {
        match self {
            Self::Concrete(t) => &t.bonds[index],
            Self::Query(q) => q.bonds()[index].bond(),
        }
    }
    fn neighbor_count(self, atom: usize) -> usize {
        match self {
            Self::Concrete(t) => t.adjacency.neighbors_of(atom).len(),
            Self::Query(q) => q.adjacency()[atom].len(),
        }
    }
    fn neighbor(self, atom: usize, position: usize) -> (usize, BondId) {
        match self {
            Self::Concrete(t) => {
                let n = t.adjacency.neighbors_of(atom)[position];
                (n.atom_index, n.bond)
            }
            Self::Query(q) => {
                let (other, bond) = q.adjacency()[atom][position];
                (other, BondId::new(bond))
            }
        }
    }
}

pub(crate) trait NeighborSource {
    fn atom_count(&self) -> usize;
    fn neighbor_count(&self, atom: usize) -> usize;
    fn neighbor_at(&self, atom: usize, position: usize) -> usize;
}

impl NeighborSource for TopologyBlock {
    fn atom_count(&self) -> usize {
        self.atoms.len()
    }

    fn neighbor_count(&self, atom: usize) -> usize {
        self.adjacency.neighbors_of(atom).len()
    }

    fn neighbor_at(&self, atom: usize, position: usize) -> usize {
        self.adjacency.neighbors_of(atom)[position].atom_index
    }
}

impl NeighborSource for QueryGraph {
    fn atom_count(&self) -> usize {
        self.num_atoms()
    }

    fn neighbor_count(&self, atom: usize) -> usize {
        self.adjacency()[atom].len()
    }

    fn neighbor_at(&self, atom: usize, position: usize) -> usize {
        self.adjacency()[atom][position].0
    }
}

pub(crate) fn connected_components_from_source(
    source: &impl NeighborSource,
) -> ConnectedComponents {
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags
    // RDKit✔️❌: unsigned int getMolFrags(const ROMol &mol, INT_VECT &mapping) {
    // RDKit✔️❌:   unsigned int natms = mol.getNumAtoms();
    // RDKit✔️❌:   mapping.resize(natms);
    // RDKit✔️❌:   return natms ? boost::connected_components(mol.getTopology(), &mapping[0])
    // RDKit✔️❌:                : 0;
    // RDKit✔️❌: };
    // END RDKIT CPP FUNCTION MolOps::getMolFrags(INT_VECT&)
    // The canonical paired result also owns atom groups required by the
    // companion overload. That output adds allocations compared with the
    // mapping-only source overload; no allocation-equivalence claim is made.
    // BEGIN REACHED BOOST connected_components / recorder / depth_first_search
    // Source defining body connected_components.hpp::components_recorder
    // Boost✔️✔️:     template < class ComponentsMap >
    // Boost✔️✔️:     class components_recorder : public dfs_visitor<>
    // Boost✔️✔️:     {
    // Boost✔️✔️:         typedef typename property_traits< ComponentsMap >::value_type comp_type;
    // Boost✔️✔️:
    // Boost✔️✔️:     public:
    // Boost✔️✔️:         components_recorder(ComponentsMap c, comp_type& c_count)
    // Boost✔️✔️:         : m_component(c), m_count(c_count)
    // Boost✔️✔️:         {
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:         template < class Vertex, class Graph > void start_vertex(Vertex, Graph&)
    // Boost✔️✔️:         {
    // Boost✔️✔️:             if (m_count == (std::numeric_limits< comp_type >::max)())
    // Boost✔️✔️:                 m_count = 0; // start counting components at zero
    // Boost✔️✔️:             else
    // Boost✔️✔️:                 ++m_count;
    // Boost✔️✔️:         }
    // Boost✔️✔️:         template < class Vertex, class Graph >
    // Boost✔️✔️:         void discover_vertex(Vertex u, Graph&)
    // Boost✔️✔️:         {
    // Boost✔️✔️:             put(m_component, u, m_count);
    // Boost✔️✔️:         }
    // Boost✔️✔️:
    // Boost✔️✔️:     protected:
    // Boost✔️✔️:         ComponentsMap m_component;
    // Boost✔️✔️:         comp_type& m_count;
    // Boost✔️✔️:     };
    // Source defining body connected_components.hpp::connected_components
    // Boost✔️✔️: template < class Graph, class ComponentMap >
    // Boost✔️✔️: inline typename property_traits< ComponentMap >::value_type
    // Boost✔️✔️: connected_components(const Graph& g,
    // Boost✔️✔️:     ComponentMap c BOOST_GRAPH_ENABLE_IF_MODELS_PARM(
    // Boost✔️✔️:         Graph, vertex_list_graph_tag))
    // Boost✔️✔️: {
    // Boost✔️✔️:     if (num_vertices(g) == 0)
    // Boost✔️✔️:         return 0;
    // Boost✔️✔️:
    // Boost✔️✔️:     typedef typename graph_traits< Graph >::vertex_descriptor Vertex;
    // Boost✔️✔️:     BOOST_CONCEPT_ASSERT((WritablePropertyMapConcept< ComponentMap, Vertex >));
    // Boost✔️✔️:     // typedef typename boost::graph_traits<Graph>::directed_category directed;
    // Boost✔️✔️:     // BOOST_STATIC_ASSERT((boost::is_same<directed, undirected_tag>::value));
    // Boost✔️✔️:
    // Boost✔️✔️:     typedef typename property_traits< ComponentMap >::value_type comp_type;
    // Boost✔️✔️:     // c_count initialized to "nil" (with nil represented by (max)())
    // Boost✔️✔️:     comp_type c_count((std::numeric_limits< comp_type >::max)());
    // Boost✔️✔️:     detail::components_recorder< ComponentMap > vis(c, c_count);
    // Boost✔️✔️:     depth_first_search(g, visitor(vis));
    // Boost✔️✔️:     return c_count + 1;
    // Boost✔️✔️: }
    // Source defining body depth_first_search.hpp::depth_first_visit_impl
    // Boost✔️✔️:     template < class IncidenceGraph, class DFSVisitor, class ColorMap,
    // Boost✔️✔️:         class TerminatorFunc >
    // Boost✔️✔️:     void depth_first_visit_impl(const IncidenceGraph& g,
    // Boost✔️✔️:         typename graph_traits< IncidenceGraph >::vertex_descriptor u,
    // Boost✔️✔️:         DFSVisitor& vis, ColorMap color, TerminatorFunc func = TerminatorFunc())
    // Boost✔️✔️:     {
    // Boost✔️✔️:         BOOST_CONCEPT_ASSERT((IncidenceGraphConcept< IncidenceGraph >));
    // Boost✔️✔️:         BOOST_CONCEPT_ASSERT((DFSVisitorConcept< DFSVisitor, IncidenceGraph >));
    // Boost✔️✔️:         typedef
    // Boost✔️✔️:             typename graph_traits< IncidenceGraph >::vertex_descriptor Vertex;
    // Boost✔️✔️:         typedef typename graph_traits< IncidenceGraph >::edge_descriptor Edge;
    // Boost✔️✔️:         BOOST_CONCEPT_ASSERT((ReadWritePropertyMapConcept< ColorMap, Vertex >));
    // Boost✔️✔️:         typedef typename property_traits< ColorMap >::value_type ColorValue;
    // Boost✔️✔️:         BOOST_CONCEPT_ASSERT((ColorValueConcept< ColorValue >));
    // Boost✔️✔️:         typedef color_traits< ColorValue > Color;
    // Boost✔️✔️:         typedef typename graph_traits< IncidenceGraph >::out_edge_iterator Iter;
    // Boost✔️✔️:         typedef std::pair< Vertex,
    // Boost✔️✔️:             std::pair< boost::optional< Edge >, std::pair< Iter, Iter > > >
    // Boost✔️✔️:             VertexInfo;
    // Boost✔️✔️:
    // Boost✔️✔️:         boost::optional< Edge > src_e;
    // Boost✔️✔️:         Iter ei, ei_end;
    // Boost✔️✔️:         std::vector< VertexInfo > stack;
    // Boost✔️✔️:
    // Boost✔️✔️:         // Possible optimization for vector
    // Boost✔️✔️:         // stack.reserve(num_vertices(g));
    // Boost✔️✔️:
    // Boost✔️✔️:         put(color, u, Color::gray());
    // Boost✔️✔️:         vis.discover_vertex(u, g);
    // Boost✔️✔️:         boost::tie(ei, ei_end) = out_edges(u, g);
    // Boost✔️✔️:         if (func(u, g))
    // Boost✔️✔️:         {
    // Boost✔️✔️:             // If this vertex terminates the search, we push empty range
    // Boost✔️✔️:             stack.push_back(std::make_pair(u,
    // Boost✔️✔️:                 std::make_pair(boost::optional< Edge >(),
    // Boost✔️✔️:                     std::make_pair(ei_end, ei_end))));
    // Boost✔️✔️:         }
    // Boost✔️✔️:         else
    // Boost✔️✔️:         {
    // Boost✔️✔️:             stack.push_back(std::make_pair(u,
    // Boost✔️✔️:                 std::make_pair(
    // Boost✔️✔️:                     boost::optional< Edge >(), std::make_pair(ei, ei_end))));
    // Boost✔️✔️:         }
    // Boost✔️✔️:         while (!stack.empty())
    // Boost✔️✔️:         {
    // Boost✔️✔️:             VertexInfo& back = stack.back();
    // Boost✔️✔️:             u = back.first;
    // Boost✔️✔️:             src_e = back.second.first;
    // Boost✔️✔️:             boost::tie(ei, ei_end) = back.second.second;
    // Boost✔️✔️:             stack.pop_back();
    // Boost✔️✔️:             // finish_edge has to be called here, not after the
    // Boost✔️✔️:             // loop. Think of the pop as the return from a recursive call.
    // Boost✔️✔️:             if (src_e)
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 call_finish_edge(vis, src_e.get(), g);
    // Boost✔️✔️:             }
    // Boost✔️✔️:             while (ei != ei_end)
    // Boost✔️✔️:             {
    // Boost✔️✔️:                 Vertex v = target(*ei, g);
    // Boost✔️✔️:                 vis.examine_edge(*ei, g);
    // Boost✔️✔️:                 ColorValue v_color = get(color, v);
    // Boost✔️✔️:                 if (v_color == Color::white())
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     vis.tree_edge(*ei, g);
    // Boost✔️✔️:                     src_e = *ei;
    // Boost✔️✔️:                     stack.push_back(std::make_pair(u,
    // Boost✔️✔️:                         std::make_pair(src_e, std::make_pair(++ei, ei_end))));
    // Boost✔️✔️:                     u = v;
    // Boost✔️✔️:                     put(color, u, Color::gray());
    // Boost✔️✔️:                     vis.discover_vertex(u, g);
    // Boost✔️✔️:                     boost::tie(ei, ei_end) = out_edges(u, g);
    // Boost✔️✔️:                     if (func(u, g))
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         ei = ei_end;
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                 }
    // Boost✔️✔️:                 else
    // Boost✔️✔️:                 {
    // Boost✔️✔️:                     if (v_color == Color::gray())
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         vis.back_edge(*ei, g);
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     else
    // Boost✔️✔️:                     {
    // Boost✔️✔️:                         vis.forward_or_cross_edge(*ei, g);
    // Boost✔️✔️:                     }
    // Boost✔️✔️:                     call_finish_edge(vis, *ei, g);
    // Boost✔️✔️:                     ++ei;
    // Boost✔️✔️:                 }
    // Boost✔️✔️:             }
    // Boost✔️✔️:             put(color, u, Color::black());
    // Boost✔️✔️:             vis.finish_vertex(u, g);
    // Boost✔️✔️:         }
    // Boost✔️✔️:     }
    // Boost✔️✔️:
    // Source defining body depth_first_search.hpp::depth_first_search
    // Boost✔️✔️: template < class VertexListGraph, class DFSVisitor, class ColorMap >
    // Boost✔️✔️: void depth_first_search(const VertexListGraph& g, DFSVisitor vis,
    // Boost✔️✔️:     ColorMap color,
    // Boost✔️✔️:     typename graph_traits< VertexListGraph >::vertex_descriptor start_vertex)
    // Boost✔️✔️: {
    // Boost✔️✔️:     typedef typename graph_traits< VertexListGraph >::vertex_descriptor Vertex;
    // Boost✔️✔️:     BOOST_CONCEPT_ASSERT((DFSVisitorConcept< DFSVisitor, VertexListGraph >));
    // Boost✔️✔️:     typedef typename property_traits< ColorMap >::value_type ColorValue;
    // Boost✔️✔️:     typedef color_traits< ColorValue > Color;
    // Boost✔️✔️:
    // Boost✔️✔️:     typename graph_traits< VertexListGraph >::vertex_iterator ui, ui_end;
    // Boost✔️✔️:     for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui)
    // Boost✔️✔️:     {
    // Boost✔️✔️:         Vertex u = implicit_cast< Vertex >(*ui);
    // Boost✔️✔️:         put(color, u, Color::white());
    // Boost✔️✔️:         vis.initialize_vertex(u, g);
    // Boost✔️✔️:     }
    // Boost✔️✔️:
    // Boost✔️✔️:     if (start_vertex != detail::get_default_starting_vertex(g))
    // Boost✔️✔️:     {
    // Boost✔️✔️:         vis.start_vertex(start_vertex, g);
    // Boost✔️✔️:         detail::depth_first_visit_impl(
    // Boost✔️✔️:             g, start_vertex, vis, color, detail::nontruth2());
    // Boost✔️✔️:     }
    // Boost✔️✔️:
    // Boost✔️✔️:     for (boost::tie(ui, ui_end) = vertices(g); ui != ui_end; ++ui)
    // Boost✔️✔️:     {
    // Boost✔️✔️:         Vertex u = implicit_cast< Vertex >(*ui);
    // Boost✔️✔️:         ColorValue u_color = get(color, u);
    // Boost✔️✔️:         if (u_color == Color::white())
    // Boost✔️✔️:         {
    // Boost✔️✔️:             vis.start_vertex(u, g);
    // Boost✔️✔️:             detail::depth_first_visit_impl(
    // Boost✔️✔️:                 g, u, vis, color, detail::nontruth2());
    // Boost✔️✔️:         }
    // Boost✔️✔️:     }
    // Boost✔️✔️: }
    // END REACHED BOOST connected_components / recorder / depth_first_search
    // vecS vertices start in index order; each white root increments the
    // component counter, and discovery records that counter. An explicit
    // (vertex, next-out-edge) stack preserves native DFS adjacency order.
    // The recorder has no finish/edge callback effects, so the component
    // label itself represents both gray and black without a second color map.
    // O(V+E) discovery and O(V) stack/labels, no graph/query/atom copies.
    // Group rows are the companion getMolFrags overload's output, assembled
    // only after discovery in atom-index order using dense component labels.
    let mut atom_to_component = vec![usize::MAX; source.atom_count()];
    let mut components = Vec::new();
    for start in 0..source.atom_count() {
        if atom_to_component[start] != usize::MAX {
            continue;
        }
        let component = components.len();
        let mut stack = vec![(start, 0_usize)];
        atom_to_component[start] = component;
        while let Some((atom, next_neighbor)) = stack.last_mut() {
            if *next_neighbor == source.neighbor_count(*atom) {
                stack.pop();
                continue;
            }
            let neighbor = source.neighbor_at(*atom, *next_neighbor);
            *next_neighbor += 1;
            if atom_to_component[neighbor] == usize::MAX {
                atom_to_component[neighbor] = component;
                stack.push((neighbor, 0));
            }
        }
        components.push(Vec::new());
    }
    // BEGIN RDKIT CPP FUNCTION MolOps::getMolFrags(VECT_INT_VECT&)
    // RDKit✔️🔝: unsigned int getMolFrags(const ROMol &mol, VECT_INT_VECT &frags) {
    // RDKit✔️🔝:   frags.clear();
    // RDKit✔️🔝:   INT_VECT mapping;
    // RDKit✔️🔝:   getMolFrags(mol, mapping);
    // RDKit✔️🔝:
    // RDKit✔️🔝:   INT_INT_VECT_MAP comMap;
    // RDKit✔️🔝:   for (unsigned int i = 0; i < mol.getNumAtoms(); i++) {
    // RDKit✔️🔝:     int mi = mapping[i];
    // RDKit✔️🔝:     if (comMap.find(mi) == comMap.end()) {
    // RDKit✔️🔝:       INT_VECT comp;
    // RDKit✔️🔝:       comMap[mi] = comp;
    // RDKit✔️🔝:     }
    // RDKit✔️🔝:     comMap[mi].push_back(i);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:
    // RDKit✔️🔝:   for (INT_INT_VECT_MAP_CI mci = comMap.begin(); mci != comMap.end(); mci++) {
    // RDKit✔️🔝:     frags.push_back((*mci).second);
    // RDKit✔️🔝:   }
    // RDKit✔️🔝:   return rdcast<unsigned int>(frags.size());
    // RDKit✔️🔝: }
    // END RDKIT CPP FUNCTION MolOps::getMolFrags(VECT_INT_VECT&)
    // Behavior: DFS assigns dense labels 0..K in ascending-root order. Source
    // comMap iteration therefore equals this vector's index order, and its
    // push_back loop equals ascending atom enumeration below. Fresh detached
    // groups represent the source-cleared output without any old entries.
    // Complexity improvement: dense component IDs replace O(log K) tree-map
    // lookups/node allocations with O(1) indexing, preserving all member and
    // fragment ordering. One O(V) pass, required output only, no atom copies.
    for (atom, component) in atom_to_component.iter().copied().enumerate() {
        components[component].push(AtomId::new(atom));
    }
    ConnectedComponents {
        atom_to_component,
        components,
    }
}

pub fn connected_components(topology: &TopologyBlock) -> Result<ConnectedComponents, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    Ok(connected_components_from_source(topology))
}

/// Label query components using the canonical graph algorithm without
/// materializing element-only topology or discarding query carrier identities.
#[doc(hidden)]
pub fn query_connected_components(query: &QueryGraph) -> Result<ConnectedComponents, PathError> {
    query.validate().map_err(PathError::QueryGraph)?;
    Ok(connected_components_from_source(query))
}

pub fn shortest_path(
    topology: &TopologyBlock,
    begin: AtomId,
    end: AtomId,
) -> Result<Vec<AtomId>, PathError> {
    // BEGIN RDKIT CPP FUNCTION MolOps::getShortestPath
    // RDKit✔️✔️: INT_LIST getShortestPath(const ROMol &mol, int aid1, int aid2) {
    // RDKit✔️✔️:   int nats = mol.getNumAtoms();
    // RDKit✔️✔️:   RANGE_CHECK(0, aid1, nats - 1);
    // RDKit✔️✔️:   RANGE_CHECK(0, aid2, nats - 1);
    // RDKit✔️✔️:   CHECK_INVARIANT(aid1 != aid2, "");
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_required_atom(topology, begin, "begin")?;
    validate_required_atom(topology, end, "end")?;
    if begin == end {
        return Err(PathError::EqualShortestPathEndpoints { atom: begin });
    }
    // RDKit✔️✔️:   INT_VECT pred(nats, -1);  // set all atoms to unprocessed state
    // RDKit✔️✔️:   pred[aid1] = -2;          // marks begin
    // RDKit✔️✔️:   pred[aid2] = -3;          // marks end
    // RDKit✔️✔️:   std::deque<int> bfsQ;
    // RDKit✔️✔️:   bfsQ.push_back(aid1);
    let mut predecessor = vec![None; topology.atoms.len()];
    let mut visited = vec![false; topology.atoms.len()];
    let mut queue = VecDeque::from([begin.index()]);
    visited[begin.index()] = true;
    let mut found = false;
    // RDKit✔️✔️:   while ((!done) && (bfsQ.size() > 0)) {
    // RDKit✔️✔️:     int curAid = bfsQ.front();
    // RDKit✔️✔️:     boost::tie(nbrIdx, endNbrs) =
    // RDKit✔️✔️:         mol.getAtomNeighbors(mol.getAtomWithIdx(curAid));
    while !found {
        let Some(current) = queue.pop_front() else {
            break;
        };
        for neighbor in topology.adjacency.neighbors_of(current) {
            let next = neighbor.atom_index;
            if next == end.index() {
                predecessor[next] = Some(current);
                found = true;
                break;
            }
            if !visited[next] {
                visited[next] = true;
                predecessor[next] = Some(current);
                queue.push_back(next);
            }
        }
    }
    // RDKit✔️✔️:   INT_LIST res;
    // RDKit✔️✔️:   if (done) {
    // RDKit✔️✔️:     int prev = aid2;
    // RDKit✔️✔️:     res.push_back(aid2);
    // RDKit✔️✔️:     while (!done) {
    // RDKit✔️✔️:       prev = pred[prev];
    // RDKit✔️✔️:       if (prev != aid1) { res.push_front(prev); }
    // RDKit✔️✔️:       else { done = true; }
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:     res.push_front(aid1);
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION MolOps::getShortestPath
    if !found {
        return Ok(Vec::new());
    }
    let mut reversed = vec![end];
    let mut current = end.index();
    while current != begin.index() {
        current = predecessor[current].expect("a found BFS path has predecessors");
        reversed.push(AtomId::new(current));
    }
    reversed.reverse();
    Ok(reversed)
}

pub fn all_paths_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &PathSearchParams,
) -> Result<Vec<GraphPath>, PathError> {
    Ok(
        all_paths_in_range(topology, target_length, target_length, params)?
            .remove(&target_length)
            .unwrap_or_default(),
    )
}

pub fn all_paths_in_range(
    topology: &TopologyBlock,
    lower_length: usize,
    upper_length: usize,
    params: &PathSearchParams,
) -> Result<BTreeMap<usize, Vec<GraphPath>>, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_range(lower_length, upper_length)?;
    let distances = params
        .only_shortest_paths
        .then(|| crate::matrices::unweighted_distance_steps(topology))
        .transpose()
        .map_err(PathError::Matrix)?;
    all_paths_from_graph(
        PathGraphAccess::Concrete(topology),
        lower_length,
        upper_length,
        params,
        distances.as_deref(),
    )
}

/// Internal detached query entry for source bond paths with real query atomic numbers.
#[doc(hidden)]
pub fn query_bond_paths_in_range(
    query: &QueryGraph,
    lower_length: usize,
    upper_length: usize,
    params: &SubgraphSearchParams,
) -> Result<BTreeMap<usize, Vec<GraphPath>>, PathError> {
    query.validate().map_err(PathError::QueryGraph)?;
    let path_params = PathSearchParams {
        representation: PathRepresentation::Bonds,
        use_hydrogens: params.use_hydrogens,
        rooted_at_atom: params.rooted_at_atom,
        only_shortest_paths: false,
        ignore_atoms: params.ignore_atoms,
    };
    all_paths_from_graph(
        PathGraphAccess::Query(query),
        lower_length,
        upper_length,
        &path_params,
        None,
    )
}

fn all_paths_from_graph(
    graph: PathGraphAccess<'_>,
    lower_length: usize,
    upper_length: usize,
    params: &PathSearchParams,
    distances: Option<&[usize]>,
) -> Result<BTreeMap<usize, Vec<GraphPath>>, PathError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::findAllPathsOfLengthsMtoN (Release_2026_03_6)
    // RDKit❗✔️: INT_PATH_LIST_MAP
    // RDKit❗✔️: findAllPathsOfLengthsMtoN(const ROMol &mol, unsigned int lowerLen,
    // RDKit❗✔️:                           unsigned int upperLen, bool useBonds, bool useHs,
    // RDKit❗✔️:                           int rootedAtAtom, bool onlyShortestPaths,
    // RDKit❗✔️:                           boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //  We can't be clever here and just use the bond adjacency matrix
    // RDKit❗✔️:   //  to solve this problem when useBonds is true.  This is because
    // RDKit❗✔️:   //  the bond adjacency matrices for the molecules C1CC1 and CC(C)C
    // RDKit❗✔️:   //  are indistinguishable.  In the second case, t-butane (and
    // RDKit❗✔️:   //  anything else with a T junction), we'll get some subgraphs mixed
    // RDKit❗✔️:   //  in with the paths.  So we have to construct paths of atoms and
    // RDKit❗✔️:   //  then convert them into bond paths.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   PRECONDITION(lowerLen <= upperLen, "");
    // RDKit❗✔️:
    // RDKit❗✔️:   // the molecule owns the distance matrix pointer (if we need to get it)
    // RDKit❗✔️:   double *distMat = onlyShortestPaths ? MolOps::getDistanceMat(mol) : nullptr;
    // RDKit❗✔️:   int *adjMat, dim;
    // RDKit❗✔️:   dim = mol.getNumAtoms();
    // RDKit❗✔️:   adjMat = new int[dim * dim];
    // RDKit❗✔️:   memset((void *)adjMat, 0, dim * dim * sizeof(int));
    // RDKit❗✔️:
    // RDKit❗✔️:   if (!distMat) {
    // RDKit❗✔️:     // generate the adjacency matrix by hand by looping over the bonds
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       Atom *beg = bond->getBeginAtom();
    // RDKit❗✔️:       Atom *end = bond->getEndAtom();
    // RDKit❗✔️:       // check for H, which we might be skipping
    // RDKit❗✔️:       if (useHs || (beg->getAtomicNum() != 1 && end->getAtomicNum() != 1)) {
    // RDKit❗✔️:         adjMat[beg->getIdx() * dim + end->getIdx()] = 1;
    // RDKit❗✔️:         adjMat[end->getIdx() * dim + beg->getIdx()] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     // if we have the distance matrix, we can just loop over that:
    // RDKit❗✔️:     for (auto i = 0; i < dim; ++i) {
    // RDKit❗✔️:       for (auto j = i + 1; j < dim; ++j) {
    // RDKit❗✔️:         if (fabs(distMat[i * dim + j] - 1) < 1e-4) {
    // RDKit❗✔️:           adjMat[i * dim + j] = 1;
    // RDKit❗✔️:           adjMat[j * dim + i] = 1;
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // if we're using bonds, we'll need to find paths of length N+1,
    // RDKit❗✔️:   // then convert them
    // RDKit❗✔️:   if (useBonds) {
    // RDKit❗✔️:     ++lowerLen;
    // RDKit❗✔️:     ++upperLen;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // find the paths themselves
    // RDKit❗✔️:   INT_PATH_LIST_MAP atomPaths = Subgraphs::pathFinderHelper(
    // RDKit❗✔️:       adjMat, dim, lowerLen, upperLen, rootedAtAtom, distMat, ignoreAtoms);
    // RDKit❗✔️:
    // RDKit❗✔️:   // clean up the adjacency matrix
    // RDKit❗✔️:   delete[] adjMat;
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_PATH_LIST_MAP res;
    // RDKit❗✔️:
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //--------------------------------------------------------
    // RDKit❗✔️:   // loop through all the paths we have and make sure that there are
    // RDKit❗✔️:   // no duplicates (duplicate = contains identical bond indices)
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //  We need to use the bond paths for this duplicate finding
    // RDKit❗✔️:   //  because, in rings, there can be many paths which share atom
    // RDKit❗✔️:   //  indices but which have different bond compositions. For example,
    // RDKit❗✔️:   //  there is only one "atom unique" path of length 5 bonds (6 atoms)
    // RDKit❗✔️:   //  through a 6-ring, but there are six bond paths.
    // RDKit❗✔️:   //
    // RDKit❗✔️:   if (!useBonds && lowerLen >= 1) {
    // RDKit❗✔️:     res[1] = atomPaths[1];
    // RDKit❗✔️:   }
    // RDKit❗✔️:   if (useBonds || upperLen > 1) {
    // RDKit❗✔️:     for (unsigned int i = lowerLen; i <= upperLen; ++i) {
    // RDKit❗✔️:       if (i <= 1) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:
    // RDKit❗✔️:       std::vector<boost::dynamic_bitset<>> invars;
    // RDKit❗✔️:
    // RDKit❗✔️:       for (PATH_LIST::const_iterator vivI = atomPaths[i].begin();
    // RDKit❗✔️:            vivI != atomPaths[i].end(); ++vivI) {
    // RDKit❗✔️:         boost::dynamic_bitset<> invar(mol.getNumBonds());
    // RDKit❗✔️:         const PATH_TYPE &resi = *vivI;
    // RDKit❗✔️:         PATH_TYPE locV;
    // RDKit❗✔️:         locV.reserve(i);
    // RDKit❗✔️:         for (unsigned int j = 0; j < i - 1; j++) {
    // RDKit❗✔️:           const Bond *bond = mol.getBondBetweenAtoms(resi[j], resi[j + 1]);
    // RDKit❗✔️:           locV.push_back(bond->getIdx());
    // RDKit❗✔️:           invar.set(bond->getIdx());
    // RDKit❗✔️:         }
    // RDKit❗✔️:         if (std::find(invars.begin(), invars.end(), invar) == invars.end()) {
    // RDKit❗✔️:           invars.push_back(invar);
    // RDKit❗✔️:           if (useBonds) {
    // RDKit❗✔️:             res[i - 1].push_back(locV);
    // RDKit❗✔️:           } else {
    // RDKit❗✔️:             res[i].push_back(resi);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::findAllPathsOfLengthsMtoN
    validate_range(lower_length, upper_length)?;
    validate_ignored_atom_mask(params.ignore_atoms, graph.atom_count())?;

    let adjacency =
        atom_adjacency_matrix(graph, params.only_shortest_paths || params.use_hydrogens);
    let (atom_lower, atom_upper) = match params.representation {
        PathRepresentation::Bonds => (
            lower_length
                .checked_add(1)
                .ok_or(PathError::LengthOverflow {
                    length: lower_length,
                })?,
            upper_length
                .checked_add(1)
                .ok_or(PathError::LengthOverflow {
                    length: upper_length,
                })?,
        ),
        PathRepresentation::Atoms => (lower_length, upper_length),
    };
    let atom_paths = path_finder_helper(
        &adjacency,
        graph.atom_count(),
        atom_lower,
        atom_upper,
        params.rooted_at_atom,
        distances,
        params.ignore_atoms,
    );

    let mut result = BTreeMap::new();
    if params.representation == PathRepresentation::Atoms && atom_lower >= 1 {
        result.insert(
            1,
            atom_paths
                .get(&1)
                .into_iter()
                .flatten()
                .map(|row| GraphPath::Atoms(row.iter().copied().map(AtomId::new).collect()))
                .collect(),
        );
    }
    if params.representation == PathRepresentation::Bonds || atom_upper > 1 {
        for length in atom_lower..=atom_upper {
            if length <= 1 {
                continue;
            }
            let mut seen_bond_sets: Vec<Vec<bool>> = Vec::new();
            for atom_path in atom_paths.get(&length).into_iter().flatten() {
                let mut bond_set = vec![false; graph.bond_count()];
                let mut bond_path = Vec::with_capacity(length - 1);
                for pair in atom_path.windows(2) {
                    let bond = graph_bond_between(graph, pair[0], pair[1])
                        .expect("validated adjacency path must have a bond");
                    bond_set[bond.index()] = true;
                    bond_path.push(bond);
                }
                if seen_bond_sets.contains(&bond_set) {
                    continue;
                }
                seen_bond_sets.push(bond_set);
                match params.representation {
                    PathRepresentation::Bonds => result
                        .entry(length - 1)
                        .or_insert_with(Vec::new)
                        .push(GraphPath::Bonds(bond_path)),
                    PathRepresentation::Atoms => result
                        .entry(length)
                        .or_insert_with(Vec::new)
                        .push(GraphPath::Atoms(
                            atom_path.iter().copied().map(AtomId::new).collect(),
                        )),
                }
            }
        }
    }
    Ok(result)
}

pub fn all_subgraphs_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &SubgraphSearchParams,
) -> Result<Vec<Vec<BondId>>, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_ignored_atom_mask(params.ignore_atoms, topology.atoms.len())?;
    if target_length == 0 {
        return Ok(Vec::new());
    }
    let neighbors = bond_neighbor_map(PathGraphAccess::Concrete(topology), params.use_hydrogens);
    Ok(all_subgraphs_of_length_from_neighbors(
        topology,
        &neighbors,
        target_length,
        params.rooted_at_atom,
        params.ignore_atoms,
    ))
}

pub fn all_subgraphs_in_range(
    topology: &TopologyBlock,
    lower_length: usize,
    upper_length: usize,
    params: &SubgraphSearchParams,
) -> Result<BTreeMap<usize, Vec<Vec<BondId>>>, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    all_subgraphs_from_graph(
        PathGraphAccess::Concrete(topology),
        lower_length,
        upper_length,
        params,
    )
}

/// Internal query subgraph entry; no concrete atom conversion or query copy.
#[doc(hidden)]
pub fn query_subgraphs_in_range(
    query: &QueryGraph,
    lower_length: usize,
    upper_length: usize,
    params: &SubgraphSearchParams,
) -> Result<BTreeMap<usize, Vec<Vec<BondId>>>, PathError> {
    query.validate().map_err(PathError::QueryGraph)?;
    all_subgraphs_from_graph(
        PathGraphAccess::Query(query),
        lower_length,
        upper_length,
        params,
    )
}

fn all_subgraphs_from_graph(
    graph: PathGraphAccess<'_>,
    lower_length: usize,
    upper_length: usize,
    params: &SubgraphSearchParams,
) -> Result<BTreeMap<usize, Vec<Vec<BondId>>>, PathError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::findAllSubgraphsOfLengthsMtoN (Release_2026_03_6)
    // RDKit❗✔️: INT_PATH_LIST_MAP findAllSubgraphsOfLengthsMtoN(
    // RDKit❗✔️:     const ROMol &mol, unsigned int lowerLen, unsigned int upperLen, bool useHs,
    // RDKit❗✔️:     int rootedAtAtom, boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   PRECONDITION(lowerLen <= upperLen, "");
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == mol.getNumAtoms(),
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   boost::dynamic_bitset<> forbidden(mol.getNumBonds());
    // RDKit❗✔️:   // if there are any ignore atoms, mark any bonds involving them as forbidden
    // RDKit❗✔️:   if (ignoreAtoms) {
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       if (ignoreAtoms->test(bond->getBeginAtomIdx()) ||
    // RDKit❗✔️:           ignoreAtoms->test(bond->getEndAtomIdx())) {
    // RDKit❗✔️:         forbidden[bond->getIdx()] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   INT_INT_VECT_MAP nbrs;
    // RDKit❗✔️:   Subgraphs::getNbrsList(mol, useHs, nbrs);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Start path at each bond
    // RDKit❗✔️:   INT_PATH_LIST_MAP res;
    // RDKit❗✔️:   for (unsigned int idx = lowerLen; idx <= upperLen; idx++) {
    // RDKit❗✔️:     PATH_LIST ordern;
    // RDKit❗✔️:     res[idx] = ordern;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // start paths at each bond:
    // RDKit❗✔️:   for (auto nbi = nbrs.begin(); nbi != nbrs.end(); nbi++) {
    // RDKit❗✔️:     int i = (*nbi).first;
    // RDKit❗✔️:     if (forbidden[i]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // if we're only returning paths rooted at a particular atom, check now
    // RDKit❗✔️:     // that this bond involves that atom:
    // RDKit❗✔️:     if (rootedAtAtom >= 0 &&
    // RDKit❗✔️:         mol.getBondWithIdx(i)->getBeginAtomIdx() !=
    // RDKit❗✔️:             static_cast<unsigned int>(rootedAtAtom) &&
    // RDKit❗✔️:         mol.getBondWithIdx(i)->getEndAtomIdx() !=
    // RDKit❗✔️:             static_cast<unsigned int>(rootedAtAtom)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // don't come back to this bond in the later subgraphs
    // RDKit❗✔️:     forbidden[i] = 1;
    // RDKit❗✔️:
    // RDKit❗✔️:     // start the recursive path building with the current bond
    // RDKit❗✔️:     PATH_TYPE spath;
    // RDKit❗✔️:     spath.clear();
    // RDKit❗✔️:     spath.push_back(i);
    // RDKit❗✔️:
    // RDKit❗✔️:     // neighbors of this bond are the next candidates
    // RDKit❗✔️:     INT_VECT cands = nbrs[i];
    // RDKit❗✔️:
    // RDKit❗✔️:     // now call the recursive function
    // RDKit❗✔️:     // little bit different from the python version
    // RDKit❗✔️:     // the result list of paths is passed as a reference, instead of on the fly
    // RDKit❗✔️:     // appending
    // RDKit❗✔️:     Subgraphs::recurseWalkRange(nbrs, spath, cands, lowerLen, upperLen,
    // RDKit❗✔️:                                 forbidden, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   nbrs.clear();
    // RDKit❗✔️:   return res;  // FIX : need some verbose testing code here
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::findAllSubgraphsOfLengthsMtoN
    validate_range(lower_length, upper_length)?;
    validate_ignored_atom_mask(params.ignore_atoms, graph.atom_count())?;
    let neighbors = bond_neighbor_map(graph, params.use_hydrogens);
    let mut result = (lower_length..=upper_length)
        .map(|length| (length, Vec::new()))
        .collect::<BTreeMap<_, _>>();
    if upper_length == 0 {
        return Ok(result);
    }
    let mut forbidden = ignored_bonds(graph, params.ignore_atoms);
    for (&start, adjacent) in &neighbors {
        if forbidden[start] {
            continue;
        }
        if !graph_root_allows_bond(graph, params.rooted_at_atom, start) {
            continue;
        }
        forbidden[start] = true;
        recurse_walk_range(
            &neighbors,
            vec![start],
            adjacent.clone(),
            lower_length,
            upper_length,
            forbidden.clone(),
            &mut result,
        );
    }
    Ok(result)
}

pub fn unique_subgraphs_of_length(
    topology: &TopologyBlock,
    target_length: usize,
    params: &UniqueSubgraphParams,
) -> Result<Vec<Vec<BondId>>, PathError> {
    // BEGIN RDKIT CPP FUNCTION RDKit::findUniqueSubgraphsOfLengthN (Release_2026_03_6)
    // RDKit❗✔️: PATH_LIST findUniqueSubgraphsOfLengthN(const ROMol &mol, unsigned int targetLen,
    // RDKit❗✔️:                                        bool useHs, bool useBO, int rootedAtAtom,
    // RDKit❗✔️:                                        boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   // start by finding all subgraphs, then uniquify
    // RDKit❗✔️:   PATH_LIST allSubgraphs = findAllSubgraphsOfLengthN(mol, targetLen, useHs,
    // RDKit❗✔️:                                                      rootedAtAtom, ignoreAtoms);
    // RDKit❗✔️:   PATH_LIST res = Subgraphs::uniquifyPaths(mol, allSubgraphs, useBO);
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::findUniqueSubgraphsOfLengthN
    topology.validate().map_err(PathError::InvalidTopology)?;
    if let Some(extra) = &params.extra_atom_invariants
        && extra.len() != topology.atoms.len()
    {
        return Err(PathError::ExtraInvariantLength {
            actual: extra.len(),
            expected: topology.atoms.len(),
        });
    }
    let all = all_subgraphs_of_length(
        topology,
        target_length,
        &SubgraphSearchParams {
            use_hydrogens: params.use_hydrogens,
            rooted_at_atom: params.rooted_at_atom,
            ignore_atoms: params.ignore_atoms,
        },
    )?;
    let mut result = Vec::new();
    let mut seen = Vec::new();
    for path in all {
        let discriminator = path_discriminator(
            topology,
            &path,
            params.use_bond_orders,
            params.extra_atom_invariants.as_deref(),
        )?;
        if !seen.contains(&discriminator) {
            seen.push(discriminator);
            result.push(path);
        }
    }
    Ok(result)
}

pub fn atom_environment(
    topology: &TopologyBlock,
    radius: usize,
    root: AtomId,
    params: &AtomEnvironmentParams,
) -> Result<AtomEnvironment, PathError> {
    // BEGIN RDKIT CPP FUNCTION findAtomEnvironmentOfRadiusN
    // RDKit✔️✔️: if (rootedAtAtom >= mol.getNumAtoms()) {
    // RDKit✔️✔️:   throw ValueErrorException("bad atom index");
    // RDKit✔️✔️: }
    topology.validate().map_err(PathError::InvalidTopology)?;
    validate_required_atom(topology, root, "root")?;
    let mut atom_distances = vec![None; topology.atoms.len()];
    atom_distances[root.index()] = Some(0);
    if radius == 0 {
        return Ok(AtomEnvironment {
            bonds: Vec::new(),
            atom_distances,
        });
    }
    // RDKit✔️✔️: std::list<std::pair<int, int>> nbrStack;
    let mut layer = VecDeque::new();
    for neighbor in topology.adjacency.neighbors_of(root.index()) {
        if params.use_hydrogens || topology.atoms[neighbor.atom_index].atomic_number() != 1 {
            layer.push_back((root.index(), neighbor.bond));
        }
    }
    let mut bonds_in = vec![false; topology.bonds.len()];
    let mut bonds = Vec::new();
    let mut completed_layers = 0usize;
    for depth in 0..radius {
        if layer.is_empty() {
            break;
        }
        let mut next_layer = VecDeque::new();
        while let Some((start_atom, bond_id)) = layer.pop_front() {
            if bonds_in[bond_id.index()] {
                continue;
            }
            bonds_in[bond_id.index()] = true;
            bonds.push(bond_id);
            let bond = &topology.bonds[bond_id.index()];
            let other = if bond.begin().index() == start_atom {
                bond.end().index()
            } else {
                bond.begin().index()
            };
            let distance = depth + 1;
            atom_distances[other] =
                Some(atom_distances[other].map_or(distance, |old| old.min(distance)));
            if depth < radius - 1 {
                for neighbor in topology.adjacency.neighbors_of(other) {
                    if !bonds_in[neighbor.bond.index()]
                        && (params.use_hydrogens
                            || topology.atoms[neighbor.atom_index].atomic_number() != 1)
                    {
                        next_layer.push_back((other, neighbor.bond));
                    }
                }
            }
        }
        layer = next_layer;
        completed_layers += 1;
    }
    // RDKit✔️✔️: if (i != radius && enforceSize) {
    // RDKit✔️✔️:   res.clear();
    // RDKit✔️✔️:   if (atomMap) { atomMap->clear(); }
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION findAtomEnvironmentOfRadiusN
    if completed_layers != radius && params.enforce_radius {
        bonds.clear();
        atom_distances.fill(None);
    }
    Ok(AtomEnvironment {
        bonds,
        atom_distances,
    })
}

pub fn bond_ids_from_atom_path(
    topology: &TopologyBlock,
    atoms: &[AtomId],
) -> Result<Vec<BondId>, PathError> {
    // BEGIN RDKIT CPP FUNCTION bondListFromAtomList
    // RDKit✔️✔️: PATH_TYPE bondListFromAtomList(const ROMol &mol, const PATH_TYPE &atomIds) {
    // RDKit✔️✔️:   PATH_TYPE bids;
    // RDKit✔️✔️:   unsigned int natms = atomIds.size();
    // RDKit✔️✔️:   if (natms <= 1) { return bids; }
    topology.validate().map_err(PathError::InvalidTopology)?;
    for (position, atom) in atoms.iter().copied().enumerate() {
        if atom.index() >= topology.atoms.len() {
            return Err(PathError::AtomPathOutOfRange {
                position,
                atom,
                atom_count: topology.atoms.len(),
            });
        }
    }
    let mut result = Vec::new();
    // RDKit✔️✔️:   for (unsigned int i = 0; i < natms; i++) {
    // RDKit✔️✔️:     for (unsigned int j = i + 1; j < natms; j++) {
    // RDKit✔️✔️:       const Bond *bnd = mol.getBondBetweenAtoms(atomIds[i], atomIds[j]);
    // RDKit✔️✔️:       if (bnd) { bids.push_back(bnd->getIdx()); }
    for i in 0..atoms.len() {
        for j in i + 1..atoms.len() {
            if let Some(bond) = bond_between(topology, atoms[i].index(), atoms[j].index()) {
                result.push(bond);
            }
        }
    }
    // RDKit✔️✔️:   return bids;
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION bondListFromAtomList
    Ok(result)
}

pub fn subtopology_from_path(
    topology: &TopologyBlock,
    bonds: &[BondId],
    params: &SubtopologyParams,
) -> Result<SubtopologyResult, PathError> {
    topology.validate().map_err(PathError::InvalidTopology)?;
    // BEGIN RDKIT CPP FUNCTION getSubsetInfo/copySelectedAtomsAndBonds
    // RDKit✔️✔️: } else if (options.method == SubsetMethod::BONDS) {
    // RDKit✔️✔️:   for (const auto &bond_idx : path) {
    // RDKit✔️✔️:     if (bond_idx < num_bonds) {
    // RDKit✔️✔️:       selectedBonds.set(bond_idx);
    // RDKit✔️✔️:       const auto &bnd = mol.getBondWithIdx(bond_idx);
    // RDKit✔️✔️:       selectedAtoms.set(bnd->getBeginAtomIdx());
    // RDKit✔️✔️:       selectedAtoms.set(bnd->getEndAtomIdx());
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    let mut selected_bonds = vec![false; topology.bonds.len()];
    let mut selected_atoms = vec![false; topology.atoms.len()];
    for bond_id in bonds.iter().copied() {
        let Some(bond) = topology.bonds.get(bond_id.index()) else {
            continue;
        };
        selected_bonds[bond_id.index()] = true;
        selected_atoms[bond.begin().index()] = true;
        selected_atoms[bond.end().index()] = true;
    }

    let mut atom_old_to_new = vec![None; topology.atoms.len()];
    let mut atom_new_to_old = Vec::new();
    let mut stereo_atom_mapping = BTreeMap::new();
    let mut copied_atoms = Vec::new();
    // RDKit✔️✔️: for (const auto &ref_atom : reference_mol.atoms()) {
    // RDKit✔️✔️:   if (!selectedAtoms[ref_atom->getIdx()]) { continue; }
    // RDKit✔️✔️:   extracted_atom->clearComputedProps();
    for atom in &topology.atoms {
        if !selected_atoms[atom.id().index()] {
            continue;
        }
        let new_id = AtomId::new(copied_atoms.len());
        atom_old_to_new[atom.id().index()] = Some(new_id);
        stereo_atom_mapping.insert(atom.id(), new_id);
        atom_new_to_old.push(Some(atom.id()));
        let mut copied = atom.clone().with_id(new_id);
        copied.clear_computed_props()?;
        copied_atoms.push(copied);
    }

    let mut bond_old_to_new = vec![None; topology.bonds.len()];
    let mut bond_new_to_old = Vec::new();
    let mut copied_bonds = Vec::new();
    // RDKit✔️✔️: for (const auto &ref_bond : reference_mol.bonds()) {
    // RDKit✔️✔️:   if (!selectedBonds[ref_bond->getIdx()]) { continue; }
    for bond in &topology.bonds {
        if !selected_bonds[bond.id().index()] {
            continue;
        }
        let new_id = BondId::new(copied_bonds.len());
        let begin =
            atom_old_to_new[bond.begin().index()].expect("selected bond endpoints are selected");
        let end =
            atom_old_to_new[bond.end().index()].expect("selected bond endpoints are selected");
        // RDKit✔️✔️: handleBondStereo(*extracted_bond, *ref_bond, reference_mol, atomMapping);
        // The sole .6 Subset::handleBondStereo owner in fragments.rs carries
        // its complete source body and getOtherAtomIdx body. Reuse it for
        // pathToSubmol's copyMolSubset route before remapping bond endpoints.
        let (stereo, stereo_atoms) =
            crate::fragments::handle_subset_bond_stereo(topology, bond, &stereo_atom_mapping);
        let mut copied = bond.clone();
        copied.clear_computed_props()?;
        copied = copied.remapped(new_id, begin, end, stereo_atoms);
        copied.set_stereo(stereo)?;
        bond_old_to_new[bond.id().index()] = Some(new_id);
        bond_new_to_old.push(Some(bond.id()));
        copied_bonds.push(copied);
    }

    let mapping = TopologyMapping {
        atoms: AtomMapping {
            old_to_new: atom_old_to_new,
            new_to_old: atom_new_to_old,
        },
        bonds: BondMapping {
            old_to_new: bond_old_to_new,
            new_to_old: bond_new_to_old,
        },
    };
    mapping
        .validate_for_counts(
            topology.atoms.len(),
            copied_atoms.len(),
            topology.bonds.len(),
            copied_bonds.len(),
        )
        .map_err(PathError::Mapping)?;

    for atom in &mut copied_atoms {
        let carrier = mapping.atoms.new_to_old[atom.id().index()]
            .expect("copied atoms always have a source row");
        atom.remap_template_attachment_order(&mapping.atoms.old_to_new)
            .map_err(|source| PathError::TemplateAttachmentRemap { carrier, source })?;
    }

    let substance_groups = remap_selected_substance_groups(topology, &mapping);
    let stereo_groups = remap_selected_stereo_groups(topology, &mapping)?;

    let subgraph = if params.copy_as_query {
        // RDKit✔️✔️: std::unique_ptr<Atom> extracted_atom{
        // RDKit✔️✔️:     options.copyAsQuery ? new QueryAtom(*ref_atom) : ref_atom->copy()};
        // RDKit✔️✔️: std::unique_ptr<Bond> extracted_bond{
        // RDKit✔️✔️:     options.copyAsQuery ? new QueryBond(*ref_bond) : ref_bond->copy()};
        let query_atoms = copied_atoms
            .into_iter()
            .map(|atom| {
                let predicate =
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number()));
                QueryAtom::from_parts(atom, predicate)
            })
            .collect();
        let query_bonds = copied_bonds
            .into_iter()
            .map(|bond| {
                let predicate = if bond.order() == cosmolkit_model::BondOrder::Unspecified {
                    QueryNode::predicate(BondQueryPredicate::Any)
                } else {
                    QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
                };
                QueryBond::from_parts(bond, predicate)
            })
            .collect();
        let graph = QueryGraph::from_parts(
            query_atoms,
            query_bonds,
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            stereo_groups,
        )
        .map_err(PathError::QueryGraph)?;
        DetachedPathSubgraph::Query {
            graph,
            substance_groups,
        }
    } else {
        let concrete = TopologyBlock::try_from_parts(
            copied_atoms,
            copied_bonds,
            substance_groups,
            stereo_groups,
        )
        .map_err(|error| match error {
            cosmolkit_model::TopologyValidationError::StereoGroup(cause) => {
                PathError::StereoGroup(cause)
            }
            error => PathError::InvalidSubsetTopology(error),
        })?;
        DetachedPathSubgraph::Concrete(concrete)
    };
    // END RDKIT CPP FUNCTION getSubsetInfo/copySelectedAtomsAndBonds
    Ok(SubtopologyResult { subgraph, mapping })
}

fn validate_required_atom(
    topology: &TopologyBlock,
    atom: AtomId,
    role: &'static str,
) -> Result<(), PathError> {
    if atom.index() >= topology.atoms.len() {
        Err(PathError::AtomOutOfRange {
            role,
            atom,
            atom_count: topology.atoms.len(),
        })
    } else {
        Ok(())
    }
}

fn validate_range(lower: usize, upper: usize) -> Result<(), PathError> {
    if lower > upper {
        Err(PathError::InvalidLengthRange { lower, upper })
    } else {
        Ok(())
    }
}

fn atom_adjacency_matrix(graph: PathGraphAccess<'_>, include_hydrogens: bool) -> Vec<bool> {
    // RDKit✔️✔️: for (bondIt = mol.beginBonds(); bondIt != mol.endBonds(); bondIt++) {
    // RDKit✔️✔️:   Atom *beg = (*bondIt)->getBeginAtom();
    // RDKit✔️✔️:   Atom *end = (*bondIt)->getEndAtom();
    // RDKit✔️✔️:   if (useHs || (beg->getAtomicNum() != 1 && end->getAtomicNum() != 1)) {
    // RDKit✔️✔️:     adjMat[beg->getIdx() * dim + end->getIdx()] = 1;
    // RDKit✔️✔️:     adjMat[end->getIdx() * dim + beg->getIdx()] = 1;
    // RDKit✔️✔️:   }
    let dimension = graph.atom_count();
    let mut adjacency = vec![false; dimension.saturating_mul(dimension)];
    for index in 0..graph.bond_count() {
        let bond = graph.bond(index);
        let begin = bond.begin().index();
        let end = bond.end().index();
        if include_hydrogens || (graph.atomic_number(begin) != 1 && graph.atomic_number(end) != 1) {
            adjacency[begin * dimension + end] = true;
            adjacency[end * dimension + begin] = true;
        }
    }
    adjacency
}

fn path_finder_helper(
    adjacency: &[bool],
    dimension: usize,
    minimum_length: usize,
    maximum_length: usize,
    root: Option<AtomId>,
    distances: Option<&[usize]>,
    ignore_atoms: Option<&[bool]>,
) -> BTreeMap<usize, Vec<Vec<usize>>> {
    // BEGIN RDKIT CPP FUNCTION RDKit::Subgraphs::pathFinderHelper (Release_2026_03_6)
    // RDKit❗✔️: INT_PATH_LIST_MAP
    // RDKit❗✔️: pathFinderHelper(int *adjMat, unsigned int dim, unsigned int minLen,
    // RDKit❗✔️:                  unsigned int maxLen, int rootedAtAtom, double *distMat,
    // RDKit❗✔️:                  boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit❗✔️:   PRECONDITION(minLen <= maxLen, "bad lengths provided");
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == dim,
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   // finds all paths of length N using an adjacency matrix,
    // RDKit❗✔️:   //  which is constructed elsewhere
    // RDKit❗✔️:   INT_PATH_LIST_MAP res;
    // RDKit❗✔️:   PATH_LIST paths;
    // RDKit❗✔️:   paths.clear();
    // RDKit❗✔️:
    // RDKit❗✔️:   if (rootedAtAtom < 0) {
    // RDKit❗✔️:     // start a path at each possible index
    // RDKit❗✔️:     for (unsigned int i = 0; i < dim; i++) {
    // RDKit❗✔️:       if (ignoreAtoms && ignoreAtoms->test(i)) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       PATH_TYPE tPath;
    // RDKit❗✔️:       tPath.push_back(i);
    // RDKit❗✔️:       paths.push_back(tPath);
    // RDKit❗✔️:     }
    // RDKit❗✔️:   } else if (rootedAtAtom < static_cast<int>(dim) &&
    // RDKit❗✔️:              (!ignoreAtoms || !ignoreAtoms->test(rootedAtAtom))) {
    // RDKit❗✔️:     // only start a path at the atom of interest:
    // RDKit❗✔️:     PATH_TYPE tPath;
    // RDKit❗✔️:     tPath.push_back(rootedAtAtom);
    // RDKit❗✔️:     paths.push_back(tPath);
    // RDKit❗✔️:   } else {
    // RDKit❗✔️:     return res;
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // and build them up one index at a time:
    // RDKit❗✔️:   for (unsigned int length = 1; length < maxLen; length++) {
    // RDKit❗✔️:     // extend each path:
    // RDKit❗✔️:     if (length >= minLen) {
    // RDKit❗✔️:       res[length] = paths;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     paths = extendPaths(adjMat, dim, paths, maxLen, distMat, ignoreAtoms);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   res[maxLen] = paths;
    // RDKit❗✔️:
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::Subgraphs::pathFinderHelper
    let mut paths = match root {
        None => (0..dimension)
            .filter(|&atom| !ignore_atoms.is_some_and(|mask| mask[atom]))
            .map(|atom| vec![atom])
            .collect(),
        Some(root) if root.index() < dimension => {
            if ignore_atoms.is_some_and(|mask| mask[root.index()]) {
                return BTreeMap::new();
            }
            vec![vec![root.index()]]
        }
        Some(_) => return BTreeMap::new(),
    };
    let mut result = BTreeMap::new();
    for length in 1..maximum_length {
        if length >= minimum_length {
            result.insert(length, paths.clone());
        }
        paths = extend_paths(
            adjacency,
            dimension,
            &paths,
            maximum_length,
            distances,
            ignore_atoms,
        );
    }
    result.insert(maximum_length, paths);
    result
}

fn extend_paths(
    adjacency: &[bool],
    dimension: usize,
    paths: &[Vec<usize>],
    allow_ring_closures: usize,
    distances: Option<&[usize]>,
    ignore_atoms: Option<&[bool]>,
) -> Vec<Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION RDKit::Subgraphs::extendPaths (Release_2026_03_6)
    // RDKit❗✔️: PATH_LIST
    // RDKit❗✔️: extendPaths(int *adjMat, unsigned int dim, const PATH_LIST &paths,
    // RDKit❗✔️:             int allowRingClosures, double *distMat,
    // RDKit❗✔️:             boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   PRECONDITION(adjMat, "no matrix");
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == dim,
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   //
    // RDKit❗✔️:   //  extend each of the currently active paths by adding
    // RDKit❗✔️:   //   a single adjacent index to the end of each
    // RDKit❗✔️:   //
    // RDKit❗✔️:   PATH_LIST res;
    // RDKit❗✔️:   PATH_LIST::const_iterator path;
    // RDKit❗✔️:   for (path = paths.begin(); path != paths.end(); ++path) {
    // RDKit❗✔️:     unsigned int endIdx = (*path)[path->size() - 1];
    // RDKit❗✔️:     unsigned int iTab = endIdx * dim;
    // RDKit❗✔️:     for (unsigned int otherIdx = 0; otherIdx < dim; otherIdx++) {
    // RDKit❗✔️:       if (ignoreAtoms && ignoreAtoms->test(otherIdx)) {
    // RDKit❗✔️:         continue;
    // RDKit❗✔️:       }
    // RDKit❗✔️:       if (adjMat[iTab + otherIdx] == 1) {
    // RDKit❗✔️:         if (distMat &&
    // RDKit❗✔️:             distMat[path->front() * dim + otherIdx] - path->size() < -0.001) {
    // RDKit❗✔️:           continue;
    // RDKit❗✔️:         }
    // RDKit❗✔️:         // test 1: make sure the new atom is not already
    // RDKit❗✔️:         //   in the path
    // RDKit❗✔️:         auto loc =
    // RDKit❗✔️:             std::find(path->begin(), path->end(), static_cast<int>(otherIdx));
    // RDKit❗✔️:         // The two conditions for adding the atom are:
    // RDKit❗✔️:         //   1) it's not there already
    // RDKit❗✔️:         //   2) it's there, but ring closures are allowed and this
    // RDKit❗✔️:         //      will be the last addition to the path.
    // RDKit❗✔️:         if (loc == path->end()) {
    // RDKit❗✔️:           // the easy case
    // RDKit❗✔️:           // PATH_TYPE newPath=*path;
    // RDKit❗✔️:           // newPath.push_back(otherIdx);
    // RDKit❗✔️:           // res.push_back(newPath);
    // RDKit❗✔️:           res.push_back(*path);
    // RDKit❗✔️:           res.rbegin()->push_back(otherIdx);
    // RDKit❗✔️:         } else if (allowRingClosures > 2 &&
    // RDKit❗✔️:                    static_cast<int>(path->size()) == allowRingClosures - 1) {
    // RDKit❗✔️:           // We *might* be adding the atom, but we need to make sure
    // RDKit❗✔️:           // that we're not just duplicating the second to last
    // RDKit❗✔️:           // element of the path:
    // RDKit❗✔️:           auto rIt = path->rbegin();
    // RDKit❗✔️:           rIt++;
    // RDKit❗✔️:           if (*rIt != static_cast<int>(otherIdx)) {
    // RDKit❗✔️:             // PATH_TYPE newPath=*path;
    // RDKit❗✔️:             // newPath.push_back(otherIdx);
    // RDKit❗✔️:             // res.push_back(newPath);
    // RDKit❗✔️:             res.push_back(*path);
    // RDKit❗✔️:             res.rbegin()->push_back(otherIdx);
    // RDKit❗✔️:           }
    // RDKit❗✔️:         }
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::Subgraphs::extendPaths
    let mut result = Vec::new();
    for path in paths {
        let end = *path.last().expect("path finder never stores an empty path");
        for other in 0..dimension {
            if ignore_atoms.is_some_and(|mask| mask[other]) {
                continue;
            }
            if !adjacency[end * dimension + other] {
                continue;
            }
            if distances.is_some_and(|matrix| matrix[path[0] * dimension + other] < path.len()) {
                continue;
            }
            if !path.contains(&other) {
                let mut extended = path.clone();
                extended.push(other);
                result.push(extended);
            } else if allow_ring_closures > 2 && path.len() == allow_ring_closures - 1 {
                let penultimate = path[path.len() - 2];
                if penultimate != other {
                    let mut extended = path.clone();
                    extended.push(other);
                    result.push(extended);
                }
            }
        }
    }
    result
}

fn bond_neighbor_map(
    graph: PathGraphAccess<'_>,
    use_hydrogens: bool,
) -> BTreeMap<usize, Vec<usize>> {
    // BEGIN RDKIT CPP FUNCTION getNbrsList
    // RDKit✔️✔️: for (int i = 0; i < nAtoms; i++) {
    // RDKit✔️✔️:   const Atom *atom = mol.getAtomWithIdx(i);
    // RDKit✔️✔️:   if (useHs || atom->getAtomicNum() != 1) {
    // RDKit✔️✔️:     while (bIt1 != end) {
    // RDKit✔️✔️:       const Bond *bond1 = mol[*bIt1];
    let mut result = BTreeMap::<usize, Vec<usize>>::new();
    for atom_index in 0..graph.atom_count() {
        if !use_hydrogens && graph.atomic_number(atom_index) == 1 {
            continue;
        }
        let atom_bonds = graph.neighbor_count(atom_index);
        for first in 0..atom_bonds {
            let (first_atom, first_bond) = graph.neighbor(atom_index, first);
            if !use_hydrogens && graph.atomic_number(first_atom) == 1 {
                continue;
            }
            result.entry(first_bond.index()).or_default();
            for second in 0..atom_bonds {
                let (second_atom, second_bond) = graph.neighbor(atom_index, second);
                if first_bond != second_bond
                    && (use_hydrogens || graph.atomic_number(second_atom) != 1)
                {
                    result
                        .get_mut(&first_bond.index())
                        .expect("entry was inserted")
                        .push(second_bond.index());
                }
            }
        }
    }
    // END RDKIT CPP FUNCTION getNbrsList
    result
}

fn all_subgraphs_of_length_from_neighbors(
    topology: &TopologyBlock,
    neighbors: &BTreeMap<usize, Vec<usize>>,
    target_length: usize,
    root: Option<AtomId>,
    ignore_atoms: Option<&[bool]>,
) -> Vec<Vec<BondId>> {
    // BEGIN RDKIT CPP FUNCTION RDKit::findAllSubgraphsOfLengthN (Release_2026_03_6)
    // RDKit❗✔️: PATH_LIST findAllSubgraphsOfLengthN(const ROMol &mol, unsigned int targetLen,
    // RDKit❗✔️:                                     bool useHs, int rootedAtAtom,
    // RDKit❗✔️:                                     boost::dynamic_bitset<> *ignoreAtoms) {
    // RDKit❗✔️:   /*********************************************
    // RDKit❗✔️:     FIX: Lots of issues here:
    // RDKit❗✔️:     - pathListType is defined as a container of "pathType", should it be a
    // RDKit❗✔️:   container
    // RDKit❗✔️:     of "pointers to pathtype"
    // RDKit❗✔️:     - to make few things clear it might be useful to typedef a
    // RDKit❗✔️:   "subgraphListType" even if it is exactly same as the "pathListType", just to
    // RDKit❗✔️:   not confuse between path vs. subgraph definitions
    // RDKit❗✔️:     - To make it consistent with the python version of this function in
    // RDKit❗✔️:   "subgraph.py"
    // RDKit❗✔️:     it return a "list of paths" instead of a "list of list of paths" (see
    // RDKit❗✔️:     "GetPathsUpTolength" in "molgraphs.cpp")
    // RDKit❗✔️:   ****************************************************************************/
    // RDKit❗✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == mol.getNumAtoms(),
    // RDKit❗✔️:                "bad ignoreAtoms size");
    // RDKit❗✔️:   boost::dynamic_bitset<> forbidden(mol.getNumBonds());
    // RDKit❗✔️:   // if there are any ignore atoms, mark any bonds involving them as forbidden
    // RDKit❗✔️:   if (ignoreAtoms) {
    // RDKit❗✔️:     for (const auto bond : mol.bonds()) {
    // RDKit❗✔️:       if (ignoreAtoms->test(bond->getBeginAtomIdx()) ||
    // RDKit❗✔️:           ignoreAtoms->test(bond->getEndAtomIdx())) {
    // RDKit❗✔️:         forbidden[bond->getIdx()] = 1;
    // RDKit❗✔️:       }
    // RDKit❗✔️:     }
    // RDKit❗✔️:   }
    // RDKit❗✔️:
    // RDKit❗✔️:   // this should be the only dependence on mol object:
    // RDKit❗✔️:   INT_INT_VECT_MAP nbrs;
    // RDKit❗✔️:   Subgraphs::getNbrsList(mol, useHs, nbrs);
    // RDKit❗✔️:
    // RDKit❗✔️:   // Start path at each bond
    // RDKit❗✔️:   PATH_LIST res;
    // RDKit❗✔️:
    // RDKit❗✔️:   // start paths at each bond:
    // RDKit❗✔️:   for (auto nbi = nbrs.begin(); nbi != nbrs.end(); ++nbi) {
    // RDKit❗✔️:     int i = (*nbi).first;
    // RDKit❗✔️:     if (forbidden[i]) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:     auto bi = mol.getBondWithIdx(i);
    // RDKit❗✔️:
    // RDKit❗✔️:     // if we're only returning paths rooted at a particular atom, check now
    // RDKit❗✔️:     // that this bond involves that atom:
    // RDKit❗✔️:     if (rootedAtAtom >= 0 &&
    // RDKit❗✔️:         bi->getBeginAtomIdx() != static_cast<unsigned int>(rootedAtAtom) &&
    // RDKit❗✔️:         bi->getEndAtomIdx() != static_cast<unsigned int>(rootedAtAtom)) {
    // RDKit❗✔️:       continue;
    // RDKit❗✔️:     }
    // RDKit❗✔️:
    // RDKit❗✔️:     // don't come back to this bond in the later subgraphs
    // RDKit❗✔️:     forbidden[i] = 1;
    // RDKit❗✔️:
    // RDKit❗✔️:     // start the recursive path building with the current bond
    // RDKit❗✔️:     PATH_TYPE spath;
    // RDKit❗✔️:     spath.clear();
    // RDKit❗✔️:     spath.push_back(i);
    // RDKit❗✔️:
    // RDKit❗✔️:     // neighbors of this bond are the next candidates
    // RDKit❗✔️:     INT_VECT cands = nbrs[i];
    // RDKit❗✔️:
    // RDKit❗✔️:     // now call the recursive function
    // RDKit❗✔️:     // little bit different from the python version
    // RDKit❗✔️:     // the result list of paths is passed as a reference, instead of on the fly
    // RDKit❗✔️:     // appending
    // RDKit❗✔️:     Subgraphs::recurseWalk(nbrs, spath, cands, targetLen, forbidden, res);
    // RDKit❗✔️:   }
    // RDKit❗✔️:   nbrs.clear();
    // RDKit❗✔️:   return res;
    // RDKit❗✔️: }
    // END RDKIT CPP FUNCTION RDKit::findAllSubgraphsOfLengthN
    let mut forbidden = ignored_bonds(PathGraphAccess::Concrete(topology), ignore_atoms);
    let mut raw_result = Vec::new();
    for (&start, adjacent) in neighbors {
        if forbidden[start] {
            continue;
        }
        if !root_allows_bond(topology, root, start) {
            continue;
        }
        forbidden[start] = true;
        recurse_walk(
            neighbors,
            vec![start],
            adjacent.clone(),
            target_length,
            forbidden.clone(),
            &mut raw_result,
        );
    }
    raw_result
        .into_iter()
        .map(|row| row.into_iter().map(BondId::new).collect())
        .collect()
}

fn recurse_walk(
    neighbors: &BTreeMap<usize, Vec<usize>>,
    path: Vec<usize>,
    mut candidates: Vec<usize>,
    target_length: usize,
    mut forbidden: Vec<bool>,
    result: &mut Vec<Vec<usize>>,
) {
    // BEGIN RDKIT CPP FUNCTION recurseWalk
    // RDKit✔️✔️: if (spath.size() == targetLen) { res.push_back(spath); return; }
    // RDKit✔️✔️: if (spath.size() > targetLen) { return; }
    if path.len() == target_length {
        result.push(path);
        return;
    }
    if path.len() > target_length {
        return;
    }
    // RDKit✔️✔️: while (cands.size() != 0) {
    // RDKit✔️✔️:   int next = cands.back(); cands.pop_back();
    while let Some(next) = candidates.pop() {
        if forbidden[next] {
            continue;
        }
        forbidden[next] = true;
        let mut stack = candidates.clone();
        if let Some(next_neighbors) = neighbors.get(&next) {
            for &bond in next_neighbors {
                if !forbidden[bond] {
                    stack.push(bond);
                }
            }
        }
        let mut next_path = path.clone();
        next_path.push(next);
        recurse_walk(
            neighbors,
            next_path,
            stack,
            target_length,
            forbidden.clone(),
            result,
        );
    }
    // END RDKIT CPP FUNCTION recurseWalk
}

fn recurse_walk_range(
    neighbors: &BTreeMap<usize, Vec<usize>>,
    path: Vec<usize>,
    mut candidates: Vec<usize>,
    lower_length: usize,
    upper_length: usize,
    mut forbidden: Vec<bool>,
    result: &mut BTreeMap<usize, Vec<Vec<BondId>>>,
) {
    // BEGIN RDKIT CPP FUNCTION recurseWalkRange
    // RDKit✔️✔️: unsigned int nsize = spath.size();
    // RDKit✔️✔️: if ((nsize >= lowerLen) && (nsize <= upperLen)) {
    // RDKit✔️✔️:   res[nsize].push_back(spath);
    // RDKit✔️✔️: }
    let length = path.len();
    if length >= lower_length && length <= upper_length {
        result
            .entry(length)
            .or_default()
            .push(path.iter().copied().map(BondId::new).collect());
    }
    if length >= upper_length {
        return;
    }
    while let Some(next) = candidates.pop() {
        if forbidden[next] {
            continue;
        }
        forbidden[next] = true;
        let mut stack = candidates.clone();
        if let Some(next_neighbors) = neighbors.get(&next) {
            for &bond in next_neighbors {
                if !forbidden[bond] {
                    stack.push(bond);
                }
            }
        }
        let mut next_path = path.clone();
        next_path.push(next);
        recurse_walk_range(
            neighbors,
            next_path,
            stack,
            lower_length,
            upper_length,
            forbidden.clone(),
            result,
        );
    }
    // END RDKIT CPP FUNCTION recurseWalkRange
}

fn root_allows_bond(topology: &TopologyBlock, root: Option<AtomId>, bond: usize) -> bool {
    graph_root_allows_bond(PathGraphAccess::Concrete(topology), root, bond)
}
fn graph_root_allows_bond(graph: PathGraphAccess<'_>, root: Option<AtomId>, bond: usize) -> bool {
    let Some(root) = root else {
        return true;
    };
    let bond = graph.bond(bond);
    bond.begin() == root || bond.end() == root
}
fn bond_between(topology: &TopologyBlock, begin: usize, end: usize) -> Option<BondId> {
    graph_bond_between(PathGraphAccess::Concrete(topology), begin, end)
}
fn graph_bond_between(graph: PathGraphAccess<'_>, begin: usize, end: usize) -> Option<BondId> {
    (0..graph.neighbor_count(begin)).find_map(|position| {
        let (other, bond) = graph.neighbor(begin, position);
        (other == end).then_some(bond)
    })
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PathDiscriminator(u32, usize, usize);

fn path_discriminator(
    topology: &TopologyBlock,
    path: &[BondId],
    use_bond_orders: bool,
    extra_atom_invariants: Option<&[u32]>,
) -> Result<PathDiscriminator, PathError> {
    // BEGIN RDKIT CPP FUNCTION calcPathDiscriminators
    // RDKit✔️✔️: std::vector<int32_t> atomsUsed(mol.getNumAtoms(), -1);
    // RDKit✔️✔️: std::vector<const Atom *> atoms;
    // RDKit✔️✔️: std::vector<uint32_t> pathDegrees;
    let mut atom_positions = vec![None; topology.atoms.len()];
    let mut atoms = Vec::new();
    let mut degrees = Vec::<u32>::new();
    for bond_id in path {
        let bond = &topology.bonds[bond_id.index()];
        for atom_id in [bond.begin(), bond.end()] {
            match atom_positions[atom_id.index()] {
                Some(position) => degrees[position] += 1,
                None => {
                    atom_positions[atom_id.index()] = Some(atoms.len());
                    atoms.push(atom_id);
                    degrees.push(1);
                }
            }
        }
    }
    // RDKit✔️✔️: uint32_t invar = atom->getAtomicNum();
    // RDKit✔️✔️: hash_combine(invar, pathDegrees[i]);
    // RDKit✔️✔️: hash_combine(invar, atom->getFormalCharge());
    // RDKit✔️✔️: int deltaMass = static_cast<int>(
    // RDKit✔️✔️:     atom->getMass() - PeriodicTable::getTable()->getAtomicWeight(atom->getAtomicNum()));
    // RDKit✔️✔️: hash_combine(invar, deltaMass);
    let mut invariants = Vec::with_capacity(atoms.len());
    for (position, atom_id) in atoms.iter().copied().enumerate() {
        let atom = &topology.atoms[atom_id.index()];
        let mut invariant = u32::from(atom.atomic_number());
        hash_combine(&mut invariant, degrees[position]);
        hash_combine(&mut invariant, i32::from(atom.formal_charge()) as u32);
        let isotope_mass =
            atomic_mass(atom.element(), atom.isotope()).map_err(PathError::PeriodicTable)?;
        let average_mass = atomic_mass(atom.element(), None).map_err(PathError::PeriodicTable)?;
        hash_combine(&mut invariant, (isotope_mass - average_mass) as i32 as u32);
        if atom.is_aromatic() {
            hash_combine(&mut invariant, 1);
        }
        if let Some(extra) = extra_atom_invariants {
            hash_combine(&mut invariant, extra[atom_id.index()]);
        }
        invariants.push(invariant);
    }
    // RDKit✔️✔️: unsigned int nCycles = path.size() / 2 + 1;
    // RDKit✔️✔️: for (unsigned int cycle = 0; cycle < nCycles; ++cycle) {
    for _ in 0..path.len() / 2 + 1 {
        let mut local = vec![Vec::<u32>::new(); atoms.len()];
        for bond_id in path {
            let bond = &topology.bonds[bond_id.index()];
            let begin = atom_positions[bond.begin().index()].expect("path atom was recorded");
            let end = atom_positions[bond.end().index()].expect("path atom was recorded");
            let mut begin_value = invariants[begin];
            let mut end_value = invariants[end];
            if use_bond_orders {
                hash_combine(&mut begin_value, bond.order().rdkit_code() as u32);
                hash_combine(&mut end_value, bond.order().rdkit_code() as u32);
            }
            local[begin].push(end_value);
            local[end].push(begin_value);
        }
        for (invariant, neighbors) in invariants.iter_mut().zip(&mut local) {
            neighbors.sort_unstable();
            *invariant = hash_range(neighbors);
        }
    }
    invariants.sort_unstable();
    Ok(PathDiscriminator(
        hash_range(&invariants),
        path.len(),
        atoms.len(),
    ))
    // END RDKIT CPP FUNCTION calcPathDiscriminators
}

fn hash_combine(seed: &mut u32, value: u32) {
    // RDKit✔️✔️: seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    *seed ^= value
        .wrapping_add(0x9e37_79b9)
        .wrapping_add((*seed).wrapping_shl(6))
        .wrapping_add(*seed >> 2);
}

fn hash_range(values: &[u32]) -> u32 {
    // RDKit✔️✔️: std::hash_result_t seed = 0;
    // RDKit✔️✔️: for (; first != last; ++first) { hash_combine(seed, *first); }
    let mut seed = 0;
    for &value in values {
        hash_combine(&mut seed, value);
    }
    seed
}

fn remap_selected_substance_groups(
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Vec<SubstanceGroup> {
    let mut selected = topology
        .substance_groups
        .iter()
        .map(|group| {
            group.can_remap_without_parent(&mapping.atoms.old_to_new, &mapping.bonds.old_to_new)
        })
        .collect::<Vec<_>>();
    loop {
        let mut changed = false;
        for (index, group) in topology.substance_groups.iter().enumerate() {
            if selected[index]
                && group
                    .parent()
                    .is_some_and(|parent| !selected.get(parent.index()).copied().unwrap_or(false))
            {
                selected[index] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    let mut sgroup_map = vec![None; topology.substance_groups.len()];
    let mut next = 0;
    for (old, keep) in selected.iter().copied().enumerate() {
        if keep {
            sgroup_map[old] = Some(SubstanceGroupId::new(next));
            next += 1;
        }
    }
    topology
        .substance_groups
        .iter()
        .enumerate()
        .filter(|(index, _)| selected[*index])
        .filter_map(|(index, group)| {
            group.remapped(
                sgroup_map[index].expect("selected SGroup has a new id"),
                &mapping.atoms.old_to_new,
                &mapping.bonds.old_to_new,
                &sgroup_map,
            )
        })
        .collect()
}

fn remap_selected_stereo_groups(
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Result<Vec<StereoGroup>, cosmolkit_model::StereoGroupError> {
    topology
        .stereo_groups
        .iter()
        .map(|group| {
            // RDKit✔️✔️: return objects.empty() ||
            // RDKit✔️✔️:        std::any_of(objects.begin(), objects.end(), [&](auto &object) {
            // RDKit✔️✔️:          return selected_indices[object->getIdx()];
            // RDKit✔️✔️:        });
            let atom_side_selected = group.atoms().is_empty()
                || group
                    .atoms()
                    .iter()
                    .any(|atom| mapping.atoms.old_to_new[atom.index()].is_some());
            let bond_side_selected = group.bonds().is_empty()
                || group
                    .bonds()
                    .iter()
                    .any(|bond| mapping.bonds.old_to_new[bond.index()].is_some());
            if !atom_side_selected || !bond_side_selected {
                return Ok(None);
            }
            let atoms = group
                .atoms()
                .iter()
                .filter_map(|atom| mapping.atoms.old_to_new[atom.index()])
                .collect();
            let bonds = group
                .bonds()
                .iter()
                .filter_map(|bond| mapping.bonds.old_to_new[bond.index()])
                .collect();
            let remapped = StereoGroup::new(group.kind(), atoms, bonds)?;
            let remapped = match group.id() {
                Some(id) => remapped.with_id(id),
                None => remapped,
            };
            // RDKit✔️✔️: extracted_stereo_groups.back().setWriteId(stereo_group.getWriteId());
            Ok(Some(remapped.with_write_id(group.write_id())))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(|rows| rows.into_iter().flatten().collect())
}

#[cfg(test)]
mod cf3d_frag_f01_tests {
    use super::{ConnectedComponents, connected_components};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, TopologyBlock,
    };
    use cosmolkit_types::Element;

    fn topology(atom_count: usize, edges: &[(usize, usize)]) -> TopologyBlock {
        let atoms = (0..atom_count)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = edges
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
            .expect("fixed component topology is valid")
    }

    #[test]
    fn cf3d_frag_f01_empty_topology_has_no_labels_or_components() {
        assert_eq!(
            connected_components(&TopologyBlock::default()).unwrap(),
            ConnectedComponents {
                atom_to_component: Vec::new(),
                components: Vec::new(),
            }
        );
    }

    #[test]
    fn cf3d_frag_f01_single_atom_is_one_component() {
        assert_eq!(
            connected_components(&topology(1, &[])).unwrap(),
            ConnectedComponents {
                atom_to_component: vec![0],
                components: vec![vec![AtomId::new(0)]],
            }
        );
    }

    #[test]
    fn cf3d_frag_f01_connected_cycle_is_one_source_ordered_component() {
        let cycle = topology(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]);
        assert_eq!(
            connected_components(&cycle).unwrap(),
            ConnectedComponents {
                atom_to_component: vec![0, 0, 0, 0],
                components: vec![vec![
                    AtomId::new(0),
                    AtomId::new(1),
                    AtomId::new(2),
                    AtomId::new(3),
                ]],
            }
        );
    }

    #[test]
    fn cf3d_frag_f01_interleaved_disconnected_rows_keep_source_membership() {
        let disconnected = topology(5, &[(0, 2), (1, 3)]);
        assert_eq!(
            connected_components(&disconnected).unwrap(),
            ConnectedComponents {
                atom_to_component: vec![0, 1, 0, 1, 2],
                components: vec![
                    vec![AtomId::new(0), AtomId::new(2)],
                    vec![AtomId::new(1), AtomId::new(3)],
                    vec![AtomId::new(4)],
                ],
            }
        );
    }

    #[test]
    fn cf3d_frag_f01_components_follow_ascending_first_source_row() {
        let disconnected = topology(7, &[(5, 6), (0, 3), (2, 4)]);
        assert_eq!(
            connected_components(&disconnected).unwrap(),
            ConnectedComponents {
                atom_to_component: vec![0, 1, 2, 0, 2, 3, 3],
                components: vec![
                    vec![AtomId::new(0), AtomId::new(3)],
                    vec![AtomId::new(1)],
                    vec![AtomId::new(2), AtomId::new(4)],
                    vec![AtomId::new(5), AtomId::new(6)],
                ],
            }
        );
    }
}

#[cfg(test)]
mod cf3d_sgids_core_2_tests {
    use super::{DetachedPathSubgraph, SubtopologyParams, subtopology_from_path};
    use cosmolkit_model::{
        Atom, AtomId, AtomSpec, Bond, BondId, BondOrder, BondSpec, StereoGroup, StereoGroupKind,
        TopologyBlock,
    };
    use cosmolkit_types::Element;

    #[test]
    fn cf3d_sgids_core_2_subset_path_preserves_ids_and_source_group_order() {
        let atoms = (0..4)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect();
        let bonds = [(0, 1), (2, 3)]
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect();
        let stereo_groups = vec![
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(0), AtomId::new(2)],
                vec![BondId::new(0), BondId::new(1)],
            )
            .expect("valid distinct stereo members")
            .with_id(17)
            .with_write_id(9),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(2)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members")
            .with_id(41)
            .with_write_id(5),
            StereoGroup::new(
                StereoGroupKind::And,
                vec![AtomId::new(1), AtomId::new(3)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members")
            .with_id(0),
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(0), AtomId::new(2)],
                vec![],
            )
            .expect("valid distinct stereo members")
            .with_write_id(12),
        ];
        let source_groups = stereo_groups.clone();
        let topology = TopologyBlock::try_from_parts(atoms, bonds, vec![], stereo_groups)
            .expect("fixed source topology and group references are valid");

        let result =
            subtopology_from_path(&topology, &[BondId::new(0)], &SubtopologyParams::default())
                .expect("selected source bond forms a valid subtopology");
        let DetachedPathSubgraph::Concrete(subgraph) = result.subgraph else {
            panic!("default subtopology result is concrete");
        };
        assert_eq!(
            result.mapping.atoms.old_to_new,
            vec![Some(AtomId::new(0)), Some(AtomId::new(1)), None, None]
        );
        assert_eq!(
            result.mapping.bonds.old_to_new,
            vec![Some(BondId::new(0)), None]
        );

        let groups = &subgraph.stereo_groups;
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].kind(), StereoGroupKind::Or);
        assert_eq!(groups[0].id(), Some(17));
        assert_eq!(groups[0].write_id(), 9);
        assert_eq!(groups[0].atoms(), &[AtomId::new(0)]);
        assert_eq!(groups[0].bonds(), &[BondId::new(0)]);

        assert_eq!(groups[1].kind(), StereoGroupKind::And);
        assert_eq!(groups[1].id(), Some(0));
        assert_eq!(groups[1].write_id(), 0);
        assert_eq!(groups[1].atoms(), &[AtomId::new(1)]);
        assert_eq!(groups[1].bonds(), &[BondId::new(0)]);

        assert_eq!(groups[2].kind(), StereoGroupKind::Absolute);
        assert_eq!(groups[2].id(), None);
        assert_eq!(groups[2].write_id(), 12);
        assert_eq!(groups[2].atoms(), &[AtomId::new(0)]);
        assert!(groups[2].bonds().is_empty());
        assert_eq!(topology.stereo_groups, source_groups);
    }
}

fn validate_ignored_atom_mask(mask: Option<&[bool]>, expected: usize) -> Result<(), PathError> {
    // RDKit✔️✔️:   PRECONDITION(!ignoreAtoms || ignoreAtoms->size() == mol.getNumAtoms(),
    // RDKit✔️✔️:                "bad ignoreAtoms size");
    if let Some(mask) = mask
        && mask.len() != expected
    {
        return Err(PathError::IgnoredAtomMaskLength {
            actual: mask.len(),
            expected,
        });
    }
    Ok(())
}

fn ignored_bonds(graph: PathGraphAccess<'_>, mask: Option<&[bool]>) -> Vec<bool> {
    // RDKit✔️❌:   // if there are any ignore atoms, mark any bonds involving them as forbidden
    // RDKit✔️❌:   if (ignoreAtoms) {
    // RDKit✔️❌:     for (const auto bond : mol.bonds()) {
    // RDKit✔️❌:       if (ignoreAtoms->test(bond->getBeginAtomIdx()) ||
    // RDKit✔️❌:           ignoreAtoms->test(bond->getEndAtomIdx())) {
    // RDKit✔️❌:         forbidden[bond->getIdx()] = 1;
    // RDKit✔️❌:       }
    // RDKit✔️❌:     }
    // RDKit✔️❌:   }

    // Source loops exactly one molecule bond pass when ignoreAtoms is present.
    // Bitmap has one byte per entry rather than packed bits: same O(B) work,
    // larger storage is retained and must not claim packed-memory parity.
    let mut forbidden = vec![false; graph.bond_count()];
    if let Some(mask) = mask {
        for index in 0..graph.bond_count() {
            let bond = graph.bond(index);
            if mask[bond.begin().index()] || mask[bond.end().index()] {
                forbidden[index] = true;
            }
        }
    }
    forbidden
}
