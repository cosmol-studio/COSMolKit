//! SMARTS/search algorithms over detached query and topology values.

mod cx_lowering;
mod generic_groups;
mod matcher;
mod query_behavior;
mod query_graph_behavior;
mod smarts_parse;
mod smarts_write;
mod target;

pub use cx_lowering::{CxQueryLoweringError, apply_cx_to_query_graph};
pub use matcher::{
    AtomCoordsMatchFunctor, ExtraAtomCheck, ExtraBondCheck, ExtraFinalCheck, QueryInput,
    SubstructMatchError, SubstructMatchOverload, SubstructMatchParams,
    SubstructMatchParamsJsonError, SubstructMatchResult, check_substruct_match_overload_support,
    get_substruct_match, get_substruct_matches, get_substruct_matches_with_params,
    has_substruct_match, substruct_match_params_to_json, try_get_substruct_matches_with_params,
    try_get_substruct_matches_with_params_and_context, update_substruct_match_params_from_json,
};
pub use query_behavior::{
    QUERY_SCAN_MAGIC_VALUE, QueryConstructionError, QueryMatchContext, SmartsParseError,
    atom_matches_query, atom_matches_query_with_context, atom_predicate_matches,
    atom_predicate_matches_with_context, atom_query_has_magic_value, bond_matches_query,
    bond_matches_query_with_context, bond_predicate_matches, bond_predicate_matches_with_context,
    build_query_match_context, convert_complex_name_to_query, is_atom_aromatic,
    make_single_or_aromatic_bond_query, query_bond_min_ring_size, query_is_bond_in_ring,
};
pub use smarts_parse::{SmartsParseParams, compile_query_fixture, parse_smarts};
pub use smarts_write::{
    SmartsWriteError, SmartsWriteParams, atom_to_smarts, bond_to_smarts, query_atom_to_smarts,
    query_bond_to_smarts, query_graph_fragment_to_cx_smarts, query_graph_fragment_to_smarts,
    query_graph_to_cx_smarts, query_graph_to_smarts, write_smarts,
};
pub use target::{SearchTarget, SearchTargetAccess};

use cosmolkit_model::{AtomId, BondId, CoordinateBlock, TopologyBlock};
pub use cosmolkit_model::{
    AtomQueryPredicate, BondQueryPredicate, QueryAtom, QueryBond, QueryGraph, QueryGraphError,
    QueryNode,
};
pub use cosmolkit_types::BondDirection;

/// Errors raised while compiling a query graph into a reusable match plan.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryCompileError {
    #[error("query graph is invalid: {0}")]
    InvalidGraph(String),
}

/// A reusable, detached query execution plan.
///
/// The plan owns the query value and a deterministic atom visitation order.
/// It contains no matcher state, live molecule, runtime cache, or operation
/// capability.  Repeated target matching can therefore share one immutable
/// plan without coupling this crate to the `cosmolkit` runtime.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledQuery {
    query: QueryGraph,
    atom_order: Vec<usize>,
    compiled_graph: matcher::CompiledQueryGraph,
}

impl CompiledQuery {
    /// Compile a query graph into a deterministic plan.
    pub fn compile(query: QueryGraph) -> Result<Self, QueryCompileError> {
        query
            .validate()
            .map_err(|error| QueryCompileError::InvalidGraph(error.to_string()))?;

        let compiled_graph = matcher::compile_query_graph(&query);
        let atom_order = matcher::compile_query_order_from_graph(&compiled_graph);

        Ok(Self {
            query,
            atom_order,
            compiled_graph,
        })
    }

    #[must_use]
    pub fn query(&self) -> &QueryGraph {
        &self.query
    }

    #[must_use]
    pub fn atom_order(&self) -> &[usize] {
        &self.atom_order
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.query.num_atoms()
    }

    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.query.num_bonds()
    }

    /// Match this immutable plan against a detached topology without
    /// rebuilding the query-side order or adjacency view.
    pub fn matches(&self, topology: &TopologyBlock) -> Result<Vec<MatchResult>, MatchError> {
        topology
            .validate()
            .map_err(|error| MatchError::InvalidTarget(error.to_string()))?;
        let coordinates = CoordinateBlock::default();
        let target = SearchTarget::new(topology, &coordinates, &topology.stereo_groups, None, None);
        self.matches_target(&target)
    }

    /// Match against the complete detached target view, including conformers
    /// and optional precomputed ring/valence assignments.
    pub fn matches_target(
        &self,
        target: &SearchTarget<'_>,
    ) -> Result<Vec<MatchResult>, MatchError> {
        matcher::get_substruct_matches_with_compiled_query(
            target,
            &self.query,
            &SubstructMatchParams::default(),
            &self.compiled_graph,
            &self.atom_order,
        )
        .map_err(MatchError::from)
    }
}

/// Compile a detached query graph into a reusable execution plan.
pub fn compile_query(query: &QueryGraph) -> Result<CompiledQuery, QueryCompileError> {
    CompiledQuery::compile(query.clone())
}

/// Behaviour facade for an existing query value.
///
/// Construction remains a module-level parser concern; this operator only
/// exposes behaviour that interprets an already-owned `QueryGraph`.
#[derive(Debug, Clone, Copy)]
pub struct QueryGraphOperator<'a> {
    inner: &'a QueryGraph,
}

impl<'a> QueryGraphOperator<'a> {
    #[must_use]
    pub const fn new(inner: &'a QueryGraph) -> Self {
        Self { inner }
    }

    #[must_use]
    pub const fn inner(self) -> &'a QueryGraph {
        self.inner
    }

    pub fn compile(self) -> Result<CompiledQuery, QueryCompileError> {
        compile_query(self.inner)
    }

    pub fn matches(self, topology: &TopologyBlock) -> Result<Vec<MatchResult>, MatchError> {
        match_query(self.inner, topology)
    }

    pub fn to_smarts(self, params: &SmartsWriteParams) -> Result<String, SmartsWriteError> {
        write_smarts(self.inner, params)
    }

    pub fn atom_to_smarts(
        self,
        atom_id: AtomId,
        params: &SmartsWriteParams,
    ) -> Result<String, SmartsWriteError> {
        atom_to_smarts(self.inner, atom_id, params)
    }

    pub fn bond_to_smarts(self, bond_id: BondId) -> Result<String, SmartsWriteError> {
        bond_to_smarts(self.inner, bond_id)
    }
}

/// A query-to-target mapping produced by a matcher.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchResult {
    pub atom_mapping: Vec<usize>,
    pub bond_mapping: Vec<usize>,
}

impl MatchResult {
    #[must_use]
    pub fn new(atom_mapping: Vec<usize>, bond_mapping: Vec<usize>) -> Self {
        Self {
            atom_mapping,
            bond_mapping,
        }
    }
}

/// Result metadata for a maximum-common-substructure search.
#[derive(Debug, Clone, PartialEq)]
pub struct McsResult {
    pub query: QueryGraph,
    pub atom_count: usize,
    pub bond_count: usize,
    pub completed: bool,
}

impl McsResult {
    #[must_use]
    pub fn new(query: QueryGraph, atom_count: usize, bond_count: usize, completed: bool) -> Self {
        Self {
            query,
            atom_count,
            bond_count,
            completed,
        }
    }
}

/// Errors raised while matching a detached query against a detached topology.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MatchError {
    #[error("query graph is invalid: {0}")]
    InvalidQuery(String),
    #[error("target topology is invalid: {0}")]
    InvalidTarget(String),
    #[error("unsupported atom query predicate: {0}")]
    UnsupportedAtomPredicate(&'static str),
    #[error("unsupported bond query predicate: {0}")]
    UnsupportedBondPredicate(&'static str),
    #[error(transparent)]
    Substruct(#[from] SubstructMatchError),
}

/// Match a detached query against a detached topology.
pub fn match_query(
    query: &QueryGraph,
    topology: &TopologyBlock,
) -> Result<Vec<MatchResult>, MatchError> {
    query
        .validate()
        .map_err(|error| MatchError::InvalidQuery(error.to_string()))?;
    let plan = compile_query(query).map_err(|error| MatchError::InvalidQuery(error.to_string()))?;
    plan.matches(topology)
}

/// Match a detached query against a complete detached target view.
pub fn match_query_target(
    query: &QueryGraph,
    target: &SearchTarget<'_>,
) -> Result<Vec<MatchResult>, MatchError> {
    query
        .validate()
        .map_err(|error| MatchError::InvalidQuery(error.to_string()))?;
    let plan = compile_query(query).map_err(|error| MatchError::InvalidQuery(error.to_string()))?;
    plan.matches_target(target)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cosmolkit_model::{AdjacencyList, Atom, AtomId, AtomSpec, Bond, BondId, BondSpec};
    use cosmolkit_types::{BondOrder, Element};

    fn query() -> QueryGraph {
        let atoms = vec![
            cosmolkit_model::QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C)),
            cosmolkit_model::QueryAtom::new(AtomId::new(1), AtomSpec::new(Element::C)),
            cosmolkit_model::QueryAtom::new(AtomId::new(2), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![
            cosmolkit_model::QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            cosmolkit_model::QueryBond::new(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Single),
            ),
        ];
        QueryGraph::from_parts(
            atoms,
            bonds,
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap()
    }

    #[test]
    fn compile_query_orders_constrained_atoms_first() {
        let plan = compile_query(&query()).unwrap();
        assert_eq!(plan.atom_order(), &[1, 0, 2]);
        assert_eq!(plan.num_atoms(), 3);
        assert_eq!(plan.num_bonds(), 2);
    }

    #[test]
    fn match_query_uses_detached_topology_and_returns_index_mappings() {
        let query = query();
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(2), AtomSpec::new(Element::O)),
            Atom::from_spec(AtomId::new(3), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(
                    AtomId::new(0),
                    AtomId::new(1),
                    cosmolkit_types::BondOrder::Single,
                ),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(
                    AtomId::new(1),
                    AtomId::new(2),
                    cosmolkit_types::BondOrder::Single,
                ),
            ),
            Bond::from_spec(
                BondId::new(2),
                BondSpec::new(
                    AtomId::new(1),
                    AtomId::new(3),
                    cosmolkit_types::BondOrder::Single,
                ),
            ),
        ];
        let topology = TopologyBlock {
            atoms,
            bonds: bonds.clone(),
            adjacency: AdjacencyList::from_topology(4, &bonds),
            ..TopologyBlock::default()
        };
        let matches = match_query(&query, &topology).unwrap();
        assert_eq!(matches.len(), 2);
        assert!(matches.iter().any(|value| value.atom_mapping == [0, 1, 2]));
        assert!(matches.iter().any(|value| value.atom_mapping == [3, 1, 2]));
        assert!(matches.iter().all(|value| value.bond_mapping[1] == 1));
    }

    #[test]
    fn compiled_query_reuses_its_planned_order() {
        let query = query();
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::O)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let topology = TopologyBlock {
            atoms,
            bonds: bonds.clone(),
            adjacency: AdjacencyList::from_topology(2, &bonds),
            ..TopologyBlock::default()
        };
        let plan = compile_query(&query).expect("compile");
        assert!(plan.matches(&topology).expect("match").is_empty());
        assert_eq!(plan.atom_order(), &[1, 0, 2]);
    }

    #[test]
    fn unsupported_predicates_fail_closed() {
        let atoms = vec![cosmolkit_model::QueryAtom::from_parts(
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            QueryNode::predicate(AtomQueryPredicate::UnsupportedFeature("recursive SMARTS")),
        )];
        let query = QueryGraph::from_parts(
            atoms,
            Vec::new(),
            Default::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .unwrap();
        let topology = TopologyBlock {
            atoms: vec![Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C))],
            adjacency: AdjacencyList::from_topology(1, &[]),
            ..TopologyBlock::default()
        };
        assert!(matches!(
            match_query(&query, &topology),
            Err(MatchError::Substruct(SubstructMatchError::Unsupported {
                branch: "recursive SMARTS",
                rdkit_function: "QueryAtom::Match",
            }))
        ));
    }

    fn cyclopropane_topology() -> TopologyBlock {
        let atoms = (0..3)
            .map(|index| Atom::from_spec(AtomId::new(index), AtomSpec::new(Element::C)))
            .collect::<Vec<_>>();
        let bonds = [(0, 1), (1, 2), (2, 0)]
            .into_iter()
            .enumerate()
            .map(|(index, (begin, end))| {
                Bond::from_spec(
                    BondId::new(index),
                    BondSpec::new(AtomId::new(begin), AtomId::new(end), BondOrder::Single),
                )
            })
            .collect::<Vec<_>>();
        TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        }
    }

    #[test]
    fn complete_writer_round_trips_ring_queries() {
        let query = parse_smarts("C1CC1", &SmartsParseParams::default()).expect("parse ring");
        let written = write_smarts(&query, &SmartsWriteParams::default()).expect("write ring");
        let reparsed = parse_smarts(&written, &SmartsParseParams::default()).expect("reparse ring");
        assert_eq!(reparsed.num_atoms(), 3);
        assert_eq!(reparsed.num_bonds(), 3);
    }

    #[test]
    fn complete_matcher_handles_ring_and_range_predicates() {
        let topology = cyclopropane_topology();
        for smarts in ["[R]", "[r3]", "[D{2-3}]"] {
            let query = parse_smarts(smarts, &SmartsParseParams::default()).expect("parse query");
            let matches = match_query(&query, &topology).expect("match query");
            assert_eq!(matches.len(), 3, "{smarts}");
        }
    }

    #[test]
    fn complete_matcher_applies_generic_group_labels() {
        let mut query = parse_smarts("C*", &SmartsParseParams::default()).expect("parse query");
        query.atoms_mut()[1]
            .atom_mut()
            .set_prop("_QueryAtomGenericLabel", "ALK");
        let params = SubstructMatchParams {
            use_generic_matchers: true,
            ..SubstructMatchParams::default()
        };

        let target = |second: Element| {
            let atoms = vec![
                Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
                Atom::from_spec(AtomId::new(1), AtomSpec::new(second)),
            ];
            let bonds = vec![Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )];
            TopologyBlock {
                adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
                atoms,
                bonds,
                ..TopologyBlock::default()
            }
        };

        let carbon = target(Element::C);
        let oxygen = target(Element::O);
        let coordinates = CoordinateBlock::default();
        let carbon_target =
            SearchTarget::new(&carbon, &coordinates, &carbon.stereo_groups, None, None);
        let oxygen_target =
            SearchTarget::new(&oxygen, &coordinates, &oxygen.stereo_groups, None, None);
        assert_eq!(
            get_substruct_matches_with_params(&carbon_target, &query, &params).len(),
            1
        );
        assert!(get_substruct_matches_with_params(&oxygen_target, &query, &params).is_empty());
    }

    #[test]
    fn complete_matcher_handles_recursive_smarts() {
        let atoms = vec![
            Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C)),
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::C)),
        ];
        let bonds = vec![Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )];
        let topology = TopologyBlock {
            adjacency: AdjacencyList::from_topology(atoms.len(), &bonds),
            atoms,
            bonds,
            ..TopologyBlock::default()
        };
        let query =
            parse_smarts("[$(C-C)]", &SmartsParseParams::default()).expect("parse recursive");
        assert_eq!(match_query(&query, &topology).expect("match").len(), 2);
    }
}
