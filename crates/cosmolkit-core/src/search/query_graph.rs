//! First-class SMARTS query graph model.
//!
//! A query is a graph with query predicates attached to atoms and bonds.  It
//! is intentionally distinct from a concrete [`Molecule`]: query predicates
//! describe matching semantics and are not chemical state. Matching consumes
//! this graph directly; no query-to-molecule projection is part of the
//! public query contract.

use super::query::{AtomQueryPredicate, BondQueryPredicate, QueryNode};
use crate::{Atom, Bond, Conformer2D, Conformer3D, Molecule, StereoGroup};
use std::collections::BTreeMap;

/// Operations that interpret a query graph.
///
/// `QueryGraph` is the canonical query value.  Parsing, rendering, matching,
/// compilation, and fingerprint integration are search-domain behaviour and
/// are intentionally grouped here instead of being implemented by the value
/// type itself.  The operator borrows the graph and therefore cannot replace
/// or mutate the query value while an operation is running.
pub struct QueryGraphOperator<'a> {
    inner: &'a QueryGraph,
}

impl std::fmt::Debug for QueryGraphOperator<'_> {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter
            .debug_struct("QueryGraphOperator")
            .field("inner", &self.inner)
            .finish()
    }
}

impl<'a> QueryGraphOperator<'a> {
    #[must_use]
    pub const fn new(inner: &'a QueryGraph) -> Self {
        Self { inner }
    }

    #[must_use]
    pub const fn inner(&self) -> &'a QueryGraph {
        self.inner
    }

    pub fn to_smarts(&self, params: &crate::SmilesWriteParams) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_to_smarts(self.inner, params)
    }

    pub fn to_cx_smarts(&self, params: &crate::SmilesWriteParams) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_to_cx_smarts(self.inner, params)
    }

    pub fn fragment_to_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_fragment_to_smarts(self.inner, params, atom_ids, bond_ids)
    }

    pub fn fragment_to_cx_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        crate::search::smarts_write::query_graph_fragment_to_cx_smarts(self.inner, params, atom_ids, bond_ids)
    }

    pub fn pattern_fingerprint(
        &self,
        params: &crate::PatternFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        let molecule = self
            .inner
            .to_molecule()
            .map_err(|error| crate::FingerprintError::Pattern {
                reason: error.to_string(),
            })?;
        crate::fingerprint::pattern_fingerprint(&molecule, params)
    }

    #[must_use]
    pub fn compile(&self) -> CompiledQuery {
        CompiledQuery::compile(self.inner.clone())
    }

    pub fn matches(&self, molecule: &Molecule) -> Vec<MatchResult> {
        crate::search::substruct::get_substruct_matches(molecule, self.inner)
    }
}

/// An atom in a [`QueryGraph`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryAtom {
    atom: Atom,
    predicate: QueryNode<AtomQueryPredicate>,
}

impl QueryAtom {
    pub(crate) fn new(atom: Atom, predicate: QueryNode<AtomQueryPredicate>) -> Self {
        Self { atom, predicate }
    }

    /// Concrete atom attributes used by predicates and match callbacks.
    #[must_use]
    pub fn atom(&self) -> &Atom {
        &self.atom
    }

    /// The complete recursive atom predicate tree.
    #[must_use]
    pub fn predicate(&self) -> &QueryNode<AtomQueryPredicate> {
        &self.predicate
    }

    /// Pattern atom index.
    #[must_use]
    pub fn index(&self) -> usize {
        self.atom.id().index()
    }

    #[must_use]
    pub fn id(&self) -> crate::AtomId {
        self.atom.id()
    }

    #[must_use]
    pub fn atom_map(&self) -> Option<u32> {
        self.atom.atom_map()
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.atom.prop(key)
    }

    #[must_use]
    pub fn query(&self) -> Option<&QueryNode<AtomQueryPredicate>> {
        Some(&self.predicate)
    }

    pub(crate) fn atom_mut(&mut self) -> &mut Atom {
        &mut self.atom
    }
}

impl std::ops::Deref for QueryAtom {
    type Target = Atom;

    fn deref(&self) -> &Self::Target {
        &self.atom
    }
}

/// A bond in a [`QueryGraph`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryBond {
    bond: Bond,
    predicate: QueryNode<BondQueryPredicate>,
}

impl QueryBond {
    pub(crate) fn new(bond: Bond, predicate: QueryNode<BondQueryPredicate>) -> Self {
        Self { bond, predicate }
    }

    /// Concrete bond attributes used by predicates and match callbacks.
    #[must_use]
    pub fn bond(&self) -> &Bond {
        &self.bond
    }

    /// The complete recursive bond predicate tree.
    #[must_use]
    pub fn predicate(&self) -> &QueryNode<BondQueryPredicate> {
        &self.predicate
    }

    /// Bond endpoints in query-graph atom-index space.
    #[must_use]
    pub fn endpoints(&self) -> (usize, usize) {
        (self.bond.begin().index(), self.bond.end().index())
    }

    #[must_use]
    pub fn id(&self) -> crate::BondId {
        self.bond.id()
    }

    #[must_use]
    pub fn begin(&self) -> crate::AtomId {
        self.bond.begin()
    }

    #[must_use]
    pub fn end(&self) -> crate::AtomId {
        self.bond.end()
    }

    #[must_use]
    pub fn query(&self) -> Option<&QueryNode<BondQueryPredicate>> {
        Some(&self.predicate)
    }

    pub(crate) fn bond_mut(&mut self) -> &mut Bond {
        &mut self.bond
    }
}

impl std::ops::Deref for QueryBond {
    type Target = Bond;

    fn deref(&self) -> &Self::Target {
        &self.bond
    }
}

/// A first-class SMARTS/MCS query graph.
#[derive(Debug, Clone, PartialEq)]
pub struct QueryGraph {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<(usize, usize)>>,
    props: BTreeMap<String, String>,
    conformers_2d: Vec<Conformer2D>,
    conformers_3d: Vec<Conformer3D>,
    stereo_groups: Vec<StereoGroup>,
}

impl QueryGraph {
    /// Borrow the search-domain behaviour facade for this query value.
    #[must_use]
    pub const fn operator(&self) -> QueryGraphOperator<'_> {
        QueryGraphOperator::new(self)
    }

    /// Parse SMARTS into a first-class query graph.
    pub fn from_smarts(smarts: &str, params: &crate::SmartsParseParams) -> Result<Self, crate::SmartsParseError> {
        super::smarts_parse::mol_from_smarts(smarts, params)
    }

    /// Build a query graph from a query-bearing construction state. Existing
    /// predicates are preserved; ordinary atoms and bonds receive the same
    /// exact-attribute defaults used by the matcher boundary.
    pub(crate) fn from_query_molecule(molecule: Molecule) -> Result<Self, QueryGraphError> {
        let atoms = molecule
            .atoms()
            .iter()
            .map(|atom| {
                let predicate = atom
                    .query()
                    .cloned()
                    .unwrap_or_else(|| QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())));
                QueryAtom::new(atom.clone(), predicate)
            })
            .collect::<Vec<_>>();
        let bonds = molecule
            .bonds()
            .iter()
            .map(|bond| {
                let predicate = bond.query().cloned().unwrap_or_else(|| {
                    if bond.order() == crate::BondOrder::Unspecified {
                        QueryNode::predicate(BondQueryPredicate::Any)
                    } else {
                        QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
                    }
                });
                Ok(QueryBond::new(bond.clone(), predicate))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for (bond_index, bond) in bonds.iter().enumerate() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if begin >= adjacency.len() || end >= adjacency.len() {
                return Err(QueryGraphError::InvalidBondEndpoint(bond_index));
            }
            adjacency[begin].push((end, bond_index));
            adjacency[end].push((begin, bond_index));
        }
        let mut props = molecule.properties().props().clone();
        if let Some(name) = molecule.properties().name() {
            props.insert("_Name".to_owned(), name.to_owned());
        }
        Ok(Self {
            atoms,
            bonds,
            adjacency,
            props,
            conformers_2d: molecule
                .coordinates_2d()
                .map_or_else(Vec::new, |coords| vec![Conformer2D::new(0, coords.to_vec())]),
            conformers_3d: molecule.conformers_3d().to_vec(),
            stereo_groups: molecule.stereo_groups().to_vec(),
        })
    }

    /// Build a concrete-attribute query graph for APIs that use a molecule as
    /// an exact structural pattern. This is a one-way construction into the
    /// query model; no query graph can be converted back into a `Molecule`.
    pub(crate) fn from_concrete_molecule(molecule: &Molecule) -> Result<Self, QueryGraphError> {
        let atoms = molecule
            .atoms()
            .iter()
            .map(|atom| {
                let mut concrete = atom.clone();
                crate::search::query::replace_atom_with_query_atom(&mut concrete);
                let predicate = concrete
                    .query()
                    .cloned()
                    .expect("replace_atom_with_query_atom installs a predicate");
                concrete.set_query(None);
                QueryAtom::new(concrete, predicate)
            })
            .collect::<Vec<_>>();
        let bonds = molecule
            .bonds()
            .iter()
            .map(|bond| {
                QueryBond::new(
                    bond.clone(),
                    if bond.order() == crate::BondOrder::Unspecified {
                        QueryNode::predicate(BondQueryPredicate::Any)
                    } else {
                        QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
                    },
                )
            })
            .collect::<Vec<_>>();
        let mut props = molecule.properties().props().clone();
        if let Some(name) = molecule.properties().name() {
            props.insert("_Name".to_owned(), name.to_owned());
        }
        Self::from_parts(
            atoms,
            bonds,
            props,
            molecule
                .coordinates_2d()
                .map_or_else(Vec::new, |coords| vec![Conformer2D::new(0, coords.to_vec())]),
            molecule.conformers_3d().to_vec(),
            molecule.stereo_groups().to_vec(),
        )
    }

    pub(crate) fn from_parts(
        atoms: Vec<QueryAtom>,
        bonds: Vec<QueryBond>,
        props: BTreeMap<String, String>,
        conformers_2d: Vec<Conformer2D>,
        conformers_3d: Vec<Conformer3D>,
        stereo_groups: Vec<StereoGroup>,
    ) -> Result<Self, QueryGraphError> {
        let mut adjacency = vec![Vec::new(); atoms.len()];
        for (bond_index, bond) in bonds.iter().enumerate() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            if begin >= adjacency.len() || end >= adjacency.len() {
                return Err(QueryGraphError::InvalidBondEndpoint(bond_index));
            }
            adjacency[begin].push((end, bond_index));
            adjacency[end].push((begin, bond_index));
        }
        Ok(Self {
            atoms,
            bonds,
            adjacency,
            props,
            conformers_2d,
            conformers_3d,
            stereo_groups,
        })
    }

    /// Materialize the graph into the legacy internal molecule workspace.
    ///
    /// This is deliberately crate-private and preserves the query-bearing
    /// atom/bond records. It exists only for source-backed helpers that have
    /// not yet been lowered to operate directly on `QueryGraph`; it is not a
    /// public conversion or a live-molecule mutation path.
    pub(crate) fn to_molecule(&self) -> Result<Molecule, crate::InvariantError> {
        let topology = crate::model::molecule::TopologyBlock {
            atoms: self.atoms.iter().map(|atom| atom.atom.clone()).collect(),
            bonds: self.bonds.iter().map(|bond| bond.bond.clone()).collect(),
            adjacency: crate::AdjacencyList::from_topology(
                self.atoms.len(),
                &self.bonds.iter().map(|bond| bond.bond.clone()).collect::<Vec<_>>(),
            ),
            substance_groups: Vec::new(),
            stereo_groups: self.stereo_groups.clone(),
        };
        let coordinates = crate::model::molecule::CoordinateBlock {
            conformers_2d: self.conformers_2d.clone(),
            conformers_3d: self.conformers_3d.clone(),
            source_coordinate_dim: if !self.conformers_3d.is_empty() {
                Some(crate::CoordinateDimension::ThreeD)
            } else if !self.conformers_2d.is_empty() {
                Some(crate::CoordinateDimension::TwoD)
            } else {
                None
            },
        };
        let mut properties = crate::MoleculeProperties::default();
        for (key, value) in &self.props {
            if key == "_Name" {
                properties = properties.with_name(value.clone());
            } else {
                properties = properties.with_prop(key.clone(), value.clone());
            }
        }
        Molecule::from_blocks(topology, coordinates, properties)
    }

    /// Number of query atoms.
    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    /// Number of query bonds.
    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// Query atoms in stable graph order.
    #[must_use]
    pub fn atoms(&self) -> &[QueryAtom] {
        &self.atoms
    }

    /// Return one query atom by graph index.
    #[must_use]
    pub fn atom(&self, index: usize) -> Option<&QueryAtom> {
        self.atoms.get(index)
    }

    pub(crate) fn atom_mut(&mut self, index: usize) -> Option<&mut Atom> {
        self.atoms.get_mut(index).map(|query_atom| &mut query_atom.atom)
    }

    /// Query bonds in stable graph order.
    #[must_use]
    pub fn bonds(&self) -> &[QueryBond] {
        &self.bonds
    }

    /// Return one query bond by graph index.
    #[must_use]
    pub fn bond(&self, index: usize) -> Option<&QueryBond> {
        self.bonds.get(index)
    }

    #[must_use]
    pub(crate) fn adjacency(&self) -> &[Vec<(usize, usize)>] {
        &self.adjacency
    }

    #[must_use]
    pub(crate) fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub(crate) fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    /// Query name parsed from the SMARTS name field, if present.
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.prop("_Name")
    }

    pub(crate) fn with_name(mut self, name: impl Into<String>) -> Self {
        self.props.insert("_Name".to_owned(), name.into());
        self
    }

    #[must_use]
    pub(crate) fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub(crate) fn conformers_3d(&self) -> &[Conformer3D] {
        &self.conformers_3d
    }

    #[must_use]
    pub(crate) fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.conformers_2d.first().map(Conformer2D::coordinates)
    }

    pub(crate) fn with_2d_coordinate_block(mut self, coords: Vec<[f64; 2]>) -> Result<Self, String> {
        if coords.len() != self.num_atoms() {
            return Err("query graph coordinate count does not match atom count".to_owned());
        }
        self.conformers_2d = vec![Conformer2D::new(0, coords)];
        Ok(self)
    }

    #[must_use]
    pub(crate) fn stereo_groups(&self) -> &[StereoGroup] {
        &self.stereo_groups
    }

    pub(crate) fn has_directional_bonds(&self) -> bool {
        self.bonds
            .iter()
            .any(|bond| bond.bond().direction() != crate::BondDirection::None)
    }

    /// Remove parser-only properties after the source parser has completed.
    /// This is the query-graph equivalent of RDKit CleanupAfterParsing for
    /// state emitted by the SMARTS grammar itself.
    pub(crate) fn cleanup_parser_state(&mut self) {
        for atom in &mut self.atoms {
            atom.atom_mut().clear_prop("_RingClosures");
            atom.atom_mut().clear_prop(crate::notation::smiles::SMILES_START_PROP);
        }
        for bond in &mut self.bonds {
            bond.bond_mut()
                .clear_prop(crate::notation::smiles::UNSPECIFIED_ORDER_PROP);
            bond.bond_mut()
                .clear_prop(crate::notation::smiles::CXSMILES_BOND_IDX_PROP);
        }
    }

    /// Render this query using the canonical SMARTS writer.
    pub fn to_smarts(&self, params: &crate::SmilesWriteParams) -> Result<String, crate::SmartsWriteError> {
        self.operator().to_smarts(params)
    }

    pub fn pattern_fingerprint(
        &self,
        params: &crate::PatternFingerprintParams,
    ) -> Result<crate::Fingerprint, crate::FingerprintError> {
        self.operator().pattern_fingerprint(params)
    }

    /// Render a CXSMARTS query.
    pub fn to_cx_smarts(&self, params: &crate::SmilesWriteParams) -> Result<String, crate::SmartsWriteError> {
        self.operator().to_cx_smarts(params)
    }

    /// Render a selected query fragment as SMARTS.
    pub fn fragment_to_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        self.operator().fragment_to_smarts(params, atom_ids, bond_ids)
    }

    /// Render a selected query fragment as CXSMARTS.
    pub fn fragment_to_cx_smarts(
        &self,
        params: &crate::SmilesWriteParams,
        atom_ids: &[crate::AtomId],
        bond_ids: Option<&[crate::BondId]>,
    ) -> Result<String, crate::SmartsWriteError> {
        self.operator().fragment_to_cx_smarts(params, atom_ids, bond_ids)
    }

    /// Compile this graph for repeated matching.
    #[must_use]
    pub fn compile(self) -> CompiledQuery {
        CompiledQuery::compile(self)
    }

    /// Match this query graph against a concrete molecule.
    pub fn matches(&self, molecule: &Molecule) -> Vec<MatchResult> {
        self.operator().matches(molecule)
    }
}

/// Query graph construction failed because a required predicate was absent.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryGraphError {
    #[error("query atom {0} has no atom predicate")]
    MissingAtomPredicate(usize),
    #[error("query bond {0} has no bond predicate")]
    MissingBondPredicate(usize),
    #[error("query bond {0} has an invalid atom endpoint")]
    InvalidBondEndpoint(usize),
}

/// A compiled query ready for repeated matching.
#[derive(Debug, Clone, PartialEq)]
pub struct CompiledQuery {
    graph: QueryGraph,
}

impl CompiledQuery {
    /// Compile a query graph without changing its predicate semantics.
    #[must_use]
    pub fn compile(graph: QueryGraph) -> Self {
        Self { graph }
    }

    /// Access the source query graph.
    #[must_use]
    pub fn graph(&self) -> &QueryGraph {
        &self.graph
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.graph.num_atoms()
    }

    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.graph.num_bonds()
    }

    /// Match the compiled query against a concrete molecule.
    pub fn matches(&self, molecule: &Molecule) -> Vec<MatchResult> {
        crate::search::substruct::get_substruct_matches(molecule, &self.graph)
    }
}

impl MatchResult {
    /// Create a match result from query-to-target atom and bond mappings.
    #[must_use]
    pub fn new(atom_mapping: Vec<usize>, bond_mapping: Vec<usize>) -> Self {
        Self {
            atom_mapping,
            bond_mapping,
        }
    }
}

impl McsResult {
    /// Construct MCS result metadata for a completed or canceled search.
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

/// A match mapping from query graph indices to target molecule indices.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatchResult {
    pub atom_mapping: Vec<usize>,
    pub bond_mapping: Vec<usize>,
}

/// Result metadata for a maximum-common-substructure search.
#[derive(Debug, Clone, PartialEq)]
pub struct McsResult {
    pub query: QueryGraph,
    pub atom_count: usize,
    pub bond_count: usize,
    pub completed: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smarts_parser_produces_first_class_query_graph() {
        let graph = crate::mol_from_smarts("[#6]-[#8]", &crate::SmartsParseParams::default()).expect("SMARTS query");
        assert_eq!(graph.num_atoms(), 2);
        assert_eq!(graph.num_bonds(), 1);
        assert!(matches!(
            graph.atoms()[0].predicate(),
            QueryNode::Predicate(AtomQueryPredicate::AtomicNumber(6))
        ));
    }

    #[test]
    fn compiled_query_matches_concrete_molecule() {
        let graph = crate::mol_from_smarts("[#6]-[#8]", &crate::SmartsParseParams::default()).expect("SMARTS query");
        let target = Molecule::from_smiles("CCO").expect("target molecule");
        let atom_count = graph.num_atoms();
        let matches = graph.compile().matches(&target);
        assert!(!matches.is_empty());
        assert_eq!(matches[0].atom_mapping.len(), atom_count);
    }
}
