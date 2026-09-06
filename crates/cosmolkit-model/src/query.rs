//! Query values used by SMARTS, MCS, and substructure algorithms.
//!
//! This module contains only query data and local graph validation. Parsing,
//! matching, serialization, and compilation belong to `cosmolkit-core` (or
//! another domain implementation crate); query data is never lowered back to a
//! concrete `Molecule`.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Hybridization};

use crate::{
    Atom, AtomId, Bond, BondId, Conformer2D, Conformer3D, CoordinateBlock,
    CoordinateValidationError, StereoGroup,
};

/// A recursive Boolean query tree over a predicate type.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryNode<T> {
    Predicate(T),
    And(Vec<QueryNode<T>>),
    Or(Vec<QueryNode<T>>),
    Xor(Vec<QueryNode<T>>),
    Not(Box<QueryNode<T>>),
}

impl<T> QueryNode<T> {
    #[must_use]
    pub fn predicate(predicate: T) -> Self {
        Self::Predicate(predicate)
    }

    #[must_use]
    pub fn and(children: Vec<Self>) -> Self {
        Self::And(children)
    }

    #[must_use]
    pub fn or(children: Vec<Self>) -> Self {
        Self::Or(children)
    }

    #[must_use]
    pub fn xor(children: Vec<Self>) -> Self {
        Self::Xor(children)
    }

    #[must_use]
    pub fn not(child: Self) -> Self {
        Self::Not(Box::new(child))
    }

    /// Append a child to a composite node.
    #[doc(hidden)]
    pub fn add_child(&mut self, child: Self) {
        match self {
            Self::And(children) | Self::Or(children) | Self::Xor(children) => children.push(child),
            Self::Predicate(_) | Self::Not(_) => {
                unreachable!("only child-vector query nodes accept children")
            }
        }
    }

    /// Toggle the canonical outer negation used by the source query merger.
    #[doc(hidden)]
    pub fn set_negation(&mut self, negated: bool) {
        match (negated, matches!(self, Self::Not(_))) {
            (true, false) => {
                let child = std::mem::replace(self, Self::And(Vec::new()));
                *self = Self::Not(Box::new(child));
            }
            (false, true) => {
                let Self::Not(child) = std::mem::replace(self, Self::And(Vec::new())) else {
                    unreachable!()
                };
                *self = *child;
            }
            _ => {}
        }
    }

    #[doc(hidden)]
    pub fn is_negated(&self) -> bool {
        matches!(self, Self::Not(_))
    }
}

/// Supported RDKit atom data functions for typed range-query leaves.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomRangeDataFunction {
    ExplicitDegree,
    NonHydrogenDegree,
    TotalDegree,
    TotalValence,
    NumAtomRings,
    NumHeteroatomNeighbors,
    NumAliphaticHeteroatomNeighbors,
    MinRingSize,
    RingBondCount,
    ImplicitHydrogenCount,
    FormalCharge,
    NegativeFormalCharge,
    AtomRingSize {
        lower: i32,
        upper: i32,
        lower_open: bool,
        upper_open: bool,
    },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AtomRangeBounds {
    LessEqual(i32),
    GreaterEqual(i32),
    Inclusive {
        lower: i32,
        upper: i32,
        lower_open: bool,
        upper_open: bool,
    },
}

/// Canonical typed representation of RDKit's `ATOM_RANGE_QUERY`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomRangeQuery {
    bounds: AtomRangeBounds,
    data_function: AtomRangeDataFunction,
}

impl AtomRangeQuery {
    #[must_use]
    pub const fn new(bounds: AtomRangeBounds, data_function: AtomRangeDataFunction) -> Self {
        Self {
            bounds,
            data_function,
        }
    }

    #[must_use]
    pub const fn bounds(&self) -> AtomRangeBounds {
        self.bounds
    }

    #[must_use]
    pub const fn data_function(&self) -> AtomRangeDataFunction {
        self.data_function
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn writer_parts(&self) -> (AtomRangeBounds, AtomRangeDataFunction) {
        (self.bounds, self.data_function)
    }
}

/// Atom-level SMARTS / MolFile query predicates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AtomQueryPredicate {
    Any,
    AtomicNumber(u8),
    AtomType { atomic_number: u8, aromatic: bool },
    AtomicNumberIn(Vec<u8>),
    AtomicNumberNotIn(Vec<u8>),
    FormalCharge(i8),
    NegativeFormalCharge(i8),
    NumRadicalElectrons(u8),
    HasChiralTag,
    MissingChiralTag,
    Isotope(u16),
    HydrogenCount(u8),
    HasImplicitHydrogen,
    ImplicitHydrogenCount(u8),
    ImplicitHydrogenCountLessEqual(u8),
    ImplicitValence(i32),
    ExplicitValence(i32),
    ExplicitDegree(u8),
    ExplicitDegreeLessEqual(u8),
    NonHydrogenDegree(u32),
    NonHydrogenDegreeLessEqual(u32),
    NonHydrogenDegreeGreaterEqual(u32),
    HeavyAtomDegree(u32),
    NumHeteroatomNeighbors(u8),
    HasHeteroatomNeighbors,
    NumAliphaticHeteroatomNeighbors(u8),
    HasAliphaticHeteroatomNeighbors,
    RingBondCount(u32),
    RingBondCountLessEqual(u8),
    HasRingBond,
    IsBridgehead,
    IsAromatic(bool),
    IsUnsaturated,
    RecursiveSmarts(RecursiveStructureQuery),
    HasProperty(String),
    PropertyValue { name: String, value: String },
    RGroupLabel(u32),
    MolFileAlias(String),
    HybridizationMatch(Hybridization),
    TotalDegree(u8),
    TotalDegreeLessEqual(u8),
    TotalDegreeGreaterEqual(u8),
    TotalValence(u8),
    TotalValenceLessEqual(u8),
    TotalValenceGreaterEqual(u8),
    InRing,
    NumAtomRings(i32),
    InRingOfSize(u8),
    InRingOfSizeLessEqual(u8),
    InRingOfSizeGreaterEqual(u8),
    SmallestRingSize(u8),
    SmallestRingSizeLessEqual(u8),
    SmallestRingSizeGreaterEqual(u8),
    Mass(u16),
    ChiralTagMatch(ChiralTag),
    ChiralPermutationMatch(u32),
    DegreeLessEqual(u8),
    DegreeGreaterEqual(u8),
    Range(AtomRangeQuery),
    UnsupportedFeature(&'static str),
}

/// Bond-level SMARTS / MolFile query predicates.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BondQueryPredicate {
    Any,
    Order(BondOrder),
    OrderIn(Vec<BondOrder>),
    IsAromatic(bool),
    IsInRing(bool),
    Direction(BondDirection),
    Stereo(BondStereo),
    HasStereo,
    IsConjugated,
    NumRingBonds(i32),
    InRingOfSize(i32),
    MinRingSize(i32),
    NumRingBondsGreaterEqual(u8),
    NumRingBondsLessEqual(u8),
    MolFileQueryCode(u32),
    HasProperty(String),
    PropertyValue { name: String, value: String },
    UnsupportedFeature(&'static str),
}

/// Recursive SMARTS query data.
#[derive(Debug, Clone, PartialEq)]
pub struct RecursiveStructureQuery {
    query_mol: Option<Box<QueryGraph>>,
    source_smarts: Option<String>,
    atom_indices: BTreeSet<i32>,
    serial_number: u32,
}

impl RecursiveStructureQuery {
    #[must_use]
    pub fn new() -> Self {
        Self {
            query_mol: None,
            source_smarts: None,
            atom_indices: BTreeSet::new(),
            serial_number: 0,
        }
    }

    #[must_use]
    pub fn from_query_graph(query: QueryGraph, serial_number: u32) -> Self {
        Self {
            query_mol: Some(Box::new(query)),
            source_smarts: None,
            atom_indices: BTreeSet::new(),
            serial_number,
        }
    }

    /// Compatibility name for source-shaped recursive-query construction.
    #[must_use]
    pub fn from_molecule(query: QueryGraph, serial_number: u32) -> Self {
        Self::from_query_graph(query, serial_number)
    }

    #[doc(hidden)]
    #[must_use]
    pub fn atom_index(atom: &Atom) -> i32 {
        atom.id().index() as i32
    }

    #[must_use]
    pub fn from_smarts(smarts: impl Into<String>, query: QueryGraph, serial_number: u32) -> Self {
        let mut value = Self::from_query_graph(query, serial_number);
        value.source_smarts = Some(smarts.into());
        value
    }

    #[must_use]
    pub fn query_mol(&self) -> Option<&QueryGraph> {
        self.query_mol.as_deref()
    }

    #[doc(hidden)]
    pub fn set_query_mol(&mut self, query: QueryGraph) {
        self.query_mol = Some(Box::new(query));
    }

    #[doc(hidden)]
    pub fn query_mol_mut(&mut self) -> Option<&mut QueryGraph> {
        self.query_mol.as_deref_mut()
    }

    #[must_use]
    pub fn source_smarts(&self) -> Option<&str> {
        self.source_smarts.as_deref()
    }

    #[doc(hidden)]
    pub fn insert_atom_index(&mut self, index: i32) {
        self.atom_indices.insert(index);
    }

    #[doc(hidden)]
    #[must_use]
    pub fn contains_atom_index(&self, index: i32) -> bool {
        self.atom_indices.contains(&index)
    }

    #[doc(hidden)]
    #[must_use]
    pub fn copy_query(&self) -> Self {
        self.clone()
    }

    #[must_use]
    pub const fn serial_number(&self) -> u32 {
        self.serial_number
    }
}

impl Default for RecursiveStructureQuery {
    fn default() -> Self {
        Self::new()
    }
}

impl Eq for RecursiveStructureQuery {}

/// A query atom combines concrete atom attributes with a query predicate tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryAtom {
    atom: Atom,
    predicate: QueryNode<AtomQueryPredicate>,
}

impl QueryAtom {
    #[must_use]
    pub fn new(id: AtomId, spec: crate::AtomSpec) -> Self {
        let predicate = QueryNode::predicate(AtomQueryPredicate::AtomicNumber(
            spec.element().atomic_number(),
        ));
        Self::from_parts(Atom::from_spec(id, spec), predicate)
    }

    #[must_use]
    pub fn from_parts(atom: Atom, predicate: QueryNode<AtomQueryPredicate>) -> Self {
        Self { atom, predicate }
    }

    #[must_use]
    pub fn atom(&self) -> &Atom {
        &self.atom
    }

    #[doc(hidden)]
    pub fn atom_mut(&mut self) -> &mut Atom {
        &mut self.atom
    }

    #[must_use]
    pub fn predicate(&self) -> &QueryNode<AtomQueryPredicate> {
        &self.predicate
    }

    #[doc(hidden)]
    pub fn predicate_mut(&mut self) -> &mut QueryNode<AtomQueryPredicate> {
        &mut self.predicate
    }

    #[doc(hidden)]
    pub fn set_predicate(&mut self, predicate: QueryNode<AtomQueryPredicate>) {
        self.predicate = predicate;
    }

    #[must_use]
    pub fn index(&self) -> usize {
        self.atom.id().index()
    }

    #[must_use]
    pub fn id(&self) -> AtomId {
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
}

/// A query bond combines concrete bond attributes with a query predicate tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryBond {
    bond: Bond,
    predicate: QueryNode<BondQueryPredicate>,
}

impl QueryBond {
    #[must_use]
    pub fn new(id: BondId, spec: crate::BondSpec) -> Self {
        let predicate = if spec.order() == BondOrder::Unspecified {
            QueryNode::predicate(BondQueryPredicate::Any)
        } else {
            QueryNode::predicate(BondQueryPredicate::Order(spec.order()))
        };
        Self::from_parts(Bond::from_spec(id, spec), predicate)
    }

    #[must_use]
    pub fn from_parts(bond: Bond, predicate: QueryNode<BondQueryPredicate>) -> Self {
        Self { bond, predicate }
    }

    #[must_use]
    pub fn bond(&self) -> &Bond {
        &self.bond
    }

    #[doc(hidden)]
    pub fn bond_mut(&mut self) -> &mut Bond {
        &mut self.bond
    }

    #[must_use]
    pub fn predicate(&self) -> &QueryNode<BondQueryPredicate> {
        &self.predicate
    }

    #[doc(hidden)]
    pub fn predicate_mut(&mut self) -> &mut QueryNode<BondQueryPredicate> {
        &mut self.predicate
    }

    #[doc(hidden)]
    pub fn set_predicate(&mut self, predicate: QueryNode<BondQueryPredicate>) {
        self.predicate = predicate;
    }

    #[must_use]
    pub fn endpoints(&self) -> (usize, usize) {
        (self.bond.begin().index(), self.bond.end().index())
    }

    #[must_use]
    pub fn id(&self) -> BondId {
        self.bond.id()
    }

    #[must_use]
    pub fn begin(&self) -> AtomId {
        self.bond.begin()
    }

    #[must_use]
    pub fn end(&self) -> AtomId {
        self.bond.end()
    }
}

/// First-class SMARTS/MCS query graph value.
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
    #[must_use]
    pub fn from_parts(
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
        let graph = Self {
            atoms,
            bonds,
            adjacency,
            props,
            conformers_2d,
            conformers_3d,
            stereo_groups,
        };
        graph.validate()?;
        Ok(graph)
    }

    /// Validate local query-graph structure without consulting a live
    /// molecule or runtime cache.
    pub fn validate(&self) -> Result<(), QueryGraphError> {
        for (position, atom) in self.atoms.iter().enumerate() {
            if atom.id() != AtomId::new(position) {
                return Err(QueryGraphError::AtomIdMismatch {
                    position,
                    id: atom.id(),
                });
            }
        }
        for (position, bond) in self.bonds.iter().enumerate() {
            if bond.id() != BondId::new(position) {
                return Err(QueryGraphError::BondIdMismatch {
                    position,
                    id: bond.id(),
                });
            }
            for atom in [bond.begin(), bond.end()] {
                if atom.index() >= self.atoms.len() {
                    return Err(QueryGraphError::InvalidBondEndpoint(position));
                }
            }
            if let Some([begin, end]) = bond.bond().stereo_atoms()
                && (begin.index() >= self.atoms.len() || end.index() >= self.atoms.len())
            {
                return Err(QueryGraphError::StereoAtomOutOfRange {
                    bond: bond.id(),
                    begin,
                    end,
                    atom_count: self.atoms.len(),
                });
            }
        }

        let coordinates = CoordinateBlock {
            conformers_2d: self.conformers_2d.clone(),
            conformers_3d: self.conformers_3d.clone(),
            source_coordinate_dim: None,
        };
        coordinates
            .validate_for_atom_count(self.atoms.len())
            .map_err(QueryGraphError::CoordinateValidation)?;

        for group in &self.stereo_groups {
            for atom in group.atoms() {
                if atom.index() >= self.atoms.len() {
                    return Err(QueryGraphError::StereoGroupAtomOutOfRange {
                        atom: *atom,
                        atom_count: self.atoms.len(),
                    });
                }
            }
            for bond in group.bonds() {
                if bond.index() >= self.bonds.len() {
                    return Err(QueryGraphError::StereoGroupBondOutOfRange {
                        bond: *bond,
                        bond_count: self.bonds.len(),
                    });
                }
            }
        }

        let mut expected_adjacency = vec![Vec::new(); self.atoms.len()];
        for (bond_index, bond) in self.bonds.iter().enumerate() {
            let begin = bond.begin().index();
            let end = bond.end().index();
            expected_adjacency[begin].push((end, bond_index));
            expected_adjacency[end].push((begin, bond_index));
        }
        if self.adjacency != expected_adjacency {
            return Err(QueryGraphError::AdjacencyMismatch);
        }
        Ok(())
    }

    #[must_use]
    pub fn num_atoms(&self) -> usize {
        self.atoms.len()
    }

    #[must_use]
    pub fn num_bonds(&self) -> usize {
        self.bonds.len()
    }

    #[must_use]
    pub fn atoms(&self) -> &[QueryAtom] {
        &self.atoms
    }

    #[must_use]
    pub fn atom(&self, index: usize) -> Option<&QueryAtom> {
        self.atoms.get(index)
    }

    #[doc(hidden)]
    pub fn atoms_mut(&mut self) -> &mut [QueryAtom] {
        &mut self.atoms
    }

    #[doc(hidden)]
    pub fn atom_mut(&mut self, index: usize) -> Option<&mut QueryAtom> {
        self.atoms.get_mut(index)
    }

    #[must_use]
    pub fn bonds(&self) -> &[QueryBond] {
        &self.bonds
    }

    #[must_use]
    pub fn bond(&self, index: usize) -> Option<&QueryBond> {
        self.bonds.get(index)
    }

    #[doc(hidden)]
    pub fn bonds_mut(&mut self) -> &mut [QueryBond] {
        &mut self.bonds
    }

    #[must_use]
    pub fn adjacency(&self) -> &[Vec<(usize, usize)>] {
        &self.adjacency
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    /// Remove a parser/runtime property from this detached query value.
    #[doc(hidden)]
    pub fn clear_prop(&mut self, key: &str) {
        self.props.remove(key);
    }

    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.prop("_Name")
    }

    #[doc(hidden)]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.props.insert("_Name".to_owned(), name.into());
        self
    }

    #[doc(hidden)]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[doc(hidden)]
    pub fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        self.props.insert(key.into(), value.into());
    }

    #[must_use]
    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.conformers_3d
    }

    #[doc(hidden)]
    pub fn add_conformer_3d(&mut self, conformer: Conformer3D) -> Result<(), String> {
        if conformer.coordinates().len() != self.num_atoms() {
            return Err("query graph coordinate count does not match atom count".to_owned());
        }
        self.conformers_3d.push(conformer);
        Ok(())
    }

    #[must_use]
    pub fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.conformers_2d.first().map(Conformer2D::coordinates)
    }

    #[doc(hidden)]
    pub fn with_2d_coordinate_block(mut self, coords: Vec<[f64; 2]>) -> Result<Self, String> {
        if coords.len() != self.num_atoms() {
            return Err("query graph coordinate count does not match atom count".to_owned());
        }
        self.conformers_2d = vec![Conformer2D::new(0, coords)];
        Ok(self)
    }

    #[must_use]
    pub fn stereo_groups(&self) -> &[StereoGroup] {
        &self.stereo_groups
    }

    #[doc(hidden)]
    pub fn add_stereo_group(&mut self, group: StereoGroup) {
        self.stereo_groups.push(group);
    }
}

/// Query graph construction failed because a graph-local constraint was invalid.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryGraphError {
    #[error("query atom at position {position} has id {id}, expected {position}")]
    AtomIdMismatch { position: usize, id: AtomId },
    #[error("query bond at position {position} has id {id}, expected {position}")]
    BondIdMismatch { position: usize, id: BondId },
    #[error("query atom {0} has no atom predicate")]
    MissingAtomPredicate(usize),
    #[error("query bond {0} has no bond predicate")]
    MissingBondPredicate(usize),
    #[error("query bond {0} has an invalid atom endpoint")]
    InvalidBondEndpoint(usize),
    #[error(
        "query bond {bond} has stereo atom references {begin}-{end} outside {atom_count} atoms"
    )]
    StereoAtomOutOfRange {
        bond: BondId,
        begin: AtomId,
        end: AtomId,
        atom_count: usize,
    },
    #[error("query stereo group references atom {atom} outside {atom_count} atoms")]
    StereoGroupAtomOutOfRange { atom: AtomId, atom_count: usize },
    #[error("query stereo group references bond {bond} outside {bond_count} bonds")]
    StereoGroupBondOutOfRange { bond: BondId, bond_count: usize },
    #[error("query graph coordinate validation failed: {0}")]
    CoordinateValidation(CoordinateValidationError),
    #[error("query graph adjacency does not match its atom and bond rows")]
    AdjacencyMismatch,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{AtomSpec, BondSpec, Element, StereoGroupKind};

    fn carbon(id: usize) -> QueryAtom {
        QueryAtom::new(AtomId::new(id), AtomSpec::new(Element::C))
    }

    #[test]
    fn from_parts_rejects_noncanonical_query_atom_ids() {
        let result = QueryGraph::from_parts(
            vec![carbon(1)],
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        );

        assert!(matches!(
            result,
            Err(QueryGraphError::AtomIdMismatch {
                position: 0,
                id
            }) if id == AtomId::new(1)
        ));
    }

    #[test]
    fn from_parts_rejects_coordinate_and_stereo_group_misalignment() {
        let coordinate_result = QueryGraph::from_parts(
            vec![carbon(0)],
            Vec::new(),
            BTreeMap::new(),
            vec![Conformer2D::new(0, Vec::new())],
            Vec::new(),
            Vec::new(),
        );
        assert!(matches!(
            coordinate_result,
            Err(QueryGraphError::CoordinateValidation(
                CoordinateValidationError::RowCount { .. }
            ))
        ));

        let group_result = QueryGraph::from_parts(
            vec![carbon(0)],
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            vec![StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(1)],
                Vec::new(),
            )],
        );
        assert!(matches!(
            group_result,
            Err(QueryGraphError::StereoGroupAtomOutOfRange {
                atom,
                atom_count: 1
            }) if atom == AtomId::new(1)
        ));
    }

    #[test]
    fn validate_rejects_stale_query_adjacency_after_detached_editing() {
        let graph = QueryGraph {
            atoms: vec![carbon(0), carbon(1)],
            bonds: vec![QueryBond::new(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            )],
            adjacency: vec![Vec::new(), Vec::new()],
            props: BTreeMap::new(),
            conformers_2d: Vec::new(),
            conformers_3d: Vec::new(),
            stereo_groups: Vec::new(),
        };

        assert_eq!(graph.validate(), Err(QueryGraphError::AdjacencyMismatch));
    }
}
