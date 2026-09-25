//! Query values used by SMARTS, MCS, and substructure algorithms.
//!
//! This module contains only query data and local graph validation. SMARTS
//! parsing, writing, matching, serialization, and compilation belong uniquely
//! to `cosmolkit-search`; query data is never lowered back to a concrete
//! `Molecule`.

use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element, Hybridization};

use crate::atom::AtomProperties;
use crate::{
    Atom, AtomId, AtomPropertyError, Bond, BondId, Conformer2D, Conformer3D, CoordinateBlock,
    CoordinateValidationError, MappingValidationError, StereoGroup, TemplateAttachmentOrder,
    TemplateAttachmentOrderError, TopologyBlock, TopologyMapping,
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
        // BEGIN RDKIT CPP FUNCTION Queries::Query::addChild
        // RDKit✔️✔️: //! adds a child to our list of children
        // RDKit✔️✔️: void addChild(CHILD_TYPE child) { this->d_children.push_back(child); }
        // END RDKIT CPP FUNCTION Queries::Query::addChild
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
        // BEGIN RDKIT CPP FUNCTION Queries::Query::setNegation
        // RDKit✔️✔️: //! sets whether or not we are negated
        // RDKit✔️✔️: void setNegation(bool what) { this->df_negate = what; }
        // END RDKIT CPP FUNCTION Queries::Query::setNegation
        // `Not` stores the same Boolean outer state in the Rust sum type.
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
        // BEGIN RDKIT CPP FUNCTION RangeQuery::RangeQuery
        // RDKit✔️✔️: //! construct and set the lower and upper bounds
        // RDKit✔️✔️: RangeQuery(MatchFuncArgType lower, MatchFuncArgType upper)
        // RDKit✔️✔️:     : d_upper(upper), d_lower(lower), df_upperOpen(true), df_lowerOpen(true) {
        // RDKit✔️✔️:   this->df_negate = false;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RangeQuery::RangeQuery
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
    FormalCharge(i32),
    NegativeFormalCharge(i32),
    NumRadicalElectrons(u8),
    HasChiralTag,
    MissingChiralTag,
    Isotope(i32),
    HydrogenCount(i32),
    HasImplicitHydrogen,
    ImplicitHydrogenCount(i32),
    ImplicitHydrogenCountLessEqual(u8),
    ImplicitValence(i32),
    ExplicitValence(i32),
    ExplicitDegree(i32),
    ExplicitDegreeLessEqual(u8),
    NonHydrogenDegree(u32),
    NonHydrogenDegreeLessEqual(u32),
    NonHydrogenDegreeGreaterEqual(u32),
    HeavyAtomDegree(u32),
    NumHeteroatomNeighbors(i32),
    HasHeteroatomNeighbors,
    NumAliphaticHeteroatomNeighbors(i32),
    HasAliphaticHeteroatomNeighbors,
    RingBondCount(i32),
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
    TotalDegree(i32),
    TotalDegreeLessEqual(u8),
    TotalDegreeGreaterEqual(u8),
    TotalValence(i32),
    TotalValenceLessEqual(u8),
    TotalValenceGreaterEqual(u8),
    InRing,
    NumAtomRings(i32),
    InRingOfSize(i32),
    InRingOfSizeLessEqual(u8),
    InRingOfSizeGreaterEqual(u8),
    SmallestRingSize(i32),
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
#[derive(Debug, PartialEq)]
pub struct RecursiveStructureQuery {
    query_graph: Option<Box<QueryGraph>>,
    source_smarts: Option<String>,
    atom_indices: BTreeSet<i32>,
    serial_number: u32,
}

impl RecursiveStructureQuery {
    #[must_use]
    pub fn new() -> Self {
        // BEGIN RDKIT CPP FUNCTION RecursiveStructureQuery::RecursiveStructureQuery
        // RDKit✔️✔️: RecursiveStructureQuery() : Queries::SetQuery<int, Atom const *, true>() {
        // RDKit✔️✔️:   setDataFunc(getAtIdx);
        // RDKit✔️✔️:   setDescription("RecursiveStructure");
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::RecursiveStructureQuery
        Self {
            query_graph: None,
            source_smarts: None,
            atom_indices: BTreeSet::new(),
            serial_number: 0,
        }
    }

    #[must_use]
    pub fn from_query_graph(query_graph: QueryGraph, serial_number: u32) -> Self {
        // BEGIN RDKIT CPP FUNCTION RecursiveStructureQuery::RecursiveStructureQuery
        // RDKit✔️✔️: RecursiveStructureQuery(ROMol const *query, unsigned int serialNumber = 0)
        // RDKit✔️✔️:     : Queries::SetQuery<int, Atom const *, true>(),
        // RDKit✔️✔️:       d_serialNumber(serialNumber) {
        // RDKit✔️✔️:   setQueryMol(query);
        // RDKit✔️✔️:   setDataFunc(getAtIdx);
        // RDKit✔️✔️:   setDescription("RecursiveStructure");
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::RecursiveStructureQuery
        Self {
            query_graph: Some(Box::new(query_graph)),
            source_smarts: None,
            atom_indices: BTreeSet::new(),
            serial_number,
        }
    }

    #[must_use]
    pub fn with_source_smarts(mut self, smarts: impl Into<String>) -> Self {
        self.source_smarts = Some(smarts.into());
        self
    }

    #[must_use]
    pub fn query_graph(&self) -> Option<&QueryGraph> {
        // BEGIN RDKIT CPP FUNCTION RecursiveStructureQuery::getQueryMol
        // RDKit✔️✔️: //! returns a pointer to our query molecule
        // RDKit✔️✔️: ROMol const *getQueryMol() const { return dp_queryMol.get(); }
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::getQueryMol
        self.query_graph.as_deref()
    }

    #[doc(hidden)]
    pub fn set_query_graph(&mut self, query_graph: QueryGraph) {
        // BEGIN RDKIT CPP FUNCTION RecursiveStructureQuery::setQueryMol
        // RDKit✔️✔️: void setQueryMol(ROMol const *query) { dp_queryMol.reset(query); }
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::setQueryMol
        self.query_graph = Some(Box::new(query_graph));
    }

    #[doc(hidden)]
    pub fn query_graph_mut(&mut self) -> Option<&mut QueryGraph> {
        self.query_graph.as_deref_mut()
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

    #[must_use]
    pub const fn serial_number(&self) -> u32 {
        self.serial_number
    }
}

impl Clone for RecursiveStructureQuery {
    fn clone(&self) -> Self {
        // BEGIN RDKIT CPP FUNCTION RecursiveStructureQuery::copy
        // RDKit✔️✔️: RecursiveStructureQuery *res = new RecursiveStructureQuery();
        // RDKit✔️✔️: res->dp_queryMol.reset(new ROMol(*dp_queryMol, true));
        // RDKit✔️✔️: for (i = d_set.begin(); i != d_set.end(); i++) {
        // RDKit✔️✔️:   res->insert(*i);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: res->d_serialNumber = d_serialNumber;
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::copy
        Self {
            query_graph: self.query_graph.clone(),
            source_smarts: self.source_smarts.clone(),
            atom_indices: self.atom_indices.clone(),
            serial_number: self.serial_number,
        }
    }
}

impl Default for RecursiveStructureQuery {
    fn default() -> Self {
        Self::new()
    }
}

impl Eq for RecursiveStructureQuery {}

/// Distinguishes source queries from ordinary carriers in uniform query storage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QueryPredicateOrigin {
    Explicit,
    CarrierDerived,
}

/// Chemical identity carried by one query atom.
///
/// Ordinary [`Atom`] values remain `Element`-constrained. Query carriers may
/// retain a source numeric atomic number that has no corresponding Element.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum QueryAtomIdentity {
    Element(Element),
    AtomicNumber(u8),
}

impl QueryAtomIdentity {
    /// Store representable atomic numbers canonically as Elements.
    #[must_use]
    pub const fn from_atomic_number(atomic_number: u8) -> Self {
        match Element::from_atomic_number(atomic_number) {
            Some(element) => Self::Element(element),
            None => Self::AtomicNumber(atomic_number),
        }
    }

    #[must_use]
    pub const fn atomic_number(self) -> u8 {
        match self {
            Self::Element(element) => element.atomic_number(),
            Self::AtomicNumber(atomic_number) => atomic_number,
        }
    }

    #[must_use]
    pub const fn element(self) -> Option<Element> {
        match self {
            Self::Element(element) => Some(element),
            Self::AtomicNumber(_) => None,
        }
    }

    const fn canonicalized(self) -> Self {
        match self {
            Self::Element(element) => Self::Element(element),
            Self::AtomicNumber(atomic_number) => Self::from_atomic_number(atomic_number),
        }
    }
}

/// A query carrier cannot be represented as an ordinary Element-only Atom.
#[derive(Debug, Clone, Copy, PartialEq, Eq, thiserror::Error)]
pub enum QueryAtomConversionError {
    #[error("query atom {atom} has non-Element atomic number {atomic_number}")]
    NonElementAtomicNumber { atom: AtomId, atomic_number: u8 },
}

/// A query atom combines carrier attributes, a predicate tree, and its origin.
///
/// Equality compares stored representation, including predicate origin. An
/// explicit query and an ordinary carrier with identical attributes and trees
/// are unequal because their source `hasQuery()` behavior differs. Equality
/// does not establish chemical or query-matching equivalence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryAtom {
    id: AtomId,
    identity: QueryAtomIdentity,
    properties: AtomProperties,
    predicate: QueryNode<AtomQueryPredicate>,
    predicate_origin: QueryPredicateOrigin,
}

impl QueryAtom {
    #[must_use]
    pub fn new(id: AtomId, spec: crate::AtomSpec) -> Self {
        let (element, properties) = spec.into_query_parts();
        let predicate =
            QueryNode::predicate(AtomQueryPredicate::AtomicNumber(element.atomic_number()));
        Self::from_properties(
            id,
            QueryAtomIdentity::Element(element),
            properties,
            predicate,
            QueryPredicateOrigin::Explicit,
        )
    }

    #[must_use]
    pub fn from_parts(atom: Atom, predicate: QueryNode<AtomQueryPredicate>) -> Self {
        let (id, element, properties) = atom.into_query_parts();
        Self::from_properties(
            id,
            QueryAtomIdentity::Element(element),
            properties,
            predicate,
            QueryPredicateOrigin::Explicit,
        )
    }

    /// Construct the uniform query carrier used internally when Molfile input
    /// contains a mixture of ordinary and query atoms.
    #[doc(hidden)]
    #[must_use]
    pub fn from_carrier_parts(atom: Atom, predicate: QueryNode<AtomQueryPredicate>) -> Self {
        let (id, element, properties) = atom.into_query_parts();
        Self::from_properties(
            id,
            QueryAtomIdentity::Element(element),
            properties,
            predicate,
            QueryPredicateOrigin::CarrierDerived,
        )
    }

    /// Construct a query carrier with independent typed identity and predicate.
    #[must_use]
    pub fn from_identity_parts(
        id: AtomId,
        identity: QueryAtomIdentity,
        predicate: QueryNode<AtomQueryPredicate>,
    ) -> Self {
        // BEGIN RDKIT CPP FUNCTION Atom::Atom(unsigned int)
        // RDKit❗✔️: Atom::Atom(unsigned int num) : RDProps() {
        // RDKit❗✔️:   d_atomicNum = num;
        // RDKit❗✔️:   initAtom();
        // RDKit❗✔️: };
        // END RDKIT CPP FUNCTION Atom::Atom(unsigned int)
        // The typed identity preserves the source u8 value while canonicalizing
        // values with an Element; the predicate remains the caller's separate AST.
        Self::from_properties(
            id,
            identity,
            AtomProperties::new(),
            predicate,
            QueryPredicateOrigin::Explicit,
        )
    }

    #[doc(hidden)]
    pub fn with_identity(mut self, identity: QueryAtomIdentity) -> Self {
        self.identity = identity.canonicalized();
        self
    }

    #[doc(hidden)]
    pub fn with_id(mut self, id: AtomId) -> Self {
        self.id = id;
        self
    }

    fn from_properties(
        id: AtomId,
        identity: QueryAtomIdentity,
        properties: AtomProperties,
        predicate: QueryNode<AtomQueryPredicate>,
        predicate_origin: QueryPredicateOrigin,
    ) -> Self {
        Self {
            id,
            identity: identity.canonicalized(),
            properties,
            predicate,
            predicate_origin,
        }
    }

    #[must_use]
    pub const fn identity(&self) -> QueryAtomIdentity {
        self.identity
    }

    #[must_use]
    pub const fn atomic_number(&self) -> u8 {
        self.identity.atomic_number()
    }

    #[must_use]
    pub const fn element(&self) -> Option<Element> {
        self.identity.element()
    }

    /// Convert at an explicit Element-only boundary, preserving all common
    /// carrier state when conversion is representable.
    pub fn try_to_atom(&self) -> Result<Atom, QueryAtomConversionError> {
        let Some(element) = self.identity.element() else {
            return Err(QueryAtomConversionError::NonElementAtomicNumber {
                atom: self.id,
                atomic_number: self.atomic_number(),
            });
        };
        Ok(Atom::from_query_parts(
            self.id,
            element,
            self.properties.clone(),
        ))
    }

    #[must_use]
    pub fn predicate(&self) -> &QueryNode<AtomQueryPredicate> {
        &self.predicate
    }

    #[doc(hidden)]
    pub fn predicate_mut(&mut self) -> &mut QueryNode<AtomQueryPredicate> {
        self.predicate_origin = QueryPredicateOrigin::Explicit;
        &mut self.predicate
    }

    #[doc(hidden)]
    pub fn set_predicate(&mut self, predicate: QueryNode<AtomQueryPredicate>) {
        self.predicate = predicate;
        self.predicate_origin = QueryPredicateOrigin::Explicit;
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn predicate_is_carrier_derived(&self) -> bool {
        matches!(self.predicate_origin, QueryPredicateOrigin::CarrierDerived)
    }

    #[must_use]
    pub const fn index(&self) -> usize {
        self.id.index()
    }

    #[must_use]
    pub const fn id(&self) -> AtomId {
        self.id
    }

    #[must_use]
    pub const fn formal_charge(&self) -> i8 {
        self.properties.formal_charge
    }

    #[must_use]
    pub const fn explicit_hydrogens(&self) -> u8 {
        self.properties.explicit_hydrogens
    }

    #[must_use]
    pub const fn chiral_tag(&self) -> ChiralTag {
        self.properties.chiral_tag
    }

    #[must_use]
    pub const fn chiral_permutation(&self) -> Option<u32> {
        self.properties.chiral_permutation
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.properties.unknown_stereo
    }

    #[must_use]
    pub const fn mol_parity(&self) -> Option<i32> {
        self.properties.mol_parity
    }

    #[must_use]
    pub const fn mol_inversion_flag(&self) -> Option<i32> {
        self.properties.mol_inversion_flag
    }

    #[must_use]
    pub const fn implicit_hydrogen(&self) -> bool {
        self.properties.implicit_hydrogen
    }

    #[must_use]
    pub fn tracked_isotopic_hydrogens(&self) -> &[u16] {
        &self.properties.tracked_isotopic_hydrogens
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.properties.is_aromatic
    }

    #[must_use]
    pub const fn isotope(&self) -> Option<u16> {
        self.properties.isotope
    }

    #[must_use]
    pub const fn atom_map(&self) -> Option<u32> {
        self.properties.atom_map
    }

    #[must_use]
    pub const fn no_implicit(&self) -> bool {
        self.properties.no_implicit
    }

    #[must_use]
    pub const fn radical_electrons(&self) -> u8 {
        self.properties.radical_electrons
    }

    #[must_use]
    pub const fn hybridization(&self) -> Hybridization {
        self.properties.hybridization
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.properties.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.properties.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.properties.computed_props.contains(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        &self.properties.computed_props
    }

    #[must_use]
    pub const fn pdb_residue_info(&self) -> Option<&crate::AtomPdbResidueInfo> {
        self.properties.pdb_residue_info.as_ref()
    }

    #[must_use]
    pub const fn template_attachment_order(&self) -> Option<&TemplateAttachmentOrder> {
        self.properties.template_attachment_order.as_ref()
    }

    #[doc(hidden)]
    pub fn set_chiral_tag(&mut self, value: ChiralTag) {
        self.properties.chiral_tag = value;
    }

    #[doc(hidden)]
    pub fn set_chiral_permutation(&mut self, value: Option<u32>) {
        self.properties.chiral_permutation = value;
    }

    #[doc(hidden)]
    pub fn set_unknown_stereo(&mut self, value: bool) {
        self.properties.unknown_stereo = value;
    }

    #[doc(hidden)]
    pub fn set_mol_parity(&mut self, value: Option<i32>) {
        self.properties.mol_parity = value;
    }

    #[doc(hidden)]
    pub fn set_mol_inversion_flag(&mut self, value: Option<i32>) {
        self.properties.mol_inversion_flag = value;
    }

    #[doc(hidden)]
    pub fn set_implicit_hydrogen(&mut self, value: bool) {
        self.properties.implicit_hydrogen = value;
    }

    #[doc(hidden)]
    pub fn set_tracked_isotopic_hydrogens(&mut self, value: Vec<u16>) {
        self.properties.tracked_isotopic_hydrogens = value;
    }

    #[doc(hidden)]
    pub fn set_aromatic(&mut self, value: bool) {
        self.properties.is_aromatic = value;
    }

    #[doc(hidden)]
    pub fn set_formal_charge(&mut self, value: i8) {
        self.properties.formal_charge = value;
    }

    #[doc(hidden)]
    pub fn set_explicit_hydrogens(&mut self, value: u8) {
        self.properties.explicit_hydrogens = value;
    }

    #[doc(hidden)]
    pub fn set_isotope(&mut self, value: Option<u16>) {
        self.properties.isotope = value.filter(|value| *value != 0);
    }

    #[doc(hidden)]
    pub fn set_atom_map(&mut self, value: Option<u32>) {
        self.properties.atom_map = value;
    }

    #[doc(hidden)]
    pub fn set_no_implicit(&mut self, value: bool) {
        self.properties.no_implicit = value;
    }

    #[doc(hidden)]
    pub fn set_radical_electrons(&mut self, value: u8) {
        self.properties.radical_electrons = value;
    }

    #[doc(hidden)]
    pub fn set_hybridization(&mut self, value: Hybridization) {
        self.properties.hybridization = value;
    }

    #[doc(hidden)]
    pub fn set_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), AtomPropertyError> {
        self.properties.set_prop(key, value)
    }

    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), AtomPropertyError> {
        self.properties.set_computed_prop(key, value)
    }

    #[doc(hidden)]
    pub fn clear_prop(&mut self, key: &str) {
        self.properties.clear_prop(key);
    }

    #[doc(hidden)]
    pub fn clear_computed_props(&mut self) {
        self.properties.clear_computed_props();
    }

    #[doc(hidden)]
    pub fn set_pdb_residue_info(&mut self, value: Option<crate::AtomPdbResidueInfo>) {
        self.properties.pdb_residue_info = value;
    }

    #[doc(hidden)]
    pub fn set_template_attachment_order(&mut self, value: Option<TemplateAttachmentOrder>) {
        self.properties.template_attachment_order = value;
    }

    #[doc(hidden)]
    pub fn remap_template_attachment_order(
        &mut self,
        old_to_new: &[Option<AtomId>],
    ) -> Result<(), TemplateAttachmentOrderError> {
        self.properties.remap_template_attachment_order(old_to_new)
    }
}

/// A query bond combines concrete bond attributes with a query predicate tree.
/// Equality includes predicate origin, as for [`QueryAtom`]; it is not a
/// chemical or query-matching equivalence test.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct QueryBond {
    bond: Bond,
    predicate: QueryNode<BondQueryPredicate>,
    predicate_origin: QueryPredicateOrigin,
}

/// Structural errors while borrowing or remapping typed query rows alongside
/// a detached concrete topology.
#[doc(hidden)]
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryStateError {
    #[error(
        "query atom row {position} (id {atom}) has non-Element atomic number {atomic_number}, which cannot align with concrete topology"
    )]
    NonElementAtomIdentity {
        position: usize,
        atom: AtomId,
        atomic_number: u8,
    },
    #[error("query atom state has {actual} rows, expected {expected}")]
    AtomCount { actual: usize, expected: usize },
    #[error("query bond state has {actual} rows, expected {expected}")]
    BondCount { actual: usize, expected: usize },
    #[error("query atom row {position} has id {actual:?}, expected {expected:?}")]
    AtomId {
        position: usize,
        actual: AtomId,
        expected: AtomId,
    },
    #[error("query bond row {position} has id {actual:?}, expected {expected:?}")]
    BondId {
        position: usize,
        actual: BondId,
        expected: BondId,
    },
    #[error(
        "query bond row {position} has endpoints {actual:?}, expected topology endpoints {expected:?}"
    )]
    BondEndpoints {
        position: usize,
        actual: (AtomId, AtomId),
        expected: (AtomId, AtomId),
    },
    #[error("query-state mapping is invalid: {0}")]
    Mapping(#[from] MappingValidationError),
    #[error("query-state mapping appends {entity} row {position}, which has no source query row")]
    AppendedRow {
        entity: &'static str,
        position: usize,
    },
    #[error("query-state mapping supplied {actual} appended {entity} rows, expected {expected}")]
    AppendedRowCount {
        entity: &'static str,
        actual: usize,
        expected: usize,
    },
}

/// A validated, non-owning view of typed query identity and predicates aligned
/// one-for-one with a detached concrete topology.
///
/// This is an overlay, not a current carrier snapshot. Alignment covers row
/// counts, IDs and bond endpoints only. After chemistry transforms, carrier
/// attributes must be read from the current `TopologyBlock`, never from the
/// query rows retained here. Predicates and origins remain available without
/// exposing those potentially stale carriers.
///
/// ```compile_fail
/// fn cannot_read_old_atoms(state: cosmolkit_model::QueryStateRef<'_>) {
///     let _ = state.atoms();
/// }
/// ```
/// ```compile_fail
/// fn cannot_read_old_bonds(state: cosmolkit_model::QueryStateRef<'_>) {
///     let _ = state.bonds();
/// }
/// ```
#[doc(hidden)]
#[derive(Debug, Clone, Copy)]
pub struct QueryStateRef<'a> {
    atoms: &'a [QueryAtom],
    bonds: &'a [QueryBond],
}

impl<'a> QueryStateRef<'a> {
    pub fn try_for_topology(
        atoms: &'a [QueryAtom],
        bonds: &'a [QueryBond],
        topology: &TopologyBlock,
    ) -> Result<Self, QueryStateError> {
        if atoms.len() != topology.atoms.len() {
            return Err(QueryStateError::AtomCount {
                actual: atoms.len(),
                expected: topology.atoms.len(),
            });
        }
        if bonds.len() != topology.bonds.len() {
            return Err(QueryStateError::BondCount {
                actual: bonds.len(),
                expected: topology.bonds.len(),
            });
        }
        for (position, query) in atoms.iter().enumerate() {
            if let QueryAtomIdentity::AtomicNumber(atomic_number) = query.identity() {
                return Err(QueryStateError::NonElementAtomIdentity {
                    position,
                    atom: query.id(),
                    atomic_number,
                });
            }
        }
        for (position, (query, carrier)) in atoms.iter().zip(&topology.atoms).enumerate() {
            if query.id() != carrier.id() {
                return Err(QueryStateError::AtomId {
                    position,
                    actual: query.id(),
                    expected: carrier.id(),
                });
            }
        }
        for (position, (query, carrier)) in bonds.iter().zip(&topology.bonds).enumerate() {
            if query.id() != carrier.id() {
                return Err(QueryStateError::BondId {
                    position,
                    actual: query.id(),
                    expected: carrier.id(),
                });
            }
            let actual = (query.begin(), query.end());
            let expected = (carrier.begin(), carrier.end());
            if actual != expected {
                return Err(QueryStateError::BondEndpoints {
                    position,
                    actual,
                    expected,
                });
            }
        }
        Ok(Self { atoms, bonds })
    }

    /// Revalidate row alignment against the current topology without exposing
    /// carrier snapshots or requiring chemistry attributes to remain unchanged.
    pub fn validate_for_topology(self, topology: &TopologyBlock) -> Result<(), QueryStateError> {
        Self::try_for_topology(self.atoms, self.bonds, topology).map(|_| ())
    }

    #[must_use]
    pub fn atom_has_query(self, atom: AtomId) -> bool {
        // BEGIN RDKIT CPP FUNCTION QueryAtom::hasQuery
        // RDKit✔️✔️: // This method can be used to distinguish query atoms from standard atoms:
        // RDKit✔️✔️: bool hasQuery() const override { return dp_query != nullptr; }
        // END RDKIT CPP FUNCTION QueryAtom::hasQuery
        // Behavior review: Explicit rows model source QueryAtom values; carrier-derived
        // rows model ordinary Atom values lifted only for uniform Rust storage.
        // Complexity review: one checked slice lookup and one enum comparison are O(1),
        // matching the source null-pointer test without allocation or cloning.
        !self.atoms[atom.index()].predicate_is_carrier_derived()
    }

    #[must_use]
    pub fn bond_has_query(self, bond: BondId) -> bool {
        // BEGIN RDKIT CPP FUNCTION QueryBond::hasQuery
        // RDKit✔️✔️: // This method can be used to distinguish query bonds from standard bonds
        // RDKit✔️✔️: bool hasQuery() const override { return dp_query != nullptr; }
        // END RDKIT CPP FUNCTION QueryBond::hasQuery
        // Behavior review: Explicit and carrier-derived rows retain the same source
        // dynamic-type distinction as atom rows.
        // Complexity review: this is an allocation-free O(1) indexed test.
        !self.bonds[bond.index()].predicate_is_carrier_derived()
    }

    #[must_use]
    pub fn atom_predicate(self, atom: AtomId) -> &'a QueryNode<AtomQueryPredicate> {
        self.atoms[atom.index()].predicate()
    }

    #[must_use]
    pub fn bond_predicate(self, bond: BondId) -> &'a QueryNode<BondQueryPredicate> {
        self.bonds[bond.index()].predicate()
    }
}

/// Remap canonical typed query rows with the same validated topology mapping
/// used for every other atom- and bond-indexed detached block.
#[doc(hidden)]
pub fn remap_query_rows(
    state: QueryStateRef<'_>,
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
) -> Result<(Vec<QueryAtom>, Vec<QueryBond>), QueryStateError> {
    remap_query_rows_with_appended(state, topology, mapping, &[], &[])
}

/// Transport existing predicate/origin rows and explicitly supplied appended
/// rows onto the *current* topology carriers. Generic remapping never infers
/// what a new query means: the owning append algorithm must provide one row
/// for each `None` in the topology mapping, in new-row order.
#[doc(hidden)]
pub fn remap_query_rows_with_appended(
    state: QueryStateRef<'_>,
    topology: &TopologyBlock,
    mapping: &TopologyMapping,
    appended_atoms: &[QueryAtom],
    appended_bonds: &[QueryBond],
) -> Result<(Vec<QueryAtom>, Vec<QueryBond>), QueryStateError> {
    mapping.validate_for_counts(
        state.atoms.len(),
        topology.atoms.len(),
        state.bonds.len(),
        topology.bonds.len(),
    )?;

    let mut next_appended_atom = 0;
    let mut atoms = Vec::with_capacity(topology.atoms.len());
    for (position, (carrier, old)) in topology
        .atoms
        .iter()
        .zip(mapping.atoms().new_to_old())
        .enumerate()
    {
        let source = if let Some(old) = old {
            &state.atoms[old.index()]
        } else {
            let source =
                appended_atoms
                    .get(next_appended_atom)
                    .ok_or(QueryStateError::AppendedRow {
                        entity: "atom",
                        position,
                    })?;
            next_appended_atom += 1;
            if source.id() != carrier.id() {
                return Err(QueryStateError::AtomId {
                    position,
                    actual: source.id(),
                    expected: carrier.id(),
                });
            }
            if old.is_none()
                && let QueryAtomIdentity::AtomicNumber(atomic_number) = source.identity()
            {
                return Err(QueryStateError::NonElementAtomIdentity {
                    position,
                    atom: source.id(),
                    atomic_number,
                });
            }
            source
        };
        let mut mapped = QueryAtom::from_parts(carrier.clone(), source.predicate.clone());
        mapped.predicate_origin = source.predicate_origin;
        atoms.push(mapped);
    }
    if next_appended_atom != appended_atoms.len() {
        return Err(QueryStateError::AppendedRowCount {
            entity: "atom",
            actual: appended_atoms.len(),
            expected: next_appended_atom,
        });
    }

    let mut next_appended_bond = 0;
    let mut bonds = Vec::with_capacity(topology.bonds.len());
    for (position, (carrier, old)) in topology
        .bonds
        .iter()
        .zip(mapping.bonds().new_to_old())
        .enumerate()
    {
        let source = if let Some(old) = old {
            &state.bonds[old.index()]
        } else {
            let source =
                appended_bonds
                    .get(next_appended_bond)
                    .ok_or(QueryStateError::AppendedRow {
                        entity: "bond",
                        position,
                    })?;
            next_appended_bond += 1;
            if source.id() != carrier.id() {
                return Err(QueryStateError::BondId {
                    position,
                    actual: source.id(),
                    expected: carrier.id(),
                });
            }
            let actual = (source.begin(), source.end());
            let expected = (carrier.begin(), carrier.end());
            if actual != expected {
                return Err(QueryStateError::BondEndpoints {
                    position,
                    actual,
                    expected,
                });
            }
            source
        };
        bonds.push(QueryBond {
            bond: carrier.clone(),
            predicate: source.predicate.clone(),
            predicate_origin: source.predicate_origin,
        });
    }
    if next_appended_bond != appended_bonds.len() {
        return Err(QueryStateError::AppendedRowCount {
            entity: "bond",
            actual: appended_bonds.len(),
            expected: next_appended_bond,
        });
    }

    QueryStateRef::try_for_topology(&atoms, &bonds, topology)?;
    Ok((atoms, bonds))
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
        Self {
            bond,
            predicate,
            predicate_origin: QueryPredicateOrigin::Explicit,
        }
    }

    /// Construct the uniform query carrier used internally when Molfile input
    /// contains a mixture of ordinary and query bonds.
    #[doc(hidden)]
    #[must_use]
    pub fn from_carrier_parts(bond: Bond, predicate: QueryNode<BondQueryPredicate>) -> Self {
        Self {
            bond,
            predicate,
            predicate_origin: QueryPredicateOrigin::CarrierDerived,
        }
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
        self.predicate_origin = QueryPredicateOrigin::Explicit;
        &mut self.predicate
    }

    #[doc(hidden)]
    pub fn set_predicate(&mut self, predicate: QueryNode<BondQueryPredicate>) {
        self.predicate = predicate;
        self.predicate_origin = QueryPredicateOrigin::Explicit;
    }

    #[doc(hidden)]
    #[must_use]
    pub const fn predicate_is_carrier_derived(&self) -> bool {
        matches!(self.predicate_origin, QueryPredicateOrigin::CarrierDerived)
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
///
/// `PartialEq` compares stored representation, including atom/bond predicate
/// origins and metadata. It is not graph isomorphism or matching equivalence.
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
            if let Some(order) = atom.template_attachment_order() {
                order
                    .validate_for_atom_count(self.atoms.len())
                    .map_err(|source| QueryGraphError::TemplateAttachmentOrder {
                        atom: atom.id(),
                        source,
                    })?;
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

    /// Clone the complete detached coordinate carrier for a domain owner that
    /// must apply an index-changing transform and rebuild this query value.
    #[doc(hidden)]
    #[must_use]
    pub fn coordinate_block(
        &self,
        source_coordinate_dim: Option<crate::CoordinateDimension>,
    ) -> CoordinateBlock {
        CoordinateBlock {
            conformers_2d: self.conformers_2d.clone(),
            conformers_3d: self.conformers_3d.clone(),
            source_coordinate_dim,
        }
    }

    #[doc(hidden)]
    pub fn add_conformer_3d(&mut self, conformer: Conformer3D) -> Result<(), QueryGraphError> {
        if conformer.coordinates().len() != self.num_atoms() {
            return Err(QueryGraphError::CoordinateValidation(
                CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: conformer.id(),
                    rows: conformer.coordinates().len(),
                    atom_count: self.num_atoms(),
                },
            ));
        }
        self.conformers_3d.push(conformer);
        Ok(())
    }

    #[must_use]
    pub fn coordinates_2d(&self) -> Option<&[[f64; 2]]> {
        self.conformers_2d.first().map(Conformer2D::coordinates)
    }

    #[doc(hidden)]
    pub fn with_2d_coordinate_block(
        mut self,
        coords: Vec<[f64; 2]>,
    ) -> Result<Self, QueryGraphError> {
        if coords.len() != self.num_atoms() {
            return Err(QueryGraphError::CoordinateValidation(
                CoordinateValidationError::RowCount {
                    dimension: "2D",
                    conformer: 0,
                    rows: coords.len(),
                    atom_count: self.num_atoms(),
                },
            ));
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
    #[error("query atom {atom} has invalid template attachment order: {source}")]
    TemplateAttachmentOrder {
        atom: AtomId,
        source: TemplateAttachmentOrderError,
    },
    #[error("query bond at position {position} has id {id}, expected {position}")]
    BondIdMismatch { position: usize, id: BondId },
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
