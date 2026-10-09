//! Query values used by SMARTS, MCS, and substructure algorithms.
//!
//! This module contains only query data and local graph validation. SMARTS
//! parsing, writing, matching, serialization, and compilation belong uniquely
//! to `cosmolkit-search`; query data is never lowered back to a concrete
//! `Molecule`.

use crate::PropertyText;
use std::collections::{BTreeMap, BTreeSet};

use cosmolkit_types::{BondDirection, BondOrder, BondStereo, ChiralTag, Element, Hybridization};

use crate::atom::AtomProperties;
use crate::sgroup::{SubstanceGroupValidationError, validate_substance_groups};
use crate::{
    Atom, AtomId, AtomPropertyError, Bond, BondId, Conformer2D, Conformer3D, CoordinateBlock,
    CoordinateValidationError, MappingValidationError, StereoGroup, SubstanceGroup,
    TemplateAttachmentOrder, TemplateAttachmentOrderError, TopologyBlock, TopologyMapping,
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
    source_smarts: Option<PropertyText>,
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
    pub fn with_source_smarts(mut self, smarts: impl Into<PropertyText>) -> Self {
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
    pub fn source_smarts(&self) -> Option<&PropertyText> {
        self.source_smarts.as_ref()
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
        // RDKit❗✔️:   Queries::Query<int, Atom const *, true> *copy() const override {
        // RDKit❗✔️:     RecursiveStructureQuery *res = new RecursiveStructureQuery();
        // RDKit❗✔️:     res->dp_queryMol.reset(new ROMol(*dp_queryMol, true));
        // RDKit❗✔️:
        // RDKit❗✔️:     std::set<int>::const_iterator i;
        // RDKit❗✔️:     for (i = d_set.begin(); i != d_set.end(); i++) {
        // RDKit❗✔️:       res->insert(*i);
        // RDKit❗✔️:     }
        // RDKit❗✔️:     res->setNegation(getNegation());
        // RDKit❗✔️:     res->d_description = d_description;
        // RDKit❗✔️:     res->d_serialNumber = d_serialNumber;
        // RDKit❗✔️:     return res;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   unsigned int getSerialNumber() const { return d_serialNumber; }
        // RDKit❗✔️:
        // END RDKIT CPP FUNCTION RecursiveStructureQuery::copy
        // Native recursive copy calls the quick ROMol constructor. Keep the
        // modeled null query graph absent; native dereferences that pointer in
        // copy, an invalid-state boundary still requiring final reconciliation.
        // No complete graph clone before removing metadata or coordinates.
        Self {
            query_graph: self
                .query_graph
                .as_ref()
                .map(|graph| Box::new(graph.source_copy(true, -1))),
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
    /// Return transported Atom member effects without changing query identity,
    /// origin or predicate. This is detached common-field transport only.
    #[doc(hidden)]
    pub fn replace_source_carrier_members_from(&mut self, source: &Atom) {
        assert_eq!(self.id(), source.id(), "source carrier row identity");
        assert_eq!(
            self.atomic_number(),
            source.atomic_number(),
            "source carrier atomic identity"
        );
        self.properties = source.source_common_properties().clone();
    }

    /// Borrow property records using the source private/computed include flags.
    #[doc(hidden)]
    pub fn property_records(
        &self,
        include_private: bool,
        include_computed: bool,
    ) -> Result<impl Iterator<Item = (&PropertyText, &crate::PropertyValue)> + '_, AtomPropertyError>
    {
        self.properties
            .props
            .filtered_ordered(include_private, include_computed)
            .map_err(AtomPropertyError::from)
    }

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

    /// Borrow a query and its detached ordered dictionary as disjoint fields.
    #[doc(hidden)]
    pub fn predicate_and_properties_mut(
        &mut self,
    ) -> (&QueryNode<AtomQueryPredicate>, &mut crate::PropertyStore) {
        // RDKit❗✔️:   QUERYATOM_QUERY *getQuery() const override { return dp_query; }
        // RDKit❗✔️:   Dict &getDict() { return d_props; }
        // This only splits detached field borrows; it grants no live molecule,
        // query mutation or commit authority, and performs no clone/allocation.
        (&self.predicate, &mut self.properties.props)
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
    pub fn props(&self) -> &BTreeMap<PropertyText, crate::PropertyValue> {
        self.properties.props.values()
    }

    #[must_use]
    pub fn prop(&self, key: impl AsRef<[u8]>) -> Option<&crate::PropertyValue> {
        self.properties.props.get(key.as_ref())
    }

    /// Read a required detached property without conversion.
    #[doc(hidden)]
    pub fn prop_required(
        &self,
        key: impl AsRef<[u8]>,
    ) -> Result<&crate::PropertyValue, crate::MissingPropertyError> {
        self.properties.props.get_required(key)
    }

    #[must_use]
    pub fn is_prop_computed(
        &self,
        key: impl AsRef<[u8]>,
    ) -> Result<bool, crate::PropertyValueError> {
        self.properties.props.is_computed(key)
    }

    #[must_use]
    pub fn computed_prop_names(
        &self,
    ) -> Result<Option<&[PropertyText]>, crate::PropertyValueError> {
        self.properties.props.computed_names()
    }

    /// Copy the canonical source dictionary and its typed property carriers.
    /// Atom members (including temporary flags and cached valences), identity,
    /// predicate and monomer info remain those of the receiving atom.
    #[doc(hidden)]
    pub fn replace_source_properties_from(&mut self, source: &Self) {
        // BEGIN RDKIT CPP FUNCTION RDProps::updateProps
        // RDKit✔️❌: void updateProps(const RDProps &source, bool preserveExisting = false) {
        // RDKit✔️❌:     d_props.update(source.getDict(), preserveExisting);
        // RDKit✔️❌:   }
        // END RDKIT CPP FUNCTION RDProps::updateProps
        // RDKit✔️✔️: const Dict &getDict() const { return d_props; }
        // BEGIN RDKIT CPP REUSED HELPER Dict::update
        // RDKit✔️❌:   void update(const Dict &other, bool preserveExisting = false) {
        // RDKit✔️❌:     if (!preserveExisting) {
        // RDKit✔️❌:       *this = other;
        // RDKit✔️❌:     } else {
        // RDKit✔️❌:       if (other._hasNonPodData) {
        // RDKit✔️❌:         _hasNonPodData = true;
        // RDKit✔️❌:       }
        // RDKit✔️❌:       for (const auto &opair : other._data) {
        // RDKit✔️❌:         Pair *target = nullptr;
        // RDKit✔️❌:         for (auto &dpair : _data) {
        // RDKit✔️❌:           if (dpair.key == opair.key) {
        // RDKit✔️❌:             target = &dpair;
        // RDKit✔️❌:             break;
        // RDKit✔️❌:           }
        // RDKit✔️❌:         }
        // RDKit✔️❌:
        // RDKit✔️❌:         if (!target) {
        // RDKit✔️❌:           // need to create blank entry and copy
        // RDKit✔️❌:           _data.push_back(Pair(opair.key));
        // RDKit✔️❌:           copy_rdvalue(_data.back().val, opair.val);
        // RDKit✔️❌:         } else {
        // RDKit✔️❌:           // just copy
        // RDKit✔️❌:           copy_rdvalue(target->val, opair.val);
        // RDKit✔️❌:         }
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        // RDKit✔️❌:   }
        // END RDKIT CPP REUSED HELPER Dict::update
        // Behavior: this ReplaceAtom source call site uses preserveExisting=false.
        // Delegate the entire ordered byte dictionary (including its actual
        // __computedProps tagged entry) to PropertyStore::update_from, whose
        // independent SF373 comparison implements both source branches. Do not
        // inspect/cast that reserved value during a dictionary copy: wrong-kind
        // payloads are copied unchanged, just as source ordinary entries are.
        // Typed property carriers below are source dictionary facts. Receiving
        // Atom members, cached valences, predicate, identity and monomer info
        // are not assigned by RDProps::updateProps and remain untouched.
        // Complexity: one deep dictionary copy plus typed property payloads.
        // Source false-branch copying is linear in entries and bytes. The
        // canonical tree/order representation adds tree allocations and an
        // additional owning key index; record this explicit performance gap.
        self.properties
            .props
            .update_from(&source.properties.props, false);
        self.properties.chiral_permutation = source.properties.chiral_permutation;
        self.properties.unknown_stereo = source.properties.unknown_stereo;
        self.properties.mol_parity = source.properties.mol_parity;
        self.properties.mol_inversion_flag = source.properties.mol_inversion_flag;
        self.properties.implicit_hydrogen = source.properties.implicit_hydrogen;
        self.properties.tracked_isotopic_hydrogens =
            source.properties.tracked_isotopic_hydrogens.clone();
        self.properties.atom_map = source.properties.atom_map;
        self.properties.template_attachment_order =
            source.properties.template_attachment_order.clone();
    }

    /// Read detached source storage without executing a chemistry getter.
    #[doc(hidden)]
    pub const fn source_valence_facts(&self) -> crate::SourceAtomValenceFacts {
        self.properties.source_valence_facts
    }

    #[doc(hidden)]
    pub fn set_source_valence_facts(&mut self, facts: crate::SourceAtomValenceFacts) {
        self.properties.source_valence_facts = facts;
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
        key: impl Into<PropertyText>,
        value: impl Into<crate::PropertyValue>,
    ) -> Result<(), AtomPropertyError> {
        self.properties.set_prop(key, value)
    }

    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<crate::PropertyValue>,
    ) -> Result<(), AtomPropertyError> {
        self.properties.set_computed_prop(key, value)
    }

    #[doc(hidden)]
    pub fn clear_prop(&mut self, key: impl AsRef<[u8]>) -> Result<(), AtomPropertyError> {
        self.properties.clear_prop(key)
    }

    #[doc(hidden)]
    pub fn clear_computed_props(&mut self) -> Result<(), AtomPropertyError> {
        self.properties.clear_computed_props()
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
        let mut bond = carrier.clone();
        // Carrier attributes and source query rows may both carry the incoming
        // identity; retain only the established query-row root in this value.
        // A query-bearing carrier clone temporarily copies its tree before it
        // is dropped here, an additional O(query-tree) cost versus old carriers.
        let _ = bond.take_query();
        bonds.push(QueryBond {
            bond,
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
        // RDKit✔️✔️: QueryBond::QueryBond(BondType bT) : Bond(bT) {
        // RDKit✔️✔️:   if (bT != Bond::UNSPECIFIED) {
        // RDKit✔️✔️:     dp_query = makeBondOrderEqualsQuery(bT);
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     dp_query = makeBondNullQuery();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️: };
        let mut bond = Bond::from_spec(id, spec);
        let predicate = bond.take_query().unwrap_or_else(|| {
            if bond.order() == BondOrder::Unspecified {
                QueryNode::predicate(BondQueryPredicate::Any)
            } else {
                QueryNode::predicate(BondQueryPredicate::Order(bond.order()))
            }
        });
        Self::from_parts(bond, predicate)
    }

    #[must_use]
    pub fn from_parts(mut bond: Bond, predicate: QueryNode<BondQueryPredicate>) -> Self {
        // RDKit✔️✔️:   void setQuery(QUERYBOND_QUERY *what) override {
        // RDKit✔️✔️:     // free up any existing query (Issue255):
        // RDKit✔️✔️:     delete dp_query;
        // RDKit✔️✔️:     dp_query = what;
        // RDKit✔️✔️:   }
        // QueryGraph keeps its existing typed root as the sole source identity;
        // the detached carrier contains attributes, never a second query tree.
        let _ = bond.take_query();
        Self {
            bond,
            predicate,
            predicate_origin: QueryPredicateOrigin::Explicit,
        }
    }

    /// Construct a uniform carrier, preserving an incoming explicit query.
    #[doc(hidden)]
    #[must_use]
    pub fn from_carrier_parts(mut bond: Bond, predicate: QueryNode<BondQueryPredicate>) -> Self {
        // RDKit✔️✔️:   bool hasQuery() const override { return dp_query != nullptr; }
        // RDKit✔️✔️:   QUERYBOND_QUERY *getQuery() const override { return dp_query; }
        match bond.take_query() {
            Some(query) => Self {
                bond,
                predicate: query,
                predicate_origin: QueryPredicateOrigin::Explicit,
            },
            None => Self {
                bond,
                predicate,
                predicate_origin: QueryPredicateOrigin::CarrierDerived,
            },
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

    /// Move a detached query bond while retaining its sole predicate identity.
    #[doc(hidden)]
    pub fn remapped(
        self,
        id: BondId,
        begin: AtomId,
        end: AtomId,
        stereo_atoms: Option<[AtomId; 2]>,
    ) -> Self {
        Self {
            bond: self.bond.remapped(id, begin, end, stereo_atoms),
            ..self
        }
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
#[derive(Debug, PartialEq)]
pub struct QueryGraph {
    atoms: Vec<QueryAtom>,
    bonds: Vec<QueryBond>,
    adjacency: Vec<Vec<(usize, usize)>>,
    props: crate::property_value::PropertyStore,
    conformers_2d: Vec<Conformer2D>,
    conformers_3d: Vec<Conformer3D>,
    source_conformer_order: Option<Vec<crate::CoordinateDimension>>,
    stereo_groups: Vec<StereoGroup>,
    substance_groups: Vec<SubstanceGroup>,
    source_ring_info: crate::SourceRingInfo,
}

impl Clone for QueryGraph {
    fn clone(&self) -> Self {
        // Default ROMol copy constructor supplies quickCopy=false, confId=-1.
        // The unique copy implementation below also owns recursive quickCopy.
        self.source_copy(false, -1)
    }
}

impl QueryGraph {
    /// Borrow the actual source cache without finding or normalizing rings.
    #[doc(hidden)]
    pub fn source_ring_info(&self) -> &crate::SourceRingInfo {
        // RDKit❗✔️:   RingInfo *getRingInfo() const { return dp_ringInfo; }
        &self.source_ring_info
    }

    /// Transport an explicit owner-produced source cache into this detached graph.
    #[doc(hidden)]
    pub fn replace_source_ring_info(&mut self, source: crate::SourceRingInfo) {
        // RDKit❗✔️:   RingInfo &operator=(const RingInfo &other) = default;
        // This move is the typed equivalent of transferring the supplied cache;
        // no graph scan, finding, renumbering, validation repair or inference.
        self.source_ring_info = source;
    }

    /// Exact ordered dictionary projection, including raw computed-list data.
    #[doc(hidden)]
    pub fn source_molecule_properties(&self) -> crate::MoleculeProperties {
        crate::MoleculeProperties::from_source_dict(self.props.clone())
    }
    /// Return only source molecule property effects from a detached algorithm.
    #[doc(hidden)]
    pub fn replace_source_molecule_properties(&mut self, source: &crate::MoleculeProperties) {
        self.props = source.source_dict().clone();
    }

    // Source copy options are an internal model prerequisite, not a new
    // public molecule transform. RecursiveStructureQuery::copy uses quickCopy.
    fn source_copy(&self, quick_copy: bool, conf_id: i32) -> Self {
        // BEGIN RDKIT CPP FUNCTION ROMol::initFromOther
        // RDKit❗✔️: void ROMol::initFromOther(const ROMol &other, bool quickCopy, int confId) {
        // RDKit❗✔️:   if (this == &other) {
        // RDKit❗✔️:     return;
        // RDKit❗✔️:   }
        // RDKit❗✔️:   numBonds = 0;
        // RDKit❗✔️:   // std::cerr<<"    init from other: "<<this<<" "<<&other<<std::endl;
        // RDKit❗✔️:   // copy over the atoms
        // RDKit❗✔️:   for (const auto oatom : other.atoms()) {
        // RDKit❗✔️:     constexpr bool updateLabel = false;
        // RDKit❗✔️:     constexpr bool takeOwnership = true;
        // RDKit❗✔️:     addAtom(oatom->copy(), updateLabel, takeOwnership);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // and the bonds:
        // RDKit❗✔️:   for (const auto obond : other.bonds()) {
        // RDKit❗✔️:     addBond(obond->copy(), true);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // ring information
        // RDKit❗✔️:   delete dp_ringInfo;
        // RDKit❗✔️:   if (other.dp_ringInfo) {
        // RDKit❗✔️:     dp_ringInfo = new RingInfo(*(other.dp_ringInfo));
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     dp_ringInfo = new RingInfo();
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // enhanced stereochemical information
        // RDKit❗✔️:   d_stereo_groups.clear();
        // RDKit❗✔️:   for (auto &otherGroup : other.d_stereo_groups) {
        // RDKit❗✔️:     std::vector<Atom *> atoms;
        // RDKit❗✔️:     for (auto &otherAtom : otherGroup.getAtoms()) {
        // RDKit❗✔️:       atoms.push_back(getAtomWithIdx(otherAtom->getIdx()));
        // RDKit❗✔️:     }
        // RDKit❗✔️:     std::vector<Bond *> bonds;
        // RDKit❗✔️:     for (auto &otherBond : otherGroup.getBonds()) {
        // RDKit❗✔️:       bonds.push_back(getBondWithIdx(otherBond->getIdx()));
        // RDKit❗✔️:     }
        // RDKit❗✔️:     d_stereo_groups.emplace_back(otherGroup.getGroupType(), std::move(atoms),
        // RDKit❗✔️:                                  std::move(bonds), otherGroup.getReadId());
        // RDKit❗✔️:     d_stereo_groups.back().setWriteId(otherGroup.getWriteId());
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (other.dp_delAtoms) {
        // RDKit❗✔️:     dp_delAtoms.reset(new boost::dynamic_bitset<>(*other.dp_delAtoms));
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     dp_delAtoms.reset(nullptr);
        // RDKit❗✔️:   }
        // RDKit❗✔️:   if (other.dp_delBonds) {
        // RDKit❗✔️:     dp_delBonds.reset(new boost::dynamic_bitset<>(*other.dp_delBonds));
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     dp_delBonds.reset(nullptr);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   if (!quickCopy) {
        // RDKit❗✔️:     // copy conformations
        // RDKit❗✔️:     for (const auto &conf : other.d_confs) {
        // RDKit❗✔️:       if (confId < 0 || rdcast<int>(conf->getId()) == confId) {
        // RDKit❗✔️:         this->addConformer(new Conformer(*conf));
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     // Copy sgroups
        // RDKit❗✔️:     for (const auto &sg : getSubstanceGroups(other)) {
        // RDKit❗✔️:       addSubstanceGroup(*this, sg);
        // RDKit❗✔️:     }
        // RDKit❗✔️:
        // RDKit❗✔️:     d_props = other.d_props;
        // RDKit❗✔️:
        // RDKit❗✔️:     // Bookmarks should be copied as well:
        // RDKit❗✔️:     for (auto abmI : other.d_atomBookmarks) {
        // RDKit❗✔️:       for (const auto *aptr : abmI.second) {
        // RDKit❗✔️:         setAtomBookmark(getAtomWithIdx(aptr->getIdx()), abmI.first);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:     for (auto bbmI : other.d_bondBookmarks) {
        // RDKit❗✔️:       for (const auto *bptr : bbmI.second) {
        // RDKit❗✔️:         setBondBookmark(getBondWithIdx(bptr->getIdx()), bbmI.first);
        // RDKit❗✔️:       }
        // RDKit❗✔️:     }
        // RDKit❗✔️:   } else {
        // RDKit❗✔️:     d_props.reset();
        // RDKit❗✔️:     STR_VECT computed;
        // RDKit❗✔️:     d_props.setVal(RDKit::detail::computedPropName, computed);
        // RDKit❗✔️:   }
        // RDKit❗✔️:
        // RDKit❗✔️:   // std::cerr<<"---------    done init from other: "<<this<<"
        // RDKit❗✔️:   // "<<&other<<std::endl;
        // RDKit❗✔️: }
        // END RDKIT CPP FUNCTION ROMol::initFromOther
        // Behavioral boundary: a newly returned detached graph has stable dense
        // IDs instead of native owning pointers; self-aliasing assignment is
        // unrepresentable. QueryGraph retains the represented RingInfo state; pending deletion masks and
        // bookmark storage, so those independent source states remain unmodeled.
        // Native conformer IDs are unsigned32; for nonnegative conf_id equality
        // with that value is equivalent to release rdcast<int>(id)==confId.
        // Wide model IDs, native debug cast checks and missing mixed-dimension
        // append order are retained representation gaps, not fabricated state.
        // Cost: linear deep copy of represented atoms, bonds, query trees,
        // properties/groups and selected coordinate rows, like native copying.
        // The existing adjacency cache and explicit conformer-order transport
        // add O(V+E+C) scalar copies; no whole-graph clone followed by erasure,
        // rank sorting, property conversion, validation or unselected coords copy.
        let atoms = self.atoms.clone();
        let bonds = self.bonds.clone();
        let adjacency = self.adjacency.clone();
        let stereo_groups = self.stereo_groups.clone();
        let (conformers_2d, conformers_3d, source_conformer_order, substance_groups, props) =
            if quick_copy {
                (
                    Vec::new(),
                    Vec::new(),
                    None,
                    Vec::new(),
                    crate::property_value::PropertyStore::from_records([(
                        PropertyText::from("__computedProps"),
                        crate::PropertyValue::StringVector(Vec::new()),
                    )]),
                )
            } else {
                let selected = |id: usize| conf_id < 0 || id == conf_id as usize;
                let conformers_2d = self
                    .conformers_2d
                    .iter()
                    .filter(|conf| selected(conf.id()))
                    .cloned()
                    .collect();
                let conformers_3d = self
                    .conformers_3d
                    .iter()
                    .filter(|conf| selected(conf.id()))
                    .cloned()
                    .collect();
                let source_conformer_order = self.source_conformer_order.as_ref().map(|order| {
                    let mut two_d = self.conformers_2d.iter();
                    let mut three_d = self.conformers_3d.iter();
                    order
                        .iter()
                        .copied()
                        .filter(|dimension| {
                            let id = match dimension {
                                crate::CoordinateDimension::TwoD => two_d
                                    .next()
                                    .expect("stored conformer order indexes its 2D collection")
                                    .id(),
                                crate::CoordinateDimension::ThreeD => three_d
                                    .next()
                                    .expect("stored conformer order indexes its 3D collection")
                                    .id(),
                            };
                            selected(id)
                        })
                        .collect()
                });
                (
                    conformers_2d,
                    conformers_3d,
                    source_conformer_order,
                    self.substance_groups.clone(),
                    self.props.clone(),
                )
            };
        Self {
            atoms,
            bonds,
            adjacency,
            props,
            conformers_2d,
            conformers_3d,
            source_conformer_order,
            stereo_groups,
            substance_groups,
            source_ring_info: self.source_ring_info.clone(),
        }
    }

    /// Borrow property records using the source private/computed include flags.
    #[doc(hidden)]
    pub fn property_records(
        &self,
        include_private: bool,
        include_computed: bool,
    ) -> Result<
        impl Iterator<Item = (&PropertyText, &crate::PropertyValue)> + '_,
        crate::MoleculePropertyError,
    > {
        self.props
            .filtered_ordered(include_private, include_computed)
            .map_err(crate::MoleculePropertyError::from)
    }

    #[must_use]
    pub fn from_parts(
        atoms: Vec<QueryAtom>,
        bonds: Vec<QueryBond>,
        props: impl IntoIterator<Item = (PropertyText, crate::PropertyValue)>,
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
            props: crate::property_value::PropertyStore::from_records(props),
            conformers_2d,
            conformers_3d,
            source_conformer_order: None,
            stereo_groups,
            substance_groups: Vec::new(),
            source_ring_info: crate::SourceRingInfo::default(),
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
            source_conformer_order: self.source_conformer_order.clone(),
        };
        coordinates
            .validate_for_atom_count(self.atoms.len())
            .map_err(QueryGraphError::CoordinateValidation)?;

        validate_substance_groups(&self.substance_groups, self.atoms.len(), self.bonds.len())?;

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
            group.validate_members()?;
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

    /// Source degree of an atom owned by this detached graph, if its row exists.
    #[doc(hidden)]
    pub fn try_atom_degree(&self, atom: AtomId) -> Option<u32> {
        // RDKit❗✔️: unsigned int Atom::getDegree() const {
        // RDKit❗✔️:   return dp_mol ? getOwningMol().getAtomDegree(this) : 0;
        // RDKit❗✔️: }
        // RDKit❗✔️: unsigned int ROMol::getAtomDegree(const Atom *at) const {
        // RDKit❗✔️:   PRECONDITION(at, "no atom");
        // RDKit❗✔️:   PRECONDITION(&at->getOwningMol() == this,
        // RDKit❗✔️:                "atom not associated with this molecule");
        // RDKit❗✔️:   return rdcast<unsigned int>(boost::out_degree(at->getIdx(), d_graph));
        // RDKit❗✔️: };
        // A QueryGraph always represents its attached atoms. Missing row data
        // stays absent; callers propagate their own structural bounds context.
        self.atoms.get(atom.index())?;
        self.adjacency
            .get(atom.index())
            .map(|neighbors| neighbors.len() as u32)
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
    pub fn props(&self) -> &BTreeMap<PropertyText, crate::PropertyValue> {
        self.props.values()
    }

    #[must_use]
    /// Borrow the same dictionary's records in actual source insertion order.
    #[doc(hidden)]
    pub fn ordered_props(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyText, &crate::PropertyValue)> + '_ {
        self.props.ordered()
    }

    pub fn prop(&self, key: impl AsRef<[u8]>) -> Option<&crate::PropertyValue> {
        self.props.get(key.as_ref())
    }

    /// Read a required detached property without conversion.
    #[doc(hidden)]
    pub fn prop_required(
        &self,
        key: impl AsRef<[u8]>,
    ) -> Result<&crate::PropertyValue, crate::MissingPropertyError> {
        self.props.get_required(key)
    }

    /// Remove a parser/runtime property from this detached query value.
    #[doc(hidden)]
    pub fn clear_prop(
        &mut self,
        key: impl AsRef<[u8]>,
    ) -> Result<(), crate::MoleculePropertyError> {
        self.props
            .clear(key)
            .map_err(crate::MoleculePropertyError::from)
    }

    /// Borrow the source's ordered computed-list fact, including presence.
    #[doc(hidden)]
    pub fn computed_prop_names(
        &self,
    ) -> Result<Option<&[PropertyText]>, crate::PropertyValueError> {
        self.props.computed_names()
    }

    /// Store an ordinary detached value and the source computed membership.
    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<crate::PropertyValue>,
    ) -> Result<(), crate::MoleculePropertyError> {
        self.props
            .set_computed(key.into(), value.into())
            .map_err(crate::MoleculePropertyError::from)
    }

    #[must_use]
    pub fn name(&self) -> Result<Option<&PropertyText>, crate::PropertyValueError> {
        self.prop("_Name")
            .map(crate::PropertyValue::as_string)
            .transpose()
    }

    #[doc(hidden)]
    pub fn with_name(mut self, name: impl Into<PropertyText>) -> Self {
        self.props
            .set("_Name".into(), crate::PropertyValue::String(name.into()))
            .expect("fixed nonempty ordinary property key");
        self
    }

    #[doc(hidden)]
    pub fn with_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<crate::PropertyValue>,
    ) -> Result<Self, crate::MoleculePropertyError> {
        self.set_prop(key, value)?;
        Ok(self)
    }

    #[doc(hidden)]
    pub fn set_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<crate::PropertyValue>,
    ) -> Result<(), crate::MoleculePropertyError> {
        self.props
            .set(key.into(), value.into())
            .map_err(crate::MoleculePropertyError::from)
    }

    /// Borrow canonical 2D conformers without cloning the coordinate carrier.
    #[doc(hidden)]
    #[must_use]
    pub fn conformers_2d(&self) -> &[Conformer2D] {
        &self.conformers_2d
    }

    #[must_use]
    pub fn conformers_3d(&self) -> &[Conformer3D] {
        &self.conformers_3d
    }

    /// Borrow the actual source-front conformer through the sole MODEL selector.
    #[doc(hidden)]
    pub fn first_source_conformer(
        &self,
    ) -> Result<Option<crate::CoordinateSourceConformer<'_>>, crate::CoordinateValidationError>
    {
        crate::coordinates::first_source_conformer_from_parts(
            &self.conformers_2d,
            &self.conformers_3d,
            self.source_conformer_order.as_deref(),
        )
    }

    /// Borrow every actual source conformer in physical append order.
    #[doc(hidden)]
    pub fn source_conformers(
        &self,
    ) -> Result<Vec<crate::CoordinateSourceConformer<'_>>, crate::CoordinateValidationError> {
        // RDKit❗❌:   inline ConstConformerIterator beginConformers() const {
        // RDKit❗❌:     return d_confs.begin();
        // RDKit❗❌:   }
        // RDKit❗❌:   inline ConstConformerIterator endConformers() const { return d_confs.end(); }
        // These are actual detached storage-order facts, not conformer ID,
        // is3D, dimensional preference, import provenance or numerical checks.
        // Cost ❌: the split dimensional carrier materializes C borrowed refs
        // versus Native's O(1) begin/end iterator pair; no coordinates clone.
        use crate::{
            CoordinateDimension as Dim, CoordinateSourceConformer as Row,
            CoordinateValidationError as Error,
        };
        let mut result = Vec::with_capacity(self.conformers_2d.len() + self.conformers_3d.len());
        if let Some(order) = self.source_conformer_order.as_deref() {
            let (mut two, mut three) = (0, 0);
            for dim in order {
                match dim {
                    Dim::TwoD => {
                        let row = self
                            .conformers_2d
                            .get(two)
                            .ok_or(Error::MissingSourceConformerOrder)?;
                        result.push(Row::TwoD(row));
                        two += 1;
                    }
                    Dim::ThreeD => {
                        let row = self
                            .conformers_3d
                            .get(three)
                            .ok_or(Error::MissingSourceConformerOrder)?;
                        result.push(Row::ThreeD(row));
                        three += 1;
                    }
                }
            }
            if two != self.conformers_2d.len() || three != self.conformers_3d.len() {
                return Err(Error::SourceConformerOrder {
                    two_d: two,
                    three_d: three,
                    expected_two_d: self.conformers_2d.len(),
                    expected_three_d: self.conformers_3d.len(),
                });
            }
        } else if self.conformers_2d.is_empty() {
            result.extend(self.conformers_3d.iter().map(Row::ThreeD));
        } else if self.conformers_3d.is_empty() {
            result.extend(self.conformers_2d.iter().map(Row::TwoD));
        } else {
            return Err(Error::MissingSourceConformerOrder);
        }
        Ok(result)
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
            source_conformer_order: self.source_conformer_order.clone(),
        }
    }

    #[doc(hidden)]
    pub fn set_source_conformer_order(
        &mut self,
        order: Option<Vec<crate::CoordinateDimension>>,
    ) -> Result<(), QueryGraphError> {
        if let Some(order) = &order {
            let two_d = order
                .iter()
                .filter(|d| **d == crate::CoordinateDimension::TwoD)
                .count();
            let three_d = order.len() - two_d;
            if two_d != self.conformers_2d.len() || three_d != self.conformers_3d.len() {
                return Err(QueryGraphError::CoordinateValidation(
                    CoordinateValidationError::SourceConformerOrder {
                        two_d,
                        three_d,
                        expected_two_d: self.conformers_2d.len(),
                        expected_three_d: self.conformers_3d.len(),
                    },
                ));
            }
        }
        self.source_conformer_order = order;
        Ok(())
    }

    #[doc(hidden)]
    #[must_use]
    pub fn source_conformer_order(&self) -> Option<&[crate::CoordinateDimension]> {
        self.source_conformer_order.as_deref()
    }

    #[doc(hidden)]
    pub fn add_conformer_3d(&mut self, conformer: Conformer3D) -> Result<(), QueryGraphError> {
        // RDKit❗❌: unsigned int ROMol::addConformer(Conformer *conf, bool assignId) {
        // RDKit❗❌:   PRECONDITION(conf, "bad conformer");
        // RDKit❗❌:   PRECONDITION(conf->getNumAtoms() == this->getNumAtoms(),
        // RDKit❗❌:                "Number of atom mismatch");
        // RDKit❗❌:   if (assignId) {
        // RDKit❗❌:     int maxId = -1;
        // RDKit❗❌:     for (auto cptr : d_confs) {
        // RDKit❗❌:       maxId = std::max((int)(cptr->getId()), maxId);
        // RDKit❗❌:     }
        // RDKit❗❌:     maxId++;
        // RDKit❗❌:     conf->setId((unsigned int)maxId);
        // RDKit❗❌:   }
        // RDKit❗❌:   conf->setOwningMol(this);
        // RDKit❗❌:   CONFORMER_SPTR nConf(conf);
        // RDKit❗❌:   d_confs.push_back(nConf);
        // RDKit❗❌:   return conf->getId();
        // RDKit❗❌: }
        // This canonical detached boundary specializes assignId=false.
        // Split dimensional storage additionally records actual append order.

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
        crate::coordinates::record_source_append(
            &mut self.source_conformer_order,
            self.conformers_2d.len(),
            self.conformers_3d.len(),
            crate::CoordinateDimension::ThreeD,
        )
        .map_err(QueryGraphError::CoordinateValidation)?;
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
        if let Some(order) = &mut self.source_conformer_order {
            order.retain(|dimension| *dimension != crate::CoordinateDimension::TwoD);
        }
        crate::coordinates::record_source_append(
            &mut self.source_conformer_order,
            0,
            self.conformers_3d.len(),
            crate::CoordinateDimension::TwoD,
        )
        .map_err(QueryGraphError::CoordinateValidation)?;
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

    /// Borrow existing detached groups without validation, merging or reordering.
    #[doc(hidden)]
    pub fn stereo_groups_mut(&mut self) -> &mut [StereoGroup] {
        &mut self.stereo_groups
    }
}

/// Replace a detached query graph's enhanced-stereo groups after validating
/// every atom and bond reference against the graph's current rows.
pub fn replace_query_stereo_groups(
    graph: &mut QueryGraph,
    groups: Vec<StereoGroup>,
) -> Result<(), QueryGraphError> {
    for group in &groups {
        for atom in group.atoms() {
            if atom.index() >= graph.atoms.len() {
                return Err(QueryGraphError::StereoGroupAtomOutOfRange {
                    atom: *atom,
                    atom_count: graph.atoms.len(),
                });
            }
        }
        for bond in group.bonds() {
            if bond.index() >= graph.bonds.len() {
                return Err(QueryGraphError::StereoGroupBondOutOfRange {
                    bond: *bond,
                    bond_count: graph.bonds.len(),
                });
            }
        }
        group.validate_members()?;
    }

    let checked = crate::merge_absolute_stereo_groups(groups)?;
    graph.stereo_groups = checked;
    Ok(())
}

/// Borrow the ordered typed substance groups owned by a detached query graph.
#[must_use]
pub fn query_substance_groups(graph: &QueryGraph) -> &[SubstanceGroup] {
    // RDKit✔️✔️: std::vector<SubstanceGroup> &getSubstanceGroups(ROMol &mol) {
    // RDKit✔️✔️:   return mol.d_sgroups;
    // RDKit✔️✔️: }
    // RDKit✔️✔️: const std::vector<SubstanceGroup> &getSubstanceGroups(const ROMol &mol) {
    // RDKit✔️✔️:   return mol.d_sgroups;
    // RDKit✔️✔️: }
    // Behavior: borrow the actual ordered collection, including every typed
    // payload and sparse role record. Read consumers observe its identity and
    // order without sorting, filtering, normalizing, or reconstructing groups.
    // Detached mutation uses the existing validated replacement boundary;
    // this read adapter does not expose a second mutable storage authority.
    // Complexity: O(1) borrow, no allocation or group/property clone.
    &graph.substance_groups
}

/// Replace only the ordered typed substance groups owned by a detached query graph.
pub fn replace_query_substance_groups(
    graph: &mut QueryGraph,
    groups: Vec<SubstanceGroup>,
) -> Result<(), QueryGraphError> {
    // RDKit✔️❌: unsigned int addSubstanceGroup(ROMol &mol, SubstanceGroup sgroup) {
    // RDKit✔️❌:   sgroup.setOwningMol(&mol);
    // RDKit✔️❌:
    // RDKit✔️❌:   auto &&sgroups = getSubstanceGroups(mol);
    // RDKit✔️❌:   unsigned int id = sgroups.size();
    // RDKit✔️❌:
    // RDKit✔️❌:   sgroups.push_back(std::move(sgroup));
    // RDKit✔️❌:
    // RDKit✔️❌:   return id;
    // RDKit✔️❌: }
    // RDKit✔️✔️: void SubstanceGroup::setOwningMol(ROMol *mol) {
    // RDKit✔️✔️:   PRECONDITION(mol, "owning molecule is nullptr");
    // RDKit✔️✔️:
    // RDKit✔️✔️:   dp_mol = mol;
    // RDKit✔️✔️: }
    // Behavior: the existing detached append callers assign the pre-append
    // collection length as dense identity, append the complete group, then
    // transfer that ordered collection here. Installation preserves every
    // payload; ownership is the graph containment, with no ROMol pointer.
    // The source helper does not set an "index" property: that write remains
    // in each source caller. Validate detached references before installation.
    // Complexity: validation scans group references and existing append
    // callers copy the prior SGroup collection, materially more work than
    // source amortized O(1) append. No molecule-wide clone is performed.
    validate_substance_groups(&groups, graph.atoms.len(), graph.bonds.len())?;
    graph.substance_groups = groups;
    Ok(())
}

/// Query graph construction failed because a graph-local constraint was invalid.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum QueryGraphError {
    #[error("{0}")]
    StereoGroup(#[from] crate::StereoGroupError),

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
    #[error("query substance group at position {position} has id {id:?}, expected {position}")]
    SubstanceGroupIdMismatch {
        position: usize,
        id: crate::SubstanceGroupId,
    },
    #[error(
        "query substance group {sgroup:?} references atom {atom}, out of range for {atom_count} atoms"
    )]
    SubstanceGroupAtomOutOfRange {
        sgroup: crate::SubstanceGroupId,
        atom: AtomId,
        atom_count: usize,
    },
    #[error(
        "query substance group {sgroup:?} references bond {bond}, out of range for {bond_count} bonds"
    )]
    SubstanceGroupBondOutOfRange {
        sgroup: crate::SubstanceGroupId,
        bond: BondId,
        bond_count: usize,
    },
    #[error("query substance group {sgroup:?} has parent {parent:?} out of range")]
    SubstanceGroupParentOutOfRange {
        sgroup: crate::SubstanceGroupId,
        parent: crate::SubstanceGroupId,
    },
    #[error("query graph coordinate validation failed: {0}")]
    CoordinateValidation(CoordinateValidationError),
    #[error("query graph adjacency does not match its atom and bond rows")]
    AdjacencyMismatch,
}

impl From<SubstanceGroupValidationError> for QueryGraphError {
    fn from(error: SubstanceGroupValidationError) -> Self {
        match error {
            SubstanceGroupValidationError::IdMismatch { position, id } => {
                Self::SubstanceGroupIdMismatch { position, id }
            }
            SubstanceGroupValidationError::AtomOutOfRange {
                sgroup,
                atom,
                atom_count,
            } => Self::SubstanceGroupAtomOutOfRange {
                sgroup,
                atom,
                atom_count,
            },
            SubstanceGroupValidationError::BondOutOfRange {
                sgroup,
                bond,
                bond_count,
            } => Self::SubstanceGroupBondOutOfRange {
                sgroup,
                bond,
                bond_count,
            },
            SubstanceGroupValidationError::ParentOutOfRange { sgroup, parent } => {
                Self::SubstanceGroupParentOutOfRange { sgroup, parent }
            }
        }
    }
}

/// Borrow query-atom properties in the same source insertion order as atom properties.
pub fn ordered_query_atom_properties(
    atom: &QueryAtom,
) -> impl ExactSizeIterator<Item = (&PropertyText, &crate::PropertyValue)> + '_ {
    // RDKit✔️✔️: for (const auto &item : _data) {
    // RDKit✔️✔️:   res.push_back(item.key);
    // RDKit✔️✔️: }
    // Dict::keys order is shared by ordinary and query atoms. Borrow the one
    // canonical store; do not clone a carrier or build a second property map.
    atom.properties.props.ordered()
}

#[cfg(test)]
mod tests {
    // These original fixtures contain UTF-8 literals. Decode only their
    // borrowed test projection; raw-byte controls assert as_bytes directly.
    // Invalid bytes fail this assertion, never change a chemistry outcome.
    fn fixture_text(value: &crate::PropertyText) -> &str {
        std::str::from_utf8(value.as_bytes()).expect("unchanged UTF-8 fixture bytes")
    }

    use super::*;
    use crate::{
        AtomSpec, BondSpec, Element, PropertyValue, SGroupAttachPoint, SGroupBondRole,
        SGroupBracket, SGroupBracketStyle, SGroupCState, SGroupConnection, SGroupData,
        SGroupDisplay, StereoGroupKind, SubstanceGroupId, SubstanceGroupKind,
    };

    fn carbon(id: usize) -> QueryAtom {
        QueryAtom::new(AtomId::new(id), AtomSpec::new(Element::C))
    }

    fn stereo_test_graph() -> QueryGraph {
        let mut atom0 = QueryAtom::from_identity_parts(
            AtomId::new(0),
            QueryAtomIdentity::Element(Element::C),
            QueryNode::predicate(AtomQueryPredicate::FormalCharge(-1)),
        );
        atom0
            .set_prop("atom-note", "preserved")
            .expect("valid atom property");
        let atom1 = QueryAtom::from_carrier_parts(
            Atom::from_spec(AtomId::new(1), AtomSpec::new(Element::N)),
            QueryNode::predicate(AtomQueryPredicate::Any),
        );
        let bond = QueryBond::from_carrier_parts(
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
            ),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );

        QueryGraph::from_parts(
            vec![atom0, atom1],
            vec![bond],
            BTreeMap::from([("graph-note".into(), crate::PropertyValue::from("preserved"))]),
            vec![
                Conformer2D::new(9, vec![[0.0, 1.0], [2.0, 3.0]]).with_prop("frame", "first"),
                Conformer2D::new(11, vec![[4.0, 5.0], [6.0, 7.0]]).with_prop("frame", "second"),
            ],
            vec![
                Conformer3D::new(17, vec![[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]], true)
                    .with_prop("frame", "three-dimensional"),
                Conformer3D::new(19, vec![[6.0, 7.0, 8.0], [9.0, 10.0, 11.0]], false)
                    .with_prop("frame", "two-dimensional"),
            ],
            vec![
                StereoGroup::new(
                    StereoGroupKind::And,
                    vec![AtomId::new(0)],
                    vec![BondId::new(0)],
                )
                .expect("valid distinct stereo members")
                .with_id(7),
            ],
        )
        .expect("valid stereo test graph")
    }

    #[test]
    fn query_sgroups_empty_on_new_graph() {
        let graph = QueryGraph::from_parts(
            Vec::new(),
            Vec::new(),
            BTreeMap::new(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
        )
        .expect("empty detached query graph is valid");

        assert!(query_substance_groups(&graph).is_empty());
    }

    #[test]
    fn query_sgroups_ordered_replace_preserves_typed_and_graph_state() {
        let mut graph = stereo_test_graph();
        let original = graph.clone();
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(7)
                .with_external_id(31)
                .with_atoms(vec![AtomId::new(0), AtomId::new(0)])
                .with_label("ordered-first"),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Superatom)
                .with_parent(SubstanceGroupId::new(0))
                .with_atoms(vec![AtomId::new(1)])
                .with_label("ordered-second"),
        ];

        replace_query_substance_groups(&mut graph, groups.clone())
            .expect("typed query SGroups are replaceable");

        assert_eq!(query_substance_groups(&graph), groups);
        assert_eq!(graph.clone(), graph);
        assert_eq!(graph.atoms, original.atoms);
        assert_eq!(graph.bonds, original.bonds);
        assert_eq!(graph.adjacency, original.adjacency);
        assert_eq!(graph.props, original.props);
        assert_eq!(graph.stereo_groups, original.stereo_groups);
        assert_eq!(graph.conformers_2d, original.conformers_2d);
        assert_eq!(graph.conformers_3d, original.conformers_3d);
        assert_eq!(graph.atoms[0].predicate(), original.atoms[0].predicate());
        assert_eq!(
            graph.atoms[0].predicate_is_carrier_derived(),
            original.atoms[0].predicate_is_carrier_derived()
        );
        assert_eq!(graph.atoms[1].predicate(), original.atoms[1].predicate());
        assert_eq!(
            graph.atoms[1].predicate_is_carrier_derived(),
            original.atoms[1].predicate_is_carrier_derived()
        );
        assert_eq!(graph.bonds[0].predicate(), original.bonds[0].predicate());
        assert_eq!(
            graph.bonds[0].predicate_is_carrier_derived(),
            original.bonds[0].predicate_is_carrier_derived()
        );
    }

    #[test]
    fn query_sgroups_explicit_transport_survives_graph_reconstruction() {
        let mut original = stereo_test_graph();
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_rdkit_sequence_id(23)
                .with_external_id(41)
                .with_atoms(vec![AtomId::new(1), AtomId::new(0), AtomId::new(1)])
                .with_label("transported query data"),
        ];
        replace_query_substance_groups(&mut original, groups.clone())
            .expect("source query SGroups are structurally valid");

        // Rebuilding an existing graph keeps from_parts' new-graph semantics,
        // then explicitly transports the typed groups into the reconstruction.
        let mut rebuilt = QueryGraph::from_parts(
            original.atoms.clone(),
            original.bonds.clone(),
            original
                .ordered_props()
                .map(|(key, value)| (key.clone(), value.clone())),
            original.conformers_2d.clone(),
            original.conformers_3d.clone(),
            original.stereo_groups.clone(),
        )
        .expect("existing non-SGroup query state reconstructs");
        replace_query_substance_groups(&mut rebuilt, query_substance_groups(&original).to_vec())
            .expect("typed SGroups explicitly transport through reconstruction");

        assert_eq!(query_substance_groups(&rebuilt), groups);
        assert_eq!(rebuilt, original);
        assert_eq!(rebuilt.atoms, original.atoms);
        assert_eq!(rebuilt.bonds, original.bonds);
        assert_eq!(rebuilt.adjacency, original.adjacency);
        assert_eq!(rebuilt.props, original.props);
        assert_eq!(rebuilt.stereo_groups, original.stereo_groups);
        assert_eq!(rebuilt.conformers_2d, original.conformers_2d);
        assert_eq!(rebuilt.conformers_3d, original.conformers_3d);
        assert_eq!(
            rebuilt.atoms[0].predicate_is_carrier_derived(),
            original.atoms[0].predicate_is_carrier_derived()
        );
        assert_eq!(
            rebuilt.atoms[1].predicate_is_carrier_derived(),
            original.atoms[1].predicate_is_carrier_derived()
        );
        assert_eq!(
            rebuilt.bonds[0].predicate_is_carrier_derived(),
            original.bonds[0].predicate_is_carrier_derived()
        );
    }

    #[test]
    fn query_sgroups_preserve_dat_polymer_metadata_and_ordered_references() {
        let mut graph = stereo_test_graph();
        let original_atoms = graph.atoms.clone();
        let original_bonds = graph.bonds.clone();
        let original_adjacency = graph.adjacency.clone();
        let original_props = graph.props.clone();
        let original_stereo_groups = graph.stereo_groups.clone();
        let original_conformers_2d = graph.conformers_2d.clone();
        let original_conformers_3d = graph.conformers_3d.clone();
        let bracket = SGroupBracket::new([[1.25, 2.5, 3.75], [4.5, 5.25, 6.75], [7.0, 8.5, 9.25]]);
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_external_id(40)
                .with_atoms(vec![AtomId::new(0), AtomId::new(0)])
                .with_data(SGroupData {
                    field_name: Some("FIELD".into()),
                    field_type: Some("T".into()),
                    field_info: Some("typed data".into()),
                    field_display: Some("display specification".into()),
                    units: Some("ppm".into()),
                    query_type: Some("Q".into()),
                    query_op: Some("OP".into()),
                    values: vec!["first".into(), "second".into()],
                })
                .with_data_field("field row one")
                .with_data_field("field row two")
                .with_prop("custom", "retained")
                .unwrap(),
            SubstanceGroup::new(
                SubstanceGroupId::new(1),
                SubstanceGroupKind::StructuralRepeatUnit,
            )
            .with_rdkit_sequence_id(73)
            .with_external_id(91)
            .with_parent(SubstanceGroupId::new(0))
            .with_atoms(vec![AtomId::new(1), AtomId::new(0), AtomId::new(1)])
            .with_parent_atoms(vec![AtomId::new(0), AtomId::new(0)])
            .with_bonds(vec![BondId::new(0)])
            .with_bond_role(BondId::new(0), SGroupBondRole::Contained)
            .with_head_crossing_bonds(vec![BondId::new(0), BondId::new(0)])
            .with_crossing_bond_correspondence(vec![BondId::new(0), BondId::new(0)])
            .with_cstates(vec![
                SGroupCState::new(BondId::new(0), [0.25, 0.5, 0.75]),
                SGroupCState::new(BondId::new(0), [1.25, 1.5, 1.75]),
            ])
            .with_display(SGroupDisplay {
                brackets: vec![bracket],
                field_position: Some([10.5, 11.75]),
                display_tag: Some("polymer-display".into()),
            })
            .with_bracket_style(SGroupBracketStyle::Parenthesis)
            .with_connection(SGroupConnection::HeadToTail)
            .with_label("repeat unit")
            .with_subtype("SRU")
            .with_expansion_state("expanded")
            .with_class("polymer-class")
            .with_component_number(6)
            .with_data_field("polymer field"),
        ];

        replace_query_substance_groups(&mut graph, groups.clone())
            .expect("DAT and polymer group state is structurally valid");

        assert_eq!(query_substance_groups(&graph), groups);
        assert_eq!(
            query_substance_groups(&graph)[0].kind(),
            &SubstanceGroupKind::Data
        );
        assert_eq!(
            query_substance_groups(&graph)[1].kind(),
            &SubstanceGroupKind::StructuralRepeatUnit
        );
        assert_eq!(
            query_substance_groups(&graph)[1].parent(),
            Some(SubstanceGroupId::new(0))
        );
        assert_eq!(
            query_substance_groups(&graph)[1].atoms(),
            &[AtomId::new(1), AtomId::new(0), AtomId::new(1)]
        );
        assert_eq!(
            query_substance_groups(&graph)[1].crossing_bond_correspondence(),
            &[BondId::new(0), BondId::new(0)]
        );
        assert_eq!(
            query_substance_groups(&graph)[1].head_crossing_bonds(),
            &[BondId::new(0), BondId::new(0)]
        );
        assert_eq!(
            query_substance_groups(&graph)[1].cstates(),
            &[
                SGroupCState::new(BondId::new(0), [0.25, 0.5, 0.75]),
                SGroupCState::new(BondId::new(0), [1.25, 1.5, 1.75]),
            ]
        );
        assert_eq!(
            query_substance_groups(&graph)[1]
                .display()
                .expect("typed polymer display")
                .brackets,
            vec![bracket]
        );
        assert_eq!(
            query_substance_groups(&graph)[0]
                .data()
                .unwrap()
                .values
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            ["first", "second"]
        );
        assert_eq!(
            query_substance_groups(&graph)[0]
                .data_fields()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["field row one", "field row two"]
        );
        assert_eq!(
            query_substance_groups(&graph)[1]
                .data_fields()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["polymer field"]
        );
        assert_eq!(graph.atoms, original_atoms);
        assert_eq!(graph.bonds, original_bonds);
        assert_eq!(graph.adjacency, original_adjacency);
        assert_eq!(graph.props, original_props);
        assert_eq!(graph.stereo_groups, original_stereo_groups);
        assert_eq!(graph.conformers_2d, original_conformers_2d);
        assert_eq!(graph.conformers_3d, original_conformers_3d);
        assert_eq!(graph.conformers_2d.len(), 2);
        assert_eq!(graph.conformers_3d.len(), 2);
        assert_eq!(graph.conformers_3d[0].id(), 17);
        assert!(graph.conformers_3d[0].is_3d());
        assert_eq!(graph.conformers_3d[1].id(), 19);
        assert!(!graph.conformers_3d[1].is_3d());
    }

    fn assert_query_sgroup_replacement_is_atomic(
        invalid: SubstanceGroup,
        expected: QueryGraphError,
    ) {
        let mut graph = stereo_test_graph();
        let existing = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_atoms(vec![AtomId::new(0)]),
        ];
        replace_query_substance_groups(&mut graph, existing).expect("valid existing SGroup state");
        let original = graph.clone();

        assert_eq!(
            replace_query_substance_groups(&mut graph, vec![invalid]),
            Err(expected)
        );
        assert_eq!(graph, original);
    }

    #[test]
    fn query_sgroups_replacement_rejects_each_invalid_reference_atomically() {
        let id = SubstanceGroupId::new(0);
        let invalid_atom = AtomId::new(2);
        let invalid_bond = BondId::new(1);

        assert_query_sgroup_replacement_is_atomic(
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data),
            QueryGraphError::SubstanceGroupIdMismatch {
                position: 0,
                id: SubstanceGroupId::new(1),
            },
        );
        for invalid in [
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_atoms(vec![invalid_atom]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_parent_atoms(vec![invalid_atom]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_attach_points(vec![
                SGroupAttachPoint {
                    atom: invalid_atom,
                    leaving_atom: None,
                    label: None,
                    order: None,
                },
            ]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_attach_points(vec![
                SGroupAttachPoint {
                    atom: AtomId::new(0),
                    leaving_atom: Some(invalid_atom),
                    label: None,
                    order: None,
                },
            ]),
        ] {
            assert_query_sgroup_replacement_is_atomic(
                invalid,
                QueryGraphError::SubstanceGroupAtomOutOfRange {
                    sgroup: id,
                    atom: invalid_atom,
                    atom_count: 2,
                },
            );
        }
        for invalid in [
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_bonds(vec![invalid_bond]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_cstates(vec![SGroupCState::new(invalid_bond, [1.0, 2.0, 3.0])]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_head_crossing_bonds(vec![invalid_bond]),
            SubstanceGroup::new(id, SubstanceGroupKind::Data)
                .with_crossing_bond_correspondence(vec![invalid_bond]),
        ] {
            assert_query_sgroup_replacement_is_atomic(
                invalid,
                QueryGraphError::SubstanceGroupBondOutOfRange {
                    sgroup: id,
                    bond: invalid_bond,
                    bond_count: 1,
                },
            );
        }
        assert_query_sgroup_replacement_is_atomic(
            SubstanceGroup::new(id, SubstanceGroupKind::Data).with_parent(SubstanceGroupId::new(1)),
            QueryGraphError::SubstanceGroupParentOutOfRange {
                sgroup: id,
                parent: SubstanceGroupId::new(1),
            },
        );
    }

    #[test]
    fn query_sgroups_valid_parent_cycles_are_not_rejected() {
        let groups = vec![
            SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                .with_parent(SubstanceGroupId::new(1)),
            SubstanceGroup::new(SubstanceGroupId::new(1), SubstanceGroupKind::Data)
                .with_parent(SubstanceGroupId::new(0)),
        ];
        let mut graph = stereo_test_graph();

        replace_query_substance_groups(&mut graph, groups.clone())
            .expect("in-range parent references are structurally valid");

        assert_eq!(query_substance_groups(&graph), groups);
        graph
            .validate()
            .expect("parent cycles are not a model rule");
    }

    #[test]
    fn query_stereo_replace_preserves_order_duplicates_and_graph_state() {
        let mut graph = stereo_test_graph();
        let original = graph.clone();
        let replacement = vec![
            StereoGroup::new(
                StereoGroupKind::Or,
                vec![AtomId::new(1), AtomId::new(0)],
                vec![BondId::new(0)],
            )
            .expect("valid distinct stereo members")
            .with_id(12),
            StereoGroup::new(
                StereoGroupKind::Absolute,
                vec![AtomId::new(0), AtomId::new(1)],
                Vec::new(),
            )
            .expect("valid distinct stereo members"),
        ];

        replace_query_stereo_groups(&mut graph, replacement.clone())
            .expect("valid ordered stereo replacement");

        assert_eq!(graph.stereo_groups, replacement);
        assert_eq!(graph.atoms, original.atoms);
        assert_eq!(graph.bonds, original.bonds);
        assert_eq!(graph.adjacency, original.adjacency);
        assert_eq!(graph.props, original.props);
        assert_eq!(graph.conformers_2d, original.conformers_2d);
        assert_eq!(graph.conformers_3d, original.conformers_3d);
        assert_eq!(graph.atoms[0].predicate(), original.atoms[0].predicate());
        assert!(!graph.atoms[0].predicate_is_carrier_derived());
        assert_eq!(graph.atoms[1].predicate(), original.atoms[1].predicate());
        assert!(graph.atoms[1].predicate_is_carrier_derived());
        assert!(graph.bonds[0].predicate_is_carrier_derived());
        assert_eq!(graph.conformers_2d[0].id(), 9);
        assert_eq!(
            graph.conformers_2d[0]
                .props()
                .get(b"frame".as_slice())
                .map(fixture_text)
                .unwrap(),
            "first"
        );
        assert_eq!(graph.conformers_2d[1].id(), 11);
        assert_eq!(graph.conformers_3d[0].id(), 17);
        assert!(graph.conformers_3d[0].is_3d());
        assert_eq!(
            graph.conformers_3d[0]
                .props()
                .get(b"frame".as_slice())
                .map(fixture_text)
                .unwrap(),
            "three-dimensional"
        );
        assert_eq!(graph.conformers_3d[1].id(), 19);
        assert!(!graph.conformers_3d[1].is_3d());
        assert_eq!(
            graph.conformers_3d[1]
                .props()
                .get(b"frame".as_slice())
                .map(fixture_text)
                .unwrap(),
            "two-dimensional"
        );
    }

    #[test]
    fn query_stereo_replace_accepts_empty_replacement() {
        let mut graph = stereo_test_graph();
        let original = graph.clone();

        replace_query_stereo_groups(&mut graph, Vec::new()).expect("empty replacement");

        assert!(graph.stereo_groups.is_empty());
        assert_eq!(graph.atoms, original.atoms);
        assert_eq!(graph.bonds, original.bonds);
        assert_eq!(graph.adjacency, original.adjacency);
        assert_eq!(graph.props, original.props);
        assert_eq!(graph.conformers_2d, original.conformers_2d);
        assert_eq!(graph.conformers_3d, original.conformers_3d);
    }

    #[test]
    fn query_stereo_replace_rejects_invalid_atom_without_mutation() {
        let mut graph = stereo_test_graph();
        let original = graph.clone();

        let result = replace_query_stereo_groups(
            &mut graph,
            vec![
                StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(2)], Vec::new())
                    .expect("valid distinct stereo members"),
            ],
        );

        assert_eq!(
            result,
            Err(QueryGraphError::StereoGroupAtomOutOfRange {
                atom: AtomId::new(2),
                atom_count: 2,
            })
        );
        assert_eq!(graph, original);
    }

    #[test]
    fn query_stereo_replace_rejects_invalid_bond_without_mutation() {
        let mut graph = stereo_test_graph();
        let original = graph.clone();

        let result = replace_query_stereo_groups(
            &mut graph,
            vec![
                StereoGroup::new(StereoGroupKind::Absolute, Vec::new(), vec![BondId::new(1)])
                    .expect("valid distinct stereo members"),
            ],
        );

        assert_eq!(
            result,
            Err(QueryGraphError::StereoGroupBondOutOfRange {
                bond: BondId::new(1),
                bond_count: 1,
            })
        );
        assert_eq!(graph, original);
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
            vec![
                StereoGroup::new(StereoGroupKind::Absolute, vec![AtomId::new(1)], Vec::new())
                    .expect("valid distinct stereo members"),
            ],
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
            props: crate::property_value::PropertyStore::default(),
            conformers_2d: Vec::new(),
            conformers_3d: Vec::new(),
            source_conformer_order: None,
            stereo_groups: Vec::new(),
            substance_groups: Vec::new(),
            source_ring_info: crate::SourceRingInfo::default(),
        };

        assert_eq!(graph.validate(), Err(QueryGraphError::AdjacencyMismatch));
    }

    #[test]
    fn typed_property_transport_query_row_remap_and_deletion_preserve_types_and_inputs() {
        let atoms = vec![
            Atom::from_spec(
                AtomId::new(0),
                AtomSpec::new(Element::C)
                    .with_prop("removed", PropertyValue::String("atom-0".into()))
                    .unwrap(),
            ),
            Atom::from_spec(
                AtomId::new(1),
                AtomSpec::new(Element::N)
                    .with_prop("first", PropertyValue::Int(7))
                    .unwrap()
                    .with_computed_prop("computed", PropertyValue::Double(-0.0))
                    .unwrap()
                    .with_prop("first", PropertyValue::Bool(true))
                    .unwrap(),
            ),
            Atom::from_spec(
                AtomId::new(2),
                AtomSpec::new(Element::O)
                    .with_prop("last", PropertyValue::Bool(false))
                    .unwrap(),
            ),
        ];
        let bonds = vec![
            Bond::from_spec(
                BondId::new(0),
                BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
                    .with_prop("removed", PropertyValue::Int(1))
                    .unwrap(),
            ),
            Bond::from_spec(
                BondId::new(1),
                BondSpec::new(AtomId::new(1), AtomId::new(2), BondOrder::Double)
                    .with_prop("bond", PropertyValue::Double(2.5))
                    .unwrap()
                    .with_computed_prop("bond-computed", PropertyValue::Bool(true))
                    .unwrap(),
            ),
        ];
        let topology = TopologyBlock::try_from_parts(atoms, bonds, Vec::new(), Vec::new()).unwrap();
        let query_atoms = topology
            .atoms
            .iter()
            .map(|atom| {
                QueryAtom::from_parts(
                    atom.clone(),
                    QueryNode::predicate(AtomQueryPredicate::AtomicNumber(atom.atomic_number())),
                )
            })
            .collect::<Vec<_>>();
        let query_bonds = topology
            .bonds
            .iter()
            .map(|bond| {
                QueryBond::from_parts(
                    bond.clone(),
                    QueryNode::predicate(BondQueryPredicate::Order(bond.order())),
                )
            })
            .collect::<Vec<_>>();
        assert_eq!(
            query_atoms[1].prop("computed"),
            Some(&PropertyValue::Double(-0.0))
        );
        assert!(query_atoms[1].is_prop_computed("computed").unwrap());
        assert_eq!(
            query_bonds[1].bond().prop("bond-computed"),
            Some(&PropertyValue::Bool(true))
        );
        assert!(
            query_bonds[1]
                .bond()
                .is_prop_computed("bond-computed")
                .unwrap()
        );
        let topology_before = topology.clone();
        let query_atoms_before = query_atoms.clone();
        let query_bonds_before = query_bonds.clone();

        let state = QueryStateRef::try_for_topology(&query_atoms, &query_bonds, &topology).unwrap();
        let mut edit = topology.begin_batch_edit().unwrap();
        edit.remove_atom(AtomId::new(0)).unwrap();
        let (edited, mapping) = edit.finish().unwrap();
        let (mapped_atoms, mapped_bonds) = remap_query_rows(state, &edited, &mapping).unwrap();

        assert_eq!(topology, topology_before);
        assert_eq!(query_atoms, query_atoms_before);
        assert_eq!(query_bonds, query_bonds_before);
        // RWMol::commitBatchEdit calls clearComputedProps(true) after removals.
        assert_eq!(mapped_atoms.len(), 2);
        assert_eq!(mapped_bonds.len(), 1);
        assert_eq!(
            mapped_atoms[0].prop("first"),
            Some(&PropertyValue::Bool(true))
        );
        assert_eq!(mapped_atoms[0].prop("computed"), None);
        assert!(!mapped_atoms[0].is_prop_computed("computed").unwrap());
        assert_eq!(
            mapped_atoms[0]
                .properties
                .props
                .ordered_keys()
                .iter()
                .map(fixture_text)
                .collect::<Vec<_>>(),
            &["first", "__computedProps"]
        );
        assert_eq!(
            mapped_atoms[1].prop("last"),
            Some(&PropertyValue::Bool(false))
        );
        assert_eq!(
            mapped_bonds[0].bond().prop("bond"),
            Some(&PropertyValue::Double(2.5))
        );
        assert_eq!(mapped_bonds[0].bond().prop("bond-computed"), None);
        assert!(
            !mapped_bonds[0]
                .bond()
                .is_prop_computed("bond-computed")
                .unwrap()
        );
        assert!(
            mapped_atoms
                .iter()
                .all(|atom| atom.prop("removed").is_none())
        );
        assert!(
            mapped_bonds
                .iter()
                .all(|bond| bond.bond().prop("removed").is_none())
        );
    }
    #[test]
    fn cf3d_flags_atom_query_carrier_roundtrip_preserves_temporary_word() {
        let mut carrier = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
        carrier.set_temporary_flags(u64::MAX);
        carrier
            .set_computed_prop("derived", "value")
            .expect("valid computed property");

        let predicate = QueryNode::predicate(AtomQueryPredicate::Any);
        let explicit = QueryAtom::from_parts(carrier.clone(), predicate.clone());
        let explicit_roundtrip = explicit
            .clone()
            .with_id(AtomId::new(3))
            .try_to_atom()
            .expect("element query converts to atom");
        assert_eq!(explicit_roundtrip.id(), AtomId::new(3));
        assert_eq!(explicit_roundtrip.temporary_flags(), u64::MAX);

        let mut derived = QueryAtom::from_carrier_parts(carrier, predicate);
        derived
            .clear_computed_props()
            .expect("original fixture property clear succeeds");
        let derived_roundtrip = derived
            .try_to_atom()
            .expect("carrier-derived query converts to atom");
        assert_eq!(derived_roundtrip.temporary_flags(), u64::MAX);
        assert_eq!(derived_roundtrip.prop("derived"), None);

        let different_flags = QueryAtom::from_parts(
            {
                let mut atom = Atom::from_spec(AtomId::new(0), AtomSpec::new(Element::C));
                atom.set_temporary_flags(1);
                atom
            },
            QueryNode::predicate(AtomQueryPredicate::Any),
        );
        assert_ne!(
            explicit, different_flags,
            "equality remains representation-based"
        );
    }

    #[test]
    fn cf3d_flags_bond_query_carrier_roundtrip_preserves_temporary_word() {
        let mut carrier = Bond::from_spec(
            BondId::new(0),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        );
        carrier.set_temporary_flags(u64::MAX);
        carrier
            .set_computed_prop("derived", "value")
            .expect("valid computed property");

        let predicate = QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single));
        let explicit = QueryBond::from_parts(carrier.clone(), predicate.clone());
        assert_eq!(explicit.bond().temporary_flags(), u64::MAX);
        assert_eq!(explicit.clone().bond().temporary_flags(), u64::MAX);

        let mut derived = QueryBond::from_carrier_parts(carrier, predicate);
        derived
            .bond_mut()
            .clear_computed_props()
            .expect("original fixture property clear succeeds");
        assert_eq!(derived.bond().temporary_flags(), u64::MAX);
        assert_eq!(derived.bond().prop("derived"), None);
    }
}

#[cfg(test)]
mod source_romol_copy_complete_tests {
    use super::*;
    use crate::{
        AtomSpec, CoordinateDimension, PropertyValue, StereoGroupKind, SubstanceGroupId,
        SubstanceGroupKind,
    };

    fn fixture() -> QueryGraph {
        let atoms = vec![QueryAtom::new(AtomId::new(0), AtomSpec::new(Element::C))];
        let mut graph = QueryGraph::from_parts(
            atoms,
            vec![],
            [
                (PropertyText::from("first"), PropertyValue::UInt(u32::MAX)),
                (
                    PropertyText::from_bytes(b"raw\0\xff"),
                    PropertyValue::String(PropertyText::from_bytes(b"bytes\xff\0")),
                ),
            ],
            vec![
                Conformer2D::new(7, vec![[1.0, 2.0]]),
                Conformer2D::new(8, vec![[3.0, 4.0]]),
            ],
            vec![
                Conformer3D::new(7, vec![[5.0, 6.0, 7.0]], false),
                Conformer3D::new(u32::MAX as usize, vec![[8.0, 9.0, 10.0]], true),
            ],
            vec![
                StereoGroup::new(StereoGroupKind::And, vec![AtomId::new(0)], vec![])
                    .expect("valid distinct stereo members")
                    .with_id(19)
                    .with_write_id(27),
            ],
        )
        .unwrap();
        graph.atoms[0]
            .set_prop("atom-note", PropertyValue::UInt(13))
            .unwrap();
        graph
            .set_source_conformer_order(Some(vec![
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
                CoordinateDimension::ThreeD,
                CoordinateDimension::TwoD,
            ]))
            .unwrap();
        replace_query_substance_groups(
            &mut graph,
            vec![
                SubstanceGroup::new(SubstanceGroupId::new(0), SubstanceGroupKind::Data)
                    .with_atoms(vec![AtomId::new(0)])
                    .with_label("group"),
            ],
        )
        .unwrap();
        graph
    }

    #[test]
    fn default_full_copy_preserves_ordered_payloads_and_detaches_storage() {
        let graph = fixture();
        let mut copy = graph.clone();
        assert_eq!(copy, graph);
        assert_eq!(
            copy.property_records(true, true)
                .unwrap()
                .map(|(key, _)| key.as_bytes().to_vec())
                .collect::<Vec<_>>(),
            [b"first".to_vec(), b"raw\0\xff".to_vec()]
        );
        assert_eq!(copy.stereo_groups[0].id(), Some(19));
        assert_eq!(copy.stereo_groups[0].write_id(), 27);
        assert_ne!(copy.atoms.as_ptr(), graph.atoms.as_ptr());
        assert_ne!(
            copy.conformers_3d[0].coordinates().as_ptr(),
            graph.conformers_3d[0].coordinates().as_ptr()
        );
        copy.atoms[0]
            .set_prop("atom-note", PropertyValue::UInt(99))
            .unwrap();
        copy.conformers_3d[0].coordinates_mut()[0][0] = 123.0;
        assert_eq!(
            graph.atoms[0].prop("atom-note"),
            Some(&PropertyValue::UInt(13))
        );
        assert_eq!(graph.conformers_3d[0].coordinates()[0][0], 5.0);
    }

    #[test]
    fn quick_copy_retains_structure_and_enhanced_groups_but_resets_source_metadata() {
        let graph = fixture();
        let quick = graph.source_copy(true, 7);
        assert_eq!(quick.atoms, graph.atoms);
        assert_eq!(quick.bonds, graph.bonds);
        assert_eq!(quick.adjacency, graph.adjacency);
        assert_eq!(quick.stereo_groups, graph.stereo_groups);
        assert!(quick.conformers_2d.is_empty());
        assert!(quick.conformers_3d.is_empty());
        assert!(quick.substance_groups.is_empty());
        assert_eq!(quick.source_conformer_order, None);
        assert_eq!(quick.props.values().len(), 1);
        assert_eq!(
            quick.prop("__computedProps"),
            Some(&PropertyValue::StringVector(vec![]))
        );
        assert!(quick.prop("first").is_none());
        assert!(quick.prop(b"raw\0\xff").is_none());
        assert_eq!(
            quick.atoms[0].prop("atom-note"),
            Some(&PropertyValue::UInt(13))
        );
        assert_eq!(graph.conformers_2d.len(), 2);
        assert_eq!(graph.substance_groups.len(), 1);
    }

    #[test]
    fn selected_conformer_copy_keeps_all_matching_ids_and_native_append_order() {
        let graph = fixture();
        let selected = graph.source_copy(false, 7);
        assert_eq!(selected.conformers_2d, vec![graph.conformers_2d[0].clone()]);
        assert_eq!(selected.conformers_3d, vec![graph.conformers_3d[0].clone()]);
        assert_eq!(
            selected.source_conformer_order,
            Some(vec![CoordinateDimension::ThreeD, CoordinateDimension::TwoD])
        );
        assert_eq!(selected.props, graph.props);
        assert_eq!(selected.substance_groups, graph.substance_groups);
        assert!(!selected.conformers_3d[0].is_3d());
        for id in [0, 99, i32::MAX] {
            let absent = graph.source_copy(false, id);
            assert!(absent.conformers_2d.is_empty());
            assert!(absent.conformers_3d.is_empty());
            assert_eq!(absent.source_conformer_order, Some(vec![]));
        }
        for id in [-1, i32::MIN] {
            assert_eq!(graph.source_copy(false, id), graph);
        }
    }

    #[test]
    fn recursive_copy_reaches_quick_graph_constructor_and_preserves_source_set() {
        let mut query = RecursiveStructureQuery::from_query_graph(fixture(), u32::MAX);
        query.insert_atom_index(-1);
        query.insert_atom_index(0);
        let mut copied = query.clone();
        let inner = copied.query_graph().unwrap();
        assert_eq!(inner, &query.query_graph().unwrap().source_copy(true, -1));
        assert!(inner.conformers_2d.is_empty());
        assert!(inner.conformers_3d.is_empty());
        assert!(inner.substance_groups.is_empty());
        assert_eq!(
            inner.prop("__computedProps"),
            Some(&PropertyValue::StringVector(vec![]))
        );
        assert_eq!(copied.serial_number(), u32::MAX);
        assert!(copied.contains_atom_index(-1));
        assert!(copied.contains_atom_index(0));
        copied.query_graph_mut().unwrap().atoms[0]
            .set_prop("atom-note", PropertyValue::UInt(99))
            .unwrap();
        assert_eq!(
            query.query_graph().unwrap().atoms[0].prop("atom-note"),
            Some(&PropertyValue::UInt(13))
        );
        assert_eq!(query.query_graph().unwrap().conformers_2d.len(), 2);
    }
}

#[cfg(test)]
mod source_ring_info_tests {
    use super::*;
    use crate::SourceRingInfo;
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            (0..3)
                .map(|i| {
                    QueryAtom::from_identity_parts(
                        AtomId::new(i),
                        QueryAtomIdentity::from_atomic_number(if i % 2 == 0 { 0 } else { 119 }),
                        QueryNode::predicate(AtomQueryPredicate::AtomicNumber(6)),
                    )
                })
                .collect(),
            [(0, 1), (1, 2), (2, 0)]
                .into_iter()
                .enumerate()
                .map(|(i, (a, b))| {
                    QueryBond::new(
                        BondId::new(i),
                        crate::BondSpec::new(AtomId::new(a), AtomId::new(b), BondOrder::Single),
                    )
                })
                .collect(),
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    fn cache() -> SourceRingInfo {
        SourceRingInfo {
            initialized: true,
            find_type: 2,
            atom_members: vec![vec![0]; 3],
            bond_members: vec![vec![0]; 3],
            atom_rings: vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]],
            bond_rings: vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]],
            atom_ring_families: vec![vec![AtomId::new(0), AtomId::new(1), AtomId::new(2)]],
            bond_ring_families: vec![vec![BondId::new(0), BondId::new(1), BondId::new(2)]],
            relevant_cycle_count: Some(1),
            fused_rings: vec![vec![false]],
            num_fused_bonds: vec![0],
        }
    }
    #[test]
    fn newly_constructed_cyclic_query_keeps_native_uninitialized_ring_state() {
        let q = graph();
        assert_eq!(q.source_ring_info(), &SourceRingInfo::default());
        assert!(!q.source_ring_info().initialized);
        assert_eq!(q.source_ring_info().find_type, 3);
        assert!(q.source_ring_info().atom_members.is_empty());
    }
    #[test]
    fn ordinary_and_quick_source_copies_retain_all_actual_cache_fields() {
        let mut q = graph();
        q.set_prop("private-work", "value").unwrap();
        q.replace_source_ring_info(cache());
        let ordinary = q.source_copy(false, -1);
        let quick = q.source_copy(true, -1);
        assert_eq!(ordinary.source_ring_info(), &cache());
        assert_eq!(quick.source_ring_info(), &cache());
        assert!(ordinary.prop("private-work").is_some());
        assert!(quick.prop("private-work").is_none());
    }
    #[test]
    fn source_cache_storage_detaches_without_reconstructing_or_repairing_memberships() {
        let mut q = graph();
        q.replace_source_ring_info(cache());
        let before = q.clone();
        let mut copy = q.clone();
        let mut state = copy.source_ring_info().clone();
        state.atom_members.push(vec![]);
        state.bond_members.push(vec![]);
        state.relevant_cycle_count = Some(7);
        copy.replace_source_ring_info(state);
        assert_eq!(q, before);
        assert_eq!(q.source_ring_info().atom_members.len(), 3);
        assert_eq!(copy.source_ring_info().atom_members.len(), 4);
        assert_eq!(copy.source_ring_info().relevant_cycle_count, Some(7));
    }
}
#[cfg(test)]
mod all_source_conformers_tests {
    use super::*;
    use crate::{CoordinateSourceConformer as Row, CoordinateValidationError as Error};
    fn graph() -> QueryGraph {
        QueryGraph::from_parts(
            vec![QueryAtom::new(
                AtomId::new(0),
                crate::AtomSpec::new(crate::Element::C),
            )],
            vec![],
            [],
            vec![],
            vec![],
            vec![],
        )
        .unwrap()
    }
    #[test]
    fn empty_and_single_dimension_keep_duplicate_ids_and_append_order() {
        let mut q = graph();
        assert!(q.source_conformers().unwrap().is_empty());
        for id in [9, 2, 9] {
            q.add_conformer_3d(Conformer3D::new(id, vec![[id as f64, 0., 0.]], false))
                .unwrap();
        }
        q.set_source_conformer_order(None).unwrap();
        let rows = q.source_conformers().unwrap();
        let ids: Vec<_> = rows
            .iter()
            .map(|r| match r {
                Row::ThreeD(c) => c.id(),
                _ => panic!(),
            })
            .collect();
        assert_eq!(ids, [9, 2, 9]);
        match rows[1] {
            Row::ThreeD(c) => assert!(std::ptr::eq(c, &q.conformers_3d()[1])),
            _ => panic!(),
        }
    }
    #[test]
    fn mixed_physical_order_preserves_flags_and_borrowed_rows() {
        let mut q = graph();
        q.add_conformer_3d(Conformer3D::new(7, vec![[1., 2., 3.]], false))
            .unwrap();
        q = q.with_2d_coordinate_block(vec![[4., 5.]]).unwrap();
        q.add_conformer_3d(Conformer3D::new(7, vec![[6., 7., 8.]], true))
            .unwrap();
        let rows = q.source_conformers().unwrap();
        assert_eq!(rows.len(), 3);
        match rows[0] {
            Row::ThreeD(c) => {
                assert!(!c.is_3d());
                assert!(std::ptr::eq(c, &q.conformers_3d()[0]));
            }
            _ => panic!(),
        }
        match rows[1] {
            Row::TwoD(c) => assert!(std::ptr::eq(c, &q.conformers_2d()[0])),
            _ => panic!(),
        }
        match rows[2] {
            Row::ThreeD(c) => assert!(c.is_3d()),
            _ => panic!(),
        }
    }
    #[test]
    fn absent_mixed_order_propagates_without_dimensional_preference() {
        let mut q = graph().with_2d_coordinate_block(vec![[1., 2.]]).unwrap();
        q.add_conformer_3d(Conformer3D::new(1, vec![[1., 2., 3.]], true))
            .unwrap();
        q.set_source_conformer_order(None).unwrap();
        assert!(matches!(
            q.source_conformers(),
            Err(Error::MissingSourceConformerOrder)
        ));
        assert_eq!(q.conformers_2d().len(), 1);
        assert_eq!(q.conformers_3d().len(), 1);
    }
    #[test]
    fn raw_nonfinite_values_are_borrowed_without_validation_or_rewrite() {
        let bits = 0x7ff8_0000_0000_0042;
        let mut q = graph();
        q.add_conformer_3d(Conformer3D::new(
            3,
            vec![[f64::from_bits(bits), f64::INFINITY, f64::NEG_INFINITY]],
            false,
        ))
        .unwrap();
        match q.source_conformers().unwrap()[0] {
            Row::ThreeD(c) => {
                assert_eq!(c.coordinates()[0][0].to_bits(), bits);
                assert_eq!(c.coordinates()[0][1], f64::INFINITY);
                assert_eq!(c.coordinates()[0][2], f64::NEG_INFINITY);
            }
            _ => panic!(),
        }
    }
}
