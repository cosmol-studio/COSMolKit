// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

use crate::{AtomId, BondQueryPredicate, QueryNode};

/// Stable bond table index.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BondId(usize);

impl BondId {
    #[must_use]
    pub const fn new(index: usize) -> Self {
        Self(index)
    }

    #[must_use]
    pub const fn index(self) -> usize {
        self.0
    }
}

impl fmt::Display for BondId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// COSMolKit core bond order.
///
/// RDKit aliases and integer codes belong in a compatibility module, not here.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondOrder {
    Null,
    Single,
    Double,
    Triple,
    Quadruple,
    Quintuple,
    Hextuple,
    OneAndHalf,
    TwoAndHalf,
    ThreeAndHalf,
    FourAndHalf,
    FiveAndHalf,
    Aromatic,
    Ionic,
    Dative,
    DativeOne,
    DativeLeft,
    DativeRight,
    Hydrogen,
    ThreeCenter,
    Other,
    Zero,
    Unspecified,
}

impl BondOrder {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::Null | Self::Unspecified => 0,
            Self::Single => 1,
            Self::Double => 2,
            Self::Triple => 3,
            Self::Quadruple => 4,
            Self::Quintuple => 5,
            Self::Hextuple => 6,
            Self::OneAndHalf => 7,
            Self::TwoAndHalf => 8,
            Self::ThreeAndHalf => 9,
            Self::FourAndHalf => 10,
            Self::FiveAndHalf => 11,
            Self::Aromatic => 12,
            Self::Ionic => 13,
            Self::Dative => 14,
            Self::DativeOne => 15,
            Self::DativeLeft => 16,
            Self::DativeRight => 17,
            Self::Hydrogen => 18,
            Self::ThreeCenter => 19,
            Self::Other => 20,
            Self::Zero => 21,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::Null | Self::Unspecified => "UNSPECIFIED",
            Self::Single => "SINGLE",
            Self::Double => "DOUBLE",
            Self::Triple => "TRIPLE",
            Self::Quadruple => "QUADRUPLE",
            Self::Quintuple => "QUINTUPLE",
            Self::Hextuple => "HEXTUPLE",
            Self::OneAndHalf => "ONEANDAHALF",
            Self::TwoAndHalf => "TWOANDAHALF",
            Self::ThreeAndHalf => "THREEANDAHALF",
            Self::FourAndHalf => "FOURANDAHALF",
            Self::FiveAndHalf => "FIVEANDAHALF",
            Self::Aromatic => "AROMATIC",
            Self::Ionic => "IONIC",
            Self::Dative => "DATIVE",
            Self::DativeOne => "DATIVEONE",
            Self::DativeLeft => "DATIVEL",
            Self::DativeRight => "DATIVER",
            Self::Hydrogen => "HYDROGEN",
            Self::ThreeCenter => "THREECENTER",
            Self::Other => "OTHER",
            Self::Zero => "ZERO",
        }
    }

    #[must_use]
    pub fn from_rdkit_name(name: &str) -> Option<Self> {
        match name {
            "UNSPECIFIED" | "ZERO" => Some(Self::Unspecified),
            "SINGLE" => Some(Self::Single),
            "DOUBLE" => Some(Self::Double),
            "TRIPLE" => Some(Self::Triple),
            "QUADRUPLE" => Some(Self::Quadruple),
            "QUINTUPLE" => Some(Self::Quintuple),
            "HEXTUPLE" => Some(Self::Hextuple),
            "ONEANDAHALF" => Some(Self::OneAndHalf),
            "TWOANDAHALF" => Some(Self::TwoAndHalf),
            "THREEANDAHALF" => Some(Self::ThreeAndHalf),
            "FOURANDAHALF" => Some(Self::FourAndHalf),
            "FIVEANDAHALF" => Some(Self::FiveAndHalf),
            "AROMATIC" => Some(Self::Aromatic),
            "IONIC" => Some(Self::Ionic),
            "DATIVE" => Some(Self::Dative),
            "DATIVEONE" => Some(Self::DativeOne),
            "DATIVEL" => Some(Self::DativeLeft),
            "DATIVER" => Some(Self::DativeRight),
            "HYDROGEN" => Some(Self::Hydrogen),
            "THREECENTER" => Some(Self::ThreeCenter),
            "OTHER" => Some(Self::Other),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondDirection {
    None,
    BeginWedge,
    BeginDash,
    EndUpRight,
    EndDownRight,
    EitherDouble,
    Unknown,
}

impl BondDirection {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::None => 0,
            Self::BeginWedge => 1,
            Self::BeginDash => 2,
            Self::EndUpRight => 3,
            Self::EndDownRight => 4,
            Self::EitherDouble => 5,
            Self::Unknown => 6,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::BeginWedge => "BEGINWEDGE",
            Self::BeginDash => "BEGINDASH",
            Self::EndUpRight => "ENDUPRIGHT",
            Self::EndDownRight => "ENDDOWNRIGHT",
            Self::EitherDouble => "EITHERDOUBLE",
            Self::Unknown => "UNKNOWN",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BondStereo {
    None,
    Any,
    Z,
    E,
    Cis,
    Trans,
    AtropCw,
    AtropCcw,
}

impl BondStereo {
    #[must_use]
    pub const fn python_code(self) -> i64 {
        match self {
            Self::None => 0,
            Self::Any => 1,
            Self::Z => 2,
            Self::E => 3,
            Self::Cis => 4,
            Self::Trans => 5,
            Self::AtropCw => 6,
            Self::AtropCcw => 7,
        }
    }

    #[must_use]
    pub const fn rdkit_name(self) -> &'static str {
        match self {
            Self::None => "NONE",
            Self::Any => "ANY",
            Self::Z => "Z",
            Self::E => "E",
            Self::Cis => "CIS",
            Self::Trans => "TRANS",
            Self::AtropCw => "ATROP_CW",
            Self::AtropCcw => "ATROP_CCW",
        }
    }
}

/// Bond construction payload. Builders assign `BondId`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BondSpec {
    begin: AtomId,
    end: AtomId,
    order: BondOrder,
    is_aromatic: bool,
    is_conjugated: bool,
    direction: BondDirection,
    stereo: BondStereo,
    stereo_atoms: Option<[AtomId; 2]>,
    unknown_stereo: bool,
    query: Option<QueryNode<BondQueryPredicate>>,
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
}

impl BondSpec {
    #[must_use]
    pub const fn new(begin: AtomId, end: AtomId, order: BondOrder) -> Self {
        Self {
            begin,
            end,
            order,
            is_aromatic: false,
            is_conjugated: false,
            direction: BondDirection::None,
            stereo: BondStereo::None,
            stereo_atoms: None,
            unknown_stereo: false,
            query: None,
            props: BTreeMap::new(),
            computed_props: BTreeSet::new(),
        }
    }

    #[must_use]
    pub const fn begin(&self) -> AtomId {
        self.begin
    }

    #[must_use]
    pub const fn end(&self) -> AtomId {
        self.end
    }

    #[must_use]
    pub const fn order(&self) -> BondOrder {
        self.order
    }

    #[must_use]
    pub const fn with_order(mut self, order: BondOrder) -> Self {
        self.order = order;
        self
    }

    #[must_use]
    pub const fn direction(&self) -> BondDirection {
        self.direction
    }

    #[must_use]
    pub const fn stereo(&self) -> BondStereo {
        self.stereo
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    #[must_use]
    pub const fn is_conjugated(&self) -> bool {
        self.is_conjugated
    }

    #[must_use]
    pub const fn with_aromatic(mut self, is_aromatic: bool) -> Self {
        self.is_aromatic = is_aromatic;
        self
    }

    #[must_use]
    pub const fn with_conjugated(mut self, is_conjugated: bool) -> Self {
        self.is_conjugated = is_conjugated;
        self
    }

    #[must_use]
    pub const fn with_direction(mut self, direction: BondDirection) -> Self {
        self.direction = direction;
        self
    }

    #[must_use]
    pub const fn with_stereo(mut self, stereo: BondStereo) -> Self {
        self.stereo = stereo;
        self
    }

    #[must_use]
    pub const fn with_stereo_atoms(mut self, begin_ref: AtomId, end_ref: AtomId) -> Self {
        self.stereo_atoms = Some([begin_ref, end_ref]);
        self
    }

    #[must_use]
    pub const fn without_stereo_atoms(mut self) -> Self {
        self.stereo_atoms = None;
        self
    }

    #[must_use]
    pub const fn stereo_atoms(&self) -> Option<[AtomId; 2]> {
        self.stereo_atoms
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn with_unknown_stereo(mut self, unknown_stereo: bool) -> Self {
        self.unknown_stereo = unknown_stereo;
        self
    }

    #[must_use]
    pub fn with_query(mut self, query: QueryNode<BondQueryPredicate>) -> Self {
        self.query = Some(query);
        self
    }

    #[must_use]
    pub fn without_query(mut self) -> Self {
        self.query = None;
        self
    }

    #[must_use]
    pub const fn query(&self) -> Option<&QueryNode<BondQueryPredicate>> {
        self.query.as_ref()
    }

    #[must_use]
    pub fn with_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.props.insert(key.into(), value.into());
        self
    }

    #[must_use]
    pub fn with_computed_prop(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        let key = key.into();
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
        self
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.computed_props.contains(key)
    }

    pub(crate) fn remapped_endpoints(
        &self,
        begin: AtomId,
        end: AtomId,
        stereo_atoms: Option<[AtomId; 2]>,
    ) -> Self {
        let mut spec = self.clone();
        spec.begin = begin;
        spec.end = end;
        spec.stereo_atoms = stereo_atoms;
        spec
    }
}

/// Immutable bond record owned by `Molecule`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Bond {
    id: BondId,
    begin: AtomId,
    end: AtomId,
    order: BondOrder,
    is_aromatic: bool,
    is_conjugated: bool,
    direction: BondDirection,
    stereo: BondStereo,
    stereo_atoms: Option<[AtomId; 2]>,
    unknown_stereo: bool,
    query: Option<QueryNode<BondQueryPredicate>>,
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
}

impl Bond {
    pub(crate) fn from_spec(id: BondId, spec: BondSpec) -> Self {
        Self {
            id,
            begin: spec.begin,
            end: spec.end,
            order: spec.order,
            is_aromatic: spec.is_aromatic,
            is_conjugated: spec.is_conjugated,
            direction: spec.direction,
            stereo: spec.stereo,
            stereo_atoms: spec.stereo_atoms,
            unknown_stereo: spec.unknown_stereo,
            query: spec.query,
            props: spec.props,
            computed_props: spec.computed_props,
        }
    }

    #[allow(dead_code)]
    pub(crate) fn remapped(
        mut self,
        id: BondId,
        begin: AtomId,
        end: AtomId,
        stereo_atoms: Option<[AtomId; 2]>,
    ) -> Self {
        self.id = id;
        self.begin = begin;
        self.end = end;
        self.stereo_atoms = stereo_atoms;
        self
    }

    #[must_use]
    pub const fn id(&self) -> BondId {
        self.id
    }

    pub(crate) fn set_id_for_construction(&mut self, id: BondId) {
        self.id = id;
    }

    #[must_use]
    pub const fn begin(&self) -> AtomId {
        self.begin
    }

    #[must_use]
    pub const fn end(&self) -> AtomId {
        self.end
    }

    #[must_use]
    pub const fn order(&self) -> BondOrder {
        self.order
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        self.is_aromatic
    }

    #[must_use]
    pub const fn is_conjugated(&self) -> bool {
        self.is_conjugated
    }

    #[must_use]
    pub const fn direction(&self) -> BondDirection {
        self.direction
    }

    #[must_use]
    pub const fn stereo(&self) -> BondStereo {
        self.stereo
    }

    #[must_use]
    pub const fn stereo_atoms(&self) -> Option<[AtomId; 2]> {
        self.stereo_atoms
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
    }

    #[must_use]
    pub const fn query(&self) -> Option<&QueryNode<BondQueryPredicate>> {
        self.query.as_ref()
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, String> {
        &self.props
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&str> {
        self.props.get(key).map(String::as_str)
    }

    /// Returns whether a property is registered as computed state.
    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.computed_props.contains(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        &self.computed_props
    }

    /// Returns the modern CIP descriptor persisted on this bond, if present.
    pub fn cip_descriptor(
        &self,
    ) -> Result<Option<crate::CipDescriptor>, crate::CipDescriptorError> {
        crate::chemistry::cip::descriptor_from_property(self.prop("_CIPCode"))
    }

    #[allow(dead_code)]
    pub(crate) fn set_order(&mut self, order: BondOrder) {
        self.order = order;
    }

    #[allow(dead_code)]
    pub(crate) fn set_endpoints(&mut self, begin: AtomId, end: AtomId) {
        self.begin = begin;
        self.end = end;
    }

    #[allow(dead_code)]
    pub(crate) fn set_aromatic(&mut self, is_aromatic: bool) {
        self.is_aromatic = is_aromatic;
    }

    #[allow(dead_code)]
    pub(crate) fn set_conjugated(&mut self, is_conjugated: bool) {
        self.is_conjugated = is_conjugated;
    }

    #[allow(dead_code)]
    pub(crate) fn set_direction(&mut self, direction: BondDirection) {
        self.direction = direction;
    }

    #[allow(dead_code)]
    pub(crate) fn set_stereo(&mut self, stereo: BondStereo) {
        // BEGIN RDKIT CPP FUNCTION Bond::setStereo
        // RDKit✔️✔️: void setStereo(BondStereo what) {
        // RDKit✔️✔️:   PRECONDITION(((what != STEREOCIS && what != STEREOTRANS) ||
        // RDKit✔️✔️:                 getStereoAtoms().size() == 2),
        // RDKit✔️✔️:                "Stereo atoms should be specified before specifying CIS/TRANS "
        // RDKit✔️✔️:                "bond stereochemistry")
        // RDKit✔️✔️:   d_stereo = what;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Bond::setStereo
        debug_assert!(
            !matches!(stereo, BondStereo::Cis | BondStereo::Trans) || self.stereo_atoms.is_some(),
            "stereo atoms should be specified before CIS/TRANS bond stereochemistry"
        );
        self.stereo = stereo;
    }

    #[allow(dead_code)]
    pub(crate) fn set_stereo_atoms(&mut self, stereo_atoms: Option<[AtomId; 2]>) {
        self.stereo_atoms = stereo_atoms;
    }

    #[allow(dead_code)]
    pub(crate) fn set_unknown_stereo(&mut self, unknown_stereo: bool) {
        self.unknown_stereo = unknown_stereo;
    }

    #[allow(dead_code)]
    pub(crate) fn set_query(&mut self, query: Option<QueryNode<BondQueryPredicate>>) {
        self.query = query;
    }

    #[allow(dead_code)]
    pub(crate) fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        // RDKit✔️✔️: d_props.setVal(key, val);
        // A non-computed write does not remove an existing computed marker.
        self.props.insert(key.into(), value.into());
    }

    #[allow(dead_code)]
    pub(crate) fn set_computed_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        // RDKit✔️🔝: if (computed) {
        // RDKit✔️🔝:   STR_VECT compLst;
        // RDKit✔️🔝:   getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝:   if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️🔝:     compLst.emplace_back(key);
        // RDKit✔️🔝:     d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝:   }
        // RDKit✔️🔝: }
        // RDKit✔️🔝: d_props.setVal(key, val);
        // The ordered set preserves membership semantics while replacing the
        // source vector's linear duplicate scan with logarithmic insertion.
        let key = key.into();
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
    }

    #[allow(dead_code)]
    pub(crate) fn clear_prop(&mut self, key: &str) {
        // RDKit✔️🔝: auto svi = std::find(compLst.begin(), compLst.end(), key);
        // RDKit✔️🔝: if (svi != compLst.end()) {
        // RDKit✔️🔝:   compLst.erase(svi);
        // RDKit✔️🔝:   d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️🔝: }
        // RDKit✔️🔝: d_props.clearVal(key);
        // BTreeSet removal preserves the source transition with logarithmic
        // lookup instead of the source vector's linear search and erase.
        self.props.remove(key);
        self.computed_props.remove(key);
    }

    pub(crate) fn clear_computed_props(&mut self) {
        // RDKit✔️🔝: for (const auto &key : compLst) {
        // RDKit✔️🔝:   d_props.clearVal(key);
        // RDKit✔️🔝: }
        // Moving the set avoids the source vector copy while preserving exact
        // membership-based clearing.
        for key in std::mem::take(&mut self.computed_props) {
            self.props.remove(&key);
        }
    }
}
