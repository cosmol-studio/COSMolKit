// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

use crate::AtomId;

pub use cosmolkit_types::{BondDirection, BondOrder, BondStereo};

/// Stable bond-table index.
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
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{}", self.0)
    }
}

/// A detached bond value violates a bond-local source constraint.
#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum BondValueError {
    #[error("CIS/TRANS bond stereo requires two reference atoms")]
    StereoAtomsRequired,
    #[error("bond property key cannot be empty")]
    EmptyPropertyKey,
}

fn validate_property_key(key: &str) -> Result<(), BondValueError> {
    // BEGIN RDKIT CPP FUNCTION RDProps::setProp empty-key precondition
    // RDKit✔️✔️: if(key.empty()) {
    // RDKit✔️✔️:   throw ValueErrorException("Cannot set property with empty key");
    // RDKit✔️✔️: }
    // END RDKIT CPP FUNCTION RDProps::setProp empty-key precondition
    if key.is_empty() {
        Err(BondValueError::EmptyPropertyKey)
    } else {
        Ok(())
    }
}

fn validate_stereo_references(
    stereo: BondStereo,
    stereo_atoms: Option<[AtomId; 2]>,
) -> Result<(), BondValueError> {
    if matches!(stereo, BondStereo::Cis | BondStereo::Trans) && stereo_atoms.is_none() {
        Err(BondValueError::StereoAtomsRequired)
    } else {
        Ok(())
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
    pub fn with_prop(
        mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<Self, BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        self.props.insert(key, value.into());
        Ok(self)
    }

    #[must_use]
    pub fn with_computed_prop(
        mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<Self, BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
        Ok(self)
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

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        &self.computed_props
    }

    pub fn validate(&self) -> Result<(), BondValueError> {
        validate_stereo_references(self.stereo, self.stereo_atoms)
    }

    pub fn remapped_endpoints(
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
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
}

impl Bond {
    pub fn from_spec(id: BondId, spec: BondSpec) -> Self {
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
            props: spec.props,
            computed_props: spec.computed_props,
        }
    }

    pub fn validate(&self) -> Result<(), BondValueError> {
        validate_stereo_references(self.stereo, self.stereo_atoms)
    }

    #[doc(hidden)]
    pub fn remapped(
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
        // RDKit✔️✔️: unsigned int getIdx() const { return d_index; }
        self.id
    }

    #[doc(hidden)]
    pub fn set_id_for_construction(&mut self, id: BondId) {
        // RDKit✔️✔️: void setIdx(unsigned int index) { d_index = index; }
        self.id = id;
    }

    #[must_use]
    pub const fn begin(&self) -> AtomId {
        // RDKit✔️✔️: unsigned int getBeginAtomIdx() const { return d_beginAtomIdx; }
        self.begin
    }

    #[must_use]
    pub const fn end(&self) -> AtomId {
        // RDKit✔️✔️: unsigned int getEndAtomIdx() const { return d_endAtomIdx; }
        self.end
    }

    #[must_use]
    pub const fn order(&self) -> BondOrder {
        // RDKit✔️✔️: BondType getBondType() const { return static_cast<BondType>(d_bondType); }
        self.order
    }

    #[must_use]
    pub const fn is_aromatic(&self) -> bool {
        // RDKit✔️✔️: bool getIsAromatic() const { return df_isAromatic; }
        self.is_aromatic
    }

    #[must_use]
    pub const fn is_conjugated(&self) -> bool {
        // RDKit✔️✔️: bool getIsConjugated() const { return df_isConjugated; }
        self.is_conjugated
    }

    #[must_use]
    pub const fn direction(&self) -> BondDirection {
        // RDKit✔️✔️: BondDir getBondDir() const { return static_cast<BondDir>(d_dirTag); }
        self.direction
    }

    #[must_use]
    pub const fn stereo(&self) -> BondStereo {
        // RDKit✔️✔️: BondStereo getStereo() const { return static_cast<BondStereo>(d_stereo); }
        self.stereo
    }

    #[must_use]
    pub const fn stereo_atoms(&self) -> Option<[AtomId; 2]> {
        // RDKit✔️✔️: const INT_VECT &getStereoAtoms() const {
        // RDKit✔️✔️:   if (!dp_stereoAtoms) {
        // RDKit✔️✔️:     const_cast<Bond *>(this)->dp_stereoAtoms = new INT_VECT();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return *dp_stereoAtoms;
        // RDKit✔️✔️: }
        // The empty source vector is projected as `None`.
        self.stereo_atoms
    }

    #[must_use]
    pub const fn unknown_stereo(&self) -> bool {
        self.unknown_stereo
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
        crate::cip::descriptor_from_property(self.prop("_CIPCode"))
    }

    #[doc(hidden)]
    pub fn set_order(&mut self, order: BondOrder) {
        // RDKit✔️✔️: void setBondType(BondType bT) { d_bondType = bT; }
        self.order = order;
    }

    #[doc(hidden)]
    pub fn set_endpoints(&mut self, begin: AtomId, end: AtomId) {
        self.begin = begin;
        self.end = end;
    }

    #[doc(hidden)]
    pub fn set_aromatic(&mut self, is_aromatic: bool) {
        // RDKit✔️✔️: void setIsAromatic(bool what) { df_isAromatic = what; }
        self.is_aromatic = is_aromatic;
    }

    #[doc(hidden)]
    pub fn set_conjugated(&mut self, is_conjugated: bool) {
        // RDKit✔️✔️: void setIsConjugated(bool what) { df_isConjugated = what; }
        self.is_conjugated = is_conjugated;
    }

    #[doc(hidden)]
    pub fn set_direction(&mut self, direction: BondDirection) {
        // RDKit✔️✔️: void setBondDir(BondDir what) { d_dirTag = what; }
        self.direction = direction;
    }

    #[doc(hidden)]
    pub fn set_stereo(&mut self, stereo: BondStereo) -> Result<(), BondValueError> {
        // BEGIN RDKIT CPP FUNCTION Bond::setStereo
        // RDKit✔️✔️: void setStereo(BondStereo what) {
        // RDKit✔️✔️:   PRECONDITION(((what != STEREOCIS && what != STEREOTRANS) ||
        // RDKit✔️✔️:                 getStereoAtoms().size() == 2),
        // RDKit✔️✔️:                "Stereo atoms should be specified before specifying CIS/TRANS "
        // RDKit✔️✔️:                "bond stereochemistry")
        // RDKit✔️✔️:   d_stereo = what;
        // RDKit✔️✔️: }
        // END RDKIT CPP FUNCTION Bond::setStereo
        validate_stereo_references(stereo, self.stereo_atoms)?;
        self.stereo = stereo;
        Ok(())
    }

    #[doc(hidden)]
    pub fn set_stereo_atoms(&mut self, stereo_atoms: Option<[AtomId; 2]>) {
        self.stereo_atoms = stereo_atoms;
    }

    #[doc(hidden)]
    pub fn set_unknown_stereo(&mut self, unknown_stereo: bool) {
        self.unknown_stereo = unknown_stereo;
    }

    #[doc(hidden)]
    pub fn set_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        // A non-computed write does not remove an existing computed marker.
        self.props.insert(key, value.into());
        Ok(())
    }

    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<String>,
    ) -> Result<(), BondValueError> {
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
        validate_property_key(&key)?;
        self.props.insert(key.clone(), value.into());
        self.computed_props.insert(key);
        Ok(())
    }

    #[doc(hidden)]
    pub fn clear_prop(&mut self, key: &str) {
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

    #[doc(hidden)]
    pub fn clear_computed_props(&mut self) {
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
