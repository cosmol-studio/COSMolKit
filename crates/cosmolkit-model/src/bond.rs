// RDKit marker convention defined in dev/source_reproduction_protocol.md.

use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

use crate::AtomId;
use crate::PropertyValue;
use crate::property_value::PropertyStore;

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
    properties: PropertyStore,
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
            properties: PropertyStore::new(),
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
        value: impl Into<PropertyValue>,
    ) -> Result<Self, BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        self.properties.set(key, value.into());
        Ok(self)
    }

    #[must_use]
    pub fn with_computed_prop(
        mut self,
        key: impl Into<String>,
        value: impl Into<PropertyValue>,
    ) -> Result<Self, BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        self.properties.set_computed(key, value.into());
        Ok(self)
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<String, PropertyValue> {
        self.properties.values()
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&PropertyValue> {
        self.properties.get(key)
    }

    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.properties.is_computed(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        self.properties.computed_names()
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
    properties: PropertyStore,
    // BEGIN RDKIT CPP FUNCTION Bond::Bond(const Bond &) temporary flags
    // RDKit❗✔️: d_flags = other.d_flags;
    // END RDKIT CPP FUNCTION Bond::Bond(const Bond &) temporary flags
    // Bond's existing derived Clone copies this typed word. It is separate
    // from ordinary/computed properties, which may be cleared independently.
    temporary_flags: u64,
}

impl Bond {
    pub fn from_spec(id: BondId, spec: BondSpec) -> Self {
        // BEGIN RDKIT CPP MEMBER Bond::d_flags default
        // RDKit✔️✔️: std::uint64_t d_flags = 0;
        // END RDKIT CPP MEMBER Bond::d_flags default
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
            properties: spec.properties,
            temporary_flags: 0,
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

    /// Returns the source-compatible temporary bond flag word.
    #[doc(hidden)]
    #[must_use]
    pub const fn temporary_flags(&self) -> u64 {
        // BEGIN RDKIT CPP FUNCTION Bond::getFlags
        // RDKit✔️✔️: std::uint64_t getFlags() const { return d_flags; }
        // END RDKIT CPP FUNCTION Bond::getFlags
        self.temporary_flags
    }

    /// Sets the source-compatible temporary bond flag word.
    #[doc(hidden)]
    pub fn set_temporary_flags(&mut self, flags: u64) {
        // BEGIN RDKIT CPP FUNCTION Bond::setFlags
        // RDKit✔️✔️: void setFlags(std::uint64_t flags) { d_flags = flags; }
        // END RDKIT CPP FUNCTION Bond::setFlags
        self.temporary_flags = flags;
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
    pub fn props(&self) -> &BTreeMap<String, PropertyValue> {
        self.properties.values()
    }

    #[must_use]
    pub fn prop(&self, key: &str) -> Option<&PropertyValue> {
        self.properties.get(key)
    }

    pub const fn order_code(&self) -> i64 {
        self.order().rdkit_code()
    }
    pub const fn order_name(&self) -> &'static str {
        self.order().rdkit_name()
    }
    pub const fn direction_code(&self) -> i64 {
        self.direction().rdkit_code()
    }
    pub const fn direction_name(&self) -> &'static str {
        self.direction().rdkit_name()
    }
    pub const fn stereo_code(&self) -> i64 {
        self.stereo().rdkit_code()
    }
    pub const fn stereo_name(&self) -> &'static str {
        self.stereo().rdkit_name()
    }

    /// Return the persisted modern CIP neighbor ranking without suppressing
    /// malformed stored values or changing any rank/order.
    pub fn cip_neighbor_order(&self) -> Result<Option<Vec<u32>>, crate::CipDescriptorError> {
        crate::cip::neighbor_order_from_property(self.prop("_CIPNeighborOrder"))
    }

    /// Returns whether a property is registered as computed state.
    #[must_use]
    pub fn is_prop_computed(&self, key: &str) -> bool {
        self.properties.is_computed(key)
    }

    #[must_use]
    pub fn computed_prop_names(&self) -> &BTreeSet<String> {
        self.properties.computed_names()
    }

    /// Returns the modern CIP descriptor persisted on this bond, if present.
    pub fn cip_descriptor(
        &self,
    ) -> Result<Option<crate::CipDescriptor>, crate::CipDescriptorError> {
        // Reading a present property using the wrong kind must propagate its
        // typed conversion error; it is not an absent descriptor.
        let value = self
            .prop("_CIPCode")
            .map(|value| value.as_string())
            .transpose()?;
        crate::cip::descriptor_from_property(value)
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
        value: impl Into<PropertyValue>,
    ) -> Result<(), BondValueError> {
        let key = key.into();
        validate_property_key(&key)?;
        // RDKit✔️✔️: d_props.setVal(key, val);
        // A non-computed write does not remove an existing computed marker.
        self.properties.set(key, value.into());
        Ok(())
    }

    #[doc(hidden)]
    pub fn set_computed_prop(
        &mut self,
        key: impl Into<String>,
        value: impl Into<PropertyValue>,
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
        self.properties.set_computed(key, value.into());
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
        self.properties.clear(key);
    }

    #[doc(hidden)]
    pub fn clear_computed_props(&mut self) {
        // RDKit✔️🔝: for (const auto &key : compLst) {
        // RDKit✔️🔝:   d_props.clearVal(key);
        // RDKit✔️🔝: }
        // Moving the set avoids the source vector copy while preserving exact
        // membership-based clearing.
        self.properties.clear_computed();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{BondQueryPredicate, QueryBond, QueryNode};

    fn property_order(bond: &Bond) -> Vec<&str> {
        bond.properties
            .ordered_keys()
            .iter()
            .map(String::as_str)
            .collect()
    }

    #[test]
    fn typed_property_transport_bond_copy_and_query_carrier_preserve_values_and_order() {
        let spec = BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single)
            .with_prop("first", PropertyValue::String("seven".to_owned()))
            .unwrap()
            .with_computed_prop("computed", PropertyValue::Bool(false))
            .unwrap()
            .with_prop("first", PropertyValue::Int(7))
            .unwrap();
        assert_eq!(spec.properties.ordered_keys(), &["first", "computed"]);
        assert_eq!(spec.prop("first"), Some(&PropertyValue::Int(7)));
        assert!(spec.is_prop_computed("computed"));

        let bond = Bond::from_spec(BondId::new(0), spec);
        let source = bond.clone();
        let mut query = QueryBond::from_parts(
            bond.clone(),
            QueryNode::predicate(BondQueryPredicate::Order(BondOrder::Single)),
        );
        assert_eq!(query.bond(), &source);
        assert_eq!(property_order(query.bond()), vec!["first", "computed"]);
        assert_eq!(query.bond().prop("first"), Some(&PropertyValue::Int(7)));
        assert_eq!(
            query.bond().prop("computed"),
            Some(&PropertyValue::Bool(false))
        );
        assert!(query.bond().is_prop_computed("computed"));

        query.bond_mut().clear_prop("first");
        query
            .bond_mut()
            .set_prop("first", PropertyValue::Double(1.25))
            .unwrap();
        assert_eq!(property_order(query.bond()), vec!["computed", "first"]);
        assert_eq!(
            query.bond().prop("first"),
            Some(&PropertyValue::Double(1.25))
        );
        assert_eq!(bond, source);
        assert_eq!(property_order(&bond), vec!["first", "computed"]);
    }
}

#[cfg(test)]
mod flags_tests {
    use super::*;

    fn single_bond(id: usize) -> Bond {
        Bond::from_spec(
            BondId::new(id),
            BondSpec::new(AtomId::new(0), AtomId::new(1), BondOrder::Single),
        )
    }

    #[test]
    fn cf3d_flags_bond_zero_default_every_bit_and_fixed_masks() {
        assert_eq!(single_bond(0).temporary_flags(), 0);

        for bit in 0..u64::BITS {
            let mut bond = single_bond(bit as usize);
            let expected = 1_u64 << bit;
            bond.set_temporary_flags(expected);
            assert_eq!(bond.temporary_flags(), expected, "bit {bit}");
        }

        for expected in [u64::MAX, 0xAAAA_AAAA_AAAA_AAAA, 0x5555_5555_5555_5555] {
            let mut bond = single_bond(0);
            bond.set_temporary_flags(expected);
            assert_eq!(bond.temporary_flags(), expected);
        }
    }

    #[test]
    fn cf3d_flags_bond_clone_remap_clear_and_equality_preserve_word() {
        let mut source = single_bond(0);
        source.set_temporary_flags(u64::MAX);
        source
            .set_prop("ordinary", "kept")
            .expect("valid ordinary property");
        source
            .set_computed_prop("derived", "cleared")
            .expect("valid computed property");

        let cloned = source.clone();
        assert_eq!(cloned, source);
        assert_eq!(cloned.temporary_flags(), u64::MAX);

        let remapped =
            source
                .clone()
                .remapped(BondId::new(7), AtomId::new(3), AtomId::new(4), None);
        assert_eq!(remapped.id(), BondId::new(7));
        assert_eq!(remapped.begin(), AtomId::new(3));
        assert_eq!(remapped.end(), AtomId::new(4));
        assert_eq!(remapped.temporary_flags(), u64::MAX);

        let mut cleared = source.clone();
        cleared.clear_computed_props();
        assert_eq!(cleared.temporary_flags(), u64::MAX);
        assert_eq!(cleared.prop("derived"), None);
        assert_eq!(
            cleared.prop("ordinary"),
            Some(&crate::PropertyValue::String("kept".to_owned()))
        );

        let mut changed_clone = source.clone();
        changed_clone.set_temporary_flags(1);
        assert_eq!(source.temporary_flags(), u64::MAX, "clone is isolated");
        assert_ne!(
            changed_clone, source,
            "equality remains representation-based"
        );
    }
}
