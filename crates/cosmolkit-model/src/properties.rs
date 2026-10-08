//! Molecule-level property values shared by parsers and algorithms.

use crate::PropertyText;
use std::collections::BTreeMap;

use crate::{AtomId, BondId, PropertyValue};

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
pub enum MoleculePropertyError {
    #[error("property key must not be empty")]
    EmptyKey,
    #[error("computed property list has the wrong value kind: {0}")]
    ComputedListKind(crate::PropertyValueError),
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MoleculeProperties {
    name: Option<PropertyText>,
    sdf_data_fields: Vec<(PropertyText, PropertyText)>,
    sdf_property_lists: Vec<SdfPropertyList>,
    props: crate::property_value::PropertyStore,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfPropertyListTarget {
    Atom,
    Bond,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfPropertyList {
    target: SdfPropertyListTarget,
    name: PropertyText,
    values: Vec<Option<PropertyValue>>,
}

impl SdfPropertyList {
    #[must_use]
    pub fn new(
        target: SdfPropertyListTarget,
        name: impl Into<PropertyText>,
        values: Vec<Option<PropertyValue>>,
    ) -> Self {
        Self {
            target,
            name: name.into(),
            values,
        }
    }

    #[must_use]
    pub const fn target(&self) -> SdfPropertyListTarget {
        self.target
    }

    #[must_use]
    pub fn name(&self) -> &PropertyText {
        &self.name
    }

    #[must_use]
    pub fn values(&self) -> &[Option<PropertyValue>] {
        &self.values
    }
}

impl MoleculeProperties {
    pub(crate) fn source_dict(&self) -> &crate::property_value::PropertyStore {
        &self.props
    }
    pub(crate) fn from_source_dict(props: crate::property_value::PropertyStore) -> Self {
        Self {
            props,
            ..Self::default()
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
        MoleculePropertyError,
    > {
        self.props
            .filtered_ordered(include_private, include_computed)
            .map_err(MoleculePropertyError::from)
    }

    #[must_use]
    pub fn name(&self) -> Option<&PropertyText> {
        self.name.as_ref()
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<PropertyText>) -> Self {
        self.name = Some(name.into());
        self
    }

    #[must_use]
    pub fn sdf_data_fields(&self) -> &[(PropertyText, PropertyText)] {
        &self.sdf_data_fields
    }

    #[must_use]
    pub fn sdf_property_lists(&self) -> &[SdfPropertyList] {
        &self.sdf_property_lists
    }

    #[must_use]
    pub fn props(&self) -> &BTreeMap<PropertyText, PropertyValue> {
        self.props.values()
    }

    #[must_use]
    /// Borrow the canonical property records in source insertion order.
    #[doc(hidden)]
    pub fn ordered_props(
        &self,
    ) -> impl ExactSizeIterator<Item = (&PropertyText, &PropertyValue)> + '_ {
        self.props.ordered()
    }

    pub fn prop(&self, key: impl AsRef<[u8]>) -> Option<&PropertyValue> {
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

    /// Returns whether a property is registered as computed state.
    #[must_use]
    pub fn is_prop_computed(
        &self,
        key: impl AsRef<[u8]>,
    ) -> Result<bool, crate::PropertyValueError> {
        self.props.is_computed(key)
    }

    #[must_use]
    pub fn computed_prop_names(
        &self,
    ) -> Result<Option<&[PropertyText]>, crate::PropertyValueError> {
        self.props.computed_names()
    }

    pub fn with_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<Self, MoleculePropertyError> {
        self.set_prop(key, value)?;
        Ok(self)
    }

    pub fn with_computed_prop(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<Self, MoleculePropertyError> {
        self.set_computed_prop(key, value)?;
        Ok(self)
    }

    #[must_use]
    pub fn with_sdf_data_field(
        mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyText>,
    ) -> Self {
        self.sdf_data_fields.push((key.into(), value.into()));
        self
    }

    #[must_use]
    pub fn with_sdf_property_list(mut self, property_list: SdfPropertyList) -> Self {
        self.sdf_property_lists.push(property_list);
        self
    }

    pub fn set_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<(), MoleculePropertyError> {
        self.props
            .set(key.into(), value.into())
            .map_err(MoleculePropertyError::from)
    }

    pub fn set_computed_prop(
        &mut self,
        key: impl Into<PropertyText>,
        value: impl Into<PropertyValue>,
    ) -> Result<(), MoleculePropertyError> {
        self.props
            .set_computed(key.into(), value.into())
            .map_err(MoleculePropertyError::from)
    }

    /// Retain native computed-list bookkeeping for a typed transient value
    /// owned by an algorithm and consumed before its detached result returns.
    /// The caller must clear this name through clear_prop when consuming it.
    #[doc(hidden)]
    pub fn register_transient_computed_name(
        &mut self,
        key: impl Into<PropertyText>,
    ) -> Result<(), MoleculePropertyError> {
        // RDKit✔️❌:     if (computed) {
        // RDKit✔️❌:       STR_VECT compLst;
        // RDKit✔️❌:       getPropIfPresent(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       if (std::find(compLst.begin(), compLst.end(), key) == compLst.end()) {
        // RDKit✔️❌:         compLst.emplace_back(key);
        // RDKit✔️❌:         d_props.setVal(RDKit::detail::computedPropName, compLst);
        // RDKit✔️❌:       }
        // RDKit✔️❌:     }
        self.props
            .register_transient_computed_name(&key.into())
            .map_err(MoleculePropertyError::from)
    }

    pub fn clear_prop(&mut self, key: impl AsRef<[u8]>) -> Result<(), MoleculePropertyError> {
        self.props.clear(key).map_err(MoleculePropertyError::from)
    }

    pub fn clear_computed_props(&mut self) -> Result<(), MoleculePropertyError> {
        self.props
            .clear_computed()
            .map_err(MoleculePropertyError::from)
    }

    pub fn remap_topology(
        &mut self,
        atom_new_to_old: &[Option<AtomId>],
        bond_new_to_old: &[Option<BondId>],
    ) {
        self.sdf_property_lists = self
            .sdf_property_lists
            .iter()
            .map(|property_list| property_list.remapped_topology(atom_new_to_old, bond_new_to_old))
            .collect();
    }
}

impl SdfPropertyList {
    fn remapped_topology(
        &self,
        atom_new_to_old: &[Option<AtomId>],
        bond_new_to_old: &[Option<BondId>],
    ) -> Self {
        let values = match self.target {
            SdfPropertyListTarget::Atom => atom_new_to_old
                .iter()
                .map(|old_row| {
                    old_row.and_then(|row| self.values.get(row.index()).cloned().flatten())
                })
                .collect(),
            SdfPropertyListTarget::Bond => bond_new_to_old
                .iter()
                .map(|old_row| {
                    old_row.and_then(|row| self.values.get(row.index()).cloned().flatten())
                })
                .collect(),
        };
        Self {
            target: self.target,
            name: self.name.clone(),
            values,
        }
    }
}

impl From<crate::property_value::PropertyStoreError> for MoleculePropertyError {
    fn from(error: crate::property_value::PropertyStoreError) -> Self {
        match error {
            crate::property_value::PropertyStoreError::EmptyKey => Self::EmptyKey,
            crate::property_value::PropertyStoreError::ComputedListKind(source) => {
                Self::ComputedListKind(source)
            }
        }
    }
}
