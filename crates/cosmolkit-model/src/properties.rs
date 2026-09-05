//! Molecule-level property values shared by parsers and algorithms.

use std::collections::{BTreeMap, BTreeSet};

use crate::{AtomId, BondId};

#[derive(Debug, Clone, PartialEq, Default)]
pub struct PropertyStore {
    pub name: Option<String>,
    pub sdf_data_fields: Vec<(String, String)>,
    pub sdf_property_lists: Vec<SdfPropertyList>,
    pub props: std::collections::BTreeMap<String, String>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct MoleculeProperties {
    name: Option<String>,
    sdf_data_fields: Vec<(String, String)>,
    sdf_property_lists: Vec<SdfPropertyList>,
    props: BTreeMap<String, String>,
    computed_props: BTreeSet<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SdfPropertyListTarget {
    Atom,
    Bond,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SdfPropertyList {
    target: SdfPropertyListTarget,
    name: String,
    values: Vec<Option<String>>,
}

impl SdfPropertyList {
    #[must_use]
    pub fn new(target: SdfPropertyListTarget, name: impl Into<String>, values: Vec<Option<String>>) -> Self {
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
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub fn values(&self) -> &[Option<String>] {
        &self.values
    }
}

impl MoleculeProperties {
    #[must_use]
    pub fn name(&self) -> Option<&str> {
        self.name.as_deref()
    }

    #[must_use]
    pub fn with_name(mut self, name: impl Into<String>) -> Self {
        self.name = Some(name.into());
        self
    }

    #[must_use]
    pub fn sdf_data_fields(&self) -> &[(String, String)] {
        &self.sdf_data_fields
    }

    #[must_use]
    pub fn sdf_property_lists(&self) -> &[SdfPropertyList] {
        &self.sdf_property_lists
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
    pub fn with_sdf_data_field(mut self, key: impl Into<String>, value: impl Into<String>) -> Self {
        self.sdf_data_fields.push((key.into(), value.into()));
        self
    }

    #[must_use]
    pub fn with_sdf_property_list(mut self, property_list: SdfPropertyList) -> Self {
        self.sdf_property_lists.push(property_list);
        self
    }

    #[allow(dead_code)]
    pub fn set_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
        // RDKit✔️✔️: d_props.setVal(key, val);
        // A non-computed write does not remove an existing computed marker.
        self.props.insert(key.into(), value.into());
    }

    #[allow(dead_code)]
    pub fn set_computed_prop(&mut self, key: impl Into<String>, value: impl Into<String>) {
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

    pub fn remap_topology(&mut self, atom_new_to_old: &[Option<AtomId>], bond_new_to_old: &[Option<BondId>]) {
        self.sdf_property_lists = self
            .sdf_property_lists
            .iter()
            .map(|property_list| property_list.remapped_topology(atom_new_to_old, bond_new_to_old))
            .collect();
    }
}

impl SdfPropertyList {
    fn remapped_topology(&self, atom_new_to_old: &[Option<AtomId>], bond_new_to_old: &[Option<BondId>]) -> Self {
        let values = match self.target {
            SdfPropertyListTarget::Atom => atom_new_to_old
                .iter()
                .map(|old_row| old_row.and_then(|row| self.values.get(row.index()).cloned().flatten()))
                .collect(),
            SdfPropertyListTarget::Bond => bond_new_to_old
                .iter()
                .map(|old_row| old_row.and_then(|row| self.values.get(row.index()).cloned().flatten()))
                .collect(),
        };
        Self {
            target: self.target,
            name: self.name.clone(),
            values,
        }
    }
}
