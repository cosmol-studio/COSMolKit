use crate::chemical_values::{EnumValue, enum_value};
use crate::property_values::{PropertyValue, PropertyValueError, decode_source_text};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

enum_value!(SdfPropertyListTarget, { Atom, Bond });
#[derive(FromPyObject, IntoPyObject)]
pub enum SdfPropertyListTargetInput {
    Enum(SdfPropertyListTarget),
    Name(String),
}
impl SdfPropertyListTargetInput {
    fn into_rust(self) -> PyResult<ck::SdfPropertyListTarget> {
        match self {
            Self::Enum(value) => Ok(value.into()),
            Self::Name(name) => SdfPropertyListTarget::from_name(&name).map(Into::into),
        }
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Named, ordered per-atom or per-bond property values, with explicit missing entries.
    pub struct SdfPropertyList {
        pub(crate) inner: ck::SdfPropertyList,
    }
}
#[pymethods]
impl SdfPropertyList {
    #[new]
    fn new(
        target: SdfPropertyListTargetInput,
        name: String,
        values: Vec<Option<PropertyValue>>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::SdfPropertyList::new(
                target.into_rust()?,
                name,
                values
                    .into_iter()
                    .map(|value| value.map(|value| value.inner))
                    .collect(),
            ),
        })
    }
    #[getter]
    /// Whether the list assigns values to atoms or bonds.
    fn target(&self) -> SdfPropertyListTarget {
        self.inner.target().into()
    }
    #[getter]
    /// Stored property-list name or molecular name; an absent molecular name is None.
    fn name(&self, py: Python<'_>) -> PyResult<String> {
        decode_source_text(py, self.inner.name())
    }
    #[getter]
    /// Independent property-value list in target-index order, preserving missing entries as None.
    fn values(&self) -> Vec<Option<PropertyValue>> {
        self.inner
            .values()
            .iter()
            .map(|value| value.clone().map(|inner| PropertyValue { inner }))
            .collect()
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Detached molecular metadata, SDF fields and tagged properties. Returned containers are independent copies.
    pub struct MoleculeProperties {
        pub(crate) inner: ck::MoleculeProperties,
    }
}
#[pymethods]
impl MoleculeProperties {
    #[new]
    fn new() -> Self {
        Self {
            inner: ck::MoleculeProperties::default(),
        }
    }
    #[getter]
    /// Stored property-list name or molecular name; an absent molecular name is None.
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .name()
            .map(|text| decode_source_text(py, text))
            .transpose()
    }
    #[getter]
    /// Ordered SDF name/value pairs, preserving repeated names.
    fn sdf_data_fields(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        self.inner
            .sdf_data_fields()
            .iter()
            .map(|(key, value)| Ok((decode_source_text(py, key)?, decode_source_text(py, value)?)))
            .collect()
    }
    #[getter]
    /// Independent copies of SDF per-atom or per-bond property lists.
    fn sdf_property_lists(&self) -> Vec<SdfPropertyList> {
        self.inner
            .sdf_property_lists()
            .iter()
            .cloned()
            .map(|inner| SdfPropertyList { inner })
            .collect()
    }
    #[getter]
    /// Independent mapping of property names to explicitly tagged values.
    fn props(&self, py: Python<'_>) -> PyResult<BTreeMap<String, PropertyValue>> {
        self.inner
            .props()
            .iter()
            .map(|(key, value)| {
                Ok((
                    decode_source_text(py, key)?,
                    PropertyValue {
                        inner: value.clone(),
                    },
                ))
            })
            .collect()
    }
    /// Look up a property by name; return None when absent without coercing its stored type.
    fn prop(&self, key: &str) -> Option<PropertyValue> {
        self.inner
            .prop(key)
            .cloned()
            .map(|inner| PropertyValue { inner })
    }
    /// Whether the property is marked as computed; malformed computed-property metadata raises PropertyValueError.
    fn is_prop_computed(&self, key: &str) -> PyResult<bool> {
        self.inner
            .is_prop_computed(key)
            .map_err(PropertyValueError::from_rust)
    }
    #[getter]
    /// Names recorded as computed properties, or None when no such metadata exists.
    fn computed_prop_names(&self, py: Python<'_>) -> PyResult<Option<Vec<String>>> {
        self.inner
            .computed_prop_names()
            .map_err(PropertyValueError::from_rust)?
            .map(|names| {
                names
                    .iter()
                    .map(|text| decode_source_text(py, text))
                    .collect()
            })
            .transpose()
    }
}
