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

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct SdfPropertyList {
    pub(crate) inner: ck::SdfPropertyList,
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
    fn target(&self) -> SdfPropertyListTarget {
        self.inner.target().into()
    }
    fn name(&self, py: Python<'_>) -> PyResult<String> {
        decode_source_text(py, self.inner.name())
    }
    fn values(&self) -> Vec<Option<PropertyValue>> {
        self.inner
            .values()
            .iter()
            .map(|value| value.clone().map(|inner| PropertyValue { inner }))
            .collect()
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct MoleculeProperties {
    pub(crate) inner: ck::MoleculeProperties,
}
#[pymethods]
impl MoleculeProperties {
    #[new]
    fn new() -> Self {
        Self {
            inner: ck::MoleculeProperties::default(),
        }
    }
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .name()
            .map(|text| decode_source_text(py, text))
            .transpose()
    }
    fn sdf_data_fields(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        self.inner
            .sdf_data_fields()
            .iter()
            .map(|(key, value)| Ok((decode_source_text(py, key)?, decode_source_text(py, value)?)))
            .collect()
    }
    fn sdf_property_lists(&self) -> Vec<SdfPropertyList> {
        self.inner
            .sdf_property_lists()
            .iter()
            .cloned()
            .map(|inner| SdfPropertyList { inner })
            .collect()
    }
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
    fn prop(&self, key: &str) -> Option<PropertyValue> {
        self.inner
            .prop(key)
            .cloned()
            .map(|inner| PropertyValue { inner })
    }
    fn is_prop_computed(&self, key: &str) -> PyResult<bool> {
        self.inner
            .is_prop_computed(key)
            .map_err(PropertyValueError::from_rust)
    }
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
