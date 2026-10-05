//! Read-only projections of existing canonical property values.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::collections::BTreeMap;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    Double,
    Bool,
}
impl From<ck::PropertyValueKind> for PropertyValueKind {
    fn from(value: ck::PropertyValueKind) -> Self {
        match value {
            ck::PropertyValueKind::String => Self::String,
            ck::PropertyValueKind::Int => Self::Int,
            ck::PropertyValueKind::UInt => Self::UInt,
            ck::PropertyValueKind::IntVector => Self::IntVector,
            ck::PropertyValueKind::Double => Self::Double,
            ck::PropertyValueKind::Bool => Self::Bool,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PropertyValue {
    inner: ck::PropertyValue,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PropertyValue {
    fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
    fn as_string(&self, py: Python<'_>) -> PyResult<&str> {
        self.inner
            .as_string()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    fn as_int(&self, py: Python<'_>) -> PyResult<i32> {
        self.inner
            .as_int()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    fn as_uint(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .as_uint()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    fn as_int_vector(&self, py: Python<'_>) -> PyResult<Vec<i32>> {
        self.inner
            .as_int_vector()
            .map(<[i32]>::to_vec)
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    fn as_double(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .as_double()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    fn as_bool(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .as_bool()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum SdfPropertyListTarget {
    Atom,
    Bond,
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfPropertyList {
    inner: ck::SdfPropertyList,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfPropertyList {
    fn target(&self) -> SdfPropertyListTarget {
        match self.inner.target() {
            ck::SdfPropertyListTarget::Atom => SdfPropertyListTarget::Atom,
            ck::SdfPropertyListTarget::Bond => SdfPropertyListTarget::Bond,
        }
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn values(&self) -> Vec<Option<PropertyValue>> {
        self.inner
            .values()
            .iter()
            .map(|v| v.clone().map(|inner| PropertyValue { inner }))
            .collect()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MoleculeProperties {
    pub(crate) inner: ck::MoleculeProperties,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeProperties {
    fn name(&self) -> Option<&str> {
        self.inner.name()
    }
    fn props(&self) -> BTreeMap<String, String> {
        self.inner.props().clone()
    }
    fn prop(&self, key: &str) -> Option<&str> {
        self.inner.prop(key)
    }
    fn is_prop_computed(&self, key: &str) -> bool {
        self.inner.is_prop_computed(key)
    }
    fn computed_prop_names(&self) -> Vec<String> {
        self.inner.computed_prop_names().iter().cloned().collect()
    }
    fn sdf_data_fields(&self) -> Vec<(String, String)> {
        self.inner.sdf_data_fields().to_vec()
    }
    fn sdf_property_lists(&self) -> Vec<SdfPropertyList> {
        self.inner
            .sdf_property_lists()
            .iter()
            .cloned()
            .map(|inner| SdfPropertyList { inner })
            .collect()
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<PropertyValueKind>()?;
    module.add_class::<PropertyValue>()?;
    module.add_class::<MoleculeProperties>()?;
    module.add_class::<SdfPropertyList>()?;
    module.add_class::<SdfPropertyListTarget>()?;
    Ok(())
}
