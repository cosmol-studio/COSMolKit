//! Read-only projections of existing canonical property values.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::collections::BTreeMap;

/// The modeled source value kinds supported by atom and bond properties.
///
/// Declared values: ``String``, ``Int``, ``UInt``, ``IntVector``, ``StringVector``, ``Double``, ``Bool``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    StringVector,
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
            ck::PropertyValueKind::StringVector => Self::StringVector,
            ck::PropertyValueKind::Double => Self::Double,
            ck::PropertyValueKind::Bool => Self::Bool,
        }
    }
}
/// Typed property value; as_* accessors require the corresponding value kind and raise PropertyValueError on a mismatch.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct PropertyValue {
    pub(crate) inner: ck::PropertyValue,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl PropertyValue {
    /// Return the PropertyValueKind discriminant; it determines the valid as_* accessor.
    fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
    /// Return the string value; a different PropertyValueKind raises PropertyValueError.
    fn as_string(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .as_string()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, text))
    }
    /// Return the signed integer value; a different PropertyValueKind raises PropertyValueError.
    fn as_int(&self, py: Python<'_>) -> PyResult<i32> {
        self.inner
            .as_int()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    /// Return the unsigned integer value; a different PropertyValueKind raises PropertyValueError.
    fn as_uint(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .as_uint()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    /// Return a list of signed integer values; a different PropertyValueKind raises PropertyValueError.
    fn as_int_vector(&self, py: Python<'_>) -> PyResult<Vec<i32>> {
        self.inner
            .as_int_vector()
            .map(<[i32]>::to_vec)
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    /// Return the floating-point value; a different PropertyValueKind raises PropertyValueError.
    fn as_double(&self, py: Python<'_>) -> PyResult<f64> {
        self.inner
            .as_double()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    /// Return the boolean value; a different PropertyValueKind raises PropertyValueError.
    fn as_bool(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .as_bool()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
}
/// Graph entity to which an interpreted SDF property list applies: Atom or Bond.
///
/// Declared values: ``Atom``, ``Bond``.
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum SdfPropertyListTarget {
    Atom,
    Bond,
}
/// Interpreted SDF atom/bond property list. Entries follow target graph order; None marks an explicit missing value.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SdfPropertyList {
    inner: ck::SdfPropertyList,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SdfPropertyList {
    /// Target identifier of this template attachment or SDF property list.
    fn target(&self) -> SdfPropertyListTarget {
        match self.inner.target() {
            ck::SdfPropertyListTarget::Atom => SdfPropertyListTarget::Atom,
            ck::SdfPropertyListTarget::Bond => SdfPropertyListTarget::Bond,
        }
    }
    /// Stored name of this value.
    fn name(&self, py: Python<'_>) -> PyResult<String> {
        crate::canonical_sdf::decode_source_text(py, self.inner.name())
    }
    /// Return typed values in target atom/bond order; None denotes an explicit missing value.
    fn values(&self) -> Vec<Option<PropertyValue>> {
        self.inner
            .values()
            .iter()
            .map(|v| v.clone().map(|inner| PropertyValue { inner }))
            .collect()
    }
}
/// Read-only molecule name, typed properties, computed-property markers and ordered SDF data fields.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MoleculeProperties {
    pub(crate) inner: ck::MoleculeProperties,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MoleculeProperties {
    /// Stored name of this value.
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .name()
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .transpose()
    }
    /// Stored typed properties; returned values do not provide mutable access to the owning molecule.
    fn props(&self, py: Python<'_>) -> PyResult<BTreeMap<String, PropertyValue>> {
        self.inner
            .props()
            .iter()
            .map(|(key, value)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    PropertyValue {
                        inner: value.clone(),
                    },
                ))
            })
            .collect()
    }
    /// Return the typed property for the requested key, or None if absent.
    fn prop(&self, key: &str) -> Option<PropertyValue> {
        self.inner
            .prop(key)
            .cloned()
            .map(|inner| PropertyValue { inner })
    }
    /// Returns whether a property is registered as computed state.
    fn is_prop_computed(&self, py: Python<'_>, key: &str) -> PyResult<bool> {
        self.inner
            .is_prop_computed(key)
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))
    }
    /// Property names marked as computed.
    fn computed_prop_names(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        let names = self
            .inner
            .computed_prop_names()
            .map_err(|e| crate::canonical_atom_bond::property_pyerr(py, e))?;
        names
            .unwrap_or_default()
            .iter()
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .collect()
    }
    /// Return ordered SDF (name, value) string pairs, preserving duplicate names.
    fn sdf_data_fields(&self, py: Python<'_>) -> PyResult<Vec<(String, String)>> {
        self.inner
            .sdf_data_fields()
            .iter()
            .map(|(key, value)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    crate::canonical_sdf::decode_source_text(py, value)?,
                ))
            })
            .collect()
    }
    /// Interpreted atom/bond SDF property lists in stored order.
    fn sdf_property_lists(&self) -> Vec<SdfPropertyList> {
        self.inner
            .sdf_property_lists()
            .iter()
            .cloned()
            .map(|inner| SdfPropertyList { inner })
            .collect()
    }
}

pyo3::create_exception!(
    cosmolkit,
    PropertyStringError,
    pyo3::exceptions::PyValueError,
    "A molecular property could not be converted to the requested textual representation."
);

/// Return the property's lexical text as a Python string without changing its stored value or type.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn property_value_to_text(py: Python<'_>, value: &PropertyValue) -> PyResult<String> {
    let text = ck::property_value_to_text(&value.inner)
        .map_err(|source| property_string_pyerr(py, source))?;
    crate::canonical_sdf::decode_source_text(py, &text)
}

pub(crate) fn property_string_pyerr(py: Python<'_>, source: ck::PropertyStringError) -> PyErr {
    let error = crate::canonical_values::annotate(
        py,
        PropertyStringError::new_err(source.to_string()),
        "properties",
        "UnsupportedKind",
        &source,
    );
    let value_kind = PropertyValueKind::from(source.kind());
    let payload = || -> PyResult<()> {
        error.value(py).setattr("error_kind", "UnsupportedKind")?;
        error.value(py).setattr("value_kind", value_kind)?;
        // Keep the registered class accessor; the generic annotation's string
        // category is retained separately rather than shadowing kind().
        error.value(py).delattr("kind")
    };
    match payload() {
        Ok(()) => error,
        Err(error) => error,
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    crate::canonical_error_accessors::attach(
        module.py().get_type::<PropertyStringError>().as_any(),
        &[("kind", "value_kind")],
    )?;
    module.add(
        "PropertyStringError",
        module.py().get_type::<PropertyStringError>(),
    )?;
    module.add_function(wrap_pyfunction!(property_value_to_text, module)?)?;
    module.add_class::<PropertyValueKind>()?;
    module.add_class::<PropertyValue>()?;
    module.add_class::<MoleculeProperties>()?;
    module.add_class::<SdfPropertyList>()?;
    module.add_class::<SdfPropertyListTarget>()?;
    Ok(())
}
