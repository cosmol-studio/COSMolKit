use crate::chemical_values::enum_value;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(PropertyValueKind, { String, Int, UInt, IntVector, StringVector, Double, Bool });

crate::error_framework::error_projection! {record
#[doc = r#"{"role": "error_record"}"#]
pub struct MissingPropertyError {
    inner: ck::MissingPropertyError,
}
}
crate::error_framework::error_projection! {methods impl MissingPropertyError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    /// Name of the required property that was not present.
    fn key(&self, py: Python<'_>) -> PyResult<String> {
        decode_source_text(py, self.inner.key())
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl MissingPropertyError {
    pub(crate) fn from_rust(inner: ck::MissingPropertyError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role": "error_record"}"#]
pub struct PropertyValueError {
    inner: ck::PropertyValueError,
}
}
crate::error_framework::error_projection! {methods impl PropertyValueError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "property"
    }
    #[getter]
    /// Property value kind required by the attempted typed read.
    fn expected(&self) -> PropertyValueKind {
        self.inner.expected().into()
    }
    #[getter]
    /// Property value kind actually stored.
    fn actual(&self) -> PropertyValueKind {
        self.inner.actual().into()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl PropertyValueError {
    pub(crate) fn from_rust(inner: ck::PropertyValueError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Explicitly tagged molecular property value. Typed reads require the matching stored variant and do not coerce other kinds.
    pub struct PropertyValue {
        pub(crate) inner: ck::PropertyValue,
    }
}
#[pymethods]
impl PropertyValue {
    // Enum variant construction is representation, not numeric coercion. Each
    // entry constructs exactly the named canonical Rust discriminant.
    #[staticmethod]
    #[pyo3(name = "String")]
    /// Create a String property value from Python text.
    fn string(value: String) -> Self {
        Self {
            inner: ck::PropertyValue::String(value.into()),
        }
    }
    #[staticmethod]
    #[pyo3(name = "Int")]
    /// Create a signed 32-bit integer property value.
    fn int(value: i32) -> Self {
        Self {
            inner: ck::PropertyValue::Int(value),
        }
    }
    #[staticmethod]
    #[pyo3(name = "UInt")]
    /// Create an unsigned 32-bit integer property value.
    fn uint(value: u32) -> Self {
        Self {
            inner: ck::PropertyValue::UInt(value),
        }
    }
    #[staticmethod]
    #[pyo3(name = "IntVector")]
    /// Create a vector of signed 32-bit integer property values.
    fn int_vector(value: Vec<i32>) -> Self {
        Self {
            inner: ck::PropertyValue::IntVector(value),
        }
    }
    #[staticmethod]
    #[pyo3(name = "StringVector")]
    /// Create a vector of string property values.
    fn string_vector(value: Vec<String>) -> Self {
        Self {
            inner: ck::PropertyValue::StringVector(value.into_iter().map(Into::into).collect()),
        }
    }
    #[staticmethod]
    #[pyo3(name = "Double")]
    /// Create a double-precision floating-point property value.
    fn double(value: f64) -> Self {
        Self {
            inner: ck::PropertyValue::Double(value),
        }
    }
    #[staticmethod]
    #[pyo3(name = "Bool")]
    /// Create a Boolean property value.
    fn bool(value: bool) -> Self {
        Self {
            inner: ck::PropertyValue::Bool(value),
        }
    }
    #[getter]
    /// Stored property variant, without numeric or textual coercion.
    fn kind(&self) -> PropertyValueKind {
        self.inner.kind().into()
    }
    /// Read a String property as text; reject a different stored kind.
    fn as_string(&self, py: Python<'_>) -> PyResult<String> {
        let text = self
            .inner
            .as_string()
            .map_err(PropertyValueError::from_rust)?;
        decode_source_text(py, text)
    }
    /// Read a signed integer property without coercing another stored kind.
    fn as_int(&self) -> PyResult<i32> {
        self.inner.as_int().map_err(PropertyValueError::from_rust)
    }
    /// Read an unsigned integer property without coercing another stored kind.
    fn as_uint(&self) -> PyResult<u32> {
        self.inner.as_uint().map_err(PropertyValueError::from_rust)
    }
    /// Read an integer-vector property as an independent list; reject a different kind.
    fn as_int_vector(&self) -> PyResult<Vec<i32>> {
        self.inner
            .as_int_vector()
            .map(<[i32]>::to_vec)
            .map_err(PropertyValueError::from_rust)
    }
    /// Read a Double value; raise PropertyValueError for a different stored variant.
    fn as_double(&self) -> PyResult<f64> {
        self.inner
            .as_double()
            .map_err(PropertyValueError::from_rust)
    }
    /// Read a Bool value; raise PropertyValueError for a different stored variant.
    fn as_bool(&self) -> PyResult<bool> {
        self.inner.as_bool().map_err(PropertyValueError::from_rust)
    }
}

pub(crate) fn decode_source_text(py: Python<'_>, text: &ck::PropertyText) -> PyResult<String> {
    // Same strict UTF-8 boundary as the production Python projection. Invalid
    // source bytes remain a UnicodeDecodeError, never replacement characters.
    pyo3::types::PyBytes::new(py, text.as_bytes())
        .call_method1("decode", ("utf-8", "strict"))?
        .extract()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn actual_missing_property_retains_record_key_and_source_absence() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let atom = ck::Atom::from_spec(ck::AtomId::new(0), ck::AtomSpec::new(ck::Element::C));
            let native = atom.prop_required("tracking\0id").unwrap_err();
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = MissingPropertyError::from_rust(native);
            assert!(error.is_instance(py, <MissingPropertyError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(
                error.value(py).getattr("key")?.extract::<String>()?,
                "tracking\0id"
            );
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.cause(py).is_none());
            assert!(
                error
                    .value(py)
                    .setattr("key", "other")
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            Ok(())
        })
        .unwrap();
    }
}

/// Counted native byte text is exposed as Python text only by strict decoding.
/// This shares the existing explicit UTF-8 boundary; no lossy replacement.
pub(crate) struct PropertyTextField;
crate::configuration::type_projection! {
    impl Container for PropertyTextField {
        type Native = cosmolkit::PropertyText;
        const TYPE: pyo3::inspect::PyStaticExpr = pyo3::type_hint_identifier!("builtins", "str");
        const INPUT: pyo3::inspect::PyStaticExpr = Self::TYPE;
        fn extract(value: pyo3::Borrowed<'_, '_, pyo3::PyAny>) -> pyo3::PyResult<Self::Native> {
            Ok(value.extract::<String>()?.into())
        }
        fn output<'py>(
            py: pyo3::Python<'py>,
            value: &Self::Native,
        ) -> pyo3::PyResult<pyo3::Bound<'py, pyo3::PyAny>> {
            use pyo3::IntoPyObjectExt;
            crate::property_values::decode_source_text(py, value)?.into_bound_py_any(py)
        }
    }
}
crate::error_framework::error_projection!(AtomPropertyReadError, AtomPropertyReadErrorKind, ck::AtomPropertyReadError, "property", {
    Operation(_), Type(_), Utf8(_),
}, cause = |error: &ck::AtomPropertyReadError| Some(match error {
    ck::AtomPropertyReadError::Operation(source) => crate::operation_errors::OperationError::from_rust(source.clone()),
    ck::AtomPropertyReadError::Type(source) => PropertyValueError::from_rust(source.clone()),
    ck::AtomPropertyReadError::Utf8(source) => pyo3::exceptions::PyUnicodeError::new_err(source.to_string()),
}), methods = {
    #[getter]
    /// Number of valid UTF-8 bytes before the decoding failure, or None for non-decoding errors.
    fn valid_up_to(&self) -> Option<usize> {
        match &self.inner { ck::AtomPropertyReadError::Utf8(e) => Some(e.valid_up_to()), _ => None }
    }
    #[getter]
    /// Length of the invalid UTF-8 sequence, or None for incomplete input or a different error kind.
    fn error_len(&self) -> Option<usize> {
        match &self.inner { ck::AtomPropertyReadError::Utf8(e) => e.error_len(), _ => None }
    }
});
