//! String transport for identifier constructors; never a chemistry parser.
use pyo3::{
    exceptions::PyTypeError,
    prelude::*,
    types::{PyBytes, PyString},
};

/// The Rust string facade accepts UTF-8 text. Keep existing Python strings
/// borrowed; decode bytes strictly once without lossy replacement or coercion.
pub(crate) struct TextInput<'py>(Bound<'py, PyString>);

impl TextInput<'_> {
    pub(crate) fn as_text(&self) -> PyResult<std::borrow::Cow<'_, str>> {
        self.0.extract()
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for TextInput<'py> {
    type Error = PyErr;
    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        // RDKit✔️✔️: python::extract<std::string> ex(input);
        // RDKit✔️✔️: if (ex.check()) {
        // RDKit✔️✔️:   return ex();
        // RDKit✔️✔️: }
        // rdmolfiles.cpp:55–59; InChI uses the same std::string boundary.
        // Scope: UTF-8 identifiers supported by the Rust &str facade. Python
        // str is retained without copying; bytes use CPython's strict decoder.
        if let Ok(text) = value.cast::<PyString>() {
            return Ok(Self(text.to_owned()));
        }
        if let Ok(bytes) = value.cast::<PyBytes>() {
            return bytes
                .call_method1("decode", ("utf-8", "strict"))?
                .cast_into::<PyString>()
                .map(Self)
                .map_err(Into::into);
        }
        Err(PyTypeError::new_err("expected str or UTF-8 bytes"))
    }
}

#[cfg(feature = "stubgen")]
impl pyo3_stub_gen::PyStubType for TextInput<'_> {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo::builtin("str") | pyo3_stub_gen::TypeInfo::builtin("bytes")
    }
}

#[cfg(all(test, feature = "python-embed-tests"))]
mod tests {
    use super::*;
    #[test]
    fn text_extraction_borrows_strings_and_preserves_counted_utf8() {
        Python::initialize();
        Python::attach(|py| {
            let original = PyString::new(py, "CCO name α\0suffix");
            let text: TextInput<'_> = original.extract().unwrap();
            assert!(text.0.is(&original));
            assert_eq!(text.as_text().unwrap(), "CCO name α\0suffix");
            let bytes = PyBytes::new(py, b"CCO name \xce\xb1\0suffix");
            let decoded: TextInput<'_> = bytes.extract().unwrap();
            assert_eq!(decoded.as_text().unwrap(), text.as_text().unwrap());
            let invalid = PyBytes::new(py, b"CCO\xff")
                .extract::<TextInput<'_>>()
                .err()
                .unwrap();
            assert!(invalid.is_instance_of::<pyo3::exceptions::PyUnicodeDecodeError>(py));
            assert!(
                py.None()
                    .bind(py)
                    .extract::<TextInput<'_>>()
                    .err()
                    .unwrap()
                    .is_instance_of::<PyTypeError>(py)
            );
        });
    }
}
