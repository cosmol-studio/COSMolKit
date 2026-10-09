//! Python filesystem-protocol extraction, separate from pure path helpers.
use pyo3::prelude::*;

/// Text-valued filesystem inputs for facade methods whose paths are UTF-8.
/// Use Python's filesystem protocol, never an arbitrary object's `str()`.
pub(crate) struct TextPath(String);

impl TextPath {
    pub(crate) fn as_str(&self) -> &str {
        &self.0
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for TextPath {
    type Error = PyErr;

    fn extract(value: Borrowed<'a, 'py, PyAny>) -> PyResult<Self> {
        let path = value
            .py()
            .import("os")?
            .getattr("fspath")?
            .call1((value,))?;
        // String extraction rejects bytes, including byte-valued PathLike.
        path.extract::<String>().map(Self)
    }
}

#[cfg(feature = "stubgen")]
impl pyo3_stub_gen::PyStubType for TextPath {
    fn type_output() -> pyo3_stub_gen::TypeInfo {
        pyo3_stub_gen::TypeInfo::builtin("str")
            | pyo3_stub_gen::TypeInfo::with_module("os.PathLike", "os".into())
    }
}
