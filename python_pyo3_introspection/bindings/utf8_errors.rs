//! Concrete standard UTF-8 failures used by native IO and stereo errors.
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module = "cosmolkit", extends=PyValueError)]
pub struct Utf8Error {
    inner: std::str::Utf8Error,
}

#[pymethods]
impl Utf8Error {
    #[getter]
    fn domain(&self) -> &'static str {
        "std"
    }
    #[getter]
    fn valid_up_to(&self) -> usize {
        self.inner.valid_up_to()
    }
    #[getter]
    fn error_len(&self) -> Option<usize> {
        self.inner.error_len()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

impl Utf8Error {
    pub(crate) fn from_rust(inner: std::str::Utf8Error) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module = "cosmolkit", extends=PyValueError)]
pub struct FromUtf8Error {
    inner: std::string::FromUtf8Error,
}

#[pymethods]
impl FromUtf8Error {
    #[getter]
    fn domain(&self) -> &'static str {
        "std"
    }
    fn as_bytes<'py>(&self, py: Python<'py>) -> Bound<'py, PyBytes> {
        PyBytes::new(py, self.inner.as_bytes())
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn utf8_error(&self, py: Python<'_>) -> PyResult<Py<Utf8Error>> {
        Utf8Error::from_rust(self.inner.utf8_error())
            .value(py)
            .extract()
            .map_err(Into::into)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

impl FromUtf8Error {
    pub(crate) fn from_rust(inner: std::string::FromUtf8Error) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::error::Error;

    #[test]
    fn real_owned_utf8_failures_keep_bytes_and_payload_not_cause() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            for (bytes, error_len) in [(vec![0, 255], Some(1_usize)), (vec![0, 0xe2, 0x82], None)] {
                let native = String::from_utf8(bytes.clone()).unwrap_err();
                assert!(native.source().is_none());
                assert!(native.utf8_error().source().is_none());
                let display = native.to_string();
                let error = FromUtf8Error::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(error.cause(py).is_none());
                assert_eq!(
                    error
                        .value(py)
                        .call_method0("as_bytes")?
                        .cast::<PyBytes>()?
                        .as_bytes(),
                    bytes
                );
                let detail = error.value(py).getattr("utf8_error")?;
                assert!(detail.is_instance_of::<Utf8Error>());
                assert_eq!(detail.getattr("valid_up_to")?.extract::<usize>()?, 1);
                assert_eq!(
                    detail.getattr("error_len")?.extract::<Option<usize>>()?,
                    error_len
                );
                assert_eq!(detail.str()?.to_str()?, display);
                assert!(detail.getattr("__cause__")?.is_none());
                assert!(
                    detail
                        .setattr("error_len", 3)
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
                );
                assert!(
                    error
                        .value(py)
                        .setattr("utf8_error", py.None())
                        .unwrap_err()
                        .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
                );
                assert_eq!(
                    error.value(py).getattr("domain")?.extract::<String>()?,
                    "std"
                );
            }
            Ok(())
        })
        .unwrap();
    }
}
