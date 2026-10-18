use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum QueryCompileErrorKind {
    InvalidGraph,
}

#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct QueryCompileError {
    inner: ck::QueryCompileError,
}
#[pymethods]
impl QueryCompileError {
    #[getter]
    fn variant(&self) -> QueryCompileErrorKind {
        match &self.inner {
            ck::QueryCompileError::InvalidGraph(_) => QueryCompileErrorKind::InvalidGraph,
        }
    }
    #[getter]
    fn message(&self) -> &str {
        match &self.inner {
            ck::QueryCompileError::InvalidGraph(message) => message,
        }
    }
    #[getter]
    fn domain(&self) -> &'static str {
        "search"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl QueryCompileError {
    pub(crate) fn from_rust(inner: ck::QueryCompileError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(error) => PyErr::from_value(error.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum AvalonEngineErrorKind {
    InvalidArguments,
    AvalonConversion,
}

#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct AvalonEngineError {
    inner: ck::AvalonEngineError,
}
#[pymethods]
impl AvalonEngineError {
    #[getter]
    fn variant(&self) -> AvalonEngineErrorKind {
        match &self.inner {
            ck::AvalonEngineError::InvalidArguments { reason: _ } => {
                AvalonEngineErrorKind::InvalidArguments
            }
            ck::AvalonEngineError::AvalonConversion { reason: _ } => {
                AvalonEngineErrorKind::AvalonConversion
            }
        }
    }
    #[getter]
    fn reason(&self) -> &str {
        match &self.inner {
            ck::AvalonEngineError::InvalidArguments { reason } => reason,
            ck::AvalonEngineError::AvalonConversion { reason } => reason,
        }
    }
    #[getter]
    fn domain(&self) -> &'static str {
        "fingerprint"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl AvalonEngineError {
    pub(crate) fn from_rust(inner: ck::AvalonEngineError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(error) => PyErr::from_value(error.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn registered_string_errors_preserve_payload_and_native_validation() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::QueryCompileError::InvalidGraph("owned\0context".to_owned());
            let expected = native.to_string();
            let error = QueryCompileError::from_rust(native);
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                "owned\0context"
            );
            assert_eq!(error.value(py).str()?.to_str()?, expected);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<QueryCompileErrorKind>()?
                    == QueryCompileErrorKind::InvalidGraph
            );
            assert!(error.cause(py).is_none());
            let native = ck::AvalonFingerprintParams {
                n_bits: 0,
                ..Default::default()
            }
            .validate()
            .unwrap_err();
            let expected = native.to_string();
            let error = AvalonEngineError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<AvalonEngineErrorKind>()?
                    == AvalonEngineErrorKind::InvalidArguments
            );
            assert_eq!(
                error.value(py).getattr("reason")?.extract::<String>()?,
                "Avalon n_bits must be at least 8"
            );
            assert_eq!(error.value(py).str()?.to_str()?, expected);
            assert!(error.cause(py).is_none());
            let error = AvalonEngineError::from_rust(ck::AvalonEngineError::AvalonConversion {
                reason: "owned conversion context".to_owned(),
            });
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<AvalonEngineErrorKind>()?
                    == AvalonEngineErrorKind::AvalonConversion
            );
            assert_eq!(
                error.value(py).getattr("reason")?.extract::<String>()?,
                "owned conversion context"
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
