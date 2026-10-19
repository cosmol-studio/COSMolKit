//! Concrete allocation failures; no unstable kind inference or message parsing.
use crate::error_framework::error_projection;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

// Like the existing ParseFloatError projection, this typed standard leaf owns
// its actual error and exposes only stable public information. kind() is an
// unstable Rust API, not a field to reconstruct from Display or Debug.
#[doc = r#"{"role":"error_leaf"}"#]
#[pyclass(module = "cosmolkit", extends=PyValueError)]
pub struct TryReserveError {
    inner: std::collections::TryReserveError,
}

#[pymethods]
impl TryReserveError {
    #[getter]
    fn domain(&self) -> &'static str {
        "std"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

impl TryReserveError {
    pub(crate) fn from_rust(inner: std::collections::TryReserveError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

error_projection!(retained UnsignedStreamArrayError, UnsignedStreamArrayErrorKind, ck::UnsignedStreamArrayError, "property", {
    MissingFirstValue, Allocation(_),
}, cause = |inner: &crate::retained_error::RetainedError<ck::UnsignedStreamArrayError>| match inner.as_ref() {
    ck::UnsignedStreamArrayError::MissingFirstValue => None,
    ck::UnsignedStreamArrayError::Allocation(value) => Some(TryReserveError::from_rust(value.clone())),
}, source = Option<Py<TryReserveError>>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_capacity_failure_keeps_exact_allocation_leaf_without_exhausting_memory() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            // u8 capacity exceeds isize::MAX, so this rejects the layout before
            // an allocator request. Do not try to exhaust available memory.
            let mut values = Vec::<u8>::new();
            let native_leaf = values.try_reserve(usize::MAX).unwrap_err();
            assert_eq!(values.capacity(), 0);
            assert!(native_leaf.source().is_none());
            let leaf_display = native_leaf.to_string();
            let native = ck::UnsignedStreamArrayError::Allocation(native_leaf);
            assert!(
                native
                    .source()
                    .unwrap()
                    .is::<std::collections::TryReserveError>()
            );
            let drawing = ck::DrawingError::VariationArray(native);
            assert!(
                drawing
                    .source()
                    .unwrap()
                    .is::<ck::UnsignedStreamArrayError>()
            );
            let ck::DrawingError::VariationArray(native) = drawing else {
                unreachable!()
            };
            let display = native.to_string();
            let error = UnsignedStreamArrayError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<TryReserveError>());
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(cause.str()?.to_str()?, leaf_display);
            assert_eq!(cause.getattr("domain")?.extract::<String>()?, "std");
            assert!(cause.getattr("__cause__")?.is_none());
            assert!(!cause.hasattr("kind")? && !cause.hasattr("variant")?);
            assert!(cause.setattr("domain", "different").is_err());
            assert!(error.value(py).setattr("source", py.None()).is_err());

            let native = ck::UnsignedStreamArrayError::MissingFirstValue;
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = UnsignedStreamArrayError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<UnsignedStreamArrayErrorKind>()?
                    == UnsignedStreamArrayErrorKind::MissingFirstValue
            );
            Ok(())
        })
        .unwrap();
    }
}
