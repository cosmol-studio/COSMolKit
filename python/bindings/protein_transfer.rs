//! Viewer transport delegates to the canonical Rust Protein API.
use crate::error_framework::{ErrorObject, error_projection};
use crate::protein_values::{Protein, ProteinProjectionError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::{error::Error, sync::Arc};

error_projection!(ProteinTransferError, ProteinTransferErrorKind, ck::ProteinTransferError, "bio", {
    VersionMismatch { expected: String, actual: String },
    InvalidData { message: String },
    Projection(_),
}, cause = |inner: &ck::ProteinTransferError| match inner {
    ck::ProteinTransferError::Projection(error) => Some(ProteinProjectionError::from_rust(error.clone())),
    _ => None,
}, source = Option<ErrorObject<ProteinProjectionError>>);

#[pymethods]
impl Protein {
    /// :meta private:
    ///
    /// Transfer the viewer hierarchy, atom values, coordinates and structural annotations.
    fn to_transfer_binary<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyBytes>> {
        self.inner
            .to_transfer_binary()
            .map(|data| pyo3::types::PyBytes::new(py, &data))
            .map_err(ProteinTransferError::from_rust)
    }

    #[staticmethod]
    /// :meta private:
    ///
    /// Restore an independently owned viewer value from the exact same CK version.
    fn from_transfer_binary(data: &[u8]) -> PyResult<Self> {
        ck::Protein::from_transfer_binary(data)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(ProteinTransferError::from_rust)
    }
}
