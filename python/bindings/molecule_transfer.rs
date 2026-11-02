//! Same-version value transfer delegates entirely to the canonical Rust API.
use crate::binary_archive::PickleError;
use crate::error_framework::{ErrorObject, error_projection};
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(MoleculeTransferError, MoleculeTransferErrorKind, ck::MoleculeTransferError, "serialization", {
    VersionMismatch { expected: String, actual: String }, Codec(_),
}, cause = |inner: &ck::MoleculeTransferError| match inner {
    ck::MoleculeTransferError::Codec(error) => Some(PickleError::from_rust(error.clone())),
    ck::MoleculeTransferError::VersionMismatch { .. } => None,
}, source = Option<ErrorObject<PickleError>>);

#[pymethods]
impl Molecule {
    /// :meta private:
    ///
    /// Encode an owned value for another CK module of exactly the same version.
    /// Preserves topology, metadata and every coordinate set, but not derived
    /// caches. Does not change this molecule or replace the persistent archive.
    fn to_transfer_binary<'py>(
        &self,
        py: Python<'py>,
    ) -> PyResult<Bound<'py, pyo3::types::PyBytes>> {
        self.inner
            .to_transfer_binary()
            .map(|data| pyo3::types::PyBytes::new(py, &data))
            .map_err(MoleculeTransferError::from_rust)
    }

    #[staticmethod]
    /// :meta private:
    ///
    /// Restore a molecular value owned by this CK module from transfer bytes.
    /// Rejects any different complete CK version before molecular decoding.
    fn from_transfer_binary(data: &[u8]) -> PyResult<Self> {
        ck::Molecule::from_transfer_binary(data)
            .map(|inner| Self { inner })
            .map_err(MoleculeTransferError::from_rust)
    }
}
