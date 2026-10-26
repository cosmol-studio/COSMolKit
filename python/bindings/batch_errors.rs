//! Ordered real batch failures and reports; preserve the native dynamic source owner.
use crate::{
    error_framework::{Source, error_projection},
    io_source_errors::{ErrorRef, standard_source},
    retained_error::RetainedError,
};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::{error::Error, path::PathBuf};
pub(crate) type BatchSource=Source<crate::error_framework::ErrorObject<crate::operation_errors::SmilesError>,Source<crate::error_framework::ErrorObject<crate::operation_errors::OperationError>,Source<crate::error_framework::ErrorObject<crate::molecular_io_errors::MolecularIoError>,Source<crate::error_framework::ErrorObject<crate::sdf_records::SdfError>,Source<crate::error_framework::ErrorObject<crate::molecular_io_errors::MolWriteError>,Source<crate::error_framework::ErrorObject<crate::drawing_errors::DrawingError>,Source<crate::error_framework::ErrorObject<crate::drawing_errors::DrawingWriteError>,Source<crate::error_framework::ErrorObject<crate::batch_errors::BatchImageError>,Source<crate::error_framework::ErrorObject<crate::smiles_writer::SmilesWriteError>,Source<crate::error_framework::ErrorObject<crate::conformer_execution::ConformerRunError>,Source<crate::error_framework::ErrorObject<crate::morgan_errors::AtomPairReadError>,Source<crate::error_framework::ErrorObject<crate::morgan_errors::MorganReadError>,Source<crate::error_framework::ErrorObject<crate::morgan_errors::TopologicalTorsionReadError>,Source<crate::error_framework::ErrorObject<crate::morgan_errors::AtomPairError>,Source<crate::error_framework::ErrorObject<crate::mol_post_errors::MolPostError>,Source<crate::error_framework::ErrorObject<crate::sdf_read_errors::SdfReadError>,Source<crate::error_framework::ErrorObject<crate::molecular_io_errors::SdfWriteError>,Source<crate::error_framework::ErrorObject<crate::pattern_fingerprint_errors::PatternFingerprintError>,Source<crate::error_framework::ErrorObject<crate::layered_fingerprint_errors::LayeredFingerprintError>,Source<crate::error_framework::ErrorObject<crate::batch_errors::BatchValidationError>,crate::io_source_errors::StandardSource>>>>>>>>>>>>>>>>>>>>;
fn project(node: ErrorRef) -> PyErr {
    if node.error().is::<ck::SmilesError>() {
        return crate::operation_errors::SmilesError::from_ref(node);
    }
    if node.error().is::<ck::OperationError>() {
        return crate::operation_errors::OperationError::from_ref(node);
    }
    if node.error().is::<ck::MolecularIoError>() {
        return crate::molecular_io_errors::MolecularIoError::from_ref(node);
    }
    if node.error().is::<ck::SdfError>() {
        return crate::sdf_records::SdfError::from_ref(node);
    }
    if node.error().is::<ck::MolWriteError>() {
        return crate::molecular_io_errors::MolWriteError::from_ref(node);
    }
    if node.error().is::<ck::DrawingError>() {
        return crate::drawing_errors::DrawingError::from_ref(node);
    }
    if node.error().is::<ck::DrawingWriteError>() {
        return crate::drawing_errors::DrawingWriteError::from_ref(node);
    }
    if node.error().is::<ck::SmilesWriteError>() {
        return crate::smiles_writer::SmilesWriteError::from_ref(node);
    }
    if node.error().is::<ck::ConformerRunError>() {
        return crate::conformer_execution::ConformerRunError::from_ref(node);
    }
    if node.error().is::<ck::AtomPairReadError>() {
        return crate::morgan_errors::AtomPairReadError::from_ref(node);
    }
    if node.error().is::<ck::MorganReadError>() {
        return crate::morgan_errors::MorganReadError::from_ref(node);
    }
    if node.error().is::<ck::TopologicalTorsionReadError>() {
        return crate::morgan_errors::TopologicalTorsionReadError::from_ref(node);
    }
    if node.error().is::<ck::AtomPairError>() {
        return crate::morgan_errors::AtomPairError::from_ref(node);
    }
    if node.error().is::<ck::BatchImageError>() {
        return BatchImageError::from_ref(node);
    }
    if node.error().is::<ck::BatchValidationError>() {
        return BatchValidationError::from_ref(node);
    }
    if let Some(v) = node.error().downcast_ref::<ck::MolPostError>() {
        return crate::mol_post_errors::MolPostError::from_rust(v.clone());
    }
    if let Some(v) = node.error().downcast_ref::<ck::SdfReadError>() {
        return crate::sdf_read_errors::SdfReadError::from_rust(v.clone());
    }
    if let Some(v) = node.error().downcast_ref::<ck::SdfWriteError>() {
        return crate::molecular_io_errors::SdfWriteError::from_rust(v.clone());
    }
    if let Some(v) = node.error().downcast_ref::<ck::PatternFingerprintError>() {
        return crate::pattern_fingerprint_errors::PatternFingerprintError::from_rust(v.clone());
    }
    if let Some(v) = node.error().downcast_ref::<ck::LayeredFingerprintError>() {
        return crate::layered_fingerprint_errors::LayeredFingerprintError::from_borrowed(v);
    }
    standard_source(node)
}
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct BatchError {
    pub(crate) inner: RetainedError<ck::BatchError>,
}
}
impl BatchError {
    pub(crate) fn from_rust(inner: ck::BatchError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        let cause = node.source().map(project);
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(v) => {
                    let e = crate::error_framework::from_payload(v);
                    e.set_cause(py, cause);
                    e
                }
                Err(e) => e,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl BatchError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "batch"
    }
    #[getter]
    /// Zero-based input record index at which the batch operation failed.
    fn index(&self) -> usize {
        self.inner.index
    }
    #[getter]
    /// Canonical name of the failed batch operation.
    fn operation(&self) -> &'static str {
        self.inner.operation
    }
    #[getter]
    /// Native message describing this record's failure.
    fn message(&self) -> &str {
        &self.inner.message
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<BatchSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct BatchValidationError {
    inner: RetainedError<ck::BatchValidationError>,
}
}
impl BatchValidationError {
    pub(crate) fn from_rust(inner: ck::BatchValidationError) -> PyErr {
        Self::from_ref(RetainedError::new(inner).node().clone())
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            match Py::new(
                py,
                Self {
                    inner: RetainedError::from_ref(node),
                },
            ) {
                Ok(v) => crate::error_framework::from_payload(v),
                Err(e) => e,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl BatchValidationError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "batch"
    }
    #[getter]
    /// Number of invalid batch records.
    fn errors(&self) -> usize {
        self.inner.errors
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    /// Return a copy of the per-record BatchError list, retaining record indices
    /// and the native causes. Editing the list does not modify this exception.
    fn record_errors(&self, py: Python<'_>) -> PyResult<Vec<crate::error_framework::ErrorObject<BatchError>>> {
        self.inner
            .record_errors
            .iter()
            .map(|e| Ok(BatchError::from_rust(e.clone()).value(py).extract()?))
            .collect()
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
error_projection!(retained BatchImageError,BatchImageErrorKind,ck::BatchImageError,"batch",{
InvalidFormat(_),Write(_),
},cause=|inner:&RetainedError<ck::BatchImageError>|inner.node().source().map(crate::drawing_errors::DrawingWriteError::from_ref),
source=Option<crate::error_framework::ErrorObject<crate::drawing_errors::DrawingWriteError>>,methods={
#[getter]#[doc = "Format involved in the failed conversion. Applies to BatchImageError.InvalidFormat. None when this payload is absent."] fn format(&self)->Option<&str>{match self.inner.as_ref(){ck::BatchImageError::InvalidFormat(v)=>Some(v),_=>None}}
});
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Ordered batch export totals and retained per-record failures.
    pub struct BatchExportReport {
        pub(crate) inner: ck::BatchExportReport,
    }
}
#[pymethods]
impl BatchExportReport {
    #[getter]
    /// Total number of records processed by the export.
    fn total(&self) -> usize {
        self.inner.total()
    }
    #[getter]
    /// Number of records successfully written.
    fn success(&self) -> usize {
        self.inner.success()
    }
    #[getter]
    /// Number of records that failed to export.
    fn failed(&self) -> usize {
        self.inner.failed()
    }
    #[getter]
    /// Independent list of per-record export errors, preserving structured causes.
    fn errors(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<crate::error_framework::ErrorObject<BatchError>>> {
        self.inner
            .errors()
            .iter()
            .map(|v| Ok(BatchError::from_rust(v.clone()).value(py).extract()?))
            .collect()
    }
    /// Write the native export report to the supplied filesystem path.
    fn write_report(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        py.detach(|| self.inner.write_report(&path))
            .map_err(BatchValidationError::from_rust)
    }
}
