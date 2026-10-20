//! Ordered real batch failures and reports; preserve the native dynamic source owner.
use crate::{
    error_framework::{Source, error_projection},
    io_source_errors::{ErrorRef, standard_source},
    retained_error::RetainedError,
};
use cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
use std::{error::Error, path::PathBuf};
pub(crate) type BatchSource=Source<Py<crate::operation_errors::SmilesError>,Source<Py<crate::operation_errors::OperationError>,Source<Py<crate::molecular_io_errors::MolecularIoError>,Source<Py<crate::sdf_records::SdfError>,Source<Py<crate::molecular_io_errors::MolWriteError>,Source<Py<crate::drawing_errors::DrawingError>,Source<Py<crate::drawing_errors::DrawingWriteError>,Source<Py<crate::batch_errors::BatchImageError>,Source<Py<crate::smiles_writer::SmilesWriteError>,Source<Py<crate::conformer_execution::ConformerRunError>,Source<Py<crate::morgan_errors::AtomPairReadError>,Source<Py<crate::morgan_errors::MorganReadError>,Source<Py<crate::morgan_errors::TopologicalTorsionReadError>,Source<Py<crate::morgan_errors::AtomPairError>,Source<Py<crate::mol_post_errors::MolPostError>,Source<Py<crate::sdf_read_errors::SdfReadError>,Source<Py<crate::molecular_io_errors::SdfWriteError>,Source<Py<crate::pattern_fingerprint_errors::PatternFingerprintError>,Source<Py<crate::layered_fingerprint_errors::LayeredFingerprintError>,Source<Py<crate::batch_errors::BatchValidationError>,crate::io_source_errors::StandardSource>>>>>>>>>>>>>>>>>>>>;
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
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct BatchError {
    pub(crate) inner: RetainedError<ck::BatchError>,
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
                    let e = PyErr::from_value(v.into_bound(py).into_any());
                    e.set_cause(py, cause);
                    e
                }
                Err(e) => e,
            }
        })
    }
}
#[pymethods]
impl BatchError {
    #[getter]
    fn domain(&self) -> &'static str {
        "batch"
    }
    #[getter]
    fn index(&self) -> usize {
        self.inner.index
    }
    #[getter]
    fn operation(&self) -> &'static str {
        self.inner.operation
    }
    #[getter]
    fn message(&self) -> &str {
        &self.inner.message
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<BatchSource>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="cosmolkit",extends=PyValueError)]
pub struct BatchValidationError {
    inner: RetainedError<ck::BatchValidationError>,
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
                Ok(v) => PyErr::from_value(v.into_bound(py).into_any()),
                Err(e) => e,
            }
        })
    }
}
#[pymethods]
impl BatchValidationError {
    #[getter]
    fn domain(&self) -> &'static str {
        "batch"
    }
    #[getter]
    fn errors(&self) -> usize {
        self.inner.errors
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn record_errors(&self, py: Python<'_>) -> PyResult<Vec<Py<BatchError>>> {
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
error_projection!(retained BatchImageError,BatchImageErrorKind,ck::BatchImageError,"batch",{
InvalidFormat(_),Write(_),
},cause=|inner:&RetainedError<ck::BatchImageError>|inner.node().source().map(crate::drawing_errors::DrawingWriteError::from_ref),
source=Option<Py<crate::drawing_errors::DrawingWriteError>>,methods={
#[getter]fn format(&self)->Option<&str>{match self.inner.as_ref(){ck::BatchImageError::InvalidFormat(v)=>Some(v),_=>None}}
});
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    pub struct BatchExportReport {
        pub(crate) inner: ck::BatchExportReport,
    }
}
#[pymethods]
impl BatchExportReport {
    #[getter]
    fn total(&self) -> usize {
        self.inner.total()
    }
    #[getter]
    fn success(&self) -> usize {
        self.inner.success()
    }
    #[getter]
    fn failed(&self) -> usize {
        self.inner.failed()
    }
    #[getter]
    fn errors(&self, py: Python<'_>) -> PyResult<Vec<Py<BatchError>>> {
        self.inner
            .errors()
            .iter()
            .map(|v| Ok(BatchError::from_rust(v.clone()).value(py).extract()?))
            .collect()
    }
    fn write_report(&self, py: Python<'_>, path: PathBuf) -> PyResult<()> {
        py.detach(|| self.inner.write_report(&path))
            .map_err(BatchValidationError::from_rust)
    }
}
