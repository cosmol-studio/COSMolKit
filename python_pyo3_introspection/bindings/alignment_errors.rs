//! Exact error closure for the complete read-only alignment query family.
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::{ErrorRef, IoError};
use crate::match_errors::{MatchSource, SubstructMatchError};
use crate::model_errors::CoordinateValidationError;
use crate::query_model_errors::QueryGraphError;
use crate::smarts_errors::{SmartsParseError, SmartsSource};
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
use std::sync::Arc;

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum ThreadCountErrorKind {
    UndefinedSignedNegation,
    ObservationIo,
    InvalidCpuList,
    AutomaticObservationUnavailable,
}
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct ThreadCountError {
    node: ErrorRef,
}
impl ThreadCountError {
    fn native(&self) -> &ck::ThreadCountError {
        self.node
            .error()
            .downcast_ref()
            .expect("native thread count error")
    }
    pub(crate) fn from_ref(node: ErrorRef) -> PyErr {
        Python::attach(|py| {
            let cause = node.source().map(IoError::from_ref);
            match Py::new(py, Self { node }) {
                Ok(value) => {
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, cause);
                    error
                }
                Err(error) => error,
            }
        })
    }
    #[cfg(test)]
    fn from_rust(inner: ck::ThreadCountError) -> PyErr {
        Self::from_ref(ErrorRef::new(Arc::new(inner)))
    }
}
#[pymethods]
impl ThreadCountError {
    #[getter]
    fn domain(&self) -> &'static str {
        "threads"
    }
    #[getter]
    fn variant(&self) -> ThreadCountErrorKind {
        match self.native() {
            ck::ThreadCountError::UndefinedSignedNegation => {
                ThreadCountErrorKind::UndefinedSignedNegation
            }
            ck::ThreadCountError::ObservationIo(_) => ThreadCountErrorKind::ObservationIo,
            ck::ThreadCountError::InvalidCpuList {
                segment: _,
                reason: _,
            } => ThreadCountErrorKind::InvalidCpuList,
            ck::ThreadCountError::AutomaticObservationUnavailable { platform: _ } => {
                ThreadCountErrorKind::AutomaticObservationUnavailable
            }
        }
    }
    #[getter]
    fn segment(&self) -> Option<usize> {
        match self.native() {
            ck::ThreadCountError::InvalidCpuList { segment, .. } => Some(*segment),
            _ => None,
        }
    }
    #[getter]
    fn reason(&self) -> Option<&'static str> {
        match self.native() {
            ck::ThreadCountError::InvalidCpuList { reason, .. } => Some(reason),
            _ => None,
        }
    }
    #[getter]
    fn platform(&self) -> Option<&'static str> {
        match self.native() {
            ck::ThreadCountError::AutomaticObservationUnavailable { platform } => Some(platform),
            _ => None,
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Option<Py<IoError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.native().to_string()
    }
}

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module="fixed_projection._native", extends=PyValueError)]
pub struct ThreadSelectionError {
    inner: Arc<ck::ThreadSelectionError>,
}
impl ThreadSelectionError {
    pub(crate) fn from_rust(inner: ck::ThreadSelectionError) -> PyErr {
        let inner = Arc::new(inner);
        Python::attach(|py| {
            let root: Arc<dyn Error + Send + Sync> = inner.clone();
            let cause = ThreadCountError::from_ref(
                ErrorRef::new(root)
                    .source()
                    .expect("native thread selection cause"),
            );
            match Py::new(py, Self { inner }) {
                Ok(value) => {
                    let error = PyErr::from_value(value.into_bound(py).into_any());
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}
#[pymethods]
impl ThreadSelectionError {
    #[getter]
    fn domain(&self) -> &'static str {
        "alignment"
    }
    #[getter]
    #[doc = r#"{"role":"error_source"}"#]
    fn source(slf: Bound<'_, Self>) -> PyResult<Py<ThreadCountError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

type AlignmentSource = Source<Py<ThreadCountError>, Source<MatchSource, SmartsSource>>;
error_projection!(AlignmentError, AlignmentErrorKind, ck::AlignmentError, "alignment", {
    InvalidCoordinates(_), Matching(_), QueryGraph(_), QueryParse(_), ThreadSelection(_),
    NoConformers, ConformerNotFound { id: i32 }, EmptyAtomMap,
    ProbeAtomOutOfRange { index, atom_count }, ReferenceAtomOutOfRange { index, atom_count },
    WeightCountMismatch { map_len: usize, weight_len: usize }, NonPositiveWeight { index },
    NoSubstructureMatch, TerminalGroupSymmetrization { message }, NumericalPrecondition { message },
    WorkerTerminated,
}, cause = |inner: &ck::AlignmentError| Python::attach(|py| match inner {
    // Every wrapper is native #[error(transparent)]: preserve the child's
    // source, retaining the wrapped error separately as a typed payload.
    ck::AlignmentError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::AlignmentError::Matching(value) => SubstructMatchError::from_rust(value.clone()).cause(py),
    ck::AlignmentError::QueryGraph(value) => QueryGraphError::from_rust(value.clone()).cause(py),
    ck::AlignmentError::QueryParse(value) => SmartsParseError::from_rust(value.clone()).cause(py),
    ck::AlignmentError::ThreadSelection(value) => ThreadSelectionError::from_rust(value.clone()).cause(py),
    ck::AlignmentError::NoConformers | ck::AlignmentError::ConformerNotFound { .. }
    | ck::AlignmentError::EmptyAtomMap | ck::AlignmentError::ProbeAtomOutOfRange { .. }
    | ck::AlignmentError::ReferenceAtomOutOfRange { .. } | ck::AlignmentError::WeightCountMismatch { .. }
    | ck::AlignmentError::NonPositiveWeight { .. } | ck::AlignmentError::NoSubstructureMatch
    | ck::AlignmentError::TerminalGroupSymmetrization { .. } | ck::AlignmentError::NumericalPrecondition { .. }
    | ck::AlignmentError::WorkerTerminated => None,
}), source = Option<AlignmentSource>, methods = {
    #[getter]
    fn index(&self) -> Option<usize> {
        match self.inner { ck::AlignmentError::ProbeAtomOutOfRange { index, .. }
        | ck::AlignmentError::ReferenceAtomOutOfRange { index, .. } | ck::AlignmentError::NonPositiveWeight { index } => Some(index), _ => None }
    }
    #[getter]
    fn atom_count(&self) -> Option<usize> {
        match self.inner { ck::AlignmentError::ProbeAtomOutOfRange { atom_count, .. }
        | ck::AlignmentError::ReferenceAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None }
    }
    #[getter]
    fn message(&self) -> Option<&'static str> {
        match self.inner { ck::AlignmentError::TerminalGroupSymmetrization { message }
        | ck::AlignmentError::NumericalPrecondition { message } => Some(message), _ => None }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinates(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> {
        match &self.inner { ck::AlignmentError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn matching(&self, py: Python<'_>) -> PyResult<Option<Py<SubstructMatchError>>> {
        match &self.inner { ck::AlignmentError::Matching(value) => SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_graph(&self, py: Python<'_>) -> PyResult<Option<Py<QueryGraphError>>> {
        match &self.inner { ck::AlignmentError::QueryGraph(value) => QueryGraphError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_parse(&self, py: Python<'_>) -> PyResult<Option<Py<SmartsParseError>>> {
        match &self.inner { ck::AlignmentError::QueryParse(value) => SmartsParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn thread_selection(&self, py: Python<'_>) -> PyResult<Option<Py<ThreadSelectionError>>> {
        match &self.inner { ck::AlignmentError::ThreadSelection(value) => ThreadSelectionError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::attachment_values::TemplateAttachmentOrderError;
    #[test]
    fn actual_alignment_error_payloads_and_transparent_sources_remain_typed_and_owned() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::AlignmentError::InvalidCoordinates(ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 0,
                    rows: 1,
                    atom_count: 2,
                });
            assert!(native.source().is_none());
            let error = AlignmentError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("coordinates")?
                    .is_instance_of::<CoordinateValidationError>()
            );
            let native =
                ck::AlignmentError::QueryGraph(ck::QueryGraphError::TemplateAttachmentOrder {
                    atom: ck::AtomId::new(1),
                    source: ck::TemplateAttachmentOrderError::Empty,
                });
            let error = AlignmentError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance_of::<TemplateAttachmentOrderError>()
            );
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is(&error.value(py).getattr("__cause__")?)
            );
            assert!(
                error
                    .value(py)
                    .getattr("query_graph")?
                    .is_instance_of::<QueryGraphError>()
            );
            let error = AlignmentError::from_rust(ck::AlignmentError::NumericalPrecondition {
                message: "native\0context",
            });
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                "native\0context"
            );
            assert!(error.value(py).setattr("message", "changed").is_err());
            for native in [
                ck::ThreadCountError::UndefinedSignedNegation,
                ck::ThreadCountError::InvalidCpuList {
                    segment: 3,
                    reason: "native\0reason",
                },
                ck::ThreadCountError::AutomaticObservationUnavailable {
                    platform: "native-platform",
                },
            ] {
                let display = native.to_string();
                let error = ThreadCountError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                assert!(error.value(py).getattr("source")?.is_none());
                assert!(error.value(py).setattr("segment", 9).is_err());
            }
            let error = ThreadCountError::from_rust(ck::ThreadCountError::ObservationIo(
                std::io::Error::from_raw_os_error(13),
            ));
            let retained = error.value(py).getattr("source")?.unbind();
            assert!(retained.bind(py).is_instance_of::<IoError>());
            assert!(retained.bind(py).is(&error.value(py).getattr("__cause__")?));
            drop(error);
            assert_eq!(
                retained
                    .bind(py)
                    .getattr("raw_os_error")?
                    .extract::<i32>()?,
                13
            );
            Ok(())
        })
        .unwrap();
    }
}
