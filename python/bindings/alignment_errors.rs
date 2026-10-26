//! Exact error closure for the complete read-only alignment query family.
use crate::error_framework::{Source, error_projection};
use crate::io_source_errors::{ErrorRef, IoError};
use crate::match_errors::{MatchSource, SubstructMatchError};
use crate::model_errors::CoordinateValidationError;
use crate::query_model_errors::QueryGraphError;
use crate::smarts_errors::{SmartsParseError, SmartsSource};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;
use std::sync::Arc;

crate::enums::binding_enum!(ThreadCountErrorKind, {
    UndefinedSignedNegation,
    ObservationIo,
    InvalidCpuList,
    AutomaticObservationUnavailable,
});
crate::error_framework::error_projection! {record
pub struct ThreadCountError {
    node: ErrorRef,
}
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
                    let error = crate::error_framework::from_payload(value);
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
crate::error_framework::error_projection! {methods impl ThreadCountError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "threads"
    }
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
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
    #[doc = "CPU-list segment containing the invalid value. Present for `ThreadCountError.InvalidCpuList`; None when unavailable."] fn segment(&self) -> Option<usize> {
        match self.native() {
            ck::ThreadCountError::InvalidCpuList { segment, .. } => Some(*segment),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Native explanation for rejecting the input. Present for `ThreadCountError.InvalidCpuList`; None when unavailable."] fn reason(&self) -> Option<&'static str> {
        match self.native() {
            ck::ThreadCountError::InvalidCpuList { reason, .. } => Some(reason),
            _ => None,
        }
    }
    #[getter]
    #[doc = "Platform associated with native thread selection. Present for `ThreadCountError.AutomaticObservationUnavailable`; None when unavailable."] fn platform(&self) -> Option<&'static str> {
        match self.native() {
            ck::ThreadCountError::AutomaticObservationUnavailable { platform } => Some(platform),
            _ => None,
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<Option<crate::error_framework::ErrorObject<IoError>>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.native().to_string()
    }
}
}

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct ThreadSelectionError {
    inner: Arc<ck::ThreadSelectionError>,
}
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
                    let error = crate::error_framework::from_payload(value);
                    error.set_cause(py, Some(cause));
                    error
                }
                Err(error) => error,
            }
        })
    }
}
crate::error_framework::error_projection! {methods impl ThreadSelectionError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "alignment"
    }
    #[staticmethod]
    #[doc = r#"{"role":"error_source"}"#]
    /// Immediate native cause, also available as __cause__; None when the native error has no cause.
    fn source(slf: Bound<'_, pyo3::exceptions::PyBaseException>) -> PyResult<crate::error_framework::ErrorObject<ThreadCountError>> {
        Ok(slf.getattr("__cause__")?.extract()?)
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

type AlignmentSource = Source<
    crate::error_framework::ErrorObject<ThreadCountError>,
    Source<MatchSource, SmartsSource>,
>;
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
    #[doc = "Index reported by the failing native check. Applies to AlignmentError.ProbeAtomOutOfRange, AlignmentError.ReferenceAtomOutOfRange, AlignmentError.NonPositiveWeight. None when this payload is absent."] fn index(&self) -> Option<usize> {
        match self.inner { ck::AlignmentError::ProbeAtomOutOfRange { index, .. }
        | ck::AlignmentError::ReferenceAtomOutOfRange { index, .. } | ck::AlignmentError::NonPositiveWeight { index } => Some(index), _ => None }
    }
    #[getter]
    #[doc = "Atom count used when validating the reported atom reference. Applies to AlignmentError.ProbeAtomOutOfRange, AlignmentError.ReferenceAtomOutOfRange. None when this payload is absent."] fn atom_count(&self) -> Option<usize> {
        match self.inner { ck::AlignmentError::ProbeAtomOutOfRange { atom_count, .. }
        | ck::AlignmentError::ReferenceAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None }
    }
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to AlignmentError.TerminalGroupSymmetrization, AlignmentError.NumericalPrecondition. None when this payload is absent."] fn message(&self) -> Option<&'static str> {
        match self.inner { ck::AlignmentError::TerminalGroupSymmetrization { message }
        | ck::AlignmentError::NumericalPrecondition { message } => Some(message), _ => None }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured CoordinateValidationError payload of AlignmentError.InvalidCoordinates; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn coordinates(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CoordinateValidationError>>> {
        match &self.inner { ck::AlignmentError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured SubstructMatchError payload of AlignmentError.Matching; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn matching(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<SubstructMatchError>>> {
        match &self.inner { ck::AlignmentError::Matching(value) => SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured QueryGraphError payload of AlignmentError.QueryGraph; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_graph(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryGraphError>>> {
        match &self.inner { ck::AlignmentError::QueryGraph(value) => QueryGraphError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured SmartsParseError payload of AlignmentError.QueryParse; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_parse(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<SmartsParseError>>> {
        match &self.inner { ck::AlignmentError::QueryParse(value) => SmartsParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None) }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ThreadSelectionError payload of AlignmentError.ThreadSelection; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn thread_selection(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ThreadSelectionError>>> {
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
                    .is_instance(<CoordinateValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
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
                    .is_instance(<TemplateAttachmentOrderError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
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
                    .is_instance(<QueryGraphError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
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
            assert!(retained.bind(py).is_instance(<IoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
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
