//! Exact XYZ reader failures, including the actual source-free standard leaf.
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_leaf"}"#]
pub struct ParseFloatError {
    inner: std::num::ParseFloatError,
}
}
crate::error_framework::error_projection! {methods impl ParseFloatError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "std"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl ParseFloatError {
    fn from_rust(inner: std::num::ParseFloatError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

error_projection!(XyzReadError, XyzReadErrorKind, ck::XyzReadError, "io", {
    EmptyBlock, AtomCount { value }, UnexpectedEof,
    MissingCoordinates { line }, Coordinate { value, line, source },
    AtomSymbol { message: String }, Topology(_), Coordinates(_), MoleculeProperty(_),
}, cause = |inner: &ck::XyzReadError| match inner {
    ck::XyzReadError::Coordinate { source, .. } => Some(ParseFloatError::from_rust(source.clone())),
    ck::XyzReadError::Topology(source) => Some(TopologyValidationError::from_rust(source.clone())),
    ck::XyzReadError::Coordinates(source) => Some(CoordinateValidationError::from_rust(source.clone())),
    ck::XyzReadError::MoleculeProperty(source) => Some(MoleculePropertyError::from_rust(source.clone())),
    ck::XyzReadError::EmptyBlock | ck::XyzReadError::AtomCount { .. } | ck::XyzReadError::UnexpectedEof
    | ck::XyzReadError::MissingCoordinates { .. } | ck::XyzReadError::AtomSymbol { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<ParseFloatError>, Source<crate::error_framework::ErrorObject<TopologyValidationError>, Source<crate::error_framework::ErrorObject<CoordinateValidationError>, crate::error_framework::ErrorObject<MoleculePropertyError>>>>>, methods = {
    #[getter]
    #[doc = "Original value rejected by the operation. Present for `XyzReadError.AtomCount`, `XyzReadError.Coordinate`; None when unavailable."] fn value(&self) -> Option<&str> { match &self.inner {
        ck::XyzReadError::AtomCount { value } | ck::XyzReadError::Coordinate { value, .. } => Some(value), _ => None,
    } }
    #[getter]
    #[doc = "Source line associated with the failure. Present for `XyzReadError.MissingCoordinates`, `XyzReadError.Coordinate`; None when unavailable."] fn line(&self) -> Option<usize> { match self.inner {
        ck::XyzReadError::MissingCoordinates { line } | ck::XyzReadError::Coordinate { line, .. } => Some(line), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_xyz_coordinate_failure_keeps_concrete_standard_leaf() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let root = ck::Molecule::from_xyz_block("1\nname\nC Q 0 0\n").unwrap_err();
            let ck::MolecularIoError::XyzRead(native) = root else {
                panic!("expected native XYZ read error")
            };
            let source = native.source().unwrap();
            assert!(source.is::<std::num::ParseFloatError>());
            assert!(source.source().is_none());
            let source_display = source.to_string();
            let display = native.to_string();
            let error = XyzReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(error.value(py).getattr("value")?.extract::<String>()?, "Q");
            let source = error.cause(py).unwrap();
            assert!(
                source.is_instance(
                    py,
                    <ParseFloatError as crate::error_framework::ExceptionPayload>::exception_type(
                        py
                    )
                    .expect("error backend")
                    .as_any()
                )
            );
            assert!(error.value(py).getattr("source")?.is(source.value(py)));
            assert_eq!(source.value(py).str()?.to_str()?, source_display);
            assert_eq!(
                source.value(py).getattr("domain")?.extract::<String>()?,
                "std"
            );
            assert!(source.cause(py).is_none());
            assert!(!source.value(py).hasattr("kind")?);
            assert!(
                source
                    .value(py)
                    .setattr("domain", "other")
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            Ok(())
        })
        .unwrap();
    }
}
