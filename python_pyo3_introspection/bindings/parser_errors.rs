//! Concrete shared parser-carrier errors, preserving the native cause levels.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::AtomPropertyError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(ParserCarrierError, ParserCarrierErrorKind, ck::ParserCarrierError, "parser_carrier", {
    PropertyString(_), AtomProperty(_), Model(_),
}, cause = |inner: &ck::ParserCarrierError| match inner {
    ck::ParserCarrierError::PropertyString(value) => Some(PropertyStringError::from_rust(*value)),
    ck::ParserCarrierError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::ParserCarrierError::Model(_) => None,
}, source = Option<Source<Py<PropertyStringError>, Py<AtomPropertyError>>>, methods = {
    #[getter]
    fn model(&self) -> Option<&str> { match &self.inner {
        ck::ParserCarrierError::Model(value) => Some(value), _ => None,
    } }
});

error_projection!(ParserStereoOrderError, ParserStereoOrderErrorKind, ck::ParserStereoOrderError, "parser_stereo_order", {
    BondIndexOutOfRange { bond_index: usize, bond_count: usize },
    ProbeLengthMismatch { incident_count: usize, probe_count: usize },
    ProbeMemberMissing { probe_index: usize },
}, cause = |_inner: &ck::ParserStereoOrderError| None::<PyErr>);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::property_values::PropertyValueKind;

    #[test]
    fn parser_errors_retain_typed_nontransparent_cause_and_native_row_payloads() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::ParserCarrierError::PropertyString(ck::PropertyStringError::UnsupportedKind {
                    kind: ck::PropertyValueKind::Bool,
                });
            assert!(native.source().unwrap().is::<ck::PropertyStringError>());
            let message = native.to_string();
            let error = ParserCarrierError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance_of::<PropertyStringError>());
            assert!(
                cause.call_method0("kind")?.extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(cause.getattr("__cause__")?.is_none());

            let native = ck::ParserStereoOrderError::ProbeLengthMismatch {
                incident_count: 3,
                probe_count: 2,
            };
            assert!(native.source().is_none());
            let message = native.to_string();
            let error = ParserStereoOrderError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(
                error
                    .value(py)
                    .getattr("incident_count")?
                    .extract::<usize>()?,
                3
            );
            assert_eq!(
                error.value(py).getattr("probe_count")?.extract::<usize>()?,
                2
            );
            assert!(error.value(py).setattr("probe_count", 3).is_err());
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
