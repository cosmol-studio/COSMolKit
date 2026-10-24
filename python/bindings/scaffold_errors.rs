//! Native scaffold failures without changing transparent source boundaries.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::matrix_errors::MatrixError;
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::TopologyValidationError;
use crate::property_numeric_errors::{PropertyUIntReadError, UIntLexicalReadError};
use crate::ring_errors::RingFindingError;
use crate::sanitize_leaf_errors::RadicalError;
use crate::stereo_errors::{LegacySource, LegacyStereoError};
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(ScaffoldError, ScaffoldErrorKind, ck::ScaffoldError, "scaffold", {
    Topology(_), Edit(_), Valence(_), Rings(_), Radicals(_), Stereo(_),
    Property(_), AtomProperty(_), BondProperty(_), UIntProperty(_), TextProperty(_), Matrix(_),
    MissingImplicitHydrogens { atom: usize => |value: &ck::AtomId| value.index() },
    ValenceRowCount { field: &'static str, actual: usize, expected: usize },
    MissingRingInfo, MatrixSize { atoms: usize },
}, cause = |inner: &ck::ScaffoldError| Python::attach(|py| match inner {
    ck::ScaffoldError::Topology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Edit(value) => TopologyEditError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Rings(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Radicals(value) => RadicalError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Stereo(value) => LegacyStereoError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Property(value) => MoleculePropertyError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::BondProperty(value) => BondValueError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::UIntProperty(value) => PropertyUIntReadError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::TextProperty(value) => PropertyStringError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::Matrix(value) => MatrixError::from_rust(value.clone()).cause(py),
    ck::ScaffoldError::MissingImplicitHydrogens { .. } | ck::ScaffoldError::ValenceRowCount { .. }
    | ck::ScaffoldError::MissingRingInfo | ck::ScaffoldError::MatrixSize { .. } => None,
}), source = Option<Source<LegacySource, crate::error_framework::ErrorObject<UIntLexicalReadError>>>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::ScaffoldError::Topology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn edit(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyEditError>>> { match &self.inner {
        ck::ScaffoldError::Edit(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::ScaffoldError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn rings(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>> { match &self.inner {
        ck::ScaffoldError::Rings(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn radicals(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RadicalError>>> { match &self.inner {
        ck::ScaffoldError::Radicals(value) => RadicalError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn stereo(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<LegacyStereoError>>> { match &self.inner {
        ck::ScaffoldError::Stereo(value) => LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn property_error(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MoleculePropertyError>>> { match &self.inner {
        ck::ScaffoldError::Property(value) => MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atom_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtomPropertyError>>> { match &self.inner {
        ck::ScaffoldError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn bond_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<BondValueError>>> { match &self.inner {
        ck::ScaffoldError::BondProperty(value) => BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn uint_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyUIntReadError>>> { match &self.inner {
        ck::ScaffoldError::UIntProperty(value) => PropertyUIntReadError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn text_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> { match &self.inner {
        ck::ScaffoldError::TextProperty(value) => PropertyStringError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn matrix(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MatrixError>>> { match &self.inner {
        ck::ScaffoldError::Matrix(value) => MatrixError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scaffold_errors_keep_transparent_causes_full_payloads_and_native_width() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::ScaffoldError::UIntProperty(ck::PropertyUIntReadError::Lexical {
                value: "invalid".into(),
                source: ck::UIntLexicalReadError::Character {
                    position: 3,
                    byte: b'x',
                },
            });
            assert!(native.source().unwrap().is::<ck::UIntLexicalReadError>());
            let message = native.to_string();
            let error = ScaffoldError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance(<UIntLexicalReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(cause.getattr("position")?.extract::<usize>()?, 3);
            assert!(
                error
                    .value(py)
                    .getattr("uint_property")?
                    .is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            let native = ck::ScaffoldError::Edit(ck::TopologyEditError::InvalidSource(
                ck::TopologyValidationError::AdjacencyMismatch,
            ));
            assert!(native.source().is_none());
            let error = ScaffoldError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("edit")?
                    .is_instance(<TopologyEditError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            let error = ScaffoldError::from_rust(ck::ScaffoldError::ValenceRowCount {
                field: "implicit_hydrogens",
                actual: usize::MAX,
                expected: 3,
            });
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<ScaffoldErrorKind>()?
                    == ScaffoldErrorKind::ValenceRowCount
            );
            assert_eq!(
                error.value(py).getattr("field")?.extract::<String>()?,
                "implicit_hydrogens"
            );
            assert_eq!(
                error.value(py).getattr("actual")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(error.value(py).getattr("expected")?.extract::<usize>()?, 3);
            assert!(error.value(py).setattr("actual", 0).is_err());
            Ok(())
        })
        .unwrap();
    }
}
