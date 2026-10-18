//! The actual transparent SMILES stereo-stage wrappers and their native values.
use crate::double_stereo_errors::DoubleBondStereoError;
use crate::error_framework::error_projection;
use crate::model_edit_errors::MoleculePropertyError;
use crate::model_errors::CoordinateValidationError;
use crate::ring_errors::RingFindingError;
use crate::stereo_errors::{LegacySource, LegacyStereoError};
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(SmilesStereoError, SmilesStereoErrorKind, ck::SmilesStereoError, "smiles_stereo", {
    Properties(_), SelectedConformerMissing { id: usize },
    ValenceRows { field: &'static str, actual, expected },
    RingRows { dimension: &'static str, actual, expected },
    Coordinates(_), Directions(_), Assignment(_), Rings(_), Valence(_),
}, cause = |inner: &ck::SmilesStereoError| Python::attach(|py| match inner {
    ck::SmilesStereoError::Properties(value) => MoleculePropertyError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::Coordinates(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::Directions(value) => DoubleBondStereoError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::Assignment(value) => LegacyStereoError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::Rings(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::SmilesStereoError::SelectedConformerMissing { .. } | ck::SmilesStereoError::ValenceRows { .. }
    | ck::SmilesStereoError::RingRows { .. } => None,
}), source = Option<LegacySource>, methods = {
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::SmilesStereoError::ValenceRows { actual, .. }
        | ck::SmilesStereoError::RingRows { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::SmilesStereoError::ValenceRows { expected, .. }
        | ck::SmilesStereoError::RingRows { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn properties(&self, py: Python<'_>) -> PyResult<Option<Py<MoleculePropertyError>>> { match &self.inner {
        ck::SmilesStereoError::Properties(value) => MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinates(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> { match &self.inner {
        ck::SmilesStereoError::Coordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn directions(&self, py: Python<'_>) -> PyResult<Option<Py<DoubleBondStereoError>>> { match &self.inner {
        ck::SmilesStereoError::Directions(value) => DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn assignment(&self, py: Python<'_>) -> PyResult<Option<Py<LegacyStereoError>>> { match &self.inner {
        ck::SmilesStereoError::Assignment(value) => LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn rings(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> { match &self.inner {
        ck::SmilesStereoError::Rings(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::SmilesStereoError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smiles_stereo_transparency_keeps_value_and_actual_deeper_cause_distinct() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::SmilesStereoError::Directions(ck::DoubleBondStereoError::NumericPropertyRead {
                    atom: Some(ck::AtomId::new(2)),
                    bond: None,
                    property: "_UnknownStereo",
                    source: ck::PropertyIntReadError::UnsignedOverflow { value: u32::MAX },
                });
            assert!(native.source().unwrap().is::<ck::PropertyIntReadError>());
            let message = native.to_string();
            let error = SmilesStereoError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(cause.getattr("value")?.extract::<u32>()?, u32::MAX);
            let payload = error.value(py).getattr("directions")?;
            assert!(payload.is_instance_of::<DoubleBondStereoError>());
            assert_eq!(payload.getattr("atom")?.extract::<usize>()?, 2);
            assert!(payload.getattr("bond")?.is_none());
            assert!(payload.setattr("atom", 3).is_err());
            Ok(())
        })
        .unwrap();
    }
}
