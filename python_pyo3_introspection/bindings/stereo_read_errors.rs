//! Actual readonly stereo-query failures and their concrete native causes.
use crate::canonical_errors::PropertyStringError;
use crate::cip_atom_code_errors::CipLabelerError;
use crate::error_framework::{Source, error_projection};
use crate::model_errors::TopologyValidationError;
use crate::potential_stereo_errors::PotentialStereoError;
use crate::ring_errors::RingFindingError;
use crate::stereo_errors::LegacySource;
use crate::utf8_errors::FromUtf8Error;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type StereoReadSource = Source<
    LegacySource,
    Source<
        Py<PotentialStereoError>,
        Source<Py<CipLabelerError>, Source<Py<PropertyStringError>, Py<FromUtf8Error>>>,
    >,
>;

error_projection!(StereoReadError, StereoReadErrorKind, ck::StereoReadError, "stereo", {
    InvalidTopology(_), Valence(_), Rings(_), PotentialStereo(_), CipLabeler(_),
    PropertyString(_), CipLabelEncoding { atom: usize => |value: &ck::AtomId| value.index(), source },
}, cause = |inner: &ck::StereoReadError| Python::attach(|py| match inner {
    ck::StereoReadError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::StereoReadError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::StereoReadError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::StereoReadError::PotentialStereo(value) => Some(PotentialStereoError::from_rust(value.clone())),
    ck::StereoReadError::CipLabeler(value) => Some(CipLabelerError::from_rust(value.clone())),
    ck::StereoReadError::PropertyString(value) => Some(PropertyStringError::from_rust(*value)),
    ck::StereoReadError::CipLabelEncoding { source, .. } => Some(FromUtf8Error::from_rust(source.clone())),
}), source = Option<StereoReadSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::StereoReadError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::property_values::PropertyValueKind;

    #[test]
    fn stereo_read_errors_preserve_native_cause_depth_encoding_and_topology_payload() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::StereoReadError::CipLabeler(ck::CipLabelerError::InvalidPropertyKind {
                    atom: 3,
                    property: "_CIPCode",
                    kind: ck::PropertyValueKind::Bool,
                });
            assert!(native.source().unwrap().is::<ck::CipLabelerError>());
            let message = native.to_string();
            let error = StereoReadError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance_of::<CipLabelerError>());
            assert_eq!(cause.getattr("atom")?.extract::<usize>()?, 3);
            assert!(
                cause.getattr("kind")?.extract::<PropertyValueKind>()? == PropertyValueKind::Bool
            );

            let bytes = vec![0, 255];
            let encoding = String::from_utf8(bytes.clone()).unwrap_err();
            let native = ck::StereoReadError::CipLabelEncoding {
                atom: ck::AtomId::new(usize::MAX),
                source: encoding,
            };
            assert!(native.source().unwrap().is::<std::string::FromUtf8Error>());
            let error = StereoReadError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<StereoReadErrorKind>()?
                    == StereoReadErrorKind::CipLabelEncoding
            );
            assert_eq!(
                error.value(py).getattr("atom")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).setattr("atom", 0).is_err());
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<FromUtf8Error>());
            assert_eq!(cause.call_method0("as_bytes")?.extract::<Vec<u8>>()?, bytes);
            assert_eq!(
                cause
                    .getattr("utf8_error")?
                    .getattr("valid_up_to")?
                    .extract::<usize>()?,
                1
            );

            let native = ck::StereoReadError::InvalidTopology(
                ck::TopologyValidationError::AdjacencyMismatch,
            );
            assert!(native.source().is_none());
            let error = StereoReadError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(
                error
                    .value(py)
                    .getattr("topology")?
                    .is_instance_of::<TopologyValidationError>()
            );
            Ok(())
        })
        .unwrap();
    }
}
