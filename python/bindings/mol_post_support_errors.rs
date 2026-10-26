//! Exact attachment and bond-direction dependencies of native IO finalization.
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::HydrogenError;
use crate::model_edit_errors::{AtomPropertyError, QueryStateError};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::stereo_errors::{AtropSource, AtropisomerError};
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type AttachmentSource = Source<
    crate::error_framework::ErrorObject<TopologyValidationError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        Source<
            crate::error_framework::ErrorObject<QueryStateError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                Source<
                    crate::error_framework::ErrorObject<TopologyEditError>,
                    Source<
                        crate::error_framework::ErrorObject<ValenceError>,
                        crate::error_framework::ErrorObject<HydrogenError>,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(AttachmentExpansionError, AttachmentExpansionErrorKind, ck::AttachmentExpansionError, "io", {
    Topology(_), Coordinates(_), ValueCount { actual: usize, expected: usize },
    Query(_), Property(_), Edit(_), Valence(_), Placement(_),
}, cause = |inner: &ck::AttachmentExpansionError| match inner {
    ck::AttachmentExpansionError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Query(value) => Some(QueryStateError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Property(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Edit(value) => Some(TopologyEditError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::AttachmentExpansionError::Placement(value) => Some(HydrogenError::from_rust(value.clone())),
    ck::AttachmentExpansionError::ValueCount { .. } => None,
}, source = Option<AttachmentSource>);

error_projection!(BondDirectionStereoError, BondDirectionStereoErrorKind, ck::BondDirectionStereoError, "stereo", {
    InvalidState(_), Valence(_), Atropisomer(_),
}, cause = |inner: &ck::BondDirectionStereoError| Python::attach(|py| match inner {
    ck::BondDirectionStereoError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::BondDirectionStereoError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).cause(py),
    ck::BondDirectionStereoError::InvalidState(_) => None,
}), source = Option<AtropSource>, methods = {
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to BondDirectionStereoError.InvalidState. None when this payload is absent."] fn message(&self) -> Option<&str> { match &self.inner {
        ck::BondDirectionStereoError::InvalidState(value) => Some(value), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of BondDirectionStereoError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> {
        match &self.inner {
            ck::BondDirectionStereoError::Valence(value) => ValenceError::from_rust(value.clone())
                .value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured AtropisomerError payload of BondDirectionStereoError.Atropisomer; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn atropisomer(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtropisomerError>>> {
        match &self.inner {
            ck::BondDirectionStereoError::Atropisomer(value) => AtropisomerError::from_rust(value.clone())
                .value(py).extract().map(Some).map_err(Into::into),
            _ => Ok(None),
        }
    }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_attachment_and_direction_errors_keep_counts_wrapped_values_and_true_causes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::AttachmentExpansionError::ValueCount {
                actual: usize::MAX,
                expected: 0,
            };
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = AttachmentExpansionError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("actual")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("actual", 0).is_err());

            let native = ck::AttachmentExpansionError::Placement(
                ck::HydrogenError::InvalidTopology(ck::TopologyValidationError::AdjacencyMismatch),
            );
            assert!(native.source().unwrap().is::<ck::HydrogenError>());
            let error = AttachmentExpansionError::from_rust(native);
            let source = error.value(py).getattr("source")?;
            assert!(source.is_instance(<HydrogenError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            // Hydrogen's InvalidTopology is a payload, not a native source.
            assert!(source.getattr("source")?.is_none());

            for (native, field) in [
                (
                    ck::BondDirectionStereoError::Valence(ck::ValenceError::InvalidTopology {
                        source: ck::TopologyValidationError::AdjacencyMismatch,
                    }),
                    "valence",
                ),
                (
                    ck::BondDirectionStereoError::Atropisomer(
                        ck::AtropisomerError::InvalidTopology {
                            source: ck::TopologyValidationError::AdjacencyMismatch,
                        },
                    ),
                    "atropisomer",
                ),
            ] {
                assert!(native.source().unwrap().is::<ck::TopologyValidationError>());
                let display = native.to_string();
                let error = BondDirectionStereoError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let source = error.value(py).getattr("source")?;
                assert!(source.is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                let payload = error.value(py).getattr(field)?;
                assert!(!payload.is_none());
                assert!(
                    payload
                        .getattr("source")?
                        .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
                );
                assert!(error.value(py).setattr(field, py.None()).is_err());
            }
            let native = ck::BondDirectionStereoError::InvalidState("state\0message".into());
            assert!(native.source().is_none());
            let error = BondDirectionStereoError::from_rust(native);
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                "state\0message"
            );
            assert!(error.value(py).getattr("source")?.is_none());
            assert!(error.value(py).setattr("message", "other").is_err());
            Ok(())
        })
        .unwrap();
    }
}
