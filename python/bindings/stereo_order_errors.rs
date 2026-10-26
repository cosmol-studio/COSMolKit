//! Tetrahedral ordering failures with native causes and detached typed payloads.
use crate::chemical_values::ChiralTag;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{AtomPropertyError, MappingValidationError};
use crate::model_errors::TopologyValidationError;
use crate::property_numeric_errors::PropertyUIntReadError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(StereoOrderError, StereoOrderErrorKind, ck::StereoOrderError, "stereo_order", {
    ChiralPermutationRead(_), ChiralPermutationWrite(_),
    ConflictingChiralPermutation { atom: usize => |value: &ck::AtomId| value.index(), typed: u32, raw: u32 },
    PermutationLength { reference: usize, probe: usize }, MissingProbeValue { reference_position: usize },
    BondIndexSourceWidth { bond_index: usize }, InvalidTopology(_), InvalidMapping(_),
    CenterOutOfRange { center, atom_count: usize }, BondOutOfRange { bond, bond_count: usize },
    BondNotIncident { bond, center }, TooManyLigands { actual: usize, maximum: usize },
    DuplicateBondLigand { position: usize, bond }, MultipleImplicitLigands { first: usize, second: usize },
    UnsupportedTetrahedralTag { tag: ChiralTag => |value: &ck::ChiralTag| (*value).into() },
    RemovedCenter { center }, RemovedLigand { bond },
    ImplicitReplacementOutOfRange { bond, old_bond_count: usize },
    ImplicitReplacementRetained { bond, mapped: usize => |value: &ck::BondId| value.index() },
    ImplicitReplacementNotInSource { bond },
}, cause = |inner: &ck::StereoOrderError| match inner {
    ck::StereoOrderError::ChiralPermutationRead(value) => Some(PropertyUIntReadError::from_rust(value.clone())),
    ck::StereoOrderError::ChiralPermutationWrite(value) => Some(AtomPropertyError::from_rust(value.clone())),
    // The native topology/mapping tuple fields have neither #[source] nor
    // #[from]. Preserve them as payloads without inventing causal edges.
    ck::StereoOrderError::InvalidTopology(_) | ck::StereoOrderError::InvalidMapping(_)
    | ck::StereoOrderError::ConflictingChiralPermutation { .. } | ck::StereoOrderError::PermutationLength { .. }
    | ck::StereoOrderError::MissingProbeValue { .. } | ck::StereoOrderError::BondIndexSourceWidth { .. }
    | ck::StereoOrderError::CenterOutOfRange { .. } | ck::StereoOrderError::BondOutOfRange { .. }
    | ck::StereoOrderError::BondNotIncident { .. } | ck::StereoOrderError::TooManyLigands { .. }
    | ck::StereoOrderError::DuplicateBondLigand { .. } | ck::StereoOrderError::MultipleImplicitLigands { .. }
    | ck::StereoOrderError::UnsupportedTetrahedralTag { .. } | ck::StereoOrderError::RemovedCenter { .. }
    | ck::StereoOrderError::RemovedLigand { .. } | ck::StereoOrderError::ImplicitReplacementOutOfRange { .. }
    | ck::StereoOrderError::ImplicitReplacementRetained { .. } | ck::StereoOrderError::ImplicitReplacementNotInSource { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<PropertyUIntReadError>, crate::error_framework::ErrorObject<AtomPropertyError>>>, methods = {
    #[getter]
    #[doc = "Stereochemical center reported by the native check. Applies to StereoOrderError.CenterOutOfRange, StereoOrderError.BondNotIncident, StereoOrderError.RemovedCenter. None when this payload is absent."] fn center(&self) -> Option<usize> { match self.inner {
        ck::StereoOrderError::CenterOutOfRange { center, .. }
        | ck::StereoOrderError::BondNotIncident { center, .. }
        | ck::StereoOrderError::RemovedCenter { center } => Some(center.index()), _ => None,
    } }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to StereoOrderError.BondOutOfRange, StereoOrderError.BondNotIncident, StereoOrderError.DuplicateBondLigand, StereoOrderError.RemovedLigand, StereoOrderError.ImplicitReplacementOutOfRange, StereoOrderError.ImplicitReplacementRetained, StereoOrderError.ImplicitReplacementNotInSource. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::StereoOrderError::BondOutOfRange { bond, .. }
        | ck::StereoOrderError::BondNotIncident { bond, .. }
        | ck::StereoOrderError::DuplicateBondLigand { bond, .. }
        | ck::StereoOrderError::RemovedLigand { bond }
        | ck::StereoOrderError::ImplicitReplacementOutOfRange { bond, .. }
        | ck::StereoOrderError::ImplicitReplacementRetained { bond, .. }
        | ck::StereoOrderError::ImplicitReplacementNotInSource { bond } => Some(bond.index()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of StereoOrderError.InvalidTopology; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::StereoOrderError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured MappingValidationError payload of StereoOrderError.InvalidMapping; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn mapping(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MappingValidationError>>> { match &self.inner {
        ck::StereoOrderError::InvalidMapping(value) => MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;
    use crate::property_numeric_errors::PropertyUIntReadErrorKind;

    #[test]
    fn stereo_order_projection_preserves_payloads_and_actual_cause_identity() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::StereoOrderError::ChiralPermutationRead(ck::PropertyUIntReadError::Negative {
                    value: -7,
                });
            assert!(native.source().is_some());
            let display = native.to_string();
            let error = StereoOrderError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<StereoOrderErrorKind>()?
                    == StereoOrderErrorKind::ChiralPermutationRead
            );
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                source
                    .getattr("variant")?
                    .extract::<PropertyUIntReadErrorKind>()?
                    == PropertyUIntReadErrorKind::Negative
            );
            assert_eq!(source.getattr("value")?.extract::<i64>()?, -7);
            assert!(source.getattr("__cause__")?.is_none());

            let native = ck::StereoOrderError::ImplicitReplacementRetained {
                bond: ck::BondId::new(2),
                mapped: ck::BondId::new(4),
            };
            assert!(native.source().is_none());
            let error = StereoOrderError::from_rust(native);
            assert_eq!(error.value(py).getattr("bond")?.extract::<usize>()?, 2);
            assert_eq!(error.value(py).getattr("mapped")?.extract::<usize>()?, 4);
            assert!(error.value(py).setattr("mapped", 8).is_err());
            assert!(error.cause(py).is_none());
            let native = ck::StereoOrderError::InvalidTopology(
                ck::TopologyValidationError::AdjacencyMismatch,
            );
            assert!(native.source().is_none());
            let error = StereoOrderError::from_rust(native);
            assert!(error.cause(py).is_none());
            assert!(
                error
                    .value(py)
                    .getattr("topology")?
                    .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            Ok(())
        })
        .unwrap();
    }
}
