//! Distinct native MCS candidate, search-state and progress error values.
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::AtomPropertyError;
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::stereo_order_errors::StereoOrderError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(SearchMcsError, SearchMcsErrorKind, ck::SearchMcsError, "search", {
    TooFewInputs { count }, ThresholdAboveOne,
    AtomOutOfRange { side, atom }, MissingConformer { side },
    CoordinateOutOfRange { side, atom }, MissingRingInfo { side },
    MissingValence { side }, ValenceOutOfRange { side, atom },
    BondOutOfRange { side, bond }, SeedExcludedBondOutOfRange { bond, count },
    SeedBondAlreadyExcluded { bond }, SeedAtomMappingMissing { atom },
    TargetAtomMappingMissing { atom }, SeedAtomNotBondEndpoint { atom, bond },
    SeedNewAtomMissing { atom },
    SeedRemainingBondUnderflow { remaining, frontier: usize },
    SeedRemainingAtomUnderflow { remaining, added: usize },
}, cause = |_: &ck::SearchMcsError| None::<PyErr>, methods = {
    #[getter]
    #[doc = "Count reported by the native failure. Applies to SearchMcsError.TooFewInputs, SearchMcsError.SeedExcludedBondOutOfRange. None when this payload is absent."] fn count(&self) -> Option<usize> { match self.inner {
        ck::SearchMcsError::TooFewInputs { count }
        | ck::SearchMcsError::SeedExcludedBondOutOfRange { count, .. } => Some(count), _ => None,
    } }
    #[getter]
    #[doc = "Probe/reference or source/target side identified by the native check. Applies to SearchMcsError.AtomOutOfRange, SearchMcsError.MissingConformer, SearchMcsError.CoordinateOutOfRange, SearchMcsError.MissingRingInfo, SearchMcsError.MissingValence, SearchMcsError.ValenceOutOfRange, SearchMcsError.BondOutOfRange. None when this payload is absent."] fn side(&self) -> Option<&'static str> { match self.inner {
        ck::SearchMcsError::AtomOutOfRange { side, .. } | ck::SearchMcsError::MissingConformer { side }
        | ck::SearchMcsError::CoordinateOutOfRange { side, .. } | ck::SearchMcsError::MissingRingInfo { side }
        | ck::SearchMcsError::MissingValence { side } | ck::SearchMcsError::ValenceOutOfRange { side, .. }
        | ck::SearchMcsError::BondOutOfRange { side, .. } => Some(side), _ => None,
    } }
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to SearchMcsError.AtomOutOfRange, SearchMcsError.CoordinateOutOfRange, SearchMcsError.ValenceOutOfRange, SearchMcsError.SeedAtomMappingMissing, SearchMcsError.TargetAtomMappingMissing, SearchMcsError.SeedAtomNotBondEndpoint, SearchMcsError.SeedNewAtomMissing. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::SearchMcsError::AtomOutOfRange { atom, .. } | ck::SearchMcsError::CoordinateOutOfRange { atom, .. }
        | ck::SearchMcsError::ValenceOutOfRange { atom, .. } | ck::SearchMcsError::SeedAtomMappingMissing { atom }
        | ck::SearchMcsError::TargetAtomMappingMissing { atom } | ck::SearchMcsError::SeedAtomNotBondEndpoint { atom, .. }
        | ck::SearchMcsError::SeedNewAtomMissing { atom } => Some(atom), _ => None,
    } }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to SearchMcsError.BondOutOfRange, SearchMcsError.SeedExcludedBondOutOfRange, SearchMcsError.SeedBondAlreadyExcluded, SearchMcsError.SeedAtomNotBondEndpoint. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::SearchMcsError::BondOutOfRange { bond, .. } | ck::SearchMcsError::SeedExcludedBondOutOfRange { bond, .. }
        | ck::SearchMcsError::SeedBondAlreadyExcluded { bond } | ck::SearchMcsError::SeedAtomNotBondEndpoint { bond, .. } => Some(bond), _ => None,
    } }
    #[getter]
    #[doc = "Remaining input or entries reported by the native failure. Applies to SearchMcsError.SeedRemainingBondUnderflow, SearchMcsError.SeedRemainingAtomUnderflow. None when this payload is absent."] fn remaining(&self) -> Option<usize> { match self.inner {
        ck::SearchMcsError::SeedRemainingBondUnderflow { remaining, .. }
        | ck::SearchMcsError::SeedRemainingAtomUnderflow { remaining, .. } => Some(remaining), _ => None,
    } }
});

error_projection!(McsProgressError, McsProgressErrorKind, ck::McsProgressError, "search", {
    Callback { message: String },
}, cause = |_: &ck::McsProgressError| None::<PyErr>);

error_projection!(McsError, McsErrorKind, ck::McsError, "search", {
    TargetTableCount { targets, tables: usize },
    ThresholdCountOutOfRange { threshold: usize, targets },
    MatchTableOutOfRange { kind, row: usize, column: usize },
    RingMembershipMissing { bond }, RingOutOfRange { ring: usize },
    MappedBondMissing { side: &'static str, begin, end }, StereoOrder(_),
    TooManyNewBonds { count: usize }, InitialSeedParse { message },
    InitialSeedMatch { message }, InitialSeedBondMissing { bond },
    ResultValueOutOfRange { kind, value: usize }, ResultQueryGraph { message },
    ResultSmarts { message }, SeedReconstructionBondMissing { begin, end },
    ResultContextMissing, Progress(_), State(_),
}, cause = |inner: &ck::McsError| Python::attach(|py| match inner {
    ck::McsError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).cause(py),
    ck::McsError::Progress(_) | ck::McsError::State(_)
    | ck::McsError::TargetTableCount { .. } | ck::McsError::ThresholdCountOutOfRange { .. }
    | ck::McsError::MatchTableOutOfRange { .. } | ck::McsError::RingMembershipMissing { .. }
    | ck::McsError::RingOutOfRange { .. } | ck::McsError::MappedBondMissing { .. }
    | ck::McsError::TooManyNewBonds { .. } | ck::McsError::InitialSeedParse { .. }
    | ck::McsError::InitialSeedMatch { .. } | ck::McsError::InitialSeedBondMissing { .. }
    | ck::McsError::ResultValueOutOfRange { .. } | ck::McsError::ResultQueryGraph { .. }
    | ck::McsError::ResultSmarts { .. } | ck::McsError::SeedReconstructionBondMissing { .. }
    | ck::McsError::ResultContextMissing => None,
}), source = Option<Source<crate::error_framework::ErrorObject<PropertyUIntReadError>, crate::error_framework::ErrorObject<AtomPropertyError>>>, methods = {
    #[getter]
    #[doc = "Number of targets associated with the failed count check, or None when not applicable."] fn targets(&self) -> Option<usize> { match self.inner {
        ck::McsError::TargetTableCount { targets, .. }
        | ck::McsError::ThresholdCountOutOfRange { targets, .. } => Some(targets), _ => None,
    } }
    #[getter]
    #[doc = "Native failure-specific category or value kind. Applies to McsError.MatchTableOutOfRange, McsError.ResultValueOutOfRange. None when this payload is absent."] fn kind(&self) -> Option<&'static str> { match self.inner {
        ck::McsError::MatchTableOutOfRange { kind, .. }
        | ck::McsError::ResultValueOutOfRange { kind, .. } => Some(kind), _ => None,
    } }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to McsError.RingMembershipMissing, McsError.InitialSeedBondMissing. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::McsError::RingMembershipMissing { bond }
        | ck::McsError::InitialSeedBondMissing { bond } => Some(bond), _ => None,
    } }
    #[getter]
    #[doc = "Beginning endpoint reported by the native check. Applies to McsError.MappedBondMissing, McsError.SeedReconstructionBondMissing. None when this payload is absent."] fn begin(&self) -> Option<usize> { match self.inner {
        ck::McsError::MappedBondMissing { begin, .. }
        | ck::McsError::SeedReconstructionBondMissing { begin, .. } => Some(begin), _ => None,
    } }
    #[getter]
    #[doc = "Ending endpoint reported by the native check. Applies to McsError.MappedBondMissing, McsError.SeedReconstructionBondMissing. None when this payload is absent."] fn end(&self) -> Option<usize> { match self.inner {
        ck::McsError::MappedBondMissing { end, .. }
        | ck::McsError::SeedReconstructionBondMissing { end, .. } => Some(end), _ => None,
    } }
    #[getter]
    #[doc = "Detailed native diagnostic message. Applies to McsError.InitialSeedParse, McsError.InitialSeedMatch, McsError.ResultQueryGraph, McsError.ResultSmarts. None when this payload is absent."] fn message(&self) -> Option<&str> { match &self.inner {
        ck::McsError::InitialSeedParse { message } | ck::McsError::InitialSeedMatch { message }
        | ck::McsError::ResultQueryGraph { message } | ck::McsError::ResultSmarts { message } => Some(message), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured StereoOrderError payload of McsError.StereoOrder; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn stereo_order(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<StereoOrderError>>> { match &self.inner {
        ck::McsError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured McsProgressError payload of McsError.Progress; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn progress(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<McsProgressError>>> { match &self.inner {
        ck::McsError::Progress(value) => McsProgressError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured SearchMcsError payload of McsError.State; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn state(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<SearchMcsError>>> { match &self.inner {
        ck::McsError::State(value) => SearchMcsError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_mcs_failure_preserves_state_progress_kind_fields_and_transparent_causes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles("CCO").unwrap();
            let before = molecule.to_binary().unwrap();
            let native = ck::maximum_common_substructure(&[&molecule]).unwrap_err();
            assert!(matches!(
                native,
                ck::McsError::State(ck::SearchMcsError::TooFewInputs { count: 1 })
            ));
            assert!(native.source().is_none());
            let display = native.to_string();
            let error = McsError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(error.value(py).getattr("source")?.is_none());
            let state = error.value(py).getattr("state")?;
            assert!(state.is_instance(<SearchMcsError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(state.getattr("count")?.extract::<usize>()?, 1);
            assert!(error.value(py).setattr("state", py.None()).is_err());
            assert_eq!(molecule.to_binary().unwrap(), before);

            let native = ck::McsError::StereoOrder(ck::StereoOrderError::ChiralPermutationRead(
                ck::PropertyUIntReadError::Negative { value: -7 },
            ));
            assert!(native.source().unwrap().is::<ck::PropertyUIntReadError>());
            let error = McsError::from_rust(native);
            let source = error.value(py).getattr("source")?;
            assert!(source.is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert_eq!(source.getattr("value")?.extract::<i64>()?, -7);
            assert!(
                error
                    .value(py)
                    .getattr("stereo_order")?
                    .is_instance(<StereoOrderError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            let native = ck::McsError::Progress(ck::McsProgressError::Callback {
                message: "callback\0context".into(),
            });
            assert!(native.source().is_none());
            let error = McsError::from_rust(native);
            assert!(error.value(py).getattr("source")?.is_none());
            let progress = error.value(py).getattr("progress")?;
            assert!(progress.is_instance(<McsProgressError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(
                progress.getattr("message")?.extract::<String>()?,
                "callback\0context"
            );

            let error = McsError::from_rust(ck::McsError::MatchTableOutOfRange {
                kind: "target",
                row: usize::MAX,
                column: 2,
            });
            assert_eq!(
                error.value(py).getattr("kind")?.extract::<String>()?,
                "target"
            );
            assert_eq!(
                error.value(py).getattr("row")?.extract::<usize>()?,
                usize::MAX
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<McsErrorKind>()?
                    == McsErrorKind::MatchTableOutOfRange
            );
            assert!(error.value(py).setattr("kind", "other").is_err());
            let error = SearchMcsError::from_rust(ck::SearchMcsError::SeedAtomNotBondEndpoint {
                atom: usize::MAX,
                bond: 3,
            });
            assert_eq!(
                error.value(py).getattr("atom")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(error.value(py).getattr("bond")?.extract::<usize>()?, 3);
            assert!(error.value(py).setattr("atom", 0).is_err());
            assert!(error.value(py).getattr("__cause__")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
