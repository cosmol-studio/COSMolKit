//! Concrete sanitization leaf errors, forwarding the native topology causes.
use crate::error_framework::error_projection;
use crate::model_errors::TopologyValidationError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(AdjustHsError, AdjustHsErrorKind, ck::AdjustHsError, "adjust_hs", {
    InvalidTopology { source }, ValenceAssignmentLength { field, actual: usize, expected: usize },
    InvalidValenceRow { atom, field, value: i32 },
    ExplicitHydrogenOverflow { atom, original_explicit_hydrogens: u8, original_implicit_valence: i32, recalculated_implicit_valence: i32 },
    Valence(_),
}, cause = |inner: &ck::AdjustHsError| Python::attach(|py| match inner {
    ck::AdjustHsError::InvalidTopology { source } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::AdjustHsError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::AdjustHsError::ValenceAssignmentLength { .. } | ck::AdjustHsError::InvalidValenceRow { .. }
    | ck::AdjustHsError::ExplicitHydrogenOverflow { .. } => None,
}), source = Option<crate::error_framework::ErrorObject<TopologyValidationError>>, methods = {
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to AdjustHsError.InvalidValenceRow, AdjustHsError.ExplicitHydrogenOverflow. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::AdjustHsError::InvalidValenceRow { atom, .. }
        | ck::AdjustHsError::ExplicitHydrogenOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    #[doc = "Field whose content or type caused the failure. Applies to AdjustHsError.ValenceAssignmentLength, AdjustHsError.InvalidValenceRow. None when this payload is absent."] fn field(&self) -> Option<&'static str> { match self.inner {
        ck::AdjustHsError::ValenceAssignmentLength { field, .. }
        | ck::AdjustHsError::InvalidValenceRow { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of AdjustHsError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::AdjustHsError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(RadicalError, RadicalErrorKind, ck::RadicalError, "radicals", {
    InvalidTopology { source }, PeriodicTableLookup { atomic_number: u8, field: &'static str },
    Valence(_), RadicalCountOutOfRange { atom: usize => |value: &ck::AtomId| value.index(), count: i32 },
}, cause = |inner: &ck::RadicalError| Python::attach(|py| match inner {
    ck::RadicalError::InvalidTopology { source } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::RadicalError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::RadicalError::PeriodicTableLookup { .. } | ck::RadicalError::RadicalCountOutOfRange { .. } => None,
}), source = Option<crate::error_framework::ErrorObject<TopologyValidationError>>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of RadicalError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::RadicalError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(HybridizationError, HybridizationErrorKind, ck::HybridizationError, "hybridization", {
    ConjugatedBondAssignmentLength { actual, expected }, InvalidTopology(_),
    ValenceAssignmentLength { field, actual, expected }, InvalidValenceRow { atom, field, value: i32 },
    IntegerOverflow { atom, field }, Valence(_),
}, cause = |inner: &ck::HybridizationError| Python::attach(|py| match inner {
    ck::HybridizationError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::HybridizationError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::HybridizationError::ConjugatedBondAssignmentLength { .. } | ck::HybridizationError::ValenceAssignmentLength { .. }
    | ck::HybridizationError::InvalidValenceRow { .. } | ck::HybridizationError::IntegerOverflow { .. } => None,
}), source = Option<crate::error_framework::ErrorObject<TopologyValidationError>>, methods = {
    #[getter]
    #[doc = "Observed value reported by the native check. Applies to HybridizationError.ConjugatedBondAssignmentLength, HybridizationError.ValenceAssignmentLength. None when this payload is absent."] fn actual(&self) -> Option<usize> { match self.inner {
        ck::HybridizationError::ConjugatedBondAssignmentLength { actual, .. }
        | ck::HybridizationError::ValenceAssignmentLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    #[doc = "Required value reported by the native check. Applies to HybridizationError.ConjugatedBondAssignmentLength, HybridizationError.ValenceAssignmentLength. None when this payload is absent."] fn expected(&self) -> Option<usize> { match self.inner {
        ck::HybridizationError::ConjugatedBondAssignmentLength { expected, .. }
        | ck::HybridizationError::ValenceAssignmentLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = "Field whose content or type caused the failure. Applies to HybridizationError.ValenceAssignmentLength, HybridizationError.InvalidValenceRow, HybridizationError.IntegerOverflow. None when this payload is absent."] fn field(&self) -> Option<&'static str> { match self.inner {
        ck::HybridizationError::ValenceAssignmentLength { field, .. }
        | ck::HybridizationError::InvalidValenceRow { field, .. }
        | ck::HybridizationError::IntegerOverflow { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to HybridizationError.InvalidValenceRow, HybridizationError.IntegerOverflow. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::HybridizationError::InvalidValenceRow { atom, .. }
        | ck::HybridizationError::IntegerOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of HybridizationError.InvalidTopology; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::HybridizationError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of HybridizationError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::HybridizationError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanitize_leaf_projection_keeps_native_payloads_and_real_source_edges() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let topology = ck::TopologyValidationError::AdjacencyMismatch;
            let native_hs = ck::AdjustHsError::InvalidTopology {
                source: topology.clone(),
            };
            let native_radicals = ck::RadicalError::InvalidTopology {
                source: topology.clone(),
            };
            let native_hybrid = ck::HybridizationError::InvalidTopology(topology);
            let messages = [
                native_hs.to_string(),
                native_radicals.to_string(),
                native_hybrid.to_string(),
            ];
            let errors = [
                AdjustHsError::from_rust(native_hs),
                RadicalError::from_rust(native_radicals),
                HybridizationError::from_rust(native_hybrid),
            ];
            for (error, expected) in errors.iter().zip(messages) {
                assert_eq!(error.value(py).str()?.to_str()?, expected);
                let source = error.value(py).getattr("source")?;
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert!(source.is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
                assert!(source.getattr("__cause__")?.is_none());
            }
            let native = ck::AdjustHsError::ExplicitHydrogenOverflow {
                atom: ck::AtomId::new(2),
                original_explicit_hydrogens: 250,
                original_implicit_valence: 6,
                recalculated_implicit_valence: 0,
            };
            assert!(native.source().is_none());
            let error = AdjustHsError::from_rust(native);
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 2);
            assert_eq!(
                error
                    .value(py)
                    .getattr("original_explicit_hydrogens")?
                    .extract::<u8>()?,
                250
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("original_implicit_valence")?
                    .extract::<i32>()?,
                6
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("recalculated_implicit_valence")?
                    .extract::<i32>()?,
                0
            );
            assert!(error.cause(py).is_none());
            assert!(error.value(py).setattr("atom", 8).is_err());
            Ok(())
        })
        .unwrap();
    }
}
