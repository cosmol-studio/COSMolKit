//! Concrete prepared-context and substructure matcher failures.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::AtomPropertyError;
use crate::model_errors::TopologyValidationError;
use crate::path_errors::PeriodicTableError;
use crate::property_numeric_errors::{PropertyIntReadError, PropertyUIntReadError};
use crate::ring_errors::RingFindingError;
use crate::stereo_order_errors::StereoOrderError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

error_projection!(QueryMatchContextError, QueryMatchContextErrorKind, ck::QueryMatchContextError, "search", {
    InvalidTopology(_), ValenceRows { field, expected, actual }, Valence(_),
    ValencePrecondition { atom, field, getter }, RingPrecondition { getter, atom, bond },
    UninitializedRings, RingMembershipRows { field, expected, actual }, Rings(_),
}, cause = |inner: &ck::QueryMatchContextError| Python::attach(|py| match inner {
    ck::QueryMatchContextError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::QueryMatchContextError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::QueryMatchContextError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::QueryMatchContextError::ValenceRows { .. } | ck::QueryMatchContextError::ValencePrecondition { .. }
    | ck::QueryMatchContextError::RingPrecondition { .. } | ck::QueryMatchContextError::UninitializedRings
    | ck::QueryMatchContextError::RingMembershipRows { .. } => None,
}), source = Option<Source<crate::error_framework::ErrorObject<TopologyValidationError>, crate::error_framework::ErrorObject<RingFindingError>>>, methods = {
    #[getter]
    #[doc = "Field whose content or type caused the failure. Applies to QueryMatchContextError.ValenceRows, QueryMatchContextError.ValencePrecondition, QueryMatchContextError.RingMembershipRows. None when this payload is absent."] fn field(&self) -> Option<&'static str> { match self.inner {
        ck::QueryMatchContextError::ValenceRows { field, .. }
        | ck::QueryMatchContextError::ValencePrecondition { field, .. }
        | ck::QueryMatchContextError::RingMembershipRows { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = "Accessor involved in the reported failure. Applies to QueryMatchContextError.ValencePrecondition, QueryMatchContextError.RingPrecondition. None when this payload is absent."] fn getter(&self) -> Option<&'static str> { match self.inner {
        ck::QueryMatchContextError::ValencePrecondition { getter, .. }
        | ck::QueryMatchContextError::RingPrecondition { getter, .. } => Some(getter), _ => None,
    } }
    #[getter]
    #[doc = "Required value reported by the native check. Applies to QueryMatchContextError.ValenceRows, QueryMatchContextError.RingMembershipRows. None when this payload is absent."] fn expected(&self) -> Option<usize> { match self.inner {
        ck::QueryMatchContextError::ValenceRows { expected, .. }
        | ck::QueryMatchContextError::RingMembershipRows { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = "Observed value reported by the native check. Applies to QueryMatchContextError.ValenceRows, QueryMatchContextError.RingMembershipRows. None when this payload is absent."] fn actual(&self) -> Option<usize> { match self.inner {
        ck::QueryMatchContextError::ValenceRows { actual, .. }
        | ck::QueryMatchContextError::RingMembershipRows { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to QueryMatchContextError.ValencePrecondition, QueryMatchContextError.RingPrecondition. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::QueryMatchContextError::ValencePrecondition { atom, .. } => Some(atom),
        ck::QueryMatchContextError::RingPrecondition { atom, .. } => atom, _ => None,
    } }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to QueryMatchContextError.RingPrecondition. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::QueryMatchContextError::RingPrecondition { bond, .. } => bond, _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured ValenceError payload of QueryMatchContextError.Valence; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::QueryMatchContextError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

pub(crate) type MatchSource = Source<
    crate::error_framework::ErrorObject<PropertyUIntReadError>,
    Source<
        crate::error_framework::ErrorObject<AtomPropertyError>,
        Source<
            crate::error_framework::ErrorObject<PropertyIntReadError>,
            Source<
                crate::error_framework::ErrorObject<TopologyValidationError>,
                crate::error_framework::ErrorObject<RingFindingError>,
            >,
        >,
    >,
>;

error_projection!(MatchError, MatchErrorKind, ck::MatchError, "search", {
    InvalidQuery(_), InvalidTarget(_), UnsupportedAtomPredicate(_), UnsupportedBondPredicate(_), Substruct(_),
}, cause = |inner: &ck::MatchError| Python::attach(|py| match inner {
    ck::MatchError::Substruct(value) => SubstructMatchError::from_rust(value.clone()).cause(py),
    ck::MatchError::InvalidQuery(_) | ck::MatchError::InvalidTarget(_)
    | ck::MatchError::UnsupportedAtomPredicate(_) | ck::MatchError::UnsupportedBondPredicate(_) => None,
}), source = Option<MatchSource>, methods = {
    #[getter]
    #[doc = "Explanation of an invalid query, or None for another failure."] fn invalid_query(&self) -> Option<String> { match &self.inner {
        ck::MatchError::InvalidQuery(value) => Some(value.clone()), _ => None,
    } }
    #[getter]
    #[doc = "Explanation of an invalid target, or None for another failure."] fn invalid_target(&self) -> Option<String> { match &self.inner {
        ck::MatchError::InvalidTarget(value) => Some(value.clone()), _ => None,
    } }
    #[getter]
    #[doc = "Unsupported atom predicate name; None for other failures."] fn unsupported_atom_predicate(&self) -> Option<&'static str> { match &self.inner {
        ck::MatchError::UnsupportedAtomPredicate(value) => Some(*value), _ => None,
    } }
    #[getter]
    #[doc = "Unsupported bond predicate name; None for other failures."] fn unsupported_bond_predicate(&self) -> Option<&'static str> { match &self.inner {
        ck::MatchError::UnsupportedBondPredicate(value) => Some(*value), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured SubstructMatchError payload of MatchError.Substruct; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn substruct(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<SubstructMatchError>>> { match &self.inner {
        ck::MatchError::Substruct(value) => SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(SubstructMatchError, SubstructMatchErrorKind, ck::SubstructMatchError, "search", {
    Unsupported { branch: &'static str, rdkit_function: &'static str },
    FinalCheckMappingLength { query_mapping: usize, target_mapping: usize, query_atoms: usize },
    FinalCheckMappingIndex { side, position: usize, index, atom_count },
    FinalCheckInvariant { invariant: &'static str, query_atom: usize },
    FinalCheckBondEndpoint { side, endpoint: &'static str, index, atom_count },
    FinalCheckMissingBond { query_bond: usize, begin: usize, end: usize },
    PeriodicTable(_), StereoOrder(_), PropertyString(_),
    PropertyInteger { property: &'static str, source }, QueryContext(_), UninitializedRecursiveQuery,
}, cause = |inner: &ck::SubstructMatchError| Python::attach(|py| match inner {
    ck::SubstructMatchError::PeriodicTable(value) => PeriodicTableError::from_rust(value.clone()).cause(py),
    ck::SubstructMatchError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).cause(py),
    ck::SubstructMatchError::PropertyString(value) => PropertyStringError::from_rust(*value).cause(py),
    ck::SubstructMatchError::PropertyInteger { source, .. } => Some(PropertyIntReadError::from_rust(source.clone())),
    ck::SubstructMatchError::QueryContext(value) => QueryMatchContextError::from_rust(value.clone()).cause(py),
    ck::SubstructMatchError::Unsupported { .. } | ck::SubstructMatchError::FinalCheckMappingLength { .. }
    | ck::SubstructMatchError::FinalCheckMappingIndex { .. } | ck::SubstructMatchError::FinalCheckInvariant { .. }
    | ck::SubstructMatchError::FinalCheckBondEndpoint { .. } | ck::SubstructMatchError::FinalCheckMissingBond { .. }
    | ck::SubstructMatchError::UninitializedRecursiveQuery => None,
}), source = Option<MatchSource>, methods = {
    #[getter]
    #[doc = "Probe/reference or source/target side identified by the native check. Applies to SubstructMatchError.FinalCheckMappingIndex, SubstructMatchError.FinalCheckBondEndpoint. None when this payload is absent."] fn side(&self) -> Option<&'static str> { match self.inner {
        ck::SubstructMatchError::FinalCheckMappingIndex { side, .. }
        | ck::SubstructMatchError::FinalCheckBondEndpoint { side, .. } => Some(side), _ => None,
    } }
    #[getter]
    #[doc = "Index reported by the failing native check. Applies to SubstructMatchError.FinalCheckMappingIndex, SubstructMatchError.FinalCheckBondEndpoint. None when this payload is absent."] fn index(&self) -> Option<usize> { match self.inner {
        ck::SubstructMatchError::FinalCheckMappingIndex { index, .. }
        | ck::SubstructMatchError::FinalCheckBondEndpoint { index, .. } => Some(index), _ => None,
    } }
    #[getter]
    #[doc = "Atom count used when validating the reported atom reference. Applies to SubstructMatchError.FinalCheckMappingIndex, SubstructMatchError.FinalCheckBondEndpoint. None when this payload is absent."] fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::SubstructMatchError::FinalCheckMappingIndex { atom_count, .. }
        | ck::SubstructMatchError::FinalCheckBondEndpoint { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured PeriodicTableError payload of SubstructMatchError.PeriodicTable; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn periodic_table(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PeriodicTableError>>> { match &self.inner {
        ck::SubstructMatchError::PeriodicTable(value) => PeriodicTableError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured StereoOrderError payload of SubstructMatchError.StereoOrder; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn stereo_order(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<StereoOrderError>>> { match &self.inner {
        ck::SubstructMatchError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured PropertyStringError payload of SubstructMatchError.PropertyString; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn property_string(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PropertyStringError>>> { match &self.inner {
        ck::SubstructMatchError::PropertyString(value) => PropertyStringError::from_rust(*value).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured QueryMatchContextError payload of SubstructMatchError.QueryContext; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_context(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryMatchContextError>>> { match &self.inner {
        ck::SubstructMatchError::QueryContext(value) => QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matcher_projection_preserves_transparent_context_and_full_width_rows() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::SubstructMatchError::QueryContext(
                ck::QueryMatchContextError::Valence(ck::ValenceError::InvalidTopology {
                    source: ck::TopologyValidationError::AdjacencyMismatch,
                }),
            );
            assert!(native.source().unwrap().is::<ck::TopologyValidationError>());
            let display = native.to_string();
            let error = SubstructMatchError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let context = error.value(py).getattr("query_context")?;
            assert!(context.is_instance(<QueryMatchContextError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(context.getattr("valence")?.is_instance(<ValenceError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);

            let native = ck::SubstructMatchError::FinalCheckMappingIndex {
                side: "target",
                position: 2,
                index: usize::MAX,
                atom_count: 3,
            };
            assert!(native.source().is_none());
            let error = SubstructMatchError::from_rust(native);
            assert_eq!(
                error.value(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                error.value(py).getattr("side")?.extract::<String>()?,
                "target"
            );
            assert!(error.value(py).setattr("index", 0).is_err());
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
