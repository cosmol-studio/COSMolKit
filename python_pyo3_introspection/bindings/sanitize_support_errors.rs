//! Native sanitization support errors; transparent payloads are not cause edges.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::canonical_errors::CanonicalRankError;
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::{MappingValidationError, MoleculePropertyError, QueryStateError};
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use crate::ring_errors::RingFindingError;
use crate::selectors::AromaticityModel;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

type AromaticitySource = Source<
    Py<TopologyValidationError>,
    Source<
        Py<MappingValidationError>,
        Source<
            Py<MoleculePropertyError>,
            Source<Py<StereoGroupError>, Py<TemplateAttachmentOrderError>>,
        >,
    >,
>;
type CleanupSource = Source<
    Py<TopologyValidationError>,
    Source<
        Py<MoleculePropertyError>,
        Source<Py<StereoGroupError>, Py<TemplateAttachmentOrderError>>,
    >,
>;

error_projection!(AromaticityError, AromaticityErrorKind, ck::AromaticityError, "aromaticity", {
    InvalidTopology(_), InvalidQueryState(_), RingInfoNotInitialized,
    RingInfoDimensionMismatch { ring_atom_count: usize, ring_bond_count: usize, topology_atom_count: usize, topology_bond_count: usize },
    RingTableLengthMismatch { atom_rows: usize, bond_rows: usize },
    RingRowLengthMismatch { ring_index, atom_count, bond_count },
    RingAtomOutOfRange { ring_index, atom, atom_count }, RingBondOutOfRange { ring_index, bond, bond_count },
    AtomOutOfRange { atom, atom_count }, BondOutOfRange { bond, bond_count },
    ValenceAssignmentLength { field, actual, expected }, InvalidValenceRow { atom, field, value: i32 },
    ExpectedRingBondNotFound { begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
    UnsupportedModel { model, detail }, UnsupportedState { field, detail }, IntegerOverflow { field },
    UnexpectedTopologyShape { input_atoms: usize, input_bonds: usize, output_atoms: usize, output_bonds: usize },
    AtomIdentityChanged { row, expected, actual }, BondIdentityChanged { row, expected, actual },
    AromaticRingCountOutOfRange { model, actual, maximum: usize }, RingFinding(_), Valence(_), Kekulize(_),
}, cause = |inner: &ck::AromaticityError| Python::attach(|py| match inner {
    ck::AromaticityError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::AromaticityError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).cause(py),
    ck::AromaticityError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::AromaticityError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::AromaticityError::Kekulize(value) => KekulizeError::from_rust(value.clone()).cause(py),
    ck::AromaticityError::RingInfoNotInitialized | ck::AromaticityError::RingInfoDimensionMismatch { .. }
    | ck::AromaticityError::RingTableLengthMismatch { .. } | ck::AromaticityError::RingRowLengthMismatch { .. }
    | ck::AromaticityError::RingAtomOutOfRange { .. } | ck::AromaticityError::RingBondOutOfRange { .. }
    | ck::AromaticityError::AtomOutOfRange { .. } | ck::AromaticityError::BondOutOfRange { .. }
    | ck::AromaticityError::ValenceAssignmentLength { .. } | ck::AromaticityError::InvalidValenceRow { .. }
    | ck::AromaticityError::ExpectedRingBondNotFound { .. } | ck::AromaticityError::UnsupportedModel { .. }
    | ck::AromaticityError::UnsupportedState { .. } | ck::AromaticityError::IntegerOverflow { .. }
    | ck::AromaticityError::UnexpectedTopologyShape { .. } | ck::AromaticityError::AtomIdentityChanged { .. }
    | ck::AromaticityError::BondIdentityChanged { .. } | ck::AromaticityError::AromaticRingCountOutOfRange { .. } => None,
}), source = Option<AromaticitySource>, methods = {
    #[getter]
    fn ring_index(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::RingRowLengthMismatch { ring_index, .. }
        | ck::AromaticityError::RingAtomOutOfRange { ring_index, .. }
        | ck::AromaticityError::RingBondOutOfRange { ring_index, .. } => Some(ring_index), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::RingRowLengthMismatch { atom_count, .. }
        | ck::AromaticityError::RingAtomOutOfRange { atom_count, .. }
        | ck::AromaticityError::AtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn bond_count(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::RingRowLengthMismatch { bond_count, .. }
        | ck::AromaticityError::RingBondOutOfRange { bond_count, .. }
        | ck::AromaticityError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::RingAtomOutOfRange { atom, .. }
        | ck::AromaticityError::AtomOutOfRange { atom, .. }
        | ck::AromaticityError::InvalidValenceRow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::RingBondOutOfRange { bond, .. }
        | ck::AromaticityError::BondOutOfRange { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    fn field(&self) -> Option<&'static str> { match self.inner {
        ck::AromaticityError::ValenceAssignmentLength { field, .. }
        | ck::AromaticityError::InvalidValenceRow { field, .. }
        | ck::AromaticityError::UnsupportedState { field, .. }
        | ck::AromaticityError::IntegerOverflow { field } => Some(field), _ => None,
    } }
    #[getter]
    fn detail(&self) -> Option<&'static str> { match self.inner {
        ck::AromaticityError::UnsupportedModel { detail, .. }
        | ck::AromaticityError::UnsupportedState { detail, .. } => Some(detail), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::ValenceAssignmentLength { expected, .. } => Some(expected),
        ck::AromaticityError::AtomIdentityChanged { expected, .. } => Some(expected.index()),
        ck::AromaticityError::BondIdentityChanged { expected, .. } => Some(expected.index()), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::ValenceAssignmentLength { actual, .. }
        | ck::AromaticityError::AromaticRingCountOutOfRange { actual, .. } => Some(actual),
        ck::AromaticityError::AtomIdentityChanged { actual, .. } => Some(actual.index()),
        ck::AromaticityError::BondIdentityChanged { actual, .. } => Some(actual.index()), _ => None,
    } }
    #[getter]
    fn row(&self) -> Option<usize> { match self.inner {
        ck::AromaticityError::AtomIdentityChanged { row, .. }
        | ck::AromaticityError::BondIdentityChanged { row, .. } => Some(row), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn model(&self) -> Option<AromaticityModel> { match self.inner {
        ck::AromaticityError::UnsupportedModel { model, .. }
        | ck::AromaticityError::AromaticRingCountOutOfRange { model, .. } => Some(model.into()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::AromaticityError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_state(&self, py: Python<'_>) -> PyResult<Option<Py<QueryStateError>>> { match &self.inner {
        ck::AromaticityError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> { match &self.inner {
        ck::AromaticityError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::AromaticityError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kekulize(&self, py: Python<'_>) -> PyResult<Option<Py<KekulizeError>>> { match &self.inner {
        ck::AromaticityError::Kekulize(value) => KekulizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(ConjugationError, ConjugationErrorKind, ck::ConjugationError, "conjugation", {
    InvalidTopology(_), ValenceAssignmentLength { field, actual: usize, expected: usize },
    InvalidValenceRow { atom, field, value: i32 }, AtomOutOfRange { atom, atom_count: usize },
    BondOutOfRange { bond: usize => |value: &ck::BondId| value.index(), bond_count: usize },
    IntegerOverflow { field }, Aromaticity(_), Valence(_),
}, cause = |inner: &ck::ConjugationError| Python::attach(|py| match inner {
    ck::ConjugationError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::ConjugationError::Aromaticity(value) => AromaticityError::from_rust(value.clone()).cause(py),
    ck::ConjugationError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::ConjugationError::ValenceAssignmentLength { .. } | ck::ConjugationError::InvalidValenceRow { .. }
    | ck::ConjugationError::AtomOutOfRange { .. } | ck::ConjugationError::BondOutOfRange { .. }
    | ck::ConjugationError::IntegerOverflow { .. } => None,
}), source = Option<AromaticitySource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::ConjugationError::InvalidValenceRow { atom, .. }
        | ck::ConjugationError::AtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn field(&self) -> Option<&'static str> { match self.inner {
        ck::ConjugationError::ValenceAssignmentLength { field, .. }
        | ck::ConjugationError::InvalidValenceRow { field, .. }
        | ck::ConjugationError::IntegerOverflow { field } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::ConjugationError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn aromaticity(&self, py: Python<'_>) -> PyResult<Option<Py<AromaticityError>>> { match &self.inner {
        ck::ConjugationError::Aromaticity(value) => AromaticityError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::ConjugationError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(CleanupError, CleanupErrorKind, ck::CleanupError, "cleanup", {
    InvalidTopology(_), Valence(_), CanonicalRank(_),
}, cause = |inner: &ck::CleanupError| Python::attach(|py| match inner {
    ck::CleanupError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::CleanupError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::CleanupError::CanonicalRank(value) => CanonicalRankError::from_rust(value.clone()).cause(py),
}), source = Option<CleanupSource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::CleanupError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::CleanupError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn canonical_rank(&self, py: Python<'_>) -> PyResult<Option<Py<CanonicalRankError>>> { match &self.inner {
        ck::CleanupError::CanonicalRank(value) => CanonicalRankError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(PropertyCacheError, PropertyCacheErrorKind, ck::PropertyCacheError, "property_cache", {
    Valence(_),
}, cause = |inner: &ck::PropertyCacheError| match inner {
    ck::PropertyCacheError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
}, source = Option<Py<ValenceError>>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Py<ValenceError>> { match &self.inner {
        ck::PropertyCacheError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map_err(Into::into),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanitize_support_errors_keep_transparent_payloads_and_nontransparent_causes() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::ValenceError::InvalidTopology {
                source: ck::TopologyValidationError::AdjacencyMismatch,
            };
            let wrapped = ck::PropertyCacheError::Valence(native.clone());
            assert!(wrapped.source().unwrap().is::<ck::ValenceError>());
            let error = PropertyCacheError::from_rust(wrapped);
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance_of::<ValenceError>());
            assert!(
                source
                    .getattr("source")?
                    .is_instance_of::<TopologyValidationError>()
            );

            // Unlike PropertyCacheError, these native variants are transparent.
            let cleanup = ck::CleanupError::Valence(native.clone());
            let conjugation = ck::ConjugationError::Valence(native.clone());
            let aromaticity = ck::AromaticityError::Valence(native);
            assert!(
                cleanup
                    .source()
                    .unwrap()
                    .is::<ck::TopologyValidationError>()
            );
            assert!(
                conjugation
                    .source()
                    .unwrap()
                    .is::<ck::TopologyValidationError>()
            );
            assert!(
                aromaticity
                    .source()
                    .unwrap()
                    .is::<ck::TopologyValidationError>()
            );
            for error in [
                CleanupError::from_rust(cleanup),
                ConjugationError::from_rust(conjugation),
                AromaticityError::from_rust(aromaticity),
            ] {
                let source = error.value(py).getattr("source")?;
                assert!(source.is_instance_of::<TopologyValidationError>());
                assert!(source.is(&error.value(py).getattr("__cause__")?));
                assert!(
                    error
                        .value(py)
                        .getattr("valence")?
                        .is_instance_of::<ValenceError>()
                );
            }
            let native = ck::AromaticityError::AtomIdentityChanged {
                row: 2,
                expected: ck::AtomId::new(4),
                actual: ck::AtomId::new(8),
            };
            let message = native.to_string();
            let error = AromaticityError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(error.value(py).getattr("expected")?.extract::<usize>()?, 4);
            assert_eq!(error.value(py).getattr("actual")?.extract::<usize>()?, 8);
            assert_eq!(error.value(py).getattr("row")?.extract::<usize>()?, 2);
            assert!(error.cause(py).is_none());
            assert!(error.value(py).setattr("actual", 9).is_err());
            Ok(())
        })
        .unwrap();
    }
}
