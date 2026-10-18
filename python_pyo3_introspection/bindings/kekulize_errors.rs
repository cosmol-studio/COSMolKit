//! The native kekulization error and its precise payload/source distinction.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::canonical_errors::CanonicalRankError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{MappingValidationError, MoleculePropertyError, QueryStateError};
use crate::model_errors::{StereoGroupError, TopologyValidationError};
use crate::ring_errors::RingFindingError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type KekulizeSource = Source<
    Py<TopologyValidationError>,
    Source<
        Py<MappingValidationError>,
        Source<
            Py<MoleculePropertyError>,
            Source<Py<StereoGroupError>, Py<TemplateAttachmentOrderError>>,
        >,
    >,
>;

error_projection!(KekulizeError, KekulizeErrorKind, ck::KekulizeError, "kekulize", {
    InvalidTopology(_), AtomSelectionLength { expected, actual }, BondSelectionLength { expected, actual },
    CandidateAtomOutOfRange { atom, atom_count: usize }, DuplicateCandidateAtom { atom },
    MatchingStateLength { field, expected, actual }, DuplicateDoneAtom { atom },
    MissingCandidateBond { begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
    MissingBacktrackAnchor { atom }, QuestionSubsetOverflow { questions: usize, bit_width: u32 },
    AromaticAtomOutsideRing { atom },
    NotKekulizable { problem_atoms: Vec<usize> => |values: &Vec<ck::AtomId>| values.iter().map(|value| value.index()).collect() },
    PostconditionValenceMismatch { atom, before: i32, after: i32 },
    UnsupportedQueryState { bond: usize => |value: &ck::BondId| value.index(), detail: &'static str },
    InvalidQueryState(_), IntegerOverflow { atom, field }, RingFinding(_), Valence(_), CanonicalRank(_),
}, cause = |inner: &ck::KekulizeError| Python::attach(|py| match inner {
    ck::KekulizeError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    // Transparent wrappers forward the wrapped error's real source, not the
    // wrapped value itself. That value remains available as a typed payload.
    ck::KekulizeError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).cause(py),
    ck::KekulizeError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::KekulizeError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::KekulizeError::CanonicalRank(value) => CanonicalRankError::from_rust(value.clone()).cause(py),
    ck::KekulizeError::AtomSelectionLength { .. } | ck::KekulizeError::BondSelectionLength { .. }
    | ck::KekulizeError::CandidateAtomOutOfRange { .. } | ck::KekulizeError::DuplicateCandidateAtom { .. }
    | ck::KekulizeError::MatchingStateLength { .. } | ck::KekulizeError::DuplicateDoneAtom { .. }
    | ck::KekulizeError::MissingCandidateBond { .. } | ck::KekulizeError::MissingBacktrackAnchor { .. }
    | ck::KekulizeError::QuestionSubsetOverflow { .. } | ck::KekulizeError::AromaticAtomOutsideRing { .. }
    | ck::KekulizeError::NotKekulizable { .. } | ck::KekulizeError::PostconditionValenceMismatch { .. }
    | ck::KekulizeError::UnsupportedQueryState { .. } | ck::KekulizeError::IntegerOverflow { .. } => None,
}), source = Option<KekulizeSource>, methods = {
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::KekulizeError::AtomSelectionLength { expected, .. }
        | ck::KekulizeError::BondSelectionLength { expected, .. }
        | ck::KekulizeError::MatchingStateLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::KekulizeError::AtomSelectionLength { actual, .. }
        | ck::KekulizeError::BondSelectionLength { actual, .. }
        | ck::KekulizeError::MatchingStateLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::KekulizeError::CandidateAtomOutOfRange { atom, .. }
        | ck::KekulizeError::DuplicateCandidateAtom { atom }
        | ck::KekulizeError::DuplicateDoneAtom { atom }
        | ck::KekulizeError::MissingBacktrackAnchor { atom }
        | ck::KekulizeError::AromaticAtomOutsideRing { atom }
        | ck::KekulizeError::PostconditionValenceMismatch { atom, .. }
        | ck::KekulizeError::IntegerOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn field(&self) -> Option<&'static str> { match self.inner {
        ck::KekulizeError::MatchingStateLength { field, .. }
        | ck::KekulizeError::IntegerOverflow { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::KekulizeError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_state(&self, py: Python<'_>) -> PyResult<Option<Py<QueryStateError>>> { match &self.inner {
        ck::KekulizeError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> { match &self.inner {
        ck::KekulizeError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::KekulizeError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn canonical_rank(&self, py: Python<'_>) -> PyResult<Option<Py<CanonicalRankError>>> { match &self.inner {
        ck::KekulizeError::CanonicalRank(value) => CanonicalRankError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_kekulization_failure_retains_atom_rows_and_transparent_payload() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let molecule = ck::Molecule::from_smiles_with_params(
                "c1cccc1",
                &ck::SmilesParseParams {
                    sanitize: false,
                    ..Default::default()
                },
            )
            .unwrap();
            let ck::OperationError::Tautomer(ck::TautomerRunError::Kekulize(native)) =
                molecule.enumerate_tautomers().unwrap_err()
            else {
                panic!("native kekulization error")
            };
            let ck::KekulizeError::NotKekulizable { problem_atoms } = &native else {
                panic!("native unresolved atoms")
            };
            let atoms = problem_atoms
                .iter()
                .map(|id| id.index())
                .collect::<Vec<_>>();
            assert!(!atoms.is_empty());
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = KekulizeError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<KekulizeErrorKind>()?
                    == KekulizeErrorKind::NotKekulizable
            );
            assert_eq!(
                error
                    .value(py)
                    .getattr("problem_atoms")?
                    .extract::<Vec<usize>>()?,
                atoms
            );
            assert!(
                error
                    .value(py)
                    .setattr("problem_atoms", Vec::<usize>::new())
                    .is_err()
            );
            assert!(error.cause(py).is_none());
            assert!(error.value(py).getattr("source")?.is_none());

            let wrapped = ck::KekulizeError::InvalidQueryState(ck::QueryStateError::AtomCount {
                actual: 2,
                expected: 3,
            });
            assert!(wrapped.source().is_none());
            let error = KekulizeError::from_rust(wrapped);
            assert!(error.cause(py).is_none());
            let value = error.value(py).getattr("query_state")?;
            assert!(value.is_instance_of::<QueryStateError>());
            assert_eq!(value.getattr("actual")?.extract::<usize>()?, 2);
            assert_eq!(value.getattr("expected")?.extract::<usize>()?, 3);
            assert!(value.getattr("__cause__")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
