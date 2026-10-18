//! Exact native tautomer failures for borrowed views and owned snapshots.
use crate::canonical_errors::PropertyStringError;
use crate::error_framework::{Source, error_projection};
use crate::kekulize_errors::KekulizeError;
use crate::match_errors::{QueryMatchContextError, SubstructMatchError};
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::ring_errors::RingFindingError;
use crate::sanitize_errors::{SanitizeError, SanitizeSource};
use crate::smiles_parse_errors::{SmilesParseError, SmilesSource};
use crate::stereo_errors::{LegacySource, LegacyStereoError};
use crate::string_errors::QueryCompileError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type TautomerSource = Source<SmilesSource, Source<LegacySource, SanitizeSource>>;

error_projection!(TautomerRunError, TautomerRunErrorKind, ck::TautomerRunError, "tautomer", {
    Smiles(_), Valence(_), Rings(_), Kekulize(_), Sanitize(_), LegacyStereo(_),
    QueryCompile(_), Match(_), MatchContext(_), Coordinates(_), Topology(_),
    AtomProperty(_), MoleculeProperty(_), BondValue(_), PropertyString(_),
    AtomMappingCount { expected, actual }, BondMappingCount { expected, actual },
    AtomCountMismatch { expected, actual }, BondCountMismatch { expected, actual },
    AtomOutOfRange { atom, atom_count: usize }, BondOutOfRange { bond, bond_count: usize },
    EmptyTransformMatch, EditCount { field: &'static str, expected, actual },
    HydrogenCountOutOfRange { atom, count: u32 }, FormalChargeOutOfRange { atom, charge: i32 },
    MissingMappedBond { bond }, MissingValence { atom }, NoCanonicalTautomer,
    SourceCacheCommit, ScoreCacheTargetMismatch, Control(_), Callback(_),
}, cause = |inner: &ck::TautomerRunError| Python::attach(|py| match inner {
    ck::TautomerRunError::Smiles(value) => SmilesParseError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Rings(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Kekulize(value) => KekulizeError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Sanitize(value) => SanitizeError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::LegacyStereo(value) => LegacyStereoError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::QueryCompile(value) => QueryCompileError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Match(value) => SubstructMatchError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::MatchContext(value) => QueryMatchContextError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Coordinates(value) => CoordinateValidationError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::Topology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::MoleculeProperty(value) => MoleculePropertyError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::BondValue(value) => BondValueError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::PropertyString(value) => PropertyStringError::from_rust(value.clone()).cause(py),
    ck::TautomerRunError::AtomMappingCount { .. }
    | ck::TautomerRunError::BondMappingCount { .. }
    | ck::TautomerRunError::AtomCountMismatch { .. }
    | ck::TautomerRunError::BondCountMismatch { .. }
    | ck::TautomerRunError::AtomOutOfRange { .. }
    | ck::TautomerRunError::BondOutOfRange { .. }
    | ck::TautomerRunError::EmptyTransformMatch
    | ck::TautomerRunError::EditCount { .. }
    | ck::TautomerRunError::HydrogenCountOutOfRange { .. }
    | ck::TautomerRunError::FormalChargeOutOfRange { .. }
    | ck::TautomerRunError::MissingMappedBond { .. }
    | ck::TautomerRunError::MissingValence { .. }
    | ck::TautomerRunError::NoCanonicalTautomer
    | ck::TautomerRunError::SourceCacheCommit
    | ck::TautomerRunError::ScoreCacheTargetMismatch
    | ck::TautomerRunError::Control(_)
    | ck::TautomerRunError::Callback(_) => None,
}), source = Option<TautomerSource>, methods = {
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::TautomerRunError::AtomMappingCount { expected, .. }
        | ck::TautomerRunError::BondMappingCount { expected, .. }
        | ck::TautomerRunError::AtomCountMismatch { expected, .. }
        | ck::TautomerRunError::BondCountMismatch { expected, .. }
        | ck::TautomerRunError::EditCount { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::TautomerRunError::AtomMappingCount { actual, .. }
        | ck::TautomerRunError::BondMappingCount { actual, .. }
        | ck::TautomerRunError::AtomCountMismatch { actual, .. }
        | ck::TautomerRunError::BondCountMismatch { actual, .. }
        | ck::TautomerRunError::EditCount { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::TautomerRunError::AtomOutOfRange { atom, .. } => Some(atom),
        ck::TautomerRunError::HydrogenCountOutOfRange { atom, .. }
        | ck::TautomerRunError::FormalChargeOutOfRange { atom, .. }
        | ck::TautomerRunError::MissingValence { atom } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::TautomerRunError::BondOutOfRange { bond, .. }
        | ck::TautomerRunError::MissingMappedBond { bond } => Some(bond), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn smiles(&self, py: Python<'_>) -> PyResult<Option<Py<SmilesParseError>>> { match &self.inner {
        ck::TautomerRunError::Smiles(value) => SmilesParseError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::TautomerRunError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn rings(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> { match &self.inner {
        ck::TautomerRunError::Rings(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kekulize(&self, py: Python<'_>) -> PyResult<Option<Py<KekulizeError>>> { match &self.inner {
        ck::TautomerRunError::Kekulize(value) => KekulizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn sanitize(&self, py: Python<'_>) -> PyResult<Option<Py<SanitizeError>>> { match &self.inner {
        ck::TautomerRunError::Sanitize(value) => SanitizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn legacy_stereo(&self, py: Python<'_>) -> PyResult<Option<Py<LegacyStereoError>>> { match &self.inner {
        ck::TautomerRunError::LegacyStereo(value) => LegacyStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_compile(&self, py: Python<'_>) -> PyResult<Option<Py<QueryCompileError>>> { match &self.inner {
        ck::TautomerRunError::QueryCompile(value) => QueryCompileError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn match_error(&self, py: Python<'_>) -> PyResult<Option<Py<SubstructMatchError>>> { match &self.inner {
        ck::TautomerRunError::Match(value) => SubstructMatchError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn match_context(&self, py: Python<'_>) -> PyResult<Option<Py<QueryMatchContextError>>> { match &self.inner {
        ck::TautomerRunError::MatchContext(value) => QueryMatchContextError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinates(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> { match &self.inner {
        ck::TautomerRunError::Coordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::TautomerRunError::Topology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atom_property(&self, py: Python<'_>) -> PyResult<Option<Py<AtomPropertyError>>> { match &self.inner {
        ck::TautomerRunError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn molecule_property(&self, py: Python<'_>) -> PyResult<Option<Py<MoleculePropertyError>>> { match &self.inner {
        ck::TautomerRunError::MoleculeProperty(value) => MoleculePropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn bond_value(&self, py: Python<'_>) -> PyResult<Option<Py<BondValueError>>> { match &self.inner {
        ck::TautomerRunError::BondValue(value) => BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn property_string(&self, py: Python<'_>) -> PyResult<Option<Py<PropertyStringError>>> { match &self.inner {
        ck::TautomerRunError::PropertyString(value) => PropertyStringError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    fn control(&self) -> Option<&str> { match &self.inner {
        ck::TautomerRunError::Control(value) => Some(value), _ => None,
    } }
    #[getter]
    fn callback(&self) -> Option<&str> { match &self.inner {
        ck::TautomerRunError::Callback(value) => Some(value), _ => None,
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tautomer_errors_preserve_transparent_causes_wrapped_values_and_context() {
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
            let native = match molecule.enumerate_tautomers().unwrap_err() {
                ck::OperationError::Tautomer(value) => value,
                error => panic!("unexpected native error: {error:?}"),
            };
            let atoms = match &native {
                ck::TautomerRunError::Kekulize(ck::KekulizeError::NotKekulizable {
                    problem_atoms,
                }) => problem_atoms
                    .iter()
                    .map(|atom| atom.index())
                    .collect::<Vec<_>>(),
                error => panic!("unexpected tautomer error: {error:?}"),
            };
            let message = native.to_string();
            assert!(native.source().is_none());
            let error = TautomerRunError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<TautomerRunErrorKind>()?
                    == TautomerRunErrorKind::Kekulize
            );
            assert!(error.value(py).getattr("source")?.is_none());
            let payload = error.value(py).getattr("kekulize")?;
            assert!(payload.is_instance_of::<KekulizeError>());
            assert_eq!(
                payload.getattr("problem_atoms")?.extract::<Vec<usize>>()?,
                atoms
            );
            assert!(
                payload
                    .setattr("problem_atoms", Vec::<usize>::new())
                    .is_err()
            );

            // Transparent sanitize exposes its actual cause, not an invented
            // sanitize level; the full wrapped value remains available too.
            let native = ck::TautomerRunError::Sanitize(ck::SanitizeError::Properties {
                stage: ck::SanitizeStage::Properties,
                source: ck::PropertyCacheError::Valence(ck::ValenceError::InvalidTopology {
                    source: ck::TopologyValidationError::AdjacencyMismatch,
                }),
            });
            assert!(native.source().unwrap().is::<ck::PropertyCacheError>());
            let error = TautomerRunError::from_rust(native);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance_of::<crate::sanitize_support_errors::PropertyCacheError>());
            assert!(
                error
                    .value(py)
                    .getattr("sanitize")?
                    .is_instance_of::<SanitizeError>()
            );

            let error = TautomerRunError::from_rust(ck::TautomerRunError::FormalChargeOutOfRange {
                atom: ck::AtomId::new(7),
                charge: i32::MIN,
            });
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 7);
            assert_eq!(
                error.value(py).getattr("charge")?.extract::<i32>()?,
                i32::MIN
            );
            let error = TautomerRunError::from_rust(ck::TautomerRunError::Callback(
                "original callback detail".to_owned(),
            ));
            assert_eq!(
                error.value(py).getattr("callback")?.extract::<String>()?,
                "original callback detail"
            );
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
