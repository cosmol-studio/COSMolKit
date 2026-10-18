//! Exact hydrogen error fields; only four native variants have source edges.
use crate::error_framework::{Source, error_projection};
use crate::metadata_values::SdfPropertyListTarget;
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
    QueryStateError,
};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::property_values::decode_source_text;
use crate::sanitize_errors::SanitizeError;
use crate::stereo_order_errors::StereoOrderError;
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type HydrogenSource = Source<
    Py<StereoOrderError>,
    Source<Py<AtomPropertyError>, Source<Py<BondValueError>, Py<MoleculePropertyError>>>,
>;

error_projection!(HydrogenError, HydrogenErrorKind, ck::HydrogenError, "hydrogens", {
    StereoOrder(_), BondProperty(_), AtomProperty(_), InvalidTopology(_), InvalidCoordinates(_),
    TopologyEdit(_), InvalidMapping(_), InvalidQueryState(_), InvalidProperty(_),
    InvalidPropertyList { target, name, expected_rows: usize, actual_rows: usize },
    Unsupported { operation: &'static str, reason }, OnlyOnAtomOutOfRange { atom, atom_count: usize },
    InvalidAdditionPlan { addition, reason }, InvalidRemovalCandidate { position: usize, atom, reason },
    MissingSourceState { stage: &'static str, state: &'static str },
    ValenceAssignmentLength { field: &'static str, expected: usize, actual: usize },
    ExplicitHydrogenOverflow { atom, current: u8 },
    InvalidStereoTransition { bond: usize => |value: &ck::BondId| value.index(), reason },
    CoordinatePlacement { addition, conformer: usize, dimension: &'static str, reason }, Valence(_), Sanitize(_),
}, cause = |inner: &ck::HydrogenError| match inner {
    ck::HydrogenError::StereoOrder(value) => Some(StereoOrderError::from_rust(value.clone())),
    ck::HydrogenError::BondProperty(value) => Some(BondValueError::from_rust(value.clone())),
    ck::HydrogenError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::HydrogenError::InvalidProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::HydrogenError::InvalidTopology(_) | ck::HydrogenError::InvalidCoordinates(_)
    | ck::HydrogenError::TopologyEdit(_) | ck::HydrogenError::InvalidMapping(_)
    | ck::HydrogenError::InvalidQueryState(_) | ck::HydrogenError::InvalidPropertyList { .. }
    | ck::HydrogenError::Unsupported { .. } | ck::HydrogenError::OnlyOnAtomOutOfRange { .. }
    | ck::HydrogenError::InvalidAdditionPlan { .. } | ck::HydrogenError::InvalidRemovalCandidate { .. }
    | ck::HydrogenError::MissingSourceState { .. } | ck::HydrogenError::ValenceAssignmentLength { .. }
    | ck::HydrogenError::ExplicitHydrogenOverflow { .. } | ck::HydrogenError::InvalidStereoTransition { .. }
    | ck::HydrogenError::CoordinatePlacement { .. } | ck::HydrogenError::Valence(_)
    | ck::HydrogenError::Sanitize(_) => None,
}, source = Option<HydrogenSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::HydrogenError::OnlyOnAtomOutOfRange { atom, .. }
        | ck::HydrogenError::InvalidRemovalCandidate { atom, .. }
        | ck::HydrogenError::ExplicitHydrogenOverflow { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn reason(&self) -> Option<&'static str> { match self.inner {
        ck::HydrogenError::Unsupported { reason, .. } | ck::HydrogenError::InvalidAdditionPlan { reason, .. }
        | ck::HydrogenError::InvalidRemovalCandidate { reason, .. } | ck::HydrogenError::InvalidStereoTransition { reason, .. }
        | ck::HydrogenError::CoordinatePlacement { reason, .. } => Some(reason), _ => None,
    } }
    #[getter]
    fn addition(&self) -> Option<usize> { match self.inner {
        ck::HydrogenError::InvalidAdditionPlan { addition, .. } => addition,
        ck::HydrogenError::CoordinatePlacement { addition, .. } => Some(addition), _ => None,
    } }
    #[getter]
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> { match &self.inner {
        ck::HydrogenError::InvalidPropertyList { name, .. } => decode_source_text(py, name).map(Some), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn target(&self) -> Option<SdfPropertyListTarget> { match self.inner {
        ck::HydrogenError::InvalidPropertyList { target, .. } => Some(target.into()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::HydrogenError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinates(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> { match &self.inner {
        ck::HydrogenError::InvalidCoordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology_edit(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyEditError>>> { match &self.inner {
        ck::HydrogenError::TopologyEdit(value) => TopologyEditError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn mapping(&self, py: Python<'_>) -> PyResult<Option<Py<MappingValidationError>>> { match &self.inner {
        ck::HydrogenError::InvalidMapping(value) => MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_state(&self, py: Python<'_>) -> PyResult<Option<Py<QueryStateError>>> { match &self.inner {
        ck::HydrogenError::InvalidQueryState(value) => QueryStateError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<Py<ValenceError>>> { match &self.inner {
        ck::HydrogenError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn sanitize(&self, py: Python<'_>) -> PyResult<Option<Py<SanitizeError>>> { match &self.inner {
        ck::HydrogenError::Sanitize(value) => SanitizeError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_hydrogen_rejection_and_source_free_sanitize_payload_keep_their_shape() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Molecule::from_smiles("C").unwrap();
            let native = source
                .with_hydrogens_with_params(&ck::AddHsParams {
                    only_on_atoms: Some(vec![ck::AtomId::new(9)]),
                    ..Default::default()
                })
                .unwrap_err();
            let ck::OperationError::Hydrogen(native) = native else {
                panic!("native hydrogen rejection");
            };
            assert!(matches!(
                native,
                ck::HydrogenError::OnlyOnAtomOutOfRange { .. }
            ));
            let message = native.to_string();
            let error = HydrogenError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, message);
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 9);
            assert_eq!(
                error.value(py).getattr("atom_count")?.extract::<usize>()?,
                1
            );
            assert!(error.value(py).setattr("atom", 0).is_err());
            assert!(error.cause(py).is_none());

            let native = ck::HydrogenError::Sanitize(ck::SanitizeError::InvalidOperations {
                bits: 1 << 20,
                unknown_bits: 1 << 20,
            });
            assert!(native.source().is_none());
            let error = HydrogenError::from_rust(native);
            assert!(error.cause(py).is_none());
            let payload = error.value(py).getattr("sanitize")?;
            assert!(payload.is_instance_of::<SanitizeError>());
            assert_eq!(payload.getattr("bits")?.extract::<u32>()?, 1 << 20);

            let native = ck::HydrogenError::AtomProperty(ck::AtomPropertyError::EmptyKey);
            assert!(native.source().unwrap().is::<ck::AtomPropertyError>());
            let error = HydrogenError::from_rust(native);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<AtomPropertyError>());
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            Ok(())
        })
        .unwrap();
    }
}
