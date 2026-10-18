//! Complete native wedge/traversal dependencies for the SMARTS writer.
use crate::chemical_values::ChiralTag;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::AtomPropertyError;
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::potential_stereo_errors::{PotentialSource, PotentialStereoError};
use crate::property_numeric_errors::PropertyUIntReadError;
use crate::property_values::{PropertyValueKind, decode_source_text};
use crate::ring_errors::{RingFindingError, RingSource};
use crate::smiles_parse_errors::SmilesParseError;
use crate::stereo_errors::{AtropisomerError, LegacyStereoError, StereoError, StereoSource};
use crate::stereo_order_errors::StereoOrderError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type WedgeSource = Source<
    Source<StereoSource, PotentialSource>,
    Source<
        RingSource,
        Source<
            Py<TopologyValidationError>,
            Source<
                Py<CoordinateValidationError>,
                Source<Py<PropertyUIntReadError>, Py<AtomPropertyError>>,
            >,
        >,
    >,
>;

error_projection!(WedgeError, WedgeErrorKind, ck::WedgeError, "stereo", {
    SourcePrecondition { message: &'static str }, SourceStateIndex { state, index, count: usize },
    SourceIndexWidth { state, index }, SourceSignedScoreOverflow { atom, operation: &'static str },
    UnsignedRankOverflow { atom, property, value }, Topology(_), Coordinates(_),
    BondOutOfRange { bond, bond_count: usize }, NonSingleBond { bond }, CenterNotIncident { bond, center },
    UnsupportedCenterTag { center, tag: ChiralTag => |value: &ck::ChiralTag| (*value).into() },
    Geometry(_), StereoOrder(_), PotentialStereo(_), RingFinding(_), Atropisomer(_),
    InvalidStereoRank { atom, property, value }, InvalidStereoRankType { atom, property, kind },
}, cause = |inner: &ck::WedgeError| Python::attach(|py| match inner {
    ck::WedgeError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::WedgeError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::WedgeError::Geometry(value) => StereoError::from_rust(value.clone()).cause(py),
    ck::WedgeError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).cause(py),
    ck::WedgeError::PotentialStereo(value) => PotentialStereoError::from_rust(value.clone()).cause(py),
    ck::WedgeError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::WedgeError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).cause(py),
    ck::WedgeError::SourcePrecondition { .. } | ck::WedgeError::SourceStateIndex { .. }
    | ck::WedgeError::SourceIndexWidth { .. } | ck::WedgeError::SourceSignedScoreOverflow { .. }
    | ck::WedgeError::UnsignedRankOverflow { .. } | ck::WedgeError::BondOutOfRange { .. }
    | ck::WedgeError::NonSingleBond { .. } | ck::WedgeError::CenterNotIncident { .. }
    | ck::WedgeError::UnsupportedCenterTag { .. } | ck::WedgeError::InvalidStereoRank { .. }
    | ck::WedgeError::InvalidStereoRankType { .. } => None,
}), source = Option<WedgeSource>, methods = {
    #[getter]
    fn state(&self) -> Option<&'static str> { match self.inner {
        ck::WedgeError::SourceStateIndex { state, .. } | ck::WedgeError::SourceIndexWidth { state, .. } => Some(state), _ => None,
    } }
    #[getter]
    fn index(&self) -> Option<usize> { match self.inner {
        ck::WedgeError::SourceStateIndex { index, .. } | ck::WedgeError::SourceIndexWidth { index, .. } => Some(index), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::WedgeError::SourceSignedScoreOverflow { atom, .. } | ck::WedgeError::UnsignedRankOverflow { atom, .. }
        | ck::WedgeError::InvalidStereoRank { atom, .. } | ck::WedgeError::InvalidStereoRankType { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::WedgeError::UnsignedRankOverflow { property, .. } | ck::WedgeError::InvalidStereoRank { property, .. }
        | ck::WedgeError::InvalidStereoRankType { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    fn value(&self, py: Python<'_>) -> PyResult<Option<Source<u32, String>>> { match &self.inner {
        ck::WedgeError::UnsignedRankOverflow { value, .. } => Ok(Some(Source::First(*value))),
        ck::WedgeError::InvalidStereoRank { value, .. } => decode_source_text(py, value).map(|text| Some(Source::Second(text))), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::WedgeError::InvalidStereoRankType { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::WedgeError::BondOutOfRange { bond, .. } | ck::WedgeError::NonSingleBond { bond }
        | ck::WedgeError::CenterNotIncident { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    fn center(&self) -> Option<usize> { match self.inner {
        ck::WedgeError::CenterNotIncident { center, .. } | ck::WedgeError::UnsupportedCenterTag { center, .. } => Some(center.index()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<Py<TopologyValidationError>>> { match &self.inner {
        ck::WedgeError::Topology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn coordinates(&self, py: Python<'_>) -> PyResult<Option<Py<CoordinateValidationError>>> { match &self.inner {
        ck::WedgeError::Coordinates(value) => CoordinateValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn geometry(&self, py: Python<'_>) -> PyResult<Option<Py<StereoError>>> { match &self.inner {
        ck::WedgeError::Geometry(value) => StereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn stereo_order(&self, py: Python<'_>) -> PyResult<Option<Py<StereoOrderError>>> { match &self.inner {
        ck::WedgeError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn potential_stereo(&self, py: Python<'_>) -> PyResult<Option<Py<PotentialStereoError>>> { match &self.inner {
        ck::WedgeError::PotentialStereo(value) => PotentialStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<Py<RingFindingError>>> { match &self.inner {
        ck::WedgeError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atropisomer(&self, py: Python<'_>) -> PyResult<Option<Py<AtropisomerError>>> { match &self.inner {
        ck::WedgeError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(SmartsTraversalError, SmartsTraversalErrorKind, ck::SmartsTraversalError, "search", {
    Valence(_), Stereo(_), Traversal(_),
}, cause = |inner: &ck::SmartsTraversalError| match inner {
    ck::SmartsTraversalError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::SmartsTraversalError::Stereo(value) => Some(LegacyStereoError::from_rust(value.as_ref().clone())),
    ck::SmartsTraversalError::Traversal(value) => Some(SmilesParseError::from_rust(value.clone())),
}, source = Option<Source<Py<ValenceError>, Source<Py<LegacyStereoError>, Py<SmilesParseError>>>>);

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model_errors::StereoGroupError;
    use std::sync::Arc;

    #[test]
    fn wedge_and_smarts_traversal_keep_exact_causes_typed_payloads_and_rank_values() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::WedgeError::Geometry(ck::StereoError::StereoGroup(
                ck::StereoGroupError::DuplicateAtom,
            ));
            assert!(native.source().unwrap().is::<ck::StereoGroupError>());
            let display = native.to_string();
            let error = WedgeError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance_of::<StereoGroupError>());
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            let payload = error.value(py).getattr("geometry")?;
            assert!(payload.is_instance_of::<StereoError>());
            assert!(error.value(py).setattr("geometry", py.None()).is_err());
            let error = WedgeError::from_rust(ck::WedgeError::Coordinates(
                ck::CoordinateValidationError::RowCount {
                    dimension: "3D",
                    conformer: 3,
                    rows: 1,
                    atom_count: 2,
                },
            ));
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance_of::<CoordinateValidationError>()
            );
            let error = WedgeError::from_rust(ck::WedgeError::UnsignedRankOverflow {
                atom: ck::AtomId::new(usize::MAX),
                property: "_CIPRank",
                value: u32::MAX,
            });
            assert_eq!(
                error.value(py).getattr("atom")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(
                error.value(py).getattr("value")?.extract::<u32>()?,
                u32::MAX
            );
            let error = WedgeError::from_rust(ck::WedgeError::InvalidStereoRank {
                atom: ck::AtomId::new(0),
                property: "_CIPRank",
                value: "bad\0rank".into(),
            });
            assert_eq!(
                error.value(py).getattr("value")?.extract::<String>()?,
                "bad\0rank"
            );
            assert!(error.value(py).getattr("source")?.is_none());
            let error = WedgeError::from_rust(ck::WedgeError::InvalidStereoRankType {
                atom: ck::AtomId::new(0),
                property: "_CIPRank",
                kind: ck::PropertyValueKind::Bool,
            });
            assert!(
                error
                    .value(py)
                    .getattr("kind")?
                    .extract::<PropertyValueKind>()?
                    == PropertyValueKind::Bool
            );
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<WedgeErrorKind>()?
                    == WedgeErrorKind::InvalidStereoRankType
            );
            let cases = [
                (
                    ck::SmartsTraversalError::Valence(
                        ck::ValenceError::ExplicitValenceCacheNotInitialized {
                            atom: ck::AtomId::new(0),
                        },
                    ),
                    "ValenceError",
                ),
                (
                    ck::SmartsTraversalError::Stereo(Arc::new(
                        ck::LegacyStereoError::InvalidTopology(
                            ck::TopologyValidationError::AdjacencyMismatch,
                        ),
                    )),
                    "LegacyStereoError",
                ),
                (
                    ck::SmartsTraversalError::Traversal(ck::SmilesParseError::Syntax {
                        offset: 3,
                        message: "unclosed branch".into(),
                    }),
                    "SmilesParseError",
                ),
            ];
            for (native, class) in cases {
                assert!(native.source().is_some());
                let display = native.to_string();
                let error = SmartsTraversalError::from_rust(native);
                assert_eq!(error.value(py).str()?.to_str()?, display);
                let cause = error.value(py).getattr("source")?;
                assert_eq!(
                    cause
                        .getattr("__class__")?
                        .getattr("__name__")?
                        .extract::<String>()?,
                    class
                );
                assert!(cause.is(&error.value(py).getattr("__cause__")?));
            }
            Ok(())
        })
        .unwrap();
    }
}
