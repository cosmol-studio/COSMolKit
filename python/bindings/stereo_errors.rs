//! Native stereo, axial-stereo and legacy failures without invented cause edges.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::canonical_errors::CanonicalRankError;
use crate::chemical_values::BondStereo;
use crate::double_stereo_errors::DoubleBondStereoError;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::molecular_hash::CipRankError;
use crate::potential_stereo_errors::PotentialStereoError;
use crate::property_numeric_errors::{PropertyIntReadError, PropertyUIntReadError};
use crate::ring_errors::RingFindingError;
use crate::sanitize_leaf_errors::HybridizationError;
use crate::sanitize_support_errors::ConjugationError;
use crate::stereo_order_errors::StereoOrderError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

pub(crate) type AtropSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<ValenceError>,
        Source<
            crate::error_framework::ErrorObject<ConjugationError>,
            Source<
                crate::error_framework::ErrorObject<HybridizationError>,
                Source<
                    crate::error_framework::ErrorObject<RingFindingError>,
                    Source<
                        crate::error_framework::ErrorObject<TopologyValidationError>,
                        Source<
                            crate::error_framework::ErrorObject<CoordinateValidationError>,
                            crate::error_framework::ErrorObject<BondValueError>,
                        >,
                    >,
                >,
            >,
        >,
    >,
>;
pub(crate) type StereoSource = Source<
    AtropSource,
    Source<
        crate::error_framework::ErrorObject<StereoOrderError>,
        Source<
            crate::error_framework::ErrorObject<AtomPropertyError>,
            crate::error_framework::ErrorObject<PropertyIntReadError>,
        >,
    >,
>;

error_projection!(StereoError, StereoErrorKind, ck::StereoError, "stereo", {
    Atropisomer(_), StereoGroup(_), InvalidTopology(_), DuplicateConformerId { id: usize },
    ConformerNotFound { requested: i32 }, ConformerAtomCountMismatch { conformer: usize, rows: usize, atom_count },
    InvalidValence { field, actual: usize, atom_count },
    InvalidValenceValue { field, atom: usize => |value: &ck::AtomId| value.index(), value: i32 },
    ZeroLengthVector { center: usize => |value: &ck::AtomId| value.index(), neighbor: usize => |value: &ck::AtomId| value.index() },
    StereoOrder(_), AtomProperty(_), HydrogenCount(_),
    UnknownStereoProperty { bond: usize => |value: &ck::BondId| value.index(), source },
}, cause = |inner: &ck::StereoError| Python::attach(|py| match inner {
    ck::StereoError::Atropisomer(value) => AtropisomerError::from_rust(value.as_ref().clone()).cause(py),
    ck::StereoError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::StereoError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::StereoError::StereoOrder(value) => Some(StereoOrderError::from_rust(value.clone())),
    ck::StereoError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::StereoError::UnknownStereoProperty { source, .. } => Some(PropertyIntReadError::from_rust(source.clone())),
    ck::StereoError::HydrogenCount(_) | ck::StereoError::DuplicateConformerId { .. }
    | ck::StereoError::ConformerNotFound { .. } | ck::StereoError::ConformerAtomCountMismatch { .. }
    | ck::StereoError::InvalidValence { .. } | ck::StereoError::InvalidValenceValue { .. }
    | ck::StereoError::ZeroLengthVector { .. } => None,
}), source = Option<StereoSource>, methods = {
    #[getter]
    fn field(&self) -> Option<&'static str> { match self.inner {
        ck::StereoError::InvalidValence { field, .. }
        | ck::StereoError::InvalidValenceValue { field, .. } => Some(field), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::StereoError::ConformerAtomCountMismatch { atom_count, .. }
        | ck::StereoError::InvalidValence { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atropisomer(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtropisomerError>>> { match &self.inner {
        ck::StereoError::Atropisomer(value) => AtropisomerError::from_rust(value.as_ref().clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn hydrogen_count(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::StereoError::HydrogenCount(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

error_projection!(AtropisomerError, AtropisomerErrorKind, ck::AtropisomerError, "stereo", {
    StereoGroup(_), SourcePrecondition { message: &'static str },
    CandidateTraversalOutOfRange { bond, bond_count }, CandidateTraversalDuplicate { bond },
    CandidateTraversalUnexpected { bond }, CandidateTraversalMissing { bond },
    Valence(_), Conjugation(_), Hybridization(_), Normalization(_), CarrierCount { atom, count: usize },
    RingFinding(_), InvalidTopology { source }, InvalidCoordinates { source },
    ConformerNotThreeDimensional { conformer: usize }, RingInfoNotSssr,
    RingAtomOutOfRange { atom, atom_count }, RingBondOutOfRange { bond, bond_count },
    RingAtomRowCount { actual, expected }, RingBondRowCount { actual, expected },
    HybridizationAssignmentLength { actual, expected }, BondUpdate { bond, source },
    StereoGroupAtomOutOfRange { atom, atom_count }, StereoGroupBondOutOfRange { bond, bond_count },
    AssignmentBondOutOfRange { bond, bond_count }, AxialBondOutOfRange { bond, bond_count },
    InvalidAtropisomerStereo { bond, stereo: BondStereo => |value: &ck::BondStereo| (*value).into() },
}, cause = |inner: &ck::AtropisomerError| match inner {
    ck::AtropisomerError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::AtropisomerError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::AtropisomerError::Conjugation(value) => Some(ConjugationError::from_rust(value.clone())),
    ck::AtropisomerError::Hybridization(value) => Some(HybridizationError::from_rust(value.clone())),
    ck::AtropisomerError::RingFinding(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::AtropisomerError::InvalidTopology { source } => Some(TopologyValidationError::from_rust(source.clone())),
    ck::AtropisomerError::InvalidCoordinates { source } => Some(CoordinateValidationError::from_rust(source.clone())),
    ck::AtropisomerError::BondUpdate { source, .. } => Some(BondValueError::from_rust(source.clone())),
    ck::AtropisomerError::Normalization(_) | ck::AtropisomerError::SourcePrecondition { .. }
    | ck::AtropisomerError::CandidateTraversalOutOfRange { .. } | ck::AtropisomerError::CandidateTraversalDuplicate { .. }
    | ck::AtropisomerError::CandidateTraversalUnexpected { .. } | ck::AtropisomerError::CandidateTraversalMissing { .. }
    | ck::AtropisomerError::CarrierCount { .. } | ck::AtropisomerError::ConformerNotThreeDimensional { .. }
    | ck::AtropisomerError::RingInfoNotSssr | ck::AtropisomerError::RingAtomOutOfRange { .. }
    | ck::AtropisomerError::RingBondOutOfRange { .. } | ck::AtropisomerError::RingAtomRowCount { .. }
    | ck::AtropisomerError::RingBondRowCount { .. } | ck::AtropisomerError::HybridizationAssignmentLength { .. }
    | ck::AtropisomerError::StereoGroupAtomOutOfRange { .. } | ck::AtropisomerError::StereoGroupBondOutOfRange { .. }
    | ck::AtropisomerError::AssignmentBondOutOfRange { .. } | ck::AtropisomerError::AxialBondOutOfRange { .. }
    | ck::AtropisomerError::InvalidAtropisomerStereo { .. } => None,
}, source = Option<AtropSource>, methods = {
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::CandidateTraversalOutOfRange { bond, .. }
        | ck::AtropisomerError::CandidateTraversalDuplicate { bond }
        | ck::AtropisomerError::CandidateTraversalUnexpected { bond }
        | ck::AtropisomerError::CandidateTraversalMissing { bond }
        | ck::AtropisomerError::RingBondOutOfRange { bond, .. }
        | ck::AtropisomerError::BondUpdate { bond, .. }
        | ck::AtropisomerError::StereoGroupBondOutOfRange { bond, .. }
        | ck::AtropisomerError::AssignmentBondOutOfRange { bond, .. }
        | ck::AtropisomerError::AxialBondOutOfRange { bond, .. }
        | ck::AtropisomerError::InvalidAtropisomerStereo { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::CarrierCount { atom, .. } | ck::AtropisomerError::RingAtomOutOfRange { atom, .. }
        | ck::AtropisomerError::StereoGroupAtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn bond_count(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::CandidateTraversalOutOfRange { bond_count, .. }
        | ck::AtropisomerError::RingBondOutOfRange { bond_count, .. }
        | ck::AtropisomerError::StereoGroupBondOutOfRange { bond_count, .. }
        | ck::AtropisomerError::AssignmentBondOutOfRange { bond_count, .. }
        | ck::AtropisomerError::AxialBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::RingAtomOutOfRange { atom_count, .. }
        | ck::AtropisomerError::StereoGroupAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn actual(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::RingAtomRowCount { actual, .. }
        | ck::AtropisomerError::RingBondRowCount { actual, .. }
        | ck::AtropisomerError::HybridizationAssignmentLength { actual, .. } => Some(actual), _ => None,
    } }
    #[getter]
    fn expected(&self) -> Option<usize> { match self.inner {
        ck::AtropisomerError::RingAtomRowCount { expected, .. }
        | ck::AtropisomerError::RingBondRowCount { expected, .. }
        | ck::AtropisomerError::HybridizationAssignmentLength { expected, .. } => Some(expected), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn normalization(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<StereoError>>> { match &self.inner {
        // The native tuple has no #[source]/#[from]; preserve it as payload.
        ck::AtropisomerError::Normalization(value) => StereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

pub(crate) type LegacySource = Source<
    StereoSource,
    Source<
        crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
        Source<
            crate::error_framework::ErrorObject<MappingValidationError>,
            Source<
                crate::error_framework::ErrorObject<PropertyUIntReadError>,
                Source<
                    crate::error_framework::ErrorObject<CanonicalRankError>,
                    Source<
                        crate::error_framework::ErrorObject<DoubleBondStereoError>,
                        crate::error_framework::ErrorObject<MoleculePropertyError>,
                    >,
                >,
            >,
        >,
    >,
>;

error_projection!(LegacyStereoError, LegacyStereoErrorKind, ck::LegacyStereoError, "stereo", {
    StereoGroup(_), InvalidTopology(_), CipRank(_), DoubleBond(_), StereoOrder(_), Valence(_),
    AtomProperty(_), BondValue(_), PotentialStereo(_), RingFinding(_), Atropisomer(_),
}, cause = |inner: &ck::LegacyStereoError| Python::attach(|py| match inner {
    ck::LegacyStereoError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::LegacyStereoError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::CipRank(value) => CipRankError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::DoubleBond(value) => DoubleBondStereoError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::Valence(value) => ValenceError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::BondValue(value) => BondValueError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::PotentialStereo(value) => PotentialStereoError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::RingFinding(value) => RingFindingError::from_rust(value.clone()).cause(py),
    ck::LegacyStereoError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).cause(py),
}), source = Option<LegacySource>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn topology(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::LegacyStereoError::InvalidTopology(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn cip_rank(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<CipRankError>>> { match &self.inner {
        ck::LegacyStereoError::CipRank(value) => CipRankError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn double_bond(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<DoubleBondStereoError>>> { match &self.inner {
        ck::LegacyStereoError::DoubleBond(value) => DoubleBondStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn stereo_order(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<StereoOrderError>>> { match &self.inner {
        ck::LegacyStereoError::StereoOrder(value) => StereoOrderError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn valence(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<ValenceError>>> { match &self.inner {
        ck::LegacyStereoError::Valence(value) => ValenceError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atom_property(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtomPropertyError>>> { match &self.inner {
        ck::LegacyStereoError::AtomProperty(value) => AtomPropertyError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn bond_value(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<BondValueError>>> { match &self.inner {
        ck::LegacyStereoError::BondValue(value) => BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn potential_stereo(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<PotentialStereoError>>> { match &self.inner {
        ck::LegacyStereoError::PotentialStereo(value) => PotentialStereoError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn ring_finding(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<RingFindingError>>> { match &self.inner {
        ck::LegacyStereoError::RingFinding(value) => RingFindingError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn atropisomer(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<AtropisomerError>>> { match &self.inner {
        ck::LegacyStereoError::Atropisomer(value) => AtropisomerError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stereo_transparency_and_normalization_payload_follow_actual_native_source() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::StereoError::Atropisomer(Box::new(ck::AtropisomerError::InvalidTopology {
                    source: ck::TopologyValidationError::AdjacencyMismatch,
                }));
            assert!(native.source().unwrap().is::<ck::TopologyValidationError>());
            let display = native.to_string();
            let error = StereoError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                error
                    .value(py)
                    .getattr("atropisomer")?
                    .is_instance(<AtropisomerError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            let native = ck::AtropisomerError::Normalization(ck::StereoError::ZeroLengthVector {
                center: ck::AtomId::new(1),
                neighbor: ck::AtomId::new(3),
            });
            assert!(native.source().is_none());
            let error = AtropisomerError::from_rust(native);
            assert!(error.cause(py).is_none());
            let normalization = error.value(py).getattr("normalization")?;
            assert!(normalization.is_instance(<StereoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(normalization.getattr("center")?.extract::<usize>()?, 1);
            assert_eq!(normalization.getattr("neighbor")?.extract::<usize>()?, 3);
            assert!(normalization.setattr("center", 4).is_err());

            let native =
                ck::LegacyStereoError::DoubleBond(ck::DoubleBondStereoError::NumericPropertyRead {
                    atom: None,
                    bond: Some(ck::BondId::new(2)),
                    property: "_UnknownStereo",
                    source: ck::PropertyIntReadError::UnsignedOverflow { value: u32::MAX },
                });
            assert!(native.source().unwrap().is::<ck::PropertyIntReadError>());
            let error = LegacyStereoError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("source")?
                    .is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            let payload = error.value(py).getattr("double_bond")?;
            assert!(payload.is_instance(<DoubleBondStereoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(payload.getattr("bond")?.extract::<usize>()?, 2);
            assert!(payload.getattr("atom")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
