//! Detached edit failures retain actual causes and source-free typed payloads.
use crate::attachment_values::TemplateAttachmentOrderError;
use crate::chemical_values::enum_value;
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MappingValidationError, MoleculePropertyError,
};
use crate::model_errors::{CoordinateValidationError, StereoGroupError, TopologyValidationError};
use crate::property_values::decode_source_text;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

enum_value!(BondEndPointsParseErrorKind, { InvalidArgument, OutOfRange, UnrepresentableEndpoint });
pub(crate) type EditSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<CoordinateValidationError>,
        Source<
            crate::error_framework::ErrorObject<MoleculePropertyError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                crate::error_framework::ErrorObject<TemplateAttachmentOrderError>,
            >,
        >,
    >,
>;

error_projection!(TopologyEditError, TopologyEditErrorKind, ck::TopologyEditError, "topology_edit", {
    StereoGroup(_), Coordinate(_), MoleculeProperty(_), AtomProperty(_),
    SourceGroupParentCycle { group: usize => |value: &ck::SubstanceGroupId| value.index() },
    InvalidSource(_), AtomOutOfRange { atom, atom_count }, BondOutOfRange { bond, bond_count: usize },
    BondEndPointsParse { bond, token_index: usize, token, kind },
    IncompleteSourceBatchMasks, NullBondAtom,
    DuplicateBond { begin: usize => |value: &ck::AtomId| value.index(), end: usize => |value: &ck::AtomId| value.index() },
    InvalidBond(_), PermutationLength { actual: usize, expected: usize },
    PermutationAtomOutOfRange { position, atom, atom_count }, PermutationDuplicateAtom { position, atom },
    InvalidMapping(_), TemplateAttachmentRemap { carrier: usize => |value: &ck::AtomId| value.index(), source }, InvalidResult(_),
}, cause = |inner: &ck::TopologyEditError| match inner {
    ck::TopologyEditError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::TopologyEditError::Coordinate(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::TopologyEditError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::TopologyEditError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::TopologyEditError::TemplateAttachmentRemap { source, .. } => Some(TemplateAttachmentOrderError::from_rust(source.clone())),
    ck::TopologyEditError::SourceGroupParentCycle { .. } | ck::TopologyEditError::InvalidSource(_)
    | ck::TopologyEditError::AtomOutOfRange { .. } | ck::TopologyEditError::BondOutOfRange { .. }
    | ck::TopologyEditError::BondEndPointsParse { .. } | ck::TopologyEditError::IncompleteSourceBatchMasks
    | ck::TopologyEditError::NullBondAtom | ck::TopologyEditError::DuplicateBond { .. }
    | ck::TopologyEditError::InvalidBond(_) | ck::TopologyEditError::PermutationLength { .. }
    | ck::TopologyEditError::PermutationAtomOutOfRange { .. } | ck::TopologyEditError::PermutationDuplicateAtom { .. }
    | ck::TopologyEditError::InvalidMapping(_) | ck::TopologyEditError::InvalidResult(_) => None,
}, source = Option<EditSource>, methods = {
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to TopologyEditError.AtomOutOfRange, TopologyEditError.PermutationAtomOutOfRange, TopologyEditError.PermutationDuplicateAtom. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::TopologyEditError::AtomOutOfRange { atom, .. }
        | ck::TopologyEditError::PermutationAtomOutOfRange { atom, .. }
        | ck::TopologyEditError::PermutationDuplicateAtom { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    #[doc = "Atom count used when validating the reported atom reference. Applies to TopologyEditError.AtomOutOfRange, TopologyEditError.PermutationAtomOutOfRange. None when this payload is absent."] fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::TopologyEditError::AtomOutOfRange { atom_count, .. }
        | ck::TopologyEditError::PermutationAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    #[doc = "Bond identifier implicated in the failure. Applies to TopologyEditError.BondOutOfRange, TopologyEditError.BondEndPointsParse. None when this payload is absent."] fn bond(&self) -> Option<usize> { match self.inner {
        ck::TopologyEditError::BondOutOfRange { bond, .. }
        | ck::TopologyEditError::BondEndPointsParse { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    #[doc = "Position reported by the native failure. Applies to TopologyEditError.PermutationAtomOutOfRange, TopologyEditError.PermutationDuplicateAtom. None when this payload is absent."] fn position(&self) -> Option<usize> { match self.inner {
        ck::TopologyEditError::PermutationAtomOutOfRange { position, .. }
        | ck::TopologyEditError::PermutationDuplicateAtom { position, .. } => Some(position), _ => None,
    } }
    #[getter]
    #[doc = "Input token rejected by the parser. Applies to TopologyEditError.BondEndPointsParse. None when this payload is absent."] fn token(&self, py: Python<'_>) -> PyResult<Option<String>> { match &self.inner {
        ck::TopologyEditError::BondEndPointsParse { token, .. } => decode_source_text(py, token).map(Some), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to TopologyEditError.BondEndPointsParse. None when this payload is absent."] fn kind(&self) -> Option<BondEndPointsParseErrorKind> { match self.inner {
        ck::TopologyEditError::BondEndPointsParse { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of TopologyEditError.InvalidSource; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn invalid_source(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::TopologyEditError::InvalidSource(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured TopologyValidationError payload of TopologyEditError.InvalidResult; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn invalid_result(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<TopologyValidationError>>> { match &self.inner {
        ck::TopologyEditError::InvalidResult(value) => TopologyValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured BondValueError payload of TopologyEditError.InvalidBond; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn invalid_bond(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<BondValueError>>> { match &self.inner {
        ck::TopologyEditError::InvalidBond(value) => BondValueError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured MappingValidationError payload of TopologyEditError.InvalidMapping; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn invalid_mapping(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<MappingValidationError>>> { match &self.inner {
        ck::TopologyEditError::InvalidMapping(value) => MappingValidationError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into), _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn topology_edit_error_keeps_actual_kind_separate_from_variant_and_cause() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::TopologyEditError::BondEndPointsParse {
                bond: ck::BondId::new(2),
                token_index: 3,
                token: "x.2".into(),
                kind: ck::BondEndPointsParseErrorKind::InvalidArgument,
            };
            let display = native.to_string();
            assert!(native.source().is_none());
            let error = TopologyEditError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<TopologyEditErrorKind>()?
                    == TopologyEditErrorKind::BondEndPointsParse
            );
            assert!(
                error
                    .value(py)
                    .getattr("kind")?
                    .extract::<BondEndPointsParseErrorKind>()?
                    == BondEndPointsParseErrorKind::InvalidArgument
            );
            assert_eq!(error.value(py).getattr("bond")?.extract::<usize>()?, 2);
            assert_eq!(
                error.value(py).getattr("token_index")?.extract::<usize>()?,
                3
            );
            assert_eq!(
                error.value(py).getattr("token")?.extract::<String>()?,
                "x.2"
            );
            assert!(error.value(py).setattr("token_index", 0).is_err());
            assert!(error.cause(py).is_none());

            let native = ck::TopologyEditError::InvalidSource(
                ck::TopologyValidationError::AdjacencyMismatch,
            );
            assert!(native.source().is_none());
            let error = TopologyEditError::from_rust(native);
            assert!(
                error
                    .value(py)
                    .getattr("invalid_source")?
                    .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(error.cause(py).is_none());
            let native = ck::TopologyEditError::AtomProperty(ck::AtomPropertyError::EmptyKey);
            assert!(native.source().unwrap().is::<ck::AtomPropertyError>());
            let error = TopologyEditError::from_rust(native);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is_instance(<AtomPropertyError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            Ok(())
        })
        .unwrap();
    }
}
