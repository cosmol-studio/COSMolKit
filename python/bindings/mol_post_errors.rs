//! Native Molfile finalization errors; no parsing or processing in the binding.
use crate::canonical_errors::PropertyStringError;
use crate::double_stereo_errors::DoubleBondStereoError;
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::HydrogenError;
use crate::model_edit_errors::{
    AtomPropertyError, BondValueError, MoleculePropertyError, QueryAtomConversionError,
    QueryStateError,
};
use crate::model_errors::{CoordinateValidationError, TopologyValidationError};
use crate::mol_post_support_errors::{AttachmentExpansionError, BondDirectionStereoError};
use crate::property_numeric_errors::PropertyIntReadError;
use crate::property_values::{PropertyValueError, PropertyValueKind, decode_source_text};
use crate::query_model_errors::QueryGraphError;
use crate::ring_errors::RingFindingError;
use crate::sanitize_errors::SanitizeError;
use crate::stereo_errors::{AtropisomerError, LegacyStereoError, StereoError};
use crate::topology_edit_errors::TopologyEditError;
use crate::valence_errors::ValenceError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type PropertySource = Source<
    crate::error_framework::ErrorObject<PropertyIntReadError>,
    Source<
        crate::error_framework::ErrorObject<PropertyStringError>,
        Source<
            crate::error_framework::ErrorObject<PropertyValueError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                crate::error_framework::ErrorObject<MoleculePropertyError>,
            >,
        >,
    >,
>;
type ModelSource = Source<
    crate::error_framework::ErrorObject<QueryAtomConversionError>,
    Source<
        crate::error_framework::ErrorObject<QueryGraphError>,
        Source<
            crate::error_framework::ErrorObject<QueryStateError>,
            Source<
                crate::error_framework::ErrorObject<TopologyValidationError>,
                Source<
                    crate::error_framework::ErrorObject<TopologyEditError>,
                    Source<
                        crate::error_framework::ErrorObject<BondValueError>,
                        crate::error_framework::ErrorObject<CoordinateValidationError>,
                    >,
                >,
            >,
        >,
    >,
>;
type StereoSource = Source<
    crate::error_framework::ErrorObject<StereoError>,
    Source<
        crate::error_framework::ErrorObject<DoubleBondStereoError>,
        Source<
            crate::error_framework::ErrorObject<BondDirectionStereoError>,
            Source<
                crate::error_framework::ErrorObject<AtropisomerError>,
                crate::error_framework::ErrorObject<LegacyStereoError>,
            >,
        >,
    >,
>;
type ChemistrySource = Source<
    crate::error_framework::ErrorObject<ValenceError>,
    Source<
        crate::error_framework::ErrorObject<RingFindingError>,
        Source<
            crate::error_framework::ErrorObject<SanitizeError>,
            Source<
                crate::error_framework::ErrorObject<HydrogenError>,
                crate::error_framework::ErrorObject<AttachmentExpansionError>,
            >,
        >,
    >,
>;
type ProcessingSource =
    Source<Source<PropertySource, ModelSource>, Source<StereoSource, ChemistrySource>>;

error_projection!(MolProcessingError, MolProcessingErrorKind, ck::MolProcessingError, "io", {
    Valence(_), PropertyInt(_), PropertyString(_), Rings(_), Stereo(_), DoubleBondStereo(_),
    BondDirectionStereo(_), Atropisomer(_), Sanitize(_), Hydrogen(_), LegacyStereo(_),
    AttachmentExpansion(_), QueryAtom(_), QueryGraph(_), QueryState(_), Topology(_),
    TopologyEdit(_), PropertyValue(_), AtomProperty(_), BondValue(_), Coordinates(_), MoleculeProperty(_),
}, cause = |inner: &ck::MolProcessingError| match inner {
    ck::MolProcessingError::Valence(value) => Some(ValenceError::from_rust(value.clone())),
    ck::MolProcessingError::PropertyInt(value) => Some(PropertyIntReadError::from_rust(value.clone())),
    ck::MolProcessingError::PropertyString(value) => Some(PropertyStringError::from_rust(*value)),
    ck::MolProcessingError::Rings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::MolProcessingError::Stereo(value) => Some(StereoError::from_rust(value.clone())),
    ck::MolProcessingError::DoubleBondStereo(value) => Some(DoubleBondStereoError::from_rust(value.clone())),
    ck::MolProcessingError::BondDirectionStereo(value) => Some(BondDirectionStereoError::from_rust(value.clone())),
    ck::MolProcessingError::Atropisomer(value) => Some(AtropisomerError::from_rust(value.clone())),
    ck::MolProcessingError::Sanitize(value) => Some(SanitizeError::from_rust(value.clone())),
    ck::MolProcessingError::Hydrogen(value) => Some(HydrogenError::from_rust(value.clone())),
    ck::MolProcessingError::LegacyStereo(value) => Some(LegacyStereoError::from_rust(value.clone())),
    ck::MolProcessingError::AttachmentExpansion(value) => Some(AttachmentExpansionError::from_rust(value.clone())),
    ck::MolProcessingError::QueryAtom(value) => Some(QueryAtomConversionError::from_rust(value.clone())),
    ck::MolProcessingError::QueryGraph(value) => Some(QueryGraphError::from_rust(value.clone())),
    ck::MolProcessingError::QueryState(value) => Some(QueryStateError::from_rust(value.clone())),
    ck::MolProcessingError::Topology(value) => Some(TopologyValidationError::from_rust(value.clone())),
    ck::MolProcessingError::TopologyEdit(value) => Some(TopologyEditError::from_rust(value.clone())),
    ck::MolProcessingError::PropertyValue(value) => Some(PropertyValueError::from_rust(value.clone())),
    ck::MolProcessingError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::MolProcessingError::BondValue(value) => Some(BondValueError::from_rust(value.clone())),
    ck::MolProcessingError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::MolProcessingError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
}, source = Option<ProcessingSource>);

#[derive(IntoPyObject)]
pub enum PostValue {
    Unsigned(u32),
    Text(String),
}

error_projection!(MolPostError, MolPostErrorKind, ck::MolPostError, "io", {
    MissingCapability(_), UnsignedPropertyOverflow { atom, property, value },
    InvalidPropertyKind { atom, property, kind }, AttachmentValue { atom, value },
    SGroupPropertyKind { group: usize => |value: &ck::SubstanceGroupId| value.index(), property, kind },
    AttachmentExpansion(_), Representation(_), DataFieldNumber { field: &'static str, value },
    QueryAtomConversion(_), Processing(_),
}, cause = |inner: &ck::MolPostError| match inner {
    ck::MolPostError::AttachmentExpansion(value) | ck::MolPostError::Processing(value) => Some(MolProcessingError::from_rust(value.clone())),
    ck::MolPostError::QueryAtomConversion(_) | ck::MolPostError::MissingCapability(_)
    | ck::MolPostError::UnsignedPropertyOverflow { .. } | ck::MolPostError::InvalidPropertyKind { .. }
    | ck::MolPostError::AttachmentValue { .. } | ck::MolPostError::SGroupPropertyKind { .. }
    | ck::MolPostError::Representation(_) | ck::MolPostError::DataFieldNumber { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<MolProcessingError>>, methods = {
    #[getter]
    #[doc = "Capability required by the failing native operation. Applies to MolPostError.MissingCapability. None when this payload is absent."] fn capability(&self) -> Option<&'static str> { match self.inner { ck::MolPostError::MissingCapability(value) => Some(value), _ => None } }
    #[getter]
    #[doc = "Representation the native operation could not support. Applies to MolPostError.Representation. None when this payload is absent."] fn representation(&self) -> Option<&'static str> { match self.inner { ck::MolPostError::Representation(value) => Some(value), _ => None } }
    #[getter]
    #[doc = "Atom identifier implicated in the failure. Applies to MolPostError.UnsignedPropertyOverflow, MolPostError.InvalidPropertyKind, MolPostError.AttachmentValue. None when this payload is absent."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { atom, .. } | ck::MolPostError::InvalidPropertyKind { atom, .. }
        | ck::MolPostError::AttachmentValue { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    #[doc = "Property name involved in the failure. Applies to MolPostError.UnsignedPropertyOverflow, MolPostError.InvalidPropertyKind, MolPostError.SGroupPropertyKind. None when this payload is absent."] fn property(&self) -> Option<&'static str> { match self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { property, .. } | ck::MolPostError::InvalidPropertyKind { property, .. }
        | ck::MolPostError::SGroupPropertyKind { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Native failure-specific category or value kind. Applies to MolPostError.InvalidPropertyKind, MolPostError.SGroupPropertyKind. None when this payload is absent."] fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::MolPostError::InvalidPropertyKind { kind, .. } | ck::MolPostError::SGroupPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    #[doc = "Value rejected or reported by the native operation. Applies to MolPostError.UnsignedPropertyOverflow, MolPostError.AttachmentValue, MolPostError.DataFieldNumber. None when this payload is absent."] fn value(&self, py: Python<'_>) -> PyResult<Option<PostValue>> { match &self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { value, .. } => Ok(Some(PostValue::Unsigned(*value))),
        ck::MolPostError::AttachmentValue { value, .. } | ck::MolPostError::DataFieldNumber { value, .. } =>
            decode_source_text(py, value).map(|value| Some(PostValue::Text(value))),
        _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    #[doc = "Structured QueryAtomConversionError payload of MolPostError.QueryAtomConversion; None for other variants. Inspect __cause__ separately for the actual exception chain."] fn query_atom_conversion(&self, py: Python<'_>) -> PyResult<Option<crate::error_framework::ErrorObject<QueryAtomConversionError>>> { match &self.inner {
        ck::MolPostError::QueryAtomConversion(value) => QueryAtomConversionError::from_rust(value.clone()).value(py).extract().map(Some).map_err(Into::into),
        _ => Ok(None),
    } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_processing_and_post_errors_keep_explicit_levels_and_typed_fields() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native =
                ck::MolPostError::AttachmentExpansion(ck::MolProcessingError::AttachmentExpansion(
                    ck::AttachmentExpansionError::Placement(ck::HydrogenError::InvalidTopology(
                        ck::TopologyValidationError::AdjacencyMismatch,
                    )),
                ));
            assert!(native.source().unwrap().is::<ck::MolProcessingError>());
            let display = native.to_string();
            let error = MolPostError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            let processing = error.value(py).getattr("source")?;
            assert!(processing.is_instance(<MolProcessingError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(processing.is(&error.value(py).getattr("__cause__")?));
            let attachment = processing.getattr("source")?;
            assert!(attachment.is_instance(<AttachmentExpansionError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            let hydrogen = attachment.getattr("source")?;
            assert!(hydrogen.is_instance(<HydrogenError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(hydrogen.getattr("source")?.is_none());

            let native = ck::MolProcessingError::BondDirectionStereo(
                ck::BondDirectionStereoError::Valence(ck::ValenceError::InvalidTopology {
                    source: ck::TopologyValidationError::AdjacencyMismatch,
                }),
            );
            assert!(
                native
                    .source()
                    .unwrap()
                    .is::<ck::BondDirectionStereoError>()
            );
            let error = MolProcessingError::from_rust(native);
            let direction = error.value(py).getattr("source")?;
            assert!(direction.is_instance(<BondDirectionStereoError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                direction
                    .getattr("source")?
                    .is_instance(<TopologyValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert!(
                direction
                    .getattr("valence")?
                    .is_instance(<ValenceError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );

            for native in [
                ck::MolPostError::UnsignedPropertyOverflow {
                    atom: ck::AtomId::new(usize::MAX),
                    property: "key\0name",
                    value: u32::MAX,
                },
                ck::MolPostError::InvalidPropertyKind {
                    atom: ck::AtomId::new(usize::MAX),
                    property: "key\0name",
                    kind: ck::PropertyValueKind::Bool,
                },
                ck::MolPostError::AttachmentValue {
                    atom: ck::AtomId::new(usize::MAX),
                    value: "value\0text".into(),
                },
                ck::MolPostError::SGroupPropertyKind {
                    group: ck::SubstanceGroupId::new(usize::MAX),
                    property: "key\0name",
                    kind: ck::PropertyValueKind::IntVector,
                },
                ck::MolPostError::DataFieldNumber {
                    field: "field\0name",
                    value: "value\0text".into(),
                },
                ck::MolPostError::QueryAtomConversion(
                    ck::QueryAtomConversionError::NonElementAtomicNumber {
                        atom: ck::AtomId::new(usize::MAX),
                        atomic_number: 255,
                    },
                ),
                ck::MolPostError::MissingCapability("capability\0name"),
                ck::MolPostError::Representation("representation\0detail"),
            ] {
                assert!(native.source().is_none());
                let display = native.to_string();
                let error = MolPostError::from_rust(native);
                let value = error.value(py);
                assert_eq!(value.str()?.to_str()?, display);
                assert!(value.getattr("source")?.is_none());
                for field in [
                    "atom",
                    "property",
                    "kind",
                    "group",
                    "field",
                    "value",
                    "query_atom_conversion",
                    "capability",
                    "representation",
                ] {
                    let payload = value.getattr(field)?;
                    if !payload.is_none() {
                        assert!(value.setattr(field, payload).is_err());
                    }
                }
                match value.getattr("variant")?.extract::<MolPostErrorKind>()? {
                    MolPostErrorKind::UnsignedPropertyOverflow => {
                        assert_eq!(value.getattr("atom")?.extract::<usize>()?, usize::MAX);
                        assert_eq!(value.getattr("value")?.extract::<u32>()?, u32::MAX);
                        assert_eq!(value.getattr("property")?.extract::<String>()?, "key\0name");
                    }
                    MolPostErrorKind::InvalidPropertyKind => assert!(
                        value.getattr("kind")?.extract::<PropertyValueKind>()?
                            == PropertyValueKind::Bool
                    ),
                    MolPostErrorKind::SGroupPropertyKind => {
                        assert_eq!(value.getattr("group")?.extract::<usize>()?, usize::MAX)
                    }
                    MolPostErrorKind::AttachmentValue | MolPostErrorKind::DataFieldNumber => {
                        assert_eq!(value.getattr("value")?.extract::<String>()?, "value\0text")
                    }
                    MolPostErrorKind::QueryAtomConversion => assert!(
                        value
                            .getattr("query_atom_conversion")?
                            .is_instance(<QueryAtomConversionError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
                    ),
                    MolPostErrorKind::MissingCapability => assert_eq!(
                        value.getattr("capability")?.extract::<String>()?,
                        "capability\0name"
                    ),
                    MolPostErrorKind::Representation => assert_eq!(
                        value.getattr("representation")?.extract::<String>()?,
                        "representation\0detail"
                    ),
                    _ => unreachable!(),
                }
            }
            Ok(())
        })
        .unwrap();
    }
}
