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
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;

type PropertySource = Source<
    Py<PropertyIntReadError>,
    Source<
        Py<PropertyStringError>,
        Source<Py<PropertyValueError>, Source<Py<AtomPropertyError>, Py<MoleculePropertyError>>>,
    >,
>;
type ModelSource = Source<
    Py<QueryAtomConversionError>,
    Source<
        Py<QueryGraphError>,
        Source<
            Py<QueryStateError>,
            Source<
                Py<TopologyValidationError>,
                Source<
                    Py<TopologyEditError>,
                    Source<Py<BondValueError>, Py<CoordinateValidationError>>,
                >,
            >,
        >,
    >,
>;
type StereoSource = Source<
    Py<StereoError>,
    Source<
        Py<DoubleBondStereoError>,
        Source<Py<BondDirectionStereoError>, Source<Py<AtropisomerError>, Py<LegacyStereoError>>>,
    >,
>;
type ChemistrySource = Source<
    Py<ValenceError>,
    Source<
        Py<RingFindingError>,
        Source<Py<SanitizeError>, Source<Py<HydrogenError>, Py<AttachmentExpansionError>>>,
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
}, source = Option<Py<MolProcessingError>>, methods = {
    #[getter]
    fn capability(&self) -> Option<&'static str> { match self.inner { ck::MolPostError::MissingCapability(value) => Some(value), _ => None } }
    #[getter]
    fn representation(&self) -> Option<&'static str> { match self.inner { ck::MolPostError::Representation(value) => Some(value), _ => None } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { atom, .. } | ck::MolPostError::InvalidPropertyKind { atom, .. }
        | ck::MolPostError::AttachmentValue { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { property, .. } | ck::MolPostError::InvalidPropertyKind { property, .. }
        | ck::MolPostError::SGroupPropertyKind { property, .. } => Some(property), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn kind(&self) -> Option<PropertyValueKind> { match self.inner {
        ck::MolPostError::InvalidPropertyKind { kind, .. } | ck::MolPostError::SGroupPropertyKind { kind, .. } => Some(kind.into()), _ => None,
    } }
    #[getter]
    fn value(&self, py: Python<'_>) -> PyResult<Option<PostValue>> { match &self.inner {
        ck::MolPostError::UnsignedPropertyOverflow { value, .. } => Ok(Some(PostValue::Unsigned(*value))),
        ck::MolPostError::AttachmentValue { value, .. } | ck::MolPostError::DataFieldNumber { value, .. } =>
            decode_source_text(py, value).map(|value| Some(PostValue::Text(value))),
        _ => Ok(None),
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn query_atom_conversion(&self, py: Python<'_>) -> PyResult<Option<Py<QueryAtomConversionError>>> { match &self.inner {
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
            assert!(processing.is_instance_of::<MolProcessingError>());
            assert!(processing.is(&error.value(py).getattr("__cause__")?));
            let attachment = processing.getattr("source")?;
            assert!(attachment.is_instance_of::<AttachmentExpansionError>());
            let hydrogen = attachment.getattr("source")?;
            assert!(hydrogen.is_instance_of::<HydrogenError>());
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
            assert!(direction.is_instance_of::<BondDirectionStereoError>());
            assert!(
                direction
                    .getattr("source")?
                    .is_instance_of::<TopologyValidationError>()
            );
            assert!(
                direction
                    .getattr("valence")?
                    .is_instance_of::<ValenceError>()
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
                            .is_instance_of::<QueryAtomConversionError>()
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
