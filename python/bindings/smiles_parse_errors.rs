//! Complete native SMILES parser/writer error values and concrete cause types.
use crate::canonical_errors::{PropertyStringError, RequiredPropertyStringError};
use crate::error_framework::{Source, error_projection};
use crate::hydrogen_errors::HydrogenError;
use crate::kekulize_errors::KekulizeError;
use crate::model_edit_errors::{AtomPropertyError, BondValueError, MoleculePropertyError};
use crate::model_errors::{CoordinateValidationError, StereoGroupError};
use crate::notation_params::CxCoordinateSelection;
use crate::parser_errors::{ParserCarrierError, ParserStereoOrderError};
use crate::potential_stereo_errors::PotentialStereoError;
use crate::property_numeric_errors::{
    PropertyIntReadError, PropertyUIntReadError, PropertyULongReadError,
};
use crate::property_values::PropertyValueError;
use crate::ring_errors::RingFindingError;
use crate::sanitize_errors::SanitizeError;
use crate::smiles_stereo_errors::SmilesStereoError;
use crate::stereo_errors::{AtropisomerError, LegacyStereoError};
use crate::stereo_order_errors::StereoOrderError;
use crate::topology_edit_errors::TopologyEditError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

type ModelSource = Source<
    crate::error_framework::ErrorObject<StereoGroupError>,
    Source<
        crate::error_framework::ErrorObject<TopologyEditError>,
        Source<
            crate::error_framework::ErrorObject<CoordinateValidationError>,
            Source<
                crate::error_framework::ErrorObject<AtomPropertyError>,
                Source<
                    crate::error_framework::ErrorObject<BondValueError>,
                    crate::error_framework::ErrorObject<MoleculePropertyError>,
                >,
            >,
        >,
    >,
>;
type PropertySource = Source<
    crate::error_framework::ErrorObject<PropertyIntReadError>,
    Source<
        crate::error_framework::ErrorObject<PropertyULongReadError>,
        Source<
            crate::error_framework::ErrorObject<PropertyUIntReadError>,
            Source<
                crate::error_framework::ErrorObject<PropertyValueError>,
                Source<
                    crate::error_framework::ErrorObject<RequiredPropertyStringError>,
                    crate::error_framework::ErrorObject<PropertyStringError>,
                >,
            >,
        >,
    >,
>;
type StereoSource = Source<
    crate::error_framework::ErrorObject<LegacyStereoError>,
    Source<
        crate::error_framework::ErrorObject<PotentialStereoError>,
        Source<
            crate::error_framework::ErrorObject<ParserStereoOrderError>,
            Source<
                crate::error_framework::ErrorObject<StereoOrderError>,
                Source<
                    crate::error_framework::ErrorObject<SmilesStereoError>,
                    crate::error_framework::ErrorObject<AtropisomerError>,
                >,
            >,
        >,
    >,
>;
type ChemistrySource = Source<
    crate::error_framework::ErrorObject<RingFindingError>,
    Source<
        crate::error_framework::ErrorObject<ParserCarrierError>,
        Source<
            crate::error_framework::ErrorObject<KekulizeError>,
            Source<
                crate::error_framework::ErrorObject<SanitizeError>,
                crate::error_framework::ErrorObject<HydrogenError>,
            >,
        >,
    >,
>;
pub(crate) type SmilesSource =
    Source<Source<ModelSource, PropertySource>, Source<StereoSource, ChemistrySource>>;

error_projection!(SmilesParseError, SmilesParseErrorKind, ck::SmilesParseError, "smiles", {
    StereoGroup(_), TopologyEdit(_), CxLinkAtomOutOfRange { atom, atom_count },
    CxLinkOrderOutOfRange { index, count }, CxZeroBondOutOfRange { bond, bond_count },
    CxTypedBondOutOfRange { bond, bond_count }, CxAtomPropertyAtomOutOfRange { atom, atom_count },
    CxCoordinateAtomOutOfRange { atom, atom_count }, WriterCanonicalInvariant(_),
    WriterInt(_), WriterULong(_), WriterRings(_), WriterLegacyStereo(_), WriterPotentialStereo(_),
    WriterParserStereoOrder(_), TraversalTooManyOpenRings, TraversalStateIndex { state: &'static str, index, count },
    WriterStereoOrder(_), WriterPropertyKind(_), Coordinates(_), WriterNumeric(_),
    WriterRequiredProperty { atom, source }, WriterPropertyList { atom, source },
    AtomProperty(_), BondProperty(_), MoleculeProperty(_), PropertyKind(_), ParserCarrier(_),
    Unsupported { token: String => |value: &char| value.to_string(), offset },
    Syntax { offset, message: String }, Atom { offset, text: String }, UnclosedRing { index },
    Cx(_), UnsupportedCx(_), UnsupportedWriter(_), CanonicalRank(_), WriterValence(_),
    WriterKekulize(_), WriterProperty(_), WriterStereo(_), ParserSanitize(_),
    ParserRemoveHydrogens(_), ParserStereo(_), ParserAtropisomer(_), WriterStereoBond(_),
    AmbiguousCoordinateSelection { two_d_count: usize, three_d_count: usize }, MissingCoordinateSelection { selection },
    WriterRootAtomOutOfRange { atom_index: usize, atom_count }, Model(_), EmptyReplacementKey, ReplacementCycle { key: String },
}, cause = |inner: &ck::SmilesParseError| match inner {
    ck::SmilesParseError::StereoGroup(value) => Some(StereoGroupError::from_rust(*value)),
    ck::SmilesParseError::TopologyEdit(value) => Some(TopologyEditError::from_rust(value.clone())),
    ck::SmilesParseError::WriterInt(value) => Some(PropertyIntReadError::from_rust(value.clone())),
    ck::SmilesParseError::WriterULong(value) => Some(PropertyULongReadError::from_rust(value.clone())),
    ck::SmilesParseError::WriterRings(value) => Some(RingFindingError::from_rust(value.clone())),
    ck::SmilesParseError::WriterLegacyStereo(value) => Some(LegacyStereoError::from_rust(value.clone())),
    ck::SmilesParseError::WriterPotentialStereo(value) => Some(PotentialStereoError::from_rust(value.clone())),
    ck::SmilesParseError::WriterParserStereoOrder(value) => Some(ParserStereoOrderError::from_rust(value.clone())),
    ck::SmilesParseError::WriterStereoOrder(value) => Some(StereoOrderError::from_rust(value.clone())),
    ck::SmilesParseError::WriterPropertyKind(value) | ck::SmilesParseError::PropertyKind(value) => Some(PropertyValueError::from_rust(value.clone())),
    ck::SmilesParseError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone())),
    ck::SmilesParseError::WriterNumeric(value) => Some(PropertyUIntReadError::from_rust(value.clone())),
    ck::SmilesParseError::WriterRequiredProperty { source, .. } => Some(RequiredPropertyStringError::from_rust(source.clone())),
    ck::SmilesParseError::WriterPropertyList { source, .. } => Some(AtomPropertyError::from_rust(source.clone())),
    ck::SmilesParseError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone())),
    ck::SmilesParseError::BondProperty(value) | ck::SmilesParseError::WriterStereoBond(value) => Some(BondValueError::from_rust(value.clone())),
    ck::SmilesParseError::MoleculeProperty(value) => Some(MoleculePropertyError::from_rust(value.clone())),
    ck::SmilesParseError::ParserCarrier(value) => Some(ParserCarrierError::from_rust(value.clone())),
    ck::SmilesParseError::WriterKekulize(value) => Some(KekulizeError::from_rust(value.clone())),
    ck::SmilesParseError::WriterProperty(value) => Some(PropertyStringError::from_rust(*value)),
    ck::SmilesParseError::ParserSanitize(value) => Some(SanitizeError::from_rust(value.clone())),
    ck::SmilesParseError::ParserRemoveHydrogens(value) => Some(HydrogenError::from_rust(value.clone())),
    ck::SmilesParseError::ParserStereo(value) => Some(SmilesStereoError::from_rust(value.clone())),
    ck::SmilesParseError::ParserAtropisomer(value) => Some(AtropisomerError::from_rust(value.clone())),
    ck::SmilesParseError::CxLinkAtomOutOfRange { .. } | ck::SmilesParseError::CxLinkOrderOutOfRange { .. }
    | ck::SmilesParseError::CxZeroBondOutOfRange { .. } | ck::SmilesParseError::CxTypedBondOutOfRange { .. }
    | ck::SmilesParseError::CxAtomPropertyAtomOutOfRange { .. } | ck::SmilesParseError::CxCoordinateAtomOutOfRange { .. }
    | ck::SmilesParseError::WriterCanonicalInvariant(_) | ck::SmilesParseError::TraversalTooManyOpenRings
    | ck::SmilesParseError::TraversalStateIndex { .. } | ck::SmilesParseError::Unsupported { .. }
    | ck::SmilesParseError::Syntax { .. } | ck::SmilesParseError::Atom { .. } | ck::SmilesParseError::UnclosedRing { .. }
    | ck::SmilesParseError::Cx(_) | ck::SmilesParseError::UnsupportedCx(_) | ck::SmilesParseError::UnsupportedWriter(_)
    | ck::SmilesParseError::CanonicalRank(_) | ck::SmilesParseError::WriterValence(_) | ck::SmilesParseError::WriterStereo(_)
    | ck::SmilesParseError::AmbiguousCoordinateSelection { .. } | ck::SmilesParseError::MissingCoordinateSelection { .. }
    | ck::SmilesParseError::WriterRootAtomOutOfRange { .. } | ck::SmilesParseError::Model(_)
    | ck::SmilesParseError::EmptyReplacementKey | ck::SmilesParseError::ReplacementCycle { .. } => None,
}, source = Option<SmilesSource>, methods = {
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxLinkAtomOutOfRange { atom, .. } => Some(atom),
        ck::SmilesParseError::CxAtomPropertyAtomOutOfRange { atom, .. }
        | ck::SmilesParseError::CxCoordinateAtomOutOfRange { atom, .. }
        | ck::SmilesParseError::WriterRequiredProperty { atom, .. }
        | ck::SmilesParseError::WriterPropertyList { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn atom_count(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxLinkAtomOutOfRange { atom_count, .. }
        | ck::SmilesParseError::CxAtomPropertyAtomOutOfRange { atom_count, .. }
        | ck::SmilesParseError::CxCoordinateAtomOutOfRange { atom_count, .. }
        | ck::SmilesParseError::WriterRootAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None,
    } }
    #[getter]
    fn bond(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxZeroBondOutOfRange { bond, .. }
        | ck::SmilesParseError::CxTypedBondOutOfRange { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[getter]
    fn bond_count(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxZeroBondOutOfRange { bond_count, .. }
        | ck::SmilesParseError::CxTypedBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None,
    } }
    #[getter]
    fn index(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxLinkOrderOutOfRange { index, .. }
        | ck::SmilesParseError::TraversalStateIndex { index, .. } => Some(index),
        ck::SmilesParseError::UnclosedRing { index } => Some(index as usize), _ => None,
    } }
    #[getter]
    fn count(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::CxLinkOrderOutOfRange { count, .. }
        | ck::SmilesParseError::TraversalStateIndex { count, .. } => Some(count), _ => None,
    } }
    #[getter]
    fn offset(&self) -> Option<usize> { match self.inner {
        ck::SmilesParseError::Unsupported { offset, .. } | ck::SmilesParseError::Syntax { offset, .. }
        | ck::SmilesParseError::Atom { offset, .. } => Some(offset), _ => None,
    } }
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn selection(&self) -> Option<CxCoordinateSelection> { match self.inner {
        ck::SmilesParseError::MissingCoordinateSelection { selection } => Some(CxCoordinateSelection { inner: selection }), _ => None,
    } }
    #[getter]
    fn writer_canonical_invariant(&self) -> Option<&'static str> { match self.inner {
        ck::SmilesParseError::WriterCanonicalInvariant(value) => Some(value), _ => None,
    } }
    #[getter]
    fn cx(&self) -> Option<&str> { match &self.inner { ck::SmilesParseError::Cx(value) => Some(value), _ => None } }
    #[getter]
    fn unsupported_cx(&self) -> Option<&'static str> { match self.inner { ck::SmilesParseError::UnsupportedCx(value) => Some(value), _ => None } }
    #[getter]
    fn unsupported_writer(&self) -> Option<&'static str> { match self.inner { ck::SmilesParseError::UnsupportedWriter(value) => Some(value), _ => None } }
    #[getter]
    fn canonical_rank(&self) -> Option<&str> { match &self.inner { ck::SmilesParseError::CanonicalRank(value) => Some(value), _ => None } }
    #[getter]
    fn writer_valence(&self) -> Option<&str> { match &self.inner { ck::SmilesParseError::WriterValence(value) => Some(value), _ => None } }
    #[getter]
    fn writer_stereo(&self) -> Option<&str> { match &self.inner { ck::SmilesParseError::WriterStereo(value) => Some(value), _ => None } }
    #[getter]
    fn model(&self) -> Option<&str> { match &self.inner { ck::SmilesParseError::Model(value) => Some(value), _ => None } }
});

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn native_smiles_syntax_and_writer_errors_keep_structured_context() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let ck::SmilesError::Parse(native) = ck::Molecule::from_smiles("CC(").unwrap_err()
            else {
                panic!("native parser failure");
            };
            let ck::SmilesParseError::Syntax { offset, message } = &native else {
                panic!("native unclosed-branch syntax failure");
            };
            let offset = *offset;
            let message = message.clone();
            let display = native.to_string();
            let error = SmilesParseError::from_rust(native);
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("offset")?.extract::<usize>()?,
                offset
            );
            assert_eq!(
                error.value(py).getattr("message")?.extract::<String>()?,
                message
            );
            assert!(error.value(py).setattr("offset", 0).is_err());
            assert!(error.cause(py).is_none());

            let native = ck::SmilesParseError::WriterNumeric(ck::PropertyUIntReadError::Negative {
                value: -7,
            });
            assert!(native.source().unwrap().is::<ck::PropertyUIntReadError>());
            let error = SmilesParseError::from_rust(native);
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(source.getattr("value")?.extract::<i32>()?, -7);

            let native = ck::SmilesParseError::MissingCoordinateSelection {
                selection: ck::CxCoordinateSelection::ThreeD { id: 7 },
            };
            let error = SmilesParseError::from_rust(native);
            let selection = error.value(py).getattr("selection")?;
            assert!(selection.is_instance_of::<CxCoordinateSelection>());
            assert_eq!(selection.getattr("id")?.extract::<usize>()?, 7);
            assert!(error.cause(py).is_none());
            Ok(())
        })
        .unwrap();
    }
}
