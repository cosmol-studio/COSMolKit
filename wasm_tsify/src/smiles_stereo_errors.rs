//! Source-defined SMILES stereo stage wrappers.
use crate::{
    double_stereo_errors,
    errors::{enum_error, exception},
    legacy_stereo_errors,
    model_errors::{CoordinateValidationError, MoleculePropertyError},
    ring_errors, valence_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError | StereoOrderError | AtomPropertyError | PropertyIntReadError | TemplateAttachmentOrderError | MappingValidationError | PropertyUIntReadError | CanonicalRankError | DoubleBondStereoError | MoleculePropertyError"
    )]
    #[derive(Clone)]
    pub type SmilesStereoCause;
    #[wasm_bindgen(typescript_type = "MoleculePropertyError")]
    pub type SmilesStereoProperties;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type SmilesStereoCoordinates;
    #[wasm_bindgen(typescript_type = "DoubleBondStereoError")]
    pub type SmilesStereoDirections;
    #[wasm_bindgen(typescript_type = "LegacyStereoError")]
    pub type SmilesStereoAssignment;
    #[wasm_bindgen(typescript_type = "RingFindingError")]
    pub type SmilesStereoRings;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type SmilesStereoValence;
}
pub(crate) fn cause(value: &ck::SmilesStereoError) -> Option<SmilesStereoCause> {
    match value {
        ck::SmilesStereoError::Properties(_) | ck::SmilesStereoError::Coordinates(_) => None,
        ck::SmilesStereoError::Directions(value) => {
            double_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::SmilesStereoError::Assignment(value) => {
            legacy_stereo_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::SmilesStereoError::Rings(value) => {
            ring_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::SmilesStereoError::Valence(value) => {
            valence_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::SmilesStereoError::SelectedConformerMissing { .. }
        | ck::SmilesStereoError::ValenceRows { .. }
        | ck::SmilesStereoError::RingRows { .. } => None,
    }
}
enum_error!(SmilesStereoError, SmilesStereoErrorKind, ck::SmilesStereoError, {
    Properties => ck::SmilesStereoError::Properties(_), SelectedConformerMissing => ck::SmilesStereoError::SelectedConformerMissing { .. },
    ValenceRows => ck::SmilesStereoError::ValenceRows { .. }, RingRows => ck::SmilesStereoError::RingRows { .. }, Coordinates => ck::SmilesStereoError::Coordinates(_),
    Directions => ck::SmilesStereoError::Directions(_), Assignment => ck::SmilesStereoError::Assignment(_), Rings => ck::SmilesStereoError::Rings(_), Valence => ck::SmilesStereoError::Valence(_),
}, {
    #[wasm_bindgen(getter)] pub fn id(&self) -> Option<usize> { match *self.inner { ck::SmilesStereoError::SelectedConformerMissing { id } => Some(id), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::SmilesStereoError::ValenceRows { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn dimension(&self) -> Option<String> { match *self.inner { ck::SmilesStereoError::RingRows { dimension, .. } => Some(dimension.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::SmilesStereoError::ValenceRows { actual, .. } | ck::SmilesStereoError::RingRows { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::SmilesStereoError::ValenceRows { expected, .. } | ck::SmilesStereoError::RingRows { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn properties(&self) -> Option<SmilesStereoProperties> { match self.inner.as_ref() { ck::SmilesStereoError::Properties(value) => Some(MoleculePropertyError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn coordinates(&self) -> Option<SmilesStereoCoordinates> { match self.inner.as_ref() { ck::SmilesStereoError::Coordinates(value) => Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn directions(&self) -> Option<SmilesStereoDirections> { match self.inner.as_ref() { ck::SmilesStereoError::Directions(value) => Some(double_stereo_errors::DoubleBondStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn assignment(&self) -> Option<SmilesStereoAssignment> { match self.inner.as_ref() { ck::SmilesStereoError::Assignment(value) => Some(legacy_stereo_errors::LegacyStereoError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn rings(&self) -> Option<SmilesStereoRings> { match self.inner.as_ref() { ck::SmilesStereoError::Rings(value) => Some(ring_errors::RingFindingError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn valence(&self) -> Option<SmilesStereoValence> { match self.inner.as_ref() { ck::SmilesStereoError::Valence(value) => Some(valence_errors::ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = SmilesStereoCause, cause = cause);
