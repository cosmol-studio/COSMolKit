//! Stereo operations preserve native transparent versus explicit source edges.
use crate::{
    atropisomer_errors,
    errors::{enum_error, exception},
    model_errors::{AtomPropertyError, StereoGroupError, TopologyValidationError},
    property_errors::PropertyIntReadError,
    stereo_order_errors::StereoOrderError,
    valence_errors::ValenceError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "StereoGroupError | ValenceError | ConjugationError | HybridizationError | RingFindingError | TopologyValidationError | CoordinateValidationError | BondValueError | StereoOrderError | AtomPropertyError | PropertyIntReadError"
    )]
    #[derive(Clone)]
    pub type StereoCause;
    #[wasm_bindgen(typescript_type = "AtropisomerError")]
    pub type StereoAtropisomer;
    #[wasm_bindgen(typescript_type = "ValenceError")]
    pub type StereoHydrogenCount;
}
pub(crate) fn cause(value: &ck::StereoError) -> Option<StereoCause> {
    match value {
        ck::StereoError::Atropisomer(value) => {
            atropisomer_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::StereoError::StereoGroup(value) => {
            Some(StereoGroupError::from_rust(*value).unchecked_into())
        }
        ck::StereoError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::StereoError::StereoOrder(value) => {
            Some(StereoOrderError::from_rust(value.clone()).unchecked_into())
        }
        ck::StereoError::AtomProperty(value) => {
            Some(AtomPropertyError::from_rust(value.clone()).unchecked_into())
        }
        ck::StereoError::UnknownStereoProperty { source, .. } => {
            Some(PropertyIntReadError::from_rust(source.clone()).unchecked_into())
        }
        ck::StereoError::HydrogenCount(_)
        | ck::StereoError::DuplicateConformerId { .. }
        | ck::StereoError::ConformerNotFound { .. }
        | ck::StereoError::ConformerAtomCountMismatch { .. }
        | ck::StereoError::InvalidValence { .. }
        | ck::StereoError::InvalidValenceValue { .. }
        | ck::StereoError::ZeroLengthVector { .. } => None,
    }
}
enum_error!(StereoError, StereoErrorKind, ck::StereoError, {
    Atropisomer => ck::StereoError::Atropisomer(_), StereoGroup => ck::StereoError::StereoGroup(_), InvalidTopology => ck::StereoError::InvalidTopology(_),
    DuplicateConformerId => ck::StereoError::DuplicateConformerId { .. }, ConformerNotFound => ck::StereoError::ConformerNotFound { .. },
    ConformerAtomCountMismatch => ck::StereoError::ConformerAtomCountMismatch { .. }, InvalidValence => ck::StereoError::InvalidValence { .. },
    InvalidValenceValue => ck::StereoError::InvalidValenceValue { .. }, ZeroLengthVector => ck::StereoError::ZeroLengthVector { .. },
    StereoOrder => ck::StereoError::StereoOrder(_), AtomProperty => ck::StereoError::AtomProperty(_), HydrogenCount => ck::StereoError::HydrogenCount(_),
    UnknownStereoProperty => ck::StereoError::UnknownStereoProperty { .. },
}, {
    #[wasm_bindgen(getter)] pub fn id(&self) -> Option<usize> { match *self.inner { ck::StereoError::DuplicateConformerId { id } => Some(id), _ => None } }
    #[wasm_bindgen(getter)] pub fn requested(&self) -> Option<i32> { match *self.inner { ck::StereoError::ConformerNotFound { requested } => Some(requested), _ => None } }
    #[wasm_bindgen(getter)] pub fn conformer(&self) -> Option<usize> { match *self.inner { ck::StereoError::ConformerAtomCountMismatch { conformer, .. } => Some(conformer), _ => None } }
    #[wasm_bindgen(getter)] pub fn rows(&self) -> Option<usize> { match *self.inner { ck::StereoError::ConformerAtomCountMismatch { rows, .. } => Some(rows), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::StereoError::ConformerAtomCountMismatch { atom_count, .. } | ck::StereoError::InvalidValence { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::StereoError::InvalidValence { field, .. } | ck::StereoError::InvalidValenceValue { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::StereoError::InvalidValence { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::StereoError::InvalidValenceValue { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i32> { match *self.inner { ck::StereoError::InvalidValenceValue { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn center(&self) -> Option<usize> { match *self.inner { ck::StereoError::ZeroLengthVector { center, .. } => Some(center.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn neighbor(&self) -> Option<usize> { match *self.inner { ck::StereoError::ZeroLengthVector { neighbor, .. } => Some(neighbor.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::StereoError::UnknownStereoProperty { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn atropisomer(&self) -> Option<StereoAtropisomer> { match self.inner.as_ref() { ck::StereoError::Atropisomer(value) => Some(atropisomer_errors::AtropisomerError::from_rust(value.as_ref().clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = hydrogenCount)] pub fn hydrogen_count(&self) -> Option<StereoHydrogenCount> { match self.inner.as_ref() { ck::StereoError::HydrogenCount(value) => Some(ValenceError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = StereoCause, cause = cause);
