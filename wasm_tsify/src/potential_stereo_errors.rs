//! Potential-stereo diagnostics expose the canonical source and payload types.
use crate::{
    canonical_errors::CanonicalRankError,
    double_stereo_errors::DoubleBondStereoError,
    enums::BondOrder,
    errors::{enum_error, exception},
    model_errors::{BondValueError, PropertyValueKind, TopologyValidationError},
    ring_errors::RingFindingError,
    stereo_order_errors::StereoOrderError,
    valence_errors::ValenceError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "RingFindingError | CanonicalRankError | TopologyValidationError | ValenceError | StereoOrderError | DoubleBondStereoError | BondValueError"
    )]
    #[derive(Clone)]
    pub type PotentialStereoCause;
}
pub(crate) fn cause(value: &ck::PotentialStereoError) -> Option<PotentialStereoCause> {
    match value {
        ck::PotentialStereoError::Ring(value) => {
            Some(RingFindingError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::Rank(value) => {
            Some(CanonicalRankError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::Valence(value) => {
            Some(ValenceError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::StereoOrder(value) => {
            Some(StereoOrderError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::DoubleStereo(value) => {
            Some(DoubleBondStereoError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::BondValue(value) => {
            Some(BondValueError::from_rust(value.clone()).unchecked_into())
        }
        ck::PotentialStereoError::InvalidPropertyKind { .. }
        | ck::PotentialStereoError::InvalidRingStereoReference { .. }
        | ck::PotentialStereoError::EmptyRingStereoReferences { .. }
        | ck::PotentialStereoError::InvalidValence { .. }
        | ck::PotentialStereoError::InvalidValenceValue { .. }
        | ck::PotentialStereoError::AtomOutOfRange { .. }
        | ck::PotentialStereoError::InvalidRingInfo { .. }
        | ck::PotentialStereoError::InvalidAtomDegree { .. }
        | ck::PotentialStereoError::InvalidBondDegree { .. }
        | ck::PotentialStereoError::InvalidStereoReferences { .. }
        | ck::PotentialStereoError::UnsupportedBondOrder { .. }
        | ck::PotentialStereoError::InvalidChiralPermutation { .. }
        | ck::PotentialStereoError::AtropisomerDependencyUnavailable { .. }
        | ck::PotentialStereoError::RefinementDidNotConverge { .. } => None,
    }
}
enum_error!(PotentialStereoError, PotentialStereoErrorKind, ck::PotentialStereoError, {
    InvalidPropertyKind => ck::PotentialStereoError::InvalidPropertyKind { .. }, InvalidRingStereoReference => ck::PotentialStereoError::InvalidRingStereoReference { .. },
    EmptyRingStereoReferences => ck::PotentialStereoError::EmptyRingStereoReferences { .. }, Ring => ck::PotentialStereoError::Ring(_), Rank => ck::PotentialStereoError::Rank(_),
    InvalidTopology => ck::PotentialStereoError::InvalidTopology(_), InvalidValence => ck::PotentialStereoError::InvalidValence { .. },
    InvalidValenceValue => ck::PotentialStereoError::InvalidValenceValue { .. }, AtomOutOfRange => ck::PotentialStereoError::AtomOutOfRange { .. },
    InvalidRingInfo => ck::PotentialStereoError::InvalidRingInfo { .. }, InvalidAtomDegree => ck::PotentialStereoError::InvalidAtomDegree { .. },
    InvalidBondDegree => ck::PotentialStereoError::InvalidBondDegree { .. }, InvalidStereoReferences => ck::PotentialStereoError::InvalidStereoReferences { .. },
    UnsupportedBondOrder => ck::PotentialStereoError::UnsupportedBondOrder { .. }, InvalidChiralPermutation => ck::PotentialStereoError::InvalidChiralPermutation { .. },
    AtropisomerDependencyUnavailable => ck::PotentialStereoError::AtropisomerDependencyUnavailable { .. }, RefinementDidNotConverge => ck::PotentialStereoError::RefinementDidNotConverge { .. },
    Valence => ck::PotentialStereoError::Valence(_), StereoOrder => ck::PotentialStereoError::StereoOrder(_), DoubleStereo => ck::PotentialStereoError::DoubleStereo(_), BondValue => ck::PotentialStereoError::BondValue(_),
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::PotentialStereoError::InvalidPropertyKind { atom, .. } | ck::PotentialStereoError::InvalidRingStereoReference { atom, .. } | ck::PotentialStereoError::EmptyRingStereoReferences { atom }
        | ck::PotentialStereoError::InvalidValenceValue { atom, .. } | ck::PotentialStereoError::AtomOutOfRange { atom, .. } | ck::PotentialStereoError::InvalidAtomDegree { atom, .. }
        | ck::PotentialStereoError::InvalidChiralPermutation { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::PotentialStereoError::InvalidBondDegree { bond, .. } | ck::PotentialStereoError::InvalidStereoReferences { bond, .. }
        | ck::PotentialStereoError::UnsupportedBondOrder { bond, .. } | ck::PotentialStereoError::AtropisomerDependencyUnavailable { bond } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn property(&self) -> Option<String> { match *self.inner { ck::PotentialStereoError::InvalidPropertyKind { property, .. } => Some(property.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::PotentialStereoError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<f64> { match *self.inner {
        ck::PotentialStereoError::InvalidRingStereoReference { value, .. } | ck::PotentialStereoError::InvalidValenceValue { value, .. } => Some(f64::from(value)),
        ck::PotentialStereoError::InvalidRingInfo { value, .. } => Some(value as f64), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::PotentialStereoError::InvalidValence { field, .. } | ck::PotentialStereoError::InvalidValenceValue { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::InvalidValence { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::InvalidValence { atom_count, .. } | ck::PotentialStereoError::AtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn reason(&self) -> Option<String> { match *self.inner { ck::PotentialStereoError::InvalidRingInfo { reason, .. } | ck::PotentialStereoError::InvalidStereoReferences { reason, .. } => Some(reason.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn row(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::InvalidRingInfo { row, .. } => Some(row), _ => None } }
    #[wasm_bindgen(getter)] pub fn limit(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::InvalidRingInfo { limit, .. } => Some(limit), _ => None } }
    #[wasm_bindgen(getter)] pub fn degree(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::InvalidAtomDegree { degree, .. } | ck::PotentialStereoError::InvalidBondDegree { degree, .. } => Some(degree), _ => None } }
    #[wasm_bindgen(getter)] pub fn endpoint(&self) -> Option<String> { match *self.inner { ck::PotentialStereoError::InvalidBondDegree { endpoint, .. } => Some(endpoint.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::PotentialStereoError::UnsupportedBondOrder { order, .. } => Some(order.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn permutation(&self) -> Option<u32> { match *self.inner { ck::PotentialStereoError::InvalidChiralPermutation { permutation, .. } => Some(permutation), _ => None } }
    #[wasm_bindgen(getter)] pub fn iterations(&self) -> Option<usize> { match *self.inner { ck::PotentialStereoError::RefinementDidNotConverge { iterations } => Some(iterations), _ => None } }
}, source = PotentialStereoCause, cause = cause);
