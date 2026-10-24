//! Double-bond stereo errors expose native context without reinterpretation.
use crate::{
    enums::{BondDirection, BondOrder, BondStereo},
    errors::{enum_error, exception},
    model_errors::{
        BondValueError, CoordinateValidationError, PropertyValueKind, TopologyValidationError,
    },
    property_errors::PropertyIntReadError,
    ring_errors::RingFindingError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(
        typescript_type = "RingFindingError | PropertyIntReadError | TopologyValidationError | BondValueError | CoordinateValidationError"
    )]
    #[derive(Clone)]
    pub type DoubleStereoCause;
}
pub(crate) fn cause(value: &ck::DoubleBondStereoError) -> Option<DoubleStereoCause> {
    match value {
        ck::DoubleBondStereoError::Rings(value) => {
            Some(RingFindingError::from_rust(value.clone()).unchecked_into())
        }
        ck::DoubleBondStereoError::NumericPropertyRead { source, .. } => {
            Some(PropertyIntReadError::from_rust(source.clone()).unchecked_into())
        }
        ck::DoubleBondStereoError::InvalidTopology(value) => {
            Some(TopologyValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::DoubleBondStereoError::BondValue(value) => {
            Some(BondValueError::from_rust(value.clone()).unchecked_into())
        }
        ck::DoubleBondStereoError::InvalidConformer(value) => {
            Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into())
        }
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { .. }
        | ck::DoubleBondStereoError::InvalidPropertyKind { .. }
        | ck::DoubleBondStereoError::BondOutOfRange { .. }
        | ck::DoubleBondStereoError::AtomOutOfRange { .. }
        | ck::DoubleBondStereoError::NotDoubleBond { .. }
        | ck::DoubleBondStereoError::UndefinedStereo { .. }
        | ck::DoubleBondStereoError::UnsupportedStereo { .. }
        | ck::DoubleBondStereoError::InvalidEndpointDegree { .. }
        | ck::DoubleBondStereoError::InvalidDirection { .. }
        | ck::DoubleBondStereoError::RankCount { .. }
        | ck::DoubleBondStereoError::RingInfoNotInitialized
        | ck::DoubleBondStereoError::RingRowCount { .. }
        | ck::DoubleBondStereoError::StereoReferenceRequired { .. }
        | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { .. }
        | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { .. }
        | ck::DoubleBondStereoError::ControllingAtomMismatch { .. }
        | ck::DoubleBondStereoError::ControllingSlotOverflow { .. } => None,
    }
}
enum_error!(DoubleBondStereoError, DoubleBondStereoErrorKind, ck::DoubleBondStereoError, {
    Rings => ck::DoubleBondStereoError::Rings(_), UnsignedPropertyOverflow => ck::DoubleBondStereoError::UnsignedPropertyOverflow { .. },
    InvalidPropertyKind => ck::DoubleBondStereoError::InvalidPropertyKind { .. }, NumericPropertyRead => ck::DoubleBondStereoError::NumericPropertyRead { .. },
    InvalidTopology => ck::DoubleBondStereoError::InvalidTopology(_), BondValue => ck::DoubleBondStereoError::BondValue(_),
    BondOutOfRange => ck::DoubleBondStereoError::BondOutOfRange { .. }, AtomOutOfRange => ck::DoubleBondStereoError::AtomOutOfRange { .. },
    NotDoubleBond => ck::DoubleBondStereoError::NotDoubleBond { .. }, UndefinedStereo => ck::DoubleBondStereoError::UndefinedStereo { .. },
    UnsupportedStereo => ck::DoubleBondStereoError::UnsupportedStereo { .. }, InvalidEndpointDegree => ck::DoubleBondStereoError::InvalidEndpointDegree { .. },
    InvalidDirection => ck::DoubleBondStereoError::InvalidDirection { .. }, RankCount => ck::DoubleBondStereoError::RankCount { .. }, RingInfoNotInitialized => ck::DoubleBondStereoError::RingInfoNotInitialized,
    RingRowCount => ck::DoubleBondStereoError::RingRowCount { .. }, InvalidConformer => ck::DoubleBondStereoError::InvalidConformer(_),
    StereoReferenceRequired => ck::DoubleBondStereoError::StereoReferenceRequired { .. }, StereoReferenceNotNeighbor => ck::DoubleBondStereoError::StereoReferenceNotNeighbor { .. },
    StereoReferenceIsOppositeEndpoint => ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { .. },
    ControllingAtomMismatch => ck::DoubleBondStereoError::ControllingAtomMismatch { .. }, ControllingSlotOverflow => ck::DoubleBondStereoError::ControllingSlotOverflow { .. },
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner {
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { atom, .. } | ck::DoubleBondStereoError::InvalidPropertyKind { atom, .. } | ck::DoubleBondStereoError::NumericPropertyRead { atom, .. } => atom.map(ck::AtomId::index),
        ck::DoubleBondStereoError::AtomOutOfRange { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner {
        ck::DoubleBondStereoError::UnsignedPropertyOverflow { bond, .. } | ck::DoubleBondStereoError::InvalidPropertyKind { bond, .. } | ck::DoubleBondStereoError::NumericPropertyRead { bond, .. } => bond.map(ck::BondId::index),
        ck::DoubleBondStereoError::BondOutOfRange { bond, .. } | ck::DoubleBondStereoError::NotDoubleBond { bond, .. } | ck::DoubleBondStereoError::UndefinedStereo { bond, .. }
        | ck::DoubleBondStereoError::UnsupportedStereo { bond, .. } | ck::DoubleBondStereoError::InvalidEndpointDegree { bond, .. } | ck::DoubleBondStereoError::StereoReferenceRequired { bond, .. }
        | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { bond, .. } | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { bond, .. }
        | ck::DoubleBondStereoError::ControllingAtomMismatch { bond, .. } | ck::DoubleBondStereoError::ControllingSlotOverflow { bond, .. } => Some(bond.index()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn property(&self) -> Option<String> { match *self.inner { ck::DoubleBondStereoError::UnsignedPropertyOverflow { property, .. } | ck::DoubleBondStereoError::InvalidPropertyKind { property, .. } | ck::DoubleBondStereoError::NumericPropertyRead { property, .. } => Some(property.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<u32> { match *self.inner { ck::DoubleBondStereoError::UnsignedPropertyOverflow { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn kind(&self) -> Option<PropertyValueKind> { match *self.inner { ck::DoubleBondStereoError::InvalidPropertyKind { kind, .. } => Some(kind.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::AtomOutOfRange { atom_count, .. } | ck::DoubleBondStereoError::RankCount { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::BondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::RankCount { actual, .. } | ck::DoubleBondStereoError::RingRowCount { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::RingRowCount { expected, .. } => Some(expected), _ => None } }
    #[wasm_bindgen(getter)] pub fn dimension(&self) -> Option<String> { match *self.inner { ck::DoubleBondStereoError::RingRowCount { dimension, .. } => Some(dimension.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::DoubleBondStereoError::NotDoubleBond { order, .. } => Some(order.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn stereo(&self) -> Option<BondStereo> { match *self.inner { ck::DoubleBondStereoError::UndefinedStereo { stereo, .. } | ck::DoubleBondStereoError::UnsupportedStereo { stereo, .. } | ck::DoubleBondStereoError::StereoReferenceRequired { stereo, .. } => Some(stereo.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn endpoint(&self) -> Option<String> { match *self.inner {
        ck::DoubleBondStereoError::InvalidEndpointDegree { endpoint, .. } | ck::DoubleBondStereoError::StereoReferenceNotNeighbor { endpoint, .. }
        | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { endpoint, .. } | ck::DoubleBondStereoError::ControllingAtomMismatch { endpoint, .. }
        | ck::DoubleBondStereoError::ControllingSlotOverflow { endpoint, .. } => Some(endpoint.into()), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn degree(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::InvalidEndpointDegree { degree, .. } | ck::DoubleBondStereoError::ControllingSlotOverflow { degree, .. } => Some(degree), _ => None } }
    #[wasm_bindgen(getter)] pub fn reference(&self) -> Option<usize> { match *self.inner { ck::DoubleBondStereoError::StereoReferenceNotNeighbor { reference, .. } | ck::DoubleBondStereoError::StereoReferenceIsOppositeEndpoint { reference, .. } | ck::DoubleBondStereoError::ControllingAtomMismatch { reference, .. } => Some(reference.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn direction(&self) -> Option<BondDirection> { match *self.inner { ck::DoubleBondStereoError::InvalidDirection { direction } => Some(direction.into()), _ => None } }
}, source = DoubleStereoCause, cause = cause);
