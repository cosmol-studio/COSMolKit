//! Shared parser-carrier failures, with the native cause edges.
use crate::{
    errors::{enum_error, exception},
    model_errors::AtomPropertyError,
    property_errors::PropertyStringError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "PropertyStringError | AtomPropertyError")]
    #[derive(Clone)]
    pub type ParserCarrierCause;
}
enum_error!(ParserCarrierError, ParserCarrierErrorKind, ck::ParserCarrierError, {
    PropertyString => ck::ParserCarrierError::PropertyString(_),
    AtomProperty => ck::ParserCarrierError::AtomProperty(_),
    Model => ck::ParserCarrierError::Model(_),
}, {
    #[wasm_bindgen(getter)] pub fn model(&self) -> Option<String> { match self.inner.as_ref() { ck::ParserCarrierError::Model(value) => Some(value.clone()), _ => None } }
}, source = ParserCarrierCause, cause = |value: &ck::ParserCarrierError| match value {
    ck::ParserCarrierError::PropertyString(value) => Some(PropertyStringError::from_rust(*value).unchecked_into()),
    ck::ParserCarrierError::AtomProperty(value) => Some(AtomPropertyError::from_rust(value.clone()).unchecked_into()),
    ck::ParserCarrierError::Model(_) => None,
});

enum_error!(ParserStereoOrderError, ParserStereoOrderErrorKind, ck::ParserStereoOrderError, {
    BondIndexOutOfRange => ck::ParserStereoOrderError::BondIndexOutOfRange { .. },
    ProbeLengthMismatch => ck::ParserStereoOrderError::ProbeLengthMismatch { .. },
    ProbeMemberMissing => ck::ParserStereoOrderError::ProbeMemberMissing { .. },
}, {
    #[wasm_bindgen(getter, js_name = bondIndex)] pub fn bond_index(&self) -> Option<usize> { match *self.inner { ck::ParserStereoOrderError::BondIndexOutOfRange { bond_index, .. } => Some(bond_index), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::ParserStereoOrderError::BondIndexOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = incidentCount)] pub fn incident_count(&self) -> Option<usize> { match *self.inner { ck::ParserStereoOrderError::ProbeLengthMismatch { incident_count, .. } => Some(incident_count), _ => None } }
    #[wasm_bindgen(getter, js_name = probeCount)] pub fn probe_count(&self) -> Option<usize> { match *self.inner { ck::ParserStereoOrderError::ProbeLengthMismatch { probe_count, .. } => Some(probe_count), _ => None } }
    #[wasm_bindgen(getter, js_name = probeIndex)] pub fn probe_index(&self) -> Option<usize> { match *self.inner { ck::ParserStereoOrderError::ProbeMemberMissing { probe_index } => Some(probe_index), _ => None } }
});
