//! Detached query-row diagnostics.
use crate::{
    errors::{enum_error, exception},
    model_errors::MappingValidationError,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "number | readonly [number, number]")]
    pub type QueryIndexPayload;
    #[wasm_bindgen(typescript_type = "MappingValidationError")]
    #[derive(Clone)]
    pub type QueryStateCause;
}
fn endpoints(value: (ck::AtomId, ck::AtomId)) -> QueryIndexPayload {
    let pair = js_sys::Array::of2(
        &JsValue::from(value.0.index()),
        &JsValue::from(value.1.index()),
    );
    js_sys::Object::freeze(&pair).unchecked_into()
}
pub(crate) fn cause(value: &ck::QueryStateError) -> Option<QueryStateCause> {
    match value {
        ck::QueryStateError::Mapping(source) => {
            Some(MappingValidationError::from_rust(source.clone()).unchecked_into())
        }
        ck::QueryStateError::NonElementAtomIdentity { .. }
        | ck::QueryStateError::AtomCount { .. }
        | ck::QueryStateError::BondCount { .. }
        | ck::QueryStateError::AtomId { .. }
        | ck::QueryStateError::BondId { .. }
        | ck::QueryStateError::BondEndpoints { .. }
        | ck::QueryStateError::AppendedRow { .. }
        | ck::QueryStateError::AppendedRowCount { .. } => None,
    }
}
enum_error!(QueryStateError, QueryStateErrorKind, ck::QueryStateError, {
    NonElementAtomIdentity => ck::QueryStateError::NonElementAtomIdentity { .. }, AtomCount => ck::QueryStateError::AtomCount { .. },
    BondCount => ck::QueryStateError::BondCount { .. }, AtomId => ck::QueryStateError::AtomId { .. }, BondId => ck::QueryStateError::BondId { .. },
    BondEndpoints => ck::QueryStateError::BondEndpoints { .. }, Mapping => ck::QueryStateError::Mapping(_),
    AppendedRow => ck::QueryStateError::AppendedRow { .. }, AppendedRowCount => ck::QueryStateError::AppendedRowCount { .. },
}, {
    #[wasm_bindgen(getter)] pub fn position(&self) -> Option<usize> { match *self.inner {
        ck::QueryStateError::NonElementAtomIdentity { position, .. } | ck::QueryStateError::AtomId { position, .. } | ck::QueryStateError::BondId { position, .. }
        | ck::QueryStateError::BondEndpoints { position, .. } | ck::QueryStateError::AppendedRow { position, .. } => Some(position), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::QueryStateError::NonElementAtomIdentity { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomicNumber)] pub fn atomic_number(&self) -> Option<u8> { match *self.inner { ck::QueryStateError::NonElementAtomIdentity { atomic_number, .. } => Some(atomic_number), _ => None } }
    #[wasm_bindgen(getter)] pub fn entity(&self) -> Option<String> { match *self.inner { ck::QueryStateError::AppendedRow { entity, .. } | ck::QueryStateError::AppendedRowCount { entity, .. } => Some(entity.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<QueryIndexPayload> { match *self.inner {
        ck::QueryStateError::AtomCount { actual, .. } | ck::QueryStateError::BondCount { actual, .. } | ck::QueryStateError::AppendedRowCount { actual, .. } => Some(JsValue::from(actual).unchecked_into()),
        ck::QueryStateError::AtomId { actual, .. } => Some(JsValue::from(actual.index()).unchecked_into()),
        ck::QueryStateError::BondId { actual, .. } => Some(JsValue::from(actual.index()).unchecked_into()),
        ck::QueryStateError::BondEndpoints { actual, .. } => Some(endpoints(actual)), _ => None,
    } }
    #[wasm_bindgen(getter)] pub fn expected(&self) -> Option<QueryIndexPayload> { match *self.inner {
        ck::QueryStateError::AtomCount { expected, .. } | ck::QueryStateError::BondCount { expected, .. } | ck::QueryStateError::AppendedRowCount { expected, .. } => Some(JsValue::from(expected).unchecked_into()),
        ck::QueryStateError::AtomId { expected, .. } => Some(JsValue::from(expected.index()).unchecked_into()),
        ck::QueryStateError::BondId { expected, .. } => Some(JsValue::from(expected.index()).unchecked_into()),
        ck::QueryStateError::BondEndpoints { expected, .. } => Some(endpoints(expected)), _ => None,
    } }
}, source = QueryStateCause, cause = cause);
