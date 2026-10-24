//! Legacy CIP ranking diagnostics retain exact signed integer widths.
use crate::{
    enums::BondOrder,
    errors::{enum_error, exception},
    model_errors::TopologyValidationError,
    query_errors,
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "MappingValidationError")]
    #[derive(Clone)]
    pub type CipRankCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type CipRankTopology;
    #[wasm_bindgen(typescript_type = "QueryStateError")]
    pub type CipRankQuery;
}
pub(crate) fn cause(value: &ck::CipRankError) -> Option<CipRankCause> {
    match value {
        ck::CipRankError::InvalidQueryState(value) => {
            query_errors::cause(value).map(JsCast::unchecked_into)
        }
        ck::CipRankError::InvalidTopology(_)
        | ck::CipRankError::ValenceRowCount { .. }
        | ck::CipRankError::NegativeImplicitHydrogen { .. }
        | ck::CipRankError::AtomMapOutOfRange { .. }
        | ck::CipRankError::InvariantCount { .. }
        | ck::CipRankError::InvariantOutOfRange { .. }
        | ck::CipRankError::TooManyNeighbors { .. }
        | ck::CipRankError::UnsupportedBondOrder { .. } => None,
    }
}
enum_error!(CipRankError, CipRankErrorKind, ck::CipRankError, {
    InvalidTopology => ck::CipRankError::InvalidTopology(_), InvalidQueryState => ck::CipRankError::InvalidQueryState(_),
    ValenceRowCount => ck::CipRankError::ValenceRowCount { .. }, NegativeImplicitHydrogen => ck::CipRankError::NegativeImplicitHydrogen { .. },
    AtomMapOutOfRange => ck::CipRankError::AtomMapOutOfRange { .. }, InvariantCount => ck::CipRankError::InvariantCount { .. },
    InvariantOutOfRange => ck::CipRankError::InvariantOutOfRange { .. }, TooManyNeighbors => ck::CipRankError::TooManyNeighbors { .. },
    UnsupportedBondOrder => ck::CipRankError::UnsupportedBondOrder { .. },
}, {
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::CipRankError::NegativeImplicitHydrogen { atom, .. } | ck::CipRankError::AtomMapOutOfRange { atom, .. } | ck::CipRankError::InvariantOutOfRange { atom, .. } | ck::CipRankError::TooManyNeighbors { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn value(&self) -> Option<i64> { match *self.inner { ck::CipRankError::NegativeImplicitHydrogen { value, .. } => Some(i64::from(value)), ck::CipRankError::InvariantOutOfRange { value, .. } => Some(value), _ => None } }
    #[wasm_bindgen(getter)] pub fn field(&self) -> Option<String> { match *self.inner { ck::CipRankError::ValenceRowCount { field, .. } => Some(field.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn actual(&self) -> Option<usize> { match *self.inner { ck::CipRankError::ValenceRowCount { actual, .. } | ck::CipRankError::InvariantCount { actual, .. } => Some(actual), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::CipRankError::ValenceRowCount { atom_count, .. } | ck::CipRankError::InvariantCount { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = mapNumber)] pub fn map_number(&self) -> Option<u32> { match *self.inner { ck::CipRankError::AtomMapOutOfRange { map_number, .. } => Some(map_number), _ => None } }
    #[wasm_bindgen(getter)] pub fn degree(&self) -> Option<usize> { match *self.inner { ck::CipRankError::TooManyNeighbors { degree, .. } => Some(degree), _ => None } }
    #[wasm_bindgen(getter, js_name = maximumSupported)] pub fn maximum_supported(&self) -> Option<usize> { match *self.inner { ck::CipRankError::TooManyNeighbors { maximum_supported, .. } => Some(maximum_supported), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::CipRankError::UnsupportedBondOrder { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::CipRankError::UnsupportedBondOrder { order, .. } => Some(order.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<CipRankTopology> { match self.inner.as_ref() { ck::CipRankError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter, js_name = queryState)] pub fn query_state(&self) -> Option<CipRankQuery> { match self.inner.as_ref() { ck::CipRankError::InvalidQueryState(value) => Some(query_errors::QueryStateError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = CipRankCause, cause = cause);
