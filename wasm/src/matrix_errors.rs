//! Matrix errors with retained transparent topology/coordinate payloads.
use crate::{
    enums::BondOrder,
    errors::{enum_error, exception},
    model_errors::{self, CoordinateValidationError, TopologyValidationError},
};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "StereoGroupError | TemplateAttachmentOrderError")]
    #[derive(Clone)]
    pub type MatrixCause;
    #[wasm_bindgen(typescript_type = "TopologyValidationError")]
    pub type MatrixTopology;
    #[wasm_bindgen(typescript_type = "CoordinateValidationError")]
    pub type MatrixCoordinates;
}
pub(crate) fn cause(value: &ck::MatrixError) -> Option<MatrixCause> {
    match value {
        ck::MatrixError::InvalidTopology(value) => {
            model_errors::topology_cause(value).map(JsCast::unchecked_into)
        }
        ck::MatrixError::InvalidCoordinates(_)
        | ck::MatrixError::ActiveBondsWithoutAtoms
        | ck::MatrixError::ActiveAtomOutOfRange { .. }
        | ck::MatrixError::DuplicateActiveAtom { .. }
        | ck::MatrixError::ActiveBondOutOfRange { .. }
        | ck::MatrixError::DuplicateActiveBond { .. }
        | ck::MatrixError::ActiveBondEndpointMissing { .. }
        | ck::MatrixError::UnsupportedBondOrder { .. }
        | ck::MatrixError::No3dConformer
        | ck::MatrixError::ConformerNotFound { .. }
        | ck::MatrixError::MatrixDimensionOverflow { .. } => None,
    }
}
enum_error!(MatrixError, MatrixErrorKind, ck::MatrixError, {
    InvalidTopology => ck::MatrixError::InvalidTopology(_), InvalidCoordinates => ck::MatrixError::InvalidCoordinates(_), ActiveBondsWithoutAtoms => ck::MatrixError::ActiveBondsWithoutAtoms,
    ActiveAtomOutOfRange => ck::MatrixError::ActiveAtomOutOfRange { .. }, DuplicateActiveAtom => ck::MatrixError::DuplicateActiveAtom { .. },
    ActiveBondOutOfRange => ck::MatrixError::ActiveBondOutOfRange { .. }, DuplicateActiveBond => ck::MatrixError::DuplicateActiveBond { .. },
    ActiveBondEndpointMissing => ck::MatrixError::ActiveBondEndpointMissing { .. }, UnsupportedBondOrder => ck::MatrixError::UnsupportedBondOrder { .. },
    No3dConformer => ck::MatrixError::No3dConformer, ConformerNotFound => ck::MatrixError::ConformerNotFound { .. }, MatrixDimensionOverflow => ck::MatrixError::MatrixDimensionOverflow { .. },
}, {
    #[wasm_bindgen(getter)] pub fn position(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ActiveAtomOutOfRange { position, .. } | ck::MatrixError::ActiveBondOutOfRange { position, .. } => Some(position), _ => None } }
    #[wasm_bindgen(getter)] pub fn atom(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ActiveAtomOutOfRange { atom, .. } | ck::MatrixError::DuplicateActiveAtom { atom, .. } | ck::MatrixError::ActiveBondEndpointMissing { atom, .. } => Some(atom.index()), _ => None } }
    #[wasm_bindgen(getter)] pub fn bond(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ActiveBondOutOfRange { bond, .. } | ck::MatrixError::DuplicateActiveBond { bond, .. } | ck::MatrixError::ActiveBondEndpointMissing { bond, .. } | ck::MatrixError::UnsupportedBondOrder { bond, .. } => Some(bond.index()), _ => None } }
    #[wasm_bindgen(getter, js_name = atomCount)] pub fn atom_count(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ActiveAtomOutOfRange { atom_count, .. } => Some(atom_count), _ => None } }
    #[wasm_bindgen(getter, js_name = bondCount)] pub fn bond_count(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ActiveBondOutOfRange { bond_count, .. } => Some(bond_count), _ => None } }
    #[wasm_bindgen(getter, js_name = firstPosition)] pub fn first_position(&self) -> Option<usize> { match *self.inner { ck::MatrixError::DuplicateActiveAtom { first_position, .. } | ck::MatrixError::DuplicateActiveBond { first_position, .. } => Some(first_position), _ => None } }
    #[wasm_bindgen(getter, js_name = secondPosition)] pub fn second_position(&self) -> Option<usize> { match *self.inner { ck::MatrixError::DuplicateActiveAtom { second_position, .. } | ck::MatrixError::DuplicateActiveBond { second_position, .. } => Some(second_position), _ => None } }
    #[wasm_bindgen(getter)] pub fn endpoint(&self) -> Option<String> { match *self.inner { ck::MatrixError::ActiveBondEndpointMissing { endpoint, .. } => Some(endpoint.into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn order(&self) -> Option<BondOrder> { match *self.inner { ck::MatrixError::UnsupportedBondOrder { order, .. } => Some(order.into()), _ => None } }
    #[wasm_bindgen(getter, js_name = conformerId)] pub fn conformer_id(&self) -> Option<usize> { match *self.inner { ck::MatrixError::ConformerNotFound { conformer_id } => Some(conformer_id), _ => None } }
    #[wasm_bindgen(getter)] pub fn dimension(&self) -> Option<usize> { match *self.inner { ck::MatrixError::MatrixDimensionOverflow { dimension } => Some(dimension), _ => None } }
    #[wasm_bindgen(getter)] pub fn topology(&self) -> Option<MatrixTopology> { match self.inner.as_ref() { ck::MatrixError::InvalidTopology(value) => Some(TopologyValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
    #[wasm_bindgen(getter)] pub fn coordinates(&self) -> Option<MatrixCoordinates> { match self.inner.as_ref() { ck::MatrixError::InvalidCoordinates(value) => Some(CoordinateValidationError::from_rust(value.clone()).unchecked_into()), _ => None } }
}, source = MatrixCause, cause = cause);
