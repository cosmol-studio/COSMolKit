//! Precise native causes for internal Protein viewer transport.
use crate::{bio_operations::ProteinProjectionError, errors::enum_error};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};

#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ProteinProjectionError")]
    #[derive(Clone)]
    pub type ProjectionCause;
}
enum_error!(ProteinTransferError, ProteinTransferErrorKind, ck::ProteinTransferError, {
    VersionMismatch => ck::ProteinTransferError::VersionMismatch { .. },
    InvalidData => ck::ProteinTransferError::InvalidData { .. },
    Projection => ck::ProteinTransferError::Projection(_),
}, {
    #[wasm_bindgen(getter)]
    pub fn expected(&self) -> Option<String> { match self.inner.as_ref() {
        ck::ProteinTransferError::VersionMismatch { expected, .. } => Some(expected.clone()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn actual(&self) -> Option<String> { match self.inner.as_ref() {
        ck::ProteinTransferError::VersionMismatch { actual, .. } => Some(actual.clone()), _ => None,
    } }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> Option<String> { match self.inner.as_ref() {
        ck::ProteinTransferError::InvalidData { message } => Some(message.clone()), _ => None,
    } }
}, source = ProjectionCause, cause = |owner: &crate::storage::ReadonlyValue<ck::ProteinTransferError>| match owner.as_ref() {
    ck::ProteinTransferError::Projection(_) => Some(ProteinProjectionError::from_owner(owner.child(|error| match error {
        ck::ProteinTransferError::Projection(error) => error,
        _ => unreachable!("immutable protein transfer cause"),
    })).unchecked_into()), _ => None,
});
