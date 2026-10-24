//! Registered molecule fragment operations retain native ordering and ownership.
use crate::{Molecule, fields::Field, operation_errors::OperationError, storage::NativeField};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
crate::fields::value_array!(ck::Molecule, Molecules, "ReadonlyArray<Molecule>");
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fragments)]
    pub fn fragments(&self) -> Result<Molecules, JsValue> {
        Vec::<ck::Molecule>::get(NativeField::owned(
            self.0.fragments().map_err(OperationError::from_rust)?,
        ))
    }
    #[wasm_bindgen(js_name=largestFragment)]
    pub fn largest_fragment(&self) -> Result<Self, JsValue> {
        self.0
            .largest_fragment()
            .map(Self)
            .map_err(OperationError::from_rust)
    }
}
