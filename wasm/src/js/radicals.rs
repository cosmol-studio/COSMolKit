//! Radical assignment preserves canonical value and mutation semantics.
use crate::Molecule;
use crate::alignment_values::operation_error;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAssignedRadicals)]
    pub fn with_assigned_radicals(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_assigned_radicals()
        self.inner
            .with_assigned_radicals()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=assignRadicals)]
    pub fn assign_radicals_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assign_radicals_()
        self.inner
            .assign_radicals_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
