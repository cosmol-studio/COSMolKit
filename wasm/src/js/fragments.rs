//! JavaScript ownership/error conversion only; chemistry stays in the facade.
use crate::Molecule;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    /// Sanitized components in source order, including conformers; self is unchanged.
    pub fn fragments(&self) -> Result<Vec<Molecule>, JsValue> {
        self.inner.fragments().map(|values| values.into_iter().map(|inner| Self {
            inner: std::sync::Arc::new(inner),
        }).collect()).map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
    /// Most atoms, last on ties. Empty input throws OperationError.
    #[wasm_bindgen(js_name = largestFragment)]
    pub fn largest_fragment(&self) -> Result<Self, JsValue> {
        self.inner.largest_fragment().map(|inner| Self { inner: std::sync::Arc::new(inner) })
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e| e))
    }
}
