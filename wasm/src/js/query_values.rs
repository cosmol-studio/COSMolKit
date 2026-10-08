//! Foundational detached query value; parsing and matching require search.
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct QueryGraph {
    pub(crate) inner: ck::QueryGraph,
}

#[wasm_bindgen]
impl QueryGraph {
    #[wasm_bindgen(js_name=numAtoms)]
    pub fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    #[wasm_bindgen(js_name=numBonds)]
    pub fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn name(&self) -> Result<JsValue, JsValue> {
        self.inner
            .prop("_Name")
            .map(|value| {
                value
                    .as_string()
                    .map_err(|e| crate::property_values::property_error(&e).unwrap_or_else(|e| e))
                    .and_then(crate::host_values::text)
            })
            .transpose()
            .map(|value| value.map_or(JsValue::NULL, JsValue::from))
    }
}
