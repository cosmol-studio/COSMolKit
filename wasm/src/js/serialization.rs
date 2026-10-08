//! Uint8Array transport without an alternative archive codec.
use crate::Molecule;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toBinary)]
    pub fn to_binary(&self) -> Result<Vec<u8>, JsValue> {
        // COSMolKit❗✔️: self.inner.to_binary()
        self.inner
            .to_binary()
            .map_err(|e| crate::pickle_errors::pickle_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromBinary)]
    pub fn from_binary(
        #[wasm_bindgen(unchecked_param_type = "Uint8Array")] data: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::Molecule::from_binary(data)
        let data = data
            .dyn_ref::<js_sys::Uint8Array>()
            .ok_or_else(|| JsValue::from(js_sys::TypeError::new("data must be a Uint8Array")))?;
        cosmolkit_wasm::Molecule::from_binary(&data.to_vec())
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| crate::pickle_errors::pickle_error(&e).unwrap_or_else(|e| e))
    }
}
