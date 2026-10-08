//! Reuse full shared KekulizeParams and error vocabulary for molecule operations.
use crate::Molecule;
use crate::alignment_values::operation_error;
use crate::host_values::type_error;
use crate::transform_parameters::KekulizeParams;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(
    inline_js = "export function visitKekulizeParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid KekulizeParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitKekulizeParams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&KekulizeParams)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withKekulizedBonds)]
    pub fn with_kekulized_bonds(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_kekulized_bonds()
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withKekulizedBondsWithParams)]
    pub fn with_kekulized_bonds_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "KekulizeParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_kekulized_bonds_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &KekulizeParams| {
            result = Some(
                self.inner
                    .with_kekulized_bonds_with_params(&p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=kekulizeBonds)]
    pub fn kekulize_bonds_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.kekulize_bonds_()
        self.inner
            .kekulize_bonds_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=kekulizeBondsWithParams)]
    pub fn kekulize_bonds_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "KekulizeParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.kekulize_bonds_with_params_(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &KekulizeParams| {
            result = Some(
                self.inner
                    .kekulize_bonds_with_params_(&p.inner)
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
