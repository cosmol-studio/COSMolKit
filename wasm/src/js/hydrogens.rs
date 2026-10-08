//! Full value/in-place hydrogen operations sharing complete batch parameter values.
use crate::Molecule;
use crate::alignment_values::operation_error;
use crate::transform_parameters::{AddHsParams, RemoveHsParams};
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(
    inline_js = "export function visitHydrogenAddParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid AddHsParams',{cause});}} export function visitHydrogenRemoveParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid RemoveHsParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitHydrogenAddParams)]
    fn visit_add_params(
        value: &JsValue,
        visit: &mut dyn FnMut(&AddHsParams),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitHydrogenRemoveParams)]
    fn visit_remove_params(
        value: &JsValue,
        visit: &mut dyn FnMut(&RemoveHsParams),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withHydrogens)]
    pub fn with_hydrogens(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_hydrogens()
        self.inner
            .with_hydrogens()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withHydrogensWithParams)]
    pub fn with_hydrogens_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "AddHsParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_hydrogens_with_params(&params.inner)
        let mut result = None;
        visit_add_params(&params, &mut |params: &AddHsParams| {
            result = Some(
                self.inner
                    .with_hydrogens_with_params(&params.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
    #[wasm_bindgen(js_name=withoutHydrogens)]
    pub fn without_hydrogens(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.without_hydrogens()
        self.inner
            .without_hydrogens()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withoutHydrogensWithParams)]
    pub fn without_hydrogens_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "RemoveHsParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.without_hydrogens_with_params(&params.inner)
        let mut result = None;
        visit_remove_params(&params, &mut |params: &RemoveHsParams| {
            result = Some(
                self.inner
                    .without_hydrogens_with_params(&params.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
    #[wasm_bindgen(js_name=addHydrogens)]
    pub fn add_hydrogens_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.add_hydrogens_()
        self.inner
            .add_hydrogens_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=addHydrogensWithParams)]
    pub fn add_hydrogens_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "AddHsParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.add_hydrogens_with_params_(&params.inner)
        let mut result = None;
        visit_add_params(&params, &mut |params: &AddHsParams| {
            result = Some(
                self.inner
                    .add_hydrogens_with_params_(&params.inner)
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
    #[wasm_bindgen(js_name=removeHydrogens)]
    pub fn remove_hydrogens_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.remove_hydrogens_()
        self.inner
            .remove_hydrogens_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=removeHydrogensWithParams)]
    pub fn remove_hydrogens_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "RemoveHsParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.remove_hydrogens_with_params_(&params.inner)
        let mut result = None;
        visit_remove_params(&params, &mut |params: &RemoveHsParams| {
            result = Some(
                self.inner
                    .remove_hydrogens_with_params_(&params.inner)
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| crate::host_values::type_error("params"))?
    }
}
