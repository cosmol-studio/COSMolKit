//! Valence preparation is a typed projection of the canonical operations.
use crate::Molecule;
use crate::alignment_values::operation_error;
use crate::host_values::{bool_value, integer};
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum ValenceModel {
    RdkitLike = 0,
}

#[wasm_bindgen]
pub struct ValenceParams {
    inner: ck::ValenceParams,
}

#[wasm_bindgen]
impl ValenceParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "ValenceModel")] model: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] strict: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::ValenceParams::default();
        if !model.is_undefined() {
            integer(&model, "model", 0.0, 0.0)?;
            inner.model = ck::ValenceModel::RdkitLike;
        }
        if !strict.is_undefined() {
            inner.strict = bool_value(&strict, "strict")?;
        }
        Ok(Self { inner })
    }

    #[wasm_bindgen(getter)]
    pub fn model(&self) -> ValenceModel {
        match self.inner.model {
            ck::ValenceModel::RdkitLike => ValenceModel::RdkitLike,
        }
    }

    #[wasm_bindgen(getter)]
    pub fn strict(&self) -> bool {
        self.inner.strict
    }
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = withAssignedValence)]
    pub fn with_assigned_valence(&self) -> Result<Self, JsValue> {
        self.inner
            .with_assigned_valence()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }

    #[wasm_bindgen(js_name = withAssignedValenceWithParams)]
    pub fn with_assigned_valence_with_params(
        &self,
        params: &ValenceParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_assigned_valence_with_params(&params.inner)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }

    #[wasm_bindgen(js_name = assignValence)]
    pub fn assign_valence_(&self) -> Result<(), JsValue> {
        self.inner
            .assign_valence_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }

    #[wasm_bindgen(js_name = assignValenceWithParams)]
    pub fn assign_valence_with_params_(&self, params: &ValenceParams) -> Result<(), JsValue> {
        self.inner
            .assign_valence_with_params_(&params.inner)
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
