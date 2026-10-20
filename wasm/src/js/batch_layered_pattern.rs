//! All six Layered/Pattern batch methods delegate to the canonical owner.
use crate::batch_boundary::{MoleculeBatch, batch_validation_error};
use crate::batch_queries::BatchQueryParams;
use crate::fingerprint_values::Fingerprint;
use crate::layered_pattern_values::{
    LayeredFingerprintParams, LayeredFingerprintResult, PatternFingerprintParams,
};
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=fingerprintLayered,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_layered(&self, #[wasm_bindgen(unchecked_optional_param_type = "LayeredFingerprintParams | LayeredFingerprintOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_layered()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })

            };
        }
        let configuration_0 = crate::layered_pattern_values::LayeredFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_layered_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_layered_with_params(
        &self,
        options: &LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_layered_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutput,unchecked_return_type="(LayeredFingerprintResult | null)[]")]
    pub fn fingerprint_layered_with_output(&self, #[wasm_bindgen(unchecked_optional_param_type = "LayeredFingerprintParams | LayeredFingerprintOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_layered_with_output()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            LayeredFingerprintResult { inner }.into()
                        })
                    })
                    .collect()
            })

            };
        }
        let configuration_0 = crate::layered_pattern_values::LayeredFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_layered_with_output_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputWithParams,unchecked_return_type="(LayeredFingerprintResult | null)[]")]
    pub fn fingerprint_layered_with_output_with_params(
        &self,
        options: &LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_layered_with_output_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            LayeredFingerprintResult { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintPattern,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_pattern(&self, #[wasm_bindgen(unchecked_optional_param_type = "PatternFingerprintParams | PatternFingerprintParamsOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_pattern()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })

            };
        }
        let configuration_0 = crate::layered_pattern_values::PatternFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_pattern_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintPatternWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_pattern_with_params(
        &self,
        options: &PatternFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_pattern_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
}
