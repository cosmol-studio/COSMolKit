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
    #[wasm_bindgen(js_name=fingerprintLayeredList,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_layered_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_layered_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintLayeredListWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_layered_list_with_params(
        &self,
        options: &LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_layered_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputList,unchecked_return_type="(LayeredFingerprintResult | null)[]")]
    pub fn fingerprint_layered_with_output_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_layered_with_output_list()
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
    }
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputListWithParams,unchecked_return_type="(LayeredFingerprintResult | null)[]")]
    pub fn fingerprint_layered_with_output_list_with_params(
        &self,
        options: &LayeredFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_layered_with_output_list_with_params(&options.inner, p)
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
    #[wasm_bindgen(js_name=fingerprintPatternList,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_pattern_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_pattern_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintPatternListWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_pattern_list_with_params(
        &self,
        options: &PatternFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_pattern_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
}
