//! Full concrete and query Layered/Pattern public API projections.
use crate::Molecule;
use crate::fingerprint_values::Fingerprint;
use crate::layered_pattern_errors::{layered_error, pattern_error};
use crate::layered_pattern_values::{
    LayeredFingerprintParams, LayeredFingerprintResult, PatternFingerprintParams,
};
use crate::query_construction::QueryGraph;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=layeredFingerprint)]
    pub fn layered_fingerprint(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint(
        self.inner
            .layered_fingerprint()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=layeredFingerprintWithParams)]
    pub fn layered_fingerprint_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_params(
        self.inner
            .layered_fingerprint_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=layeredFingerprintWithOutput)]
    pub fn layered_fingerprint_with_output(&self) -> Result<LayeredFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output(
        self.inner
            .layered_fingerprint_with_output()
            .map(|inner| LayeredFingerprintResult { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=layeredFingerprintWithOutputWithParams)]
    pub fn layered_fingerprint_with_output_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<LayeredFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output_with_params(
        self.inner
            .layered_fingerprint_with_output_with_params(&params.inner)
            .map(|inner| LayeredFingerprintResult { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen(js_name=layeredQueryFingerprintWithParams)]
pub fn layered_query_fingerprint_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::layered_query_fingerprint_with_params(&query.inner, &params.inner)
    ck::layered_query_fingerprint_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen(js_name=layeredQueryFingerprintWithOutputWithParams)]
pub fn layered_query_fingerprint_with_output_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, JsValue> {
    // COSMolKit❗✔️: ck::layered_query_fingerprint_with_output_with_params(&query.inner, &params.inner)
    ck::layered_query_fingerprint_with_output_with_params(&query.inner, &params.inner)
        .map(|inner| LayeredFingerprintResult { inner })
        .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=patternFingerprint)]
    pub fn pattern_fingerprint(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .pattern_fingerprint(
        self.inner
            .pattern_fingerprint()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=patternFingerprintWithParams)]
    pub fn pattern_fingerprint_with_params(
        &self,
        params: &PatternFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .pattern_fingerprint_with_params(
        self.inner
            .pattern_fingerprint_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen(js_name=patternQueryFingerprint)]
pub fn pattern_query_fingerprint(query: &QueryGraph) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::pattern_query_fingerprint(&query.inner)
    ck::pattern_query_fingerprint(&query.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen(js_name=patternQueryFingerprintWithParams)]
pub fn pattern_query_fingerprint_with_params(
    query: &QueryGraph,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::pattern_query_fingerprint_with_params(&query.inner, &params.inner)
    ck::pattern_query_fingerprint_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
}
