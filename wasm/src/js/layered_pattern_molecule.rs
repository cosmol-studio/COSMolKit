//! Full concrete and query Layered/Pattern public API projections.
use crate::Molecule;
use crate::fingerprint_values::Fingerprint;
use crate::layered_pattern_errors::{layered_error, pattern_error};
use crate::layered_pattern_values::{
    LayeredFingerprintParams, LayeredFingerprintResult, PatternFingerprintParams,
};
use crate::query_values::QueryGraph;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(inline_js = "export function visitLayeredParams(value,visit){visit(value);}")]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitLayeredParams)]
    fn visit_params(value: &JsValue, visit: &mut dyn FnMut(&LayeredFingerprintParams)) -> Result<(), JsValue>;
}
impl LayeredFingerprintParams {
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        if visit_params(value, &mut |params: &LayeredFingerprintParams| {
            inner = Some(params.inner.clone());
        }).is_ok() {
            return inner.map(|inner| Self { inner }).ok_or_else(|| js_sys::TypeError::new("invalid LayeredFingerprintParams").into());
        }
        Self::from_js_options(value)
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayered)]
    pub fn fingerprint_layered(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "LayeredFingerprintParams | LayeredFingerprintOptions")]
        params: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint(
        let configured = (!params.is_undefined()).then(|| LayeredFingerprintParams::from_configuration(&params)).transpose()?;
        configured.as_ref().map_or_else(
            || self.inner.fingerprint_layered(),
            |params| self.inner.fingerprint_layered_with_params(&params.inner),
        )
            .map(|inner| Fingerprint { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayeredWithParams)]
    pub fn fingerprint_layered_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_params(
        self.inner
            .fingerprint_layered_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutput)]
    pub fn fingerprint_layered_with_output(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "LayeredFingerprintParams | LayeredFingerprintOptions")]
        params: JsValue,
    ) -> Result<LayeredFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output(
        let configured = (!params.is_undefined()).then(|| LayeredFingerprintParams::from_configuration(&params)).transpose()?;
        configured.as_ref().map_or_else(
            || self.inner.fingerprint_layered_with_output(),
            |params| self.inner.fingerprint_layered_with_output_with_params(&params.inner),
        )
            .map(|inner| LayeredFingerprintResult { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintLayeredWithOutputWithParams)]
    pub fn fingerprint_layered_with_output_with_params(
        &self,
        params: &LayeredFingerprintParams,
    ) -> Result<LayeredFingerprintResult, JsValue> {
        // COSMolKit❗✔️: .layered_fingerprint_with_output_with_params(
        self.inner
            .fingerprint_layered_with_output_with_params(&params.inner)
            .map(|inner| LayeredFingerprintResult { inner })
            .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen(js_name=fingerprintLayeredQueryWithParams)]
pub fn fingerprint_layered_query_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::layered_query_fingerprint_with_params(&query.inner, &params.inner)
    ck::fingerprint_layered_query_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen(js_name=fingerprintLayeredQueryWithOutputWithParams)]
pub fn fingerprint_layered_query_with_output_with_params(
    query: &QueryGraph,
    params: &LayeredFingerprintParams,
) -> Result<LayeredFingerprintResult, JsValue> {
    // COSMolKit❗✔️: ck::layered_query_fingerprint_with_output_with_params(&query.inner, &params.inner)
    ck::fingerprint_layered_query_with_output_with_params(&query.inner, &params.inner)
        .map(|inner| LayeredFingerprintResult { inner })
        .map_err(|e| layered_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintPattern)]
    pub fn fingerprint_pattern(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .pattern_fingerprint(
        self.inner
            .fingerprint_pattern()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintPatternWithParams)]
    pub fn fingerprint_pattern_with_params(
        &self,
        params: &PatternFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .pattern_fingerprint_with_params(
        self.inner
            .fingerprint_pattern_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen(js_name=fingerprintPatternQuery)]
pub fn fingerprint_pattern_query(query: &QueryGraph) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::pattern_query_fingerprint(&query.inner)
    ck::fingerprint_pattern_query(&query.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
}
#[wasm_bindgen(js_name=fingerprintPatternQueryWithParams)]
pub fn fingerprint_pattern_query_with_params(
    query: &QueryGraph,
    params: &PatternFingerprintParams,
) -> Result<Fingerprint, JsValue> {
    // COSMolKit❗✔️: ck::pattern_query_fingerprint_with_params(&query.inner, &params.inner)
    ck::fingerprint_pattern_query_with_params(&query.inner, &params.inner)
        .map(|inner| Fingerprint { inner })
        .map_err(|e| pattern_error(&e).unwrap_or_else(|e| e))
}
