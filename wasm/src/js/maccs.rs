//! MACCS values, all registered calls and complete source error context.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::fingerprint_values::Fingerprint;
use crate::host_values::usize_value;
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MaccsFingerprintParams {
    inner: ck::MaccsFingerprintParams,
}
#[wasm_bindgen]
impl MaccsFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] n_bits: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::MaccsFingerprintParams { n_bits },
        let mut inner = ck::MaccsFingerprintParams::default();
        if !n_bits.is_undefined() {
            inner.n_bits = usize_value(&n_bits, "nBits")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> usize {
        // COSMolKit❗✔️: self.inner.n_bits
        self.inner.n_bits
    }
}
#[wasm_bindgen]
pub struct MaccsFingerprintError {
    kind: String,
    message: String,
    cause: JsValue,
    option: Option<String>,
    reason: Option<String>,
    bit: Option<usize>,
}
#[wasm_bindgen]
impl MaccsFingerprintError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "Fingerprint".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn option(&self) -> JsValue {
        self.option
            .as_deref()
            .map_or(JsValue::NULL, JsValue::from_str)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn reason(&self) -> JsValue {
        self.reason
            .as_deref()
            .map_or(JsValue::NULL, JsValue::from_str)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn bit(&self) -> JsValue {
        self.bit
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
    }
}
pub(crate) fn maccs_error(source: &ck::MaccsFingerprintError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: object.setattr("domain", "Fingerprint")?;
    use ck::MaccsFingerprintError as E;
    let kind = match source {
        E::UnsupportedOption { .. } => "UnsupportedOption",
        E::MissingPattern { .. } => "MissingPattern",
        E::Topology(_) => "Topology",
        E::Rings(_) => "Rings",
        E::Valence(_) => "Valence",
        E::Paths(_) => "Paths",
        E::Smarts(_) => "Smarts",
        E::QueryCompile(_) => "QueryCompile",
        E::QueryContext(_) => "QueryContext",
        E::Match(_) => "Match",
        E::Value(_) => "Value",
    };
    let (option, reason, bit) = match source {
        E::UnsupportedOption { option, reason } => {
            (Some((*option).to_owned()), Some((*reason).to_owned()), None)
        }
        E::MissingPattern { bit } => (None, None, Some(*bit)),
        _ => (None, None, None),
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("MaccsFingerprintError");
    let e: JsValue = e.into();
    set(&e, "domain", "Fingerprint".into())?;
    set(&e, "kind", kind.into())?;
    if let Some(v) = &option {
        set(&e, "option", v.as_str().into())?;
    }
    if let Some(v) = &reason {
        set(&e, "reason", v.as_str().into())?;
    }
    if let Some(v) = bit {
        set(&e, "bit", JsValue::from_f64(v as f64))?;
    }
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        MaccsFingerprintError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            option,
            reason,
            bit,
        }
        .into(),
    )?;
    Ok(e)
}
fn error(e: ck::MaccsFingerprintError) -> JsValue {
    maccs_error(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=maccsFingerprint)]
    pub fn maccs_fingerprint(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .maccs_fingerprint()
        self.inner
            .maccs_fingerprint()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=maccsFingerprintRaw)]
    pub fn maccs_fingerprint_raw(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .maccs_fingerprint_raw()
        self.inner
            .maccs_fingerprint_raw()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=maccsFingerprintWithParams)]
    pub fn maccs_fingerprint_with_params(
        &self,
        params: &MaccsFingerprintParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .maccs_fingerprint_with_params(&params.inner)
        self.inner
            .maccs_fingerprint_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
}
