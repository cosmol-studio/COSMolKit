//! Canonical Layered and Pattern errors retain the Python domain spelling.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct LayeredFingerprintError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl LayeredFingerprintError {
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
}
pub(crate) fn layered_error(source: &ck::LayeredFingerprintError) -> Result<JsValue, JsValue> {
    use ck::LayeredFingerprintError as E;
    let kind = match source {
        E::InvalidArguments { .. } => "InvalidArguments",
        E::Topology(_) => "Topology",
        E::Query(_) => "Query",
        E::Rings(_) => "Rings",
        E::Paths(_) => "Paths",
        E::Value(_) => "Value",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("LayeredFingerprintError");
    let error: JsValue = error.into();
    set(&error, "domain", "Fingerprint".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        LayeredFingerprintError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    if let E::InvalidArguments { reason } = source {
        set(&error, "reason", (*reason).into())?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct PatternFingerprintError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl PatternFingerprintError {
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
}
pub(crate) fn pattern_error(source: &ck::PatternFingerprintError) -> Result<JsValue, JsValue> {
    use ck::PatternFingerprintError as E;
    let kind = match source {
        E::EmptyFingerprint => "EmptyFingerprint",
        E::InvalidArguments { .. } => "InvalidArguments",
        E::BitLengthMismatch { .. } => "BitLengthMismatch",
        E::Topology(_) => "Topology",
        E::Query(_) => "Query",
        E::QueryCarrier(_) => "QueryCarrier",
        E::Rings(_) => "Rings",
        E::Smarts(_) => "Smarts",
        E::QueryCompile(_) => "QueryCompile",
        E::QueryContext(_) => "QueryContext",
        E::Match(_) => "Match",
        E::Value(_) => "Value",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("PatternFingerprintError");
    let error: JsValue = error.into();
    set(&error, "domain", "Fingerprint".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        PatternFingerprintError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    if let E::InvalidArguments { reason } = source {
        set(&error, "reason", (*reason).into())?;
    }
    if let E::BitLengthMismatch { left, right } = source {
        set(&error, "left", JsValue::from_f64(*left as f64))?;
        set(&error, "right", JsValue::from_f64(*right as f64))?;
    }
    Ok(error)
}
