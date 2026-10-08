//! Canonical embedding parameter errors and original nested execution failures.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct ConformerError {
    kind: String,
    message: String,
    detail: Option<String>,
}
#[wasm_bindgen]
impl ConformerError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "conformer".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn detail(&self) -> JsValue {
        self.detail
            .as_ref()
            .map_or(JsValue::NULL, |v| v.as_str().into())
    }
}
pub(crate) fn params_error(source: &ck::ConformerError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: pub enum ConformerError {
    let (kind, detail) = match source {
        ck::ConformerError::Unsupported => ("Unsupported", None),
        ck::ConformerError::CannotNormalizeZeroLengthVector => {
            ("CannotNormalizeZeroLengthVector", None)
        }
        ck::ConformerError::GenerationFailed(v) => ("GenerationFailed", Some(v.clone())),
        ck::ConformerError::InvalidEmbedParametersJson(v) => {
            ("InvalidEmbedParametersJson", Some(v.clone()))
        }
        ck::ConformerError::WasmImplicitClockSeedUnsupported => {
            ("WasmImplicitClockSeedUnsupported", None)
        }
    };
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("ConformerError");
    let e: JsValue = e.into();
    set(&e, "domain", "conformer".into())?;
    set(&e, "kind", kind.into())?;
    set(
        &e,
        "detail",
        ConformerError {
            kind: kind.into(),
            message: source.to_string(),
            detail,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct ConformerRunError {
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl ConformerRunError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "conformer".into()
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
pub(crate) fn run_error(source: &ck::ConformerRunError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: fmt::Display::fmt(self.source().unwrap(), f)
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("ConformerRunError");
    let e: JsValue = e.into();
    set(&e, "domain", "conformer".into())?;
    set(
        &e,
        "detail",
        ConformerRunError {
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
