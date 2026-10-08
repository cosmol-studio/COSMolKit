//! Source filesystem and image errors are observable, never converted to success.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BatchImageError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl BatchImageError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "batch".into()
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
pub(crate) fn batch_image_error(source: &ck::BatchImageError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::BatchImageError::InvalidFormat(..) => "InvalidFormat",
        ck::BatchImageError::Write(..) => "Write",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BatchImageError");
    let error: JsValue = error.into();
    set(&error, "domain", "batch".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        BatchImageError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if let ck::BatchImageError::InvalidFormat(format) = source {
        set(&error, "format", JsValue::from_str(format))?;
    }
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct DrawingWriteError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl DrawingWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "drawing".into()
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
pub(crate) fn drawing_write_error(source: &ck::DrawingWriteError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::DrawingWriteError::Drawing(..) => "Drawing",
        ck::DrawingWriteError::Io { .. } => "Io",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("DrawingWriteError");
    let error: JsValue = error.into();
    set(&error, "domain", "drawing".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        DrawingWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if let ck::DrawingWriteError::Io { path, source } = source {
        set(
            &error,
            "filename",
            JsValue::from_str(&path.to_string_lossy()),
        )?;
        set(
            &error,
            "errno",
            source
                .raw_os_error()
                .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64)),
        )?;
    }
    set(&error, "cause", cause)?;
    Ok(error)
}
pub(crate) fn io_error(source: &std::io::Error) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("IoError");
    let error: JsValue = error.into();
    set(&error, "domain", "io".into())?;
    set(
        &error,
        "kind",
        JsValue::from_str(&format!("{:?}", source.kind())),
    )?;
    set(
        &error,
        "errno",
        source
            .raw_os_error()
            .map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64)),
    )?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
