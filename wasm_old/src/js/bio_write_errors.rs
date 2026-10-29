//! Typed BIO writer variants, original paths and complete error cause chains.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
fn field(context: &JsValue, name: &str) -> JsValue {
    js_sys::Reflect::get(context, &name.into()).unwrap_or(JsValue::NULL)
}
#[wasm_bindgen]
pub struct BioPdbWriteError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioPdbWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
    pub fn chain(&self) -> JsValue {
        field(&self.context, "chain")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn length(&self) -> JsValue {
        field(&self.context, "length")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn serial(&self) -> JsValue {
        field(&self.context, "serial")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn from(&self) -> JsValue {
        field(&self.context, "from")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn value(&self) -> JsValue {
        field(&self.context, "value")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn field(&self) -> JsValue {
        field(&self.context, "field")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn detail(&self) -> JsValue {
        field(&self.context, "detail")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn path(&self) -> JsValue {
        field(&self.context, "path")
    }
}
pub(crate) fn pdb_error(source: &ck::BioPdbWriteError) -> Result<JsValue, JsValue> {
    let context: JsValue = js_sys::Object::new().into();
    set(&context, "chain", JsValue::NULL)?;
    set(&context, "length", JsValue::NULL)?;
    set(&context, "serial", JsValue::NULL)?;
    set(&context, "from", JsValue::NULL)?;
    set(&context, "value", JsValue::NULL)?;
    set(&context, "field", JsValue::NULL)?;
    set(&context, "detail", JsValue::NULL)?;
    set(&context, "path", JsValue::NULL)?;
    let kind = match source {
        ck::BioPdbWriteError::InvalidStructure(_) => "InvalidStructure",
        ck::BioPdbWriteError::ChainNameTooLong { chain, length } => {
            set(&context, "chain", chain.as_str().into())?;
            set(&context, "length", (*length as f64).into())?;
            "ChainNameTooLong"
        }
        ck::BioPdbWriteError::NegativeSerial { serial } => {
            set(&context, "serial", (*serial).into())?;
            "NegativeSerial"
        }
        ck::BioPdbWriteError::SerialIncrementOverflow { from } => {
            set(&context, "from", (*from).into())?;
            "SerialIncrementOverflow"
        }
        ck::BioPdbWriteError::SerialOffsetOverflow { serial } => {
            set(&context, "serial", (*serial).into())?;
            "SerialOffsetOverflow"
        }
        ck::BioPdbWriteError::NegativeBase36Value { value } => {
            set(&context, "value", (*value).into())?;
            "NegativeBase36Value"
        }
        ck::BioPdbWriteError::UnrepresentableField { field, detail } => {
            set(&context, "field", (*field).into())?;
            set(&context, "detail", detail.as_str().into())?;
            "UnrepresentableField"
        }
        ck::BioPdbWriteError::Io { path, .. } => {
            set(&context, "path", path.to_string_lossy().as_ref().into())?;
            "Io"
        }
    };
    // COSMolKit❗✔️: annotate(py, error, kind, &source)
    // Source-defined categories and causal chain; one scalar/context copy per error.
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioPdbWriteError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(&e, "chain", field(&context, "chain"))?;
    set(&e, "length", field(&context, "length"))?;
    set(&e, "serial", field(&context, "serial"))?;
    set(&e, "from", field(&context, "from"))?;
    set(&e, "value", field(&context, "value"))?;
    set(&e, "field", field(&context, "field"))?;
    set(&e, "path", field(&context, "path"))?;
    set(
        &e,
        "detail",
        BioPdbWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct BioMmcifWriteError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioMmcifWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
    pub fn field(&self) -> JsValue {
        field(&self.context, "field")
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn path(&self) -> JsValue {
        field(&self.context, "path")
    }
}
pub(crate) fn mmcif_error(source: &ck::BioMmcifWriteError) -> Result<JsValue, JsValue> {
    let context: JsValue = js_sys::Object::new().into();
    set(&context, "field", JsValue::NULL)?;
    set(&context, "path", JsValue::NULL)?;
    let kind = match source {
        ck::BioMmcifWriteError::Structure(_) => "Structure",
        ck::BioMmcifWriteError::Cif(_) => "Cif",
        ck::BioMmcifWriteError::InvalidText { field } => {
            set(&context, "field", (*field).into())?;
            "InvalidText"
        }
        ck::BioMmcifWriteError::Io(_) => "Io",
        ck::BioMmcifWriteError::FileWrite { path, .. } => {
            set(&context, "path", path.to_string_lossy().as_ref().into())?;
            "FileWrite"
        }
    };
    // COSMolKit❗✔️: annotate(py, error, kind, &source)
    // Source-defined categories and causal chain; one scalar/context copy per error.
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("BioMmcifWriteError");
    let e: JsValue = e.into();
    set(&e, "domain", "bio".into())?;
    set(&e, "kind", kind.into())?;
    set(&e, "field", field(&context, "field"))?;
    set(&e, "path", field(&context, "path"))?;
    set(
        &e,
        "detail",
        BioMmcifWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
