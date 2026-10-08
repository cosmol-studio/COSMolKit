//! Source-defined XYZ and MOL2 failures retain their fields and recursive causes.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct XyzReadError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl XyzReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
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
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"value".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=line,unchecked_return_type="number | null")]
    pub fn field_1(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"line".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn xyz_read_error(source: &ck::XyzReadError) -> Result<JsValue, JsValue> {
    let fields = js_sys::Object::new();
    set(&fields, "value", JsValue::NULL)?;
    set(&fields, "line", JsValue::NULL)?;
    use ck::XyzReadError as E;
    let kind = match source {
        E::EmptyBlock => "EmptyBlock",
        E::UnexpectedEof => "UnexpectedEof",
        E::AtomCount { value } => {
            set(&fields, "value", value.into())?;
            "AtomCount"
        }
        E::MissingCoordinates { line } => {
            set(&fields, "line", (*line as u32).into())?;
            "MissingCoordinates"
        }
        E::Coordinate { value, line, .. } => {
            set(&fields, "value", value.into())?;
            set(&fields, "line", (*line as u32).into())?;
            "Coordinate"
        }
        E::AtomSymbol { .. } => "AtomSymbol",
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
        E::MoleculeProperty(_) => "MoleculeProperty",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("XyzReadError");
    let e: JsValue = e.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let v = js_sys::Reflect::get(&fields, &key)?;
        if !v.is_null() {
            set(&e, &key.as_string().unwrap(), v)?;
        }
    }
    set(
        &e,
        "detail",
        XyzReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct XyzWriteError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl XyzWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
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
    #[wasm_bindgen(getter,js_name=id,unchecked_return_type="number | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"id".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn xyz_write_error(source: &ck::XyzWriteError) -> Result<JsValue, JsValue> {
    let fields = js_sys::Object::new();
    set(&fields, "id", JsValue::NULL)?;
    use ck::XyzWriteError as E;
    let kind = match source {
        E::Property(..) => "Property",
        E::ConformerNotFound { id } => {
            set(&fields, "id", (*id as u32).into())?;
            "ConformerNotFound"
        }
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("XyzWriteError");
    let e: JsValue = e.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let v = js_sys::Reflect::get(&fields, &key)?;
        if !v.is_null() {
            set(&e, &key.as_string().unwrap(), v)?;
        }
    }
    set(
        &e,
        "detail",
        XyzWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct Mol2ReadError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl Mol2ReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
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
    #[wasm_bindgen(getter,js_name=feature,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"feature".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=parseDetail,unchecked_return_type="string | null")]
    pub fn field_1(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"parseDetail".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn mol2_read_error(source: &ck::Mol2ReadError) -> Result<JsValue, JsValue> {
    let fields = js_sys::Object::new();
    set(&fields, "feature", JsValue::NULL)?;
    set(&fields, "parseDetail", JsValue::NULL)?;
    use ck::Mol2ReadError as E;
    let kind = match source {
        E::PropertyString(..) => "PropertyString",
        E::Parse(value) => {
            set(&fields, "parseDetail", value.into())?;
            "Parse"
        }
        E::Unsupported { feature } => {
            set(&fields, "feature", (*feature).into())?;
            "Unsupported"
        }
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
        E::AtomProperty(_) => "AtomProperty",
        E::MoleculeProperty(_) => "MoleculeProperty",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("Mol2ReadError");
    let e: JsValue = e.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let v = js_sys::Reflect::get(&fields, &key)?;
        if !v.is_null() {
            set(&e, &key.as_string().unwrap(), v)?;
        }
    }
    set(
        &e,
        "detail",
        Mol2ReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct Mol2PostError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl Mol2PostError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "io".into()
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
    #[wasm_bindgen(getter,js_name=stage,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"stage".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn mol2_post_error(source: &ck::Mol2PostError) -> Result<JsValue, JsValue> {
    let fields = js_sys::Object::new();
    set(&fields, "stage", JsValue::NULL)?;
    let kind = source.stage;
    set(&fields, "stage", source.stage.into())?;
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("Mol2PostError");
    let e: JsValue = e.into();
    set(&e, "domain", "io".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    for key in js_sys::Object::keys(&fields).iter() {
        let v = js_sys::Reflect::get(&fields, &key)?;
        if !v.is_null() {
            set(&e, &key.as_string().unwrap(), v)?;
        }
    }
    set(
        &e,
        "detail",
        Mol2PostError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
