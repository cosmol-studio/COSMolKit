//! Complete canonical PickleError variants and Python payload semantics.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct PickleError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl PickleError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "serialization".into()
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
    #[wasm_bindgen(getter,js_name=version,unchecked_return_type="number | null")]
    pub fn field_0(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"version".into())
    }
    #[wasm_bindgen(getter,js_name=major,unchecked_return_type="number | null")]
    pub fn field_1(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"major".into())
    }
    #[wasm_bindgen(getter,js_name=minor,unchecked_return_type="number | null")]
    pub fn field_2(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"minor".into())
    }
    #[wasm_bindgen(getter,js_name=section,unchecked_return_type="number | null")]
    pub fn field_3(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"section".into())
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="number | null")]
    pub fn field_4(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"expected".into())
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn field_5(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"actual".into())
    }
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="number | null")]
    pub fn field_6(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"value".into())
    }
    #[wasm_bindgen(getter,js_name=typeName,unchecked_return_type="string | null")]
    pub fn field_7(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"typeName".into())
    }
    #[wasm_bindgen(getter,js_name=count,unchecked_return_type="number | null")]
    pub fn field_8(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"count".into())
    }
}
pub(crate) fn pickle_error(source: &ck::PickleError) -> Result<JsValue, JsValue> {
    use ck::PickleError as E;
    let fields = js_sys::Object::new();
    set(&fields, "version", JsValue::NULL)?;
    set(&fields, "major", JsValue::NULL)?;
    set(&fields, "minor", JsValue::NULL)?;
    set(&fields, "section", JsValue::NULL)?;
    set(&fields, "expected", JsValue::NULL)?;
    set(&fields, "actual", JsValue::NULL)?;
    set(&fields, "value", JsValue::NULL)?;
    set(&fields, "typeName", JsValue::NULL)?;
    set(&fields, "count", JsValue::NULL)?;
    let mut message = source.to_string();
    let kind = match source {
        E::UnexpectedEof => "UnexpectedEof",
        E::UnsupportedVersion(version) => {
            set(&fields, "version", (*version).into())?;
            "UnsupportedVersion"
        }
        E::UnsupportedArchiveVersion { major, minor } => {
            set(&fields, "major", (*major).into())?;
            set(&fields, "minor", (*minor).into())?;
            "UnsupportedArchiveVersion"
        }
        E::UnsupportedSectionVersion { section, version } => {
            set(&fields, "section", (*section).into())?;
            set(&fields, "version", (*version).into())?;
            "UnsupportedSectionVersion"
        }
        E::MissingRequiredSection(section) => {
            set(&fields, "section", (*section).into())?;
            "MissingRequiredSection"
        }
        E::DuplicateSection(section) => {
            set(&fields, "section", (*section).into())?;
            "DuplicateSection"
        }
        E::UnknownRequiredSection(section) => {
            set(&fields, "section", (*section).into())?;
            "UnknownRequiredSection"
        }
        E::DataLengthMismatch { expected, actual } => {
            set(&fields, "expected", (*expected as u32).into())?;
            set(&fields, "actual", (*actual as u32).into())?;
            "DataLengthMismatch"
        }
        E::InvalidEnumValue { value, type_name } => {
            set(&fields, "value", (*value).into())?;
            set(&fields, "typeName", (*type_name).into())?;
            "InvalidEnumValue"
        }
        E::InvalidArchive(value) => {
            message = value.clone();
            "InvalidArchive"
        }
        E::InvalidMolecule(value) => {
            message = value.clone();
            "InvalidMolecule"
        }
        E::TooManyAtoms(count) => {
            set(&fields, "count", (*count as u32).into())?;
            "TooManyAtoms"
        }
        E::TooManyBonds(count) => {
            set(&fields, "count", (*count as u32).into())?;
            "TooManyBonds"
        }
        E::StringTooLong(count) => {
            set(&fields, "count", (*count as u32).into())?;
            "StringTooLong"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&message);
    error.set_name("PickleError");
    set(&error, "domain", "serialization".into())?;
    set(&error, "kind", kind.into())?;
    set(&error, "cause", cause.clone())?;
    for key in js_sys::Object::keys(&fields).iter() {
        let key = key.as_string().expect("own field key");
        set(
            &error,
            &key,
            js_sys::Reflect::get(&fields, &key.clone().into())?,
        )?;
    }
    set(
        &error,
        "detail",
        PickleError {
            kind: kind.into(),
            message,
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(error.into())
}
