//! Complete canonical MolWriteError vocabulary and payloads.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MolWriteError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl MolWriteError {
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
    #[wasm_bindgen(getter,js_name=subset,unchecked_return_type="string | null")]
    pub fn field_0(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"subset".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=valueDetail,unchecked_return_type="string | null")]
    pub fn field_1(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"valueDetail".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=twoD,unchecked_return_type="number | null")]
    pub fn field_2(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"twoD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=threeD,unchecked_return_type="number | null")]
    pub fn field_3(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"threeD".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=dimension,unchecked_return_type="CoordinateDimension | null")]
    pub fn field_4(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"dimension".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=id,unchecked_return_type="number | null")]
    pub fn field_5(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"id".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn capability(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"capability".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn mol_write_error(source: &ck::MolWriteError) -> Result<JsValue, JsValue> {
    let fields = js_sys::Object::new();
    set(&fields, "subset", JsValue::NULL)?;
    set(&fields, "valueDetail", JsValue::NULL)?;
    set(&fields, "twoD", JsValue::NULL)?;
    set(&fields, "threeD", JsValue::NULL)?;
    set(&fields, "dimension", JsValue::NULL)?;
    set(&fields, "id", JsValue::NULL)?;
    use ck::MolWriteError as E;
    let kind = match source {
        E::UnsignedProperty(..) => "UnsignedProperty",
        E::UnsupportedSubset(v) => {
            set(&fields, "subset", (*v).into())?;
            "UnsupportedSubset"
        }
        E::Value(v) => {
            set(&fields, "valueDetail", v.into())?;
            "Value"
        }
        E::AmbiguousCoordinates { two_d, three_d } => {
            set(&fields, "twoD", (*two_d as u32).into())?;
            set(&fields, "threeD", (*three_d as u32).into())?;
            "AmbiguousCoordinates"
        }
        E::MissingCoordinate { dimension, id } => {
            set(
                &fields,
                "dimension",
                JsValue::from(match dimension {
                    ck::CoordinateDimension::TwoD => 0u32,
                    ck::CoordinateDimension::ThreeD => 1u32,
                }),
            )?;
            set(&fields, "id", (*id as u32).into())?;
            "MissingCoordinate"
        }
        E::CoordinateDimensionMismatch => "CoordinateDimensionMismatch",
        E::QueryGraph(_) => "QueryGraph",
        E::QueryAtom(_) => "QueryAtom",
        E::QueryState(_) => "QueryState",
        #[cfg(feature = "cap-search")]
        E::QuerySmarts(_) => "QuerySmarts",
        E::Valence(_) => "Valence",
        E::Kekulize(_) => "Kekulize",
        E::Atropisomer(_) => "Atropisomer",
        E::Wedge(_) => "Wedge",
        // Avalon enables the writer's internal coordinate generator even when
        // public depiction APIs are disabled; its concrete errors still exist.
        #[cfg(any(feature = "cap-depict", feature = "cap-fingerprints"))]
        E::Depict(_) => "Depict",
        E::MissingCapability(capability) => {
            set(&fields, "capability", (*capability).into())?;
            "MissingCapability"
        }
        E::Topology(_) => "Topology",
        E::Coordinates(_) => "Coordinates",
        E::Property(_) => "Property",
        E::Detached(_) => "Detached",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("MolWriteError");
    let e: JsValue = error.into();
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
        MolWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
