//! Complete detached property values used by finalized records.
use crate::alignment_values::set;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PropertyValueKind {
    String,
    Int,
    UInt,
    IntVector,
    Double,
    Bool,
    StringVector,
}
pub(crate) fn kind(value: ck::PropertyValueKind) -> u32 {
    match value {
        ck::PropertyValueKind::String => 0,
        ck::PropertyValueKind::Int => 1,
        ck::PropertyValueKind::UInt => 2,
        ck::PropertyValueKind::IntVector => 3,
        ck::PropertyValueKind::Double => 4,
        ck::PropertyValueKind::Bool => 5,
        ck::PropertyValueKind::StringVector => 6,
    }
}
#[wasm_bindgen]
pub struct PropertyValue {
    pub(crate) inner: ck::PropertyValue,
}
#[wasm_bindgen]
impl PropertyValue {
    #[wasm_bindgen(unchecked_return_type = "PropertyValueKind")]
    pub fn kind(&self) -> u32 {
        kind(self.inner.kind())
    }
    #[wasm_bindgen(js_name=asString)]
    pub fn as_string(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.as_string()
        self.inner
            .as_string()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
            .and_then(crate::host_values::text)
    }
    #[wasm_bindgen(js_name=asInt)]
    pub fn as_int(&self) -> Result<i32, JsValue> {
        // COSMolKit❗✔️: self.inner.as_int()
        self.inner
            .as_int()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=asUint)]
    pub fn as_uint(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: self.inner.as_uint()
        self.inner
            .as_uint()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=asIntVector,unchecked_return_type="number[]")]
    pub fn as_int_vector(&self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.as_int_vector()
        self.inner
            .as_int_vector()
            .map(|v| {
                v.iter()
                    .map(|v| JsValue::from(*v))
                    .collect::<js_sys::Array>()
                    .into()
            })
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=asDouble)]
    pub fn as_double(&self) -> Result<f64, JsValue> {
        // COSMolKit❗✔️: self.inner.as_double()
        self.inner
            .as_double()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=asBool)]
    pub fn as_bool(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.as_bool()
        self.inner
            .as_bool()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=asStringVector,unchecked_return_type="string[]")]
    pub fn as_string_vector(&self) -> Result<JsValue, JsValue> {
        self.inner
            .as_string_vector()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))?
            .iter()
            .map(|value| crate::host_values::text(value).map(JsValue::from))
            .collect::<Result<js_sys::Array, JsValue>>()
            .map(Into::into)
    }
}
#[wasm_bindgen]
pub struct PropertyValueError {
    inner: ck::PropertyValueError,
}
#[wasm_bindgen]
impl PropertyValueError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "property".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        "KindMismatch".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        JsValue::NULL
    }
    #[wasm_bindgen(unchecked_return_type = "PropertyValueKind")]
    pub fn expected(&self) -> u32 {
        kind(self.inner.expected())
    }
    #[wasm_bindgen(unchecked_return_type = "PropertyValueKind")]
    pub fn actual(&self) -> u32 {
        kind(self.inner.actual())
    }
}
pub(crate) fn property_error(source: &ck::PropertyValueError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("PropertyValueError");
    let e: JsValue = error.into();
    set(&e, "domain", "property".into())?;
    set(&e, "kind", "KindMismatch".into())?;
    set(&e, "expected", kind(source.expected()).into())?;
    set(&e, "actual", kind(source.actual()).into())?;
    set(
        &e,
        "detail",
        PropertyValueError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum SdfPropertyListTarget {
    Atom,
    Bond,
}
#[wasm_bindgen]
pub struct SdfPropertyList {
    pub(crate) inner: ck::SdfPropertyList,
}
#[wasm_bindgen]
impl SdfPropertyList {
    #[wasm_bindgen(unchecked_return_type = "SdfPropertyListTarget")]
    pub fn target(&self) -> u32 {
        match self.inner.target() {
            ck::SdfPropertyListTarget::Atom => 0,
            ck::SdfPropertyListTarget::Bond => 1,
        }
    }
    pub fn name(&self) -> Result<String, JsValue> {
        crate::host_values::text(self.inner.name())
    }
    #[wasm_bindgen(unchecked_return_type = "(PropertyValue | null)[]")]
    pub fn values(&self) -> JsValue {
        self.inner
            .values()
            .iter()
            .map(|v| {
                v.as_ref().map_or(JsValue::NULL, |inner| {
                    PropertyValue {
                        inner: inner.clone(),
                    }
                    .into()
                })
            })
            .collect::<js_sys::Array>()
            .into()
    }
}
#[wasm_bindgen]
pub struct MoleculeProperties {
    pub(crate) inner: ck::MoleculeProperties,
}
pub(crate) fn fields(values: &[(ck::PropertyText, ck::PropertyText)]) -> Result<JsValue, JsValue> {
    values
        .iter()
        .map(|(k, v)| {
            let a = js_sys::Array::new();
            a.push(&crate::host_values::text(k)?.into());
            a.push(&crate::host_values::text(v)?.into());
            Ok(JsValue::from(a))
        })
        .collect::<Result<js_sys::Array, JsValue>>()
        .map(Into::into)
}
#[wasm_bindgen]
impl MoleculeProperties {
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn name(&self) -> Result<JsValue, JsValue> {
        crate::host_values::optional_text(self.inner.name())
    }
    #[wasm_bindgen(unchecked_return_type = "Map<string, PropertyValue>")]
    pub fn props(&self) -> Result<js_sys::Map, JsValue> {
        let m = js_sys::Map::new();
        for (k, v) in self.inner.ordered_props() {
            m.set(
                &crate::host_values::text(k)?.into(),
                &PropertyValue { inner: v.clone() }.into(),
            );
        }
        Ok(m)
    }
    #[wasm_bindgen(unchecked_return_type = "PropertyValue | null")]
    pub fn prop(&self, key: &str) -> JsValue {
        self.inner
            .prop(key)
            .map_or(JsValue::NULL, |v| PropertyValue { inner: v.clone() }.into())
    }
    #[wasm_bindgen(js_name=isPropComputed)]
    pub fn is_prop_computed(&self, key: &str) -> Result<bool, JsValue> {
        self.inner
            .is_prop_computed(key)
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=computedPropNames,unchecked_return_type="string[]")]
    pub fn computed_prop_names(&self) -> Result<JsValue, JsValue> {
        let names = self
            .inner
            .computed_prop_names()
            .map_err(|e| property_error(&e).unwrap_or_else(|e| e))?;
        names
            .unwrap_or_default()
            .iter()
            .map(|v| crate::host_values::text(v).map(JsValue::from))
            .collect::<Result<js_sys::Array, JsValue>>()
            .map(Into::into)
    }
    #[wasm_bindgen(js_name=sdfDataFields,unchecked_return_type="[string, string][]")]
    pub fn sdf_data_fields(&self) -> Result<JsValue, JsValue> {
        fields(self.inner.sdf_data_fields())
    }
    #[wasm_bindgen(js_name=sdfPropertyLists,unchecked_return_type="SdfPropertyList[]")]
    pub fn sdf_property_lists(&self) -> JsValue {
        self.inner
            .sdf_property_lists()
            .iter()
            .cloned()
            .map(|inner| JsValue::from(SdfPropertyList { inner }))
            .collect::<js_sys::Array>()
            .into()
    }
}
