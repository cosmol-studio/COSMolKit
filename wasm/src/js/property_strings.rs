//! Source property strings retain absent values and full owner-defined spelling.
use crate::Molecule;
use crate::alignment_values::set;
use crate::host_values::usize_value;
use crate::property_values::kind;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

fn native_property(value: &JsValue) -> Result<ck::PropertyValue, JsValue> {
    if let Some(v) = value.as_bool() { return Ok(ck::PropertyValue::Bool(v)); }
    if let Some(v) = value.as_string() { return Ok(v.into()); }
    if let Some(v) = value.as_f64() {
        if v.is_finite() && v.fract() == 0.0 && !(v == 0.0 && v.is_sign_negative()) {
            if v >= i32::MIN as f64 && v <= i32::MAX as f64 { return Ok(ck::PropertyValue::Int(v as i32)); }
            return crate::host_values::u32_value(value,"atom property integer").map(ck::PropertyValue::UInt);
        }
        return Ok(ck::PropertyValue::Double(v));
    }
    if js_sys::Array::is_array(value) {
        let array = js_sys::Array::from(value);
        if array.iter().all(|v| v.as_string().is_some()) {
            return Ok(ck::PropertyValue::StringVector(array.iter().map(|v| v.as_string().unwrap().into()).collect()));
        }
        return array.iter().map(|v| crate::host_values::i32_value(&v,"atom property list item")).collect::<Result<Vec<_>,_>>().map(ck::PropertyValue::IntVector);
    }
    Err(crate::host_values::type_error("atom property value"))
}
fn property_to_js(value: &ck::PropertyValue) -> Result<JsValue,JsValue> {
    Ok(match value {
        ck::PropertyValue::Bool(v) => (*v).into(),
        ck::PropertyValue::Int(v) => (*v).into(),
        ck::PropertyValue::UInt(v) => (*v).into(),
        ck::PropertyValue::Double(v) => (*v).into(),
        ck::PropertyValue::String(v) => crate::host_values::text(v)?.into(),
        ck::PropertyValue::IntVector(v) => v.iter().map(|v| JsValue::from(*v)).collect::<js_sys::Array>().into(),
        ck::PropertyValue::StringVector(v) => v.iter().map(|s| crate::host_values::text(s).map(JsValue::from)).collect::<Result<js_sys::Array,_>>()?.into(),
    })
}
#[wasm_bindgen]
pub struct PropertyStringError {
    inner: ck::PropertyStringError,
}
#[wasm_bindgen]
impl PropertyStringError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "property".into()
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
    pub fn kind(&self) -> u32 {
        kind(self.inner.kind())
    }
}
pub(crate) fn string_error(source: &ck::PropertyStringError) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("PropertyStringError");
    let e: JsValue = e.into();
    set(&e, "domain", "property".into())?;
    set(&e, "kind", kind(source.kind()).into())?;
    set(&e, "detail", PropertyStringError { inner: *source }.into())?;
    Ok(e)
}
#[wasm_bindgen]
impl Molecule {
    /// Read typed user metadata; a missing key returns null.
    #[wasm_bindgen(js_name=atomProperty,unchecked_return_type="boolean | number | string | number[] | string[] | null")]
    pub fn atom_property(&self, #[wasm_bindgen(unchecked_param_type="number")] atom:JsValue,key:&str)->Result<JsValue,JsValue> {
        let value = self.inner.atom_property(usize_value(&atom,"atom")?,key)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e|e))?;
        value.as_ref().map(property_to_js).transpose().map(|v|v.unwrap_or(JsValue::NULL))
    }
    /// Return a new molecule; the receiver and its shared copies are unchanged.
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=withAtomProperty)]
    pub fn with_atom_property(&self, #[wasm_bindgen(unchecked_param_type="number")] atom:JsValue,key:&str,#[wasm_bindgen(unchecked_param_type="boolean | number | string | number[] | string[]")] value:JsValue)->Result<Molecule,JsValue> {
        self.inner.with_atom_property(usize_value(&atom,"atom")?,key,&native_property(&value)?)
            .map(|inner|Self{inner:std::sync::Arc::new(inner)})
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e|e))
    }
    /// Set metadata in place; COW keeps other molecule copies unchanged.
    #[cfg(feature = "cap-transforms")]
    #[wasm_bindgen(js_name=setAtomProperty)]
    pub fn set_atom_property_(&self, #[wasm_bindgen(unchecked_param_type="number")] atom:JsValue,key:&str,#[wasm_bindgen(unchecked_param_type="boolean | number | string | number[] | string[]")] value:JsValue)->Result<(),JsValue> {
        self.inner.set_atom_property_(usize_value(&atom,"atom")?,key,&native_property(&value)?)
            .map_err(|e| crate::alignment_values::operation_error(&e).unwrap_or_else(|e|e))
    }
    #[wasm_bindgen(js_name=atomPropertyString,unchecked_return_type="string | null")]
    pub fn atom_property_string(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.atom_property_string(id,key)
        self.inner
            .atom_property_string(usize_value(&id, "id")?, key)
            .map_err(|e| string_error(&e).unwrap_or_else(|e| e))
            .and_then(|v| crate::host_values::optional_text(v.as_ref()))
    }
    #[wasm_bindgen(js_name=bondPropertyString,unchecked_return_type="string | null")]
    pub fn bond_property_string(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] id: JsValue,
        key: &str,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.bond_property_string(id,key)
        self.inner
            .bond_property_string(usize_value(&id, "id")?, key)
            .map_err(|e| string_error(&e).unwrap_or_else(|e| e))
            .and_then(|v| crate::host_values::optional_text(v.as_ref()))
    }
}
