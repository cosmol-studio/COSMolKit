//! Private checked host-value conversion shared by language projections.
use js_sys::{Array, ArrayBuffer};
use wasm_bindgen::prelude::*;
pub(crate) fn text(value: &cosmolkit_wasm::rust::PropertyText) -> Result<String, JsValue> {
    std::str::from_utf8(value.as_bytes())
        .map(str::to_owned)
        .map_err(|error| js_sys::TypeError::new(&format!("binding text encoding: {error}")).into())
}
pub(crate) fn optional_text(
    value: Option<&cosmolkit_wasm::rust::PropertyText>,
) -> Result<JsValue, JsValue> {
    value
        .map(text)
        .transpose()
        .map(|value| value.map_or(JsValue::NULL, JsValue::from))
}
pub(crate) fn type_error(name: &str) -> JsValue {
    js_sys::TypeError::new(&format!("invalid {name}")).into()
}
pub(crate) fn integer(value: &JsValue, name: &str, min: f64, max: f64) -> Result<f64, JsValue> {
    let n = value.as_f64().ok_or_else(|| type_error(name))?;
    if !n.is_finite() || n.fract() != 0.0 || n < min || n > max {
        return Err(
            js_sys::RangeError::new(&format!("{name} is outside its integer range")).into(),
        );
    }
    Ok(n)
}
pub(crate) fn i32_value(value: &JsValue, name: &str) -> Result<i32, JsValue> {
    Ok(integer(value, name, i32::MIN as f64, i32::MAX as f64)? as i32)
}
pub(crate) fn u32_value(value: &JsValue, name: &str) -> Result<u32, JsValue> {
    Ok(integer(value, name, 0.0, u32::MAX as f64)? as u32)
}
pub(crate) fn usize_value(value: &JsValue, name: &str) -> Result<usize, JsValue> {
    Ok(integer(value, name, 0.0, usize::MAX as f64)? as usize)
}
pub(crate) fn bool_value(value: &JsValue, name: &str) -> Result<bool, JsValue> {
    value.as_bool().ok_or_else(|| type_error(name))
}
pub(crate) fn sequence(value: &JsValue, name: &str) -> Result<Array, JsValue> {
    if !Array::is_array(value) && !ArrayBuffer::is_view(value) {
        return Err(type_error(name));
    }
    Ok(Array::from(value))
}
