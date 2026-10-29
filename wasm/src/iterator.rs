//! Standard synchronous iterator result transport; values retain their owner.
use wasm_bindgen::JsValue;
pub(crate) fn result(value: Option<JsValue>) -> Result<JsValue, JsValue> {
    let result = js_sys::Object::new();
    js_sys::Reflect::set(&result, &"done".into(), &value.is_none().into())?;
    js_sys::Reflect::set(
        &result,
        &"value".into(),
        &value.unwrap_or(JsValue::UNDEFINED),
    )?;
    Ok(result.into())
}
