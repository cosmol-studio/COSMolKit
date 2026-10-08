//! Source property strings retain absent values and full owner-defined spelling.
use crate::Molecule;
use crate::alignment_values::set;
use crate::host_values::usize_value;
use crate::property_values::kind;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
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
