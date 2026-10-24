use crate::{
    atoms::Element,
    errors::{error_projection, exception},
};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
error_projection!(PeriodicTableError, ck::PeriodicTableError, {
    #[wasm_bindgen(getter)]
    pub fn element(&self) -> Result<Element, JsValue> {
        match *self.inner {
            ck::PeriodicTableError::UnknownIsotope { element, .. } => Ok(Element(element)),
            _ => {
                Err(js_sys::TypeError::new("unrecognized native PeriodicTableError variant").into())
            }
        }
    }
    #[wasm_bindgen(getter)]
    pub fn isotope(&self) -> Result<u16, JsValue> {
        match *self.inner {
            ck::PeriodicTableError::UnknownIsotope { isotope, .. } => Ok(isotope),
            _ => {
                Err(js_sys::TypeError::new("unrecognized native PeriodicTableError variant").into())
            }
        }
    }
});
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum PeriodicTableErrorKind {
    UnknownIsotope,
}
#[wasm_bindgen]
impl PeriodicTableError {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> Result<PeriodicTableErrorKind, JsValue> {
        match *self.inner {
            ck::PeriodicTableError::UnknownIsotope { .. } => {
                Ok(PeriodicTableErrorKind::UnknownIsotope)
            }
            _ => {
                Err(js_sys::TypeError::new("unrecognized native PeriodicTableError variant").into())
            }
        }
    }
}
