//! Typed XYZ/MOL2 parameters and thin canonical molecule transport calls.
use crate::Molecule;
use crate::host_values::{bool_value, type_error, u32_value, usize_value};
use crate::io_errors::io_error;
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum Mol2Type {
    Corina,
}
#[wasm_bindgen]
pub struct Mol2ReadParams {
    inner: ck::Mol2ReadParams,
}
#[wasm_bindgen]
impl Mol2ReadParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "Mol2Type")] variant: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] cleanup_substructures: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::Mol2ReadParams::default();
        if !sanitize.is_undefined() {
            inner.sanitize = bool_value(&sanitize, "sanitize")?;
        }
        if !remove_hs.is_undefined() {
            inner.remove_hs = bool_value(&remove_hs, "removeHs")?;
        }
        if !cleanup_substructures.is_undefined() {
            inner.cleanup_substructures =
                bool_value(&cleanup_substructures, "cleanupSubstructures")?;
        }
        if !variant.is_undefined() {
            inner.variant = match u32_value(&variant, "variant")? {
                0 => ck::Mol2Type::Corina,
                _ => return Err(js_sys::RangeError::new("invalid Mol2Type").into()),
            };
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[wasm_bindgen(getter,js_name=removeHs)]
    pub fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Mol2Type")]
    pub fn variant(&self) -> u32 {
        match self.inner.variant {
            ck::Mol2Type::Corina => 0,
        }
    }
    #[wasm_bindgen(getter,js_name=cleanupSubstructures)]
    pub fn cleanup_substructures(&self) -> bool {
        self.inner.cleanup_substructures
    }
}
#[wasm_bindgen]
pub struct XyzWriteParams {
    inner: ck::XyzWriteParams,
}
#[wasm_bindgen]
impl XyzWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] precision: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::XyzWriteParams::default();
        if !conformer_id.is_undefined() {
            inner.conformer_id = if conformer_id.is_null() {
                None
            } else {
                Some(usize_value(&conformer_id, "conformerId")?)
            };
        }
        if !precision.is_undefined() {
            inner.precision = u32_value(&precision, "precision")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id
            .map_or(JsValue::NULL, |id| JsValue::from(id as u32))
    }
    #[wasm_bindgen(getter)]
    pub fn precision(&self) -> u32 {
        self.inner.precision
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMol2Params(v,f){try{f(v);}catch(cause){throw new TypeError('invalid Mol2ReadParams',{cause});}} export function visitXyzParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid XyzWriteParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitMol2Params)]
    fn visit_mol2(value: &JsValue, visit: &mut dyn FnMut(&Mol2ReadParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitXyzParams)]
    fn visit_xyz(value: &JsValue, visit: &mut dyn FnMut(&XyzWriteParams)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromXyzBlock)]
    pub fn from_xyz_block(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_xyz_block(text)
        cosmolkit_wasm::Molecule::from_xyz_block(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readXyz)]
    pub fn read_xyz(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_xyz(text)
        cosmolkit_wasm::Molecule::read_xyz(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromMol2)]
    pub fn from_mol2(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_mol2(text)
        cosmolkit_wasm::Molecule::from_mol2(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=readMol2)]
    pub fn read_mol2(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_mol2(text)
        cosmolkit_wasm::Molecule::read_mol2(text)
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fromMol2WithParams)]
    pub fn from_mol2_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "Mol2ReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::from_mol2_with_params(text,&params.inner)
        let mut result = None;
        visit_mol2(&params, &mut |p: &Mol2ReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::from_mol2_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=readMol2WithParams)]
    pub fn read_mol2_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "Mol2ReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::Molecule::read_mol2_with_params(text,&params.inner)
        let mut result = None;
        visit_mol2(&params, &mut |p: &Mol2ReadParams| {
            result = Some(
                cosmolkit_wasm::Molecule::read_mol2_with_params(text, &p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=toXyz)]
    pub fn to_xyz(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_xyz()
        self.inner
            .to_xyz()
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeXyz)]
    pub fn write_xyz(&self, path: &str) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_xyz(path)
        self.inner
            .write_xyz(path)
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=toXyzWithParams)]
    pub fn to_xyz_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "XyzWriteParams")] params: JsValue,
    ) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_xyz_with_params(&p.inner)
        let mut result = None;
        visit_xyz(&params, &mut |p: &XyzWriteParams| {
            result = Some(
                self.inner
                    .to_xyz_with_params(&p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=writeXyzWithParams)]
    pub fn write_xyz_with_params(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "XyzWriteParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.write_xyz_with_params(path,&p.inner)
        let mut result = None;
        visit_xyz(&params, &mut |p: &XyzWriteParams| {
            result = Some(
                self.inner
                    .write_xyz_with_params(path, &p.inner)
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
