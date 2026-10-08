//! InChI language transport; algorithms remain in the public Rust facade.
use crate::{
    Molecule,
    alignment_values::set,
    host_values::{bool_value, type_error},
};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum InchiErrorKind {
    AllocationFailed,
    UnsupportedState,
    InvalidInput,
    InvalidSourceOutput,
    SanitizeFailed,
    Toolkit,
    SourcePort,
}

fn kind(value: ck::InchiErrorKind) -> InchiErrorKind {
    match value {
        ck::InchiErrorKind::AllocationFailed => InchiErrorKind::AllocationFailed,
        ck::InchiErrorKind::UnsupportedState => InchiErrorKind::UnsupportedState,
        ck::InchiErrorKind::InvalidInput => InchiErrorKind::InvalidInput,
        ck::InchiErrorKind::InvalidSourceOutput => InchiErrorKind::InvalidSourceOutput,
        ck::InchiErrorKind::SanitizeFailed => InchiErrorKind::SanitizeFailed,
        ck::InchiErrorKind::Toolkit => InchiErrorKind::Toolkit,
        ck::InchiErrorKind::SourcePort => InchiErrorKind::SourcePort,
    }
}

#[wasm_bindgen]
pub struct InchiError {
    inner: ck::InchiError,
}
#[wasm_bindgen]
impl InchiError {
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> InchiErrorKind {
        kind(self.inner.kind)
    }
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "inchi".into()
    }
    #[wasm_bindgen(getter)]
    pub fn operation(&self) -> String {
        self.inner.operation.into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn detail(&self) -> String {
        self.inner.detail.clone()
    }
}

fn error(source: ck::InchiError) -> JsValue {
    let kind = kind(source.kind);
    let exception = js_sys::Error::new(&source.to_string());
    exception.set_name("InchiError");
    let value: JsValue = exception.into();
    let fields = || -> Result<(), JsValue> {
        set(&value, "domain", "inchi".into())?;
        set(&value, "kind", (kind as u32).into())?;
        set(&value, "operation", source.operation.into())?;
        set(&value, "detail", InchiError { inner: source }.into())?;
        Ok(())
    };
    match fields() {
        Ok(()) => value,
        Err(error) => error,
    }
}

#[wasm_bindgen]
pub struct InchiReadParams {
    inner: ck::InchiReadParams,
}
#[wasm_bindgen]
impl InchiReadParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hydrogens: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::InchiReadParams {
                sanitize: if sanitize.is_undefined() {
                    true
                } else {
                    bool_value(&sanitize, "sanitize")?
                },
                remove_hydrogens: if remove_hydrogens.is_undefined() {
                    true
                } else {
                    bool_value(&remove_hydrogens, "removeHydrogens")?
                },
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[wasm_bindgen(setter, js_name=sanitize)]
    pub fn set_sanitize(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.sanitize = bool_value(&value, "sanitize")?;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=removeHydrogens)]
    pub fn remove_hydrogens(&self) -> bool {
        self.inner.remove_hydrogens
    }
    #[wasm_bindgen(setter, js_name=removeHydrogens)]
    pub fn set_remove_hydrogens(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.remove_hydrogens = bool_value(&value, "removeHydrogens")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct InchiWriteParams {
    inner: ck::InchiWriteParams,
}
#[wasm_bindgen]
impl InchiWriteParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "string")] options: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::InchiWriteParams {
                options: if options.is_undefined() {
                    String::new()
                } else {
                    options.as_string().ok_or_else(|| type_error("options"))?
                },
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn options(&self) -> String {
        self.inner.options.clone()
    }
    #[wasm_bindgen(setter, js_name=options)]
    pub fn set_options(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "string")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.options = value.as_string().ok_or_else(|| type_error("options"))?;
        Ok(())
    }
}

// Borrow parameter instances without consuming or freeing reusable objects.
#[wasm_bindgen(
    inline_js = "export function visitInchiRead(v,f){f(v);} export function visitInchiWrite(v,f){f(v);}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitInchiRead)]
    fn visit_read(
        value: &JsValue,
        callback: &mut dyn FnMut(&InchiReadParams),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch, js_name=visitInchiWrite)]
    fn visit_write(
        value: &JsValue,
        callback: &mut dyn FnMut(&InchiWriteParams),
    ) -> Result<(), JsValue>;
}

fn options(value: &JsValue, allowed: &[&str]) -> Result<(), JsValue> {
    if !value.is_object() || value.is_null() || js_sys::Array::is_array(value) {
        return Err(type_error("InChI options"));
    }
    for key in js_sys::Object::keys(&js_sys::Object::from(value.clone())) {
        if !key
            .as_string()
            .is_some_and(|k| allowed.contains(&k.as_str()))
        {
            return Err(type_error("unknown InChI option"));
        }
    }
    Ok(())
}
fn read_params(value: &JsValue) -> Result<ck::InchiReadParams, JsValue> {
    if value.is_undefined() {
        return Ok(Default::default());
    }
    let mut result = None;
    if visit_read(value, &mut |p: &InchiReadParams| result = Some(p.inner)).is_ok()
        && let Some(result) = result
    {
        return Ok(result);
    }
    options(value, &["sanitize", "removeHydrogens"])?;
    let mut result = ck::InchiReadParams::default();
    let sanitize = js_sys::Reflect::get(value, &"sanitize".into())?;
    let remove = js_sys::Reflect::get(value, &"removeHydrogens".into())?;
    if !sanitize.is_undefined() {
        result.sanitize = bool_value(&sanitize, "sanitize")?;
    }
    if !remove.is_undefined() {
        result.remove_hydrogens = bool_value(&remove, "removeHydrogens")?;
    }
    Ok(result)
}
fn write_params(value: &JsValue) -> Result<ck::InchiWriteParams, JsValue> {
    if value.is_undefined() {
        return Ok(Default::default());
    }
    let mut result = None;
    if visit_write(value, &mut |p: &InchiWriteParams| {
        result = Some(p.inner.clone())
    })
    .is_ok()
        && let Some(result) = result
    {
        return Ok(result);
    }
    options(value, &["options"])?;
    let text = js_sys::Reflect::get(value, &"options".into())?;
    Ok(ck::InchiWriteParams {
        options: if text.is_undefined() {
            String::new()
        } else {
            text.as_string().ok_or_else(|| type_error("options"))?
        },
    })
}

#[wasm_bindgen(js_name=inchiToInchiKey)]
pub fn inchi_to_inchi_key(text: &str) -> Result<String, JsValue> {
    ck::inchi_to_inchi_key(text).map_err(error)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fromInchi)]
    pub fn from_inchi(
        text: &str,
        #[wasm_bindgen(
            unchecked_optional_param_type = "InchiReadParams | { sanitize?: boolean; removeHydrogens?: boolean }"
        )]
        params: JsValue,
    ) -> Result<Self, JsValue> {
        cosmolkit_wasm::Molecule::from_inchi_with_params(text, &read_params(&params)?)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fromInchiWithParams)]
    pub fn from_inchi_with_params(text: &str, params: &InchiReadParams) -> Result<Self, JsValue> {
        cosmolkit_wasm::Molecule::from_inchi_with_params(text, &params.inner)
            .map(|inner| Self {
                inner: inner.into(),
            })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=toInchi)]
    pub fn to_inchi(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "InchiWriteParams | { options?: string }")]
        params: JsValue,
    ) -> Result<String, JsValue> {
        self.inner
            .to_inchi_with_params(&write_params(&params)?)
            .map_err(error)
    }
    #[wasm_bindgen(js_name=toInchiWithParams)]
    pub fn to_inchi_with_params(&self, params: &InchiWriteParams) -> Result<String, JsValue> {
        self.inner
            .to_inchi_with_params(&params.inner)
            .map_err(error)
    }
    #[wasm_bindgen(js_name=toInchiKey)]
    pub fn to_inchi_key(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "InchiWriteParams | { options?: string }")]
        params: JsValue,
    ) -> Result<String, JsValue> {
        self.inner
            .to_inchi_key_with_params(&write_params(&params)?)
            .map_err(error)
    }
    #[wasm_bindgen(js_name=toInchiKeyWithParams)]
    pub fn to_inchi_key_with_params(&self, params: &InchiWriteParams) -> Result<String, JsValue> {
        self.inner
            .to_inchi_key_with_params(&params.inner)
            .map_err(error)
    }
}
