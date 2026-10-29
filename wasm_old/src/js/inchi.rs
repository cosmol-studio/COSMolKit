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
#[cosmolkit_wasm::javascript_options("InchiReadParamsOptions")]
#[wasm_bindgen]
impl InchiReadParams {
    #[wasm_bindgen(js_name=setSanitize)]
    pub fn set_sanitize_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_sanitize(value)
    }
    #[wasm_bindgen(js_name=setRemoveHs)]
    pub fn set_remove_hs_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_remove_hs(value)
    }
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] sanitize: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] remove_hs: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::InchiReadParams {
                sanitize: if sanitize.is_undefined() {
                    true
                } else {
                    bool_value(&sanitize, "sanitize")?
                },
                remove_hs: if remove_hs.is_undefined() {
                    true
                } else {
                    bool_value(&remove_hs, "removeHs")?
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
    #[wasm_bindgen(getter, js_name=removeHs)]
    pub fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[wasm_bindgen(setter, js_name=removeHs)]
    pub fn set_remove_hs(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner.remove_hs = bool_value(&value, "removeHs")?;
        Ok(())
    }
}

#[wasm_bindgen]
pub struct InchiWriteParams {
    inner: ck::InchiWriteParams,
}
#[cosmolkit_wasm::javascript_options("InchiWriteParamsOptions")]
#[wasm_bindgen]
impl InchiWriteParams {
    #[wasm_bindgen(js_name=setOptions)]
    pub fn set_options_method(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "string")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_options(value)
    }
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

// Both call forms use the same constructor-derived configuration boundary.
fn read_params(value: &JsValue) -> Result<ck::InchiReadParams, JsValue> {
    Ok(InchiReadParams::from_js_configuration(value)?.inner)
}
fn write_params(value: &JsValue) -> Result<ck::InchiWriteParams, JsValue> {
    Ok(InchiWriteParams::from_js_configuration(value)?.inner)
}

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum InchiDiagnosticLevel {
    Warning,
    Error,
}
#[wasm_bindgen]
pub struct InchiResult {
    inner: ck::InchiResult,
}
#[wasm_bindgen]
impl InchiResult {
    #[wasm_bindgen(getter, js_name=inchi)]
    pub fn inchi(&self) -> String {
        self.inner.inchi().to_owned()
    }
    #[wasm_bindgen(getter, js_name=returnCode)]
    pub fn return_code(&self) -> i32 {
        self.inner.return_code()
    }
    #[wasm_bindgen(getter, js_name=message)]
    pub fn message(&self) -> String {
        self.inner.message().to_owned()
    }
    #[wasm_bindgen(getter, js_name=log)]
    pub fn log(&self) -> String {
        self.inner.log().to_owned()
    }
    #[wasm_bindgen(getter, js_name=auxInfo)]
    pub fn aux_info(&self) -> String {
        self.inner.aux_info().to_owned()
    }
    #[wasm_bindgen(getter, js_name=diagnostics)]
    pub fn diagnostics(&self) -> Vec<InchiDiagnostic> {
        self.inner
            .diagnostics()
            .iter()
            .cloned()
            .map(|inner| InchiDiagnostic { inner })
            .collect()
    }
}
#[wasm_bindgen]
pub struct InchiDiagnostic {
    inner: ck::InchiDiagnostic,
}
#[wasm_bindgen]
impl InchiDiagnostic {
    #[wasm_bindgen(getter, js_name=level)]
    pub fn level(&self) -> InchiDiagnosticLevel {
        match self.inner.level() {
            ck::InchiDiagnosticLevel::Warning => InchiDiagnosticLevel::Warning,
            ck::InchiDiagnosticLevel::Error => InchiDiagnosticLevel::Error,
        }
    }
    #[wasm_bindgen(getter, js_name=message)]
    pub fn message(&self) -> String {
        self.inner.message().to_owned()
    }
}

#[wasm_bindgen(js_name=inchiToKey)]
pub fn inchi_to_key(text: &str) -> Result<String, JsValue> {
    ck::inchi_to_key(text).map_err(error)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=toInchiResult)]
    pub fn to_inchi_result(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "InchiWriteParams | InchiWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<InchiResult, JsValue> {
        self.inner
            .to_inchi_result_with_params(&write_params(&params)?)
            .map(|inner| InchiResult { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=toInchiResultWithParams)]
    pub fn to_inchi_result_with_params(
        &self,
        params: &InchiWriteParams,
    ) -> Result<InchiResult, JsValue> {
        self.inner
            .to_inchi_result_with_params(&params.inner)
            .map(|inner| InchiResult { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fromInchi)]
    pub fn from_inchi(
        text: &str,
        #[wasm_bindgen(unchecked_optional_param_type = "InchiReadParams | InchiReadParamsOptions")]
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
        #[wasm_bindgen(
            unchecked_optional_param_type = "InchiWriteParams | InchiWriteParamsOptions"
        )]
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
        #[wasm_bindgen(
            unchecked_optional_param_type = "InchiWriteParams | InchiWriteParamsOptions"
        )]
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
