//! SMILES text/error transport. No chemistry or source defaults live here.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct SmilesStereoError {
    inner: ck::SmilesStereoError,
}
#[wasm_bindgen]
impl SmilesStereoError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "smiles".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        stereo_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn id(&self) -> JsValue {
        match self.inner {
            ck::SmilesStereoError::SelectedConformerMissing { id } => (id as u32).into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn actual(&self) -> JsValue {
        match self.inner {
            ck::SmilesStereoError::ValenceRows { actual, .. }
            | ck::SmilesStereoError::RingRows { actual, .. } => (actual as u32).into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn expected(&self) -> JsValue {
        match self.inner {
            ck::SmilesStereoError::ValenceRows { expected, .. }
            | ck::SmilesStereoError::RingRows { expected, .. } => (expected as u32).into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn field(&self) -> JsValue {
        match self.inner {
            ck::SmilesStereoError::ValenceRows { field, .. } => field.into(),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "string | null")]
    pub fn dimension(&self) -> JsValue {
        match self.inner {
            ck::SmilesStereoError::RingRows { dimension, .. } => dimension.into(),
            _ => JsValue::NULL,
        }
    }
}
fn stereo_kind(source: &ck::SmilesStereoError) -> &'static str {
    match source {
        ck::SmilesStereoError::Properties(_) => "Properties",
        ck::SmilesStereoError::SelectedConformerMissing { .. } => "SelectedConformerMissing",
        ck::SmilesStereoError::ValenceRows { .. } => "ValenceRows",
        ck::SmilesStereoError::RingRows { .. } => "RingRows",
        ck::SmilesStereoError::Coordinates(_) => "Coordinates",
        ck::SmilesStereoError::Directions(_) => "Directions",
        ck::SmilesStereoError::Assignment(_) => "Assignment",
        ck::SmilesStereoError::Rings(_) => "Rings",
        ck::SmilesStereoError::Valence(_) => "Valence",
    }
}
pub(crate) fn smiles_stereo_error(source: &ck::SmilesStereoError) -> Result<JsValue, JsValue> {
    let detail = SmilesStereoError {
        inner: source.clone(),
    };
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SmilesStereoError");
    let error: JsValue = error.into();
    set(&error, "domain", "smiles".into())?;
    set(&error, "kind", detail.kind().into())?;
    let cause = detail.cause()?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    set(&error, "id", detail.id())?;
    set(&error, "actual", detail.actual())?;
    set(&error, "expected", detail.expected())?;
    set(&error, "field", detail.field())?;
    set(&error, "dimension", detail.dimension())?;
    set(&error, "detail", detail.into())?;
    Ok(error)
}

#[wasm_bindgen]
pub struct SmilesWriteError {
    kind: String,
    message: String,
    cause: JsValue,
}

#[wasm_bindgen]
impl SmilesWriteError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "smiles".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}

pub(crate) fn smiles_write_error(source: &ck::SmilesWriteError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::SmilesWriteError::Write(..) => "Write",
        ck::SmilesWriteError::Fragment(..) => "Fragment",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SmilesWriteError");
    let error: JsValue = error.into();
    set(&error, "domain", "smiles".into())?;
    set(&error, "kind", kind.into())?;
    if !cause.is_null() {
        set(&error, "cause", cause.clone())?;
    }
    set(
        &error,
        "detail",
        SmilesWriteError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = toSmiles)]
    pub fn to_smiles(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmilesWriteParams | SmilesWriteParamsOptions"
        )]
        params: JsValue,
    ) -> Result<String, JsValue> {
        let params = crate::smiles_parameters::SmilesWriteParams::from_js_configuration(&params)?;
        self.to_smiles_with_params(&params)
    }
    #[wasm_bindgen(js_name = toSmilesWithParams)]
    pub fn to_smiles_with_params(
        &self,
        params: &crate::smiles_parameters::SmilesWriteParams,
    ) -> Result<String, JsValue> {
        self.inner
            .to_smiles_with_params(&params.inner)
            .map_err(|error| smiles_write_error(&error).unwrap_or_else(|error| error))
            .and_then(|value| crate::host_values::text(&value))
    }
}

#[wasm_bindgen]
pub struct SmilesError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl SmilesError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "smiles".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn smiles_error(source: &ck::SmilesError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::SmilesError::Parse(..) => "Parse",
        ck::SmilesError::Hydrogen(..) => "Hydrogen",
        ck::SmilesError::Sanitize(..) => "Sanitize",
        ck::SmilesError::Stereo(..) => "Stereo",
        ck::SmilesError::Construction(..) => "Construction",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("SmilesError");
    let error: JsValue = error.into();
    set(&error, "domain", "smiles".into())?;
    set(&error, "kind", kind.into())?;
    if !cause.is_null() {
        set(&error, "cause", cause.clone())?;
    }
    set(
        &error,
        "detail",
        SmilesError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(error)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name = fromSmilesWithSanitize)]
    pub fn from_smiles_with_sanitize(
        smiles: &str,
        #[wasm_bindgen(unchecked_param_type = "boolean")] sanitize: JsValue,
    ) -> Result<Molecule, JsValue> {
        let sanitize = crate::host_values::bool_value(&sanitize, "sanitize")?;
        cosmolkit_wasm::Molecule::from_smiles_with_sanitize(smiles, sanitize)
            .map(|inner| Molecule {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| smiles_error(&e).unwrap_or_else(|e| e))
    }

    #[wasm_bindgen(js_name = fromSmiles)]
    pub fn from_smiles(
        smiles: &str,
        #[wasm_bindgen(
            unchecked_optional_param_type = "SmilesParseParams | SmilesParseParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Molecule, JsValue> {
        let params = crate::smiles_parameters::SmilesParseParams::from_js_configuration(&params)?;
        Self::from_smiles_with_params(smiles, &params)
    }
    #[wasm_bindgen(js_name = fromSmilesWithParams)]
    pub fn from_smiles_with_params(
        smiles: &str,
        params: &crate::smiles_parameters::SmilesParseParams,
    ) -> Result<Molecule, JsValue> {
        cosmolkit_wasm::Molecule::from_smiles_with_params(smiles, &params.inner)
            .map(|inner| Molecule {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| smiles_error(&e).unwrap_or_else(|e| e))
    }
}
