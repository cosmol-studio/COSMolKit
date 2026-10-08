//! SMILES text/error transport. No chemistry or source defaults live here.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

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
    pub fn to_smiles(&self) -> Result<String, JsValue> {
        self.inner
            .to_smiles()
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
    pub fn from_smiles(smiles: &str) -> Result<Molecule, JsValue> {
        cosmolkit_wasm::Molecule::from_smiles(smiles)
            .map(|inner| Molecule {
                inner: std::sync::Arc::new(inner),
            })
            .map_err(|e| smiles_error(&e).unwrap_or_else(|e| e))
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
