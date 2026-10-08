//! Public fingerprint source errors retain their exact causal chain.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AtomPairReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl AtomPairReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
#[wasm_bindgen]
pub struct FingerprintPreparationError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl FingerprintPreparationError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
#[wasm_bindgen]
pub struct FingerprintJsonError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl FingerprintJsonError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
#[cfg(feature = "batch")]
#[wasm_bindgen]
pub struct BatchFingerprintOutputError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[cfg(feature = "batch")]
#[wasm_bindgen]
impl BatchFingerprintOutputError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn atom_pair_error(source: &ck::AtomPairReadError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::AtomPairReadError::Preparation(_) => "Preparation",
        ck::AtomPairReadError::Generator(_) => "Generator",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("AtomPairReadError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        AtomPairReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
pub(crate) fn preparation_error(
    source: &ck::FingerprintPreparationError,
) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::FingerprintPreparationError::MissingPreparedValence => "MissingPreparedValence",
        ck::FingerprintPreparationError::RingPreparation(_) => "RingPreparation",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("FingerprintPreparationError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        FingerprintPreparationError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
pub(crate) fn json_error(source: &ck::FingerprintJsonError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::FingerprintJsonError::Parse(_) => "Parse",
        ck::FingerprintJsonError::Invalid(_) => "Invalid",
        ck::FingerprintJsonError::UnsupportedComponent { .. } => "UnsupportedComponent",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("FingerprintJsonError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        FingerprintJsonError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    if let ck::FingerprintJsonError::UnsupportedComponent {
        component,
        source_type,
    } = source
    {
        set(&error, "component", (*component).into())?;
        set(&error, "sourceType", source_type.as_str().into())?;
    }
    Ok(error)
}
#[cfg(feature = "batch")]
pub(crate) fn output_error(source: &ck::BatchFingerprintOutputError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::BatchFingerprintOutputError::MissingAdditionalOutput { .. } => {
            "MissingAdditionalOutput"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BatchFingerprintOutputError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        BatchFingerprintOutputError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    let ck::BatchFingerprintOutputError::MissingAdditionalOutput { fingerprint_kind } = source;
    set(&error, "fingerprintKind", (*fingerprint_kind).into())?;
    Ok(error)
}

#[wasm_bindgen]
pub struct MorganReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl MorganReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}

pub(crate) fn morgan_error(source: &ck::MorganReadError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::MorganReadError::Preparation(_) => "Preparation",
        ck::MorganReadError::Generator(_) => "Generator",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("MorganReadError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        MorganReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}

#[wasm_bindgen]
pub struct TopologicalTorsionReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl TopologicalTorsionReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "fingerprints".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        self.kind.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}

pub(crate) fn topological_torsion_error(
    source: &ck::TopologicalTorsionReadError,
) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::TopologicalTorsionReadError::Preparation(_) => "Preparation",
        ck::TopologicalTorsionReadError::Generator(_) => "Generator",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("TopologicalTorsionReadError");
    let error: JsValue = error.into();
    set(&error, "domain", "fingerprints".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        TopologicalTorsionReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
