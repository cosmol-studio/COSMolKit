//! Canonical descriptor error categories, exact context and real source chains.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct DescriptorReadError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl DescriptorReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "descriptors".into()
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
pub(crate) fn descriptor_read_error(source: &ck::DescriptorReadError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: pub enum DescriptorReadError {
    let kind = match source {
        ck::DescriptorReadError::MissingPreparedValence => "MissingPreparedValence",
        ck::DescriptorReadError::CachePoisoned => "CachePoisoned",
        ck::DescriptorReadError::MissingInitializedRings => "MissingInitializedRings",
        ck::DescriptorReadError::Algorithm { .. } => "Algorithm",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("DescriptorReadError");
    let e: JsValue = e.into();
    set(&e, "domain", "descriptors".into())?;
    set(&e, "kind", kind.into())?;
    set(
        &e,
        "detail",
        DescriptorReadError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&e, "cause", cause)?;
    }
    Ok(e)
}
#[wasm_bindgen]
pub struct DescriptorError {
    inner: ck::DescriptorError,
}
#[wasm_bindgen]
impl DescriptorError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "descriptors".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        descriptor_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter,js_name=function,unchecked_return_type="string | null")]
    pub fn function(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::TopologyEdit { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::InvalidConnectivityPath { function, .. } => {
                JsValue::from_str(function)
            }
            ck::DescriptorError::InvalidTopology { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Path { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::MissingComputedScalar { function, .. } => {
                JsValue::from_str(function)
            }
            ck::DescriptorError::MissingLabuteHydrogens { function, .. } => {
                JsValue::from_str(function)
            }
            ck::DescriptorError::MissingLabuteAsa { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::MissingCrippenMrContributions { function, .. } => {
                JsValue::from_str(function)
            }
            ck::DescriptorError::MissingCrippenMr { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Hydrogens { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::CountOverflow { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Valence { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::InvalidValenceRows { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Unsupported { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Ring { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Search { function, .. } => JsValue::from_str(function),
            ck::DescriptorError::Stereo { function, .. } => JsValue::from_str(function),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=field,unchecked_return_type="string | null")]
    pub fn field(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MissingFinalHydrogenState { field, .. } => {
                JsValue::from_str(field)
            }
            ck::DescriptorError::CrippenParamNumeric { field, .. } => JsValue::from_str(field),
            ck::DescriptorError::InvalidCrippenOptionalRows { field, .. } => {
                JsValue::from_str(field)
            }
            ck::DescriptorError::CountOverflow { field, .. } => JsValue::from_str(field),
            ck::DescriptorError::InvalidValenceRows { field, .. } => JsValue::from_str(field),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=expectedRows,unchecked_return_type="number | null")]
    pub fn expected_rows(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::InvalidConnectivityPath { expected_rows, .. } => {
                JsValue::from_f64(*expected_rows as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=actualRows,unchecked_return_type="number | null")]
    pub fn actual_rows(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::InvalidConnectivityPath { actual_rows, .. } => {
                actual_rows.map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64))
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::InvalidHallKierContributionRows { actual, .. } => {
                JsValue::from_f64(*actual as f64)
            }
            ck::DescriptorError::InvalidCrippenOptionalRows { actual, .. } => {
                JsValue::from_f64(*actual as f64)
            }
            ck::DescriptorError::InvalidValenceRows { actual, .. } => {
                JsValue::from_f64(*actual as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=minimum,unchecked_return_type="number | null")]
    pub fn minimum(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::InvalidHallKierContributionRows { minimum, .. } => {
                JsValue::from_f64(*minimum as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=includeSulfurPhosphorus,unchecked_return_type="boolean | null")]
    pub fn include_sulfur_phosphorus(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MissingComputedScalar {
                include_sulfur_phosphorus,
                ..
            } => JsValue::from_bool(*include_sulfur_phosphorus),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=contribsLen,unchecked_return_type="number | null")]
    pub fn contribs_len(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MismatchedBinArrays { contribs_len, .. } => {
                JsValue::from_f64(*contribs_len as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=binPropLen,unchecked_return_type="number | null")]
    pub fn bin_prop_len(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MismatchedBinArrays { bin_prop_len, .. } => {
                JsValue::from_f64(*bin_prop_len as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=binsLen,unchecked_return_type="number | null")]
    pub fn bins_len(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MismatchedBinArrays { bins_len, .. } => {
                JsValue::from_f64(*bins_len as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=cell,unchecked_return_type="string | null")]
    pub fn cell(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::CrippenParamNumeric { cell, .. } => JsValue::from_str(cell),
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=expected,unchecked_return_type="number | null")]
    pub fn expected(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::InvalidCrippenOptionalRows { expected, .. } => {
                JsValue::from_f64(*expected as f64)
            }
            ck::DescriptorError::InvalidValenceRows { expected, .. } => {
                JsValue::from_f64(*expected as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=row,unchecked_return_type="number | null")]
    pub fn row(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::MissingCrippenDefaultPattern { row, .. } => {
                JsValue::from_f64(*row as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=detail,unchecked_return_type="string | null")]
    pub fn detail(&self) -> JsValue {
        match &self.inner {
            ck::DescriptorError::Unsupported { detail, .. } => JsValue::from_str(detail),
            _ => JsValue::NULL,
        }
    }
}
fn descriptor_error_kind(source: &ck::DescriptorError) -> &'static str {
    use ck::DescriptorError as E;
    match source {
        E::TopologyEdit { .. } => "TopologyEdit",
        E::MissingFinalHydrogenState { .. } => "MissingFinalHydrogenState",
        E::InvalidConnectivityPath { .. } => "InvalidConnectivityPath",
        E::InvalidTopology { .. } => "InvalidTopology",
        E::Path { .. } => "Path",
        E::InvalidHallKierContributionRows { .. } => "InvalidHallKierContributionRows",
        E::MissingComputedScalar { .. } => "MissingComputedScalar",
        E::MissingLabuteHydrogens { .. } => "MissingLabuteHydrogens",
        E::MissingLabuteAsa { .. } => "MissingLabuteAsa",
        E::MismatchedBinArrays { .. } => "MismatchedBinArrays",
        E::CrippenParamNumeric { .. } => "CrippenParamNumeric",
        E::MissingCrippenMrContributions { .. } => "MissingCrippenMrContributions",
        E::MissingCrippenMr { .. } => "MissingCrippenMr",
        E::InvalidCrippenOptionalRows { .. } => "InvalidCrippenOptionalRows",
        E::MissingCrippenDefaultPattern { .. } => "MissingCrippenDefaultPattern",
        E::Hydrogens { .. } => "Hydrogens",
        E::CountOverflow { .. } => "CountOverflow",
        E::Valence { .. } => "Valence",
        E::InvalidValenceRows { .. } => "InvalidValenceRows",
        E::Unsupported { .. } => "Unsupported",
        E::Ring { .. } => "Ring",
        E::Search { .. } => "Search",
        E::Stereo { .. } => "Stereo",
    }
}
pub(crate) fn descriptor_error(source: &ck::DescriptorError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: pub enum DescriptorError {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("DescriptorError");
    let error: JsValue = error.into();
    set(&error, "domain", "descriptors".into())?;
    set(&error, "kind", descriptor_error_kind(source).into())?;
    let detail = DescriptorError {
        inner: source.clone(),
    };
    set(&error, "function", detail.function())?;
    set(&error, "field", detail.field())?;
    set(&error, "expectedRows", detail.expected_rows())?;
    set(&error, "actualRows", detail.actual_rows())?;
    set(&error, "actual", detail.actual())?;
    set(&error, "minimum", detail.minimum())?;
    set(
        &error,
        "includeSulfurPhosphorus",
        detail.include_sulfur_phosphorus(),
    )?;
    set(&error, "contribsLen", detail.contribs_len())?;
    set(&error, "binPropLen", detail.bin_prop_len())?;
    set(&error, "binsLen", detail.bins_len())?;
    set(&error, "cell", detail.cell())?;
    set(&error, "expected", detail.expected())?;
    set(&error, "row", detail.row())?;
    set(&error, "detail", detail.into())?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
