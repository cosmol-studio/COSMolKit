//! Native drawing exceptions retain the facade enum and source chain.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct DrawingError {
    kind: String,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl DrawingError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "drawing".into()
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
pub(crate) fn drawing_error(source: &ck::DrawingError) -> Result<JsValue, JsValue> {
    use ck::DrawingError as E;
    let kind = match source {
        E::SvgTextProjection(..) => "SvgTextProjection",
        E::PropertyString(..) => "PropertyString",
        E::VariationArray(..) => "VariationArray",
        E::VariationIndex { .. } => "VariationIndex",
        E::DataFieldDouble { .. } => "DataFieldDouble",
        E::Property(..) => "Property",
        E::Topology(..) => "Topology",
        E::Coordinates(..) => "Coordinates",
        E::Mapping(..) => "Mapping",
        E::Kekulize(..) => "Kekulize",
        E::Hydrogen(..) => "Hydrogen",
        E::Wedge(..) => "Wedge",
        E::Valence(..) => "Valence",
        E::CoordinateGeneration(..) => "CoordinateGeneration",
        E::SvgParse(..) => "SvgParse",
        E::PngEncode(..) => "PngEncode",
        E::StateRows { .. } => "StateRows",
        E::HydrogenAppend { .. } => "HydrogenAppend",
        E::InvalidDimensions { .. } => "InvalidDimensions",
        E::PixmapAllocation { .. } => "PixmapAllocation",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("DrawingError");
    let error: JsValue = error.into();
    set(&error, "domain", "drawing".into())?;
    set(&error, "kind", kind.into())?;
    set(
        &error,
        "detail",
        DrawingError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
        }
        .into(),
    )?;
    match source {
        E::InvalidDimensions { width, height } | E::PixmapAllocation { width, height } => {
            set(&error, "width", JsValue::from_f64(*width as f64))?;
            set(&error, "height", JsValue::from_f64(*height as f64))?;
        }
        E::StateRows {
            field,
            actual,
            expected,
        } => {
            set(&error, "field", JsValue::from_str(field))?;
            set(&error, "actual", JsValue::from_f64(*actual as f64))?;
            set(&error, "expected", JsValue::from_f64(*expected as f64))?;
        }
        E::HydrogenAppend { row, reason } => {
            set(
                &error,
                "row",
                row.map_or(JsValue::NULL, |v| JsValue::from_f64(v as f64)),
            )?;
            set(&error, "reason", JsValue::from_str(reason))?;
        }
        _ => {}
    }
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
