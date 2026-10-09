//! Typed tuple and error projection of the canonical chiral-center query.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::host_values::bool_value;
use cosmolkit_wasm::rust as ck;
use std::error::Error as _;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct StereoReadError {
    kind: String,
    message: String,
    cause: JsValue,
}

#[wasm_bindgen]
impl StereoReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "stereo".into()
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

fn read_error(source: &ck::StereoReadError) -> Result<JsValue, JsValue> {
    let kind = match source {
        ck::StereoReadError::InvalidTopology(_) => "InvalidTopology",
        ck::StereoReadError::Valence(_) => "Valence",
        ck::StereoReadError::Rings(_) => "Rings",
        ck::StereoReadError::PotentialStereo(_) => "PotentialStereo",
        ck::StereoReadError::CipLabeler(_) => "CipLabeler",
        ck::StereoReadError::PropertyString(_) => "PropertyString",
        ck::StereoReadError::CipLabelEncoding { .. } => "CipLabelEncoding",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("StereoReadError");
    let error: JsValue = error.into();
    set(&error, "domain", "stereo".into())?;
    set(&error, "kind", kind.into())?;
    if !cause.is_null() {
        set(&error, "cause", cause.clone())?;
    }
    set(
        &error,
        "detail",
        StereoReadError {
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
    #[wasm_bindgen(js_name = findChiralCenters, unchecked_return_type = "Array<[number, string]>")]
    pub fn find_chiral_centers(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_unassigned: JsValue,
    ) -> Result<JsValue, JsValue> {
        let include_unassigned = if include_unassigned.is_undefined() {
            false
        } else {
            bool_value(&include_unassigned, "includeUnassigned")?
        };
        let rows = self
            .inner
            .find_chiral_centers(include_unassigned)
            .map_err(|error| read_error(&error).unwrap_or_else(|error| error))?;
        let output = js_sys::Array::new();
        for (atom, label) in rows {
            let row = js_sys::Array::new();
            row.push(&JsValue::from_f64(atom as f64));
            row.push(&label.into());
            output.push(&row);
        }
        Ok(output.into())
    }
}
