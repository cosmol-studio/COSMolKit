//! Stable template and layout error vocabulary from the public facade.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct Coordinate2DTemplateError {
    inner: ck::Coordinate2DTemplateError,
}
#[wasm_bindgen]
impl Coordinate2DTemplateError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "depict".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        template_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DTemplateError::InvalidDefaultRow { index, .. } => {
                JsValue::from_f64(*index as f64)
            }
            ck::Coordinate2DTemplateError::InvalidTopology { index, .. } => {
                JsValue::from_f64(*index as f64)
            }
            ck::Coordinate2DTemplateError::NonElementIdentity { index, .. } => {
                JsValue::from_f64(*index as f64)
            }
            ck::Coordinate2DTemplateError::RingInitialization { index, .. } => {
                JsValue::from_f64(*index as f64)
            }
            ck::Coordinate2DTemplateError::ConnectedComponents { index, .. } => {
                JsValue::from_f64(*index as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=path,unchecked_return_type="string | null")]
    pub fn path(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DTemplateError::ExternalOpen { path, .. } => JsValue::from_str(path),
            ck::Coordinate2DTemplateError::ExternalRead { path, .. } => JsValue::from_str(path),
            ck::Coordinate2DTemplateError::ExternalInvalidSmarts { path, .. } => {
                JsValue::from_str(path)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=line,unchecked_return_type="number | null")]
    pub fn line(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DTemplateError::ExternalInvalidSmarts { line, .. } => {
                JsValue::from_f64(*line as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=row,unchecked_return_type="string | null")]
    pub fn row(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DTemplateError::MissingCoordinates { row, .. } => JsValue::from_str(row),
            ck::Coordinate2DTemplateError::ThreeDimensionalCoordinates { row, .. } => {
                JsValue::from_str(row)
            }
            ck::Coordinate2DTemplateError::MultipleFragments { row, .. } => JsValue::from_str(row),
            ck::Coordinate2DTemplateError::NotRingSystem { row, .. } => JsValue::from_str(row),
            _ => JsValue::NULL,
        }
    }
}
fn template_error_kind(source: &ck::Coordinate2DTemplateError) -> &'static str {
    use ck::Coordinate2DTemplateError as E;
    match source {
        E::InvalidDefaultRow { .. } => "InvalidDefaultRow",
        E::InvalidTopology { .. } => "InvalidTopology",
        E::NonElementIdentity { .. } => "NonElementIdentity",
        E::RingInitialization { .. } => "RingInitialization",
        E::ConnectedComponents { .. } => "ConnectedComponents",
        E::ExternalOpen { .. } => "ExternalOpen",
        E::ExternalRead { .. } => "ExternalRead",
        E::ExternalInvalidSmarts { .. } => "ExternalInvalidSmarts",
        E::MissingCoordinates { .. } => "MissingCoordinates",
        E::ThreeDimensionalCoordinates { .. } => "ThreeDimensionalCoordinates",
        E::MultipleFragments { .. } => "MultipleFragments",
        E::NotRingSystem { .. } => "NotRingSystem",
    }
}
pub(crate) fn template_error(source: &ck::Coordinate2DTemplateError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: pub enum Coordinate2DTemplateError {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("Coordinate2DTemplateError");
    let error: JsValue = error.into();
    set(&error, "domain", "depict".into())?;
    set(&error, "kind", template_error_kind(source).into())?;
    let detail = Coordinate2DTemplateError {
        inner: source.clone(),
    };
    set(&error, "index", detail.index())?;
    set(&error, "path", detail.path())?;
    set(&error, "line", detail.line())?;
    set(&error, "row", detail.row())?;
    set(&error, "detail", detail.into())?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct Coordinate2DLayoutError {
    inner: ck::Coordinate2DLayoutError,
}
#[wasm_bindgen]
impl Coordinate2DLayoutError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "depict".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        layout_error_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::AtomAlreadyEmbedded { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::AtomNotEmbedded { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::EmptyAttachment { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::InvalidRankProperty { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { atom, .. } => {
                JsValue::from_f64(*atom as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::AtomIndexOutOfRange { atom_count, .. } => {
                JsValue::from_f64(*atom_count as f64)
            }
            ck::Coordinate2DLayoutError::AtomCountTooLarge { atom_count, .. } => {
                JsValue::from_f64(*atom_count as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=count,unchecked_return_type="number | null")]
    pub fn count(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::NotEnoughEmbeddedNeighbors { count, .. } => {
                JsValue::from_f64(*count as f64)
            }
            ck::Coordinate2DLayoutError::GeometryNotEnoughNeighbors { count, .. } => {
                JsValue::from_f64(*count as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=centre,unchecked_return_type="number | null")]
    pub fn centre(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::NonTetrahedralNoLigand { centre, .. } => {
                JsValue::from_f64(*centre as f64)
            }
            ck::Coordinate2DLayoutError::NonTetrahedralLigandOverflow { centre, .. } => {
                JsValue::from_f64(*centre as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn bond(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::CisTransBondInvalid { bond, .. } => {
                JsValue::from_f64(*bond as f64)
            }
            ck::Coordinate2DLayoutError::CollisionBondInvalid { bond, .. } => {
                JsValue::from_f64(*bond as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=first,unchecked_return_type="number | null")]
    pub fn first(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::UndefinedSamplingDistance { first, .. } => {
                JsValue::from_f64(*first as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=second,unchecked_return_type="number | null")]
    pub fn second(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::UndefinedSamplingDistance { second, .. } => {
                JsValue::from_f64(*second as f64)
            }
            _ => JsValue::NULL,
        }
    }
    #[wasm_bindgen(getter,js_name=key,unchecked_return_type="string | null")]
    pub fn key(&self) -> JsValue {
        match &self.inner {
            ck::Coordinate2DLayoutError::InvalidRankProperty { key, .. } => JsValue::from_str(key),
            _ => JsValue::NULL,
        }
    }
}
fn layout_error_kind(source: &ck::Coordinate2DLayoutError) -> &'static str {
    use ck::Coordinate2DLayoutError as E;
    match source {
        E::AtomIndexOutOfRange { .. } => "AtomIndexOutOfRange",
        E::AtomAlreadyEmbedded { .. } => "AtomAlreadyEmbedded",
        E::AtomNotEmbedded { .. } => "AtomNotEmbedded",
        E::NotEnoughEmbeddedNeighbors { .. } => "NotEnoughEmbeddedNeighbors",
        E::CoincidentPoints => "CoincidentPoints",
        E::InvalidAngle => "InvalidAngle",
        E::NoCommonAtoms => "NoCommonAtoms",
        E::MismatchedTopology => "MismatchedTopology",
        E::EmptyAttachment { .. } => "EmptyAttachment",
        E::NonTetrahedralNoLigand { .. } => "NonTetrahedralNoLigand",
        E::NonTetrahedralLigandOverflow { .. } => "NonTetrahedralLigandOverflow",
        E::CisTransBondInvalid { .. } => "CisTransBondInvalid",
        E::CollisionBondInvalid { .. } => "CollisionBondInvalid",
        E::UndefinedSamplingDistance { .. } => "UndefinedSamplingDistance",
        E::AtomCountTooLarge { .. } => "AtomCountTooLarge",
        E::InvalidRankProperty { .. } => "InvalidRankProperty",
        E::GeometryNotEnoughNeighbors { .. } => "GeometryNotEnoughNeighbors",
        E::TemplateMatch(..) => "TemplateMatch",
        E::GraphPath(..) => "GraphPath",
        E::GraphDistance(..) => "GraphDistance",
    }
}
pub(crate) fn layout_error(source: &ck::Coordinate2DLayoutError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: pub enum Coordinate2DLayoutError {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("Coordinate2DLayoutError");
    let error: JsValue = error.into();
    set(&error, "domain", "depict".into())?;
    set(&error, "kind", layout_error_kind(source).into())?;
    let detail = Coordinate2DLayoutError {
        inner: source.clone(),
    };
    set(&error, "atom", detail.atom())?;
    set(&error, "atomCount", detail.atom_count())?;
    set(&error, "count", detail.count())?;
    set(&error, "centre", detail.centre())?;
    set(&error, "bond", detail.bond())?;
    set(&error, "first", detail.first())?;
    set(&error, "second", detail.second())?;
    set(&error, "key", detail.key())?;
    set(&error, "detail", detail.into())?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
