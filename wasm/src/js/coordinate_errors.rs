//! Structured canonical coordinate errors retained across the WASM boundary.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct CoordinateInputError {
    inner: ck::CoordinateInputError,
}
#[wasm_bindgen]
impl CoordinateInputError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "coordinates".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        input_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn input_kind(e: &ck::CoordinateInputError) -> &'static str {
    use ck::CoordinateInputError as E;
    match e {
        E::UnknownZPolicy { .. } => "UnknownZPolicy",
        E::RowCount { .. } => "RowCount",
        E::Shape { .. } => "Shape",
        E::NonFinite { .. } => "NonFinite",
        E::NonZeroZ { .. } => "NonZeroZ",
        E::ThreeColumnsForbidden => "ThreeColumnsForbidden",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::ConformerIdOverflow { .. } => "ConformerIdOverflow",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
    }
}
pub(crate) fn input_error(source: &ck::CoordinateInputError) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("CoordinateInputError");
    let e: JsValue = e.into();
    set(&e, "domain", "coordinates".into())?;
    set(&e, "kind", input_kind(source).into())?;
    use ck::CoordinateInputError as E;
    match source {
        E::UnknownZPolicy { value } => set(&e, "value", value.as_str().into())?,
        E::RowCount {
            dimension,
            expected,
            actual,
        } => {
            set(&e, "dimension", (*dimension).into())?;
            set(&e, "expected", (*expected as f64).into())?;
            set(&e, "actual", (*actual as f64).into())?;
        }
        E::Shape {
            dimension,
            row,
            columns,
            expected,
        } => {
            set(&e, "dimension", (*dimension).into())?;
            set(&e, "row", (*row as f64).into())?;
            set(&e, "columns", (*columns as f64).into())?;
            set(&e, "expected", (*expected).into())?;
        }
        E::NonFinite {
            dimension,
            row,
            column,
            value,
        } => {
            set(&e, "dimension", (*dimension).into())?;
            set(&e, "row", (*row as f64).into())?;
            set(&e, "column", (*column as f64).into())?;
            set(&e, "value", (*value).into())?;
        }
        E::NonZeroZ { row, z } => {
            set(&e, "row", (*row as f64).into())?;
            set(&e, "z", (*z).into())?;
        }
        E::ConformerNotFound {
            conformer_id,
            count,
        } => {
            set(&e, "conformerId", (*conformer_id as f64).into())?;
            set(&e, "count", (*count as f64).into())?;
        }
        E::ConformerIdOverflow { max_id } => set(&e, "maxId", (*max_id as f64).into())?,
        E::ThreeColumnsForbidden | E::InvalidCoordinates(..) => {}
    }
    if let Some(cause) = source.source() {
        set(&e, "cause", source_error(cause)?)?;
    }
    set(
        &e,
        "detail",
        CoordinateInputError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
pub struct Coordinate3DReadError {
    inner: ck::Coordinate3DReadError,
}
#[wasm_bindgen]
impl Coordinate3DReadError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "coordinates".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        "ConformerNotFound".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> usize {
        match self.inner {
            ck::Coordinate3DReadError::ConformerNotFound { conformer_id, .. } => conformer_id,
        }
    }
    #[wasm_bindgen(getter)]
    pub fn count(&self) -> usize {
        match self.inner {
            ck::Coordinate3DReadError::ConformerNotFound { count, .. } => count,
        }
    }
}
pub(crate) fn read_error(source: &ck::Coordinate3DReadError) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("Coordinate3DReadError");
    let e: JsValue = e.into();
    let detail = Coordinate3DReadError {
        inner: source.clone(),
    };
    set(&e, "domain", "coordinates".into())?;
    set(&e, "kind", "ConformerNotFound".into())?;
    set(&e, "conformerId", (detail.conformer_id() as f64).into())?;
    set(&e, "count", (detail.count() as f64).into())?;
    set(&e, "detail", detail.into())?;
    Ok(e)
}
#[wasm_bindgen]
pub struct TransformError {
    inner: ck::TransformError,
}
#[wasm_bindgen]
impl TransformError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "transform".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        transform_kind(&self.inner).into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
}
fn transform_kind(e: &ck::TransformError) -> &'static str {
    use ck::TransformError as E;
    match e {
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::WeightCount { .. } => "WeightCount",
        E::NonFiniteWeight { .. } => "NonFiniteWeight",
        E::InvalidWeightSum { .. } => "InvalidWeightSum",
        E::NoSelectedAtoms => "NoSelectedAtoms",
        E::NonFinitePoint { .. } => "NonFinitePoint",
        E::NonFiniteTransform { .. } => "NonFiniteTransform",
        E::EigenDecompositionDidNotConverge => "EigenDecompositionDidNotConverge",
        E::No3dConformer => "No3dConformer",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::AtomsNotBonded { .. } => "AtomsNotBonded",
        E::RingBondNotMovable { .. } => "RingBondNotMovable",
        E::BothAngleBondsInRing { .. } => "BothAngleBondsInRing",
        E::CoincidentCoordinates { .. } => "CoincidentCoordinates",
        E::UndefinedRotationAxis => "UndefinedRotationAxis",
        E::InvalidTargetValue { .. } => "InvalidTargetValue",
    }
}
pub(crate) fn transform_error(source: &ck::TransformError) -> Result<JsValue, JsValue> {
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("TransformError");
    let e: JsValue = e.into();
    set(&e, "domain", "transform".into())?;
    set(&e, "kind", transform_kind(source).into())?;
    use ck::TransformError as E;
    match source {
        E::AtomOutOfRange {
            role,
            atom,
            atom_count,
        } => {
            set(&e, "role", (*role).into())?;
            set(&e, "atomIndex", (atom.index() as f64).into())?;
            set(&e, "atomCount", (*atom_count as f64).into())?;
        }
        E::ConformerNotFound { conformer_id } => {
            set(&e, "conformerId", (*conformer_id as f64).into())?
        }
        E::NonFinitePoint { role, axis } => {
            set(&e, "role", (*role).into())?;
            set(&e, "axis", (*axis).into())?;
        }
        E::WeightCount { actual, required } => {
            set(&e, "actual", (*actual as f64).into())?;
            set(&e, "required", (*required as f64).into())?;
        }
        E::NonFiniteWeight { atom, value } => {
            set(&e, "atomIndex", (atom.index() as f64).into())?;
            set(&e, "value", (*value).into())?;
        }
        E::InvalidWeightSum { sum } => set(&e, "sum", (*sum).into())?,
        E::NonFiniteTransform { row, column } => {
            set(&e, "row", (*row as f64).into())?;
            set(&e, "column", (*column as f64).into())?;
        }
        E::AtomsNotBonded {
            first_role,
            first,
            second_role,
            second,
        }
        | E::CoincidentCoordinates {
            first_role,
            first,
            second_role,
            second,
        } => {
            set(&e, "firstRole", (*first_role).into())?;
            set(&e, "first", (first.index() as f64).into())?;
            set(&e, "secondRole", (*second_role).into())?;
            set(&e, "second", (second.index() as f64).into())?;
        }
        E::RingBondNotMovable { bond } => set(&e, "bondIndex", (bond.index() as f64).into())?,
        E::BothAngleBondsInRing {
            first_bond,
            second_bond,
        } => {
            set(&e, "firstBond", (first_bond.index() as f64).into())?;
            set(&e, "secondBond", (second_bond.index() as f64).into())?;
        }
        E::InvalidTargetValue { quantity, value } => {
            set(&e, "quantity", (*quantity).into())?;
            set(&e, "value", (*value).into())?;
        }
        E::InvalidTopology(..)
        | E::InvalidCoordinates(..)
        | E::NoSelectedAtoms
        | E::EigenDecompositionDidNotConverge
        | E::No3dConformer
        | E::UndefinedRotationAxis => {}
    }
    if let Some(cause) = source.source() {
        set(&e, "cause", source_error(cause)?)?;
    }
    set(
        &e,
        "detail",
        TransformError {
            inner: source.clone(),
        }
        .into(),
    )?;
    Ok(e)
}
