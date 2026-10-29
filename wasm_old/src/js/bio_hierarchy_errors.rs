//! Structured BIO hierarchy errors, shared with validation and selection projections.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct BioStructureError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioStructureError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="number | null")]
    pub fn value(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"value".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=start,unchecked_return_type="number | null")]
    pub fn start(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"start".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=len,unchecked_return_type="number | null")]
    pub fn len(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"len".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=tableLen,unchecked_return_type="number | null")]
    pub fn table_len(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"tableLen".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=table,unchecked_return_type="string | null")]
    pub fn table(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"table".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=expectedStart,unchecked_return_type="number | null")]
    pub fn expected_start(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"expectedStart".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=actualStart,unchecked_return_type="number | null")]
    pub fn actual_start(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"actualStart".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=covered,unchecked_return_type="number | null")]
    pub fn covered(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"covered".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=index,unchecked_return_type="number | null")]
    pub fn index(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"index".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"atomCount".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=coordinateCount,unchecked_return_type="number | null")]
    pub fn coordinate_count(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"coordinateCount".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=entityId,unchecked_return_type="number | null")]
    pub fn entity_id(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"entityId".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=subchain,unchecked_return_type="string | null")]
    pub fn subchain(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"subchain".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=operation,unchecked_return_type="string | null")]
    pub fn operation(&self) -> JsValue {
        js_sys::Reflect::get(&self.context, &"operation".into())
            .map(|value| {
                if value.is_undefined() {
                    JsValue::NULL
                } else {
                    value
                }
            })
            .unwrap_or(JsValue::NULL)
    }
}
fn make_biostructureerror(
    source: &ck::BioStructureError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BioStructureError");
    let error: JsValue = error.into();
    set(&error, "domain", "bio".into())?;
    set(&error, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        let value = js_sys::Reflect::get(&context, &key)?;
        js_sys::Reflect::set(&error, &key, &value)?;
    }
    set(
        &error,
        "detail",
        BioStructureError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct BioRowModelError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioRowModelError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
fn make_biorowmodelerror(
    source: &ck::BioRowModelError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BioRowModelError");
    let error: JsValue = error.into();
    set(&error, "domain", "bio".into())?;
    set(&error, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        let value = js_sys::Reflect::get(&context, &key)?;
        js_sys::Reflect::set(&error, &key, &value)?;
    }
    set(
        &error,
        "detail",
        BioRowModelError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct BioRowChainError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioRowChainError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
fn make_biorowchainerror(
    source: &ck::BioRowChainError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BioRowChainError");
    let error: JsValue = error.into();
    set(&error, "domain", "bio".into())?;
    set(&error, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        let value = js_sys::Reflect::get(&context, &key)?;
        js_sys::Reflect::set(&error, &key, &value)?;
    }
    set(
        &error,
        "detail",
        BioRowChainError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
#[wasm_bindgen]
pub struct BioRowTraverseError {
    kind: String,
    message: String,
    cause: JsValue,
    context: JsValue,
}
#[wasm_bindgen]
impl BioRowTraverseError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "bio".into()
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
fn make_biorowtraverseerror(
    source: &ck::BioRowTraverseError,
    kind: &str,
    context: JsValue,
) -> Result<JsValue, JsValue> {
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BioRowTraverseError");
    let error: JsValue = error.into();
    set(&error, "domain", "bio".into())?;
    set(&error, "kind", kind.into())?;
    for key in js_sys::Object::keys(&js_sys::Object::from(context.clone())).iter() {
        let value = js_sys::Reflect::get(&context, &key)?;
        js_sys::Reflect::set(&error, &key, &value)?;
    }
    set(
        &error,
        "detail",
        BioRowTraverseError {
            kind: kind.into(),
            message: source.to_string(),
            cause: cause.clone(),
            context,
        }
        .into(),
    )?;
    if !cause.is_null() {
        set(&error, "cause", cause)?;
    }
    Ok(error)
}
pub(crate) fn structure_error(source: &ck::BioStructureError) -> Result<JsValue, JsValue> {
    let context: JsValue = js_sys::Object::new().into();
    let kind = match source {
        ck::BioStructureError::RowIndexTooLarge { value } => {
            set(&context, "value", ((*value) as f64).into())?;
            "RowIndexTooLarge"
        }
        ck::BioStructureError::RowSpanOverflow { start, len } => {
            set(&context, "start", (*start).into())?;
            set(&context, "len", (*len).into())?;
            "RowSpanOverflow"
        }
        ck::BioStructureError::RowSpanOutOfBounds {
            start,
            len,
            table_len,
        } => {
            set(&context, "start", (*start).into())?;
            set(&context, "len", (*len).into())?;
            set(&context, "tableLen", ((*table_len) as f64).into())?;
            "RowSpanOutOfBounds"
        }
        ck::BioStructureError::TableTooLarge { table, len } => {
            set(&context, "table", (*table).into())?;
            set(&context, "len", ((*len) as f64).into())?;
            "TableTooLarge"
        }
        ck::BioStructureError::NonContiguousSpan {
            table,
            expected_start,
            actual_start,
        } => {
            set(&context, "table", (*table).into())?;
            set(&context, "expectedStart", (*expected_start).into())?;
            set(&context, "actualStart", (*actual_start).into())?;
            "NonContiguousSpan"
        }
        ck::BioStructureError::IncompleteCoverage {
            table,
            covered,
            table_len,
        } => {
            set(&context, "table", (*table).into())?;
            set(&context, "covered", (*covered).into())?;
            set(&context, "tableLen", ((*table_len) as f64).into())?;
            "IncompleteCoverage"
        }
        ck::BioStructureError::ParentMismatch { table, index } => {
            set(&context, "table", (*table).into())?;
            set(&context, "index", (*index).into())?;
            "ParentMismatch"
        }
        ck::BioStructureError::RowReferenceOutOfBounds {
            table,
            index,
            table_len,
        } => {
            set(&context, "table", (*table).into())?;
            set(&context, "index", (*index).into())?;
            set(&context, "tableLen", ((*table_len) as f64).into())?;
            "RowReferenceOutOfBounds"
        }
        ck::BioStructureError::CoordinateCountMismatch {
            atom_count,
            coordinate_count,
        } => {
            set(&context, "atomCount", ((*atom_count) as f64).into())?;
            set(
                &context,
                "coordinateCount",
                ((*coordinate_count) as f64).into(),
            )?;
            "CoordinateCountMismatch"
        }
        ck::BioStructureError::EntitySubchainMismatch {
            entity_id,
            subchain,
        } => {
            set(&context, "entityId", (entity_id.value()).into())?;
            set(&context, "subchain", (subchain.as_str()).into())?;
            "EntitySubchainMismatch"
        }
        ck::BioStructureError::EmptyResidueSpan { operation } => {
            set(&context, "operation", (*operation).into())?;
            "EmptyResidueSpan"
        }
        ck::BioStructureError::ImpossibleCrystalAngle => "ImpossibleCrystalAngle",
        ck::BioStructureError::AtomNotFound => "AtomNotFound",
    };
    make_biostructureerror(source, kind, context)
}
pub(crate) fn row_model_error(source: &ck::BioRowModelError) -> Result<JsValue, JsValue> {
    make_biorowmodelerror(
        source,
        match source {
            ck::BioRowModelError::MissingModelNumber => "MissingModelNumber",
        },
        js_sys::Object::new().into(),
    )
}
pub(crate) fn row_chain_error(source: &ck::BioRowChainError) -> Result<JsValue, JsValue> {
    make_biorowchainerror(
        source,
        match source {
            ck::BioRowChainError::MissingCanonicalChainName => "MissingCanonicalChainName",
        },
        js_sys::Object::new().into(),
    )
}
pub(crate) fn row_traverse_error(source: &ck::BioRowTraverseError) -> Result<JsValue, JsValue> {
    make_biorowtraverseerror(
        source,
        match source {
            ck::BioRowTraverseError::Model(_) => "Model",
            ck::BioRowTraverseError::Chain(_) => "Chain",
        },
        js_sys::Object::new().into(),
    )
}
