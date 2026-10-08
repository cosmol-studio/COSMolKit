//! Exact MatrixError variants, fields and recursive source causes.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MatrixError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl MatrixError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "matrices".into()
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
    #[wasm_bindgen(getter,js_name=position,unchecked_return_type="number | null")]
    pub fn field_0(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"position".into())
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn field_1(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atom".into())
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn field_2(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atomCount".into())
    }
    #[wasm_bindgen(getter,js_name=firstPosition,unchecked_return_type="number | null")]
    pub fn field_3(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"firstPosition".into())
    }
    #[wasm_bindgen(getter,js_name=secondPosition,unchecked_return_type="number | null")]
    pub fn field_4(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"secondPosition".into())
    }
    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn field_5(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"bond".into())
    }
    #[wasm_bindgen(getter,js_name=bondCount,unchecked_return_type="number | null")]
    pub fn field_6(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"bondCount".into())
    }
    #[wasm_bindgen(getter,js_name=endpoint,unchecked_return_type="string | null")]
    pub fn field_7(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"endpoint".into())
    }
    #[wasm_bindgen(getter,js_name=order,unchecked_return_type="BondOrder | null")]
    pub fn field_8(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"order".into())
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn field_9(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"conformerId".into())
    }
    #[wasm_bindgen(getter,js_name=dimension,unchecked_return_type="number | null")]
    pub fn field_10(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"dimension".into())
    }
}
pub(crate) fn matrix_error(source: &ck::MatrixError) -> Result<JsValue, JsValue> {
    use ck::MatrixError as E;
    let fields = js_sys::Object::new();
    set(&fields, "position", JsValue::NULL)?;
    set(&fields, "atom", JsValue::NULL)?;
    set(&fields, "atomCount", JsValue::NULL)?;
    set(&fields, "firstPosition", JsValue::NULL)?;
    set(&fields, "secondPosition", JsValue::NULL)?;
    set(&fields, "bond", JsValue::NULL)?;
    set(&fields, "bondCount", JsValue::NULL)?;
    set(&fields, "endpoint", JsValue::NULL)?;
    set(&fields, "order", JsValue::NULL)?;
    set(&fields, "conformerId", JsValue::NULL)?;
    set(&fields, "dimension", JsValue::NULL)?;
    let kind = match source {
        E::InvalidTopology(..) => "InvalidTopology",
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::ActiveBondsWithoutAtoms => "ActiveBondsWithoutAtoms",
        E::ActiveAtomOutOfRange {
            position,
            atom,
            atom_count,
        } => {
            set(&fields, "position", JsValue::from(*position as u32))?;
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(&fields, "atomCount", JsValue::from(*atom_count as u32))?;
            "ActiveAtomOutOfRange"
        }
        E::DuplicateActiveAtom {
            atom,
            first_position,
            second_position,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(
                &fields,
                "firstPosition",
                JsValue::from(*first_position as u32),
            )?;
            set(
                &fields,
                "secondPosition",
                JsValue::from(*second_position as u32),
            )?;
            "DuplicateActiveAtom"
        }
        E::ActiveBondOutOfRange {
            position,
            bond,
            bond_count,
        } => {
            set(&fields, "position", JsValue::from(*position as u32))?;
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(&fields, "bondCount", JsValue::from(*bond_count as u32))?;
            "ActiveBondOutOfRange"
        }
        E::DuplicateActiveBond {
            bond,
            first_position,
            second_position,
        } => {
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(
                &fields,
                "firstPosition",
                JsValue::from(*first_position as u32),
            )?;
            set(
                &fields,
                "secondPosition",
                JsValue::from(*second_position as u32),
            )?;
            "DuplicateActiveBond"
        }
        E::ActiveBondEndpointMissing {
            bond,
            endpoint,
            atom,
        } => {
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(&fields, "endpoint", JsValue::from(*endpoint))?;
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            "ActiveBondEndpointMissing"
        }
        E::UnsupportedBondOrder { bond, order } => {
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(&fields, "order", JsValue::from(*order as u32))?;
            "UnsupportedBondOrder"
        }
        E::No3dConformer => "No3dConformer",
        E::ConformerNotFound { conformer_id } => {
            set(&fields, "conformerId", JsValue::from(*conformer_id as u32))?;
            "ConformerNotFound"
        }
        E::MatrixDimensionOverflow { dimension } => {
            set(&fields, "dimension", JsValue::from(*dimension as u32))?;
            "MatrixDimensionOverflow"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("MatrixError");
    set(&error, "domain", "matrices".into())?;
    set(&error, "kind", kind.into())?;
    set(&error, "cause", cause.clone())?;
    for key in js_sys::Object::keys(&fields).iter() {
        let key = key.as_string().expect("own field key");
        set(
            &error,
            &key,
            js_sys::Reflect::get(&fields, &key.clone().into())?,
        )?;
    }
    set(
        &error,
        "detail",
        MatrixError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(error.into())
}
