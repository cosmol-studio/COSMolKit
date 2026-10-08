//! Hash values retain all 64 bits and canonical typed failure fields.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::host_values::{sequence, u32_value};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MoleculeHashError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl MoleculeHashError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "molecular_hash".into()
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
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"actual".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"atomCount".into()).unwrap_or(JsValue::NULL)
    }
}
#[wasm_bindgen]
pub struct CipRankError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl CipRankError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "cip_ranking".into()
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
    #[wasm_bindgen(getter,js_name=field,unchecked_return_type="string | null")]
    pub fn field(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"field".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=actual,unchecked_return_type="number | null")]
    pub fn actual(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"actual".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn atom_count(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"atomCount".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn atom(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"atom".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="number | bigint | null")]
    pub fn value(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"value".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=mapNumber,unchecked_return_type="number | null")]
    pub fn map_number(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"mapNumber".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=degree,unchecked_return_type="number | null")]
    pub fn degree(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"degree".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=maximumSupported,unchecked_return_type="number | null")]
    pub fn maximum_supported(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"maximumSupported".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn bond(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"bond".into()).unwrap_or(JsValue::NULL)
    }
    #[wasm_bindgen(getter,js_name=order,unchecked_return_type="number | null")]
    pub fn order(&self) -> JsValue {
        js_sys::Reflect::get(&self.fields, &"order".into()).unwrap_or(JsValue::NULL)
    }
}
pub(crate) fn hash_error(source: &ck::MoleculeHashError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: use ck::MoleculeHashError as E;
    use ck::MoleculeHashError as E;
    let kind = match source {
        E::EmptyMolecule => "EmptyMolecule",
        E::MissingPreparedValence => "MissingPreparedValence",
        E::CipRanks(_) => "CipRanks",
        E::InvalidTopology(_) => "InvalidTopology",
        E::RankCount { .. } => "RankCount",
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("MoleculeHashError");
    let e: JsValue = e.into();
    let fields = js_sys::Object::new();
    for name in ["actual", "atomCount"] {
        set(&fields, name, JsValue::NULL)?;
    }
    if let E::RankCount { actual, atom_count } = source {
        set(&fields, "actual", (*actual as u32).into())?;
        set(&fields, "atomCount", (*atom_count as u32).into())?;
        set(&e, "actual", (*actual as u32).into())?;
        set(&e, "atomCount", (*atom_count as u32).into())?;
    }
    set(&e, "domain", "molecular_hash".into())?;
    set(&e, "kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        MoleculeHashError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
pub(crate) fn cip_error(source: &ck::CipRankError) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: use ck::CipRankError as E;
    use ck::CipRankError as E;
    let fields = js_sys::Object::new();
    for name in [
        "field",
        "actual",
        "atomCount",
        "atom",
        "value",
        "mapNumber",
        "degree",
        "maximumSupported",
        "bond",
        "order",
    ] {
        set(&fields, name, JsValue::NULL)?;
    }
    let kind = match source {
        E::InvalidTopology(_) => "InvalidTopology",
        E::InvalidQueryState(_) => "InvalidQueryState",
        E::ValenceRowCount {
            field,
            actual,
            atom_count,
        } => {
            set(&fields, "field", (*field).into())?;
            set(&fields, "actual", (*actual as u32).into())?;
            set(&fields, "atomCount", (*atom_count as u32).into())?;
            "ValenceRowCount"
        }
        E::NegativeImplicitHydrogen { atom, value } => {
            set(&fields, "atom", (atom.index() as u32).into())?;
            set(&fields, "value", (*value).into())?;
            "NegativeImplicitHydrogen"
        }
        E::AtomMapOutOfRange { atom, map_number } => {
            set(&fields, "atom", (atom.index() as u32).into())?;
            set(&fields, "mapNumber", (*map_number).into())?;
            "AtomMapOutOfRange"
        }
        E::InvariantCount { actual, atom_count } => {
            set(&fields, "actual", (*actual as u32).into())?;
            set(&fields, "atomCount", (*atom_count as u32).into())?;
            "InvariantCount"
        }
        E::InvariantOutOfRange { atom, value } => {
            set(&fields, "atom", (atom.index() as u32).into())?;
            set(&fields, "value", JsValue::from(*value))?;
            "InvariantOutOfRange"
        }
        E::TooManyNeighbors {
            atom,
            degree,
            maximum_supported,
        } => {
            set(&fields, "atom", (atom.index() as u32).into())?;
            set(&fields, "degree", (*degree as u32).into())?;
            set(
                &fields,
                "maximumSupported",
                (*maximum_supported as u32).into(),
            )?;
            "TooManyNeighbors"
        }
        E::UnsupportedBondOrder { bond, order } => {
            set(&fields, "bond", (bond.index() as u32).into())?;
            set(&fields, "order", order.rdkit_code().into())?;
            "UnsupportedBondOrder"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("CipRankError");
    let e: JsValue = e.into();
    set(&e, "domain", "cip_ranking".into())?;
    set(&e, "kind", kind.into())?;
    for key in js_sys::Object::keys(&fields).iter() {
        let value = js_sys::Reflect::get(&fields, &key)?;
        if !value.is_null() {
            set(&e, &key.as_string().unwrap(), value)?;
        }
    }
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        CipRankError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=molecularHash)]
    pub fn molecular_hash(&self) -> Result<u64, JsValue> {
        // COSMolKit❗✔️: self.inner.molecular_hash()
        self.inner
            .molecular_hash()
            .map_err(|e| hash_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=molecularHashWithRanks)]
    pub fn molecular_hash_with_ranks(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] ranks: JsValue,
    ) -> Result<u64, JsValue> {
        // COSMolKit❗✔️: self.inner.molecular_hash_with_ranks(ranks)
        let ranks = sequence(&ranks, "ranks")?
            .iter()
            .map(|v| u32_value(&v, "ranks[]"))
            .collect::<Result<Vec<_>, _>>()?;
        self.inner
            .molecular_hash_with_ranks(&ranks)
            .map_err(|e| hash_error(&e).unwrap_or_else(|e| e))
    }
}
