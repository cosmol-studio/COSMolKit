//! Complete canonical ValenceError vocabulary and source fields.
use crate::alignment_values::{set, source_error};
use cosmolkit_wasm::rust as ck;
use std::error::Error as RustError;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct ValenceError {
    kind: String,
    message: String,
    cause: JsValue,
    fields: js_sys::Object,
}
#[wasm_bindgen]
impl ValenceError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "valence".into()
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
    #[wasm_bindgen(getter,js_name=atom,unchecked_return_type="number | null")]
    pub fn field_0(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atom".into())
    }
    #[wasm_bindgen(getter,js_name=explicitValence,unchecked_return_type="number | null")]
    pub fn field_1(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"explicitValence".into())
    }
    #[wasm_bindgen(getter,js_name=physicalBonds,unchecked_return_type="number | null")]
    pub fn field_2(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"physicalBonds".into())
    }
    #[wasm_bindgen(getter,js_name=atomicNumber,unchecked_return_type="number | null")]
    pub fn field_3(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atomicNumber".into())
    }
    #[wasm_bindgen(getter,js_name=formalCharge,unchecked_return_type="number | null")]
    pub fn field_4(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"formalCharge".into())
    }
    #[wasm_bindgen(getter,js_name=phase,unchecked_return_type="string | null")]
    pub fn field_5(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"phase".into())
    }
    #[wasm_bindgen(getter,js_name=calculated,unchecked_return_type="number | null")]
    pub fn field_6(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"calculated".into())
    }
    #[wasm_bindgen(getter,js_name=reason,unchecked_return_type="string | null")]
    pub fn field_7(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"reason".into())
    }
    #[wasm_bindgen(getter,js_name=atomCount,unchecked_return_type="number | null")]
    pub fn field_8(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"atomCount".into())
    }
    #[wasm_bindgen(getter,js_name=neighborAtom,unchecked_return_type="number | null")]
    pub fn field_9(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"neighborAtom".into())
    }
    #[wasm_bindgen(getter,js_name=bond,unchecked_return_type="number | null")]
    pub fn field_10(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"bond".into())
    }
    #[wasm_bindgen(getter,js_name=bondCount,unchecked_return_type="number | null")]
    pub fn field_11(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"bondCount".into())
    }
    #[wasm_bindgen(getter,js_name=begin,unchecked_return_type="number | null")]
    pub fn field_12(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"begin".into())
    }
    #[wasm_bindgen(getter,js_name=end,unchecked_return_type="number | null")]
    pub fn field_13(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"end".into())
    }
    #[wasm_bindgen(getter,js_name=value,unchecked_return_type="number | null")]
    pub fn field_14(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"value".into())
    }
    #[wasm_bindgen(getter,js_name=field,unchecked_return_type="string | null")]
    pub fn field_15(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"field".into())
    }
    #[wasm_bindgen(getter,js_name=explicit,unchecked_return_type="number | null")]
    pub fn field_16(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"explicit".into())
    }
    #[wasm_bindgen(getter,js_name=implicit,unchecked_return_type="number | null")]
    pub fn field_17(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"implicit".into())
    }
    #[wasm_bindgen(getter,js_name=neighborHydrogens,unchecked_return_type="number | null")]
    pub fn field_18(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"neighborHydrogens".into())
    }
    #[wasm_bindgen(getter,js_name=order,unchecked_return_type="BondOrder | null")]
    pub fn field_19(&self) -> Result<JsValue, JsValue> {
        js_sys::Reflect::get(&self.fields, &"order".into())
    }
}
pub(crate) fn valence_error(source: &ck::ValenceError) -> Result<JsValue, JsValue> {
    use ck::ValenceError as E;
    let fields = js_sys::Object::new();
    set(&fields, "atom", JsValue::NULL)?;
    set(&fields, "explicitValence", JsValue::NULL)?;
    set(&fields, "physicalBonds", JsValue::NULL)?;
    set(&fields, "atomicNumber", JsValue::NULL)?;
    set(&fields, "formalCharge", JsValue::NULL)?;
    set(&fields, "phase", JsValue::NULL)?;
    set(&fields, "calculated", JsValue::NULL)?;
    set(&fields, "reason", JsValue::NULL)?;
    set(&fields, "atomCount", JsValue::NULL)?;
    set(&fields, "neighborAtom", JsValue::NULL)?;
    set(&fields, "bond", JsValue::NULL)?;
    set(&fields, "bondCount", JsValue::NULL)?;
    set(&fields, "begin", JsValue::NULL)?;
    set(&fields, "end", JsValue::NULL)?;
    set(&fields, "value", JsValue::NULL)?;
    set(&fields, "field", JsValue::NULL)?;
    set(&fields, "explicit", JsValue::NULL)?;
    set(&fields, "implicit", JsValue::NULL)?;
    set(&fields, "neighborHydrogens", JsValue::NULL)?;
    set(&fields, "order", JsValue::NULL)?;
    let kind = match source {
        E::PiElectronExplicitValenceCacheNotInitialized { atom } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            "PiElectronExplicitValenceCacheNotInitialized"
        }
        E::PiElectronInvariant {
            atom,
            explicit_valence,
            physical_bonds,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(
                &fields,
                "explicitValence",
                JsValue::from(*explicit_valence as u32),
            )?;
            set(
                &fields,
                "physicalBonds",
                JsValue::from(*physical_bonds as u32),
            )?;
            "PiElectronInvariant"
        }
        E::InvalidValence {
            atom,
            atomic_number,
            formal_charge,
            phase,
            calculated,
            reason,
            message,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(
                &fields,
                "atomicNumber",
                JsValue::from(*atomic_number as u32),
            )?;
            set(&fields, "formalCharge", JsValue::from(*formal_charge))?;
            set(&fields, "phase", JsValue::from(format!("{phase:?}")))?;
            set(
                &fields,
                "calculated",
                calculated.map_or(JsValue::NULL, JsValue::from),
            )?;
            set(&fields, "reason", JsValue::from(*reason))?;
            let _ = message;
            "InvalidValence"
        }
        E::InvalidTopology { source } => {
            let _ = source;
            "InvalidTopology"
        }
        E::AtomOutOfRange { atom, atom_count } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(&fields, "atomCount", JsValue::from(*atom_count as u32))?;
            "AtomOutOfRange"
        }
        E::AdjacencyAtomOutOfRange {
            atom,
            neighbor_atom,
            atom_count,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(
                &fields,
                "neighborAtom",
                JsValue::from(*neighbor_atom as u32),
            )?;
            set(&fields, "atomCount", JsValue::from(*atom_count as u32))?;
            "AdjacencyAtomOutOfRange"
        }
        E::AdjacencyBondOutOfRange {
            atom,
            bond,
            bond_count,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(&fields, "bondCount", JsValue::from(*bond_count as u32))?;
            "AdjacencyBondOutOfRange"
        }
        E::AdjacencyEndpointMismatch {
            atom,
            neighbor_atom,
            bond,
            begin,
            end,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(
                &fields,
                "neighborAtom",
                JsValue::from(*neighbor_atom as u32),
            )?;
            set(&fields, "bond", JsValue::from(bond.index() as u32))?;
            set(&fields, "begin", JsValue::from(begin.index() as u32))?;
            set(&fields, "end", JsValue::from(end.index() as u32))?;
            "AdjacencyEndpointMismatch"
        }
        E::InvalidExplicitValenceInput { atom, value } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(&fields, "value", JsValue::from(*value))?;
            "InvalidExplicitValenceInput"
        }
        E::PeriodicTableLookup {
            atomic_number,
            field,
        } => {
            set(
                &fields,
                "atomicNumber",
                JsValue::from(*atomic_number as u32),
            )?;
            set(&fields, "field", JsValue::from(*field))?;
            "PeriodicTableLookup"
        }
        E::ExplicitValenceCacheNotInitialized { atom } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            "ExplicitValenceCacheNotInitialized"
        }
        E::ImplicitValenceCacheNotInitialized { atom } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            "ImplicitValenceCacheNotInitialized"
        }
        E::HydrogenCountOverflow {
            atom,
            explicit,
            implicit,
            neighbor_hydrogens,
        } => {
            set(&fields, "atom", JsValue::from(atom.index() as u32))?;
            set(&fields, "explicit", JsValue::from(*explicit as u32))?;
            set(&fields, "implicit", JsValue::from(*implicit as u32))?;
            set(
                &fields,
                "neighborHydrogens",
                JsValue::from(*neighbor_hydrogens as u32),
            )?;
            "HydrogenCountOverflow"
        }
        E::BadBondType { bond, order } => {
            set(
                &fields,
                "bond",
                bond.map_or(JsValue::NULL, |id| JsValue::from(id.index() as u32)),
            )?;
            set(&fields, "order", JsValue::from(*order as u32))?;
            "BadBondType"
        }
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("ValenceError");
    set(&error, "domain", "valence".into())?;
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
        ValenceError {
            kind: kind.into(),
            message: source.to_string(),
            cause,
            fields,
        }
        .into(),
    )?;
    Ok(error.into())
}
