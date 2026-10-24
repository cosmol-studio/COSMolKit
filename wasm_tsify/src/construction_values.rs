//! Detached construction facts. Native builders own indexing and validation.
use crate::{
    atoms::{Element, ElementInput},
    enums::BondOrderInput,
    fields::Field,
    input,
};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct AtomSpec(pub(crate) ck::AtomSpec);
crate::fields::value_field!(ck::AtomSpec, AtomSpec, AtomSpecInput, "AtomSpec");
#[wasm_bindgen]
impl AtomSpec {
    #[wasm_bindgen(constructor)]
    pub fn new(element: ElementInput) -> Result<Self, JsValue> {
        Element::with_ref(&element, |v| Ok(Self(ck::AtomSpec::new(*v))))
    }
    #[wasm_bindgen(js_name=withFormalCharge)]
    pub fn with_formal_charge(&self, value: Ts<input::I8Input>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_formal_charge(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=withExplicitHydrogens)]
    pub fn with_explicit_hydrogens(&self, value: Ts<input::U8Input>) -> Result<Self, JsValue> {
        Ok(Self(
            self.0.clone().with_explicit_hydrogens(input::read(value)?),
        ))
    }
    #[wasm_bindgen(js_name=withAtomMap)]
    pub fn with_atom_map(&self, value: Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_atom_map(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=withIsotope)]
    pub fn with_isotope(&self, value: Ts<input::U16Input>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_isotope(input::read(value)?)))
    }
    #[wasm_bindgen(js_name=withNoImplicit)]
    pub fn with_no_implicit(&self, value: Ts<input::BoolInput>) -> Result<Self, JsValue> {
        Ok(Self(self.0.clone().with_no_implicit(input::read(value)?)))
    }
}
#[wasm_bindgen]
pub struct BondSpec(pub(crate) ck::BondSpec);
crate::fields::value_field!(ck::BondSpec, BondSpec, BondSpecInput, "BondSpec");
#[wasm_bindgen]
impl BondSpec {
    #[wasm_bindgen(constructor)]
    pub fn new(
        begin: Ts<input::IndexInput>,
        end: Ts<input::IndexInput>,
        order: BondOrderInput,
    ) -> Result<Self, JsValue> {
        Ok(Self(ck::BondSpec::new(
            ck::AtomId::new(input::read(begin)?),
            ck::AtomId::new(input::read(end)?),
            ck::BondOrder::input(order)?,
        )))
    }
}
