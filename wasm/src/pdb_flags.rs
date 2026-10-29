use crate::{fields::Field, input};
use cosmolkit as ck;
use tsify::Ts;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct PdbReadFlags(ck::PdbReadFlags);
crate::fields::value_field!(
    ck::PdbReadFlags,
    PdbReadFlags,
    PdbReadFlagsInput,
    "PdbReadFlags"
);
#[wasm_bindgen]
impl PdbReadFlags {
    #[wasm_bindgen(getter,js_name=NONE)]
    pub fn none() -> Self {
        Self(ck::PdbReadFlags::NONE)
    }
    #[wasm_bindgen(getter,js_name=KEEP_ALTERNATE_AND_PSEUDO_ATOMS)]
    pub fn keep_alternate_and_pseudo_atoms() -> Self {
        Self(ck::PdbReadFlags::KEEP_ALTERNATE_AND_PSEUDO_ATOMS)
    }
    #[wasm_bindgen(getter,js_name=STANDARD_RESIDUE_BOND_ORDERS)]
    pub fn standard_residue_bond_orders() -> Self {
        Self(ck::PdbReadFlags::STANDARD_RESIDUE_BOND_ORDERS)
    }
    #[wasm_bindgen(js_name=fromBits)]
    pub fn from_bits(bits: Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::PdbReadFlags::from_bits(input::read::<u32>(bits)?)))
    }
    #[wasm_bindgen(getter)]
    pub fn bits(&self) -> u32 {
        self.0.bits()
    }
    #[wasm_bindgen(getter,js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
    pub fn contains(&self, other: PdbReadFlagsInput) -> Result<bool, JsValue> {
        Self::with_ref(&other, |v| Ok(self.0.contains(*v)))
    }
    pub fn union(&self, other: PdbReadFlagsInput) -> Result<Self, JsValue> {
        Self::with_ref(&other, |v| Ok(Self(self.0 | *v)))
    }
    pub fn intersection(&self, other: PdbReadFlagsInput) -> Result<Self, JsValue> {
        Self::with_ref(&other, |v| Ok(Self(self.0 & *v)))
    }
}
