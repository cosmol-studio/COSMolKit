//! Immutable canonical descriptor outputs and the source option vocabulary.
use crate::host_values::u32_value;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum RotatableBondsOptions {
    Default,
    NonStrict,
    Strict,
    StrictLinkages,
}
pub(crate) fn rotatable_option(v: &JsValue) -> Result<ck::RotatableBondsOptions, JsValue> {
    // COSMolKit❗✔️: Default,
    // COSMolKit❗✔️: NonStrict,
    // COSMolKit❗✔️: Strict,
    // COSMolKit❗✔️: StrictLinkages,
    Ok(match u32_value(v, "RotatableBondsOptions")? {
        0 => ck::RotatableBondsOptions::Default,
        1 => ck::RotatableBondsOptions::NonStrict,
        2 => ck::RotatableBondsOptions::Strict,
        3 => ck::RotatableBondsOptions::StrictLinkages,
        _ => return Err(js_sys::RangeError::new("invalid RotatableBondsOptions").into()),
    })
}
#[wasm_bindgen]
pub struct CrippenTotals {
    pub(crate) inner: ck::CrippenTotals,
}
#[wasm_bindgen]
impl CrippenTotals {
    #[wasm_bindgen(getter)]
    pub fn logp(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.logp
        self.inner.logp
    }
    #[wasm_bindgen(getter,js_name=molarRefractivity)]
    pub fn molar_refractivity(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.molar_refractivity
        self.inner.molar_refractivity
    }
}
#[wasm_bindgen]
pub struct LabuteAsaContributions {
    pub(crate) inner: ck::LabuteAsaContributions,
}
#[wasm_bindgen]
impl LabuteAsaContributions {
    #[wasm_bindgen(getter)]
    pub fn asa(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.asa
        self.inner.asa
    }
    #[wasm_bindgen(getter,js_name=hydrogenContribution)]
    pub fn hydrogen_contribution(&self) -> f64 {
        // COSMolKit❗✔️: self.inner.hydrogen_contribution
        self.inner.hydrogen_contribution
    }
    #[wasm_bindgen(getter,js_name=atomContributions)]
    pub fn atom_contributions(&self) -> Vec<f64> {
        // COSMolKit❗✔️: self.inner.atom_contributions.clone()
        self.inner.atom_contributions.clone()
    }
}
