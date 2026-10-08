//! Complete ring policies and canonical operations, without cache or algorithm authority.
use crate::Molecule;
use crate::alignment_values::operation_error;
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct RingSearchParams {
    inner: ck::RingSearchParams,
}
#[wasm_bindgen]
impl RingSearchParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_dative_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_hydrogen_bonds: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::RingSearchParams::default();
        if !include_dative_bonds.is_undefined() {
            inner.include_dative_bonds = bool_value(&include_dative_bonds, "includeDativeBonds")?;
        }
        if !include_hydrogen_bonds.is_undefined() {
            inner.include_hydrogen_bonds =
                bool_value(&include_hydrogen_bonds, "includeHydrogenBonds")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=includeDativeBonds)]
    pub fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[wasm_bindgen(getter,js_name=includeHydrogenBonds)]
    pub fn include_hydrogen_bonds(&self) -> bool {
        self.inner.include_hydrogen_bonds
    }
}
#[wasm_bindgen(
    inline_js = "export function visitRingSearchParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid RingSearchParams',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitRingSearchParams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&RingSearchParams)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=withAssignedRings)]
    pub fn with_assigned_rings(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_assigned_rings()
        self.inner
            .with_assigned_rings()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=assignRings)]
    pub fn assign_rings_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assign_rings_()
        self.inner
            .assign_rings_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withAssignedRingFamilies)]
    pub fn with_assigned_ring_families(&self) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_assigned_ring_families()
        self.inner
            .with_assigned_ring_families()
            .map(|inner| Self {
                inner: Arc::new(inner),
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=assignRingFamilies)]
    pub fn assign_ring_families_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assign_ring_families_()
        self.inner
            .assign_ring_families_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withAssignedRingFamiliesWithParams)]
    pub fn with_assigned_ring_families_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "RingSearchParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: self.inner.with_assigned_ring_families_with_params(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &RingSearchParams| {
            result = Some(
                self.inner
                    .with_assigned_ring_families_with_params(&p.inner)
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=assignRingFamiliesWithParams)]
    pub fn assign_ring_families_with_params_(
        &self,
        #[wasm_bindgen(unchecked_param_type = "RingSearchParams")] params: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assign_ring_families_with_params_(&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &RingSearchParams| {
            result = Some(
                self.inner
                    .assign_ring_families_with_params_(&p.inner)
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
