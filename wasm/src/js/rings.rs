//! Complete ring policies and canonical operations, without cache or algorithm authority.
use crate::Molecule;
use crate::alignment_values::operation_error;
use crate::host_values::*;
use cosmolkit_wasm::rust as ck;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
/// Owned symmetrized SSSR result, independent of subsequent molecule edits.
#[wasm_bindgen]
pub struct RingSearchResult {
    inner: ck::RingSearchResult,
}
#[wasm_bindgen]
impl RingSearchResult {
    #[wasm_bindgen(getter, js_name=atomRings, unchecked_return_type="number[][]")]
    pub fn atom_rings(&self) -> js_sys::Array {
        self.inner.atom_rings().iter().map(|ring| {
            ring.iter().map(|id| JsValue::from(id.index() as u32)).collect::<js_sys::Array>()
        }).collect()
    }
    #[wasm_bindgen(getter, js_name=bondRings, unchecked_return_type="number[][]")]
    pub fn bond_rings(&self) -> js_sys::Array {
        self.inner.bond_rings().iter().map(|ring| {
            ring.iter().map(|id| JsValue::from(id.index() as u32)).collect::<js_sys::Array>()
        }).collect()
    }
    #[wasm_bindgen(getter, js_name=numRings)]
    pub fn num_rings(&self) -> usize { self.inner.num_rings() }
}
#[wasm_bindgen]
pub struct RingSearchParams {
    inner: ck::RingSearchParams,
}
#[cosmolkit_wasm::javascript_options("RingSearchParamsOptions")]
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
    #[wasm_bindgen(js_name=symmetrizedSssr)]
    pub fn symmetrized_sssr(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type="RingSearchParams | RingSearchParamsOptions")]
        params: JsValue,
    ) -> Result<RingSearchResult, JsValue> {
        if params.is_undefined() {
            return self.inner.symmetrized_sssr().map(|inner| RingSearchResult { inner })
                .map_err(|e| operation_error(&e).unwrap_or_else(|e| e));
        }
        let params = RingSearchParams::from_js_configuration(&params)?;
        self.inner.symmetrized_sssr_with_params(&params.inner)
            .map(|inner| RingSearchResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=symmetrizedSssrWithParams)]
    pub fn symmetrized_sssr_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type="RingSearchParams")] params: JsValue,
    ) -> Result<RingSearchResult, JsValue> {
        let mut result = None;
        visit_params(&params, &mut |p: &RingSearchParams| {
            result = Some(self.inner.symmetrized_sssr_with_params(&p.inner)
                .map(|inner| RingSearchResult { inner })
                .map_err(|e| operation_error(&e).unwrap_or_else(|e| e)));
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
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
    #[wasm_bindgen(js_name=assignRings_)]
    pub fn assign_rings_(&self) -> Result<(), JsValue> {
        // COSMolKit❗✔️: self.inner.assign_rings_()
        self.inner
            .assign_rings_()
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withAssignedRingFamilies)]
    pub fn with_assigned_ring_families(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "RingSearchParams | RingSearchParamsOptions"
        )]
        params: JsValue,
    ) -> Result<Self, JsValue> {
        if params.is_undefined() {
            return {
                // COSMolKit❗✔️: self.inner.with_assigned_ring_families()
                self.inner
                    .with_assigned_ring_families()
                    .map(|inner| Self {
                        inner: Arc::new(inner),
                    })
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
            };
        }
        let configuration_0 = crate::rings::RingSearchParams::from_js_configuration(&params)?;
        self.with_assigned_ring_families_with_params(configuration_0.into())
    }
    #[wasm_bindgen(js_name=assignRingFamilies_)]
    pub fn assign_ring_families_(
        &self,
        #[wasm_bindgen(
            unchecked_optional_param_type = "RingSearchParams | RingSearchParamsOptions"
        )]
        params: JsValue,
    ) -> Result<(), JsValue> {
        if params.is_undefined() {
            return {
                // COSMolKit❗✔️: self.inner.assign_ring_families_()
                self.inner
                    .assign_ring_families_()
                    .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
            };
        }
        let configuration_0 = crate::rings::RingSearchParams::from_js_configuration(&params)?;
        self.assign_ring_families_with_params_(configuration_0.into())
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
    #[wasm_bindgen(js_name=assignRingFamiliesWithParams_)]
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
