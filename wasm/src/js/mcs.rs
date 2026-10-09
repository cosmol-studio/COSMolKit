//! MCS host transport and live configuration views; no chemistry implementation.
use crate::{Molecule, host_values::*, query_values::QueryGraph};
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use std::{cell::RefCell, rc::Rc, sync::Arc};
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum McsAtomComparator {
    Any,
    Elements,
    Isotopes,
    AnyHeavyAtom,
}
impl McsAtomComparator {
    fn core(self) -> ck::McsAtomComparator {
        match self {
            Self::Any => ck::McsAtomComparator::AtomCompareAny,
            Self::Elements => ck::McsAtomComparator::AtomCompareElements,
            Self::Isotopes => ck::McsAtomComparator::AtomCompareIsotopes,
            Self::AnyHeavyAtom => ck::McsAtomComparator::AtomCompareAnyHeavyAtom,
        }
    }
    fn from_core(value: ck::McsAtomComparator) -> Self {
        match value {
            ck::McsAtomComparator::AtomCompareAny => Self::Any,
            ck::McsAtomComparator::AtomCompareElements => Self::Elements,
            ck::McsAtomComparator::AtomCompareIsotopes => Self::Isotopes,
            ck::McsAtomComparator::AtomCompareAnyHeavyAtom => Self::AnyHeavyAtom,
        }
    }
    fn from_value(value: &JsValue) -> Result<Self, JsValue> {
        match u32_value(value, "McsAtomComparator")? {
            0 => Ok(Self::Any),
            1 => Ok(Self::Elements),
            2 => Ok(Self::Isotopes),
            3 => Ok(Self::AnyHeavyAtom),
            _ => Err(type_error("McsAtomComparator")),
        }
    }
}
#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum McsBondComparator {
    Any,
    Order,
    OrderExact,
}
impl McsBondComparator {
    fn core(self) -> ck::McsBondComparator {
        match self {
            Self::Any => ck::McsBondComparator::BondCompareAny,
            Self::Order => ck::McsBondComparator::BondCompareOrder,
            Self::OrderExact => ck::McsBondComparator::BondCompareOrderExact,
        }
    }
    fn from_core(value: ck::McsBondComparator) -> Self {
        match value {
            ck::McsBondComparator::BondCompareAny => Self::Any,
            ck::McsBondComparator::BondCompareOrder => Self::Order,
            ck::McsBondComparator::BondCompareOrderExact => Self::OrderExact,
        }
    }
    fn from_value(value: &JsValue) -> Result<Self, JsValue> {
        match u32_value(value, "McsBondComparator")? {
            0 => Ok(Self::Any),
            1 => Ok(Self::Order),
            2 => Ok(Self::OrderExact),
            _ => Err(type_error("McsBondComparator")),
        }
    }
}

#[wasm_bindgen]
pub struct McsAtomCompareParameters {
    inner: Rc<RefCell<ck::McsParameters>>,
}
#[cosmolkit_wasm::javascript_options("McsAtomCompareParametersOptions")]
#[wasm_bindgen]
impl McsAtomCompareParameters {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_valences: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_chiral_tag: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_formal_charge: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ring_matches_ring_only: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] complete_rings_only: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_isotope: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_distance: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::McsParameters::default();
        if !match_valences.is_undefined() {
            inner.atom_compare_parameters.match_valences =
                bool_value(&match_valences, "matchValences")?;
        }
        if !match_chiral_tag.is_undefined() {
            inner.atom_compare_parameters.match_chiral_tag =
                bool_value(&match_chiral_tag, "matchChiralTag")?;
        }
        if !match_formal_charge.is_undefined() {
            inner.atom_compare_parameters.match_formal_charge =
                bool_value(&match_formal_charge, "matchFormalCharge")?;
        }
        if !ring_matches_ring_only.is_undefined() {
            inner.atom_compare_parameters.ring_matches_ring_only =
                bool_value(&ring_matches_ring_only, "ringMatchesRingOnly")?;
        }
        if !complete_rings_only.is_undefined() {
            inner.atom_compare_parameters.complete_rings_only =
                bool_value(&complete_rings_only, "completeRingsOnly")?;
        }
        if !match_isotope.is_undefined() {
            inner.atom_compare_parameters.match_isotope =
                bool_value(&match_isotope, "matchIsotope")?;
        }
        if !max_distance.is_undefined() {
            inner.atom_compare_parameters.max_distance = max_distance
                .as_f64()
                .ok_or_else(|| type_error("maxDistance"))?;
        }
        Ok(Self {
            inner: Rc::new(RefCell::new(inner)),
        })
    }
    #[wasm_bindgen(getter, js_name=matchValences)]
    pub fn match_valences(&self) -> bool {
        self.inner.borrow().atom_compare_parameters.match_valences
    }
    #[wasm_bindgen(setter, js_name=matchValences)]
    pub fn set_match_valences(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchValences")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .match_valences = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchChiralTag)]
    pub fn match_chiral_tag(&self) -> bool {
        self.inner.borrow().atom_compare_parameters.match_chiral_tag
    }
    #[wasm_bindgen(setter, js_name=matchChiralTag)]
    pub fn set_match_chiral_tag(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchChiralTag")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .match_chiral_tag = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchFormalCharge)]
    pub fn match_formal_charge(&self) -> bool {
        self.inner
            .borrow()
            .atom_compare_parameters
            .match_formal_charge
    }
    #[wasm_bindgen(setter, js_name=matchFormalCharge)]
    pub fn set_match_formal_charge(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchFormalCharge")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .match_formal_charge = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=ringMatchesRingOnly)]
    pub fn ring_matches_ring_only(&self) -> bool {
        self.inner
            .borrow()
            .atom_compare_parameters
            .ring_matches_ring_only
    }
    #[wasm_bindgen(setter, js_name=ringMatchesRingOnly)]
    pub fn set_ring_matches_ring_only(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "ringMatchesRingOnly")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .ring_matches_ring_only = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=completeRingsOnly)]
    pub fn complete_rings_only(&self) -> bool {
        self.inner
            .borrow()
            .atom_compare_parameters
            .complete_rings_only
    }
    #[wasm_bindgen(setter, js_name=completeRingsOnly)]
    pub fn set_complete_rings_only(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "completeRingsOnly")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .complete_rings_only = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchIsotope)]
    pub fn match_isotope(&self) -> bool {
        self.inner.borrow().atom_compare_parameters.match_isotope
    }
    #[wasm_bindgen(setter, js_name=matchIsotope)]
    pub fn set_match_isotope(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchIsotope")?;
        self.inner
            .borrow_mut()
            .atom_compare_parameters
            .match_isotope = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=maxDistance)]
    pub fn max_distance(&self) -> f64 {
        self.inner.borrow().atom_compare_parameters.max_distance
    }
    #[wasm_bindgen(setter, js_name=maxDistance)]
    pub fn set_max_distance(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = value.as_f64().ok_or_else(|| type_error("maxDistance"))?;
        self.inner.borrow_mut().atom_compare_parameters.max_distance = value;
        Ok(())
    }
}
impl McsAtomCompareParameters {
    fn snapshot(&self) -> ck::McsAtomCompareParameters {
        self.inner.borrow().atom_compare_parameters.clone()
    }
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        if visit_mcsatomcompareparameters(value, &mut |params: &McsAtomCompareParameters| {
            inner = Some(params.inner.borrow().clone())
        })
        .is_ok()
        {
            return inner
                .map(|inner| Self {
                    inner: Rc::new(RefCell::new(inner)),
                })
                .ok_or_else(|| type_error("params"));
        }
        Self::from_js_options(value)
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMcsAtomCompareParameters(v,f){try{f(v);}catch(cause){throw new TypeError('invalid McsAtomCompareParameters',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitMcsAtomCompareParameters)]
    fn visit_mcsatomcompareparameters(
        v: &JsValue,
        f: &mut dyn FnMut(&McsAtomCompareParameters),
    ) -> Result<(), JsValue>;
}

#[wasm_bindgen]
pub struct McsBondCompareParameters {
    inner: Rc<RefCell<ck::McsParameters>>,
}
#[cosmolkit_wasm::javascript_options("McsBondCompareParametersOptions")]
#[wasm_bindgen]
impl McsBondCompareParameters {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] ring_matches_ring_only: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] complete_rings_only: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_fused_rings: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        match_fused_rings_strict: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] match_stereo: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::McsParameters::default();
        if !ring_matches_ring_only.is_undefined() {
            inner.bond_compare_parameters.ring_matches_ring_only =
                bool_value(&ring_matches_ring_only, "ringMatchesRingOnly")?;
        }
        if !complete_rings_only.is_undefined() {
            inner.bond_compare_parameters.complete_rings_only =
                bool_value(&complete_rings_only, "completeRingsOnly")?;
        }
        if !match_fused_rings.is_undefined() {
            inner.bond_compare_parameters.match_fused_rings =
                bool_value(&match_fused_rings, "matchFusedRings")?;
        }
        if !match_fused_rings_strict.is_undefined() {
            inner.bond_compare_parameters.match_fused_rings_strict =
                bool_value(&match_fused_rings_strict, "matchFusedRingsStrict")?;
        }
        if !match_stereo.is_undefined() {
            inner.bond_compare_parameters.match_stereo = bool_value(&match_stereo, "matchStereo")?;
        }
        Ok(Self {
            inner: Rc::new(RefCell::new(inner)),
        })
    }
    #[wasm_bindgen(getter, js_name=ringMatchesRingOnly)]
    pub fn ring_matches_ring_only(&self) -> bool {
        self.inner
            .borrow()
            .bond_compare_parameters
            .ring_matches_ring_only
    }
    #[wasm_bindgen(setter, js_name=ringMatchesRingOnly)]
    pub fn set_ring_matches_ring_only(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "ringMatchesRingOnly")?;
        self.inner
            .borrow_mut()
            .bond_compare_parameters
            .ring_matches_ring_only = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=completeRingsOnly)]
    pub fn complete_rings_only(&self) -> bool {
        self.inner
            .borrow()
            .bond_compare_parameters
            .complete_rings_only
    }
    #[wasm_bindgen(setter, js_name=completeRingsOnly)]
    pub fn set_complete_rings_only(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "completeRingsOnly")?;
        self.inner
            .borrow_mut()
            .bond_compare_parameters
            .complete_rings_only = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchFusedRings)]
    pub fn match_fused_rings(&self) -> bool {
        self.inner
            .borrow()
            .bond_compare_parameters
            .match_fused_rings
    }
    #[wasm_bindgen(setter, js_name=matchFusedRings)]
    pub fn set_match_fused_rings(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchFusedRings")?;
        self.inner
            .borrow_mut()
            .bond_compare_parameters
            .match_fused_rings = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchFusedRingsStrict)]
    pub fn match_fused_rings_strict(&self) -> bool {
        self.inner
            .borrow()
            .bond_compare_parameters
            .match_fused_rings_strict
    }
    #[wasm_bindgen(setter, js_name=matchFusedRingsStrict)]
    pub fn set_match_fused_rings_strict(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchFusedRingsStrict")?;
        self.inner
            .borrow_mut()
            .bond_compare_parameters
            .match_fused_rings_strict = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=matchStereo)]
    pub fn match_stereo(&self) -> bool {
        self.inner.borrow().bond_compare_parameters.match_stereo
    }
    #[wasm_bindgen(setter, js_name=matchStereo)]
    pub fn set_match_stereo(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "matchStereo")?;
        self.inner.borrow_mut().bond_compare_parameters.match_stereo = value;
        Ok(())
    }
}
impl McsBondCompareParameters {
    fn snapshot(&self) -> ck::McsBondCompareParameters {
        self.inner.borrow().bond_compare_parameters.clone()
    }
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        if visit_mcsbondcompareparameters(value, &mut |params: &McsBondCompareParameters| {
            inner = Some(params.inner.borrow().clone())
        })
        .is_ok()
        {
            return inner
                .map(|inner| Self {
                    inner: Rc::new(RefCell::new(inner)),
                })
                .ok_or_else(|| type_error("params"));
        }
        Self::from_js_options(value)
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMcsBondCompareParameters(v,f){try{f(v);}catch(cause){throw new TypeError('invalid McsBondCompareParameters',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitMcsBondCompareParameters)]
    fn visit_mcsbondcompareparameters(
        v: &JsValue,
        f: &mut dyn FnMut(&McsBondCompareParameters),
    ) -> Result<(), JsValue>;
}

#[wasm_bindgen]
pub struct McsParameters {
    inner: Rc<RefCell<ck::McsParameters>>,
}
#[cosmolkit_wasm::javascript_options("McsParametersOptions")]
#[wasm_bindgen]
impl McsParameters {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] store_all: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] maximize_bonds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] timeout: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] verbose: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "McsAtomCompareParameters | McsAtomCompareParametersOptions"
        )]
        atom_compare_parameters: JsValue,
        #[wasm_bindgen(
            unchecked_optional_param_type = "McsBondCompareParameters | McsBondCompareParametersOptions"
        )]
        bond_compare_parameters: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "McsAtomComparator")]
        atom_comparator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "McsBondComparator")]
        bond_comparator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string")] initial_seed: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::McsParameters::default();
        if !store_all.is_undefined() {
            inner.store_all = bool_value(&store_all, "storeAll")?;
        }
        if !maximize_bonds.is_undefined() {
            inner.maximize_bonds = bool_value(&maximize_bonds, "maximizeBonds")?;
        }
        if !threshold.is_undefined() {
            inner.threshold = threshold.as_f64().ok_or_else(|| type_error("threshold"))?;
        }
        if !timeout.is_undefined() {
            inner.timeout = u32_value(&timeout, "timeout")?;
        }
        if !verbose.is_undefined() {
            inner.verbose = bool_value(&verbose, "verbose")?;
        }
        if !atom_compare_parameters.is_undefined() {
            inner.atom_compare_parameters =
                McsAtomCompareParameters::from_configuration(&atom_compare_parameters)?.snapshot();
        }
        if !bond_compare_parameters.is_undefined() {
            inner.bond_compare_parameters =
                McsBondCompareParameters::from_configuration(&bond_compare_parameters)?.snapshot();
        }
        if !atom_comparator.is_undefined() {
            inner.atom_comparator = McsAtomComparator::from_value(&atom_comparator)?.core();
        }
        if !bond_comparator.is_undefined() {
            inner.bond_comparator = McsBondComparator::from_value(&bond_comparator)?.core();
        }
        if !initial_seed.is_undefined() {
            inner.initial_seed = initial_seed
                .as_string()
                .ok_or_else(|| type_error("initialSeed"))?;
        }
        Ok(Self {
            inner: Rc::new(RefCell::new(inner)),
        })
    }
    #[wasm_bindgen(getter, js_name=storeAll)]
    pub fn store_all(&self) -> bool {
        self.inner.borrow().store_all
    }
    #[wasm_bindgen(setter, js_name=storeAll)]
    pub fn set_store_all(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "storeAll")?;
        self.inner.borrow_mut().store_all = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=maximizeBonds)]
    pub fn maximize_bonds(&self) -> bool {
        self.inner.borrow().maximize_bonds
    }
    #[wasm_bindgen(setter, js_name=maximizeBonds)]
    pub fn set_maximize_bonds(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "maximizeBonds")?;
        self.inner.borrow_mut().maximize_bonds = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=threshold)]
    pub fn threshold(&self) -> f64 {
        self.inner.borrow().threshold
    }
    #[wasm_bindgen(setter, js_name=threshold)]
    pub fn set_threshold(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = value.as_f64().ok_or_else(|| type_error("threshold"))?;
        self.inner.borrow_mut().threshold = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=timeout)]
    pub fn timeout(&self) -> u32 {
        self.inner.borrow().timeout
    }
    #[wasm_bindgen(setter, js_name=timeout)]
    pub fn set_timeout(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = u32_value(&value, "timeout")?;
        self.inner.borrow_mut().timeout = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=verbose)]
    pub fn verbose(&self) -> bool {
        self.inner.borrow().verbose
    }
    #[wasm_bindgen(setter, js_name=verbose)]
    pub fn set_verbose(&mut self, value: JsValue) -> Result<(), JsValue> {
        let value = bool_value(&value, "verbose")?;
        self.inner.borrow_mut().verbose = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=atomCompareParameters)]
    pub fn atom_compare_parameters(&self) -> McsAtomCompareParameters {
        McsAtomCompareParameters {
            inner: self.inner.clone(),
        }
    }
    #[wasm_bindgen(setter, js_name=atomCompareParameters)]
    pub fn set_atom_compare_parameters(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "McsAtomCompareParameters | McsAtomCompareParametersOptions"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        let value = McsAtomCompareParameters::from_configuration(&value)?.snapshot();
        self.inner.borrow_mut().atom_compare_parameters = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=bondCompareParameters)]
    pub fn bond_compare_parameters(&self) -> McsBondCompareParameters {
        McsBondCompareParameters {
            inner: self.inner.clone(),
        }
    }
    #[wasm_bindgen(setter, js_name=bondCompareParameters)]
    pub fn set_bond_compare_parameters(
        &mut self,
        #[wasm_bindgen(
            unchecked_param_type = "McsBondCompareParameters | McsBondCompareParametersOptions"
        )]
        value: JsValue,
    ) -> Result<(), JsValue> {
        let value = McsBondCompareParameters::from_configuration(&value)?.snapshot();
        self.inner.borrow_mut().bond_compare_parameters = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=atomComparator)]
    pub fn atom_comparator(&self) -> McsAtomComparator {
        McsAtomComparator::from_core(self.inner.borrow().atom_comparator)
    }
    #[wasm_bindgen(setter, js_name=atomComparator)]
    pub fn set_atom_comparator(&mut self, value: McsAtomComparator) -> Result<(), JsValue> {
        let value = value.core();
        self.inner.borrow_mut().atom_comparator = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=bondComparator)]
    pub fn bond_comparator(&self) -> McsBondComparator {
        McsBondComparator::from_core(self.inner.borrow().bond_comparator)
    }
    #[wasm_bindgen(setter, js_name=bondComparator)]
    pub fn set_bond_comparator(&mut self, value: McsBondComparator) -> Result<(), JsValue> {
        let value = value.core();
        self.inner.borrow_mut().bond_comparator = value;
        Ok(())
    }
    #[wasm_bindgen(getter, js_name=initialSeed)]
    pub fn initial_seed(&self) -> String {
        self.inner.borrow().initial_seed.clone()
    }
    #[wasm_bindgen(setter, js_name=initialSeed)]
    pub fn set_initial_seed(&mut self, value: String) -> Result<(), JsValue> {
        let value = value;
        self.inner.borrow_mut().initial_seed = value;
        Ok(())
    }
}
impl McsParameters {
    fn snapshot(&self) -> ck::McsParameters {
        self.inner.borrow().clone()
    }
    fn from_configuration(value: &JsValue) -> Result<Self, JsValue> {
        let mut inner = None;
        if visit_mcsparameters(value, &mut |params: &McsParameters| {
            inner = Some(params.inner.borrow().clone())
        })
        .is_ok()
        {
            return inner
                .map(|inner| Self {
                    inner: Rc::new(RefCell::new(inner)),
                })
                .ok_or_else(|| type_error("params"));
        }
        Self::from_js_options(value)
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMcsParameters(v,f){try{f(v);}catch(cause){throw new TypeError('invalid McsParameters',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitMcsParameters)]
    fn visit_mcsparameters(v: &JsValue, f: &mut dyn FnMut(&McsParameters)) -> Result<(), JsValue>;
}

#[wasm_bindgen]
pub struct McsResult {
    inner: ck::McsResult,
}
#[wasm_bindgen]
impl McsResult {
    #[wasm_bindgen(getter, unchecked_return_type = "QueryGraph | null")]
    pub fn query(&self) -> JsValue {
        self.inner
            .query
            .clone()
            .map_or(JsValue::NULL, |inner| QueryGraph { inner }.into())
    }
    #[wasm_bindgen(getter, js_name=atomCount)]
    pub fn atom_count(&self) -> usize {
        self.inner.atom_count
    }
    #[wasm_bindgen(getter, js_name=bondCount)]
    pub fn bond_count(&self) -> usize {
        self.inner.bond_count
    }
    #[wasm_bindgen(getter)]
    pub fn completed(&self) -> bool {
        self.inner.completed
    }
    #[wasm_bindgen(getter)]
    pub fn smarts(&self) -> Result<String, JsValue> {
        text(&self.inner.smarts)
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Map<string, QueryGraph>")]
    pub fn degenerate(&self) -> Result<Map, JsValue> {
        let results = Map::new();
        for (key, query) in &self.inner.degenerate {
            results.set(
                &text(key)?.into(),
                &QueryGraph {
                    inner: query.clone(),
                }
                .into(),
            );
        }
        Ok(results)
    }
}
#[wasm_bindgen]
pub struct McsError {
    inner: ck::McsError,
}
#[wasm_bindgen]
impl McsError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "mcs".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        match &self.inner {
            ck::McsError::State(_) => "State",
            ck::McsError::Progress(_) => "Progress",
            ck::McsError::StereoOrder(_) => "StereoOrder",
            ck::McsError::TargetTableCount { .. } => "TargetTableCount",
            ck::McsError::ThresholdCountOutOfRange { .. } => "ThresholdCountOutOfRange",
            ck::McsError::MatchTableOutOfRange { .. } => "MatchTableOutOfRange",
            ck::McsError::RingMembershipMissing { .. } => "RingMembershipMissing",
            ck::McsError::RingOutOfRange { .. } => "RingOutOfRange",
            ck::McsError::MappedBondMissing { .. } => "MappedBondMissing",
            ck::McsError::TooManyNewBonds { .. } => "TooManyNewBonds",
            ck::McsError::InitialSeedParse { .. } => "InitialSeedParse",
            ck::McsError::InitialSeedMatch { .. } => "InitialSeedMatch",
            ck::McsError::InitialSeedBondMissing { .. } => "InitialSeedBondMissing",
            ck::McsError::ResultValueOutOfRange { .. } => "ResultValueOutOfRange",
            ck::McsError::ResultQueryGraph { .. } => "ResultQueryGraph",
            ck::McsError::ResultSmarts { .. } => "ResultSmarts",
            ck::McsError::SeedReconstructionBondMissing { .. } => "SeedReconstructionBondMissing",
            ck::McsError::ResultContextMissing => "ResultContextMissing",
        }
        .into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
}
fn mcs_error(inner: ck::McsError) -> JsValue {
    let detail = McsError { inner };
    let error = js_sys::Error::new(&detail.message());
    error.set_name("McsError");
    let set = |key: &str, value: JsValue| js_sys::Reflect::set(&error, &key.into(), &value);
    if let Err(error) = set("domain", detail.domain().into())
        .and_then(|_| set("kind", detail.kind().into()))
        .and_then(|_| set("detail", detail.into()))
    {
        return error;
    }
    error.into()
}
#[wasm_bindgen(
    inline_js = "export function visitMcsMolecule(v,f){try{f(v);}catch(cause){throw new TypeError('invalid Molecule',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch, js_name=visitMcsMolecule)]
    fn visit_molecule(v: &JsValue, f: &mut dyn FnMut(&Molecule)) -> Result<(), JsValue>;
}
fn inputs(value: &JsValue) -> Result<Vec<Arc<cosmolkit_wasm::Molecule>>, JsValue> {
    sequence(value, "inputs")?
        .iter()
        .map(|value| {
            let mut input = None;
            visit_molecule(&value, &mut |mol: &Molecule| {
                input = Some(mol.inner.clone())
            })?;
            input.ok_or_else(|| type_error("input molecule"))
        })
        .collect()
}
#[wasm_bindgen(js_name=maximumCommonSubstructure)]
pub fn maximum_common_substructure(
    #[wasm_bindgen(unchecked_param_type = "Molecule[]")] molecules: JsValue,
    #[wasm_bindgen(unchecked_optional_param_type="McsParameters | McsParametersOptions")] params: JsValue,
) -> Result<McsResult, JsValue> {
    let inputs = inputs(&molecules)?;
    let references: Vec<_> = inputs.iter().map(|mol| &**mol).collect();
    let result = if params.is_undefined() {
        cosmolkit_wasm::maximum_common_substructure(&references)
    } else {
        cosmolkit_wasm::maximum_common_substructure_with_params(
            &references,
            &McsParameters::from_configuration(&params)?.snapshot(),
        )
    };
    result.map(|inner| McsResult { inner }).map_err(mcs_error)
}
#[wasm_bindgen(js_name=maximumCommonSubstructureWithParams)]
pub fn maximum_common_substructure_with_params(
    #[wasm_bindgen(unchecked_param_type = "Molecule[]")] molecules: JsValue,
    params: &McsParameters,
) -> Result<McsResult, JsValue> {
    let inputs = inputs(&molecules)?;
    let references: Vec<_> = inputs.iter().map(|mol| &**mol).collect();
    cosmolkit_wasm::maximum_common_substructure_with_params(&references, &params.snapshot())
        .map(|inner| McsResult { inner })
        .map_err(mcs_error)
}
