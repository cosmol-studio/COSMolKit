//! Complete immutable Morgan generator/provider/call configuration transport.
use crate::atom_pair_parameters::{
    AtomPairAtomInvariantsGenerator, optional_array, optional_u32_array,
};
use crate::host_values::{bool_value, i32_value, sequence, type_error, u32_value};
use crate::query_construction::QueryGraph;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MorganParams {
    pub(crate) inner: ck::MorganParams,
}
#[wasm_bindgen]
impl MorganParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] radius: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_bond_types: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_ring_membership: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] only_nonzero_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        include_redundant_environments: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] count_simulation: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        count_bounds: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bits_per_feature: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_fingerprint_values.rs::MorganParams::new:
        //         // Configuration transport only: validation belongs to canonical callers.
        // O(fields + supplied array length), detached configurations only.
        let mut inner = ck::MorganParams::default();
        if !radius.is_undefined() {
            inner.radius = u32_value(&radius, "radius")?;
        }
        if !include_chirality.is_undefined() {
            inner.include_chirality = bool_value(&include_chirality, "includeChirality")?;
        }
        if !use_bond_types.is_undefined() {
            inner.use_bond_types = bool_value(&use_bond_types, "useBondTypes")?;
        }
        if !include_ring_membership.is_undefined() {
            inner.include_ring_membership =
                bool_value(&include_ring_membership, "includeRingMembership")?;
        }
        if !only_nonzero_invariants.is_undefined() {
            inner.only_nonzero_invariants =
                bool_value(&only_nonzero_invariants, "onlyNonzeroInvariants")?;
        }
        if !include_redundant_environments.is_undefined() {
            inner.include_redundant_environments = bool_value(
                &include_redundant_environments,
                "includeRedundantEnvironments",
            )?;
        }
        if !fp_size.is_undefined() {
            inner.fp_size = u32_value(&fp_size, "fpSize")?;
        }
        if !count_simulation.is_undefined() {
            inner.count_simulation = bool_value(&count_simulation, "countSimulation")?;
        }
        if let Some(v) = optional_u32_array(&count_bounds, "countBounds")? {
            inner.count_bounds = v;
        }
        if !bits_per_feature.is_undefined() {
            inner.bits_per_feature = u32_value(&bits_per_feature, "bitsPerFeature")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=radius)]
    pub fn radius(&self) -> u32 {
        self.inner.radius
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[wasm_bindgen(getter,js_name=useBondTypes)]
    pub fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types
    }
    #[wasm_bindgen(getter,js_name=includeRingMembership)]
    pub fn include_ring_membership(&self) -> bool {
        self.inner.include_ring_membership
    }
    #[wasm_bindgen(getter,js_name=onlyNonzeroInvariants)]
    pub fn only_nonzero_invariants(&self) -> bool {
        self.inner.only_nonzero_invariants
    }
    #[wasm_bindgen(getter,js_name=includeRedundantEnvironments)]
    pub fn include_redundant_environments(&self) -> bool {
        self.inner.include_redundant_environments
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=countSimulation)]
    pub fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    #[wasm_bindgen(getter,js_name=countBounds,unchecked_return_type="number[]")]
    pub fn count_bounds(&self) -> Array {
        self.inner
            .count_bounds
            .iter()
            .copied()
            .map(JsValue::from)
            .collect()
    }
    #[wasm_bindgen(getter,js_name=bitsPerFeature)]
    pub fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> String {
        self.inner.info_string()
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
    #[wasm_bindgen(js_name=withJson)]
    pub fn with_json(&self, json: &str) -> Result<Self, JsValue> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::fingerprint_source_errors::json_error(&e).unwrap_or_else(|e| e))
    }
}
#[wasm_bindgen]
pub struct MorganInvariants {
    inner: ck::MorganInvariants,
}
#[wasm_bindgen]
impl MorganInvariants {
    pub fn connectivity() -> Self {
        Self {
            inner: ck::MorganInvariants::Connectivity,
        }
    }
    pub fn features() -> Self {
        Self {
            inner: ck::MorganInvariants::Features,
        }
    }
}
#[wasm_bindgen(
    inline_js = "export function visitMorganParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MorganParams',{cause});}} export function visitMorganInvariants(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MorganInvariants',{cause});}} export function visitMorganQuery(v,f){try{f(v);}catch(cause){throw new TypeError('invalid QueryGraph',{cause});}} export function visitMorganAtom(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MorganAtomInvariantsGenerator',{cause});}} export function visitMorganBond(v,f){try{f(v);}catch(cause){throw new TypeError('invalid MorganBondInvariantsGenerator',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitMorganParams)]
    fn visit_params(value: &JsValue, visit: &mut dyn FnMut(&MorganParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMorganInvariants)]
    fn visit_invariants(
        value: &JsValue,
        visit: &mut dyn FnMut(&MorganInvariants),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMorganQuery)]
    fn visit_query(value: &JsValue, visit: &mut dyn FnMut(&QueryGraph)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMorganAtom)]
    fn visit_atom(
        value: &JsValue,
        visit: &mut dyn FnMut(&MorganAtomInvariantsGenerator),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitMorganBond)]
    fn visit_bond(
        value: &JsValue,
        visit: &mut dyn FnMut(&MorganBondInvariantsGenerator),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct MorganFingerprintParams {
    pub(crate) inner: ck::MorganFingerprintParams,
}
#[wasm_bindgen]
impl MorganFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "MorganParams | null")] generator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms:JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_bond_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganInvariants | null")]
        invariants: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::MorganFingerprintParams::default();
        if !generator.is_null() && !generator.is_undefined() {
            visit_params(&generator, &mut |v: &MorganParams| {
                inner.generator = v.inner.clone()
            })?;
        }
        if !invariants.is_null() && !invariants.is_undefined() {
            visit_invariants(&invariants, &mut |v: &MorganInvariants| {
                inner.invariants = v.inner.clone()
            })?;
        }
        inner.from_atoms = optional_u32_array(&from_atoms, "fromAtoms")?;
        inner.ignore_atoms = optional_u32_array(&ignore_atoms, "ignoreAtoms")?;
        inner.custom_atom_invariants =
            optional_u32_array(&custom_atom_invariants, "customAtomInvariants")?;
        inner.custom_bond_invariants =
            optional_u32_array(&custom_bond_invariants, "customBondInvariants")?;
        if !conformer_id.is_undefined() {
            inner.conformer_id = i32_value(&conformer_id, "conformerId")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn generator(&self) -> MorganParams {
        MorganParams {
            inner: self.inner.generator.clone(),
        }
    }
    #[wasm_bindgen(getter)]
    pub fn invariants(&self) -> MorganInvariants {
        MorganInvariants {
            inner: self.inner.invariants.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        optional_array(self.inner.from_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=ignoreAtoms,unchecked_return_type="number[] | null")]
    pub fn ignore_atoms(&self) -> JsValue {
        optional_array(self.inner.ignore_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customAtomInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_atom_invariants(&self) -> JsValue {
        optional_array(self.inner.custom_atom_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customBondInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_bond_invariants(&self) -> JsValue {
        optional_array(self.inner.custom_bond_invariants.as_deref())
    }
}
#[wasm_bindgen]
pub struct MorganAtomInvariantsGenerator {
    pub(crate) inner: ck::MorganAtomInvariantsGenerator,
}
#[wasm_bindgen]
impl MorganAtomInvariantsGenerator {
    pub fn connectivity(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_ring_membership: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::MorganAtomInvariantsGenerator::connectivity(
                if include_ring_membership.is_undefined() {
                    true
                } else {
                    bool_value(&include_ring_membership, "includeRingMembership")?
                },
            ),
        })
    }
    pub fn features(
        #[wasm_bindgen(unchecked_optional_param_type = "QueryGraph[] | null")] patterns: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_fingerprint_values.rs::MorganAtomInvariantsGenerator::features:
        //                             .map(|q| q.inner.clone())
        // Capture detached query configuration only; no live molecule copies or matching here.
        let patterns = if patterns.is_null() || patterns.is_undefined() {
            None
        } else {
            let mut result = Vec::new();
            for v in sequence(&patterns, "patterns")?.iter() {
                visit_query(&v, &mut |q: &QueryGraph| result.push(q.inner.clone()))?;
            }
            Some(result)
        };
        Ok(Self {
            inner: ck::MorganAtomInvariantsGenerator::features(patterns),
        })
    }
    #[wasm_bindgen(js_name=atomPair)]
    pub fn atom_pair(generator: &AtomPairAtomInvariantsGenerator) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::atom_pair(generator.inner),
        }
    }
}
pub(crate) fn optional_atom(
    value: &JsValue,
) -> Result<Option<ck::MorganAtomInvariantsGenerator>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    let mut out = None;
    visit_atom(value, &mut |v: &MorganAtomInvariantsGenerator| {
        out = Some(v.inner.clone())
    })?;
    out.map(Some).ok_or_else(|| type_error("atomInvariants"))
}
pub(crate) fn optional_bond(
    value: &JsValue,
) -> Result<Option<ck::MorganBondInvariantsGenerator>, JsValue> {
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    let mut out = None;
    visit_bond(value, &mut |v: &MorganBondInvariantsGenerator| {
        out = Some(v.inner)
    })?;
    out.map(Some).ok_or_else(|| type_error("bondInvariants"))
}
#[wasm_bindgen]
pub struct MorganBondInvariantsGenerator {
    pub(crate) inner: ck::MorganBondInvariantsGenerator,
}
#[wasm_bindgen]
impl MorganBondInvariantsGenerator {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_bond_types: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(use_bond_types, include_chirality)
    }
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_bond_types: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::MorganBondInvariantsGenerator::new(
                if use_bond_types.is_undefined() {
                    true
                } else {
                    bool_value(&use_bond_types, "useBondTypes")?
                },
                if include_chirality.is_undefined() {
                    false
                } else {
                    bool_value(&include_chirality, "includeChirality")?
                },
            ),
        })
    }
    #[wasm_bindgen(js_name=useBondTypes)]
    pub fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types()
    }
    #[wasm_bindgen(js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        self.inner.include_chirality()
    }
}
#[wasm_bindgen]
pub struct MorganCallParams {
    pub(crate) inner: ck::MorganCallParams,
}
#[wasm_bindgen]
impl MorganCallParams {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms:JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_bond_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id,
        )
    }
    #[wasm_bindgen(js_name=new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms:JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_bond_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::MorganCallParams::new(
                optional_u32_array(&from_atoms, "fromAtoms")?,
                optional_u32_array(&ignore_atoms, "ignoreAtoms")?,
                optional_u32_array(&custom_atom_invariants, "customAtomInvariants")?,
                optional_u32_array(&custom_bond_invariants, "customBondInvariants")?,
                if conformer_id.is_undefined() {
                    -1
                } else {
                    i32_value(&conformer_id, "conformerId")?
                },
            ),
        })
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        optional_array(self.inner.from_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=ignoreAtoms,unchecked_return_type="number[] | null")]
    pub fn ignore_atoms(&self) -> JsValue {
        optional_array(self.inner.ignore_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customAtomInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_atom_invariants(&self) -> JsValue {
        optional_array(self.inner.custom_atom_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customBondInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_bond_invariants(&self) -> JsValue {
        optional_array(self.inner.custom_bond_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
}
