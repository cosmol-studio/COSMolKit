//! Complete immutable torsion configurations, with canonical defaults and JSON.
use crate::atom_pair_parameters::{
    AtomPairAtomInvariantsGenerator, optional_array, optional_u32_array,
};
use crate::host_values::{bool_value, i32_value, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct TopologicalTorsionParams {
    pub(crate) inner: ck::TopologicalTorsionParams,
}
#[wasm_bindgen]
impl TopologicalTorsionParams {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] torsion_atom_count: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] only_shortest_paths: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] count_simulation: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bits_per_feature: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        count_bounds: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(
            torsion_atom_count,
            only_shortest_paths,
            include_chirality,
            count_simulation,
            fp_size,
            bits_per_feature,
            count_bounds,
        )
    }
    fn new(
        torsion_atom_count: JsValue,
        only_shortest_paths: JsValue,
        include_chirality: JsValue,
        count_simulation: JsValue,
        fp_size: JsValue,
        bits_per_feature: JsValue,
        count_bounds: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::TopologicalTorsionParams {
        let mut inner = ck::TopologicalTorsionParams::default();
        if !torsion_atom_count.is_undefined() {
            inner.torsion_atom_count = u32_value(&torsion_atom_count, "torsionAtomCount")?;
        }
        if !only_shortest_paths.is_undefined() {
            inner.only_shortest_paths = bool_value(&only_shortest_paths, "onlyShortestPaths")?;
        }
        if !include_chirality.is_undefined() {
            inner.include_chirality = bool_value(&include_chirality, "includeChirality")?;
        }
        if !count_simulation.is_undefined() {
            inner.count_simulation = bool_value(&count_simulation, "countSimulation")?;
        }
        if !fp_size.is_undefined() {
            inner.fp_size = u32_value(&fp_size, "fpSize")?;
        }
        if !bits_per_feature.is_undefined() {
            inner.bits_per_feature = u32_value(&bits_per_feature, "bitsPerFeature")?;
        }
        if let Some(v) = optional_u32_array(&count_bounds, "countBounds")? {
            inner.count_bounds = v;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=torsionAtomCount)]
    pub fn torsion_atom_count(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.torsion_atom_count
        self.inner.torsion_atom_count
    }
    #[wasm_bindgen(getter,js_name=onlyShortestPaths)]
    pub fn only_shortest_paths(&self) -> bool {
        // COSMolKit❗✔️: self.inner.only_shortest_paths
        self.inner.only_shortest_paths
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        // COSMolKit❗✔️: self.inner.include_chirality
        self.inner.include_chirality
    }
    #[wasm_bindgen(getter,js_name=countSimulation)]
    pub fn count_simulation(&self) -> bool {
        // COSMolKit❗✔️: self.inner.count_simulation
        self.inner.count_simulation
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.fp_size
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=bitsPerFeature)]
    pub fn bits_per_feature(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.bits_per_feature
        self.inner.bits_per_feature
    }
    #[wasm_bindgen(getter,js_name=countBounds,unchecked_return_type="number[]")]
    pub fn count_bounds(&self) -> Array {
        // COSMolKit❗✔️: self.inner.count_bounds
        self.inner
            .count_bounds
            .iter()
            .copied()
            .map(JsValue::from)
            .collect()
    }
}
#[wasm_bindgen(
    inline_js = "export function visitTorsionParams(value,visit){try { visit(value); } catch(cause) { throw new TypeError('invalid AtomPair parameter', {cause}); }} export function visitTorsionInvariants(value,visit){try { visit(value); } catch(cause) { throw new TypeError('invalid AtomPair parameter', {cause}); }}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitTorsionParams)]
    fn visit_params(
        value: &JsValue,
        visit: &mut dyn FnMut(&TopologicalTorsionParams),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitTorsionInvariants)]
    fn visit_invariants(
        value: &JsValue,
        visit: &mut dyn FnMut(&AtomPairAtomInvariantsGenerator),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct TopologicalTorsionFingerprintParams {
    pub(crate) inner: ck::TopologicalTorsionFingerprintParams,
}
#[wasm_bindgen]
impl TopologicalTorsionFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionParams | null")]
        generator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_bond_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AtomPairAtomInvariantsGenerator | null")]
        atom_invariants_generator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        use_legacy_stereo_perception: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::TopologicalTorsionFingerprintParams::default();
        if !generator.is_null() && !generator.is_undefined() {
            let mut copied = None;
            visit_params(&generator, &mut |value: &TopologicalTorsionParams| {
                copied = Some(value.inner.clone())
            })?;
            inner.generator = copied.ok_or_else(|| type_error("generator"))?;
        }
        if !atom_invariants_generator.is_null() && !atom_invariants_generator.is_undefined() {
            let mut copied = None;
            visit_invariants(
                &atom_invariants_generator,
                &mut |value: &AtomPairAtomInvariantsGenerator| copied = Some(value.inner),
            )?;
            inner.atom_invariants_generator = copied;
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
        if !use_legacy_stereo_perception.is_undefined() {
            inner.use_legacy_stereo_perception =
                bool_value(&use_legacy_stereo_perception, "useLegacyStereoPerception")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn generator(&self) -> TopologicalTorsionParams {
        TopologicalTorsionParams {
            inner: self.inner.generator.clone(),
        }
    }
    #[wasm_bindgen(getter,js_name=atomInvariantsGenerator,unchecked_return_type="AtomPairAtomInvariantsGenerator | null")]
    pub fn atom_invariants_generator(&self) -> JsValue {
        self.inner
            .atom_invariants_generator
            .map_or(JsValue::NULL, |inner| {
                AtomPairAtomInvariantsGenerator { inner }.into()
            })
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[wasm_bindgen(getter,js_name=useLegacyStereoPerception)]
    pub fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
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
pub struct TopologicalTorsionCallParams {
    pub(crate) inner: ck::TopologicalTorsionCallParams,
}
#[wasm_bindgen]
impl TopologicalTorsionCallParams {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")]from_atoms:JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_bond_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        use_legacy_stereo_perception: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id,
            use_legacy_stereo_perception,
        )
    }
    fn new(
        from_atoms: JsValue,
        ignore_atoms: JsValue,
        custom_atom_invariants: JsValue,
        custom_bond_invariants: JsValue,
        conformer_id: JsValue,
        use_legacy_stereo_perception: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::TopologicalTorsionCallParams {
        let mut inner = ck::TopologicalTorsionCallParams::default();
        inner.from_atoms = optional_u32_array(&from_atoms, "fromAtoms")?;
        inner.ignore_atoms = optional_u32_array(&ignore_atoms, "ignoreAtoms")?;
        inner.custom_atom_invariants =
            optional_u32_array(&custom_atom_invariants, "customAtomInvariants")?;
        inner.custom_bond_invariants =
            optional_u32_array(&custom_bond_invariants, "customBondInvariants")?;
        if !conformer_id.is_undefined() {
            inner.conformer_id = i32_value(&conformer_id, "conformerId")?;
        }
        if !use_legacy_stereo_perception.is_undefined() {
            inner.use_legacy_stereo_perception =
                bool_value(&use_legacy_stereo_perception, "useLegacyStereoPerception")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.from_atoms
        optional_array(self.inner.from_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=ignoreAtoms,unchecked_return_type="number[] | null")]
    pub fn ignore_atoms(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.ignore_atoms
        optional_array(self.inner.ignore_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customAtomInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_atom_invariants(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.custom_atom_invariants
        optional_array(self.inner.custom_atom_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customBondInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_bond_invariants(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.custom_bond_invariants
        optional_array(self.inner.custom_bond_invariants.as_deref())
    }
    #[wasm_bindgen(getter,js_name=conformerId)]
    pub fn conformer_id(&self) -> i32 {
        // COSMolKit❗✔️: self.inner.conformer_id
        self.inner.conformer_id
    }
    #[wasm_bindgen(getter,js_name=useLegacyStereoPerception)]
    pub fn use_legacy_stereo_perception(&self) -> bool {
        // COSMolKit❗✔️: self.inner.use_legacy_stereo_perception
        self.inner.use_legacy_stereo_perception
    }
}
#[wasm_bindgen]
pub struct LegacyTopologicalTorsionParams {
    pub(crate) inner: ck::LegacyTopologicalTorsionParams,
}
#[wasm_bindgen]
impl LegacyTopologicalTorsionParams {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] torsion_atom_count: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bits_per_entry: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")]from_atoms:JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(
            torsion_atom_count,
            include_chirality,
            fp_size,
            bits_per_entry,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
        )
    }
    #[wasm_bindgen(js_name=new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] torsion_atom_count: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bits_per_entry: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        ignore_atoms: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        custom_atom_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::LegacyTopologicalTorsionParams {
        let mut inner = ck::LegacyTopologicalTorsionParams::default();
        if !torsion_atom_count.is_undefined() {
            inner.torsion_atom_count = u32_value(&torsion_atom_count, "torsionAtomCount")?;
        }
        if !include_chirality.is_undefined() {
            inner.include_chirality = bool_value(&include_chirality, "includeChirality")?;
        }
        if !fp_size.is_undefined() {
            inner.fp_size = u32_value(&fp_size, "fpSize")?;
        }
        if !bits_per_entry.is_undefined() {
            inner.bits_per_entry = u32_value(&bits_per_entry, "bitsPerEntry")?;
        }
        inner.from_atoms = optional_u32_array(&from_atoms, "fromAtoms")?;
        inner.ignore_atoms = optional_u32_array(&ignore_atoms, "ignoreAtoms")?;
        inner.custom_atom_invariants =
            optional_u32_array(&custom_atom_invariants, "customAtomInvariants")?;
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=torsionAtomCount)]
    pub fn torsion_atom_count(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.torsion_atom_count
        self.inner.torsion_atom_count
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        // COSMolKit❗✔️: self.inner.include_chirality
        self.inner.include_chirality
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.fp_size
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=bitsPerEntry)]
    pub fn bits_per_entry(&self) -> u32 {
        // COSMolKit❗✔️: self.inner.bits_per_entry
        self.inner.bits_per_entry
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.from_atoms
        optional_array(self.inner.from_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=ignoreAtoms,unchecked_return_type="number[] | null")]
    pub fn ignore_atoms(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.ignore_atoms
        optional_array(self.inner.ignore_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=customAtomInvariants,unchecked_return_type="number[] | null")]
    pub fn custom_atom_invariants(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.custom_atom_invariants
        optional_array(self.inner.custom_atom_invariants.as_deref())
    }
}

#[wasm_bindgen]
impl TopologicalTorsionParams {
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
