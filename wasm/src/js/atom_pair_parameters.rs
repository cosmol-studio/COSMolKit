//! Complete immutable AtomPair configuration transport through cosmolkit.
use crate::fingerprint_values::u32_array;
use crate::host_values::{bool_value, i32_value, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
pub(crate) fn optional_u32_array(value: &JsValue, name: &str) -> Result<Option<Vec<u32>>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        u32_array(value, name).map(Some)
    }
}
pub(crate) fn optional_array(value: Option<&[u32]>) -> JsValue {
    value.map_or(JsValue::NULL, |v| {
        v.iter()
            .copied()
            .map(JsValue::from)
            .collect::<Array>()
            .into()
    })
}
#[wasm_bindgen]
pub struct AtomPairParams {
    pub(crate) inner: ck::AtomPairParams,
}
#[wasm_bindgen]
impl AtomPairParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] min_distance: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_distance: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_2d: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] count_simulation: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] bits_per_feature: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        count_bounds: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_fingerprint_values.rs::AtomPairParams::new:
        //                 count_bounds: count_bounds
        //                     .unwrap_or_else(|| ck::AtomPairParams::default().count_bounds),
        // O(fields + supplied arrays); no generator validation or chemistry duplicated.
        let mut inner = ck::AtomPairParams::default();
        if !min_distance.is_undefined() {
            inner.min_distance = u32_value(&min_distance, "minDistance")?;
        }
        if !max_distance.is_undefined() {
            inner.max_distance = u32_value(&max_distance, "maxDistance")?;
        }
        if !include_chirality.is_undefined() {
            inner.include_chirality = bool_value(&include_chirality, "includeChirality")?;
        }
        if !use_2d.is_undefined() {
            inner.use_2d = bool_value(&use_2d, "use2d")?;
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
        if let Some(bounds) = optional_u32_array(&count_bounds, "countBounds")? {
            inner.count_bounds = bounds;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=minDistance)]
    pub fn min_distance(&self) -> u32 {
        self.inner.min_distance
    }
    #[wasm_bindgen(getter,js_name=maxDistance)]
    pub fn max_distance(&self) -> u32 {
        self.inner.max_distance
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[wasm_bindgen(getter,js_name=use2d)]
    pub fn use_2d(&self) -> bool {
        self.inner.use_2d
    }
    #[wasm_bindgen(getter,js_name=countSimulation)]
    pub fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=bitsPerFeature)]
    pub fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
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
pub struct AtomPairAtomInvariantsGenerator {
    pub(crate) inner: ck::AtomPairAtomInvariantsGenerator,
}
#[wasm_bindgen]
impl AtomPairAtomInvariantsGenerator {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] include_chirality: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        topological_torsion_correction: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::AtomPairAtomInvariantsGenerator {
                include_chirality: if include_chirality.is_undefined() {
                    false
                } else {
                    bool_value(&include_chirality, "includeChirality")?
                },
                topological_torsion_correction: if topological_torsion_correction.is_undefined() {
                    false
                } else {
                    bool_value(
                        &topological_torsion_correction,
                        "topologicalTorsionCorrection",
                    )?
                },
            },
        })
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[wasm_bindgen(getter,js_name=topologicalTorsionCorrection)]
    pub fn topological_torsion_correction(&self) -> bool {
        self.inner.topological_torsion_correction
    }
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> String {
        self.inner.info_string()
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> String {
        self.inner.to_json()
    }
}
#[wasm_bindgen(
    inline_js = "export function visitAtomPairParams(value,visit){try { visit(value); } catch(cause) { throw new TypeError('invalid AtomPair parameter', {cause}); }} export function visitAtomPairInvariants(value,visit){try { visit(value); } catch(cause) { throw new TypeError('invalid AtomPair parameter', {cause}); }}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitAtomPairParams)]
    fn visit_params(value: &JsValue, visit: &mut dyn FnMut(&AtomPairParams))
    -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitAtomPairInvariants)]
    fn visit_invariants(
        value: &JsValue,
        visit: &mut dyn FnMut(&AtomPairAtomInvariantsGenerator),
    ) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct AtomPairFingerprintParams {
    pub(crate) inner: ck::AtomPairFingerprintParams,
}
#[wasm_bindgen]
impl AtomPairFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "AtomPairParams | null")] generator: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms:JsValue,
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
        let mut inner = ck::AtomPairFingerprintParams::default();
        if !generator.is_null() && !generator.is_undefined() {
            let mut copied = None;
            visit_params(&generator, &mut |value: &AtomPairParams| {
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
    pub fn generator(&self) -> AtomPairParams {
        AtomPairParams {
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
