//! Full Layered and Pattern configuration/result transport.
use crate::atom_pair_parameters::{optional_array, optional_u32_array};
use crate::fingerprint_values::Fingerprint;
use crate::host_values::{bool_value, u32_value, usize_value};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct LayeredFingerprintLayers {
    inner: ck::LayeredFingerprintLayers,
}
#[wasm_bindgen]
impl LayeredFingerprintLayers {
    #[wasm_bindgen(js_name=fromBitsRetain)]
    pub fn from_bits_retain(
        #[wasm_bindgen(unchecked_param_type = "number")] bits: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::LayeredFingerprintLayers::from_bits_retain(u32_value(&bits, "bits")?),
        })
    }
    pub fn bits(&self) -> u32 {
        self.inner.bits()
    }
}
#[wasm_bindgen(
    inline_js = "export function visitLayeredMask(value,visit){try{visit(value);}catch(cause){throw new TypeError('invalid Fingerprint',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitLayeredMask)]
    fn visit_mask(value: &JsValue, visit: &mut dyn FnMut(&Fingerprint)) -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct LayeredFingerprintParams {
    pub(crate) inner: ck::LayeredFingerprintParams,
}
#[wasm_bindgen]
impl LayeredFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] layers: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] min_path: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_path: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] fp_size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number[] | Uint32Array | null")]
        atom_counts: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "Fingerprint | null")]
        set_only_bits: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] branched_paths: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type="number[] | Uint32Array | null")] from_atoms:JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: canonical_layered.rs::LayeredFingerprintParams::new:
        //                 layers: ck::LayeredFingerprintLayers::from_bits_retain(layers),
        //                 set_only_bits: set_only_bits.map(|value| value.inner.clone()),
        // Preserve source flag bits and detached mask ownership; no paths enumerated here.
        let mut inner = ck::LayeredFingerprintParams::default();
        if !layers.is_undefined() {
            inner.layers =
                ck::LayeredFingerprintLayers::from_bits_retain(u32_value(&layers, "layers")?);
        }
        if !min_path.is_undefined() {
            inner.min_path = u32_value(&min_path, "minPath")?;
        }
        if !max_path.is_undefined() {
            inner.max_path = u32_value(&max_path, "maxPath")?;
        }
        if !fp_size.is_undefined() {
            inner.fp_size = u32_value(&fp_size, "fpSize")?;
        }
        inner.atom_counts = optional_u32_array(&atom_counts, "atomCounts")?;
        inner.from_atoms = optional_u32_array(&from_atoms, "fromAtoms")?;
        if !set_only_bits.is_null() && !set_only_bits.is_undefined() {
            visit_mask(&set_only_bits, &mut |value: &Fingerprint| {
                inner.set_only_bits = Some(value.inner.clone())
            })?;
        }
        if !branched_paths.is_undefined() {
            inner.branched_paths = bool_value(&branched_paths, "branchedPaths")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter)]
    pub fn layers(&self) -> u32 {
        self.inner.layers.bits()
    }
    #[wasm_bindgen(getter,js_name=minPath)]
    pub fn min_path(&self) -> u32 {
        self.inner.min_path
    }
    #[wasm_bindgen(getter,js_name=maxPath)]
    pub fn max_path(&self) -> u32 {
        self.inner.max_path
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[wasm_bindgen(getter,js_name=atomCounts,unchecked_return_type="number[] | null")]
    pub fn atom_counts(&self) -> JsValue {
        optional_array(self.inner.atom_counts.as_deref())
    }
    #[wasm_bindgen(getter,js_name=fromAtoms,unchecked_return_type="number[] | null")]
    pub fn from_atoms(&self) -> JsValue {
        optional_array(self.inner.from_atoms.as_deref())
    }
    #[wasm_bindgen(getter,js_name=setOnlyBits,unchecked_return_type="Fingerprint | null")]
    pub fn set_only_bits(&self) -> JsValue {
        self.inner
            .set_only_bits
            .as_ref()
            .map_or(JsValue::NULL, |inner| {
                Fingerprint {
                    inner: inner.clone(),
                }
                .into()
            })
    }
    #[wasm_bindgen(getter,js_name=branchedPaths)]
    pub fn branched_paths(&self) -> bool {
        self.inner.branched_paths
    }
}
#[wasm_bindgen]
pub struct PatternFingerprintParams {
    pub(crate) inner: ck::PatternFingerprintParams,
}
#[wasm_bindgen]
impl PatternFingerprintParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] n_bits: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] tautomeric: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::PatternFingerprintParams::default();
        if !n_bits.is_undefined() {
            inner.n_bits = usize_value(&n_bits, "nBits")?;
        }
        if !tautomeric.is_undefined() {
            inner.tautomeric = bool_value(&tautomeric, "tautomeric")?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> usize {
        self.inner.n_bits
    }
    #[wasm_bindgen(getter)]
    pub fn tautomeric(&self) -> bool {
        self.inner.tautomeric
    }
}
#[wasm_bindgen]
pub struct LayeredFingerprintResult {
    pub(crate) inner: ck::LayeredFingerprintResult,
}
#[wasm_bindgen]
impl LayeredFingerprintResult {
    pub fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[wasm_bindgen(js_name=atomCounts,unchecked_return_type="number[] | null")]
    pub fn atom_counts(&self) -> JsValue {
        optional_array(self.inner.atom_counts())
    }
}
