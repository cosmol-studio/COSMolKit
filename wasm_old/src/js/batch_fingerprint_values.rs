//! Read-only canonical batch fingerprint snapshots, copied at the language boundary.
use crate::fingerprint_values::Fingerprint;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BatchFingerprintAdditionalOutput {
    pub(crate) inner: ck::BatchFingerprintAdditionalOutput,
}
#[wasm_bindgen]
impl BatchFingerprintAdditionalOutput {
    #[wasm_bindgen(js_name=atomCounts,unchecked_return_type="number[] | null")]
    pub fn atom_counts(&self) -> JsValue {
        self.inner.atom_counts().map_or(JsValue::NULL, |v| {
            v.iter()
                .copied()
                .map(JsValue::from)
                .collect::<Array>()
                .into()
        })
    }
    #[wasm_bindgen(js_name=atomToBits,unchecked_return_type="bigint[][] | null")]
    pub fn atom_to_bits(&self) -> JsValue {
        self.inner.atom_to_bits().map_or(JsValue::NULL, |v| {
            v.iter()
                .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                .collect::<Array>()
                .into()
        })
    }
    #[wasm_bindgen(js_name=bitInfoMap,unchecked_return_type="Map<bigint, [number, number][]> | null")]
    pub fn bit_info_map(&self) -> JsValue {
        self.inner.bit_info_map().map_or(JsValue::NULL, |v| {
            let map = Map::new();
            for (&key, rows) in v {
                let rows = rows
                    .iter()
                    .map(|&(a, b)| Array::of2(&a.into(), &b.into()))
                    .collect::<Array>();
                map.set(&key.into(), &rows);
            }
            map.into()
        })
    }
    #[wasm_bindgen(js_name=bitPaths,unchecked_return_type="Map<bigint, number[][]> | null")]
    pub fn bit_paths(&self) -> JsValue {
        self.inner.bit_paths().map_or(JsValue::NULL, |v| {
            let map = Map::new();
            for (&key, rows) in v {
                let rows = rows
                    .iter()
                    .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                    .collect::<Array>();
                map.set(&key.into(), &rows);
            }
            map.into()
        })
    }
    #[wasm_bindgen(js_name=atomsPerBit,unchecked_return_type="Map<bigint, number[][]> | null")]
    pub fn atoms_per_bit(&self) -> JsValue {
        self.inner.atoms_per_bit().map_or(JsValue::NULL, |v| {
            let map = Map::new();
            for (&key, rows) in v {
                let rows = rows
                    .iter()
                    .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                    .collect::<Array>();
                map.set(&key.into(), &rows);
            }
            map.into()
        })
    }
}
#[wasm_bindgen]
pub struct BatchFingerprintOutput {
    pub(crate) inner: ck::BatchFingerprintOutput,
}
#[wasm_bindgen]
impl BatchFingerprintOutput {
    pub fn fingerprint(&self) -> Fingerprint {
        Fingerprint {
            inner: self.inner.fingerprint().clone(),
        }
    }
    #[wasm_bindgen(js_name=additionalOutput)]
    pub fn additional_output(&self) -> Result<BatchFingerprintAdditionalOutput, JsValue> {
        // COSMolKit❗✔️: canonical_batch_fingerprint_values.rs:
        //         self.inner.additional_output().map(|inner| BatchFingerprintAdditionalOutput { inner: inner.clone() })
        self.inner
            .additional_output()
            .map(|inner| BatchFingerprintAdditionalOutput {
                inner: inner.clone(),
            })
            .map_err(|e| crate::fingerprint_source_errors::output_error(&e).unwrap_or_else(|e| e))
    }
}
