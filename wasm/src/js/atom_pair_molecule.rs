//! Scalar AtomPair calls borrow canonical parameters and mutable additional output.
use crate::Molecule;
use crate::atom_pair_parameters::AtomPairFingerprintParams;
use crate::fingerprint_additional::with_output;
use crate::fingerprint_source_errors::atom_pair_error;
use crate::fingerprint_values::{
    Fingerprint, SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintAtomPair)]
    pub fn fingerprint_atom_pair(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_fingerprint(
        self.inner
            .fingerprint_atom_pair()
            .map(|inner| Fingerprint { inner })
            .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairWithParams)]
    pub fn fingerprint_atom_pair_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_atom_pair_with_params(&params.inner, output)
                .map(|inner| Fingerprint { inner })
                .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
        })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparse)]
    pub fn fingerprint_atom_pair_sparse(&self) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_sparse_fingerprint(
        self.inner
            .fingerprint_atom_pair_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseWithParams)]
    pub fn fingerprint_atom_pair_sparse_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_sparse_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_atom_pair_sparse_with_params(&params.inner, output)
                .map(|inner| SparseBitFingerprint { inner })
                .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
        })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairCount)]
    pub fn fingerprint_atom_pair_count(&self) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .atom_pair_count_fingerprint(
        self.inner
            .fingerprint_atom_pair_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairCountWithParams)]
    pub fn fingerprint_atom_pair_count_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .atom_pair_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_atom_pair_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint32 { inner })
                .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
        })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCount)]
    pub fn fingerprint_atom_pair_sparse_count(&self) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_sparse_count_fingerprint(
        self.inner
            .fingerprint_atom_pair_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCountWithParams)]
    pub fn fingerprint_atom_pair_sparse_count_with_params(
        &self,
        params: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .atom_pair_sparse_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_atom_pair_sparse_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint { inner })
                .map_err(|e| atom_pair_error(&e).unwrap_or_else(|e| e))
        })
    }
}
