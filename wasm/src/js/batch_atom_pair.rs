//! All ten AtomPair batch methods delegate to the public batch owner.
use crate::atom_pair_parameters::AtomPairFingerprintParams;
use crate::batch_boundary::{MoleculeBatch, batch_validation_error};
use crate::batch_fingerprint_values::BatchFingerprintOutput;
use crate::batch_queries::BatchQueryParams;
use crate::fingerprint_values::{
    Fingerprint, SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32,
};
use crate::host_values::bool_value;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=fingerprintAtomPairList,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_atom_pair_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_atom_pair_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairListWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_atom_pair_list_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCountList,unchecked_return_type="(SparseCountFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_count_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_atom_pair_sparse_count_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            SparseCountFingerprint { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCountListWithParams,unchecked_return_type="(SparseCountFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_count_list_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_sparse_count_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            SparseCountFingerprint { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairCountList,unchecked_return_type="(SparseCountFingerprint32 | null)[]")]
    pub fn fingerprint_atom_pair_count_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_atom_pair_count_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            SparseCountFingerprint32 { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairCountListWithParams,unchecked_return_type="(SparseCountFingerprint32 | null)[]")]
    pub fn fingerprint_atom_pair_count_list_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_count_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            SparseCountFingerprint32 { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseBitsList,unchecked_return_type="(SparseBitFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_bits_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_atom_pair_sparse_bits_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| SparseBitFingerprint { inner }.into())
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseBitsListWithParams,unchecked_return_type="(SparseBitFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_bits_list_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_sparse_bits_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| SparseBitFingerprint { inner }.into())
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutputList,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_atom_pair_with_output_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_atom_pair_with_output_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            BatchFingerprintOutput { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutputListWithParams,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_atom_pair_with_output_list_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "boolean")] collect_additional_output: JsValue,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let collect = bool_value(&collect_additional_output, "collectAdditionalOutput")?;
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_with_output_list_with_params(&options.inner, collect, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| {
                            BatchFingerprintOutput { inner }.into()
                        })
                    })
                    .collect()
            })
    }
}
