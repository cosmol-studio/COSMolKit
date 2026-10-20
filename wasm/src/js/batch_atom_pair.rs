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
    #[wasm_bindgen(js_name=fingerprintAtomPair,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_atom_pair(&self, #[wasm_bindgen(unchecked_optional_param_type = "AtomPairFingerprintParams | AtomPairFingerprintParamsOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_atom_pair()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })

            };
        }
        let configuration_0 = crate::atom_pair_parameters::AtomPairFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_atom_pair_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_atom_pair_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCount,unchecked_return_type="(SparseCountFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_count(&self, #[wasm_bindgen(unchecked_optional_param_type = "AtomPairFingerprintParams | AtomPairFingerprintParamsOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_atom_pair_sparse_count()
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

            };
        }
        let configuration_0 = crate::atom_pair_parameters::AtomPairFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_atom_pair_sparse_count_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseCountWithParams,unchecked_return_type="(SparseCountFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_count_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_sparse_count_with_params(&options.inner, p)
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
    #[wasm_bindgen(js_name=fingerprintAtomPairCount,unchecked_return_type="(SparseCountFingerprint32 | null)[]")]
    pub fn fingerprint_atom_pair_count(&self, #[wasm_bindgen(unchecked_optional_param_type = "AtomPairFingerprintParams | AtomPairFingerprintParamsOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_atom_pair_count()
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

            };
        }
        let configuration_0 = crate::atom_pair_parameters::AtomPairFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_atom_pair_count_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairCountWithParams,unchecked_return_type="(SparseCountFingerprint32 | null)[]")]
    pub fn fingerprint_atom_pair_count_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_count_with_params(&options.inner, p)
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
    #[wasm_bindgen(js_name=fingerprintAtomPairSparse,unchecked_return_type="(SparseBitFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse(&self, #[wasm_bindgen(unchecked_optional_param_type = "AtomPairFingerprintParams | AtomPairFingerprintParamsOptions")] options: JsValue, #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue) -> Result<Array, JsValue> {
        if options.is_undefined() && params.is_undefined() {
            return {
        self.inner
            .fingerprint_atom_pair_sparse()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map_or(JsValue::NULL, |inner| SparseBitFingerprint { inner }.into())
                    })
                    .collect()
            })

            };
        }
        let configuration_0 = crate::atom_pair_parameters::AtomPairFingerprintParams::from_js_configuration(&options)?;
        let configuration_1 = crate::batch_queries::BatchQueryParams::from_js_configuration(&params)?;
        self.fingerprint_atom_pair_sparse_with_params(&configuration_0, &configuration_1)
    }
    #[wasm_bindgen(js_name=fingerprintAtomPairSparseWithParams,unchecked_return_type="(SparseBitFingerprint | null)[]")]
    pub fn fingerprint_atom_pair_sparse_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_sparse_with_params(&options.inner, p)
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
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutput,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_atom_pair_with_output(&self,
        #[wasm_bindgen(unchecked_optional_param_type = "AtomPairFingerprintParams | AtomPairFingerprintParamsOptions")] options: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "BatchQueryParams | BatchQueryParamsOptions")] params: JsValue,
    ) -> Result<Array, JsValue> {
        if !options.is_undefined() || !params.is_undefined() {
            let options = AtomPairFingerprintParams::from_js_configuration(&options)?;
            let params = BatchQueryParams::from_js_configuration(&params)?;
            return self.fingerprint_atom_pair_with_output_with_params(&options, true.into(), &params);
        }
        self.inner
            .fingerprint_atom_pair_with_output()
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
    #[wasm_bindgen(js_name=fingerprintAtomPairWithOutputWithParams,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_atom_pair_with_output_with_params(
        &self,
        options: &AtomPairFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "boolean")] collect_additional_output: JsValue,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let collect = bool_value(&collect_additional_output, "collectAdditionalOutput")?;
        params
            .execute(|p| {
                self.inner
                    .fingerprint_atom_pair_with_output_with_params(&options.inner, collect, p)
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
