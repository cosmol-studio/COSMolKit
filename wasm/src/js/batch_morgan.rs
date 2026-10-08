//! Complete six Morgan batch adapters through native public forwarding.
use crate::batch_boundary::{MoleculeBatch, batch_validation_error};
use crate::batch_fingerprint_values::BatchFingerprintOutput;
use crate::batch_queries::BatchQueryParams;
use crate::fingerprint_values::Fingerprint;
use crate::host_values::bool_value;
use crate::morgan_parameters::{
    MorganCallParams, MorganFingerprintParams, MorganParams, optional_atom, optional_bond,
};
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=fingerprintMorganList,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_morgan_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_morgan_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| v.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintMorganListWithParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_morgan_list_with_params(
        &self,
        options: &MorganFingerprintParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| {
                self.inner
                    .fingerprint_morgan_list_with_params(&options.inner, p)
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| v.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintMorganWithOutputList,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_morgan_with_output_list(&self) -> Result<Array, JsValue> {
        self.inner
            .fingerprint_morgan_with_output_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| {
                        v.map_or(JsValue::NULL, |inner| {
                            BatchFingerprintOutput { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintMorganWithOutputListWithParams,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_morgan_with_output_list_with_params(
        &self,
        options: &MorganFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "boolean")] collect_additional_output: JsValue,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let collect = bool_value(&collect_additional_output, "collectAdditionalOutput")?;
        params
            .execute(|p| {
                self.inner.fingerprint_morgan_with_output_list_with_params(
                    &options.inner,
                    collect,
                    p,
                )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| {
                        v.map_or(JsValue::NULL, |inner| {
                            BatchFingerprintOutput { inner }.into()
                        })
                    })
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintMorganListWithGeneratorParams,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprint_morgan_list_with_generator_params(
        &self,
        options: &MorganParams,
        #[wasm_bindgen(unchecked_param_type = "MorganAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_param_type = "MorganBondInvariantsGenerator | null")]
        bond_invariants: JsValue,
        call: &MorganCallParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let atom = optional_atom(&atom_invariants)?;
        let bond = optional_bond(&bond_invariants)?;
        params
            .execute(|p| {
                self.inner.fingerprint_morgan_list_with_generator_params(
                    &options.inner,
                    atom.as_ref(),
                    bond.as_ref(),
                    &call.inner,
                    p,
                )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| v.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                    .collect()
            })
    }
    #[wasm_bindgen(js_name=fingerprintMorganWithOutputListWithGeneratorParams,unchecked_return_type="(BatchFingerprintOutput | null)[]")]
    pub fn fingerprint_morgan_with_output_list_with_generator_params(
        &self,
        options: &MorganParams,
        #[wasm_bindgen(unchecked_param_type = "MorganAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_param_type = "MorganBondInvariantsGenerator | null")]
        bond_invariants: JsValue,
        call: &MorganCallParams,
        #[wasm_bindgen(unchecked_param_type = "boolean")] collect_additional_output: JsValue,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let atom = optional_atom(&atom_invariants)?;
        let bond = optional_bond(&bond_invariants)?;
        let collect = bool_value(&collect_additional_output, "collectAdditionalOutput")?;
        params
            .execute(|p| {
                self.inner
                    .fingerprint_morgan_with_output_list_with_generator_params(
                        &options.inner,
                        atom.as_ref(),
                        bond.as_ref(),
                        &call.inner,
                        collect,
                        p,
                    )
            })
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| {
                        v.map_or(JsValue::NULL, |inner| {
                            BatchFingerprintOutput { inner }.into()
                        })
                    })
                    .collect()
            })
    }
}
