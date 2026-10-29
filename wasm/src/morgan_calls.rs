//! Molecular Morgan calls borrow both configuration and optional output state.
use crate::{
    Molecule,
    fingerprint_output::{FingerprintAdditionalOutput, FingerprintAdditionalOutputInput},
    morgan_configuration::{MorganFingerprintParams, MorganFingerprintParamsInput},
    morgan_errors::MorganReadError,
};
use wasm_bindgen::prelude::*;
macro_rules! morgan_call {
    ($method:ident,$js:ident,$native:ident,$output:path) => {
        #[wasm_bindgen]
        impl Molecule {
            #[wasm_bindgen(js_name=$js)]
            pub fn $method(
                &self,
                params: Option<MorganFingerprintParamsInput>,
                additional_output: Option<FingerprintAdditionalOutputInput>,
            ) -> Result<$output, JsValue> {
                MorganFingerprintParams::with_input(params.as_ref(), |params| {
                    FingerprintAdditionalOutput::with_optional(
                        additional_output.as_ref(),
                        |output| {
                            self.0
                                .$native(params, output)
                                .map($output)
                                .map_err(MorganReadError::from_rust)
                        },
                    )
                })
            }
        }
        paste::paste! {
            #[wasm_bindgen]
            impl Molecule {
                #[wasm_bindgen(js_name=[<$js WithParams>])]
                pub fn $native(&self, params: MorganFingerprintParamsInput, additional_output: Option<FingerprintAdditionalOutputInput>) -> Result<$output, JsValue> {
                    self.$method(Some(params), additional_output)
                }
            }
        }
    };
}
morgan_call!(
    fingerprint_morgan,
    fingerprintMorgan,
    fingerprint_morgan_with_params,
    crate::fingerprint_values::Fingerprint
);
morgan_call!(
    fingerprint_morgan_sparse,
    fingerprintMorganSparse,
    fingerprint_morgan_sparse_with_params,
    crate::sparse_fingerprints::SparseBitFingerprint
);
morgan_call!(
    fingerprint_morgan_count,
    fingerprintMorganCount,
    fingerprint_morgan_count_with_params,
    crate::sparse_fingerprints::SparseCountFingerprint32
);
morgan_call!(
    fingerprint_morgan_sparse_count,
    fingerprintMorganSparseCount,
    fingerprint_morgan_sparse_count_with_params,
    crate::sparse_fingerprints::SparseCountFingerprint
);
