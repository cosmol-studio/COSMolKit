//! Read-only native fingerprint results and immutable provenance values.
use crate::{errors::enum_error, fingerprint_output::*, fingerprint_values::Fingerprint};
use cosmolkit as ck;
use wasm_bindgen::{JsCast, prelude::*};
enum_error!(BatchFingerprintOutputError,BatchFingerprintOutputErrorKind,ck::BatchFingerprintOutputError,{
    MissingAdditionalOutput=>ck::BatchFingerprintOutputError::MissingAdditionalOutput{..}
},{
    #[wasm_bindgen(getter,js_name=fingerprintKind)]
    pub fn fingerprint_kind(&self)->String { match self.inner.as_ref(){ck::BatchFingerprintOutputError::MissingAdditionalOutput{fingerprint_kind}=> (*fingerprint_kind).into()} }
});
#[wasm_bindgen]
pub struct BatchFingerprintOutput(pub(crate) ck::BatchFingerprintOutput);
crate::fields::value_field!(
    ck::BatchFingerprintOutput,
    BatchFingerprintOutput,
    BatchFingerprintOutputInput,
    "BatchFingerprintOutput"
);
crate::fields::optional_value_array!(
    ck::BatchFingerprintOutput,
    BatchFingerprintOutput,
    BatchFingerprintOutputs,
    "ReadonlyArray<BatchFingerprintOutput | undefined>"
);
#[wasm_bindgen]
impl BatchFingerprintOutput {
    #[wasm_bindgen(getter)]
    pub fn fingerprint(&self) -> Fingerprint {
        Fingerprint(self.0.fingerprint().clone())
    }
    #[wasm_bindgen(getter,js_name=additionalOutput)]
    pub fn additional_output(&self) -> Result<BatchFingerprintAdditionalOutput, JsValue> {
        self.0
            .additional_output()
            .map(|v| BatchFingerprintAdditionalOutput(v.clone()))
            .map_err(BatchFingerprintOutputError::from_rust)
    }
}
#[wasm_bindgen]
pub struct BatchFingerprintAdditionalOutput(ck::BatchFingerprintAdditionalOutput);
#[wasm_bindgen]
impl BatchFingerprintAdditionalOutput {
    #[wasm_bindgen(getter,js_name=atomCounts)]
    pub fn atom_counts(&self) -> Option<AtomCounts> {
        self.0
            .atom_counts()
            .map(|v| array(v.iter().copied().map(JsValue::from)).unchecked_into())
    }
    #[wasm_bindgen(getter,js_name=atomToBits)]
    pub fn atom_to_bits(&self) -> Option<AtomToBits> {
        self.0.atom_to_bits().map(|v| {
            array(
                v.iter()
                    .map(|r| array(r.iter().copied().map(JsValue::from))),
            )
            .unchecked_into()
        })
    }
    #[wasm_bindgen(getter,js_name=bitInfoMap)]
    pub fn bit_info_map(&self) -> Result<Option<BitInfoMap>, JsValue> {
        self.0.bit_info_map().map(info_map).transpose()
    }
    #[wasm_bindgen(getter,js_name=bitPaths)]
    pub fn bit_paths(&self) -> Result<Option<BitPaths>, JsValue> {
        self.0.bit_paths().map(paths).transpose()
    }
    #[wasm_bindgen(getter,js_name=atomsPerBit)]
    pub fn atoms_per_bit(&self) -> Result<Option<BitPaths>, JsValue> {
        self.0.atoms_per_bit().map(paths).transpose()
    }
}
