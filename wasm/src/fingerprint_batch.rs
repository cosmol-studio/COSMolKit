//! Dense collection transport; shape validation is the native owner's query.
use crate::{
    errors::{enum_error, exception},
    fingerprint_values::Fingerprint,
    input,
};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
crate::fields::optional_value_array!(
    ck::Fingerprint,
    Fingerprint,
    FingerprintArray,
    "ReadonlyArray<Fingerprint | null | undefined>"
);
enum_error!(FingerprintBatchError,FingerprintBatchErrorKind,ck::FingerprintBatchError,{
    MissingRow=>ck::FingerprintBatchError::MissingRow{..},WidthMismatch=>ck::FingerprintBatchError::WidthMismatch{..}
},{
    #[wasm_bindgen(getter)] pub fn index(&self)->usize{match *self.inner{ck::FingerprintBatchError::MissingRow{index}|ck::FingerprintBatchError::WidthMismatch{index,..}=>index}}
    #[wasm_bindgen(getter)] pub fn actual(&self)->Option<u32>{match *self.inner{ck::FingerprintBatchError::WidthMismatch{actual,..}=>Some(actual),_=>None}}
    #[wasm_bindgen(getter)] pub fn expected(&self)->Option<u32>{match *self.inner{ck::FingerprintBatchError::WidthMismatch{expected,..}=>Some(expected),_=>None}}
});
#[wasm_bindgen]
pub struct FingerprintBatch(pub(crate) ck::FingerprintBatch);
#[wasm_bindgen]
impl FingerprintBatch {
    #[wasm_bindgen(constructor)]
    pub fn new(values: FingerprintArray) -> Result<Self, JsValue> {
        Ok(Self(ck::FingerprintBatch::new(values.native_values()?)))
    }
    #[wasm_bindgen(getter)]
    pub fn len(&self) -> usize {
        self.0.len()
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> u32 {
        self.0.n_bits()
    }
    pub fn get(&self, index: tsify::Ts<input::IndexInput>) -> Result<Option<Fingerprint>, JsValue> {
        self.0
            .values()
            .get(input::read::<usize>(index)?)
            .map(|value| value.clone().map(Fingerprint))
            .ok_or_else(|| js_sys::RangeError::new("fingerprint index is out of range").into())
    }
    #[wasm_bindgen(js_name=toUint8Array)]
    pub fn to_uint8_array(&self) -> Result<Vec<u8>, JsValue> {
        self.0
            .validate_matrix()
            .map_err(FingerprintBatchError::from_rust)?;
        self.0
            .values()
            .iter()
            .flat_map(|value| value.as_ref().into_iter())
            .flat_map(|value| (0..value.n_bits()).map(|bit| value.get_bit(bit).map(u8::from)))
            .collect::<Result<Vec<_>, _>>()
            .map_err(crate::fingerprint_errors::FingerprintError::from_rust)
    }
}
