//! Bit-vector operations forward to the native value without temporary copies.
use crate::{fields::Field, fingerprint_errors::FingerprintError, input};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct Fingerprint(pub(crate) ck::Fingerprint);
crate::fields::value_field!(
    ck::Fingerprint,
    Fingerprint,
    FingerprintInput,
    "Fingerprint"
);
#[wasm_bindgen]
impl Fingerprint {
    #[wasm_bindgen(js_name=fromOnBits)]
    pub fn from_on_bits(
        n_bits: tsify::Ts<input::U32Input>,
        on_bits: <Vec<u32> as Field>::Input,
    ) -> Result<Fingerprint, JsValue> {
        ck::Fingerprint::from_on_bits(input::read::<u32>(n_bits)?, Vec::<u32>::input(on_bits)?)
            .map(Self)
            .map_err(FingerprintError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> u32 {
        self.0.n_bits()
    }
    #[wasm_bindgen(js_name=getBit)]
    pub fn get_bit(&self, bit: tsify::Ts<input::U32Input>) -> Result<bool, JsValue> {
        self.0
            .get_bit(input::read::<u32>(bit)?)
            .map_err(FingerprintError::from_rust)
    }
    #[wasm_bindgen(getter,js_name=onBits)]
    pub fn on_bits(&self) -> Vec<u32> {
        self.0.on_bits()
    }
    pub fn tanimoto(&self, other: FingerprintInput) -> Result<f64, JsValue> {
        Self::with_ref(&other, |other| {
            self.0.tanimoto(other).map_err(FingerprintError::from_rust)
        })
    }
    #[wasm_bindgen(js_name=toUint8Array)]
    pub fn to_uint8_array(&self) -> Result<Vec<u8>, JsValue> {
        (0..self.0.n_bits())
            .map(|bit| self.0.get_bit(bit).map(u8::from))
            .collect::<Result<Vec<_>, _>>()
            .map_err(FingerprintError::from_rust)
    }
}
