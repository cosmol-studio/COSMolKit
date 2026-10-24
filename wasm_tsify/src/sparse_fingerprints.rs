//! Sparse fingerprint values; arithmetic remains entirely in the Rust facade.
use crate::{fields::Field, fingerprint_errors::FingerprintError, input};
use cosmolkit as ck;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub struct SparseBitFingerprint(pub(crate) ck::SparseBitFingerprint);
crate::fields::optional_value_array!(
    ck::SparseBitFingerprint,
    SparseBitFingerprint,
    SparseBitFingerprintArray,
    "ReadonlyArray<SparseBitFingerprint | undefined>"
);
crate::fields::optional_value_array!(
    ck::SparseCountFingerprint,
    SparseCountFingerprint,
    SparseCountFingerprintArray,
    "ReadonlyArray<SparseCountFingerprint | undefined>"
);
crate::fields::optional_value_array!(
    ck::SparseCountFingerprint32,
    SparseCountFingerprint32,
    SparseCountFingerprint32Array,
    "ReadonlyArray<SparseCountFingerprint32 | undefined>"
);
crate::fields::value_field!(
    ck::SparseBitFingerprint,
    SparseBitFingerprint,
    SparseBitFingerprintInput,
    "SparseBitFingerprint"
);
#[wasm_bindgen]
impl SparseBitFingerprint {
    #[wasm_bindgen(constructor)]
    pub fn new(n_bits: tsify::Ts<input::U32Input>) -> Result<Self, JsValue> {
        Ok(Self(ck::SparseBitFingerprint::new(input::read::<u32>(
            n_bits,
        )?)))
    }
    #[wasm_bindgen(getter,js_name=nBits)]
    pub fn n_bits(&self) -> u32 {
        self.0.n_bits()
    }
    #[wasm_bindgen(getter,js_name=onBits)]
    pub fn on_bits(&self) -> Result<<Vec<i32> as Field>::Output, JsValue> {
        <Vec<i32> as Field>::get(crate::storage::NativeField::owned(self.0.on_bits()))
    }
}

// Both native count widths provide the same operations. Only the index ABI
// differs: u64 stays bigint; u32 is a checked number.
macro_rules! sparse_counts {
    ($name:ident,$input:ident,$elements:ident,$native:ty,$index:ty,$wire:ty,$ts:literal,$class:literal) => {
        #[wasm_bindgen]
        pub struct $name(pub(crate) $native);
        crate::fields::value_field!($native, $name, $input, $class);
        #[wasm_bindgen]
        extern "C" {
            #[wasm_bindgen(typescript_type=$ts)]
            pub type $elements;
        }
        #[wasm_bindgen]
        impl $name {
            #[wasm_bindgen(constructor)]
            pub fn new(length: tsify::Ts<$wire>) -> Result<Self, JsValue> {
                Ok(Self(<$native>::new(input::read::<$index>(length)?)))
            }
            #[wasm_bindgen(getter)]
            pub fn length(&self) -> $index {
                self.0.length()
            }
            pub fn value(&self, index: tsify::Ts<$wire>) -> Result<i32, JsValue> {
                self.0
                    .value(input::read::<$index>(index)?)
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(js_name=setValue)]
            pub fn set_value(
                &mut self,
                index: tsify::Ts<$wire>,
                value: tsify::Ts<input::I32Input>,
            ) -> Result<(), JsValue> {
                self.0
                    .set_value(input::read::<$index>(index)?, input::read::<i32>(value)?)
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(getter,js_name=nonzeroElements)]
            pub fn nonzero_elements(&self) -> Result<$elements, JsValue> {
                let map = js_sys::Map::new();
                for (&index, &count) in self.0.nonzero_elements() {
                    map.set(&JsValue::from(index), &JsValue::from(count));
                }
                Ok(crate::fields::readonly_map(&map)?.unchecked_into())
            }
            #[wasm_bindgen(js_name=totalValue)]
            pub fn total_value(
                &self,
                use_abs: Option<tsify::Ts<input::BoolInput>>,
            ) -> Result<i32, JsValue> {
                self.0
                    .total_value(
                        use_abs
                            .map(input::read::<bool>)
                            .transpose()?
                            .unwrap_or(false),
                    )
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(js_name=fuzzyAnd)]
            pub fn fuzzy_and(&self, other: $input) -> Result<Self, JsValue> {
                Self::with_ref(&other, |other| {
                    self.0
                        .fuzzy_and(other)
                        .map(Self)
                        .map_err(FingerprintError::from_rust)
                })
            }
            #[wasm_bindgen(js_name=fuzzyOr)]
            pub fn fuzzy_or(&self, other: $input) -> Result<Self, JsValue> {
                Self::with_ref(&other, |other| {
                    self.0
                        .fuzzy_or(other)
                        .map(Self)
                        .map_err(FingerprintError::from_rust)
                })
            }
            #[wasm_bindgen(js_name=withAdded)]
            pub fn with_added(&self, other: $input) -> Result<Self, JsValue> {
                Self::with_ref(&other, |other| {
                    self.0
                        .with_added(other)
                        .map(Self)
                        .map_err(FingerprintError::from_rust)
                })
            }
            #[wasm_bindgen(js_name=withSubtracted)]
            pub fn with_subtracted(&self, other: $input) -> Result<Self, JsValue> {
                Self::with_ref(&other, |other| {
                    self.0
                        .with_subtracted(other)
                        .map(Self)
                        .map_err(FingerprintError::from_rust)
                })
            }
            #[wasm_bindgen(js_name=withAddedScalar)]
            pub fn with_added_scalar(
                &self,
                value: tsify::Ts<input::I32Input>,
            ) -> Result<Self, JsValue> {
                self.0
                    .with_added_scalar(input::read::<i32>(value)?)
                    .map(Self)
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(js_name=withSubtractedScalar)]
            pub fn with_subtracted_scalar(
                &self,
                value: tsify::Ts<input::I32Input>,
            ) -> Result<Self, JsValue> {
                self.0
                    .with_subtracted_scalar(input::read::<i32>(value)?)
                    .map(Self)
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(js_name=withMultipliedScalar)]
            pub fn with_multiplied_scalar(
                &self,
                value: tsify::Ts<input::I32Input>,
            ) -> Result<Self, JsValue> {
                self.0
                    .with_multiplied_scalar(input::read::<i32>(value)?)
                    .map(Self)
                    .map_err(FingerprintError::from_rust)
            }
            #[wasm_bindgen(js_name=withDividedScalar)]
            pub fn with_divided_scalar(
                &self,
                value: tsify::Ts<input::I32Input>,
            ) -> Result<Self, JsValue> {
                self.0
                    .with_divided_scalar(input::read::<i32>(value)?)
                    .map(Self)
                    .map_err(FingerprintError::from_rust)
            }
        }
    };
}
sparse_counts!(
    SparseCountFingerprint,
    SparseCountFingerprintInput,
    SparseCountElements,
    ck::SparseCountFingerprint,
    u64,
    input::U64Input,
    "ReadonlyMap<bigint, number>",
    "SparseCountFingerprint"
);
sparse_counts!(
    SparseCountFingerprint32,
    SparseCountFingerprint32Input,
    SparseCountElements32,
    ck::SparseCountFingerprint32,
    u32,
    input::U32Input,
    "ReadonlyMap<number, number>",
    "SparseCountFingerprint32"
);
