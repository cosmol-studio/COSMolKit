//! Detached value transport; every operation delegates to the public facade.
use crate::fingerprint_errors::throw_fingerprint_error as error;
use crate::host_values::{bool_value, i32_value, sequence, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use wasm_bindgen::prelude::*;
pub(crate) fn u64_value(value: &JsValue, name: &str) -> Result<u64, JsValue> {
    if !value.is_bigint() {
        return Err(type_error(name));
    }
    u64::try_from(value.clone())
        .map_err(|_| js_sys::RangeError::new(&format!("{name} is outside u64 range")).into())
}
pub(crate) fn u32_array(value: &JsValue, name: &str) -> Result<Vec<u32>, JsValue> {
    sequence(value, name)?
        .iter()
        .map(|v| u32_value(&v, name))
        .collect()
}
#[wasm_bindgen]
pub struct Fingerprint {
    pub(crate) inner: ck::Fingerprint,
}
#[wasm_bindgen]
impl Fingerprint {
    #[wasm_bindgen(js_name=fromOnBits)]
    pub fn from_on_bits(
        #[wasm_bindgen(unchecked_param_type = "number")] n_bits: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] on_bits: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs::from_on_bits:
        //         ck::Fingerprint::from_on_bits(n_bits, on_bits)
        //             .map(|inner| Self { inner })
        ck::Fingerprint::from_on_bits(u32_value(&n_bits, "nBits")?, u32_array(&on_bits, "onBits")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=nBits)]
    pub fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    #[wasm_bindgen(js_name=onBits,unchecked_return_type="number[]")]
    pub fn on_bits(&self) -> Array {
        self.inner
            .on_bits()
            .into_iter()
            .map(JsValue::from)
            .collect()
    }
    pub fn tanimoto(&self, other: &Fingerprint) -> Result<f64, JsValue> {
        self.inner.tanimoto(&other.inner).map_err(error)
    }
}
#[wasm_bindgen]
pub struct SparseBitFingerprint {
    pub(crate) inner: ck::SparseBitFingerprint,
}
#[wasm_bindgen]
impl SparseBitFingerprint {
    #[wasm_bindgen(js_name=nBits)]
    pub fn n_bits(&self) -> u32 {
        self.inner.n_bits()
    }
    #[wasm_bindgen(js_name=onBits,unchecked_return_type="number[]")]
    pub fn on_bits(&self) -> Array {
        self.inner
            .on_bits()
            .into_iter()
            .map(JsValue::from)
            .collect()
    }
}

#[wasm_bindgen]
pub struct SparseCountFingerprint {
    pub(crate) inner: ck::SparseCountFingerprint,
}
#[wasm_bindgen]
impl SparseCountFingerprint {
    #[wasm_bindgen(js_name=new)]
    pub fn new_value(
        #[wasm_bindgen(unchecked_param_type = "bigint")] length: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        Ok(Self {
            inner: ck::SparseCountFingerprint::new(u64_value(&length, "length")?),
        })
    }
    pub fn length(&self) -> u64 {
        self.inner.length()
    }
    pub fn value(
        &self,
        #[wasm_bindgen(unchecked_param_type = "bigint")] index: JsValue,
    ) -> Result<i32, JsValue> {
        self.inner.value(u64_value(&index, "index")?).map_err(error)
    }
    #[wasm_bindgen(js_name=setValue)]
    pub fn set_value(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "bigint")] index: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_value(u64_value(&index, "index")?, i32_value(&value, "value")?)
            .map_err(error)
    }
    #[wasm_bindgen(js_name=nonzeroElements,unchecked_return_type="Map<bigint, number>")]
    pub fn nonzero_elements(&self) -> Map {
        let map = Map::new();
        for (&key, &value) in self.inner.nonzero_elements() {
            map.set(&key.into(), &value.into());
        }
        map
    }
    #[wasm_bindgen(js_name=totalValue)]
    pub fn total_value(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_abs: JsValue,
    ) -> Result<i32, JsValue> {
        self.inner
            .total_value(if use_abs.is_undefined() {
                false
            } else {
                bool_value(&use_abs, "useAbs")?
            })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fuzzyAnd)]
    pub fn fuzzy_and(
        &self,
        other: &SparseCountFingerprint,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.fuzzy_and(&other.inner)
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fuzzyOr)]
    pub fn fuzzy_or(
        &self,
        other: &SparseCountFingerprint,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.fuzzy_or(&other.inner)
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withAdded)]
    pub fn with_added(
        &self,
        other: &SparseCountFingerprint,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_added(&other.inner)
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withSubtracted)]
    pub fn with_subtracted(
        &self,
        other: &SparseCountFingerprint,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_subtracted(&other.inner)
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withAddedScalar)]
    pub fn with_added_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_added_scalar(value)
        self.inner
            .with_added_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withSubtractedScalar)]
    pub fn with_subtracted_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_subtracted_scalar(value)
        self.inner
            .with_subtracted_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withMultipliedScalar)]
    pub fn with_multiplied_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_multiplied_scalar(value)
        self.inner
            .with_multiplied_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withDividedScalar)]
    pub fn with_divided_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_divided_scalar(value)
        self.inner
            .with_divided_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
}

#[wasm_bindgen]
pub struct SparseCountFingerprint32 {
    pub(crate) inner: ck::SparseCountFingerprint32,
}
#[wasm_bindgen]
impl SparseCountFingerprint32 {
    #[wasm_bindgen(js_name=new)]
    pub fn new_value(
        #[wasm_bindgen(unchecked_param_type = "number")] length: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        Ok(Self {
            inner: ck::SparseCountFingerprint32::new(u32_value(&length, "length")?),
        })
    }
    pub fn length(&self) -> u32 {
        self.inner.length()
    }
    pub fn value(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<i32, JsValue> {
        self.inner.value(u32_value(&index, "index")?).map_err(error)
    }
    #[wasm_bindgen(js_name=setValue)]
    pub fn set_value(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.inner
            .set_value(u32_value(&index, "index")?, i32_value(&value, "value")?)
            .map_err(error)
    }
    #[wasm_bindgen(js_name=nonzeroElements,unchecked_return_type="Map<number, number>")]
    pub fn nonzero_elements(&self) -> Map {
        let map = Map::new();
        for (&key, &value) in self.inner.nonzero_elements() {
            map.set(&key.into(), &value.into());
        }
        map
    }
    #[wasm_bindgen(js_name=totalValue)]
    pub fn total_value(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")] use_abs: JsValue,
    ) -> Result<i32, JsValue> {
        self.inner
            .total_value(if use_abs.is_undefined() {
                false
            } else {
                bool_value(&use_abs, "useAbs")?
            })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fuzzyAnd)]
    pub fn fuzzy_and(
        &self,
        other: &SparseCountFingerprint32,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.fuzzy_and(&other.inner)
        self.inner
            .fuzzy_and(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fuzzyOr)]
    pub fn fuzzy_or(
        &self,
        other: &SparseCountFingerprint32,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.fuzzy_or(&other.inner)
        self.inner
            .fuzzy_or(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withAdded)]
    pub fn with_added(
        &self,
        other: &SparseCountFingerprint32,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_added(&other.inner)
        self.inner
            .with_added(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withSubtracted)]
    pub fn with_subtracted(
        &self,
        other: &SparseCountFingerprint32,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_subtracted(&other.inner)
        self.inner
            .with_subtracted(&other.inner)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withAddedScalar)]
    pub fn with_added_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_added_scalar(value)
        self.inner
            .with_added_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withSubtractedScalar)]
    pub fn with_subtracted_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_subtracted_scalar(value)
        self.inner
            .with_subtracted_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withMultipliedScalar)]
    pub fn with_multiplied_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_multiplied_scalar(value)
        self.inner
            .with_multiplied_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=withDividedScalar)]
    pub fn with_divided_scalar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: canonical_values.rs: self.inner.with_divided_scalar(value)
        self.inner
            .with_divided_scalar(i32_value(&value, "value")?)
            .map(|inner| Self { inner })
            .map_err(error)
    }
}
