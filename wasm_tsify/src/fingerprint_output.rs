//! Native output owner; allocation and collection semantics are not rebuilt.
use cosmolkit as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct FingerprintAdditionalOutput(
    pub(crate) crate::storage::NativeField<ck::FingerprintAdditionalOutput>,
);
crate::storage::mutable_owner!(
    FingerprintAdditionalOutput,
    ck::FingerprintAdditionalOutput,
    FingerprintAdditionalOutputInput,
    "FingerprintAdditionalOutput"
);
impl FingerprintAdditionalOutput {
    pub(crate) fn with_optional<R>(
        value: Option<&FingerprintAdditionalOutputInput>,
        apply: impl FnOnce(Option<&mut ck::FingerprintAdditionalOutput>) -> Result<R, JsValue>,
    ) -> Result<R, JsValue> {
        match value {
            Some(value) => Self::with_mut(value, |value| apply(Some(value))),
            None => apply(None),
        }
    }
}
#[wasm_bindgen]
extern "C" {
    #[wasm_bindgen(typescript_type = "ReadonlyArray<number>")]
    pub type AtomCounts;
    #[wasm_bindgen(typescript_type = "ReadonlyArray<ReadonlyArray<bigint>>")]
    pub type AtomToBits;
    #[wasm_bindgen(
        typescript_type = "ReadonlyMap<bigint, ReadonlyArray<readonly [number, number]>>"
    )]
    pub type BitInfoMap;
    #[wasm_bindgen(typescript_type = "ReadonlyMap<bigint, ReadonlyArray<ReadonlyArray<number>>>")]
    pub type BitPaths;
}
pub(crate) fn array(values: impl IntoIterator<Item = JsValue>) -> JsValue {
    let output = values.into_iter().collect::<js_sys::Array>();
    js_sys::Object::freeze(&output);
    output.into()
}
pub(crate) fn paths(
    values: &std::collections::BTreeMap<u64, Vec<Vec<i32>>>,
) -> Result<BitPaths, JsValue> {
    let output = js_sys::Map::new();
    for (&key, rows) in values {
        output.set(
            &JsValue::from(key),
            &array(
                rows.iter()
                    .map(|row| array(row.iter().copied().map(JsValue::from))),
            ),
        );
    }
    Ok(crate::fields::readonly_map(&output)?.unchecked_into())
}
pub(crate) fn info_map(
    values: &std::collections::BTreeMap<u64, Vec<(u32, u32)>>,
) -> Result<BitInfoMap, JsValue> {
    let output = js_sys::Map::new();
    for (&key, rows) in values {
        output.set(
            &JsValue::from(key),
            &array(rows.iter().map(|&(a, b)| array([a.into(), b.into()]))),
        );
    }
    Ok(crate::fields::readonly_map(&output)?.unchecked_into())
}
#[wasm_bindgen]
impl FingerprintAdditionalOutput {
    #[wasm_bindgen(constructor)]
    pub fn new() -> Self {
        Self(crate::storage::NativeField::owned(
            ck::FingerprintAdditionalOutput::new(),
        ))
    }
    pub fn default() -> Self {
        Self(crate::storage::NativeField::owned(
            ck::FingerprintAdditionalOutput::default(),
        ))
    }
    #[wasm_bindgen(js_name=allocateAtomCounts)]
    pub fn allocate_atom_counts(&self) -> Result<(), JsValue> {
        self.0
            .write(ck::FingerprintAdditionalOutput::allocate_atom_counts)
    }
    #[wasm_bindgen(js_name=allocateAtomToBits)]
    pub fn allocate_atom_to_bits(&self) -> Result<(), JsValue> {
        self.0
            .write(ck::FingerprintAdditionalOutput::allocate_atom_to_bits)
    }
    #[wasm_bindgen(js_name=allocateBitInfoMap)]
    pub fn allocate_bit_info_map(&self) -> Result<(), JsValue> {
        self.0
            .write(ck::FingerprintAdditionalOutput::allocate_bit_info_map)
    }
    #[wasm_bindgen(js_name=allocateBitPaths)]
    pub fn allocate_bit_paths(&self) -> Result<(), JsValue> {
        self.0
            .write(ck::FingerprintAdditionalOutput::allocate_bit_paths)
    }
    #[wasm_bindgen(js_name=allocateAtomsPerBit)]
    pub fn allocate_atoms_per_bit(&self) -> Result<(), JsValue> {
        self.0
            .write(ck::FingerprintAdditionalOutput::allocate_atoms_per_bit)
    }
    #[wasm_bindgen(getter,js_name=atomCounts)]
    pub fn atom_counts(&self) -> Result<Option<AtomCounts>, JsValue> {
        self.0.read(|native| {
            native
                .atom_counts()
                .map(|values| array(values.iter().copied().map(JsValue::from)).unchecked_into())
        })
    }
    #[wasm_bindgen(getter,js_name=atomToBits)]
    pub fn atom_to_bits(&self) -> Result<Option<AtomToBits>, JsValue> {
        self.0.read(|native| {
            native.atom_to_bits().map(|values| {
                array(
                    values
                        .iter()
                        .map(|row| array(row.iter().copied().map(JsValue::from))),
                )
                .unchecked_into()
            })
        })
    }
    #[wasm_bindgen(getter,js_name=bitInfoMap)]
    pub fn bit_info_map(&self) -> Result<Option<BitInfoMap>, JsValue> {
        self.0
            .read(|native| native.bit_info_map().map(info_map))?
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=bitPaths)]
    pub fn bit_paths(&self) -> Result<Option<BitPaths>, JsValue> {
        self.0
            .read(|native| native.bit_paths().map(paths))?
            .transpose()
    }
    #[wasm_bindgen(getter,js_name=atomsPerBit)]
    pub fn atoms_per_bit(&self) -> Result<Option<BitPaths>, JsValue> {
        self.0
            .read(|native| native.atoms_per_bit().map(paths))?
            .transpose()
    }
}
