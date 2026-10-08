//! Unique mutable output and detached public read snapshots.
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Map};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct FingerprintAdditionalOutput {
    pub(crate) inner: std::cell::RefCell<ck::FingerprintAdditionalOutput>,
}
#[wasm_bindgen]
impl FingerprintAdditionalOutput {
    #[wasm_bindgen(constructor)]
    pub fn constructor() -> Self {
        // COSMolKit❗✔️: inner: ck::FingerprintAdditionalOutput::new(),
        Self {
            inner: ck::FingerprintAdditionalOutput::new().into(),
        }
    }
    #[wasm_bindgen(js_name=new)]
    pub fn new_value() -> Self {
        // COSMolKit❗✔️: inner: ck::FingerprintAdditionalOutput::new(),
        Self {
            inner: ck::FingerprintAdditionalOutput::new().into(),
        }
    }
    pub fn default() -> Self {
        // COSMolKit❗✔️: inner: ck::FingerprintAdditionalOutput::default(),
        Self {
            inner: ck::FingerprintAdditionalOutput::default().into(),
        }
    }

    #[wasm_bindgen(js_name=atomCounts,unchecked_return_type="number[] | null")]
    pub fn atom_counts(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.atom_counts()
        self.inner
            .borrow()
            .atom_counts()
            .map_or(JsValue::NULL, |v| {
                v.iter()
                    .copied()
                    .map(JsValue::from)
                    .collect::<Array>()
                    .into()
            })
    }
    #[wasm_bindgen(js_name=atomToBits,unchecked_return_type="bigint[][] | null")]
    pub fn atom_to_bits(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.atom_to_bits()
        self.inner
            .borrow()
            .atom_to_bits()
            .map_or(JsValue::NULL, |v| {
                v.iter()
                    .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                    .collect::<Array>()
                    .into()
            })
    }
    #[wasm_bindgen(js_name=bitInfoMap,unchecked_return_type="Map<bigint, [number, number][]> | null")]
    pub fn bit_info_map(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.bit_info_map()
        self.inner
            .borrow()
            .bit_info_map()
            .map_or(JsValue::NULL, |v| {
                let map = Map::new();
                for (&key, rows) in v {
                    let rows = rows
                        .iter()
                        .map(|&(a, b)| Array::of2(&a.into(), &b.into()))
                        .collect::<Array>();
                    map.set(&key.into(), &rows);
                }
                map.into()
            })
    }
    #[wasm_bindgen(js_name=bitPaths,unchecked_return_type="Map<bigint, number[][]> | null")]
    pub fn bit_paths(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.bit_paths()
        self.inner.borrow().bit_paths().map_or(JsValue::NULL, |v| {
            let map = Map::new();
            for (&key, rows) in v {
                let rows = rows
                    .iter()
                    .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                    .collect::<Array>();
                map.set(&key.into(), &rows);
            }
            map.into()
        })
    }
    #[wasm_bindgen(js_name=atomsPerBit,unchecked_return_type="Map<bigint, number[][]> | null")]
    pub fn atoms_per_bit(&self) -> JsValue {
        // COSMolKit❗✔️: self.inner.atoms_per_bit()
        self.inner
            .borrow()
            .atoms_per_bit()
            .map_or(JsValue::NULL, |v| {
                let map = Map::new();
                for (&key, rows) in v {
                    let rows = rows
                        .iter()
                        .map(|row| row.iter().copied().map(JsValue::from).collect::<Array>())
                        .collect::<Array>();
                    map.set(&key.into(), &rows);
                }
                map.into()
            })
    }
}
#[wasm_bindgen]
impl FingerprintAdditionalOutput {
    #[wasm_bindgen(js_name=allocateAtomCounts)]
    pub fn allocate_atom_counts(&mut self) {
        // COSMolKit❗✔️: self.inner.allocate_atom_counts();
        self.inner.get_mut().allocate_atom_counts();
    }
    #[wasm_bindgen(js_name=allocateAtomToBits)]
    pub fn allocate_atom_to_bits(&mut self) {
        // COSMolKit❗✔️: self.inner.allocate_atom_to_bits();
        self.inner.get_mut().allocate_atom_to_bits();
    }
    #[wasm_bindgen(js_name=allocateBitInfoMap)]
    pub fn allocate_bit_info_map(&mut self) {
        // COSMolKit❗✔️: self.inner.allocate_bit_info_map();
        self.inner.get_mut().allocate_bit_info_map();
    }
    #[wasm_bindgen(js_name=allocateBitPaths)]
    pub fn allocate_bit_paths(&mut self) {
        // COSMolKit❗✔️: self.inner.allocate_bit_paths();
        self.inner.get_mut().allocate_bit_paths();
    }
    #[wasm_bindgen(js_name=allocateAtomsPerBit)]
    pub fn allocate_atoms_per_bit(&mut self) {
        // COSMolKit❗✔️: self.inner.allocate_atoms_per_bit();
        self.inner.get_mut().allocate_atoms_per_bit();
    }
}
#[wasm_bindgen(
    inline_js = "export function visitFingerprintAdditionalOutput(value,visit){try{visit(value);}catch(cause){throw new TypeError('invalid FingerprintAdditionalOutput',{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitFingerprintAdditionalOutput)]
    fn visit_output(
        value: &JsValue,
        visit: &mut dyn FnMut(&FingerprintAdditionalOutput),
    ) -> Result<(), JsValue>;
}
pub(crate) fn with_output<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&mut ck::FingerprintAdditionalOutput>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() {
        return run(None);
    }
    if value.is_undefined() {
        return Err(crate::host_values::type_error("additionalOutput"));
    }
    let mut run = Some(run);
    let mut result = None;
    visit_output(value, &mut |output: &FingerprintAdditionalOutput| {
        if let Some(run) = run.take() {
            result = Some(run(Some(&mut output.inner.borrow_mut())));
        }
    })?;
    result.ok_or_else(|| crate::host_values::type_error("additionalOutput"))?
}
