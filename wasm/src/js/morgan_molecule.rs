//! All registered scalar Morgan calls, persistent generators and shared settings.
use crate::Molecule;
use crate::fingerprint_additional::with_output;
use crate::fingerprint_source_errors::morgan_error;
use crate::fingerprint_values::{
    Fingerprint, SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32, u32_array,
};
use crate::host_values::{bool_value, i32_value, sequence, type_error, u32_value};
use crate::morgan_parameters::{
    MorganAtomInvariantsGenerator, MorganBondInvariantsGenerator, MorganCallParams,
    MorganFingerprintParams, MorganParams,
};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
fn error(e: ck::MorganReadError) -> JsValue {
    morgan_error(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintMorganWithGenerator)]
    pub fn fingerprint_morgan_with_generator(
        &self,
        generator: &MorganFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganCallParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_morgan_with_generator(&generator.inner, params, output)
                    .map(|inner| Fingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganCountWithGenerator)]
    pub fn fingerprint_morgan_count_with_generator(
        &self,
        generator: &MorganFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganCallParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .morgan_count_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_morgan_count_with_generator(&generator.inner, params, output)
                    .map(|inner| SparseCountFingerprint32 { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparseWithGenerator)]
    pub fn fingerprint_morgan_sparse_with_generator(
        &self,
        generator: &MorganFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganCallParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_morgan_sparse_with_generator(&generator.inner, params, output)
                    .map(|inner| SparseBitFingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparseCountWithGenerator)]
    pub fn fingerprint_morgan_sparse_count_with_generator(
        &self,
        generator: &MorganFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganCallParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_morgan_sparse_count_with_generator(
                        &generator.inner,
                        params,
                        output,
                    )
                    .map(|inner| SparseCountFingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparseCount)]
    pub fn fingerprint_morgan_sparse_count(&self) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint(
        self.inner
            .fingerprint_morgan_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparseCountWithParams)]
    pub fn fingerprint_morgan_sparse_count_with_params(
        &self,
        params: &MorganFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_morgan_sparse_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparse)]
    pub fn fingerprint_morgan_sparse(&self) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint(
        self.inner
            .fingerprint_morgan_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintMorganSparseWithParams)]
    pub fn fingerprint_morgan_sparse_with_params(
        &self,
        params: &MorganFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_sparse_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_morgan_sparse_with_params(&params.inner, output)
                .map(|inner| SparseBitFingerprint { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorganCount)]
    pub fn fingerprint_morgan_count(&self) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .morgan_count_fingerprint(
        self.inner
            .fingerprint_morgan_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintMorganCountWithParams)]
    pub fn fingerprint_morgan_count_with_params(
        &self,
        params: &MorganFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .morgan_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_morgan_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint32 { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintMorgan)]
    pub fn fingerprint_morgan(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_fingerprint(
        self.inner
            .fingerprint_morgan()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintMorganWithParams)]
    pub fn fingerprint_morgan_with_params(
        &self,
        params: &MorganFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .morgan_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_morgan_with_params(&params.inner, output)
                .map(|inner| Fingerprint { inner })
                .map_err(error)
        })
    }
}
#[wasm_bindgen(
    inline_js = "export function visitScalarMorganparams(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid MorganParams\",{cause});}} export function visitScalarMorganatom(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid MorganAtomInvariantsGenerator\",{cause});}} export function visitScalarMorganbond(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid MorganBondInvariantsGenerator\",{cause});}} export function visitScalarMorgancall(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid MorganCallParams\",{cause});}} export function visitScalarMorganMolecule(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid Molecule\",{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitScalarMorganparams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&MorganParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarMorganatom)]
    fn visit_atom(
        v: &JsValue,
        f: &mut dyn FnMut(&MorganAtomInvariantsGenerator),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarMorganbond)]
    fn visit_bond(
        v: &JsValue,
        f: &mut dyn FnMut(&MorganBondInvariantsGenerator),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarMorgancall)]
    fn visit_call(v: &JsValue, f: &mut dyn FnMut(&MorganCallParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarMorganMolecule)]
    fn visit_molecule(v: &JsValue, f: &mut dyn FnMut(&Molecule)) -> Result<(), JsValue>;
}
fn with_params<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::MorganParams>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_params(value, &mut |v: &MorganParams| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("MorganParams"))?
}
fn with_atom<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::MorganAtomInvariantsGenerator>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_atom(value, &mut |v: &MorganAtomInvariantsGenerator| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("MorganAtomInvariantsGenerator"))?
}
fn with_bond<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::MorganBondInvariantsGenerator>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_bond(value, &mut |v: &MorganBondInvariantsGenerator| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("MorganBondInvariantsGenerator"))?
}
fn with_call<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::MorganCallParams>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_call(value, &mut |v: &MorganCallParams| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("MorganCallParams"))?
}
#[wasm_bindgen]
pub struct MorganFingerprintGenerator {
    inner: ck::MorganFingerprintGenerator,
}
#[wasm_bindgen]
impl MorganFingerprintGenerator {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "MorganParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganBondInvariantsGenerator | null")]
        bond_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(params, atom_invariants, bond_invariants)
    }
    #[wasm_bindgen(js_name=new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "MorganParams | null")] params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "MorganBondInvariantsGenerator | null")]
        bond_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::MorganFingerprintGenerator::new(
        with_params(&params, |params| {
            with_atom(&atom_invariants, |atom| {
                with_bond(&bond_invariants, |bond| {
                    ck::MorganFingerprintGenerator::new(params, atom, bond)
                        .map(|inner| Self { inner })
                        .map_err(error)
                })
            })
        })
    }
    #[wasm_bindgen(js_name=fromJson)]
    pub fn from_json(json: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::MorganFingerprintGenerator::from_json(json)
        ck::MorganFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    pub fn settings(&self) -> MorganSettings {
        // COSMolKit❗✔️: inner: self.inner.settings(),
        MorganSettings {
            inner: self.inner.settings(),
        }
    }
    #[wasm_bindgen(js_name=infoString)]
    pub fn info_string(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.info_string()
        self.inner.info_string().map_err(error)
    }
    #[wasm_bindgen(js_name=toJson)]
    pub fn to_json(&self) -> Result<String, JsValue> {
        // COSMolKit❗✔️: self.inner.to_json()
        self.inner
            .to_json()
            .map_err(error)
            .and_then(|value| crate::host_values::text(&value))
    }
    #[wasm_bindgen(js_name=fingerprints,unchecked_return_type="(Fingerprint | null)[]")]
    pub fn fingerprints(
        &self,
        #[wasm_bindgen(unchecked_param_type = "(Molecule | null)[]")] molecules: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: self.inner.fingerprints(&rows, num_threads)
        let mut handles: Vec<Option<Arc<cosmolkit_wasm::Molecule>>> = Vec::new();
        for value in sequence(&molecules, "molecules")?.iter() {
            if value.is_null() {
                handles.push(None);
            } else {
                visit_molecule(&value, &mut |m: &Molecule| {
                    handles.push(Some(Arc::clone(&m.inner)))
                })?;
            }
        }
        let rows = handles.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
        cosmolkit_wasm::morgan_generator_fingerprints(
            &self.inner,
            &rows,
            if num_threads.is_undefined() {
                1
            } else {
                i32_value(&num_threads, "numThreads")?
            },
        )
        .map_err(error)
        .map(|values| {
            values
                .into_iter()
                .map(|v| v.map_or(JsValue::NULL, |inner| Fingerprint { inner }.into()))
                .collect()
        })
    }
    #[wasm_bindgen(js_name=counts,unchecked_return_type="(SparseCountFingerprint32 | null)[]")]
    pub fn counts(
        &self,
        #[wasm_bindgen(unchecked_param_type = "(Molecule | null)[]")] molecules: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: self.inner.counts(&rows, num_threads)
        let mut handles: Vec<Option<Arc<cosmolkit_wasm::Molecule>>> = Vec::new();
        for value in sequence(&molecules, "molecules")?.iter() {
            if value.is_null() {
                handles.push(None);
            } else {
                visit_molecule(&value, &mut |m: &Molecule| {
                    handles.push(Some(Arc::clone(&m.inner)))
                })?;
            }
        }
        let rows = handles.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
        cosmolkit_wasm::morgan_generator_counts(
            &self.inner,
            &rows,
            if num_threads.is_undefined() {
                1
            } else {
                i32_value(&num_threads, "numThreads")?
            },
        )
        .map_err(error)
        .map(|values| {
            values
                .into_iter()
                .map(|v| {
                    v.map_or(JsValue::NULL, |inner| {
                        SparseCountFingerprint32 { inner }.into()
                    })
                })
                .collect()
        })
    }
    #[wasm_bindgen(js_name=sparseFingerprints,unchecked_return_type="(SparseBitFingerprint | null)[]")]
    pub fn sparse_fingerprints(
        &self,
        #[wasm_bindgen(unchecked_param_type = "(Molecule | null)[]")] molecules: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: self.inner.sparse_fingerprints(&rows, num_threads)
        let mut handles: Vec<Option<Arc<cosmolkit_wasm::Molecule>>> = Vec::new();
        for value in sequence(&molecules, "molecules")?.iter() {
            if value.is_null() {
                handles.push(None);
            } else {
                visit_molecule(&value, &mut |m: &Molecule| {
                    handles.push(Some(Arc::clone(&m.inner)))
                })?;
            }
        }
        let rows = handles.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
        cosmolkit_wasm::morgan_generator_sparse_fingerprints(
            &self.inner,
            &rows,
            if num_threads.is_undefined() {
                1
            } else {
                i32_value(&num_threads, "numThreads")?
            },
        )
        .map_err(error)
        .map(|values| {
            values
                .into_iter()
                .map(|v| v.map_or(JsValue::NULL, |inner| SparseBitFingerprint { inner }.into()))
                .collect()
        })
    }
    #[wasm_bindgen(js_name=sparseCounts,unchecked_return_type="(SparseCountFingerprint | null)[]")]
    pub fn sparse_counts(
        &self,
        #[wasm_bindgen(unchecked_param_type = "(Molecule | null)[]")] molecules: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
    ) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: self.inner.sparse_counts(&rows, num_threads)
        let mut handles: Vec<Option<Arc<cosmolkit_wasm::Molecule>>> = Vec::new();
        for value in sequence(&molecules, "molecules")?.iter() {
            if value.is_null() {
                handles.push(None);
            } else {
                visit_molecule(&value, &mut |m: &Molecule| {
                    handles.push(Some(Arc::clone(&m.inner)))
                })?;
            }
        }
        let rows = handles.iter().map(|m| m.as_deref()).collect::<Vec<_>>();
        cosmolkit_wasm::morgan_generator_sparse_counts(
            &self.inner,
            &rows,
            if num_threads.is_undefined() {
                1
            } else {
                i32_value(&num_threads, "numThreads")?
            },
        )
        .map_err(error)
        .map(|values| {
            values
                .into_iter()
                .map(|v| {
                    v.map_or(JsValue::NULL, |inner| {
                        SparseCountFingerprint { inner }.into()
                    })
                })
                .collect()
        })
    }
}
#[wasm_bindgen]
pub struct MorganSettings {
    inner: ck::MorganSettings,
}
#[wasm_bindgen]
impl MorganSettings {
    #[wasm_bindgen(getter,js_name=radius)]
    pub fn radius(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: self.inner.radius()
        self.inner.radius().map_err(error)
    }
    #[wasm_bindgen(js_name=setRadius)]
    pub fn set_radius(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_radius(value)
        self.inner
            .set_radius(u32_value(&value, "radius")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=radius)]
    pub fn assign_radius(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_radius(value)
    }
    #[wasm_bindgen(getter,js_name=onlyNonzeroInvariants)]
    pub fn only_nonzero_invariants(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.only_nonzero_invariants()
        self.inner.only_nonzero_invariants().map_err(error)
    }
    #[wasm_bindgen(js_name=setOnlyNonzeroInvariants)]
    pub fn set_only_nonzero_invariants(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_only_nonzero_invariants(value)
        self.inner
            .set_only_nonzero_invariants(bool_value(&value, "onlyNonzeroInvariants")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=onlyNonzeroInvariants)]
    pub fn assign_only_nonzero_invariants(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_only_nonzero_invariants(value)
    }
    #[wasm_bindgen(getter,js_name=includeRedundantEnvironments)]
    pub fn include_redundant_environments(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.include_redundant_environments()
        self.inner.include_redundant_environments().map_err(error)
    }
    #[wasm_bindgen(js_name=setIncludeRedundantEnvironments)]
    pub fn set_include_redundant_environments(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_include_redundant_environments(value)
        self.inner
            .set_include_redundant_environments(bool_value(&value, "includeRedundantEnvironments")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=includeRedundantEnvironments)]
    pub fn assign_include_redundant_environments(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_include_redundant_environments(value)
    }
    #[wasm_bindgen(getter,js_name=includeChirality)]
    pub fn include_chirality(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.include_chirality()
        self.inner.include_chirality().map_err(error)
    }
    #[wasm_bindgen(js_name=setIncludeChirality)]
    pub fn set_include_chirality(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_include_chirality(value)
        self.inner
            .set_include_chirality(bool_value(&value, "includeChirality")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=includeChirality)]
    pub fn assign_include_chirality(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_include_chirality(value)
    }
    #[wasm_bindgen(getter,js_name=countSimulation)]
    pub fn count_simulation(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.count_simulation()
        self.inner.count_simulation().map_err(error)
    }
    #[wasm_bindgen(js_name=setCountSimulation)]
    pub fn set_count_simulation(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_count_simulation(value)
        self.inner
            .set_count_simulation(bool_value(&value, "countSimulation")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=countSimulation)]
    pub fn assign_count_simulation(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_count_simulation(value)
    }
    #[wasm_bindgen(getter,js_name=fpSize)]
    pub fn fp_size(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: self.inner.fp_size()
        self.inner.fp_size().map_err(error)
    }
    #[wasm_bindgen(js_name=setFpSize)]
    pub fn set_fp_size(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_fp_size(value)
        self.inner
            .set_fp_size(u32_value(&value, "fpSize")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=fpSize)]
    pub fn assign_fp_size(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_fp_size(value)
    }
    #[wasm_bindgen(getter,js_name=bitsPerFeature)]
    pub fn bits_per_feature(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: self.inner.bits_per_feature()
        self.inner.bits_per_feature().map_err(error)
    }
    #[wasm_bindgen(js_name=setBitsPerFeature)]
    pub fn set_bits_per_feature(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_bits_per_feature(value)
        self.inner
            .set_bits_per_feature(u32_value(&value, "bitsPerFeature")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=bitsPerFeature)]
    pub fn assign_bits_per_feature(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_bits_per_feature(value)
    }
    #[wasm_bindgen(getter,js_name=countBounds,unchecked_return_type="number[]")]
    pub fn count_bounds(&self) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: self.inner.count_bounds()
        self.inner
            .count_bounds()
            .map_err(error)
            .map(|rows| rows.into_iter().map(JsValue::from).collect())
    }
    #[wasm_bindgen(js_name=setCountBounds)]
    pub fn set_count_bounds(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_count_bounds(value)
        self.inner
            .set_count_bounds(u32_array(&value, "countBounds")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=countBounds)]
    pub fn assign_count_bounds(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_count_bounds(value)
    }
    pub fn params(&self) -> Result<MorganParams, JsValue> {
        // COSMolKit❗✔️: .map(|inner| MorganParams { inner })
        self.inner
            .params()
            .map(|inner| MorganParams { inner })
            .map_err(error)
    }
}
