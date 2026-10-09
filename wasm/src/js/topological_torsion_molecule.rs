//! All registered scalar TopologicalTorsion calls, persistent generators and shared settings.
use crate::Molecule;
use crate::atom_pair_parameters::AtomPairAtomInvariantsGenerator;
use crate::fingerprint_additional::with_output;
use crate::fingerprint_source_errors::topological_torsion_error;
use crate::fingerprint_values::{
    Fingerprint, SparseBitFingerprint, SparseCountFingerprint, SparseCountFingerprint32, u32_array,
};
use crate::host_values::{bool_value, i32_value, sequence, type_error, u32_value};
use crate::topological_torsion_parameters::{
    TopologicalTorsionCallParams, TopologicalTorsionFingerprintParams, TopologicalTorsionParams,
};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::sync::Arc;
use wasm_bindgen::prelude::*;
fn error(e: ck::TopologicalTorsionReadError) -> JsValue {
    topological_torsion_error(&e).unwrap_or_else(|e| e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionWithGenerator)]
    pub fn fingerprint_topological_torsion_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionCallParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_topological_torsion_with_generator(
                        &generator.inner,
                        params,
                        output,
                    )
                    .map(|inner| Fingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionCountWithGenerator)]
    pub fn fingerprint_topological_torsion_count_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionCallParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_topological_torsion_count_with_generator(
                        &generator.inner,
                        params,
                        output,
                    )
                    .map(|inner| SparseCountFingerprint32 { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseWithGenerator)]
    pub fn fingerprint_topological_torsion_sparse_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionCallParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_topological_torsion_sparse_with_generator(
                        &generator.inner,
                        params,
                        output,
                    )
                    .map(|inner| SparseBitFingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseCountWithGenerator)]
    pub fn fingerprint_topological_torsion_sparse_count_with_generator(
        &self,
        generator: &TopologicalTorsionFingerprintGenerator,
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionCallParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "FingerprintAdditionalOutput | null")]
        output: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint_with_generator(
        let output = if output.is_undefined() {
            JsValue::NULL
        } else {
            output
        };
        with_call(&params, |params| {
            with_output(&output, |output| {
                self.inner
                    .fingerprint_topological_torsion_sparse_count_with_generator(
                        &generator.inner,
                        params,
                        output,
                    )
                    .map(|inner| SparseCountFingerprint { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseCount)]
    pub fn fingerprint_topological_torsion_sparse_count(
        &self,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint(
        self.inner
            .fingerprint_topological_torsion_sparse_count()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseCountWithParams)]
    pub fn fingerprint_topological_torsion_sparse_count_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_topological_torsion_sparse_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparse)]
    pub fn fingerprint_topological_torsion_sparse(&self) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint(
        self.inner
            .fingerprint_topological_torsion_sparse()
            .map(|inner| SparseBitFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseWithParams)]
    pub fn fingerprint_topological_torsion_sparse_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseBitFingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_sparse_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_topological_torsion_sparse_with_params(&params.inner, output)
                .map(|inner| SparseBitFingerprint { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionCount)]
    pub fn fingerprint_topological_torsion_count(
        &self,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint(
        self.inner
            .fingerprint_topological_torsion_count()
            .map(|inner| SparseCountFingerprint32 { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionCountWithParams)]
    pub fn fingerprint_topological_torsion_count_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<SparseCountFingerprint32, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_count_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_topological_torsion_count_with_params(&params.inner, output)
                .map(|inner| SparseCountFingerprint32 { inner })
                .map_err(error)
        })
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsion)]
    pub fn fingerprint_topological_torsion(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint(
        self.inner
            .fingerprint_topological_torsion()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionWithParams)]
    pub fn fingerprint_topological_torsion_with_params(
        &self,
        params: &TopologicalTorsionFingerprintParams,
        #[wasm_bindgen(unchecked_param_type = "FingerprintAdditionalOutput | null")]
        additional_output: JsValue,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_fingerprint_with_params(
        with_output(&additional_output, |output| {
            self.inner
                .fingerprint_topological_torsion_with_params(&params.inner, output)
                .map(|inner| Fingerprint { inner })
                .map_err(error)
        })
    }
}
#[wasm_bindgen(
    inline_js = "export function visitScalarTopologicalTorsionparams(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid TopologicalTorsionParams\",{cause});}} export function visitScalarTopologicalTorsionatom(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid AtomPairAtomInvariantsGenerator\",{cause});}} export function visitScalarTopologicalTorsioncall(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid TopologicalTorsionCallParams\",{cause});}} export function visitScalarTopologicalTorsionMolecule(v,f){try{f(v);}catch(cause){throw new TypeError(\"invalid Molecule\",{cause});}}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitScalarTopologicalTorsionparams)]
    fn visit_params(
        v: &JsValue,
        f: &mut dyn FnMut(&TopologicalTorsionParams),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarTopologicalTorsionatom)]
    fn visit_atom(
        v: &JsValue,
        f: &mut dyn FnMut(&AtomPairAtomInvariantsGenerator),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarTopologicalTorsioncall)]
    fn visit_call(
        v: &JsValue,
        f: &mut dyn FnMut(&TopologicalTorsionCallParams),
    ) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitScalarTopologicalTorsionMolecule)]
    fn visit_molecule(v: &JsValue, f: &mut dyn FnMut(&Molecule)) -> Result<(), JsValue>;
}
fn with_params<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::TopologicalTorsionParams>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_params(value, &mut |v: &TopologicalTorsionParams| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("TopologicalTorsionParams"))?
}
fn with_atom<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::AtomPairAtomInvariantsGenerator>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_atom(value, &mut |v: &AtomPairAtomInvariantsGenerator| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("AtomPairAtomInvariantsGenerator"))?
}
fn with_call<T>(
    value: &JsValue,
    run: impl FnOnce(Option<&ck::TopologicalTorsionCallParams>) -> Result<T, JsValue>,
) -> Result<T, JsValue> {
    if value.is_null() || value.is_undefined() {
        return run(None);
    }
    let mut run = Some(run);
    let mut out = None;
    visit_call(value, &mut |v: &TopologicalTorsionCallParams| {
        if let Some(run) = run.take() {
            out = Some(run(Some(&v.inner)));
        }
    })?;
    out.ok_or_else(|| type_error("TopologicalTorsionCallParams"))?
}
#[wasm_bindgen]
pub struct TopologicalTorsionFingerprintGenerator {
    inner: ck::TopologicalTorsionFingerprintGenerator,
}
#[wasm_bindgen]
impl TopologicalTorsionFingerprintGenerator {
    #[wasm_bindgen(constructor)]
    pub fn construct(
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AtomPairAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        Self::new(params, atom_invariants)
    }
    #[wasm_bindgen(js_name=new)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "TopologicalTorsionParams | null")]
        params: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "AtomPairAtomInvariantsGenerator | null")]
        atom_invariants: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::TopologicalTorsionFingerprintGenerator::new(
        with_params(&params, |params| {
            with_atom(&atom_invariants, |atom| {
                ck::TopologicalTorsionFingerprintGenerator::new(params, atom.copied())
                    .map(|inner| Self { inner })
                    .map_err(error)
            })
        })
    }
    #[wasm_bindgen(js_name=fromJson)]
    pub fn from_json(json: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: ck::TopologicalTorsionFingerprintGenerator::from_json(json)
        ck::TopologicalTorsionFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    pub fn settings(&self) -> TopologicalTorsionSettings {
        // COSMolKit❗✔️: inner: self.inner.settings(),
        TopologicalTorsionSettings {
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
        self.inner.to_json().map_err(error)
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
        cosmolkit_wasm::topological_torsion_generator_fingerprints(
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
        cosmolkit_wasm::topological_torsion_generator_counts(
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
        cosmolkit_wasm::topological_torsion_generator_sparse_fingerprints(
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
        cosmolkit_wasm::topological_torsion_generator_sparse_counts(
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
pub struct TopologicalTorsionSettings {
    inner: ck::TopologicalTorsionSettings,
}
#[wasm_bindgen]
impl TopologicalTorsionSettings {
    #[wasm_bindgen(getter,js_name=torsionAtomCount)]
    pub fn torsion_atom_count(&self) -> Result<u32, JsValue> {
        // COSMolKit❗✔️: self.inner.torsion_atom_count()
        self.inner.torsion_atom_count().map_err(error)
    }
    #[wasm_bindgen(js_name=setTorsionAtomCount)]
    pub fn set_torsion_atom_count(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_torsion_atom_count(value)
        self.inner
            .set_torsion_atom_count(u32_value(&value, "torsionAtomCount")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=torsionAtomCount)]
    pub fn assign_torsion_atom_count(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "number")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_torsion_atom_count(value)
    }
    #[wasm_bindgen(getter,js_name=onlyShortestPaths)]
    pub fn only_shortest_paths(&self) -> Result<bool, JsValue> {
        // COSMolKit❗✔️: self.inner.only_shortest_paths()
        self.inner.only_shortest_paths().map_err(error)
    }
    #[wasm_bindgen(js_name=setOnlyShortestPaths)]
    pub fn set_only_shortest_paths(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        // COSMolKit❗✔️: .set_only_shortest_paths(value)
        self.inner
            .set_only_shortest_paths(bool_value(&value, "onlyShortestPaths")?)
            .map_err(error)
    }
    #[wasm_bindgen(setter,js_name=onlyShortestPaths)]
    pub fn assign_only_shortest_paths(
        &mut self,
        #[wasm_bindgen(unchecked_param_type = "boolean")] value: JsValue,
    ) -> Result<(), JsValue> {
        self.set_only_shortest_paths(value)
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
    pub fn params(&self) -> Result<TopologicalTorsionParams, JsValue> {
        // COSMolKit❗✔️: .map(|inner| TopologicalTorsionParams { inner })
        self.inner
            .params()
            .map(|inner| TopologicalTorsionParams { inner })
            .map_err(error)
    }
}

use crate::topological_torsion_parameters::LegacyTopologicalTorsionParams;

#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseCountLegacy)]
    pub fn fingerprint_topological_torsion_sparse_count_legacy(
        &self,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_sparse_count_fingerprint(
        self.inner
            .fingerprint_topological_torsion_sparse_count_legacy()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionSparseCountLegacyWithParams)]
    pub fn fingerprint_topological_torsion_sparse_count_legacy_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_sparse_count_fingerprint_with_params(
        self.inner
            .fingerprint_topological_torsion_sparse_count_legacy_with_params(&params.inner)
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionCountLegacy)]
    pub fn fingerprint_topological_torsion_count_legacy(
        &self,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_count_fingerprint(
        self.inner
            .fingerprint_topological_torsion_count_legacy()
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionCountLegacyWithParams)]
    pub fn fingerprint_topological_torsion_count_legacy_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<SparseCountFingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_count_fingerprint_with_params(
        self.inner
            .fingerprint_topological_torsion_count_legacy_with_params(&params.inner)
            .map(|inner| SparseCountFingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionLegacy)]
    pub fn fingerprint_topological_torsion_legacy(&self) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_fingerprint(
        self.inner
            .fingerprint_topological_torsion_legacy()
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fingerprintTopologicalTorsionLegacyWithParams)]
    pub fn fingerprint_topological_torsion_legacy_with_params(
        &self,
        params: &LegacyTopologicalTorsionParams,
    ) -> Result<Fingerprint, JsValue> {
        // COSMolKit❗✔️: .legacy_topological_torsion_fingerprint_with_params(
        self.inner
            .fingerprint_topological_torsion_legacy_with_params(&params.inner)
            .map(|inner| Fingerprint { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=topologicalTorsionIds,unchecked_return_type="bigint[]")]
    pub fn topological_torsion_ids(&self) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_ids(
        self.inner
            .topological_torsion_ids()
            .map(|values| values.into_iter().map(JsValue::from).collect())
            .map_err(error)
    }
    #[wasm_bindgen(js_name=topologicalTorsionIdsWithParams,unchecked_return_type="bigint[]")]
    pub fn topological_torsion_ids_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] torsion_atom_count: JsValue,
    ) -> Result<Array, JsValue> {
        // COSMolKit❗✔️: .topological_torsion_ids_with_params(
        self.inner
            .topological_torsion_ids_with_params(u32_value(
                &torsion_atom_count,
                "torsionAtomCount",
            )?)
            .map(|values| values.into_iter().map(JsValue::from).collect())
            .map_err(error)
    }
}
