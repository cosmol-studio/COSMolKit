//! Full MMFF values and errors; all chemistry delegates to cosmolkit.
use crate::Molecule;
use crate::alignment_values::{operation_error, set, source_error};
use crate::host_values::{bool_value, i32_value, type_error};
use crate::uff::source_conformer_id;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::{error::Error as RustError, sync::Arc};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct MmffEvaluationParams {
    inner: ck::MmffEvaluationParams,
}
#[wasm_bindgen]
impl MmffEvaluationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "string")] mmff_variant: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] non_bonded_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::MmffEvaluationParams {
        let mut inner = ck::MmffEvaluationParams::default();
        if !mmff_variant.is_undefined() {
            inner.mmff_variant = mmff_variant
                .as_string()
                .ok_or_else(|| type_error("mmffVariant"))?;
        }
        if !non_bonded_threshold.is_undefined() {
            inner.non_bonded_threshold = non_bonded_threshold
                .as_f64()
                .ok_or_else(|| type_error("nonBondedThreshold"))?;
        }
        if !conformer_id.is_undefined() {
            inner.conformer_id = source_conformer_id(&conformer_id)?;
        }
        if !ignore_interfragment_interactions.is_undefined() {
            inner.ignore_interfragment_interactions = bool_value(
                &ignore_interfragment_interactions,
                "ignoreInterfragmentInteractions",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=mmffVariant)]
    pub fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[wasm_bindgen(getter,js_name=nonBondedThreshold)]
    pub fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id
            .map_or(JsValue::NULL, |v| JsValue::from(v as u32))
    }
    #[wasm_bindgen(getter,js_name=ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[wasm_bindgen]
pub struct MmffOptimizationParams {
    inner: ck::MmffOptimizationParams,
}
#[wasm_bindgen]
impl MmffOptimizationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "string")] mmff_variant: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] non_bonded_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::MmffOptimizationParams {
        let mut inner = ck::MmffOptimizationParams::default();
        if !mmff_variant.is_undefined() {
            inner.mmff_variant = mmff_variant
                .as_string()
                .ok_or_else(|| type_error("mmffVariant"))?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = i32_value(&max_iterations, "maxIterations")?;
        }
        if !non_bonded_threshold.is_undefined() {
            inner.non_bonded_threshold = non_bonded_threshold
                .as_f64()
                .ok_or_else(|| type_error("nonBondedThreshold"))?;
        }
        if !conformer_id.is_undefined() {
            inner.conformer_id = source_conformer_id(&conformer_id)?;
        }
        if !ignore_interfragment_interactions.is_undefined() {
            inner.ignore_interfragment_interactions = bool_value(
                &ignore_interfragment_interactions,
                "ignoreInterfragmentInteractions",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=mmffVariant)]
    pub fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[wasm_bindgen(getter,js_name=maxIterations)]
    pub fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(getter,js_name=nonBondedThreshold)]
    pub fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
    }
    #[wasm_bindgen(getter,js_name=conformerId,unchecked_return_type="number | null")]
    pub fn conformer_id(&self) -> JsValue {
        self.inner
            .conformer_id
            .map_or(JsValue::NULL, |v| JsValue::from(v as u32))
    }
    #[wasm_bindgen(getter,js_name=ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[wasm_bindgen]
pub struct MmffConformerOptimizationParams {
    inner: ck::MmffConformerOptimizationParams,
}
#[wasm_bindgen]
impl MmffConformerOptimizationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string")] mmff_variant: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] non_bonded_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::MmffConformerOptimizationParams {
        let mut inner = ck::MmffConformerOptimizationParams::default();
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = i32_value(&max_iterations, "maxIterations")?;
        }
        if !mmff_variant.is_undefined() {
            inner.mmff_variant = mmff_variant
                .as_string()
                .ok_or_else(|| type_error("mmffVariant"))?;
        }
        if !non_bonded_threshold.is_undefined() {
            inner.non_bonded_threshold = non_bonded_threshold
                .as_f64()
                .ok_or_else(|| type_error("nonBondedThreshold"))?;
        }
        if !ignore_interfragment_interactions.is_undefined() {
            inner.ignore_interfragment_interactions = bool_value(
                &ignore_interfragment_interactions,
                "ignoreInterfragmentInteractions",
            )?;
        }
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter,js_name=numThreads)]
    pub fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[wasm_bindgen(getter,js_name=maxIterations)]
    pub fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(getter,js_name=mmffVariant)]
    pub fn mmff_variant(&self) -> String {
        self.inner.mmff_variant.clone()
    }
    #[wasm_bindgen(getter,js_name=nonBondedThreshold)]
    pub fn non_bonded_threshold(&self) -> f64 {
        self.inner.non_bonded_threshold
    }
    #[wasm_bindgen(getter,js_name=ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[wasm_bindgen]
pub struct MmffEnergyGradient {
    inner: ck::MmffEnergyGradient,
}
#[wasm_bindgen]
impl MmffEnergyGradient {
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
    #[wasm_bindgen(unchecked_return_type = "number[]")]
    pub fn gradient(&self) -> Array {
        self.inner
            .gradient()
            .iter()
            .map(|v| JsValue::from(*v))
            .collect()
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeResult {
    inner: cosmolkit_wasm::MmffOptimizeMoleculeResult,
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeResult {
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    #[wasm_bindgen(js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[wasm_bindgen(js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeConfResult {
    inner: ck::MmffOptimizeMoleculeConfResult,
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfResult {
    #[wasm_bindgen(js_name=statusCode)]
    pub fn status_code(&self) -> i32 {
        self.inner.status_code()
    }
    #[wasm_bindgen(js_name=needsMore)]
    pub fn needs_more(&self) -> bool {
        self.inner.needs_more()
    }
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
}
#[wasm_bindgen]
pub struct MmffOptimizeMoleculeConfsResult {
    inner: cosmolkit_wasm::MmffOptimizeMoleculeConfsResult,
}
#[wasm_bindgen]
impl MmffOptimizeMoleculeConfsResult {
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    #[wasm_bindgen(js_name=conformerResults,unchecked_return_type="MmffOptimizeMoleculeConfResult[]")]
    pub fn conformer_results(&self) -> Array {
        self.inner
            .conformer_results()
            .iter()
            .copied()
            .map(|inner| JsValue::from(MmffOptimizeMoleculeConfResult { inner }))
            .collect()
    }
}
#[wasm_bindgen]
pub struct MmffOptimizationError {
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl MmffOptimizationError {
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn mmff_optimization_error(
    source: &ck::MmffOptimizationError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: .setattr("domain", "mmff_optimization")
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("MmffOptimizationError");
    let e: JsValue = e.into();
    set(&e, "domain", "mmff_optimization".into())?;
    set(&e, "kind", "MmffOptimization".into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        MmffOptimizationError {
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=mmffEnergyGradient,unchecked_return_type="MmffEnergyGradient | null")]
    pub fn mmff_energy_gradient(&self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.mmff_energy_gradient()
        self.inner
            .mmff_energy_gradient()
            .map(|inner| {
                inner.map_or(JsValue::NULL, |inner| {
                    JsValue::from(MmffEnergyGradient { inner })
                })
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=mmffEnergyGradientWithParams,unchecked_return_type="MmffEnergyGradient | null")]
    pub fn mmff_energy_gradient_with_params(
        &self,
        params: &MmffEvaluationParams,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.mmff_energy_gradient_with_params(&params.inner)
        self.inner
            .mmff_energy_gradient_with_params(&params.inner)
            .map(|inner| {
                inner.map_or(JsValue::NULL, |inner| {
                    JsValue::from(MmffEnergyGradient { inner })
                })
            })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withMmffOptimized)]
    pub fn with_mmff_optimized(&self) -> Result<MmffOptimizeMoleculeResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized()
        self.inner
            .with_mmff_optimized()
            .map(|inner| MmffOptimizeMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withMmffOptimizedWithParams)]
    pub fn with_mmff_optimized_with_params(
        &self,
        params: &MmffOptimizationParams,
    ) -> Result<MmffOptimizeMoleculeResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_with_params(&params.inner)
        self.inner
            .with_mmff_optimized_with_params(&params.inner)
            .map(|inner| MmffOptimizeMoleculeResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withMmffOptimizedConfs)]
    pub fn with_mmff_optimized_confs(&self) -> Result<MmffOptimizeMoleculeConfsResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_confs()
        self.inner
            .with_mmff_optimized_confs()
            .map(|inner| MmffOptimizeMoleculeConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withMmffOptimizedConfsWithParams)]
    pub fn with_mmff_optimized_confs_with_params(
        &self,
        params: &MmffConformerOptimizationParams,
    ) -> Result<MmffOptimizeMoleculeConfsResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_mmff_optimized_confs_with_params(&params.inner)
        self.inner
            .with_mmff_optimized_confs_with_params(&params.inner)
            .map(|inner| MmffOptimizeMoleculeConfsResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
