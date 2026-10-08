//! Typed UFF projection; calculation and optimization stay in cosmolkit.
use crate::Molecule;
use crate::alignment_values::{operation_error, set, source_error};
use crate::host_values::{bool_value, i32_value, type_error};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::{error::Error as RustError, sync::Arc};
use wasm_bindgen::prelude::*;
pub(crate) fn source_conformer_id(value: &JsValue) -> Result<Option<usize>, JsValue> {
    // COSMolKit❗✔️: let conformer_id = conformer_id.and_then(crate::mmff_binding::source_conformer_id);
    if value.is_null() || value.is_undefined() {
        return Ok(None);
    }
    let n = i32_value(value, "conformerId")?;
    Ok((n >= 0).then_some(n as usize))
}
#[wasm_bindgen]
pub struct UffEvaluationParams {
    inner: ck::UffEvaluationParams,
}
#[wasm_bindgen]
impl UffEvaluationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] vdw_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::UffEvaluationParams {
        let mut inner = ck::UffEvaluationParams::default();
        if !vdw_threshold.is_undefined() {
            inner.vdw_threshold = vdw_threshold
                .as_f64()
                .ok_or_else(|| type_error("vdwThreshold"))?;
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
    #[wasm_bindgen(getter,js_name=vdwThreshold)]
    pub fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
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
pub struct UffOptimizationParams {
    inner: ck::UffOptimizationParams,
}
#[wasm_bindgen]
impl UffOptimizationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] vdw_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] conformer_id: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::UffOptimizationParams {
        let mut inner = ck::UffOptimizationParams::default();
        if !max_iterations.is_undefined() {
            inner.max_iterations = i32_value(&max_iterations, "maxIterations")?;
        }
        if !vdw_threshold.is_undefined() {
            inner.vdw_threshold = vdw_threshold
                .as_f64()
                .ok_or_else(|| type_error("vdwThreshold"))?;
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
    #[wasm_bindgen(getter,js_name=maxIterations)]
    pub fn max_iterations(&self) -> i32 {
        self.inner.max_iterations
    }
    #[wasm_bindgen(getter,js_name=vdwThreshold)]
    pub fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
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
pub struct UffConformerOptimizationParams {
    inner: ck::UffConformerOptimizationParams,
}
#[wasm_bindgen]
impl UffConformerOptimizationParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number")] num_threads: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] max_iterations: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] vdw_threshold: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean")]
        ignore_interfragment_interactions: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: inner: ck::UffConformerOptimizationParams {
        let mut inner = ck::UffConformerOptimizationParams::default();
        if !num_threads.is_undefined() {
            inner.num_threads = i32_value(&num_threads, "numThreads")?;
        }
        if !max_iterations.is_undefined() {
            inner.max_iterations = i32_value(&max_iterations, "maxIterations")?;
        }
        if !vdw_threshold.is_undefined() {
            inner.vdw_threshold = vdw_threshold
                .as_f64()
                .ok_or_else(|| type_error("vdwThreshold"))?;
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
    #[wasm_bindgen(getter,js_name=vdwThreshold)]
    pub fn vdw_threshold(&self) -> f64 {
        self.inner.vdw_threshold
    }
    #[wasm_bindgen(getter,js_name=ignoreInterfragmentInteractions)]
    pub fn ignore_interfragment_interactions(&self) -> bool {
        self.inner.ignore_interfragment_interactions
    }
}
#[wasm_bindgen]
pub struct UffEnergyGradient {
    inner: ck::UffEnergyGradient,
}
#[wasm_bindgen]
impl UffEnergyGradient {
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
pub struct UffOptimizationResult {
    inner: cosmolkit_wasm::UffOptimizationResult,
}
#[wasm_bindgen]
impl UffOptimizationResult {
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
    pub fn energy(&self) -> f64 {
        self.inner.energy()
    }
}
#[wasm_bindgen]
pub struct UffConformerResult {
    inner: ck::UffConformerResult,
}
#[wasm_bindgen]
impl UffConformerResult {
    #[wasm_bindgen(js_name=conformerId)]
    pub fn conformer_id(&self) -> usize {
        self.inner.conformer_id()
    }
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
pub struct UffConformerOptimizationResult {
    inner: cosmolkit_wasm::UffConformerOptimizationResult,
}
#[wasm_bindgen]
impl UffConformerOptimizationResult {
    pub fn molecule(&self) -> Molecule {
        Molecule {
            inner: Arc::new(self.inner.molecule()),
        }
    }
    #[wasm_bindgen(js_name=conformerResults,unchecked_return_type="UffConformerResult[]")]
    pub fn conformer_results(&self) -> Array {
        self.inner
            .conformer_results()
            .iter()
            .copied()
            .map(|inner| JsValue::from(UffConformerResult { inner }))
            .collect()
    }
}
#[wasm_bindgen]
pub struct UffOptimizationErrorKind {
    inner: ck::UffOptimizationErrorKind,
}
#[wasm_bindgen]
impl UffOptimizationErrorKind {
    #[wasm_bindgen(getter)]
    pub fn variant(&self) -> String {
        match self.inner {
            ck::UffOptimizationErrorKind::MissingConformer { .. } => "MissingConformer",
            ck::UffOptimizationErrorKind::Rings => "Rings",
            ck::UffOptimizationErrorKind::Optimization => "Optimization",
            ck::UffOptimizationErrorKind::ConformerOptimization => "ConformerOptimization",
            ck::UffOptimizationErrorKind::Evaluation => "Evaluation",
        }
        .into()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "number | null")]
    pub fn requested(&self) -> JsValue {
        match self.inner {
            ck::UffOptimizationErrorKind::MissingConformer { requested } => {
                requested.map_or(JsValue::NULL, |n| JsValue::from(n as u32))
            }
            _ => JsValue::NULL,
        }
    }
    pub fn equals(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
#[wasm_bindgen]
pub struct UffOptimizationError {
    kind: ck::UffOptimizationErrorKind,
    message: String,
    cause: JsValue,
}
#[wasm_bindgen]
impl UffOptimizationError {
    pub fn kind(&self) -> UffOptimizationErrorKind {
        UffOptimizationErrorKind { inner: self.kind }
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.message.clone()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> JsValue {
        self.cause.clone()
    }
}
pub(crate) fn uff_optimization_error(
    source: &ck::UffOptimizationError,
) -> Result<JsValue, JsValue> {
    // COSMolKit❗✔️: .setattr("domain", "uff_optimization")
    let kind = UffOptimizationErrorKind {
        inner: source.kind(),
    };
    let cause = source.source().map_or(Ok(JsValue::NULL), source_error)?;
    let e = js_sys::Error::new(&source.to_string());
    e.set_name("UffOptimizationError");
    let e: JsValue = e.into();
    set(&e, "domain", "uff_optimization".into())?;
    set(&e, "kind", kind.variant().into())?;
    set(&e, "requested", kind.requested())?;
    set(&e, "_kind", kind.into())?;
    if !cause.is_null() {
        set(&e, "cause", cause.clone())?;
    }
    set(
        &e,
        "detail",
        UffOptimizationError {
            kind: source.kind(),
            message: source.to_string(),
            cause,
        }
        .into(),
    )?;
    Ok(e)
}
#[wasm_bindgen]
impl Molecule {
    #[wasm_bindgen(js_name=uffEnergyGradient)]
    pub fn uff_energy_gradient(&self) -> Result<UffEnergyGradient, JsValue> {
        // COSMolKit❗✔️: self.inner.uff_energy_gradient()
        self.inner
            .uff_energy_gradient()
            .map(|inner| UffEnergyGradient { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=uffEnergyGradientWithParams)]
    pub fn uff_energy_gradient_with_params(
        &self,
        params: &UffEvaluationParams,
    ) -> Result<UffEnergyGradient, JsValue> {
        // COSMolKit❗✔️: self.inner.uff_energy_gradient_with_params(&params.inner)
        self.inner
            .uff_energy_gradient_with_params(&params.inner)
            .map(|inner| UffEnergyGradient { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withUffOptimized)]
    pub fn with_uff_optimized(&self) -> Result<UffOptimizationResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized()
        self.inner
            .with_uff_optimized()
            .map(|inner| UffOptimizationResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withUffOptimizedWithParams)]
    pub fn with_uff_optimized_with_params(
        &self,
        params: &UffOptimizationParams,
    ) -> Result<UffOptimizationResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_with_params(&params.inner)
        self.inner
            .with_uff_optimized_with_params(&params.inner)
            .map(|inner| UffOptimizationResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withUffOptimizedConfs)]
    pub fn with_uff_optimized_confs(&self) -> Result<UffConformerOptimizationResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_confs()
        self.inner
            .with_uff_optimized_confs()
            .map(|inner| UffConformerOptimizationResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=withUffOptimizedConfsWithParams)]
    pub fn with_uff_optimized_confs_with_params(
        &self,
        params: &UffConformerOptimizationParams,
    ) -> Result<UffConformerOptimizationResult, JsValue> {
        // COSMolKit❗✔️: self.inner.with_uff_optimized_confs_with_params(&params.inner)
        self.inner
            .with_uff_optimized_confs_with_params(&params.inner)
            .map(|inner| UffConformerOptimizationResult { inner })
            .map_err(|e| operation_error(&e).unwrap_or_else(|e| e))
    }
}
