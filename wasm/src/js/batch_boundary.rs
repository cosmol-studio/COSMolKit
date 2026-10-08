//! Ordered batch transport over the sole public facade batch owner.
use crate::Molecule;
use crate::alignment_values::{set, source_error};
use crate::host_values::{bool_value, sequence, type_error, usize_value};
use crate::smiles_parameters::SmilesParseParams;
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use std::error::Error as RustError;
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Clone, Copy)]
pub enum BatchErrorMode {
    Strict,
    KeepErrors,
}
impl BatchErrorMode {
    fn core(self) -> ck::BatchErrorMode {
        match self {
            Self::Strict => ck::BatchErrorMode::Strict,
            Self::KeepErrors => ck::BatchErrorMode::KeepErrors,
        }
    }
    fn from_core(value: ck::BatchErrorMode) -> Self {
        match value {
            ck::BatchErrorMode::Strict => Self::Strict,
            ck::BatchErrorMode::KeepErrors => Self::KeepErrors,
        }
    }
}

fn optional_usize(value: &JsValue, name: &str) -> Result<Option<usize>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        usize_value(value, name).map(Some)
    }
}
fn optional_bool(value: &JsValue, name: &str) -> Result<Option<bool>, JsValue> {
    if value.is_null() || value.is_undefined() {
        Ok(None)
    } else {
        bool_value(value, name).map(Some)
    }
}

#[wasm_bindgen]
pub struct BatchParams {
    pub(crate) inner: ck::BatchParams,
}
#[wasm_bindgen]
impl BatchParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "BatchErrorMode | null")] errors: Option<
            BatchErrorMode,
        >,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] n_jobs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] progress_bar: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: ck::BatchParams {
                errors: errors
                    .map_or_else(|| ck::BatchParams::default().errors, BatchErrorMode::core),
                n_jobs: optional_usize(&n_jobs, "nJobs")?,
                progress_bar: optional_bool(&progress_bar, "progressBar")?,
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn errors(&self) -> BatchErrorMode {
        BatchErrorMode::from_core(self.inner.errors)
    }
    #[wasm_bindgen(getter, js_name = nJobs, unchecked_return_type = "number | null")]
    pub fn n_jobs(&self) -> JsValue {
        self.inner
            .n_jobs
            .map_or(JsValue::NULL, |n| JsValue::from_f64(n as f64))
    }
    #[wasm_bindgen(getter, js_name = progressBar, unchecked_return_type = "boolean | null")]
    pub fn progress_bar(&self) -> JsValue {
        self.inner
            .progress_bar
            .map_or(JsValue::NULL, JsValue::from_bool)
    }
}

#[wasm_bindgen]
pub struct BatchError {
    pub(crate) inner: ck::BatchError,
}
#[wasm_bindgen]
impl BatchError {
    pub fn index(&self) -> usize {
        self.inner.index
    }
    pub fn operation(&self) -> String {
        self.inner.operation.into()
    }
    pub fn message(&self) -> String {
        self.inner.message.clone()
    }
    #[wasm_bindgen(unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner.source().map_or(Ok(JsValue::NULL), source_error)
    }
    #[wasm_bindgen(js_name = asDict, unchecked_return_type = "[string, string][]")]
    pub fn as_dict(&self) -> Array {
        [
            ("index", self.inner.index.to_string()),
            ("operation", self.operation()),
            ("message", self.message()),
        ]
        .into_iter()
        .map(|(key, value)| {
            let pair = Array::new();
            pair.push(&key.into());
            pair.push(&value.into());
            pair
        })
        .collect()
    }
}

#[wasm_bindgen]
pub struct BatchValidationError {
    inner: ck::BatchValidationError,
}
#[wasm_bindgen]
impl BatchValidationError {
    #[wasm_bindgen(getter)]
    pub fn domain(&self) -> String {
        "batch".into()
    }
    #[wasm_bindgen(getter)]
    pub fn kind(&self) -> String {
        "Validation".into()
    }
    #[wasm_bindgen(getter)]
    pub fn message(&self) -> String {
        self.inner.to_string()
    }
    #[wasm_bindgen(getter)]
    pub fn errors(&self) -> usize {
        self.inner.errors
    }
    #[wasm_bindgen(getter, js_name = recordErrors)]
    pub fn record_errors(&self) -> Vec<BatchError> {
        self.inner
            .record_errors
            .iter()
            .cloned()
            .map(|inner| BatchError { inner })
            .collect()
    }
    #[wasm_bindgen(getter, unchecked_return_type = "never | null")]
    pub fn reason(&self) -> JsValue {
        JsValue::NULL
    }
    #[wasm_bindgen(getter, unchecked_return_type = "Error | null")]
    pub fn cause(&self) -> Result<JsValue, JsValue> {
        self.inner
            .record_errors
            .first()
            .map_or(Ok(JsValue::NULL), batch_record_error)
    }
}

pub(crate) fn batch_record_error(source: &ck::BatchError) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BatchError");
    let error: JsValue = error.into();
    set(&error, "domain", "batch".into())?;
    set(&error, "kind", "Record".into())?;
    set(&error, "index", JsValue::from_f64(source.index as f64))?;
    set(&error, "operation", source.operation.into())?;
    set(
        &error,
        "detail",
        BatchError {
            inner: source.clone(),
        }
        .into(),
    )?;
    if let Some(cause) = source.source() {
        set(&error, "cause", source_error(cause)?)?;
    }
    Ok(error)
}
pub(crate) fn batch_validation_error(
    source: &ck::BatchValidationError,
) -> Result<JsValue, JsValue> {
    let error = js_sys::Error::new(&source.to_string());
    error.set_name("BatchValidationError");
    let error: JsValue = error.into();
    set(&error, "domain", "batch".into())?;
    set(&error, "kind", "Validation".into())?;
    set(&error, "errors", JsValue::from_f64(source.errors as f64))?;
    set(&error, "reason", JsValue::NULL)?;
    let records: Array = source
        .record_errors
        .iter()
        .map(|e| JsValue::from(BatchError { inner: e.clone() }))
        .collect();
    set(&error, "recordErrors", records.into())?;
    set(
        &error,
        "detail",
        BatchValidationError {
            inner: source.clone(),
        }
        .into(),
    )?;
    if let Some(cause) = source.record_errors.first() {
        set(&error, "cause", batch_record_error(cause)?)?;
    }
    Ok(error)
}

#[wasm_bindgen]
pub struct BatchRecord {
    pub(crate) inner: cosmolkit_wasm::BatchRecord,
}
#[wasm_bindgen]
impl BatchRecord {
    pub fn molecule(value: &Molecule) -> Self {
        Self {
            inner: cosmolkit_wasm::BatchRecord::molecule(&value.inner),
        }
    }
    pub fn error(value: &BatchError) -> Self {
        Self {
            inner: cosmolkit_wasm::BatchRecord::error(&value.inner),
        }
    }
    #[wasm_bindgen(js_name = moleculeValue, unchecked_return_type = "Molecule | null")]
    pub fn molecule_value(&self) -> JsValue {
        self.inner.molecule_value().map_or(JsValue::NULL, |inner| {
            Molecule {
                inner: Arc::new(inner),
            }
            .into()
        })
    }
    #[wasm_bindgen(js_name = errorValue, unchecked_return_type = "BatchError | null")]
    pub fn error_value(&self) -> JsValue {
        self.inner
            .error_value()
            .map_or(JsValue::NULL, |inner| BatchError { inner }.into())
    }
}

// The synchronous imported iterator uses wasm-bindgen's RefFromWasmAbi
// callback support, borrowing each class for one callback. It neither reads
// raw pointers nor consumes user objects to build the owned Rust vector.
#[wasm_bindgen(
    inline_js = "export function visit_batch_records(records, visit) { for (const record of records) visit(record); }"
)]
extern "C" {
    #[wasm_bindgen(catch)]
    fn visit_batch_records(
        records: &Array,
        visit: &mut dyn FnMut(&BatchRecord),
    ) -> Result<(), JsValue>;
}

#[wasm_bindgen]
pub struct MoleculeBatch {
    pub(crate) inner: cosmolkit_wasm::MoleculeBatch,
}
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name = fromRecords)]
    pub fn from_records(
        #[wasm_bindgen(unchecked_param_type = "BatchRecord[]")] records: JsValue,
        mode: BatchErrorMode,
    ) -> Result<Self, JsValue> {
        if !Array::is_array(&records) {
            return Err(type_error("records"));
        }
        let records: &Array = records.unchecked_ref();
        let mut copied = Vec::with_capacity(records.length() as usize);
        visit_batch_records(records, &mut |record: &BatchRecord| {
            copied.push(record.inner.clone())
        })?;
        cosmolkit_wasm::MoleculeBatch::from_records(copied, mode.core())
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = fromSmilesList)]
    pub fn from_smiles_list(
        #[wasm_bindgen(unchecked_param_type = "string[]")] smiles: JsValue,
    ) -> Result<Self, JsValue> {
        let values = strings(&smiles)?;
        cosmolkit_wasm::MoleculeBatch::from_smiles_list(&values)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = fromSmilesListWithParams)]
    pub fn from_smiles_list_with_params(
        #[wasm_bindgen(unchecked_param_type = "string[]")] smiles: JsValue,
        parse: &SmilesParseParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        let values = strings(&smiles)?;
        cosmolkit_wasm::MoleculeBatch::from_smiles_list_with_params(
            &values,
            &parse.inner,
            &params.inner,
        )
        .map(|inner| Self { inner })
        .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    #[wasm_bindgen(js_name = isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    #[wasm_bindgen(js_name = validMask, unchecked_return_type = "boolean[]")]
    pub fn valid_mask(&self) -> Array {
        self.inner
            .valid_mask()
            .into_iter()
            .map(JsValue::from_bool)
            .collect()
    }
    #[wasm_bindgen(js_name = invalidMask, unchecked_return_type = "boolean[]")]
    pub fn invalid_mask(&self) -> Array {
        self.inner
            .invalid_mask()
            .into_iter()
            .map(JsValue::from_bool)
            .collect()
    }
    #[wasm_bindgen(js_name = validCount)]
    pub fn valid_count(&self) -> usize {
        self.inner.valid_count()
    }
    #[wasm_bindgen(js_name = invalidCount)]
    pub fn invalid_count(&self) -> usize {
        self.inner.invalid_count()
    }
    pub fn errors(&self) -> Vec<BatchError> {
        self.inner
            .errors()
            .into_iter()
            .map(|inner| BatchError { inner })
            .collect()
    }
    #[wasm_bindgen(js_name = parallelJobs, unchecked_return_type = "number | null")]
    pub fn parallel_jobs(&self) -> JsValue {
        self.inner
            .parallel_jobs()
            .map_or(JsValue::NULL, |n| JsValue::from_f64(n as f64))
    }
    #[wasm_bindgen(js_name = progressBar, unchecked_return_type = "boolean | null")]
    pub fn progress_bar(&self) -> JsValue {
        self.inner
            .progress_bar()
            .map_or(JsValue::NULL, JsValue::from_bool)
    }
    #[wasm_bindgen(js_name = toList, unchecked_return_type = "(Molecule | null)[]")]
    pub fn to_list(&self) -> Array {
        self.inner
            .to_list()
            .into_iter()
            .map(|v| {
                v.map_or(JsValue::NULL, |inner| {
                    Molecule {
                        inner: Arc::new(inner),
                    }
                    .into()
                })
            })
            .collect()
    }
    #[wasm_bindgen(js_name = withValidRecords)]
    pub fn with_valid_records(&self) -> Self {
        Self {
            inner: self.inner.with_valid_records(),
        }
    }
    #[wasm_bindgen(js_name = withParallelJobs)]
    pub fn with_parallel_jobs(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number | null")] n_jobs: JsValue,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_parallel_jobs(optional_usize(&n_jobs, "nJobs")?)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withProgressBar)]
    pub fn with_progress_bar(
        &self,
        #[wasm_bindgen(unchecked_param_type = "boolean | null")] progress_bar: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            inner: self
                .inner
                .with_progress_bar(optional_bool(&progress_bar, "progressBar")?),
        })
    }
}
fn strings(value: &JsValue) -> Result<Vec<String>, JsValue> {
    sequence(value, "smiles")?
        .iter()
        .map(|v| v.as_string().ok_or_else(|| type_error("SMILES element")))
        .collect()
}

#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name = sanitize)]
    pub fn sanitize(&self) -> Result<Self, JsValue> {
        self.inner
            .sanitize()
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = sanitizeWithParams)]
    pub fn sanitize_with_params(
        &self,
        options: &crate::transform_parameters::SanitizeParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .sanitize_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withHydrogens)]
    pub fn with_hydrogens(&self) -> Result<Self, JsValue> {
        self.inner
            .with_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withHydrogensWithParams)]
    pub fn with_hydrogens_with_params(
        &self,
        options: &crate::transform_parameters::AddHsParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_hydrogens_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withoutHydrogens)]
    pub fn without_hydrogens(&self) -> Result<Self, JsValue> {
        self.inner
            .without_hydrogens()
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withoutHydrogensWithParams)]
    pub fn without_hydrogens_with_params(
        &self,
        options: &crate::transform_parameters::RemoveHsParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .without_hydrogens_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withKekulizedBonds)]
    pub fn with_kekulized_bonds(&self) -> Result<Self, JsValue> {
        self.inner
            .with_kekulized_bonds()
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name = withKekulizedBondsWithParams)]
    pub fn with_kekulized_bonds_with_params(
        &self,
        options: &crate::transform_parameters::KekulizeParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_kekulized_bonds_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name = with2dCoordinates)]
    pub fn with_2d_coordinates(&self) -> Result<Self, JsValue> {
        self.inner
            .with_2d_coordinates()
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name = with2dCoordinatesWithParams)]
    pub fn with_2d_coordinates_with_params(
        &self,
        options: &crate::transform_parameters::Coordinate2DParams,
        params: &BatchParams,
    ) -> Result<Self, JsValue> {
        self.inner
            .with_2d_coordinates_with_params(&options.inner, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
}
