//! Synchronous query/callback transport over the sole public batch owner.
use crate::batch_boundary::{MoleculeBatch, batch_validation_error};
use crate::host_values::{bool_value, type_error, u32_value, usize_value};
use crate::smiles_parameters::SmilesWriteParams;
use cosmolkit_wasm::rust as ck;
use js_sys::{Array, Function};
use std::sync::{Arc, Mutex};
use wasm_bindgen::prelude::*;
#[wasm_bindgen]
pub struct BatchQueryParams {
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    progress_callback: Option<Function>,
}
#[wasm_bindgen]
impl BatchQueryParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] n_jobs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] progress_bar: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "(() => void) | null")]
        progress_callback: JsValue,
    ) -> Result<Self, JsValue> {
        Ok(Self {
            n_jobs: if n_jobs.is_null() || n_jobs.is_undefined() {
                None
            } else {
                Some(usize_value(&n_jobs, "nJobs")?)
            },
            progress_bar: if progress_bar.is_null() || progress_bar.is_undefined() {
                None
            } else {
                Some(bool_value(&progress_bar, "progressBar")?)
            },
            progress_callback: if progress_callback.is_null() || progress_callback.is_undefined() {
                None
            } else {
                Some(
                    progress_callback
                        .dyn_into::<Function>()
                        .map_err(|_| type_error("progressCallback"))?,
                )
            },
        })
    }
    #[wasm_bindgen(getter,js_name=nJobs,unchecked_return_type="number | null")]
    pub fn n_jobs(&self) -> JsValue {
        self.n_jobs
            .map_or(JsValue::NULL, |n| JsValue::from_f64(n as f64))
    }
    #[wasm_bindgen(getter,js_name=progressBar,unchecked_return_type="boolean | null")]
    pub fn progress_bar(&self) -> JsValue {
        self.progress_bar.map_or(JsValue::NULL, JsValue::from_bool)
    }
    #[wasm_bindgen(getter,js_name=progressCallback,unchecked_return_type="(() => void) | null")]
    pub fn progress_callback(&self) -> JsValue {
        self.progress_callback
            .as_ref()
            .map_or(JsValue::NULL, |f| f.clone().into())
    }
}
impl BatchQueryParams {
    pub(crate) fn execute<T>(
        &self,
        call: impl FnOnce(&ck::BatchQueryParams) -> Result<T, ck::BatchValidationError>,
    ) -> Result<T, JsValue> {
        // COSMolKit❗✔️: canonical_batch_params.rs::execute, boundary behavior:
        //         let callback_error = Arc::new(Mutex::new(None));
        //                     if let Err(source) = callback.call0(py) {
        //                         let mut first = error.lock().expect("callback error mutex");
        //                         if first.is_none() {
        //                             *first = Some(source);
        //                         }
        //                     }
        // Retain the first native JS exception while canonical work completes.
        // wasm-bindgen supplies Send/Sync for JsValue on non-atomic wasm targets;
        // no adapter unsafe implementations or thread-pool substitutions exist.
        let failure: Arc<Mutex<Option<JsValue>>> = Arc::new(Mutex::new(None));
        let progress_callback = self.progress_callback.as_ref().map(|callback| {
            let callback = callback.clone();
            let failure = Arc::clone(&failure);
            Arc::new(move || {
                if let Err(error) = callback.call0(&JsValue::UNDEFINED) {
                    let mut first = failure.lock().expect("callback error mutex");
                    if first.is_none() {
                        *first = Some(error);
                    }
                }
            }) as Arc<dyn Fn() + Send + Sync>
        });
        let execution = ck::BatchQueryParams {
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
            progress_callback,
        };
        // Python gives canonical operation errors precedence when both fail.
        let result =
            call(&execution).map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))?;
        let error = failure.lock().expect("callback error mutex").take();
        match error {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }
}
fn strings(values: Vec<Option<ck::PropertyText>>) -> Result<Array, JsValue> {
    values
        .into_iter()
        .map(|v| crate::host_values::optional_text(v.as_ref()))
        .collect()
}
fn matrices(values: Vec<Option<Vec<Vec<f64>>>>) -> Array {
    values
        .into_iter()
        .map(|matrix| {
            matrix.map_or(JsValue::NULL, |rows| {
                rows.into_iter()
                    .map(|row| row.into_iter().map(JsValue::from_f64).collect::<Array>())
                    .collect::<Array>()
                    .into()
            })
        })
        .collect()
}
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=toSmilesList,unchecked_return_type="(string | null)[]")]
    pub fn to_smiles_list(&self) -> Result<Array, JsValue> {
        self.inner
            .to_smiles_list()
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
            .and_then(strings)
    }
    #[wasm_bindgen(js_name=toSmilesListWithParams,unchecked_return_type="(string | null)[]")]
    pub fn to_smiles_list_with_params(
        &self,
        options: &SmilesWriteParams,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| self.inner.to_smiles_list_with_params(&options.inner, p))
            .and_then(strings)
    }
    #[cfg(feature = "cap-conformer")]
    #[wasm_bindgen(js_name=dgBoundsMatrixList,unchecked_return_type="(number[][] | null)[]")]
    pub fn dg_bounds_matrix_list(&self) -> Result<Array, JsValue> {
        self.inner
            .dg_bounds_matrix_list()
            .map(matrices)
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-conformer")]
    #[wasm_bindgen(js_name=dgBoundsMatrixListWithParams,unchecked_return_type="(number[][] | null)[]")]
    pub fn dg_bounds_matrix_list_with_params(
        &self,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        params
            .execute(|p| self.inner.dg_bounds_matrix_list_with_params(p))
            .map(matrices)
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=toSvgList,unchecked_return_type="(string | null)[]")]
    pub fn to_svg_list(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
    ) -> Result<Array, JsValue> {
        self.inner
            .to_svg_list(u32_value(&width, "width")?, u32_value(&height, "height")?)
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| v.map_or(JsValue::NULL, JsValue::from))
                    .collect()
            })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[cfg(feature = "cap-depict")]
    #[wasm_bindgen(js_name=toSvgListWithParams,unchecked_return_type="(string | null)[]")]
    pub fn to_svg_list_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number")] height: JsValue,
        params: &BatchQueryParams,
    ) -> Result<Array, JsValue> {
        let width = u32_value(&width, "width")?;
        let height = u32_value(&height, "height")?;
        params
            .execute(|p| self.inner.to_svg_list_with_params(width, height, p))
            .map(|values| {
                values
                    .into_iter()
                    .map(|v| v.map_or(JsValue::NULL, JsValue::from))
                    .collect()
            })
    }
}
