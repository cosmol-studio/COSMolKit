//! Complete SDF batch parameters, reads, record values, exports and iterators.
use crate::alignment_values::set;
use crate::batch_boundary::{BatchRecord, MoleculeBatch, batch_validation_error};
use crate::batch_images::BatchExportReport;
use crate::host_values::{bool_value, sequence, type_error, u32_value, usize_value};
use crate::sdf_datasets::{SdfDataset, SdfReader, SdfRecordStream};
use crate::sdf_reading::SdfReadParams;
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
fn error(e: ck::BatchValidationError) -> JsValue {
    batch_validation_error(&e).unwrap_or_else(|e| e)
}
fn mode(v: &JsValue) -> Result<ck::BatchErrorMode, JsValue> {
    match u32_value(v, "mode")? {
        0 => Ok(ck::BatchErrorMode::Strict),
        1 => Ok(ck::BatchErrorMode::KeepErrors),
        _ => Err(js_sys::RangeError::new("invalid BatchErrorMode").into()),
    }
}
fn jobs(v: &JsValue) -> Result<Option<usize>, JsValue> {
    if v.is_null() || v.is_undefined() {
        Ok(None)
    } else {
        usize_value(v, "nJobs").map(Some)
    }
}
fn optional_bool(v: &JsValue) -> Result<Option<bool>, JsValue> {
    if v.is_null() || v.is_undefined() {
        Ok(None)
    } else {
        bool_value(v, "progressBar").map(Some)
    }
}
fn optional_string(v: &JsValue) -> Result<Option<String>, JsValue> {
    if v.is_null() || v.is_undefined() {
        Ok(None)
    } else {
        v.as_string()
            .map(Some)
            .ok_or_else(|| type_error("reportPath"))
    }
}
#[wasm_bindgen]
pub struct BatchExportParams {
    inner: ck::BatchExportParams,
}
#[wasm_bindgen]
impl BatchExportParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "SdfFormat")] format: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "BatchErrorMode | null")] errors: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] n_jobs: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "boolean | null")] progress_bar: JsValue,
    ) -> Result<Self, JsValue> {
        let mut inner = ck::BatchExportParams::default();
        if !format.is_undefined() {
            inner.format = match u32_value(&format, "format")? {
                0 => ck::SdfFormat::V2000,
                1 => ck::SdfFormat::V3000,
                _ => return Err(js_sys::RangeError::new("invalid SdfFormat").into()),
            };
        }
        if !errors.is_undefined() && !errors.is_null() {
            inner.errors = Some(mode(&errors)?);
        }
        inner.n_jobs = jobs(&n_jobs)?;
        inner.progress_bar = optional_bool(&progress_bar)?;
        Ok(Self { inner })
    }
    #[wasm_bindgen(getter, unchecked_return_type = "SdfFormat")]
    pub fn format(&self) -> u32 {
        match self.inner.format {
            ck::SdfFormat::V2000 => 0,
            ck::SdfFormat::V3000 => 1,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "BatchErrorMode | null")]
    pub fn errors(&self) -> JsValue {
        self.inner.errors.map_or(JsValue::NULL, |mode| JsValue::from(match mode {
            ck::BatchErrorMode::Strict => 0u32,
            ck::BatchErrorMode::KeepErrors => 1u32,
        }))
    }
    #[wasm_bindgen(getter,js_name=nJobs,unchecked_return_type="number | null")]
    pub fn n_jobs(&self) -> JsValue {
        self.inner
            .n_jobs
            .map_or(JsValue::NULL, |n| JsValue::from(n as u32))
    }
    #[wasm_bindgen(getter,js_name=progressBar,unchecked_return_type="boolean | null")]
    pub fn progress_bar(&self) -> JsValue {
        self.inner.progress_bar.map_or(JsValue::NULL, JsValue::from)
    }
}
#[wasm_bindgen(
    inline_js = "export function visitSdfBatchRead(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SdfReadParams',{cause});}} export function visitSdfBatchExport(v,f){try{f(v);}catch(cause){throw new TypeError('invalid BatchExportParams',{cause});}} export function visitSdfDataset(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SdfDataset',{cause});}} export function attachSdfBatchIterator(v){Object.defineProperty(v,Symbol.iterator,{value:function(){return this;}});return v;}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSdfBatchRead)]
    fn visit_read(v: &JsValue, f: &mut dyn FnMut(&SdfReadParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitSdfBatchExport)]
    fn visit_export(v: &JsValue, f: &mut dyn FnMut(&BatchExportParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(catch,js_name=visitSdfDataset)]
    fn visit_dataset(v: &JsValue, f: &mut dyn FnMut(&SdfDataset)) -> Result<(), JsValue>;
    #[wasm_bindgen(js_name=attachSdfBatchIterator)]
    fn attach_iterator(v: JsValue) -> JsValue;
}
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=fromSdfRecords)]
    pub fn from_sdf_records(text: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::MoleculeBatch::from_sdf_records(text)
        cosmolkit_wasm::MoleculeBatch::from_sdf_records(text)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=readSdf)]
    pub fn read_sdf(path: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::MoleculeBatch::read_sdf(path)
        cosmolkit_wasm::MoleculeBatch::read_sdf(path)
            .map(|inner| Self { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=fromSdfRecordsWithParams)]
    pub fn from_sdf_records_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] read: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BatchErrorMode")] error_mode: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number | null")] n_jobs: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::MoleculeBatch::from_sdf_records_with_params(text,&p.inner,mode,n_jobs)
        let mode = mode(&error_mode)?;
        let n_jobs = jobs(&n_jobs)?;
        let mut result = None;
        visit_read(&read, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::MoleculeBatch::from_sdf_records_with_params(
                    text, &p.inner, mode, n_jobs,
                )
                .map(|inner| Self { inner })
                .map_err(error),
            );
        })?;
        result.ok_or_else(|| type_error("read"))?
    }
    #[wasm_bindgen(js_name=readSdfWithParams)]
    pub fn read_sdf_with_params(
        text: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] read: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BatchErrorMode")] error_mode: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number | null")] n_jobs: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] progress_bar: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::MoleculeBatch::read_sdf_with_params(text,&p.inner,mode,n_jobs,progress_bar)
        let mode = mode(&error_mode)?;
        let n_jobs = jobs(&n_jobs)?;
        let progress_bar = bool_value(&progress_bar, "progressBar")?;
        let mut result = None;
        visit_read(&read, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::MoleculeBatch::read_sdf_with_params(
                    text,
                    &p.inner,
                    mode,
                    n_jobs,
                    progress_bar,
                )
                .map(|inner| Self { inner })
                .map_err(error),
            );
        })?;
        result.ok_or_else(|| type_error("read"))?
    }
    #[wasm_bindgen(js_name=fromDatasetIndices)]
    pub fn from_dataset_indices(
        #[wasm_bindgen(unchecked_param_type = "SdfDataset")] dataset: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number[] | Uint32Array")] indices: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BatchErrorMode")] error_mode: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::MoleculeBatch::from_dataset_indices(&d.inner,&indices,mode)
        let mode = mode(&error_mode)?;
        let indices = sequence(&indices, "indices")?
            .iter()
            .map(|v| usize_value(&v, "index"))
            .collect::<Result<Vec<_>, _>>()?;
        let mut result = None;
        visit_dataset(&dataset, &mut |d: &SdfDataset| {
            result = Some(
                cosmolkit_wasm::MoleculeBatch::from_dataset_indices(&d.inner, &indices, mode)
                    .map(|inner| Self { inner })
                    .map_err(error),
            );
        })?;
        result.ok_or_else(|| type_error("dataset"))?
    }
    #[wasm_bindgen(unchecked_return_type = "BatchRecord | null")]
    pub fn get(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .get(usize_value(&index, "index")?)
            .map_or(JsValue::NULL, |inner| BatchRecord { inner }.into()))
    }
    pub fn records(&self) -> Vec<BatchRecord> {
        self.inner
            .records()
            .into_iter()
            .map(|inner| BatchRecord { inner })
            .collect()
    }
    #[wasm_bindgen(js_name=writeSdf)]
    pub fn write_sdf(&self, path: &str) -> Result<BatchExportReport, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf(path)
        self.inner
            .write_sdf(path)
            .map(|inner| BatchExportReport { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=writeSdfFiles)]
    pub fn write_sdf_files(&self, path: &str) -> Result<BatchExportReport, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_files(path)
        self.inner
            .write_sdf_files(path)
            .map(|inner| BatchExportReport { inner })
            .map_err(error)
    }
    #[wasm_bindgen(js_name=writeSdfWithParams)]
    pub fn write_sdf_with_params(
        &self,
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "BatchExportParams")] params: JsValue,
        #[wasm_bindgen(unchecked_param_type = "string | null")] report_path: JsValue,
    ) -> Result<BatchExportReport, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_with_params(path,&p.inner,report_path.as_deref())
        let report_path = optional_string(&report_path)?;
        let mut result = None;
        visit_export(&params, &mut |p: &BatchExportParams| {
            result = Some(
                self.inner
                    .write_sdf_with_params(path, &p.inner, report_path.as_deref())
                    .map(|inner| BatchExportReport { inner })
                    .map_err(error),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=writeSdfFilesWithParams)]
    pub fn write_sdf_files_with_params(
        &self,
        directory: &str,
        #[wasm_bindgen(unchecked_param_type = "BatchExportParams")] params: JsValue,
        #[wasm_bindgen(unchecked_param_type = "(string | null)[] | null")] filenames: JsValue,
        #[wasm_bindgen(unchecked_param_type = "string | null")] report_path: JsValue,
    ) -> Result<BatchExportReport, JsValue> {
        // COSMolKit❗✔️: self.inner.write_sdf_files_with_params(directory,&p.inner,filenames.as_deref(),report_path.as_deref())
        let report_path = optional_string(&report_path)?;
        let filenames = if filenames.is_null() {
            None
        } else {
            Some(
                sequence(&filenames, "filenames")?
                    .iter()
                    .map(|v| {
                        if v.is_null() {
                            Ok(None)
                        } else {
                            v.as_string()
                                .map(Some)
                                .ok_or_else(|| type_error("filename"))
                        }
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            )
        };
        let mut result = None;
        visit_export(&params, &mut |p: &BatchExportParams| {
            result = Some(
                self.inner
                    .write_sdf_files_with_params(
                        directory,
                        &p.inner,
                        filenames.as_deref(),
                        report_path.as_deref(),
                    )
                    .map(|inner| BatchExportReport { inner })
                    .map_err(error),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
pub struct SdfBatchIterator {
    inner: cosmolkit_wasm::SdfBatchIterator,
}
#[wasm_bindgen]
impl SdfBatchIterator {
    #[wasm_bindgen(js_name=nextBatch,unchecked_return_type="MoleculeBatch | null")]
    pub fn next_batch(&mut self) -> Result<JsValue, JsValue> {
        self.inner
            .next_batch()
            .map(|b| b.map_or(JsValue::NULL, |inner| MoleculeBatch { inner }.into()))
            .map_err(error)
    }
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<MoleculeBatch>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        let b = self.inner.next_batch().map_err(error)?;
        let result = js_sys::Object::new();
        set(&result, "done", JsValue::from(b.is_none()))?;
        set(
            &result,
            "value",
            b.map_or(JsValue::UNDEFINED, |inner| MoleculeBatch { inner }.into()),
        )?;
        Ok(result.into())
    }
}
#[wasm_bindgen]
pub struct SdfReaderBatchIterator {
    inner: cosmolkit_wasm::SdfReaderBatchIterator,
}
#[wasm_bindgen]
impl SdfReaderBatchIterator {
    #[wasm_bindgen(js_name=nextBatch,unchecked_return_type="MoleculeBatch | null")]
    pub fn next_batch(&mut self) -> Result<JsValue, JsValue> {
        self.inner
            .next_batch()
            .map(|b| b.map_or(JsValue::NULL, |inner| MoleculeBatch { inner }.into()))
            .map_err(error)
    }
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<MoleculeBatch>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        let b = self.inner.next_batch().map_err(error)?;
        let result = js_sys::Object::new();
        set(&result, "done", JsValue::from(b.is_none()))?;
        set(
            &result,
            "value",
            b.map_or(JsValue::UNDEFINED, |inner| MoleculeBatch { inner }.into()),
        )?;
        Ok(result.into())
    }
}
#[wasm_bindgen]
impl SdfDataset {
    #[wasm_bindgen(unchecked_return_type = "SdfBatchIterator & IterableIterator<MoleculeBatch>")]
    pub fn batches(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] size: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BatchErrorMode")] error_mode: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number | null")] n_jobs: JsValue,
        #[wasm_bindgen(unchecked_param_type = "boolean")] progress_bar: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs,progress_bar)
        let size = usize_value(&size, "size")?;
        let mode = mode(&error_mode)?;
        let n_jobs = jobs(&n_jobs)?;
        let progress_bar = bool_value(&progress_bar, "progressBar")?;
        self.inner
            .batches(size, mode, n_jobs, progress_bar)
            .map(|inner| attach_iterator(SdfBatchIterator { inner }.into()))
            .map_err(error)
    }
}
#[wasm_bindgen]
impl SdfReader {
    #[wasm_bindgen(
        unchecked_return_type = "SdfReaderBatchIterator & IterableIterator<MoleculeBatch>"
    )]
    pub fn batches(
        &self,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] size: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "BatchErrorMode")] error_mode: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number | null")] n_jobs: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs)
        let size = if size.is_undefined() {
            1024
        } else {
            usize_value(&size, "size")?
        };
        let mode = if error_mode.is_undefined() {
            ck::BatchErrorMode::Strict
        } else {
            mode(&error_mode)?
        };
        self.inner
            .batches(size, mode, jobs(&n_jobs)?)
            .map(|inner| attach_iterator(SdfReaderBatchIterator { inner }.into()))
            .map_err(error)
    }
}
#[wasm_bindgen]
impl SdfRecordStream {
    #[wasm_bindgen(
        unchecked_return_type = "SdfReaderBatchIterator & IterableIterator<MoleculeBatch>"
    )]
    pub fn batches(
        self,
        #[wasm_bindgen(unchecked_param_type = "number")] size: JsValue,
        #[wasm_bindgen(unchecked_param_type = "BatchErrorMode")] error_mode: JsValue,
        #[wasm_bindgen(unchecked_param_type = "number | null")] n_jobs: JsValue,
    ) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.batches(size,mode,n_jobs)
        self.inner
            .batches(
                usize_value(&size, "size")?,
                mode(&error_mode)?,
                jobs(&n_jobs)?,
            )
            .map(|inner| attach_iterator(SdfReaderBatchIterator { inner }.into()))
            .map_err(error)
    }
}
#[wasm_bindgen(typescript_custom_section)]
const ITERATORS: &str = "export interface SdfBatchIterator extends IterableIterator<MoleculeBatch> {}\nexport interface SdfReaderBatchIterator extends IterableIterator<MoleculeBatch> {}";
