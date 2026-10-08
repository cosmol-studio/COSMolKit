//! Full indexed and forward supplier API, with exact u64 bigint metadata.
use crate::alignment_values::set;
use crate::host_values::{type_error, usize_value};
use crate::io_errors::{io_error, sdf_error};
use crate::sdf_reading::{SdfReadParams, SdfRecord};
use cosmolkit_wasm::rust as ck;
use wasm_bindgen::prelude::*;
#[wasm_bindgen(
    inline_js = "export function visitSupplierParams(v,f){try{f(v);}catch(cause){throw new TypeError('invalid SdfReadParams',{cause});}} export function attachSdfIterator(v){Object.defineProperty(v,Symbol.iterator,{value:function(){return this;}});return v;}"
)]
extern "C" {
    #[wasm_bindgen(catch,js_name=visitSupplierParams)]
    fn visit_params(v: &JsValue, f: &mut dyn FnMut(&SdfReadParams)) -> Result<(), JsValue>;
    #[wasm_bindgen(js_name=attachSdfIterator)]
    fn attach_iterator(v: JsValue) -> JsValue;
}
#[wasm_bindgen]
pub struct SdfDataset {
    pub(crate) inner: cosmolkit_wasm::SdfDataset,
}
#[wasm_bindgen]
impl SdfDataset {
    pub fn open(path: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfDataset::open(path)
        cosmolkit_wasm::SdfDataset::open(path)
            .map(|inner| Self { inner })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=openWithParams)]
    pub fn open_with_params(
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfDataset::open_with_params(path,&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::SdfDataset::open_with_params(path, &p.inner)
                    .map(|inner| Self { inner })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
pub struct SdfRecordStream {
    pub(crate) inner: cosmolkit_wasm::SdfRecordStream,
}
#[wasm_bindgen]
impl SdfRecordStream {
    pub fn open(path: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfRecordStream::open(path)
        cosmolkit_wasm::SdfRecordStream::open(path)
            .map(|inner| Self { inner })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=openWithParams)]
    pub fn open_with_params(
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfRecordStream::open_with_params(path,&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::SdfRecordStream::open_with_params(path, &p.inner)
                    .map(|inner| Self { inner })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
pub struct SdfReader {
    pub(crate) inner: cosmolkit_wasm::SdfReader,
}
#[wasm_bindgen]
impl SdfReader {
    pub fn open(path: &str) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfReader::open(path)
        cosmolkit_wasm::SdfReader::open(path)
            .map(|inner| Self { inner })
            .map_err(|e| io_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=openWithParams)]
    pub fn open_with_params(
        path: &str,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<Self, JsValue> {
        // COSMolKit❗✔️: cosmolkit_wasm::SdfReader::open_with_params(path,&p.inner)
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                cosmolkit_wasm::SdfReader::open_with_params(path, &p.inner)
                    .map(|inner| Self { inner })
                    .map_err(|e| io_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
}
#[wasm_bindgen]
pub struct SdfRecordMetadata {
    inner: ck::SdfRecordMetadata,
}
#[wasm_bindgen]
impl SdfRecordMetadata {
    pub fn index(&self) -> usize {
        self.inner.index()
    }
    #[wasm_bindgen(js_name=byteOffset)]
    pub fn byte_offset(&self) -> u64 {
        self.inner.byte_offset()
    }
    #[wasm_bindgen(js_name=byteLen)]
    pub fn byte_len(&self) -> u64 {
        self.inner.byte_len()
    }
    #[wasm_bindgen(js_name=byteRange,unchecked_return_type="[bigint, bigint]")]
    pub fn byte_range(&self) -> JsValue {
        let (a, b) = self.inner.byte_range();
        let out = js_sys::Array::new();
        out.push(&JsValue::from(a));
        out.push(&JsValue::from(b));
        out.into()
    }
    #[wasm_bindgen(js_name=lineRange,unchecked_return_type="[number, number]")]
    pub fn line_range(&self) -> JsValue {
        let (a, b) = self.inner.line_range();
        let out = js_sys::Array::new();
        out.push(&JsValue::from(a as u32));
        out.push(&JsValue::from(b as u32));
        out.into()
    }
    #[wasm_bindgen(unchecked_return_type = "string | null")]
    pub fn title(&self) -> JsValue {
        self.inner.title().map_or(JsValue::NULL, JsValue::from)
    }
}
#[wasm_bindgen]
pub struct SdfDatasetIterator {
    inner: cosmolkit_wasm::SdfDatasetIterator,
}
#[wasm_bindgen]
impl SdfDatasetIterator {
    #[wasm_bindgen(unchecked_return_type = "IteratorResult<SdfRecord>")]
    pub fn next(&mut self) -> Result<JsValue, JsValue> {
        // COSMolKit❗✔️: self.inner.next().transpose()
        let next = self
            .inner
            .next()
            .transpose()
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))?;
        let result = js_sys::Object::new();
        set(&result, "done", JsValue::from(next.is_none()))?;
        set(
            &result,
            "value",
            next.map_or(JsValue::UNDEFINED, |inner| SdfRecord { inner }.into()),
        )?;
        Ok(result.into())
    }
}
#[wasm_bindgen]
impl SdfDataset {
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    #[wasm_bindgen(js_name=isEmpty)]
    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    pub fn path(&self) -> String {
        self.inner.path().to_string_lossy().into_owned()
    }
    #[wasm_bindgen(unchecked_return_type = "SdfRecordMetadata | null")]
    pub fn metadata(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<JsValue, JsValue> {
        Ok(self
            .inner
            .metadata(usize_value(&index, "index")?)
            .map_or(JsValue::NULL, |inner| SdfRecordMetadata { inner }.into()))
    }
    pub fn record(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<SdfRecord, JsValue> {
        self.inner
            .record(usize_value(&index, "index")?)
            .map(|inner| SdfRecord { inner })
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=recordWithParams)]
    pub fn record_with_params(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
        #[wasm_bindgen(unchecked_param_type = "SdfReadParams")] params: JsValue,
    ) -> Result<SdfRecord, JsValue> {
        let index = usize_value(&index, "index")?;
        let mut result = None;
        visit_params(&params, &mut |p: &SdfReadParams| {
            result = Some(
                self.inner
                    .record_with_params(index, &p.inner)
                    .map(|inner| SdfRecord { inner })
                    .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e)),
            );
        })?;
        result.ok_or_else(|| type_error("params"))?
    }
    #[wasm_bindgen(js_name=recordText)]
    pub fn record_text(
        &self,
        #[wasm_bindgen(unchecked_param_type = "number")] index: JsValue,
    ) -> Result<String, JsValue> {
        self.inner
            .record_text(usize_value(&index, "index")?)
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(unchecked_return_type = "SdfDatasetIterator & IterableIterator<SdfRecord>")]
    pub fn iter(&self) -> JsValue {
        attach_iterator(
            SdfDatasetIterator {
                inner: self.inner.iter(),
            }
            .into(),
        )
    }
}
#[wasm_bindgen]
impl SdfRecordStream {
    #[wasm_bindgen(js_name=nextRecord,unchecked_return_type="SdfRecord | null")]
    pub fn next_record(&mut self) -> Result<JsValue, JsValue> {
        self.inner
            .next_record()
            .map(|r| r.map_or(JsValue::NULL, |inner| SdfRecord { inner }.into()))
            .map_err(|e| sdf_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=isEnd)]
    pub fn is_end(&self) -> bool {
        self.inner.is_end()
    }
    #[wasm_bindgen(js_name=recordsConsumed)]
    pub fn records_consumed(&self) -> usize {
        self.inner.records_consumed()
    }
    #[wasm_bindgen(js_name=bytesConsumed)]
    pub fn bytes_consumed(&self) -> u64 {
        self.inner.bytes_consumed()
    }
    #[wasm_bindgen(js_name=linesConsumed)]
    pub fn lines_consumed(&self) -> usize {
        self.inner.lines_consumed()
    }
}
#[wasm_bindgen]
impl SdfReader {
    pub fn path(&self) -> String {
        self.inner.path().to_string_lossy().into_owned()
    }
    pub fn params(&self) -> SdfReadParams {
        SdfReadParams {
            inner: self.inner.params(),
        }
    }
}

#[wasm_bindgen(typescript_custom_section)]
const ITERATOR: &str = "export interface SdfDatasetIterator extends IterableIterator<SdfRecord> {}";
