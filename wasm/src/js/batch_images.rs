//! Image/report transport; filesystem effects remain in canonical Rust APIs.
use crate::batch_boundary::{BatchError, BatchParams, MoleculeBatch, batch_validation_error};
use crate::host_values::{sequence, type_error, u32_value};
use cosmolkit_wasm::rust as ck;
use js_sys::Array;
use wasm_bindgen::prelude::*;
#[cfg(feature = "cap-depict")]
#[wasm_bindgen]
pub struct BatchImageParams {
    inner: ck::BatchImageParams,
}
#[cfg(feature = "cap-depict")]
#[wasm_bindgen]
impl BatchImageParams {
    #[wasm_bindgen(constructor)]
    pub fn new(
        #[wasm_bindgen(unchecked_optional_param_type = "string")] format: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] width: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "number")] height: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "BatchParams | null")] execution: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "(string | null)[] | null")]
        filenames: JsValue,
        #[wasm_bindgen(unchecked_optional_param_type = "string | null")] report_path: JsValue,
    ) -> Result<Self, JsValue> {
        let defaults = ck::BatchImageParams::default();
        let mut selected = defaults.execution;
        if !execution.is_null() && !execution.is_undefined() {
            visit_execution(&execution, &mut |value: &BatchParams| {
                selected = value.inner
            })?;
        }
        Ok(Self {
            inner: ck::BatchImageParams {
                format: if format.is_undefined() {
                    defaults.format
                } else {
                    format.as_string().ok_or_else(|| type_error("format"))?
                },
                width: if width.is_undefined() {
                    defaults.width
                } else {
                    u32_value(&width, "width")?
                },
                height: if height.is_undefined() {
                    defaults.height
                } else {
                    u32_value(&height, "height")?
                },
                execution: selected,
                filenames: if filenames.is_null() || filenames.is_undefined() {
                    None
                } else {
                    Some(
                        sequence(&filenames, "filenames")?
                            .iter()
                            .map(|name| {
                                if name.is_null() {
                                    Ok(None)
                                } else {
                                    name.as_string()
                                        .map(Some)
                                        .ok_or_else(|| type_error("filename"))
                                }
                            })
                            .collect::<Result<_, _>>()?,
                    )
                },
                report_path: if report_path.is_null() || report_path.is_undefined() {
                    None
                } else {
                    Some(
                        report_path
                            .as_string()
                            .ok_or_else(|| type_error("reportPath"))?
                            .into(),
                    )
                },
            },
        })
    }
    #[wasm_bindgen(getter)]
    pub fn format(&self) -> String {
        self.inner.format.clone()
    }
    #[wasm_bindgen(getter)]
    pub fn width(&self) -> u32 {
        self.inner.width
    }
    #[wasm_bindgen(getter)]
    pub fn height(&self) -> u32 {
        self.inner.height
    }
    #[wasm_bindgen(getter)]
    pub fn execution(&self) -> BatchParams {
        BatchParams {
            inner: self.inner.execution,
        }
    }
    #[wasm_bindgen(getter, unchecked_return_type = "(string | null)[] | null")]
    pub fn filenames(&self) -> JsValue {
        self.inner
            .filenames
            .as_ref()
            .map_or(JsValue::NULL, |values| {
                values
                    .iter()
                    .map(|name| {
                        name.as_ref()
                            .map_or(JsValue::NULL, |v| JsValue::from_str(v))
                    })
                    .collect::<Array>()
                    .into()
            })
    }
    #[wasm_bindgen(getter,js_name=reportPath,unchecked_return_type="string | null")]
    pub fn report_path(&self) -> JsValue {
        self.inner
            .report_path
            .as_ref()
            .map_or(JsValue::NULL, |path| {
                JsValue::from_str(&path.to_string_lossy())
            })
    }
}
// The callback only copies BatchParams. Borrowing failures become argument
// TypeError values while retaining wasm-bindgen's original error as cause.
#[cfg(feature = "cap-depict")]
#[wasm_bindgen(
    inline_js = "export function visit_execution(value, visit) { try { visit(value); } catch (cause) { throw new TypeError('invalid BatchParams', { cause }); } }"
)]
extern "C" {
    #[wasm_bindgen(catch)]
    fn visit_execution(value: &JsValue, visit: &mut dyn FnMut(&BatchParams))
    -> Result<(), JsValue>;
}
#[wasm_bindgen]
pub struct BatchExportReport {
    pub(crate) inner: ck::BatchExportReport,
}
#[wasm_bindgen]
impl BatchExportReport {
    pub fn total(&self) -> usize {
        self.inner.total()
    }
    pub fn success(&self) -> usize {
        self.inner.success()
    }
    pub fn failed(&self) -> usize {
        self.inner.failed()
    }
    pub fn errors(&self) -> Vec<BatchError> {
        self.inner
            .errors()
            .iter()
            .cloned()
            .map(|inner| BatchError { inner })
            .collect()
    }
    #[wasm_bindgen(getter)]
    pub fn written(&self) -> usize {
        self.inner.written
    }
    #[wasm_bindgen(js_name=writeReport)]
    pub fn write_report(&self, path: &str) -> Result<(), JsValue> {
        self.inner
            .write_report(std::path::Path::new(path))
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
}
#[cfg(feature = "cap-depict")]
#[wasm_bindgen]
impl MoleculeBatch {
    #[wasm_bindgen(js_name=writeImages)]
    pub fn write_images(&self, directory: &str) -> Result<BatchExportReport, JsValue> {
        self.inner
            .write_images(directory)
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
    #[wasm_bindgen(js_name=writeImagesWithParams)]
    pub fn write_images_with_params(
        &self,
        directory: &str,
        options: &BatchImageParams,
    ) -> Result<BatchExportReport, JsValue> {
        self.inner
            .write_images_with_params(directory, &options.inner)
            .map(|inner| BatchExportReport { inner })
            .map_err(|e| batch_validation_error(&e).unwrap_or_else(|e| e))
    }
}
