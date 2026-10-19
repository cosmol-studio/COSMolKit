//! Native indexed and streaming suppliers. Cursor movement remains in Rust.
use crate::{
    batch_errors::BatchValidationError,
    extension::MoleculeBatch,
    io_batch_params::BatchErrorModeInput,
    molecular_io_errors::MolecularIoError,
    sdf_records::{SdfError, SdfRecord},
    selectors::SdfReadParams,
};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::{fs::File, io::BufReader, path::PathBuf};

#[pyclass(module = "cosmolkit", frozen)]
pub struct SdfDataset {
    pub(crate) inner: ck::SdfDataset,
}
#[pymethods]
impl SdfDataset {
    #[pyo3(signature=(index,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    fn record(
        &self,
        py: Python<'_>,
        index: usize,
        params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<
            Option<crate::selectors::SdfCoordinateModeInput>,
        >,
    ) -> PyResult<SdfRecord> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || strict_parsing.is_supplied()
            || expand_attachment_points.is_supplied()
            || process_property_lists.is_supplied()
            || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        let params = params.as_deref();
        if keyword_form {
            self._record_keywords(
                py,
                index,
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )
        } else {
            self._record_params(py, index, params)
        }
    }
    #[doc = r#"{"role":"overload","target":"SdfDataset.record"}"#]
    #[pyo3(signature=(index,params=None))]
    fn _record_params(
        &self,
        py: Python<'_>,
        index: usize,
        params: Option<&SdfReadParams>,
    ) -> PyResult<SdfRecord> {
        py.detach(|| match params {
            Some(p) => self.inner.record_with_params(index, &p.inner),
            None => self.inner.record(index),
        })
        .map(SdfRecord::from_rust)
        .map_err(SdfError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"SdfDataset.record"}"#]
    #[pyo3(signature=(index,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    fn _record_keywords(
        &self,
        py: Python<'_>,
        index: usize,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<SdfRecord> {
        let p = crate::molecular_io_calls::sdf_keywords(
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )?;
        self._record_params(py, index, Some(&p))
    }
    fn record_with_params(
        &self,
        py: Python<'_>,
        index: usize,
        params: &SdfReadParams,
    ) -> PyResult<SdfRecord> {
        self._record_params(py, index, Some(params))
    }

    #[staticmethod]
    #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    fn open(
        py: Python<'_>,
        path: PathBuf,
        params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<
            Option<crate::selectors::SdfCoordinateModeInput>,
        >,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || strict_parsing.is_supplied()
            || expand_attachment_points.is_supplied()
            || process_property_lists.is_supplied()
            || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        let params = params.as_deref();
        if keyword_form {
            Self::_open_keywords(
                py,
                path,
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )
        } else {
            Self::_open_params(py, path, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfDataset.open"}"#]
    #[pyo3(signature=(path,params=None))]
    fn _open_params(
        py: Python<'_>,
        path: PathBuf,
        params: Option<&SdfReadParams>,
    ) -> PyResult<Self> {
        let path = crate::molecular_io_calls::path_text(&path)?;
        py.detach(|| match params {
            Some(p) => ck::SdfDataset::open_with_params(path, &p.inner),
            None => ck::SdfDataset::open(path),
        })
        .map(|inner| Self { inner })
        .map_err(MolecularIoError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfDataset.open"}"#]
    #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    fn _open_keywords(
        py: Python<'_>,
        path: PathBuf,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let p = crate::molecular_io_calls::sdf_keywords(
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )?;
        Self::_open_params(py, path, Some(&p))
    }
    #[staticmethod]
    fn open_with_params(py: Python<'_>, path: PathBuf, params: &SdfReadParams) -> PyResult<Self> {
        Self::_open_params(py, path, Some(params))
    }

    fn len(&self) -> usize {
        self.inner.len()
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn path(&self) -> PathBuf {
        self.inner.path().to_path_buf()
    }
    fn metadata(&self, index: usize) -> Option<crate::detached_values::SdfRecordMetadata> {
        self.inner
            .metadata(index)
            .cloned()
            .map(|inner| crate::detached_values::SdfRecordMetadata { inner })
    }
    fn record_text(&self, py: Python<'_>, index: usize) -> PyResult<String> {
        py.detach(|| self.inner.record_text(index))
            .map_err(SdfError::from_rust)
    }
    fn iter(&self) -> SdfDatasetIterator {
        SdfDatasetIterator {
            inner: self.inner.iter(),
        }
    }
    fn __iter__(&self) -> SdfDatasetIterator {
        self.iter()
    }
    #[pyo3(signature=(size,mode=None,n_jobs=None,progress_bar=false))]
    fn batches(
        &self,
        py: Python<'_>,
        size: usize,
        mode: Option<BatchErrorModeInput>,
        n_jobs: Option<usize>,
        progress_bar: bool,
    ) -> PyResult<SdfBatchIterator> {
        let mode = mode
            .map(BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        py.detach(|| self.inner.batches(size, mode, n_jobs, progress_bar))
            .map(|inner| SdfBatchIterator { inner })
            .map_err(BatchValidationError::from_rust)
    }
}
#[pyclass(module = "cosmolkit")]
pub struct SdfDatasetIterator {
    inner: ck::SdfDatasetIterator,
}
#[pymethods]
impl SdfDatasetIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        py.detach(|| self.inner.next())
            .transpose()
            .map(|v| v.map(SdfRecord::from_rust))
            .map_err(SdfError::from_rust)
    }
}
#[pyclass(module = "cosmolkit")]
pub struct SdfRecordStream {
    inner: Option<ck::SdfRecordStream<BufReader<File>>>,
}
impl SdfRecordStream {
    fn active(&self) -> PyResult<&ck::SdfRecordStream<BufReader<File>>> {
        self.inner.as_ref().ok_or_else(Self::transferred)
    }
    fn transferred() -> PyErr {
        pyo3::exceptions::PyRuntimeError::new_err(
            "SdfRecordStream ownership has been transferred to batches; the stream is no longer usable",
        )
    }
}
#[pymethods]
impl SdfRecordStream {
    #[staticmethod]
    #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    fn open(
        py: Python<'_>,
        path: PathBuf,
        params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<
            Option<crate::selectors::SdfCoordinateModeInput>,
        >,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || strict_parsing.is_supplied()
            || expand_attachment_points.is_supplied()
            || process_property_lists.is_supplied()
            || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        let params = params.as_deref();
        if keyword_form {
            Self::_open_keywords(
                py,
                path,
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )
        } else {
            Self::_open_params(py, path, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfRecordStream.open"}"#]
    #[pyo3(signature=(path,params=None))]
    fn _open_params(
        py: Python<'_>,
        path: PathBuf,
        params: Option<&SdfReadParams>,
    ) -> PyResult<Self> {
        let path = crate::molecular_io_calls::path_text(&path)?;
        py.detach(|| match params {
            Some(p) => ck::SdfRecordStream::open_with_params(path, &p.inner),
            None => ck::SdfRecordStream::open(path),
        })
        .map(|inner| Self { inner: Some(inner) })
        .map_err(MolecularIoError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfRecordStream.open"}"#]
    #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    fn _open_keywords(
        py: Python<'_>,
        path: PathBuf,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let p = crate::molecular_io_calls::sdf_keywords(
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )?;
        Self::_open_params(py, path, Some(&p))
    }
    #[staticmethod]
    fn open_with_params(py: Python<'_>, path: PathBuf, params: &SdfReadParams) -> PyResult<Self> {
        Self::_open_params(py, path, Some(params))
    }

    fn next_record(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        let inner = self.inner.as_mut().ok_or_else(Self::transferred)?;
        py.detach(|| inner.next_record())
            .map(|v| v.map(SdfRecord::from_rust))
            .map_err(SdfError::from_rust)
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<PyRef<'_, Self>> {
        slf.active()?;
        Ok(slf)
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        self.next_record(py)
    }
    fn is_end(&self) -> PyResult<bool> {
        Ok(self.active()?.is_end())
    }
    fn records_consumed(&self) -> PyResult<usize> {
        Ok(self.active()?.records_consumed())
    }
    fn bytes_consumed(&self) -> PyResult<u64> {
        Ok(self.active()?.bytes_consumed())
    }
    fn lines_consumed(&self) -> PyResult<usize> {
        Ok(self.active()?.lines_consumed())
    }
    /// Transfer the current stream to a batch iterator. This stream becomes
    /// unusable even if native batch-parameter validation returns an error.
    #[pyo3(signature=(size,mode=None,n_jobs=None))]
    fn batches(
        &mut self,
        py: Python<'_>,
        size: usize,
        mode: Option<BatchErrorModeInput>,
        n_jobs: Option<usize>,
    ) -> PyResult<SdfReaderBatchIterator> {
        self.active()?;
        let mode = mode
            .map(BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        let inner = self.inner.take().ok_or_else(Self::transferred)?;
        py.detach(|| inner.batches(size, mode, n_jobs))
            .map(|inner| SdfReaderBatchIterator { inner })
            .map_err(BatchValidationError::from_rust)
    }
}
#[pyclass(module = "cosmolkit", frozen)]
pub struct SdfReader {
    inner: ck::SdfReader,
}
#[pymethods]
impl SdfReader {
    #[staticmethod]
    #[pyo3(signature=(path,params=crate::arguments::Argument::Missing,*,sanitize=crate::arguments::Argument::Missing,remove_hs=crate::arguments::Argument::Missing,strict_parsing=crate::arguments::Argument::Missing,expand_attachment_points=crate::arguments::Argument::Missing,process_property_lists=crate::arguments::Argument::Missing,coordinate_mode=crate::arguments::Argument::Missing))]
    fn open(
        py: Python<'_>,
        path: PathBuf,
        params: crate::arguments::Argument<Option<PyRef<'_, SdfReadParams>>>,
        sanitize: crate::arguments::Argument<Option<bool>>,
        remove_hs: crate::arguments::Argument<Option<bool>>,
        strict_parsing: crate::arguments::Argument<Option<bool>>,
        expand_attachment_points: crate::arguments::Argument<Option<bool>>,
        process_property_lists: crate::arguments::Argument<Option<bool>>,
        coordinate_mode: crate::arguments::Argument<
            Option<crate::selectors::SdfCoordinateModeInput>,
        >,
    ) -> PyResult<Self> {
        let keyword_form = sanitize.is_supplied()
            || remove_hs.is_supplied()
            || strict_parsing.is_supplied()
            || expand_attachment_points.is_supplied()
            || process_property_lists.is_supplied()
            || coordinate_mode.is_supplied();
        crate::arguments::configuration_arguments!(params; sanitize, remove_hs, strict_parsing, expand_attachment_points, process_property_lists, coordinate_mode);
        let params = params.as_deref();
        if keyword_form {
            Self::_open_keywords(
                py,
                path,
                sanitize,
                remove_hs,
                strict_parsing,
                expand_attachment_points,
                process_property_lists,
                coordinate_mode,
            )
        } else {
            Self::_open_params(py, path, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfReader.open"}"#]
    #[pyo3(signature=(path,params=None))]
    fn _open_params(
        py: Python<'_>,
        path: PathBuf,
        params: Option<&SdfReadParams>,
    ) -> PyResult<Self> {
        let path = crate::molecular_io_calls::path_text(&path)?;
        py.detach(|| match params {
            Some(p) => ck::SdfReader::open_with_params(path, &p.inner),
            None => ck::SdfReader::open(path),
        })
        .map(|inner| Self { inner })
        .map_err(MolecularIoError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"SdfReader.open"}"#]
    #[pyo3(signature=(path,*,sanitize=None,remove_hs=None,strict_parsing=None,expand_attachment_points=None,process_property_lists=None,coordinate_mode=None))]
    fn _open_keywords(
        py: Python<'_>,
        path: PathBuf,
        sanitize: Option<bool>,
        remove_hs: Option<bool>,
        strict_parsing: Option<bool>,
        expand_attachment_points: Option<bool>,
        process_property_lists: Option<bool>,
        coordinate_mode: Option<crate::selectors::SdfCoordinateModeInput>,
    ) -> PyResult<Self> {
        let p = crate::molecular_io_calls::sdf_keywords(
            sanitize,
            remove_hs,
            strict_parsing,
            expand_attachment_points,
            process_property_lists,
            coordinate_mode,
        )?;
        Self::_open_params(py, path, Some(&p))
    }
    #[staticmethod]
    fn open_with_params(py: Python<'_>, path: PathBuf, params: &SdfReadParams) -> PyResult<Self> {
        Self::_open_params(py, path, Some(params))
    }

    fn path(&self) -> PathBuf {
        self.inner.path().to_path_buf()
    }
    fn params(&self) -> SdfReadParams {
        SdfReadParams {
            inner: *self.inner.params(),
        }
    }
    #[pyo3(signature=(size,mode=None,n_jobs=None))]
    fn batches(
        &self,
        py: Python<'_>,
        size: usize,
        mode: Option<BatchErrorModeInput>,
        n_jobs: Option<usize>,
    ) -> PyResult<SdfReaderBatchIterator> {
        let mode = mode
            .map(BatchErrorModeInput::into_rust)
            .transpose()?
            .unwrap_or_default();
        py.detach(|| self.inner.batches(size, mode, n_jobs))
            .map(|inner| SdfReaderBatchIterator { inner })
            .map_err(BatchValidationError::from_rust)
    }
}
#[pyclass(module = "cosmolkit")]
pub struct SdfBatchIterator {
    inner: ck::SdfBatchIterator,
}
#[pymethods]
impl SdfBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn next_batch(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        py.detach(|| self.inner.next_batch())
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(BatchValidationError::from_rust)
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.next_batch(py)
    }
}
#[pyclass(module = "cosmolkit")]
pub struct SdfReaderBatchIterator {
    inner: ck::SdfReaderBatchIterator<BufReader<File>>,
}
#[pymethods]
impl SdfReaderBatchIterator {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn next_batch(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        py.detach(|| self.inner.next_batch())
            .map(|v| v.map(|inner| MoleculeBatch { inner }))
            .map_err(BatchValidationError::from_rust)
    }
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<MoleculeBatch>> {
        self.next_batch(py)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stream_batches_transfer_cursor_ownership_and_native_validation_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
                .join("artifacts")
                .join(format!("stream-transfer-{}", std::process::id()));
            std::fs::create_dir_all(&directory).unwrap();
            let path = directory.join("source.sdf");
            let text: String = ["C", "CO", "CCO"]
                .into_iter()
                .map(|smiles| {
                    ck::Molecule::from_smiles(smiles)
                        .unwrap()
                        .to_sdf_with_params(&ck::MolBlockWriteParams {
                            include_coordinates: false,
                            ..Default::default()
                        })
                        .unwrap()
                })
                .collect();
            std::fs::write(&path, text).unwrap();
            let open = || {
                py.get_type::<SdfRecordStream>()
                    .call_method1("open", (&path,))
            };
            let assert_transferred = |stream: &Bound<'_, pyo3::types::PyAny>| -> PyResult<()> {
                for method in [
                    "next_record",
                    "__next__",
                    "__iter__",
                    "is_end",
                    "records_consumed",
                    "bytes_consumed",
                    "lines_consumed",
                ] {
                    let error = stream.call_method0(method).unwrap_err();
                    assert!(error.is_instance_of::<pyo3::exceptions::PyRuntimeError>(py));
                    assert!(error.to_string().contains("ownership has been transferred"));
                }
                let error = stream.call_method1("batches", (1,)).unwrap_err();
                assert!(error.is_instance_of::<pyo3::exceptions::PyRuntimeError>(py));
                Ok(())
            };

            let stream = open()?;
            let alias = stream.clone();
            assert!(!stream.call_method0("next_record")?.is_none());
            assert_eq!(
                stream
                    .call_method0("records_consumed")?
                    .extract::<usize>()?,
                1
            );
            let mut native = ck::SdfRecordStream::open(path.to_str().unwrap()).unwrap();
            native.next_record().unwrap();
            let mut expected = native
                .batches(1, ck::BatchErrorMode::default(), Some(1))
                .unwrap();
            let iterator = stream.call_method1("batches", (1, "strict", 1))?;
            assert_transferred(&stream)?;
            assert_transferred(&alias)?;
            drop(stream);
            drop(alias);
            let mut count = 0;
            for batch in iterator.try_iter()? {
                let batch = batch?.extract::<Py<MoleculeBatch>>()?;
                let actual = batch.borrow(py);
                let expected = expected.next_batch().unwrap().unwrap();
                assert_eq!(actual.inner.records().len(), expected.records().len());
                for (actual, expected) in actual.inner.records().iter().zip(expected.records()) {
                    match (actual, expected) {
                        (
                            ck::BatchRecord::Molecule(actual),
                            ck::BatchRecord::Molecule(expected),
                        ) => {
                            assert_eq!(actual, expected);
                        }
                        _ => panic!("valid SDF fixture must contain molecule records"),
                    }
                }
                count += 1;
            }
            assert_eq!(count, 2);
            assert!(iterator.call_method0("next_batch")?.is_none());
            assert!(expected.next_batch().unwrap().is_none());

            for (size, jobs) in [(0, None), (1, Some(0))] {
                let stream = open()?;
                let native = ck::SdfRecordStream::open(path.to_str().unwrap()).unwrap();
                let expected = native
                    .batches(size, ck::BatchErrorMode::default(), jobs)
                    .err()
                    .unwrap();
                let error = stream
                    .call_method1("batches", (size, py.None(), jobs))
                    .unwrap_err();
                assert!(error.is_instance_of::<BatchValidationError>(py));
                assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
                assert_transferred(&stream)?;
            }

            // Python argument conversion precedes the consuming Rust call.
            let stream = open()?;
            assert!(stream.call_method1("batches", ("invalid",)).is_err());
            assert_eq!(
                stream
                    .call_method0("records_consumed")?
                    .extract::<usize>()?,
                0
            );
            assert!(!stream.call_method0("next_record")?.is_none());
            Ok(())
        })
        .unwrap();
    }
}
