//! Forward-only supplier projection; indexed and batch suppliers share canonical_batch.
use crate::canonical_sdf::{SdfReadParams, SdfRecord};
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};

/// Incremental SDF record stream with read position and record/byte/line counters.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct SdfRecordStream {
    inner: Option<ck::SdfRecordStream<std::io::BufReader<std::fs::File>>>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl SdfRecordStream {
    /// Open the input SDF path for indexed or streaming access; does not interpret a text block as a path.
    #[staticmethod]
    fn open(py: Python<'_>, path: &str) -> PyResult<Self> {
        ck::SdfRecordStream::open(path)
            .map(|inner| Self { inner: Some(inner) })
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Open the input SDF path for indexed or streaming access; does not interpret a text block as a path. Uses the supplied configuration object.
    #[staticmethod]
    fn open_with_params(py: Python<'_>, path: &str, params: &SdfReadParams) -> PyResult<Self> {
        ck::SdfRecordStream::open_with_params(path, &params.inner)
            .map(|inner| Self { inner: Some(inner) })
            .map_err(|e| crate::canonical_molecular_io::error_pyerr(py, e))
    }
    /// Read the next SDF record, or None at end of input; malformed input remains an error.
    fn next_record(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        self.inner
            .as_mut()
            .ok_or_else(stream_transferred)?
            .next_record()
            .map(|r| r.map(|inner| SdfRecord { inner }))
            .map_err(|e| crate::canonical_sdf::sdf_pyerr(py, e))
    }
    /// Return whether the SDF stream has reached end of input.
    fn is_end(&self) -> PyResult<bool> {
        Ok(self.inner.as_ref().ok_or_else(stream_transferred)?.is_end())
    }
    /// Number of records consumed from this stream.
    fn records_consumed(&self) -> PyResult<usize> {
        Ok(self
            .inner
            .as_ref()
            .ok_or_else(stream_transferred)?
            .records_consumed())
    }
    /// Number of source bytes consumed from this stream.
    fn bytes_consumed(&self) -> PyResult<u64> {
        Ok(self
            .inner
            .as_ref()
            .ok_or_else(stream_transferred)?
            .bytes_consumed())
    }
    /// Number of source lines consumed from this stream.
    fn lines_consumed(&self) -> PyResult<usize> {
        Ok(self
            .inner
            .as_ref()
            .ok_or_else(stream_transferred)?
            .lines_consumed())
    }
    /// Return an iterator that reads records in bounded batches.
    #[pyo3(signature=(size,mode,n_jobs))]
    fn batches(
        &mut self,
        py: Python<'_>,
        size: usize,
        mode: &Bound<'_, PyAny>,
        n_jobs: Option<usize>,
    ) -> PyResult<crate::canonical_batch::SdfReaderBatchIterator> {
        let mode = crate::canonical_batch::error_mode(Some(mode))?;
        self.inner
            .take()
            .ok_or_else(stream_transferred)?
            .batches(size, mode, n_jobs)
            .map(|inner| crate::canonical_batch::SdfReaderBatchIterator { inner })
            .map_err(|e| crate::canonical_batch::batch_error(py, e))
    }
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    #[gen_stub(override_return_type(type_repr = "SdfRecord"))]
    fn __next__(&mut self, py: Python<'_>) -> PyResult<Option<SdfRecord>> {
        self.next_record(py)
    }
}

fn stream_transferred() -> PyErr {
    pyo3::exceptions::PyRuntimeError::new_err("SDF stream has been transferred to a batch iterator")
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<SdfRecordStream>()
}
