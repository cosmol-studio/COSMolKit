//! Frozen transport values for the canonical typed batch API.
use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::sync::{Arc, Mutex};

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BatchParams {
    pub(crate) inner: ck::BatchParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BatchExportParams {
    pub(crate) inner: ck::BatchExportParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchExportParams {
    #[new]
    #[pyo3(signature=(*,format=None,errors=None,n_jobs=None,progress_bar=None))]
    fn new(
        format: Option<&str>,
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BatchExportParams {
                format: crate::canonical_batch::write_params(format, true, true)?.format,
                errors: crate::canonical_batch::error_mode(errors)?,
                n_jobs: crate::canonical_batch::n_jobs(n_jobs)?,
                progress_bar,
            },
        })
    }
    #[getter]
    fn format(&self) -> &'static str {
        match self.inner.format {
            ck::SdfFormat::V2000 => "v2000",
            ck::SdfFormat::V3000 => "v3000",
        }
    }
    #[getter]
    fn errors(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        BatchParams {
            inner: ck::BatchParams {
                errors: self.inner.errors,
                n_jobs: self.inner.n_jobs,
                progress_bar: self.inner.progress_bar,
            },
        }
        .errors(py)
    }
    #[getter]
    fn n_jobs(&self) -> Option<usize> {
        self.inner.n_jobs
    }
    #[getter]
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchParams {
    #[new]
    #[pyo3(signature=(*,errors=None,n_jobs=None,progress_bar=None))]
    fn new(
        errors: Option<&Bound<'_, PyAny>>,
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: ck::BatchParams {
                errors: crate::canonical_batch::error_mode(errors)?,
                n_jobs,
                progress_bar,
            },
        })
    }
    #[getter]
    fn errors(&self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        let name = match self.inner.errors {
            ck::BatchErrorMode::Strict => "RAISE",
            ck::BatchErrorMode::KeepErrors => "KEEP",
        };
        Ok(py
            .import("cosmolkit")?
            .getattr("BatchErrorMode")?
            .getattr(name)?
            .unbind())
    }
    #[getter]
    fn n_jobs(&self) -> Option<usize> {
        self.inner.n_jobs
    }
    #[getter]
    fn progress_bar(&self) -> Option<bool> {
        self.inner.progress_bar
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BatchQueryParams {
    n_jobs: Option<usize>,
    progress_bar: Option<bool>,
    progress_callback: Option<Py<PyAny>>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchQueryParams {
    #[new]
    #[pyo3(signature=(*,n_jobs=None,progress_bar=None,progress_callback=None))]
    fn new(
        n_jobs: Option<usize>,
        progress_bar: Option<bool>,
        progress_callback: Option<Py<PyAny>>,
        py: Python<'_>,
    ) -> PyResult<Self> {
        if progress_callback
            .as_ref()
            .is_some_and(|callback| !callback.bind(py).is_callable())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "progress_callback must be callable",
            ));
        }
        Ok(Self {
            n_jobs,
            progress_bar,
            progress_callback,
        })
    }
    #[getter]
    fn n_jobs(&self) -> Option<usize> {
        self.n_jobs
    }
    #[getter]
    fn progress_bar(&self) -> Option<bool> {
        self.progress_bar
    }
    #[getter]
    fn progress_callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.progress_callback
            .as_ref()
            .map(|callback| callback.clone_ref(py))
    }
}
impl BatchQueryParams {
    pub(crate) fn execute<T: Send>(
        &self,
        py: Python<'_>,
        call: impl FnOnce(&ck::BatchQueryParams) -> Result<T, ck::BatchValidationError> + Send,
    ) -> PyResult<T> {
        // Callback transport only: the canonical facade owns row order/tick count.
        // Python exceptions remain observable and are returned after canonical work;
        // the first callback failure is retained, with no silent fallback or panic.
        let callback_error = Arc::new(Mutex::new(None));
        let progress_callback = self.progress_callback.as_ref().map(|callback| {
            let callback = callback.clone_ref(py);
            let error = Arc::clone(&callback_error);
            Arc::new(move || {
                Python::attach(|py| {
                    if let Err(source) = callback.call0(py) {
                        let mut first = error.lock().expect("callback error mutex");
                        if first.is_none() {
                            *first = Some(source);
                        }
                    }
                })
            }) as Arc<dyn Fn() + Send + Sync>
        });
        let params = ck::BatchQueryParams {
            n_jobs: self.n_jobs,
            progress_bar: self.progress_bar,
            progress_callback,
        };
        // Worker callbacks acquire Python's GIL. Release the caller's GIL for
        // canonical Rust work so parallel callbacks cannot deadlock its owner.
        let result = py
            .detach(|| call(&params))
            .map_err(|source| crate::canonical_batch::batch_error(py, source))?;
        let error = callback_error.lock().expect("callback error mutex").take();
        match error {
            Some(error) => Err(error),
            None => Ok(result),
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct BatchImageParams {
    pub(crate) inner: ck::BatchImageParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl BatchImageParams {
    #[new]
    #[pyo3(signature=(*,format="png",width=300,height=300,execution=None,filenames=None,report_path=None))]
    fn new(
        format: &str,
        width: u32,
        height: u32,
        execution: Option<&BatchParams>,
        filenames: Option<Vec<Option<String>>>,
        report_path: Option<String>,
    ) -> Self {
        Self {
            inner: ck::BatchImageParams {
                format: format.into(),
                width,
                height,
                execution: execution.map(|value| value.inner).unwrap_or_default(),
                filenames,
                report_path: report_path.map(Into::into),
            },
        }
    }
    #[getter]
    fn format(&self) -> String {
        self.inner.format.clone()
    }
    #[getter]
    fn width(&self) -> u32 {
        self.inner.width
    }
    #[getter]
    fn height(&self) -> u32 {
        self.inner.height
    }
    #[getter]
    fn execution(&self) -> BatchParams {
        BatchParams {
            inner: self.inner.execution,
        }
    }
    #[getter]
    fn filenames(&self) -> Option<Vec<Option<String>>> {
        self.inner.filenames.clone()
    }
    #[getter]
    fn report_path(&self) -> Option<String> {
        self.inner
            .report_path
            .as_ref()
            .map(|value| value.to_string_lossy().into_owned())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<BatchExportParams>()?;
    module.add_class::<BatchParams>()?;
    module.add_class::<BatchQueryParams>()?;
    module.add_class::<BatchImageParams>()?;
    Ok(())
}
