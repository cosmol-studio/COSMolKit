//! InChI transport only; algorithms remain in the canonical facade.
use ::cosmolkit as ck;
use pyo3::{
    exceptions::{PyTypeError, PyValueError},
    prelude::*,
};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};

pyo3::create_exception!(cosmolkit, InchiError, PyValueError);
#[cosmolkit_macros::python_enum]
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum InchiErrorKind {
    AllocationFailed,
    UnsupportedState,
    InvalidInput,
    InvalidSourceOutput,
    SanitizeFailed,
    Toolkit,
    SourcePort,
}
pub(crate) fn error(py: Python<'_>, source: ck::InchiError) -> PyErr {
    let kind = match source.kind {
        ck::InchiErrorKind::AllocationFailed => InchiErrorKind::AllocationFailed,
        ck::InchiErrorKind::UnsupportedState => InchiErrorKind::UnsupportedState,
        ck::InchiErrorKind::InvalidInput => InchiErrorKind::InvalidInput,
        ck::InchiErrorKind::InvalidSourceOutput => InchiErrorKind::InvalidSourceOutput,
        ck::InchiErrorKind::SanitizeFailed => InchiErrorKind::SanitizeFailed,
        ck::InchiErrorKind::Toolkit => InchiErrorKind::Toolkit,
        ck::InchiErrorKind::SourcePort => InchiErrorKind::SourcePort,
    };
    let exception = InchiError::new_err(source.to_string());
    let fields = || -> PyResult<()> {
        let value = exception.value(py);
        value.setattr("domain", "inchi")?;
        value.setattr("kind", Py::new(py, kind)?)?;
        value.setattr("operation", source.operation)?;
        value.setattr("detail", source.detail)?;
        Ok(())
    };
    match fields() {
        Ok(()) => exception,
        Err(error) => error,
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct InchiReadParams {
    pub(crate) inner: ck::InchiReadParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl InchiReadParams {
    #[new]
    #[pyo3(signature=(sanitize=true, remove_hs=true))]
    fn new(sanitize: bool, remove_hs: bool) -> Self {
        Self {
            inner: ck::InchiReadParams::new(sanitize, remove_hs),
        }
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[setter]
    fn set_sanitize(&mut self, value: bool) {
        self.inner.sanitize = value;
    }
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[setter]
    fn set_remove_hs(&mut self, value: bool) {
        self.inner.remove_hs = value;
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct InchiWriteParams {
    pub(crate) inner: ck::InchiWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl InchiWriteParams {
    #[new]
    #[pyo3(signature=(options=String::new()))]
    fn new(options: String) -> Self {
        Self {
            inner: ck::InchiWriteParams::new(options),
        }
    }
    #[getter]
    fn options(&self) -> &str {
        &self.inner.options
    }
    #[setter]
    fn set_options(&mut self, value: String) {
        self.inner.options = value;
    }
}
pub(crate) fn write_params(
    params: Option<&InchiWriteParams>,
    options: Option<String>,
) -> PyResult<ck::InchiWriteParams> {
    if params.is_some() && options.is_some() {
        return Err(PyTypeError::new_err(
            "params and options are mutually exclusive",
        ));
    }
    Ok(match params {
        Some(p) => p.inner.clone(),
        None => ck::InchiWriteParams {
            options: options.unwrap_or_default(),
        },
    })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn inchi_to_key(py: Python<'_>, inchi: &str) -> PyResult<String> {
    ck::inchi_to_key(inchi).map_err(|e| error(py, e))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<InchiReadParams>()?;
    module.add_class::<InchiWriteParams>()?;
    module.add_class::<InchiErrorKind>()?;
    module.add("InchiError", module.py().get_type::<InchiError>())?;
    module.add_function(wrap_pyfunction!(inchi_to_key, module)?)?;
    Ok(())
}
