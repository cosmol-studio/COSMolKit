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

pyo3::create_exception!(
    cosmolkit,
    InchiError,
    PyValueError,
    "InChI conversion or key generation failed; diagnostic results expose engine messages and status."
);
/// Stable category for failures at the toolkit-neutral InChI boundary.
///
/// Declared values: ``AllocationFailed``, ``UnsupportedState``, ``InvalidInput``, ``InvalidSourceOutput``, ``SanitizeFailed``, ``Toolkit``, ``SourcePort``.
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

/// Writable configuration for InChI parsing and chemical preprocessing.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct InchiReadParams {
    pub(crate) inner: ck::InchiReadParams,
}
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl InchiReadParams {
    /// Configure InChI parsing and chemical preprocessing; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(sanitize=true, remove_hs=true))]
    fn new(sanitize: bool, remove_hs: bool) -> Self {
        Self {
            inner: ck::InchiReadParams::new(sanitize, remove_hs),
        }
    }
    /// Return a new result that will perform the selected chemical sanitization stages. The source molecule is unchanged.
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize
    }
    #[setter]
    fn set_sanitize(&mut self, value: bool) {
        self.inner.sanitize = value;
    }
    /// Whether removable explicit hydrogens are removed during input conversion.
    #[getter]
    fn remove_hs(&self) -> bool {
        self.inner.remove_hs
    }
    #[setter]
    fn set_remove_hs(&mut self, value: bool) {
        self.inner.remove_hs = value;
    }
}
/// Writable configuration for InChI/InChIKey generation options.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct InchiWriteParams {
    pub(crate) inner: ck::InchiWriteParams,
}
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl InchiWriteParams {
    /// Configure InChI/InChIKey generation options; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(options=String::new()))]
    fn new(options: String) -> Self {
        Self {
            inner: ck::InchiWriteParams::new(options),
        }
    }
    /// InChI option string passed to the InChI engine; an empty string selects its default options.
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
/// Convert InChI text directly into its InChIKey without constructing a molecule; invalid input raises InchiError.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn inchi_to_key(py: Python<'_>, inchi: crate::text_input::TextInput<'_>) -> PyResult<String> {
    ck::inchi_to_key(&inchi.as_text()?).map_err(|e| error(py, e))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<InchiReadParams>()?;
    module.add_class::<InchiWriteParams>()?;
    module.add_class::<InchiErrorKind>()?;
    module.add("InchiError", module.py().get_type::<InchiError>())?;
    module.add_function(wrap_pyfunction!(inchi_to_key, module)?)?;
    Ok(())
}
