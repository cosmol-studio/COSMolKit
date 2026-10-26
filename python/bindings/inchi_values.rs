use crate::chemical_values::enum_value;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;

enum_value!(InchiDiagnosticLevel, { Warning, Error });
#[pyclass(module = "cosmolkit", frozen)]
/// InChI output with AuxInfo, native status, message, log and structured diagnostics.
pub struct InchiResult {
    pub(crate) inner: ck::InchiResult,
}
#[pymethods]
impl InchiResult {
    #[getter]
    /// Generated InChI identifier.
    fn inchi(&self) -> &str {
        self.inner.inchi()
    }
    #[getter]
    /// Native InChI return status code.
    fn return_code(&self) -> i32 {
        self.inner.return_code()
    }
    #[getter]
    /// Message emitted by the native InChI calculation.
    fn message(&self) -> &str {
        self.inner.message()
    }
    #[getter]
    /// Native InChI diagnostic log.
    fn log(&self) -> &str {
        self.inner.log()
    }
    #[getter]
    /// AuxInfo text accompanying the InChI result.
    fn aux_info(&self) -> &str {
        self.inner.aux_info()
    }
    #[getter]
    /// Independent list of structured warning/error diagnostics.
    fn diagnostics(&self) -> Vec<InchiDiagnostic> {
        self.inner
            .diagnostics()
            .iter()
            .cloned()
            .map(|inner| InchiDiagnostic { inner })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}
#[pyclass(module = "cosmolkit", frozen)]
/// One warning or error message produced by the InChI calculation.
pub struct InchiDiagnostic {
    pub(crate) inner: ck::InchiDiagnostic,
}
#[pymethods]
impl InchiDiagnostic {
    #[getter]
    /// Severity of the InChI diagnostic.
    fn level(&self) -> InchiDiagnosticLevel {
        self.inner.level().into()
    }
    #[getter]
    /// Message emitted by the native InChI calculation.
    fn message(&self) -> &str {
        self.inner.message()
    }
    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }
}

enum_value!(InchiErrorKind, {
    AllocationFailed, UnsupportedState, InvalidInput, InvalidSourceOutput,
    SanitizeFailed, Toolkit, SourcePort
});

crate::error_framework::error_projection! {record
pub struct InchiError {
    inner: ck::InchiError,
}
}
crate::error_framework::error_projection! {methods impl InchiError {
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
    fn variant(&self) -> InchiErrorKind {
        self.inner.kind.into()
    }
    #[getter]
    /// Native InChI operation reporting the failure.
    fn operation(&self) -> &'static str {
        self.inner.operation
    }
    #[getter]
    /// Native explanation of the InChI failure.
    fn detail(&self) -> &str {
        &self.inner.detail
    }
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "inchi"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl InchiError {
    pub(crate) fn from_rust(inner: ck::InchiError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

#[pyfunction]
/// Convert an InChI string to its InChIKey; invalid input retains the native structured error.
pub fn inchi_to_key(inchi: &str) -> PyResult<String> {
    ck::inchi_to_key(inchi).map_err(InchiError::from_rust)
}

#[pyfunction]
#[pyo3(signature = (text, params=crate::arguments::Argument::Missing, *, sanitize=crate::arguments::Argument::Missing, remove_hs=crate::arguments::Argument::Missing))]
/// Parse an InChI string into a molecular value with the selected sanitization and hydrogen-removal policy.
pub fn mol_from_inchi(
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, crate::params::InchiReadParams>>>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    remove_hs: crate::arguments::Argument<Option<bool>>,
) -> PyResult<Molecule> {
    let keyword_form = sanitize.is_supplied() || remove_hs.is_supplied();
    crate::arguments::configuration_arguments!(params; sanitize, remove_hs);
    let params = params.as_deref();
    if keyword_form {
        _mol_from_inchi_keywords(text, sanitize, remove_hs)
    } else {
        _mol_from_inchi_params(text, params)
    }
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_inchi"}"#]
#[pyo3(signature = (text, params=None))]
pub fn _mol_from_inchi_params(
    text: &str,
    params: Option<&crate::params::InchiReadParams>,
) -> PyResult<Molecule> {
    Molecule::_from_inchi_params(text, params)
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"mol_from_inchi"}"#]
#[pyo3(signature = (text, *, sanitize=None, remove_hs=None))]
pub fn _mol_from_inchi_keywords(
    text: &str,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
) -> PyResult<Molecule> {
    Molecule::_from_inchi_keywords(text, sanitize, remove_hs)
}
