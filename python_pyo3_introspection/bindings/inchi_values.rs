use crate::chemical_values::enum_value;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

enum_value!(InchiErrorKind, {
    AllocationFailed, UnsupportedState, InvalidInput, InvalidSourceOutput,
    SanitizeFailed, Toolkit, SourcePort
});

#[pyclass(module="cosmolkit", extends=PyValueError)]
pub struct InchiError {
    inner: ck::InchiError,
}
#[pymethods]
impl InchiError {
    #[getter]
    fn variant(&self) -> InchiErrorKind {
        self.inner.kind.into()
    }
    #[getter]
    fn operation(&self) -> &'static str {
        self.inner.operation
    }
    #[getter]
    fn detail(&self) -> &str {
        &self.inner.detail
    }
    #[getter]
    fn domain(&self) -> &'static str {
        "inchi"
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl InchiError {
    pub(crate) fn from_rust(inner: ck::InchiError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[pyfunction]
pub fn inchi_to_key(inchi: &str) -> PyResult<String> {
    ck::inchi_to_key(inchi).map_err(InchiError::from_rust)
}

#[pyfunction]
#[pyo3(signature = (text, params=crate::arguments::Argument::Missing, *, sanitize=crate::arguments::Argument::Missing, remove_hs=crate::arguments::Argument::Missing))]
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
