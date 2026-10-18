use crate::chemical_values::enum_value;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

enum_value!(InchiErrorKind, {
    AllocationFailed, UnsupportedState, InvalidInput, InvalidSourceOutput,
    SanitizeFailed, Toolkit, SourcePort
});

#[pyclass(module="fixed_projection._native", extends=PyValueError)]
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
#[pyo3(signature = (text, params=None, *, sanitize=None, remove_hs=None))]
pub fn mol_from_inchi(
    text: &str,
    params: Option<&crate::params::InchiReadParams>,
    sanitize: Option<bool>,
    remove_hs: Option<bool>,
) -> PyResult<Molecule> {
    if let Some(params) = params {
        if sanitize.is_some() || remove_hs.is_some() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        _mol_from_inchi_params(text, Some(params))
    } else {
        _mol_from_inchi_keywords(text, sanitize, remove_hs)
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
