use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum BioSelectionParseErrorKind {
    Syntax,
    SeqidRange,
}

#[pyclass(module = "fixed_projection._native", extends=PyValueError)]
pub struct BioSelectionParseError {
    inner: ck::BioSelectionParseError,
}

#[pymethods]
impl BioSelectionParseError {
    #[getter]
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    fn variant(&self) -> BioSelectionParseErrorKind {
        match self.inner {
            ck::BioSelectionParseError::Syntax(_) => BioSelectionParseErrorKind::Syntax,
            ck::BioSelectionParseError::SeqidRange(_) => BioSelectionParseErrorKind::SeqidRange,
        }
    }
    #[getter]
    fn cid(&self) -> Option<&str> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => Some(value.cid()),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    #[getter]
    fn pos(&self) -> Option<usize> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => Some(value.pos()),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    #[getter]
    fn info(&self) -> Option<&str> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => value.info(),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

impl BioSelectionParseError {
    fn from_rust(inner: ck::BioSelectionParseError) -> PyErr {
        // This Rust enum deliberately declares no Error::source. Its public
        // syntax accessors are payloads, not invented exception causes.
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct BioSelection {
    pub(crate) inner: ck::BioSelection,
}

#[pymethods]
impl BioSelection {
    #[staticmethod]
    fn from_cid(cid: &str) -> PyResult<Self> {
        ck::BioSelection::from_cid(cid)
            .map(|inner| Self { inner })
            .map_err(BioSelectionParseError::from_rust)
    }
    fn to_cid(&self) -> String {
        self.inner.to_cid()
    }
}
