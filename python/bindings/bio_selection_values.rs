use cosmolkit as ck;
use pyo3::prelude::*;

crate::enums::binding_enum!(BioSelectionParseErrorKind, {
    Syntax,
    SeqidRange,
});

crate::error_framework::error_projection! {record
pub struct BioSelectionParseError {
    inner: ck::BioSelectionParseError,
}
}

crate::error_framework::error_projection! {methods impl BioSelectionParseError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "bio"
    }
    #[getter]
    /// Exact native error category. Variant-specific payload fields are None when that variant is not active.
    fn variant(&self) -> BioSelectionParseErrorKind {
        match self.inner {
            ck::BioSelectionParseError::Syntax(_) => BioSelectionParseErrorKind::Syntax,
            ck::BioSelectionParseError::SeqidRange(_) => BioSelectionParseErrorKind::SeqidRange,
        }
    }
    #[getter]
    #[doc = "CID expression that failed syntax parsing; None for sequence-ID range failures."] fn cid(&self) -> Option<&str> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => Some(value.cid()),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    #[getter]
    #[doc = "Syntax failure position in the CID expression; None for sequence-ID range failures."] fn pos(&self) -> Option<usize> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => Some(value.pos()),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    #[getter]
    #[doc = "Optional syntax-parser explanation; None for sequence-ID range failures."] fn info(&self) -> Option<&str> {
        match &self.inner {
            ck::BioSelectionParseError::Syntax(value) => value.info(),
            ck::BioSelectionParseError::SeqidRange(_) => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}

impl BioSelectionParseError {
    fn from_rust(inner: ck::BioSelectionParseError) -> PyErr {
        // This Rust enum deliberately declares no Error::source. Its public
        // syntax accessors are payloads, not invented exception causes.
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}

#[pyclass(module = "cosmolkit", frozen)]
/// Parsed CID selector for models, chains, residues and atoms in a biological structure.
pub struct BioSelection {
    pub(crate) inner: ck::BioSelection,
}

#[pymethods]
impl BioSelection {
    #[staticmethod]
    /// Parse a CID expression into a reusable structure selector.
    fn from_cid(cid: &str) -> PyResult<Self> {
        ck::BioSelection::from_cid(cid)
            .map(|inner| Self { inner })
            .map_err(BioSelectionParseError::from_rust)
    }
    /// Serialize this selection as a CID expression.
    fn to_cid(&self) -> String {
        self.inner.to_cid()
    }
}
