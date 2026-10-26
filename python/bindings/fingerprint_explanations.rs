use super::errors::AtomCodeExplanationError;
use cosmolkit as ck;
use pyo3::prelude::*;

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen)]
    /// Decoded element, branch, pi-electron and optional chirality components of an atom-pair atom code.
    pub struct AtomCodeExplanation {
        inner: ck::AtomCodeExplanation,
    }
}

#[pymethods]
impl AtomCodeExplanation {
    #[staticmethod]
    #[pyo3(signature = (code, branch_subtract=0, include_chirality=false))]
    /// Decode an atom-pair atom code, optionally adjusting branch count and interpreting chirality bits.
    fn from_code(code: u64, branch_subtract: i64, include_chirality: bool) -> PyResult<Self> {
        ck::AtomCodeExplanation::from_code(code, branch_subtract, include_chirality)
            .map(|inner| Self { inner })
            .map_err(AtomCodeExplanationError::from_rust)
    }
    #[getter]
    /// Element symbol decoded from the atom code.
    fn symbol(&self) -> &'static str {
        self.inner.symbol()
    }
    #[getter]
    /// Branch count decoded from the atom code after the requested subtraction.
    fn branch_count(&self) -> u32 {
        self.inner.branch_count()
    }
    #[getter]
    /// Pi-electron count encoded in the atom code.
    fn pi_electrons(&self) -> u32 {
        self.inner.pi_electrons()
    }
    #[getter]
    /// Decoded chirality label, or None when chirality is not included.
    fn chirality(&self) -> Option<&'static str> {
        self.inner.chirality()
    }
}

#[pyfunction]
#[pyo3(signature = (score, size=4))]
/// Decode a topological-torsion path score into the canonical atom-code components for the specified path size.
pub fn explain_path_score(score: u64, size: usize) -> Vec<(&'static str, u32, u32)> {
    ck::explain_path_score(score, size)
}
