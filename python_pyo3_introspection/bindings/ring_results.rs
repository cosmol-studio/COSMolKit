use pyo3::prelude::*;

crate::configuration::type_projection! {
    /// Owned symmetrized SSSR rows; indices refer to the queried molecular version.
    #[pyclass(module = "cosmolkit._native", frozen, from_py_object)]
    #[derive(Clone)]
    pub struct RingSearchResult {
        pub(crate) inner: cosmolkit::RingSearchResult,
    }
}

#[pymethods]
impl RingSearchResult {
    #[getter]
    fn atom_rings(&self) -> Vec<Vec<usize>> {
        self.inner.atom_rings().iter().map(|ring| ring.iter().map(|id| id.index()).collect()).collect()
    }
    #[getter]
    fn bond_rings(&self) -> Vec<Vec<usize>> {
        self.inner.bond_rings().iter().map(|ring| ring.iter().map(|id| id.index()).collect()).collect()
    }
    #[getter]
    fn num_rings(&self) -> usize {
        self.inner.num_rings()
    }
}
