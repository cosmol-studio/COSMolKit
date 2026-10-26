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
    /// Ordered atom-index cycles for the symmetrized SSSR, returned as independent lists.
    fn atom_rings(&self) -> Vec<Vec<usize>> {
        self.inner
            .atom_rings()
            .iter()
            .map(|ring| ring.iter().map(|id| id.index()).collect())
            .collect()
    }
    #[getter]
    /// Ordered bond-index cycles corresponding to atom_rings, returned as independent lists.
    fn bond_rings(&self) -> Vec<Vec<usize>> {
        self.inner
            .bond_rings()
            .iter()
            .map(|ring| ring.iter().map(|id| id.index()).collect())
            .collect()
    }
    #[getter]
    /// Number of rings in this symmetrized SSSR result.
    fn num_rings(&self) -> usize {
        self.inner.num_rings()
    }
}
