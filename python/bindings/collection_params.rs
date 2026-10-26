use cosmolkit as ck;
use pyo3::prelude::*;

ck::project_configuration!(ck::LayeredFingerprintParams);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Explicit correspondence between one probe atom and one reference atom for alignment.
    pub struct AlignmentAtomMap {
        pub(crate) inner: ck::AlignmentAtomMap,
    }
}
#[pymethods]
impl AlignmentAtomMap {
    #[new]
    fn py_new(probe_atom: usize, reference_atom: usize) -> Self {
        Self::new(probe_atom, reference_atom)
    }
    #[staticmethod]
    /// Construct a correspondence using zero-based probe and reference atom indices.
    fn new(probe_atom: usize, reference_atom: usize) -> Self {
        Self {
            inner: ck::AlignmentAtomMap::new(probe_atom, reference_atom),
        }
    }
    #[getter]
    /// Zero-based atom index in the probe molecule.
    fn probe_atom(&self) -> usize {
        self.inner.probe_atom
    }
    #[getter]
    /// Zero-based atom index in the reference molecule.
    fn reference_atom(&self) -> usize {
        self.inner.reference_atom
    }
}

ck::project_configuration!(ck::AlignmentParameters);
