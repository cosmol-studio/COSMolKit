use crate::native::Fingerprint;
use cosmolkit as ck;
use pyo3::prelude::*;

impl crate::parameters::NativeScalar for ck::LayeredFingerprintLayers {
    type Repr = u32;
    fn from_scalar(value: u32) -> Self {
        Self::from_bits_retain(value)
    }
    fn scalar(self) -> u32 {
        self.bits()
    }
}
impl crate::parameters::Configuration for Fingerprint {
    type Native = ck::Fingerprint;
    fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
        Ok(Self { inner })
    }
    fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
        Ok(self.inner.clone())
    }
}
ck::project_configuration!(ck::LayeredFingerprintParams);

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
#[derive(Clone)]
pub struct AlignmentAtomMap {
    pub(crate) inner: ck::AlignmentAtomMap,
}
#[pymethods]
impl AlignmentAtomMap {
    #[staticmethod]
    fn new(probe_atom: usize, reference_atom: usize) -> Self {
        Self {
            inner: ck::AlignmentAtomMap::new(probe_atom, reference_atom),
        }
    }
    #[getter]
    fn probe_atom(&self) -> usize {
        self.inner.probe_atom
    }
    #[getter]
    fn reference_atom(&self) -> usize {
        self.inner.reference_atom
    }
}

impl crate::parameters::Configuration for AlignmentAtomMap {
    type Native = ck::AlignmentAtomMap;
    fn from_native(_py: Python<'_>, inner: Self::Native) -> PyResult<Self> {
        Ok(Self { inner })
    }
    fn to_native(&self, _py: Python<'_>) -> PyResult<Self::Native> {
        Ok(self.inner)
    }
}
ck::project_configuration!(ck::AlignmentParameters);
