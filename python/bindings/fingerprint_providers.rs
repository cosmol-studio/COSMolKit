use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

ck::project_configuration!(ck::AtomPairAtomInvariantsGenerator, methods {
    /// Return the native invariant-generator settings summary.
    fn info_string(&self, py: Python<'_>) -> PyResult<String> { self.inner.read(py, |native| native.info_string()) }
    /// Serialize the invariant-generator settings using the native JSON format.
    fn to_json(&self, py: Python<'_>) -> PyResult<String> { self.inner.read(py, |native| native.to_json()) }
});

ck::project_configuration!(ck::MorganBondInvariantsGenerator);

#[pyclass(module = "cosmolkit")]
/// Mutable native fingerprint provenance collector. Allocate the desired outputs before passing it to fingerprint generation.
pub struct FingerprintAdditionalOutput {
    pub(crate) inner: ck::FingerprintAdditionalOutput,
}

#[pymethods]
impl FingerprintAdditionalOutput {
    #[new]
    fn construct() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::new(),
        }
    }
    #[staticmethod]
    /// Create an empty additional-output collector; enable individual outputs before fingerprint generation.
    fn new() -> Self {
        Self::construct()
    }
    #[staticmethod]
    /// Create a collector using the native default allocation policy.
    fn default() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::default(),
        }
    }
    /// Enable collection of per-atom contribution counts.
    fn allocate_atom_counts(&mut self) {
        self.inner.allocate_atom_counts();
    }
    /// Enable collection of contributing bit identifiers per atom.
    fn allocate_atom_to_bits(&mut self) {
        self.inner.allocate_atom_to_bits();
    }
    /// Enable collection of generator-specific bit provenance tuples.
    fn allocate_bit_info_map(&mut self) {
        self.inner.allocate_bit_info_map();
    }
    /// Enable collection of contributing paths per bit.
    fn allocate_bit_paths(&mut self) {
        self.inner.allocate_bit_paths();
    }
    /// Enable collection of participating atom sets per bit.
    fn allocate_atoms_per_bit(&mut self) {
        self.inner.allocate_atoms_per_bit();
    }
    #[getter]
    /// Per-atom feature contribution counts, or None when not collected.
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
    #[getter]
    /// Contributing fingerprint bit identifiers for each atom, or None when not collected.
    fn atom_to_bits(&self) -> Option<Vec<Vec<u64>>> {
        self.inner.atom_to_bits().map(<[Vec<u64>]>::to_vec)
    }
    #[getter]
    /// Bit identifiers mapped to native feature descriptors (center/radius or atom-pair endpoints, depending on the generator), or None when not collected.
    fn bit_info_map(&self) -> Option<BTreeMap<u64, Vec<(u32, u32)>>> {
        self.inner.bit_info_map().cloned()
    }
    #[getter]
    /// Bit identifiers mapped to contributing native paths, or None when not collected.
    fn bit_paths(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.bit_paths().cloned()
    }
    #[getter]
    /// Bit identifiers mapped to participating atom sets, or None when not collected.
    fn atoms_per_bit(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.atoms_per_bit().cloned()
    }
}
