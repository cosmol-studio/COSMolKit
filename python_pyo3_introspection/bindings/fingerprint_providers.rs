use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

ck::project_configuration!(ck::AtomPairAtomInvariantsGenerator, methods {
    fn info_string(&self) -> String { self.inner.info_string() }
    fn to_json(&self) -> String { self.inner.to_json() }
});

ck::project_configuration!(ck::MorganBondInvariantsGenerator);

#[pyclass(module = "cosmolkit")]
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
    fn new() -> Self {
        Self::construct()
    }
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::default(),
        }
    }
    fn allocate_atom_counts(&mut self) {
        self.inner.allocate_atom_counts();
    }
    fn allocate_atom_to_bits(&mut self) {
        self.inner.allocate_atom_to_bits();
    }
    fn allocate_bit_info_map(&mut self) {
        self.inner.allocate_bit_info_map();
    }
    fn allocate_bit_paths(&mut self) {
        self.inner.allocate_bit_paths();
    }
    fn allocate_atoms_per_bit(&mut self) {
        self.inner.allocate_atoms_per_bit();
    }
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
    fn atom_to_bits(&self) -> Option<Vec<Vec<u64>>> {
        self.inner.atom_to_bits().map(<[Vec<u64>]>::to_vec)
    }
    fn bit_info_map(&self) -> Option<BTreeMap<u64, Vec<(u32, u32)>>> {
        self.inner.bit_info_map().cloned()
    }
    fn bit_paths(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.bit_paths().cloned()
    }
    fn atoms_per_bit(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.atoms_per_bit().cloned()
    }
}
