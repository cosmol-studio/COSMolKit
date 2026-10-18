use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(module = "fixed_projection._native")]
pub struct AtomPairAtomInvariantsGenerator {
    pub(crate) inner: ck::AtomPairAtomInvariantsGenerator,
}

#[pymethods]
impl AtomPairAtomInvariantsGenerator {
    #[new]
    #[pyo3(signature = (*, include_chirality=ck::AtomPairAtomInvariantsGenerator::default().include_chirality, topological_torsion_correction=ck::AtomPairAtomInvariantsGenerator::default().topological_torsion_correction))]
    fn construct(include_chirality: bool, topological_torsion_correction: bool) -> Self {
        Self {
            inner: ck::AtomPairAtomInvariantsGenerator {
                include_chirality,
                topological_torsion_correction,
            },
        }
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[setter]
    fn set_include_chirality(&mut self, value: bool) {
        self.inner.include_chirality = value;
    }
    #[getter]
    fn topological_torsion_correction(&self) -> bool {
        self.inner.topological_torsion_correction
    }
    #[setter]
    fn set_topological_torsion_correction(&mut self, value: bool) {
        self.inner.topological_torsion_correction = value;
    }
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct MorganBondInvariantsGenerator {
    pub(crate) inner: ck::MorganBondInvariantsGenerator,
}

#[pymethods]
impl MorganBondInvariantsGenerator {
    #[new]
    #[pyo3(signature = (*, use_bond_types=ck::MorganBondInvariantsGenerator::default().use_bond_types(), include_chirality=ck::MorganBondInvariantsGenerator::default().include_chirality()))]
    fn construct(use_bond_types: bool, include_chirality: bool) -> Self {
        Self {
            inner: ck::MorganBondInvariantsGenerator::new(use_bond_types, include_chirality),
        }
    }
    #[staticmethod]
    #[pyo3(signature = (use_bond_types=ck::MorganBondInvariantsGenerator::default().use_bond_types(), include_chirality=ck::MorganBondInvariantsGenerator::default().include_chirality()))]
    fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self::construct(use_bond_types, include_chirality)
    }
    #[getter]
    fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types()
    }
    #[setter]
    fn set_use_bond_types(&mut self, value: bool) {
        self.inner = ck::MorganBondInvariantsGenerator::new(value, self.inner.include_chirality());
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality()
    }
    #[setter]
    fn set_include_chirality(&mut self, value: bool) {
        self.inner = ck::MorganBondInvariantsGenerator::new(self.inner.use_bond_types(), value);
    }
}

#[pyclass(module = "fixed_projection._native")]
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
