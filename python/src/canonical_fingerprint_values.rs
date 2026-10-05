//! Detached fingerprint value projections through the public facade.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganParams {
    pub(crate) inner: ck::MorganParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganParams {
    #[new]
    #[pyo3(signature = (*, radius=3, include_chirality=false, use_bond_types=true, include_ring_membership=true, only_nonzero_invariants=false, include_redundant_environments=false, fp_size=2048, count_simulation=false, count_bounds=None, bits_per_feature=1))]
    fn new(
        radius: u32,
        include_chirality: bool,
        use_bond_types: bool,
        include_ring_membership: bool,
        only_nonzero_invariants: bool,
        include_redundant_environments: bool,
        fp_size: u32,
        count_simulation: bool,
        count_bounds: Option<Vec<u32>>,
        bits_per_feature: u32,
    ) -> Self {
        // Configuration transport only: validation belongs to canonical callers.
        // An explicit empty vector remains distinct from the source default.
        Self {
            inner: ck::MorganParams {
                radius,
                include_chirality,
                use_bond_types,
                include_ring_membership,
                only_nonzero_invariants,
                include_redundant_environments,
                fp_size,
                count_simulation,
                count_bounds: count_bounds
                    .unwrap_or_else(|| ck::MorganParams::default().count_bounds),
                bits_per_feature,
            },
        }
    }

    #[getter]
    fn radius(&self) -> u32 {
        self.inner.radius
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[getter]
    fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types
    }
    #[getter]
    fn include_ring_membership(&self) -> bool {
        self.inner.include_ring_membership
    }
    #[getter]
    fn only_nonzero_invariants(&self) -> bool {
        self.inner.only_nonzero_invariants
    }
    #[getter]
    fn include_redundant_environments(&self) -> bool {
        self.inner.include_redundant_environments
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
}

/// One uniquely owned canonical metadata value, with detached read results.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct AdditionalOutput {
    pub(crate) inner: ck::AdditionalOutput,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AdditionalOutput {
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: ck::AdditionalOutput::default(),
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

    // Only borrowed results are copied for Python ownership. The canonical
    // AdditionalOutput stays unique, and Option preserves unallocated states.
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

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<MorganParams>()?;
    module.add_class::<AdditionalOutput>()?;
    Ok(())
}
