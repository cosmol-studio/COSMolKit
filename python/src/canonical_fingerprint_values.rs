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

/// Detached invariant selection, interpreted by the canonical Rust owner.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganInvariants {
    pub(crate) inner: ck::MorganInvariants,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganInvariants {
    #[staticmethod]
    fn connectivity() -> Self {
        Self {
            inner: ck::MorganInvariants::Connectivity,
        }
    }

    #[staticmethod]
    fn features() -> Self {
        Self {
            inner: ck::MorganInvariants::Features,
        }
    }
}

/// Owned call options; optional empty vectors retain their source meaning.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganFingerprintParams {
    pub(crate) inner: ck::MorganFingerprintParams,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganFingerprintParams {
    #[new]
    #[pyo3(signature = (*, generator=None, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None, conformer_id=-1, invariants=None))]
    fn new(
        generator: Option<&MorganParams>,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
        invariants: Option<&MorganInvariants>,
    ) -> Self {
        Self {
            inner: ck::MorganFingerprintParams {
                generator: generator
                    .map(|value| value.inner.clone())
                    .unwrap_or_default(),
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                invariants: invariants
                    .map(|value| value.inner.clone())
                    .unwrap_or(ck::MorganInvariants::Connectivity),
            },
        }
    }

    #[getter]
    fn generator(&self) -> MorganParams {
        MorganParams {
            inner: self.inner.generator.clone(),
        }
    }

    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[getter]
    fn invariants(&self) -> MorganInvariants {
        MorganInvariants {
            inner: self.inner.invariants.clone(),
        }
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
    #[new]
    fn new() -> Self {
        Self::default()
    }

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
    fn __repr__(&self) -> String {
        format!(
            "AdditionalOutput(atom_to_bits={}, bit_info_map={}, bit_paths={}, atom_counts={}, atoms_per_bit={})",
            self.inner.atom_to_bits().is_some(),
            self.inner.bit_info_map().is_some(),
            self.inner.bit_paths().is_some(),
            self.inner.atom_counts().is_some(),
            self.inner.atoms_per_bit().is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairParams {
    pub(crate) inner: ck::AtomPairParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairParams {
    #[new]
    #[pyo3(signature = (*, min_distance=1, max_distance=30, include_chirality=false, use_2d=true, count_simulation=true, fp_size=2048, bits_per_feature=1, count_bounds=None))]
    fn new(
        min_distance: u32,
        max_distance: u32,
        include_chirality: bool,
        use_2d: bool,
        count_simulation: bool,
        fp_size: u32,
        bits_per_feature: u32,
        count_bounds: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::AtomPairParams {
                min_distance,
                max_distance,
                include_chirality,
                use_2d,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds: count_bounds
                    .unwrap_or_else(|| ck::AtomPairParams::default().count_bounds),
            },
        }
    }
    #[getter]
    fn min_distance(&self) -> u32 {
        self.inner.min_distance
    }
    #[getter]
    fn max_distance(&self) -> u32 {
        self.inner.max_distance
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[getter]
    fn use_2d(&self) -> bool {
        self.inner.use_2d
    }
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairAtomInvariantsGenerator {
    pub(crate) inner: ck::AtomPairAtomInvariantsGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairAtomInvariantsGenerator {
    #[new]
    #[pyo3(signature = (*, include_chirality=false, topological_torsion_correction=false))]
    fn new(include_chirality: bool, topological_torsion_correction: bool) -> Self {
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
    #[getter]
    fn topological_torsion_correction(&self) -> bool {
        self.inner.topological_torsion_correction
    }
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    fn __repr__(&self) -> String {
        format!(
            "AtomPairAtomInvariantsGenerator({})",
            self.inner.info_string()
        )
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairFingerprintParams {
    pub(crate) inner: ck::AtomPairFingerprintParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairFingerprintParams {
    #[new]
    #[pyo3(signature = (*, generator=None, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None, conformer_id=-1, atom_invariants_generator=None, use_legacy_stereo_perception=true))]
    fn new(
        generator: Option<&AtomPairParams>,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
        atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
        use_legacy_stereo_perception: bool,
    ) -> Self {
        Self {
            inner: ck::AtomPairFingerprintParams {
                generator: generator
                    .map(|value| value.inner.clone())
                    .unwrap_or_default(),
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator: atom_invariants_generator.map(|value| value.inner),
                use_legacy_stereo_perception,
            },
        }
    }
    #[getter]
    fn generator(&self) -> AtomPairParams {
        AtomPairParams {
            inner: self.inner.generator.clone(),
        }
    }
    #[getter]
    fn atom_invariants_generator(&self) -> Option<AtomPairAtomInvariantsGenerator> {
        self.inner
            .atom_invariants_generator
            .map(|inner| AtomPairAtomInvariantsGenerator { inner })
    }
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[getter]
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalTorsionParams {
    pub(crate) inner: ck::TopologicalTorsionParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionParams {
    #[new]
    #[pyo3(signature = (*, torsion_atom_count=4, only_shortest_paths=false, include_chirality=false, count_simulation=true, fp_size=2048, bits_per_feature=1, count_bounds=None))]
    fn new(
        torsion_atom_count: u32,
        only_shortest_paths: bool,
        include_chirality: bool,
        count_simulation: bool,
        fp_size: u32,
        bits_per_feature: u32,
        count_bounds: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::TopologicalTorsionParams {
                torsion_atom_count,
                only_shortest_paths,
                include_chirality,
                count_simulation,
                fp_size,
                bits_per_feature,
                count_bounds: count_bounds
                    .unwrap_or_else(|| ck::TopologicalTorsionParams::default().count_bounds),
            },
        }
    }
    #[getter]
    fn torsion_atom_count(&self) -> u32 {
        self.inner.torsion_atom_count
    }
    #[getter]
    fn only_shortest_paths(&self) -> bool {
        self.inner.only_shortest_paths
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalTorsionFingerprintParams {
    pub(crate) inner: ck::TopologicalTorsionFingerprintParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionFingerprintParams {
    #[new]
    #[pyo3(signature = (*, generator=None, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None, custom_bond_invariants=None, conformer_id=-1, atom_invariants_generator=None, use_legacy_stereo_perception=true))]
    fn new(
        generator: Option<&TopologicalTorsionParams>,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
        atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
        use_legacy_stereo_perception: bool,
    ) -> Self {
        Self {
            inner: ck::TopologicalTorsionFingerprintParams {
                generator: generator
                    .map(|value| value.inner.clone())
                    .unwrap_or_default(),
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                atom_invariants_generator: atom_invariants_generator.map(|value| value.inner),
                use_legacy_stereo_perception,
            },
        }
    }
    #[getter]
    fn generator(&self) -> TopologicalTorsionParams {
        TopologicalTorsionParams {
            inner: self.inner.generator.clone(),
        }
    }
    #[getter]
    fn atom_invariants_generator(&self) -> Option<AtomPairAtomInvariantsGenerator> {
        self.inner
            .atom_invariants_generator
            .map(|inner| AtomPairAtomInvariantsGenerator { inner })
    }
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    #[getter]
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<TopologicalTorsionParams>()?;
    module.add_class::<TopologicalTorsionFingerprintParams>()?;
    module.add_class::<AtomPairParams>()?;
    module.add_class::<AtomPairFingerprintParams>()?;
    module.add_class::<AtomPairAtomInvariantsGenerator>()?;
    module.add_class::<MorganParams>()?;
    module.add_class::<MorganInvariants>()?;
    module.add_class::<MorganFingerprintParams>()?;
    module.add_class::<AdditionalOutput>()?;
    Ok(())
}
