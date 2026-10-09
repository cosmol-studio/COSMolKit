//! Detached fingerprint value projections through the public facade.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganParams {
    pub(crate) inner: ck::MorganParams,
}

#[cosmolkit_macros::python_configuration]
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
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_values::fingerprint_json_pyerr(py, error))
    }
}

/// Detached invariant selection, interpreted by the canonical Rust owner.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(PartialEq)]
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
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganFingerprintParams {
    pub(crate) inner: ck::MorganFingerprintParams,
}

#[cosmolkit_macros::python_configuration]
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
pub(crate) struct FingerprintAdditionalOutput {
    pub(crate) inner: ck::FingerprintAdditionalOutput,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FingerprintAdditionalOutput {
    #[new]
    fn py_new() -> Self {
        Self::new()
    }
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::new(),
        }
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

    // Only borrowed results are copied for Python ownership. The canonical
    // FingerprintAdditionalOutput stays unique, and Option preserves unallocated states.
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
            "FingerprintAdditionalOutput(atom_to_bits={}, bit_info_map={}, bit_paths={}, atom_counts={}, atoms_per_bit={})",
            self.inner.atom_to_bits().is_some(),
            self.inner.bit_info_map().is_some(),
            self.inner.bit_paths().is_some(),
            self.inner.atom_counts().is_some(),
            self.inner.atoms_per_bit().is_some()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct AtomPairParams {
    pub(crate) inner: ck::AtomPairParams,
}
#[cosmolkit_macros::python_configuration]
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
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_values::fingerprint_json_pyerr(py, error))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct AtomPairAtomInvariantsGenerator {
    pub(crate) inner: ck::AtomPairAtomInvariantsGenerator,
}
#[cosmolkit_macros::python_configuration]
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
#[pyclass(module = "cosmolkit")]
pub(crate) struct AtomPairFingerprintParams {
    pub(crate) inner: ck::AtomPairFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
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
#[pyclass(module = "cosmolkit")]
pub(crate) struct TopologicalTorsionParams {
    pub(crate) inner: ck::TopologicalTorsionParams,
}
#[cosmolkit_macros::python_configuration]
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
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_values::fingerprint_json_pyerr(py, error))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TopologicalTorsionFingerprintParams {
    pub(crate) inner: ck::TopologicalTorsionFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
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
    module.add_class::<TopologicalTorsionFingerprintGenerator>()?;
    module.add_class::<TopologicalTorsionSettings>()?;
    module.add_class::<TopologicalTorsionCallParams>()?;
    module.add_class::<TopologicalTorsionParams>()?;
    module.add_class::<AtomCodeExplanation>()?;
    module.add_class::<AtomPairsParameters>()?;
    module.add_class::<AtomPairAtomCodeResult>()?;
    module.add(
        "AtomCodeExplanationError",
        module.py().get_type::<AtomCodeExplanationError>(),
    )?;
    module.add_class::<LegacyTopologicalTorsionParams>()?;
    module.add_class::<TopologicalTorsionFingerprintParams>()?;
    module.add_class::<AtomPairParams>()?;
    module.add_class::<AtomPairFingerprintParams>()?;
    module.add_class::<AtomPairAtomInvariantsGenerator>()?;
    module.add_class::<MorganParams>()?;
    module.add_class::<MorganAtomInvariantsGenerator>()?;
    module.add_class::<MorganBondInvariantsGenerator>()?;
    module.add_class::<MorganFingerprintGenerator>()?;
    module.add_class::<MorganSettings>()?;
    module.add_class::<MorganCallParams>()?;
    module.add_class::<MorganInvariants>()?;
    module.add_class::<MorganFingerprintParams>()?;
    module.add_class::<FingerprintAdditionalOutput>()?;
    Ok(())
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TopologicalTorsionCallParams {
    pub(crate) inner: ck::TopologicalTorsionCallParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionCallParams {
    #[new]
    #[pyo3(signature=(*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,conformer_id=-1,use_legacy_stereo_perception=true))]
    fn new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
        use_legacy_stereo_perception: bool,
    ) -> Self {
        Self {
            inner: ck::TopologicalTorsionCallParams {
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
                use_legacy_stereo_perception,
            },
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
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalTorsionFingerprintGenerator {
    pub(crate) inner: ck::TopologicalTorsionFingerprintGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionFingerprintGenerator {
    #[new]
    #[pyo3(signature=(*,params=None,atom_invariants_generator=None))]
    fn py_new(
        py: Python<'_>,
        params: Option<&TopologicalTorsionParams>,
        atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
    ) -> PyResult<Self> {
        Self::new(py, params, atom_invariants_generator)
    }
    #[staticmethod]
    #[pyo3(signature=(*,params=None,atom_invariants_generator=None))]
    fn new(
        py: Python<'_>,
        params: Option<&TopologicalTorsionParams>,
        atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
    ) -> PyResult<Self> {
        ck::TopologicalTorsionFingerprintGenerator::new(
            params.map(|p| &p.inner),
            atom_invariants_generator.map(|i| i.inner),
        )
        .map(|inner| Self { inner })
        .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[staticmethod]
    fn from_json(py: Python<'_>, json: &str) -> PyResult<Self> {
        ck::TopologicalTorsionFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings {
            inner: self.inner.settings(),
        }
    }
    fn info_string(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .info_string()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    fn to_json(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "TopologicalTorsionFingerprintGenerator({})",
            self.info_string(py)?
        ))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::Fingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        self.inner
            .fingerprints(&rows, num_threads)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| crate::canonical_values::Fingerprint { inner }))
                    .collect()
            })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn sparse_fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseBitFingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        self.inner
            .sparse_fingerprints(&rows, num_threads)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseBitFingerprint { inner })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseCountFingerprint32>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        self.inner
            .counts(&rows, num_threads)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseCountFingerprint32 {
                            inner,
                        })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn sparse_counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseCountFingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        self.inner
            .sparse_counts(&rows, num_threads)
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseCountFingerprint { inner })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TopologicalTorsionSettings {
    inner: ck::TopologicalTorsionSettings,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionSettings {
    #[pyo3(name = "set_torsion_atom_count")]
    fn set_torsion_atom_count_method(&mut self, value: u32) -> PyResult<()> {
        self.set_torsion_atom_count(value)
    }
    #[pyo3(name = "set_only_shortest_paths")]
    fn set_only_shortest_paths_method(&mut self, value: bool) -> PyResult<()> {
        self.set_only_shortest_paths(value)
    }
    #[pyo3(name = "set_include_chirality")]
    fn set_include_chirality_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_chirality(value)
    }
    #[pyo3(name = "set_count_simulation")]
    fn set_count_simulation_method(&mut self, value: bool) -> PyResult<()> {
        self.set_count_simulation(value)
    }
    #[pyo3(name = "set_fp_size")]
    fn set_fp_size_method(&mut self, value: u32) -> PyResult<()> {
        self.set_fp_size(value)
    }
    #[pyo3(name = "set_bits_per_feature")]
    fn set_bits_per_feature_method(&mut self, value: u32) -> PyResult<()> {
        self.set_bits_per_feature(value)
    }
    #[pyo3(name = "set_count_bounds")]
    fn set_count_bounds_method(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.set_count_bounds(value)
    }
    #[getter]
    fn torsion_atom_count(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .torsion_atom_count()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_torsion_atom_count(&mut self, value: u32) -> PyResult<()> {
        self.inner.set_torsion_atom_count(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn only_shortest_paths(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .only_shortest_paths()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_only_shortest_paths(&mut self, value: bool) -> PyResult<()> {
        self.inner.set_only_shortest_paths(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn include_chirality(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .include_chirality()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_include_chirality(&mut self, value: bool) -> PyResult<()> {
        self.inner.set_include_chirality(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn count_simulation(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .count_simulation()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_count_simulation(&mut self, value: bool) -> PyResult<()> {
        self.inner.set_count_simulation(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn fp_size(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .fp_size()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_fp_size(&mut self, value: u32) -> PyResult<()> {
        self.inner.set_fp_size(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn bits_per_feature(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .bits_per_feature()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_bits_per_feature(&mut self, value: u32) -> PyResult<()> {
        self.inner.set_bits_per_feature(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    #[getter]
    fn count_bounds(&self, py: Python<'_>) -> PyResult<Vec<u32>> {
        self.inner
            .count_bounds()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    #[setter]
    fn set_count_bounds(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.inner.set_count_bounds(value).map_err(|e| {
            Python::attach(|py| crate::canonical_values::topological_torsion_pyerr(py, e))
        })
    }
    fn params(&self, py: Python<'_>) -> PyResult<TopologicalTorsionParams> {
        self.inner
            .params()
            .map(|inner| TopologicalTorsionParams { inner })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        let p = self
            .inner
            .params()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))?;
        Ok(format!(
            "TopologicalTorsionSettings(torsion_atom_count={}, fp_size={}, count_simulation={}, include_chirality={}, only_shortest_paths={})",
            p.torsion_atom_count,
            p.fp_size,
            p.count_simulation,
            p.include_chirality,
            p.only_shortest_paths
        ))
    }
}

/// Owned, immutable legacy torsion settings and call selections.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct LegacyTopologicalTorsionParams {
    pub(crate) inner: ck::LegacyTopologicalTorsionParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LegacyTopologicalTorsionParams {
    #[new]
    #[pyo3(signature = (*, torsion_atom_count=4, include_chirality=false, fp_size=2048, bits_per_entry=4, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None))]
    fn py_new(
        torsion_atom_count: u32,
        include_chirality: bool,
        fp_size: u32,
        bits_per_entry: u32,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self::new(
            torsion_atom_count,
            include_chirality,
            fp_size,
            bits_per_entry,
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
        )
    }
    #[staticmethod]
    #[pyo3(signature = (*, torsion_atom_count=4, include_chirality=false, fp_size=2048, bits_per_entry=4, from_atoms=None, ignore_atoms=None, custom_atom_invariants=None))]
    fn new(
        torsion_atom_count: u32,
        include_chirality: bool,
        fp_size: u32,
        bits_per_entry: u32,
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::LegacyTopologicalTorsionParams::new(
                torsion_atom_count,
                include_chirality,
                fp_size,
                bits_per_entry,
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
            ),
        }
    }
    #[getter]
    fn torsion_atom_count(&self) -> u32 {
        self.inner.torsion_atom_count
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[getter]
    fn bits_per_entry(&self) -> u32 {
        self.inner.bits_per_entry
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
}

pyo3::create_exception!(
    cosmolkit,
    AtomCodeExplanationError,
    pyo3::exceptions::PyKeyError
);
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairAtomCodeResult {
    pub(crate) inner: ck::AtomPairAtomCodeResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairAtomCodeResult {
    #[getter]
    fn code(&self) -> u32 {
        self.inner.code
    }
    #[getter]
    fn molecule(&self) -> crate::drawing_binding::Molecule {
        crate::drawing_binding::Molecule {
            inner: self.inner.molecule.clone(),
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomCodeExplanation {
    pub(crate) inner: ck::AtomCodeExplanation,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomCodeExplanation {
    #[staticmethod]
    #[pyo3(signature = (code, branch_subtract=0, include_chirality=false))]
    fn from_code(
        code: u64,
        branch_subtract: i64,
        include_chirality: bool,
        py: Python<'_>,
    ) -> PyResult<Self> {
        ck::AtomCodeExplanation::from_code(code, branch_subtract, include_chirality)
            .map(|inner| Self { inner })
            .map_err(|error| {
                let ck::AtomCodeExplanationError::UnknownChirality { code } = error;
                let exception = AtomCodeExplanationError::new_err(code);
                let object = exception.value(py);
                for (name, value) in [("domain", "Fingerprint"), ("kind", "UnknownChirality")] {
                    if let Err(attribute_error) = object.setattr(name, value) {
                        return attribute_error;
                    }
                }
                if let Err(attribute_error) = object.setattr("code", code) {
                    return attribute_error;
                }
                exception
            })
    }
    fn symbol(&self) -> &'static str {
        self.inner.symbol()
    }
    fn branch_count(&self) -> u32 {
        self.inner.branch_count()
    }
    fn pi_electrons(&self) -> u32 {
        self.inner.pi_electrons()
    }
    fn chirality(&self) -> Option<&'static str> {
        self.inner.chirality()
    }
}

// Registered persistent Morgan projections. All computation delegates to ck.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganAtomInvariantsGenerator {
    pub(crate) inner: ck::MorganAtomInvariantsGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganAtomInvariantsGenerator {
    #[staticmethod]
    #[pyo3(signature=(include_ring_membership=true))]
    fn connectivity(include_ring_membership: bool) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::connectivity(include_ring_membership),
        }
    }
    #[staticmethod]
    #[pyo3(signature=(patterns=None))]
    fn features(
        py: Python<'_>,
        patterns: Option<Vec<Py<crate::canonical_search::QueryGraph>>>,
    ) -> PyResult<Self> {
        let patterns = patterns
            .map(|patterns| {
                patterns
                    .iter()
                    .map(|query| {
                        query
                            .try_borrow(py)
                            .map(|q| q.inner.clone())
                            .map_err(PyErr::from)
                    })
                    .collect::<PyResult<Vec<_>>>()
            })
            .transpose()?;
        Ok(Self {
            inner: ck::MorganAtomInvariantsGenerator::features(patterns),
        })
    }
    #[staticmethod]
    fn atom_pair(generator: &AtomPairAtomInvariantsGenerator) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::atom_pair(generator.inner),
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganBondInvariantsGenerator {
    pub(crate) inner: ck::MorganBondInvariantsGenerator,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganBondInvariantsGenerator {
    #[new]
    #[pyo3(signature=(*,use_bond_types=true,include_chirality=false))]
    fn py_new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self::new(use_bond_types, include_chirality)
    }
    #[staticmethod]
    #[pyo3(signature=(*,use_bond_types=true,include_chirality=false))]
    fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self {
            inner: ck::MorganBondInvariantsGenerator::new(use_bond_types, include_chirality),
        }
    }
    #[getter]
    fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types()
    }
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganCallParams {
    pub(crate) inner: ck::MorganCallParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganCallParams {
    #[new]
    #[pyo3(signature=(*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,conformer_id=-1))]
    fn py_new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
    ) -> Self {
        Self::new(
            from_atoms,
            ignore_atoms,
            custom_atom_invariants,
            custom_bond_invariants,
            conformer_id,
        )
    }
    #[staticmethod]
    #[pyo3(signature=(*,from_atoms=None,ignore_atoms=None,custom_atom_invariants=None,custom_bond_invariants=None,conformer_id=-1))]
    fn new(
        from_atoms: Option<Vec<u32>>,
        ignore_atoms: Option<Vec<u32>>,
        custom_atom_invariants: Option<Vec<u32>>,
        custom_bond_invariants: Option<Vec<u32>>,
        conformer_id: i32,
    ) -> Self {
        Self {
            inner: ck::MorganCallParams::new(
                from_atoms,
                ignore_atoms,
                custom_atom_invariants,
                custom_bond_invariants,
                conformer_id,
            ),
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
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganFingerprintGenerator {
    pub(crate) inner: ck::MorganFingerprintGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganFingerprintGenerator {
    #[new]
    #[pyo3(signature=(*,params=None,atom_invariants=None,bond_invariants=None))]
    fn py_new(
        py: Python<'_>,
        params: Option<&MorganParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> PyResult<Self> {
        Self::new(py, params, atom_invariants, bond_invariants)
    }
    #[staticmethod]
    #[pyo3(signature=(*,params=None,atom_invariants=None,bond_invariants=None))]
    fn new(
        py: Python<'_>,
        params: Option<&MorganParams>,
        atom_invariants: Option<&MorganAtomInvariantsGenerator>,
        bond_invariants: Option<&MorganBondInvariantsGenerator>,
    ) -> PyResult<Self> {
        ck::MorganFingerprintGenerator::new(
            params.map(|p| &p.inner),
            atom_invariants.map(|p| &p.inner),
            bond_invariants.map(|p| &p.inner),
        )
        .map(|inner| Self { inner })
        .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[staticmethod]
    fn from_json(py: Python<'_>, json: &str) -> PyResult<Self> {
        ck::MorganFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    fn settings(&self) -> MorganSettings {
        MorganSettings {
            inner: self.inner.settings(),
        }
    }
    fn info_string(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .info_string()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    fn to_json(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_json()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    fn __repr__(&self, py: Python<'_>) -> PyResult<String> {
        Ok(format!(
            "MorganFingerprintGenerator({})",
            self.info_string(py)?
        ))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::Fingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        // RDKit❗✔️:   {
        // RDKit❗✔️:     NOGIL gil;
        // RDKit❗✔️:     fps = std::move(func(tmols, numThreads));
        // RDKit❗✔️:   }
        // Borrow guards live until all source workers join; detach invokes only
        // Rust facade types and keeps every None slot in its original position.
        py.detach(|| self.inner.fingerprints(&rows, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| value.map(|inner| crate::canonical_values::Fingerprint { inner }))
                    .collect()
            })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseCountFingerprint32>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        // RDKit❗✔️:   {
        // RDKit❗✔️:     NOGIL gil;
        // RDKit❗✔️:     fps = std::move(func(tmols, numThreads));
        // RDKit❗✔️:   }
        // Borrow guards live until all source workers join; detach invokes only
        // Rust facade types and keeps every None slot in its original position.
        py.detach(|| self.inner.counts(&rows, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseCountFingerprint32 {
                            inner,
                        })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn sparse_fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseBitFingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        // RDKit❗✔️:   {
        // RDKit❗✔️:     NOGIL gil;
        // RDKit❗✔️:     fps = std::move(func(tmols, numThreads));
        // RDKit❗✔️:   }
        // Borrow guards live until all source workers join; detach invokes only
        // Rust facade types and keeps every None slot in its original position.
        py.detach(|| self.inner.sparse_fingerprints(&rows, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseBitFingerprint { inner })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn sparse_counts(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<Vec<Option<crate::canonical_values::SparseCountFingerprint>>> {
        let borrowed = molecules
            .iter()
            .map(|m| m.as_ref().map(|m| m.try_borrow(py)).transpose())
            .collect::<Result<Vec<_>, _>>()?;
        let rows = borrowed
            .iter()
            .map(|m| m.as_ref().map(|m| &m.inner))
            .collect::<Vec<_>>();
        // RDKit❗✔️:   {
        // RDKit❗✔️:     NOGIL gil;
        // RDKit❗✔️:     fps = std::move(func(tmols, numThreads));
        // RDKit❗✔️:   }
        // Borrow guards live until all source workers join; detach invokes only
        // Rust facade types and keeps every None slot in its original position.
        py.detach(|| self.inner.sparse_counts(&rows, num_threads))
            .map(|values| {
                values
                    .into_iter()
                    .map(|value| {
                        value.map(|inner| crate::canonical_values::SparseCountFingerprint { inner })
                    })
                    .collect()
            })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganSettings {
    inner: ck::MorganSettings,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganSettings {
    #[pyo3(name = "set_radius")]
    fn set_radius_method(&mut self, value: u32) -> PyResult<()> {
        self.set_radius(value)
    }
    #[pyo3(name = "set_only_nonzero_invariants")]
    fn set_only_nonzero_invariants_method(&mut self, value: bool) -> PyResult<()> {
        self.set_only_nonzero_invariants(value)
    }
    #[pyo3(name = "set_include_redundant_environments")]
    fn set_include_redundant_environments_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_redundant_environments(value)
    }
    #[pyo3(name = "set_include_chirality")]
    fn set_include_chirality_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_chirality(value)
    }
    #[pyo3(name = "set_count_simulation")]
    fn set_count_simulation_method(&mut self, value: bool) -> PyResult<()> {
        self.set_count_simulation(value)
    }
    #[pyo3(name = "set_fp_size")]
    fn set_fp_size_method(&mut self, value: u32) -> PyResult<()> {
        self.set_fp_size(value)
    }
    #[pyo3(name = "set_bits_per_feature")]
    fn set_bits_per_feature_method(&mut self, value: u32) -> PyResult<()> {
        self.set_bits_per_feature(value)
    }
    #[pyo3(name = "set_count_bounds")]
    fn set_count_bounds_method(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.set_count_bounds(value)
    }
    #[getter]
    fn radius(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .radius()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_radius(&mut self, value: u32) -> PyResult<()> {
        self.inner
            .set_radius(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn only_nonzero_invariants(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .only_nonzero_invariants()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_only_nonzero_invariants(&mut self, value: bool) -> PyResult<()> {
        self.inner
            .set_only_nonzero_invariants(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn include_redundant_environments(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .include_redundant_environments()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_include_redundant_environments(&mut self, value: bool) -> PyResult<()> {
        self.inner
            .set_include_redundant_environments(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn include_chirality(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .include_chirality()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_include_chirality(&mut self, value: bool) -> PyResult<()> {
        self.inner
            .set_include_chirality(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn count_simulation(&self, py: Python<'_>) -> PyResult<bool> {
        self.inner
            .count_simulation()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_count_simulation(&mut self, value: bool) -> PyResult<()> {
        self.inner
            .set_count_simulation(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn fp_size(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .fp_size()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_fp_size(&mut self, value: u32) -> PyResult<()> {
        self.inner
            .set_fp_size(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn bits_per_feature(&self, py: Python<'_>) -> PyResult<u32> {
        self.inner
            .bits_per_feature()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_bits_per_feature(&mut self, value: u32) -> PyResult<()> {
        self.inner
            .set_bits_per_feature(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    #[getter]
    fn count_bounds(&self, py: Python<'_>) -> PyResult<Vec<u32>> {
        self.inner
            .count_bounds()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    #[setter]
    fn set_count_bounds(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.inner
            .set_count_bounds(value)
            .map_err(|e| Python::attach(|py| crate::canonical_values::morgan_pyerr(py, e)))
    }
    fn params(&self, py: Python<'_>) -> PyResult<MorganParams> {
        self.inner
            .params()
            .map(|inner| MorganParams { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairsParameters;
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairsParameters {
    #[staticmethod]
    fn version() -> &'static str {
        ck::AtomPairsParameters::version()
    }
    #[staticmethod]
    fn num_type_bits() -> u32 {
        ck::AtomPairsParameters::num_type_bits()
    }
    #[staticmethod]
    fn num_pi_bits() -> u32 {
        ck::AtomPairsParameters::num_pi_bits()
    }
    #[staticmethod]
    fn num_branch_bits() -> u32 {
        ck::AtomPairsParameters::num_branch_bits()
    }
    #[staticmethod]
    fn num_chiral_bits() -> u32 {
        ck::AtomPairsParameters::num_chiral_bits()
    }
    #[staticmethod]
    fn code_size() -> u32 {
        ck::AtomPairsParameters::code_size()
    }
    #[staticmethod]
    fn num_path_bits() -> u32 {
        ck::AtomPairsParameters::num_path_bits()
    }
    #[staticmethod]
    fn max_path_length() -> u32 {
        ck::AtomPairsParameters::max_path_length()
    }
    #[staticmethod]
    fn num_atom_pair_fingerprint_bits() -> u32 {
        ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
    }
    #[staticmethod]
    fn atom_types() -> Vec<u32> {
        ck::AtomPairsParameters::atom_types()
    }
}
