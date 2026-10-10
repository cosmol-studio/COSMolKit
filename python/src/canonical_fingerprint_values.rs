//! Detached fingerprint value projections through the public facade.

use ::cosmolkit as ck;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;

/// Writable configuration for Morgan fingerprint generation.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct MorganParams {
    pub(crate) inner: ck::MorganParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganParams {
    /// Configure Morgan fingerprint generation; omitted fields use the defaults shown in the signature.
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

    /// Morgan environment radius, measured in graph bonds.
    #[getter]
    fn radius(&self) -> u32 {
        self.inner.radius
    }
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    /// Whether bond types contribute to fingerprint invariants.
    #[getter]
    fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types
    }
    /// Whether ring membership contributes to Morgan connectivity invariants.
    #[getter]
    fn include_ring_membership(&self) -> bool {
        self.inner.include_ring_membership
    }
    /// Whether atoms with zero invariants are excluded as environment centers.
    #[getter]
    fn only_nonzero_invariants(&self) -> bool {
        self.inner.only_nonzero_invariants
    }
    /// Whether redundant Morgan environments are retained.
    #[getter]
    fn include_redundant_environments(&self) -> bool {
        self.inner.include_redundant_environments
    }
    /// Number of bins/bits in the folded fingerprint.
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    /// Whether occurrence counts are represented using threshold bits.
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    /// Occurrence-count thresholds used by count simulation.
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
    /// Number of hashed bit positions generated per feature.
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    /// Return a JSON string containing this configuration.
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    /// Return an independent configuration updated from the supplied JSON text; invalid JSON/options raise an error.
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
    /// Select atom invariants based on chemical graph connectivity.
    #[staticmethod]
    fn connectivity() -> Self {
        Self {
            inner: ck::MorganInvariants::Connectivity,
        }
    }

    /// Select pharmacophore/chemical-feature atom invariants for Morgan fingerprinting.
    #[staticmethod]
    fn features() -> Self {
        Self {
            inner: ck::MorganInvariants::Features,
        }
    }
}

/// Owned call options; optional empty vectors retain their source meaning.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct MorganFingerprintParams {
    pub(crate) inner: ck::MorganFingerprintParams,
}

#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganFingerprintParams {
    /// Configure Morgan generator selection and per-call atom/bond invariants; omitted fields use the defaults shown in the signature.
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

    /// Fingerprint generator configuration; nested edits update the parent parameter object.
    #[getter]
    fn generator(&self) -> MorganParams {
        MorganParams {
            inner: self.inner.generator.clone(),
        }
    }

    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    /// Caller-supplied bond invariants in bond-index order.
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
    /// Stored 3D conformer identifier used by this operation; not its position in the conformer list.
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    /// Morgan atom-invariant selection, such as connectivity or feature invariants.
    #[getter]
    fn invariants(&self) -> MorganInvariants {
        MorganInvariants {
            inner: self.inner.invariants.clone(),
        }
    }
}

/// Optional fingerprint metadata collector.
///
/// Allocate the desired atom-count, atom-to-bit or bit-environment sinks before a
/// fingerprint call. Unallocated sinks return None. The collector does not own
/// the returned fingerprint.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct FingerprintAdditionalOutput {
    pub(crate) inner: ck::FingerprintAdditionalOutput,
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl FingerprintAdditionalOutput {
    /// Construct a FingerprintAdditionalOutput value from the supplied inputs.
    #[new]
    fn py_new() -> Self {
        Self::new()
    }
    /// Construct a FingerprintAdditionalOutput value from the supplied inputs.
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::new(),
        }
    }

    /// Construct a FingerprintAdditionalOutput value from the supplied inputs.
    #[staticmethod]
    fn default() -> Self {
        Self {
            inner: ck::FingerprintAdditionalOutput::default(),
        }
    }

    /// Enable collection of atom counts for subsequent fingerprint calls using this collector.
    fn allocate_atom_counts(&mut self) {
        self.inner.allocate_atom_counts();
    }
    /// Enable collection of atom to bits for subsequent fingerprint calls using this collector.
    fn allocate_atom_to_bits(&mut self) {
        self.inner.allocate_atom_to_bits();
    }
    /// Enable collection of bit info map for subsequent fingerprint calls using this collector.
    fn allocate_bit_info_map(&mut self) {
        self.inner.allocate_bit_info_map();
    }
    /// Enable collection of bit paths for subsequent fingerprint calls using this collector.
    fn allocate_bit_paths(&mut self) {
        self.inner.allocate_bit_paths();
    }
    /// Enable collection of atoms per bit for subsequent fingerprint calls using this collector.
    fn allocate_atoms_per_bit(&mut self) {
        self.inner.allocate_atoms_per_bit();
    }

    // Only borrowed results are copied for Python ownership. The canonical
    // FingerprintAdditionalOutput stays unique, and Option preserves unallocated states.
    /// Return atom-indexed feature participation counts, or None if this sink was not allocated.
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts().map(<[u32]>::to_vec)
    }
    /// Return atom-indexed lists of contributed fingerprint bits, or None if this sink was not allocated.
    fn atom_to_bits(&self) -> Option<Vec<Vec<u64>>> {
        self.inner.atom_to_bits().map(<[Vec<u64>]>::to_vec)
    }
    /// Return bit-to-environment mapping with atom/radius pairs, or None if this sink was not allocated.
    fn bit_info_map(&self) -> Option<BTreeMap<u64, Vec<(u32, u32)>>> {
        self.inner.bit_info_map().cloned()
    }
    /// Return bit-to-path mapping of contributing bond-index sequences, or None if this sink was not allocated.
    fn bit_paths(&self) -> Option<BTreeMap<u64, Vec<Vec<i32>>>> {
        self.inner.bit_paths().cloned()
    }
    /// Return bit-to-environment mapping of contributing atom-index sequences, or None if this sink was not allocated.
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

/// Writable configuration for atom-pair fingerprint generation.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct AtomPairParams {
    pub(crate) inner: ck::AtomPairParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairParams {
    /// Configure atom-pair fingerprint generation; omitted fields use the defaults shown in the signature.
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
    /// Minimum atom-pair separation included in fingerprinting.
    #[getter]
    fn min_distance(&self) -> u32 {
        self.inner.min_distance
    }
    /// Maximum separation for atom-pair fingerprints or MCS coordinate matching.
    #[getter]
    fn max_distance(&self) -> u32 {
        self.inner.max_distance
    }
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    /// Whether atom-pair distances are graph distances rather than 3D geometric distances.
    #[getter]
    fn use_2d(&self) -> bool {
        self.inner.use_2d
    }
    /// Whether occurrence counts are represented using threshold bits.
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    /// Number of bins/bits in the folded fingerprint.
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    /// Number of hashed bit positions generated per feature.
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    /// Occurrence-count thresholds used by count simulation.
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    /// Return a JSON string containing this configuration.
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    /// Return an independent configuration updated from the supplied JSON text; invalid JSON/options raise an error.
    fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_values::fingerprint_json_pyerr(py, error))
    }
}
/// Atom-pair atom-code generator with chirality and torsion-correction configuration.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct AtomPairAtomInvariantsGenerator {
    pub(crate) inner: ck::AtomPairAtomInvariantsGenerator,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairAtomInvariantsGenerator {
    /// Construct a AtomPairAtomInvariantsGenerator value from the supplied inputs.
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
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    /// Whether atom-pair atom codes include the topological-torsion branch correction.
    #[getter]
    fn topological_torsion_correction(&self) -> bool {
        self.inner.topological_torsion_correction
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    /// Return a JSON string containing this configuration.
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
/// Writable configuration for atom-pair generator selection and per-call invariants.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct AtomPairFingerprintParams {
    pub(crate) inner: ck::AtomPairFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairFingerprintParams {
    /// Configure atom-pair generator selection and per-call invariants; omitted fields use the defaults shown in the signature.
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
    /// Fingerprint generator configuration; nested edits update the parent parameter object.
    #[getter]
    fn generator(&self) -> AtomPairParams {
        AtomPairParams {
            inner: self.inner.generator.clone(),
        }
    }
    /// Atom-invariant generator used instead of the default invariant calculation.
    #[getter]
    fn atom_invariants_generator(&self) -> Option<AtomPairAtomInvariantsGenerator> {
        self.inner
            .atom_invariants_generator
            .map(|inner| AtomPairAtomInvariantsGenerator { inner })
    }
    /// Stored 3D conformer identifier used by this operation; not its position in the conformer list.
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    /// Whether legacy stereochemistry perception is used by the fingerprint operation.
    #[getter]
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    /// Caller-supplied bond invariants in bond-index order.
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
}
/// Writable configuration for topological-torsion fingerprint generation.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct TopologicalTorsionParams {
    pub(crate) inner: ck::TopologicalTorsionParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionParams {
    /// Configure topological-torsion fingerprint generation; omitted fields use the defaults shown in the signature.
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
    /// Number of atoms in each topological torsion path.
    #[getter]
    fn torsion_atom_count(&self) -> u32 {
        self.inner.torsion_atom_count
    }
    /// Whether only torsion paths that are shortest between their endpoints are retained.
    #[getter]
    fn only_shortest_paths(&self) -> bool {
        self.inner.only_shortest_paths
    }
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    /// Whether occurrence counts are represented using threshold bits.
    #[getter]
    fn count_simulation(&self) -> bool {
        self.inner.count_simulation
    }
    /// Number of bins/bits in the folded fingerprint.
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    /// Number of hashed bit positions generated per feature.
    #[getter]
    fn bits_per_feature(&self) -> u32 {
        self.inner.bits_per_feature
    }
    /// Occurrence-count thresholds used by count simulation.
    #[getter]
    fn count_bounds(&self) -> Vec<u32> {
        self.inner.count_bounds.clone()
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self) -> String {
        self.inner.info_string()
    }
    /// Return a JSON string containing this configuration.
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    /// Return an independent configuration updated from the supplied JSON text; invalid JSON/options raise an error.
    fn with_json(&self, py: Python<'_>, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|error| crate::canonical_values::fingerprint_json_pyerr(py, error))
    }
}
/// Writable configuration for topological-torsion generator selection and per-call invariants.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct TopologicalTorsionFingerprintParams {
    pub(crate) inner: ck::TopologicalTorsionFingerprintParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionFingerprintParams {
    /// Configure topological-torsion generator selection and per-call invariants; omitted fields use the defaults shown in the signature.
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
    /// Fingerprint generator configuration; nested edits update the parent parameter object.
    #[getter]
    fn generator(&self) -> TopologicalTorsionParams {
        TopologicalTorsionParams {
            inner: self.inner.generator.clone(),
        }
    }
    /// Atom-invariant generator used instead of the default invariant calculation.
    #[getter]
    fn atom_invariants_generator(&self) -> Option<AtomPairAtomInvariantsGenerator> {
        self.inner
            .atom_invariants_generator
            .map(|inner| AtomPairAtomInvariantsGenerator { inner })
    }
    /// Stored 3D conformer identifier used by this operation; not its position in the conformer list.
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    /// Whether legacy stereochemistry perception is used by the fingerprint operation.
    #[getter]
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    /// Caller-supplied bond invariants in bond-index order.
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

/// Writable configuration for per-call topological-torsion atom selection and invariants.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct TopologicalTorsionCallParams {
    pub(crate) inner: ck::TopologicalTorsionCallParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionCallParams {
    /// Configure per-call topological-torsion atom selection and invariants; omitted fields use the defaults shown in the signature.
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
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    /// Caller-supplied bond invariants in bond-index order.
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
    /// Stored 3D conformer identifier used by this operation; not its position in the conformer list.
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
    /// Whether legacy stereochemistry perception is used by the fingerprint operation.
    #[getter]
    fn use_legacy_stereo_perception(&self) -> bool {
        self.inner.use_legacy_stereo_perception
    }
}
/// Reusable topological-torsion fingerprint generator with a live settings view and bit/count, dense/sparse output methods.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TopologicalTorsionFingerprintGenerator {
    pub(crate) inner: ck::TopologicalTorsionFingerprintGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionFingerprintGenerator {
    /// Construct a TopologicalTorsionFingerprintGenerator value from the supplied inputs.
    #[new]
    #[pyo3(signature=(*,params=None,atom_invariants_generator=None))]
    fn py_new(
        py: Python<'_>,
        params: Option<&TopologicalTorsionParams>,
        atom_invariants_generator: Option<&AtomPairAtomInvariantsGenerator>,
    ) -> PyResult<Self> {
        Self::new(py, params, atom_invariants_generator)
    }
    /// Construct a TopologicalTorsionFingerprintGenerator value from the supplied inputs.
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
    /// Construct a generator from serialized JSON settings; invalid JSON/options raise an error.
    #[staticmethod]
    fn from_json(py: Python<'_>, json: &str) -> PyResult<Self> {
        ck::TopologicalTorsionFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    /// Return a live settings view; edits affect subsequent calls on this generator.
    fn settings(&self) -> TopologicalTorsionSettings {
        TopologicalTorsionSettings {
            inner: self.inner.settings(),
        }
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .info_string()
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
    }
    /// Return a JSON string containing this configuration.
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
    /// Generate fixed-width bit fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<crate::fingerprint_numpy::FingerprintBatch> {
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
            .map_err(|e| crate::canonical_values::topological_torsion_pyerr(py, e))
            .and_then(|values| crate::fingerprint_numpy::FingerprintBatch::from_values(py, values))
    }
    /// Generate sparse bit fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
    /// Generate folded count fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
    /// Generate sparse count fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
/// Live topological-torsion generator settings. Changes affect subsequent generator calls.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TopologicalTorsionSettings {
    inner: ck::TopologicalTorsionSettings,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TopologicalTorsionSettings {
    /// Update the owning generator setting. Number of atoms in each topological torsion path.
    #[pyo3(name = "set_torsion_atom_count")]
    fn set_torsion_atom_count_method(&mut self, value: u32) -> PyResult<()> {
        self.set_torsion_atom_count(value)
    }
    /// Update the owning generator setting. Whether only torsion paths that are shortest between their endpoints are retained.
    #[pyo3(name = "set_only_shortest_paths")]
    fn set_only_shortest_paths_method(&mut self, value: bool) -> PyResult<()> {
        self.set_only_shortest_paths(value)
    }
    /// Update the owning generator setting. Whether stereochemical information contributes to fingerprint features.
    #[pyo3(name = "set_include_chirality")]
    fn set_include_chirality_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_chirality(value)
    }
    /// Update the owning generator setting. Whether occurrence counts are represented using threshold bits.
    #[pyo3(name = "set_count_simulation")]
    fn set_count_simulation_method(&mut self, value: bool) -> PyResult<()> {
        self.set_count_simulation(value)
    }
    /// Update the owning generator setting. Number of bins/bits in the folded fingerprint.
    #[pyo3(name = "set_fp_size")]
    fn set_fp_size_method(&mut self, value: u32) -> PyResult<()> {
        self.set_fp_size(value)
    }
    /// Update the owning generator setting. Number of hashed bit positions generated per feature.
    #[pyo3(name = "set_bits_per_feature")]
    fn set_bits_per_feature_method(&mut self, value: u32) -> PyResult<()> {
        self.set_bits_per_feature(value)
    }
    /// Update the owning generator setting. Occurrence-count thresholds used by count simulation.
    #[pyo3(name = "set_count_bounds")]
    fn set_count_bounds_method(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.set_count_bounds(value)
    }
    /// Number of atoms in each topological torsion path.
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
    /// Whether only torsion paths that are shortest between their endpoints are retained.
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
    /// Whether stereochemical information contributes to fingerprint features.
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
    /// Whether occurrence counts are represented using threshold bits.
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
    /// Number of bins/bits in the folded fingerprint.
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
    /// Number of hashed bit positions generated per feature.
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
    /// Occurrence-count thresholds used by count simulation.
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
    /// Return the parameter values associated with this result or settings view.
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
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct LegacyTopologicalTorsionParams {
    pub(crate) inner: ck::LegacyTopologicalTorsionParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl LegacyTopologicalTorsionParams {
    /// Configure topological-torsion vector entry points; omitted fields use the defaults shown in the signature.
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
    /// Configure topological-torsion vector entry points; omitted fields use the defaults shown in the signature.
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
    /// Number of atoms in each topological torsion path.
    #[getter]
    fn torsion_atom_count(&self) -> u32 {
        self.inner.torsion_atom_count
    }
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality
    }
    /// Number of bins/bits in the folded fingerprint.
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    /// Number of hashed bits generated per torsion entry.
    #[getter]
    fn bits_per_entry(&self) -> u32 {
        self.inner.bits_per_entry
    }
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
}

pyo3::create_exception!(
    cosmolkit,
    AtomCodeExplanationError,
    pyo3::exceptions::PyKeyError,
    "An atom-pair or torsion code could not be decoded into the requested explanation."
);
/// Atom-pair code calculation result containing the code and the associated molecule value.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairAtomCodeResult {
    pub(crate) inner: ck::AtomPairAtomCodeResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairAtomCodeResult {
    /// Stored residue/atom code; its interpretation is defined by the owning value type.
    #[getter]
    fn code(&self) -> u32 {
        self.inner.code
    }
    /// Return the molecule produced by this operation; the input molecule remains independently owned.
    #[getter]
    fn molecule(&self) -> crate::drawing_binding::Molecule {
        crate::drawing_binding::Molecule {
            inner: self.inner.molecule.clone(),
        }
    }
}
/// Decoded atom-pair atom-code fields: element symbol, branch count, pi-electron count and chirality.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomCodeExplanation {
    pub(crate) inner: ck::AtomCodeExplanation,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomCodeExplanation {
    /// Decode an atom-pair atom code into its element, branching, pi-electron and chirality fields.
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
    /// Chemical element symbol represented by this value.
    fn symbol(&self) -> &'static str {
        self.inner.symbol()
    }
    /// Branch-count field decoded from an atom-pair atom code.
    fn branch_count(&self) -> u32 {
        self.inner.branch_count()
    }
    /// Pi-electron field decoded from an atom-pair atom code.
    fn pi_electrons(&self) -> u32 {
        self.inner.pi_electrons()
    }
    /// Chirality field decoded from an atom-pair atom code.
    fn chirality(&self) -> Option<&'static str> {
        self.inner.chirality()
    }
}

// Registered persistent Morgan projections. All computation delegates to ck.
/// Immutable, independently captured atom-invariant provider configuration.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganAtomInvariantsGenerator {
    pub(crate) inner: ck::MorganAtomInvariantsGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganAtomInvariantsGenerator {
    /// Select atom invariants based on chemical graph connectivity.
    #[staticmethod]
    #[pyo3(signature=(include_ring_membership=true))]
    fn connectivity(include_ring_membership: bool) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::connectivity(include_ring_membership),
        }
    }
    /// Select pharmacophore/chemical-feature atom invariants for Morgan fingerprinting.
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
    /// Select atom-pair atom-code invariants with the requested chirality setting.
    #[staticmethod]
    fn atom_pair(generator: &AtomPairAtomInvariantsGenerator) -> Self {
        Self {
            inner: ck::MorganAtomInvariantsGenerator::atom_pair(generator.inner),
        }
    }
}
/// Immutable explicit bond-provider flags, captured independently of live settings.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct MorganBondInvariantsGenerator {
    pub(crate) inner: ck::MorganBondInvariantsGenerator,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganBondInvariantsGenerator {
    /// Construct a MorganBondInvariantsGenerator value from the supplied inputs.
    #[new]
    #[pyo3(signature=(*,use_bond_types=true,include_chirality=false))]
    fn py_new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self::new(use_bond_types, include_chirality)
    }
    /// Construct a MorganBondInvariantsGenerator value from the supplied inputs.
    #[staticmethod]
    #[pyo3(signature=(*,use_bond_types=true,include_chirality=false))]
    fn new(use_bond_types: bool, include_chirality: bool) -> Self {
        Self {
            inner: ck::MorganBondInvariantsGenerator::new(use_bond_types, include_chirality),
        }
    }
    /// Whether bond types contribute to fingerprint invariants.
    #[getter]
    fn use_bond_types(&self) -> bool {
        self.inner.use_bond_types()
    }
    /// Whether stereochemical information contributes to fingerprint features.
    #[getter]
    fn include_chirality(&self) -> bool {
        self.inner.include_chirality()
    }
}
/// Writable configuration for per-call Morgan atom selection and invariants.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct MorganCallParams {
    pub(crate) inner: ck::MorganCallParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganCallParams {
    /// Configure per-call Morgan atom selection and invariants; omitted fields use the defaults shown in the signature.
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
    /// Configure per-call Morgan atom selection and invariants; omitted fields use the defaults shown in the signature.
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
    /// Atom indices used as fingerprint starting centers; an omitted list uses all eligible atoms.
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    /// Atom indices excluded from fingerprint feature enumeration.
    #[getter]
    fn ignore_atoms(&self) -> Option<Vec<u32>> {
        self.inner.ignore_atoms.clone()
    }
    /// Caller-supplied atom invariants in atom-index order.
    #[getter]
    fn custom_atom_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_atom_invariants.clone()
    }
    /// Caller-supplied bond invariants in bond-index order.
    #[getter]
    fn custom_bond_invariants(&self) -> Option<Vec<u32>> {
        self.inner.custom_bond_invariants.clone()
    }
    /// Stored 3D conformer identifier used by this operation; not its position in the conformer list.
    #[getter]
    fn conformer_id(&self) -> i32 {
        self.inner.conformer_id
    }
}
/// Reusable Morgan fingerprint generator. settings() is a live view; later edits affect subsequent calls, not previously returned fingerprints.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct MorganFingerprintGenerator {
    pub(crate) inner: ck::MorganFingerprintGenerator,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganFingerprintGenerator {
    /// Construct a MorganFingerprintGenerator value from the supplied inputs.
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
    /// Construct a MorganFingerprintGenerator value from the supplied inputs.
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
    /// Construct a generator from serialized JSON settings; invalid JSON/options raise an error.
    #[staticmethod]
    fn from_json(py: Python<'_>, json: &str) -> PyResult<Self> {
        ck::MorganFingerprintGenerator::from_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Return a live settings view; edits affect subsequent calls on this generator.
    fn settings(&self) -> MorganSettings {
        MorganSettings {
            inner: self.inner.settings(),
        }
    }
    /// Return a readable description of the generator settings.
    fn info_string(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .info_string()
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
    /// Return a JSON string containing this configuration.
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
    /// Generate fixed-width bit fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
    #[pyo3(signature=(molecules,*,num_threads=1))]
    fn fingerprints(
        &self,
        py: Python<'_>,
        molecules: Vec<Option<Py<crate::drawing_binding::Molecule>>>,
        num_threads: i32,
    ) -> PyResult<crate::fingerprint_numpy::FingerprintBatch> {
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
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
            .and_then(|values| crate::fingerprint_numpy::FingerprintBatch::from_values(py, values))
    }
    /// Generate folded count fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
    /// Generate sparse bit fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
    /// Generate sparse count fingerprints for the supplied molecules in input order using this generator and the per-call parameters.
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
/// Live Morgan generator settings. Property assignment and set_* methods update the owning generator.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MorganSettings {
    inner: ck::MorganSettings,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MorganSettings {
    /// Update the owning generator setting. Morgan environment radius, measured in graph bonds.
    #[pyo3(name = "set_radius")]
    fn set_radius_method(&mut self, value: u32) -> PyResult<()> {
        self.set_radius(value)
    }
    /// Update the owning generator setting. Whether atoms with zero invariants are excluded as environment centers.
    #[pyo3(name = "set_only_nonzero_invariants")]
    fn set_only_nonzero_invariants_method(&mut self, value: bool) -> PyResult<()> {
        self.set_only_nonzero_invariants(value)
    }
    /// Update the owning generator setting. Whether redundant Morgan environments are retained.
    #[pyo3(name = "set_include_redundant_environments")]
    fn set_include_redundant_environments_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_redundant_environments(value)
    }
    /// Update the owning generator setting. Whether stereochemical information contributes to fingerprint features.
    #[pyo3(name = "set_include_chirality")]
    fn set_include_chirality_method(&mut self, value: bool) -> PyResult<()> {
        self.set_include_chirality(value)
    }
    /// Update the owning generator setting. Whether occurrence counts are represented using threshold bits.
    #[pyo3(name = "set_count_simulation")]
    fn set_count_simulation_method(&mut self, value: bool) -> PyResult<()> {
        self.set_count_simulation(value)
    }
    /// Update the owning generator setting. Number of bins/bits in the folded fingerprint.
    #[pyo3(name = "set_fp_size")]
    fn set_fp_size_method(&mut self, value: u32) -> PyResult<()> {
        self.set_fp_size(value)
    }
    /// Update the owning generator setting. Number of hashed bit positions generated per feature.
    #[pyo3(name = "set_bits_per_feature")]
    fn set_bits_per_feature_method(&mut self, value: u32) -> PyResult<()> {
        self.set_bits_per_feature(value)
    }
    /// Update the owning generator setting. Occurrence-count thresholds used by count simulation.
    #[pyo3(name = "set_count_bounds")]
    fn set_count_bounds_method(&mut self, value: Vec<u32>) -> PyResult<()> {
        self.set_count_bounds(value)
    }
    /// Morgan environment radius, measured in graph bonds.
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
    /// Whether atoms with zero invariants are excluded as environment centers.
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
    /// Whether redundant Morgan environments are retained.
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
    /// Whether stereochemical information contributes to fingerprint features.
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
    /// Whether occurrence counts are represented using threshold bits.
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
    /// Number of bins/bits in the folded fingerprint.
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
    /// Number of hashed bit positions generated per feature.
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
    /// Occurrence-count thresholds used by count simulation.
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
    /// Return the parameter values associated with this result or settings view.
    fn params(&self, py: Python<'_>) -> PyResult<MorganParams> {
        self.inner
            .params()
            .map(|inner| MorganParams { inner })
            .map_err(|e| crate::canonical_values::morgan_pyerr(py, e))
    }
}

/// Read-only constants defining atom-pair code bit widths, supported atom types and fingerprint index space.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct AtomPairsParameters;
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl AtomPairsParameters {
    /// Return the version identifier recorded by this value or implementation.
    #[staticmethod]
    fn version() -> &'static str {
        ck::AtomPairsParameters::version()
    }
    /// Bit width used for the atom-type field of atom-pair codes.
    #[staticmethod]
    fn num_type_bits() -> u32 {
        ck::AtomPairsParameters::num_type_bits()
    }
    /// Bit width used for the pi-electron field of atom-pair codes.
    #[staticmethod]
    fn num_pi_bits() -> u32 {
        ck::AtomPairsParameters::num_pi_bits()
    }
    /// Bit width used for the branch-count field of atom-pair codes.
    #[staticmethod]
    fn num_branch_bits() -> u32 {
        ck::AtomPairsParameters::num_branch_bits()
    }
    /// Bit width used for the chirality field of atom-pair codes.
    #[staticmethod]
    fn num_chiral_bits() -> u32 {
        ck::AtomPairsParameters::num_chiral_bits()
    }
    /// Total bit width of the packed atom-pair atom code.
    #[staticmethod]
    fn code_size() -> u32 {
        ck::AtomPairsParameters::code_size()
    }
    /// Bit width used for atom-pair path lengths.
    #[staticmethod]
    fn num_path_bits() -> u32 {
        ck::AtomPairsParameters::num_path_bits()
    }
    /// Maximum path length representable by the atom-pair code format.
    #[staticmethod]
    fn max_path_length() -> u32 {
        ck::AtomPairsParameters::max_path_length()
    }
    /// Logical atom-pair fingerprint index-space size.
    #[staticmethod]
    fn num_atom_pair_fingerprint_bits() -> u32 {
        ck::AtomPairsParameters::num_atom_pair_fingerprint_bits()
    }
    /// Element atomic numbers included in the atom-pair atom-type table.
    #[staticmethod]
    fn atom_types() -> Vec<u32> {
        ck::AtomPairsParameters::atom_types()
    }
}
