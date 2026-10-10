//! Thin native projections of the canonical conformer facade.
use crate::drawing_binding::Molecule;
use ::cosmolkit as ck;

use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use std::collections::BTreeMap;
/// Writable configuration for distance-geometry 3D conformer embedding.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct EmbedParams {
    pub(crate) inner: ck::EmbedParams,
}
// Preset-only fields and supplied bounds are not constructor keywords: retain
// them when changing one public setting rather than rebuilding the preset.
#[cosmolkit_macros::python_configuration(fieldwise)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl EmbedParams {
    /// Configure distance-geometry 3D conformer embedding; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature = (*, max_iterations=0, num_threads=1, random_seed=-1, clear_conformers=true, use_random_coords=false, box_size_mult=2.0, rand_neg_eig=true, num_zero_fail=1, coord_map=None, optimizer_force_tol=1e-3, ignore_smoothing_failures=false, enforce_chirality=true, use_exp_torsion_angle_prefs=false, use_basic_knowledge=false, verbose=false, basin_thresh=5.0, prune_rms_thresh=-1.0, only_heavy_atoms_for_rms=true, et_version=2, embed_fragments_separately=true, use_small_ring_torsions=false, use_macrocycle_torsions=false, use_macrocycle14config=false, timeout=0, cpci=None, force_trans_amides=true, use_symmetry_for_pruning=true, bounds_mat_force_scaling=1.0, track_failures=false, enable_sequential_random_seeds=false, symmetrize_conjugated_terminal_groups_for_pruning=true))]
    #[pyo3(
        text_signature = "(*, max_iterations=0, num_threads=1, random_seed=-1, clear_conformers=True, use_random_coords=False, box_size_mult=2.0, rand_neg_eig=True, num_zero_fail=1, coord_map=None, optimizer_force_tol=1e-3, ignore_smoothing_failures=False, enforce_chirality=True, use_exp_torsion_angle_prefs=False, use_basic_knowledge=False, verbose=False, basin_thresh=5.0, prune_rms_thresh=-1.0, only_heavy_atoms_for_rms=True, et_version=2, embed_fragments_separately=True, use_small_ring_torsions=False, use_macrocycle_torsions=False, use_macrocycle14config=False, timeout=0, cpci=None, force_trans_amides=True, use_symmetry_for_pruning=True, bounds_mat_force_scaling=1.0, track_failures=False, enable_sequential_random_seeds=False, symmetrize_conjugated_terminal_groups_for_pruning=True)"
    )]
    fn py_new(
        max_iterations: u32,
        num_threads: i32,
        random_seed: i32,
        clear_conformers: bool,
        use_random_coords: bool,
        box_size_mult: f64,
        rand_neg_eig: bool,
        num_zero_fail: u32,
        coord_map: Option<BTreeMap<i32, [f64; 3]>>,
        optimizer_force_tol: f64,
        ignore_smoothing_failures: bool,
        enforce_chirality: bool,
        use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool,
        verbose: bool,
        basin_thresh: f64,
        prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool,
        et_version: u32,
        embed_fragments_separately: bool,
        use_small_ring_torsions: bool,
        use_macrocycle_torsions: bool,
        use_macrocycle14config: bool,
        timeout: u32,
        cpci: Option<BTreeMap<(u32, u32), f64>>,
        force_trans_amides: bool,
        use_symmetry_for_pruning: bool,
        bounds_mat_force_scaling: f64,
        track_failures: bool,
        enable_sequential_random_seeds: bool,
        symmetrize_conjugated_terminal_groups_for_pruning: bool,
    ) -> Self {
        // Only a new detached configuration is initialized here. Generation
        // borrows this frozen value; source execution progress lives in reports.
        let mut inner = ck::EmbedParams::new();
        inner.max_iterations = max_iterations;
        inner.num_threads = num_threads;
        inner.random_seed = random_seed;
        inner.clear_conformers = clear_conformers;
        inner.use_random_coords = use_random_coords;
        inner.box_size_mult = box_size_mult;
        inner.rand_neg_eig = rand_neg_eig;
        inner.num_zero_fail = num_zero_fail;
        inner.coord_map = coord_map;
        inner.optimizer_force_tol = optimizer_force_tol;
        inner.ignore_smoothing_failures = ignore_smoothing_failures;
        inner.enforce_chirality = enforce_chirality;
        inner.use_exp_torsion_angle_prefs = use_exp_torsion_angle_prefs;
        inner.use_basic_knowledge = use_basic_knowledge;
        inner.verbose = verbose;
        inner.basin_thresh = basin_thresh;
        inner.prune_rms_thresh = prune_rms_thresh;
        inner.only_heavy_atoms_for_rms = only_heavy_atoms_for_rms;
        inner.et_version = et_version;
        inner.embed_fragments_separately = embed_fragments_separately;
        inner.use_small_ring_torsions = use_small_ring_torsions;
        inner.use_macrocycle_torsions = use_macrocycle_torsions;
        inner.use_macrocycle14config = use_macrocycle14config;
        inner.timeout = timeout;
        inner.cpci = cpci;
        inner.force_trans_amides = force_trans_amides;
        inner.use_symmetry_for_pruning = use_symmetry_for_pruning;
        inner.bounds_mat_force_scaling = bounds_mat_force_scaling;
        inner.track_failures = track_failures;
        inner.enable_sequential_random_seeds = enable_sequential_random_seeds;
        inner.symmetrize_conjugated_terminal_groups_for_pruning =
            symmetrize_conjugated_terminal_groups_for_pruning;
        Self { inner }
    }
    /// Configure distance-geometry 3D conformer embedding; omitted fields use the defaults shown in the signature.
    #[staticmethod]
    fn new() -> Self {
        Self {
            inner: ck::EmbedParams::new(),
        }
    }
    /// Return distance-geometry embedding parameters without experimental torsion or basic-knowledge constraints.
    #[staticmethod]
    fn dg() -> Self {
        Self {
            inner: ck::EmbedParams::dg(),
        }
    }
    /// Return distance-geometry parameters with basic chemical knowledge constraints.
    #[staticmethod]
    fn kdg() -> Self {
        Self {
            inner: ck::EmbedParams::kdg(),
        }
    }
    /// Return distance-geometry parameters with experimental torsion preferences.
    #[staticmethod]
    fn etdg() -> Self {
        Self {
            inner: ck::EmbedParams::etdg(),
        }
    }
    /// Return distance-geometry parameters using version 2 experimental torsion preferences.
    #[staticmethod]
    fn etdg_v2() -> Self {
        Self {
            inner: ck::EmbedParams::etdg_v2(),
        }
    }
    /// Return embedding parameters combining experimental torsions and basic chemical knowledge.
    #[staticmethod]
    fn etkdg() -> Self {
        Self {
            inner: ck::EmbedParams::etkdg(),
        }
    }
    /// Return ETKDG embedding parameters using version 2 torsion preferences.
    #[staticmethod]
    fn etkdg_v2() -> Self {
        Self {
            inner: ck::EmbedParams::etkdg_v2(),
        }
    }
    /// Return ETKDGv3 embedding parameters, including its macrocycle geometry configuration.
    #[staticmethod]
    fn etkdg_v3() -> Self {
        Self {
            inner: ck::EmbedParams::etkdg_v3(),
        }
    }
    /// Return ETKDGv3 parameters with additional small-ring torsion preferences.
    #[staticmethod]
    fn sr_etkdg_v3() -> Self {
        Self {
            inner: ck::EmbedParams::sr_etkdg_v3(),
        }
    }
    /// Maximum iterations allowed by the optimizer or embedding algorithm.
    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations.clone()
    }
    /// Requested worker count; interpretation of zero follows the corresponding operation.
    #[getter]
    fn num_threads(&self) -> i32 {
        self.inner.num_threads.clone()
    }
    /// Embedding random seed; a negative value requests nondeterministic seeding.
    #[getter]
    fn random_seed(&self) -> i32 {
        self.inner.random_seed.clone()
    }
    /// Whether existing 3D conformers are cleared before embedding.
    #[getter]
    fn clear_conformers(&self) -> bool {
        self.inner.clear_conformers.clone()
    }
    /// Whether embedding starts from random Cartesian coordinates instead of distance-matrix eigenvectors.
    #[getter]
    fn use_random_coords(&self) -> bool {
        self.inner.use_random_coords.clone()
    }
    /// Scale controlling the initial random-coordinate box.
    #[getter]
    fn box_size_mult(&self) -> f64 {
        self.inner.box_size_mult.clone()
    }
    /// Whether negative eigenvalues are handled by randomizing the corresponding coordinates.
    #[getter]
    fn rand_neg_eig(&self) -> bool {
        self.inner.rand_neg_eig.clone()
    }
    /// Number of zero eigenvalues that causes embedding to fail.
    #[getter]
    fn num_zero_fail(&self) -> u32 {
        self.inner.num_zero_fail.clone()
    }
    /// Fixed 3D coordinates indexed by atom index for embedding.
    #[getter]
    fn coord_map(&self) -> Option<BTreeMap<i32, [f64; 3]>> {
        self.inner.coord_map.clone()
    }
    /// Force tolerance for the distance-geometry optimizer.
    #[getter]
    fn optimizer_force_tol(&self) -> f64 {
        self.inner.optimizer_force_tol.clone()
    }
    /// Whether embedding may proceed after bounds-smoothing failure.
    #[getter]
    fn ignore_smoothing_failures(&self) -> bool {
        self.inner.ignore_smoothing_failures.clone()
    }
    /// Whether embedding enforces specified stereochemistry.
    #[getter]
    fn enforce_chirality(&self) -> bool {
        self.inner.enforce_chirality.clone()
    }
    /// Whether experimental torsion-angle preferences are used.
    #[getter]
    fn use_exp_torsion_angle_prefs(&self) -> bool {
        self.inner.use_exp_torsion_angle_prefs.clone()
    }
    /// Whether basic chemical geometry constraints are used during embedding.
    #[getter]
    fn use_basic_knowledge(&self) -> bool {
        self.inner.use_basic_knowledge.clone()
    }
    /// Whether diagnostic output is enabled.
    #[getter]
    fn verbose(&self) -> bool {
        self.inner.verbose.clone()
    }
    /// Energy threshold used when considering distinct embedding basins.
    #[getter]
    fn basin_thresh(&self) -> f64 {
        self.inner.basin_thresh.clone()
    }
    /// RMSD threshold in angstroms for pruning similar conformers; a negative value disables pruning.
    #[getter]
    fn prune_rms_thresh(&self) -> f64 {
        self.inner.prune_rms_thresh.clone()
    }
    /// Whether conformer RMSD pruning excludes hydrogen atoms.
    #[getter]
    fn only_heavy_atoms_for_rms(&self) -> bool {
        self.inner.only_heavy_atoms_for_rms.clone()
    }
    /// Version of the experimental torsion-angle preference parameter set.
    #[getter]
    fn et_version(&self) -> u32 {
        self.inner.et_version.clone()
    }
    /// Whether disconnected fragments are embedded independently.
    #[getter]
    fn embed_fragments_separately(&self) -> bool {
        self.inner.embed_fragments_separately.clone()
    }
    /// Whether small-ring torsion preferences are included.
    #[getter]
    fn use_small_ring_torsions(&self) -> bool {
        self.inner.use_small_ring_torsions.clone()
    }
    /// Whether macrocycle torsion preferences are included.
    #[getter]
    fn use_macrocycle_torsions(&self) -> bool {
        self.inner.use_macrocycle_torsions.clone()
    }
    /// Whether macrocycle-specific 1-4 distance bounds are used.
    #[getter]
    fn use_macrocycle14config(&self) -> bool {
        self.inner.use_macrocycle14config.clone()
    }
    /// Time limit in seconds for the corresponding search or embedding operation.
    #[getter]
    fn timeout(&self) -> u32 {
        self.inner.timeout.clone()
    }
    /// Custom pairwise Coulomb-like interaction strengths keyed by atom-index pairs.
    #[getter]
    fn cpci(&self) -> Option<BTreeMap<(u32, u32), f64>> {
        self.inner.cpci.clone()
    }
    /// Whether amide bonds are constrained toward trans geometry during embedding.
    #[getter]
    fn force_trans_amides(&self) -> bool {
        self.inner.force_trans_amides.clone()
    }
    /// Whether symmetry-equivalent atom mappings are considered during RMSD pruning.
    #[getter]
    fn use_symmetry_for_pruning(&self) -> bool {
        self.inner.use_symmetry_for_pruning.clone()
    }
    /// Scale applied to distance-bound constraint forces.
    #[getter]
    fn bounds_mat_force_scaling(&self) -> f64 {
        self.inner.bounds_mat_force_scaling.clone()
    }
    /// Whether counts of embedding failure causes are recorded.
    #[getter]
    fn track_failures(&self) -> bool {
        self.inner.track_failures.clone()
    }
    /// Counts of recorded embedding failures by failure-cause category.
    #[getter]
    fn failures(&self) -> Vec<u32> {
        self.inner.failures.clone()
    }
    /// Whether successive embeddings derive seeds sequentially from the configured seed.
    #[getter]
    fn enable_sequential_random_seeds(&self) -> bool {
        self.inner.enable_sequential_random_seeds.clone()
    }
    /// Whether conjugated terminal groups are treated as symmetry-equivalent when pruning conformers.
    #[getter]
    fn symmetrize_conjugated_terminal_groups_for_pruning(&self) -> bool {
        self.inner
            .symmetrize_conjugated_terminal_groups_for_pruning
            .clone()
    }
    /// Return a JSON string containing this configuration.
    fn to_json(&self) -> String {
        self.inner.to_json()
    }
    /// Return an independent configuration updated from the supplied JSON text; invalid JSON/options raise an error.
    fn with_json(&self, json: &str) -> PyResult<Self> {
        self.inner
            .with_json(json)
            .map(|inner| Self { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
    fn __repr__(&self) -> String {
        format!(
            "EmbedParams(random_seed={}, num_threads={}, prune_rms_thresh={}, clear_conformers={})",
            self.inner.random_seed,
            self.inner.num_threads,
            self.inner.prune_rms_thresh,
            self.inner.clear_conformers
        )
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<EmbedParams>()?;
    module.add_class::<EmbedMoleculeResult>()?;
    module.add_class::<EmbedMultipleConfsResult>()?;
    Ok(())
}

/// Single-conformer embedding result: molecule, updated parameters, conformer ID and success flag.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct EmbedMoleculeResult {
    pub(crate) inner: ck::EmbedMoleculeResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl EmbedMoleculeResult {
    /// Return the molecule produced by this operation; the input molecule remains independently owned.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    /// Return the parameter values associated with this result or settings view.
    fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
    /// Identifier of the generated conformer, or the operation failure sentinel.
    fn conf_id(&self) -> i32 {
        self.inner.conf_id()
    }
    /// Whether the embedding attempt produced a conformer.
    fn ok(&self) -> bool {
        self.inner.ok()
    }
    fn __repr__(&self) -> String {
        format!(
            "EmbedMoleculeResult(conf_id={}, ok={})",
            self.inner.conf_id(),
            self.inner.ok()
        )
    }
}

/// Multi-conformer embedding result containing generated IDs, requested/generated counts and updated parameters.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct EmbedMultipleConfsResult {
    pub(crate) inner: ck::EmbedMultipleConfsResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl EmbedMultipleConfsResult {
    /// Return the molecule produced by this operation; the input molecule remains independently owned.
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    /// Return the parameter values associated with this result or settings view.
    fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
    /// Identifiers of the successfully generated conformers in generation order.
    fn conf_ids(&self) -> Vec<i32> {
        self.inner.conf_ids().to_vec()
    }
    /// Number of successfully generated conformers.
    fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
    /// Number of conformers requested from the embedding operation.
    fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    fn __repr__(&self) -> String {
        format!(
            "EmbedMultipleConfsResult(requested_num_confs={}, generated_count={})",
            self.inner.requested_num_confs(),
            self.inner.generated_count()
        )
    }
}
