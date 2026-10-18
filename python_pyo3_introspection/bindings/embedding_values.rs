use crate::chemical_values::enum_value;
use crate::native::Molecule;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

enum_value!(CoordinateDimension, { TwoD, ThreeD });

#[pyclass(module = "fixed_projection._native", eq, eq_int, from_py_object)]
#[derive(Clone, Copy, PartialEq)]
pub enum ConformerErrorKind {
    Unsupported,
    CannotNormalizeZeroLengthVector,
    GenerationFailed,
    InvalidEmbedParametersJson,
}

#[pyclass(module="fixed_projection._native", extends=pyo3::exceptions::PyValueError)]
pub struct ConformerError {
    inner: ck::ConformerError,
}
#[pymethods]
impl ConformerError {
    #[getter]
    fn variant(&self) -> ConformerErrorKind {
        match &self.inner {
            ck::ConformerError::Unsupported => ConformerErrorKind::Unsupported,
            ck::ConformerError::CannotNormalizeZeroLengthVector => {
                ConformerErrorKind::CannotNormalizeZeroLengthVector
            }
            ck::ConformerError::GenerationFailed(_) => ConformerErrorKind::GenerationFailed,
            ck::ConformerError::InvalidEmbedParametersJson(_) => {
                ConformerErrorKind::InvalidEmbedParametersJson
            }
        }
    }
    #[getter]
    fn domain(&self) -> &'static str {
        "conformer"
    }
    // These unnamed tuple payloads contain source messages, not nested errors.
    #[getter]
    fn message(&self) -> Option<&str> {
        match &self.inner {
            ck::ConformerError::GenerationFailed(message)
            | ck::ConformerError::InvalidEmbedParametersJson(message) => Some(message),
            ck::ConformerError::Unsupported
            | ck::ConformerError::CannotNormalizeZeroLengthVector => None,
        }
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl ConformerError {
    pub(crate) fn from_rust(inner: ck::ConformerError) -> PyErr {
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(error) => PyErr::from_value(error.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}

#[pyclass(module = "fixed_projection._native", from_py_object)]
#[derive(Clone)]
pub struct EmbedParams {
    pub(crate) inner: ck::EmbedParams,
}
macro_rules! embedding_parameters {
    (fields {$($field:ident: $ty:ty),+ $(,)?}, presets {$($preset:ident),+ $(,)?}) => {
        paste::paste! {
            #[pymethods]
            impl EmbedParams {
                #[new]
                #[pyo3(signature=(*, $($field=ck::EmbedParams::default().$field),+))]
                fn construct($($field: $ty),+) -> Self {
                    // The Rust owner has private bounds storage and nonprojected
                    // callback/diagnostic state. Its real constructor initializes
                    // those members; the canonical configuration exposes only
                    // the following declared fields, with no replacement policy.
                    let mut inner = ck::EmbedParams::new();
                    $(inner.$field = $field;)+
                    Self { inner }
                }
                #[staticmethod] fn new() -> Self { Self { inner: ck::EmbedParams::new() } }
                $(#[staticmethod] fn $preset() -> Self { Self { inner: ck::EmbedParams::$preset() } })+
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
                fn to_json(&self) -> String { self.inner.to_json() }
                fn with_json(&self, json: &str) -> PyResult<Self> {
                    self.inner.with_json(json).map(|inner| Self { inner }).map_err(ConformerError::from_rust)
                }
                #[getter]
                fn failures(&self) -> Vec<u32> { self.inner.failures().clone() }
            }
        }
    };
}
embedding_parameters! {
    fields {
        max_iterations: u32, num_threads: i32, random_seed: i32, clear_conformers: bool,
        use_random_coords: bool, box_size_mult: f64, rand_neg_eig: bool, num_zero_fail: u32,
        coord_map: Option<BTreeMap<i32,[f64;3]>>, optimizer_force_tol: f64,
        ignore_smoothing_failures: bool, enforce_chirality: bool, use_exp_torsion_angle_prefs: bool,
        use_basic_knowledge: bool, verbose: bool, basin_thresh: f64, prune_rms_thresh: f64,
        only_heavy_atoms_for_rms: bool, et_version: u32, embed_fragments_separately: bool,
        use_small_ring_torsions: bool, use_macrocycle_torsions: bool, use_macrocycle14config: bool,
        timeout: u32, cpci: Option<BTreeMap<(u32,u32),f64>>, force_trans_amides: bool,
        use_symmetry_for_pruning: bool, bounds_mat_force_scaling: f64, track_failures: bool,
        enable_sequential_random_seeds: bool, symmetrize_conjugated_terminal_groups_for_pruning: bool
    },
    presets { dg, kdg, etdg, etdg_v2, etkdg, etkdg_v2, etkdg_v3, sr_etkdg_v3 }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct EmbedMoleculeResult {
    pub(crate) inner: ck::EmbedMoleculeResult,
}
#[pymethods]
impl EmbedMoleculeResult {
    #[new]
    fn new(molecule: &Molecule, conf_id: i32, params: &EmbedParams) -> Self {
        Self {
            inner: ck::EmbedMoleculeResult {
                molecule: molecule.inner.clone(),
                conf_id,
                params: params.inner.clone(),
            },
        }
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conf_id(&self) -> i32 {
        self.inner.conf_id()
    }
    fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
    fn ok(&self) -> bool {
        self.inner.ok()
    }
}
#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct EmbedMultipleConfsResult {
    pub(crate) inner: ck::EmbedMultipleConfsResult,
}
#[pymethods]
impl EmbedMultipleConfsResult {
    #[new]
    fn new(
        molecule: &Molecule,
        conf_ids: Vec<i32>,
        requested_num_confs: u32,
        params: &EmbedParams,
    ) -> Self {
        Self {
            inner: ck::EmbedMultipleConfsResult {
                molecule: molecule.inner.clone(),
                conf_ids,
                requested_num_confs,
                params: params.inner.clone(),
            },
        }
    }
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    fn conf_ids(&self) -> &[i32] {
        self.inner.conf_ids()
    }
    fn requested_num_confs(&self) -> u32 {
        self.inner.requested_num_confs()
    }
    fn generated_count(&self) -> usize {
        self.inner.generated_count()
    }
    fn params(&self) -> EmbedParams {
        EmbedParams {
            inner: self.inner.params().clone(),
        }
    }
}

#[cfg(test)]
mod diagnostic_tests {
    use super::*;

    #[test]
    fn actual_embedding_diagnostics_survive_result_parameter_projection() {
        let molecule = ck::Molecule::from_smiles("CCO").unwrap();
        let mut options = ck::EmbedParams::etkdg_v3();
        options.random_seed = 17;
        options.track_failures = true;
        let result = molecule
            .with_3d_conformer_result_with_params(&options)
            .unwrap();
        let expected = result.params().failures().clone();
        assert!(!expected.is_empty());
        assert!(options.failures().is_empty());
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let result = Py::new(py, EmbedMoleculeResult { inner: result })?;
            let options = result.bind(py).call_method0("params")?;
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            assert!(options.getattr("track_failures")?.extract::<bool>()?);
            assert!(
                options
                    .setattr("failures", vec![1_u32])
                    .unwrap_err()
                    .is_instance_of::<pyo3::exceptions::PyAttributeError>(py)
            );
            let copy = options.getattr("failures")?;
            copy.call_method1("append", (u32::MAX,))?;
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            drop(result);
            assert_eq!(
                options.getattr("failures")?.extract::<Vec<u32>>()?,
                expected
            );
            Ok(())
        })
        .unwrap();
    }
}
