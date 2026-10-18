use crate::native::Fingerprint;
use cosmolkit as ck;
use pyo3::prelude::*;

#[pyclass(module = "fixed_projection._native")]
pub struct LayeredFingerprintParams {
    pub(crate) inner: ck::LayeredFingerprintParams,
}
#[pymethods]
impl LayeredFingerprintParams {
    #[new]
    #[pyo3(signature = (*, layers=ck::LayeredFingerprintParams::default().layers.bits(), min_path=ck::LayeredFingerprintParams::default().min_path, max_path=ck::LayeredFingerprintParams::default().max_path, fp_size=ck::LayeredFingerprintParams::default().fp_size, atom_counts=None, set_only_bits=None, branched_paths=ck::LayeredFingerprintParams::default().branched_paths, from_atoms=None))]
    fn new(
        layers: u32,
        min_path: u32,
        max_path: u32,
        fp_size: u32,
        atom_counts: Option<Vec<u32>>,
        set_only_bits: Option<&Fingerprint>,
        branched_paths: bool,
        from_atoms: Option<Vec<u32>>,
    ) -> Self {
        Self {
            inner: ck::LayeredFingerprintParams {
                layers: ck::LayeredFingerprintLayers::from_bits_retain(layers),
                min_path,
                max_path,
                fp_size,
                atom_counts,
                set_only_bits: set_only_bits.map(|value| value.inner.clone()),
                branched_paths,
                from_atoms,
            },
        }
    }
    #[getter]
    fn layers(&self) -> u32 {
        self.inner.layers.bits()
    }
    #[setter]
    fn set_layers(&mut self, value: u32) {
        self.inner.layers = ck::LayeredFingerprintLayers::from_bits_retain(value);
    }
    #[getter]
    fn min_path(&self) -> u32 {
        self.inner.min_path
    }
    #[setter]
    fn set_min_path(&mut self, value: u32) {
        self.inner.min_path = value;
    }
    #[getter]
    fn max_path(&self) -> u32 {
        self.inner.max_path
    }
    #[setter]
    fn set_max_path(&mut self, value: u32) {
        self.inner.max_path = value;
    }
    #[getter]
    fn fp_size(&self) -> u32 {
        self.inner.fp_size
    }
    #[setter]
    fn set_fp_size(&mut self, value: u32) {
        self.inner.fp_size = value;
    }
    #[getter]
    fn atom_counts(&self) -> Option<Vec<u32>> {
        self.inner.atom_counts.clone()
    }
    #[setter]
    fn set_atom_counts(&mut self, value: Option<Vec<u32>>) {
        self.inner.atom_counts = value;
    }
    #[getter]
    fn set_only_bits(&self) -> Option<Fingerprint> {
        self.inner
            .set_only_bits
            .clone()
            .map(|inner| Fingerprint { inner })
    }
    #[setter]
    fn set_set_only_bits(&mut self, value: Option<&Fingerprint>) {
        self.inner.set_only_bits = value.map(|value| value.inner.clone());
    }
    #[getter]
    fn branched_paths(&self) -> bool {
        self.inner.branched_paths
    }
    #[setter]
    fn set_branched_paths(&mut self, value: bool) {
        self.inner.branched_paths = value;
    }
    #[getter]
    fn from_atoms(&self) -> Option<Vec<u32>> {
        self.inner.from_atoms.clone()
    }
    #[setter]
    fn set_from_atoms(&mut self, value: Option<Vec<u32>>) {
        self.inner.from_atoms = value;
    }
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
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

#[pyclass(module = "fixed_projection._native")]
pub struct AlignmentParameters {
    pub(crate) inner: ck::AlignmentParameters,
}
#[pymethods]
impl AlignmentParameters {
    #[new]
    #[pyo3(signature = (*, probe_conformer_id=ck::AlignmentParameters::default().probe_conformer_id, reference_conformer_id=ck::AlignmentParameters::default().reference_conformer_id, atom_map=None, weights=None, reflect=ck::AlignmentParameters::default().reflect, max_iterations=ck::AlignmentParameters::default().max_iterations))]
    fn construct(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_map: Option<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            inner: ck::AlignmentParameters::new(
                probe_conformer_id,
                reference_conformer_id,
                atom_map.map(|maps| maps.into_iter().map(|value| value.inner).collect()),
                weights,
                reflect,
                max_iterations,
            ),
        }
    }
    #[staticmethod]
    #[pyo3(signature = (probe_conformer_id=ck::AlignmentParameters::default().probe_conformer_id, reference_conformer_id=ck::AlignmentParameters::default().reference_conformer_id, atom_map=None, weights=None, reflect=ck::AlignmentParameters::default().reflect, max_iterations=ck::AlignmentParameters::default().max_iterations))]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_map: Option<Vec<AlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self::construct(
            probe_conformer_id,
            reference_conformer_id,
            atom_map,
            weights,
            reflect,
            max_iterations,
        )
    }
    #[getter]
    fn probe_conformer_id(&self) -> i32 {
        self.inner.probe_conformer_id
    }
    #[setter]
    fn set_probe_conformer_id(&mut self, value: i32) {
        self.inner.probe_conformer_id = value;
    }
    #[getter]
    fn reference_conformer_id(&self) -> i32 {
        self.inner.reference_conformer_id
    }
    #[setter]
    fn set_reference_conformer_id(&mut self, value: i32) {
        self.inner.reference_conformer_id = value;
    }
    #[getter]
    fn atom_map(&self) -> Option<Vec<AlignmentAtomMap>> {
        self.inner.atom_map.as_ref().map(|maps| {
            maps.iter()
                .copied()
                .map(|inner| AlignmentAtomMap { inner })
                .collect()
        })
    }
    #[setter]
    fn set_atom_map(&mut self, value: Option<Vec<AlignmentAtomMap>>) {
        self.inner.atom_map = value.map(|maps| maps.into_iter().map(|value| value.inner).collect());
    }
    #[getter]
    fn weights(&self) -> Option<Vec<f64>> {
        self.inner.weights.clone()
    }
    #[setter]
    fn set_weights(&mut self, value: Option<Vec<f64>>) {
        self.inner.weights = value;
    }
    #[getter]
    fn reflect(&self) -> bool {
        self.inner.reflect
    }
    #[setter]
    fn set_reflect(&mut self, value: bool) {
        self.inner.reflect = value;
    }
    #[getter]
    fn max_iterations(&self) -> u32 {
        self.inner.max_iterations
    }
    #[setter]
    fn set_max_iterations(&mut self, value: u32) {
        self.inner.max_iterations = value;
    }
}
