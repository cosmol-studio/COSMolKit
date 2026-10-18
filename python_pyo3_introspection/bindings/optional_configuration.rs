use crate::construction_values::ReactionCoordinateSelection;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyclass(module = "fixed_projection._native")]
pub struct ReactionRunParams {
    pub(crate) inner: ck::ReactionRunParams,
}
#[pymethods]
impl ReactionRunParams {
    #[new]
    #[pyo3(signature=(*, max_products=ck::ReactionRunParams::default().max_products, coordinate_selections=None, copy_atom_properties=ck::ReactionRunParams::default().copy_atom_properties))]
    fn new(
        max_products: u32,
        coordinate_selections: Option<Vec<ReactionCoordinateSelection>>,
        copy_atom_properties: bool,
    ) -> Self {
        let selections = coordinate_selections
            .map(|values| values.into_iter().map(|value| value.inner).collect())
            .unwrap_or_else(|| ck::ReactionRunParams::default().coordinate_selections);
        Self {
            inner: ck::ReactionRunParams::new(max_products, selections, copy_atom_properties),
        }
    }
    #[getter]
    fn max_products(&self) -> u32 {
        self.inner.max_products
    }
    #[setter]
    fn set_max_products(&mut self, value: u32) {
        self.inner.max_products = value;
    }
    #[getter]
    fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> {
        self.inner
            .coordinate_selections
            .iter()
            .copied()
            .map(|inner| ReactionCoordinateSelection { inner })
            .collect()
    }
    #[setter]
    fn set_coordinate_selections(&mut self, value: Option<Vec<ReactionCoordinateSelection>>) {
        self.inner.coordinate_selections = value
            .map(|values| values.into_iter().map(|value| value.inner).collect())
            .unwrap_or_else(|| ck::ReactionRunParams::default().coordinate_selections);
    }
    #[getter]
    fn copy_atom_properties(&self) -> bool {
        self.inner.copy_atom_properties
    }
    #[setter]
    fn set_copy_atom_properties(&mut self, value: bool) {
        self.inner.copy_atom_properties = value;
    }
}

#[pyclass(module = "fixed_projection._native")]
pub struct Coordinate2DParams {
    pub(crate) inner: ck::Coordinate2DParams,
}
#[pymethods]
impl Coordinate2DParams {
    #[new]
    #[pyo3(signature=(*, coordinate_map=None, canonical_orientation=ck::Coordinate2DParams::default().canonical_orientation, clear_existing_2d=ck::Coordinate2DParams::default().clear_existing_2d, flips_per_sample=ck::Coordinate2DParams::default().flips_per_sample, samples=ck::Coordinate2DParams::default().samples, sample_seed=ck::Coordinate2DParams::default().sample_seed, permute_degree_four=ck::Coordinate2DParams::default().permute_degree_four, force_rdkit=ck::Coordinate2DParams::default().force_rdkit, use_ring_templates=ck::Coordinate2DParams::default().use_ring_templates, use_path_angle_expansion=ck::Coordinate2DParams::default().use_path_angle_expansion))]
    fn new(
        coordinate_map: Option<BTreeMap<usize, [f64; 2]>>,
        canonical_orientation: bool,
        clear_existing_2d: bool,
        flips_per_sample: u32,
        samples: u32,
        sample_seed: i32,
        permute_degree_four: bool,
        force_rdkit: bool,
        use_ring_templates: bool,
        use_path_angle_expansion: bool,
    ) -> Self {
        Self {
            inner: ck::Coordinate2DParams {
                coordinate_map: coordinate_map
                    .unwrap_or_else(|| ck::Coordinate2DParams::default().coordinate_map),
                canonical_orientation,
                clear_existing_2d,
                flips_per_sample,
                samples,
                sample_seed,
                permute_degree_four,
                force_rdkit,
                use_ring_templates,
                use_path_angle_expansion,
            },
        }
    }
    #[getter]
    fn coordinate_map(&self) -> BTreeMap<usize, [f64; 2]> {
        self.inner.coordinate_map.clone()
    }
    #[setter]
    fn set_coordinate_map(&mut self, value: Option<BTreeMap<usize, [f64; 2]>>) {
        self.inner.coordinate_map =
            value.unwrap_or_else(|| ck::Coordinate2DParams::default().coordinate_map);
    }
    #[getter]
    fn canonical_orientation(&self) -> bool {
        self.inner.canonical_orientation
    }
    #[setter]
    fn set_canonical_orientation(&mut self, value: bool) {
        self.inner.canonical_orientation = value;
    }
    #[getter]
    fn clear_existing_2d(&self) -> bool {
        self.inner.clear_existing_2d
    }
    #[setter]
    fn set_clear_existing_2d(&mut self, value: bool) {
        self.inner.clear_existing_2d = value;
    }
    #[getter]
    fn flips_per_sample(&self) -> u32 {
        self.inner.flips_per_sample
    }
    #[setter]
    fn set_flips_per_sample(&mut self, value: u32) {
        self.inner.flips_per_sample = value;
    }
    #[getter]
    fn samples(&self) -> u32 {
        self.inner.samples
    }
    #[setter]
    fn set_samples(&mut self, value: u32) {
        self.inner.samples = value;
    }
    #[getter]
    fn sample_seed(&self) -> i32 {
        self.inner.sample_seed
    }
    #[setter]
    fn set_sample_seed(&mut self, value: i32) {
        self.inner.sample_seed = value;
    }
    #[getter]
    fn permute_degree_four(&self) -> bool {
        self.inner.permute_degree_four
    }
    #[setter]
    fn set_permute_degree_four(&mut self, value: bool) {
        self.inner.permute_degree_four = value;
    }
    #[getter]
    fn force_rdkit(&self) -> bool {
        self.inner.force_rdkit
    }
    #[setter]
    fn set_force_rdkit(&mut self, value: bool) {
        self.inner.force_rdkit = value;
    }
    #[getter]
    fn use_ring_templates(&self) -> bool {
        self.inner.use_ring_templates
    }
    #[setter]
    fn set_use_ring_templates(&mut self, value: bool) {
        self.inner.use_ring_templates = value;
    }
    #[getter]
    fn use_path_angle_expansion(&self) -> bool {
        self.inner.use_path_angle_expansion
    }
    #[setter]
    fn set_use_path_angle_expansion(&mut self, value: bool) {
        self.inner.use_path_angle_expansion = value;
    }
}
