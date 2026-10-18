use crate::collection_params::AlignmentAtomMap;
use cosmolkit as ck;
use pyo3::prelude::*;

macro_rules! map_parameters {
    ($name:ident, prefix {$($prefix:ident: $prefix_ty:ty),* $(,)?}, {$($field:ident: $ty:ty),+ $(,)?}) => {
        #[pyclass(module="fixed_projection._native")]
        pub struct $name { pub(crate) inner: ck::$name }
        paste::paste! {
            #[pymethods]
            impl $name {
                #[new]
                #[pyo3(signature=(*, atom_maps=None, $($prefix=ck::$name::default().$prefix,)* $($field=ck::$name::default().$field),+))]
                fn construct(atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>, $($prefix: $prefix_ty,)* $($field: $ty),+) -> Self {
                    Self { inner: ck::$name {
                        atom_maps: atom_maps.map(|maps| maps.into_iter().map(|row| row.into_iter().map(|value| value.inner).collect()).collect()).unwrap_or_else(|| ck::$name::default().atom_maps),
                        $($prefix,)* $($field),+
                    } }
                }
                #[staticmethod]
                #[pyo3(signature=($($prefix=ck::$name::default().$prefix,)* atom_maps=None, $($field=ck::$name::default().$field),+))]
                fn new($($prefix: $prefix_ty,)* atom_maps: Option<Vec<Vec<AlignmentAtomMap>>>, $($field: $ty),+) -> Self {
                    let atom_maps = atom_maps.map(|maps| maps.into_iter().map(|row| row.into_iter().map(|value| value.inner).collect()).collect()).unwrap_or_else(|| ck::$name::default().atom_maps);
                    Self { inner: ck::$name::new($($prefix,)* atom_maps, $($field),+) }
                }
                #[getter]
                fn atom_maps(&self) -> Vec<Vec<AlignmentAtomMap>> {
                    self.inner.atom_maps.iter().map(|row| row.iter().copied().map(|inner| AlignmentAtomMap { inner }).collect()).collect()
                }
                #[setter]
                fn set_atom_maps(&mut self, value: Option<Vec<Vec<AlignmentAtomMap>>>) {
                    self.inner.atom_maps = value.map(|maps| maps.into_iter().map(|row| row.into_iter().map(|value| value.inner).collect()).collect()).unwrap_or_else(|| ck::$name::default().atom_maps);
                }
                $(#[getter] fn $prefix(&self) -> $prefix_ty { self.inner.$prefix.clone() }
                #[setter] fn [<set_ $prefix>](&mut self, value: $prefix_ty) { self.inner.$prefix = value; })*
                $(#[getter] fn $field(&self) -> $ty { self.inner.$field.clone() }
                #[setter] fn [<set_ $field>](&mut self, value: $ty) { self.inner.$field = value; })+
            }
        }
    };
}
map_parameters!(BestAlignmentParameters, prefix { probe_conformer_id: i32, reference_conformer_id: i32 }, {
    weights: Option<Vec<f64>>,
    reflect: bool, max_iterations: u32, max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool, ignore_hydrogens: bool, num_threads: i32,
});
map_parameters!(CoordinateRmsdParameters, prefix { probe_conformer_id: i32, reference_conformer_id: i32 }, {
    weights: Option<Vec<f64>>,
    max_matches: i32, symmetrize_conjugated_terminal_groups: bool,
});
map_parameters!(AllConformerRmsdParameters, prefix {}, {
    weights: Option<Vec<f64>>, max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool, ignore_hydrogens: bool, num_threads: i32,
});

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct AlignmentTransform {
    pub(crate) inner: ck::AlignmentTransform,
}
#[pymethods]
impl AlignmentTransform {
    #[new]
    fn new(matrix: [[f64; 4]; 4]) -> Self {
        Self {
            inner: ck::AlignmentTransform { matrix },
        }
    }
    fn matrix(&self) -> [[f64; 4]; 4] {
        *self.inner.matrix()
    }
}
#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct AlignmentResult {
    pub(crate) inner: ck::AlignmentResult,
}
#[pymethods]
impl AlignmentResult {
    #[new]
    fn new(rmsd: f64, transform: AlignmentTransform, atom_map: Vec<AlignmentAtomMap>) -> Self {
        Self {
            inner: ck::AlignmentResult {
                rmsd,
                transform: transform.inner,
                atom_map: atom_map.into_iter().map(|value| value.inner).collect(),
            },
        }
    }
    fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
    fn transform(&self) -> AlignmentTransform {
        AlignmentTransform {
            inner: *self.inner.transform(),
        }
    }
    fn atom_map(&self) -> Vec<AlignmentAtomMap> {
        self.inner
            .atom_map()
            .iter()
            .copied()
            .map(|inner| AlignmentAtomMap { inner })
            .collect()
    }
}
#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone, Copy)]
pub struct ConformerRmsd {
    pub(crate) inner: ck::ConformerRmsd,
}
#[pymethods]
impl ConformerRmsd {
    #[new]
    fn new(probe_conformer_id: usize, reference_conformer_id: usize, rmsd: f64) -> Self {
        Self {
            inner: ck::ConformerRmsd {
                probe_conformer_id,
                reference_conformer_id,
                rmsd,
            },
        }
    }
    fn probe_conformer_id(&self) -> usize {
        self.inner.probe_conformer_id()
    }
    fn reference_conformer_id(&self) -> usize {
        self.inner.reference_conformer_id()
    }
    fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
}
#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct ConformerAlignmentReport {
    pub(crate) inner: ck::ConformerAlignmentReport,
}
#[pymethods]
impl ConformerAlignmentReport {
    #[new]
    fn new(rmsds: Vec<f64>) -> Self {
        Self {
            inner: ck::ConformerAlignmentReport { rmsds },
        }
    }
    fn rmsds(&self) -> &[f64] {
        self.inner.rmsds()
    }
}
