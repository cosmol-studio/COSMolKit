use crate::collection_params::AlignmentAtomMap;
use cosmolkit as ck;
use pyo3::prelude::*;

ck::project_configuration!(ck::BestAlignmentParameters);
ck::project_configuration!(ck::CoordinateRmsdParameters);
ck::project_configuration!(ck::AllConformerRmsdParameters);

#[pyclass(module = "cosmolkit", frozen, from_py_object)]
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
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
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
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
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
#[pyclass(module = "cosmolkit", frozen, from_py_object)]
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
