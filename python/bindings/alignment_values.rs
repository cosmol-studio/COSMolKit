use crate::collection_params::AlignmentAtomMap;
use cosmolkit as ck;
use pyo3::prelude::*;

ck::project_configuration!(ck::BestAlignmentParameters);
ck::project_configuration!(ck::CoordinateRmsdParameters);
ck::project_configuration!(ck::AllConformerRmsdParameters);

crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// Readonly homogeneous 4-by-4 transform mapping aligned probe coordinates into the reference frame.
    pub struct AlignmentTransform {
        pub(crate) inner: ck::AlignmentTransform,
    }
}
#[pymethods]
impl AlignmentTransform {
    #[new]
    fn new(matrix: [[f64; 4]; 4]) -> Self {
        Self {
            inner: ck::AlignmentTransform { matrix },
        }
    }
    #[getter]
    /// Homogeneous 4-by-4 transformation matrix.
    fn matrix(&self) -> [[f64; 4]; 4] {
        *self.inner.matrix()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Alignment result retaining RMSD, transformation and probe/reference atom correspondence.
    pub struct AlignmentResult {
        pub(crate) inner: ck::AlignmentResult,
    }
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
    #[getter]
    /// Root-mean-square deviation for the associated alignment or conformer pair.
    fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
    #[getter]
    /// Transform determined by the alignment.
    fn transform(&self) -> AlignmentTransform {
        AlignmentTransform {
            inner: *self.inner.transform(),
        }
    }
    #[getter]
    /// Independent list of probe-to-reference atom correspondences used for alignment.
    fn atom_map(&self) -> Vec<AlignmentAtomMap> {
        self.inner
            .atom_map()
            .iter()
            .copied()
            .map(|inner| AlignmentAtomMap { inner })
            .collect()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone, Copy)]
    /// RMSD value associated with explicit probe and reference conformer identifiers.
    pub struct ConformerRmsd {
        pub(crate) inner: ck::ConformerRmsd,
    }
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
    #[getter]
    /// Identifier of the probe conformer.
    fn probe_conformer_id(&self) -> usize {
        self.inner.probe_conformer_id()
    }
    #[getter]
    /// Identifier of the reference conformer.
    fn reference_conformer_id(&self) -> usize {
        self.inner.reference_conformer_id()
    }
    #[getter]
    /// Root-mean-square deviation for the associated alignment or conformer pair.
    fn rmsd(&self) -> f64 {
        self.inner.rmsd()
    }
}
crate::configuration::type_projection! {
    #[pyclass(module = "cosmolkit", frozen, from_py_object)]
    #[derive(Clone)]
    /// Ordered RMSD results returned by native conformer alignment.
    pub struct ConformerAlignmentReport {
        pub(crate) inner: ck::ConformerAlignmentReport,
    }
}
#[pymethods]
impl ConformerAlignmentReport {
    #[new]
    fn new(rmsds: Vec<f64>) -> Self {
        Self {
            inner: ck::ConformerAlignmentReport { rmsds },
        }
    }
    #[getter]
    /// RMSD values in the native alignment result order.
    fn rmsds(&self) -> &[f64] {
        self.inner.rmsds()
    }
}
