//! Alignment projections of the canonical public facade. No chemistry lives here.
use ::cosmolkit as ck;
use pyo3::{exceptions::PyValueError, prelude::*};
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pymethods};
use pyo3_stub_gen_derive::remove_gen_stub;

fn check_projected_weight_count(
    weights: Option<&[f64]>,
    expected: usize,
) -> Result<(), ck::AlignmentError> {
    // RDKit✔️✔️:   if (wtsVec) {
    // RDKit✔️✔️:     if (wtsVec->size() != nAtms) {
    // RDKit✔️✔️:       throw_value_error("Incorrect number of weights specified");
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    if let Some(weights) = weights {
        if weights.len() != expected {
            return Err(ck::AlignmentError::WeightCountMismatch {
                map_len: expected,
                weight_len: weights.len(),
            });
        }
    }
    Ok(())
}

// Source-wrapper projections preserve the parameter object's original [] fields.
fn projected_weights(weights: &Option<Vec<f64>>) -> Option<Vec<f64>> {
    // RDKit✔️✔️: RDNumeric::DoubleVector *translateDoubleSeq(const python::object &doubleSeq) {
    // RDKit✔️✔️:   PySequenceHolder<double> doubles(doubleSeq);
    // RDKit✔️✔️:   unsigned int nDoubles = doubles.size();
    // RDKit✔️✔️:   RDNumeric::DoubleVector *doubleVec;
    // RDKit✔️✔️:   doubleVec = nullptr;
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   if (nDoubles > 0) {
    // RDKit✔️✔️:     doubleVec = new RDNumeric::DoubleVector(nDoubles);
    // RDKit✔️✔️:     for (i = 0; i < nDoubles; ++i) {
    // RDKit✔️✔️:       doubleVec->setVal(i, doubles[i]);
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return doubleVec;
    // RDKit✔️✔️: }
    weights
        .as_ref()
        .filter(|values| !values.is_empty())
        .cloned()
}
fn projected_atom_map(map: &Option<Vec<PyAlignmentAtomMap>>) -> Option<Vec<ck::AlignmentAtomMap>> {
    // RDKit✔️✔️: std::vector<std::pair<int, int>> *translateAtomMap(
    // RDKit✔️✔️:     const python::object &atomMap) {
    // RDKit✔️✔️:   PySequenceHolder<python::object> pyAtomMap(atomMap);
    // RDKit✔️✔️:   std::vector<std::pair<int, int>> *res;
    // RDKit✔️✔️:   res = nullptr;
    // RDKit✔️✔️:   unsigned int i;
    // RDKit✔️✔️:   unsigned int n = pyAtomMap.size();
    // RDKit✔️✔️:   if (n > 0) {
    // RDKit✔️✔️:     res = new std::vector<std::pair<int, int>>;
    // RDKit✔️✔️:     for (i = 0; i < n; ++i) {
    // RDKit✔️✔️:       PySequenceHolder<int> item(pyAtomMap[i]);
    // RDKit✔️✔️:       if (item.size() != 2) {
    // RDKit✔️✔️:         delete res;
    // RDKit✔️✔️:         res = nullptr;
    // RDKit✔️✔️:         throw_value_error("Incorrect format for an atomMap");
    // RDKit✔️✔️:       }
    // RDKit✔️✔️:       res->push_back(std::pair<int, int>(item[0], item[1]));
    // RDKit✔️✔️:     }
    // RDKit✔️✔️:   }
    // RDKit✔️✔️:   return res;
    // RDKit✔️✔️: }
    map.as_ref()
        .filter(|values| !values.is_empty())
        .map(|values| values.iter().map(Into::into).collect())
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "AlignmentAtomMap",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "A probe-to-reference atom pair used by molecular alignment."]
pub(crate) struct PyAlignmentAtomMap {
    probe_atom: usize,
    reference_atom: usize,
}

impl From<ck::AlignmentAtomMap> for PyAlignmentAtomMap {
    fn from(value: ck::AlignmentAtomMap) -> Self {
        Self {
            probe_atom: value.probe_atom,
            reference_atom: value.reference_atom,
        }
    }
}

impl From<&PyAlignmentAtomMap> for ck::AlignmentAtomMap {
    fn from(value: &PyAlignmentAtomMap) -> Self {
        Self {
            probe_atom: value.probe_atom,
            reference_atom: value.reference_atom,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentAtomMap {
    #[new]
    fn new(probe_atom: usize, reference_atom: usize) -> Self {
        Self {
            probe_atom,
            reference_atom,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentAtomMap(probe_atom={}, reference_atom={})",
            self.probe_atom, self.reference_atom
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "AlignmentParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters for an explicit or first-match molecular alignment."]
pub(crate) struct PyAlignmentParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_map: Option<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
}

impl PyAlignmentParameters {
    pub(crate) fn wrapper_parameters(
        &self,
        atom_count: usize,
    ) -> Result<ck::AlignmentParameters, ck::AlignmentError> {
        // RDKit✔️✔️: PyObject *getMolAlignTransform(const ROMol &prbMol, const ROMol &refMol,
        // RDKit✔️✔️:                                int prbCid = -1, int refCid = -1,
        // RDKit✔️✔️:                                python::object atomMap = python::list(),
        // RDKit✔️✔️:                                python::object weights = python::list(),
        // RDKit✔️✔️:                                bool reflect = false,
        // RDKit✔️✔️:                                unsigned int maxIters = 50) {
        // RDKit✔️✔️:   std::unique_ptr<MatchVectType> aMap(translateAtomMap(atomMap));
        // RDKit✔️✔️:   unsigned int nAtms;
        // RDKit✔️✔️:   if (aMap) {
        // RDKit✔️✔️:     nAtms = aMap->size();
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     nAtms = prbMol.getNumAtoms();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
        // RDKit✔️✔️:   if (wtsVec) {
        // RDKit✔️✔️:     if (wtsVec->size() != nAtms) {
        // RDKit✔️✔️:       throw_value_error("Incorrect number of weights specified");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   RDGeom::Transform3D trans;
        // RDKit✔️✔️:   double rmsd;
        // RDKit✔️✔️:   {
        // RDKit✔️✔️:     NOGIL gil;
        // RDKit✔️✔️:     rmsd = MolAlign::getAlignmentTransform(prbMol, refMol, trans, prbCid,
        // RDKit✔️✔️:                                            refCid, aMap.get(), wtsVec.get(),
        // RDKit✔️✔️:                                            reflect, maxIters);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   return generateRmsdTransMatchPyTuple(rmsd, trans);
        // RDKit✔️✔️: }
        // RDKit✔️✔️: double AlignMolecule(ROMol &prbMol, const ROMol &refMol, int prbCid = -1,
        // RDKit✔️✔️:                      int refCid = -1, python::object atomMap = python::list(),
        // RDKit✔️✔️:                      python::object weights = python::list(),
        // RDKit✔️✔️:                      bool reflect = false, unsigned int maxIters = 50) {
        // RDKit✔️✔️:   std::unique_ptr<MatchVectType> aMap(translateAtomMap(atomMap));
        // RDKit✔️✔️:   unsigned int nAtms;
        // RDKit✔️✔️:   if (aMap) {
        // RDKit✔️✔️:     nAtms = aMap->size();
        // RDKit✔️✔️:   } else {
        // RDKit✔️✔️:     nAtms = prbMol.getNumAtoms();
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   std::unique_ptr<RDNumeric::DoubleVector> wtsVec(translateDoubleSeq(weights));
        // RDKit✔️✔️:   if (wtsVec) {
        // RDKit✔️✔️:     if (wtsVec->size() != nAtms) {
        // RDKit✔️✔️:       throw_value_error("Incorrect number of weights specified");
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:
        // RDKit✔️✔️:   double rmsd;
        // RDKit✔️✔️:   {
        // RDKit✔️✔️:     NOGIL gil;
        // RDKit✔️✔️:     rmsd = MolAlign::alignMol(prbMol, refMol, prbCid, refCid, aMap.get(),
        // RDKit✔️✔️:                               wtsVec.get(), reflect, maxIters);
        // RDKit✔️✔️:   }
        // RDKit✔️✔️:   return rmsd;
        // RDKit✔️✔️: }
        let params = self.core_parameters();
        let expected = params.atom_map.as_ref().map_or(atom_count, Vec::len);
        check_projected_weight_count(params.weights.as_deref(), expected)?;
        Ok(params)
    }

    pub(crate) fn core_parameters(&self) -> ck::AlignmentParameters {
        ck::AlignmentParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_map: projected_atom_map(&self.atom_map),
            weights: projected_weights(&self.weights),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_map=None, weights=None, reflect=false, max_iterations=50))]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_map: Option<Vec<PyAlignmentAtomMap>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_map,
            weights,
            reflect,
            max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "BestAlignmentParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters for source-compatible best molecular alignment and RMSD."]
pub(crate) struct PyBestAlignmentParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
    num_threads: i32,
}

impl PyBestAlignmentParameters {
    pub(crate) fn wrapper_parameters(
        &self,
    ) -> Result<ck::BestAlignmentParameters, ck::AlignmentError> {
        // RDKit✔️✔️:   pyBestAlignmentParams(int maxMatches_, bool symmetrizeTerminalGroups_,
        // RDKit✔️✔️:                         bool ignoreHs_, int numThreads_, python::object map_,
        // RDKit✔️✔️:                         python::object weights_)
        // RDKit✔️✔️:       : BestAlignmentParams{maxMatches_, symmetrizeTerminalGroups_, ignoreHs_,
        // RDKit✔️✔️:                             numThreads_, std::vector<MatchVectType>(),
        // RDKit✔️✔️:                             nullptr} {
        // RDKit✔️✔️:     unsigned int nAtms = 0;
        // RDKit✔️✔️:     if (map_ != python::object()) {
        // RDKit✔️✔️:       map = translateAtomMapSeq(map_);
        // RDKit✔️✔️:       if (!map.empty()) {
        // RDKit✔️✔️:         nAtms = map.front().size();
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     wtsVec.reset(translateDoubleSeq(weights_));
        // RDKit✔️✔️:     if (wtsVec) {
        // RDKit✔️✔️:       if (!map.empty() && wtsVec->size() != nAtms) {
        // RDKit✔️✔️:         throw_value_error("Incorrect number of weights specified");
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       weights = wtsVec.get();
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        let params = self.core_parameters();
        if let Some(first) = params.atom_maps.first() {
            check_projected_weight_count(params.weights.as_deref(), first.len())?;
        }
        Ok(params)
    }

    pub(crate) fn core_parameters(&self) -> ck::BestAlignmentParameters {
        ck::BestAlignmentParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: projected_weights(&self.weights),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
            ignore_hydrogens: self.ignore_hydrogens,
            num_threads: self.num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyBestAlignmentParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_maps=None, weights=None, reflect=false, max_iterations=50, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true, ignore_hydrogens=true, num_threads=1))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            reflect,
            max_iterations,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "AllConformerRmsdParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters accepted by source-compatible all-conformer best RMSD."]
pub(crate) struct PyAllConformerRmsdParameters {
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
    ignore_hydrogens: bool,
    num_threads: i32,
}

impl PyAllConformerRmsdParameters {
    pub(crate) fn wrapper_parameters(
        &self,
    ) -> Result<ck::AllConformerRmsdParameters, ck::AlignmentError> {
        // RDKit✔️✔️:   pyBestAlignmentParams(int maxMatches_, bool symmetrizeTerminalGroups_,
        // RDKit✔️✔️:                         bool ignoreHs_, int numThreads_, python::object map_,
        // RDKit✔️✔️:                         python::object weights_)
        // RDKit✔️✔️:       : BestAlignmentParams{maxMatches_, symmetrizeTerminalGroups_, ignoreHs_,
        // RDKit✔️✔️:                             numThreads_, std::vector<MatchVectType>(),
        // RDKit✔️✔️:                             nullptr} {
        // RDKit✔️✔️:     unsigned int nAtms = 0;
        // RDKit✔️✔️:     if (map_ != python::object()) {
        // RDKit✔️✔️:       map = translateAtomMapSeq(map_);
        // RDKit✔️✔️:       if (!map.empty()) {
        // RDKit✔️✔️:         nAtms = map.front().size();
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:
        // RDKit✔️✔️:     wtsVec.reset(translateDoubleSeq(weights_));
        // RDKit✔️✔️:     if (wtsVec) {
        // RDKit✔️✔️:       if (!map.empty() && wtsVec->size() != nAtms) {
        // RDKit✔️✔️:         throw_value_error("Incorrect number of weights specified");
        // RDKit✔️✔️:       }
        // RDKit✔️✔️:       weights = wtsVec.get();
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:   }
        let params = self.core_parameters();
        if let Some(first) = params.atom_maps.first() {
            check_projected_weight_count(params.weights.as_deref(), first.len())?;
        }
        Ok(params)
    }

    pub(crate) fn core_parameters(&self) -> ck::AllConformerRmsdParameters {
        ck::AllConformerRmsdParameters {
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: projected_weights(&self.weights),
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
            ignore_hydrogens: self.ignore_hydrogens,
            num_threads: self.num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAllConformerRmsdParameters {
    #[new]
    #[pyo3(signature = (atom_maps=None, weights=None, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true, ignore_hydrogens=true, num_threads=1))]
    fn new(
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
        ignore_hydrogens: bool,
        num_threads: i32,
    ) -> Self {
        Self {
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
            ignore_hydrogens,
            num_threads,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "CoordinateRmsdParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters for RMSD measurement in the existing coordinate frame."]
pub(crate) struct PyCoordinateRmsdParameters {
    probe_conformer_id: i32,
    reference_conformer_id: i32,
    atom_maps: Vec<Vec<PyAlignmentAtomMap>>,
    weights: Option<Vec<f64>>,
    max_matches: i32,
    symmetrize_conjugated_terminal_groups: bool,
}

impl PyCoordinateRmsdParameters {
    pub(crate) fn core_parameters(&self) -> ck::CoordinateRmsdParameters {
        ck::CoordinateRmsdParameters {
            probe_conformer_id: self.probe_conformer_id,
            reference_conformer_id: self.reference_conformer_id,
            atom_maps: self
                .atom_maps
                .iter()
                .map(|map| map.iter().map(Into::into).collect())
                .collect(),
            weights: projected_weights(&self.weights),
            max_matches: self.max_matches,
            symmetrize_conjugated_terminal_groups: self.symmetrize_conjugated_terminal_groups,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyCoordinateRmsdParameters {
    #[new]
    #[pyo3(signature = (probe_conformer_id=-1, reference_conformer_id=-1, atom_maps=None, weights=None, max_matches=1_000_000, symmetrize_conjugated_terminal_groups=true))]
    fn new(
        probe_conformer_id: i32,
        reference_conformer_id: i32,
        atom_maps: Option<Vec<Vec<PyAlignmentAtomMap>>>,
        weights: Option<Vec<f64>>,
        max_matches: i32,
        symmetrize_conjugated_terminal_groups: bool,
    ) -> Self {
        Self {
            probe_conformer_id,
            reference_conformer_id,
            atom_maps: atom_maps.unwrap_or_default(),
            weights,
            max_matches,
            symmetrize_conjugated_terminal_groups,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "ConformerAlignmentParameters",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone)]
#[doc = "Parameters for aligning selected or all conformers of one molecule."]
pub(crate) struct PyConformerAlignmentParameters {
    atom_indices: Option<Vec<usize>>,
    conformer_ids: Option<Vec<usize>>,
    weights: Option<Vec<f64>>,
    reflect: bool,
    max_iterations: u32,
}

impl PyConformerAlignmentParameters {
    pub(crate) fn core_parameters(&self) -> ck::ConformerAlignmentParameters {
        ck::ConformerAlignmentParameters {
            atom_indices: self.atom_indices.clone(),
            conformer_ids: self.conformer_ids.clone(),
            weights: projected_weights(&self.weights),
            reflect: self.reflect,
            max_iterations: self.max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerAlignmentParameters {
    #[new]
    #[pyo3(signature = (atom_indices=None, conformer_ids=None, weights=None, reflect=false, max_iterations=50))]
    fn new(
        atom_indices: Option<Vec<usize>>,
        conformer_ids: Option<Vec<usize>>,
        weights: Option<Vec<f64>>,
        reflect: bool,
        max_iterations: u32,
    ) -> Self {
        Self {
            atom_indices,
            conformer_ids,
            weights,
            reflect,
            max_iterations,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", name = "AlignmentTransform", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlignmentTransform {
    matrix: [[f64; 4]; 4],
}

impl From<ck::AlignmentTransform> for PyAlignmentTransform {
    fn from(value: ck::AlignmentTransform) -> Self {
        Self {
            matrix: value.matrix,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentTransform {
    #[doc = "Return the source-oriented 4x4 homogeneous transform matrix."]
    fn matrix(&self) -> Vec<Vec<f64>> {
        self.matrix.iter().map(|row| row.to_vec()).collect()
    }

    fn __repr__(&self) -> String {
        "AlignmentTransform(matrix=4x4)".to_string()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", name = "AlignmentResult", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyAlignmentResult {
    inner: ck::AlignmentResult,
}

impl From<ck::AlignmentResult> for PyAlignmentResult {
    fn from(inner: ck::AlignmentResult) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyAlignmentResult {
    #[doc = "Return the aligned root-mean-square deviation."]
    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    #[doc = "Return the transform mapping probe coordinates onto the reference."]
    fn transform(&self) -> PyAlignmentTransform {
        self.inner.transform.into()
    }

    #[doc = "Return the selected probe-to-reference atom map."]
    fn atom_map(&self) -> Vec<PyAlignmentAtomMap> {
        self.inner
            .atom_map
            .iter()
            .copied()
            .map(Into::into)
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "AlignmentResult(rmsd={}, mapped_atoms={})",
            self.inner.rmsd,
            self.inner.atom_map.len()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", name = "ConformerRmsd", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PyConformerRmsd {
    inner: ck::ConformerRmsd,
}

impl From<ck::ConformerRmsd> for PyConformerRmsd {
    fn from(inner: ck::ConformerRmsd) -> Self {
        Self { inner }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerRmsd {
    fn probe_conformer_id(&self) -> usize {
        self.inner.probe_conformer_id
    }

    fn reference_conformer_id(&self) -> usize {
        self.inner.reference_conformer_id
    }

    fn rmsd(&self) -> f64 {
        self.inner.rmsd
    }

    fn __repr__(&self) -> String {
        format!(
            "ConformerRmsd(probe_conformer_id={}, reference_conformer_id={}, rmsd={})",
            self.inner.probe_conformer_id, self.inner.reference_conformer_id, self.inner.rmsd
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(
    module = "cosmolkit",
    name = "ConformerAlignmentReport",
    skip_from_py_object
)]
#[derive(Clone)]
pub(crate) struct PyConformerAlignmentReport {
    rmsds: Vec<f64>,
}

impl From<ck::ConformerAlignmentReport> for PyConformerAlignmentReport {
    fn from(value: ck::ConformerAlignmentReport) -> Self {
        Self { rmsds: value.rmsds }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), remove_gen_stub)]
#[pymethods]
impl PyConformerAlignmentReport {
    #[doc = "Return RMS values in the source conformer-selection order."]
    fn rmsds(&self) -> Vec<f64> {
        self.rmsds.clone()
    }

    fn __repr__(&self) -> String {
        format!("ConformerAlignmentReport(rmsds={})", self.rmsds.len())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("AlignmentError", module.py().get_type::<AlignmentError>())?;
    module.add_class::<PyAlignmentAtomMap>()?;
    module.add_class::<PyAlignmentParameters>()?;
    module.add_class::<PyBestAlignmentParameters>()?;
    module.add_class::<PyAllConformerRmsdParameters>()?;
    module.add_class::<PyCoordinateRmsdParameters>()?;
    module.add_class::<PyConformerAlignmentParameters>()?;
    module.add_class::<PyAlignmentTransform>()?;
    module.add_class::<PyAlignmentResult>()?;
    module.add_class::<PyConformerRmsd>()?;
    module.add_class::<PyConformerAlignmentReport>()?;
    Ok(())
}

pyo3::create_exception!(cosmolkit, AlignmentError, PyValueError);
pub(crate) fn alignment_pyerr(py: Python<'_>, source: ck::AlignmentError) -> PyErr {
    use ck::AlignmentError as E;
    let kind = match &source {
        E::InvalidCoordinates(..) => "InvalidCoordinates",
        E::Matching(..) => "Matching",
        E::QueryGraph(..) => "QueryGraph",
        E::QueryParse(..) => "QueryParse",
        E::ThreadSelection(..) => "ThreadSelection",
        E::NoConformers => "NoConformers",
        E::ConformerNotFound { .. } => "ConformerNotFound",
        E::EmptyAtomMap => "EmptyAtomMap",
        E::ProbeAtomOutOfRange { .. } => "ProbeAtomOutOfRange",
        E::ReferenceAtomOutOfRange { .. } => "ReferenceAtomOutOfRange",
        E::WeightCountMismatch { .. } => "WeightCountMismatch",
        E::NonPositiveWeight { .. } => "NonPositiveWeight",
        E::NoSubstructureMatch => "NoSubstructureMatch",
        E::TerminalGroupSymmetrization { .. } => "TerminalGroupSymmetrization",
        E::NumericalPrecondition { .. } => "NumericalPrecondition",
        E::WorkerTerminated => "WorkerTerminated",
    };
    let error = AlignmentError::new_err(source.to_string());
    let attributes = || -> PyResult<()> {
        let v = error.value(py);
        v.setattr("domain", "alignment")?;
        v.setattr("kind", kind)?;
        match &source {
            E::ConformerNotFound { id } => v.setattr("id", *id)?,
            E::ProbeAtomOutOfRange { index, atom_count }
            | E::ReferenceAtomOutOfRange { index, atom_count } => {
                v.setattr("index", *index)?;
                v.setattr("atom_count", *atom_count)?;
            }
            E::WeightCountMismatch {
                map_len,
                weight_len,
            } => {
                v.setattr("map_len", *map_len)?;
                v.setattr("weight_len", *weight_len)?;
            }
            E::NonPositiveWeight { index } => v.setattr("index", *index)?,
            E::TerminalGroupSymmetrization { message } | E::NumericalPrecondition { message } => {
                v.setattr("message", *message)?
            }
            _ => {}
        }
        Ok(())
    };
    if let Err(attribute_error) = attributes() {
        return attribute_error;
    }
    error.set_cause(
        py,
        std::error::Error::source(&source)
            .map(|cause| crate::drawing_binding::source_pyerr(py, cause)),
    );
    error
}
