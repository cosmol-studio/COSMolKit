//! Thin language projections of canonical TAU configuration and validated results.
use crate::drawing_binding::Molecule;
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::sync::{Arc, Mutex};

pyo3::create_exception!(cosmolkit, TautomerRunError, PyValueError);
pyo3::create_exception!(cosmolkit, TautomerCatalogError, PyValueError);

pub(crate) fn run_pyerr(py: Python<'_>, source: &ck::TautomerRunError) -> PyErr {
    use ck::TautomerRunError as E;
    let kind = match source {
        E::Smiles(..) => "Smiles",
        E::Valence(..) => "Valence",
        E::Rings(..) => "Rings",
        E::Kekulize(..) => "Kekulize",
        E::Sanitize(..) => "Sanitize",
        E::LegacyStereo(..) => "LegacyStereo",
        E::QueryCompile(..) => "QueryCompile",
        E::Match(..) => "Match",
        E::MatchContext(..) => "MatchContext",
        E::Topology(..) => "Topology",
        E::AtomProperty(..) => "AtomProperty",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::BondValue(..) => "BondValue",
        E::PropertyString(..) => "PropertyString",
        E::AtomMappingCount { .. } => "AtomMappingCount",
        E::BondMappingCount { .. } => "BondMappingCount",
        E::AtomCountMismatch { .. } => "AtomCountMismatch",
        E::BondCountMismatch { .. } => "BondCountMismatch",
        E::AtomOutOfRange { .. } => "AtomOutOfRange",
        E::BondOutOfRange { .. } => "BondOutOfRange",
        E::EmptyTransformMatch => "EmptyTransformMatch",
        E::EditCount { .. } => "EditCount",
        E::HydrogenCountOutOfRange { .. } => "HydrogenCountOutOfRange",
        E::FormalChargeOutOfRange { .. } => "FormalChargeOutOfRange",
        E::MissingMappedBond { .. } => "MissingMappedBond",
        E::MissingValence { .. } => "MissingValence",
        E::NoCanonicalTautomer => "NoCanonicalTautomer",
        E::Control(..) => "Control",
        E::Callback(..) => "Callback",
    };
    crate::canonical_values::annotate(
        py,
        TautomerRunError::new_err(source.to_string()),
        "tautomer",
        kind,
        source,
    )
}
fn catalog_pyerr(py: Python<'_>, source: ck::TautomerCatalogError) -> PyErr {
    use ck::TautomerCatalogError as E;
    let kind = match &source {
        E::BadInputFile { .. } => "BadInputFile",
        E::BadStreamContents(..) => "BadStreamContents",
        E::InvalidUtf8(..) => "InvalidUtf8",
        E::Transform(..) => "Transform",
        E::TransformIndexOutOfRange { .. } => "TransformIndexOutOfRange",
        E::DeserializationUnderConstruction => "DeserializationUnderConstruction",
    };
    crate::canonical_values::annotate(
        py,
        TautomerCatalogError::new_err(source.to_string()),
        "tautomer_catalog",
        kind,
        &source,
    )
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum TautomerEnumerationStatus {
    Completed = 0,
    MaxTautomersReached = 1,
    MaxTransformsReached = 2,
    Canceled = 3,
}
impl From<ck::TautomerEnumerationStatus> for TautomerEnumerationStatus {
    fn from(value: ck::TautomerEnumerationStatus) -> Self {
        match value {
            ck::TautomerEnumerationStatus::Completed => Self::Completed,
            ck::TautomerEnumerationStatus::MaxTautomersReached => Self::MaxTautomersReached,
            ck::TautomerEnumerationStatus::MaxTransformsReached => Self::MaxTransformsReached,
            ck::TautomerEnumerationStatus::Canceled => Self::Canceled,
        }
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct TautomerScoreTerm {
    inner: ck::TautomerScoreTerm,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScoreTerm {
    #[new]
    fn py_new(name: String, smarts: String, score: i32) -> Self {
        Self::new(name, smarts, score)
    }
    #[staticmethod]
    fn new(name: String, smarts: String, score: i32) -> Self {
        Self {
            inner: ck::TautomerScoreTerm::new(name, smarts, score),
        }
    }
    fn name(&self) -> &str {
        self.inner.name()
    }
    fn smarts(&self) -> &str {
        self.inner.smarts()
    }
    fn score(&self) -> i32 {
        self.inner.score()
    }
    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct TautomerScoreParams {
    pub(crate) inner: ck::TautomerScoreParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScoreParams {
    #[new]
    #[pyo3(signature=(terms=None))]
    fn new(terms: Option<Vec<PyRef<'_, TautomerScoreTerm>>>) -> Self {
        Self {
            inner: ck::TautomerScoreParams {
                terms: terms.map(|terms| terms.iter().map(|term| term.inner.clone()).collect()),
            },
        }
    }
    #[getter]
    fn terms(&self) -> Option<Vec<TautomerScoreTerm>> {
        self.inner.terms.as_ref().map(|terms| {
            terms
                .iter()
                .cloned()
                .map(|inner| TautomerScoreTerm { inner })
                .collect()
        })
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TautomerScore {
    pub(crate) inner: ck::TautomerScore,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScore {
    fn ring(&self) -> i32 {
        self.inner.ring()
    }
    fn substructure(&self) -> i32 {
        self.inner.substructure()
    }
    fn hetero_hydrogen(&self) -> i32 {
        self.inner.hetero_hydrogen()
    }
    fn total(&self) -> i32 {
        self.inner.total()
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TautomerMoleculeView {
    inner: ck::TautomerMoleculeView<'static>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerMoleculeView {
    fn properties(&self) -> crate::canonical_property_values::MoleculeProperties {
        crate::canonical_property_values::MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }
    fn to_owned(&self) -> Self {
        Self {
            inner: self.inner.to_owned(),
        }
    }
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    fn atoms(&self) -> Vec<crate::canonical_atom_bond::Atom> {
        let metadata = self.inner.atom_metadata();
        self.inner
            .atoms()
            .iter()
            .enumerate()
            .map(|(i, atom)| crate::canonical_atom_bond::Atom {
                inner: atom.clone(),
                degree: self
                    .inner
                    .atom_degree(ck::AtomId::new(i))
                    .expect("existing atom"),
                metadata: metadata
                    .as_ref()
                    .map(|rows| rows[i].clone())
                    .map_err(Clone::clone),
            })
            .collect()
    }
    fn bonds(&self) -> Vec<crate::canonical_atom_bond::Bond> {
        self.inner
            .bonds()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
            .collect()
    }
    fn atom(&self, index: usize) -> Option<crate::canonical_atom_bond::Atom> {
        let inner = self.inner.atom(ck::AtomId::new(index))?.clone();
        Some(crate::canonical_atom_bond::Atom {
            inner,
            degree: self.inner.atom_degree(ck::AtomId::new(index))?,
            metadata: self.inner.atom_metadata().map(|rows| rows[index].clone()),
        })
    }
    fn bond(&self, index: usize) -> Option<crate::canonical_atom_bond::Bond> {
        self.inner
            .bond(ck::BondId::new(index))
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
    }
    fn atom_degree(&self, index: usize) -> Option<usize> {
        self.inner.atom_degree(ck::AtomId::new(index))
    }
    fn atom_metadata(
        &self,
        py: Python<'_>,
    ) -> PyResult<Vec<crate::canonical_atom_bond::AtomMetadata>> {
        self.inner
            .atom_metadata()
            .map(|rows| {
                rows.into_iter()
                    .map(|inner| crate::canonical_atom_bond::AtomMetadata { inner })
                    .collect()
            })
            .map_err(|error| crate::canonical_atom_bond::valence_pyerr(py, error))
    }
    fn to_smiles(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_smiles()
            .map_err(|error| run_pyerr(py, &error))
    }
    fn tautomer_score(&self, py: Python<'_>) -> PyResult<TautomerScore> {
        self.inner
            .tautomer_score()
            .map(|inner| TautomerScore { inner })
            .map_err(|error| run_pyerr(py, &error))
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TautomerProgress {
    inner: ck::TautomerProgress<'static>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerProgress {
    fn to_owned(&self) -> Self {
        Self {
            inner: self.inner.to_owned(),
        }
    }
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn __len__(&self) -> usize {
        self.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status().into()
    }
    fn num_transforms(&self) -> u32 {
        self.inner.num_transforms()
    }
    fn modified_atoms(&self) -> Vec<usize> {
        self.inner
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn modified_bonds(&self) -> Vec<usize> {
        self.inner
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn entries(&self) -> Vec<(String, TautomerMoleculeView)> {
        self.inner
            .entries()
            .map(|(key, value)| {
                (
                    key.to_owned(),
                    TautomerMoleculeView {
                        inner: value.to_owned(),
                    },
                )
            })
            .collect()
    }
}

// Each invocation owns its error transport. Shared parameter values never mix
// errors from concurrent or recursive Python callbacks.
type CallbackFailure = Arc<Mutex<Option<PyErr>>>;
struct PyCallback {
    callable: Py<PyAny>,
    failure: CallbackFailure,
}
struct PyScorer {
    callable: Py<PyAny>,
    failure: CallbackFailure,
}
fn retain_failure(failure: &CallbackFailure, error: PyErr) -> ck::TautomerRunError {
    let detail = error.to_string();
    let mut slot = failure.lock().expect("callback error mutex");
    if slot.is_none() {
        *slot = Some(error)
    }
    ck::TautomerRunError::Callback(detail)
}
impl ck::TautomerEnumerationCallback for PyCallback {
    fn should_continue(
        &self,
        source: ck::TautomerMoleculeView<'_>,
        progress: ck::TautomerProgress<'_>,
    ) -> Result<bool, ck::TautomerRunError> {
        Python::attach(|py| {
            let result = (|| -> PyResult<bool> {
                let source = Py::new(
                    py,
                    TautomerMoleculeView {
                        inner: source.to_owned(),
                    },
                )?;
                let progress = Py::new(
                    py,
                    TautomerProgress {
                        inner: progress.to_owned(),
                    },
                )?;
                // RDKit✔️❌:     return getCallbackOverride()(boost::ref(mol), boost::ref(pyRes));
                // Python retention requires one detached read-only source snapshot;
                // progress owns its entries, as the pinned wrapper copies pyRes.
                self.callable
                    .bind(py)
                    .call1((source, progress))?
                    .is_truthy()
            })();
            result.map_err(|error| retain_failure(&self.failure, error))
        })
    }
}
impl ck::TautomerScorer for PyScorer {
    fn score(&self, molecule: ck::TautomerMoleculeView<'_>) -> Result<i32, ck::TautomerRunError> {
        // RDKit✔️❌:   int operator()(const ROMol &m) {
        // RDKit✔️❌:     return python::extract<int>(dp_obj(boost::ref(m)));
        // RDKit✔️❌:   }
        // The same callable and signed-int extraction are preserved. A detached
        // snapshot adds O(V+E) copying so Python can retain its read-only view.
        Python::attach(|py| {
            let result = (|| -> PyResult<i32> {
                let value = Py::new(
                    py,
                    TautomerMoleculeView {
                        inner: molecule.to_owned(),
                    },
                )?;
                self.callable.bind(py).call1((value,))?.extract()
            })();
            result.map_err(|error| retain_failure(&self.failure, error))
        })
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object)]
pub(crate) struct TautomerParams {
    pub(crate) inner: ck::TautomerParams,
    callback: Option<Py<PyAny>>,
    scorer: Option<Py<PyAny>>,
}
impl TautomerParams {
    fn from_inner(inner: ck::TautomerParams) -> Self {
        Self {
            inner,
            callback: None,
            scorer: None,
        }
    }
    fn cloned(&self, py: Python<'_>) -> Self {
        Self {
            inner: self.inner.clone(),
            callback: self.callback.as_ref().map(|x| x.clone_ref(py)),
            scorer: self.scorer.as_ref().map(|x| x.clone_ref(py)),
        }
    }
    fn invocation(&self, py: Python<'_>) -> (ck::TautomerParams, CallbackFailure) {
        let failure = Arc::new(Mutex::new(None));
        let mut params = self.inner.clone();
        params.set_callback(self.callback.as_ref().map(|callable| {
            Arc::new(PyCallback {
                callable: callable.clone_ref(py),
                failure: failure.clone(),
            }) as Arc<dyn ck::TautomerEnumerationCallback>
        }));
        params.set_scorer(self.scorer.as_ref().map(|callable| {
            Arc::new(PyScorer {
                callable: callable.clone_ref(py),
                failure: failure.clone(),
            }) as Arc<dyn ck::TautomerScorer>
        }));
        (params, failure)
    }
}
fn callable(py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<Option<Py<PyAny>>> {
    if value
        .as_ref()
        .is_some_and(|value| !value.bind(py).is_callable())
    {
        return Err(PyTypeError::new_err("expected a callable or None"));
    }
    Ok(value)
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerParams {
    #[new]
    #[pyo3(signature=(*,max_tautomers=1000,max_transforms=1000,remove_sp3_stereo=true,remove_bond_stereo=true,remove_isotopic_hydrogens=true,reassign_stereo=true,callback=None,scorer=None,score_params=None))]
    fn new(
        py: Python<'_>,
        max_tautomers: u32,
        max_transforms: u32,
        remove_sp3_stereo: bool,
        remove_bond_stereo: bool,
        remove_isotopic_hydrogens: bool,
        reassign_stereo: bool,
        callback: Option<Py<PyAny>>,
        scorer: Option<Py<PyAny>>,
        score_params: Option<PyRef<'_, TautomerScoreParams>>,
    ) -> PyResult<Self> {
        let mut value = Self::from_inner(
            ck::TautomerParams::default()
                .with_max_tautomers(max_tautomers)
                .with_max_transforms(max_transforms)
                .with_remove_sp3_stereo(remove_sp3_stereo)
                .with_remove_bond_stereo(remove_bond_stereo)
                .with_remove_isotopic_hydrogens(remove_isotopic_hydrogens)
                .with_reassign_stereo(reassign_stereo),
        );
        value.callback = callable(py, callback)?;
        value.scorer = callable(py, scorer)?;
        if let Some(params) = score_params {
            value.inner.score_params = params.inner.clone()
        }
        Ok(value)
    }
    #[staticmethod]
    fn v1(py: Python<'_>) -> PyResult<Self> {
        ck::TautomerParams::v1()
            .map(Self::from_inner)
            .map_err(|e| catalog_pyerr(py, e))
    }
    #[staticmethod]
    fn from_transform_data(
        py: Python<'_>,
        data: Vec<(String, String, String, String)>,
    ) -> PyResult<Self> {
        let borrowed: Vec<_> = data
            .iter()
            .map(|(a, b, c, d)| (a.as_str(), b.as_str(), c.as_str(), d.as_str()))
            .collect();
        ck::TautomerParams::from_transform_data(&borrowed)
            .map(Self::from_inner)
            .map_err(|e| catalog_pyerr(py, e))
    }
    #[staticmethod]
    fn from_transform_file(py: Python<'_>, path: String) -> PyResult<Self> {
        ck::TautomerParams::from_transform_file(path)
            .map(Self::from_inner)
            .map_err(|e| catalog_pyerr(py, e))
    }
    fn transform_count(&self) -> usize {
        self.inner.transform_count()
    }
    fn callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.callback.as_ref().map(|x| x.clone_ref(py))
    }
    fn set_callback(&mut self, py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<()> {
        self.callback = callable(py, value)?;
        Ok(())
    }
    fn scorer(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.scorer.as_ref().map(|x| x.clone_ref(py))
    }
    fn set_scorer(&mut self, py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<()> {
        self.scorer = callable(py, value)?;
        Ok(())
    }
    #[getter]
    fn score_params(&self) -> TautomerScoreParams {
        TautomerScoreParams {
            inner: self.inner.score_params.clone(),
        }
    }
    #[setter]
    fn set_score_params(&mut self, value: &TautomerScoreParams) {
        self.inner.score_params = value.inner.clone()
    }
    fn max_tautomers(&self) -> u32 {
        self.inner.max_tautomers()
    }
    fn set_max_tautomers(&mut self, value: u32) {
        self.inner.set_max_tautomers(value)
    }
    fn with_max_tautomers(&self, py: Python<'_>, value: u32) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_max_tautomers(value);
        result
    }
    fn max_transforms(&self) -> u32 {
        self.inner.max_transforms()
    }
    fn set_max_transforms(&mut self, value: u32) {
        self.inner.set_max_transforms(value)
    }
    fn with_max_transforms(&self, py: Python<'_>, value: u32) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_max_transforms(value);
        result
    }
    fn remove_sp3_stereo(&self) -> bool {
        self.inner.remove_sp3_stereo()
    }
    fn set_remove_sp3_stereo(&mut self, value: bool) {
        self.inner.set_remove_sp3_stereo(value)
    }
    fn with_remove_sp3_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_sp3_stereo(value);
        result
    }
    fn remove_bond_stereo(&self) -> bool {
        self.inner.remove_bond_stereo()
    }
    fn set_remove_bond_stereo(&mut self, value: bool) {
        self.inner.set_remove_bond_stereo(value)
    }
    fn with_remove_bond_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_bond_stereo(value);
        result
    }
    fn remove_isotopic_hydrogens(&self) -> bool {
        self.inner.remove_isotopic_hydrogens()
    }
    fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        self.inner.set_remove_isotopic_hydrogens(value)
    }
    fn with_remove_isotopic_hydrogens(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_isotopic_hydrogens(value);
        result
    }
    fn reassign_stereo(&self) -> bool {
        self.inner.reassign_stereo()
    }
    fn set_reassign_stereo(&mut self, value: bool) {
        self.inner.set_reassign_stereo(value)
    }
    fn with_reassign_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_reassign_stereo(value);
        result
    }
}
fn invocation(
    py: Python<'_>,
    params: Option<&TautomerParams>,
) -> (ck::TautomerParams, CallbackFailure) {
    match params {
        Some(params) => params.invocation(py),
        None => (ck::TautomerParams::default(), Arc::new(Mutex::new(None))),
    }
}
fn operation_result<T>(
    py: Python<'_>,
    result: Result<T, ck::OperationError>,
    failure: CallbackFailure,
) -> PyResult<T> {
    if let Some(error) = failure.lock().expect("callback error mutex").take() {
        return Err(error);
    }
    result.map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
}
pub(crate) fn enumerate(
    py: Python<'_>,
    molecule: &Molecule,
    params: Option<&TautomerParams>,
) -> PyResult<TautomerEnumeration> {
    let (params, failure) = invocation(py, params);
    operation_result(
        py,
        molecule.inner.enumerate_tautomers_with_params(&params),
        failure,
    )
    .map(|inner| TautomerEnumeration { inner })
}
pub(crate) fn canonical(
    py: Python<'_>,
    molecule: &Molecule,
    params: Option<&TautomerParams>,
) -> PyResult<Molecule> {
    let (params, failure) = invocation(py, params);
    operation_result(
        py,
        molecule.inner.canonical_tautomer_with_params(&params),
        failure,
    )
    .map(|inner| Molecule { inner })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TautomerEnumeration {
    inner: ck::TautomerEnumeration,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl TautomerEnumeration {
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn __len__(&self) -> usize {
        self.len()
    }
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status().into()
    }
    fn modified_atoms(&self) -> Vec<usize> {
        self.inner
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn modified_bonds(&self) -> Vec<usize> {
        self.inner
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    fn canonical_smiles(&self) -> Vec<String> {
        self.inner
            .canonical_smiles()
            .into_iter()
            .map(str::to_owned)
            .collect()
    }
    fn get(&self, index: usize) -> Option<Molecule> {
        self.inner
            .get(index)
            .cloned()
            .map(|inner| Molecule { inner })
    }
    fn __getitem__(&self, index: isize) -> PyResult<Molecule> {
        // RDKit✔️✔️:     if (pos < 0) {
        // RDKit✔️✔️:       pos += size();
        // RDKit✔️✔️:     }
        // RDKit✔️✔️:     if (pos < 0 || pos >= size()) {
        // RDKit✔️✔️:       PyErr_SetString(PyExc_IndexError, "index out of bounds");
        let index = if index < 0 {
            index + self.len() as isize
        } else {
            index
        };
        if index < 0 {
            return Err(PyIndexError::new_err("index out of bounds"));
        }
        self.get(index as usize)
            .ok_or_else(|| PyIndexError::new_err("index out of bounds"))
    }
    fn entries(&self) -> Vec<(String, Molecule)> {
        self.inner
            .entries()
            .map(|(key, inner)| {
                (
                    key.to_owned(),
                    Molecule {
                        inner: inner.clone(),
                    },
                )
            })
            .collect()
    }
    #[gen_stub(override_return_type(type_repr="typing.Iterator[Molecule]",imports=("typing")))]
    fn iter<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        let molecules = self
            .inner
            .iter()
            .cloned()
            .map(|inner| Py::new(py, Molecule { inner }))
            .collect::<PyResult<Vec<_>>>()?;
        PyList::new(py, molecules)?.call_method0("__iter__")
    }
    #[gen_stub(override_return_type(type_repr="typing.Iterator[Molecule]",imports=("typing")))]
    fn __iter__<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        self.iter(py)
    }
    fn canonical_tautomer(&self, py: Python<'_>) -> PyResult<Molecule> {
        self.inner
            .canonical_tautomer()
            .map(|inner| Molecule { inner })
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
    }
    fn canonical_tautomer_with_params(
        &self,
        py: Python<'_>,
        params: &TautomerParams,
    ) -> PyResult<Molecule> {
        let (params, failure) = params.invocation(py);
        operation_result(
            py,
            self.inner.canonical_tautomer_with_params(&params),
            failure,
        )
        .map(|inner| Molecule { inner })
    }
}
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn default_tautomer_score_terms() -> Vec<TautomerScoreTerm> {
    ck::default_tautomer_score_terms()
        .iter()
        .cloned()
        .map(|inner| TautomerScoreTerm { inner })
        .collect()
}
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn canonical_tautomer_from_molecules(
    py: Python<'_>,
    #[gen_stub(override_type(type_repr="typing.Iterable[Molecule]",imports=("typing")))]
    molecules: &Bound<'_, PyAny>,
) -> PyResult<Molecule> {
    canonical_tautomer_from_iterable(py, molecules, None)
}
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pyfunction]
fn canonical_tautomer_from_molecules_with_params(
    py: Python<'_>,
    #[gen_stub(override_type(type_repr="typing.Iterable[Molecule]",imports=("typing")))]
    molecules: &Bound<'_, PyAny>,
    params: &TautomerParams,
) -> PyResult<Molecule> {
    canonical_tautomer_from_iterable(py, molecules, Some(params))
}
fn canonical_tautomer_from_iterable(
    py: Python<'_>,
    molecules: &Bound<'_, PyAny>,
    params: Option<&TautomerParams>,
) -> PyResult<Molecule> {
    let values = molecules
        .try_iter()?
        .map(|value| -> PyResult<ck::Molecule> {
            let value = value?;
            let molecule = value.extract::<PyRef<'_, Molecule>>()?;
            Ok(molecule.inner.clone())
        })
        .collect::<PyResult<Vec<_>>>()?;
    let (params, failure) = invocation(py, params);
    operation_result(
        py,
        ck::canonical_tautomer_from_molecules_with_params(&values, &params),
        failure,
    )
    .map(|inner| Molecule { inner })
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "TautomerRunError",
        module.py().get_type::<TautomerRunError>(),
    )?;
    module.add(
        "TautomerCatalogError",
        module.py().get_type::<TautomerCatalogError>(),
    )?;
    module.add_class::<TautomerEnumerationStatus>()?;
    module.add_class::<TautomerScoreTerm>()?;
    module.add_class::<TautomerScoreParams>()?;
    module.add_class::<TautomerScore>()?;
    module.add_class::<TautomerParams>()?;
    module.add_class::<TautomerEnumeration>()?;
    module.add_class::<TautomerMoleculeView>()?;
    module.add_class::<TautomerProgress>()?;
    module.add_function(wrap_pyfunction!(default_tautomer_score_terms, module)?)?;
    module.add_function(wrap_pyfunction!(canonical_tautomer_from_molecules, module)?)?;
    module.add_function(wrap_pyfunction!(
        canonical_tautomer_from_molecules_with_params,
        module
    )?)?;
    Ok(())
}
