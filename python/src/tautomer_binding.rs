//! Thin language projections of canonical TAU configuration and validated results.
use crate::drawing_binding::Molecule;
use ::cosmolkit as ck;
use pyo3::exceptions::{PyIndexError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyList;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pymethods};
use std::sync::{Arc, Mutex};

pyo3::create_exception!(
    cosmolkit,
    TautomerRunError,
    PyValueError,
    "Tautomer enumeration, scoring or canonicalization failed."
);
pyo3::create_exception!(
    cosmolkit,
    TautomerCatalogError,
    PyValueError,
    "A tautomer transformation or scoring catalog could not be loaded or interpreted."
);

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
        E::Coordinates(..) => "Coordinates",
        E::SourceCacheCommit => "SourceCacheCommit",
        E::ScoreCacheTargetMismatch => "ScoreCacheTargetMismatch",
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

/// Completion state of a tautomer-enumeration run.
///
/// Declared values: ``Completed``, ``MaxTautomersReached``, ``MaxTransformsReached``, ``Canceled``.
#[cosmolkit_macros::python_enum]
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

/// Named SMARTS pattern and integer score used by tautomer ranking.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct TautomerScoreTerm {
    inner: ck::TautomerScoreTerm,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScoreTerm {
    /// Construct a TautomerScoreTerm value from the supplied inputs.
    #[new]
    fn py_new(name: String, smarts: String, score: i32) -> Self {
        Self::new(name, smarts, score)
    }
    /// Construct a TautomerScoreTerm value from the supplied inputs.
    #[staticmethod]
    fn new(name: String, smarts: String, score: i32) -> Self {
        Self {
            inner: ck::TautomerScoreTerm::new(name, smarts, score),
        }
    }
    /// Stored name of this value.
    fn name(&self) -> &str {
        self.inner.name()
    }
    /// SMARTS text defining the query or scoring pattern.
    fn smarts(&self) -> &str {
        self.inner.smarts()
    }
    /// Integer weight contributed by a matching tautomer scoring term.
    fn score(&self) -> i32 {
        self.inner.score()
    }
    fn __eq__(&self, py: Python<'_>, other: &Bound<'_, PyAny>) -> PyResult<Py<PyAny>> {
        match other.extract::<PyRef<'_, Self>>() {
            Ok(other) => Ok((self.inner == other.inner)
                .into_pyobject(py)?
                .to_owned()
                .into_any()
                .unbind()),
            Err(_) => Ok(py.NotImplemented()),
        }
    }
}
/// Writable configuration for tautomer SMARTS scoring terms.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
#[derive(Clone)]
pub(crate) struct TautomerScoreParams {
    pub(crate) inner: ck::TautomerScoreParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScoreParams {
    /// Configure tautomer SMARTS scoring terms; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(terms=None))]
    fn new(terms: Option<Vec<PyRef<'_, TautomerScoreTerm>>>) -> Self {
        Self {
            inner: ck::TautomerScoreParams {
                terms: terms.map(|terms| terms.iter().map(|term| term.inner.clone()).collect()),
            },
        }
    }
    /// Named SMARTS scoring terms used to rank tautomers.
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
/// Tautomer ranking contributions for ring patterns, SMARTS terms and heteroatom hydrogens, plus their total.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct TautomerScore {
    pub(crate) inner: ck::TautomerScore,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerScore {
    /// Ring-pattern contribution to the tautomer score.
    fn ring(&self) -> i32 {
        self.inner.ring()
    }
    /// SMARTS substructure contribution to the tautomer score.
    fn substructure(&self) -> i32 {
        self.inner.substructure()
    }
    /// Heteroatom hydrogen contribution to the tautomer score.
    fn hetero_hydrogen(&self) -> i32 {
        self.inner.hetero_hydrogen()
    }
    /// Total number of processed items.
    fn total(&self) -> i32 {
        self.inner.total()
    }
}

/// Read-only molecule view valid during a tautomer callback. Use to_owned() for an independent Molecule that outlives the callback.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TautomerMoleculeView {
    inner: ck::TautomerMoleculeView<'static>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerMoleculeView {
    /// Return a detached snapshot of the stored properties.
    fn properties(&self) -> crate::canonical_property_values::MoleculeProperties {
        crate::canonical_property_values::MoleculeProperties {
            inner: self.inner.properties().clone(),
        }
    }
    /// Return an owned snapshot that remains usable independently of this borrowed callback view.
    fn to_owned(&self) -> Self {
        Self {
            inner: self.inner.to_owned(),
        }
    }
    /// Number of atoms in the graph or selected structure.
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    /// Number of bonds in the graph.
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    /// Return atom rows in stored graph/hierarchy order.
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
    /// Return bond rows in graph order.
    fn bonds(&self) -> Vec<crate::canonical_atom_bond::Bond> {
        self.inner
            .bonds()
            .iter()
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
            .collect()
    }
    /// Return the read-only atom at the zero-based graph index, or None if the index is out of range.
    fn atom(&self, index: usize) -> Option<crate::canonical_atom_bond::Atom> {
        let inner = self.inner.atom(ck::AtomId::new(index))?.clone();
        Some(crate::canonical_atom_bond::Atom {
            inner,
            degree: self.inner.atom_degree(ck::AtomId::new(index))?,
            metadata: self.inner.atom_metadata().map(|rows| rows[index].clone()),
        })
    }
    /// Return the read-only bond at the zero-based graph index, or None if the index is out of range.
    fn bond(&self, index: usize) -> Option<crate::canonical_atom_bond::Bond> {
        self.inner
            .bond(ck::BondId::new(index))
            .cloned()
            .map(|inner| crate::canonical_atom_bond::Bond { inner })
    }
    /// Degree from the validated detached adjacency, independent of valence errors.
    fn atom_degree(&self, index: usize) -> Option<usize> {
        self.inner.atom_degree(ck::AtomId::new(index))
    }
    /// Canonical metadata query through its unique foundational owner.
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
    /// Return a SMILES string using the selected writer settings; this does not modify the molecule.
    fn to_smiles(&self, py: Python<'_>) -> PyResult<String> {
        self.inner
            .to_smiles()
            .map_err(|error| run_pyerr(py, &error))
            .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
    }
    /// Return ring, SMARTS-pattern and heteroatom-hydrogen score contributions for this tautomer.
    fn tautomer_score(&mut self, py: Python<'_>) -> PyResult<TautomerScore> {
        self.inner
            .tautomer_score()
            .map(|inner| TautomerScore { inner })
            .map_err(|error| run_pyerr(py, &error))
    }
}
/// Read-only tautomer callback progress view. Use to_owned() to retain a snapshot after the callback.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TautomerProgress {
    inner: ck::TautomerProgress<'static>,
    score_entries: Vec<(ck::PropertyText, Py<TautomerMoleculeView>)>,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerProgress {
    /// Return an owned snapshot that remains usable independently of this borrowed callback view.
    fn to_owned(&self, py: Python<'_>) -> PyResult<Self> {
        let score_entries = self
            .score_entries
            .iter()
            .map(|(key, value)| {
                let snapshot = value.bind(py).try_borrow()?.inner.to_owned();
                Ok((
                    key.clone(),
                    Py::new(py, TautomerMoleculeView { inner: snapshot })?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            inner: self.inner.to_owned(),
            score_entries,
        })
    }
    /// Return the number of stored entries.
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn __len__(&self) -> usize {
        self.len()
    }
    /// Return whether there are no stored entries.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    /// Return the current enumeration termination/progress status.
    fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status().into()
    }
    /// Number of tautomer transformations attempted so far.
    fn num_transforms(&self) -> u32 {
        self.inner.num_transforms()
    }
    /// Indices of atoms participating in enumerated tautomer transformations.
    fn modified_atoms(&self) -> Vec<usize> {
        self.inner
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    /// Indices of bonds participating in enumerated tautomer transformations.
    fn modified_bonds(&self) -> Vec<usize> {
        self.inner
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    /// Return retained entries in their stored order.
    fn entries(&mut self, py: Python<'_>) -> PyResult<Vec<(String, Py<TautomerMoleculeView>)>> {
        self.score_entries
            .iter()
            .map(|(key, value)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    value.clone_ref(py),
                ))
            })
            .collect()
    }
}
impl TautomerProgress {
    fn from_live(py: Python<'_>, progress: &mut ck::TautomerProgress<'_>) -> PyResult<Self> {
        let inner = progress.to_owned();
        let score_entries = progress
            .entries()
            .map(|(key, value)| {
                Ok((
                    key.clone(),
                    Py::new(
                        py,
                        TautomerMoleculeView {
                            inner: value.to_owned(),
                        },
                    )?,
                ))
            })
            .collect::<PyResult<Vec<_>>>()?;
        Ok(Self {
            inner,
            score_entries,
        })
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
        mut source: ck::TautomerMoleculeView<'_>,
        mut progress: ck::TautomerProgress<'_>,
    ) -> Result<bool, ck::TautomerRunError> {
        // BEGIN RDKIT CPP FUNCTION PyTautomerEnumeratorCallback::operator()
        // RDKit❗❌:   bool operator()(
        // RDKit❗❌:       const ROMol &mol,
        // RDKit❗❌:       const MolStandardize::TautomerEnumeratorResult &res) override {
        // RDKit❗❌:     PyTautomerEnumeratorResult pyRes(res);
        // RDKit❗❌:     return getCallbackOverride()(boost::ref(mol), boost::ref(pyRes));
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION PyTautomerEnumeratorCallback::operator()
        // Python may retain invocation snapshots. Only completed score-cache
        // writes are transported back into the actual restricted loans.
        // Source result-map copies share ROMol pointers; detached snapshots
        // alone would lose the .6 scoreRings state change.
        Python::attach(|py| {
            let source_value = Py::new(
                py,
                TautomerMoleculeView {
                    inner: source.to_owned(),
                },
            )
            .map_err(|error| retain_failure(&self.failure, error))?;
            let progress_value = TautomerProgress::from_live(py, &mut progress)
                .and_then(|value| Py::new(py, value))
                .map_err(|error| retain_failure(&self.failure, error))?;
            let result = self
                .callable
                .bind(py)
                .call1((source_value.clone_ref(py), progress_value.clone_ref(py)))
                .and_then(|value| value.is_truthy())
                .map_err(|error| retain_failure(&self.failure, error));
            {
                let snapshot = source_value
                    .bind(py)
                    .try_borrow()
                    .map_err(|error| retain_failure(&self.failure, error.into()))?;
                source.retain_score_cache_from(&snapshot.inner)?;
            }
            let snapshot = progress_value
                .bind(py)
                .try_borrow()
                .map_err(|error| retain_failure(&self.failure, error.into()))?;
            let mut actual = progress.entries();
            if actual.len() != snapshot.score_entries.len() {
                return Err(ck::TautomerRunError::ScoreCacheTargetMismatch);
            }
            for ((key, mut target), (saved_key, value)) in
                actual.by_ref().zip(&snapshot.score_entries)
            {
                if key != saved_key {
                    return Err(ck::TautomerRunError::ScoreCacheTargetMismatch);
                }
                let saved = value
                    .bind(py)
                    .try_borrow()
                    .map_err(|error| retain_failure(&self.failure, error.into()))?;
                target.retain_score_cache_from(&saved.inner)?;
            }
            result
        })
    }
}
impl ck::TautomerScorer for PyScorer {
    fn score(
        &self,
        mut molecule: ck::TautomerMoleculeView<'_>,
    ) -> Result<i32, ck::TautomerRunError> {
        // BEGIN RDKIT CPP FUNCTION pyobjFunctor::operator()
        // RDKit❗❌:   int operator()(const ROMol &m) {
        // RDKit❗❌:     return python::extract<int>(dp_obj(boost::ref(m)));
        // RDKit❗❌:   }
        // END RDKIT CPP FUNCTION pyobjFunctor::operator()
        // Preserve one scorer call and signed extraction. The existing graph
        // snapshot cost is explicit; actual cache publication precedes a later
        // Python/extraction error and never recomputes the score.
        Python::attach(|py| {
            let value = Py::new(
                py,
                TautomerMoleculeView {
                    inner: molecule.to_owned(),
                },
            )
            .map_err(|error| retain_failure(&self.failure, error))?;
            let result = self
                .callable
                .bind(py)
                .call1((value.clone_ref(py),))
                .and_then(|result| result.extract::<i32>())
                .map_err(|error| retain_failure(&self.failure, error));
            let snapshot = value
                .bind(py)
                .try_borrow()
                .map_err(|error| retain_failure(&self.failure, error.into()))?;
            molecule.retain_score_cache_from(&snapshot.inner)?;
            result
        })
    }
}

/// Writable configuration for tautomer enumeration and canonical selection.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", skip_from_py_object, dict, weakref)]
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
#[cosmolkit_macros::python_configuration(existing_setters)]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl TautomerParams {
    /// Configure tautomer enumeration and canonical selection; omitted fields use the defaults shown in the signature.
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
    /// Return tautomer enumeration parameters using the v1 transform rule set.
    #[staticmethod]
    fn v1(py: Python<'_>) -> PyResult<Self> {
        ck::TautomerParams::v1()
            .map(Self::from_inner)
            .map_err(|e| catalog_pyerr(py, e))
    }
    /// Construct tautomer enumeration parameters from explicitly supplied transformation rule text.
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
    /// Load tautomer transformation rules from the supplied filesystem path.
    #[staticmethod]
    fn from_transform_file(py: Python<'_>, path: String) -> PyResult<Self> {
        ck::TautomerParams::from_transform_file(path)
            .map(Self::from_inner)
            .map_err(|e| catalog_pyerr(py, e))
    }
    /// Number of transformation rules loaded in the enumerator.
    fn transform_count(&self) -> usize {
        self.inner.transform_count()
    }
    /// Optional tautomer progress callback; return False to stop enumeration and retain the corresponding status.
    #[getter]
    fn callback(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.callback.as_ref().map(|x| x.clone_ref(py))
    }
    /// Replace the tautomer progress callback; callback failures propagate rather than being ignored.
    fn set_callback(&mut self, py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<()> {
        self.callback = callable(py, value)?;
        Ok(())
    }
    /// Optional callback computing a tautomer ranking score from a read-only molecule view.
    #[getter]
    fn scorer(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.scorer.as_ref().map(|x| x.clone_ref(py))
    }
    /// Replace the tautomer ranking callback.
    fn set_scorer(&mut self, py: Python<'_>, value: Option<Py<PyAny>>) -> PyResult<()> {
        self.scorer = callable(py, value)?;
        Ok(())
    }
    /// SMARTS scoring terms used when no custom scorer is supplied.
    #[getter]
    fn score_params(&self) -> TautomerScoreParams {
        TautomerScoreParams {
            inner: self.inner.score_params.clone(),
        }
    }
    #[setter]
    fn set_score_params(&mut self, value: Option<&TautomerScoreParams>) {
        self.inner.score_params = value.map(|params| params.inner.clone()).unwrap_or_default()
    }
    /// Maximum number of distinct tautomers retained during enumeration.
    #[getter]
    fn max_tautomers(&self) -> u32 {
        self.inner.max_tautomers()
    }
    /// Update this configuration setting. Maximum number of distinct tautomers retained during enumeration.
    fn set_max_tautomers(&mut self, value: u32) {
        self.inner.set_max_tautomers(value)
    }
    /// Return a new configuration with this setting. Maximum number of distinct tautomers retained during enumeration.
    fn with_max_tautomers(&self, py: Python<'_>, value: u32) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_max_tautomers(value);
        result
    }
    /// Maximum number of tautomer transformations attempted.
    #[getter]
    fn max_transforms(&self) -> u32 {
        self.inner.max_transforms()
    }
    /// Update this configuration setting. Maximum number of tautomer transformations attempted.
    fn set_max_transforms(&mut self, value: u32) {
        self.inner.set_max_transforms(value)
    }
    /// Return a new configuration with this setting. Maximum number of tautomer transformations attempted.
    fn with_max_transforms(&self, py: Python<'_>, value: u32) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_max_transforms(value);
        result
    }
    /// Whether stereo labels at potentially modified tetrahedral centers are removed.
    #[getter]
    fn remove_sp3_stereo(&self) -> bool {
        self.inner.remove_sp3_stereo()
    }
    /// Update this configuration setting. Whether stereo labels at potentially modified tetrahedral centers are removed.
    fn set_remove_sp3_stereo(&mut self, value: bool) {
        self.inner.set_remove_sp3_stereo(value)
    }
    /// Return a new configuration with this setting. Whether stereo labels at potentially modified tetrahedral centers are removed.
    fn with_remove_sp3_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_sp3_stereo(value);
        result
    }
    /// Whether stereo labels on potentially modified bonds are removed.
    #[getter]
    fn remove_bond_stereo(&self) -> bool {
        self.inner.remove_bond_stereo()
    }
    /// Update this configuration setting. Whether stereo labels on potentially modified bonds are removed.
    fn set_remove_bond_stereo(&mut self, value: bool) {
        self.inner.set_remove_bond_stereo(value)
    }
    /// Return a new configuration with this setting. Whether stereo labels on potentially modified bonds are removed.
    fn with_remove_bond_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_bond_stereo(value);
        result
    }
    /// Whether isotope labels on hydrogens involved in tautomer transforms are removed.
    #[getter]
    fn remove_isotopic_hydrogens(&self) -> bool {
        self.inner.remove_isotopic_hydrogens()
    }
    /// Update this configuration setting. Whether isotope labels on hydrogens involved in tautomer transforms are removed.
    fn set_remove_isotopic_hydrogens(&mut self, value: bool) {
        self.inner.set_remove_isotopic_hydrogens(value)
    }
    /// Return a new configuration with this setting. Whether isotope labels on hydrogens involved in tautomer transforms are removed.
    fn with_remove_isotopic_hydrogens(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_remove_isotopic_hydrogens(value);
        result
    }
    /// Whether stereochemistry is reassigned after tautomer transformations.
    #[getter]
    fn reassign_stereo(&self) -> bool {
        self.inner.reassign_stereo()
    }
    /// Update this configuration setting. Whether stereochemistry is reassigned after tautomer transformations.
    fn set_reassign_stereo(&mut self, value: bool) {
        self.inner.set_reassign_stereo(value)
    }
    /// Return a new configuration with this setting. Whether stereochemistry is reassigned after tautomer transformations.
    fn with_reassign_stereo(&self, py: Python<'_>, value: bool) -> Self {
        let mut result = self.cloned(py);
        result.inner.set_reassign_stereo(value);
        result
    }

    #[setter(max_tautomers)]
    fn assign_max_tautomers(&mut self, value: u32) -> PyResult<()> {
        self.set_max_tautomers(value);
        Ok(())
    }
    #[setter(max_transforms)]
    fn assign_max_transforms(&mut self, value: u32) -> PyResult<()> {
        self.set_max_transforms(value);
        Ok(())
    }
    #[setter(remove_sp3_stereo)]
    fn assign_remove_sp3_stereo(&mut self, value: bool) -> PyResult<()> {
        self.set_remove_sp3_stereo(value);
        Ok(())
    }
    #[setter(remove_bond_stereo)]
    fn assign_remove_bond_stereo(&mut self, value: bool) -> PyResult<()> {
        self.set_remove_bond_stereo(value);
        Ok(())
    }
    #[setter(remove_isotopic_hydrogens)]
    fn assign_remove_isotopic_hydrogens(&mut self, value: bool) -> PyResult<()> {
        self.set_remove_isotopic_hydrogens(value);
        Ok(())
    }
    #[setter(reassign_stereo)]
    fn assign_reassign_stereo(&mut self, value: bool) -> PyResult<()> {
        self.set_reassign_stereo(value);
        Ok(())
    }
    #[setter(callback)]
    fn assign_callback(slf: &Bound<'_, Self>, value: Option<Py<PyAny>>) -> PyResult<()> {
        slf.try_borrow_mut()?.set_callback(slf.py(), value)
    }
    #[setter(scorer)]
    fn assign_scorer(slf: &Bound<'_, Self>, value: Option<Py<PyAny>>) -> PyResult<()> {
        slf.try_borrow_mut()?.set_scorer(slf.py(), value)
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
/// Ordered tautomer enumeration result with termination status, modified graph indices and canonical selection methods.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct TautomerEnumeration {
    inner: ck::TautomerEnumeration,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[cfg_attr(not(feature = "stubgen"), pyo3_stub_gen_derive::remove_gen_stub)]
#[pymethods]
impl TautomerEnumeration {
    /// Return the number of stored entries.
    fn len(&self) -> usize {
        self.inner.len()
    }
    fn __len__(&self) -> usize {
        self.len()
    }
    /// Return whether there are no stored entries.
    fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
    /// Return why enumeration finished, including completion, limits or callback cancellation.
    fn status(&self) -> TautomerEnumerationStatus {
        self.inner.status().into()
    }
    /// Indices of atoms participating in enumerated tautomer transformations.
    fn modified_atoms(&self) -> Vec<usize> {
        self.inner
            .modified_atoms()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    /// Indices of bonds participating in enumerated tautomer transformations.
    fn modified_bonds(&self) -> Vec<usize> {
        self.inner
            .modified_bonds()
            .iter()
            .map(|id| id.index())
            .collect()
    }
    /// Canonical SMILES strings for retained tautomers in enumeration order.
    fn canonical_smiles(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        self.inner
            .canonical_smiles()
            .into_iter()
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .collect()
    }
    /// Return the retained tautomer at the supplied zero-based index.
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
    /// Return retained tautomer molecules in enumeration order.
    fn entries(&self, py: Python<'_>) -> PyResult<Vec<(String, Molecule)>> {
        self.inner
            .entries()
            .map(|(key, inner)| {
                Ok((
                    crate::canonical_sdf::decode_source_text(py, key)?,
                    Molecule {
                        inner: inner.clone(),
                    },
                ))
            })
            .collect()
    }
    /// Return an iterator over entries in stored order.
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
    /// Return the highest-ranked canonical tautomer using the configured scoring rules; leave the source unchanged.
    fn canonical_tautomer(&self, py: Python<'_>) -> PyResult<Molecule> {
        self.inner
            .canonical_tautomer()
            .map(|inner| Molecule { inner })
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, error))
    }
    /// Return the highest-ranked canonical tautomer using the configured scoring rules; leave the source unchanged. Uses the supplied configuration object.
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
/// Return the built-in named SMARTS tautomer scoring terms.
#[cfg_attr(feature = "stubgen", pyo3_stub_gen::derive::gen_stub_pyfunction)]
#[pyfunction]
fn default_tautomer_score_terms() -> Vec<TautomerScoreTerm> {
    ck::default_tautomer_score_terms()
        .iter()
        .cloned()
        .map(|inner| TautomerScoreTerm { inner })
        .collect()
}
/// Choose the canonical tautomer from the supplied molecules using the configured score terms.
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
/// Choose the canonical tautomer from the supplied molecules using the configured score terms. Uses the supplied configuration object.
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
    // Preserve input order and repeated handles while lending shared values.
    // The facade owns the isolated COW candidates used by scoring.
    let values = molecules
        .try_iter()?
        .map(|value| -> PyResult<Py<Molecule>> { Ok(value?.extract::<Py<Molecule>>()?) })
        .collect::<PyResult<Vec<_>>>()?;
    let (params, failure) = invocation(py, params);
    let host_failure = failure.clone();
    let with_candidate =
        |index: usize, action: &mut dyn FnMut(&ck::Molecule) -> Result<(), ck::OperationError>| {
            let value = values[index].bind(py).try_borrow().map_err(|error| {
                ck::OperationError::Tautomer(retain_failure(&host_failure, error.into()))
            })?;
            action(&value.inner)
        };
    operation_result(
        py,
        ck::canonical_tautomer_from_molecule_hosts_with_params(
            values.len(),
            with_candidate,
            &params,
        ),
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

#[cfg(all(test, feature = "python-embed-tests"))]
mod recovery_search06_python_cache_tests {
    use super::*;
    use std::ffi::CString;
    #[test]
    fn recovery_search06_python_iterable_repeated_object_preserves_original() {
        Python::initialize();
        Python::attach(|py| {
            let source = ck::Molecule::from_smiles_with_params(
                "c1ccccc1",
                &ck::SmilesParseParams {
                    sanitize: false,
                    remove_hs: false,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(matches!(
                source.num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
            let value = Py::new(py, Molecule { inner: source }).unwrap();
            let inputs = PyList::new(py, [value.clone_ref(py), value.clone_ref(py)]).unwrap();
            canonical_tautomer_from_iterable(py, inputs.as_any(), None).unwrap();
            assert!(matches!(
                value.bind(py).borrow().inner.num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
        });
    }
    #[test]
    fn recovery_search06_python_scorer_error_preserves_all_inputs() {
        Python::initialize();
        Python::attach(|py| {
            let code=CString::new("def score(m):\n    try:\n        m.tautomer_score()\n    except ValueError:\n        pass\n    raise ValueError('after score')\n").unwrap();
            let filename = CString::new("search06-local.rs-inline").unwrap();
            let name = CString::new("search06_local").unwrap();
            let module =
                PyModule::from_code(py, code.as_c_str(), filename.as_c_str(), name.as_c_str())
                    .unwrap();
            let callable = module.getattr("score").unwrap().unbind();
            let mut params = ck::TautomerParams::default();
            let failure: CallbackFailure = Default::default();
            params.set_scorer(Some(Arc::new(PyScorer {
                callable,
                failure: failure.clone(),
            })));
            let cold = || {
                ck::Molecule::from_smiles_with_params(
                    "c1ccccc1",
                    &ck::SmilesParseParams {
                        sanitize: false,
                        remove_hs: false,
                        ..Default::default()
                    },
                )
                .unwrap()
            };
            let mut inputs = [cold(), cold()];
            let result = ck::canonical_tautomer_from_molecules_with_params(&mut inputs, &params);
            assert!(matches!(
                result,
                Err(ck::OperationError::Tautomer(
                    ck::TautomerRunError::Callback(_)
                ))
            ));
            assert!(failure.lock().unwrap().is_some());
            assert!(matches!(
                inputs[0].num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
            assert!(matches!(
                inputs[1].num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
        });
    }
    #[test]
    fn recovery_search06_python_callback_error_preserves_source_and_retains_view() {
        Python::initialize();
        Python::attach(|py| {
            let code=CString::new("saved = []\ndef callback(source, progress):\n    try:\n        source.tautomer_score()\n    except ValueError:\n        pass\n    first = progress.entries()\n    again = progress.entries()\n    assert first[0][1] is again[0][1]\n    try:\n        first[0][1].tautomer_score()\n    except ValueError:\n        pass\n    saved.append(first[0][1])\n    raise ValueError('after source and entry scores')\n").unwrap();
            let filename = CString::new("search06-callback.rs-inline").unwrap();
            let name = CString::new("search06_callback").unwrap();
            let module =
                PyModule::from_code(py, code.as_c_str(), filename.as_c_str(), name.as_c_str())
                    .unwrap();
            let callable = module.getattr("callback").unwrap().unbind();
            let failure: CallbackFailure = Default::default();
            let mut params = ck::TautomerParams::default();
            params.set_callback(Some(Arc::new(PyCallback {
                callable,
                failure: failure.clone(),
            })));
            let mut source = ck::Molecule::from_smiles_with_params(
                "CC(=O)c1ccccc1",
                &ck::SmilesParseParams {
                    sanitize: false,
                    remove_hs: false,
                    ..Default::default()
                },
            )
            .unwrap();
            assert!(matches!(
                source.num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
            let result = source.enumerate_tautomers_with_params(&params);
            assert!(matches!(
                result,
                Err(ck::OperationError::Tautomer(
                    ck::TautomerRunError::Callback(_)
                ))
            ));
            assert!(failure.lock().unwrap().is_some());
            assert!(matches!(
                source.num_rings(),
                Err(ck::DescriptorReadError::MissingInitializedRings)
            ));
            assert_eq!(module.getattr("saved").unwrap().len().unwrap(), 1);
        });
    }
}
