//! Canonical detached query and search projections; chemistry lives behind cosmolkit.
use ::cosmolkit as ck;
use pyo3::exceptions::{PyRuntimeError, PyTypeError, PyValueError};
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

/// A query atom is not an Element-only Atom: wildcard/OR queries can have Z=0.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct QueryAtom {
    inner: ck::QueryAtom,
    degree: usize,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl QueryAtom {
    fn id(&self) -> usize {
        self.inner.id().index()
    }
    fn atomic_number(&self) -> u8 {
        self.inner.atomic_number()
    }
    fn formal_charge(&self) -> i8 {
        self.inner.formal_charge()
    }
    fn explicit_hydrogens(&self) -> u8 {
        self.inner.explicit_hydrogens()
    }
    fn isotope(&self) -> Option<u16> {
        self.inner.isotope()
    }
    fn atom_map(&self) -> Option<u32> {
        self.inner.atom_map()
    }
    fn is_aromatic(&self) -> bool {
        self.inner.is_aromatic()
    }
    fn no_implicit(&self) -> bool {
        self.inner.no_implicit()
    }
    fn radical_electrons(&self) -> u8 {
        self.inner.radical_electrons()
    }
    fn degree(&self) -> usize {
        self.degree
    }
    fn __repr__(&self) -> String {
        format!(
            "QueryAtom(id={}, atomic_number={})",
            self.id(),
            self.atomic_number()
        )
    }
}

pyo3::create_exception!(cosmolkit, SmartsParseError, PyValueError);
pyo3::create_exception!(cosmolkit, SmartsWriteError, PyValueError);
pyo3::create_exception!(cosmolkit, SubstructMatchError, PyValueError);
pyo3::create_exception!(cosmolkit, QueryCompileError, PyValueError);
pyo3::create_exception!(cosmolkit, MatchError, PyValueError);

pub(crate) fn parse_pyerr(py: Python<'_>, source: ck::SmartsParseError) -> PyErr {
    use ck::SmartsParseError as E;
    let kind = match &source {
        E::StereoGroup(_) => "StereoGroup",
        E::MissingRecursiveQueryGraph => "MissingRecursiveQueryGraph",
        E::CxLowering(_) => "CxLowering",
        E::QueryGraph(_) => "QueryGraph",
        E::ParserCarrier(_) => "ParserCarrier",
        E::AtomProperty(_) => "AtomProperty",
        E::BondProperty(_) => "BondProperty",
        E::MoleculeProperty(_) => "MoleculeProperty",
        E::UnclosedBracket(_) => "UnclosedBracket",
        E::UnexpectedCharacter { .. } => "UnexpectedCharacter",
        E::UnexpectedEnd(_) => "UnexpectedEnd",
        E::InvalidAtomPrimitive { .. } => "InvalidAtomPrimitive",
        E::UnclosedParenthesis(_) => "UnclosedParenthesis",
        E::UnbalancedRingClosure(_) => "UnbalancedRingClosure",
        E::SelfRingClosure { .. } => "SelfRingClosure",
        E::DuplicateRingBond { .. } => "DuplicateRingBond",
        E::CxSmiles(_) => "CxSmiles",
        E::Parse(_) => "Parse",
        E::UnsupportedFeature(_) => "UnsupportedFeature",
        E::TemplateAttachmentRemap { .. } => "TemplateAttachmentRemap",
    };
    let error = crate::canonical_values::annotate(
        py,
        SmartsParseError::new_err(source.to_string()),
        "search",
        kind,
        &source,
    );
    let attributes = || -> PyResult<()> {
        let value = error.value(py);
        match source {
            E::UnclosedBracket(position) | E::UnclosedParenthesis(position) => {
                value.setattr("position", position)?
            }
            E::UnexpectedCharacter {
                position,
                character,
                context,
            } => {
                value.setattr("position", position)?;
                value.setattr("character", char::from(character).to_string())?;
                value.setattr("context", context)?;
            }
            E::InvalidAtomPrimitive { position, detail } => {
                value.setattr("position", position)?;
                value.setattr("detail", detail)?;
            }
            E::UnbalancedRingClosure(number) => value.setattr("ring", number)?,
            E::SelfRingClosure { ring, atom } => {
                value.setattr("ring", ring)?;
                value.setattr("atom", atom)?;
            }
            E::DuplicateRingBond {
                ring,
                begin_atom,
                end_atom,
            } => {
                value.setattr("ring", ring)?;
                value.setattr("begin_atom", begin_atom)?;
                value.setattr("end_atom", end_atom)?;
            }
            E::UnsupportedFeature(feature) => value.setattr("feature", feature)?,
            E::TemplateAttachmentRemap { carrier, .. } => value.setattr("carrier", carrier)?,
            _ => (),
        }
        Ok(())
    };
    match attributes() {
        Ok(()) => error,
        Err(e) => e,
    }
}

pub(crate) fn write_pyerr(py: Python<'_>, source: ck::SmartsWriteError) -> PyErr {
    use ck::SmartsWriteError as E;
    let kind = match &source {
        E::StereoGroup(_) => "StereoGroup",
        E::Traversal(_) => "Traversal",
        E::CxCoordinates(_) => "CxCoordinates",
        E::CxRingInfo(_) => "CxRingInfo",
        E::CxRingAtomOrderIndex { .. } => "CxRingAtomOrderIndex",
        E::CxRingStereoReferenceMissing { .. } => "CxRingStereoReferenceMissing",
        E::CxWedge(_) => "CxWedge",
        E::CxBondConfigAtropMissingCarriers { .. } => "CxBondConfigAtropMissingCarriers",
        E::CxAtomPropertyOutput(_) => "CxAtomPropertyOutput",
        E::CxMissingConformer => "CxMissingConformer",
        E::CxCoordinateSource(_) => "CxCoordinateSource",
        E::CxCoordinateOutput(_) => "CxCoordinateOutput",
        E::CxAtomPropertyWrite { .. } => "CxAtomPropertyWrite",
        E::CxCoordinateStorage { .. } => "CxCoordinateStorage",
        E::CxSourceBondOutOfRange { .. } => "CxSourceBondOutOfRange",
        E::CxMoleculePropertyUInt { .. } => "CxMoleculePropertyUInt",
        E::CxSgroupVectorCast { .. } => "CxSgroupVectorCast",
        E::CxSgroupPropertyWrite { .. } => "CxSgroupPropertyWrite",
        E::CxStereoGroup(_) => "CxStereoGroup",
        E::CxSourceAtomOutOfRange { .. } => "CxSourceAtomOutOfRange",
        E::MoleculePropertyWrite(_) => "MoleculePropertyWrite",
        E::SourceAtomCount { .. } => "SourceAtomCount",
        E::AtomPropertyWrite { .. } => "AtomPropertyWrite",
        E::CanonicalTraversal(_) => "CanonicalTraversal",
        E::UnwritableBondQuery { .. } => "UnwritableBondQuery",
        E::SourceAtomToLeftIndex { .. } => "SourceAtomToLeftIndex",
        E::SourceBondBeginIndex { .. } => "SourceBondBeginIndex",
        E::AtomMapInt { .. } => "AtomMapInt",
        E::AtomTypeAtomicNumber { .. } => "AtomTypeAtomicNumber",
        E::ChargeMagnitudeOverflow { .. } => "ChargeMagnitudeOverflow",
        E::PropertyValue(_) => "PropertyValue",
        E::CxRequiredProperty { .. } => "CxRequiredProperty",
        E::CxPropertyList { .. } => "CxPropertyList",
        E::CxCoordinateSelectionArity { .. } => "CxCoordinateSelectionArity",
        E::CxMissingOutputOrder { .. } => "CxMissingOutputOrder",
        E::CxOutputOrderPropertyType { .. } => "CxOutputOrderPropertyType",
        E::CxAtomPropertyKind { .. } => "CxAtomPropertyKind",
        E::CxAtomPropertyUInt { .. } => "CxAtomPropertyUInt",
        E::CxCoordinateSelection { .. } => "CxCoordinateSelection",
        E::CxComposition(_) => "CxComposition",
        E::CxOutputOrder { .. } => "CxOutputOrder",
        E::CxRowCount { .. } => "CxRowCount",
        E::CxBondPropertyUInt { .. } => "CxBondPropertyUInt",
        E::CxAtomPropertyInt { .. } => "CxAtomPropertyInt",
        E::CxSgroupPropertyUInt { .. } => "CxSgroupPropertyUInt",
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::Property(_) => "Property",
        E::Valence(_) => "Valence",
        E::InvalidGraph(_) => "InvalidGraph",
        E::QueryGraphTraversalUnsupported { .. } => "QueryGraphTraversalUnsupported",
        E::OrAboveAndBelowAnd => "OrAboveAndBelowAnd",
        E::UnknownCombination { .. } => "UnknownCombination",
        E::MissingRecursiveQueryMolecule => "MissingRecursiveQueryMolecule",
        E::SourceBondDirection { .. } => "SourceBondDirection",
        E::UnsupportedBondQuery { .. } => "UnsupportedBondQuery",
        E::UnsupportedAtomQuery { .. } => "UnsupportedAtomQuery",
        E::CompositeChildCount { .. } => "CompositeChildCount",
        E::XorComposite => "XorComposite",
        E::QueryGraphCxExtensionsUnsupported { .. } => "QueryGraphCxExtensionsUnsupported",
        E::RootedAtomOutOfRange { .. } => "RootedAtomOutOfRange",
        E::EmptyAtomSelection => "EmptyAtomSelection",
        E::EmptyBondSelection => "EmptyBondSelection",
        E::FragmentAtomOutOfRange { .. } => "FragmentAtomOutOfRange",
        E::FragmentBondOutOfRange { .. } => "FragmentBondOutOfRange",
        E::BondAtomNotEndpoint { .. } => "BondAtomNotEndpoint",
    };
    crate::canonical_values::annotate(
        py,
        SmartsWriteError::new_err(source.to_string()),
        "search",
        kind,
        &source,
    )
}

pub(crate) fn substruct_pyerr(py: Python<'_>, source: ck::SubstructMatchError) -> PyErr {
    let kind = match &source {
        ck::SubstructMatchError::Unsupported { .. } => "Unsupported",
        ck::SubstructMatchError::FinalCheckMappingLength { .. } => "FinalCheckMappingLength",
        ck::SubstructMatchError::FinalCheckMappingIndex { .. } => "FinalCheckMappingIndex",
        ck::SubstructMatchError::FinalCheckInvariant { .. } => "FinalCheckInvariant",
        ck::SubstructMatchError::FinalCheckBondEndpoint { .. } => "FinalCheckBondEndpoint",
        ck::SubstructMatchError::FinalCheckMissingBond { .. } => "FinalCheckMissingBond",
        ck::SubstructMatchError::StereoOrder(_) => "StereoOrder",
        ck::SubstructMatchError::PeriodicTable(_) => "PeriodicTable",
        ck::SubstructMatchError::PropertyString(_) => "PropertyString",
        ck::SubstructMatchError::PropertyInteger { .. } => "PropertyInteger",
        ck::SubstructMatchError::QueryContext(_) => "QueryContext",
    };
    let error = crate::canonical_values::annotate(
        py,
        SubstructMatchError::new_err(source.to_string()),
        "search",
        kind,
        &source,
    );
    if let ck::SubstructMatchError::Unsupported {
        branch,
        rdkit_function,
    } = source
    {
        if let Err(e) = error
            .value(py)
            .setattr("branch", branch)
            .and_then(|()| error.value(py).setattr("rdkit_function", rdkit_function))
        {
            return e;
        }
    }
    error
}

pub(crate) fn match_pyerr(py: Python<'_>, source: ck::MatchError) -> PyErr {
    use ck::MatchError as E;
    let kind = match &source {
        E::InvalidQuery(_) => "InvalidQuery",
        E::InvalidTarget(_) => "InvalidTarget",
        E::UnsupportedAtomPredicate(_) => "UnsupportedAtomPredicate",
        E::UnsupportedBondPredicate(_) => "UnsupportedBondPredicate",
        E::Substruct(_) => "Substruct",
    };
    crate::canonical_values::annotate(
        py,
        MatchError::new_err(source.to_string()),
        "search",
        kind,
        &source,
    )
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct QueryGraph {
    pub(crate) inner: ck::QueryGraph,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl QueryGraph {
    #[staticmethod]
    fn from_smarts(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::search::from_smarts(text)
            .map(|inner| Self { inner })
            .map_err(|error| parse_pyerr(py, error))
    }

    #[staticmethod]
    fn from_smarts_with_params(
        py: Python<'_>,
        text: &str,
        params: &SmartsParseParams,
    ) -> PyResult<Self> {
        ck::search::from_smarts_with_params(text, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|error| parse_pyerr(py, error))
    }

    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    fn name(&self, py: Python<'_>) -> PyResult<Option<String>> {
        self.inner
            .name()
            .map_err(|error| crate::canonical_atom_bond::property_pyerr(py, error))?
            .map(|text| crate::canonical_sdf::decode_source_text(py, text))
            .transpose()
    }
    fn __len__(&self) -> usize {
        self.inner.num_atoms()
    }
    fn __repr__(&self) -> String {
        format!(
            "QueryGraph(atoms={}, bonds={})",
            self.inner.num_atoms(),
            self.inner.num_bonds()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct CompiledQuery {
    pub(crate) inner: ck::CompiledQuery,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl CompiledQuery {
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    fn atom_order(&self) -> Vec<usize> {
        self.inner.atom_order().to_vec()
    }
    fn query(&self) -> QueryGraph {
        QueryGraph {
            inner: self.inner.query().clone(),
        }
    }
    fn __repr__(&self) -> String {
        format!(
            "CompiledQuery(atoms={}, bonds={})",
            self.inner.num_atoms(),
            self.inner.num_bonds()
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct MatchResult {
    pub(crate) inner: ck::MatchResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl MatchResult {
    fn atom_mapping(&self) -> Vec<usize> {
        self.inner.atom_mapping.clone()
    }
    fn bond_mapping(&self) -> Vec<usize> {
        self.inner.bond_mapping.clone()
    }
    fn atom_pairs(&self) -> Vec<(usize, usize)> {
        // RDKit✔️✔️: std::for_each(matches.begin(), matches.end(), [res, &matches](const auto &pair) {
        // RDKit✔️✔️:   PyObject *pyPair = PyTuple_New(2);
        // RDKit✔️✔️:   PyTuple_SetItem(pyPair, 0, PyInt_FromLong(pair.first));
        // RDKit✔️✔️:   PyTuple_SetItem(pyPair, 1, PyInt_FromLong(pair.second));
        // RDKit✔️✔️:   PyTuple_SetItem(res, &pair - &matches.front(), pyPair);
        // RDKit✔️✔️: });
        // Dense query-indexed mappings already carry source order. One linear
        // pair projection allocates one pair per row, like the source wrapper.
        self.inner
            .atom_mapping
            .iter()
            .copied()
            .enumerate()
            .collect()
    }
    fn __repr__(&self) -> String {
        format!(
            "MatchResult(atom_mapping={:?}, bond_mapping={:?})",
            self.inner.atom_mapping, self.inner.bond_mapping
        )
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SmartsParseParams {
    pub(crate) inner: ck::SmartsParseParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmartsParseParams {
    #[new]
    #[pyo3(signature = (*, allow_cxsmiles=true, strict_cxsmiles=true, parse_name=true, merge_hs=false, skip_cleanup=false, debug_parse=false, replacements=None))]
    fn new(
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
        parse_name: bool,
        merge_hs: bool,
        skip_cleanup: bool,
        debug_parse: bool,
        replacements: Option<BTreeMap<String, String>>,
    ) -> Self {
        Self {
            inner: ck::SmartsParseParams {
                allow_cxsmiles,
                strict_cxsmiles,
                parse_name,
                merge_hs,
                skip_cleanup,
                debug_parse,
                replacements: replacements
                    .unwrap_or_default()
                    .into_iter()
                    .map(|(key, value)| (key.into(), value.into()))
                    .collect(),
            },
        }
    }
    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles
    }
    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles
    }
    #[getter]
    fn parse_name(&self) -> bool {
        self.inner.parse_name
    }
    #[getter]
    fn merge_hs(&self) -> bool {
        self.inner.merge_hs
    }
    #[getter]
    fn skip_cleanup(&self) -> bool {
        self.inner.skip_cleanup
    }
    #[getter]
    fn debug_parse(&self) -> bool {
        self.inner.debug_parse
    }
    #[getter]
    fn replacements(&self) -> PyResult<BTreeMap<String, String>> {
        self.inner
            .replacements
            .iter()
            .map(|(key, value)| {
                let key = std::str::from_utf8(key.as_bytes())
                    .map_err(|error| PyValueError::new_err(error.to_string()))?;
                let value = std::str::from_utf8(value.as_bytes())
                    .map_err(|error| PyValueError::new_err(error.to_string()))?;
                Ok((key.to_owned(), value.to_owned()))
            })
            .collect()
    }
    fn __repr__(&self) -> String {
        format!("SmartsParseParams({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SmartsWriteParams {
    pub(crate) inner: ck::SmartsWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmartsWriteParams {
    #[new]
    #[pyo3(signature = (*, include_atom_maps=true, isomeric_smiles=true, include_dative_bonds=true, rooted_at_atom=None))]
    fn new(
        include_atom_maps: bool,
        isomeric_smiles: bool,
        include_dative_bonds: bool,
        rooted_at_atom: Option<usize>,
    ) -> Self {
        Self {
            inner: ck::SmartsWriteParams {
                include_atom_maps,
                isomeric_smiles,
                include_dative_bonds,
                rooted_at_atom,
            },
        }
    }
    #[getter]
    fn include_atom_maps(&self) -> bool {
        self.inner.include_atom_maps
    }
    #[getter]
    fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles
    }
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds
    }
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom
    }
    fn __repr__(&self) -> String {
        format!("SmartsWriteParams({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct SubstructMatchParams {
    pub(crate) inner: ck::SubstructMatchParams,
    final_match: Option<Py<PyAny>>,
    atom_match: Option<Py<PyAny>>,
    bond_match: Option<Py<PyAny>>,
}

// Errors are per invocation, not attached to reusable params or thread-local
// state. Never hold this mutex while executing arbitrary/reentrant Python.
struct CallbackErrors(Mutex<Option<PyErr>>);
impl CallbackErrors {
    fn call(&self, run: impl FnOnce(Python<'_>) -> PyResult<bool>) -> bool {
        let mut slot = match self.0.lock() {
            Ok(slot) => slot,
            Err(poisoned) => {
                let mut slot = poisoned.into_inner();
                if slot.is_none() {
                    *slot = Some(PyRuntimeError::new_err(
                        "match callback error state poisoned",
                    ));
                }
                return false;
            }
        };
        if slot.is_some() {
            return false;
        }
        drop(slot);
        match Python::attach(run) {
            Ok(value) => value,
            Err(error) => {
                slot = self
                    .0
                    .lock()
                    .unwrap_or_else(|poisoned| poisoned.into_inner());
                if slot.is_none() {
                    *slot = Some(error);
                }
                false
            }
        }
    }
}

impl SubstructMatchParams {
    pub(crate) fn from_inner(inner: ck::SubstructMatchParams) -> Self {
        Self {
            inner,
            final_match: None,
            atom_match: None,
            bond_match: None,
        }
    }

    pub(crate) fn has_python_callbacks(&self) -> bool {
        self.final_match.is_some() || self.atom_match.is_some() || self.bond_match.is_some()
    }

    pub(crate) fn with_callbacks<T>(
        &self,
        py: Python<'_>,
        molecule: &ck::Molecule,
        run: impl FnOnce(&ck::SubstructMatchParams) -> Result<T, ck::SubstructMatchError>,
    ) -> PyResult<T> {
        let mut params = self.inner.clone();
        let errors = Arc::new(CallbackErrors(Mutex::new(None)));
        if let Some(callback) = &self.final_match {
            // RDKit✔️❌: bool operator()(const ROMol &m, std::span<const unsigned int> match) {
            // RDKit✔️❌:   // grab the GIL
            // RDKit✔️❌:   PyGILStateHolder h;
            // RDKit✔️❌:   // boost::python doesn't handle std::span, so we need to convert the span to
            // RDKit✔️❌:   // a vector before calling into python:
            // RDKit✔️❌:   std::vector<unsigned int> matchVec(match.begin(), match.end());
            // RDKit✔️❌:   return python::extract<bool>(dp_obj(boost::ref(m), boost::ref(matchVec)));
            // RDKit✔️❌: }
            // Same O(query atoms) mapping copy; unlike Boost's borrowed wrapper,
            // each callback owns a cheap COW molecule snapshot and error slot.
            let callback = callback.clone_ref(py);
            let molecule = molecule.clone();
            let errors = Arc::clone(&errors);
            params.extra_final_check = Some(Arc::new(move |_, ids| {
                errors.call(|py| {
                    let target = Py::new(
                        py,
                        crate::drawing_binding::Molecule::from_inner(molecule.clone()),
                    )?;
                    callback
                        .call1(py, (target, ids.to_vec()))?
                        .extract::<bool>(py)
                })
            }));
        }
        if let Some(callback) = &self.atom_match {
            // RDKit✔️❌: bool operator()(const T &a1, const T &a2) {
            // RDKit✔️❌:   // grab the GIL
            // RDKit✔️❌:   PyGILStateHolder h;
            // RDKit✔️❌:   return python::extract<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
            // RDKit✔️❌: }
            // Binding snapshots preserve wildcard query identities, but owning
            // their property values costs more than borrowed Boost wrappers.
            let callback = callback.clone_ref(py);
            let errors = Arc::clone(&errors);
            let molecule = molecule.clone();
            let metadata = molecule.atom_metadata(false);
            params.extra_atom_check = Some(Arc::new(move |graph, query, _target, atom| {
                errors.call(|py| {
                    let query = Py::new(
                        py,
                        QueryAtom {
                            inner: query.clone(),
                            degree: graph
                                .adjacency()
                                .get(query.id().index())
                                .map_or(0, Vec::len),
                        },
                    )?;
                    let target = Py::new(
                        py,
                        crate::canonical_atom_bond::Atom {
                            inner: atom.clone(),
                            degree: molecule
                                .topology()
                                .adjacency
                                .neighbors_of(atom.id().index())
                                .len(),
                            metadata: metadata
                                .as_ref()
                                .map(|rows| rows[atom.id().index()].clone())
                                .map_err(Clone::clone),
                        },
                    )?;
                    callback.call1(py, (query, target))?.extract::<bool>(py)
                })
            }));
        }
        if let Some(callback) = &self.bond_match {
            // RDKit✔️❌: bool operator()(const T &a1, const T &a2) {
            // RDKit✔️❌:   // grab the GIL
            // RDKit✔️❌:   PyGILStateHolder h;
            // RDKit✔️❌:   return python::extract<bool>(dp_obj(boost::ref(a1), boost::ref(a2)));
            // RDKit✔️❌: }
            // Bond specialization: detached owning values, not live storage.
            let callback = callback.clone_ref(py);
            let errors = Arc::clone(&errors);
            params.extra_bond_check = Some(Arc::new(move |query, target| {
                errors.call(|py| {
                    let query = Py::new(
                        py,
                        crate::canonical_atom_bond::Bond {
                            inner: query.clone(),
                        },
                    )?;
                    let target = Py::new(
                        py,
                        crate::canonical_atom_bond::Bond {
                            inner: target.clone(),
                        },
                    )?;
                    callback.call1(py, (query, target))?.extract::<bool>(py)
                })
            }));
        }
        let result = run(&params);
        let mut stored = errors
            .0
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        if let Some(error) = stored.take() {
            return Err(error);
        }
        result.map_err(|error| substruct_pyerr(py, error))
    }
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstructMatchParams {
    #[new]
    #[pyo3(signature = (*, max_matches=1000, uniquify=true, use_chirality=false, use_enhanced_stereo=false, specified_stereo_query_matches_unspecified=false, use_query_query_matches=false, recursion_possible=true, max_recursive_matches=1000, num_threads=1, aromatic_matches_conjugated=false, aromatic_matches_single_or_double=false, atom_properties=None, bond_properties=None, extra_atom_check_overrides_default_check=false, extra_bond_check_overrides_default_check=false, use_generic_matchers=false, final_match=None, atom_match=None, bond_match=None))]
    fn new(
        py: Python<'_>,
        max_matches: usize,
        uniquify: bool,
        use_chirality: bool,
        use_enhanced_stereo: bool,
        specified_stereo_query_matches_unspecified: bool,
        use_query_query_matches: bool,
        recursion_possible: bool,
        max_recursive_matches: usize,
        num_threads: i32,
        aromatic_matches_conjugated: bool,
        aromatic_matches_single_or_double: bool,
        atom_properties: Option<Vec<String>>,
        bond_properties: Option<Vec<String>>,
        extra_atom_check_overrides_default_check: bool,
        extra_bond_check_overrides_default_check: bool,
        use_generic_matchers: bool,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Callable[[Molecule, typing.Sequence[builtins.int]], builtins.bool]]", imports=("typing", "builtins")))]
        final_match: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Callable[[QueryAtom, Atom], builtins.bool]]", imports=("typing", "builtins")))]
        atom_match: Option<Py<PyAny>>,
        #[gen_stub(override_type(type_repr="typing.Optional[typing.Callable[[Bond, Bond], builtins.bool]]", imports=("typing", "builtins")))]
        bond_match: Option<Py<PyAny>>,
    ) -> PyResult<Self> {
        for (name, callback) in [
            ("final_match", &final_match),
            ("atom_match", &atom_match),
            ("bond_match", &bond_match),
        ] {
            if callback
                .as_ref()
                .is_some_and(|value| !value.bind(py).is_callable())
            {
                return Err(PyTypeError::new_err(format!(
                    "{name} must be callable or None"
                )));
            }
        }
        Ok(Self {
            inner: ck::SubstructMatchParams {
                max_matches,
                uniquify,
                use_chirality,
                use_enhanced_stereo,
                specified_stereo_query_matches_unspecified,
                use_query_query_matches,
                recursion_possible,
                max_recursive_matches,
                num_threads,
                aromatic_matches_conjugated,
                aromatic_matches_single_or_double,
                atom_properties: atom_properties.unwrap_or_default(),
                bond_properties: bond_properties.unwrap_or_default(),
                extra_atom_check_overrides_default_check,
                extra_bond_check_overrides_default_check,
                use_generic_matchers,
                ..Default::default()
            },
            final_match,
            atom_match,
            bond_match,
        })
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Optional[typing.Callable[[Molecule, typing.Sequence[builtins.int]], builtins.bool]]", imports=("typing", "builtins")))]
    fn final_match(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.final_match.as_ref().map(|value| value.clone_ref(py))
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Optional[typing.Callable[[QueryAtom, Atom], builtins.bool]]", imports=("typing", "builtins")))]
    fn atom_match(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.atom_match.as_ref().map(|value| value.clone_ref(py))
    }
    #[getter]
    #[gen_stub(override_return_type(type_repr="typing.Optional[typing.Callable[[Bond, Bond], builtins.bool]]", imports=("typing", "builtins")))]
    fn bond_match(&self, py: Python<'_>) -> Option<Py<PyAny>> {
        self.bond_match.as_ref().map(|value| value.clone_ref(py))
    }
    #[gen_stub(skip)]
    fn __traverse__(&self, visit: pyo3::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        visit.call(&self.final_match)?;
        visit.call(&self.atom_match)?;
        visit.call(&self.bond_match)
    }
    #[gen_stub(skip)]
    fn __clear__(&mut self) {
        self.final_match = None;
        self.atom_match = None;
        self.bond_match = None;
    }
    #[getter]
    fn max_matches(&self) -> usize {
        self.inner.max_matches
    }
    #[getter]
    fn uniquify(&self) -> bool {
        self.inner.uniquify
    }
    #[getter]
    fn use_chirality(&self) -> bool {
        self.inner.use_chirality
    }
    #[getter]
    fn use_enhanced_stereo(&self) -> bool {
        self.inner.use_enhanced_stereo
    }
    #[getter]
    fn specified_stereo_query_matches_unspecified(&self) -> bool {
        self.inner.specified_stereo_query_matches_unspecified
    }
    #[getter]
    fn use_query_query_matches(&self) -> bool {
        self.inner.use_query_query_matches
    }
    #[getter]
    fn recursion_possible(&self) -> bool {
        self.inner.recursion_possible
    }
    #[getter]
    fn max_recursive_matches(&self) -> usize {
        self.inner.max_recursive_matches
    }
    #[getter]
    fn num_threads(&self) -> i32 {
        self.inner.num_threads
    }
    #[getter]
    fn aromatic_matches_conjugated(&self) -> bool {
        self.inner.aromatic_matches_conjugated
    }
    #[getter]
    fn aromatic_matches_single_or_double(&self) -> bool {
        self.inner.aromatic_matches_single_or_double
    }
    #[getter]
    fn atom_properties(&self) -> Vec<String> {
        self.inner.atom_properties.clone()
    }
    #[getter]
    fn bond_properties(&self) -> Vec<String> {
        self.inner.bond_properties.clone()
    }
    #[getter]
    fn extra_atom_check_overrides_default_check(&self) -> bool {
        self.inner.extra_atom_check_overrides_default_check
    }
    #[getter]
    fn extra_bond_check_overrides_default_check(&self) -> bool {
        self.inner.extra_bond_check_overrides_default_check
    }
    #[getter]
    fn use_generic_matchers(&self) -> bool {
        self.inner.use_generic_matchers
    }
    fn __repr__(&self) -> String {
        format!("SubstructMatchParams({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smarts(py: Python<'_>, text: &str) -> PyResult<QueryGraph> {
    ck::parse_smarts(text)
        .map(|inner| QueryGraph { inner })
        .map_err(|e| parse_pyerr(py, e))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smarts_with_params(
    py: Python<'_>,
    text: &str,
    params: &SmartsParseParams,
) -> PyResult<QueryGraph> {
    ck::parse_smarts_with_params(text, &params.inner)
        .map(|inner| QueryGraph { inner })
        .map_err(|e| parse_pyerr(py, e))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn compile_query(py: Python<'_>, query: &QueryGraph) -> PyResult<CompiledQuery> {
    ck::compile_query(&query.inner)
        .map(|inner| CompiledQuery { inner })
        .map_err(|e| {
            crate::canonical_values::annotate(
                py,
                QueryCompileError::new_err(e.to_string()),
                "search",
                "InvalidGraph",
                &e,
            )
        })
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn write_smarts(
    py: Python<'_>,
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> PyResult<String> {
    ck::write_smarts(&query.inner, &params.inner)
        .map_err(|e| write_pyerr(py, e))
        .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn write_cx_smarts(
    py: Python<'_>,
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> PyResult<String> {
    ck::write_cx_smarts(&query.inner, &params.inner)
        .map_err(|e| write_pyerr(py, e))
        .and_then(|text| crate::canonical_sdf::decode_source_text(py, &text))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<QueryAtom>()?;
    module.add_class::<QueryGraph>()?;
    module.add_class::<CompiledQuery>()?;
    module.add_class::<MatchResult>()?;
    module.add_class::<SmartsParseParams>()?;
    module.add_class::<SmartsWriteParams>()?;
    module.add_class::<SubstructMatchParams>()?;
    module.add(
        "SmartsParseError",
        module.py().get_type::<SmartsParseError>(),
    )?;
    module.add(
        "SmartsWriteError",
        module.py().get_type::<SmartsWriteError>(),
    )?;
    module.add(
        "SubstructMatchError",
        module.py().get_type::<SubstructMatchError>(),
    )?;
    module.add(
        "QueryCompileError",
        module.py().get_type::<QueryCompileError>(),
    )?;
    module.add("MatchError", module.py().get_type::<MatchError>())?;
    module.add_function(wrap_pyfunction!(parse_smarts, module)?)?;
    module.add_function(wrap_pyfunction!(parse_smarts_with_params, module)?)?;
    module.add_function(wrap_pyfunction!(compile_query, module)?)?;
    module.add_function(wrap_pyfunction!(write_smarts, module)?)?;
    module.add_function(wrap_pyfunction!(write_cx_smarts, module)?)?;
    Ok(())
}
