//! Canonical detached query and search projections; chemistry lives behind cosmolkit.
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{gen_stub_pyclass, gen_stub_pyfunction, gen_stub_pymethods};
use std::collections::BTreeMap;

pyo3::create_exception!(cosmolkit, SmartsParseError, PyValueError);
pyo3::create_exception!(cosmolkit, SmartsWriteError, PyValueError);
pyo3::create_exception!(cosmolkit, SubstructMatchError, PyValueError);
pyo3::create_exception!(cosmolkit, QueryCompileError, PyValueError);
pyo3::create_exception!(cosmolkit, MatchError, PyValueError);

pub(crate) fn parse_pyerr(py: Python<'_>, source: ck::SmartsParseError) -> PyErr {
    use ck::SmartsParseError as E;
    let kind = match &source {
        E::UnclosedBracket(_) => "UnclosedBracket",
        E::UnexpectedCharacter { .. } => "UnexpectedCharacter",
        E::UnexpectedEnd(_) => "UnexpectedEnd",
        E::InvalidAtomPrimitive { .. } => "InvalidAtomPrimitive",
        E::UnclosedParenthesis(_) => "UnclosedParenthesis",
        E::UnbalancedRingClosure(_) => "UnbalancedRingClosure",
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
                value.setattr("character", character.to_string())?;
                value.setattr("context", context)?;
            }
            E::InvalidAtomPrimitive { position, detail } => {
                value.setattr("position", position)?;
                value.setattr("detail", detail)?;
            }
            E::UnbalancedRingClosure(number) => value.setattr("ring", number)?,
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
        E::InvalidPropertyKind { .. } => "InvalidPropertyKind",
        E::Property(_) => "Property",
        E::InvalidGraph(_) => "InvalidGraph",
        E::QueryGraphTraversalUnsupported { .. } => "QueryGraphTraversalUnsupported",
        E::OrAboveAndBelowAnd => "OrAboveAndBelowAnd",
        E::UnknownCombination { .. } => "UnknownCombination",
        E::MissingRecursiveQueryMolecule => "MissingRecursiveQueryMolecule",
        E::UnsupportedBondDirection { .. } => "UnsupportedBondDirection",
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
        ck::SubstructMatchError::PeriodicTable(_) => "PeriodicTable",
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
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct QueryGraph {
    pub(crate) inner: ck::QueryGraph,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl QueryGraph {
    fn num_atoms(&self) -> usize {
        self.inner.num_atoms()
    }
    fn num_bonds(&self) -> usize {
        self.inner.num_bonds()
    }
    fn name(&self) -> Option<String> {
        self.inner.prop("_Name").map(str::to_owned)
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
#[pyclass(module = "cosmolkit", frozen)]
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
#[pyclass(module = "cosmolkit", frozen)]
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
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SmartsParseParams {
    pub(crate) inner: ck::SmartsParseParams,
}
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
                replacements: replacements.unwrap_or_default(),
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
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements.clone()
    }
    fn __repr__(&self) -> String {
        format!("SmartsParseParams({:?})", self.inner)
    }
}

#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SmartsWriteParams {
    pub(crate) inner: ck::SmartsWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SmartsWriteParams {
    #[new]
    #[pyo3(signature = (*, include_atom_maps=true, do_isomeric_smiles=true, include_dative_bonds=true, rooted_at_atom=None))]
    fn new(
        include_atom_maps: bool,
        do_isomeric_smiles: bool,
        include_dative_bonds: bool,
        rooted_at_atom: Option<usize>,
    ) -> Self {
        Self {
            inner: ck::SmartsWriteParams {
                include_atom_maps,
                do_isomeric_smiles,
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
    fn do_isomeric_smiles(&self) -> bool {
        self.inner.do_isomeric_smiles
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
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct SubstructMatchParams {
    pub(crate) inner: ck::SubstructMatchParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl SubstructMatchParams {
    #[new]
    #[pyo3(signature = (*, max_matches=1000, uniquify=true, use_chirality=false, use_enhanced_stereo=false, specified_stereo_query_matches_unspecified=false, use_query_query_matches=false, recursion_possible=true, max_recursive_matches=1000, num_threads=1, aromatic_matches_conjugated=false, aromatic_matches_single_or_double=false, atom_properties=None, bond_properties=None, extra_atom_check_overrides_default_check=false, extra_bond_check_overrides_default_check=false, use_generic_matchers=false))]
    fn new(
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
    ) -> Self {
        Self {
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
        }
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
    ck::search::parse_smarts(text)
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
    ck::search::parse_smarts_with_params(text, &params.inner)
        .map(|inner| QueryGraph { inner })
        .map_err(|e| parse_pyerr(py, e))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn compile_query(py: Python<'_>, query: &QueryGraph) -> PyResult<CompiledQuery> {
    ck::search::compile_query(&query.inner)
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
    ck::search::write_smarts(&query.inner, &params.inner).map_err(|e| write_pyerr(py, e))
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn write_cx_smarts(
    py: Python<'_>,
    query: &QueryGraph,
    params: &SmartsWriteParams,
) -> PyResult<String> {
    ck::search::write_cx_smarts(&query.inner, &params.inner).map_err(|e| write_pyerr(py, e))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
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
    let search = PyModule::new(module.py(), "cosmolkit.search")?;
    search.add_function(wrap_pyfunction!(parse_smarts, &search)?)?;
    search.add_function(wrap_pyfunction!(parse_smarts_with_params, &search)?)?;
    search.add_function(wrap_pyfunction!(compile_query, &search)?)?;
    search.add_function(wrap_pyfunction!(write_smarts, &search)?)?;
    search.add_function(wrap_pyfunction!(write_cx_smarts, &search)?)?;
    module.add("search", &search)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("cosmolkit.search", &search)?;
    Ok(())
}
