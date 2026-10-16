//! SMILES configuration, module calls and an indexed molecule collection.
use crate::indexing;
use crate::native::Molecule;
use crate::parameters::parameters;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::collections::BTreeMap;

parameters!(SmilesParseParams, ck::SmilesParseParams, {
    sanitize: bool,
    allow_cxsmiles: bool,
    strict_cxsmiles: bool,
    parse_name: bool,
    remove_hs: bool,
    skip_cleanup: bool,
    debug_parse: bool,
    replacements: BTreeMap<String, String>,
});

#[pyfunction]
#[pyo3(signature = (text, params=None, *, sanitize=None, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, remove_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
pub fn mol_from_smiles_configured(
    text: &str,
    params: Option<&SmilesParseParams>,
    sanitize: Option<bool>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    remove_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
) -> PyResult<Molecule> {
    // Only these explicitly bound keywords are supported. The generator does
    // not infer call options from the fields or constructor of SmilesParseParams.
    if let Some(params) = params {
        if sanitize.is_some()
            || allow_cxsmiles.is_some()
            || strict_cxsmiles.is_some()
            || parse_name.is_some()
            || remove_hs.is_some()
            || skip_cleanup.is_some()
            || debug_parse.is_some()
            || replacements.is_some()
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        _mol_from_smiles_configured_params(text, Some(params))
    } else {
        _mol_from_smiles_configured_keywords(
            text,
            sanitize,
            allow_cxsmiles,
            strict_cxsmiles,
            parse_name,
            remove_hs,
            skip_cleanup,
            debug_parse,
            replacements,
        )
    }
}

#[pyfunction]
#[doc = r#"{"role": "overload", "target": "mol_from_smiles_configured"}"#]
#[pyo3(signature = (text, params=None))]
pub fn _mol_from_smiles_configured_params(
    text: &str,
    params: Option<&SmilesParseParams>,
) -> PyResult<Molecule> {
    let defaults;
    let effective = if let Some(params) = params {
        &params.inner
    } else {
        defaults = ck::SmilesParseParams::default();
        &defaults
    };
    ck::Molecule::from_smiles_with_params(text, effective)
        .map(|inner| Molecule { inner })
        .map_err(crate::errors::SmilesError::from_rust)
}

#[pyfunction]
#[doc = r#"{"role": "overload", "target": "mol_from_smiles_configured"}"#]
#[pyo3(signature = (text, *, sanitize=None, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, remove_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
pub fn _mol_from_smiles_configured_keywords(
    text: &str,
    sanitize: Option<bool>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    remove_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
) -> PyResult<Molecule> {
    let mut value = ck::SmilesParseParams::default();
    if let Some(option) = sanitize {
        value.sanitize = option;
    }
    if let Some(option) = allow_cxsmiles {
        value.allow_cxsmiles = option;
    }
    if let Some(option) = strict_cxsmiles {
        value.strict_cxsmiles = option;
    }
    if let Some(option) = parse_name {
        value.parse_name = option;
    }
    if let Some(option) = remove_hs {
        value.remove_hs = option;
    }
    if let Some(option) = skip_cleanup {
        value.skip_cleanup = option;
    }
    if let Some(option) = debug_parse {
        value.debug_parse = option;
    }
    if let Some(option) = replacements {
        value.replacements = option;
    }
    _mol_from_smiles_configured_params(text, Some(&SmilesParseParams { inner: value }))
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct MoleculeBatch {
    values: Vec<Option<Py<Molecule>>>,
}

#[pymethods]
impl MoleculeBatch {
    #[new]
    fn new(values: Vec<Option<Py<Molecule>>>) -> Self {
        Self { values }
    }
    fn __len__(&self) -> usize {
        self.values.len()
    }
    fn __getitem__(
        &self,
        py: Python<'_>,
        index: indexing::Index<'_>,
    ) -> PyResult<indexing::Indexed<Option<Py<Molecule>>, Py<Self>>> {
        indexing::get(py, index, &self.values, |values| Self { values })
    }
}
