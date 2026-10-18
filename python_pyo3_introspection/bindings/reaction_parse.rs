//! SMIRKS parsing projects the Rust facade; no parser or chemistry here.
use crate::error_framework::{Source, error_projection};
use crate::params::ReactionParseParams;
use crate::property_values::decode_source_text;
use crate::query_model_errors::QueryGraphError;
use crate::reaction::Reaction;
use crate::reaction_validation::ReactionValidationError;
use crate::reaction_values::ReactionRole;
use crate::smarts_errors::{CxQueryLoweringError, SmartsParseError};
use crate::smiles_parse_errors::SmilesParseError;
use crate::stereo_order_errors::StereoOrderError;
use crate::utf8_errors::Utf8Error;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::error::Error;

#[doc = r#"{"role":"error_record"}"#]
#[pyclass(module = "fixed_projection._native", extends=PyValueError)]
pub struct CxParseError {
    inner: ck::CxParseError,
}
#[pymethods]
impl CxParseError {
    #[getter]
    fn domain(&self) -> &'static str {
        "cx"
    }
    #[getter]
    fn offset(&self) -> usize {
        self.inner.offset
    }
    #[getter]
    fn message(&self) -> &str {
        &self.inner.message
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
impl CxParseError {
    fn from_rust(inner: ck::CxParseError) -> PyErr {
        let ck::CxParseError {
            offset: _,
            message: _,
        } = &inner;
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => PyErr::from_value(value.into_bound(py).into_any()),
            Err(error) => error,
        })
    }
}
type ParseSource = Source<
    Source<
        Py<SmartsParseError>,
        Source<Py<SmilesParseError>, Source<Py<Utf8Error>, Py<QueryGraphError>>>,
    >,
    Source<
        Py<ReactionValidationError>,
        Source<Py<CxParseError>, Source<Py<CxQueryLoweringError>, Py<StereoOrderError>>>,
    >,
>;
error_projection!(ReactionParseError, ReactionParseErrorKind, ck::ReactionParseError, "reaction", {
    ParserCleanup { role, template, source },
    TemplateAtomBounds { role, template, atom, atom_count: usize },
    StereoDegreeArithmetic { reactant: u32, product },
    StereoInversionFlag { product, atom, flag: Option<i32> },
    TemplateProperty { source }, Separators { count }, MultiStep { count },
    ComponentBounds { start: usize, end: usize }, OffsetOverflow,
    Smarts { role, template, text, source }, Smiles { role, template, text, source },
    ComponentEncoding { role, template, text, source }, Model { role, template, source },
    AgentFragments { source }, CxParse { role, template, start_atom, start_bond, source },
    CxLowering { role, template, start_atom, start_bond, source }, StereoOrder { product, atom, source },
}, cause = |inner: &ck::ReactionParseError| match inner {
    ck::ReactionParseError::ParserCleanup { source, .. } | ck::ReactionParseError::Smarts { source, .. }
    | ck::ReactionParseError::AgentFragments { source } => Some(SmartsParseError::from_rust(source.clone())),
    ck::ReactionParseError::Smiles { source, .. } => Some(SmilesParseError::from_rust(source.clone())),
    ck::ReactionParseError::ComponentEncoding { source, .. } => Some(Utf8Error::from_rust(*source)),
    ck::ReactionParseError::Model { source, .. } => Some(QueryGraphError::from_rust(source.clone())),
    ck::ReactionParseError::TemplateProperty { source } => Some(ReactionValidationError::from_rust(source.clone())),
    ck::ReactionParseError::CxParse { source, .. } => Some(CxParseError::from_rust(source.clone())),
    ck::ReactionParseError::CxLowering { source, .. } => Some(CxQueryLoweringError::from_rust(source.clone())),
    ck::ReactionParseError::StereoOrder { source, .. } => Some(StereoOrderError::from_rust(source.clone())),
    ck::ReactionParseError::TemplateAtomBounds { .. } | ck::ReactionParseError::StereoDegreeArithmetic { .. }
    | ck::ReactionParseError::StereoInversionFlag { .. } | ck::ReactionParseError::Separators { .. }
    | ck::ReactionParseError::MultiStep { .. } | ck::ReactionParseError::ComponentBounds { .. }
    | ck::ReactionParseError::OffsetOverflow => None,
}, source = Option<ParseSource>, methods = {
    #[getter]
    fn role(&self) -> Option<ReactionRole> { match self.inner {
        ck::ReactionParseError::ParserCleanup { role, .. } | ck::ReactionParseError::TemplateAtomBounds { role, .. }
        | ck::ReactionParseError::Smarts { role, .. } | ck::ReactionParseError::Smiles { role, .. }
        | ck::ReactionParseError::ComponentEncoding { role, .. } | ck::ReactionParseError::Model { role, .. }
        | ck::ReactionParseError::CxParse { role, .. } | ck::ReactionParseError::CxLowering { role, .. } => Some(role.into()), _ => None,
    } }
    #[getter]
    fn template(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::ParserCleanup { template, .. } | ck::ReactionParseError::TemplateAtomBounds { template, .. }
        | ck::ReactionParseError::Smarts { template, .. } | ck::ReactionParseError::Smiles { template, .. }
        | ck::ReactionParseError::ComponentEncoding { template, .. } | ck::ReactionParseError::Model { template, .. }
        | ck::ReactionParseError::CxParse { template, .. } | ck::ReactionParseError::CxLowering { template, .. } => Some(template), _ => None,
    } }
    #[getter]
    fn atom(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::TemplateAtomBounds { atom, .. } | ck::ReactionParseError::StereoInversionFlag { atom, .. }
        | ck::ReactionParseError::StereoOrder { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    fn product(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::StereoDegreeArithmetic { product, .. } => Some(product as usize),
        ck::ReactionParseError::StereoInversionFlag { product, .. } | ck::ReactionParseError::StereoOrder { product, .. } => Some(product), _ => None,
    } }
    #[getter]
    fn count(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::Separators { count } | ck::ReactionParseError::MultiStep { count } => Some(count), _ => None,
    } }
    #[getter]
    fn start_atom(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::CxParse { start_atom, .. } | ck::ReactionParseError::CxLowering { start_atom, .. } => Some(start_atom), _ => None,
    } }
    #[getter]
    fn start_bond(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::CxParse { start_bond, .. } | ck::ReactionParseError::CxLowering { start_bond, .. } => Some(start_bond), _ => None,
    } }
    #[getter]
    fn text(&self, py: Python<'_>) -> PyResult<Option<String>> { match &self.inner {
        ck::ReactionParseError::Smarts { text, .. } | ck::ReactionParseError::Smiles { text, .. }
        | ck::ReactionParseError::ComponentEncoding { text, .. } => decode_source_text(py, text).map(Some), _ => Ok(None),
    } }
});

pub(crate) fn keywords(
    use_smiles: Option<bool>,
    sanitize: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
) -> ReactionParseParams {
    let defaults = ck::ReactionParseParams::default();
    ReactionParseParams {
        inner: ck::ReactionParseParams {
            use_smiles: use_smiles.unwrap_or(defaults.use_smiles),
            sanitize: sanitize.unwrap_or(defaults.sanitize),
            replacements: replacements.unwrap_or(defaults.replacements),
            allow_cxsmiles: allow_cxsmiles.unwrap_or(defaults.allow_cxsmiles),
            strict_cxsmiles: strict_cxsmiles.unwrap_or(defaults.strict_cxsmiles),
        },
    }
}

#[pyfunction]
#[pyo3(signature=(text, params=None, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
pub fn parse_smirks(
    text: &str,
    params: Option<&ReactionParseParams>,
    use_smiles: Option<bool>,
    sanitize: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
) -> PyResult<Reaction> {
    let configured = use_smiles.is_some()
        || sanitize.is_some()
        || replacements.is_some()
        || allow_cxsmiles.is_some()
        || strict_cxsmiles.is_some();
    if params.is_some() && configured {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    if configured {
        _parse_smirks_keywords(
            text,
            use_smiles,
            sanitize,
            replacements,
            allow_cxsmiles,
            strict_cxsmiles,
        )
    } else {
        _parse_smirks_params(text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"parse_smirks"}"#]
#[pyo3(signature=(text, params=None))]
pub fn _parse_smirks_params(
    text: &str,
    params: Option<&ReactionParseParams>,
) -> PyResult<Reaction> {
    match params {
        Some(params) => ck::parse_smirks_with_params(text, &params.inner),
        None => ck::parse_smirks(text),
    }
    .map(Reaction::from_rust)
    .map_err(ReactionParseError::from_rust)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"parse_smirks"}"#]
#[pyo3(signature=(text, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
pub fn _parse_smirks_keywords(
    text: &str,
    use_smiles: Option<bool>,
    sanitize: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
) -> PyResult<Reaction> {
    let params = keywords(
        use_smiles,
        sanitize,
        replacements,
        allow_cxsmiles,
        strict_cxsmiles,
    );
    _parse_smirks_params(text, Some(&params))
}
#[pyfunction]
pub fn parse_smirks_with_params(text: &str, params: &ReactionParseParams) -> PyResult<Reaction> {
    ck::parse_smirks_with_params(text, &params.inner)
        .map(Reaction::from_rust)
        .map_err(ReactionParseError::from_rust)
}

#[pyfunction]
#[pyo3(signature=(text, params=None, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
pub fn reaction_from_smirks(
    text: &str,
    params: Option<&ReactionParseParams>,
    use_smiles: Option<bool>,
    sanitize: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
) -> PyResult<Reaction> {
    let configured = use_smiles.is_some()
        || sanitize.is_some()
        || replacements.is_some()
        || allow_cxsmiles.is_some()
        || strict_cxsmiles.is_some();
    if params.is_some() && configured {
        return Err(pyo3::exceptions::PyTypeError::new_err(
            "Params and keyword configuration are mutually exclusive",
        ));
    }
    if configured {
        _reaction_from_smirks_keywords(
            text,
            use_smiles,
            sanitize,
            replacements,
            allow_cxsmiles,
            strict_cxsmiles,
        )
    } else {
        _reaction_from_smirks_params(text, params)
    }
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"reaction_from_smirks"}"#]
#[pyo3(signature=(text, params=None))]
pub fn _reaction_from_smirks_params(
    text: &str,
    params: Option<&ReactionParseParams>,
) -> PyResult<Reaction> {
    match params {
        Some(params) => ck::Reaction::from_smirks_with_params(text, &params.inner),
        None => ck::reaction_from_smirks(text),
    }
    .map(Reaction::from_rust)
    .map_err(ReactionParseError::from_rust)
}
#[pyfunction]
#[doc = r#"{"role":"overload","target":"reaction_from_smirks"}"#]
#[pyo3(signature=(text, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
pub fn _reaction_from_smirks_keywords(
    text: &str,
    use_smiles: Option<bool>,
    sanitize: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
) -> PyResult<Reaction> {
    let params = keywords(
        use_smiles,
        sanitize,
        replacements,
        allow_cxsmiles,
        strict_cxsmiles,
    );
    _reaction_from_smirks_params(text, Some(&params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::query_compiler::QueryGraph;
    use pyo3::types::{PyDict, PyModule};

    fn exact(actual: &Reaction, expected: &ck::Reaction) {
        assert_eq!(
            actual.inner().reactant_templates(),
            expected.reactant_templates()
        );
        assert_eq!(
            actual.inner().product_templates(),
            expected.product_templates()
        );
        assert_eq!(actual.inner().agent_templates(), expected.agent_templates());
        assert_eq!(actual.inner().is_initialized(), expected.is_initialized());
        assert_eq!(
            actual.inner().implicit_properties(),
            expected.implicit_properties()
        );
    }

    #[test]
    fn complete_smirks_family_forwards_all_forms_templates_and_concrete_parse_context() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "reaction_parse_test")?;
            module.add_class::<Reaction>()?;
            module.add_class::<ReactionParseParams>()?;
            module.add_function(wrap_pyfunction!(parse_smirks, &module)?)?;
            module.add_function(wrap_pyfunction!(parse_smirks_with_params, &module)?)?;
            module.add_function(wrap_pyfunction!(reaction_from_smirks, &module)?)?;
            let class = py.get_type::<Reaction>();
            for (text, native_params) in [
                ("[C:1]>>[O:1]", ck::ReactionParseParams::default()),
                (
                    "CC>[Na+]>CO",
                    ck::ReactionParseParams {
                        use_smiles: true,
                        sanitize: true,
                        ..Default::default()
                    },
                ),
                (
                    "{X}>>{X}",
                    ck::ReactionParseParams {
                        replacements: BTreeMap::from([("{X}".into(), "[C:1]".into())]),
                        ..Default::default()
                    },
                ),
            ] {
                let expected = ck::parse_smirks_with_params(text, &native_params).unwrap();
                let params = Py::new(
                    py,
                    ReactionParseParams {
                        inner: native_params.clone(),
                    },
                )?;
                let kwargs = PyDict::new(py);
                kwargs.set_item("use_smiles", native_params.use_smiles)?;
                kwargs.set_item("sanitize", native_params.sanitize)?;
                kwargs.set_item("replacements", &native_params.replacements)?;
                kwargs.set_item("allow_cxsmiles", native_params.allow_cxsmiles)?;
                kwargs.set_item("strict_cxsmiles", native_params.strict_cxsmiles)?;
                let mut values = vec![
                    class.call_method1("from_smirks", (text, params.clone_ref(py)))?,
                    class.call_method1("from_smirks_with_params", (text, params.clone_ref(py)))?,
                    class.call_method("from_smirks", (text,), Some(&kwargs))?,
                    module
                        .getattr("parse_smirks")?
                        .call1((text, params.clone_ref(py)))?,
                    module
                        .getattr("parse_smirks_with_params")?
                        .call1((text, params.clone_ref(py)))?,
                    module
                        .getattr("parse_smirks")?
                        .call((text,), Some(&kwargs))?,
                    module
                        .getattr("reaction_from_smirks")?
                        .call1((text, params.clone_ref(py)))?,
                    module
                        .getattr("reaction_from_smirks")?
                        .call((text,), Some(&kwargs))?,
                ];
                if native_params == ck::ReactionParseParams::default() {
                    values.push(class.call_method1("from_smirks", (text,))?);
                    values.push(module.getattr("parse_smirks")?.call1((text,))?);
                    values.push(module.getattr("reaction_from_smirks")?.call1((text,))?);
                }
                for value in &values {
                    exact(&*value.extract::<PyRef<Reaction>>()?, &expected);
                }
                for callable in [
                    class.getattr("from_smirks")?,
                    module.getattr("parse_smirks")?,
                    module.getattr("reaction_from_smirks")?,
                ] {
                    assert!(
                        callable
                            .call((text, params.clone_ref(py)), Some(&kwargs))
                            .unwrap_err()
                            .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                    );
                }
                let retained = values[0]
                    .call_method1("reactant_template", (0,))?
                    .extract::<Py<QueryGraph>>()?;
                let expected_query = expected.reactant_template(0).unwrap().clone();
                drop(values);
                drop(params);
                assert_eq!(retained.borrow(py).inner(), &expected_query);
            }
            let partial = PyDict::new(py);
            partial.set_item("use_smiles", true)?;
            let actual = class.call_method("from_smirks", ("CC>>CO",), Some(&partial))?;
            let expected = ck::Reaction::from_smirks_with_params(
                "CC>>CO",
                &ck::ReactionParseParams {
                    use_smiles: true,
                    ..Default::default()
                },
            )
            .unwrap();
            exact(&*actual.extract::<PyRef<Reaction>>()?, &expected);

            for callable in [
                class.getattr("from_smirks")?,
                module.getattr("parse_smirks")?,
                module.getattr("reaction_from_smirks")?,
            ] {
                let error = callable.call1(("CC",)).unwrap_err();
                assert!(error.is_instance_of::<ReactionParseError>(py));
                assert!(
                    error
                        .value(py)
                        .getattr("variant")?
                        .extract::<ReactionParseErrorKind>()?
                        == ReactionParseErrorKind::Separators
                );
                assert_eq!(error.value(py).getattr("count")?.extract::<usize>()?, 0);
                assert!(error.value(py).getattr("source")?.is_none());
                assert!(error.value(py).setattr("count", 3).is_err());
                let error = callable.call1(("[C:1]>>[",)).unwrap_err();
                assert!(error.is_instance_of::<ReactionParseError>(py));
                let value = error.value(py);
                assert!(
                    value
                        .getattr("source")?
                        .is_instance_of::<SmartsParseError>()
                );
                assert!(value.getattr("source")?.is(&value.getattr("__cause__")?));
                assert!(value.getattr("role")?.extract::<ReactionRole>()? == ReactionRole::Product);
                assert_eq!(value.getattr("template")?.extract::<usize>()?, 0);
                assert_eq!(value.getattr("text")?.extract::<String>()?, "[");
            }
            let native = ck::ReactionParseError::CxParse {
                role: ck::ReactionRole::Agent,
                template: usize::MAX,
                start_atom: usize::MAX,
                start_bond: 7,
                source: ck::CxParseError::new(usize::MAX, "message\0text"),
            };
            assert!(native.source().unwrap().is::<ck::CxParseError>());
            let display = native.to_string();
            let error = ReactionParseError::from_rust(native);
            let value = error.value(py);
            assert_eq!(value.str()?.to_str()?, display);
            assert_eq!(value.getattr("start_atom")?.extract::<usize>()?, usize::MAX);
            assert_eq!(value.getattr("template")?.extract::<usize>()?, usize::MAX);
            let source = value.getattr("source")?;
            assert!(source.is_instance_of::<CxParseError>());
            assert_eq!(source.getattr("offset")?.extract::<usize>()?, usize::MAX);
            assert_eq!(
                source.getattr("message")?.extract::<String>()?,
                "message\0text"
            );
            assert!(source.getattr("__cause__")?.is_none());
            assert!(source.setattr("offset", 0).is_err());
            Ok(())
        })
        .unwrap();
    }
}
