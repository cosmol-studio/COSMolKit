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
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::error::Error;

crate::error_framework::error_projection! {record
#[doc = r#"{"role":"error_record"}"#]
pub struct CxParseError {
    inner: ck::CxParseError,
}
}
crate::error_framework::error_projection! {methods impl CxParseError {
    #[getter]
    /// Native subsystem reporting this failure, preserved independently of the human-readable message.
    fn domain(&self) -> &'static str {
        "cx"
    }
    #[getter]
    #[doc = "Byte offset of the CXSMILES parsing failure."] fn offset(&self) -> usize {
        self.inner.offset
    }
    #[getter]
    #[doc = "Original CXSMILES parser explanation."] fn message(&self) -> &str {
        &self.inner.message
    }
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
}
impl CxParseError {
    fn from_rust(inner: ck::CxParseError) -> PyErr {
        let ck::CxParseError {
            offset: _,
            message: _,
        } = &inner;
        Python::attach(|py| match Py::new(py, Self { inner }) {
            Ok(value) => crate::error_framework::from_payload(value),
            Err(error) => error,
        })
    }
}
type ParseSource = Source<
    Source<
        crate::error_framework::ErrorObject<SmartsParseError>,
        Source<
            crate::error_framework::ErrorObject<SmilesParseError>,
            Source<
                crate::error_framework::ErrorObject<Utf8Error>,
                crate::error_framework::ErrorObject<QueryGraphError>,
            >,
        >,
    >,
    Source<
        crate::error_framework::ErrorObject<ReactionValidationError>,
        Source<
            crate::error_framework::ErrorObject<CxParseError>,
            Source<
                crate::error_framework::ErrorObject<CxQueryLoweringError>,
                crate::error_framework::ErrorObject<StereoOrderError>,
            >,
        >,
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
    #[doc = "Role of the affected reaction template or transformation participant. Present for `ReactionParseError.ParserCleanup`, `ReactionParseError.TemplateAtomBounds`, `ReactionParseError.Smarts`, `ReactionParseError.Smiles`, `ReactionParseError.ComponentEncoding`, `ReactionParseError.Model`, `ReactionParseError.CxParse`, `ReactionParseError.CxLowering`; None when unavailable."] fn role(&self) -> Option<ReactionRole> { match self.inner {
        ck::ReactionParseError::ParserCleanup { role, .. } | ck::ReactionParseError::TemplateAtomBounds { role, .. }
        | ck::ReactionParseError::Smarts { role, .. } | ck::ReactionParseError::Smiles { role, .. }
        | ck::ReactionParseError::ComponentEncoding { role, .. } | ck::ReactionParseError::Model { role, .. }
        | ck::ReactionParseError::CxParse { role, .. } | ck::ReactionParseError::CxLowering { role, .. } => Some(role.into()), _ => None,
    } }
    #[getter]
    #[doc = "Index of the affected reaction template. Present for `ReactionParseError.ParserCleanup`, `ReactionParseError.TemplateAtomBounds`, `ReactionParseError.Smarts`, `ReactionParseError.Smiles`, `ReactionParseError.ComponentEncoding`, `ReactionParseError.Model`, `ReactionParseError.CxParse`, `ReactionParseError.CxLowering`; None when unavailable."] fn template(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::ParserCleanup { template, .. } | ck::ReactionParseError::TemplateAtomBounds { template, .. }
        | ck::ReactionParseError::Smarts { template, .. } | ck::ReactionParseError::Smiles { template, .. }
        | ck::ReactionParseError::ComponentEncoding { template, .. } | ck::ReactionParseError::Model { template, .. }
        | ck::ReactionParseError::CxParse { template, .. } | ck::ReactionParseError::CxLowering { template, .. } => Some(template), _ => None,
    } }
    #[getter]
    #[doc = "Atom index involved in the failure. Present for `ReactionParseError.TemplateAtomBounds`, `ReactionParseError.StereoInversionFlag`, `ReactionParseError.StereoOrder`; None when unavailable."] fn atom(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::TemplateAtomBounds { atom, .. } | ck::ReactionParseError::StereoInversionFlag { atom, .. }
        | ck::ReactionParseError::StereoOrder { atom, .. } => Some(atom.index()), _ => None,
    } }
    #[getter]
    #[doc = "Product-template index associated with the failure. Present for `ReactionParseError.StereoDegreeArithmetic`, `ReactionParseError.StereoInversionFlag`, `ReactionParseError.StereoOrder`; None when unavailable."] fn product(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::StereoDegreeArithmetic { product, .. } => Some(product as usize),
        ck::ReactionParseError::StereoInversionFlag { product, .. } | ck::ReactionParseError::StereoOrder { product, .. } => Some(product), _ => None,
    } }
    #[getter]
    #[doc = "Count observed when the failure occurred. Present for `ReactionParseError.Separators`, `ReactionParseError.MultiStep`; None when unavailable."] fn count(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::Separators { count } | ck::ReactionParseError::MultiStep { count } => Some(count), _ => None,
    } }
    #[getter]
    #[doc = "Starting atom offset for the reaction template. Present for `ReactionParseError.CxParse`, `ReactionParseError.CxLowering`; None when unavailable."] fn start_atom(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::CxParse { start_atom, .. } | ck::ReactionParseError::CxLowering { start_atom, .. } => Some(start_atom), _ => None,
    } }
    #[getter]
    #[doc = "Starting bond offset for the reaction template. Present for `ReactionParseError.CxParse`, `ReactionParseError.CxLowering`; None when unavailable."] fn start_bond(&self) -> Option<usize> { match self.inner {
        ck::ReactionParseError::CxParse { start_bond, .. } | ck::ReactionParseError::CxLowering { start_bond, .. } => Some(start_bond), _ => None,
    } }
    #[getter]
    #[doc = "Original input text retained by the error. Present for `ReactionParseError.Smarts`, `ReactionParseError.Smiles`, `ReactionParseError.ComponentEncoding`; None when unavailable."] fn text(&self, py: Python<'_>) -> PyResult<Option<String>> { match &self.inner {
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
) -> ck::ReactionParseParams {
    let defaults = ck::ReactionParseParams::default();
    ck::ReactionParseParams {
        use_smiles: use_smiles.unwrap_or(defaults.use_smiles),
        sanitize: sanitize.unwrap_or(defaults.sanitize),
        replacements: replacements.unwrap_or(defaults.replacements),
        allow_cxsmiles: allow_cxsmiles.unwrap_or(defaults.allow_cxsmiles),
        strict_cxsmiles: strict_cxsmiles.unwrap_or(defaults.strict_cxsmiles),
    }
}

#[pyfunction]
#[pyo3(signature=(text, params=crate::arguments::Argument::Missing, *, use_smiles=crate::arguments::Argument::Missing, sanitize=crate::arguments::Argument::Missing, replacements=crate::arguments::Argument::Missing, allow_cxsmiles=crate::arguments::Argument::Missing, strict_cxsmiles=crate::arguments::Argument::Missing))]
/// Parse reaction SMIRKS into a Reaction using the selected native parsing parameters.
pub fn parse_smirks(
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, ReactionParseParams>>>,
    use_smiles: crate::arguments::Argument<Option<bool>>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    replacements: crate::arguments::Argument<Option<BTreeMap<String, String>>>,
    allow_cxsmiles: crate::arguments::Argument<Option<bool>>,
    strict_cxsmiles: crate::arguments::Argument<Option<bool>>,
) -> PyResult<Reaction> {
    let keyword_form = use_smiles.is_supplied()
        || sanitize.is_supplied()
        || replacements.is_supplied()
        || allow_cxsmiles.is_supplied()
        || strict_cxsmiles.is_supplied();
    crate::arguments::configuration_arguments!(params; use_smiles, sanitize, replacements, allow_cxsmiles, strict_cxsmiles);
    let params = params.as_deref();
    if keyword_form {
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
    let params = Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
    match params {
        Some(params) => ck::parse_smirks_with_params(text, &params),
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
    ck::parse_smirks_with_params(text, &params)
        .map(Reaction::from_rust)
        .map_err(ReactionParseError::from_rust)
}
#[pyfunction]
/// Parse reaction SMIRKS into a Reaction with explicit ReactionParseParams.
pub fn parse_smirks_with_params(text: &str, params: &ReactionParseParams) -> PyResult<Reaction> {
    _parse_smirks_params(text, Some(params))
}

#[pyfunction]
#[pyo3(signature=(text, params=crate::arguments::Argument::Missing, *, use_smiles=crate::arguments::Argument::Missing, sanitize=crate::arguments::Argument::Missing, replacements=crate::arguments::Argument::Missing, allow_cxsmiles=crate::arguments::Argument::Missing, strict_cxsmiles=crate::arguments::Argument::Missing))]
/// Construct a Reaction from SMIRKS using the same canonical parser as parse_smirks.
pub fn reaction_from_smirks(
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, ReactionParseParams>>>,
    use_smiles: crate::arguments::Argument<Option<bool>>,
    sanitize: crate::arguments::Argument<Option<bool>>,
    replacements: crate::arguments::Argument<Option<BTreeMap<String, String>>>,
    allow_cxsmiles: crate::arguments::Argument<Option<bool>>,
    strict_cxsmiles: crate::arguments::Argument<Option<bool>>,
) -> PyResult<Reaction> {
    let keyword_form = use_smiles.is_supplied()
        || sanitize.is_supplied()
        || replacements.is_supplied()
        || allow_cxsmiles.is_supplied()
        || strict_cxsmiles.is_supplied();
    crate::arguments::configuration_arguments!(params; use_smiles, sanitize, replacements, allow_cxsmiles, strict_cxsmiles);
    let params = params.as_deref();
    if keyword_form {
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
    let params = Python::attach(|py| params.map(|value| value.to_rust(py)).transpose())?;
    match params {
        Some(params) => ck::Reaction::from_smirks_with_params(text, &params),
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
    ck::Reaction::from_smirks_with_params(text, &params)
        .map(Reaction::from_rust)
        .map_err(ReactionParseError::from_rust)
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
                    ReactionParseParams::from_rust(py, native_params.clone())?,
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
                assert!(error.is_instance(py, <ReactionParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
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
                assert!(error.is_instance(py, <ReactionParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
                let value = error.value(py);
                assert!(
                    value
                        .getattr("source")?
                        .is_instance(<SmartsParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
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
            assert!(source.is_instance(<CxParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
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
