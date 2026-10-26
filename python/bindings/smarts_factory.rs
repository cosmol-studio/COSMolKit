//! Typed constructor routing; the sole native parser owns all SMARTS behavior.
use crate::notation_params::SmartsParseParams;
use crate::query_compiler::QueryGraph;
use crate::smarts_errors::SmartsParseError;
use cosmolkit as ck;
use pyo3::exceptions::PyTypeError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[pyfunction]
#[pyo3(signature = (text, params=crate::arguments::Argument::Missing, *, allow_cxsmiles=crate::arguments::Argument::Missing, strict_cxsmiles=crate::arguments::Argument::Missing, parse_name=crate::arguments::Argument::Missing, merge_hs=crate::arguments::Argument::Missing, skip_cleanup=crate::arguments::Argument::Missing, debug_parse=crate::arguments::Argument::Missing, replacements=crate::arguments::Argument::Missing))]
/// Parse SMARTS into a QueryGraph. Supply SmartsParseParams or keyword options, not both; syntax failures retain structured diagnostics.
pub fn parse_smarts(
    text: &str,
    params: crate::arguments::Argument<Option<PyRef<'_, SmartsParseParams>>>,
    allow_cxsmiles: crate::arguments::Argument<Option<bool>>,
    strict_cxsmiles: crate::arguments::Argument<Option<bool>>,
    parse_name: crate::arguments::Argument<Option<bool>>,
    merge_hs: crate::arguments::Argument<Option<bool>>,
    skip_cleanup: crate::arguments::Argument<Option<bool>>,
    debug_parse: crate::arguments::Argument<Option<bool>>,
    replacements: crate::arguments::Argument<Option<BTreeMap<String, String>>>,
) -> PyResult<QueryGraph> {
    let keyword_form = allow_cxsmiles.is_supplied()
        || strict_cxsmiles.is_supplied()
        || parse_name.is_supplied()
        || merge_hs.is_supplied()
        || skip_cleanup.is_supplied()
        || debug_parse.is_supplied()
        || replacements.is_supplied();
    crate::arguments::configuration_arguments!(params; allow_cxsmiles, strict_cxsmiles, parse_name, merge_hs, skip_cleanup, debug_parse, replacements);
    let params = params.as_deref();
    if keyword_form {
        _parse_smarts_keywords(
            text,
            allow_cxsmiles,
            strict_cxsmiles,
            parse_name,
            merge_hs,
            skip_cleanup,
            debug_parse,
            replacements,
        )
    } else {
        _parse_smarts_params(text, params)
    }
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"parse_smarts"}"#]
#[pyo3(signature = (text, params=None))]
pub fn _parse_smarts_params(
    text: &str,
    params: Option<&SmartsParseParams>,
) -> PyResult<QueryGraph> {
    let result = match params {
        Some(params) => Python::attach(|py| {
            ck::parse_smarts_with_params(text, &params.to_rust(py)?)
                .map_err(SmartsParseError::from_rust)
        })
        .map(QueryGraph::from_rust),
        None => ck::parse_smarts(text)
            .map(QueryGraph::from_rust)
            .map_err(SmartsParseError::from_rust),
    };
    result
}

#[pyfunction]
#[doc = r#"{"role":"overload","target":"parse_smarts"}"#]
#[pyo3(signature = (text, *, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, merge_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
pub fn _parse_smarts_keywords(
    text: &str,
    allow_cxsmiles: Option<bool>,
    strict_cxsmiles: Option<bool>,
    parse_name: Option<bool>,
    merge_hs: Option<bool>,
    skip_cleanup: Option<bool>,
    debug_parse: Option<bool>,
    replacements: Option<BTreeMap<String, String>>,
) -> PyResult<QueryGraph> {
    let mut inner = ck::SmartsParseParams::default();
    if let Some(value) = allow_cxsmiles {
        inner.allow_cxsmiles = value;
    }
    if let Some(value) = strict_cxsmiles {
        inner.strict_cxsmiles = value;
    }
    if let Some(value) = parse_name {
        inner.parse_name = value;
    }
    if let Some(value) = merge_hs {
        inner.merge_hs = value;
    }
    if let Some(value) = skip_cleanup {
        inner.skip_cleanup = value;
    }
    if let Some(value) = debug_parse {
        inner.debug_parse = value;
    }
    if let Some(value) = replacements {
        inner.replacements = value
            .into_iter()
            .map(|(key, value)| (key.into(), value.into()))
            .collect();
    }
    ck::parse_smarts_with_params(text, &inner)
        .map(QueryGraph::from_rust)
        .map_err(SmartsParseError::from_rust)
}

#[pyfunction]
/// Parse SMARTS into a QueryGraph with explicit SmartsParseParams.
pub fn parse_smarts_with_params(text: &str, params: &SmartsParseParams) -> PyResult<QueryGraph> {
    _parse_smarts_params(text, Some(params))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::property_numeric_errors::PropertyUIntReadError;
    use crate::smarts_errors::{CxQueryLoweringError, SmartsParseErrorKind};
    use pyo3::types::{PyDict, PyModule};
    use std::error::Error;

    #[test]
    fn actual_smarts_factories_keep_one_parser_and_typed_context() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = PyModule::new(py, "smarts_factory_test")?;
            module.add_class::<QueryGraph>()?;
            module.add_class::<SmartsParseParams>()?;
            module.add_function(wrap_pyfunction!(parse_smarts, &module)?)?;
            let text = "[C][H]";
            let inner = ck::SmartsParseParams {
                merge_hs: true,
                ..Default::default()
            };
            let expected = ck::parse_smarts_with_params(text, &inner).unwrap();
            let params = Py::new(py, SmartsParseParams::from_rust(py, inner)?)?;
            let kwargs = PyDict::new(py);
            kwargs.set_item("merge_hs", true)?;
            let from_keywords = module
                .getattr("parse_smarts")?
                .call((text,), Some(&kwargs))?;
            assert_eq!(
                from_keywords.extract::<PyRef<'_, QueryGraph>>()?.inner(),
                &expected
            );
            let from_params = module
                .getattr("parse_smarts")?
                .call1((text, params.clone_ref(py)))?;
            assert_eq!(
                from_params.extract::<PyRef<'_, QueryGraph>>()?.inner(),
                &expected
            );
            let class = module.getattr("QueryGraph")?;
            let from_class = class.call_method("from_smarts", (text,), Some(&kwargs))?;
            assert_eq!(
                from_class.extract::<PyRef<'_, QueryGraph>>()?.inner(),
                &expected
            );
            let default = module.getattr("parse_smarts")?.call1(("CO",))?;
            assert_eq!(
                default.extract::<PyRef<'_, QueryGraph>>()?.inner(),
                &ck::parse_smarts("CO").unwrap()
            );
            assert!(
                module
                    .getattr("parse_smarts")?
                    .call((text, params), Some(&kwargs))
                    .unwrap_err()
                    .is_instance_of::<PyTypeError>(py)
            );

            let expected = ck::parse_smarts("[C").unwrap_err();
            assert!(expected.source().is_none());
            let error = module.getattr("parse_smarts")?.call1(("[C",)).unwrap_err();
            assert!(error.is_instance(py, <SmartsParseError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            assert!(error.cause(py).is_none());

            let error = SmartsParseError::from_rust(ck::SmartsParseError::UnclosedBracket(7));
            assert!(
                error
                    .value(py)
                    .getattr("variant")?
                    .extract::<SmartsParseErrorKind>()?
                    == SmartsParseErrorKind::UnclosedBracket
            );
            assert_eq!(error.value(py).getattr("position")?.extract::<usize>()?, 7);
            assert!(error.value(py).setattr("position", 0).is_err());
            let native = ck::SmartsParseError::CxLowering(ck::CxQueryLoweringError::Numeric(
                ck::PropertyUIntReadError::Negative { value: -7 },
            ));
            assert!(native.source().unwrap().is::<ck::CxQueryLoweringError>());
            let error = SmartsParseError::from_rust(native);
            let source = error.value(py).getattr("source")?;
            assert!(source.is(&error.value(py).getattr("__cause__")?));
            assert!(source.is_instance(<CxQueryLoweringError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert!(
                source
                    .getattr("source")?
                    .is_instance(<PropertyUIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?
            );
            assert_eq!(
                source
                    .getattr("source")?
                    .getattr("value")?
                    .extract::<i32>()?,
                -7
            );
            Ok(())
        })
        .unwrap();
    }
}
