//! Native reaction validation and initialization failures.
use crate::error_framework::{Source, error_projection};
use crate::model_edit_errors::AtomPropertyError;
use crate::property_numeric_errors::PropertyIntReadError;
use crate::reaction_values::{ReactionRole, ReactionValidationReport};
use cosmolkit as ck;
use pyo3::prelude::*;
use std::error::Error;

#[cfg(test)]
mod tests {
    use super::*;
    use pyo3::ffi::c_str;

    #[test]
    fn reaction_validation_and_initialization_forward_native_reports_and_value_state() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let module = pyo3::types::PyModule::new(py, "reaction_validation_test")?;
            module.add_class::<crate::reaction::Reaction>()?;
            module.add_class::<crate::query_compiler::QueryGraph>()?;
            module.add_class::<crate::params::ReactionValidationParams>()?;
            module.add_class::<ReactionValidationError>()?;
            module.add_class::<ReactionInitializationError>()?;
            let scope = pyo3::types::PyDict::new(py);
            scope.set_item("ck", &module)?;
            py.run(c_str!(r#"
query = ck.QueryGraph.from_smarts("[C:1]")
reaction = ck.Reaction.from_templates([query], [query], [])
params = ck.ReactionValidationParams(silent=True)
assert reaction.is_initialized() is False
for report in (reaction.validate(), reaction.validate(params), reaction.validate(silent=True), reaction.validate_with_params(params)):
    assert report.is_valid and report.num_errors == 0
for initialized in (reaction.with_initialized(), reaction.with_initialized(params), reaction.with_initialized(silent=True), reaction.with_initialized_with_params(params)):
    assert initialized.is_initialized() is True
    assert initialized.num_reactant_templates() == reaction.num_reactant_templates()
    assert reaction.is_initialized() is False
for method in (reaction.validate, reaction.with_initialized):
    try:
        method(params, silent=False)
    except TypeError:
        pass
    else:
        raise AssertionError("mixed configuration accepted")
empty = ck.Reaction()
report = empty.validate(silent=True)
assert report.is_valid is False
try:
    empty.with_initialized(silent=True)
except ck.ReactionInitializationError as error:
    assert error.report.num_errors == report.num_errors
    assert error.source is None and error.__cause__ is None
    try:
        error.report = None
    except AttributeError:
        pass
    else:
        raise AssertionError("invalid report must remain readonly")
else:
    raise AssertionError("invalid reaction initialization succeeded")
assert empty.is_initialized() is False
"#), Some(&scope), None)?;

            let mut atom = ck::QueryAtom::new(ck::AtomId::new(0), ck::AtomSpec::new(ck::Element::C));
            atom.set_prop("molAtomMapNumber", ck::PropertyValue::Bool(true)).unwrap();
            let graph = ck::QueryGraph::from_parts(vec![atom], vec![], [], vec![], vec![], vec![]).unwrap();
            let native = ck::Reaction::from_templates(vec![graph.clone()], vec![graph], vec![]).unwrap();
            let expected = native.validate_with_params(&ck::ReactionValidationParams { silent: true }).unwrap_err();
            let message = expected.to_string();
            let reaction = Py::new(py, crate::reaction::Reaction::from_rust(native))?;
            let keywords = pyo3::types::PyDict::new(py);
            keywords.set_item("silent", true)?;
            let error = reaction.bind(py).call_method("validate", (), Some(&keywords)).unwrap_err();
            assert!(error.is_instance(py, <ReactionValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert_eq!(error.value(py).str()?.to_str()?, message);
            let cause = error.value(py).getattr("source")?;
            assert!(cause.is(&error.value(py).getattr("__cause__")?));
            assert!(cause.is_instance(<PropertyIntReadError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any())?);
            assert_eq!(error.value(py).getattr("atom")?.extract::<usize>()?, 0);
            assert!(error.value(py).setattr("property", "other").is_err());
            let error = reaction.bind(py).call_method("with_initialized", (), Some(&keywords)).unwrap_err();
            assert!(error.is_instance(py, <ReactionInitializationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(error.cause(py).unwrap().is_instance(py, <ReactionValidationError as crate::error_framework::ExceptionPayload>::exception_type(py).expect("error backend").as_any()));
            assert!(error.value(py).getattr("report")?.is_none());
            Ok(())
        }).unwrap();
    }
}

error_projection!(ReactionValidationError, ReactionValidationErrorKind, ck::ReactionValidationError, "reaction", {
    Property { role, template, atom, property, source },
    MapOverflow { role, template, atom, property, value: u32 },
    AtomBounds { role, template, atom, atom_count: usize },
    MissingReactingAtom { template, atom, map: i32 },
    QueryArithmetic { template, atom },
    Annotation { template, atom, property, source },
}, cause = |inner: &ck::ReactionValidationError| match inner {
    ck::ReactionValidationError::Property { source, .. } => Some(PropertyIntReadError::from_rust(source.clone())),
    ck::ReactionValidationError::Annotation { source, .. } => Some(AtomPropertyError::from_rust(source.clone())),
    ck::ReactionValidationError::MapOverflow { .. } | ck::ReactionValidationError::AtomBounds { .. }
    | ck::ReactionValidationError::MissingReactingAtom { .. } | ck::ReactionValidationError::QueryArithmetic { .. } => None,
}, source = Option<Source<crate::error_framework::ErrorObject<PropertyIntReadError>, crate::error_framework::ErrorObject<AtomPropertyError>>>, methods = {
    #[getter]
    fn role(&self) -> Option<ReactionRole> { match self.inner {
        ck::ReactionValidationError::Property { role, .. }
        | ck::ReactionValidationError::MapOverflow { role, .. }
        | ck::ReactionValidationError::AtomBounds { role, .. } => Some(role.into()), _ => None,
    } }
    #[getter]
    fn template(&self) -> usize { match self.inner {
        ck::ReactionValidationError::Property { template, .. }
        | ck::ReactionValidationError::MapOverflow { template, .. }
        | ck::ReactionValidationError::AtomBounds { template, .. }
        | ck::ReactionValidationError::MissingReactingAtom { template, .. }
        | ck::ReactionValidationError::QueryArithmetic { template, .. }
        | ck::ReactionValidationError::Annotation { template, .. } => template,
    } }
    #[getter]
    fn atom(&self) -> usize { match self.inner {
        ck::ReactionValidationError::Property { atom, .. }
        | ck::ReactionValidationError::MapOverflow { atom, .. }
        | ck::ReactionValidationError::AtomBounds { atom, .. }
        | ck::ReactionValidationError::MissingReactingAtom { atom, .. }
        | ck::ReactionValidationError::QueryArithmetic { atom, .. }
        | ck::ReactionValidationError::Annotation { atom, .. } => atom.index(),
    } }
    #[getter]
    fn property(&self) -> Option<&'static str> { match self.inner {
        ck::ReactionValidationError::Property { property, .. }
        | ck::ReactionValidationError::MapOverflow { property, .. }
        | ck::ReactionValidationError::Annotation { property, .. } => Some(property), _ => None,
    } }
});

error_projection!(ReactionInitializationError, ReactionInitializationErrorKind, ck::ReactionInitializationError, "reaction", {
    Validation { source }, Invalid { report },
}, cause = |inner: &ck::ReactionInitializationError| match inner {
    ck::ReactionInitializationError::Validation { source } => Some(ReactionValidationError::from_rust(source.clone())),
    ck::ReactionInitializationError::Invalid { .. } => None,
}, source = Option<crate::error_framework::ErrorObject<ReactionValidationError>>, methods = {
    #[getter]
    #[doc = r#"{"role":"error_payload"}"#]
    fn report(&self) -> Option<ReactionValidationReport> { match &self.inner {
        ck::ReactionInitializationError::Invalid { report } => Some(ReactionValidationReport { inner: report.clone() }), _ => None,
    } }
});
