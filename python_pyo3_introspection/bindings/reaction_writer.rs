//! Reaction serialization calls and actual retained native errors.
use crate::error_framework::{Source, error_projection};
use crate::path_errors::PathError;
use crate::reaction_values::ReactionRole;
use crate::retained_error::RetainedError;
use crate::smarts_writer::SmartsWriteError;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::error::Error;
error_projection!(retained ReactionWriteError,ReactionWriteErrorKind,ck::ReactionWriteError,"reaction",{
Cx(_),Template{role,template,source},Connectivity{role,template,source},CoordinateSelectionArity{expected:usize,actual:usize},
},cause=|inner:&RetainedError<ck::ReactionWriteError>|match inner.as_ref(){
ck::ReactionWriteError::Cx(_)|ck::ReactionWriteError::Template{..}=>Some(SmartsWriteError::from_ref(inner.node().source().expect("actual SMARTS writer source"))),
ck::ReactionWriteError::Connectivity{source,..}=>Some(PathError::from_rust(source.clone())),
ck::ReactionWriteError::CoordinateSelectionArity{..}=>None,
},source=Option<Source<Py<SmartsWriteError>,Py<PathError>>>,methods={
#[getter]fn role(&self)->Option<ReactionRole>{match self.inner.as_ref(){ck::ReactionWriteError::Template{role,..}|ck::ReactionWriteError::Connectivity{role,..}=>Some((*role).into()),_=>None}}
#[getter]fn template(&self)->Option<usize>{match self.inner.as_ref(){ck::ReactionWriteError::Template{template,..}|ck::ReactionWriteError::Connectivity{template,..}=>Some(*template),_=>None}}
});
pub(crate) fn keywords(
    canonical: Option<bool>,
    isomeric_smiles: Option<bool>,
    rooted_at_atom: Option<usize>,
    include_dative_bonds: Option<bool>,
    include_cx: Option<bool>,
    cx_fields: Option<crate::notation_params::CxSmilesFieldsInput>,
    coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
) -> PyResult<crate::search_reaction_params::ReactionWriteParams> {
    let defaults = ck::ReactionWriteParams::default();
    crate::search_reaction_params::ReactionWriteParams::new(
        canonical.unwrap_or(defaults.canonical),
        isomeric_smiles.unwrap_or(defaults.isomeric_smiles),
        rooted_at_atom,
        include_dative_bonds.unwrap_or(defaults.include_dative_bonds),
        include_cx.unwrap_or(defaults.include_cx),
        cx_fields,
        coordinate_selections,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reaction::Reaction;
    use pyo3::types::PyDict;
    #[test]
    fn reaction_writer_family_matches_native_options_and_preserves_actual_rich_sources() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let source = ck::Reaction::from_smirks("[C:1].[O:2]>[N]>[C:1][O:2]").unwrap();
            let templates = (
                source.reactant_templates().to_vec(),
                source.agent_templates().to_vec(),
                source.product_templates().to_vec(),
            );
            let initialized = source.is_initialized();
            let actual = Py::new(py, Reaction::from_rust(source.clone()))?;
            let config = ck::ReactionWriteParams {
                canonical: true,
                ..Default::default()
            };
            let params = Py::new(
                py,
                crate::search_reaction_params::ReactionWriteParams {
                    inner: config.clone(),
                },
            )?;
            let kw = PyDict::new(py);
            kw.set_item("canonical", true)?;
            for (name, alias, default, configured) in [
                (
                    "to_smirks",
                    "to_smirks_with_params",
                    source.to_smirks().unwrap(),
                    source.to_smirks_with_params(&config).unwrap(),
                ),
                (
                    "to_cx_smirks",
                    "to_cx_smirks_with_params",
                    source.to_cx_smirks().unwrap(),
                    source.to_cx_smirks_with_params(&config).unwrap(),
                ),
            ] {
                for form in 0..4 {
                    let result = match form {
                        0 => actual.bind(py).call_method0(name)?,
                        1 => actual
                            .bind(py)
                            .call_method1(name, (params.clone_ref(py),))?,
                        2 => actual.bind(py).call_method(name, (), Some(&kw))?,
                        _ => actual
                            .bind(py)
                            .call_method1(alias, (params.clone_ref(py),))?,
                    }
                    .extract::<String>()?;
                    assert_eq!(
                        result.as_bytes(),
                        if form == 0 {
                            default.as_bytes()
                        } else {
                            configured.as_bytes()
                        }
                    );
                }
                let error = actual
                    .bind(py)
                    .call_method(name, (params.clone_ref(py),), Some(&kw))
                    .unwrap_err();
                assert!(error.is_instance_of::<pyo3::exceptions::PyTypeError>(py));
            }
            assert_eq!(actual.borrow(py).inner().is_initialized(), initialized);
            assert_eq!(actual.borrow(py).inner().reactant_templates(), templates.0);
            assert_eq!(actual.borrow(py).inner().agent_templates(), templates.1);
            assert_eq!(actual.borrow(py).inner().product_templates(), templates.2);
            let bad = ck::ReactionWriteParams {
                rooted_at_atom: Some(usize::MAX),
                ..Default::default()
            };
            let expected = source.to_smirks_with_params(&bad).unwrap_err();
            let kw = PyDict::new(py);
            kw.set_item("rooted_at_atom", usize::MAX)?;
            let error = actual
                .bind(py)
                .call_method("to_smirks", (), Some(&kw))
                .unwrap_err();
            assert!(error.is_instance_of::<ReactionWriteError>(py));
            assert_eq!(error.value(py).str()?.to_str()?, expected.to_string());
            let ck::ReactionWriteError::Template { role, template, .. } = expected else {
                panic!("native template error expected")
            };
            assert!(error.value(py).getattr("role")?.extract::<ReactionRole>()? == role.into());
            assert_eq!(
                error.value(py).getattr("template")?.extract::<usize>()?,
                template
            );
            let child = error.cause(py).unwrap();
            assert!(child.is_instance_of::<SmartsWriteError>(py));
            assert!(error.value(py).getattr("source")?.is(child.value(py)));
            assert!(error.value(py).setattr("template", 0).is_err());
            assert_eq!(
                child.value(py).getattr("atom")?.extract::<usize>()?,
                usize::MAX
            );
            let rich = ReactionWriteError::from_rust(ck::ReactionWriteError::Cx(
                ck::SmartsWriteError::UnsupportedAtomQuery {
                    predicate: ck::AtomQueryPredicate::AtomicNumberNotIn(vec![6, 8]),
                },
            ));
            let child = rich.cause(py).unwrap();
            let payload = child.value(py).getattr("predicate")?.unbind();
            drop(rich);
            drop(child);
            drop(actual);
            assert_eq!(
                payload.bind(py).getattr("value")?.extract::<Vec<u8>>()?,
                vec![6, 8]
            );
            assert!(payload.bind(py).setattr("value", vec![7]).is_err());
            let connectivity =
                ReactionWriteError::from_rust(ck::ReactionWriteError::Connectivity {
                    role: ck::ReactionRole::Product,
                    template: usize::MAX,
                    source: ck::PathError::InvalidLengthRange { lower: 3, upper: 2 },
                });
            assert!(
                connectivity
                    .cause(py)
                    .unwrap()
                    .is_instance_of::<PathError>(py)
            );
            assert_eq!(
                connectivity
                    .value(py)
                    .getattr("template")?
                    .extract::<usize>()?,
                usize::MAX
            );
            let arity =
                ReactionWriteError::from_rust(ck::ReactionWriteError::CoordinateSelectionArity {
                    expected: 4,
                    actual: 1,
                });
            assert!(arity.cause(py).is_none());
            assert_eq!(arity.value(py).getattr("expected")?.extract::<usize>()?, 4);
            assert_eq!(arity.value(py).getattr("actual")?.extract::<usize>()?, 1);
            Ok(())
        })
        .unwrap();
    }
}
