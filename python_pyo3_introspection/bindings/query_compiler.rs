//! Immutable native query transport and the single facade compiler.
use crate::notation_params::SmartsParseParams;
use crate::string_errors::QueryCompileError;
use cosmolkit as ck;
use pyo3::prelude::*;
use std::sync::Arc;

#[derive(Clone)]
enum QueryOwner {
    Graph(Arc<ck::QueryGraph>),
    Compiled(Arc<ck::CompiledQuery>),
    McsPrimary(Arc<ck::McsResult>),
    McsDegenerate(Arc<ck::McsResult>, ck::PropertyText),
    ReactionTemplate(crate::reaction::ReactionOwner, ck::ReactionRole, usize),
    RemovedTemplate(Arc<ck::ReactionTemplateRemoval>, usize),
    WriterRecursive(crate::retained_error::RetainedError<ck::SmartsWriteError>),
    SdfRecord(Arc<ck::SdfRecord>),
}

#[pyclass(module = "fixed_projection._native", frozen, from_py_object)]
#[derive(Clone)]
pub struct QueryGraph {
    owner: QueryOwner,
}

impl QueryGraph {
    pub(crate) fn from_rust(inner: ck::QueryGraph) -> Self {
        Self {
            owner: QueryOwner::Graph(Arc::new(inner)),
        }
    }
    pub(crate) fn from_sdf_record(owner: Arc<ck::SdfRecord>) -> Self {
        Self {
            owner: QueryOwner::SdfRecord(owner),
        }
    }
    pub(crate) fn inner(&self) -> &ck::QueryGraph {
        match &self.owner {
            QueryOwner::Graph(owner) => owner,
            QueryOwner::Compiled(owner) => owner.query(),
            QueryOwner::McsPrimary(owner) => owner.query.as_ref().expect("native MCS query view"),
            QueryOwner::McsDegenerate(owner, key) => &owner.degenerate[key],
            QueryOwner::ReactionTemplate(owner, role, index) => match role {
                ck::ReactionRole::Reactant => &owner.inner().reactant_templates()[*index],
                ck::ReactionRole::Product => &owner.inner().product_templates()[*index],
                ck::ReactionRole::Agent => &owner.inner().agent_templates()[*index],
            },
            QueryOwner::RemovedTemplate(owner, index) => &owner.removed_templates()[*index],
            QueryOwner::SdfRecord(owner) => owner.query_graph().expect("native SDF query view"),
            QueryOwner::WriterRecursive(owner) => {
                let ck::SmartsWriteError::UnsupportedAtomQuery {
                    predicate: ck::AtomQueryPredicate::RecursiveSmarts(query),
                } = owner.as_ref()
                else {
                    unreachable!("native recursive writer payload")
                };
                query.query_graph().expect("native recursive graph view")
            }
        }
    }
    pub(crate) fn from_mcs_query(owner: Arc<ck::McsResult>) -> Option<Self> {
        owner.query.as_ref()?;
        Some(Self {
            owner: QueryOwner::McsPrimary(owner),
        })
    }
    pub(crate) fn from_mcs_degenerate(owner: Arc<ck::McsResult>, key: ck::PropertyText) -> Self {
        Self {
            owner: QueryOwner::McsDegenerate(owner, key),
        }
    }
    pub(crate) fn from_reaction_template(
        owner: crate::reaction::ReactionOwner,
        role: ck::ReactionRole,
        index: usize,
    ) -> Self {
        Self {
            owner: QueryOwner::ReactionTemplate(owner, role, index),
        }
    }
    pub(crate) fn from_removed_template(
        owner: Arc<ck::ReactionTemplateRemoval>,
        index: usize,
    ) -> Self {
        Self {
            owner: QueryOwner::RemovedTemplate(owner, index),
        }
    }
    pub(crate) fn from_writer_recursive(
        owner: crate::retained_error::RetainedError<ck::SmartsWriteError>,
    ) -> Self {
        Self {
            owner: QueryOwner::WriterRecursive(owner),
        }
    }
}

#[pymethods]
impl QueryGraph {
    #[staticmethod]
    #[pyo3(signature = (text, params=None, *, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, merge_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
    fn from_smarts(
        text: &str,
        params: Option<&SmartsParseParams>,
        allow_cxsmiles: Option<bool>,
        strict_cxsmiles: Option<bool>,
        parse_name: Option<bool>,
        merge_hs: Option<bool>,
        skip_cleanup: Option<bool>,
        debug_parse: Option<bool>,
        replacements: Option<std::collections::BTreeMap<String, String>>,
    ) -> PyResult<Self> {
        let keywords = [
            allow_cxsmiles,
            strict_cxsmiles,
            parse_name,
            merge_hs,
            skip_cleanup,
            debug_parse,
        ]
        .iter()
        .any(Option::is_some)
            || replacements.is_some();
        if params.is_some() && keywords {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and keyword configuration are mutually exclusive",
            ));
        }
        if keywords {
            Self::_from_smarts_keywords(
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
            Self::_from_smarts_params(text, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"QueryGraph.from_smarts"}"#]
    #[pyo3(signature = (text, params=None))]
    fn _from_smarts_params(text: &str, params: Option<&SmartsParseParams>) -> PyResult<Self> {
        crate::smarts_factory::_parse_smarts_params(text, params)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"QueryGraph.from_smarts"}"#]
    #[pyo3(signature = (text, *, allow_cxsmiles=None, strict_cxsmiles=None, parse_name=None, merge_hs=None, skip_cleanup=None, debug_parse=None, replacements=None))]
    fn _from_smarts_keywords(
        text: &str,
        allow_cxsmiles: Option<bool>,
        strict_cxsmiles: Option<bool>,
        parse_name: Option<bool>,
        merge_hs: Option<bool>,
        skip_cleanup: Option<bool>,
        debug_parse: Option<bool>,
        replacements: Option<std::collections::BTreeMap<String, String>>,
    ) -> PyResult<Self> {
        crate::smarts_factory::_parse_smarts_keywords(
            text,
            allow_cxsmiles,
            strict_cxsmiles,
            parse_name,
            merge_hs,
            skip_cleanup,
            debug_parse,
            replacements,
        )
    }
    #[staticmethod]
    fn from_smarts_with_params(text: &str, params: &SmartsParseParams) -> PyResult<Self> {
        Self::_from_smarts_params(text, Some(params))
    }
    fn num_atoms(&self) -> usize {
        self.inner().num_atoms()
    }
    fn num_bonds(&self) -> usize {
        self.inner().num_bonds()
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct CompiledQuery {
    inner: Arc<ck::CompiledQuery>,
}

impl CompiledQuery {
    pub(crate) fn inner(&self) -> &ck::CompiledQuery {
        &self.inner
    }
}

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
            owner: QueryOwner::Compiled(Arc::clone(&self.inner)),
        }
    }
}

#[pyfunction]
pub fn compile_query(query: &QueryGraph) -> PyResult<CompiledQuery> {
    ck::compile_query(query.inner())
        .map(|inner| CompiledQuery {
            inner: Arc::new(inner),
        })
        .map_err(QueryCompileError::from_rust)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn actual_compiler_preserves_native_graph_order_and_owned_view() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let graph = ck::parse_smarts("C(=O)O").unwrap();
            let expected = ck::compile_query(&graph).unwrap();
            let query = QueryGraph::from_rust(graph);
            let compiled = compile_query(&query)?;
            assert_eq!(compiled.num_atoms(), expected.num_atoms());
            assert_eq!(compiled.num_bonds(), expected.num_bonds());
            assert_eq!(compiled.atom_order(), expected.atom_order());
            let view = compiled.query();
            let QueryOwner::Compiled(owner) = &view.owner else {
                panic!("expected native owner view")
            };
            assert!(Arc::ptr_eq(owner, &compiled.inner));
            assert!(std::ptr::eq(view.inner(), compiled.inner.query()));
            drop(query);
            drop(compiled);
            let view = Py::new(py, view)?;
            assert_eq!(
                view.bind(py)
                    .call_method0("num_atoms")?
                    .extract::<usize>()?,
                expected.num_atoms()
            );
            assert_eq!(
                view.bind(py)
                    .call_method0("num_bonds")?
                    .extract::<usize>()?,
                expected.num_bonds()
            );
            let module = PyModule::new(py, "compiler_probe")?;
            module.add_function(wrap_pyfunction!(compile_query, &module)?)?;
            let compiled = module.getattr("compile_query")?.call1((view,))?;
            assert_eq!(
                compiled
                    .call_method0("atom_order")?
                    .extract::<Vec<usize>>()?,
                expected.atom_order()
            );
            Ok(())
        })
        .unwrap();
    }
}
