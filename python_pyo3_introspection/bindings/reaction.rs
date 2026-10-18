//! Reaction template ownership and model operations, forwarded to the facade.
use crate::error_framework::error_projection;
use crate::params::ReactionParseParams;
use crate::params::{ReactionTemplateRemovalParams, ReactionValidationParams};
use crate::query_compiler::QueryGraph;
use crate::query_model_errors::QueryGraphError;
use crate::reaction_parse::ReactionParseError;
use crate::reaction_validation::{ReactionInitializationError, ReactionValidationError};
use crate::reaction_values::ReactionRole;
use crate::reaction_values::ReactionValidationReport;
use crate::search_callbacks::SubstructMatchParams;
use cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;
use std::error::Error;
use std::sync::Arc;
use std::sync::Mutex;

error_projection!(ReactionModelError, ReactionModelErrorKind, ck::ReactionModelError, "reaction", {
    Template { role, template: usize, source }, TemplateIndex { role, index: usize, count: usize },
}, cause = |inner: &ck::ReactionModelError| match inner {
    ck::ReactionModelError::Template { source, .. } => Some(QueryGraphError::from_rust(source.clone())),
    ck::ReactionModelError::TemplateIndex { .. } => None,
}, source = Option<Py<QueryGraphError>>, methods = {
    #[getter]
    fn role(&self) -> ReactionRole {
        match self.inner { ck::ReactionModelError::Template { role, .. } | ck::ReactionModelError::TemplateIndex { role, .. } => role.into() }
    }
});

#[pyclass(module = "fixed_projection._native")]
pub struct Reaction {
    owner: ReactionOwner,
    match_owner: Mutex<Option<Py<SubstructMatchParams>>>,
}
#[derive(Clone)]
pub(crate) enum ReactionOwner {
    Value(Arc<ck::Reaction>),
    Removal(Arc<ck::ReactionTemplateRemoval>),
}
impl ReactionOwner {
    pub(crate) fn inner(&self) -> &ck::Reaction {
        match self {
            Self::Value(inner) => inner,
            Self::Removal(inner) => inner.reaction(),
        }
    }

    pub(crate) fn inner_mut(&mut self) -> &mut ck::Reaction {
        if let Self::Removal(owner) = self {
            if Arc::get_mut(owner).is_none() {
                *self = Self::Value(Arc::new(owner.reaction().clone()));
            }
        }
        match self {
            Self::Value(owner) => Arc::make_mut(owner),
            Self::Removal(owner) => {
                &mut Arc::get_mut(owner).expect("unique removal owner").reaction
            }
        }
    }
}
impl Reaction {
    pub(crate) fn from_rust(inner: ck::Reaction) -> Self {
        Self {
            owner: ReactionOwner::Value(Arc::new(inner)),
            match_owner: Mutex::new(None),
        }
    }
    pub(crate) fn inner(&self) -> &ck::Reaction {
        self.owner.inner()
    }
    pub(crate) fn native_mut(&mut self) -> &mut ck::Reaction {
        self.owner.inner_mut()
    }
    pub(crate) fn copy_match_owner(&self) -> Option<Py<SubstructMatchParams>> {
        Python::attach(|py| {
            self.match_owner
                .lock()
                .expect("reaction configuration lock")
                .as_ref()
                .map(|value| value.clone_ref(py))
        })
    }
    fn with_native(&self, inner: ck::Reaction) -> Self {
        Self {
            owner: ReactionOwner::Value(Arc::new(inner)),
            match_owner: Mutex::new(self.copy_match_owner()),
        }
    }
    fn with_removal(&self, inner: ck::ReactionTemplateRemoval) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: Arc::new(inner),
            match_owner: Mutex::new(self.copy_match_owner()),
        }
    }
    fn template_view(&self, role: ck::ReactionRole, index: usize) -> QueryGraph {
        QueryGraph::from_reaction_template(self.owner.clone(), role, index)
    }
}

#[pymethods]
impl Reaction {
    #[pyo3(signature=(reactants,params=None,*,max_products=None,coordinate_selections=None,copy_atom_properties=None))]
    fn run(
        &mut self,
        py: Python<'_>,
        reactants: Vec<Py<crate::native::Molecule>>,
        params: Option<&crate::optional_configuration::ReactionRunParams>,
        max_products: Option<u32>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
        copy_atom_properties: Option<bool>,
    ) -> PyResult<Vec<Vec<crate::native::Molecule>>> {
        if params.is_some()
            && (max_products.is_some()
                || coordinate_selections.is_some()
                || copy_atom_properties.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "params and keyword configuration are mutually exclusive",
            ));
        }
        match params {
            Some(params) => self._run_params(py, reactants, Some(params)),
            None => self._run_keywords(
                py,
                reactants,
                max_products,
                coordinate_selections,
                copy_atom_properties,
            ),
        }
    }
    #[doc = r#"{"role":"overload","target":"Reaction.run"}"#]
    #[pyo3(signature=(reactants,params=None))]
    fn _run_params(
        &mut self,
        py: Python<'_>,
        reactants: Vec<Py<crate::native::Molecule>>,
        params: Option<&crate::optional_configuration::ReactionRunParams>,
    ) -> PyResult<Vec<Vec<crate::native::Molecule>>> {
        let defaults;
        let params = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::ReactionRunParams::default();
                &defaults
            }
        };
        crate::reaction_calls::run(self, py, &reactants, params)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.run"}"#]
    #[pyo3(signature=(reactants,*,max_products=None,coordinate_selections=None,copy_atom_properties=None))]
    fn _run_keywords(
        &mut self,
        py: Python<'_>,
        reactants: Vec<Py<crate::native::Molecule>>,
        max_products: Option<u32>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
        copy_atom_properties: Option<bool>,
    ) -> PyResult<Vec<Vec<crate::native::Molecule>>> {
        let params = crate::reaction_calls::keywords(
            max_products,
            coordinate_selections,
            copy_atom_properties,
        );
        self._run_params(py, reactants, Some(&params))
    }
    fn set_match_params(&mut self, py: Python<'_>, params: &SubstructMatchParams) -> PyResult<()> {
        let owner = Py::new(py, params.with_native(py, params.inner.clone()))?;
        self.owner.inner_mut().set_match_params(&params.inner);
        *self
            .match_owner
            .lock()
            .expect("reaction configuration lock") = Some(owner);
        Ok(())
    }
    #[pyo3(signature=(params=None,*,canonical=None,isomeric_smiles=None,rooted_at_atom=None,include_dative_bonds=None,include_cx=None,cx_fields=None,coordinate_selections=None))]
    fn to_smirks(
        &self,
        py: Python<'_>,
        params: Option<&crate::search_reaction_params::ReactionWriteParams>,
        canonical: Option<bool>,
        isomeric_smiles: Option<bool>,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: Option<bool>,
        include_cx: Option<bool>,
        cx_fields: Option<crate::notation_params::CxSmilesFieldsInput>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
    ) -> PyResult<String> {
        if params.is_some()
            && (canonical.is_some()
                || isomeric_smiles.is_some()
                || rooted_at_atom.is_some()
                || include_dative_bonds.is_some()
                || include_cx.is_some()
                || cx_fields.is_some()
                || coordinate_selections.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "params and keyword configuration are mutually exclusive",
            ));
        }
        if params.is_some()
            || (canonical.is_none()
                && isomeric_smiles.is_none()
                && rooted_at_atom.is_none()
                && include_dative_bonds.is_none()
                && include_cx.is_none()
                && cx_fields.is_none()
                && coordinate_selections.is_none())
        {
            self._writer_to_smirks_params(py, params)
        } else {
            self._writer_to_smirks_keywords(
                py,
                canonical,
                isomeric_smiles,
                rooted_at_atom,
                include_dative_bonds,
                include_cx,
                cx_fields,
                coordinate_selections,
            )
        }
    }
    fn to_smirks_with_params(
        &self,
        py: Python<'_>,
        params: &crate::search_reaction_params::ReactionWriteParams,
    ) -> PyResult<String> {
        let text = py
            .detach(|| self.inner().to_smirks_with_params(&params.inner))
            .map_err(crate::reaction_writer::ReactionWriteError::from_rust)?;
        crate::property_values::decode_source_text(py, &text)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.to_smirks"}"#]
    #[pyo3(signature=(params=None))]
    fn _writer_to_smirks_params(
        &self,
        py: Python<'_>,
        params: Option<&crate::search_reaction_params::ReactionWriteParams>,
    ) -> PyResult<String> {
        let text = match params {
            Some(params) => py.detach(|| self.inner().to_smirks_with_params(&params.inner)),
            None => py.detach(|| self.inner().to_smirks()),
        }
        .map_err(crate::reaction_writer::ReactionWriteError::from_rust)?;
        crate::property_values::decode_source_text(py, &text)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.to_smirks"}"#]
    #[pyo3(signature=(*,canonical=None,isomeric_smiles=None,rooted_at_atom=None,include_dative_bonds=None,include_cx=None,cx_fields=None,coordinate_selections=None))]
    fn _writer_to_smirks_keywords(
        &self,
        py: Python<'_>,
        canonical: Option<bool>,
        isomeric_smiles: Option<bool>,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: Option<bool>,
        include_cx: Option<bool>,
        cx_fields: Option<crate::notation_params::CxSmilesFieldsInput>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
    ) -> PyResult<String> {
        let params = crate::reaction_writer::keywords(
            canonical,
            isomeric_smiles,
            rooted_at_atom,
            include_dative_bonds,
            include_cx,
            cx_fields,
            coordinate_selections,
        )?;
        self._writer_to_smirks_params(py, Some(&params))
    }

    #[pyo3(signature=(params=None,*,canonical=None,isomeric_smiles=None,rooted_at_atom=None,include_dative_bonds=None,include_cx=None,cx_fields=None,coordinate_selections=None))]
    fn to_cx_smirks(
        &self,
        py: Python<'_>,
        params: Option<&crate::search_reaction_params::ReactionWriteParams>,
        canonical: Option<bool>,
        isomeric_smiles: Option<bool>,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: Option<bool>,
        include_cx: Option<bool>,
        cx_fields: Option<crate::notation_params::CxSmilesFieldsInput>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
    ) -> PyResult<String> {
        if params.is_some()
            && (canonical.is_some()
                || isomeric_smiles.is_some()
                || rooted_at_atom.is_some()
                || include_dative_bonds.is_some()
                || include_cx.is_some()
                || cx_fields.is_some()
                || coordinate_selections.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "params and keyword configuration are mutually exclusive",
            ));
        }
        if params.is_some()
            || (canonical.is_none()
                && isomeric_smiles.is_none()
                && rooted_at_atom.is_none()
                && include_dative_bonds.is_none()
                && include_cx.is_none()
                && cx_fields.is_none()
                && coordinate_selections.is_none())
        {
            self._writer_to_cx_smirks_params(py, params)
        } else {
            self._writer_to_cx_smirks_keywords(
                py,
                canonical,
                isomeric_smiles,
                rooted_at_atom,
                include_dative_bonds,
                include_cx,
                cx_fields,
                coordinate_selections,
            )
        }
    }
    fn to_cx_smirks_with_params(
        &self,
        py: Python<'_>,
        params: &crate::search_reaction_params::ReactionWriteParams,
    ) -> PyResult<String> {
        let text = py
            .detach(|| self.inner().to_cx_smirks_with_params(&params.inner))
            .map_err(crate::reaction_writer::ReactionWriteError::from_rust)?;
        crate::property_values::decode_source_text(py, &text)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.to_cx_smirks"}"#]
    #[pyo3(signature=(params=None))]
    fn _writer_to_cx_smirks_params(
        &self,
        py: Python<'_>,
        params: Option<&crate::search_reaction_params::ReactionWriteParams>,
    ) -> PyResult<String> {
        let text = match params {
            Some(params) => py.detach(|| self.inner().to_cx_smirks_with_params(&params.inner)),
            None => py.detach(|| self.inner().to_cx_smirks()),
        }
        .map_err(crate::reaction_writer::ReactionWriteError::from_rust)?;
        crate::property_values::decode_source_text(py, &text)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.to_cx_smirks"}"#]
    #[pyo3(signature=(*,canonical=None,isomeric_smiles=None,rooted_at_atom=None,include_dative_bonds=None,include_cx=None,cx_fields=None,coordinate_selections=None))]
    fn _writer_to_cx_smirks_keywords(
        &self,
        py: Python<'_>,
        canonical: Option<bool>,
        isomeric_smiles: Option<bool>,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: Option<bool>,
        include_cx: Option<bool>,
        cx_fields: Option<crate::notation_params::CxSmilesFieldsInput>,
        coordinate_selections: Option<Vec<crate::construction_values::ReactionCoordinateSelection>>,
    ) -> PyResult<String> {
        let params = crate::reaction_writer::keywords(
            canonical,
            isomeric_smiles,
            rooted_at_atom,
            include_dative_bonds,
            include_cx,
            cx_fields,
            coordinate_selections,
        )?;
        self._writer_to_cx_smirks_params(py, Some(&params))
    }
    fn match_params(&self, py: Python<'_>) -> SubstructMatchParams {
        let native = self.inner().match_params().clone();
        match self
            .match_owner
            .lock()
            .expect("reaction configuration lock")
            .as_ref()
        {
            Some(owner) => owner.borrow(py).with_native(py, native),
            None => SubstructMatchParams::from_native(native),
        }
    }
    fn with_match_params(&self, py: Python<'_>, params: &SubstructMatchParams) -> PyResult<Self> {
        let inner = self
            .inner()
            .with_match_params(&params.inner)
            .map_err(ReactionModelError::from_rust)?;
        let owner = Py::new(py, params.with_native(py, inner.match_params().clone()))?;
        Ok(Self {
            owner: ReactionOwner::Value(Arc::new(inner)),
            match_owner: Mutex::new(Some(owner)),
        })
    }
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(value) = self
            .match_owner
            .lock()
            .expect("reaction configuration lock")
            .as_ref()
        {
            visit.call(value)?;
        }
        Ok(())
    }
    fn __clear__(&self) {
        let previous = self
            .match_owner
            .lock()
            .expect("reaction configuration lock")
            .take();
        drop(previous);
    }
    #[staticmethod]
    #[pyo3(signature=(text, params=None, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
    fn from_smirks(
        text: &str,
        params: Option<&ReactionParseParams>,
        use_smiles: Option<bool>,
        sanitize: Option<bool>,
        replacements: Option<BTreeMap<String, String>>,
        allow_cxsmiles: Option<bool>,
        strict_cxsmiles: Option<bool>,
    ) -> PyResult<Self> {
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
            Self::_from_smirks_keywords(
                text,
                use_smiles,
                sanitize,
                replacements,
                allow_cxsmiles,
                strict_cxsmiles,
            )
        } else {
            Self::_from_smirks_params(text, params)
        }
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Reaction.from_smirks"}"#]
    #[pyo3(signature=(text, params=None))]
    fn _from_smirks_params(text: &str, params: Option<&ReactionParseParams>) -> PyResult<Self> {
        match params {
            Some(params) => ck::Reaction::from_smirks_with_params(text, &params.inner),
            None => ck::Reaction::from_smirks(text),
        }
        .map(Self::from_rust)
        .map_err(ReactionParseError::from_rust)
    }
    #[staticmethod]
    #[doc = r#"{"role":"overload","target":"Reaction.from_smirks"}"#]
    #[pyo3(signature=(text, *, use_smiles=None, sanitize=None, replacements=None, allow_cxsmiles=None, strict_cxsmiles=None))]
    fn _from_smirks_keywords(
        text: &str,
        use_smiles: Option<bool>,
        sanitize: Option<bool>,
        replacements: Option<BTreeMap<String, String>>,
        allow_cxsmiles: Option<bool>,
        strict_cxsmiles: Option<bool>,
    ) -> PyResult<Self> {
        let params = crate::reaction_parse::keywords(
            use_smiles,
            sanitize,
            replacements,
            allow_cxsmiles,
            strict_cxsmiles,
        );
        Self::_from_smirks_params(text, Some(&params))
    }
    #[staticmethod]
    fn from_smirks_with_params(text: &str, params: &ReactionParseParams) -> PyResult<Self> {
        ck::Reaction::from_smirks_with_params(text, &params.inner)
            .map(Self::from_rust)
            .map_err(ReactionParseError::from_rust)
    }
    #[pyo3(signature = (params=None, *, silent=None))]
    fn validate(
        &self,
        params: Option<&ReactionValidationParams>,
        silent: Option<bool>,
    ) -> PyResult<ReactionValidationReport> {
        if params.is_some() && silent.is_some() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and keyword configuration are mutually exclusive",
            ));
        }
        if silent.is_some() {
            self._validate_keywords(silent)
        } else {
            self._validate_params(params)
        }
    }
    #[doc = r#"{"role":"overload","target":"Reaction.validate"}"#]
    #[pyo3(signature = (params=None))]
    fn _validate_params(
        &self,
        params: Option<&ReactionValidationParams>,
    ) -> PyResult<ReactionValidationReport> {
        match params {
            Some(params) => self.inner().validate_with_params(&params.inner),
            None => self.inner().validate(),
        }
        .map(|inner| ReactionValidationReport { inner })
        .map_err(ReactionValidationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.validate"}"#]
    #[pyo3(signature = (*, silent=None))]
    fn _validate_keywords(&self, silent: Option<bool>) -> PyResult<ReactionValidationReport> {
        let params = ReactionValidationParams {
            inner: ck::ReactionValidationParams {
                silent: silent.unwrap_or_else(|| ck::ReactionValidationParams::default().silent),
            },
        };
        self._validate_params(Some(&params))
    }
    fn validate_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> PyResult<ReactionValidationReport> {
        self._validate_params(Some(params))
    }

    #[pyo3(signature = (params=None, *, silent=None))]
    fn with_initialized(
        &self,
        params: Option<&ReactionValidationParams>,
        silent: Option<bool>,
    ) -> PyResult<Self> {
        if params.is_some() && silent.is_some() {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Params and keyword configuration are mutually exclusive",
            ));
        }
        if silent.is_some() {
            self._with_initialized_keywords(silent)
        } else {
            self._with_initialized_params(params)
        }
    }
    #[doc = r#"{"role":"overload","target":"Reaction.with_initialized"}"#]
    #[pyo3(signature = (params=None))]
    fn _with_initialized_params(
        &self,
        params: Option<&ReactionValidationParams>,
    ) -> PyResult<Self> {
        match params {
            Some(params) => self.inner().with_initialized_with_params(&params.inner),
            None => self.inner().with_initialized(),
        }
        .map(|inner| self.with_native(inner))
        .map_err(ReactionInitializationError::from_rust)
    }
    #[doc = r#"{"role":"overload","target":"Reaction.with_initialized"}"#]
    #[pyo3(signature = (*, silent=None))]
    fn _with_initialized_keywords(&self, silent: Option<bool>) -> PyResult<Self> {
        let params = ReactionValidationParams {
            inner: ck::ReactionValidationParams {
                silent: silent.unwrap_or_else(|| ck::ReactionValidationParams::default().silent),
            },
        };
        self._with_initialized_params(Some(&params))
    }
    fn with_initialized_with_params(&self, params: &ReactionValidationParams) -> PyResult<Self> {
        self._with_initialized_params(Some(params))
    }

    #[new]
    fn new() -> Self {
        Self::from_rust(ck::Reaction::new())
    }
    #[staticmethod]
    fn from_templates(
        py: Python<'_>,
        reactants: Vec<Py<QueryGraph>>,
        products: Vec<Py<QueryGraph>>,
        agents: Vec<Py<QueryGraph>>,
    ) -> PyResult<Self> {
        let owned = |values: Vec<Py<QueryGraph>>| {
            values
                .into_iter()
                .map(|value| value.borrow(py).inner().clone())
                .collect()
        };
        ck::Reaction::from_templates(owned(reactants), owned(products), owned(agents))
            .map(Self::from_rust)
            .map_err(ReactionModelError::from_rust)
    }
    fn num_reactant_templates(&self) -> usize {
        self.inner().num_reactant_templates()
    }
    fn num_product_templates(&self) -> usize {
        self.inner().num_product_templates()
    }
    fn num_agent_templates(&self) -> usize {
        self.inner().num_agent_templates()
    }
    fn reactant_templates(&self) -> Vec<QueryGraph> {
        self.inner()
            .reactant_templates()
            .iter()
            .enumerate()
            .map(|(index, _)| self.template_view(ck::ReactionRole::Reactant, index))
            .collect()
    }
    fn product_templates(&self) -> Vec<QueryGraph> {
        self.inner()
            .product_templates()
            .iter()
            .enumerate()
            .map(|(index, _)| self.template_view(ck::ReactionRole::Product, index))
            .collect()
    }
    fn agent_templates(&self) -> Vec<QueryGraph> {
        self.inner()
            .agent_templates()
            .iter()
            .enumerate()
            .map(|(index, _)| self.template_view(ck::ReactionRole::Agent, index))
            .collect()
    }
    fn reactant_template(&self, index: usize) -> PyResult<QueryGraph> {
        self.inner()
            .reactant_template(index)
            .map(|_| self.template_view(ck::ReactionRole::Reactant, index))
            .map_err(ReactionModelError::from_rust)
    }
    fn product_template(&self, index: usize) -> PyResult<QueryGraph> {
        self.inner()
            .product_template(index)
            .map(|_| self.template_view(ck::ReactionRole::Product, index))
            .map_err(ReactionModelError::from_rust)
    }
    fn agent_template(&self, index: usize) -> PyResult<QueryGraph> {
        self.inner()
            .agent_template(index)
            .map(|_| self.template_view(ck::ReactionRole::Agent, index))
            .map_err(ReactionModelError::from_rust)
    }
    fn with_reactant_template(&self, template: &QueryGraph) -> PyResult<Self> {
        self.inner()
            .with_reactant_template(template.inner().clone())
            .map(|inner| self.with_native(inner))
            .map_err(ReactionModelError::from_rust)
    }
    fn with_product_template(&self, template: &QueryGraph) -> PyResult<Self> {
        self.inner()
            .with_product_template(template.inner().clone())
            .map(|inner| self.with_native(inner))
            .map_err(ReactionModelError::from_rust)
    }
    fn with_agent_template(&self, template: &QueryGraph) -> PyResult<Self> {
        self.inner()
            .with_agent_template(template.inner().clone())
            .map(|inner| self.with_native(inner))
            .map_err(ReactionModelError::from_rust)
    }
    fn is_initialized(&self) -> bool {
        self.inner().is_initialized()
    }
    fn implicit_properties(&self) -> bool {
        self.inner().implicit_properties()
    }
    fn with_implicit_properties(&self, enabled: bool) -> Self {
        self.with_native(self.inner().with_implicit_properties(enabled))
    }
    fn without_agents(&self) -> ReactionTemplateRemoval {
        self.with_removal(self.inner().without_agents())
    }

    #[pyo3(signature = (params=None, *, threshold_unmapped_atoms=None, move_to_agent_templates=None))]
    fn without_unmapped_reactants(
        &self,
        params: Option<&ReactionTemplateRemovalParams>,
        threshold_unmapped_atoms: Option<f64>,
        move_to_agent_templates: Option<bool>,
    ) -> PyResult<ReactionTemplateRemoval> {
        if params.is_some()
            && (threshold_unmapped_atoms.is_some() || move_to_agent_templates.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        Ok(if params.is_some() {
            self._without_unmapped_reactants_params(params)
        } else {
            self._without_unmapped_reactants_keywords(
                threshold_unmapped_atoms,
                move_to_agent_templates,
            )
        })
    }
    #[doc = r#"{"role":"overload","target":"Reaction.without_unmapped_reactants"}"#]
    #[pyo3(signature = (params=None))]
    fn _without_unmapped_reactants_params(
        &self,
        params: Option<&ReactionTemplateRemovalParams>,
    ) -> ReactionTemplateRemoval {
        self.with_removal(match params {
            Some(params) => self
                .inner()
                .without_unmapped_reactants_with_params(&params.inner),
            None => self.inner().without_unmapped_reactants(),
        })
    }
    #[doc = r#"{"role":"overload","target":"Reaction.without_unmapped_reactants"}"#]
    #[pyo3(signature = (*, threshold_unmapped_atoms=None, move_to_agent_templates=None))]
    fn _without_unmapped_reactants_keywords(
        &self,
        threshold_unmapped_atoms: Option<f64>,
        move_to_agent_templates: Option<bool>,
    ) -> ReactionTemplateRemoval {
        let defaults = ck::ReactionTemplateRemovalParams::default();
        let params = ReactionTemplateRemovalParams {
            inner: ck::ReactionTemplateRemovalParams {
                threshold_unmapped_atoms: threshold_unmapped_atoms
                    .unwrap_or(defaults.threshold_unmapped_atoms),
                move_to_agent_templates: move_to_agent_templates
                    .unwrap_or(defaults.move_to_agent_templates),
            },
        };
        self._without_unmapped_reactants_params(Some(&params))
    }
    fn without_unmapped_reactants_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        self._without_unmapped_reactants_params(Some(params))
    }

    #[pyo3(signature = (params=None, *, threshold_unmapped_atoms=None, move_to_agent_templates=None))]
    fn without_unmapped_products(
        &self,
        params: Option<&ReactionTemplateRemovalParams>,
        threshold_unmapped_atoms: Option<f64>,
        move_to_agent_templates: Option<bool>,
    ) -> PyResult<ReactionTemplateRemoval> {
        if params.is_some()
            && (threshold_unmapped_atoms.is_some() || move_to_agent_templates.is_some())
        {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "parameter object and configuration keywords are mutually exclusive",
            ));
        }
        Ok(if params.is_some() {
            self._without_unmapped_products_params(params)
        } else {
            self._without_unmapped_products_keywords(
                threshold_unmapped_atoms,
                move_to_agent_templates,
            )
        })
    }
    #[doc = r#"{"role":"overload","target":"Reaction.without_unmapped_products"}"#]
    #[pyo3(signature = (params=None))]
    fn _without_unmapped_products_params(
        &self,
        params: Option<&ReactionTemplateRemovalParams>,
    ) -> ReactionTemplateRemoval {
        self.with_removal(match params {
            Some(params) => self
                .inner()
                .without_unmapped_products_with_params(&params.inner),
            None => self.inner().without_unmapped_products(),
        })
    }
    #[doc = r#"{"role":"overload","target":"Reaction.without_unmapped_products"}"#]
    #[pyo3(signature = (*, threshold_unmapped_atoms=None, move_to_agent_templates=None))]
    fn _without_unmapped_products_keywords(
        &self,
        threshold_unmapped_atoms: Option<f64>,
        move_to_agent_templates: Option<bool>,
    ) -> ReactionTemplateRemoval {
        let defaults = ck::ReactionTemplateRemovalParams::default();
        let params = ReactionTemplateRemovalParams {
            inner: ck::ReactionTemplateRemovalParams {
                threshold_unmapped_atoms: threshold_unmapped_atoms
                    .unwrap_or(defaults.threshold_unmapped_atoms),
                move_to_agent_templates: move_to_agent_templates
                    .unwrap_or(defaults.move_to_agent_templates),
            },
        };
        self._without_unmapped_products_params(Some(&params))
    }
    fn without_unmapped_products_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        self._without_unmapped_products_params(Some(params))
    }
}

#[pyclass(module = "fixed_projection._native", frozen)]
pub struct ReactionTemplateRemoval {
    inner: Arc<ck::ReactionTemplateRemoval>,
    match_owner: Mutex<Option<Py<SubstructMatchParams>>>,
}
#[pymethods]
impl ReactionTemplateRemoval {
    #[getter]
    fn reaction(&self) -> Reaction {
        let match_owner = Python::attach(|py| {
            self.match_owner
                .lock()
                .expect("reaction configuration lock")
                .as_ref()
                .map(|value| value.clone_ref(py))
        });
        Reaction {
            owner: ReactionOwner::Removal(Arc::clone(&self.inner)),
            match_owner: Mutex::new(match_owner),
        }
    }
    #[getter]
    fn removed_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .removed_templates()
            .iter()
            .enumerate()
            .map(|(index, _)| QueryGraph::from_removed_template(Arc::clone(&self.inner), index))
            .collect()
    }
    fn __traverse__(&self, visit: pyo3::gc::PyVisit<'_>) -> Result<(), pyo3::PyTraverseError> {
        if let Some(value) = self
            .match_owner
            .lock()
            .expect("reaction configuration lock")
            .as_ref()
        {
            visit.call(value)?;
        }
        Ok(())
    }
    fn __clear__(&self) {
        let previous = self
            .match_owner
            .lock()
            .expect("reaction configuration lock")
            .take();
        drop(previous);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matching_configuration_keeps_native_values_real_callbacks_and_collectible_owners() {
        use pyo3::ffi::c_str;
        use pyo3::types::PyDict;
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let native = ck::Reaction::from_smirks("[C:1]>>[C:1]").unwrap();
            let native_default = native.match_params().max_matches;
            let native_configured = native.with_match_params(&ck::SubstructMatchParams { max_matches: 2, num_threads: 2, ..Default::default() }).unwrap();
            let reaction = Py::new(py, Reaction::from_rust(native))?;
            let scope = PyDict::new(py);
            scope.set_item("reaction", reaction)?;
            scope.set_item("Params", py.get_type::<SubstructMatchParams>())?;
            scope.set_item("mol", Py::new(py, crate::native::Molecule { inner: ck::Molecule::from_smiles("CCC").unwrap() })?)?;
            scope.set_item("query", Py::new(py, QueryGraph::from_rust(ck::parse_smarts("C").unwrap()))?)?;
            scope.set_item("native_default", native_default)?;
            scope.set_item("native_configured", native_configured.match_params().max_matches)?;
            py.run(c_str!(r#"
import gc
import weakref
seen = []
def accept(target, indices):
    assert target is not mol
    assert target.to_binary() == mol.to_binary()
    seen.append(tuple(indices))
    return True
def atom(query_atom, atom):
    return True
def bond(query_bond, bond):
    return True
params = Params(max_matches=2, num_threads=2, final_match=accept, atom_match=atom, bond_match=bond)
configured = reaction.with_match_params(params)
assert configured.match_params().max_matches == native_configured
assert reaction.match_params().max_matches == native_default
params.max_matches = 1
params.final_match = None
assert configured.match_params().max_matches == 2
assert configured.match_params().final_match is accept
copied = configured.match_params()
copied.max_matches = 9
copied.final_match = None
assert configured.match_params().max_matches == 2
assert configured.match_params().final_match is accept
values = [configured.with_implicit_properties(False),
          configured.with_reactant_template(configured.reactant_template(0)),
          configured.with_product_template(configured.product_template(0)),
          configured.with_agent_template(configured.reactant_template(0)),
          configured.with_initialized(silent=True)]
removals = [configured.without_agents(), configured.without_unmapped_reactants(), configured.without_unmapped_products()]
values += [item.reaction for item in removals]
del configured, params, copied
gc.collect()
for value in values:
    stored = value.match_params()
    assert stored.max_matches == 2 and stored.num_threads == 2
    assert stored.final_match is accept and stored.atom_match is atom and stored.bond_match is bond
    seen.clear()
    assert len(mol.substruct_matches(query, stored)) == 2
    assert len(seen) == 2
marker = LookupError('retained reaction configuration callback')
def fail(target, indices):
    raise marker
stored.final_match = fail
try:
    mol.substruct_matches(query, stored)
except LookupError as error:
    assert error is marker
else:
    raise AssertionError('original callback exception lost')
assert values[-1].match_params().final_match is accept
class Callback:
    def __call__(self, target, indices):
        return True
def make_cycle(removal):
    callback = Callback()
    params = Params(final_match=callback)
    configured = reaction.with_match_params(params)
    if removal:
        owner = configured.without_agents()
        view = owner.reaction.reactant_template(0)
    else:
        owner = configured
        view = configured.reactant_template(0)
    callback.owner = owner
    return weakref.ref(callback), view
for removal in (False, True):
    retained, view = make_cycle(removal)
    gc.collect()
    assert retained() is None
    assert view.num_atoms() == 1
"#), Some(&scope), None)?;
            Ok(())
        }).unwrap();
    }

    #[test]
    fn native_removal_keeps_result_owners_and_configured_native_rules() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let graph = |text| ck::parse_smarts(text).unwrap();
            let native = ck::Reaction::from_templates(
                vec![graph("[C:1]"), graph("O"), graph("NN")],
                vec![graph("[C:1]"), graph("Cl")],
                vec![graph("[Na+]")],
            )
            .unwrap();
            let reaction = Reaction::from_rust(native.clone());
            let native_agents = native.without_agents();
            let removed = reaction.without_agents();
            let result_view = removed.reaction();
            let removed_view = removed.removed_templates().remove(0);
            assert!(std::ptr::eq(result_view.inner(), removed.inner.reaction()));
            assert!(std::ptr::eq(
                removed_view.inner(),
                &removed.inner.removed_templates()[0]
            ));
            assert_eq!(
                result_view.num_agent_templates(),
                native_agents.reaction().num_agent_templates()
            );
            assert_eq!(reaction.num_agent_templates(), 1);
            drop(removed);
            assert_eq!(result_view.num_reactant_templates(), 3);
            assert_eq!(removed_view.inner().num_atoms(), 1);
            for move_to_agent_templates in [true, false] {
                let params = ReactionTemplateRemovalParams {
                    inner: ck::ReactionTemplateRemovalParams {
                        threshold_unmapped_atoms: 0.5,
                        move_to_agent_templates,
                    },
                };
                let expected = native.without_unmapped_reactants_with_params(&params.inner);
                let via_params = reaction.without_unmapped_reactants(Some(&params), None, None)?;
                let via_keywords = reaction.without_unmapped_reactants(
                    None,
                    Some(0.5),
                    Some(move_to_agent_templates),
                )?;
                for actual in [
                    via_params,
                    via_keywords,
                    reaction.without_unmapped_reactants_with_params(&params),
                ] {
                    assert_eq!(
                        actual.reaction().num_reactant_templates(),
                        expected.reaction().num_reactant_templates()
                    );
                    assert_eq!(
                        actual.reaction().num_agent_templates(),
                        expected.reaction().num_agent_templates()
                    );
                    assert_eq!(
                        actual
                            .removed_templates()
                            .iter()
                            .map(|query| query.inner().num_atoms())
                            .collect::<Vec<_>>(),
                        expected
                            .removed_templates()
                            .iter()
                            .map(ck::QueryGraph::num_atoms)
                            .collect::<Vec<_>>()
                    );
                }
                let expected = native.without_unmapped_products_with_params(&params.inner);
                for actual in [
                    reaction.without_unmapped_products(Some(&params), None, None)?,
                    reaction.without_unmapped_products(
                        None,
                        Some(0.5),
                        Some(move_to_agent_templates),
                    )?,
                    reaction.without_unmapped_products_with_params(&params),
                ] {
                    assert_eq!(
                        actual.reaction().num_product_templates(),
                        expected.reaction().num_product_templates()
                    );
                    assert_eq!(
                        actual.reaction().num_agent_templates(),
                        expected.reaction().num_agent_templates()
                    );
                    assert_eq!(
                        actual.removed_templates().len(),
                        expected.removed_templates().len()
                    );
                }
                assert!(
                    reaction
                        .without_unmapped_products(Some(&params), Some(0.5), None)
                        .err()
                        .unwrap()
                        .is_instance_of::<pyo3::exceptions::PyTypeError>(py)
                );
            }
            assert_eq!(reaction.num_reactant_templates(), 3);
            assert_eq!(reaction.num_product_templates(), 2);
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn native_template_operations_preserve_order_views_and_errors() {
        Python::initialize();
        Python::attach(|py| -> PyResult<()> {
            let reactant = ck::parse_smarts("[C:1]O").unwrap();
            let product = ck::parse_smarts("[C:1]=O").unwrap();
            let agent = ck::parse_smarts("[Na+]").unwrap();
            let native = ck::Reaction::from_templates(
                vec![reactant.clone()],
                vec![product.clone()],
                vec![agent.clone()],
            )
            .unwrap();
            let reaction = Reaction::from_templates(
                py,
                vec![Py::new(py, QueryGraph::from_rust(reactant))?],
                vec![Py::new(py, QueryGraph::from_rust(product))?],
                vec![Py::new(py, QueryGraph::from_rust(agent))?],
            )?;
            assert_eq!(
                reaction.num_reactant_templates(),
                native.num_reactant_templates()
            );
            assert_eq!(
                reaction.num_product_templates(),
                native.num_product_templates()
            );
            assert_eq!(reaction.num_agent_templates(), native.num_agent_templates());
            let reactant_view = reaction.reactant_template(0)?;
            let product_view = reaction.product_templates().remove(0);
            let agent_view = reaction.agent_template(0)?;
            assert!(std::ptr::eq(
                reactant_view.inner(),
                &reaction.inner().reactant_templates()[0]
            ));
            assert!(std::ptr::eq(
                product_view.inner(),
                &reaction.inner().product_templates()[0]
            ));
            assert!(std::ptr::eq(
                agent_view.inner(),
                &reaction.inner().agent_templates()[0]
            ));
            let extended = reaction
                .with_reactant_template(&reactant_view)?
                .with_product_template(&product_view)?
                .with_agent_template(&agent_view)?;
            assert_eq!(
                (
                    extended.num_reactant_templates(),
                    extended.num_product_templates(),
                    extended.num_agent_templates()
                ),
                (2, 2, 2)
            );
            assert_eq!(reaction.num_reactant_templates(), 1);
            assert_eq!(
                reaction
                    .with_implicit_properties(true)
                    .implicit_properties(),
                true
            );
            assert_eq!(reaction.implicit_properties(), native.implicit_properties());
            let native_error = native.product_template(usize::MAX).unwrap_err();
            let display = native_error.to_string();
            let error = reaction.product_template(usize::MAX).err().unwrap();
            assert_eq!(error.value(py).str()?.to_str()?, display);
            assert_eq!(
                error.value(py).getattr("index")?.extract::<usize>()?,
                usize::MAX
            );
            assert_eq!(error.value(py).getattr("count")?.extract::<usize>()?, 1);
            assert!(
                error.value(py).getattr("role")?.extract::<ReactionRole>()?
                    == ReactionRole::Product
            );
            assert!(error.cause(py).is_none());
            drop(reaction);
            assert_eq!(reactant_view.inner().num_atoms(), 2);
            assert_eq!(product_view.inner().num_bonds(), 1);
            assert_eq!(agent_view.inner().num_atoms(), 1);
            let root = ck::ReactionModelError::Template {
                role: ck::ReactionRole::Agent,
                template: 7,
                source: ck::QueryGraphError::StereoGroup(ck::StereoGroupError::DuplicateAtom),
            };
            assert!(root.source().unwrap().is::<ck::QueryGraphError>());
            let error = ReactionModelError::from_rust(root);
            let cause = error.cause(py).unwrap();
            assert!(cause.is_instance_of::<QueryGraphError>(py));
            assert!(error.value(py).getattr("source")?.is(cause.value(py)));
            assert!(
                cause
                    .cause(py)
                    .unwrap()
                    .is_instance_of::<crate::model_errors::StereoGroupError>(py)
            );
            let payload = ck::QueryGraphError::CoordinateValidation(
                ck::CoordinateValidationError::Missing3DConformer { id: 9 },
            );
            assert!(payload.source().is_none());
            let error = QueryGraphError::from_rust(payload);
            assert!(error.cause(py).is_none());
            assert_eq!(
                error
                    .value(py)
                    .getattr("coordinate_validation")?
                    .getattr("id")?
                    .extract::<usize>()?,
                9
            );
            Ok(())
        })
        .unwrap();
    }
}
