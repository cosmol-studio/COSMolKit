//! Typed reaction projections. Algorithms and operation contracts live in cosmolkit.
use crate::canonical_search::{QueryGraph, SubstructMatchParams};
use crate::canonical_smiles_writer::CxSmilesFields;
use crate::drawing_binding::Molecule;
use ::cosmolkit as ck;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
#[cfg(feature = "stubgen")]
use pyo3_stub_gen::derive::{
    gen_stub_pyclass, gen_stub_pyclass_enum, gen_stub_pyfunction, gen_stub_pymethods,
};
use std::collections::BTreeMap;
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReactionRole {
    Reactant,
    Product,
    Agent,
}
impl From<ck::ReactionRole> for ReactionRole {
    fn from(value: ck::ReactionRole) -> Self {
        match value {
            ck::ReactionRole::Reactant => Self::Reactant,
            ck::ReactionRole::Product => Self::Product,
            ck::ReactionRole::Agent => Self::Agent,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReactionValidationSeverity {
    Warning,
    Error,
}
impl From<ck::ReactionValidationSeverity> for ReactionValidationSeverity {
    fn from(value: ck::ReactionValidationSeverity) -> Self {
        match value {
            ck::ReactionValidationSeverity::Warning => Self::Warning,
            ck::ReactionValidationSeverity::Error => Self::Error,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass_enum)]
#[pyclass(module = "cosmolkit", eq, eq_int, frozen, from_py_object)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReactionValidationIssueKind {
    MissingReactants,
    MissingProducts,
    DuplicateReactantMap,
    UnmappedReactant,
    DuplicateProductMap,
    MissingProductMapReactant,
    UnmappedProduct,
    UnmappedReactantMaps,
    MultipleCharge,
    MultipleHydrogenCount,
    MultipleMass,
    MultipleIsotope,
}
impl From<ck::ReactionValidationIssueKind> for ReactionValidationIssueKind {
    fn from(value: ck::ReactionValidationIssueKind) -> Self {
        match value {
            ck::ReactionValidationIssueKind::MissingReactants => Self::MissingReactants,
            ck::ReactionValidationIssueKind::MissingProducts => Self::MissingProducts,
            ck::ReactionValidationIssueKind::DuplicateReactantMap => Self::DuplicateReactantMap,
            ck::ReactionValidationIssueKind::UnmappedReactant => Self::UnmappedReactant,
            ck::ReactionValidationIssueKind::DuplicateProductMap => Self::DuplicateProductMap,
            ck::ReactionValidationIssueKind::MissingProductMapReactant => {
                Self::MissingProductMapReactant
            }
            ck::ReactionValidationIssueKind::UnmappedProduct => Self::UnmappedProduct,
            ck::ReactionValidationIssueKind::UnmappedReactantMaps => Self::UnmappedReactantMaps,
            ck::ReactionValidationIssueKind::MultipleCharge => Self::MultipleCharge,
            ck::ReactionValidationIssueKind::MultipleHydrogenCount => Self::MultipleHydrogenCount,
            ck::ReactionValidationIssueKind::MultipleMass => Self::MultipleMass,
            ck::ReactionValidationIssueKind::MultipleIsotope => Self::MultipleIsotope,
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionCoordinateSelection {
    pub(crate) inner: ck::ReactionCoordinateSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionCoordinateSelection {
    #[staticmethod]
    fn auto() -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::auto(),
        }
    }
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::two_d(id),
        }
    }
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::three_d(id),
        }
    }
    #[getter]
    fn id(&self) -> Option<usize> {
        self.inner.id()
    }
    #[getter]
    fn is_auto(&self) -> bool {
        self.inner.is_auto()
    }
    #[getter]
    fn is_2d(&self) -> bool {
        self.inner.is_2d()
    }
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionParseParams {
    pub(crate) inner: ck::ReactionParseParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionParseParams {
    #[new]
    #[pyo3(signature=(*,use_smiles=false,sanitize=false,replacements=None,allow_cxsmiles=true,strict_cxsmiles=true))]
    fn new(
        use_smiles: bool,
        sanitize: bool,
        replacements: Option<BTreeMap<String, String>>,
        allow_cxsmiles: bool,
        strict_cxsmiles: bool,
    ) -> Self {
        Self {
            inner: ck::ReactionParseParams::new(
                use_smiles,
                sanitize,
                replacements.unwrap_or_default(),
                allow_cxsmiles,
                strict_cxsmiles,
            ),
        }
    }
    #[getter]
    fn use_smiles(&self) -> bool {
        self.inner.use_smiles()
    }
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize()
    }
    #[getter]
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements().clone()
    }
    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles()
    }
    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionValidationParams {
    pub(crate) inner: ck::ReactionValidationParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationParams {
    #[new]
    #[pyo3(signature=(*,silent=false))]
    fn new(silent: bool) -> Self {
        Self {
            inner: ck::ReactionValidationParams::new(silent),
        }
    }
    #[getter]
    fn silent(&self) -> bool {
        self.inner.silent()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionSingleRunParams {
    pub(crate) inner: ck::ReactionSingleRunParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionSingleRunParams {
    #[new]
    #[pyo3(signature=(*,coordinate_selection=None))]
    fn new(coordinate_selection: Option<&ReactionCoordinateSelection>) -> Self {
        Self {
            inner: ck::ReactionSingleRunParams::new(
                coordinate_selection.map_or(ck::ReactionCoordinateSelection::Auto, |v| v.inner),
            ),
        }
    }
    #[getter]
    fn coordinate_selection(&self) -> ReactionCoordinateSelection {
        ReactionCoordinateSelection {
            inner: self.inner.coordinate_selection(),
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionRunParams {
    pub(crate) inner: ck::ReactionRunParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionRunParams {
    #[new]
    #[pyo3(signature=(*,max_products=1000,coordinate_selections=None))]
    fn new(
        max_products: u32,
        coordinate_selections: Option<Vec<PyRef<'_, ReactionCoordinateSelection>>>,
    ) -> Self {
        Self {
            inner: ck::ReactionRunParams::new(
                max_products,
                coordinate_selections
                    .unwrap_or_default()
                    .iter()
                    .map(|v| v.inner)
                    .collect(),
            ),
        }
    }
    #[getter]
    fn max_products(&self) -> u32 {
        self.inner.max_products()
    }
    #[getter]
    fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| ReactionCoordinateSelection { inner: *v })
            .collect()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionApplyParams {
    pub(crate) inner: ck::ReactionApplyParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionApplyParams {
    #[new]
    #[pyo3(signature=(*,remove_unmatched_atoms=true))]
    fn new(remove_unmatched_atoms: bool) -> Self {
        Self {
            inner: ck::ReactionApplyParams::new(remove_unmatched_atoms),
        }
    }
    #[getter]
    fn remove_unmatched_atoms(&self) -> bool {
        self.inner.remove_unmatched_atoms()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionTemplateRemovalParams {
    pub(crate) inner: ck::ReactionTemplateRemovalParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionTemplateRemovalParams {
    #[new]
    #[pyo3(signature=(*,threshold_unmapped_atoms=0.2,move_to_agent_templates=true))]
    fn new(threshold_unmapped_atoms: f64, move_to_agent_templates: bool) -> Self {
        Self {
            inner: ck::ReactionTemplateRemovalParams::new(
                threshold_unmapped_atoms,
                move_to_agent_templates,
            ),
        }
    }
    #[getter]
    fn threshold_unmapped_atoms(&self) -> f64 {
        self.inner.threshold_unmapped_atoms()
    }
    #[getter]
    fn move_to_agent_templates(&self) -> bool {
        self.inner.move_to_agent_templates()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionWriteParams {
    pub(crate) inner: ck::ReactionWriteParams,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionWriteParams {
    #[new]
    #[pyo3(signature=(*,canonical=false,do_isomeric_smiles=true,rooted_at_atom=None,include_dative_bonds=true,include_cx=false,cx_fields=None,coordinate_selections=None))]
    fn new(
        canonical: bool,
        do_isomeric_smiles: bool,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: bool,
        include_cx: bool,
        cx_fields: Option<&CxSmilesFields>,
        coordinate_selections: Option<Vec<PyRef<'_, ReactionCoordinateSelection>>>,
    ) -> Self {
        Self {
            inner: ck::ReactionWriteParams::new(
                canonical,
                do_isomeric_smiles,
                rooted_at_atom,
                include_dative_bonds,
                include_cx,
                cx_fields.map_or(ck::CxSmilesFields::ALL, |v| v.inner),
                coordinate_selections
                    .unwrap_or_default()
                    .iter()
                    .map(|v| v.inner)
                    .collect(),
            ),
        }
    }
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical()
    }
    #[getter]
    fn do_isomeric_smiles(&self) -> bool {
        self.inner.do_isomeric_smiles()
    }
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom()
    }
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds()
    }
    #[getter]
    fn include_cx(&self) -> bool {
        self.inner.include_cx()
    }
    #[getter]
    fn cx_fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.cx_fields(),
        }
    }
    #[getter]
    fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| ReactionCoordinateSelection { inner: *v })
            .collect()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Reaction {
    pub(crate) inner: ck::Reaction,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Reaction {
    fn run(
        &mut self,
        py: Python<'_>,
        reactants: Vec<PyRef<'_, Molecule>>,
        params: &ReactionRunParams,
    ) -> PyResult<Vec<Vec<Molecule>>> {
        let inputs: Vec<_> = reactants.iter().map(|value| &value.inner).collect();
        self.inner
            .run(&inputs, &params.inner)
            .map(product_sets)
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, &error))
    }
    #[new]
    fn new() -> Self {
        Self {
            inner: ck::Reaction::new(),
        }
    }
    #[staticmethod]
    fn from_smirks(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Reaction::from_smirks(text)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(py, &e))
    }
    #[staticmethod]
    fn from_smirks_with_params(
        py: Python<'_>,
        text: &str,
        params: &ReactionParseParams,
    ) -> PyResult<Self> {
        ck::Reaction::from_smirks_with_params(text, &params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(py, &e))
    }
    #[staticmethod]
    fn from_templates(
        py: Python<'_>,
        reactants: Vec<PyRef<'_, QueryGraph>>,
        products: Vec<PyRef<'_, QueryGraph>>,
        agents: Vec<PyRef<'_, QueryGraph>>,
    ) -> PyResult<Self> {
        ck::Reaction::from_templates(
            reactants.iter().map(|v| v.inner.clone()).collect(),
            products.iter().map(|v| v.inner.clone()).collect(),
            agents.iter().map(|v| v.inner.clone()).collect(),
        )
        .map(|inner| Self { inner })
        .map_err(|e| model_error(py, &e))
    }
    fn num_reactant_templates(&self) -> usize {
        self.inner.num_reactant_templates()
    }
    fn reactant_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .reactant_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    fn reactant_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .reactant_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    fn with_reactant_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_reactant_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    fn num_product_templates(&self) -> usize {
        self.inner.num_product_templates()
    }
    fn product_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .product_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    fn product_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .product_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    fn with_product_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_product_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    fn num_agent_templates(&self) -> usize {
        self.inner.num_agent_templates()
    }
    fn agent_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .agent_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    fn agent_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .agent_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    fn with_agent_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_agent_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    fn is_initialized(&self) -> bool {
        self.inner.is_initialized()
    }
    fn implicit_properties(&self) -> bool {
        self.inner.implicit_properties()
    }
    fn with_implicit_properties(&self, enabled: bool) -> Self {
        Self {
            inner: self.inner.with_implicit_properties(enabled),
        }
    }
    fn match_params(&self) -> SubstructMatchParams {
        SubstructMatchParams {
            inner: self.inner.match_params().clone(),
        }
    }
    fn with_match_params(&self, py: Python<'_>, params: &SubstructMatchParams) -> PyResult<Self> {
        self.inner
            .with_match_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    fn to_smirks(&self, py: Python<'_>) -> PyResult<String> {
        let v = self.inner.to_smirks().map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    fn to_smirks_with_params(
        &self,
        py: Python<'_>,
        params: &ReactionWriteParams,
    ) -> PyResult<String> {
        let v = self
            .inner
            .to_smirks_with_params(&params.inner)
            .map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    fn to_cx_smirks(&self, py: Python<'_>) -> PyResult<String> {
        let v = self.inner.to_cx_smirks().map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    fn to_cx_smirks_with_params(
        &self,
        py: Python<'_>,
        params: &ReactionWriteParams,
    ) -> PyResult<String> {
        let v = self
            .inner
            .to_cx_smirks_with_params(&params.inner)
            .map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    fn validate(&self, py: Python<'_>) -> PyResult<ReactionValidationReport> {
        self.inner
            .validate()
            .map(|inner| ReactionValidationReport { inner })
            .map_err(|e| validation_error(py, &e))
    }
    fn validate_with_params(
        &self,
        py: Python<'_>,
        params: &ReactionValidationParams,
    ) -> PyResult<ReactionValidationReport> {
        self.inner
            .validate_with_params(&params.inner)
            .map(|inner| ReactionValidationReport { inner })
            .map_err(|e| validation_error(py, &e))
    }
    fn with_initialized(&self, py: Python<'_>) -> PyResult<Reaction> {
        self.inner
            .with_initialized()
            .map(|inner| Reaction { inner })
            .map_err(|e| initialization_error(py, &e))
    }
    fn with_initialized_with_params(
        &self,
        py: Python<'_>,
        params: &ReactionValidationParams,
    ) -> PyResult<Reaction> {
        self.inner
            .with_initialized_with_params(&params.inner)
            .map(|inner| Reaction { inner })
            .map_err(|e| initialization_error(py, &e))
    }
    fn without_agents(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_agents(),
        }
    }
    fn without_unmapped_reactants(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_unmapped_reactants(),
        }
    }
    fn without_unmapped_reactants_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self
                .inner
                .without_unmapped_reactants_with_params(&params.inner),
        }
    }
    fn without_unmapped_products(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_unmapped_products(),
        }
    }
    fn without_unmapped_products_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self
                .inner
                .without_unmapped_products_with_params(&params.inner),
        }
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionTemplateRemoval {
    inner: ck::ReactionTemplateRemoval,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionTemplateRemoval {
    #[getter]
    fn reaction(&self) -> Reaction {
        Reaction {
            inner: self.inner.reaction().clone(),
        }
    }
    #[getter]
    fn removed_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .removed_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionApplyResult {
    inner: ck::ReactionApplyResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionApplyResult {
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    #[getter]
    fn changed(&self) -> bool {
        self.inner.changed()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionValidationReport {
    inner: ck::ReactionValidationReport,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationReport {
    #[getter]
    fn warnings(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .warnings()
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    #[getter]
    fn errors(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .errors()
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    #[getter]
    fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
    #[getter]
    fn num_warnings(&self) -> usize {
        self.inner.num_warnings()
    }
    #[getter]
    fn num_errors(&self) -> usize {
        self.inner.num_errors()
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionValidationIssue {
    inner: ck::ReactionValidationIssue,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationIssue {
    #[getter]
    fn kind(&self) -> ReactionValidationIssueKind {
        self.inner.kind().into()
    }
    #[getter]
    fn severity(&self) -> ReactionValidationSeverity {
        self.inner.severity().into()
    }
    #[getter]
    fn role(&self) -> Option<ReactionRole> {
        self.inner.role().map(Into::into)
    }
    #[getter]
    fn template(&self) -> Option<usize> {
        self.inner.template()
    }
    #[getter]
    fn atom(&self) -> Option<usize> {
        self.inner.atom().map(|v| v.index())
    }
    #[getter]
    fn map(&self) -> Option<i32> {
        self.inner.map()
    }
    #[getter]
    fn maps(&self) -> Vec<i32> {
        self.inner.maps().to_vec()
    }
    #[getter]
    fn detail(&self) -> String {
        self.inner.detail().to_owned()
    }
}
pyo3::create_exception!(cosmolkit, ReactionModelError, PyValueError);
pub(crate) fn model_error(py: Python<'_>, source: &ck::ReactionModelError) -> PyErr {
    use ck::ReactionModelError as E;
    let kind = match source {
        E::Template { .. } => "Template",
        E::TemplateIndex { .. } => "TemplateIndex",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionModelError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::Template { role, template, .. } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
            }
            E::TemplateIndex {
                role, index, count, ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("index", index)?;
                exception_value.setattr("count", count)?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionParseError, PyValueError);
pub(crate) fn parse_error(py: Python<'_>, source: &ck::ReactionParseError) -> PyErr {
    use ck::ReactionParseError as E;
    let kind = match source {
        E::ParserCleanup { .. } => "ParserCleanup",
        E::TemplateAtomBounds { .. } => "TemplateAtomBounds",
        E::StereoDegreeArithmetic { .. } => "StereoDegreeArithmetic",
        E::StereoInversionFlag { .. } => "StereoInversionFlag",
        E::TemplateProperty { .. } => "TemplateProperty",
        E::Separators { .. } => "Separators",
        E::MultiStep { .. } => "MultiStep",
        E::ComponentBounds { .. } => "ComponentBounds",
        E::OffsetOverflow => "OffsetOverflow",
        E::Smarts { .. } => "Smarts",
        E::Smiles { .. } => "Smiles",
        E::ComponentEncoding { .. } => "ComponentEncoding",
        E::Model { .. } => "Model",
        E::AgentFragments { .. } => "AgentFragments",
        E::CxParse { .. } => "CxParse",
        E::CxLowering { .. } => "CxLowering",
        E::StereoOrder { .. } => "StereoOrder",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionParseError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::ParserCleanup { role, template, .. } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
            }
            E::TemplateAtomBounds {
                role,
                template,
                atom,
                atom_count,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("atom_count", atom_count)?;
            }
            E::StereoDegreeArithmetic {
                reactant, product, ..
            } => {
                exception_value.setattr("reactant", reactant)?;
                exception_value.setattr("product", product)?;
            }
            E::StereoInversionFlag {
                product,
                atom,
                flag,
                ..
            } => {
                exception_value.setattr("product", product)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("flag", flag)?;
            }
            E::Separators { count, .. } => {
                exception_value.setattr("count", count)?;
            }
            E::MultiStep { count, .. } => {
                exception_value.setattr("count", count)?;
            }
            E::ComponentBounds { start, end, .. } => {
                exception_value.setattr("start", start)?;
                exception_value.setattr("end", end)?;
            }
            E::Smarts {
                role,
                template,
                text,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value
                    .setattr("text", crate::canonical_sdf::decode_source_text(py, text)?)?;
            }
            E::Smiles {
                role,
                template,
                text,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value
                    .setattr("text", crate::canonical_sdf::decode_source_text(py, text)?)?;
            }
            E::ComponentEncoding {
                role,
                template,
                text,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value
                    .setattr("text", crate::canonical_sdf::decode_source_text(py, text)?)?;
            }
            E::Model { role, template, .. } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
            }
            E::CxParse {
                role,
                template,
                start_atom,
                start_bond,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("start_atom", start_atom)?;
                exception_value.setattr("start_bond", start_bond)?;
            }
            E::CxLowering {
                role,
                template,
                start_atom,
                start_bond,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("start_atom", start_atom)?;
                exception_value.setattr("start_bond", start_bond)?;
            }
            E::StereoOrder { product, atom, .. } => {
                exception_value.setattr("product", product)?;
                exception_value.setattr("atom", atom.index())?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionRunError, PyValueError);
pub(crate) fn run_error(py: Python<'_>, source: &ck::ReactionRunError) -> PyErr {
    use ck::ReactionRunError as E;
    let kind = match source {
        E::Product { .. } => "Product",
        E::CoordinateSelectionArity { .. } => "CoordinateSelectionArity",
        E::Initialization(..) => "Initialization",
        E::NeedsInitialization => "NeedsInitialization",
        E::ReactantArity { .. } => "ReactantArity",
        E::ReactantTemplateIndex { .. } => "ReactantTemplateIndex",
        E::Matching { .. } => "Matching",
        E::MatchAtomIndex { .. } => "MatchAtomIndex",
        E::CombinationLevel { .. } => "CombinationLevel",
        E::CombinationSize { .. } => "CombinationSize",
        E::EmptyCombinationLevels => "EmptyCombinationLevels",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionRunError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    // Transparent Rust errors still retain the concrete wrapper at the language boundary.
    match source {
        E::Initialization(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        _ => (),
    }
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::Product { set, template, .. } => {
                exception_value.setattr("set", set)?;
                exception_value.setattr("template", template)?;
            }
            E::CoordinateSelectionArity {
                expected, actual, ..
            } => {
                exception_value.setattr("expected", expected)?;
                exception_value.setattr("actual", actual)?;
            }
            E::ReactantArity {
                expected, actual, ..
            } => {
                exception_value.setattr("expected", expected)?;
                exception_value.setattr("actual", actual)?;
            }
            E::ReactantTemplateIndex { index, count, .. } => {
                exception_value.setattr("index", index)?;
                exception_value.setattr("count", count)?;
            }
            E::Matching {
                reactant, template, ..
            } => {
                exception_value.setattr("reactant", reactant)?;
                exception_value.setattr("template", template)?;
            }
            E::MatchAtomIndex {
                reactant,
                template,
                atom,
                atom_count,
                ..
            } => {
                exception_value.setattr("reactant", reactant)?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom)?;
                exception_value.setattr("atom_count", atom_count)?;
            }
            E::CombinationLevel { level, count, .. } => {
                exception_value.setattr("level", level)?;
                exception_value.setattr("count", count)?;
            }
            E::CombinationSize {
                expected, actual, ..
            } => {
                exception_value.setattr("expected", expected)?;
                exception_value.setattr("actual", actual)?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionApplyError, PyValueError);
pub(crate) fn apply_error(py: Python<'_>, source: &ck::ReactionApplyError) -> PyErr {
    use ck::ReactionApplyError as E;
    let kind = match source {
        E::SourceText(..) => "SourceText",
        E::Initialization(..) => "Initialization",
        E::ApplicabilityArity { .. } => "ApplicabilityArity",
        E::AddsProductAtom { .. } => "AddsProductAtom",
        E::Matching(..) => "Matching",
        E::Product(..) => "Product",
        E::Edit(..) => "Edit",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionApplyError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    // Transparent Rust errors still retain the concrete wrapper at the language boundary.
    match source {
        E::SourceText(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Initialization(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Matching(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Product(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Edit(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        _ => (),
    }
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::ApplicabilityArity {
                reactants,
                products,
                ..
            } => {
                exception_value.setattr("reactants", reactants)?;
                exception_value.setattr("products", products)?;
            }
            E::AddsProductAtom { atom, .. } => {
                exception_value.setattr("atom", atom.index())?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionProductError, PyValueError);
pub(crate) fn product_error(py: Python<'_>, source: &ck::ReactionProductError) -> PyErr {
    use ck::ReactionProductError as E;
    let kind = match source {
        E::Coordinate(..) => "Coordinate",
        E::TopologyEdit(..) => "TopologyEdit",
        E::Adjacency(..) => "Adjacency",
        E::RingFinding(..) => "RingFinding",
        E::SourceUInt(..) => "SourceUInt",
        E::Valence(..) => "Valence",
        E::StereoOrder(..) => "StereoOrder",
        E::DoubleBondStereo(..) => "DoubleBondStereo",
        E::CoordinateSelection(..) => "CoordinateSelection",
        E::CarrierIdentity(..) => "CarrierIdentity",
        E::MoleculeProperty(..) => "MoleculeProperty",
        E::AtomProperty(..) => "AtomProperty",
        E::BondValue(..) => "BondValue",
        E::TemplateProperty(..) => "TemplateProperty",
        E::StereoGetter { .. } => "StereoGetter",
        E::PropertyInt { .. } => "PropertyInt",
        E::PropertyUInt { .. } => "PropertyUInt",
        E::MissingProperty { .. } => "MissingProperty",
        E::Invariant { .. } => "Invariant",
        E::RowOverflow { .. } => "RowOverflow",
        E::QueryReactantBond { .. } => "QueryReactantBond",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionProductError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    // Transparent Rust errors still retain the concrete wrapper at the language boundary.
    match source {
        E::Coordinate(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::TopologyEdit(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Adjacency(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::RingFinding(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::SourceUInt(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::Valence(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::StereoOrder(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::DoubleBondStereo(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::CoordinateSelection(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::CarrierIdentity(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::MoleculeProperty(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::AtomProperty(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::BondValue(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        E::TemplateProperty(cause) => {
            error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause)))
        }
        _ => (),
    }
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::StereoGetter { atom, .. } => {
                exception_value.setattr("atom", atom.index())?;
            }
            E::PropertyInt { atom, key, .. } => {
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("key", key)?;
            }
            E::PropertyUInt { atom, key, .. } => {
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("key", key)?;
            }
            E::MissingProperty { atom, key, .. } => {
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("key", key)?;
            }
            E::Invariant {
                stage,
                detail,
                reactant_atom,
                product_atom,
                bond,
                ..
            } => {
                exception_value.setattr("stage", stage)?;
                exception_value.setattr("detail", detail)?;
                exception_value.setattr("reactant_atom", reactant_atom)?;
                exception_value.setattr("product_atom", product_atom)?;
                exception_value.setattr("bond", bond.map(|v| v.index()))?;
            }
            E::RowOverflow { kind, index, .. } => {
                exception_value.setattr("row_kind", kind)?;
                exception_value.setattr("index", index)?;
            }
            E::QueryReactantBond { bond, .. } => {
                exception_value.setattr("bond", bond.index())?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionWriteError, PyValueError);
pub(crate) fn write_error(py: Python<'_>, source: &ck::ReactionWriteError) -> PyErr {
    use ck::ReactionWriteError as E;
    let kind = match source {
        E::Cx(..) => "Cx",
        E::Template { .. } => "Template",
        E::Connectivity { .. } => "Connectivity",
        E::CoordinateSelectionArity { .. } => "CoordinateSelectionArity",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionWriteError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    // Transparent Rust errors still retain the concrete wrapper at the language boundary.
    match source {
        E::Cx(cause) => error.set_cause(py, Some(crate::canonical_values::source_pyerr(py, cause))),
        _ => (),
    }
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::Template { role, template, .. } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
            }
            E::Connectivity { role, template, .. } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
            }
            E::CoordinateSelectionArity {
                expected, actual, ..
            } => {
                exception_value.setattr("expected", expected)?;
                exception_value.setattr("actual", actual)?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionValidationError, PyValueError);
pub(crate) fn validation_error(py: Python<'_>, source: &ck::ReactionValidationError) -> PyErr {
    use ck::ReactionValidationError as E;
    let kind = match source {
        E::Property { .. } => "Property",
        E::MapOverflow { .. } => "MapOverflow",
        E::AtomBounds { .. } => "AtomBounds",
        E::MissingReactingAtom { .. } => "MissingReactingAtom",
        E::QueryArithmetic { .. } => "QueryArithmetic",
        E::Annotation { .. } => "Annotation",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionValidationError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::Property {
                role,
                template,
                atom,
                property,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("property", property)?;
            }
            E::MapOverflow {
                role,
                template,
                atom,
                property,
                value,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("property", property)?;
                exception_value.setattr("value", value)?;
            }
            E::AtomBounds {
                role,
                template,
                atom,
                atom_count,
                ..
            } => {
                exception_value.setattr("role", ReactionRole::from(*role))?;
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("atom_count", atom_count)?;
            }
            E::MissingReactingAtom {
                template,
                atom,
                map,
                ..
            } => {
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("map", map)?;
            }
            E::QueryArithmetic { template, atom, .. } => {
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
            }
            E::Annotation {
                template,
                atom,
                property,
                ..
            } => {
                exception_value.setattr("template", template)?;
                exception_value.setattr("atom", atom.index())?;
                exception_value.setattr("property", property)?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
pyo3::create_exception!(cosmolkit, ReactionInitializationError, PyValueError);
pub(crate) fn initialization_error(
    py: Python<'_>,
    source: &ck::ReactionInitializationError,
) -> PyErr {
    use ck::ReactionInitializationError as E;
    let kind = match source {
        E::Validation { .. } => "Validation",
        E::Invalid { .. } => "Invalid",
    };
    let error = crate::canonical_values::annotate(
        py,
        ReactionInitializationError::new_err(source.to_string()),
        "reaction",
        kind,
        source,
    );
    let attrs = || -> PyResult<()> {
        let exception_value = error.value(py);
        match source {
            E::Invalid { report, .. } => {
                exception_value.setattr(
                    "report",
                    ReactionValidationReport {
                        inner: report.clone(),
                    },
                )?;
            }
            _ => (),
        };
        Ok(())
    };
    match attrs() {
        Ok(()) => error,
        Err(e) => e,
    }
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smirks(py: Python<'_>, text: &str) -> PyResult<Reaction> {
    Reaction::from_smirks(py, text)
}
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smirks_with_params(
    py: Python<'_>,
    text: &str,
    params: &ReactionParseParams,
) -> PyResult<Reaction> {
    Reaction::from_smirks_with_params(py, text, params)
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ReactionRole>()?;
    module.add_class::<ReactionValidationSeverity>()?;
    module.add_class::<ReactionValidationIssueKind>()?;
    module.add_class::<ReactionCoordinateSelection>()?;
    module.add_class::<ReactionParseParams>()?;
    module.add_class::<ReactionValidationParams>()?;
    module.add_class::<ReactionSingleRunParams>()?;
    module.add_class::<ReactionRunParams>()?;
    module.add_class::<ReactionApplyParams>()?;
    module.add_class::<ReactionTemplateRemovalParams>()?;
    module.add_class::<ReactionWriteParams>()?;
    module.add_class::<Reaction>()?;
    module.add_class::<ReactionTemplateRemoval>()?;
    module.add_class::<ReactionApplyResult>()?;
    module.add_class::<ReactionValidationReport>()?;
    module.add_class::<ReactionValidationIssue>()?;
    module.add(
        "ReactionModelError",
        module.py().get_type::<ReactionModelError>(),
    )?;
    module.add(
        "ReactionParseError",
        module.py().get_type::<ReactionParseError>(),
    )?;
    module.add(
        "ReactionRunError",
        module.py().get_type::<ReactionRunError>(),
    )?;
    module.add(
        "ReactionApplyError",
        module.py().get_type::<ReactionApplyError>(),
    )?;
    module.add(
        "ReactionProductError",
        module.py().get_type::<ReactionProductError>(),
    )?;
    module.add(
        "ReactionWriteError",
        module.py().get_type::<ReactionWriteError>(),
    )?;
    module.add(
        "ReactionValidationError",
        module.py().get_type::<ReactionValidationError>(),
    )?;
    module.add(
        "ReactionInitializationError",
        module.py().get_type::<ReactionInitializationError>(),
    )?;
    module.add_function(wrap_pyfunction!(parse_smirks, module)?)?;
    module.add_function(wrap_pyfunction!(parse_smirks_with_params, module)?)?;
    Ok(())
}
pub(crate) fn product_sets(sets: Vec<Vec<ck::Molecule>>) -> Vec<Vec<Molecule>> {
    sets.into_iter()
        .map(|set| set.into_iter().map(|inner| Molecule { inner }).collect())
        .collect()
}
pub(crate) fn apply_result(inner: ck::ReactionApplyResult) -> ReactionApplyResult {
    ReactionApplyResult { inner }
}
