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
/// Source template role, retaining source vector order.
///
/// Declared values: ``Reactant``, ``Product``, ``Agent``.
#[cosmolkit_macros::python_enum]
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
/// Severity of a reaction validation issue: warning or error.
///
/// Declared values: ``Warning``, ``Error``.
#[cosmolkit_macros::python_enum]
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
/// Typed category of reaction template or atom-map validation failure.
///
/// Declared values: ``MissingReactants``, ``MissingProducts``, ``DuplicateReactantMap``, ``UnmappedReactant``, ``DuplicateProductMap``, ``MissingProductMapReactant``, ``UnmappedProduct``, ``UnmappedReactantMaps``, ``MultipleCharge``, ``MultipleHydrogenCount``, ``MultipleMass``, ``MultipleIsotope``.
#[cosmolkit_macros::python_enum]
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
/// Explicit reaction coordinate selection: automatic, stored 2D, or a specified 3D conformer ID.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen, eq)]
#[derive(PartialEq)]
pub(crate) struct ReactionCoordinateSelection {
    pub(crate) inner: ck::ReactionCoordinateSelection,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionCoordinateSelection {
    /// Construct a coordinate selection for automatic coordinate resolution, rejecting an ambiguous selection.
    #[staticmethod]
    fn auto() -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::auto(),
        }
    }
    /// Construct a coordinate selection for the separate stored 2D conformer.
    #[staticmethod]
    fn two_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::two_d(id),
        }
    }
    /// Construct a coordinate selection for the stored 3D conformer with the supplied ID.
    #[staticmethod]
    fn three_d(id: usize) -> Self {
        Self {
            inner: ck::ReactionCoordinateSelection::three_d(id),
        }
    }
    /// Requested stored 2D or 3D conformer ID, or None for automatic selection.
    #[getter]
    fn id(&self) -> Option<usize> {
        self.inner.id()
    }
    /// Return whether this selection uses automatic coordinate resolution.
    #[getter]
    fn is_auto(&self) -> bool {
        self.inner.is_auto()
    }
    /// Return whether this selection uses stored 2D coordinates.
    #[getter]
    fn is_2d(&self) -> bool {
        self.inner.is_2d()
    }
    /// Whether the conformer/input is designated three-dimensional.
    #[getter]
    fn is_3d(&self) -> bool {
        self.inner.is_3d()
    }
}
/// Writable configuration for reaction SMIRKS parsing.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionParseParams {
    pub(crate) inner: ck::ReactionParseParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionParseParams {
    /// Configure reaction SMIRKS parsing; omitted fields use the defaults shown in the signature.
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
    /// Whether reaction components are parsed as SMILES rather than SMARTS.
    #[getter]
    fn use_smiles(&self) -> bool {
        self.inner.use_smiles()
    }
    /// Apply to a new molecule and return the result: perform the selected chemical sanitization stages. The source molecule is unchanged.
    #[getter]
    fn sanitize(&self) -> bool {
        self.inner.sanitize()
    }
    /// Text substitutions applied before parsing.
    #[getter]
    fn replacements(&self) -> BTreeMap<String, String> {
        self.inner.replacements().clone()
    }
    /// Whether a CX extension following the graph notation is parsed.
    #[getter]
    fn allow_cxsmiles(&self) -> bool {
        self.inner.allow_cxsmiles()
    }
    /// Whether malformed CX extension data is rejected.
    #[getter]
    fn strict_cxsmiles(&self) -> bool {
        self.inner.strict_cxsmiles()
    }
}
/// Writable configuration for reaction validation diagnostics.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionValidationParams {
    pub(crate) inner: ck::ReactionValidationParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationParams {
    /// Configure reaction validation diagnostics; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,silent=false))]
    fn new(silent: bool) -> Self {
        Self {
            inner: ck::ReactionValidationParams::new(silent),
        }
    }
    /// Whether reaction validation suppresses diagnostic output.
    #[getter]
    fn silent(&self) -> bool {
        self.inner.silent()
    }
}
/// Writable configuration for single-reactant reaction coordinate selection.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionSingleRunParams {
    pub(crate) inner: ck::ReactionSingleRunParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionSingleRunParams {
    /// Configure single-reactant reaction coordinate selection; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,coordinate_selection=None))]
    fn new(coordinate_selection: Option<&ReactionCoordinateSelection>) -> Self {
        Self {
            inner: ck::ReactionSingleRunParams::new(
                coordinate_selection.map_or(ck::ReactionCoordinateSelection::Auto, |v| v.inner),
            ),
        }
    }
    /// Explicit selection of stored 2D coordinates, a 3D conformer, or automatic resolution.
    #[getter]
    fn coordinate_selection(&self) -> ReactionCoordinateSelection {
        ReactionCoordinateSelection {
            inner: self.inner.coordinate_selection(),
        }
    }
}
/// Writable configuration for multi-reactant reaction execution.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionRunParams {
    pub(crate) inner: ck::ReactionRunParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionRunParams {
    /// Configure multi-reactant reaction execution; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,max_products=1000,coordinate_selections=None,copy_atom_properties=false))]
    fn new(
        max_products: u32,
        coordinate_selections: Option<Vec<PyRef<'_, ReactionCoordinateSelection>>>,
        copy_atom_properties: bool,
    ) -> Self {
        Self {
            inner: ck::ReactionRunParams::new(
                max_products,
                coordinate_selections
                    .unwrap_or_default()
                    .iter()
                    .map(|v| v.inner)
                    .collect(),
                copy_atom_properties,
            ),
        }
    }
    /// RDKit default 1000; zero leaves matching and combinations unlimited.
    #[getter]
    fn max_products(&self) -> u32 {
        self.inner.max_products()
    }
    /// CK extension: fill missing product user properties from atom origins.
    /// False preserves native behavior; it never clears existing properties.
    #[getter]
    fn copy_atom_properties(&self) -> bool {
        self.inner.copy_atom_properties()
    }
    /// Empty selects Auto for each input; otherwise input arity must match.
    #[getter]
    fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| ReactionCoordinateSelection { inner: *v })
            .collect()
    }
}
/// Writable configuration for single-reactant reaction application.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionApplyParams {
    pub(crate) inner: ck::ReactionApplyParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionApplyParams {
    /// Configure single-reactant reaction application; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,remove_unmatched_atoms=true))]
    fn new(remove_unmatched_atoms: bool) -> Self {
        Self {
            inner: ck::ReactionApplyParams::new(remove_unmatched_atoms),
        }
    }
    /// Whether atoms not matched by the reaction template are removed during application.
    #[getter]
    fn remove_unmatched_atoms(&self) -> bool {
        self.inner.remove_unmatched_atoms()
    }
}
/// Writable configuration for filtering unmapped reaction templates.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionTemplateRemovalParams {
    pub(crate) inner: ck::ReactionTemplateRemovalParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionTemplateRemovalParams {
    /// Configure filtering unmapped reaction templates; omitted fields use the defaults shown in the signature.
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
    /// Unmapped-atom fraction above which a reaction template is removed.
    #[getter]
    fn threshold_unmapped_atoms(&self) -> f64 {
        self.inner.threshold_unmapped_atoms()
    }
    /// Whether removed reactant/product templates are retained as agent templates.
    #[getter]
    fn move_to_agent_templates(&self) -> bool {
        self.inner.move_to_agent_templates()
    }
}
/// Writable configuration for reaction SMIRKS/CX output.
///
/// Set fields in the constructor or assign them afterward. Omitted values use the
/// documented constructor defaults; invalid assignments leave the previous value unchanged.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", dict, weakref)]
pub(crate) struct ReactionWriteParams {
    pub(crate) inner: ck::ReactionWriteParams,
}
#[cosmolkit_macros::python_configuration]
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionWriteParams {
    /// Configure reaction SMIRKS/CX output; omitted fields use the defaults shown in the signature.
    #[new]
    #[pyo3(signature=(*,canonical=false,isomeric_smiles=true,rooted_at_atom=None,include_dative_bonds=true,include_cx=false,cx_fields=None,coordinate_selections=None))]
    fn new(
        canonical: bool,
        isomeric_smiles: bool,
        rooted_at_atom: Option<usize>,
        include_dative_bonds: bool,
        include_cx: bool,
        cx_fields: Option<&CxSmilesFields>,
        coordinate_selections: Option<Vec<PyRef<'_, ReactionCoordinateSelection>>>,
    ) -> Self {
        Self {
            inner: ck::ReactionWriteParams::new(
                canonical,
                isomeric_smiles,
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
    /// Whether canonical atom traversal is used for output.
    #[getter]
    fn canonical(&self) -> bool {
        self.inner.canonical()
    }
    /// Whether isotope and stereochemical information is included in the output notation.
    #[getter]
    fn isomeric_smiles(&self) -> bool {
        self.inner.isomeric_smiles()
    }
    /// Atom index at which output traversal starts, or None for the default traversal.
    #[getter]
    fn rooted_at_atom(&self) -> Option<usize> {
        self.inner.rooted_at_atom()
    }
    /// Whether dative bonds are included in the requested graph operation/output.
    #[getter]
    fn include_dative_bonds(&self) -> bool {
        self.inner.include_dative_bonds()
    }
    /// Whether reaction output includes supported CX annotations.
    #[getter]
    fn include_cx(&self) -> bool {
        self.inner.include_cx()
    }
    /// Bit mask selecting CX annotation fields.
    #[getter]
    fn cx_fields(&self) -> CxSmilesFields {
        CxSmilesFields {
            inner: self.inner.cx_fields(),
        }
    }
    /// Per template in original reactant, agent, product insertion order.
    #[getter]
    fn coordinate_selections(&self) -> Vec<ReactionCoordinateSelection> {
        self.inner
            .coordinate_selections()
            .iter()
            .map(|v| ReactionCoordinateSelection { inner: *v })
            .collect()
    }
}
/// Owned reaction template with reactant, agent and product QueryGraph values.
///
/// Construct from SMIRKS or explicit templates. run() accepts the ordered reactant
/// list and returns product sets without mutating those reactants. Template edits
/// use value-returning methods.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit")]
pub(crate) struct Reaction {
    pub(crate) inner: ck::Reaction,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl Reaction {
    /// Run this reaction on the ordered reactant Molecule list. Return a list of product sets, each a list of Molecule values in product-template order. Reactants are unchanged; params and keyword options are mutually exclusive.
    #[pyo3(signature = (reactants, params=None))]
    fn run(
        &mut self,
        py: Python<'_>,
        reactants: Vec<PyRef<'_, Molecule>>,
        params: Option<&ReactionRunParams>,
    ) -> PyResult<Vec<Vec<Molecule>>> {
        let defaults;
        let params = match params {
            Some(params) => &params.inner,
            None => {
                defaults = ck::ReactionRunParams::default();
                &defaults
            }
        };
        let inputs: Vec<_> = reactants.iter().map(|value| &value.inner).collect();
        self.inner
            .run(&inputs, params)
            .map(product_sets)
            .map_err(|error| crate::drawing_binding::operation_pyerr(py, &error))
    }
    /// Construct a Reaction value from the supplied inputs.
    #[new]
    fn new() -> Self {
        Self {
            inner: ck::Reaction::new(),
        }
    }
    /// Parse reaction SMIRKS text into a Reaction; parsing failures raise ReactionParseError.
    #[staticmethod]
    fn from_smirks(py: Python<'_>, text: &str) -> PyResult<Self> {
        ck::Reaction::from_smirks(text)
            .map(|inner| Self { inner })
            .map_err(|e| parse_error(py, &e))
    }
    /// Parse reaction SMIRKS text into a Reaction; parsing failures raise ReactionParseError. Uses the supplied configuration object.
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
    /// Construct a Reaction from reactant, product and agent query templates.
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
    /// Return the number of reactant query templates.
    fn num_reactant_templates(&self) -> usize {
        self.inner.num_reactant_templates()
    }
    /// Return reactant QueryGraph template(s) in stored template order; indices are zero-based.
    fn reactant_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .reactant_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    /// Return reactant QueryGraph template(s) in stored template order; indices are zero-based.
    fn reactant_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .reactant_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    /// Return a new Reaction with the supplied reactant QueryGraph template; leave the source reaction unchanged.
    fn with_reactant_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_reactant_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    /// Return the number of product query templates.
    fn num_product_templates(&self) -> usize {
        self.inner.num_product_templates()
    }
    /// Return product QueryGraph template(s) in stored template order; indices are zero-based.
    fn product_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .product_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    /// Return product QueryGraph template(s) in stored template order; indices are zero-based.
    fn product_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .product_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    /// Return a new Reaction with the supplied product QueryGraph template; leave the source reaction unchanged.
    fn with_product_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_product_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    /// Return the number of agent query templates.
    fn num_agent_templates(&self) -> usize {
        self.inner.num_agent_templates()
    }
    /// Return agent QueryGraph template(s) in stored template order; indices are zero-based.
    fn agent_templates(&self) -> Vec<QueryGraph> {
        self.inner
            .agent_templates()
            .iter()
            .cloned()
            .map(|inner| QueryGraph { inner })
            .collect()
    }
    /// Return agent QueryGraph template(s) in stored template order; indices are zero-based.
    fn agent_template(&self, py: Python<'_>, index: usize) -> PyResult<QueryGraph> {
        self.inner
            .agent_template(index)
            .map(|v| QueryGraph { inner: v.clone() })
            .map_err(|e| model_error(py, &e))
    }
    /// Return a new Reaction with the supplied agent QueryGraph template; leave the source reaction unchanged.
    fn with_agent_template(&self, py: Python<'_>, template: &QueryGraph) -> PyResult<Self> {
        self.inner
            .with_agent_template(template.inner.clone())
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    /// Return whether reactant matchers have been initialized.
    fn is_initialized(&self) -> bool {
        self.inner.is_initialized()
    }
    /// Return whether unspecified product atom properties may be inherited from reactants.
    fn implicit_properties(&self) -> bool {
        self.inner.implicit_properties()
    }
    /// Return a new Reaction with the requested product-property inheritance policy.
    fn with_implicit_properties(&self, enabled: bool) -> Self {
        Self {
            inner: self.inner.with_implicit_properties(enabled),
        }
    }
    /// Return the substructure-matching configuration used to find reactant template matches.
    fn match_params(&self) -> SubstructMatchParams {
        SubstructMatchParams::from_inner(self.inner.match_params().clone())
    }
    /// Return a new Reaction using the supplied SubstructMatchParams for reactant template matching.
    fn with_match_params(&self, py: Python<'_>, params: &SubstructMatchParams) -> PyResult<Self> {
        if params.has_python_callbacks() {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(
                "Python match callbacks are supported by Molecule substructure matching; persistent Reaction callbacks are not yet projected",
            ));
        }
        self.inner
            .with_match_params(&params.inner)
            .map(|inner| Self { inner })
            .map_err(|e| model_error(py, &e))
    }
    /// Return reaction SMIRKS text with reactants, agents and products; leave the Reaction unchanged.
    fn to_smirks(&self, py: Python<'_>) -> PyResult<String> {
        let v = self.inner.to_smirks().map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    /// Return reaction SMIRKS text with reactants, agents and products; leave the Reaction unchanged. Uses the supplied configuration object.
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
    /// Return reaction SMIRKS text with supported CX annotations.
    fn to_cx_smirks(&self, py: Python<'_>) -> PyResult<String> {
        let v = self.inner.to_cx_smirks().map_err(|e| write_error(py, &e))?;
        crate::canonical_sdf::decode_source_text(py, &v)
    }
    /// Return reaction SMIRKS text with supported CX annotations. Uses the supplied configuration object.
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
    /// Return a ReactionValidationReport with warnings and errors for template and atom-map consistency; do not silently repair the reaction.
    fn validate(&self, py: Python<'_>) -> PyResult<ReactionValidationReport> {
        self.inner
            .validate()
            .map(|inner| ReactionValidationReport { inner })
            .map_err(|e| validation_error(py, &e))
    }
    /// Return a ReactionValidationReport with warnings and errors for template and atom-map consistency; do not silently repair the reaction. Uses the supplied configuration object.
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
    /// Return an initialized Reaction with compiled reactant matchers; leave the source Reaction unchanged.
    fn with_initialized(&self, py: Python<'_>) -> PyResult<Reaction> {
        self.inner
            .with_initialized()
            .map(|inner| Reaction { inner })
            .map_err(|e| initialization_error(py, &e))
    }
    /// Return an initialized Reaction with compiled reactant matchers; leave the source Reaction unchanged. Uses the supplied configuration object.
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
    /// Return a new Reaction with all agent templates removed.
    fn without_agents(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_agents(),
        }
    }
    /// Return a ReactionTemplateRemoval containing the new reaction and removed reactant templates, according to the unmapped-atom threshold.
    fn without_unmapped_reactants(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_unmapped_reactants(),
        }
    }
    /// Return a ReactionTemplateRemoval containing the new reaction and removed reactant templates, according to the unmapped-atom threshold. Uses the supplied configuration object.
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
    /// Return a ReactionTemplateRemoval containing the new reaction and removed product templates, according to the unmapped-atom threshold.
    fn without_unmapped_products(&self) -> ReactionTemplateRemoval {
        ReactionTemplateRemoval {
            inner: self.inner.without_unmapped_products(),
        }
    }
    /// Return a ReactionTemplateRemoval containing the new reaction and removed product templates, according to the unmapped-atom threshold. Uses the supplied configuration object.
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
/// New reaction plus the query templates removed by a template-filtering operation.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionTemplateRemoval {
    inner: ck::ReactionTemplateRemoval,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionTemplateRemoval {
    /// Reaction value produced by this result.
    #[getter]
    fn reaction(&self) -> Reaction {
        Reaction {
            inner: self.inner.reaction().clone(),
        }
    }
    /// QueryGraph templates removed by the filtering operation.
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
/// Result of single-reactant reaction application: the resulting molecule and whether it changed.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionApplyResult {
    inner: ck::ReactionApplyResult,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionApplyResult {
    /// Return the molecule produced by this operation; the input molecule remains independently owned.
    #[getter]
    fn molecule(&self) -> Molecule {
        Molecule {
            inner: self.inner.molecule().clone(),
        }
    }
    /// Whether reaction application changed the molecule.
    #[getter]
    fn changed(&self) -> bool {
        self.inner.changed()
    }
}
/// Reaction-template validation warnings and errors; is_valid reports whether errors are absent.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionValidationReport {
    inner: ck::ReactionValidationReport,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationReport {
    /// Reaction validation warning issues in report order.
    #[getter]
    fn warnings(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .warnings()
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    /// Return the reaction-template validation issues classified as errors.
    #[getter]
    fn errors(&self) -> Vec<ReactionValidationIssue> {
        self.inner
            .errors()
            .iter()
            .cloned()
            .map(|inner| ReactionValidationIssue { inner })
            .collect()
    }
    /// Whether this result satisfies its domain validity conditions.
    #[getter]
    fn is_valid(&self) -> bool {
        self.inner.is_valid()
    }
    /// Number of reaction validation warnings.
    #[getter]
    fn num_warnings(&self) -> usize {
        self.inner.num_warnings()
    }
    /// Number of reaction validation errors.
    #[getter]
    fn num_errors(&self) -> usize {
        self.inner.num_errors()
    }
}
/// Structured reaction validation issue, with severity, template role/index, atom/map context and detail.
#[cfg_attr(feature = "stubgen", gen_stub_pyclass)]
#[pyclass(module = "cosmolkit", frozen)]
pub(crate) struct ReactionValidationIssue {
    inner: ck::ReactionValidationIssue,
}
#[cfg_attr(feature = "stubgen", gen_stub_pymethods)]
#[pymethods]
impl ReactionValidationIssue {
    /// Classification/discriminant of this value as defined by its owning type.
    #[getter]
    fn kind(&self) -> ReactionValidationIssueKind {
        self.inner.kind().into()
    }
    /// Validation issue severity.
    #[getter]
    fn severity(&self) -> ReactionValidationSeverity {
        self.inner.severity().into()
    }
    /// Reactant, product or agent template role associated with this issue.
    #[getter]
    fn role(&self) -> Option<ReactionRole> {
        self.inner.role().map(Into::into)
    }
    /// Zero-based template index associated with this issue, when available.
    #[getter]
    fn template(&self) -> Option<usize> {
        self.inner.template()
    }
    /// Atom index associated with the validation issue, when available.
    #[getter]
    fn atom(&self) -> Option<usize> {
        self.inner.atom().map(|v| v.index())
    }
    /// Atom-map number associated with this issue, when available.
    #[getter]
    fn map(&self) -> Option<i32> {
        self.inner.map()
    }
    /// Atom-map numbers associated with this issue.
    #[getter]
    fn maps(&self) -> Vec<i32> {
        self.inner.maps().to_vec()
    }
    /// Human-readable context for this issue.
    #[getter]
    fn detail(&self) -> String {
        self.inner.detail().to_owned()
    }
}
pyo3::create_exception!(
    cosmolkit,
    ReactionModelError,
    PyValueError,
    "The reaction templates or their references are structurally invalid."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionParseError,
    PyValueError,
    "Reaction SMARTS, SMILES or reaction-block text could not be parsed."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionRunError,
    PyValueError,
    "A reaction could not be executed with the supplied reactants and options."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionApplyError,
    PyValueError,
    "A reaction could not be applied to the molecule."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionProductError,
    PyValueError,
    "A reaction product could not be constructed or finalized."
);
pub(crate) fn product_error(py: Python<'_>, source: &ck::ReactionProductError) -> PyErr {
    use ck::ReactionProductError as E;
    let kind = match source {
        E::StereoGroup(_) => "StereoGroup",
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
pyo3::create_exception!(
    cosmolkit,
    ReactionWriteError,
    PyValueError,
    "The reaction could not be serialized in the requested notation or format."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionValidationError,
    PyValueError,
    "Reaction-template validation could not be completed."
);
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
pyo3::create_exception!(
    cosmolkit,
    ReactionInitializationError,
    PyValueError,
    "Reaction templates could not be initialized for execution."
);
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
/// Parse reaction SMIRKS text into a Reaction using the canonical reaction parser.
#[cfg_attr(feature = "stubgen", gen_stub_pyfunction)]
#[pyfunction]
fn parse_smirks(py: Python<'_>, text: &str) -> PyResult<Reaction> {
    Reaction::from_smirks(py, text)
}
/// Parse reaction SMIRKS text into a Reaction using the canonical reaction parser. Uses the supplied configuration object.
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
