//! Reaction result transport; chemistry remains in the detached reaction owner.

use crate::{
    Molecule, PropertyText, ReactionInitializationError, ReactionModelError, ReactionParseError,
    ReactionParseParams, ReactionTemplateRemovalParams, ReactionValidationError,
    ReactionValidationParams, ReactionValidationReport, ReactionWriteError, ReactionWriteParams,
};
use cosmolkit_model::QueryGraph;
use cosmolkit_reaction::SubstructMatchParams;

/// Public reaction object. Detached template storage and chemistry belong to
/// `cosmolkit-reaction`; live molecule execution belongs to this facade.
#[derive(Debug, Clone, Default)]
pub struct Reaction {
    inner: cosmolkit_reaction::Reaction,
}

impl Reaction {
    pub fn new() -> Self {
        Self::default()
    }
    pub fn from_smirks(text: &str) -> Result<Self, ReactionParseError> {
        parse_smirks(text)
    }
    pub fn from_smirks_with_params(
        text: &str,
        params: &ReactionParseParams,
    ) -> Result<Self, ReactionParseError> {
        parse_smirks_with_params(text, params)
    }
    pub fn from_templates(
        reactants: Vec<QueryGraph>,
        products: Vec<QueryGraph>,
        agents: Vec<QueryGraph>,
    ) -> Result<Self, ReactionModelError> {
        cosmolkit_reaction::Reaction::from_templates(reactants, products, agents)
            .map(|inner| Self { inner })
    }
    pub fn reactant_templates(&self) -> &[QueryGraph] {
        self.inner.reactant_templates()
    }
    pub fn product_templates(&self) -> &[QueryGraph] {
        self.inner.product_templates()
    }
    pub fn agent_templates(&self) -> &[QueryGraph] {
        self.inner.agent_templates()
    }
    pub fn reactant_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        self.inner.reactant_template(index)
    }
    pub fn product_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        self.inner.product_template(index)
    }
    pub fn agent_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        self.inner.agent_template(index)
    }
    pub fn num_reactant_templates(&self) -> usize {
        self.inner.num_reactant_templates()
    }
    pub fn num_product_templates(&self) -> usize {
        self.inner.num_product_templates()
    }
    pub fn num_agent_templates(&self) -> usize {
        self.inner.num_agent_templates()
    }
    pub fn with_reactant_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        self.inner
            .with_reactant_template(template)
            .map(|inner| Self { inner })
    }
    pub fn with_product_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        self.inner
            .with_product_template(template)
            .map(|inner| Self { inner })
    }
    pub fn with_agent_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        self.inner
            .with_agent_template(template)
            .map(|inner| Self { inner })
    }
    pub fn is_initialized(&self) -> bool {
        self.inner.is_initialized()
    }
    pub fn implicit_properties(&self) -> bool {
        self.inner.implicit_properties()
    }
    pub fn with_implicit_properties(&self, enabled: bool) -> Self {
        Self {
            inner: self.inner.with_implicit_properties(enabled),
        }
    }
    pub fn match_params(&self) -> &SubstructMatchParams {
        self.inner.match_params()
    }
    pub fn with_match_params(
        &self,
        params: &SubstructMatchParams,
    ) -> Result<Self, ReactionModelError> {
        self.inner
            .with_match_params(params)
            .map(|inner| Self { inner })
    }
    pub fn to_smirks(&self) -> Result<PropertyText, ReactionWriteError> {
        self.inner.to_smirks()
    }
    pub fn to_smirks_with_params(
        &self,
        params: &ReactionWriteParams,
    ) -> Result<PropertyText, ReactionWriteError> {
        self.inner.to_smirks_with_params(params)
    }
    pub fn to_cx_smirks(&self) -> Result<PropertyText, ReactionWriteError> {
        self.inner.to_cx_smirks()
    }
    pub fn to_cx_smirks_with_params(
        &self,
        params: &ReactionWriteParams,
    ) -> Result<PropertyText, ReactionWriteError> {
        self.inner.to_cx_smirks_with_params(params)
    }
    pub fn validate(&self) -> Result<ReactionValidationReport, ReactionValidationError> {
        self.inner.validate()
    }
    pub fn validate_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<ReactionValidationReport, ReactionValidationError> {
        self.inner.validate_with_params(params)
    }
    pub fn with_initialized(&self) -> Result<Self, ReactionInitializationError> {
        self.inner.with_initialized().map(|inner| Self { inner })
    }
    pub fn with_initialized_with_params(
        &self,
        params: &ReactionValidationParams,
    ) -> Result<Self, ReactionInitializationError> {
        self.inner
            .with_initialized_with_params(params)
            .map(|inner| Self { inner })
    }
    pub fn without_agents(&self) -> ReactionTemplateRemoval {
        self.inner.without_agents().into()
    }
    pub fn without_unmapped_reactants(&self) -> ReactionTemplateRemoval {
        self.inner.without_unmapped_reactants().into()
    }
    pub fn without_unmapped_reactants_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        self.inner
            .without_unmapped_reactants_with_params(params)
            .into()
    }
    pub fn without_unmapped_products(&self) -> ReactionTemplateRemoval {
        self.inner.without_unmapped_products().into()
    }
    pub fn without_unmapped_products_with_params(
        &self,
        params: &ReactionTemplateRemovalParams,
    ) -> ReactionTemplateRemoval {
        self.inner
            .without_unmapped_products_with_params(params)
            .into()
    }
    pub(crate) fn detached_mut(&mut self) -> &mut cosmolkit_reaction::Reaction {
        &mut self.inner
    }
}

pub fn parse_smirks(text: &str) -> Result<Reaction, ReactionParseError> {
    cosmolkit_reaction::parse_smirks(text).map(|inner| Reaction { inner })
}
pub fn parse_smirks_with_params(
    text: &str,
    params: &ReactionParseParams,
) -> Result<Reaction, ReactionParseError> {
    cosmolkit_reaction::parse_smirks_with_params(text, params).map(|inner| Reaction { inner })
}

#[derive(Debug, Clone)]
pub struct ReactionTemplateRemoval {
    pub reaction: Reaction,
    pub removed_templates: Vec<QueryGraph>,
}
impl ReactionTemplateRemoval {
    pub fn reaction(&self) -> &Reaction {
        &self.reaction
    }
    pub fn removed_templates(&self) -> &[QueryGraph] {
        &self.removed_templates
    }
}
impl From<cosmolkit_reaction::ReactionTemplateRemoval> for ReactionTemplateRemoval {
    fn from(value: cosmolkit_reaction::ReactionTemplateRemoval) -> Self {
        Self {
            reaction: Reaction {
                inner: value.reaction,
            },
            removed_templates: value.removed_templates,
        }
    }
}

/// The source bool and the molecule finalized by the sole operation runtime.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionApplyResult {
    pub(crate) molecule: Molecule,
    pub(crate) changed: bool,
}

impl ReactionApplyResult {
    pub fn molecule(&self) -> &Molecule {
        &self.molecule
    }
    pub fn changed(&self) -> bool {
        self.changed
    }
}

impl From<(Molecule, bool)> for ReactionApplyResult {
    fn from((molecule, changed): (Molecule, bool)) -> Self {
        // Preserve the algorithm's bool, including a true result with no
        // graph difference. Conversion owns the already finalized molecule;
        // it cannot validate, commit, or infer chemistry from graph equality.
        Self { molecule, changed }
    }
}

#[cfg(feature = "cap-reaction")]
pub(crate) fn assemble_product_sets(
    molecules: Vec<Molecule>,
    lengths: Vec<usize>,
) -> Result<Vec<Vec<Molecule>>, crate::OperationError> {
    let mut remaining = molecules.into_iter();
    let mut sets = Vec::with_capacity(lengths.len());
    for length in lengths {
        if length > remaining.len() {
            return Err(crate::OperationError::InvalidAlgorithmResult {
                operation: "reaction_products",
                field: "product set length",
                actual: length,
                expected: remaining.len(),
            });
        }
        sets.push(remaining.by_ref().take(length).collect());
    }
    if remaining.len() != 0 {
        return Err(crate::OperationError::InvalidAlgorithmResult {
            operation: "reaction_products",
            field: "ungrouped products",
            actual: remaining.len(),
            expected: 0,
        });
    }
    Ok(sets)
}
