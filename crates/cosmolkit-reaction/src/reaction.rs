use crate::{ReactionModelError, ReactionRole};
use cosmolkit_model::QueryGraph;
use cosmolkit_search::SubstructMatchParams;
use std::collections::BTreeMap;

/// Immutable detached reaction value. Construction does not chemically
/// validate or initialize its templates. Preparation returns a new value.
#[derive(Debug)]
pub struct Reaction {
    pub(crate) reactants: Vec<QueryGraph>,
    pub(crate) products: Vec<QueryGraph>,
    pub(crate) agents: Vec<QueryGraph>,
    pub(crate) needs_init: bool,
    pub(crate) implicit_properties: bool,
    pub(crate) match_params: SubstructMatchParams,
    pub(crate) properties: BTreeMap<String, String>,
}

impl Clone for Reaction {
    fn clone(&self) -> Self {
        // RDKit❗✔️: void copy(const ChemicalReaction &other) {
        // RDKit❗✔️:     RDProps::operator=(other);
        // RDKit❗✔️:     df_needsInit = other.df_needsInit;
        // RDKit❗✔️:     df_implicitProperties = other.df_implicitProperties;
        // RDKit❗✔️:     m_reactantTemplates.clear();
        // RDKit❗✔️:     m_reactantTemplates.reserve(other.m_reactantTemplates.size());
        // RDKit❗✔️:     for (ROMOL_SPTR reactant_template : other.m_reactantTemplates) {
        // RDKit❗✔️:       m_reactantTemplates.emplace_back(new RWMol(*reactant_template));
        // RDKit❗✔️:     }
        // RDKit❗✔️:     m_productTemplates.clear();
        // RDKit❗✔️:     m_productTemplates.reserve(other.m_productTemplates.size());
        // RDKit❗✔️:     for (ROMOL_SPTR product_template : other.m_productTemplates) {
        // RDKit❗✔️:       m_productTemplates.emplace_back(new RWMol(*product_template));
        // RDKit❗✔️:     }
        // RDKit❗✔️:     m_agentTemplates.clear();
        // RDKit❗✔️:     m_agentTemplates.reserve(other.m_agentTemplates.size());
        // RDKit❗✔️:     for (ROMOL_SPTR agent_template : other.m_agentTemplates) {
        // RDKit❗✔️:       m_agentTemplates.emplace_back(new RWMol(*agent_template));
        // RDKit❗✔️:     }
        // RDKit❗✔️:     d_substructParams = other.d_substructParams;
        Self {
            reactants: self.reactants.clone(),
            products: self.products.clone(),
            agents: self.agents.clone(),
            needs_init: self.needs_init,
            implicit_properties: self.implicit_properties,
            match_params: self.match_params.clone(),
            properties: self.properties.clone(),
        }
    }
}

impl Default for Reaction {
    fn default() -> Self {
        Self::new()
    }
}

impl Reaction {
    pub fn new() -> Self {
        // RDKit❗✔️: bool df_needsInit{true};
        // RDKit❗✔️:   bool df_implicitProperties{false};
        // RDKit❗✔️:   MOL_SPTR_VECT m_reactantTemplates, m_productTemplates, m_agentTemplates;
        // RDKit❗✔️:   SubstructMatchParameters d_substructParams;
        Self {
            reactants: Vec::new(),
            products: Vec::new(),
            agents: Vec::new(),
            needs_init: true,
            implicit_properties: false,
            match_params: SubstructMatchParams::default(),
            properties: BTreeMap::new(),
        }
    }

    pub fn from_templates(
        reactants: Vec<QueryGraph>,
        products: Vec<QueryGraph>,
        agents: Vec<QueryGraph>,
    ) -> Result<Self, ReactionModelError> {
        // RDKit❗✔️: ChemicalReaction() : RDProps() {}
        // CK structural boundary: each canonical query validates before storage;
        // this does not execute ChemicalReaction::validate or change needsInit.
        for (role, templates) in [
            (ReactionRole::Reactant, &reactants),
            (ReactionRole::Product, &products),
            (ReactionRole::Agent, &agents),
        ] {
            for (index, graph) in templates.iter().enumerate() {
                validate_template(graph, role, index)?;
            }
        }
        Ok(Self {
            reactants,
            products,
            agents,
            ..Self::new()
        })
    }

    pub fn reactant_templates(&self) -> &[QueryGraph] {
        // RDKit❗✔️: const MOL_SPTR_VECT &getReactants() const {
        // RDKit❗✔️:     return this->m_reactantTemplates;
        &self.reactants
    }

    pub fn reactant_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        // RDKit❗✔️: const MOL_SPTR_VECT &getReactants() const {
        // RDKit❗✔️:     return this->m_reactantTemplates;
        self.reactants
            .get(index)
            .ok_or(ReactionModelError::TemplateIndex {
                role: ReactionRole::Reactant,
                index,
                count: self.reactants.len(),
            })
    }

    pub fn num_reactant_templates(&self) -> usize {
        // RDKit❗✔️: unsigned int getNumReactantTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_reactantTemplates.size());
        self.reactants.len()
    }

    pub fn with_reactant_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // RDKit❗❌: unsigned int addReactantTemplate(ROMOL_SPTR mol) {
        // RDKit❗❌:     this->df_needsInit = true;
        // RDKit❗❌:     this->m_reactantTemplates.push_back(mol);
        // RDKit❗❌:     return rdcast<unsigned int>(this->m_reactantTemplates.size());
        validate_template(&template, ReactionRole::Reactant, self.reactants.len())?;
        // Approved immutable value projection copies the reaction before append.
        let mut result = self.clone();
        result.needs_init = true;
        result.reactants.push(template);
        Ok(result)
    }

    pub fn product_templates(&self) -> &[QueryGraph] {
        // RDKit❗✔️: const MOL_SPTR_VECT &getProducts() const { return this->m_productTemplates;
        &self.products
    }

    pub fn product_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        // RDKit❗✔️: const MOL_SPTR_VECT &getProducts() const { return this->m_productTemplates;
        self.products
            .get(index)
            .ok_or(ReactionModelError::TemplateIndex {
                role: ReactionRole::Product,
                index,
                count: self.products.len(),
            })
    }

    pub fn num_product_templates(&self) -> usize {
        // RDKit❗✔️: unsigned int getNumProductTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_productTemplates.size());
        self.products.len()
    }

    pub fn with_product_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // RDKit❗❌: unsigned int addProductTemplate(ROMOL_SPTR mol) {
        // RDKit❗❌:     this->m_productTemplates.push_back(mol);
        // RDKit❗❌:     return rdcast<unsigned int>(this->m_productTemplates.size());
        validate_template(&template, ReactionRole::Product, self.products.len())?;
        // Approved immutable value projection copies the reaction before append.
        let mut result = self.clone();
        result.products.push(template);
        Ok(result)
    }

    pub fn agent_templates(&self) -> &[QueryGraph] {
        // RDKit❗✔️: const MOL_SPTR_VECT &getAgents() const { return this->m_agentTemplates; }
        &self.agents
    }

    pub fn agent_template(&self, index: usize) -> Result<&QueryGraph, ReactionModelError> {
        // RDKit❗✔️: const MOL_SPTR_VECT &getAgents() const { return this->m_agentTemplates; }
        self.agents
            .get(index)
            .ok_or(ReactionModelError::TemplateIndex {
                role: ReactionRole::Agent,
                index,
                count: self.agents.len(),
            })
    }

    pub fn num_agent_templates(&self) -> usize {
        // RDKit❗✔️: unsigned int getNumAgentTemplates() const {
        // RDKit❗✔️:     return rdcast<unsigned int>(this->m_agentTemplates.size());
        self.agents.len()
    }

    pub fn with_agent_template(&self, template: QueryGraph) -> Result<Self, ReactionModelError> {
        // RDKit❗❌: unsigned int addAgentTemplate(ROMOL_SPTR mol) {
        // RDKit❗❌:     this->m_agentTemplates.push_back(mol);
        // RDKit❗❌:     return rdcast<unsigned int>(this->m_agentTemplates.size());
        validate_template(&template, ReactionRole::Agent, self.agents.len())?;
        // Approved immutable value projection copies the reaction before append.
        let mut result = self.clone();
        result.agents.push(template);
        Ok(result)
    }

    pub fn is_initialized(&self) -> bool {
        // RDKit❗✔️: bool isInitialized() const { return !df_needsInit; }
        !self.needs_init
    }

    pub fn implicit_properties(&self) -> bool {
        // RDKit❗✔️: bool getImplicitPropertiesFlag() const { return df_implicitProperties; }
        self.implicit_properties
    }

    pub fn with_implicit_properties(&self, enabled: bool) -> Self {
        // RDKit❗❌: void setImplicitPropertiesFlag(bool val) { df_implicitProperties = val; }
        let mut result = self.clone();
        result.implicit_properties = enabled;
        result
    }

    pub fn match_params(&self) -> &SubstructMatchParams {
        // RDKit❗✔️: const SubstructMatchParameters &getSubstructParams() const {
        // RDKit❗✔️:     return d_substructParams;
        &self.match_params
    }

    pub fn with_match_params(
        &self,
        params: &SubstructMatchParams,
    ) -> Result<Self, ReactionModelError> {
        // RDKit❗❌: SubstructMatchParameters &getSubstructParams() { return d_substructParams; }
        // Source exposes direct configuration assignment without added constraints.
        let mut result = self.clone();
        result.match_params = params.clone();
        Ok(result)
    }
}

fn validate_template(
    graph: &QueryGraph,
    role: ReactionRole,
    template: usize,
) -> Result<(), ReactionModelError> {
    graph
        .validate()
        .map_err(|source| ReactionModelError::Template {
            role,
            template,
            source,
        })
}
